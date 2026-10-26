"""Offset-preserving global DFTHC mathematical refit. No circuit construction."""
from __future__ import annotations
import argparse
import dataclasses
import hashlib
import importlib.util
import json
from pathlib import Path
import time

import numpy as np
import scipy
import scipy.linalg as la
import scipy.optimize as opt
from scipy.special import softmax

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent
LOCAL=HERE/'local'
SHAPE=(16,32,20)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def load_helpers():
    spec=importlib.util.spec_from_file_location('gram_study',PROJECT/'dfthc_gram/study.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def pack(matrix):
    n=matrix.shape[-1];i,j=np.tril_indices(n)
    return matrix[...,i,j]*np.where(i==j,1.,np.sqrt(2.))

def unpack(vector,n):
    i,j=np.tril_indices(n);values=vector/np.where(i==j,1.,np.sqrt(2.))
    out=np.zeros(vector.shape[:-1]+(n,n));out[...,i,j]=values;out[...,j,i]=values
    return out

def projector_rows(u):return pack(u[...,None]*u[...,None,:])

def polar(z):
    if z.ndim!=2 or z.shape[0]<z.shape[1] or not np.isrealobj(z) or not np.isfinite(z).all():
        raise ValueError('Polar input must be finite, tall and two-dimensional.')
    left,values,vt=la.svd(z,full_matrices=False,check_finite=False)
    if values[-1]<1e-10*max(values[0],1.):raise ValueError('Polar frame is rank deficient.')
    return left@vt,values,vt.T

def polar_pullback(frame,gradient,values,right):
    a=frame.T@gradient;rhs=a-a.T
    k=right@((right.T@rhs@right)/(values[:,None]+values[None,:]))@right.T
    inverse=(right/values)@right.T
    return (gradient-frame@a)@inverse+frame@k

def gram_loss(x,target,target_squared):
    """Rows of x are symmetric-pair factors; target is symmetric."""
    gram=x@x.T;tx=x@target
    squared=float(np.sum(gram*gram)-2*np.sum(tx*x)+target_squared)
    if squared < -1e-9:raise ArithmeticError('Unresolved global residual cancellation.')
    return max(squared,0.),4*(gram@x-tx)

class GlobalModel:
    def __init__(self,target,offset,tau,shape=SHAPE):
        if (offset.ndim!=2 or offset.shape[0]!=offset.shape[1] or
            not np.isrealobj(target) or not np.isrealobj(offset) or
            not np.isscalar(tau) or not np.isrealobj(tau) or not np.isfinite(tau)):
            raise ValueError('Global model requires real matrices and a finite offset norm.')
        r,b,c=shape;n=offset.shape[0]
        if target.shape!=(n*(n+1)//2,)*2 or c<2 or c-1>b or r*b<n or tau<=0:
            raise ValueError('Invalid global model dimensions or offset norm.')
        if not np.isfinite(target).all() or not np.allclose(target,target.T):
            raise ValueError('Target must be finite and symmetric.')
        if not np.isfinite(offset).all() or not np.allclose(offset,offset.T):
            raise ValueError('Offset must be finite and symmetric.')
        self.l=la.cholesky(-offset,lower=False,check_finite=True)
        self.target=target;self.offset=offset;self.tau=float(tau)
        self.r,self.b,self.c,self.n=r,b,c,n
        self.d=target.shape[0];self.target_squared=float(np.sum(target**2))
        self.scale=max(self.target_squared,1.)
        self.z_size=r*b*n;self.w_size=r*b*(c-1)

    def split(self,raw):
        if raw.shape!=(self.z_size+self.w_size+self.r,):raise ValueError('Parameter length.')
        z=raw[:self.z_size].reshape(self.r*self.b,self.n)
        w=raw[self.z_size:self.z_size+self.w_size].reshape(self.r,self.b,self.c-1)
        gamma=raw[-self.r:]
        return z,w,gamma

    def materialize(self,raw):
        z,free,gamma=self.split(raw);frame,values,right=polar(z)
        s=(frame@self.l).reshape(self.r,self.b,self.n)
        norm=np.linalg.norm(s,axis=-1)
        if np.min(norm)<1e-12:raise ValueError('Zero projector row.')
        u=s/norm[...,None];probability=softmax(gamma)
        alpha=np.sqrt(self.tau*probability)
        if np.min(alpha)<1e-12:raise ValueError('Unresolved offset allocation.')
        w=np.concatenate([(-norm**2/alpha[:,None])[...,None],free],axis=2)
        v=np.zeros((self.r,self.c));v[:,0]=alpha
        return u,w,v,(frame,values,right,s,norm,probability)

    def value_gradient(self,raw):
        u,w,v,cache=self.materialize(raw)
        frame,values,right,s,norm,probability=cache;alpha=v[:,0]
        p=projector_rows(u);x=np.einsum('rbc,rbd->rcd',w,p,optimize=True)
        squared,gx=gram_loss(x.reshape(-1,self.d),self.target,self.target_squared)
        gx=gx.reshape(self.r,self.c,self.d)
        gw=np.einsum('rbd,rcd->rbc',p,gx[:,1:],optimize=True)
        packed=np.einsum('rbc,rcd->rbd',w[:,:,1:],gx[:,1:],optimize=True)
        gu=2*np.einsum('rbij,rbj->rbi',unpack(packed,self.n),u,optimize=True)
        gs=(gu-u*np.sum(gu*u,axis=-1,keepdims=True))/norm[...,None]
        gs-=2*np.einsum('rbi,rij->rbj',s,unpack(gx[:,0],self.n),optimize=True)/alpha[:,None,None]
        gf=gs.reshape(-1,self.n)@self.l.T
        gz=polar_pullback(frame,gf,values,right)
        ga=-np.sum(gx[:,0]*x[:,0],axis=1)/alpha
        gg=(alpha*ga-probability*np.sum(alpha*ga))/2
        gradient=np.concatenate([gz.ravel(),gw.ravel(),gg])
        return squared/self.scale,gradient/self.scale

    def direct_error(self,raw):
        u,w,_,_=self.materialize(raw)
        x=np.einsum('rbc,rbd->rcd',w,projector_rows(u),optimize=True).reshape(-1,self.d)
        return float(np.linalg.norm(x.T@x-self.target))

class FixedModel:
    """Decoded directions and offset coefficients fixed; all free cores refit globally."""
    def __init__(self,model,u,coefficients):
        self.model=model;self.u=u;self.p=projector_rows(u);self.coefficients=coefficients

    def materialize(self,raw):
        m=self.model;free=raw[:-m.r].reshape(m.r,m.b,m.c-1)
        prob=softmax(raw[-m.r:]);alpha=np.sqrt(m.tau*prob)
        w=np.concatenate([(self.coefficients/alpha[:,None])[...,None],free],axis=2)
        v=np.zeros((m.r,m.c));v[:,0]=alpha
        return w,v,prob

    def value_gradient(self,raw):
        m=self.model;w,v,prob=self.materialize(raw)
        x=np.einsum('rbc,rbd->rcd',w,self.p,optimize=True)
        squared,gx=gram_loss(x.reshape(-1,m.d),m.target,m.target_squared)
        gx=gx.reshape(m.r,m.c,m.d)
        gw=np.einsum('rbd,rcd->rbc',self.p,gx[:,1:],optimize=True)
        alpha=v[:,0];ga=-np.sum(gx[:,0]*x[:,0],axis=1)/alpha
        gg=(alpha*ga-prob*np.sum(alpha*ga))/2
        return squared/m.scale,np.concatenate([gw.ravel(),gg])/m.scale

def core_initialize(model,frame,sweeps=4):
    r,b,c,n=model.r,model.b,model.c,model.n
    raw=np.concatenate([frame.ravel(),np.zeros(model.w_size+r)])
    u,w,_,_=model.materialize(raw);p=projector_rows(u)
    inverses=[];q=[]
    for block in p:
        ev,vec=la.eigh(block@block.T)
        if ev[0]<1e-12:raise ValueError('Initial projector dictionary is singular.')
        inverse=(vec/np.sqrt(ev))@vec.T
        inverses.append(inverse);q.append(block.T@inverse)
    q=np.stack(q);inverses=np.stack(inverses)
    tq=np.einsum('rdi,de->rie',q,model.target,optimize=True)
    projected=np.einsum('rie,rej->rij',tq,q,optimize=True)
    x=np.einsum('rbc,rbd->rcd',w,p,optimize=True)
    history=[]
    for _ in range(sweeps):
        for block in range(r):
            current=x.reshape(-1,model.d)@q[block]
            own=x[block,1:]@q[block]
            residual=projected[block]-current.T@current+own.T@own
            ev,vec=la.eigh((residual+residual.T)/2)
            indices=np.arange(b-1,max(b-c,-1),-1)
            factors=vec[:,indices]*np.sqrt(np.maximum(ev[indices],0))
            w[block,:,1:]=inverses[block]@factors
            x[block]=w[block].T@p[block]
        squared,_=gram_loss(x.reshape(-1,model.d),model.target,model.target_squared)
        history.append(float(np.sqrt(squared)))
    if np.any(np.diff(history)>1e-8):raise ArithmeticError('Core update increased residual.')
    raw[model.z_size:model.z_size+model.w_size]=w[:,:,1:].ravel()
    return raw,history

def fit(model,raw,iterations,label):
    history=[];last=[None]
    def objective(x):
        value,gradient=model.value_gradient(x);last[:]=[(x.copy(),value)]
        return value,gradient
    scale=model.scale if isinstance(model,GlobalModel) else model.model.scale
    def callback(x):
        cached=last[0]
        value=cached[1] if cached is not None and np.array_equal(cached[0],x) else objective(x)[0]
        history.append(float(np.sqrt(max(value*scale,0))))
        if len(history)%50==0:print(json.dumps({'fit':label,'iteration':len(history),'error':history[-1]}),flush=True)
    started=time.monotonic()
    allocation_count=model.r if isinstance(model,GlobalModel) else model.model.r
    bounds=[(None,None)]*(len(raw)-allocation_count)+[(-4.,4.)]*allocation_count
    result=opt.minimize(objective,raw,method='L-BFGS-B',jac=True,callback=callback,bounds=bounds,
        options={'maxiter':iterations,'maxfun':20000,'maxls':30,'maxcor':12,'ftol':1e-14,'gtol':1e-9})
    value,gradient=model.value_gradient(result.x)
    report={'iteration_cap':iterations,'iterations':int(result.nit),'evaluations':int(result.nfev),
        'optimizer_success':bool(result.success),'termination':str(result.message),
        'error_history':history,'contraction_error':float(np.sqrt(max(value*scale,0))),
        'gradient_infinity_norm':float(np.max(np.abs(gradient))),
        'allocation_logit_bounds':[-4.,4.],
        'elapsed_seconds':time.monotonic()-started,'global_optimality_claim':False}
    return result.x,report

def seed_frames(model,u0,w0,pair):
    helper=load_helpers();r=model.r;b=model.b
    x0=np.einsum('rbc,rbd->rcd',w0,projector_rows(u0),optimize=True)
    cross=(x0.reshape(-1,model.d)@x0.reshape(-1,model.d).T).reshape(65,30,65,30)
    residual=np.sum(cross**2,axis=(1,3));chosen=[]
    for _ in range(r):
        diagonal=np.diag(residual).copy();diagonal[chosen]=-np.inf
        index=int(diagonal.argmax());pivot=residual[index,index]
        if pivot<=1e-15:raise ValueError('Unresolved release-block seed.')
        chosen.append(index);residual-=np.outer(residual[:,index],residual[index,:])/pivot
        residual=(residual+residual.T)/2
    release,_=helper.greedy_dictionary(u0[chosen],w0[chosen],b)
    release_frame=polar(release.reshape(r*b,model.n))[0]
    ev,vec=la.eigh(pair,subset_by_index=[model.d-r,model.d-1])
    spectral=[]
    for column in vec[:,::-1].T:
        matrix=unpack(column,model.n);values,vectors=la.eigh(matrix)
        keep=np.argsort(np.abs(values))[-b:]
        spectral.extend(vectors[:,keep].T)
    spectral_frame=polar(np.array(spectral))[0]
    return [('released',release_frame),('spectral',spectral_frame)]

def official_candidate(helper,verifier,maths,inputs,model,u,w,v,beta,name,two_body_words=None):
    tree=inputs.params['tree_edges'];qm=1<<(beta-1)
    words2=helper.encode(u.reshape(-1,model.n),tree,beta) if two_body_words is None else two_body_words
    decoded=maths.reconstruct_u_vectors(words2/qm,tree).reshape(u.shape)
    params={**inputs.params,'U':decoded,'W':w,'H2_I':v+w.sum(axis=1)}
    nw,nv,ns,nb=verifier.native_inputs(params)
    if np.max(np.abs(nv-v))>1e-10:raise ArithmeticError('Native-offset serialization changed offsets.')
    candidate=dataclasses.replace(inputs,params=params,w=nw,W=nv,S=ns,beta_1=nb)
    x=maths.reconstruct_x(decoded,nw)
    h1=maths.build_h1_prime(inputs.ref['h1_thc'],ns,nb,inputs.ref['electrons'],nv,x)
    _,vectors=la.eigh(h1);words1=helper.encode(vectors.T,tree,beta)
    words=np.concatenate([words1,words2])
    metrics=verifier.evaluate_words(candidate,words,tree,beta)
    # Official metric code may report NaN normalization outside the SOS domain.
    metrics={k:(None if isinstance(value,(float,np.floating)) and not np.isfinite(value) else value)
             for k,value in metrics.items()}
    params['angle_words']=words;params['tree_edges']=tree
    np.savez(LOCAL/(name+'.npz'),**params)
    resource=None if metrics['lambda_eff'] is None else helper.resources(verifier,SHAPE,metrics['lambda_eff'],beta)
    contraction=np.einsum('rc,rcij->ij',nv,x,optimize=True)
    return {'name':name,'shape':list(SHAPE),'precision_bits':beta,'metrics':metrics,
        'metric_gates_pass':False if metrics['lambda_eff'] is None else bool(verifier.passes_accuracy(metrics)),
        'offset_matrix_drift':float(np.linalg.norm(contraction-model.offset)),
        'offset_norm_drift':float(abs(np.sum(nv**2)-model.tau)),
        'private_candidate_sha256':sha(LOCAL/(name+'.npz')),
        'whole_task_estimate':resource,'valid_competition_score':False},decoded

def gauge_value_gradient(raw,w,v):
    r,b,c=w.shape;k=c-1;rotations=[];pullbacks=[]
    for block in raw.reshape(r,k,k):
        rotation,values,right=polar(block)
        rotations.append(rotation);pullbacks.append((values,right))
    rotations=np.stack(rotations)
    rotated=np.concatenate([w[:,:,:1],w[:,:,1:]@rotations],axis=2)
    total=v+rotated.sum(axis=1)
    norms=np.abs(rotated).sum(axis=1)+np.abs(total)
    value=float(np.sum(norms**2)/4)
    gw=norms[:,None,:]*(np.sign(rotated)+np.sign(total[:,None,:]))/2
    go=w[:,:,1:].swapaxes(-1,-2)@gw[:,:,1:]
    gradient=np.stack([polar_pullback(q,g,*cache) for q,g,cache in zip(rotations,go,pullbacks)])
    return value,gradient.ravel(),rotated

def optimize_gauge(w,v):
    k=w.shape[2]-1;start=np.broadcast_to(np.eye(k),(len(w),k,k)).copy().ravel()
    initial=gauge_value_gradient(start,w,v)[0]
    def objective(raw):
        value,gradient,_=gauge_value_gradient(raw,w,v)
        return value/max(initial,1.),gradient/max(initial,1.)
    result=opt.minimize(objective,start,jac=True,method='L-BFGS-B',
        options={'maxiter':150,'maxfun':5000,'maxls':30,'ftol':1e-12,'gtol':1e-8})
    value,_,rotated=gauge_value_gradient(result.x,w,v)
    if value>initial+1e-8:raise ArithmeticError('Copy gauge increased normalization.')
    return rotated,{'iteration_cap':150,'iterations':int(result.nit),
        'optimizer_success':bool(result.success),'termination':str(result.message),
        'lambda_h2_before':initial,'lambda_h2_after':value,
        'global_optimality_claim':False}

def prepare(args):
    helper=load_helpers();verifier,maths,inputs,u,pair,provenance=helper.load(args.data)
    LOCAL.mkdir(exist_ok=True)
    x=maths.reconstruct_x(u,inputs.w);offset=np.einsum('rc,rcij->ij',inputs.W,x,optimize=True)
    tau=float(np.sum(inputs.W**2));model=GlobalModel(pair,offset,tau)
    seeds=seed_frames(model,u,inputs.w,pair)
    np.savez(LOCAL/'prepared.npz',offset=offset,tau=tau,**{name:frame for name,frame in seeds})
    out={'status':'PINNED_INPUTS','shape':list(SHAPE),'provenance':provenance,
        'helper_sha256':sha(PROJECT/'dfthc_gram/study.py'),
        'source_sha256':sha(Path(__file__)),'private_arrays_sha256':sha(LOCAL/'prepared.npz'),
        'offset_eigenvalue_range':[float(la.eigvalsh(offset)[0]),float(la.eigvalsh(offset)[-1])],
        'offset_norm_squared':tau,'offset_matrix_frobenius':float(np.linalg.norm(offset)),
        'baseline_cost':float(json.loads((PROJECT/'dfthc_c30/verified.json').read_text())['score']),
        'recorded_score_context':helper.score_context()}
    write_json(HERE/'inputs.json',out);print(json.dumps(out,indent=2),flush=True)

def run(args):
    helper=load_helpers();verifier,maths,inputs,_,pair,provenance=helper.load(args.data)
    meta=json.loads((HERE/'inputs.json').read_text());wl=json.loads((HERE/'symbolic.json').read_text())
    if meta['provenance']!=provenance or meta['source_sha256']!=sha(Path(__file__)):
        raise ValueError('Stale input audit.')
    if wl['status']!='PASS' or not all(sha(HERE/n)==value for n,value in wl['source_sha256'].items()):
        raise ValueError('WL symbolic gate is stale or failing.')
    if sha(LOCAL/'prepared.npz')!=meta['private_arrays_sha256']:raise ValueError('Prepared data changed.')
    with np.load(LOCAL/'prepared.npz',allow_pickle=False) as data:
        model=GlobalModel(pair,data['offset'],float(data['tau']))
        seeds=[(name,data[name].copy()) for name in ['released','spectral']]
    floating=[]
    for name,frame in seeds:
        start,initial=core_initialize(model,frame)
        print(json.dumps({'fit':name,'initialized_errors':initial}),flush=True)
        raw,report=fit(model,start,1200,name)
        u,w,v,cache=model.materialize(raw)
        x=maths.reconstruct_x(u,w);contraction=np.einsum('rc,rcij->ij',v,x,optimize=True)
        direct=model.direct_error(raw)
        if abs(direct-report['contraction_error'])>1e-7:raise ArithmeticError('Global loss mismatch.')
        offset_drift=float(np.linalg.norm(contraction-model.offset));tau_drift=float(abs(np.sum(v*v)-model.tau))
        if max(offset_drift,tau_drift)>1e-8:raise ArithmeticError('Floating offset contract failed.')
        np.savez(LOCAL/('float-'+name+'.npz'),raw=raw,u=u,w=w,v=v)
        h1=maths.build_h1_prime(inputs.ref['h1_thc'],inputs.S,inputs.beta_1,
                              inputs.ref['electrons'],v,x)
        _,one=la.eigh(h1)
        float_metrics=maths.evaluate_dfthc_vectors(np.concatenate([one.T,u.reshape(-1,model.n)]),
            w,v,inputs.S,inputs.beta_1,inputs.ref,energy_upper=inputs.energy_upper)
        floating.append({'name':name,'initial_core_sweeps':initial,'fit':report,
            'actual_pair_error':direct,'offset_matrix_drift':offset_drift,'offset_norm_drift':tau_drift,
            'unquantized_metrics':float_metrics,
            'private_sha256':sha(LOCAL/('float-'+name+'.npz'))})
        print(json.dumps({'floating':name,'actual_pair_error':direct,'offset_matrix_drift':offset_drift}),flush=True)
        write_json(HERE/'floating-progress.json',floating)
    selected=min(floating,key=lambda row:row['actual_pair_error'])['name']
    with np.load(LOCAL/('float-'+selected+'.npz'),allow_pickle=False) as data:
        raw=data['raw'].copy()
    u,w,v,cache=model.materialize(raw);coefficients=-(cache[4]**2)
    initial=np.concatenate([w[:,:,1:].ravel(),raw[-model.r:]])
    candidates=[]
    for beta in [10,12,14,16]:
        tree=inputs.params['tree_edges'];words=helper.encode(u.reshape(-1,model.n),tree,beta)
        decoded=maths.reconstruct_u_vectors(words/(1<<(beta-1)),tree).reshape(u.shape)
        fixed=FixedModel(model,decoded,coefficients)
        fitted,report=fit(fixed,initial,300,'decoded-'+str(beta))
        nw,nv,_=fixed.materialize(fitted)
        result,check_decoded=official_candidate(helper,verifier,maths,inputs,model,decoded,nw,nv,beta,
            'global-beta'+str(beta),two_body_words=words)
        if not np.allclose(check_decoded,decoded,rtol=0,atol=1e-13):raise ArithmeticError('Angle re-encoding changed decoded vectors.')
        result['fit']=report;candidates.append(result)
        print(json.dumps({'candidate':result['name'],'metrics':result['metrics'],'cost':result['whole_task_estimate']}),flush=True)
        write_json(HERE/'candidate-progress.json',candidates)
    best=min(candidates,key=lambda row:row['metrics']['error_2'])
    with np.load(LOCAL/(best['name']+'.npz'),allow_pickle=False) as data:
        bw=data['W'];bu=data['U'];bv=data['H2_I']-bw.sum(axis=1);best_words=data['angle_words'][model.n:]
    gauged,gauge_report=optimize_gauge(bw,bv)
    before=np.einsum('rbc,rbd->rcd',bw,projector_rows(bu)).reshape(-1,model.d)
    after=np.einsum('rbc,rbd->rcd',gauged,projector_rows(bu)).reshape(-1,model.d)
    invariance=float(np.linalg.norm(before.T@before-after.T@after))
    if invariance>1e-9:raise ArithmeticError('Copy gauge changed the pair tensor.')
    gauge_report['pair_tensor_drift']=invariance
    offset_gauge_drift=float(np.linalg.norm(np.einsum('rc,rcd->d',bv,
        after.reshape(model.r,model.c,model.d)-before.reshape(model.r,model.c,model.d))))
    if offset_gauge_drift>1e-10:raise ArithmeticError('Copy gauge changed the offset contraction.')
    gauge_report['offset_contraction_drift']=offset_gauge_drift
    gauged_result,_=official_candidate(helper,verifier,maths,inputs,model,bu,gauged,bv,
        best['precision_bits'],'global-gauged-beta'+str(best['precision_bits']),two_body_words=best_words)
    gauged_result['gauge_fit']=gauge_report;candidates.append(gauged_result)
    result={'status':'FEASIBILITY_ONLY_NO_CIRCUITS','shape':list(SHAPE),'provenance':provenance,
        'source_sha256':{name:sha(HERE/name) for name in ['refit.py','GramCores.wl','GlobalRefit.wl','symbolic.wls','symbolic.json','inputs.json']},
        'helper_sha256':sha(PROJECT/'dfthc_gram/study.py'),
        'versions':{'numpy':np.__version__,'scipy':scipy.__version__},
        'floating_fits':floating,'selected_floating_fit':selected,'candidates':candidates,
        'metric_passes':sum(row['metric_gates_pass'] for row in candidates),
        'baseline_cost':meta['baseline_cost'],'recorded_score_context':meta['recorded_score_context'],
        'new_competition_score':None,'circuit_implementation':False,'committed':False,'pushed':False}
    write_json(HERE/'results.json',result)
    print(json.dumps({'status':result['status'],'metric_passes':result['metric_passes']}),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('stage',choices=['prepare','run'])
    parser.add_argument('--data',type=Path,required=True);args=parser.parse_args()
    {'prepare':prepare,'run':run}[args.stage](args)
