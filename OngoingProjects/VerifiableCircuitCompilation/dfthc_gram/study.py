"""Pinned DFTHC Gram-core study: numerical feasibility, never circuit building."""
from __future__ import annotations
import argparse
import dataclasses
import hashlib
import importlib
import json
from pathlib import Path
import sys

import numpy as np
import scipy
import scipy.linalg as la
import scipy.optimize as opt

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent
LOCAL=HERE/'local'
PRIME=2147483647
GUARD=1e-8
SCORE_REPORTS={
    'invsqrt':('results/baseline-verified.json','compensated_quadratic/verified.json'),
    'DQI':('other_tasks/dqi-baseline.json','dqi_gf8/verified.json'),
    'DFTHC':('other_tasks/dfthcsos-baseline.json','dfthc_c30/verified.json')}
SHAPES=[(65,60,30),(65,59,30),(65,56,30),(65,48,30),(65,32,30),
        (65,60,29),(65,60,24),(65,60,16),(32,32,16),(20,48,24),
        (16,32,20),(10,32,30),(8,32,32)]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def score_context():
    """Derive the unchanged context from recorded PASS reports, without reruns."""
    tracks={};reports={}
    for track,names in SCORE_REPORTS.items():
        rows=[]
        for name in names:
            row=json.loads((PROJECT/name).read_text())
            if row['status']!='PASS' or not np.isfinite(row['score']) or row['score']<=0:
                raise ValueError('Invalid recorded score report: '+name)
            rows.append(row);reports[name]=sha(PROJECT/name)
        baseline,current=[row['score'] for row in rows]
        tracks[track]={'baseline_cost':baseline,'current_cost':current,
                       'improvement':baseline/current,'report_paths':list(names)}
    combined=float(np.prod([row['improvement'] for row in tracks.values()])**(1/3))
    return {'report_sha256':reports,'tracks':tracks,'combined_local_estimate':combined,
            'formula':'geometric mean of the three baseline/current cost ratios',
            'other_tracks_reverified_in_this_study':False,'new_competition_score':False}

def load(data):
    pinned=json.loads((PROJECT/'dfthc_feasibility/inputs.json').read_text())
    for name,digest in pinned['input_sha256'].items():
        if sha(data/name)!=digest:raise ValueError('Official dependency changed: '+name)
    sys.path.insert(0,str(data/'verifiers'))
    verifier=importlib.import_module('dfthcsos.verify')
    maths=importlib.import_module('dfthcsos.scoreinfra.dfthc_math')
    source=PROJECT/'dfthc_c30/local/submission/submission.npz'
    build=json.loads((PROJECT/'dfthc_c30/build.json').read_text())
    if sha(source)!=build['factorization_sha256']:raise ValueError('C30 factorization changed.')
    with np.load(source,allow_pickle=False) as archive:
        params={k:archive[k].copy() for k in archive.files}
    inputs=verifier.load_inputs(params)
    inputs=dataclasses.replace(inputs,params={**inputs.params,
        'angle_words':params['angle_words'],'tree_edges':params['tree_edges']})
    vectors=maths.reconstruct_u_vectors(params['angle_words']/512,params['tree_edges'])
    u=vectors[76:].reshape(65,60,76)
    u=u/np.linalg.norm(u,axis=2,keepdims=True)
    pair_path=PROJECT/'dfthc_feasibility/local/shifted-pairs.f64'
    if sha(pair_path)!=pinned['matrix_sha256']:raise ValueError('Shifted target changed.')
    pair=np.fromfile(pair_path,dtype='<f8').reshape(2926,2926)
    provenance={'official':pinned['input_sha256'],'c30_factorization':sha(source),
        'shifted_pair_matrix':sha(pair_path),
        'incumbent_report':sha(PROJECT/'dfthc_c30/verified.json'),
        'recorded_score_reports':score_context()['report_sha256']}
    return verifier,maths,inputs,u,pair,provenance

def gram(u,v):return (u@v.swapaxes(-1,-2))**2

def weights_for(u,u0,w0):
    g=gram(u,u);rhs=gram(u,u0)@w0
    # Independent per-block dense solves; singular dictionaries are rejected.
    values=np.linalg.eigvalsh(g)
    if float(values.min())<1e-12:raise ValueError('Singular or unresolved dictionary Gram.')
    return np.linalg.solve(g,rhs)

def factor_loss(u,w,u0,w0):
    own=np.sum((w@w.swapaxes(-1,-2))*gram(u,u),axis=(-2,-1))
    cross=np.sum((w@w0.swapaxes(-1,-2))*gram(u,u0),axis=(-2,-1))
    target=np.sum((w0@w0.swapaxes(-1,-2))*gram(u0,u0),axis=(-2,-1))
    return own-2*cross+target

def gradient(u,w,u0,w0):
    s=u@u.swapaxes(-1,-2);t=u@u0.swapaxes(-1,-2)
    return 4*((w@w.swapaxes(-1,-2)*s)@u-(w@w0.swapaxes(-1,-2)*t)@u0)

def greedy_dictionary(u0,w0,count):
    if type(count)is not int or not 1<=count<=u0.shape[1]:raise ValueError('Dictionary count.')
    indices=[]
    for u,w in zip(u0,w0):
        residual=gram(u,u).copy();chosen=[]
        for _ in range(count):
            diagonal=np.diag(residual)
            gain=np.sum((residual@w)**2,axis=1)/np.maximum(diagonal,1e-30)
            gain[chosen]=-np.inf
            index=int(gain.argmax());pivot=float(residual[index,index])
            if pivot<1e-12:raise ValueError('No independent projector pivot.')
            chosen.append(index)
            residual-=np.outer(residual[:,index],residual[index,:])/pivot
            residual=(residual+residual.T)/2
        indices.append(chosen)
    selected=u0[np.arange(len(u0))[:,None],np.array(indices)].copy()
    return selected,np.array(indices)

def optimize_dictionary(u,u0,w0,steps):
    w=weights_for(u,u0,w0);initial=factor_loss(u,w,u0,w0)
    history=[float(initial.sum())];accepted=np.zeros(len(u),dtype=int)
    for _ in range(steps):
        old=factor_loss(u,w,u0,w0);g=gradient(u,w,u0,w0)
        tangent=g-np.sum(g*u,axis=2,keepdims=True)*u
        norm2=np.sum(tangent**2,axis=(1,2))
        rate=1/(1+np.sqrt(norm2));pending=norm2>1e-28
        for _ in range(16):
            trial=u-rate[:,None,None]*tangent
            trial/=np.linalg.norm(trial,axis=2,keepdims=True)
            tw=weights_for(trial,u0,w0);loss=factor_loss(trial,tw,u0,w0)
            take=pending&(loss<=old-1e-4*rate*norm2+1e-14)
            u[take]=trial[take];w[take]=tw[take];accepted[take]+=1
            pending&=~take
            if not pending.any():break
            rate[pending]*=.5
        history.append(float(factor_loss(u,w,u0,w0).sum()))
    if np.any(np.diff(history)>1e-9):raise ArithmeticError('Variable projection increased loss.')
    return u,w,{'steps':steps,'accepted_block_steps':int(accepted.sum()),
        'initial_squared_factor_error':float(initial.sum()),
        'final_squared_factor_error':history[-1],'loss_history':history,
        'global_optimality_claim':False}

def encode(vectors,tree,beta=10):
    reduced=np.array(vectors,copy=True,dtype=float)
    reduced/=np.linalg.norm(reduced,axis=1,keepdims=True)
    angles=np.empty((len(reduced),75))
    for j,(p,c) in reversed(list(enumerate(tree))):
        radius=np.hypot(reduced[:,p],reduced[:,c])
        sign=np.where(np.abs(reduced[:,p])<1e-12,np.where(reduced[:,c]>0,-1,1),np.where(reduced[:,p]<0,-1,1))
        angles[:,j]=np.where(radius<1e-12,0,np.arctan2(reduced[:,c]*sign,reduced[:,p]*sign))
        reduced[:,p]=radius*sign;reduced[:,c]=0
    return np.rint(angles*((1<<(beta-1))/np.pi)).astype(np.int64)%(1<<beta)

def evaluate(verifier,maths,inputs,u,w,v,tree,name,u0=None,w0=None,fixed_words=None):
    r,b,c=w.shape
    two_words=encode(u.reshape(r*b,76),tree) if fixed_words is None else fixed_words[76:]
    decoded=maths.reconstruct_u_vectors(two_words/512,tree).reshape(r,b,76)
    if u0 is not None:w=weights_for(decoded,u0,w0)
    params={**inputs.params,'U':decoded,'W':w,'H2_I':v+w.sum(axis=1)}
    native_w,native_v,native_s,native_beta=verifier.native_inputs(params)
    if float(np.max(np.abs(native_v-v)))>1e-11:raise ArithmeticError('Native-offset serialization drift.')
    w,v=native_w,native_v
    candidate=dataclasses.replace(inputs,params=params,w=w,W=v,S=native_s,beta_1=native_beta)
    x=maths.reconstruct_x(decoded,w)
    h1=maths.build_h1_prime(inputs.ref['h1_thc'],inputs.S,inputs.beta_1,
                          inputs.ref['electrons'],v,x)
    if fixed_words is None:
        _,columns=np.linalg.eigh(h1);one_words=encode(columns.T,tree)
    else:one_words=fixed_words[:76]
    words=np.concatenate([one_words,two_words])
    metrics=verifier.evaluate_words(candidate,words,tree,10)
    if not np.isfinite(list(metrics.values())).all():raise ArithmeticError('Nonfinite metrics.')
    np.savez(LOCAL/(name+'.npz'),u=decoded,w=w,v=v,words=words)
    return {'name':name,'shape':[r,b,c],'metrics':metrics,
        'metric_gates_pass':bool(verifier.passes_accuracy(metrics)),
        'precision':10,'one_body_rows_recomputed':fixed_words is None,
        'all_angle_words_unchanged':fixed_words is not None and bool(np.array_equal(words,fixed_words)),
        'weights_refitted_after_quantization':u0 is not None,
        'private_candidate_sha256':sha(LOCAL/(name+'.npz')),
        'verification_scope':'OFFICIAL_DECODED_WORD_METRICS_ONLY_NO_CIRCUITS'},decoded,w

def resources(verifier,shape,lam,beta=10):
    r,b,c=shape;a=76+r*b;bits=(a-1).bit_length()
    h=verifier.heuristic_costs(76,r,b,c,beta);rows=[]
    for j in range(1,76):
        aq=2*bits-1+j*beta
        if aq>128:continue
        at=2*((75+j-1)//j)*(a-2)
        t=at+h['toffoli'];q=h['persistent_qubits']+j*beta+max(h['temporary_qubits'],bits-1)
        rows.append({'shape':list(shape),'beta':beta,'batch':j,'angle_qubits':aq,
            'angle_toffoli':at,'toffoli':t,'qubits':q,'heuristic':h,
            'assumed_lambda_eff':lam,'conditional_cost':lam*t*q})
    return min(rows,key=lambda row:row['conditional_cost'])

def coupled_value_gradient(raw,cross,a,scale):
    r,c=cross.shape[:2]
    z=raw.reshape(r,c);norm=np.linalg.norm(z,axis=1,keepdims=True)
    if float(norm.min())<1e-14:raise ArithmeticError('Zero removal direction.')
    z=z/norm
    products=np.einsum('rcsd,sd->rcs',cross,z,optimize=True)
    gy=np.einsum('rc,rcs->rs',z,products,optimize=True)
    az=np.einsum('rcd,rd->rc',a,z,optimize=True)
    value=float(np.sum(gy**2)-2*np.sum(z*az))
    grad=4*(np.einsum('rs,rcs->rc',gy,products,optimize=True)-az)
    grad=(grad-z*np.sum(z*grad,axis=1,keepdims=True))/norm
    return value/scale,grad.ravel()/scale

def offset_preserving_compression(u,w,v,pair):
    """Coupled full-target rank-one removals in the offset-orthogonal spheres."""
    r,b,c=w.shape;i,j=np.tril_indices(76);scale=np.where(i==j,1.,np.sqrt(2.))
    p=u[:,:,i]*u[:,:,j]*scale
    x=np.einsum('rbc,rbp->rcp',w,p,optimize=True)
    flat=x.reshape(r*c,-1);residual=flat.T@flat-pair
    baseline_squared=float(np.sum(residual**2))
    basis=np.stack([la.null_space(row[None,:]) for row in v])
    if basis.shape!=(r,c,c-1):raise ValueError('Offset nullspace dimension unresolved.')
    xs=np.einsum('rck,rcp->rkp',basis,x,optimize=True)
    xsflat=xs.reshape(r*(c-1),-1)
    cross=(xsflat@xsflat.T).reshape(r,c-1,r,c-1)
    a=np.einsum('rcp,rdp->rcd',(xsflat@residual).reshape(xs.shape),xs,optimize=True)
    covariance=xs@xs.swapaxes(-1,-2)
    _,eigenvectors=np.linalg.eigh(covariance);initial=eigenvectors[:,:,0]
    def value_gradient(raw):return coupled_value_gradient(raw,cross,a,baseline_squared)
    before=value_gradient(initial.ravel())[0]
    fit=opt.minimize(value_gradient,initial.ravel(),method='L-BFGS-B',jac=True,
        options={'maxiter':200,'maxls':30,'ftol':1e-14,'gtol':1e-9})
    z=fit.x.reshape(r,c-1);z/=np.linalg.norm(z,axis=1,keepdims=True)
    directions=np.einsum('rck,rk->rc',basis,z)
    complements=[]
    for direction in directions:
        pivot=int(np.abs(direction).argmax())
        direction=direction*np.sign(direction[pivot]);axis=np.eye(c)[pivot]
        skew=np.outer(direction,axis)-np.outer(axis,direction)
        rotation=np.eye(c)+skew+skew@skew/(1+direction[pivot])
        complements.append(np.delete(rotation,pivot,axis=1))
    q=np.stack(complements)
    # K = [q q^T/tau + S, q; q^T, tau]. Only S is rescaled.
    # The resulting two-body tensor is affine in 65 nonnegative scales.
    vh=v/np.linalg.norm(v,axis=1,keepdims=True)
    parallel=np.einsum('rc,rcp->rp',vh,x,optimize=True)
    cov=np.einsum('rck,rcp->rkp',q,x,optimize=True)-np.einsum(
        'rk,rp->rkp',np.einsum('rck,rc->rk',q,vh),parallel)
    cf=cov.reshape(r*(c-1),-1)
    gg=(cf@cf.T).reshape(r,c-1,r,c-1)
    normal=np.sum(gg**2,axis=(1,3))
    target=pair-parallel.T@parallel
    rhs=np.einsum('rcp,rcp->r',(cf@target).reshape(cov.shape),cov,optimize=True)
    upper=la.cholesky((normal+normal.T)/2,lower=False)
    target_ls=la.solve_triangular(upper.T,rhs,lower=True)
    scales,nnls_residual=opt.nnls(upper,target_ls,maxiter=500)
    pv=vh[:,:,None]*vh[:,None,:]
    transform=pv+np.sqrt(scales)[:,None,None]*(np.eye(c)-pv)
    nw=w@transform@q;nv=np.einsum('rck,rc->rk',q,v)
    after_x=np.einsum('rbc,rbp->rcp',nw,p,optimize=True).reshape(r*(c-1),-1)
    after_error=float(np.linalg.norm(after_x.T@after_x-pair))
    kkt=normal@scales-rhs
    kkt_error=max(float(np.max(np.abs(kkt[scales>1e-10]))),
                  float(np.max(np.maximum(-kkt[scales<=1e-10],0),initial=0)))
    offset_drift=float(np.max(np.abs(w@v[:,:,None]-nw@nv[:,:,None])))
    norm_drift=float(np.max(np.abs(np.sum(v*v,axis=1)-np.sum(nv*nv,axis=1))))
    if max(offset_drift,norm_drift)>1e-10:raise ArithmeticError('Offset-preservation invariant failed.')
    return nw,nv,{'iterations':int(fit.nit),'optimizer_success':bool(fit.success),
        'initial_predicted_error':float(np.sqrt(baseline_squared*(1+before))),
        'before_schur_rescale_error':float(np.sqrt(baseline_squared*(1+fit.fun))),
        'final_predicted_error':after_error,
        'schur_scale_range':[float(scales.min()),float(scales.max())],
        'schur_nnls_kkt_error':kkt_error,'schur_nnls_residual':float(nnls_residual),
        'baseline_pair_residual':float(np.sqrt(baseline_squared)),
        'max_offset_coefficient_drift':offset_drift,'max_offset_norm_drift':norm_drift,
        'global_optimality_claim':False}

def span_bound(u,pair,rank_cap):
    rows=u.reshape(-1,76);i,j=np.tril_indices(76);scale=np.where(i==j,1.,np.sqrt(2.))
    p=(rows[:,i]*rows[:,j]*scale).T
    q,r=la.qr(p,mode='economic',check_finite=True)
    # Pivot-free QR is used only where its diagonal is numerically resolved.
    if len(np.diag(r))>2926 or np.min(np.abs(np.diag(r)))<1e-10:
        raise ValueError('Dictionary span requires rank refinement.')
    compressed=q.T@pair@q
    ev,vec=la.eigh((compressed+compressed.T)/2,check_finite=True)
    order=np.argsort(ev)[::-1];keep=order[:min(rank_cap,len(ev))]
    keep=keep[ev[keep]>0]
    basis=q@vec[:,keep]
    fitted=(basis*ev[keep])@basis.T
    error=float(np.linalg.norm(pair-fitted))
    return {'dictionary_projectors':len(rows),'span_dimension':q.shape[1],
        'relaxed_rank_cap':rank_cap,'error_lower_bound':error,
        'numerical_guard':GUARD,'orthogonality_residual':float(np.linalg.norm(q.T@q-np.eye(q.shape[1]))),
        'screen':'REJECT_THIS_FIXED_DICTIONARY' if error-GUARD>.020 else 'NOT_RULED_OUT',
        'scope':'Global target; fixed projector directions; all block constraints relaxed. Not a bound on reoptimized dictionaries.'}

def prepare(args):
    verifier,maths,inputs,u,pair,provenance=load(args.data)
    LOCAL.mkdir(exist_ok=True)
    baseline=verifier.evaluate_words(inputs,inputs.params['angle_words'],inputs.params['tree_edges'],10)
    incumbent=json.loads((PROJECT/'dfthc_c30/verified.json').read_text())
    for k,v in baseline.items():
        if isinstance(v,bool):assert v==incumbent['metrics'][k]
        else:assert abs(v-incumbent['metrics'][k])<1e-10,k
    # Numerically instantiate the Gram representation, retaining raw arrays privately.
    write_json(LOCAL/'arrays.json',{'u':u.tolist(),'w':inputs.w.tolist(),'v':inputs.W.tolist()})
    meta={'status':'PINNED_INPUTS','provenance':provenance,'private_arrays_sha256':sha(LOCAL/'arrays.json'),
        'prepare_source_sha256':sha(Path(__file__)),'shape':[65,60,30],
        'baseline_metrics':baseline,'prime':PRIME,'target_frobenius':float(np.linalg.norm(pair)),
        'baseline_cost':incumbent['score'],'baseline_lambda_eff':incumbent['lambda_eff']}
    write_json(HERE/'inputs.json',meta)
    print(json.dumps({k:v for k,v in meta.items() if k not in ['provenance','baseline_metrics']},indent=2),flush=True)

def run(args):
    verifier,maths,inputs,u0,pair,provenance=load(args.data)
    meta=json.loads((HERE/'inputs.json').read_text())
    if meta['provenance']!=provenance or meta['prepare_source_sha256']!=sha(Path(__file__)):
        raise ValueError('Stale input audit.')
    wl=json.loads((HERE/'symbolic.json').read_text())
    if wl['status']!='PASS' or not all(wl['symbolic']['checks'].values()):raise ValueError('WL gate failed.')
    g0=gram(u0,u0);w0=inputs.w;v0=inputs.W;tree=inputs.params['tree_edges']
    candidates=[];private={}
    for b in [59,56,48,32]:
        u,indices=greedy_dictionary(u0,w0,b)
        if b in [48,32]:u,w,fit=optimize_dictionary(u,u0,w0,50)
        else:w=weights_for(u,u0,w0);fit={'steps':0,'squared_factor_error':float(factor_loss(u,w,u0,w0).sum())}
        name='dictionary-B'+str(b)
        result,decoded,weights=evaluate(verifier,maths,inputs,u,w,v0,tree,name,u0,w0)
        result['fit']=fit;result['initial_dictionary_indices_sha256']=hashlib.sha256(indices.tobytes()).hexdigest()
        candidates.append(result);private[name]=(decoded,weights)
        print(json.dumps({'candidate':name,'metrics':result['metrics'],'metric_gates_pass':result['metric_gates_pass']}),flush=True)
    covariance=w0.swapaxes(1,2)@g0@w0
    values,vectors=np.linalg.eigh(covariance)
    for c in [29,24,16]:
        rotation=vectors[:,:,::-1][:,:,:c]
        w=w0@rotation;v=np.einsum('rkc,rk->rc',rotation,v0)
        result,_,_=evaluate(verifier,maths,inputs,u0,w,v,tree,'copies-C'+str(c))
        candidates.append(result)
        print(json.dumps({'candidate':result['name'],'metrics':result['metrics'],'metric_gates_pass':result['metric_gates_pass']}),flush=True)
    nw,nv,fit=offset_preserving_compression(u0,w0,v0,pair)
    result,_,_=evaluate(verifier,maths,inputs,u0,nw,nv,tree,'offset-preserving-C29',
        fixed_words=inputs.params['angle_words'])
    result['fit']=fit
    actual_error=result['metrics']['error_2']
    if abs(actual_error-fit['final_predicted_error'])>1e-9:raise ArithmeticError('Coupled residual model disagrees with official metrics.')
    candidates.append(result)
    print(json.dumps({'candidate':result['name'],'metrics':result['metrics'],'metric_gates_pass':result['metric_gates_pass'],'fit':fit}),flush=True)
    bounds={}
    for name in ['dictionary-B32']:
        bounds[name]=span_bound(private[name][0],pair,65*30)
        print(json.dumps({'fixed_dictionary_bound':bounds[name]}),flush=True)
    eig=la.eigvalsh(pair)
    def tail(k):return float(np.sqrt(np.sum(eig[eig<0]**2)+np.sum(np.sort(eig[eig>0])[::-1][k:]**2)))
    shapes=[]
    for shape in SHAPES:
        rank=min(2926,shape[0]*min(shape[1:]));lower=tail(rank)
        shapes.append({'shape':list(shape),'rank_cap':rank,'psd_error_lower_bound':lower,
            'rank_screen':'REJECT' if lower-GUARD>.020 else 'NOT_RULED_OUT',
            'reference_loader':resources(verifier,shape,meta['baseline_lambda_eff']),
            'normalization_is_assumed_not_fitted':True})
    for candidate in candidates:
        candidate['reference_loader']=resources(verifier,candidate['shape'],candidate['metrics']['lambda_eff'])
        if candidate['all_angle_words_unchanged']:
            old=json.loads((PROJECT/'dfthc_c30/verified.json').read_text())
            h=verifier.heuristic_costs(76,*candidate['shape'],10)
            tt=old['resources']['angles_toffoli']+h['toffoli']
            qq=h['persistent_qubits']+30+max(h['temporary_qubits'],old['resources']['angles_temporary_qubits'])
            candidate['unchanged_loader_estimate']={'toffoli':tt,'qubits':qq,'angle_qubits':128,
                'cost':candidate['metrics']['lambda_eff']*tt*qq,
                'scope':'Existing verified ANGLE family unchanged; new factors have word metrics only, no new full circuit verification.'}
        candidate['valid_competition_score']=False
    context=score_context()
    if context['report_sha256']!=provenance['recorded_score_reports']:
        raise ValueError('Recorded score reports changed during the study.')
    out={'status':'FEASIBILITY_ONLY_NO_CIRCUITS','source_sha256':{n:sha(HERE/n) for n in
         ['study.py','GramCores.wl','run.wls','inputs.json','symbolic.json']},
         'provenance':provenance,'versions':{'numpy':np.__version__,'scipy':scipy.__version__},
         'baseline_cost':meta['baseline_cost'],'baseline_lambda_eff':meta['baseline_lambda_eff'],
         'candidates':candidates,'fixed_dictionary_bounds':bounds,'shape_estimates':shapes,
         'feasible_word_metric_candidates':sum(row['metric_gates_pass'] for row in candidates),
         'numerical_guard':GUARD,'recorded_score_context':context,
         'current_combined_local_estimate':context['combined_local_estimate'],
         'new_competition_score':None,'circuit_implementation':False,'pushed':False}
    write_json(HERE/'results.json',out)
    print(json.dumps({'status':out['status'],'candidates':len(candidates),'metric_passes':out['feasible_word_metric_candidates']}),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('stage',choices=['prepare','run'])
    parser.add_argument('--data',type=Path,required=True);args=parser.parse_args()
    {'prepare':prepare,'run':run}[args.stage](args)
