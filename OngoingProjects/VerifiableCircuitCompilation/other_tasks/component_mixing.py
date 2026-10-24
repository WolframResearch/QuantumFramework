"""Bounded orthogonal component-mixing feasibility; no angle circuit changes."""

import argparse
import dataclasses
import hashlib
import importlib
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

HERE=Path(__file__).resolve().parent


def pair_objective(a,b,theta):
    c,s=np.cos(theta),np.sin(theta)
    return float(np.abs(c*a+s*b).sum()**2+np.abs(-s*a+c*b).sum()**2)


def best_pair(a,b):
    """All sign boundaries and stationary points in one pi/2 period.

    Floating evaluation instantiates the WL-derived piecewise formula; it is
    not a formal floating-point optimality certificate.
    """
    a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
    if a.ndim!=1 or a.shape!=b.shape or not len(a) or not np.isfinite([a,b]).all():
        raise ValueError("Expected nonempty finite vectors with equal shapes.")
    period=math.pi/2
    active=(a!=0)|(b!=0)
    boundaries=np.unique(np.r_[0.,np.mod(np.arctan2(b[active],a[active]),period),period])
    candidates=list(boundaries)
    for left,right in zip(boundaries[:-1],boundaries[1:]):
        mid=(left+right)/2;c,s=math.cos(mid),math.sin(mid)
        sx,sy=np.sign(c*a+s*b),np.sign(-s*a+c*b)
        aa,bb,cc,dd=sx@a,sx@b,sy@a,sy@b
        cosine=(aa*aa+dd*dd-bb*bb-cc*cc)/2
        sine=aa*bb-cc*dd
        stationary=math.atan2(sine,cosine)/2
        candidates.extend(stationary+k*period for k in range(-2,3)
                          if left<stationary+k*period<right)
    values=np.array([pair_objective(a,b,t) for t in candidates])
    k=int(values.argmin())
    return float(candidates[k]),float(values[k])


def mix_pairs(augmented,rounds=2,pairs_per_group=60):
    mixed=augmented.copy();accepted=0;attempted=0
    for _ in range(rounds):
        for group in mixed:
            norms=np.abs(group).sum(axis=0)
            products=np.sign(group).T@group
            derivative=2*(norms[:,None]*products-norms[None,:]*products.T)
            i,j=np.triu_indices(group.shape[1],1)
            order=np.argsort(-np.abs(derivative[i,j]),kind="stable")[:pairs_per_group]
            for k in order:
                p,q=int(i[k]),int(j[k]);a,b=group[:,p].copy(),group[:,q].copy()
                old=pair_objective(a,b,0);theta,new=best_pair(a,b);attempted+=1
                if new<old-1e-12*max(1.,old):
                    c,s=math.cos(theta),math.sin(theta)
                    group[:,p],group[:,q]=c*a+s*b,-s*a+c*b
                    accepted+=1
    return mixed,{"rounds":rounds,"pairs_per_group":pairs_per_group,
                  "attempted_pair_rotations":attempted,"accepted_pair_rotations":accepted}


def spectral_mix(augmented):
    return np.array([group@np.linalg.svd(group,full_matrices=False)[2].T for group in augmented])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data",type=Path,required=True)
    args=parser.parse_args();started=time.monotonic()
    sys.path.insert(0,str(args.data/"verifiers"))
    verifier=importlib.import_module("dfthcsos.verify")
    source=args.data/"sample_submission/dfthcsos/submission.npz"
    with np.load(source,allow_pickle=False) as archive:
        params={k:archive[k].copy() for k in archive.files}
    inputs=verifier.load_inputs(params)
    words,edges=params["angle_words"],params["tree_edges"]
    brot=12
    base=verifier.evaluate_words(inputs,words,edges,brot)
    profile=json.loads((HERE/"dfthcsos-baseline.json").read_text())
    if profile["status"]!="PASS" or not math.isclose(base["lambda_eff"],profile["lambda_eff"],rel_tol=1e-12):
        raise ValueError("Baseline decoded words do not match the full verifier record.")
    augmented=np.concatenate([inputs.w,params["H2_I"][:,None,:]],axis=1)
    paired,search=mix_pairs(augmented)
    controls={"identity":augmented,"spectral_control":spectral_mix(augmented),"pair_coordinate":paired}
    records=[];best_name="identity";best_cost=profile["score"];best_params=params
    gram=np.einsum("ric,rjc->rij",augmented,augmented)
    for name,array in controls.items():
        updated={k:v.copy() for k,v in params.items()}
        updated["W"],updated["H2_I"]=array[:,:-1,:].copy(),array[:,-1,:].copy()
        native_W=updated["H2_I"]-updated["W"].sum(axis=1)
        candidate=dataclasses.replace(inputs,params=updated,w=updated["W"],W=native_W)
        metrics=verifier.evaluate_words(candidate,words,edges,brot)
        score=metrics["lambda_eff"]*profile["toffoli_count"]*profile["num_qubits"]
        passed=verifier.passes_accuracy(metrics)
        error=np.linalg.norm(np.einsum("ric,rjc->rij",array,array)-gram)/np.linalg.norm(gram)
        record={"name":name,"status":"OFFICIAL_METRICS_PASS" if passed else "OFFICIAL_METRICS_FAIL",
                "metrics":metrics,"score_from_unchanged_circuit_resources":score,
                "improvement":profile["score"]/score,"relative_augmented_gram_error":float(error)}
        records.append(record)
        if passed and score<best_cost:
            best_name,best_cost,best_params=name,score,updated
    local=HERE/"local";local.mkdir(exist_ok=True)
    np.savez(local/"mixed-feasibility.npz",**best_params)
    # Full official verification uses the unchanged supplied circuit family.
    full=verifier.verify(args.data/"sample_submission/dfthcsos",best_params,
                         tree_edges=edges,salt="qf-other-tasks-feasibility-v1")
    report={"status":"FEASIBILITY_ONLY","source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
            "analysis_source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "symbolic_sha256":hashlib.sha256((HERE/"symbolic.json").read_bytes()).hexdigest(),
            "search":search,"candidates":records,"best_name":best_name,
            "angle_words_unchanged":True,"factorization_dimensions_unchanged":True,
            "full_official_best_result":full,"best_improvement":profile["score"]/full["score"],
            "local_candidate":"local/mixed-feasibility.npz",
            "local_candidate_sha256":hashlib.sha256((local/"mixed-feasibility.npz").read_bytes()).hexdigest(),
            "elapsed_seconds":time.monotonic()-started}
    (HERE/"mixing-results.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2),flush=True)


if __name__=="__main__":
    main()
