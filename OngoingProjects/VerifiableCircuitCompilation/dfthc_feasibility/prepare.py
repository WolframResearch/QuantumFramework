"""Prepare the official shifted pair matrix; audit released numeric inputs.

All large arrays are kept under ignored local/. No optimization is performed.
"""

import argparse
import dataclasses
import hashlib
import importlib
import json
from pathlib import Path
import sys

import numpy as np
import scipy
import scipy.linalg

HERE = Path(__file__).resolve().parent
SHAPES = [(65,60,60),(65,60,30),(15,57,19),(16,57,19),
          (20,48,24),(24,48,24),(32,32,16),(32,48,24),
          (40,48,32),(48,48,32),(56,48,40),(65,48,30)]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pack_matrix(tensor):
    n = tensor.shape[0]
    i,j = np.tril_indices(n)
    scale = np.where(i == j, 1., np.sqrt(2.))
    return tensor[i[:,None],j[:,None],i[None,:],j[None,:]] * scale[:,None] * scale[None,:]


def resource_row(verifier, shape, beta, batch, lambda_eff):
    r,b,c = shape
    addresses = 76+r*b
    selection = int(np.ceil(np.log2(addresses)))
    angle_temp = selection-1
    angle_q = selection+batch*beta+angle_temp
    h = verifier.heuristic_costs(76,r,b,c,beta)
    angle_t = 2*int(np.ceil(75/batch))*(addresses-2)
    total_t = angle_t+h['toffoli']
    total_q = h['persistent_qubits']+batch*beta+max(h['temporary_qubits'],angle_temp)
    return {'shape':list(shape),'beta':beta,'batch':batch,'addresses':addresses,
            'angle_qubits':angle_q,'angle_toffoli':angle_t,
            'heuristic':h,'toffoli':total_t,'qubits':total_q,
            'conditional_cost':lambda_eff*total_t*total_q}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--released',type=Path,required=True)
    args = parser.parse_args()
    sys.path.insert(0,str(args.data/'verifiers'))
    verifier = importlib.import_module('dfthcsos.verify')
    math = importlib.import_module('dfthcsos.scoreinfra.dfthc_math')
    source = args.data/'sample_submission/dfthcsos/submission.npz'
    with np.load(source,allow_pickle=False) as z:
        params = {k:z[k].copy() for k in z.files}
    inputs = verifier.load_inputs(params)
    n = inputs.num_orbitals
    target = np.array(inputs.ref['tei'],dtype=float,copy=True)
    identity = np.eye(n)
    target += (np.einsum('pq,rs->pqrs',inputs.S,identity)
               +np.einsum('pq,rs->pqrs',identity,inputs.S))/2
    symmetry = {'within_pair':float(np.max(np.abs(target-target.swapaxes(0,1)))),
                'pair_exchange':float(np.max(np.abs(target-target.transpose(2,3,0,1))))}
    if max(symmetry.values()) != 0:
        raise ValueError('Target lacks exact pair symmetries; revise packing before screening.')
    pair = pack_matrix(target)
    local = HERE/'local';local.mkdir(exist_ok=True)
    pair.astype('<f8').tofile(local/'shifted-pairs.f64')
    vectors = math.reconstruct_u_vectors(params['angle_words']/2**11,params['tree_edges'])
    x = math.reconstruct_x(vectors[n:].reshape(65,60,n),inputs.w)
    reconstructed = math.reconstruct_h2(x)
    full_residual = float(np.linalg.norm(reconstructed-target))
    packed_residual = float(np.linalg.norm(pack_matrix(reconstructed)-pair))
    del reconstructed,x,target
    incumbent = json.loads((HERE.parent/'other_tasks/mixing-results.json').read_text())['full_official_best_result']
    if abs(packed_residual-full_residual)>1e-10:
        raise ValueError('Pair packing did not preserve the official full Frobenius residual.')
    # Independent SciPy driver. WL computes its own spectrum from the binary matrix.
    values = scipy.linalg.eigvalsh(pair,driver='evd',check_finite=True)
    np.save(local/'scipy-eigenvalues.npy',values)
    (local/'scipy-eigenvalues.json').write_text(json.dumps(values.tolist()))
    audits = []
    research = args.released/'starter_kit/baselines/FeMoCo76o113e_65_60_60_ep50000_seed0.npz'
    with np.load(research,allow_pickle=False) as z:
        keys = list(z.files)
        research_params = {k:z[k] for k in ['U','W','H2_I','H1_I','Hsym_lower'] if k in z}
    audits.append({'name':'official_research_baseline','sha256':sha(research),
                   'keys':keys,'shape':list(research_params['W'].shape),
                   'matches_sample_parameters':{k:bool(np.array_equal(v,params[k])) for k,v in research_params.items()},
                   'note':'Stored training diagnostics are not competition scores.'})
    public = args.released/'femoco_public_words.npz'
    with np.load(public,allow_pickle=False) as z:
        words,w,h2i = z['words'],z['W'],z['H2_I']
    public_vectors = math.reconstruct_u_vectors(words/2**9,params['tree_edges'])
    updated = {**params,'U':public_vectors[n:].reshape(65,60,n),'W':w,'H2_I':h2i}
    candidate = dataclasses.replace(inputs,params=updated,w=w,W=h2i-w.sum(axis=1))
    metrics = verifier.evaluate_words(candidate,words,params['tree_edges'],10)
    audits.append({'name':'public_C30_10bit','sha256':sha(public),'shape':list(w.shape),
                   'word_shape':list(words.shape),'word_range':[int(words.min()),int(words.max())],
                   'metrics':{k:(v if not isinstance(v,float) or np.isfinite(v) else None) for k,v in metrics.items()},
                   'official_metric_gates_pass':bool(verifier.passes_accuracy(metrics)),
                   'full_circuit_verification':'NOT_RUN_AUDIT_ONLY'})
    rows = [resource_row(verifier,shape,beta,j,incumbent['lambda_eff'])
            for shape in SHAPES for beta in [10,11,12,14,16] for j in range(1,76)
            if 2*int(np.ceil(np.log2(76+shape[0]*shape[1])))-1+j*beta<=128]
    provenance = {str(p.relative_to(args.data)):sha(p) for p in [source,
        args.data/'verifiers/dfthcsos/FeMoCo76o113e.fcidump',
        args.data/'verifiers/dfthcsos/FeMoCo76o113e.fixed_shifts.npz',
        args.data/'verifiers/dfthcsos/verify.py',
        args.data/'verifiers/dfthcsos/scoreinfra/costs_prepare.py',
        args.data/'verifiers/dfthcsos/scoreinfra/dfthc_math.py']}
    report = {'status':'AUDIT_AND_FEASIBILITY_ONLY','versions':{'numpy':np.__version__,'scipy':scipy.__version__},
              'input_sha256':provenance,'matrix_sha256':sha(local/'shifted-pairs.f64'),
              'pair_dimension':len(pair),'symmetry_max_error':symmetry,
              'target_frobenius':float(np.linalg.norm(pair)),
              'baseline_full_residual':full_residual,'baseline_packed_residual':packed_residual,
              'incumbent_cost':incumbent['score'],'reference_lambda_eff':incumbent['lambda_eff'],
              'shapes':[list(s) for s in SHAPES], 'released_numeric_audits':audits,'cost_rows':rows}
    (HERE/'inputs.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['cost_rows','input_sha256']},indent=2),flush=True)


if __name__ == '__main__':
    main()
