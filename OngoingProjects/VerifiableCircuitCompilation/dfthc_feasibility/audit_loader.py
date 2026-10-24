"""Profile the released signed-gauge four-address loader in memory.

The external Apache-2.0 source is not vendored. This audit does not write a
submission or run the full circuit verifier. Counts are static projections.
"""

import argparse
import dataclasses
import hashlib
import importlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE_HASHES = {
    'gauge_qrom.py':'e7d8569a51cb563687f2c7593090ca605de834baa6282b878777ea7ec0b79466',
    'banked_lab.py':'a77e4819ba4f2b0b059847539b04b137179e82af5d0cea60da533e0012304e2a',
    'chemistry_refine.py':'83cc9eb7196a01fb729f1934a9fb38dd71b3b53a7331e334bd452e0b737f8539'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--words',type=Path,required=True)
    args = parser.parse_args()
    for name,digest in SOURCE_HASHES.items():
        if hashlib.sha256((args.source/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('Released source hash mismatch: '+name)
    sys.path[:0] = [str(args.data/'verifiers'),str(args.source)]
    verifier = importlib.import_module('dfthcsos.verify')
    maths = importlib.import_module('dfthcsos.scoreinfra.dfthc_math')
    gauge = importlib.import_module('gauge_qrom')
    banks = importlib.import_module('banked_lab')
    with np.load(args.data/'sample_submission/dfthcsos/submission.npz',allow_pickle=False) as z:
        params = {k:z[k].copy() for k in z.files}
    inputs = verifier.load_inputs(params)
    if hashlib.sha256(args.words.read_bytes()).hexdigest() != '9afd0ca5e90dd1740fe713e4ac1c26736184ea9bc4d4b57a6eeefafd608ee193':
        raise ValueError('Public word table hash mismatch.')
    with np.load(args.words,allow_pickle=False) as z:
        words,w,h2i = z['words'].astype(np.int64),z['W'],z['H2_I']
    data = {**params,'angle_words':words,'W':w,'H2_I':h2i,
            'U':maths.reconstruct_u_vectors(words[76:]/512,params['tree_edges']).reshape(65,60,76)}
    signed,evidence = gauge.canonicalize(SimpleNamespace(math=maths),data,10)
    candidate = dataclasses.replace(inputs,params=signed,w=w,W=h2i-w.sum(axis=1))
    metrics = verifier.evaluate_words(candidate,signed['angle_words'],signed['tree_edges'],10)
    if not verifier.passes_accuracy(metrics):
        raise ValueError('Released signed-gauge words failed official metric gates.')
    built = [gauge.build(delta,return_report=True) for delta in banks.batch_deltas(signed['angle_words'],3)]
    circuits,profiles = zip(*built)
    resources = verifier.resource_costs(candidate,circuits,10,2*sum(p['static_toffoli'] for p in profiles))
    projected = verifier.score_from_metrics(metrics,resources)
    incumbent = json.loads((HERE/'inputs.json').read_text())['incumbent_cost']
    report = {'status':'RELEASED_SOURCE_AUDIT_STATIC_RESOURCE_PROJECTION',
              'source_sha256':SOURCE_HASHES,'gauge_evidence':evidence,
              'metrics':metrics,'metric_gates_pass':True,'resource_projection':resources,
              'conditional_cost':projected,'conditional_improvement_over_incumbent':incumbent/projected,
              'circuit_profiles':profiles,'circuit_count':len(circuits),
              'full_circuit_verification':'NOT_RUN','submission_written':False,
              'scope':'No column-gauge optimization or factor training; public factors fixed.'}
    (HERE/'released-loader.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='circuit_profiles'},indent=2),flush=True)


if __name__ == '__main__':
    main()
