"""Unmodified full official benchmark for the compensated quadratic circuit."""
import argparse
import json
from pathlib import Path
import platform
import time

import kickmix as km
import numpy as np

import generator as gen
from verify_quadratic import load_verifier, load_comparison, SALT, VERIFIER_SHA256

EVIDENCE=['CompensatedQuadratic.wl','run.wls','coefficients.json','generator.py',
          'screen.py','screening.json','build.json','verify.py']
COMPARISONS=['results/optimized-quadratic-verified.json','results/baseline-verified.json',
             'dqi_gf8/verified.json','dfthc_c30/verified.json',
             'other_tasks/dqi-baseline.json','other_tasks/dfthcsos-baseline.json']


def input_hashes():
    return {**{name:gen.digest(gen.HERE/name) for name in EVIDENCE},
            **{'../'+name:gen.digest(gen.PROJECT/name) for name in COMPARISONS+
               ['quadratic.py','optimized_quadratic.py','verify_quadratic.py',
                'feasibility/B-quadratic-192.json']}}


def check_inputs():
    gen.checked_rows()
    screen=json.loads((gen.HERE/'screening.json').read_text())
    build=json.loads((gen.HERE/'build.json').read_text())
    if (screen['generator_sha256']!=gen.digest(gen.HERE/'generator.py') or
            screen['screen_source_sha256']!=gen.digest(gen.HERE/'screen.py') or
            screen['coefficients_sha256']!=gen.digest(gen.HERE/'coefficients.json') or
            screen['verifier_sha256']!=VERIFIER_SHA256 or
            screen['standalone_default_matches_selected'] is not True or
            screen['model']['status']!='INTEGER_MODEL_PASS' or
            screen['model']['checked_inputs']!=3_133_452 or
            screen['model']['phase_failures']!=0 or
            screen['model']['partial_accumulator_violations']!=0 or
            screen['model']['minimum_margin_units']<=0 or
            len(screen['variants'])!=4 or
            any(not all(row['boundary'][key] is True for key in
                ['input_restored','scratch_clean','exact_phase_model_match','register_fit','priced_counts_match'])
                for row in screen['variants'])):
        raise ValueError('Integer screening or compiled boundary evidence is stale/failing.')
    if (build['generator_sha256']!=screen['generator_sha256'] or
            build['coefficients_sha256']!=screen['coefficients_sha256'] or
            build['circuit_sha256']!=screen['circuit_sha256'] or
            {key:build[key] for key in screen['selected_options']}!=screen['selected_options'] or
            any(gen.digest(gen.PROJECT/name)!=digest for name,digest in build['dependency_sha256'].items()) or
            gen.digest(gen.HERE/'local'/build['artifact'])!=build['circuit_sha256']):
        raise ValueError('Build, circuit or implementation dependency changed.')
    return screen,build


def compare(actual,baseline):
    read=lambda name:json.loads((gen.PROJECT/name).read_text())
    old=read('results/optimized-quadratic-verified.json')
    dqi,df=read('dqi_gf8/verified.json'),read('dfthc_c30/verified.json')
    dqi_base,df_base=read('other_tasks/dqi-baseline.json'),read('other_tasks/dfthcsos-baseline.json')
    if any(row['status']!='PASS' for row in [old,dqi,df,dqi_base,df_base]):
        raise ValueError('A comparison report is not verified.')
    if (old['verified_shots']!=200_540_928 or dqi['verified_shots']!=102_400 or
            df['verified_shots']!=254_464):
        raise ValueError('Unexpected comparison verification coverage.')
    for name,digest in old['input_sha256'].items():
        if gen.digest(gen.PROJECT/name)!=digest:
            raise ValueError('An incumbent input no longer matches its full report.')
    ratios={'invsqrt':baseline['score']/actual,'DQI':dqi_base['score']/dqi['score'],
            'DFTHC':df_base['score']/df['score']}
    combined=float(np.prod(list(ratios.values()))**(1/3))
    previous=(baseline['score']/old['score']*ratios['DQI']*ratios['DFTHC'])**(1/3)
    return dict(baseline_cost=baseline['score'],incumbent_cost=old['score'],
        baseline_improvement=baseline['score']/actual,incumbent_improvement=old['score']/actual,
        incumbent_cost_reduction_percent=100*(1-actual/old['score']),
        baseline_cost_reduction_percent=100*(1-actual/baseline['score']),
        track_improvements=ratios,combined_local_estimate=combined,
        previous_combined_local_estimate=previous,
        combined_improvement_percent=100*(combined/previous-1),
        reference_score=8.43991,reference_score_ratio=combined/8.43991,
        other_tracks_reverified_in_this_run=False,new_kaggle_submission=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    screen,build=check_inputs();inputs=input_hashes()
    verifier,official=load_verifier(args.data);baseline,_=load_comparison(args.data)
    spec=verifier.ProblemSpecInvSqrt.for_n_p(10)
    if (spec.width,spec.max_x,spec.precision_bits,spec.shots_per_x,spec.num_shots)!=(22,3_133_452,24,64,200_540_928):
        raise ValueError('Full official benchmark specification changed.')
    artifact=gen.HERE/'local'/build['artifact'];circuit=km.Circuit.from_file(artifact)
    metadata=dict(design=build['design'],build=build,input_sha256=inputs,
        verifier_sha256=VERIFIER_SHA256,official_verifier_modified=False,salt=SALT,
        python=platform.python_version(),numpy=np.__version__,kickmix=km.__version__,
        protocol=dict(width=spec.width,max_x=spec.max_x,precision_bits=spec.precision_bits,
            shots_per_x=spec.shots_per_x,verified_shots=spec.num_shots,beta=str(spec.beta),
            batch_size=verifier.DEFAULT_BATCH_SIZE,max_workers=args.workers))
    report_path=gen.HERE/'verified.json'
    report_path.write_text(json.dumps({**metadata,'status':'RUNNING'},indent=2)+'\n')
    print(json.dumps(dict(status='RUNNING',shots=spec.num_shots,num_qubits=build['num_qubits'],
                        toffoli_count=build['toffoli_count'],circuit_sha256=build['circuit_sha256'])),flush=True)
    started=time.monotonic()
    try:
        result=verifier.verify(circuit,spec,salt=SALT,max_workers=args.workers)
        if (result['status']!='PASS' or result['verified_shots']!=spec.num_shots or
                result['circuit_sha256']!=build['circuit_sha256'] or
                result['num_qubits']!=build['num_qubits'] or
                gen.digest(artifact)!=build['circuit_sha256'] or
                gen.digest(official)!=VERIFIER_SHA256 or inputs!=input_hashes()):
            raise RuntimeError('Full verification coverage, artifact identity or input stability failed.')
        report={**metadata,**result,'comparison':compare(result['score'],baseline),
            'model_estimated_cost':screen['variants'][0]['estimated_cost_from_exact_model_margin'],
            'elapsed_seconds':time.monotonic()-started}
        report_path.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({key:report[key] for key in ['status','score','num_qubits','toffoli_count',
              'toffoli_equivalent','rotation_count','verified_shots','comparison','elapsed_seconds']},indent=2),flush=True)
    except Exception as error:
        report_path.write_text(json.dumps({**metadata,'status':'FAIL','failure_type':type(error).__name__,
            'elapsed_seconds':time.monotonic()-started},indent=2)+'\n')
        raise


if __name__=='__main__':main()
