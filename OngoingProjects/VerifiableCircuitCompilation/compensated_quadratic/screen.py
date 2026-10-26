"""Exhaustive integer model and compiled boundary probes for compensation."""
import argparse
import json
import math
from pathlib import Path

import kickmix as km
import numpy as np

import generator as gen
from verify_quadratic import load_verifier, VERIFIER_SHA256

MAX_X=3_133_452
VARIANTS=[dict(arithmetic='padded',recycle_coefficient=True,adder_budget=None),
          dict(arithmetic='padded',recycle_coefficient=False,adder_budget=None),
          dict(arithmetic='controlled',recycle_coefficient=True,adder_budget=None),
          dict(arithmetic='padded',recycle_coefficient=True,adder_budget=130)]


def model(xs):
    xs=np.asarray(xs,dtype=np.int64)
    if xs.ndim!=1 or not ((xs>=1)&(xs<=MAX_X)).all():
        raise ValueError('Expected a one-dimensional list of legal positive inputs.')
    shifts=(22-np.searchsorted(1<<np.arange(23),xs,side='right'))//2
    normalized=xs<<(2*shifts)
    address=(normalized>>14)-64
    table=np.asarray(gen.quantized_rows(),dtype=np.int64)
    lengths=[int(max(table[:,j])).bit_length() for j in range(3)]
    values=[table[address,j].copy() for j in range(3)]
    bad=np.zeros(len(xs),dtype=bool)
    for j in (1,0):
        for m in range(7,21):
            shift=gen.FRACTIONAL_BITS[j]-gen.FRACTIONAL_BITS[j+1]-m
            term=values[j+1]<<shift if shift>=0 else values[j+1]>>(-shift)
            values[j]-=((normalized>>(20-m))&1)*term
            bad|=(values[j]<=0)|(values[j]>=2**lengths[j])
    return (values[0]<<shifts).astype(np.uint64)<<np.uint64(24),bad


def exhaustive(verifier):
    xs=np.arange(1,MAX_X+1,dtype=np.uint64)
    target,tolerance=verifier.reference_values(verifier.ProblemSpecInvSqrt.for_n_p(10),xs)
    fail=capacity=0;worst=-1.;worst_x=0;margin=1<<64
    for start in range(0,MAX_X,65536):
        sl=slice(start,min(start+65536,MAX_X))
        phase,bad=model(xs[sl])
        used=np.abs((phase-target[sl]).view(np.int64)).astype(np.uint64)
        ratios=used/tolerance[sl];i=int(ratios.argmax())
        if float(ratios[i])>worst:worst=float(ratios[i]);worst_x=int(xs[sl][i])
        fail+=int((used>tolerance[sl]).sum());capacity+=int(bad.sum())
        margin=min(margin,int((tolerance[sl].astype(np.int64)-used.astype(np.int64)).min()))
    if fail or capacity or margin<=0:
        raise AssertionError('Exhaustive model failed phase, capacity or synthesis margin.')
    return dict(status='INTEGER_MODEL_PASS',checked_inputs=MAX_X,phase_failures=fail,
        partial_accumulator_violations=capacity,minimum_margin_units=margin,
        worst_error_over_tolerance=worst,worst_x=worst_x)


def boundary_sample(circuit,verifier):
    inputs={1,MAX_X}
    for shift in range(11):
        for address in range(64,257):
            boundary=(address<<14)>>(2*shift)
            inputs.update(boundary+offset for offset in (-1,0,1) if 1<=boundary+offset<=MAX_X)
    xs=np.tile(np.array(sorted(inputs),dtype=np.uint64),4)
    xreg,scratch=verifier.validate_registers(circuit,verifier.ProblemSpecInvSqrt.for_n_p(10))
    sim=km.Simulator(batch_size=len(xs),counting_operations=True,seed=371811)
    sim.write_across_shots(xreg,verifier.to_bits(xs,22));sim.do(circuit)
    expected,bad=model(xs)
    actual=np.array([verifier.half_turns_to_units(sim.read_shot_phase(k)) for k in range(len(xs))],dtype=np.uint64)
    result=dict(unique_inputs=len(inputs),shots=len(xs),seed=371811,
        input_restored=bool(np.array_equal(verifier.from_bits(sim.read_across_shots(xreg)),xs)),
        scratch_clean=not bool(sim.read_across_shots(scratch).any()),
        exact_phase_model_match=bool(np.array_equal(actual,expected)),
        register_fit=not bool(bad.any()))
    static,measured=circuit.operation_statistics(),sim.operation_statistics()
    result['priced_counts_match']=all(a==len(xs)*b for a,b in [
        (measured.num_toffoli_like,static.num_toffoli_like),(measured.num_t_like,static.num_t_like),
        (measured.num_z_pow(min_clifford_hierarchy_level=4),static.num_z_pow(min_clifford_hierarchy_level=4))])
    if not all(result[key] is True for key in ['input_restored','scratch_clean','exact_phase_model_match','register_fit','priced_counts_match']):
        raise AssertionError('Compiled boundary, phase or cleanup probe failed.')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--data',type=Path,required=True)
    args=parser.parse_args();verifier,_=load_verifier(args.data)
    gen.checked_rows();screen=exhaustive(verifier)
    records=[]
    for options in VARIANTS:
        circuit=gen.build_circuit(**options);stats=circuit.operation_statistics()
        rotations=stats.num_z_pow(min_clifford_hierarchy_level=4)
        eps=math.pi*screen['minimum_margin_units']/(2**64*rotations)
        estimate=circuit.num_qubits*(stats.num_toffoli_like+
                   (stats.num_t_like+rotations*verifier.rotation_t_count(eps))/4)
        records.append(dict(options=options,num_qubits=circuit.num_qubits,
            toffoli_count=stats.num_toffoli_like,rotation_count=rotations,
            estimated_cost_from_exact_model_margin=estimate,
            boundary=boundary_sample(circuit,verifier)))
    records.sort(key=lambda record:record['estimated_cost_from_exact_model_margin'])
    selected=records[0]['options']
    build=gen.build_to_file(gen.HERE/'local/compensated-quadratic.kmb',**selected)
    if gen.canonical_bytes(gen.build_circuit())!=(gen.HERE/'local'/build['artifact']).read_bytes():
        raise AssertionError('Standalone default does not reproduce the selected circuit.')
    (gen.HERE/'build.json').write_text(json.dumps(build,indent=2)+'\n')
    report=dict(status='SCREENING_ONLY_NOT_OFFICIAL_VERIFICATION',model=screen,variants=records,
        standalone_default_matches_selected=True,
        selected_options=selected,generator_sha256=gen.digest(gen.HERE/'generator.py'),
        screen_source_sha256=gen.digest(Path(__file__)),
        coefficients_sha256=gen.digest(gen.HERE/'coefficients.json'),verifier_sha256=VERIFIER_SHA256,
        circuit_sha256=build['circuit_sha256'])
    (gen.HERE/'screening.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
