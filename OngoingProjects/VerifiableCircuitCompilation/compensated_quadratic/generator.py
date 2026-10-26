"""Compensated-curvature phase circuit from the frozen quadratic source.

Build requires kickmix, the parent arithmetic modules, frozen source JSON and
this directory's WL-derived coefficients.json. No competition import or
Wolfram runtime is used by the generator.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import kickmix as km

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
sys.path.insert(0, str(PROJECT))
from quadratic import DESIGN_SHA256, _normalize, canonical_bytes, load_coefficients
from optimized_quadratic import _horner_stage

FRACTIONAL_BITS = (30, 23, 10)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def quantized_rows(widths=FRACTIONAL_BITS):
    if (len(widths) != 3 or any(type(n) is not int or not 1 <= n <= 36 for n in widths)):
        raise ValueError("Expected three integer fractional widths in 1..36.")
    result = []
    h = Fraction(1, 64)
    for row in load_coefficients():
        a, b, c = (Fraction(abs(value), 2**36) for value in row)
        coarse = Fraction(round(c * 2**widths[2]), 2**widths[2])
        delta = coarse-c
        result.append((round((a+delta*h*h/8)*2**widths[0]),
                       round((b+delta*h)*2**widths[1]), int(coarse*2**widths[2])))
    return tuple(result)


def checked_rows():
    report = json.loads((HERE/'coefficients.json').read_text())
    if (report['status'] != 'DERIVATION_AND_CAPACITY_PASS'
            or report['fractional_bits'] != list(FRACTIONAL_BITS)
            or report['source_sha256'] != DESIGN_SHA256
            or report['wl_source_sha256'] != digest(HERE/'CompensatedQuadratic.wl')
            or report['runner_sha256'] != digest(HERE/'run.wls')
            or not all(value is True for value in report['derivation']['checks'].values())
            or report['certificate']['all_capacity_proved'] is not True
            or len(report['certificate']['intervals']) != 192
            or any(row['capacity'] is not True for row in report['certificate']['intervals'])
            or report['table'] != [list(row) for row in quantized_rows()]):
        raise ValueError('WL derivation, source binding or coefficient correspondence failed.')
    return quantized_rows()


def build_circuit(*, arithmetic='padded', recycle_coefficient=False, adder_budget=None):
    rows = checked_rows()
    if arithmetic not in ('padded', 'controlled'):
        raise ValueError('Unknown arithmetic implementation.')
    lengths = [max(row[j] for row in rows).bit_length() for j in range(3)]
    if min(lengths) <= 0:
        raise ValueError('Empty coefficient register.')
    builder = km.CircuitBuilder()
    if adder_budget is not None:
        builder.set_max_qubits(adder_budget)
    x = builder.create_quantum_register(22, name='x')
    flags = builder.alloc_qubits(10)
    registers = [builder.alloc_qubits(n) for n in lengths]
    table = [0]*64 + [sum(c << sum(lengths[:j]) for j,c in enumerate(row)) for row in rows]
    c2_table = [0]*64 + [row[2] for row in rows]
    address = x[14:]
    lookup_register = lambda: km.array([q for register in registers for q in register])

    def stage(index, inverse):
        _horner_stage(builder, x, registers, FRACTIONAL_BITS, index, inverse, arithmetic)

    with builder.mark('normalize'):
        _normalize(builder, x, flags, False)
    with builder.mark('coefficient_lookup'):
        builder.init_lookup(table=table, address=address, target=lookup_register())
    with builder.mark('first_horner'):
        stage(1, False)
    if recycle_coefficient:
        with builder.mark('release_c2'):
            builder.del_lookup(table=c2_table, address=address, target=registers[2])
            builder.free(registers[2])
    with builder.mark('second_horner'):
        stage(0, False)
    extension = builder.alloc_qubits(10)
    word = km.array([*registers[0], *extension])
    with builder.mark('phase'):
        for step in range(10):
            builder.left_rotate(word[step:step+lengths[0]+1], 1, control=flags[step])
        for bit, qubit in enumerate(word):
            builder.z_pow(qubit, Fraction(2)**(bit-FRACTIONAL_BITS[0]-10+1))
        for step in reversed(range(10)):
            builder.right_rotate(word[step:step+lengths[0]+1], 1, control=flags[step])
    builder.free(extension)
    with builder.mark('inverse_second_horner'):
        stage(0, True)
    if recycle_coefficient:
        with builder.mark('reload_c2'):
            registers[2] = builder.alloc_qubits(lengths[2])
            builder.init_lookup(table=c2_table, address=address, target=registers[2])
    with builder.mark('inverse_first_horner'):
        stage(1, True)
    with builder.mark('coefficient_unlookup'):
        builder.del_lookup(table=table, address=address, target=lookup_register())
    with builder.mark('inverse_normalize'):
        _normalize(builder, x, flags, True)
    for register in [flags, *registers]:
        builder.free(register)
    return builder.finish_circuit()


def build_to_file(output, **options):
    circuit = build_circuit(**options)
    binary = canonical_bytes(circuit)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(binary)
    stats = circuit.operation_statistics()
    return dict(design='compensated-quadratic-192',artifact=output.name,
        fractional_bits=list(FRACTIONAL_BITS),physical_coefficient_widths=[30,22,9],
        arithmetic=options.get('arithmetic','padded'),
        recycle_coefficient=options.get('recycle_coefficient',False),
        adder_budget=options.get('adder_budget'),circuit_sha256=hashlib.sha256(binary).hexdigest(),
        coefficients_sha256=digest(HERE/'coefficients.json'),generator_sha256=digest(Path(__file__)),
        dependency_sha256={name:digest(PROJECT/name) for name in
                          ['quadratic.py','optimized_quadratic.py','feasibility/B-quadratic-192.json']},
        num_qubits=circuit.num_qubits,toffoli_count=stats.num_toffoli_like,
        t_count=stats.num_t_like,rotation_count=stats.num_z_pow(min_clifford_hierarchy_level=4),
        reaction_depth=circuit.reaction_depth(),num_ops=len(circuit),kickmix_version=km.__version__)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'local/compensated-quadratic.kmb')
    args=parser.parse_args()
    report=build_to_file(args.output)
    (HERE/'build.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
