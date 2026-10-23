"""Reduced-width quadratic circuit, derived from the frozen 36-bit table.

The generator needs kickmix, quadratic.py and its checked-in coefficient JSON.
The competition code and Wolfram kernel are verification dependencies only.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import kickmix as km

from quadratic import (HERE, DESIGN_SHA256, TABLE_SHA256, _normalize,
                       canonical_bytes, load_coefficients)

FRACTIONAL_BITS = (30, 23, 16)


def quantized_rows(fractional_bits=FRACTIONAL_BITS):
    if (len(fractional_bits) != 3 or
            any(type(n) is not int or not 1 <= n <= 36 for n in fractional_bits)):
        raise ValueError("Expected three integer fractional widths from 1 through 36.")
    return tuple(tuple(round(Fraction(abs(c), 2**(36-f)))
                       for c, f in zip(row, fractional_bits))
                 for row in load_coefficients())


def _horner_stage(builder, x, registers, widths, stage, inverse, arithmetic):
    accumulator, factor = registers[stage], registers[stage+1]
    scratch = builder.alloc_qubits(len(accumulator)) if arithmetic == "padded" else None
    for m in range(7, 21):
        # floor(2**f_acc * (factor / 2**f_factor) * 2**-m).
        shift = m + widths[stage+1] - widths[stage]
        sources = factor[max(0, shift):]
        target = accumulator[max(0, -shift):]
        if not len(sources) or not len(target):
            continue
        add = builder.iadd if inverse else builder.isub
        if scratch is None:
            add(sources, target=target, control=x[20-m])
        else:
            products = scratch[:len(sources)]
            builder.init_and(x[20-m], sources, products)
            add(scratch[:len(target)], target=target)
            builder.del_and(x[20-m], sources, products)
    if scratch is not None:
        builder.free(scratch)


def build_circuit(fractional_bits=FRACTIONAL_BITS, *, arithmetic="padded",
                  recycle_coefficient=True, adder_budget=None):
    """Build a candidate; only the default configuration has full-run evidence."""
    rows = quantized_rows(fractional_bits)
    if arithmetic not in ("padded", "controlled"):
        raise ValueError("Unknown arithmetic implementation.")
    lengths = [max(row[j] for row in rows).bit_length() for j in range(3)]
    if not all(lengths):
        raise ValueError("A coefficient register would be empty.")
    builder = km.CircuitBuilder()
    if adder_budget is not None:
        builder.set_max_qubits(adder_budget)
    x = builder.create_quantum_register(22, name="x")
    flags = builder.alloc_qubits(10)
    registers = [builder.alloc_qubits(n) for n in lengths]
    table = [0]*64 + [sum(c << sum(lengths[:j]) for j, c in enumerate(row))
                     for row in rows]
    c2_table = [0]*64 + [row[2] for row in rows]
    address = x[14:]

    def lookup_register():
        return km.array([q for register in registers for q in register])

    def stage(index, inverse):
        _horner_stage(builder, x, registers, fractional_bits, index, inverse, arithmetic)

    with builder.mark("normalize"):
        _normalize(builder, x, flags, False)
    with builder.mark("coefficient_lookup"):
        builder.init_lookup(table=table, address=address, target=lookup_register())
    with builder.mark("first_horner"):
        stage(1, False)
    if recycle_coefficient:
        with builder.mark("release_c2"):
            builder.del_lookup(table=c2_table, address=address, target=registers[2])
            builder.free(registers[2])
    with builder.mark("second_horner"):
        stage(0, False)

    # These zeros are live only for phase extraction, not during arithmetic.
    extension = builder.alloc_qubits(10)
    word = km.array([*registers[0], *extension])
    with builder.mark("phase"):
        for step in range(10):
            builder.left_rotate(word[step:step+lengths[0]+1], 1, control=flags[step])
        for bit, qubit in enumerate(word):
            builder.z_pow(qubit, Fraction(2)**(bit-fractional_bits[0]-10+1))
        for step in reversed(range(10)):
            builder.right_rotate(word[step:step+lengths[0]+1], 1, control=flags[step])
    builder.free(extension)

    with builder.mark("inverse_second_horner"):
        stage(0, True)
    if recycle_coefficient:
        with builder.mark("reload_c2"):
            registers[2] = builder.alloc_qubits(lengths[2])
            builder.init_lookup(table=c2_table, address=address, target=registers[2])
    with builder.mark("inverse_first_horner"):
        stage(1, True)
    with builder.mark("coefficient_unlookup"):
        builder.del_lookup(table=table, address=address, target=lookup_register())
    with builder.mark("inverse_normalize"):
        _normalize(builder, x, flags, True)
    for register in [flags, *registers]:
        builder.free(register)
    return builder.finish_circuit()


def build_to_file(output: Path):
    circuit = build_circuit()
    binary = canonical_bytes(circuit)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(binary)
    stats = circuit.operation_statistics()
    return {"design": "B-quadratic-192-reduced", "artifact": output.name,
            "fractional_bits": list(FRACTIONAL_BITS),
            "design_sha256": DESIGN_SHA256, "source_table_sha256": TABLE_SHA256,
            "quantized_table_sha256": hashlib.sha256(json.dumps(
                quantized_rows(), separators=(",", ":")).encode()).hexdigest(),
            "circuit_sha256": hashlib.sha256(binary).hexdigest(),
            "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "generator_dependency_sha256": hashlib.sha256((HERE/"quadratic.py").read_bytes()).hexdigest(),
            "kickmix_version": km.__version__, "num_qubits": circuit.num_qubits,
            "toffoli_count": stats.num_toffoli_like, "t_count": stats.num_t_like,
            "rotation_count": stats.num_z_pow(min_clifford_hierarchy_level=4),
            "reaction_depth": circuit.reaction_depth(), "num_ops": len(circuit)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE/"results/optimized-quadratic.kmb")
    args = parser.parse_args()
    manifest = build_to_file(args.output)
    args.output.with_suffix(".build.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps(manifest, indent=2))
