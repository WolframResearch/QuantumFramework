"""Standalone kickmix generator for the certified, fixed B-quadratic-192 design.

Building needs kickmix and the checked-in certificate JSON only. It does not
import the competition generator, verifier, or feasibility research code.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path

import kickmix as km

HERE = Path(__file__).resolve().parent
DESIGN = HERE / "feasibility/B-quadratic-192.json"
DESIGN_SHA256 = "2b768e13c7054d0e58bc4edf7afc0ff64559685aa309ad250b275fd9fd6a4c8d"
TABLE_SHA256 = "a849cf5f46c5a39ba4c8c3b46fdaf49c28cd1bf375bd2d69d778ac6a57ca9a32"
WIDTH, SHIFT_BITS, COEF_BITS, INTERVAL_BITS = 22, 10, 36, 6


def load_coefficients(path: Path = DESIGN) -> tuple[tuple[int, ...], ...]:
    """Reject any change to the reviewed design document or its integer table."""
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != DESIGN_SHA256:
        raise ValueError("The fixed certified design document has changed.")
    design = json.loads(raw)
    expected = {"name": "B-quadratic-192", "family": "B", "degree": 2,
                "interval_bits": INTERVAL_BITS, "coef_bits": COEF_BITS,
                "delta_bits": 20, "certified": True}
    if any(design.get(k) != v for k, v in expected.items()):
        raise ValueError("Unexpected certified design parameters.")
    table = design["table"]
    digest = hashlib.sha256(json.dumps(table, separators=(",", ":")).encode()).hexdigest()
    if digest != TABLE_SHA256:
        raise ValueError("The certified integer coefficients have changed.")
    return tuple(tuple(row) for row in table)


def canonical_bytes(circuit: km.Circuit) -> bytes:
    stream = io.BytesIO()
    circuit.to_file(stream, format="kmb")
    return stream.getvalue()


def _normalize(builder, x, flags, inverse):
    """Shift x by two bits per unary flag, preserving reversibility."""
    steps = reversed(range(SHIFT_BITS)) if inverse else range(SHIFT_BITS)
    for step in steps:
        lanes = (x[2 * step::2], x[2 * step + 1::2])
        if inverse:
            for lane in lanes:
                builder.right_rotate(lane, 1, control=flags[step])
        builder.x(x[-2:])
        if inverse:
            builder.del_and(x[-2], x[-1], flags[step])
        else:
            builder.init_and(x[-2], x[-1], flags[step])
        builder.x(x[-2:])
        if not inverse:
            for lane in lanes:
                builder.left_rotate(lane, 1, control=flags[step])


def _horner(builder, x, coefficients, scratch, inverse):
    """c1 -= T(delta*c2); c0 -= T(delta*c1), with row-wise truncation."""
    pairs = [(coefficients[1], coefficients[2]), (coefficients[0], coefficients[1])]
    for accumulator, factor in reversed(pairs) if inverse else pairs:
        for fraction_bit in range(INTERVAL_BITS + 1, WIDTH - 1):
            control = x[WIDTH - 2 - fraction_bit]
            sources = km.array([factor[i] for i in range(fraction_bit, COEF_BITS)])
            products = scratch[:COEF_BITS - fraction_bit]
            builder.init_and(control, sources, products)
            if inverse:
                builder.iadd(scratch, target=accumulator)
            else:
                builder.isub(scratch, target=accumulator)
            builder.del_and(control, sources, products)


def _shift_phase_word(builder, flags, word, inverse):
    steps = reversed(range(SHIFT_BITS)) if inverse else range(SHIFT_BITS)
    for step in steps:
        lane = word[step:step + COEF_BITS + 1]
        if inverse:
            builder.right_rotate(lane, 1, control=flags[step])
        else:
            builder.left_rotate(lane, 1, control=flags[step])


def build_circuit() -> km.Circuit:
    coefficients = load_coefficients()
    # Normalized M is in [1,4). Addresses below 1 are unused but still charged.
    table = [0] * 2**INTERVAL_BITS + [
        sum(abs(value) << (COEF_BITS * degree) for degree, value in enumerate(row))
        for row in coefficients]
    builder = km.CircuitBuilder()
    x = builder.create_quantum_register(WIDTH, name="x")
    flags = builder.alloc_qubits(SHIFT_BITS)
    registers = [builder.alloc_qubits(COEF_BITS) for _ in range(3)]
    extension = builder.alloc_qubits(SHIFT_BITS)
    scratch = builder.alloc_qubits(COEF_BITS)
    address = x[WIDTH - INTERVAL_BITS - 2:]
    lookup = km.array([qubit for register in registers for qubit in register])
    phase_word = km.array([*registers[0], *extension])

    def compute(inverse):
        with builder.mark("uncompute_poly" if inverse else "compute_poly"):
            if not inverse:
                with builder.mark("normalize"):
                    _normalize(builder, x, flags, False)
                with builder.mark("coefficient_lookup"):
                    builder.init_lookup(table=table, address=address, target=lookup)
                with builder.mark("horner"):
                    _horner(builder, x, registers, scratch, False)
            with builder.mark("denormalize"):
                _shift_phase_word(builder, flags, phase_word, inverse)
            if inverse:
                with builder.mark("horner"):
                    _horner(builder, x, registers, scratch, True)
                with builder.mark("coefficient_unlookup"):
                    builder.del_lookup(table=table, address=address, target=lookup)
                with builder.mark("normalize"):
                    _normalize(builder, x, flags, True)

    compute(False)
    with builder.mark("phase_by_bits"):
        for bit, qubit in enumerate(phase_word):
            builder.z_pow(qubit, Fraction(2)**(bit - len(phase_word) + 1))
    compute(True)
    for register in [flags, *registers, extension, scratch]:
        builder.free(register)
    return builder.finish_circuit()


def build_to_file(output: Path) -> dict:
    circuit = build_circuit()
    binary = canonical_bytes(circuit)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(binary)
    statistics = circuit.operation_statistics()
    return {"design": "B-quadratic-192", "artifact": output.name,
            "design_sha256": DESIGN_SHA256, "coefficient_table_sha256": TABLE_SHA256,
            "circuit_sha256": hashlib.sha256(binary).hexdigest(),
            "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "kickmix_version": km.__version__, "num_qubits": circuit.num_qubits,
            "toffoli_count": statistics.num_toffoli_like,
            "t_count": statistics.num_t_like,
            "rotation_count": statistics.num_z_pow(min_clifford_hierarchy_level=4),
            "reaction_depth": circuit.reaction_depth(), "num_ops": len(circuit)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "results/quadratic.kmb")
    args = parser.parse_args()
    result = build_to_file(args.output)
    args.output.with_suffix(".build.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
