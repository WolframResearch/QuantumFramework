"""Resource-only invsqrt architecture probes using the official, unpatched kickmix.

This script does not produce a submission and never calls the full verifier.  It
measures compiled schematic circuits; approximation certificates remain separate
inputs and do not certify that the schematics implement those approximations.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import dataclasses
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time

import kickmix as km
import numpy as np


HERE = Path(__file__).resolve().parent
N_P = 10
WIDTH = 22
MAX_X = 3_133_452
TOFFOLI_NAMES = {"CCX", "CCZ", "CSWAP"}
STAGES = {"normalize", "select_interval", "load_coefficients", "coefficient_lookup",
          "coefficient_unlookup", "horner", "square", "multiply_m", "correct",
          "multiply_y0", "denormalize", "phase_by_bits", "phase_polynomial",
          "exponent_indicator"}


def load_official(data: Path):
    sys.path.insert(0, str(data / "verifiers"))
    from invsqrt import nr_reference_circuit as nr
    from invsqrt import verify_submission_invsqrt as verifier
    return nr, verifier


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _stage_key(stack: list[str]) -> str:
    for i, mark in enumerate(stack):
        if mark in STAGES:
            return "/".join(stack[:i + 1])
    return stack[0] if stack else "outside_tags"


def tagged_profile(circuit) -> dict:
    """Parse actual INFO_TAG nesting and reconcile disjoint stage gate counts."""
    rows = defaultdict(Counter)
    stack = []
    for line in str(circuit).splitlines():
        if line.startswith('INFO_TAG "@enter:'):
            stack.append(line[len('INFO_TAG "@enter:'):-1])
        elif line.startswith('INFO_TAG "@exit:'):
            name = line[len('INFO_TAG "@exit:'):-1]
            if not stack or stack.pop() != name:
                raise ValueError(f"Unbalanced INFO_TAG: {line}")
        else:
            key = _stage_key(stack)
            op = line.split()[0] if line.split() else ""
            row = rows[key]
            row["instruction_count_excluding_tags"] += 1
            if op in TOFFOLI_NAMES:
                row["toffoli_like"] += 1
                row[op] += 1
                if " if " in line:
                    row["classically_conditioned_toffoli_like"] += 1
            if " if " in line:
                row["classically_conditioned_instructions"] += 1
            if op == "HMR":
                row["HMR"] += 1
            if op == "Z_POW":
                half_turns = Fraction(line.split()[2])
                angle_class = ("clifford" if half_turns % Fraction(1, 2) == 0 else
                               "t" if half_turns % Fraction(1, 4) == 0 else "fine")
                row[f"z_pow_{angle_class}"] += 1
                if " if " in line:
                    row[f"classically_conditioned_z_pow_{angle_class}"] += 1
    if stack:
        raise ValueError("Unclosed INFO_TAG")
    stats = circuit.operation_statistics()
    toffolis = sum(v["toffoli_like"] for v in rows.values())
    rotations = sum(v["z_pow_fine"] for v in rows.values())
    t_angles = sum(v["z_pow_t"] for v in rows.values())
    assert toffolis == stats.num_toffoli_like
    assert rotations == stats.num_z_pow(min_clifford_hierarchy_level=4)
    assert t_angles == stats.num_t_like
    return {"stages": dict(sorted((k, dict(v)) for k, v in rows.items())),
            "sum_stage_toffoli_like": toffolis,
            "sum_stage_fine_rotations": rotations,
            "sum_stage_t_angles": t_angles,
            "classically_conditioned_toffoli_like": sum(
                v["classically_conditioned_toffoli_like"] for v in rows.values()),
            "classically_conditioned_base_t_angles": sum(
                v["classically_conditioned_z_pow_t"] for v in rows.values()),
            "classically_conditioned_fine_rotations": sum(
                v["classically_conditioned_z_pow_fine"] for v in rows.values()),
            "reconciles_operation_statistics": True,
            "counting_caveat": "Static counts include all classically conditioned operations. "
                "The official objective uses mean executed counts from simulation; these "
                "resource probes use static upper counts and have no full-verifier claim."}


def circuit_stats(circuit) -> dict:
    stats = circuit.operation_statistics()
    return {"num_qubits": circuit.num_qubits,
            "num_operations": len(circuit),
            "toffoli_like": stats.num_toffoli_like,
            "base_t_angles": stats.num_t_like,
            "fine_rotations": stats.num_z_pow(min_clifford_hierarchy_level=4),
            "reaction_depth": circuit.reaction_depth(),
            "operation_counts": dict(stats.items())}


def margin_objective(stats: dict, fraction, verifier) -> dict:
    """Price a hypothetical fraction of the smallest allowed absolute phase error."""
    fraction = Fraction(str(fraction))
    # Identical integer tolerance rule to reference_values, evaluated at max_x.
    minimum_tolerance_units = math.isqrt((1 << 130) // MAX_X) >> 25
    margin_units = fraction.numerator * minimum_tolerance_units // fraction.denominator
    record = {"margin_fraction_of_minimum_tolerance": str(fraction),
              "minimum_tolerance_units": minimum_tolerance_units,
              "assumed_margin_units": margin_units,
              "status": "CONDITIONAL_RESOURCE_ESTIMATE"}
    if margin_units <= 0:
        return {**record, "status": "NO_POSITIVE_SYNTHESIS_MARGIN"}
    rotations = stats["fine_rotations"]
    eps = math.pi * margin_units / (2**64 * rotations) if rotations else None
    synthesis_t = rotations * verifier.rotation_t_count(eps) if rotations else 0
    cost = stats["toffoli_like"] + (stats["base_t_angles"] + synthesis_t) / 4
    return {**record, "rotation_eps": eps, "synthesis_t_count": synthesis_t,
            "toffoli_equivalent_static": cost,
            "estimated_objective_Q_times_C": stats["num_qubits"] * cost}


def _coefficient_table(design: dict) -> list[int]:
    degree, r, w = (int(design[k]) for k in ("degree", "interval_bits", "coef_bits"))
    rows = design["table"]
    if len(rows) != 3 * 2**r:
        raise ValueError(f"Expected {3 * 2**r} interval rows")
    packed = [0] * 2**r
    for row in rows:
        if len(row) != degree + 1:
            raise ValueError("Coefficient count does not match degree")
        for j, value in enumerate(row):
            if (int(value) < 0) != bool(j % 2) and int(value) != 0:
                raise ValueError("Expected alternating signs, with c0 positive")
            if abs(int(value)) >= 2**w:
                raise ValueError("Coefficient magnitude does not fit its allocated register")
        packed.append(sum(abs(int(value)) << (w * j) for j, value in enumerate(row)))
    assert len(packed) == 2**(r + 2)
    return packed


def _poly_horner(builder, inverse, design, x, coefs, scratch, nr):
    r, d, w = (int(design[k]) for k in ("interval_bits", "delta_bits", "coef_bits"))
    delta = [(x[WIDTH - 2 - m], m) for m in range(r + 1, min(d, WIDTH - 2) + 1)]
    steps = [(coefs[j], coefs[j + 1]) for j in reversed(range(len(coefs) - 1))]
    for acc, mult in reversed(steps) if inverse else steps:
        for control, m in delta:
            terms = [(mult[i], i - m) for i in range(m, w)]
            if terms:
                nr._add_row(builder, inverse, control, terms, acc, scratch, subtract=True)


def build_poly_probe(design: dict, nr, *, nr_bits: int | None = None):
    """Build a resource schematic with one joint coefficient lookup and its deletion.

    Horner uses magnitudes with alternating subtraction, as in the official
    generator.  A appends the official NR stages.  B phases the direct polynomial
    with ten extra clean bits, preserving coefficient storage for uncomputation.
    """
    family = design["family"]
    if family not in ("A", "B"):
        raise ValueError(f"Unknown architecture family {family}")
    r, w, degree = (int(design[k]) for k in ("interval_bits", "coef_bits", "degree"))
    if not 0 <= r <= 19:
        raise ValueError("Address must fit normalized 22-bit input")
    table = _coefficient_table(design)
    nrb = int(nr_bits if nr_bits is not None else design.get("nr_bits", 32))
    builder = km.CircuitBuilder()
    x = builder.create_quantum_register(WIDTH, name="x")
    flags = builder.alloc_qubits(N_P)
    coefs = [builder.alloc_qubits(w) for _ in range(degree + 1)]
    address = x[WIDTH - (r + 2):]
    lookup_target = km.array([q for coef in coefs for q in coef])
    owned = [flags, *coefs]
    if family == "A":
        square = builder.alloc_qubits(nrb)
        z = builder.alloc_qubits(nrb + 2)
        y1 = builder.alloc_qubits(nrb + N_P + 1)
        scratch = builder.alloc_qubits(max(w, nrb + 2))
        owned += [square, z, y1, scratch]
        p = y1[:nrb + 1]
    else:
        extra = builder.alloc_qubits(N_P)
        y1 = km.array([*coefs[0], *extra])
        scratch = builder.alloc_qubits(w)
        owned += [extra, scratch]

    def compute(inverse: bool):
        with builder.mark("uncompute_poly" if inverse else "compute_poly"):
            if not inverse:
                with builder.mark("normalize"):
                    nr._normalize(builder, False, x, flags)
                with builder.mark("coefficient_lookup"):
                    builder.init_lookup(table=table, address=address, target=lookup_target)
                with builder.mark("horner"):
                    _poly_horner(builder, False, design, x, coefs, scratch, nr)
            if family == "A":
                stages = [
                    ("square", lambda: nr._square(builder, inverse, coefs[0], square, scratch)),
                    ("multiply_m", lambda: nr._multiply_m(builder, inverse, x, square, z, scratch)),
                    ("correct", lambda: nr._correct(builder, inverse, z)),
                    ("multiply_y0", lambda: nr._multiply_y0(builder, inverse, coefs[0], z, p, scratch)),
                    ("denormalize", lambda: nr._denormalize(builder, inverse, flags, y1))]
                for name, action in reversed(stages) if inverse else stages:
                    with builder.mark(name):
                        action()
            else:
                # len(y1)=w+10 gives wp=w-1 in the official helper, hence each
                # lane has w active bits plus the clean next shift destination.
                with builder.mark("denormalize"):
                    nr._denormalize(builder, inverse, flags, y1)
            if inverse:
                with builder.mark("horner"):
                    _poly_horner(builder, True, design, x, coefs, scratch, nr)
                with builder.mark("coefficient_unlookup"):
                    builder.del_lookup(table=table, address=address, target=lookup_target)
                with builder.mark("normalize"):
                    nr._normalize(builder, True, x, flags)

    compute(False)
    nr._phase_by_bits(builder, y1)
    compute(True)
    for reg in owned:
        builder.free(reg)
    circuit = builder.finish_circuit()
    metadata = {"family": family, "degree": degree, "interval_bits": r,
                "coef_bits": w, "delta_bits": int(design["delta_bits"]),
                "interval_count": 3 * 2**r, "lookup_entries_including_padding": len(table),
                "lookup_target_bits": len(lookup_target),
                "lookup_strategy": "one concatenated coefficient-magnitude table; init_lookup then del_lookup",
                "explicitly_allocated_qubits": WIDTH + sum(len(reg) for reg in owned),
                "phase_fractional_bits": len(y1), "resource_only": True}
    if family == "A":
        metadata["nr_precision"] = {k: nrb for k in ("square_bits", "product_bits", "result_bits")}
    return circuit, metadata


def poly_model_words(design: dict, xs, *, nr_bits=None):
    """Bit-exact integer model of the schematic, with register-capacity evidence.

    Every multiplier drops individual low partial-product bits, rather than
    rounding an exact product.  Arithmetic is reduced modulo the actual target
    width after recording any capacity violation, matching finite registers.
    """
    xs = np.asarray(xs, dtype=np.int64)
    if not ((1 <= xs) & (xs <= MAX_X)).all():
        raise ValueError("Model inputs outside the benchmark domain")
    r, w, d = (int(design[k]) for k in ("interval_bits", "coef_bits", "delta_bits"))
    nrb = int(nr_bits if nr_bits is not None else design.get("nr_bits", 32))
    bitlength = np.searchsorted(np.left_shift(1, np.arange(WIDTH + 1)), xs, side="right")
    s = (WIDTH - bitlength) // 2
    xnorm = xs << (2 * s)
    addresses = (xnorm >> (WIDTH - 2 - r)) - 2**r
    table = np.asarray(design["table"], dtype=np.int64)
    coefs = [np.abs(table[addresses, j]) for j in range(int(design["degree"]) + 1)]
    violations = np.zeros(len(xs), dtype=bool)
    bounds = {}

    def fit(value, width, label, positive=False):
        bad = (value < (1 if positive else 0)) | (value >= 2**width)
        violations[:] |= bad
        bounds[label] = {"bits": width, "minimum_word": int(value.min()),
                         "maximum_word": int(value.max()),
                         "minimum_headroom_words": int((2**width - 1 - value).min()),
                         "capacity_violations": int(bad.sum())}
        return value & (2**width - 1)

    for j in reversed(range(int(design["degree"]))):
        product = np.zeros_like(xs)
        for m in range(r + 1, min(d, WIDTH - 2) + 1):
            product += ((xnorm >> (WIDTH - 2 - m)) & 1) * (coefs[j + 1] >> m)
        coefs[j] = fit(coefs[j] - product, w, f"horner_c{j}", positive=True)
    y0 = coefs[0]
    if design["family"] == "B":
        words = fit(y0 << s, w + N_P, "denormalized_direct")
        return words, w + N_P, violations, bounds

    square = np.zeros_like(xs)
    for k in range(1, w + 1):
        yk = (y0 >> (w - k)) & 1
        if 2 * k <= nrb:
            square += yk << (nrb - 2 * k)
        for ell in range(k + 1, w + 1):
            if k + ell - 1 <= nrb:
                square += (yk & ((y0 >> (w - ell)) & 1)) << (nrb - (k + ell - 1))
    square = fit(square, nrb, "square")
    z = np.zeros_like(xs)
    for i in range(WIDTH):
        e = i - WIDTH + 2  # product_bits == square_bits == nr_bits.
        z += ((xnorm >> i) & 1) * (square << e if e >= 0 else square >> -e)
    z = fit(z, nrb + 2, "multiply_m")
    z = fit((3 << nrb) - z, nrb + 2, "correction")
    p = np.zeros_like(xs)
    for k in range(1, w + 1):
        p += ((y0 >> (w - k)) & 1) * (z >> k)
    p = fit(p, nrb + 1, "multiply_y0")
    words = fit(p << s, nrb + N_P + 1, "denormalized_newton")
    return words, nrb + N_P + 1, violations, bounds


def screen_poly_model(design, verifier, target, tolerance, *, nr_bits=None, chunk_size=65_536):
    """Every valid x, once, in the integer model; never an official circuit PASS."""
    start_time = time.monotonic()
    worst, worst_x, minimum_margin = -1.0, 0, 1 << 64
    failed, overflow_inputs = 0, 0
    aggregate = {}
    first_overflow_x = None
    for start in range(1, MAX_X + 1, chunk_size):
        xs = np.arange(start, min(start + chunk_size, MAX_X + 1), dtype=np.int64)
        words, fractional_bits, overflow, bounds = poly_model_words(design, xs, nr_bits=nr_bits)
        overflow_inputs += int(overflow.sum())
        if overflow.any() and first_overflow_x is None:
            first_overflow_x = int(xs[np.flatnonzero(overflow)[0]])
        for name, values in bounds.items():
            if name not in aggregate:
                aggregate[name] = values
            else:
                prev = aggregate[name]
                for key in ("minimum_word", "minimum_headroom_words"):
                    prev[key] = min(prev[key], values[key])
                prev["maximum_word"] = max(prev["maximum_word"], values["maximum_word"])
                prev["capacity_violations"] += values["capacity_violations"]
        theta = words.astype(np.uint64) << np.uint64(64 - fractional_bits)
        allowed, expected = tolerance[start - 1:start - 1 + len(xs)], target[start - 1:start - 1 + len(xs)]
        used = np.abs((theta - expected).view(np.int64)).view(np.uint64)
        failed += int((used > allowed).sum())
        ratios = used / allowed
        index = int(ratios.argmax())
        if float(ratios[index]) > worst:
            worst, worst_x = float(ratios[index]), int(xs[index])
        minimum_margin = min(minimum_margin,
                             int((allowed.astype(np.int64) - used.astype(np.int64)).min()))
    status = "MODEL_PASS_CIRCUIT_UNVERIFIED" if failed == 0 and overflow_inputs == 0 and minimum_margin > 0 else "MODEL_FAIL_CIRCUIT_UNVERIFIED"
    return {"status": status, "checked_inputs": MAX_X, "failed_phase_inputs": failed,
            "inputs_with_register_capacity_violation": overflow_inputs,
            "first_capacity_violation_x": first_overflow_x,
            "worst_error_over_tolerance": worst, "worst_x": worst_x,
            "minimum_margin_units": minimum_margin,
            "fractional_phase_bits": fractional_bits, "register_ranges": aggregate,
            "elapsed_seconds": time.monotonic() - start_time,
            "limitation": "Exhaustive deterministic integer model only; full circuit and measurement branches remain unverified."}


def simulator_model_samples(circuit, design, verifier, *, nr_bits=None):
    """Small basis-state check of normalization boundaries and compiled arithmetic."""
    inputs = {1, 2, 3, 4, 7, 8, 15, 16, MAX_X}
    for k in range(WIDTH):
        inputs.update(2**k + e for e in (-1, 0, 1) if 1 <= 2**k + e <= MAX_X)
    r = int(design["interval_bits"])
    for i in np.unique(np.linspace(0, 3 * 2**r - 1, 17, dtype=int)):
        boundary = 2**20 + int(i) * 2**(20 - r)
        inputs.update(boundary + e for e in (-1, 0, 1) if 1 <= boundary + e <= MAX_X)
    xs = np.tile(np.asarray(sorted(inputs), dtype=np.uint64), 2)
    sim = km.Simulator(batch_size=len(xs), counting_operations=True, seed=90210)
    xreg = circuit.register_data["x"]
    scratch = km.array([km.q(i) for i in range(circuit.num_qubits) if km.q(i) not in set(xreg)])
    sim.write_across_shots(xreg, verifier.to_bits(xs, WIDTH))
    sim.do(circuit)
    words, f, _, _ = poly_model_words(design, xs, nr_bits=nr_bits)
    expected = words.astype(np.uint64) << np.uint64(64 - f)
    actual = np.asarray([verifier.half_turns_to_units(sim.read_shot_phase(i))
                         for i in range(len(xs))], dtype=np.uint64)
    restored = bool(np.array_equal(verifier.from_bits(sim.read_across_shots(xreg)), xs))
    clean = bool(not sim.read_across_shots(scratch).any())
    matching = bool(np.array_equal(actual, expected))
    measured = sim.operation_statistics()
    static = circuit.operation_statistics()
    counts_agree = (measured.num_toffoli_like == len(xs) * static.num_toffoli_like
                    and measured.num_t_like == len(xs) * static.num_t_like
                    and measured.num_z_pow(min_clifford_hierarchy_level=4)
                    == len(xs) * static.num_z_pow(min_clifford_hierarchy_level=4))
    if not (restored and clean and matching):
        raise AssertionError(f"Circuit/model sample mismatch for {design['name']}")
    return {"status": "SAMPLED_CIRCUIT_MODEL_AGREEMENT", "unique_inputs": len(inputs),
            "basis_shots": len(xs), "seed": 90210, "input_values": sorted(inputs),
            "x_restored": restored, "scratch_clean": clean, "phase_matches_integer_model_exactly": matching,
            "sampled_priced_operation_counts_equal_static_counts": bool(counts_agree),
            "mean_executed_toffoli_like": measured.num_toffoli_like / len(xs),
            "mean_executed_base_t_angles": measured.num_t_like / len(xs),
            "mean_executed_fine_rotations": measured.num_z_pow(min_clifford_hierarchy_level=4) / len(xs),
            "limitation": "A small seeded sample, not the official all-input repeated-shot circuit verification."}


def _prefix(builder, inverse, nr, x, flags, selectors, coefs, square, z, scratch):
    prec = nr.Precision(12, 17, 28, 27, 28)
    steps = [
        ("normalize", lambda: nr._normalize(builder, inverse, x, flags)),
        ("select_interval", lambda: nr._select_interval(builder, inverse, x, selectors)),
        ("load_coefficients", lambda: nr._load_coefficients(builder, prec, x, selectors, coefs)),
        ("horner", lambda: nr._horner(builder, inverse, prec, x, selectors, coefs, scratch)),
        ("square", lambda: nr._square(builder, inverse, coefs[0], square, scratch)),
        ("multiply_m", lambda: nr._multiply_m(builder, inverse, x, square, z, scratch)),
        ("correct", lambda: nr._correct(builder, inverse, z))]
    with builder.mark("uncompute_prefix" if inverse else "compute_prefix"):
        for name, action in reversed(steps) if inverse else steps:
            with builder.mark(name):
                action()


def _exponent_indicator(builder, s, flags, scratch):
    """XOR exact unary exponent equality into a clean scratch (self-inverse)."""
    if s == 0:
        builder.x(scratch)
        builder.cx(flags[0], scratch)
    elif s == N_P:
        builder.cx(flags[N_P - 1], scratch)
    elif 1 <= s < N_P:
        builder.x(flags[s])
        builder.ccx(flags[s - 1], flags[s], scratch)
        builder.x(flags[s])
    else:
        raise ValueError("Exponent outside 0..10")


def build_phase_probe(variant: dict, nr):
    """Propagate exported parity phases after the exact incumbent prefix."""
    builder = km.CircuitBuilder()
    x = builder.create_quantum_register(WIDTH, name="x")
    flags = builder.alloc_qubits(N_P)
    selectors = builder.alloc_qubits(3)
    coefs = [builder.alloc_qubits(17) for _ in range(4)]
    square = builder.alloc_qubits(28)
    z = builder.alloc_qubits(29)
    scratch = builder.alloc_qubits(29)
    owned = [flags, selectors, *coefs, square, z, scratch]
    regs = (x, flags, selectors, coefs, square, z, scratch)
    _prefix(builder, False, nr, *regs)
    indicator = scratch[0]  # Arithmetic scratch is clean when the prefix finishes.
    inputs = [*coefs[0], *z, indicator]
    with builder.mark("phase_polynomial"):
        for block in variant["exponent_blocks"]:
            s = int(block["s"])
            terms = block["terms"]
            has_indicator = any(46 in term["bits"] for term in terms)
            if has_indicator:
                with builder.mark("exponent_indicator"):
                    _exponent_indicator(builder, s, flags, indicator)
            for term in terms:
                bits = [int(i) for i in term["bits"]]
                if not bits or len(set(bits)) != len(bits) or min(bits) < 0 or max(bits) > 46:
                    raise ValueError(f"Invalid parity support {bits}")
                qs = [inputs[i] for i in bits]
                for control in qs[:-1]:
                    builder.cx(control, qs[-1])
                builder.z_pow(qs[-1], 2 * Fraction(str(term["turns"])))
                for control in reversed(qs[:-1]):
                    builder.cx(control, qs[-1])
            if has_indicator:
                with builder.mark("exponent_indicator"):
                    _exponent_indicator(builder, s, flags, indicator)
    _prefix(builder, True, nr, *regs)
    for reg in owned:
        builder.free(reg)
    circuit = builder.finish_circuit()
    return circuit, {"family": "C", "resource_only": True,
                     "explicitly_allocated_qubits": WIDTH + sum(len(reg) for reg in owned),
                     "phase_terms": sum(len(b["terms"]) for b in variant["exponent_blocks"]),
                     "exponent_blocks": len(variant["exponent_blocks"]),
                     "indicator_strategy": "reuse clean arithmetic scratch; exact CCX compute and exact CCX clear",
                     "phase_bit_mapping": "indices0..16=y0 LSB0..16; indices17..45=z LSB0..28; index46=exponent indicator",
                     "limitation": "fixed-s0-optimistic is diagnostic only; all-exponents retains exported terms without synthesis optimization"}


def _records(path: Path, key: str):
    if not path.exists():
        return []
    data = json.loads(path.read_text())
    return data if isinstance(data, list) else data[key]


def describe(circuit, metadata, verifier, margins, bound=None, margin_fraction=None):
    stats = circuit_stats(circuit)
    result = {**metadata, **stats, "tagged_profile": tagged_profile(circuit),
              "objective_margin_sweep": [margin_objective(stats, f, verifier) for f in margins]}
    if bound is not None:
        result["exported_certified_relative_bound"] = str(bound)
        if margin_fraction is None:
            margin_fraction = 1 - Fraction(str(bound)) * 2**24
    if margin_fraction is not None:
        result["objective_at_exported_bound"] = margin_objective(stats, margin_fraction, verifier)
        result["objective_at_exported_bound"]["caveat"] = (
            "Uses an imported approximation/arithmetic bound as a hypothetical synthesis "
            "margin. This resource circuit has not been independently shown to implement that bound.")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True,
                        help="Path to the separately downloaded competition bundle")
    parser.add_argument("--designs", type=Path, default=HERE / "designs.json")
    parser.add_argument("--phase-terms", type=Path, default=HERE / "phase-terms.json")
    parser.add_argument("--output", type=Path, default=HERE / "resource-results.json")
    parser.add_argument("--margin-fractions", nargs="+", default=["1/2", "1/4", "1/10", "1/100"])
    args = parser.parse_args()
    nr, verifier = load_official(args.data)
    started = time.monotonic()
    sample = args.data / "sample_submission/invsqrt.kmb"
    incumbent = HERE.parent / "results/candidate.kmb"
    files = {"sample_submission/invsqrt.kmb": sample,
             "../results/candidate.kmb": incumbent,
             "verifiers/invsqrt/nr_reference_circuit.py": args.data / "verifiers/invsqrt/nr_reference_circuit.py",
             "verifiers/invsqrt/verify_submission_invsqrt.py": args.data / "verifiers/invsqrt/verify_submission_invsqrt.py"}
    result = {"status": "RESOURCE_ONLY", "python": sys.version, "platform": platform.platform(),
              "kickmix_version": km.__version__, "kickmix_module": Path(km.__file__).name,
              "kickmix_binary_sha256": sha256(Path(km.__file__)),
              "probe_source_sha256": sha256(Path(__file__)),
              "source_sha256": {label: sha256(p) for label, p in files.items()},
              "scope": "No production candidate, no full verification, no architecture selected.",
              "reference_and_incumbent": {}, "polynomial_probes": [], "phase_probes": []}
    for name, path, precision in (
            ("official_sample", sample, dataclasses.asdict(nr.Precision())),
            ("incumbent", incumbent, dataclasses.asdict(nr.Precision(12, 17, 28, 27, 28)))):
        c = km.Circuit.from_file(path)
        result["reference_and_incumbent"][name] = describe(
            c, {"file": "sample_submission/invsqrt.kmb" if name == "official_sample" else "../results/candidate.kmb",
                "precision": precision, "resource_only": False},
            verifier, args.margin_fractions)
    designs = _records(args.designs, "designs")
    if designs:
        target, tolerance = verifier.reference_values(
            verifier.ProblemSpecInvSqrt.for_n_p(N_P), np.arange(1, MAX_X + 1, dtype=np.uint64))
    for design in designs:
        widths = [int(design.get("nr_bits", 32))] if design["family"] == "A" else [None]
        if design["family"] == "A" and widths[0] != 28:
            widths.append(28)
        for nrb in widths:
            c, metadata = build_poly_probe(design, nr, nr_bits=nrb)
            metadata["name"] = design["name"] + (f"-nr{nrb}" if nrb is not None else "")
            metadata["imported_design"] = {k: v for k, v in design.items() if k not in ("table", "intervals")}
            bound, margin = design.get("certified_relative_bound"), design.get("margin_fraction")
            if nrb == 28 and nrb != int(design.get("nr_bits", 32)):
                bound, margin = None, None
                metadata["precision_status"] = "UNCERTIFIED_WIDTH_VARIANT; no arithmetic guarantee transferred"
            probe = describe(c, metadata, verifier, args.margin_fractions, bound, margin)
            probe["simulator_model_samples"] = simulator_model_samples(
                c, design, verifier, nr_bits=nrb)
            probe["full_domain_integer_model"] = screen_poly_model(
                design, verifier, target, tolerance, nr_bits=nrb)
            model = probe["full_domain_integer_model"]
            if model["status"] == "MODEL_PASS_CIRCUIT_UNVERIFIED":
                model_fraction = Fraction(model["minimum_margin_units"], int(tolerance[-1]))
                probe["objective_at_model_margin"] = margin_objective(probe, model_fraction, verifier)
                probe["objective_at_model_margin"]["caveat"] = "Model-checked synthesis margin; circuit remains unverified."
            result["polynomial_probes"].append(probe)
            print(json.dumps({"name": probe["name"], "model_status": model["status"],
                              "worst_error_over_tolerance": model["worst_error_over_tolerance"],
                              "minimum_margin_units": model["minimum_margin_units"],
                              "overflow_inputs": model["inputs_with_register_capacity_violation"]}), flush=True)
    for variant in _records(args.phase_terms, "variants"):
        c, metadata = build_phase_probe(variant, nr)
        metadata["name"] = variant["name"]
        metadata["imported_variant"] = {k: v for k, v in variant.items() if k != "exponent_blocks"}
        result["phase_probes"].append(describe(
            c, metadata, verifier, args.margin_fractions,
            variant.get("certified_relative_bound"), variant.get("margin_fraction")))
    result["input_sha256"] = {p.name: sha256(p) for p in (args.designs, args.phase_terms) if p.exists()}
    result["pending_inputs"] = [str(p) for p in (args.designs, args.phase_terms) if not p.exists()]
    result["elapsed_seconds"] = time.monotonic() - started
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    compact = {"output": str(args.output), "pending_inputs": result["pending_inputs"],
               "circuits": [{"name": k, "Q": v["num_qubits"], "Toffoli": v["toffoli_like"],
                             "fine_rotations": v["fine_rotations"]}
                            for k, v in result["reference_and_incumbent"].items()] +
                           [{"name": p["name"], "Q": p["num_qubits"], "Toffoli": p["toffoli_like"],
                             "fine_rotations": p["fine_rotations"]}
                            for p in result["polynomial_probes"] + result["phase_probes"]]}
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
