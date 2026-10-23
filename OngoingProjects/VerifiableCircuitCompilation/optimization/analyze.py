"""Exhaustive integer screening and compiled circuit comparison, not a full PASS."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import kickmix as km
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import optimized_quadratic as opt
from verify_quadratic import load_verifier, VERIFIER_SHA256

MAX_X = 3_133_452


def model(xs, widths):
    xs = np.asarray(xs, dtype=np.int64)
    if not ((1 <= xs) & (xs <= MAX_X)).all():
        raise ValueError("Inputs outside the benchmark domain.")
    shifts = (22 - np.searchsorted(1 << np.arange(23), xs, side="right")) // 2
    normalized = xs << (2*shifts)
    addresses = (normalized >> 14) - 64
    rows = np.asarray(opt.quantized_rows(widths), dtype=np.int64)[addresses]
    values = [rows[:, j].copy() for j in range(3)]
    lengths = [max(row[j] for row in opt.quantized_rows(widths)).bit_length()
               for j in range(3)]
    violations = np.zeros(len(xs), dtype=bool)
    ranges = {}
    for j in (1, 0):
        # Mathematical word arithmetic, with no quantum register manipulation.
        summands = []
        for m in range(7, 21):
            shift = widths[j] - widths[j+1] - m
            row = values[j+1] << shift if shift >= 0 else values[j+1] >> -shift
            summands.append(((normalized >> (20-m)) & 1) * row)
        values[j] -= np.sum(summands, axis=0)
        bad = (values[j] <= 0) | (values[j] >= 2**lengths[j])
        violations |= bad
        ranges[f"c{j}"] = [int(values[j].min()), int(values[j].max())]
    return values[0] << shifts, widths[0]+10, violations, ranges


def screen(widths, target, tolerance):
    worst, worst_x, margin, failures, capacity = -1., None, 1 << 64, 0, 0
    for start in range(1, MAX_X+1, 65_536):
        xs = np.arange(start, min(start+65_536, MAX_X+1))
        word, f, bad, _ = model(xs, widths)
        phase = word.astype(np.uint64) << np.uint64(64-f)
        expected, allowed = target[start-1:start-1+len(xs)], tolerance[start-1:start-1+len(xs)]
        used = np.abs((phase-expected).view(np.int64)).astype(np.uint64)
        ratios = used/allowed
        k = int(ratios.argmax())
        if ratios[k] > worst:
            worst, worst_x = float(ratios[k]), int(xs[k])
        margin = min(margin, int((allowed.astype(np.int64)-used.astype(np.int64)).min()))
        failures += int((used > allowed).sum())
        capacity += int(bad.sum())
    return {"status": "MODEL_PASS_CIRCUIT_UNVERIFIED" if not (failures or capacity) and margin>0
            else "MODEL_FAIL_CIRCUIT_UNVERIFIED", "checked_inputs": MAX_X,
            "failed_phase_inputs": failures, "capacity_violations": capacity,
            "worst_error_over_tolerance": worst, "worst_x": worst_x,
            "minimum_margin_units": margin}


def sample(circuit, widths, verifier):
    # All reachable interval transitions, normalization transitions and endpoints.
    inputs = {1, MAX_X}
    for exponent in range(11):
        for address in range(64, 257):
            boundary = (address << 14) >> (2*exponent)
            inputs.update(boundary+e for e in (-1, 0, 1)
                          if 1 <= boundary+e <= MAX_X)
    xs = np.tile(np.array(sorted(inputs), dtype=np.uint64), 4)
    sim = km.Simulator(batch_size=len(xs), counting_operations=True, seed=371811)
    xreg, scratch = verifier.validate_registers(circuit, verifier.ProblemSpecInvSqrt.for_n_p(10))
    sim.write_across_shots(xreg, verifier.to_bits(xs, 22))
    sim.do(circuit)
    words, f, bad, _ = model(xs, widths)
    expected = words.astype(np.uint64) << np.uint64(64-f)
    actual = np.array([verifier.half_turns_to_units(sim.read_shot_phase(k))
                       for k in range(len(xs))], dtype=np.uint64)
    restored = np.array_equal(verifier.from_bits(sim.read_across_shots(xreg)), xs)
    clean = not sim.read_across_shots(scratch).any()
    matching = np.array_equal(actual, expected)
    measured, static = sim.operation_statistics(), circuit.operation_statistics()
    counts_agree = all(a == len(xs)*b for a, b in (
        (measured.num_toffoli_like, static.num_toffoli_like),
        (measured.num_t_like, static.num_t_like),
        (measured.num_z_pow(min_clifford_hierarchy_level=4), static.num_z_pow(min_clifford_hierarchy_level=4))))
    if not (restored and clean and matching and not bad.any() and counts_agree):
        raise AssertionError("Compiled circuit disagrees with model, cleanup or resource counts.")
    return {"unique_inputs": len(inputs), "basis_shots": len(xs), "seed": 371811,
            "x_restored": bool(restored), "scratch_clean": bool(clean),
            "phase_matches_integer_model_exactly": bool(matching),
            "sampled_priced_counts_equal_static_counts": counts_agree}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    args = parser.parse_args()
    verifier, _ = load_verifier(args.data)
    target, tolerance = verifier.reference_values(verifier.ProblemSpecInvSqrt.for_n_p(10),
                                                 np.arange(1, MAX_X+1, dtype=np.uint64))
    certificates = json.loads((Path(__file__).parent/"width-certificates.json").read_text())
    records = []
    for certificate in certificates["certificates"]:
        widths = tuple(certificate["fractional_bits"])
        if not certificate["certified"]:
            records.append({"fractional_bits": widths, "status": "SYMBOLIC_CERTIFICATE_REJECTED"})
            continue
        screening = screen(widths, target, tolerance)
        if screening["status"] != "MODEL_PASS_CIRCUIT_UNVERIFIED":
            raise AssertionError("Symbolically certified widths failed the integer model.")
        variants = []
        for arithmetic, recycle in (("padded", False), ("padded", True), ("controlled", True)):
            for budget in (None, 140, 150):
                try:
                    circuit = opt.build_circuit(widths, arithmetic=arithmetic,
                        recycle_coefficient=recycle, adder_budget=budget)
                except ValueError:
                    continue
                stats = circuit.operation_statistics()
                rotations = stats.num_z_pow(min_clifford_hierarchy_level=4)
                eps = math.pi*screening["minimum_margin_units"]/(2**64*rotations)
                cost = circuit.num_qubits*(stats.num_toffoli_like +
                    (stats.num_t_like + rotations*verifier.rotation_t_count(eps))/4)
                variants.append({"arithmetic": arithmetic,"recycle_coefficient": recycle,
                    "adder_budget": budget,"num_qubits": circuit.num_qubits,
                    "toffoli_count": stats.num_toffoli_like,"rotation_count": rotations,
                    "estimated_cost_from_model_margin": cost})
        records.append({"fractional_bits": widths,"screening": screening,
                        "variants": sorted(variants, key=lambda r:r["estimated_cost_from_model_margin"])})
    compiled = opt.build_circuit()
    samples = sample(compiled, opt.FRACTIONAL_BITS, verifier)
    report = {"status": "SCREENING_ONLY_NOT_FULL_VERIFICATION", "verifier_sha256": VERIFIER_SHA256,
              "analysis_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "width_certificates_sha256": hashlib.sha256((Path(__file__).parent/"width-certificates.json").read_bytes()).hexdigest(),
              "candidates": records,"selected_fractional_bits": opt.FRACTIONAL_BITS,
              "selected_sample": samples}
    (Path(__file__).parent/"screening.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"best_each": [{"widths": r["fractional_bits"], "best": r["variants"][0]}
                        for r in records if "variants" in r], "sample": samples}, indent=2))


if __name__ == "__main__":
    main()
