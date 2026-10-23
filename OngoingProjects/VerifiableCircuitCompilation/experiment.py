"""Local inverse-square-root experiments using the unmodified competition verifier.

Competition data stays outside this repository. Model screening is explicitly separate
from circuit verification; only verify() can produce an official-contract PASS.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import sys
import time

import kickmix as km
import numpy as np

HERE = Path(__file__).resolve().parent
MAX_X = 3_133_452


def load_competition(data: Path):
    sys.path.insert(0, str(data / "verifiers"))
    from invsqrt import nr_reference_circuit as nr
    from invsqrt import verify_submission_invsqrt as verifier
    return nr, verifier


def provenance(data: Path) -> dict:
    paths = ["sample_submission/invsqrt.kmb",
             "verifiers/invsqrt/nr_reference_circuit.py",
             "verifiers/invsqrt/verify_submission_invsqrt.py"]
    return {"python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "kickmix": km.__version__,
            "data_sha256": {p: hashlib.sha256((data / p).read_bytes()).hexdigest()
                            for p in paths}}


def save(path: Path, result: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + "\n")


def screen_full_domain(nr, verifier, precision, chunk_size=65_536) -> dict:
    """Integer phase/tolerance comparison at every valid x, before circuit simulation."""
    spec = verifier.ProblemSpecInvSqrt.for_n_p(10)
    minimum_margin, worst_ratio, worst_x, failed = 1 << 64, -1.0, 0, 0
    fractional_bits = precision.output_bits(10)
    for start in range(1, MAX_X + 1, chunk_size):
        xs = np.arange(start, min(start + chunk_size, MAX_X + 1), dtype=np.uint64)
        words = nr.fixed_point_model(10, precision, xs).astype(np.uint64)
        # Our searched widths all have F <= 40, so this shift exactly represents
        # the dyadic model phase in the verifier's uint64 units, modulo one turn.
        theta = words << np.uint64(64 - fractional_bits)
        target, tolerance = verifier.reference_values(spec, xs)
        used = np.abs((theta - target).view(np.int64)).view(np.uint64)
        bad = used > tolerance
        failed += int(bad.sum())
        ratios = used / tolerance
        index = int(ratios.argmax())
        if ratios[index] > worst_ratio:
            worst_ratio, worst_x = float(ratios[index]), int(xs[index])
        # Signed Python integers avoid uint64 underflow on rejected candidates.
        minimum_margin = min(minimum_margin,
                             int((tolerance.astype(np.int64) - used.astype(np.int64)).min()))
    return {"status": "MODEL_PASS" if failed == 0 and minimum_margin > 0 else "MODEL_FAIL",
            "checked_inputs": MAX_X, "failed_inputs": failed,
            "worst_error_over_tolerance": worst_ratio, "worst_x": worst_x,
            "minimum_margin_units": minimum_margin}


def estimated_cost(circuit, margin_units: int, verifier) -> dict:
    stats = circuit.operation_statistics()
    rotations = stats.num_z_pow(min_clifford_hierarchy_level=4)
    eps = math.pi * margin_units / 2**64 / rotations if rotations else None
    t_synthesis = rotations * verifier.rotation_t_count(eps) if rotations else 0
    c = stats.num_toffoli_like + (stats.num_t_like + t_synthesis) / 4
    return {"num_qubits": circuit.num_qubits, "toffoli_count": stats.num_toffoli_like,
            "rotation_count": rotations, "toffoli_equivalent": c,
            "model_estimated_score": circuit.num_qubits * c}


def search(nr, verifier, output: Path) -> dict:
    sample = np.unique(np.r_[np.arange(1, 1025),
                            np.linspace(1, MAX_X, 8192, dtype=np.int64),
                            2**np.arange(22), np.maximum(1, 2**np.arange(22) - 1)])
    survivors = []
    tested = 0
    for widths in itertools.product(range(10, 15), range(14, 18),
                                    range(25, 29), range(25, 29), range(26, 30)):
        precision = nr.Precision(*widths)
        tested += 1
        try:
            error = nr.relative_error(10, precision, sample)
            if float(error.max()) >= 2**-24:
                continue
            circuit = nr.build_nr_circuit(10, prec=precision)
        except (OverflowError, ValueError):
            continue
        # Gate cost is a first ordering only. Full-domain margin determines the
        # actual rotation-synthesis cost below.
        stat = circuit.operation_statistics()
        proxy = circuit.num_qubits * (stat.num_toffoli_like + 8 * stat.num_z_pow(min_clifford_hierarchy_level=4))
        survivors.append((proxy, precision))
    full_results = []
    best = None
    for _, precision in sorted(survivors, key=lambda p: p[0])[:24]:
        result = screen_full_domain(nr, verifier, precision)
        result["precision"] = dataclasses.asdict(precision)
        if result["status"] == "MODEL_PASS":
            circuit = nr.build_nr_circuit(10, prec=precision)
            result.update(estimated_cost(circuit, result["minimum_margin_units"], verifier))
            if best is None or result["model_estimated_score"] < best["model_estimated_score"]:
                best = result
                circuit.to_file(output / "candidate.kmb", format="kmb")
        full_results.append(result)
        save(output / "search.json", {"sampled_configurations": tested,
             "sample_survivors": len(survivors), "fully_checked": full_results, "best": best})
    if best is None:
        raise RuntimeError("No candidate passed full-domain model screening.")
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["baseline", "search", "candidate"])
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=HERE / "results")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--salt", default="qf-prototype-local-validation-v1")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    nr, verifier = load_competition(args.data)
    started = time.monotonic()
    metadata = provenance(args.data)
    if args.mode == "search":
        best = search(nr, verifier, args.output)
        save(args.output / "candidate-model.json", {**metadata, **best,
             "elapsed_seconds": time.monotonic() - started})
        print(json.dumps(best, indent=2), flush=True)
        return
    path = (args.data / "sample_submission/invsqrt.kmb" if args.mode == "baseline"
            else args.output / "candidate.kmb")
    circuit = km.Circuit.from_file(path)
    result = verifier.verify(circuit, verifier.ProblemSpecInvSqrt.for_n_p(10),
                             salt=args.salt, max_workers=args.workers)
    save(args.output / f"{args.mode}-verified.json", {**metadata, **result,
         "salt": args.salt, "elapsed_seconds": time.monotonic() - started})
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
