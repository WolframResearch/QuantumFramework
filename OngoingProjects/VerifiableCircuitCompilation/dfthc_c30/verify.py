"""Run the unmodified official verifier on the standalone C30 family."""

import argparse
import importlib
import json
from pathlib import Path
import platform
import time

import kickmix as km
import numpy as np

from generator import HERE, digest, official_modules


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--submission", type=Path, default=HERE / "local/standalone")
    parser.add_argument("--report", type=Path, default=HERE / "verified.json")
    args = parser.parse_args()
    verifier, _ = official_modules(args.data)
    circuits_module = importlib.import_module("dfthcsos.scoreinfra.circuits")
    build = json.loads((HERE / "build.json").read_text())
    for profile in build["profiles"]:
        if digest(args.submission / f"angle_k_{profile['k']}.kmb") != profile["sha256"]:
            raise ValueError("Circuit hash differs from build evidence.")
    if digest(args.submission / "submission.npz") != build["factorization_sha256"]:
        raise ValueError("Factorization hash differs from build evidence.")
    salt = "qf-standalone-c30-full-v1"
    started = time.monotonic()
    result = verifier.verify(args.submission, salt=salt)
    (args.submission.parent / "official-result.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if result["verified_shots"] != 3976 * 64 or result["num_circuits"] != 25:
        raise AssertionError("Unexpected full verification coverage.")
    prior = HERE.parent
    baseline = json.loads((prior / "other_tasks/dfthcsos-baseline.json").read_text())["score"]
    incumbent = json.loads((prior / "other_tasks/mixing-results.json").read_text())["full_official_best_result"]["score"]
    invsqrt = 3.9385957568764147
    dqi = (json.loads((prior / "other_tasks/dqi-baseline.json").read_text())["score"] /
           json.loads((prior / "dqi_gf8/verified.json").read_text())["score"])
    result.update({"salt": salt, "shots_per_address": circuits_module.SHOTS_PER_ADDRESS,
                   "passes": ["forward", "reverse"], "official_verifier_modified": False,
                   "python": platform.python_version(), "numpy": np.__version__, "kickmix": km.__version__,
                   "cost_improvement_over_baseline": baseline / result["score"],
                   "cost_improvement_over_incumbent": incumbent / result["score"],
                   "cost_reduction_percent_over_incumbent": 100 * (1 - result["score"] / incumbent),
                   "combined_local_estimate": (invsqrt * dqi * baseline / result["score"]) ** (1 / 3),
                   "new_kaggle_submission": False,
                   "build_sha256": digest(HERE / "build.json"),
                   "source_sha256": {p.name: digest(p) for p in (HERE / "generator.py", HERE / "loader.py", Path(__file__))},
                   "official_dependency_sha256": {str(p.relative_to(args.data)): digest(p) for p in
                       (args.data / "verifiers/dfthcsos/verify.py", args.data / "verifiers/dfthcsos/scoreinfra/circuits.py")},
                   "elapsed_seconds": time.monotonic() - started})
    args.report.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
