"""Run the unmodified official full verifier for the frozen quadratic circuit."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import platform
import sys
import time

import kickmix as km
import numpy as np

from quadratic import HERE, DESIGN, DESIGN_SHA256, build_to_file

SALT = "qf-prototype-local-validation-v1"
VERIFIER_SHA256 = "ddef520dec6b68d1cf0a6f328ddd36d79eb84570b9192cf3f7acaf073f06733a"


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_verifier(data):
    source = data / "verifiers/invsqrt/verify_submission_invsqrt.py"
    if file_hash(source) != VERIFIER_SHA256:
        raise ValueError("The downloaded official verifier differs from the reviewed version.")
    sys.path.insert(0, str(data / "verifiers"))
    verifier = importlib.import_module("invsqrt.verify_submission_invsqrt")
    if Path(verifier.__file__).resolve() != source.resolve():
        raise ValueError("An unexpected verifier module was already imported.")
    return verifier, source


def load_comparison(data):
    baseline = json.loads((HERE / "results/baseline-verified.json").read_text())
    incumbent = json.loads((HERE / "results/candidate-verified.json").read_text())
    baseline_hash = file_hash(data / "sample_submission/invsqrt.kmb")
    for record in (baseline, incumbent):
        if (record["status"] != "PASS" or record["verified_shots"] != 200_540_928
                or record["data_sha256"]["verifiers/invsqrt/verify_submission_invsqrt.py"] != VERIFIER_SHA256
                or record["salt"] != SALT):
            raise ValueError("Comparison report does not match the full verification protocol.")
    if baseline["circuit_sha256"] != baseline_hash:
        raise ValueError("The supplied baseline circuit has changed.")
    return baseline, incumbent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "results")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    verifier, verifier_source = load_verifier(args.data)
    baseline, incumbent = load_comparison(args.data)
    spec = verifier.ProblemSpecInvSqrt.for_n_p(10)
    if (spec.width, spec.max_x, spec.precision_bits, spec.shots_per_x, spec.num_shots) != (
            22, 3_133_452, 24, 64, 200_540_928):
        raise ValueError("The full benchmark specification has changed.")
    args.output.mkdir(parents=True, exist_ok=True)
    artifact = args.output / "quadratic.kmb"
    manifest = build_to_file(artifact)
    artifact.with_suffix(".build.json").write_text(json.dumps(manifest, indent=2) + "\n")
    circuit = km.Circuit.from_file(artifact)
    report_path = args.output / "quadratic-verified.json"
    metadata = {"design": manifest["design"], "build": manifest,
                "verifier_sha256": VERIFIER_SHA256,
                "runner_sha256": file_hash(Path(__file__)),
                "python": sys.version, "platform": platform.platform(),
                "numpy": np.__version__, "salt": SALT,
                "protocol": {"width": spec.width, "max_x": spec.max_x,
                             "precision_bits": spec.precision_bits,
                             "shots_per_x": spec.shots_per_x, "verified_shots": spec.num_shots,
                             "beta": str(spec.beta), "batch_size": verifier.DEFAULT_BATCH_SIZE,
                             "max_workers": args.workers},
                "official_verifier_modified": False}
    report_path.write_text(json.dumps({**metadata, "status": "RUNNING"}, indent=2) + "\n")
    print(json.dumps({"status": "RUNNING", "shots": spec.num_shots, "build": manifest}), flush=True)
    started = time.monotonic()
    try:
        result = verifier.verify(circuit, spec, salt=SALT, max_workers=args.workers)
        if (file_hash(DESIGN) != DESIGN_SHA256 or file_hash(verifier_source) != VERIFIER_SHA256
                or file_hash(artifact) != manifest["circuit_sha256"]
                or file_hash(HERE / "quadratic.py") != manifest["generator_sha256"]
                or file_hash(Path(__file__)) != metadata["runner_sha256"]):
            raise RuntimeError("A verification input or executable changed during the run.")
        if result["circuit_sha256"] != manifest["circuit_sha256"]:
            raise RuntimeError("Verifier circuit hash differs from the generated artifact.")
        actual = result["score"]
        comparison = {"baseline_cost": baseline["score"], "incumbent_cost": incumbent["score"],
                      "baseline_improvement": baseline["score"] / actual,
                      "incumbent_improvement": incumbent["score"] / actual,
                      "baseline_cost_reduction_percent": 100 * (1 - actual / baseline["score"]),
                      "incumbent_cost_reduction_percent": 100 * (1 - actual / incumbent["score"]),
                      "aggregate_if_other_tracks_at_baseline": (baseline["score"] / actual)**(1/3),
                      "baseline_circuit_sha256": baseline["circuit_sha256"],
                      "incumbent_circuit_sha256": incumbent["circuit_sha256"]}
        report = {**metadata, **result, "comparison": comparison,
                  "elapsed_seconds": time.monotonic() - started}
        report_path.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2), flush=True)
    except Exception as error:
        report_path.write_text(json.dumps({**metadata, "status": "FAIL",
            "failure_type": type(error).__name__, "elapsed_seconds": time.monotonic()-started}, indent=2) + "\n")
        raise


if __name__ == "__main__":
    main()
