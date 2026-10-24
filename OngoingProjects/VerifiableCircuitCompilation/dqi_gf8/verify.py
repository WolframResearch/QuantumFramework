"""Run the unmodified full official DQI verifier on the generated circuit."""

import argparse
import importlib
import json
from pathlib import Path
import sys
import time

import kickmix as km

if __package__:
    from .gf8_backend import digest
else:
    from gf8_backend import digest

HERE = Path(__file__).resolve().parent
SALT = "qf-other-tasks-feasibility-v1"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--circuit", type=Path, default=HERE / "local/dqi-gf8.kmb")
    parser.add_argument("--build", type=Path, default=HERE / "build.json")
    parser.add_argument("--report", type=Path, default=HERE / "verified.json")
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    baseline_path = HERE.parent / "other_tasks/dqi-baseline.json"
    baseline = json.loads(baseline_path.read_text())
    build = json.loads(args.build.read_text())
    if (build["m"], build["n"], build["degree"]) != (4095, 70, 12):
        raise ValueError("The full report requires the unchanged benchmark parameters.")
    for name, expected in build["source_sha256"].items():
        if digest(HERE / name) != expected:
            raise ValueError("The built source binding changed: " + name)
    if digest(args.circuit) != build["circuit_sha256"]:
        raise ValueError("The candidate circuit differs from the build record.")
    inputs = baseline["input_sha256"]
    for name, expected in inputs.items():
        if digest(args.data / name) != expected:
            raise ValueError("An official input changed from the baseline: " + name)
    bindings = {p.name: digest(p) for p in
                (Path(__file__), HERE / "generator.py", HERE / "gf8_backend.py",
                 HERE / "run-certificates.wls", HERE / "certificates.json", args.build)}
    metadata = {"coverage": "DEFAULT_FULL_OFFICIAL_DQI", "salt": SALT,
                "official_code_modified": False, "input_sha256": inputs,
                "source_sha256": bindings, "baseline_sha256": digest(baseline_path),
                "kickmix_version": km.__version__, "workers": args.workers}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps({**metadata, "status": "RUNNING"}, indent=2) + "\n")
    print(json.dumps({"status": "RUNNING", "qubits": build["num_qubits"],
                      "static_toffoli_count": build["static_toffoli_count"]}), flush=True)
    started = time.monotonic()
    try:
        sys.path.insert(0, str(args.data / "verifiers"))
        official = importlib.import_module("dqi.verify_submission_dqi")
        result = official.verify(km.Circuit.from_file(args.circuit),
                                 official.ProblemSpecDQI.benchmark_instance(),
                                 salt=SALT, max_workers=args.workers)
        for name, expected in inputs.items():
            if digest(args.data / name) != expected:
                raise RuntimeError("An official input changed during verification.")
        for name, expected in bindings.items():
            path = args.build if name == args.build.name else HERE / name
            if digest(path) != expected:
                raise RuntimeError("A bound source changed during verification.")
        if result["circuit_sha256"] != build["circuit_sha256"]:
            raise RuntimeError("The verifier's canonical circuit hash differs.")
        improvement = baseline["score"] / result["score"]
        report = {**metadata, **result,
                  "comparison": {"baseline_cost": baseline["score"],
                                 "improvement": improvement,
                                 "cost_reduction_percent": 100 * (1 - 1 / improvement),
                                 "baseline_qubits": baseline["num_qubits"],
                                 "baseline_mean_toffoli": baseline["toffoli_count"]},
                  "elapsed_seconds": time.monotonic() - started}
        args.report.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2), flush=True)
    except Exception as error:
        args.report.write_text(json.dumps({**metadata, "status": "FAIL",
                              "failure_type": type(error).__name__,
                              "elapsed_seconds": time.monotonic() - started}, indent=2) + "\n")
        raise


if __name__ == "__main__":
    main()
