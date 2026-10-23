"""Run the unmodified, default full verifier for the reduced quadratic circuit."""

import argparse
import json
from pathlib import Path
import platform
import sys
import time

import kickmix as km
import numpy as np

import optimized_quadratic as opt
from quadratic import DESIGN_SHA256, HERE
from verify_quadratic import SALT, VERIFIER_SHA256, file_hash, load_comparison, load_verifier


def evidence_hashes():
    names = ["quadratic.py", "optimized_quadratic.py", "verify_quadratic.py",
             "verify_optimized_quadratic.py", "feasibility/B-quadratic-192.json",
             "optimization/WidthOptimization.wl", "optimization/run-widths.wls",
             "optimization/verify-widths.wls", "optimization/width-certificates.json",
             "optimization/wl-checks.json", "optimization/analyze.py", "optimization/screening.json"]
    return {name: file_hash(HERE/name) for name in names}


def check_evidence():
    certificates = json.loads((HERE/"optimization/width-certificates.json").read_text())
    selected = next(c for c in certificates["certificates"]
                    if c["fractional_bits"] == list(opt.FRACTIONAL_BITS))
    if (not selected["certified"] or len(selected["intervals"]) != 192 or
            any(c["phase"] is not True or c["capacity"] is not True for c in selected["intervals"]) or
            selected["table"] != [list(row) for row in opt.quantized_rows()] or
            certificates["source_sha256"] != DESIGN_SHA256 or
            certificates["certificate_source_sha256"] != file_hash(HERE/"optimization/WidthOptimization.wl")):
        raise ValueError("Symbolic certificate does not match the generator.")
    checks = json.loads((HERE/"optimization/wl-checks.json").read_text())
    if (checks["failed"] != 0 or checks["succeeded"] != checks["tests"] or
            checks["verifier_source_sha256"] != file_hash(HERE/"optimization/verify-widths.wls")):
        raise ValueError("WL checks are failing or stale.")
    screen = json.loads((HERE/"optimization/screening.json").read_text())
    selected_screen = next(c for c in screen["candidates"]
                           if c["fractional_bits"] == list(opt.FRACTIONAL_BITS))
    if (screen["analysis_source_sha256"] != file_hash(HERE/"optimization/analyze.py") or
            screen["width_certificates_sha256"] != file_hash(HERE/"optimization/width-certificates.json") or
            screen["selected_fractional_bits"] != list(opt.FRACTIONAL_BITS) or
            selected_screen["screening"]["status"] != "MODEL_PASS_CIRCUIT_UNVERIFIED" or
            selected_screen["screening"]["checked_inputs"] != 3_133_452 or
            selected_screen["screening"]["minimum_margin_units"] <= 0 or
            not all(screen["selected_sample"][k] is True for k in (
                "x_restored","scratch_clean","phase_matches_integer_model_exactly",
                "sampled_priced_counts_equal_static_counts"))):
        raise ValueError("Integer screening evidence is stale.")
    return evidence_hashes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE/"results")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    inputs = check_evidence()
    verifier, verifier_source = load_verifier(args.data)
    baseline, _ = load_comparison(args.data)
    incumbent = json.loads((HERE/"results/quadratic-verified.json").read_text())
    if (incumbent["status"] != "PASS" or incumbent["verified_shots"] != 200_540_928 or
            incumbent["verifier_sha256"] != VERIFIER_SHA256 or incumbent["salt"] != SALT or
            incumbent["build"]["generator_sha256"] != file_hash(HERE/"quadratic.py")):
        raise ValueError("The quadratic incumbent report does not match the comparison protocol.")
    spec = verifier.ProblemSpecInvSqrt.for_n_p(10)
    if (spec.width, spec.max_x, spec.precision_bits, spec.shots_per_x, spec.num_shots) != (
            22, 3_133_452, 24, 64, 200_540_928):
        raise ValueError("The full benchmark specification has changed.")
    args.output.mkdir(parents=True, exist_ok=True)
    artifact = args.output/"optimized-quadratic.kmb"
    manifest = opt.build_to_file(artifact)
    artifact.with_suffix(".build.json").write_text(json.dumps(manifest, indent=2)+"\n")
    circuit = km.Circuit.from_file(artifact)
    report_path = args.output/"optimized-quadratic-verified.json"
    metadata = {"design": manifest["design"], "build": manifest, "input_sha256": inputs,
                "verifier_sha256": VERIFIER_SHA256, "salt": SALT,
                "python": sys.version, "platform": platform.platform(), "numpy": np.__version__,
                "official_verifier_modified": False,
                "protocol": {"width": spec.width, "max_x": spec.max_x,
                             "precision_bits": spec.precision_bits,"shots_per_x": spec.shots_per_x,
                             "verified_shots": spec.num_shots, "beta": str(spec.beta),
                             "batch_size": verifier.DEFAULT_BATCH_SIZE, "max_workers": args.workers}}
    report_path.write_text(json.dumps({**metadata,"status":"RUNNING"}, indent=2)+"\n")
    print(json.dumps({"status":"RUNNING","shots":spec.num_shots,"build":manifest}), flush=True)
    started = time.monotonic()
    try:
        result = verifier.verify(circuit, spec, salt=SALT, max_workers=args.workers)
        if (inputs != evidence_hashes() or file_hash(verifier_source) != VERIFIER_SHA256 or
                file_hash(artifact) != manifest["circuit_sha256"] or
                result["circuit_sha256"] != manifest["circuit_sha256"]):
            raise RuntimeError("Verification input changed or circuit hash disagrees.")
        if result["status"] != "PASS" or result["verified_shots"] != spec.num_shots:
            raise RuntimeError("The full verifier did not pass every required shot.")
        actual = result["score"]
        if not actual < incumbent["score"]:
            raise RuntimeError("The verified candidate does not improve the incumbent.")
        comparison = {"baseline_cost":baseline["score"], "quadratic_incumbent_cost":incumbent["score"],
                      "baseline_improvement":baseline["score"]/actual,
                      "quadratic_incumbent_improvement":incumbent["score"]/actual,
                      "baseline_cost_reduction_percent":100*(1-actual/baseline["score"]),
                      "quadratic_incumbent_cost_reduction_percent":100*(1-actual/incumbent["score"]),
                      "baseline_circuit_sha256":baseline["circuit_sha256"],
                      "quadratic_incumbent_circuit_sha256":incumbent["circuit_sha256"],
                      "aggregate_if_other_tracks_at_baseline":(baseline["score"]/actual)**(1/3)}
        report = {**metadata, **result, "comparison": comparison,
                  "elapsed_seconds": time.monotonic()-started}
        report_path.write_text(json.dumps(report, indent=2)+"\n")
        print(json.dumps(report, indent=2), flush=True)
    except Exception as error:
        report_path.write_text(json.dumps({**metadata,"status":"FAIL",
            "failure_type":type(error).__name__,"elapsed_seconds":time.monotonic()-started}, indent=2)+"\n")
        raise


if __name__ == "__main__":
    main()
