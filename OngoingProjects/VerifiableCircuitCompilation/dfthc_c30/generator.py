"""Build C30 circuits from the hash-pinned public factors, without training."""

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys
import time

import numpy as np

from loader import batch_deltas, build

HERE = Path(__file__).resolve().parent
PUBLIC_SHA256 = "9afd0ca5e90dd1740fe713e4ac1c26736184ea9bc4d4b57a6eeefafd608ee193"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_digest(array):
    a = np.ascontiguousarray(array)
    return hashlib.sha256(a.dtype.str.encode() + str(a.shape).encode() + a.tobytes()).hexdigest()


def official_modules(data):
    pinned = json.loads((HERE.parent / "dfthc_feasibility/inputs.json").read_text())["input_sha256"]
    for relative, expected in pinned.items():
        if digest(data / relative) != expected:
            raise ValueError("Competition dependency hash mismatch: " + relative)
    sys.path.insert(0, str(data / "verifiers"))
    return (importlib.import_module("dfthcsos.verify"),
            importlib.import_module("dfthcsos.scoreinfra.dfthc_math"))


def signed_gauge(words, tree, maths):
    """Released gauge: change only vector signs, with a checked 1e-11 drift cap."""
    vectors = maths.reconstruct_u_vectors(words / 512, tree)
    reduced = vectors.copy()
    angles = np.empty_like(words, dtype=np.float64)
    for j, (p, q) in reversed(list(enumerate(tree))):
        norm = np.hypot(reduced[:, p], reduced[:, q])
        sign = np.where(np.abs(reduced[:, p]) < 1e-12,
                        np.where(reduced[:, q] > 0, -1, 1),
                        np.where(reduced[:, p] < 0, -1, 1))
        angles[:, j] = np.where(norm < 1e-12, 0,
                               np.arctan2(reduced[:, q] * sign, reduced[:, p] * sign))
        reduced[:, p] = np.where(norm < 1e-12, 0, norm * sign)
        reduced[:, q] = 0
    signed = np.rint(angles * (512 / np.pi)).astype(np.int64) % 1024
    decoded = maths.reconstruct_u_vectors(signed / 512, tree)
    signs = np.where(np.sum(vectors * decoded, axis=1) < 0, -1, 1)
    drift = float(np.max(np.abs(decoded * signs[:, None] - vectors)))
    if drift > 1e-11 or np.any(((signed >> 9) ^ (signed >> 8)) & 1):
        raise ArithmeticError("Released sign gauge failed its invariant.")
    return signed, decoded, {"max_vector_drift_up_to_sign": drift,
                             "max_unit_norm_error": float(np.max(np.abs(np.sum(decoded**2, axis=1) - 1))),
                             "only_row_sign_gauge": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--public-factors", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "local/standalone")
    parser.add_argument("--report", type=Path, default=HERE / "build.json")
    args = parser.parse_args()
    started = time.monotonic()
    if digest(args.public_factors) != PUBLIC_SHA256:
        raise ValueError("Public factor file changed from the reviewed release.")
    verifier, maths = official_modules(args.data)
    with np.load(args.data / "sample_submission/dfthcsos/submission.npz", allow_pickle=False) as z:
        baseline = {k: z[k].copy() for k in z.files}
    with np.load(args.public_factors, allow_pickle=False) as z:
        words, weights, h2 = z["words"].astype(np.int64), z["W"].copy(), z["H2_I"].copy()
    signed, vectors, gauge = signed_gauge(words, baseline["tree_edges"], maths)
    params = {"U": vectors[76:].reshape(65, 60, 76), "W": weights, "H2_I": h2,
              "H1_I": baseline["H1_I"], "Hsym_lower": baseline["Hsym_lower"],
              "tree_edges": baseline["tree_edges"], "angle_words": signed}
    inputs = verifier.load_inputs(params)
    metrics = verifier.evaluate_words(inputs, signed, params["tree_edges"], 10)
    if not verifier.passes_accuracy(metrics):
        raise ValueError("Fixed public factors failed official metric gates.")
    if args.output.exists() and list(args.output.iterdir()):
        raise FileExistsError("Build into a new or empty output directory.")
    args.output.mkdir(parents=True, exist_ok=True)
    np.savez(args.output / "submission.npz", **params)
    signed.astype("<u2").tofile(args.output.parent / "words.u16")
    profiles = []
    for k, delta in enumerate(batch_deltas(signed, 3)):
        circuit, profile = build(delta, return_report=True)
        path = args.output / f"angle_k_{k}.kmb"
        circuit.to_file(path, format="kmb")
        profiles.append({"k": k, **profile, "sha256": digest(path), "bytes": path.stat().st_size,
                         "num_ops": len(circuit), "reaction_depth": circuit.reaction_depth()})
    report = {"status": "BUILT_NOT_FULLY_VERIFIED", "public_factor_sha256": PUBLIC_SHA256,
              "factor_shape": [65, 60, 30], "factor_training": False,
              "fixed_array_sha256": {k: array_digest(params[k]) for k in
                                     ("W", "H2_I", "H1_I", "Hsym_lower", "tree_edges")},
              "weights_and_offsets_exactly_released": np.array_equal(params["W"], weights) and np.array_equal(params["H2_I"], h2),
              "gauge": gauge, "metrics": metrics, "metric_gates_pass": True,
              "words_sha256": digest(args.output.parent / "words.u16"),
              "factorization_sha256": digest(args.output / "submission.npz"),
              "circuit_count": len(profiles), "profiles": profiles,
              "source_sha256": {p.name: digest(p) for p in (HERE / "loader.py", Path(__file__))},
              "elapsed_seconds": time.monotonic() - started}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "profiles"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
