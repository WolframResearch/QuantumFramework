"""Build the external reference DQI decoder with exact GF(8) arithmetic."""

import argparse
import importlib.util
import json
from pathlib import Path
import time

import kickmix as km

if __package__:
    from .gf8_backend import GF8Builder, digest
else:
    from gf8_backend import GF8Builder, digest

HERE = Path(__file__).resolve().parent
BASELINE = HERE.parent / "other_tasks/dqi-baseline.json"


def load_reference(data):
    relative = "verifiers/dqi/dqi_reference_circuit.py"
    reference = Path(data) / relative
    expected = json.loads(BASELINE.read_text())["input_sha256"][relative]
    if digest(reference) != expected:
        raise ValueError("The external reference changed from the profiled baseline.")
    spec = importlib.util.spec_from_file_location("qf_external_dqi_reference", reference)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build(data, *, m=4095, n=70):
    if not (isinstance(m, int) and isinstance(n, int) and 2 <= n < m <= 4095):
        raise ValueError("Expected integer parameters 2 <= n < m <= 4095.")
    reference = load_reference(data)
    builder = GF8Builder()
    syndromes = [builder.create_quantum_register(12, name=f"S{i}") for i in range(n)]
    outcomes = [builder.create_classical_register(12, name=f"E{j}") for j in range(m)]
    reference.ProosAndZalkaDecoder(n, km.GF2Field(12)).build_syndrome_decoding_mbuc(
        builder, syndromes=syndromes, error_mbuc=outcomes,
        n_iters=reference.default_n_iters(n))
    return builder.finish_circuit(), dict(builder.counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "local/dqi-gf8.kmb")
    parser.add_argument("--report", type=Path, default=HERE / "build.json")
    parser.add_argument("--m", type=int, default=4095)
    parser.add_argument("--n", type=int, default=70)
    args = parser.parse_args()
    started = time.monotonic()
    circuit, counts = build(args.data, m=args.m, n=args.n)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    circuit.to_file(args.output, format="kmb")
    stats = circuit.operation_statistics()
    record = {
        "status": "BUILT_NOT_VERIFIED", "m": args.m, "n": args.n, "degree": 12,
        "num_qubits": circuit.num_qubits, "static_toffoli_count": stats.num_toffoli_like,
        "num_ops": len(circuit), "num_bits": circuit.num_bits,
        "operation_counts": dict(stats), "arithmetic_calls": counts,
        "circuit_sha256": digest(args.output), "circuit_bytes": args.output.stat().st_size,
        "source_sha256": {p.name: digest(p) for p in
                          (Path(__file__), HERE / "gf8_backend.py", HERE / "certificates.json")},
        "reference_sha256": json.loads(BASELINE.read_text())["input_sha256"]["verifiers/dqi/dqi_reference_circuit.py"],
        "elapsed_seconds": time.monotonic() - started,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2), flush=True)


if __name__ == "__main__":
    main()
