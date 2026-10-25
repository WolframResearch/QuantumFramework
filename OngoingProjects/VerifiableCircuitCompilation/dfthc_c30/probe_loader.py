"""Check arbitrary-target XOR, cleanup and involution beyond official chaining."""

import argparse
import json
from pathlib import Path

import kickmix as km
import numpy as np

from generator import HERE, digest
from loader import batch_deltas, build, mobius_transform


def check_family(circuits, deltas, repeats, seed):
    rows = len(deltas[0])
    addresses = np.tile(np.arange(rows), repeats)
    selection = circuits[0].register_data["selection"]
    word_regs = [circuits[0].register_data[f"word_{j}"] for j in range(3)]
    address_bits = ((addresses[None, :] >> np.arange(len(selection))[:, None]) & 1).astype(bool)
    rng = np.random.default_rng(seed)
    sim = km.Simulator(batch_size=len(addresses), seed=seed, counting_operations=True)
    sim.use_same_registers_as(circuits[0])
    checks = 0
    for circuit, delta in zip(circuits, deltas, strict=True):
        dirty = rng.integers(0, 1024, size=(len(addresses), 3), dtype=np.int64)
        sim.clear_for_shot()
        sim.write_across_shots(selection, address_bits)
        for j, register in enumerate(word_regs):
            sim.write_across_shots(register, ((dirty[:, j][None, :] >> np.arange(10)[:, None]) & 1).astype(bool))
        public = {q.id for q in selection} | {q.id for register in word_regs for q in register}
        work = km.array([km.q(j) for j in range(circuit.num_qubits) if j not in public])
        for expected in (dirty ^ delta[addresses], dirty):
            sim.do(circuit)
            actual = np.stack([(np.int64(1) << np.arange(10)) @ sim.read_across_shots(register).astype(np.int64)
                               for register in word_regs], axis=1)
            if not np.array_equal(actual, expected):
                raise AssertionError("Dirty-target XOR or involution failed.")
            if (np.any(sim.read_phase_flipped_across_shots()) or
                    not np.array_equal(sim.read_across_shots(selection), address_bits) or
                    np.any(sim.read_across_shots(work))):
                raise AssertionError("Phase, selection or workspace invariant failed.")
            checks += len(addresses)
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--submission", type=Path, default=HERE / "local/standalone")
    args = parser.parse_args()
    with np.load(args.submission / "submission.npz", allow_pickle=False) as z:
        words = z["angle_words"]
    circuits = [km.Circuit.from_file(args.submission / f"angle_k_{k}.kmb", format="kmb") for k in range(25)]
    full_checks = check_family(circuits, batch_deltas(words, 3), 4, 897235)
    synthetic_checks = 0
    rng = np.random.default_rng(987521)
    for rows in (16, 20, 32, 128):
        for kind in ("zero", "constant", "nonconstant"):
            table = np.zeros((rows, 3), dtype=np.int64)
            if kind == "constant":
                table[:] = [255, 768, 1023]
            elif kind == "nonconstant":
                table = rng.integers(-256, 256, size=(rows, 3)) % 1024
            synthetic_checks += check_family([build(table)], [table], 4, rows)
    invalid = [np.zeros((0, 3), dtype=int), np.zeros((16, 2), dtype=int),
               np.zeros((16, 3)), np.full((16, 3), -1), np.full((16, 3), 256),
               np.full((16, 3), 1024), np.zeros((18, 3), dtype=int)]
    for table in invalid:
        try:
            build(table)
        except ValueError:
            pass
        else:
            raise AssertionError("Malformed input accepted.")
    for table in ([1, 2, 3], [], [1] * 6):
        try:
            mobius_transform(table)
        except ValueError:
            pass
        else:
            raise AssertionError("Malformed ANF table accepted.")
    for j in (0, -1, True, 1.5):
        try:
            batch_deltas(words, j)
        except ValueError:
            pass
        else:
            raise AssertionError("Malformed batch size accepted.")
    report = {"status": "PASS", "seed": 897235, "released_dirty_target_shots_per_address": 4,
              "released_circuit_applications_checked": full_checks,
              "synthetic_circuit_applications_checked": synthetic_checks,
              "synthetic_rows": [16, 20, 32, 128], "synthetic_regimes": ["zero", "constant", "nonconstant"],
              "invalid_cases_rejected": len(invalid) + 3 + 4,
              "phase_errors": 0, "workspace_errors": 0, "selection_errors": 0, "value_errors": 0,
              "source_sha256": {p.name: digest(p) for p in (HERE / "loader.py", Path(__file__))}}
    (HERE / "probes.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
