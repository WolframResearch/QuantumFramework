"""Fast contract and model/circuit cross-checks; full verification is experiment.py."""
import dataclasses
import json
import os
from pathlib import Path
import unittest

import kickmix as km
import numpy as np

from experiment import HERE, load_competition


class PrototypeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = Path(os.environ["QF_INVSQRT_DATA"])
        cls.nr, cls.verifier = load_competition(cls.data)
        cls.spec = dataclasses.replace(cls.verifier.ProblemSpecInvSqrt.for_n_p(10), max_x=4)

    def test_missing_phase_is_rejected(self):
        b = km.CircuitBuilder()
        b.create_quantum_register(22, name="x")
        with self.assertRaisesRegex(AssertionError, "phase error"):
            self.verifier.verify(b.finish_circuit(), self.spec, max_workers=1)

    def test_changed_input_is_rejected(self):
        b = km.CircuitBuilder()
        x = b.create_quantum_register(22, name="x")
        b.x(x[0])
        with self.assertRaisesRegex(AssertionError, "x was not restored"):
            self.verifier.verify(b.finish_circuit(), self.spec, max_workers=1)

    def test_dirty_scratch_is_rejected(self):
        b = km.CircuitBuilder()
        b.create_quantum_register(22, name="x")
        b.x(b.alloc_qubits(1)[0])
        with self.assertRaisesRegex(AssertionError, "Scratch qubit"):
            self.verifier.verify(b.finish_circuit(), self.spec, max_workers=1)

    def test_candidate_compute_matches_integer_model(self):
        widths = json.loads((HERE / "results/candidate-model.json").read_text())["precision"]
        precision = self.nr.Precision(**widths)
        circuit = self.nr.build_nr_circuit(10, prec=precision, compute_only=True)
        rng = np.random.default_rng(20261022)
        xs = np.r_[np.arange(1, 33), np.array([1 << k for k in range(22)]),
                   rng.integers(1, 3_133_453, size=202)].astype(np.uint64)
        sim = km.Simulator(batch_size=len(xs), seed=71)
        sim.write_across_shots(circuit.register_data["x"], self.verifier.to_bits(xs, 22))
        sim.do(circuit)
        observed = self.verifier.from_bits(sim.read_across_shots(circuit.register_data["y1"]))
        expected = self.nr.fixed_point_model(10, precision, xs)
        np.testing.assert_array_equal(observed, expected)

    def test_candidate_resources_improve(self):
        baseline = km.Circuit.from_file(self.data / "sample_submission/invsqrt.kmb")
        candidate = km.Circuit.from_file(HERE / "results/candidate.kmb")
        self.assertLess(candidate.num_qubits, baseline.num_qubits)
        self.assertLess(candidate.operation_statistics().num_toffoli_like,
                        baseline.operation_statistics().num_toffoli_like)


if __name__ == "__main__":
    unittest.main()
