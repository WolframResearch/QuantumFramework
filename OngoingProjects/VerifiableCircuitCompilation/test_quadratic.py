"""Compiler-equivalence and frozen-coefficient checks before full verification."""

import dataclasses
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import kickmix as km

from feasibility import resource_probe as probe
from quadratic import (DESIGN, DESIGN_SHA256, HERE, TABLE_SHA256,
                       build_circuit, build_to_file, canonical_bytes, load_coefficients)


class QuadraticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = Path(os.environ["QF_INVSQRT_DATA"])
        cls.nr, cls.verifier = probe.load_official(cls.data)
        cls.design = json.loads(DESIGN.read_text())
        cls.circuit = build_circuit()

    def test_exact_frozen_coefficients(self):
        coefficients = load_coefficients()
        self.assertEqual(len(coefficients), 192)
        self.assertEqual(coefficients, tuple(tuple(row) for row in self.design["table"]))
        self.assertEqual(hashlib.sha256(DESIGN.read_bytes()).hexdigest(), DESIGN_SHA256)
        self.assertEqual(hashlib.sha256(json.dumps(coefficients, separators=(",", ":")).encode()).hexdigest(), TABLE_SHA256)
        with tempfile.TemporaryDirectory() as directory:
            corrupted = json.loads(DESIGN.read_text())
            corrupted["table"][0][0] += 1
            path = Path(directory) / "changed.json"
            path.write_text(json.dumps(corrupted))
            with self.assertRaisesRegex(ValueError, "fixed certified design"):
                load_coefficients(path)

    def test_byte_identical_to_reviewed_feasibility_circuit(self):
        reference, _ = probe.build_poly_probe(self.design, self.nr)
        self.assertEqual(canonical_bytes(self.circuit), canonical_bytes(reference))

    def test_expected_resource_counts_and_registers(self):
        statistics = self.circuit.operation_statistics()
        self.assertEqual(self.circuit.num_qubits, 220)
        self.assertEqual(statistics.num_toffoli_like, 4449)
        self.assertEqual(statistics.num_t_like, 1)
        self.assertEqual(statistics.num_z_pow(min_clifford_hierarchy_level=4), 43)
        self.assertEqual(list(self.circuit.register_data), ["x"])
        self.assertEqual(len(self.circuit.register_data["x"]), 22)

    def test_boundaries_match_independent_integer_model(self):
        report = probe.simulator_model_samples(self.circuit, self.design, self.verifier)
        self.assertTrue(report["phase_matches_integer_model_exactly"])
        self.assertTrue(report["scratch_clean"])
        self.assertTrue(report["x_restored"])
        self.assertGreaterEqual(report["unique_inputs"], 93)

    def test_official_verifier_small_domain(self):
        spec = dataclasses.replace(self.verifier.ProblemSpecInvSqrt.for_n_p(10), max_x=1024)
        report = self.verifier.verify(self.circuit, spec, salt="quadratic-small-domain-test", max_workers=1)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["verified_shots"], 65_536)

    def test_serialized_artifact_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "quadratic.kmb"
            report = build_to_file(path)
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), report["circuit_sha256"])
            self.assertEqual(canonical_bytes(km.Circuit.from_file(path)), canonical_bytes(self.circuit))

    def test_build_without_competition_or_research_modules(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "feasibility").mkdir()
            shutil.copy2(HERE / "quadratic.py", root / "quadratic.py")
            shutil.copy2(DESIGN, root / "feasibility/B-quadratic-192.json")
            output = root / "standalone.kmb"
            subprocess.run([sys.executable, "-I", str(root / "quadratic.py"), "--output", str(output)],
                           check=True, capture_output=True, text=True)
            self.assertEqual(output.read_bytes(), canonical_bytes(self.circuit))


if __name__ == "__main__":
    unittest.main()
