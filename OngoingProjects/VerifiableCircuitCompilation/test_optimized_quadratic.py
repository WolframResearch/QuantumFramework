"""Rounded certificate correspondence, lifetime scheduling and official checks."""

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
import numpy as np

import optimized_quadratic as opt
from optimization.analyze import model, sample
from quadratic import HERE, DESIGN, canonical_bytes
from verify_quadratic import load_verifier
from verify_optimized_quadratic import check_evidence


class ReducedQuadraticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.verifier, _ = load_verifier(Path(os.environ["QF_INVSQRT_DATA"]))
        cls.circuit = opt.build_circuit()

    def test_exact_correspondence_to_symbolic_certificate(self):
        certs = json.loads((HERE/"optimization/width-certificates.json").read_text())["certificates"]
        cert = next(c for c in certs if c["fractional_bits"] == list(opt.FRACTIONAL_BITS))
        self.assertTrue(cert["certified"])
        self.assertEqual(cert["table"], [list(row) for row in opt.quantized_rows()])
        self.assertEqual(cert["register_bits"], [30,22,15])
        self.assertEqual(len(check_evidence()), 12)

    def test_compiled_resources(self):
        stats = self.circuit.operation_statistics()
        self.assertEqual((self.circuit.num_qubits, stats.num_toffoli_like, stats.num_t_like,
                          stats.num_z_pow(min_clifford_hierarchy_level=4)), (142,3460,1,37))
        self.assertEqual(list(self.circuit.register_data), ["x"])
        self.assertEqual(len(self.circuit.register_data["x"]), 22)

    def test_all_boundary_samples_and_measurement_cleanup(self):
        report = sample(self.circuit, opt.FRACTIONAL_BITS, self.verifier)
        self.assertGreaterEqual(report["unique_inputs"], 4095)
        self.assertTrue(report["scratch_clean"])

    def test_retained_coefficient_and_controlled_arithmetic_match_model(self):
        for arithmetic, recycle in (("padded",False),("controlled",True)):
            with self.subTest(arithmetic=arithmetic,recycle=recycle):
                circuit = opt.build_circuit(arithmetic=arithmetic,recycle_coefficient=recycle)
                sample(circuit, opt.FRACTIONAL_BITS, self.verifier)

    def test_official_default_shots_small_domain(self):
        spec = dataclasses.replace(self.verifier.ProblemSpecInvSqrt.for_n_p(10),max_x=4096)
        result = self.verifier.verify(self.circuit,spec,salt="reduced-quadratic-unit-test",max_workers=1)
        self.assertEqual(result["status"],"PASS")
        self.assertEqual(result["verified_shots"],262_144)

    def test_too_narrow_phase_precision_is_detected(self):
        xs = np.arange(1,4097,dtype=np.uint64)
        word,f,_,_ = model(xs,(24,20,14))
        phase = word.astype(np.uint64) << np.uint64(64-f)
        expected,allowed = self.verifier.reference_values(self.verifier.ProblemSpecInvSqrt.for_n_p(10),xs)
        used = np.abs((phase-expected).view(np.int64)).astype(np.uint64)
        self.assertGreater(int((used>allowed).sum()),0)

    def test_input_validation(self):
        for widths in ((0,23,16),(30.,23,16),(30,23),(37,23,16)):
            with self.assertRaises(ValueError):opt.quantized_rows(widths)
        with self.assertRaises(ValueError):model([0],opt.FRACTIONAL_BITS)
        with self.assertRaises(ValueError):opt.build_circuit(arithmetic="unknown")

    def test_binary_and_isolated_generator(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/"feasibility").mkdir()
            for name in ("quadratic.py","optimized_quadratic.py"):
                shutil.copy2(HERE/name,root/name)
            shutil.copy2(DESIGN,root/"feasibility/B-quadratic-192.json")
            artifact = root/"standalone.kmb"
            # -I removes cwd from sys.path, so explicitly permit only this bundle.
            subprocess.run([sys.executable,"-I","-c",
                "import sys,runpy;sys.path.insert(0,sys.argv[1]);sys.argv=sys.argv[2:];runpy.run_path(sys.argv[0],run_name='__main__')",
                str(root),str(root/"optimized_quadratic.py"),"--output",str(artifact)],
                check=True,capture_output=True,text=True)
            self.assertEqual(artifact.read_bytes(),canonical_bytes(self.circuit))
            self.assertEqual(canonical_bytes(km.Circuit.from_file(artifact)),artifact.read_bytes())

    def test_incumbent_source_unchanged(self):
        incumbent = json.loads((HERE/"results/quadratic-verified.json").read_text())
        self.assertEqual(hashlib.sha256((HERE/"quadratic.py").read_bytes()).hexdigest(),
                         incumbent["build"]["generator_sha256"])


if __name__ == "__main__":
    unittest.main()
