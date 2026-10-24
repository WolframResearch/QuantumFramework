"""Independent reference checks for arithmetic and component-mixing probes."""

import hashlib
import json
import math
from pathlib import Path
import unittest

import numpy as np

from other_tasks import component_mixing as mixing
from other_tasks import field_probes as fields

HERE=Path(__file__).resolve().parent


class OtherTasksTests(unittest.TestCase):
    def test_piecewise_pair_search_against_dense_reference(self):
        rng=np.random.default_rng(5813)
        cases=[(np.zeros(7),np.zeros(7)),(np.ones(7),np.ones(7)),
               (np.array([1.,1.]),np.array([1.,-1.]))]
        cases.extend((rng.normal(size=7),rng.normal(size=7)) for _ in range(12))
        theta=np.linspace(0,math.pi/2,4097)
        for a,b in cases:
            angle,value=mixing.best_pair(a,b)
            c,s=np.cos(theta)[:,None],np.sin(theta)[:,None]
            grid=np.abs(c*a+s*b).sum(axis=1)**2+np.abs(-s*a+c*b).sum(axis=1)**2
            self.assertLessEqual(value,float(grid.min())+1e-10)
            self.assertAlmostEqual(value,mixing.pair_objective(a,b,angle),places=10)

    def test_mixing_preserves_gram_and_improves_or_retains_objective(self):
        rng=np.random.default_rng(3141);a=rng.normal(size=(3,5,6))
        b,metadata=mixing.mix_pairs(a,rounds=2,pairs_per_group=6)
        self.assertEqual(metadata["attempted_pair_rotations"],36)
        self.assertTrue(np.allclose(a@a.transpose(0,2,1),b@b.transpose(0,2,1),atol=1e-12,rtol=1e-12))
        self.assertLessEqual(float(np.square(np.abs(b).sum(axis=1)).sum()),
                             float(np.square(np.abs(a).sum(axis=1)).sum())+1e-12)

    def test_spectral_control_preserves_gram(self):
        a=np.random.default_rng(871).normal(size=(2,8,5));b=mixing.spectral_mix(a)
        self.assertTrue(np.allclose(a@a.transpose(0,2,1),b@b.transpose(0,2,1),atol=1e-12,rtol=1e-12))

    def test_invalid_pair_inputs(self):
        for a,b in (([],[]),([1],[1,2]),([float("nan")],[1])):
            with self.assertRaises(ValueError):mixing.best_pair(a,b)

    def test_both_official_baselines_pass_full_protocols(self):
        d=json.loads((HERE/"dqi-baseline.json").read_text())
        f=json.loads((HERE/"dfthcsos-baseline.json").read_text())
        self.assertEqual((d["status"],d["verified_shots"],d["log2_delta"]),("PASS",102400,-128.))
        self.assertEqual(d["value_errors"]+d["phase_errors"],0)
        self.assertGreater(d["repetition_factor"],1)
        self.assertEqual((f["status"],f["verified_shots"],f["num_circuits"]),("PASS",254464,10))
        self.assertEqual(f["resources"]["angles_toffoli"],79480)

    def test_arithmetic_probes_reproduce_resources_and_cleanup(self):
        symbolic=json.loads((HERE/"symbolic.json").read_text())
        plan=next(p for p in symbolic["tower_plans"] if p["subfield_degree"]==3)
        circuit=fields.multiplication_circuit(plan,True)
        self.assertEqual((circuit.num_qubits,circuit.operation_statistics().num_toffoli_like),(59,54))
        self.assertTrue(fields.check_circuit(circuit)["scratch_clean"])
        circuit=fields.inverse_phase_circuit(plan)
        self.assertEqual((circuit.num_qubits,circuit.operation_statistics().num_toffoli_like),(82,210))
        self.assertTrue(fields.check_circuit(circuit,True)["phase_correct"])

    def test_report_source_bindings(self):
        sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
        symbolic=json.loads((HERE/"symbolic.json").read_text())
        for name,digest in symbolic["source_sha256"].items():self.assertEqual(sha(HERE/name),digest)
        probes=json.loads((HERE/"field-results.json").read_text())
        self.assertEqual(probes["source_sha256"],sha(HERE/"field_probes.py"))
        self.assertEqual(probes["symbolic_sha256"],sha(HERE/"symbolic.json"))
        mixed=json.loads((HERE/"mixing-results.json").read_text())
        self.assertEqual(mixed["analysis_source_sha256"],sha(HERE/"component_mixing.py"))
        self.assertEqual(mixed["symbolic_sha256"],sha(HERE/"symbolic.json"))
        self.assertEqual(mixed["full_official_best_result"]["status"],"PASS")
        self.assertGreater(mixed["best_improvement"],1)

    def test_inverse_phase_zero_and_multibit_masks(self):
        symbolic=json.loads((HERE/"symbolic.json").read_text())
        plan=next(p for p in symbolic["tower_plans"] if p["subfield_degree"]==3)
        for mask in (0,0xA55):
            circuit=fields.inverse_phase_circuit(plan,mask)
            self.assertTrue(fields.check_circuit(circuit,True,mask)["phase_correct"])


if __name__=="__main__":unittest.main()
