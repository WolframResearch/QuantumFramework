"""Native cleanup, conditional phases and decoder integration checks."""

import importlib
import os
from pathlib import Path
import sys
import unittest

import kickmix as km
import numpy as np

from .gf8_backend import GF8Builder, GF, linear_gates
from .generator import build
from other_tasks.field_probes import packed, unpack, check_circuit


class IntegrationTests(unittest.TestCase):
    def test_dirty_target_and_quantum_control(self):
        for controlled in (False, True):
            builder = GF8Builder()
            a, b, target = (builder.create_quantum_register(12, name=name)
                            for name in ("a", "b", "target"))
            control = builder.create_quantum_register(1, name="control")[0] if controlled else True
            builder.ixor_gf2_mul(a, b, target=target, control=control)
            circuit = builder.finish_circuit()
            self.assertEqual(circuit.operation_statistics().num_toffoli_like, 54 if controlled else 42)
            self.assertEqual(check_circuit(circuit)["shots"], 8192)

    def test_native_inverse_cleanup_with_all_masks(self):
        builder = GF8Builder()
        source = builder.create_quantum_register(12, name="source")
        mask = builder.create_classical_register(12, name="mask")
        one = builder.alloc_qubits(12)
        builder.x(one[0])
        output, chain = builder.init_gf2_inverse_with_scaffold(source)
        builder.gf2_phase_by_product(output, one, mask=mask, field=GF)
        builder.del_gf2_inverse_with_scaffold(source, target=output, scaffold=chain, field=GF)
        builder.x(one[0])
        builder.free(one + output + chain)
        circuit = builder.finish_circuit()
        self.assertEqual(len(chain), 48)
        self.assertEqual(circuit.operation_statistics().num_toffoli_like, 210)
        values = np.arange(4096, dtype=np.uint64)
        masks = np.random.default_rng(58).integers(0, 4096, len(values), dtype=np.uint64)
        sim = km.Simulator(batch_size=len(values), seed=613)
        sim.write_across_shots(source, packed(values))
        sim.write_across_shots(mask, packed(masks))
        sim.do(circuit)
        self.assertTrue(np.array_equal(unpack(sim.read_across_shots(source)), values))
        self.assertTrue(np.array_equal(unpack(sim.read_across_shots(mask)), masks))
        expected = [bool((GF.invert(int(a)) & int(m)).bit_count() % 2) for a, m in zip(values, masks)]
        self.assertTrue(np.array_equal(sim.read_phase_flipped_across_shots(), expected))
        scratch = km.array([km.q(i) for i in range(circuit.num_qubits) if km.q(i) not in set(source)])
        self.assertFalse(sim.read_across_shots(scratch).any())

    def test_conditional_division_with_native_cleanup(self):
        builder = GF8Builder()
        a, b, answer = (builder.create_quantum_register(12, name=name)
                        for name in ("a", "b", "answer"))
        mask = builder.create_classical_register(12, name="mask")
        enabled = builder.create_classical_register(1, name="enabled")
        builder.push_condition(enabled[0])
        quotient, chain = builder.init_gf2_div_with_scaffold(a, b)
        builder.cx(quotient, answer)
        builder.gf2_phase_by_product(quotient, a, mask=mask, field=GF)
        builder.del_gf2_div_with_scaffold(a, b, target=quotient, scaffold=chain, field=GF)
        builder.free(quotient + chain)
        builder.pop_condition()
        circuit = builder.finish_circuit()
        self.assertEqual(len(chain), 60)
        self.assertEqual(circuit.operation_statistics().num_toffoli_like, 252)
        rng = np.random.default_rng(18)
        divisors = np.tile(np.arange(4096, dtype=np.uint64), 2)
        numerators = rng.integers(0, 4096, len(divisors), dtype=np.uint64)
        masks = rng.integers(0, 4096, len(divisors), dtype=np.uint64)
        active = np.repeat([0, 1], 4096).astype(np.uint64)
        sim = km.Simulator(batch_size=len(divisors), seed=513)
        for register, values in ((a, numerators), (b, divisors), (mask, masks)):
            sim.write_across_shots(register, packed(values))
        sim.write_across_shots(enabled, packed(active, 1))
        sim.do(circuit)
        expected = np.array([GF.mul(int(x), GF.invert(int(y))) for x, y in zip(numerators, divisors)], dtype=np.uint64) * active
        for register, values in ((a, numerators), (b, divisors), (answer, expected), (mask, masks), (enabled, active)):
            self.assertTrue(np.array_equal(unpack(sim.read_across_shots(register)), values))
        phases = [bool((GF.mul(int(q), int(x)) & int(m)).bit_count() % 2)
                  for q, x, m in zip(expected, numerators, masks)]
        self.assertTrue(np.array_equal(sim.read_phase_flipped_across_shots(), phases))
        public = set(a + b + answer)
        scratch = km.array([km.q(i) for i in range(circuit.num_qubits) if km.q(i) not in public])
        self.assertFalse(sim.read_across_shots(scratch).any())

    def test_frobenius_maps(self):
        builder = GF8Builder()
        for power in builder.frobenius_gates:
            with self.subTest(power=power):
                values = np.arange(4096, dtype=np.uint64)
                probe = GF8Builder()
                target = probe.create_quantum_register(12, name="target")
                probe._frobenius(power, target)
                sim = km.Simulator(batch_size=len(values), seed=18)
                sim.write_across_shots(target, packed(values))
                sim.do(probe.finish_circuit())
                expected = np.array([GF.frobenius(int(x), power) for x in values])
                self.assertTrue(np.array_equal(unpack(sim.read_across_shots(target)), expected))

    def test_false_control_and_invalid_inputs(self):
        builder = GF8Builder()
        a, b, target = (builder.create_quantum_register(12, name=name) for name in ("a", "b", "target"))
        builder.ixor_gf2_mul(a, b, target=target, control=False)
        self.assertEqual(len(builder.finish_circuit()), 0)
        with self.assertRaises(ValueError):
            builder.ixor_gf2_mul(a, b, target=target, field=km.GF2Field(3))
        with self.assertRaises(ValueError):
            builder.init_gf2_inverse_with_scaffold(a[:11])
        with self.assertRaises(ValueError):
            builder.init_gf2_div_with_scaffold(a, b, target="wrong")
        with self.assertRaises(ValueError):
            linear_gates(np.zeros((12, 12), dtype=np.uint8))
        with self.assertRaises(ValueError):
            build(Path("unused"), m=2, n=4)

    def test_reduced_decoder_official_protocol(self):
        data = Path(os.environ["QF_INVSQRT_DATA"])
        circuit, calls = build(data, m=15, n=4)
        sys.path.insert(0, str(data / "verifiers"))
        verifier = importlib.import_module("dqi.verify_submission_dqi")
        spec = verifier.ProblemSpecDQI.benchmark_instance(m=15, n=4, num_shots=2048)
        report = verifier.verify(circuit, spec, salt="qf-gf8-integration-small-v1", max_workers=1)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["value_errors"], 0)
        self.assertEqual(report["phase_errors"], 0)
        self.assertGreater(calls["division_scaffolds"], 0)
        self.assertEqual(calls["inverse_scaffolds"], 15)


if __name__ == "__main__":
    unittest.main()
