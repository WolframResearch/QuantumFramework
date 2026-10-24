"""Native-basis arithmetic adapter for the external DQI reference decoder.

Only the three arithmetic entry points used by that decoder are replaced.
Other builder calls, including native zero-Toffoli cleanup, are delegated.
"""

from collections import Counter
import hashlib
import json
from pathlib import Path

import kickmix as km
import numpy as np

HERE = Path(__file__).resolve().parent
FEASIBILITY = HERE.parent / "other_tasks"
GF = km.GF2Field(12)
GF8 = km.GF2Field(3)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def linear_gates(matrix):
    """Synthesize an invertible binary map by reversed row elimination."""
    reduced = np.array(matrix, dtype=np.uint8, copy=True)
    if reduced.shape != (12, 12) or not np.isin(reduced, [0, 1]).all():
        raise ValueError("Expected a 12 by 12 binary matrix.")
    gates = []
    for column in range(12):
        pivot = next((row for row in range(column, 12) if reduced[row, column]), None)
        if pivot is None:
            raise ValueError("The binary map must be invertible.")
        if pivot != column:
            reduced[[column, pivot]] = reduced[[pivot, column]]
            gates.append(("swap", column, pivot))
        for row in range(12):
            if row != column and reduced[row, column]:
                reduced[row] ^= reduced[column]
                gates.append(("cx", column, row))
    return tuple(reversed(gates))


def xor_edges(matrix):
    return tuple((source, target) for target, row in enumerate(matrix)
                 for source, value in enumerate(row) if value)


class GF8Builder:
    """Adapter for degree 12 and native modulus 4179; not a general builder API."""

    def __init__(self):
        symbolic = json.loads((FEASIBILITY / "symbolic.json").read_text())
        certificate = json.loads((HERE / "certificates.json").read_text())
        if (certificate["status"] != "PASS" or
                not all(certificate["checks"].values()) or
                certificate["source_sha256"] != digest(HERE / "run-certificates.wls") or
                certificate["field_source_sha256"] != digest(FEASIBILITY / "FieldArithmetic.wl") or
                certificate["symbolic_sha256"] != digest(FEASIBILITY / "symbolic.json")):
            raise ValueError("Regenerate the exact WL certificates before compiling.")
        plan = next(p for p in symbolic["tower_plans"] if p["subfield_degree"] == 3)
        if not all(plan[k] for k in ("basis_certificate", "interpolation_certificate", "tensor_certificate")):
            raise ValueError("The native-basis multiplication certificate failed.")
        if GF.modulus != 4179 or not symbolic["six_product"]["certificate"]:
            raise ValueError("The native field or GF(8) identity does not match.")
        self.raw = km.CircuitBuilder()
        self.counts = Counter()
        self.evaluations = tuple(xor_edges(m) for m in plan["evaluation_maps"])
        self.outputs = tuple(xor_edges(m) for m in plan["output_maps"])
        six = symbolic["six_product"]
        self.forms = tuple(tuple(i for i, value in enumerate(row) if value)
                           for row in six["input_forms"])
        self.product_outputs = tuple(tuple(i for i, row in enumerate(six["output_forms"]) if row[term])
                                     for term in range(6))
        self.frobenius_gates = {int(p): linear_gates(m) for p, m in certificate["frobenius_maps"].items()}

    def __getattr__(self, name):
        return getattr(self.raw, name)

    @staticmethod
    def _field(field, *registers):
        field = GF if field is None else field
        if field.degree != 12 or field.modulus != 4179 or any(len(r) != 12 for r in registers):
            raise ValueError("This adapter requires native GF(4096) registers of width 12.")

    def _xor(self, edges, source, target):
        if edges:
            self.raw.cx(km.array([source[a] for a, _ in edges]),
                        km.array([target[b] for _, b in edges]))

    def _frobenius(self, power, target):
        for gate, a, b in self.frobenius_gates[power % 12]:
            getattr(self.raw, gate)(target[a], target[b])

    def _multiply(self, lhs, rhs, target):
        self.counts["gf8_multiplications"] += 1
        builder = self.raw
        with builder.mark("gf8_mul"):
            left, right, product = (builder.alloc_qubits(3) for _ in range(3))
            temp = builder.alloc_qubits(1)[0]
            builder.reset(left + right)
            for evaluate, output in zip(self.evaluations, self.outputs):
                self._xor(evaluate, lhs, left)
                self._xor(evaluate, rhs, right)
                builder.reset(product)
                for form, targets in zip(self.forms, self.product_outputs):
                    first, *rest = form
                    for bit in rest:
                        builder.cx(left[bit], left[first])
                        builder.cx(right[bit], right[first])
                    builder.init_and(left[first], right[first], temp)
                    for bit in targets:
                        builder.cx(temp, product[bit])
                    builder.del_and(left[first], right[first], temp)
                    for bit in reversed(rest):
                        builder.cx(right[bit], right[first])
                        builder.cx(left[bit], left[first])
                self._xor(output, product, target)
                builder.del_gf2_mul(left, right, target=product, field=GF8)
                self._xor(evaluate, rhs, right)
                self._xor(evaluate, lhs, left)
            builder.free(km.array([temp]) + product + right + left)

    def ixor_gf2_mul(self, lhs, rhs, *, target, field=None, control=True):
        self._field(field, lhs, rhs, target)
        if control is False:
            return
        if control is True:
            self._multiply(lhs, rhs, target)
        else:
            self.counts["controlled_multiplications"] += 1
            product = self.raw.alloc_qubits(12)
            self.raw.reset(product)
            self._multiply(lhs, rhs, product)
            self.raw.ccx(control, product, target)
            self.raw.del_gf2_mul(lhs, rhs, target=product, field=GF)
            self.raw.free(product)

    def _inverse(self, source, target, scaffold):
        """Same ordered beta(2), beta(4), beta(8), beta(10) scaffold as native."""
        self.raw.reset(target + scaffold)
        nodes = {1: source, **{k: scaffold[i*12:(i+1)*12]
                              for i, k in enumerate((2, 4, 8, 10))}}
        for previous, current in ((1, 2), (2, 4), (4, 8)):
            self.raw.cx(nodes[previous], target)
            self._frobenius(previous, target)
            self._multiply(nodes[previous], target, nodes[current])
            self._frobenius(-previous, target)
            self.raw.cx(nodes[previous], target)
        self._frobenius(2, nodes[8])
        self._multiply(nodes[8], nodes[2], nodes[10])
        self._frobenius(1, nodes[10])
        self._multiply(nodes[10], source, target)
        self._frobenius(1, target)

    def init_gf2_inverse_with_scaffold(self, source, *, target="alloc", field=None):
        self._field(field, source)
        if isinstance(target, str):
            if target != "alloc":
                raise ValueError("Expected target='alloc' or a native register.")
            target = self.raw.alloc_qubits(12)
        self._field(field, target)
        self.counts["inverse_scaffolds"] += 1
        scaffold = self.raw.alloc_qubits(48)
        self._inverse(source, target, scaffold)
        return target, scaffold

    def init_gf2_div_with_scaffold(self, lhs, rhs, *, target="alloc", field=None):
        self._field(field, lhs, rhs)
        if isinstance(target, str):
            if target != "alloc":
                raise ValueError("Expected target='alloc' or a native register.")
            target = self.raw.alloc_qubits(12)
        self._field(field, target)
        self.counts["division_scaffolds"] += 1
        self.raw.reset(target)
        scaffold = self.raw.alloc_qubits(60)
        inverse, chain = scaffold[:12], scaffold[12:]
        self._inverse(rhs, inverse, chain)
        self._multiply(lhs, inverse, target)
        return target, scaffold
