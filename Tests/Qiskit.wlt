BeginTestSection["QuantumQASM"]

(* This suite needs a working qiskit Python session. It SELF-GUARDS: when qiskit is
   unavailable (e.g. CI, which runs the build-paclet action, not RunTests.wls, and has
   no Python) every Python-touching test trivially passes. The `qtest` helper holds its
   arguments and only evaluates them when qiskitOK is True. *)

qiskitOK = TrueQ @ Quiet @ Check[Head[QuantumCircuitOperator["Bell"]["Qiskit"]] === QiskitCircuit, False];

SetAttributes[qtest, HoldAll];
qtest[actual_, expected_, id_] := VerificationTest[If[qiskitOK, actual, "skip"], If[qiskitOK, expected, "skip"], TestID -> id];

pe = Wolfram`QuantumFramework`PackageScope`PythonEvaluate;

qc    = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}];   (* Bell + measure *)
qcU   = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}];             (* unitary, no measure *)
cmap  = {{1, 0}, {1, 3}, {2, 0}, {2, 3}};
tspec = <|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> cmap|>;

(* Ground truth: raw qiskit transpile + dumps in the SAME session, for string cross-check. *)
gtQASM[bytes_, optsAssoc_, version_] := Block[{pb = bytes, po = optsAssoc, pv = version},
    pe[Context[pb], "
import pickle
from qiskit import transpile, qasm2, qasm3
from qiskit.transpiler import CouplingMap
c = pickle.loads(<* pb *>)
kw = dict(<* po *>)
if kw.get('coupling_map') is not None:
    kw['coupling_map'] = CouplingMap([tuple(e) for e in kw['coupling_map']])
out = transpile(c, **kw) if kw else c
(qasm3.dumps(out) if <* pv *> == 3 else qasm2.dumps(out))
"]
];

(* Correctness: transpiling to a universal basis preserves the unitary (up to global phase). *)
unitEquivQ[bytes_, basis_] := Block[{pb = bytes, pbasis = basis},
    pe[Context[pb], "
import pickle
from qiskit import transpile
from qiskit.quantum_info import Operator
c = pickle.loads(<* pb *>)
out = transpile(c, basis_gates=<* pbasis *>, optimization_level=1)
bool(Operator(out).equiv(Operator(c)))
"]
];


(* ---- Layer 1: cross-check QuantumQASM output against raw qiskit, string-for-string ---- *)

qtest[QuantumQASM[qc], gtQASM[qc["Qiskit"]["Bytes"], <||>, 3], "faithful-qasm3-matches-qiskit"]

qtest[QuantumQASM[qc, "Version" -> 2], gtQASM[qc["Qiskit"]["Bytes"], <||>, 2], "faithful-qasm2-matches-qiskit"]

qtest[
    QuantumQASM[qc, {"x", "sx", "rz", "cz"}, "OptimizationLevel" -> 1, "SeedTranspiler" -> 1234],
    gtQASM[qc["Qiskit"]["Bytes"], <|"basis_gates" -> {"x", "sx", "rz", "cz"}, "optimization_level" -> 1, "seed_transpiler" -> 1234|>, 3],
    "basis-only-matches-qiskit"
]

qtest[
    QuantumQASM[qc, {"x", "sx", "rz", "cz"}, "CouplingMap" -> cmap, "OptimizationLevel" -> 3, "SeedTranspiler" -> 1234],
    gtQASM[qc["Qiskit"]["Bytes"], <|"basis_gates" -> {"x", "sx", "rz", "cz"}, "coupling_map" -> cmap, "optimization_level" -> 3, "seed_transpiler" -> 1234|>, 3],
    "native-routed-matches-qiskit"
]


(* ---- Layer 2: semantic correctness ---- *)

qtest[unitEquivQ[qcU["Qiskit"]["Bytes"], {"x", "sx", "rz", "cz"}], True, "transpile-preserves-unitary"]

qtest[
    With[{s = QuantumQASM[qc, {"x", "sx", "rz", "cz"}, "CouplingMap" -> cmap, "OptimizationLevel" -> 3]},
        StringContainsQ[s, "cz $1, $0"] && ! StringContainsQ[s, "u2"] && ! StringContainsQ[s, "u3"]],
    True, "native-directed-cz-no-u2u3"
]

qtest[
    With[{s = QuantumQASM[qc]}, StringContainsQ[s, "h q[0]"] && StringContainsQ[s, "cx q[0], q[1]"]],
    True, "faithful-generic-gates"
]

qtest[StringStartsQ[QuantumQASM[qc, "Version" -> 2], "OPENQASM 2.0"], True, "qasm2-prefix"]

qtest[StringStartsQ[QuantumQASM[qc], "OPENQASM 3.0"], True, "qasm3-prefix"]


(* ---- Layer 3: design, consistency, downvalues ---- *)

qtest[
    QuantumQASM[qc, {"x", "sx", "rz", "cz"}, "CouplingMap" -> cmap, "OptimizationLevel" -> 1, "SeedTranspiler" -> 7] ===
        QuantumQASM[qc, {"x", "sx", "rz", "cz"}, "couplingMap" -> cmap, "optimization_level" -> 1, "seed_transpiler" -> 7],
    True, "casing-variants-identical"
]

qtest[StringQ @ QuantumQASM[qc, {"x", "sx", "rz", "cz"}, "seed_transpiler" -> 7], True, "snake-case-passthrough"]

qtest[StringStartsQ[QuantumQASM[qc["Qiskit"], {"x", "sx", "rz", "cz"}, "CouplingMap" -> cmap], "OPENQASM 3.0"], True, "from-qiskitcircuit"]

qtest[
    With[{f = QuantumQASM[qc, {"x", "sx", "rz", "cz"}, "NoSuchParam" -> 5]},
        FailureQ[f] && MemberQ[f["Unknown"], "no_such_param"]],
    True, "unknown-option-failure"
]

qtest[FailureQ @ QuantumQASM[qc, {"x", "cz"}], True, "non-universal-basis-failure"]

qtest[FailureQ @ QuantumQASM[qc, {"x", "sx", "rz", "cz"}, "CouplingMap" -> cmap, "coupling_map" -> cmap], True, "duplicate-key-failure"]

qtest[FailureQ @ QuantumQASM[qc, {"x", "sx", "rz", "cz"}, "BasisGates" -> {"x"}], True, "basis-given-twice-failure"]


(* ---- QiskitTarget ---- *)

qtest[Head @ QiskitTarget[tspec], QiskitTarget, "target-head"]

qtest[QiskitTarget[tspec]["NumQubits"], 4, "target-numqubits"]

qtest[QiskitTarget[tspec]["OperationNames"], {"cz", "measure", "reset", "rz", "sx", "x"}, "target-operationnames"]

qtest[QiskitTarget[tspec]["CouplingMap"], cmap, "target-couplingmap"]

qtest[Head @ ToBoxes @ QiskitTarget[tspec], InterpretationBox, "target-makeboxes"]

qtest[FailureQ @ QiskitTarget[<|"BasisGates" -> {"x"}, "Bogus" -> 1|>], True, "target-bad-key-failure"]

(* a MEASURED circuit transpiled via a QiskitTarget built from a bare unitary basis must
   succeed — proves the automatic measure/reset add. *)
qtest[StringStartsQ[QuantumQASM[qc, "Target" -> QiskitTarget[tspec]], "OPENQASM 3.0"], True, "measured-via-target"]

qtest[StringStartsQ[QuantumQASM[qc, "Target" -> tspec], "OPENQASM 3.0"], True, "inline-target-spec"]

EndTestSection[]
