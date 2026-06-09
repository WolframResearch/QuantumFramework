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

qtest[QuantumQASM[qc, "Version" -> 3], gtQASM[qc["Qiskit"]["Bytes"], <||>, 3], "faithful-qasm3-matches-qiskit"]

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
    With[{s = QuantumQASM[qc, "Version" -> 3]}, StringContainsQ[s, "h q[0]"] && StringContainsQ[s, "cx q[0], q[1]"]],
    True, "qiskit-faithful-named-gates"
]

qtest[StringStartsQ[QuantumQASM[qc, "Version" -> 2], "OPENQASM 2.0"], True, "qasm2-prefix"]

(* bare QuantumQASM (no args) is the native, dependency-free WL emission *)
qtest[StringStartsQ[QuantumQASM[qc], "OPENQASM 3.0"], True, "native-default-prefix"]

(* a Provider / Backend request engages the hardware-provider export path (transpile against
   a real backend) in qiskitQASM, which lives in Qiskit.m and is reached cross-file from the
   QuantumQASM hub: the guard pins that cross-file reach so the branch cannot silently fall
   back to an unevaluated symbol. "Backend" -> Automatic uses the local Aer/BasicProvider, so
   no credentials are needed; the measured qc avoids the save_statevector instruction that a
   pure-unitary circuit on Aer would emit (qasm3 cannot export it). *)
qtest[StringStartsQ[QuantumQASM[qc, "Backend" -> Automatic], "OPENQASM 3.0;"], True, "provider-backend-path-exports-qasm3"]


(* ---- qubit-only guard: OpenQASM models qubits only, so any qudit dimension > 2 is a
   clear Failure rather than an unevaluated emitter call or an opaque ConfirmBy error.
   These are pure-native checks (no qiskit needed), so they use VerificationTest directly. *)

VerificationTest[
    FailureQ @ QuantumQASM[QuantumOperator["X"[3]]],
    True, TestID -> "qudit-operator-failure"
]

VerificationTest[
    FailureQ @ QuantumQASM[QuantumOperator["X"[3]], "Simple"],
    True, TestID -> "qudit-operator-simple-failure"
]

VerificationTest[
    FailureQ @ QuantumQASM[QuantumCircuitOperator[{"X"[3] -> 1}]],
    True, TestID -> "qudit-circuit-failure"
]

VerificationTest[
    FailureQ @ QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "X"[3] -> 2}], "WL"],
    True, TestID -> "mixed-qubit-qudit-circuit-failure"
]

VerificationTest[
    QuantumQASM[QuantumCircuitOperator[{"X"[3] -> 1}]]["NonQubitDimensions"],
    {3}, TestID -> "qudit-failure-reports-dimensions"
]

(* the guard must not reject genuine qubit operators/circuits *)
VerificationTest[
    StringQ @ QuantumQASM[QuantumOperator["X"]],
    True, TestID -> "qubit-operator-still-emits"
]


(* ---- joint measurements: a joint multi-qubit measurement carries a single classical
   outcome register whose dimension is the product of the measured qubits' dimensions (2^n
   for n qubits). That register is emitted as OpenQASM bits, not a quantum wire, so the guard
   tests the quantum wires of each operator (positive orders) and a joint measurement on
   qubits exports exactly like the equivalent separate single-qubit measurements. These are
   pure-native checks (no qiskit needed). *)

VerificationTest[
    StringQ @ QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1, 2}}]],
    True, TestID -> "joint-measurement-exports"
]

VerificationTest[
    QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1, 2}}]] ===
        QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}]],
    True, TestID -> "joint-measurement-matches-separate"
]

(* a genuine non-qubit quantum wire is still rejected when a joint measurement also adds a
   2^n outcome register: only the qutrit dimension is reported, not the classical register. *)
VerificationTest[
    QuantumQASM[QuantumCircuitOperator[{"H"[3] -> 1, "H" -> 2, "H" -> 3, {1, 2, 3}}]]["NonQubitDimensions"],
    {3}, TestID -> "qutrit-with-joint-measurement-reports-only-qutrit"
]

(* the wire dimensions are read per operator, so a non-qubit wire whose dimension is masked
   in the circuit's aggregate OutputDimensions by a downstream measurement is still caught. *)
VerificationTest[
    QuantumQASM[QuantumCircuitOperator[{"H"[3] -> 1, "H" -> 2, {1, 2, 3}}]]["NonQubitDimensions"],
    {3}, TestID -> "qutrit-masked-by-measurement-still-caught"
]


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


(* ---- IBM sample decoding: classical-bit -> ascending-qubit reordering ----

   IBM / qiskit return each shot in classical-bit order, which follows the measurement's
   emission (target) order; QuantumMeasurement reports outcomes in ascending qubit order. A
   non-ascending (permuted) measurement, or any backend whose transpiled measure -> clbit map
   is not the identity, must be reordered with the per-clbit original-qubit map captured at
   submission time. These are pure-WL checks (no qiskit / no credentials needed) so they run
   in CI and lock the decode against silent mis-ordering. *)

reorder = Wolfram`QuantumFramework`PackageScope`qiskitReorderByQubit;

VerificationTest[reorder[{1, 0, 1}, {0, 1, 2}], {1, 0, 1}, TestID -> "reorder-identity"]
VerificationTest[reorder[{1, 1, 0}, {2, 0, 1}], {1, 0, 1}, TestID -> "reorder-permuted"]
VerificationTest[reorder[{1, 1, 0}, {2, 0}],    {1, 1, 0}, TestID -> "reorder-length-mismatch-noop"]
VerificationTest[reorder[{1, 1, 0}, Automatic], {1, 1, 0}, TestID -> "reorder-no-map-noop"]

(* synthetic Completed IBMJob carrying raw hex samples + the captured map: the permuted
   measurement {3,1,2} of X@1, X@3 (qubit values q1=1, q2=0, q3=1) is read into classical
   bits in target order [q3,q1,q2] = [1,1,0] (integer 0x3), and the map {2,0,1} must reorder
   it to the ascending-qubit outcome {1,0,1} that QuantumCircuitOperator[...][] gives. *)
ibmRawJob[map_, hex_, nbits_] := Wolfram`QuantumFramework`IBMJob[<|
    "ID" -> "regression", "Status" -> "Completed", "Backend" -> "sim", "MeasuredQubits" -> map,
    "Raw" -> <|"Results" -> <|"results" -> {<|"data" -> <|"c" -> <|"samples" -> {hex}, "num_bits" -> nbits|>|>|>}|>|>
|>]

VerificationTest[
    Keys @ Normal @ ibmRawJob[{2, 0, 1}, "0x3", 3]["Counts"],
    {{1, 0, 1}},
    TestID -> "ibm-decode-permuted-reordered"
]

VerificationTest[
    First @ Keys @ Normal @ ibmRawJob[{2, 0, 1}, "0x3", 3]["Counts"] ===
        Replace[First @ Keys @ Select[QuantumCircuitOperator[{"X" -> 1, "X" -> 3, {3, 1, 2}}][]["Probabilities"], # > 1/2 &], QuditName[v_, ___] :> v],
    True,
    TestID -> "ibm-decode-matches-exact-distribution"
]

(* without a map the decode falls back to classical-bit order: the old, unreordered outcome *)
VerificationTest[
    First @ Keys @ Normal @ ibmRawJob[Automatic, "0x3", 3]["Counts"],
    {1, 1, 0},
    TestID -> "ibm-decode-no-map-is-clbit-order"
]

(* identity map (ascending full measurement) is unchanged: q1=1,q2=0,q3=1 -> 0x5 -> {1,0,1} *)
VerificationTest[
    First @ Keys @ Normal @ ibmRawJob[{0, 1, 2}, "0x5", 3]["Counts"],
    {1, 0, 1},
    TestID -> "ibm-decode-identity-unchanged"
]

EndTestSection[]
