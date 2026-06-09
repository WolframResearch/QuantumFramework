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


(* ---- estimator primitive: observable layout + the IBMJob result accessors ----

   IBMJobSubmit[..., "Primitive" -> "estimator", "Observable" -> "ZZZ"] used to fail on hardware
   with IBM error 1501 ("number of qubits of the circuit (156) does not match ... the observable
   (3)"): the ISA circuit is transpiled to the backend's full physical register, but the bare
   Pauli string was submitted unchanged at logical width. The fix delegates submission to qiskit's
   EstimatorV2, which lays out AND widens the observable to the device register via the transpiled
   circuit's layout. The estimator result is the expectation value, decoded by qiskit (not bit
   samples), so the sampler-only accessors are inert and ExpectationValue(s) / StandardErrors are
   the live ones. *)

(* core mechanism, in the qiskit session (no IBM auth): a 3-qubit observable is widened to the
   device width by apply_layout on the transpiled circuit's layout *)
obsLayoutWidth[deviceQubits_] := Block[{pn = deviceQubits}, pe[Context[pn], "
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import SparsePauliOp
from qiskit.providers.fake_provider import GenericBackendV2
be = GenericBackendV2(num_qubits=<* pn *>, seed=1)
c = QuantumCircuit(3); c.h(0); c.cx(0, 1); c.cx(1, 2)
isa = transpile(c, be, optimization_level=1)
int(SparsePauliOp('ZZZ').apply_layout(isa.layout).num_qubits)
"]];
qtest[obsLayoutWidth[7], 7, "estimator-observable-widened-to-device"]

(* an estimator IBMJob exposes the qiskit-decoded expectation value + standard error; the value
   cached at Refresh under "EstimatorResult" is read without any service round-trip *)
estJob = Wolfram`QuantumFramework`IBMJob[<|
    "ID" -> "est", "Status" -> "Completed", "Backend" -> "sim", "Primitive" -> "estimator",
    "EstimatorResult" -> <|"EVs" -> {0.42}, "Stds" -> {0.05}|>
|>];

VerificationTest[estJob["ExpectationValue"], 0.42, TestID -> "estimator-expectation-value"]
VerificationTest[estJob["ExpectationValues"], {0.42}, TestID -> "estimator-expectation-values"]
VerificationTest[estJob["StandardErrors"], {0.05}, TestID -> "estimator-standard-errors"]
VerificationTest[estJob[], 0.42, TestID -> "estimator-default-is-expectation-value"]
VerificationTest[estJob["Primitive"], "estimator", TestID -> "estimator-primitive-label"]

(* sampler-only accessors are not applicable to an estimator job and say so, instead of decoding
   the (absent) bit samples into garbage *)
VerificationTest[estJob["Counts"], Missing["NotApplicable", "estimator"], TestID -> "estimator-counts-not-applicable"]
VerificationTest[estJob["Measurement"], Missing["NotApplicable", "estimator"], TestID -> "estimator-measurement-not-applicable"]
VerificationTest[estJob["Samples"], Missing["NotApplicable", "estimator"], TestID -> "estimator-samples-not-applicable"]

(* a multi-observable estimator keeps the full expectation-value list (no single-value unwrap) *)
VerificationTest[
    Wolfram`QuantumFramework`IBMJob[<|
        "ID" -> "est2", "Status" -> "Completed", "Primitive" -> "estimator",
        "EstimatorResult" -> <|"EVs" -> {0.1, -0.2}, "Stds" -> {0.01, 0.02}|>
    |>]["ExpectationValue"],
    {0.1, -0.2},
    TestID -> "estimator-multi-observable-expectation-values"
]

(* the estimator summary box renders without recursion / messages *)
VerificationTest[FreeQ[ToBoxes[estJob], TerminatedEvaluation], True, TestID -> "estimator-box-no-recursion"]

(* a job with no "Primitive" key defaults to sampler, so the existing sampler decode is unchanged *)
VerificationTest[
    Wolfram`QuantumFramework`IBMJob[<|"ID" -> "s", "Backend" -> "sim"|>]["Primitive"],
    "sampler",
    TestID -> "primitive-defaults-to-sampler"
]


(* ---- QiskitCircuit constructor: OpenQASM-string import, bad-argument guard, and the
   summary-box recursion regression ----

   A QiskitCircuit wraps a pickled qiskit circuit (a ByteArray). Three behaviors are pinned:

   (1) an OpenQASM string is imported into a real handle, so QiskitCircuit[QuantumQASM[qco]]
       round-trips a circuit (the import reuses QuantumCircuitOperator[source]["Qiskit"]);
   (2) a non-OpenQASM string is a clean $Failed with a message, never a live object;
   (3) a handle whose argument is not a ByteArray must NOT recurse when displayed. The
       summary-box formatter is restricted to QiskitCircuit[_ByteArray]; a malformed handle
       falls back to default boxes instead of re-entering MakeBoxes forever (which previously
       produced TerminatedEvaluation["RecursionLimit"] when an accessor re-embedded the object). *)

(* tolerance-aware unitary equivalence, up to global phase; atol matched to QuantumQASM's
   ~6-significant-figure angle literals, so the truncated angles in the QASM text do not
   register as a structural difference. *)
qcEquivTol[ba_, bb_] := Block[{pa = ba, pb = bb}, pe[Context[pa], "
import pickle
from qiskit.quantum_info import Operator
a = pickle.loads(<* pa *>); b = pickle.loads(<* pb *>)
bool(Operator(a).equiv(Operator(b), rtol=1e-4, atol=1e-4))
"]];

(* (1) the reported form QiskitCircuit[QuantumQASM[QCO[...]]] builds a real ByteArray handle *)
qtest[Head @ QiskitCircuit[QuantumQASM[QuantumCircuitOperator["Toffoli"]]], QiskitCircuit, "qasm-string-builds-handle"]

qtest[MatchQ[QiskitCircuit[QuantumQASM[qcU]], QiskitCircuit[_ByteArray]], True, "qasm-string-handle-is-bytearray"]

(* the round trip preserves the unitary (Toffoli and Bell), up to global phase and QASM precision *)
qtest[
    qcEquivTol[QuantumCircuitOperator["Toffoli"]["Qiskit"]["Bytes"], QiskitCircuit[QuantumQASM[QuantumCircuitOperator["Toffoli"]]]["Bytes"]],
    True, "qasm-string-roundtrip-toffoli-unitary"
]

qtest[qcEquivTol[qcU["Qiskit"]["Bytes"], QiskitCircuit[QuantumQASM[qcU]]["Bytes"]], True, "qasm-string-roundtrip-bell-unitary"]

(* OpenQASM 2.0 source imports as well as the default 3.0 *)
qtest[MatchQ[QiskitCircuit[QuantumQASM[qcU, "Version" -> 2]], QiskitCircuit[_ByteArray]], True, "qasm2-string-imports"]

(* displaying a freshly imported handle uses the summary box without recursing and without
   messages (the icon falls back to a graphic whenever $quantumServiceIcon is not an Image,
   so Show never receives a non-graphics symbol) *)
qtest[FreeQ[ToBoxes @ QiskitCircuit[QuantumQASM[qcU]], TerminatedEvaluation], True, "qasm-string-handle-displays"]

(* (2) a non-OpenQASM string is rejected cleanly with a message - pure WL, no qiskit needed *)
VerificationTest[
    QiskitCircuit["not an openqasm program"],
    $Failed, {QiskitCircuit::badarg},
    TestID -> "non-qasm-string-fails-with-message"
]

(* (3) recursion regression: a non-ByteArray handle must format without blowing the stack.
   Pure WL, no qiskit needed; would TerminatedEvaluation["RecursionLimit"] before the fix. *)
VerificationTest[
    Block[{$RecursionLimit = 200}, FreeQ[ToBoxes @ QiskitCircuit[12345], TerminatedEvaluation]],
    True, TestID -> "nonbytearray-handle-no-recursion"
]

(* the bad-argument guard keys on a literal String type, so it does not fire on pattern
   expressions: QiskitCircuit[...] patterns stay inert and match normally - pure WL *)
VerificationTest[
    MatchQ[QiskitCircuit[ByteArray[{1, 2, 3}]], QiskitCircuit[_ByteArray]],
    True, TestID -> "bytearray-handle-matches-pattern"
]

VerificationTest[
    Length @ Cases[{QiskitCircuit[ByteArray[{1}]], $Failed, QiskitCircuit[ByteArray[{2}]]}, QiskitCircuit[_ByteArray]],
    2, TestID -> "qiskitcircuit-pattern-not-shadowed-by-guard"
]


(* ================================================================================================
   Layer 5: QiskitCircuit <-> QuantumCircuitOperator gate-level fidelity.

   The earlier layers exercise the QASM-string direction and byte-level unitary equivalence, but
   never the QiskitCircuit -> QuantumCircuitOperator decode (qc_to_QuantumCircuitOperator). That is
   exactly where named operators are reconstructed, and where a decode could silently emit a
   degenerate Label -> None operator (controlled / parametrized gates routed through dead list
   forms) or a wrong-dimension operator ("SX" became a 2-qubit S(x)X in 2.0.0; sqrt(X) is "V"),
   while every existing test stayed green. This layer pins, per gate and per accessor:
     (a) the decoded circuit reproduces the source unitary up to global phase  (opEquivQ),
     (b) no decoded operator is a degenerate Label -> None state              (noNoneQ),
     (c) named structure survives (controlled gates stay controlled shortcuts),
     (d) the QiskitCircuit accessors report correct outcomes.
   ================================================================================================ *)

(* matrix equivalence up to a global phase and float tolerance *)
matEquivQ[a_, b_, tol_ : 1.*^-7] := Block[{m1 = N @ Normal @ a, m2 = N @ Normal @ b, idx, r},
    If[Dimensions[m1] =!= Dimensions[m2], Return[False]];
    idx = FirstPosition[m2, x_ /; Abs[x] > tol, None, {2}];
    If[idx === None, Return[Max[Abs[Flatten[m1]]] < tol]];
    r = Extract[m1, idx] / Extract[m2, idx];
    Max[Abs[Flatten[m1 - r m2]]] < tol
];
opEquivQ[x_, y_] := matEquivQ[x["Sort"]["Matrix"], y["Sort"]["Matrix"]];
qcoRT[qco_] := QuantumCircuitOperator[QiskitCircuit[qco]];                 (* decode round trip *)
noNoneQ[qco_] := FreeQ[#["Label"] & /@ qco["Flatten"]["Operators"], None]; (* no degenerate op *)
(* a clean decode that both preserves the unitary and leaves no degenerate operator *)
decodeCleanQ[qco_] := With[{b = qcoRT[qco]}, opEquivQ[b, qco] && noNoneQ[b]];

(* ---- 5a: single-qubit named gates decode to the right unitary, no degenerate operator ---- *)
qtest[decodeCleanQ[QuantumCircuitOperator[{"X" -> {1}}]], True, "decode-X"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"Y" -> {1}}]], True, "decode-Y"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"Z" -> {1}}]], True, "decode-Z"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"H" -> {1}}]], True, "decode-H"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"S" -> {1}}]], True, "decode-S"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"T" -> {1}}]], True, "decode-T"]
qtest[decodeCleanQ[QuantumCircuitOperator[{SuperDagger["T"] -> {1}}]], True, "decode-Tdg"]
qtest[decodeCleanQ[QuantumCircuitOperator[{SuperDagger["S"] -> {1}}]], True, "decode-Sdg"]

(* sqrt(X): qiskit sx must decode to QF "V" (sqrt(X)), not "SX" (which is a 2-qubit S(x)X in 2.0.0) *)
qtest[decodeCleanQ[QuantumCircuitOperator[{"V" -> {1}}]], True, "decode-V-sqrtX"]
qtest[QuantumCircuitOperator[QiskitCircuit[QuantumCircuitOperator[{"V" -> {1}}]]]["Sort"]["Dimensions"], {2, 2}, "decode-V-is-single-qubit"]
qtest[decodeCleanQ[QuantumCircuitOperator[{SuperDagger["V"] -> {1}}]], True, "decode-Vdg-sqrtXdg"]

(* ---- 5b: parametrized single-qubit gates decode through the "Name"[args] call form ---- *)
qtest[decodeCleanQ[QuantumCircuitOperator[{"RX"[0.7] -> {1}}]], True, "decode-RX"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"RY"[1.1] -> {1}}]], True, "decode-RY"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"RZ"[0.5] -> {1}}]], True, "decode-RZ"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"P"[0.9] -> {1}}]], True, "decode-P"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"U"[0.3, 0.4, 0.5] -> {1}}]], True, "decode-U3"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"U2"[0.3, 0.4] -> {1}}]], True, "decode-U2"]
(* the symptom that exposed the list-form rot: a None-labeled operator must never appear *)
qtest[
    noNoneQ @ QuantumCircuitOperator[QiskitCircuit[QuantumCircuitOperator[{"RX"[0.7] -> {1}, "RZ"[0.5] -> {1}, "U"[0.3, 0.4, 0.5] -> {2}}]]],
    True, "decode-parametrized-no-none-label"
]

(* ---- 5c: SWAP and controlled gates (1-control, open control, multi-control, controlled-param) ---- *)
qtest[decodeCleanQ[QuantumCircuitOperator[{"SWAP" -> {1, 2}}]], True, "decode-SWAP"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"CNOT" -> {1, 2}}]], True, "decode-CNOT"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"CZ" -> {1, 2}}]], True, "decode-CZ"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"C0"["NOT" -> 2] -> {1, 2}}]], True, "decode-open-control-C0NOT"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"C"["R"[0.7, "X"] -> 2] -> {1, 2}}]], True, "decode-controlled-RX"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"C"["V" -> 2] -> {1, 2}}]], True, "decode-controlled-V"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"C"["Z" -> 3, {1}, {2}] -> {1, 2, 3}}]], True, "decode-mixed-control-CCZ"]

(* ---- 5d: the originally reported case, end to end ---- *)
qtest[decodeCleanQ[QuantumCircuitOperator["Toffoli"]], True, "decode-Toffoli-clean"]
(* the exact reported symptom: QuantumShortcut yields clean shortcuts, no degenerate element *)
qtest[FreeQ[QuantumShortcut[qcoRT[QuantumCircuitOperator["Toffoli"]]], None], True, "decode-Toffoli-no-none-shortcut"]
(* named structure survives: the standard Toffoli decomposition has exactly six controlled-NOTs,
   each a "C"[...] shortcut rather than a raw unitary block *)
qtest[Length @ Cases[QuantumShortcut[qcoRT[QuantumCircuitOperator["Toffoli"]]], "C"[___]], 6, "decode-Toffoli-six-controlled-nots"]

(* ---- 5e: composite circuits ---- *)
qtest[decodeCleanQ[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]], True, "decode-Bell"]
qtest[decodeCleanQ[QuantumCircuitOperator["Fourier"[4]]], True, "decode-QFT4"]
qtest[decodeCleanQ[QuantumCircuitOperator[{"H" -> 1, "RX"[0.4] -> 2, "CNOT" -> {1, 2}, SuperDagger["T"] -> 2, "CZ" -> {2, 3}}]], True, "decode-mixed-circuit"]
qtest[opEquivQ[qcoRT[QuantumCircuitOperator[{"Permutation"[Cycles[{{1, 2, 3}}]] -> {1, 2, 3}}]], QuantumCircuitOperator[{"Permutation"[Cycles[{{1, 2, 3}}]] -> {1, 2, 3}}]], True, "decode-Permutation"]

(* ---- 5f: measurement decode structure (a measured circuit decodes to a QuantumMeasurementOperator) ---- *)
qtest[
    ! FreeQ[QuantumCircuitOperator[QiskitCircuit[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}]]]["Elements"], _QuantumMeasurementOperator],
    True, "decode-measurement-structure"
]

(* ---- 5g: QiskitCircuit accessors report correct outcomes ---- *)
qtest[QiskitCircuit[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["Qubits"], 2, "accessor-qubits"]
qtest[QiskitCircuit[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["Depth"], 2, "accessor-depth"]
qtest[Lookup[QiskitCircuit[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["Ops"], "cx"], 1, "accessor-ops-cx-count"]
(* the qiskit-side unitary (Operator) equals the source unitary up to global phase *)
qtest[matEquivQ[QiskitCircuit[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["QuantumOperator"]["Matrix"], QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]["Sort"]["Matrix"]], True, "accessor-quantumoperator-unitary"]
qtest[matEquivQ[QiskitCircuit[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["Matrix"], QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]["Sort"]["Matrix"]], True, "accessor-matrix-unitary"]
(* Decompose keeps the unitary invariant *)
qtest[matEquivQ[QiskitCircuit[QuantumCircuitOperator["Toffoli"]]["Decompose"]["Matrix"], QuantumCircuitOperator["Toffoli"]["Sort"]["Matrix"]], True, "accessor-decompose-unitary"]

(* ---- 5h: scale / outcome on a larger circuit (8x8 ... 64x64 unitaries) ---- *)
qtest[decodeCleanQ[QuantumCircuitOperator["Fourier"[6]]], True, "decode-QFT6-scale"]

EndTestSection[]
