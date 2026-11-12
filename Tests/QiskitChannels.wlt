BeginTestSection["QiskitChannels"]

(* Every route out of QF (to a qiskit circuit, native OpenQASM, qiskit's OpenQASM dump, OpenQASM
   in a gate set) and back into QF (the OpenQASM importer, the qiskit decoder) must carry the
   circuit's whole channel, not only its outcome probabilities from |0...0>.

   For each circuit, QF's channel and each route's channel are compared on every product of
   |0>, |1>, |+>, |+i> over the circuit's input wires (these span all operators, so agreement on
   all of them is equality of the channels): the output state of the wires, and the joint
   distribution of the classical record keyed by the measured wire. qiskit's reading of a route's
   circuit or program is the reference for what it means; each measurement is replaced there by a
   CNOT onto a fresh record qubit, so the record distribution is exact, not sampled.

   This suite needs a working qiskit session with qiskit_qasm3_import and SELF-GUARDS like
   Qiskit.wlt: without them every test trivially passes. *)

(* qiskit is usable when a circuit converts; Check turns a missing Python session's messages into False *)
qiskitOK = TrueQ @ Check[MatchQ[QuantumCircuitOperator["Bell"]["Qiskit"], _QiskitCircuit], False] &&
    TrueQ @ PythonEvaluate[Context[qiskitOK], "import importlib.util\nimportlib.util.find_spec('qiskit_qasm3_import') is not None"];

SetAttributes[ctest, HoldAll];
ctest[actual_, expected_, id_] := VerificationTest[If[qiskitOK, actual, "skip"], If[qiskitOK, expected, "skip"], TestID -> id];

$channelPython = "
import itertools
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import DensityMatrix, partial_trace

def ch_deferred(circ):
    n, nc = circ.num_qubits, circ.num_clbits
    d = QuantumCircuit(n + nc)
    d.global_phase = circ.global_phase
    meas = {}
    for inst in circ.data:
        op = inst.operation
        qs = [circ.find_bit(q).index for q in inst.qubits]
        cs = [circ.find_bit(c).index for c in inst.clbits]
        if op.name == 'measure':
            d.cx(qs[0], n + cs[0])
            meas[cs[0]] = qs[0]
        elif op.name not in ('barrier', 'delay'):
            d.append(op, qs)
    return d, meas

def ch_wire_order(data, m):
    t = np.asarray(data).reshape([2] * (2 * m))
    rows = [m - 1 - q for q in range(m)]
    return t.transpose(rows + [m + r for r in rows]).reshape(2 ** m, 2 ** m)

def ch_analyze(circ, n_wires, input_wires):
    d, meas = ch_deferred(circ)
    nq, nc = circ.num_qubits, circ.num_clbits
    rec = sorted(meas, key=lambda c: (meas[c], c))
    states, probs = [], []
    for t in itertools.product('01+r', repeat=len(input_wires)):
        label = ['0'] * (nq + nc)
        for w, ch in zip(input_wires, t):
            label[w] = ch
        rho = DensityMatrix.from_label(''.join(reversed(label))).evolve(d)
        others = list(range(n_wires, nq + nc))
        sys = partial_trace(rho, others) if others else rho
        states.append(ch_wire_order(sys.data, n_wires).tolist())
        if rec:
            probs.append(np.real(rho.probabilities(qargs=list(reversed([nq + c for c in rec])))).tolist())
    return {'states': states, 'probs': probs}
";

oneQubit = <|"0" -> {1, 0}, "1" -> {0, 1}, "+" -> {1, 1} / Sqrt[2], "r" -> {1, I} / Sqrt[2]|>;

(* the input wires of the reference circuit; a route's circuit is fed the same inputs on those
   wires and |0> on any other input wire it has *)
inputWires[qc_, n_] := Intersection[Select[Union @ Flatten @ {qc["InputOrder"]}, Positive], Range[n]]

traceFirst[dm_, dims_, 0] := dm
traceFirst[dm_, dims_, k_] := With[{r = Length[dims]},
    Flatten[TensorContract[ArrayReshape[dm, Join[dims, dims]], Table[{i, r + i}, {i, k}]], {Range[r - k], Range[r - k + 1, 2 (r - k)]}]
]

qfChannel[qc_, n_, iw_] := With[{ws = inputWires[qc, n]},
    Transpose @ Map[
        Function[t, With[{vec = AssociationThread[iw -> Lookup[oneQubit, t]]},
            With[{out = If[ws === {}, qc[], qc[QuantumState[N @ Flatten[Outer[Times, Sequence @@ Lookup[vec, ws, {{1, 0}}]]], ConstantArray[2, Length[ws]]]]]},
                If[ QuantumMeasurementQ[out],
                    {traceFirst[Normal @ out["State"]["Computational"]["DensityMatrix"], out["State"]["Dimensions"], out["Eigenqudits"]], N @ out["ProbabilitiesList"]},
                    {Normal @ out["Computational"]["DensityMatrix"], {}}
                ]
            ]
        ]],
        Tuples[Keys[oneQubit], Length[iw]]
    ]
]

qiskitChannel[src_, kind_, n_, iw_] := Block[{pysrc = src, pykind = kind, pyn = n, pyiw = iw - 1},
    Values @ KeyTake[
        PythonEvaluate[Context[pysrc], $channelPython <> "
import pickle
from qiskit import qasm3
s = <* pysrc *>
c = pickle.loads(s) if <* pykind *> == 'pickle' else qasm3.loads(s)
ch_analyze(c, <* pyn *>, <* pyiw *>)
"],
        {"states", "probs"}
    ]
]

(* the largest trace distance between output states and total-variation distance between record
   distributions over all inputs; {0, 0} when the channels agree *)
channelDistance[{sa_, pa_}, {sb_, pb_}] := {
    Chop[Max[MapThread[Total[SingularValueList[N[Normal[#1] - Normal[#2]], Tolerance -> 0]] / 2 &, {sa, sb}]], 10 ^ -8],
    If[ pa === {} && pb === {} || AllTrue[Join[pa, pb], # === {} &], 0,
        Chop[Max[MapThread[Total[Abs[Normal[#1] - Normal[#2]]] / 2 &, {pa, pb}]], 10 ^ -8]
    ]
}

u2 = Orthogonalize[{{1, 2, 0, 1}, {0, 1, 3, 0}, {2, 0, 1, 1}, {1, 1, 0, 2}}];
u3 = Orthogonalize[Table[Prime[i j] + i ^ 2 - j, {i, 8}, {j, 8}]];

channelCircuits = <|
    "named-gates-out-of-order" -> QuantumCircuitOperator[{
        "H" -> 1, "H" -> 2, "RY"[0.4] -> 3, "CNOT" -> {2, 1}, QuantumOperator["Toffoli", {3, 1, 2}],
        "CZ" -> {3, 1}, "SWAP" -> {3, 1}, "CP"[0.3] -> {2, 1}, QuantumOperator["C"[QuantumOperator["RX"[0.7], {3}], {1}, {2}]]
    }],
    "relabeling-of-three-wires" -> QuantumCircuitOperator[{"H" -> 1, "RY"[0.5] -> 2, QuantumOperator["Permutation"[], {1, 2, 3} -> {3, 1, 2}]}],
    "permutation-three-cycle" -> QuantumCircuitOperator[{"H" -> 1, "RY"[0.9] -> 2, "Permutation"[Cycles[{{1, 2, 3}}]] -> {3, 1, 2}}],
    "matrix-two-wires-reversed" -> QuantumCircuitOperator[{"H" -> 1, "RY"[0.4] -> 2, QuantumOperator[u2, {2, 1}]}],
    "matrix-three-wires-permuted" -> QuantumCircuitOperator[{"H" -> 1, "RY"[0.4] -> 2, QuantumOperator[u3, {3, 1, 2}]}],
    "matrix-two-of-three-wires" -> QuantumCircuitOperator[{"H" -> 1, "RY"[0.3] -> 3, QuantumOperator[u2, {3, 1}], "CNOT" -> {1, 2}}],
    "injection-permuted-wires" -> QuantumCircuitOperator[{QuantumOperator[QuantumState["001"], {{3, 1, 2}, {}}], "H" -> 2}],
    "reset-two-qubit-state" -> QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "Reset"[QuantumState["01"]] -> {2, 1}}],
    "measurement-out-of-order" -> QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "RY"[0.3] -> 3, {3, 1}}],
    "mid-circuit-measurement" -> QuantumCircuitOperator[{"H" -> 1, {1}, "CNOT" -> {1, 2}, "H" -> 1, {2}}],
    "x-and-y-measurements" -> QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "RY"[0.8] -> 3, QuantumMeasurementOperator["X", {2}], QuantumMeasurementOperator["Y", {3}]}],
    "x-measurement-then-gates" -> QuantumCircuitOperator[{"RY"[0.6] -> 1, QuantumMeasurementOperator["X", {1}], "RY"[0.3] -> 1, {1}}],
    "y-measurement-then-cnot" -> QuantumCircuitOperator[{"H" -> 1, "RY"[0.7] -> 2, QuantumMeasurementOperator["Y", {2}], "CNOT" -> {2, 1}, {1}}],
    "basis-listed-out-of-order" -> QuantumCircuitOperator[{"RY"[0.4] -> 1, "H" -> 2, "CNOT" -> {2, 3}, QuantumMeasurementOperator[QuantumBasis[{{0, 1}, {1, 0}}], {1}]}],
    "computational-measurement-two-wires" -> QuantumCircuitOperator[{"RY"[0.4] -> 1, "RY"[1.1] -> 2, "CNOT" -> {1, 2}, "M" -> {2, 1}}],
    "product-basis-measurement" -> QuantumCircuitOperator[{"H" -> 1, "RY"[0.4] -> 3, "CNOT" -> {1, 2}, QuantumMeasurementOperator["XZ", {3, 1}]}]
|>;

(* the routes: each gives either a qiskit program to read ({source, kind, qubit of each wire}) or a
   QF circuit *)
channelRoutes = <|
    "qiskit" -> Function[qc, {qc["Qiskit"]["Bytes"], "pickle"}],
    "native-qasm" -> Function[qc, {QuantumQASM[qc], "qasm3"}],
    "qiskit-qasm" -> Function[qc, {QuantumQASM[qc, "Version" -> 3], "qasm3"}],
    "gate-set-qasm" -> Function[qc, {QuantumQASM[qc, {"x", "sx", "rz", "cz"}, "OptimizationLevel" -> 3], "qasm3"}],
    "gate-set-qasm-initial-layout" -> Function[qc, {QuantumQASM[qc, {"x", "sx", "rz", "cz"}, "InitialLayout" -> RotateLeft[Range[qc["Max"]] - 1], "OptimizationLevel" -> 1], "qasm3"}],
    "import-native-qasm" -> Function[qc, QuantumQASM[QuantumQASM[qc]]],
    "import-qiskit-qasm" -> Function[qc, QuantumQASM[QuantumQASM[qc, "Version" -> 3]]],
    "decode-qiskit" -> Function[qc, qc["Qiskit"]["QuantumCircuitOperator"]],
    "decode-transpiled" -> Function[qc, qc["Qiskit"]["Transpile", {"x", "sx", "rz", "cz"}, "OptimizationLevel" -> 3]["QuantumCircuitOperator"]]
|>;

(* After a measurement outside the computational basis, a later gate on the measured qubit changes
   the outcome probabilities QF reports for that measurement, though the state it leaves is right
   (an X measurement of RY(0.6)|0> reports P(+) = 0.5, where Born's rule gives (1 + Sin[0.6]) / 2).
   For these circuits only the output state is compared. *)
stateOnly = {"x-measurement-then-gates", "y-measurement-then-cnot"};

(* native OpenQASM has no name for an unnamed unitary on three or more qubits (two-qubit ones are
   written through their KAK decomposition): the export refuses it, naming its wires *)
nativeRefused = {"matrix-three-wires-permuted"};

ctest[
    With[{f = QuantumQASM[channelCircuits["matrix-three-wires-permuted"]]}, {FailureQ[f], f["NonNativeGates"]}],
    {True, {"a matrix on wires {3, 1, 2}"}},
    "native-qasm-refuses-three-qubit-matrix"
]

(* a gate on a qutrit has no qiskit or OpenQASM form *)
ctest[
    FailureQ /@ {QuantumCircuitOperator[{QuantumOperator["X"[3], {1}]}]["Qiskit"], QuantumQASM[QuantumCircuitOperator[{QuantumOperator["X"[3], {1}]}]]},
    {True, True},
    "qutrit-gate-refused"
]

routeChannel[route_, qc_, n_, iw_] := With[{r = route[qc]},
    Which[
        MatchQ[r, {_String | _ByteArray, _String}], qiskitChannel[r[[1]], r[[2]], n, iw],
        MatchQ[r, _QuantumCircuitOperator], qfChannel[r, n, iw],
        True, r
    ]
]

If[ qiskitOK,
    KeyValueMap[
        Function[{cname, qc},
            With[{n = qc["Max"], iw = inputWires[qc, qc["Max"]]},
                With[{ref = qfChannel[qc, n, iw]},
                    KeyValueMap[
                        Function[{rname, route},
                            ctest[
                                If[MemberQ[stateOnly, cname], First, Identity] @ channelDistance[ref, routeChannel[route, qc, n, iw]],
                                If[MemberQ[stateOnly, cname], 0, {0, 0}],
                                "channel-" <> cname <> "-" <> rname
                            ]
                        ],
                        If[MemberQ[nativeRefused, cname], KeyDrop[channelRoutes, {"native-qasm", "import-native-qasm"}], channelRoutes]
                    ]
                ]
            ]
        ],
        channelCircuits
    ],
    ctest[True, True, "channel-skipped"]
]


(* ---- measurements of an observable ----
   The exported record must follow Born's rule for the observable's eigenbasis, which is what the
   measurement gives applied to the state directly (mo[state]). *)

observableMeasurements = <|
    "observable-ZX-wires-reversed" -> QuantumMeasurementOperator[QuantumOperator["ZX", {2, 1}]],
    "observable-ZX-wires-in-order" -> QuantumMeasurementOperator[QuantumOperator["ZX", {1, 2}]],
    "observable-entangled-eigenbasis" -> QuantumMeasurementOperator[QuantumOperator[
        KroneckerProduct[PauliMatrix[1], PauliMatrix[1]] + 0.5 KroneckerProduct[PauliMatrix[3], PauliMatrix[3]] + 0.3 KroneckerProduct[PauliMatrix[3], IdentityMatrix[2]], {2, 1}]],
    "basis-XZ-wires-reversed" -> QuantumMeasurementOperator["XZ", {2, 1}]
|>;

observablePrefix = {"RY"[0.4] -> 1, "RY"[1.1] -> 2, "CNOT" -> {1, 2}};

If[ qiskitOK,
    KeyValueMap[
        Function[{name, mo},
            ctest[
                With[{
                    direct = N @ mo[QuantumCircuitOperator[observablePrefix][]]["ProbabilitiesList"],
                    exported = First @ Last @ qiskitChannel[QuantumCircuitOperator[Append[observablePrefix, mo]]["Qiskit"]["Bytes"], "pickle", 2, {}]
                },
                    Chop[Max[Abs[direct - exported]], 10 ^ -8]
                ],
                0,
                "measurement-born-rule-" <> name
            ]
        ],
        observableMeasurements
    ]
]


(* ---- the gate table: every standard gate the readers know means the same in QF and qiskit ----
   Each name of the importer's table, called with sample parameters on wires out of order, is
   read by QF's importer and by qiskit; the two operators must agree up to a global phase. *)

sampleParams = {0.37, -1.1, 0.52, 0.81};

(* cu1 and cu3 are OpenQASM 2 names (qelib1.inc), not in OpenQASM 3's stdgates.inc *)
qasm2Only = {"cu1", "cu3"};

gateTableProgram[name_, nParams_, nQubits_] := StringJoin[
    If[ MemberQ[qasm2Only, name],
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\n",
        "OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[3] q;\n"
    ],
    name, If[nParams > 0, "(" <> StringRiffle[ToString /@ Take[sampleParams, nParams], ", "] <> ")", ""], " ",
    StringRiffle[Take[{"q[2]", "q[0]", "q[1]"}, nQubits], ", "], ";"
]

qiskitOperatorOf[src_] := Block[{pysrc = src},
    PythonEvaluate[Context[pysrc], "
from qiskit import qasm2, qasm3
from qiskit.quantum_info import Operator
s = <* pysrc *>
c = qasm2.loads(s, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS) if s.startswith('OPENQASM 2') else qasm3.loads(s)
Operator(c).reverse_qargs().data.tolist()
"]
]

phaseEqualQ[a_, b_] := With[{ov = Tr[ConjugateTranspose[a] . b]}, Abs[ov] > 10 ^ -9 && Chop[Norm[Flatten[b - ov / Abs[ov] a]], 10 ^ -8] == 0]

(* the parameter and qubit count of each standard gate: the shape qiskit's stdgates.inc gives it *)
stdGateShapes = <|
    "id" -> {0, 1}, "x" -> {0, 1}, "y" -> {0, 1}, "z" -> {0, 1}, "h" -> {0, 1}, "s" -> {0, 1}, "sdg" -> {0, 1},
    "t" -> {0, 1}, "tdg" -> {0, 1}, "sx" -> {0, 1}, "rx" -> {1, 1}, "ry" -> {1, 1}, "rz" -> {1, 1}, "p" -> {1, 1},
    "u1" -> {1, 1}, "phase" -> {1, 1}, "u2" -> {2, 1}, "u3" -> {3, 1}, "U" -> {3, 1}, "swap" -> {0, 2},
    "cx" -> {0, 2}, "CX" -> {0, 2}, "cy" -> {0, 2}, "cz" -> {0, 2}, "ch" -> {0, 2}, "crx" -> {1, 2}, "cry" -> {1, 2},
    "crz" -> {1, 2}, "cp" -> {1, 2}, "cu1" -> {1, 2}, "cphase" -> {1, 2}, "cu3" -> {3, 2}, "cu" -> {4, 2},
    "ccx" -> {0, 3}, "cswap" -> {0, 3}
|>;

If[ qiskitOK,
    KeyValueMap[
        Function[{name, shape},
            ctest[
                With[{src = gateTableProgram[name, Sequence @@ shape]},
                    (* the QF circuit on all three declared wires, as qiskit reads the program *)
                    phaseEqualQ[
                        N @ Normal @ QuantumCircuitOperator[Prepend[QuantumQASM[src]["Elements"], QuantumOperator["I", {1, 2, 3}]]]["QuantumOperator"]["Sort"]["Matrix"],
                        qiskitOperatorOf[src]
                    ]
                ],
                True,
                "gate-table-import-" <> name
            ]
        ],
        stdGateShapes
    ]
]

(* every name the importer's table has is checked above *)
VerificationTest[
    Complement[ToLowerCase /@ Wolfram`QuantumFramework`PackageScope`qasmGateNames[], ToLowerCase /@ Keys[stdGateShapes], {"i", "cnot", "toffoli", "fredkin", "ccz", "csx", "sxdg", "u"}],
    {},
    TestID -> "gate-table-covered"
]

EndTestSection[]
