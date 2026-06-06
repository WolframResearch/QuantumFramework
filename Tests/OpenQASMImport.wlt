BeginTestSection["OpenQASMImport"]

(* Native Wolfram-Language OpenQASM importer. This suite is PURE WL: it needs no qiskit /
   Python session, so it runs everywhere (a strict improvement over the qiskit-coupled
   QuantumQASM export tests). The importer is exercised through the public surface:
   QuantumQASM[source] and QuantumCircuitOperator[source]. *)

(* matrix-equivalence helpers *)
mnorm[a_, b_] := Norm[Flatten[Normal[N[a["Matrix"]]] - Normal[N[b["Matrix"]]]]];
exactQ[a_, b_] := TrueQ[Chop[mnorm[a, b]] == 0];
closeQ[a_, b_] := TrueQ[mnorm[a, b] < 1.*^-4];

bell = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}];

(* ---- basic parse: v2 and v3, with measurement ---- *)

VerificationTest[
    Head @ QuantumQASM["OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\ncreg c[2];\nh q[0];\ncx q[0],q[1];\nmeasure q[0] -> c[0];\nmeasure q[1] -> c[1];"],
    QuantumCircuitOperator,
    TestID -> "v2-bell-measure-parses"
]

VerificationTest[
    Head @ QuantumQASM["OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[2] q;\nbit[2] c;\nh q[0];\ncx q[0], q[1];\nc[0] = measure q[0];\nc[1] = measure q[1];"],
    QuantumCircuitOperator,
    TestID -> "v3-bell-measure-parses"
]

(* ---- v2 / v3 parity: same Bell unitary in both dialects ---- *)

VerificationTest[
    exactQ[
        QuantumQASM["OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\nh q[0];\ncx q[0],q[1];"],
        QuantumQASM["OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[2] q;\nh q[0];\ncx q[0], q[1];"]
    ],
    True,
    TestID -> "v2-v3-parity"
]

(* ---- exact round-trip on hand-written standard QASM ---- *)

VerificationTest[
    exactQ[
        QuantumQASM["OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[2] q;\nh q[0];\ncx q[0], q[1];\nrz(pi/3) q[1];"],
        QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "RZ"[Pi/3] -> 2}]
    ],
    True,
    TestID -> "exact-roundtrip-standard"
]

(* ---- emitter round-trip (tolerance: the WL emitter rounds angles to ~6 digits) ---- *)

VerificationTest[
    closeQ[QuantumCircuitOperator[bell["QASM"]], bell],
    True,
    TestID -> "emitter-roundtrip"
]

(* ---- custom gate definition -> named subcircuit ---- *)

VerificationTest[
    exactQ[
        QuantumQASM["OPENQASM 3.0;\ninclude \"stdgates.inc\";\ngate mygate a, b { h a; cx a, b; }\nqubit[2] q;\nmygate q[0], q[1];"],
        bell
    ],
    True,
    TestID -> "custom-gate-def"
]

(* ---- named controlled gates ---- *)

VerificationTest[
    exactQ[QuantumQASM["OPENQASM 3.0;\nqubit[2] q;\ncx q[0],q[1];"], QuantumCircuitOperator[{"CNOT" -> {1, 2}}]],
    True,
    TestID -> "named-cx"
]

VerificationTest[
    exactQ[
        QuantumQASM["OPENQASM 3.0;\nqubit[3] q;\nccx q[0],q[1],q[2];"],
        QuantumCircuitOperator[{QuantumOperator["C"["X", {1, 2}, {}], {3}]}]
    ],
    True,
    TestID -> "named-ccx-toffoli"
]

(* ---- gate modifiers ---- *)

VerificationTest[
    exactQ[
        QuantumQASM["OPENQASM 3.0;\nqubit[2] q;\ncx q[0],q[1];"],
        QuantumQASM["OPENQASM 3.0;\nqubit[2] q;\nctrl @ x q[0],q[1];"]
    ],
    True,
    TestID -> "ctrl-modifier-equals-named"
]

VerificationTest[
    exactQ[
        QuantumQASM["OPENQASM 3.0;\nqubit[1] q;\ninv @ s q[0];"],
        QuantumQASM["OPENQASM 3.0;\nqubit[1] q;\nsdg q[0];"]
    ],
    True,
    TestID -> "inv-modifier-equals-sdg"
]

VerificationTest[
    Length @ QuantumQASM["OPENQASM 3.0;\nqubit[1] q;\npow(2) @ x q[0];"]["Operators"],
    2,
    TestID -> "pow-modifier-repeats"
]

(* ---- multi-register ordering ---- *)

VerificationTest[
    QuantumQASM["OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg a[2];\nqreg b[1];\ncx a[1], b[0];"]["Operators"][[1]]["Order"],
    {{2, 3}, {2, 3}},
    TestID -> "multi-register-order"
]

(* ---- reset / barrier / gphase parse ---- *)

VerificationTest[
    Head @ QuantumQASM["OPENQASM 3.0;\nqubit[2] q;\nh q[0];\nbarrier q;\nreset q[1];\ngphase(pi/4);"],
    QuantumCircuitOperator,
    TestID -> "reset-barrier-gphase"
]

(* ---- parametrized gate with a pi expression ---- *)

VerificationTest[
    QuantumQASM["OPENQASM 3.0;\nqubit[1] q;\nrx(pi/2) q[0];"]["Operators"][[1]]["Label"],
    Subscript["R", "X"][Pi/2],
    TestID -> "param-pi-expression"
]

(* ---- failure categories ---- *)

VerificationTest[
    QuantumQASM["OPENQASM 3.0;\nqubit[1] q;\nfoo q[0];"]["Category"],
    "UnknownGate",
    TestID -> "fail-unknown-gate"
]

VerificationTest[
    QuantumQASM["OPENQASM 3.0;\nqubit[1] q;\nfor int i in [0:2] { h q[0]; }"]["Category"],
    "Unsupported",
    TestID -> "fail-unsupported-construct"
]

VerificationTest[
    QuantumQASM["OPENQASM 3.0;\nqubit[1] q;\nh q[5];"]["Category"],
    "Range",
    TestID -> "fail-out-of-range"
]

VerificationTest[
    QuantumQASM["OPENQASM 3.0;\nqubit[1] q;\nrx(Run[\"x\"]) q[0];"]["Category"],
    "BadExpression",
    TestID -> "fail-injection-refused"
]

VerificationTest[
    FailureQ @ QuantumQASM["OPENQASM 3.0;\nqubit[1] q;\ncx q[0];"],
    True,
    TestID -> "fail-arity"
]

(* ---- named-circuit guard regression: the QASM overload must not shadow named circuits ---- *)

VerificationTest[
    Head @ QuantumCircuitOperator["Bell"],
    QuantumCircuitOperator,
    TestID -> "named-circuit-still-works"
]

VerificationTest[
    FailureQ @ Quiet @ QuantumCircuitOperator["DefinitelyNotACircuit"],
    True,
    TestID -> "invalid-name-still-fails"
]

(* ---- constructor surface ---- *)

VerificationTest[
    exactQ[QuantumCircuitOperator["OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[2] q;\nh q[0];\ncx q[0], q[1];"], bell],
    True,
    TestID -> "constructor-from-string"
]

VerificationTest[
    With[{tmp = FileNameJoin[{$TemporaryDirectory, "qf-openqasm-test.qasm"}]},
        Export[tmp, "OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[2] q;\nh q[0];\ncx q[0], q[1];", "Text"];
        exactQ[QuantumCircuitOperator[File[tmp]], bell]
    ],
    True,
    TestID -> "constructor-from-file"
]

(* ---- hub delegation: properties go through QuantumQASM ---- *)

VerificationTest[
    bell["QASM"] === QuantumQASM[bell, "WL"],
    True,
    TestID -> "circuit-qasm-delegates-to-hub"
]

VerificationTest[
    QuantumOperator["H", {1}]["QASM"] === QuantumQASM[QuantumOperator["H", {1}]],
    True,
    TestID -> "operator-qasm-delegates-to-hub"
]

VerificationTest[
    QuantumOperator["H", {1}]["SimpleQASM"] === QuantumQASM[QuantumOperator["H", {1}], "Simple"],
    True,
    TestID -> "operator-simpleqasm-delegates-to-hub"
]

(* the native "WL" export mode is dependency-free (no qiskit) and starts with the header *)
VerificationTest[
    StringStartsQ[QuantumQASM[bell, "WL"], "OPENQASM 3.0"],
    True,
    TestID -> "wl-export-dependency-free"
]

(* ---- back-compat: deprecated ImportQASMCircuit still imports (now returns a QCO) ---- *)

VerificationTest[
    exactQ[
        Quiet[ImportQASMCircuit["OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[2] q;\nh q[0];\ncx q[0], q[1];"], ImportQASMCircuit::deprecated],
        bell
    ],
    True,
    TestID -> "import-qasm-circuit-deprecated-shim"
]

EndTestSection[]
