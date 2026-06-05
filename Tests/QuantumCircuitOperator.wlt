BeginTestSection["QuantumCircuitOperator - constructors"]

VerificationTest[QuantumCircuitOperator["GHZ"[3]]["Arity"], 3, TestID -> "GHZ-3"]

VerificationTest[QuantumCircuitOperator["Bell"]["Arity"], 2, TestID -> "Bell-bare"]

VerificationTest[QuantumCircuitOperator["Toffoli"]["Arity"], 3, TestID -> "Toffoli-bare"]

VerificationTest[QuantumCircuitOperator["Fredkin"]["Arity"], 3, TestID -> "Fredkin-bare"]

VerificationTest[QuantumCircuitOperator["Fourier"[3]]["Arity"], 3, TestID -> "Fourier-3"]

VerificationTest[QuantumCircuitOperator["InverseFourier"[3]]["Arity"], 3, TestID -> "InverseFourier-3"]

(* sequence-of-gates list (semantic composition - kept as List) *)
VerificationTest[
    QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]["Arity"],
    2,
    TestID -> "Sequence-of-gates-list"
]

(* Pauli-string circuit shortcut *)
VerificationTest[QuantumCircuitOperator["XYZ"]["Arity"], 3, TestID -> "PauliString-XYZ"]

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - basic action"]

(* Bell circuit takes |00> -> Bell pair *)
VerificationTest[
    Chop @ Norm[QuantumCircuitOperator["Bell"][QuantumState["00"]]["StateVector"] - QuantumState["Bell"]["StateVector"]],
    0,
    TestID -> "Bell-circuit-on-00"
]

VerificationTest[Length @ QuantumCircuitOperator["GHZ"[3]]["Elements"], 3, TestID -> "GHZ-gate-count"]

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - named circuit fidelity"]

(* These cover the class of v2.0 list-form regressions where a named circuit
   silently built a wrong operator: the Multiplexer collapsed every input list
   to arity 1, the Fourier circuit used a legacy {"PhaseShift", k} arg that
   v2.0 interprets as a tensor product, and the Graph circuit defaulted to
   the legacy {"C", "1"} controlled-gate spec. The fixes are observable as
   non-trivial arity and as matrix equivalence with the corresponding
   QuantumOperator. *)

(* Fourier circuit must match QuantumOperator["Fourier"[d]] *)
fourierMatDiff[n_] := Chop @ Norm @ Flatten[
    N @ (QuantumOperator @ QuantumCircuitOperator["Fourier"[n]])["Sort"]["Matrix"]
    - N @ QuantumOperator["Fourier"[2^n]]["Matrix"]
]

VerificationTest[fourierMatDiff[2], 0, TestID -> "Fourier-2-matches-QO"]
VerificationTest[fourierMatDiff[3], 0, TestID -> "Fourier-3-matches-QO"]
VerificationTest[fourierMatDiff[4], 0, TestID -> "Fourier-4-matches-QO"]

(* Multiplexer arity scales with the number of ops *)
VerificationTest[QuantumCircuitOperator["Multiplexer"["X", "Y"]]["Arity"], 2, TestID -> "Multiplexer-2-arity"]
VerificationTest[QuantumCircuitOperator["Multiplexer"["X", "Y", "Z"]]["Arity"], 3, TestID -> "Multiplexer-3-arity"]
VerificationTest[QuantumCircuitOperator["Multiplexer"["X", "Y", "Z", "H"]]["Arity"], 3, TestID -> "Multiplexer-4-arity"]
VerificationTest[QuantumCircuitOperator["Multiplexer"["X", "Y", "Z", "H", "S"]]["Arity"], 4, TestID -> "Multiplexer-5-arity"]
VerificationTest[QuantumCircuitOperator["Multiplexer"["I", "X"]]["Arity"], 2, TestID -> "Multiplexer-with-I"]

(* Multiplexer must produce a unitary *)
VerificationTest[
    (QuantumOperator @ QuantumCircuitOperator["Multiplexer"["X", "Y", "Z", "H"]])["UnitaryQ"],
    True,
    TestID -> "Multiplexer-unitary"
]

(* Graph circuit arity tracks vertex count *)
VerificationTest[QuantumCircuitOperator["Graph"[CompleteGraph[5]]]["Arity"], 5, TestID -> "Graph-K5-arity"]
VerificationTest[QuantumCircuitOperator["Graph"[PathGraph[Range[4]]]]["Arity"], 4, TestID -> "Graph-P4-arity"]
VerificationTest[QuantumCircuitOperator["Graph"[CycleGraph[3]]]["Arity"], 3, TestID -> "Graph-C3-arity"]
VerificationTest[QuantumCircuitOperator["Graph"[CompleteGraph[3], 2]]["Arity"], 5, TestID -> "Graph-K3-offset-2"]

(* Graph element count: 1 H per vertex + 1 CNOT per edge.
   K5 has 5 vertices, 10 edges -> 15 elements. *)
VerificationTest[Length @ QuantumCircuitOperator["Graph"[CompleteGraph[5]]]["Elements"], 15, TestID -> "Graph-K5-element-count"]

(* Constructing every named circuit emits no QF-level messages (no failprop,
   undefprop, or InvalidName from the construction itself). *)
quietBuild[expr_] := Module[{r, msgs = {}},
    Internal`HandlerBlock[
        {"Message", Function[m,
            With[{s = ToString[InputForm[m]]},
                If[ StringContainsQ[s, "failprop" | "undefprop" | "InvalidName"],
                    AppendTo[msgs, m]
                ]
            ]
        ]},
        r = expr
    ];
    {Head[r], Length[msgs]}
]

VerificationTest[quietBuild[QuantumCircuitOperator["Graph"[CompleteGraph[5]]]], {QuantumCircuitOperator, 0}, TestID -> "Graph-K5-quiet-build"]
VerificationTest[quietBuild[QuantumCircuitOperator["Fourier"[3]]], {QuantumCircuitOperator, 0}, TestID -> "Fourier-3-quiet-build"]
VerificationTest[quietBuild[QuantumCircuitOperator["Multiplexer"["X", "Y", "Z", "H", "S"]]], {QuantumCircuitOperator, 0}, TestID -> "Multiplexer-quiet-build"]

(* The Grover formula-default rule used to feed Nothing as the third positional
   arg of "PhaseOracle"[formula, varSpec, Nothing], which doesn't vanish in
   function-call arg position and missed every PhaseOracle dispatch pattern -
   the call fell back to the Grover rule itself and recursed to the depth
   limit. With the recursion fixed and an invalidArgs catch-all in place, each
   of these constructions must return a QuantumCircuitOperator. *)

VerificationTest[Head @ QuantumCircuitOperator["Grover"[]], QuantumCircuitOperator, TestID -> "Grover-bare"]
VerificationTest[Head @ QuantumCircuitOperator["GroverPhase"[]], QuantumCircuitOperator, TestID -> "GroverPhase-bare"]
VerificationTest[Head @ QuantumCircuitOperator["Grover0"[]], QuantumCircuitOperator, TestID -> "Grover0-bare"]
VerificationTest[Head @ QuantumCircuitOperator["GroverPhase0"[]], QuantumCircuitOperator, TestID -> "GroverPhase0-bare"]

(* Grover with an explicit oracle / formula goes through the QuantumFrameworkOperatorQ
   dispatcher (line 121), whose RHS previously used a {nameString, order, gate}
   list interpreted post-v2.0 as a sequence-of-gates rather than the intended
   single-gate call. *)

VerificationTest[
    Head @ QuantumCircuitOperator["Grover"[QuantumOperator["X", {1}]]],
    QuantumCircuitOperator,
    TestID -> "Grover-with-op"
]

VerificationTest[
    Head @ QuantumCircuitOperator["GroverPhase"[BooleanFunction[2^6, 3]]],
    QuantumCircuitOperator,
    TestID -> "GroverPhase-with-formula"
]

(* DeutschJozsa[Phase]Oracle previously used the same {nameString, formula}
   list-form post-stripping "DeutschJozsa". The resulting circuit was dispatched
   as two separate gates and the downstream DeutschJozsa[] / Deutsch[] circuits
   inherited the broken width. *)

VerificationTest[
    Head @ QuantumCircuitOperator["DeutschJozsaPhaseOracle"[]],
    QuantumCircuitOperator,
    TestID -> "DeutschJozsaPhaseOracle-bare"
]

VerificationTest[
    Head @ QuantumCircuitOperator["DeutschJozsaBooleanOracle"[]],
    QuantumCircuitOperator,
    TestID -> "DeutschJozsaBooleanOracle-bare"
]

VerificationTest[
    Head @ QuantumCircuitOperator["DeutschJozsa"[]],
    QuantumCircuitOperator,
    TestID -> "DeutschJozsa-bare"
]

VerificationTest[
    Head @ QuantumCircuitOperator["DeutschPhase"[]],
    QuantumCircuitOperator,
    TestID -> "DeutschPhase-bare"
]

VerificationTest[
    Head @ QuantumCircuitOperator["Deutsch"[]],
    QuantumCircuitOperator,
    TestID -> "Deutsch-bare"
]

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - shortcut roundtrip"]

(* QuantumShortcut compresses a circuit to its named-shorthand form;
   re-feeding it through the QuantumCircuitOperator constructor reconstructs
   an equivalent circuit. *)

VerificationTest[
    QuantumCircuitOperator @ QuantumShortcut[QuantumCircuitOperator["Bell"]] == QuantumCircuitOperator["Bell"],
    True,
    TestID -> "Roundtrip-Bell"
]

VerificationTest[
    QuantumCircuitOperator @ QuantumShortcut[QuantumCircuitOperator["GHZ"[3]]] == QuantumCircuitOperator["GHZ"[3]],
    True,
    TestID -> "Roundtrip-GHZ-3"
]

VerificationTest[
    QuantumCircuitOperator @ QuantumShortcut[QuantumCircuitOperator["Toffoli"]] == QuantumCircuitOperator["Toffoli"],
    True,
    TestID -> "Roundtrip-Toffoli"
]

VerificationTest[
    QuantumCircuitOperator @ QuantumShortcut[QuantumCircuitOperator["Fredkin"]] == QuantumCircuitOperator["Fredkin"],
    True,
    TestID -> "Roundtrip-Fredkin"
]

(* hand-built circuit *)
VerificationTest[
    With[{qc = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "Phase"[Pi/4] -> 2}]},
        QuantumCircuitOperator @ QuantumShortcut[qc] == qc
    ],
    True,
    TestID -> "Roundtrip-handbuilt"
]

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - Qiskit export"]

(* Qiskit export must hand Python a tuple (name, args), not a WL function call.
   The bug surfaced as `'Measure'[(1,)]` on the Python side; lock it down. *)

VerificationTest[
    Quiet @ Head @ QuantumCircuitOperator[{{1}}]["Qiskit"],
    QiskitCircuit,
    TestID -> "Qiskit-measurement-from-list"
]

VerificationTest[
    Quiet @ Head @ QuantumCircuitOperator["Bell"]["Qiskit"],
    QiskitCircuit,
    TestID -> "Qiskit-Bell"
]

VerificationTest[
    Quiet @ Head @ QuantumCircuitOperator[{"U"[Pi/3, Pi/4, Pi/5] -> 1}]["Qiskit"],
    QiskitCircuit,
    TestID -> "Qiskit-U3"
]

(* Last-resort invalidArgs catch-all: a known name with unmatched call shape
   returns Failure["InvalidArguments"] rather than recursing or falling
   through unevaluated. This is the safety net that would have turned the
   GroverPhase[] recursion-limit hang into a clear failure. *)

VerificationTest[
    Quiet @ QuantumCircuitOperator["Fourier"[-1]],
    Failure["InvalidArguments", _],
    SameTest -> MatchQ,
    TestID -> "InvalidArgs-Fourier"
]

VerificationTest[
    Quiet @ QuantumCircuitOperator["GHZ"["not-an-integer"]],
    Failure["InvalidArguments", _],
    SameTest -> MatchQ,
    TestID -> "InvalidArgs-GHZ"
]

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - failure"]

VerificationTest[
    Quiet @ QuantumCircuitOperator["NotACircuit"[2]],
    Failure["InvalidName", _],
    SameTest -> MatchQ,
    TestID -> "InvalidName-call-form"
]

VerificationTest[
    Quiet @ QuantumCircuitOperator["NotACircuit"],
    Failure["InvalidName", _],
    SameTest -> MatchQ,
    TestID -> "InvalidName-bare"
]

EndTestSection[]
