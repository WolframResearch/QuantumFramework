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


shouldBeUnitary[circuit_] := MatchQ[QuantumOperator @ circuit, _QuantumOperator] && Chop @ QuantumOperator[circuit]["UnitaryQ"]
SetAttributes[shouldBeUnitary, HoldFirst]

VerificationTest[shouldBeUnitary[QuantumCircuitOperator["GrayOracle"[# &]]], True, TestID -> "Unitary-GrayOracle-default"]
VerificationTest[shouldBeUnitary[QuantumCircuitOperator["GrayOracle"[# &, 3, Automatic, "RX"]]], True, TestID -> "Unitary-GrayOracle-RX"]
VerificationTest[shouldBeUnitary[QuantumCircuitOperator["GrayOracle"[# &, 3, Automatic, "RY"]]], True, TestID -> "Unitary-GrayOracle-RY"]
VerificationTest[shouldBeUnitary[QuantumCircuitOperator["GrayOracle"[# &, 3, Automatic, "RZ"]]], True, TestID -> "Unitary-GrayOracle-RZ"]
VerificationTest[shouldBeUnitary[QuantumCircuitOperator["BooleanOracleR"[BooleanFunction[1, 3]]]], True, TestID -> "Unitary-BooleanOracleR"]
VerificationTest[shouldBeUnitary[QuantumCircuitOperator["BooleanOracle"[BooleanFunction[1, 3]]]], True, TestID -> "Unitary-BooleanOracle"]
VerificationTest[shouldBeUnitary[QuantumCircuitOperator["PhaseOracle"[BooleanFunction[1, 3]]]], True, TestID -> "Unitary-PhaseOracle"]
VerificationTest[shouldBeUnitary[QuantumCircuitOperator["BernsteinVaziraniOracle"[]]], True, TestID -> "Unitary-BernsteinVaziraniOracle"]
VerificationTest[shouldBeUnitary[QuantumCircuitOperator["GroverDiffusion"[{1, 2, 3}]]], True, TestID -> "Unitary-GroverDiffusion"]
VerificationTest[shouldBeUnitary[QuantumCircuitOperator["GroverPhaseDiffusion"[{1, 2, 3}]]], True, TestID -> "Unitary-GroverPhaseDiffusion"]
VerificationTest[shouldBeUnitary[QuantumCircuitOperator["Toffoli"]], True, TestID -> "Unitary-Toffoli"]
VerificationTest[shouldBeUnitary[QuantumCircuitOperator["Fredkin"]], True, TestID -> "Unitary-Fredkin"]
VerificationTest[shouldBeUnitary[QuantumCircuitOperator["Magic"[]]], True, TestID -> "Unitary-Magic"]

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


BeginTestSection["QuantumCircuitOperator - state preparation round-trip"]

(* QuantumCircuitOperator["QuantumState"[qs]][] must reproduce qs (the |0...0> -> qs
   map). Guards the multiplexer-dagger angle negation: a state with any nonzero
   RY/RZ angle is wrong if the dagger is not applied. Exact-equality round-trips
   first (including global phase), then fidelity for numeric / random states. *)

(* exact-equality round-trips: Norm of the state-vector difference is exactly 0 *)
VerificationTest[
    Chop @ Norm[QuantumCircuitOperator["QuantumState"[QuantumState["1"]]][]["StateVector"] - QuantumState["1"]["StateVector"]],
    0,
    TestID -> "StatePreparation-roundtrip-ket1"
]

VerificationTest[
    Chop @ Norm[QuantumCircuitOperator["QuantumState"[QuantumState["010"]]][]["StateVector"] - QuantumState["010"]["StateVector"]],
    0,
    TestID -> "StatePreparation-roundtrip-ket010"
]

VerificationTest[
    Chop @ Norm[QuantumCircuitOperator["QuantumState"[QuantumState["GHZ"[3]]]][]["StateVector"] - QuantumState["GHZ"[3]]["StateVector"]],
    0,
    TestID -> "StatePreparation-roundtrip-GHZ3"
]

VerificationTest[
    Chop @ Norm[QuantumCircuitOperator["QuantumState"[QuantumState["W"[3]]]][]["StateVector"] - QuantumState["W"[3]]["StateVector"]],
    0,
    TestID -> "StatePreparation-roundtrip-W3"
]

(* global phase is preserved, not just the ray *)
VerificationTest[
    With[{qs = Exp[I 6/5] QuantumState["+"]},
        Chop @ Norm[N @ QuantumCircuitOperator["QuantumState"[qs]][]["StateVector"] - N @ qs["StateVector"]]
    ],
    0,
    TestID -> "StatePreparation-roundtrip-global-phase"
]

(* fidelity round-trips for numeric / random states across dimensions *)
VerificationTest[
    Block[{qs = QuantumState[Normalize[{1, 2 I, -3, 0, 0.5, -1.5 I, 2, 1 + I}]]},
        Chop[1 - QuantumDistance[QuantumCircuitOperator["QuantumState"[qs]][], qs, "Fidelity"]]
    ],
    1,
    SameTest -> (Abs[#1 - #2] < 10^-8 &),
    TestID -> "StatePreparation-roundtrip-explicit-vector"
]

VerificationTest[
    Table[
        With[{qs = (SeedRandom[n]; QuantumState["RandomPure"[n]])},
            Chop[1 - QuantumDistance[QuantumCircuitOperator["QuantumState"[qs]][], qs, "Fidelity"]]
        ],
        {n, 1, 5}
    ],
    {1, 1, 1, 1, 1},
    SameTest -> (Max[Abs[#1 - #2]] < 10^-8 &),
    TestID -> "StatePreparation-roundtrip-random-1to5"
]

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - StatePreparation block-diagonal (heterogeneous qudits)"]

(* The block-diagonal algorithm (Method -> "BlockDiagonal") prepares an arbitrary pure
   state on a register whose qudits may have different dimensions. Workflow under test:
   generate a random state vector for a dimension profile, build the circuit, run qc[]
   (no input -> |0...0>), and compare with the original vector (norm, exact amplitudes,
   fidelity). prepCirc forces "BlockDiagonal" so every profile, including all-qubit ones,
   exercises this algorithm rather than the Automatic multiplexer route. *)

heteroDims = {{2, 3, 2}, {4, 2, 3}, {3, 3, 3}, {5, 3}, {2, 2, 2, 2}, {2, 5}, {7}};

svGen[dims_, seed_] := (SeedRandom[seed]; Normalize[RandomComplex[{-1 - I, 1 + I}, Times @@ dims]]);

prepCirc[dims_, seed_] := QuantumCircuitOperator["StatePreparation"[QuantumState[svGen[dims, seed], QuantumBasis[dims]], Method -> "BlockDiagonal"]];

(* qc[] preserves the norm *)
VerificationTest[
    Table[Chop[Norm[prepCirc[d, 11][]["StateVector"]] - 1], {d, heteroDims}],
    ConstantArray[0, Length[heteroDims]],
    SameTest -> (Max[Abs[#1 - #2]] < 10^-8 &),
    TestID -> "StatePreparation-blockdiag-norm"
]

(* qc[] reproduces the original amplitude vector exactly *)
VerificationTest[
    Table[Chop[Norm[prepCirc[d, 11][]["StateVector"] - svGen[d, 11]]], {d, heteroDims}],
    ConstantArray[0, Length[heteroDims]],
    SameTest -> (Max[Abs[#1 - #2]] < 10^-8 &),
    TestID -> "StatePreparation-blockdiag-vector-roundtrip"
]

(* fidelity with the target is 1 (QuantumDistance "Fidelity" returns a distance, 0 when identical) *)
VerificationTest[
    Table[Chop[1 - QuantumDistance[prepCirc[d, 11][], QuantumState[svGen[d, 11], QuantumBasis[d]], "Fidelity"]], {d, heteroDims}],
    ConstantArray[1, Length[heteroDims]],
    SameTest -> (Max[Abs[#1 - #2]] < 10^-8 &),
    TestID -> "StatePreparation-blockdiag-fidelity"
]

(* the prepared state keeps the requested per-qudit dimensions *)
VerificationTest[
    Table[prepCirc[d, 11][]["Dimensions"], {d, heteroDims}],
    heteroDims,
    TestID -> "StatePreparation-blockdiag-dims-preserved"
]

(* the compiled circuit operator is unitary *)
VerificationTest[
    Table[(QuantumOperator @ prepCirc[d, 11])["UnitaryQ"], {d, heteroDims}],
    ConstantArray[True, Length[heteroDims]],
    TestID -> "StatePreparation-blockdiag-unitary"
]

(* one block-diagonal gate per qudit (count the flattened operators) *)
VerificationTest[
    Table[Length[prepCirc[d, 11]["Flatten"]["Operators"]], {d, heteroDims}],
    Length /@ heteroDims,
    TestID -> "StatePreparation-blockdiag-gate-count"
]

(* explicit |0...0> input matches the qc[] form *)
VerificationTest[
    With[{d = {2, 3, 2}},
        Chop @ Norm[
            prepCirc[d, 11][QuantumState[SparseArray[{1 -> 1}, Times @@ d], QuantumBasis[d]]]["StateVector"]
            - prepCirc[d, 11][]["StateVector"]
        ]
    ],
    0,
    SameTest -> (Abs[#1 - #2] < 10^-8 &),
    TestID -> "StatePreparation-blockdiag-explicit-input-matches"
]

(* a state given in a non-computational basis is prepared via its Computational form *)
VerificationTest[
    With[{qs = QuantumState["+", QuantumBasis["X"]]},
        QuantumCircuitOperator["StatePreparation"[qs, Method -> "BlockDiagonal"]][] == qs["Computational"]
    ],
    True,
    TestID -> "StatePreparation-blockdiag-noncomputational-basis"
]

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - StatePreparation method routing"]

(* Automatic routes all-qubit states through the multiplexer and any qudit state through
   the block-diagonal algorithm, with no explicit Method. "QuantumState" and
   "StatePreparation" are the same entry point; "StatePrep" has been removed. *)

(* single qudit (dim != 2) round-trips under Automatic *)
VerificationTest[
    With[{qs = QuantumState[svGen[{5}, 3], QuantumBasis[{5}]]},
        Chop @ Norm[QuantumCircuitOperator["StatePreparation"[qs]][]["StateVector"] - svGen[{5}, 3]]
    ],
    0,
    SameTest -> (Abs[#1 - #2] < 10^-8 &),
    TestID -> "StatePreparation-automatic-single-qudit"
]

(* heterogeneous multi-qudit round-trips under Automatic (no Method) *)
VerificationTest[
    Table[
        With[{qs = QuantumState[svGen[d, 4], QuantumBasis[d]]},
            Chop @ Norm[QuantumCircuitOperator["StatePreparation"[qs]][]["StateVector"] - svGen[d, 4]]
        ],
        {d, {{3, 2}, {2, 3, 2}, {3, 3, 3}}}
    ],
    {0, 0, 0},
    SameTest -> (Max[Abs[#1 - #2]] < 10^-8 &),
    TestID -> "StatePreparation-automatic-heterogeneous"
]

(* "QuantumState" and "StatePreparation" name the same entry point *)
VerificationTest[
    With[{qs = QuantumState["GHZ"[3]]},
        InputForm[QuantumCircuitOperator["QuantumState"[qs]]] === InputForm[QuantumCircuitOperator["StatePreparation"[qs]]]
    ],
    True,
    TestID -> "StatePreparation-QuantumState-alias"
]

(* Automatic on an all-qubit state uses the multiplexer: more gates than block-diagonal *)
VerificationTest[
    With[{qs = QuantumState[svGen[{2, 2, 2}, 9], QuantumBasis[{2, 2, 2}]]},
        Length[QuantumCircuitOperator["StatePreparation"[qs]]["Flatten"]["Operators"]] >
            Length[QuantumCircuitOperator["StatePreparation"[qs, Method -> "BlockDiagonal"]]["Flatten"]["Operators"]]
    ],
    True,
    TestID -> "StatePreparation-multiplexer-vs-blockdiagonal-gatecount"
]

(* "Multiplexer" forced on a qudit is qubit-only: returns a Failure with ::qubitonly *)
VerificationTest[
    FailureQ @ QuantumCircuitOperator["StatePreparation"[QuantumState[svGen[{3}, 1], QuantumBasis[{3}]], Method -> "Multiplexer"]],
    True,
    {QuantumCircuitOperator::qubitonly},
    TestID -> "StatePreparation-multiplexer-qudit-failure"
]

(* an unrecognized Method returns a Failure with ::badmethod *)
VerificationTest[
    FailureQ @ QuantumCircuitOperator["StatePreparation"[QuantumState["GHZ"[2]], Method -> "Nonsense"]],
    True,
    {QuantumCircuitOperator::badmethod},
    TestID -> "StatePreparation-badmethod-failure"
]

(* a mixed state returns a Failure and emits QuantumCircuitOperator::nonpure *)
VerificationTest[
    FailureQ @ QuantumCircuitOperator["StatePreparation"[QuantumState["RandomMixed", {2, 3}]]],
    True,
    {QuantumCircuitOperator::nonpure},
    TestID -> "StatePreparation-nonpure-failure"
]

(* the no-argument forms build a circuit (3-qubit uniform superposition) *)
VerificationTest[Head @ QuantumCircuitOperator["StatePreparation"[]], QuantumCircuitOperator, TestID -> "StatePreparation-noarg"]

VerificationTest[Head @ QuantumCircuitOperator["QuantumState"[]], QuantumCircuitOperator, TestID -> "QuantumState-noarg"]

(* "StatePrep" has been removed: it no longer matches a constructor rule *)
VerificationTest[
    FailureQ @ QuantumCircuitOperator["StatePrep"[QuantumState["GHZ"[3]]]],
    True,
    {QuantumCircuitOperator::invalidName},
    TestID -> "StatePrep-removed"
]

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - Graph constructor"]

(* gh#2: QuantumCircuitOperator[g_Graph] used to call the now-removed
   GraphTensorNetwork and return an inert TensorNetworkQuantumCircuit[...].
   It must evaluate to a real QuantumCircuitOperator after the rename to
   ToTensorNetworkGraph. *)
VerificationTest[
    Head @ QuantumCircuitOperator[CycleGraph[4]],
    QuantumCircuitOperator,
    TestID -> "Graph-CycleGraph-head"
]

VerificationTest[
    QuantumCircuitOperator[CycleGraph[4]]["GateCount"],
    4,
    TestID -> "Graph-CycleGraph-gateCount"
]

VerificationTest[
    Head @ QuantumCircuitOperator[PathGraph[Range[3]]],
    QuantumCircuitOperator,
    TestID -> "Graph-PathGraph-head"
]

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - contraction path"]

(* gh#4: TensorNetworkCompile now exposes the contraction-path argument
   through the "Path" option, reachable via
   Method -> {"TensorNetwork", "Path" -> spec}. Greedy and Optimal paths
   must produce results numerically equal to the Automatic default. *)
With[{
    qc = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "H" -> 2}],
    psi0 = QuantumState["00"]
},
    VerificationTest[
        Chop @ Normal @ qc[psi0, Method -> {"TensorNetwork", "Path" -> "Greedy"}]["StateVector"],
        Chop @ Normal @ qc[psi0]["StateVector"],
        TestID -> "ContractionPath-Greedy-eq-default"
    ];

    VerificationTest[
        Chop @ Normal @ qc[psi0, Method -> {"TensorNetwork", "Path" -> "Optimal"}]["StateVector"],
        Chop @ Normal @ qc[psi0]["StateVector"],
        TestID -> "ContractionPath-Optimal-eq-default"
    ];

    VerificationTest[
        Head @ qc[psi0, Method -> {"TensorNetwork", "Path" -> "Greedy"}],
        QuantumState,
        TestID -> "ContractionPath-Greedy-returns-state"
    ]
]

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - Schrodinger method"]

(* Method -> "Schrodinger" folds gates one at a time; each per-gate apply
   re-enters quantumCircuitApply, so its options must travel inside the
   Method spec. VerificationTest fails on any emitted message, so these
   also guard against OptionValue::nodef from a stray top-level option. *)
With[{
    qc = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}],
    psi0 = QuantumState["Register"[2]]
},
    VerificationTest[
        Normal @ qc[psi0, Method -> "Schrodinger"]["StateVector"],
        {1/Sqrt[2], 0, 0, 1/Sqrt[2]},
        TestID -> "Schrodinger-Bell-stateVector"
    ];

    VerificationTest[
        Head @ qc[psi0, Method -> "Schrodinger"],
        QuantumState,
        TestID -> "Schrodinger-returns-state"
    ];

    VerificationTest[
        Normal @ qc[psi0, Method -> "Schrodinger"]["StateVector"],
        Normal @ qc[psi0]["StateVector"],
        TestID -> "Schrodinger-eq-TensorNetwork-default"
    ]
]

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - operator placement"]

(* An operator rule `op -> order` must mean the same whether written bare or
   inside the operator list. For a permutation (a swap, whose output and input
   orders differ) the two forms used to diverge: the bare rule reordered the
   circuit correctly while the list form collapsed the swap to the identity. *)

VerificationTest[
    QuantumCircuitOperator[QuantumOperator["Permutation" -> {2, 1}] -> {2, 3}] ==
        QuantumCircuitOperator[{QuantumOperator["Permutation" -> {2, 1}] -> {2, 3}}],
    True,
    TestID -> "Placement-rule-vs-list-permutation"
]

VerificationTest[
    QuantumCircuitOperator[QuantumOperator["CNOT"] -> {2, 3}] ==
        QuantumCircuitOperator[{QuantumOperator["CNOT"] -> {2, 3}}],
    True,
    TestID -> "Placement-rule-vs-list-CNOT"
]

(* The placed permutation is a genuine SWAP on the target qudits, not identity. *)
VerificationTest[
    Normal @ QuantumCircuitOperator[{QuantumOperator["Permutation" -> {2, 1}] -> {2, 3}}][
            "CircuitOperator"]["OrderedMatrixRepresentation"],
    {{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}},
    TestID -> "Placement-permutation-is-SWAP"
]

(* Acting on |0,1,0>, a swap on qudits {2, 3} yields |0,0,1> for both forms. *)
With[{
    in = QuantumTensorProduct[QuantumState["0"], QuantumState["1"], QuantumState["0"]],
    out = QuantumTensorProduct[QuantumState["0"], QuantumState["0"], QuantumState["1"]]
},
    VerificationTest[
        QuantumCircuitOperator[QuantumOperator["Permutation" -> {2, 1}] -> {2, 3}][in] == out,
        True,
        TestID -> "Placement-swap-action-rule"
    ];
    VerificationTest[
        QuantumCircuitOperator[{QuantumOperator["Permutation" -> {2, 1}] -> {2, 3}}][in] == out,
        True,
        TestID -> "Placement-swap-action-list"
    ]
]

EndTestSection[]
