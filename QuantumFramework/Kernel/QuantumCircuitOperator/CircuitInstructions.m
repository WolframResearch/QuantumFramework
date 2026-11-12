Package["Wolfram`QuantumFramework`"]

PackageScope["circuitInstructions"]
PackageScope["instructionFailure"]



(* The one lowering of a circuit for export. Both exporters (the qiskit circuit builder and the
   native OpenQASM writer) read this list and only print it: neither reorders wires, changes a
   measurement basis or decides what a reset means on its own.

   Each instruction is one of
     {"Gate", qo}            a unitary operator on the wires qo["InputOrder"] (equal to its output order)
     {"Reset", w}            discard qubit w and leave |0> there
     {"Measure", w, k}       measure qubit w in the computational basis into classical bit k (from 0)
     {"Barrier", wires}      a barrier ({} spans the whole circuit)

   What the other circuit elements become:
     a measurement in the basis B on wires T: the unitary taking B's outcome vectors to the
       computational basis, then a computational measurement of each wire of T. QF numbers the
       outcomes by the measured wires in increasing order, and the rotation is written in that
       order, so classical bits read back by wire give QF's outcome. The inverse rotation after the
       measurement leaves the qubits in the outcome's basis vector, as QF leaves them.
     "Reset"[state] on wires w, and a state injected on wires w (an operator with no input):
       a reset of each wire of w, then the gates preparing state from |0...0>.
     a permutation of wires: the SWAP gates that compose to it.
   Anything with no such sequence (a channel, a gate on qudits that are not qubits, a non-unitary
   or wire-moving operator, a mixed or symbolic reset state, a measurement that is symbolic, of
   qudits, or with degenerate or non-projective outcomes) is a Failure naming it. *)

instructionFailure[template_String, params_List : {}] :=
    Failure["QuantumCircuitInstructions", <|"MessageTemplate" -> template, "MessageParameters" -> params|>]

circuitInstructions[qco_QuantumCircuitOperator] := Enclose @ With[{
    raw = Catenate[Confirm @ elementInstructions[#] & /@ qco["Flatten"]["Elements"]]
},
    (* classical bits are numbered in the order the measurements happen *)
    With[{positions = Flatten @ Position[raw, {"Measure", _}, {1}, Heads -> False]},
        ReplacePart[raw, MapIndexed[#1 -> Append[raw[[#1]], First[#2] - 1] &, positions]]
    ]
]

elementInstructions[barrier_ ? BarrierQ] := {{"Barrier", Cases[barrier, _Integer ? Positive, Infinity]}}

elementInstructions[qmo_QuantumMeasurementOperator] := measurementInstructions[qmo]

elementInstructions[qc_QuantumChannel] :=
    instructionFailure["The channel `1` has no gate sequence; only unitary gates, resets and measurements do.", {qc["Label"]}]

elementInstructions[qo_QuantumOperator] := Which[
    MatchQ[qo["Label"], "Reset"[_]] && qo["OutputOrder"] === qo["InputOrder"],
        stateInjectionInstructions[resetState[qo], qo["Sort"]["OutputOrder"], qo["Label"]],
    qo["InputDimension"] == 1 && qo["OutputDimension"] >= 2,
        stateInjectionInstructions[qo["Sort"]["State"], qo["Sort"]["OutputOrder"], qo["Label"]],
    qo["Dimensions"] === {},
        {{"Gate", qo}},
    ! AllTrue[qo["Dimensions"], # == 2 &],
        instructionFailure["The operator `1` acts on qudits that are not qubits, which have no gate sequence.", {qo["Label"]}],
    qo["OutputOrder"] =!= qo["InputOrder"] && Sort[qo["OutputOrder"]] === Sort[qo["InputOrder"]] && qo["Sort"]["OutputOrder"] === qo["Sort"]["InputOrder"],
        elementInstructions[qo["Sort"]],
    qo["OutputOrder"] =!= qo["InputOrder"],
        instructionFailure["The operator `1` maps wires `2` to wires `3`; only operators acting on the wires they read have a gate sequence.", {qo["Label"], qo["InputOrder"], qo["OutputOrder"]}],
    MatchQ[qo["Label"], "\[Pi]"[__]],
        permutationInstructions[qo],
    MatrixQ[qo["Matrix"], NumericQ] && ! UnitaryMatrixQ[N @ Normal @ qo["Matrix"], Tolerance -> 10 ^ -8],
        instructionFailure["The operator `1` is not unitary, so it has no gate sequence.", {qo["Label"]}],
    True,
        {{"Gate", qo}}
]

elementInstructions[el_] := instructionFailure["The circuit element `1` has no gate sequence.", {el}]


(* -- a pure qubit state left on wires (a reset to it, or its injection) -- *)

zeroStateQ[state_] := With[{v = N @ Normal[state["StateVector"]]},
    TrueQ[Chop[Norm[v - UnitVector[Length[v], 1]]] == 0]
]

stateInjectionInstructions[state_, wires_List, label_] := Enclose @ Which[
    ! AllTrue[state["Dimensions"], # == 2 &],
        instructionFailure["`1` leaves a state that is not a qubit state, which has no gate sequence.", {label}],
    ! state["NumericQ"],
        instructionFailure["`1` leaves a symbolic state, which has no gate sequence; give the state numeric amplitudes.", {label}],
    ! state["PureStateQ"],
        instructionFailure["`1` leaves a mixed state, which has no gate sequence.", {label}],
    True,
        Join[
            {"Reset", #} & /@ wires,
            If[ zeroStateQ[state],
                {},
                (* the wires are now |0...0>: the preparation circuit's |0> kets add nothing *)
                Catenate @ Map[
                    If[ #["InputDimension"] == 1 && #["OutputDimension"] >= 2 && zeroStateQ[#["Sort"]["State"]], {}, Confirm @ elementInstructions[#]] &,
                    (* the state keeps its label: a named state has a named preparation *)
                    QuantumCircuitOperator[QuantumCircuitOperator[state], wires]["Flatten"]["Operators"]
                ]
            ]
        ]
]


(* -- a permutation of wires: SWAPs, checked against the operator's own matrix -- *)

(* the qubit permutation sigma with M |...1 at wire j...> = |...1 at wire sigma[j]...> *)
permutationOfMatrix[m_, k_] := Table[
    With[{col = m[[All, 2 ^ (k - j) + 1]]},
        With[{pos = FirstPosition[Chop[N[Abs[col]]], 1., {0}, {1}][[1]]},
            If[pos == 0, 0, k - Floor[Log2[pos - 1]]]
        ]
    ],
    {j, k}
]

(* the matrix of SWAP on positions a, b of k qubits *)
swapMatrix[k_, {a_, b_}] := With[{idx = Range[0, 2 ^ k - 1]},
    With[{img = FromDigits[ReplacePart[#, {a -> #[[b]], b -> #[[a]]}], 2] & /@ IntegerDigits[idx, 2, k]},
        Normal @ SparseArray[Thread[Transpose[{img + 1, idx + 1}] -> 1], {2 ^ k, 2 ^ k}]
    ]
]

permutationInstructions[qo_] := Enclose @ With[{
    sorted = qo["Sort"], k = Length[qo["InputOrder"]]
},
    With[{
        m = N @ Normal @ sorted["Matrix"],
        wires = sorted["InputOrder"],
        sigma = permutationOfMatrix[N @ Normal @ sorted["Matrix"], k]
    },
        ConfirmAssert[Sort[sigma] === Range[k]];
        (* a cycle (c1 c2 ... cn) of sigma (c1 -> c2 -> ... -> cn -> c1): SWAP(c1, c2), ...,
           SWAP(c1, cn) in time order; the product is checked against the matrix, so a wrong
           convention cannot pass silently *)
        With[{swaps = Catenate[Function[c, {First[c], #} & /@ Rest[c]] /@ First[PermutationCycles[sigma]]]},
            ConfirmAssert[
                Chop[Norm[Flatten[Fold[swapMatrix[k, #2] . #1 &, IdentityMatrix[2 ^ k], swaps] - m]]] == 0
            ];
            {"Gate", QuantumOperator["SWAP", wires[[#]]]} & /@ swaps
        ]
    ]
]


(* -- measurements -- *)

(* the unit vector v with P == v v^dagger, for a rank-1 projector P *)
projectorVector[p_] := With[{j = First @ Ordering[N @ Re @ Diagonal[p], -1]},
    p[[All, j]] / Sqrt[p[[j, j]]]
]

measurementInstructions[qmo_] := Enclose @ With[{
    targets = qmo["Target"], k = Length[qmo["Target"]]
},
    ConfirmAssert[AllTrue[qmo["TargetDimensions"], # == 2 &],
        instructionFailure["A measurement of qudits that are not qubits has no gate sequence."]];
    With[{povm = N @ Normal[#] & /@ qmo["POVMElements"]},
        ConfirmAssert[AllTrue[povm, MatrixQ[#, NumericQ] &],
            instructionFailure["The measurement `1` has symbolic outcome operators, which have no gate sequence.", {qmo["Label"]}]];
        ConfirmAssert[Length[povm] === 2 ^ k && AllTrue[povm, MatrixRank[#, Tolerance -> 10 ^ -8] == 1 &],
            instructionFailure[
                "The measurement `1` does not have one outcome per basis state of its `2` qubits (its outcomes are degenerate or it is not projective), so a computational measurement after a change of basis cannot reproduce it.",
                {qmo["Label"], k}
            ]];
        (* QF writes each outcome operator over the measured wires in increasing order and numbers
           outcome i by the bits of i over those wires. The rotation takes outcome vector i to that
           computational state; its inverse after the measurement leaves the qubits in outcome
           vector i, as QF does. *)
        With[{rotation = Conjugate[projectorVector /@ povm]},
            ConfirmAssert[UnitaryMatrixQ[rotation, Tolerance -> 10 ^ -8],
                instructionFailure["The outcome vectors of the measurement `1` are not orthonormal, so it is not a basis measurement.", {qmo["Label"]}]];
            If[ AllTrue[povm, DiagonalMatrixQ[Chop[#]] &],
                (* a measurement in the computational basis: which outcome QF reports for each
                   computational state depends on the kind of measurement (an observable or "M"
                   reports the qubits' values, a basis lists its own order), so it is read off QF's
                   own evaluation in a circuit, and the bits are relabeled to match *)
                computationalMeasurementInstructions[qmo, Sort[targets], targets],
                Join[
                    {"Gate", #} & /@ basisChangeOperators[rotation, Sort[targets]],
                    {"Measure", #} & /@ targets,
                    {"Gate", #} & /@ basisChangeOperators[ConjugateTranspose[rotation], Sort[targets]]
                ]
            ]
        ]
    ]
]

(* outcome[[b + 1]]: the outcome QF reports in a circuit when the measured qubits (increasing
   wires) are in computational state b *)
computationalOutcomes[qmo_, wires_] := With[{k = Length[wires]},
    Table[
        First @ Ordering[
            N @ QuantumCircuitOperator[{QuantumOperator[QuantumState[UnitVector[2 ^ k, b + 1], ConstantArray[2, k]], {wires, {}}], qmo}][]["ProbabilitiesList"],
            -1
        ] - 1,
        {b, 0, 2 ^ k - 1}
    ]
]

computationalMeasurementInstructions[qmo_, wires_, targets_] := Enclose @ With[{
    outcome = computationalOutcomes[qmo, wires], k = Length[wires]
},
    ConfirmAssert[Sort[outcome] === Range[0, 2 ^ k - 1]];
    If[ outcome === Range[0, 2 ^ k - 1],
        {"Measure", #} & /@ targets,
        (* the permutation taking computational state b to the state of its reported outcome *)
        With[{perm = Normal @ SparseArray[Thread[Transpose[{outcome + 1, Range[2 ^ k]}] -> 1], {2 ^ k, 2 ^ k}]},
            Join[
                {"Gate", #} & /@ basisChangeOperators[perm, wires],
                {"Measure", #} & /@ targets,
                {"Gate", #} & /@ basisChangeOperators[Transpose[perm], wires]
            ]
        ]
    ]
]

(* the basis change as single-qubit operators when it is a tensor product of them, otherwise one
   operator on all the wires *)
basisChangeOperators[u_, wires_] := Replace[kroneckerFactors[u, Length[wires]], {
    factors_List :> MapThread[QuantumOperator[#1, {#2}] &, {factors, wires}],
    _ :> {QuantumOperator[u, wires]}
}]

(* {A1, ..., Ak} with u == A1 (x) ... (x) Ak (2x2 factors), or Missing; the split of the first
   qubit is an operator-Schmidt decomposition, which has rank 1 exactly when it is a product *)
kroneckerFactors[u_, 1] := {u}
kroneckerFactors[u_, k_] := With[{
    r = Flatten[Transpose[ArrayReshape[u, {2, 2 ^ (k - 1), 2, 2 ^ (k - 1)}], {1, 3, 2, 4}], {{1, 2}, {3, 4}}]
},
    With[{svd = SingularValueDecomposition[r]},
        With[{s = Diagonal[svd[[2]]]},
            If[ Length[s] > 1 && Chop[s[[2]], 10 ^ -8] != 0,
                Missing["NotProduct"],
                With[{
                    a = ArrayReshape[svd[[1]][[All, 1]] Sqrt[s[[1]]], {2, 2}],
                    b = ArrayReshape[Conjugate[svd[[3]][[All, 1]]] Sqrt[s[[1]]], {2 ^ (k - 1), 2 ^ (k - 1)}]
                },
                    Replace[kroneckerFactors[b, k - 1], rest_List :> Prepend[rest, a]]
                ]
            ]
        ]
    ]
]
