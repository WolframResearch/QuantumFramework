Package["Wolfram`QuantumFramework`"]

PackageScope[stabilizerMeasurementCircuitQ]
PackageScope[stabilizerMeasurementApply]
PackageScope[stabilizerSamplingOptionsQ]
PackageScope[stabilizerSamplingWithoutMeasurement]



(* ============================================================================ *)
(* Circuits with computational-basis measurements under Method -> "Stabilizer". *)
(*                                                                              *)
(* Exact mode returns <|record -> PauliStabilizer, ...|>: every record the      *)
(* circuit can produce, each with its normalized post-measurement tableau.      *)
(* Sampling mode ("Shots" -> n) returns n records drawn with their exact        *)
(* probabilities, or with "States" -> True n rules record -> tableau. A record has *)
(* one bit per measured wire, sorted by wire and, for repeated measurements of  *)
(* one wire, in time order: the order of the outcomes of the default method.    *)
(*                                                                              *)
(* With a pure stabilizer input, Clifford gates and computational measurements, *)
(* every branch has the same generators and differs only in signs, and whether *)
(* ps["M", q] is random reads only the generators. So a measurement is random   *)
(* in every branch or in none: the first branch tells whether the count doubles,*)
(* and every record has probability 1 / Length[result].                         *)
(* ============================================================================ *)


stabilizerMeasurementCircuitQ[qco_QuantumCircuitOperator] := stabilizerHasMeasurementQ[qco["Operators"]]

stabilizerHasMeasurementQ[ops_List] := AnyTrue[ops, stabilizerHasMeasurementQ]
stabilizerHasMeasurementQ[_QuantumMeasurementOperator] := True
stabilizerHasMeasurementQ[qco_QuantumCircuitOperator] := stabilizerHasMeasurementQ[qco["Operators"]]
stabilizerHasMeasurementQ[_] := False


stabilizerMeasurementFailure[reason_String, template_String, params_List, extra_Association : <||>] :=
    Failure["StabilizerMeasurement", Join[<|"MessageTemplate" -> template, "MessageParameters" -> params, "Reason" -> reason|>, extra]]

(* Confirm and Enclose share one tag, so a Confirm in a helper reaches the     *)
(* Enclose of stabilizerMeasurementApply. Enclose wraps a confirmed Failure in  *)
(* a "ConfirmationFailed" one; the original is returned.                       *)
$stabilizerMeasurementTag = "StabilizerMeasurement";

stabilizerConfirm[expr_] := Confirm[expr, Null, $stabilizerMeasurementTag]

stabilizerReturnFailure[failure_Failure] := Replace[failure["Expression"], Except[_Failure] :> failure]


(* ============================================================================ *)
(* A computational-basis measurement on wires T: its outcome operators, read in *)
(* the computational basis on the sorted wires, are the projectors |b><b| for   *)
(* b = 0 ... 2^|T| - 1, in that order. Outcome b is then both the b-th outcome  *)
(* of the default method (whatever its label: 0/1, or +/- for a named basis)   *)
(* and the bits of b on the sorted wires; an instrument listing the projectors  *)
(* in another order is refused. Building the 2^|T| operators grows as 4^|T|,   *)
(* so an instrument on more than $stabilizerMatrixMaxQubits wires is refused   *)
(* before they are built; one measurement per wire is the same instrument.     *)
(* The verdict depends on the instrument's stored operator and the relative    *)
(* order of its wires, not on which wires they are, so it is remembered under  *)
(* that key: one classification serves QuantumMeasurementOperator[{q}] on      *)
(* every wire q. The instrument's wires must be listed sorted or in an order   *)
(* that is its own inverse ({2, 1}, {3, 2, 1}, ...): for an order such as      *)
(* {2, 3, 1}, "Operators" applies the inverse of the wire permutation, so its  *)
(* projectors cannot be read on the sorted wires.                               *)
(* ============================================================================ *)

$stabilizerMeasurementMemo = <||>;

stabilizerComputationalWires[qmo_QuantumMeasurementOperator] := With[{wires = Sort[qmo["InputOrder"]]},
    Which[
        ! AllTrue[qmo["InputDimensions"], # == 2 &], Missing["NotQubit"],
        Length[wires] > $stabilizerMatrixMaxQubits, Missing["TooLarge"],
        With[{order = Ordering[qmo["InputOrder"]]}, order[[order]] =!= Range[Length[order]]], Missing["WireOrder"],
        True, With[{key = stabilizerMeasurementKey[qmo]},
            If[ Lookup[$stabilizerMeasurementMemo, key, stabilizerRemember[$stabilizerMeasurementMemo, key, stabilizerComputationalQ[qmo, wires]]],
                wires,
                Missing["NotComputationalMeasurement"]
            ]
        ]
    ]
]

(* The stored operator (qmo[[1]], its stored state qmo[[1, 1]]: the "Operator" *)
(* and "State" properties return converted copies) with every wire replaced by *)
(* its rank among the wires.                                                    *)
stabilizerMeasurementKey[qmo_QuantumMeasurementOperator] := With[{
    rank = AssociationThread[# -> Range[Length[#]]] & @ Union[qmo["InputOrder"], qmo["Target"], Select[qmo[[1]]["OutputOrder"], Positive]]
},
    Hash[{qmo[[1, 1]], Replace[qmo[[1]]["InputOrder"], rank, {1}], Replace[qmo[[1]]["OutputOrder"], rank, {1}], Replace[qmo["Target"], rank, {1}]}]
]

stabilizerComputationalQ[qmo_QuantumMeasurementOperator, wires_List] := Sort[qmo["Target"]] === wires &&
    With[{ops = Values[qmo["Operators"]]},
        AllTrue[ops, Sort[#["InputOrder"]] === wires && Sort[#["OutputOrder"]] === wires &] &&
        Length[ops] == 2 ^ Length[wires] &&
        (stabilizerProjectorIndex[#["Sort"]["MatrixRepresentation"]] & /@ ops) === Range[2 ^ Length[wires]]
    ]

(* k when the matrix is |k><k| (1-based, computational basis), else 0; read *)
(* from its nonzero entries only.                                          *)
stabilizerProjectorIndex[m_ ? SquareMatrixQ] := With[{sa = SparseArray[m]},
    If[ Length[sa["NonzeroPositions"]] == 1 && Equal @@ First[sa["NonzeroPositions"]] &&
            stabilizerArrayZeroQ[sa["NonzeroValues"] - 1],
        First[First[sa["NonzeroPositions"]]],
        0
    ]
]
stabilizerProjectorIndex[_] := 0

stabilizerShortcutName[op_] := Replace[QuantumShortcut[op], {one_} :> one]


(* ============================================================================ *)
(* Execution steps: the flattened circuit split into runs of unitary elements, *)
(* each a list of compiled Clifford gates, and single measurements, each a list *)
(* of sorted wires. Every element is checked here, before execution starts.     *)
(* ============================================================================ *)

stabilizerMeasurementSteps[qco_QuantumCircuitOperator] :=
    stabilizerMeasurementStep /@ Split[qco["Flatten"]["Operators"],
        ! (MatchQ[#1, _QuantumMeasurementOperator] || MatchQ[#2, _QuantumMeasurementOperator]) &]

stabilizerMeasurementStep[{qmo_QuantumMeasurementOperator}] := Replace[stabilizerComputationalWires[qmo], {
    wires_List :> "Measure" -> wires,
    Missing["WireOrder"] :> stabilizerMeasurementFailure["WireOrder",
        "`1` lists its wires in an order Method -> \"Stabilizer\" does not read; list them sorted.", {stabilizerShortcutName[qmo]}],
    Missing["NotQubit"] :> stabilizerMeasurementFailure["NotQubit",
        "`1` acts on a wire that is not a qubit.", {stabilizerShortcutName[qmo]}],
    Missing["TooLarge"] :> stabilizerMeasurementFailure["TooLarge",
        "A measurement on `1` wires has more outcomes than Method -> \"Stabilizer\" reads (at most `2` wires); measure the wires one at a time.",
        {Length[qmo["InputOrder"]], $stabilizerMatrixMaxQubits}],
    _ :> stabilizerMeasurementFailure["NotComputationalMeasurement",
        "`1` is not a computational-basis measurement; Method -> \"Stabilizer\" measures in the computational basis only.",
        {stabilizerShortcutName[qmo]}]
}]

(* A wire below 1 is a measurement record (wire 0, -1, ...): a gate controlled *)
(* by a record is not supported.                                              *)
stabilizerMeasurementStep[run_List] /; Min[Catenate[Through[run["InputOrder"]]]] < 1 :=
    stabilizerMeasurementFailure["RecordWire",
        "`1` acts on wire `2`, a measurement record; Method -> \"Stabilizer\" does not support gates controlled by records.",
        With[{op = SelectFirst[run, Min[#["InputOrder"]] < 1 &]}, {stabilizerShortcutName[op], Min[op["InputOrder"]]}]]

stabilizerMeasurementStep[run_List] := With[{plan = stabilizerExecutionPlan[QuantumCircuitOperator[run]]},
    If[ ListQ[plan["Gates"]],
        "Gates" -> plan["Gates"],
        stabilizerMeasurementFailure["NotClifford",
            "`1` cannot be applied by Method -> \"Stabilizer\" in a circuit with measurements: `2`.",
            stabilizerSpecRefusal[SelectFirst[plan["Specs"], ! ListQ[encodeStabilizerGates[{#}]] &]]]
    ]
]

stabilizerSpecRefusal[spec : Missing[_, _QuantumOperator]] := stabilizerRefusal[spec]
stabilizerSpecRefusal[spec_] := {spec, "it is not a Clifford unitary"}


(* ============================================================================ *)
(* Input: a concrete pure stabilizer state on at least the circuit's wires.     *)
(* ============================================================================ *)

stabilizerMeasurementInput[qco_, Automatic] := PauliStabilizer[Max[qco["Arity"], qco["Max"]]]

stabilizerMeasurementInput[qco_, ps_PauliStabilizer] := Which[
    ! ConcretePauliStabilizerQ[ps],
        stabilizerMeasurementFailure["SymbolicSigns",
            "The input tableau has symbolic signs; Method -> \"Stabilizer\" with measurements needs numeric signs.", {}],
    qco["Max"] > ps["Qubits"], ps["PadRight", qco["Max"]],
    True, ps
]

stabilizerMeasurementInput[qco_, qs_QuantumState] := Which[
    ! AllTrue[qs["OutputDimensions"], # == 2 &],
        stabilizerMeasurementFailure["NotQubit", "The input state has wires that are not qubits.", {}],
    qs["PureStateQ"] && qs["InputDimensions"] === {},
        With[{ps = PauliStabilizer[If[qs["VectorQ"], qs, qs["VectorState"]]["Normalized"]]},
            If[ConcretePauliStabilizerQ[ps], stabilizerMeasurementInput[qco, ps], stabilizerNotStabilizerInput[]]
        ],
    True, stabilizerNotStabilizerInput[]
]

stabilizerNotStabilizerInput[] := stabilizerMeasurementFailure["NotStabilizerState",
    "The input state is not a pure qubit stabilizer state.", {}]


(* Sampling options are checked before any other work. *)
stabilizerSamplingOptions[shots_, states_] := Which[
    ! MatchQ[shots, None | (_Integer ? NonNegative)],
        stabilizerMeasurementFailure["InvalidOption", "\"Shots\" must be a nonnegative integer, not `1`.", {shots}],
    ! BooleanQ[states],
        stabilizerMeasurementFailure["InvalidOption", "\"States\" must be True or False, not `1`.", {states}],
    states && shots === None,
        stabilizerMeasurementFailure["StatesWithoutShots", "\"States\" -> True needs \"Shots\".", {}],
    True, {shots, states}
]

(* A circuit without measurements refuses "Shots" or "States" given with a      *)
(* value other than its default, and an invalid "MaxBranches".                  *)
stabilizerSamplingOptionsQ[opts___] := With[{given = Flatten[{opts}]},
    ! MatchQ[FilterRules[given, "Shots" | "States"], {(("Shots" -> None) | ("States" -> False)) ...}] ||
    ! stabilizerMaxBranchesQ[Lookup[given, "MaxBranches", Automatic]]
]

stabilizerSamplingWithoutMeasurement[opts___] := With[{max = Lookup[Flatten[{opts}], "MaxBranches", Automatic]},
    If[ stabilizerMaxBranchesQ[max],
        stabilizerMeasurementFailure["ShotsWithoutMeasurement",
            "\"Shots\" and \"States\" apply to circuits with measurements; this circuit has none.", {}],
        stabilizerMaxBranches[max, 1]
    ]
]


stabilizerMaxBranchesQ[max_] := MatchQ[max, Automatic | _Integer ? Positive]

stabilizerMaxBranches[Automatic, n_Integer] := Floor[2 ^ 30 / (1500 + 140 n)]
stabilizerMaxBranches[max_Integer ? Positive, _] := max
stabilizerMaxBranches[max_, _] := stabilizerMeasurementFailure["InvalidOption",
    "\"MaxBranches\" must be a positive integer or Automatic, not `1`.", {max}]


(* ============================================================================ *)
(* The branch loop. State: {measured wires in time order, <|bits -> tableau|>}. *)
(* ============================================================================ *)

stabilizerBranchStep[_][{wires_, branches_}, "Gates" -> gates_] :=
    {wires, Map[stabilizerConfirm @ applyCompiledFold[#, gates] &, branches]}

stabilizerBranchStep[max_][state_, "Measure" -> qubits_] := Fold[stabilizerMeasureWire[max], state, qubits]

(* The first branch's outcomes tell whether the count doubles; the limit is   *)
(* checked against that count before any other branch is measured.           *)
stabilizerMeasureWire[max_][{wires_, branches_}, q_] := With[{first = stabilizerConfirm @ First[branches]["M", q]},
    With[{count = Length[branches] Length[first]},
        If[ count > max,
            stabilizerConfirm @ stabilizerMeasurementFailure["MaxBranches",
                "Measuring wire `1` would make `2` branches, more than \"MaxBranches\" -> `3`.", {q, count, max},
                <|"Count" -> count|>]
        ];
        {
            Append[wires, q],
            Join @@ KeyValueMap[
                Function[{key, outcomes}, KeyMap[Append[key, #] &, outcomes]],
                Prepend[stabilizerConfirm @ #["M", q] & /@ Rest[branches], First[Keys[branches]] -> first]
            ]
        }
    ]
]

(* Positions of the bits sorted by wire, time order within a wire. *)
stabilizerRecordOrder[wires_List] := Ordering[Transpose[{wires, Range[Length[wires]]}]]


(* ============================================================================ *)
(* Sampling: one tableau per shot. At a measurement, ps["M", q] gives one or two *)
(* outcomes; a random one is drawn with RandomInteger[] from the caller's random *)
(* stream. Records and states come from the same draws, so "States" -> True and *)
(* "States" -> False give the same records for the same seed.                  *)
(* ============================================================================ *)

stabilizerSampleStep[{bits_, ps_}, "Gates" -> gates_] := {bits, stabilizerConfirm @ applyCompiledFold[ps, gates]}

stabilizerSampleStep[state_, "Measure" -> qubits_] := Fold[stabilizerSampleWire, state, qubits]

stabilizerSampleWire[{bits_, ps_}, q_] := With[{outcomes = stabilizerConfirm @ ps["M", q]},
    With[{bit = If[Length[outcomes] == 1, First[Keys[outcomes]], RandomInteger[]]},
        {Append[bits, bit], outcomes[bit]}
    ]
]

(* The gates before the first measurement run once; every shot starts after them. *)
stabilizerSample[init_, steps_List, shots_Integer, states_] := With[{
    prefix = TakeWhile[steps, MatchQ["Gates" -> _]],
    order = stabilizerRecordOrder[Catenate[Cases[steps, ("Measure" -> wires_) :> wires]]]
},
    With[{start = Fold[stabilizerSampleStep, {{}, init}, prefix], rest = Drop[steps, Length[prefix]]},
        Table[
            With[{shot = Fold[stabilizerSampleStep, start, rest]},
                If[states, First[shot][[order]] -> Last[shot], First[shot][[order]]]
            ],
            shots
        ]
    ]
]


stabilizerMeasurementApply[qco_QuantumCircuitOperator, qs_, opts : OptionsPattern[PauliStabilizerApply]] := Enclose[
    With[{unknown = Complement[Keys[Flatten[{opts}]], Keys[Options[PauliStabilizerApply]]]},
        If[ unknown =!= {},
            stabilizerConfirm @ stabilizerMeasurementFailure["UnknownOption", "Unknown Method -> \"Stabilizer\" options: `1`.", {unknown}]
        ]
    ];
    With[{values = OptionValue[PauliStabilizerApply, {opts}, {"Compress", "Shots", "States", "MaxBranches"}]},
        If[ First[values] =!= None,
            stabilizerConfirm @ stabilizerMeasurementFailure["InvalidOption", "\"Compress\" applies to circuits without measurements.", {}]
        ];
        With[{
            sampling = stabilizerConfirm @ stabilizerSamplingOptions[values[[2]], values[[3]]],
            init = stabilizerConfirm @ stabilizerMeasurementInput[qco, qs]
        },
        With[{
            max = stabilizerConfirm @ stabilizerMaxBranches[values[[4]], init["Qubits"]],
            steps = stabilizerConfirm /@ stabilizerMeasurementSteps[qco]
        },
            If[ First[sampling] === None,
                With[{result = Fold[stabilizerBranchStep[max], {{}, <|{} -> init|>}, steps]},
                    KeyMap[Part[#, stabilizerRecordOrder[First[result]]] &, Last[result]]
                ],
                stabilizerSample[init, steps, First[sampling], Last[sampling]]
            ]
        ]]
    ],
    stabilizerReturnFailure,
    $stabilizerMeasurementTag
]
