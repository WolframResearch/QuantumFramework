Package["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`TensorNetworks`"]

PackageExport["QuantumCircuitOperator"]

PackageScope["QuantumCircuitOperatorQ"]
PackageScope["BarrierQ"]
PackageScope["injectionResets"]
PackageScope["FromCircuitOperatorShorthand"]


QuantumCircuitOperator::invalidName = "`1` is not a recognized QuantumCircuitOperator constructor"

QuantumCircuitOperator::invalidArgs = "QuantumCircuitOperator constructor `1` did not match any rule"

QuantumCircuitOperator::nonpure = "State preparation requires a pure state, but `1` is not pure"


BarrierQ[barrier_] := MatchQ[barrier, "Barrier" | "Barrier"[_ ? orderQ, ___] | "Barrier"[Span[_Integer, _Integer | All], ___]]

quantumCircuitOperatorQ[QuantumCircuitOperator[data : HoldPattern[Association[(_String -> _) ...]]]] :=
    KeyExistsQ[data, "Elements"] && AllTrue[Lookup[data, "Elements"], BarrierQ[#] || QuantumFrameworkOperatorQ[#] &]

QuantumCircuitOperatorQ[___] := False

QuantumCircuitOperatorQ[qco_QuantumCircuitOperator] := System`Private`HoldValidQ[qco] || quantumCircuitOperatorQ[Unevaluated[qco]]

QuantumCircuitOperatorQ[___] := False


qco_QuantumCircuitOperator /; System`Private`HoldNotValidQ[qco] && quantumCircuitOperatorQ[Unevaluated[qco]] := System`Private`HoldSetValid[qco]


Options[QuantumCircuitOperator] := Join[{"Parameters" -> {}, "Picture" -> Automatic}, Options[CircuitDraw]]

(* constructors *)

circuitNamePattern = name_String | name_String[___] /; MemberQ[$QuantumCircuitOperatorNames, name]

FromCircuitOperatorShorthand[barrier_ ? BarrierQ] := barrier
FromCircuitOperatorShorthand[{}] := QuantumCircuitOperator[]
FromCircuitOperatorShorthand[qc_ ? QuantumCircuitOperatorQ] := qc
FromCircuitOperatorShorthand[qc_ ? QuantumCircuitOperatorQ -> order_ ? orderQ] := QuantumCircuitOperator[qc, order]
FromCircuitOperatorShorthand[arg : circuitNamePattern] := QuantumCircuitOperator[arg]
FromCircuitOperatorShorthand[(arg : circuitNamePattern) -> args_] := QuantumCircuitOperator[FromCircuitOperatorShorthand[arg], args]
FromCircuitOperatorShorthand[(arg : circuitNamePattern) -> order_ ? orderQ -> args_] := QuantumCircuitOperator[FromCircuitOperatorShorthand[arg], order, args]
FromCircuitOperatorShorthand[op : Except[_ ? QuantumCircuitOperatorQ, _ ? QuantumFrameworkOperatorQ]] := op
FromCircuitOperatorShorthand[arg_] := Enclose @ Replace[Confirm @ FromOperatorShorthand[arg], ops_List :> QuantumCircuitOperator[ops]]



QuantumCircuitOperator[f_Failure, ___] := f

QuantumCircuitOperator[arg_, order_ ? orderQ, args___] := QuantumCircuitOperator[QuantumCircuitOperator[arg, args], order]

QuantumCircuitOperator[arg_ -> order_ ? orderQ, args___] := QuantumCircuitOperator[QuantumCircuitOperator[arg, args], order]

(* Options[CircuitDraw] joins about 300 drawing and Graphics options, built anew on every call. SetOptions
   changes their values but not their names, so the names are collected once, on first use, to filter the
   options a circuit keeps; this assumes Options[CircuitDraw] is not reassigned after the paclet loads. *)
circuitDrawOptionNames[] := circuitDrawOptionNames[] = Keys[Options[CircuitDraw]]

QuantumCircuitOperator[operators_ ? ListQ, opts : OptionsPattern[]] := Enclose @ With[{parameters = ToList @ OptionValue["Parameters"], picture = OptionValue["Picture"]},
    QuantumCircuitOperator[With[{
        ops = (Confirm[If[BarrierQ[#] || parameters === {} && MatchQ[picture, Automatic | None], #,
            Head[#][#, "Parameters" -> Join[#["ParameterSpec"], parameters], "Picture" -> Replace[picture, Automatic | None :> #["Picture"]]]] & @ FromCircuitOperatorShorthand[#]] & /@ operators)
    },
        <|
            "Elements" -> ops,
            "Label" -> Replace[Composition @@ Reverse @ Through[DeleteCases[ops, _ ? BarrierQ]["Label"]], Identity -> None],
            DeleteDuplicatesBy[First] @ FilterRules[{opts}, circuitDrawOptionNames[]]
        |>
    ]]
]


QuantumCircuitOperator[op : Except[_ ? QuantumCircuitOperatorQ | _ ? ListQ, _ ? QuantumFrameworkOperatorQ], args___] := QuantumCircuitOperator[{op}, args]

QuantumCircuitOperator[op_ ? QuantumCircuitOperatorQ, opts : OptionsPattern[]] := QuantumCircuitOperator[op["Elements"], opts, op["Options"]]

QuantumCircuitOperator[qco_ ? QuantumCircuitOperatorQ | {qco_ ? QuantumCircuitOperatorQ}] := qco

QuantumCircuitOperator[params: Except[{Alternatives @@ $QuantumCircuitOperatorNames, ___}, _List]] :=
    Enclose @ QuantumCircuitOperator[ConfirmBy[QuantumOperator[#], QuantumOperatorQ] & @@ Replace[#, param : Except[_List] :> {param}] & /@ params]

QuantumCircuitOperator[g_ ? GraphQ, opts___] := TensorNetworkQuantumCircuit[ToTensorNetworkGraph[g, opts, Method -> "Random"]]

QuantumCircuitOperator[arg : Except[_List | _Association | _File] /; ! (StringQ[arg] && (qasmStringQ[arg] || qasmFileQ[arg]))] := QuantumCircuitOperator[{arg}]

QuantumCircuitOperator[operators_, label_, opts : OptionsPattern[]] := QuantumCircuitOperator[ToList[operators], "Label" -> label, opts]

QuantumCircuitOperator[] := QuantumCircuitOperator[{}]


(* composition *)

(qco_QuantumCircuitOperator ? QuantumCircuitOperatorQ)[op_ ? QuantumFrameworkOperatorQ] :=
    QuantumCircuitOperator[Prepend[qco["Elements"], op], Replace[{qco["Label"], op["Label"]}, {{None, _} | {_, None} -> None, {lbl1_, lbl2_} :> lbl1[lbl2]}]]

Options[quantumCircuitApply] := {Method -> Automatic}

(* A state injected onto wires that already carry qudits (an operator with no input, on a wire an
   earlier element or the circuit's input uses) replaces those qudits: they are discarded and the
   state is left there, the reset of the wires to the state. A circuit on fixed wires has no other
   place for the displaced qudits. Injections onto fresh wires, or onto wires an earlier element
   ended, are left as they are. An injection onto some wires in use and some fresh ones first
   puts |0> on its fresh wires, so that the reset acts on qudits that are there and the circuit's
   input wires stay what they were. An injection of other dimensions than the qudits it displaces
   traces those qudits out and then leaves the state on the freed wires. *)
injectionResets[qco_QuantumCircuitOperator] := With[{elements = qco["Flatten"]["Elements"]},
    If[ NoneTrue[elements, injectionQ],
        qco,
        QuantumCircuitOperator[
            First @ Fold[
                Function[{acc, el}, With[{live = acc[[2]]},
                    If[ injectionQ[el] && IntersectingQ[el["OutputOrder"], Keys[live]],
                        injectionReset[acc, el["Sort"]],
                        {Append[acc[[1]], el], liveWires[live, el]}
                    ]
                ]],
                {{}, KeySelect[AssociationThread[qco["InputOrder"], qco["InputDimensions"]], Positive]},
                elements
            ],
            qco["Options"]
        ]
    ]
]

(* the elements replacing an injection el (sorted) onto some live wires, and the live wires after *)
injectionReset[{elements_, live_}, el_] := With[{
    wires = el["OutputOrder"], dims = el["OutputDimensions"], used = KeyExistsQ[live, #] & /@ el["OutputOrder"]
},
    {
        Join[
            elements,
            If[ Lookup[live, Pick[wires, used]] === Pick[dims, used],
                Join[
                    If[ And @@ used,
                        {},
                        {QuantumOperator[QuantumState["Register"[Pick[dims, used, False]]], {Pick[wires, used, False], {}}]}
                    ],
                    {QuantumOperator["Reset"[el["State"]], wires]}
                ],
                Append[QuantumOperator["Trace"[live[#]], {#}] & /@ Pick[wires, used], el]
            ]
        ],
        Join[live, AssociationThread[wires, dims]]
    }
]

injectionQ[el_] := QuantumOperatorQ[el] && el["InputDimension"] == 1 && el["OutputDimension"] >= 2

(* the wires carrying a qudit after an element, with their dimensions: it takes its input wires
   and leaves its output wires, so a wire it ends (a post-selection, a discard) is free again *)
liveWires[live_, el_ ? BarrierQ] := live
liveWires[live_, el_] := Join[KeyDrop[live, el["InputOrder"]], KeySelect[AssociationThread[el["OutputOrder"], el["OutputDimensions"]], Positive]]

quantumCircuitApply[qco_QuantumCircuitOperator, qs_QuantumState, opts : OptionsPattern[]] /; qco["InputDimensions"][[Ordering[qco["InputOrder"]]]] == qs["OutputDimensions"] := Replace[
    OptionValue[Method],
    {
        "Schrodinger" :> Block[{n = 1, m = qco["Gates"]},
            Progress`EvaluateWithProgress[
                Fold[(n++; FullSimplify @ #2[#1, Method -> {"TensorNetwork", "Computational" -> False}]) &, qs, injectionResets[qco]["NormalOperators"]],
                <|"Text" -> "Folding gates", "Progress" :> n / m, "Percentage" :> n / m, "ElapsedTime" -> Automatic, "RemainingTime" -> Automatic|>
            ]
        ],
        Automatic | "TensorNetwork" | {"TensorNetwork", subOpts___} :> TensorNetworkApply[injectionResets[qco]["Flatten"], qs, FilterRules[{subOpts}, Options[TensorNetworkApply]]],
        "QuEST" :> QuESTApply[qco, qs],
        "Qiskit" | {"Qiskit", subOpts___} :> qco["Qiskit"][qs, subOpts],
        "Stabilizer" :> PauliStabilizerApply[qco, qs],
        {"Stabilizer", subOpts___} :> PauliStabilizerApply[qco, qs, subOpts],
        _ -> $Failed
    }
]

QuantumCircuitOperator::dim = "Circuit expecting dimensions `1`, but the state has dimensions `2`."

quantumCircuitApply[qco_QuantumCircuitOperator, qs_QuantumState, OptionsPattern[]] :=
    (Message[QuantumCircuitOperator::dim, qco["InputDimensions"], qs["OutputDimensions"]]; $Failed)


(qco_QuantumCircuitOperator ? QuantumCircuitOperatorQ)[qs_QuantumState ? QuantumStateQ, opts : OptionsPattern[quantumCircuitApply]] :=
    quantumCircuitApply[
        qco /* QuantumCircuitOperator[
            Join[
                MapIndexed[
                    With[{idx = #1 - Min[qco["FullInputOrder"]] + 1}, If[idx <= qs["OutputQudits"], "I"[qs["Output"]["Extract", {idx}]], "I"] -> #1] &,
                    qco["FreeOrder"]
                ],
                MapIndexed[
                    If[#2[[1]] <= qs["OutputQudits"], "I"[qs["Output"]["Extract", #2 + Length[qco["FullInputOrder"]]]], "I"] -> #1 &,
                    Max[qco["Order"]] + Range[Max[0, qs["OutputQudits"] - Length[qco["FullInputOrder"]]]]
                ]
            ]
        ],
        If[# === {}, qs, QuantumTensorProduct[qs, QuantumState["Register"[#]]]] & @
            ConstantArray[2, Max[0, Length[qco["FullInputOrder"]] - qs["OutputQudits"]]],
        opts
    ]

(qco_QuantumCircuitOperator ? QuantumCircuitOperatorQ)[opts : OptionsPattern[]] := With[
    {method = OptionValue[quantumCircuitApply, {opts}, Method]},
    Switch[method,
    "Stabilizer",
    PauliStabilizerApply[qco, Automatic],
    {"Stabilizer", ___},
    PauliStabilizerApply[qco, Automatic, Sequence @@ Rest[method]],
    "QuEST",
    QuESTApply[qco, QuantumState["Register"[qco["InputDimensions"]]]],
    _,
    (QuantumCircuitOperator[
        Join[
            MapThread[QuantumState["Register", #2, "Label" -> Ket[{"0"}]] -> {#1} &, {qco["InputOrder"], qco["InputDimensions"]}],
            Map[QuantumState["Register", 2, "Label" -> Ket[{"0"}]] -> {#1} &, Complement[qco["FreeOrder"], qco["InputOrder"]]]
        ]
    ] /* qco)[QuantumState[1, 1], opts]
]
]

(qco_QuantumCircuitOperator ? QuantumCircuitOperatorQ)[qm_QuantumMeasurement, opts : OptionsPattern[]] :=
    QuantumMeasurement[qco[qm["QuantumOperator"]]["QuantumOperator", opts]]

op_QuantumMeasurementOperator[qco_QuantumCircuitOperator ? QuantumCircuitOperatorQ] :=
    QuantumCircuitOperator[Append[qco["Elements"], op], Replace[{op["Label"], qco["Label"]}, {{None, _} | {_, None} -> None, {lbl1_, lbl2_} :> lbl1[lbl2]}]]


(* A composition joins its circuits' elements, which are valid already, so the circuit is built
   from them directly; an argument that is a circuit is used as it is. *)
compositionCircuit[qc_ ? QuantumCircuitOperatorQ] := qc

compositionCircuit[arg_] := QuantumCircuitOperator[arg]

QuantumCircuitOperator /: comp : Composition[___, _QuantumCircuitOperator ? QuantumCircuitOperatorQ, ___] := Enclose @
With[{ops = ConfirmBy[compositionCircuit[#], QuantumCircuitOperatorQ] & /@ List @@ Unevaluated[comp]},
    QuantumCircuitOperator[<|
        "Elements" -> Flatten[Through[Reverse[ops]["Elements"]], 1],
        "Label" -> If[# === {} || MemberQ[#, None], None, Composition @@ #] & @ Through[ops["Label"]]
    |>]
]

QuantumCircuitOperator /: comp : RightComposition[___, _QuantumCircuitOperator ? QuantumCircuitOperatorQ, ___] := Enclose @
With[{ops = ConfirmBy[compositionCircuit[#], QuantumCircuitOperatorQ] & /@ List @@ Unevaluated[comp]},
    QuantumCircuitOperator[<|
        "Elements" -> Flatten[Through[ops["Elements"]], 1],
        "Label" -> If[# === {} || MemberQ[#, None], None, Composition @@ Reverse @ #] & @ Through[ops["Label"]]
    |>]
]


(* dagger *)

SuperDagger[qco_QuantumCircuitOperator] ^:= qco["Dagger"]

SuperStar[qco_QuantumCircuitOperator] ^:= qco["Conjugate"]

Transpose[qco_QuantumCircuitOperator] ^:= qco["Transpose"]


(* equality *)

QuantumCircuitOperator /: Equal[left___, qco_QuantumCircuitOperator, right___] :=
    Equal @@ (If[QuantumCircuitOperatorQ[#], #["CircuitOperator"], #] & /@ {left, qco, right})

QuantumCircuitOperator /: Unequal[left___, qco_QuantumCircuitOperator, right___] := ! Equal[left, qco, right]


(* part *)

Part[qco_QuantumCircuitOperator, i_, rest___] ^:= With[{elem = ToList[qco["Elements"][[i]]]}, If[QuantumCircuitOperatorQ[elem], elem[[rest]], QuantumCircuitOperator[elem, qco["Options"]][[rest]]]]

Part[qco_QuantumCircuitOperator] ^:= qco


(* reorder *)

QuantumCircuitOperator[qc_ ? QuantumCircuitOperatorQ, order_ ? orderQ] := QuantumCircuitOperator[qc, {order, order}]

QuantumCircuitOperator[qc_ ? QuantumCircuitOperatorQ, order : {_ ? orderQ, _ ? orderQ}] := Block[{
    outOrder = Select[qc["FullOutputOrder"], Positive],
    inOrder = Select[qc["FullInputOrder"], Positive],
    length,
    outRepl, inRepl
},
    length = Max[Length[outOrder], Length[inOrder]];
    outRepl = Thread[Sort @ Take[Join[outOrder, inOrder], UpTo[length]] -> Take[Join[order[[1]], order[[2]], Drop[outOrder, UpTo[Length[order[[1]]]]]], UpTo[length]]];
    inRepl = Thread[Sort @ Take[Join[inOrder, outOrder], UpTo[length]] -> Take[Join[order[[2]], order[[1]], Drop[inOrder, UpTo[Length[order[[2]]]]]], UpTo[length]]];
    {outRepl, inRepl} = {Join[outRepl, inRepl], Join[inRepl, outRepl]};
    QuantumCircuitOperator[
        Which[
            BarrierQ[#], # /. inRepl,
            True, With[{newOrder = {#["FullOutputOrder"] /. outRepl, #["FullInputOrder"] /. inRepl}},
                If[QuantumOperatorQ[#], #["Reorder", newOrder], Head[#][#, newOrder]]
            ]
        ] & /@ qc["Elements"],
        qc["Options"]
    ]
]

(* simplify *)

Scan[
    (Symbol[#][qco_QuantumCircuitOperator, args___] ^:= qco[#, args]) &,
    {"Simplify", "FullSimplify", "Chop", "ComplexExpand"}
]


(* parameterization *)

(qco_QuantumCircuitOperator ? QuantumCircuitOperatorQ)[rules_ ? AssociationQ] :=
    QuantumCircuitOperator[
        If[BarrierQ[#] || #["ParameterArity"] == 0, #, #[KeyTake[rules, #["Parameters"]]]] & /@ qco["Elements"],
        qco["Label"] /. rules
    ] /; ContainsOnly[Keys[rules], qco["Parameters"]]

(qco_QuantumCircuitOperator ? QuantumCircuitOperatorQ)[parameters : PatternSequence[p : Except[_Association], ___]] /; ! MemberQ[QuantumCircuitOperator["Properties"], p] /; Length[{parameters}] <= Length[qco["Parameters"]] :=
    qco[AssociationThread[Take[qco["Parameters"], UpTo[Length[{parameters}]]], {parameters}]]

