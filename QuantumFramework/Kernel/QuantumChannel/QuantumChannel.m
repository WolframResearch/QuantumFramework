Package["Wolfram`QuantumFramework`"]

PackageExport["QuantumChannel"]

PackageScope["QuantumChannelQ"]


QuantumChannel::invalidName = "`1` is not a recognized QuantumChannel constructor"
QuantumChannel::invalidArgs = "QuantumChannel constructor `1` did not match any rule"
QuantumChannel::emptyKraus = "An empty list does not define a channel; provide at least one Kraus operator."


quantumChannelQ[QuantumChannel[qo_QuantumOperator /; QuantumOperatorQ[Unevaluated[qo]]]] :=
    qo["OutputQudits"] - qo["InputQudits"] >= 0 && AllTrue[Join[qo["InputOrder"], Drop[qo["OutputOrder"], qo["OutputQudits"] - qo["InputQudits"]]], Positive]

quantumChannelQ[___] := False

QuantumChannelQ[qc_QuantumChannel] := System`Private`HoldValidQ[qc] || quantumChannelQ[Unevaluated[qc]]

QuantumChannelQ[___] := False

qc_QuantumChannel /; System`Private`HoldNotValidQ[qc] && quantumChannelQ[Unevaluated[qc]] := System`Private`HoldSetValid[qc]


(* A channel's operator-sum runs over a nonempty Kraus set, so an empty list is
   not a channel. Reject it before the general list path lets {} reach MapAt and
   surface an internal part error. *)
QuantumChannel[{}, ___] := (
    Message[QuantumChannel::emptyKraus];
    Failure["EmptyKraus", <|"MessageTemplate" :> QuantumChannel::emptyKraus|>]
)

(* A one-element Kraus list is the deterministic channel rho -> M rho M^dagger:
   an isometry with no environment to trace. Route it to the single-operator
   constructor; the general list path below would hand this lone operator a
   dimension-1 environment qudit that collapses, leaving the output wire on the
   non-positive environment label and the channel invalid. *)
QuantumChannel[{opArg_}, args___] := Enclose @ QuantumChannel[ConfirmBy[QuantumOperator[opArg, args], QuantumOperatorQ]["Computational"]]

(* each Kraus operator has to build: one that fails stops the channel with its Failure *)
QuantumChannel[opArgs_List, args___] := Enclose @ Block[{ops = ConfirmBy[QuantumOperator[#, args], QuantumOperatorQ]["Computational"] & /@ opArgs, order, inputDims, outputDims},
    order = Union @@@ Thread[Through[ops["Order"]]];
    inputDims = Merge[AssociationThread[#["InputOrder"], #["InputDimensions"]] & /@ ops, Identity];
    ConfirmAssert[AllTrue[inputDims, Apply[Equal]]];
    outputDims = Merge[AssociationThread[#["OutputOrder"], #["OutputDimensions"]] & /@ ops, Identity];
    ConfirmAssert[AllTrue[outputDims, Apply[Equal]]];
    inputDims = Values[inputDims][[All, 1]];
    outputDims = Values[outputDims][[All, 1]];
    QuantumChannel @ QuantumOperator[
        StackQuantumOperators[#["OrderedOutput", order[[1]], QuditBasis[outputDims]]["OrderedInput", order[[2]], QuditBasis[inputDims]["Dual"]] & /@ ops],
        MapAt[Prepend[0], order, {1}]
    ]
]

QuantumChannel[qo_ ? QuantumOperatorQ] /; ! qo["SortedQ"] := QuantumChannel[qo["Sort"]]

QuantumChannel[qc_ ? QuantumChannelQ, order : _ ? orderQ | Automatic] := qc["Reorder", {Join[Select[qc["FullOutputOrder"], NonPositive], order], order}]

QuantumChannel[qc_ ? QuantumChannelQ, order : {_ ? orderQ | Automatic, _ ? orderQ | Automatic}] := qc["Reorder", order]

QuantumChannel[qc_ ? QuantumChannelQ, args___] := QuantumChannel[QuantumOperator[qc["QuantumOperator"], args]]

QuantumChannel[qm : _ ? QuantumMeasurementOperatorQ | _ ? QuantumMeasurementQ] := QuantumChannel[QuantumOperator[qm["POVM"]]]


(qc_QuantumChannel ? QuantumChannelQ)[qm_QuantumMeasurement, args___] := QuantumMeasurement @ QuantumCircuitOperator[{qm, qc}][args]

(qc_QuantumChannel ? QuantumChannelQ)[op_ ? QuantumFrameworkOperatorQ] := QuantumCircuitOperator[{op, qc}]["QuantumOperator"]

(qc_QuantumChannel ? QuantumChannelQ)[args___] := QuantumCircuitOperator[qc][args]

(* (qc1_QuantumChannel ? QuantumChannelQ)[qc2_ ? QuantumChannelQ] := Enclose @ Module[{
    top, bottom, traceQudits, result
},
    top = qc1["SortOutput"];
    bottom = qc2["SortOutput"];

    traceQudits = qc1["TraceQudits"] + qc2["TraceQudits"];
    top = QuantumOperator[top,
        {Join[1 - Drop[Reverse[Range[traceQudits]], qc2["TraceQudits"]], Drop[top["OutputOrder"], qc1["TraceQudits"]]], top["InputOrder"]}
    ];
    bottom = QuantumOperator[bottom,
        {Join[1 - Take[Reverse[Range[traceQudits]], qc2["TraceQudits"]], Drop[bottom["OutputOrder"], qc2["TraceQudits"]]], bottom["InputOrder"]}
    ];
    result = top[bottom]["SortOutput"];
    QuantumChannel[result]
] *)


(* equality *)

QuantumChannel /: Equal[qc__QuantumChannel] := Equal @@ (#["QuantumOperator"] & /@ {qc})

QuantumChannel /: Unequal[qc__QuantumChannel] := ! Equal[qc]


(* dagger *)

SuperDagger[qc_QuantumChannel] ^:= qc["Adjoint"]

Transpose[qc_QuantumChannel] ^:= qc["Adjoint"]["Conjugate"]


(* simplify *)

Scan[
    (Symbol[#][qc_QuantumChannel, args___] ^:= qc[#, args]) &,
    {"Simplify", "FullSimplify", "Chop", "ComplexExpand"}
]


(* parameterization *)

(qc_QuantumChannel ? QuantumChannelQ)[ps : PatternSequence[p : Except[_Association], ___]] /; ! MemberQ[QuantumChannel["Properties"], p] && Length[{ps}] <= qc["ParameterArity"] :=
    qc[AssociationThread[Take[qc["Parameters"], UpTo[Length[{ps}]]], {ps}]]

(qc_QuantumChannel ? QuantumChannelQ)[rules_ ? AssociationQ] /; ContainsOnly[Keys[rules], qc["Parameters"]] :=
    QuantumChannel[qc["Operator"][rules]]

