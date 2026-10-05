(* ::Package:: *)

Package["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`Arrays`"]

PackageImport["Wolfram`QuantumFramework`Gates`"]


PackageExport["QuantumOperator"]

PackageScope["QuantumOperatorQ"]
PackageScope["StackQuantumOperators"]
PackageScope["$QuantumOperatorBroadcastLimit"]
PackageScope["matrixMapAmplitudes"]
PackageScope["lazyMatrixMapAmplitudesQ"]
PackageScope["padQuantumOperators"]
PackageScope["padJumpOperators"]
PackageScope["jumpRatesQ"]


(* What the matrix constructors below accept: a rank-2 array container of any
   tier.  MatrixQ was the test, and it is narrower than the constructors need -
   it is False for a NumericArray, for a MatrixSymbol and for a deferred
   contraction, all of which are perfectly good operator matrices.  Anything it
   rejected fell through to the Diagonal catch-all at the foot of this file,
   which reads a non-atomic argument as a SCALAR eigenvalue: an entire matrix
   became HoldForm[m] * IdentityMatrix, silently. *)

matrixContainerQ[m_] := ArrayContainerQ[m] && MatchQ[ArrayDimensions[m], {_Integer, _Integer}]


QuantumOperator::invalidInputOrder = "input order should be a list of distinct input qudit positions"
QuantumOperator::invalidOutputOrder = "output order should be a list of distinct output qudit positions"
QuantumOperator::invalidName = "`1` is not a recognized QuantumOperator constructor"
QuantumOperator::invalidArgs = "QuantumOperator constructor `1` did not match any rule"
QuantumOperator::broadcast = "broadcasting a `1`-qudit operator over an order of length `2` would build a tensor power of dimension `3`, exceeding `4`; supply an operator and order of matching size"
QuantumOperator::rates = "the rates `1` do not fit the jump operators: give as many rates as there are jump operators, `2`, or a `2` by `2` rate matrix"


(* largest implied dimension the order-driven multiplicity broadcast may materialize *)
$QuantumOperatorBroadcastLimit = 2^24

quantumOperatorQ[QuantumOperator[qs_QuantumState /; QuantumStateQ[Unevaluated[qs]], {_ ? orderQ, _ ? orderQ}]] := True

quantumOperatorQ[___] := False


QuantumOperatorQ[qo_QuantumOperator] := System`Private`HoldValidQ[qo] || quantumOperatorQ[Unevaluated[qo]]

QuantumOperatorQ[___] := False


qo_QuantumOperator /; System`Private`HoldNotValidQ[qo] && quantumOperatorQ[Unevaluated[qo]] := System`Private`HoldSetValid[qo]


(* constructors *)

SetAttributes[QuantumOperator, NHoldRest]

QuantumOperator[arg : _ ? QuantumStateQ, order : {_ ? orderQ, _ ? orderQ}, opts__] :=
    Enclose @ QuantumOperator[ConfirmBy[QuantumState[arg, opts], QuantumStateQ], order]

QuantumOperator[qs_ ? QuantumStateQ, order : {outputOrder_ ? orderQ, inputOrder_ ? orderQ}] /;
    qs["Qudits"] == Length[outputOrder] + Length[inputOrder] && (qs["OutputQudits"] != Length[outputOrder] || qs["InputQudits"] != Length[inputOrder]) :=
    QuantumOperator[qs["SplitDual", Length[outputOrder]], order]

QuantumOperator[qs_ ? QuantumStateQ] :=
    QuantumOperator[
        qs,
        If[qs["InputQudits"] > qs["OutputQudits"], {Automatic, Range[qs["InputQudits"]]}, {Range[qs["OutputQudits"]], Automatic}]
    ]

QuantumOperator[arg : _ ? QuantumStateQ, outputOrder : _ ? orderQ | Automatic, inputOrder : _ ? orderQ | Automatic, opts___] :=
    QuantumOperator[QuantumState[arg, opts], {outputOrder, inputOrder}]

QuantumOperator[qs : _ ? QuantumStateQ, autoOrder : _ ? orderQ | Automatic | {Automatic, Automatic}, opts___] := With[{
    order = Replace[autoOrder, Automatic | {Automatic, Automatic} :> Range[qs["OutputQudits"]]]
},
    QuantumOperator[qs, If[Length[order] == qs["OutputQudits"], {order, Automatic}, {Automatic, order}], opts]
]


QuantumOperator[qs_ ? QuantumStateQ, {Automatic, order_ ? orderQ}, opts___] :=
    QuantumOperator[
        qs,
        {
            # - Min[#, 1] + 1 & @ Reverse @ Take[Join[Reverse @ order, Min[order] - Range[Max[1, qs["OutputQudits"] - qs["InputQudits"]]]], UpTo @ qs["OutputQudits"]],
            order
        },
        opts
    ]

QuantumOperator[qs_ ? QuantumStateQ, {order_ ? orderQ, Automatic}, opts___] :=
    QuantumOperator[
        qs,
        {
            order,
            # - Min[#, 1] + 1 & @ Reverse @ Take[Join[Reverse @ order, Min[order] - Range[Max[1, qs["InputQudits"] - qs["OutputQudits"]]]], UpTo @ qs["InputQudits"]]
        },
        opts
    ]

QuantumOperator[qs_ ? QuantumStateQ, opts : PatternSequence[Except[{_ ? orderQ, _ ? orderQ}], ___]] := QuantumOperator[
    QuantumOperator[qs],
    QuantumBasis[qs["Basis"], opts]
]

QuantumOperator[qs_ ? QuantumStateQ, {outputOrder_ ? orderQ, inputOrder_}] /; qs["OutputDimension"] == 1 && Length[outputOrder] > 0 :=
    QuantumOperator[qs, {{}, inputOrder}]

QuantumOperator[qs_ ? QuantumStateQ, {outputOrder_, inputOrder_ ? orderQ}] /; qs["InputDimension"] == 1 && Length[inputOrder] > 0 :=
    QuantumOperator[qs, {outputOrder, {}}]


QuantumOperator[qb : _QuditBasis | _QuantumBasis, opts___] := QuantumOperator[QuantumState[qb], opts]


QuantumOperator[tensor_ ? TensorQ /; TensorRank[tensor] > 2, order : _ ? autoOrderQ : Automatic, args___, opts : OptionsPattern[]] := Block[{
    dimensions = TensorDimensions[tensor],
    outputOrder,
    basis,
    inputDimension, outputDimension
},
    basis = QuantumBasis[args];
    If[ Times @@ dimensions === basis["Dimension"],
        dimensions = basis["Dimensions"]
    ];
    outputOrder = Replace[order, {
        Automatic :> Range[Min[basis["OutputQudits"], Length[dimensions]]],
        {o_ ? orderQ, _} :> o,
        {Automatic, o_ ? orderQ} :> Range[Max[Length[dimensions] - Length[o], 0]]}
    ];
    {outputDimension, inputDimension} = Times @@@ TakeDrop[dimensions, Length[outputOrder]];
    If[ basis["OutputDimension"] != outputDimension,
        basis = QuantumBasis[basis, "Output" -> QuditBasis[dimensions[[;; Length[outputOrder]]]]]
    ];
    If[ basis["InputDimension"] != inputDimension,
        basis = QuantumBasis[basis, "Input" -> QuditBasis[dimensions[[Length[outputOrder] + 1 ;;]]]["Dual"]]
    ];
    QuantumOperator[
        ArrayReshape[tensor, {outputDimension, inputDimension}],
        order,
        basis,
        opts
    ]
]

QuantumOperator[assoc_Association, order : (_ ? orderQ) : {1}, args___, opts : OptionsPattern[]] := Enclose @ Module[{
    quditBasis,
    basis,
    tensorDimensions
},
    quditBasis = QuditBasis[
        Association @ Catenate @ MapIndexed[
            With[{counts = #1, i = First[#2]}, MapIndexed[{QuditName[#1], i} -> UnitVector[Length[counts], First[#2]] &, Keys @ counts]] &,
            Counts /@ Transpose[ConfirmBy[List @@@ Keys[assoc], MatrixQ]]
        ],
        args
    ];
    ConfirmAssert[Length[assoc] > 0 && Equal @@ TensorDimensions /@ assoc];
    ConfirmAssert[EvenQ @ quditBasis["Qudits"]];
    basis = QuantumBasis[
        "Output" -> QuantumPartialTrace[quditBasis, Range[quditBasis["Qudits"] / 2 + 1, quditBasis["Qudits"]]],
        "Input" -> QuantumPartialTrace[quditBasis, Range[quditBasis["Qudits"] / 2]]["Dual"]
    ];
    tensorDimensions = TensorDimensions @ First[assoc];
    QuantumOperator[
        ArrayReshape[Lookup[KeyMap[QuditName[List @@ #] &, assoc], quditBasis["Names"], 0], Join[basis["Dimensions"], tensorDimensions]],
        Join[Complement[Range[Length[tensorDimensions]], order], order],
        basis,
        opts
    ]
]


QuantumOperator::invalidState = "invalid state specification";


(* d as n qudit dimensions: equal when d is an n-th power, otherwise its smallest
   divisors first and the rest last *)
splitDimension[1, _] := {1}

splitDimension[d_, 1] := {d}

splitDimension[d_, n_] := With[{root = Round[d ^ (1 / n)]},
    If[ root ^ n == d,
        ConstantArray[root, n],
        With[{factors = primeFactors[d]}, Append[Take[factors, n - 1], Times @@ Drop[factors, n - 1]]]
    ]
]

(* the fewest copies of a matrix that fill its orders: each order length divides into at
   most as many qudits as its side's dimension has prime factors *)
matrixBroadcastMultiplicity[dims_, orders_] := With[{
    sides = Select[Transpose[{dims, Length /@ orders}], First[#] > 1 &]
},
    If[ sides === {},
        1,
        SelectFirst[
            Divisors[GCD @@ sides[[All, 2]]],
            k |-> AllTrue[sides, #[[2]] / k <= Length[primeFactors[#[[1]]]] &]
        ]
    ]
]

(* With no basis given, the order says how many qudits each side of the matrix spans, and
   a matrix too small for its order is broadcast over it, as a named operator is. A side
   of dimension 1 spans none. *)
QuantumOperator[matrix_ ? matrixContainerQ, order : _ ? orderQ | {_ ? orderQ | Automatic, _ ? orderQ | Automatic}, opts : OptionsPattern[]] := Enclose @ Module[{
    dims = ArrayDimensions[matrix], orders, multiplicity, op
},
    orders = MapThread[
        Which[#1 == 1, {}, MatchQ[#2, Automatic | {}], Range[Length[primeFactors[#1]]], True, #2] &,
        {dims, Replace[order, o_ ? orderQ :> {o, o}]}
    ];
    multiplicity = matrixBroadcastMultiplicity[dims, orders];
    If[ MissingQ[multiplicity],
        (* no number of copies fills both orders: each side spans as many qudits as it can *)
        orders = MapThread[Take[#2, UpTo[Length[primeFactors[#1]]]] &, {dims, orders}];
        multiplicity = 1
    ];
    op = ConfirmBy[
        QuantumOperator[matrix, QuantumBasis @@ MapThread[splitDimension[#1, Length[#2] / multiplicity] &, {dims, orders}]],
        QuantumOperatorQ
    ];
    Which[
        multiplicity == 1,
        QuantumOperator[op["State"], orders, opts],
        op["Dimension"] ^ multiplicity > $QuantumOperatorBroadcastLimit,
        Message[QuantumOperator::broadcast,
            Max[op["OutputQudits"], op["InputQudits"]], Max[Length /@ orders],
            op["Dimension"] ^ multiplicity, $QuantumOperatorBroadcastLimit
        ];
        $Failed,
        True,
        QuantumOperator[{op, multiplicity}, orders, opts]
    ]
]

QuantumOperator[matrix_ ? matrixContainerQ, order : _ ? autoOrderQ, args___, opts : OptionsPattern[]] := Enclose @ Module[{
    op = ConfirmBy[QuantumOperator[matrix, args], QuantumOperatorQ],
    newOutputOrder, newInputOrder
},
    {newOutputOrder, newInputOrder} = Replace[order, {
        out_ ? orderQ /; Length[out] == op["InputQudits"] :> {out, out},
        out_ ? orderQ :> {out, op["InputOrder"]},
        Automatic :> op["Order"],
        {out : _ ? orderQ | Automatic, in : _ ? orderQ | Automatic} :> {Replace[out, Automatic :> op["OutputOrder"]], Replace[in, Automatic :> op["InputOrder"]]}
    }];
    QuantumOperator[op["State"]["Split", Length[newOutputOrder]], {newOutputOrder, newInputOrder}, opts]
]

QuantumOperator[matrix_ ? matrixContainerQ, opts : OptionsPattern[]] := QuantumOperator[matrix, QuantumBasis[primeFactors[#1], primeFactors[#2]] & @@ ArrayDimensions[matrix], opts]

QuantumOperator[matrix_ ? matrixContainerQ, args__, opts : OptionsPattern[]] := Block[{
    outMultiplicity, inMultiplicity, result
},
    result = Enclose @ Module[{newMatrix = matrix, outputs, inputs,
        basis, newOutputQuditBasis, newInputQuditBasis, state},
        {outputs, inputs} = ArrayDimensions[newMatrix];
        basis = ConfirmBy[QuantumBasis[args], QuantumBasisQ, "Invalid basis"];
        If[ basis["Dimension"] =!= outputs * inputs,
            If[ basis["InputDimension"] == 1,
                basis = QuantumBasis[basis, "Input" -> basis["Output"]["Dual"]]
            ];
            outMultiplicity = If[basis["OutputDimension"] === 1, Indeterminate, Log[basis["OutputDimension"], outputs]];
            inMultiplicity = If[basis["InputDimension"] === 1, Indeterminate, Log[basis["InputDimension"], inputs]];
            If[
                IntegerQ[outMultiplicity] && IntegerQ[inMultiplicity],
                (* multiply existing basis *)

                basis = QuantumBasis[basis,
                    "Output" -> QuditBasis[basis["Output"], outMultiplicity],
                    "Input" -> QuditBasis[basis["Input"], inMultiplicity]
                ],

                (* add one extra qudit *)
                newOutputQuditBasis = QuantumTensorProduct[basis["Output"], QuditBasis[Ceiling[outputs / basis["OutputDimension"]] /. 1 -> Sequence[]]];
                newInputQuditBasis = QuantumTensorProduct[basis["Input"], QuditBasis[Ceiling[inputs / basis["InputDimension"]] /. 1 -> Sequence[]]];

                newMatrix = kroneckerProduct[
                    newMatrix,
                    With[{outs = Ceiling[newOutputQuditBasis["Dimension"] / outputs], ins = Ceiling[newInputQuditBasis["Dimension"] / inputs]},
                        Replace[{} -> {{1}}] @ identityMatrix[Max[outs, ins]][[;; outs, ;; ins]]
                    ]
                ];
                basis = QuantumBasis[basis,
                    "Output" -> newOutputQuditBasis,
                    "Input" -> If[newInputQuditBasis["DualQ"], newInputQuditBasis, newInputQuditBasis["Dual"]]
                ];
            ]
        ];
        state = ConfirmBy[
            QuantumState[
                ArrayVector[newMatrix],
                basis
            ],
            QuantumStateQ,
            Message[QuantumOperator::invalidState]
        ];
        QuantumOperator[state, {Range[basis["FullOutputQudits"]], Range[basis["FullInputQudits"]]}, opts]
    ];
    result /; !FailureQ[Unevaluated @ result]
]

QuantumOperator[array_ ? NumericArrayQ, args___] := QuantumOperator[Normal @ array, args]


QuantumOperator[n_Integer, args___] := QuantumOperator["PhaseShift"[n], args]


QuantumOperator[Labeled[arg_, label_], opts___] := QuantumOperator[arg, opts, "Label" -> label]

QuantumOperator[{}, opts___] := QuantumOperator[QuantumState[{}, 0], opts]

QuantumOperator[arg_List, opts___] := QuantumOperator[QuantumState[Flatten[arg]], opts]

QuantumOperator[arg_, order1 : _ ? orderQ -> order2 : _ ? orderQ, opts___] :=
    QuantumOperator[arg, {order2, order1}, opts]


(* Mutation *)

QuantumOperator[qo_ ? QuantumOperatorQ, order_ ? orderQ] /;
    Length[order] == Length[Union[qo["FullOutputOrder"], qo["FullInputOrder"]]] :=
    With[{repl = Thread[Union[qo["FullOutputOrder"], qo["FullInputOrder"]] -> order]},
        QuantumOperator[qo, {qo["FullOutputOrder"] /. repl, qo["FullInputOrder"] /. repl}]
    ]

QuantumOperator[qo_ ? QuantumOperatorQ, order : (_ ? orderQ | Automatic)] :=
    QuantumOperator[qo, {order, order}]

QuantumOperator[qo_ ? QuantumOperatorQ, outputOrder : (_ ? orderQ | Automatic), inputOrder : (_ ? orderQ | Automatic)] :=
    QuantumOperator[qo, {outputOrder, inputOrder}]

QuantumOperator[qo_ ? QuantumOperatorQ, order : {order1 : _ ? orderQ | Automatic, order2 : _ ? orderQ | Automatic}, opts : OptionsPattern[]] := Block[{
    outputOrder = Replace[order1, Automatic -> qo["OutputOrder"]],
    inputOrder = Replace[order2, Automatic -> qo["InputOrder"]],
    outputQudits = Max[qo["FullOutputQudits"], 1],
    inputQudits = Max[qo["FullInputQudits"], 1]
},
    With[{
        multiplicity = LCM[Quotient[Length[outputOrder], outputQudits], Quotient[Length[inputOrder], inputQudits]]
    },
        If[ qo["Dimension"] ^ multiplicity > $QuantumOperatorBroadcastLimit,
            Message[QuantumOperator::broadcast,
                Max[outputQudits, inputQudits], Max[Length[outputOrder], Length[inputOrder]],
                qo["Dimension"] ^ multiplicity, $QuantumOperatorBroadcastLimit
            ];
            $Failed,
            QuantumOperator[{qo, multiplicity}, {outputOrder, inputOrder}, opts]
        ]
    ] /;
        Length[outputOrder] > outputQudits && Divisible[Length[outputOrder], outputQudits] &&
        Length[inputOrder] > inputQudits && Divisible[Length[inputOrder], inputQudits]
]


QuantumOperator[qo_ ? QuantumOperatorQ, order : {_ ? orderQ | Automatic, _ ? orderQ | Automatic}] := qo["Reorder", order, True]

QuantumOperator[qo_ ? QuantumOperatorQ, order1 : _ ? orderQ | Automatic -> order2 : _ ? orderQ | Automatic, opts___] :=
    QuantumOperator[qo, {order2, order1}, opts]

QuantumOperator[qo_ ? QuantumOperatorQ, opts : PatternSequence[Except[_ ? QuantumBasisQ], ___],
    outputOrder : (_ ? orderQ | Automatic), inputOrder : (_ ? orderQ | Automatic)] := Enclose @
    QuantumOperator[qo, {outputOrder, inputOrder}, ConfirmBy[QuantumBasis[qo["Basis"], opts], QuantumBasisQ]]

QuantumOperator[qo_ ? QuantumOperatorQ, order : (_ ? orderQ | Automatic), opts : PatternSequence[Except[_ ? QuantumBasisQ], ___]] := Enclose @
    QuantumOperator[qo, {order, order}, ConfirmBy[QuantumBasis[qo["Basis"], opts], QuantumBasisQ]]

QuantumOperator[qo_ ? QuantumOperatorQ, order : _ ? autoOrderQ, opts : PatternSequence[Except[_ ? QuantumBasisQ], ___]] := Enclose @
    QuantumOperator[qo, order, ConfirmBy[QuantumBasis[qo["Basis"], opts], QuantumBasisQ]]

QuantumOperator[qo_ ? QuantumOperatorQ, opts : PatternSequence[Except[_ ? autoOrderQ | _ ? QuantumBasisQ | _ ? QuantumOperatorQ], ___]] := Enclose @
    QuantumOperator[qo, qo["Order"], ConfirmBy[QuantumBasis[qo["Basis"], opts], QuantumBasisQ]]

QuantumOperator[{qo : _ ? QuantumOperatorQ, multiplicity_Integer ? Positive}] := QuantumOperator[{qo, multiplicity}, Range[multiplicity]]

QuantumOperator[{qo : _ ? QuantumOperatorQ, multiplicity_Integer ? Positive}, opts__] :=
    QuantumOperator[QuantumTensorProduct @ Table[qo, multiplicity], opts]

QuantumOperator[qc_ ? QuantumCircuitOperatorQ, opts___] := QuantumOperator[qc["QuantumOperator"], opts]

QuantumOperator[q : _ ? QuantumChannelQ | _ ? QuantumMeasurementOperatorQ | _ ? QuantumMeasurementQ, opts___] := QuantumOperator[q["Operator"], opts]

(* A Failure from a nested constructor propagates as-is; the diagonal catch-all
   below would otherwise absorb it as a scalar eigenvalue. *)
QuantumOperator[failure_Failure, ___] := failure

QuantumOperator[x : Except[_ ? QuantumStateQ | _ ? QuantumOperatorQ | _ ? QuantumCircuitOperatorQ | _ ? QuantumGateQ], args___] := Enclose @
    ConfirmBy[QuantumOperator["Diagonal"[If[AtomQ[x], x, HoldForm[x]]], args], QuantumOperatorQ]


(* change of basis *)

QuantumOperator[qo_ ? QuantumOperatorQ, name_ ? nameQ, opts___] := QuantumOperator[qo, QuantumBasis[name], opts]

QuantumOperator[qo_ ? QuantumOperatorQ, qb_ ? QuantumBasisQ, opts___] := QuantumOperator[qo, qo["Order"], qb, opts]

QuantumOperator[qo_ ? QuantumOperatorQ, order : _ ? autoOrderQ, qb_ ? QuantumBasisQ, opts___] :=
Enclose @ Block[{
    newBasis
},
    If[qo["Basis"] == qb && qo["Order"] === order, Return[QuantumOperator[QuantumState[qo["State"]["State"], qb], order]]];

    newBasis = If[
        qb["InputDimension"] == 1 && qo["InputDimension"] > 1,
        QuantumBasis[qb, "Input" -> qb["Output"]["Dual"]],
        QuantumBasis[qb]
    ];

    newBasis = QuantumBasis[qb,
        "Output" -> Confirm @ QuditBasis[qo["Output"], newBasis["Output"]],
        "Input" -> Confirm @ QuditBasis[qo["Input"], newBasis["Input"]],
        opts
    ];

    ConfirmAssert[qo["Dimension"] == newBasis["Dimension"], "Basis dimensions are inconsistent"];
    QuantumOperator[
        ConfirmBy[
            QuantumState[
                qo["State"],
                newBasis
            ],
            QuantumStateQ,
            Message[QuantumOperator::invalidState]
        ],
        order
    ]
]

QuantumOperator[qo_ ? QuantumOperatorQ,
    outputOrder : _ ? orderQ | Automatic : Automatic, inputOrder : _ ? orderQ | Automatic : Automatic, qb_ ? QuantumBasisQ, opts___] :=
    QuantumOperator[qo, {outputOrder, inputOrder}, qb, opts]

(* composition *)

(qo_QuantumOperator ? QuantumOperatorQ)[qb_ ? QuantumBasisQ] := With[{states = qo /@ qb["BasisStates"]},
    QuantumBasis[
        AssociationThread[First[states]["Names"], Through[states["StateVector"]]],
        "Label" -> Replace[qo["Label"] @* qb["Label"], Identity | None @* _ | _ @* None -> None],
        qb["Options"]
    ]
]

QuantumOperator::incompatiblePictures = "Pictures `` and `` are incompatible with this operation"

(qo_QuantumOperator ? QuantumOperatorQ)[qs_ ? QuantumStateQ, opts___] := QuantumCircuitOperator[qo][qs, opts]

(* The direct contraction rules below glue qo1's input wires to qo2's output wires with an
   identity factor, which is faithful only when both sides carry the same basis (a dual flag
   apart) on every shared wire. A mismatched basis falls through to the circuit fallback, which
   changes basis through the computational basis, so rep(qo1 @ qo2) = rep(qo1) . rep(qo2) always. *)
composableBasesQ[qo1_, qo2_] := With[{shared = Intersection[qo1["InputOrder"], qo2["OutputOrder"]]},
    shared === {} ||
        (* fully-aligned wires compare as whole bases, skipping the per-qudit walk *)
        qo1["InputOrder"] === shared === qo2["OutputOrder"] &&
            With[{in = qo1["Input"], out = qo2["Output"]}, in === out || in["Dual"] === out] ||
        And @@ MapThread[
            #1 === #2 || #1["Dual"] === #2 &,
            {
                qo1["Input"]["Decompose"][[Lookup[AssociationThread[qo1["InputOrder"], Range[Length[qo1["InputOrder"]]]], shared]]],
                qo2["Output"]["Decompose"][[Lookup[AssociationThread[qo2["OutputOrder"], Range[Length[qo2["OutputOrder"]]]], shared]]]
            }
        ]
]

(qo1_QuantumOperator ? QuantumOperatorQ)[qo2_QuantumOperator ? QuantumOperatorQ] /; (
    (qo1["MatrixQ"] || qo2["MatrixQ"]) && qo1["Picture"] === qo2["Picture"] &&
    Intersection[qo1["OutputOrder"], Complement[qo2["OutputOrder"], qo1["InputOrder"]]] === {} &&
    Intersection[Complement[qo1["InputOrder"], qo2["OutputOrder"]], qo2["InputOrder"]] === {} &&
    composableBasesQ[qo1, qo2]
) := With[{shift = Max[qo1["Order"], qo2["Order"]] + 1},
    qo1["Bend", shift][qo2["Bend", shift]]["Unbend"]
]


(qo1_QuantumOperator ? QuantumOperatorQ)[qo2_QuantumOperator ? QuantumOperatorQ] /; (
    qo1["VectorQ"] && qo2["VectorQ"] && qo1["Picture"] === qo2["Picture"] &&
    Intersection[qo1["OutputOrder"], Complement[qo2["OutputOrder"], qo1["InputOrder"]]] === {} &&
    Intersection[Complement[qo1["InputOrder"], qo2["OutputOrder"]], qo2["InputOrder"]] === {} &&
    composableBasesQ[qo1, qo2]
) := Module[{
    s1, s2, q1OutOrder, q1InOrder, q2OutOrder, q2InOrder,
    q1OutQ, q1InQ, q2OutQ, q2InQ,
    shared, posQ1In, posQ2Out, sharedQ1Idx, sharedQ2Idx,
    keepQ1InIdx, keepQ2OutIdx,
    contractPairs, resultTensor, resultOutOrder, resultInOrder,
    resultOutDecompose, resultInDecompose
},
    s1 = qo1["Sort"]; s2 = qo2["Sort"];
    q1OutOrder = s1["OutputOrder"]; q1InOrder = s1["InputOrder"];
    q2OutOrder = s2["OutputOrder"]; q2InOrder = s2["InputOrder"];
    q1OutQ = Length[q1OutOrder]; q1InQ = Length[q1InOrder];
    q2OutQ = Length[q2OutOrder]; q2InQ = Length[q2InOrder];

    shared = Intersection[q1InOrder, q2OutOrder];
    posQ1In = AssociationThread[q1InOrder, Range[q1InQ]];
    posQ2Out = AssociationThread[q2OutOrder, Range[q2OutQ]];
    sharedQ1Idx = Lookup[posQ1In, shared];
    sharedQ2Idx = Lookup[posQ2Out, shared];
    keepQ1InIdx = Complement[Range[q1InQ], sharedQ1Idx];
    keepQ2OutIdx = Complement[Range[q2OutQ], sharedQ2Idx];

    contractPairs = Transpose[{q1OutQ + sharedQ1Idx, q1OutQ + q1InQ + sharedQ2Idx}];
    (* Two array containers contract through ArrayContract, which never forms the
       outer product: TensorProduct[t1, t2] of two rank-4 operator tensors is rank
       8 with d^8 entries, materialized in full before any index is contracted.
       ArrayContract also joins the operand tiers, where a raw TensorProduct would
       treat a symbolic or lazy operand as a SCALAR and multiply, the unevaluated
       expression then reading downstream as a single amplitude.

       The operands go to ArrayContract as an inactive tensor product, which names
       them one by one; a plain list of them would read as a single ragged array.
       A scalar operand is not a container and stays on the TensorProduct branch,
       and so does the case with nothing to contract, which is a plain tensor
       product. *)
    resultTensor = With[{t1 = s1["StateTensor"], t2 = s2["StateTensor"]},
        If[ ArrayContainerQ[t1] && ArrayContainerQ[t2] && Length[contractPairs] > 0,
            ArrayContract[Inactive[TensorProduct][t1, t2], contractPairs],
            If[Length[contractPairs] == 0,
                TensorProduct[t1, t2],
                TensorContract[TensorProduct[t1, t2], contractPairs]
            ]
        ]
    ];

    (* After contraction, surviving axes appear in tensor order:
         [q1.outputs] [q1.inputs kept] [q2.outputs kept] [q2.inputs]
       Permute so all outputs come first, then all inputs:
         [q1.outputs] [q2.outputs kept] [q1.inputs kept] [q2.inputs] *)
    With[{
        n1Out = q1OutQ, n1In = Length[keepQ1InIdx],
        n2Out = Length[keepQ2OutIdx], n2In = q2InQ
    },
        (* ArrayTranspose rather than Transpose, so a non-explicit result stays
           a container: Transpose on a symbolic node reaches the expression
           tree.  For an explicit tensor the two are identical. *)
        resultTensor = ArrayTranspose[resultTensor,
            FindPermutation @ Join[
                Range[n1Out],
                Range[n1Out + n1In + 1, n1Out + n1In + n2Out],
                Range[n1Out + 1, n1Out + n1In],
                Range[n1Out + n1In + n2Out + 1, n1Out + n1In + n2Out + n2In]
            ]
        ]
    ];

    resultOutOrder = Join[q1OutOrder, q2OutOrder[[keepQ2OutIdx]]];
    resultInOrder = Join[q1InOrder[[keepQ1InIdx]], q2InOrder];
    resultOutDecompose = Join[s1["Output"]["Decompose"], s2["Output"]["Decompose"][[keepQ2OutIdx]]];
    resultInDecompose = Join[s1["Input"]["Decompose"][[keepQ1InIdx]], s2["Input"]["Decompose"]];

    QuantumOperator[
        QuantumState[
            ArrayVector @ resultTensor,
            QuantumBasis[
                "Output" -> If[resultOutDecompose === {}, QuditBasis[], QuantumTensorProduct @@ resultOutDecompose],
                "Input"  -> If[resultInDecompose === {}, QuditBasis[], QuantumTensorProduct @@ resultInDecompose],
                "Label" -> qo1["Label"] @* qo2["Label"],
                "Picture" -> qo1["Picture"],
                "ParameterSpec" -> MergeParameterSpecs[qo1, qo2]
            ]
        ],
        {resultOutOrder, resultInOrder}
    ]
]


(qo1_QuantumOperator ? QuantumOperatorQ)[qo2_ ? QuantumOperatorQ] := QuantumOperator @ QuantumCircuitOperator[{qo2, qo1}]

orderDuplicates[xs_List] := Block[{next = Function[{ys, y}, If[MemberQ[ys, y], next[ys, y + 1], y]]}, Fold[Append[#1, next[#1, #2]] &, {}, xs]]


(* the measurement's "POVM" form carries its dilation together with the targets that go with it *)
(qo_QuantumOperator ? QuantumOperatorQ)[qmo_ ? QuantumMeasurementOperatorQ] := With[{povm = qmo["POVM"]}, {op = qo @ povm["QuantumOperator"]},
    If[ContainsAll[Select[op["OutputOrder"], Positive], Select[qmo["OutputOrder"], Positive]] && op["VectorQ"], QuantumMeasurementOperator[op, povm["Targets"]], op]
]

(qo_QuantumOperator ? QuantumOperatorQ)[qc_ ? QuantumChannelQ] := With[{op = qo @ qc["QuantumOperator"]},
    If[ContainsAll[Select[op["OutputOrder"], Positive], Select[qc["OutputOrder"], Positive]] && op["VectorQ"], QuantumChannel[op], op]
]

(qo_QuantumOperator ? QuantumOperatorQ)[qm_ ? QuantumMeasurementQ] :=
    If[QuantumMeasurementOperatorQ[#], QuantumMeasurement[#["Sort"]], #] & @ qo[qm["QuantumOperator"]]

(qo_QuantumOperator ? QuantumOperatorQ)[qco_QuantumCircuitOperator ? QuantumCircuitOperatorQ] :=
    QuantumCircuitOperator[Append[qco["Operators"], qo]]


expandQuditBasis[qb_QuditBasis, order1_ ? orderQ, order2_ ? orderQ, defaultDim_Integer : 2] := Enclose @ (
    ConfirmAssert[Length[order1] == qb["Qudits"]];
    QuantumTensorProduct[order2 /. Append[Thread[order1 -> qb["Decompose"]], _Integer -> QuditBasis[defaultDim]]]
)


QuantumOperator /: HoldPattern[Plus[ops__QuantumOperator]] /; Length[{ops}] > 1 := Fold[addQuantumOperators, {ops}]

QuantumOperator /: HoldPattern[Plus[x : Except[_QuantumOperator], qo_QuantumOperator]] := With[{op = qo["Sort"]},
    (* Add x on the diagonal by adding x times the identity operator on op's own
       output space, reusing the operator sum.  Building the identity in the
       computational basis and letting the sum reconcile keeps this correct for a
       density-matrix-type operator, whose stored matrix is the d^2 x d^2
       superoperator rather than the d x d name-dimensioned matrix. *)
    With[{id = With[{v = QuantumOperator[IdentityMatrix[Times @@ op["OutputDimensions"]], op["Order"], QuantumBasis[op["OutputDimensions"]]]},
                    If[op["StateType"] === "Matrix", v["ToMatrix"], v]]},
        QuantumOperator[op + x id, "Label" -> If[op["Label"] === None, None, x + op["Label"]]]
    ]
]

(* The amplitudes an operator stores for the matrix mat: the flattened matrix for a
   vector-type state, the doubled layout for a matrix-type one. *)
operatorAmplitudes[mat_, vectorQ_, nameDimensions_, dimension_] := If[ vectorQ,
    Flatten[mat],
    ArrayReshape[
        Transpose[ArrayReshape[mat, Join[nameDimensions, nameDimensions]], 2 <-> 3],
        {dimension, dimension}
    ]
]

matrixOperator[op_QuantumOperator, mat_, opts___] := QuantumOperator[
    QuantumState[
        operatorAmplitudes[mat, op["VectorQ"], op["MatrixNameDimensions"], op["Dimension"]],
        op["Basis"]
    ],
    op["Order"],
    opts
]

(* A scalar base with a square-operator exponent is the matrix exponential of the
   stored matrix, base^op = MatrixExp[Log[base] op], the same stored-basis reading as
   MatrixExp[op], which stays regular at parameter values where eigenvalues collide.
   A zero base has no logarithm, and 0^op is the limit of base^op as base -> 0
   (zeroBasePower): for a number operator n, 0^n is the vacuum projector, as in the
   zero-temperature thermal state and the total-loss channel. A base that is a
   declared parameter reaches the same limit when its value is zero on the lazy route
   of matrixMapOperator, where the base is read only when the values arrive;
   parameters that are not plain symbols substitute into the closed form instead. The
   general rule further down reads Power arguments matrix-first (MatrixPower) and
   would compute op^base. The base is never an array container or a quantum object.
   An inexact base is zero when it equals 0, an exact number when zeroBaseQ finds it
   zero (without the conventions that read a tiny nonzero number as zero), and a
   symbolic one when PossibleZeroQ does. *)
scalarPowerBaseQ[base_] :=
    FreeQ[base,
        _ ? ArrayContainerQ | _ ? QuantumFrameworkOperatorQ |
        _QuantumState | _QuantumMeasurement | _QuantumBasis | _QuditBasis | _QuditName |
        _Failure | _String | _List
    ]

scalarBasePower[base_ /; Which[! NumericQ[base], PossibleZeroQ[base], InexactNumberQ[base], base == 0, True, zeroBaseQ[base]], mat_] :=
    zeroBasePower[If[InexactNumberQ[base], N[mat], mat]]

scalarBasePower[base_, mat_] := matrixExponential[Log[base] mat]

QuantumOperator /: Power[base_ ? scalarPowerBaseQ, qo_QuantumOperator] /; TrueQ[qo["SquareQ"]] :=
    matrixMapOperator[scalarBasePower[base, #] &, qo, Power[base, #] &]

(* The generic NumericFunction rule below would match Exp[qo] directly with f = Exp,
   so the built-in rewrite of Exp into Power[E, ...] never runs on an operator
   argument. Delegating to the Power normal form makes the spellings of the
   operator exponential one object: the matrix exponential of the stored matrix
   for a square operator, a loud MatrixPower failure otherwise. *)
QuantumOperator /: Exp[qo_QuantumOperator] := E ^ qo

QuantumOperator /: f_Symbol[left : Except[_QuantumOperator] ..., qo_QuantumOperator, right : Except[_QuantumOperator | OptionsPattern[]] ..., opts : OptionsPattern[]] /; MemberQ[Attributes[f], NumericFunction] :=
    matrixMapOperator[matrixFunction[f, #, {left}, {right}, opts] &, qo, f[left, #, right] &]

QuantumOperator /: MatrixExp[qo_QuantumOperator] := matrixMapOperator[matrixExponential, qo, Exp]

(* The operator g(M) for the matrix M of qo in sorted order, labelled with labelF of
   its label: evaluated now, or kept lazy in the parameters (below). Where g fails,
   the result is the Failure that names the reason. *)
matrixMapOperator[g_, qo_, labelF_] := Enclose[
    With[
        {op = qo["Sort"]},
        {label = If[op["Label"] === None, None, labelF[op["Label"]]]},
        If[ lazyMatrixMapQ[op],
            parametricMatrixMapOperator[g, qo, op, label],
            matrixOperator[op, ConfirmBy[Confirm[g[op["Matrix"]]], matrixContainerQ, Defer[label]], "Label" -> label]
        ]
    ],
    innermostFailure
]

(* A failed Confirm returns a Failure wrapping the Failure it confirmed; the innermost
   one names the reason. *)
innermostFailure[failure_] := Replace[failure["Expression"], {inner_ ? FailureQ :> innermostFailure[inner], _ :> failure}]

(* A matrix map g of an operator with declared parameters (f of the matrix for a
   NumericFunction f, the matrix exponential for Exp, MatrixExp and base^op) keeps g
   unapplied until values arrive: the amplitudes are a lazy Function of the
   parameters whose body applies g to the matrix with the values in place.
   Substituting into a closed form of g(M) instead divides by differences of
   eigenvalues, so it is Indeterminate where the values make two of them collide
   (sin M = sin(r)/r M at r = 0) and loses accuracy near such values; applying g to
   the substituted matrix is regular there. Reading the amplitudes gives the closed
   form (the Function at its own parameters), and a partial substitution curries
   the Function. A scalar multiple stays lazy; a sum or a product with another
   operator, or the action on a state, is built from the closed form, so values at
   a collision fail there and values near one lose accuracy without a message.

   The parameters can be the variables of a Function only when they are plain
   unprotected symbols: an indexed th[1] or a Subscript is not, and those operators
   keep the closed form, as do non-square ones, which have no matrix function. *)
lazyMatrixMapQ[op_] := op["ParameterArity"] > 0 && TrueQ[op["SquareQ"]] &&
    MatchQ[op["Parameters"], {__Symbol}] && AllTrue[op["Parameters"], ! MemberQ[Attributes[#], Protected] &]

(* Amplitudes built by this route: a Function of the parameters whose body is
   matrixMapAmplitudes. They are read and substituted through that body. *)
lazyMatrixMapAmplitudesQ[amplitudes_] := MatchQ[amplitudes, HoldPattern[Function[_Symbol | {__Symbol}, _matrixMapAmplitudes]]]

(* The amplitudes of g of the matrix, or a Failure when g is not defined there (Log
   at a zero eigenvalue), so that a substitution reaching that value fails. *)
matrixMapAmplitudes[g_, mat_, vectorQ_, nameDimensions_, dimension_] := Enclose[
    operatorAmplitudes[
        ConfirmBy[Confirm[g[ConfirmBy[mat, ! valuelessEntriesQ[#] &]]], MatrixQ],
        vectorQ, nameDimensions, dimension
    ],
    innermostFailure
]

(* The matrix of the amplitudes amps in the parameter-free basis qb, in the sorted
   order the matrix map works in. *)
amplitudesMatrix[amps_, qb_, order_] := Normal[QuantumOperator[QuantumState[amps, qb], order]["Sort"]["Matrix"]]

(* The matrix the body applies g to, held so that the parameters in it stay free.
   For an operator qo that is itself such a lazy result (g1 of an operator, then
   g2), it is the inner body read as a matrix, so the values reach the innermost
   matrix first rather than the closed form of g1. It is read from qo as given,
   before sorting its order, since sorting rebuilds the amplitudes from the closed
   form; the sorting happens inside the body instead. Otherwise it is the sorted
   matrix: for a diagonal matrix the List of its diagonal entries, made into a sparse
   diagonal matrix only after the substitution, and for any other the matrix as a
   List. A SparseArray is atomic, so neither Function application nor ReplaceAll
   would reach the parameters inside one. The basis written into the body carries no
   parameter specification, which a substitution would otherwise overwrite. *)
heldOperatorMatrix[qo_, op_] /; lazyMatrixMapAmplitudesQ[qo["State"]["State"]] := With[{
    qb = QuantumBasis[qo["Basis"], "ParameterSpec" -> {}],
    order = qo["Order"]
},
    Replace[Extract[qo["State"]["State"], {2}, Hold], Hold[body_] :> Hold[amplitudesMatrix[body, qb, order]]]
]

heldOperatorMatrix[_, op_] := heldMatrix[op["Matrix"]]

heldMatrix[mat_ ? diagonalMatrixQ] := With[{d = Normal[Diagonal[mat]]}, Hold[DiagonalMatrix[d, TargetStructure -> "Sparse"]]]
heldMatrix[mat_] := With[{m = Normal[mat]}, Hold[m]]

(* Function binds its first argument, so With cannot write the parameters there;
   Apply puts them in front of the held body instead. The shape is declared so that
   Wolfram`Arrays` does not probe the body; Arrays keeps that declaration, keyed on
   the Function, for the rest of the session. Evaluating the closed form here fails
   construction when g is undefined for every value of the parameters, and caches
   it for later reads. *)
parametricMatrixMapOperator[g_, qo_, op_, label_] := Enclose @ With[{
    vectorQ = op["VectorQ"],
    nameDimensions = op["MatrixNameDimensions"],
    dimension = op["Dimension"]
},
    With[{
        amplitudes = Function @@ Prepend[
            Replace[heldOperatorMatrix[qo, op], Hold[mat_] :> Hold[matrixMapAmplitudes[g, mat, vectorQ, nameDimensions, dimension]]],
            op["Parameters"]
        ]
    },
        ConfirmBy[lazyClosedForm[amplitudes], ArrayQ];
        ArrayDeclareShape[amplitudes, ArrayDimensions[op["State"]["State"]]];
        QuantumOperator[QuantumState[amplitudes, op["Basis"]], op["Order"], "Label" -> label]
    ]
]

(* The exponential e^M of the operator acting on the state, as Exp[qo][qs] gives: a
   pure state is multiplied by e^M, a mixed state rho becomes e^M rho e^(M^dagger),
   and a superoperator acts on the density vector. An operator on any other qudits
   than exactly the state's, with their dimensions, goes through Exp[qo][qs], which
   extends the operator by the identity on qudits it does not act on, extends the
   state to qudits it lacks, and fails on a dimension mismatch. *)
QuantumOperator /: MatrixExp[qo_QuantumOperator, qs_QuantumState] /; ! wholeRegisterQ[qo, qs] := Exp[qo][qs]

QuantumOperator /: MatrixExp[qo_QuantumOperator, qs_QuantumState] := Enclose @ With[{op = qo["Sort"]},
    QuantumState[
        ConfirmBy[exponentialAction[op, QuantumState[qs, QuantumBasis[op["Input"], qs["Input"]]]], ArrayQ],
        QuantumBasis[
            op["Output"],
            "Label" -> If[op["Label"] === None || qs["Label"] === None, None, Exp[op["Label"]][qs["Label"]]]
        ]
    ]
]

wholeRegisterQ[qo_, qs_] := Sort[Thread[{qo["InputOrder"], qo["InputDimensions"]}]] === Thread[{Range[qs["OutputQudits"]], qs["OutputDimensions"]}]

exponentialAction[op_ ? (#["VectorQ"] &), qs_ ? (#["VectorQ"] &)] := matrixExponential[op["Matrix"], qs["StateVector"]]
exponentialAction[op_ ? (#["VectorQ"] &), qs_] := With[{u = matrixExponential[op["Matrix"]]}, u . qs["DensityMatrix"] . ConjugateTranspose[u]]
exponentialAction[op_, qs_] := ReshapeArray[matrixExponential[op["ToMatrix"]["Matrix"], qs["DensityVector"]], {#, #} & @ op["OutputDimension"]]


(* operators padded with identities to the qudits any of them acts on, each qudit in the
   basis of the last operator acting on it *)
padQuantumOperators[ops : {__QuantumOperator}] := Module[{
    orderInput, orderOutput
},
    orderInput = With[{
            order = Union @@ Through[ops["InputOrder"]],
            qbMap = Association[Thread[#["InputOrder"] -> #["Input"]["Decompose"]] & /@ ops]
        },
            {"OrderedInput", order, QuantumTensorProduct[order /. qbMap]}
    ];
    orderOutput = With[{
            order = Union @@ Through[ops["OutputOrder"]],
            qbMap = Association[Thread[#["OutputOrder"] -> #["Output"]["Decompose"]] & /@ ops]
        },
            {"OrderedOutput", order, QuantumTensorProduct[order /. qbMap]}
    ];
    ((# @@ orderInput)["Sort"] @@ orderOutput)["Sort"] & /@ ops
]

(* {Hamiltonian, jump operators} with each jump operator on some of the qudits acting as
   the identity on the rest: on the Hamiltonian's own order when it covers them all, else
   on the union of every order. Without a Hamiltonian the jump operators share their union. *)
padJumpOperators[h_ ? QuantumOperatorQ, ls : {___ ? QuantumOperatorQ}] := Which[
    AllTrue[ls, #["Order"] === h["Order"] &],
    {h, ls},
    AllTrue[ls, SubsetQ[h["OutputOrder"], #["OutputOrder"]] && SubsetQ[h["InputOrder"], #["InputOrder"]] &],
    {h, If[ #["Order"] === h["Order"],
            #,
            #["OrderedInput", h["InputOrder"], h["Input"]]["OrderedOutput", h["OutputOrder"], h["Output"]]
        ] & /@ ls},
    True,
    With[{ops = padQuantumOperators[Prepend[ls, h]]}, {First[ops], Rest[ops]}]
]

padJumpOperators[None, ls : {___ ? QuantumOperatorQ}] := {None, If[SameQ @@ Through[ls["Order"]], ls, padQuantumOperators[ls]]}

padJumpOperators[h_, ls_] := {h, ls}

(* rates for a list of jump operators: none (every rate is 1), one rate per jump operator,
   or a square matrix of them (a Kossakowski matrix) *)
jumpRatesQ[ls_List, gammas_] :=
    gammas === {} ||
    VectorQ[gammas] && Length[gammas] == Length[ls] ||
    MatrixQ[gammas] && Dimensions[gammas] === {Length[ls], Length[ls]}

addQuantumOperators[qo1_QuantumOperator ? QuantumOperatorQ, qo2_QuantumOperator ? QuantumOperatorQ] := Enclose @ Module[{
    ordered1, ordered2
},
    {ordered1, ordered2} = padQuantumOperators[{qo1, qo2}];
    ConfirmAssert[ordered1["Dimensions"] == ordered2["Dimensions"]];
    QuantumOperator[
        QuantumState[
            addQuantumStates[ordered1["State"], ordered2["State"]],
            "Label" -> If[ordered1["Label"] === None || ordered2["Label"] === None, None, ordered1["Label"] + ordered2["Label"]],
            "ParameterSpec" -> MergeParameterSpecs[ordered1, ordered2]
        ],
        ordered1["Order"]
    ]
]


(* differentiation *)

QuantumOperator /: D[op : _QuantumOperator, args___] := QuantumOperator[D[op["State"], args], op["Order"]]


(* dagger *)

SuperDagger[qo_QuantumOperator] ^:= qo["Dagger"]

SuperStar[qo_QuantumOperator] ^:= qo["Conjugate"]

Transpose[qo_QuantumOperator, args___] ^:= qo["Transpose", args]

Inverse[qo_QuantumOperator] ^:= qo ^ -1


(* Trace *)

QuantumOperator /: Tr[qo_QuantumOperator] := Tr @ qo["Matrix"]


(* commutator *)

QuantumOperator /: Commutator[a_QuantumOperator, b_QuantumOperator] := a@b - b@a


(* simplify *)

Scan[
    (Symbol[#][qo_QuantumOperator, args___] ^:= qo[#, args]) &,
    {"Simplify", "FullSimplify", "Chop", "ComplexExpand"}
]


(* join *)

QuantumOperator[qo_ ? QuantumOperatorQ] := qo

QuantumOperator[qo__QuantumOperator ? QuantumOperatorQ] := QuantumOperator["Multiplexer"[qo]]


(* equality *)

QuantumOperator /: Equal[qo__QuantumOperator] :=
    Equal @@ (#["Picture"] & /@ {qo}) && And @@ Thread[Equal @@ (Chop @ SetPrecisionNumeric @ ArrayVector @ #["Sort"]["MatrixRepresentation"] & /@
        If[Or @@ Through[{qo}["MatrixQ"]], Through[{qo}["ToMatrix"]], {qo}])]

QuantumOperator /: Unequal[qo__QuantumOperator] := ! Equal[qo]


(* conversion *)

QuantumOperator[obj : _QuantumMeasurementOperator | _QuantumMeasurement | _QuantumChannel | _QuantumCircuitOperator, opts___] :=
    QuantumOperator[obj["QuantumOperator"], opts]


(* parameterization *)

(qo_QuantumOperator ? QuantumOperatorQ)[opts___] := QuantumCircuitOperator[qo][opts]

(qo_QuantumOperator ? QuantumOperatorQ)[ps : PatternSequence[p : Except[_Association], ___]] /; ! MemberQ[QuantumOperator["Properties"], p] && Length[{ps}] <= qo["ParameterArity"] :=
    qo[AssociationThread[Take[qo["Parameters"], UpTo[Length[{ps}]]], {ps}]]

(qo_QuantumOperator ? QuantumOperatorQ)[rules_ ? AssociationQ] /; ContainsOnly[Keys[rules], qo["Parameters"]] :=
    Enclose[QuantumOperator[Confirm[qo["State"][rules]], qo["Order"]], innermostFailure]


(* *)

StackQuantumOperators[ops : {_ ? QuantumOperatorQ ..}, name_ : "\[ScriptCapitalE]"] := Block[{
    basis = First[ops]["Basis"],
    order = MapAt[Prepend[#, Min[#] - 1] &, First[ops]["Order"], {1}]
},
    basis = QuantumBasis[basis,
        "Output" -> QuantumTensorProduct[QuditBasis[Subscript[name, #] & /@ Range @ Length @ ops], basis["Output"]],
        "Input" -> basis["Input"]
    ];
    QuantumOperator[
        QuantumOperator[
            SparseArray[#["MatrixRepresentation"] & /@ ops],
            QuantumBasis[basis["OutputDimensions"], basis["InputDimensions"]]
        ],
        order,
        basis
    ]
]
