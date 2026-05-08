(* ::Package:: *)

Package["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`QuantumFramework`Gates`"]


PackageExport["QuantumOperator"]

PackageScope["QuantumOperatorQ"]
PackageScope["StackQuantumOperators"]


QuantumOperator::invalidInputOrder = "input order should be a list of distinct input qudit positions"
QuantumOperator::invalidOutputOrder = "output order should be a list of distinct output qudit positions"
QuantumOperator::invalidName = "`1` is not a recognized QuantumOperator constructor"

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


QuantumOperator[matrix_ ? MatrixQ, order : _ ? autoOrderQ, args___, opts : OptionsPattern[]] := Enclose @ Module[{
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

QuantumOperator[matrix_ ? MatrixQ, opts : OptionsPattern[]] := QuantumOperator[matrix, QuantumBasis[primeFactors[#1], primeFactors[#2]] & @@ Dimensions[matrix], opts]

QuantumOperator[matrix_ ? MatrixQ, args__, opts : OptionsPattern[]] := Block[{
    outMultiplicity, inMultiplicity, result
},
    result = Enclose @ Module[{newMatrix = matrix, outputs, inputs,
        basis, newOutputQuditBasis, newInputQuditBasis, state},
        {outputs, inputs} = Dimensions[newMatrix];
        basis = ConfirmBy[QuantumBasis[args], QuantumBasisQ, "Invalid basis"];
        If[ basis["Dimension"] =!= outputs * inputs,
            If[ basis["InputDimension"] == 1,
                basis = QuantumBasis[basis, "Input" -> basis["Output"]["Dual"]]
            ];
            outMultiplicity = Quiet @ Log[basis["OutputDimension"], outputs];
            inMultiplicity = Quiet @ Log[basis["InputDimension"], inputs];
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
                Flatten[newMatrix],
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
    QuantumOperator[
        {qo, LCM[Quotient[Length[outputOrder], outputQudits], Quotient[Length[inputOrder], inputQudits]]},
        {outputOrder, inputOrder},
        opts
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

(qo1_QuantumOperator ? QuantumOperatorQ)[qo2_QuantumOperator ? QuantumOperatorQ] /; (
    (qo1["MatrixQ"] || qo2["MatrixQ"]) && qo1["Picture"] === qo2["Picture"] &&
    Intersection[qo1["OutputOrder"], Complement[qo2["OutputOrder"], qo1["InputOrder"]]] === {} &&
    Intersection[Complement[qo1["InputOrder"], qo2["OutputOrder"]], qo2["InputOrder"]] === {}
) := With[{shift = Max[qo1["Order"], qo2["Order"]] + 1},
    qo1["Bend", shift][qo2["Bend", shift]]["Unbend"]
]


(qo1_QuantumOperator ? QuantumOperatorQ)[qo2_QuantumOperator ? QuantumOperatorQ] /; (
    qo1["VectorQ"] && qo2["VectorQ"] && qo1["Picture"] === qo2["Picture"] &&
    Intersection[qo1["OutputOrder"], Complement[qo2["OutputOrder"], qo1["InputOrder"]]] === {} &&
    Intersection[Complement[qo1["InputOrder"], qo2["OutputOrder"]], qo2["InputOrder"]] === {}
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
    resultTensor = If[Length[contractPairs] == 0,
        TensorProduct[s1["StateTensor"], s2["StateTensor"]],
        TensorContract[TensorProduct[s1["StateTensor"], s2["StateTensor"]], contractPairs]
    ];

    (* After contraction, surviving axes appear in tensor order:
         [q1.outputs] [q1.inputs kept] [q2.outputs kept] [q2.inputs]
       Permute so all outputs come first, then all inputs:
         [q1.outputs] [q2.outputs kept] [q1.inputs kept] [q2.inputs] *)
    With[{
        n1Out = q1OutQ, n1In = Length[keepQ1InIdx],
        n2Out = Length[keepQ2OutIdx], n2In = q2InQ
    },
        resultTensor = Transpose[resultTensor,
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
            SparseArrayFlatten @ resultTensor,
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


(qo_QuantumOperator ? QuantumOperatorQ)[qmo_ ? QuantumMeasurementOperatorQ] := With[{op = qo @ qmo["SuperOperator"]},
    If[ContainsAll[Select[op["OutputOrder"], Positive], Select[qmo["OutputOrder"], Positive]] && op["VectorQ"], QuantumMeasurementOperator[op, qmo["Targets"]], op]
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
    QuantumOperator[
        Plus[DiagonalMatrix[ConstantArray[x, Min[op["MatrixNameDimensions"]], SparseArray], 0, op["MatrixNameDimensions"]], op["Matrix"]],
        op["Order"],
        op["Basis"],
        "Label" -> If[op["Label"] === None, None, x + op["Label"]]
    ]
]

matrixOperator[op_QuantumOperator, mat_, opts___] := QuantumOperator[
    QuantumState[
        If[ op["VectorQ"],
            Flatten[mat],
            ArrayReshape[
                Transpose[ArrayReshape[mat, Join[#, #] & @ op["MatrixNameDimensions"]], 2 <-> 3],
                {#, #} & @ op["Dimension"]
            ]
        ],
        op["Basis"]
    ],
    op["Order"],
    opts
]

QuantumOperator /: f_Symbol[left : Except[_QuantumOperator] ..., qo_QuantumOperator, right : Except[_QuantumOperator | OptionsPattern[]] ..., opts : OptionsPattern[]] /; MemberQ[Attributes[f], NumericFunction] := Enclose @ With[{
    op = qo["Sort"]
},
    matrixOperator[
        op,
        ConfirmBy[matrixFunction[f, op["Matrix"], {left}, {right}, opts], MatrixQ],
        "Label" -> If[op["Label"] === None, None, f[left, op["Label"], right]]
    ]
]

QuantumOperator /: MatrixExp[qo_QuantumOperator] := Enclose @ With[{
    op = qo["Sort"]
},
    matrixOperator[
        op,
        ConfirmBy[MatrixExp[op["Matrix"]], MatrixQ],
        "Label" -> If[op["Label"] === None, None, Exp[op["Label"]]]
    ]
]

QuantumOperator /: MatrixExp[qo_QuantumOperator, qs_QuantumState] := Enclose @ With[{op = qo["Sort"]},
    QuantumState[
        If[ op["VectorQ"] && qs["VectorQ"],
            MatrixExp[op["Matrix"], QuantumState[qs, QuantumBasis[op["Input"], qs["Input"]]]["StateVector"]],
            ArrayReshape[
                MatrixExp[op["ToMatrix"]["Matrix"], QuantumState[qs, QuantumBasis[op["Input"], qs["Input"]]]["DensityVector"]],
                {#, #} & @ op["OutputDimension"]
            ]
        ],
        QuantumBasis[
            op["Output"],
            "Label" -> If[op["Label"] === None || qs["Label"] === None, None, Exp[op["Label"]][qs["Label"]]]
        ]
    ]
]


addQuantumOperators[qo1_QuantumOperator ? QuantumOperatorQ, qo2_QuantumOperator ? QuantumOperatorQ] := Enclose @ Module[{
    orderInput, orderOutput,
    ordered1, ordered2
},
    orderInput = With[{
            order = Union[qo1["InputOrder"], qo2["InputOrder"]],
            qbMap = Association[Thread[qo1["InputOrder"] -> qo1["Input"]["Decompose"]], Thread[qo2["InputOrder"] -> qo2["Input"]["Decompose"]]]
        },
            {"OrderedInput", order, QuantumTensorProduct[order /. qbMap]}
    ];
    orderOutput = With[{
            order = Union[qo1["OutputOrder"], qo2["OutputOrder"]],
            qbMap = Association[Thread[qo1["OutputOrder"] -> qo1["Output"]["Decompose"]], Thread[qo2["OutputOrder"] -> qo2["Output"]["Decompose"]]]
        },
            {"OrderedOutput", order, QuantumTensorProduct[order /. qbMap]}
    ];
    ordered1 = ((qo1 @@ orderInput)["Sort"] @@ orderOutput)["Sort"];
    ordered2 = ((qo2 @@ orderInput)["Sort"] @@ orderOutput)["Sort"];
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
    Equal @@ (#["Picture"] & /@ {qo}) && And @@ Thread[Equal @@ (Chop @ SetPrecisionNumeric @ SparseArrayFlatten @ #["Sort"]["MatrixRepresentation"] & /@
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
    QuantumOperator[qo["State"][rules], qo["Order"]]


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
