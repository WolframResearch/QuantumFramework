Package["Wolfram`QuantumFramework`"]

PackageScope["UnitaryEulerAngles"]
PackageScope["UnitaryEulerAnglesWithPhase"]
PackageScope["TwoQubitKAK"]
PackageScope["TwoQubitKAKSymbolic"]
PackageScope["$MagicBasis"]



$QuantumOperatorProperties = {
    "Order",
    "InputOrder", "OutputOrder", "ControlOrder", "TargetOrder",
    "MatrixRepresentation", "Matrix",
    "TensorRepresentation", "Tensor",
    "Ordered", "OrderedInput", "OrderedOutput",
    "SortInput", "SortOutput", "Sort", "SortedQ",
    "ReverseOutput", "ReverseInput", "Reverse",
    "Shift",
    "OrderedMatrixRepresentation", "OrderedMatrix",
    "OrderedTensorRepresentation", "OrderedTensor",
    "Arity", "MaxArity", "FullArity", "TargetArity",
    "Range", "FullInputOrder", "FullOutputOrder", "FullOrder", "InputQuditOrder", "OutputQuditOrder", "QuditOrder",
    "OutputOrderQuditMapping", "InputOrderQuditMapping",
    "FirstOutputQudit", "LastOutputQudit", "FirstInputQudit", "LastInputQudit", "InputOrderQuditMapping",
    "HermitianQ", "UnitaryQ", "Eigenvalues", "Eigenvectors", "Eigensystem", "Projectors",
    "Transpose", "ConjugateTranspose",
    "UnstackOutput", "UnstackInput",
    "QuantumOperator", "Operator",
    "Computational", "Diagonalize",
    "Dagger", "Dual",
    "TraceNorm",
    "PauliDecompose",
    "KAK", "Decompose",
    "CircuitDiagram",
    "OrderedFormula",
    "Simplify", "FullSimplify", "Chop", "ComplexExpand"
};

QuantumOperator["Properties"] := Union @ Join[$QuantumOperatorProperties, Complement[QuantumState["Properties"], {
    "BlochCartesianCoordinates", "BlochSphericalCoordinates", "BlochPlot"
}]];

QuantumOperatorProp[qo_, "Properties"] := Union @ $QuantumOperatorProperties

QuantumOperatorProp[qo_, "AllProperties"] := Union @ Join[QuantumOperator["Properties"], Complement[qo["State"]["AllProperties"], {
    "BlochCartesianCoordinates", "BlochSphericalCoordinates", "BlochPlot"
}]]


qo_QuantumOperator["ValidQ"] := QuantumOperatorQ[qo]


QuantumOperator::undefprop = "property `` is undefined for this state";

QuantumOperator::failprop = "property `` have failed with ``";


(* basic getters *)

QuantumOperatorProp[QuantumOperator[state_, _], "State"] := state

(* Direct structural getter: the basis lives in the wrapped state. Without this, qo["Basis"]
   falls through all ~124 QuantumOperatorProp clauses to the generic delegation, whose guard
   recomputes AllProperties on both operator and state and intersects them. Because the
   cache-eligibility guard reads qo["Basis"]["ParameterArity"] on every property access, that
   fall-through taxed every read; the direct getter removes it. *)
QuantumOperatorProp[QuantumOperator[state_, _], "Basis"] := state["Basis"]

QuantumOperatorProp[QuantumOperator[_, order : {_, _}], "Order"] := order

QuantumOperatorProp[QuantumOperator[_, {_, inputOrder_}], "InputOrder"] := inputOrder

QuantumOperatorProp[QuantumOperator[_, {outputOrder_, _}], "OutputOrder"] := outputOrder


(qo_QuantumOperator[prop_ ? propQ, args___]) /; QuantumOperatorQ[qo] := With[{
    result = QuantumOperatorProp[qo, prop, args]
    },
    If[ TrueQ[$QuantumFrameworkPropCache] &&
        ! MemberQ[{"Properties", "AllProperties", "State", "Basis", "Order", "InputOrder", "OutputOrder"}, prop] &&
        QuantumOperatorProp[qo, "Basis"]["ParameterArity"] == 0,
        cacheProperty[QuantumOperatorProp[qo, prop, args], result],
        result
    ] /;
        (!FailureQ[Unevaluated @ result] || Message[QuantumOperator::failprop, prop, result]) &&
        (!MatchQ[result, _QuantumOperatorProp] || Message[QuantumOperator::undefprop, prop])
]


(* computed properties *)

QuantumOperatorProp[qo_, "Arity"] := Length @ qo["InputOrder"]

QuantumOperatorProp[qo_, "Range"] := Max[qo["InputOrder"]] - Min[qo["InputOrder"]] + 1

QuantumOperatorProp[qo_, "MaxArity"] := Max[qo["InputQudits"], qo["Arity"]]

QuantumOperatorProp[qo_, "FullArity"] := Max[qo["InputQudits"], qo["Range"]]

QuantumOperatorProp[qo_, "FullInputOrder"] := Which[qo["InputDimension"] > 1,
    Take[
        If[MatchQ[qo["InputOrder"], {___, _ ? NonPositive, ___}], Identity, # - Min[#, 1] + 1 &] @
            If[ Length[qo["InputOrder"]] > 0,
                Join[Complement[Range[Max[qo["InputOrder"]] - qo["InputQudits"] + 1, Max[qo["InputOrder"]]], qo["InputOrder"]], qo["InputOrder"]],
                Range[qo["InputQudits"]]
            ],
        - qo["InputQudits"]
    ],
    qo["InputDimension"] == 0,
    {1},
    True,
    {}
]

QuantumOperatorProp[qo_, "FullOutputOrder"] := Which[qo["OutputDimension"] > 1,
    Take[
        If[ MatchQ[qo["OutputOrder"], {___, _ ? NonPositive, ___}], Identity, # - Min[#, 1] + 1 &] @
            If[ Length[qo["OutputOrder"]] > 0,
                Join[Complement[Range[Max[qo["OutputOrder"]] - qo["OutputQudits"] + 1, Max[qo["OutputOrder"]]], qo["OutputOrder"]], qo["OutputOrder"]],
                Range[qo["OutputQudits"]]
            ],
        - qo["OutputQudits"]
    ],
    qo["OutputDimension"] == 0,
    {1},
    True,
    {}
]

QuantumOperatorProp[qo_, "FullOrder"] := {qo["FullOutputOrder"], qo["FullInputOrder"]}

QuantumOperatorProp[qo_, "SetFullOutputOrder"] := QuantumOperator[qo, {qo["FullOutputOrder"], qo["InputOrder"]}]

QuantumOperatorProp[qo_, "SetFullInputOrder"] := QuantumOperator[qo, {qo["OutputOrder"], qo["FullInputOrder"]}]

QuantumOperatorProp[qo_, "SetFullOrder"] := QuantumOperator[qo, {qo["FullOutputOrder"], qo["FullInputOrder"]}]


QuantumOperatorProp[qo_, "ControlOrder1"] :=
    FirstCase[qo["Label"], Subscript["C", _][control1_, ___] | Interpretation[_, Subscript["C", _][control1_, ___]] :> control1, {}, {0}]

QuantumOperatorProp[qo_, "ControlOrder0"] :=
    FirstCase[qo["Label"], Subscript["C", _][_, control0_]  | Interpretation[_, Subscript["C", _][_, control0_]] :> control0, {}, {0}]

QuantumOperatorProp[qo_, "ControlOrder"] :=
    FirstCase[qo["Label"], Subscript["C", _][control1_, control0_ : {}] | Interpretation[_, Subscript["C", _][control1_, control0_ : {}]] :> Join[control1, control0], {}, {0}]

QuantumOperatorProp[qo_, "TargetOrder"] := Enclose[DeleteCases[qo["InputOrder"], Alternatives @@ Confirm @ qo["ControlOrder"]], qo["InputOrder"] &]

QuantumOperatorProp[qo_, "ControlArity"] := Length @ qo["ControlOrder"]

QuantumOperatorProp[qo_, "TargetArity"] := Length @ qo["TargetOrder"]

QuantumOperatorProp[qo_, "FirstInputQudit"] := Min @ qo["FullInputOrder"]

QuantumOperatorProp[qo_, "LastInputQudit"] := Max @ qo["FullInputOrder"]

QuantumOperatorProp[qo_, "FirstOutputQudit"] := Min @ qo["FullOutputOrder"]

QuantumOperatorProp[qo_, "LastOutputQudit"] := Max @ qo["FullOutputOrder"]

QuantumOperatorProp[qo_, "InputQuditOrder"] := qo["InputOrder"] - Min[qo["FullInputOrder"]] + 1

QuantumOperatorProp[qo_, "OutputQuditOrder"] := qo["OutputOrder"] - Min[qo["FullOutputOrder"]] + 1

QuantumOperatorProp[qo_, "QuditOrder"] := {qo["OutputQuditOrder"], qo["InputQuditOrder"]}

QuantumOperatorProp[qo_, "InputOrderQuditMapping"] := Thread[qo["FullInputOrder"] -> Range[qo["InputQudits"]]]

QuantumOperatorProp[qo_, "OutputOrderQuditMapping"] := Thread[qo["FullOutputOrder"] -> Range[qo["OutputQudits"]]]

QuantumOperatorProp[qo_, "OutputOrderDimensions"] := AssociationMap[
    qo["OutputDimensions"][[Replace[#, qo["OutputOrderQuditMapping"]]]] &,
    qo["FullOutputOrder"]
]

QuantumOperatorProp[qo_, "InputOrderDimensions"] := AssociationMap[
    qo["InputDimensions"][[Replace[#, qo["InputOrderQuditMapping"]]]] &,
    qo["FullInputOrder"]
]

QuantumOperatorProp[qo_, "TargetDimensions"] := qo["TargetOrder"] /. qo["InputOrderDimensions"]


QuantumOperatorProp[qo_, "SquareQ"] := qo["OutputDimension"] == qo["InputDimension"]

QuantumOperatorProp[qo_, "Tensor"] := qo["Sort"]["StateTensor"]

QuantumOperatorProp[qo_, "TensorRepresentation"] := qo["Sort"]["State"]["TensorRepresentation"]


QuantumOperatorProp[qo_, "Matrix"] := qo["Sort"]["StateMatrix"]

QuantumOperatorProp[qo_, "MatrixRepresentation"] := qo["Computational"]["Matrix"]

QuantumOperatorProp[qo_, "Operator"] := qo["Amplitudes"]

QuantumOperatorProp[qo_, "Formula"] := qo["Sort"]["State"]["Formula"]


QuantumOperatorProp[qo_, "QuantumOperator"] := qo

QuantumOperatorProp[qo_, name : "Computational" | "SchmidtBasis" | "SpectralBasis"] := QuantumOperator[qo["State"][name], qo["Order"]]


QuantumOperatorProp[qo_, "OrderedMatrix"] := qo["Ordered"]["Matrix"]

QuantumOperatorProp[qo_, "OrderedMatrixRepresentation"] := qo["Ordered"]["MatrixRepresentation"]

QuantumOperatorProp[qo_, "OrderedTensor"] := qo["Ordered"]["Tensor"]

QuantumOperatorProp[qo_, "OrderedTensorRepresentation"] := qo["Ordered"]["TensorRepresentation"]


QuantumOperatorProp[qo_, "MatrixQuantumState"] /; qo["OutputDimension"] == qo["InputDimension"] := QuantumState[qo["Matrix"], qo["Output"]]

QuantumOperatorProp[qo_, prop : "MatrixState" | "VectorState" | "ToMatrix" | "ToVector" | "Vector"] := QuantumOperator[qo["State"][prop], qo["Order"]]

QuantumOperatorProp[qo_, "PermuteInput", perm_Cycles] := QuantumOperator[
    qo["State"]["PermuteInput", perm],
    qo["Order"]
]

QuantumOperatorProp[qo_, "PermuteOutput", perm_Cycles] := QuantumOperator[
    qo["State"]["PermuteOutput", perm],
    qo["Order"]
]

QuantumOperatorProp[qo_, "OrderInputExtra", inputSource_ ? orderQ, inputTarget_ ? orderQ] := Module[{
    extra = DeleteCases[inputTarget, Alternatives @@ inputSource],
    inputPerm, outputTarget
},
    (* assume extra qudits are on the right *)
    inputPerm = FindPermutation[Join[inputSource, extra], inputTarget];
    outputTarget = qo["FullOutputOrder"] /. Thread[Permute[qo["FullInputOrder"], inputPerm] -> inputTarget];
    QuantumOperator[
        qo["PermuteInput", inputPerm],
        {outputTarget, Sort @ inputTarget}
    ]
]

QuantumOperatorProp[qo_, "OrderOutputExtra", outputSource_ ? orderQ, outputTarget_ ? orderQ] := Module[{
    extra = DeleteCases[outputTarget, Alternatives @@ outputSource],
    outputPerm, inputTarget
},
    (* assume extra qudits are on the right *)
    outputPerm = FindPermutation[Join[outputSource, extra], outputTarget];
    inputTarget = qo["FullInputOrder"] /. Thread[Permute[qo["FullOutputOrder"], outputPerm] -> outputTarget];
    QuantumOperator[
        qo["PermuteOutput", outputPerm],
        {Sort @ outputTarget, inputTarget}
    ]
]

QuantumOperatorProp[qo_, "SortInput"] := If[
    OrderedQ[qo["InputOrder"]],
    qo,
    QuantumOperator[
        qo[
            "PermuteInput",
            InversePermutation @ FindPermutation[Ordering @ Ordering @ qo["FullInputOrder"]]
        ]["State"],
        {qo["OutputOrder"], Sort @ qo["InputOrder"]}
    ]
]

QuantumOperatorProp[qo_, "SortOutput"] := If[
    OrderedQ[qo["OutputOrder"]],
    qo,
    QuantumOperator[
        qo[
            "PermuteOutput",
            InversePermutation @ FindPermutation[Ordering @ Ordering @ qo["FullOutputOrder"]]
        ]["State"],
        {Sort @ qo["OutputOrder"], qo["InputOrder"]}
    ]
]


QuantumOperatorProp[qo_, "Sort"] := QuantumOperator[
    qo["SortOutput"]["SortInput"],
    "Label" -> collectLabel @ Replace[qo["Label"], {
        Subscript["C", subLabel_][controls__] :>
            Subscript["C", sortLabel[subLabel, qo["TargetOrder"]]][controls],
        Subscript["R", subLabel_][angle_] :>
            Subscript["R", sortLabel[subLabel, qo["InputOrder"]]][angle],
        "\[Pi]"[perm__] :> "\[Pi]" @@ {perm}[[Ordering[qo["InputOrder"]]]],
        subLabel_  :> sortLabel[subLabel, qo["InputOrder"]]
    }]
]

QuantumOperatorProp[qo_, "SortedQ"] := AllTrue[qo["Order"], OrderedQ]


QuantumOperatorProp[qo_, "ReverseOutput"] := QuantumOperator[qo["State"], {Reverse @ qo["OutputOrder"], qo["InputOrder"]}]

QuantumOperatorProp[qo_, "ReverseInput"] := QuantumOperator[qo["State"], {qo["OutputOrder"], Reverse @ qo["InputOrder"]}]

QuantumOperatorProp[qo_, "Reverse"] := QuantumOperator[qo["State"], Reverse /@ qo["Order"]]


QuantumOperatorProp[qo_, "Ordered" | "OrderedInput"] := qo["OrderedInput", Sort @ qo["FullInputOrder"]]

QuantumOperatorProp[qo_, "OrderedOutput"] := qo["OrderedOutput", Sort @ qo["FullOutputOrder"]]

QuantumOperatorProp[qo_, "OrderedOutput", from_Integer, to_Integer] /; qo["OutputQudits"] == 1 :=
    qo["OrderedOutput", Range[from, to], QuditBasis[qo["Output"], to - from + 1]]

QuantumOperatorProp[qo_, "Ordered" | "OrderedInput", from_Integer, to_Integer] /; qo["InputQudits"] == 1 :=
    qo["OrderedInput", Range[from, to], QuditBasis[qo["Input"], to - from + 1]]

QuantumOperatorProp[qo_, prop : "Ordered" | "OrderedInput" | "OrderedOutput", from_Integer, to_Integer] :=
    qo[prop, Range[from, to]]

QuantumOperatorProp[qo_, "Ordered" | "OrderedInput", order_ ? orderQ] :=
    qo["OrderedInput", order, QuantumTensorProduct[
        order /. Join[
            Thread[qo["FullInputOrder"] -> qo["Input"]["Decompose"]],
            Thread[qo["FullOutputOrder"] -> qo["Output"]["Decompose"]],
            # -> QuditBasis[2] & /@ Complement[order, Join @@ qo["FullOrder"]]
        ]
    ]]

QuantumOperatorProp[qo_, "OrderedOutput", order_ ? orderQ] :=
    qo["OrderedOutput", order, QuantumTensorProduct[
            order /. Join[
            Thread[qo["FullOutputOrder"] -> qo["Output"]["Decompose"]],
            Thread[qo["FullInputOrder"] -> qo["Input"]["Decompose"]],
            # -> QuditBasis[2] & /@ Complement[order, Join @@ qo["FullOrder"]]
        ]
    ]]

QuantumOperatorProp[qo_, "Ordered", from_Integer, to_Integer, qb_ ? QuditBasisQ] := qo["Ordered", Range[from, to], qb]

QuantumOperatorProp[qo_, "Ordered", order_ ? orderQ, qb_ ? QuditBasisQ] :=
    qo["OrderedInput", order, qb]

QuantumOperatorProp[qo_, "Ordered" | "OrderedInput" | "OrderedOutput", {}, ___] := qo

QuantumOperatorProp[qo_, "OrderedInput", qb_ ? QuditBasisQ] := qo["OrderedInput", qo["FullInputOrder"], qb]

QuantumOperatorProp[qo_, "OrderedInput", order_ ? orderQ, qb_ ? QuditBasisQ] := Enclose @ Block[{
    arity, pos
},
    ConfirmAssert[ContainsAll[order, qo["FullInputOrder"]], "Given order should contain all operator order qudits"];
    arity = Length[order];
    pos = Catenate @ Lookup[PositionIndex[order], qo["FullInputOrder"], Nothing];
    ConfirmAssert[arity <= qb["Qudits"], "Order size should be less than or equal to number of qudits"];
    QuantumOperator[
        If[ arity > qo["InputQudits"],
            QuantumTensorProduct[
                QuantumOperator[qo, {qo["OutputOrder"], qo["FullInputOrder"]}, "Input" -> qb["Extract", pos]],
                With[{iqb = qb["Delete", pos]["Dual"]},
                    QuantumOperator["Identity"[iqb], Max[qo["LastOutputQudit"], qo["LastInputQudit"]] + Range @ iqb["Qudits"]]
                ]
            ],
            QuantumOperator[qo, "Input" -> qb["Extract", pos]]
        ]["OrderInputExtra", qo["FullInputOrder"], order],
        "Label" -> qo["Label"]
    ]
]

QuantumOperatorProp[qo_, "OrderedOutput", qb_ ? QuditBasisQ] := qo["OrderedOutput", qo["FullOutputOrder"], qb]

QuantumOperatorProp[qo_, "OrderedOutput", order_ ? orderQ, qb_ ? QuditBasisQ] := Enclose @ Block[{
    arity, pos
},
    ConfirmAssert[ContainsAll[order, qo["FullOutputOrder"]], "Given order should contain all operator order qudits"];
    arity = Length[order];
    pos = Catenate @ Lookup[PositionIndex[order], qo["FullOutputOrder"], Nothing];
    ConfirmAssert[arity <= qb["Qudits"], "Order size should be less than or equal to number of qudits"];
    QuantumOperator[
        If[ arity > qo["OutputQudits"],
            QuantumTensorProduct[
                QuantumOperator[qo, {qo["FullOutputOrder"], qo["InputOrder"]}, "Output" -> qb["Extract", pos]],
                With[{iqb = qb["Delete", pos]},
                    QuantumOperator["Identity"[iqb], Max[qo["LastOutputQudit"], qo["LastInputQudit"]] + Range @ iqb["Qudits"]]
                ]
            ],
            QuantumOperator[qo, "Output" -> qb["Extract", pos]]
        ]["OrderOutputExtra", qo["FullOutputOrder"], order],
        "Label" -> qo["Label"]
    ]
]

QuantumOperatorProp[qo_, "OrderedFormula", OptionsPattern[]] /; qo["State"]["DegenerateStateQ"] := 0

QuantumOperatorProp[qo_, "OrderedFormula", OptionsPattern[]] /; qo["State"]["EmptyStateQ"] := ""

QuantumOperatorProp[qo_, "OrderedFormula", OptionsPattern["Normalize" -> False]] := With[{s = qo["State"]["Pure"]},
    With[{
        v = SparseArray @ s[If[TrueQ[OptionValue["Normalize"]], "NormalizedStateVector", "StateVector"]],
        d = s["InputDimension"],
        order = Join @@ qo["Order"]
    },
        Dot[
            Map[
                Replace[qn_QuditName :> With[{names = qn["Name"]},
                    With[{qudits = #["Qudits"] & /@ names},
                        QuditName @@ MapThread[
                            QuditName[Thread[Subscript[Flatten[{#1["Name"]}], OverHat /@ #2]], "Dual" -> #1["DualQ"]] &,
                            {names, TakeList[order, qudits]}
                        ] /; Total[qudits] == Length[order]
                    ] /; MatchQ[names, {__QuditName}]
                ]],
                With[{pos = Catenate @ v["ExplicitPositions"]}, s["Names", Thread[{Quotient[pos - 1, d] + 1, Mod[pos - 1, d] + 1}]]]
            ],
            v["ExplicitValues"]
        ]
    ]
]

QuantumOperatorProp[qo_, "Reorder", order : {_ ? orderQ | Automatic, _ ? orderQ | Automatic}, controlQ_ : False] := Block[{
    output = Replace[order[[1]], Automatic :> qo["FullOutputOrder"]], input = Replace[order[[2]], Automatic :> qo["FullInputOrder"]],
    inputRepl, outputRepl
},
    inputRepl =
        Thread[
            Take[
                If[ TrueQ[controlQ],
                    Join[#, Complement[qo["FullInputOrder"], #]] & @ Join[Sort @ qo["ControlOrder"], qo["TargetOrder"]],
                    Take[qo["FullInputOrder"], UpTo[Length[qo["FullInputOrder"]]]]
                ],
                UpTo[Length[input]]
            ] ->
            Take[input, UpTo[Length[qo["FullInputOrder"]]]]
        ];
    outputRepl =
        Thread[
            Take[qo["FullOutputOrder"], UpTo[Length[output]]] ->
            Take[output, UpTo[Length[qo["FullOutputOrder"]]]]];
    QuantumOperator[
        qo["State"],
        {qo["FullOutputOrder"] /. outputRepl, qo["FullInputOrder"] /. inputRepl},
        "Label" -> ReplaceAll[qo["Label"],
            Subscript["C", name_][c1_, c0_] :> RuleCondition[Subscript["C", name][c1 /. inputRepl, c0 /. inputRepl]]
        ]
    ]
]

QuantumOperatorProp[qo_, "Reorder", order_ ? orderQ, args___] := qo["Reorder", {order, Automatic}, args]


QuantumOperatorProp[qo_, "Shift", n : _Integer : 1] := qo["Reorder", qo["Order"] /. k_Integer :> k + n]



QuantumOperatorProp[qo_, "UnstackOutput", n_Integer : 1] /; 1 <= n <= qo["OutputQudits"] :=
    QuantumOperator[#, {Drop[qo["OutputOrder"], {n}], qo["InputOrder"]}] & /@ qo["State"]["UnstackOutput", n]

QuantumOperatorProp[qo_, "UnstackInput", n_Integer : 1] /; 1 <= n <= qo["InputQudits"] :=
    QuantumOperator[#, {qo["OutputOrder"], Drop[qo["InputOrder"], {n}]}] & /@ qo["State"]["UnstackInput", n]


QuantumOperatorProp[qo_, "HermitianQ"] := HermitianMatrixQ[qo["Matrix"]]

QuantumOperatorProp[qo_, "UnitaryQ"] := UnitaryMatrixQ[N[qo["Matrix"]]]

QuantumOperatorProp[qo_, "Eigenvalues", opts___] /; qo["SquareQ"] := eigenvalues[qo["MatrixRepresentation"], opts, "Sort" -> False, "Normalize" -> True]

QuantumOperatorProp[qo_, "Eigenvectors", opts___] /; qo["SquareQ"] := eigenvectors[qo["MatrixRepresentation"], opts, "Sort" -> False, "Normalize" -> True]

QuantumOperatorProp[qo_, "Eigensystem", opts___] /; qo["SquareQ"] := eigensystem[qo["MatrixRepresentation"], opts, "Sort" -> False, "Normalize" -> True]

QuantumOperatorProp[qo_, "Eigenbasis", opts___] /; qo["SquareQ"] := With[{vecs = qo["Eigenvectors", opts, "Sort" -> False]},
    QuantumBasis[AssociationThread[Symbol["\[FormalV]" <> ToString[#]] & /@ Range[Length[vectors]], vectors]]
]

QuantumOperatorProp[qo_, "SplitBasis"] := With[{
    output = Thread[{qo["FullOutputOrder"], qo["Output"]["Decompose"]}],
    input = Thread[{qo["FullInputOrder"], qo["Input"]["Decompose"]}]
},
    {
        QuantumBasis[
            "Output" -> QuantumTensorProduct @ Select[output, NonPositive[#[[1]]] &][[All, 2]],
            "Input" -> QuantumTensorProduct @ Select[input, NonPositive[#[[1]]] &][[All, 2]]
        ],
        QuantumBasis[
            "Output" -> QuantumTensorProduct @ Select[output, Positive[#[[1]]] &][[All, 2]],
            "Input" -> QuantumTensorProduct @ Select[input, Positive[#[[1]]] &][[All, 2]]
        ]
    }
]

QuantumOperatorProp[qo_, "Projectors", opts___] /; qo["SquareQ"] := projector /@ SparseArray @ Chop @ qo["Eigenvectors", opts]

QuantumOperatorProp[qo_, "Diagonalize", opts___] /; qo["SquareQ"] := Block[{vectors, values},
    {values, vectors} = qo["Eigensystem", opts, "Sort" -> True];
    QuantumOperator[
        DiagonalMatrix[values],
        Take[#, UpTo[1]] & /@ qo["Order"],
        QuantumBasis[AssociationThread[Subscript["s", #] & /@ Range[Length[values]], vectors], qo["Basis"]["Options"]]
    ]
]

QuantumOperatorProp[qo_, "Transpose"] := With[{qudits = Min[qo["OutputQudits"], qo["InputQudits"]]},
    qo["Transpose", Thread[{Take[qo["OutputOrder"], qudits], Take[qo["InputOrder"], qudits]}]]
]

QuantumOperatorProp[qo_, "Transpose", order : {(List | Rule)[_Integer, _Integer]...}] := Block[{
    outputMap = MapAt[1, Rule @@@ order, {All, 2}],
    inputMap = MapAt[0, Rule @@@ Reverse /@ order, {All, 2}],
    map, out, in, qudits
},
    map = Join[Replace[qo["FullOutputOrder"], Append[outputMap, o_ :> 0[o]], {1}], Replace[qo["FullInputOrder"], Append[inputMap, i_ :> 1[i]], {1}]];
    {out, in} = {Cases[map, 0[_]], Cases[map, 1[_]]};
    qudits = Join[
        Replace[Keys[outputMap], Append[qo["OutputOrderQuditMapping"], _ -> Nothing], {1}],
        Replace[Keys[inputMap], Append[qo["InputOrderQuditMapping"], _ -> Nothing], {1}] + qo["OutputQudits"]
    ];
    QuantumOperator[
        QuantumState[qo["State"], qo["Basis"]]["Permute", FindPermutation[map, Join[out, in]]],
        Map[First, {out, in}, {2}]
    ]
]


QuantumOperatorProp[qo_, prop : "Dagger" | "ConjugateTranspose" | "Inverse"] := simplifyLabel @ QuantumOperator[
    qo["State"][prop], {qo["InputOrder"], qo["OutputOrder"]}]

QuantumOperatorProp[qo_, prop : "Conjugate" | "Dual"] := QuantumOperator[qo["State"][prop], qo["Order"]]

QuantumOperatorProp[qo_, "Bend", shift : _Integer ? Positive : Automatic] :=
    simplifyLabel @ If[ qo["MatrixQ"],
        QuantumOperator[qo["State"]["Bend"], Join[#, # + Replace[shift, Automatic :> Max[qo["Order"]] - Min[#] + 1]] & /@ qo["Order"]],
        QuantumTensorProduct[qo, qo["Conjugate"]["Shift", Replace[shift, Automatic :> Max[qo["Order"]]]]]
    ]

QuantumOperatorProp[qo_, "Unbend"] :=
    If[ qo["MatrixQ"],
        qo,
        QuantumOperator[qo["Sort"]["State"]["Unbend"], Take[Sort[#], Length[#] / 2] & /@ qo["Order"]]
    ]

QuantumOperatorProp[qo_, "Double"] := With[{state = qo["State"]["Double"]},
    QuantumOperator[
        QuantumState[
            state,
            QuantumBasis[
                "Output" -> QuantumTensorProduct[Times @@@ Partition[state["Output"]["Decompose"], 2]],
                "Input" -> QuantumTensorProduct[Times @@@ Partition[state["Input"]["Decompose"], 2]],
                "Label" -> With[{label = Replace[qo["Label"], Interpretation[_, label_] :> label]}, Interpretation[Style[label, Bold], label]],
                state["Basis"]["Options"]
            ]
        ],
        qo["Order"]
    ]
]

QuantumOperatorProp[qo_, "Undouble"] := Enclose @ QuantumOperator[
    QuantumState[qo["State"], QuantumBasis[
        Splice[{#, #}] & /@ ConfirmBy[Sqrt[qo["OutputDimensions"]], AllTrue[IntegerQ]],
        Splice[{#, #}] & /@ ConfirmBy[Sqrt[qo["InputDimensions"]], AllTrue[IntegerQ]],
        "Label" -> Replace[qo["Label"], Interpretation[_, label_] :> label]
    ]]["Undouble"],
    qo["Order"]
]

QuantumOperatorProp[qo_, "TensorReverseInput", order_ ? orderQ] :=
    QuantumOperator[qo["State"]["TensorReverseInput", order /. qo["InputOrderQuditMapping"]], qo["Order"]]

QuantumOperatorProp[qo_, "TensorReverseOutput", order_ ? orderQ] :=
    QuantumOperator[qo["State"]["TensorReverseOutput", order /. qo["OutputOrderQuditMapping"]], qo["Order"]]

QuantumOperatorProp[qo_, "Tr" | "Trace" | "TraceNorm" | "Norm"] := Tr[qo["Matrix"]]

QuantumOperatorProp[qo_, "Normalize"] := qo / Abs[qo["Norm"]]

QuantumOperatorProp[qo_, "PrimeBasis", outputQudit : _Integer | Automatic : Automatic, inputQudit : _Integer | Automatic : Automatic, opts___] :=
With[{state = qo["State"]["PrimeBasis"]},
    QuantumOperator[state, {
        If[ outputQudit === Automatic,
            # - Min[#, 1] + 1 &[Max[qo["OutputOrder"]] - Reverse @ Range[state["OutputQudits"]] + 1],
            outputQudit + Range[state["OutputQudits"]] - 1
        ],
        Which[
            inputQudit === Automatic && outputQudit === Automatic,
            # - Min[#, 1] + 1 &[Max[qo["InputOrder"]] - Reverse @ Range[state["InputQudits"]] + 1],
            inputQudit === Automatic && outputQudit =!= Automatic,
            outputQudit + Range[state["InputQudits"]] - 1,
            True,
            inputQudit + Range[state["InputQudits"]] - 1
        ]
    },
    opts,
    "Label" -> qo["Label"]
    ]
]

QuantumOperatorProp[qo_, prop : "Simplify" | "FullSimplify" | "Chop" | "ComplexExpand" | "EigenPrune", args___] :=
    Enclose @ QuantumOperator[ConfirmBy[qo["State"][prop, args], QuantumStateQ], qo["Order"]]


(* evolution *)

QuantumOperatorProp[qo_, "EvolutionOperator", args___] := QuantumEvolve[qo, None, args]

QuantumOperatorProp[qo_, "NEvolutionOperator", args___] /; qo["ParameterArity"] == 1 := Module[{
    parameter = First @ qo["Parameters"],
    parameterSpec = First @ qo["ParameterSpec"],
    initialParameter = First @ qo["InitialParameters"],
    leftEquations, rightEquations, initialEquations, equations
},
        leftEquations = Flatten[
            (1 / I) qo["MatrixRepresentation"] . Partition[
                    Subscript["u", #][parameter] & /@ Range[qo["OutputDimension"] ^ 2],
                    qo["OutputDimension"]
                ]
        ];
        rightEquations = Subscript["u", #]'[parameter] & /@ Range[qo["OutputDimension"] ^ 2];
        equations = Equal @@@ Transpose[{rightEquations, leftEquations}];
        initialEquations = Equal @@@
            Transpose @ {
                Subscript["u", #][initialParameter] & /@ Range[qo["OutputDimension"] ^ 2],
                Flatten[IdentityMatrix[qo["OutputDimension"]]]
            };
        QuantumOperator[
            Partition[Flatten @
                NDSolveValue[
                    Join[equations, initialEquations],
                    Subscript["u", #][parameter] & /@ Range[qo["OutputDimension"] ^ 2],
                    Evaluate @ parameterSpec,
                    args
                ],
                qo["OutputDimension"]
            ],
            qo["Order"],
            "ParameterSpec" -> First @ qo["ParameterSpec"]
        ]
    ]

QuantumOperatorProp[qo_, "EigenvaluePlot", args___] /; qo["ParameterArity"] == 1 :=
    ReImPlot[Evaluate @ qo["Eigenvalues"], Evaluate @ First @ qo["ParameterSpec"],
        args
    ]


QuantumOperatorProp[qo_, "PauliDecompose"] /; qo["InputQudits"] == qo["OutputQudits"] && MatchQ[qo["Dimensions"], {2 ..}] := With[{
    n = qo["InputQudits"],
    mat = {{0, 0, 1 / 2, 1 / 2}, {1 / 2, I / 2, 0, 0}, {1 / 2, - I / 2, 0, 0}, {0, 0, - 1 / 2, 1 / 2}}
},
    KeyMap[StringJoin @ Replace[#, {1 -> "X", 2 -> "Y", 3 -> "Z", 4 -> "I"}, {1}] &] @ Association @ Most @ ArrayRules[
        Nest[
            Transpose[# . mat, RotateRight[Range[n]]] &,
            Flatten[ArrayReshape[qo["Matrix"], ConstantArray[2, 2 n]], Array[{#, # + n} &, n]],
            n
        ]
    ]
]

QuantumOperatorProp[qo_, "PauliDecompose"] /; qo["InputDimensions"] == qo["OutputDimensions"] := Enclose @ Block[{
    dims = qo["InputDimensions"],
    targetValues = Flatten[qo["MatrixRepresentation"]],
    paulis, ops, params, values
},
	paulis = Tuples[{"I", "X", "Y", "Z"}, Length[dims]];
	ops = Map[MapIndexed[QuantumOperator, Thread[{#, dims}]] &, paulis];
	params = \[FormalT] @@@ paulis;
	values = Flatten @ Total @ MapThread[#1 kroneckerProduct @@ (#["Matrix"] & /@ #2) &, {params, ops}];
	DeleteCases[AssociationThread[ops, First @ ConfirmMatch[SolveValues[Thread[values == targetValues], params], {__}]], 0]
]


UnitaryEulerAngles[_, b_, c_, d_] := With[{theta = Simplify[2 ArcSin[Abs[b]]]},
    If[NumericQ[theta] && theta == 0, {0, 0, Simplify[Mod[Arg[d], 2 Pi]]}, {Simplify[Mod[Arg[c], 2 Pi]], theta, Simplify[Mod[Arg[b] - Pi, 2 Pi]]}]
]
UnitaryEulerAngles[u_ ? SquareMatrixQ] /; Dimensions[u] === {2, 2} := UnitaryEulerAngles[u[[1, 1]], u[[1, 2]] , u[[2, 1]], u[[2, 2]]]

UnitaryEulerAngles[qo_QuantumOperator] := UnitaryEulerAngles[qo["MatrixRepresentation"]]

UnitaryEulerAnglesWithPhase[u_ ? SquareMatrixQ] /; Dimensions[u] === {2, 2} := With[{a = u[[1, 1]], b = u[[1, 2]], c = u[[2, 1]], d = u[[2, 2]]},
    With[{phase = Simplify[Arg[a]], theta = Simplify[2 ArcSin[Abs[b]]]},
        {If[NumericQ[theta] && theta == 0, {0, 0, Simplify[Mod[Arg[d] - phase, 2 Pi]]}, {Simplify[Mod[Arg[c] - phase, 2 Pi]], theta, Simplify[Mod[Arg[b] - Pi - phase, 2 Pi]]}], phase}
    ]
]

UnitaryEulerAnglesWithPhase[mat_ ? MatrixQ] /; Times @@ Dimensions[mat] == 4 := UnitaryEulerAnglesWithPhase[ArrayReshape[mat, {2, 2}]]

UnitaryEulerAnglesWithPhase[qo_QuantumOperator] := UnitaryEulerAnglesWithPhase[qo["MatrixRepresentation"]]

UnitaryAngles[arg_] := UnitaryEulerAngles[arg][[{2, 1, 3}]]

UnitaryAnglesWithPhase[arg_] := MapAt[#[[{2, 1, 3}]] &, UnitaryEulerAnglesWithPhase[arg], {1}]


QuantumOperatorProp[qo_, "EulerAngles"] /; qo["VectorQ"] && qo["Dimension"] == 4 := First[UnitaryEulerAnglesWithPhase[qo["MatrixRepresentation"]]]

QuantumOperatorProp[qo_, "EulerAnglesWithPhase"] /; qo["VectorQ"] && qo["Dimension"] == 4 := UnitaryEulerAnglesWithPhase[qo["MatrixRepresentation"]]

QuantumOperatorProp[qo_, "ZYZ"] /; qo["VectorQ"] && qo["Dimension"] == 4 := Enclose @ Module[{angles, phase},
    {angles, phase} = UnitaryAnglesWithPhase[qo["MatrixRepresentation"]];
    If[ NumericQ[phase] && phase == 0,
        QuantumOperator["U" @@ angles, qo["Order"]],
        QuantumCircuitOperator[{QuantumOperator["GlobalPhase"[phase]], QuantumOperator["U" @@ angles, qo["Order"]]}]
    ]
]


(* ::Section:: Two-qubit Cartan (KAK) decomposition *)

(* Magic (Bell) basis: conjugation by it maps local SU(2) (x) SU(2) to real SO(4). *)
$MagicBasis = (1 / Sqrt[2]) {{1, 0, 0, I}, {0, I, 1, 0}, {0, I, -1, 0}, {1, 0, 0, -I}};

(* Eigenvalue signs of XX, YY, ZZ on the magic-basis Bell states (rows = Bell, cols = X/Y/Z). *)
$kakBellSigns := $kakBellSigns = With[{
    mbd = ConjugateTranspose[$MagicBasis],
    pauliPair = KroneckerProduct[PauliMatrix[#], PauliMatrix[#]] &
},
    Transpose[Chop @ Re @ Table[Diagonal[mbd . pauliPair[k] . $MagicBasis], {k, 3}]]
]

(* Split a 4x4 local operator L = KroneckerProduct[a, b] into its two single-qubit factors
   via the reshuffle (realignment) trick: the reshuffled matrix is rank one. *)
kakKroneckerFactor[mat_] := Module[{reshuffled, u, s, v},
    reshuffled = ArrayReshape[
        Table[mat[[2 (i1 - 1) + i2, 2 (j1 - 1) + j2]], {i1, 2}, {j1, 2}, {i2, 2}, {j2, 2}],
        {4, 4}
    ];
    {u, s, v} = SingularValueDecomposition[reshuffled];
    {
        ArrayReshape[Sqrt[s[[1, 1]]] u[[All, 1]], {2, 2}],
        ArrayReshape[Sqrt[s[[1, 1]]] Conjugate[v[[All, 1]]], {2, 2}]
    }
]

(* Off-diagonal mass of proj^T . m2 . proj: zero iff proj diagonalizes the complex-symmetric m2. *)
kakOffDiagonal[proj_, m2_] := With[{d = Transpose[proj] . m2 . proj},
    Norm[Flatten[d - DiagonalMatrix[Diagonal[d]]]]
]

(* Real-orthonormal simultaneous diagonalizer of the commuting pair {Re[m2], Im[m2]}.
   m2 is complex-symmetric and unitary, so a single real orthogonal basis diagonalizes it; we get it
   from Eigensystem of the real-symmetric combination Re[m2] + sep Im[m2]. A *fixed* separator merges
   two distinct m2 eigenvalues whenever their phases are mirror-symmetric about ArcTan[sep] (a
   measure-zero resonance, but reachable: the c1 = ArcTan[Sqrt[2]]/2 class with nontrivial local
   factors fails at sep = Sqrt[2]). Fast-path Sqrt[2.]; on the rare resonance, fall back to the
   separator that best diagonalizes m2. *)
kakEigenbasis[m2_] := Module[{seps = {Sqrt[2.], Sqrt[3.], Sqrt[5.], N[Pi], N[E]}, first},
    first = Transpose[Eigensystem[Re[m2] + First[seps] Im[m2]][[2]]];
    If[ kakOffDiagonal[first, m2] < 10.^-9,
        first,
        First @ MinimalBy[
            Table[Transpose[Eigensystem[Re[m2] + sep Im[m2]][[2]]], {sep, Rest[seps]}],
            kakOffDiagonal[#, m2] &
        ]
    ]
]

(* Core: factor a 4x4 unitary u as u = e^{i phi} (a1 (x) a2) . exp(i (c1 XX + c2 YY + c3 ZZ)) . (b1 (x) b2).
   Returns <|"GlobalPhase", "Left" -> {a1, a2}, "Canonical" -> {c1, c2, c3}, "Right" -> {b1, b2}|>. *)
TwoQubitKAK[u_ ? SquareMatrixQ] /; Dimensions[u] === {4, 4} := Module[{
    mb = $MagicBasis, mbd = ConjugateTranspose[$MagicBasis], signs = $kakBellSigns,
    su, uMagic, m2, proj, k1, k2, theta, coords, corr, left, right
},
    su = u / Det[u] ^ (1/4);
    uMagic = mbd . su . mb;
    m2 = Transpose[uMagic] . uMagic;
    (* commuting real-symmetric pair Re[m2], Im[m2] share a real orthonormal eigenbasis;
       kakEigenbasis picks a separator that actually diagonalizes m2 (guards the resonance) *)
    proj = kakEigenbasis[m2];
    If[Re[Det[proj]] < 0, proj[[All, 1]] = - proj[[All, 1]]];        (* K2 in SO(4) *)
    k2 = proj;
    theta = Arg[Diagonal[Transpose[k2] . m2 . k2]] / 2;
    k1 = uMagic . k2 . DiagonalMatrix[Exp[- I theta]];
    If[Re[Det[k1]] < 0, k1[[All, 1]] = - k1[[All, 1]]; theta[[1]] += Pi];   (* K1 in SO(4) *)
    coords = Transpose[signs] . theta / 4;
    corr = mb . DiagonalMatrix[Exp[I (theta - signs . coords)]] . mbd;       (* local residual *)
    left = kakKroneckerFactor[mb . k1 . mbd . corr];
    right = kakKroneckerFactor[mb . Transpose[k2] . mbd];
    <|"GlobalPhase" -> Arg[Det[u] ^ (1/4)], "Left" -> left, "Canonical" -> coords, "Right" -> right|>
]


(* Pivot-anchored Kronecker factor: read A (x) B off a nonzero 2x2 block and a nonzero pivot inside it.
   Exact and symbolic-safe (no SingularValueDecomposition, which returns a Sqrt[0] on degenerate
   symbolic locals such as the local factors of a controlled-phase gate). *)
kakKroneckerFactorSymbolic[mat_, simp_] := Module[{blocks, a0, c0, b, r0, s0},
    blocks = Table[mat[[2 (a - 1) + 1 ;; 2 a, 2 (c - 1) + 1 ;; 2 c]], {a, 2}, {c, 2}];
    {a0, c0} = First @ Position[blocks, blk_ /; AnyTrue[Flatten[blk], ! PossibleZeroQ[#] &], {2}, Heads -> False];
    b = blocks[[a0, c0]];
    {r0, s0} = First @ Position[b, e_ /; ! PossibleZeroQ[e], {2}, Heads -> False];
    {Table[simp[blocks[[a, c]][[r0, s0]] / b[[r0, s0]]], {a, 2}, {c, 2}], b}
]

(* Continuous angle of a unit-modulus symbolic expression. Arg gives the principal value
   (Arg[E^{-I theta}] is theta mod 2 Pi, not theta), which leaves the magic phases unable to cancel and
   makes K1 nonlocal for gates such as RXX. Reading the exponent off the Exp form instead recovers the
   continuous angle exactly. Valid here because every input (the diagonalized m2 eigenvalues and
   Det[u]^(1/4)) is unit modulus under the real-parameter assumption. *)
kakUnwrapAngle[z_] := - I PowerExpand[Log[Simplify[TrigToExp[z]]]]

(* Structured-symbolic KAK: the same magic-basis factorization in exact arithmetic. ComplexExpand makes
   Re/Im split once the parameters are declared real, an exact Sqrt[2] separates the eigenvalues,
   Orthogonalize forces a real-orthonormal eigenbasis, kakUnwrapAngle reads off the continuous angles,
   and FullSimplify[..., assum] reduces to closed form. Returns the same association as TwoQubitKAK (for
   a generic dense symbolic matrix the eigensolve is a Root blowup; the caller bounds it). *)
TwoQubitKAKSymbolic[u_ ? SquareMatrixQ, assum_] /; Dimensions[u] === {4, 4} := Module[{
    mb = $MagicBasis, mbd = ConjugateTranspose[$MagicBasis], signs = $kakBellSigns,
    simp, su, uMagic, m2, proj, k1, theta, coords, corr
},
    simp = FullSimplify[#, assum] &;
    su = u / Det[u] ^ (1/4);
    uMagic = simp @ ComplexExpand[mbd . su . mb];
    m2 = simp @ ComplexExpand[Transpose[uMagic] . uMagic];
    proj = simp @ ComplexExpand @ Transpose[Orthogonalize[
        Eigensystem[ComplexExpand[Re[m2]] + Sqrt[2] ComplexExpand[Im[m2]]][[2]]]];
    If[TrueQ[simp[Re[Det[proj]]] < 0], proj[[All, 1]] = - proj[[All, 1]]];
    theta = simp[kakUnwrapAngle[Diagonal[Transpose[proj] . m2 . proj]] / 2];
    k1 = simp @ ComplexExpand[uMagic . proj . DiagonalMatrix[Exp[- I theta]]];
    If[TrueQ[simp[Re[Det[k1]]] < 0], k1[[All, 1]] = - k1[[All, 1]]; theta[[1]] += Pi];
    coords = simp[Transpose[signs] . theta / 4];
    corr = simp @ ComplexExpand[mb . DiagonalMatrix[Exp[I (theta - signs . coords)]] . mbd];
    <|
        "GlobalPhase" -> simp @ kakUnwrapAngle[Det[u] ^ (1/4)],
        "Left" -> kakKroneckerFactorSymbolic[simp @ ComplexExpand[mb . k1 . mbd . corr], simp],
        "Canonical" -> coords,
        "Right" -> kakKroneckerFactorSymbolic[simp @ ComplexExpand[mb . Transpose[proj] . mbd], simp]
    |>
]

(* Wall-clock budget for the symbolic eigensolve before bailing to the undecomposed-operator path. *)
$kakSymbolicTimeBudget = 60;

(* exp(i c P (x) P) gadget as a circuit fragment on {o1, o2}: a CNOT, an Rz on the target, a CNOT,
   conjugated into the X / Y / Z basis. Single-qubit pieces are emitted as "U" gates (phaseless);
   the dropped global phases are recovered once, at the end, by projection. *)
$kakRotZ = MatrixExp[-I # / 2 PauliMatrix[3]] &;
$kakBasisX = (PauliMatrix[1] + PauliMatrix[3]) / Sqrt[2];                      (* Hadamard *)
$kakBasisY = MatrixExp[-I Pi / 4 PauliMatrix[1]];                            (* Rx[Pi/2] up to phase *)

(* Use the with-phase angles: "U" @@ UnitaryAngles[mat] alone drops Arg[mat[[1,1]]], which is a
   *relative* (not global) phase error for gates like Rz. The with-phase angles represent mat up to
   a single global phase, which folds into the final projected GlobalPhase. Chop first: a near-diagonal
   factor (e.g. the identity local of X (x) I) carries tiny off-diagonal noise that fails the exact
   theta == 0 test in UnitaryEulerAngles and sends Arg of numerical zeros into the angles. *)
kakUGate[mat_, order_] := QuantumOperator["U" @@ First[UnitaryAnglesWithPhase[Chop[mat]]], order]

kakInteractionGadget[c_, basis_, {o1_, o2_}] := If[TrueQ[Chop[c] == 0] || PossibleZeroQ[c], {},
    {
        If[basis === IdentityMatrix[2], Nothing, kakUGate[basis, {o1}]],
        If[basis === IdentityMatrix[2], Nothing, kakUGate[basis, {o2}]],
        QuantumOperator["CX", {o1, o2}],
        kakUGate[$kakRotZ[-2 c], {o2}],
        QuantumOperator["CX", {o1, o2}],
        If[basis === IdentityMatrix[2], Nothing, kakUGate[ConjugateTranspose[basis], {o1}]],
        If[basis === IdentityMatrix[2], Nothing, kakUGate[ConjugateTranspose[basis], {o2}]]
    }
]

(* Symbolic emission: single-qubit pieces go in as direct matrix gates rather than "U" angles.
   Re-extracting Euler angles from an Arg / ArcTan-laden symbolic factor is unreliable; a matrix gate
   is exact and preserves the factor's own phase, so the whole circuit equals the operator up to the
   single global phase Arg[Det[u]^(1/4)] returned by the worker, with no projection needed. *)
kakInteractionGadgetSymbolic[c_, basis_, {o1_, o2_}] := If[PossibleZeroQ[c], {},
    {
        If[basis === IdentityMatrix[2], Nothing, QuantumOperator[basis, {o1}]],
        If[basis === IdentityMatrix[2], Nothing, QuantumOperator[basis, {o2}]],
        QuantumOperator["CX", {o1, o2}],
        QuantumOperator[$kakRotZ[-2 c], {o2}],
        QuantumOperator["CX", {o1, o2}],
        If[basis === IdentityMatrix[2], Nothing, QuantumOperator[ConjugateTranspose[basis], {o1}]],
        If[basis === IdentityMatrix[2], Nothing, QuantumOperator[ConjugateTranspose[basis], {o2}]]
    }
]

QuantumCircuitOperator::kaksymbolic = "KAK could not decompose `1` symbolically (a generic dense symbolic 2-qubit unitary is a Root blowup); it is returned undecomposed. Substitute numeric parameters, or pass a structured / parametric gate.";

QuantumOperatorProp[qo_, "KAK"] /; qo["Dimension"] == 16 && MatchQ[qo["Dimensions"], {2, 2, 2, 2}] :=
Enclose @ Module[{
    mat = Normal @ qo["MatrixRepresentation"], numericQ, vars, assum, data, c, a, b, o1, o2, order,
    gates, circuit, circMat, residualPhase
},
    numericQ = MatrixQ[Quiet @ N @ mat, NumericQ];
    order = qo["Order"][[1]];
    {o1, o2} = order;
    If[ numericQ,
    (* ---- numeric path: "U"-angle gadgets, global phase recovered by projection ---- *)
        data = TwoQubitKAK[N @ mat];
        {c, a, b} = {data["Canonical"], data["Left"], data["Right"]};
        (* circuit order applies first -> last: right locals, ZZ, YY, XX gadgets, left locals *)
        gates = Flatten @ {
            kakUGate[b[[1]], {o1}], kakUGate[b[[2]], {o2}],
            kakInteractionGadget[c[[3]], IdentityMatrix[2], {o1, o2}],
            kakInteractionGadget[c[[2]], $kakBasisY, {o1, o2}],
            kakInteractionGadget[c[[1]], $kakBasisX, {o1, o2}],
            kakUGate[a[[1]], {o1}], kakUGate[a[[2]], {o2}]
        };
        circMat = Normal @ (QuantumOperator @ QuantumCircuitOperator[gates])["MatrixRepresentation"];
        residualPhase = Arg[Tr[ConjugateTranspose[N @ circMat] . N @ mat] / 4];
        QuantumCircuitOperator[
            If[Chop[residualPhase] == 0, gates, Prepend[gates, QuantumOperator["GlobalPhase"[residualPhase]]]],
            "Label" -> "KAK"
        ]
    ,
    (* ---- structured-symbolic path: direct matrix gates, global phase = worker's gp (no projection).
       A generic dense symbolic eigensolve is a Root blowup, so the worker is time-bounded and bails. ---- *)
        (* the free parameters: actual non-System symbols (Variables misses these inside Cos / E^...) *)
        vars = DeleteDuplicates @ Cases[mat, _Symbol ? (Context[#] =!= "System`" &), Infinity];
        assum = Element[vars, Reals];
        data = TimeConstrained[TwoQubitKAKSymbolic[mat, assum], $kakSymbolicTimeBudget, $Failed];
        If[ ! AssociationQ[data],
            Message[QuantumCircuitOperator::kaksymbolic, qo["Label"]];
            Return[QuantumCircuitOperator[{qo}, "Label" -> "KAK"], Module]
        ];
        {c, a, b} = {data["Canonical"], data["Left"], data["Right"]};
        gates = Flatten @ {
            QuantumOperator[b[[1]], {o1}], QuantumOperator[b[[2]], {o2}],
            kakInteractionGadgetSymbolic[c[[3]], IdentityMatrix[2], {o1, o2}],
            kakInteractionGadgetSymbolic[c[[2]], $kakBasisY, {o1, o2}],
            kakInteractionGadgetSymbolic[c[[1]], $kakBasisX, {o1, o2}],
            QuantumOperator[a[[1]], {o1}], QuantumOperator[a[[2]], {o2}]
        };
        residualPhase = data["GlobalPhase"];
        circuit = QuantumCircuitOperator[
            If[PossibleZeroQ[residualPhase], gates, Prepend[gates, QuantumOperator["GlobalPhase"[residualPhase]]]],
            "Label" -> "KAK"
        ];
        (* Safety self-check: symbolic eigendecomposition under a generic real assumption is not always
           reliable (the eigenvectors, hence the local factors, can come out wrong even when the
           coordinates are right). Verify the circuit reproduces the operator at a probe substitution;
           if it does not, return the operator undecomposed rather than a silently wrong circuit. *)
        Module[{probeSub, checkErr},
            probeSub = Thread[vars -> Table[Sqrt[Prime[k + 1]] / 5, {k, Max[Length[vars], 1]}]];
            checkErr = Chop @ Norm[Flatten[
                N[Normal[(QuantumOperator @ circuit)["MatrixRepresentation"]] /. probeSub] - N[mat /. probeSub]], Infinity];
            If[ TrueQ[checkErr < 10.^-6],
                circuit,
                Message[QuantumCircuitOperator::kaksymbolic, qo["Label"]];
                QuantumCircuitOperator[{qo}, "Label" -> "KAK"]
            ]
        ]
    ]
]

(* Umbrella synthesis into the gate set: 1-qubit -> ZYZ, 2-qubit -> KAK. Larger operators fall
   through to the standard undefined-property path (the n>=3 QSD slot is future work). *)
QuantumOperatorProp[qo_, "Decompose", gateset : {___String} : {"U", "CX"}] /;
    qo["Dimension"] == 4 && qo["VectorQ"] := qo["ZYZ"]

QuantumOperatorProp[qo_, "Decompose", gateset : {___String} : {"U", "CX"}] /;
    qo["Dimension"] == 16 && MatchQ[qo["Dimensions"], {2, 2, 2, 2}] := qo["KAK"]

(* OpenQASM emit workers: the public ["QASM"]/["SimpleQASM"] downvalues route through the
   QuantumQASM hub, which calls these. *)

(* The single-qubit-unitary branch fires only on a genuine 1-qubit gate. SquareQ
   (OutputDimension == InputDimension) is required because non-square state-injection
   tensors share the {2, 2} leg signature: a 2-output/0-input Cup also reports
   Dimensions === {2, 2} but has a 4x1 matrix, so without this guard UnitaryAngles is
   handed a non-square matrix and the emit fails opaquely. Non-square operators fall
   through to the catch-all, which emits the "// Unimplemented" marker that
   qasmEmitCircuit turns into a clean QuantumQASM failure naming the gate. *)
(* Format a real for an OpenQASM literal: full machine precision (NumberForm truncates to 6
   digits, capping round-trip fidelity, and renders small/large magnitudes as a "x10"
   superscript that wraps across lines, which is not a valid literal). InputForm gives a flat
   decimal; "*^" becomes the OpenQASM "e" exponent (the importer parses both). *)
qasmReal[x_] := StringReplace[ToString[N[x], InputForm], "*^" -> "e"]

(* Decompose a 1-qubit gate as U(angles) up to a global phase. UnitaryAnglesWithPhase, not
   UnitaryAngles: the latter folds Arg[mat[[1,1]]] into the phi/lambda angles, which is a
   diagonal Z-conjugation (mat == D.U(UnitaryAngles).D, D = diag(E^(I a/2), E^(-I a/2))),
   NOT a global phase, so U(UnitaryAngles[mat]) is a different operator whenever mat[[1,1]]
   is complex. The leftover phase here is a true global phase and is dropped (callers are all
   uncontrolled single-qubit contexts; the controlled path emits it as a controlled gphase). *)
qasmEmitSimple[qo_] /; qo["SquareQ"] && qo["Dimensions"] === {2, 2} := Enclose @ With[{
    angles = qasmReal /@ ConfirmBy[First @ UnitaryAnglesWithPhase[qo["MatrixRepresentation"]], AllTrue[NumericQ]]
},
    StringTemplate["U(``, ``, ``)"][Sequence @@ angles]
]

qasmEmitSimple[qo_] :=
    Replace[qo["Label"], {
        "SWAP" | "\[Pi]"[_, _] :> "swap",
        Superscript[label_, CircleTimes[_]] :> qasmEmitSimple[QuantumOperator[label]],
        label_ :> StringTemplate["// Unimplemented QASM for operator with label: `` ----"][label]
    }]

QuantumOperatorProp[qo_, "SimpleQASM"] := QuantumQASM[qo, "Simple"]

QuantumOperatorProp[qo_, "TargetOperator"] := Module[{control1, control0, n, m},
    control1 = qo["ControlOrder1"];
    control0 = qo["ControlOrder0"];
    n = Length[control1];
    m = Length[control0];
    If[n + m == 0, Return[qo]];
    QuantumOperator[
        QuantumTensorProduct[
            QuantumOperator[QuantumState["Register"[n, 2 ^ n - 1]]["Dagger"], {{}, control1}],
            QuantumOperator[QuantumState["Register"[m, 0]]["Dagger"], {{}, control0}]
        ] @ qo @
        QuantumTensorProduct[
            QuantumOperator[QuantumState["Register"[n, 2 ^ n - 1]], {control1, {}}],
            QuantumOperator[QuantumState["Register"[m, 0]], {control0, {}}]
        ],
        "Label" -> Replace[qo["Label"], Subscript["C", label_][__] :> label]
    ]
]

qasmEmitOperator[qo_] /; qo["ControlOrder"] =!= {} && MatchQ[qo["TargetOrder"], {_}] :=
    Replace[qo["Label"], {
        Subscript["C", _][control1_, control0_] :>
            With[{n = Length[control1], m = Length[control0]},
                Module[{proj, mat, angles, phase, ctrlPrefix, ctrlQubits, gateLine, phaseLine},
                    proj =
                        QuantumTensorProduct[{
                            If[ n > 0,
                                QuantumOperator[QuantumState["Register"[n, 2 ^ n - 1]]["Dagger"], {{}, control1}],
                                Nothing
                            ],
                            If[ m > 0,
                                QuantumOperator[QuantumState["Register"[m, 0]]["Dagger"], {{}, control0}],
                                Nothing
                            ]
                        }] @ qo @
                        QuantumTensorProduct[{
                            If[ n > 0,
                                QuantumOperator[QuantumState["Register"[n, 2 ^ n - 1]], {control1, {}}],
                                Nothing
                            ],
                            If[ m > 0,
                                QuantumOperator[QuantumState["Register"[m, 0]], {control0, {}}],
                                Nothing
                            ]
                        }];
                    mat = proj["MatrixRepresentation"];
                    If[ ! MatchQ[Dimensions[mat], {2, 2}], Return[$Failed, Module]];
                    (* Split the 1-qubit target into a U(angles) carrying no global phase plus the
                       leftover phase: mat == E^(I phase) . U(angles). A bare U cannot carry that
                       phase (its (1,1) entry is real). Uncontrolled it is unobservable, but under
                       control it is a physical relative phase, so emit it as a controlled gphase on
                       the control qubits alongside the controlled U. *)
                    {angles, phase} = UnitaryAnglesWithPhase[mat];
                    ctrlPrefix = StringTemplate["ctrl(``) @ negctrl(``) @ "][n, m];
                    ctrlQubits = StringRiffle[Map[StringTemplate["q[``]"], Join[control1, control0] - 1], " "];
                    gateLine = ctrlPrefix <> StringTemplate["U(``, ``, ``)"][Sequence @@ (qasmReal /@ angles)] <>
                        " " <> StringRiffle[Map[StringTemplate["q[``]"], Join[control1, control0, qo["TargetOrder"]] - 1], " "] <> ";";
                    phaseLine = If[ TrueQ[Chop[N[phase]] == 0],
                        Nothing,
                        ctrlPrefix <> "gphase(" <> qasmReal[phase] <> ") " <> ctrlQubits <> ";"
                    ];
                    StringRiffle[DeleteCases[{gateLine, phaseLine}, Nothing], "\n"]
                ]
            ],
        _ :> $Failed
        }
    ]

(* a global phase is a 0-qudit (1x1) operator e^(i phi): OpenQASM 3 gphase(phi) *)
qasmEmitOperator[qo_] /; qo["Dimensions"] === {} :=
    "gphase(" <> qasmReal[Arg[qo["MatrixRepresentation"][[1, 1]]]] <> ");"

(* a qudit permutation decomposes into SWAP gates: each cycle (c1 c2 ... ck) becomes the
   transpositions (c1 c2)(c1 c3) ... (c1 ck), applied left to right (verified by matrix). *)
qasmEmitOperator[qo_] /; MatchQ[qo["Label"], "\[Pi]"[__]] := With[{
    order = qo["InputOrder"],
    transpositions = Catenate[
        Function[cyc, {First[cyc], #} & /@ Rest[cyc]] /@ PermutationCycles[List @@ qo["Label"]][[1]]
    ]
},
    StringRiffle[
        ("swap q[" <> ToString[order[[#[[1]]]] - 1] <> "] q[" <> ToString[order[[#[[2]]]] - 1] <> "];") & /@ transpositions,
        "\n"
    ]
]

qasmEmitOperator[qo_] /; MatchQ[qo["Dimensions"], {2 ..}] :=
    qasmEmitSimple[qo] <> " " <> StringRiffle[Map[StringTemplate["q[``]"], qo["InputOrder"] - 1], " "] <> ";"

QuantumOperatorProp[qo_, "QASM"] := QuantumQASM[qo]


QuantumOperatorProp[qo_, "CircuitDiagram", opts___] := QuantumCircuitOperator[qo]["Diagram", opts]


(* state properties *)

QuantumOperatorProp[qo_, args : PatternSequence[prop_String, ___] | PatternSequence[{prop_String, ___}, ___]] /;
    MemberQ[Intersection[qo["State"]["AllProperties"], qo["AllProperties"]], prop] := qo["State"][args]

