Package["Wolfram`QuantumFramework`"]

PackageScope["roundoffChopOperator"]



$QuantumMeasurementOperatorProperties = {
    "QuantumOperator", "Targets", "Target", "TargetCount",
    "TargetIndex",
    "Operator", "Basis", "MatrixRepresentation", "POVMElements",
    "OrderedMatrixRepresentation", "OrderedPOVMElements",
    "Arity", "Eigenqudits", "Dimensions", "Order", "HermitianQ", "UnitaryQ",
    "Eigenbasis", "Eigenvalues", "EigenvalueVectors", "Eigenvectors",
    "Eigendimensions", "Eigendimension",
    "StateDimensions", "StateDimension",
    "TargetDimensions", "TargetDimension",
    "StateQudits", "TargetBasis", "StateBasis", "CanonicalBasis", "Canonical",
    "ProjectionQ", "POVMQ",
    "SuperOperator", "POVM",
    "Shift"
};


QuantumMeasurementOperator["Properties"] := Union @ Join[
    $QuantumMeasurementOperatorProperties,
    $QuantumOperatorProperties
]

qmo_QuantumMeasurementOperator["ValidQ"] := QuantumMeasurementOperatorQ[qmo]


QuantumMeasurementOperator::undefprop = "QuantumMeasurementOperator property `` is undefined for this operator";


$QuantumMeasurementOperatorPreventCache = {"Properties", "AllProperties", "QuantumOperator", "Operator", "SuperOperator", "DiscardExtraQudits", "Basis"}

(qmo_QuantumMeasurementOperator[prop_ ? propQ, args___]) /; QuantumMeasurementOperatorQ[qmo] := With[{
    result = QuantumMeasurementOperatorProp[qmo, prop, args]
    },
    If[ TrueQ[$QuantumFrameworkPropCache] &&
        ! MemberQ[$QuantumMeasurementOperatorPreventCache, prop] &&
        QuantumMeasurementOperatorProp[qmo, "Basis"]["ParameterArity"] == 0,
        cacheProperty[QuantumMeasurementOperatorProp[qmo, prop, args], result],
        result
    ] /; !FailureQ[Unevaluated @ result] && (!MatchQ[result, _QuantumMeasurementOperatorProp] || Message[QuantumMeasurementOperator::undefprop, prop])
]

CacheProperty[QuantumMeasurementOperator][args___, value_] := PrependTo[
    DownValues[QuantumMeasurementOperatorProp],
    HoldPattern[QuantumMeasurementOperatorProp[args]] :> value
]

QuantumMeasurementOperatorProp[qmo_, "Properties"] := Union @ $QuantumMeasurementOperatorProperties

QuantumMeasurementOperatorProp[qmo_, "AllProperties"] :=
    DeleteDuplicates @ Join[QuantumMeasurementOperator["Properties"], qmo["Operator"]["AllProperties"]]


(* getters *)

QuantumMeasurementOperatorProp[_[op_, _], "Operator" | "QuantumOperator"] := op

(* direct getter so qmo["Basis"] does not fall through to the AllProperties-intersection delegation *)
QuantumMeasurementOperatorProp[_[op_, _], "Basis"] := op["Basis"]

QuantumMeasurementOperatorProp[_[_, targets_], "Targets"] := targets

QuantumMeasurementOperatorProp[qmo_, "Target" | "TargetOrder"] := Join @@ qmo["Targets"]

QuantumMeasurementOperatorProp[qmo_, "Arity" | "TargetCount"] := Length[qmo["Target"]]

QuantumMeasurementOperatorProp[qmo_, "Eigenorder"] := Replace[Select[qmo["FullOutputOrder"], NonPositive], {} -> {0}]

QuantumMeasurementOperatorProp[qmo_, "Eigenindex"] :=
    Catenate @ Position[qmo["FullOutputOrder"], _ ? NonPositive, {1}]

QuantumMeasurementOperatorProp[qmo_, "TargetIndex"] :=
    Catenate @ Lookup[PositionIndex[qmo["FullOutputOrder"]], qmo["Target"]]

QuantumMeasurementOperatorProp[qmo_, "TargetDimensions"] :=
    Part[qmo["OutputDimensions"], qmo["TargetIndex"]]

QuantumMeasurementOperatorProp[qmo_, "TargetDimension"] := Times @@ qmo["TargetDimensions"]

QuantumMeasurementOperatorProp[qmo_, "ExtraQudits"] := Count[qmo["OutputOrder"], _ ? NonPositive]

QuantumMeasurementOperatorProp[qmo_, "Eigenqudits"] := Max[qmo["ExtraQudits"], 1]

QuantumMeasurementOperatorProp[qmo_, "Eigendimensions"] :=
    If[qmo["ExtraQudits"] > 0, qmo["OutputDimensions"][[qmo["Eigenindex"]]], {Times @@ qmo["OutputDimensions"][[qmo["TargetIndex"]]]}]

QuantumMeasurementOperatorProp[qmo_, "Eigendimension"] := Times @@ qmo["Eigendimensions"]

QuantumMeasurementOperatorProp[qmo_, "Eigenbasis"] := With[{povm = qmo["POVM"]}, povm["Output"]["Extract", povm["Eigenindex"]]]

QuantumMeasurementOperatorProp[qmo_, "Eigenvalues"] := qmo["Eigenbasis"]["Names"]

QuantumMeasurementOperatorProp[qmo_, "Eigenvectors"] := qmo["Eigenbasis"]["Elements"]

QuantumMeasurementOperatorProp[qmo_, "EigenvalueVectors"] := Replace[Normal /@ qmo["Eigenvalues"], {Interpretation[_, {v_, _}] :> Replace[v, _List :> Splice[v]], v_ :> Ket[{v}]}, {2}]

QuantumMeasurementOperatorProp[qmo_, "StateQudits"] := qmo["OutputQudits"] - qmo["ExtraQudits"]

QuantumMeasurementOperatorProp[qmo_, "StateDimensions"] := Drop[qmo["Dimensions"], qmo["ExtraQudits"]]

QuantumMeasurementOperatorProp[qmo_, "StateDimension"] := Times @@ qmo["StateDimensions"]

QuantumMeasurementOperatorProp[qmo_, "TargetBasis"] := qmo["Output"]["Extract", qmo["TargetIndex"]]

QuantumMeasurementOperatorProp[qmo_, "StateBasis"] :=
    QuantumBasis[qmo["Basis"], "Output" -> Last @ qmo["Output"]["Split", qmo["ExtraQudits"]], "Input" -> qmo["Input"]]

QuantumMeasurementOperatorProp[qmo_, "CanonicalBasis"] :=
    QuantumBasis[qmo["Basis"], "Output" -> QuantumTensorProduct[qmo["TargetBasis"]["Reverse"], qmo["StateBasis"]["Output"]], "Input" -> qmo["Input"]]


canonicalEigenPermutation[qmo_] := Block[{accumIndex = PositionIndex[FoldList[Times, qmo["TargetDimensions"]]]},
	FindPermutation @ Catenate[
        Reverse /@ TakeList[
            Range[qmo["TargetCount"]],
            Reverse @ Differences @ Prepend[0] @ Catenate @ Lookup[accumIndex, FoldList[Times, Reverse[qmo["Eigendimensions"]]]]
        ]
    ]
]


(* A projective measurement is built from the eigenvectors of its operator's
   "MatrixRepresentation", whose tensor factors follow the sorted qudit order, and its target
   names the measured qudits. The rules below for these properties build the measurement in
   that sorted order and then attach the operator's stored order or the target list to the
   result, which is right only when both are sorted. So an operator stored on an unsorted
   order, or a target listed out of order, is measured as the sorted operator on the sorted
   target: the same observable on the same qudits. This rule has to stay above those rules,
   which would otherwise answer first. *)
QuantumMeasurementOperatorProp[qmo_,
    prop : "SuperOperator" | "POVM" | "Canonical" | "DiscardExtraQudits" | "ReverseEigenQudits" | "Double" | "Conjugate" | "Dual" | "Unbend",
    args___
] /; qmo["ProjectionQ"] && ! (qmo["QuantumOperator"]["SortedQ"] && OrderedQ[qmo["Target"]]) :=
    QuantumMeasurementOperator[qmo["QuantumOperator"]["Sort"], {Sort @ qmo["Target"]}][prop, args]



QuantumMeasurementOperatorProp[qmo_, "Canonical", OptionsPattern[{"Reverse" -> True, "CanonicalBasis" -> True}]] /; qmo["Eigendimension"] == qmo["TargetDimension"] := With[{
    basis = qmo["CanonicalBasis"],
    perm = canonicalEigenPermutation[qmo]
    (* perm = EchoLabel[{qmo["Targets"], canonicalEigenPermutation[qmo]}] @ FindPermutation[Join @@ Sort /@ qmo["Targets"], ReverseSort[qmo["Target"]]] *)
},
    QuantumMeasurementOperator[
        QuantumOperator[
            QuantumState[
                QuantumState[
                    qmo["SuperOperator"]["State"],
                    QuantumBasis[Join[Permute[Reverse[qmo["TargetDimensions"]], InversePermutation[perm]], qmo["StateBasis"]["OutputDimensions"]], basis["InputDimensions"]]
                ]["PermuteOutput", perm],
                basis
            ]["PermuteOutput", PermutationProduct[
                FindPermutation[Reverse[qmo["Target"]], ReverseSort[qmo["Target"]]],
                If[TrueQ[OptionValue["Reverse"]], FindPermutation[Reverse[Range[qmo["TargetCount"]]]], Cycles[{}]]
            ]],
            {Join[Range[- qmo["TargetCount"] + 1, 0], DeleteCases[qmo["OutputOrder"], _ ? NonPositive]], qmo["InputOrder"]}
        ],
        Sort @ qmo["Target"]
    ] /; TrueQ[OptionValue["CanonicalBasis"]]
]

QuantumMeasurementOperatorProp[qmo_, "Canonical", OptionsPattern[{"Reverse" -> True, "CanonicalBasis" -> False}]] /; Length[qmo["Eigenorder"]] == qmo["TargetCount"] := QuantumMeasurementOperator[
        QuantumOperator[
            qmo["SuperOperator"]["State"]["PermuteOutput", PermutationProduct[
                FindPermutation[Reverse[qmo["Target"]], ReverseSort[qmo["Target"]]],
                If[TrueQ[OptionValue["Reverse"]], FindPermutation[Reverse[Range[qmo["TargetCount"]]]], Cycles[{}]]
            ]],
            {Join[Range[- qmo["TargetCount"] + 1, 0], DeleteCases[qmo["OutputOrder"], _ ? NonPositive]], qmo["InputOrder"]}
        ],
        Sort @ qmo["Target"]
    ]

(* TODO: sort target by modifying eigenbasis *)
QuantumMeasurementOperatorProp[qmo_, "Canonical", OptionsPattern[{"Reverse" -> True, "CanonicalBasis" -> False}]] := QuantumMeasurementOperator[
        QuantumOperator[
            qmo["SuperOperator"]["State"]["PermuteOutput",
                If[TrueQ[OptionValue["Reverse"]], FindPermutation[Reverse[Range[qmo["Eigenqudits"]]]], Cycles[{}]]
            ],
            {Join[Range[- qmo["Eigenqudits"] + 1, 0], DeleteCases[qmo["OutputOrder"], _ ? NonPositive]], qmo["InputOrder"]}
        ],
        qmo["Target"]
    ]


QuantumMeasurementOperatorProp[qmo_, "SortTarget"] := qmo["Canonical", "Reverse" -> False, "CanonicalBasis" -> False]

QuantumMeasurementOperatorProp[qmo_, "Sort", args___] := QuantumMeasurementOperator[qmo["QuantumOperator"]["Sort", args], qmo["Target"]]["SortTarget"]


QuantumMeasurementOperatorProp[qmo_, "ReverseEigenQudits"] := QuantumMeasurementOperator[
    QuantumOperator[
        qmo["SuperOperator"]["State"]["PermuteOutput", FindPermutation[Reverse[Range[qmo["Eigenqudits"]]]]],
        qmo["Order"]
    ],
    qmo["Targets"]
]


QuantumMeasurementOperatorProp[qmo_, "Type"] := Which[
    Count[qmo["OutputOrder"], _ ? NonPositive] == 0 && qmo["OutputDimensions"] == qmo["InputDimensions"],
    "Projection",
    Count[qmo["OutputOrder"], _ ? NonPositive] > 0,
    "POVM",
    True,
    "Unknown"
]

QuantumMeasurementOperatorProp[qmo_, "ProjectionQ"] := qmo["Type"] === "Projection"

QuantumMeasurementOperatorProp[qmo_, "POVMQ"] := qmo["Type"] === "POVM"

QuantumMeasurementOperatorProp[qmo_, "POVMElements"] := # / Mean[Diagonal[Total[#]]] & @ (# . # & /@ Through[Values[qmo["Operators"]]["MatrixRepresentation"]])

QuantumMeasurementOperatorProp[qmo_, "Operators"] /; qmo["POVMQ"] :=
    AssociationThread[qmo["Eigenvalues"], QuantumOperator[#, {Drop[qmo["OutputOrder"], qmo["Eigenqudits"]], qmo["InputOrder"]}, qmo["StateBasis"]] & /@
        ArrayReshape[qmo["Tensor"], Catenate @ MapAt[{Times @@ #} &, {1}] @ TakeDrop[qmo["Dimensions"], qmo["Eigenqudits"]]]
    ]

(* A projective measurement acts through its "POVM" form, whose Kraus operators are the
   projectors onto the orthonormal eigenbasis "SuperOperator" builds, one per outcome
   under the labels of "Eigenvalues"; so its operators are those. *)
QuantumMeasurementOperatorProp[qmo_, "Operators"] /; qmo["ProjectionQ"] := qmo["POVM"]["Operators"]

(* The partial trace leaves roundoff where the traced operator has zero entries, and inside a
   degenerate eigenspace the eigenvectors Eigensystem returns would follow it. Each entry of the
   traced operator sums as many entries of the operator of qmo as the traced qudits have
   dimensions, so its roundoff is measured against that operator before the trace: an entry of
   the traced state below 100 N 10^-p times the largest entry of the operator's state, for the
   N x N operator at precision p, or below 10^-10 times it where that is smaller, as it is
   below about 12 digits, is written 0 (roundoffChop), and the traced operator keeps its
   basis. That writes 0 for all of a traced operator that cancels, as X (x) Z does on its first
   qudit in any frame, where roundoff alone leaves numbers. The threshold scales with the
   operator, so a relative structure below it is written 0 at every scale: the off-diagonal
   entries of 10^6 (1 + 10^-15 X) are, as those of 1 + 10^-15 X are, and it is measured in the
   computational basis. A threshold that does not scale would also write 0 for entries that
   are not roundoff: all of an operator whose entries are all small, and the off-diagonal
   entries of 1 + 10^-11 X, whose eigenvectors |+> and |-> would give way to the computational
   basis that Eigensystem returns for the identity. Numbers above the threshold are kept, even
   when roundoff accumulated in building the operator put them there, as in an observable
   evolved over a long time: inside a twofold eigenspace they then pick the eigenvectors. A
   traced operator whose state is not an array of numbers gets Chop, of its state and of its
   basis. *)
roundoffChopOperator[traced_, qmo_] := With[
    {qo = traced["QuantumOperator"], operator = qmo["QuantumOperator"]},
    {state = qo["State"], reference = operator["State"]["State"]},
    If[ ArrayQ[state["State"], _, NumericQ] && ArrayQ[reference, _, NumericQ],
        QuantumMeasurementOperator[
            QuantumOperator[QuantumState[roundoffChop[state["State"], operator["OutputDimension"], reference], state["Basis"]], qo["Order"]],
            traced["Target"]
        ],
        Chop[traced]
    ]
]

QuantumMeasurementOperatorProp[qmo_, "SuperOperator", defaultEigenvalues_ : Automatic] := Module[{
    trace,
    traceQudits,
    traceDimension,
    tracedOperator,
    eigenvalues, eigenvectors, projectors,
    eigenBasis, outputBasis, inputBasis, operator
},
    trace = DeleteCases[qmo["FullInputOrder"], Alternatives @@ qmo["Target"]];
    traceQudits = trace - Min[qmo["FullInputOrder"]] + 1;
    traceDimension = Times @@ Lookup[qmo["QuantumOperator"]["InputOrderDimensions"], trace];
    If[
        ! qmo["ProjectionQ"],

        qmo["Operator"],

        tracedOperator = roundoffChopOperator[Simplify @ QuantumPartialTrace[qmo, trace], qmo];

        (* "Orthogonalize" repairs the degenerate-eigenspace basis Eigensystem returns:
           without it the eigen-projectors do not resolve the identity and the outcome
           probabilities of a degenerate observable are silently wrong. *)
        {eigenvalues, eigenvectors} = profile["Eigensystem"] @ Simplify @ tracedOperator["Eigensystem", "Sort" -> True, "Orthogonalize" -> True];
        (* The observable measured on the target is the partial trace over the other qudits
           divided by their dimension, so an operator that acts as A on the target and as the
           identity elsewhere is measured as A. The division scales the eigenvalues and leaves
           the eigenvectors alone. *)
        eigenvalues = PadRight[Replace[defaultEigenvalues, Automatic -> eigenvalues / traceDimension], Length[eigenvectors], 0];
        projectors = projector /@ eigenvectors;

        eigenBasis = QuditBasis[
            MapIndexed[
                Interpretation[Tooltip[Style[Subscript["\[ScriptCapitalE]", #1], Bold], StringTemplate["Eigenvalue ``"][First @ #2]], {#1, #2}] &,
                eigenvalues
            ],
            eigenvectors
        ];

        outputBasis = QuantumPartialTrace[qmo["Output"], Catenate @ Position[qmo["FullOutputOrder"], Alternatives @@ qmo["Target"]]];
        inputBasis = QuantumPartialTrace[qmo["Input"], Catenate @ Position[qmo["FullInputOrder"], Alternatives @@ qmo["Target"]]];

        (* construct *)
        operator = QuantumOperator[
            SparseArray @ Map[kroneckerProduct[IdentityMatrix[traceDimension, SparseArray], #] &, projectors],
            QuantumBasis[
                "Output" -> QuantumTensorProduct[
                    eigenBasis,
                    QuditBasis[outputBasis["Dimensions"]],
                    QuditBasis[tracedOperator["OutputDimensions"]]
                ],
                "Input" -> QuditBasis @ Join[inputBasis["Dimensions"], tracedOperator["InputDimensions"]]
            ]
        ];

        (* change back basis *)
        operator = profile["basis change"] @ QuantumOperator[
            operator,
            {{0, 1}, qmo["InputOrder"]},
            QuantumBasis[
                "Output" -> QuantumTensorProduct[
                    eigenBasis,
                    outputBasis,
                    tracedOperator["Output"]
                ],
                "Input" -> QuantumTensorProduct[inputBasis, tracedOperator["Input"]],
                "Label" -> "Eigen"[qmo["Label"]],
                "ParameterSpec" -> qmo["ParameterSpec"]
            ]
        ];

        (* permute and set order *)
        Simplify @ QuantumOperator[
            operator[
                "PermuteOutput", InversePermutation @ FindPermutation[Prepend[1 + Join[traceQudits, qmo["Target"] - Min[qmo["InputOrder"]] + 1], 1]]
            ][
                "PermuteInput", InversePermutation @ FindPermutation[Join[traceQudits, qmo["Target"] - Min[qmo["InputOrder"]] + 1]]
            ],
            {Prepend[Sort @ qmo["OutputOrder"], 0], Sort @ qmo["InputOrder"]}
        ]
    ]
]

QuantumMeasurementOperatorProp[qmo_, "POVM", args___] := QuantumMeasurementOperator[qmo["SuperOperator", args], qmo["Targets"]]

QuantumMeasurementOperatorProp[qmo_, "QASM"] := StringRiffle[MapIndexed[StringTemplate["c[``] = measure q[``];"][First[#2] - 1, #1 - 1] &, qmo["Target"]], "\n"]


QuantumMeasurementOperatorProp[qmo_, "Shift", n : _Integer : 1] :=
    QuantumMeasurementOperator[QuantumOperator[qmo]["Reorder", qmo["Order"] /. k_Integer ? Positive :> k + n], qmo["Target"] + n]

QuantumMeasurementOperatorProp[qmo_, "Bend", autoShift : _Integer ? Positive : Automatic] := With[{
    shift = Replace[autoShift, Automatic :> Max[qmo["Order"]]],
    target = qmo["Target"]
},
    If[ qmo["POVMQ"],
        QuantumMeasurementOperator[QuantumOperator[QuantumChannel[qmo]["Bend", shift]], Join[target, target - Min[target] + 1 + shift]],
        QuantumMeasurementOperator[QuantumOperator[qmo]["Bend", shift], Join[target, target - Min[target] + 1 + shift]]
    ]
]

QuantumMeasurementOperatorProp[qmo_, prop : "Conjugate" | "Dual" | "Unbend"] :=
    QuantumMeasurementOperator[qmo["SuperOperator"][prop], qmo["Target"]]


QuantumMeasurementOperatorProp[qmo_, "DiscardExtraQudits"] := QuantumOperator[
    Fold[
        #2[#1] &,
        qmo["SuperOperator"],
        With[{picture = qmo["Picture"]},
            MapThread[
                QuantumOperator[With[{d = Sqrt[#1["Dimension"]]},
                    If[picture === "PhaseSpace" && IntegerQ[d], "Double"["Trace"[d]], "Trace"[#1]]], {#2}
                ] &,
                {qmo["Eigenbasis"]["Decompose"], qmo["Eigenorder"]}
            ]
        ]
    ],
    qmo["InputOrder"] -> qmo["TargetOrder"],
    "Label" -> "Measurement"[qmo["Label"]],
    qmo["Basis"]["Options"]
]


QuantumMeasurementOperatorProp[qmo_, "CircuitDiagram", opts___] :=
    QuantumCircuitOperator[qmo]["Diagram", opts]


(* operator properties *)

(* Changing the underlying operator to the computational basis strips the eigenvalue-tagged
   pointer labels, which silently corrupts measurement statistics for non-Z bases.
   When the eigenbasis is not already computational, leave the measurement untouched. *)
QuantumMeasurementOperatorProp[qmo_, "Computational"] /; ! TrueQ[qmo["Eigenbasis"]["ComputationalQ"]] := qmo

QuantumMeasurementOperatorProp[qmo_, prop : "Ordered" | "SortOutput" | "SortInput" | "Computational" | "Simplify" | "FullSimplify" | "Chop" | "ComplexExpand" | "Reorder", args___] :=
    QuantumMeasurementOperator[qmo["QuantumOperator"][prop, args], qmo["Target"]]

QuantumMeasurementOperatorProp[qmo_, prop : "Dagger", args___] :=
    qmo["SuperOperator"][prop, args]

QuantumMeasurementOperatorProp[qmo_, prop : "Double", args___] :=
    QuantumMeasurementOperator[qmo["SuperOperator"][prop, args], qmo["Target"]]


QuantumMeasurementOperatorProp[qmo_, args : PatternSequence[prop_String, ___]] /;
    MemberQ[Intersection[qmo["Operator"]["AllProperties"], qmo["AllProperties"]], prop] := qmo["Operator"][args]

