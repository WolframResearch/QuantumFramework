Package["Wolfram`QuantumFramework`"]



$QuantumMeasurementOperatorProperties = {
    "QuantumOperator", "Target", "Targets",
    "TargetIndex",
    "Operator", "Basis", "MatrixRepresentation", "POVMElements",
    "OrderedMatrixRepresentation", "OrderedPOVMElements",
    "Arity", "Eigenqudits", "Dimensions", "Order", "HermitianQ", "UnitaryQ", "Eigenvalues", "Eigenvectors",
    "ProjectionQ", "POVMQ",
    "SuperOperator", "POVM"
};


QuantumMeasurementOperator["Properties"] := Union @ Join[
    $QuantumMeasurementOperatorProperties,
    $QuantumOperatorProperties
]

qmo_QuantumMeasurementOperator["ValidQ"] := QuantumMeasurementOperatorQ[qmo]


QuantumMeasurementOperator::undefprop = "QuantumMeasurementOperator property `` is undefined for this operator";


(qmo_QuantumMeasurementOperator[prop_ ? propQ, args___]) /; QuantumMeasurementOperatorQ[qmo] := With[{
    result = QuantumMeasurementOperatorProp[qmo, prop, args]
    },
    If[TrueQ[$QuantumFrameworkPropCache], QuantumMeasurementOperatorProp[qmo, prop, args] = result, result]
        /; !FailureQ[Unevaluated @ result] && (!MatchQ[result, _QuantumMeasurementOperatorProp] || Message[QuantumMeasurementOperator::undefprop, prop])
]

QuantumMeasurementOperatorProp[qmo_, "Properties"] :=
    DeleteDuplicates @ Join[QuantumMeasurementOperator["Properties"], qmo["Operator"]["Properties"]]


(* getters *)

QuantumMeasurementOperatorProp[_[op_, _], "Operator" | "QuantumOperator"] := op

QuantumMeasurementOperatorProp[_[_, target_], "Target"] := target

QuantumMeasurementOperatorProp[qmo_, "Targets"] := Length[qmo["Target"]]

QuantumMeasurementOperatorProp[qmo_, "TargetIndex"] := Catenate @ Lookup[PositionIndex[qmo["FullInputOrder"]], qmo["Target"]]

QuantumMeasurementOperatorProp[qmo_, "TargetDimensions"] :=
    Extract[qmo["InputDimensions"], Position[qmo["FullInputOrder"], Alternatives @@ qmo["Target"]] ]

QuantumMeasurementOperatorProp[qmo_, "TargetDimension"] := Times @@ qmo["TargetDimensions"]

QuantumMeasurementOperatorProp[qmo_, "Eigenqudits"] := If[qmo["POVMQ"],
    Enclose[
        First @ Confirm @ FirstPosition[
            FoldList[Times, qmo["OutputDimensions"]],
            qmo["TargetDimension"]
        ],
        1 &
    ],
    1
]

QuantumMeasurementOperatorProp[qmo_, "Type"] := Which[
    qmo["OutputDimension"] == qmo["InputDimension"],
    "Projection",
    qmo["OutputDimension"] > qmo["InputDimension"],
    "POVM",
    True,
    "Destructive"
]

QuantumMeasurementOperatorProp[qmo_, "ProjectionQ"] := qmo["Type"] === "Projection"

QuantumMeasurementOperatorProp[qmo_, "POVMQ"] := qmo["Type"] === "POVM"

QuantumMeasurementOperatorProp[qmo_, "POVMElements"] := If[qmo["POVMQ"], qmo["Tensor"], qmo["Projectors"]]

QuantumMeasurementOperatorProp[qmo_, "OrderedPOVMElements"] := If[qmo["POVMQ"],
    qmo["OrderedTensor"],
    projector /@ qmo["OrderedMatrix"]
]

QuantumMeasurementOperatorProp[qmo_, "Operators"] := If[qmo["POVMQ"],
    AssociationThread[Range[0, Length[qmo["Tensor"]] - 1], QuantumOperator[#, QuantumBasis["Output" -> qmo["Basis"]["Input"]], qmo["Order"]] & /@ qmo["Tensor"]],
    AssociationThread[Eigenvalues[qmo["Matrix"]], QuantumOperator[projector @ #, qmo["Basis"], qmo["Order"]] & /@ Eigenvectors[qmo["OrderedMatrix"]]]
]

QuantumMeasurementOperatorProp[qmo_, "SuperOperator"] := Module[{
    trace,
    traceQudits,
    tracedOperator,
    eigenvalues, eigenvectors, projectors,
    eigenBasis, outputBasis, inputBasis, operator
},
    trace = Complement[qmo["InputOrder"], qmo["Target"]];
    traceQudits = trace - Min[qmo["InputOrder"]] + 1;

    If[
        qmo["POVMQ"],

        qmo["Operator"],

        tracedOperator = QuantumPartialTrace[
            qmo,
            If[qmo["POVMQ"], {# + qmo["OutputQudits"] - qmo["InputQudits"], #} & /@ trace, trace]
        ];

        {eigenvalues, eigenvectors} = profile["Eigensystem"] @ tracedOperator["Eigensystem", "Sort" -> True];
        projectors = projector /@ eigenvectors;

        eigenBasis = QuditBasis[
            MapIndexed[
                Interpretation[Tooltip[Style[#, Bold], StringTemplate["Eigenvalue ``"][First @ #2]], {#1, #2}] &,
                eigenvalues
            ],
            eigenvectors
        ];

        outputBasis = QuantumPartialTrace[qmo["Output"], Catenate @ Position[qmo["OutputOrder"], Alternatives @@ qmo["Target"]]];
        inputBasis = QuantumPartialTrace[qmo["Input"], Catenate @ Position[qmo["InputOrder"], Alternatives @@ qmo["Target"]]];

        (* construct *)
        operator = QuantumOperator[
            SparseArray @ Map[kroneckerProduct[IdentityMatrix[Times @@ qmo["InputDimensions"][[traceQudits]], SparseArray], #] &, projectors],

            QuantumBasis[
                "Output" -> QuantumTensorProduct[
                    eigenBasis,
                    QuditBasis[outputBasis["Dimensions"]],
                    QuditBasis[tracedOperator["OutputDimensions"]]
                ],
                "Input" -> QuantumTensorProduct[QuditBasis[inputBasis["Dimensions"]], QuditBasis[tracedOperator["InputDimensions"]]]
            ]
        ];

        (* change back basis *)
        operator = profile["basis change"] @ QuantumOperator[
            operator,
            QuantumBasis[
                "Output" -> QuantumTensorProduct[
                    eigenBasis,
                    outputBasis,
                    tracedOperator["Output"]
                ],
                "Input" -> QuantumTensorProduct[inputBasis, tracedOperator["Input"]],
                "Label" -> qmo["Label"],
                "ParameterSpec" -> qmo["ParameterSpec"]
            ]
        ];

        (* permute and set order *)
        QuantumOperator[
            operator[
                {"PermuteOutput", InversePermutation @ FindPermutation[Prepend[1 + Join[traceQudits, qmo["Target"] - Min[qmo["InputOrder"]] + 1], 1]]}
            ][
                {"PermuteInput", InversePermutation @ FindPermutation[Join[traceQudits, qmo["Target"] - Min[qmo["InputOrder"]] + 1]]}
            ],
            {Automatic, qmo["InputOrder"]}
        ]
    ]
]

QuantumMeasurementOperatorProp[qmo_, "POVM"] := QuantumMeasurementOperator[qmo["SuperOperator"], qmo["Target"]]

(* operator properties *)

QuantumMeasurementOperatorProp[qmo_, prop : "Ordered" | {"Ordered", __} | "Sort" | "Computational"] :=
    QuantumMeasurementOperator[qmo["QuantumOperator"][prop], qmo["Target"]]

QuantumMeasurementOperatorProp[qmo_, args : PatternSequence[prop_String, ___] | PatternSequence[{prop_String, ___}, ___]] /;
    MemberQ[Intersection[qmo["Operator"]["Properties"], qmo["Properties"]], prop] := qmo["Operator"][args]

