Package["Wolfram`QuantumFramework`"]

PackageExport["QuantumState"]

PackageScope["quantumStateQ"]
PackageScope["QuantumStateQ"]
PackageScope["addQuantumStates"]


QuantumState::invalidName = "`1` is not a recognized QuantumState constructor"


quantumStateQ[QuantumState[state_SparseArray ? stateQ, qb_QuantumBasis /; QuantumBasisQ[Unevaluated[qb]]]] := Length[state] === qb["Dimension"]

quantumStateQ[___] := False


QuantumStateQ[qs_QuantumState] := System`Private`HoldValidQ[qs] || quantumStateQ[Unevaluated[qs]]

QuantumStateQ[___] := False


qs_QuantumState /; System`Private`HoldNotValidQ[qs] && quantumStateQ[Unevaluated[qs]] := System`Private`HoldSetValid[qs]


(* basis argument input *)

QuantumState[state_ ? stateQ] := QuantumState[state, QuantumBasis[primeFactors[Length[state]]]]

QuantumState[state_ ? stateQ, basisArgs__] /; ! QuantumBasisQ[basisArgs] := Enclose @ Block[{
    basis, multiplicity
},
    basis = ConfirmBy[QuantumBasis[basisArgs], QuantumBasisQ];
    multiplicity = basisMultiplicity[Length[state], basis["Dimension"]];
    basis = ConfirmBy[QuantumBasis[basis, multiplicity], QuantumBasisQ];
    QuantumState[
        PadRight[state, Table[basis["Dimension"], TensorRank[state]]],
        basis
    ]
]


(* association input *)

QuantumState[state_ ? AssociationQ, basisArgs___] /; VectorQ[Values[state]] := Enclose @ Module[{
    basis = ConfirmBy[QuantumBasis[basisArgs], QuantumBasisQ], multiplicity},
    multiplicity = basisMultiplicity[Length[state], basis["Dimension"]];
    basis = ConfirmBy[QuantumBasis[basis, multiplicity], QuantumBasisQ];
    ConfirmAssert[ContainsOnly[QuditName /@ Keys[state], basis["ElementNames"]], "Association keys and basis names don't match"];
    QuantumState[
        Values @ KeyMap[QuditName, state][[Key /@ basis["ElementNames"]]] /. _Missing -> 0,
        basis
    ]
]


(* eigenvalues input *)

QuantumState["Eigenvalues" -> eigenvalues_ ? VectorQ, basisArgs___] := With[{
    basis = QuantumBasis[basisArgs]
},
    QuantumState[
        Total @ MapThread[#1 #2 &, {eigenvalues, basis["Projectors"]}],
        basis
    ] /; Length[eigenvalues] == basis["Dimension"]
]

(* conversion *)

QuantumState[obj : _QuantumOperator | _QuantumMeasurementOperator | _QuantumMeasurement | _QuantumChannel | _QuantumCircuitOperator, opts___] :=
    QuantumState[obj["State"], opts]


(* active basis transform *)

QuantumState[qb_QuditBasis, opts___] := QuantumState[QuantumBasis[qb], opts]

QuantumState[qb_ ? QuantumBasisQ, opts___] := Enclose @ If[
    qb["Picture"] === "PhaseSpace",

    With[{dims = Sqrt[qb["Dimensions"]], n = qb["Qudits"]},
        ConfirmAssert[AllTrue[dims, IntegerQ]];
        QuantumState[
            ArrayReshape[
                Transpose[
                    ArrayReshape[Inverse[qb["Matrix"]], Join[#, #] & @ Catenate[{#, #} & /@ dims]],
                    If[ n > 1,
                        PermutationProduct[
                            Cycles[NestList[# + 2 &, {2, 2 n + 1 }, n - 1]],
                            Cycles[NestList[# + 4 &, {2, 3}, n - 1]]
                        ],
                        Cycles[{{2, 3}}]
                    ]
                ],
                {#, #} & @ qb["Dimension"]
            ],
            QuantumBasis[dims, dims, opts, qb["Options"]]
        ]
    ],

    QuantumState[Flatten[Inverse[qb["Matrix"]]], QuantumBasis[qb["QuditBasis"], qb["Dimensions"], opts, qb["Options"]]]
]


(* number *)

QuantumState[x : Except[_ ? QuantumStateQ | _ ? stateQ], basisArgs___] := QuantumState[{x}, QuantumBasis[basisArgs]]


(* expand basis *)

QuantumState[state : Except[_ ? QuantumStateQ], args : Except[_ ? QuantumBasisQ]] :=
    Enclose @ QuantumState[state, ConfirmBy[QuantumBasis[args], QuantumBasisQ]]

QuantumState[state_ ? stateQ, basis_ ? QuantumBasisQ] := QuantumState[
    state,
    QuantumTensorProduct[basis, QuantumBasis[Max[2, Length[state] - basis["Dimension"]]]]
] /; Length[state] > basis["Dimension"] > 0


(* pad state *)

QuantumState[state_ ? stateQ, basis_ ? QuantumBasisQ] := QuantumState[
    PadRight[state, Table[basis["Dimension"], TensorRank[state]]],
    basis
] /; Length[state] < basis["Dimension"]


(* Mutation *)

QuantumState[state_ ? stateQ, basis_ ? QuantumBasisQ] /; !SparseArrayQ[state] && state =!= {} :=
    Enclose @ QuantumState[ConfirmBy[SparseArray[state, Length[state]], SparseArrayQ], basis]

QuantumState[qs_ ? QuantumStateQ, args : Except[_ ? QuantumBasisQ, Except[Alternatives @@ $QuantumBasisPictures, _ ? nameQ | _Integer]]] :=
    Enclose @ QuantumState[qs, ConfirmBy[QuantumBasis[args], QuantumBasisQ]]

QuantumState[qs_ ? QuantumStateQ, args : PatternSequence[Except[_ ? QuantumBasisQ | _ ? QuantumStateQ], ___]] :=
    Enclose @ QuantumState[qs, ConfirmBy[QuantumBasis[qs["Basis"], args], QuantumBasisQ]]


(* change of basis *)

QuantumState[qs_ ? QuantumStateQ, newBasis_ ? QuantumBasisQ] /; ! newBasis["SortedQ"] := QuantumState[qs, newBasis["Sort"]]

QuantumState[qs_ ? QuantumStateQ, newBasis_ ? QuantumBasisQ] /; qs["Basis"] == newBasis || qs["ComputationalQ"] && newBasis["ComputationalQ"] := QuantumState[qs["State"], newBasis]

QuantumState[qs_ ? QuantumStateQ, newBasis_ ? QuantumBasisQ] /; qs["Dimension"] == newBasis["Dimension"] :=
    Enclose[Which[
        qs["Dimension"] == 0,
        QuantumState[qs["State"], newBasis],
        qs["VectorQ"],
        QuantumState[
            SparseArrayFlatten @ ConfirmQuiet[
                Dot[
                    MatrixInverse[newBasis["Output"]["ReducedMatrix"]],
                    qs["Output"]["ReducedMatrix"] . qs["StateMatrix"] . MatrixInverse[qs["Input"]["ReducedMatrix"]],
                    newBasis["Input"]["ReducedMatrix"]
                ],
                Dot::dotsh
            ],
            newBasis
        ],
        qs["MatrixQ"],
        QuantumState[
            ConfirmQuiet[
                Dot[
                    MatrixInverse[newBasis["ReducedMatrix"]],
                    qs["Basis"]["ReducedMatrix"] . qs["DensityMatrix"] . MatrixInverse[qs["Basis"]["ReducedMatrix"]],
                    newBasis["ReducedMatrix"]
                ],
                Dot::dotsh
            ],
            newBasis
        ],
        True,
        $Failed
    ],
    (ReleaseHold[#["HeldMessageCall"]]; Throw[$Failed]) &
]

QuantumState[qs_ ? QuantumStateQ, newBasis_ ? QuantumBasisQ] := Switch[
    qs["StateType"],
    "Vector",
    QuantumState[PadRight[qs["State"], newBasis["Dimension"]], newBasis],
    "Matrix",
    QuantumState[PadRight[qs["State"], {newBasis["Dimension"], newBasis["Dimension"]}], newBasis]
]

QuantumState[qs_ ? QuantumStateQ, newBasis_ ? QuantumBasisQ, opts__] :=
    QuantumState[qs, QuantumBasis[newBasis, opts]]


(* equality *)

QuantumState /: Equal[qs__QuantumState] :=
    Which[
        And @@ (#["VectorStateQ"] & /@ {qs}),
        Thread[Equal @@ (Chop @ SetPrecisionNumeric[#["Computational"]["CanonicalStateVector"]] & /@ {qs})],
        And @@ (#["MatrixStateQ"] & /@ {qs}),
        Thread[Equal @@ (Chop @ SetPrecisionNumeric[SparseArrayFlatten @ #["NormalizedMatrixRepresentation"]] & /@ {qs})],
        True,
        Thread[Equal @@ Through[{qs}["MatrixState"]]]
    ]

QuantumState /: Unequal[qs__QuantumState] := ! Equal[qs]


(* numeric function *)

QuantumState /: f_Symbol[left : Except[_QuantumState] ..., qs_QuantumState, right : Except[_QuantumState | OptionsPattern[]] ..., opts : OptionsPattern[]] /; MemberQ[Attributes[f], NumericFunction] :=
    Enclose @ QuantumState[
        ConfirmBy[If[MatchQ[f, Times | Conjugate] && qs["VectorQ"], f[left, qs["StateVector"], right], matrixFunction[f, qs["DensityMatrix"], {left}, {right}, opts]], stateQ],
        qs["Basis"]
    ]

(* Trace *)

QuantumState /: Tr[qs_QuantumState] := Tr @ qs["DensityMatrix"]


(* addition *)

addQuantumStates[qs1_QuantumState ? QuantumStateQ, qs2_QuantumState ? QuantumStateQ] /; qs1["Dimension"] == qs2["Dimension"] :=
    QuantumState[
        If[ qs1["StateType"] === qs2["StateType"] === "Vector",
            qs1["StateVector"] + QuantumState[qs2, qs1["Basis"]]["StateVector"],
            qs1["DensityMatrix"] + QuantumState[qs2, qs1["Basis"]]["DensityMatrix"]
        ],
        qs1["Basis"],
        "ParameterSpec" -> MergeParameterSpecs[qs1, qs2]
    ]

QuantumState /: HoldPattern[Plus[states__QuantumState]] /; Length[{states}] > 1 :=
    If[ Equal @@ (#["Dimension"] & /@ {states}), Fold[addQuantumStates, {states}],
        Failure["QuantumState", <|"MessageTemplate" -> "Incompatible dimensions"|>]
    ]


(* multiplication *)

multiplyQuantumStates[qs1_QuantumState, qs2_QuantumState] /; qs1["Dimension"] == qs2["Dimension"] :=
    QuantumState[
        QuantumState[
            If[ qs1["StateType"] === qs2["StateType"] === "Vector",
                qs1["VectorRepresentation"] * qs2["VectorRepresentation"],
                qs1["MatrixRepresentation"] * ArrayReshape[qs2["MatrixRepresentation"], Dimensions @ qs1["MatrixRepresentation"]]
            ],
            QuantumBasis[qs1["Dimensions"]]
        ],
        qs1["Basis"],
        "ParameterSpec" -> MergeParameterSpecs[qs1, qs2]
    ]

QuantumState /: HoldPattern[Times[states : _QuantumState ? QuantumStateQ ...]] :=
    If[ Equal @@ (#["Dimension"] & /@ {states}), Fold[multiplyQuantumStates, {states}],
        Failure["QuantumState", <|"MessageTemplate" -> "Incompatible dimensions"|>]
    ]


(* differentiation *)

QuantumState /: D[qs : _QuantumState, args___] := QuantumState[D[qs["State"], args], qs["Basis"]]


(* duals *)

SuperDagger[qs_QuantumState] ^:= qs["Dagger"]

SuperStar[qs_QuantumState] ^:= qs["Conjugate"]

Transpose[qs_QuantumState, args___] ^:= qs["Transpose", args]

Inverse[qs_QuantumState] ^:= qs ^ -1


(* simplify *)

Scan[
    (Symbol[#][qs_QuantumState, args___] ^:= qs[#, args]) &,
    {"Simplify", "FullSimplify", "Chop", "ComplexExpand"}
]


(* join *)

QuantumState[qs_QuantumState ? QuantumStateQ] := qs

QuantumState[qs__QuantumState ? QuantumStateQ] := QuantumState[
    If[ And @@ (#["PureStateQ"] & /@ {qs}),
        SparseArrayFlatten[blockDiagonalMatrix[#["StateMatrix"] & /@ {qs}]],
        blockDiagonalMatrix[#["DensityMatrix"] & /@ {qs}]
    ],
    Plus @@ (#["Basis"] & /@ {qs}),
    "Label" -> CirclePlus @@ (#["Label"] & /@ {qs}),
    "ParameterSpec" -> MergeParameterSpecs[qs1, qs2]
]


(* composition *)


(top_QuantumState ? QuantumStateQ)[(bot_QuantumState ? QuantumStateQ)] := (QuantumOperator[top] @ QuantumOperator[bot])["Sort"]["State"]


(qs1_QuantumState ? QuantumStateQ)[(qs2_QuantumState ? QuantumStateQ)] /; (
    qs1["InputQudits"] >= 1 && qs2["OutputQudits"] >= 1 && (
        (qs1["InputQudits"] <= qs2["OutputQudits"] && Take[qs2["OutputDimensions"], qs1["InputQudits"]] === qs1["InputDimensions"]) ||
        (qs2["OutputQudits"] <= qs1["InputQudits"] && Take[qs1["InputDimensions"], qs2["OutputQudits"]] === qs2["OutputDimensions"])
    ) && (
        qs1["OutputDimension"] * qs2["InputDimension"] == 0 ||
        (
            qs1["VectorQ"] && qs1["Picture"] === "Schrodinger" && qs2["Picture"] === "Schrodinger" &&
            (qs2["VectorQ"] || qs1["InputQudits"] <= qs2["OutputQudits"])
        ) || qs1["InputQudits"] == qs2["OutputQudits"]
    )
) := Block[{
    kMin = Min[qs1["InputQudits"], qs2["OutputQudits"]],
    inDim = qs2["InputDimension"],
    restOut, restOutDim, restIn, restInDim,
    q1, q2, q1Mat, q1ConjMat, liftedQ1, liftedQ1Dag, q2Mat, liftedQ2, state, outBasis, inBasis
},
    restOut = qs2["Output"]["Decompose"][[kMin + 1 ;;]];
    restOutDim = Times @@ Drop[qs2["OutputDimensions"], kMin];
    restIn = qs1["Input"]["Decompose"][[kMin + 1 ;;]];
    restInDim = Times @@ Drop[qs1["InputDimensions"], kMin];
    Which[
        qs1["OutputDimension"] * inDim == 0,
            QuantumState[{},
                QuantumBasis[
                    "Output" -> qs1["Output"], "Input" -> qs2["Input"],
                    "Label" -> qs1["Label"] @* qs2["Label"],
                    "Picture" -> If[MemberQ[{qs1["Picture"], qs2["Picture"]}, "PhaseSpace"], "PhaseSpace", qs1["Picture"]],
                    "ParameterSpec" -> MergeParameterSpecs[qs1, qs2]
                ]
            ],
        qs1["VectorQ"] && qs1["Picture"] === "Schrodinger" && qs2["Picture"] === "Schrodinger",
            q1 = qs1["Computational"]; q2 = qs2["Computational"];
            q1Mat = q1["Matrix"];
            liftedQ1 = If[restOutDim == 1, q1Mat, KroneckerProduct[q1Mat, IdentityMatrix[restOutDim, SparseArray]]];
            state = If[TrueQ[q2["VectorQ"]],
                q2Mat = q2["Matrix"];
                liftedQ2 = If[restInDim == 1, q2Mat, KroneckerProduct[q2Mat, IdentityMatrix[restInDim, SparseArray]]];
                SparseArrayFlatten[liftedQ1 . liftedQ2],
                q1ConjMat = ConjugateTranspose[q1Mat];
                liftedQ1Dag = If[restOutDim == 1, q1ConjMat, KroneckerProduct[q1ConjMat, IdentityMatrix[restOutDim, SparseArray]]];
                If[inDim == 1,
                    liftedQ1 . q2["DensityMatrix"] . liftedQ1Dag,
                    KroneckerProduct[liftedQ1, IdentityMatrix[inDim, SparseArray]] . q2["DensityMatrix"] . KroneckerProduct[liftedQ1Dag, IdentityMatrix[inDim, SparseArray]]
                ]
            ];
            outBasis = Which[
                qs1["OutputDimension"] == 1 && restOut === {}, qs1["Output"],
                qs1["OutputDimension"] == 1, QuantumTensorProduct @ restOut,
                restOut === {}, qs1["Output"],
                True, QuantumTensorProduct @ Join[qs1["Output"]["Decompose"], restOut]
            ];
            inBasis = Which[
                inDim == 1 && restIn === {}, qs2["Input"],
                inDim == 1, QuantumTensorProduct @ restIn,
                restIn === {}, qs2["Input"],
                True, QuantumTensorProduct @ Join[qs2["Input"]["Decompose"], restIn]
            ];
            QuantumState[state,
                QuantumBasis[
                    "Output" -> outBasis, "Input" -> inBasis,
                    "Label" -> qs1["Label"] @* qs2["Label"],
                    "ParameterSpec" -> MergeParameterSpecs[qs1, qs2]
                ]
            ],
        True,
            q1 = qs1["Computational"]["Double"]; q2 = qs2["Computational"]["Double"];
            state = ArrayReshape[
                q1["StateMatrix"] . q2["StateMatrix"],
                Table[qs1["OutputDimension"] inDim, 2]
            ];
            With[{
                s = QuantumState[state, "Output" -> QuditBasis @ qs1["OutputDimensions"], "Input" -> QuditBasis @ qs2["InputDimensions"]],
                b = QuantumBasis[
                    "Output" -> qs1["Output"], "Input" -> qs2["Input"],
                    "Label" -> qs1["Label"] @* qs2["Label"],
                    "Picture" -> If[MemberQ[{qs1["Picture"], qs2["Picture"]}, "PhaseSpace"], "PhaseSpace", qs1["Picture"]],
                    "ParameterSpec" -> MergeParameterSpecs[qs1, qs2]
                ]
            },
                QuantumState[s, b]
            ]
    ]
]


(* parameterization *)

(qs_QuantumState ? QuantumStateQ)[ps : PatternSequence[p : Except[_Association], ___]] /; ! MemberQ[QuantumState["Properties"], p] && Length[{ps}] <= qs["ParameterArity"] :=
    qs[AssociationThread[Take[qs["Parameters"], UpTo[Length[{ps}]]], {ps}]]

(qs_QuantumState ? QuantumStateQ)[rules_ ? AssociationQ] /; ContainsOnly[Keys[rules], qs["Parameters"]] :=
    QuantumState[
        Map[ReplaceAll[rules], qs["State"], {If[qs["VectorQ"], 1, 2]}],
        qs["Basis"][rules]
    ]

