Package["Wolfram`QuantumFramework`"]

PackageExport["QuantumState"]

PackageScope["QuantumStateQ"]


QuantumState::inState = "is invalid";

QuantumState::inBasis = "has invalid basis";

QuantumState::incompatible = "is incompatible with its basis";


QuantumStateQ[QuantumState[state_, basis_]] :=
    (stateQ[state] || (Message[QuantumState::inState]; False)) &&
    (QuantumBasisQ[basis] || (Message[QuantumState::inBasis]; False)) &&
    (Length[state] === basis["Dimension"] || (Message[QuantumState::incompatible]; False))

QuantumStateQ[___] := False


(* basis argument input *)

QuantumState[state_ ? stateQ, basisArgs___] /; !QuantumBasisQ[basisArgs] := Enclose @ Module[{
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


(* expand basis *)

QuantumState[state : Except[_ ? QuantumStateQ], args : Except[_ ? QuantumBasisQ]] :=
    Enclose @ QuantumState[state, ConfirmBy[QuantumBasis[args], QuantumBasisQ]]

QuantumState[state_ ? stateQ, basis_ ? QuantumBasisQ] := QuantumState[
    state,
    QuantumTensorProduct[basis, QuantumBasis[Max[2, Length[state] - basis["Dimension"]]]]
] /; Length[state] > basis["Dimension"]


(* pad state *)

QuantumState[state_ ? stateQ, basis_ ? QuantumBasisQ] := QuantumState[
    PadRight[state, Table[basis["Dimension"], TensorRank[state]]],
    basis
] /; Length[state] < basis["Dimension"]


(* Mutation *)

QuantumState[state_ ? stateQ, basis_ ? QuantumBasisQ] /; !SparseArrayQ[state] :=
    QuantumState[SparseArray[state], basis]


QuantumState[qs_ ? QuantumStateQ, args : Except[_ ? QuantumBasisQ, Except[Alternatives @@ $QuantumBasisPictures, _ ? nameQ]]] :=
    Enclose @ QuantumState[qs, ConfirmBy[QuantumBasis[args], QuantumBasisQ]]

QuantumState[qs_ ? QuantumStateQ, args : PatternSequence[Except[_ ? QuantumBasisQ], ___]] :=
    Enclose @ QuantumState[qs, ConfirmBy[QuantumBasis[qs["Basis"], args], QuantumBasisQ]]


(* change of basis *)

QuantumState[qs_ ? QuantumStateQ, newBasis_ ? QuantumBasisQ] /; ! newBasis["SortedQ"] := QuantumState[qs, newBasis["Sort"]]

QuantumState[qs_ ? QuantumStateQ, newBasis_ ? QuantumBasisQ] /; qs["Basis"] == newBasis := QuantumState[qs["State"], newBasis]

QuantumState[qs_ ? QuantumStateQ, newBasis_ ? QuantumBasisQ] /; qs["ElementDimension"] == newBasis["ElementDimension"] := Switch[
    qs["StateType"],
    "Vector",
    QuantumState[
        Simplify @ Flatten[
            PseudoInverse[newBasis["OutputMatrix"]] . (qs["OutputMatrix"] . qs["StateMatrix"] . PseudoInverse[qs["InputMatrix"]]) . newBasis["InputMatrix"]
            (* PseudoInverse[qs["OutputMatrix"]] . qs["StateMatrix"] . newBasis["OutputMatrix"] *)
        ],
        newBasis
    ],
    "Matrix",
    QuantumState[
        PseudoInverse[newBasis["Matrix"]] . (qs["Basis"]["Matrix"] . qs["DensityMatrix"] . PseudoInverse[qs["Basis"]["Matrix"]]) . newBasis["Matrix"],
        newBasis
    ]
]


(* renew basis *)

QuantumState[qs_ ? QuantumStateQ] := qs["Computational"]

(*QuantumState[qs_ ? QuantumStateQ, args__] := With[{
    newBasis = QuantumBasis[qs["Basis"], args]},
    If[ qs["Basis"] === newBasis,
        qs,
        QuantumState[qs["State"], newBasis]
    ]
]
*)

(* equality *)

QuantumState /: Equal[qs : _QuantumState ...] :=
    Equal @@ (#["Picture"] & /@ {qs}) && Equal @@ (#["NormalizedMatrixRepresentation"] & /@ {qs})

(* addition *)

addQuantumStates[qs1_QuantumState ? QuantumStateQ, qs2_QuantumState ? QuantumStateQ] /; qs1["Dimension"] == qs2["Dimension"] :=
    QuantumState[
        QuantumState[
            If[ qs1["StateType"] === qs2["StateType"] === "Vector",
                qs1["VectorRepresentation"] + qs2["VectorRepresentation"],
                qs1["MatrixRepresentation"] + qs2["MatrixRepresentation"]
            ],
            QuantumBasis[qs1["Dimensions"]]
        ],
        qs1["Basis"]
    ]

QuantumState /: Plus[states : _QuantumState...] := Fold[addQuantumStates, {states}]


(* multiplication *)

QuantumState /: (qs1_QuantumState ? QuantumStateQ) * (qs2_QuantumState ? QuantumStateQ) /; qs1["Dimension"] == qs2["Dimension"] :=
    QuantumState[
        QuantumState[
            If[ qs1["StateType"] === qs2["StateType"] === "Vector",
                qs1["VectorRepresentation"] * qs2["VectorRepresentation"],
                qs1["MatrixRepresentation"] * ArrayReshape[qs2["MatrixRepresentation"], Dimensions @ qs1["MatrixRepresentation"]]
            ],
            QuantumBasis[qs1["Dimensions"]]
        ],
        qs1["Basis"]
    ]

QuantumState /: (x : (_ ? NumericQ) | _Symbol) * (qs_QuantumState ? QuantumStateQ) :=
    QuantumState[
        x qs["State"],
        qs["Basis"]
    ]

QuantumState /: f_Symbol[left : _ ? NumericQ ..., qs_QuantumState, right : _ ? NumericQ ...] /; MemberQ[Attributes[f], NumericFunction] :=
    Enclose @ QuantumState[ConfirmBy[MatrixFunction[f[left, #, right] &, qs["DensityMatrix"]], MatrixQ], qs["Basis"]]


(* composition *)


(qs1_QuantumState ? QuantumStateQ)[(qs2_QuantumState ? QuantumStateQ)] /; qs1["Input"] == qs2["Output"] := QuantumState[
    Simplify @ Flatten[qs1["Pure"]["StateMatrix"] . qs2["Pure"]["StateMatrix"]],
    QuantumBasis[
        "Output" -> qs1["Pure"]["Output"], "Input" -> qs2["Pure"]["Input"],
        "Label" -> qs1["Label"] @* qs2["Label"],
        "ParameterSpec" -> Join[qs2["ParameterSpec"], qs1["ParameterSpec"]]
    ]
]

(qs1_QuantumState ? QuantumStateQ)[(qs2_QuantumState ? QuantumStateQ)] /; qs1["InputDimension"] == qs2["OutputDimension"] :=
    QuantumState[
        QuantumState[
            Simplify @ Flatten[qs1["PureMatrix"] . qs2["PureMatrix"]],
            QuantumBasis["Output" -> QuditBasis[qs1["OutputDimensions"]], "Input" -> QuditBasis[qs2["InputDimensions"]]]
        ],
        QuantumBasis[
            "Output" -> qs1["Output"], "Input" -> qs2["Input"],
            "Label" -> qs1["Label"] @* qs2["Label"],
            "ParameterSpec" -> Join[qs2["ParameterSpec"], qs1["ParameterSpec"]]
        ]
    ]


(* parameterization *)

(qs_QuantumState ? QuantumStateQ)[ps___] /; Length[{ps}] <= qs["ParameterArity"] :=
    QuantumState[
        Map[ReplaceAll[Thread[Take[qs["Parameters"], Length[{ps}]] -> {ps}]], qs["State"], {-1}],
        QuantumBasis[qs["Basis"], "ParameterSpec" -> Drop[qs["ParameterSpec"], Length[{ps}]]]
    ]
