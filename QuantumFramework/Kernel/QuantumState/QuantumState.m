Package["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`Arrays`"]

PackageExport["QuantumState"]

PackageScope["quantumStateQ"]
PackageScope["QuantumStateQ"]
PackageScope["addQuantumStates"]
PackageScope["valueOutsideFacts"]
PackageScope["domainFailure"]
PackageScope["$checkValues"]
PackageScope["basisRecord"]
PackageScope["basisDeclared"]
PackageScope["valuesAgainst"]


QuantumState::invalidName = "`1` is not a recognized QuantumState constructor"
QuantumState::invalidArgs = "QuantumState constructor `1` did not match any rule"
QuantumState::padded = "state of length `1` is zero-padded to dimension `2` (`3` qudits); use a length that is a power of the qudit dimension, or \"Register\"[n] for an n-qudit register"


(* Validity is a shape contract, not a container contract: the amplitudes may
   be held in any array container whose leading dimension matches the basis. *)

quantumStateQ[QuantumState[state_ ? stateQ, qb_QuantumBasis /; QuantumBasisQ[Unevaluated[qb]]]] := First[ArrayDimensions[state]] === qb["Dimension"]

quantumStateQ[___] := False


QuantumStateQ[qs_QuantumState] := System`Private`HoldValidQ[qs] || quantumStateQ[Unevaluated[qs]]

QuantumStateQ[___] := False


qs_QuantumState /; System`Private`HoldNotValidQ[qs] && quantumStateQ[Unevaluated[qs]] := System`Private`HoldSetValid[qs]


(* basis argument input *)

QuantumState[state_ ? stateQ] := QuantumState[state, QuantumBasis[primeFactors[First @ ArrayDimensions[state]]]]

QuantumState[state_ ? stateQ, basisArgs__] /; ! QuantumBasisQ[basisArgs] := Enclose @ Block[{
    basis, multiplicity, dim = First @ ArrayDimensions[state]
},
    basis = ConfirmBy[QuantumBasis[basisArgs], QuantumBasisQ];
    multiplicity = basisMultiplicity[dim, basis["Dimension"]];
    basis = ConfirmBy[QuantumBasis[basis, multiplicity], QuantumBasisQ];
    (* a multi-qudit pad means the qudit count was inferred by rounding up,
       not given exactly; the single-qudit amplitude-prefix pad stays silent *)
    If[ multiplicity > 1 && dim =!= basis["Dimension"],
        Message[QuantumState::padded, dim, basis["Dimension"], basis["Qudits"]]
    ];
    QuantumState[
        If[dim === basis["Dimension"], state, padState[state, basis["Dimension"]]],
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


(* "SymbolDomains" -> Reals on a state or operator: once the amplitudes are known, every symbol in them gets
   the fact *)
QuantumState[state_ ? stateQ, basis : QuantumBasis[data_Association] ? QuantumBasisQ] /; KeyExistsQ[Lookup[data, "SymbolDomains", <||>], All] :=
    QuantumState[state, QuantumBasis[<|data, "SymbolDomains" -> allDomainsRecord[data["SymbolDomains"], amplitudeValues[state]]|>]]

(* An amplitude that needs a sign the facts do not give, Sqrt[1 - q] with q a declared parameter and nothing
   more, gets one: the range it leaves for its symbol joins the facts, with one message (assumedSignRecord).
   A symbol with any stated domain, the user's, a role's, or one given by "SymbolDomains" -> Reals, keeps it
   as given: an explicit Reals may mean values outside a channel's physical range on purpose. *)
QuantumState[state_ ? stateQ, basis : QuantumBasis[data_Association] ? QuantumBasisQ] /;
    data["ParameterSpec"] =!= {} && ! FreeQ[amplitudeValues[state], Power[_, _Rational] | _Log] :=
    With[{r = assumedSignRecord[Lookup[data, "SymbolDomains", <||>], Complement[basis["Parameters"], Keys[Lookup[data, "SymbolDomains", <||>]]], amplitudeValues[state]]},
        (
            Message[QuantumState::assumed, Replace[Last[r], (h : Greater | GreaterEqual | Less | LessEqual)[a_, 0] :> a, {1}], And @@ Last[r]];
            QuantumState[state, QuantumBasis[<|data, "SymbolDomains" -> First[r]|>]]
        ) /; Last[r] =!= {}
    ]

amplitudeValues[s_SparseArray] := s["NonzeroValues"]
amplitudeValues[s_List] := Flatten[s]
amplitudeValues[_] := {}


(* expand basis *)

QuantumState[state : Except[_ ? QuantumStateQ], args : Except[_ ? QuantumBasisQ]] :=
    Enclose @ QuantumState[state, ConfirmBy[QuantumBasis[args], QuantumBasisQ]]

QuantumState[state_ ? stateQ, basis_ ? QuantumBasisQ] := QuantumState[
    state,
    QuantumTensorProduct[basis, QuantumBasis[Max[2, First @ ArrayDimensions[state] - basis["Dimension"]]]]
] /; First[ArrayDimensions[state]] > basis["Dimension"] > 0


(* pad state *)

(* Zero-padding to the basis dimension is defined only where the zeros can be
   written: a lazy or symbolic container has no elements to pad, so an
   under-dimensioned one is left for validity to reject rather than being
   grown into something whose amplitudes do not exist. *)

padState[state_ ? ArrayExplicitQ, dim_Integer] := PadRight[state, Table[dim, ArrayRank[state]]]

QuantumState[state_ ? stateQ, basis_ ? QuantumBasisQ] := QuantumState[
    padState[state, basis["Dimension"]],
    basis
] /; ArrayExplicitQ[state] && First[ArrayDimensions[state]] < basis["Dimension"]


(* Mutation *)

(* A list of amplitudes is normalized to a SparseArray on intake; any other
   container is stored as it arrived. *)

QuantumState[state_List ? stateQ, basis_ ? QuantumBasisQ] /; state =!= {} :=
    Enclose @ QuantumState[ConfirmBy[SparseArray[state, Length[state]], SparseArrayQ], basis]

QuantumState[qs_ ? QuantumStateQ, args : Except[_ ? QuantumBasisQ, Except[Alternatives @@ $QuantumBasisPictures, _ ? nameQ | _Integer]]] :=
    Enclose @ QuantumState[qs, ConfirmBy[QuantumBasis[args], QuantumBasisQ]]

QuantumState[qs_ ? QuantumStateQ, args : PatternSequence[Except[_ ? QuantumBasisQ | _ ? QuantumStateQ], ___]] :=
    Enclose @ QuantumState[qs, ConfirmBy[QuantumBasis[qs["Basis"], args], QuantumBasisQ]]


(* change of basis *)

(* row-side transform of a doubled (density-matrix) state: the flattened (output, input)
   index transforms as Out (x) (In^-1)^T, the vec image of Out . S . In^-1, so a dual
   input qudit enters through the inverse transpose (conjugated for a unitary basis);
   with no input qudits this is the basis ReducedMatrix itself *)
doubledReducedMatrix[qb_] := If[ qb["InputDimension"] === 1,
    qb["ReducedMatrix"],
    KroneckerProduct[qb["Output"]["ReducedMatrix"], Transpose[MatrixInverse[qb["Input"]["ReducedMatrix"]]]]
]

(* The facts about the symbols do not depend on the basis, so a change of basis joins the facts of the old
   basis with those of the new one, as any merge does: a symbol named by both lies in both sets of domains,
   and one whose domains meet nowhere, or whose facts were taken away, keeps none. *)
QuantumState[qs_ ? QuantumStateQ, newBasis_ ? QuantumBasisQ] /; symbolDomainsOf[qs] =!= <||> && symbolDomainsOf[qs] =!= symbolDomainsOf[newBasis] := With[{
    domains = MergeSymbolDomains[qs["Basis"], newBasis]
},
    QuantumState[qs, QuantumBasis[newBasis, "SymbolDomains" -> domains]] /; ! sameSymbolDomainsQ[domains, symbolDomainsOf[newBasis]]
]

sameSymbolDomainsQ[r1_, r2_] := KeySort[Sort /@ r1] === KeySort[Sort /@ r2]

QuantumState[qs_ ? QuantumStateQ, newBasis_ ? QuantumBasisQ] /; ! newBasis["SortedQ"] := QuantumState[qs, newBasis["Sort"]]

QuantumState[qs_ ? QuantumStateQ, newBasis_ ? QuantumBasisQ] /; qs["Basis"] == newBasis || qs["ComputationalQ"] && newBasis["ComputationalQ"] := QuantumState[qs["State"], newBasis]

QuantumState[qs_ ? QuantumStateQ, newBasis_ ? QuantumBasisQ] /; qs["Dimension"] == newBasis["Dimension"] :=
    Enclose[Which[
        qs["Dimension"] == 0,
        QuantumState[qs["State"], newBasis],
        qs["VectorQ"],
        QuantumState[
            (* TODO: ConfirmQuiet hides Dot::dotsh from shape mismatches; precompute shapes via Dimensions and assert before Dot *)
            ArrayVector @ ConfirmQuiet[
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
        With[{oldMatrix = doubledReducedMatrix[qs["Basis"]], newMatrix = doubledReducedMatrix[newBasis]},
            QuantumState[
                (* TODO: ConfirmQuiet hides Dot::dotsh from shape mismatches; precompute shapes via Dimensions and assert before Dot *)
                ConfirmQuiet[
                    Dot[
                        MatrixInverse[newMatrix],
                        oldMatrix . qs["DensityMatrix"] . MatrixInverse[oldMatrix],
                        newMatrix
                    ],
                    Dot::dotsh
                ],
                newBasis
            ]
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
        Thread[Equal @@ (Chop @ SetPrecisionNumeric[ArrayVector @ #["NormalizedMatrixRepresentation"]] & /@ {qs})],
        True,
        Thread[Equal @@ Through[{qs}["MatrixState"]]]
    ]

QuantumState /: Unequal[qs__QuantumState] := ! Equal[qs]


(* numeric function *)

QuantumState /: (f_Symbol ? numericFunctionQ)[left : Except[_QuantumState] ..., qs_QuantumState, right : Except[_QuantumState | OptionsPattern[]] ..., opts : OptionsPattern[]] :=
    Enclose @ QuantumState[
        ConfirmBy[If[MatchQ[f, Times | Conjugate] && qs["VectorQ"], f[left, qs["StateVector"], right], matrixFunction[f, qs["DensityMatrix"], {left}, {right}, opts]], stateQ],
        qs["Basis"]
    ]

(* Trace *)

QuantumState /: Tr[qs_QuantumState] := Tr @ qs["DensityMatrix"]


(* addition *)

(* qs2 in the basis of qs1 already carries the facts of both *)
addQuantumStates[qs1_QuantumState ? QuantumStateQ, qs2_QuantumState ? QuantumStateQ] /; qs1["Dimension"] == qs2["Dimension"] :=
    With[{q2 = QuantumState[qs2, qs1["Basis"]]},
        QuantumState[
            If[ qs1["StateType"] === qs2["StateType"] === "Vector",
                qs1["StateVector"] + q2["StateVector"],
                qs1["DensityMatrix"] + q2["DensityMatrix"]
            ],
            qs1["Basis"],
            "ParameterSpec" -> MergeParameterSpecs[qs1, qs2],
            "SymbolDomains" -> symbolDomainsOf[q2]
        ]
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
        "ParameterSpec" -> MergeParameterSpecs[qs1, qs2],
        "SymbolDomains" -> MergeSymbolDomains[qs1, qs2]
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
        ArrayVector[blockDiagonalMatrix[#["StateMatrix"] & /@ {qs}]],
        blockDiagonalMatrix[#["DensityMatrix"] & /@ {qs}]
    ],
    Plus @@ (#["Basis"] & /@ {qs}),
    "Label" -> CirclePlus @@ (#["Label"] & /@ {qs}),
    "ParameterSpec" -> MergeParameterSpecs[qs],
    "SymbolDomains" -> MergeSymbolDomains[qs]
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
    q1, q2, q1Mat, q1ConjMat, liftedQ1, liftedQ1Dag, q2Mat, liftedQ2, state, outBasis, inBasis, resBasis
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
                    "ParameterSpec" -> MergeParameterSpecs[qs1, qs2],
        "SymbolDomains" -> MergeSymbolDomains[qs1, qs2]
                ]
            ],
        qs1["VectorQ"] && qs1["Picture"] === "Schrodinger" && qs2["Picture"] === "Schrodinger",
            q1 = qs1["Computational"]; q2 = qs2["Computational"];
            q1Mat = q1["Matrix"];
            liftedQ1 = If[restOutDim == 1, q1Mat, KroneckerProduct[q1Mat, IdentityMatrix[restOutDim, SparseArray]]];
            state = If[TrueQ[q2["VectorQ"]],
                q2Mat = q2["Matrix"];
                liftedQ2 = If[restInDim == 1, q2Mat, KroneckerProduct[q2Mat, IdentityMatrix[restInDim, SparseArray]]];
                ArrayVector[liftedQ1 . liftedQ2],
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
            resBasis = QuantumBasis[
                "Output" -> outBasis, "Input" -> inBasis,
                "Label" -> qs1["Label"] @* qs2["Label"],
                "ParameterSpec" -> MergeParameterSpecs[qs1, qs2],
        "SymbolDomains" -> MergeSymbolDomains[qs1, qs2]
            ];
            If[ resBasis["ComputationalQ"],
                QuantumState[state, resBasis],
                (* the product data is in the computational basis; change it into the basis the result carries *)
                QuantumState[
                    QuantumState[state,
                        "Output" -> QuditBasis[resBasis["OutputDimensions"]],
                        "Input" -> QuditBasis[resBasis["InputDimensions"]]
                    ],
                    resBasis
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
                    "ParameterSpec" -> MergeParameterSpecs[qs1, qs2],
        "SymbolDomains" -> MergeSymbolDomains[qs1, qs2]
                ]
            },
                QuantumState[s, b]
            ]
    ]
]


(* parameterization *)

(qs_QuantumState ? QuantumStateQ)[ps : PatternSequence[p : Except[_Association], ___]] /; ! MemberQ[QuantumState["Properties"], p] && Length[{ps}] <= qs["ParameterArity"] :=
    qs[AssociationThread[Take[qs["Parameters"], UpTo[Length[{ps}]]], {ps}]]

(* Substitution goes through the container rather than over its elements: a
   lazy container binds all its parameters in one whole-array evaluation, where
   an element-wise map would call the underlying function once per amplitude. *)

(* Values that leave an amplitude without a value give a Failure rather than a state
   with such amplitudes: Indeterminate where a closed form divides by a difference of
   eigenvalues that the values make vanish, or where they fall outside every case of
   a Piecewise closed form, as in that of 0^op. *)
(qs_QuantumState ? QuantumStateQ)[rules_ ? AssociationQ] /; ContainsOnly[Keys[rules], qs["Parameters"]] :=
    Replace[valueOutsideFacts[qs, rules], {
        bad_List :> domainFailure[QuantumState, bad],
        None :> Enclose @ With[{state = Confirm[substituteAmplitudes[qs["State"], Normal[rules]]]},
            ConfirmAssert[! newlyValuelessQ[qs["State"], state], "The values leave an amplitude without a value."];
            QuantumState[state, qs["Basis"][rules]]
        ]
    }]

(* A value is checked against the facts its symbol carries before it is put in: a symbol QF takes as real (a
   declared parameter, a gate angle, a symbol declared Reals) refuses a value with an imaginary part, and a
   range, an integer or a positive domain refuses a value outside it. Every step taken before the substitution
   may have used the fact (Conjugate[t] -> t, |e^(I t)|^2 -> 1), so a value that breaks it would give answers
   that are silently wrong. Only numbers are checked: a symbolic value is taken as given. Rounding is allowed:
   an imaginary part, or a step outside a range, below 10^-10 of the value's size passes. The result is the
   first violation, {symbol, value, domain}, or None. *)
$checkValues = True

valueOutsideFacts[_, _, ___] /; ! TrueQ[$checkValues] := None

(* exempt lists symbols whose real requirement is waived (declared Complexes elsewhere in a circuit); their ranges
   are still checked. Values that are not numbers cost nothing: nothing is read from the basis for them. *)
valueOutsideFacts[x_, rules_, exempt_ : {}] := With[{record = basisRecord[x]},
    If[ record === <||> && AllTrue[Values[rules], ! NumericQ[#] || Im[#] == 0 &],
        None,
        valuesAgainst[record, basisDeclared[x], rules, exempt]
    ]
]

(* the check itself, on a record and a list of declared parameters. Real numbers given to symbols with no range,
   the common case, return at once *)
valuesAgainst[record_, declared_, rules_, exempt_] := With[{numeric = Select[Normal[rules], NumericQ @* Last]},
    Which[
        numeric === {}, None,
        record === <||> && AllTrue[numeric, Im[Last[#]] == 0 &], None,
        True, FirstCase[
            numeric,
            (s_ -> v_) :> With[{domains = Join[
                If[realSymbolQ[record, declared, s] && ! MemberQ[exempt, s], {Reals}, {}],
                DeleteCases[Lookup[record, Key[s], {}], Reals | Complexes]
            ]},
                With[{bad = SelectFirst[domains, ! valueInDomainQ[v, #] &, None]}, {s, v, bad} /; bad =!= None]
            ],
            None
        ]
    ]
]

(* real as realDomainSymbols counts it: declared or recorded, and not recorded as Complexes only *)
realSymbolQ[record_, declared_, s_] := With[{d = Lookup[record, Key[s], None]},
    If[ d === None, MemberQ[declared, s], ! MatchQ[d, {Complexes} | {}]]
]

(* the record and the declared parameters, read from the data of the basis of a state, operator, channel or
   measurement operator, without property calls *)
basisOf[qb_QuantumBasis] := qb
basisOf[QuantumState[_, qb_QuantumBasis]] := qb
basisOf[(QuantumOperator | QuantumChannel | QuantumMeasurementOperator)[x : _QuantumState | _QuantumOperator, ___]] := basisOf[x]
basisOf[x_] := x["Basis"]

basisRecord[x_] := symbolDomainsOf[basisOf[x]]
basisDeclared[x_] := declaredOf[basisOf[x]]
declaredOf[QuantumBasis[data_Association]] := Replace[Lookup[data, "ParameterSpec", {}], spec : {__List} :> spec[[All, 1]]]
declaredOf[qb_] := qb["Parameters"]

$domainTolerance = 10 ^ -10

nearlyRealQ[v_] := Im[v] == 0 || Precision[v] =!= Infinity && Abs[Im[v]] <= $domainTolerance Max[1, Abs[v]]

valueInDomainQ[v_, Reals] := nearlyRealQ[v]
valueInDomainQ[v_, Integers] := nearlyRealQ[v] && If[Precision[v] === Infinity, IntegerQ[v], Abs[Re[v] - Round[Re[v]]] <= $domainTolerance Max[1, Abs[v]]]
valueInDomainQ[v_, Interval[{lo_, hi_}]] := nearlyRealQ[v] && With[{x = Re[v], tol = If[Precision[v] === Infinity, 0, $domainTolerance Max[1, Abs[v]]]},
    (lo === - Infinity || x >= lo - tol) && (hi === Infinity || x <= hi + tol)
]
valueInDomainQ[v_, d_Interval] := AnyTrue[List @@ d, valueInDomainQ[v, Interval[#]] &]
valueInDomainQ[v_, NonNegativeReals] := valueInDomainQ[v, Interval[{0, Infinity}]]
valueInDomainQ[v_, PositiveReals] := valueInDomainQ[v, Interval[{0, Infinity}]] && ! (Precision[v] === Infinity && v == 0)
valueInDomainQ[v_, NonPositiveReals] := valueInDomainQ[v, Interval[{- Infinity, 0}]]
valueInDomainQ[v_, NegativeReals] := valueInDomainQ[v, Interval[{- Infinity, 0}]] && ! (Precision[v] === Infinity && v == 0)
valueInDomainQ[v_, NonNegativeIntegers] := valueInDomainQ[v, Integers] && valueInDomainQ[v, Interval[{0, Infinity}]]
valueInDomainQ[v_, PositiveIntegers] := valueInDomainQ[v, Integers] && valueInDomainQ[v, Interval[{1, Infinity}]]
valueInDomainQ[v_, NonPositiveIntegers] := valueInDomainQ[v, Integers] && valueInDomainQ[v, Interval[{- Infinity, 0}]]
valueInDomainQ[v_, NegativeIntegers] := valueInDomainQ[v, Integers] && valueInDomainQ[v, Interval[{- Infinity, -1}]]
valueInDomainQ[v_, Rationals] := nearlyRealQ[v] && (Precision[v] =!= Infinity || Element[v, Rationals] =!= False)
valueInDomainQ[_, _] := True

(* the message is issued under the head of the object the values were given to, as the other QF messages are,
   and the same Failure is returned *)
domainFailure[head_, {s_, v_, d_}] := (
    Message[MessageName[head, "domain"], s, v, domainText[d]];
    Failure["DomainValue", <|"MessageTemplate" :> MessageName[head, "domain"], "MessageParameters" -> {s, v, domainText[d]}|>]
)

domainText[Reals] := "is real here (a declared parameter or a gate angle is taken as real)"
domainText[Integers] := "is an integer here"
domainText[d_Interval] := "lies in " <> ToString[d, InputForm] <> " here"
domainText[d_] := "lies in " <> ToString[d] <> " here"

Scan[
    (MessageName[#, "domain"] = "the value `2` cannot be given to `1`, which `3`; to allow it, give `1` a wider domain with \"SymbolDomains\" when the object is built, Complexes for complex values.") &,
    {QuantumState, QuantumOperator, QuantumChannel, QuantumMeasurementOperator, QuantumCircuitOperator}
]

(* True when an amplitude that had a value before the substitution has none after it,
   read amplitude by amplitude. Before and after it, an amplitude that has not chosen a
   case of a Piecewise has a value, whatever its cases hold, and one without a value
   outside any Piecewise has none. *)
newlyValuelessQ[before_, after_] := valuelessEntriesQ[after] && If[
    ArrayQ[before] && Dimensions[before] === Dimensions[after],
    Or @@ MapThread[! valuelessOutsidePiecewiseQ[#1] && valuelessOutsidePiecewiseQ[#2] &, Flatten /@ Normal /@ {before, after}],
    AnyTrue[Flatten[{Normal[after]}], valuelessOutsidePiecewiseQ]
]

valuelessOutsidePiecewiseQ[x_] := ! FreeQ[x /. HoldPattern[Piecewise[___]] :> 0, Indeterminate | Undefined | _DirectedInfinity]

(* The amplitudes of f[qo] for an operator with declared parameters, a Function of
   the parameters whose body is matrixMapAmplitudes, are substituted through
   the unevaluated body: binding every parameter calls it, binding some writes the
   values into the held body and keeps the rest as its parameters. ArrayReplaceAll
   curries through the body evaluated at the parameters instead, which is the
   closed form the Function exists to avoid substituting into. *)
substituteAmplitudes[state_ ? lazyMatrixMapAmplitudesQ, rules_List] := With[{
    parameters = Flatten[{First[state]}]
},
    With[{free = DeleteCases[parameters, Alternatives @@ Keys[rules]]},
        If[ free === {},
            state @@ Replace[parameters, rules, {1}],
            With[{curried = Function @@ Prepend[Extract[state, {2}, Hold] /. rules, free]},
                ArrayDeclareShape[curried, ArrayDimensions[state]];
                curried
            ]
        ]
    ]
]

substituteAmplitudes[state_, rules_] := ArrayReplaceAll[state, rules]

