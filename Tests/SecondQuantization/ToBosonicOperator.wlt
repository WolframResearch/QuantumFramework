Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

$v = FieldVariables[];
$vab = FieldVariables[{\[FormalA], \[FormalB]}];


BeginTestSection["ToBosonicOperator - single mode"]

(* a^dagger a (number operator) via ** matches direct composition *)
VerificationTest[
    With[{a = AnnihilationOperator[6]},
        ToBosonicOperator[SuperDagger[\[FormalA]] ** \[FormalA], $v, 6]["Matrix"] == (SuperDagger[a] @ a)["Matrix"]
    ],
    True,
    TestID -> "ToBosonicOperator-NumberOperator"
]

(* Only the creator appears; it realizes as the adjoint of its mode *)
VerificationTest[
    ToBosonicOperator[SuperDagger[\[FormalA]], $v, 6]["Matrix"] == AnnihilationOperator[6]["Dagger"]["Matrix"],
    True,
    TestID -> "ToBosonicOperator-OnlyDagger"
]

(* (a^dagger)^2 a^2 via repeated ** matches direct repeated composition *)
VerificationTest[
    With[{a = AnnihilationOperator[]},
        ToBosonicOperator[SuperDagger[\[FormalA]] ** SuperDagger[\[FormalA]] ** \[FormalA] ** \[FormalA], $v, $FockSize]["Matrix"] ==
            (SuperDagger[a] @ SuperDagger[a] @ a @ a)["Matrix"]
    ],
    True,
    TestID -> "ToBosonicOperator-RepeatedComposition"
]

(* Explicit size overrides $FockSize *)
VerificationTest[
    ToBosonicOperator[SuperDagger[\[FormalA]] ** \[FormalA], $v, 8]["Matrix"] ==
        (SuperDagger[AnnihilationOperator[8]] @ AnnihilationOperator[8])["Matrix"],
    True,
    TestID -> "ToBosonicOperator-ExplicitSize"
]

(* A normal-ordered chain (as returned by BosonicNormalOrder) using GeneralizedPower notation *)
VerificationTest[
    With[{a = AnnihilationOperator[6]},
        ToBosonicOperator[
            GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] **
                GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2],
            $v,
            6
        ]["Matrix"] == (SuperDagger[a] @ SuperDagger[a] @ a @ a)["Matrix"]
    ],
    True,
    TestID -> "ToBosonicOperator-GeneralizedPowerChain"
]

(* Passing an already-realized QuantumOperator through (no formal symbols present) is a no-op *)
VerificationTest[
    With[{alpha = 0.3 + 0.1 I, a = AnnihilationOperator[6]},
        With[{concreteGen = alpha SuperDagger[a] - Conjugate[alpha] a},
            ToBosonicOperator[Exp[concreteGen], $v, 6]["Matrix"] == Exp[concreteGen]["Matrix"]
        ]
    ],
    True,
    TestID -> "ToBosonicOperator-PassThroughConcrete"
]

EndTestSection[]


BeginTestSection["ToBosonicOperator - two mode"]

(* a1^dagger a2 + a2^dagger a1 (beamsplitter-like coupling), built from field variables pre-bound to
   local names via Module, as FieldVariables[] is typically used throughout this codebase *)
VerificationTest[
    Module[{a1 = AnnihilationOperator[6, {1}], a2 = AnnihilationOperator[6, {2}], v1, v1dag, v2, v2dag},
        {{v1dag, v2dag}, {v1, v2}} = FieldVariables[{1, 2}];
        ToBosonicOperator[v1dag ** v2 + v2dag ** v1, FieldVariables[{1, 2}], 6]["Matrix"] ==
            (SuperDagger[a1] @ a2 + SuperDagger[a2] @ a1)["Matrix"]
    ],
    True,
    TestID -> "ToBosonicOperator-TwoModeCoupling"
]

VerificationTest[
    Module[{a1 = AnnihilationOperator[4, {1}], a2 = AnnihilationOperator[4, {2}]},
        ToBosonicOperator[\[FormalA] ** \[FormalB] + SuperDagger[\[FormalA]], $vab, 4]["Matrix"] ==
            (a1 @ a2 + SuperDagger[a1])["Matrix"]
    ],
    True,
    TestID -> "ToBosonicOperator-ProductPlusBareOperator"
]

VerificationTest[
    Module[{a1 = AnnihilationOperator[4, {1}], a2 = AnnihilationOperator[4, {2}]},
        Chop[Norm[Flatten[
            ToBosonicOperator[Exp[0.1 (\[FormalA] ** \[FormalB] + SuperDagger[\[FormalA]])], $vab, 4]["Matrix"] -
                Exp[0.1 (a1 @ a2 + SuperDagger[a1])]["Matrix"]
        ]]]
    ],
    0,
    TestID -> "ToBosonicOperator-ProductPlusBareOperatorInExp"
]

VerificationTest[
    With[{a = AnnihilationOperator[6]},
        ToBosonicOperator[
            Evaluate @ BosonicNormalOrder[\[FormalA] ** SuperDagger[\[FormalA]] ** \[FormalA], $v],
            $v,
            6
        ]["Matrix"] == (a + SuperDagger[a] @ a @ a)["Matrix"]
    ],
    True,
    TestID -> "ToBosonicOperator-NormalOrderedPolynomial"
]

EndTestSection[]


BeginTestSection["ToBosonicOperator - exponentials and other wrapping functions"]

(* Weak-order displacement exponential should match MatrixExp and DisplacementOperator directly.
   Regression test: Exp[symbolicExpr] eagerly canonicalizes to Power[E, symbolicExpr] before field vars
   are substituted for real operators, which without care dispatches as a matrix power of E instead of
   the matrix exponential. *)
VerificationTest[
    Module[{alpha = 0.3 + 0.1 I, a = AnnihilationOperator[6]},
        Chop[Norm[Flatten[
            ToBosonicOperator[Exp[alpha SuperDagger[\[FormalA]] - Conjugate[alpha] \[FormalA]], $v, 6]["Matrix"] -
                MatrixExp[alpha SuperDagger[a] - Conjugate[alpha] a]["Matrix"]
        ]]]
    ],
    0,
    TestID -> "ToBosonicOperator-WeakDisplacementExponential"
]

(* Independent cross-check against the dedicated DisplacementOperator constructor *)
VerificationTest[
    Module[{alpha = 0.3 + 0.1 I},
        Chop[Norm[Flatten[
            ToBosonicOperator[Exp[alpha SuperDagger[\[FormalA]] - Conjugate[alpha] \[FormalA]], $v, 6]["Matrix"] -
                DisplacementOperator[alpha, 6, "Ordering" -> "Weak"]["Matrix"]
        ]]]
    ],
    0,
    TestID -> "ToBosonicOperator-MatchesDisplacementOperator"
]

(* Exp built up from local variables (v = field var, vdag = SuperDagger[v] computed beforehand), the
   natural way FieldVariables[] results get used throughout this codebase *)
VerificationTest[
    Module[{alpha = 0.3 + 0.1 I, a = AnnihilationOperator[6], v = \[FormalA], vdag},
        vdag = SuperDagger[v];
        Chop[Norm[Flatten[
            ToBosonicOperator[Exp[alpha vdag - Conjugate[alpha] v], {{vdag}, {v}}, 6]["Matrix"] -
                Exp[alpha SuperDagger[a] - Conjugate[alpha] a]["Matrix"]
        ]]]
    ],
    0,
    TestID -> "ToBosonicOperator-ExpWithIntermediateVariables"
]

(* Cos of a Hermitian combination, checking any wrapping function (not just Exp) realizes correctly *)
VerificationTest[
    With[{a = AnnihilationOperator[6]},
        Chop[Norm[Flatten[
            ToBosonicOperator[Cos[SuperDagger[\[FormalA]] + \[FormalA]], $v, 6]["Matrix"] -
                Cos[SuperDagger[a] + a]["Matrix"]
        ]]]
    ],
    0,
    TestID -> "ToBosonicOperator-CosOfOperator"
]

(* A plain scalar base raised to an operator power, not going through Exp at all *)
VerificationTest[
    With[{a = AnnihilationOperator[6]},
        Chop[Norm[Flatten[
            ToBosonicOperator[2^(SuperDagger[\[FormalA]] + \[FormalA]), $v, 6]["Matrix"] -
                MatrixExp[Log[2] (SuperDagger[a] + a)]["Matrix"]
        ]]]
    ],
    0,
    TestID -> "ToBosonicOperator-GeneralScalarBasePower"
]

(* A squared operator in an exponent, as the squeeze orderings return it: the GeneralizedPower inside
   the exponent is realized before the exponential, or MatrixExp is left holding it. *)
VerificationTest[
    With[{a = Normal[AnnihilationOperator[5]["Matrix"]]},
        Simplify[Normal[ToBosonicOperator[Exp[z GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]], $v, 5]["Matrix"]] -
            MatrixExp[z Transpose[a] . Transpose[a]]]
    ],
    ConstantArray[0, {5, 5}],
    TestID -> "ToBosonicOperator-SquareInSymbolicExponent"
]

(* Ordered outputs can carry a commutative product of operator powers inside **, here c^N d^N; with
   symbolic c and d QF does not multiply the two operators, so the product is composed instead. *)
VerificationTest[
    With[{a = Normal[AnnihilationOperator[4]["Matrix"]], n = SuperDagger[\[FormalA]] ** \[FormalA]},
        Simplify[Normal[ToBosonicOperator[Exp[SuperDagger[\[FormalA]]] ** (c^n d^n), $v, 4]["Matrix"]] -
            MatrixExp[Transpose[a]] . DiagonalMatrix[(c d)^Range[0, 3]]]
    ],
    ConstantArray[0, {4, 4}],
    TestID -> "ToBosonicOperator-CommutingProductOfOperatorPowers"
]

EndTestSection[]


BeginTestSection["ToBosonicOperator - ordered products"]

$exactBlock[gen_, big_, modes_] :=
    With[{a = N @ Normal[AnnihilationOperator[big]["Matrix"]]},
        With[{ops = Table[Fold[KroneckerProduct, ReplacePart[ConstantArray[IdentityMatrix[big], modes], j -> a]], {j, modes}]},
            With[{idx = If[modes === 1, Range[4], Flatten[Table[big i + j + 1, {i, 0, 2}, {j, 0, 2}]]]},
                MatrixExp[gen[ops, ConjugateTranspose /@ ops]][[idx, idx]]
            ]
        ]
    ];

$block[op_, size_, modes_] :=
    With[{idx = If[modes === 1, Range[4], Flatten[Table[size i + j + 1, {i, 0, 2}, {j, 0, 2}]]]},
        N[Normal[op["Matrix"]]][[idx, idx]]
    ];

VerificationTest[
    Max[Abs[N[Normal[ToBosonicOperator[BosonicExpOrder[Exp[1/2 SuperDagger[\[FormalA]] ** \[FormalA]], $v], $v, 6]["Matrix"]]] -
        DiagonalMatrix[N[Exp[Range[0, 5]/2]]]]] < 10^-12,
    True,
    TestID -> "ToBosonicOperator-NormalOrderedNumberExponential"
]

VerificationTest[
    Normal[ToBosonicOperator[NormalOrdered[Exp[-SuperDagger[\[FormalA]] ** \[FormalA]], $v], $v, 5]["Matrix"]],
    Normal[SparseArray[{1, 1} -> 1, {5, 5}]],
    TestID -> "ToBosonicOperator-NormalOrderedVacuumProjector"
]

VerificationTest[
    Max[Abs[$block[ToBosonicOperator[BosonicExpOrder[Exp[0.1 SuperDagger[\[FormalA]] ** \[FormalA] ** SuperDagger[\[FormalA]]], $v], $v, 6], 6, 1] -
        $exactBlock[0.1 #2[[1]] . #1[[1]] . #2[[1]] &, 40, 1]]] < 10^-12,
    True,
    TestID -> "ToBosonicOperator-NormalOrderedBlasiak"
]

VerificationTest[
    Max[Abs[$block[ToBosonicOperator[BosonicExpOrder[Exp[0.3 (SuperDagger[\[FormalA]] ** \[FormalB] - SuperDagger[\[FormalB]] ** \[FormalA])], $vab], $vab, 6], 6, 2] -
        $exactBlock[0.3 (#2[[1]] . #1[[2]] - #2[[2]] . #1[[1]]) &, 14, 2]]] < 10^-12,
    True,
    TestID -> "ToBosonicOperator-NormalOrderedBeamSplitter"
]

VerificationTest[
    Max[Abs[$block[ToBosonicOperator[BosonicExpOrder[
            Exp[0.2 SuperDagger[\[FormalA]] ** \[FormalB] + 0.1 SuperDagger[\[FormalA]] ** \[FormalA] + 0.05 SuperDagger[\[FormalB]] ** \[FormalA]], $vab], $vab, 6], 6, 2] -
        $exactBlock[0.2 #2[[1]] . #1[[2]] + 0.1 #2[[1]] . #1[[1]] + 0.05 #2[[2]] . #1[[1]] &, 14, 2]]] < 10^-12,
    True,
    TestID -> "ToBosonicOperator-NormalOrderedPassiveMixing"
]

VerificationTest[
    Max[Abs[$block[ToBosonicOperator[BosonicExpOrder[
            Exp[0.1 (GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] - GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]) +
                0.15 (SuperDagger[\[FormalA]] ** \[FormalB] - SuperDagger[\[FormalB]] ** \[FormalA])], $vab], $vab, 6], 6, 2] -
        $exactBlock[0.1 (#2[[1]] . #2[[1]] - #1[[1]] . #1[[1]]) + 0.15 (#2[[1]] . #1[[2]] - #2[[2]] . #1[[1]]) &, 14, 2]]] < 10^-10,
    True,
    TestID -> "ToBosonicOperator-NormalOrderedSqueezeAndMixing"
]

VerificationTest[
    Max[Abs[$block[ToBosonicOperator[BosonicExpOrder[Exp[0.3 (SuperDagger[\[FormalA]] ** \[FormalB] - SuperDagger[\[FormalB]] ** \[FormalA])], $vab,
            "Ordering" -> "Antinormal"], $vab, 12], 12, 2] -
        $exactBlock[0.3 (#2[[1]] . #1[[2]] - #2[[2]] . #1[[1]]) &, 14, 2]]] < 10^-8,
    True,
    TestID -> "ToBosonicOperator-AntinormalOrderedBeamSplitter"
]

EndTestSection[]
