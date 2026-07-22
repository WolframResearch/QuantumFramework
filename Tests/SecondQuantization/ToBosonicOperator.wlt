BeginTestSection["ToBosonicOperator - single mode"]

(* a^dagger a (number operator) via ** matches direct composition *)
VerificationTest[
    With[{a = AnnihilationOperator[6]},
        ToBosonicOperator[SuperDagger[\[FormalA]] ** \[FormalA], 6]["Matrix"] == (SuperDagger[a] @ a)["Matrix"]
    ],
    True,
    TestID -> "ToBosonicOperator-NumberOperator"
]

(* Only the dagger appears; automatic var detection must still find the mode *)
VerificationTest[
    ToBosonicOperator[SuperDagger[\[FormalA]], 6]["Matrix"] == AnnihilationOperator[6]["Dagger"]["Matrix"],
    True,
    TestID -> "ToBosonicOperator-OnlyDagger"
]

(* (a^dagger)^2 a^2 via repeated ** matches direct repeated composition *)
VerificationTest[
    With[{a = AnnihilationOperator[]},
        ToBosonicOperator[SuperDagger[\[FormalA]] ** SuperDagger[\[FormalA]] ** \[FormalA] ** \[FormalA], $FockSize]["Matrix"] ==
            (SuperDagger[a] @ SuperDagger[a] @ a @ a)["Matrix"]
    ],
    True,
    TestID -> "ToBosonicOperator-RepeatedComposition"
]

(* Explicit size overrides $FockSize *)
VerificationTest[
    ToBosonicOperator[SuperDagger[\[FormalA]] ** \[FormalA], {\[FormalA]}, 8]["Matrix"] ==
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
            ToBosonicOperator[Exp[concreteGen], 6]["Matrix"] == Exp[concreteGen]["Matrix"]
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
        {v1, v1dag, v2, v2dag} = FieldVariables[{1, 2}];
        ToBosonicOperator[v1dag ** v2 + v2dag ** v1, {v1, v2}, 6]["Matrix"] ==
            (SuperDagger[a1] @ a2 + SuperDagger[a2] @ a1)["Matrix"]
    ],
    True,
    TestID -> "ToBosonicOperator-TwoModeCoupling"
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
            ToBosonicOperator[Exp[alpha SuperDagger[\[FormalA]] - Conjugate[alpha] \[FormalA]], 6]["Matrix"] -
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
            ToBosonicOperator[Exp[alpha SuperDagger[\[FormalA]] - Conjugate[alpha] \[FormalA]], 6]["Matrix"] -
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
            ToBosonicOperator[Exp[alpha vdag - Conjugate[alpha] v], {v}, 6]["Matrix"] -
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
            ToBosonicOperator[Cos[SuperDagger[\[FormalA]] + \[FormalA]], 6]["Matrix"] -
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
            ToBosonicOperator[2^(SuperDagger[\[FormalA]] + \[FormalA]), 6]["Matrix"] -
                MatrixExp[Log[2] (SuperDagger[a] + a)]["Matrix"]
        ]]]
    ],
    0,
    TestID -> "ToBosonicOperator-GeneralScalarBasePower"
]

EndTestSection[]
