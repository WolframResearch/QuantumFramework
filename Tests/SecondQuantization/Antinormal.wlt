
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

(* Shared scaffolding: a battery of monomials/polynomials to check invariants over, and a
   structural predicate asserting that no creator precedes an annihilator in any monomial *)
$antiVars = {\[FormalA], SuperDagger[\[FormalA]]};

$antiBattery = {
    SuperDagger[\[FormalA]] ** \[FormalA],
    SuperDagger[\[FormalA]] ** SuperDagger[\[FormalA]] ** \[FormalA] ** \[FormalA],
    SuperDagger[\[FormalA]] ** \[FormalA] ** SuperDagger[\[FormalA]],
    GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 3] ** \[FormalA],
    \[FormalA] ** SuperDagger[\[FormalA]] ** \[FormalA] ** SuperDagger[\[FormalA]],
    SuperDagger[\[FormalA]] ** \[FormalA] + 3 SuperDagger[\[FormalA]] ** SuperDagger[\[FormalA]] ** \[FormalA],
    2 GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] ** \[FormalA] - 5 SuperDagger[\[FormalA]],
    -2 \[FormalA] ** SuperDagger[\[FormalA]]
};

$antiOrderedQ[expr_, vars_] := AllTrue[
    If[Head[expr] === Plus, List @@ expr, {expr}],
    Function[term,
        OrderedQ @ Map[
            Boole[! FreeQ[#, SuperDagger]] &,
            Select[
                Replace[
                    FirstCase[term, _NonCommutativeMultiply, term, {0, Infinity}],
                    {NonCommutativeMultiply[f__] :> {f}, f_ :> {f}}
                ],
                ! FreeQ[#, Alternatives @@ vars] &
            ]
        ]
    ]
];


BeginTestSection["BosonicAntiNormalOrder - single mode"]

(* a^dagger a = a a^dagger - 1 *)
VerificationTest[
    BosonicAntiNormalOrder[SuperDagger[\[FormalA]] ** \[FormalA], $antiVars],
    -1 + \[FormalA] ** SuperDagger[\[FormalA]],
    TestID -> "AntiNormal-NumberOperator"
]

(* (a^dagger)^2 a^2 = a^2 (a^dagger)^2 - 4 a a^dagger + 2 *)
VerificationTest[
    BosonicAntiNormalOrder[
        SuperDagger[\[FormalA]] ** SuperDagger[\[FormalA]] ** \[FormalA] ** \[FormalA], $antiVars],
    2 - 4 \[FormalA] ** SuperDagger[\[FormalA]] +
        GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2] **
            GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2],
    TestID -> "AntiNormal-SquaredNumberOperator"
]

(* a^dagger a a^dagger = a (a^dagger)^2 - a^dagger *)
VerificationTest[
    BosonicAntiNormalOrder[
        SuperDagger[\[FormalA]] ** \[FormalA] ** SuperDagger[\[FormalA]], $antiVars],
    \[FormalA] ** GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] -
        SuperDagger[\[FormalA]],
    TestID -> "AntiNormal-Sandwich"
]

(* An already anti-ordered monomial is a fixed point *)
VerificationTest[
    BosonicAntiNormalOrder[\[FormalA] ** SuperDagger[\[FormalA]], $antiVars],
    \[FormalA] ** SuperDagger[\[FormalA]],
    TestID -> "AntiNormal-Idempotent"
]

(* Structural invariant: no creator ever precedes an annihilator in the result *)
VerificationTest[
    AllTrue[$antiBattery, $antiOrderedQ[BosonicAntiNormalOrder[#, $antiVars], $antiVars] &],
    True,
    TestID -> "AntiNormal-StructurallyAntiOrdered"
]

(* Algebraic invariant: re-normal-ordering the anti-ordered form recovers the normal form,
   i.e. the reduction rewrites the operator without changing it *)
VerificationTest[
    AllTrue[$antiBattery,
        Simplify[
            BosonicNormalOrder[BosonicAntiNormalOrder[#, $antiVars], $antiVars] -
                BosonicNormalOrder[#, $antiVars]
        ] === 0 &],
    True,
    TestID -> "AntiNormal-RoundTripThroughNormalOrder"
]

EndTestSection[]


BeginTestSection["BosonicAntiNormalOrder - two modes"]

(* Modes commute with each other, only each mode's own pair reorders:
   a^dagger a b^dagger b = a b a^dagger b^dagger - a a^dagger - b b^dagger + 1 *)
VerificationTest[
    BosonicAntiNormalOrder[
        SuperDagger[\[FormalA]] ** \[FormalA] ** SuperDagger[\[FormalB]] ** \[FormalB],
        {\[FormalA], SuperDagger[\[FormalA]], \[FormalB], SuperDagger[\[FormalB]]}],
    1 - \[FormalA] ** SuperDagger[\[FormalA]] - \[FormalB] ** SuperDagger[\[FormalB]] +
        \[FormalA] ** \[FormalB] ** SuperDagger[\[FormalA]] ** SuperDagger[\[FormalB]],
    TestID -> "AntiNormal-TwoModeNumberProduct"
]

(* Mixed two-mode polynomial: check the operator is preserved *)
VerificationTest[
    Module[{vars = {\[FormalA], SuperDagger[\[FormalA]], \[FormalB], SuperDagger[\[FormalB]]}, x},
        x = 2 SuperDagger[\[FormalA]] ** \[FormalA] +
            SuperDagger[\[FormalA]] ** SuperDagger[\[FormalB]] ** \[FormalA] ** \[FormalB];
        Simplify[
            BosonicNormalOrder[BosonicAntiNormalOrder[x, vars], vars] - BosonicNormalOrder[x, vars]
        ] === 0
    ],
    True,
    TestID -> "AntiNormal-TwoModeRoundTrip"
]

EndTestSection[]


BeginTestSection["BosonicAntiNormalOrder - scalars"]

(* Symbolic coefficients declared through "Scalars" must survive the reduction *)
VerificationTest[
    BosonicAntiNormalOrder[
        \[Alpha] SuperDagger[\[FormalA]] ** \[FormalA] +
            \[Beta] SuperDagger[\[FormalA]] ** SuperDagger[\[FormalA]] ** \[FormalA] ** \[FormalA],
        $antiVars,
        "Scalars" -> {\[Alpha], \[Beta]}],
    -\[Alpha] + 2 \[Beta] + \[Alpha] \[FormalA] ** SuperDagger[\[FormalA]] -
        4 \[Beta] \[FormalA] ** SuperDagger[\[FormalA]] +
        \[Beta] GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2] **
            GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2],
    TestID -> "AntiNormal-ScalarCoefficients"
]

EndTestSection[]


BeginTestSection["BosonicAntiNormalOrder - methods"]

(* "Blasiak" reaches anti-normal order by conjugating the normal-ordering formula with the
   automorphism a -> a^dagger, a^dagger -> -a; it must agree with the Grobner reduction *)
VerificationTest[
    AllTrue[$antiBattery,
        Simplify[
            BosonicAntiNormalOrder[#, $antiVars, Method -> "Blasiak"] -
                BosonicAntiNormalOrder[#, $antiVars]
        ] === 0 &],
    True,
    TestID -> "AntiNormal-BlasiakMatchesGrobner"
]

(* The conjugating automorphism introduces a sign on odd-degree creator blocks, so the Blasiak
   path must accept a monomial carrying a scalar coefficient *)
VerificationTest[
    Simplify[
        BosonicAntiNormalOrder[
            GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] ** \[FormalA],
            $antiVars, Method -> "Blasiak"] -
        (\[FormalA] ** GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] -
            2 SuperDagger[\[FormalA]])
    ] === 0,
    True,
    TestID -> "AntiNormal-BlasiakSignedMonomial"
]

(* Same requirement seen directly on BosonicNormalOrder: a coefficient in front of a monomial
   is not part of the ladder word and must be carried through the Blasiak reduction *)
VerificationTest[
    Simplify[
        BosonicNormalOrder[-2 \[FormalA] ** SuperDagger[\[FormalA]], $antiVars, Method -> "Blasiak"] -
            BosonicNormalOrder[-2 \[FormalA] ** SuperDagger[\[FormalA]], $antiVars]
    ] === 0,
    True,
    TestID -> "AntiNormal-BlasiakCoefficientCarried"
]


VerificationTest[
    BosonicAntiNormalOrder[SuperDagger[\[FormalA]] ** \[FormalA], $antiVars, Method -> "NoSuchMethod"],
    $Failed,
    {BosonicAntiNormalOrder::unknownMethod},
    TestID -> "AntiNormal-UnknownMethod"
]

EndTestSection[]


BeginTestSection["BosonicAntiNormalOrder - matrix realization"]

(* Independent numeric check against truncated Fock matrices: the anti-ordered polynomial must
   act identically to the original.  Terms are realized separately so that the c-number part
   becomes a multiple of the identity, and only a sub-block is compared because a truncated
   ladder operator differs from the exact one near the top of the space. *)
VerificationTest[
    Module[{size = 12, block = 8, am, adm, toMatrix, x, anti},
        am = Normal[AnnihilationOperator[size]["Matrix"]];
        adm = ConjugateTranspose[am];
        toMatrix[e_] := Total @ Map[
            If[ FreeQ[#, \[FormalA]],
                # IdentityMatrix[size],
                # /. {
                    GeneralizedPower[NonCommutativeMultiply, op_, p_] :>
                        MatrixPower[op /. {SuperDagger[\[FormalA]] -> adm, \[FormalA] -> am}, p],
                    SuperDagger[\[FormalA]] -> adm,
                    \[FormalA] -> am,
                    NonCommutativeMultiply -> Dot
                }
            ] &,
            If[Head[e] === Plus, List @@ e, {e}]
        ];
        x = SuperDagger[\[FormalA]] ** SuperDagger[\[FormalA]] ** \[FormalA] ** \[FormalA];
        anti = BosonicAntiNormalOrder[x, $antiVars];
        toMatrix[x][[;; block, ;; block]] == toMatrix[anti][[;; block, ;; block]]
    ],
    True,
    TestID -> "AntiNormal-MatrixRealization"
]

EndTestSection[]
