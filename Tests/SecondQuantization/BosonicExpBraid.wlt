Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

$a = \[FormalA];
$ad = SuperDagger[\[FormalA]];
$n = $ad ** $a;
$disp = \[Alpha] $ad - Conjugate[\[Alpha]] $a;

exponentOf[Exp[x_]] := x;
exponentOf[_] := $Failed;
sameGeneratorQ[res_, expected_] := FullSimplify[exponentOf[res] - expected] === 0;

(* Factor-wise comparison: two products of exponentials agree when their
   generators do. Structural === is too strict, since FullSimplify may leave
   one side factored (\[Alpha](x+y)) and the other expanded (\[Alpha]x+\[Alpha]y). *)
sameBraidQ[res_, expected_] :=
    Length[res] === Length[expected] &&
    AllTrue[
        Transpose[{List @@ res, List @@ expected}],
        FullSimplify[exponentOf[#[[1]]] - exponentOf[#[[2]]]] === 0 &
    ];

(* Truncated-oscillator realization via the paclet's own realizer, which maps a
   c-number summand to a multiple of the identity - matrix - scalar would
   instead subtract the scalar from every entry. N first, so MatrixExp never
   sees an exact matrix. *)
matrixOf[expr_, num_, size_ : 60] :=
    With[{e = N[expr /. num]},
        ToBosonicOperator[e, {$a}, size]["MatrixRepresentation"]
    ];
braidHoldsQ[input_, out_, num_, keep_ : 12] :=
    With[{blk = Normal[#][[;; keep, ;; keep]] &},
        Max[Abs[blk[matrixOf[input, num]] - blk[matrixOf[out, num]]]] < 10.^-10
    ];


BeginTestSection["BosonicExpBraid - number-operator generator"]

(* The identity that motivated the function:
   e^(\[Beta] n) D(\[Alpha]) = e^(\[Alpha] e^\[Beta] a\[Dagger] - \[Alpha]\[Conjugate] e^-\[Beta] a) e^(\[Beta] n). *)
VerificationTest[
    BosonicExpBraid[Exp[\[Beta] $n] ** Exp[$disp]],
    Exp[\[Alpha] Exp[\[Beta]] $ad - Conjugate[\[Alpha]] Exp[-\[Beta]] $a] ** Exp[\[Beta] $n],
    {},
    TestID -> "BosonicExpBraid-NumberPastDisplacement"
]

VerificationTest[
    braidHoldsQ[
        Exp[\[Beta] $n] ** Exp[$disp],
        BosonicExpBraid[Exp[\[Beta] $n] ** Exp[$disp]],
        {\[Beta] -> 0.3 + 0.2 I, \[Alpha] -> 0.5 - 0.4 I}
    ],
    True,
    {},
    TestID -> "BosonicExpBraid-NumberPastDisplacementNumeric"
]

(* Any partner is allowed, including GeneralizedPower and unordered products.
   Each a\[Dagger] picks up e^\[Lambda], each a picks up e^-\[Lambda]. *)
VerificationTest[
    BosonicExpBraid[Exp[\[Lambda] $n] ** Exp[\[Xi] $ad ** $ad]],
    Exp[\[Xi] Exp[2 \[Lambda]] $ad ** $ad] ** Exp[\[Lambda] $n],
    {},
    TestID -> "BosonicExpBraid-NumberPastQuadraticPartner"
]

VerificationTest[
    BosonicExpBraid[
        Exp[\[Lambda] $n] ** Exp[\[Xi] GeneralizedPower[NonCommutativeMultiply, $a, 3]]],
    Exp[\[Xi] Exp[-3 \[Lambda]] GeneralizedPower[NonCommutativeMultiply, $a, 3]] ** Exp[\[Lambda] $n],
    {},
    TestID -> "BosonicExpBraid-NumberPastGeneralizedPower"
]

(* a**a\[Dagger] = n + 1 differs from n by a central term, so the ad-action is the
   same and this partner comes back untouched - it is NOT normal-ordered. *)
VerificationTest[
    BosonicExpBraid[Exp[\[Lambda] $n] ** Exp[\[Xi] $a ** $ad]],
    Exp[\[Xi] $a ** $ad] ** Exp[\[Lambda] $n],
    {},
    TestID -> "BosonicExpBraid-DoesNotOrderThePartner"
]

EndTestSection[]


BeginTestSection["BosonicExpBraid - linear and quadratic generators"]

(* Ad_(\[Alpha] a\[Dagger]) shifts a by -\[Alpha] and fixes a\[Dagger]. *)
VerificationTest[
    BosonicExpBraid[Exp[\[Alpha] $ad] ** Exp[\[Lambda] $n]],
    Exp[\[Lambda] $ad ** ($a - \[Alpha])] ** Exp[\[Alpha] $ad],
    {},
    TestID -> "BosonicExpBraid-CreationPastNumber"
]

VerificationTest[
    braidHoldsQ[
        Exp[\[Alpha] $ad] ** Exp[\[Lambda] $n],
        BosonicExpBraid[Exp[\[Alpha] $ad] ** Exp[\[Lambda] $n]],
        {\[Lambda] -> -0.35 + 0.1 I, \[Alpha] -> 0.4 - 0.2 I}
    ],
    True,
    {},
    TestID -> "BosonicExpBraid-CreationPastNumberNumeric"
]

(* Squeeze generator: Ad is the Bogoliubov rotation
   a\[Dagger] -> Cosh[r] a\[Dagger] - Sinh[r] a at r = 1/2. *)
VerificationTest[
    sameBraidQ[
        BosonicExpBraid[
            Exp[(1/4) (GeneralizedPower[NonCommutativeMultiply, $ad, 2] -
                GeneralizedPower[NonCommutativeMultiply, $a, 2])] ** Exp[\[Alpha] $ad]],
        Exp[\[Alpha] (Cosh[1/2] $ad - Sinh[1/2] $a)] **
            Exp[(1/4) (GeneralizedPower[NonCommutativeMultiply, $ad, 2] -
                GeneralizedPower[NonCommutativeMultiply, $a, 2])]
    ],
    True,
    {},
    TestID -> "BosonicExpBraid-SqueezePastCreation"
]

VerificationTest[
    braidHoldsQ[
        Exp[(1/4) ($ad ** $ad - $a ** $a)] ** Exp[\[Alpha] $ad],
        BosonicExpBraid[Exp[(1/4) ($ad ** $ad - $a ** $a)] ** Exp[\[Alpha] $ad]],
        {\[Alpha] -> 0.3 - 0.25 I}
    ],
    True,
    {},
    TestID -> "BosonicExpBraid-SqueezePastCreationNumeric"
]

EndTestSection[]


BeginTestSection["BosonicExpBraid - direction"]

(* Left transport: e^A e^B = e^B e^(Ad_-B A). *)
VerificationTest[
    BosonicExpBraid[Exp[\[Beta] $n] ** Exp[\[Alpha] $ad], "Direction" -> "Left"],
    Exp[\[Alpha] $ad] ** Exp[\[Beta] $ad ** ($a + \[Alpha])],
    {},
    TestID -> "BosonicExpBraid-LeftDirection"
]

VerificationTest[
    braidHoldsQ[
        Exp[\[Beta] $n] ** Exp[\[Alpha] $ad],
        BosonicExpBraid[Exp[\[Beta] $n] ** Exp[\[Alpha] $ad], "Direction" -> "Left"],
        {\[Beta] -> 0.2 - 0.3 I, \[Alpha] -> 0.45 + 0.15 I}
    ],
    True,
    {},
    TestID -> "BosonicExpBraid-LeftDirectionNumeric"
]

VerificationTest[
    BosonicExpBraid[Exp[\[Beta] $n] ** Exp[$disp], "Direction" -> "Sideways"],
    $Failed,
    {BosonicExpBraid::direction},
    TestID -> "BosonicExpBraid-RejectsUnknownDirection"
]

EndTestSection[]


BeginTestSection["BosonicExpBraid - stays in its lane"]

(* Commuting exponentials are swapped, never merged - merging belongs to
   BosonicExpSimplify. *)
VerificationTest[
    BosonicExpBraid[Exp[\[Beta] $n] ** Exp[\[Gamma] $n]],
    Exp[\[Gamma] $n] ** Exp[\[Beta] $n],
    {},
    TestID -> "BosonicExpBraid-DoesNotMergeCommutingPair"
]

(* An exponential and its inverse are swapped, never cancelled. *)
VerificationTest[
    BosonicExpBraid[Exp[\[Beta] $n] ** Exp[-\[Beta] $n]],
    Exp[-\[Beta] $n] ** Exp[\[Beta] $n],
    {},
    TestID -> "BosonicExpBraid-DoesNotCancelInversePair"
]

(* A generator outside the recognized families transports nothing, silently. *)
VerificationTest[
    BosonicExpBraid[Exp[\[Lambda] $ad ** $a ** $a] ** Exp[\[Alpha] $ad]],
    Exp[\[Lambda] $ad ** $a ** $a] ** Exp[\[Alpha] $ad],
    {},
    TestID -> "BosonicExpBraid-UnrecognizedGeneratorIsInert"
]

(* Mixed families in one generator (number + linear) are not claimed. *)
VerificationTest[
    BosonicExpBraid[Exp[\[Lambda] $n + \[Alpha] $ad] ** Exp[\[Gamma] $a]],
    Exp[\[Lambda] $n + \[Alpha] $ad] ** Exp[\[Gamma] $a],
    {},
    TestID -> "BosonicExpBraid-MixedFamilyIsInert"
]

EndTestSection[]


BeginTestSection["BosonicExpBraid - multimode"]

(* Per-mode number generators transport each mode with its own scaling. *)
VerificationTest[
    Module[{b = \[FormalB], bd = SuperDagger[\[FormalB]]},
        BosonicExpBraid[
            Exp[\[Lambda] $n + \[Mu] bd ** b] ** Exp[\[Alpha] $ad ** b]]
    ],
    Module[{b = \[FormalB], bd = SuperDagger[\[FormalB]]},
        Exp[\[Alpha] Exp[\[Lambda] - \[Mu]] $ad ** b] ** Exp[\[Lambda] $n + \[Mu] bd ** b]
    ],
    {},
    TestID -> "BosonicExpBraid-MultimodeNumberGenerator"
]

(* A cross-mode generator (beam splitter) is not a recognized family. *)
VerificationTest[
    Module[{b = \[FormalB], bd = SuperDagger[\[FormalB]]},
        BosonicExpBraid[Exp[\[Theta] $ad ** b] ** Exp[\[Alpha] $ad]]
    ],
    Module[{b = \[FormalB]},
        Exp[\[Theta] $ad ** b] ** Exp[\[Alpha] $ad]
    ],
    {},
    TestID -> "BosonicExpBraid-CrossModeGeneratorIsInert"
]

EndTestSection[]
