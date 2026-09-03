Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

(* Local abbreviations. \[FormalA] is the annihilation operator, its
   SuperDagger the creation operator, n the number operator. *)
$a = \[FormalA];
$ad = SuperDagger[\[FormalA]];
$n = $ad ** $a;

(* Displacement generator \[Alpha] a\[Dagger] - Conjugate[\[Alpha]] a. *)
$disp = \[Alpha] $ad - Conjugate[\[Alpha]] $a;

(* Exponent of a single Exp, for comparing generators coefficient-wise:
   two linear combinations of the ladder operators are equal iff their
   difference collects to 0 under Plus/Times. *)
exponentOf[Exp[x_]] := x;
exponentOf[_] := $Failed;
sameGeneratorQ[res_, expected_] :=
    FullSimplify[exponentOf[res] - expected] === 0;


BeginTestSection["BosonicExpSimplify - adjacent pair reduction"]

(* [\[Alpha] a\[Dagger], \[Beta] a] = -\[Alpha]\[Beta] is a c-number, so the Weyl/BCH shortcut applies. *)
VerificationTest[
    BosonicExpSimplify[Exp[\[Alpha] $ad] ** Exp[\[Beta] $a]],
    Exp[-\[Alpha] \[Beta]/2] Exp[\[Alpha] $ad + \[Beta] $a],
    {},
    TestID -> "BosonicExpSimplify-WeylPairCentralCommutator"
]

(* Commuting exponents merge outright. Plus does not factor \[Beta] x + \[Gamma] x on its
   own, so compare the generators rather than their printed form. *)
VerificationTest[
    sameGeneratorQ[
        BosonicExpSimplify[Exp[\[Beta] $n] ** Exp[\[Gamma] $n]],
        (\[Beta] + \[Gamma]) $n
    ],
    True,
    {},
    TestID -> "BosonicExpSimplify-CommutingPairMerges"
]

(* An exponential times its inverse must collapse to the identity, not to a
   two-factor product. This is the pre-pass that In[77] needed. *)
VerificationTest[
    BosonicExpSimplify[Exp[\[Beta] $n] ** Exp[-\[Beta] $n]],
    1,
    {},
    TestID -> "BosonicExpSimplify-InversePairCancels"
]

(* Same cancellation buried inside a longer product: the pair rule has to see
   an adjacent window, not the whole expression. *)
VerificationTest[
    BosonicExpSimplify[
        Exp[\[Alpha] $ad] ** Exp[\[Beta] $n] ** Exp[-\[Beta] $n] ** Exp[\[Gamma] $a]],
    Exp[-\[Alpha] \[Gamma]/2] Exp[\[Alpha] $ad + \[Gamma] $a],
    {},
    TestID -> "BosonicExpSimplify-InversePairInsideLongerProduct"
]

EndTestSection[]


BeginTestSection["BosonicExpSimplify - conjugation"]

(* e^(\[Beta] n) D(\[Alpha]) e^(-\[Beta] n) = D'(\[Alpha] e^\[Beta]): the ad-action closes on span{a, a\[Dagger]}. *)
VerificationTest[
    sameGeneratorQ[
        FullSimplify @ BosonicExpSimplify[
            Exp[\[Beta] $n] ** Exp[$disp] ** Exp[-\[Beta] $n]],
        \[Alpha] Exp[\[Beta]] $ad - Conjugate[\[Alpha]] Exp[-\[Beta]] $a
    ],
    True,
    {},
    TestID -> "BosonicExpSimplify-ConjugationRotatesDisplacement"
]

(* Same identity at \[Beta] = -I \[Omega] t, where the result is a genuine unitary
   displacement D(\[Alpha] e^(-I \[Omega] t)) - the Heisenberg-picture rotation. *)
VerificationTest[
    sameGeneratorQ[
        FullSimplify[
            BosonicExpSimplify[
                Exp[\[Beta] $n] ** Exp[$disp] ** Exp[-\[Beta] $n]] /. \[Beta] -> -I \[Omega] t],
        \[Alpha] Exp[-I \[Omega] t] $ad - Conjugate[\[Alpha]] Exp[I \[Omega] t] $a
    ],
    True,
    {},
    TestID -> "BosonicExpSimplify-ConjugationPhaseRotation"
]

(* Several factors inside one sandwich conjugate independently. Here the
   grouped middle a**a has no closed ad-action (\[Lambda] would be operator-valued),
   so only the per-factor distribution can finish it:
   e^(\[Alpha] a\[Dagger]) a a e^(-\[Alpha] a\[Dagger]) = (a - \[Alpha])(a - \[Alpha]). *)
VerificationTest[
    BosonicExpSimplify[Exp[\[Alpha] $ad] ** $a ** $a ** Exp[-\[Alpha] $ad]],
    ($a - \[Alpha]) ** ($a - \[Alpha]),
    {},
    TestID -> "BosonicExpSimplify-ConjugationDistributesOverFactors"
]

(* A middle that does close as a whole is handled in one step, since a
   Flat pattern lets the single-factor rule absorb the run:
   e^(\[Beta] n) a\[Dagger] a e^(-\[Beta] n) = a\[Dagger] a. *)
VerificationTest[
    BosonicExpSimplify[Exp[\[Beta] $n] ** $ad ** $a ** Exp[-\[Beta] $n]],
    $ad ** $a,
    {},
    TestID -> "BosonicExpSimplify-ConjugationOfNumberOperatorIsFixed"
]

EndTestSection[]


BeginTestSection["BosonicExpSimplify - four-factor gap"]

(* In[77]: the trailing Exp[\[Beta] n] cancels the closing Exp[-\[Beta] n], so this is
   In[76] with a redundant factor. It must reduce to the same two-factor answer,
   and must emit no messages while doing it. *)
VerificationTest[
    sameGeneratorQ[
        FullSimplify @ First @ BosonicExpSimplify[
            Exp[\[Beta] $n] ** Exp[$disp] ** Exp[-\[Beta] $n] ** Exp[\[Beta] $n]],
        \[Alpha] Exp[\[Beta]] $ad - Conjugate[\[Alpha]] Exp[-\[Beta]] $a
    ],
    True,
    {},
    TestID -> "BosonicExpSimplify-FourFactorFirstFactorIsRotatedDisplacement"
]

VerificationTest[
    Last @ BosonicExpSimplify[
        Exp[\[Beta] $n] ** Exp[$disp] ** Exp[-\[Beta] $n] ** Exp[\[Beta] $n]],
    Exp[\[Beta] $n],
    {},
    TestID -> "BosonicExpSimplify-FourFactorSecondFactorIsNumberExponential"
]

VerificationTest[
    Length @ BosonicExpSimplify[
        Exp[\[Beta] $n] ** Exp[$disp] ** Exp[-\[Beta] $n] ** Exp[\[Beta] $n]],
    2,
    {},
    TestID -> "BosonicExpSimplify-FourFactorReducesToTwoFactors"
]

EndTestSection[]


BeginTestSection["BosonicExpSimplify - declines silently"]

(* [\[Beta] n, \[Alpha] a\[Dagger] - \[Alpha]\[Conjugate] a] is operator-valued, so no closed form is
   claimed. The input must come back untouched and, crucially, with no
   NonCommutativePolynomialReduce::ncpv leaking out of the internals. *)
VerificationTest[
    BosonicExpSimplify[Exp[\[Beta] $n] ** Exp[$disp]],
    Exp[\[Beta] $n] ** Exp[$disp],
    {},
    TestID -> "BosonicExpSimplify-NonCentralPairReturnsInputSilently"
]

VerificationTest[
    BosonicExpSimplify[Exp[\[Beta] $n] ** Exp[$disp] ** Exp[\[Beta] $n]],
    Exp[\[Beta] $n] ** Exp[$disp] ** Exp[\[Beta] $n],
    {},
    TestID -> "BosonicExpSimplify-NonInverseSandwichReturnsInputSilently"
]

EndTestSection[]


BeginTestSection["BosonicExpSimplify - polynomial past exponential"]

(* a e^(\[Alpha] a\[Dagger]) = e^(\[Alpha] a\[Dagger]) (a + \[Alpha]). *)
VerificationTest[
    BosonicExpSimplify[$a ** Exp[\[Alpha] $ad]],
    Exp[\[Alpha] $ad] ** (\[Alpha] + $a),
    {},
    TestID -> "BosonicExpSimplify-PolynomialPastCreationExponential"
]

EndTestSection[]


BeginTestSection["BosonicExpSimplify - numerical validation"]

(* Truncated-oscillator check of the four-factor result: compare the top block
   of the symbolic answer against the direct matrix product. *)
VerificationTest[
    Module[{d = 60, keep = 12, av, ad, nOp, num, expo, lhs, rhs, blk},
        av = SparseArray[Band[{1, 2}] -> Sqrt[Range[d - 1]], {d, d}];
        ad = Transpose[av];
        nOp = ad . av;
        num = {\[Beta] -> 0.3 + 0.2 I, \[Alpha] -> 0.5 - 0.4 I};
        (* The trailing factor cancels the closing one, so the input equals
           e^(\[Beta] n) D(\[Alpha]) as a matrix product. *)
        lhs = MatrixExp[(\[Beta] /. num) nOp] .
            MatrixExp[(\[Alpha] /. num) ad - Conjugate[\[Alpha] /. num] av];
        (* Scalars before matrices, so MatrixExp only sees a numeric matrix. *)
        expo = FullSimplify @ exponentOf @ First @ BosonicExpSimplify[
            Exp[\[Beta] $n] ** Exp[$disp] ** Exp[-\[Beta] $n] ** Exp[\[Beta] $n]];
        rhs = MatrixExp[expo /. num /. {$ad -> ad, $a -> av}] .
            MatrixExp[(\[Beta] /. num) nOp];
        blk = Normal[#][[;; keep, ;; keep]] &;
        Max[Abs[blk[lhs] - blk[rhs]]] < 10.^-10
    ],
    True,
    {},
    TestID -> "BosonicExpSimplify-FourFactorNumericalValidation"
]

EndTestSection[]


BeginTestSection["BosonicExpSimplify - Assumptions"]

(* Squeeze generator 1/2 Conjugate[z] a a - 1/2 z a\[Dagger] a\[Dagger]. *)
$squeeze[z_] := Exp[1/2 Conjugate[z] $a ** $a - 1/2 z $ad ** $ad];

(* The su(1,1) closure runs Simplify internally; without assumptions it cannot
   reduce Conjugate[r E^(I \[Theta])] and leaves the rotation angle unresolved. Passing
   the assumptions in gives the Bogoliubov form in a single call. *)
VerificationTest[
    BosonicExpSimplify[
        $squeeze[-r E^(I \[Theta])] ** $a ** $squeeze[r E^(I \[Theta])],
        Assumptions -> r > 0 && \[Theta] \[Element] Reals
    ],
    $a Cosh[r] - E^(I \[Theta]) Sinh[r] $ad,
    {},
    TestID -> "BosonicExpSimplify-AssumptionsSqueezedAnnihilation"
]

VerificationTest[
    BosonicExpSimplify[
        $squeeze[-r E^(I \[Theta])] ** $ad ** $squeeze[r E^(I \[Theta])],
        Assumptions -> r > 0 && \[Theta] \[Element] Reals
    ],
    $ad Cosh[r] - E^(-I \[Theta]) Sinh[r] $a,
    {},
    TestID -> "BosonicExpSimplify-AssumptionsSqueezedCreation"
]

(* The option must survive the explicit-vars overload alongside "Scalars". *)
VerificationTest[
    BosonicExpSimplify[
        $squeeze[-r E^(I \[Theta])] ** $a ** $squeeze[r E^(I \[Theta])],
        {$a, $ad},
        "Scalars" -> {r, \[Theta]},
        Assumptions -> r > 0 && \[Theta] \[Element] Reals
    ],
    $a Cosh[r] - E^(I \[Theta]) Sinh[r] $ad,
    {},
    TestID -> "BosonicExpSimplify-AssumptionsWithExplicitVars"
]

(* Assumptions are additive only: with none supplied the result must still be
   the unassumed closed form, equal to the assumed one under Simplify. *)
VerificationTest[
    Simplify[
        BosonicExpSimplify[$squeeze[-r E^(I \[Theta])] ** $a ** $squeeze[r E^(I \[Theta])]] -
            ($a Cosh[r] - E^(I \[Theta]) Sinh[r] $ad),
        r > 0 && \[Theta] \[Element] Reals
    ],
    0,
    {},
    TestID -> "BosonicExpSimplify-NoAssumptionsUnchanged"
]

EndTestSection[]
