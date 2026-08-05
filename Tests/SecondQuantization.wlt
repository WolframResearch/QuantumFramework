(* Regression tests for the SecondQuantization context, which had none: its
   operators exercise heavy composition - a displacement is a MatrixExp of
   ladder operators - so it doubles as a sweep of operator application. *)

Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

BeginTestSection["SecondQuantization"]

SetFockSpaceSize[8]

VerificationTest[
    {$FockSize, FockState[2]["Dimension"], Normal[FockState[2]["ProbabilitiesList"]]},
    {8, 8, UnitVector[8, 3]},
    TestID -> "SQ-Fock-state-basics"
]

VerificationTest[
    Chop[Normal[(AnnihilationOperator[] @ FockState[2])["StateVector"]] -
        Sqrt[2] Normal[FockState[1]["StateVector"]]],
    ConstantArray[0, 8],
    TestID -> "SQ-annihilation-lowers-with-sqrt-n"
]

(* A coherent state has mean photon number and variance both |alpha|^2. *)
VerificationTest[
    With[{cs = CoherentState[][0.5], n = AnnihilationOperator[]["Dagger"] @ AnnihilationOperator[]},
        {
            Chop[(cs["Dagger"] @ n @ cs)["Scalar"] - 0.25, 1*^-6],
            Chop[OperatorVariance[cs, n] - 0.25, 1*^-6]
        }
    ],
    {0, 0},
    TestID -> "SQ-coherent-state-Poisson-statistics"
]

VerificationTest[
    Chop[(DisplacementOperator[0.3] @ FockState[0])["Norm"] - 1, 1*^-4],
    0,
    TestID -> "SQ-displacement-preserves-the-norm"
]

VerificationTest[
    Chop[G2Coherence[CoherentState[][0.5]] - 1, 1*^-2],
    0,
    TestID -> "SQ-coherent-g2-is-one"
]

(* BosonicVEV[expr] auto-detects the field variables and must agree with the
   explicit two-argument form.  <0|a a^dag|0> = 1 because [a, a^dag] = 1, and
   <0|a^dag a|0> = 0 because a annihilates the vacuum. *)
VerificationTest[
    With[{a = \[FormalA], ad = SuperDagger[\[FormalA]]},
        {
            BosonicVEV[a ** ad],
            BosonicVEV[a ** ad] === BosonicVEV[a ** ad, {a, ad}],
            BosonicVEV[ad ** a]
        }
    ],
    {1, True, 0},
    TestID -> "SQ-BosonicVEV-single-arg-auto-detects-variables"
]

(* Phase-space functions at the origin alpha = 0. The closed-form kernels carry a
   base^(m - n) (Wigner) or base^n base^m (Husimi) factor; on the diagonal the
   exponent is 0, so the origin value is the limit base^0 = 1. W(0) = (2/pi) <parity>
   and Q(0) = rho00/pi: the vacuum gives 2/pi and 1/pi, Fock |1> gives -2/pi (parity
   -1) and 0. *)

VerificationTest[
    WignerFunction[FockState[0, 8], 0],
    2/Pi,
    TestID -> "SQ-Wigner-origin-vacuum-is-2overPi"
]

VerificationTest[
    HusimiQFunction[FockState[0, 8], 0],
    1/Pi,
    TestID -> "SQ-Husimi-origin-vacuum-is-1overPi"
]

VerificationTest[
    WignerFunction[FockState[1, 8], 0],
    -2/Pi,
    TestID -> "SQ-Wigner-origin-Fock1-is-negative-parity"
]

VerificationTest[
    HusimiQFunction[FockState[1, 8], 0],
    0,
    TestID -> "SQ-Husimi-origin-Fock1-vanishes"
]

(* The {x, p} quadrature entry points route through the alpha form, so they too are
   finite at the origin. *)
VerificationTest[
    WignerFunction[FockState[0, 8], {0, 0}],
    1/Pi,
    TestID -> "SQ-Wigner-origin-xp-entry"
]

VerificationTest[
    HusimiQFunction[FockState[0, 8], {0, 0}],
    1/(2 Pi),
    TestID -> "SQ-Husimi-origin-xp-entry"
]

(* Away from the origin the closed forms are untouched: numeric agreement with the
   Gaussian vacuum W(alpha) = (2/pi) e^{-2 |alpha|^2}, and the symbolic form itself. *)
VerificationTest[
    Chop[WignerFunction[FockState[0, 20], 0.5] - (2/Pi) Exp[-2 (0.5)^2], 1*^-10],
    0,
    TestID -> "SQ-Wigner-offorigin-numeric-preserved"
]

VerificationTest[
    Simplify[WignerFunction[FockState[0, 20], \[FormalAlpha]] - (2/Pi) Exp[-2 Abs[\[FormalAlpha]]^2]],
    0,
    TestID -> "SQ-Wigner-symbolic-form-preserved"
]

EndTestSection[]
