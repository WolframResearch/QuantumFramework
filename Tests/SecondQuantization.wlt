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
   explicit two-argument form.  Single mode: [a, a^dag] = 1 gives <0|a a^dag|0> = 1,
   and a annihilates the vacuum so <0|a^dag a|0> = 0. *)
VerificationTest[
    With[{a = \[FormalA], ad = SuperDagger[\[FormalA]]},
        {
            BosonicVEV[a ** ad],
            BosonicVEV[a ** ad] === BosonicVEV[a ** ad, {a, ad}],
            BosonicVEV[ad ** a]
        }
    ],
    {1, True, 0},
    TestID -> "SQ-BosonicVEV-single-mode-auto-detect"
]

(* Auto-detection does real work only with more than one mode.  Independent modes
   commute, so <0|a a^dag b b^dag|0> = 1, while a fully normal-ordered product has
   no c-number part: <0|a^dag a b^dag b|0> = 0.  Single-arg must match two-arg. *)
VerificationTest[
    With[{a = \[FormalA], ad = SuperDagger[\[FormalA]], b = \[FormalB], bd = SuperDagger[\[FormalB]]},
        {
            BosonicVEV[a ** ad ** b ** bd],
            BosonicVEV[a ** ad ** b ** bd] === BosonicVEV[a ** ad ** b ** bd, {a, ad, b, bd}],
            BosonicVEV[ad ** a ** bd ** b]
        }
    ],
    {1, True, 0},
    TestID -> "SQ-BosonicVEV-multimode-auto-detect"
]

(* Wick's theorem for the vacuum: <0|a^n (a^dag)^n|0> = n!.  Exercises the
   auto-detect path on a family, single-arg agreeing with two-arg.  n is kept
   concrete: the reduction does not accept a symbolic exponent. *)
VerificationTest[
    With[{a = \[FormalA], ad = SuperDagger[\[FormalA]]},
        With[{chains = Table[NonCommutativeMultiply @@ Join[ConstantArray[a, n], ConstantArray[ad, n]], {n, 5}]},
            {
                BosonicVEV /@ chains,
                (BosonicVEV[#] === BosonicVEV[#, {a, ad}] &) /@ chains
            }
        ]
    ],
    {{1, 2, 6, 24, 120}, {True, True, True, True, True}},
    TestID -> "SQ-BosonicVEV-Wick-n-factorial"
]

(* A symbolic scalar coefficient passes straight through the vacuum projection:
   <0| lambda (a a^dag) |0> = lambda. *)
VerificationTest[
    With[{a = \[FormalA], ad = SuperDagger[\[FormalA]]},
        {
            BosonicVEV[\[Lambda] (a ** ad)],
            BosonicVEV[\[Lambda] (a ** ad)] === BosonicVEV[\[Lambda] (a ** ad), {a, ad}]
        }
    ],
    {\[Lambda], True},
    TestID -> "SQ-BosonicVEV-symbolic-coefficient"
]

(* Independent route: the symbolic normal-ordering result matches the finite-cutoff
   Fock-space matrix element <0|a a^dag|0> built from the ladder operators. *)
VerificationTest[
    With[{a = \[FormalA], ad = SuperDagger[\[FormalA]]},
        BosonicVEV[a ** ad] ===
            (FockState[0]["Dagger"] @ (AnnihilationOperator[] @ AnnihilationOperator[]["Dagger"]) @ FockState[0])["Scalar"]
    ],
    True,
    TestID -> "SQ-BosonicVEV-matches-Fock-matrix-element"
]

(* Independent modes decouple.  [a, b] = [a, b^dag] = 0, so the VEV factorizes
   across modes and cross-mode interleaving leaves it unchanged: <0|a a^dag b b^dag|0>
   = <0|a a^dag|0> <0|b b^dag|0>, the interleaved <0|a b a^dag b^dag|0> = 1, and the
   reordered <0|a^dag b a b^dag|0> = 0. *)
VerificationTest[
    With[{a = \[FormalA], ad = SuperDagger[\[FormalA]], b = \[FormalB], bd = SuperDagger[\[FormalB]]},
        {
            BosonicVEV[a ** ad ** b ** bd] === BosonicVEV[a ** ad] BosonicVEV[b ** bd],
            BosonicVEV[a ** b ** ad ** bd],
            BosonicVEV[ad ** b ** a ** bd]
        }
    ],
    {True, 1, 0},
    TestID -> "SQ-BosonicVEV-modes-decouple"
]

(* Number-conservation selection rule: the vacuum expectation vanishes unless the
   creation and annihilation counts balance.  Unbalanced monomials give 0, while
   the balanced a a a^dag a^dag follows Wick, <0|a^2 (a^dag)^2|0> = 2. *)
VerificationTest[
    With[{a = \[FormalA], ad = SuperDagger[\[FormalA]]},
        {BosonicVEV[a ** a ** ad], BosonicVEV[a ** ad ** ad], BosonicVEV[a ** a ** ad ** ad]}
    ],
    {0, 0, 2},
    TestID -> "SQ-BosonicVEV-number-selection-rule"
]

(* A c-number has no field operators, so it is its own vacuum expectation:
   <0|c|0> = c. The single-argument auto-detect path finds an empty variable list
   and must return the c-number rather than recurse to the recursion limit; the
   explicit empty-list form BosonicVEV[c, {}] routes through the same rule. *)
VerificationTest[
    {BosonicVEV[7], BosonicVEV[\[Lambda]], BosonicVEV[7, {}]},
    {7, \[Lambda], 7},
    TestID -> "SQ-BosonicVEV-cnumber-is-itself"
]

(* Normal-ordering with no field variables is the identity on the c-number: with
   nothing to reduce, the reducer is never invoked and the constant passes through
   unchanged and message-free. *)
VerificationTest[
    {BosonicNormalOrder[7, {}], BosonicNormalOrder[\[Lambda], {}]},
    {7, \[Lambda]},
    TestID -> "SQ-BosonicNormalOrder-empty-vars-is-identity"
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

(* The origin value is the mean parity for ANY single-mode state,
   W(0) = (2/pi) Sum rho_nn (-1)^n, and Q(0) = rho00/pi. A superposition exercises
   the off-diagonal elements, which vanish at the origin, and the sign of the
   parity, neither of which a single Fock ket reaches: (|0>+|2>)/Sqrt[2] adds two
   even parities (rho02 drops out) while (|0>+|1>)/Sqrt[2] cancels to 0. *)
VerificationTest[
    WignerFunction[QuantumState[{1, 0, 1, 0, 0, 0, 0, 0}/Sqrt[2], 8], 0],
    2/Pi,
    TestID -> "SQ-Wigner-origin-superposition-even-parity"
]

VerificationTest[
    WignerFunction[QuantumState[{1, 1, 0, 0, 0, 0, 0, 0}/Sqrt[2], 8], 0],
    0,
    TestID -> "SQ-Wigner-origin-superposition-parity-cancels"
]

VerificationTest[
    HusimiQFunction[QuantumState[{1, 1, 0, 0, 0, 0, 0, 0}/Sqrt[2], 8], 0],
    1/(2 Pi),
    TestID -> "SQ-Husimi-origin-superposition-vacuum-weight"
]

(* All three SymbolicForm branches of WignerFunction are finite at the origin. *)
VerificationTest[
    WignerFunction[FockState[0, 8], 0, SymbolicForm -> "Wirtinger"],
    2/Pi,
    TestID -> "SQ-Wigner-origin-Wirtinger-form"
]

VerificationTest[
    Activate[WignerFunction[FockState[0, 8], 0, SymbolicForm -> "LaguerreForm"]],
    2/Pi,
    TestID -> "SQ-Wigner-origin-LaguerreForm"
]

(* The origin value is the alpha -> 0 limit: a tiny displacement reproduces it. *)
VerificationTest[
    Chop[WignerFunction[FockState[0, 20], 1.*^-6] - 2/Pi, 1.*^-8],
    0,
    TestID -> "SQ-Wigner-origin-is-the-limit"
]

(* Independent check of the fixed origin value against the numeric-grid route,
   which uses the Laguerre recurrence rather than the closed-form kernel. *)
VerificationTest[
    Chop[WignerFunction[FockState[0, 8], {0, 0}] -
        WignerRepresentation[FockState[0, 8], {-2, 2}, {-2, 2}][0, 0], 1.*^-4],
    0,
    TestID -> "SQ-Wigner-origin-matches-grid-route"
]

(* Husimi vacuum keeps its Gaussian closed form Q(alpha) = (1/pi) e^{-|alpha|^2}. *)
VerificationTest[
    Simplify[HusimiQFunction[FockState[0, 20], \[FormalAlpha]] - (1/Pi) Exp[-Abs[\[FormalAlpha]]^2]],
    0,
    TestID -> "SQ-Husimi-symbolic-form-preserved"
]

(* Multimode states fall outside the single-mode definition and guard to $Failed. *)
VerificationTest[
    WignerFunction[QuantumTensorProduct[FockState[0, 4], FockState[1, 4]], 0],
    $Failed,
    {WignerFunction::multimode},
    TestID -> "SQ-Wigner-multimode-guards-to-Failed"
]

(* Away from the origin a complex coherence separates alpha from its conjugate.
   For psi = (|0> + I |1>)/Sqrt[2] the Husimi function is Q(alpha) = (1/pi) |<alpha|psi>|^2;
   at alpha = 2/5 + 3/10 I this equals (37/(40 pi)) e^{-1/4}, distinct from the mirror
   Q(Conjugate[alpha]), so it pins the bra/ket ordering of the Husimi kernel. *)
VerificationTest[
    HusimiQFunction[QuantumState[{1, I, 0, 0, 0, 0, 0, 0}/Sqrt[2], 8], 2/5 + 3/10 I],
    (37/(40 Pi)) Exp[-1/4],
    TestID -> "SQ-Husimi-complex-coherence-pins-bra-ket-order"
]

(* A rank-2 mixed state with a complex coherence carries the fix past pure states:
   the origin gives the mean parity (2/pi) Sum rho_nn (-1)^n and rho00/pi, and
   off-origin the Husimi value is (1/pi) <alpha|rho|alpha>. *)
VerificationTest[
    With[{rho = QuantumState[{{2/3, -I/3, 0, 0}, {I/3, 1/3, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}, 4]},
        {WignerFunction[rho, 0], HusimiQFunction[rho, 0], HusimiQFunction[rho, 2/5 + 3/10 I]}
    ],
    {2/(3 Pi), 2/(3 Pi), (19/(20 Pi)) Exp[-1/4]},
    TestID -> "SQ-mixed-state-origin-and-complex-coherence"
]

(* General symbolic Wigner beyond the vacuum: the one-photon Fock state has
   W_1(alpha) = (2/pi)(-1) L_1(4|alpha|^2) e^{-2|alpha|^2} with L_1(x) = 1 - x. *)
VerificationTest[
    Simplify[WignerFunction[FockState[1, 20], \[FormalAlpha]] -
        (2/Pi) (-1) (1 - 4 Abs[\[FormalAlpha]]^2) Exp[-2 Abs[\[FormalAlpha]]^2]],
    0,
    TestID -> "SQ-Wigner-symbolic-Fock1-form"
]

EndTestSection[]
