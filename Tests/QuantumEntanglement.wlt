(* ::Package:: *)

BeginTestSection["QuantumEntanglement"]


(* ========== Setup ========== *)

conc = QuantumEntanglementMonotone[#, "Concurrence"] &;

bell    = QuantumState["PhiPlus"];                                              (* (|00>+|11>)/Sqrt2 *)
bellDM  = QuantumState[{{1, 0, 0, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 0, 0, 1}} / 2, {2, 2}];
prod    = QuantumState["00"];
sep     = QuantumState[DiagonalMatrix[{2, 1, 3, 2} / 8], {2, 2}];              (* diagonal => separable *)

werner[p_] := p {{1, 0, 0, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 0, 0, 1}} / 2 + (1 - p) / 4 IdentityMatrix[4];

(* Symbolic families with known closed-form concurrence. These are the regression guards:
   the ConcurrenceVector bug was symbolic-only (every numeric value was already correct), because
   Fold[Subtract] @ SingularValueList assumed a decreasing singular-value order that WL guarantees
   for numeric matrices but not for symbolic ones. A numeric test cannot catch it; a symbolic one can. *)
rhoLinear = {{q/2, 0, 0, q/2}, {0, 1 - q, 0, 0}, {0, 0, 0, 0}, {q/2, 0, 0, q/2}};   (* C = q            *)
rhoClamp  = {{(1 - q)/2, 0, 0, 0}, {0, q/2, q/2, 0}, {0, q/2, q/2, 0}, {0, 0, 0, (1 - q)/2}}; (* C = Max[0, 2q-1] *)

(* Independent textbook Wootters via a DIFFERENT route (eigenvalues of rho.rho~, not QF's
   singular values of Sqrt[rho].Sqrt[rho~]); no QuantumFramework symbols. *)
wootters[r_] := Module[{Y2 = KroneckerProduct[PauliMatrix[2], PauliMatrix[2]], ev},
    ev = Sqrt[Max[#, 0] & /@ ReverseSort[Re @ Eigenvalues[r . Y2 . Conjugate[r] . Y2]]];
    Max[0, ev[[1]] - ev[[2]] - ev[[3]] - ev[[4]]]
];


(* ========== Concurrence: pure states (closed form) ========== *)

(* Bell state is maximally entangled *)
VerificationTest[
    conc[bell],
    1,
    TestID -> "Concurrence-BellPure"
]

(* product state is separable *)
VerificationTest[
    conc[prod],
    0,
    TestID -> "Concurrence-ProductZero"
]

(* cos t |00> + sin t |11> has concurrence |sin 2t| *)
VerificationTest[
    FullSimplify[conc[QuantumState[{Cos[t], 0, 0, Sin[t]}, {2, 2}]], 0 < t < Pi/2],
    Sin[2 t],
    TestID -> "Concurrence-PureClosedForm"
]

(* pure-state concurrence equals the reduced-purity form Sqrt[2 (1 - Tr rhoA^2)] *)
VerificationTest[
    FullSimplify[
        conc[QuantumState[{Cos[t], 0, 0, Sin[t]}, {2, 2}]] -
            Sqrt[2 (1 - Tr[MatrixPower[Normal @ QuantumPartialTrace[QuantumState[{Cos[t], 0, 0, Sin[t]}, {2, 2}], {2}]["DensityMatrix"], 2]])],
        0 < t < Pi/2
    ],
    0,
    TestID -> "Concurrence-PureReducedPurityIdentity"
]


(* ========== Concurrence: SYMBOLIC mixed states (regression for the SingularValueList ordering bug) ========== *)

(* q |Phi+><Phi+| + (1-q) |01><01| has C = q. Returned 0 before the 2 Max - Total fix. *)
VerificationTest[
    FullSimplify[conc[QuantumState[rhoLinear, {2, 2}]], 0 < q < 1],
    q,
    TestID -> "Concurrence-SymbolicMixed-Linear"
]

VerificationTest[
    Limit[FullSimplify[conc[QuantumState[rhoLinear, {2, 2}]], 0 < q < 1], q -> 0],
    0,
    TestID -> "Concurrence-SymbolicMixed-SeparableLimit"
]

VerificationTest[
    Limit[FullSimplify[conc[QuantumState[rhoLinear, {2, 2}]], 0 < q < 1], q -> 1],
    1,
    TestID -> "Concurrence-SymbolicMixed-BellLimit"
]

(* Bell-diagonal q |psi+><psi+| + (1-q)/2 (|00><00|+|11><11|) has C = Max[0, 2q-1]:
   the symbolic result must fire the separability clamp, giving 2q-1 above q=1/2 and exactly 0 below. *)
VerificationTest[
    FullSimplify[conc[QuantumState[rhoClamp, {2, 2}]], 1/2 < q < 1],
    2 q - 1,
    TestID -> "Concurrence-SymbolicClamp-Entangled"
]

VerificationTest[
    FullSimplify[conc[QuantumState[rhoClamp, {2, 2}]], 0 < q < 1/2],
    0,
    TestID -> "Concurrence-SymbolicClamp-Separable"
]


(* ========== Concurrence: independent Wootters cross-check (different computational route) ========== *)

VerificationTest[
    SeedRandom[20260803];
    Max @ Table[
        With[{r = QuantumState["RandomMixed", {2, 2}]},
            Abs[conc[r] - wootters[Normal @ r["DensityMatrix"]]]
        ],
        {12}
    ] < 10^-8,
    True,
    TestID -> "Concurrence-IndependentWootters"
]


(* ========== Concurrence: physics invariants ========== *)

(* bounded in [0, 1] *)
VerificationTest[
    SeedRandom[7];
    AllTrue[Table[conc[QuantumState["RandomMixed", {2, 2}]], {12}], 0 <= # <= 1 &],
    True,
    TestID -> "Concurrence-Bounded"
]

(* invariant under local unitaries U1 (x) U2: entanglement cannot change under local operations *)
VerificationTest[
    SeedRandom[3];
    Module[{r = Normal @ QuantumState["RandomMixed", {2, 2}]["DensityMatrix"], u},
        u = KroneckerProduct[MatrixExp[I 0.7 PauliMatrix[2]], MatrixExp[I 1.1 PauliMatrix[1]]];
        Chop[conc[QuantumState[u . r . ConjugateTranspose[u], {2, 2}]] - conc[QuantumState[r, {2, 2}]]]
    ],
    0,
    TestID -> "Concurrence-LocalUnitaryInvariant"
]


(* ========== Werner state: entanglement sudden death ========== *)

(* C(Werner p) = Max[0, (3p-1)/2]: exactly 0 at and below the p = 1/3 threshold, positive above *)
VerificationTest[
    Chop[conc[QuantumState[N @ werner[1/2], {2, 2}]] - 1/4],
    0,
    TestID -> "Werner-Above"
]

VerificationTest[
    conc[QuantumState[N @ werner[3/10], {2, 2}]],
    0,
    TestID -> "Werner-BelowThreshold"
]

VerificationTest[
    {
        QuantumEntangledQ[QuantumState[N @ werner[4/5], {2, 2}], Automatic, "Concurrence"],
        QuantumEntangledQ[QuantumState[N @ werner[1/5], {2, 2}], Automatic, "Concurrence"]
    },
    {True, False},
    TestID -> "Werner-EntangledQ-TracksConcurrence"
]


(* ========== Higher dimensional (d > 2) ========== *)

(* maximally-entangled qutrit pair, as a density matrix so it routes through ConcurrenceVector:
   I-concurrence Sqrt[2 (d-1)/d] = 2/Sqrt[3] for d = 3. The exact input is kept exact, so compare the
   difference numerically (an exact-form difference need not be structurally zero). *)
VerificationTest[
    Chop[
        N[conc[QuantumState[QuantumState[Normalize[Flatten[IdentityMatrix[3]]], {3, 3}]["DensityMatrix"], {3, 3}]] -
            2 / Sqrt[3]]
    ],
    0,
    TestID -> "Concurrence-QutritMaxEntangled"
]


(* ========== Higher dimensional (d > 2) MIXED: I-concurrence semantics + conditioning guards ==========

   For d > 2 the value is the Rungta-Buzek-Caves-Hillery-Milburn I-concurrence (PRA 64, 042315 (2001)):
   exact = Sqrt[2 (1 - Tr rhoA^2)] for a pure state of any dimension, and a LOWER BOUND on the
   convex-roof concurrence (which has no closed form) for a mixed state. The mixed d1 d2 > 4 path reads
   the singular values off Eigenvalues of rho.rho-tilde, keeping native precision (exact input stays
   exact); the retired operator square-root route leaked RowReduce::luc on ill-conditioned machine
   reductions and blew up to unsimplified Root objects on exact input. *)

(* I-concurrence of a normalized bipartite d x d vector via the reduced-state purity: an independent
   reference that does not touch ConcurrenceVector, exact for pure states of any dimension. *)
iConcurrence[v_, d_] := Sqrt[2 Max[0, Re[1 - Tr[MatrixPower[# . ConjugateTranspose[#] &[ArrayReshape[v, {d, d}]], 2]]]]];

(* pure qutrit pair (not maximally entangled), handed in as a density matrix so it routes through the
   mixed ConcurrenceVector path: C = Sqrt[2 (1 - Tr rhoA^2)] *)
VerificationTest[
    SeedRandom[314];
    With[{psi = Normalize[RandomComplex[{-1 - I, 1 + I}, 9]]},
        Abs[conc[QuantumState[QuantumState[psi, {3, 3}]["DensityMatrix"], {3, 3}]] - iConcurrence[psi, 3]] < 10^-6
    ],
    True,
    TestID -> "Concurrence-QutritPure-IConcurrenceIdentity"
]

(* ill-conditioned d>2 mixed (seed 99, smallest reduced eigenvalue ~3*10^-3) used to leak RowReduce::luc
   from the internal MatrixPower square root; Check returns $Failed if any message fires, so a re-leak
   fails the test. The value must be a clean, real, non-negative number. *)
VerificationTest[
    SeedRandom[99];
    Check[
        With[{r = conc[QuantumState["RandomMixed", {3, 3}]]}, NumericQ[r] && r >= 0],
        $Failed
    ],
    True,
    TestID -> "Concurrence-QutritMixed-NoMessageLeak"
]

(* exact-rational isotropic state 1/2 |Omega><Omega| + 1/2 I/9: the retired operator square-root route
   blew up into a huge unsimplified Root expression here (the exact matrix square roots never simplify),
   while the eigenvalue route keeps it exact AND tractable, returning the closed form 2/(3 Sqrt[3]).
   Assert it stays exact (Precision Infinity), equals the closed form, and computes without blowing up. *)
VerificationTest[
    With[{rho = 1/2 Outer[Times, #, Conjugate[#]] &[Normalize[{1, 0, 0, 0, 1, 0, 0, 0, 1}]] + IdentityMatrix[9] / 18},
        With[{r = TimeConstrained[conc[QuantumState[rho, {3, 3}]], 20, $TimedOut]},
            r =!= $TimedOut && Precision[r] === Infinity && Chop[N[r - 2 / (3 Sqrt[3])]] == 0
        ]
    ],
    True,
    TestID -> "Concurrence-QutritIsotropic-Exact"
]

(* lower-bound semantics: the I-concurrence never exceeds Sum p_i C_I(psi_i) of an explicit pure
   decomposition (an upper bound on the convex roof), and is non-negative *)
VerificationTest[
    SeedRandom[2718];
    Module[{ps = {2, 3, 5} / 10, psis, rho, cQF, bound},
        psis = Table[Normalize[RandomComplex[{-1 - I, 1 + I}, 9]], {3}];
        rho = Sum[ps[[i]] Outer[Times, psis[[i]], Conjugate[psis[[i]]]], {i, 3}];
        cQF = conc[QuantumState[rho, {3, 3}]];
        bound = Total[MapThread[#1 iConcurrence[#2, 3] &, {ps, psis}]];
        0 <= cQF <= bound + 10^-6
    ],
    True,
    TestID -> "Concurrence-QutritMixed-LowerBound"
]

(* the numeric d>2 eigenvalue route reproduces the operator-square-root (Wootters) route it replaced,
   computed here through a message-free manual Hermitian square root: both are the same RBCHM lower
   bound and must agree (this is the tight value guard the 2-qubit-only Wootters cross-check cannot give) *)
VerificationTest[
    SeedRandom[2];
    Module[{rho = N @ Normal @ QuantumState["RandomMixed", {3, 3}]["DensityMatrix"], hSqrt, yGen, ref},
        hSqrt[m_] := With[{es = Eigensystem[m]},
            Transpose[es[[2]]] . DiagonalMatrix[Sqrt[Clip[Re[es[[1]]], {0, Infinity}]]] . Conjugate[es[[2]]]
        ];
        yGen[n_] := Catenate @ Table[SparseArray[{{j, k} -> -I, {k, j} -> I}, {n, n}], {k, 2, n}, {j, k - 1}];
        ref = Norm @ Flatten @ Outer[
            Function[{ya, yb},
                With[{o = Normal @ KroneckerProduct[ya, yb]},
                    With[{sv = SingularValueList[hSqrt[rho] . hSqrt[o . Conjugate[rho] . o]]},
                        Max[0, 2 Max[sv] - Total[sv]]
                    ]
                ]
            ],
            yGen[3], yGen[3], 1
        ];
        Abs[conc[QuantumState[rho, {3, 3}]] - ref] < 10^-8
    ],
    True,
    TestID -> "Concurrence-QutritMixed-MatchesOperatorSqrtRoute"
]

(* the d>2 mixed I-concurrence is FRAME-DEPENDENT: unlike the two-qubit Wootters concurrence it is NOT
   invariant under local unitaries U1 (x) U2 (it rises or falls toward the true, LU-invariant convex-roof
   value it bounds). It stays a valid lower bound in every local frame, so a local unitary never pushes
   it above the decomposition bound Sum p_i C_I(psi_i) of the original pure ensemble. c0 > 0 keeps the
   check non-vacuous (the base state is genuinely entangled). *)
VerificationTest[
    SeedRandom[2718];
    Module[{ps = {2, 3, 5} / 10, psis, rho, bound, u, c0, c1},
        psis = Table[Normalize[RandomComplex[{-1 - I, 1 + I}, 9]], {3}];
        rho = Sum[ps[[i]] Outer[Times, psis[[i]], Conjugate[psis[[i]]]], {i, 3}];
        bound = Total[MapThread[#1 iConcurrence[#2, 3] &, {ps, psis}]];
        u = KroneckerProduct[
            MatrixExp[I (# + ConjugateTranspose[#] &[RandomComplex[{-1 - I, 1 + I}, {3, 3}]])],
            MatrixExp[I (# + ConjugateTranspose[#] &[RandomComplex[{-1 - I, 1 + I}, {3, 3}]])]
        ];
        c0 = conc[QuantumState[rho, {3, 3}]];
        c1 = conc[QuantumState[u . rho . ConjugateTranspose[u], {3, 3}]];
        (* c0 > 0: base state genuinely entangled; Abs[c0 - c1] > 0: the value actually moves under the
           local unitary (frame-dependence, which an LU-invariant reimplementation would break); c1 <= bound:
           it stays a valid lower bound on the LU-invariant convex roof in the rotated frame *)
        c0 > 10^-3 && Abs[c0 - c1] > 10^-6 && 0 <= c1 <= bound + 10^-6
    ],
    True,
    TestID -> "Concurrence-QutritMixed-FrameDependentLowerBound"
]

(* separable product d>2 mixed rhoA (x) rhoB => 0 *)
VerificationTest[
    SeedRandom[9];
    With[{
        rhoA = # / Tr[#] &[# . ConjugateTranspose[#] &[RandomComplex[{-1 - I, 1 + I}, {3, 3}]]],
        rhoB = # / Tr[#] &[# . ConjugateTranspose[#] &[RandomComplex[{-1 - I, 1 + I}, {3, 3}]]]
    },
        Chop[conc[QuantumState[KroneckerProduct[rhoA, rhoB], {3, 3}]], 10^-6]
    ],
    0,
    TestID -> "Concurrence-QutritProduct-Zero"
]


(* ========== Cross-monotone consistency ========== *)

(* Bell state: negativity 1/2, entanglement entropy 1 bit *)
VerificationTest[
    QuantumEntanglementMonotone[bell, "Negativity"],
    1/2,
    TestID -> "Negativity-Bell"
]

VerificationTest[
    QuantityMagnitude @ QuantumEntanglementMonotone[bell, "EntanglementEntropy"],
    1,
    TestID -> "EntanglementEntropy-Bell"
]

(* separable state: every monotone vanishes *)
VerificationTest[
    Chop[{conc[sep], QuantumEntanglementMonotone[sep, "Negativity"]}],
    {0, 0},
    TestID -> "Separable-AllMonotonesZero"
]

(* Bell prepared as a density matrix (mixed path) agrees with the pure (vector) path *)
VerificationTest[
    Chop[conc[bellDM] - 1],
    0,
    TestID -> "Concurrence-BellDensityPathMatchesVector"
]


(* ========== EntanglementEntropy: SYMBOLIC pure states ========== *)


ent = QuantityMagnitude @ QuantumEntanglementMonotone[#, "EntanglementEntropy"] &;

(* the standard entangler: RY(2t) then CNOT carries |00> to Cos[t] |00> + Sin[t] |11> *)
psiTheta = QuantumCircuitOperator[{"RY"[2 \[Theta]] -> 1, "CNOT" -> {1, 2}}][prod];

psiQutrit = {Cos[\[Theta]], Sin[\[Theta]] Cos[\[Phi]], 
   Sin[\[Theta]] Sin[\[Phi]]} . (QuantumState[#, {3, 3}] & /@ {"00", "11", "22"})

(* cos^2 and sin^2 weights*)
VerificationTest[
 Simplify[
  ent[psiTheta] - (-Cos[\[Theta]]^2 Log2[Cos[\[Theta]]^2] - 
     Sin[\[Theta]]^2 Log2[Sin[\[Theta]]^2]), 
  0 < \[Theta] < \[Pi]/2], 0, 
 TestID -> "EntanglementEntropy-SymbolicPure-ClosedForm"]

(* the same state as a density matrix takes the mixed branch, which was always right; the two must agree *)
VerificationTest[
 Simplify[
  ent[QuantumState[
     psiTheta[
      "MatrixState"]]] - (-Cos[\[Theta]]^2 Log2[Cos[\[Theta]]^2] - 
     Sin[\[Theta]]^2 Log2[Sin[\[Theta]]^2]), 
  0 < \[Theta] < \[Pi]/2], 0, 
 TestID -> "EntanglementEntropy-SymbolicVectorMatchesDensity"]

(* not qubit-specific: three symbolic Schmidt weights on a qutrit pair *)
VerificationTest[
 Simplify[
  ent[psiQutrit] - (-Cos[\[Theta]]^2 Log2[
       Cos[\[Theta]]^2] - (Sin[\[Theta]] Cos[\[Phi]])^2 Log2[(Sin[\
\[Theta]] Cos[\[Phi]])^2] - (Sin[\[Theta]] Sin[\[Phi]])^2 Log2[(Sin[\
\[Theta]] Sin[\[Phi]])^2]),
  0 < \[Theta] < Pi/2 && 0 < \[Phi] < 2 Pi], 0,
 TestID -> "EntanglementEntropy-SymbolicQutritPure-ClosedForm"]

(* the numeric side of the guard, on the same vector branch: weights 3/4 and 1/4. EntanglementEntropy-Bell
   pins that branch only at the symmetric 1/2, 1/2 point, where an ordering or normalization slip cancels. *)
psiUneven = QuantumCircuitOperator[{"RY"[Pi/3] -> 1, "CNOT" -> {1, 2}}][prod];

VerificationTest[
 Simplify[ent[psiUneven] - (2 - 3 Log2[3]/4)], 0,
 TestID -> "EntanglementEntropy-ExactUnequalWeights"]


(* ========== Normalization contract: every monotone reads the direction, not the scale ==========
   Entanglement is a property of the ray, so an input whose trace (mixed) or vector norm (pure) is not
   1 must give the same value as its normalized form. The concurrence path read 0 on a scaled pure state
   until its bipartition was normalized like the sibling monotones already were; these guard the shared
   contract. *)

ubell = QuantumState[2 bell["StateVector"]];   (* unnormalized Bell, vector norm 2 *)

(* pure: the reduced-purity route clamps a scaled 2 (1 - Purity) to 0 without normalization *)
VerificationTest[
    QuantumEntanglementMonotone[ubell, "Concurrence"],
    1,
    TestID -> "Concurrence-UnnormalizedPureScaleInvariant"
]

VerificationTest[
    QuantumEntanglementMonotone[ubell, "ConcurrenceVector"],
    QuantumEntanglementMonotone[bell, "ConcurrenceVector"],
    TestID -> "ConcurrenceVector-ScaleInvariant"
]

(* mixed: a Werner(4/5) scaled by 3 (trace 3) reads the same as the normalized state, and every mixed-state
   monotone in the family agrees on the scaled vs normalized input *)
VerificationTest[
    With[{w = QuantumState[N @ werner[4/5], {2, 2}], wu = QuantumState[3 N @ werner[4/5], {2, 2}]},
        Chop[
            Function[m, QuantumEntanglementMonotone[wu, m] - QuantumEntanglementMonotone[w, m]] /@
                {"Concurrence", "Negativity", "LogNegativity", "Realignment"}
        ]
    ],
    {0, 0, 0, 0},
    TestID -> "Monotones-ScaleInvariantAcrossFamily"
]

(* entanglement entropy (pure only) reads the same for the scaled and normalized pure state *)
VerificationTest[
    QuantityMagnitude @ QuantumEntanglementMonotone[ubell, "EntanglementEntropy"],
    1,
    TestID -> "EntanglementEntropy-UnnormalizedPureScaleInvariant"
]


(* ========== Negativity / LogNegativity / Realignment on the Bell state (value checks) ========== *)

VerificationTest[
    QuantumEntanglementMonotone[bell, "LogNegativity"],
    1,
    TestID -> "LogNegativity-Bell"
]

(* CCNR realigned Bell has trace norm 2, so the criterion value is 2 - 1 = 1 *)
VerificationTest[
    Chop[QuantumEntanglementMonotone[bell, "Realignment"] - 1],
    0,
    TestID -> "Realignment-Bell"
]


(* ========== Separability criteria: PPT (negativity) and CCNR (realignment) are complementary ==========
   Each QuantumEntanglementMonotone-backed criterion is one-directional: a positive value certifies
   entanglement, a non-positive one does not certify separability. For 2x2 and 2x3 the negativity/PPT
   value is exact (necessary and sufficient), while realignment/CCNR is only sufficient; in higher
   dimensions realignment catches positive-partial-transpose (bound) entangled states that negativity
   misses. Neither criterion dominates the other. These fix the method-explicit facts; the choice of
   the QuantumEntangledQ default over them is a separate, deliberate decision. *)

(* A PSD two-qubit entangled state (negativity > 0, exact here) that CCNR fails to flag: a rank-deficient
   real state with 0.5% maximally-mixed noise sits below the realignment threshold while still entangled. *)
ccnrMissRho = With[{
    r0 = {{46/1000, 33/1000, -166/1000, -22/1000}, {33/1000, 30/1000, -117/1000, -59/1000},
          {-166/1000, -117/1000, 604/1000, 56/1000}, {-22/1000, -59/1000, 56/1000, 320/1000}}
}, N[(1 - 1/200) r0 + (1/200) IdentityMatrix[4] / 4]];

VerificationTest[
    {
        Min @ Re @ Eigenvalues @ ccnrMissRho >= 0,                                          (* genuinely PSD *)
        QuantumEntanglementMonotone[QuantumState[ccnrMissRho, {2, 2}], "Negativity"] > 0,   (* entangled *)
        QuantumEntanglementMonotone[QuantumState[ccnrMissRho, {2, 2}], "Realignment"] <= 0  (* CCNR blind *)
    },
    {True, True, True},
    TestID -> "Criterion-CCNR-MissesWhatPPTCatches-2qubit"
]

VerificationTest[
    {
        QuantumEntangledQ[QuantumState[ccnrMissRho, {2, 2}], Automatic, "Negativity"],
        QuantumEntangledQ[QuantumState[ccnrMissRho, {2, 2}], Automatic, "Realignment"]
    },
    {True, False},
    TestID -> "Criterion-EntangledQ-NegativityCatches-RealignmentMisses"
]

(* The reverse blind spot: the Horodecki 3x3 bound entangled state (P. Horodecki, Phys. Lett. A 232, 333
   (1997)) is PPT, so negativity is 0 and misses it, while realignment/CCNR detects it. *)
horodeckiRho = With[{a = 1/2}, (1/(8 a + 1)) {
    {a, 0, 0, 0, a, 0, 0, 0, a}, {0, a, 0, 0, 0, 0, 0, 0, 0}, {0, 0, a, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, a, 0, 0, 0, 0, 0}, {a, 0, 0, 0, a, 0, 0, 0, a}, {0, 0, 0, 0, 0, a, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, (1 + a)/2, 0, Sqrt[1 - a^2]/2}, {0, 0, 0, 0, 0, 0, 0, a, 0},
    {a, 0, 0, 0, a, 0, Sqrt[1 - a^2]/2, 0, (1 + a)/2}
}];

VerificationTest[
    {
        Chop[QuantumEntanglementMonotone[QuantumState[N @ horodeckiRho, {3, 3}], "Negativity"]],
        QuantumEntanglementMonotone[QuantumState[N @ horodeckiRho, {3, 3}], "Realignment"] > 0
    },
    {0, True},
    TestID -> "Criterion-Horodecki-PPTMisses-CCNRCatches"
]

VerificationTest[
    {
        QuantumEntangledQ[QuantumState[N @ horodeckiRho, {3, 3}], Automatic, "Negativity"],
        QuantumEntangledQ[QuantumState[N @ horodeckiRho, {3, 3}], Automatic, "Realignment"]
    },
    {False, True},
    TestID -> "Criterion-EntangledQ-Horodecki-RealignmentCatches-NegativityMisses"
]

(* The dimension-aware default certifies both examples: negativity for the 2-qubit state (exact there),
   the negativity-or-realignment union for the Horodecki 3x3 bound entangled state. *)
VerificationTest[
    {
        QuantumEntangledQ[QuantumState[ccnrMissRho, {2, 2}]],
        QuantumEntangledQ[QuantumState[N @ horodeckiRho, {3, 3}]]
    },
    {True, True},
    TestID -> "EntangledQ-DimensionAwareDefault-CertifiesBoth"
]

(* a separable state stays False under the default in both dimensions *)
VerificationTest[
    {
        QuantumEntangledQ[sep],
        QuantumEntangledQ[QuantumState[N[IdentityMatrix[9] / 9], {3, 3}]]
    },
    {False, False},
    TestID -> "EntangledQ-DimensionAwareDefault-SeparableFalse"
]


(* ========== EntanglementEntropy is guarded to pure inputs ========== *)

(* a genuinely mixed state is not an entanglement-entropy case: it returns Indeterminate with a message,
   since the reduced von Neumann entropy there counts classical ignorance, not entanglement *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[N[IdentityMatrix[4] / 4], {2, 2}], "EntanglementEntropy"],
    Indeterminate,
    {QuantumEntanglementMonotone::mixedentropy},
    TestID -> "EntanglementEntropy-MixedGuarded"
]

(* a pure state supplied as a density matrix is still computed (the reduced-state branch), no message *)
VerificationTest[
    QuantityMagnitude @ QuantumEntanglementMonotone[bellDM, "EntanglementEntropy"],
    1,
    TestID -> "EntanglementEntropy-PureDensityMatrixStillComputes"
]

(* RenyiEntropy shares the pure-only semantics, so it is guarded the same way *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[N[IdentityMatrix[4] / 4], {2, 2}], "RenyiEntropy"],
    Indeterminate,
    {QuantumEntanglementMonotone::mixedentropy},
    TestID -> "RenyiEntropy-MixedGuarded"
]

(* a pure state still computes: for equal Schmidt weights every Renyi order equals the von Neumann value,
   so the alpha = 1/2 Renyi entanglement entropy of the Bell state is 1 bit, from both the vector and the
   density-matrix form *)
VerificationTest[
    Chop[{QuantumEntanglementMonotone[bell, "RenyiEntropy"], QuantumEntanglementMonotone[bellDM, "RenyiEntropy"]} - 1],
    {0, 0},
    TestID -> "RenyiEntropy-PureStillComputes"
]


(* ========== Non-physical (non-PSD) input warning ========== *)

(* a Hermitian 2-qubit matrix with a negative eigenvalue is not a physical state; a monotone still returns
   a number but warns, since the value is meaningless there *)
VerificationTest[
    QuantumEntanglementMonotone[
        QuantumState[{{0.5, 0, 0, 0.7}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0.7, 0, 0, 0.5}}, {2, 2}],
        "Negativity"
    ],
    _ ? NumericQ,
    {QuantumEntanglementMonotone::notphysical},
    SameTest -> MatchQ,
    TestID -> "Monotone-NonPSD-Warns"
]

(* a physical mixed state does not warn *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[N @ werner[4/5], {2, 2}], "Negativity"],
    _ ? NumericQ,
    {},
    SameTest -> MatchQ,
    TestID -> "Monotone-PSD-NoWarn"
]

(* A monotone computes from its input divided by the trace (a state vector divided by its norm), so the
   warning is judged on that matrix: a non-physical matrix warns and a physical one is silent at any scale at
   which the kernel can decide the sign of the trace, and at which the machine entries and a state vector's
   squared norm stay in the normal range. A matrix that is not Hermitian is not a state even with a
   positive spectrum, and neither is one whose trace is 0, negative, or not real. Each warning test makes
   one call, so the declared message is pinned to that call. *)

nonHermitian = KroneckerProduct[#, #] &[{{1/2, 1}, {0, 1/2}}];   (* trace 1, the eigenvalue 1/4 four times *)
pureDM = Outer[Times, #, Conjugate[#]] &[Normalize[{1, 2, 3, 4 I}]];   (* a pure state as a density matrix *)

(* the non-physical matrix above, scaled by 10^-9, still warns and gives its unscaled negativity *)
VerificationTest[
    Chop[
        QuantumEntanglementMonotone[
            QuantumState[1.*^-9 {{0.5, 0, 0, 0.7}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0.7, 0, 0, 0.5}}, {2, 2}],
            "Negativity"
        ] - 7/10
    ],
    0,
    {QuantumEntanglementMonotone::notphysical},
    TestID -> "Monotone-NonPSD-Scaled-Warns"
]

(* a matrix that is not Hermitian warns though its eigenvalues are positive *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[nonHermitian, {2, 2}], "Negativity"],
    _ ? NumericQ,
    {QuantumEntanglementMonotone::notphysical},
    SameTest -> MatchQ,
    TestID -> "Monotone-NonHermitian-Warns"
]

(* in machine numbers and scaled by 10^-9 it still warns: every entry of the scaled matrix lies below the
   10^-8 tolerance of the Hermiticity test, which therefore runs on the matrix divided by its trace *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[1.*^-9 N[nonHermitian], {2, 2}], "Negativity"],
    _ ? NumericQ,
    {QuantumEntanglementMonotone::notphysical},
    SameTest -> MatchQ,
    TestID -> "Monotone-NonHermitian-Scaled-Warns"
]

(* a zero trace cannot be divided out: the zero vector warns *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[{0, 0, 0, 0}, {2, 2}], "Negativity"],
    _ ? NumericQ,
    {QuantumEntanglementMonotone::notphysical},
    SameTest -> MatchQ,
    TestID -> "Monotone-ZeroTrace-Warns"
]

(* a negative trace warns at any scale: -10^-9 times the Werner state has no eigenvalue below -10^-8, and the
   division by its trace hands the monotone the Werner state itself, so it returns the Werner negativity *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[-10^-9 werner[4/5], {2, 2}], "Negativity"],
    7/20,
    {QuantumEntanglementMonotone::notphysical},
    TestID -> "Monotone-NegativeTrace-Warns"
]

(* a physical state is silent at any scale: the Werner state times 10^-9, and a pure state times 10^9, whose
   zero eigenvalues come out with round-off proportional to the scale, beyond an absolute -10^-8 floor *)
VerificationTest[
    {
        QuantumEntanglementMonotone[QuantumState[1.*^-9 N @ werner[4/5], {2, 2}], "Negativity"],
        QuantumEntanglementMonotone[QuantumState[1.*^9 N[pureDM], {2, 2}], "Negativity"]
    },
    {_ ? NumericQ, _ ? NumericQ},
    {},
    SameTest -> MatchQ,
    TestID -> "Monotone-PSD-AnyScale-NoWarn"
]

(* an algebraic trace whose sign Positive cannot decide is reduced with RootReduce, and the state is physical:
   the density matrix of the ket (1 + I Sqrt[2]) |00> + |11> has the trace (1 + I Sqrt[2]) (1 - I Sqrt[2])
   + 1, which is 4, and gives the negativity of the normalized ket without a warning *)
VerificationTest[
    Chop[
        N @ QuantumEntanglementMonotone[
            QuantumState[QuantumState[{1 + I Sqrt[2], 0, 0, 1}, {2, 2}]["DensityMatrix"], {2, 2}],
            "Negativity"
        ] - Sqrt[3]/4
    ],
    0,
    {},
    TestID -> "Monotone-PSD-UndecidedTrace-NoWarn"
]

(* an asymmetry of 10^-17 opposite a structural zero of a SparseArray is round-off, not a defect *)
VerificationTest[
    Chop[
        QuantumEntanglementMonotone[
            QuantumState[N[werner[4/5]] + SparseArray[{{2, 3} -> 1.*^-17}, {4, 4}], {2, 2}],
            "Negativity"
        ] - 7/20
    ],
    0,
    {},
    TestID -> "Monotone-PSD-SparseRoundoff-NoWarn"
]

(* a ket is judged in the computational basis: {1, 0, -1, 0}, in a basis whose two elements for the first
   qubit are both |0>, is the zero vector there, and warns *)
VerificationTest[
    QuantumEntanglementMonotone[
        QuantumState[{1, 0, -1, 0}, QuantumTensorProduct[QuantumBasis[{{1, 0}, {1, 0}}], QuantumBasis[2]]],
        "Negativity"
    ],
    _ ? NumericQ,
    {QuantumEntanglementMonotone::notphysical},
    SameTest -> MatchQ,
    TestID -> "Monotone-ZeroKet-DependentBasis-Warns"
]

(* Exact numbers on which Equal or Positive fails. zR and zR3 are 0 and trig2 is 2: Positive cannot decide
   the sign of zR or trig2, and Equal cannot decide zR3 == 0 until RootReduce reduces it. QuantumState stores
   a ket with several entries zR in a SparseArray whose default element is zR, on which Norm returns a wrong
   value or crashes the kernel. The sign of a trace scaled by trig2 stays undecided, and the matrix is then
   judged as given. tA is positive, of order 10^-200, so 1/2 + I tA is not real. These pin the helper's
   verdict, which decides the warning, without the messages the monotones' own arithmetic on such entries
   gives. *)

zR = Sqrt[3 + 2 Sqrt[2]] - 1 - Sqrt[2];
zR3 = Expand[(Sqrt[3 + 2 Sqrt[2]] - Sqrt[2])^2] - 1;
tA = Sqrt[2 + 10^-200] - Sqrt[2];
trig2 = 2 (Cos[1] + I Sin[1]) (Cos[1] - I Sin[1]);
notPhysicalQ = Wolfram`QuantumFramework`PackageScope`numericStateNotPSDQ;

(* a ket whose only nonzero entry is Equal-undecidable zR3 is the zero vector, and warns *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[{zR3, 0, 0, 0}, {2, 2}], "Negativity"],
    _ ? NumericQ,
    {QuantumEntanglementMonotone::notphysical},
    SameTest -> MatchQ,
    TestID -> "Monotone-ExactZeroKet-Warns"
]

(* an exact rational is compared as it is: a ket whose squared norm is 10^-400, below the machine range, is
   physical, and neither warns nor raises General::munfl *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[{10^-200, 0, 0, 0}, {2, 2}], "Negativity"],
    _ ? NumericQ,
    {},
    SameTest -> MatchQ,
    TestID -> "Monotone-TinyExactKet-NoWarn"
]

VerificationTest[
    {
        notPhysicalQ[QuantumState[{zR3, 0}]],
        notPhysicalQ[QuantumState[{zR, zR, zR, 1}, {2, 2}]],
        notPhysicalQ[QuantumState[{{1/2 + I tA, 0}, {0, 1/2}}]],
        notPhysicalQ[QuantumState[trig2 werner[4/5], {2, 2}]],
        notPhysicalQ[QuantumState[trig2 {{1/2, 0, 0, 7/10}, {0, 0, 0, 0}, {0, 0, 0, 0}, {7/10, 0, 0, 1/2}}, {2, 2}]]
    },
    {True, False, True, False, True},
    TestID -> "Helper-UndecidableExactEntries"
]

EndTestSection[]
