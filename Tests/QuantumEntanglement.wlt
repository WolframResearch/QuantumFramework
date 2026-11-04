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
    {QuantumEntanglementMonotone::realignment, QuantumEntanglementMonotone::realignment},
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
    {QuantumEntanglementMonotone::realignment},
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
    {QuantumEntanglementMonotone::realignment},
    TestID -> "Criterion-CCNR-MissesWhatPPTCatches-2qubit"
]

VerificationTest[
    {
        QuantumEntangledQ[QuantumState[ccnrMissRho, {2, 2}], Automatic, "Negativity"],
        QuantumEntangledQ[QuantumState[ccnrMissRho, {2, 2}], Automatic, "Realignment"]
    },
    {True, Indeterminate},
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
    {QuantumEntanglementMonotone::realignment},
    TestID -> "Criterion-Horodecki-PPTMisses-CCNRCatches"
]

VerificationTest[
    {
        QuantumEntangledQ[QuantumState[N @ horodeckiRho, {3, 3}], Automatic, "Negativity"],
        QuantumEntangledQ[QuantumState[N @ horodeckiRho, {3, 3}], Automatic, "Realignment"]
    },
    {Indeterminate, True},
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
    Chop[QuantityMagnitude[{QuantumEntanglementMonotone[bell, "RenyiEntropy"], QuantumEntanglementMonotone[bellDM, "RenyiEntropy"]}] - 1],
    {0, 0},
    TestID -> "RenyiEntropy-PureStillComputes"
]


(* ========== Non-physical input refused ========== *)

(* a Hermitian 2-qubit matrix with a negative eigenvalue is not a physical state, and a monotone refuses it *)
VerificationTest[
    QuantumEntanglementMonotone[
        QuantumState[{{0.5, 0, 0, 0.7}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0.7, 0, 0, 0.5}}, {2, 2}],
        "Negativity"
    ],
    _Failure,
    {QuantumEntanglementMonotone::notphysical},
    SameTest -> MatchQ,
    TestID -> "Monotone-NonPSD-Refused"
]

(* a physical mixed state is read without a message *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[N @ werner[4/5], {2, 2}], "Negativity"],
    _ ? NumericQ,
    {},
    SameTest -> MatchQ,
    TestID -> "Monotone-PSD-NoWarn"
]

(* A monotone computes from its input divided by the trace (a state vector divided by its norm), and the
   matrix is judged relative to its own scale: a non-physical matrix is refused and a physical one is read
   without a message at any scale, a machine one for a machine matrix and a rational one for a matrix of exact
   rational or complex rational entries. A matrix that is not Hermitian is not a state even with a positive
   spectrum, and neither is one whose trace is 0, negative, or not real. Each refusal test makes one call, so
   the declared message is pinned to that call. *)

nonHermitian = KroneckerProduct[#, #] &[{{1/2, 1}, {0, 1/2}}];   (* trace 1, the eigenvalue 1/4 four times *)
pureDM = Outer[Times, #, Conjugate[#]] &[Normalize[{1, 2, 3, 4 I}]];   (* a pure state as a density matrix *)

(* the non-physical matrix above, scaled by 10^-9, is refused as well *)
VerificationTest[
    QuantumEntanglementMonotone[
        QuantumState[1.*^-9 {{0.5, 0, 0, 0.7}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0.7, 0, 0, 0.5}}, {2, 2}],
        "Negativity"
    ],
    _Failure,
    {QuantumEntanglementMonotone::notphysical},
    SameTest -> MatchQ,
    TestID -> "Monotone-NonPSD-Scaled-Refused"
]

(* a matrix that is not Hermitian is refused though its eigenvalues are positive *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[nonHermitian, {2, 2}], "Negativity"],
    _Failure,
    {QuantumEntanglementMonotone::notphysical},
    SameTest -> MatchQ,
    TestID -> "Monotone-NonHermitian-Refused"
]

(* in machine numbers and scaled by 10^-9 it is refused as well: every entry of the scaled matrix lies below
   10^-8, and the Hermiticity test reads it relative to its largest entry *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[1.*^-9 N[nonHermitian], {2, 2}], "Negativity"],
    _Failure,
    {QuantumEntanglementMonotone::notphysical},
    SameTest -> MatchQ,
    TestID -> "Monotone-NonHermitian-Scaled-Refused"
]

(* a zero trace cannot be divided out: the zero vector is refused *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[{0, 0, 0, 0}, {2, 2}], "Negativity"],
    _Failure,
    {QuantumEntanglementMonotone::zero},
    SameTest -> MatchQ,
    TestID -> "Monotone-ZeroKet-Refused"
]

(* a negative trace is refused at any scale: -10^-9 times the Werner state has no eigenvalue below -10^-8, but
   every one of its nonzero eigenvalues is negative *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[-10^-9 werner[4/5], {2, 2}], "Negativity"],
    _Failure,
    {QuantumEntanglementMonotone::notphysical},
    SameTest -> MatchQ,
    TestID -> "Monotone-NegativeTrace-Refused"
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

(* 10^-300 times the Werner state lies in the normal range of machine numbers, and its trace is judged
   against its scale without forming a number below that range: the monotone reads it without a message and
   raises no General::munfl *)
VerificationTest[
    Chop[QuantumEntanglementMonotone[QuantumState[1.*^-300 N @ werner[4/5], {2, 2}], "Negativity"] - 7/20],
    0,
    {},
    TestID -> "Monotone-PSD-TinyMachineScale-NoMessage"
]

(* a trace whose sign Positive cannot decide is not taken for a negative or non-real one: the density matrix
   of the ket (1 + I Sqrt[2]) |00> + |11> has the trace (1 + I Sqrt[2]) (1 - I Sqrt[2]) + 1, which is 4,
   and gives the negativity of the normalized ket without a message *)
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

(* {1, 0, -1, 0}, in a basis whose two elements for the first qubit are both |0>, is refused: the elements of
   its basis are not independent *)
VerificationTest[
    QuantumEntanglementMonotone[
        QuantumState[{1, 0, -1, 0}, QuantumTensorProduct[QuantumBasis[{{1, 0}, {1, 0}}], QuantumBasis[2]]],
        "Negativity"
    ],
    _Failure,
    {QuantumEntanglementMonotone::dependent},
    SameTest -> MatchQ,
    TestID -> "Monotone-DependentBasis-Refused"
]

(* the entropies, as every monotone, refuse a single qubit, which has no split, before anything is computed: a
   qubit with the eigenvalue -1/5 gets that one message, not the refusal of a non-physical state *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[{{0.5, 0.7}, {0.7, 0.5}}], "EntanglementEntropy"],
    _Failure,
    {QuantumEntanglementMonotone::nosplit},
    SameTest -> MatchQ,
    TestID -> "Monotone-EntanglementEntropy-SingleQudit-Refused"
]

VerificationTest[
    QuantumEntanglementMonotone[QuantumState[{{0.5, 0.7}, {0.7, 0.5}}], "RenyiEntropy"],
    _Failure,
    {QuantumEntanglementMonotone::nosplit},
    SameTest -> MatchQ,
    TestID -> "Monotone-RenyiEntropy-SingleQudit-Refused"
]

(* Exact numbers on which Equal or Positive fails. zR and zR3 are 0 and trig2 is 2: Positive cannot decide
   the sign of zR or trig2, and Equal cannot decide zR3 == 0 until RootReduce reduces it. QuantumState stores
   the ket {zR, zR, zR, 1} in a SparseArray whose default element is zR, and the helper reads it as the
   nonzero vector it is. The sign of a trace scaled by trig2 stays undecided, and the matrix is then judged
   as given. tA is positive, of order 10^-200, so the trace 1 + I tA is not real. These pin the verdict of
   QuantumDistance's helper, which decides its notphysical warning, without the messages the monotones' own
   arithmetic on such entries gives. *)

zR = Sqrt[3 + 2 Sqrt[2]] - 1 - Sqrt[2];
zR3 = Expand[(Sqrt[3 + 2 Sqrt[2]] - Sqrt[2])^2] - 1;
tA = Sqrt[2 + 10^-200] - Sqrt[2];
trig2 = 2 (Cos[1] + I Sin[1]) (Cos[1] - I Sin[1]);
notPhysicalQ = Wolfram`QuantumFramework`PackageScope`numericStateNotPSDQ;

(* a ket whose only nonzero entry is Equal-undecidable zR3 is the zero vector, and is refused *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[{zR3, 0, 0, 0}, {2, 2}], "Negativity"],
    _Failure,
    {QuantumEntanglementMonotone::zero},
    SameTest -> MatchQ,
    TestID -> "Monotone-ExactZeroKet-Refused"
]

(* an exact rational is compared as it is: a ket whose squared norm is 10^-400, below the machine range, is
   physical, is read without a message and raises no General::munfl *)
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

(* a matrix with a symbolic entry is judged through its trace alone, and a machine trace against the norm of
   its numeric entries, which scales with the matrix: the state is physical at 10^-9 as at 1, and a trace of
   10^-12 against numeric entries of order 1 counts as 0, so that matrix is not physical *)
VerificationTest[
    With[{m = {{0.5, 0, 0, x}, {0, 0, 0, 0}, {0, 0, 0, 0}, {Conjugate[x], 0, 0, 0.5}}},
        {
            notPhysicalQ[QuantumState[m, {2, 2}]],
            notPhysicalQ[QuantumState[1.*^-9 m, {2, 2}]],
            notPhysicalQ[QuantumState[ReplacePart[m, {4, 4} -> -0.5 + 1.*^-12], {2, 2}]]
        }
    ],
    {False, False, True},
    TestID -> "Helper-SymbolicEntry-AnyScale"
]

(* an exact trace is compared with 0 and 1 without converting it to a machine number: the trace
   Sqrt[2] 10^-400, and the squared norm of the ket {2^(1/4) 10^-200, 0}, lie below the machine range, and
   raise no General::munfl *)
VerificationTest[
    {
        notPhysicalQ[QuantumState[Sqrt[2] 10^-400 {{1/2, 0}, {0, 1/2}}]],
        notPhysicalQ[QuantumState[{2^(1/4) 10^-200, 0}]]
    },
    {False, False},
    {},
    TestID -> "Helper-TinyExactTrace-NoMessage"
]

(* arbitrary-precision entries keep their precision in the scale their trace is judged against: in 30-digit
   numbers and multiplied by 10^-400, below the range of machine numbers, the matrix with the eigenvalue
   -1/5 is not physical and the Werner state is, with no General::munfl *)
VerificationTest[
    {
        notPhysicalQ[QuantumState[10^-400 N[{{1/2, 0, 0, 7/10}, {0, 0, 0, 0}, {0, 0, 0, 0}, {7/10, 0, 0, 1/2}}, 30], {2, 2}]],
        notPhysicalQ[QuantumState[10^-400 N[werner[4/5], 30], {2, 2}]]
    },
    {True, False},
    {},
    TestID -> "Helper-ArbitraryPrecisionTinyScale"
]


(* ========== Concurrence: amplitude damping of an entangled pair, symbolic in the rates ========== *)

(* Sin[th] |00> + Cos[th] |11> under independent amplitude damping of each qubit, rates gA and gB, with
   QF's "-" = |1><0|, so |0> decays to |1>. With pX = Exp[-gX tau] the pair stays the X state
       rho00,00 = Sin^2 pA pB, rho01,01 = Sin^2 pA (1 - pB), rho10,10 = Sin^2 (1 - pA) pB,
       rho11,11 = Cos^2 + Sin^2 (1 - pA) (1 - pB), rho00,11 = Sin Cos Sqrt[pA pB],
   whose concurrence 2 Max[0, |rho00,11| - Sqrt[rho01,01 rho10,10]] is (Yu and Eberly, PRL 93, 140404 (2004))
       C = 2 Sqrt[pA pB] |Sin[th]| Max[0, |Cos[th]| - |Sin[th]| Sqrt[(1 - pA) (1 - pB)]].
   QuantumEvolve has no jump-only form (a None Hamiltonian returns unevaluated), so the Hamiltonian is the
   zero operator. *)
dampedPairArgs = Sequence[
    QuantumOperator[0 IdentityMatrix[4]],
    {QuantumOperator["-", {1}], QuantumOperator["-", {2}]} -> {\[Gamma]A, \[Gamma]B},
    QuantumState[{Sin[\[Theta]], 0, 0, Cos[\[Theta]]}],
    \[Tau]
];
dampedPairAssumptions = \[Gamma]A >= 0 && \[Gamma]B >= 0 && \[Tau] >= 0 && 0 < \[Theta] < Pi / 2;
dampedPairClosedForm[th_, gA_, gB_, tau_] := With[{pA = Exp[-gA tau], pB = Exp[-gB tau]},
    2 Sqrt[pA pB] Abs[Sin[th]] Max[0, Abs[Cos[th]] - Abs[Sin[th]] Sqrt[(1 - pA) (1 - pB)]]
];
realWithinQ[tol_] := VectorQ[#1, NumericQ] && VectorQ[#1, Im[#] == 0 &] && Max[Abs[#1 - #2]] < tol &;

(* the evolved state, under real nonnegative rates, is the X state above: unit trace, Hermitian and
   positive, the coherence bounded by the populations it joins, |rho00,11|^2 <= rho00,00 rho11,11, and
   decaying at (gA + gB) / 2, the mean of the two damping rates *)
VerificationTest[
    With[{
        rho = Normal @ QuantumEvolve[dampedPairArgs, Assumptions -> dampedPairAssumptions]["DensityMatrix"],
        pA = Exp[-\[Gamma]A \[Tau]], pB = Exp[-\[Gamma]B \[Tau]], s = Sin[\[Theta]], c = Cos[\[Theta]]
    },
        FullSimplify[
            {
                rho - {
                    {s^2 pA pB, 0, 0, s c Sqrt[pA pB]},
                    {0, s^2 pA (1 - pB), 0, 0},
                    {0, 0, s^2 (1 - pA) pB, 0},
                    {s c Sqrt[pA pB], 0, 0, c^2 + s^2 (1 - pA) (1 - pB)}
                },
                Tr[rho],
                rho - ConjugateTranspose[rho],
                Abs[rho[[1, 4]]]^2 <= rho[[1, 1]] rho[[4, 4]] && rho[[2, 2]] >= 0 && rho[[3, 3]] >= 0
            },
            dampedPairAssumptions
        ]
    ],
    {ConstantArray[0, {4, 4}], 1, ConstantArray[0, {4, 4}], True},
    {},
    TestID -> "Concurrence-DampedPair-XState"
]

(* only qubit A decays: no population reaches |01>, the coherence falls as Sqrt[pA], and the concurrence
   reduces symbolically to Sin[2 th] Exp[-gA tau / 2], a decay with no death *)
VerificationTest[
    FullSimplify[
        conc[QuantumEvolve[dampedPairArgs, Assumptions -> dampedPairAssumptions]] /. \[Gamma]B -> 0,
        dampedPairAssumptions
    ],
    Sin[2 \[Theta]] Exp[-\[Gamma]A \[Tau] / 2],
    {},
    TestID -> "Concurrence-DampedPair-OneSidedDecay-Symbolic"
]

(* The concurrence built with no assumptions and evaluated afterwards, as in a plot over the parameters.
   When a rate is zero an eigenvalue of rho.rho-tilde is exactly zero; machine numbers substituted into the
   closed form leave it as round-off of either sign, and an unclipped square root turns it imaginary: Max
   then cannot compare, and the value stays unevaluated. Tolerance 10^-6: a round-off zero under a square
   root is of order Sqrt[eps]. *)
dampedPairConcurrence = conc @ QuantumEvolve[dampedPairArgs];
dampedPairAt[th_, gA_, gB_, tau_] := dampedPairConcurrence /. {\[Theta] -> th, \[Gamma]A -> gA, \[Gamma]B -> gB, \[Tau] -> tau};
dampedPairGrid = Tuples[{N[{Pi/6, Pi/4, Pi/3, 5 Pi/12}], {0.3, 1.3, 3.1}}];
dampedPairExactGrid = Tuples[{{Pi/6, Pi/3}, {3/10, 31/10}}];

(* one qubit decays, the other is untouched: C = |Sin[2 th]| Exp[-g tau / 2] on either side *)
VerificationTest[
    Join[dampedPairAt[#1, 1., 0., #2] & @@@ dampedPairGrid, dampedPairAt[#1, 0., 0.7, #2] & @@@ dampedPairGrid],
    Join[Abs[Sin[2 #1]] Exp[-#2 / 2] & @@@ dampedPairGrid, Abs[Sin[2 #1]] Exp[-0.7 #2 / 2] & @@@ dampedPairGrid],
    {},
    SameTest -> realWithinQ[10^-6],
    TestID -> "Concurrence-DampedPair-OneSidedDecay-Machine"
]

(* no decay on either side: the evolution is the identity and the pure state keeps C = |Sin[2 th]| *)
VerificationTest[
    dampedPairAt[#1, 0., 0., #2] & @@@ dampedPairGrid,
    Abs[Sin[2 #1]] & @@@ dampedPairGrid,
    {},
    SameTest -> realWithinQ[10^-6],
    TestID -> "Concurrence-DampedPair-NoDecay-Machine"
]

(* equal rates: with more weight on the decaying |00>, |Sin[th]| > |Cos[th]|, the populations of |01> and
   |10> outgrow the coherence and the entanglement dies at the finite time 1 - Exp[-tau] = |Cot[th]|, after
   which C is exactly 0; with more weight on |11> it only decays. The times are 0.9, 1 and 1.1 death times
   of th = 5 Pi/12, and the same times for th = Pi/12, which never dies. *)
dampedPairDeathTimes = -{0.9, 1, 1.1} Log[1 - Cot[5 Pi/12.]];
VerificationTest[
    dampedPairAt[#1, 1., 1., #2] & @@@ Tuples[{N[{Pi/12, 5 Pi/12}], dampedPairDeathTimes}],
    dampedPairClosedForm[#1, 1., 1., #2] & @@@ Tuples[{N[{Pi/12, 5 Pi/12}], dampedPairDeathTimes}],
    {},
    SameTest -> realWithinQ[10^-6],
    TestID -> "Concurrence-DampedPair-SuddenDeath-Machine"
]

(* across the death time of th = 5 Pi/12 the concurrence is positive before, zero to round-off at, and
   exactly 0 after *)
VerificationTest[
    With[{c = dampedPairAt[5 Pi/12., 1., 1., #] & /@ dampedPairDeathTimes},
        {c[[1]] > 0, Abs[c[[2]]] < 10^-12, c[[3]] == 0}
    ],
    {True, True, True},
    {},
    TestID -> "Concurrence-DampedPair-DeathTime-Machine"
]

(* at the exact death time of th = 5 Pi/12, Exp[-tau] = Sqrt[3] - 1, the coherence equals the geometric mean
   of the |01> and |10> populations, |rho00,11| = Sqrt[rho01,01 rho10,10] = (Sqrt[3] - 1) / 4, and the
   concurrence is exactly 0. The state is reduced before the concurrence is taken: the closed form built
   without assumptions holds this zero as nested radicals that N cannot certify (N::meprec). *)
VerificationTest[
    conc @ QuantumState[
        RootReduce @ FullSimplify[Normal[QuantumEvolve[dampedPairArgs]["DensityMatrix"]] /.
            {\[Theta] -> 5 Pi/12, \[Gamma]A -> 1, \[Gamma]B -> 1, \[Tau] -> -Log[Sqrt[3] - 1]}],
        {2, 2}
    ],
    0,
    {},
    TestID -> "Concurrence-DampedPair-DeathTime-Exact"
]

(* the symbolic result is the Yu-Eberly closed form at unequal nonzero rates, exact and read to 20 digits.
   The two-rate form, equal rates included, is compared at exact points rather than by FullSimplify: its
   Wootters values are nested radicals, Sqrt[a +- 2 E^(g tau) Cos[th] Sqrt[...]], that FullSimplify does
   not denest. *)
VerificationTest[
    N[dampedPairAt[#1, 1, 3/10, #2] & @@@ dampedPairExactGrid, 20],
    N[dampedPairClosedForm[#1, 1, 3/10, #2] & @@@ dampedPairExactGrid, 20],
    {},
    SameTest -> realWithinQ[10^-15],
    TestID -> "Concurrence-DampedPair-ExactClosedForm"
]

(* long times: a surviving pair, |Cos[th]| > |Sin[th]|, loses entanglement as the product of the two
   coherences, C ~ 2 Sin[th] (Cos[th] - Sin[th]) Exp[-(gA + gB) tau / 2]. The relative correction is led
   by the slower rate: Sqrt[(1 - pA) (1 - pB)] ~ 1 - (pA + pB) / 2 gives
   C / C_asymptotic - 1 ~ Sin[th] / (2 (Cos[th] - Sin[th])) Exp[-gB tau] for gB < gA. Exact at th = Pi/6,
   gA = 1, gB = 3/10, tau = 60, where the next order is below 10^-8 of the coefficient. *)
VerificationTest[
    N[(dampedPairAt[Pi/6, 1, 3/10, 60] / (2 Sin[Pi/6] (Cos[Pi/6] - Sin[Pi/6]) Exp[-39]) - 1) Exp[18], 20],
    N[Sin[Pi/6] / (2 (Cos[Pi/6] - Sin[Pi/6])), 20],
    {},
    SameTest -> (Abs[#1 - #2] < 10^-8 &),
    TestID -> "Concurrence-DampedPair-LongTime"
]

(* ========== Monotones of an approximate state known to a given accuracy ========== *)

(* A numerical solution of the master equation is known only to the solver's error, about 10^-8 here, far above the
   round-off of its machine entries. At a zero rate the exact state has a zero eigenvalue (the |00>,|11> block is
   rank one when only qubit A decays), which that error pushes to either sign. Tolerance -> e accepts a state with no
   eigenvalue below -e (Weyl's inequality: a state within e of a density matrix has none) and reads the monotone on
   the nearest density matrix. The method, its references and the error each monotone keeps are in the comments of
   toleratedState and stateToleranceCut in QuantumEntanglement.m. *)
numericDampedPair[th_, gA_, gB_] := QuantumEvolve[
    QuantumOperator[0 IdentityMatrix[4]],
    {QuantumOperator["-", {1}], QuantumOperator["-", {2}]} -> {gA, gB},
    QuantumState[{Sin[th], 0, 0, Cos[th]}],
    {\[FormalT], 0, 4}
];
numericDampedPairAt[th_, gA_, gB_, tau_] := numericDampedPair[th, gA, gB] /. \[FormalT] -> tau;
numericZeroRateGrid = Tuples[{N[{Pi/8, Pi/3}], {0.3, 1.3, 3.}}];
(* the X state's partial transpose has the block {{rho01,01, rho00,11}, {rho00,11, rho10,10}}, whose smaller
   eigenvalue is the only one that can be negative *)
xStateNegativity[th_, gA_, gB_, tau_] := With[{pA = Exp[-gA tau], pB = Exp[-gB tau], s = Sin[th], c = Cos[th]},
    With[{a = s^2 pA (1 - pB), b = s^2 (1 - pA) pB, z = s c Sqrt[pA pB]}, Max[0, Sqrt[((a - b) / 2)^2 + z^2] - (a + b) / 2]]
];

(* the solver leaves the zero eigenvalue at -2.4 10^-9 at th = Pi/3, gA = 1, gB = 0, tau = 1.3: the default, the
   round-off of the entries, refuses it, an accuracy below that eigenvalue refuses it, one above accepts it *)
VerificationTest[
    With[{st = numericDampedPairAt[Pi/3., 1., 0., 1.3]},
        {
            FailureQ[QuantumEntanglementMonotone[st, "Concurrence"]],
            FailureQ[QuantumEntanglementMonotone[st, "Concurrence", Tolerance -> 10^-10]],
            NumericQ[QuantumEntanglementMonotone[st, "Concurrence", Tolerance -> 10^-8]]
        }
    ],
    {True, True, True},
    {QuantumEntanglementMonotone::notphysical, QuantumEntanglementMonotone::notphysical},
    TestID -> "Tolerance-WeylCut"
]

(* only qubit A decays: the concurrence of the numerical state is Sin[2 th] Exp[-gA tau / 2] and its negativity the
   X-state value, to within 10^-6. The scale is the solver's: its entries are good to about 10^-8, and the correction,
   of the size of the most negative eigenvalue, at most about 5 10^-8 on this grid, is far below the 10^-6 tested. This
   test checks the physics of the accepted state; the
   projection itself is checked on the Bell-diagonal state below. *)
VerificationTest[
    With[{e = 10^-6},
        {
            Max[Abs[(QuantumEntanglementMonotone[numericDampedPairAt[#1, 1., 0., #2], "Concurrence", Tolerance -> e] - Sin[2 #1] Exp[-#2 / 2]) & @@@ numericZeroRateGrid]] < 10^-6,
            Max[Abs[(QuantumEntanglementMonotone[numericDampedPairAt[#1, 1., 0., #2], "Negativity", Tolerance -> e] - xStateNegativity[#1, 1., 0., #2]) & @@@ numericZeroRateGrid]] < 10^-6
        }
    ],
    {True, True},
    {},
    TestID -> "Tolerance-NumericOneSidedDecay"
]

(* a state with no negative eigenvalue is not moved: at two nonzero rates the smallest eigenvalue is far above the
   solver's error, and the value with an accuracy given is the value without one *)
VerificationTest[
    With[{st = numericDampedPairAt[Pi/3., 1., 0.3, 1.3]},
        QuantumEntanglementMonotone[st, #, Tolerance -> 10^-6] - QuantumEntanglementMonotone[st, #] & /@ {"Concurrence", "Negativity"}
    ],
    {0., 0.},
    {},
    TestID -> "Tolerance-ValidStateUnchanged"
]

(* the correction is the Hilbert-Schmidt-nearest density matrix (Smolin, Gambetta and Smith, PRL 108, 070502
   (2012)): on the Bell-diagonal state with weights
   {7/10, 2/10, 1/10 + d, -d}, the negative weight is set to 0 and the others lowered by d / 3, so the concurrence
   2 lambda_max - 1 is 2/5 - 2 d / 3 and the negativity lambda_max - 1/2 is 1/5 - d / 3; clipping and rescaling would
   give 2/5 - 7 d / 5 *)
VerificationTest[
    With[{d = 10^-7, bell = Normal @ QuantumState[#]["DensityMatrix"] & /@ {"PhiPlus", "PhiMinus", "PsiPlus", "PsiMinus"}},
        With[{st = QuantumState[N[{7/10, 2/10, 1/10 + d, -d} . bell], {2, 2}]},
            {
                QuantumEntanglementMonotone[st, "Concurrence", Tolerance -> 10^-6] - (2/5 - 2 d / 3),
                QuantumEntanglementMonotone[st, "Negativity", Tolerance -> 10^-6] - (1/5 - d / 3)
            }
        ]
    ],
    {0, 0},
    {},
    SameTest -> (Max[Abs[#1]] < 10^-12 &),
    TestID -> "Tolerance-NearestStateProjection"
]

(* a matrix with an eigenvalue far below the accuracy is not a state and stays refused; an exact matrix is exact,
   so an accuracy does not admit its negative eigenvalue; an accuracy that is not a nonnegative number is refused *)
VerificationTest[
    {
        FailureQ @ QuantumEntanglementMonotone[QuantumState[N[{{1/2, 0, 0, 1/2}, {0, -1/100, 0, 0}, {0, 0, 1/100, 0}, {1/2, 0, 0, 1/2}}], {2, 2}], Tolerance -> 10^-6],
        FailureQ @ QuantumEntanglementMonotone[QuantumState[{{1/2, 0, 0, 1/2}, {0, -10^-9, 0, 0}, {0, 0, 10^-9, 0}, {1/2, 0, 0, 1/2}}, {2, 2}], Tolerance -> 10^-6],
        FailureQ @ QuantumEntanglementMonotone[QuantumState["PhiPlus"], Tolerance -> -1]
    },
    {True, True, True},
    {QuantumEntanglementMonotone::notphysical, QuantumEntanglementMonotone::notphysical, QuantumEntanglementMonotone::tolerance},
    TestID -> "Tolerance-Refusals"
]

(* the accuracy covers the anti-Hermitian part too: the Bell state plus an antisymmetric error of operator norm 0.3 e
   is within e of a density matrix, refused at the round-off of its entries and accepted at Tolerance -> e, its
   negativity that of the Bell state to within the error *)
VerificationTest[
    With[{e = 10^-6, st = QuantumState[N[Normal @ QuantumState["PhiPlus"]["DensityMatrix"] + 3/10 10^-6 SparseArray[{{2, 3} -> 1, {3, 2} -> -1}, {4, 4}]], {2, 2}]},
        {
            FailureQ[QuantumEntanglementMonotone[st, "Negativity"]],
            Abs[QuantumEntanglementMonotone[st, "Negativity", Tolerance -> e] - 1/2] < e
        }
    ],
    {True, True},
    {QuantumEntanglementMonotone::notphysical},
    TestID -> "Tolerance-AntiHermitianError"
]

(* a rule that is not an option is refused with the arguments, as a split *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState["PhiPlus"], "Negativity", Foo -> 1],
    Failure["InvalidArguments", <|"MessageTemplate" :> QuantumEntanglementMonotone::args, "MessageParameters" -> {{"Negativity", Foo -> 1}}|>],
    {QuantumEntanglementMonotone::args},
    TestID -> "Tolerance-UnknownOption"
]


(* ========== Entanglement of an approximate state known to a given accuracy ========== *)

(* With equal rates and th = Pi/3 the pair Sin[th] |00> + Cos[th] |11> loses its entanglement at the finite time
   1 - Exp[-tau] = Cot[th] (Yu and Eberly): there the smallest eigenvalue of the partial transpose crosses 0, and the
   solver's error, about 10^-7 at that time, gives it either sign. With Tolerance -> 10^-6 the answer holds for every
   state within that accuracy: entangled well before, separable well after (2 x 2, a positive partial transpose is
   separable), and Indeterminate at the death time itself. *)
VerificationTest[
    With[{ev = numericDampedPair[Pi/3., 1., 1.], td = - Log[1 - Cot[Pi/3.]]},
        QuantumEntangledQ[ev /. \[FormalT] -> #, Tolerance -> 10^-6] & /@ {0.5, td - 0.01, td, td + 0.01, 1.5}
    ],
    {True, True, Indeterminate, False, False},
    {},
    TestID -> "Tolerance-EntangledQ-SuddenDeath"
]

(* The isotropic qutrit state p |Phi><Phi| + (1 - p) I / 9 at p = 1/2 has smallest partial-transpose eigenvalue -1/9
   and realignment value ||R||_1 - 1 = 2/3. An error e moves the first by up to 3 e and the second by up to 9 e: at
   e = 0.05 the partial transpose no longer decides, the realignment that Automatic reads beyond 2 x 3 still proves
   entanglement, and at e = 0.1 neither does. *)
VerificationTest[
    With[{iso = QuantumState[N[Outer[Times, Flatten[IdentityMatrix[3]], Flatten[IdentityMatrix[3]]] / 6 + IdentityMatrix[9] / 18], {3, 3}]},
        {
            QuantumEntangledQ[iso, Tolerance -> 10^-6],
            QuantumEntangledQ[iso, Tolerance -> 0.05],
            QuantumEntangledQ[iso, "Negativity", Tolerance -> 0.05],
            QuantumEntangledQ[iso, Tolerance -> 0.1]
        }
    ],
    {True, True, Indeterminate, Indeterminate},
    {},
    TestID -> "Tolerance-EntangledQ-IsotropicQutrits"
]

(* The isotropic qutrit state at p = 1/20 is at Hilbert-Schmidt distance (1/20) Sqrt[8/9] from I/9, inside the ball
   of radius 1 / Sqrt[72] in which every state is separable (Gurvits and Barnum): separable for e = 10^-6, whose
   9 x 9 neighborhood stays inside, and Indeterminate for e = 0.03, whose neighborhood reaches past it. *)
VerificationTest[
    With[{iso = QuantumState[N[Outer[Times, Flatten[IdentityMatrix[3]], Flatten[IdentityMatrix[3]]] / 60 + 19 IdentityMatrix[9] / 180], {3, 3}]},
        {QuantumEntangledQ[iso, Tolerance -> 10^-6], QuantumEntangledQ[iso, Tolerance -> 0.03]}
    ],
    {False, Indeterminate},
    {},
    TestID -> "Tolerance-EntangledQ-PurityBall"
]

(* Two qubits of GHZ hold (|00><00| + |11><11|) / 2, classically correlated and separable, but with a partial
   transpose of eigenvalue 0: adding e (|00><11| + |11><00|) gives it the eigenvalue -e, an entangled state within any
   accuracy. So without Tolerance it is separable, with any Tolerance Indeterminate, while the cut of one qubit from
   the other two is entangled either way. *)
VerificationTest[
    With[{g = QuantumState[N @ Normal @ QuantumState["GHZ"[3]]["DensityMatrix"], {2, 2, 2}]},
        {
            QuantumEntangledQ[g, {{1}, {2}}],
            QuantumEntangledQ[g, {{1}, {2}}, Tolerance -> 10^-6],
            QuantumEntangledQ[g, {{1}, {2, 3}}, Tolerance -> 10^-6]
        }
    ],
    {False, Indeterminate, True},
    {},
    TestID -> "Tolerance-EntangledQ-ClassicalCorrelations"
]

(* an accuracy that is not a nonnegative number, a method other than the partial-transpose tests, and an unknown
   option are refused; an exact state has no error to bound and is answered as without Tolerance *)
VerificationTest[
    {
        FailureQ @ QuantumEntangledQ[QuantumState["PhiPlus"], Tolerance -> -1],
        FailureQ @ QuantumEntangledQ[QuantumState[N @ Normal @ QuantumState["PhiPlus"]["DensityMatrix"], {2, 2}], "Concurrence", Tolerance -> 10^-6],
        FailureQ @ QuantumEntangledQ[QuantumState["PhiPlus"], Foo -> 1],
        QuantumEntangledQ[QuantumState["PhiPlus"], Tolerance -> 1/10]
    },
    {True, True, True, True},
    {QuantumEntangledQ::tolerance, QuantumEntangledQ::tolmethod, QuantumEntangledQ::args},
    TestID -> "Tolerance-EntangledQ-Refusals"
]


(* ========== Monotones of a symbolic state under assumptions on its parameters ========== *)

(* The Werner state p |Psi-><Psi-| + (1 - p) I / 4 is entangled exactly for p > 1/3, with concurrence (3 p - 1) / 2
   and negativity (3 p - 1) / 4 there and 0 below. Without assumptions the concurrence is a Max over square roots of
   conjugated p; the assumptions decide every comparison, so the value is the threshold law itself, and the same
   assumptions given through $Assumptions act as the option does. *)
VerificationTest[
    Block[{p},
        With[{w = QuantumState[p Outer[Times, {0, 1, -1, 0}, {0, 1, -1, 0}] / 2 + (1 - p) IdentityMatrix[4] / 4, {2, 2}]},
            Expand[{
                QuantumEntanglementMonotone[w, "Concurrence", Assumptions -> 1/3 < p <= 1],
                QuantumEntanglementMonotone[w, "Concurrence", Assumptions -> 0 <= p < 1/3],
                QuantumEntanglementMonotone[w, "Negativity", Assumptions -> 1/3 < p <= 1],
                Block[{$Assumptions = 1/3 < p <= 1}, QuantumEntanglementMonotone[w, "Concurrence"]]
            } - {(3 p - 1) / 2, 0, (3 p - 1) / 4, (3 p - 1) / 2}]
        ]
    ],
    {0, 0, 0, 0},
    {},
    TestID -> "Assumptions-WernerThreshold"
]

(* Cos[th] |00> + Sin[th] |11> with 0 < th < Pi/2 has Schmidt coefficients Cos[th] and Sin[th]: concurrence
   2 Cos[th] Sin[th] = Sin[2 th], negativity Cos[th] Sin[th], half of it as for every pure pair of qubits, and
   entanglement entropy the binary entropy of Cos[th]^2. Under the assumptions every conjugate and absolute value of
   th is resolved. *)
VerificationTest[
    Block[{th},
        With[{ps = QuantumState[{Cos[th], 0, 0, Sin[th]}, {2, 2}], asm = 0 < th < Pi/2},
            With[{v = {
                QuantumEntanglementMonotone[ps, "Concurrence", Assumptions -> asm],
                QuantumEntanglementMonotone[ps, "Negativity", Assumptions -> asm],
                QuantityMagnitude @ QuantumEntanglementMonotone[ps, "EntanglementEntropy", Assumptions -> asm]
            }},
                {
                    FreeQ[v, Conjugate | Abs],
                    Simplify[v - {Sin[2 th], Sin[2 th] / 2, - Cos[th]^2 Log2[Cos[th]^2] - Sin[th]^2 Log2[Sin[th]^2]}, asm]
                }
            ]
        ]
    ],
    {True, {0, 0, 0}},
    {},
    TestID -> "Assumptions-PureStateMonotones"
]

(* Only qubit A decays, with "-" = |1><0| on it: the pair Sin[th] |00> + Cos[th] |11> keeps its |01>, |10>
   populations at 0, its coherence falls as Exp[-g t / 2], and the concurrence is Sin[2 th] Exp[-g t / 2] with no
   sudden death. With the assumptions the monotone returns that closed form directly. *)
VerificationTest[
    Block[{g, t, th},
        QuantumEntanglementMonotone[
            QuantumEvolve[QuantumOperator[0 IdentityMatrix[4]], {QuantumOperator["-", {1}], QuantumOperator["-", {2}]} -> {g, 0}, QuantumState[{Sin[th], 0, 0, Cos[th]}], t],
            "Concurrence",
            Assumptions -> 0 < th < Pi/2 && g > 0 && t > 0
        ] - Sin[2 th] Exp[- g t / 2]
    ],
    0,
    {},
    TestID -> "Assumptions-OneSidedDamping"
]

(* a state of numbers has no parameters for the assumptions to act on: its value is the one without them *)
VerificationTest[
    Block[{p},
        With[{w = QuantumState[N[4/5 Outer[Times, {0, 1, -1, 0}, {0, 1, -1, 0}] / 2 + 1/5 IdentityMatrix[4] / 4], {2, 2}]},
            QuantumEntanglementMonotone[w, "Concurrence", Assumptions -> p > 0] === QuantumEntanglementMonotone[w, "Concurrence"]
        ]
    ],
    True,
    {},
    TestID -> "Assumptions-NumericStateUnchanged"
]

(* assumptions that no parameter value meets are refused: Simplify under them may return any value *)
VerificationTest[
    Block[{p},
        With[{w = QuantumState[p Outer[Times, {0, 1, -1, 0}, {0, 1, -1, 0}] / 2 + (1 - p) IdentityMatrix[4] / 4, {2, 2}]},
            FailureQ /@ {
                QuantumEntanglementMonotone[w, "Concurrence", Assumptions -> p > 1 && p < 0],
                QuantumEntanglementMonotone[w, "Concurrence", Assumptions -> False]
            }
        ]
    ],
    {True, True},
    {QuantumEntanglementMonotone::assumptions, QuantumEntanglementMonotone::assumptions},
    TestID -> "Assumptions-Contradictory"
]


(* ========== Negativity of an exact pure state ========== *)

(* The partial transpose of a pure pair with Schmidt coefficients s_i has eigenvalues s_i^2 and +- s_i s_j, so the
   negativity is ((Sum[s_i])^2 - 1) / 2 and the log-negativity 2 Log2[Sum[s_i]]. The maximally entangled pair of
   qutrits has s_i = 1 / Sqrt[3]: negativity 1 and log-negativity Log2[3], as a state vector and as a density matrix. *)
VerificationTest[
    With[{v = Flatten[IdentityMatrix[3]] / Sqrt[3]},
        Flatten @ Table[
            QuantumEntanglementMonotone[QuantumState[s, {3, 3}], m],
            {s, {v, Outer[Times, v, v]}}, {m, {"Negativity", "LogNegativity"}}
        ]
    ],
    {1, Log2[3], 1, Log2[3]},
    {},
    TestID -> "ExactPureNegativity-MaximallyEntangledQutrits"
]

(* a product of two qutrit states with nested-radical amplitudes has one Schmidt coefficient, 1, and negativity
   exactly 0; through the 9 x 9 partial transpose the exact trace norm does not finish in two minutes *)
VerificationTest[
    TimeConstrained[
        With[{v = Normalize @ Flatten @ KroneckerProduct[{1, 2^(1/3), I}, {Sqrt[2], 1 - I, 3^(1/4)}]},
            QuantumEntanglementMonotone[QuantumState[Outer[Times, v, Conjugate[v]], {3, 3}], "Negativity"]
        ],
        10
    ],
    0,
    {},
    TestID -> "ExactPureNegativity-NestedRadicalProduct"
]

(* an entangled exact pair of a qutrit and a qubit, part A the larger: the exact value agrees to 30 digits with
   ((Total of the absolute eigenvalues of the partial transpose) - 1) / 2 computed at 50 digits *)
VerificationTest[
    With[{v = Normalize[{1, 2, I, 3^(1/4), 0, 1 - I}]},
        With[{rho = Outer[Times, v, Conjugate[v]]},
            Abs[
                N[QuantumEntanglementMonotone[QuantumState[rho, {3, 2}], "Negativity"], 40] -
                (Total[Abs[Eigenvalues[N[Flatten[Transpose[ArrayReshape[rho, {3, 2, 3, 2}], {1, 4, 3, 2}], {{1, 2}, {3, 4}}], 50]]]] - 1) / 2
            ] < 10^-30
        ]
    ],
    True,
    {},
    TestID -> "ExactPureNegativity-QutritQubit"
]

EndTestSection[]
