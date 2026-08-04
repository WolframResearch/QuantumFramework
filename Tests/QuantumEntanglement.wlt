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

EndTestSection[]
