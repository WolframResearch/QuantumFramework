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
   I-concurrence Sqrt[2 (d-1)/d] = 2/Sqrt[3] for d = 3 *)
VerificationTest[
    Chop[
        conc[QuantumState[QuantumState[Normalize[Flatten[IdentityMatrix[3]]], {3, 3}]["DensityMatrix"], {3, 3}]] -
            2 / Sqrt[3]
    ],
    0,
    TestID -> "Concurrence-QutritMaxEntangled"
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


EndTestSection[]
