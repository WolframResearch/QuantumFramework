(* ::Package:: *)

BeginTestSection["QuantumEntanglementSpec"]


(* Behavior tests for QuantumEntangledQ and QuantumEntanglementMonotone, one or more per rule of the
   entanglement spec (ENTANGLEMENT-SPEC-DRAFT.md, approved 2026-10-18). Test IDs follow its test list:
   A the state, B the split, C the method, D the True/False answers, E numbers, F measure values.
   Every reference value is computed below in plain Wolfram Language, without QuantumFramework.
   A test that expects "a message and a Failure" requires at least one message ({__}); the message names
   are pinned when the fix that issues them is written. *)


(* ========== References (plain Wolfram Language) ========== *)

proj[v_] := Outer[Times, v, Conjugate[v]];
rayDM[v_] := proj[v] / (Conjugate[v] . v);
blochDM[{x_, y_, z_}] := (IdentityMatrix[2] + x PauliMatrix[1] + y PauliMatrix[2] + z PauliMatrix[3]) / 2;

(* rho_{(ij),(kl)} reshaped to t[i, j, k, l] *)
ptB[rho_, {dA_, dB_}] := ArrayReshape[Transpose[ArrayReshape[rho, {dA, dB, dA, dB}], {1, 4, 3, 2}], {dA dB, dA dB}];
realign[rho_, {dA_, dB_}] := ArrayReshape[Transpose[ArrayReshape[rho, {dA, dB, dA, dB}], {1, 3, 2, 4}], {dA^2, dB^2}];
negRef[rho_, d_] := (Total[Abs[Eigenvalues[ptB[rho, d]]]] - 1) / 2;
ccnrRef[rho_, d_] := Total[SingularValueList[realign[rho, d]]] - 1;
trB[rho_, {dA_, dB_}] := TensorContract[ArrayReshape[rho, {dA, dB, dA, dB}], {{2, 4}}];
trA[rho_, {dA_, dB_}] := TensorContract[ArrayReshape[rho, {dA, dB, dA, dB}], {{1, 3}}];
entropyBits[m_] := With[{p = Select[Re @ Eigenvalues[m], # > 10^-13 &]}, - Total[p Log2[p]]];

(* textbook Wootters concurrence through the eigenvalues of rho.rho-tilde *)
woottersRef[r_] := With[{y2 = KroneckerProduct[PauliMatrix[2], PauliMatrix[2]]},
    With[{ev = Sqrt[Clip[ReverseSort[Re @ Eigenvalues[r . y2 . Conjugate[r] . y2]], {0, Infinity}]]},
        Max[0, ev[[1]] - ev[[2]] - ev[[3]] - ev[[4]]]
    ]
];

(* pure-state concurrence Sqrt[2 (1 - Tr rhoA^2)] from the coefficient matrix *)
pureConcRef[v_, {dA_, dB_}] := With[{m = ArrayReshape[v, {dA, dB}]},
    Sqrt[2 (1 - Re @ Tr[MatrixPower[m . ConjugateTranspose[m], 2]])]
];

(* J = S(A) - Sum_b p_b S(A|b), with B measured in its computational basis, and discord D = I - J, in bits *)
discordRefs[rho_, {dA_, dB_}] := With[{
    sA = entropyBits[trB[rho, {dA, dB}]], sB = entropyBits[trA[rho, {dA, dB}]], sAB = entropyBits[rho],
    branches = Table[With[{pb = KroneckerProduct[IdentityMatrix[dA], proj[UnitVector[dB, b]]]}, pb . rho . pb], {b, dB}]
},
    With[{j = sA - Total[If[Re[Tr[#]] > 10^-13, Re[Tr[#]] entropyBits[trB[# / Tr[#], {dA, dB}]], 0] & /@ branches]},
        {j, sA + sB - sAB - j}
    ]
];

randomMixed[d_] := With[{m = RandomComplex[{-1 - I, 1 + I}, {d, d}]}, # / Tr[#] & [m . ConjugateTranspose[m]]];
randomKet[d_] := Normalize[RandomComplex[{-1 - I, 1 + I}, d]];
haarU[d_] := RandomVariate[CircularUnitaryMatrixDistribution[d]];
(* an orthonormal frame per qudit: the columns of each unitary *)
localFrame[us_List] := QuantumBasis[QuantumTensorProduct @@
    (QuditBasis[AssociationThread[Subscript["e", #] & /@ Range[Length[#]], Transpose[#]]] & /@ us)];

within[x_, y_, tol_] := NumericQ[N[x]] && NumericQ[N[y]] && Abs[N[x] - N[y]] < tol;


(* ========== States ========== *)

bell = QuantumState["PhiPlus"];
bellY = QuantumState[bell, QuantumBasis["PauliY", 2]];
ghz3 = QuantumState["GHZ"[3]];
ghz4 = QuantumState["GHZ"[4]];
w3 = QuantumState["W"[3]];
w3m = Normal @ w3["DensityMatrix"];
g3m = Normal @ ghz3["DensityMatrix"];
rhoDegM = (w3m + g3m) / 2;            (* spectrum {1/2, 1/2, 0, ...}: a repeated nonzero eigenvalue *)
rhoNDM = (w3m + 2 g3m) / 3;           (* spectrum {2/3, 1/3, 0, ...} *)
rhoDeg = QuantumState[rhoDegM, {2, 2, 2}];
rhoND = QuantumState[rhoNDM, {2, 2, 2}];
inFrame[qs_, name_] := QuantumState[qs, QuantumBasis[name, 3]];
skew = QuditBasis[<|"a" -> {1, 0}, "b" -> {1, 1} / Sqrt[2]|>];
skewFrame = QuantumBasis[QuantumTensorProduct[skew, skew, skew]];
mm4 = QuantumState[IdentityMatrix[4] / 4, {2, 2}];
classical = QuantumState[DiagonalMatrix[{1, 0, 0, 1}] / 2, {2, 2}];
singlet = {0, 1, -1, 0} / Sqrt[2];
werner[p_] := p proj[singlet] + (1 - p) IdentityMatrix[4] / 4;     (* entangled exactly for p > 1/3 *)
zR3 = Expand[(Sqrt[3 + 2 Sqrt[2]] - Sqrt[2])^2] - 1;               (* exactly 0 *)

(* separable by construction: mixtures of product states *)
exactSep22 = 1/3 KroneckerProduct[blochDM[{1/2, 0, 1/3}], blochDM[{0, 1/5, 1/2}]] +
    2/3 KroneckerProduct[blochDM[{1/7, 1/3, 0}], blochDM[{1/2, 1/4, -1/3}]];
exactSep23 = 1/4 KroneckerProduct[blochDM[{1/2, 0, 1/3}], rayDM[{1, 2, 0}]] +
    3/4 KroneckerProduct[blochDM[{0, 1/3, -1/2}], rayDM[{1, 1, 3}]];
exactSep33 = 1/2 KroneckerProduct[rayDM[{1, 2, 0}], rayDM[{0, 1, 1}]] +
    1/2 KroneckerProduct[rayDM[{1, 1, 1}], rayDM[{2, 0, 1}]];
transcendentalSep22 = 1/2 KroneckerProduct[blochDM[{Sin[1]/2, 0, Cos[1]/3}], blochDM[{0, 1/5, 1/2}]] +
    1/2 KroneckerProduct[blochDM[{1/7, Cos[2]/3, 0}], blochDM[{1/2, 1/4, -1/3}]];

(* P. Horodecki, Phys. Lett. A 232, 333 (1997): PPT entangled 3x3 (a = 1/2, realignment detects it) and
   2x4 (b = 1/2, realignment does not) *)
horodecki33 = With[{a = 1/2}, (1/(8 a + 1)) {
    {a, 0, 0, 0, a, 0, 0, 0, a}, {0, a, 0, 0, 0, 0, 0, 0, 0}, {0, 0, a, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, a, 0, 0, 0, 0, 0}, {a, 0, 0, 0, a, 0, 0, 0, a}, {0, 0, 0, 0, 0, a, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, (1 + a)/2, 0, Sqrt[1 - a^2]/2}, {0, 0, 0, 0, 0, 0, 0, a, 0},
    {a, 0, 0, 0, a, 0, Sqrt[1 - a^2]/2, 0, (1 + a)/2}
}];
horodecki24 = With[{b = 1/2}, (1/(7 b + 1)) {
    {b, 0, 0, 0, 0, b, 0, 0}, {0, b, 0, 0, 0, 0, b, 0}, {0, 0, b, 0, 0, 0, 0, b}, {0, 0, 0, b, 0, 0, 0, 0},
    {0, 0, 0, 0, (1 + b)/2, 0, 0, Sqrt[1 - b^2]/2}, {b, 0, 0, 0, 0, b, 0, 0}, {0, b, 0, 0, 0, 0, b, 0},
    {0, 0, b, 0, Sqrt[1 - b^2]/2, 0, 0, (1 + b)/2}
}];

(* Breuer family on N x N (N even): lambda |psi0><psi0| + (1 - lambda) 2 P_sym / (N (N + 1)), PPT entangled for
   0 < lambda <= 1/(N + 2) (Breuer 2006), and realignment-blind for lambda <= 1/N *)
breuerState[lambda_, n_] := With[{
    psi = Normal @ SparseArray[Table[{(k - 1) n + (n + 1 - k)} -> (-1)^(k - 1) / Sqrt[n], {k, n}], {n^2}],
    swap = Normal @ SparseArray[Flatten @ Table[{(i - 1) n + j, (j - 1) n + i} -> 1, {i, n}, {j, n}], {n^2, n^2}]
},
    lambda proj[psi] + (1 - lambda) (IdentityMatrix[n^2] + swap) / (n (n + 1))
];

(* the issue 23 family: a mixture of three product projectors with one tiny weight r *)
productProjector[{d1_, d2_}] := proj[Flatten @ KroneckerProduct[randomKet[d1], randomKet[d2]]];
nearBoundarySeparable[dims_, r_] := With[{ps = Table[productProjector[dims], 3]}, (1 - r) / 2 (ps[[1]] + ps[[2]]) + r ps[[3]]];

productDMs[] := (SeedRandom[101]; Catenate @ Table[QuantumState[productProjector[dims], dims], {dims, {{2, 3}, {3, 3}, {4, 4}}}, {50}]);
productKets[] := (SeedRandom[202]; Catenate @ Table[
    QuantumState[Flatten @ KroneckerProduct[randomKet[dims[[1]]], randomKet[dims[[2]]]], dims], {dims, {{2, 2}, {2, 3}, {3, 3}}}, {50}]);
nearBoundaryStates[] := (SeedRandom[5005]; Flatten[Table[
    QuantumState[nearBoundarySeparable[dims, r], dims], {dims, {{2, 3}, {3, 3}}}, {r, {10^-9, 10^-8, 10^-7}}, {15}], 2]);

countTrue[states_, method_] := Count[QuantumEntangledQ[#, Automatic, method] & /@ states, True];
countTrueDefault[states_] := Count[QuantumEntangledQ /@ states, True];


(* ========== A. The state (spec 2.1) ========== *)

(* a mixed state stored in another orthonormal basis is the same state *)
VerificationTest[
    With[{v = TimeConstrained[QuantumEntanglementMonotone[inFrame[rhoND, "PauliX"], Automatic, "Negativity"], 60, $TimedOut]},
        within[v, negRef[N @ rhoNDM, {2, 4}], 10^-10]
    ],
    True,
    {___},
    TestID -> "A1-MixedStateInPauliXBasis-Negativity"
]

VerificationTest[
    With[{v = TimeConstrained[QuantumEntanglementMonotone[inFrame[rhoDeg, "PauliX"], Automatic, "Negativity"], 60, $TimedOut]},
        within[v, negRef[N @ rhoDegM, {2, 4}], 10^-10]
    ],
    True,
    {___},
    TestID -> "A2-DegenerateMixedStateInPauliXBasis-Negativity"
]

VerificationTest[
    SeedRandom[3];
    Max @ Table[
        With[{qs = QuantumState[QuantumState[N @ rhoDegM, {2, 2, 2}], localFrame[{haarU[2], haarU[2], haarU[2]}]]},
            Abs[QuantumEntanglementMonotone[qs, Automatic, "Negativity"] - negRef[N @ rhoDegM, {2, 4}]]
        ],
        200
    ] < 10^-8,
    True,
    {___},
    TestID -> "A3-MachineMixedStateInRandomFrames-Negativity"
]

VerificationTest[
    {QuantumEntanglementMonotone[bellY, "Negativity"], Chop @ N @ QuantumEntanglementMonotone[bellY, "Realignment"]},
    {1/2, 1.},
    {___},
    TestID -> "A4-BellInPauliYBasis-Values"
]

VerificationTest[
    TimeConstrained[
        With[{s = inFrame[rhoDeg, "PauliX"]["SchmidtBasis", 2]},
            FailureQ[s] || With[{m = Normal @ s["Computational"]["DensityMatrix"]},
                Max[Abs[N @ Flatten[m / Tr[m] - rhoDegM]]] < 10^-10
            ]
        ],
        60, $TimedOut
    ],
    True,
    {___},
    TestID -> "A5-MixedSchmidtBasis-IsTheState"
]

VerificationTest[
    FailureQ @ QuantumEntanglementMonotone[QuantumState[{zR3, 0, 0, 0}, {2, 2}], "Negativity"],
    True,
    {__},
    TestID -> "A6-ZeroNormKet-Failure"
]

VerificationTest[
    FailureQ @ QuantumEntanglementMonotone[
        QuantumState[{1, 0, -1, 0}, QuantumTensorProduct[QuantumBasis[{{1, 0}, {1, 0}}], QuantumBasis[2]]],
        "Negativity"
    ],
    True,
    {__},
    TestID -> "A7-DependentBasis-Failure"
]

VerificationTest[
    With[{qs = QuantumState[QuantumState[N @ rhoNDM, {2, 2, 2}], skewFrame]},
        Append[
            FailureQ @ QuantumEntanglementMonotone[qs, Automatic, #] & /@ {"Negativity", "LogNegativity", "Realignment", "Concurrence"},
            FailureQ @ QuantumEntangledQ[qs]
        ]
    ],
    {True, True, True, True, True},
    {__},
    TestID -> "A8-NonOrthonormalBasis-Failure"
]

(* the flat-list split of a mixed state stored in another basis: qubit 2 against qubits 1 and 3 *)
swap12[rho_] := ArrayReshape[Transpose[ArrayReshape[rho, ConstantArray[2, 6]], {2, 1, 3, 5, 4, 6}], {8, 8}];

VerificationTest[
    With[{v = TimeConstrained[QuantumEntanglementMonotone[inFrame[rhoND, "PauliY"], {2}, "Negativity"], 60, $TimedOut]},
        within[v, negRef[N @ swap12[rhoNDM], {2, 4}], 10^-10]
    ],
    True,
    {___},
    TestID -> "A9-MixedStateInPauliYBasis-FlatListSplit"
]

VerificationTest[
    With[{v = TimeConstrained[QuantumEntanglementMonotone[inFrame[rhoND, "PauliY"], Automatic, "Realignment"], 60, $TimedOut]},
        within[v, ccnrRef[N @ rhoNDM, {2, 4}], 10^-10]
    ],
    True,
    {___},
    TestID -> "A10-MixedStateInPauliYBasis-Realignment"
]

(* one 4-level qudit holding a mixed two-qubit state, stored in the Fourier[4] basis, split 2 x 2 *)
VerificationTest[
    With[{v = TimeConstrained[
            QuantumEntanglementMonotone[QuantumState[QuantumState[werner[2/3], 4], QuantumBasis["Fourier"[4]]], Automatic, "Negativity"],
            60, $TimedOut]},
        within[v, negRef[N @ werner[2/3], {2, 2}], 10^-10]
    ],
    True,
    {___},
    TestID -> "A11-MixedQuditInFourierBasis-Negativity"
]

(* one 32-level qudit holding a machine mixed state, stored in the Fourier[32] basis, split 4 x 8: deciding that
   the frame can be changed must not take the exact rank of the 32 x 32 Fourier matrix *)
VerificationTest[
    SeedRandom[3];
    With[{rho = randomMixed[32]},
        With[{v = TimeConstrained[
                QuantumEntanglementMonotone[QuantumState[QuantumState[rho, 32], QuantumBasis["Fourier"[32]]], 4, "Negativity"],
                20, $TimedOut]},
            within[v, negRef[rho, {4, 8}], 10^-10]
        ]
    ],
    True,
    TestID -> "A12-MachineMixedQuditInFourier32Basis-Negativity"
]

(* two qutrits holding a machine mixed state, stored in the Fourier[3] (x) Fourier[3] basis: the d > 2 concurrence
   bound changes under local unitaries, so it must be read in the computational frame, as for the same state stored there *)
VerificationTest[
    SeedRandom[17];
    With[{rho = randomMixed[9]},
        within[
            QuantumEntanglementMonotone[QuantumState[QuantumState[rho, {3, 3}], QuantumBasis[QuantumTensorProduct[QuditBasis["Fourier"[3]], QuditBasis["Fourier"[3]]]]], "Concurrence"],
            QuantumEntanglementMonotone[QuantumState[rho, {3, 3}], "Concurrence"],
            10^-10
        ]
    ],
    True,
    TestID -> "A13-MixedTwoQutritsInFourierBasis-Concurrence"
]

(* the exact isotropic two-qutrit state (p = 1/2) stored in the Fourier[3] (x) Fourier[3] basis: the change of basis
   leaves long sums of cube roots of unity, which must be simplified before an exact trace norm *)
VerificationTest[
    With[{iso = With[{phi = Flatten[IdentityMatrix[3]] / Sqrt[3]}, proj[phi] / 2 + IdentityMatrix[9] / 18]},
        within[
            TimeConstrained[
                QuantumEntanglementMonotone[QuantumState[QuantumState[iso, {3, 3}], QuantumBasis[QuantumTensorProduct[QuditBasis["Fourier"[3]], QuditBasis["Fourier"[3]]]]], "Negativity"],
                30, $TimedOut],
            negRef[N @ iso, {3, 3}],
            10^-10
        ]
    ],
    True,
    TestID -> "A14-ExactMixedTwoQutritsInFourierBasis-Negativity"
]


(* ========== B. The split (spec 2.3) ========== *)

VerificationTest[
    QuantumEntangledQ[ghz3, {{1}, {3}}],
    False,
    TestID -> "B1-GHZ3-PairSplit-Separable"
]

(* a Bell pair on qubits 1 and 3, |0> on 2 and 4: the flat list {1, 2} cuts the pair, {1, 3} keeps it inside one part *)
VerificationTest[
    With[{qs = QuantumState[Normal @ SparseArray[{1 -> 1 / Sqrt[2], 11 -> 1 / Sqrt[2]}, 16], {2, 2, 2, 2}]},
        {QuantumEntangledQ[qs, {1, 2}], QuantumEntangledQ[qs, {1, 3}]}
    ],
    {True, False},
    TestID -> "B2-FlatListSplitFollowsQuditOrder"
]

VerificationTest[
    QuantumEntangledQ[QuantumState[{1, 0, 0, 1} / Sqrt[2], 4], 2],
    True,
    TestID -> "B3-SingleQuditIntegerSplit"
]

VerificationTest[
    {FailureQ @ QuantumEntangledQ[ghz3, {{1}, {4}}], FailureQ @ QuantumEntanglementMonotone[ghz3, {{1}, {4}}, "Concurrence"]},
    {True, True},
    {__},
    TestID -> "B4-IndexOutOfRange-Failure"
]

VerificationTest[
    {FailureQ @ QuantumEntangledQ[bell, {1, 1}], FailureQ @ QuantumEntanglementMonotone[bell, {1, 1}, "Negativity"]},
    {True, True},
    {__},
    TestID -> "B5-RepeatedIndex-Failure"
]

VerificationTest[
    FailureQ @ QuantumEntangledQ[ghz4, {{1, 2}, {2, 3}}],
    True,
    {__},
    TestID -> "B6-OverlappingGroups-Failure"
]

VerificationTest[
    FailureQ @ QuantumEntangledQ[ghz3, {{}, {1}}],
    True,
    {__},
    TestID -> "B7-EmptyGroup-Failure"
]

VerificationTest[
    {FailureQ @ QuantumEntangledQ[bell, 0], FailureQ @ QuantumEntangledQ[bell, -1]},
    {True, True},
    {__},
    TestID -> "B8-ZeroOrNegativeSplit-Failure"
]

VerificationTest[
    FailureQ @ QuantumEntangledQ[bell, 3],
    True,
    {__},
    TestID -> "B9-IntegerOnTwoQudits-Failure"
]

VerificationTest[
    {
        FailureQ @ QuantumEntangledQ[QuantumState["0"]],
        FailureQ @ QuantumEntangledQ[QuantumState[{1, 0, 0}]],
        FailureQ @ QuantumEntanglementMonotone[QuantumState["0"], "Negativity"]
    },
    {True, True, True},
    {__},
    TestID -> "B10-SingleQuditNoSplit-Failure"
]


(* ========== C. The method argument (spec 2.4) ========== *)

VerificationTest[
    {QuantumEntangledQ[bell, "Negativity"], QuantumEntangledQ[bell, Automatic, "Negativity"]},
    {True, True},
    TestID -> "C1-MethodInSecondSlot"
]

VerificationTest[
    QuantumEntanglementMonotone[bell, "Negativity"],
    1/2,
    TestID -> "C2-MonotoneMethodInSecondSlot"
]

VerificationTest[
    {FailureQ @ QuantumEntangledQ[bell, Automatic, "PPT"], FailureQ @ QuantumEntanglementMonotone[bell, "PPT"]},
    {True, True},
    {__},
    TestID -> "C3-UnknownMethod-Failure"
]

VerificationTest[
    FailureQ @ QuantumEntangledQ[bell, Automatic, #] & /@ {"MutualInformationI", "MutualInformationJ", "Discord"},
    {True, True, True},
    {__},
    TestID -> "C4-CorrelationMeasuresRejected"
]

VerificationTest[
    {QuantumEntangledQ[bell, Automatic, "EntanglementEntropy"], QuantumEntangledQ[QuantumState["0+"], Automatic, "EntanglementEntropy"]},
    {True, False},
    TestID -> "C5a-EntropyMethodDecidesPureStates"
]

VerificationTest[
    QuantumEntangledQ[mm4, Automatic, "EntanglementEntropy"],
    Indeterminate,
    {__},
    TestID -> "C5b-EntropyMethodOnMixedState"
]

VerificationTest[
    {
        QuantumEntangledQ[bell, {"RenyiEntropy", 1/2}],
        QuantumEntangledQ[bell, Automatic, {"RenyiEntropy", 1/2}],
        N @ QuantityMagnitude @ QuantumEntanglementMonotone[bell, {"RenyiEntropy", 1/2}],
        N @ QuantityMagnitude @ QuantumEntanglementMonotone[bell, Automatic, {"RenyiEntropy", 1/2}]
    },
    {True, True, 1., 1.},
    TestID -> "C6-RenyiListMethodInBothSlots"
]


(* ========== D. QuantumEntangledQ answers (spec 3) ========== *)

VerificationTest[
    {QuantumEntangledQ[bell], QuantumEntangledQ[QuantumState["0+"]]},
    {True, False},
    TestID -> "D1-BellAndProduct"
]

VerificationTest[
    {QuantumEntangledQ[QuantumState[werner[1/2], {2, 2}]], QuantumEntangledQ[QuantumState[werner[1/3], {2, 2}]]},
    {True, False},
    TestID -> "D2-ExactWernerAboveAndAtThreshold"
]

VerificationTest[
    QuantumEntangledQ[QuantumState[KroneckerProduct[DiagonalMatrix[{1, 2}] / 3, DiagonalMatrix[{1, 1, 2}] / 4], {2, 3}]],
    False,
    TestID -> "D3-Separable2x3"
]

VerificationTest[
    QuantumEntangledQ[QuantumState[N @ horodecki33, {3, 3}]],
    True,
    TestID -> "D4-Horodecki3x3-RealignmentDetects"
]

VerificationTest[
    QuantumEntangledQ[QuantumState[N @ breuerState[3/20, 4], {4, 4}]],
    Indeterminate,
    {___},
    TestID -> "D5-Breuer4x4-Undecided"
]

VerificationTest[
    QuantumEntangledQ[QuantumState[N @ horodecki24, {2, 4}]],
    Indeterminate,
    {___},
    TestID -> "D6-Horodecki2x4-Undecided"
]

VerificationTest[
    QuantumEntangledQ[QuantumState[N[IdentityMatrix[9] / 9], {3, 3}]],
    False,
    TestID -> "D7-MaximallyMixed3x3-PurityBound"
]

VerificationTest[
    QuantumEntangledQ[QuantumState[N @ exactSep33, {3, 3}]],
    Indeterminate,
    {___},
    TestID -> "D8-Separable3x3AbovePurityBound-Undecided"
]

VerificationTest[
    QuantumEntangledQ[QuantumState[Normalize @ Flatten @ KroneckerProduct[{1, 2, 0}, {0, 1, 1}], 3]],
    False,
    TestID -> "D9-PureProduct3x3"
]

VerificationTest[
    TimeConstrained[QuantumEntangledQ[QuantumState[exactSep23, {2, 3}]], 60, $TimedOut],
    False,
    TestID -> "D10-ExactSeparable2x3-Decided"
]

VerificationTest[
    MatchQ[TimeConstrained[QuantumEntangledQ[QuantumState[exactSep33, {3, 3}]], 60, $TimedOut], True | False | Indeterminate],
    True,
    {___},
    TestID -> "D11-ExactSeparable3x3-NeverUnevaluated"
]


(* ========== E. Numbers (spec 4) ========== *)

VerificationTest[
    countTrue[productDMs[], "Concurrence"],
    0,
    {___},
    TestID -> "E1-MachineProductDensityMatrices-Concurrence"
]

VerificationTest[
    countTrue[productKets[], "Concurrence"],
    0,
    {___},
    TestID -> "E2-MachineProductKets-Concurrence"
]

VerificationTest[
    countTrue[nearBoundaryStates[], "Concurrence"],
    0,
    {___},
    TestID -> "E3-NearBoundarySeparableMixtures-Concurrence"
]

VerificationTest[
    With[{states = Join[productDMs[], productKets[], nearBoundaryStates[]]},
        {countTrue[states, "Negativity"], countTrue[states, "Realignment"], countTrueDefault[states]}
    ],
    {0, 0, 0},
    {___},
    TestID -> "E4-SameStatesOtherMethods"
]

VerificationTest[
    QuantumEntangledQ[QuantumState[N @ werner[1/3 + 10^-6], {2, 2}]],
    True,
    TestID -> "E5-MachineWernerClearlyEntangled"
]

VerificationTest[
    QuantumEntangledQ[QuantumState[werner[1/3 + 10^-12], {2, 2}]],
    True,
    TestID -> "E6-ExactWernerBarelyEntangled"
]

VerificationTest[
    QuantumEntangledQ[QuantumState[werner[1/3 + 10.^-12], {2, 2}]],
    False,
    TestID -> "E7-MachineWernerWithinRoundoff"
]

VerificationTest[
    With[{t = AbsoluteTiming[TimeConstrained[QuantumEntangledQ[QuantumState[exactSep22, {2, 2}]], 60, $TimedOut]]},
        {Last[t], First[t] < 15}
    ],
    {False, True},
    TestID -> "E8-ExactSeparable2x2-RadicalZero"
]

VerificationTest[
    With[{t = AbsoluteTiming[TimeConstrained[QuantumEntangledQ[QuantumState[transcendentalSep22, {2, 2}]], 60, $TimedOut]]},
        {MatchQ[Last[t], False | Indeterminate], First[t] < 15}
    ],
    {True, True},
    {___},
    TestID -> "E9-ExactTranscendentalSeparable-Bounded"
]

VerificationTest[
    Block[{p}, QuantumEntangledQ[QuantumState[werner[p], {2, 2}]]],
    Indeterminate,
    {__},
    TestID -> "E10-SymbolicWerner-Indeterminate"
]

VerificationTest[
    Block[{phi}, QuantumEntangledQ[QuantumState[{1, 0, 0, Exp[I phi]} / Sqrt[2]]]],
    True,
    {__},
    TestID -> "E11-SymbolicMaximallyEntangledFamily-True"
]

VerificationTest[
    Block[{theta}, QuantumEntangledQ[QuantumState[{Cos[theta], Sin[theta], 0, 0}]]],
    False,
    {__},
    TestID -> "E12-SymbolicProductFamily-False"
]

VerificationTest[
    Block[{a, b, c},
        With[{v = TimeConstrained[
                QuantumEntanglementMonotone[QuantumState[{0, a, b, 0, c, 0, 0, 0}, {2, 2, 2}], {{1}, {2}}, "Concurrence"],
                120, $TimedOut],
            point = Thread[{a, b, c} -> Normalize[{0.3 + 0.2 I, 0.5 - 0.1 I, 0.4 + 0.6 I}]]},
            within[v /. point, 2 Abs[b c] /. point, 10^-10]
        ]
    ],
    True,
    TestID -> "E13-SymbolicPairConcurrence-Evaluates"
]

VerificationTest[
    TimeConstrained[QuantumEntanglementMonotone[QuantumState[exactSep22, {2, 2}], "Concurrence"], 60, $TimedOut],
    0,
    TestID -> "E14-ExactMixedTwoQubitConcurrence"
]

(* a 30-digit pure 3-qubit state whose pair {1, 2} is 4/5 |Phi+><Phi+| + 1/5 |01><01|, concurrence 4/5: the pair is
   mixed, and the zero entries its partial trace leaves must not make it read as pure *)
VerificationTest[
    With[{v = N[Sqrt[4/5] Flatten[KroneckerProduct[{1, 0, 0, 1} / Sqrt[2], {1, 1} / Sqrt[2]]] + Sqrt[1/5] Flatten[KroneckerProduct[{0, 1, 0, 0}, {1, -1} / Sqrt[2]]], 30]},
        within[QuantumEntanglementMonotone[QuantumState[v, {2, 2, 2}], {{1}, {2}}, "Concurrence"], 4/5, 10^-20]
    ],
    True,
    TestID -> "E15-HighPrecisionMixedPairConcurrence"
]

(* near-pure separable two-qubit mixtures (1 - e) P1 + e P2 of product projectors, machine numbers, e from 3e-14 to 1e-9 *)
VerificationTest[
    SeedRandom[606];
    countTrue[Flatten @ Table[QuantumState[(1 - e) productProjector[{2, 2}] + e productProjector[{2, 2}], {2, 2}], {e, {3*^-14, 1*^-13, 1*^-12, 1*^-11, 1*^-10, 1*^-9}}, {35}], "Concurrence"],
    0,
    TestID -> "E16-NearPureSeparableTwoQubits-Concurrence"
]

(* near-pure Bell-diagonal states (1 - e) Phi+ + e Phi-, concurrence 1 - 2 e: mixed by less than the purity test can
   see, they must not be read as pure; machine numbers with e = 1e-14, and 6-digit numbers with e = 1e-5 *)
VerificationTest[
    With[{phiP = proj[{1, 0, 0, 1} / Sqrt[2]], phiM = proj[{1, 0, 0, -1} / Sqrt[2]]},
        {
            within[QuantumEntanglementMonotone[QuantumState[N[(1 - 10^-14) phiP + 10^-14 phiM], {2, 2}], "Concurrence"], 1 - 2 10^-14, 5 10^-15],
            within[QuantumEntanglementMonotone[QuantumState[SetPrecision[(1 - 10^-5) phiP + 10^-5 phiM, 6], {2, 2}], "Concurrence"], 1 - 2 10^-5, 5 10^-6]
        }
    ],
    {True, True},
    TestID -> "E17-NearPureMixedStates-Concurrence"
]

(* an exact separable mixture of three pure product states with rational and Gaussian-rational amplitudes: its
   Wootters combination is a sum of square roots of cubic roots that cancels to exactly zero *)
VerificationTest[
    With[{k = {{1, 2}, {3, -1}, {2, 1 + I}}, l = {{1, 1}, {1, -2}, {2 - I, 1}}},
        QuantumEntanglementMonotone[
            QuantumState[Total[MapThread[#3 proj[Normalize[Flatten @ KroneckerProduct[#1, #2]]] &, {k, l, {1/2, 1/3, 1/6}}]], {2, 2}],
            "Concurrence"
        ]
    ],
    0,
    TestID -> "E18-ExactThreeProductMixture-Concurrence"
]


(* ========== F. Measure values (spec 5) ========== *)

VerificationTest[
    SeedRandom[7];
    Max @ Catenate @ Table[
        With[{rho = randomMixed[Times @@ dims]},
            With[{qs = QuantumState[rho, dims], n = negRef[rho, dims]},
                {Abs[QuantumEntanglementMonotone[qs, "Negativity"] - n], Abs[QuantumEntanglementMonotone[qs, "LogNegativity"] - Log2[2 n + 1]]}
            ]
        ],
        {dims, {{2, 3}, {3, 3}, {2, 4}}}, {5}
    ] < 10^-10,
    True,
    TestID -> "F1-NegativityMatchesReference"
]

VerificationTest[
    Chop @ N[{QuantumEntanglementMonotone[ghz3, {{1}, {2}}, "Concurrence"], QuantumEntanglementMonotone[w3, {{1}, {3}}, "Concurrence"]} - {0, 2/3}],
    {0, 0},
    TestID -> "F2-PairConcurrenceOfPureStates"
]

VerificationTest[
    SeedRandom[11];
    With[{v = Flatten @ KroneckerProduct[randomKet[2^9], {1, 0}]},
        within[
            TimeConstrained[QuantumEntanglementMonotone[QuantumState[v, ConstantArray[2, 10]], {Range[4], Range[5, 9]}, "Concurrence"], 120, $TimedOut],
            pureConcRef[v[[1 ;; ;; 2]], {16, 32}],
            10^-10
        ]
    ],
    True,
    TestID -> "F3-LargePurePairConcurrence"
]

(* the same with exact amplitudes: the pure pair must not fall to the exact mixed-state route *)
VerificationTest[
    With[{v = Flatten @ KroneckerProduct[Normalize[Table[Mod[k^2, 7] + 1, {k, 2^7}]], {1, 0}]},
        within[
            TimeConstrained[QuantumEntanglementMonotone[QuantumState[v, ConstantArray[2, 8]], {Range[3], Range[4, 7]}, "Concurrence"], 60, $TimedOut],
            pureConcRef[v[[1 ;; ;; 2]], {8, 16}],
            10^-10
        ]
    ],
    True,
    TestID -> "F3b-ExactLargePurePairConcurrence"
]

(* a rank-3 entangled 2 x 3 mixed state given as machine numbers: its zero eigenvalues come out as round-off, which
   must not enter the value; reference: the generator-pair bound at 60 digits *)
VerificationTest[
    With[{
        rho = # / Tr[#] & [proj[UnitVector[6, 1] + UnitVector[6, 6]] + proj[UnitVector[6, 2] + UnitVector[6, 3]] / 2 + proj[UnitVector[6, 1] + UnitVector[6, 2] + I UnitVector[6, 5]] / 5],
        gens = Function[n, Catenate @ Table[SparseArray[{{j, k} -> -I, {k, j} -> I}, {n, n}], {k, 2, n}, {j, k - 1}]]
    },
        within[
            QuantumEntanglementMonotone[QuantumState[N @ rho, {2, 3}], "Concurrence"],
            With[{r = N[rho, 60]}, Norm @ Flatten @ Outer[
                Function[{a, b}, With[{o = Normal @ KroneckerProduct[a, b]},
                    With[{l = Sqrt[Clip[Re @ Eigenvalues[r . o . Conjugate[r] . o], {0, Infinity}]]}, Max[0, 2 Max[l] - Total[l]]]]],
                gens[2], gens[3], 1]],
            10^-13
        ]
    ],
    True,
    TestID -> "F15-RankDeficientMixedConcurrence"
]

VerificationTest[
    SeedRandom[13];
    Max @ Table[With[{rho = randomMixed[4]}, Abs[QuantumEntanglementMonotone[QuantumState[rho, {2, 2}], "Concurrence"] - woottersRef[rho]]], {10}] < 10^-8,
    True,
    TestID -> "F4-TwoQubitConcurrenceMatchesWootters"
]

VerificationTest[
    QuantumEntanglementMonotone[mm4, "Realignment"],
    -1/2,
    {__},
    TestID -> "F5-RealignmentRawValueWithWarning"
]

VerificationTest[
    QuantumEntangledQ[QuantumState[N @ horodecki33, {3, 3}], Automatic, "Realignment"],
    True,
    TestID -> "F6-RealignmentInsideEntangledQIsSilent"
]

VerificationTest[
    SeedRandom[17];
    Max @ Append[
        Table[With[{rho = randomMixed[Times @@ dims]}, Abs[QuantumEntanglementMonotone[QuantumState[rho, dims], "Realignment"] - ccnrRef[rho, dims]]],
            {dims, {{3, 2}, {2, 4}, {3, 3}}}, {3}] // Catenate,
        With[{rho = randomMixed[4]},
            Abs[QuantumEntanglementMonotone[QuantumState[QuantumState[rho, {2, 2}], QuantumBasis["PauliY", 2]], "Realignment"] - ccnrRef[rho, {2, 2}]]]
    ] < 10^-10,
    True,
    {___},
    TestID -> "F7-RealignmentMatchesReference"
]

VerificationTest[
    FailureQ @ QuantumEntanglementMonotone[ghz3, {{1}, {4}}, "Realignment"],
    True,
    {__},
    TestID -> "F8-RealignmentInvalidSplit-Failure"
]

VerificationTest[
    QuantumEntanglementMonotone[bell, Automatic, {"RenyiEntropy", 1}],
    Quantity[1, "Bits"],
    TestID -> "F9-RenyiOrderOneIsEntanglementEntropy"
]

VerificationTest[
    MatchQ[QuantumEntanglementMonotone[bell, Automatic, #], Quantity[_, "Bits"]] & /@
        {"EntanglementEntropy", "RenyiEntropy", "MutualInformationI", "MutualInformationJ", "Discord"},
    {True, True, True, True, True},
    TestID -> "F10-EntropiesInBits"
]

VerificationTest[
    SeedRandom[19];
    Max @ Catenate @ Table[
        With[{rho = randomMixed[Times @@ dims]},
            With[{qs = QuantumState[rho, dims], ref = discordRefs[rho, dims]},
                Abs[N @ QuantityMagnitude[{QuantumEntanglementMonotone[qs, "MutualInformationJ"], QuantumEntanglementMonotone[qs, "Discord"]}] - ref]
            ]
        ],
        {dims, {{2, 3}, {3, 3}}}, {3}
    ] < 10^-10,
    True,
    TestID -> "F11-DiscordWithQutritB"
]

VerificationTest[
    {FailureQ @ QuantumEntanglementMonotone[QuantumState["0"], "Discord"], FailureQ @ QuantumEntanglementMonotone[QuantumState["0"], "MutualInformationI"]},
    {True, True},
    {__},
    TestID -> "F12-CorrelationMeasuresOnOneQudit-Failure"
]

VerificationTest[
    SeedRandom[23];
    With[{rho = randomMixed[9], u = KroneckerProduct[haarU[3], IdentityMatrix[3]]},
        With[{qs = QuantumState[rho, {3, 3}], qsU = QuantumState[u . rho . ConjugateTranspose[u], {3, 3}]},
            {
                Chop[QuantumEntanglementMonotone[bell] - 1],
                Chop[QuantumEntanglementMonotone[QuantumState[N @ werner[4/5], {2, 2}]] - 7/10],
                Chop[QuantumEntanglementMonotone[qs] - QuantumEntanglementMonotone[qs, "Negativity"]],
                Chop[QuantumEntanglementMonotone[qsU] - QuantumEntanglementMonotone[qs]]
            }
        ]
    ],
    {0, 0, 0, 0},
    TestID -> "F13-DefaultMethod"
]

VerificationTest[
    QuantumEntanglementMonotone[classical, "MutualInformationI"],
    Quantity[1, "Bits"],
    TestID -> "F14-MutualInformationOfClassicalCorrelation"
]


EndTestSection[]
