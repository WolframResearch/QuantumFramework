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

localUnitary[rho_, us_List] := With[{u = KroneckerProduct @@ us}, u . rho . ConjugateTranspose[u]];

(* the Horodecki states of horodecki33 and horodecki24 at any a, b in (0, 1), all PPT entangled *)
horodecki33Family[a_] := (1/(8 a + 1)) {
    {a, 0, 0, 0, a, 0, 0, 0, a}, {0, a, 0, 0, 0, 0, 0, 0, 0}, {0, 0, a, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, a, 0, 0, 0, 0, 0}, {a, 0, 0, 0, a, 0, 0, 0, a}, {0, 0, 0, 0, 0, a, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, (1 + a)/2, 0, Sqrt[1 - a^2]/2}, {0, 0, 0, 0, 0, 0, 0, a, 0},
    {a, 0, 0, 0, a, 0, Sqrt[1 - a^2]/2, 0, (1 + a)/2}
};
horodecki24Family[b_] := (1/(7 b + 1)) {
    {b, 0, 0, 0, 0, b, 0, 0}, {0, b, 0, 0, 0, 0, b, 0}, {0, 0, b, 0, 0, 0, 0, b}, {0, 0, 0, b, 0, 0, 0, 0},
    {0, 0, 0, 0, (1 + b)/2, 0, 0, Sqrt[1 - b^2]/2}, {b, 0, 0, 0, 0, b, 0, 0}, {0, b, 0, 0, 0, 0, b, 0},
    {0, 0, b, 0, Sqrt[1 - b^2]/2, 0, 0, (1 + b)/2}
};
(* Bennett et al., PRL 82, 5385 (1999): the normalized projector onto the complement of the Tiles UPB, rank 4 *)
tilesState = With[{vs = {
        Flatten @ KroneckerProduct[{1, 0, 0}, {1, -1, 0}], Flatten @ KroneckerProduct[{0, 0, 1}, {0, 1, -1}],
        Flatten @ KroneckerProduct[{1, -1, 0}, {0, 0, 1}], Flatten @ KroneckerProduct[{0, 1, -1}, {1, 0, 0}],
        Flatten @ KroneckerProduct[{1, 1, 1}, {1, 1, 1}]}},
    (IdentityMatrix[9] - Total[rayDM /@ vs]) / 4
];
(* a state of d1 x d2 qudits placed in the first levels of D1 x D2 qudits: the same entanglement, larger local spaces *)
embedState[rho_, {d1_, d2_}, {e1_, e2_}] := With[{idx = Flatten @ Table[(i - 1) e2 + j, {i, d1}, {j, d2}]},
    Normal @ SparseArray[Flatten @ MapIndexed[{idx[[First[#2]]], idx[[Last[#2]]]} -> #1 &, rho, {2}], {e1 e2, e1 e2}]
];
productMixture[dims_, weights_] := Total @ MapThread[#1 KroneckerProduct[randomMixed[First[dims]], randomMixed[Last[dims]]] &, {weights}];
exactProduct[dims_] := KroneckerProduct @@ (With[{m = RandomInteger[{-3, 3}, {#, #}] + I RandomInteger[{-3, 3}, {#, #}]}, # / Tr[#] &[m . ConjugateTranspose[m]]] & /@ dims);

skewQubit = QuditBasis[<|"a" -> {1, 0}, "b" -> {1, 1} / Sqrt[2]|>];
skewQutrit = QuditBasis[<|"a" -> {1, 0, 0}, "b" -> {1, 1, 0} / Sqrt[2], "c" -> {1, 1, 1} / Sqrt[3]|>];
skewLocalFrame[dims_] := QuantumBasis[QuantumTensorProduct @@ Replace[dims, {2 -> skewQubit, 3 -> skewQutrit}, {1}]];
namedFrames[dims_] := {
    QuantumBasis[QuantumTensorProduct @@ (QuditBasis["Fourier"[#]] & /@ dims)],
    QuantumBasis[QuantumTensorProduct @@ (If[# == 2, QuditBasis["PauliY"], QuditBasis["Fourier"[#]]] & /@ dims)]
};
randomFrames[dims_, n_] := Table[localFrame[haarU /@ dims], n];
inLocalFrame[state_, dims_, frame_] := QuantumState[QuantumState[state, dims], frame];
methodValues[qs_] := QuantumEntanglementMonotone[qs, Automatic, #] & /@ {"Negativity", "Realignment"};


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
            ! MatchQ[s, _QuantumState] || With[{m = Normal @ s["Computational"]["DensityMatrix"]},
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

(* a density matrix changed into a frame that is not orthonormal is the same state: every method reads its
   computational form, which is rhoNDM itself *)
VerificationTest[
    With[{qs = QuantumState[QuantumState[N @ rhoNDM, {2, 2, 2}], skewFrame], ref = QuantumState[N @ rhoNDM, {2, 2, 2}]},
        {
            within[QuantumEntanglementMonotone[qs, Automatic, "Negativity"], negRef[N @ rhoNDM, {2, 4}], 10^-10],
            within[QuantumEntanglementMonotone[qs, Automatic, "LogNegativity"], Log2[2 negRef[N @ rhoNDM, {2, 4}] + 1], 10^-10],
            within[QuantumEntanglementMonotone[qs, Automatic, "Realignment"], ccnrRef[N @ rhoNDM, {2, 4}], 10^-10],
            within[QuantumEntanglementMonotone[qs, Automatic, "Concurrence"], QuantumEntanglementMonotone[ref, Automatic, "Concurrence"], 10^-10],
            QuantumEntangledQ[qs] === QuantumEntangledQ[ref] === True
        }
    ],
    {True, True, True, True, True},
    {___},
    TestID -> "A8-NonOrthonormalBasis-ReadsTheState"
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


(* mixed and pure states, exact and machine, stored in a product of local orthonormal frames, give the
   plain-WL values of the state itself *)
VerificationTest[
    SeedRandom[41];
    With[{cases = {
            {N @ werner[2/3], {2, 2}}, {werner[2/3], {2, 2}}, {randomMixed[6], {2, 3}}, {randomMixed[9], {3, 3}},
            {proj[randomKet[6]], {2, 3}}, {proj[Normalize[{1, 2, 0, 3}]], {2, 2}}
        }},
        Max @ Flatten @ Map[
            Function[c, With[{rho = c[[1]], dims = c[[2]]},
                Table[
                    Abs[N @ methodValues[inLocalFrame[rho, dims, frame]] - N @ {negRef[rho, dims], ccnrRef[rho, dims]}],
                    {frame, Join[namedFrames[dims], randomFrames[dims, 3]]}
                ]
            ]],
            cases
        ] < 10^-9
    ],
    True,
    {___},
    TestID -> "A15-OrthonormalFrames-SameValues"
]

(* kets in orthonormal frames, whole and split into a pair after a partial trace *)
VerificationTest[
    SeedRandom[43];
    With[{v = randomKet[8], frames = Join[namedFrames[{2, 2, 2}], randomFrames[{2, 2, 2}, 3]]},
        With[{pair = trB[proj[v], {4, 2}]},
            Max @ Flatten @ Table[
                With[{qs = inLocalFrame[v, {2, 2, 2}, frame]},
                    {
                        Abs[N @ QuantumEntanglementMonotone[qs, {{1}, {2}}, "Negativity"] - negRef[pair, {2, 2}]],
                        Abs[N @ QuantumEntanglementMonotone[qs, Automatic, "Negativity"] - negRef[proj[v], {2, 4}]],
                        Abs[N @ QuantumEntanglementMonotone[qs, Automatic, "Concurrence"] - pureConcRef[v, {2, 4}]]
                    }
                ],
                {frame, frames}
            ] < 10^-9
        ]
    ],
    True,
    {___},
    TestID -> "A16-OrthonormalFrames-KetAndItsPairs"
]

(* the same in frames that are not orthonormal: {|0>, |+>} per qubit and a skew qutrit frame. The ket with
   coefficients (1, 0, 0, 1)/Sqrt[2] in {|0>, |+>} (x) {|0>, |+>} is (|00> + |++>) normalized, with
   negativity 1/6 *)
VerificationTest[
    With[{qs = QuantumState[Normalize[{1, 0, 0, 1}], skewLocalFrame[{2, 2}]]},
        {QuantumEntanglementMonotone[qs, "Negativity"], QuantumEntanglementMonotone[qs, "Realignment"], QuantumEntanglementMonotone[qs, "Concurrence"]}
    ],
    {1/6, 1/3, 1/3},
    {___},
    SameTest -> (Max[Abs[N[#1 - #2]]] < 10^-12 &),
    TestID -> "A17-NonOrthonormalFrame-Ket"
]

VerificationTest[
    SeedRandom[47];
    With[{v = randomKet[8], w = randomKet[12]},
        Max @ {
            Abs[N @ QuantumEntanglementMonotone[QuantumState[QuantumState[v, {2, 2, 2}], skewLocalFrame[{2, 2, 2}]], {{1}, {2}}, "Negativity"] - negRef[trB[proj[v], {4, 2}], {2, 2}]],
            Abs[N @ QuantumEntanglementMonotone[QuantumState[QuantumState[v, {2, 2, 2}], skewLocalFrame[{2, 2, 2}]], {{1}, {3}}, "Concurrence"] -
                With[{r = TensorContract[ArrayReshape[proj[v], {2, 2, 2, 2, 2, 2}], {{2, 5}}]}, QuantumEntanglementMonotone[QuantumState[ArrayReshape[r, {4, 4}], {2, 2}], "Concurrence"]]],
            Abs[N @ QuantumEntanglementMonotone[QuantumState[QuantumState[w, {2, 3, 2}], skewLocalFrame[{2, 3, 2}]], {{1}, {2}}, "Negativity"] - negRef[trB[proj[w], {6, 2}], {2, 3}]],
            Abs[N @ QuantumEntanglementMonotone[QuantumState[QuantumState[w, {2, 3, 2}], skewLocalFrame[{2, 3, 2}]], {1, 3}, "Negativity"] -
                negRef[ArrayReshape[Transpose[ArrayReshape[proj[w], {2, 3, 2, 2, 3, 2}], {1, 3, 2, 4, 6, 5}], {12, 12}], {4, 3}]]
        } < 10^-9
    ],
    True,
    {___},
    TestID -> "A18-NonOrthonormalFrame-PairsOfKets"
]

(* a density matrix changed into a frame that is not orthonormal, the content of A8: the stored matrix is
   S^-1 rho S, and its computational form is rho itself *)
VerificationTest[
    SeedRandom[53];
    With[{rho = randomMixed[8]},
        With[{qs = QuantumState[QuantumState[rho, {2, 2, 2}], skewLocalFrame[{2, 2, 2}]]},
            {
                within[QuantumEntanglementMonotone[qs, Automatic, "Negativity"], negRef[rho, {2, 4}], 10^-9],
                within[QuantumEntanglementMonotone[qs, {{1}, {3}}, "Negativity"], negRef[ArrayReshape[TensorContract[ArrayReshape[rho, ConstantArray[2, 6]], {{2, 5}}], {4, 4}], {2, 2}], 10^-9],
                within[QuantumEntanglementMonotone[qs, Automatic, "Realignment"], ccnrRef[rho, {2, 4}], 10^-9]
            }
        ]
    ],
    {True, True, True},
    {___},
    TestID -> "A19-NonOrthonormalFrame-ChangedDensityMatrix"
]

(* a state held in the phase-space picture has no density matrix in the computational basis: refused *)
VerificationTest[
    With[{states = {QuantumPhaseSpaceTransform[QuantumState["PhiPlus"]], QuantumPhaseSpaceTransform[QuantumState[werner[2/3], {2, 2}]]}},
        Flatten @ Map[{FailureQ[QuantumEntanglementMonotone[#, "Negativity"]], FailureQ[QuantumEntanglementMonotone[#, "Realignment"]], FailureQ[QuantumEntangledQ[#]]} &, states]
    ],
    ConstantArray[True, 6],
    {__},
    TestID -> "A20-PhaseSpacePicture-Refused"
]

(* a mixed state stored in a frame whose elements are dependent has no density matrix in the computational
   basis either *)
VerificationTest[
    With[{qs = QuantumState[werner[2/3], QuantumTensorProduct[QuantumBasis[{{1, 0}, {1, 0}}], QuantumBasis[2]]]},
        {FailureQ[QuantumEntanglementMonotone[qs, "Negativity"]], FailureQ[QuantumEntangledQ[qs]]}
    ],
    {True, True},
    {__},
    TestID -> "A21-MixedStateInDependentFrame-Refused"
]


(* a mixed state gives no state other than itself: the property refuses it, or returns the same state *)
VerificationTest[
    SeedRandom[59];
    With[{w3m = Normal @ QuantumState["W"[3]]["DensityMatrix"], g3m = Normal @ QuantumState["GHZ"[3]]["DensityMatrix"]},
        With[{cases = {(w3m + g3m) / 2, (w3m + 2 g3m) / 3, N[(w3m + g3m) / 2], randomMixed[8]}},
            Flatten @ Table[
                With[{s = TimeConstrained[inLocalFrame[rho, {2, 2, 2}, frame]["SchmidtBasis", 2], 30, $TimedOut]},
                    ! MatchQ[s, _QuantumState | $TimedOut] ||
                        MatchQ[s, _QuantumState] && Max[Abs[N @ Flatten[Normal[s["Computational"]["DensityMatrix"]] / Tr[Normal[s["Computational"]["DensityMatrix"]]] - rho]]] < 10^-10
                ],
                {rho, cases}, {frame, {QuantumBasis[{2, 2, 2}], QuantumBasis["PauliX", 3], QuantumBasis["PauliY", 3]}}
            ]
        ]
    ],
    ConstantArray[True, 12],
    {___},
    TestID -> "A22-MixedSchmidtBasis-SameStateOrRefused"
]

(* a pure state held as a density matrix has the Schmidt basis of its state vector: the same state, with at
   most min(dA, dB) nonzero amplitudes *)
VerificationTest[
    SeedRandom[61];
    With[{cases = {{Normalize[{1, 2, 0, 3}], {2, 2}}, {randomKet[6], {2, 3}}, {Normalize[{1, 0, 0, 0, 1, 0, 0, 0, 1}], {3, 3}}}},
        Map[
            Function[c, With[{s = QuantumState[proj[c[[1]]], c[[2]]]["SchmidtBasis", First[c[[2]]]]},
                MatchQ[s, _QuantumState] &&
                    Max[Abs[N @ Flatten[Normal[s["Computational"]["DensityMatrix"]] - proj[c[[1]]]]]] < 10^-10 &&
                    Count[Chop[N @ Normal @ s["StateVector"]], Except[0 | 0.]] <= Min[c[[2]]]
            ]],
            cases
        ]
    ],
    {True, True, True},
    TestID -> "A23-PureDensityMatrixSchmidtBasis"
]

(* a pair split traces out the qudits outside the pair, which needs the density matrix in the computational basis:
   a state held in the phase-space picture, or a density matrix in a frame whose elements are dependent, is refused
   for a pair split as it is for a split of the whole state *)
VerificationTest[
    With[{states = {
        QuantumPhaseSpaceTransform[w3],
        QuantumState[KroneckerProduct[proj[{1, 0}], proj[{1, 0, 0, 1} / Sqrt[2]]], QuantumTensorProduct[QuantumBasis[{{1, 0}, {1, 0}}], QuantumBasis[2], QuantumBasis[2]]]
    }},
        Flatten @ Map[
            {
                FailureQ[QuantumEntanglementMonotone[#, {{1}, {2}}, "Negativity"]],
                FailureQ[QuantumEntanglementMonotone[#, {{2}, {3}}, "Concurrence"]],
                FailureQ[QuantumEntangledQ[#, {{2}, {3}}]]
            } &,
            states
        ]
    ],
    ConstantArray[True, 6],
    {__},
    TestID -> "A24-PairSplitWithoutDensityMatrix-Refused"
]

(* the same for a flat list of qudits, which permutes the qudits before the split *)
VerificationTest[
    With[{bellPS = QuantumPhaseSpaceTransform[bell], w3PS = QuantumPhaseSpaceTransform[w3]},
        Flatten @ Map[
            {FailureQ[QuantumEntanglementMonotone[#[[1]], #[[2]], "Negativity"]], FailureQ[QuantumEntangledQ[#[[1]], #[[2]]]]} &,
            {{bellPS, {1}}, {bellPS, {2}}, {w3PS, {2}}, {w3PS, {1, 3}}}
        ]
    ],
    ConstantArray[True, 8],
    {__},
    TestID -> "A25-FlatListSplitWithoutDensityMatrix-Refused"
]

(* a pure density matrix of any trace has the Schmidt basis of the same state, exact and machine *)
VerificationTest[
    With[{rho = proj[{1, 2, 0, 3}]},
        {
            Simplify[Normal @ QuantumState[rho, {2, 2}]["SchmidtBasis"]["Computational"]["DensityMatrix"] - rho] === ConstantArray[0, {4, 4}],
            Max[Abs[Normal @ QuantumState[N @ rho, {2, 2}]["SchmidtBasis"]["Computational"]["DensityMatrix"] - rho]] < 10^-12
        }
    ],
    {True, True},
    TestID -> "A26-UnnormalizedPureDensityMatrixSchmidtBasis"
]

(* an exact pure density matrix whose purity, read in another frame, is a sum of radicals equal to 1 *)
VerificationTest[
    With[{rho = 3 proj[Normalize[{1, Sqrt[2], I, 1 + Sqrt[3]}]]},
        Max[Abs[N[Normal @ QuantumState[QuantumState[rho, {2, 2}], QuantumBasis[#, 2]]["SchmidtBasis", 2]["Computational"]["DensityMatrix"] - rho]]] < 10^-12 & /@ {"PauliX", "Fourier"}
    ],
    {True, True},
    TestID -> "A27-ExactPureDensityMatrixInOtherFrameSchmidtBasis"
]

(* the concurrence vector refuses a state with no density matrix in the computational basis, as every method does *)
VerificationTest[
    FailureQ[QuantumEntanglementMonotone[#, Automatic, "ConcurrenceVector"]] & /@ {
        QuantumPhaseSpaceTransform[bell],
        QuantumState[werner[2/3], QuantumTensorProduct[QuantumBasis[{{1, 0}, {1, 0}}], QuantumBasis[2]]]
    },
    {True, True},
    {__},
    TestID -> "A28-ConcurrenceVectorWithoutDensityMatrix-Failure"
]

(* a symbolic density matrix has no state vector valid for all its parameter values: no Schmidt basis *)
VerificationTest[
    Block[{t}, MatchQ[QuantumState[proj[{Cos[t], 0, 0, Sin[t]}], {2, 2}]["SchmidtBasis"], _QuantumState]],
    False,
    {QuantumState::failprop},
    TestID -> "A29-SymbolicDensityMatrixSchmidtBasis-Refused"
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

(* rank 2, at most the larger of its marginal ranks (2 and 2), with a positive partial transpose: separable *)
VerificationTest[
    QuantumEntangledQ[QuantumState[N @ exactSep33, {3, 3}]],
    False,
    {___},
    TestID -> "D8-Separable3x3AbovePurityBound-LowRank"
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


(* states that are entangled with a positive partial transpose never read as separable, exact or machine,
   in any local frame, and in larger local spaces *)
VerificationTest[
    SeedRandom[67];
    With[{cases = Join[
            {horodecki33Family[#], {3, 3}} & /@ {1/4, 1/2, 3/4},
            {horodecki24Family[#], {2, 4}} & /@ {1/4, 1/2, 3/4},
            {breuerState[#, 4], {4, 4}} & /@ {1/20, 1/10, 3/20},
            {{tilesState, {3, 3}}, {embedState[horodecki33Family[1/2], {3, 3}, {3, 8}], {3, 8}}, {embedState[horodecki24Family[1/2], {2, 4}, {2, 8}], {2, 8}}}
        ]},
        Flatten @ Map[
            Function[c, With[{rho = c[[1]], dims = c[[2]]},
                MatchQ[TimeConstrained[QuantumEntangledQ[#], 30, $TimedOut], True | Indeterminate] & /@ {
                    QuantumState[N @ rho, dims], QuantumState[rho, dims],
                    QuantumState[N @ localUnitary[rho, haarU /@ dims], dims]
                }
            ]],
            cases
        ]
    ],
    ConstantArray[True, 36],
    {___},
    TestID -> "D12-PPTEntangledFamily-NeverFalse"
]

(* separable states inside the purity bound Tr rho^2 <= 1/(dA dB - 1), exact and machine, read as separable:
   (1 - e) I/d + e sigma with sigma a product state has purity (1 - e^2)/d + e^2 Tr sigma^2, within the bound
   for e <= 1/(d - 1) *)
VerificationTest[
    SeedRandom[71];
    Flatten @ Table[
        With[{d = Times @@ dims},
            With[{
                exact = (1 - 1/d) IdentityMatrix[d] / d + (1/d) exactProduct[dims],
                machine = With[{s = productMixture[dims, {1/2, 1/2}]}, N[1 - 1/d] IdentityMatrix[d] / d + N[1/d] s]
            },
                {
                    TimeConstrained[QuantumEntangledQ[QuantumState[exact, dims]], 30, $TimedOut],
                    QuantumEntangledQ[QuantumState[machine, dims]]
                }
            ]
        ],
        {dims, {{3, 3}, {2, 4}, {3, 4}, {4, 4}}}
    ],
    ConstantArray[False, 8],
    TestID -> "D13-SeparableInsidePurityBound-False"
]

(* separable states that no separability proof covers (beyond 2 x 3, full rank, purity above the bound)
   are not read as separable: three pure product states plus a tenth of a full-rank product state *)
VerificationTest[
    SeedRandom[73];
    Table[
        With[{d = Times @@ dims},
            With[{rho = 9/10 Total[Table[KroneckerProduct[proj[randomKet[First[dims]]], proj[randomKet[Last[dims]]]], 3]] / 3 +
                    1/10 KroneckerProduct[randomMixed[First[dims]], randomMixed[Last[dims]]]},
                QuantumEntangledQ[QuantumState[rho, dims]]
            ]
        ],
        {dims, {{3, 3}, {2, 4}, {3, 4}}}
    ],
    {Indeterminate, Indeterminate, Indeterminate},
    {___},
    TestID -> "D14-SeparableWithoutProof-NotFalse"
]

(* exact separable mixtures of product states are decided, never left unevaluated, and without messages *)
VerificationTest[
    SeedRandom[79];
    Table[
        With[{rho = (exactProduct[dims] + 2 exactProduct[dims]) / 3},
            TimeConstrained[QuantumEntangledQ[QuantumState[rho, dims]], 30, $TimedOut]
        ],
        {dims, {{2, 2}, {2, 3}, {3, 2}}}
    ],
    {False, False, False},
    {},
    TestID -> "D15-ExactSeparableSmallDimensions-False"
]

VerificationTest[
    SeedRandom[83];
    Table[
        With[{rho = (exactProduct[dims] + 2 exactProduct[dims]) / 3},
            MatchQ[TimeConstrained[QuantumEntangledQ[QuantumState[rho, dims]], 30, $TimedOut], False | Indeterminate]
        ],
        {dims, {{3, 3}, {2, 4}, {4, 4}}}
    ],
    {True, True, True},
    {},
    TestID -> "D16-ExactSeparableLargerDimensions-Decided"
]

(* exact states entangled by a tiny margin are detected *)
VerificationTest[
    {
        QuantumEntangledQ[QuantumState[werner[1/3 + 10^-3], {2, 2}]],
        QuantumEntangledQ[QuantumState[werner[1/3 + 10^-20], {2, 2}]],
        QuantumEntangledQ[QuantumState[With[{phi = Flatten[IdentityMatrix[3]] / Sqrt[3]}, (1/4 + 10^-15) proj[phi] + (3/4 - 10^-15) IdentityMatrix[9] / 9], {3, 3}]],
        QuantumEntangledQ[QuantumState[werner[1/3 - 10^-20], {2, 2}]]
    },
    {True, True, True, False},
    {},
    TestID -> "D17-ExactNearThreshold"
]

(* entangled states with a negative partial transpose read as entangled, mixed and pure *)
VerificationTest[
    SeedRandom[89];
    Flatten @ Table[
        With[{v = randomKet[Times @@ dims]},
            With[{rho = 9/10 proj[v] + 1/10 IdentityMatrix[Times @@ dims] / Times @@ dims},
                {QuantumEntangledQ[QuantumState[rho, dims]], QuantumEntangledQ[QuantumState[v, dims]]}
            ]
        ],
        {dims, {{2, 2}, {2, 3}, {3, 3}, {2, 4}, {4, 4}}}
    ],
    ConstantArray[True, 10],
    TestID -> "D18-NPTFamily-True"
]

(* Tiles in 3 x 4 plus a weight e on each |i><i| (x) |4><4|: positive partial transpose, rank 7, marginal ranks 3
   and 4, and entangled, since projecting B onto its first three levels gives Tiles back. Eigenvalues of size e, near
   round-off, never prove it separable: one dropped from rho while kept in a marginal, where three of them add, would
   make the low-rank proof *)
VerificationTest[
    With[{noisy = Function[e, (1 - 3 e) embedState[tilesState, {3, 3}, {3, 4}] + e Total[Table[KroneckerProduct[proj[UnitVector[3, i]], proj[UnitVector[4, 4]]], {i, 3}]]]},
        Flatten @ Table[
            MatchQ[QuantumEntangledQ[QuantumState[#, {3, 4}]], True | Indeterminate] & /@ {N @ noisy[e], noisy[e], SetPrecision[noisy[e], 30]},
            {e, {10^-6, 10^-9, 5 10^-11, 3 10^-11, 10^-13}}
        ]
    ],
    ConstantArray[True, 15],
    {___},
    TestID -> "D19-PPTEntangledNearRoundoff-NeverFalse"
]

(* Tiles and the Horodecki 3 x 3 state with weight e in the first three levels of 4 x 4 and 1 - e on |44>: positive
   partial transpose, entangled, and within e of a pure product state. Given with p digits, where e is far above
   round-off, it is never read as pure or of low rank, and the realignment proves it entangled once its value is
   above the 10^-10 cut *)
VerificationTest[
    With[{corner = Function[{sigma, e}, e embedState[sigma, {3, 3}, {4, 4}] + (1 - e) proj[UnitVector[16, 16]]]},
        {
            Flatten @ Table[
                MatchQ[QuantumEntangledQ[QuantumState[SetPrecision[corner[sigma, pe[[2]]], pe[[1]]], {4, 4}]], True | Indeterminate],
                {sigma, {tilesState, horodecki33Family[1/2]}},
                {pe, {{12, 10^-10}, {10, 10^-8}, {8, 10^-6}, {6, 10^-4}, {5, 10^-3}}}
            ],
            QuantumEntangledQ[QuantumState[SetPrecision[corner[tilesState, 10^-4], 6], {4, 4}], Automatic, #] & /@ {Automatic, "Realignment"}
        }
    ],
    {ConstantArray[True, 10], {True, True}},
    {___},
    TestID -> "D20-PPTEntangledAtLowPrecision-NeverFalse"
]

(* the zero ket has no entanglement to read: never True, by Automatic or by the concurrence *)
VerificationTest[
    Flatten @ Table[(QuantumEntangledQ[QuantumState[ConstantArray[0, Times @@ dims], dims], Automatic, m] =!= True), {dims, {{2, 2}, {3, 3}}}, {m, {Automatic, "Concurrence"}}],
    {True, True, True, True},
    {___},
    TestID -> "D21-ZeroKet-NeverTrue"
]

(* the pure-state cut on approximate input follows the precision: a 30-digit state 5 10^-27 away from |00><00|,
   an impurity above the cut n 10^(2 - p) = 1.6 10^-27, is not read as pure *)
VerificationTest[
    With[{phi = Flatten[IdentityMatrix[4]] / 2},
        QuantumEntangledQ[QuantumState[SetPrecision[(1 - 5 10^-27) proj[UnitVector[16, 1]] + 5 10^-27 proj[phi], 30], {4, 4}]] =!= False
    ],
    True,
    {___},
    TestID -> "D22-ApproximatePureCutFollowsPrecision"
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


(* symbolic input, real parameters over the values for which the input is a density matrix: a family inside
   the purity bound for every real p is separable, and the isotropic family changes at p = 1/4 *)
VerificationTest[
    Block[{p},
        {
            QuantumEntangledQ[QuantumState[IdentityMatrix[9] / 9 + p / (36 (1 + p^2)) KroneckerProduct[DiagonalMatrix[{1, -1, 0}], DiagonalMatrix[{1, 0, -1}]], {3, 3}]],
            QuantumEntangledQ[QuantumState[With[{phi = Flatten[IdentityMatrix[3]] / Sqrt[3]}, p proj[phi] + (1 - p) IdentityMatrix[9] / 9], {3, 3}]]
        }
    ],
    {False, Indeterminate},
    {__},
    TestID -> "E19-SymbolicFamiliesBeyond2x3"
]

(* the parameters of a symbolic state are its outermost non-numeric parts, indexed ones included *)
VerificationTest[
    Block[{x, p},
        {
            QuantumEntangledQ[QuantumState[KroneckerProduct[DiagonalMatrix[{x[1], 1 - x[1]}], DiagonalMatrix[{x[2], 1 - x[2]}]], {2, 2}]],
            QuantumEntangledQ[QuantumState[KroneckerProduct[DiagonalMatrix[{Subscript[p, 1], 1 - Subscript[p, 1]}], DiagonalMatrix[{Subscript[p, 2], 1 - Subscript[p, 2]}]], {2, 2}]]
        }
    ],
    {False, False},
    {QuantumEntangledQ::realparameters, QuantumEntangledQ::realparameters},
    TestID -> "E20-SymbolicIndexedParameters"
]

(* machine coefficients, whose zero imaginary parts 0. are no nonreal constants, and a family whose realignment value
   has none in closed form: decided or Indeterminate, with no message from Resolve *)
VerificationTest[
    Block[{a, p, q},
        {
            QuantumEntangledQ[QuantumState[N @ KroneckerProduct[DiagonalMatrix[{a, 1 - a}], DiagonalMatrix[{1/3, 2/3}]], {2, 2}]],
            With[{r1 = KroneckerProduct[proj[{1, 0, 0}], proj[{0, 1, 0}]], r2 = proj[Flatten[IdentityMatrix[3]] / Sqrt[3]]},
                QuantumEntangledQ[QuantumState[p r2 + q r1 + (1 - p - q) IdentityMatrix[9] / 9, {3, 3}]]
            ]
        }
    ],
    {False, Indeterminate},
    {QuantumEntangledQ::realparameters, QuantumEntangledQ::realparameters},
    TestID -> "E21-SymbolicNonrealConstants-NoResolveMessage"
]

(* a dense three-parameter family on 3 x 3: every symbolic step has a time limit *)
VerificationTest[
    SeedRandom[21];
    With[{r = Table[With[{m = RandomInteger[{-3, 3}, {9, 9}] + I RandomInteger[{-3, 3}, {9, 9}]}, # / Tr[#] & [m . ConjugateTranspose[m]]], 3]},
        Block[{a, b, c},
            With[{t = AbsoluteTiming[QuantumEntangledQ[QuantumState[a r[[1]] + b r[[2]] + c r[[3]] + (1 - a - b - c) IdentityMatrix[9] / 9, {3, 3}]]]},
                {Last[t], First[t] < 60}
            ]
        ]
    ],
    {Indeterminate, True},
    {QuantumEntangledQ::realparameters},
    TestID -> "E22-SymbolicDenseFamily-Bounded"
]

(* the exact Werner state 10^-46 away from p = 1/3, the margin a difference of transcendental numbers: decided on
   both sides with no message, the same in any order of evaluation *)
VerificationTest[
    With[{d = Abs[Pi - Rationalize[N[Pi, 60], 10^-45]]},
        QuantumEntangledQ[QuantumState[werner[#], {2, 2}]] & /@ {1/3 - d, 1/3 + d, 1/3 - d, 1/3 + d}
    ],
    {False, True, False, True},
    {},
    TestID -> "E23-ExactTranscendentalTinyMargin-Decided"
]

(* an algebraic margin of 10^-300, beyond what $MaxExtraPrecision reads by default: decided on both sides *)
VerificationTest[
    With[{d = Abs[Sqrt[2] - Rationalize[N[Sqrt[2], 400], 10^-300]]},
        QuantumEntangledQ[QuantumState[werner[#], {2, 2}]] & /@ {1/3 - d, 1/3 + d}
    ],
    {False, True},
    {},
    TestID -> "E24-ExactAlgebraicTinyMargin-Decided"
]

(* an exact ket held in the PauliX frame, whose coefficients carry an exact zero written as radicals: no N::meprec *)
VerificationTest[
    With[{qs = QuantumState[QuantumState[Normalize[{1, 2, 3, 4}], {2, 2}], QuantumBasis["PauliX", 2]]},
        QuantumEntangledQ[qs, Automatic, #] & /@ {Automatic, "Negativity", "Concurrence", "EntanglementEntropy"}
    ],
    {True, True, True, True},
    {},
    TestID -> "E25-ExactKetInPauliXFrame-NoMessage"
]

(* symbolic families with machine coefficients reach Resolve as exact numbers: no Resolve message *)
VerificationTest[
    Block[{x},
        {
            QuantumEntangledQ[QuantumState[N @ With[{phi = Flatten[IdentityMatrix[3]] / Sqrt[3]}, x proj[phi] + (1 - x) IdentityMatrix[9] / 9], {3, 3}]],
            QuantumEntangledQ[QuantumState[N @ werner[x], {2, 2}], Automatic, "Concurrence"]
        }
    ],
    {Indeterminate, Indeterminate},
    {QuantumEntangledQ::realparameters, QuantumEntangledQ::realparameters},
    TestID -> "E26-SymbolicMachineCoefficients-NoResolveMessage"
]

(* a machine family whose coherence c makes it entangled for every c != 0, with negativities far below and just
   above 10^-10: no coefficient is dropped as round-off, so it is never read as separable *)
VerificationTest[
    Block[{c},
        QuantumEntangledQ[QuantumState[(1 - e) proj[{1, 0, 0, 0}] + (e / 2) {{1, 0, 0, c}, {0, 0, 0, 0}, {0, 0, 0, 0}, {c, 0, 0, 1}}, {2, 2}]] /.
            e -> # & /@ {1.*^-5, 1.*^-9}
    ],
    {Indeterminate, Indeterminate},
    {QuantumEntangledQ::realparameters, QuantumEntangledQ::realparameters},
    TestID -> "E27-SymbolicMachineSmallCoherence-NotFalse"
]

(* algebraic margins of 10^-300 of high degree, sums of cube roots, fifth roots and square roots of primes: decided by
   the long readings, within the time limit *)
VerificationTest[
    With[{ds = Abs[# - Rationalize[N[#, 400], 10^-300]] & /@ {Total[Prime[Range[6]] ^ (1/3)], Total[Prime[Range[4]] ^ (1/5)], Total[Sqrt[Prime[Range[12]]]]}},
        With[{t = AbsoluteTiming[Flatten[{QuantumEntangledQ[QuantumState[werner[1/3 - #], {2, 2}]], QuantumEntangledQ[QuantumState[werner[1/3 + #], {2, 2}]]} & /@ ds]]},
            {Last[t], First[t] < 10}
        ]
    ],
    {{False, True, False, True, False, True}, True},
    {},
    TestID -> "E28-ExactHighDegreeAlgebraicMargins-Decided"
]

(* an exact ket held in the Fourier[3] (x) Fourier[3] frame: every method reads it without a message *)
VerificationTest[
    With[{qs = QuantumState[QuantumState[Normalize[{1, 2, 0, 0, 1, 3, 1, 0, 2}], {3, 3}], QuantumBasis[QuantumTensorProduct[QuditBasis["Fourier"[3]], QuditBasis["Fourier"[3]]]]]},
        With[{t = AbsoluteTiming[{
                QuantumEntanglementMonotone[qs, Automatic, "EntanglementEntropy"], QuantumEntanglementMonotone[qs, Automatic, "Negativity"],
                QuantumEntanglementMonotone[qs, Automatic, "Realignment"], QuantumEntangledQ[qs]
            }]},
            {Last[t][[4]], First[t] < 30}
        ]
    ],
    {True, True},
    {},
    TestID -> "E29-ExactKetInFourierFrame-NoMessage"
]

(* a symbolic ket held in the PauliX frame keeps its coefficients as they are, and its values come as fast *)
VerificationTest[
    Block[{a, b, c, d},
        With[{qs = QuantumState[QuantumState[{a, b, c, d}, {2, 2}], QuantumBasis["PauliX", 2]]},
            First[AbsoluteTiming[QuantumEntanglementMonotone[qs, Automatic, #]]] < 20 & /@ {"Negativity", "Realignment"}
        ]
    ],
    {True, True},
    TestID -> "E30-SymbolicKetInPauliXFrame-Bounded"
]

(* the concurrence of the machine zero ket is 0, with no division by its zero Schmidt weights *)
VerificationTest[
    QuantumEntanglementMonotone[QuantumState[{0., 0., 0., 0.}, {2, 2}], "Concurrence"],
    0,
    {QuantumEntanglementMonotone::notphysical},
    TestID -> "E31-MachineZeroKetConcurrence"
]

(* kets held in the Fourier[7] (x) Fourier[2] frame, whose coefficients in the computational frame are sums of 7th
   roots of unity: the product ket reads as a product, the Bell ket as maximally entangled, with no message *)
VerificationTest[
    With[{frame = QuantumBasis[QuantumTensorProduct[QuditBasis["Fourier"[7]], QuditBasis["Fourier"[2]]]]},
        With[{product = QuantumState[QuantumState[UnitVector[14, 1], {7, 2}], frame], bellLike = QuantumState[QuantumState[Normalize[UnitVector[14, 1] + UnitVector[14, 14]], {7, 2}], frame]},
            {
                {QuantumEntanglementMonotone[product, "Concurrence"], QuantumEntanglementMonotone[product, Automatic, "Negativity"], QuantumEntangledQ[product]},
                N @ {QuantumEntanglementMonotone[bellLike, Automatic, "Negativity"], QuantumEntanglementMonotone[bellLike, Automatic, "Realignment"], QuantumEntangledQ[bellLike]}
            }
        ]
    ],
    {{0, 0, False}, {0.5, 1., True}},
    {},
    TestID -> "E32-ExactKetsInFourier7Frame-NoMessage"
]

(* a pure product ket of cube roots, whose partial-transpose sign is an exact zero of high algebraic degree: read as
   0 by FullSimplify before RootReduce takes the time limit *)
VerificationTest[
    With[{t = AbsoluteTiming[QuantumEntangledQ[QuantumState[Flatten[KroneckerProduct[{1, 2^(1/3), 3^(1/3)}, {1, 5^(1/3), 7^(1/3)}]], {3, 3}]]]},
        {Last[t], First[t] < 10}
    ],
    {False, True},
    {},
    TestID -> "E33-ExactHighDegreeZero-Decided"
]

(* a symbolic family with machine coefficients far below machine underflow, entangled for every c != 0: its numbers
   are read exactly, so nothing underflows and it is never read as separable *)
VerificationTest[
    Block[{c},
        QuantumEntangledQ[QuantumState[DiagonalMatrix[{0, 1, 1, 0}] / 2 + # c SparseArray[{{2, 3} -> 1, {3, 2} -> 1}, {4, 4}], {2, 2}]] & /@ {1.*^-170, 1.*^-300}
    ],
    {Indeterminate, Indeterminate},
    {QuantumEntangledQ::realparameters, QuantumEntangledQ::realparameters},
    TestID -> "E34-SymbolicTinyMachineCoefficients-NoUnderflow"
]

(* an exact zero built from a constant whose 5000-digit value is slow, Khinchin's: FullSimplify reads it as 0 before
   that reading is tried *)
VerificationTest[
    With[{t = AbsoluteTiming[QuantumEntangledQ[QuantumState[{Cos[Khinchin], 0, Sin[Khinchin], 0}, {2, 2}]]]},
        {Last[t], First[t] < 5}
    ],
    {False, True},
    {},
    TestID -> "E35-ExactZeroOfSlowConstant-Decided"
]

(* the named realignment on an exact ket of cube roots, whose exact singular values take a minute: read within the
   time limit, never True on this product state *)
VerificationTest[
    With[{t = AbsoluteTiming[QuantumEntangledQ[QuantumState[Flatten[KroneckerProduct[{1, 2^(1/3), 3^(1/3)}, {1, 5^(1/3), 7^(1/3)}]], {3, 3}], Automatic, "Realignment"]]},
        {MatchQ[Last[t], False | Indeterminate], First[t] < 15}
    ],
    {True, True},
    {},
    TestID -> "E36-NamedRealignmentOnExactInput-Bounded"
]

(* a symbolic state with exact coefficients of 10^-400: no physicality test on it, so no underflow message *)
VerificationTest[
    Block[{c}, QuantumEntangledQ[QuantumState[DiagonalMatrix[{0, 1, 1, 0}] / 2 + 10^-400 c {{0, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 0}}, {2, 2}]]],
    Indeterminate,
    {QuantumEntangledQ::realparameters},
    TestID -> "E37-SymbolicTinyExactCoefficients-NoUnderflow"
]

(* the parameter message comes with every symbolic input, also when normalization or the partial trace removes the
   parameter from the two-part state *)
VerificationTest[
    Block[{a, t},
        {
            QuantumEntangledQ[QuantumState[(1 + a^2) proj[{1, 0, 0, 1} / Sqrt[2]], {2, 2}]],
            QuantumEntangledQ[QuantumState[Flatten[KroneckerProduct[{1, 0}, {1, 0}, {Cos[t], Sin[t]}]], {2, 2, 2}], {{1}, {2}}]
        }
    ],
    {True, False},
    {QuantumEntangledQ::realparameters, QuantumEntangledQ::realparameters},
    TestID -> "E38-ParameterMessageOnEverySymbolicInput"
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

(* the default is unchanged by local unitaries, mixed and pure, for several seeds *)
VerificationTest[
    Max @ Flatten @ Table[
        SeedRandom[seed];
        With[{rho = randomMixed[Times @@ dims], v = randomKet[Times @@ dims], us = haarU /@ dims},
            {
                Abs[QuantumEntanglementMonotone[QuantumState[rho, dims]] - QuantumEntanglementMonotone[QuantumState[localUnitary[rho, us], dims]]],
                Abs[QuantumEntanglementMonotone[QuantumState[v, dims]] - QuantumEntanglementMonotone[QuantumState[(KroneckerProduct @@ us) . v, dims]]]
            }
        ],
        {seed, {1, 2, 3, 4, 5}}, {dims, {{3, 3}, {2, 4}, {3, 4}}}
    ] < 10^-9,
    True,
    TestID -> "F16-DefaultMonotone-LocalUnitaryInvariant"
]

(* concurrence for a pure two-part state and for two qubits, negativity for a mixed two-part state beyond two
   qubits, decided on the state after the split and the partial trace *)
VerificationTest[
    SeedRandom[97];
    With[{v = randomKet[27], rho = randomMixed[9], w = randomKet[9]},
        Chop[{
            QuantumEntanglementMonotone[QuantumState["W"[3]], {{1}, {3}}] - QuantumEntanglementMonotone[QuantumState["W"[3]], {{1}, {3}}, "Concurrence"],
            QuantumEntanglementMonotone[QuantumState[w, {3, 3}]] - pureConcRef[w, {3, 3}],
            QuantumEntanglementMonotone[QuantumState[proj[w], {3, 3}]] - pureConcRef[w, {3, 3}],
            QuantumEntanglementMonotone[QuantumState[rho, {3, 3}]] - negRef[rho, {3, 3}],
            QuantumEntanglementMonotone[QuantumState[v, {3, 3, 3}], {{1}, {2}}] - negRef[ArrayReshape[TensorContract[ArrayReshape[proj[v], ConstantArray[3, 6]], {{3, 6}}], {9, 9}], {3, 3}]
        }, 10^-9]
    ],
    {0, 0, 0, 0, 0},
    TestID -> "F17-DefaultMonotone-DecidedOnTheTwoPartState"
]

(* the mutual informations and the discord are positive on separable states too, so they never prove entanglement *)
VerificationTest[
    (QuantumEntangledQ[classical, Automatic, #] =!= True) & /@ {"MutualInformationI", "MutualInformationJ", "Discord"},
    {True, True, True},
    {___},
    TestID -> "F18-CorrelationMeasuresNeverProveEntanglement"
]

(* the entanglement entropy decides a pure state held as a density matrix, as it does a ket *)
VerificationTest[
    QuantumEntangledQ[QuantumState[#, {2, 2}], Automatic, "EntanglementEntropy"] & /@ {proj[{1, 0, 0, 1} / Sqrt[2]], proj[{1, 1, 0, 0} / Sqrt[2]]},
    {True, False},
    TestID -> "F19-EntropyMethodDecidesPureDensityMatrices"
]

(* the correlation measures refuse a state with no density matrix in the computational basis, as every method does *)
VerificationTest[
    With[{states = {QuantumPhaseSpaceTransform[bell], QuantumState[werner[2/3], QuantumTensorProduct[QuantumBasis[{{1, 0}, {1, 0}}], QuantumBasis[2]]]}},
        Flatten @ Outer[
            FailureQ[TimeConstrained[QuantumEntanglementMonotone[#1, Automatic, #2], 30, $TimedOut]] &,
            states,
            {"MutualInformationI", "MutualInformationJ", "Discord"}
        ]
    ],
    ConstantArray[True, 6],
    {__},
    TestID -> "F20-CorrelationMeasuresWithoutDensityMatrix-Failure"
]


EndTestSection[]
