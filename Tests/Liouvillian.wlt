BeginTestSection["Liouvillian - Hamiltonian relationship"]

(* For purely-coherent dynamics (no Lindblad operators), the Hamiltonian
   super-operator and the Liouvillian super-operator are related by
   ℋ == I ℒ. *)

VerificationTest[
    Block[{d = {3, "X"[2]}, h, hh, ll},
        SeedRandom[42];
        h = QuantumOperator["RandomHermitian"[d]];
        hh = QuantumOperator["Hamiltonian"[h, {}, {}]];
        ll = QuantumOperator["Liouvillian"[h, {}, {}]];
        hh == I ll
    ],
    True,
    TestID -> "ℋ-eq-I-ℒ-no-decoherence"
];

(* The same structural identity holds with Lindblad terms present. *)
VerificationTest[
    Block[{d = 2, h, ls, γs, hh, ll},
        SeedRandom[7];
        h = QuantumOperator["RandomHermitian"[d]];
        ls = {QuantumOperator["RandomHermitian"[d]]};
        γs = {0.3};
        hh = QuantumOperator["Hamiltonian"[h, ls, γs]];
        ll = QuantumOperator["Liouvillian"[h, ls, γs]];
        hh == I ll
    ],
    True,
    TestID -> "ℋ-eq-I-ℒ-with-Lindblad"
];

EndTestSection[]


BeginTestSection["Liouvillian - evolution equivalence"]

(* All eight evolution forms must yield the same final density matrix at t=1.
   The tolerance is loose because some routes go through numerical ODE
   integration while others use direct matrix exponentiation. *)

cmpEvolve[ref_, other_] := Chop @ Norm @ Flatten[N @ ref["DensityMatrix"] - N @ other["DensityMatrix"]]

With[{tol = 1*^-4},
    Block[{d, n, h, ls, γs, ρ, hh, ll, fA, fB, fC, fD, fE, fF, fG, fH},
        SeedRandom[1];
        d = 2;
        n = 0;
        h  = QuantumOperator["RandomHermitian"[d]];
        ls = Table[QuantumOperator["RandomHermitian"[d]], n];
        γs = RandomReal[1, n];
        ρ  = QuantumState["RandomMixed", d];
        hh = QuantumOperator["Hamiltonian"[h, ls, γs]];
        ll = QuantumOperator["Liouvillian"[h, ls, γs]];

        fA = Chop @ QuantumEvolve[hh, ρ, {t, 0, 2}][1];
        fB = Chop @ QuantumEvolve[hh, ρ, t][1];
        fC = Chop[QuantumEvolve[hh, None, {t, 0, 2}][1][ρ]];
        fD = Chop[QuantumEvolve[hh, None, t][1][ρ]];
        fE = Chop[Exp[-I hh][ρ]];
        fF = Chop[Exp[ll][ρ]];
        fG = Chop @ QuantumEvolve[h, ls -> γs, ρ, {t, 0, 2}][1];
        fH = Chop[QuantumEvolve[h, ls -> γs, None, {t, 0, 2}][1][ρ]];

        VerificationTest[cmpEvolve[fA, fB] < tol, True, TestID -> "Evolve-state-symbolic-t"];
        VerificationTest[cmpEvolve[fA, fC] < tol, True, TestID -> "Evolve-op-then-state"];
        VerificationTest[cmpEvolve[fA, fD] < tol, True, TestID -> "Evolve-op-symbolic-then-state"];
        VerificationTest[cmpEvolve[fA, fE] < tol, True, TestID -> "Evolve-Exp-Iℋ"];
        VerificationTest[cmpEvolve[fA, fF] < tol, True, TestID -> "Evolve-Exp-ℒ"];
        VerificationTest[cmpEvolve[fA, fG] < tol, True, TestID -> "Evolve-bare-Lindblad-form"];
        VerificationTest[cmpEvolve[fA, fH] < tol, True, TestID -> "Evolve-bare-Lindblad-op-then-state"];
    ]
]

EndTestSection[]


BeginTestSection["Liouvillian - Kossakowski matrix rates"]

(* A symbolic rate matrix must build: the zero-rate skip in the matrix branch
   applies only to provably zero entries, so undecidable comparisons like
   γ == 0 fall through to the general term instead of leaving an If behind. *)
VerificationTest[
    Block[{σm, id, l1, l2, h, ll},
        σm = QuantumOperator[{{0, 1}, {0, 0}}];
        id = QuantumOperator["I"];
        l1 = QuantumTensorProduct[σm, id];
        l2 = QuantumTensorProduct[id, σm];
        h = QuantumOperator["ZI"] + QuantumOperator["IZ"];
        ll = QuantumOperator["Liouvillian"[h, {l1, l2}, {{\[FormalG], \[FormalC]}, {\[FormalC], \[FormalG]}}]];
        MatchQ[ll, _QuantumOperator] && FreeQ[ll["Matrix"], If]
    ],
    True,
    TestID -> "Kossakowski-symbolic-rates-build"
];

(* QuantumEvolve accepts a rate matrix; a diagonal one must reproduce the
   vector-rate evolution. *)
VerificationTest[
    Block[{σm, id, l1, l2, h, ρ, fVec, fMat},
        SeedRandom[3];
        σm = QuantumOperator[{{0, 1}, {0, 0}}];
        id = QuantumOperator["I"];
        l1 = QuantumTensorProduct[σm, id];
        l2 = QuantumTensorProduct[id, σm];
        h = 0.35 (QuantumOperator["ZI"] + QuantumOperator["IZ"]);
        ρ = QuantumState["RandomMixed", {2, 2}];
        fVec = QuantumEvolve[h, {l1, l2} -> {0.8, 0.5}, ρ, {t, 0, 1}][1];
        fMat = QuantumEvolve[h, {l1, l2} -> {{0.8, 0.}, {0., 0.5}}, ρ, {t, 0, 1}][1];
        Norm[Flatten[N[fVec["DensityMatrix"]] - N[fMat["DensityMatrix"]]]] < 1*^-4
    ],
    True,
    TestID -> "Kossakowski-diagonal-equals-vector-rates"
];

(* Full off-diagonal rate matrix against a hand-written master equation
   dρ/dt = -i[H,ρ] + Σ_jk β_jk (L_j ρ L_k† - {L_k† L_j, ρ}/2). *)
VerificationTest[
    Block[{σm, id, l1, l2, lm, h, hm, β, ρ, f, rhs, ref},
        SeedRandom[3];
        σm = QuantumOperator[{{0, 1}, {0, 0}}];
        id = QuantumOperator["I"];
        l1 = QuantumTensorProduct[σm, id];
        l2 = QuantumTensorProduct[id, σm];
        lm = Normal /@ {l1["Matrix"], l2["Matrix"]};
        h = 0.35 (QuantumOperator["ZI"] + QuantumOperator["IZ"]);
        hm = Normal[h["Matrix"]];
        β = {{0.8, 0.3}, {0.3, 0.5}};
        ρ = QuantumState["RandomMixed", {2, 2}];
        f = QuantumEvolve[h, {l1, l2} -> β, ρ, {t, 0, 1}][1];
        rhs[m_] := -I (hm . m - m . hm) + Sum[
            β[[j, k]] (lm[[j]] . m . ConjugateTranspose[lm[[k]]] -
                (ConjugateTranspose[lm[[k]]] . lm[[j]] . m + m . ConjugateTranspose[lm[[k]]] . lm[[j]]) / 2),
            {j, 2}, {k, 2}];
        ref = NDSolveValue[{\[FormalR]'[t] == rhs[\[FormalR][t]], \[FormalR][0] == Normal[ρ["DensityMatrix"]]}, \[FormalR][1], {t, 0, 1}];
        Norm[Flatten[N[f["Computational"]["DensityMatrix"]] - ref]] < 1*^-4
    ],
    True,
    TestID -> "Kossakowski-matrix-rates-vs-master-equation"
];

EndTestSection[]


BeginTestSection["QuantumEvolve - DSolve options"]

(* DSolve-only options (notably Assumptions) must be forwarded to DSolveValue
   on the symbolic branch instead of being filtered against NDSolveValue. *)
VerificationTest[
    Block[{captured = $Failed, dsvOpts = Options[DSolveValue]},
        Block[{DSolveValue},
            Options[DSolveValue] = dsvOpts;
            DSolveValue[_, ret_, _, o___] := (captured = Flatten[{o}]; ret);
            QuantumEvolve[QuantumOperator["Z"], {QuantumOperator[{{0, 1}, {0, 0}}]} -> {\[FormalG]},
                QuantumState["1"], t, "ReturnSolution" -> True, Assumptions -> \[FormalG] > 0]
        ];
        MemberQ[captured, Assumptions -> \[FormalG] > 0]
    ],
    True,
    TestID -> "Evolve-Assumptions-reaches-DSolveValue"
];

EndTestSection[]


BeginTestSection["Liouvillian - phase-space (Wootters)"]

With[{d = 2},
    Block[{h, ll, mat, sf, tf, sopS, sopT, daggerLeft, daggerRight, denom},
        SeedRandom[1];
        h = QuantumOperator["RandomHermitian"[d]];
        ll = QuantumOperator["Liouvillian"[h, {}, {}]];

        (* Dagger commutes with the Wootters phase-space transform of Exp[ℒ]. *)
        daggerLeft = QuantumPhaseSpaceTransform[Exp[ll]["Dagger"], "Wootters"[d]]["StateVector"];
        daggerRight = QuantumPhaseSpaceTransform[Exp[ll], "Wootters"[d]]["Dagger"]["StateVector"];
        VerificationTest[
            Chop @ Norm @ Flatten[daggerLeft - daggerRight],
            0,
            SameTest -> (Chop[#1 - #2] === 0 &),
            TestID -> "Phase-space-Dagger-commutes"
        ];

        mat = Chop[Exp[QuantumPhaseSpaceTransform[ll, "Wootters"[d]]]]["Matrix"];
        denom = 1 - d^2 Min[mat];
        sf = (mat - Min[mat]) / denom;
        tf = (1 / denom) IdentityMatrix[d^2] + (1 - 1 / denom) (1 / d^2);

        (* Doubly-stochastic-style decomposition: sf . Inverse[tf] reproduces mat. *)
        VerificationTest[
            Chop @ Norm @ Flatten[sf . Inverse[tf] - mat],
            0,
            SameTest -> (Chop[#1 - #2] === 0 &),
            TestID -> "S-InvT-equals-M"
        ];

        (* The factors commute at d = 2. *)
        VerificationTest[
            Chop @ Norm @ Flatten[sf . Inverse[tf] - Inverse[tf] . sf],
            0,
            SameTest -> (Chop[#1 - #2] === 0 &),
            TestID -> "S-InvT-commute-d2"
        ];

        (* Operators reconstructed from sf and Inverse[tf] compose to Exp[QPST[ℒ]]. *)
        sopS = QuantumOperator[QuantumState[Flatten @ sf, QuantumPhaseSpaceTransform[ll, "Wootters"[d]]["Basis"]]];
        sopT = QuantumOperator[QuantumState[Flatten @ Inverse[tf], QuantumPhaseSpaceTransform[ll, "Wootters"[d]]["Basis"]]];
        VerificationTest[
            Chop[sopS[sopT] == Exp[QuantumPhaseSpaceTransform[ll, "Wootters"[d]]]],
            True,
            TestID -> "S-T-compose-equals-Exp-QPST-ℒ"
        ];
    ]
]

EndTestSection[]
