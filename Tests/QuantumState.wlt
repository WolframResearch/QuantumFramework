BeginTestSection["QuantumState - constructors"]

VerificationTest[QuantumState[]["Dimension"], 2, TestID -> "Empty-default-0"]

VerificationTest[QuantumState["Bell"]["Dimension"], 4, TestID -> "Bell-bare"]

VerificationTest[QuantumState["GHZ"]["Dimension"], 8, TestID -> "GHZ-default"]

VerificationTest[QuantumState["GHZ"[4]]["Dimension"], 16, TestID -> "GHZ-4"]

VerificationTest[QuantumState["W"[3]]["Dimension"], 8, TestID -> "W-3"]

VerificationTest[QuantumState["Werner"]["Dimension"], 4, TestID -> "Werner-default"]

VerificationTest[QuantumState["Werner"[1/3, 2]]["Dimension"], 4, TestID -> "Werner-call"]

VerificationTest[QuantumState["Dicke"[4, 2]]["Dimension"], 16, TestID -> "Dicke-4-2"]

VerificationTest[QuantumState["BlochVector"[{0, 1, 0}]]["Dimension"], 2, TestID -> "BlochVector"]

VerificationTest[QuantumState["UniformSuperposition"[3]]["Dimension"], 8, TestID -> "UniformSuperposition-3"]

VerificationTest[QuantumState["UniformMixture"[2]]["Dimension"], 4, TestID -> "UniformMixture-2"]

VerificationTest[QuantumState["Plus"]["Dimension"], 2, TestID -> "Plus-bare"]

VerificationTest[QuantumState["+"]["Dimension"], 2, TestID -> "+-shorthand"]

VerificationTest[QuantumState["101"]["Dimension"], 8, TestID -> "Digit-string-101"]

VerificationTest[QuantumState[{1/Sqrt[2], 1/Sqrt[2]}]["Dimension"], 2, TestID -> "Explicit-vector"]

EndTestSection[]


BeginTestSection["QuantumState - properties"]

VerificationTest[QuantumState["Bell"]["Qudits"], 2, TestID -> "Bell-Qudits"]

VerificationTest[QuantumState["Bell"]["Purity"] == 1, True, TestID -> "Bell-PurePurity"]

VerificationTest[Chop[QuantumState["UniformMixture"[2]]["Purity"] - 1/4], 0, TestID -> "UnifMix-Purity"]

VerificationTest[QuantumState["Bell"]["MatrixRepresentation"] // Dimensions, {4, 4}, TestID -> "Bell-MatrixDim"]

EndTestSection[]


BeginTestSection["QuantumState - bipartite operations"]

VerificationTest[
    QuantumPartialTrace[QuantumState["Bell"], {1}]["Dimension"],
    2,
    TestID -> "Bell-PartialTrace-1"
]

VerificationTest[
    Chop[QuantityMagnitude @ UnitConvert[QuantumPartialTrace[QuantumState["Bell"], {1}]["VonNeumannEntropy"], "Bits"] - 1],
    0,
    TestID -> "Bell-EntropyHalf"
]

EndTestSection[]


BeginTestSection["QuantumState - composition fast path"]

(* Each test compares the qs1[qs2] dispatch against the QuantumOperator route
   to verify the fast path matches the reference semantics. *)

slowState[a_, b_] := (QuantumOperator[a] @ QuantumOperator[b])["Sort"]["State"]
densityDiff[a_, b_] := Chop @ Norm @ Flatten[N @ a["DensityMatrix"] - N @ b["DensityMatrix"]]

With[{phiBra = QuantumState[{1/Sqrt[2], 1/Sqrt[2]}]["Dagger"], psi2 = QuantumState[Normalize[{1, 2, 3, 4}], {2, 2}]},
    VerificationTest[densityDiff[phiBra[psi2], slowState[phiBra, psi2]], 0, TestID -> "FastPath-bra-on-ket-prefix"];
    VerificationTest[phiBra[psi2]["VectorQ"], True, TestID -> "FastPath-bra-on-ket-vector"];
    VerificationTest[Dimensions @ phiBra[psi2]["State"], {2}, TestID -> "FastPath-bra-on-ket-dim"];
]

With[{phiBra = QuantumState[{1/Sqrt[2], 1/Sqrt[2]}]["Dagger"], rho2 = QuantumState[Normalize[{1, 2, 3, 4}], {2, 2}]["MatrixState"]},
    VerificationTest[densityDiff[phiBra[rho2], slowState[phiBra, rho2]], 0, TestID -> "FastPath-bra-on-mixed-prefix"];
    VerificationTest[phiBra[rho2]["MatrixQ"], True, TestID -> "FastPath-bra-on-mixed-matrixQ"];
]

With[{
    hOp = QuantumState[QuantumOperator["H"]],
    ket = QuantumState[Normalize[{1, 1}]],
    ket2 = QuantumState[Normalize[{1, 2, 3, 4}], {2, 2}]
},
    VerificationTest[densityDiff[hOp[ket], slowState[hOp, ket]], 0, TestID -> "FastPath-op-on-ket-1q"];
    VerificationTest[densityDiff[hOp[ket["MatrixState"]], slowState[hOp, ket["MatrixState"]]], 0, TestID -> "FastPath-op-on-mixed-1q"];
    VerificationTest[densityDiff[hOp[ket2], slowState[hOp, ket2]], 0, TestID -> "FastPath-op-on-ket-prefix"];
    VerificationTest[densityDiff[hOp[ket2["MatrixState"]], slowState[hOp, ket2["MatrixState"]]], 0, TestID -> "FastPath-op-on-mixed-prefix"];
]

VerificationTest[
    Block[{rp, rm},
        SeedRandom[7];
        rp = QuantumState["RandomPure"];
        rm = QuantumState["RandomMixed"[2]];
        Chop[rp["Dagger"][rm]["Norm"] - slowState[rp["Dagger"], rm]["Norm"]]
    ],
    0,
    TestID -> "FastPath-random-bra-mixed-norm"
]

(* Symmetric direction: qs1 has more input qudits than qs2 has output qudits.
   qs2 fills only a prefix of qs1's input legs; qs1 keeps the residual inputs. *)

With[{cnotOp = QuantumState[QuantumOperator["CNOT"]], ket = QuantumState[Normalize[{1, 1}]]},
    VerificationTest[densityDiff[cnotOp[ket], slowState[cnotOp, ket]], 0, TestID -> "FastPath-symmetric-cnot-on-ket-densitymatrix"];
    VerificationTest[cnotOp[ket]["OutputDimensions"], {2, 2}, TestID -> "FastPath-symmetric-cnot-on-ket-out"];
    VerificationTest[cnotOp[ket]["InputDimensions"], {2}, TestID -> "FastPath-symmetric-cnot-on-ket-in"];
]

With[{cnotOp = QuantumState[QuantumOperator["CNOT"]], hOp = QuantumState[QuantumOperator["H"]]},
    VerificationTest[densityDiff[cnotOp[hOp], slowState[cnotOp, hOp]], 0, TestID -> "FastPath-symmetric-cnot-on-h"];
]

(* PhaseSpace falls through to the QuantumOperator route (fast path is
   Schrodinger-only); verify the result still matches the slow path. *)

With[{phiBra = QuantumState[{1/Sqrt[2], 1/Sqrt[2]}]["Dagger"], psiW = QuantumState[Normalize[{1, 0, 0, 1}], {2, 2}, "Picture" -> "PhaseSpace"]},
    VerificationTest[densityDiff[phiBra[psiW], slowState[phiBra, psiW]], 0, TestID -> "FastPath-phasespace-bra-on-wigner"];
    VerificationTest[phiBra[psiW]["Picture"], "Schrodinger", TestID -> "FastPath-phasespace-picture"];
]

With[{
    psiW1 = QuantumState[Normalize[{1, 0, 0, 1}], {2, 2}, "Picture" -> "PhaseSpace"],
    psiW2 = QuantumState[Normalize[{1, 1}], {2}, "Picture" -> "PhaseSpace"]["Dagger"]
},
    VerificationTest[
        Chop[psiW2[psiW1]["Norm"] - slowState[psiW2, psiW1]["Norm"]],
        0,
        TestID -> "FastPath-phasespace-pair-norm"
    ];
]

EndTestSection[]


BeginTestSection["QuantumState - failure"]

VerificationTest[
    Quiet @ QuantumState["NotAnActualName"["bar"]],
    Failure["InvalidName", _],
    SameTest -> MatchQ,
    TestID -> "InvalidName-call-form"
]

EndTestSection[]
