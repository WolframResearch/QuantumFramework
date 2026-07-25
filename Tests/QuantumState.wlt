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


BeginTestSection["QuantumState - named states with a basis or option tail"]

(* The fixed-vector named states take the same tail as every other constructor
   here: a basis specification followed by options, applied to the state the
   name builds. The amplitudes are the ones the bare name produces, tagged with
   the requested basis, so a dimension-preserving tail leaves the state alone. *)

VerificationTest[Normal @ QuantumState["0", QuantumBasis[2]]["StateVector"], {1, 0}, {}, TestID -> "NamedTail-0"]

VerificationTest[Normal @ QuantumState["1", QuantumBasis[2]]["StateVector"], {0, 1}, {}, TestID -> "NamedTail-1"]

VerificationTest[Normal @ QuantumState["Plus", QuantumBasis[2]]["StateVector"], {1/Sqrt[2], 1/Sqrt[2]}, {}, TestID -> "NamedTail-Plus"]

VerificationTest[Normal @ QuantumState["Minus", QuantumBasis[2]]["StateVector"], {1/Sqrt[2], -1/Sqrt[2]}, {}, TestID -> "NamedTail-Minus"]

VerificationTest[Normal @ QuantumState["Left", QuantumBasis[2]]["StateVector"], {1/Sqrt[2], -I/Sqrt[2]}, {}, TestID -> "NamedTail-Left"]

VerificationTest[Normal @ QuantumState["Right", QuantumBasis[2]]["StateVector"], {1/Sqrt[2], I/Sqrt[2]}, {}, TestID -> "NamedTail-Right"]

VerificationTest[Normal @ QuantumState["PhiPlus", QuantumBasis[2]]["StateVector"], {1/Sqrt[2], 0, 0, 1/Sqrt[2]}, {}, TestID -> "NamedTail-PhiPlus"]

VerificationTest[Normal @ QuantumState["PhiMinus", QuantumBasis[2]]["StateVector"], {1/Sqrt[2], 0, 0, -1/Sqrt[2]}, {}, TestID -> "NamedTail-PhiMinus"]

VerificationTest[Normal @ QuantumState["PsiPlus", QuantumBasis[2]]["StateVector"], {0, 1/Sqrt[2], 1/Sqrt[2], 0}, {}, TestID -> "NamedTail-PsiPlus"]

VerificationTest[Normal @ QuantumState["PsiMinus", QuantumBasis[2]]["StateVector"], {0, 1/Sqrt[2], -1/Sqrt[2], 0}, {}, TestID -> "NamedTail-PsiMinus"]

(* A basis name is as good a tail as a QuantumBasis, and the empty call form is
   the literal shape the name dispatcher forwards to. *)
VerificationTest[
    Normal @ QuantumState["Plus", "PauliBasis"]["StateVector"],
    {1/Sqrt[2], 1/Sqrt[2], 0, 0},
    {},
    TestID -> "NamedTail-Plus-PauliBasis"
]

VerificationTest[
    Normal @ QuantumState["Plus"[], QuantumBasis[2]]["StateVector"],
    {1/Sqrt[2], 1/Sqrt[2]},
    {},
    TestID -> "NamedTail-Plus-empty-call-form"
]

(* Whatever the tail, the state stays a unit vector. *)
VerificationTest[
    Union @ Flatten @ Outer[
        QuantumState[#1, #2]["Norm"] &,
        {"0", "1", "Plus", "Minus", "Left", "Right", "PhiPlus", "PhiMinus", "PsiPlus", "PsiMinus"},
        {QuantumBasis[2], "PauliBasis", QuantumBasis[2, 2]}
    ],
    {1},
    {},
    TestID -> "NamedTail-normalization"
]

(* The name's own label is a default the tail can override, so the state stays
   recognizable through a change of basis but an explicit option still wins. *)
VerificationTest[
    Union @ Map[
        QuantumState[#, QuantumBasis[2]]["Label"] === QuantumState[#]["Label"] &,
        {"0", "1", "Plus", "Minus", "Left", "Right", "PhiPlus", "PhiMinus", "PsiPlus", "PsiMinus"}
    ],
    {True},
    {},
    TestID -> "NamedTail-label-survives"
]

VerificationTest[QuantumState["Plus", QuantumBasis[2], "Label" -> "q"]["Label"], "q", {}, TestID -> "NamedTail-label-overridable"]

(* "0" and "1" reach the digit string rule, their aliases reach the name
   dispatcher. The two routes have to land on the same object. *)
VerificationTest[QuantumState["Zero", QuantumBasis[2]] === QuantumState["0", QuantumBasis[2]], True, {}, TestID -> "NamedTail-Zero-matches-digit-route"]

VerificationTest[QuantumState["Zero", "PauliBasis"] === QuantumState["0", "PauliBasis"], True, {}, TestID -> "NamedTail-Zero-matches-digit-route-Pauli"]

VerificationTest[QuantumState["Up", QuantumBasis[2]] === QuantumState["0", QuantumBasis[2]], True, {}, TestID -> "NamedTail-Up-matches-digit-route"]

VerificationTest[QuantumState["One", QuantumBasis[2]] === QuantumState["1", QuantumBasis[2]], True, {}, TestID -> "NamedTail-One-matches-digit-route"]

(* The bare forms are what they were before the tail was accepted. *)
VerificationTest[
    Normal[QuantumState[#]["StateVector"]] & /@
        {"0", "1", "Plus", "Minus", "Left", "Right", "PhiPlus", "PhiMinus", "PsiPlus", "PsiMinus"},
    {{1, 0}, {0, 1}, {1/Sqrt[2], 1/Sqrt[2]}, {1/Sqrt[2], -1/Sqrt[2]}, {1/Sqrt[2], -I/Sqrt[2]}, {1/Sqrt[2], I/Sqrt[2]},
     {1/Sqrt[2], 0, 0, 1/Sqrt[2]}, {1/Sqrt[2], 0, 0, -1/Sqrt[2]}, {0, 1/Sqrt[2], 1/Sqrt[2], 0}, {0, 1/Sqrt[2], -1/Sqrt[2], 0}},
    {},
    TestID -> "NamedTail-bare-forms-unchanged"
]

(* A trailing integer is a qudit dimension, like everywhere else in this file.
   The qudit count lives inside the name instead, as "Plus"[n], and that form
   is unaffected by the tail the base rules now accept. *)
VerificationTest[Normal @ QuantumState["Plus"[2]]["StateVector"], {1/2, 1/2, 1/2, 1/2}, {}, TestID -> "NamedTail-Plus-2-unchanged"]

VerificationTest[QuantumState["Plus"[3]]["Dimension"], 8, {}, TestID -> "NamedTail-Plus-3-unchanged"]

VerificationTest[Normal @ QuantumState["Left"[2]]["StateVector"], {1/2, -I/2, -I/2, -1/2}, {}, TestID -> "NamedTail-Left-2-unchanged"]

VerificationTest[QuantumState["PhiPlus"[2]]["Dimension"], 16, {}, TestID -> "NamedTail-PhiPlus-2-unchanged"]

VerificationTest[QuantumState["Plus", 2]["Dimension"], 2, {}, TestID -> "NamedTail-trailing-integer-is-a-dimension"]

VerificationTest[QuantumState["Plus"[2], QuantumBasis[2]] === QuantumState["Plus"[2]], True, {}, TestID -> "NamedTail-Plus-2-accepts-tail"]

(* Names that carry their own arguments already threaded the tail; they still do. *)
VerificationTest[QuantumState["GHZ", QuantumBasis[2]]["Dimension"], 8, {}, TestID -> "NamedTail-GHZ-still-threads"]

(* Accepting a tail must not turn an unmatched call shape into a silent success. *)
VerificationTest[
    QuantumState["Plus"["bad-arg"]],
    Failure["InvalidArguments", _],
    {QuantumState::invalidArgs},
    SameTest -> MatchQ,
    TestID -> "NamedTail-bad-call-shape-still-fails"
]

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


BeginTestSection["QuantumState - padding"]

(* A vector whose length is not a power of the qudit dimension is zero-padded
   into the next larger qudit space; that changes the qudit count, so it warns. *)
VerificationTest[
    QuantumState[Table[0, {12}], 2]["Dimension"],
    16,
    {QuantumState::padded},
    TestID -> "Pad-nonpower-vector-warns"
]

VerificationTest[
    QuantumState[Table[0, {12}], 2]["Norm"],
    0,
    {QuantumState::padded},
    TestID -> "Pad-zero-vector-norm0"
]

VerificationTest[
    QuantumState[ConstantArray[0, {12, 12}], 2]["Dimension"],
    16,
    {QuantumState::padded},
    TestID -> "Pad-nonpower-matrix-warns"
]

(* Exact power of the qudit dimension: no padding, no message. *)
VerificationTest[
    QuantumState[Normalize @ Table[1, {16}], 2]["Dimension"],
    16,
    TestID -> "Pad-exact-power-silent"
]

(* Single-qudit amplitude-prefix pad (multiplicity 1) stays silent; the
   tensor-network initial-state prepend relies on QuantumState[{1}, dim]. *)
VerificationTest[
    Normal @ QuantumState[{1}, 2]["StateVector"],
    {1, 0},
    TestID -> "Pad-prefix-single-qudit-silent"
]

EndTestSection[]


BeginTestSection["QuantumState - failure"]

VerificationTest[
    QuantumState["NotAnActualName"["bar"]],
    Failure["InvalidName", _],
    {QuantumState::invalidName},
    SameTest -> MatchQ,
    TestID -> "InvalidName-call-form"
]

VerificationTest[
    QuantumState["GHZ"["bad-arg"]],
    Failure["InvalidArguments", _],
    {QuantumState::invalidArgs},
    SameTest -> MatchQ,
    TestID -> "InvalidArgs-GHZ"
]

(* A bare unrecognized name used to fall through to the scalar-amplitude rule
   and build the 1-dimensional state {"NotAnActualName"}, so that ["Norm"]
   returned Abs["NotAnActualName"] with no message at all. It has to fail the
   same way the call form does. *)
VerificationTest[
    QuantumState["NotAnActualName"],
    Failure["InvalidName", _],
    {QuantumState::invalidName},
    SameTest -> MatchQ,
    TestID -> "InvalidName-bare-form"
]

(* The name is unrecognized whatever the argument tail says. *)
VerificationTest[
    QuantumState["NotAnActualName", QuantumBasis[2]],
    Failure["InvalidName", _],
    {QuantumState::invalidName},
    SameTest -> MatchQ,
    TestID -> "InvalidName-bare-form-with-basis"
]

(* Near misses are what a user actually types: a truncation, a wrong case, and
   a name belonging to a sibling constructor. One test each, because three
   copies of the same message in a single test trips General::stop. *)
VerificationTest[
    FailureQ @ QuantumState["Bel"],
    True,
    {QuantumState::invalidName},
    TestID -> "InvalidName-truncated"
]

VerificationTest[
    FailureQ @ QuantumState["ghz"],
    True,
    {QuantumState::invalidName},
    TestID -> "InvalidName-wrong-case"
]

VerificationTest[
    FailureQ @ QuantumState["Hadamard"],
    True,
    {QuantumState::invalidName},
    TestID -> "InvalidName-belongs-to-sibling"
]

(* Nothing symbolic may reach a numeric property. Asserting NumericQ is False
   would not discriminate, since Abs["NotAnActualName"] is not numeric either;
   the property has to come back Missing, which only a Failure produces. *)
VerificationTest[
    MatchQ[QuantumState["NotAnActualName"]["Norm"], _Missing],
    True,
    {QuantumState::invalidName},
    TestID -> "InvalidName-no-symbolic-norm"
]

(* Registered names stay untouched, and stay silent. *)
VerificationTest[
    QuantumState["Plus"]["Norm"],
    1,
    {},
    TestID -> "ValidName-unaffected-Plus"
]

VerificationTest[
    {QuantumState["Bell"]["Dimension"], QuantumState["Register"[3]]["Qudits"]},
    {4, 3},
    {},
    TestID -> "ValidName-unaffected-Bell-Register"
]

(* The string forms that bypass $QuantumStateNames must survive the guard: the
   zero-qudit empty state, basis-element digit strings (in any qudit dimension),
   computational and Pauli-eigenstate letter sequences, and the "+i"/"-i"
   aliases, none of which appear in the registry. Only "" is a bare literal and
   so ordered ahead of the guard on its own; the rest are conditioned or
   patterned rules that stay reachable only by preceding it in NamedStates.m. *)
VerificationTest[
    {
        QuantumState[""]["Dimension"],
        QuantumState["0101"]["Qudits"],
        QuantumState["2", 3]["AmplitudesList"],
        QuantumState["+0L"]["Qudits"],
        QuantumState["-i"]["AmplitudesList"],
        QuantumState["+i"]["AmplitudesList"]
    },
    {1, 4, {0, 0, 1}, 3, {1/Sqrt[2], -I/Sqrt[2]}, {1/Sqrt[2], I/Sqrt[2]}},
    {},
    TestID -> "Shorthand-strings-not-swallowed"
]

(* The digit and letter-sequence rules take an argument tail, so they compete
   with the invalidName arm at the same arity, and neither pattern outranks the
   other: they win only by preceding it in NamedStates.m. These pin that file
   order; the letter-sequence tail is exercised nowhere else in the suite. *)
VerificationTest[
    {
        QuantumState["0101", QuantumBasis[16]]["Qudits"],
        QuantumState["2", 3]["AmplitudesList"],
        QuantumState["12", 3]["Qudits"],
        QuantumState["+", QuantumBasis[2]]["Norm"]
    },
    {4, {0, 0, 1}, 2, 1},
    {},
    TestID -> "Shorthand-strings-with-a-tail"
]

(* The "+i"/"-i" aliases forward their tail to the named state, so the tail
   reaches the basis instead of slipping past into the scalar-amplitude rule,
   where the string itself would become the amplitude and the norm would come
   back as Abs["+i"] with no message. The norm pins that: a leak cannot be 1. *)
VerificationTest[
    {
        QuantumState["+i", QuantumBasis[2]] === QuantumState["Right", QuantumBasis[2]],
        QuantumState["+i", QuantumBasis[2]]["Norm"]
    },
    {True, 1},
    {},
    TestID -> "Alias-with-tail-forwards-not-leaks"
]

(* A two-level named state given a qutrit basis embeds on the first two levels,
   silently, the way any short amplitude prefix does. *)
VerificationTest[
    {
        QuantumState["-i", 3] === QuantumState["Left", 3],
        Normal @ QuantumState["-i", 3]["StateVector"],
        QuantumState["-i", 3]["Norm"]
    },
    {True, {1/Sqrt[2], -I/Sqrt[2], 0}, 1},
    {},
    TestID -> "Alias-with-tail-forwards-not-leaks-qutrit"
]

(* The empty state is a zero-qudit form and takes no basis, so with a tail it is
   as unrecognized as any other string rather than an Abs[""] amplitude. *)
VerificationTest[
    QuantumState["", QuantumBasis[2]],
    Failure["InvalidName", _],
    {QuantumState::invalidName},
    SameTest -> MatchQ,
    TestID -> "Empty-state-with-tail-is-unrecognized"
]

(* "Properties" is a query on the symbol itself, not a state name. *)
VerificationTest[
    MemberQ[QuantumState["Properties"], "VonNeumannEntropy"],
    True,
    {},
    TestID -> "Properties-query-not-a-name"
]

EndTestSection[]


BeginTestSection["QuantumState - cross-basis composition"]

(* Shared fixtures: diag(1,-1) tagged in the PauliX eigenbasis is sigma_x; the
   mixed state is a generic qubit density matrix; the qutrit diagonal carries a
   symbolic phase in the Fourier[3] frame. *)
xFrameDiag = QuantumOperator[DiagonalMatrix[{1, -1}], QuantumBasis["PauliX"]]["State"];
mixedRho = QuantumState[{{3/4, I/4}, {-I/4, 1/4}}];
fourierPhaseDiag = QuantumOperator[DiagonalMatrix[{1, Exp[I \[FormalTheta]], Exp[2 I \[FormalTheta]]}], QuantumBasis["Fourier"[3]]]["State"];

(* The qs1[qs2] fast path multiplies computational matrices; the result's data
   and tagged basis must agree. diag(1,-1) tagged in the PauliX eigenbasis is
   sigma_x, so acting on |0> must give |1> in every basis-aware read. *)

VerificationTest[
    Normal @ xFrameDiag[QuantumState[{1, 0}]]["Computational"]["StateVector"],
    {0, 1},
    TestID -> "CrossBasis-XFrameDiag-OnKet0"
]

(* The composed state stays tagged in the non-computational frame; only the
   data is rebased into it. *)
VerificationTest[
    xFrameDiag[QuantumState[{1, 0}]]["Output"]["ComputationalQ"],
    False,
    TestID -> "CrossBasis-FramePreserved"
]

(* Symbolic phase: diag(1, E^(I t)) in the X frame acts on |0> as
   (1/2) {1 + E^(I t), 1 - E^(I t)} computationally. *)
VerificationTest[
    FullSimplify[
        Normal @ QuantumOperator[DiagonalMatrix[{1, Exp[I \[FormalTheta]]}], QuantumBasis["PauliX"]]["State"][QuantumState[{1, 0}]]["Computational"]["StateVector"] -
            {1 + Exp[I \[FormalTheta]], 1 - Exp[I \[FormalTheta]]}/2
    ],
    {0, 0},
    TestID -> "CrossBasis-SymbolicPhase"
]

(* Zero-phase limit: the identity tagged in the X frame is a pure relabel. *)
VerificationTest[
    Normal @ QuantumOperator[IdentityMatrix[2], QuantumBasis["PauliX"]]["State"][QuantumState[{1, 0}]]["Computational"]["StateVector"],
    {1, 0},
    TestID -> "CrossBasis-IdentityRelabel"
]

(* Both operands tagged in the same non-computational frame: diag(1, I) in the
   X frame fixes |x0> = |+>, whose computational vector is {1,1}/Sqrt[2]. *)
VerificationTest[
    Normal @ QuantumOperator[DiagonalMatrix[{1, I}], QuantumBasis["PauliX"]]["State"][QuantumState[{1, 0}, QuantumBasis["PauliX"]]]["Computational"]["StateVector"],
    {1, 1}/Sqrt[2],
    TestID -> "CrossBasis-SameXFrame"
]

VerificationTest[
    QuantumOperator[DiagonalMatrix[{1, I}], QuantumBasis["PauliX"]]["State"][QuantumState[{1, 0}, QuantumBasis["PauliX"]]]["Norm"],
    1,
    TestID -> "CrossBasis-NormPreserved"
]

(* A non-computational qudit beyond the contracted prefix: X on qudit 1 of
   |0> (x) |x0> leaves the X-tagged rest qudit's data consistent. *)
VerificationTest[
    Normal @ QuantumOperator["X"]["State"][
        QuantumState[{1, 0, 0, 0}, QuantumBasis["Output" -> QuantumTensorProduct[QuditBasis[2], QuditBasis["PauliX"]]]]
    ]["Computational"]["StateVector"],
    {0, 0, 1, 1}/Sqrt[2],
    TestID -> "CrossBasis-NonComputationalRest"
]

(* Mixed sub-case: sigma_x . rho . sigma_x with an X-tagged qs1. *)
VerificationTest[
    Normal @ xFrameDiag[mixedRho]["Computational"]["DensityMatrix"],
    {{1/4, -I/4}, {I/4, 3/4}},
    TestID -> "CrossBasis-MixedConjugation"
]

(* Complex frame: diag(1,-1) tagged in the PauliY eigenbasis is sigma_y, whose
   action on |0> is {0, I}; a dropped conjugate in the rebase would flip it. *)
VerificationTest[
    Normal @ QuantumOperator[DiagonalMatrix[{1, -1}], QuantumBasis["PauliY"]]["State"][QuantumState[{1, 0}]]["Computational"]["StateVector"],
    {0, I},
    TestID -> "CrossBasis-YFrameComplexConjugation"
]

(* Two different non-computational frames: sigma_x (X-tagged) acting on the
   Y-frame ket |y0> = {1, I}/Sqrt[2]. *)
VerificationTest[
    Normal @ xFrameDiag[QuantumState[{1, 0}, QuantumBasis["PauliY"]]]["Computational"]["StateVector"],
    {I, 1}/Sqrt[2],
    TestID -> "CrossBasis-DifferentFrames"
]

(* Qudit symbolic case: a symbolic-phase diagonal in the Fourier[3] frame
   applied to |0>, checked against the closed form sum_k d_k <f_k|0> |f_k>
   built from the frame's own elements. *)
VerificationTest[
    With[{
        diag = {1, Exp[I \[FormalTheta]], Exp[2 I \[FormalTheta]]},
        els = Normal /@ QuantumBasis["Fourier"[3]]["Output"]["Elements"],
        res = Normal @ fourierPhaseDiag[QuantumState[{1, 0, 0}, 3]]["Computational"]["StateVector"]
    },
        FullSimplify @ ComplexExpand[res - Total[MapThread[#1 Conjugate[#2[[1]]] #2 &, {diag, els}]]]
    ],
    {0, 0, 0},
    TestID -> "CrossBasis-QutritFourierSymbolic"
]

(* The rebased Fourier-frame data of that state stays normalized for symbolic
   real phase. *)
VerificationTest[
    FullSimplify @ ComplexExpand[# . Conjugate[#]] & @
        Normal @ fourierPhaseDiag[QuantumState[{1, 0, 0}, 3]]["StateVector"],
    1,
    TestID -> "CrossBasis-QutritSymbolicNorm"
]

(* Operator on operator, both X-tagged: the result carries a non-computational
   INPUT frame, so the input-side rebase is exercised; computationally this is
   sigma_x . (H diag(1,I) H). *)
VerificationTest[
    With[{res = xFrameDiag[QuantumOperator[DiagonalMatrix[{1, I}], QuantumBasis["PauliX"]]["State"]]},
        {Normal @ res["Computational"]["Matrix"], res["Input"]["ComputationalQ"]}
    ],
    {{{(1 - I)/2, (1 + I)/2}, {(1 + I)/2, (1 - I)/2}}, False},
    TestID -> "CrossBasis-OperatorOnOperator-InputFrame"
]

(* Invariants of the mixed-branch result: unit trace, Hermiticity, positivity. *)
VerificationTest[
    With[{m = Normal @ xFrameDiag[mixedRho]["Computational"]["DensityMatrix"]},
        {Tr[m], HermitianMatrixQ[m], PositiveSemidefiniteMatrixQ[m]}
    ],
    {1, True, True},
    TestID -> "CrossBasis-MixedInvariants"
]

(* Mixed state with a non-computational rest qudit: X on qudit 1 of
   rho = (|0,x0><0,x0| + |1,x1><1,x1|)/2 tagged comp (x) PauliX gives
   (|1,x0><1,x0| + |0,x1><0,x1|)/2, with all density-matrix invariants intact. *)
VerificationTest[
    With[{m = Normal @ QuantumOperator["X"]["State"][
            QuantumState[DiagonalMatrix[{1/2, 0, 0, 1/2}], QuantumBasis["Output" -> QuantumTensorProduct[QuditBasis[2], QuditBasis["PauliX"]]]]
        ]["Computational"]["DensityMatrix"]},
        {m, Tr[m], HermitianMatrixQ[m], PositiveSemidefiniteMatrixQ[m]}
    ],
    {{{1/4, -1/4, 0, 0}, {-1/4, 1/4, 0, 0}, {0, 0, 1/4, 1/4}, {0, 0, 1/4, 1/4}}, 1, True, True},
    TestID -> "CrossBasis-MixedNonComputationalRest"
]

(* An unnormalized input passes through faithfully: the rebase must not
   silently renormalize. *)
VerificationTest[
    xFrameDiag[QuantumState[{2, 0}]]["Norm"],
    2,
    TestID -> "CrossBasis-UnnormalizedFaithful"
]

(* The fast path agrees with the general operator route on stored data, frame
   tag, and computational read (concrete and symbolic phase). *)
VerificationTest[
    With[{
        fastC = xFrameDiag[QuantumState[{1, 0}]],
        genC = (QuantumOperator[xFrameDiag] @ QuantumOperator[QuantumState[{1, 0}]])["Sort"]["State"],
        fastS = QuantumOperator[DiagonalMatrix[{1, Exp[I \[FormalTheta]]}], QuantumBasis["PauliX"]]["State"][QuantumState[{1, 0}]],
        genS = (QuantumOperator[QuantumOperator[DiagonalMatrix[{1, Exp[I \[FormalTheta]]}], QuantumBasis["PauliX"]]["State"]] @ QuantumOperator[QuantumState[{1, 0}]])["Sort"]["State"]
    },
        {
            Normal @ fastC["State"] === Normal @ genC["State"],
            fastC["Output"]["ComputationalQ"] === genC["Output"]["ComputationalQ"] === False,
            FullSimplify[Normal @ fastS["Computational"]["StateVector"] - Normal @ genS["Computational"]["StateVector"]]
        }
    ],
    {True, True, {0, 0}},
    TestID -> "CrossBasis-FastMatchesGeneralRoute"
]

(* Independent numeric reference for the composition law. *)
VerificationTest[
    N @ Normal @ QuantumOperator[DiagonalMatrix[{1, Exp[I 0.7]}], QuantumBasis["PauliX"]]["State"][QuantumState[{1, 0}]]["Computational"]["StateVector"],
    ({{1., 1.}, {1., -1.}}/Sqrt[2.]) . DiagonalMatrix[{1., Exp[I 0.7]}] . ({{1., 1.}, {1., -1.}}/Sqrt[2.]) . {1., 0.},
    SameTest -> (Norm[#1 - #2] < 1*^-12 &),
    TestID -> "CrossBasis-NumericReference"
]

(* Computational compositions are untouched by the frame handling. *)
VerificationTest[
    Normal @ QuantumOperator["X"]["State"][QuantumState[{1, 0}]]["StateVector"],
    {0, 1},
    TestID -> "CrossBasis-ComputationalUnchanged"
]

VerificationTest[
    QuantumState[{1, 0}, QuantumBasis["PauliX"]]["Dagger"][QuantumState[{1, 0}]]["Scalar"],
    1/Sqrt[2],
    TestID -> "CrossBasis-BraKetScalar"
]

EndTestSection[]


BeginTestSection["QuantumState - matrix-route change of basis"]

(* Shared fixtures: CNOT dressed in the PauliY (x) PauliY frame is the canonical
   operator-shaped state whose INPUT basis has complex elements; the PauliX twin
   is the real-frame control; the qutrit clock in the Fourier[3] frame covers
   d > 2. For a pure operator-shaped state the vector route (StateVector) and the
   matrix route (MatrixState -> DensityMatrix) must land on the same computational
   object: rho_comp === w . w^dagger exactly. *)
yyCnot = QuantumOperator[QuantumOperator["CNOT"], QuantumBasis[{"PauliY", "PauliY"}]]["State"];
xxCnot = QuantumOperator[QuantumOperator["CNOT"], QuantumBasis[{"PauliX", "PauliX"}]]["State"];
fourierClock = QuantumOperator[QuantumOperator["Z"[3]], QuantumBasis["Fourier"[3]]]["State"];

(* The vector route is the ground truth: rebasing must not change the abstract
   operator, so the computational statevector is vec(CNOT) exactly. *)
VerificationTest[
    Simplify[Normal @ yyCnot["Computational"]["StateVector"] - Flatten[Normal @ QuantumOperator["CNOT"]["Matrix"]]],
    ConstantArray[0, 16],
    TestID -> "MatrixRoute-VectorRouteGroundTruth-YY"
]

(* Cross-route agreement on a complex input frame: the dual input leg of the
   doubled object must enter conjugated, matching vec(Out . S . In^-1). *)
VerificationTest[
    With[{w = Normal @ yyCnot["Computational"]["StateVector"]},
        Simplify[Normal @ yyCnot["MatrixState"]["Computational"]["DensityMatrix"] - KroneckerProduct[w, Conjugate[w]]]
    ],
    ConstantArray[0, {16, 16}],
    TestID -> "MatrixRoute-CrossRouteAgreement-YY"
]

(* Same law on a qutrit with a complex Fourier input frame. *)
VerificationTest[
    With[{w = Normal @ fourierClock["Computational"]["StateVector"]},
        Simplify[Normal @ fourierClock["MatrixState"]["Computational"]["DensityMatrix"] - KroneckerProduct[w, Conjugate[w]]]
    ],
    ConstantArray[0, {9, 9}],
    TestID -> "MatrixRoute-CrossRouteAgreement-Fourier3"
]

(* Real input frames are blind to the dual conjugation; the PauliX twin pins the
   already-correct behavior. *)
VerificationTest[
    With[{w = Normal @ xxCnot["Computational"]["StateVector"]},
        Simplify[Normal @ xxCnot["MatrixState"]["Computational"]["DensityMatrix"] - KroneckerProduct[w, Conjugate[w]]]
    ],
    ConstantArray[0, {16, 16}],
    TestID -> "MatrixRoute-RealFrameRegression-XX"
]

(* Round trip: computational and back into the tagged frame must reproduce the
   stored density matrix, so both directions of the transform stay inverse to
   each other. *)
VerificationTest[
    With[{ms = yyCnot["MatrixState"]},
        Simplify[Normal @ QuantumState[ms["Computational"], ms["Basis"]]["DensityMatrix"] - Normal @ ms["DensityMatrix"]]
    ],
    ConstantArray[0, {16, 16}],
    TestID -> "MatrixRoute-RoundTrip-YY"
]

(* Input-rank-0 mixed states have no dual leg: a frame-diagonal mixture reads
   computationally as the same mixture of the frame's own projectors. *)
VerificationTest[
    With[{
        m = Normal @ QuantumState[DiagonalMatrix[{1/3, 2/3}], QuantumBasis["PauliY"]]["Computational"]["DensityMatrix"],
        els = Normal /@ QuantumBasis["PauliY"]["Output"]["Elements"]
    },
        Simplify[m - (KroneckerProduct[els[[1]], Conjugate[els[[1]]]]/3 + 2 KroneckerProduct[els[[2]], Conjugate[els[[2]]]]/3)]
    ],
    ConstantArray[0, {2, 2}],
    TestID -> "MatrixRoute-MixedInputRankZero-YFrame"
]

EndTestSection[]
