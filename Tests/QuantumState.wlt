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

(* A basis that can hold the amplitudes leaves the state a unit vector, whether
   it matches the dimension, raises it, or is a frame other than the
   computational one. *)
VerificationTest[
    Union @ Flatten @ Outer[
        QuantumState[#1, #2]["Norm"] &,
        {"0", "1", "Plus", "Minus", "Left", "Right", "PhiPlus", "PhiMinus", "PsiPlus", "PsiMinus"},
        {QuantumBasis[2], "PauliBasis", QuantumBasis[2, 2], QuantumBasis["PauliX"]}
    ],
    {1},
    {},
    TestID -> "NamedTail-normalization"
]

(* The tail says what the amplitudes are coefficients of; it does not rotate
   them into the frame. In a dimension-matched frame that distinction is
   visible, where a larger basis would hide it behind zero-padding: the name
   keeps its coefficient pattern, so QuantumState["Plus", PauliX] carries
   {1,1}/Sqrt[2] against the PauliX elements rather than the {1,0} a change of
   basis would give. The digit string rule reads its tail the same way, which
   is what makes the two routes into "0" agree on a well-formed tail; they part
   company on a malformed one, where each rejects in its own way. *)
VerificationTest[
    {
        Normal @ QuantumState["Plus", QuantumBasis["PauliX"]]["StateVector"],
        Normal @ QuantumState["0", QuantumBasis["PauliX"]]["StateVector"],
        QuantumState["Zero", QuantumBasis["PauliX"]] === QuantumState["0", QuantumBasis["PauliX"]]
    },
    {{1/Sqrt[2], 1/Sqrt[2]}, {1, 0}, True},
    {},
    TestID -> "NamedTail-tags-coefficients-not-a-change-of-basis"
]

(* These names denote kets of at least two levels. A basis with fewer levels
   would be filled by dropping amplitudes, returning an unnormalized state, and
   a basis carrying an input leg would return an operator shaped object; both
   are rejected. One test per rejected tail, because several copies of one
   message in a single test trip General::stop. *)
VerificationTest[
    QuantumState["Plus", 1],
    Failure["InvalidArguments", _],
    {QuantumState::invalidArgs},
    SameTest -> MatchQ,
    TestID -> "NamedTail-one-level-basis-rejected"
]

VerificationTest[
    QuantumState["Plus", 0],
    Failure["InvalidArguments", _],
    {QuantumState::invalidArgs},
    SameTest -> MatchQ,
    TestID -> "NamedTail-zero-level-basis-rejected"
]

VerificationTest[
    QuantumState["PsiMinus", QuantumBasis[1]],
    Failure["InvalidArguments", _],
    {QuantumState::invalidArgs},
    SameTest -> MatchQ,
    TestID -> "NamedTail-one-level-basis-rejected-two-qubit-name"
]

(* Without this the state comes back with a dual leg and dimension 6. *)
VerificationTest[
    QuantumState["Plus", QuantumBasis[QuditBasis[2], QuditBasis[3]]],
    Failure["InvalidArguments", _],
    {QuantumState::invalidArgs},
    SameTest -> MatchQ,
    TestID -> "NamedTail-input-leg-rejected"
]

(* A tail that is no basis at all reports the same way, rather than surfacing
   the internal predicate that turned it down. *)
VerificationTest[
    QuantumState["Plus", QuantumState["Bell"]],
    Failure["InvalidArguments", _],
    {QuantumState::invalidArgs},
    SameTest -> MatchQ,
    TestID -> "NamedTail-non-basis-tail-rejected"
]

VerificationTest[
    QuantumState["Minus", 3/2],
    Failure["InvalidArguments", _],
    {QuantumState::invalidArgs},
    SameTest -> MatchQ,
    TestID -> "NamedTail-non-integer-tail-rejected"
]

(* A basis with room to spare embeds rather than truncating, so the norm holds
   even when the dimension is not a power of the qudit dimension. *)
VerificationTest[
    {QuantumState["Plus", 3]["Norm"], Normal @ QuantumState["Plus", 3]["StateVector"]},
    {1, {1/Sqrt[2], 1/Sqrt[2], 0}},
    {},
    TestID -> "NamedTail-larger-basis-embeds"
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

(* The coefficients above are the ones handed in, so they read the same under
   either convention. What discriminates them is the state those coefficients
   describe: read computationally, (|+> + |->)/Sqrt[2] is |0>, which points
   along +z where the bare name points along +x. So the name survives the tail
   as a coefficient pattern, not as a physical state, and a caller who wants
   |+> written in the X frame wants QuantumState[QuantumState["Plus"], basis]
   instead, which is the other branch of the ambiguity and gives {1, 0}. *)
VerificationTest[
    {
        Normal @ QuantumState["Plus", QuantumBasis["PauliX"]]["Computational"]["StateVector"],
        QuantumState["Plus", QuantumBasis["PauliX"]]["BlochCartesianCoordinates"],
        QuantumState["Plus"]["BlochCartesianCoordinates"],
        Normal @ QuantumState[QuantumState["Plus"], QuantumBasis["PauliX"]]["StateVector"]
    },
    {{1, 0}, {0, 0, 1}, {1, 0, 0}, {1, 0}},
    {},
    TestID -> "NamedTail-tagged-state-is-not-the-named-eigenvector"
]

(* A tail that keeps the two-qubit structure keeps the entanglement the Bell
   names are named for: one bit across the cut. A tail that does not is where
   the name stops describing the state, so the same read on a mixed-dimension
   basis returns a product state. *)
VerificationTest[
    QuantumPartialTrace[QuantumState["PhiPlus", QuantumBasis[2]], {1}]["VonNeumannEntropy"],
    Quantity[1, "Bits"],
    {},
    TestID -> "NamedTail-Bell-name-keeps-one-bit-across-the-cut"
]

VerificationTest[
    {
        QuantumState["PhiPlus", QuantumBasis[{2, 3}]]["Qudits"],
        QuantumPartialTrace[QuantumState["PhiPlus", QuantumBasis[{2, 3}]], {1}]["VonNeumannEntropy"],
        QuantumState["PhiPlus", QuantumBasis[4]]["Qudits"]
    },
    {2, Quantity[0, "Bits"], 1},
    {},
    TestID -> "NamedTail-Bell-name-loses-entanglement-off-two-qubits"
]

(* The guard reads a dimension rather than assuming a qubit, so the one-qudit
   names embed into any d >= 2 with the norm intact. Checked away from d = 2,
   and against a non-computational qutrit frame, because the qubit tails are the
   ones that stay quiet. "Fourier"[3] is the qutrit frame; QuantumBasis levels a
   trailing integer at three Fourier qubits instead. *)
VerificationTest[
    Union @ Flatten @ Outer[
        QuantumState[#1, #2]["Norm"] &,
        {"0", "1", "Plus", "Minus", "Left", "Right"},
        {3, 4, 5, QuantumBasis["Fourier"[3]], QuantumBasis["Fourier"[5]]}
    ],
    {1},
    {},
    TestID -> "NamedTail-qudit-dimensions-normalized"
]

(* Four amplitudes against d levels is a trichotomy in how d divides 4, and it
   is the structure a caller has to know before reading a Bell name's tail. For
   d = 2 the four amplitudes fill two qubits, which is the only branch where the
   name still describes the state. For d dividing 4 above that, and for every
   d > 4, they fill a single qudit, so there is no cut to be entangled across.
   d = 3 is the odd one, since 3 divides neither 4 nor any power of 2 below it:
   the four amplitudes grow to d^2 = 9 across two qutrits and the rest is padded,
   which warns. It is pinned separately below so this table stays silent. *)
VerificationTest[
    Table[{QuantumState["PhiPlus", d]["Qudits"], QuantumState["PhiPlus", d]["Dimension"]}, {d, {2, 4, 5, 6, 7}}],
    {{2, 4}, {1, 4}, {1, 5}, {1, 6}, {1, 7}},
    {},
    TestID -> "NamedTail-four-amplitudes-against-d-levels"
]

VerificationTest[
    Union @ Flatten @ Table[QuantumState[n, d]["Norm"], {n, {"PhiPlus", "PhiMinus", "PsiPlus", "PsiMinus"}}, {d, {2, 4, 5, 6, 7}}],
    {1},
    {},
    TestID -> "NamedTail-Bell-names-normalized-across-d"
]

(* The third branch: three levels hold four amplitudes only by taking two
   qutrits and padding the other five slots, which is the one accepted tail that
   warns. The amplitudes land in the two-qutrit corner, so the state is
   (|00> + |11>)/Sqrt[2] read in d = 3, not a qutrit Bell state. *)
VerificationTest[
    With[{q = QuantumState["PhiPlus", 3]}, {q["Qudits"], q["Dimension"], Normal @ q["StateVector"]}],
    {2, 9, {1/Sqrt[2], 0, 0, 1/Sqrt[2], 0, 0, 0, 0, 0}},
    {QuantumState::padded},
    TestID -> "NamedTail-four-amplitudes-against-three-levels"
]

(* A counted tail widens the basis one qudit per count, and the label widens
   with it, so the name is worn by every factor while the amplitudes fill only
   the first. Pinned because it is silent, normalized and a ket, and so invisible
   to every other check here. *)
VerificationTest[
    With[{q = QuantumState["Plus", 2, 3]},
        {q["Label"], Normal @ q["StateVector"], q == QuantumState["Plus"[3]]}],
    {"+"^CircleTimes[3], {1/Sqrt[2], 1/Sqrt[2], 0, 0, 0, 0, 0, 0}, False},
    {},
    TestID -> "NamedTail-counted-tail-widens-the-label-not-the-amplitudes"
]

(* Each one-qudit name is an eigenvector of a Pauli axis, and each Bell name is
   maximally entangled; a tail that preserves the structure preserves both. Only
   "Plus" and "PhiPlus" were spot-checked elsewhere. *)
VerificationTest[
    QuantumState[#, QuantumBasis[2]]["BlochCartesianCoordinates"] & /@ {"0", "1", "Plus", "Minus", "Left", "Right"},
    {{0, 0, 1}, {0, 0, -1}, {1, 0, 0}, {-1, 0, 0}, {0, -1, 0}, {0, 1, 0}},
    {},
    TestID -> "NamedTail-one-qudit-names-keep-their-axis"
]

VerificationTest[
    QuantumPartialTrace[QuantumState[#, QuantumBasis[2]], {1}]["VonNeumannEntropy"] & /@
        {"PhiPlus", "PhiMinus", "PsiPlus", "PsiMinus"},
    ConstantArray[Quantity[1, "Bits"], 4],
    {},
    TestID -> "NamedTail-all-four-Bell-names-keep-one-bit"
]

(* The two routes into "0" agree on a well-formed tail whatever its dimension,
   which is the property worth pinning; they are not interchangeable on a
   malformed one, since only the spelled-out name reaches the guard. *)
VerificationTest[
    Union @ Map[
        QuantumState["Zero", #] === QuantumState["0", #] &,
        {QuantumBasis[2], "PauliBasis", 2, 3, 4, QuantumBasis[2, 2], QuantumBasis["PauliX"], QuantumBasis["Fourier", 3]}
    ],
    {True},
    {},
    TestID -> "NamedTail-routes-agree-across-dimensions"
]

VerificationTest[
    {
        MatchQ[QuantumState["Zero", QuantumBasis[QuditBasis[2], QuditBasis[3]]], _Failure],
        QuantumState["0", QuantumBasis[QuditBasis[2], QuditBasis[3]]]["InputQudits"],
        QuantumState["Zero", QuantumBasis[3]]["Dimension"]
    },
    {True, 1, 3},
    {QuantumState::invalidArgs},
    TestID -> "NamedTail-routes-diverge-on-a-malformed-tail"
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

(* The registry the invalidName guard tests names against. *)
stateNames = Wolfram`QuantumFramework`PackageScope`$QuantumStateNames;

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

(* A bare unrecognized name used to fall through to the scalar-amplitude rule,
   which wrapped it as a single amplitude and let the default qubit basis pad it,
   giving the two-dimensional state {"NotAnActualName", 0} whose ["Norm"] was
   Abs["NotAnActualName"], with no message at all. It has to fail the same way
   the call form does. *)
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

(* Near misses are what a user actually types: truncations, wrong case, names
   belonging to a sibling constructor, a property name, and strings that look
   like the shorthands without being one. Repeating one message symbol past the
   third occurrence appends General::stop, so the expected list names it once. *)
VerificationTest[
    Select[
        {
            "Bel", "GHZZ", "ghz", "plus", "bell",
            "Hadamard", "Fourier", "PauliBasis",
            "VonNeumannEntropy", "AllProperties",
            "3.5", "0x", "+-i", "-I", "2i", " GHZ", "Left ", "l", "r"
        },
        ! FailureQ @ QuantumState[#] &
    ],
    {},
    {QuantumState::invalidName, QuantumState::invalidName, QuantumState::invalidName, General::stop},
    TestID -> "InvalidName-near-misses"
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

(* Registered names stay untouched, and stay silent. The norm of a single name is
   left to the registry sweeps below, which subsume it: a Failure's ["Norm"] is
   Missing and so would show up there. What those sweeps do not see is the shape
   of the state, so that is what this pins. *)
VerificationTest[
    {QuantumState["Bell"]["Dimension"], QuantumState["Register"[3]]["Qudits"]},
    {4, 3},
    {},
    TestID -> "ValidName-unaffected-Bell-Register"
]

(* Every registered name builds a state rather than reaching the guard. Select
   rather than Union so a failure names the offender. *)
VerificationTest[
    Select[stateNames, FailureQ @ QuantumState[#] &],
    {},
    {},
    TestID -> "ValidName-sweep-registry-none-fail"
]

(* The sweep above cannot catch a name the guard swallows, since the guard fires
   only on non-members: it is the registry that defines the exemption. The drift
   that IS reachable is a name implemented as a "Name"[...] rule but never added
   to $QuantumStateNames, whose bare form the guard then rejects while its call
   form works. Both directions have to stay empty. *)
VerificationTest[
    With[{
        implemented = Union @ Cases[
            Cases[DownValues[QuantumState], HoldPattern[QuantumState[h_[___], ___]] :> HoldForm[h], Infinity],
            _String, Infinity
        ],
        registered = stateNames
    },
        {Complement[implemented, registered], Complement[registered, implemented]}
    ],
    {{}, {}},
    {},
    TestID -> "ValidName-implemented-and-registered-agree"
]

(* And each name builds a unit vector exactly, not merely a non-Failure. The only
   two exceptions are irreducible: both Random* names draw machine-precision
   numbers, so no exact answer exists to compare against. Naming them keeps a
   third name that silently goes numeric visible. *)
VerificationTest[
    Select[stateNames, QuantumState[#]["Norm"] =!= 1 &],
    {"RandomPure", "RandomMixed"},
    {},
    TestID -> "ValidName-sweep-registry-exactly-normalized"
]

(* The two random ones are still normalized to machine precision. *)
VerificationTest[
    Max @ Abs[(QuantumState[#]["Norm"] & /@ {"RandomPure", "RandomMixed"}) - 1] < 10^-15,
    True,
    {},
    TestID -> "ValidName-random-names-are-float-not-wrong"
]

(* The Werner state is the one registry entry that is a family rather than a
   single state, so it is where the constructor can be held symbolically: p stays
   free and the density matrix stays exact. Trace is 1 for every p, and
   Hermiticity holds exactly when p is real, as a mixing weight must be. The
   purity is the closed form 1 - 2p + 4p^2/3, whose three landmarks pin the
   family: 1 at p = 0, the pure singlet; its minimum 1/4 = 1/d at p = 3/4, the
   maximally mixed state; and 1/3 at p = 1. *)
VerificationTest[
    With[{rho = Normal @ QuantumState["Werner"[\[FormalP]]]["DensityMatrix"]},
        {
            Simplify @ Tr[rho],
            Assuming[\[FormalP] \[Element] Reals, Simplify[rho - ConjugateTranspose[rho]]],
            Simplify @ Tr[rho . rho],
            Simplify[Tr[rho . rho] /. \[FormalP] -> {0, 3/4, 1}],
            Minimize[{1 - 2 \[FormalP] + 4 \[FormalP]^2/3, 0 <= \[FormalP] <= 1}, \[FormalP]]
        }
    ],
    {
        1,
        ConstantArray[0, {4, 4}],
        1 - 2 \[FormalP] + 4 \[FormalP]^2/3,
        {1, 1/4, 1/3},
        {1/4, {\[FormalP] -> 3/4}}
    },
    {},
    TestID -> "ValidName-Werner-is-symbolic-in-p"
]

(* The endpoints as objects, not just as purities. *)
VerificationTest[
    {QuantumState["Werner"[0]]["PureStateQ"], QuantumState["Werner"[1]]["PureStateQ"]},
    {True, False},
    {},
    TestID -> "ValidName-Werner-endpoints"
]

(* An independent reading of the same invariant: "Norm" is the constructor's own
   property, so on its own it cannot catch a wrong "Norm". Total probability is no
   use as the second reading, being normalized by construction and so identically
   1 whatever the state (QuantumState[{3,4}] has norm 5 and total probability 1).
   The stored data is the honest cross-check: the Euclidean length of the state
   vector for a ket, the trace of the density matrix for a mixture. *)
VerificationTest[
    Max @ Map[
        Function[name,
            With[{q = QuantumState[name]},
                Abs[q["Norm"] - If[q["PureStateQ"],
                    Norm @ Flatten @ Normal @ q["StateVector"],
                    Sqrt @ Tr @ Normal @ q["DensityMatrix"]
                ]]
            ]
        ],
        stateNames
    ] < 10^-15,
    True,
    {},
    TestID -> "ValidName-norm-agrees-with-stored-data"
]

(* The mixtures in the registry are density matrices, so normalization is only one
   of their invariants: each must also be Hermitian and positive semidefinite with
   unit trace, or it is not a state at all. *)
VerificationTest[
    Select[
        Select[stateNames, ! QuantumState[#]["PureStateQ"] &],
        With[{m = Normal @ QuantumState[#]["DensityMatrix"]},
            ! (HermitianMatrixQ[m] && PositiveSemidefiniteMatrixQ[m] && Chop[Tr[m] - 1] === 0)
        ] &
    ],
    {},
    {},
    TestID -> "ValidName-mixtures-are-density-matrices"
]

(* Names x argument tails, the two-dimensional space the guard actually ranges
   over: no combination reaches the guard. Norm is deliberately not asserted here,
   since a tail is not norm-preserving in general (see below). *)
VerificationTest[
    Select[
        Tuples[{
            stateNames,
            {{}, {QuantumBasis[2]}, {"PauliBasis"}, {QuantumBasis[2], "Label" -> "z"}}
        }],
        FailureQ @ QuantumState[First[#], Sequence @@ Last[#]] &
    ],
    {},
    {},
    TestID -> "ValidName-sweep-registry-times-tails"
]

(* A tail is NOT norm-preserving in general, so the sweeps above must not be read
   as a claim about arbitrary tails. A name whose dimension exceeds the basis is
   TRUNCATED, not padded, and silently: the 3-qubit Dicke state carries 3
   amplitudes of 1/Sqrt[3] across 8 components, only the two in positions 2 and 3
   survive into the 4-dimensional Pauli basis, and the norm drops to Sqrt[2/3]
   with no message. Against this particular 4-dimensional tail only "Dicke" and
   "Graph" change dimension at all; the count is a property of the tail, not of
   the registry, and the qutrit case below shows it moving. *)
VerificationTest[
    {
        Normal @ QuantumState["Dicke"]["StateVector"],
        Normal @ QuantumState["Dicke", "PauliBasis"]["StateVector"],
        QuantumState["Dicke", "PauliBasis"]["Norm"],
        Select[
            stateNames,
            QuantumState[#, QuantumBasis[2]]["Dimension"] =!= QuantumState[#]["Dimension"] &
        ]
    },
    {
        {0, 1/Sqrt[3], 1/Sqrt[3], 0, 1/Sqrt[3], 0, 0, 0},
        {0, 1/Sqrt[3], 1/Sqrt[3], 0},
        Sqrt[2/3],
        {"Dicke", "Graph"}
    },
    {},
    TestID -> "ValidName-tail-truncation-is-not-norm-preserving"
]

(* The sweeps above range over qubit tails, where a name's own dimension is a
   power of the basis dimension. A qutrit tail breaks that: the four two-qubit
   Bell names hold 4 amplitudes, which no power of 3 matches, so they are padded
   out to 9 and cease to be the state they name. QuantumState::padded says so,
   which is why these cannot join the silent sweeps. *)
VerificationTest[
    Select[stateNames, Check[QuantumState[#, 3]; False, True] &],
    {"PsiPlus", "PsiMinus", "PhiPlus", "PhiMinus"},
    {QuantumState::padded, QuantumState::padded, QuantumState::padded, General::stop},
    TestID -> "ValidName-qutrit-tail-pads-the-two-qubit-names"
]

(* The zero-subsystem boundary of the same dispatch: an empty register is the
   one-dimensional trivial state, and the arm must not intercept it. *)
VerificationTest[
    {#["Dimension"], #["Norm"]} & /@ {
        QuantumState["Register"[0]], QuantumState["GHZ"[0]], QuantumState["W"[0]],
        QuantumState["RandomPure"[0]], QuantumState["UniformSuperposition"[0]],
        QuantumState["RandomMixed"[0]], QuantumState["UniformMixture"[0]]
    },
    ConstantArray[{1, 1}, 7],
    {},
    TestID -> "ValidName-zero-subsystem-boundary"
]

(* An unrecognized name is unrecognized through the dagger route too, which sits
   above the arm and re-enters dispatch: without a Confirm it returned the
   Missing from a property query on the Failure instead of the Failure. *)
VerificationTest[
    {FailureQ @ QuantumState[SuperDagger["NotAnActualName"]], QuantumState[SuperDagger["Plus"]]["Norm"]},
    {True, 1},
    {QuantumState::invalidName},
    TestID -> "InvalidName-through-SuperDagger"
]

(* The independent reference for the family claim: each constructor that owns an
   invalidName message reports an unknown bare name with its OWN, not merely with
   some message mentioning the string. The messages are asserted here rather than
   captured and discarded, so the expected list carries all five. QuditBasis wraps
   its Failure in a ConfirmationFailed from the basis it builds, which is why the
   message and not the Failure tag is what the family shares. *)
VerificationTest[
    Head /@ {
        QuantumState["NotAnActualName"], QuantumOperator["NotAnActualName"],
        QuantumChannel["NotAnActualName"], QuantumCircuitOperator["NotAnActualName"],
        QuditBasis["NotAnActualName"]
    },
    ConstantArray[Failure, 5],
    {
        QuantumState::invalidName, QuantumOperator::invalidName,
        QuantumChannel::invalidName, QuantumCircuitOperator::invalidName,
        QuditBasis::invalidName
    },
    TestID -> "InvalidName-family-agreement"
]

(* The string forms that reach a state without appearing in $QuantumStateNames:
   the zero-qudit empty state, basis-element digit strings (in any qudit
   dimension), computational and Pauli-eigenstate letter sequences, and the
   "+i"/"-i" aliases. The arm-last file order asserted in NamedStates.m is what
   keeps all but "" reachable. *)
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

(* Two spellings of the zero-qudit state that disagree on the amplitude, which is
   why the norm sweeps above exclude "": QuantumState[""] is built from an empty
   amplitude list and so has norm 0, while "Register"[0] is built from {1} and is
   normalized. Both carry zero qudits, so it is the amplitude and not the qudit
   count that decides. *)
VerificationTest[
    {
        {QuantumState[""]["Qudits"], QuantumState[""]["Norm"], Normal @ QuantumState[""]["StateVector"]},
        {QuantumState["Register"[0]]["Qudits"], QuantumState["Register"[0]]["Norm"], Normal @ QuantumState["Register"[0]]["StateVector"]}
    },
    {{0, 0, {0}}, {0, 1, {1}}},
    {},
    TestID -> "Shorthand-empty-state-has-zero-norm"
]

(* Same forms carrying an argument tail, which is where they compete with the arm
   at equal arity. Nothing else in the suite gives the letter-sequence rule a
   tail. *)
VerificationTest[
    {
        QuantumState["0101", QuantumBasis[16]]["Qudits"],
        QuantumState["12", 3]["Qudits"],
        QuantumState["+0L", QuantumBasis[8]]["Norm"]
    },
    {4, 2, 1},
    {},
    TestID -> "Shorthand-strings-with-a-tail"
]

(* Any sequence over the eigenstate letters is one qudit per character, so a
   two-letter mix is a two-qubit product state and not a near-miss name. It is the
   tensor product of its characters as a state; the two differ only in label,
   since the sequence names itself where the tensor product spells out its
   factors. *)
VerificationTest[
    {
        QuantumState["L+"]["Qudits"],
        QuantumState["L+"]["Norm"],
        Normal @ QuantumState["L+"]["StateVector"],
        Normal @ QuantumState["L+"]["StateVector"] ===
            Normal @ QuantumTensorProduct[QuantumState["L"], QuantumState["+"]]["StateVector"],
        QuantumState["L+"]["Label"]
    },
    {2, 1, {1/2, 1/2, -I/2, -I/2}, True, "L+"},
    {},
    TestID -> "Shorthand-mixed-letter-pair-is-a-product-state"
]

(* Each one-character or "+i"/"-i" spelling and the name it stands for. *)
aliasSpellings = {{"+", "Plus"}, {"-", "Minus"}, {"L", "Left"}, {"R", "Right"}, {"+i", "Right"}, {"-i", "Left"}};

(* Bare, every alias is the same object as the name it spells. Select rather than
   Union so a failure names the spelling that broke. *)
VerificationTest[
    Select[aliasSpellings, QuantumState[First[#]] =!= QuantumState[Last[#]] &],
    {},
    {},
    TestID -> "Shorthand-alias-agrees-with-name-bare"
]

(* And with a tail too, which is where the one-character spellings used to part
   company: they reached the letter-sequence rule and came back labelled "I", the
   identity, so QuantumState["+", QuantumBasis[2]] was the plus state wearing the
   wrong name. Gating that rule on a length above one makes every spelling the
   same object, label included. *)
VerificationTest[
    Cases[
        Map[QuantumState[#, QuantumBasis[2]] &, aliasSpellings, {2}],
        {a_, b_} /; a =!= b || a["Norm"] =!= 1
    ],
    {},
    {},
    TestID -> "Shorthand-alias-agrees-with-name-under-a-tail"
]

(* The label is what the letter-sequence route used to overwrite, so it is worth
   naming rather than folding into the object comparison above. *)
VerificationTest[
    Map[QuantumState[First[#], QuantumBasis[2]]["Label"] &, aliasSpellings],
    {"+", "-", "L", "R", "R", "L"},
    {},
    TestID -> "Shorthand-alias-keeps-the-state-label-under-a-tail"
]

(* A genuine sequence keeps its own name under a tail for the same reason, the
   rule now passing the string down as a default label the way the digit rule
   two lines above it always did. Untailed it was already right; with a basis it
   used to come back as that basis's label. *)
VerificationTest[
    {
        QuantumState["+0L", QuantumBasis[8]]["Label"],
        QuantumState["0101", QuantumBasis[16]]["Label"],
        QuantumState["+0L", QuantumBasis[8]]["Norm"]
    },
    {"+0L", "0101", 1},
    {},
    TestID -> "Shorthand-sequence-keeps-its-label-under-a-tail"
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

(* "" is a good name carrying no levels, so a basis is a bad ARGUMENT to it, not
   a bad name: it reports the way the other named states report a tail they cannot
   carry, rather than contradicting the tests above that build QuantumState[""]. A
   tail used to slip past the alias rules entirely and normalize Abs[""]. *)
VerificationTest[
    QuantumState["", QuantumBasis[2]],
    Failure["InvalidArguments", _],
    {QuantumState::invalidArgs},
    SameTest -> MatchQ,
    TestID -> "Empty-state-with-tail-is-a-bad-argument"
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
