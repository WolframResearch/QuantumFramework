(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Operational.wlt

   The code presented as the objects quantum theory uses for error correction:
   the encoding isometry, the codewords, the code-space projector, the syndrome
   instrument, the recovery channels, the effective logical channel, and the
   Knill-Laflamme matrix.

   Two tests carry the file.

   QEC-Operational-logical-channel-agrees-with-the-rate: the logical channel is
   built from the coset engine, not by composing 4^n x 4^n superoperators, and
   the claim that nothing is lost is exactly 1 - q_I == QECLogicalErrorRate,
   symbolically, on the Steane code.

   QEC-Operational-encoder-speaks-the-code's-logical-basis: V^dagger Zbar_j V = Z_j
   and V^dagger Xbar_j V = X_j for the code's own reported logical operators, on
   every small named code.  The encoding circuit alone fails this on the 5-qubit
   code (its |0> output is the -1 eigenstate of Zbar); the isometry is
   canonicalized, and this test is what keeps it so.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

(* The row functions and helpers that the public QECPauli and QECCode["CSS", ...]
   stand on, in PackageScope since step 5 of the API redesign. *)
qecPauliVector = Symbol["Wolfram`QuantumFramework`QEC`PackageScope`pauliVector"];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";

qecPauliMatrix = Symbol[qecScope <> "pauliRowMatrix"];
qecProjectors  = Symbol[qecScope <> "codeSyndromeProjectors"];

bitFlip = QECCode["BitFlipCode"];
five = QECCode["5QubitCode"];
steane = QECCode["SteaneCode"];
shor = QECCode["ShorCode"];
twoLogical = QECCode["DistanceTwo", 4];

opMatrix[p_String] := Normal[qecPauliMatrix[qecPauliVector[p]]]
kron[l_List] := Fold[KroneckerProduct, l]
onLogical[k_, j_, s_] := kron @ ReplacePart[ConstantArray[IdentityMatrix[2], k], j -> s]
zeroQ[m_] := AllTrue[Flatten[Chop[Simplify[m]]], # === 0 &]
encoderMatrix[c_] := Normal[c["Encoder"]["MatrixRepresentation"]]

(* V^dagger P V for each logical P of the given kind, against the one-qubit Pauli s
   on logical qubit j. *)
logicalActionQ[c_, kind_String, s_] := With[{v = encoderMatrix[c], k = c["LogicalQubits"]},
    And @@ Table[
        zeroQ[ConjugateTranspose[v] . opMatrix[c[kind][[j]]] . v - onLogical[k, j, s]],
        {j, k}
    ]
]


(* ============================================================================
   The isometry, the codewords, the projector
   ============================================================================ *)

VerificationTest[
    Table[With[{v = encoderMatrix[c]}, zeroQ[ConjugateTranspose[v] . v - IdentityMatrix[2^c["LogicalQubits"]]]],
        {c, {bitFlip, five, steane, twoLogical}}],
    {True, True, True, True},
    TestID -> "QEC-Operational-encoder-is-an-isometry"
]

VerificationTest[
    Table[{logicalActionQ[c, "LogicalZ", PauliMatrix[3]], logicalActionQ[c, "LogicalX", PauliMatrix[1]]},
        {c, {bitFlip, five, steane, shor, twoLogical}}],
    ConstantArray[{True, True}, 5],
    TestID -> "QEC-Operational-encoder-speaks-the-code's-logical-basis"
]

(* The same projector reached two ways: from the isometry, and as the product of
   (1 + g)/2 over the generators, which never sees the encoder. *)
VerificationTest[
    Table[
        zeroQ[Normal[c["Codespace"]["MatrixRepresentation"]] - Normal[First[qecProjectors[First[c]]]]],
        {c, {bitFlip, five, steane}}
    ],
    {True, True, True},
    TestID -> "QEC-Operational-codespace-is-the-stabilizer-projector"
]

VerificationTest[
    {Length[#], Union[Head /@ #]} & /@ {bitFlip["Codewords"], twoLogical["Codewords"]},
    {{2, {QuantumState}}, {4, {QuantumState}}},
    TestID -> "QEC-Operational-one-codeword-per-logical-basis-state"
]

(* |0_L> is the +1 eigenstate of Zbar, |1_L> the -1: the order is the logical basis
   order, which is what makes "Codewords" and "Encoder" the same object. *)
VerificationTest[
    With[{z = opMatrix[First[five["LogicalZ"]]]},
        Simplify[Conjugate[#] . z . #] & /@ (Normal[#["StateVector"]] & /@ five["Codewords"])
    ],
    {1, -1},
    TestID -> "QEC-Operational-codewords-are-in-Zbar-order"
]


(* ============================================================================
   Syndrome instrument and recovery
   ============================================================================ *)

VerificationTest[
    Head[bitFlip["SyndromeMeasurement"]],
    QuantumMeasurementOperator,
    TestID -> "QEC-Operational-syndrome-measurement-is-a-measurement-operator"
]

(* Outcome i is the syndrome IntegerDigits[i - 1, 2, m]: X on the middle qubit trips
   both checks, {1, 1}, the fourth outcome, with certainty. *)
VerificationTest[
    Values[Normal[bitFlip["SyndromeMeasurement"][QuantumOperator["IXI"][First[bitFlip["Codewords"]]]]["Probabilities"]]],
    {0, 0, 0, 1},
    TestID -> "QEC-Operational-syndrome-outcome-order-matches-Syndrome"
]

VerificationTest[
    With[{psi = QuantumOperator["IXI"][First[bitFlip["Codewords"]]]},
        zeroQ[
            Normal[bitFlip["Recovery", bitFlip["Syndrome", "IXI"]][QuantumState[psi["DensityMatrix"]]]["DensityMatrix"]] -
            Normal[First[bitFlip["Codewords"]]["DensityMatrix"]]
        ]
    ],
    True,
    TestID -> "QEC-Operational-recovery-undoes-a-correctable-error"
]

VerificationTest[
    QECSyndromeCircuit[bitFlip]["SyndromeMeasurement"] === bitFlip["SyndromeMeasurement"],
    True,
    TestID -> "QEC-Operational-syndrome-circuit-states-the-instrument-it-implements"
]

(* The recovery is one Kraus operator on all three wires.  A one-element Kraus list
   with an order is read by the engine as a channel whose environment replaces the
   first wire, so this pins the route that avoids it. *)
VerificationTest[
    bitFlip["Recovery", {1, 0}]["Order"],
    {{1, 2, 3}, {1, 2, 3}},
    TestID -> "QEC-Operational-recovery-acts-on-every-wire"
]


(* ============================================================================
   The logical channel
   ============================================================================ *)

VerificationTest[
    Simplify /@ bitFlip["LogicalPauliProbabilities", QECNoiseModel["BitFlip", qecP]],
    <|"I" -> Simplify[(1 - qecP)^2 (1 + 2 qecP)], "X" -> Simplify[(3 - 2 qecP) qecP^2], "Y" -> 0, "Z" -> 0|>,
    TestID -> "QEC-Operational-bit-flip-logical-channel-is-exact"
]

VerificationTest[
    Simplify[1 - steane["LogicalPauliProbabilities", QECNoiseModel["Depolarizing", qecP]]["I"] -
        QECLogicalErrorRate[steane, QECNoiseModel["Depolarizing", qecP]]],
    0,
    TestID -> "QEC-Operational-logical-channel-agrees-with-the-rate"
]

VerificationTest[
    Simplify[Total[five["LogicalPauliProbabilities", QECNoiseModel["Depolarizing", qecP]]]],
    1,
    TestID -> "QEC-Operational-logical-probabilities-sum-to-one"
]

VerificationTest[
    With[{ch = steane["LogicalChannel", QECNoiseModel["Depolarizing", 1/10]],
          q = steane["LogicalPauliProbabilities", QECNoiseModel["Depolarizing", 1/10]]},
        zeroQ[Normal[ch[QuantumState["0"]]["DensityMatrix"]] - DiagonalMatrix[{q["I"] + q["Z"], q["X"] + q["Y"]}]]
    ],
    True,
    TestID -> "QEC-Operational-logical-channel-acts-as-its-probabilities-say"
]

(* Noise that only ever produces the identity leaves one Kraus operator; the channel
   must still act on the logical wire, not on an environment in its place. *)
VerificationTest[
    bitFlip["LogicalChannel", QECNoiseModel["BitFlip", 0]]["Order"],
    {{1}, {1}},
    TestID -> "QEC-Operational-noiseless-logical-channel-is-the-identity-on-its-wire"
]

VerificationTest[
    Simplify[1 - QECCode["Repetition", 3]["LogicalPauliProbabilities", QuantumChannel["BitFlip"[qecP]]]["I"]],
    Simplify[3 qecP^2 - 2 qecP^3],
    TestID -> "QEC-Operational-logical-channel-accepts-a-QuantumChannel"
]

VerificationTest[
    steane["LogicalChannel", QECNoiseModel["Circuit", 1/100]],
    Missing["NotAvailable", "Circuit"],
    {QECCode::nolevel},
    TestID -> "QEC-Operational-logical-channel-is-code-capacity-only"
]


(* ============================================================================
   Knill-Laflamme
   ============================================================================ *)

VerificationTest[
    bitFlip["KnillLaflammeMatrix", {"III", "XII", "IXI", "IIX"}],
    IdentityMatrix[4],
    TestID -> "QEC-Operational-bit-flip-corrects-single-X"
]

VerificationTest[
    bitFlip["KnillLaflammeMatrix", {"XII", "IXX"}],
    Missing["NotCorrectable", <|"Errors" -> {"XII", "IXX"}, "Product" -> "XXX"|>],
    TestID -> "QEC-Operational-bit-flip-names-the-logical-that-breaks-it"
]

VerificationTest[
    {steane["CorrectableQ", "Weight"[1]], steane["CorrectableQ", "Weight"[2]], five["CorrectableQ", "Weight"[1]]},
    {True, False, True},
    TestID -> "QEC-Operational-correctable-weight-matches-the-distance"
]

(* The symplectic answer against the definition: <W_i| Ea^dagger Eb |W_j> = h_ab delta_ij
   computed densely from the codewords, on a set that includes a stabilizer (ZZI, so
   h is not diagonal). *)
VerificationTest[
    With[{errs = {"III", "XII", "ZZI"}, ws = Normal[#["StateVector"]] & /@ bitFlip["Codewords"]},
        Table[
            Conjugate[ws[[i]]] . ConjugateTranspose[opMatrix[ea]] . opMatrix[eb] . ws[[j]],
            {i, 2}, {j, 2}, {ea, errs}, {eb, errs}
        ] === {{bitFlip["KnillLaflammeMatrix", errs], ConstantArray[0, {3, 3}]},
               {ConstantArray[0, {3, 3}], bitFlip["KnillLaflammeMatrix", errs]}}
    ],
    True,
    TestID -> "QEC-Operational-Knill-Laflamme-matches-its-dense-definition"
]

VerificationTest[
    bitFlip["KnillLaflammeMatrix", {"XX"}],
    $Failed,
    {QECCode::klerrors},
    TestID -> "QEC-Operational-Knill-Laflamme-refuses-a-wrong-size-error"
]


(* ============================================================================
   Noise and channels
   ============================================================================ *)

VerificationTest[
    QECNoiseModel[QuantumChannel["Depolarizing"[qecQ]]]["Probabilities"],
    {1 - 3 qecQ / 4, qecQ / 4, qecQ / 4, qecQ / 4},
    TestID -> "QEC-Operational-noise-from-a-named-channel"
]

VerificationTest[
    QECNoiseModel[QuantumChannel[{Sqrt[1 - qecU - qecW] IdentityMatrix[2], Sqrt[qecU] PauliMatrix[1], Sqrt[qecW] PauliMatrix[3]}]]["Probabilities"],
    {1 - qecU - qecW, qecU, 0, qecW},
    TestID -> "QEC-Operational-noise-from-Kraus-operators"
]

VerificationTest[
    QECNoiseModel[QuantumChannel["AmplitudeDamping"[qecG]]],
    $Failed,
    {QECNoiseModel::notpauli},
    TestID -> "QEC-Operational-non-Pauli-channel-is-refused"
]

(* A general Pauli model had no channel before; now it does, and it reads back. *)
VerificationTest[
    QECNoiseModel[QECNoiseModel[<|"X" -> qecU, "Z" -> qecW|>]["QuantumChannel"]]["Probabilities"],
    {1 - qecU - qecW, qecU, 0, qecW},
    TestID -> "QEC-Operational-general-Pauli-noise-round-trips-through-a-channel"
]

VerificationTest[
    Simplify[QECLogicalErrorRate[bitFlip, QuantumChannel["BitFlip"[qecP]]] -
        QECLogicalErrorRate[bitFlip, QECNoiseModel["BitFlip", qecP]]],
    0,
    TestID -> "QEC-Operational-rate-accepts-a-QuantumChannel"
]


(* ============================================================================
   The dense guard, and the listing
   ============================================================================ *)

VerificationTest[
    Block[{$QECDenseQubitLimit = 5}, steane["Encoder"]],
    $Failed,
    {QECCode::dense},
    TestID -> "QEC-Operational-dense-objects-respect-the-limit"
]

(* The cheap objects are not guarded: the logical channel and the Knill-Laflamme
   matrix never form a 2^n vector. *)
VerificationTest[
    Block[{$QECDenseQubitLimit = 5}, steane["CorrectableQ", "Weight"[1]]],
    True,
    TestID -> "QEC-Operational-Knill-Laflamme-is-not-dense"
]

VerificationTest[
    SubsetQ[steane["Properties"], {"Encoder", "Codewords", "Codespace", "SyndromeMeasurement", "Recovery",
        "LogicalChannel", "LogicalPauliProbabilities", "KnillLaflammeMatrix", "CorrectableQ"}],
    True,
    TestID -> "QEC-Operational-properties-are-listed"
]
