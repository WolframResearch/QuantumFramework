(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Transversal.wlt

   Transversal gates: which Cliffords a code performs by acting qubit by qubit.

   Got26 sec. 11.2 gives the criterion -- U applied to every qubit separately is
   a gate gadget exactly when it maps the stabilizer to itself -- and sec. 11.3
   works it out for the 7-qubit code.  The load-bearing test here is
   QEC-Transversal-S-on-the-7-qubit-code-is-the-logical-Sdg, because it is the
   one a sign-blind conjugation gets wrong: Y^7 = -Ybar (eq. 11.22), so the
   transversal S performs the logical S^dagger and not the logical S.  Its
   companion is QEC-Transversal-matrices-match-the-engine, which pins the
   convention every sign below is read in.

   The negative cases carry as much as the positive ones.  On the 5-qubit code
   neither H nor S nor the transversal CNOT is a gate gadget at all, and yet the
   code is not gateless: the cyclic Clifford SH is transversal on it, and exactly
   twelve of the twenty-four one-qubit Cliffords are.  A scan that reported
   "none" would be as wrong as one that reported "all".
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";

qecAction     = Symbol[qecScope <> "transversalAction"];
qecMatrix     = Symbol[qecScope <> "transversalMatrix"];
qecConjugate  = Symbol[qecScope <> "transversalConjugate"];
qecImages     = Symbol[qecScope <> "transversalGeneratorImages"];
qecEngineGates = Symbol[qecScope <> "instructionEngineGates"];
qecEncoding   = Symbol[qecScope <> "codeEncodingGates"];
qecApplyGates = Symbol[qecScope <> "applyGates"];
qecCodeData   = Symbol[qecScope <> "codeData"];
qecOneQubitOps = Symbol[qecScope <> "$circuitOneQubitOps"];

steane = QECCode["SteaneCode"];
five = QECCode["5QubitCode"];
four = QECCode["DistanceTwo", 4];

(* Two matrices are the same gate when they differ by a global phase, which is all
   a gate gadget ever determines. *)
qecSamePhaseQ[m1_, m2_] := With[{p = Simplify[m1 . ConjugateTranspose[m2]]},
    Simplify[p - p[[1, 1]] IdentityMatrix[Length[p]]] === ConstantArray[0, Dimensions[p]]
]

(* The action table read as strings: the image of X, of Z, of Y. *)
qecActionStrings[gate_] := QECPauliString /@ Lookup[qecAction[gate, 1], {{1, 0}, {0, 1}, {1, 1}}]


(* ============================================================================
   The action tables, and the convention they are read in
   ============================================================================ *)

(* Not typed in: computed from the gate's matrix.  These are the engine's matrices,
   compared here so that a drift in either convention fails loudly rather than
   quietly changing every sign in the file. *)
VerificationTest[
    AllTrue[
        {"H", "S", "V", "X", "Y", "Z", "CNOT", "CZ", "SWAP"},
        qecSamePhaseQ[qecMatrix[#], Normal[QuantumOperator[#]["MatrixRepresentation"]]] &
    ],
    True,
    TestID -> "QEC-Transversal-matrices-match-the-engine"
]

(* H: X <-> Z and Y -> -Y. *)
VerificationTest[
    qecActionStrings["H"],
    {"Z", "X", "-Y"},
    TestID -> "QEC-Transversal-H-exchanges-X-and-Z"
]

(* S: X -> Y, Z -> Z, Y -> -X.  The sign on Y is the one that matters later. *)
VerificationTest[
    qecActionStrings["S"],
    {"Y", "Z", "-X"},
    TestID -> "QEC-Transversal-S-sends-X-to-Y"
]

(* Sdg is S run backwards, and the two differ ONLY in phases -- which is exactly
   why the frame propagator, which drops phases, cannot tell them apart. *)
VerificationTest[
    {qecActionStrings["Sdg"], qecActionStrings[{"S", "Sdg"}]},
    {{"-Y", "Z", "X"}, {"X", "Z", "Y"}},
    TestID -> "QEC-Transversal-Sdg-undoes-S"
]

(* V is Sqrt[X] in the engine's convention, so it fixes X and sends Z to -Y. *)
VerificationTest[
    qecActionStrings["V"],
    {"X", "-Y", "Z"},
    TestID -> "QEC-Transversal-V-is-the-square-root-of-X"
]

(* A list is a sequence in time order: H then S is not S then H. *)
VerificationTest[
    {qecActionStrings[{"H", "S"}], qecActionStrings[{"S", "H"}]},
    {{"Z", "Y", "X"}, {"-Y", "X", "-Z"}},
    TestID -> "QEC-Transversal-a-list-is-a-sequence-in-time-order"
]


(* ============================================================================
   The 7-qubit code (sec. 11.3)
   ============================================================================ *)

(* Eq. 11.22, on its own, because everything below turns on it: Y^7 is MINUS the
   logical Y, since Xbar Zbar = i Y^7 and Ybar = i Xbar Zbar. *)
VerificationTest[
    {QECPauliString[QECPauliProduct["XXXXXXX", "ZZZZZZZ"]],
     QECPauliString[QECPauliProduct[{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, QECPauliProduct["XXXXXXX", "ZZZZZZZ"]]]},
    {"iYYYYYYY", "-YYYYYYY"},
    TestID -> "QEC-Transversal-Y-to-the-seven-is-minus-the-logical-Y"
]

(* Eq. 11.18-11.19: the transversal Hadamard is the logical Hadamard. *)
VerificationTest[
    With[{g = QECTransversalGate[steane, "H"]},
        {g["TransversalQ"], g["LogicalGate"], Values[g["LogicalAction"]]}
    ],
    {True, "H", {"Z", "X"}},
    TestID -> "QEC-Transversal-H-on-the-7-qubit-code-is-the-logical-H"
]

(* THE test, eq. 11.20-11.22.  Xbar -> -Ybar and Zbar -> Zbar is the action of
   S^dagger, not of S: the transversal S performs the INVERSE of the gate it is
   made of.  A conjugation that dropped signs would report "S" here. *)
VerificationTest[
    With[{g = QECTransversalGate[steane, "S"]},
        {g["TransversalQ"], g["LogicalGate"], Values[g["LogicalAction"]]}
    ],
    {True, "Sdg", {"-Y", "Z"}},
    TestID -> "QEC-Transversal-S-on-the-7-qubit-code-is-the-logical-Sdg"
]

(* Being a gate gadget is the other half: every generator has to come back a
   stabilizer WITH ITS SIGN, or the gate takes the state out of the code space.
   All six do, because every element of this stabilizer has weight 4. *)
VerificationTest[
    With[{g = QECTransversalGate[steane, "S"]},
        Normal[g["StabilizerImages"][All, "StabilizerQ"]]
    ],
    ConstantArray[True, 6],
    TestID -> "QEC-Transversal-S-keeps-every-generator-a-stabilizer"
]

(* Sec. 11.3's conclusion: U^7 is a gadget for EVERY one-qubit Clifford, all
   twenty-four of them, which is what makes this code the chapter's favourite. *)
VerificationTest[
    Length[QECTransversalGate[steane, All]],
    24,
    TestID -> "QEC-Transversal-the-7-qubit-code-admits-every-one-qubit-Clifford"
]

(* The book's own explanation of the S^dagger: the transversal U performs the
   logical U*, the complex conjugate.  Checked for every one-qubit gate, up to the
   global phase a gate gadget never fixes. *)
VerificationTest[
    AllTrue[
        {"I", "X", "Y", "Z", "H", "S", "Sdg", "V", "Vdg"},
        qecSamePhaseQ[
            qecMatrix[QECTransversalGate[steane, #]["LogicalGate"]],
            Conjugate[qecMatrix[#]]
        ] &
    ],
    True,
    TestID -> "QEC-Transversal-the-logical-gate-is-the-complex-conjugate"
]

(* Eq. 11.23-11.26, now with the signs the phase-blind version could not see. *)
VerificationTest[
    With[{g = QECTransversalGate[steane, "CNOT"]},
        {g["TransversalQ"], g["LogicalGate"], Values[g["LogicalAction"]]}
    ],
    {True, "CNOT", {"XX", "IX", "ZI", "ZZ"}},
    TestID -> "QEC-Transversal-CNOT-between-two-blocks-is-the-logical-CNOT"
]

(* And the other two two-block gates, which come free with the same symmetry. *)
VerificationTest[
    {QECTransversalGate[steane, "CZ"]["LogicalGate"],
     QECTransversalGate[steane, "SWAP"]["LogicalGate"]},
    {"CZ", "SWAP"},
    TestID -> "QEC-Transversal-CZ-and-SWAP-are-also-transversal-on-the-7-qubit-code"
]

(* The engine's opinion, on the state rather than on the algebra: transversal H
   applied to the encoded |0> leaves the state in the code space -- every generator
   still has expectation 1 -- and turns Zbar = +1 into Xbar = +1, which is the
   logical Hadamard acting on the encoded qubit. *)
VerificationTest[
    Module[{data = qecCodeData[steane], state},
        state = qecApplyGates[
            qecApplyGates[PauliStabilizer[7], qecEncoding[data]],
            qecEngineGates[QECTransversalGate[steane, "H"]["Instructions"]]
        ];
        {AllTrue[steane["Generators"], state["Expectation", #] === 1 &],
         state["Expectation", First[steane["LogicalOperators"]["X"]]] === 1}
    ],
    {True, True},
    TestID -> "QEC-Transversal-H-takes-the-encoded-zero-to-the-encoded-plus"
]


(* ============================================================================
   The 5-qubit code: different symmetry, different gates
   ============================================================================ *)

(* Neither H nor S is a gadget here, and a gate that is not one has no logical
   action to report -- the images are still there to look at. *)
VerificationTest[
    With[{g = QECTransversalGate[five, "H"]},
        {g["TransversalQ"], g["LogicalAction"], g["LogicalGate"]}
    ],
    {False, Missing["NotTransversal"], Missing["NotTransversal"]},
    TestID -> "QEC-Transversal-H-is-not-a-gadget-on-the-5-qubit-code"
]

(* Sec. 11.4: transversal CNOT forces the generators to split into X-only and
   Z-only, which is what it means to be CSS.  Register.wl found this modulo sign;
   this is the same conclusion with the phases carried. *)
VerificationTest[
    {QECTransversalGate[five, "CNOT"]["TransversalQ"], five["CSSQ"]},
    {False, False},
    TestID -> "QEC-Transversal-CNOT-needs-a-CSS-code"
]

(* And yet the code is not gateless: the cyclic Clifford X -> Y -> Z -> X is
   transversal on it, and exactly half the one-qubit Cliffords are. *)
VerificationTest[
    With[{scan = QECTransversalGate[five, All]},
        {QECTransversalGate[five, {"S", "H"}]["TransversalQ"],
         Length[scan],
         SubsetQ[Keys[scan], {"I", "X", "Y", "Z"}]}
    ],
    {True, 12, True},
    TestID -> "QEC-Transversal-the-5-qubit-code-admits-the-cyclic-Clifford"
]


(* ============================================================================
   More than one logical qubit in a block
   ============================================================================ *)

(* On the [[4,2,2]] code of sec. 11.2 the transversal H acts on BOTH logical
   qubits at once -- it is a Hadamard on each followed by a swap of the two -- so
   there is no one-qubit gate to name it with, and the action is the answer. *)
VerificationTest[
    With[{g = QECTransversalGate[four, "H"]},
        {g["TransversalQ"], Values[g["LogicalAction"]], g["LogicalGate"]}
    ],
    {True, {"IZ", "ZI", "IX", "XI"}, Missing["NotNamed"]},
    TestID -> "QEC-Transversal-H-on-the-4-qubit-code-swaps-the-logical-pair"
]


(* ============================================================================
   The circuit
   ============================================================================ *)

(* Qubit by qubit, which is the definition; a sequence is gate-major, so its depth
   is the length of the sequence and not the number of qubits. *)
VerificationTest[
    With[{g = QECTransversalGate[steane, {"H", "S"}]},
        {Take[g["Instructions"], 2], Length[g["Instructions"]], g["Depth"], g["Qubits"]}
    ],
    {{{"H", 1}, {"H", 2}}, 14, 2, 7},
    TestID -> "QEC-Transversal-instructions-are-qubit-by-qubit"
]

(* A two-block gate is emitted on the register layout of Register.wl: qubit q of
   block 1 with qubit q of block 2. *)
VerificationTest[
    {QECTransversalGate[steane, "CNOT"]["Instructions"],
     QECTransversalGate[steane, "I"]["Instructions"]},
    {Table[{"CNOT", q, q + 7}, {q, 7}], {}},
    TestID -> "QEC-Transversal-two-block-instructions-use-the-register-layout"
]

(* Sdg had to join the circuit vocabulary for this file, and the engine has no
   inverse-S of its own, so it is emitted as three S -- the same treatment Vdg
   already had. *)
VerificationTest[
    {MemberQ[qecOneQubitOps, "Sdg"], qecEngineGates[{{"Sdg", 3}}]},
    {True, {"S" -> 3, "S" -> 3, "S" -> 3}},
    TestID -> "QEC-Transversal-Sdg-is-three-S-for-the-engine"
]


(* ============================================================================
   Refusals
   ============================================================================ *)

VerificationTest[
    QECTransversalGate[steane, "T"],
    $Failed,
    {QECTransversalGate::gate},
    TestID -> "QEC-Transversal-refuses-a-gate-it-cannot-build"
]

VerificationTest[
    QECTransversalGate[steane, "H"]["Nonsense"],
    Missing["NotFound", "Nonsense"],
    {QECTransversalGate::noprop},
    TestID -> "QEC-Transversal-unknown-property-is-refused"
]
