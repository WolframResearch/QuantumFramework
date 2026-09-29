(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Surface.wlt

   Step 5 of the API redesign: the public surface after the renames.

   Three things are pinned.  QECPauli is one object whose properties are the seven
   former verbs, and returns what they returned.  The renamed heads and the
   properties that replace free functions give the same answers under the new
   spelling.  And every retired name still works, returns the same value, and
   warns exactly once per session -- the one-release deprecation the plan
   promises.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";
qecPauliVector = Symbol[qecScope <> "pauliVector"];
qecPauliString = Symbol[qecScope <> "pauliString"];
qecPauliProduct = Symbol[qecScope <> "pauliProduct"];
qecHammingMatrix = Symbol[qecScope <> "classicalHammingMatrix"];

steane = QECCode["SteaneCode"];
five = QECCode["5QubitCode"];


(* ============================================================================
   QECPauli
   ============================================================================ *)

VerificationTest[
    QECPauli["XZZXI"],
    QECPauli[{1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0}],
    TestID -> "QEC-Surface-Pauli-canonical-form-is-the-row"
]

VerificationTest[
    With[{p = QECPauli["-iXY"]},
        {p["String"], p["Vector"], p["Weight"], p["Phase"], p["Qubits"], p["Support"], p["HermitianQ"]}
    ],
    {"-iXY", {1, 1, 0, 1, 3}, 2, 3, 2, {1, 2}, False},
    TestID -> "QEC-Surface-Pauli-properties"
]

VerificationTest[
    {QECPauli[{1, 0, 5}], QECPauli[QECPauli["X"]]},
    {QECPauli[{1, 0, 1}], QECPauli["X"]},
    TestID -> "QEC-Surface-Pauli-normalizes-phase-and-is-idempotent"
]

VerificationTest[
    {(QECPauli["X"] ** QECPauli["Z"])["String"], QECPauli["XX"]["Product", "ZZ", "YY"]["String"]},
    {"-iY", "-II"},
    TestID -> "QEC-Surface-Pauli-product-carries-the-Z4-phase"
]

VerificationTest[
    {QECPauli["XZZXI"]["CommuteQ", "IXZZX"], QECPauli["X"]["CommuteQ", QECPauli["Z"]]},
    {True, False},
    TestID -> "QEC-Surface-Pauli-commutation"
]

(* The object and the verbs it replaces agree on every Pauli of two qubits. *)
VerificationTest[
    With[{all = Flatten[Outer[StringJoin, {"I", "X", "Y", "Z"}, {"I", "X", "Y", "Z"}]]},
        And @@ Flatten @ Table[
            QECPauli[a]["Product", b]["Vector"] === qecPauliProduct[a, b],
            {a, all}, {b, all}
        ]
    ],
    True,
    TestID -> "QEC-Surface-Pauli-product-matches-the-row-algebra"
]

VerificationTest[
    {Normal[QECPauli["Y"]["Matrix"]], Head[QECPauli["XZ"]["QuantumOperator"]]},
    {{{0, -I}, {I, 0}}, QuantumOperator},
    TestID -> "QEC-Surface-Pauli-as-a-matrix-and-an-operator"
]

VerificationTest[
    QECPauli[{"XI", "ZZ"}],
    {QECPauli["XI"], QECPauli["ZZ"]},
    TestID -> "QEC-Surface-Pauli-maps-over-a-list"
]

(* An object goes wherever a string or a row went. *)
VerificationTest[
    {five["Syndrome", QECPauli["XIIII"]] === five["Syndrome", "XIIII"], QECPauliQ[QECPauli["XX"]]},
    {True, True},
    TestID -> "QEC-Surface-Pauli-is-accepted-where-a-string-is"
]

VerificationTest[
    QECPauli["XQ"],
    $Failed,
    {QECPauli::invalid},
    TestID -> "QEC-Surface-Pauli-refuses-a-bad-string"
]

VerificationTest[
    QECPauli[{1, 2, 0}],
    $Failed,
    {QECPauli::badrow},
    TestID -> "QEC-Surface-Pauli-refuses-a-bad-row"
]


(* ============================================================================
   The renamed heads and the new properties
   ============================================================================ *)

VerificationTest[
    Head[QECFaultTolerantCircuit[{{"R", 1}, {"H", 1}, {"M", 1}}, steane]],
    QECFaultTolerantCircuit,
    TestID -> "QEC-Surface-fault-tolerant-circuit-head"
]

VerificationTest[
    StringQ[QECStim[QECCode["BitFlipCode"], QECNoiseModel["Circuit", 1/100], 1]],
    True,
    TestID -> "QEC-Surface-Stim-writer"
]

VerificationTest[
    With[{bf = QECCode["BitFlipCode"], noise = QECNoiseModel["Circuit", 1/100]},
        QECDetectorModel[bf, noise, 2]["StimString"] === QECStim[bf, noise, 2]
    ],
    True,
    TestID -> "QEC-Surface-detector-model-gives-its-Stim-source"
]

VerificationTest[
    {Head[QECCode["Catalog"]], MemberQ[Normal[Keys[QECCode["Catalog"]]], "SteaneCode"]},
    {Dataset, True},
    TestID -> "QEC-Surface-catalog-is-asked-of-the-head"
]

VerificationTest[
    {steane["RemoveQubit"] === QECRemoveQubit[steane],
     steane["RemoveQubit", 3] === QECRemoveQubit[steane, 3],
     QECCode["BitFlipCode"]["Concatenate", QECCode["PhaseFlipCode"]] ===
         QECConcatenate[QECCode["BitFlipCode"], QECCode["PhaseFlipCode"]],
     With[{c4 = QECCode[{"XXXX", "ZZZZ"}]}, c4["Paste", c4, 1, 1] === QECPasteCodes[c4, c4, 1, 1]]},
    {True, True, True, True},
    TestID -> "QEC-Surface-constructions-as-properties-agree-with-the-functions"
]

VerificationTest[
    SubsetQ[steane["Properties"], {"Concatenate", "RemoveQubit", "Paste"}],
    True,
    TestID -> "QEC-Surface-construction-properties-are-listed"
]

VerificationTest[
    QECCode["CSS", {{0, 0, 0, 1, 1, 1, 1}, {0, 1, 1, 0, 0, 1, 1}, {1, 0, 1, 0, 1, 0, 1}}]["Parameters"],
    {7, 1, 3},
    TestID -> "QEC-Surface-Hamming-matrix-is-an-ordinary-matrix"
]


(* ============================================================================
   The retired names: same value, one warning
   ============================================================================ *)

VerificationTest[
    QECPauliVector["XZ"],
    qecPauliVector["XZ"],
    {QECPauliVector::qecdeprecated},
    TestID -> "QEC-Surface-alias-PauliVector-works-and-warns"
]

(* Second use in the same session: the same value and no message. *)
VerificationTest[
    QECPauliVector["ZX"],
    qecPauliVector["ZX"],
    TestID -> "QEC-Surface-alias-warns-once"
]

VerificationTest[
    QECPauliString[{1, 1, 0}],
    "Y",
    {QECPauliString::qecdeprecated},
    TestID -> "QEC-Surface-alias-PauliString-works-and-warns"
]

VerificationTest[
    QECPauliWeight["XIZ"],
    2,
    {QECPauliWeight::qecdeprecated},
    TestID -> "QEC-Surface-alias-PauliWeight-works-and-warns"
]

VerificationTest[
    QECPauliCommuteQ["X", "Z"],
    False,
    {QECPauliCommuteQ::qecdeprecated},
    TestID -> "QEC-Surface-alias-PauliCommuteQ-works-and-warns"
]

VerificationTest[
    QECPauliProduct["X", "Z"],
    {1, 1, 3},
    {QECPauliProduct::qecdeprecated},
    TestID -> "QEC-Surface-alias-PauliProduct-works-and-warns"
]

VerificationTest[
    QECPauliPhase["-iX"],
    3,
    {QECPauliPhase::qecdeprecated},
    TestID -> "QEC-Surface-alias-PauliPhase-works-and-warns"
]

VerificationTest[
    Head[QECFaultTolerant[{{"R", 1}, {"M", 1}}, steane]],
    QECFaultTolerantCircuit,
    {QECFaultTolerant::qecdeprecated},
    TestID -> "QEC-Surface-alias-FaultTolerant-works-and-warns"
]

VerificationTest[
    QECStimCircuit[QECCode["BitFlipCode"]],
    QECStim[QECCode["BitFlipCode"]],
    {QECStimCircuit::qecdeprecated},
    TestID -> "QEC-Surface-alias-StimCircuit-works-and-warns"
]

VerificationTest[
    QECCodeCatalog[] === QECCode["Catalog"],
    True,
    {QECCodeCatalog::qecdeprecated},
    TestID -> "QEC-Surface-alias-CodeCatalog-works-and-warns"
]

VerificationTest[
    QECClassicalHammingMatrix[3],
    qecHammingMatrix[3],
    {QECClassicalHammingMatrix::qecdeprecated},
    TestID -> "QEC-Surface-alias-ClassicalHammingMatrix-works-and-warns"
]
