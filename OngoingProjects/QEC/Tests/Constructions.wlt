(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Constructions.wlt

   Codes built out of codes: CSS, concatenation, qubit removal, pasting
   (Gottesman thesis sec. 3.5 and 8.6, QECC book ch. 5).

   The regression anchors live here.  Two of them are exact reconstructions of
   textbook codes by an entirely different route -- the Steane code as a CSS code
   over the classical Hamming matrix, and the Shor code as the phase-flip code
   concatenated with the bit-flip code -- and a third reproduces Gottesman's own
   worked surgery on the five-qubit code.  If a refactor breaks the algebra, these
   are what notice.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

(* Transitional: the rebuilt QEC core still lives under OngoingProjects/QEC/.
   Once it moves into Kernel/QEC/ and PacletInfo.wl lists its context, this Get
   disappears and the Needs above is enough.  The repo root is found from the
   loaded paclet, which RunTests.wls points at this checkout through
   PacletDirectoryLoad, so the tests always run against the source tree. *)
Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

(* Derived properties are memoised on the code's data.  A kernel that already
   held results from an earlier load would let those outrank the definitions just
   read, so the tests would silently check the old package. *)
QECClearCache[];


(* ============================================================================
   CSS
   ============================================================================ *)

VerificationTest[
    QECCode["CSS", QECClassicalHammingMatrix[3]]["Parameters"],
    {7, 1, 3},
    TestID -> "QEC-CSS-hamming-parameters"
]

(* The Steane code exactly: the generators come out in a different order, so the
   claim is about the stabilizer group, checked with signs in both directions. *)
VerificationTest[
    Module[{css = QECCode["CSS", QECClassicalHammingMatrix[3]], steane = QECCode["SteaneCode"]},
        AllTrue[steane["GeneratorVectors"], css["StabilizerMemberQ", #] &] &&
        AllTrue[css["GeneratorVectors"], steane["StabilizerMemberQ", #] &]
    ],
    True,
    TestID -> "QEC-CSS-reproduces-steane"
]

VerificationTest[
    Sort[QECCode["CSS", QECClassicalHammingMatrix[3]]["Generators"]] === Sort[QECCode["SteaneCode"]["Generators"]],
    True,
    TestID -> "QEC-CSS-reproduces-steane-generators"
]

VerificationTest[QECCode["CSS", {{1, 1, 1, 1}}, {{1, 1, 1, 1}}]["Parameters"], {4, 2, 2}, TestID -> "QEC-CSS-four-qubit"]

VerificationTest[QECCode["CSS", {{1, 1, 1, 1}}]["Distance"], 2, TestID -> "QEC-CSS-single-matrix"]

VerificationTest[QECCode["CSS", QECClassicalHammingMatrix[3]]["CSSQ"], True, TestID -> "QEC-CSS-is-css"]

VerificationTest[Quiet[QECCode["CSS", {{1, 1, 0}}, {{1, 0, 0}}]], $Failed, TestID -> "QEC-CSS-reject-non-orthogonal"]

VerificationTest[Quiet[QECCode["CSS", {{1, 1}}, {{1, 1, 1}}]], $Failed, TestID -> "QEC-CSS-reject-dimension-mismatch"]

VerificationTest[
    Quiet[Check[QECCode["CSS", {{1, 1, 0}}, {{1, 0, 0}}], "fired", QECCode::cssorth]],
    "fired",
    TestID -> "QEC-CSS-orthogonality-message"
]


(* ============================================================================
   Concatenation
   ============================================================================ *)

VerificationTest[
    QECConcatenate[QECCode["PhaseFlipCode"], QECCode["BitFlipCode"]]["Parameters"],
    {9, 1, 3},
    TestID -> "QEC-Concatenate-shor-parameters"
]

(* The Shor code exactly, again as a statement about the group. *)
VerificationTest[
    Module[{built = QECConcatenate[QECCode["PhaseFlipCode"], QECCode["BitFlipCode"]], shor = QECCode["ShorCode"]},
        AllTrue[shor["GeneratorVectors"], built["StabilizerMemberQ", #] &] &&
        AllTrue[built["GeneratorVectors"], shor["StabilizerMemberQ", #] &]
    ],
    True,
    TestID -> "QEC-Concatenate-reproduces-shor"
]

VerificationTest[
    Sort[QECConcatenate[QECCode["PhaseFlipCode"], QECCode["BitFlipCode"]]["Generators"]] === Sort[QECCode["ShorCode"]["Generators"]],
    True,
    TestID -> "QEC-Concatenate-reproduces-shor-generators"
]

(* [[n1,k,d1]] outer with [[n2,1,d2]] inner gives [[n1 n2, k, d1 d2]]. *)
VerificationTest[
    QECConcatenate[QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}], QECCode["BitFlipCode"]]["Parameters"],
    {12, 2, 2},
    TestID -> "QEC-Concatenate-multi-logical"
]

VerificationTest[
    QECConcatenate[QECCode["BitFlipCode"], QECCode["BitFlipCode"]]["Parameters"],
    {9, 1, 1},
    TestID -> "QEC-Concatenate-repetition-squared"
]

VerificationTest[
    Quiet[QECConcatenate[QECCode["PhaseFlipCode"], QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}]]],
    $Failed,
    TestID -> "QEC-Concatenate-reject-inner-with-many-logicals"
]

(* A Y in an outer generator translates to i Xbar Zbar, which must come back
   Hermitian: the resulting code must be constructible at all, and its generator
   phases must be real. *)
VerificationTest[
    Module[{outer = QECCode[{"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"}], built},
        built = QECConcatenate[outer, QECCode["BitFlipCode"]];
        built =!= $Failed && SubsetQ[{0, 2}, DeleteDuplicates[built["Phases"]]]
    ],
    True,
    TestID -> "QEC-Concatenate-Y-stays-hermitian"
]

(* d1 d2 is a lower bound on the concatenated distance, not an equality: the
   bit-flip code has quantum distance 1 but corrects bit flips at weight 3, and
   concatenating it under the [[5,1,3]] code buys distance 5, not 3.  The witness
   below realises it, and weights 1 to 4 are ruled out exhaustively by the search. *)
VerificationTest[
    QECConcatenate[QECCode[{"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"}], QECCode["BitFlipCode"]]["Parameters"],
    {15, 1, 5},
    TestID -> "QEC-Concatenate-five-with-bitflip"
]

VerificationTest[
    Module[{code = QECConcatenate[QECCode["5QubitCode"], QECCode["BitFlipCode"]]},
        QECPauliWeight[code["MinimumWeightLogical"]] === code["Distance"] &&
        code["LogicalPauliQ", code["MinimumWeightLogical"]]
    ],
    True,
    TestID -> "QEC-Concatenate-distance-witness"
]

(* Where both factors have a genuine quantum distance, the product is attained. *)
VerificationTest[
    QECConcatenate[QECCode["PhaseFlipCode"], QECCode["BitFlipCode"]]["Distance"],
    3,
    TestID -> "QEC-Concatenate-shor-distance"
]


(* ============================================================================
   Qubit removal
   ============================================================================ *)

VerificationTest[QECRemoveQubit[QECCode["5QubitCode"]]["Parameters"], {4, 2, 2}, TestID -> "QEC-Remove-five-parameters"]

(* Gottesman's own worked example: removing a qubit from the [[5,1,3]] code
   leaves the [[4,2,2]] code with these two generators. *)
VerificationTest[
    QECRemoveQubit[QECCode["5QubitCode"]]["Generators"],
    {"XZZX", "YXXY"},
    TestID -> "QEC-Remove-five-generators"
]

VerificationTest[QECRemoveQubit[QECCode["SteaneCode"]]["Parameters"], {6, 2, 2}, TestID -> "QEC-Remove-steane"]

VerificationTest[Head[QECRemoveQubit[QECCode["5QubitCode"], 3]], QECCode, TestID -> "QEC-Remove-specific-qubit"]

(* [[n,k,d]] becomes [[n-1,k+1,d-1]] whichever qubit is taken. *)
VerificationTest[
    Module[{code = QECCode["5QubitCode"]},
        DeleteDuplicates[
            Table[Quiet[QECRemoveQubit[code, q]], {q, 5}] /. c_QECCode :> c["Parameters"]
        ]
    ],
    {{4, 2, 2}},
    TestID -> "QEC-Remove-any-qubit-same-parameters"
]

VerificationTest[Quiet[QECRemoveQubit[QECCode["BitFlipCode"]]], $Failed, TestID -> "QEC-Remove-reject-no-pivot"]

VerificationTest[
    Quiet[Check[QECRemoveQubit[QECCode["BitFlipCode"]], "fired", QECRemoveQubit::nopivot]],
    "fired",
    TestID -> "QEC-Remove-nopivot-message"
]


(* ============================================================================
   Pasting
   ============================================================================ *)

VerificationTest[
    QECPasteCodes[
        QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}],
        QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}], 1, 1
    ]["Parameters"],
    {8, 5, 2},
    TestID -> "QEC-Paste-parameters"
]

VerificationTest[
    QECPasteCodes[
        QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}],
        QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}], 1, 1
    ]["Generators"],
    {"XXXXIIII", "IIIIXXXX", "ZZZZZZZZ"},
    TestID -> "QEC-Paste-generators"
]

VerificationTest[
    Quiet[QECPasteCodes[QECCode["BitFlipCode"], QECCode["SteaneCode"], 1, 1]],
    $Failed,
    TestID -> "QEC-Paste-reject-unequal-pairing"
]

VerificationTest[
    Quiet[QECPasteCodes[QECCode["BitFlipCode"], QECCode["PhaseFlipCode"], 5, 1]],
    $Failed,
    TestID -> "QEC-Paste-reject-bad-solo-count"
]

VerificationTest[
    Quiet[Check[QECPasteCodes[QECCode["BitFlipCode"], QECCode["SteaneCode"], 1, 1], "fired", QECPasteCodes::pairing]],
    "fired",
    TestID -> "QEC-Paste-pairing-message"
]

(* Pasting is a tensor construction: the two blocks keep their own qubits. *)
VerificationTest[
    QECPasteCodes[QECCode["BitFlipCode"], QECCode["PhaseFlipCode"], 2, 2]["Qubits"],
    6,
    TestID -> "QEC-Paste-qubit-count"
]
