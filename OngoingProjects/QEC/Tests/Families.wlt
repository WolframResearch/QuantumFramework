(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Families.wlt

   Named codes and parametric families, all reached through the QECCode object
   itself: QECCode["SteaneCode"], QECCode["Repetition", 5], QECCode["Hamming", 4].
   The scalable families of roadmap item 4 join this list as further names.
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
   Named codes
   ============================================================================ *)

VerificationTest[QECCode["5QubitCode"]["Parameters"], {5, 1, 3}, TestID -> "QEC-Families-five-parameters"]

VerificationTest[
    QECCode["5QubitCode"]["Generators"],
    {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"},
    TestID -> "QEC-Families-five-generators"
]

VerificationTest[QECCode["SteaneCode"]["Parameters"], {7, 1, 3}, TestID -> "QEC-Families-steane-parameters"]

VerificationTest[QECCode["ShorCode"]["Parameters"], {9, 1, 3}, TestID -> "QEC-Families-shor-parameters"]

VerificationTest[QECCode["BitFlipCode"]["Generators"], {"ZZI", "IZZ"}, TestID -> "QEC-Families-bitflip-generators"]

VerificationTest[QECCode["PhaseFlipCode"]["Generators"], {"XXI", "IXX"}, TestID -> "QEC-Families-phaseflip-generators"]

VerificationTest[Quiet[QECCode["NoSuchCode"]], $Failed, TestID -> "QEC-Families-unknown-name"]

VerificationTest[
    Quiet[Check[QECCode["NoSuchCode"], "fired", QECCode::unknown]],
    "fired",
    TestID -> "QEC-Families-unknown-name-message"
]


(* ============================================================================
   Repetition codes
   ============================================================================ *)

VerificationTest[
    QECCode["Repetition", 3]["Generators"],
    QECCode["BitFlipCode"]["Generators"],
    TestID -> "QEC-Families-repetition-three-is-bitflip"
]

VerificationTest[
    QECCode["PhaseRepetition", 3]["Generators"],
    QECCode["PhaseFlipCode"]["Generators"],
    TestID -> "QEC-Families-phase-repetition-three-is-phaseflip"
]

VerificationTest[QECCode["Repetition", 5]["Parameters"], {5, 1, 1}, TestID -> "QEC-Families-repetition-five"]

VerificationTest[Quiet[QECCode["Repetition", 1]], $Failed, TestID -> "QEC-Families-repetition-too-small"]

(* The repetition code detects bit flips and is blind to phase flips, which is why
   its quantum distance is 1 whatever n is. *)
VerificationTest[
    Table[QECCode["Repetition", n]["Distance"], {n, 2, 6}],
    {1, 1, 1, 1, 1},
    TestID -> "QEC-Families-repetition-distance-is-one"
]

VerificationTest[
    Table[QECCode["Repetition", n]["Syndrome", StringJoin[ReplacePart[ConstantArray["I", n], 1 -> "X"]]], {n, 3, 5}],
    {{1, 0}, {1, 0, 0}, {1, 0, 0, 0}},
    TestID -> "QEC-Families-repetition-detects-bit-flips"
]


(* ============================================================================
   The [[n, n-2, 2]] family
   ============================================================================ *)

VerificationTest[QECCode["DistanceTwo", 4]["Parameters"], {4, 2, 2}, TestID -> "QEC-Families-distance-two-4"]

VerificationTest[QECCode["DistanceTwo", 6]["Parameters"], {6, 4, 2}, TestID -> "QEC-Families-distance-two-6"]

VerificationTest[Quiet[QECCode["DistanceTwo", 5]], $Failed, TestID -> "QEC-Families-distance-two-odd"]

VerificationTest[Quiet[QECCode["DistanceTwo", 2]], $Failed, TestID -> "QEC-Families-distance-two-too-small"]

(* It saturates the quantum Singleton bound, n - k >= 2(d - 1). *)
VerificationTest[
    AllTrue[{4, 6, 8},
        With[{p = QECCode["DistanceTwo", #]["Parameters"]},
            p[[1]] - p[[2]] === 2 (p[[3]] - 1)
        ] &
    ],
    True,
    TestID -> "QEC-Families-distance-two-saturates-singleton"
]


(* ============================================================================
   The Hamming CSS family
   ============================================================================ *)

VerificationTest[Dimensions[QECClassicalHammingMatrix[3]], {3, 7}, TestID -> "QEC-Families-hamming-matrix-shape"]

VerificationTest[
    QECClassicalHammingMatrix[3],
    {{0, 0, 0, 1, 1, 1, 1}, {0, 1, 1, 0, 0, 1, 1}, {1, 0, 1, 0, 1, 0, 1}},
    TestID -> "QEC-Families-hamming-matrix"
]

VerificationTest[
    Sort[QECCode["Hamming", 3]["Generators"]] === Sort[QECCode["SteaneCode"]["Generators"]],
    True,
    TestID -> "QEC-Families-hamming-three-is-steane"
]

VerificationTest[QECCode["Hamming", 4]["Parameters"], {15, 7, 3}, TestID -> "QEC-Families-hamming-four"]

VerificationTest[Quiet[QECCode["Hamming", 2]], $Failed, TestID -> "QEC-Families-hamming-too-small"]

(* [[2^r - 1, 2^r - 1 - 2r, 3]]. *)
VerificationTest[
    Table[Take[QECCode["Hamming", r]["Parameters"], 2], {r, 3, 5}],
    Table[{2^r - 1, 2^r - 1 - 2 r}, {r, 3, 5}],
    TestID -> "QEC-Families-hamming-parameters-formula"
]


(* ============================================================================
   The catalog
   ============================================================================ *)

VerificationTest[Head[QECCodeCatalog[]], Dataset, TestID -> "QEC-Families-catalog-head"]

VerificationTest[
    KeyExistsQ[Normal[QECCodeCatalog[]], "SteaneCode"],
    True,
    TestID -> "QEC-Families-catalog-contains-steane"
]

VerificationTest[
    AllTrue[Values[Normal[QECCodeCatalog[]]], Keys[#] === {"n", "k", "d"} &],
    True,
    TestID -> "QEC-Families-catalog-columns"
]

VerificationTest[
    Normal[QECCodeCatalog[]]["ShorCode"],
    <|"n" -> 9, "k" -> 1, "d" -> 3|>,
    TestID -> "QEC-Families-catalog-shor-row"
]
