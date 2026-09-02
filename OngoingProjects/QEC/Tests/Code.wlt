(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Code.wlt

   Building the code object and reading it back: validation of the generators,
   the direct properties, and the two-way traffic with PauliStabilizer.
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

qecFive = {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"};


(* ============================================================================
   Construction
   ============================================================================ *)

VerificationTest[Head[QECCode[{"ZZI", "IZZ"}]], QECCode, TestID -> "QEC-Code-head"]

VerificationTest[QECCode[{"ZZI", "IZZ"}]["Parameters"], {3, 1, 1}, TestID -> "QEC-Code-parameters-bitflip"]

VerificationTest[QECCode[qecFive]["Parameters"], {5, 1, 3}, TestID -> "QEC-Code-parameters-five"]

VerificationTest[
    QECCode[{"IIIXXXX", "XIXIXIX", "IXXIIXX", "IIIZZZZ", "ZIZIZIZ", "IZZIIZZ"}]["Parameters"],
    {7, 1, 3},
    TestID -> "QEC-Code-parameters-steane"
]

VerificationTest[
    QECCode[{StringRepeat["X", 8], StringRepeat["Z", 8]}]["Parameters"],
    {8, 6, 2},
    TestID -> "QEC-Code-parameters-distance-two"
]

VerificationTest[QECCode[{"ZZI", "IZZ"}]["Qubits"], 3, TestID -> "QEC-Code-qubits"]

VerificationTest[QECCode[{"ZZI", "IZZ"}]["StabilizerCount"], 2, TestID -> "QEC-Code-stabilizer-count"]

VerificationTest[QECCode[{"ZZI", "IZZ"}]["LogicalQubits"], 1, TestID -> "QEC-Code-logical-qubits"]

VerificationTest[QECCode[{"ZZI", "IZZ"}]["Generators"], {"ZZI", "IZZ"}, TestID -> "QEC-Code-generators"]

VerificationTest[
    QECCode[{"ZZI", "IZZ"}]["CheckMatrix"],
    {{0, 0, 0, 1, 1, 0}, {0, 0, 0, 0, 1, 1}},
    TestID -> "QEC-Code-check-matrix"
]

VerificationTest[QECCode[{"-ZI", "IZ"}]["Signs"], {-1, 1}, TestID -> "QEC-Code-signs"]

VerificationTest[QECCode[{"-ZI", "IZ"}]["Phases"], {2, 0}, TestID -> "QEC-Code-phases"]

VerificationTest[QECCode[{"-ZI", "IZ"}]["Generators"], {"-ZI", "IZ"}, TestID -> "QEC-Code-signed-generators"]

(* Rows are accepted as readily as strings, and give the same code. *)
VerificationTest[
    QECCode[QECPauliVector /@ {"ZZI", "IZZ"}] === QECCode[{"ZZI", "IZZ"}],
    True,
    TestID -> "QEC-Code-from-rows"
]


(* ============================================================================
   Validation
   ============================================================================ *)

VerificationTest[Quiet[QECCode[{"XXX", "ZZZ"}]], $Failed, TestID -> "QEC-Code-reject-noncommuting"]

VerificationTest[Quiet[QECCode[{"ZZI", "IZZ", "ZIZ"}]], $Failed, TestID -> "QEC-Code-reject-dependent"]

VerificationTest[Quiet[QECCode[{"ZQI", "IZZ"}]], $Failed, TestID -> "QEC-Code-reject-bad-letter"]

VerificationTest[Quiet[QECCode[{"ZZI", "IZZZ"}]], $Failed, TestID -> "QEC-Code-reject-ragged"]

VerificationTest[Quiet[QECCode[{"ZI", "IZ", "ZZ"}]], $Failed, TestID -> "QEC-Code-reject-overcomplete"]

(* Each rejection issues its own message, and the non-commuting one names the
   offending pair.  Check fires only for the message it is given, so these assert
   which message came out, not merely that something failed. *)
VerificationTest[
    Quiet[Check[QECCode[{"XXX", "ZZZ"}], "fired", QECCode::noncomm]],
    "fired",
    TestID -> "QEC-Code-noncommuting-message"
]

VerificationTest[
    Quiet[Check[QECCode[{"ZZI", "IZZ", "ZIZ"}], "fired", QECCode::dep]],
    "fired",
    TestID -> "QEC-Code-dependent-message"
]

VerificationTest[
    Quiet[Check[QECCode[{"ZI", "IZ", "ZZ"}], "fired", QECCode::overcomplete]],
    "fired",
    TestID -> "QEC-Code-overcomplete-message"
]


(* ============================================================================
   Properties
   ============================================================================ *)

VerificationTest[
    AllTrue[QECCode[{"ZZI", "IZZ"}]["Properties"], StringQ],
    True,
    TestID -> "QEC-Code-properties-list"
]

VerificationTest[
    MemberQ[QECCode[{"ZZI", "IZZ"}]["Properties"], "Distance"],
    True,
    TestID -> "QEC-Code-properties-contains-distance"
]

VerificationTest[
    Quiet[QECCode[{"ZZI", "IZZ"}]["NoSuchKey"]],
    Missing["NotFound", "NoSuchKey"],
    TestID -> "QEC-Code-unknown-property"
]

VerificationTest[
    Quiet[Check[QECCode[{"ZZI", "IZZ"}]["NoSuchKey"], "fired", QECCode::noprop]],
    "fired",
    TestID -> "QEC-Code-unknown-property-message"
]

VerificationTest[
    QECCode[{"ZZI", "IZZ"}]["GeneratorVectors"],
    QECPauliVector /@ {"ZZI", "IZZ"},
    TestID -> "QEC-Code-generator-vectors"
]


(* ============================================================================
   Traffic with the engine
   ============================================================================ *)

VerificationTest[
    Head[QECCode[{"ZZI", "IZZ"}]["PauliStabilizer"]],
    PauliStabilizer,
    TestID -> "QEC-Code-to-pauli-stabilizer"
]

VerificationTest[
    Normal[QECCode[{"ZI", "IZ"}]["State"]["StateVector"]],
    {1, 0, 0, 0},
    TestID -> "QEC-Code-state-vector"
]

VerificationTest[
    QECCode[{"ZI", "IZ"}]["PauliStabilizer"]["Expectation", "ZZ"],
    1,
    TestID -> "QEC-Code-expectation"
]

VerificationTest[
    QECCode[PauliStabilizer[2]]["Parameters"],
    {2, 0, Infinity},
    TestID -> "QEC-Code-from-pauli-stabilizer"
]

VerificationTest[
    QECCode[PauliStabilizer[2]["X", 1]]["Signs"],
    {-1, 1},
    TestID -> "QEC-Code-from-pauli-stabilizer-signs"
]

(* A code built from a state's stabilizers has no logical qubits, so its distance
   is Infinity: there is nothing left to protect. *)
VerificationTest[
    QECCode[PauliStabilizer[3]["H", 1]["CNOT", 1, 2]]["LogicalQubits"],
    0,
    TestID -> "QEC-Code-state-has-no-logical-qubits"
]

VerificationTest[
    Module[{code = QECCode[{"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"}]},
        AllTrue[code["Generators"], code["PauliStabilizer"]["Expectation", #] === 1 &]
    ],
    True,
    TestID -> "QEC-Code-stabilizer-roundtrip"
]
