(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Syndrome.wlt

   Syndromes, the lookup decoder, and the correction cycle.

   The syndrome is cross-checked against the engine: measuring a check on the
   damaged state must give the same bit as the symplectic product does.  That is
   the oracle the package cannot fake, because it goes through an actual
   stabilizer state rather than through our own linear algebra.
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

qecWeightOne[n_] := QECPauliString /@ Wolfram`QuantumFramework`QEC`PackageScope`weightOneVectors[n];


(* ============================================================================
   Syndromes
   ============================================================================ *)

VerificationTest[QECCode["BitFlipCode"]["Syndrome", "XII"], {1, 0}, TestID -> "QEC-Syndrome-bitflip-X1"]

VerificationTest[QECCode["BitFlipCode"]["Syndrome", "IXI"], {1, 1}, TestID -> "QEC-Syndrome-bitflip-X2"]

VerificationTest[QECCode["BitFlipCode"]["Syndrome", "IIX"], {0, 1}, TestID -> "QEC-Syndrome-bitflip-X3"]

VerificationTest[QECCode["BitFlipCode"]["Syndrome", "III"], {0, 0}, TestID -> "QEC-Syndrome-identity"]

VerificationTest[QECCode["BitFlipCode"]["Syndrome", "ZII"], {0, 0}, TestID -> "QEC-Syndrome-undetectable"]

VerificationTest[QECCode[qecFive]["Syndrome", "IXIII"], {1, 0, 0, 0}, TestID -> "QEC-Syndrome-five-qubit"]

(* A sign on the error cannot change which checks it anticommutes with. *)
VerificationTest[
    QECCode["BitFlipCode"]["Syndrome", "-XII"],
    QECCode["BitFlipCode"]["Syndrome", "XII"],
    TestID -> "QEC-Syndrome-sign-irrelevant"
]

VerificationTest[Quiet[QECCode["BitFlipCode"]["Syndrome", "XIII"]], $Failed, TestID -> "QEC-Syndrome-wrong-size"]

VerificationTest[Quiet[QECCode["BitFlipCode"]["Syndrome", "XQI"]], $Failed, TestID -> "QEC-Syndrome-bad-letter"]

(* The five-qubit code is perfect: its fifteen weight-one errors have fifteen
   distinct nonzero syndromes, which with the trivial one fill all sixteen. *)
VerificationTest[
    Module[{code = QECCode[qecFive], syndromes},
        syndromes = code["Syndrome", #] & /@ qecWeightOne[5];
        Length[DeleteDuplicates[syndromes]] === 15 && FreeQ[syndromes, {0, 0, 0, 0}]
    ],
    True,
    TestID -> "QEC-Syndrome-five-qubit-distinct"
]

(* Cross-check against the engine: the syndrome bit is (1 - <check>)/2 measured on
   the damaged state. *)
VerificationTest[
    Module[{code = QECCode[qecFive], errors, ours, oracle},
        errors = Flatten[Table[{p, q}, {q, 1, 5}, {p, {"X", "Y", "Z"}}], 1];
        ours = code["Syndrome", StringJoin[ReplacePart[ConstantArray["I", 5], #[[2]] -> #[[1]]]]] & /@ errors;
        oracle = Function[e,
            (1 - (PauliStabilizer["5QubitCode"][e[[1]], e[[2]]]["Expectation", #] & /@ qecFive)) / 2
        ] /@ errors;
        ours === oracle
    ],
    True,
    TestID -> "QEC-Syndrome-vs-engine"
]

VerificationTest[
    Module[{code = QECCode["SteaneCode"], errors, ours, oracle},
        errors = Flatten[Table[{p, q}, {q, 1, 7}, {p, {"X", "Y", "Z"}}], 1];
        ours = code["Syndrome", StringJoin[ReplacePart[ConstantArray["I", 7], #[[2]] -> #[[1]]]]] & /@ errors;
        oracle = Function[e,
            (1 - (PauliStabilizer["SteaneCode"][e[[1]], e[[2]]]["Expectation", #] & /@ code["Generators"])) / 2
        ] /@ errors;
        ours === oracle
    ],
    True,
    TestID -> "QEC-Syndrome-vs-engine-steane"
]

VerificationTest[Length[QECCode["BitFlipCode"]["SyndromeTable"]], 9, TestID -> "QEC-Syndrome-table-size"]

VerificationTest[QECCode["BitFlipCode"]["SyndromeTable"]["XII"], {1, 0}, TestID -> "QEC-Syndrome-table-entry"]

VerificationTest[QECCode["BitFlipCode"]["SyndromeTable"]["ZII"], {0, 0}, TestID -> "QEC-Syndrome-table-trivial"]

VerificationTest[QECCode[qecFive]["SyndromeTable"]["IXIII"], {1, 0, 0, 0}, TestID -> "QEC-Syndrome-table-five"]


(* ============================================================================
   The decoder
   ============================================================================ *)

(* The bit-flip code has d = 1 and so guarantees nothing, but its weight-one table
   is exactly what it is for: correcting single X errors.  The decoder's reach is
   Max[t, 1] for that reason. *)
VerificationTest[
    QECPauliString /@ QECCode["BitFlipCode"]["Decoder"],
    <|{0, 0} -> "III", {1, 0} -> "XII", {1, 1} -> "IXI", {0, 1} -> "IIX"|>,
    TestID -> "QEC-Syndrome-decoder-bitflip"
]

VerificationTest[Length[QECCode[qecFive]["Decoder"]], 16, TestID -> "QEC-Syndrome-decoder-five-complete"]

VerificationTest[QECCode[qecFive]["Decode", {1, 0, 0, 0}], "IXIII", TestID -> "QEC-Syndrome-decode-five"]

VerificationTest[QECCode[qecFive]["Decode", {0, 0, 0, 0}], "IIIII", TestID -> "QEC-Syndrome-decode-trivial"]

VerificationTest[QECCode["BitFlipCode"]["Decode", {0, 0}], "III", TestID -> "QEC-Syndrome-decode-bitflip-trivial"]

(* Every weight-one error is recovered from its own syndrome. *)
VerificationTest[
    Module[{code = QECCode[qecFive]},
        AllTrue[qecWeightOne[5], code["Decode", code["Syndrome", #]] === # &]
    ],
    True,
    TestID -> "QEC-Syndrome-decode-inverts-weight-one"
]

VerificationTest[
    MissingQ[QECCode["SteaneCode"]["Decode", {1, 1, 0, 0, 0, 1}]],
    True,
    TestID -> "QEC-Syndrome-decode-unknown"
]

VerificationTest[Quiet[QECCode["BitFlipCode"]["Decode", {1, 0, 1}]], $Failed, TestID -> "QEC-Syndrome-decode-wrong-length"]

VerificationTest[Quiet[QECCode["BitFlipCode"]["Decode", {1, 2}]], $Failed, TestID -> "QEC-Syndrome-decode-not-binary"]

(* The table reaches the weight the code can actually correct.  On a distance-5
   code that is weight two, which the prototype's weight-one-only table could not
   reach: it declared correctable errors undecodable. *)
VerificationTest[
    Wolfram`QuantumFramework`QEC`PackageScope`codeCorrectableWeight[First[QECCode["SteaneCode"]]],
    1,
    TestID -> "QEC-Syndrome-correctable-weight-steane"
]

VerificationTest[
    Wolfram`QuantumFramework`QEC`PackageScope`codeCorrectableWeight[First[QECCode["BitFlipCode"]]],
    0,
    TestID -> "QEC-Syndrome-correctable-weight-bitflip"
]

VerificationTest[
    Module[{code = QECCode["Repetition", 5]},
        {code["Distance"], Wolfram`QuantumFramework`QEC`PackageScope`codeCorrectableWeight[First[code]]}
    ],
    {1, 0},
    TestID -> "QEC-Syndrome-correctable-weight-repetition"
]

(* Every stored correction has least weight for its syndrome. *)
VerificationTest[
    Module[{code = QECCode[qecFive], decoder},
        decoder = code["Decoder"];
        AllTrue[Normal[decoder],
            With[{syn = First[#], corr = Last[#]},
                QECPauliWeight[corr] === Min[QECPauliWeight /@ Select[
                    Prepend[Wolfram`QuantumFramework`QEC`PackageScope`weightOneVectors[5], QECPauliVector["IIIII"]],
                    code["Syndrome", #] === syn &
                ]]
            ] &
        ]
    ],
    True,
    TestID -> "QEC-Syndrome-decoder-minimum-weight"
]


VerificationTest[
    Wolfram`QuantumFramework`QEC`PackageScope`codeDecoderReach[First[QECCode["BitFlipCode"]]],
    1,
    TestID -> "QEC-Syndrome-decoder-reach-floor"
]

(* Asking for more reach than the default gives a bigger table, and the weight-one
   entries of the default table survive in it unchanged. *)
VerificationTest[
    Module[{code = QECCode["SteaneCode"], small, big},
        small = code["Decoder"];
        big = code["Decoder", 2];
        Length[big] > Length[small] && AllTrue[Keys[small], big[Key[#]] === small[Key[#]] &]
    ],
    True,
    TestID -> "QEC-Syndrome-decoder-deeper-reach"
]


(* ============================================================================
   The correction cycle
   ============================================================================ *)

VerificationTest[QECCode["BitFlipCode"]["CorrectionCycle", "IXI"]["Success"], True, TestID -> "QEC-Cycle-bitflip-success"]

VerificationTest[QECCode["BitFlipCode"]["CorrectionCycle", "IXI"]["Correction"], "IXI", TestID -> "QEC-Cycle-bitflip-correction"]

VerificationTest[QECCode["BitFlipCode"]["CorrectionCycle", "IXI"]["Residual"], "III", TestID -> "QEC-Cycle-bitflip-residual"]

VerificationTest[QECCode["BitFlipCode"]["CorrectionCycle", "IXI"]["Outcome"], "Corrected", TestID -> "QEC-Cycle-bitflip-outcome"]

VerificationTest[QECCode[qecFive]["CorrectionCycle", "IIIII"]["Success"], True, TestID -> "QEC-Cycle-no-error"]

VerificationTest[
    AllTrue[qecWeightOne[5], QECCode[qecFive]["CorrectionCycle", #]["Success"] &],
    True,
    TestID -> "QEC-Cycle-all-weight-one"
]

VerificationTest[
    AllTrue[qecWeightOne[7], QECCode["SteaneCode"]["CorrectionCycle", #]["Success"] &],
    True,
    TestID -> "QEC-Cycle-all-weight-one-steane"
]

(* An undetectable error on the bit-flip code is a logical error, not an
   undecodable syndrome.  The prototype reported both as Success -> False with
   nothing to tell them apart. *)
VerificationTest[QECCode["BitFlipCode"]["CorrectionCycle", "ZII"]["Success"], False, TestID -> "QEC-Cycle-logical-failure"]

VerificationTest[
    QECCode["BitFlipCode"]["CorrectionCycle", "ZII"]["Outcome"],
    "LogicalError",
    TestID -> "QEC-Cycle-logical-outcome"
]

VerificationTest[
    QECCode["BitFlipCode"]["CorrectionCycle", "ZII"]["Syndrome"],
    {0, 0},
    TestID -> "QEC-Cycle-logical-has-trivial-syndrome"
]

VerificationTest[QECCode[qecFive]["CorrectionCycle", "XXIII"]["Success"], False, TestID -> "QEC-Cycle-weight-two-fails"]

VerificationTest[QECCode[qecFive]["CorrectionCycle", "XXIII"]["Correction"], "IIIZI", TestID -> "QEC-Cycle-weight-two-correction"]

VerificationTest[
    QECCode[qecFive]["CorrectionCycle", "XXIII"]["Outcome"],
    "LogicalError",
    TestID -> "QEC-Cycle-weight-two-outcome"
]

VerificationTest[
    QECCode["SteaneCode"]["CorrectionCycle", "IXIIZII"]["Outcome"],
    "Undecodable",
    TestID -> "QEC-Cycle-undecodable-outcome"
]

VerificationTest[
    MissingQ[QECCode["SteaneCode"]["CorrectionCycle", "IXIIZII"]["Correction"]],
    True,
    TestID -> "QEC-Cycle-undecodable-correction"
]

(* The residual is the actual product of error and correction, phase included. *)
VerificationTest[
    Module[{code = QECCode[qecFive]},
        AllTrue[qecWeightOne[5],
            code["CorrectionCycle", #]["Residual"] === QECPauliString[QECPauliProduct[#, code["CorrectionCycle", #]["Correction"]]] &
        ]
    ],
    True,
    TestID -> "QEC-Cycle-residual-is-the-product"
]

(* Success means the residual is a stabilizer up to the unobservable global phase. *)
VerificationTest[
    Module[{code = QECCode["SteaneCode"]},
        AllTrue[qecWeightOne[7],
            With[{result = code["CorrectionCycle", #]},
                result["Success"] === code["StabilizerMemberQ", StringDelete[result["Residual"], StartOfString ~~ ("-i" | "-" | "i")]]
            ] &
        ]
    ],
    True,
    TestID -> "QEC-Cycle-success-means-stabilizer"
]


(* ============================================================================
   The cycle on an actual state
   ============================================================================ *)

VerificationTest[
    QECCode["BitFlipCode"]["PhysicalCorrectionCycle", "IXI"]["Success"],
    True,
    TestID -> "QEC-Physical-bitflip-success"
]

VerificationTest[
    QECCode["BitFlipCode"]["PhysicalCorrectionCycle", "IXI"]["Correction"],
    "IXI",
    TestID -> "QEC-Physical-bitflip-correction"
]

(* A residual logical Z fixes the fiducial codeword, so the physical cycle counts
   it as restored where the symplectic cycle calls it a logical error. *)
VerificationTest[
    QECCode["BitFlipCode"]["PhysicalCorrectionCycle", "ZII"]["Success"],
    True,
    TestID -> "QEC-Physical-logical-Z-restores-fiducial-state"
]

VerificationTest[
    {QECCode["BitFlipCode"]["PhysicalCorrectionCycle", "ZII"]["Success"],
     QECCode["BitFlipCode"]["CorrectionCycle", "ZII"]["Success"]},
    {True, False},
    TestID -> "QEC-Physical-and-symplectic-disagree-as-designed"
]

VerificationTest[
    AllTrue[qecWeightOne[5], QECCode["5QubitCode"]["PhysicalCorrectionCycle", #]["Success"] &],
    True,
    TestID -> "QEC-Physical-five-qubit"
]

VerificationTest[
    AllTrue[qecWeightOne[3],
        QECCode["BitFlipCode"]["PhysicalCorrectionCycle", #]["Syndrome"] === QECCode["BitFlipCode"]["Syndrome", #] &
    ],
    True,
    TestID -> "QEC-Physical-syndromes-agree"
]


(* ============================================================================
   The quantum Hamming bound
   ============================================================================ *)

VerificationTest[QECCode[qecFive]["PerfectQ"], True, TestID -> "QEC-Perfect-five-qubit"]

VerificationTest[QECCode["BitFlipCode"]["PerfectQ"], False, TestID -> "QEC-Perfect-bitflip"]

VerificationTest[QECCode["SteaneCode"]["PerfectQ"], False, TestID -> "QEC-Perfect-steane"]

VerificationTest[
    QECCode[{StringRepeat["X", 8], StringRepeat["Z", 8]}]["PerfectQ"],
    False,
    TestID -> "QEC-Perfect-distance-two"
]
