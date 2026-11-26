(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Style.wlt

   The layer spoken the framework's way (November 2026).
   - Circuits the layer hands out are circuits QuantumCircuitOperator accepts as they
     are: "Reset" -> q, "Measurement" -> q, SuperDagger["S"] -> q, and a verification
     readout as a QuantumMeasurementOperator labelled "Herald".
   - Objects applied to what they act on: dec[syndrome], code[QECPauli[...]].
   - Equal compares what an object represents, not how it was written.
   - QuantumTensorProduct joins codes, as it joins every other framework object.
   - The framework's PauliStabilizer names for the textbook codes are accepted.
   - A detector model's Stim source is the property "Stim".
   - Every object has a predicate (QECCodeQ, ...), "ValidQ" and "AllProperties"; an
     unknown property gives ::undefprop and stays unevaluated; a refusal is a Failure;
     memoisation follows $QuantumFrameworkPropCache.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

steane = QECCode["SteaneCode"];
surface3 = QECCode["Surface"[3]];


(* ============================================================================
   Circuits
   ============================================================================ *)

VerificationTest[
    QECCode["BitFlipCode"]["SyndromeCircuit"]["Instructions"][[{1, 2, 4}]],
    {"Reset" -> 4, "CNOT" -> {1, 4}, "Measurement" -> 4},
    TestID -> "QEC-Style-reset-and-measurement-names"
]

VerificationTest[
    QECTransversalGate[steane, SuperDagger["S"]]["Instructions"][[1]],
    SuperDagger["S"] -> 1,
    TestID -> "QEC-Style-dagger-in-and-out"
]

(* Every list the layer hands out is a QuantumCircuitOperator as it stands. *)
VerificationTest[
    Head /@ {
        QuantumCircuitOperator[QECCode["SteaneCode"]["SyndromeCircuit"]["Instructions"]],
        QuantumCircuitOperator[QECCatState[4]["Instructions"]],
        QuantumCircuitOperator[QECFaultTolerantCircuit[{"Reset" -> 1, "H" -> 1, "Measurement" -> 1}, steane]["Instructions"]]
    },
    {QuantumCircuitOperator, QuantumCircuitOperator, QuantumCircuitOperator},
    TestID -> "QEC-Style-instructions-are-framework-circuits"
]

VerificationTest[
    With[{h = Cases[QECCatState[4]["Instructions"], _QuantumMeasurementOperator]},
        {Length[h], First[h]["Label"], QECCatState[4]["GateCounts"]["Herald"]}
    ],
    {1, "Herald", 1},
    TestID -> "QEC-Style-heralds-are-labelled-measurements"
]

VerificationTest[
    QECFaultTolerantCircuit[{"R" -> 1, "H" -> 1}, steane][[1]],
    "InvalidArguments",
    {QECFaultTolerantCircuit::oldform},
    TestID -> "QEC-Style-short-names-refused"
]


(* ============================================================================
   Objects applied to what they act on
   ============================================================================ *)

VerificationTest[
    With[{dec = QECDecoder[surface3, QECNoiseModel["Depolarizing"[1/100]], "Method" -> "Matching"]},
        dec[{0, 0, 0, 0, 1, 0, 0, 0}] === dec["Decode", {0, 0, 0, 0, 1, 0, 0, 0}]
    ],
    True,
    TestID -> "QEC-Style-decoder-applied-to-a-syndrome"
]

VerificationTest[
    steane[QECPauli["XIIIIII"]] === steane["Syndrome", "XIIIIII"],
    True,
    TestID -> "QEC-Style-code-applied-to-a-Pauli"
]


(* ============================================================================
   Equality
   ============================================================================ *)

VerificationTest[
    {QECCode["CSS"[{{0, 0, 0, 1, 1, 1, 1}, {0, 1, 1, 0, 0, 1, 1}, {1, 0, 1, 0, 1, 0, 1}}]] == steane,
     steane == QECCode["7QubitCode"],
     QECCode["Surface"[3]] == QECCode["Toric"[3]],
     QECCode[{"ZZI", "IZZ"}] == QECCode[{"ZZI", "ZIZ"}],
     QECCode[{"ZZI", "IZZ"}] == QECCode[{"-ZZI", "IZZ"}]},
    {True, True, False, True, False},
    TestID -> "QEC-Style-codes-equal-by-their-group"
]

VerificationTest[
    {QECNoiseModel["Depolarizing"[p]] == QECNoiseModel[<|"X" -> p/3, "Y" -> p/3, "Z" -> p/3|>],
     QECNoiseModel["BitFlip"[p]] == QECNoiseModel["PhaseFlip"[p]],
     QECPauli["XZ"] == QECPauli[{1, 0, 0, 1, 0}]},
    {True, False, True},
    TestID -> "QEC-Style-noise-and-paulis-equal-by-content"
]


(* ============================================================================
   Tensor product, names, Stim
   ============================================================================ *)

VerificationTest[
    {QuantumTensorProduct[steane, QECCode["5QubitCode"]]["Parameters"],
     QuantumTensorProduct[QECCode["Repetition"[3]], QECCode["Repetition"[3]], QECCode["Repetition"[3]]]["Qubits"]},
    {{12, 2, 3}, 9},
    TestID -> "QEC-Style-tensor-product-of-codes"
]

VerificationTest[
    {QECCode["9QubitCode"] == QECCode["ShorCode"], MemberQ[$QECCodeNames, "7QubitCode"]},
    {True, True},
    TestID -> "QEC-Style-pauli-stabilizer-names"
]

VerificationTest[
    With[{dem = QECDetectorModel[QECCode["BitFlipCode"], QECNoiseModel["Circuit"[1/100]], 1]},
        StringQ[dem["Stim"]]
    ],
    True,
    TestID -> "QEC-Style-stim-is-a-property"
]


(* ============================================================================
   Objects like the framework's: predicates, ValidQ, undefprop, Failure, cache
   ============================================================================ *)

VerificationTest[
    {QECCodeQ[steane], QECCodeQ[QECCode[<|"Qubits" -> 3|>]], QECCodeQ[3],
     QECNoiseModelQ[QECNoiseModel["Depolarizing"[p]]], QECDecoderQ[QECDecoder[steane]],
     QECCatStateQ[QECCatState[4]], steane["ValidQ"]},
    {True, False, False, True, True, True, True},
    TestID -> "QEC-Style-predicates"
]

VerificationTest[
    Complement[steane["AllProperties"], steane["Properties"]],
    {"AllProperties", "Properties", "ValidQ"},
    TestID -> "QEC-Style-all-properties"
]

VerificationTest[
    MatchQ[QECNoiseModel["Depolarizing"[p]]["Nonsense"], _QECNoiseModel["Nonsense"]],
    True,
    {QECNoiseModel::undefprop},
    TestID -> "QEC-Style-unknown-property-left-unevaluated"
]

VerificationTest[
    {#[[1]], FailureQ[#]} & @ QECCode[{"ZI", "IZ", "ZZ"}],
    {"overcomplete", True},
    {QECCode::overcomplete},
    TestID -> "QEC-Style-refusal-is-a-failure"
]

VerificationTest[
    Module[{f = Symbol["Wolfram`QuantumFramework`QEC`PackageScope`codeMinimumLogical"], n},
        QECClearCache[];
        n = Length[DownValues[Evaluate[f]]];
        Block[{Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = False}, QECCode["SteaneCode"]["Distance"]];
        {Length[DownValues[Evaluate[f]]] - n, QECCode["SteaneCode"]["Distance"]; Length[DownValues[Evaluate[f]]] - n}
    ],
    {0, 1},
    TestID -> "QEC-Style-cache-follows-the-framework"
]

VerificationTest[
    {#[[1]], FailureQ[#]} & @ QECDecoder[steane, QECNoiseModel["Depolarizing"[1/100]], "Method" -> "Matching"],
    {"notgraphlike", True},
    {QECDecoder::notgraphlike},
    TestID -> "QEC-Style-matching-refusal-is-a-failure"
]

VerificationTest[
    {#[[1]], FailureQ[#]} & @ steane["LogicalChannel", "noise"],
    {"noisespec", True},
    {QECCode::noisespec},
    TestID -> "QEC-Style-a-bad-noise-is-a-failure"
]
