(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECCodeCatalog]
PackageExport[QECClassicalHammingMatrix]


(* ============================================================================ *)
(* Named codes and code families.                                               *)
(*                                                                              *)
(* Everything here is a constructor for the same QECCode object, reached through *)
(* the object itself rather than through a separate exported symbol per family:  *)
(* QECCode["SteaneCode"], QECCode["Repetition", 5], QECCode["CSS", hx, hz].      *)
(* The scalable families of roadmap item 4 (SurfaceCode, ToricCode, bivariate    *)
(* bicycle) slot in here as further named constructors without touching the      *)
(* rest of the package.                                                         *)
(* ============================================================================ *)


(* ---- fixed named codes ---- *)

(* PauliStabilizer's named codes give n stabilizers, which fix one codeword; the
   code is the group they generate minus the last one, which is the logical Z
   pinning that particular codeword. *)
$QECNamedCodes := $QECNamedCodes = <|
    "BitFlipCode" -> {"ZZI", "IZZ"},
    "PhaseFlipCode" -> {"XXI", "IXX"},
    "ShorCode" -> {
        "ZZIIIIIII", "IZZIIIIII", "IIIZZIIII", "IIIIZZIII", "IIIIIIZZI", "IIIIIIIZZ",
        "XXXXXXIII", "IIIXXXXXX"
    },
    "5QubitCode" -> Most[Wolfram`QuantumFramework`PauliStabilizer["5QubitCode"]["Stabilizers"]],
    "SteaneCode" -> Most[Wolfram`QuantumFramework`PauliStabilizer["SteaneCode"]["Stabilizers"]]
|>

QECCode::unknown = "`1` is not a known code or code family. Named codes: `2`. Families: `3`.";

$QECFamilies = {"Repetition", "PhaseRepetition", "DistanceTwo", "Hamming", "CSS"};

QECCode[name_String] := If[
    KeyExistsQ[$QECNamedCodes, name],
    QECCode[$QECNamedCodes[name]],
    Message[QECCode::unknown, name, StringRiffle[Keys[$QECNamedCodes], ", "], StringRiffle[$QECFamilies, ", "]];
    $Failed
]


(* ---- parametric families ---- *)

QECCode::repsize = "The repetition code needs at least 2 qubits.";
QECCode::evensize = "The [[n, n-2, 2]] family needs an even number of qubits, at least 4.";
QECCode::hammingsize = "The Hamming CSS construction needs r >= 3 for the classical code to be weakly self-dual.";

(* Neighbouring-pair checks: ZZ for bit flips, XX for phase flips (thesis sec. 2.2). *)
repetitionRows[n_Integer, letter_ : "Z"] := Table[
    QECPauliVector[StringJoin[ReplacePart[ConstantArray["I", n], {i -> letter, i + 1 -> letter}]]],
    {i, n - 1}
]

QECCode["Repetition", n_Integer] := If[
    n >= 2, QECCode[repetitionRows[n, "Z"]], Message[QECCode::repsize]; $Failed
]

QECCode["PhaseRepetition", n_Integer] := If[
    n >= 2, QECCode[repetitionRows[n, "X"]], Message[QECCode::repsize]; $Failed
]

(* The [[n, n-2, 2]] family of thesis sec. 8.1, saturating the quantum Singleton bound. *)
QECCode["DistanceTwo", n_Integer] := If[
    EvenQ[n] && n >= 4,
    QECCode[{StringRepeat["X", n], StringRepeat["Z", n]}],
    Message[QECCode::evensize]; $Failed
]

(* Columns are the nonzero binary words of length r: the classical [2^r-1, 2^r-1-r, 3]
   Hamming code. *)
QECClassicalHammingMatrix::usage = "QECClassicalHammingMatrix[r] gives the parity-check matrix of the classical Hamming code with r checks.";

QECClassicalHammingMatrix[r_Integer ? Positive] := Transpose[Table[IntegerDigits[j, 2, r], {j, 2^r - 1}]]

(* The CSS code built on it, [[2^r-1, 2^r-1-2r, 3]], which is the Steane code at r = 3. *)
QECCode["Hamming", r_Integer] := If[
    r >= 3, QECCode["CSS", QECClassicalHammingMatrix[r]], Message[QECCode::hammingsize]; $Failed
]


(* ---- the catalog ---- *)

QECCodeCatalog::usage = "QECCodeCatalog[] gives a Dataset of the available named codes and code families with their [[n, k, d]] parameters.";

$QECCatalogExamples = {
    {"Repetition", 5}, {"PhaseRepetition", 5}, {"DistanceTwo", 6}, {"Hamming", 4}
};

QECCodeCatalog[] := Dataset @ Association @ Join[
    Table[
        name -> AssociationThread[{"n", "k", "d"}, QECCode[name]["Parameters"]],
        {name, Keys[$QECNamedCodes]}
    ],
    Table[
        StringTemplate["QECCode[\"``\", ``]"][First[spec], Last[spec]] ->
            AssociationThread[{"n", "k", "d"}, (QECCode @@ spec)["Parameters"]],
        {spec, $QECCatalogExamples}
    ]
]
