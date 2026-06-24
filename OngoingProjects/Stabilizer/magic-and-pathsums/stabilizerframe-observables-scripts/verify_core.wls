#!/usr/bin/env wolframscript
(* Core verification: frameExpectation / frameNorm2 / frameInnerProduct vs dense *)
(* over random Clifford+T circuits. Ends with a result Association.              *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Stabilizer/stabilizerframe-observables-scripts/framealgebra.wls"];

SeedRandom[2026];

(* dense Pauli-string matrix, qubit 1 = most significant *)
pauliMat[str_String] := Module[{sign, body, mats},
    {sign, body} = If[StringStartsQ[str, "-"], {-1, StringDrop[str, 1]}, {1, str}];
    mats = Characters[body] /. {"I" -> PauliMatrix[0],
        "X" -> PauliMatrix[1], "Y" -> PauliMatrix[2], "Z" -> PauliMatrix[3]};
    sign If[Length[mats] == 1, First[mats], KroneckerProduct @@ mats]
];
denseExpect[f_, str_] := Module[{v = Normal @ f["StateVector"]},
    Conjugate[v] . (pauliMat[str] . v)
];

(* random circuit spec list (arrow forms) with >= 1 T, first T makes the frame. *)
randomSpecs[n_, nClifford_, nT_] := Module[{gateset, cliff, ts},
    gateset = If[n >= 2, {"H", "S", "V", "X", "Y", "Z", "CNOT", "CZ", "SWAP"}, {"H", "S", "V", "X", "Y", "Z"}];
    cliff = Table[Module[{g = RandomChoice[gateset]},
        If[MemberQ[{"CNOT", "CZ", "SWAP"}, g], g -> RandomSample[Range[n], 2], g -> RandomInteger[{1, n}]]],
        nClifford];
    ts = Table["T" -> RandomInteger[{1, n}], nT];
    Join[Take[cliff, Min[2, Length[cliff]]], ts, Drop[cliff, Min[2, Length[cliff]]]]
];

buildFrame[n_, nClifford_, nT_] := Fold[#1[#2] &, PauliStabilizer[n], randomSpecs[n, nClifford, nT]];
coherentFrameQ[f_] := Head[f] === StabilizerFrame && KeyExistsQ[First[f], "Paulis"];
allPauliStrings[n_] := StringJoin /@ Tuples[{"I", "X", "Y", "Z"}, n];

(* ---- (1) expectation vs dense across many random frames -------------------- *)
expCases = DeleteCases[Flatten[Table[
    Module[{f = buildFrame[n, 4, nT], strs, errs},
        If[! coherentFrameQ[f], Nothing,
            strs = If[n <= 3, allPauliStrings[n], RandomSample[allPauliStrings[n], 12]];
            errs = Table[Abs[N[frameExpectation[f, s] - denseExpect[f, s]]], {s, strs}];
            {n, nT, f["Length"], Max[errs]}]],
    {n, 1, 4}, {nT, 1, 3}, {trial, 1, 6}], 2], Nothing];
expWorst = Max[expCases[[All, 4]]];

(* ---- (2) norm + inner product vs dense (g shares f's reference) ------------ *)
ipCases = DeleteCases[Flatten[Table[
    Module[{f = buildFrame[n, 4, nT], g, vF, vG},
        If[! coherentFrameQ[f], Nothing,
            g = f["T", RandomInteger[{1, n}]];
            vF = Normal @ f["StateVector"]; vG = Normal @ g["StateVector"];
            {n, nT,
             Abs[N[frameNorm2[f] - Conjugate[vF] . vF]],
             Abs[N[frameInnerProduct[f, g] - Conjugate[vF] . vG]]}]],
    {n, 1, 4}, {nT, 1, 2}, {trial, 1, 4}], 2], Nothing];
normWorst = Max[ipCases[[All, 3]]];
ipWorst = Max[ipCases[[All, 4]]];

(* ---- (3) Born consistency: P+ + P- == 1, probabilities real and >= 0 ------- *)
bornCases = DeleteCases[Flatten[Table[
    Module[{f = buildFrame[n, 4, 2], born, marg},
        If[! coherentFrameQ[f], Nothing,
            born = Chop[N @ framePauliBorn[f, StringJoin[ConstantArray["Z", n]]]];
            marg = Chop[N @ frameZmarginal[f, 1]];
            {n, Abs[Total[born] - 1], Abs[Total[marg] - 1], Re[Min[born]], Re[Min[marg]],
             Max[Abs[Im[born]]], Max[Abs[Im[marg]]]}]],
    {n, 1, 4}, {trial, 1, 5}], 1], Nothing];
bornSumWorst = Max[bornCases[[All, 2]], bornCases[[All, 3]]];
bornMinProb = Min[bornCases[[All, 4]], bornCases[[All, 5]]];
bornImWorst = Max[bornCases[[All, 6]], bornCases[[All, 7]]];

WriteString[$Output, ToString[#, InputForm], "\n"] & @ <|
    "expectation_cases" -> Length[expCases],
    "expectation_worst_abs_err" -> expWorst,
    "norm_worst_abs_err" -> normWorst,
    "innerproduct_cases" -> Length[ipCases],
    "innerproduct_worst_abs_err" -> ipWorst,
    "born_cases" -> Length[bornCases],
    "born_sum_to_one_worst" -> bornSumWorst,
    "born_min_probability" -> bornMinProb,
    "born_imag_residue_worst" -> bornImWorst
|>
