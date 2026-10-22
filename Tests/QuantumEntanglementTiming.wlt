(* ::Package:: *)

BeginTestSection["QuantumEntanglementTiming"]


(* Speed rule of the entanglement spec (ENTANGLEMENT-SPEC-DRAFT.md, section 6). Each case is timed as the
   minimum of 7 runs with ClearSystemCache[] before each and the QF property cache off, so a repeated run
   cannot read a cached result. The baseline times were recorded once on 96f48907 into
   QuantumEntanglementBaseline.wl, next to this file.
   - a case whose baseline answer is correct may take at most 1.25 times its baseline time (+1 ms);
   - a case whose baseline answer is wrong or does not finish has its own limit, "Limit" in the baseline file. *)


(* ========== Cases ========== *)

proj[v_] := Outer[Times, v, Conjugate[v]];
randomMixed[d_] := With[{m = RandomComplex[{-1 - I, 1 + I}, {d, d}]}, # / Tr[#] & [m . ConjugateTranspose[m]]];
randomKet[d_] := Normalize[RandomComplex[{-1 - I, 1 + I}, d]];
singlet = {0, 1, -1, 0} / Sqrt[2];
werner[p_] := p proj[singlet] + (1 - p) IdentityMatrix[4] / 4;
blochDM[{x_, y_, z_}] := (IdentityMatrix[2] + x PauliMatrix[1] + y PauliMatrix[2] + z PauliMatrix[3]) / 2;
exactSep22 = 1/3 KroneckerProduct[blochDM[{1/2, 0, 1/3}], blochDM[{0, 1/5, 1/2}]] +
    2/3 KroneckerProduct[blochDM[{1/7, 1/3, 0}], blochDM[{1/2, 1/4, -1/3}]];
w3m = Normal @ QuantumState["W"[3]]["DensityMatrix"];
g3m = Normal @ QuantumState["GHZ"[3]]["DensityMatrix"];
isotropic3 = With[{phi = Flatten[IdentityMatrix[3]] / Sqrt[3]}, proj[phi] / 2 + IdentityMatrix[9] / 18];   (* entangled: p = 1/2 > 1/4 *)

SeedRandom[31];
ket10 = Flatten @ KroneckerProduct[randomKet[2^9], {1, 0}];      (* qubits 1-9 pure, qubit 10 in |0> *)
ket11 = Flatten @ KroneckerProduct[randomKet[2^10], {1, 0}];
ket10random = randomKet[2^10];
ket32 = randomKet[32];
mixed24 = randomMixed[8];
mixed33 = randomMixed[9];
mixed44 = randomMixed[16];

shippedFile = FileNameJoin[{ParentDirectory[First[PacletFind["Wolfram/QuantumFramework"]]["Location"]], "Tests", "QuantumEntanglement.wlt"}];

timingCases = <|
    "G1" -> Hold[TestReport[shippedFile]],
    "G2" -> Hold[QuantumEntanglementMonotone[QuantumState[ket10, ConstantArray[2, 10]], {Range[4], Range[5, 9]}, "Concurrence"]],
    "G3" -> Hold[QuantumEntanglementMonotone[QuantumState[ket11, ConstantArray[2, 11]], {Range[5], Range[6, 10]}, "Concurrence"]],
    "G4" -> Hold[QuantumEntanglementMonotone[QuantumState[QuantumState[ket32, 32], QuantumBasis["Fourier"[32]]], Automatic, "Negativity"]],
    "G5" -> Hold[With[{qs = QuantumState[QuantumState[Normalize[{1, 2, 3, 4}], {2, 2}], QuantumBasis["PauliX", 2]]},
        {QuantumEntanglementMonotone[qs, "EntanglementEntropy"], QuantumEntanglementMonotone[qs, "Concurrence"]}]],
    "G6" -> Hold[{QuantumEntangledQ[QuantumState["PhiPlus"]], QuantumEntangledQ[QuantumState[N @ werner[4/5], {2, 2}]]}],
    "G7" -> Hold[{
        QuantumEntangledQ[QuantumState[mixed24, {2, 4}]], QuantumEntangledQ[QuantumState[mixed33, {3, 3}]], QuantumEntangledQ[QuantumState[mixed44, {4, 4}]],
        QuantumEntanglementMonotone[QuantumState[mixed33, {3, 3}], "Negativity"], QuantumEntanglementMonotone[QuantumState[mixed44, {4, 4}], "Negativity"]}],
    "G8" -> Hold[QuantumEntanglementMonotone[QuantumState[mixed33, {3, 3}], "Concurrence"]],
    "G9a" -> Hold[QuantumEntangledQ[QuantumState[ket10random, 2], {{1}, {2}}]],
    "G9b" -> Hold[QuantumEntangledQ[QuantumState[ket10random, 2], {1}]],
    "G9c" -> Hold[QuantumEntangledQ[QuantumState[ket10random, 2], {Range[5], Range[6, 10]}]],
    "G10" -> Hold[{QuantumEntanglementMonotone[QuantumState["GHZ"[3]], {{1}, {2}}, "Concurrence"], QuantumEntanglementMonotone[QuantumState["W"[3]], {{1}, {3}}, "Concurrence"]}],
    "G11" -> Hold[QuantumEntanglementMonotone[QuantumState[QuantumState[(w3m + g3m) / 2, {2, 2, 2}], QuantumBasis["PauliX", 3]], Automatic, "Negativity"]],
    "G12" -> Hold[Block[{a, b, c}, QuantumEntanglementMonotone[QuantumState[{0, a, b, 0, c, 0, 0, 0}, {2, 2, 2}], {{1}, {2}}, "Concurrence"]]],
    "G13" -> Hold[QuantumEntanglementMonotone[QuantumState[exactSep22, {2, 2}], "Concurrence"]],
    "G14" -> Hold[QuantumEntangledQ[QuantumState[QuantumState[isotropic3, {3, 3}], QuantumBasis[QuantumTensorProduct[QuditBasis["Fourier"[3]], QuditBasis["Fourier"[3]]]]]]],
    "G15" -> Hold[QuantumEntanglementMonotone[QuantumState[proj[Normalize @ Flatten @ KroneckerProduct[{1, 2^(1/3), I}, {Sqrt[2], 1 - I, 3^(1/4)}]], {3, 3}], "Negativity"]]
|>;

SetAttributes[timeMin, HoldFirst];
timeMin[Hold[expr_], n_] := Block[{Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = False},
    Min @ Table[ClearSystemCache[]; First @ AbsoluteTiming[expr], n]
];

baselineFile = FileNameJoin[{ParentDirectory[First[PacletFind["Wolfram/QuantumFramework"]]["Location"]], "Tests", "QuantumEntanglementBaseline.wl"}];
baseline = If[FileExistsQ[baselineFile], Get[baselineFile], <||>];
limitFor[id_] := With[{b = Lookup[baseline, id, <||>]},
    If[Lookup[b, "BaselineAnswer"] === "correct", Max[1.25 b["Seconds"], b["Seconds"] + 0.001], Lookup[b, "Limit", 0]]
];
(* the seven runs together may take at most seven times the limit, so a case that never finishes fails instead of hanging *)
withinLimit[id_] := With[{h = timingCases[id], limit = limitFor[id]},
    TimeConstrained[timeMin[h, 7], 7 limit + 1, Infinity] <= limit
];


(* ========== Tests ========== *)

VerificationTest[withinLimit["G1"], True, {___}, TestID -> "G1-ShippedTestFile"]
VerificationTest[withinLimit["G2"], True, {___}, TestID -> "G2-PurePair10Qubits"]
VerificationTest[withinLimit["G3"], True, {___}, TestID -> "G3-PurePair11Qubits"]
VerificationTest[withinLimit["G4"], True, {___}, TestID -> "G4-PureStateInFourier32Basis"]
VerificationTest[withinLimit["G5"], True, {___}, TestID -> "G5-ExactPureTwoQubitsInPauliXBasis"]
VerificationTest[withinLimit["G6"], True, {___}, TestID -> "G6-BellAndWernerDefault"]
VerificationTest[withinLimit["G7"], True, {___}, TestID -> "G7-RandomMixedDefaultAndNegativity"]
VerificationTest[withinLimit["G8"], True, {___}, TestID -> "G8-RandomMixed3x3Concurrence"]
VerificationTest[withinLimit["G9a"], True, {___}, TestID -> "G9a-Random10QubitsPairSplit"]
VerificationTest[withinLimit["G9b"], True, {___}, TestID -> "G9b-Random10QubitsOneVsRest"]
VerificationTest[withinLimit["G9c"], True, {___}, TestID -> "G9c-Random10QubitsFiveVsFive"]
VerificationTest[withinLimit["G10"], True, {___}, TestID -> "G10-GHZAndWPairConcurrence"]
VerificationTest[withinLimit["G11"], True, {___}, TestID -> "G11-ExactMixedStateInPauliXBasis"]
VerificationTest[withinLimit["G12"], True, {___}, TestID -> "G12-SymbolicPairConcurrence"]
VerificationTest[withinLimit["G13"], True, {___}, TestID -> "G13-ExactMixedTwoQubitConcurrence"]
VerificationTest[withinLimit["G14"], True, {___}, TestID -> "G14-ExactIsotropicInFourierBasis"]
VerificationTest[withinLimit["G15"], True, {___}, TestID -> "G15-ExactNestedRadicalProduct"]


EndTestSection[]
