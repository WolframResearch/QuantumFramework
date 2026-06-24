#!/usr/bin/env wolframscript
(* Scaling: rank-based frameExpectation (chi^2 * poly(n)) vs dense densification *)
(* (2^n). Plus an n-independent correctness anchor needing no dense reference.   *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Stabilizer/stabilizerframe-observables-scripts/framealgebra.wls"];

SeedRandom[11];

pauliMat[str_String] := Module[{mats = Characters[str] /. {"I" -> PauliMatrix[0],
        "X" -> PauliMatrix[1], "Y" -> PauliMatrix[2], "Z" -> PauliMatrix[3]}},
    If[Length[mats] == 1, First[mats], KroneckerProduct @@ mats]];

(* fixed-T frame: a depth-n random Clifford layer + nT T-gates (chi = 2^nT). *)
buildFixedT[n_, nT_] := Module[{cliff, ts},
    cliff = Table[Module[{g = RandomChoice[{"H", "S", "V", "CNOT", "CZ"}]},
        If[MemberQ[{"CNOT", "CZ"}, g], g -> RandomSample[Range[n], 2], g -> RandomInteger[{1, n}]]], n];
    ts = Table["T" -> RandomInteger[{1, n}], nT];
    Fold[#1[#2] &, PauliStabilizer[n], Join[Take[cliff, 2], ts, Drop[cliff, 2]]]];

obsString[n_] := StringJoin @ ReplacePart[ConstantArray["I", n], {1 -> "X", 2 -> "Z"}];

minTime[expr_Hold] := Min @ Table[ClearSystemCache[]; First @ AbsoluteTiming @ ReleaseHold[expr], 3];

(* rank-only scaling (no dense), few-T many-qubit *)
rankRows = Table[
    Module[{f = buildFixedT[n, 3], s, t},
        s = obsString[n];
        t = Min @ Table[ClearSystemCache[]; First @ AbsoluteTiming @ frameExpectation[f, s], 3];
        {n, f["Length"], Round[1000 t, 0.01]}],     (* {n, chi, rank ms} *)
    {n, {4, 6, 8, 10, 12, 16, 20, 30, 40}}];

(* dense densification cost for the same observable (capped where 2^n is sane) *)
denseRows = Table[
    Module[{f = buildFixedT[n, 3], s, M, t},
        s = obsString[n];
        M = pauliMat[s];
        t = Min @ Table[ClearSystemCache[];
            First @ AbsoluteTiming @ Module[{v = Normal @ f["StateVector"]}, Conjugate[v] . (M . v)], 2];
        {n, Round[1000 t, 0.01]}],
    {n, {4, 6, 8, 10, 12, 14}}];

(* n-independent correctness anchor: qubits 1,2,3 = T|+>, rest |0>. No dense.    *)
(* <Y1 Y2 Y3> = (1/Sqrt2)^3, <Z1> = 0, <Z4> = +1, norm = 1.                      *)
anchorRows = Table[
    Module[{f = Fold[#1[#2] &, PauliStabilizer[n],
            {"H" -> 1, "H" -> 2, "H" -> 3, "T" -> 1, "T" -> 2, "T" -> 3}],
        sYYY, sZ1, sZ4},
        sYYY = StringJoin @ ReplacePart[ConstantArray["I", n], {1 -> "Y", 2 -> "Y", 3 -> "Y"}];
        sZ1 = StringJoin @ ReplacePart[ConstantArray["I", n], 1 -> "Z"];
        sZ4 = StringJoin @ ReplacePart[ConstantArray["I", n], 4 -> "Z"];
        {n, f["Length"],
         Chop[N[frameExpectation[f, sYYY] - (1/Sqrt[2])^3]],
         Chop[N[frameExpectation[f, sZ1]]],
         Chop[N[frameExpectation[f, sZ4] - 1]],
         Chop[N[frameNorm2[f] - 1]]}],
    {n, {6, 10, 20, 40, 60}}];
anchorWorst = Max[Abs[Flatten[anchorRows[[All, 3 ;; 6]]]]];

WriteString[$Output, ToString[#, InputForm], "\n"] & @ <|
    "rank_rows_n_chi_ms" -> rankRows,
    "dense_rows_n_ms" -> denseRows,
    "anchor_rows_n_chi_dYYY_Z1_dZ4_dnorm" -> anchorRows,
    "anchor_worst_abs_err" -> anchorWorst
|>
