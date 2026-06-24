#!/usr/bin/env wolframscript
(* Clean densification baseline: cost of materializing the 2^n state vector,      *)
(* the unavoidable first step to read ANY observable off a frame today, vs the    *)
(* rank-path frameExpectation. Few-T (chi = 8) frames, varying n.                 *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Stabilizer/stabilizerframe-observables-scripts/framealgebra.wls"];

SeedRandom[11];
buildFixedT[n_, nT_] := Module[{cliff, ts},
    cliff = Table[Module[{g = RandomChoice[{"H", "S", "V", "CNOT", "CZ"}]},
        If[MemberQ[{"CNOT", "CZ"}, g], g -> RandomSample[Range[n], 2], g -> RandomInteger[{1, n}]]], n];
    ts = Table["T" -> RandomInteger[{1, n}], nT];
    Fold[#1[#2] &, PauliStabilizer[n], Join[Take[cliff, 2], ts, Drop[cliff, 2]]]];
obsString[n_] := StringJoin @ ReplacePart[ConstantArray["I", n], {1 -> "X", 2 -> "Z"}];

rows = Table[
    Module[{f = buildFixedT[n, 3], s = obsString[n], tRank, tDense},
        tRank = Min @ Table[ClearSystemCache[]; First @ AbsoluteTiming @ frameExpectation[f, s], 3];
        tDense = Min @ Table[ClearSystemCache[]; First @ AbsoluteTiming[Normal @ f["StateVector"]], 2];
        {n, f["Length"], Round[1000 tRank, 0.01], Round[1000 tDense, 0.01]}],
    {n, {4, 6, 8, 10, 12, 14}}];

WriteString[$Output, ToString[#, InputForm], "\n"] & @ <|
    "rows_n_chi_rankMs_densifyMs" -> rows
|>
