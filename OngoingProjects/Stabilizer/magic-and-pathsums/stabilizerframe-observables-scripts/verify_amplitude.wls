#!/usr/bin/env wolframscript
(* Amplitude / computational-Born verification vs dense, on coherent frames.     *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Stabilizer/stabilizerframe-observables-scripts/framealgebra.wls"];

SeedRandom[7];

randomSpecs[n_, nClifford_, nT_] := Module[{gateset, cliff, ts},
    gateset = If[n >= 2, {"H", "S", "V", "X", "Y", "Z", "CNOT", "CZ", "SWAP"}, {"H", "S", "V", "X", "Y", "Z"}];
    cliff = Table[Module[{g = RandomChoice[gateset]},
        If[MemberQ[{"CNOT", "CZ", "SWAP"}, g], g -> RandomSample[Range[n], 2], g -> RandomInteger[{1, n}]]],
        nClifford];
    ts = Table["T" -> RandomInteger[{1, n}], nT];
    Join[Take[cliff, Min[2, Length[cliff]]], ts, Drop[cliff, Min[2, Length[cliff]]]]];
buildFrame[n_, nClifford_, nT_] := Fold[#1[#2] &, PauliStabilizer[n], randomSpecs[n, nClifford, nT]];
coherentFrameQ[f_] := Head[f] === StabilizerFrame && KeyExistsQ[First[f], "Paulis"];
idx[x_] := FromDigits[x, 2] + 1;

(* (0) sanity: reference single-amplitude ratio vs dense ref vector ----------- *)
refRatioCases = DeleteCases[Flatten[Table[
    Module[{f = buildFrame[n, 4, 1], ref, vref, x0, errs},
        If[! coherentFrameQ[f], Nothing,
            ref = First[f]["Components"][[1, 2]];
            vref = Normal @ ref["State"]["StateVector"];
            x0 = stabSupportPoint[ref];
            errs = Table[Module[{x = IntegerDigits[b, 2, n], pred, act},
                pred = stabRatio[ref, x, x0] vref[[idx[x0]]];
                act = vref[[idx[x]]];
                Abs[N[pred - act]]], {b, 0, 2^n - 1}];
            Max[errs]]],
    {n, 1, 5}, {trial, 1, 4}], 1], Nothing];
refRatioWorst = Max[refRatioCases];

(* (1) full complex amplitude <x|psi> vs frame["StateVector"] ----------------- *)
ampCases = DeleteCases[Flatten[Table[
    Module[{f = buildFrame[n, 4, nT], ref, vref, x0, anchor, v, errs},
        If[! coherentFrameQ[f], Nothing,
            ref = First[f]["Components"][[1, 2]];
            vref = Normal @ ref["State"]["StateVector"];
            x0 = stabSupportPoint[ref];
            anchor = vref[[idx[x0]]];            (* global-phase pin (kernel: CH-form) *)
            v = Normal @ f["StateVector"];
            errs = Table[Module[{x = IntegerDigits[b, 2, n]},
                Abs[N[frameAmplitude[f, x, anchor, x0] - v[[idx[x]]]]]], {b, 0, 2^n - 1}];
            {n, nT, Max[errs]}]],
    {n, 1, 4}, {nT, 1, 3}, {trial, 1, 3}], 2], Nothing];
ampWorst = Max[ampCases[[All, 3]]];

(* (2) computational Born P(x) vs |amp|^2, and Sum_x P(x) == norm2 ------------ *)
probCases = DeleteCases[Flatten[Table[
    Module[{f = buildFrame[n, 4, nT], v, probs, errs, total, n2},
        If[! coherentFrameQ[f], Nothing,
            v = Normal @ f["StateVector"];
            probs = Table[frameCompProb[f, IntegerDigits[b, 2, n]], {b, 0, 2^n - 1}];
            errs = Table[Abs[N[probs[[b + 1]] - Abs[v[[b + 1]]]^2]], {b, 0, 2^n - 1}];
            total = Total[probs];
            n2 = frameNorm2[f];
            {n, nT, Max[errs], Abs[N[total - n2]], Re[Min[N[probs]]]}]],
    {n, 1, 4}, {nT, 1, 2}, {trial, 1, 3}], 2], Nothing];
probWorst = Max[probCases[[All, 3]]];
probTotalWorst = Max[probCases[[All, 4]]];
probMin = Min[probCases[[All, 5]]];

WriteString[$Output, ToString[#, InputForm], "\n"] & @ <|
    "ref_ratio_cases" -> Length[refRatioCases], "ref_ratio_worst_err" -> refRatioWorst,
    "amplitude_cases" -> Length[ampCases], "amplitude_worst_err" -> ampWorst,
    "compprob_cases" -> Length[probCases], "compprob_worst_err" -> probWorst,
    "compprob_total_vs_norm2_worst" -> probTotalWorst, "compprob_min" -> probMin
|>
