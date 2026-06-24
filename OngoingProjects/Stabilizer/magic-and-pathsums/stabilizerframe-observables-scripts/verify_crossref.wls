#!/usr/bin/env wolframscript
(* Distinct-reference inner product + shortcoming-3 demonstration.                *)
(*  (1) frameInnerProductGeneralExact (relating-Pauli reduction) == dense.        *)
(*  (2) shared-reference: general == poly shared-ref specialization.              *)
(*  (3) per-pair magnitude via closed form == |dense per-pair| (poly magnitude).  *)
(*  (4) shortcoming 3: naive Sum_i c_i s_i (stored components) != coherent state; *)
(*      the relating-Pauli reconstruction does match.                            *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Stabilizer/stabilizerframe-observables-scripts/framealgebra.wls"];

SeedRandom[99];

randomSpecs[n_, nC_, nT_] := Module[{gs, cl, ts},
    gs = If[n >= 2, {"H", "S", "V", "X", "Y", "Z", "CNOT", "CZ", "SWAP"}, {"H", "S", "V", "X", "Y", "Z"}];
    cl = Table[Module[{g = RandomChoice[gs]},
        If[MemberQ[{"CNOT", "CZ", "SWAP"}, g], g -> RandomSample[Range[n], 2], g -> RandomInteger[{1, n}]]], nC];
    ts = Table["T" -> RandomInteger[{1, n}], nT];
    Join[Take[cl, Min[2, Length[cl]]], ts, Drop[cl, Min[2, Length[cl]]]]];
buildFrame[n_, nC_, nT_] := Fold[#1[#2] &, PauliStabilizer[n], randomSpecs[n, nC, nT]];
coherentFrameQ[f_] := Head[f] === StabilizerFrame && KeyExistsQ[First[f], "Paulis"];
denseInner[fA_, fB_] := Conjugate[Normal @ fA["StateVector"]] . Normal @ fB["StateVector"];

(* (1)+(2) distinct- and shared-reference inner product *)
genCases = DeleteCases[Flatten[Table[
    Module[{fA = buildFrame[n, 4, nT], fB, shared},
        If[! coherentFrameQ[fA], Nothing,
            fB = fA["H", RandomInteger[{1, n}]]["S", RandomInteger[{1, n}]];
            shared = First[fA]["Components"][[1, 2]] === First[fB]["Components"][[1, 2]];
            {n, nT, shared,
             Abs[N[frameInnerProductGeneralExact[fA, fB] - denseInner[fA, fB]]]}]],
    {n, 1, 4}, {nT, 1, 2}, {trial, 1, 4}], 2], Nothing];
genWorst = Max[genCases[[All, 4]]];
allDistinct = AllTrue[genCases, #[[3]] === False &];

sharedConsist = DeleteCases[Flatten[Table[
    Module[{f = buildFrame[n, 4, nT], g},
        If[! coherentFrameQ[f], Nothing,
            g = f["T", RandomInteger[{1, n}]];
            Abs[N[frameInnerProductGeneralExact[f, g] - frameInnerProduct[f, g]]]]],
    {n, 1, 4}, {nT, 1, 2}, {trial, 1, 3}], 2], Nothing];
sharedConsistWorst = Max[sharedConsist];

(* (3) per-pair magnitude: closed form == |dense per-pair| *)
magCases = DeleteCases[Flatten[Table[
    Module[{fA = buildFrame[n, 3, nT], fB, refA, refB, aP, bP, nn, errs},
        If[! coherentFrameQ[fA], Nothing,
            fB = fA["H", RandomInteger[{1, n}]];
            refA = First[fA]["Components"][[1, 2]]; refB = First[fB]["Components"][[1, 2]];
            nn = refA["Qubits"];
            aP = sympFromRelating /@ First[fA]["Paulis"]; bP = sympFromRelating /@ First[fB]["Paulis"];
            errs = Flatten @ Table[Module[{Q = sympTimes[sympDagger[aP[[i]]], bP[[j]]]},
                Abs[N[pairMagClosed[refA, Q, refB] - Abs[pairAmpDense[refA, Q, refB, nn]]]]],
                {i, Length[aP]}, {j, Length[bP]}];
            Max[errs]]],
    {n, 2, 4}, {nT, 1, 2}, {trial, 1, 3}], 2], Nothing];
magWorst = Max[magCases];

(* (4) shortcoming 3: naive stored-component sum vs coherent reconstruction *)
trapCases = DeleteCases[Flatten[Table[
    Module[{f = buildFrame[n, 4, nT], coherent, naive, relPauli, ref, n0, v1},
        If[! coherentFrameQ[f], Nothing,
            coherent = Normal @ f["StateVector"];
            naive = Total[#1 Normal[#2["State"]["StateVector"]] & @@@ First[f]["Components"]];
            ref = First[f]["Components"][[1, 2]]; n0 = ref["Qubits"]; v1 = ref["State"]["StateVector"];
            relPauli = Total @ MapThread[
                #1[[1]] Normal[pauliApplyVec[sympFromRelating[#2], v1, n0]] &,
                {First[f]["Components"], First[f]["Paulis"]}];
            {n, nT,
             Max[Abs[N[naive - coherent]]],        (* naive: expected NONZERO (the trap) *)
             Max[Abs[N[relPauli - coherent]]]}]],   (* relating-Pauli: expected ~0 *)
    {n, 1, 4}, {nT, 1, 2}, {trial, 1, 3}], 2], Nothing];
naiveTrapMax = Max[trapCases[[All, 3]]];
relPauliMax = Max[trapCases[[All, 4]]];

WriteString[$Output, ToString[#, InputForm], "\n"] & @ <|
    "general_inner_cases" -> Length[genCases],
    "general_all_distinct_reference" -> allDistinct,
    "general_inner_worst_err" -> genWorst,
    "shared_consistency_worst" -> sharedConsistWorst,
    "perpair_magnitude_closedform_worst" -> magWorst,
    "shortcoming3_naive_stored_sum_worst (TRAP, expected large)" -> naiveTrapMax,
    "shortcoming3_relating_pauli_reconstruction_worst" -> relPauliMax
|>
