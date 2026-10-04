PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]]; Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
P[l_, x_] := (Print["@@ ", StringPadRight[l, 58], " ", ToString[Short[x, 2], InputForm]]);
kerrH[n_] := Module[{a = AnnihilationOperator[n], ad}, ad = a["Dagger"];
   {-0.5 (ad @ a) + 0.25 (ad @ ad @ a @ a) + 1.0 (a + ad), a}];
n = 40; {H, L} = kerrH[n]; psi0 = FockState[0, n];
{eqs, ret, par} = QuantumEvolve[H, {L} -> {0.2}, psi0, {\[FormalT], 0, 10}, "ReturnEquations" -> True];
t = First[par];
P["1 NDSolveValue, plain", First @ AbsoluteTiming[NDSolveValue[eqs, ret, par]]];
P["2 + StepMonitor :> (time = t)", First @ AbsoluteTiming[Module[{time}, NDSolveValue[eqs, ret, par, StepMonitor :> (time = t)]]]];
P["3 + Progress`EvaluateWithProgress around 2", First @ AbsoluteTiming[Module[{time = 0.},
   Progress`EvaluateWithProgress[NDSolveValue[eqs, ret, par, StepMonitor :> (time = t)],
     <|"Text" -> "Integrating", "Progress" :> time/10, "Percentage" :> time/10, "ElapsedTime" -> Automatic, "RemainingTime" -> Automatic|>]]]];
P["4 QuantumEvolve, ReturnSolution -> True (no state built)", First @ AbsoluteTiming[QuantumEvolve[H, {L} -> {0.2}, psi0, {\[FormalT], 0, 10}, "ReturnSolution" -> True]]];
P["5 QuantumEvolve, full", First @ AbsoluteTiming[QuantumEvolve[H, {L} -> {0.2}, psi0, {\[FormalT], 0, 10}]]];
P["6 QuantumEvolve, MergeInterpolatingFunctions -> False", First @ AbsoluteTiming[QuantumEvolve[H, {L} -> {0.2}, psi0, {\[FormalT], 0, 10}, "MergeInterpolatingFunctions" -> False]]];
(* is the overhead a fixed wall-clock cost after a threshold? same check at N=30 *)
{H3, L3} = kerrH[30];
P["N=30: 3-style progress-wrapped solve", Module[{e, r, p, time = 0.},
   {e, r, p} = QuantumEvolve[H3, {L3} -> {0.2}, FockState[0, 30], {\[FormalT], 0, 10}, "ReturnEquations" -> True];
   First @ AbsoluteTiming[Progress`EvaluateWithProgress[NDSolveValue[e, r, p, StepMonitor :> (time = First[p])], <|"Progress" :> time/10|>]]]];
