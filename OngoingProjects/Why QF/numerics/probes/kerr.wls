PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]]; Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
P[l_, x_] := (Print["@@ ", StringPadRight[l, 52], " ", ToString[x, InputForm]]);
kerrH[n_] := Module[{a = AnnihilationOperator[n], ad}, ad = a["Dagger"];
   {-0.5 (ad @ a) + 0.25 (ad @ ad @ a @ a) + 1.0 (a + ad), a}];
Do[Module[{H, L, te, rho, steps = 0},
    {H, L} = kerrH[n];
    {te, rho} = AbsoluteTiming[QuantumEvolve[H, {L} -> {0.2}, FockState[0, n], {\[FormalT], 0, 10}]];
    P["Kerr Lindblad N=" <> ToString[n] <> " evolve (s)", te]],
  {n, {30, 35, 40, 45, 50}}];
Module[{H, L, eqs, t1, t2},
  {H, L} = kerrH[40];
  {t1, eqs} = AbsoluteTiming[QuantumEvolve[H, {L} -> {0.2}, FockState[0, 40], {\[FormalT], 0, 10}, "ReturnEquations" -> True]];
  P["N=40 building the equations only (s)", t1];
  P["  RHS leaf count (size of the ODE system expression)", LeafCount[eqs[[1]]]];
  {t2, _} = AbsoluteTiming[NDSolveValue[eqs[[1]], eqs[[2]], eqs[[3]]]];
  P["  NDSolveValue on those equations directly (s)", t2]];
