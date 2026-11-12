PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
b[e_, t_:300] := TimeConstrained[MemoryConstrained[e, 6*^9, "OOM"], t, "TIMEOUT"];
P[l_, x_] := Print["@@ ", StringPadRight[l, 44], " ", ToString[x, InputForm]];
(* driven-dissipative Kerr oscillator (QuTiP's standard): H = -D n + K/2 a+a+aa + F(a+a+), L = sqrt(k) a *)
kerrH[n_] := Module[{a = AnnihilationOperator[n], ad}, ad = a["Dagger"];
   {-0.5 (ad @ a) + 0.25 (ad @ ad @ a @ a) + 1.0 (a + ad), a}];
Do[Module[{H, L, tb, te, rho, nT},
    {tb, {H, L}} = AbsoluteTiming[kerrH[n]];
    {te, rho} = AbsoluteTiming[b[QuantumEvolve[H, {L} -> {0.2}, FockState[0, n], {\[FormalT], 0, 10}]]];
    nT = If[MatchQ[rho, _QuantumState],
       b[Chop @ Re @ Tr[Normal[(rho /. \[FormalT] -> 10)["DensityMatrix"]] . Normal[(L["Dagger"] @ L)["Matrix"]]], 60], rho];
    P["Kerr N=" <> ToString[n] <> " build / evolve / <n>(10)", {Round[tb, 0.01], If[NumericQ[te], Round[te, 0.01], te], nT}]],
  {n, {10, 20, 30, 40}}];

Print[];
(* steady state from the Liouvillian null space, no dedicated solver needed *)
Module[{n = 20, H, L, liou, m, ns, rhoSS, tS},
  {H, L} = kerrH[n];
  {tS, ns} = AbsoluteTiming[b[Module[{lv = QuantumOperator["Hamiltonian"[H, {L}, {{0.2}}]]},
       NullSpace[N @ Normal @ lv["MatrixRepresentation"]]], 120]];
  P["steady state via Liouvillian null space N=20 (s)", Round[tS, 0.01]];
  P["  null-space dimension", If[ListQ[ns], Length[ns], ns]]];

Print[];
(* arbitrary precision: Rabi P0(1) = Cos[1]^2 to 30 digits *)
Module[{psi, p, ex = Cos[1]^2},
  psi = b[QuantumEvolve[QuantumOperator["X"], QuantumState["0"], {\[FormalT], 0, 1},
        WorkingPrecision -> 40, PrecisionGoal -> 30, AccuracyGoal -> 30], 120];
  p = If[MatchQ[psi, _QuantumState], Abs[Normal[psi["StateVector"] /. \[FormalT] -> 1][[1]]]^2, psi];
  P["WorkingPrecision 40: P0(1)", p];
  P["  correct digits vs Cos[1]^2", If[NumericQ[p], Round[-Log10[Abs[p - ex]/ex]], "n/a"]];
  P["  machine-precision run correct digits", Module[{q = Abs[Normal[QuantumEvolve[QuantumOperator["X"], QuantumState["0"], {\[FormalT], 0, 1}]["StateVector"] /. \[FormalT] -> 1][[1]]]^2},
       Round[-Log10[Abs[q - ex]/ex]]]]];
