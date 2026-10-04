PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
b[e_, t_:300] := TimeConstrained[MemoryConstrained[e, 6*^9, "OOM"], t, "TIMEOUT"];
P[l_, x_] := Print["@@ ", StringPadRight[l, 46], " ", ToString[x, InputForm]];
alpha = -2 Pi 0.2;
leak[T_, lam_ ? NumericQ] := Module[{bb = AnnihilationOperator[4], bd, sig = T/4., A, ox, oy, H, psi, v},
  bd = bb["Dagger"];
  A = Pi/NIntegrate[Exp[-(t - T/2)^2/(2 sig^2)], {t, 0, T}];
  ox = A Exp[-(\[FormalT] - T/2)^2/(2 sig^2)];
  oy = -lam D[ox, \[FormalT]]/alpha;
  H = (alpha/2) (bd @ bd @ bb @ bb) + (ox/2) (bb + bd) + (oy/2) (I (bd - bb));
  psi = QuantumEvolve[H, FockState[0, 4], {\[FormalT], 0, T}];
  v = Normal[psi["StateVector"] /. \[FormalT] -> T];
  Abs[v[[3]]]^2];
{tm, best} = AbsoluteTiming[b[Quiet @ FindMinimum[leak[8, l], {l, 0.5, 1.0}], 300]];
P["FindMinimum over DRAG lambda at T=8ns", best];
P["  seconds", Round[tm, 0.1]];
Print[];
Module[{N = 12, a, ad, sm, H, rho, tm2},
  a = AnnihilationOperator[N, {2}]; ad = a["Dagger"];
  sm = QuantumOperator[{{0, 1}, {0, 0}}, {1}];
  H = 2 Pi 0.25 ((sm["Dagger"] @ a) + (sm @ ad));
  {tm2, rho} = AbsoluteTiming[b[QuantumEvolve[H, {a} -> {0.05}, QuantumTensorProduct[QuantumState["1"], FockState[0, N]], {\[FormalT], 0, 20}]]];
  P["damped JC qubit x cavity(12) Lindblad (s)", Round[tm2, 0.01]];
  P["  Head / dims", If[MatchQ[rho, _QuantumState], rho["Dimensions"], Short[rho, 1]]];
  If[MatchQ[rho, _QuantumState],
    P["  P(excited) at t = 0, 1, 2, 4", Table[Chop@Re@Total[Normal[(rho /. \[FormalT] -> t)["DensityMatrix"]][[N + 1 ;; 2 N, N + 1 ;; 2 N]] // Diagonal], {t, {0, 1, 2, 4}}]]]];
Print[];
ev = b[Quiet @ QuantumEvolve[QuantumOperator["X"], QuantumState["0"], {\[FormalT], 0, 3},
   "AdditionalEquations" -> {WhenEvent[Abs[Indexed[\[FormalS][\[FormalT]], 1]]^2 < 1/2, "StopIntegration"]}], 60];
P["WhenEvent via AdditionalEquations, Head", Head[ev]];
P["  last time in solution (expect Pi/4 = 0.7854)",
   Quiet @ Check[Max[Cases[ev, InterpolatingFunction[dom_, ___] :> dom[[1, 2]], Infinity]], "n/a"]];
