PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
b[e_, t_:300] := TimeConstrained[MemoryConstrained[e, 6*^9, "OOM"], t, "TIMEOUT"];
P[l_, x_] := Print["@@ ", StringPadRight[l, 50], " ", ToString[x, InputForm]];
N0 = 12;
a2 = AnnihilationOperator[N0, {2}]; sm = QuantumOperator[{{0, 1}, {0, 0}}, {1}];
H = 2 Pi 0.25 ((sm["Dagger"] @ a2) + (sm @ a2["Dagger"]));
psi0 = QuantumTensorProduct[QuantumState["1"], FockState[0, N0]];
P["H order / dims", {H["Order"], H["Dimensions"]}];
P["jump a on {2}: order / dims", {a2["Order"], a2["Dimensions"]}];
P["jump a on {2}: Matrix head", Head[a2["Matrix"]]];
(* the same jump operator written on the full space *)
aFull = QuantumTensorProduct[QuantumOperator[IdentityMatrix[2], {1}], a2];
P["jump I(x)a: order / dims", {aFull["Order"], aFull["Dimensions"]}];
Do[Module[{tm, rho},
   {tm, rho} = AbsoluteTiming[b[Quiet @ QuantumEvolve[H, {e[[2]]} -> {0.05}, psi0, {\[FormalT], 0, 20}]]];
   P[e[[1]] <> ": head, seconds", {Head[rho], Round[tm, 0.01]}];
   If[MatchQ[rho, _QuantumState],
     P["  P(excited) at t = 0, 1, 2, 4", Table[Chop @ Re @ Total @ Diagonal[Normal[(rho /. \[FormalT] -> t)["DensityMatrix"]][[N0 + 1 ;; 2 N0, N0 + 1 ;; 2 N0]]], {t, {0, 1, 2, 4}}]]]],
  {e, {{"jump on subsystem {2}", a2}, {"jump extended to full space", aFull}}}];
Print[];
ev = b[Quiet @ QuantumEvolve[QuantumOperator["X"], QuantumState["0"], {\[FormalT], 0, 3},
   "AdditionalEquations" -> {WhenEvent[Abs[Indexed[\[FormalS][\[FormalT]], 1]]^2 < 1/2, "StopIntegration"]}], 60];
ifs = Cases[{ev["State"], ev["StateVector"], ev}, _InterpolatingFunction, Infinity, Heads -> True];
P["event: InterpolatingFunctions found", Length[ifs]];
P["  their domains (expect end at Pi/4 = 0.7854)", Union[#["Domain"] & /@ ifs]];
noev = b[Quiet @ QuantumEvolve[QuantumOperator["X"], QuantumState["0"], {\[FormalT], 0, 3}], 60];
P["  control without the event: domains", Union[#["Domain"] & /@ Cases[{noev["State"], noev}, _InterpolatingFunction, Infinity, Heads -> True]]];
