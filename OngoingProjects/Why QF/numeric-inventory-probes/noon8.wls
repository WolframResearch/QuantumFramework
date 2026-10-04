PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]]; Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
P[l_, x_] := (Print["@@ ", StringPadRight[l, 52], " ", ToString[Short[x, 2], InputForm]]);
qfi[rho_, g_] := Module[{l, v, m},
  {l, v} = Eigensystem[(rho + ConjugateTranspose[rho])/2]; l = Chop[Re[l]];
  m = Conjugate[v] . g . Transpose[v];
  2 Total[Flatten @ Table[If[l[[i]] + l[[j]] > 1*^-12, (l[[i]] - l[[j]])^2/(l[[i]] + l[[j]]) Abs[m[[i, j]]]^2, 0],
     {i, Length[l]}, {j, Length[l]}]]];
noon[n_, eta_] := Module[{d = n + 1, a1, a2, id, a1F, a2F, psi, H, rho, r, n2},
  id = QuantumOperator[IdentityMatrix[d], {1}];
  a1 = AnnihilationOperator[d, {1}]; a2 = AnnihilationOperator[d, {2}];
  a1F = QuantumTensorProduct[a1, QuantumOperator[IdentityMatrix[d], {2}]];
  a2F = QuantumTensorProduct[QuantumOperator[IdentityMatrix[d], {1}], a2];
  psi = QuantumState[(UnitVector[d^2, n d + 1] + UnitVector[d^2, n + 1])/Sqrt[2.], {d, d}];
  H = 0. (a1F["Dagger"] @ a1F);
  rho = QuantumEvolve[H, {a1F, a2F} -> {1., 1.}, psi, {\[FormalT], 0, -Log[eta]}];
  r = Normal[(rho /. \[FormalT] -> -Log[eta])["DensityMatrix"]];
  n2 = Normal[(a2F["Dagger"] @ a2F)["Matrix"]];
  {qfi[r, n2], n^2 eta^n}];
P["retry lossy NOON n=8 eta=0.9 {QFI, n^2 eta^n}", First @ {noon[8, 0.9]}];
P["  memory in use after (GB)", N[MemoryInUse[]/2^30]];
