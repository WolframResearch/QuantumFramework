p[l_, e_] := WriteString[$Output, l <> " = " <> ToString[InputForm[e]] <> "\n"];
\[Sigma] = Table[PauliMatrix[j], {j, 3}];
reduceA[\[Rho]_] := TensorContract[ArrayReshape[\[Rho], {2, 2, 2, 2}], {{2, 4}}];
reduceB[\[Rho]_] := TensorContract[ArrayReshape[\[Rho], {2, 2, 2, 2}], {{1, 3}}];
partialT[\[Rho]_] := ArrayReshape[Transpose[ArrayReshape[\[Rho], {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}];
vne[\[Rho]_] := -Total[# Log2[#] & /@ Select[Re@Eigenvalues[\[Rho]], # > 10^-12 &]];

(* D1 tensor products *)
p["D1-state", KroneckerProduct[{1, 0}, {0, 1}]];
p["D1-op", KroneckerProduct[PauliMatrix[1], PauliMatrix[3]]];

(* D2 Bell / GHZ / W *)
bell = {1, 0, 0, 1}/Sqrt[2];
ghz = {1, 0, 0, 0, 0, 0, 0, 1}/Sqrt[2];
w = {0, 1, 1, 0, 1, 0, 0, 0}/Sqrt[3];
p["D2-bell", bell]; p["D2-ghz", ghz]; p["D2-w", w];

(* D3 reduced density matrices *)
\[Rho]bell = KroneckerProduct[bell, Conjugate@bell];
p["D3-reduceA", reduceA[\[Rho]bell]];

(* D4 entanglement via Schmidt rank *)
schmidtRank[\[Psi]_, {dA_, dB_}] := Length@Select[SingularValueList[ArrayReshape[\[Psi], {dA, dB}]], # > 10^-10 &];
p["D4-bell", schmidtRank[bell, {2, 2}]];
p["D4-prod", schmidtRank[KroneckerProduct[{1, 0}, {0, 1}], {2, 2}]];

(* D5 entropy, concurrence, negativity *)
p["D5-entropy", vne[reduceA[\[Rho]bell]]];
concurrence[\[Psi]_] := Abs[Conjugate[\[Psi]] . KroneckerProduct[PauliMatrix[2], PauliMatrix[2]] . Conjugate[\[Psi]]];
p["D5-conc-bell", concurrence[bell]];
p["D5-conc-prod", concurrence[KroneckerProduct[{1, 0}, {0, 1}]]];
p["D5-neg-bell", Total@Abs@Select[Eigenvalues[partialT[\[Rho]bell]], # < 0 &]];

(* D6 Schmidt coefficients *)
p["D6-schmidt", SingularValueList[ArrayReshape[bell, {2, 2}]]];

(* D7 no-signalling: local X on qubit 1 leaves reduced state of qubit 2 unchanged *)
XI = KroneckerProduct[PauliMatrix[1], IdentityMatrix[2]];
p["D7-before", reduceB[\[Rho]bell]];
p["D7-after", reduceB[XI . \[Rho]bell . ConjugateTranspose[XI]]];
p["D7-equal", reduceB[\[Rho]bell] == reduceB[XI . \[Rho]bell . ConjugateTranspose[XI]]];
