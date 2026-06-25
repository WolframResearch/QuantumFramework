p[l_, e_] := WriteString[$Output, l <> " = " <> ToString[InputForm[e]] <> "\n"];
\[Sigma] = Table[PauliMatrix[j], {j, 3}];
densityMatrix[r_] := 1/2 (IdentityMatrix[2] + r . \[Sigma]);

(* C1 projective measurement of Z on |+> *)
\[Rho] = densityMatrix[{1, 0, 0}];
proj = TensorProduct[#, Conjugate@#] & /@ {{1, 0}, {0, 1}};
p["C1-probs", Tr[# . \[Rho]] & /@ proj];
p["C1-eigs", Eigenvalues[PauliMatrix[3]]];
p["C1-mean", Tr[PauliMatrix[3] . \[Rho]]];

(* C2 measure |0> in the X basis *)
xpm = {{1, 1}, {1, -1}}/Sqrt[2];
p["C2-probs", Abs[Conjugate[xpm] . {1, 0}]^2];

(* C3 post-measurement (collapsed) branches *)
\[Psi] = {1, 1}/Sqrt[2];
P0 = TensorProduct[{1, 0}, {1, 0}]; P1 = TensorProduct[{0, 1}, {0, 1}];
p["C3-branches", {P0 . \[Psi], P1 . \[Psi]}];

(* C4 SIC-POVM (tetrahedron) *)
verts = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}}/Sqrt[3];
effects = (1/4) (IdentityMatrix[2] + # . \[Sigma]) & /@ verts;
p["C4-complete", Total[effects] == IdentityMatrix[2]];
p["C4-psd", And @@ (AllTrue[Eigenvalues[#], # >= 0 &] & /@ effects)];
p["C4-probs", Tr[# . densityMatrix[{0, 0, 1}]] & /@ effects];
p["C4-total", Total[Tr[# . densityMatrix[{0, 0, 1}]] & /@ effects]];

(* C5 finite-shot statistics vs Born rule *)
SeedRandom[42];
p["C5-counts", Counts@RandomChoice[{1/2, 1/2} -> {0, 1}, 2000]];
p["C5-theory", N@{1/2, 1/2}];

(* C6 correlated measurement of a Bell pair *)
bell = {1, 0, 0, 1}/Sqrt[2];
p["C6-joint", Abs[bell]^2];
