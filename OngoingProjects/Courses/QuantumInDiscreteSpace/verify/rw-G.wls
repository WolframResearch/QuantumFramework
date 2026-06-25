p[l_, e_] := WriteString[$Output, l <> " = " <> ToString[InputForm[e]] <> "\n"];
ketd[d_, ks__] := Fold[Flatten@KroneckerProduct[#1, #2] &, UnitVector[d, # + 1] & /@ {ks}];

(* G1 clock & shift (Weyl-Heisenberg), d=3 *)
\[Omega] = Exp[2 Pi I/3];
Xd = RotateRight[IdentityMatrix[3]];
Zd = DiagonalMatrix[\[Omega]^Range[0, 2]];
p["G1-X", Xd];
p["G1-Z", Zd];
p["G1-weyl", Zd . Xd == \[Omega] Xd . Zd];

(* G2 qutrit SUM gate *)
sum = Sum[KroneckerProduct[KroneckerProduct[UnitVector[3, a + 1], UnitVector[3, a + 1]], MatrixPower[Xd, a]], {a, 0, 2}];
p["G2-sum|10>", sum . ketd[3, 1, 0]];

(* G3 spin-1 algebra *)
Jz = {{1, 0, 0}, {0, 0, 0}, {0, 0, -1}};
Jx = {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}}/Sqrt[2];
Jy = {{0, -I, 0}, {I, 0, -I}, {0, I, 0}}/Sqrt[2];
p["G3-Jz-eig", Eigenvalues[Jz]];
p["G3-comm", Commutator[Jx, Jy, Dot] == I Jz];
p["G3-casimir", Jx . Jx + Jy . Jy + Jz . Jz];

(* G4 qutrit Fourier / MUB *)
F = Table[\[Omega]^(j k), {j, 0, 2}, {k, 0, 2}]/Sqrt[3];
p["G4-MUB", Simplify[Abs[F[[1]]]^2]];
p["G4-unitary", Simplify[ConjugateTranspose[F] . F] == IdentityMatrix[3]];

(* G5 KCBS contextuality *)
ct2 = Cos[Pi/5]/(1 + Cos[Pi/5]); ct = Sqrt[ct2]; st = Sqrt[1 - ct2];
lvec[k_] := {ct, st Sin[4 Pi k/5], st Cos[4 Pi k/5]};
p["G5-ortho", FullSimplify[lvec[0] . lvec[1]]];
p["G5-S", FullSimplify[Sum[Abs[Conjugate[{1, 0, 0}] . lvec[k]]^2, {k, 0, 4}]]];
