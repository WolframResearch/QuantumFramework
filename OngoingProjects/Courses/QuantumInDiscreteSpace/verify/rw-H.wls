p[l_, e_] := WriteString[$Output, l <> " = " <> ToString[InputForm[e]] <> "\n"];
\[Sigma] = Table[PauliMatrix[j], {j, 3}];
densityMatrix[r_] := 1/2 (IdentityMatrix[2] + r . \[Sigma]);
partialT[\[Rho]_] := ArrayReshape[Transpose[ArrayReshape[\[Rho], {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}];

(* H1 Kraus channel: amplitude damping on |1> *)
g = 3/10;
K0 = {{1, 0}, {0, Sqrt[1 - g]}}; K1 = {{0, Sqrt[g]}, {0, 0}};
channel[\[Rho]_, kr_] := Sum[k . \[Rho] . ConjugateTranspose[k], {k, kr}];
p["H1-AD|1>", channel[{{0, 0}, {0, 1}}, {K0, K1}]];

(* H2 Choi + CP (Kraus-based) *)
choi = Sum[KroneckerProduct[Flatten[k], Conjugate@Flatten[k]], {k, {K0, K1}}];
p["H2-choi", choi];
p["H2-CP", PositiveSemidefiniteMatrixQ[N@choi]];

(* H3 named noise as Kraus *)
bitflip[p_] := {Sqrt[1 - p] IdentityMatrix[2], Sqrt[p] PauliMatrix[1]};
depol[p_] := Prepend[Sqrt[p/3] # & /@ Table[PauliMatrix[j], {j, 3}], Sqrt[1 - p] IdentityMatrix[2]];
p["H3-bitflip|0>", channel[{{1, 0}, {0, 0}}, bitflip[1/4]]];
p["H3-depol|+>", channel[densityMatrix[{1, 0, 0}], depol[1/2]] // Simplify];

(* H4 trace preservation + Stinespring dilation *)
p["H4-TP", Simplify[ConjugateTranspose[K0] . K0 + ConjugateTranspose[K1] . K1]];
stinespring = ArrayFlatten[{{K0}, {K1}}];    (* isometry: system -> system (x) env *)
p["H4-isometry", Simplify[ConjugateTranspose[stinespring] . stinespring == IdentityMatrix[2]]];

(* H5/H6 Lindblad master equation via DSolve (T1 decay -> ground state) *)
L = {{0, 1}, {0, 0}};
diss[\[Rho]_] := L . \[Rho] . ConjugateTranspose[L] - 1/2 (ConjugateTranspose[L] . L . \[Rho] + \[Rho] . ConjugateTranspose[L] . L);
\[Rho]t = DSolveValue[{\[Rho]'[t] == diss[\[Rho][t]], \[Rho][0] == {{0, 0}, {0, 1}}}, \[Rho][t] \[Element] Matrices[{2, 2}], t];
p["H5-rho(t)", \[Rho]t // Simplify];
p["H5-pop(1)", (\[Rho]t /. t -> 1)[[2, 2]]];
p["H6-steady", \[Rho]t /. t -> Infinity];

(* H7 entanglement sudden death (both qubits amplitude-damped) *)
conc[v_, gg_] := Module[{k0, k1, rho, kr},
   k0 = {{1, 0}, {0, Sqrt[1 - gg]}}; k1 = {{0, Sqrt[gg]}, {0, 0}};
   rho = KroneckerProduct[v, Conjugate@v];
   kr = Flatten[Table[KroneckerProduct[a, b], {a, {k0, k1}}, {b, {k0, k1}}], 1];
   -Total@Select[Re@Eigenvalues@partialT[Sum[m . rho . ConjugateTranspose[m], {m, kr}]], # < 0 &]];
p["H7-phiplus", Table[conc[{1, 0, 0, 1}/Sqrt[2], gg], {gg, {0, 0.3, 0.5, 0.9}}] // N];
p["H7-weighted", Table[conc[{Sqrt[1/5], 0, 0, Sqrt[4/5]}, gg], {gg, {0, 0.3, 0.5, 0.6}}] // Chop];
