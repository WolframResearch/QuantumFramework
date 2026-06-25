p[l_, e_] := WriteString[$Output, l <> " = " <> ToString[InputForm[e]] <> "\n"];
\[Sigma] = Table[PauliMatrix[j], {j, 3}];
densityMatrix[r_] := 1/2 (IdentityMatrix[2] + r . \[Sigma]);
reduceA[\[Rho]_] := TensorContract[ArrayReshape[\[Rho], {2, 2, 2, 2}], {{2, 4}}];
partialT[\[Rho]_] := ArrayReshape[Transpose[ArrayReshape[\[Rho], {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}];
vne[\[Rho]_] := -Total[# Log2[#] & /@ Select[Re@Eigenvalues[\[Rho]], # > 10^-12 &]];

(* E1 ensemble + purity *)
\[Rho] = 1/2 densityMatrix[{0, 0, 1}] + 1/2 densityMatrix[{1, 0, 0}];
p["E1-rho", \[Rho]];
p["E1-purity", Tr[\[Rho] . \[Rho]]];
p["E1-pure", Tr[densityMatrix[{0, 0, 1}] . densityMatrix[{0, 0, 1}]]];
p["E1-maxmix", Tr[(IdentityMatrix[2]/2) . (IdentityMatrix[2]/2)]];

(* E2 von Neumann entropy *)
p["E2-pure", vne[densityMatrix[{0, 0, 1}]]];
p["E2-maxmix", vne[IdentityMatrix[2]/2]];
p["E2-bell", vne[reduceA[KroneckerProduct[{1, 0, 0, 1}/Sqrt[2], {1, 0, 0, 1}/Sqrt[2]]]]];

(* E3 Werner threshold (standard convention) *)
singlet = {0, 1, -1, 0}/Sqrt[2];
werner[q_] := q KroneckerProduct[singlet, singlet] + (1 - q) IdentityMatrix[4]/4;
neg[\[Rho]_] := Total@Abs@Select[Eigenvalues[partialT[\[Rho]]], # < 0 &];
p["E3-neg", FullSimplify[neg[werner[q]], 0 <= q <= 1]];
p["E3-threshold", Reduce[neg[werner[q]] > 0 && 0 <= q <= 1, q]];

(* E4 fidelity + trace distance *)
\[Psi]0 = {1, 0}; \[Psi]p = {1, 1}/Sqrt[2];
overlap = Conjugate[\[Psi]0] . \[Psi]p;
p["E4-rootF-F", {Abs[overlap], Abs[overlap]^2}];
\[Rho]0 = TensorProduct[\[Psi]0, \[Psi]0]; \[Rho]p = TensorProduct[\[Psi]p, Conjugate@\[Psi]p];
p["E4-tracedist", 1/2 Total@Abs@Eigenvalues[\[Rho]0 - \[Rho]p]];

(* E5 purification *)
purify[\[Rho]_] := With[{e = Eigensystem[\[Rho]]}, Total@MapThread[Sqrt[#1] Flatten@KroneckerProduct[#2, #2] &, e]];
p["E5-purified", purify[IdentityMatrix[2]/2]];
p["E5-reduced", reduceA[KroneckerProduct[purify[IdentityMatrix[2]/2], purify[IdentityMatrix[2]/2]]]];
