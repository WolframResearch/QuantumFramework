p[l_, e_] := WriteString[$Output, l <> " = " <> ToString[InputForm[e]] <> "\n"];
\[Sigma] = Table[PauliMatrix[j], {j, 3}];
densityMatrix[r_] := 1/2 (IdentityMatrix[2] + r . \[Sigma]);
vne[\[Rho]_] := -Total[# Log2[#] & /@ Select[Re@Eigenvalues[\[Rho]], # > 10^-12 &]];

(* L1 no-cloning *)
ov = Conjugate[{1, 0}] . ({1, 1}/Sqrt[2]);
p["L1", {ov, ov^2, ov == ov^2}];

(* L3 Holevo *)
rhoEns = 1/2 densityMatrix[{0, 0, 1}] + 1/2 densityMatrix[{1, 0, 0}];
p["L3-holevo", N@vne[rhoEns]];

(* L4 GHZ / Mermin *)
ghz = {1, 0, 0, 0, 0, 0, 0, 1}/Sqrt[2];
X = PauliMatrix[1]; Y = PauliMatrix[2];
ev[ops_] := Re[Conjugate[ghz] . KroneckerProduct @@ ops . ghz];
p["L4", {ev[{X, X, X}], ev[{X, Y, Y}], ev[{Y, X, Y}], ev[{Y, Y, X}]}];

(* Part M: truncated Fock space *)
n = 12;
a = DiagonalMatrix[Sqrt@Range[n - 1], 1];
fock[k_] := UnitVector[n, k + 1];

(* M1 ladder *)
p["M1-a|2>", Chop[a . fock[2] - Sqrt[2] fock[1]]];
p["M1-N|3>", Chop[(ConjugateTranspose[a] . a) . fock[3] - 3 fock[3]]];

(* M2 coherent state *)
al = 1.2;
disp[\[Alpha]_] := MatrixExp[\[Alpha] ConjugateTranspose[a] - Conjugate[\[Alpha]] a];
coh = disp[al] . fock[0];
p["M2-eigen", Chop[(a . coh - al coh)][[1 ;; 6]]];
p["M2-meanN", Chop[Re[Conjugate[coh] . (ConjugateTranspose[a] . a) . coh]]];

(* M3 Poisson *)
p["M3-probs", Chop[Abs[coh]^2][[1 ;; 6]]];
p["M3-poisson", Table[N@PDF[PoissonDistribution[al^2], k], {k, 0, 5}]];

(* M4 squeezing *)
squeeze[\[Xi]_] := MatrixExp[1/2 (Conjugate[\[Xi]] a . a - \[Xi] ConjugateTranspose[a] . ConjugateTranspose[a])];
quadX = (a + ConjugateTranspose[a])/Sqrt[2];
var[\[Psi]_, op_] := Re[Conjugate[\[Psi]] . op . op . \[Psi] - (Conjugate[\[Psi]] . op . \[Psi])^2];
p["M4-vacuum", var[fock[0], quadX]];
p["M4-squeezed", var[squeeze[1.0] . fock[0], quadX]];

(* M5 g2(0) *)
g2rho[\[Rho]_] := Re[Tr[\[Rho] . ConjugateTranspose[a] . ConjugateTranspose[a] . a . a]]/Re[Tr[\[Rho] . ConjugateTranspose[a] . a]]^2;
g2vec[\[Psi]_] := g2rho[KroneckerProduct[\[Psi], Conjugate@\[Psi]]];
thermal = With[{x = 0.5}, (1 - x) DiagonalMatrix[x^Range[0, n - 1]]];
thermal = thermal/Tr[thermal];
p["M5-g2", {g2rho[thermal], g2vec[coh], g2vec[fock[1]]}];
