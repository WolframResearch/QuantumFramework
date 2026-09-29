out[label_, x_] := WriteString[$Output, label, " => ", ToString[x, InputForm], "\n"];
SetAttributes[tm, HoldFirst];
tm[e_] := Min[Table[ClearSystemCache[]; First[AbsoluteTiming[e;]], 3]];
(* Tolerance semantics: absolute or relative? *)
out["Tol 1e-9 abs test, off 1e-10, diag 1,2", DiagonalMatrixQ[{{1., 1.*^-10}, {0., 2.}}, Tolerance -> 1.*^-9]];
out["Tol 1e-9, off 1e-10, diag 1e6", DiagonalMatrixQ[{{1.*^6, 1.*^-10}, {0., 2.*^6}}, Tolerance -> 1.*^-9]];
out["Tol 1e-9, off 1e-4, diag 1e6", DiagonalMatrixQ[{{1.*^6, 1.*^-4}, {0., 2.*^6}}, Tolerance -> 1.*^-9]];
out["Tol 1e-9, off 1e-10, diag 1e-6", DiagonalMatrixQ[{{1.*^-6, 1.*^-10}, {0., 2.*^-6}}, Tolerance -> 1.*^-9]];
out["Tol 0 on 0.+0.I entries", DiagonalMatrixQ[{{1. + 0. I, 0. + 0. I}, {0. + 0. I, 2. I}}, Tolerance -> 0]];
out["Tol 0 on -0.", DiagonalMatrixQ[{{1., -0.}, {0., 2.}}, Tolerance -> 0]];
(* resonant coupling beside a large entry *)
eps = 1.*^-9; t = 1.*^6;
M = {{1.*^6, 0., 0.}, {0., 0., eps}, {0., eps, 0.}};
exact = MatrixExp[-I t {{0., eps}, {eps, 0.}}];  (* 2x2 block alone: rotation by eps t = 1e-3 *)
out["resonant: default DMQ", DiagonalMatrixQ[M]];
out["resonant: exact block mixing |U23|", Abs[exact[[1, 2]]]];
out["resonant: MatrixExp full |U23| and block error", With[{u = MatrixExp[-I t M]}, {Abs[u[[2, 3]]], Max[Abs[u[[2 ;;, 2 ;;]] - exact]]}]];
out["resonant: entrywise block error", Max[Abs[DiagonalMatrix[Exp[-I t Diagonal[M]]][[2 ;;, 2 ;;]] - exact]]];
out["relative size eps / max|M|", eps / 1.*^6];
(* strict detection on real QF matrices *)
PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "base", "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
zz[n_] := With[{J = RandomReal[{-1, 1}, {n, n}], z = 1 - 2 Tuples[{0, 1}, n]},
    Total[Flatten[Table[J[[i, j]] z[[All, i]] z[[All, j]], {i, n}, {j, i + 1, n}], 1]]];
SeedRandom[7];
liou[n_, h_] := QuantumOperator["Liouvillian"[
    QuantumOperator[SparseArray[Band[{1, 1}] -> h], Range[n]],
    Table[QuantumOperator[SparseArray[KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2], n], k -> {{1, 0}, {0, -1}}]], Range[n]], {k, n}],
    ConstantArray[0.2, n]]];
L5 = (0.7 liou[5, zz[5]])["Sort"]["Matrix"];
out["Liouvillian n=5 strict DMQ", {DiagonalMatrixQ[L5, Tolerance -> 0], tm[DiagonalMatrixQ[L5, Tolerance -> 0]], tm[DiagonalMatrixQ[L5]]}];
I12 = (-I 0.3 QuantumOperator[SparseArray[Band[{1, 1}] -> zz[12]], Range[12]])["Sort"]["Matrix"];
out["Ising n=12 strict DMQ", {DiagonalMatrixQ[I12, Tolerance -> 0], tm[DiagonalMatrixQ[I12, Tolerance -> 0]], tm[DiagonalMatrixQ[I12]]}];
rz = QuantumOperator["PauliZ"[3]]["Sort"]["Matrix"];
out["PauliZ[3] clock matrix strict/default", {Normal[rz], DiagonalMatrixQ[rz, Tolerance -> 0], DiagonalMatrixQ[rz]}];
zzop = QuantumCircuitOperator[{QuantumOperator["Z", {1}], QuantumOperator["Z", {2}]}]["QuantumOperator", Method -> "Schrodinger"]["Sort"]["Matrix"];
out["ZZ via circuit strict", {DiagonalMatrixQ[zzop, Tolerance -> 0], Head[zzop]}];
dense4096 = RandomComplex[{-1 - I, 1 + I}, {4096, 4096}];
out["strict vs default cost, dense random 4096", {tm[DiagonalMatrixQ[dense4096, Tolerance -> 0]], tm[DiagonalMatrixQ[dense4096]]}];
dd = DiagonalMatrix[RandomReal[1, 4096]];
out["strict vs default cost, dense diag 4096", {tm[DiagonalMatrixQ[dd, Tolerance -> 0]], tm[DiagonalMatrixQ[dd]]}];
