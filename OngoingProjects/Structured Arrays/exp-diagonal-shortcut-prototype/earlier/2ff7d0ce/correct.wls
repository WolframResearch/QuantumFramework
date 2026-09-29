(* Correctness battery for the diagonal exponential. Run as
     wolframscript -file correct.wls base|proto
   Lines tagged SAME are printed as a hash of the full InputForm so that the base and
   prototype runs can be diffed; lines tagged ERR print a numeric error against an
   independent reference; lines tagged MSG print the messages a case raised. *)
which = Last[$ScriptCommandLine];
PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], which, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
out[tag_, label_, x_] := WriteString[$Output, tag, "\t", label, "\t", ToString[x, InputForm], "\n"];
SetAttributes[{same, msgs}, HoldRest];
same[label_, e_] := out["SAME", label, Hash[Normal[e]]];
msgs[label_, e_] := Module[{m = {}}, Internal`HandlerBlock[{"Message", (AppendTo[m, #] &)}, e]; out["MSG", label, Length[m]]];
mat[qo_] := Normal[qo["Matrix"]];
err[a_, b_] := Max[Abs[N[Normal[a] - Normal[b]]]];

zz[n_] := With[{J = RandomReal[{-1, 1}, {n, n}], z = 1 - 2 Tuples[{0, 1}, n]},
    Total[Flatten[Table[J[[i, j]] z[[All, i]] z[[All, j]], {i, n}, {j, i + 1, n}], 1]]];
SeedRandom[11]; h8 = zz[8];
ising8 = QuantumOperator[SparseArray[Band[{1, 1}] -> h8], Range[8]];

(* 1. numeric diagonal, three spellings, against dense MatrixExp; unitarity *)
u = Exp[-I 0.3 ising8];
out["ERR", "Ising8 Exp vs dense MatrixExp", err[mat[u], MatrixExp[-I 0.3 Normal[DiagonalMatrix[h8]]]]];
out["ERR", "Ising8 Exp vs entrywise", err[Diagonal[mat[u]], Exp[-I 0.3 h8]]];
out["ERR", "Ising8 unitarity", err[mat[u] . ConjugateTranspose[mat[u]], IdentityMatrix[256]]];
out["ERR", "Ising8 E^ vs Exp vs MatrixExp", {err[mat[E ^ (-I 0.3 ising8)], mat[u]], err[mat[MatrixExp[-I 0.3 ising8]], mat[u]]}];
out["ERR", "Ising8 2^qo vs MatrixExp[Log2]", err[mat[2 ^ (-I 0.3 ising8)], MatrixExp[Log[2] (-I 0.3) DiagonalMatrix[h8]]]];
out["HEAD", "Ising8 result matrix head", Head[u["Matrix"]]];
out["ERR", "Ising8 theta->0 limit is identity", err[mat[Exp[0. ising8]], IdentityMatrix[256]]];

(* 2. exact and symbolic diagonal *)
exact = QuantumOperator[DiagonalMatrix[{I Pi, Log[2], 0, -I Pi / 2}], {1, 2}];
same["exact diag Exp", Exp[exact]["Matrix"]];
out["VAL", "exact diag Exp", mat[Exp[exact]]];
sym = QuantumOperator[DiagonalMatrix[{a, b, c, d}], {1, 2}];
same["symbolic diag Exp", Exp[sym]["Matrix"]];
out["VAL", "symbolic diag Exp", mat[Exp[sym]]];
same["symbolic diag x^qo", (x ^ sym)["Matrix"]];
same["symbolic MatrixExp[qo]", MatrixExp[sym]["Matrix"]];
same["symbolic Exp label", Exp[sym]["Label"]];

(* 3. named operators built through the exponential *)
same["RZ[theta]", QuantumOperator["RZ"[\[Theta]]]["Matrix"]];
same["RZ[0.4]", QuantumOperator["RZ"[0.4]]["Matrix"]];
same["RZ[theta, 3]", QuantumOperator["RZ"[\[Theta], 3]]["Matrix"]];
same["R[theta, ZZ]", QuantumOperator["R"[\[Theta], "ZZ"]]["Matrix"]];
same["R[theta, ZZZ]", QuantumOperator["R"[\[Theta], "ZZZ"]]["Matrix"]];
same["R[0.7, ZIZ]", QuantumOperator["R"[0.7, "ZIZ"]]["Matrix"]];
same["RX[theta]", QuantumOperator["RX"[\[Theta]]]["Matrix"]];
same["R[theta, XX]", QuantumOperator["R"[\[Theta], "XX"]]["Matrix"]];
same["RZ label", QuantumOperator["RZ"[\[Theta]]]["Label"]];

(* 4. parametrized diagonal: substitutions, collisions *)
hp = QuantumOperator[SparseArray[Band[{1, 1}] -> -I \[Gamma] h8], Range[8], "Parameters" -> {\[Gamma]}];
up = Exp[hp];
out["HEAD", "param Exp amplitudes head", Head[up["State"]["State"]]];
out["ERR", "param Exp at 0.3 vs numeric Exp", err[mat[up[0.3]], mat[u]]];
out["ERR", "param Exp at 0 is identity", err[mat[up[0]], IdentityMatrix[256]]];
same["param Exp closed form", up["Matrix"]];
coll = QuantumOperator[DiagonalMatrix[{\[Alpha], \[Beta]}], "Parameters" -> {\[Alpha], \[Beta]}];
out["VAL", "collision Exp at a=b=1", mat[Exp[coll][1, 1]]];
out["VAL", "collision Sin at a=b=0", mat[Sin[coll][0, 0]]];
out["VAL", "collision Exp[Sin] at a=b=0", mat[Exp[Sin[coll]][0, 0]]];
out["VAL", "collision MatrixExp at a=b=2", mat[MatrixExp[coll][2, 2]]];
offd = QuantumOperator[{{\[Alpha], \[Beta]}, {\[Beta], -\[Alpha]}}, "Parameters" -> {\[Alpha], \[Beta]}];
out["VAL", "non-diagonal param Exp at 0,0", mat[Exp[offd][0, 0]]];
out["ERR", "non-diagonal param Exp at 1/3,1/5", err[mat[Exp[offd][1/3, 1/5]], MatrixExp[{{1/3, 1/5}, {1/5, -1/3}}]]];
vanish = QuantumOperator[{{\[Alpha], \[Beta] \[Alpha]}, {0, 1}}, "Parameters" -> {\[Alpha], \[Beta]}];
out["ERR", "off-diagonal vanishing only at some values", err[mat[Exp[vanish][0.5, 0.2]], MatrixExp[{{0.5, 0.1}, {0, 1}}]]];

(* 5. action on a state *)
SeedRandom[3]; v8 = Normalize[RandomComplex[{-1 - I, 1 + I}, 256]];
out["ERR", "MatrixExp[qo, qs] vs dense", err[MatrixExp[-I 0.3 ising8, QuantumState[v8]]["StateVector"], MatrixExp[-I 0.3 DiagonalMatrix[h8], v8]]];

(* 6. matrix-type: dephasing Liouvillian, n = 3 *)
liou[n_, h_] := QuantumOperator["Liouvillian"[
    QuantumOperator[SparseArray[Band[{1, 1}] -> h], Range[n]],
    Table[QuantumOperator[SparseArray[KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2], n], k -> {{1, 0}, {0, -1}}]], Range[n]], {k, n}],
    ConstantArray[0.2, n]]];
SeedRandom[5]; h3 = zz[3]; L3 = liou[3, h3];
out["ERR", "Liouvillian Exp vs dense MatrixExp", err[mat[Exp[0.7 L3]], MatrixExp[0.7 Normal[L3["Matrix"]]]]];
rho0 = QuantumState["RandomMixed"[3]];
rhot = Exp[0.7 L3][rho0];
out["ERR", "Liouvillian trace preserved", Abs[Tr[rhot["DensityMatrix"]] - 1]];
out["ERR", "Liouvillian populations conserved", err[Diagonal[rhot["DensityMatrix"]], Diagonal[rho0["DensityMatrix"]]]];
out["ERR", "Liouvillian coherence decay closed form", With[{x = 1, y = 6, dist = HammingDistance[IntegerDigits[0, 2, 3], IntegerDigits[5, 2, 3]]},
    Abs[rhot["DensityMatrix"][[x, y]] - rho0["DensityMatrix"][[x, y]] Exp[-I (h3[[x]] - h3[[y]]) 0.7 - 2 0.2 dist 0.7]]]];
out["ERR", "Liouvillian MatrixExp[qo, qs] vs Exp[qo][qs]", err[MatrixExp[0.7 L3, rho0]["DensityMatrix"], rhot["DensityMatrix"]]];

(* 7. non-diagonal controls: must be unchanged *)
tfi4 = QuantumOperator[Total[Table[KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2], 4], k -> {{0, 1}, {1, 0}}], {k, 4}]], Range[4]];
same["TFI4 Exp", Exp[-I 0.3 tfi4]["Matrix"]];
SeedRandom[9]; herm = With[{m = RandomComplex[{-1 - I, 1 + I}, {16, 16}]}, (m + ConjugateTranspose[m]) / 2];
same["random Hermitian Exp", Exp[-I 0.3 QuantumOperator[herm, {1, 2, 3, 4}]]["Matrix"]];
same["exact non-diagonal Exp", Exp[QuantumOperator[{{0, 1}, {1, 0}}]]["Matrix"]];
same["symbolic non-diagonal Exp", Exp[QuantumOperator[{{a, b}, {b, a}}]]["Matrix"]];
nearOff = QuantumOperator[{{1., 1.*^-12}, {1.*^-12, 2.}}];
same["coupling 1e-12 (not diagonal) Exp", Exp[-I nearOff]["Matrix"]];

(* 8. roundoff-level coupling (counts as diagonal): error against MatrixExp *)
tiny = {{1., 1.*^-15}, {1.*^-15, 2.}};
out["ERR", "coupling 1e-15 dropped vs MatrixExp", err[mat[Exp[-I QuantumOperator[tiny]]], MatrixExp[-I tiny]]];
bigSpread = {{1.*^6, 0, 0}, {0, 0., 1.*^-9}, {0, 1.*^-9, 1.*^-3}};
out["BOOL", "1e-9 coupling beside a 1e6 entry counts as diagonal", DiagonalMatrixQ[bigSpread]];
out["ERR", "1e-9 coupling beside 1e6, t = 1: shortcut vs MatrixExp", err[mat[Exp[-I QuantumOperator[bigSpread, 3]]], MatrixExp[-I bigSpread]]];
out["ERR", "1e-9 coupling beside 1e6, t = 1e6: shortcut vs exact eigen route", With[{t = 1.*^6},
    err[mat[Exp[-I t QuantumOperator[bigSpread, 3]]][[2 ;;, 2 ;;]], MatrixExp[-I t bigSpread[[2 ;;, 2 ;;]]]]]];

(* 9. underflow and edge cases *)
msgs["underflow Exp[-2000 qo]", out["VAL", "underflow Exp[-2000 qo]", mat[Exp[-2000. QuantumOperator[DiagonalMatrix[{1., 0.}]]]]]];
msgs["complex underflow", out["VAL", "complex underflow", mat[Exp[QuantumOperator[DiagonalMatrix[{-2000. + 3. I, 1. I}]]]]]];
msgs["Ising8 Exp messages", Exp[-I 0.3 ising8]];
out["VAL", "zero operator Exp", mat[Exp[QuantumOperator[SparseArray[{}, {4, 4}], {1, 2}]]]];
out["VAL", "non-square E^qo head", Head[E ^ QuantumOperator[{{1, 0}, {0, 1}, {0, 0}, {0, 0}}]]];
same["non-square E^qo", E ^ QuantumOperator[{{1, 0}, {0, 1}, {0, 0}, {0, 0}}]];
out["VAL", "arbitrary precision diag Exp precision", Precision[mat[Exp[QuantumOperator[DiagonalMatrix[N[{1, 2}, 30]]]]]]];

(* 10. resonant roundoff coupling beside a large entry, and the same through f[qo] *)
res = {{1.*^6, 0., 0.}, {0., 0., 1.*^-9}, {0., 1.*^-9, 0.}};
resBlock = MatrixExp[-I 1.*^6 {{0., 1.*^-9}, {1.*^-9, 0.}}];
out["ERR", "resonant Exp block error, t = 1e6", err[mat[Exp[-I 1.*^6 QuantumOperator[res, 3]]][[2 ;;, 2 ;;]], resBlock]];
out["ERR", "resonant Cos[qo] block error, t = 1e6", err[mat[Cos[1.*^6 QuantumOperator[res, 3]]][[2 ;;, 2 ;;]], MatrixFunction[Cos, 1.*^6 {{0., 1.*^-9}, {1.*^-9, 0.}}]]];
out["BOOL", "symbolic hidden zero strict", DiagonalMatrixQ[{{a, Sin[y]^2 + Cos[y]^2 - 1}, {0, b}}, Tolerance -> 0]];
same["symbolic hidden-zero Exp", Exp[QuantumOperator[{{a, Sin[y]^2 + Cos[y]^2 - 1}, {0, b}}]]["Matrix"]];
out["PACKED", "underflow result packed", Developer`PackedArrayQ[Normal[Exp[-2000. QuantumOperator[DiagonalMatrix[{1., 0.}]]]["Matrix"]]]];
