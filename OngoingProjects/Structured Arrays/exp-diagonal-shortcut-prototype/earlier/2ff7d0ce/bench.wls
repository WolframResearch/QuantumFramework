(* Timing harness: wolframscript -file bench.wls base|proto group
   Every time is the best of k runs after ClearSystemCache[] (k = 3 below 5 s, one run
   above), wrapped in TimeConstrained; "TIMEOUT" marks a case that did not finish. *)
{which, group} = Take[$ScriptCommandLine, -2];
PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], which, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
out[label_, x_] := WriteString[$Output, which, "\t", group, "\t", label, "\t", ToString[x, InputForm], "\n"];
SetAttributes[{tm, once}, HoldFirst];
once[e_, limit_] := (ClearSystemCache[]; TimeConstrained[Round[First[AbsoluteTiming[e;]], 0.0001], limit, "TIMEOUT"]);
tm[e_, limit_ : 300] := With[{t1 = once[e, limit]}, If[t1 === "TIMEOUT" || t1 > 5, t1, Min[t1, once[e, limit], once[e, limit]]]];

zz[n_] := With[{J = RandomReal[{-1, 1}, {n, n}], z = 1 - 2 Tuples[{0, 1}, n]},
    Total[Flatten[Table[J[[i, j]] z[[All, i]] z[[All, j]], {i, n}, {j, i + 1, n}], 1]]];
zOn[n_, k_] := QuantumOperator[SparseArray[KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2, SparseArray], n], k -> SparseArray[{{1, 1} -> 1, {2, 2} -> -1}]]], Range[n]];
xSum[n_] := QuantumOperator[Total[Table[KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2, SparseArray], n], k -> SparseArray[{{1, 2} -> 1., {2, 1} -> 1.}]], {k, n}]], Range[n]];

Switch[group,

"ising",
Scan[
    With[{n = #}, SeedRandom[7]; Module[{h = zz[n], qo, v, r},
        qo = QuantumOperator[SparseArray[Band[{1, 1}] -> h], Range[n]];
        out["Exp[-I 0.3 qo] n=" <> ToString[n], tm[r = Exp[-I 0.3 qo], 400]];
        If[r =!= Null && Head[r] === QuantumOperator,
            out["  error vs entrywise n=" <> ToString[n], Max[Abs[Normal[Diagonal[r["Matrix"]]] - Exp[-I 0.3 h]]]]];
        SeedRandom[1]; v = QuantumState[Normalize[RandomComplex[{-1 - I, 1 + I}, 2^n]]];
        out["MatrixExp[-I 0.3 qo, qs] n=" <> ToString[n], tm[MatrixExp[-I 0.3 qo, v], 400]];
    ]] &,
    {12, 13, 14, 15, 16}
],

"liouvillian",
Scan[
    With[{n = #}, SeedRandom[7]; Module[{h = zz[n], L, H, r, m},
        H = QuantumOperator[SparseArray[Band[{1, 1}] -> h], Range[n]];
        out["build Liouvillian n=" <> ToString[n], tm[L = QuantumOperator["Liouvillian"[H, Table[zOn[n, k], {k, n}], ConstantArray[0.2, n]]], 600]];
        out["  Head stored state / Head Matrix n=" <> ToString[n], {Head[L["State"]["State"]], tm[m = (0.7 L)["Sort"]["Matrix"]], Head[m]}];
        out["Exp[0.7 L] n=" <> ToString[n], tm[r = Exp[0.7 L], 400]];
        If[Head[r] === QuantumOperator,
            out["  error vs entrywise n=" <> ToString[n], Max[Abs[Normal[Diagonal[r["Matrix"]]] - Exp[Normal[Diagonal[m]]]]]]];
    ]] &,
    If[which === "proto2", {4, 5, 6, 7}, {4, 5, 6}]
],

"qaoa",
Scan[
    With[{n = #}, SeedRandom[7]; Module[{h = zz[n], hp, U, pts25 = Subdivide[0.05, 1.25, 24], pts100 = Subdivide[0.05, 1.25, 99], one, byte, limit},
        hp = QuantumOperator[SparseArray[Band[{1, 1}] -> -I \[Gamma] h], Range[n], "Parameters" -> {\[Gamma]}];
        out["param build Exp n=" <> ToString[n], tm[U = Exp[hp], 600]];
        If[Head[U] === QuantumOperator,
            out["  ByteCount n=" <> ToString[n], ByteCount[U]];
            out["  one point n=" <> ToString[n], one = tm[U[0.3]["Matrix"], 600]];
            out["  error at 0.3 n=" <> ToString[n], Max[Abs[Normal[Diagonal[U[0.3]["Matrix"]]] - Exp[-I 0.3 h]]]];
            If[NumberQ[one] && one < 2,
                out["  25 points n=" <> ToString[n], tm[U[#]["Matrix"] & /@ pts25, 600]];
                out["  100 points n=" <> ToString[n], tm[U[#]["Matrix"] & /@ pts100, 600]]
            ]
        ]
    ]] &,
    If[which === "base", {8, 10, 12}, {8, 10, 12, 14, 16}]
],

"controls",
Module[{tfi, herm8, herm10, near, sym, r},
    SeedRandom[3];
    tfi = xSum[10];
    out["TFI n=10 Exp", tm[r = Exp[-I 0.3 tfi]]];
    herm8 = QuantumOperator[With[{m = RandomComplex[{-1 - I, 1 + I}, {256, 256}]}, (m + ConjugateTranspose[m]) / 2], Range[8]];
    out["random Hermitian n=8 Exp", tm[Exp[-I 0.3 herm8]]];
    herm10 = QuantumOperator[With[{m = RandomComplex[{-1 - I, 1 + I}, {1024, 1024}]}, (m + ConjugateTranspose[m]) / 2], Range[10]];
    out["random Hermitian n=10 Exp", tm[Exp[-I 0.3 herm10]]];
    near = QuantumOperator[SparseArray[Band[{1, 1}] -> zz[10]] + SparseArray[{{1, 1024} -> 1.*^-12, {1024, 1} -> 1.*^-12}, {1024, 1024}], Range[10]];
    out["diagonal + 1e-12 corner coupling n=10 Exp", tm[Exp[-I 0.3 near]]];
    sym = QuantumOperator[{{a, b, 0, 0}, {b, a, 0, 0}, {0, 0, c, d}, {0, 0, d, c}}, {1, 2}];
    out["symbolic block 4x4 Exp", tm[Exp[sym]]];
    out["RX[theta] build", tm[QuantumOperator["RX"[\[Theta]]]]];
    out["R[theta, XX] build", tm[QuantumOperator["R"[\[Theta], "XX"]]]];
    out["RZ[theta] build", tm[QuantumOperator["RZ"[\[Theta]]]]];
    out["R[theta, ZZ] build", tm[QuantumOperator["R"[\[Theta], "ZZ"]]]];
    out["R[0.3, ZZZZ] build", tm[QuantumOperator["R"[0.3, "ZZZZ"]]]];
    out["Trotter 4-qubit ZZ+X order 2 x 5 build", tm[QuantumCircuitOperator["Trotterization"[{QuantumOperator["ZZ", {1, 2}], QuantumOperator["ZZ", {2, 3}], QuantumOperator["ZZ", {3, 4}], xSum[4]}, 2, 5]]]];
],

"detect",
Module[{m, sa, tfi, sym},
    SeedRandom[3];
    m = RandomComplex[{-1 - I, 1 + I}, {1024, 1024}];
    sa = SparseArray[Band[{1, 1}] -> RandomReal[1, 2^16]];
    tfi = xSum[12]["Matrix"];
    sym = Table[Sum[Subscript[x, i, j, k], {k, 3}], {i, 64}, {j, 64}];
    out["strict test, dense random 1024 (non-diagonal)", tm[DiagonalMatrixQ[m, Tolerance -> 0]]];
    out["strict test, SparseArray diagonal 2^16", tm[DiagonalMatrixQ[sa, Tolerance -> 0]]];
    out["strict test, sparse TFI n=12 (non-diagonal)", tm[DiagonalMatrixQ[tfi, Tolerance -> 0]]];
    out["strict test, symbolic dense 64 (non-diagonal)", tm[DiagonalMatrixQ[sym, Tolerance -> 0]]];
    out["strict test, dense diagonal 4096", With[{d = Normal[SparseArray[Band[{1, 1}] -> RandomComplex[1, 4096]]]}, tm[DiagonalMatrixQ[d, Tolerance -> 0]]]];
]
]
