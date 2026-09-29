(* Follow-ons, base paclet only: powers, spectrum and QuantumEvolve of a diagonal operator,
   each against the entrywise closed form. *)
PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "base", "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
out[label_, x_] := WriteString[$Output, label, "\t", ToString[x, InputForm], "\n"];
SetAttributes[{tm, once}, HoldFirst];
once[e_, limit_] := (ClearSystemCache[]; TimeConstrained[Round[First[AbsoluteTiming[e;]], 0.0001], limit, "TIMEOUT"]);
tm[e_, limit_ : 300] := With[{t1 = once[e, limit]}, If[t1 === "TIMEOUT" || t1 > 5, t1, Min[t1, once[e, limit], once[e, limit]]]];
zz[n_] := With[{J = RandomReal[{-1, 1}, {n, n}], z = 1 - 2 Tuples[{0, 1}, n]},
    Total[Flatten[Table[J[[i, j]] z[[All, i]] z[[All, j]], {i, n}, {j, i + 1, n}], 1]]];

(* MatrixPower qo^x on a diagonal: principal branch entry by entry *)
Scan[With[{n = #}, SeedRandom[7]; Module[{h = zz[n], qo, r},
    qo = QuantumOperator[SparseArray[Band[{1, 1}] -> h], Range[n]];
    out["qo^0.5 n=" <> ToString[n], tm[r = qo ^ 0.5]];
    out["  qo^0.5 error vs entrywise principal branch n=" <> ToString[n], Max[Abs[Normal[Diagonal[r["Matrix"]]] - h ^ 0.5]]];
    out["  entrywise time n=" <> ToString[n], tm[SparseArray[Band[{1, 1}] -> h ^ 0.5]]];
    out["qo^3 n=" <> ToString[n], tm[qo ^ 3]];
]] &, {10, 12}];

(* Eigenvalues / Eigensystem *)
Scan[With[{n = #}, SeedRandom[7]; Module[{h = zz[n], qo, ev},
    qo = QuantumOperator[SparseArray[Band[{1, 1}] -> h], Range[n]];
    out["Eigenvalues n=" <> ToString[n], tm[ev = qo["Eigenvalues"], 600]];
    out["  order: equals h sorted by decreasing Abs n=" <> ToString[n], If[ListQ[ev], ev === h[[Reverse[OrderingBy[h, Abs]]]], ev]];
    out["  order: equals h in diagonal order n=" <> ToString[n], If[ListQ[ev], ev === h, ev]];
    out["  Eigenvalues[structured] n=" <> ToString[n], tm[Eigenvalues[DiagonalMatrix[h, TargetStructure -> "Structured"]]]];
    out["  Eigenvalues[structured] sorted equal n=" <> ToString[n], If[ListQ[ev], Sort[Eigenvalues[DiagonalMatrix[h, TargetStructure -> "Structured"]]] == Sort[ev], ev]];
]] &, {8, 10}];
out["Eigensystem of structured DiagonalMatrix (WL 15.0)", Eigensystem[DiagonalMatrix[{1., 3., 4.}, TargetStructure -> "Structured"]]];
out["Eigensystem of SparseArray diagonal", Normal /@ Eigensystem[SparseArray[Band[{1, 1}] -> {1., 3., 4.}]]];
out["Eigensystem of exact degenerate diagonal", Eigensystem[DiagonalMatrix[{2, 1, 2}]]];

(* QuantumEvolve with a time-independent diagonal H *)
Scan[With[{n = #}, SeedRandom[7]; Module[{h = zz[n], H, psi0, psi, v},
    H = QuantumOperator[SparseArray[Band[{1, 1}] -> h], Range[n]];
    SeedRandom[2]; v = Normalize[RandomComplex[{-1 - I, 1 + I}, 2^n]]; psi0 = QuantumState[v];
    out["QuantumEvolve numeric {t,0,1} n=" <> ToString[n], tm[psi = QuantumEvolve[H, psi0, {\[FormalT], 0, 1}], 600]];
    out["  error at t=1 vs entrywise closed form n=" <> ToString[n], If[Head[psi] === QuantumState, Max[Abs[Normal[psi[1]["StateVector"]] - Exp[-I h] v]], psi]];
    out["  closed form time n=" <> ToString[n], tm[Exp[-I h] v]];
]] &, {6, 8, 10}];
Module[{H, psi},
    H = QuantumOperator[SparseArray[Band[{1, 1}] -> {1, -1, 2, 0}], {1, 2}];
    out["QuantumEvolve symbolic t, 2 qubits", tm[psi = QuantumEvolve[H, QuantumState["++"], \[FormalT]], 300]];
    out["  symbolic result", If[Head[psi] === QuantumState, Normal[psi["StateVector"]], psi]];
];
