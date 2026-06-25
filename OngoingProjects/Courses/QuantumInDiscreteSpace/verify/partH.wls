PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
SetAttributes[rec, HoldRest];
rec[id_, e_] := Module[{v = e}, WriteString[$Output, id <> " = " <> ToString[InputForm[v]] <> "\n"]; v];
SetAttributes[safe, HoldFirst]; safe[e_] := Check[e, "ERR"];

(* H1 Kraus channel: amplitude damping on |1> *)
g = 3/10;
K0 = {{1, 0}, {0, Sqrt[1 - g]}}; K1 = {{0, Sqrt[g]}, {0, 0}};
ad = QuantumChannel[{K0, K1}];
rec["H1-AD|1>", ad[QuantumState["1"]]["DensityMatrix"] // Normal // safe];
rec["H1-AD|+>", ad[QuantumState["+"]]["DensityMatrix"] // Normal // safe];

(* H2 Choi matrix from Kraus + complete positivity (Choi PSD) *)
vec[K_] := Flatten[K];
choi = Sum[KroneckerProduct[vec[K], Conjugate[vec[K]]], {K, {K0, K1}}];
rec["H2-choi", choi // Normal];
rec["H2-CP-PSD", PositiveSemidefiniteMatrixQ[N[choi]]];

(* H3 named noise channels *)
rec["H3-bitflip|0>", QuantumChannel["BitFlip"[1/4]][QuantumState["0"]]["DensityMatrix"] // Normal // safe];
rec["H3-depol|+>", QuantumChannel["Depolarizing"[1/2]][QuantumState["+"]]["DensityMatrix"] // Normal // safe];

(* H4 trace preservation + Stinespring dilation *)
rec["H4-TP-sumKdK", Simplify[ConjugateTranspose[K0] . K0 + ConjugateTranspose[K1] . K1]];
rec["H4-stinespring-dims", ad["QuantumOperator"]["MatrixRepresentation"] // Dimensions // safe];

(* H5/H6 Lindblad master equation (rule form {L}->{g}): T1 decay of |1> -> steady |0> *)
Lminus = QuantumOperator[{{0, 1}, {0, 0}}];
H0 = QuantumOperator[{{0, 0}, {0, 0}}];
sol = QuantumEvolve[H0, {Lminus} -> {1}, QuantumState["1"], {t, 0, 5}] // safe;
rec["H5-pop-t0", Re[(sol[0]["DensityMatrix"] // Normal)[[2, 2]]]];
rec["H5-pop-t1", Re[(sol[1]["DensityMatrix"] // Normal)[[2, 2]]]];
rec["H6-steady-t5", Chop[Re[(sol[5]["DensityMatrix"] // Normal)]]];

(* H7 entanglement sudden death: amplitude damping on BOTH qubits *)
conc[vec_, gg_] := Module[{k0, k1, rho, kr, rho2},
   k0 = {{1, 0}, {0, Sqrt[1 - gg]}}; k1 = {{0, Sqrt[gg]}, {0, 0}};
   rho = QuantumState[vec]["DensityMatrix"] // Normal;
   kr = Flatten[Table[KroneckerProduct[a, b], {a, {k0, k1}}, {b, {k0, k1}}], 1];
   rho2 = Sum[m . rho . ConjugateTranspose[m], {m, kr}];
   N@QuantumEntanglementMonotone[QuantumState[rho2], {1}, "Concurrence"]];
phiplus = {1, 0, 0, 1}/Sqrt[2];
weighted = {Sqrt[1/5], 0, 0, Sqrt[4/5]};
rec["H7-PhiPlus-g{0,.3,.5,.9}", Chop@Table[conc[phiplus, gg], {gg, {0, 0.3, 0.5, 0.9}}]];
rec["H7-weighted-g{0,.3,.5,.6}", Chop@Table[conc[weighted, gg], {gg, {0, 0.3, 0.5, 0.6}}]];
