PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
b[e_, ___] := e;
SetAttributes[b, HoldFirst];
P[l_, x_] := Print["@@ ", StringPadRight[l, 52], " ", ToString[x, InputForm]];

(* 1. Two-qubit tomography from simulated shots: linear inversion, maximum likelihood, Bayesian *)
SeedRandom[7];
target = QuantumState[N @ {Cos[Pi/8], 0, 0, Sin[Pi/8]}];   (* partially entangled *)
mops = QuantumMeasurementOperator /@ (StringJoin /@ Tuples[{"X", "Y", "Z"}, 2]);
P["9 two-qubit Pauli settings build", Head /@ Take[mops, 2]];
Do[Module[{t1, data, t2, est},
   {t1, data} = AbsoluteTiming[b[QuantumMeasurementSimulation[target, mops, shots], 120]];
   P["simulate 9 x " <> ToString[shots] <> " shots (s)", Round[t1, 0.01]];
   {t2, est} = AbsoluteTiming[b[QuantumStateEstimate[data], 300]];
   P["  QuantumStateEstimate (s), head", {Round[t2, 0.01], Head[est]}];
   If[! StringQ[est],
     P["  properties", Quiet @ Check[est["Properties"], "?"]];
     P["  fidelity of MLE to target", Quiet @ Check[Chop @ QuantumDistance[est["MaximumLikelihoodState"], target, "Fidelity"], "?"]]]],
  {shots, {1000, 10000}}];

Print[];
(* 2. Variational ground state, 6-qubit transverse-field Ising chain, hardware-efficient ansatz:
   the model is symbolic, the optimum is found numerically, compared with exact diagonalization *)
n = 6; hx = 1.0;
Hising = Total[Table[-QuantumOperator["Z", {i}] @ QuantumOperator["Z", {i + 1}], {i, n - 1}]] +
   Total[Table[-hx QuantumOperator["X", {i}], {i, n}]];
e0 = Min @ Eigenvalues[N @ Normal @ Hising["Matrix"]];
P["exact ground energy (6 sites, open chain)", e0];
pars = Array[\[FormalTheta], 2 n];
ansatz = QuantumCircuitOperator @ Join[
   Table["RY"[\[FormalTheta][i]] -> i, {i, n}],
   Table["CNOT" -> {i, i + 1}, {i, n - 1}],
   Table["RY"[\[FormalTheta][n + i]] -> i, {i, n}]];
energy[th_ ? (VectorQ[#, NumericQ] &)] := Re @ Module[{psi = (ansatz /. Thread[pars -> th])[]},
   (psi["Dagger"] @ Hising @ psi)["Scalar"]];
{tE, eval} = AbsoluteTiming[b[energy[ConstantArray[0.1, 2 n]], 60]];
P["one energy evaluation (s)", Round[tE, 0.001]];
SeedRandom[3];
{tF, fm} = AbsoluteTiming[b[Quiet @ FindMinimum[energy[pars], Transpose[{pars, RandomReal[{-0.3, 0.3}, 2 n]}]], 600]];
P["FindMinimum VQE (s), energy", {Round[tF, 0.1], If[ListQ[fm], First[fm], fm]}];
P["  relative error vs exact", If[ListQ[fm], Abs[(First[fm] - e0)/e0], "n/a"]];

Print[];
(* 3. Fubini-Study metric of the ansatz state, the input to quantum natural gradient *)
{tG, g} = AbsoluteTiming[b[Quiet @ FubiniStudyMetricTensor[ansatz[], "Matrix", "Parameters" -> pars], 300]];
P["FubiniStudyMetricTensor 12 parameters (s), dims", {Round[tG, 0.01], Quiet@Check[Dimensions[g], Head[g]]}];
