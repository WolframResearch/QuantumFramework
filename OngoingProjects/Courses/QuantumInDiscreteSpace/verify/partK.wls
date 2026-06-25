PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`QuantumOptimization`"];
SetAttributes[rec, HoldRest];
rec[id_, e_] := Module[{v = e}, WriteString[$Output, id <> " = " <> ToString[InputForm[v]] <> "\n"]; v];
SetAttributes[safe, HoldFirst]; safe[e_] := Check[e, "ERR"];

(* K1 Deutsch-Jozsa: constant vs balanced in one query (Z on q1 = balanced f(x)=x1) *)
djLayers[oracle_] := QuantumCircuitOperator[Join[Table["H" -> q, {q, 3}], oracle, Table["H" -> q, {q, 3}]]];
rec["K1-constant-P(000)", djLayers[{}][QuantumState["000"]]["ProbabilitiesList"][[1]] // safe];
rec["K1-balanced-P(000)", djLayers[{"Z" -> 1}][QuantumState["000"]]["ProbabilitiesList"][[1]] // safe];

(* K2 Bernstein-Vazirani: recover hidden s=101 in one query (Z on s-bits) *)
bv = QuantumCircuitOperator[Join[Table["H" -> q, {q, 3}], {"Z" -> 1, "Z" -> 3}, Table["H" -> q, {q, 3}]]];
rec["K2-BV-output", bv[QuantumState["000"]]["StateVector"] // Normal // Chop // safe];

(* K3 Simon: 2 input + 1 output qubit, f(x)=x1 XOR x2, hidden s=11; input register -> y.s=0 *)
simon = QuantumCircuitOperator[{"H" -> 1, "H" -> 2, "CNOT" -> {1, 3}, "CNOT" -> {2, 3}, "H" -> 1, "H" -> 2}];
rec["K3-input-probs", QuantumPartialTrace[simon[QuantumState["000"]], {3}]["ProbabilitiesList"] // Chop // safe];

(* K4 Grover: amplitude amplification, manual G = D.O for n=3 marking |101> *)
nq = 3; NN = 2^nq; tgt = 6;                 (* |101> is index 6 (1-based) *)
Omat = IdentityMatrix[NN]; Omat[[tgt, tgt]] = -1;
sv = ConstantArray[1/Sqrt[NN], NN];
Dmat = 2 Outer[Times, sv, sv] - IdentityMatrix[NN];
Gop = QuantumOperator[Dmat . Omat];
psi0 = QuantumState["UniformSuperposition"[3]];
grov = Nest[Gop[#] &, psi0, 2];
rec["K4-grover-P(101)", grov["ProbabilitiesList"][[tgt]] // safe];
rec["K4-grover-Pmax", Max[grov["ProbabilitiesList"]] // safe];

(* K5 QFT: Fourier circuit equals the DFT; QFT|0...0> is uniform *)
qft = QuantumCircuitOperator["Fourier"[3]] // safe;
rec["K5-qft|000>", QuantumCircuitOperator["Fourier"[3]][QuantumState["000"]]["ProbabilitiesList"] // Chop // safe];
rec["K5-qft-unitary", Chop[(QuantumCircuitOperator["Fourier"[3]]["QuantumOperator"]["MatrixRepresentation"] // Normal)[[1]]] // safe];

(* K6 phase-estimation readout core: inverse QFT maps a phase-kicked register to binary(phi) *)
phi = 3/8;
kicked = QuantumState[Table[Exp[2 Pi I phi k], {k, 0, 7}]/Sqrt[8]];
iqft = QuantumCircuitOperator["Fourier"[3]]["QuantumOperator"]["Dagger"];
rec["K6-readout-probs", iqft[kicked]["ProbabilitiesList"] // Chop // safe];
rec["K6-readout-index", (First[Ordering[iqft[kicked]["ProbabilitiesList"] // Chop, -1]] - 1) // safe];

(* K7 coherent matrix-on-vector (HHL primitive): Pauli-decompose A, multiplexer, ancilla *)
A = {{1, 1/2}, {1/2, 1}}; xvec = {3/5, 4/5};
sig = QuantumOperator[A]["PauliDecompose"];
mp = QuantumCircuitOperator["Multiplexer" @@ Keys[sig]];
anc = QuantumState[Sqrt[Values[sig]]];
qcAx = QuantumCircuitOperator[{anc, QuantumState[xvec] -> {2}, mp, anc["Conjugate"]["Dagger"]}];
rec["K7-pauli-decomp", sig // safe];
rec["K7-coherent-Ax", qcAx[]["AmplitudesList"] // safe];
rec["K7-classical-Ax", A . xvec];

(* K8 VQE: ground-state energy of H = Z via RY(theta) ansatz *)
ans = QuantumCircuitOperator[{"RY"[th] -> 1}, "Parameters" -> {th}];
phi = ans[QuantumState["0"]];
cost[t_] := Re[Tr[QuantumOperator["Z"]["Matrix"] . (ans[<|th -> t|>][QuantumState["0"]]["DensityMatrix"])]];
rec["K8-VQE-min", NMinimize[{cost[t], -Pi <= t <= Pi}, t][[1]] // safe];

(* K9 QAOA / Max-Cut on a triangle: exact max cut from cost Hamiltonian *)
edges = {{1, 2}, {2, 3}, {1, 3}};
zz[{i_, j_}] := QuantumOperator[StringJoin[Table[Which[q == i || q == j, "Z", True, "I"], {q, 3}]]];
Hc = Total[(1/2)(QuantumOperator["III"]["Matrix"] - zz[#]["Matrix"]) & /@ edges];
rec["K9-maxcut", Max[Eigenvalues[Hc // Normal] // Re] // safe];

(* K10 teleportation by deferred measurement: q3 ends in the input state *)
psiIn = QuantumState[{Cos[Pi/5], Sin[Pi/5] Exp[I Pi/3]}];
init = QuantumTensorProduct[psiIn, QuantumState["PhiPlus"]];
tele = QuantumCircuitOperator[{"CNOT" -> {1, 2}, "H" -> 1, "CNOT" -> {2, 3}, "CZ" -> {1, 3}}];
out = tele[init];
q3 = QuantumPartialTrace[out, {1, 2}];
rec["K10-teleport-fidelity", N@QuantumSimilarity[q3, psiIn] // Chop // safe];
rec["K10-q3-purity", N@q3["Purity"] // safe];

(* K11 superdense coding: 2 bits through 1 qubit; all 4 messages decode *)
encode[g_] := QuantumCircuitOperator[Join[{"H" -> 1, "CNOT" -> {1, 2}}, g, {"CNOT" -> {1, 2}, "H" -> 1}]];
decode[g_] := encode[g][QuantumState["00"]]["ProbabilitiesList"];
rec["K11-msg-I(00)", decode[{}] // Chop // safe];
rec["K11-msg-X(01)", decode[{"X" -> 1}] // Chop // safe];
rec["K11-msg-Z(10)", decode[{"Z" -> 1}] // Chop // safe];
rec["K11-msg-XZ(11)", decode[{"Z" -> 1, "X" -> 1}] // Chop // safe];

(* K12 entanglement swapping: measure pair (2,3) in Bell basis -> qubits 1 and 4 entangled *)
twoPairs = QuantumTensorProduct[QuantumState["PhiPlus"], QuantumState["PhiPlus"]];
rec["K12-conc(1,4)-before", N@QuantumEntanglementMonotone[QuantumPartialTrace[twoPairs, {2, 3}], {1}, "Concurrence"] // safe];
bm = QuantumCircuitOperator[{"CNOT" -> {2, 3}, "H" -> 2}][twoPairs];
branch = QuantumMeasurementOperator[{2, 3}][bm]["States"][[1]]["Normalized"];
rec["K12-conc(1,4)-after", N@QuantumEntanglementMonotone[QuantumPartialTrace[branch, {2, 3}], {1}, "Concurrence"] // safe];
