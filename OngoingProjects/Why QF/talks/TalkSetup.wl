(* Definitions behind Talk-Slides-v3: Get this file once before the first slide. *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];

(* === 1. driven, decaying qubit === *)
\[Sigma]m = QuantumOperator[{{0, 1}, {0, 0}}];
\[Rho]0 = QuantumState[Array[\[Rho], {2, 2}]];

(* === 2. Liouvillian of a watched qubit, Mollow spectrum === *)
{id2, X, Y, Z} = Table[PauliMatrix[j], {j, 0, 3}];
liouvillian[H_, leaks_List, d_] := With[{id = IdentityMatrix[d]},
    -I (KroneckerProduct[H, id] - KroneckerProduct[id, Transpose[H]]) +
        Total[Function[c, KroneckerProduct[c, Conjugate[c]] -
            (KroneckerProduct[ConjugateTranspose[c] . c, id] +
                KroneckerProduct[id, Transpose[ConjugateTranspose[c] . c]])/2] /@ leaks]];
lower = (X - I Y)/2;
steady[H_, leaks_] := With[{v = ArrayReshape[First@NullSpace[liouvillian[H, leaks, 2]], {2, 2}]}, v/Tr[v]];
resolventSpectrum[big_, seedMat_, readout_, tones_] := With[{n = Length[big], d = Length[readout]},
    Table[2 Re@Tr[readout . ArrayReshape[If[Abs[w] < 10^-9,
        LeastSquares[Join[N[-big], {1. Flatten[IdentityMatrix[d]]}], Join[Flatten[seedMat], {0.}]],
        LinearSolve[I w IdentityMatrix[n] - big, Flatten[seedMat]]], {d, d}]], {w, tones}]];
mollow[om_, tones_] := With[{big = liouvillian[(om/2) X, {lower}, 2], rss = steady[(om/2) X, {lower}]},
    With[{dm = lower - Tr[lower . rss] id2}, resolventSpectrum[big, dm . rss, ConjugateTranspose[dm], tones]]];

(* === 3. bit-flip code with measurement records as wires === *)
zz = KroneckerProduct[PauliMatrix[3], PauliMatrix[3]];
syn[q1_, q2_] := QuantumMeasurementOperator[
    {QuantumOperator[(KroneckerProduct[id2, id2] + zz)/2, {q1, q2}],
        QuantumOperator[(KroneckerProduct[id2, id2] - zz)/2, {q1, q2}]}, {q1, q2}];
enc = QuantumCircuitOperator[{"CNOT", "CNOT" -> {1, 3}}, "Label" -> "Encode"];
qec = QuantumCircuitOperator[{enc, "RX"[2 \[Theta]] -> 2, syn[1, 2], syn[2, 3],
    "C"["X", {0}, {-1}], "C"["X" -> 2, {0, -1}, {}], "C"["X" -> 3, {-1}, {0}]}];
coded = qec[QuantumTensorProduct[QuantumState[{\[Alpha], \[Beta]}], QuantumState["00"]]];
asQEC = {\[Theta] \[Element] Reals, \[Alpha] \[Element] Reals, \[Beta] \[Element] Reals, \[Alpha]^2 + \[Beta]^2 == 1};

(* === 4. continuously watched particle: Riccati steady state, grid simulation, Kalman filter === *)
riccati = {2 Cc/m - 4 \[Lambda] Sx^2, Sp/m - m \[Omega]^2 Sx - 4 \[Lambda] Sx Cc,
    -2 m \[Omega]^2 Cc + \[Lambda] \[HBar]^2 - 4 \[Lambda] Cc^2};
positiveQ[sol_] := TrueQ@Resolve[ForAll[{m, \[Lambda], \[Omega], \[HBar]},
    m > 0 && \[Lambda] > 0 && \[Omega] >= 0 && \[HBar] > 0,
    (sol[[1]] | sol[[3]]) \[Element] Reals && sol[[1]] > 0 && sol[[3]] > 0], Reals];
{SxSS, CcSS, SpSS} = SelectFirst[FullSimplify@SolveValues[Thread[riccati == 0], {Sx, Cc, Sp}], positiveQ];
asTrap = {\[Lambda] > 0, \[Omega] > 0, m > 0, \[HBar] > 0};
sqRatio[u_] := Sqrt[2 (Sqrt[1 + u^2] - 1)]/u;

hbar = 1; mass = 1;
grid[n_, l_] := With[{dx = l/n}, <|"n" -> n, "x" -> dx Range[-n/2, n/2 - 1],
    "p" -> (2 Pi hbar/l) RotateLeft[Range[-n/2, n/2 - 1], n/2]|>];
ft[psi_] := Fourier[psi, FourierParameters -> {0, -1}];
ift[c_] := InverseFourier[c, FourierParameters -> {0, -1}];
gaussian[g_, x0_, p0_, s_] := Normalize[Exp[-(g["x"] - x0)^2/(4 s^2) + I p0 g["x"]/hbar]];
xMean[g_, psi_] := g["x"] . Abs[psi]^2;
outcome[g_, lam_, dt_][psi_] := RandomChoice[Abs[psi]^2 -> g["x"]] +
    RandomVariate[NormalDistribution[0, 1/Sqrt[4 lam dt]]];
measUpdate[g_, lam_, dt_, xb_][psi_] := Normalize[Exp[-lam dt (g["x"] - xb)^2] psi];
potentialStep[g_, V_, dt_][psi_] := Exp[-I V[g["x"]] dt/(2 hbar)] psi;
kineticStep[g_, dt_][psi_] := ift[Exp[-I g["p"]^2 dt/(2 mass hbar)] ft[psi]];
evolveStep[g_, V_, dt_] := potentialStep[g, V, dt] @* kineticStep[g, dt] @* potentialStep[g, V, dt];
trajectory[g_, V_, lam_, dt_, nt_][psi0_] := NestList[
    evolveStep[g, V, dt][measUpdate[g, lam, dt, Sow@outcome[g, lam, dt][#]][#]] &, psi0, nt];

wR = 1.; lamH = 2.; SxZP = hbar/(2 mass wR); SpZP = hbar mass wR/2;
gH = grid[256, 24.]; dtH = 0.002; ntH = 1500;
{states, {outs}} = BlockRandom[SeedRandom[21];
    Reap[trajectory[gH, 0.5 mass wR^2 #^2 &, lamH, dtH, ntH][gaussian[gH, 0., 0., Sqrt[SxZP]]]]];
xGrid = xMean[gH, #] & /@ states;

num = {\[HBar] -> hbar, m -> mass, \[Omega] -> wR, \[Lambda] -> lamH};
covDrift = Function[{sx, cc, sp}, Evaluate[riccati /. num /. {Sx -> sx, Cc -> cc, Sp -> sp}]];
filterStep[dt_][{mx_, mp_, sx_, cc_, sp_}, dY_] := With[{dWf = dY - 2 Sqrt[lamH] mx dt},
    Join[{mx + mp/mass dt + 2 Sqrt[lamH] sx dWf, mp - mass wR^2 mx dt + 2 Sqrt[lamH] cc dWf},
        {sx, cc, sp} + dt covDrift[sx, cc, sp]]];
mFilter = FoldList[filterStep[dtH], {0., 0., SxZP, 0., SpZP}, 2 Sqrt[lamH] dtH outs][[All, 1]];

(* === 5. linear optics === *)
nsMatrix = {{1 - Sqrt[2], 2^(-1/4), Sqrt[3/Sqrt[2] - 2]}, {2^(-1/4), 1/2, (1 - Sqrt[2])/2},
    {Sqrt[3/Sqrt[2] - 2], (1 - Sqrt[2])/2, Sqrt[2] - 1/2}};
nsGate = QuantumInterferometer[nsMatrix];

(* === 6. stabilizers === *)
code = PauliStabilizer["5QubitCode"];
checks = Most[code["Stabilizers"]];
syndrome[err_] := (1 - (code[err]["Expectation", #] & /@ checks))/2;
errors = Flatten[Table[{p -> q}, {p, {"X", "Y", "Z"}}, {q, 1, 5}], 1];
magicState[n_] := PauliStabilizer[n][Join[{"H" -> 1}, Table["CNOT" -> {i, i + 1}, {i, 1, n - 1}],
    {"T" -> 1, "T" -> 2, "T" -> 3}]];

(* === 7. DRAG pulse on a five-level transmon === *)
envelope = A (Pi/\[Tau]) (1 - Cos[2 Pi t/\[Tau]]);
hamiltonian[c_] := (\[Alpha]/2) (c["Dagger"] @ c["Dagger"] @ c @ c) + \[Delta] (c["Dagger"] @ c) +
    (envelope/2) (c + c["Dagger"]) - (\[Lambda] D[envelope, t]/(2 \[Alpha])) (I (c["Dagger"] - c));
anharmonicity = -2 Pi 0.2;
pulse[n_, {amplitude_, drag_, detuning_, length_}] := QuantumEvolve[
    QuantumOperator[hamiltonian[AnnihilationOperator[n]], "Parameters" -> {\[Alpha], A, \[Lambda], \[Delta], \[Tau]}][
        anharmonicity, amplitude, drag, detuning, length],
    FockState[0, n], {t, 0, length}, Method -> "ExplicitRungeKutta", PrecisionGoal -> 12, AccuracyGoal -> 14];
finalPopulations[n_, {amplitude_, drag_, detuning_, length_}] :=
    pulse[n, {amplitude, drag, detuning, length}][length]["ProbabilitiesList"];

(* === 8. boson sampling === *)
modeList[occ_List] := Catenate@MapIndexed[ConstantArray[First[#2], #1] &, occ];
bosonProb[w_, s_, t_] := Abs[Permanent[w[[modeList[t], modeList[s]]]]]^2/Times @@ Factorial[Join[s, t]];
distProb[w_, s_, t_] := Permanent[(Abs[w]^2)[[modeList[t], modeList[s]]]]/Times @@ Factorial[t];
fermionProb[w_, s_, t_] := Abs[Det[w[[modeList[t], modeList[s]]]]]^2;
Wbs = QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2]["TransferMatrix"];
occupations[n_] := PadRight[#, 2 n] & /@ <|"collision-free" -> ConstantArray[1, n],
    "pairs" -> Append[ConstantArray[2, Quotient[n, 2]], Mod[n, 2]], "bunched" -> {n}|>;
randomInterferometer[n_] := BlockRandom[SeedRandom[1107 + 1000 n]; QuantumInterferometer["Random"[2 n]]];
timeAmplitude[q_, s_List] := With[{t = RotateRight[s, Total[s]]},
    ClearSystemCache[]; First@AbsoluteTiming[q["Amplitude", s -> t]]];
