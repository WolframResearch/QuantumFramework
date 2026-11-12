<!-- Build: MarkdownToNotebook["Talk-Slides-50min.md", "Talk-Slides-50min.nb", "EvaluateSeparator" -> None]. The slides share one evaluation context, so later slides use definitions made on earlier ones. -->

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

# What the Wolfram Language and QuantumFramework Do for Physicists

#| style: Subtitle

Symbolic results, exact checks, and the pictures that go with them

#| style: Subsubtitle

Mohammad Bahrami, Wolfram Research

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 1. Open Quantum Systems, Solved Once for Every Qubit

A driven qubit with relaxation and dephasing, every rate and every entry of the initial state a symbol.

#| collapse: input
```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
$Assumptions = \[Gamma]1 > 0 && \[Gamma]\[Phi] >= 0 && \[CapitalOmega] \[Element] Reals && t >= 0;
\[Sigma]m = QuantumOperator[{{0, 1}, {0, 0}}];
\[Rho]0 = QuantumState[Array[\[Rho], {2, 2}]];
\[Rho]t = QuantumEvolve[\[CapitalOmega]/2 QuantumOperator["Z"],
   {Sqrt[\[Gamma]1] \[Sigma]m, Sqrt[\[Gamma]\[Phi]] QuantumOperator["Z"]}, \[Rho]0, t];
FullSimplify /@ Normal[\[Rho]t["Computational"]["DensityMatrix"]] // MatrixForm
```

Populations relax at $\gamma_1$; coherences precess at $\Omega$ and decay at $\gamma_1/2+2\gamma_\phi$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### $T_2\le2T_1$ as a theorem

A statement about all admissible rates, decided by `Reduce`.

#| collapse: input
```wl
Reduce[ForAll[{\[Gamma]1, \[Gamma]\[Phi]}, \[Gamma]1 > 0 && \[Gamma]\[Phi] >= 0,
  1/(\[Gamma]1/2 + 2 \[Gamma]\[Phi]) <= 2/\[Gamma]1]]
```

Derived from the solution, not quoted, with equality exactly when pure dephasing vanishes.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A coherent error, with measurement records as wires

Encode, rotate qubit 2 by an arbitrary angle, measure two parities onto record wires, and correct by gates controlled on those records.

#| hidden: true
```wl
$Assumptions = True;
id2 = IdentityMatrix[2]; zz = KroneckerProduct[PauliMatrix[3], PauliMatrix[3]];
syn[q1_, q2_] := QuantumMeasurementOperator[
   {QuantumOperator[(KroneckerProduct[id2, id2] + zz)/2, {q1, q2}],
    QuantumOperator[(KroneckerProduct[id2, id2] - zz)/2, {q1, q2}]}, {q1, q2}];
enc = QuantumCircuitOperator[{"CNOT", "CNOT" -> {1, 3}}, "Label" -> "Encode"];
qec = QuantumCircuitOperator[{enc, "RX"[2 \[Theta]] -> 2, syn[1, 2], syn[2, 3],
    "C"["X", {0}, {-1}], "C"["X" -> 2, {0, -1}, {}], "C"["X" -> 3, {-1}, {0}]}];
coded = qec[QuantumTensorProduct[QuantumState[{\[Alpha], \[Beta]}], QuantumState["00"]]];
```

```wl
#| collapse: input
qec["Diagram"]
```

The classical record is part of the circuit, so the whole protocol is one object.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The error is digitized and corrected exactly

Syndrome distribution and final fidelity, for a symbolic data qubit and a symbolic error angle.

#| collapse: input
```wl
as = {\[Theta] \[Element] Reals, \[Alpha] \[Element] Reals, \[Beta] \[Element] Reals, \[Alpha]^2 + \[Beta]^2 == 1};
Assuming[as, FullSimplify[ComplexExpand[
   Diagonal[Normal[QuantumPartialTrace[coded, {3, 4, 5}]["DensityMatrix"]]]]]]
```

#| collapse: input
```wl
Assuming[as, FullSimplify[ComplexExpand[QuantumSimilarity[QuantumPartialTrace[coded, {1, 2}],
    enc[QuantumTensorProduct[QuantumState[{\[Alpha], \[Beta]}], QuantumState["00"]]], "Fidelity"]]]]
```

The continuous rotation becomes "no error" with probability $\cos^2\theta$ or "qubit 2 flipped" with $\sin^2\theta$, and the correction restores every input at every angle.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 2. The Quantum Zeno Crossover as an Exceptional Point

A qubit Rabi-driven at $\Omega$ while $\hat\sigma_z$ is watched at strength $k$: the full Liouvillian spectrum, symbolic.

#| hidden: true
```wl
{id2, X, Y, Z} = Table[PauliMatrix[j], {j, 0, 3}];
liouvillian[H_, leaks_List, d_] := With[{id = IdentityMatrix[d]},
   -I (KroneckerProduct[H, id] - KroneckerProduct[id, Transpose[H]]) +
    Total[Function[c, KroneckerProduct[c, Conjugate[c]] -
        (KroneckerProduct[ConjugateTranspose[c] . c, id] +
           KroneckerProduct[id, Transpose[ConjugateTranspose[c] . c]])/2] /@ leaks]];
```

#| collapse: input
```wl
Eigenvalues@FullSimplify[liouvillian[(\[CapitalOmega]/2) X, {Sqrt[k] Z}, 2], \[CapitalOmega] > 0 && k > 0]
```

#| collapse: input
```wl
Normal@Series[-k + Sqrt[k^2 - \[CapitalOmega]^2], {k, Infinity, 1}]
```

The oscillating pair coalesces at $k=\Omega$ and splits into two real rates; the slow one is the Zeno freeze $\Omega^2/2k$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Watching harder freezes the qubit

Rabi frequency and slowest decay rate against measurement strength.

#| collapse: input
```wl
Plot[{Re@Sqrt[1 - r^2], r - Re@Sqrt[r^2 - 1]}, {r, 0, 3}, Frame -> True, GridLines -> {{1}, None},
 PlotLegends -> {"Rabi frequency / \[CapitalOmega]", "slowest decay rate / \[CapitalOmega]"},
 FrameLabel -> {"k/\[CapitalOmega]", "Liouvillian eigenvalue / \[CapitalOmega]"}]
```

Measurement speeds relaxation only up to the exceptional point; past it, it slows the dynamics down.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The Mollow triplet

Fluorescence of a resonantly driven atom by quantum regression: one linear solve per frequency on the same Liouvillian.

#| hidden: true
```wl
lower = (X - I Y)/2;
steady[H_, leaks_] := With[{v = ArrayReshape[First@NullSpace[liouvillian[H, leaks, 2]], {2, 2}]}, v/Tr[v]];
resolventSpectrum[big_, seedMat_, readout_, tones_] := With[{n = Length[big], d = Length[readout]},
   Table[2 Re@Tr[readout . ArrayReshape[If[Abs[w] < 10^-9,
       LeastSquares[Join[N[-big], {1. Flatten[IdentityMatrix[d]]}], Join[Flatten[seedMat], {0.}]],
       LinearSolve[I w IdentityMatrix[n] - big, Flatten[seedMat]]], {d, d}]], {w, tones}]];
```

#| collapse: input
```wl
mollow[om_, tones_] := With[{big = liouvillian[(om/2) X, {lower}, 2], rss = steady[(om/2) X, {lower}]},
   With[{dm = lower - Tr[lower . rss] id2}, resolventSpectrum[big, dm . rss, ConjugateTranspose[dm], tones]]];
tones = Range[-9, 9, 0.05];
{weak, strong} = mollow[#, tones] & /@ {0.6, 6.};
ListLinePlot[{Transpose[{tones, weak/Max[weak]}], Transpose[{tones, strong/Max[strong]}]},
 Frame -> True, PlotRange -> All, PlotLegends -> {"\[CapitalOmega] = 0.6\[Gamma]", "\[CapitalOmega] = 6\[Gamma]"},
 FrameLabel -> {"frequency - drive frequency", "inelastic spectrum (scaled to peak)"}]
```

Strong driving dresses the atom and splits one line into three, with sidebands at $\pm\Omega$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Fluorescence squeezed below the vacuum, exactly

The homodyne noise at zero frequency is $S(0)=1+8\Omega^2(4\Omega^2-1)/(1+2\Omega^2)^3$; minimize it over the drive.

#| collapse: input
```wl
Minimize[{1 + 8 \[CapitalOmega]^2 (4 \[CapitalOmega]^2 - 1)/(1 + 2 \[CapitalOmega]^2)^3, \[CapitalOmega] > 0}, \[CapitalOmega]]
```

Below the vacuum for every $\Omega<\gamma/2$: nonclassical light from a single atom, visible only to phase-sensitive detection.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 3. Continuous Measurement: Squeezing in Closed Form

A trapped particle whose position is watched continuously: the conditional covariances obey a Riccati equation, solved with $\lambda,\omega,m,\hbar$ symbolic.

#| collapse: input
```wl
riccati = {2 Cc/m - 4 \[Lambda] Sx^2, Sp/m - m \[Omega]^2 Sx - 4 \[Lambda] Sx Cc,
   -2 m \[Omega]^2 Cc + \[Lambda] \[HBar]^2 - 4 \[Lambda] Cc^2};
steady = SolveValues[Thread[riccati == 0], {Sx, Cc, Sp}] // FullSimplify;
positiveQ[sol_] := TrueQ@Resolve[ForAll[{m, \[Lambda], \[Omega], \[HBar]},
     m > 0 && \[Lambda] > 0 && \[Omega] >= 0 && \[HBar] > 0,
     (sol[[1]] | sol[[3]]) \[Element] Reals && sol[[1]] > 0 && sol[[3]] > 0], Reals];
{SxSS, CcSS, SpSS} = SelectFirst[steady, positiveQ];
as = {\[Lambda] > 0, \[Omega] > 0, m > 0, \[HBar] > 0};
Simplify[SxSS SpSS - CcSS^2, as]
```

The steady conditional state is pure.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Below the zero-point at every strength

In ground-state units the position variance depends only on $\mu=2\hbar\lambda/m\omega^2$; the cell proves it is below one for all $\mu$ and finds both limits.

#| collapse: input
```wl
sqRatio[u_] := Sqrt[2 (Sqrt[1 + u^2] - 1)]/u;
{Simplify[SxSS/(\[HBar]/(2 m \[Omega])) == sqRatio[2 \[HBar] \[Lambda]/(m \[Omega]^2)], as],
 Resolve[ForAll[u, u > 0, sqRatio[u] < 1], Reals],
 Asymptotic[sqRatio[u], {u, 0, 2}, Assumptions -> u >= 0],
 Asymptotic[sqRatio[u], {u, Infinity, 1}, Assumptions -> u >= 0]}
```

Weak measurement squeezes as $1-\mu^2/8$, strong as $\sqrt{2/\mu}$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What an imperfect detector costs

A detector catching a fraction $\eta$ of the signal: steady purity, and the strength needed to squeeze.

#| collapse: input
```wl
{SxE, CcE, SpE} = {SxSS, CcSS, SpSS} /. {\[Lambda] -> \[Eta] \[Lambda], \[HBar] -> \[HBar]/Sqrt[\[Eta]]};
{Simplify[\[HBar]/(2 Sqrt[SxE SpE - CcE^2]), {\[HBar] > 0, 0 < \[Eta] <= 1, \[Lambda] > 0, \[Omega] > 0, m > 0}],
 Simplify@Reduce[sqRatio[Sqrt[\[Eta]] u]/Sqrt[\[Eta]] <= 1 && u > 0 && 0 < \[Eta] < 1, u, Reals]}
```

The purity is exactly $\sqrt\eta$ whatever the trap, mass or strength, and squeezing survives past a threshold that diverges as the detector fails.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Five numbers track 256 amplitudes

The full conditional wavefunction on a grid, against a quantum Kalman filter fed only the noisy readings.

#| hidden: true
```wl
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
```

#| hidden: true
```wl
wR = 1.; lamH = 2.; SxZP = hbar/(2 mass wR); SpZP = hbar mass wR/2;
gH = grid[256, 24.]; dtH = 0.002; ntH = 1500;
{states, {outs}} = BlockRandom[SeedRandom[21];
   Reap[trajectory[gH, 0.5 mass wR^2 #^2 &, lamH, dtH, ntH][gaussian[gH, 0., 0., Sqrt[SxZP]]]]];
xGrid = xMean[gH, #] & /@ states;
```

#| collapse: input
```wl
num = {\[HBar] -> hbar, m -> mass, \[Omega] -> wR, \[Lambda] -> lamH};
covDrift = Function[{sx, cc, sp}, Evaluate[riccati /. num /. {Sx -> sx, Cc -> cc, Sp -> sp}]];
filterStep[dt_][{mx_, mp_, sx_, cc_, sp_}, dY_] := With[{dWf = dY - 2 Sqrt[lamH] mx dt},
   Join[{mx + mp/mass dt + 2 Sqrt[lamH] sx dWf, mp - mass wR^2 mx dt + 2 Sqrt[lamH] cc dWf},
    {sx, cc, sp} + dt covDrift[sx, cc, sp]]];
mFilter = FoldList[filterStep[dtH], {0., 0., SxZP, 0., SpZP}, 2 Sqrt[lamH] dtH outs][[All, 1]];
ListLinePlot[{xGrid, mFilter}, DataRange -> {0, ntH dtH},
 PlotStyle -> {Directive[Thick, Gray], Directive[Red, Dashed]},
 PlotLegends -> {"grid (256 amplitudes)", "Kalman filter (5 numbers)"}, Frame -> True,
 FrameLabel -> {"t", "\[LeftAngleBracket]x\[RightAngleBracket]"}]
```

In a quadratic potential the conditional state is a Gaussian whose center follows the record, so the wavefunction is never needed.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 4. Linear Optics with Symbolic Photon Numbers

A beam splitter acting on $|n_1,n_2\rangle$, with every photon number a symbol.

#| collapse: input
```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
bs = QuantumInterferometer[{"BS"[th, ph]}];
bs["Amplitude", {n1, n2} -> {t1, n1 + n2 - t1}]
```

A Wigner $D$-function: the beam splitter is a rotation of a spin $(n_1+n_2)/2$, and photon numbers are its magnetic quantum numbers.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Hong-Ou-Mandel for every $n$

$|n,n\rangle$ on a balanced splitter.

#| collapse: input
```wl
hom = QuantumInterferometer[{"BS"[Pi/4, 0]}]["Amplitude", {n, n} -> {k, 2 n - k}];
{hom, Table[hom /. {n -> m, k -> j}, {m, 1, 3}, {j, 1, 2 m - 1, 2}], Table[Abs[hom /. {n -> 2, k -> j}]^2, {j, 0, 4}]}
```

One formula for all $n$: every odd split is exactly zero, and $n=1$ is the HOM dip.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Squeezed pairs

Opposite squeezers on a balanced splitter, $r$ symbolic.

#| collapse: input
```wl
symbolicPair = QuantumInterferometer[{"S"[r], "S"[-r] -> 2, "BS"[Pi/4]}];
{FullSimplify[Table[symbolicPair["Probability", {0, 0} -> {n, n}], {n, 0, 2}], r > 0],
 FullSimplify[symbolicPair["Probability", {0, 0} -> {2, 0}], r > 0]}
```

```wl
#| collapse: input
symbolicPair["Diagram"]
```

Photons come only in pairs, with the geometric law $\operatorname{sech}^2r\tanh^{2n}r$: a two-mode squeezed vacuum.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The KLM nonlinear sign gate, exactly

A three-mode real unitary with one auxiliary photon, heralded on the auxiliary pattern.

#| collapse: input
```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
nsMatrix = {{1 - Sqrt[2], 2^(-1/4), Sqrt[3/Sqrt[2] - 2]}, {2^(-1/4), 1/2, (1 - Sqrt[2])/2},
   {Sqrt[3/Sqrt[2] - 2], (1 - Sqrt[2])/2, Sqrt[2] - 1/2}};
nsGate = QuantumInterferometer[nsMatrix];
{RootReduce[nsMatrix . Transpose[nsMatrix]] == IdentityMatrix[3],
 RootReduce[Table[nsGate["Amplitude", {n, 1, 0} -> {n, 1, 0}], {n, 0, 3}]]}
```

```wl
#| collapse: input
nsGate["Diagram"]
```

A sign flip on $|2\rangle$ alone, which no phase shifter can do, and the three-photon amplitude shows where the gate stops working.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Four-photon suppression needs an unbalanced splitter

The return amplitude of $|2,2\rangle$ on a splitter of angle $\theta$, and the reflectivities that kill it.

#| collapse: input
```wl
fourBS = QuantumInterferometer[{"BS"[\[Theta]]}];
fa = FullSimplify[fourBS["Amplitude", {2, 2} -> {2, 2}], \[Theta] \[Element] Reals];
{FullSimplify[fa == 1 - 6 Sin[\[Theta]]^2 Cos[\[Theta]]^2, \[Theta] \[Element] Reals],
 R /. Solve[1 - 6 R (1 - R) == 0 && 0 < R < 1, R, Reals]}
```

$1-6R(1-R)$ vanishes at $R=(3\pm\sqrt3)/6$, two settings that are mirror images under $R\leftrightarrow1-R$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 5. Stabilizers: Codes and a Thousand Qubits

The five-qubit code: the syndrome of each of the fifteen single-qubit errors.

#| collapse: input
```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
code = PauliStabilizer["5QubitCode"];
checks = Most[code["Stabilizers"]];
syndrome[err_] := (1 - (code[err]["Expectation", #] & /@ checks))/2;
errors = Flatten[Table[{p -> q}, {p, {"X", "Y", "Z"}}, {q, 1, 5}], 1];
AssociationMap[syndrome, errors]
```

#| collapse: input
```wl
Sort[syndrome /@ errors] === Sort@DeleteCases[Tuples[{0, 1}, 4], {0, 0, 0, 0}]
```

Each error has its own nonzero syndrome, none missing and none repeated: what "perfect" means.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A thousand-qubit GHZ state

Entanglement entropy across one qubit and across half the chain, with no state vector.

#| collapse: input
```wl
ghz = PauliStabilizer[1000][Prepend[Table["CNOT" -> {q - 1, q}, {q, 2, 1000}], "H"]];
{AbsoluteTiming[ghz["Entropy", {1}]], ghz["Entropy", Range[500]]}
```

One bit either way: GHZ entanglement is a single shared parity, however the chain is cut.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The price of magic

T gates on a GHZ backbone at growing width, and a single T on $|+\rangle$.

#| collapse: input
```wl
magicState[n_] := PauliStabilizer[n][Join[{"H" -> 1}, Table["CNOT" -> {i, i + 1}, {i, 1, n - 1}],
     {"T" -> 1, "T" -> 2, "T" -> 3}]];
{Table[{n, magicState[n]["Length"], 2^n}, {n, {4, 10, 20}}],
 FullSimplify[(1 + PauliStabilizer[1][{"H", "T"}]["Expectation", "X"])/2]}
```

The cost grows with the number of T gates, not qubits, and the irrational Born probability comes out exactly.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 6. Snyder's Quantum Spacetime

Coordinates that do not commute, $[x^\mu,x^\nu]=ia^2M^{\mu\nu}$, yet keep Lorentz symmetry. Declare the algebra once and test its Casimir.

#| hidden: true
```wl
eta[m_, n_] := DiagonalMatrix[{1, -1, -1, -1}][[m + 1, n + 1]];
MGen[m_, n_] := Which[m == n, 0, m < n, M[m, n], True, -M[n, m]];
gen = Join[Table[x[m], {m, 0, 3}], Flatten@Table[M[m, n], {m, 0, 3}, {n, m + 1, 3}]];
xxRel = Flatten@Table[Commutator[x[n], x[m]] + I a^2 M[m, n], {m, 0, 3}, {n, m + 1, 3}];
mxRel = Flatten@Table[Commutator[M[m, n], x[r]] - I (eta[n, r] x[m] - eta[m, r] x[n]),
    {m, 0, 3}, {n, m + 1, 3}, {r, 0, 3}];
mmRel = Replace[Subsets[Flatten[Table[{m, n}, {m, 0, 3}, {n, m + 1, 3}], 1], {2}],
   {{m_, n_}, {r_, s_}} :> Commutator[M[r, s], M[m, n]] +
     I (eta[n, r] MGen[m, s] - eta[n, s] MGen[m, r] - eta[m, r] MGen[n, s] + eta[m, s] MGen[n, r]), {1}];
alg = NonCommutativeAlgebra[<|"Generators" -> {gen},
    "CommutationRelations" -> Join[xxRel, mxRel, mmRel], "ScalarVariables" -> {a}|>];
```

#| collapse: input
```wl
C1 = Sum[eta[m, m] eta[n, n] MGen[m, n] ** MGen[m, n], {m, 0, 3}, {n, m + 1, 3}];
xSq = Sum[eta[m, m] x[m] ** x[m], {m, 0, 3}];
{Union@Table[NonCommutativeExpand[Commutator[xSq - a^2 C1, x[r]], alg], {r, 0, 3}],
 Union@Flatten@Table[NonCommutativeExpand[Commutator[xSq - a^2 C1, MGen[r, s]], alg], {r, 0, 3}, {s, r + 1, 3}],
 NonCommutativeExpand[Commutator[xSq, x[0]], alg]}
```

$x^2-a^2C_1$ is central, and the last entry is exactly the residue the $a^2$ term cancels.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A minimal length

The deformed uncertainty relation $\Delta x\ge(1+a^2\Delta p^2)/2\Delta p$.

#| collapse: input
```wl
Assuming[a > 0, Minimize[{(1 + a^2 dp^2)/(2 dp), dp > 0}, dp] // FullSimplify]
```

No event can be localized below the algebra's length $a$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 7. DRAG: Fast Gates on a Transmon

Five transmon levels, a short π pulse with a DRAG term, and a scan of the DRAG coefficient.

#| hidden: true
```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
envelope = A (Pi/\[Tau]) (1 - Cos[2 Pi t/\[Tau]]);
hamiltonian[c_] := (\[Alpha]/2) (c["Dagger"] @ c["Dagger"] @ c @ c) + \[Delta] (c["Dagger"] @ c) +
   (envelope/2) (c + c["Dagger"]) - (\[Lambda] D[envelope, t]/(2 \[Alpha])) (I (c["Dagger"] - c));
anharmonicity = -2 Pi 0.2;
finalPopulations[n_, {amplitude_, drag_, detuning_, length_}] := QuantumEvolve[
   QuantumOperator[hamiltonian[AnnihilationOperator[n]], "Parameters" -> {\[Alpha], A, \[Lambda], \[Delta], \[Tau]}][
    anharmonicity, amplitude, drag, detuning, length],
   FockState[0, n], {t, 0, length}, Method -> "ExplicitRungeKutta", PrecisionGoal -> 12, AccuracyGoal -> 14][length][
  "ProbabilitiesList"];
```

#| collapse: input
```wl
ListLogPlot[Transpose@Table[With[{p = finalPopulations[5, {1, drag, 0, 8}]},
    {{drag, Total[p[[3 ;;]]]}, {drag, 1 - p[[2]]}}], {drag, 0, 2, 0.1}],
 Joined -> True, PlotMarkers -> Automatic, PlotLegends -> {"leakage", "transfer error"},
 Frame -> True, FrameLabel -> {"\[Lambda]", "population"}]
```

Leakage is smallest near the weak-drive value, the transfer error well below it: suppressing leakage and completing the flip are different targets.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 8. Why Boson Sampling Is Hard

Bosons, classical particles and fermions on one balanced beam splitter: permanent, permanent of $|W|^2$, determinant.

#| collapse: input
```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
modeList[occ_List] := Catenate@MapIndexed[ConstantArray[First[#2], #1] &, occ];
bosonProb[w_, s_, t_] := Abs[Permanent[w[[modeList[t], modeList[s]]]]]^2/Times @@ Factorial[Join[s, t]];
distProb[w_, s_, t_] := Permanent[(Abs[w]^2)[[modeList[t], modeList[s]]]]/Times @@ Factorial[t];
fermionProb[w_, s_, t_] := Abs[Det[w[[modeList[t], modeList[s]]]]]^2;
Wbs = QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2]["TransferMatrix"];
Table[Through[{bosonProb, distProb, fermionProb}[Wbs, {1, 1}, t]], {t, {{2, 0}, {1, 1}, {0, 2}}}]
```

Bosons bunch, classical particles split, fermions antibunch: the only difference is the sign attached to exchanging two histories.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What one amplitude costs

Collision-free, paired and fully bunched inputs in random networks, one amplitude each, timed.

#| collapse: input
```wl
occupations[n_] := PadRight[#, 2 n] & /@ <|"collision-free" -> ConstantArray[1, n],
    "pairs" -> Append[ConstantArray[2, Quotient[n, 2]], Mod[n, 2]], "bunched" -> {n}|>;
randomInterferometer[n_] := BlockRandom[SeedRandom[1107 + 1000 n]; QuantumInterferometer["Random"[2 n]]];
timeAmplitude[q_, s_List] := With[{t = RotateRight[s, Total[s]]},
   ClearSystemCache[]; First@AbsoluteTiming[q["Amplitude", s -> t]]];
AssociationMap[Function[n, With[{q = randomInterferometer[n]}, timeAmplitude[q, #] & /@ occupations[n]]],
 {6, 12, 18, 20}]
```

Generic amplitudes pay the exponential cost of Ryser and Glynn while structured ones do not, which is why the hardness argument needs random networks.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Also in the Wolfram Language

### The scale of a self-gravitating quantum particle

Balance quantum spreading against self-gravity in the Schrodinger-Newton equation.

#| collapse: input
```wl
Solve[D[\[HBar]^2/(m l^2) - G m^2/l, l] == 0, l]
```

$\ell=\hbar^2/Gm^3$: astronomical for an electron, laboratory-sized for a heavy nanoparticle.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Continuous error correction: syndromes from noisy records

Four probabilities, updated by Bayes' rule from two weak stabilizer records.

#| hidden: true
```wl
{id2, X, Y, Z} = Table[PauliMatrix[j], {j, 0, 3}];
op3[o_, i_] := KroneckerProduct @@ ReplacePart[{id2, id2, id2}, i -> o];
ket[bits_] := Normal@SparseArray[FromDigits[bits, 2] + 1 -> 1, 8];
S1 = op3[Z, 1] . op3[Z, 2]; S2 = op3[Z, 2] . op3[Z, 3];
signs = {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
proj = (IdentityMatrix[8] + #[[1]] S1) . (IdentityMatrix[8] + #[[2]] S2)/4 & /@ signs;
```

#| hidden: true
```wl
measurementStep[H_, watched_, effs_, unwatched_, dt_] :=
  With[{id = IdentityMatrix[Length[H]],
    drift = I H + Total[ConjugateTranspose[#] . # & /@ Join[watched, unwatched]]/2,
    corr = Total[MapThread[#1 (#2 . #2) &, {effs, watched}]]/2},
   Function[{rho, dJ}, With[{sig = Total[MapThread[Sqrt[#2] #1 #3 &, {dJ, effs, watched}]]},
     With[{M = id - drift dt + sig + sig . sig/2 - corr dt},
      With[{top = M . rho . ConjugateTranspose[M] +
          dt Total[# . rho . ConjugateTranspose[#] & /@ unwatched] +
          dt Total[MapThread[(1 - #2) #1 . rho . ConjugateTranspose[#1] &, {watched, effs}]]},
       top/Re@Tr[top]]]]]];
measurementRecord[rho_, watched_, effs_, dt_, kick_] :=
  MapThread[Sqrt[#3] Re@Tr[(#1 + ConjugateTranspose[#1]) . rho] dt + #2 &, {watched, kick, effs}];
trajectory[rho0_, H_, watched_, effs_, unwatched_, dt_, tf_, seed_] :=
  BlockRandom[SeedRandom[seed];
   With[{step = measurementStep[H, watched, effs, unwatched, dt],
     kicks = RandomVariate[NormalDistribution[0, Sqrt[dt]], {Round[tf/dt], Length[watched]}]},
    With[{states = FoldList[step[#1, measurementRecord[#1, watched, effs, dt, #2]] &, rho0, kicks]},
     <|"states" -> states,
      "record" -> MapThread[measurementRecord[#1, watched, effs, dt, #2] &, {Most[states], kicks}]|>]]];
```

#| hidden: true
```wl
kappa = 1.; eta = 1.; gamma = 0.05; dt = 0.004;
run = trajectory[KroneckerProduct[#, #] &@ket[{0, 0, 0}], 0 IdentityMatrix[8],
   Sqrt[kappa] {S1, S2}, {eta, eta}, Sqrt[gamma] (op3[X, #] & /@ Range[3]), dt, 40., 3];
```

#| hidden: true
```wl
syndromeFilter[p0_, record_] := With[{c = 2 Sqrt[eta kappa]},
   FoldList[With[{w = (1/4 + Exp[-4 gamma dt] (#1 - 1/4)) Exp[c signs . #2]}, w/Total[w]] &,
    p0, record]];
pFilter = syndromeFilter[{1., 0., 0., 0.}, run["record"]];
```

#| collapse: input
```wl
ListLinePlot[Transpose[pFilter], DataRange -> {0, 40}, PlotRange -> {-0.05, 1.05},
 PlotLegends -> {"code", "X1 flipped", "X2 flipped", "X3 flipped"}, Frame -> True,
 FrameLabel -> {"time", "syndrome probability"}]
```

Each sample is buried in noise, yet the filter follows every bit flip.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### NDSolve integrating a monitored qubit

A stochastic equation turned into a random ordinary one and handed to `NDSolve`.

```wl
#| hidden: true
Clear[zz]
```

#| hidden: true
```wl
driftStrat[{x_, y_, z_, q_}, {gCI_, gBA_, gd_, gphi_, g1_, om_}] := With[{g2 = g1/2 + gphi + gd},
   {((gBA + gCI)/2 - g2 - gCI z^2) x - Sqrt[gBA gCI] y z,
    ((gBA + gCI)/2 - g2 - gCI z^2) y + Sqrt[gBA gCI] x z - om z,
    g1 (1 - z) + gCI z (1 - z^2) + om y,
    Sqrt[gCI] z}];
diffStrat[{x_, y_, z_, q_}, {gCI_, gBA_, __}] :=
  {-Sqrt[gCI] x z - Sqrt[gBA] y, -Sqrt[gCI] y z + Sqrt[gBA] x, Sqrt[gCI] (1 - z^2), 1};
driftIto[{x_, y_, z_}, {_, _, gd_, gphi_, g1_, om_}] := With[{g2 = g1/2 + gphi + gd},
   {-g2 x, -g2 y - om z, g1 (1 - z) + om y}];
rs = {Subscript[\[CapitalGamma], "CI"], Subscript[\[CapitalGamma], "BA"], Subscript[\[CapitalGamma], "d"],
   Subscript[\[Gamma], "\[Phi]"], Subscript[\[Gamma], "1"], Subscript[\[CapitalOmega], "x"]};
```

#| hidden: true
```wl
wienerFun[tf_, h_, seed_] := BlockRandom[SeedRandom[seed];
   With[{grid = Subdivide[0., tf, Ceiling[tf/h]]},
    Interpolation[Transpose[{grid, Prepend[Accumulate[
         RandomVariate[NormalDistribution[0, Sqrt[grid[[2]]]], Length[grid] - 1]], 0.]}],
     InterpolationOrder -> 3]]];
oneTraj[r_, init_, tf_, h_, seed_] := Module[{Wf = wienerFun[tf, h, seed], v = {xx, yy, zz, QQ}, t},
   NDSolveValue[Join[
     Thread[Through[(Derivative[1] /@ v)[t]] ==
       driftStrat[Through[v[t]], r] + diffStrat[Through[v[t]], r] Wf'[t]],
     Thread[Through[v[0]] == init]], v, {t, 0, tf}, MaxStepSize -> h/4]];
```

#| collapse: input
```wl
tr = oneTraj[{0.252, 0.1, 0.176, 1., 1., 3.}, {0., 0., 1., 0.}, 3., 1/120., 7];
Plot[{tr[[3]][t], tr[[4]][t]}, {t, 0, 3}, Frame -> True, PlotLegends -> {"z(t)", "Q(t)"}]
```

The monitored component and the readout move together because one noise drives both.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Shor's discrete logarithm on a toy elliptic curve

Elliptic-curve arithmetic from scratch, an exact simulation of the two-register quantum step, and twenty shots.

#| hidden: true
```wl
EllipticCurve[a_, b_, p_]["Add"][O, q_] := q;
EllipticCurve[a_, b_, p_]["Add"][q_, O] := q;
EllipticCurve[a_, b_, p_]["Add"][{x_, y1_}, {x_, y2_}] /; Mod[y1 + y2, p] == 0 := O;
EllipticCurve[a_, b_, p_]["Add"][{x_, y_}, {x_, y_}] := With[{s = Mod[(3 x^2 + a) PowerMod[2 y, -1, p], p]},
   With[{xr = Mod[s^2 - 2 x, p]}, {xr, Mod[s (x - xr) - y, p]}]];
EllipticCurve[a_, b_, p_]["Add"][{x1_, y1_}, {x2_, y2_}] := With[{s = Mod[(y2 - y1) PowerMod[x2 - x1, -1, p], p]},
   With[{xr = Mod[s^2 - x1 - x2, p]}, {xr, Mod[s (x1 - xr) - y1, p]}]];
EllipticCurve[a_, b_, p_]["Mult"][_, 0] := O;
EllipticCurve[a_, b_, p_]["Mult"][pt_, n_Integer] := With[{add = EllipticCurve[a, b, p]["Add"]},
   Fold[With[{d = add[#1, #1]}, If[#2 == 1, add[d, pt], d]] &, O, IntegerDigits[n, 2]]];
EllipticCurve[a_, b_, p_]["Points"] := Catenate@Table[With[{rhs = Mod[x^3 + a x + b, p]},
     Thread[{x, Select[Range[0, p - 1], Mod[#^2, p] == rhs &]}]], {x, 0, p - 1}];
```

#| hidden: true
```wl
ec = EllipticCurve[1, 1, 23];
pts = ec["Points"]; r = Length[pts] + 1;
G = SelectFirst[pts, Length@DeleteDuplicates@Table[ec["Mult"][#, k], {k, 0, r - 1}] == r &];
Q = ec["Mult"][G, 7];
{r, G, Q}
```

#| hidden: true
```wl
elts = Prepend[pts, O]; idx = AssociationThread[elts -> Range[0, r - 1]];
t = 7; nn = 2^t;
f = Table[idx[ec["Add"][ec["Mult"][G, a], ec["Mult"][Q, b]]], {a, 0, nn - 1}, {b, 0, nn - 1}];
probs = Total[Table[Abs[Fourier[Boole[Map[# == w &, f, {2}]], FourierParameters -> {0, -1}]]^2,
      {w, 0, r - 1}]]/nn^2;
```

#| collapse: input
```wl
SeedRandom[1];
samples = RandomChoice[Flatten[probs] -> Tuples[Range[0, nn - 1], 2], 20];
keyFrom[{c_, d_}] := With[{j = Mod[Round[r c/nn], r], l = Mod[Round[r d/nn], r]},
   If[CoprimeQ[j, r], Mod[l PowerMod[j, -1, r], r], Missing[]]];
cands = keyFrom /@ samples;
{cands, Select[DeleteMissing@DeleteDuplicates[cands], ec["Mult"][G, #] === Q &]}
```

Every usable shot returns the secret key, read from the interference pattern without scanning multiples of the generator.
