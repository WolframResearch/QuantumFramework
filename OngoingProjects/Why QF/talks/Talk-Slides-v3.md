<!-- All slides share one evaluation context. The initialization cell at the end (a "Code" cell, which Mathematica runs when the notebook opens) loads TalkSetup.wl; an evaluated md2nb build runs cells in document order, so it would need that cell first. -->

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

# Asking Physics Questions of the Wolfram Language

#| style: Subtitle

Eight questions answered by computation, with QuantumFramework

#| style: Subsubtitle

Mohammad Bahrami, Wolfram Research

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Eight questions

1. What does a precessing, decaying qubit do, for every rate at once? *The master equation solved symbolically.*
2. Does watching a qubit speed it up or freeze it? *A Liouvillian spectrum in closed form.*
3. Can a continuous error be fixed by a discrete correction? *Measurement records as circuit wires.*
4. Can watching a particle squeeze it below its zero-point motion? *A Riccati equation solved and proved for all parameters.*
5. What does a beam splitter do to $n$ photons, for every $n$? *Photon numbers kept as symbols.*
6. How far can exact simulation go without a state vector? *Stabilizer algebra.*
7. How short can a gate on a superconducting qubit be? *A five-level pulse model, evolved and scanned.*
8. Why are photon amplitudes expensive to compute? *Permanents and their cost.*

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 1. What does a precessing, decaying qubit do, for every rate at once?

**Idea:** solve the Lindblad equation once, with the precession frequency $\Omega$, the relaxation rate $\gamma_1$, the dephasing rate $\gamma_\phi$ and all four entries $\rho_{ij}$ of the initial state left as symbols. In the code `\[Sigma]m` is the lowering operator $|0\rangle\langle1|$ (the jump operator for relaxation) and `\[Rho]0` is the symbolic initial density matrix with entries `\[Rho][i, j]`.

```wl
Assuming[\[Gamma]1 > 0 && \[Gamma]\[Phi] >= 0 && \[CapitalOmega] \[Element] Reals && t >= 0,
 \[Rho]t = QuantumEvolve[\[CapitalOmega]/2 QuantumOperator["Z"], {Sqrt[\[Gamma]1] \[Sigma]m, Sqrt[\[Gamma]\[Phi]] QuantumOperator["Z"]}, \[Rho]0, t];
 FullSimplify /@ Normal[\[Rho]t["Computational"]["DensityMatrix"]]] // MatrixForm
```

The complete relaxation law of a qubit under amplitude damping and pure dephasing. Populations relax at $\gamma_1$ and never feel $\Omega$; coherences rotate at $\Omega$ and decay at $\gamma_1/2+2\gamma_\phi$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The same solution, drawn on the Bloch sphere

Take the Bloch vector of the symbolic solution, set every entry `\[Rho][i, j]` of the initial state to $1/2$ (the state $|+\rangle$), fix $\Omega=2\pi$, and pick three rate pairs $(\gamma_1,\gamma_\phi)$: relaxation only, relaxation with dephasing, and fast relaxation.

```wl
Show[ResourceFunction["BlochSpherePlot"][],
 ParametricPlot3D[Evaluate[Re[\[Rho]t["BlochVector"] /. \[Rho][_, _] -> 1/2 /. {\[CapitalOmega] -> 2 Pi, \[Gamma]1 -> #1, \[Gamma]\[Phi] -> #2}] & @@@ {{0.2, 0}, {0.2, 0.5}, {1.5, 0}}], {t, 0, 12},
  PlotStyle -> Thick, PlotLegends -> {"relaxation", "relaxation + dephasing", "fast relaxation"}], ImageSize -> 450]
```

Every trajectory spirals to the same pole, $|0\rangle$; dephasing only pulls the spiral inward faster, toward the axis.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### $T_2\le2T_1$ is a theorem, not a rule of thumb

With $T_1=1/\gamma_1$ and $T_2=1/(\gamma_1/2+2\gamma_\phi)$ read off the solution, ask `Reduce` whether the bound holds for every admissible rate.

```wl
Reduce[ForAll[{\[Gamma]1, \[Gamma]\[Phi]}, \[Gamma]1 > 0 && \[Gamma]\[Phi] >= 0, 1/(\[Gamma]1/2 + 2 \[Gamma]\[Phi]) <= 2/\[Gamma]1]]
```

It holds for all rates, with equality exactly when pure dephasing vanishes.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 2. Does watching a qubit speed it up or freeze it?

**Idea:** a qubit Rabi-driven at $\Omega$ while $\hat\sigma_z$ is measured at strength $k$ has a Liouvillian whose spectrum can be found exactly. In the code `liouvillian[H, {L}, d]` builds the Lindblad generator for Hamiltonian $H$ and jump operators $L$ on a $d$-level system, and `X`, `Z` are Pauli matrices.

```wl
ev = Eigenvalues@FullSimplify[liouvillian[(\[CapitalOmega]/2) X, {Sqrt[k] Z}, 2], \[CapitalOmega] > 0 && k > 0]
```

Two of the four rates, $-k\pm\sqrt{k^2-\Omega^2}$, carry the Rabi oscillation. They meet at $k=\Omega$, an exceptional point.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The two rates collide and separate

Follow the Rabi pair (`ev[[3 ;; 4]]`, the last two eigenvalues) in the complex plane as the measurement strength $k$ grows from $0$ to $3\Omega$, in units $\Omega=1$.

```wl
ListPlot[Transpose@Table[ReIm[ev[[3 ;; 4]] /. {\[CapitalOmega] -> 1, k -> kk}], {kk, 0, 3, 0.02}],
 AspectRatio -> 1, Frame -> True, FrameLabel -> {"Re \[Lambda] / \[CapitalOmega]  (minus the decay rate)", "Im \[Lambda] / \[CapitalOmega]  (frequency)"},
 PlotLegends -> {TraditionalForm[HoldForm[-k - Sqrt[k^2 - \[CapitalOmega]^2]]], TraditionalForm[HoldForm[-k + Sqrt[k^2 - \[CapitalOmega]^2]]]},
 Epilog -> {Red, PointSize[Large], Point[{-1, 0}]}, ImageSize -> 420]
```

Below $k=\Omega$ the pair moves left together, oscillating while decaying faster. At the red point the oscillation dies; beyond it one rate runs off to fast decay and the other returns toward zero.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Watching harder freezes the qubit

Expand the slower rate of the pair, $-k+\sqrt{k^2-\Omega^2}$, at strong measurement, $k\gg\Omega$.

```wl
Normal@Series[-k + Sqrt[k^2 - \[CapitalOmega]^2], {k, Infinity, 1}]
```

Now plot both features of the Rabi pair against $r=k/\Omega$, each as a positive number in units of $\Omega$: its damped Rabi frequency $\sqrt{\Omega^2-k^2}$, which exists only below the exceptional point, and the slowest decay rate, which is $k$ below it and $k-\sqrt{k^2-\Omega^2}$ above it.

```wl
Plot[{Re@Sqrt[1 - r^2], r - Re@Sqrt[r^2 - 1]}, {r, 0, 3}, Frame -> True, GridLines -> {{1}, None},
 PlotLegends -> {"damped Rabi frequency / \[CapitalOmega]", "slowest decay rate / \[CapitalOmega]"},
 FrameLabel -> {"k/\[CapitalOmega]", "frequency or rate / \[CapitalOmega]"}]
```

Measurement speeds relaxation only up to $k=\Omega$; beyond it the qubit freezes at the rate $\Omega^2/2k$, the quantum Zeno effect.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The same machinery gives the light an atom emits

Now the jump operator is spontaneous decay instead of a $\hat\sigma_z$ measurement, at rate $\gamma=1$, which sets the units of drive and detuning; `mollow[\[CapitalOmega], tones]` returns the inelastic spectrum at drive $\Omega$ on a list of detunings. Two-time correlations follow from this Liouvillian (quantum regression), so a fluorescence spectrum is one linear solve per frequency.

```wl
tones = Range[-9, 9, 0.05];
{weak, strong} = mollow[#, tones] & /@ {0.6, 6.};
ListLinePlot[{Transpose[{tones, weak/Max[weak]}], Transpose[{tones, strong/Max[strong]}]},
 Frame -> True, PlotRange -> All, PlotLegends -> {"\[CapitalOmega] = 0.6\[Gamma]", "\[CapitalOmega] = 6\[Gamma]"},
 FrameLabel -> {"detuning from the drive / \[Gamma]", "inelastic spectrum (scaled to peak)"}]
```

A strong drive dresses the atom and splits its single line into three, with sidebands at $\pm\Omega$: the Mollow triplet.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 3. Can a continuous error be fixed by a discrete correction?

**Idea:** in QuantumFramework a measurement result is a wire in the circuit. Encode one qubit in three, rotate the middle one by an arbitrary angle $2\theta$, measure two parities onto record wires, and correct with gates controlled by those records. The setup file builds this circuit as `qec`.

```wl
Show[qec["Diagram"], ImageSize -> 650]
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The rotation becomes a coin toss

Here `coded` is the output of `qec` acting on $(\alpha|0\rangle+\beta|1\rangle)|00\rangle$, and `asQEC` declares $\alpha,\beta,\theta$ real with $\alpha^2+\beta^2=1$. Tracing out the three data qubits leaves the joint distribution of the two parity records.

```wl
syndromes = Assuming[asQEC, FullSimplify[ComplexExpand[Diagonal[Normal[QuantumPartialTrace[coded, {3, 4, 5}]["DensityMatrix"]]]]]]
```

The four entries are the probabilities of the record pairs 00, 01, 10, 11. Plot the first and last against $\theta$: both parities even (no error) and both odd (qubit 2 flipped).

```wl
Plot[Evaluate@syndromes[[{1, 4}]], {\[Theta], 0, Pi}, Frame -> True, PlotLegends -> {"no error detected", "qubit 2 flipped"},
 FrameLabel -> {"\[Theta]  (qubit 2 rotated by 2\[Theta])", "probability"}, ImageSize -> 400]
```

Only two outcomes ever occur, with probabilities $\cos^2\theta$ and $\sin^2\theta$, independent of $\alpha$ and $\beta$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### And the correction is exact

Tracing out the two records leaves the corrected data qubits; compare them with the ideal encoded state, `enc` (the encoding circuit alone) applied to the same input.

```wl
Assuming[asQEC, FullSimplify[ComplexExpand[QuantumSimilarity[QuantumPartialTrace[coded, {1, 2}],
   enc[QuantumTensorProduct[QuantumState[{\[Alpha], \[Beta]}], QuantumState["00"]]], "Fidelity"]]]]
```

Fidelity 1 for every real-amplitude input and every angle. Measurement turns a continuous error into a discrete one, which is why quantum error correction works at all.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 4. Can watching a particle squeeze it below its zero-point motion?

**Idea:** for a particle in a harmonic trap whose position is measured continuously at strength $\lambda$, the uncertainties of the conditional state obey a Riccati equation. Its steady state comes out with $\lambda,\omega,m,\hbar$ symbolic: `SxSS`, `SpSS` and `CcSS` are the steady position variance $\Sigma_{xx}$, momentum variance $\Sigma_{pp}$ and position-momentum covariance $\Sigma_{xp}$; `asTrap` says all four parameters are positive; `sqRatio[\[Mu]]` is the closed-form ratio of the steady position variance to the zero-point variance $\hbar/2m\omega$.

```wl
{Simplify[SxSS SpSS - CcSS^2, asTrap], Simplify[SxSS/(\[HBar]/(2 m \[Omega])) == sqRatio[2 \[HBar] \[Lambda]/(m \[Omega]^2)], asTrap]}
```

The steady state is pure (its covariance determinant $\Sigma_{xx}\Sigma_{pp}-\Sigma_{xp}^2$ takes the minimum value $\hbar^2/4$), and its position spread depends only on one number, $\mu=2\hbar\lambda/m\omega^2$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Squeezed at every measurement strength

`Resolve` decides whether the position variance is below the zero-point value for every $\mu>0$.

```wl
Resolve[ForAll[\[Mu], \[Mu] > 0, sqRatio[\[Mu]] < 1], Reals]
```

The exact ratio against $\mu$, with its weak-measurement limit $1-\mu^2/8$ and strong-measurement limit $\sqrt{2/\mu}$ dashed.

```wl
Plot[{sqRatio[\[Mu]], 1 - \[Mu]^2/8, Sqrt[2/\[Mu]]}, {\[Mu], 0, 12}, PlotRange -> {0, 1.1}, Frame -> True,
 PlotStyle -> {Thick, Dashed, Dashed}, PlotLegends -> {"exact", Row[{"weak limit  ", TraditionalForm[HoldForm[1 - \[Mu]^2/8]]}], Row[{"strong limit  ", TraditionalForm[HoldForm[Sqrt[2/\[Mu]]]]}]},
 FrameLabel -> {"measurement strength \[Mu]", "position variance / zero-point variance"}, ImageSize -> 450]
```

Watching always squeezes the position below the zero-point value: slowly for weak measurement, as $\sqrt{2/\mu}$ for strong.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What an imperfect detector costs

A detector that catches a fraction $\eta$ of the signal weakens the information gain but not the measurement kick. In the steady solution this is the substitution $\lambda\to\eta\lambda$, $\hbar\to\hbar/\sqrt\eta$, after which the purity of the Gaussian state is $\hbar/2\sqrt{\Sigma_{xx}\Sigma_{pp}-\Sigma_{xp}^2}$.

```wl
With[{sol = {SxSS, CcSS, SpSS} /. {\[Lambda] -> \[Eta] \[Lambda], \[HBar] -> \[HBar]/Sqrt[\[Eta]]}},
 Simplify[\[HBar]/(2 Sqrt[sol[[1]] sol[[3]] - sol[[2]]^2]), {\[HBar] > 0, 0 < \[Eta] <= 1, \[Lambda] > 0, \[Omega] > 0, m > 0}]]
```

The steady purity is exactly $\sqrt\eta$, whatever the trap, the mass or the measurement strength.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The watched particle, reading by reading

The full wavefunction on a 256-point grid, updated by Bayes' rule after every noisy position reading and evolved in the trap between readings ($\hbar=m=\omega=1$, $\lambda=2$). `states` holds the wavefunction after every step on the grid `gH`, and `ntH dtH` is the total time; the plot shows every sixth step for $|x|<3$.

```wl
With[{keep = Flatten@Position[gH["x"], x_ /; Abs[x] < 3]},
 ListDensityPlot[Transpose[Abs[states[[;; ;; 6, keep]]]^2], DataRange -> {{0, ntH dtH}, {-3, 3}}, PlotRange -> All,
  ColorFunction -> "SunsetColors", AspectRatio -> 1/2.5, FrameLabel -> {"time", "position x"}, ImageSize -> 600]]
```

The packet stays narrow while its center wanders with the readings: measurement localizes it and kicks it at the same time.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Five numbers do the job of 256 amplitudes

A quantum Kalman filter (two means and three second moments, fed only the readings) against the mean of the full wavefunction: `xGrid` is $\langle x\rangle$ from the grid, `mFilter` the filter's estimate.

```wl
ListLinePlot[{xGrid, mFilter}, DataRange -> {0, ntH dtH}, PlotStyle -> {Directive[Thick, Gray], Directive[Red, Dashed]},
 PlotLegends -> {"full wavefunction (256 amplitudes)", "Kalman filter (5 numbers)"}, Frame -> True,
 FrameLabel -> {"time", "\[LeftAngleBracket]x\[RightAngleBracket]"}, ImageSize -> 550]
```

The two curves agree up to the finite time step. In a harmonic trap the conditional state stays Gaussian, so in continuous time the filter is exact, not an approximation.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 5. What does a beam splitter do to $n$ photons, for every $n$?

**Idea:** `QuantumInterferometer` keeps the photon numbers, angles and squeezing as symbols, and draws itself. Start with one beam splitter whose mixing angle $\theta$ and phase $\phi$ are symbols.

```wl
bs = QuantumInterferometer[{"BS"[\[Theta], \[Phi]]}]
```

Its diagram: two modes entering one splitter.

```wl
Show[bs["Diagram"], ImageSize -> 200]
```

The amplitude from $|n_1,n_2\rangle$ to $|t_1,n_1+n_2-t_1\rangle$, with the photon numbers $n_1,n_2,t_1$ also symbols (the splitter conserves the total).

```wl
bs["Amplitude", {n1, n2} -> {t1, n1 + n2 - t1}]
```

It is a Wigner $D$-function: a beam splitter rotates a spin of size $(n_1+n_2)/2$, and the photon numbers are its magnetic quantum numbers.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Hong-Ou-Mandel for every $n$

$n$ photons in each port of a balanced splitter, with $k$ photons leaving in mode 1: the amplitude as one formula in $n$ and $k$.

```wl
hom = QuantumInterferometer[{"BS"[Pi/4, 0]}]["Amplitude", {n, n} -> {k, 2 n - k}]
```

Its square gives the output distribution; draw it for $n=1$ to $4$, one bar for each $k$.

```wl
Row[Table[BarChart[Table[Abs[hom /. {n -> m, k -> j}]^2, {j, 0, 2 m}], ChartLabels -> Range[0, 2 m], Frame -> True,
   FrameLabel -> {"photons in mode 1", "probability"}, PlotLabel -> Row[{"n = ", m}], ImageSize -> 220], {m, 4}], Spacer[10]]
```

Every odd split is exactly zero. This holds for every $n$, not only those shown: the amplitude $d^{\,n}_{0,k-n}(\pi/2)$ is proportional to the Legendre function $P_n^{k-n}(0)$, which vanishes whenever $k$ is odd. $n=1$ is the familiar Hong-Ou-Mandel dip.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Four photons need an unbalanced splitter to cancel

The amplitude `fa` for $|2,2\rangle$ to come out as $|2,2\rangle$ through a splitter of angle $\theta$, checked against $1-6R(1-R)$ with reflectivity $R=\sin^2\theta$, and the reflectivities where it vanishes.

```wl
fa = FullSimplify[QuantumInterferometer[{"BS"[\[Theta]]}]["Amplitude", {2, 2} -> {2, 2}], \[Theta] \[Element] Reals];
{FullSimplify[fa == 1 - 6 Sin[\[Theta]]^2 Cos[\[Theta]]^2, \[Theta] \[Element] Reals], R /. Solve[1 - 6 R (1 - R) == 0 && 0 < R < 1, R, Reals]}
```

The return probability against $R$, with the two zeros marked in red.

```wl
Plot[(1 - 6 R (1 - R))^2, {R, 0, 1}, Frame -> True, FrameLabel -> {"reflectivity R", "probability |2,2\[RightAngleBracket] \[Rule] |2,2\[RightAngleBracket]"},
 Epilog -> {Red, PointSize[Large], Point[{(3 - Sqrt[3])/6, 0}], Point[{(3 + Sqrt[3])/6, 0}]}, ImageSize -> 400]
```

The amplitude is $1-6R(1-R)$, zero at $R=(3\pm\sqrt3)/6$ and not at the balanced point: the cancelling splitter changes with photon number.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A sign flip no phase shifter can make

The KLM nonlinear sign gate: `nsGate` is the three-mode network of a fixed real $3\times3$ unitary from the setup file, drawn here.

```wl
Show[nsGate["Diagram"], ImageSize -> 300]
```

The input $\{n,1,0\}$ puts $n$ signal photons in mode 1 and one helper photon in mode 2; success is heralded by finding the helper pattern $\{1,0\}$ again. The heralded amplitude for $n=0$ to $3$:

```wl
RootReduce[Table[nsGate["Amplitude", {n, 1, 0} -> {n, 1, 0}], {n, 0, 3}]]
```

$0$, $1$ and $2$ signal photons get amplitudes $\tfrac12,\tfrac12,-\tfrac12$: a sign on $|2\rangle$ alone, with success probability $1/4$ for every input. The three-photon amplitude differs, which is where the gate stops working.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Squeezed light comes in pairs

Squeeze mode 1 with parameter $r$ and mode 2 with $-r$, then mix them on a balanced splitter.

```wl
pair = QuantumInterferometer[{"S"[r], "S"[-r] -> 2, "BS"[Pi/4]}]
```

The network drawn: two squeezers followed by the splitter.

```wl
Show[pair["Diagram"], ImageSize -> 220]
```

Starting from vacuum, the probability of $j$ photons in each output mode, for $j=0,1,2$.

```wl
FullSimplify[Table[pair["Probability", {0, 0} -> {j, j}], {j, 0, 2}], r > 0]
```

Two checks: the probability of an unequal outcome, $|2,0\rangle$, and the pair law $\operatorname{sech}^2r\tanh^{2j}r$ summed over all $j$.

```wl
FullSimplify[{pair["Probability", {0, 0} -> {2, 0}], Sum[Sech[r]^2 Tanh[r]^(2 j), {j, 0, Infinity}]}, r > 0]
```

Pair probabilities follow $\operatorname{sech}^2r\tanh^{2j}r$; unequal outcomes vanish and the pair law sums to one, so the output modes always hold equal photon numbers: a two-mode squeezed vacuum.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A cat's interference fringes

`CatState[][2, 0]` is the even superposition of the coherent states $|\pm2\rangle$. Tabulate its Wigner function on the phase-space window $|x|\le6$, $|p|\le2$.

```wl
wigner = WignerRepresentation[CatState[][2, 0], {-6, 6}, {-2, 2}, "GridSize" -> 90];
```

Tabulate its Husimi function, the Wigner function smoothed over a vacuum-sized patch, on $|x|\le6$, $|p|\le3$.

```wl
husimi = HusimiQRepresentation[CatState[][2, 0], {-6, 6}, {-3, 3}, "GridSize" -> 90];
```

The Wigner function as a surface over $(x,p)$.

```wl
Plot3D[wigner[x, p], {x, -6, 6}, {p, -2, 2}, PlotRange -> All, Mesh -> None, PlotPoints -> 50, ColorFunction -> "SunsetColors", AxesLabel -> {"x", "p", None}, PlotLabel -> "Wigner", ImageSize -> 300]
```

The Husimi function on the same axes.

```wl
Plot3D[husimi[x, p], {x, -6, 6}, {p, -3, 3}, PlotRange -> All, Mesh -> None, PlotPoints -> 50, ColorFunction -> "SunsetColors", AxesLabel -> {"x", "p", None}, PlotLabel -> "Husimi", ImageSize -> 300]
```

The Wigner function has one lobe for each coherent state and fringes between them that dip below zero, the signature of the superposition; the Husimi function keeps the lobes and smooths the fringes away.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 6. How far can exact simulation go without a state vector?

**Idea:** a stabilizer state is named by the Pauli operators that fix it, and Clifford gates rewrite that list. Start with the five-qubit code: which of its four checks does each single-qubit error flip? `errors` lists the fifteen single-qubit Pauli errors, and `syndrome[e]` returns the four check outcomes (1 = flipped) for error `e`.

```wl
table = syndrome /@ errors;
```

Draw the table: one column per error, one row per check, blue where the check flips.

```wl
ArrayPlot[Transpose[table], Mesh -> True, ColorRules -> {0 -> White, 1 -> Darker[Blue]}, Frame -> True,
 FrameTicks -> {{Thread[{Range[4], {"check 1", "check 2", "check 3", "check 4"}}], None}, {Thread[{Range[15], #[[1, 1]] <> ToString[#[[1, 2]]] & /@ errors}], None}}, ImageSize -> 600]
```

Check that the fifteen patterns are exactly the fifteen nonzero four-bit strings.

```wl
Sort[table] === Sort@DeleteCases[Tuples[{0, 1}, 4], {0, 0, 0, 0}]
```

The fifteen errors flip fifteen different nonzero patterns, using every pattern exactly once: what makes the code "perfect".

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A thousand qubits, exactly

Build a 1000-qubit GHZ state as a stabilizer tableau: a Hadamard on qubit 1, then a chain of CNOTs.

```wl
ghz = PauliStabilizer[1000][Prepend[Table["CNOT" -> {q - 1, q}, {q, 2, 1000}], "H"]];
```

The entanglement entropy of qubit 1 with the rest, and the time it takes.

```wl
AbsoluteTiming[ghz["Entropy", {1}]]
```

The entanglement entropy of the first 500 qubits with the other 500.

```wl
ghz["Entropy", Range[500]]
```

One bit of entanglement for either cut, computed in a fraction of a second: GHZ entanglement is a single shared parity, however the chain is cut.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The price of non-Clifford gates

Three T gates on a GHZ backbone of growing width: `magicState[n]` is an $n$-qubit GHZ circuit followed by T on qubits 1 to 3, and its `"Length"` is the number of terms in its stabilizer frame. Compare with the $2^n$ amplitudes of the state vector.

```wl
ListLogPlot[{Table[{n, magicState[n]["Length"]}, {n, 4, 24, 2}], Table[{n, 2^n}, {n, 4, 24, 2}]}, Joined -> True,
 PlotMarkers -> Automatic, Frame -> True, FrameLabel -> {"qubits", "terms or amplitudes"},
 PlotLegends -> {"stabilizer frame terms (3 T gates)", "state-vector amplitudes"}, ImageSize -> 420]
```

One T gate on $|+\rangle=H|0\rangle$: the probability of finding $|+\rangle$ is $(1+\langle X\rangle)/2$, with $\langle X\rangle$ from the stabilizer frame.

```wl
FullSimplify[(1 + PauliStabilizer[1][{"H", "T"}]["Expectation", "X"])/2]
```

At a fixed number of T gates the stabilizer frame keeps the same number of terms as qubits are added (each term, a tableau, grows only polynomially), while the state vector doubles with each qubit; and after $TH|0\rangle$ the probability of finding $|+\rangle$ comes out as the irrational $(2+\sqrt2)/4$, exactly.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 7. How short can a gate on a superconducting qubit be?

**Idea:** a transmon is an oscillator with a small anharmonicity. Keep five levels, drive an 8 ns π pulse, and add a DRAG correction proportional to the pulse's time derivative; everything is one parametrized `QuantumOperator`. `pulse[levels, {A, \[Lambda], \[Delta], \[Tau]}]` evolves the transmon from $|0\rangle$ under a pulse of amplitude $A$ ($A=1$ is a π pulse), DRAG coefficient $\lambda$, detuning $\delta$ and length $\tau$ in ns; `s[t]` is the state at time $t$.

```wl
Row[Table[With[{s = pulse[5, {1, drag, 0, 8}]},
   ListLinePlot[Transpose@Table[s[tt]["ProbabilitiesList"][[;; 3]], {tt, 0, 8, 0.1}], DataRange -> {0, 8}, PlotRange -> {0, 1},
    Frame -> True, FrameLabel -> {"time (ns)", "population"}, PlotLabel -> If[drag == 0, "plain pulse", "with DRAG (\[Lambda] = 1)"],
    PlotLegends -> {"|0\[RightAngleBracket]", "|1\[RightAngleBracket]", "|2\[RightAngleBracket]"}, ImageSize -> 330]], {drag, {0, 1}}], Spacer[10]]
```

The plain pulse visibly pushes population into $|2\rangle$, outside the qubit; DRAG suppresses that leakage during the pulse.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Suppressing leakage is not the same as completing the flip

Scan the DRAG coefficient at 8 ns and track both errors. `finalPopulations` returns the level populations at the end of the pulse: leakage is the population above $|1\rangle$, transfer error is $1-P_{|1\rangle}$.

```wl
ListLogPlot[Transpose@Table[With[{p = finalPopulations[5, {1, drag, 0, 8}]},
    {{drag, Total[p[[3 ;;]]]}, {drag, 1 - p[[2]]}}], {drag, 0, 2, 0.1}],
 Joined -> True, PlotMarkers -> Automatic, PlotLegends -> {"leakage", "transfer error"},
 Frame -> True, FrameLabel -> {"DRAG coefficient \[Lambda]", "probability"}, ImageSize -> 420]
```

Leakage is smallest near the weak-drive value $\lambda\approx1$, but the flip is most complete well below it, so amplitude, detuning and DRAG must be calibrated together. These runs are closed-system; the tutorial's full scan under the Lindblad equation finds that the speed limit is then set by $T_1$ and $T_2$, not by leakage.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 8. Why are photon amplitudes expensive to compute?

**Idea:** the probability of an output pattern is a permanent for bosons, a permanent of $|W|^2$ for distinguishable particles, and a determinant for fermions. Compare the three on one balanced beam splitter. `bosonProb`, `distProb` and `fermionProb[W, s, t]` give the probability of output pattern $t$ from input pattern $s$ through transfer matrix $W$; `Wbs` is the balanced splitter's transfer matrix.

```wl
BarChart[Transpose@Table[Through[{bosonProb, distProb, fermionProb}[Wbs, {1, 1}, out]], {out, {{2, 0}, {1, 1}, {0, 2}}}],
 ChartLabels -> {{"bosons", "classical particles", "fermions"}, {"|2,0\[RightAngleBracket]", "|1,1\[RightAngleBracket]", "|0,2\[RightAngleBracket]"}}, Frame -> True,
 FrameLabel -> {None, "probability"}, ImageSize -> 500]
```

Bosons always leave together, classical particles split half the time, fermions never leave together. Bosons and fermions differ only in the sign of the exchange term; distinguishable particles have no exchange term at all.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What one amplitude costs

Time to compute one output amplitude in a random $2n$-mode network, for three ways of loading $n$ photons. `occupations[n]` gives the three input patterns, `randomInterferometer[n]` a fixed random network, and `timeAmplitude[q, s]` the seconds one amplitude takes from a cleared cache.

```wl
timings = AssociationMap[Function[n, With[{q = randomInterferometer[n]}, timeAmplitude[q, #] & /@ occupations[n]]], Range[4, 22, 2]];
ListLogPlot[Table[Transpose[{Keys[timings], Lookup[Values[timings], kind]}], {kind, Keys[First[timings]]}], Joined -> True, PlotMarkers -> Automatic,
 PlotLegends -> Keys[First[timings]], Frame -> True, FrameLabel -> {"photons n", "seconds per amplitude"}, ImageSize -> 450]
```

One photon per mode and photons in pairs cost the same, and that cost grows exponentially with every photon; all photons in one mode stays cheap. Computing a probability is not yet sampling: distinguishable photons also have permanent probabilities, yet are easy to sample. The hardness argument of Aaronson and Arkhipov works where modes far outnumber photons, so that almost every output is of the expensive, collision-free kind, and turns that cost into sampling hardness through complexity theory.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## What answered each question

1. **Decaying qubit:** the master equation solved once with every rate symbolic, and a bound proved by `Reduce`.
2. **Watched qubit:** an exact Liouvillian spectrum, with the Zeno freeze read off its exceptional point.
3. **Continuous error:** measurement records as circuit wires, giving an exact symbolic correction.
4. **Watched particle:** a Riccati equation solved and proved for all parameters, and a Kalman filter built on it that tracks a full wavefunction simulation.
5. **Beam splitter:** photon numbers kept as symbols, so laws hold for every $n$, and every network draws itself.
6. **Without a state vector:** stabilizer algebra, exact at a thousand qubits.
7. **Fast gate:** one parametrized pulse model, evolved and scanned.
8. **Photon amplitudes:** permanents and determinants for the three particle statistics, and the measured cost of permanent amplitudes.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Initialization

All definitions behind the slides are in one file next to this notebook.

```wl
#| style: Code
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Why QF/talks/TalkSetup.wl"]
```
