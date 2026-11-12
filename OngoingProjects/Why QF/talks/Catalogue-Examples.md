---
Template: Default
Title: "What Wolfram Language and QuantumFramework Compute Exactly: Worked Examples"
Author: Mads Bahrami
---

# What Wolfram Language and QuantumFramework Compute Exactly: Worked Examples

Each example below asks one physics question with the parameters left as symbols, and the answer comes back as a formula, not a table of numbers. Every cell runs as written on QuantumFramework at commit ed916377.

Load the framework and its optics context:

```wl
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"]
```
<!-- run time: 2.6 s, mostly loading -->

## A QAOA Layer on MaxCut: the Exact Optimum of the Energy Landscape

One layer of the quantum approximate optimization algorithm prepares $|+\rangle^{\otimes n}$, applies $e^{-i\gamma Z_iZ_j/2}$ on every edge and $e^{-i\theta X/2}$ on every vertex, and scores the result by the expected cut $\langle H_C\rangle = \sum_{(i,j)}\langle(1 - Z_iZ_j)/2\rangle$. How good is the best choice of $(\gamma, \theta)$? A circuit simulator evaluates the landscape at sampled angles and hands the result to a numerical optimizer. Here the angles stay symbols, so the landscape is a formula and its maximum is a number known exactly.

Build $\langle\psi(\gamma,\theta)|H_C|\psi(\gamma,\theta)\rangle$ for any graph, straight from its edge list:

```wl
qaoa[graph_] := With[{e = List @@@ EdgeList[graph], n = VertexCount[graph]},
  With[{ψ = QuantumCircuitOperator[{"H" -> Range[n], Splice[("R"[γ, "ZZ"] -> #) & /@ e], "RX"[θ] -> Range[n]}][]},
   FullSimplify[ψ["Dagger"][Total[(QuantumOperator[StringRepeat["I", n]] - QuantumOperator["ZZ", #])/2 & /@ e][ψ]]["Scalar"], {γ, θ} ∈ Reals]]];
qaoa[CycleGraph[3]]
```
<!-- run time: 0.79 s -->

The triangle's landscape is a trigonometric polynomial in the two angles. Maximize it exactly:

```wl
Maximize[qaoa[CycleGraph[3]], {γ, θ}]
```
<!-- run time: 0.83 s -->

The maximum is $2$, the classical maximum cut of a triangle (no two-coloring cuts all three edges of an odd cycle). In other words, one layer already finds the best cut with certainty here, and the optimal angles come back as exact closed forms, not as the output of a converged search.

The complete graph $K_4$ is harder: its classical maximum cut is $4$ of $6$ edges. Ask for the one-layer optimum:

```wl
First @ Maximize[qaoa[CompleteGraph[4]], {γ, θ}]
```
<!-- run time: 2.21 s -->

The optimum is a root of a cubic, and it lies below $4$. In other words, the exact landscape certifies that one layer cannot reach the true cut of $K_4$; a sampled optimizer could only report that it had not found a better point.

Check against the trivial circuit: at $\gamma = 0$ the cost layer does nothing, the state is an equal superposition of all colorings, and each of the six edges is cut with probability $1/2$:

```wl
Simplify[qaoa[CompleteGraph[4]] /. γ -> 0]
```
<!-- run time: 1.48 s -->

## Two Exchange-Coupled Spins: the Entanglement Entropy Is the Binary Entropy of $\cos^2 2Jt$

Two spins coupled by the Heisenberg exchange $H = J(XX + YY + ZZ)$ start in the product state $|01\rangle$. How entangled are they at time $t$? A numerical solver returns the entropy on a time grid for one value of $J$. Here $J$ and $t$ stay symbols, and the entropy comes back as a function.

Evolve $|01\rangle$ and take the reduced state of the first spin:

```wl
$Assumptions = J ∈ Reals && t ∈ Reals;
heis = QuantumEvolve[J (QuantumOperator["XX"] + QuantumOperator["YY"] + QuantumOperator["ZZ"]), QuantumState["01"], t];
FullSimplify[Normal[QuantumPartialTrace[heis, {2}]["DensityMatrix"]]]
```
<!-- run time: 0.22 s -->

The reduced state is diagonal with populations $\cos^2 2Jt$ and $\sin^2 2Jt$. The excitation swaps between the two spins at angular frequency $4J$, the singlet-triplet splitting, because $|01\rangle$ is an equal superposition of the singlet and the $m = 0$ triplet.

Compute $S(t) = -\operatorname{Tr}\rho_1\log_2\rho_1$:

```wl
sJ = FullSimplify[QuantumPartialTrace[heis, {2}]["VonNeumannEntropy"]]
```
<!-- run time: 0.61 s -->

This is $h_2(\cos^2 2Jt) = -\cos^2 2Jt\,\log_2\cos^2 2Jt - \sin^2 2Jt\,\log_2\sin^2 2Jt$, rewritten with $\csc 4Jt$ and $\tan 2Jt$. In other words, the spins become maximally entangled whenever the excitation is shared equally and return to a product state whenever it sits on one spin.

At $Jt = \pi/8$ the populations are equal, so the entropy should be exactly one bit:

```wl
sJ /. t -> π/(8 J)
```
<!-- run time: < 0.001 s -->

## Qubit Relaxation for Every Initial State and Both Rates: $T_2 \le 2T_1$ as a Theorem

A qubit with splitting $\omega$ relaxes at rate $\gamma_1$ (jump operator $\sigma_- = |0\rangle\langle 1|$) and dephases at rate $\gamma_\phi$ (jump operator $Z$). What is $\rho(t)$ for every initial state, and what bound ties the coherence time to the relaxation time? A numerical solver integrates one initial state at one pair of rates. Here the initial density matrix has four symbolic entries and both rates are symbols.

Solve the Lindblad equation from $\rho_0 = (r_{ij})$:

```wl
$Assumptions = γ1 > 0 && γφ >= 0 && ω ∈ Reals && t >= 0;
ρt = QuantumEvolve[ω/2 QuantumOperator["Z"], {QuantumOperator[{{0, 1}, {0, 0}}], QuantumOperator["Z"]} -> {γ1, γφ},
   QuantumState[Array[r, {2, 2}]], t];
MatrixForm[ρm = FullSimplify[Normal[ρt["DensityMatrix"]]]]
```
<!-- run time: 0.12 s -->

The excited population decays at $\gamma_1$ into the ground state, and the coherence precesses at $\omega$ while it decays at $\gamma_1/2 + 2\gamma_\phi$. So $T_1 = 1/\gamma_1$, and $T_2$ is read from the same matrix that covers every qubit experiment.

Read the coherence decay rate $1/T_2$ off the computed off-diagonal entry:

```wl
γT2 = FullSimplify[-Re[D[Log[ρm[[1, 2]]], t]]]
```
<!-- run time: 0.03 s -->

Ask whether $T_2 \le 2T_1$ holds for all admissible rates:

```wl
Reduce[ForAll[{γ1, γφ}, γ1 > 0 && γφ >= 0, 1/γT2 <= 2/γ1]]
```
<!-- run time: 0.001 s -->

Find where the bound is saturated:

```wl
Reduce[1/γT2 == 2/γ1 && γ1 > 0 && γφ >= 0, γφ]
```
<!-- run time: < 0.001 s -->

The bound is tight only without pure dephasing. In other words, $T_2 = 2T_1$ is the limit set by energy relaxation alone, and any dephasing pushes $T_2$ below it.

At long times every initial state should relax to $|0\rangle\langle 0|$, carrying its full trace there:

```wl
Limit[ρm, t -> Infinity, Assumptions -> γ1 > 0 && γφ >= 0 && ω ∈ Reals]
```
<!-- run time: 0.09 s -->

## Quantum Fisher Information of a Phase: $n$ for Product Probes, $n^2$ for GHZ

A phase $\phi$ is written on $n$ qubits by $R_Z(\phi)$ on each. The quantum Fisher information $F_Q$ bounds every estimator through $\operatorname{Var}\hat\phi \ge 1/F_Q$, and for a pure state it is four times the Fubini-Study metric, $F_Q = 4\left(\langle\partial_\phi\psi|\partial_\phi\psi\rangle - |\langle\psi|\partial_\phi\psi\rangle|^2\right)$. How does it grow with $n$ for a product probe and for a GHZ probe? A numerical tool estimates it by finite differences at one value of $\phi$. Here $\phi$ stays a symbol and $F_Q$ comes out exact.

Define $F_Q$ of the phase family generated from a probe state (`FubiniStudyMetricTensor` lives in the optimization context), and compare $|+\rangle^{\otimes n}$ with $(|0\cdots0\rangle + |1\cdots1\rangle)/\sqrt2$ for $n = 1$ to $5$:

```wl
Needs["Wolfram`QuantumFramework`QuantumOptimization`"];
qfi[st_] := Simplify[4 FubiniStudyMetricTensor[
     QuantumState[QuantumCircuitOperator[Table["RZ"[φ] -> k, {k, st["Qudits"]}]][st], "Parameters" -> {φ}], "Matrix"][[1, 1]], φ ∈ Reals];
Table[{n, qfi[QuantumState[StringRepeat["+", n]]], qfi[QuantumState["GHZ"[n]]]}, {n, 5}]
```
<!-- run time: 0.83 s -->

The product probe gives $F_Q = n$, the standard quantum limit, and the GHZ probe gives $n^2$, the Heisenberg limit. Both are free of $\phi$, because the generator $\tfrac12\sum_i Z_i$ commutes with the evolution it generates. In other words, $F_Q = 4\operatorname{Var}(\tfrac12\sum_i Z_i)$: each qubit of the product state contributes its own unit of variance, while the GHZ state puts all its weight on the two extreme eigenvalues $\pm n/2$, so the variance grows as $n^2$.

Cross-check with the variance formula directly, for the four-qubit GHZ state:

```wl
With[{v = Normal[QuantumState["GHZ"[4]]["StateVector"]], g = DiagonalMatrix[Total /@ (1/2 - Tuples[{0, 1}, 4])]},
 4 (v . g . g . v - (v . g . v)^2)]
```
<!-- run time: 0.001 s -->

## Coherent Information of Amplitude Damping: $h_2((1-\gamma)p) - h_2(\gamma p)$, Maximized Exactly

How much quantum information survives amplitude damping with decay probability $\gamma$? Purify the input $\rho_p = (1-p)|0\rangle\langle0| + p|1\rangle\langle1|$ with a reference qubit, send only the system through the channel, and compute the coherent information $I_c = S(B) - S(RB)$. A numerical tool evaluates this at chosen $(p, \gamma)$ and scans for the best input. Here both stay symbols.

Purify, apply the channel to the first qubit, and take the two entropies:

```wl
out = QuantumChannel["AmplitudeDamping"[g], {1}][QuantumState[DiagonalMatrix[{1 - p, p}]]["Purify"]];
ic = FullSimplify[QuantityMagnitude[QuantumPartialTrace[out, {2}]["VonNeumannEntropy"] - out["VonNeumannEntropy"]],
   0 < p < 1 && 0 < g < 1]
```
<!-- run time: 0.97 s -->

Check that this is the difference of two binary entropies, $h_2((1-\gamma)p) - h_2(\gamma p)$:

```wl
h2[x_] := -x Log[2, x] - (1 - x) Log[2, 1 - x];
FullSimplify[ic == h2[(1 - g) p] - h2[g p], 0 < p < 1 && 0 < g < 1]
```
<!-- run time: 1.81 s -->

The output keeps the surviving population $(1-\gamma)p$ and the environment receives the leaked population $\gamma p$; the coherent information is their competition. At $\gamma = 1/2$ the two are equal and $I_c$ vanishes for every input:

```wl
FullSimplify[ic /. g -> 1/2, 0 < p < 1]
```
<!-- run time: 0.01 s -->

Maximize over the input at $\gamma = 1/4$:

```wl
FullSimplify[Maximize[{ic /. g -> 1/4, 0 <= p <= 1}, p]]
```
<!-- run time: 0.20 s -->

The best input carries $\log_2(4/3)$ bits at $p = 4/9$, not at the maximally mixed $p = 1/2$. In other words, the optimal input underweights the excited state, because every excitation sent is also the only thing the environment can learn from.

## The Data-Processing Inequality Under Named Channels: `Reduce` Finds No Violating Noise Strength

The quantum relative entropy $S(\rho\|\sigma) = \operatorname{Tr}\rho(\log_2\rho - \log_2\sigma)$ measures how well $\rho$ can be told from $\sigma$, and no channel may raise it: $S(\mathcal N(\rho)\|\mathcal N(\sigma)) \le S(\rho\|\sigma)$. Take two qubit states with non-parallel Bloch vectors, so they do not commute. Does the inequality hold for every noise strength $p$ of a channel? A numerical tool checks it at sampled $p$. Here $p$ is a symbol, and `Reduce` searches the whole interval $0 < p < 1$ for a violation.

Set up the pair and the relative entropy, with the domain of $p$ supplied so the closed form simplifies:

```wl
{qr, qs} = QuantumState["BlochVector"[#]] & /@ {{0, 0, 3/5}, {2/5 Sqrt[3], 0, 2/5}};
rel[a_, b_] := FullSimplify[QuantityMagnitude @ QuantumDistance[a, b, "RelativeEntropy"], 0 < p < 1];
rel[QuantumChannel["Depolarizing"[p]][qr], QuantumChannel["Depolarizing"[p]][qs]]
```
<!-- run time: 1.02 s -->

After depolarizing noise the relative entropy is an explicit function of $p$, built from logarithms of the shrunken eigenvalues. Ask `Reduce` for any $p$ that raises it, under four named channels:

```wl
Table[c -> Reduce[rel[QuantumChannel[c[p]][qr], QuantumChannel[c[p]][qs]] > rel[qr, qs] && 0 < p < 1, p, Reals],
 {c, {"Depolarizing", "PhaseDamping", "BitFlip", "PhaseFlip"}}]
```
<!-- run time: 9.42 s -->

Every answer is `False`: no noise strength in $(0, 1)$ makes the processed pair easier to distinguish. In other words, each channel is a proof for this pair, not a sample of it.

Check that the named distance is $S(\rho\|\sigma)$ with the arguments in this order, against the matrix-logarithm definition:

```wl
With[{a = Normal[qr["DensityMatrix"]], b = Normal[qs["DensityMatrix"]]},
 N[{Tr[a . (MatrixLog[a] - MatrixLog[b])]/Log[2], rel[qr, qs]}]]
```
<!-- run time: 0.01 s -->

## SU(2) Double Cover from a Controlled Rotation: the Fringe Is $\sin^2(\alpha/4)$

A rotation by $2\pi$ returns every Bloch vector but multiplies the spinor by $-1$, so no measurement on the rotated qubit alone can see it. Put an ancilla in superposition, rotate a target by $R_X(\alpha)$ only in the ancilla's $|1\rangle$ branch, and close the interferometer with a second Hadamard. What is the probability of finding the ancilla in $|1\rangle$, as a function of $\alpha$? A circuit simulator returns it at each sampled angle. Here $\alpha$ is a symbol.

Build the circuit and read all four outcome probabilities:

```wl
crx = QuantumOperator["Controlled"["RX"[α] -> {2}, {1}]];
probs = FullSimplify[#, α ∈ Reals] & /@ QuantumCircuitOperator[{"H" -> 1, crx, "H" -> 1, {1}, {2}}][]["Probabilities"]
```
<!-- run time: 0.28 s -->

Add the two outcomes with the ancilla in $|1\rangle$:

```wl
fringe = Simplify[probs[[3]] + probs[[4]]]
```
<!-- run time: < 0.001 s -->

The fringe is $\sin^2(\alpha/4)$, with period $4\pi$ in the rotation angle. In other words, the ancilla reads out $(1 - \operatorname{Re}\langle 0|R_X(\alpha)|0\rangle)/2 = (1 - \cos(\alpha/2))/2$, the overlap of the rotated and unrotated branches, and that overlap is a spinor amplitude, which needs $4\pi$ to return.

At $\alpha = 2\pi$ the target is back where it started, yet the ancilla flips with certainty; at $4\pi$ it does not:

```wl
fringe /. {{α -> 2 π}, {α -> 4 π}}
```
<!-- run time: < 0.001 s -->

## Rabi and Landau-Zener Propagators: Elementary Functions and Parabolic-Cylinder Functions

A qubit driven by a field rotating at frequency $\omega$, $H(t) = \tfrac{\omega_0}{2}Z + \tfrac{\Omega}{2}(\cos\omega t\,X + \sin\omega t\,Y)$, and a qubit swept through an avoided crossing, $H(t) = \tfrac{vt}{2}Z + \tfrac{\Delta}{2}X$, both have Hamiltonians that do not commute with themselves at different times. What are their propagators $U(t)$? A numerical solver integrates the time-ordered evolution for one set of parameters. Here every parameter stays a symbol.

Compute the rotating-field propagator and its transition probability $|\langle 1|U(t)|0\rangle|^2$:

```wl
$Assumptions = ω0 ∈ Reals && ω ∈ Reals && Ω > 0 && t ∈ Reals;
urab = QuantumEvolve[ω0/2 QuantumOperator["Z"] + Ω/2 (Cos[ω t] QuantumOperator["X"] + Sin[ω t] QuantumOperator["Y"]), None, t];
FullSimplify[ComplexExpand[Abs[Normal[urab["Matrix"]][[2, 1]]]^2]]
```
<!-- run time: 0.48 s -->

This is the Rabi formula, with the generalized Rabi frequency $\sqrt{\Omega^2 + (\omega - \omega_0)^2}$ setting the oscillation and the detuning $\delta = \omega - \omega_0$ limiting the transfer to $\Omega^2/(\Omega^2 + \delta^2)$.

Compute the Landau-Zener propagator and list the special functions it is built from:

```wl
$Assumptions = True;
lz = Normal[QuantumEvolve[v t/2 QuantumOperator["Z"] + Δ/2 QuantumOperator["X"], None, t]["Matrix"]];
Union[Cases[lz, (h : Hypergeometric1F1 | ParabolicCylinderD | HermiteH)[___] :> h, Infinity]]
```
<!-- run time: 0.12 s -->

The propagator is built from parabolic-cylinder functions $D_\nu$, the solutions of the Weber equation that the linear sweep turns the Schrodinger equation into. Simplified, its first entry is a Kummer function:

```wl
FullSimplify[lz[[1, 1]]]
```
<!-- run time: 0.71 s -->

This is $e^{-ivt^2/4}\,{}_1F_1\!\left(\tfrac{i\Delta^2}{8v};\ \tfrac12;\ \tfrac{i}{2}vt^2\right)$: the single ratio $\Delta^2/v$ sets the index, and with it whether the sweep is adiabatic or diabatic. In other words, the whole time-resolved crossing, including the oscillations around the asymptotic Landau-Zener probability, is in one closed form.

Check that the matrix solves $i\,\dot U = H(t)\,U$, at $v = \Delta = 1$ and $t = 0.7$:

```wl
With[{m = lz /. {v -> 1, Δ -> 1}}, Chop[N[(I D[m, t] - (t/2 PauliMatrix[3] + PauliMatrix[1]/2) . m) /. t -> 0.7]]]
```
<!-- run time: 0.03 s -->

## A Beam Splitter with Symbolic Photon Numbers: the Amplitude Is a Wigner D-Function

Send $n_1$ and $n_2$ photons into the two ports of a beam splitter with mixing angle $\theta$ and phase $\varphi$. What is the amplitude to find $t_1$ and $t_2$ photons at the outputs, for every photon number at once? A Fock-space simulator truncates the space and answers one set of photon numbers at a time. Here the four photon numbers stay symbols.

Ask for $\langle t_1, t_2|U(\theta,\varphi)|n_1, n_2\rangle$ with $n_1, n_2, t_1, t_2, \theta, \varphi$ all symbolic:

```wl
QuantumInterferometer[{"BS"[θ, φ]}]["Amplitude", {n1, n2} -> {t1, t2}]
```
<!-- run time: 0.04 s -->

The answer is a single Wigner D-function times a Kronecker delta that conserves the photon number. In other words, a beam splitter is a rotation by $2\theta$ of a spin $j = (n_1+n_2)/2$, with the photon-number difference playing the role of $m$: two modes holding $N$ photons carry the spin-$N/2$ representation of SU(2) (Schwinger's construction), and the beam splitter rotates it.

Make the splitter balanced ($\theta = \pi/4$, $\varphi = 0$), send $n$ photons into each port, and ask for $k$ photons in the first output:

```wl
hom = QuantumInterferometer[{"BS"[π/4, 0]}]["Amplitude", {n, n} -> {k, 2 n - k}]
```
<!-- run time: 0.004 s -->

This is the Hong-Ou-Mandel effect for any $n$. At $n = k = 1$ it is $d^{1}_{00}(\pi/2) = \cos(\pi/2) = 0$: two single photons never leave through different ports. The same $d$-function vanishes whenever $k$ is odd, so with $n$ photons in each input, each output port only ever holds an even number of photons.

Confirm that every odd-$k$ amplitude is zero for $n = 1$ to $6$:

```wl
Union @ Flatten @ Table[hom, {n, 6}, {k, 1, 2 n - 1, 2}]
```
<!-- run time: 0.002 s -->

## Photon Loss: a Binomial Law for Symbolic Photon Numbers, and a Hong-Ou-Mandel Dip That Survives It

An attenuator with transmissivity $\eta$ lets each photon through with probability $\eta$. What is the probability that $n$ photons leave $k$ survivors, for every $n$ and $k$? And if both outputs of a balanced beam splitter lose light, does the Hong-Ou-Mandel zero survive? A density-matrix simulator answers both at one value of $\eta$ and one truncated photon number at a time. Here $\eta$, $n$ and $k$ stay symbols.

Ask for $P(k \mid n)$ through a pure-loss element with $n$, $k$ and $\eta$ symbolic:

```wl
loss = QuantumInterferometer[{"Loss"[η]}];
loss["Probability", {n} -> {k}]
```
<!-- run time: 0.09 s -->

The answer is the binomial law $\binom{n}{k}\eta^k(1-\eta)^{n-k}$: each photon is transmitted or lost independently of the others. The surviving photon number of a Fock state is therefore never sharp once $0 < \eta < 1$.

Put equal loss $\eta$ on both outputs of a balanced beam splitter and send one photon into each input:

```wl
lossyHOM = QuantumInterferometer[{"BS"[π/4], "Loss"[η] -> 1, "Loss"[η] -> 2}];
lossyHOM["Probabilities", {1, 1}]
```
<!-- run time: 0.02 s -->

The coincidence $\{1,1\}$ stays exactly 0 for every $\eta$, and each bunched outcome carries $\eta^2/2$. In other words, loss acts after the interference: the pair always leaves together, and loss only decides how many of the two photons reach the detectors. That is why equal detector inefficiency lowers the count rate but never fills the Hong-Ou-Mandel dip.

Check that the binomial law is normalized for symbolic $n$:

```wl
Sum[loss["Probability", {n} -> {k}], {k, 0, n}]
```
<!-- run time: 0.03 s -->

## Two Squeezers into a Balanced Beam Splitter: a Two-Mode Squeezed Vacuum

Squeeze the vacuum in each of two modes, with complex parameters $\xi_1$ and $\xi_2$, then mix the modes on a balanced beam splitter. When does one photon appear in each output, and what is the full photon-number distribution? Fock-space simulators such as Strawberry Fields answer at fixed numerical squeezing with a cutoff. Here the squeezing parameters stay symbols.

Ask for the coincidence amplitude $\langle 1,1|U|0,0\rangle$ with $\xi_1$, $\xi_2$ symbolic:

```wl
QuantumInterferometer[{"S"[ξ1] -> 1, "S"[ξ2] -> 2, "BS"[π/4, 0]}]["Amplitude", {0, 0} -> {1, 1}]
```
<!-- run time: 0.007 s -->

The amplitude is a difference of the two squeezers' pair amplitudes $e^{i\arg\xi}\tanh|\xi|$. Equal squeezing cancels it exactly: the splitter maps two identical squeezed vacua to two identical squeezed vacua, which hold only even photon numbers. Opposite squeezing, $\xi_2 = -\xi_1$, makes the two terms add, so the coincidence is largest there.

Take opposite squeezing $r$ and $-r$ with $r > 0$, and compute $P(m, m)$ for $m = 0, \dots, 3$:

```wl
tmsv = QuantumInterferometer[{"S"[r] -> 1, "S"[-r] -> 2, "BS"[π/4, 0]}];
Table[tmsv["Probability", {0, 0} -> {m, m}, Assumptions -> r > 0], {m, 0, 3}]
```
<!-- run time: 0.005 s -->

The two outputs always hold the same photon number, with $P(m,m) = \tanh^{2m} r/\cosh^2 r$: the two-mode squeezed vacuum, whose photons are made in pairs. Either mode alone is thermal. At $r = 0$ only $m = 0$ survives, which is the vacuum. Photon numbers in a squeezed network must be integers, so the general $m$ is read from this pattern rather than computed.

Check the pattern against the computed probabilities up to $m = 6$, and check that it sums to 1 over all $m$:

```wl
{FullSimplify[Table[tmsv["Probability", {0, 0} -> {m, m}, Assumptions -> r > 0] - Tanh[r]^(2 m)/Cosh[r]^2, {m, 0, 6}], r > 0],
 Sum[Tanh[r]^(2 m)/Cosh[r]^2, {m, 0, Infinity}, Assumptions -> r > 0]}
```
<!-- run time: 0.20 s -->

## The Three-Mode Fourier Interferometer: Exact Zeros in the Output

Send one photon into each port of a three-mode Fourier interferometer (a tritter). Which output patterns occur, and with what probability? A numerical simulator returns the forbidden patterns as numbers near $10^{-32}$, and the reader must decide whether they are zero. Built from the exact Fourier matrix, every probability is an exact rational.

Compute the full output distribution for $|1,1,1\rangle$:

```wl
QuantumInterferometer[FourierMatrix[3]]["Probabilities", {1, 1, 1}]
```
<!-- run time: 0.12 s -->

Every pattern with two photons in one mode and one in another is exactly 0, $|1,1,1\rangle$ survives at $1/3$, and the three fully bunched patterns share the rest at $2/9$ each. In other words, the cyclic symmetry of the Fourier matrix gives every history leading to a 2+1 pattern a partner with a cube-root-of-unity phase, and the three cancel. This is the suppression law that experiments use to certify that photons are indistinguishable.

Check that the distribution is normalized:

```wl
Total[QuantumInterferometer[FourierMatrix[3]]["Probabilities", {1, 1, 1}]]
```
<!-- run time: 0.10 s -->

## Cascaded Beam Splitters: Symbolic Amplitudes with Sums and Path Conditions

Chain beam splitters across three modes, and ask for the amplitude from $|n_1, n_2, n_3\rangle$ to $|t_1, t_2, t_3\rangle$ with every photon number symbolic. A single beam splitter gives one Wigner D-function. A cascade has to sum over the photon number inside the network and to know which output a photon can reach at all. Every Fock-space or permanent code needs the photon numbers as integers.

Two splitters, on modes 1, 2 and then on modes 2, 3:

```wl
QuantumInterferometer[{"BS"[θ1, ϕ1] -> {1, 2}, "BS"[θ2, ϕ2] -> {2, 3}}, 3]["Amplitude", {n1, n2, n3} -> {t1, t2, t3}]
```
<!-- run time: 0.01 s -->

The amplitude is a product of two D-functions under one condition, $n_1 + n_2 \ge t_1$. Photons that start in mode 3 meet mode 1 only through the first splitter, which they never reach, so mode 1 can hold at most the $n_1 + n_2$ photons that went through it. The Kronecker delta conserves the total.

Add a third splitter back on modes 1, 2, so photons can return to mode 1:

```wl
mesh = QuantumInterferometer[{"BS"[θ1, ϕ1] -> {1, 2}, "BS"[θ2, ϕ2] -> {2, 3}, "BS"[θ3, ϕ3] -> {1, 2}}, 3];
mesh["Amplitude", {n1, n2, 0} -> {t1, t2, t3}]
```
<!-- run time: 0.007 s -->

Now the amplitude is an inactive sum over $K_1$, the number of photons in mode 1 between the first and the last splitter, with the bound $\min(n_1+n_2,\ t_1+t_2)$ set by the photon flow. The sum stays unevaluated so that the formula holds for every photon number.

Substitute $|2,1,0\rangle \to |1,1,1\rangle$, activate the sum, and compare with the direct query at those numbers:

```wl
FullSimplify[Activate[mesh["Amplitude", {n1, n2, 0} -> {t1, t2, t3}] /. {n1 -> 2, n2 -> 1, t1 -> 1, t2 -> 1, t3 -> 1}] - mesh["Amplitude", {2, 1, 0} -> {1, 1, 1}]]
```
<!-- run time: 0.02 s -->

## The Nonlinear Sign Gate: Heralded Output for a Symbolic Input

Linear optics cannot give the two-photon component a sign flip on its own: a phase shifter multiplies $|n\rangle$ by $e^{in\phi}$, so it changes $|1\rangle$ as well. Knill, Laflamme and Milburn obtained $c_0|0\rangle + c_1|1\rangle + c_2|2\rangle \mapsto c_0|0\rangle + c_1|1\rangle - c_2|2\rangle$ by adding one auxiliary photon and one vacuum mode, and keeping only the runs in which the auxiliary detectors see $\{1, 0\}$. What comes out for an arbitrary input, and how often does it work? Here $c_0, c_1, c_2$ stay symbols.

Build the three-mode network from its real transfer matrix (Knill, Laflamme and Milburn, Nature 409, 46 (2001)):

```wl
ns = QuantumInterferometer[{
   {1 - Sqrt[2], 2^(-1/4), Sqrt[3/Sqrt[2] - 2]},
   {2^(-1/4), 1/2, (1 - Sqrt[2])/2},
   {Sqrt[3/Sqrt[2] - 2], (1 - Sqrt[2])/2, Sqrt[2] - 1/2}}];
RootReduce[Table[ns["Amplitude", {m, 1, 0} -> {m, 1, 0}], {m, 0, 2}]]
```
<!-- run time: 0.18 s -->

The heralded amplitudes $\langle m,1,0|U|m,1,0\rangle$ are $1/2$, $1/2$ and $-1/2$: the same size for $m = 0, 1, 2$, with the sign flipped only for two photons.

Send the symbolic signal with the auxiliary photon and vacuum through the network (four Fock levels per mode hold every sector), and keep the signal amplitudes where the auxiliary modes read $\{1, 0\}$:

```wl
nsOut = ns[QuantumTensorProduct[QuantumState[{c0, c1, c2, 0}, 4], FockState[1, 4], FockState[0, 4]]];
RootReduce[ArrayReshape[Normal[nsOut["StateVector"]], {4, 4, 4}][[All, 2, 1]]]
```
<!-- run time: 0.06 s -->

The heralded signal is $\tfrac12(c_0, c_1, -c_2)$. In other words, the gate acts identically on every input, and its success probability is $\tfrac14(|c_0|^2 + |c_1|^2 + |c_2|^2) = \tfrac14$ whatever the input state is. That independence of the input is what lets the gate serve as a building block for a controlled-sign gate between photons.

Check that the gate is defined only up to two signal photons: the heralded amplitude for $m = 3$ is neither $1/2$ nor $-1/2$:

```wl
RootReduce[ns["Amplitude", {3, 1, 0} -> {3, 1, 0}]]
```
<!-- run time: 0.002 s -->

## From Interferometer to Fock-Space Circuit: Symbolic Parameters Carry Over

The permanent formulas never build a Fock space. When a computation needs one (a mixed input, a partial trace, a QuantumFramework measurement), `qi["CircuitOperator", d]` turns the interferometer into a `QuantumCircuitOperator` on $d$ levels per mode, with the beam-splitter angles still symbolic.

Build the Fock-space circuit of a beam splitter with symbolic $\theta$, $\varphi$ on three levels per mode:

```wl
bsCircuit = QuantumInterferometer[{"BS"[θ, ϕ]}]["CircuitOperator", 3]
```
<!-- run time: 0.02 s -->

Check it against the permanent route: $\langle 2,0|\,\text{circuit}\,|1,1\rangle$ minus the interferometer's amplitude for $\{1,1\} \to \{2,0\}$:

```wl
FullSimplify[(SuperDagger[FockState[{2, 0}, 3]] @ bsCircuit[FockState[{1, 1}, 3]])["Scalar"] - QuantumInterferometer[{"BS"[θ, ϕ]}]["Amplitude", {1, 1} -> {2, 0}], Element[{θ, ϕ}, Reals]]
```
<!-- run time: 0.14 s -->

## Bosons, Distinguishable Particles and Fermions at One Beam Splitter: Permanent, Probability Permanent, Determinant

Send two particles, one into each input of a beam splitter with transfer matrix $W$. The coincidence amplitude adds the two histories (both transmitted, both reflected). Identical bosons add them with a plus sign, which gives the permanent of $W$; identical fermions add them with a minus sign, which gives the determinant; distinguishable particles add probabilities, which gives the permanent of $|W_{jk}|^2$. All three come out as formulas in $\theta$ and $\varphi$.

Compute $\mathrm{Per}\,W$, $\det W$ and $\mathrm{Per}\,|W|^2$ for a beam splitter with symbolic $\theta$, $\varphi$:

```wl
wbs = QuantumInterferometer[{"BS"[θ, ϕ]}]["Unitary"];
FullSimplify[{Permanent[wbs], Det[wbs], Permanent[Abs[wbs]^2]}, Element[{θ, ϕ}, Reals]]
```
<!-- run time: 0.01 s -->

The boson coincidence amplitude is $\cos 2\theta$, which vanishes at the balanced splitter (Hong-Ou-Mandel). The fermion amplitude is $\det W = 1$ for every $\theta$: two identical fermions always leave through different ports. Distinguishable particles give $\cos^4\theta + \sin^4\theta$, which never drops below $1/2$. In other words, one sign under exchange separates complete bunching from complete antibunching, and the same sign separates the permanent, for which no efficient algorithm is known, from the determinant, which Gaussian elimination computes in polynomial time.

Check that the permanent is the amplitude `QuantumInterferometer` itself returns for $\{1,1\} \to \{1,1\}$:

```wl
FullSimplify[Permanent[wbs] - QuantumInterferometer[{"BS"[θ, ϕ]}]["Amplitude", {1, 1} -> {1, 1}], Element[{θ, ϕ}, Reals]]
```
<!-- run time: 0.001 s -->

## A Coherent State Under the Oscillator: It Stays Coherent, with the Phase $e^{-i\omega t/2}$

Under $H = \omega(a^\dagger a + \tfrac12)$, a coherent state should stay coherent with its label rotating, $e^{-iHt}|\alpha\rangle = e^{-i\omega t/2}\,|\alpha e^{-i\omega t}\rangle$. Can this be checked for all $\alpha$, $\omega$ and $t$ at once? QuTiP would propagate one $\alpha$ on a time grid. Here the propagator is computed with $\omega$ and $t$ symbolic and applied to the coherent state with symbolic $\alpha$, in QuantumFramework's default truncated Fock space of 16 levels.

Compute the propagator $e^{-iHt}$ and read its diagonal:

```wl
aOp = AnnihilationOperator[];
uosc = QuantumEvolve[ω (aOp["Dagger"] @ aOp + 1/2), None, t];
Diagonal[Normal[uosc["Matrix"]]]
```
<!-- run time: 2.4 s -->

The propagator is diagonal in the photon number, with phases $e^{-i\omega(n+1/2)t}$. So it multiplies the $n$-th amplitude $\alpha^n/\sqrt{n!}$ by $e^{-i\omega t/2}(e^{-i\omega t})^n$, which is exactly the coherent state with label $\alpha e^{-i\omega t}$ times the vacuum phase.

Compare amplitude by amplitude, using unnormalized coherent states so that the comparison includes the phase $e^{-i\omega t/2}$ and does not depend on how the truncated state is normalized:

```wl
Simplify[(uosc @ CoherentState["Normalized" -> False])["AmplitudesList"] -
  E^(-I ω t/2) CoherentState["Normalized" -> False][\[FormalAlpha] E^(-I ω t)]["AmplitudesList"]]
```
<!-- run time: 0.09 s -->

All 16 amplitudes agree for every $\alpha$, $\omega$, $t$. The statement holds at the 16-level cutoff, amplitude by amplitude. Because the propagator is diagonal, the truncation does not mix levels, and the same identity holds at any cutoff.

Check with the state comparison `==` on normalized states, which ignores the global phase and so tests only the rotated label:

```wl
uosc @ CoherentState[] == CoherentState[][\[FormalAlpha] E^(-I ω t)]
```
<!-- run time: 0.08 s -->

## Hong-Ou-Mandel with Partially Distinguishable Photons: the Dip Depth Is the Polarization Overlap

Two photons meet at a beam splitter with angle $\theta$. The first is horizontally polarized; the second is polarized at angle $\chi$, $\cos\chi\,|H\rangle + \sin\chi\,|V\rangle$, and the detectors do not resolve polarization. How does the coincidence probability depend on $\theta$ and $\chi$? Treat the polarizations as two extra modes, ordered $1_H, 2_H, 1_V, 2_V$, with one beam splitter acting on the horizontal pair and an identical one on the vertical pair. The input is a superposition of $|1,1,0,0\rangle$ and $|1,0,0,1\rangle$, so each output amplitude is the same superposition of two interferometer amplitudes.

Sum $|\cos\chi\,A(1100 \to T) + \sin\chi\,A(1001 \to T)|^2$ over the four outputs $T$ with one photon in each path:

```wl
pol = QuantumInterferometer[{"BS"[θ] -> {1, 2}, "BS"[θ] -> {3, 4}}, 4];
coinc = Total[Abs[Cos[χ] pol["Amplitude", {1, 1, 0, 0} -> #] + Sin[χ] pol["Amplitude", {1, 0, 0, 1} -> #]]^2 & /@
    {{1, 1, 0, 0}, {1, 0, 0, 1}, {0, 1, 1, 0}, {0, 0, 1, 1}}];
FullSimplify[coinc == Cos[θ]^4 + Sin[θ]^4 - Sin[2 θ]^2 Cos[χ]^2/2, Element[{θ, χ}, Reals]]
```
<!-- run time: 0.16 s -->

The coincidence probability is $\cos^4\theta + \sin^4\theta - \tfrac12 \sin^2 2\theta\,\cos^2\chi$: the classical value for distinguishable particles, minus an interference term weighted by the squared overlap $\cos^2\chi$ of the two polarization states.

Take the limits: identical photons ($\chi = 0$), orthogonal polarizations ($\chi = \pi/2$), and the balanced splitter ($\theta = \pi/4$):

```wl
FullSimplify[{coinc /. χ -> 0, coinc /. χ -> π/2, coinc /. θ -> π/4}, Element[{θ, χ}, Reals]]
```
<!-- run time: 0.01 s -->

Identical photons give $\cos^2 2\theta$, the Hong-Ou-Mandel result; orthogonal polarizations give $\cos^4\theta + \sin^4\theta$, the result for distinguishable particles. At the balanced splitter the coincidence is $\tfrac12\sin^2\chi$: the depth of the Hong-Ou-Mandel dip measures how indistinguishable the photons are.

## Entanglement Sudden Death: the Death Time and Its Condition from One Reduce

Two qubits start in $\cos\theta\,|00\rangle + \sin\theta\,|11\rangle$, and each decays on its own from $|1\rangle$ to $|0\rangle$ at rate $\gamma$. Every local population decays smoothly, yet the entanglement can vanish at a finite time (Yu and Eberly). When, and for which $\theta$? A numerical solver answers this on a grid of $\theta$ and $t$. Here $\theta$, $\gamma$ and $t$ stay symbols, so the answer is a formula together with the range of $\theta$ where it holds.

Solve the master equation with jump operator $\sigma_- = |0\rangle\langle 1|$ on each qubit, and no Hamiltonian (in the frame rotating with the qubits the free evolution only adds phases, which the entanglement does not see):

```wl
$Assumptions = γ > 0 && t >= 0 && 0 < θ < π/2;
lower = {{0, 1}, {0, 0}};
ρpair = QuantumEvolve[QuantumOperator[0 IdentityMatrix[4]],
    {QuantumOperator[lower, {1}], QuantumOperator[lower, {2}]} -> {γ, γ},
    QuantumState[{Cos[θ], 0, 0, Sin[θ]}], t];
MatrixForm[dm = FullSimplify[Normal[ρpair["DensityMatrix"]]]]
```
<!-- run time: 0.36 s -->

The state stays an X-state: the coherence $\rho_{14}$ between $|00\rangle$ and $|11\rangle$ decays as $e^{-\gamma t}$, while the decay of $|11\rangle$ fills the one-excitation populations $\rho_{22} = \rho_{33}$. For this X-state, where $\rho_{23} = 0$, the concurrence is $\max\left(0,\ 2(|\rho_{14}| - \sqrt{\rho_{22}\rho_{33}})\right)$.

Compute $2(|\rho_{14}| - \sqrt{\rho_{22}\rho_{33}})$:

```wl
conc = FullSimplify[2 (Abs[dm[[1, 4]]] - Sqrt[dm[[2, 2]] dm[[3, 3]]])]
```
<!-- run time: 0.13 s -->

Find when it reaches zero at a positive time:

```wl
FullSimplify[Reduce[conc == 0 && t > 0 && γ > 0 && 0 < θ < π/2, t, Reals]]
```
<!-- run time: 1.4 s -->

The entanglement dies at $t^* = \frac{1}{\gamma}\ln\frac{1}{1-\cot\theta}$, and only when $\theta > \pi/4$, that is, when the doubly excited $|11\rangle$ carries more weight than $|00\rangle$. In other words, the formula and the condition for it to exist come out of the same computation. Below $\pi/4$ the coherence always wins and the entanglement only fades asymptotically; as $\theta \to \pi/4$ from above, $t^*$ diverges.

Check against QuantumFramework's own concurrence, which knows nothing of the X-state formula: at $\theta = \pi/3$ and $\gamma = 1$ it should still be positive at $0.9\,t^*$ and vanish exactly at $t^*$:

```wl
QuantumEntanglementMonotone[QuantumState[dm /. {θ -> π/3, γ -> 1, t -> #}], {{1}, {2}}, "Concurrence"] & /@ (Log[1/(1 - Cot[π/3])] {9/10, 1}) // N
```
<!-- run time: 0.18 s -->

## A Qutrit Wigner Function Under Depolarizing Noise: the Negativity Dies at an Exact Time

A qutrit state can be written as a discrete Wigner function: nine real numbers on a $3 \times 3$ phase-space grid that sum to one. Negative entries are what separate it from a probability distribution, and their total size, the mana of Veitch et al. (2014), is a resource monotone for magic-state quantum computation. Under noise the negative entries fill in. At what time, exactly, does the last one reach zero? A numerical solver gives $W(t)$ on a time grid and the crossing to the grid spacing. Here the rate $\gamma$ and the time $t$ stay symbols.

The strange state $(|1\rangle - |2\rangle)/\sqrt{2}$ and its Wigner function on the grid:

```wl
strange = QuantumState[{0, 1, -1}/Sqrt[2]];
Normal[strange["PhaseSpace"]]
```
<!-- run time: 0.16 s -->

All the negativity sits at the origin of phase space, and $-1/3$ is the most negative value a qutrit Wigner function can take (the phase-point operators have eigenvalues $\pm 1$).

Take depolarizing noise built from the eight non-identity Weyl operators $X^a Z^b$, each a jump operator at rate $\gamma/9$, and write its Lindblad generator $G$ in the Wigner basis. List its distinct entries and its column sums:

```wl
$Assumptions = γ > 0 && t >= 0;
{x3, z3} = Normal[QuantumOperator[#]["Matrix"]] & /@ {"X"[3], "Z"[3]};
weyl = QuantumOperator[MatrixPower[x3, #1] . MatrixPower[z3, #2]] & @@@ Rest[Tuples[Range[0, 2], 2]];
G = Normal[QuantumPhaseSpaceTransform[QuantumOperator["Liouvillian"[None, weyl, ConstantArray[γ/9, 8]]]]["Matrix"]];
{Union[Flatten[G]], Total[G]}
```
<!-- run time: 0.95 s -->

The generator is real, its off-diagonal entries are all $\gamma/9$, and every column sums to zero. In other words, in this representation the Lindblad equation is a classical rate equation: each Weyl operator translates the phase-space grid, so the noise is a random walk that hops from every point to each of the other eight at rate $\gamma/9$, carrying the negative entry along with it.

Solve the master equation from the strange state and read the Wigner function at time $t$:

```wl
wt = Simplify[Normal[QuantumEvolve[QuantumOperator[0 IdentityMatrix[3]], weyl -> ConstantArray[γ/9, 8], strange, t]["PhaseSpace"]]]
```
<!-- run time: 0.24 s -->

The origin entry relaxes as $\frac19 - \frac49 e^{-\gamma t}$ toward the uniform value $\frac19$, while the other eight decay toward it from above. Find when the origin entry reaches zero:

```wl
Reduce[wt[[1, 1]] == 0 && t > 0 && γ > 0, t, Reals]
```
<!-- run time: 0.55 s -->

The negativity dies at $t^* = 2\ln 2/\gamma$, the moment the surviving fraction $e^{-\gamma t}$ of the strange state falls to $1/4$. From then on all nine entries are nonnegative, the Wigner function is a probability distribution, and the mana is zero.

Check the long-time limit, the maximally mixed state, whose Wigner function is uniform:

```wl
Limit[wt, t -> Infinity, Assumptions -> γ > 0]
```
<!-- run time: 0.06 s -->

## The Amplitude-Damping Master Equation: It Is the Kraus Channel at $1 - e^{-\gamma t}$

A qubit decays from $|1\rangle$ to $|0\rangle$ at rate $\gamma$. Integrating its master equation for a time $t$ gives a linear map on density matrices. Is that map the textbook amplitude-damping channel, with Kraus operators $\mathrm{diag}(1, \sqrt{1-\lambda})$ and $\sqrt{\lambda}\,|0\rangle\langle 1|$, and at which damping parameter $\lambda$? A numerical tool compares the two on sample states and times. Here the initial state is a general qubit density matrix with symbolic population $p$ and coherence $c$, and $\gamma$ and $t$ are symbols.

Solve the master equation with jump operator $\sigma_- = |0\rangle\langle 1|$ from the general state $\rho_0$:

```wl
$Assumptions = γ > 0 && t >= 0 && 0 <= p <= 1 && t1 >= 0 && t2 >= 0;
σm = QuantumOperator[{{0, 1}, {0, 0}}];
ρ0 = QuantumState[{{p, c}, {Conjugate[c], 1 - p}}];
me = Normal[QuantumEvolve[QuantumOperator[0 IdentityMatrix[2]], {σm} -> {γ}, ρ0, t]["DensityMatrix"]]
```
<!-- run time: 0.07 s -->

The excited population decays at $\gamma$ and the coherence at $\gamma/2$. Compare with QuantumFramework's named amplitude-damping channel at $\lambda = 1 - e^{-\gamma t}$ acting on the same symbolic state:

```wl
FullSimplify[me == Normal[QuantumChannel["AmplitudeDamping"[1 - E^(-γ t)]][ρ0]["DensityMatrix"]]]
```
<!-- run time: 0.45 s -->

In other words, the master equation integrated for time $t$ is the Kraus channel at $\lambda = 1 - e^{-\gamma t}$, for every initial state and every rate.

Since evolution for $t_1$ followed by $t_2$ is evolution for $t_1 + t_2$, composing two such channels must give the channel at the summed time. Check the composition as channel objects:

```wl
ad[s_] := QuantumChannel["AmplitudeDamping"[1 - E^(-γ s)]];
FullSimplify[Normal[(ad[t1] @ ad[t2])[ρ0]["DensityMatrix"]] == Normal[ad[t1 + t2][ρ0]["DensityMatrix"]]]
```
<!-- run time: 1.05 s -->

## Two Atoms Decaying Into a Shared Field: Rates $1 \pm r$ and an Exactly Dark Singlet

Two two-level atoms close together radiate into the same field, so their decay is correlated: the Lindblad equation carries a rate matrix $\gamma\begin{pmatrix}1 & r\\ r & 1\end{pmatrix}$ with correlation $0 \le r \le 1$ (the Kossakowski form). How fast do the singlet and triplet decay, and for which $r$ does a state stop decaying? QuTiP takes a list of collapse operators, so the rate matrix must first be diagonalized by hand, and $r$ is a number in each run. Here $r$ stays a symbol, and time is in units of $1/\gamma$.

Build the Liouvillian with the full rate matrix and ask for its spectrum:

```wl
$Assumptions = 0 <= r <= 1 && t >= 0;
{L1, L2} = {QuantumOperator[{{0, 1}, {0, 0}}, {1}], QuantumOperator[{{0, 1}, {0, 0}}, {2}]};
pair = QuantumOperator["Liouvillian"[QuantumOperator[0 IdentityMatrix[4]], {L1, L2}, {{1, r}, {r, 1}}]];
Eigenvalues[Normal[pair["Matrix"]]]
```
<!-- run time: 0.78 s -->

The spectrum contains the two collective decay rates $1 - r$ and $1 + r$, the eigenvalues of the rate matrix: the environment couples to the antisymmetric and symmetric combinations of the two atoms, not to each atom. At $r = 1$ several eigenvalues reach zero, so the set of stationary states grows.

Solve the master equation for the singlet $(|01\rangle - |10\rangle)/\sqrt2$ and the triplet $(|01\rangle + |10\rangle)/\sqrt2$, and read the mean number of excitations:

```wl
excitation[ψ_] := Simplify[Tr[Normal[QuantumEvolve[QuantumOperator[0 IdentityMatrix[4]], {L1, L2} -> {{1, r}, {r, 1}}, ψ, t]["DensityMatrix"]] . DiagonalMatrix[{0, 1, 1, 2}]]];
excitation /@ {QuantumState[{0, 1, -1, 0}/Sqrt[2]], QuantumState[{0, 1, 1, 0}/Sqrt[2]]}
```
<!-- run time: 0.70 s -->

The singlet decays at the subradiant rate $1 - r$ and the triplet at the superradiant rate $1 + r$: the same two atoms in the same environment, distinguished only by a relative sign. In other words, the two emission amplitudes interfere destructively in the singlet and constructively in the triplet.

For which $r$ is the singlet an exact stationary state, $\mathcal{L}[\rho_S] = 0$?

```wl
singlet = QuantumState[{0, 1, -1, 0}/Sqrt[2]];
Reduce[Thread[Normal[pair["Matrix"]] . Flatten[Normal[singlet["DensityMatrix"]]] == 0], r]
```
<!-- run time: 0.002 s -->

Only at perfect correlation: there the singlet is a decoherence-free state, frozen for all time, which is the $r \to 1$ end of its decay rate $1 - r$.

## A Coherent Error on the Three-Qubit Code: Measurement Records Undo It for Every Angle

Encode $a|0\rangle + b|1\rangle$ as $a|000\rangle + b|111\rangle$, hit the middle qubit with a coherent rotation $e^{-i\theta X}$ (not a Pauli error), measure the parities $Z_1Z_2$ and $Z_2Z_3$, and correct with $X$ gates conditioned on the outcomes. The textbook argument (Nielsen and Chuang, §10.2) says the syndrome measurement turns the rotation into either no error or an $X_2$ error, with weights $\cos^2\theta$ and $\sin^2\theta$. A simulator with mid-circuit measurement samples one outcome per shot and estimates the fidelity. QuantumFramework keeps each measurement outcome as a pointer qubit, so the correction is an ordinary controlled gate and the whole cycle is one symbolic state, with $a$, $b$ and $\theta$ left as symbols.

A parity measurement given as the two projectors $(1 \pm Z_iZ_j)/2$, which learns the parity and leaves superpositions inside each parity sector intact:

```wl
$Assumptions = θ ∈ Reals;
zz = KroneckerProduct[PauliMatrix[3], PauliMatrix[3]];
syndrome[i_, j_] := QuantumMeasurementOperator[{QuantumOperator[(IdentityMatrix[4] + zz)/2, {i, j}], QuantumOperator[(IdentityMatrix[4] - zz)/2, {i, j}]}, {i, j}];
```
<!-- run time: 0.0006 s -->

The circuit: encode, rotate, record the two parities on pointer wires $0$ and $-1$, then apply $X_1$, $X_2$ or $X_3$ controlled on the record patterns $(1,0)$, $(1,1)$, $(0,1)$:

```wl
qec = QuantumCircuitOperator[{"CNOT" -> {1, 2}, "CNOT" -> {1, 3}, QuantumOperator["RX"[2 θ], {2}], syndrome[1, 2], syndrome[2, 3],
    QuantumOperator["C"[QuantumOperator["X", {1}], {0}, {-1}]],
    QuantumOperator["C"[QuantumOperator["X", {2}], {0, -1}, {}]],
    QuantumOperator["C"[QuantumOperator["X", {3}], {-1}, {0}]]}];
outCode = qec[QuantumTensorProduct[QuantumState[{a, b}], QuantumState["00"]]];
```
<!-- run time: 0.37 s -->

Compare the output, wires ordered $\{-1, 0, 1, 2, 3\}$, with the product of the record state $\cos\theta\,|00\rangle - i\sin\theta\,|11\rangle$ and the codeword:

```wl
Normal[outCode["StateVector"]] == Normal[QuantumTensorProduct[QuantumState[{Cos[θ], 0, 0, -I Sin[θ]}], QuantumState[a UnitVector[8, 1] + b UnitVector[8, 8]]]["StateVector"]]
```
<!-- run time: 0.009 s -->

The output factorizes exactly. In other words, the syndrome measurement has digitized the rotation into "no error" and "$X_2$" with weights $\cos^2\theta$ and $\sin^2\theta$, the record-controlled gates have fixed both branches, and the whole error now lives in the record qubits. The code qubits are back in $a|000\rangle + b|111\rangle$ for every $\theta$, with no classical communication step anywhere.

Check that the code qubits alone, tracing out the records, are in the pure codeword:

```wl
code = a UnitVector[8, 1] + b UnitVector[8, 8];
FullSimplify[Normal[QuantumPartialTrace[outCode, {1, 2}]["DensityMatrix"]] == Outer[Times, code, Conjugate[code]]]
```
<!-- run time: 0.01 s -->

## Two Qubits Have No Bound Entanglement: an Identity and One Resolve

A bound-entangled state is entangled, yet its partial transpose is still positive, so no local protocol can distill a Bell pair from it. Can a two-qubit state be bound entangled? Sampling random states with a numerical package finds none, but a sample is evidence, not a proof. Here the statement is proved for every two-qubit density matrix, written with all 16 of its real Pauli coefficients left as symbols. The tool is Wootters' concurrence $C = \max(0, \lambda_1 - \lambda_2 - \lambda_3 - \lambda_4)$, where $\lambda_1 \ge \dots \ge \lambda_4$ are the square roots of the eigenvalues of $R = \rho\tilde\rho$, $\tilde\rho = (\sigma_y\otimes\sigma_y)\rho^*(\sigma_y\otimes\sigma_y)$, and $C = 0$ exactly on separable states.

Build $\rho = \sum_{ij} c_{ij}\,\sigma_i\otimes\sigma_j$, its partial transpose on the second qubit (transpose every $2\times2$ block), and $R$. The $\lambda_k^2$ are the eigenvalues of $R$, so $\sum\lambda_k^2 = \mathrm{Tr}R$, $\sum\lambda_k^4 = \mathrm{Tr}R^2$ and $\prod\lambda_k = \det\rho$. Test $16\det\rho^{T_B} = (\mathrm{Tr}R)^2 - 2\,\mathrm{Tr}R^2 + 8\det\rho$ in all 16 coefficients:

```wl
paulis = Flatten[Outer[KroneckerProduct, #, #, 1] &[PauliMatrix[{0, 1, 2, 3}]], 1];
ρ2q = Array[c, 16] . paulis;
ρ2qTB = ArrayFlatten[Map[Transpose, Partition[ρ2q, {2, 2}], {2}]];
R = ρ2q . KroneckerProduct[PauliMatrix[2], PauliMatrix[2]] . (Array[c, 16] . Conjugate[paulis]) . KroneckerProduct[PauliMatrix[2], PauliMatrix[2]];
Expand[16 Det[ρ2qTB] - (Tr[R]^2 - 2 Tr[R . R] + 8 Det[ρ2q])]
```
<!-- run time: 0.02 s -->

The right side is a symmetric polynomial in the $\lambda_k$. Check that it equals $\prod_k(2\lambda_k - \sum_j\lambda_j)$:

```wl
With[{λ = Array[l, 4]}, Expand[Times @@ (2 λ - Total[λ]) - (Total[λ^2]^2 - 2 Total[λ^4] + 8 Times @@ λ)]]
```
<!-- run time: 0.0002 s -->

So $16\det\rho^{T_B} = \prod_k(2\lambda_k - \sum_j\lambda_j)$ for every two-qubit state. In other words, the sign of the partial transpose's determinant is a statement about the same four numbers that make up the concurrence, and the square roots have dropped out.

A positive partial transpose makes that determinant nonnegative. Ask whether, for every ordered quadruple $\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge \lambda_4 \ge 0$ with a nonnegative product, the concurrence is clipped to zero:

```wl
With[{λ = Array[l, 4]}, Resolve[ForAll[λ, GreaterEqual @@ Append[λ, 0] && Times @@ (2 λ - Total[λ]) >= 0, First[λ] - Total[Rest[λ]] <= 0], Reals]]
```
<!-- run time: 0.04 s -->

`True`: a positive partial transpose forces $C = 0$, which is separability. A two-qubit state is either entangled with a negative partial transpose or separable, so there is nothing left to be bound. For $2\times4$ and $3\times3$ systems a positive partial transpose no longer implies separability, and bound entangled states exist (Horodecki).

Cross-check on a Bell state, where $\lambda = (1, 0, 0, 0)$ makes the product $-1$:

```wl
16 Det[Normal[QuantumPartialTranspose[QuantumState["PhiPlus"], {2}]["DensityMatrix"]]]
```
<!-- run time: 0.03 s -->

## Robertson-Schrödinger for Every Qubit State and Every Pair of Observables: the Slack Factorizes

For two observables $A$, $B$ on a state $\rho$, the Robertson-Schrödinger relation says $\sigma_A^2\sigma_B^2 \ge \left(\tfrac12\langle\{A,B\}\rangle - \langle A\rangle\langle B\rangle\right)^2 + \left(\tfrac1{2i}\langle[A,B]\rangle\right)^2$. For a qubit, every observable is $a_0 I + \vec a\cdot\vec\sigma$ and the $a_0 I$ part drops out of both sides, so the general case is $A = \vec a\cdot\vec\sigma$, $B = \vec b\cdot\vec\sigma$ on $\rho = \tfrac12(I + \vec r\cdot\vec\sigma)$: nine real symbols. A numerical check draws random states and observables; here the inequality is decided for all of them at once.

Compute the slack, left side minus right side:

```wl
slackRS = With[{ρq = (IdentityMatrix[2] + {r1, r2, r3} . PauliMatrix[{1, 2, 3}])/2, A = {a1, a2, a3} . PauliMatrix[{1, 2, 3}], B = {b1, b2, b3} . PauliMatrix[{1, 2, 3}]},
   With[{ev = Tr[# . ρq] &},
    Expand[(ev[A . A] - ev[A]^2) (ev[B . B] - ev[B]^2) - (ev[(A . B + B . A)/2] - ev[A] ev[B])^2 - (ev[A . B - B . A]/(2 I))^2]]];
Factor[slackRS]
```
<!-- run time: 0.002 s -->

The slack is $|\vec a\times\vec b|^2\,(1 - |\vec r|^2)$: how far the two observables are from commuting, times how mixed the state is. It vanishes for a pure state or a commuting pair, and only then.

Ask `Reduce` for any physical Bloch vector and any two observables that violate the relation, working on the expanded polynomial in all nine variables:

```wl
Reduce[slackRS < 0 && r1^2 + r2^2 + r3^2 <= 1, {r1, r2, r3, a1, a2, a3, b1, b2, b3}, Reals]
```
<!-- run time: 0.45 s -->

`False`: no violation exists anywhere in the nine-dimensional space. In other words, the relation is a theorem about the polynomial, decided by quantifier elimination, not a sampled observation.

Cross-check that every pure state ($|\vec r| = 1$) saturates the relation:

```wl
Simplify[slackRS /. r3 -> Sqrt[1 - r1^2 - r2^2]]
```
<!-- run time: 0.0003 s -->

## Discriminating Two Qubit States: Reduce Finds the Helstrom Bound Over Every Measurement

One of two known states arrives, $\rho_1$ with prior $p$ or $\rho_2$ with prior $1-p$, and one measurement must say which. The measurement is a pair of effects $E_1 + E_2 = I$, $0 \le E_2 \le I$, and the error is $P_{\mathrm{err}} = p\,\mathrm{Tr}[\rho_1E_2] + (1-p)\,\mathrm{Tr}[\rho_2E_1] = (1-p) + \mathrm{Tr}[\Lambda E_2]$ with $\Lambda = p\rho_1 - (1-p)\rho_2$. What is the smallest error over all measurements? A numerical approach solves a semidefinite program for one pair of states. Here the minimum over every effect is derived, with the spectrum of $\Lambda$ left symbolic.

The set of effects does not change under a unitary, so work in the eigenbasis of $\Lambda$, $\Lambda = \mathrm{diag}(l_1, l_2)$, and take the most general effect $E_2 = \begin{pmatrix} a & b + ic \\ b - ic & e\end{pmatrix}$ with $E_2 \ge 0$ and $I - E_2 \ge 0$. Ask for every $k$ that bounds $\mathrm{Tr}[\Lambda E_2]$ from below for all such effects:

```wl
Reduce[ForAll[{a, b, c, e}, 0 <= a <= 1 && 0 <= e <= 1 && a e >= b^2 + c^2 && (1 - a) (1 - e) >= b^2 + c^2,
   Tr[DiagonalMatrix[{l1, l2}] . {{a, b + I c}, {b - I c, e}}] >= k], k, Reals]
```
<!-- run time: 0.04 s -->

The best bound is the sum of the negative eigenvalues of $\Lambda$, reached by the projector onto its negative eigenspace. With $\mathrm{Tr}\Lambda = 2p - 1$ that sum is $\tfrac12(2p - 1 - \|\Lambda\|_1)$, so $P_{\mathrm{err}} \ge \tfrac12(1 - \|\Lambda\|_1)$ for every measurement, with equality for that projector: the Helstrom bound, derived rather than quoted. When both eigenvalues are positive the best effect is $E_2 = 0$: always guess $\rho_1$, because the prior outweighs anything a measurement could reveal (and symmetrically for two negative eigenvalues).

Cross-check against QuantumFramework at equal priors, $|0\rangle$ against a pure state $60^\circ$ away on the Bloch sphere: the minimum error $\tfrac12$ plus the sum of the negative eigenvalues of $\Lambda$ should equal $(1 - D)/2$, with $D$ the trace distance computed by `QuantumDistance`:

```wl
With[{qa = QuantumState["0"], qb = QuantumState["BlochVector"[{Sin[π/3], 0, Cos[π/3]}]]},
 {1/2 + Total[Select[Eigenvalues[Normal[(qa["DensityMatrix"] - qb["DensityMatrix"])/2]], Negative]], (1 - QuantumDistance[qa, qb, "Trace"])/2}]
```
<!-- run time: 0.01 s -->

## A Continuously Measured Trapped Particle: Squeezing, Purity and the Efficiency Threshold in Closed Form

Watch the position of a particle in a harmonic trap continuously, with measurement rate $\lambda$ and detector efficiency $\eta$. Does the conditional state become narrower in position than the ground state, and for which measurement strengths? A trajectory simulation answers this at one set of numbers. Here the trap and the detector stay symbols. Because the Hamiltonian is quadratic, a Gaussian conditional state stays Gaussian, and its three covariances $\Sigma_x$, $C$, $\Sigma_p$ obey a deterministic Riccati equation; the detector efficiency scales the information gain by $\eta$ while the back-action keeps the full rate $\lambda$.

Find the steady state of $\dot\Sigma_x = 2C/m - 4\eta\lambda\Sigma_x^2$, $\dot C = \Sigma_p/m - m\omega^2\Sigma_x - 4\eta\lambda\Sigma_x C$, $\dot\Sigma_p = -2m\omega^2 C + \lambda\hbar^2 - 4\eta\lambda C^2$, keeping the one root that `Resolve` proves has $\Sigma_x, \Sigma_p > 0$ for every positive $m, \lambda, \omega, \hbar$ and $0 < \eta \le 1$:

```wl
riccatiη = {2 cxp/m - 4 η λ vx^2, vp/m - m ω^2 vx - 4 η λ vx cxp, -2 m ω^2 cxp + λ ℏ^2 - 4 η λ cxp^2};
{vxS, cxpS, vpS} = SelectFirst[SolveValues[Thread[riccatiη == 0], {vx, cxp, vp}],
   Resolve[ForAll[{m, λ, ω, ℏ, η}, m > 0 && λ > 0 && ω > 0 && ℏ > 0 && 0 < η <= 1,
      (#[[1]] | #[[3]]) ∈ Reals && #[[1]] > 0 && #[[3]] > 0], Reals] &];
vxS
```
<!-- run time: 0.14 s -->

The purity of a Gaussian state is $\hbar/\left(2\sqrt{\Sigma_x\Sigma_p - C^2}\right)$:

```wl
Simplify[ℏ/(2 Sqrt[vxS vpS - cxpS^2]), m > 0 && λ > 0 && ω > 0 && ℏ > 0 && 0 < η <= 1]
```
<!-- run time: 0.004 s -->

The purity is $\sqrt\eta$, whatever the trap, the mass or the measurement rate. In other words, the light that escapes detection still kicks the particle, and the only trace it leaves on the steady state is this one number.

Write the steady position variance in units of the zero-point variance $\hbar/2m\omega$, with the dimensionless strength $\mu = 2\hbar\lambda/m\omega^2$:

```wl
ratio = FullSimplify[vxS/(ℏ/(2 m ω)) /. λ -> μ m ω^2/(2 ℏ), m > 0 && μ > 0 && ω > 0 && ℏ > 0 && 0 < η <= 1]
```
<!-- run time: 0.06 s -->

Ask for every efficiency and strength at which the conditional state is narrower than the ground state:

```wl
Reduce[ratio < 1 && μ > 0 && 0 < η <= 1, {η, μ}, Reals]
```
<!-- run time: 0.006 s -->

A perfect detector squeezes the position below the zero-point at every strength $\mu > 0$, however gently it watches. An imperfect one does so only above the threshold $\mu^* = 2\sqrt{1-\eta}/\eta^{3/2}$, which falls to zero as $\eta \to 1$ and diverges as $\eta \to 0$: as long as some light is caught, watching hard enough always beats the vacuum.

Cross-check the weak-measurement limit. As $\mu \to 0$ the information gain and the back-action heating vanish together, so the state should settle to a thermal state of purity $\sqrt\eta$, whose position variance is the zero-point variance divided by its purity:

```wl
Asymptotic[ratio, {μ, 0, 2}, Assumptions -> μ > 0 && 0 < η <= 1]
```
<!-- run time: 0.05 s -->

The variance tends to $1/\sqrt\eta$ times the zero-point variance, as for a thermal state of purity $\sqrt\eta$, and at $\eta = 1$ the squeezing starts as $1 - \mu^2/8$.

## A Monitored Qubit Stays in the Bloch Ball: Ito's Lemma and One Reduce

A qubit (a transmon read out by homodyne detection of $\sigma_z$) follows a stochastic Bloch equation: relaxation $\gamma_1$, pure dephasing $\gamma_\phi$, measurement dephasing $\Gamma_d$, a Rabi drive $\Omega$, and a measurement record that kicks the state through an informative rate $\Gamma_{CI}$ and a phase back-action rate $\Gamma_{BA}$. The state is physical only while $r = |\vec r| \le 1$. Does every trajectory stay in the ball, for every detector efficiency $\eta = (\Gamma_{CI} + \Gamma_{BA})/2\Gamma_d \le 1$ and every choice of rates? A simulation checks one set of rates and a finite ensemble. Here Ito's lemma produces the equation for $r^2$, and `Reduce` decides its sign on the sphere.

Define the Ito process with $r^2$ as its output; converting it to standard form applies Ito's lemma and adds $r^2$ as a fourth state with its own drift and noise:

```wl
blochSDE = ItoProcess[{
    {-(γ1/2 + γφ + Γd) x[t], -(γ1/2 + γφ + Γd) y[t] - Ω z[t], γ1 (1 - z[t]) + Ω y[t]},
    {{-Sqrt[ΓCI] x[t] z[t] - Sqrt[ΓBA] y[t]}, {-Sqrt[ΓCI] y[t] z[t] + Sqrt[ΓBA] x[t]}, {Sqrt[ΓCI] (1 - z[t]^2)}},
    x[t]^2 + y[t]^2 + z[t]^2}, {{x, y, z}, {x0, y0, z0}}, t];
{r2Drift, r2Noise} = With[{std = ItoProcess[blochSDE]}, {std["Drift"][[4]], std["Diffusion"][[4, 1]]}];
onSphere = {x[t] -> Sqrt[1 - zs^2] Cos[φs], y[t] -> Sqrt[1 - zs^2] Sin[φs], z[t] -> zs};
{noiseOnSphere, driftOnSphere} = Factor[Simplify[{r2Noise, r2Drift} /. onSphere]]
```
<!-- run time: 0.02 s -->

The noise on $r^2$ vanishes on the sphere, so a measurement kick only slides the state along it. The drift is $(1-z)\left[(1+z)(\Gamma_{CI} + \Gamma_{BA} - 2\Gamma_d - 2\gamma_\phi) - \gamma_1(1-z)\right]$, with no dependence on the drive or the azimuth.

Ask whether any point of the sphere, any efficiency $0 \le \eta \le 1$ and any nonnegative rates make that drift point outward:

```wl
Reduce[(driftOnSphere /. ΓBA -> 2 η Γd - ΓCI) > 0 && -1 <= zs <= 1 && 0 <= η <= 1 &&
   ΓCI >= 0 && 2 η Γd - ΓCI >= 0 && Γd >= 0 && γφ >= 0 && γ1 >= 0, {zs, η, ΓCI, Γd, γφ, γ1}, Reals]
```
<!-- run time: 0.006 s -->

`False`. On the sphere the noise is tangent and the drift points inward or along it, so the exact flow never leaves the ball. In other words, positivity of the conditional density matrix is a property of the equation, not of a particular run.

Cross-check that the bound is sharp: drop $\eta \le 1$ and ask for which efficiencies some point of the sphere has an outward drift:

```wl
Reduce[Exists[zs, -1 < zs < 1, (driftOnSphere /. ΓBA -> 2 η Γd - ΓCI) > 0] && Γd > 0 && γφ >= 0 && γ1 >= 0, η, Reals]
```
<!-- run time: 0.006 s -->

An outward drift appears exactly when $\eta > 1 + \gamma_\phi/\Gamma_d$: the record would carry more information than the measurement dephasing pays for, and the equation no longer describes a density matrix.

## The Collapse Master Equation in Spherical Coordinates: Ito's Lemma Shows That Measurement Purifies

A qubit with Hamiltonian $H = \frac{\omega}{2}\sigma_x$ is measured continuously along $\sigma_z$ at rate $\gamma$. Its state obeys the Ito stochastic master equation $d\rho = \left(-i[H,\rho] + \gamma(\sigma_z\rho\sigma_z - \rho)\right)dt + \sqrt{\gamma}\left(\sigma_z\rho + \rho\sigma_z - 2\langle\sigma_z\rangle\rho\right)dW$. What does the length $r$ of the Bloch vector do along one run? Changing variables in an Ito equation is not the ordinary chain rule: Ito's lemma adds a term from the second derivatives and the noise. A trajectory solver such as QuTiP's `smesolve` only samples runs. Here `ItoProcess` holds the equation as an object, and its `"Derivative"` property applies Ito's lemma to any function of the state.

Write the master equation as an Ito process for the Bloch vector $(x, y, z)$, reading the drift and the noise off $\rho = (1 + \vec r\cdot\vec\sigma)/2$:

```wl
σs = Table[PauliMatrix[k], {k, 3}]; ρb = (IdentityMatrix[2] + {x, y, z} . σs)/2;
bl[m_] := Simplify[Table[Tr[m . s], {s, σs}]]; hx = ω/2 PauliMatrix[1]; sz = PauliMatrix[3];
sme = ItoProcess[{bl[-I (hx . ρb - ρb . hx) + γ (sz . ρb . sz - ρb)],
    List /@ bl[Sqrt[γ] (sz . ρb + ρb . sz - Tr[sz . ρb + ρb . sz] ρb)]}, {{x, y, z}, {x0, y0, z0}}, t]
```
<!-- run time: 0.01 s -->

The drift holds the measurement dephasing $-2\gamma$ of $x$ and $y$ and the rotation by $\omega$ about $x$. The noise vanishes at the poles $z = \pm 1$, the eigenstates of $\sigma_z$.

Apply Ito's lemma to $r = \sqrt{x^2+y^2+z^2}$ and write the result in spherical coordinates $(r, \theta, \varphi)$:

```wl
toSph = Thread[{x[t], y[t], z[t]} -> r {Sin[θ] Cos[φ], Sin[θ] Sin[φ], Cos[θ]}];
dr = FullSimplify[sme["Derivative", Sqrt[x[t]^2 + y[t]^2 + z[t]^2]] /. toSph, r > 0 && 0 < θ < π && γ > 0]
```
<!-- run time: 0.10 s -->

Both the drift and the noise of $r$ carry the factor $r^2 - 1$. At $r = 1$ both vanish, so a pure state stays pure on every run; for $r < 1$ the drift $2\gamma(1 - r^2)\sin^2\theta/r$ is positive, so a mixed state purifies. The Hamiltonian drops out, since a rotation does not change the length. In other words, the purification is the Ito term at work: the ordinary chain rule alone gives $-2\gamma r\sin^2\theta$, the shrinking of the averaged Lindblad state, and the second-derivative term of the noise reverses its sign.

Check against Ito's product rule, $d(r^2) = 2r\,dr + (dr)^2$, applied to the drift and noise of $r$ just found:

```wl
FullSimplify[(sme["Derivative", x[t]^2 + y[t]^2 + z[t]^2][[1, 1]] /. toSph) == 2 r dr[[1, 1]] + dr[[2, 1, 1]]^2, r > 0 && 0 < θ < π && γ > 0]
```
<!-- run time: 0.02 s -->

Get all three spherical equations from the same call, then sample one run that starts mixed at $r = 0.5$, with $\gamma = 1$ and $\omega = 0.3$, and plot $r(t)$:

```wl
sph = Simplify[sme["Derivative", ToSphericalCoordinates[{x[t], y[t], z[t]}]] /. toSph, r > 0 && Sin[θ] > 0 && γ > 0];
ListLinePlot[RandomFunction[ItoProcess[sph /. {ω -> 0.3, γ -> 1}, {{r, θ, φ}, {0.5, π/3, 0}}, t],
   {0, 4, 0.001}]["PathComponent", 1], PlotRange -> {0, 1.05}, AxesLabel -> {"γ t", "r"}]
```
<!-- run time: 0.54 s -->

The radius climbs to the surface of the Bloch ball within a few $1/\gamma$ and stays there, while the angles keep moving under the drive and the record.

## A Hamiltonian with White Noise: Converting Stratonovich to Ito Produces the Dephasing Drift

A qubit has the random Hamiltonian $H = \frac{\omega}{2}\sigma_x + \sqrt{\gamma}\,\xi(t)\,\sigma_z$, with $\xi$ white noise. Because white noise is the limit of smooth noise, the equation of motion is a Stratonovich equation, and each run is unitary. The average over runs, however, is read off the Ito form, whose drift contains an extra term $\frac12 (b\cdot\nabla) b$ from the noise $b$. Here the conversion is done on the symbolic equation.

Write the noise-free drift $-i[H_0, \rho]$ and the Stratonovich noise $-i\sqrt{\gamma}[\sigma_z, \rho]$ in Bloch coordinates:

```wl
pauli = Table[PauliMatrix[k], {k, 3}]; blochOf[m_] := Simplify[Table[Tr[m . s], {s, pauli}]];
ρv = (IdentityMatrix[2] + {x, y, z} . pauli)/2; ham = ω/2 PauliMatrix[1];
{drift0, noiseB} = {blochOf[-I (ham . ρv - ρv . ham)], blochOf[-I Sqrt[γ] (PauliMatrix[3] . ρv - ρv . PauliMatrix[3])]}
```
<!-- run time: 0.001 s -->

The drift has no $\gamma$ in it, and the noise only rotates the Bloch vector about $z$, so each run keeps its length.

Convert the Stratonovich process to Ito form and read its drift:

```wl
ito = ItoProcess[StratonovichProcess[{drift0, List /@ noiseB}, {{x, y, z}, {x0, y0, z0}}, t]];
ito["Drift"]
```
<!-- run time: 0.001 s -->

The conversion adds $-2\gamma x$ and $-2\gamma y$: dephasing at rate $2\gamma$ in the plane perpendicular to the noise axis. In other words, every single run is a pure rotation, yet the runs rotate by different random angles about $z$, and the Ito drift records how their average shrinks.

Check that this drift is the Bloch form of the Lindblad equation $\dot\rho = -i[H_0,\rho] + \gamma(\sigma_z\rho\sigma_z - \rho)$:

```wl
blochOf[-I (ham . ρv - ρv . ham) + γ (PauliMatrix[3] . ρv . PauliMatrix[3] - ρv)] == (ito["Drift"] /. Thread[{x[t], y[t], z[t]} -> {x, y, z}])
```
<!-- run time: 0.001 s -->

## A Laser with a Diffusing Phase: the Exact Mean Field and the Non-Stationary Two-Time Correlation

The field of a laser is $E(t) = \cos(\omega t + \phi_t)$, and its phase diffuses, $d\phi = \sqrt{2\kappa}\,dW$. What are the average field and the correlation $\langle E(t)E(t+\tau)\rangle$? Averaging sampled traces gives these up to sampling error at chosen $\kappa$. Here the field is a `TransformedProcess` of a Wiener process, and its moments come out exactly with $\omega$, $\kappa$, $t$ and $\tau$ as symbols.

Build the field process and take its mean, $\langle E(t)\rangle$:

```wl
field = TransformedProcess[Cos[ω t + ϕ[t]], ϕ \[Distributed] WienerProcess[0, Sqrt[2 κ]], t];
Mean[field[t]]
```
<!-- run time: 0.07 s -->

The mean field oscillates at $\omega$ and decays as $e^{-\kappa t}$, although no trace ever loses amplitude. The decay is pure dephasing among runs: their phases spread with variance $2\kappa t$ and the average cancels.

Compute $\langle E(t)E(t+\tau)\rangle$ for $\tau > 0$ as the covariance function plus the product of the means:

```wl
twoTime = FullSimplify[ComplexExpand[CovarianceFunction[field, t, t + τ] + Mean[field[t]] Mean[field[t + τ]]], t > 0 && τ > 0 && ω ∈ Reals]
```
<!-- run time: 0.83 s -->

The correlation is $\frac12 e^{-\kappa\tau}\cos\omega\tau + \frac12 e^{-\kappa(4t+\tau)}\cos\omega(2t+\tau)$. The first term depends only on the lag. The second depends on $t$ itself, because every run starts from the same phase: it carries the sum $\phi_t + \phi_{t+\tau}$, whose variance $2\kappa(4t+\tau)$ grows with $t$. In other words, the field is not stationary until it has forgotten its initial phase, on the time scale $1/4\kappa$.

Check that at long times only the stationary term survives:

```wl
Limit[twoTime, t -> ∞, Assumptions -> κ > 0 && τ > 0 && ω ∈ Reals]
```
<!-- run time: 0.06 s -->

## Gaussian Kraus Operators of a Weak Measurement: the Commutator of Two Kicks and an Exact Purity Law

A detector that reads a Pauli observable $A$ over a short window, with noise variance $D$, returns a reading $\bar I$ and updates the state by the Kraus operator $M(\bar I) = e^{\bar I A/2D}$, followed by normalization (Korotkov's quantum Bayesian update). With two detectors on $\sigma_z$ and $\sigma_x$, does the order of the two updates matter, and what happens to the purity? A trajectory simulation answers one reading at a time. Here the readings stay symbols.

Compute the commutator of a $\sigma_z$ update with reading $a$ and a $\sigma_x$ update with reading $b$:

```wl
kraus[ib_, A_] := MatrixExp[ib/(2 var) A];
FullSimplify[kraus[a, PauliMatrix[3]] . kraus[b, PauliMatrix[1]] - kraus[b, PauliMatrix[1]] . kraus[a, PauliMatrix[3]], {a, b} ∈ Reals && var > 0]
```
<!-- run time: 0.04 s -->

The commutator is $2i\sinh(a/2D)\sinh(b/2D)\,\sigma_y$ for every pair of readings. It inherits $[\sigma_z, \sigma_x] = 2i\sigma_y$ from the observables, and it is first order in each reading, so the two orders agree only as the window shrinks.

Take one full step, a half $\sigma_z$ update on either side of an update along an axis at angle $\vartheta$ in the $xz$-plane, and compute its determinant:

```wl
mStep = kraus[a/2, PauliMatrix[3]] . kraus[b, Cos[ϑ] PauliMatrix[3] + Sin[ϑ] PauliMatrix[1]] . kraus[a/2, PauliMatrix[3]];
FullSimplify[Det[mStep], {a, b, ϑ} ∈ Reals && var > 0]
```
<!-- run time: 0.006 s -->

Each factor is the exponential of a traceless matrix, so the step is unimodular for every reading and every angle. Since the determinant is multiplicative and $\det\rho = (1 - r^2)/4$, the normalized update gives $1 - r'^2 = (1 - r^2)/(\mathrm{Tr}\,M\rho M)^2$ exactly. In other words, a pure state stays pure (as it must under any single Kraus operator), and for a mixed state each step divides the mixedness $1 - r^2$ by the square of the weight $\mathrm{Tr}\,M\rho M$ of the reading, whatever the readings and the axis.

Check that the readings form a POVM: weighted by $w(\bar I) = \mathcal N(0, D)\,e^{-1/2D}$, the operators $M(\bar I)^2$ integrate to the identity:

```wl
Integrate[PDF[NormalDistribution[0, Sqrt[var]], ib] Exp[-1/(2 var)] kraus[ib, PauliMatrix[3]] . kraus[ib, PauliMatrix[3]], {ib, -∞, ∞}, Assumptions -> var > 0]
```
<!-- run time: 0.60 s -->

## The Quantum Kalman-Bucy Filter as an ItoProcess: Purity Is a Noiseless Attractor

A free particle of mass $m$ is watched continuously in position at strength $\lambda$. A Gaussian conditional state stays Gaussian, so five numbers describe it: the means $\langle x\rangle, \langle p\rangle$ and the covariances $V_x, C, V_p$. Only the means feel the noise; the covariances follow a Riccati equation (the quantum Kalman-Bucy filter). The purity of a Gaussian state is $\hbar/(2\sqrt{V_xV_p - C^2})$. Does a run keep it? A general-purpose integrator of the wavefunction equation drifts in norm and must renormalize at every step. Here the filter is an `ItoProcess`, and Ito's lemma answers for $V_xV_p - C^2$ with every parameter symbolic.

Define the filter and apply Ito's lemma to $P = V_xV_p - C^2$:

```wl
kb = ItoProcess[{{mp/mass, 0, 2 cv/mass - 4 λ vx^2, vp/mass - 4 λ vx cv, λ ℏ^2 - 4 λ cv^2},
    {{2 Sqrt[λ] vx}, {2 Sqrt[λ] cv}, {0}, {0}, {0}}}, {{mx, mp, vx, cv, vp}, {0, 0, vx0, 0, vp0}}, t];
Factor[kb["Derivative", vx[t] vp[t] - cv[t]^2]]
```
<!-- run time: 0.001 s -->

The noise term of $P$ is zero and its drift is $-4\lambda V_x\,(P - \hbar^2/4)$. The record never touches the purity, and the pure value $P = \hbar^2/4$ is a fixed point. Because $V_x > 0$, it attracts: a mixed Gaussian purifies at the rate $4\lambda V_x$, the same on every run, and numerical error in $P$ is damped rather than accumulated.

Check with one sampled run that starts mixed ($V_x = V_p = 1$, $C = 0$, with $\hbar = m = \lambda = 1$) and read $P$ at the end, with no renormalization step:

```wl
With[{v = Last[RandomFunction[kb /. {mass -> 1, λ -> 1, ℏ -> 1, vx0 -> 1, vp0 -> 1}, {0, 6, 0.01}]["Values"]]}, v[[3]] v[[5]] - v[[4]]^2]
```
<!-- run time: 0.002 s -->

## A Gaussian Pulse on a Qubit: the Propagator Is a Rotation by the Pulse Area

Drive a qubit on resonance with a Gaussian envelope, $H(t) = \tfrac12\Omega_0 e^{-(t-t_0)^2/2\tau^2}\sigma_x$. Which gate has it made by time $t$? A numerical integrator returns one $2\times 2$ matrix for each choice of $\Omega_0$, $t_0$, $\tau$ and $t$. Here the unknown is the whole propagator, a $2\times 2$ matrix function, and the equation $\dot U = -iH(t)U$ is solved once with every parameter symbolic.

Solve $\dot U = -\tfrac{i}{2}\Omega_0 e^{-(t-t_0)^2/2\tau^2}\sigma_x U$ with $U(0) = \mathbb{1}$, declaring $U(t)$ a $2\times 2$ matrix:

```wl
upulse = DSolveValue[{u'[t] == -(I/2) Ω0 Exp[-(t - t0)^2/(2 τ^2)] PauliMatrix[1] . u[t], u[0] == IdentityMatrix[2]}, u[t] ∈ Matrices[{2, 2}], t]
```
<!-- run time: 0.20 s -->

The entries are exponentials of error functions. The error function appears because it is the running integral of the Gaussian, which suggests that the envelope enters only through the area swept out so far.

Compute that area, $\theta(t) = \int_0^t \Omega(s)\,ds$:

```wl
area = Integrate[Ω0 Exp[-(s - t0)^2/(2 τ^2)], {s, 0, t}, Assumptions -> τ > 0]
```
<!-- run time: 0.66 s -->

Test whether the solved propagator is exactly the rotation $e^{-i\theta(t)\sigma_x/2}$:

```wl
FullSimplify[upulse == MatrixExp[-I area/2 PauliMatrix[1]], τ > 0]
```
<!-- run time: 0.09 s -->

The shape of the pulse has washed out: only its area survives, so a Gaussian and a square pulse of equal area make the same gate. In other words, this is the area theorem, and it holds because $H(t)$ points along $\sigma_x$ at every instant, so the Hamiltonians at different times commute and the time-ordered exponential becomes an ordinary one. Add a detuning term $\Delta\sigma_z$ and that commutation fails, and the shape matters again.

Compare with a numerical integration at $\Omega_0 = 3$, $t_0 = 1$, $\tau = 1/4$, $t = 2$:

```wl
Max[Abs[(upulse /. {Ω0 -> 3, t0 -> 1, τ -> 1/4, t -> 2}) - NDSolveValue[{v'[s] == -(I/2) 3 Exp[-(s - 1)^2/(2 (1/4)^2)] PauliMatrix[1] . v[s], v[0] == IdentityMatrix[2]}, v[2], {s, 0, 2}]]]
```
<!-- run time: 0.007 s -->

The difference sits at the integrator's tolerance.

## Ehrenfest and Mandelstam-Tamm for Any H(t): Identities Without Solving for the State

For a qubit with a completely arbitrary Hamiltonian $H(t) = h_0(t)\mathbb{1} + h_x(t)\sigma_x + h_y(t)\sigma_y + h_z(t)\sigma_z$, does $\tfrac{d}{dt}\langle A\rangle = i\langle[H(t), A]\rangle$ hold, and does the state move through Hilbert space at speed $\Delta E(t)$? A simulator can only check these on one trajectory for one chosen drive. Here the four drive functions stay undetermined, the state is never solved for, and the Schrödinger equation enters only as a rule for the time derivatives ($\hbar = 1$).

Write the ket as $(c_1(t), c_2(t))$ and the bra as $(b_1(t), b_2(t))$, with $\dot{|\psi\rangle} = -iH|\psi\rangle$ and $\dot{\langle\psi|} = i\langle\psi|H$. Then differentiate $\langle A\rangle$ for an arbitrary constant matrix $A$ and subtract $i\langle[H, A]\rangle$:

```wl
hE = h0[t] IdentityMatrix[2] + hx[t] PauliMatrix[1] + hy[t] PauliMatrix[2] + hz[t] PauliMatrix[3];
ket = {c1[t], c2[t]}; bra = {b1[t], b2[t]};
schr = Join[Thread[D[ket, t] -> -I hE . ket], Thread[D[bra, t] -> I bra . hE]];
obsA = Array[a, {2, 2}];
Simplify[(D[bra . obsA . ket, t] /. schr) - I bra . (hE . obsA - obsA . hE) . ket]
```
<!-- run time: 0.0007 s -->

The residual is zero for every drive and every observable. In other words, Ehrenfest's theorem is a consequence of the equation of motion alone, not of any particular solution, so it holds even for drives whose propagator has no closed form.

Now the Fubini-Study speed squared, $\langle\dot\psi|\dot\psi\rangle - |\langle\psi|\dot\psi\rangle|^2$, against the energy variance $\langle H^2\rangle - \langle H\rangle^2$:

```wl
speed2 = (D[bra, t] . D[ket, t] - (bra . D[ket, t]) (D[bra, t] . ket)) /. schr;
Simplify[speed2 - (bra . hE . hE . ket - (bra . hE . ket)^2)]
```
<!-- run time: 0.0004 s -->

The state moves through projective Hilbert space at exactly $\Delta E(t)$ at every instant (Anandan and Aharonov). The Mandelstam-Tamm bound follows as geometry: the distance between start and end, $\arccos|\langle\psi_0|\psi_t\rangle|$, cannot exceed the length of the path, $\int_0^t \Delta E\,dt'$.

The identity part $h_0(t)$ only adds a global phase, so it must not move the state. On normalized states the variance does not depend on it:

```wl
Simplify[D[bra . hE . hE . ket - (bra . hE . ket)^2, h0[t]], bra . ket == 1]
```
<!-- run time: 0.003 s -->

## The Jaynes-Cummings Propagator: Exact for Every Photon Number n

A two-level atom couples to one cavity mode, $H = \Delta\sigma_z/2 + g(\sigma_+ a + \sigma_- a^\dagger)$ in the frame rotating at the cavity frequency. The interaction only exchanges one excitation between the atom and the field, so $|e, n\rangle$ couples to $|g, n+1\rangle$ and nothing else. What is the propagator for every photon number? A Fock-space solver truncates the field at some cutoff and answers for each $n$ below it. Here $n$, $g$, $\Delta$ and $t$ stay symbols.

Exponentiate the $2\times 2$ block on $\{|e,n\rangle, |g,n+1\rangle\}$, whose off-diagonal element is $\langle e,n|\sigma_+ a|g,n+1\rangle = g\sqrt{n+1}$:

```wl
jc = {{Δ/2, g Sqrt[n + 1]}, {g Sqrt[n + 1], -Δ/2}};
ujc = FullSimplify[MatrixExp[-I t jc], g > 0 && n >= 0 && Δ ∈ Reals]
```
<!-- run time: 0.66 s -->

The propagator rotates at the generalized Rabi frequency $\Omega_n = \sqrt{\Delta^2 + 4g^2(n+1)}$. The probability that an atom starting in $|e, n\rangle$ has emitted its photon is $|U_{21}|^2$:

```wl
Simplify[Abs[ujc[[2, 1]]]^2, g > 0 && n >= 0 && t > 0 && Δ ∈ Reals]
```
<!-- run time: 0.05 s -->

On resonance the atom empties completely, at a frequency that grows as $\sqrt{n+1}$; detuned, it never empties, with a ceiling set by $4g^2(n+1)/\Omega_n^2$. The $+1$ under the square root is spontaneous emission: even at $n = 0$ the vacuum makes the atom oscillate.

The dressed states split by the difference of the block's eigenvalues; on resonance:

```wl
Simplify[Differences[Eigenvalues[jc /. Δ -> 0]], g > 0 && n >= 0]
```
<!-- run time: 0.002 s -->

The splitting is $2g\sqrt{n+1}$: the vacuum Rabi splitting at $n = 0$, and in general a splitting that grows as the square root of the photon number, which is the graininess of the field that a classical drive does not have.

Check against the full QuantumFramework Hamiltonian $g(\sigma_+ a + \sigma_- a^\dagger)$ with the field truncated at 6 levels: its propagator, restricted to $\{|e,2\rangle, |g,3\rangle\}$, should be the block at $n = 2$:

```wl
hq = g QuantumOperator[{{0, 1}, {0, 0}}, {1}] @ AnnihilationOperator[6, {2}];
ufull = MatrixExp[-I t Normal[(hq + hq["Dagger"])["Matrix"]]];
FullSimplify[ufull[[{3, 10}, {3, 10}]] == (ujc /. {Δ -> 0, n -> 2}), g > 0]
```
<!-- run time: 0.13 s -->

## Syndrome Populations of a Repetition Code: the Closed Equation from a General ρ

Three qubits carry the bit-flip code, with stabilizers $Z_1Z_2$ and $Z_2Z_3$, and each qubit flips at rate $\gamma$. The four syndrome sectors have populations $p_j = \mathrm{Tr}[\Pi_j\rho]$. Do these four numbers obey a closed equation of their own, and what is it? A simulator would integrate the master equation for some initial states and fit the populations. Here the question is put to a completely general $8\times 8$ matrix $\rho$ with 64 free entries, and the answer is solved for rather than proposed.

Set up the stabilizers, the four sector projectors $\Pi_{(s_1,s_2)} = \tfrac14(1 + s_1Z_1Z_2)(1 + s_2Z_2Z_3)$ and a general $\rho$, and read off the rank of each projector:

```wl
{id2, sx, sz} = PauliMatrix /@ {0, 1, 3};
op[o_, q_] := KroneckerProduct @@ ReplacePart[{id2, id2, id2}, q -> o];
stabs = {op[sz, 1] . op[sz, 2], op[sz, 2] . op[sz, 3]};
proj[{s1_, s2_}] := (IdentityMatrix[8] + s1 stabs[[1]]) . (IdentityMatrix[8] + s2 stabs[[2]])/4;
sectors = {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}}; rho = Array[r, {8, 8}];
Tr[proj[#]] & /@ sectors
```
<!-- run time: 0.0005 s -->

For each sector, take the trace of $\Pi_j$ against the bit-flip dissipator $\sum_k (X_k\rho X_k - \rho)$, and ask `SolveAlways` for coefficients $c_s$ that make it equal to $\sum_s c_s p_s$ for every one of the 64 entries of $\rho$:

```wl
inPops[tr_] := (c /@ sectors) . (p /@ sectors) /. First[SolveAlways[tr == (c /@ sectors) . (Tr[proj[#] . rho] & /@ sectors), Flatten[rho]]];
drift = (γ inPops[Sum[Tr[proj[#] . (op[sx, k] . rho . op[sx, k] - rho)], {k, 3}]] &) /@ sectors
```
<!-- run time: 0.01 s -->

A solution exists, and that existence is the closure: 64 numbers go in, and the answer depends on only four combinations of them. Each sector gains weight at rate $\gamma$ from each of the other three and loses it at $3\gamma$, because each single flip carries a sector onto exactly one other.

Add back $4\gamma p_j$ to each entry:

```wl
Simplify[drift + 4 γ (p /@ sectors)]
```
<!-- run time: 0.002 s -->

Every entry is $\gamma$ times the total population, so with $\sum_s p_s = 1$ the equation is $\dot p_j = \gamma(1 - 4p_j)$: every sector relaxes to $1/4$ at rate $4\gamma$.

The closure is a property of the dynamics, not of the method. Add a Hamiltonian term that does not commute with the stabilizers, $-i[X_1, \rho]$, and ask the same question:

```wl
SolveAlways[Tr[proj[{1, 1}] . (-I (op[sx, 1] . rho - rho . op[sx, 1]))] == (c /@ sectors) . (Tr[proj[#] . rho] & /@ sectors), Flatten[rho]]
```
<!-- run time: 0.002 s -->

No coefficients exist: the populations now couple to coherences between sectors, and four numbers are no longer enough.

## Continuous Position Measurement: the Conditional Variances Carry No Noise

A free particle is watched continuously, so its conditional state obeys a stochastic equation in which an observable $O$ picks up a noise term $\sqrt{\lambda}\,\langle\{x - \langle x\rangle, O\}\rangle\,dW$. The means $\langle x\rangle$ and $\langle p\rangle$ wander randomly. Do the variances $V_x$, $V_p$ and the covariance $C$ wander too? A trajectory simulation would answer by sampling many runs. Here the canonical commutator is handed to the kernel's non-commutative algebra, a Gaussian state supplies the moments, and the noise coefficients are computed exactly.

Declare the algebra generated by $x$ and $p$ with $px = xp - i\hbar$; `NonCommutativeExpand` then writes every product with $x$ to the left of $p$:

```wl
xpAlg = NonCommutativeAlgebra[<|"Generators" -> {{xo, po}}, "CommutationRelations" -> {po ** xo - xo ** po + I ℏ}|>];
NonCommutativeExpand[po ** po ** xo, xpAlg]
```
<!-- run time: 0.0001 s -->

For a Gaussian state with means $\mu_x, \mu_p$, variances $V_x, V_p$ and symmetric covariance $C$, the ordered generating function is $\langle e^{\alpha x}e^{\beta p}\rangle = \exp\!\big(\alpha\mu_x + \beta\mu_p + \tfrac12(\alpha^2V_x + 2\alpha\beta C + \beta^2V_p) + \tfrac{i\hbar}{2}\alpha\beta\big)$, so each ordered monomial $x^ap^b$ has expectation $\partial_\alpha^a\partial_\beta^b$ of it at zero. Build that expectation and test it on $\langle xp\rangle$ and $\langle px\rangle$:

```wl
gauss = Exp[α mx + β mp + (α^2 Vx + 2 α β Cc + β^2 Vp)/2 + I ℏ α β/2];
mean[o_] := With[{q = NonCommutativeExpand[o, xpAlg] /. {GeneralizedPower[NonCommutativeMultiply, g_, k_] :> g^k, NonCommutativeMultiply -> Times}},
   Total[(#2 D[gauss, {α, #1[[1]]}, {β, #1[[2]]}] /. {α -> 0, β -> 0}) & @@@ CoefficientRules[q, {xo, po}]]];
{mean[xo ** po], mean[po ** xo]}
```
<!-- run time: 0.0006 s -->

The two orderings differ by exactly $i\hbar$, the commutator, while their average is the classical covariance.

Now the $dW$ coefficient of each second moment, with the Itô correction for the products of means: $dV_x$ takes the noise of $\langle x^2\rangle$ minus $2\mu_x$ times the noise of $\langle x\rangle$, and likewise for $C$ and $V_p$:

```wl
noise[o_] := Sqrt[λ] (mean[xo ** o + o ** xo] - 2 mx mean[o]);
Simplify[{noise[xo ** xo] - 2 mx noise[xo], noise[(xo ** po + po ** xo)/2] - mx noise[po] - mp noise[xo], noise[po ** po] - 2 mp noise[po]}]
```
<!-- run time: 0.002 s -->

All three noise coefficients vanish. In other words, every run puts the wave packet in a different place, but all runs give it the same width at the same time: the conditional variances obey ordinary differential equations. The noise in $dV_x$ is proportional to the third cumulant of $x$, and a Gaussian state has none. The same holds for $C$ and $V_p$: every combination that enters the noise is symmetrically ordered, so the commutator drops out and only the third cumulants, all zero, remain.

A Gaussian's third moment has no content beyond its first two, $\langle x^3\rangle = \mu_x^3 + 3\mu_xV_x$:

```wl
mean[xo ** xo ** xo]
```
<!-- run time: 0.0005 s -->

## Symmetry Proofs at Symbolic Dimension: Tensor Identities with No Components

An antisymmetric tensor contracted with a symmetric one gives zero; the Ricci tensor is symmetric; the curl of a gradient vanishes. Each follows from index symmetries alone. A numerical library can check them only on one random instance in one fixed dimension. Here the tensors are declared by their symmetries in dimension $d$, never given entries, and `TensorReduce` brings each expression to a canonical form.

Contract an antisymmetric $A$ with a symmetric $S$, $\sum_{ij}A_{ij}S_{ij}$, in dimension $d$:

```wl
TensorReduce[TensorContract[TensorProduct[aT, sT], {{1, 3}, {2, 4}}], Assumptions -> {aT ∈ Arrays[{d, d}, Reals, Antisymmetric[{1, 2}]], sT ∈ Arrays[{d, d}, Reals, Symmetric[{1, 2}]]}]
```
<!-- run time: 0.001 s -->

Declare the Riemann tensor by its three single-term symmetries (antisymmetry in each pair, symmetry under pair exchange) and a symmetric metric, form the Ricci tensor $R_{bd} = R_{abcd}\,g^{ac}$, and compare it with its transpose:

```wl
riemSym = {{Cycles[{{1, 2}}], -1}, {Cycles[{{3, 4}}], -1}, {Cycles[{{1, 3}, {2, 4}}], 1}};
asmR = {rT ∈ Arrays[{d, d, d, d}, Reals, riemSym], gT ∈ Arrays[{d, d}, Reals, Symmetric[{1, 2}]]};
ricci = TensorContract[TensorProduct[rT, gT], {{1, 5}, {3, 6}}];
TensorReduce[ricci, Assumptions -> asmR] === TensorReduce[TensorTranspose[ricci, {2, 1}], Assumptions -> asmR]
```
<!-- run time: 0.002 s -->

The Ricci tensor is symmetric in every dimension, which is what lets it sit on the left of Einstein's equations next to the symmetric stress tensor.

Contract the Levi-Civita tensor with the Hessian $H_{jk} = \partial_j\partial_k f$, $(\nabla\times\nabla f)_i = \varepsilon_{ijk}H_{jk}$:

```wl
TensorReduce[TensorContract[TensorProduct[epsT, hT], {{2, 4}, {3, 5}}], Assumptions -> {epsT ∈ Arrays[{3, 3, 3}, Reals, Antisymmetric[{1, 2, 3}]], hT ∈ Arrays[{3, 3}, Reals, Symmetric[{1, 2}]]}]
```
<!-- run time: 0.0007 s -->

The result is a zero vector of length 3, not the scalar 0: the canonical form keeps the shape, so it can sit inside a larger expression.

The canonicalizer does not answer 0 indiscriminately. Contract an antisymmetric tensor with itself, $\sum_{ij}A_{ij}A_{ij}$, which is a sum of squares:

```wl
TensorReduce[TensorContract[TensorProduct[aT, aT], {{1, 3}, {2, 4}}], Assumptions -> {aT ∈ Arrays[{d, d}, Reals, Antisymmetric[{1, 2}]]}]
```
<!-- run time: 0.0006 s -->

## The Riemann Identities Derived from the Young Projector, First Bianchi Included

The Riemann tensor has three single-term symmetries and one identity of a different kind, the first Bianchi identity $R_{abcd} + R_{acdb} + R_{adbc} = 0$, which sets a sum of three permuted copies to zero. `TensorReduce` works only with symmetries of the form $R = \pm R^\sigma$, so it cannot reach Bianchi from the declared pair symmetries:

```wl
TensorReduce[rB + TensorTranspose[rB, {1, 3, 4, 2}] + TensorTranspose[rB, {1, 4, 2, 3}], Assumptions -> rB ∈ Arrays[{d, d, d, d}, Reals, {{Cycles[{{1, 2}}], -1}, {Cycles[{{3, 4}}], -1}, {Cycles[{{1, 3}, {2, 4}}], 1}}]]
```
<!-- run time: 0.0007 s -->

Can all four identities instead be found, without listing them, from the one fact that the Riemann tensor lives in the $\{2,2\}$ irreducible representation of $S_4$? Project a random integer tensor onto that representation with the TensorNetworks Young projector, and collect its 24 slot-permuted copies:

```wl
Needs["Wolfram`TensorNetworks`Symmetry`"];
SeedRandom[7];
riem = YoungProject[RandomInteger[{-9, 9}, {4, 4, 4, 4}], YoungTableau[{{1, 3}, {2, 4}}]];
perms4 = Permutations[Range[4]];
copies = Flatten[Transpose[riem, #]] & /@ perms4;
MatrixRank[copies]
```
<!-- run time: 0.006 s -->

The 24 copies span only a 2-dimensional space, the dimension of the $\{2,2\}$ representation, so they obey 22 independent linear relations. Every element $\sum_\sigma c_\sigma\sigma$ of the group algebra that annihilates the tensor is an identity of the Riemann tensor. Search the three-term sums over permutations that fix slot 1:

```wl
Select[Subsets[Flatten[Position[perms4, {1, _, _, _}]], {3}], Total[copies[[#]]] == 0 Range[256] &] /. i_Integer :> perms4[[i]]
```
<!-- run time: 0.001 s -->

The first relation found is the cyclic sum over slots 2, 3, 4: the first Bianchi identity, produced by the projector and not supplied. The second is the same identity after the antisymmetry in the last pair. A search over two-term sums $\sigma_a \pm \sigma_b$ returns the pair antisymmetries and the pair exchange in the same way.

Count the components: the image of the $\{2,2\}$ projector on all $4^4$ tensors should have the textbook dimension of the Riemann tensor in four dimensions, $d^2(d^2-1)/12 = 20$, one less than the 21 that the pair symmetries alone allow:

```wl
MatrixRank[Flatten[YoungProject[ArrayReshape[#, {4, 4, 4, 4}], YoungTableau[{{1, 3}, {2, 4}}]]] & /@ IdentityMatrix[256]]
```
<!-- run time: 0.18 s -->

## Spectrum Estimation from n Copies: Schur-Sampling Probabilities in Closed Form

Given $n$ copies of an unknown qubit state $\rho$, measure which irreducible representation of $S_n$ (which Young diagram $\lambda$) the joint state lies in. The outcome probability is $\Pr(\lambda) = \mathrm{Tr}[P_\lambda\,\rho^{\otimes n}]$, and Keyl and Werner showed that $\lambda/n$ concentrates on the spectrum of $\rho$. What are these probabilities for an arbitrary mixed state? A numerical code would evaluate them for one $\rho$ at a time. Here the Bloch vector stays symbolic.

Build the orthogonal projector onto each $\lambda$-block of four qubits from the Young projectors of its standard tableaux, $P_\lambda = MM^+$ with $M = \sum_T P_T$, and check that the three blocks are Hermitian, idempotent and complete:

```wl
basis4 = ArrayReshape[#, {2, 2, 2, 2}] & /@ IdentityMatrix[16];
isoProj[λ_] := With[{tabs = Select[YoungTableau[TakeList[#, λ]] & /@ Permutations[Range[4]], StandardTableauQ]},
   With[{m = Transpose[Function[ket, Flatten[Total[YoungProject[ket, #] & /@ tabs]]] /@ basis4]}, m . PseudoInverse[m]]];
shapes = {{4}, {3, 1}, {2, 2}}; projs = isoProj /@ shapes;
{Tr /@ projs, Total[projs] == IdentityMatrix[16], And @@ (# . # == # && # == ConjugateTranspose[#] & /@ projs)}
```
<!-- run time: 0.07 s -->

The block dimensions are $(2J+1)$ times the multiplicity of total spin $J$: the shape $\lambda = (\lambda_1, \lambda_2)$ is total spin $J = (\lambda_1 - \lambda_2)/2$.

Take $\rho = \tfrac12(\mathbb{1} + \vec b\cdot\vec\sigma)$ with a symbolic Bloch vector, and compute $\Pr(\lambda) = \mathrm{Tr}[P_\lambda\rho^{\otimes 4}]$:

```wl
ρq = (IdentityMatrix[2] + bx PauliMatrix[1] + by PauliMatrix[2] + bz PauliMatrix[3])/2;
schur = Factor[Tr[# . KroneckerProduct @@ ConstantArray[ρq, 4]] & /@ projs]
```
<!-- run time: 0.009 s -->

The probabilities depend on $\vec b$ only through $|\vec b|^2$, that is, only on the eigenvalues $(1 \pm |\vec b|)/2$: the Young diagram carries the spectrum and nothing about the eigenbasis. A pure state ($|\vec b| = 1$) puts all its weight on the symmetric diagram $\{4\}$, so the measurement says nothing beyond "pure"; it is the mixedness that spreads weight onto $\{3,1\}$ and $\{2,2\}$. The maximally mixed state gives the block dimensions divided by 16.

Compare with the closed form $\Pr(\lambda) = f^\lambda s_\lambda(a, b)$, with $f^\lambda$ the dimension of the $S_4$ representation and $s_\lambda(a,b) = (ab)^{\lambda_2}\,(a^{\lambda_1-\lambda_2+1} - b^{\lambda_1-\lambda_2+1})/(a - b)$ the Schur polynomial of the eigenvalues:

```wl
With[{a = (1 + rb)/2, b = (1 - rb)/2}, Simplify[(schur /. bz -> Sqrt[rb^2 - bx^2 - by^2]) - (TableauDimension[DeleteCases[{##}, 0]] (a b)^#2 (a^(#1 - #2 + 1) - b^(#1 - #2 + 1))/(a - b) & @@@ PadRight[shapes])]]
```
<!-- run time: 0.0006 s -->

## A Stabilizer Measurement Kept Symbolic: Random Outcomes as Bits, Determined Outcomes as Functions of Them

Prepare two qubits in $|+\rangle$ and copy their parity onto a third with two CNOTs, then measure all three in the $Z$ basis. Two outcomes are fair coins; the third is fixed by them. A stabilizer simulator such as Stim samples shots, and the dependence shows up only as a pattern in the samples. Here the measurement is not sampled at all: each random outcome becomes a formal bit, and each determined outcome is stored as a function of those bits.

Prepare the state and read its stabilizer generators:

```wl
parity = PauliStabilizer[3][{"H" -> 1, "H" -> 2, "CNOT" -> {1, 3}, "CNOT" -> {2, 3}}];
parity["Stabilizers"]
```
<!-- run time: 0.44 s -->

The generator $Z_1 Z_2 Z_3$ stabilizes the state, so the three $Z$ outcomes cannot be independent.

Measure all three qubits symbolically, and read the generators and the stored outcomes:

```wl
ms = parity["SymbolicMeasure", {1, 2, 3}];
{ms["Stabilizers"], ms["Outcomes"]}
```
<!-- run time: 0.002 s -->

Only two formal bits $s_1, s_2$ are allocated, one per genuinely random qubit, and each enters a generator as the sign $(-1)^{s_k} = 1 - 2 s_k$. The third outcome is not a new coin: it is stored as a Boolean expression in $s_1$ and $s_2$. In other words, the measured state is still one object, and the correlation the conserved $Z_1Z_2Z_3$ forces is held as a formula.

Tabulate the stored third outcome over the four values of $(s_1, s_2)$:

```wl
bits = Union @ Cases[ms["Phase"], \[FormalS][_], Infinity];
Flatten[Table[{a, b, First[Values[ms["Outcomes"]]] /. Thread[bits -> {a, b}]}, {a, 0, 1}, {b, 0, 1}], 1]
```
<!-- run time: 0.0002 s -->

The third bit is $s_1 \oplus s_2$ in every row, the syndrome of the simplest parity check.

Cross-check with Stim on the same circuit: in every one of 1000 sampled shots the third bit equals the XOR of the first two:

```wl
ExternalEvaluate["Python", "import stim\nc = stim.Circuit('H 0 1\\nCX 0 2 1 2\\nM 0 1 2')\ns = c.compile_sampler().sample(1000)\nbool(((s[:, 0] ^ s[:, 1]) == s[:, 2]).all())"]
```
<!-- run time: 2.8 s, mostly starting Python -->

## A Repetition Code under a QuantumFramework Channel: the Logical Error Rate Is a Polynomial

Encode one bit in three, let each physical qubit flip with probability $p$, and decode by majority vote. What is the probability that the logical bit is wrong, as a function of $p$? A Monte Carlo decoder estimates it at a list of $p$ values. The QEC layer under development beside QuantumFramework (QECCore, not yet in the paclet) returns it as a polynomial, with the noise given as a QuantumFramework channel.

Load the QEC layer from its source folder:

```wl
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/QEC/QECCore/QECCore.wl"];
MemberQ[$ContextPath, "Wolfram`QuantumFramework`QEC`"]
```
<!-- run time: 2.9 s, loading -->

Ask for the logical error rate of the three-qubit repetition code under `QuantumChannel["BitFlip"[p]]`:

```wl
rate[n_] := Expand[QECLogicalErrorRate[QECCode["Repetition", n], QuantumChannel["BitFlip"[p]]]];
rate[3]
```
<!-- run time: 0.03 s -->

The rate is $3p^2 - 2p^3 = 3p^2(1-p) + p^3$: the decoder fails exactly when two or three qubits flip. The linear term is absent, so for small $p$ the code beats a bare qubit.

Find where the five-qubit code beats the three-qubit code:

```wl
Reduce[rate[5] < rate[3] && 0 < p < 1, p]
```
<!-- run time: 0.03 s -->

Below $p = 1/2$ the longer code is better; above it, worse. In other words, the threshold of the repetition code under bit flips is exactly one half, the crossing of two polynomials rather than of two sampled curves.

Cross-check from QuantumFramework primitives alone: apply the bit-flip channel to each qubit of $|000\rangle$ and add the probabilities of the outcomes with two or more ones:

```wl
noisy = Fold[QuantumChannel["BitFlip"[p], {#2}][#1] &, QuantumState["000"], {1, 2, 3}];
Simplify[Total[Pick[noisy["ProbabilitiesList"], Total /@ Tuples[{0, 1}, 3], 2 | 3]] == rate[3], 0 < p < 1]
```
<!-- run time: 0.84 s -->

## A Drive Pulse as an Expression: a DRAG Sweep by Substitution and a Program Edited by Pattern

A transmon is driven by a microwave pulse, and DRAG adds a second quadrature proportional to the time derivative of the first, with a coefficient $\beta$ found by sweeping it. In a numeric pulse SDK the pulse is an array of samples and the sweep is a loop that regenerates the array. With QUALink the pulse stays an expression in $t$ until the compiler samples it, each sweep value is one substitution into that expression, and the compiled program is itself an expression that the ordinary pattern tools query and edit. Nothing below touches a device.

Load QUALink (TransmonLab, its optional physics backend, goes on the path too):

```wl
PacletDirectoryLoad[{"/Users/mohammadb/Documents/GitHub/mbahram/TransmonLab/TransmonLab", "/Users/mohammadb/Documents/GitHub/mbahram/QUALink/QUALink"}];
Needs["WolframQuantumComputation`QUALink`"]
```
<!-- run time: 0.22 s, loading -->

Write the in-phase envelope as a Gaussian and the quadrature as $-\beta$ times its derivative:

```wl
envelope = Exp[-(t - 20)^2/50];
drag = {envelope, -β D[envelope, t]}
```
<!-- run time: 0.0001 s -->

Compile a program that sweeps $\beta$ over three values, and ask its intermediate representation where the drive is played:

```wl
dragProg = QUAProgram["Qubit" -> {5.105*^9, drag, {t, 0, 40}}, "Resonator" -> 6.012*^9, "β" -> {-1, 0, 1}, "Averages" -> 100];
Cases[dragProg["IR"], KeyValuePattern["Head" -> "SwitchPlay"], Infinity]
```
<!-- run time: 0.004 s -->

The compiler found $\beta$ free in the envelope, substituted each value into the expression, sampled one pulse per value, and turned the sweep into a single switch over those pulses inside the program. The program is an association, so the question "where is the drive played" is one `Cases` call rather than a search through generated text.

Check that the samples are the expression we wrote: at $\beta = 1$ the ratio of quadrature to in-phase samples must be $(t-20)/25$ on the 1 ns grid:

```wl
wf = dragProg["Waveforms"];
Max[Abs[wf["qubit_drive_wf_3_q"]/wf["qubit_drive_wf_3"] - Table[(t - 20)/25, {t, 0, 39}]]] < 10^-12
```
<!-- run time: 0.0001 s -->

Once the sweep has picked $\beta = 1$, pin the program to that pulse by rewriting the switch into a single play with `MapAt`:

```wl
pinned = MapAt[<|"Head" -> "Play", "Pulse" -> Last[#Pulses], "Element" -> #Element, "AmpExpr" -> None, "DurationExpr" -> None, "Chirp" -> None|> &,
   dragProg["IR"], Position[dragProg["IR"], KeyValuePattern["Head" -> "SwitchPlay"]]];
Cases[pinned, KeyValuePattern["Head" -> "Play" | "SwitchPlay"], Infinity]
```
<!-- run time: 0.0001 s -->

Every point of the loop now plays the $\beta = 1$ pulse, the readout and the averaging loop are untouched, and the original program is unchanged, since `MapAt` returns a new expression.

## The CZ Gate from Two Transmons: the {11, 02} Block Derived from the Duffing Hamiltonian

A flux-tunable transmon pair makes its CZ gate by bringing $|11\rangle$ into resonance with the non-computational state $|02\rangle$. Why only $|11\rangle$, and how long must the flux be held? Take two three-level Duffing oscillators with frequencies $\omega_j$, anharmonicities $\alpha_j$ and an exchange coupling $g$, all symbols, and let the algebra reduce the nine-level problem to the $2\times 2$ block that is the gate. A numerical tool diagonalizes the nine-level matrix at chosen parameter values; here the block comes out with its parameters intact.

Build the pair Hamiltonian $\sum_j \left[\omega_j \hat n_j + \frac{\alpha_j}{2}\hat n_j(\hat n_j - 1)\right] + g(\hat b_1^\dagger \hat b_2 + \hat b_2^\dagger \hat b_1)$ on three levels per transmon:

```wl
b = DiagonalMatrix[Sqrt[{1, 2}], 1];
{b1, b2} = {KroneckerProduct[b, IdentityMatrix[3]], KroneckerProduct[IdentityMatrix[3], b]};
duffing[bj_, ω_, α_] := With[{nj = Transpose[bj] . bj}, ω nj + α/2 nj . (nj - IdentityMatrix[9])];
hPair = duffing[b1, ω1, α1] + duffing[b2, ω2, α2] + g (Transpose[b1] . b2 + Transpose[b2] . b1);
With[{nTot = Transpose[b1] . b1 + Transpose[b2] . b2}, hPair . nTot == nTot . hPair]
```
<!-- run time: 0.001 s -->

The total excitation number commutes with the Hamiltonian, so the two-excitation states $|11\rangle, |02\rangle, |20\rangle$ evolve only among themselves. The flux pulse tunes $|11\rangle$ to $|02\rangle$ and leaves $|20\rangle$ detuned, so the gate lives in a two-state block.

Extract the block on $\{|11\rangle, |02\rangle\}$, with energies measured from $E_{02}$:

```wl
pos = 1 + FromDigits[#, 3] & /@ {{1, 1}, {0, 2}};
blk = Simplify[hPair[[pos, pos]] - hPair[[pos[[2]], pos[[2]]]] IdentityMatrix[2]]
```
<!-- run time: 0.001 s -->

The diagonal is the detuning $\delta = \omega_1 - \omega_2 - \alpha_2$, and the coupling is $\sqrt 2\, g$, not $g$: lowering transmon 1 and raising transmon 2 from one to two quanta picks up $\sqrt 2$ from $\hat b^\dagger|1\rangle = \sqrt 2|2\rangle$. In other words, $|11\rangle$ has a resonant partner in the second excited level, while the one-excitation states $|01\rangle, |10\rangle$ are coupled by $g$ but detuned by $\alpha_2$, so they only shift by order $g^2/\alpha_2$. The two-level block neglects that shift; that asymmetry is what makes the phase conditional.

Flux-tune to the crossing ($\omega_1 = \omega_2 + \alpha_2$) and evolve for a time $t$:

```wl
cz = FullSimplify[MatrixExp[-I (blk /. ω1 -> ω2 + α2) t]]
```
<!-- run time: 0.02 s -->

At resonance the block is $\sqrt 2\,g\,\sigma_x$, so the population swaps $|11\rangle \to |02\rangle \to |11\rangle$ at the rate $\sqrt 2\,g$. After one full swap, $t = \pi/(\sqrt 2\, g)$, the amplitude to remain in $|11\rangle$ is $\cos\pi = -1$ for every $g$: the $2\pi$ rotation sign of a spinor, which within this block only $|11\rangle$ acquires.

Check the $|11\rangle$ sign against QuantumFramework's CZ gate, with the three other diagonal entries set to 1 as the two-level block assumes:

```wl
DiagonalMatrix[{1, 1, 1, cz[[1, 1]] /. t -> π/(Sqrt[2] g)}] == Normal[QuantumOperator["CZ"]["Matrix"]]
```
<!-- run time: 0.002 s -->

## The Quantum Zeno Effect: an Exceptional Point of the Liouvillian at $k = \Omega$

A qubit driven at Rabi frequency $\Omega$ is measured continuously along $\sigma_z$ at strength $k$, which on average adds the dephasing $k(\sigma_z\rho\sigma_z - \rho)$. Weak watching only damps the Rabi oscillation; strong watching freezes the qubit (the quantum Zeno effect). Where does one regime turn into the other? A numerical tool computes the spectrum at sampled $k$ and sees the crossover as a kink in a plot. Here $k$ and $\Omega$ stay symbols.

Build the Liouvillian and ask for its eigenvalues:

```wl
zeno = Normal[QuantumOperator["Liouvillian"[Ω/2 QuantumOperator["X"], {QuantumOperator["Z"]}, {k}]]["Matrix"]];
Eigenvalues[zeno]
```
<!-- run time: 0.28 s -->

The eigenvalue $0$ is the steady state and $-2k$ is the decay of the coherence that the drive does not touch. The other two, $-k \pm \sqrt{k^2 - \Omega^2}$, carry the physics: for $k < \Omega$ they are complex, a Rabi oscillation that the measurement damps; for $k > \Omega$ they are real, two overdamped rates. They meet at $k = \Omega$.

Count independent eigenvectors at the meeting point and, for comparison, at $k = 2\Omega$:

```wl
MatrixRank[Eigenvectors[zeno /. k -> #]] & /@ {Ω, 2 Ω}
```
<!-- run time: 0.006 s -->

At $k = \Omega$ the two eigenvectors merge as well as the eigenvalues, so the Liouvillian cannot be diagonalized there. In other words, the crossover is an exceptional point, not a mere degeneracy.

Expand the slower of the two real rates for strong measurement:

```wl
Series[-k + Sqrt[k^2 - Ω^2], {k, Infinity, 1}]
```
<!-- run time: 0.007 s -->

The slowest relaxation rate falls as $\Omega^2/(2k)$: past the exceptional point, measuring harder makes the qubit move more slowly, which is the Zeno effect.
