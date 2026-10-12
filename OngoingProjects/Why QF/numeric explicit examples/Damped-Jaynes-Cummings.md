---
Template: Default
Title: "When does the photon come back? Jaynes-Cummings with dissipation"
---

# When does the photon come back? Jaynes-Cummings with dissipation

A single two-level atom sits inside a cavity that holds one electromagnetic mode. The atom can give its excitation to the field as a photon and take it back again, and in a closed cavity this exchange is a clean oscillation: the Rabi oscillation. Real cavities are not closed. The atom can decay into modes outside the cavity, and the cavity mirror leaks photons to the outside world. Both losses drain the excitation that the oscillation is passing back and forth.

So the question is: with both channels open, what happens to the exchange? Two things, and this example computes both. First, an atom that starts high up the ladder, with several photons already in the cavity, shows a Rabi oscillation under a decaying envelope: the exchange keeps going while the excitation slowly bleeds away. Second, an atom that starts with a single quantum, so that only one photon can ever be in play, splits into three regimes set by how the coupling compares with the loss rates: the excitation either oscillates under decay, or stops oscillating at a critical ratio, or never oscillates at all and the atom simply decays faster than it would in free space.

## The model: two loss channels draining one exchange

Let $\lvert e\rangle$ and $\lvert g\rangle$ be the excited and ground atomic states, and $\lvert n\rangle$ the cavity Fock states with $n$ photons. The atom is the first subsystem and the field the second. The atomic operators are $\hat\sigma_-=\lvert g\rangle\langle e\rvert$ and $\hat\sigma_+=\hat\sigma_-^\dagger$, and the field operator obeys $\hat a\lvert n\rangle=\sqrt n\,\lvert n-1\rangle$. In the rotating-wave approximation the Jaynes-Cummings Hamiltonian is

$$
\hat H=\hbar\nu\,\hat a^\dagger\hat a
+\tfrac12\hbar\omega\,\hat\sigma_z
+\hbar g\left(\hat\sigma_-\hat a^\dagger+\hat\sigma_+\hat a\right),
\qquad
\hat\sigma_z=\lvert e\rangle\langle e\rvert-\lvert g\rangle\langle g\rvert.
$$

The interaction turns $\lvert e,n\rangle$ into $\lvert g,n+1\rangle$ and back, so it conserves the total number of excitations $\hat a^\dagger\hat a+\hat\sigma_+\hat\sigma_-$. In other words, every photon the field gains is a quantum the atom has lost, and that sum stays fixed under the coupling. We take the resonant case $\nu=\omega$. On resonance the two states $\lvert e,n\rangle$ and $\lvert g,n+1\rangle$ that the coupling mixes carry the same free energy, so the first two terms contribute only an overall phase to each excitation manifold and the excited-state population is set entirely by the coupling and the losses.

Dissipation enters through a Lindblad master equation with two jump operators, atomic decay and cavity photon loss,

$$
\frac{d\hat\rho}{dt}=-\frac{i}{\hbar}[\hat H,\hat\rho]
+\mathcal D[\hat L_1]\hat\rho+\mathcal D[\hat L_2]\hat\rho,
\qquad
\hat L_1=\sqrt{\gamma}\,\hat\sigma_-\otimes\hat I,\quad
\hat L_2=\sqrt{\kappa}\,\hat I\otimes\hat a,
$$

with $\mathcal D[\hat L]\hat\rho=\hat L\hat\rho\hat L^\dagger-\tfrac12(\hat L^\dagger\hat L\hat\rho+\hat\rho\hat L^\dagger\hat L)$. Here $\gamma$ is the atomic decay rate and $\kappa$ the cavity photon-loss rate. The atomic jump $\hat\sigma_-$ lowers the atom without touching the field; the cavity jump $\hat a$ removes a photon without touching the atom.

Load the checked-out QF paclet and its two contexts:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/"]
```

Load the main context:

```wl
Needs["Wolfram`QuantumFramework`"]
```

Expose the bosonic (second-quantization) context, which carries `AnnihilationOperator` and `FockState`:

```wl
Needs["Wolfram`QuantumFramework`SecondQuantization`"]
```

The bosonic constructors truncate the cavity to `$FockSize` levels. Twelve is comfortable for an initial state with a couple of photons:

```wl
SetFockSpaceSize[12]
```

Build photon annihilation on the second subsystem; this same operator is also the cavity-loss jump:

```wl
a = AnnihilationOperator[{2}]
```

Use the labels $\lvert0\rangle_a=\lvert e\rangle$ and $\lvert1\rangle_a=\lvert g\rangle$, the excited-first convention. In that convention QF's named `"-"` is exactly $\hat\sigma_-=\lvert g\rangle\langle e\rvert$ on the first subsystem, and `"+"` is its adjoint $\hat\sigma_+$:

```wl
σminus = QuantumOperator["-"]
```

The atomic raising operator:

```wl
σplus = QuantumOperator["+"]
```

For QF operators, `x @ y` composes them, aligning the subsystems each addresses and carrying the untouched subsystem through. Form the dimensionless interaction $\hat\sigma_-\hat a^\dagger+\hat\sigma_+\hat a$; `SuperDagger` takes the operator adjoint:

```wl
h = σminus @ SuperDagger[a] + σplus @ a
```

Confirm that the interaction is Hermitian:

```wl
h["HermitianQ"]
```

## Damped Rabi oscillations: exchange under a decaying envelope

Start the atom excited with two photons already in the cavity, $\lvert e,2\rangle$. The atom is `QuantumState["0"]` in the excited-first convention, and `FockState[2]` is the two-photon field:

```wl
ψ0 = QuantumTensorProduct[QuantumState["0"], FockState[2]]
```

Fix the resonant frequencies and a coupling much larger than the losses we will apply, so that many oscillations fit inside the decay time:

```wl
{ν, ω, gc} = {10., 10., 15.}
```

Assemble the full Jaynes-Cummings Hamiltonian:

```wl
H = ν SuperDagger[a] @ a + 1/2 ω QuantumOperator["Z"] + gc h
```

The excited-state probability is the $\lvert e\rangle\langle e\rvert$ entry of the atomic reduced density matrix. In QF the listed subsystem is the one discarded, so trace out the field, subsystem 2, and read the first diagonal entry:

```wl
excitedProbability[evolved_] := Re @ QuantumPartialTrace[evolved, {2}]["DensityMatrix"][[1, 1]]
```

Choose three dissipation strengths with the atom and cavity losing at equal rates $\gamma=\kappa$. QF's `"Hamiltonian"` form takes the generator, the jump operators and their rates, and builds the full Lindblad generator; no explicit $i$ is needed, and each jump operator auto-pads with the identity on the subsystem it leaves alone. Evolve the initial state over a fixed window for each rate:

```wl
trajectories = Table[
    QuantumEvolve[QuantumOperator["Hamiltonian"[H, {σminus, a}, {rate, rate}]], ψ0, {t, 0, 2}],
    {rate, {0.8, 0.4, 0.2}}];
```

Each trajectory is a state parametrized by the time $t$. Plot the excited-state probability of all three:

```wl
ListLinePlot[
    Table[Transpose[{#, (excitedProbability[traj] /. t -> #)}] &[N[Range[0, 200]/100]], {traj, trajectories}],
    PlotRange -> {0, 1}, Frame -> True, FrameLabel -> {"t", "P(e)"},
    PlotLegends -> {"\[Gamma] = \[Kappa] = 0.8", "\[Gamma] = \[Kappa] = 0.4", "\[Gamma] = \[Kappa] = 0.2"},
    PlotLabel -> "Damped Rabi oscillations", ImageSize -> 560, AspectRatio -> 1/2.4
]
```
As one can see, the excited-state probability oscillates under an envelope that decays toward zero, and as the rate grows the envelope collapses sooner. The oscillation runs at $2g\sqrt{n+1}=2g\sqrt3$, the Rabi frequency for the $\lvert e,2\rangle\leftrightarrow\lvert g,3\rangle$ pair: the oscillation is the coherent exchange, the envelope is the loss. In every case the excitation ends in the joint ground state $\lvert g,0\rangle$, which neither channel can leave. Loss does not change the exchange frequency here, only how long the exchange survives.

## One excitation: strong coupling, crossover and the Purcell regime

Give the system a single quantum of energy, the atom excited and the cavity empty, $\lvert e,0\rangle$:

```wl
ψe0 = QuantumTensorProduct[QuantumState["0"], FockState[0]]
```

With one excitation the problem closes on three states. The coupling connects only $\lvert e,0\rangle$ and $\lvert g,1\rangle$. Either loss sends the system to $\lvert g,0\rangle$, which is dark: it has no excitation to lose and the coupling cannot lift it back up. Because the excited atom lives only in $\lvert e,0\rangle$, and the recycling term of the master equation feeds population into $\lvert g,0\rangle$ but never back into $\lvert e,0\rangle$, the excited amplitude follows the non-Hermitian "no-jump" evolution alone. Writing the amplitudes of $\lvert e,0\rangle$ and $\lvert g,1\rangle$ as $c_e$ and $c_1$,

$$
\frac{dc_e}{dt}=-\frac{\gamma}{2}c_e-ig\,c_1,
\qquad
\frac{dc_1}{dt}=-\frac{\kappa}{2}c_1-ig\,c_e,
\qquad
P_e=\lvert c_e\rvert^2.
$$

Encode the $2\times2$ generator of these equations:

```wl
noJump = {{-γ/2, -I g}, {-I g, -κ/2}};
```

Its eigenvalues set the dynamics:

```wl
Simplify[Eigenvalues[noJump]]
```

The two eigenvalues are $\tfrac14\!\left(-\gamma-\kappa\pm\sqrt{(\gamma-\kappa)^2-16g^2}\right)$. The sign of the discriminant $(\gamma-\kappa)^2-16g^2$ splits the behavior into three regimes:

- **Strong coupling**, $g>\lvert\kappa-\gamma\rvert/4$: the square root is imaginary, the eigenvalues are a complex-conjugate pair, and $P_e$ is a vacuum Rabi oscillation under a decaying envelope: the photon is emitted and reabsorbed several times before it escapes.
- **Crossover**, $g=\lvert\kappa-\gamma\rvert/4$: the two eigenvalues merge, the oscillation disappears, and the decay is critically damped.
- **Purcell regime**, $\kappa\gg g$: the cavity empties before the photon can be reabsorbed, so the atom simply decays, but at the enhanced rate $\gamma+4g^2/\kappa$ rather than its bare rate $\gamma$.

We check all three against the full master-equation evolution below, where the $2\times2$ picture and QF's open-system solver must agree.

Set the coupling to one for this comparison, so the loss rates can be read directly against it, and rebuild the Hamiltonian:

```wl
gc = 1;
```

Rebuild with the new coupling:

```wl
H = ν SuperDagger[a] @ a + 1/2 ω QuantumOperator["Z"] + gc h;
```

The excited amplitude is the matrix exponential of the generator applied to the excited state $(1,0)$, and the probability is the squared modulus of its first entry. Build the closed form once per rate pair, with the generator kept exact, so that at the crossover, where the generator is defective, the matrix exponential stays symbolic rather than taking an ill-conditioned numerical route:

```wl
peExact[κv_, γv_] := peExact[κv, γv] = With[
    {m = noJump /. {g -> gc, κ -> Rationalize[κv], γ -> Rationalize[γv]}},
    Function[τ, Evaluate[Abs[First[MatrixExp[m τ] . {1, 0}]]^2]]];
```

Pick one loss pair in each regime, holding the atomic rate at $\gamma=0.1$ and moving the cavity rate $\kappa$ through $\{0.4,4.1,20\}$. At $g=1$ these sit on the strong, crossover and Purcell sides of $g=\lvert\kappa-\gamma\rvert/4$:

```wl
regimes = {{0.4, 0.1}, {4.1, 0.1}, {20., 0.1}};
```

Evolve the one-excitation state under the full Lindblad dynamics for each regime; the rate list pairs $\gamma$ with $\hat\sigma_-$ and $\kappa$ with $\hat a$:

```wl
evolutions = Table[
    QuantumEvolve[QuantumOperator["Hamiltonian"[H, {σminus, a}, {γκ[[2]], γκ[[1]]}]], ψe0, {t, 0, 10}],
    {γκ, regimes}];
```

Confirm that the full open-system evolution reproduces the exact two-state result. Sample both on a grid and report the largest discrepancy across all three regimes:

```wl
Max @ Table[
    Abs[(excitedProbability[evolutions[[k]]] /. t -> #) - peExact[regimes[[k, 1]], regimes[[k, 2]]][#]] & /@ N[Range[0, 100]/10],
    {k, 3}]
```

As expected, the two agree to the accuracy of the time integration: the $2\times2$ no-jump picture captures the full master-equation dynamics of the excited population, as the dark-state argument promised. Plot the three regimes together, the exact result as solid lines and the QF evolution as dots:

```wl
Show[
    Plot[Evaluate[peExact[#[[1]], #[[2]]][t] & /@ regimes], {t, 0, 10},
        PlotStyle -> {Blue, Darker[Green], Red}, PlotRange -> {0, 1}],
    ListPlot[
        Table[Transpose[{#, (excitedProbability[evolutions[[k]]] /. t -> #)}] &[N[Range[0, 40]/4]], {k, 3}],
        PlotStyle -> (Directive[#, PointSize[0.012]] & /@ {Blue, Darker[Green], Red})],
    Plot[Exp[-(0.1 + 4 gc^2/20.) t], {t, 0, 10}, PlotStyle -> Directive[Gray, Dashed]],
    Frame -> True, Axes -> False, FrameLabel -> {"t", "P(e)"},
    PlotLabel -> "One excitation: strong coupling, crossover, Purcell",
    PlotLegends -> SwatchLegend[{Blue, Darker[Green], Red, Gray}, {"strong: \[Kappa] = 0.4", "crossover: \[Kappa] = 4.1", "Purcell: \[Kappa] = 20", "exp(-(\[Gamma] + 4 g^2/\[Kappa]) t)"}],
    ImageSize -> 560
]
```
As one can see, in the strong-coupling curve the excitation sloshes between atom and cavity several times before it leaks away: the photon is reabsorbed before it escapes. At crossover the oscillation is just gone, the curve falling monotonically. In the Purcell curve the atom decays as a pure exponential, and the dashed reference $e^{-(\gamma+4g^2/\kappa)t}$ sits on top of it: a leaky cavity speeds the atom's decay above its bare rate, the Purcell effect, because the cavity carries the excitation away before the atom can take it back.

The lesson is one dimensionless comparison. When the coupling beats the loss, $g>\lvert\kappa-\gamma\rvert/4$, the excitation is shared coherently and the oscillation survives; when the loss beats the coupling, the excitation leaves before it can return and the oscillation is gone. The same master equation, read in QF through the `"Hamiltonian"` jump-operator form, produces both the damped exchange of the multi-photon ladder and the three single-excitation regimes, and in the solvable one-excitation case it agrees with the elementary $2\times2$ calculation.

<!-- wl cells verified against QuantumFramework HEAD (worktree busy-herschel-1f9c79), WL 15, 2026-10-12. -->
