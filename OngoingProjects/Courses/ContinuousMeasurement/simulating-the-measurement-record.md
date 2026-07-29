---
Template: Default
---

# Trajectories, the Record, and the Master Equation

**Simulating the stochastic Schrödinger equation of a free particle under continuous weak position measurement in the Wolfram Language. Continuous homodyne detection produces a stochastic measurement record; conditioning the state on that record unravels the Lindblad master equation into individual quantum trajectories, each a pure conditional state undergoing measurement-induced localization; and averaging the conditional states over the measurement noise recovers the unconditional, decohering master-equation state. Every simulated quantity is matched to an observable of a continuous-measurement experiment.**

Mads Bahrami (last updated: July 28, 2026)

### Setting the Stage

This is a companion to *Watching a Quantum Particle*, which derives the stochastic Schrödinger equation (SSE) of continuous position measurement. Here I take the SSE as given and answer one question: how does the simulation connect to a real experiment? The answer has three parts, and the document is built around them. Continuous homodyne detection produces a stochastic *measurement record*. Conditioning on one record yields one *quantum trajectory*, a pure conditional state whose mean diffuses and whose variance localizes. Averaging the conditional states over the measurement noise, or equivalently discarding the record, returns the *unconditional state*, which obeys the Lindblad master equation and decoheres. One trajectory is the a-posteriori state inferred in a single run; the unconditional state is the a-priori ensemble that remains when the record is not read.

Every number below is produced by a cell you can rerun and change. This is a live notebook, so evaluate the cells top to bottom and vary the measurement strength, the initial state, and the ensemble size. I assume basic quantum mechanics and working Wolfram Language; no stochastic calculus is needed beyond the defining property of the Wiener process, that its increment over a step of length $dt$ is a zero-mean Gaussian of variance $dt$. Throughout, $\hbar$, the mass $m$, and the measurement strength $\lambda$ (equivalently the measurement rate) stay symbolic wherever the algebra is symbolic, and take numerical values (`hbar`, `mass`) only in simulation.

## Part I: The Equation and the Detector's Output

A particle with Hamiltonian $\hat H$ under continuous position measurement obeys

$$
d|\psi\rangle = \left[-\tfrac{i}{\hbar}\hat H\,dt - \tfrac{\lambda}{2}\bigl(\hat x - \langle\hat x\rangle\bigr)^2 dt + \sqrt{\lambda}\,\bigl(\hat x - \langle\hat x\rangle\bigr)\,dW\right]|\psi\rangle .
$$

The three terms are the free unitary evolution under $\hat H$; the stochastic backaction term proportional to the Wiener increment $dW$, which conditions the state toward a more sharply defined position; and the nonlinear drift $-\tfrac{\lambda}{2}(\hat x-\langle\hat x\rangle)^2$, which preserves the norm of the conditional state. In this diffusive, norm-preserving form the SSE is the quantum-state-diffusion (QSD) equation of Gisin and Percival. The Wiener increment $dW$ is the detector shot noise and, equivalently, the innovation of the associated quantum filter, and the same increment appears in the measurement record output by the detector,

$$
dY_t = 2\sqrt{\lambda}\,\langle\hat x\rangle_t\,dt + dW_t .
$$

The record is the conditional mean position $\langle\hat x\rangle_t$ superposed on white noise. This is precisely the output of a continuous position detector: the homodyne photocurrent of light scattered by the particle, or the voltage of a probe following its motion. Over one time step $dt$ the detector returns one outcome $\bar x$ of a weak position measurement, the state is acted on by the measurement (Kraus) operator $\hat M(\bar x)\propto e^{-\lambda\,dt\,(\hat x-\bar x)^2}$, whose $\hat M^\dagger\hat M$ is the POVM element for outcome $\bar x$, and renormalized, and it then evolves unitarily for $dt$:

$$
|\psi(t+dt)\rangle = \frac{e^{-i\hat H\,dt/\hbar}\,\hat M(\bar x)\,|\psi(t)\rangle}{\|\hat M(\bar x)\,|\psi(t)\rangle\|} .
$$

The SSE is exactly the $dt\to 0$ limit of this conditional update, with the outcome $\bar x$ distributed according to the measurement probability density $\|\hat M(\bar x)|\psi\rangle\|^2$: it is the equation of motion for $|\psi\rangle$ conditioned on the record, and its nonlinearity, the $\langle\hat x\rangle$ terms, is exactly the cost of the renormalization by $\|\hat M(\bar x)|\psi\rangle\|$. Because $e^{-i\hat H\,dt/\hbar}$ is unitary it preserves the norm, so the split ordering used in the code, measure then evolve, requires only the single normalization shown. The measurement strength $\lambda$ sets the information rate of each weak measurement: larger $\lambda$ gives a sharper, more informative outcome and stronger measurement backaction.

### Glossary

The six symbols of Part I and their operational status. The last column states whether each is measurable, and if so how; the essential point is that only the record is measured, while the conditional mean is inferred from it.

| Symbol | Meaning | Observable, and how it is tested |
|---|---|---|
| $Y_t$ | the integrated measurement record, $dY_t = 2\sqrt{\lambda}\,\langle\hat x\rangle\,dt + dW_t$ | Directly observed: the (integrated) homodyne photocurrent. |
| $\bar x$ | one weak-measurement outcome per step | Directly observed; a single outcome has signal-to-noise below one. |
| $dW_t$ | the Wiener increment, the detector shot noise, equal to the innovation $dY_t - 2\sqrt{\lambda}\,\langle\hat x\rangle\,dt$ | Recoverable once $\langle\hat x\rangle$ is estimated: the record minus its conditional mean. |
| $\langle\hat x\rangle$ | the conditional mean position, $\int x\,\lvert\psi(x)\rvert^2\,dx$ | Not a single-shot readout; the filtered position estimate (the red curve), tested as the estimator output and used in feedback. |
| $\hat x$ | the position operator, the measured observable | Observable: its weak-measurement outcomes are what the detector samples. |
| $\hat M(\bar x)$ | the Kraus (measurement) operator, $\propto e^{-\lambda\,dt\,(\hat x-\bar x)^2}$ | Not observable; it defines the measurement update and the outcome statistics. |

Only the record $Y_t$, through its outcomes $\bar x$, is measured; the conditional mean $\langle\hat x\rangle$ is inferred from it by filtering, which is exactly what Part VIII shows a real experiment doing in real time.

## Part II: How to Simulate It

We build this on a position grid, one piece at a time. First fix units, $\hbar = m = 1$:

```wl
hbar = 1.; mass = 1.;
```

The computational grid holds the position samples and the matching momenta in the wrapped ordering `Fourier` uses (nonnegative frequencies first, the negative half rotated to the end). The two arguments are `n`, the number of grid points, and `ll`, the box length $L$; the returned association carries the spacing $dx = L/n$ under `"dx"`, the positions under `"x"`, and the momenta $p = (2\pi\hbar/L)\times\text{integers}$ under `"p"`:

```wl
grid[n_, ll_] := With[{dx = ll/n},
   <|"n" -> n, "dx" -> dx, "x" -> dx Range[-n/2, n/2 - 1],
     "p" -> (2 Pi hbar/ll) RotateLeft[Range[-n/2, n/2 - 1], n/2]|>];
```

A normalized Gaussian wave packet on grid `g`, of position width $s$ (so $V_x = s^2$), centered at $x_0$ with mean momentum $p_0$:

```wl
gaussian[g_, x0_, p0_, s_] := Normalize[Exp[-(g["x"] - x0)^2/(4 s^2) + I p0 g["x"]/hbar]];
```

Free unitary evolution of the wavefunction `ψ` on grid `g` for a step `dt` is diagonal in momentum, so it is a split-operator (Fourier) step: forward transform, multiply by the kinetic propagator $e^{-i p^2 dt/(2 m\hbar)}$, inverse transform:

```wl
unitaryStep[g_, dt_][ψ_] := InverseFourier[
   Exp[-I g["p"]^2 dt/(2 mass hbar)] Fourier[ψ, FourierParameters -> {0, -1}],
   FourierParameters -> {0, -1}];
```

The conditional mean position $\langle\hat x\rangle$ is the inner product of the grid positions with the probability density $|\psi|^2$ (arguments: grid `g`, wavefunction `ψ`):

```wl
xMean[g_, ψ_] := g["x"] . Abs[ψ]^2;
```

and the conditional position variance $V_x = \langle\hat x^2\rangle - \langle\hat x\rangle^2$ is the second moment minus the square of the first:

```wl
xVar[g_, ψ_] := g["x"]^2 . Abs[ψ]^2 - xMean[g, ψ]^2;
```

A measurement outcome is sampled from the exact outcome distribution $\|\hat M(\bar x)|\psi\rangle\|^2$: a point drawn from the density $|\psi|^2$ and convolved with detector Gaussian noise of standard deviation $1/\sqrt{4\lambda\,dt}$. The arguments are the grid `g`, the measurement strength `λ`, and the step `dt`:

```wl
drawExact[g_, λ_, dt_][ψ_] := RandomChoice[Abs[ψ]^2 -> g["x"]] +
   RandomVariate[NormalDistribution[0, 1/Sqrt[4 λ dt]]];
```

The measurement update applies the Kraus operator $\hat M(\bar x)$ for outcome `xb` and renormalizes the conditional state:

```wl
measUpdate[g_, λ_, dt_, xb_][ψ_] := Normalize[Exp[-λ dt (g["x"] - xb)^2] ψ];
```

and one full time step is the split conditional update from Part I, measurement followed by unitary evolution:

```wl
stepAt[g_, λ_, dt_, xb_][ψ_] := unitaryStep[g, dt][measUpdate[g, λ, dt, xb][ψ]];
```

One experimental run over `nt` steps draws a measurement outcome at each step, records it, and propagates the conditional state from initial state `ψ0`. `Reap` and `Sow` accumulate the measurement record while `NestList` accumulates the quantum trajectory:

```wl
expt[g_, λ_, dt_, nt_][ψ0_] := Reap[
   NestList[stepAt[g, λ, dt, Sow@drawExact[g, λ, dt][#]][#] &, ψ0, nt]];
```

## Part III: One Run Is One Experiment

Simulate one run, starting from a broad packet so the measurement-induced localization is visible. Fix the grid, the measurement strength `λR`, the step `dt`, and the number of steps `nt`, and retain the full trajectory of conditional states and the measurement record:

```wl
g = grid[256, 50.]; λR = 1.; dt = 0.005; nt = 600;
{states, {record}} = BlockRandom[SeedRandom[11]; expt[g, λR, dt, nt][gaussian[g, 0., 0., 2.]]];
```

Extract the conditional mean $\langle\hat x\rangle$ and the conditional width $\sqrt{V_x}$ along the trajectory:

```wl
xt = xMean[g, #] & /@ states; widths = Sqrt[xVar[g, #] & /@ states];
```

Check the step count, the localization of the width, and the purity of the final conditional state:

```wl
AssociationThread[{"steps", "\!\(\*SqrtBox[\"Vx\"]\): start -> end", "purity of final state"},
  {Length[record], {First@widths, Last@widths},
   Re@Tr[#.#] &@ Outer[Times, Last@states, Conjugate@Last@states]}]
```

The conditional width $\sqrt{V_x}$ relaxes from $2$ to the steady value $(\hbar/(4\lambda m))^{1/4}\approx 0.71$, the square root of the steady conditional variance $V_x = \sqrt{\hbar/(4\lambda m)} = \tfrac12$, at which measurement-induced localization balances free dispersion, and the final conditional state is exactly pure. This is the defining property of an efficiently monitored trajectory: in a single run the particle occupies a definite, if randomly placed, pure state. Continuous observation does not mix the state; it localizes it.

Compare the raw measurement record against the conditional mean inferred from it. The record is the gray point cloud; the conditional mean $\langle\hat x\rangle$ estimated from it is the red curve:

```wl
Show[
 ListPlot[Transpose[{dt Range[nt], record}],
  PlotStyle -> Directive[Gray, Opacity[0.3], PointSize[0.004]]],
 ListLinePlot[Transpose[{dt Range[0, nt], xt}], PlotStyle -> Directive[Thick, Red]],
 Frame -> True, ImageSize -> 540, AspectRatio -> 1/2, PlotRange -> All,
 FrameLabel -> {"time t", "detector reading, and the mean ⟨x⟩ you infer"},
 PlotLabel -> "One run: the noisy record (gray) and the trajectory it determines (red)"]
```

The individual outcomes scatter over $\pm 20$ while the particle remains within a few units of the origin: the signal-to-noise ratio of a single weak measurement is far below one. The conditional mean is recovered only after the measurement operator has integrated hundreds of outcomes. The record is the raw data; the trajectory is the estimate inferred from it.

It is worth stating precisely what $\langle\hat x\rangle$ is, since the detector never reports it directly. It is the first moment of the conditional state, $\langle\hat x\rangle = \int x\,|\psi(x)|^2\,dx$, the minimum-variance estimate of the position given the record so far. No single outcome is that estimate; a lone $\bar x$ has signal-to-noise below one, and $\langle\hat x\rangle$ is recovered only by propagating the full record through the conditional dynamics, which is what the trajectory does. It is worth tracking for three reasons: it is the signal carried by the record, $dY_t = 2\sqrt{\lambda}\,\langle\hat x\rangle\,dt + dW_t$; the particle is localized within $\sqrt{V_x}$ of it, so the conditional mean and variance together specify the Gaussian state; and a feedback controller can act only on the estimate it computes in real time from the photocurrent, which is exactly this one. So $\langle\hat x\rangle(t)$ is not the eigenvalue of a single projective measurement but the filtered position estimate reconstructed from the continuous record, as physical as any conditional expectation.

The conditional width relaxes along the same deterministic curve in every run:

```wl
ListLinePlot[Transpose[{dt Range[0, nt], widths}],
 Frame -> True, GridLines -> Automatic, ImageSize -> 480, AspectRatio -> 1/2,
 PlotRange -> {0, All}, FrameLabel -> {"time t", "packet width \!\(\*SqrtBox[\"Vx\"]\)"},
 PlotLabel -> "The measurement sharpens the packet to a fixed width"]
```

The conditional mean and width are moments of the full conditional density $|\psi(x,t)|^2$, which is worth displaying directly and comparing between two runs. Restricting to the central region and reading time upward: both conditional densities localize to the same steady width, but each conditional mean follows an independent stochastic path, so the density ridge differs between runs.

```wl
band10 = Flatten@Position[UnitStep[10 - Abs[g["x"]]], 1];
xr = g["x"][[{First@band10, Last@band10}]];
showDensity[st_, lbl_] := ArrayPlot[Reverse[(Abs[st]^2)[[All, band10]]],
   DataRange -> {xr, {0, nt dt}}, ColorFunction -> "SunsetColors", AspectRatio -> 1.4,
   ImageSize -> 250, Frame -> True, FrameLabel -> {"x", "t"}, PlotLabel -> lbl];
states2 = BlockRandom[SeedRandom[4]; First@expt[g, λR, dt, nt][gaussian[g, 0., 0., 2.]]];
GraphicsRow[{showDensity[states, "run 1"], showDensity[states2, "run 2"]}, ImageSize -> 540]
```

## Part IV: Many Runs Are Many Possible Experiments

A different realization of the measurement noise produces a different trajectory. A helper returns the conditional-mean path $\langle\hat x\rangle(t)$ of one run:

```wl
pathX[g_, λ_, dt_, nt_][ψ0_] := xMean[g, #] & /@
   NestList[stepAt[g, λ, dt, drawExact[g, λ, dt][#]][#] &, ψ0, nt];
```

Simulate an ensemble of forty:

```wl
manyX = BlockRandom[SeedRandom[3]; Table[pathX[g, λR, dt, nt][gaussian[g, 0., 0., 2.]], {40}]];
```

Overlay the conditional means, one faint curve per realization:

```wl
ListLinePlot[Transpose[{dt Range[0, nt], #}] & /@ manyX,
 PlotStyle -> Directive[Opacity[0.4]], Frame -> True, ImageSize -> 540, AspectRatio -> 1/2,
 FrameLabel -> {"time t", "inferred mean ⟨x⟩"},
 PlotLabel -> "Forty runs: each detector record sends the particle on its own walk"]
```

Those curves are the estimates; the underlying measurement records also differ from run to run. Each realization has its own outcome sequence $\bar x_k$ and its own integrated record $Y_t = 2\sqrt{\lambda}\int\bar x\,dt$. The left panel shows the raw outcomes for three realizations, dominated by shot noise; the right shows the integrated record for thirty, in which the conditional-mean signal appears only as a slow drift beneath the Wiener (Brownian) spread.

```wl
manyRec = BlockRandom[SeedRandom[3];
   Table[expt[g, λR, dt, nt][gaussian[g, 0., 0., 2.]][[2, 1]], {30}]];
manyY = (2 Sqrt[λR] dt Accumulate[#]) & /@ manyRec;
GraphicsRow[{
  ListPlot[Transpose[{dt Range[nt], #}] & /@ manyRec[[1 ;; 3]],
   PlotStyle -> Directive[Opacity[0.35], PointSize[0.004]], Frame -> True,
   FrameLabel -> {"time t", "reading"}, PlotLabel -> "three raw records", ImageSize -> 280],
  ListLinePlot[Transpose[{dt Range[nt], #}] & /@ manyY,
   PlotStyle -> Directive[Opacity[0.4]], Frame -> True,
   FrameLabel -> {"time t", "integrated record Y"}, PlotLabel -> "thirty integrated records", ImageSize -> 280]
  }, ImageSize -> 600]
```

Each conditional mean follows an independent path and the ensemble spreads over several units. But at any fixed time all realizations share the identical conditional width, because the conditional variance obeys a deterministic (noiseless) equation while only the mean is driven by the record. To confirm that the ensemble spread of the means is large while that of the variances vanishes, take forty final conditional states:

```wl
finals = BlockRandom[SeedRandom[3];
   Table[Last@NestList[stepAt[g, λR, dt, drawExact[g, λR, dt][#]][#] &, gaussian[g, 0., 0., 2.], nt], {40}]];
```

and compare the ensemble standard deviation of the conditional means against that of the conditional variances:

```wl
AssociationThread[{"spread of centers ⟨x⟩ across runs", "spread of widths Vx across runs"},
  {StandardDeviation[xMean[g, #] & /@ finals], StandardDeviation[xVar[g, #] & /@ finals]}]
```

The conditional means scatter by several units; the conditional variances agree to fourteen digits. Every realization yields a differently centered conditional state of the same width. This is the outcome distribution obtained by repeating the experiment.

## Part V: Averaging Gives the Master Equation

If the measurement record is discarded, or equivalently the conditional states are averaged over the measurement noise without post-selecting on their records, the description of the particle is the unconditional density matrix $\rho = E(|\psi\rangle\langle\psi|)$, where $E(\cdots)$ denotes the ensemble average over the Wiener noise. This unconditional state obeys the Lindblad master equation

$$
\frac{d\rho}{dt} = -\frac{i}{\hbar}[\hat H,\rho] - \frac{\lambda}{2}\bigl[\hat x,[\hat x,\rho]\bigr],
$$

the standard equation of spatial decoherence: the double-commutator dissipator suppresses off-diagonal coherence in the position basis at a rate set by $\lambda$. It integrates by the same operator split used for the trajectories. On the position grid the dissipator damps each element $\rho_{ij}$ by $e^{-\frac{\lambda}{2}(x_i-x_j)^2 dt}$, and the Hamiltonian part is the same free unitary evolution. Fix a grid `gS`, a step `dtS`, a horizon `tf`, and the number of steps `ntS`, starting from an initial pure state `ψ0`:

```wl
gS = grid[32, 16.]; dtS = 0.0005; tf = 0.4; ntS = Round[tf/dtS]; ψ0 = gaussian[gS, 0., 0., 1.];
```

The dissipative step is one matrix of damping factors applied to $\rho$ elementwise:

```wl
decohere = Exp[-(λR/2) Outer[(#1 - #2)^2 &, gS["x"], gS["x"]] dtS];
```

The Hamiltonian step is the same split-operator propagator as before, assembled as a matrix `uni` so it can act on both sides of $\rho$:

```wl
uni = Transpose[unitaryStep[gS, dtS] /@ IdentityMatrix[32]];
```

One master-equation step applies the dissipator, then the unitary:

```wl
lindStep[ρ_] := uni . (ρ decohere) . ConjugateTranspose[uni];
```

Propagate from the initial pure state to the horizon:

```wl
ρMaster = Nest[lindStep, Outer[Times, ψ0, Conjugate@ψ0], ntS];
```

Now average many conditional trajectories on the same grid and compare. Simulate eight hundred realizations and retain the final conditional state of each:

```wl
finalsS = BlockRandom[SeedRandom[7]; Table[
    Nest[stepAt[gS, λR, dtS, drawExact[gS, λR, dtS][#]][#] &, ψ0, ntS], {800}]];
```

The ensemble average is the unconditional density matrix $\rho = E(|\psi\rangle\langle\psi|)$:

```wl
ρAvg = Mean[Outer[Times, #, Conjugate[#]] & /@ finalsS];
```

To evaluate observables, represent $\hat x$ and $\hat p^2$ as matrices on this grid, $\hat x$ diagonal and $\hat p^2$ via the Fourier transform:

```wl
xm = DiagonalMatrix[gS["x"]];
p2m = 2 mass Transpose[InverseFourier[gS["p"]^2/(2 mass) Fourier[#, FourierParameters -> {0, -1}],
      FourierParameters -> {0, -1}] & /@ IdentityMatrix[32]];
```

A helper returns $\langle\hat x\rangle$, $\langle\hat x^2\rangle$, $\langle\hat p^2\rangle$, and the purity $\mathrm{Tr}(\rho^2)$ from any density matrix:

```wl
observables[ρ_] := {Re@Tr[ρ . xm], Re@Tr[ρ . xm . xm], Re@Tr[ρ . p2m], Re@Tr[ρ . ρ]};
```

Place the master-equation state and the trajectory ensemble average side by side:

```wl
TableForm[{observables[ρMaster], observables[ρAvg]},
 TableHeadings -> {{"master equation", "800 trajectories"}, {"⟨x⟩", "⟨x²⟩", "⟨p²⟩", "purity"}}]
```

The two agree to the Monte-Carlo sampling error of eight hundred realizations, $1/\sqrt{800}\approx4\%$. The ensemble average of the conditional states reproduces the unconditional master-equation state: this is the consistency condition of the trajectory unravelling. In the last column, the purity $\mathrm{Tr}(\rho^2)$ of the average is about $0.62$, well below one, so the unconditional state is *mixed*, yet every conditional trajectory is exactly pure. The mixed, decohering state of standard decoherence theory is precisely the ensemble average over these pure conditional trajectories, one of which is realized, and knowable, in each run.

The separation between conditional and unconditional purity develops in time. Each conditional trajectory remains pure, its $\mathrm{Tr}(\rho^2)$ fixed at $1$; the unconditional average decoheres, its purity decreasing from $1$ toward $\approx 0.62$. On one pair of axes, the constant line and the decaying curve are the entire content of decoherence as an ensemble average over pure conditional trajectories:

```wl
trajS = BlockRandom[SeedRandom[7];
   Table[NestList[stepAt[gS, λR, dtS, drawExact[gS, λR, dtS][#]][#] &, ψ0, ntS][[1 ;; ;; 25]], {300}]];
purTimes = dtS Range[0, ntS, 25];
purEns = Table[
   With[{ρ = Mean[Outer[Times, #, Conjugate[#]] & /@ trajS[[All, k]]]}, Re@Tr[ρ . ρ]],
   {k, Length@trajS[[1]]}];
ListLinePlot[{Transpose[{purTimes, purEns}], Transpose[{purTimes, ConstantArray[1., Length@purTimes]}]},
 PlotStyle -> {Directive[Thick, ColorData[97, 1]], Directive[Dashed, Gray]}, PlotRange -> {0, 1.05},
 Frame -> True, FrameLabel -> {"time t", "purity Tr(ρ²)"},
 PlotLegends -> {"average over runs", "each single run"},
 PlotLabel -> "Pure trajectories, mixed average", ImageSize -> 460]
```

Is each simulated state physical? A conditional state is pure, so $|\psi\rangle\langle\psi|$ must be positive semidefinite with eigenvalues $\{1,0,\dots,0\}$, and the unconditional average $\rho$ must be a valid density matrix: positive semidefinite, Hermitian, of unit trace. Take the smallest eigenvalue of each conditional projector, over a sample of realizations:

```wl
projMin = Min[Re@Eigenvalues[Outer[Times, #, Conjugate[#]]]] & /@ finalsS[[;; 20]];
```

and report the per-realization norm error, the smallest projector eigenvalue, and the unconditional $\rho$'s smallest eigenvalue, trace, and deviation from Hermiticity:

```wl
AssociationThread[
  {"per-run ‖ψ‖ error", "per-run |ψ⟩⟨ψ| min eigenvalue", "ρ min eigenvalue", "ρ trace", "ρ Hermiticity error"},
  {Max[Abs[# . Conjugate[#] - 1] & /@ finalsS], Min@projMin,
   Min[Re@Eigenvalues[ρAvg]], Re@Tr[ρAvg], Max@Abs[ρAvg - ConjugateTranspose[ρAvg]]}]
```

Every conditional state stays normalized to machine precision, no projector eigenvalue falls below zero by more than rounding ($\sim10^{-16}$), and the unconditional average is a valid density matrix: positive semidefinite to the same precision, Hermitian, of unit trace. The split-operator scheme keeps each state exactly physical, because it applies the measurement as an exact Kraus operator and renormalizes, never admitting a negative probability.

The third column exhibits measurement backaction heating the particle. The heating rate is exactly computable: the moments close, and the momentum variance obeys an equation whose only measurement input is a constant diffusion source,

```wl
moments = First@DSolve[{x2'[t] == xp[t]/m, xp'[t] == 2 p2[t]/m, p2'[t] == λ ℏ^2,
    x2[0] == v0, xp[0] == 0, p2[0] == q0}, {x2, xp, p2}, t];
{"⟨p²⟩(t)" -> (p2[t] /. moments), "⟨x²⟩(t)" -> Simplify[x2[t] /. moments]}
```

So $\langle\hat p^2\rangle(t) = \langle\hat p^2\rangle_0 + \lambda\hbar^2 t$: continuous position measurement heats the momentum at the constant rate $\lambda\hbar^2$, independent of the mass, the state, and the Hamiltonian. This is the measurement backaction demanded by the uncertainty principle as the cost of position information, and at these parameters it takes $\langle\hat p^2\rangle$ from $0.25$ to $0.25 + \lambda\hbar^2 t_f = 0.65$, matching the master equation above. Increasing $\lambda$ increases the backaction heating without bound; only a sufficiently fine grid can represent the resulting momenta, and a coarse grid fails silently here.

With the moment equations in hand, return to the ensemble of conditional means from Part IV. Their ensemble variance is not free: it is exactly $\langle\hat x^2\rangle(t) - V_x(t)$, the unconditional position variance minus the conditional one. This is the variance decomposition of the quantum filter, the total position variance splitting into the part resolved into the fluctuating conditional mean and the part remaining as the conditional variance, and both terms are the deterministic solutions above. The envelope around the ensemble is therefore not fitted, but read directly from the moment equations:

```wl
ψIV = gaussian[g, 0., 0., 2.];
{v0IV, q0IV} = {g["x"]^2 . Abs[ψIV]^2, g["p"]^2 . Abs[Fourier[ψIV, FourierParameters -> {0, -1}]]^2};
x2master[tt_] := x2[t] /. moments /. {m -> mass, λ -> λR, ℏ -> hbar, v0 -> v0IV, q0 -> q0IV, t -> tt};
band = Table[Sqrt[Max[x2master[k dt] - widths[[k + 1]]^2, 0.]], {k, 0, nt}];
Show[
 ListLinePlot[Transpose[{dt Range[0, nt], #}] & /@ manyX, PlotStyle -> Directive[Opacity[0.25], Thin]],
 ListLinePlot[{Transpose[{dt Range[0, nt], band}], Transpose[{dt Range[0, nt], -band}]},
  PlotStyle -> Directive[Thick, Red]],
 Frame -> True, ImageSize -> 540, AspectRatio -> 1/2,
 FrameLabel -> {"time t", "inferred mean ⟨x⟩, with the predicted band"},
 PlotLabel -> "The spread of estimates is \!\(\*SqrtBox[\"(⟨x²⟩ - Vx)\"]\), from the moment equations"]
```

The band brackets the forty realizations to their Monte-Carlo scatter. The information content of the record, how tightly the conditional mean is determined and how far it fluctuates, is fixed in advance by the same moment equations that set the backaction heating.

## Part VI: How Good Is the Simulation? A Builtin Cross-Check

The split-operator scheme is hand-written, and the way to validate it is to hand the same SSE to Wolfram Language's built-in stochastic integrator and check for agreement. `ItoProcess` represents an Ito stochastic differential equation $dx = a\,dt + b\,dW$ directly, and `RandomFunction` integrates it. The built-in cannot use the spectral propagator, so the Hamiltonian is supplied as an explicit finite-difference kinetic-energy matrix, a tridiagonal second derivative with periodic boundary conditions (arguments: grid size `nn`, spacing `dx`):

```wl
fdKinetic[nn_, dx_] := -hbar^2/(2 mass dx^2) Normal@SparseArray[
   {{i_, i_} -> -2., {i_, j_} /; Abs[i - j] == 1 -> 1., {1, nn} -> 1., {nn, 1} -> 1.}, {nn, nn}];
```

The process itself stacks the real and imaginary parts of $\psi$ into one real state vector, with the drift $a$ and the scalar-noise diffusion $b$ read directly from the SSE:

```wl
sseProcess[g_, λ_, ψ0_] := Module[
   {xs = g["x"], kmat = fdKinetic[g["n"], g["dx"]], us, vs, xav, a, drift, diff},
   us = Table[Unique["u"], g["n"]]; vs = Table[Unique["v"], g["n"]];
   xav = xs . (us^2 + vs^2); a = xs - xav;
   drift = Join[kmat . vs/hbar - (λ/2) a^2 us, -kmat . us/hbar - (λ/2) a^2 vs];
   diff = Transpose[{Sqrt[λ] Join[a us, a vs]}];
   ItoProcess[{drift, diff}, {Join[us, vs], Join[Re@ψ0, Im@ψ0]}, {t, 0}]];
```

Choose a small grid the built-in can afford, and an initial packet:

```wl
gI = grid[24, 12.]; ψI = gaussian[gI, 0., 0., 1.]; nI = gI["n"];
```

Integrate four hundred trajectories with `RandomFunction`, using the scalar-noise stochastic Runge-Kutta method:

```wl
tdIto = BlockRandom[SeedRandom[7]; RandomFunction[sseProcess[gI, λR, ψI],
    {0., tf, dtS}, 400, Method -> "StochasticRungeKuttaScalarNoise"]];
```

Reassemble each final real state vector into a complex wavefunction:

```wl
itoFinals = (#[[;; nI]] + I #[[nI + 1 ;;]]) & /@ tdIto["ValueList"][[All, -1]];
```

For comparison, represent $\hat x$ and $\hat p^2$ as matrices on this grid:

```wl
xmI = DiagonalMatrix[gI["x"]]; p2mI = 2 mass fdKinetic[nI, gI["dx"]];
```

Average the built-in trajectories into a density matrix, renormalizing each first since its norm has drifted:

```wl
itoRho = Mean[Outer[Times, #, Conjugate[#]]/(# . Conjugate[#]) & /@ itoFinals];
```

Compare the built-in's norm drift and observables against the same analytic master-equation moments:

```wl
AssociationThread[
  {"ItoProcess norm drift", "ItoProcess {⟨x²⟩, ⟨p²⟩, purity}", "master equation {⟨x²⟩, ⟨p²⟩}"},
  {Max[Abs[# . Conjugate[#] - 1] & /@ itoFinals],
   {Re@Tr[itoRho . xmI . xmI], Re@Tr[itoRho . p2mI], Re@Tr[itoRho . itoRho]},
   {x2[t], p2[t]} /. moments /. {v0 -> 1., q0 -> Re[Conjugate[ψI] . p2mI . ψI],
      m -> 1., λ -> 1., ℏ -> 1., t -> tf}}]
```

The built-in reproduces the physics: its ensemble matches the master-equation moments within Monte-Carlo scatter, and it decoheres to the same mixed purity. But it is the coarser tool here, for two reasons, both traced to the measured operator. A general-purpose integrator does not preserve the nonlinear norm invariant of the SSE, so $\|\psi\|^2$ departs from $1$ by about $10^{-4}$, against $10^{-15}$ for the split-operator scheme, and each trajectory must be renormalized before observables are evaluated. And because $\hat x$ is unbounded, the norm-preserving $(\hat x-\langle\hat x\rangle)^2$ drift is stiff: it limits the step size, a limit that tightens as the grid grows, so the built-in can afford only a small grid with a finite-difference kinetic energy that slightly underestimates the backaction heating.

Neither limitation is intrinsic to `ItoProcess`. Both follow from propagating the full state vector while the operator is unbounded; in the right coordinates the built-in is exact. A qubit trajectory simulator propagates not the state vector but the three Bloch components $\langle\hat\sigma_x\rangle,\langle\hat\sigma_y\rangle,\langle\hat\sigma_z\rangle$, a complete parametrization of a two-level state, confined to the ball $r\le 1$ that the stochastic master equation holds invariant, so positivity is enforced geometrically and no renormalization is required. The free particle admits the same reduction. Under continuous position measurement a Gaussian state remains Gaussian, so five real parameters suffice: the means $\langle\hat x\rangle,\langle\hat p\rangle$ and the covariances $V_x,C,V_p$, where only the means carry the Wiener noise and the covariances obey the deterministic Riccati flow of Part IV. This is the continuous quantum Kalman-Bucy filter. Supply its five parameters to `ItoProcess` in place of the grid:

```wl
momentProc = ItoProcess[
   {{mp[t]/mass, 0, 2 vc[t]/mass - 4 λR vx[t]^2,
     vp[t]/mass - 4 λR vx[t] vc[t], λR hbar^2 - 4 λR vc[t]^2},
    {{2 Sqrt[λR] vx[t]}, {2 Sqrt[λR] vc[t]}, {0}, {0}, {0}},
    {mx[t], mp[t], vx[t], vc[t], vp[t]}},
   {{mx, mp, vx, vc, vp}, {0., 0., 1., 0., 0.25}}, {t, 0}];
```

Normalization is now automatic, since there is no unnormalized state vector: the conditional state is Gaussian by construction. The purity invariant $V_xV_p-C^2$ replaces the norm, and its deterministic flow has an attracting fixed point at $\hbar^2/4$, so numerical error is contracted rather than accumulated. Integrate at forty times the grid step, over four thousand realizations, and read the backaction heating and the purity error:

```wl
tdM = BlockRandom[SeedRandom[7]; RandomFunction[momentProc, {0., tf, 0.02}, 4000,
    Method -> "StochasticRungeKuttaScalarNoise"]];
finM = tdM["ValueList"][[All, -1]];
AssociationThread[{"purity error", "E[⟨p²⟩](tf)"},
  {Max[Abs[#[[3]] #[[5]] - #[[4]]^2 - 0.25] & /@ finM], Mean[#[[5]] + #[[2]]^2 & /@ finM]}]
```

The purity invariant holds to about $3\times10^{-4}$ at that coarse step, and $\langle\hat p^2\rangle$ reaches the same $0.25+\lambda\hbar^2 t_f=0.65$ as the master equation, with no renormalization. The built-in is now as accurate here as for a qubit, and for the same reason: a complete, low-dimensional parametrization whose single constraint the flow preserves autonomously.

This reduction is the special case, not the general tool. The five-parameter description exists only because the free particle is a linear (Gaussian) system, so the moment hierarchy closes; introduce a potential in $\hat H$ and it does not close, the conditional state is no longer Gaussian, and only the grid split-operator scheme applies. This sets the division of labor. When the conditional state lies in a small closed parametrization, a qubit's Bloch vector or a Gaussian's moments, the built-in integrator suffices and the geometry preserves physicality. Otherwise, hand-write the step whose structure can be exploited, here the exact Kraus operator with its renormalization, which is independent of $\hat H$.

## Part VII: The Simulation-to-Experiment Dictionary

Every simulated quantity has an experimental counterpart. The correspondence is collected here.

| In the simulation | Experimental counterpart |
|---|---|
| the record, the outcomes $\bar x_k$ (or $dY_t$) | the homodyne photocurrent, a continuous noisy position readout |
| one trajectory $\lvert\psi(t)\rangle$ | the conditional state in one run, reconstructed by filtering the photocurrent |
| the conditional mean $\langle\hat x\rangle(t)$, red curve | the filtered position estimate |
| the conditional variance $V_x(t)$ | the residual (a-posteriori) position uncertainty |
| trajectory purity $=1$ | at unit efficiency, monitoring keeps the conditional state pure |
| the unconditional $\rho(t)$, purity $<1$ | the state when the record is discarded: the mixed, decohered master-equation state |
| ensemble spread of trajectories | the shot-to-shot scatter of repeated experiments |
| $\langle\hat p^2\rangle$ growing as $\lambda\hbar^2 t$ | measurement backaction, the momentum diffusion driving the standard quantum limit |

The organizing principle: the master equation is the a-priori state, averaged over every outcome the detector could return; a single trajectory is the a-posteriori state conditioned on one realized record. Filtering the record to reconstruct the conditional state is exactly what a real-time monitoring or feedback experiment does; discarding it leaves the unconditional, decohering state. The simulation provides both, and exposes the pure conditional trajectories comprising the mixed unconditional average.

## Part VIII: Comparison with Experiment

The free particle under continuous position measurement is an idealization, but its measurement-theoretic content is realized directly in cavity and levitated optomechanics, where the position of a mechanical oscillator is monitored continuously by scattered or cavity light. The correspondence is exact for the measurement, and the numbers below are taken from two experiments that reconstruct a mechanical resonator's quantum trajectory from its measurement record.

In Magrini et al. (Nature **595**, 373 (2021), [arXiv:2012.15188](https://arxiv.org/abs/2012.15188)), an optically levitated nanoparticle is monitored by position-resolving light near the Heisenberg imprecision-backaction limit, and its motional quantum state is reconstructed in real time by a Kalman filter, the classical name for the continuous quantum filter of Part VI. The filter output is the experimental analogue of the red curve: the conditional mean and variance propagated from the photocurrent. The reported conditional position uncertainty is $1.3$ times the zero-point value, and feedback on the estimate cools the oscillator to a mean occupation $n = 0.56\pm 0.02$, its motional ground state, from room temperature. Rossi et al. (Nature **563**, 53 (2018), [arXiv:1805.05087](https://arxiv.org/abs/1805.05087); and the companion Phys. Rev. Lett. **123**, 163601 (2019), "Observing and verifying the quantum trajectory of a mechanical resonator") perform the same continuous position measurement on a soft-clamped membrane at near-unit measurement efficiency, reconstruct the conditional state by causal Wiener/Kalman filtering, reach a conditional occupation $n = 0.29\pm 0.03$, and operate $9$ dB below the quantum backaction limit of sideband cooling.

The quantities carrying measured values map as follows:

| Simulated quantity | Measured value |
|---|---|
| the record $dY_t$ | the homodyne photocurrent of the scattered/cavity light |
| conditional mean and variance $\langle\hat x\rangle, V_x$ | the real-time Kalman-filter estimate and its covariance |
| steady conditional position uncertainty | $1.3\times$ zero-point (Magrini 2021) |
| conditional occupation reached by feedback | $n = 0.56\pm 0.02$ (Magrini), $n = 0.29\pm 0.03$ (Rossi) |
| measurement efficiency $\eta$ | near unity (Rossi); sub-unity $\eta$ gives conditional purity below one |
| backaction heating $\lambda\hbar^2$ | radiation-pressure backaction, $9$ dB below the sideband-cooling backaction limit (Rossi) |

The Part VI moment filter is not merely analogous to these experiments; it is the identical algorithm. Both propagate a Gaussian conditional state by the continuous quantum Kalman-Bucy equations, the means driven by the record and the covariances following a deterministic Riccati flow, and the experiments run exactly this on a field-programmable gate array in real time.

Two qualifications are needed. First, the particle here is free, $\hat H = \hat p^2/2m$, whereas these oscillators are harmonically trapped: the measurement machinery, the SSE, the record, the Kalman filter, the backaction, and the conditional-unconditional split, transfers unchanged, but the trap adds a restoring term to the drift, so the steady conditional variance is set by the ratio of the measurement rate to the mechanical and thermal rates rather than by $(\hbar/4\lambda m)^{1/4}$. Second, the experiments run at efficiency $\eta<1$: some scattered light is lost, so the reconstructed conditional state is mixed and its purity is below one, whereas the simulation at $\eta=1$ produces pure conditional trajectories. Restoring $\eta<1$ in the simulation, which replaces the pure-state SSE by a stochastic master equation for a conditional density matrix, reproduces the sub-unity conditional purity these experiments report, and is the extension flagged below.

### Where This Leaves Us

The full loop is now in runnable form: a detector model emitting a stochastic record, one realization mapping a record to a pure conditional trajectory whose mean diffuses while its variance localizes, an ensemble whose average reproduces the Lindblad master equation as a mixture of those pure conditional states, and the exact backaction heating rate $\lambda\hbar^2$ that measurement pays for position information. The most consequential parameter is the measurement strength $\lambda$: reduce it and the conditional state barely responds while the record remains shot-noise dominated; increase it and the conditional state localizes to the steady width almost immediately while the momentum heats faster than a coarse grid can represent. The natural extensions are to reduce the detection efficiency below unity, so that a single record no longer purifies the conditional state and the trajectory becomes a mixed conditional density matrix, or to add a potential to $\hat H$ and study the same record-conditioned dynamics in a trap, the setting of cavity optomechanics.
