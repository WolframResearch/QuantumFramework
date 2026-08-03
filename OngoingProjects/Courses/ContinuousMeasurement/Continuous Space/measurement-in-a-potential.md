---
Template: Default
---

# Continuous Measurement in a Potential

**A live, self-contained walk through the simulation of a quantum particle whose position is measured continuously and weakly, generalized from the free particle to any one-dimensional Hamiltonian $\hat H=\hat p^2/2m+V(\hat x)$: we write the stochastic Schrödinger equation as a Bayesian update on a position grid, recover the free particle as the $V=0$ special case, watch a harmonic trap squeeze the conditional state below the vacuum while its unconditional average heats without bound, and watch a double well, where no finite set of moments closes, localize the particle in one well and freeze its tunneling by the quantum Zeno effect.**

Mads Bahrami. Last updated July 29, 2026.

### Setting the Stage: How This Notebook Flows

This notebook covers one topic as a continuous story: continuous weak measurement of a particle's position, written first as a stochastic Schrödinger equation and its Bayesian update, then built into a split-operator simulation on a grid, then run for three Hamiltonians in turn, the free particle, the harmonic trap, and the double well, with the conditional state, the unconditional (Lindblad) average, the moment filter, the measurement backaction, and the quantum Zeno effect each read off a cell you can rerun.

I strongly believe in a computation-first narrative for learning: in a sense, if I cannot compute it, I cannot claim to understand it. So almost nothing here is asserted without a cell that computes it on the spot, and almost every claim is checked against an independent calculation, symbolic against numeric or one representation against another.

This is a live Wolfram notebook, not a static paper. Evaluate the cells top to bottom; many depend on symbols defined earlier, so a fresh kernel wants a fresh top-to-bottom pass. Read it as a movie rather than a reference: each cell sets up the next question, and the prose between cells is the connective tissue, not decoration. My suggestion is to focus on each output and its meaning before worrying about every detail of the input code. And remember that the code is yours: change the potential, the measurement strength $\lambda$, the initial state, the grid, rerun, and see what moves. That is the point.

Prerequisites are light: basic quantum mechanics (states, operators, the position and momentum representations) and enough Wolfram Language to read a one-line function. No stochastic calculus is needed beyond one fact, stated when we need it: over a step $dt$ the Wiener increment $dW$ is a zero-mean Gaussian of variance $dt$. The pace steepens gently in Part IV, where a matrix Riccati equation is solved in closed form; the two companion essays *Watching a Quantum Particle* and *Trajectories, the Record, and the Master Equation* develop the free-particle case at more length and are optional background, not requirements. Throughout, $\hbar$, the mass $m$, and $\lambda$ stay symbolic wherever the algebra is symbolic, and take numerical values (`hbar`, `mass`) only in simulation. Comments and corrections are welcome.

Let's start!

## Part I: The Stochastic Schrödinger Equation: Measurement as Bayes' Rule

Start with the equation of motion. A particle with Hamiltonian $\hat H$ whose position is measured continuously and weakly obeys a stochastic Schrödinger equation (SSE), the diffusive quantum-state-diffusion form,

$$
d|\psi\rangle = \left[-\tfrac{i}{\hbar}\hat H\,dt \;-\; \tfrac{\lambda}{2}\bigl(\hat x - \langle\hat x\rangle\bigr)^2 dt \;+\; \sqrt{\lambda}\,\bigl(\hat x - \langle\hat x\rangle\bigr)\,dW\right]|\psi\rangle .
$$

In words: the first term is ordinary unitary evolution under $\hat H$; the last is the stochastic measurement backaction, proportional to the Wiener increment $dW$, which nudges the state toward a sharper position; and the middle term is the nonlinear drift whose only job is to keep $|\psi\rangle$ normalized. The strength $\lambda$ (the measurement rate) sets how fast position information arrives. The detector's own output over the same increment is the measurement record

$$
dY_t = 2\sqrt{\lambda}\,\langle\hat x\rangle_t\,dt + dW_t ,
$$

which is the conditional mean position $\langle\hat x\rangle$ buried in white noise. Equivalently, the record is signal (the mean) plus shot noise (the same $dW$), and because that one increment appears in both the state equation and the record, the state can be reconstructed from the observed current. That reconstruction is quantum filtering, and it is what a real experiment computes in real time. We will return to the filter concretely in Part IV.

Now the crucial reframing, and the answer to "where is the update?": the SSE is the $dt\to0$ limit of a concrete two-step update, and it is the update, not the differential form, that the simulation runs. Over one step the detector returns a single weak-measurement outcome $\bar x$, a reading of the position blurred by Gaussian noise, and the state is revised by **Bayes' rule**. The reading has likelihood

$$
P(\bar x \mid x) \propto e^{-2\lambda\,dt\,(x - \bar x)^2},
$$

a Gaussian in $\bar x$ centered on the true position $x$. Because the reading is noisy, the outcome is not drawn from $|\psi(x)|^2$ directly but from the *prior predictive* $P(\bar x) = \int P(\bar x\mid x)\,|\psi(x)|^2\,dx$, the density $|\psi(x)|^2$ convolved with the detector noise. And the conditional (posterior) density is the prior reweighted by the likelihood,

$$
|\psi(x)|^2 \;\longrightarrow\; \frac{P(\bar x\mid x)\;|\psi(x)|^2}{P(\bar x)} ,
$$

which is Bayes' rule verbatim: posterior $\propto$ likelihood $\times$ prior. At the level of amplitudes this reweighting is a Kraus operator $\hat M(\bar x)\propto e^{-\lambda\,dt\,(\hat x - \bar x)^2}$, the square root of the likelihood, applied and renormalized; multiplying each amplitude by $\sqrt{P(\bar x\mid x)}$ carries the phase along. The particle then evolves unitarily for the step, so one full step is Bayes' rule followed by Hamiltonian evolution:

$$
|\psi(t+dt)\rangle = \frac{e^{-i\hat H\,dt/\hbar}\;\hat M(\bar x)\,|\psi(t)\rangle}{\bigl\|\hat M(\bar x)\,|\psi(t)\rangle\bigr\|} .
$$

In short: measure by Bayes' rule, then evolve. Notice the one fact we will lean on for the rest of the notebook: the measurement half of this update knows nothing about $\hat H$. The outcome distribution, the Bayesian reweighting, and the renormalization are identical for a free particle, a trap, or a double well; only the unitary $e^{-i\hat H\,dt/\hbar}$ carries the potential. So generalizing the simulation to any $V(x)$ will mean generalizing that one factor and leaving the measurement untouched, which is exactly what the next Part does.

## Part II: One Function Carries the Potential: A Split-Operator Step on a Grid

Part I gave the update as an equation; to run it we have to make one decision that everything else follows from. A wavefunction $\psi(x)$ is defined at every point of the infinite line, which is more than a computer can hold, so we keep it only on a finite grid: $n$ equally spaced points across a window of length $L$, at spacing $dx=L/n$. The state becomes a list of $n$ complex numbers $\psi_j=\psi(x_j)$, and every operation below is arithmetic on that list. The care is all in choosing the window and the spacing, and in carrying, next to the positions, the momenta we will need. Let me explain each piece before writing the grid.

**Why a periodic box.** Cutting space to a window of length $L$ leaves two loose ends, and we glue them: the point $x=-L/2$ is identified with $x=+L/2$, so the segment becomes a ring of circumference $L$. We choose this periodic boundary for one reason. Momentum states are the plane waves $e^{ipx/\hbar}$, and on a ring such a wave is allowed only if it closes up smoothly after going once around, which forces the momentum onto the discrete ladder

$$
p_m=\frac{2\pi\hbar}{L}\,m,\qquad m=\ldots,-1,0,1,\ldots
$$

In other words, periodicity is exactly what turns momentum from a continuum into an evenly spaced list we can store beside the positions, and it makes those plane waves an exact basis. It costs nothing as long as $|\psi|^2$ stays negligible at the edges, so the packet never feels the seam where the ends are glued.

**How to choose the length and the spacing.** These answer two independent physical questions. The length $L$ must be large enough to hold the whole wavefunction, tails included, for the entire run: the packet spreads, drifts, and (under measurement) heats, and $L$ must contain all of it, or a tail that reaches one edge reappears at the other. The spacing $dx$ must be small enough to resolve the fastest wiggle of $\psi$, which is set by its largest momentum: a component of momentum $p$ oscillates with de Broglie wavelength $\lambda=2\pi\hbar/p$, and capturing a wave needs at least two samples across it, so $dx$ small enough that the biggest momentum you expect stays below $\pi\hbar/dx$. The two ends are mirror images: the box length fixes the momentum spacing $dp=2\pi\hbar/L$ (fine momentum needs a big box), the sampling fixes the momentum reach $\pi\hbar/dx$ (high momentum needs fine sampling). So choosing the position grid already fixes the momentum grid, which is why one function returns both `x` and `p`.

**What the Fourier transform is, and why we want it.** We will keep switching between two descriptions of the same state: how much amplitude sits at each position, the list $\psi_j$, and how much sits at each momentum, the coefficients $c_m$ when we write $\psi$ as a sum of the plane waves above, $\psi(x)=\sum_m c_m\,e^{ip_mx/\hbar}$. The recipe that turns the first list into the second is the *Fourier transform*: a reversible change of basis, the fixed sum $c_m=\tfrac{1}{\sqrt n}\sum_j \psi_j\,e^{-ip_mx_j/\hbar}$ that correlates the state against each plane wave, with an inverse that rebuilds $\psi$.

**How a wavefunction sits on the grid.** After all that, representing $\psi(x)$ is the easy part: evaluate it at the grid points and normalize, $\psi\mapsto(\psi(x_1),\dots,\psi(x_n))$ with $\sum_j|\psi_j|^2=1$. That list is the state in the position description; its Fourier transform is the same state in the momentum description. First fix units so that $\hbar=m=1$:

```wl
hbar = 1.; mass = 1.;
```

Define he grid given $n$ the space and $l$ the size of the box:

```wl
grid[n_, l_] := With[{dx = l/n},
   <|"n" -> n, "dx" -> dx, "x" -> dx Range[-n/2, n/2 - 1],
     "p" -> (2 Pi hbar/l) RotateLeft[Range[-n/2, n/2 - 1], n/2]|>];
```

The grid itself returns an association with $n$, the spacing $dx$, the positions $x_j=dx\{-n/2,\dots,n/2-1\}$ (centered on the origin, so the box is $[-L/2,L/2)$), and the quantized momenta $p_m=2\pi\hbar\,m/L$. The momenta are stored rotated, nonnegative $m$ first and the negative half moved to the end, because that is the order `Fourier` returns them in, so that later one multiplication lands each amplitude on its own momentum.

A convenient initial state is a Gaussian wave packet centered at $x_0$ with mean momentum $p_0$ and position width $s$. Following the recipe just above, we write its formula on the grid points `g["x"]` and normalize:

```wl
gaussian[g_, x0_, p0_, s_] := Normalize[Exp[-(g["x"] - x0)^2/(4 s^2) + I p0 g["x"]/hbar]];
```

The two numbers we will watch most are the conditional mean position $\langle\hat x\rangle=\int x\,|\psi(x)|^2\,dx$ and the conditional variance $\Sigma_x=\langle\hat x^2\rangle-\langle\hat x\rangle^2$. On the grid each is a dot product against the probability density $|\psi|^2$:

```wl
xMean[g_, ψ_] := g["x"] . Abs[ψ]^2;
xVar[g_, ψ_] := g["x"]^2 . Abs[ψ]^2 - xMean[g, ψ]^2;
```

For the measurement step, we implement it in the two pieces: first the outcome $\bar x$, then the Bayesian update of the state.

The outcome $\bar x$ is a sample from the prior predictive $P(\bar x)$. Recall that $P(\bar x)$ is $|\psi|^2$ convolved with Gaussian noise of variance $1/(4\lambda\,dt)$. So for the sampling of $\bar x$, we can treat it as $\bar x = X + \xi$, with $X$ distributed as $|\psi(x)|^2$ and $\xi$ distributed as a Gaussian with zero mean and standard deviation $1/\sqrt{4\lambda\,dt}$:

```wl
drawExact[g_, λ_, dt_][ψ_] := RandomChoice[Abs[ψ]^2 -> g["x"]] +
   RandomVariate[NormalDistribution[0, 1/Sqrt[4 λ dt]]];
```

Second the Bayesian update itself: reweight the amplitudes by the square root of the likelihood, $e^{-\lambda\,dt\,(x-\bar x)^2}$, and renormalize. That is the whole measurement step:

```wl
measUpdate[g_, λ_, dt_, xb_][ψ_] := Normalize[Exp[-λ dt (g["x"] - xb)^2] ψ];
```

Is this really Bayes' rule? It should be: the updated density $|\psi_{\rm post}(x)|^2$ ought to equal the prior $|\psi(x)|^2$ times the likelihood $e^{-2\lambda\,dt\,(x-\bar x)^2}$, normalized. Confirm the two constructions agree, for a test state and outcome:

```wl
Module[{gDemo = grid[512, 40.], ψt, post, bayes, λ = .75, dt = .01, xb = -.5},
 ψt = gaussian[gDemo, 1., 0.5, 2.];
 post = Abs[measUpdate[gDemo, λ, dt, xb][ψt]]^2;
 bayes = Normalize[Exp[-2 λ dt (gDemo["x"] - xb)^2] Abs[ψt]^2, Total];
 Max@Abs[post - bayes]]
```

The difference is at machine level, so `measUpdate` is exactly the Bayesian posterior; the operator form and the "likelihood times prior" form are two representations of the same update.

Two small diagnostics will be handy later: the norm error (how far $\||\psi\rangle\|$ strays from $1$) and the density matrix $|\psi\rangle\langle\psi|$ built from a state vector. Define both:

```wl
normError[ψ_] := Abs[Norm[ψ] - 1];
projector[ψ_] := KroneckerProduct[ψ, Conjugate[ψ]];
```

For $\hat H = \hat p^2/2m + V(\hat x)$ the potential is diagonal in position, so its propagator $e^{-iV(\hat x)\,dt/\hbar}$ is just a pointwise phase. Define aa half-step of it:

```wl
potentialStep[g_, Vfun_, dt_][ψ_] := Exp[-I Vfun[g["x"]] dt/(2 hbar)] ψ;
```

Let's explain why the half-step. To approximate one step of the true dynamics, $U(dt)=e^{-i\hat H\,dt/\hbar}$ with $\hat H=\hat T+V$, we use the Trotter-Suzuki formula: a family of operator-splitting methods for $e^{h(A+B)}$ when $e^{hA}$ and $e^{hB}$ are each easy to compute but $A$ and $B$ do not commute. Here $A=-i\hat T/\hbar$ and $B=-iV/\hbar$ are the kinetic and potential generators and $h=dt$. The first-order member of the family of operator-splitting is the Lie-Trotter splitting, $S_1(h)=e^{hA}e^{hB}$, with local error $O(h^2)$ per step and global error $O(h)$ over a fixed total time. The symmetric second-order member is the Strang splitting, $S_2(h)=e^{hB/2}e^{hA}e^{hB/2}$, with local error $O(h^3)$ and global error $O(h^2)$. More generally, higher-order Suzuki compositions give an order-$p$ method with local error $O(h^{p+1})$ and global error $O(h^p)$, the one lost power coming from the roughly $1/h$ steps needed to cover a fixed interval. We want a global error of order $dt^2$, so we take $S_2$, the Strang splitting: with $h=dt$ it is $e^{-iV\,dt/2\hbar}\,e^{-i\hat T\,dt/\hbar}\,e^{-iV\,dt/2\hbar}$, a half potential kick around a full kinetic drift, and it is exact whenever $[A,B]=0$.

```wl
unitaryStep[g_, dt_][ψ_] := InverseFourier[
   Exp[-I g["p"]^2 dt/(2 mass hbar)] Fourier[ψ, FourierParameters -> {0, -1}],
   FourierParameters -> {0, -1}];
evolveStep[g_, Vfun_, dt_] :=
   potentialStep[g, Vfun, dt] @* unitaryStep[g, dt] @* potentialStep[g, Vfun, dt];
```

One full conditional step is now described as follows: draw an outcome, apply the Bayesian `measUpdate`, then `evolveStep`. Define the step:

```wl
stepV[g_, Vfun_, λ_, dt_, xb_][ψ_] := evolveStep[g, Vfun, dt][measUpdate[g, λ, dt, xb][ψ]];
```

A whole experimental run is that step iterated `nt` times from an initial state, with a fresh outcome drawn at each step. One runner does it: `trajectory` returns the sequence of conditional states and `Sow`s each measurement outcome as it goes, so a bare call gives just the states, while wrapping it in `Reap` hands back the states together with the record. Define it:

```wl
trajectory[g_, Vfun_, λ_, dt_, nt_][ψ0_] :=
   NestList[stepV[g, Vfun, λ, dt, Sow@drawExact[g, λ, dt][#]][#] &, ψ0, nt];
```

That is the entire engine. Everything below is these functions with a different `Vfun`, a different initial state, or a different measurement strength.

## Part III: The Free Particle: The V=0 Case, Recovered and Heated

Let's exercise the engine on the case we already understand analytically, the free particle, by passing `0 &` for the potential. The point is twofold: to see the generalized code reproduce the known free-particle physics, and to build the unconditional (master-equation) machinery we will reuse in the trap and the well.

The first known fact: under measurement a broad packet localizes to a steady conditional width $(\hbar/4\lambda m)^{1/4}$. Run one trajectory from a wide (width $2$) packet and track its conditional width over the run:

```wl
gF = grid[256, 50.]; λR = 1.; dt = 0.005; nt = 600;
statesF = BlockRandom[SeedRandom[11]; trajectory[gF, 0 &, λR, dt, nt][gaussian[gF, 0., 0., 2.]]];
widthsF = Sqrt[xVar[gF, #] & /@ statesF];
```

Plot the width and show it approaches the predicted steady value, $(\hbar/4\lambda m)^{1/4}$:

```wl
ListLinePlot[widthsF, PlotRange -> {.5, All}, DataRange -> {0, nt dt}, Frame -> True,
 FrameLabel -> {"time t", "conditional width"},
 PlotLabel -> "A broad packet localizes to the steady width (dashed)",
 Epilog -> {Dashed, Line[{{0, (hbar/(4 λR mass))^(1/4)}, {nt dt, (hbar/(4 λR mass))^(1/4)}}]}]
```

And the substantive guarantee of the split-operator step is that the norm still sits at $1$:

```wl
normError[Last@statesF]
```

At unit detection efficiency the conditional state stays pure, which is what lets us carry it as a single state vector; the state only becomes mixed when we average over the noise, which is the master equation, coming next.

To average, we need the unconditional state as a density matrix $\rho$. Discarding the record and averaging the conditional states over the Wiener noise gives the Lindblad master equation, which integrates by the same operator split: the measurement dissipator damps each element $\rho_{ij}$ by the factor $e^{-\frac{\lambda}{2}(x_i-x_j)^2 dt}$, and the Hamiltonian acts as `evolveStep` turned into a matrix. Define the damping factor first:

```wl
decohereMat[g_, λ_, dt_] := Exp[-(λ/2) Outer[(#1 - #2)^2 &, g["x"], g["x"]] dt];
```

Turn `evolveStep` into a matrix by applying it to each basis vector, and assemble one master-equation step as "damp, then evolve on both sides":

```wl
uniMat[g_, Vfun_, dt_] := Transpose[evolveStep[g, Vfun, dt] /@ IdentityMatrix[g["n"]]];
lindblad[g_, Vfun_, λ_, dt_] := With[{u = uniMat[g, Vfun, dt], dec = decohereMat[g, λ, dt]},
   Function[ρ, u . (ρ dec) . ConjugateTranspose[u]]];
```

Position observables read straight off $\rho$'s diagonal, but momentum observables like $\langle\hat p^2\rangle$ need the momentum-squared operator as a matrix, built spectrally:

```wl
p2Op[g_] := Transpose[InverseFourier[g["p"]^2 Fourier[#, FourierParameters -> {0, -1}],
      FourierParameters -> {0, -1}] & /@ IdentityMatrix[g["n"]]];
```

Now the second known fact, the measurement backaction. Propagate the free master equation from a pure packet and watch its momentum variance $\langle\hat p^2\rangle(t)$ climb to the predicted final value:

```wl
gS = grid[32, 16.];
dtS = 0.0005; tf = 0.4; ntS = Round[tf/dtS];
ψ0 = gaussian[gS, 0., 0., 1.];
ρfree = NestList[lindblad[gS, 0 &, λR, dtS], projector[ψ0], ntS];
With[{pred = Re[Conjugate[ψ0] . p2Op[gS] . ψ0] + λR hbar^2 tf},
 ListLinePlot[Re@Tr[# . p2Op[gS]] & /@ ρfree, DataRange -> {0, tf}, Frame -> True,
  FrameLabel -> {"time t", "⟨p²⟩"}, PlotLabel -> "Measurement heats momentum linearly (predicted value dashed)",
  Epilog -> {Dashed, Line[{{0, pred}, {tf, pred}}]}]]
```

The momentum variance rises exactly as $\langle\hat p^2\rangle_0+\lambda\hbar^2 t$: position measurement feeds momentum at the constant rate $\lambda\hbar^2$, the price the uncertainty principle charges for position information. That rate is state- and Hamiltonian-independent, so for the free particle, whose $\hat H$ commutes with $\hat p^2$, it is the whole story. In a potential it is not: the force couples $\langle\hat p^2\rangle$ to $\langle\hat x^2\rangle$ as the packet swings, which is why Part V follows the total energy instead. This backaction is the one piece of physics common to every case below; we meet it again competing with the trap and with tunneling.

As another example, let's consider a delocalized initial state. Take a "cat," a superposition of two packets sitting at $x=\pm6$. Build it and look at its two lobes:

```wl
gC = grid[256, 50.];
ψCat = Normalize[gaussian[gC, -6., 0., 1.] + gaussian[gC, 6., 0., 1.]];
ListLinePlot[ψCat, DataRange -> {gC["x"][[1]], gC["x"][[-1]]}, Frame -> True,
 FrameLabel -> {"x", "ψ(x)"}, PlotLabel -> "A cat state: two packets at x = ±6"]
```

Its initial variance is large, dominated by that $\pm6$ lobe separation:

```wl
xVar[gC, ψCat]
```

Now watch the measurement collapse it. Take a moderate strength $\lambda=0.1$, chosen so the two lobes coexist for a visible stretch before one wins rather than collapsing in the first few steps, and run one trajectory:

```wl
statesC = BlockRandom[SeedRandom[5]; trajectory[gC, 0 &, 0.1, 0.005, 600][ψCat]];
```

See the collapse as a surface of the conditional density $|\psi(x,t)|^2$ over position and time:

```wl
band = Flatten@Position[UnitStep[10 - Abs[gC["x"]]], 1];
ListPlot3D[(Abs[statesC]^2)[[1 ;; ;; 6, band]], PlotRange -> All,
 DataRange -> {gC["x"][[{First@band, Last@band}]], {0, 3.}},
 ColorFunction -> (ColorData["DarkRainbow"][#3] &), Mesh -> None,
 AxesLabel -> {"x", "t", "|ψ|²"}, PlotLabel -> "A cat collapses to one packet"]
```

The two ridges compete and one grows while the other fades, with no potential in sight.

That surface is the particle's actual density, which no experiment can see. What an observer sees is the record, the stream of noisy outcomes $\bar x_k$, and from it alone they must infer which lobe won. That inference is a running bet between two hypotheses, the left lobe at $x_L=-6$ against the right at $x_R=+6$. Each outcome has likelihood $\propto e^{-2\lambda\,dt(\bar x_k-x_h)^2}$ under hypothesis $h$, so the log-odds for left over right accumulate one term at a time: each term is $-2\lambda\,dt$ times $(\bar x_k+6)^2-(\bar x_k-6)^2$, which is exactly $24\,\bar x_k$, so the increment is simply $-48\lambda\,dt\,\bar x_k$. Rerun the same trajectory, keep its outcomes, and accumulate them into the posterior `PleftRec`; first look at the raw stream the detector logs:

```wl
{catStates, {catOuts}} = BlockRandom[SeedRandom[5]; Reap[trajectory[gC, 0 &, 0.1, 0.005, 600][ψCat]]];
PleftRec = LogisticSigmoid[Accumulate[-48 (0.1) (0.005) catOuts]];
ListLinePlot[catOuts, PlotRange -> All, Frame -> True,
 FrameLabel -> {"measurement step k", "outcome"}, PlotLabel -> "The raw record: signal buried in noise"]
```

That raw stream is almost pure noise: each reading scatters with standard deviation $1/\sqrt{4\lambda\,dt}$, far wider than the $\pm6$ lobe separation, so no single outcome names the lobe. The verdict lives only in the slow average that `PleftRec` accumulates, which the grid's own final $\langle\hat x\rangle$ should confirm, ending negative for the left lobe:

```wl
xMean[gC, Last@catStates]
```

So the record says *left lobe* and the grid agrees. The road there is the instructive part: the belief builds only over many steps, and can even lean the wrong way first. Plot the posterior forming over the first many steps, from the even prior $\tfrac12$:

```wl
ListLinePlot[Prepend[PleftRec, 0.5][[1 ;; 120]], DataRange -> {0, 119}, Frame -> True, ImageSize -> 500,
 AspectRatio -> 1/2, PlotRange -> {0, 1}, GridLines -> {None, {{0.5, Gray}}},
 FrameLabel -> {"measurement step k", "observer's posterior \!\(\*SubscriptBox[\(P\), \(left\)]\)"},
 PlotLabel -> "The record decides the lobe (here → 1, the left lobe won)"]
```

$P_{\rm left}$ starts at $\tfrac12$ and fluctuates, at first even favoring the wrong lobe, then converges to $1$ as the accumulating record raises the signal-to-noise ratio: the state has localized to the left lobe. The observer reaches this from the measurement record alone, never from the wavefunction.

That verdict was this record's; a different noise history would decide differently. The initial cat is symmetric, so nothing favors either lobe: the measurement forces a definite choice, but which lobe is set entirely by the noise. Overlay several independent seeded records and watch their posteriors split, some pinning to $1$ (left), some to $0$ (right):

```wl
ListLinePlot[
 Table[BlockRandom[SeedRandom[k];
    LogisticSigmoid[Accumulate[-48 (0.1) (0.005)
        Reap[trajectory[gC, 0 &, 0.1, 0.005, 100][ψCat]][[2, 1]]]]], {k, 6}],
 PlotRange -> {0, 1.1}, Frame -> True, GridLines -> {None, {{0.5, Gray}}},
 FrameLabel -> {"measurement step k", "\!\(\*SubscriptBox[\(P\), \(left\)]\)"}, PlotLabel -> "Different records, different verdicts"]
```

The runs split between the two lobes, each driven to a definite outcome by its own record. This is the measurement record doing the observer's inference, the same accumulation a real detector performs in software. Now we turn on a potential and ask what survives.

## Part IV: The Harmonic Trap: Measurement That Squeezes Below the Vacuum

A trapped particle is never still: even in its ground state its position keeps the irreducible zero-point spread $\Sigma_x^{zp}=\hbar/2m\omega$, the blur of the quantum vacuum. So here is the question this Part answers: **watch that position continuously, and can you come to know it more sharply than the vacuum itself allows? What does the state become?**

The harmonic trap is where this can be settled exactly. Because $\hat H$ is quadratic, a Gaussian conditional state stays Gaussian under both the measurement update and the evolution, so the entire state collapses to five numbers, the two means $\langle\hat x\rangle,\langle\hat p\rangle$ and the three covariances $\Sigma_x, C, \Sigma_p$, and we can solve their steady state with algebra before ever touching the grid. (This is also the real setting of cavity and levitated optomechanics, where light monitors a mechanical mode.)

Start with the covariances, and the fact that makes the trap solvable: they close among themselves. The rates for $\Sigma_x$, $C=\tfrac12\langle\{\hat x,\hat p\}\rangle-\langle\hat x\rangle\langle\hat p\rangle$, and $\Sigma_p$ depend only on each other, never on the means or the noise, so they run deterministically (a matrix Riccati flow). Write the three rates, keeping $\lambda,\omega,m,\hbar$ symbolic, and solve for the fixed point where all three vanish:

```wl
ClearAll[Σx, Cc, Σp, λ, ω, m, ℏ];
riccati = {2 Cc/m - 4 λ Σx^2, Σp/m - m ω^2 Σx - 4 λ Σx Cc, -2 m ω^2 Cc + λ ℏ^2 - 4 λ Cc^2};
steady = Solve[Thread[riccati == 0], {Σx, Cc, Σp}];
```

`Solve` returns several algebraic roots; the physical one has $\Sigma_x>0$ and $\Sigma_p>0$. Select it by a quick numeric probe and simplify the result symbolically:

```wl
posRoot[rule_] := With[{v = {Σx, Cc, Σp} /. rule /. {λ -> 1., ω -> 1., m -> 1., ℏ -> 1.}},
   Max@Abs@Im[v] < 10^-10 && Re[v[[1]]] > 0 && Re[v[[3]]] > 0];
{ΣxSS, CcSS, ΣpSS} = Simplify[{Σx, Cc, Σp} /. First@Select[steady, posRoot], {λ > 0, ω > 0, m > 0, ℏ > 0}]
```

The first entry is the steady position variance $\Sigma_x^{ss}=\dfrac{1}{2\sqrt2}\sqrt{\dfrac{-m\omega^2+\sqrt{4\hbar^2\lambda^2+m^2\omega^4}}{\lambda^2 m}}$. Two checks that it is right. First, perfect measurement keeps the state pure, and a pure Gaussian sits exactly at the uncertainty floor, so $\Sigma_x\Sigma_p-C^2$ must come out to $\hbar^2/4$:

```wl
Simplify[ΣxSS ΣpSS - CcSS^2, {λ > 0, ω > 0, m > 0, ℏ > 0}]
```

And switching off the trap, $\omega\to0$, should return the free-particle steady value $\sqrt{\hbar/4\lambda m}$:

```wl
Simplify[Limit[ΣxSS, ω -> 0], {λ > 0, m > 0, ℏ > 0}]
```

Both check out: the steady state is a pure Gaussian that reduces to the free particle when the trap is switched off. Now the payoff, the question we opened with. Set the steady variance against the vacuum value $\Sigma_x^{zp}=\hbar/2m\omega$: form the ratio $\Sigma_x^{ss}/\Sigma_x^{zp}$ and simplify:

```wl
μ = 2 ℏ λ/(m ω^2);
Simplify[ΣxSS/(ℏ/(2 m ω)), {λ > 0, ω > 0, m > 0, ℏ > 0}]
```

In the single dimensionless measurement strength $\mu=2\hbar\lambda/m\omega^2$ that is exactly the clean closed form $\sqrt{2(\sqrt{1+\mu^2}-1)}/\mu$; confirm the two are equal:

```wl
Simplify[ΣxSS/(ℏ/(2 m ω)) == Sqrt[2 (Sqrt[1 + μ^2] - 1)]/μ, {λ > 0, ω > 0, m > 0, ℏ > 0}]
```

That ratio is below $1$ for every $\mu>0$. Its two limiting regimes follow by expanding the closed form (writing $q$ for $\mu$ as a bare expansion variable):

```wl
Block[{q}, Simplify[Normal@Series[Sqrt[2 (Sqrt[1 + q^2] - 1)]/q, {q, 0, 2}], q > 0]]
```

That is the weak limit, dipping just below the zero-point value as $\mu$ grows from zero. For strong measurement, $\mu\to\infty$, the ratio itself vanishes; the telling quantity is $\sqrt\mu$ times it, which should approach a constant:

```wl
Block[{q}, Limit[Sqrt[q] Sqrt[2 (Sqrt[1 + q^2] - 1)]/q, q -> Infinity]]
```

So the ratio is $\approx1-\mu^2/8$ for weak watching and $\approx\sqrt{2/\mu}$ for strong. There is the answer to our question, stated plainly: continuous position measurement pushes the conditional variance below the vacuum at *any* strength $\mu>0$, and deeper the harder you watch. The conditional state is a pure squeezed state, narrower in position than the ground state itself, because the measurement is steadily extracting position information. Fix a strongly-measured regime, $\omega=1$ and $\lambda=2$ (so $\mu=4$), and read the position variance in units of the vacuum:

```wl
ωR = 1.; λR = 2.;
ΣxZP = hbar/(2 mass ωR); ΣpZP = hbar mass ωR/2;
harmNum = {ℏ -> hbar, m -> mass, ω -> ωR, λ -> λR};
(ΣxSS /. harmNum)/ΣxZP
```

Below $1$: squeezed below the vacuum. The price shows in the momentum, anti-squeezed above its own zero-point value by the reciprocal cost:

```wl
(ΣpSS /. harmNum)/ΣpZP
```

The position is squeezed well below the zero-point variance and the momentum is anti-squeezed above it: the backaction paid for the position information. All of this came from the closed-form moments. Does the grid engine, which assumes nothing about Gaussianity, agree? Evolve one trajectory from the ground state in the trap and read the late-time conditional variance:

```wl
Vharm = 0.5 mass ωR^2 #^2 &;
gH = grid[256, 24.]; dtH = 0.002; tfH = 8.; ntH = Round[tfH/dtH];
runH = BlockRandom[SeedRandom[11]; trajectory[gH, Vharm, λR, dtH, ntH][gaussian[gH, 0., 0., Sqrt[ΣxZP]]]];
ΣxH = xVar[gH, #] & /@ runH;
Mean[ΣxH[[Round[0.6 ntH] ;;]]]
```

That late-time grid variance should match the closed-form Riccati steady value and sit below the zero-point $\Sigma_x^{zp}=\tfrac12$; read the two reference values side by side:

```wl
{ΣxSS /. harmNum, ΣxZP}
```

The grid variance relaxes to the Riccati steady value, matching it to grid resolution and sitting well below the zero-point. The relaxation is the same in every run, because the conditional variance obeys a deterministic (noiseless) equation while only the mean is driven by the record. Plot the variance against time with the zero-point and Riccati lines drawn in:

```wl
ListLinePlot[{Transpose[{dtH Range[0, ntH], ΣxH}]},
 Frame -> True, ImageSize -> 500, AspectRatio -> 1/2, PlotRange -> {0, 0.55},
 GridLines -> {None, {{ΣxZP, Directive[Gray, Dashed]}, {ΣxSS /. harmNum, Directive[Red, Dashed]}}},
 FrameLabel -> {"time t", "conditional variance Σx"},
 PlotLabel -> "Measurement squeezes Σx below the zero-point line (gray) to the Riccati value (red)"]
```

The grid agreed, but it was never needed: the whole conditional state is those five numbers, the two means driven by the record and the three covariances running on their own. This five-number object is the quantum Kalman-Bucy filter, the algorithm an optomechanics experiment integrates in real time to track its mechanical mode from the measurement current. Hand its equations to the stochastic integrator `ItoProcess`, the drifts in place and the noise entering only through the means:

```wl
ClearAll[mx, mp, Σx, Cc, Σp];
momentProc = ItoProcess[
   {{mp[t]/mass, -mass ωR^2 mx[t], 2 Cc[t]/mass - 4 λR Σx[t]^2,
     Σp[t]/mass - mass ωR^2 Σx[t] - 4 λR Σx[t] Cc[t], -2 mass ωR^2 Cc[t] + λR hbar^2 - 4 λR Cc[t]^2},
    {{2 Sqrt[λR] Σx[t]}, {2 Sqrt[λR] Cc[t]}, {0}, {0}, {0}},
    {mx[t], mp[t], Σx[t], Cc[t], Σp[t]}},
   {{mx, mp, Σx, Cc, Σp}, {0., 0., ΣxZP, 0., ΣpZP}}, {t, 0}];
```

Integrate an ensemble of these filters and read the mean conditional variance across the ensemble together with its run-to-run spread:

```wl
tdM = BlockRandom[SeedRandom[7]; RandomFunction[momentProc, {0., tfH, 0.01}, 300,
    Method -> "StochasticRungeKuttaScalarNoise"]];
finM = tdM["ValueList"][[All, -1]];
{Mean[#[[3]] & /@ finM], StandardDeviation[#[[3]] & /@ finM]}
```

The mean sits on the Riccati value and the spread is exactly zero: every realization carries the *same* variance. Contrast the conditional mean, which spreads freely across runs because it is the part the record drives:

```wl
StandardDeviation[#[[1]] & /@ finM]
```

Grid, filter, and algebra return the same variance, and across the ensemble that variance has exactly zero spread while the mean wanders freely from run to run. That is the division of labor made visible: the covariance is fixed and deterministic, and all the randomness lives in the means, fed by the record through the gains $2\sqrt{\lambda}\,\Sigma_x$ and $2\sqrt{\lambda}\,C$. Five numbers and a 256-point grid describe the same state, because for a quadratic trap the Gaussian picture is exact. One honest caveat from the lab: we assumed perfect detection, which is what keeps the state pure and below the vacuum. A real experiment collects only a fraction $\eta<1$ of the light, mixing the conditional state and lifting its uncertainty back above the vacuum (Magrini et al. measure $1.3$ times zero-point), the extension we return to at the close.

### Closing the Loop: The Record Reconstructs the State

The filter above ran on its own manufactured noise, so it showed that the covariances agree but never that the *same* record rebuilds the *same* conditional mean. That reconstruction, the observer recovering the state from the current, is the whole point of quantum filtering and the promise Part I made here. Deliver it in two steps: first confirm the record carries exactly the noise the theory claims, then feed that record to the filter and watch it track the grid.

Run one long trajectory in the trap, this time keeping every measurement outcome $\bar x_k$ next to the conditional mean $\langle\hat x\rangle_k$ the grid computes on the spot:

```wl
{recStates, {recOuts}} = BlockRandom[SeedRandom[21]; Reap[trajectory[gH, Vharm, λR, dtH, ntH][gaussian[gH, 0., 0., Sqrt[ΣxZP]]]]];
xGrid = xMean[gH, #] & /@ recStates;
```

The record is those outcomes assembled into a current. Over one step the integrated signal is $\Delta Y_k=2\sqrt\lambda\,dt\,\bar x_k$, and its *innovation*, the part no filter could have predicted, is what is left after subtracting the expected signal $2\sqrt\lambda\,\langle\hat x\rangle\,dt$:

$$
\Delta W_k=\Delta Y_k-2\sqrt\lambda\,\langle\hat x\rangle_k\,dt=2\sqrt\lambda\,dt\,(\bar x_k-\langle\hat x\rangle_k).
$$

The theory is specific about $\Delta W$: it should be the bare Wiener increment $dW$, mean zero, variance $dt$, uncorrelated from one step to the next. The variance follows from $\bar x=\langle\hat x\rangle+\eta$ with $\eta$ of variance $1/(4\lambda\,dt)$, giving $4\lambda\,dt^2\cdot 1/(4\lambda\,dt)=dt$. First its mean and variance: the mean should read near $0$ and the variance-over-$dt$ ratio near $1$.

```wl
ΔY = 2 Sqrt[λR] dtH recOuts;
ΔW = ΔY - 2 Sqrt[λR] Most[xGrid] dtH;
{Mean[ΔW], Variance[ΔW]/dtH}
```

And successive steps should be uncorrelated; the lag-1 through lag-4 autocorrelations should all sit near zero:

```wl
Table[Correlation[Drop[ΔW, -k], Drop[ΔW, k]], {k, 4}]
```

The innovation comes out with mean consistent with zero, variance consistent with $dt$ to within the $\pm\sqrt{2/N}$ sampling scatter of a finite run, and lag autocorrelations that all sit inside the $\pm1/\sqrt N$ band of true white noise. The record is signal plus white noise, exactly as posed. (At finite $dt$ the intrinsic conditional spread lifts the variance by a further $4\lambda\,dt\,\Sigma_x$, a fraction that is below the sampling floor here and vanishes as $dt\to0$.) Two pictures make the point: scaled by $\sqrt{dt}$ the innovation histogram is a standard normal, and its autocorrelation is flat.

```wl
Show[Histogram[ΔW/Sqrt[dtH], 60, "PDF"], Plot[PDF[NormalDistribution[], z], {z, -4, 4}, PlotStyle -> Red],
 Frame -> True, ImageSize -> 460, AspectRatio -> 1/2, FrameLabel -> {"ΔW / √dt", "density"},
 PlotLabel -> "Innovations are standard-normal (red curve)"]
```

```wl
ListPlot[Table[Correlation[Drop[ΔW, -k], Drop[ΔW, k]], {k, 1, 20}], Filling -> Axis, PlotRange -> {-0.05, 0.05},
 Frame -> True, ImageSize -> 460, AspectRatio -> 1/2, GridLines -> {None, 2/Sqrt[Length[ΔW]] {1, -1}},
 FrameLabel -> {"lag", "autocorrelation"}, PlotLabel -> "No structure across lags (±2/√N band)"]
```

Now the reconstruction. Hand the filter the recorded increments $\Delta Y_k$ and nothing else, integrating the mean equations

$$
dm_x=\tfrac{m_p}{m}\,dt+2\sqrt\lambda\,\Sigma_x\,dW_f,\qquad dm_p=-m\omega^2 m_x\,dt+2\sqrt\lambda\,C\,dW_f,
$$

in which the filter forms its own innovation $dW_f=\Delta Y_k-2\sqrt\lambda\,m_x\,dt$ from its current estimate $m_x$, and the covariances $\Sigma_x,C,\Sigma_p$ follow the same deterministic Riccati flow as before, independent of the record. Write the covariance drift and one filter step:

```wl
covDrift[Σx_, Cc_, Σp_] := {2 Cc/mass - 4 λR Σx^2, Σp/mass - mass ωR^2 Σx - 4 λR Σx Cc, -2 mass ωR^2 Cc + λR hbar^2 - 4 λR Cc^2};
filterStep[{mx_, mp_, Σx_, Cc_, Σp_}, ΔYk_] := With[{dWf = ΔYk - 2 Sqrt[λR] mx dtH, cov = {Σx, Cc, Σp} + dtH covDrift[Σx, Cc, Σp]},
   {mx + mp/mass dtH + 2 Sqrt[λR] Σx dWf, mp - mass ωR^2 mx dtH + 2 Sqrt[λR] Cc dWf, cov[[1]], cov[[2]], cov[[3]]}];
```

Fold that step along the recorded increments from the ground state, and compare the filter's mean to the grid's:

```wl
mFilter = FoldList[filterStep, {0., 0., ΣxZP, 0., ΣpZP}, ΔY][[All, 1]];
Max@Abs[mFilter - xGrid]
```

That maximum deviation only means something against the range over which $\langle\hat x\rangle$ itself swings, so read the swing too:

```wl
{Min@xGrid, Max@xGrid}
```

Driven by the recorded current alone, the five-number filter reproduces the grid's conditional mean to a small fraction of the range over which $\langle\hat x\rangle$ itself swings: the residual is set by the grid spacing and the step size, not by anything statistical. The same record the grid produced, fed to a filter that never sees the wavefunction, rebuilds the same trajectory of the mean. That is quantum filtering in one line of physics: the experiment measures the current, the filter integrates five numbers in real time, and out comes the conditional state. Overlay the two and they lie on top of each other:

```wl
ListLinePlot[{Transpose[{dtH Range[0, ntH], xGrid}], Transpose[{dtH Range[0, ntH], mFilter}]},
 PlotStyle -> {Directive[Thick, GrayLevel[0.6]], Directive[Red, Dashed]}, Frame -> True, ImageSize -> 520,
 AspectRatio -> 1/2, PlotLegends -> {"grid ⟨x⟩", "filter m_x from the record"},
 FrameLabel -> {"time t", "conditional mean position"},
 PlotLabel -> "One record, two engines: the filter reconstructs the grid's ⟨x⟩(t)"]
```

## Part V: Backaction in the Trap: An Unconditional State That Never Settles

The conditional state settled into a steady squeezed pure state. It is natural to expect the unconditional state, the average over the record, to settle too, perhaps into a thermal state. It does not, and seeing why is worth a moment. The trap conserves energy on its own, but the measurement backaction keeps injecting energy the closed trap has no way to remove. Watch this in the moments: the unconditional second moments close by themselves, and their energy has a constant source. Solve the unconditional moment equations in closed form and read the energy $E(t)$:

```wl
ClearAll[XX, SS, PP];
usol = First@DSolve[{XX'[t] == 2 SS[t]/m, SS'[t] == PP[t]/m - m ω^2 XX[t],
    PP'[t] == -2 m ω^2 SS[t] + λ ℏ^2, XX[0] == x0, SS[0] == 0, PP[0] == p0}, {XX, SS, PP}, t];
energy = Simplify[PP[t]/(2 m) + (1/2) m ω^2 XX[t] /. usol]
```

It is linear in $t$, so its time derivative should be a constant, the steady heating rate:

```wl
Simplify[D[energy, t]]
```

The energy is exactly linear in time, $E(t)=E_0+\lambda\hbar^2 t/2m$, with no oscillatory part: the breathing of $\langle\hat x^2\rangle$ and $\langle\hat p^2\rangle$ against each other cancels in the sum, leaving only steady heating at $dE/dt=\lambda\hbar^2/2m$. Equivalently, the mean occupation climbs linearly, $dn/dt=\lambda\hbar/2m\omega$: the unconditional state is a mixed Gaussian whose energy and mean occupation grow without bound. With no damping and no detailed balance there is no thermodynamic temperature to assign, so it is thermal-like only in its Gaussian shape, not a thermal state. Confirm it on the grid by propagating the master equation and reading the energy and the individual moments at a few times. First two observables on $\rho$, position-squared as a diagonal expectation and the energy:

```wl
gHL = grid[96, 24.]; dtHL = 0.002;
x2Mean[g_, ρ_] := Re[Diagonal[ρ] . g["x"]^2];
Emeas[ρ_] := Re@Tr[ρ . p2Op[gHL]]/(2 mass) + 0.5 mass ωR^2 x2Mean[gHL, ρ];
```

Now propagate from the ground state and tabulate:

```wl
ρ0HL = projector[gaussian[gHL, 0., 0., Sqrt[ΣxZP]]];
tsHL = {0., 0.5, 1., 1.5, 2.};
ρsHL = FoldList[Nest[lindblad[gHL, Vharm, λR, dtHL], #1, Round[(#2[[2]] - #2[[1]])/dtHL]] &,
   ρ0HL, Partition[tsHL, 2, 1]];
TableForm[{Emeas /@ ρsHL, x2Mean[gHL, #] & /@ ρsHL, Re@Tr[#.p2Op[gHL]] & /@ ρsHL, Re@Tr[#.#] & /@ ρsHL},
 TableHeadings -> {{"E", "⟨x²⟩", "⟨p²⟩", "purity"}, tsHL}]
```

As predicted, the energy climbs in equal steps over $0\le t\le2$, a straight ramp of slope exactly $\lambda\hbar^2/2m$; both $\langle\hat x^2\rangle$ and $\langle\hat p^2\rangle$ grow, and the purity falls from $1$ as the unconditional state mixes. There is no steady state to plot, only a ramp:

```wl
ListLinePlot[Transpose[{tsHL, Emeas /@ ρsHL}],
 Frame -> True, ImageSize -> 460, AspectRatio -> 1/2, PlotRange -> {0, All}, Mesh -> All,
 FrameLabel -> {"time t", "unconditional energy E(t)"},
 PlotLabel -> "The trapped unconditional state heats linearly: dE/dt = λℏ²/2m"]
```

Note the contrast with Part IV, which is starker here than for the free particle: one realization sits forever in a fixed pure squeezed packet, while the average over realizations heats without bound. A genuine thermal *steady* state needs a mechanical damping channel to balance the heating, which is exactly why the optomechanics experiments cool with feedback rather than watch the bare measurement; restoring that damping is a one-line addition to the master equation and the natural next extension.

Before we move on, let us summarize what the two quadratic cases, free and harmonic, have taught us:

- The measurement half of the update is Bayes' rule, $|\psi(x)|^2\to P(\bar x\mid x)|\psi(x)|^2/P(\bar x)$, and it is independent of the Hamiltonian; only the unitary carries the potential, recovered exactly at $V=0$.
- Continuous position measurement heats the momentum at the state-independent rate $\lambda\hbar^2$, $\langle\hat p^2\rangle(t)=\langle\hat p^2\rangle_0+\lambda\hbar^2 t$.
- In the harmonic trap the conditional state is a pure squeezed Gaussian with $\Sigma_x^{ss}/\Sigma_x^{zp}=\sqrt{2(\sqrt{1+\mu^2}-1)}/\mu<1$, below the zero-point for every $\mu>0$.
- Because $\hat H$ is quadratic the Gaussian moments close, so the five-number Kalman-Bucy filter is exact and agrees with the 256-point grid.
- The unconditional trapped state has no steady state; it heats linearly at $dE/dt=\lambda\hbar^2/2m$.

Every one of these rested on the moments closing. The next Part removes that.

## Part VI: The Double Well: Where the Moments Fail and the Grid Wins

Introduce a quartic potential and the Gaussian description breaks. For a double well $V(x)=b(x^2-a^2)^2$ the force is cubic, so the equation for the mean momentum needs the third moment, which a non-Gaussian state does not fix from its mean and variance. See it directly by differentiating the potential. The double-well force is cubic:

```wl
ClearAll[x, a, b];
D[b (x^2 - a^2)^2, x]
```

So $d\langle\hat p\rangle/dt=-\langle V'\rangle=-4b(\langle\hat x^3\rangle-a^2\langle\hat x\rangle)$ pulls in the third moment $\langle\hat x^3\rangle$. The harmonic force, by contrast, is linear and closes the hierarchy at the first moment:

```wl
D[(1/2) m ω^2 x^2, x]
```

Where the harmonic force $m\omega^2 x$ closes the hierarchy at the first moment, the cubic force $4b(x^3-a^2x)$ pulls in $\langle\hat x^3\rangle$, which pulls in $\langle\hat x^4\rangle$, and no finite set of moments closes. There is no five-number filter here; only the grid split-operator scheme, which never assumes a shape for the conditional state, applies. Set up the well, with minima at $x=\pm a$ and a barrier between them, and the mask that measures how much probability sits in the left well. Weight the single barrier-top point $x=0$ by one half (it belongs to neither well), so a left-right symmetric state registers exactly $P_{\rm left}=\tfrac12$ instead of missing that point; `(1 - Sign[x])/2` does this in one Listable stroke, returning $1$ for $x<0$, $\tfrac12$ at $x=0$, and $0$ for $x>0$:

```wl
aa = 1.2; bb = 1.0; Vdw = bb (#^2 - aa^2)^2 &;
leftMask[g_] := (1 - Sign[g["x"]])/2;
```

We will want the well's two lowest eigenstates on two different grids, so package the diagonalization once as a helper that returns the Hamiltonian, the two lowest energies, the symmetric and antisymmetric states, and their one-well combination:

```wl
wellDoublet[g_] := Block[{ham, ev, vec, ord, sym, anti},
   ham = p2Op[g]/(2 mass) + DiagonalMatrix[Vdw[g["x"]]];
   {ev, vec} = Eigensystem[ham]; ord = Ordering[Re@ev];
   sym = Normalize[vec[[ord[[1]]]]]; anti = Normalize[vec[[ord[[2]]]]];
   <|"H" -> ham, "E" -> Re@ev[[ord[[1 ;; 2]]]], "sym" -> sym, "anti" -> anti,
     "left" -> With[{plus = Normalize[sym + anti], minus = Normalize[sym - anti]},
        If[leftMask[g] . Abs[plus]^2 >= 0.5, plus, minus]]|>];
```

Apply it on a fine grid and read the doublet, the barrier height, the tunneling period $T_{\rm tun}=2\pi\hbar/(E_1-E_0)$, and a check that the Hamiltonian is Hermitian:

```wl
gD = grid[256, 14.]; dw = wellDoublet[gD];
ψSym = dw["sym"]; ψAnti = dw["anti"]; Ttun = 2 Pi hbar/(dw["E"][[2]] - dw["E"][[1]]);
leftP[ψ_] := leftMask[gD] . Abs[ψ]^2;
dw["E"]
```

The two lowest energies form a close doublet. They sit far below the barrier height $ba^4$, and their splitting $E_1-E_0$ fixes the tunneling period $T_{\rm tun}=2\pi\hbar/(E_1-E_0)$; read the barrier and the period:

```wl
{bb aa^4, Ttun}
```

A validity check: the Hamiltonian should be Hermitian to machine precision, so its spectrum is real and the `Re` above only strips rounding:

```wl
Max@Abs[dw["H"] - ConjugateTranspose[dw["H"]]]
```

The two lowest states form a doublet split by $\Delta=E_1-E_0$, far below the barrier height, so a particle prepared in one well tunnels to the other with period $T_{\rm tun}=2\pi\hbar/\Delta$. The Hamiltonian is Hermitian to machine precision, so its spectrum is real and the `Re` above only strips rounding. Now confirm the claim that forecloses the moment filter: under measurement the conditional state is genuinely non-Gaussian. Track the skewness and excess kurtosis of $|\psi(x)|^2$ along a measured trajectory (both are zero for a Gaussian):

```wl
central[ψ_, k_] := With[{d = Abs[ψ]^2, mu = xMean[gD, ψ]}, (gD["x"] - mu)^k . d];
skew[ψ_] := central[ψ, 3]/xVar[gD, ψ]^(3/2); exkurt[ψ_] := central[ψ, 4]/xVar[gD, ψ]^2 - 3;
trajDW = BlockRandom[SeedRandom[3]; trajectory[gD, Vdw, 2., 0.002, 1500][gaussian[gD, -aa, 0., 0.45]]];
Max@Abs[skew /@ trajDW[[1 ;; ;; 100]]]
```

That largest skewness magnitude runs to order one, far from the Gaussian value of zero. The excess kurtosis does the same:

```wl
Max@Abs[exkurt /@ trajDW[[1 ;; ;; 100]]]
```

Both run to order one, far from the Gaussian value of zero: a five-number filter that forces skewness and kurtosis to zero cannot represent this state, and only the grid survives. With that established, watch the two physical effects of watching a particle in a double well. First, measurement-induced localization. Prepare the symmetric eigenstate, spread equally over both wells, and monitor it repeatedly:

The symmetric eigenstate starts balanced over the two wells; the half-weighted centre point makes that exactly $P_{\rm left}=\tfrac12$:

```wl
leftP[ψSym]
```

Now monitor it across eight measured runs at a moderate strength $\lambda=2$, chosen so the two wells coexist for a visible stretch before one wins, over a run kept short enough (to $t=1.6$) that the backaction has not yet heated the particle over the barrier. Each run should localize toward one well, $P_{\rm left}$ near $0$ or $1$:

```wl
locRuns = BlockRandom[SeedRandom[10]; Table[Last@trajectory[gD, Vdw, 2., 0.002, 800][ψSym], {8}]];
leftP /@ locRuns
```

The runs localize, each toward one well or the other, and every conditional state stays a normalized pure state (only the ensemble average, as in Part V, is mixed). Different records send the particle to different wells: the measurement resolves a delocalized state into a definite, if random, location. Two runs side by side show the density localizing, each to whichever well its record selects:

```wl
bandD = Flatten@Position[UnitStep[3.5 - Abs[gD["x"]]], 1];
showRun[seed_] := ListPlot3D[(Abs[#]^2 & /@ BlockRandom[SeedRandom[seed];
      trajectory[gD, Vdw, 2., 0.002, 800][ψSym]])[[All, bandD]], PlotRange -> All,
   DataRange -> {gD["x"][[{First@bandD, Last@bandD}]], {0, 1.6}},
   ColorFunction -> (ColorData["DarkRainbow"][#3] &), Mesh -> None, AxesLabel -> {"x", "t", "|ψ|²"},
   PlotLabel -> "run " <> ToString[seed]];
GraphicsRow[{showRun[4], showRun[3]}, ImageSize -> 500]
```

Second, the quantum Zeno effect: watching the particle suppresses the coherent tunneling between wells. The cleanest way to see it is the unconditional (Lindblad) population, which is deterministic and needs no averaging. Recall from Part III that the dissipator dephases position; here it dephases the left-right coherence at rate $2\lambda a^2$, converting coherent oscillation into slow incoherent equilibration. On a coarser grid the dense master equation can afford, take the one-well state from the doublet helper and propagate it at a few measurement strengths. Set up the grid, the state, and the two observables:

```wl
gDL = grid[128, 12.]; dtDL = 0.005;
ψLwell = wellDoublet[gDL]["left"];
leftPL[ρ_] := Re[Diagonal[ρ] . leftMask[gDL]];
EmeasL[ρ_] := Re@Tr[ρ . p2Op[gDL]]/(2 mass) + Re[Diagonal[ρ] . Vdw[gDL["x"]]];
```

Propagate over one tunneling period at $\lambda=0$ (no measurement), a moderate $\lambda=0.5$, and a strong $\lambda=4$, and tabulate the left-well population against time:

```wl
smp = Round[Ttun/dtDL/8];
zenoStates[λ_] := NestList[Nest[lindblad[gDL, Vdw, λ, dtDL], #, smp] &, projector[ψLwell], 8];
zTimes = Range[0, 8] smp dtDL; zStates = zenoStates /@ {0., 0.5, 4.};
TableForm[Transpose[Prepend[Map[leftPL, zStates, {2}], zTimes]],
 TableHeadings -> {None, {"time t", "P_left, λ=0", "λ=0.5", "λ=4"}}]
```

Read the columns. Without measurement the left-well population oscillates coherently, swinging from almost fully in the left well to almost fully in the right at the half period, where the particle has tunneled across, and back. At $\lambda=0.5$ that coherent oscillation is gone: the population relaxes monotonically toward equipartition at $0.5$ and never swings across and returns. The coherent tunneling transfer is suppressed, which is the quantum Zeno signature of a continuous position measurement: watching destroys the phase coherence the tunneling needs. But there is a second effect competing with the first, and it is the same backaction $\lambda\hbar^2$ heating from Part V. Read the energy injected over one period, for $\lambda=0,0.5,4$:

```wl
EmeasL[Last@#] & /@ zStates
```

Compare those against the barrier height $ba^4$:

```wl
bb aa^4
```

By one period the energy has already climbed past the barrier at $\lambda=0.5$, and far higher still at $\lambda=4$: the stronger the measurement, the faster the population equilibrates over the top rather than through the barrier. So watching does two things at once, both growing with $\lambda$: it dephases the wells, killing the coherent oscillation, and it heats the particle, eventually over the barrier. Plot the three populations to see both at work:

```wl
ListLinePlot[Transpose[{zTimes, leftPL /@ #}] & /@ zStates,
 PlotLegends -> {"λ=0 (unwatched)", "λ=0.5 (dephased)", "λ=4 (heated over barrier)"}, Frame -> True, ImageSize -> 500,
 AspectRatio -> 1/2, PlotRange -> {0, 1}, PlotMarkers -> Automatic,
 FrameLabel -> {"time t", "left-well population"},
 PlotLabel -> "Watching kills the coherent oscillation (λ=0.5) and, harder still, heats over the barrier (λ=4)"]
```

Holding the particle in one well would take the dephasing without the heating, and unit-efficiency position measurement cannot separate them: the one coupling delivers both. This is the honest reach of the Zeno effect under continuous position measurement, and the reason a real experiment steers with feedback rather than relying on the measurement alone.

## Part VII: The Division of Labor: When to Filter, When to Grid

Let's recall what we have built and what it tells us to reach for. Two ways to simulate a continuously monitored particle appeared, and the choice between them is set by the Hamiltonian. When the conditional state lives in a small closed parametrization, the built-in stochastic integrator suffices and the geometry keeps the state physical on its own: a qubit is three Bloch numbers in a ball, and a linear (Gaussian) system, the free particle or the harmonic trap, is five moments whose covariance flow has an attracting pure-state fixed point. Both close because the Hamiltonian is at most quadratic, and for both the Kalman-Bucy moment filter is exact, the very algorithm the optomechanics experiments run.

Otherwise, hand-write the step whose structure you can exploit. Any potential beyond quadratic breaks the moment closure, the conditional state is no longer Gaussian, and only the grid split-operator scheme applies: the measurement as an exact Bayesian (Kraus) update with its renormalization, independent of $\hat H$, and the Hamiltonian as a Strang split whose kinetic factor is spectral and whose potential factor is a pointwise phase. The double well is the smallest example that forces this choice, and it repays it with physics no moment filter can represent: the non-Gaussian conditional state, the collapse into one well, and the Zeno suppression of coherent tunneling set against the backaction heating that drives the particle over the barrier.

### Where This Leaves Us

We started from one equation, the stochastic Schrödinger equation, rewrote its measurement half as a Bayesian update, and built it into a split-operator engine in which a single function carries the potential; then we ran that engine across the free particle (recovered exactly at $V=0$, heating at $\lambda\hbar^2$), the harmonic trap (a pure squeezed conditional state below the vacuum, an unconditional state that heats without bound), and the double well (a non-Gaussian conditional state that localizes into one well and has its tunneling frozen by observation while the same backaction heats it over the barrier). The measurement was the same in every case; only the unitary changed. The natural next steps are to lower the detection efficiency below unity, so a single record no longer purifies the conditional state and the trajectory becomes a mixed conditional density matrix, or to add mechanical damping, so the trapped unconditional state finally reaches the thermal steady state the bare measurement never does. Both are small edits to the engine you now have, and the code is yours to make them.
