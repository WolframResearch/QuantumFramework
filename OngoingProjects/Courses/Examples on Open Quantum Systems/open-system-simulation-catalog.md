---
Template: Default
Title: Watching Quantum Things: A Build-It-Yourself Catalog of Open-System Simulations
Author: Mads Bahrami
---

# Watching Quantum Things: A Build-It-Yourself Catalog of Open-System Simulations

<!-- #| style: Subtitle -->
A single, self-contained tour of twenty small open quantum systems, from master equations for the unconditioned average to monitored trajectories and feedback, each built from plain matrices, computed on the spot, and shown as a picture.

<!-- #| style: Author -->
Mads Bahrami

<!-- #| style: Affiliation -->
Wolfram Research, Inc.

### Setting the Stage: How This Catalog Flows

This is a hands-on tour of open quantum systems: systems that leak energy or information to their environment, watched by a detector on the output. We build a driven atom relaxing to its steady state, light in a leaky cavity, a superposition that decoheres almost the instant it forms, an atom whose fluorescence we read first as single photons and then as a homodyne current, a mechanical resonator cooled by measurement and feedback, and a dozen more. Each becomes a few lines of code you can run and change, and each ends in a picture.

I strongly believe in a computation-first way of learning: in a sense, if I cannot compute it, I cannot claim to understand it. Real understanding shows up when you can actually do the thing, simulate it, predict it, watch it, not when you can only repeat words about it. So the rhythm here is always the same: say an idea in plain words, say it again more plainly, compute it, look at what came out, and let that raise the next question.

Two notes on style. First, the file stands on its own: no other document, textbook, or note is needed. Second, I use the standard names (density matrix, Bloch vector, master equation, dissipator, coherent state) and explain each in one plain sentence the first time it appears; the code uses the same standard names, so it reads like the physics.

The thing you are reading is a live notebook. Evaluate the cells from top to bottom, because later cells lean on earlier ones. My suggestion is to look at each output and what it means before picking apart the code that made it. And you are not locked into anything: change a number, watch a different system, reseed the randomness, and see what moves. That is the whole point.

### What You Need, and How to Read This

You need three ideas. First, the state is a **density matrix** $\rho$: an $n\times n$ Hermitian matrix ($n$ = number of levels, 2 for a qubit, more for an oscillator), with the level populations on its diagonal (summing to one) and the coherences off-diagonal. Second, the system has a Hamiltonian $\hat H$ that drives its unitary evolution when left alone. Third, one random quantity recurs below: the Wiener increment $dW$, a Gaussian random step with variance equal to the time step $dt$. Where the pace steepens I will say so.

We set $\hbar = 1$ throughout, so every rate and frequency is in units of inverse time.

Let's start!

## The Toolkit: Everything the Twenty Examples Share

Before any specific system, we build the small set of tools every example reuses. There are about ten of them, each a few lines. Once they are in hand, each of the twenty systems is just a short setup, a run, a picture, and a check.

### The density matrix, and quantities read from it

The smallest system is a qubit. Everything is built from the three Pauli matrices $X$, $Y$, $Z$ and the $2\times2$ identity:

```wl
{id2, X, Y, Z} = Table[PauliMatrix[j], {j, 0, 3}];
```

The basis states are the excited state $|e\rangle = \{1,0\}$ and the ground state $|g\rangle = \{0,1\}$; their equal superposition is $|+\rangle = (|e\rangle + |g\rangle)/\sqrt2$:

```wl
excited = {1, 0}; ground = {0, 1}; plus = {1, 1}/Sqrt[2];
```

The lowering operator $\hat\sigma_- = |g\rangle\langle e|$ takes the excited state to the ground state (an atom emitting a photon). Build it from the Pauli matrices:

```wl
lower = (X - I Y)/2;
```

Check that it sends $|e\rangle$ to $|g\rangle$:

```wl
lower . excited
```

It lands on the ground state. Its conjugate transpose is the raising operator $\hat\sigma_+$, written inline when needed.

A state vector describes only a pure state. The general state, pure or mixed, is the density matrix; for a pure state it is the outer product $\rho = |\psi\rangle\langle\psi|$:

```wl
densityMatrix[v_] := KroneckerProduct[v, Conjugate[v]];
```

Form the density matrix of the excited state:

```wl
densityMatrix[excited]
```

For $|e\rangle$ this is a single 1 in the top corner.

Two quantities summarize any qubit state. The first is the **Bloch vector** $(\langle X\rangle, \langle Y\rangle, \langle Z\rangle)$: it has length 1 for a pure state and shrinks as the state mixes. Define it:

```wl
blochVector[rho_] := Re[Tr[rho . #] & /@ {X, Y, Z}];
```

Check that $|+\rangle$ points along $X$:

```wl
blochVector[densityMatrix[plus]]
```

It points along $X$, length one.

The second quantity is the **purity** $\mathrm{Tr}(\rho^2)$: 1 for a pure state, $1/2$ for the maximally mixed qubit. Define it:

```wl
purity[rho_] := Re@Tr[rho . rho];
```

Read it for a pure state and for the maximally mixed one:

```wl
{purity[densityMatrix[plus]], purity[id2/2]}
```

The pure state gives one, the maximally mixed state one half.

Last, the **expectation value** $\langle A\rangle = \mathrm{Tr}(A\rho)$ of any operator, which predicts what a meter reads:

```wl
expectation[op_, rho_] := Re@Tr[op . rho];
```

### The two terms of the master equation

A state changes for two reasons. The first is unitary evolution under the Hamiltonian $\hat H$: the commutator term $-i[\hat H,\rho]$. Define it:

```wl
commutatorTerm[H_, rho_] := -I (H . rho - rho . H);
```

The second is a leak to the environment, one Lindblad (jump) operator $\hat c$ per channel (for a decaying atom, $\hat\sigma_-$). Its **dissipator** transfers population along the channel and damps coherences, $\mathcal{D}[\hat c]\rho = \hat c\rho\hat c^\dagger - \tfrac12(\hat c^\dagger\hat c\rho + \rho\hat c^\dagger\hat c)$:

```wl
dissipator[c_, rho_] := c . rho . ConjugateTranspose[c] -
   (ConjugateTranspose[c] . c . rho + rho . ConjugateTranspose[c] . c)/2;
```

A dissipator is trace-preserving, $\mathrm{Tr}(\mathcal{D}[\hat c]\rho) = 0$, so the total probability stays one. Confirm it for a general $\hat c$ and $\rho$:

```wl
FullSimplify[Tr[dissipator[{{c11, c12}, {c21, c22}}, {{r11, r12}, {r21, r22}}]]]
```

Zero with no assumptions.

Together the two terms give the **Lindblad master equation** $\dot\rho = \mathcal{L}\rho = -i[\hat H,\rho] + \sum_k\mathcal{D}[\hat c_k]\rho$, the evolution of an open system averaged over all measurement outcomes:

```wl
lindbladian[H_, leaks_List, rho_] := commutatorTerm[H, rho] + Total[dissipator[#, rho] & /@ leaks];
```

This is a system that leaks but is not watched: smooth and deterministic, because averaging over every detector outcome removes the randomness. The first several examples are this equation with different $\hat H$ and leaks.

### The master equation, solved two independent ways

We solve the master equation two independent ways and check they agree. The first vectorizes $\rho$ (stack its rows) and applies the **Liouvillian**, the matrix whose action is $\mathcal{L}$. In this row-stacking convention (matching WL's `Flatten`, with $\overline{\phantom{x}}$ the entrywise conjugate) it is the Kronecker form
$$\mathcal{L} = -i\big(\hat H\otimes\mathbb{1} - \mathbb{1}\otimes\hat H^{T}\big) + \sum_k\Big[\hat c_k\otimes\overline{\hat c_k} - \tfrac12\big(\hat c_k^\dagger\hat c_k\otimes\mathbb{1} + \mathbb{1}\otimes(\hat c_k^\dagger\hat c_k)^{T}\big)\Big],$$
which the code builds equivalently by feeding $\mathcal{L}$ each basis matrix and stacking the results:

```wl
liouvillian[H_, leaks_List, d_] := Module[{e}, Transpose@Table[
   Flatten@lindbladian[H, leaks, ArrayReshape[UnitVector[d^2, e], {d, d}]], {e, d^2}]];
```

The state at time $t$ is then a matrix exponential of the Liouvillian applied to the vectorized initial state, reshaped back to a matrix, $\mathrm{vec}\,\rho(t) = e^{\mathcal{L}t}\,\mathrm{vec}\,\rho_0$:

```wl
evolve[H_, leaks_List, rho0_, t_] := With[{d = Length[rho0]},
   ArrayReshape[MatrixExp[liouvillian[H, leaks, d] N[t]] . Flatten[N@rho0], {d, d}]];
```

The Liouvillian also gives the **steady state** $\rho_{\mathrm{ss}}$, satisfying $\mathcal{L}\rho_{\mathrm{ss}} = 0$: the null space of the Liouvillian (there can be more than one), normalized to unit trace:

```wl
steadyState[H_, leaks_List] := With[{ns = NullSpace[liouvillian[H, leaks, Length[H]]]},
   #/Tr[#] & /@ (ArrayReshape[#, {Length[H], Length[H]}] & /@ ns)];
```

The second way integrates the master equation directly with an ODE solver, sharing no code with the Liouvillian route:

```wl
evolveODE[H_, leaks_List, rho0_, t1_] :=
  NDSolveValue[{s'[t] == lindbladian[H, leaks, s[t]], s[0] == N@rho0}, s, {t, 0, t1}];
```

Now the cross-check. Fix a test system: a qubit driven by $X$, dephased through $Z$, starting in the ground state:

```wl
Hx = 1.0 X; oneLeak = {Sqrt[0.4] Z}; start = densityMatrix[ground];
```

Integrate its master equation once with the ODE route:

```wl
solved = evolveODE[Hx, oneLeak, start, 10.0];
```

Measure the largest disagreement between the two solvers over the whole span:

```wl
Max@Table[Norm[evolve[Hx, oneLeak, start, tt] - solved[tt], "Frobenius"], {tt, 0, 10, 0.5}]
```

The largest disagreement is at the ODE solver's tolerance: the matrix exponential and the ODE march give the same state. Either is a trustworthy reference.

### The positivity-preserving measurement step

Now put a detector on a leak and condition the state on its output. The naive update, nudging $\rho$ by the record, can push $\rho$ off the set of valid states and give a negative eigenvalue, a negative probability. See it break: take $|+\rangle$, measure $Z$ with no Hamiltonian, and feed one large noise kick into the naive update:

```wl
naiveStep = densityMatrix[plus] + dissipator[Z, densityMatrix[plus]] 0.1 +
   (Z . densityMatrix[plus] + densityMatrix[plus] . Z -
      Tr[(2 Z) . densityMatrix[plus]] densityMatrix[plus]) 0.5;
```

Read the smallest eigenvalue of the updated state:

```wl
Min@Eigenvalues[naiveStep]
```

A negative eigenvalue: the naive update has left the set of valid states. Shrinking $dt$ makes it rarer but never impossible.

The fix is the structure-preserving filter of [Rouchon and Ralph (2015)](https://arxiv.org/abs/1410.5345): write the update as a sum of Kraus terms $A\rho A^\dagger$, an operator times $\rho$ times its conjugate transpose, so the result is manifestly positive at any step size. With monitored operators $\{\hat c_k\}$ at detection efficiencies $\{\eta_k\}$, unmonitored leaks $\{\hat l_j\}$, record increments $dy_k$, the no-jump generator $K = i\hat H + \tfrac12\big(\sum_k \hat c_k^\dagger \hat c_k + \sum_j \hat l_j^\dagger \hat l_j\big)$, and the measurement kick $s = \sum_k \sqrt{\eta_k}\,dy_k\,\hat c_k$, the corrected Kraus operator is
$$M = \mathbb{1} - K\,dt + s + \tfrac12 s^2 - \tfrac{dt}{2}\sum_k \eta_k\,\hat c_k^2,$$
and the stepped state is the renormalized sum
$$\tilde\rho = M\rho M^\dagger + dt\sum_j \hat l_j\,\rho\,\hat l_j^\dagger + dt\sum_k (1-\eta_k)\,\hat c_k\,\rho\,\hat c_k^\dagger, \qquad \rho' = \frac{\tilde\rho}{\mathrm{Tr}\,\tilde\rho}.$$
The three terms of $\tilde\rho$ are the measured record, the unmonitored leaks, and the unrecorded fraction $(1-\eta_k)$ of each watched channel; each has Kraus form, so $\rho$ stays valid for any $dt$, and averaging over the $dy_k$ recovers the Lindblad master equation. Build it: it takes the fixed model pieces and returns a function mapping a state and a record to the next state:

```wl
measurementStep[H_, watched_List, effs_List, unwatched_List, dt_] :=
  With[{id = IdentityMatrix[Length[H]],
    drift = I H + Total[ConjugateTranspose[#] . # & /@ Join[watched, unwatched]]/2,
    corr = Total[MapThread[#1 (#2 . #2) &, {effs, watched}]]/2},
   Function[{rho, dy},
    Module[{sig, M, top},
     sig = Total@MapThread[Sqrt[#2] #1 #3 &, {dy, effs, watched}];
     M = id - drift dt + sig + sig . sig/2 - corr dt;
     top = M . rho . ConjugateTranspose[M] +
       dt Total[# . rho . ConjugateTranspose[#] & /@ unwatched] +
       dt Total[MapThread[(1 - #2) #1 . rho . ConjugateTranspose[#1] &, {watched, effs}]];
     top/Re@Tr[top]]]];
```

Here `watched` are the monitored channels, `effs` their detection efficiencies (1 for perfect, less if the detector misses some output), and `unwatched` the undetected leaks. Revisit the kick that broke the naive step, now through the positivity-preserving step:

```wl
fixedStep = measurementStep[0 id2, {Z}, {1.}, {}, 0.1][densityMatrix[plus], {0.5}];
```

Read the smallest eigenvalue and the trace of the stepped state:

```wl
{Min@Re@Eigenvalues[fixedStep], Re@Tr[fixedStep]}
```

The smallest eigenvalue is now zero to rounding and the trace is one: the step lands on the boundary of the state space instead of leaving it. Positivity no longer depends on the step size.

### The measurement record, and one trajectory

The homodyne record at each step is the signal $\langle\hat c + \hat c^\dagger\rangle\,dt$ plus fresh Gaussian noise. Write it:

```wl
measurementRecord[rho_, watched_List, effs_List, dt_, kick_List] :=
  MapThread[Sqrt[#3] Re@Tr[(#1 + ConjugateTranspose[#1]) . rho] dt + #2 &,
   {watched, kick, effs}];
```

`kick` is the Wiener increment $dW$, drawn fresh each step. One **quantum trajectory** is: draw $dW$, form the record, advance with the positivity step. This is the only place randomness enters. It returns the times, the states, and the record:

```wl
trajectory[rho0_, H_, watched_List, effs_List, unwatched_List, dt_, tf_, seed_] :=
  BlockRandom[SeedRandom[seed];
   Module[{n = Round[tf/dt], step, kicks, states, record},
    step = measurementStep[H, watched, effs, unwatched, dt];
    kicks = RandomVariate[NormalDistribution[0, Sqrt[dt]], {n, Length[watched]}];
    states = FoldList[
      Function[{r, dw}, step[r, measurementRecord[r, watched, effs, dt, dw]]], rho0, kicks];
    record = MapThread[measurementRecord[#1, watched, effs, dt, #2] &, {Most[states], kicks}];
    <|"times" -> dt Range[0, n], "states" -> states, "record" -> record|>]];
```

### The check that ties the catalog together: trajectories average to the master equation

The central fact, used as a check almost everywhere: averaging many trajectories recovers the master-equation solution. Averaging over the record is what tracing out the detector means. Fix the step and the horizon:

```wl
dt = 0.005; tf = 8.0;
```

Run 200 trajectories of the test qubit:

```wl
crowd = Table[
   trajectory[start, Hx, oneLeak, {1.}, {}, dt, tf, k]["states"], {k, 1, 200}];
```

Integrate the master equation as the reference:

```wl
smoothRef = evolveODE[Hx, oneLeak, start, tf];
```

Measure the largest gap between the ensemble mean and that reference:

```wl
Max@MapThread[Norm[#1 - #2, "Frobenius"] &,
   {Mean[crowd], smoothRef /@ (dt Range[0, Round[tf/dt]])}]
```

The gap is the Monte-Carlo scatter of a finite ensemble, shrinking as $1/\sqrt N$. The trajectory ensemble and the master equation share no code yet agree: a numerical validation, repeated throughout, that a stochastic example is built right.

One more fact from the same run: each **conditional** trajectory stays pure, while the ensemble average is mixed. Compare the purity of one trajectory's final state with that of the ensemble mean:

```wl
{purity[Last@trajectory[start, Hx, oneLeak, {1.}, {}, dt, tf, 1]["states"]],
 purity[Last@Mean[crowd]]}
```

The single trajectory ends pure; the average ends mixed. This is the heart of the subject: under ideal monitoring (a pure initial state, every output channel watched, unit detector efficiency, no extra unrecorded noise) decoherence does not happen to any one conditioned state, which stays pure, and mixing is entirely the result of averaging over the unrecorded outcomes. Relax any of those conditions, add inefficiency, an unmonitored channel, thermal drive, or dephasing, and the conditioned state itself becomes mixed; the pure-trajectory picture is the ideal limit, and where the examples below invoke it they are in that limit.

### The Bloch-sphere plot

For a qubit the natural picture is the Bloch sphere: the Bloch vector's tip lies on the surface for a pure state, inside for a mixed one. `BlochSpherePlot` draws the labeled sphere:

```wl
bloch = ResourceFunction["BlochSpherePlot"];
```

A small octahedron will mark where a trajectory ends; build it from its six vertices:

```wl
diamond[p_, r_] := GraphicsComplex[p + # & /@ (r {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}}),
   Polygon[{{1, 3, 5}, {3, 2, 5}, {2, 4, 5}, {4, 1, 5}, {3, 1, 6}, {2, 3, 6}, {4, 2, 6}, {1, 4, 6}}]];
```

The plotter uses one visual grammar for every trajectory: the path is a curve, a small ball marks where it starts, a diamond marks where it ends, and all three share one color per trajectory. A label names the picture, a legend names the trajectories when they mean different things, and extra marks (fixed points, targets) ride along as plain graphics primitives:

```wl
blochPlot[paths_, lbl_ : None, names_ : None, extras_ : {}] := With[{g = Show[bloch[],
     Graphics3D[MapIndexed[{ColorData[97, First@#2], Thick, Line[#1],
         Sphere[First[#1], 0.035], diamond[Last[#1], 0.06]} &, paths]],
     Graphics3D[extras], PlotLabel -> lbl, ImageSize -> 340]},
   If[names === None, g, Legended[g, LineLegend[ColorData[97, #] & /@ Range[Length[names]], names]]]];
```

Run one monitored trajectory of the test qubit:

```wl
demoRun = trajectory[{{1/2, 1/2}, {1/2, 1/2}}, Hx, oneLeak, {1.}, {}, 0.01, 3.0, 4];
```

Now visualize the monitored trajectory on the Bloch sphere, from its starting ball to its final diamond:

```wl
blochPlot[{blochVector /@ demoRun["states"]}, "one monitored trajectory, ball to diamond"]
```

The Bloch vector wanders as the record comes in; one seed is one run. Change the seed and it wanders differently; averaging many recovers the smooth relaxation checked above.

The sphere shows the shape of the path; the same run also reads as three numbers in time, one curve for each Bloch component:

```wl
ListLinePlot[Transpose[blochVector /@ demoRun["states"]],
 PlotStyle -> {ColorData[97, 1], ColorData[97, 2], ColorData[97, 3]},
 PlotLegends -> {"\[LeftAngleBracket]\[Sigma]x\[RightAngleBracket]", "\[LeftAngleBracket]\[Sigma]y\[RightAngleBracket]", "\[LeftAngleBracket]\[Sigma]z\[RightAngleBracket]"},
 Frame -> True, GridLines -> Automatic, PlotRange -> {-1.05, 1.05}, ImageSize -> 520,
 FrameLabel -> {"time", "Bloch component"}, DataRange -> MinMax@demoRun["times"],
 PlotLabel -> "the same run, component by component"]
```

All three components move. $\langle\hat\sigma_x\rangle$ starts at one and decays toward zero, since the drive turns the state about the $x$-axis and cannot rebuild what watching $\hat\sigma_z$ erases; meanwhile $\langle\hat\sigma_y\rangle$ and $\langle\hat\sigma_z\rangle$ wander as that same drive rotates them and the record tips the state toward a pole. The sphere and the curves are one trajectory seen two ways: the sphere gives the shape of the path, the curves give each direction changing in time. This component-by-component reading is how a watched qubit is usually plotted, and most later examples show just the one component the measurement makes interesting, $\langle\hat\sigma_z\rangle$ for a $\hat\sigma_z$ measurement. With the qubit toolkit in hand, on to the first systems.

### The same machinery, in standard notation

Every example below opens with a one-line statement of its problem in standard notation, so the formal problem and the plain solution sit side by side. It is only names for the machinery you just built. The state is the density matrix $\rho$; the averaged evolution is the **master equation**
$$\dot\rho = \mathcal{L}\rho = -i[\hat H,\rho] + \sum_k \mathcal{D}[\hat c_k]\rho, \qquad \hbar = 1,$$
whose Hamiltonian part $-i[\hat H,\rho]$ is `commutatorTerm` and each of whose dissipators
$$\mathcal{D}[\hat c]\rho = \hat c\,\rho\,\hat c^\dagger - \tfrac12\bigl(\hat c^\dagger\hat c\,\rho + \rho\,\hat c^\dagger\hat c\bigr)$$
is `dissipator`. For a qubit the **Bloch vector** is $(x,y,z) = (\langle\hat\sigma_x\rangle,\langle\hat\sigma_y\rangle,\langle\hat\sigma_z\rangle)$, with $\rho = \tfrac12(\mathbb{1} + x\hat\sigma_x + y\hat\sigma_y + z\hat\sigma_z)$; the lowering operator is $\hat\sigma_- = |g\rangle\langle e|$, and for an oscillator the annihilation operator is $\hat a$ with $[\hat a,\hat a^\dagger] = 1$. When the record is kept, `measurementStep` integrates a **stochastic master equation** whose measurement (diffusion) term for a monitored $\hat c$ is $\mathcal{H}[\hat c]\rho = \hat c\rho + \rho\hat c^\dagger - \langle\hat c + \hat c^\dagger\rangle\rho$, driven by $dW$ with $dW^2 = dt$; photon counting instead applies the jump $\mathcal{G}[\hat c]\rho = \hat c\rho\hat c^\dagger/\mathrm{Tr}[\hat c\rho\hat c^\dagger] - \rho$ at click times.

### The Twenty Governing Equations, at a Glance

Here is the dynamics of every example in one place. All stochastic equations are Itô equations. For a Hermitian observable $\hat M$ monitored through $\hat c=\sqrt{k}\,\hat M$, the catalog uses
$$
k\mathcal D[\hat M]\rho=-\frac{k}{2}[\hat M,[\hat M,\rho]],\qquad
dJ=2\sqrt{k}\,\langle\hat M\rangle_c\,dt+dW,
$$
so $k$ is the coefficient carried by the monitored operator $\sqrt{k}\hat M$ in the code. Equivalently, the rescaled record $dy=dJ/(2\sqrt{k})$ obeys $dy=\langle\hat M\rangle_c\,dt+dW/(2\sqrt{k})$.

1. **Pure dephasing.**
   $$
   \dot\rho=\gamma\mathcal D[\hat\sigma_z]\rho,
   \qquad
   (\dot x,\dot y,\dot z)=(-2\gamma x,-2\gamma y,0).
   $$

2. **Driven, damped atom.**
   $$
   \dot\rho=-i\left[\frac{\Omega}{2}\hat\sigma_x,\rho\right]
   +\gamma\mathcal D[\hat\sigma_-]\rho.
   $$

3. **Leaky cavity cat.**
   $$
   \dot\rho=\gamma\mathcal D[\hat a]\rho,
   \qquad
   \rho(0)=\frac{(|\alpha\rangle+|-\alpha\rangle)(\langle\alpha|+\langle-\alpha|)}
   {2\bigl(1+e^{-2|\alpha|^2}\bigr)}.
   $$

4. **Atom-made cat and its decoherence.** During the dispersive atom-cavity interaction and the subsequent storage interval, respectively,
   $$
   \dot\rho_{AC}=-i[\chi\,\hat a^\dagger\hat a\,\hat\sigma_z,\rho_{AC}],
   \qquad
   \dot\rho_C=\gamma\mathcal D[\hat a]\rho_C.
   $$

5. **Dispersive qubit readout.**
   $$
   \dot\rho=-i\bigl[\epsilon(\hat a+\hat a^\dagger)
   +\chi\,\hat a^\dagger\hat a\,\hat\sigma_z,\rho\bigr]
   +\gamma\mathcal D[\hat a]\rho,
   $$
   and, conditional on $\sigma_z=s=\pm1$, the coherent pointer obeys
   $$
   \dot\alpha_s=-i\epsilon-\left(is\chi+\frac{\gamma}{2}\right)\alpha_s.
   $$

6. **Quantum Brownian motion.** In the high-temperature Caldeira-Leggett limit,
   $$
   \dot\rho=-i\left[\frac{\hat P^2}{2M}+\frac12M\Omega^2\hat X^2,\rho\right]
   -i\gamma[\hat X,\{\hat P,\rho\}]
   -2\gamma Mk_BT[\hat X,[\hat X,\rho]].
   $$
   For the Gaussian simulation this closes as
   $$
   \dot{\bar r}=A\bar r,\qquad
   \dot\Sigma=A\Sigma+\Sigma A^T+D,\qquad
   A=\begin{pmatrix}0&1/M\\-M\Omega^2&-2\gamma\end{pmatrix},\quad
   D=\begin{pmatrix}0&0\\0&4\gamma Mk_BT\end{pmatrix}.
   $$
   The completely positive warm-oscillator comparison uses
   $$
   \dot\rho=\gamma(n_T+1)\mathcal D[\hat a]\rho
   +\gamma n_T\mathcal D[\hat a^\dagger]\rho.
   $$

7. **Photon counting.** With $\lambda_t=\operatorname{Tr}(\hat c^\dagger\hat c\rho_c)$ and $\mathbb E[dN_t\mid\rho_c]=\lambda_t\,dt$,
   $$
   d\rho_c=\bigl(-i[\hat H,\rho_c]+\mathcal D[\hat c]\rho_c\bigr)\,dt
   +\mathcal G[\hat c]\rho_c\,(dN_t-\lambda_t\,dt),
   $$
   where $\mathcal G[\hat c]\rho_c=\hat c\rho_c\hat c^\dagger/\lambda_t-\rho_c$.

8. **Homodyne detection.**
   $$
   d\rho_c=\bigl(-i[\hat H,\rho_c]+\mathcal D[\hat c]\rho_c\bigr)\,dt
   +\mathcal H[\hat c]\rho_c\,dW,
   \qquad
   dJ=\langle\hat c+\hat c^\dagger\rangle_c\,dt+dW.
   $$

9. **Heterodyne detection.** Put $\hat c_I=\hat c/\sqrt2$ and $\hat c_Q=i\hat c/\sqrt2$. Then
   $$
   d\rho_c=\bigl(-i[\hat H,\rho_c]+\mathcal D[\hat c]\rho_c\bigr)\,dt
   +\mathcal H[\hat c_I]\rho_c\,dW_I
   +\mathcal H[\hat c_Q]\rho_c\,dW_Q,
   $$
   with $dW_I^2=dW_Q^2=dt$ and $dW_IdW_Q=0$; the two records are
   $$
   dJ_j=\langle\hat c_j+\hat c_j^\dagger\rangle_c\,dt+dW_j,
   \qquad j\in\{I,Q\}.
   $$

10. **Quantum Zeno dynamics.** For the code's monitored operator $\sqrt{k}\,\hat\sigma_z$,
    $$
    d\rho_c=-i\frac{\Omega}{2}[\hat\sigma_x,\rho_c]\,dt
    +k\mathcal D[\hat\sigma_z]\rho_c\,dt
    +\sqrt{k}\,\mathcal H[\hat\sigma_z]\rho_c\,dW.
    $$
    Equivalently, the deterministic measurement term is $-\tfrac{k}{2}[\hat\sigma_z,[\hat\sigma_z,\rho_c]]\,dt$. This convention gives the Liouvillian pair $-k\pm i\sqrt{\Omega^2-k^2}$ and the exceptional point $k=\Omega$ used below.

11. **Charge qubit monitored by a QPC.** Since $\hat n_1=(\mathbb{1}-\hat\sigma_z)/2$,
    $$
    d\rho_c=-i\left[\frac{\Omega}{2}\hat\sigma_x,\rho_c\right]\,dt
    +\frac{\kappa}{4}\mathcal D[\hat\sigma_z]\rho_c\,dt
    +\frac{\sqrt{\kappa}}{2}\mathcal H[\hat\sigma_z]\rho_c\,dW,
    \qquad
    dJ=\sqrt{\kappa}\,\langle\hat\sigma_z\rangle_c\,dt+dW.
    $$
    Dropping the innovation term gives the unconditional equation $\dot\rho=-i[(\Omega/2)\hat\sigma_x,\rho]+\kappa\mathcal D[\hat n_1]\rho$.

12. **Measurement-induced localization.**
    $$
    d\rho_c=k\mathcal D[\hat\sigma_z]\rho_c\,dt
    +\sqrt{k}\,\mathcal H[\hat\sigma_z]\rho_c\,dW,
    \qquad
    dy=\langle\hat\sigma_z\rangle_c\,dt+\frac{dW}{2\sqrt{k}}.
    $$

13. **Quantum Kalman filter.** With $\hat r=(\hat x,\hat p)^T$,
    $$
    \dot\Sigma=A\Sigma+\Sigma A^T+D-\Sigma C^TC\Sigma,
    \qquad
    d\bar r_c=A\bar r_c\,dt+\Sigma C^T(dY-C\bar r_c\,dt),
    $$
    $$
    A=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\qquad
    D=\begin{pmatrix}0&0\\0&k\end{pmatrix},\qquad
    C=\begin{pmatrix}2\sqrt{k}&0\end{pmatrix},\qquad
    dY=C\bar r_c\,dt+dW.
    $$

14. **Markovian measurement feedback.** For the feedback Hamiltonian $\hat H_{\rm fb}=J(t)\hat F$,
    $$
    \dot\rho=-i[\hat H,\rho]+\mathcal D[\hat c]\rho
    -i[\hat F,\hat c\rho+\rho\hat c^\dagger]
    +\frac1\eta\mathcal D[\hat F]\rho.
    $$
    Equivalently,
    $$
    \dot\rho=-i\left[\hat H+\frac12(\hat c^\dagger\hat F+\hat F\hat c),\rho\right]
    +\mathcal D[\hat c-i\hat F]\rho+\frac{1-\eta}{\eta}\mathcal D[\hat F]\rho.
    $$
    The simulated ideal loop has $\hat c=\sqrt G\,\hat\sigma_z$, $\hat F=\sqrt G\,\hat\sigma_y$, and therefore $\dot\rho=G\mathcal D[\hat\sigma_z-i\hat\sigma_y]\rho$.

15. **Thermalization.**
    $$
    \dot\rho=-i[\omega\hat a^\dagger\hat a,\rho]
    +\gamma(n_T+1)\mathcal D[\hat a]\rho
    +\gamma n_T\mathcal D[\hat a^\dagger]\rho,
    \qquad
    n_T=(e^{\beta\omega}-1)^{-1}.
    $$

16. **Rapid purification by adaptive measurement.** For $\hat M=\sin\theta\,\hat\sigma_x+\cos\theta\,\hat\sigma_z$,
    $$
    d\rho_c=k\mathcal D[\hat M]\rho_c\,dt+\sqrt{k}\,\mathcal H[\hat M]\rho_c\,dW,
    $$
    and the ensemble-mean impurity $L=(1-|\vec a|^2)/2$ obeys
    $$
    \mathbb E[dL]=-4kL\bigl(\sin^2\theta+2L\cos^2\theta\bigr)dt.
    $$
    Feedback holds $\theta=\pi/2$, giving the deterministic equation $\dot L=-4kL$.

17. **Feedback cooling.** In the code convention, with monitored operator $\sqrt{k}\,\hat x$,
    $$
    \begin{aligned}
    d\rho_c={}&-i[\omega\hat b^\dagger\hat b-f(t)\hat x,\rho_c]\,dt
    +k\mathcal D[\hat x]\rho_c\,dt+\sqrt{\eta k}\,\mathcal H[\hat x]\rho_c\,dW\\
    &+\gamma(n_T+1)\mathcal D[\hat b]\rho_c\,dt
    +\gamma n_T\mathcal D[\hat b^\dagger]\rho_c\,dt,
    \qquad f(t)=-G\langle\hat p\rangle_c.
    \end{aligned}
    $$

18. **Linear and nonlinear trajectory bookkeepings.** The unnormalized linear state obeys
    $$
    d|\tilde\psi\rangle=\left(-i\hat H-\frac12\sum_j\hat R_j^\dagger\hat R_j\right)|\tilde\psi\rangle dt
    +\sum_j\hat R_j|\tilde\psi\rangle dW_j.
    $$
    The normalized conditional state obeys
    $$
    d\rho_c=\mathcal L\rho_c\,dt+\sum_j\mathcal H[\hat R_j]\rho_c\,d\widehat W_j.
    $$

19. **Driven atom, averaged and watched.** The general unconditioned equation is
    $$
    \dot\rho=-i\left[\frac{\Delta\omega}{2}\hat\sigma_z+\frac{\Omega}{2}\hat\sigma_x,\rho\right]
    +\gamma(n_T+1)\mathcal D[\hat\sigma_-]\rho
    +\gamma n_T\mathcal D[\hat\sigma_+]\rho
    +\gamma k_d\mathcal D[\hat\sigma_z]\rho.
    $$
    In the resonant, zero-temperature, no-extra-dephasing corner watched by ideal homodyne detection,
    $$
    d\rho_c=\left(-i\left[\frac{\Omega}{2}\hat\sigma_x,\rho_c\right]
    +\gamma\mathcal D[\hat\sigma_-]\rho_c\right)dt
    +\sqrt\gamma\,\mathcal H[\hat\sigma_-]\rho_c\,dW.
    $$

20. **Mollow spectrum.** The two-time correlations are generated by
    $$
    \dot\rho=-i\left[\frac{\Omega}{2}\hat\sigma_x,\rho\right]
    +\gamma\mathcal D[\hat\sigma_-]\rho,
    $$
    and quantum regression gives
    $$
    S_{\rm inel}(\mu)\propto\operatorname{Re}\int_0^\infty
    e^{i\mu t}\,\langle\delta\hat\sigma_+(t)\delta\hat\sigma_-(0)\rangle_{\rm ss}\,dt.
    $$

## Part One: Systems We Can Solve Symbolically

The first two systems have exact closed-form answers we can check against the code. Both are qubits, so first a convenience: read the master equation as three **Bloch equations** for $(\dot x, \dot y, \dot z)$. Write the general Bloch state $\rho = \tfrac12(\mathbb{1} + x X + y Y + z Z)$:

```wl
blochState[x_, y_, z_] := (id2 + x X + y Y + z Z)/2;
```

Feed it to the Lindbladian and extract each component's rate:

```wl
rates[H_, leaks_List] := ComplexExpand@Table[
   Tr[d . lindbladian[H, leaks, blochState[x, y, z]]], {d, {X, Y, Z}}];
```

### Pure Dephasing

**The problem.** A two-level system dephased along $\hat\sigma_z$, with no Hamiltonian:
$$\dot\rho = \gamma\,\mathcal{D}[\hat\sigma_z]\rho.$$
Since $\hat\sigma_z^\dagger\hat\sigma_z = \mathbb{1}$ this is pure dephasing, so in Bloch components $\dot x = -2\gamma x$, $\dot y = -2\gamma y$, $\dot z = 0$, giving $x(t) = x_0 e^{-2\gamma t}$, $y(t) = y_0 e^{-2\gamma t}$, $z(t) = z_0$. Both poles $|\hat\sigma_z = \pm 1\rangle$ (the $\hat\sigma_z$ eigenstates $|e\rangle$ and $|g\rangle$, not to be confused with the $\hat\sigma_x$ eigenstate the code calls `plus`) are fixed points, so there is no unique steady state: any mixture $\rho_{\mathrm{ss}} = p_+\,|e\rangle\langle e| + p_-\,|g\rangle\langle g|$ sits still forever. Here is the same problem in plain words, solved from scratch and checked against the machine.

Pure dephasing continuously measures $\hat\sigma_z$: it damps the transverse coherences and leaves the populations fixed. No Hamiltonian, one dissipator. Read the Bloch equations:

```wl
Clear[\[Gamma]]; Simplify[rates[0 id2, {Sqrt[\[Gamma]] Z}], \[Gamma] > 0]
```

The transverse components $x, y$ decay at rate $2\gamma$; $z$ is fixed. Hand those rates, exactly as the Lindbladian produced them, to the differential-equation solver as the Bloch equations, and solve from any initial Bloch vector:

```wl
Clear[x, y, z, t, x0, y0, z0];
closedDephasing[t_] = {x[t], y[t], z[t]} /. First@DSolve[Join[
     Thread[{x'[t], y'[t], z'[t]} ==
       (Simplify[rates[0 id2, {Sqrt[\[Gamma]] Z}], \[Gamma] > 0] /. {x -> x[t], y -> y[t], z -> z[t]})],
     {x[0] == x0, y[0] == y0, z[0] == z0}], {x, y, z}, t]
```

The transverse components decay exponentially, $z$ is frozen. Check the closed form against the master equation over time:

```wl
With[{gg = 0.5, a0 = {0.8, 0.4, 0.3}},
 Max@Table[Norm[(closedDephasing[tt] /. {\[Gamma] -> gg, x0 -> a0[[1]], y0 -> a0[[2]], z0 -> a0[[3]]}) -
    blochVector[evolve[0 id2, {Sqrt[gg] Z}, blochState @@ a0, tt]]], {tt, 0, 3, 0.3}]]
```

They agree to rounding.

Because $z$ never changes, there is no unique steady state: any mixture of $|e\rangle$ and $|g\rangle$ is stationary. Ask the toolkit for all the steady states of the dephasing generator:

```wl
steadyState[0 id2, {Sqrt[0.5] Z}]
```

The two poles return, $|e\rangle\langle e|$ and $|g\rangle\langle g|$, each a stationary density matrix on its own. Count them:

```wl
Length@%
```

Two, so a whole line of mixtures between them: a degenerate steady state. Every *driven* system in the catalog has a unique steady state; the degeneracy returns only where the drive is off and the dynamics is again pure dephasing, as in the undriven measurement-induced localization of Part Five.

Now visualize the dephasing behavior above: start two Bloch vectors off-axis and watch each slide straight toward the $z$-axis (ball to diamond), keeping its height, while the poles (black dots) stay fixed:

```wl
With[{gg = 0.6, starts = {{0.85, 0.25, 0.4}, {-0.6, -0.55, -0.35}}},
 blochPlot[
  Table[Table[closedDephasing[tt] /. {\[Gamma] -> gg, x0 -> s[[1]], y0 -> s[[2]], z0 -> s[[3]]}, {tt, 0, 2.5, 0.05}], {s, starts}],
  "dephasing slides each state to the z-axis", {"start (0.85, 0.25, 0.4)", "start (-0.6, -0.55, -0.35)"},
  {Black, PointSize[0.03], Point[{0, 0, 1}], Point[{0, 0, -1}]}]]
```

Each tip drifts to the $z$-axis and stops at its own height: the transverse coherence is gone, the populations untouched. This is dephasing in its purest form, decoherence with nothing else attached.

### A Driven, Damped Atom Relaxing to Steady State

**The problem.** A two-level atom driven by a laser (a drive about $\hat\sigma_x$) while it decays to its ground state:
$$\mathcal{L}\rho = -i\,\frac{\Omega}{2}[\hat\sigma_x,\rho] + \gamma\,\mathcal{D}[\hat\sigma_-]\rho.$$
In Bloch components $\dot x = -\frac{\gamma}{2}x$, $\dot y = -\Omega z - \frac{\gamma}{2}y$, $\dot z = \Omega y - \gamma(z+1)$, which spiral in to the steady state
$$x_{\mathrm{ss}} = 0,\qquad y_{\mathrm{ss}} = \frac{2\Omega\gamma}{\gamma^2 + 2\Omega^2},\qquad z_{\mathrm{ss}} = -\frac{\gamma^2}{\gamma^2 + 2\Omega^2},$$
at the rate set by the eigenvalues $\lambda_\pm = -\frac{3\gamma}{4} \pm i\tilde\Omega$, with $\tilde\Omega = \sqrt{\Omega^2 - (\gamma/4)^2}$ (a damped oscillation once $\Omega > \gamma/4$). Solve it and watch the spiral.

Now both terms at once: a laser drive about $X$ and spontaneous emission through $\hat\sigma_-$. Read the Bloch equations:

```wl
Clear[\[Gamma], \[CapitalOmega]]; Simplify[rates[(\[CapitalOmega]/2) X, {Sqrt[\[Gamma]] lower}], \[Gamma] > 0]
```

The component $x$ decays; the drive couples $y$ and $z$ while emission pulls $z$ toward the ground state. This has a unique steady state. Find it by setting all three rates to zero:

```wl
Clear[x, y, z];
steady = First@Solve[Thread[rates[(\[CapitalOmega]/2) X, {Sqrt[\[Gamma]] lower}] == 0], {x, y, z}];
```

Read the steady state in closed form:

```wl
Simplify[{x, y, z} /. steady, \[Gamma] > 0]
```

The steady state in closed form: no transverse lean, and a tilt set by the ratio of drive to emission.

How is it approached? The Bloch equations are linear, so the approach is set by the eigenvalues of their drift matrix:

```wl
Simplify[Eigenvalues[D[rates[(\[CapitalOmega]/2) X, {Sqrt[\[Gamma]] lower}], {{x, y, z}}]], \[Gamma] > 0]
```

Two are complex conjugates: a **spiral**, a ringing approach, whenever the drive exceeds a quarter of the emission rate. Confirm the closed-form steady state against the toolkit's steady-state finder:

```wl
With[{gg = 1., om = 3.},
 Norm[(Simplify[{x, y, z} /. steady] /. {\[Gamma] -> gg, \[CapitalOmega] -> om}) -
   blochVector[First@steadyState[(om/2) X, {Sqrt[gg] lower}]]]]
```

They match to rounding. From the excited state the Bloch vector will spiral in (the complex eigenvalue pair guarantees the ringing).

But a special pair of initial states skips the spiral entirely and slides straight in, and the reason sits in the drift matrix itself. Check that $\hat x = (1,0,0)$ is one of its eigenvectors, with eigenvalue $-\gamma/2$:

```wl
Simplify[D[rates[(\[CapitalOmega]/2) X, {Sqrt[\[Gamma]] lower}], {{x, y, z}}] . {1, 0, 0}, \[Gamma] > 0]
```

The result is $-\tfrac{\gamma}{2}(1,0,0)$: a deviation from the steady state leaning purely along $x$ never rotates into $y$ or $z$; it only shrinks, at rate $\gamma/2$. In other words, the $x$-line through the steady state is invariant, and motion along it is a straight slide. The pure states on that line are where it pierces the unit sphere; since $x_{\mathrm{ss}} = 0$, attach $x = \pm\sqrt{1 - y_{\mathrm{ss}}^2 - z_{\mathrm{ss}}^2}$ to the steady components:

```wl
flatPair = Simplify[With[{ss = {x, y, z} /. steady},
    {# Sqrt[1 - ss[[2]]^2 - ss[[3]]^2], ss[[2]], ss[[3]]} & /@ {1, -1}], \[CapitalOmega] > 0 && \[Gamma] > 0]
```

The square root collapses to a rational form: the pair is $\bigl(\pm 2\Omega^2,\ 2\Omega\gamma,\ -\gamma^2\bigr)/(\gamma^2 + 2\Omega^2)$. Check that both have unit Bloch length, so both are pure:

```wl
Simplify[#.# & /@ flatPair, \[CapitalOmega] > 0 && \[Gamma] > 0]
```

Both give one exactly.

This two-state family is the **physically realizable (PR) pair** of the driven atom: a two-state ensemble a monitored trajectory can stay inside. The contrast is with the eigenstates of the mixed steady state (the *orthogonal* ensemble), which are not realizable here, and the difference is visible in the drift. Read the Lindblad velocity $\dot{\vec r}$ at one endpoint of each pair, the PR endpoint against the orthogonal endpoint:

```wl
driftAt[r_, om_, gg_] := Re@Table[Tr[d . lindbladian[(om/2) X, {Sqrt[gg] lower}, blochState @@ r]], {d, {X, Y, Z}}];
prPair = N[flatPair /. {\[CapitalOmega] -> 3., \[Gamma] -> 1.}];
orthPair = {#, -#} &@Normalize@N[{x, y, z} /. steady /. {\[CapitalOmega] -> 3., \[Gamma] -> 1.}];
{driftAt[prPair[[1]], 3., 1.], driftAt[orthPair[[1]], 3., 1.]}
```

The PR endpoint's drift points purely along $\hat x$, the direction of its own chord (the two PR states differ only in $x$), so the state slides along the line joining the pair and never leaves it. The orthogonal endpoint's drift points along $\hat z$, across its chord (which lies in the $y$-$z$ plane), so the drive immediately rotates that state out of its own pair: no continuous measurement can keep the system perpetually in the orthogonal eigenbasis. That is exactly why the orthogonal ensemble is not physically realizable while the PR pair is.

Now visualize all three paths, the spiral from the excited state and the two flat slides from the derived pair, with the steady state marked. This is the *unconditional* picture, the record averaged away, so every start relaxes into the mixed steady state:

```wl
With[{gg = 1., om = 3.},
 Module[{ss, starts, path},
  ss = N[{x, y, z} /. steady /. {\[Gamma] -> gg, \[CapitalOmega] -> om}];
  starts = N[flatPair /. {\[Gamma] -> gg, \[CapitalOmega] -> om}];
  path[r0_] := With[{run = evolveODE[(om/2) X, {Sqrt[gg] lower}, r0, 6.]},
    blochVector[run[#]] & /@ Range[0, 6, 0.02]];
  blochPlot[{path[densityMatrix[excited]], path[blochState @@ starts[[1]]],
    path[blochState @@ starts[[2]]]},
   "every start lands on the steady state (black dot)",
   {"from the excited state", "robust start +", "robust start -"},
   {Black, PointSize[0.035], Point[ss]}]]]
```

The excited start rings down to the steady state; the two special starts glide straight in; all three diamonds land on the same Bloch vector (the black dot), inside the ball. A driven, damped atom's *unconditional* steady state is mixed, set by the balance of drive and emission.

## Part Two: Light in a Box, and States That Live in Two Places

The next systems are a single mode of light in a cavity. Its states span the Fock ladder $|0\rangle, |1\rangle, |2\rangle, \dots$, so the density matrix is larger, one row and column per Fock state. We truncate at a cutoff high enough that the state never reaches it, and check by raising the cutoff.

The **annihilation** and **creation** operators step down and up the ladder, $\hat a\,|n\rangle = \sqrt{n}\,|n-1\rangle$ and $\hat a^\dagger|n\rangle = \sqrt{n+1}\,|n+1\rangle$, so $\hat a$ carries $\sqrt{n}$ on its superdiagonal:

```wl
annihilation[n_] := N@DiagonalMatrix[Sqrt[Range[n - 1]], 1];
creation[n_] := ConjugateTranspose[annihilation[n]];
```

The **coherent state** $|\alpha\rangle$ (the most classical state, an eigenstate of $\hat a$) is the Fock-space sum
$$|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^{\infty}\frac{\alpha^n}{\sqrt{n!}}\,|n\rangle.$$
Build it on the coefficients $\alpha^n/\sqrt{n!}$ up to the cutoff, normalized in the truncated basis (the numerical normalization stands in for the $e^{-|\alpha|^2/2}$ prefactor):

```wl
coherentState[n_, a_] := With[{v = Table[If[k == 0 && a == 0, 1, a^k]/Sqrt[k!], {k, 0, n - 1}]}, N[v/Norm[v]]];
```

The **displacement operator** $\hat D(\alpha) = e^{\alpha\hat a^\dagger - \alpha^*\hat a}$, which shifts the vacuum to $|\alpha\rangle$:

```wl
displacement[n_, a_] := MatrixExp[a creation[n] - Conjugate[a] annihilation[n]];
```

The **Wigner function** $W(x,p)$ pictures a state on the phase-space plane: positive where the state is concentrated and, for a nonclassical state, negative in places. In the canonical $(x,p)$ normalization its displaced-parity definition is
$$W(x,p) = \frac{1}{\pi}\,\mathrm{Tr}\!\left[\hat D^\dagger(\alpha)\rho\hat D(\alpha)\hat{Π}\right], \qquad \alpha = \frac{x + i p}{\sqrt{2}}.$$

A word on the coordinate convention, kept fixed for every phase-space picture below. The axes $x$ and $p$ are the dimensionless quadratures $\hat x = (\hat a + \hat a^\dagger)/\sqrt2$, $\hat p = (\hat a - \hat a^\dagger)/(i\sqrt2)$ (so $[\hat x,\hat p] = i$ and the vacuum has $\Sigma_{xx} = \Sigma_{pp} = 1/2$), and $\alpha = (x + ip)/\sqrt2$ is the complex amplitude. The prefactor $1/\pi$ is the textbook normalization on these axes, $\int W\,dx\,dp = 1$, so the vacuum is $W(x,p) = \tfrac1\pi e^{-x^2-p^2}$.

Rather than exponentiate the displacement operator at every grid point, evaluate this trace directly in the Fock basis. Define
$$z = 2\alpha = \sqrt{2}(x + i p), \qquad b = |z|^2 = 4|\alpha|^2.$$
For the $\ell$th superdiagonal of $\rho$, define the radial Laguerre sum
$$F_\ell(b) = \sum_{m=0}^{n-\ell-1}\rho_{m,m+\ell}(-1)^m\sqrt{\frac{m!\,\ell!}{(m+\ell)!}}\,L_m^{(\ell)}(b),$$
where $L_m^{(\ell)}$ is a generalized Laguerre polynomial. The complete Wigner function is then
$$W(x,p) = \frac{1}{\pi}e^{-b/2}\,\mathrm{Re}\!\left[F_0(b) + 2\sum_{\ell=1}^{n-1}\frac{z^\ell}{\sqrt{\ell!}}F_\ell(b)\right].$$
In other words, the recurrence computes the radial function $F_\ell$ for one diagonal of $\rho$; the powers $z^\ell$ supply the angular dependence, and the Gaussian $e^{-b/2}$ supplies the envelope. Build the whole evaluation as one function: the Laguerre recurrence lives inside as a local helper, the input may be a state vector or a density matrix (a vector becomes its density matrix first), the cutoff is read off the input's size, and the grid, exact for the truncated state so the frame's far corners stay clean, is returned as an interpolating function of $(x,p)$:

```wl
wignerRepresentation[psi_?VectorQ, reach_, pts_] := wignerRepresentation[densityMatrix[psi], reach, pts];
wignerRepresentation[rho_?MatrixQ, reach_, pts_] := Module[{n = Length[rho], grid, a2, b, lag, w},
   grid = N@Subdivide[-reach, reach, pts - 1];
   a2 = Sqrt[2.] Outer[#2 + I #1 &, grid, grid]; b = Abs[a2]^2;
   lag[L_, c_] := With[{m = Length[c]},
     Which[
      m == 1, c[[1]] + 0. b,
      m == 2, c[[1]] - c[[2]] Sqrt[1./(L + 1)] (L + 1. - b),
      True, With[{yy = Fold[Function[{y, cj}, With[{j = cj[[2]]},
           {cj[[1]] - y[[2]] Sqrt[((j - 1.) (L + j - 1.))/((L + j) j)],
            y[[1]] - y[[2]] (L + 2. j - 1 - b) Sqrt[1./((L + j) j)]}]],
         {c[[-2]] + 0. b, c[[-1]] + 0. b}, Transpose[{Reverse[c[[;; -3]]], N@Range[m - 1, 2, -1]}]]},
        yy[[1]] - yy[[2]] Sqrt[1./(L + 1)] (L + 1. - b)]]];
   w = Re[Fold[lag[#2 - 1, If[#2 == 1, 1, 2] Diagonal[rho, #2 - 1]] + #1 a2/Sqrt[#2] &,
       ConstantArray[2. rho[[1, -1]], {pts, pts}], Range[n - 1, 1, -1]]] Exp[-b/2]/Pi;
   ListInterpolation[Transpose[w], {{-reach, reach}, {-reach, reach}}]];
```

Check it on the vacuum, whose Wigner function is the closed form $\tfrac{1}{\pi}e^{-x^2-p^2}$, feeding the state vector directly:

```wl
With[{w = wignerRepresentation[coherentState[8, 0], 3., 25]},
 Max@Table[Abs[w[xx, pp] - (1/Pi) Exp[-xx^2 - pp^2]], {xx, -2, 2, 0.5}, {pp, -2, 2, 0.5}]]
```

They agree to rounding on the grid's own points; between points the interpolation is accurate to the grid's resolution.

Wigner pictures are compared, across times or across states, and the comparison is only honest on one shared color scale. The plotter takes a list of Wigner interpolating functions and a label per member, finds the largest $|W|$ over the whole list, and draws every panel and the bar legend with the same diverging map: a light neutral at zero so empty phase space reads as empty, warm (orange) for positive, cold (blue) for negative, with the $W=0$ contour drawn as a thin gray line so the sign change is unmistakable. An optional color map replaces the default; it receives the signed fraction $W/\max|W|$ in $[-1,1]$:

```wl
ClearAll[wignerPlot]
wignerPlot[ws : {__InterpolatingFunction}, labels_List, color_ : Automatic] :=
  With[{m = Max[Abs[#["ValuesOnGrid"]] & /@ ws], reach = Max[Abs[ws[[1]]["Domain"]]],
    cf = If[color === Automatic,
      Blend[{{-1, RGBColor[0.15, 0.4, 0.8]}, {0, GrayLevel[0.97]}, {1, RGBColor[0.9, 0.45, 0.1]}}, #] &, color]},
   With[{raw = cf[#/m] &},
    Legended[GraphicsRow[MapThread[
       Show[DensityPlot[#1[xx, pp], {xx, -reach, reach}, {pp, -reach, reach}, PlotPoints -> 55,
          PlotRange -> All, ColorFunction -> raw, ColorFunctionScaling -> False,
          FrameLabel -> {"x", "p"}, PlotLabel -> #2],
         ContourPlot[#1[xx, pp] == 0, {xx, -reach, reach}, {pp, -reach, reach},
          ContourStyle -> Directive[Thin, Gray], PlotPoints -> 40]] &,
       {ws, labels}]],
     BarLegend[{raw, {-m, m}}, LegendMarkerSize -> 130]]]];
wignerPlot[w_InterpolatingFunction, label_, color_ : Automatic] := wignerPlot[{w}, {label}, color];
```

### A Cavity Cat, and Why It Decoheres Almost at Once

**The problem.** A single leaky mode of light (a cavity dribbling photons), with no drive:
$$\dot\rho = \gamma\,\mathcal{D}[\hat a]\rho.$$
A lone coherent state only shrinks, staying pure, $|\alpha\rangle \to |\alpha e^{-\gamma t/2}\rangle$, but a *cat* $|\alpha\rangle + |\beta\rangle$ of two well-separated amplitudes keeps its cross-term only through a coherence factor
$$|C(t)| = \exp\!\bigl[-\tfrac12\,|\alpha-\beta|^2\,(1 - e^{-\gamma t})\bigr] \;\approx\; \exp\!\bigl(-\tfrac12\,|\alpha-\beta|^2\,\gamma t\bigr)\ \text{at short time,}$$
so the coherence decays at a rate $\propto |\alpha-\beta|^2$, quadratic in the phase-space separation. Build it in a truncated Fock basis and watch the fringes vanish.

Set a high cutoff and a cavity decay rate:

```wl
topCat = 30; blankCat = IdentityMatrix[topCat]; downCat = annihilation[topCat]; \[Gamma]cat = 0.7;
```

First a single coherent state: under damping it stays coherent, its amplitude decaying as $\alpha \to \alpha e^{-\gamma t/2}$. Check that the damped state equals a fresh, smaller coherent state:

```wl
With[{later = evolve[0 blankCat, {Sqrt[\[Gamma]cat] downCat}, densityMatrix[coherentState[topCat, 2.]], 1.0]},
 Re[Conjugate[coherentState[topCat, 2. Exp[-\[Gamma]cat/2]]] . later . coherentState[topCat, 2. Exp[-\[Gamma]cat/2]]]]
```

Overlap one: a coherent state stays coherent.

Now the interesting state, a **cat**, the superposition of two coherent states $|\alpha\rangle + |{-}\alpha\rangle$:

```wl
\[Alpha] = 2.; cat = Normalize[coherentState[topCat, \[Alpha]] + coherentState[topCat, -\[Alpha]]];
```

First a rigorous check on the solver. Amplitude damping has an exact operator-sum (Kraus) solution, summed over the number of photons lost $m$,
$$\rho(t) = \sum_{m\ge 0} E_m\,\rho_0\,E_m^\dagger, \qquad E_m = \sqrt{\frac{(1-e^{-\gamma t})^{m}}{m!}}\;e^{-\gamma t\,\hat a^\dagger\hat a/2}\,\hat a^{m},$$
with $\hat a^m$ removing $m$ photons and $e^{-\gamma t\,\hat a^\dagger\hat a/2}$ shrinking the survivors; the set is trace-preserving, so it needs no renormalization. Write it:

```wl
photonLoss[rho0_, t_] := With[{shrink = DiagonalMatrix[Exp[-\[Gamma]cat t Range[0, topCat - 1]/2]], keep = 1 - Exp[-\[Gamma]cat t]},
   Sum[With[{M = keep^(m/2)/Sqrt[m!] shrink . If[m == 0, blankCat, MatrixPower[downCat, m]]},
      M . rho0 . ConjugateTranspose[M]], {m, 0, topCat - 1}]];
```

Confirm the master equation and this Kraus sum agree on the cat:

```wl
Max@Abs@Flatten[evolve[0 blankCat, {Sqrt[\[Gamma]cat] downCat}, densityMatrix[cat], 0.6] - photonLoss[densityMatrix[cat], 0.6]]
```

They agree to rounding.

Now visualize the cat's decoherence. Build the Wigner film over a short window:

```wl
catTimes = Range[0, .4, .1];
catFilm = Table[wignerRepresentation[evolve[0 blankCat, {Sqrt[\[Gamma]cat] downCat}, densityMatrix[cat], t], 4.2, 55], {t, catTimes}];
```

Draw it on one shared color scale:

```wl
wignerPlot[catFilm, "time " <> ToString[#] & /@ catTimes]
```

At first there are two Gaussian lobes (the two coherent states) with interference fringes between them: the fringes are the coherence, the off-diagonal cross-term. Across this short window the fringes fade fast while the lobes barely shrink: the cat is collapsing toward a classical mixture of the two coherent states while both are still nearly full size. Decoherence is far faster than energy loss.

How fast? The coherence $|C(t)|$ is the closed-form envelope from the problem statement, and the master equation can be made to report the same quantity on its own. Evolve the cat once over the plotted span:

```wl
catRun = evolveODE[0 blankCat, {Sqrt[\[Gamma]cat] downCat}, densityMatrix[cat], 2.];
```

Under amplitude damping the state stays in the span of its two shrinking lobes,
$$\rho(t) = p\,\bigl(|\alpha_t\rangle\langle\alpha_t| + |{-\alpha_t}\rangle\langle{-\alpha_t}|\bigr) + c\,\bigl(|\alpha_t\rangle\langle{-\alpha_t}| + |{-\alpha_t}\rangle\langle\alpha_t|\bigr), \qquad \alpha_t = \alpha e^{-\gamma t/2},$$
and the coherence is the weight ratio $|C(t)| = c/p$. The two lobes are not orthogonal (they overlap by $o = \langle\alpha_t|{-\alpha_t}\rangle$), so read $p$ and $c$ off two matrix elements of $\rho(t)$ by undoing that overlap with a two-by-two linear solve:

```wl
catCoherence[t_?NumericQ] := Module[{va, vb, o, m1, m2, p, cc},
   va = coherentState[topCat, \[Alpha] Exp[-\[Gamma]cat t/2]]; vb = coherentState[topCat, -\[Alpha] Exp[-\[Gamma]cat t/2]];
   o = Re[Conjugate[va] . vb];
   m1 = Re[Conjugate[va] . catRun[t] . va]; m2 = Re[Conjugate[va] . catRun[t] . vb];
   {p, cc} = LinearSolve[{{1 + o^2, 2 o}, {2 o, 1 + o^2}}, {m1, m2}];
   cc/p];
```

Now visualize both, the closed-form envelope and the master-equation extraction, for our chosen separation:

```wl
Plot[{Exp[-2 \[Alpha]^2 (1 - Exp[-\[Gamma]cat t])], catCoherence[t]}, {t, 0, 2},
 PlotStyle -> {ColorData[97, 1], Directive[ColorData[97, 2], Dashed]},
 PlotLegends -> {"closed form", "master equation"}, Frame -> True, GridLines -> Automatic,
 FrameLabel -> {"time", "coherence |C(t)|"}, PlotRange -> {0, 1}, ImageSize -> 460,
 PlotLabel -> "master-equation coherence lies on the closed-form envelope"]
```

The dashed master-equation extraction lies on the closed form: the envelope is not an approximation but the exact law of the simulated dynamics. The decoherence rate grows as the *square* of the separation, $|\alpha - (-\alpha)|^2 = 4|\alpha|^2$.

Every Fock-space example truncates the infinite number-state ladder at a finite cutoff, so each needs a check that the cutoff is high enough not to change the answer. Recompute the mean photon number of a decaying coherent state at cutoffs 20, 30, and 40:

```wl
With[{observable = Function[nn, expectation[creation[nn] . annihilation[nn], evolve[0 IdentityMatrix[nn],
      {Sqrt[\[Gamma]cat] annihilation[nn]}, densityMatrix[coherentState[nn, 2.]], 0.5]]]},
 observable /@ {20, 30, 40}]
```

The three cutoffs give the same mean photon number, so the truncation is fine enough not to affect the answer: 20 already holds the whole state, and raising it changes nothing.

### Making a Cat With an Atom, and Reading Its Decoherence

**The problem.** A *far-detuned* atom, its transition frequency sitting far from the cavity's so it can trade no energy with the light and only turns its phase, imprints a state-dependent phase shift under
$$\hat H_C = \chi\,\hat a^\dagger\hat a\,\hat\sigma_z,\qquad |\alpha\rangle \to |\alpha e^{\pm i\phi}\rangle,\quad \phi = \chi\tau.$$
An atom prepared in the equal superposition $(|e\rangle + |g\rangle)/\sqrt2$ therefore leaves the field in a two-phase **cat**, $|\alpha e^{-i\phi}\rangle + |\alpha e^{+i\phi}\rangle$, entangled with the atom. After the cat decoheres for a delay (the leak of the last example), a second atom probes the same field, and the two-atom correlation $\eta = p(e\,|\,e) - p(e\,|\,g)$ (the second atom's chance of coming out excited when the first did, minus when the first came out in the ground state) is nonzero only while the field keeps its coherence, so $\eta(\tau)$ tracks the surviving coherence $C(\tau)$ and falls to zero as the cat blurs into a mixture. This is the Haroche experiment; build the atom-made cat and measure how fast the correlation dies.

In plainer terms, the atom is a phase stamp, not an absorber: far off resonance it takes no photon from the field, it only rotates the field's phase by $+\phi$ or $-\phi$ according to its state, so an atom in $(|e\rangle+|g\rangle)/\sqrt2$ leaves the two phase-rotated fields superposed and entangled with it. One probe atom cannot reveal this coherence, since its own outcome averages over the two phases; only the correlation $\eta$ between two atoms reading the same field exposes it, and only while it survives. Set the cutoff and the cavity's decay rate:

```wl
topHar = 28; blankHar = IdentityMatrix[topHar]; downHar = annihilation[topHar]; \[Gamma]har = 0.5;
```

The **trace distance** $\tfrac12\mathrm{Tr}|\rho_1 - \rho_2|$ will measure how far the cat sits from a classical mixture:

```wl
traceDistance[a_, b_] := Total[SingularValueList[a - b]]/2;
```

Because the singular values of any matrix sum to its trace norm, `Total[SingularValueList[a - b]]/2` is exactly $\tfrac12\mathrm{Tr}|a-b|$; for the Hermitian difference of two states those singular values are just the magnitudes of its eigenvalues.

The two-phase cat, the equal superposition of the two rotated coherent states $2 e^{\pm i\phi}$:

```wl
catRho[\[Phi]_] := densityMatrix@Normalize[coherentState[topHar, 2. Exp[-I \[Phi]]] + coherentState[topHar, 2. Exp[I \[Phi]]]];
```

The phase $\phi$ sets the separation: small $\phi$ leaves the two coherent states close, large $\phi$ far apart.

Now watch the two cats decohere under the same leak. Evolve each under the cavity's master equation across the delays:

```wl
harTimes = {0., 0.5, 1.0};
harSmall = evolveODE[0 blankHar, {Sqrt[\[Gamma]har] downHar}, catRho[0.5], Max[harTimes]];
harLarge = evolveODE[0 blankHar, {Sqrt[\[Gamma]har] downHar}, catRho[1.3], Max[harTimes]];
```

Build one Wigner film per phase:

```wl
harFilmSmall = Table[wignerRepresentation[harSmall[t], 4.0, 45], {t, harTimes}];
harFilmLarge = Table[wignerRepresentation[harLarge[t], 4.0, 45], {t, harTimes}];
```

Draw both films on one shared color scale so the two separations are directly comparable, the small-phase triple then the large-phase triple:

```wl
wignerPlot[Join[harFilmSmall, harFilmLarge],
 Join["small \[Phi], \[Tau]=" <> ToString[#] & /@ harTimes, "large \[Phi], \[Tau]=" <> ToString[#] & /@ harTimes]]
```

Read the two triples against each other: the small-phase cat keeps its fringes across the whole window, while the large-phase cat, its two states farther apart, has collapsed to a fringeless mixture by the end. Same cavity, same leak, same clock, one color scale; only the separation differs, and the fringes die at a rate set by its square.

The second atom reads this decay as a number. Measure the surviving coherence as the trace distance between the decohering cat and the classical mixture of its two coherent states: full coherence when the cat is intact, zero once it has decohered. The mixture the cat decays toward is the equal blend of the two pointer states, each amplitude shrinking with the cavity,
$$\rho_{\mathrm{mix}}(\phi,t) = \tfrac12\big(|\alpha_+\rangle\langle\alpha_+| + |\alpha_-\rangle\langle\alpha_-|\big), \qquad \alpha_\pm(t) = 2\,e^{-\gamma t/2}\,e^{\pm i\phi},$$
with $\gamma = \gamma_{\mathrm{har}}$:

```wl
mixture[\[Phi]_, t_] := With[{sh = Exp[-\[Gamma]har t/2]},
   (densityMatrix[coherentState[topHar, 2. sh Exp[-I \[Phi]]]] + densityMatrix[coherentState[topHar, 2. sh Exp[I \[Phi]]]])/2];
```

The surviving coherence is the trace distance between the leak-evolved cat $\rho_{\mathrm{cat}}(t) = e^{\gamma t\,\mathcal{D}[\hat a]}[\rho_{\mathrm{cat}}(\phi)]$ and that mixture,
$$\mathrm{coherenceLeft}(\phi,t) = \tfrac12\,\mathrm{Tr}\big|\,\rho_{\mathrm{cat}}(t) - \rho_{\mathrm{mix}}(\phi,t)\,\big|,$$
computed over a list of delays, one master-equation run per phase:

```wl
coherenceLeft[\[Phi]_, delays_] := With[{run = evolveODE[0 blankHar, {Sqrt[\[Gamma]har] downHar}, catRho[\[Phi]], Max[delays]]},
   Table[traceDistance[run[t], mixture[\[Phi], t]], {t, delays}]];
```

Track it for the small and the large phase:

```wl
delays = Range[0, 3, 0.25];
{close, far} = {coherenceLeft[0.5, delays], coherenceLeft[1.3, delays]};
```

Now visualize the two normalized coherences against the delay:

```wl
ListLinePlot[{Transpose[{delays, close/First@close}], Transpose[{delays, far/First@far}]},
 Frame -> True, GridLines -> Automatic, PlotMarkers -> Automatic,
 PlotLegends -> {"small phase shift", "large phase shift"}, PlotRange -> All,
 FrameLabel -> {"delay before the second atom reads", "surviving coherence"},
 PlotLabel -> "cat coherence decays faster for a larger phase shift"]
```

Both fall, but the large-$\phi$ curve plunges: the same separation-squared decoherence law as before, now with the separation set by the atom's phase shift. The trace distance plotted here is a theoretical proxy for what the experiment actually reads out, the two-atom correlation $\eta$ above; both track the surviving coherence $C(\tau)$ and fall to zero together, and the faster collapse at larger $\phi$ is the experimental signature of the scaling.

### Dispersive Readout: A Qubit Measured by a Cavity

**The problem.** A switch $\hat\sigma_z$ read out by a driven, damped oscillator whose frequency it shifts, a solid-state meter:
$$\dot\rho = -i[\epsilon(\hat a + \hat a^\dagger),\rho] - i[\chi\,\hat a^\dagger\hat a\,\hat\sigma_z,\rho] + \gamma\,\mathcal{D}[\hat a]\rho.$$
The drive $\epsilon$ and damping $\gamma$ settle the oscillator into one of two switch-conditioned coherent states, a pointer pointing one way or the other, while the oscillator's own leak decoheres any which-way superposition of the switch. Build the joint switch-plus-light system and watch the pointer split in two and the switch blur.

A qubit is read out by a driven, damped cavity whose frequency it shifts (dispersive coupling). The cavity settles into one of two qubit-conditioned coherent states, the **pointer states**; the cavity field is the meter, its amplitude the pointer. Set the drive, coupling, and decay:

```wl
\[Epsilon] = 1.2; \[Chi] = 0.6; \[Gamma]disp = 1.0;
```

The pointer amplitude obeys a linear driven-damped equation, one per qubit state $s = \pm1$, namely $\dot\alpha = -i\epsilon - (is\chi + \tfrac{\gamma}{2})\alpha$ with $\alpha(0)=0$, so in closed form
$$\alpha_s(t) = \frac{-i\epsilon}{is\chi + \gamma/2}\Big(1 - e^{-(is\chi + \gamma/2)\,t}\Big),$$
settling to the mirror-image steady amplitudes $\alpha_s^{\mathrm{ss}} = -i\epsilon/(is\chi + \gamma/2)$. Solve it:

```wl
pointer[which_] := First[amp /. DSolve[{amp'[t] == -I \[Epsilon] - (I which \[Chi] + \[Gamma]disp/2) amp[t], amp[0] == 0}, amp, t]];
```

Evaluate the two pointers, one per switch state:

```wl
{up, dn} = {pointer[1], pointer[-1]};
```

Read the two late-time amplitudes:

```wl
Chop[{up[9.], dn[9.]}]
```

The two steady amplitudes are mirror images, one per qubit state.

Now visualize the pointer-state separation by drawing both trajectories from the origin to their steady points:

```wl
ListLinePlot[{Table[Sqrt[2] {Re@up[t], Im@up[t]}, {t, 0, 9, 0.03}], Table[Sqrt[2] {Re@dn[t], Im@dn[t]}, {t, 0, 9, 0.03}]},
 PlotStyle -> {ColorData[97, 1], ColorData[97, 2]}, AspectRatio -> 1, ImageSize -> 360,
 PlotLegends -> {"pointer if switch up", "pointer if switch down"}, Frame -> True, GridLines -> Automatic,
 FrameLabel -> {"x", "p"}, PlotLabel -> "the meter pointer settles to two places, one per switch state",
 Epilog -> {PointSize[0.03], ColorData[97, 1], Point[Sqrt[2] {Re@up[9.], Im@up[9.]}],
   ColorData[97, 2], Point[Sqrt[2] {Re@dn[9.], Im@dn[9.]}]}]
```

What a real readout does with those two pointers is discriminate them in the IQ plane. Detecting the cavity output (heterodyne) turns each qubit-conditioned coherent pointer into a Gaussian blob of measurement outcomes, centered at its amplitude with the shot-noise width of a coherent state. As the pointers separate the two blobs pull apart and the assignment becomes reliable. Sample the two outcome clouds at three readout times:

```wl
iqCloud[ptr_, t_, seed_] := BlockRandom[SeedRandom[seed];
   Table[Sqrt[2] {Re@ptr[t], Im@ptr[t]} + RandomVariate[NormalDistribution[0, 1], 2], 250]];
iqTimes = {0.5, 1.5, 4.};
```

Now visualize the readout separating in time, the two clouds (up in blue, down in orange) drifting apart across the IQ plane:

```wl
GraphicsRow[Table[
   ListPlot[{iqCloud[up, t, 1], iqCloud[dn, t, 2]}, PlotStyle -> {Directive[Opacity[0.5], PointSize[0.015]], Directive[Opacity[0.5], PointSize[0.015]]},
    PlotRange -> {{-6, 6}, {-6, 6}}, AspectRatio -> 1, Frame -> True, GridLines -> Automatic, ImageSize -> 260,
    FrameLabel -> {"I", "Q"}, PlotLabel -> "t = " <> ToString[t]], {t, iqTimes}], ImageSize -> 800]
```

Early the clouds overlap and the qubit state is ambiguous; by the last frame they are cleanly resolved. Turn that into numbers: the separation $|\alpha_\uparrow - \alpha_\downarrow|$ sets the single-shot assignment error $\tfrac12\mathrm{erfc}(|\alpha_\uparrow - \alpha_\downarrow|/2)$, the overlap of the two Gaussians. Track both against readout time:

```wl
sepDisp[t_] := Abs[up[t] - dn[t]];
ListLinePlot[{Table[{t, sepDisp[t]}, {t, 0, 6, 0.05}], Table[{t, 0.5 Erfc[sepDisp[t]/2]}, {t, 0, 6, 0.05}]},
 PlotStyle -> {ColorData[97, 1], ColorData[97, 3]}, Frame -> True, GridLines -> Automatic, ImageSize -> 460,
 PlotLegends -> {"pointer separation |\[Alpha]\[UpArrow]-\[Alpha]\[DownArrow]|", "assignment error 0.5 erfc(sep/2)"},
 FrameLabel -> {"readout time", "separation / error"}, PlotLabel -> "longer readout separates the pointers and lowers the error"]
```

The separation rises to its steady value while the assignment error falls toward a few percent: longer integration buys distinguishability. This is exactly the single-shot readout curve a circuit-QED experiment optimizes, the flip side of the coherence the same measurement destroys.

Measurement has a cost: as the pointer separates and reveals the qubit state, the qubit's coherence is destroyed. To see it we need the joint qubit-cavity density matrix and two helpers. First, embed a cavity operator on the joint space:

```wl
onAtom[metering_] := KroneckerProduct[IdentityMatrix[2], metering];
```

Second, the partial trace over the cavity, recovering the qubit's reduced state:

```wl
traceMode[joint_, rungs_] := With[{r = ArrayReshape[joint, {2, rungs, 2, rungs}]},
   Table[Sum[r[[a, k, b, k]], {k, rungs}], {a, 2}, {b, 2}]];
```

Build a reader for the qubit coherence $|\rho_{eg}|$ as the pointer resolves the qubit:

```wl
qubitCoherence[rungs_] := Module[{count = creation[rungs] . annihilation[rungs], H, run},
   H = \[Epsilon] onAtom[annihilation[rungs] + creation[rungs]] + \[Chi] KroneckerProduct[Z, count];
   run = evolveODE[H, {onAtom[Sqrt[\[Gamma]disp] annihilation[rungs]]},
     densityMatrix[Flatten@KroneckerProduct[plus, coherentState[rungs, 0]]], 4.];
   Function[t, 2 Abs[traceMode[run[t], rungs][[1, 2]]]]];
```

Confirm the cutoff by comparing two truncations at the same time:

```wl
reading = qubitCoherence[9];
{reading[2.], qubitCoherence[12][2.]}
```

The two cutoffs agree, so 9 is enough.

Now visualize the qubit coherence draining as the pointer states separate:

```wl
Plot[reading[t], {t, 0, 4}, Frame -> True, GridLines -> Automatic, PlotRange -> {0, 1}, ImageSize -> 420,
 FrameLabel -> {"time", "qubit coherence"},
 PlotLabel -> "the pointer measures away the qubit coherence"]
```

The coherence falls to near zero as the pointer states separate: a measurement that distinguishes the qubit states must destroy their superposition.

The cavity field, the meter itself, tells the same story in its Wigner picture, and what it shows depends on how you read the qubit. Trace the qubit out and the field is a fringeless mixture of the two pointer states as soon as they separate, since the which-way record now lives in the qubit. But keep only the records where the qubit is found in $|+\rangle$ (the basis where its coherence lives), and the conditional field is a cat of the two pointers, its fringes carrying exactly the qubit's surviving coherence. A helper for the field conditioned on finding the qubit in a given state $|q\rangle$, the qubit sandwiched off and the result renormalized,
$$\rho_{\mathrm{field}\,|\,q} = \frac{\langle q|\rho|q\rangle}{\mathrm{Tr}\,\langle q|\rho|q\rangle},$$
where $\langle q|\rho|q\rangle$ leaves a field operator and its trace is the probability $\langle q|\rho_{\mathrm{qubit}}|q\rangle$ of that outcome:

```wl
fieldGiven[joint_, rungs_, q_] := With[{r = ArrayReshape[joint, {2, rungs, 2, rungs}]},
   With[{m = Sum[Conjugate[q[[a]]] q[[b]] r[[a, All, b, All]], {a, 2}, {b, 2}]}, m/Tr[m]]];
```

Run the joint master equation once, at a cutoff wide enough for the pointer blobs:

```wl
dispRungs = 12;
dispRun = evolveODE[\[Epsilon] onAtom[annihilation[dispRungs] + creation[dispRungs]] + \[Chi] KroneckerProduct[Z, creation[dispRungs] . annihilation[dispRungs]],
   {onAtom[Sqrt[\[Gamma]disp] annihilation[dispRungs]]},
   densityMatrix[Flatten@KroneckerProduct[plus, coherentState[dispRungs, 0]]], 3.];
```

Build the Wigner film of the $|+\rangle$-conditioned field as the measurement proceeds:

```wl
dispTimes = {0., 1., 2., 3.};
dispFilm = Table[wignerRepresentation[fieldGiven[dispRun[t], dispRungs, plus], 3.8, 45], {t, dispTimes}];
```

Draw it on one shared color scale:

```wl
wignerPlot[dispFilm, "time " <> ToString[#] & /@ dispTimes]
```

The field starts as the vacuum, one blob at the origin, and drifts as the drive builds the pointer amplitude. As the two pointer components begin to separate, the conditional field passes through a faint cat: an elongated blob with a patch of genuine negativity beneath it, the fringe that the qubit's surviving coherence can still paint. It never becomes a bright cat, and that is the physics: the meter resolves its two readings exactly by carrying the which-way information away, so by the time the pointers are distinct (last panel) the coherence that would fringe them is spent, and the field is a classical mixture of two readings. Contrast the Haroche cat of the earlier section: there the dispersive atom imprinted the two phases *before* the cavity leaked, so the separation came first and the decoherence after; here they are one and the same process. This is dispersive qubit readout in circuit QED, where the pointer is a microwave tone reflected off a resonator.

## Part Three: A Heavy Particle, and the Two Faces of a Warm Bath

**The problem.** A particle of mass $M$ in a gentle trap $\tfrac12 M\Omega^2\hat X^2$, coupled to a warm bath (the Caldeira-Leggett model). In the high-temperature limit,
$$\dot\rho = -i\Bigl[\frac{\hat P^2}{2M} + \tfrac12 M\Omega^2\hat X^2,\ \rho\Bigr] - i\gamma\bigl[\hat X,\{\hat P,\rho\}\bigr] - 2\gamma M k_B T\,\bigl[\hat X,[\hat X,\rho]\bigr].$$
Friction damps the momentum at rate $2\gamma$; the double commutator $[\hat X,[\hat X,\rho]]$ decoheres a superposition of two positions a distance $\Delta x$ apart at a rate $\propto M k_B T\,\Delta x^2$, quadratic in the separation, the same law as the light-cat. For a smooth (Gaussian) blob the first and second moments of $(\hat X,\hat P)$ close into a small set of rate equations, exactly the five below.

A heavy particle in a harmonic trap $\tfrac12 M\Omega^2\hat X^2$, coupled to a warm bath (Caldeira-Leggett). Friction damps the momentum; the bath's fluctuations grow the position variance to thermal equilibrium. For a Gaussian state the first and second moments of $(\hat X,\hat P)$ close into five equations: the mean $(\langle X\rangle, \langle P\rangle)$ and the covariance matrix $\Sigma$ with entries $(\Sigma_{xx}, \Sigma_{xp}, \Sigma_{pp})$. The high-temperature master equation above holds only when $k_B T \gg \hbar\Omega$, so fix the mass, trap frequency, friction, and a temperature well above the level spacing $\hbar\Omega = 1$:

```wl
M = 1.; \[CapitalOmega]qbm = 1.; \[Gamma]qbm = 0.15; kT = 5.; Ddiff = 4 \[Gamma]qbm M kT;
```

With $k_B T = 5$ against a level spacing of $1$, the bath is genuinely warm. Written out, those five equations are
$$\dot{\langle X\rangle} = \frac{\langle P\rangle}{M}, \qquad \dot{\langle P\rangle} = -M\Omega^2\langle X\rangle - 2\gamma\langle P\rangle,$$
$$\dot\Sigma_{xx} = \frac{2\Sigma_{xp}}{M}, \qquad \dot\Sigma_{xp} = \frac{\Sigma_{pp}}{M} - M\Omega^2\Sigma_{xx} - 2\gamma\Sigma_{xp}, \qquad \dot\Sigma_{pp} = -2M\Omega^2\Sigma_{xp} - 4\gamma\Sigma_{pp} + D,$$
with diffusion coefficient $D = 4\gamma M k_B T$. The initial covariance must itself be a valid quantum state, so start from the trap's ground state, the minimum-uncertainty blob with $\Sigma_{xx} = \Sigma_{pp} = \hbar/2 = 1/2$ (so $\det\Sigma = \hbar^2/4 = 1/4$), displaced off the origin. Integrate them from there:

```wl
particle = NDSolveValue[{
   cx'[t] == cp[t]/M, cp'[t] == -M \[CapitalOmega]qbm^2 cx[t] - 2 \[Gamma]qbm cp[t],
   vx'[t] == 2 vc[t]/M, vc'[t] == vp[t]/M - M \[CapitalOmega]qbm^2 vx[t] - 2 \[Gamma]qbm vc[t],
   vp'[t] == -2 M \[CapitalOmega]qbm^2 vc[t] - 4 \[Gamma]qbm vp[t] + Ddiff,
   cx[0] == 2., cp[0] == 0., vx[0] == 0.5, vc[0] == 0., vp[0] == 0.5},
  {cx, cp, vx, vc, vp}, {t, 0, 30}];
```

Two readers, the mean and the covariance matrix:

```wl
center[t_] := {particle[[1]][t], particle[[2]][t]};
spread[t_] := {{particle[[3]][t], particle[[4]][t]}, {particle[[4]][t], particle[[5]][t]}};
```

Two checks. First, the mean spirals in from its displaced start to rest:

```wl
{center[0], center[25]}
```

It starts displaced and ends near the origin. Second, the position variance $\Sigma_{xx}$ grows from its ground-state value $\hbar/2 = 1/2$ to the thermal value $k_B T/(M\Omega^2)$, equipartition for a harmonic trap:

```wl
{particle[[3]][0], particle[[3]][25], kT/(M \[CapitalOmega]qbm^2)}
```

It lands on the thermal value, ten times the ground-state spread. Check it a second way: the steady covariance solves the Lyapunov equation, where diffusion balances friction:

```wl
drift = {{0, 1/M}, {-M \[CapitalOmega]qbm^2, -2 \[Gamma]qbm}};
LyapunovSolve[drift, -{{0, 0}, {0, Ddiff}}]
```

The Lyapunov steady state matches the integrated one.

The same covariance bookkeeping returns under monitoring in Part Five. The Kalman filter's watched oscillator obeys $\dot\Sigma = A\Sigma + \Sigma A^{T} + D - \Sigma C^{T} C\Sigma$, the Lyapunov form above plus one term the bath cannot supply, the measurement's information gain $-\Sigma C^{T} C\Sigma$: diffusion widens the blob to thermal here, watching narrows it below the ground-state spread there.

One more check the model must pass to be believed: a real quantum state obeys the Robertson-Schrodinger uncertainty relation $\det\Sigma = \Sigma_{xx}\Sigma_{pp} - \Sigma_{xp}^2 \ge \hbar^2/4 = 1/4$, so the margin $\det\Sigma - 1/4$ must stay non-negative for the whole run. The high-temperature Caldeira-Leggett generator is not of Lindblad form, so this is not automatic; watch it directly:

```wl
Plot[Det@spread[t] - 1/4, {t, 0, 25}, Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"time", "uncertainty margin  det \[CapitalSigma] - 1/4"},
 PlotLabel -> "the state stays above the uncertainty floor"]
```

The margin starts at zero, where the ground state saturates the bound, and only grows: the blob swells above the Heisenberg floor and never dips below it. Had we started from the old sub-ground blob $\Sigma_{xx} = \Sigma_{pp} = 0.1$, the margin would have opened at $0.01 - 0.25 < 0$, an invalid state the evolution could never repair; and pushed to low temperature, this same generator can drive an initially valid state below the floor, the known price of its non-Lindblad form.

Those checks read the mean and the variance as numbers; the phase-space picture puts them in one object. That object is the state's **Wigner function** $W(x,p)$, the quasiprobability the cat used in Part Two, except that for a Gaussian state no Fock-basis sum is needed: $W$ is exactly the bivariate normal set by the five moments, centered at the mean $\bar r = (\langle X\rangle, \langle P\rangle)$ and shaped by the covariance $\Sigma$,
$$W(x,p) = \frac{1}{2\pi\sqrt{\det\Sigma}}\,\exp\!\left[-\tfrac12\,(r - \bar r)^{T}\,\Sigma^{-1}\,(r - \bar r)\right], \qquad r = (x,p),$$
normalized so $\int W\,dx\,dp = 1$. Transcribe it as a function of time, reading $\bar r$ and $\Sigma$ off the solution at each $t$:

```wl
wignerGaussian[t_?NumericQ] := With[
   {r0 = center[t], iS = Inverse[spread[t]], nrm = 1/(2 Pi Sqrt[Det[spread[t]]])},
   Function[{x, p}, nrm Exp[-(({x, p} - r0) . iS . ({x, p} - r0))/2]]];
```

Trace the mean's whole path once, a faint guide to lay under the moving blob:

```wl
centerLine = ListLinePlot[Table[center[t], {t, 0, 25, 0.05}], PlotStyle -> Directive[Thin, Gray]];
```

Now set it in motion from the displaced start: the Wigner bump rides the gray path in toward the origin as friction drains the mean, and fills out to the thermal size as diffusion widens it, the red dot marking where the mean sits at each instant. Press play:

```wl
Animate[Show[
   DensityPlot[wignerGaussian[t][x, p], {x, -6, 6}, {p, -6, 6}, PlotPoints -> 60, PlotRange -> All,
    ColorFunction -> "SunsetColors", AspectRatio -> 1, Frame -> True, ImageSize -> 400,
    FrameLabel -> {"x", "p"}, PlotLabel -> "Wigner W(x,p) at t = " <> ToString[t],
    Epilog -> {Red, Point[center[t]]}], centerLine],
 {t, 0, 25, .1}, AnimationRunning -> False]
```

It stays a single Gaussian bump the whole way, never developing finer structure: the state never leaves the five-number family, the closure that let five equations describe it in the first place. Friction and diffusion are the same bath's two effects, tied by the fluctuation-dissipation relation: the friction rate $\gamma$ is set by the coupling, while the diffusion (and so the thermal spread the blob swells to) grows with temperature.

The sharper effect, and why this sits beside the cavity cat: a superposition of two positions decoheres, and a warm bath does it faster than a cold one. Here the model changes from the Caldeira-Leggett moment equations above to the thermal-oscillator Lindblad master equation (a decay channel $\hat a$ at rate $\gamma(n_T+1)$ and an excitation channel $\hat a^\dagger$ at rate $\gamma n_T$, the same warm-bath form used later for thermalization), integrated in a truncated Fock basis. Build the same coherent-state cat and evolve it under that bath. Set the oscillator cutoff and the bath's rate:

```wl
twoRung = 22; \[Gamma]warm = 0.5;
```

The surviving coherence is the trace distance between the coherently evolved cat and the classical mixture its two lobes would form, both carried through the same thermal bath $\mathcal{L}\rho = \gamma(n_T+1)\,\mathcal{D}[\hat a]\rho + \gamma n_T\,\mathcal{D}[\hat a^\dagger]\rho$:
$$\mathrm{warmCoherence}(n_T,\alpha,t) = \tfrac12\,\mathrm{Tr}\big|\,\rho_{\mathrm{cat}}(t) - \tfrac12\big[\rho_{+\alpha}(t) + \rho_{-\alpha}(t)\big]\,\big|,\quad \rho_{\mathrm{cat}}(t) = e^{\mathcal{L}t}\big[|\mathrm{cat}\rangle\langle\mathrm{cat}|\big],\ \ \rho_{\pm\alpha}(t) = e^{\mathcal{L}t}\big[|{\pm\alpha}\rangle\langle{\pm\alpha}|\big],$$
with $\alpha = \text{gap}$. It is largest while the cat is still distinct from that mixture and falls to zero as the cross-terms decohere (both thermalize to the same steady state, so the distance can only shrink). The bath's $\hat a^\dagger$ channel heats a coherent state, so the two lobes will not stay coherent and are evolved rather than written down; compute it over a list of delays, one run each for the cat and its lobes:

```wl
warmCoherence[nT_, gap_, delays_] := With[
   {fall = Sqrt[\[Gamma]warm (nT + 1)] annihilation[twoRung], climb = Sqrt[\[Gamma]warm nT] creation[twoRung],
    pair = densityMatrix@Normalize[coherentState[twoRung, gap] + coherentState[twoRung, -gap]],
    blobA = densityMatrix[coherentState[twoRung, gap]], blobB = densityMatrix[coherentState[twoRung, -gap]]},
   With[{runs = evolveODE[0 IdentityMatrix[twoRung], {fall, climb}, #, Max[delays]] & /@ {pair, blobA, blobB}},
    Table[Total@SingularValueList[runs[[1]][t] - (runs[[2]][t] + runs[[3]][t])/2]/2, {t, delays}]]];
```

Now visualize the temperature dependence by comparing a cold and a warm bath at the same separation:

```wl
warmTimes = Range[0, 1.5, 0.06];
ListLinePlot[{Transpose[{warmTimes, warmCoherence[0.15, 1.5, warmTimes]}], Transpose[{warmTimes, warmCoherence[1.5, 1.5, warmTimes]}]},
 PlotMarkers -> Automatic, Frame -> True, GridLines -> Automatic, PlotRange -> All,
 PlotLegends -> {"cold bath", "warm bath"}, FrameLabel -> {"time", "surviving coherence"},
 PlotLabel -> "a warm bath decoheres the cat faster than a cold one"]
```

The warm bath decoheres the cat far faster: the rate grows with temperature. As with the cavity cat, the rate also grows with the separation: together the decoherence rate scales as $(2\bar n + 1)$ times the *square* of the separation. The $+1$ is the zero-point (vacuum) term that decoheres even a cold cavity cat; temperature adds to it, so heating speeds decoherence but is not what starts it. For macroscopic separations the rate is enormous at any temperature: the fuller answer to why superpositions of large objects become so difficult to preserve as their separation grows. Sweep the cat size continuously at the fixed short delay, and the two effects trade off:

```wl
ListPlot[With[{sep = Range[0, 3, 0.05]}, {#, warmCoherence[1.5, #, {0.15}][[1]]} & /@ sep],
 Frame -> True, GridLines -> Automatic, PlotRange -> All, ImageSize -> 460,
 FrameLabel -> {"cat separation", "surviving coherence at fixed short delay"},
 PlotLabel -> "surviving coherence peaks at a modest cat size, then falls"]
```

The curve is not monotonic: a bigger cat starts with more coherence but decoheres faster, at a rate set by the separation squared, so at a fixed short delay the two trade off and the coherence peaks at a modest size before the quadratic decoherence overruns it. (Stop at separation 3: past there the Fock cutoff clips the coherent states.)

## Part Four: One Leak, Three Ways to Watch It

So far we mostly averaged over the record. Now we keep it. The key fact: the same master equation can be **unravelled** into trajectories in genuinely different ways, each giving different single histories that all average to the same master equation. We watch one driven, emitting atom three ways. Fix the atom, its drive, its decay channel, and the time grid:

```wl
\[CapitalOmega]atom = 3.0; \[Gamma]atom = 1.0; Hdrive = (\[CapitalOmega]atom/2) X; cAtom = Sqrt[\[Gamma]atom] lower;
dtAtom = 0.01; when = dtAtom Range[0, 600];
```

A reader for the excited-state population:

```wl
excitedPop[rho_] := (1 + blochVector[rho][[3]])/2;
```

Integrate the master equation for the population every unravelling must average back to:

```wl
smoothAtom = evolveODE[Hdrive, {cAtom}, densityMatrix[excited], 6.0];
averageChance = excitedPop[smoothAtom[#]] & /@ when;
```

### Counting the Flashes, One at a Time

**The problem.** Any leak $\dot\rho = -i[\hat H,\rho] + \mathcal{D}[\hat c]\rho$ can be *unravelled* by watching its output, each way of detecting the emitted light giving one unravelling into pure-state histories. Counting photons one click at a time gives the jump trajectory
$$d|\psi\rangle = \Bigl[dN\,\Bigl(\tfrac{\hat c}{\sqrt{\langle\hat c^\dagger\hat c\rangle}} - \mathbb{1}\Bigr) + dt\,\Bigl(\tfrac{\langle\hat c^\dagger\hat c\rangle}{2} - \tfrac{\hat c^\dagger\hat c}{2} - i\hat H\Bigr)\Bigr]|\psi\rangle,$$
where the click $dN \in \{0,1\}$ arrives at rate $\langle\hat c^\dagger\hat c\rangle\,dt$: between clicks the norm decays under the non-Hermitian $\hat H - \tfrac{i}{2}\hat c^\dagger\hat c$, and at a click the state jumps by $\hat c$. Averaging $|\psi\rangle\langle\psi|$ over histories returns the master equation exactly. Watch one counted history and check the average.

The first way, **photon counting**: detect each emitted photon. Between clicks the state evolves under the non-Hermitian $\hat H - \tfrac{i}{2}\hat c^\dagger\hat c$; at each click it jumps by $\hat c$ to the ground state. A single run is one sweep that does both jobs at once: `FoldList` advances and keeps the conditioned state at every step (a jump with probability $\langle\hat c^\dagger\hat c\rangle\,dt$, else the no-click drift), while `Sow` writes down the time whenever the jump fires. It returns the inferred state history and the raw click record together, because in a real experiment they are one run, not two. Everything that follows, one trajectory, the ensemble average, and the waiting-time statistic, reads off realizations of this single generator:

```wl
countingRun[seed_, steps_] := BlockRandom[SeedRandom[seed];
   With[{noClick = MatrixExp[-I (Hdrive - I/2 ConjugateTranspose[cAtom] . cAtom) dtAtom]},
    With[{r = Reap@FoldList[
         Function[{psi, n},
          If[RandomReal[] < dtAtom Re[Conjugate[cAtom . psi] . (cAtom . psi)],
            (Sow[n dtAtom]; Normalize[cAtom . psi]),
            Normalize[noClick . psi]]],
         excited, Range[steps]]},
     <|"states" -> First[r], "clicks" -> Flatten[Last[r]]|>]]];
```

The same population, read from a state vector:

```wl
excitedPopK[psi_] := (1 + Re[Conjugate[psi] . Z . psi])/2;
```

Record one counted history over the displayed window, and read its two parts, the conditioned states and the click times, straight off the one run:

```wl
oneRun = countingRun[3, 600]; oneCount = oneRun["states"]; oneClicks = oneRun["clicks"];
```

Now visualize one counted trajectory's excited-state population against the master equation: the trajectory drives up and down and plunges at each emission, while the average remains smooth:

```wl
ListLinePlot[{Transpose[{when, excitedPopK /@ oneCount}], Transpose[{when, averageChance}]},
 PlotStyle -> {ColorData[97, 1], Directive[Thick, Red, Dashed]}, Frame -> True, GridLines -> Automatic,
 PlotLegends -> {"one trajectory (jumps)", "master equation"}, PlotRange -> {0, 1.05}, ImageSize -> 500,
 FrameLabel -> {"time", "excited-state population"}, PlotLabel -> "photon counting: one trajectory jumps, the average is smooth"]
```

The plunges are the emitted photons, one at a time. No single run is smooth; smoothness is in the average.

The plot shows the inferred *state*; the same run also carries the *record*, the list of click times a photon counter actually writes down, already returned as `oneRun["clicks"]` and bound to `oneClicks` above. Every plunge is one entry in it, because the jump that plunges the state is the very step that sowed the time. Count the clicks in this record:

```wl
Length[oneClicks]
```

A handful of photons over the run. Now put the record and the state on one time axis with two vertical scales: the conditioned population on the left, each click marked by a red line, and the cumulative count $N(t)$ on the right, a staircase that rises by one at every click:

```wl
ListLinePlot[{Transpose[{when, excitedPopK /@ oneCount}],
   Transpose[{Riffle[oneClicks, oneClicks], Riffle[Range@Length@oneClicks - 1, Range@Length@oneClicks]}]},
 MultiaxisArrangement -> All, PlotRange -> All, Frame -> True, ImageSize -> 560,
 GridLines -> {{#, Red} & /@ oneClicks, None},
 FrameLabel -> {{"excited population", "cumulative count N(t)"}, {"time", None}},
 PlotLabel -> "clicks (red lines) fall at the population plunges; the staircase counts them"]
```

Every red line sits on a plunge: the click is the emission, and the plunge is the state collapsing to the ground state the instant the photon is registered. The staircase (right scale) is the raw experimental data; the population trajectory (left scale) is what the theorist infers from it.

Now run an ensemble of these same trajectories, carried well past the display window so each click record is long enough for the timing statistic below; our single run above is seed 3 of this same generator, the ensemble's third member seen over just the plotted window:

```wl
nCount = 200; nWaitSteps = 4000;
countRuns = countingRun[#, nWaitSteps] & /@ Range[nCount];
```

Averaging the conditioned states across the ensemble must recover the master equation. Compare the ensemble-mean excited population, over the displayed window, against the integrated reference:

```wl
Max@Abs[(excitedPop /@ Mean[Table[densityMatrix /@ Take[countRuns[[k]]["states"], Length[when]], {k, nCount}]]) - averageChance]
```

The gap shrinks as $1/\sqrt N$: the jump trajectories average back to the master equation.

The click record carries a statistic no smooth curve can: the distribution of *waiting times* between clicks. A click leaves the atom in $|g\rangle$, so the time $\tau$ to the next click has the density
$$w(\tau) = \langle\psi(\tau)|\hat c^\dagger\hat c|\psi(\tau)\rangle, \qquad |\psi(\tau)\rangle = e^{-i\hat H_{\mathrm{eff}}\tau}\,|g\rangle,\quad \hat H_{\mathrm{eff}} = \hat H - \tfrac{i}{2}\hat c^\dagger\hat c,$$
the unnormalized no-click state evolved from the ground state. Since it starts in $|g\rangle$, $w(0) = \gamma\,|\langle e|g\rangle|^2 = 0$: the atom cannot emit twice at once, it must be re-driven to the excited state first, so very short gaps are forbidden. This is antibunching. Read the density directly off the no-click evolution:

```wl
waitDensity[t_] := With[{psi = MatrixExp[-I (Hdrive - I/2 ConjugateTranspose[cAtom] . cAtom) t] . ground},
   Re[Conjugate[psi] . ConjugateTranspose[cAtom] . cAtom . psi]];
```

No separate simulation is needed. Because every click resets the atom exactly to $|g\rangle$, the gap between two consecutive clicks in any trajectory is already a draw from $w(\tau)$: it begins where $w$ begins. Pool those gaps straight from the ensemble's click records, one set of differences per run (the leading interval, which starts from the excited state rather than from $|g\rangle$, drops out on its own, since a difference needs two clicks):

```wl
waitTimes = Flatten[Differences[#["clicks"]] & /@ countRuns];
```

Overlay their histogram on the closed form:

```wl
Show[Histogram[waitTimes, {0, 6, 0.2}, "PDF", ChartStyle -> ColorData[97, 1]],
 Plot[waitDensity[t], {t, 0, 6}, PlotStyle -> Directive[Thick, Red]],
 Frame -> True, GridLines -> Automatic, ImageSize -> 480,
 FrameLabel -> {"waiting time between clicks", "density"},
 PlotLabel -> "the trajectory's own inter-click gaps land on the closed-form density"]
```

The bars land on the closed-form density: it vanishes at zero (no two clicks at once), rises to a first hump near half a Rabi period, then oscillates and decays, the re-driven atom's Rabi cycle written into the timing of the clicks. This antibunching is the hallmark of single-photon (nonclassical) light, and it lives entirely in the click timing, invisible in the smooth population.

### Homodyne: Reading One Quadrature

**The problem.** Mixing the output with a strong reference beam before detecting it (homodyne) turns the clicks into a continuous noisy current and gives a *diffusive* stochastic master equation
$$d\rho = -i[\hat H,\rho]\,dt + \mathcal{D}[\hat c]\rho\,dt + \mathcal{H}[\hat c]\rho\,dW,$$
with a real Wiener increment $dW$ ($dW^2 = dt$). The conditioned state now purifies smoothly instead of by jumps. The point of this section is the bookkeeping: what a single run gives us, which quantities are *observed* and which are *inferred*, and how every homodyne object is one of those two or is built from them.

The second way, **homodyne** detection: mix the output with a strong reference beam and measure one quadrature as a continuous noisy current. Run one homodyne trajectory of the same atom:

```wl
smoothHistory = trajectory[densityMatrix[excited], Hdrive, {cAtom}, {1.}, {}, dtAtom, 6.0, 3];
```

`trajectory` hands back exactly two things, and everything else on this page is one of them or is computed from them. The **conditional state** $\rho_c(t)$ at each time node (`smoothHistory["states"]`) is the observer's running estimate of the atom given the record so far. The **record increment** $dJ$ at each step (`smoothHistory["record"]`) is the raw output of the detector. The model tying them is
$$dJ = \underbrace{\sqrt\gamma\,\langle\hat\sigma_x\rangle}_{\text{signal}}\,dt + dW,$$
the signal read from $\rho_c$, the noise $dW$ a Wiener increment. Hold the two apart: $dJ$ is **observed**, and $\rho_c$ is **inferred**, reconstructed from $dJ$ by the filter (the positivity step from the toolkit).

Take the conditional state first. Its excited population is one smooth diffusing history, no jumps:

```wl
ListLinePlot[{Transpose[{when, excitedPop /@ smoothHistory["states"]}], Transpose[{when, averageChance}]},
 PlotStyle -> {ColorData[97, 1], Directive[Thick, Red, Dashed]}, Frame -> True, GridLines -> Automatic,
 PlotLegends -> {"one conditional history", "master equation"}, PlotRange -> {0, 1.05}, ImageSize -> 500,
 FrameLabel -> {"time", "excited-state population"}, PlotLabel -> "the conditional state diffuses smoothly, no jumps"]
```

Smooth diffusion, not the plunges of a counted trajectory, yet averaging many conditional states recovers the same master equation:

```wl
Max@Abs[(excitedPop /@ Mean[Table[trajectory[densityMatrix[excited], Hdrive, {cAtom}, {1.}, {}, dtAtom, 6.0, k]["states"], {k, 150}]]) - averageChance]
```

The gap is Monte-Carlo scatter: the conditional state is inferred from a record, and its average over all records is exactly the unconditioned master equation.

Now the record itself, and the four quantities read from it, each labeled by how it is obtained and whether it is observed or inferred.

**The record increment $dJ$** (observed). *How obtained:* the detector's output each step, here generated together with the state. *Status:* observed, the one thing that enters from the world; everything else is downstream of it. Plot it:

```wl
dJ = Flatten@smoothHistory["record"];
```

```wl
ListLinePlot[Transpose[{Most@when, dJ}],
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"time", "record increment  dJ"}, PlotLabel -> "dJ, the observed record: almost pure noise per step"]
```

It looks like noise because it is almost all noise: per step the signal $\sqrt\gamma\langle\hat\sigma_x\rangle\,dt \approx 0.005$ sits under a $dW$ of size $\sqrt{dt} = 0.1$, twenty times larger. A single increment carries almost no information about the atom.

**The integrated record $J(t) = \int_0^t dJ$** (observed). *How obtained:* sum the observed increments. *Status:* observed, the same data accumulated. Plot it against its drift, the inferred signal integrated:

```wl
signal = Sqrt[\[Gamma]atom] (blochVector[#][[1]] & /@ Most[smoothHistory["states"]]);
Jt = Prepend[Accumulate[dJ], 0.];
driftJ = Prepend[Accumulate[signal dtAtom], 0.];
```

```wl
ListLinePlot[{Transpose[{when, Jt}], Transpose[{when, driftJ}]},
 Frame -> True, GridLines -> Automatic,
 PlotLegends -> {"integrated record  J(t)", "inferred drift  \[Integral]signal dt"},
 FrameLabel -> {"time", "integrated record  J(t)"}, PlotLabel -> "the observed record integrates to a drifting path; its drift is the signal"]
```

Integrating is itself a filter: the signal accumulates coherently into a drift while the noise only grows as $\sqrt t$, so in $J(t)$ the signal is comparable to the noise instead of twenty times below it. This is the honest, window-free picture of the data. (The `signal` used for the drift is inferred, and gets its own panel below.)

**The current $dJ/dt$** (observed, singular). *How obtained:* divide the observed increment by $dt$. *Status:* observed but not a proper quantity: $dJ/dt = \text{signal} + dW/dt$, and $dW/dt$ is white noise of variance $1/dt$, which diverges as $dt\to0$. On the grid it is a $\pm 1/\sqrt{dt}$ band around a $\pm 1$ signal:

```wl
ListLinePlot[Transpose[{Most@when, dJ/dtAtom}],
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"time", "current  dJ/dt"}, PlotLabel -> "the current dJ/dt: a white-noise band with no dt\[Rule]0 limit"]
```

There is no current *curve*. The instantaneous homodyne current is a distribution, not a function, which is why any smooth current needs an arbitrary averaging window and why $J(t)$ above is the cleaner observed object.

**The signal $\sqrt\gamma\langle\hat\sigma_x\rangle$** (inferred). *How obtained:* read off the conditional state. *Status:* inferred, the conditional mean of the current, never measured directly. This is the clean thing buried in the record:

```wl
ListLinePlot[Transpose[{Most@when, signal}],
 Frame -> True, GridLines -> Automatic, PlotRange -> {-1.05, 1.05},
 FrameLabel -> {"time", "signal  \!\(\*SqrtBox[\(\[Gamma]\)]\) \[LeftAngleBracket]\[Sigma]x\[RightAngleBracket]"},
 PlotLabel -> "the signal, inferred from the conditional state, not measured"]
```

No detector emits this curve. It is what the filter reconstructs, smooth because the conditional state is smooth, and it is exactly the drift the observed $J(t)$ climbed along.

**The innovation $dW$** (inferred). *How obtained:* $dW = dJ - \text{signal}\,dt$, the observed increment minus the inferred prediction. *Status:* inferred, a residual, not a second measurement:

```wl
dW = dJ - signal dtAtom;
ListLinePlot[Transpose[{Most@when, dW}],
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"time", "innovation  dW"}, PlotLabel -> "the innovation dW = dJ - signal dt: the residual"]
```

If the filter is optimal, the residual is zero-mean white noise. Check its autocorrelation, with the $\pm 2/\sqrt N$ band white noise would sit inside:

```wl
With[{innov = dW/Sqrt[dtAtom], ci = 2/Sqrt[Length[dW]]},
 ListPlot[CorrelationFunction[innov, {1, 10}], Filling -> Axis, PlotRange -> {-0.2, 0.2},
  Frame -> True, GridLines -> {None, {{ci, Directive[Gray, Dashed]}, {-ci, Directive[Gray, Dashed]}}},
  PlotMarkers -> Automatic,
  FrameLabel -> {"lag", "innovation autocorrelation"}, PlotLabel -> "the innovation is white within confidence limits"]]
```

Every nonzero lag sits in the band: the filter has extracted every predictable part of the record, and what remains is irreducible measurement noise. This whiteness is the diffusive counterpart of the counting record's antibunched clicks, and the same optimal-filter signature the Kalman tracker will show.

So the section is one ledger. A single observed record $dJ$ (equivalently its integral $J$) enters from the world; from it the conditional state, the signal, and the innovation $dW$ are all computed. Only $dJ$ is data. Everything that looks like the atom, the smooth population, the signal curve, is inferred from it.

### Heterodyne: Reading Both Quadratures

**The problem.** Reading *both* shadows at once (heterodyne) measures both quadratures through a complex noise:
$$d\rho = -i[\hat H,\rho]\,dt + \mathcal{D}[\hat c]\rho\,dt + \mathcal{H}[dZ^*\hat c]\rho,\qquad dZ\,dZ^* = dt,\quad dZ^2 = 0.$$
It is one point in a whole *general-dyne* family running from one sharp shadow (homodyne) to two blurry ones, every member averaging to the identical master equation while giving genuinely different single histories. Show that this particular reading is just two half-strength leaks a quarter-turn apart, then watch.

The third way, **heterodyne**: measure both quadratures at once. It is equivalent to two half-strength homodyne channels a quarter-turn apart, $\hat c/\sqrt2$ and $i\hat c/\sqrt2$. First prove those two dissipators sum to the original, for any $\hat c$ and $\rho$:

```wl
With[{c = Array[Subscript[\[FormalC], ##] &, {2, 2}], 
   r = Array[Subscript[\[FormalR], ##] &, {2, 2}]}, 
  dissipator[c/Sqrt[2], r] + dissipator[I c/Sqrt[2], r] == 
   dissipator[c, r]] // FullSimplify
```

It returns True, with no assumptions: the two channels reproduce the same dissipator, so the master equation is unchanged. The trajectory explores differently. Confirm the ensemble still averages to the master equation:

```wl
Max@Abs[(excitedPop /@ Mean[Table[
      trajectory[densityMatrix[excited], Hdrive, {cAtom/Sqrt[2], I cAtom/Sqrt[2]}, {1., 1.}, {}, dtAtom, 6.0, k]["states"], {k, 150}]]) - averageChance]
```

Again the average matches.

The trajectory explores differently, so hold apart what heterodyne *observes* from what it *infers*, exactly as with homodyne. `trajectory` still returns two objects, the conditional state $\rho_c$ (inferred) and the record (observed), but the record now has *two* columns, one real current per quadrature:
$$dJ_I = \underbrace{\sqrt{\gamma/2}\,\langle\hat\sigma_x\rangle}_{\text{signal}}\,dt + dW_I,\qquad dJ_Q = \underbrace{\sqrt{\gamma/2}\,\langle\hat\sigma_y\rangle}_{\text{signal}}\,dt + dW_Q.$$
Equivalently the two real records are the real and imaginary parts of a single *complex* record $dZ$, whose increment satisfies $dZ\,dZ^* = dt$ and $dZ^2 = 0$ (the complex noise at the head of this section): two independent real noises carried in one complex number. The single new fact over homodyne is that two conjugate quadratures, $\langle\hat\sigma_x\rangle$ and $\langle\hat\sigma_y\rangle$, are read at once. Run one history and keep its record:

```wl
hetOne = trajectory[densityMatrix[excited], Hdrive, {cAtom/Sqrt[2], I cAtom/Sqrt[2]}, {1., 1.}, {}, dtAtom, 6.0, 3];
```

**The conditional state $\rho_c$** (inferred). *How obtained:* the filter's running estimate of the atom given the record so far. *Status:* inferred. Where homodyne resolved one transverse coordinate, heterodyne resolves both; plot the run's $\langle\hat\sigma_x\rangle$ and $\langle\hat\sigma_y\rangle$ together:

```wl
ListLinePlot[{Transpose[{when, blochVector[#][[1]] & /@ hetOne["states"]}], Transpose[{when, blochVector[#][[2]] & /@ hetOne["states"]}]},
 PlotStyle -> {ColorData[97, 1], ColorData[97, 2]}, PlotLegends -> {"\[LeftAngleBracket]\[Sigma]x\[RightAngleBracket]", "\[LeftAngleBracket]\[Sigma]y\[RightAngleBracket]"},
 Frame -> True, GridLines -> Automatic, PlotRange -> {-1.05, 1.05},
 FrameLabel -> {"time", "conditional transverse Bloch vector"}, PlotLabel -> "one heterodyne run resolves both quadratures"]
```

Both components diffuse smoothly, no jumps, and both are live: the run carries a full transverse Bloch vector where the homodyne run knew only one coordinate.

**The record increments $dJ_I, dJ_Q$** (observed). *How obtained:* the detector's two outputs each step. *Status:* observed, the only data that enters from the world; everything else is downstream. Read the two columns:

```wl
dJI = hetOne["record"][[All, 1]];
dJQ = hetOne["record"][[All, 2]];
```

```wl
ListLinePlot[{Transpose[{Most@when, dJI}], Transpose[{Most@when, dJQ}]},
 PlotStyle -> {ColorData[97, 1], ColorData[97, 2]}, PlotLegends -> {"I record", "Q record"},
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"time", "record increments"}, PlotLabel -> "both quadrature records: almost pure noise per step"]
```

Each is nearly pure noise, and for a sharper reason than homodyne's: each channel carries only $\sqrt{\gamma/2}$ of signal, half the power homodyne put into its one quadrature, against the same shot noise.

**The integrated record $J(t) = \int_0^t dZ$** (observed). *How obtained:* accumulate each column. *Status:* observed, the same data summed. A complex record integrates to a path in the $IQ$ plane, an observed walk that drifts along the inferred signal.

That signal is the per-channel coupling $\sqrt{\gamma/2}$ times each conditional quadrature, read off $\rho_c$ rather than measured, one channel for $\langle\hat\sigma_x\rangle$ and one for $\langle\hat\sigma_y\rangle$:

```wl
hetSigI = Sqrt[\[Gamma]atom/2] (blochVector[#][[1]] & /@ Most[hetOne["states"]]);
hetSigQ = Sqrt[\[Gamma]atom/2] (blochVector[#][[2]] & /@ Most[hetOne["states"]]);
```

The observed walk is the running sum of the two record columns, each starting from the origin:

```wl
JI = Prepend[Accumulate[dJI], 0.]; JQ = Prepend[Accumulate[dJQ], 0.];
```

The drift is that same running sum on the inferred signal, the smooth curve the walk fluctuates around:

```wl
driftI = Prepend[Accumulate[hetSigI dtAtom], 0.]; driftQ = Prepend[Accumulate[hetSigQ dtAtom], 0.];
```

Now trace the observed walk against that inferred drift:

```wl
ListLinePlot[{Transpose[{JI, JQ}], Transpose[{driftI, driftQ}]},
 PlotStyle -> {Gray, Directive[Thick, ColorData[97, 3]]}, PlotLegends -> {"integrated record", "inferred drift"},
 Frame -> True, GridLines -> Automatic, AspectRatio -> 1,
 FrameLabel -> {"I", "Q"}, PlotLabel -> "the observed record is a drifting walk in the IQ plane"]
```

The jagged path is the raw data; the smooth curve it drifts along is inferred. Integrating both quadratures at once is a 2D filter: each signal accumulates coherently into the drift while its noise only spreads as $\sqrt t$, so the walk's net displacement carries the record's information where a single increment held almost none.

**The innovations $dW_I, dW_Q$** (inferred). *How obtained:* $dW_I = dJ_I - \text{signal}_I\,dt$ and $dW_Q = dJ_Q - \text{signal}_Q\,dt$, each observed increment minus its inferred prediction. *Status:* inferred residuals. If the filter is optimal each is zero-mean white noise, and because $dZ^2 = 0$ the two are mutually independent. Check the whiteness of each channel first, against the $\pm 2/\sqrt N$ band:

```wl
dWI = dJI - hetSigI dtAtom;
dWQ = dJQ - hetSigQ dtAtom;
```

```wl
With[{ci = 2/Sqrt[Length[dWI]]},
 ListPlot[{CorrelationFunction[dWI/Sqrt[dtAtom], {1, 10}], CorrelationFunction[dWQ/Sqrt[dtAtom], {1, 10}]},
  PlotRange -> {-0.2, 0.2}, PlotMarkers -> Automatic, PlotLegends -> {"I channel", "Q channel"},
  Frame -> True, GridLines -> {None, {{ci, Directive[Gray, Dashed]}, {-ci, Directive[Gray, Dashed]}}},
  FrameLabel -> {"lag", "innovation autocorrelation"}, PlotLabel -> "each innovation channel is white"]]
```

Both correlograms hug zero within the band, the signature of white residuals: the filter has extracted every predictable part of each record.

Now the heterodyne-specific check with no homodyne counterpart, whether the two residuals are independent of each other. Read off the $2\times2$ correlation matrix of the innovation channels:

```wl
Correlation[Transpose[{dWI, dWQ}]] // MatrixForm
```

The diagonal is one and the off-diagonal sits near zero, within the same $\pm 2/\sqrt N$ band: the two residuals are not just each white but mutually independent, exactly the simulation realizing $dZ^2 = 0$. Reading both conjugate quadratures buys two independent white innovations, one per shadow.

Keeping the record splits purity into two quantities that need not agree: the purity of each conditioned state, and the purity of the state you get by averaging them. Run ninety homodyne histories of the same atom:

```wl
homRuns = Table[trajectory[densityMatrix[excited], Hdrive, {cAtom}, {1.}, {}, dtAtom, 6.0, s]["states"], {s, 90}];
```

Follow both purities in time:

```wl
ListLinePlot[{Transpose[{when, Mean[purity /@ #] & /@ Transpose[homRuns]}],
   Transpose[{when, purity[Mean[#]] & /@ Transpose[homRuns]}]},
 PlotStyle -> {ColorData[97, 1], ColorData[97, 2]}, Frame -> True, GridLines -> Automatic,
 PlotLegends -> {"mean purity of conditioned states", "purity of the mean state"},
 PlotRange -> {0, 1.05}, ImageSize -> 520, FrameLabel -> {"time", "purity"},
 PlotLabel -> "every run stays pure; only the record-blind average mixes"]
```

The conditioned purity holds at one for the whole run; the averaged state's purity falls and settles at the mixed steady value. Decoherence sits entirely in the second curve: it is not an event in any watched atom's history but the cost of discarding the record.

Counting, homodyne, heterodyne: three different measurements of one atom, three different trajectory ensembles, one shared master equation. The measurement you choose sets what you learn moment to moment; the unconditioned dynamics is the same. Between homodyne and heterodyne lies the continuous general-dyne family, all on this same average.

## Part Five: What Watching Does, and How to Act on It

Keeping the record does more than reveal the state: measurement reshapes the dynamics, and the record can be fed back to control the system. Five examples follow: a watched-frozen qubit, a charge qubit read by a detector, measurement-induced collapse, a Kalman filter tracking an oscillator, and measurement feedback (all qubits except the oscillator). First a helper to read the $z$ and $x$ components:

```wl
zOf[states_] := blochVector[#][[3]] & /@ states; xOf[states_] := blochVector[#][[1]] & /@ states;
```

### Watching Harder Freezes It

**The problem.** A driven switch continuously measured along $\hat\sigma_z$, the quantum Zeno setup. The conditioned state obeys the stochastic master equation
$$d\rho = -i\,\tfrac{\Omega}{2}[\hat\sigma_x,\rho]\,dt - \tfrac{k}{2}\,\bigl[\hat\sigma_z,[\hat\sigma_z,\rho]\bigr]\,dt + \sqrt{k}\,\bigl(\hat\sigma_z\rho + \rho\hat\sigma_z - 2\langle\hat\sigma_z\rangle\rho\bigr)\,dW,$$
with drive $\Omega$ and measurement strength $k$. When $k \ll \Omega$ the switch swings almost freely; when $k \gg \Omega$ the watching pins it near a pole and $z(t)$ becomes a random telegraph, the "watched pot" that never boils. Watch a gently-watched and a hard-watched history.

A driven qubit measured along $\hat\sigma_z$. Measure weakly ($k \ll \Omega$) and it swings almost freely; measure strongly ($k \gg \Omega$) and it is pinned near a $\hat\sigma_z$ eigenstate with occasional jumps, the quantum Zeno effect. Set the drive and the fine time grid:

```wl
\[CapitalOmega]zeno = 2.0; fastZeno = 0.004; clockZeno = fastZeno Range[0, 3000];
```

Run a weakly and a strongly measured trajectory from the same start and seed:

```wl
gentle = trajectory[densityMatrix[excited], \[CapitalOmega]zeno/2 X, {Sqrt[0.15] Z}, {1.}, {}, fastZeno, 12., 7];
hard = trajectory[densityMatrix[excited], \[CapitalOmega]zeno/2 X, {Sqrt[6.] Z}, {1.}, {}, fastZeno, 12., 7];
```

Compare the magnitude of the time-averaged $\langle\hat\sigma_z\rangle$:

```wl
{Abs@Mean[zOf@gentle["states"]], Abs@Mean[zOf@hard["states"]]}
```

The weakly measured $\langle\hat\sigma_z\rangle$ averages to near zero (it swings freely); the strongly measured one keeps a large mean (pinned).

Now visualize the weakly and strongly measured trajectories together:

```wl
ListLinePlot[{Transpose[{clockZeno, zOf@gentle["states"]}], Transpose[{clockZeno, zOf@hard["states"]}]},
 PlotStyle -> {ColorData[97, 1], Directive[ColorData[97, 2], Opacity[0.8]]}, Frame -> True, GridLines -> Automatic,
 PlotLegends -> {"weak measurement (free)", "strong measurement (pinned)"}, PlotRange -> {-1.05, 1.05},
 ImageSize -> 560, FrameLabel -> {"time", "\[LeftAngleBracket]\[Sigma]z\[RightAngleBracket]"},
 PlotLabel -> "strong measurement freezes the qubit (Zeno)"]
```

The strongly measured qubit clings to an eigenstate for long stretches, then jumps. Strong measurement freezes the dynamics, the Zeno effect: the more often you measure, the less the state can evolve. The jumps are the drive occasionally winning.

The single trace shows the freezing but cannot pin its rate. That rate is exact and sits in the master equation's own spectrum, and following it across the measurement strength tells a two-stage story: watching harder first kills the coherent Rabi oscillation, then freezes the relaxation left behind. Both stages are eigenvalues of the Liouvillian. Build it at a chosen fraction $r = k/\Omega$ of the drive and keep the nonzero ones, so the whole crossover reads against that one ratio and the exceptional point falls at $r = 1$:

```wl
zenoEigs[ratio_] := Select[Eigenvalues[liouvillian[(\[CapitalOmega]zeno/2) X, {Sqrt[ratio \[CapitalOmega]zeno] Z}, 2]], Abs[#] > 10^-9 &];
```

Below the exceptional point, at half the critical strength ($r = 1/2$), the three of them are a complex-conjugate pair and one lone real number:

```wl
zenoEigs[1/2]
```

The pair is $-k \pm i\sqrt{\Omega^2 - k^2}$ (here $k = \Omega/2$), a Rabi oscillation the measurement only damps; the lone real number is the transverse dephasing at $-2k$, off to the side and taking no part in the freezing.

The oscillation lives in the pair's imaginary part and the relaxation in the least-negative real part, so read one frequency and one rate off the spectrum at each ratio:

```wl
zenoFreq[ratio_] := Max[Abs[Im[zenoEigs[ratio]]]];
zenoRate[ratio_] := -Max[Re[zenoEigs[ratio]]];
```

Follow both across the sweep, the oscillation frequency and the slowest relaxation rate against the ratio $k/\Omega$:

```wl
ListLinePlot[
 {Table[{r, zenoFreq[r]}, {r, 0.01, 3., 0.01}], Table[{r, zenoRate[r]}, {r, 0.01, 3., 0.01}]},
 Frame -> True, GridLines -> {{{1, Directive[Gray, Dashed]}}, Automatic},
 PlotLegends -> Placed[{"Rabi oscillation  |Im \[Lambda]|", "relaxation rate  |Re \[Lambda]|"}, {Right, Top}],
 FrameLabel -> {"k/\[CapitalOmega] ratio", "Liouvillian eigenvalue  (1/time)"},
 PlotLabel -> "watching harder: the oscillation dies at k = \[CapitalOmega], the relaxation rate peaks there"]
```

The oscillation frequency slides to zero at $k = \Omega$ and stays there, the coherent swinging gone; at that same strength the relaxation rate reaches its peak and turns over, so watching harder past $\Omega$ slows the decay rather than quickening it. The pair has met on the real axis at an exceptional point, where the eigenvalues and their eigenvectors coincide and the oscillation can no longer be told apart from the decay.

Past that peak the surviving rate falls as one over the strength. Read it across a wide range on log-log, with the $\Omega^2/(2k)$ law drawn on the strong side:

```wl
zenoKs = \[CapitalOmega]zeno Table[2.^e, {e, -3, 4, 0.08}];
ListLogLogPlot[
 {Table[{k, zenoRate[k/\[CapitalOmega]zeno]}, {k, zenoKs}], Table[{k, \[CapitalOmega]zeno^2/(2 k)}, {k, Select[zenoKs, # >= \[CapitalOmega]zeno &]}]},
 Joined -> {False, True},
 PlotLegends -> Placed[{"measured rate  |Re \[Lambda]|", "\!\(\*SuperscriptBox[\(\[CapitalOmega]\), \(2\)]\)/2k  (the Zeno freeze)"}, {Right, Top}],
 Frame -> True, GridLines -> {{{\[CapitalOmega]zeno, Directive[Gray, Dashed]}}, None},
 FrameLabel -> {"measurement strength k", "slowest relaxation rate  |Re \[Lambda]|"},
 PlotLabel -> "the rate follows k up to k = \[CapitalOmega], then freezes as \!\(\*SuperscriptBox[\(\[CapitalOmega]\), \(2\)]\)/2k"]
```

Below $k = \Omega$ the rate climbs along $k$; above it the points settle onto the $\Omega^2/(2k)$ line, the leftover relaxation frozen as $1/k$, the small gap near the peak closing as the watching strengthens. Watching past the optimal strength does not pin the switch harder, it slows what motion is left, the quantum Zeno suppression.

The whole crossover is exact, read straight off the deterministic Liouvillian with no gate, no time step, and no ensemble. Its slowest eigenvalue returns at once in the next section as the width of the central peak the quantum point contact reads off its own current, so the Zeno rate here and that spectrum there are one number seen twice. Coherent Rabi giving way to overdamped freezing through a Liouvillian exceptional point is the modern reading of when watching stops helping and starts freezing ([Snizhko, Kumar, and Romito, 2020](https://arxiv.org/abs/2003.10476)), seen on a driven transmon by [Slichter and Siddiqi, 2016](https://arxiv.org/abs/1512.04006).

### A Charge Qubit Read by a Quantum Point Contact

**The problem.** The last two examples watched an abstract qubit; here is the same watching realized in a solid-state device. A charge qubit (one electron shared between two quantum dots, "dot 1" or "dot 2") is monitored by a nearby quantum point contact (QPC) whose current depends on which dot holds the electron. With interdot tunnelling $\Omega$ and no detuning, the unconditional master equation is
$$\dot\rho = -i\Bigl[\tfrac{\Omega}{2}\hat\sigma_x,\rho\Bigr] + \kappa\,\mathcal{D}[\hat n_1]\rho,\qquad \hat n_1 = \tfrac12(\mathbb{1} - \hat\sigma_z),$$
with $\hat\sigma_z = \hat n_2 - \hat n_1$ the which-dot operator and $\kappa$ the measurement-induced dephasing rate. The QPC current is a continuous weak measurement of $\hat\sigma_z$: when $\kappa \ll \Omega$ the electron tunnels coherently and the current carries an oscillation at the Rabi frequency; when $\kappa \gg \Omega$ the qubit is pinned and the current is a random telegraph. Build the conditioned dynamics, the detector current, and the current's power spectrum.

The dephasing $\hat n_1$ is a projector, and $\mathcal{D}[\hat n_1] = \tfrac14\mathcal{D}[\hat\sigma_z]$ (the identity part drops from any dissipator), so watching the QPC is exactly watching $\hat\sigma_z$ at strength $\kappa/4$: the machinery of the last two examples applies unchanged. Fix the tunnelling and the two regimes:

```wl
OmQpc = 2.0; kSlow = 0.6; kFast = 12.0; qpcLeak[\[Kappa]_] := Sqrt[\[Kappa]/4] Z;
```

The conditioned trajectory uses the toolkit's positivity-preserving step; the recorded current is $J = \sqrt{\kappa}\,\langle\hat\sigma_z\rangle + \text{shot noise}$, the QPC output. Run one weakly-measured (coherent) history and one strongly-measured (Zeno) history:

```wl
dtQpc = 0.004; tfQpc = 10.;
qpcSlow = trajectory[densityMatrix[ground], (OmQpc/2) X, {qpcLeak[kSlow]}, {1.}, {}, dtQpc, tfQpc, 4];
qpcFast = trajectory[densityMatrix[ground], (OmQpc/2) X, {qpcLeak[kFast]}, {1.}, {}, dtQpc, tfQpc, 4];
whenQpc = qpcSlow["times"];
```

The QPC hands back a single real record, $dJ = \sqrt{\kappa}\,\langle\hat\sigma_z\rangle\,dt + dW$, one continuous measurement of $\hat\sigma_z$: it is the homodyne ledger again in a solid-state skin, the record the only data and the coherent which-dot swing inferred from it, buried per step under shot noise many times larger.

Contrast the two regimes directly, the coherent swing against the strongly-watched telegraph:

```wl
ListLinePlot[{Transpose[{whenQpc, zOf@qpcSlow["states"]}], Transpose[{whenQpc, zOf@qpcFast["states"]}]},
 PlotStyle -> {ColorData[97, 1], Directive[ColorData[97, 2], Opacity[0.85]]}, Frame -> True, GridLines -> Automatic,
 PlotLegends -> {"weak (coherent tunnelling)", "strong (Zeno telegraph)"}, PlotRange -> {-1.05, 1.05}, ImageSize -> 560,
 FrameLabel -> {"time", "\[LeftAngleBracket]\[Sigma]z\[RightAngleBracket]"}, PlotLabel -> "measurement strength selects coherent oscillation or telegraph"]
```

Weak measurement leaves the electron tunnelling freely; strong measurement pins it in one dot with rare jumps.

The distinctive QPC signature is not in any single record but in the *spectrum*, because a coherent oscillation is a *frequency*: it cancels in the time-integral of the record and shows itself only in the spectrum, painting a peak at the Rabi frequency in the fluctuations of $\hat\sigma_z$ that the detector current inherits. The power spectrum is the Fourier transform of the steady-state autocorrelation $\langle\delta\hat\sigma_z(\tau)\,\delta\hat\sigma_z(0)\rangle_{\mathrm{ss}}$, where $\delta\hat\sigma_z = \hat\sigma_z - \langle\hat\sigma_z\rangle_{\mathrm{ss}}$ is the fluctuation about the steady mean. Quantum regression is the statement that this correlation runs on the *same* Liouvillian as the state: plant the seed $\delta\hat\sigma_z\,\rho_{\mathrm{ss}}$, propagate it by $e^{\mathcal L\tau}$, and read $\delta\hat\sigma_z$ back off it, so $C(\tau) = \mathrm{Tr}\!\big[\delta\hat\sigma_z\,e^{\mathcal L\tau}(\delta\hat\sigma_z\,\rho_{\mathrm{ss}})\big]$.

Fourier-transforming in $\tau$ replaces the propagator $e^{\mathcal L\tau}$ by a single matrix inverse, the resolvent $(i\omega - \mathcal L)^{-1}$, and diagonalizing the Liouvillian breaks that inverse into one simple pole per eigenvalue $\lambda_j$. An eigenvalue $\lambda_j = -\gamma_j + i\nu_j$ is one decaying, oscillating mode, and its pole is a Lorentzian centered at the frequency $\nu_j = \mathrm{Im}\,\lambda_j$ with half-width the decay rate $\gamma_j = |\mathrm{Re}\,\lambda_j|$. Each Lorentzian is weighted by $c_j\,r_j$: the coefficient $c_j$ is how much of the seed lies along mode $j$, its coordinate when the seed is written in the eigenbasis, and $r_j$ is how strongly the readout $\delta\hat\sigma_z$ sees that mode, its overlap with it. Summing the poles, with the steady mode $\lambda_j = 0$ dropped because $\delta\hat\sigma_z$ has zero steady mean,
$$S(\omega) = 2\,\mathrm{Re}\sum_{\lambda_j\neq 0}\frac{c_j\,r_j}{i\omega - \lambda_j}.$$
This is where the slowest eigenvalue of the Zeno crossover returns: it is the pole nearest the imaginary axis, so it is the narrowest Lorentzian, and its half-width is that same relaxation rate.

Build the general spectrum for any Liouvillian, seed, and readout: eigendecompose the Liouvillian, find the seed's coordinates $c_j$ by solving for them in the eigenbasis, read each mode's overlap $r_j$ with the readout, and sum the poles at every requested frequency:

```wl
regressionSpectrum[big_, seedMat_, readout_, tones_] := Module[{d = Length[readout], vals, vecs, coeffs, resu},
   {vals, vecs} = Eigensystem[Transpose[big]];
   coeffs = LinearSolve[Transpose[vecs], Flatten[seedMat]];
   resu = Table[Tr[readout . ArrayReshape[vecs[[j]], {d, d}]], {j, d^2}];
   Table[2 Re@Total@Table[If[Abs[vals[[j]]] < 10^-9, 0, coeffs[[j]] resu[[j]]/(I w - vals[[j]])], {j, d^2}], {w, tones}]];
```

`vals` are the eigenvalues $\lambda_j$; `coeffs` are the $c_j$, produced by `LinearSolve` writing the seed as a combination of the eigenmodes; `resu` are the $r_j$, each mode's trace against the readout; the last line sums $c_j r_j/(i\omega - \lambda_j)$ over the nonzero eigenvalues.

Now specialize to the QPC. Build its Liouvillian and steady state, center $\hat\sigma_z$ on its steady value to form $\delta\hat\sigma_z$, and hand the symmetrized seed $(\delta\hat\sigma_z\,\rho_{\mathrm{ss}} + \rho_{\mathrm{ss}}\,\delta\hat\sigma_z)/2$ together with that same $\delta\hat\sigma_z$ as the readout to the general function:

```wl
qpcSpectrum[\[Kappa]_, tones_] := With[{big = liouvillian[(OmQpc/2) X, {qpcLeak[\[Kappa]]}, 2], rss = First@steadyState[(OmQpc/2) X, {qpcLeak[\[Kappa]]}]},
   With[{dz = Z - Re[Tr[Z . rss]] id2}, regressionSpectrum[big, (dz . rss + rss . dz)/2, dz, tones]]];
```

Here `dz` is $\delta\hat\sigma_z$, and symmetrizing the seed keeps the correlation real and even in $\tau$, which is the physical spectrum a detector measures.

Confirm the normalization is right: the spectrum integrates to the steady variance of $\hat\sigma_z$, which is one for the maximally mixed steady state:

```wl
With[{wide = Range[-40, 40, 0.02]}, 0.02 Total[qpcSpectrum[1., wide]]/(2 Pi)]
```

It returns one, so the spectral weight is accounted for. Now sweep the measurement strength across the $\hat\sigma_z$ spectrum, from weak through the exceptional point at $\kappa = 4\Omega$ to strong:

```wl
kappaSweep = {1., 4., 6., 8., 12.}; tonesQpc = Range[-6, 6, 0.04];
sweepCols = ColorData["TemperatureMap"] /@ Subdivide[Length[kappaSweep] - 1];
ListLinePlot[Table[Transpose[{tonesQpc, qpcSpectrum[k, tonesQpc]}], {k, kappaSweep}],
 PlotStyle -> (Directive[#, Thickness[0.006]] & /@ sweepCols),
 PlotLegends -> Placed[LineLegend[sweepCols, {"\[Kappa] = 1  (weak)", "\[Kappa] = 4", "\[Kappa] = 6", "\[Kappa] = 8  (EP, \[Kappa] = 4\[CapitalOmega])", "\[Kappa] = 12  (Zeno)"}], {0.8, 0.68}],
 Frame -> True, GridLines -> {{-OmQpc, 0, OmQpc}, None}, PlotRange -> {{-6, 6}, {0, 7}}, ImageSize -> 620,
 FrameLabel -> {"frequency", "\!\(\*SubscriptBox[\(S\), \(zz\)]\)(\[Omega])"},
 PlotLabel -> "watching harder marches the Rabi peaks inward, merging them at the exceptional point"]
```

At weak watching the spectrum is split, twin peaks near $\pm\Omega$, the coherent Rabi oscillation written into the fluctuations of $\hat\sigma_z$, the solid-state echo of the Mollow triplet's sidebands. Watching harder marches those peaks inward to the damped Rabi frequency $\pm\sqrt{\Omega^2 - (\kappa/4)^2}$ and broadens them, until at $\kappa = 4\Omega$ they collide and merge at zero: this is the exceptional point of the Zeno crossover read in frequency, the eigenvalue pair meeting on the real axis seen as two peaks becoming one. Past it a single central Lorentzian remains, narrowing as the Zeno rate $2\Omega^2/\kappa$, its half-width the same slowest eigenvalue that section tracked, so the whole crossover, not just that one rate, is seen twice.

What a lab actually records is the detector *current* spectrum, the flat shot-noise floor with this $S_{zz}(\omega)$ riding on it (scaled by the measurement coupling); a fundamental result, the Korotkov-Averin bound, caps that current peak at four times the shot-noise background for a quantum-limited detector, the solid-state statement of the measurement limit that runs through this whole catalog.

Two different knobs move two different crossovers, and they should not be confused: the *measurement strength* $\kappa$ sets coherent-oscillation versus Zeno-telegraph behaviour (shown above), while the detector *transparency* sets whether the record is read as discrete electron jumps (a counting unravelling like the photon counter) or as this diffusive current. Both unravellings average to the one master equation above.

### Measurement-Induced Localization

**The problem.** The general continuous measurement of an observable (here $\hat\sigma_z$) on an *undriven* switch, written in the catalog's measurement-strength convention:
$$d\rho = k\,\mathcal D[\hat\sigma_z]\rho\,dt + \sqrt{k}\,\bigl(\hat\sigma_z\rho + \rho\hat\sigma_z - 2\langle\hat\sigma_z\rangle\rho\bigr)\,dW,\qquad dy = \langle\hat\sigma_z\rangle\,dt + \frac{dW}{2\sqrt{k}}.$$
With nothing pushing, the state's knowledge sharpens toward an eigenstate of $\hat\sigma_z$ at the rate set by $k$: each history collapses to one pole, the choice fixed only by the noise, while the average over histories stays balanced. Watch many histories make their choice.

Turn off the drive and measure $\hat\sigma_z$ on a state leaning toward $+z$, $\langle\hat\sigma_z\rangle = z_0 = 0.5$. With no Hamiltonian, measurement alone drives each run to a $\hat\sigma_z$ eigenstate, $+1$ or $-1$, and the record that drives the collapse is also what reads out which pole it chose. Fix the biased start, the measured operator, and a fine grid:

```wl
z0Loc = 0.5; x0Loc = Sqrt[1 - z0Loc^2]; measZ = Sqrt[1.5] Z;
dtLoc = 0.002; tfLoc = 2.; gridLoc = dtLoc Range[0, Round[tfLoc/dtLoc]];
```

One run reads the $z$ component throughout; run a fan of them, and beside them solve the master equation for the record-blind mean:

```wl
locZ[seed_] := zOf@trajectory[blochState[x0Loc, 0, z0Loc], 0 id2, {measZ}, {1.}, {}, dtLoc, tfLoc, seed]["states"];
fanLoc = ParallelTable[locZ[s], {s, 60}];
meanLoc = zOf[evolveODE[0 id2, {measZ}, blochState[x0Loc, 0, z0Loc], tfLoc][#] & /@ gridLoc];
```

Now visualize the fan against that conserved mean: every path is absorbed at a pole while the thick master-equation curve holds flat at $z_0$:

```wl
Show[ListLinePlot[Transpose[{gridLoc, #}] & /@ fanLoc, PlotStyle -> Directive[ColorData[97, 1], Opacity[0.12]]],
 ListLinePlot[Transpose[{gridLoc, meanLoc}], PlotStyle -> Directive[Thick, ColorData[97, 2]]],
 Frame -> True, GridLines -> Automatic, PlotRange -> {-1.05, 1.05}, ImageSize -> 560,
 FrameLabel -> {"time", "\[LeftAngleBracket]\[Sigma]z\[RightAngleBracket]"},
 PlotLabel -> "each run is absorbed at a pole; the master-equation mean holds at z0"]
```

Every run ends at a pole, yet the thick curve never leaves $z_0$: the dephasing $\mathcal{D}[\hat\sigma_z]$ leaves the populations untouched, so $\langle\hat\sigma_z\rangle$ is conserved on average. That makes $\langle\hat\sigma_z\rangle(t)$ a bounded martingale trapped in $[-1, 1]$, and since every run is absorbed at $\pm 1$, the split follows from the optional stopping theorem with nothing to fit: $\langle\hat\sigma_z\rangle(\infty)$ averages back to $z_0$, so $p_+ - p_- = z_0$ and $p_+ + p_- = 1$ give $p_\pm = (1 \pm z_0)/2$. Born is a theorem about the martingale, not a fraction the run count has to chase. The symmetric $z_0 = 0$ start, the state the earlier examples called $|+\rangle$, is the equal-weight special case.

The record is what a detector actually delivers, and it both drives the collapse and reads it out. Accumulate each run's record into the integrated readout $Y(t) = \int_0^t dJ$ and keep the pole it lands on:

```wl
readOut[seed_] := With[{run = trajectory[blochState[x0Loc, 0, z0Loc], 0 id2, {measZ}, {1.}, {}, dtLoc, tfLoc, seed]},
   {Prepend[Accumulate[Flatten@run["record"]], 0.], Sign[Last@zOf@run["states"]]}];
readData = ParallelTable[readOut[s], {s, 400}];
```

Now visualize the integrated record, one path per run, colored by the pole the run reached: each single increment is almost pure noise, but the accumulated record fans into two well-separated groups:

```wl
With[{show = readData[[;; 50]]},
 ListLinePlot[Transpose[{gridLoc, #[[1]]}] & /@ show,
  PlotStyle -> (If[#[[2]] > 0, Directive[ColorData[97, 1], Opacity[0.45]], Directive[ColorData[97, 2], Opacity[0.45]]] & /@ show),
  Frame -> True, GridLines -> Automatic, ImageSize -> 560,
  FrameLabel -> {"time", "integrated record  Y(t)"}, PlotLabel -> "the integrated record splits by the pole each run chose"]]
```

The sign of $Y(T)$ is the readout's verdict on the eigenstate; measure how often it is right as the integration window grows:

```wl
fidLoc[T_] := With[{iT = Round[T/dtLoc] + 1}, Mean[Boole[Sign[#[[1, iT]]] == #[[2]]] & /@ readData]];
ListLinePlot[Table[{T, fidLoc[T]}, {T, 0.05, tfLoc, 0.05}], PlotRange -> {0.45, 1.02},
 Frame -> True, GridLines -> Automatic, ImageSize -> 500,
 FrameLabel -> {"integration time T", "assignment fidelity"}, PlotLabel -> "longer integration reads the eigenstate out with certainty"]
```

The fidelity climbs to one: given enough record, the sign of its integral names the eigenstate for sure, a continuous measurement gone projective in the long-time limit. Two things keep this honest rather than a lookup. The eigenstate is not fixed in advance, it is chosen by the very record that reads it out, so the two groups form only as the collapse completes. And the flat-weight integral $Y$ is not the sharpest estimator: the conditional $\langle\hat\sigma_z\rangle(t)$ of the fan above is the Bayesian-optimal filter, which a uniform integral only approximates. This trajectory-read-from-its-own-record is what [Murch, Weber, Macklin, and Siddiqi](https://arxiv.org/abs/1305.7270) tracked on a superconducting qubit, reconstructing the conditioned state from the readout and confirming it by tomography.

### Quantum Kalman Filter: Tracking an Oscillator

**The problem.** Continuously measure the position of a harmonic oscillator; a Gaussian state stays Gaussian and the stochastic master equation collapses to a quantum Kalman-Bucy filter. The conditional covariance $\Sigma$ obeys a *deterministic* matrix Riccati equation, independent of the record,
$$\dot \Sigma = A \Sigma + \Sigma A^{T} + D - \Sigma\,C^{T} C\,\Sigma,$$
with the free-rotation drift, the position readout, and the measurement back-action
$$A = \begin{pmatrix}0 & 1\\-1 & 0\end{pmatrix},\qquad C = \begin{pmatrix}2\sqrt{k} & 0\end{pmatrix},\qquad D = \begin{pmatrix}0 & 0\\0 & k\end{pmatrix}.$$
In the three independent entries $(\Sigma_{xx}, \Sigma_{xp}, \Sigma_{pp})$ of $\Sigma$ this reads
$$\frac{d}{dt}\begin{pmatrix}\Sigma_{xx}\\\Sigma_{xp}\\\Sigma_{pp}\end{pmatrix} = \begin{pmatrix}2\Sigma_{xp} - 4k\,\Sigma_{xx}^{2}\\\Sigma_{pp} - \Sigma_{xx} - 4k\,\Sigma_{xx}\Sigma_{xp}\\k - 2\Sigma_{xp} - 4k\,\Sigma_{xp}^{2}\end{pmatrix},$$
the quadratic $-4k(\cdots)$ terms being the information gain $-\Sigma C^{T} C \Sigma$, which is what makes this a Riccati and not a linear Lyapunov equation. Balanced against the back-action ($+D$), that gain relaxes the uncertainty to a steady value fixed only by how hard you watch, squeezing it below the oscillator's own quietest spread; the running best-guess then follows a linear equation driven by the record. Track the shrinking uncertainty and check its standstill.

Continuously measure the position of a harmonic oscillator and you can track it with a Gaussian estimate: a conditional mean and covariance. The covariance obeys a *deterministic* Riccati equation, independent of the record, set only by the measurement strength. It is the warm particle's covariance law from Part Three, $\dot\Sigma = A\Sigma + \Sigma A^{T} + D$, now carrying the measurement's information gain $-\Sigma C^{T} C\Sigma$: diffusion widened that blob to thermal, watching squeezes this one below the ground state. Set the strength and integrate the Riccati equations from a broad start:

```wl
kKal = 1;
riccati = NDSolveValue[{vx'[t] == 2 vc[t] - 4 kKal vx[t]^2, vc'[t] == vp[t] - vx[t] - 4 kKal vx[t] vc[t],
   vp'[t] == -2 vc[t] + kKal - 4 kKal vc[t]^2, vx[0] == 3., vc[0] == 0., vp[0] == 3.}, {vx, vc, vp}, {t, 0, 6}];
```

A reader for the covariance matrix:

```wl
covariance[t_] := {{riccati[[1]][t], riccati[[2]][t]}, {riccati[[2]][t], riccati[[3]][t]}};
```

Compare the position variance at the start and at the end:

```wl
{riccati[[1]][0], riccati[[1]][6]}
```

The position variance has shrunk from its broad start below the ground-state variance $1/2$: continuous measurement conditionally squeezes it below the zero-point level. Two checks. First, the conditional state is pure, its uncertainty product at the Heisenberg floor:

```wl
Det@covariance[6.]
```

The uncertainty product is $1/4$: the conditional state is pure. Second, the steady covariance is independently the fixed point of the Riccati equation:

```wl
{N[vx /. Solve[{2 vc - 4 kKal vx^2 == 0, vp - vx - 4 kKal vx vc == 0,
      -2 vc + kKal - 4 kKal vc^2 == 0, vx > 0}, {vx, vc, vp}, Reals][[1]]], riccati[[1]][6.]}
```

The fixed point matches the late value.

The conditional state is Gaussian, so picture it as a one-sigma ellipse: the covariance $\Sigma$ around a center $c$ traces $c + \Sigma^{1/2}\{\cos s, \sin s\}$ as $s$ runs the circle. Define it:

```wl
covarianceEllipse[c_, \[CapitalSigma]_] := Table[c + Re[MatrixPower[\[CapitalSigma], 1/2]] . {Cos[s], Sin[s]}, {s, 0, 2 Pi, 0.12}];
```

Now visualize the covariance ellipse over time, broad at first and then tightening and tilting to a small steady one:

```wl
trackTimes = {0, 0.1, 0.25, 0.5, 1, 2, 4};
Legended[Graphics[{Thick, MapIndexed[{ColorData["SunsetColors"][First[#2]/Length[trackTimes]], Line@covarianceEllipse[{0, 0}, covariance[#1]]} &, trackTimes]},
  Frame -> True, GridLines -> Automatic, AspectRatio -> 1, ImageSize -> 400, FrameLabel -> {"position", "momentum"},
  PlotLabel -> "measurement shrinks the covariance to a tight, tilted, steady one"],
 SwatchLegend[ColorData["SunsetColors"][#/Length[trackTimes]] & /@ Range[Length[trackTimes]],
  "t = " <> ToString[#] & /@ trackTimes]]
```

This conditional estimate is exactly what a lab computes in real time to track a trapped particle or membrane, more precisely than its own zero-point spread.

The covariance is only half the filter. It says how sharp the estimate is; the other half is the estimate itself, the conditional *mean*, threading through the noisy record in real time. Watch the filter track. Set a fine grid and draw one record's worth of Wiener increments:

```wl
kalDt = 0.004; kalN = 1500; kalGrid = kalDt Range[0, kalN];
kalDW = BlockRandom[SeedRandom[7]; RandomVariate[NormalDistribution[0, Sqrt[kalDt]], kalN]];
```

A note on the constant $k$. The Riccati above uses the catalog's own measurement-strength convention, $k = k_{\mathrm{cat}} = 2k_{\mathrm{Jacobs}}$, so the informative terms carry $4k$ and the mean gain is $2\sqrt{k}\,V$; the record noise and the gain must then be consistent, which pins the record normalization below. What the code calls a "reference conditional mean" is the physical conditional mean of one simulated record, driven by its own innovation $dW$ with that gain. Roll it out from a displaced start:

```wl
kalTrue = FoldList[Function[{xp, jw}, With[{j = jw[[1]], dw = jw[[2]]},
     {xp[[1]] + xp[[2]] kalDt + 2 Sqrt[kKal] riccati[[1]][kalGrid[[j]]] dw,
      xp[[2]] - xp[[1]] kalDt + 2 Sqrt[kKal] riccati[[2]][kalGrid[[j]]] dw}]],
   {3., 0.}, Transpose[{Range[kalN], kalDW}]];
```

All the observer ever sees is the record, the signal $\langle x\rangle\,dt$ buried in white noise. Consistency with the $4k$-Riccati and the $2\sqrt k$ gain fixes its noise level at $dy = \langle x\rangle\,dt + dW/(2\sqrt k)$:

```wl
kalRec = MapThread[#1[[1]] kalDt + #2/(2 Sqrt[kKal]) &, {Most@kalTrue, kalDW}];
```

Now the filter. Start it from a deliberately *wrong* guess for the position and feed it only the record: at each step it forms its own innovation, the record minus its own prediction, $\widehat{dW} = 2\sqrt{k}\,(dy - \hat x\,dt)$, and corrects the estimate with the Riccati gain,
$$d\hat x = \hat p\,dt + 2\sqrt{k}\,\Sigma_{xx}(t)\,\widehat{dW}, \qquad d\hat p = -\hat x\,dt + 2\sqrt{k}\,\Sigma_{xp}(t)\,\widehat{dW},$$
the deterministic Riccati covariance $\Sigma(t)$ setting the gain:

```wl
kalFilt = FoldList[Function[{xp, jdy}, With[{j = jdy[[1]], dwf = 2 Sqrt[kKal] (jdy[[2]] - xp[[1]] kalDt)},
     {xp[[1]] + xp[[2]] kalDt + 2 Sqrt[kKal] riccati[[1]][kalGrid[[j]]] dwf,
      xp[[2]] - xp[[1]] kalDt + 2 Sqrt[kKal] riccati[[2]][kalGrid[[j]]] dwf}]],
   {-3., 0.}, Transpose[{Range[kalN], kalRec}]];
```

Read the $\pm 1\sigma$ band off the covariance:

```wl
kalBand = Sqrt[riccati[[1]][#]] & /@ kalGrid;
```

Now visualize the tracking: the reference conditional mean, and the estimate from the incorrect prior rising out of its wrong start into a shrinking $\pm 1\sigma$ band:

```wl
Legended[Show[
  ListLinePlot[{Transpose[{kalGrid, kalFilt[[All, 1]] + kalBand}], Transpose[{kalGrid, kalFilt[[All, 1]] - kalBand}]},
   Filling -> {1 -> {2}}, PlotStyle -> Directive[ColorData[97, 1], Opacity[0.15]]],
  ListLinePlot[{Transpose[{kalGrid, kalTrue[[All, 1]]}], Transpose[{kalGrid, kalFilt[[All, 1]]}]},
   PlotStyle -> {Directive[Thick, ColorData[97, 2]], Directive[Thick, ColorData[97, 1]]}],
  Frame -> True, GridLines -> Automatic, PlotRange -> All,
  FrameLabel -> {"time", "position"}, PlotLabel -> "the filter locks onto the reference conditional mean from the record alone"],
 LineLegend[{ColorData[97, 2], ColorData[97, 1]}, {"reference conditional mean", "estimate from wrong prior \[PlusMinus]\[Sigma]"}]]
```

The estimate begins wildly wrong (opposite sign) inside a wide band, then the record pulls it onto the reference conditional mean and the band tightens to the steady conditional width: a wrong prior is forgotten, because the record carries the same information whatever the observer first believed. (The reference mean is not a hidden classical position; it is the conditional mean of the simulated record, the best anyone could know.) Confirm the lock, comparing the two once the filter has settled:

```wl
{Last@kalTrue[[All, 1]], Last@kalFilt[[All, 1]]}
```

The two agree to several digits.

The engine of the filter is its own innovation, the record minus the *filter's* prediction (not the reference mean's). For an optimal filter this must be white noise, zero-mean and uncorrelated; measure its autocorrelation past the transient, with the $\pm 2/\sqrt N$ band that white noise would sit inside:

```wl
kalInnov = MapThread[2 Sqrt[kKal] (#2 - #1[[1]] kalDt) &, {Most@kalFilt, kalRec}][[Round[kalN/3] ;;]]/Sqrt[kalDt];
With[{ci = 2/Sqrt[Length[kalInnov]]},
 ListPlot[CorrelationFunction[kalInnov, {1, 10}], Filling -> Axis, PlotRange -> {-0.2, 0.2},
  Frame -> True, GridLines -> {None, {{ci, Directive[Gray, Dashed]}, {-ci, Directive[Gray, Dashed]}}},
  ImageSize -> 440, PlotMarkers -> Automatic,
  FrameLabel -> {"lag", "innovation autocorrelation"}, PlotLabel -> "the filter's innovation is white within confidence limits"]]
```

The autocorrelation is flat at zero for every nonzero lag: the filter has wrung all the structure out of the record, leaving pure noise. A whitened innovation is the working definition of an optimal filter, and it is the real-time diagnostic a lab watches to know its tracking model is right.

### Measurement Feedback: Steering With the Record

**The problem.** Feed the homodyne record back as a Hamiltonian $\hat H_{\mathrm{fb}} = J(t)\,\hat F$; averaged over the noise, the conditioned-plus-feedback evolution becomes the Wiseman-Milburn feedback master equation
$$\dot\rho = -i[\hat H,\rho] + \mathcal{D}[\hat c]\rho - i[\hat F,\ \hat c\rho + \rho\hat c^\dagger] + \tfrac{1}{\eta}\,\mathcal{D}[\hat F]\rho,$$
with $\eta$ the detection efficiency, this is just a plain leak with a shifted operator, $\mathcal{D}[\hat c - i\hat F]\rho$, plus a Hamiltonian correction. The feedback engineers the effective damping, here steering a shapeless switch onto a chosen state. Drive $\langle\hat\sigma_x\rangle$ to one with the record itself and check it against this equation.

Close the loop: feed the measurement record back as a Hamiltonian in proportion, so the record steers the state. Concretely the measured operator is $\hat c = \sqrt{G}\,\hat\sigma_z$ and the feedback generator is $\hat F = \sqrt{G}\,\hat\sigma_y$, so the loop measures $\hat\sigma_z$ and rotates about $Y$ by an amount set by the record. The shifted operator is then $\hat c - i\hat F = \sqrt{G}\,(\hat\sigma_z - i\hat\sigma_y)$, and the Hamiltonian correction $\tfrac12(\hat c^\dagger\hat F + \hat F\hat c)$ vanishes because $\hat\sigma_z$ and $\hat\sigma_y$ anticommute, so the ideal loop is the single leak $G\,\mathcal{D}[\hat\sigma_z - i\hat\sigma_y]$. Set the gain and the time grid:

```wl
Gfb = 0.8; tickFb = 0.01; track = tickFb Range[0, 200];
```

A rotation about $Y$ has a closed form, $R_y(\phi) = e^{-i\phi\hat\sigma_y} = \cos\phi\,\mathbb{1} - i\sin\phi\,\hat\sigma_y$, and it is what the generator $\hat F$ produces: the feedback unitary applied each step is $e^{-i\hat F\,dy} = R_y(\sqrt{G}\,dy)$, so $\hat F$ is the generator and $R_y$ the finite rotation, not one and the same. The step is light:

```wl
rotateY[\[Phi]_] := Cos[\[Phi]] id2 - I Sin[\[Phi]] Y;
```

One steered run: measure, then rotate by an amount set by the record just read:

```wl
steer[seed_] := BlockRandom[SeedRandom[seed];
   With[{measure = measurementStep[0 id2, {Sqrt[Gfb] Z}, {1.}, {}, tickFb]},
    FoldList[Function[{rho, dw}, With[{dy = measurementRecord[rho, {Sqrt[Gfb] Z}, {1.}, tickFb, {dw}]},
       With[{u = rotateY[Sqrt[Gfb] dy[[1]]]}, u . measure[rho, dy] . ConjugateTranspose[u]]]],
     id2/2, RandomVariate[NormalDistribution[0, Sqrt[tickFb]], 200]]]];
```

This feedback has a target: the $+X$ state is a fixed point the loop leaves alone, so from the maximally mixed state the loop should build it. Confirm $+X$ is the dark state (the shifted operator annihilates it):

```wl
Chop[(Sqrt[Gfb] (Z - I Y)) . plus]
```

The zero vector confirms it. Integrate the feedback master equation for the predicted $\langle X\rangle$:

```wl
predicted = blochVector[evolveODE[0 id2, {Sqrt[Gfb] (Z - I Y)}, id2/2, 2.][#]][[1]] & /@ track;
```

Average sixty steered runs:

```wl
steeredMean = Mean[Table[xOf@steer[s], {s, 60}]];
```

Now visualize one run, the ensemble mean, and the master-equation prediction climbing together:

```wl
ListLinePlot[{Transpose[{track, xOf@steer[2]}], Transpose[{track, steeredMean}], Transpose[{track, predicted}]},
 PlotStyle -> {Directive[ColorData[97, 1], Opacity[0.45]], Directive[Thick, ColorData[97, 2]], Directive[Dashed, Red]},
 PlotLegends -> {"one steered run", "average of many", "predicted"}, Frame -> True, GridLines -> Automatic,
 PlotRange -> {-0.35, 1.05}, ImageSize -> 520, FrameLabel -> {"time", "\[LeftAngleBracket]X\[RightAngleBracket]"},
 PlotLabel -> "feedback drives the maximally mixed state to a pure state"]
```

$\langle X\rangle$ climbs from zero to one: the maximally mixed state is driven to a pure state by its own record. This is measurement-based feedback control, the basis of feedback cooling and qubit stabilization.

The loop has a causal chain worth seeing whole: the noisy record drives a control command, which drives the state's response. Rerun one steered history, keeping the rotation the record commands at each step alongside the state:

```wl
steerRec[seed_] := BlockRandom[SeedRandom[seed]; With[{measure = measurementStep[0 id2, {Sqrt[Gfb] Z}, {1.}, {}, tickFb]},
   Rest@FoldList[Function[{prev, dw}, With[{rho = prev[[1]]}, With[{dy = measurementRecord[rho, {Sqrt[Gfb] Z}, {1.}, tickFb, {dw}]},
       With[{u = rotateY[Sqrt[Gfb] dy[[1]]]}, {u . measure[rho, dy] . ConjugateTranspose[u], Sqrt[Gfb] dy[[1]]}]]]],
     {id2/2, 0.}, RandomVariate[NormalDistribution[0, Sqrt[tickFb]], 200]]]];
feedRun = steerRec[2];
```

Now visualize the loop in time, two panels side by side: the feedback rotation the raw record commands at each step, a noise band, and the state's response $\langle\hat\sigma_x\rangle$ climbing to the target:

```wl
Row[{
  ListLinePlot[Transpose[{track[[2 ;;]], feedRun[[All, 2]]}], Frame -> True, GridLines -> Automatic,
   ImageSize -> Medium, PlotRange -> All, FrameLabel -> {"time", "feedback rotation"},
   PlotLabel -> "the command the raw record drives"], "  ",
  ListLinePlot[Transpose[{track[[2 ;;]], blochVector[#[[1]]][[1]] & /@ feedRun}], Frame -> True, GridLines -> Automatic,
   ImageSize -> Medium, PlotRange -> All, FrameLabel -> {"time", "\[LeftAngleBracket]X\[RightAngleBracket]"},
   PlotLabel -> "the response it produces"]}]
```

Command to response, left to right: the loop reads the noisy record as a rotation, a noise band, and the rotations walk the state smoothly to the target.

But a real loop never has a perfect detector, and inefficiency caps how pure the steered state can get. The Wiseman-Milburn feedback master equation carries an efficiency $\eta$ that adds an irreducible feedback-noise term $\tfrac{1-\eta}{\eta}\mathcal{D}[\hat F]$; solve its steady state against $\eta$:

```wl
steadyXeff[\[Eta]_] := blochVector[First@steadyState[0 id2, {Sqrt[Gfb] (Z - I Y), Sqrt[(1 - \[Eta])/\[Eta]] Sqrt[Gfb] Y}]][[1]];
etasFb = Range[0.3, 1, 0.05];
ListLinePlot[Transpose[{etasFb, steadyXeff /@ etasFb}], PlotMarkers -> Automatic, PlotStyle -> ColorData[97, 2], Frame -> True,
 GridLines -> Automatic, PlotRange -> {0, 1.05}, ImageSize -> 460, FrameLabel -> {"detection efficiency \[Eta]", "steady \[LeftAngleBracket]X\[RightAngleBracket]"},
 PlotLabel -> "a leaky detector cannot fully purify: feedback noise sets the floor"]
```

Steady $\langle X\rangle$ climbs to one only as $\eta \to 1$; a lossy detector leaves the loop feeding its own missing information back as noise, and the target is reached only in part. Perfect stabilization needs a perfect record, the same efficiency limit that capped the feedback cooling.

## Part Six: Settling Warm, Purifying Fast, Cooling Cold

### Thermalization to a Bath

**The problem.** A warm bath both takes and gives energy, settling a system at its own temperature. For an oscillator,
$$\dot\rho = -i[\omega\hat a^\dagger\hat a,\rho] + \gamma(n_T+1)\,\mathcal{D}[\hat a]\rho + \gamma\,n_T\,\mathcal{D}[\hat a^\dagger]\rho,\qquad n_T = \frac{1}{e^{\beta\omega} - 1},$$
where the down-rate $\gamma(n_T+1)$ beats the up-rate $\gamma n_T$ by exactly the Boltzmann factor $e^{\beta\omega}$. This *detailed balance* forces the steady state to the thermal distribution $\rho_{\mathrm{ss}}\propto e^{-\beta\omega\,\hat a^\dagger\hat a}$ with mean occupation $n_T$. Settle a cold and a hot start and confirm both land on the thermal distribution.

A warm bath both absorbs and emits, settling a system at its temperature. Model it as two dissipators: decay $\hat a$ at rate $\gamma(n_T+1)$ and excitation $\hat a^\dagger$ at rate $\gamma n_T$. Set the cutoff, the ladder, the bath occupation and its rate:

```wl
topTh = 24; a = annihilation[topTh]; count = creation[topTh] . a; nT = 1.2; \[Gamma]th = 0.5;
```

The two channels, decay and excitation, with their thermal weights:

```wl
fall = Sqrt[\[Gamma]th (nT + 1)] a; climb = Sqrt[\[Gamma]th nT] creation[topTh];
```

A reader for the mean occupation:

```wl
units[rho_] := expectation[count, rho];
```

Settle a cold (vacuum) and a hot start:

```wl
fromCold = evolveODE[0 IdentityMatrix[topTh], {fall, climb}, densityMatrix[coherentState[topTh, 0]], 16.];
fromHot = evolveODE[0 IdentityMatrix[topTh], {fall, climb}, densityMatrix[UnitVector[topTh, 9]], 16.];
```

Read both late occupations against $n_T$:

```wl
{units[fromCold[16.]], units[fromHot[16.]], nT}
```

Both reach the same mean occupation $n_T$.

Not just the mean: the full population distribution settles to the thermal (Boltzmann) shape, the geometric distribution $p_n = (1 - r)\,r^n$ with ratio $r = \frac{n_T}{n_T+1} = e^{-\beta\omega}$, each Fock level less likely than the one below by that factor. Read the settled populations and the thermal shape they should match:

```wl
pops = Re@Diagonal[fromCold[16.]][[;; 8]];
warmShape = With[{ratio = nT/(nT + 1)}, (1 - ratio) ratio^Range[0, 7]];
```

Measure the largest mismatch:

```wl
Max@Abs[pops - warmShape]
```

The settled populations match the geometric thermal distribution to a fraction of a percent. The reason is **detailed balance**: the excitation rate divided by the decay rate is exactly the Boltzmann ratio that sets the distribution:

```wl
{(\[Gamma]th nT)/(\[Gamma]th (nT + 1)), nT/(nT + 1)}
```

Now visualize both thermalization checks: the two means converging to the thermal occupation and the settled populations matching the thermal distribution:

```wl
GraphicsRow[{
  ListLinePlot[{Transpose[{Range[0, 16, 0.05], units[fromCold[#]] & /@ Range[0, 16, 0.05]}],
     Transpose[{Range[0, 16, 0.05], units[fromHot[#]] & /@ Range[0, 16, 0.05]}], {{0, nT}, {16, nT}}},
    PlotStyle -> {ColorData[97, 1], ColorData[97, 2], Directive[Dashed, Gray]},
    PlotLegends -> {"from cold", "from hot", "warm equilibrium"}, Frame -> True, GridLines -> Automatic,
    ImageSize -> 360, FrameLabel -> {"time", "mean occupation"}, PlotLabel -> "both starts reach the thermal mean"],
  BarChart[{pops, warmShape}, ChartLegends -> {"settled", "thermal"}, ChartLabels -> {Range[0, 7], None},
    Frame -> True, ImageSize -> 340, FrameLabel -> {"Fock level", "probability"}, PlotLabel -> "the settled distribution is thermal"]},
 ImageSize -> 720]
```

Detailed balance has a sharper signature than the bar chart shows: on a logarithmic axis the geometric thermal distribution is a straight line, and its constant slope is the Boltzmann ratio. Plot the settled and thermal populations on a semilog scale, and the successive ratio $p_{n+1}/p_n$ beside them:

```wl
GraphicsRow[{
  ListLogPlot[{pops, warmShape}, Joined -> True, PlotMarkers -> Automatic, PlotStyle -> {ColorData[97, 1], Directive[Dashed, Gray]},
   PlotLegends -> {"settled", "thermal"}, Frame -> True, GridLines -> Automatic, ImageSize -> 360,
   FrameLabel -> {"Fock level n", "population (log)"}, PlotLabel -> "geometric tail: a straight line on a log axis"],
  ListLinePlot[Transpose[{Range[0, 6], pops[[2 ;;]]/pops[[;; 7]]}], PlotMarkers -> Automatic, PlotRange -> {0, 1},
   Frame -> True, GridLines -> Automatic, ImageSize -> 360, PlotStyle -> ColorData[97, 1],
   Epilog -> {Directive[Dashed, Gray], Line[{{0, nT/(nT + 1)}, {6, nT/(nT + 1)}}]},
   FrameLabel -> {"Fock level n", "\!\(\*SubscriptBox[\(p\), \(n + 1\)]\)/\!\(\*SubscriptBox[\(p\), \(n\)]\)"},
   PlotLabel -> "the ratio is flat at the Boltzmann factor"]},
 ImageSize -> 720]
```

The settled populations fall along a straight line parallel to the thermal one, and their level-to-level ratio sits flat on $n_T/(n_T+1) = e^{-\beta\omega}$: every rung is less likely than the one below by exactly the Boltzmann factor, the visible fingerprint of detailed balance. Damping and thermalization are the same bath: it drains a hot system and feeds a cold one, and the balance of decay and excitation is thermal equilibrium. This is relaxation of an oscillator to its bath temperature, derived from the two rates.

### Rapid Purification by Feedback

**The problem.** How fast a measurement purifies a qubit depends on the angle $\theta$ between the measured axis and the Bloch vector. Measuring $\hat M = \sin\theta\,\hat\sigma_x + \cos\theta\,\hat\sigma_z$, the impurity $L = \tfrac12(1 - |\vec a|^2)$, with $\vec a = (x,y,z)$ the Bloch vector, falls on average as
$$\overline{dL} = -4k\,L\,\bigl(\sin^2\theta + 2L\cos^2\theta\bigr)\,dt$$
(here $\theta$ is measured from the Bloch vector and $k = k_{\mathrm{cat}} = 2k_{\mathrm{Jacobs}}$, the catalog's convention, so the crosswise rate is $4k$; in Jacobs' own $\theta$-from-the-pole and $k_{\mathrm{Jacobs}}$ it reads $-8k_{\mathrm{Jacobs}}L[1-(1-2L)\cos^2\theta]$), fastest when the measurement is held perpendicular to the Bloch vector ($\theta = \pi/2$). Feeding back a rotation that keeps it perpendicular makes the purification deterministic, $L(t) = e^{-4kt}L(0)$, while a fixed measurement, left aligned with the state it purifies, has no single rate.

One caveat sets what "optimal" means: this steering maximizes the *mean purity at a fixed time*, but it actually lengthens the *mean time* to reach a chosen purity, so the best protocol depends on which you care about (Wiseman-Ralph). Compare a fixed and a steered measurement.

Measurement purifies a qubit (Part Five), but the rate depends on the angle between the measured axis and the Bloch vector: fastest when the measurement is crosswise to the state, slowest when aligned. A fixed measurement drives the state along its own axis, where it learns slowest. Feeding back a rotation that keeps the measurement crosswise purifies fastest. Compare fixed and steered, following the mean linear entropy from the maximally mixed state. The impurity reader:

```wl
linearEntropy[rho_] := (1 - Norm[blochVector[rho]]^2)/2;
```

Set the measurement strength and the time grid:

```wl
kPur = 1.; tickPur = 0.01; span = tickPur Range[0, 120];
```

One fixed-axis run, measuring $Z$ throughout:

```wl
fixedBlur[seed_] := linearEntropy /@ trajectory[id2/2, 0 id2, {Sqrt[kPur] Z}, {1.}, {}, tickPur, 1.2, seed]["states"];
```

The crosswise direction to a given Bloch vector:

```wl
crosswise[lean_] := With[{c = Cross[lean, {0, 0, 1.}]}, If[Norm[c] < 1.*^-6, {1., 0., 0.}, Normalize[c]]];
```

One steered step: measure along the axis crosswise to the current state:

```wl
steerStep[rho_, dw_] := With[{L = Sqrt[kPur] (crosswise[blochVector[rho]] . {X, Y, Z})},
   measurementStep[0 id2, {L}, {1.}, {}, tickPur][rho, measurementRecord[rho, {L}, {1.}, tickPur, {dw}]]];
```

One steered run:

```wl
steeredBlur[seed_] := linearEntropy /@ BlockRandom[SeedRandom[seed];
    FoldList[steerStep, id2/2, RandomVariate[NormalDistribution[0, Sqrt[tickPur]], 120]]];
```

Average eighty runs of each strategy:

```wl
fixedAvg = Mean[fixedBlur /@ Range[80]]; steeredAvg = Mean[steeredBlur /@ Range[80]];
```

Compare the two at the same instant:

```wl
{fixedAvg[[61]], steeredAvg[[61]]}
```

At the same time the steered strategy has reached a lower linear entropy.

Now visualize both purification strategies on a log scale, where a straight line represents exponential decay:

```wl
ListLogPlot[{Transpose[{span, fixedAvg}], Transpose[{span, steeredAvg}]}, Joined -> True,
 PlotStyle -> {ColorData[97, 1], ColorData[97, 2]}, Frame -> True, GridLines -> Automatic,
 PlotLegends -> {"fixed measurement", "steered measurement"}, ImageSize -> 480,
 FrameLabel -> {"time", "linear entropy (log scale)"}, PlotLabel -> "feedback purifies faster, and deterministically"]
```

The steered curve is steeper and smoother (a straight line tracking its deterministic exponential, while the fixed one is ragged, with no single rate): feedback makes the mean purification faster and deterministic.

But the mean hides the caveat, so look at the distributions. Run longer ensembles of each strategy so nearly every run reaches a target purity:

```wl
purGrid = tickPur Range[0, 200]; purN = 200;
fixedBlurL[seed_] := linearEntropy /@ trajectory[id2/2, 0 id2, {Sqrt[kPur] Z}, {1.}, {}, tickPur, 2., seed]["states"];
steeredBlurL[seed_] := linearEntropy /@ BlockRandom[SeedRandom[seed]; FoldList[steerStep, id2/2, RandomVariate[NormalDistribution[0, Sqrt[tickPur]], 200]]];
fixedRunsPur = fixedBlurL /@ Range[purN]; steeredRunsPur = steeredBlurL /@ Range[purN];
```

First, that "deterministic" claim: the steered runs should each hug the exact law $L(t) = e^{-4kt}L(0)$ with $L(0) = 1/2$, while a fixed run wanders. Overlay thirty steered paths on the analytic exponential, on a log scale:

```wl
Show[ListLogPlot[Transpose[{purGrid, #}] & /@ steeredRunsPur[[;; 30]], Joined -> True, PlotStyle -> Directive[ColorData[97, 2], Opacity[0.2]]],
 LogPlot[0.5 Exp[-4 kPur t], {t, 0, 2}, PlotStyle -> Directive[Thick, Black, Dashed]],
 Frame -> True, GridLines -> Automatic, ImageSize -> 480, PlotRange -> {{0, 2}, {10^-4, 1}},
 FrameLabel -> {"time", "linear entropy (log)"}, PlotLabel -> "steered runs hug the deterministic exponential (dashed)"]
```

Each steered path clings to the dashed line: the feedback really does make the purification deterministic run by run, not just on average.

Now the objective-dependence. At a *fixed time* the steered strategy wins, its impurity distribution sitting lower and tighter than the fixed one:

```wl
Histogram[{fixedRunsPur[[All, 61]], steeredRunsPur[[All, 61]]}, {0, 0.35, 0.02}, "PDF",
 ChartLegends -> {"fixed", "steered"}, ChartStyle -> {ColorData[97, 1], ColorData[97, 2]}, Frame -> True, ImageSize -> 440,
 FrameLabel -> {"linear entropy at t = 0.6", "density"}, PlotLabel -> "at fixed time, steered is lower and tighter"]
```

But the *first-passage time* to a chosen purity tells the opposite story. Record when each run first drops below $L = 0.05$:

```wl
fptPur[run_] := With[{pos = FirstPosition[run, l_ /; l < 0.05]}, If[MissingQ[pos], 2., (First[pos] - 1) tickPur]];
{fptFixed, fptSteer} = {fptPur /@ fixedRunsPur, fptPur /@ steeredRunsPur};
Histogram[{fptFixed, fptSteer}, {0, 1.2, 0.06}, "PDF", ChartLegends -> {"fixed", "steered"},
 ChartStyle -> {ColorData[97, 1], ColorData[97, 2]}, Frame -> True, ImageSize -> 440,
 FrameLabel -> {"first-passage time to L < 0.05", "density"}, PlotLabel -> "yet the fixed measurement often reaches the target sooner"]
```

Compare the two means directly, fixed-time impurity against mean hitting time:

```wl
{{"fixed-time L", Mean[fixedRunsPur[[All, 61]]], Mean[steeredRunsPur[[All, 61]]]},
 {"mean first-passage", Mean[fptFixed], Mean[fptSteer]}}
```

The steered strategy has the lower fixed-time impurity but the *longer* mean first-passage: its determinism forbids the lucky fast records that let a fixed measurement occasionally reach the target early. Which protocol is "optimal" depends on the question, sharp average purity by a deadline or shortest typical time to a goal, and no single curve can answer both. This is the Wiseman-Ralph point, made by the distributions rather than asserted.

### Feedback Cooling of a Mechanical Oscillator

**The problem.** Continuously measure a mechanical oscillator's position (its step-down is $\hat b$) through a cavity and push against its estimated motion to cool it, against a thermal bath:
$$d\rho = -i[\omega\hat b^\dagger\hat b - f(t)\hat x,\ \rho]\,dt + k\,\mathcal D[\hat x]\rho\,dt + \sqrt{\eta k}\,\bigl(\hat x\rho + \rho\hat x - 2\langle\hat x\rangle\rho\bigr)\,dW + \gamma(n_T+1)\mathcal{D}[\hat b]\rho\,dt + \gamma n_T\mathcal{D}[\hat b^\dagger]\rho\,dt.$$
Watching sharpens the estimate, the feedback force $f \propto -G\langle\hat p\rangle$ drains the energy, and the leftover is a balance of measurement rate, efficiency $\eta$, gain $G$, and heating $\gamma n_T$. Watching alone only heats; feedback cools, but only up to an optimal gain beyond which it feeds the detector's noise back onto the oscillator. Cool a warm start, find the best gain, and sweep the detector quality.

The same measure-and-feedback loop cools a mechanical oscillator, but cooling takes more than watching: the loop must act. Three pieces enter now that the Kalman tracker did not have. A thermal bath (damping $\gamma$, occupation $n_T$) keeps trying to warm the oscillator; the continuous position measurement (strength $k$, efficiency $\eta$) tracks it; and a feedback force pushes back against the *estimated* velocity, $f = -G\,\langle\hat p\rangle_c$, draining energy. The state is Gaussian throughout, so, as in the Kalman filter, the conditional covariance is deterministic and the conditional mean is stochastic. This example stands alone, so fix its own constants:

```wl
kCool = 1.; etaCool = 1.; gammaCool = 0.1; nCool = 10.;
```

The conditional covariance is the tracking Riccati again, now carrying the bath (a drift toward the thermal spread $n_T + 1/2$) and the efficiency $\eta$ on the information-gain terms while the back-action term $+k$ stays full:

```wl
riccatiCool[\[Eta]_] := NDSolveValue[{
    vx'[t] == 2 vc[t] - 4 \[Eta] kCool vx[t]^2 - gammaCool vx[t] + gammaCool (nCool + 1/2),
    vc'[t] == vp[t] - vx[t] - 4 \[Eta] kCool vx[t] vc[t] - gammaCool vc[t],
    vp'[t] == -2 vc[t] + kCool - 4 \[Eta] kCool vc[t]^2 - gammaCool vp[t] + gammaCool (nCool + 1/2),
    vx[0] == nCool + 1/2, vc[0] == 0., vp[0] == nCool + 1/2}, {vx, vc, vp}, {t, 0, 80}];
```

Its steady value is the sharpest the observer can ever know the oscillator, the conditional floor, and its Kalman gain $2\sqrt{\eta k}\,V$ drives the estimate:

```wl
steadyCovCool[\[Eta]_] := Through[riccatiCool[\[Eta]][80.]];
{vxCool, vcCool, vpCool} = steadyCovCool[1.];
condFloorCool = (vxCool + vpCool)/2 - 0.5;
kxCool = 2 Sqrt[kCool] vxCool; kpCool = 2 Sqrt[kCool] vcCool;
```

Here is the key point, and the reason watching alone is not cooling: averaged over all records, the unconditional occupation is the conditional floor plus the spread of the estimate about it, and that spread solves a Lyapunov equation with the closed-loop drift $A_G$ (free rotation, bath damping, and the feedback $-G$ on $\hat p$) fed by the Kalman-gain noise:

```wl
nUncCool[gg_] := With[{s = LyapunovSolve[{{-gammaCool/2, 1}, {-1, -gammaCool/2 - gg}},
     -{{kxCool^2, kxCool kpCool}, {kxCool kpCool, kpCool^2}}]},
   (s[[1, 1]] + vxCool + s[[2, 2]] + vpCool)/2 - 0.5];
```

With no feedback, $G = 0$, watching does not cool at all: it *heats*, because measurement back-action kicks the momentum and nothing removes the energy. Confirm this two ways, the conditional-plus-spread bookkeeping against the plain unconditional master equation (back-action on $\Sigma_{pp}$, bath, no information gain):

```wl
{nUncCool[0.],
 With[{u = NDSolveValue[{vx'[t] == 2 vc[t] - gammaCool vx[t] + gammaCool (nCool + 1/2),
      vc'[t] == vp[t] - vx[t] - gammaCool vc[t], vp'[t] == -2 vc[t] + kCool - gammaCool vp[t] + gammaCool (nCool + 1/2),
      vx[0] == nCool + 0.5, vc[0] == 0., vp[0] == nCool + 0.5}, {vx, vp}, {t, 0, 200}]},
  (u[[1]][200.] + u[[2]][200.])/2 - 0.5]}
```

Both give the same number, above the bath's $n_T = 10$: the meter alone leaves the oscillator hotter than it found it.

Now close the loop. The conditional mean is a stochastic trajectory driven by the innovation, with the feedback force damping the momentum estimate:

```wl
Gcool = 2.; dtCool = 0.01; nStepCool = 2500; timesCool = dtCool Range[0, nStepCool];
```

At $t=0$ the oscillator is thermal, so the total covariance is $(n_T + 1/2)\mathbb{1}$. That total splits into the conditional floor the filter always carries plus the spread of the estimate about it, so the means must be drawn with the full *difference* $(n_T+\tfrac12)\mathbb{1} - V_c$, including its off-diagonal $-v_c$, to avoid double-counting the initial uncertainty and to make the initial unconditional state exactly thermal:

```wl
meanRunCool[gg_, seed_] := BlockRandom[SeedRandom[seed];
   Module[{init = RandomVariate[MultinormalDistribution[{0., 0.},
        {{nCool + 0.5 - vxCool, -vcCool}, {-vcCool, nCool + 0.5 - vpCool}}]],
     kicks = RandomVariate[NormalDistribution[0, Sqrt[dtCool]], nStepCool]},
    FoldList[Function[{xp, dw},
      {xp[[1]] + (xp[[2]] - gammaCool/2 xp[[1]]) dtCool + kxCool dw,
       xp[[2]] + (-xp[[1]] - gammaCool/2 xp[[2]] - gg xp[[2]]) dtCool + kpCool dw}],
     init, kicks]]];
```

The unconditional occupation in time is the ensemble spread of these means plus the conditional floor. Track it with the loop open and closed:

```wl
nEnsCool = 120;
nDrainCool[gg_] := Module[{ens = Table[meanRunCool[gg, s], {s, nEnsCool}]},
   Table[(Variance[ens[[All, j, 1]]] + vxCool + Variance[ens[[All, j, 2]]] + vpCool)/2 - 0.5,
    {j, 1, nStepCool + 1, 25}]];
{coolOff, coolOn} = {nDrainCool[0.], nDrainCool[Gcool]};
```

Now visualize the cooling: the occupation with the loop off (heating above the bath) and on (draining toward, but not reaching, the conditional floor):

```wl
ListLinePlot[{Transpose[{timesCool[[1 ;; ;; 25]], coolOff}], Transpose[{timesCool[[1 ;; ;; 25]], coolOn}]},
 PlotStyle -> {ColorData[97, 2], ColorData[97, 3]}, Frame -> True, GridLines -> Automatic, PlotRange -> {0, All},
 PlotLegends -> {"loop open (watching only)", "loop closed (feedback)"}, ImageSize -> 520,
 FrameLabel -> {"time", "phonon number"}, PlotLabel -> "feedback cools; watching alone heats",
 Epilog -> {Directive[Gray, Dashed], Line[{{0, nCool}, {Last@timesCool, nCool}}],
   Directive[Black, Dashed], Line[{{0, condFloorCool}, {Last@timesCool, condFloorCool}}]}]
```

The open loop climbs above the thermal line (gray); the closed loop drains fast and settles well below it, above the conditional floor (black).

To see the loop acting, take one record and plot the estimated position with the feedback force it drives:

```wl
coolRec = meanRunCool[Gcool, 4];
Column[{
  ListLinePlot[Transpose[{timesCool, coolRec[[All, 1]]}], PlotStyle -> ColorData[97, 1], Frame -> True,
   GridLines -> Automatic, ImageSize -> 520, PlotRange -> {{0, 25}, All}, FrameLabel -> {"time", "estimated position"},
   PlotLabel -> "the tracked motion rings down"],
  ListLinePlot[Transpose[{timesCool, -Gcool coolRec[[All, 2]]}], PlotStyle -> ColorData[97, 3], Frame -> True,
   GridLines -> Automatic, ImageSize -> 520, PlotRange -> {{0, 25}, All}, FrameLabel -> {"time", "feedback force f"},
   PlotLabel -> "the force the loop applies"]}]
```

The estimate rings down as the force pushes against it, the mechanical analogue of a shock absorber built from a measurement.

How hard should the loop push? Sweep the gain and read the steady occupation off the exact Lyapunov balance:

```wl
gainsCool = Range[0, 12, 0.4];
ListLinePlot[Transpose[{gainsCool, nUncCool /@ gainsCool}], PlotStyle -> ColorData[97, 3], Frame -> True,
 GridLines -> Automatic, ImageSize -> 500, PlotRange -> {0, All}, FrameLabel -> {"feedback gain G", "steady phonon number"},
 PlotLabel -> "an optimal gain: too little leaves it warm, too much feeds noise back",
 Epilog -> {Directive[Gray, Dashed], Line[{{0, nCool}, {Last@gainsCool, nCool}}],
   Directive[Black, Dashed], Line[{{0, condFloorCool}, {Last@gainsCool, condFloorCool}}]}]
```

The occupation falls steeply from the heated open-loop value, reaches a minimum near $G \approx 2$, then rises again: beyond the optimum the loop feeds the estimate's own measurement noise back onto the oscillator, the imprecision-back-action tradeoff that caps every real feedback-cooling experiment. It never reaches the conditional floor (black): that would need a noiseless actuator acting on a perfectly known state.

Finally, a better detector lowers both the floor and the achievable cold. Sweep the efficiency at the gain optimized for the ideal detector (a fixed gain, not re-optimized per efficiency, which is why the cooled curve is a conservative estimate at low $\eta$):

```wl
nUncCoolEta[\[Eta]_, gg_] := Module[{c = steadyCovCool[\[Eta]], kx, kp, s},
   kx = 2 Sqrt[\[Eta] kCool] c[[1]]; kp = 2 Sqrt[\[Eta] kCool] c[[2]];
   s = LyapunovSolve[{{-gammaCool/2, 1}, {-1, -gammaCool/2 - gg}}, -{{kx^2, kx kp}, {kx kp, kp^2}}];
   (s[[1, 1]] + c[[1]] + s[[2, 2]] + c[[3]])/2 - 0.5];
etasCool = Range[0.2, 1, 0.1];
ListLinePlot[{Transpose[{etasCool, (#[[1]] + #[[3]])/2 - 0.5 & /@ (steadyCovCool /@ etasCool)}],
   Transpose[{etasCool, nUncCoolEta[#, Gcool] & /@ etasCool}]},
 PlotStyle -> {Directive[Black, Dashed], ColorData[97, 3]}, Frame -> True, GridLines -> Automatic, PlotRange -> {0, All},
 PlotLegends -> {"conditional floor", "cooled occupation (fixed gain)"}, ImageSize -> 500,
 FrameLabel -> {"detector efficiency", "phonon number"}, PlotLabel -> "a better detector cools colder"]
```

Both curves fall as the detector improves, and the gap between them is the price of noisy actuation. This is measurement-based cooling of a nanomechanical drum or a levitated nanoparticle toward its ground state, not by contact with something colder but by watching precisely and pushing back: the same conditional covariance that told the Kalman filter how well it could track now sets how cold the feedback can pull.

## Part Seven: Two Bookkeepings, and the Colour of the Light

### Linear vs Nonlinear: One Trajectory, Two Descriptions

**The problem.** One diffusive history, written two equivalent ways. The *linear* stochastic Schrodinger equation runs an unnormalized state under plain noise,
$$d\psi = K\,\psi\,dt + \sum_j R_j\,\psi\,dW_j,\qquad K = -i\hat H - \tfrac12\sum_j R_j^\dagger R_j,$$
whose squared norm $\|\psi\|^2$ is exactly the *likelihood* of the record that drove it. The *nonlinear* master equation instead carries the normalized state an experimenter tracks, driven not by that raw noise but by the **innovation** $d\widehat W_j = dW_j - v_j\,dt$, the increment with the conditional mean it predicts removed,
$$d\rho = \mathcal{L}\rho\,dt + \sum_j\bigl(R_j\rho + \rho R_j^\dagger - v_j\rho\bigr)\,d\widehat W_j,\qquad v_j = 2\,\mathrm{Re}\,\mathrm{Tr}[R_j\rho].$$
The startling fact is that the plain average of the *unnormalized* states, with no reweighting at all, already is the master equation. Build the linear version and check.

One diffusive trajectory has two equivalent descriptions. The **nonlinear** SME (used so far) carries the normalized state and renormalizes each step. The **linear** SSE evolves an unnormalized state under plain (unconditioned) noise, never renormalizing,
$$d|\tilde\psi\rangle = \big(-i\hat H - \tfrac12\hat c^\dagger\hat c\big)|\tilde\psi\rangle\,dt + \hat c\,|\tilde\psi\rangle\,dW,$$
its squared norm $\langle\tilde\psi|\tilde\psi\rangle$ the likelihood of the record. The striking fact: the plain average of the unnormalized outer products $|\tilde\psi\rangle\langle\tilde\psi|$, unweighted, is already the master-equation $\rho$. Fix the atom and the non-Hermitian drift $K = -i\hat H - \tfrac12\hat c^\dagger\hat c$:

```wl
\[CapitalOmega]lin = 2.; cLin = lower; Knh = -I (\[CapitalOmega]lin/2) X - ConjugateTranspose[cLin] . cLin/2;
step = 0.005; steps = 600;
```

One unnormalized run under plain noise, never renormalizing:

```wl
linearRun[seed_] := BlockRandom[SeedRandom[seed];
    Fold[#1 + Knh . #1 step + cLin . #1 #2 &, excited, RandomVariate[NormalDistribution[0, Sqrt[step]], steps]]];
```

Run four hundred of them:

```wl
linearStates = linearRun /@ Range[400];
```

Average the outer products, with no reweighting at all:

```wl
linearAverage = Mean[Outer[Times, #, Conjugate[#]] & /@ linearStates];
```

Measure the gap to the master equation:

```wl
Max@Abs@Flatten[linearAverage - evolve[(\[CapitalOmega]lin/2) X, {cLin}, densityMatrix[excited], steps step]]
```

The unweighted average matches the master equation to the Monte-Carlo scatter, with no renormalization.

The squared norms are genuine likelihoods, so they average to one:

```wl
Mean[Re[Conjugate[#] . #] & /@ linearStates]
```

Different records carry different weight.

Now visualize the distribution of their likelihood weights:

```wl
Histogram[Re[Conjugate[#] . #] & /@ linearStates, 30, "PDF", Frame -> True, GridLines -> Automatic, ImageSize -> 460,
 FrameLabel -> {"likelihood weight", "density"},
 PlotLabel -> "the likelihood weights average to one"]
```

These are two faithful descriptions of the same diffusive process, related by a change of probability measure. Renormalize the state each step and feed the record back into the drift, and you get the physical (nonlinear) state an experimenter tracks in real time. Evolve the unnormalized state under plain, unconditioned noise and never renormalize, and you get the linear picture, often easier to simulate. The linear state's *direction* still fluctuates with the noise; what the weight $\|\psi\|^2$ isolates is not the randomness but the physical likelihood of each record, the Girsanov factor that turns the reference measure into the physical one.

One record carries both bookkeepings at once. Roll out a single linear run keeping its whole history:

```wl
linPathA5[seed_] := BlockRandom[SeedRandom[seed];
    FoldList[#1 + Knh . #1 step + cLin . #1 #2 &, excited, RandomVariate[NormalDistribution[0, Sqrt[step]], steps]]];
oneA5 = linPathA5[7]; gridA5 = step Range[0, steps];
```

Read the two descriptions off that one run: the normalized state (the physical trajectory you would track) and the log-weight it accumulates (its likelihood score), on aligned axes:

```wl
Column[{
  ListLinePlot[Transpose[{gridA5, (1 + Re[Conjugate[#] . Z . #])/2 &[Normalize[#]] & /@ oneA5}], PlotStyle -> ColorData[97, 1],
   Frame -> True, GridLines -> Automatic, PlotRange -> {{0, steps step}, {0, 1.05}}, ImageSize -> 520,
   FrameLabel -> {"time", "excited population (normalized)"}, PlotLabel -> "the state you track"],
  ListLinePlot[Transpose[{gridA5, Log[Re[Conjugate[#] . #]] & /@ oneA5}], PlotStyle -> ColorData[97, 3],
   Frame -> True, GridLines -> Automatic, PlotRange -> {{0, steps step}, All}, ImageSize -> 520,
   FrameLabel -> {"time", "log weight  log \[LeftDoubleBracketingBar]\[Psi]\[RightDoubleBracketingBar]\!\(\*SuperscriptBox[\(\), \(2\)]\)"}, PlotLabel -> "the likelihood it carries"]}]
```

The same record has a normalized-state reading and a running weight; neither alone is the trajectory, together they are.

The catch of the linear method is that the weights spread, so a handful of records come to dominate the average, the importance-sampling degeneracy. Measure it with the effective sample size $N_{\mathrm{eff}} = (\sum w)^2/\sum w^2$ across the ensemble over time:

```wl
essPathsA5 = linPathA5 /@ Range[400];
essA5 = Table[With[{w = (Re[Conjugate[#] . #] &@#[[j]]) & /@ essPathsA5}, (Total[w]^2/Total[w^2])/400.], {j, 1, steps + 1, 20}];
ListLinePlot[Transpose[{gridA5[[1 ;; ;; 20]], essA5}], PlotStyle -> ColorData[97, 2], Frame -> True, GridLines -> Automatic,
 PlotRange -> {0, 1.02}, ImageSize -> 480, FrameLabel -> {"time", "effective sample size / N"},
 PlotLabel -> "the linear weights degenerate: fewer records carry the average"]
```

The effective fraction falls from one toward a fraction of the ensemble: as time goes on, the unweighted linear average is carried by ever fewer high-weight records, and its Monte-Carlo error grows even though it remains unbiased. That is the price of the linear picture's convenience, and the reason the nonlinear (renormalized) equation, which keeps every trajectory equally weighted, is what a long real-time filter uses.

### The Driven Atom, Averaged and Watched

**The problem.** The fully worked driven, damped, emitting atom, both averaged and watched. In full generality its unconditional master equation carries detuning, thermal drive, and pure dephasing,
$$\mathcal{L}\rho = -i[\check H,\rho] + \gamma(n_T+1)\,\mathcal{D}[\hat\sigma_-]\rho + \gamma\,n_T\,\mathcal{D}[\hat\sigma_+]\rho + \gamma\,k_d\,\mathcal{D}[\hat\sigma_z]\rho,\qquad \check H = \tfrac{\Delta\omega}{2}\hat\sigma_z + \tfrac{\Omega}{2}\hat\sigma_x,$$
and relaxes to a mixed steady state inside the ball. We take the standard resonance-fluorescence corner, resonant ($\Delta\omega = 0$), zero temperature ($\bar n = 0$), no extra dephasing ($k_d = 0$), which reduces it to a single emission channel $R_1 = \sqrt\gamma\,\hat\sigma_-$,
$$\mathcal{L}\rho = -i\Bigl[\tfrac{\Omega}{2}\hat\sigma_x,\rho\Bigr] + \gamma\,\mathcal{D}[\hat\sigma_-]\rho.$$
Watching that one channel with unit efficiency keeps the conditioned history pure on the surface of the ball; the record-blind average is the mixed interior point. Compare the two. (Restore any of $\Delta\omega$, $\bar n$, $k_d$ and the conditioned state would itself become mixed, since a thermal or dephasing channel is then present and unmonitored.)

Pull the threads together on one worked example: a driven, emitting atom. Averaged (unconditioned), it relaxes to a mixed steady state inside the Bloch ball. Fix the atom and compute its steady Bloch vector:

```wl
\[CapitalOmega]dr = 3.; \[Gamma]dr = 1.; settled = blochVector[First@steadyState[\[CapitalOmega]dr/2 X, {Sqrt[\[Gamma]dr] lower}]];
```

Read its length (a pure state has length one):

```wl
Norm[settled]
```

The length is well short of one: the averaged state is mixed.

But a conditioned trajectory, keeping the record, stays pure on the surface of the ball, fluctuating around that mixed average. Run one:

```wl
watched = trajectory[densityMatrix[excited], \[CapitalOmega]dr/2 X, {Sqrt[\[Gamma]dr] lower}, {1.}, {}, 0.01, 8., 5];
```

Confirm it stays pure, reading the largest linear entropy along the whole run:

```wl
Max[(1 - Norm[blochVector[#]]^2)/2 & /@ watched["states"]]
```

The linear entropy stays at zero to rounding.

Now visualize the pure trajectory on the Bloch-sphere surface, with the mixed averaged state marked inside:

```wl
blochPlot[{blochVector /@ watched["states"]}, "a pure trajectory circles the mixed average (black dot)", None,
 {Black, PointSize[0.035], Point[settled]}]
```

The surface path is the atom as an observer with the record sees it; the interior point is the atom traced over the record. Neither is more true; they answer different questions. This picture, a pure trajectory circling a mixed average, is the theme of all twenty systems.

### The Mollow Triplet: Resonance-Fluorescence Spectrum

**The problem.** The spectrum (power spectrum) of the light a driven atom emits, from the record alone. The light divides into two parts: an *elastic* (coherent) line at the drive frequency, a delta function carrying the weight $|\langle\hat\sigma_-\rangle_{\mathrm{ss}}|^2$, and the *inelastic* (fluorescence) spectrum of the fluctuations $\delta\hat\sigma_- = \hat\sigma_- - \langle\hat\sigma_-\rangle_{\mathrm{ss}}$, which by quantum regression is a resolvent of the Liouvillian at the steady state,
$$S_{\mathrm{inel}}(\mu) \propto \mathrm{Re}\int_0^\infty e^{i\mu t}\,\langle\delta\hat\sigma_+(t)\,\delta\hat\sigma_-(0)\rangle_{\mathrm{ss}}\,dt.$$
For weak driving the inelastic part is a single line; for strong driving $\Omega \gg \gamma$ it splits into the three-peaked *Mollow triplet*, a central line and two sidebands at $\pm\Omega$, the dynamical Stark effect.

Homodyne detection instead reads a chosen quadrature, and its spectrum can dip *below* the shot-noise floor: squeezed, nonclassical light. Build both from the same big matrix.

One last set of quantities from the master equation, all through the eigen-based spectrum tool built for the quantum point contact (`regressionSpectrum`). It needs no shifted frequency grid, because the steady mode is dropped explicitly. The inelastic fluorescence spectrum is its correlation of $\delta\hat\sigma_-$:

```wl
mollowInelastic[\[CapitalOmega]_, \[Gamma]_, tones_] := With[{big = liouvillian[(\[CapitalOmega]/2) X, {Sqrt[\[Gamma]] lower}, 2],
    rss = First@steadyState[(\[CapitalOmega]/2) X, {Sqrt[\[Gamma]] lower}]},
   With[{dm = lower - Tr[lower . rss] id2}, regressionSpectrum[big, dm . rss, ConjugateTranspose[dm], tones]]];
```

Sweep frequencies (centered on zero, no offset needed) for a weak and a strong drive:

```wl
tones = Range[-9, 9, 0.05];
{weak, strong} = {mollowInelastic[0.6, 1., tones], mollowInelastic[6., 1., tones]};
```

Weakly driven, the inelastic spectrum is a single line at the drive frequency. Strongly driven, it splits into the **Mollow triplet**: a central line and two sidebands at $\pm\Omega$.

Now visualize the weak- and strong-drive inelastic spectra together. Their absolute scales differ by orders of magnitude (the weak line would vanish on the strong axis), so plot each normalized to its own peak to compare the *shapes*, a single line against the triplet:

```wl
ListLinePlot[{Transpose[{tones, weak/Max[weak]}], Transpose[{tones, strong/Max[strong]}]},
 PlotStyle -> {ColorData[97, 1], ColorData[97, 2]}, Frame -> True, GridLines -> {{-6, 0, 6}, None},
 PlotLegends -> {"weak drive", "strong drive"}, ImageSize -> 560, PlotRange -> All,
 FrameLabel -> {"frequency (relative to drive)", "inelastic spectrum (each / its peak)"},
 PlotLabel -> "strong driving splits the inelastic spectrum into the Mollow triplet"]
```

The sidebands sit near $\pm\Omega$: under strong driving the atom's levels are dressed into a four-level ladder, whose transitions produce the three lines. Confirm the sidebands land at the drive strength:

```wl
Select[tones[[#]] & /@ FindPeaks[strong][[All, 1]], Abs[Abs[#] - 6] < 0.6 &]
```

The peaks sit at $\pm\Omega$. Confirm too that the inelastic spectrum is nowhere negative, as a genuine power spectrum must be:

```wl
Min[strong]
```

Nonnegative to rounding.

The elastic line lives separately, a delta at zero frequency whose weight is the coherent dipole; it dominates at weak drive and all but vanishes at strong drive as the light turns fully inelastic:

```wl
elasticWeight[\[CapitalOmega]_] := With[{rss = First@steadyState[(\[CapitalOmega]/2) X, {Sqrt[1.] lower}]}, Abs[Tr[lower . rss]]^2];
{elasticWeight[0.6], elasticWeight[6.]}
```

The coherent weight falls by more than an order of magnitude from weak to strong drive. To watch the triplet being born, sweep the drive continuously and stack the inelastic spectra into a heatmap. The central peak dwarfs the sidebands at every drive, so scale each row to its own peak to keep the sidebands visible, and overlay the dressed-transition lines $\mu = \pm\Omega$ as guides:

```wl
mollowDrives = Range[0.3, 8, 0.25];
mollowHeat = (#/Max[#]) & /@ (mollowInelastic[#, 1., Range[-11, 11, 0.1]] & /@ mollowDrives);
ArrayPlot[mollowHeat, DataRange -> {{-11, 11}, {First@mollowDrives, Last@mollowDrives}}, DataReversed -> True,
 ColorFunction -> "SunsetColors", Frame -> True, AspectRatio -> 1.1, ImageSize -> 460,
 Epilog -> {Directive[White, Dashed, Opacity[0.7]], Line[{{0.3, 0.3}, {8, 8}}], Line[{{-0.3, 0.3}, {-8, 8}}]},
 FrameLabel -> {"frequency", "drive \[CapitalOmega]"}, PlotLabel -> "the Mollow triplet is born as the drive grows (rows scaled to peak; \[Mu]=\[PlusMinus]\[CapitalOmega] dashed)"]
```

The bright ridge splits into three, the outer two fanning apart along $\mu = \pm\Omega$: the dressed-state sidebands. The triplet is what an ordinary spectrometer sees.

A homodyne detector sees something an intensity spectrometer cannot: it measures one quadrature $\hat\sigma_- e^{-i\theta} + \hat\sigma_+ e^{i\theta}$, and its spectrum, normalized so the shot-noise floor is one, can fall below that floor. Build it from the measurement current fluctuations:

```wl
mollowHom[\[CapitalOmega]_, \[Gamma]_, \[Theta]_, tones_] := Module[{c = Sqrt[\[Gamma]] lower, big, rss, cth, mq, mean},
   big = liouvillian[(\[CapitalOmega]/2) X, {c}, 2]; rss = First@steadyState[(\[CapitalOmega]/2) X, {c}];
   cth = c Exp[-I \[Theta]]; mq = cth + ConjugateTranspose[cth]; mean = Re@Tr[mq . rss];
   1 + regressionSpectrum[big, cth . rss + rss . ConjugateTranspose[cth] - mean rss, mq - mean id2, tones]];
```

At weak drive, one quadrature is squeezed and the perpendicular one is not. Compare the two at $\Omega = 0.25\gamma$:

```wl
ListLinePlot[{Transpose[{tones, mollowHom[0.25, 1., Pi/2, tones]}], Transpose[{tones, mollowHom[0.25, 1., 0., tones]}]},
 PlotStyle -> {ColorData[97, 1], ColorData[97, 2]}, Frame -> True, GridLines -> {None, {{1, Directive[Gray, Dashed]}}},
 PlotLegends -> {"squeezed quadrature (\[Theta] = \[Pi]/2)", "anti-squeezed (\[Theta] = 0)"}, PlotRange -> {0.7, All}, ImageSize -> 560,
 FrameLabel -> {"frequency", "homodyne spectrum / shot noise"}, PlotLabel -> "resonance fluorescence dips below the shot-noise floor"]
```

The $\theta = \pi/2$ quadrature dips below the dashed shot-noise line while its partner rises above it: the atom's fluorescence is squeezed, its noise pushed below the vacuum in one quadrature at the cost of excess in the other. Confirm the dip is real:

```wl
Min[mollowHom[0.25, 1., Pi/2, tones]]
```

Below one: sub-shot-noise, a signature of nonclassical light that only phase-sensitive (homodyne) detection can reveal, and that the necessarily-nonnegative intensity spectrum of the triplet can never show. The whole spectroscopy, elastic line, inelastic triplet, and squeezing dip, comes from the two-time correlations of the same master equation built on the first page.

## Where This Leaves Us

From plain matrices and a few rules, we built a small laboratory for open quantum systems: the master equation for the unconditioned average, the positivity-preserving step for a kept record, and plotters for the Bloch sphere, phase space, the covariance ellipse, and spectra. On that base, twenty systems: pure dephasing, a driven damped atom, a cavity cat decohering at a rate set by the square of its separation, the Haroche experiment, dispersive qubit readout, quantum Brownian motion, the jump/homodyne/heterodyne unravellings, the quantum Zeno effect, a charge qubit read by a quantum point contact, measurement-induced collapse, a quantum Kalman filter, measurement feedback, thermalization, rapid purification, feedback cooling, the linear-vs-nonlinear trajectory, and the Mollow triplet. One theme runs through all of them: keeping the record turns the master equation's average into a pure, jagged, individual trajectory, and averaging the trajectories back recovers the master equation. Change a number in any cell and watch it move.




