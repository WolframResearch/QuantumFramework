---
Template: Default
---

# Watching a Qubit: Continuous Measurement and Quantum Trajectories

**A live, computation-first walk through the diffusive stochastic master equation of a continuously monitored qubit in the Wolfram Language: the deterministic Lindblad equation describes the ensemble that ignores the detector, while conditioning the state on a homodyne readout unravels it into individual stochastic quantum trajectories; a positivity-preserving integrator, built from plain matrices, turns each detector record into a trajectory; and averaging the trajectories over the measurement noise recovers the Lindblad state, with measurement-induced dephasing, spontaneous-emission diffusion, a two-point correlation of the readout current, inefficient detection, and a qutrit all read off cells you can rerun.**

Mads Bahrami. Last updated July 29, 2026.

### Setting the Stage: How This Notebook Flows

This notebook tells one story as a continuous sequence: what it means to watch a quantum system weakly and continuously, how that watching turns the deterministic Lindblad master equation into a stochastic master equation (SME) for the state conditioned on the detector's record, how a single positivity-preserving update rule integrates the SME from plain matrices, and how the same map, run over and over, reproduces measurement-induced dephasing, spontaneous-emission diffusion, the two-point statistics of the homodyne current, and, once the map is dimension-general and efficiency-aware, inefficient detection and a qutrit.

I strongly believe in a computation-first narrative for learning: in a sense, if I cannot compute it, I cannot claim to understand it. So almost nothing here is asserted without a cell that computes it on the spot, and almost every claim is checked against an independent calculation, an ensemble average against a master-equation solution, or a positivity test against the scheme's promise.

This is a live Wolfram notebook, not a static paper. Evaluate the cells top to bottom; each cell holds one computation and depends on symbols defined earlier, so a fresh kernel wants a fresh top-to-bottom pass. Read it as a movie rather than a reference: each cell sets up the next question, and the prose between cells is the connective tissue, not decoration. My suggestion is to focus on each output and its meaning before worrying about every detail of the input. And remember that the code is yours: change the measurement strength, the drive, the initial state, the number of trajectories, rerun, and see what moves.

The integrator, the object this notebook is really about, is a short function over plain complex matrices with no dependency on the QuantumFramework paclet. Two convenience pieces come from the Wolfram Function Repository: `ResourceFunction["LindbladSolve"]` integrates the unconditional master equation we check the trajectories against, and `ResourceFunction["BlochSpherePlot"]` draws the sphere. Both install automatically the first time you evaluate them, so the first run needs a network connection. Prerequisites are light: basic quantum mechanics (states, density matrices, the Pauli operators) and enough Wolfram Language to read a one-line function. No stochastic calculus is needed beyond one fact, stated when we first use it: over a step of length $dt$ the Wiener increment $dW$ is a zero-mean Gaussian of variance $dt$. Throughout we work in natural units with $\hbar = 1$, so Hamiltonians are measured in frequency and rates and times are reciprocal.

Let's start!

## Part I: Two Descriptions of a Monitored Qubit

Start with the object every open-system calculation begins from. A quantum system coupled to an environment, and left unobserved, has a state $\rho$ (a density matrix) that obeys the Lindblad master equation

$$
\frac{d\rho}{dt} = -i[H,\rho] + \sum_k \mathcal{D}[c_k]\rho , \qquad \mathcal{D}[c]\rho = c\,\rho\,c^\dagger - \tfrac12\{c^\dagger c,\rho\} ,
$$

where $H$ is the Hamiltonian, the $c_k$ are the jump (Lindblad) operators through which the environment acts, and $\mathcal{D}[c]$ is the dissipator (the superoperator that carries irreversible loss; the anticommutator $\{c^\dagger c,\rho\}$ keeps the trace fixed). In words: this is the equation for the state of a system whose environment is there but not being watched, and it is deterministic and smooth. Its solution is the *unconditional* state, the average over everything the environment could have done.

Now change one thing: put a detector on one of those channels. Suppose an experimenter monitors the system through a jump operator $L$, for instance by collecting the light it scatters and mixing it with a local oscillator (homodyne detection). The detector returns a continuous, noisy signal, and the state of the system, given that signal, is no longer the smooth unconditional $\rho$. It is a *conditional* state that jitters with the record. Splitting the jump operators into monitored ones $L_r$ (with detector efficiency $\eta_r$) and unmonitored ones $V_j$ (environment we cannot see), the conditional state obeys the diffusive **stochastic master equation**

$$
d\rho = -i[H,\rho]\,dt + \sum_j \mathcal{D}[V_j]\rho\,dt + \sum_r \mathcal{D}[L_r]\rho\,dt + \sum_r \sqrt{\eta_r}\;\mathcal{H}[L_r]\rho\;dW_r ,
$$

with the measurement superoperator (the innovation term) $\mathcal{H}[c]\rho = c\rho + \rho c^\dagger - \mathrm{Tr}[(c+c^\dagger)\rho]\,\rho$, and each detector emitting a readout increment

$$
dR_r = \sqrt{\eta_r}\,\mathrm{Tr}[(L_r + L_r^\dagger)\rho]\,dt + dW_r .
$$

In other words, the SME is the master equation plus one stochastic term per monitored channel: $\mathcal{H}[L_r]$ nudges the state toward agreement with the reading, weighted by the Wiener increment $dW_r$, which is the detector's shot noise. That same $dW_r$ appears in the record $dR_r$, so the record is signal (the mean $\mathrm{Tr}[(L+L^\dagger)\rho]$) buried in noise, and the state can be reconstructed from the record. The rates are folded into the operators throughout, $L_r \leftarrow \sqrt{\gamma_r}\,L_r$, so a monitored channel of rate $\gamma$ and unit coupling is passed as $\sqrt{\gamma}\,L$; this is the convention the code below expects.

Two limits are worth naming before we compute. If no channel is monitored, every $\sqrt{\eta_r}\to 0$ and the SME collapses back to the Lindblad master equation: the unwatched state. And if we run the SME for many independent noise realizations and average, the $dW_r$ terms average to zero (a Wiener increment has zero mean), so the mean conditional state obeys the master equation again. This is the single most important fact in the subject, and we will verify it directly: **one trajectory is stochastic and, at unit efficiency, pure; the average of many trajectories is the deterministic, generally mixed, master-equation state.** We will return to the readout $dR_r$ in Part VI, where it becomes the object of study.

One subtlety underlies everything that follows, so state it plainly. The only thing a detector *directly* produces is the record $R_t$: one real, noisy number per time step. The conditional state $\rho_c(t)$, and therefore its Bloch vector, its purity, and every $\langle A\rangle=\mathrm{Tr}[\rho_c A]$ we will plot from a single run, is not itself measured but *inferred*, and the SME is exactly the optimal filter that does the inferring: it is the quantum analogue of a Kalman filter, turning the record into the observer's best estimate of the state. So when a plot below shows $\langle Z\rangle(t)$ for one trajectory, read it as a quantity an experimenter reconstructs by feeding a recorded current through this same equation, not as a raw meter reading. **The record is the observable; the state, and everything computed from it, is what you infer from the record.** That is why Part VI, which analyzes the record directly, is where the essay finally touches measured data, and why the menu of things one can do with a record (estimate the state, its spectrum, the system's parameters, or a feedback signal) is worth laying out once we have one in hand.

## Part II: A Positivity-Preserving Update

A naive Euler step of the SME can push $\rho$ to a matrix with a negative eigenvalue, an unphysical state, especially at large $dt$. The scheme we use avoids this by construction. Over one step it forms an *unnormalized* conditional state $\tilde\rho_{n+1}$ and divides it by its own trace,

$$
\rho_{n+1} = \frac{\tilde\rho_{n+1}}{\mathrm{Tr}[\tilde\rho_{n+1}]} ,
$$

$$
\tilde\rho_{n+1} = M\,\rho_n\,M^\dagger + dt\sum_j V_j\,\rho_n\,V_j^\dagger + dt\sum_r (1-\eta_r)\,L_r\,\rho_n\,L_r^\dagger ,
$$

$$
M = 1 - \Bigl(iH + \tfrac12\textstyle\sum_r L_r^\dagger L_r + \tfrac12\sum_j V_j^\dagger V_j\Bigr)dt \;+\; \sum_r \sqrt{\eta_r}\,L_r\,dy_r \;+\; \tfrac12\sum_{r,s}\sqrt{\eta_r\eta_s}\,L_r L_s\,(dy_r dy_s - \delta_{rs}\,dt) ,
$$

with $dy_r = \sqrt{\eta_r}\,\mathrm{Tr}[(L_r+L_r^\dagger)\rho_n]\,dt + dW_r$. This positivity-preserving update is based on a paper of Rouchon and Ralph, [arXiv:1410.5345](https://arxiv.org/abs/1410.5345). In words: $\tilde\rho_{n+1}$ is a sum of terms of the form $A\,\rho_n\,A^\dagger$, and any such sum is positive semidefinite whenever $\rho_n$ is, so $\tilde\rho_{n+1}$ is a valid unnormalized state at *any* step size, and dividing by its trace $\mathrm{Tr}[\tilde\rho_{n+1}]$ restores unit trace. The unmonitored channels and the unmonitored fraction $(1-\eta_r)$ of each detector re-enter as ordinary dissipators, the light that leaks away unrecorded.

The paper keeps that double sum $\tfrac12\sum_{r,s}\sqrt{\eta_r\eta_s}\,L_r L_s\,(dy_r dy_s - \delta_{rs}\,dt)$; the code groups it into a single square. Writing $S = \sum_r \sqrt{\eta_r}\,L_r\,dy_r$, the scalar increments $dy_r$ factor out, so $\sum_{r,s}\sqrt{\eta_r\eta_s}\,L_r L_s\,dy_r dy_s = S^2$ identically (the operator order $L_r L_s$ is preserved), while the $\delta_{rs}$ diagonal contributes $-\tfrac12 dt\sum_r \eta_r L_r^2$. The double sum is therefore exactly $\tfrac12 S^2 - \tfrac12 dt\sum_r \eta_r L_r^2$, an algebraic identity valid for any number of monitored operators, not an approximation, so

$$
M = H_{\rm eff} + S + \tfrac12 S^2 - L_{\rm corr}, \quad H_{\rm eff} = 1 - \Bigl(iH + \tfrac12\textstyle\sum_r L_r^\dagger L_r + \tfrac12\sum_j V_j^\dagger V_j\Bigr)dt, \quad L_{\rm corr} = \tfrac12 dt\sum_r \eta_r L_r^2 ,
$$

and $L_{\rm corr}$ is the Ito correction. This grouped form is algebraically identical to the double-sum $M$ above, and it is the form the code builds. Everything that does not change from step to step ($H_{\rm eff}$, $L_{\rm corr}$) is computed once, and the per-step work is one readout, one contraction $S$, and one $M\rho M^\dagger$.

Begin with the readout increments of one step, the vector $dy$ with one entry per monitored operator. This is $\mathcal{R}$: it takes the current state, the monitored operators, the step, the Wiener draws $dw$, and the efficiencies, and returns $dy_r = \sqrt{\eta_r}\,\mathrm{Tr}[(L_r+L_r^\dagger)\rho]\,dt + dW_r$:

```wl
ClearAll[\[ScriptCapitalR], \[ScriptCapitalT]];
\[ScriptCapitalR][ρ_, Ls_, dt_, dw_, η_] :=
   MapThread[Sqrt[#3] Tr[(#1 + ConjugateTranspose[#1]) . ρ] dt + #2 &, {Ls, dw, η}];
```

Now the trajectory map itself. This is $\mathcal{T}$: given an initial density matrix `ρ0`, a Hamiltonian `H`, a list of monitored operators `Ls` (already scaled by $\sqrt{\gamma}$), optional efficiencies `η` and unmonitored operators `V`, a step `dt`, and a final time `tf`, it returns the list of conditioned density matrices, one per time point, and `Sow`s the readout vector $dy$ at every step so a `Reap` around the call recovers the record. Read it as the collapsed scheme above: precompute `heff` and `lcorr`, draw all the noise at once, and fold the per-step update `step` over the noise with `FoldList`:

```wl
\[ScriptCapitalT][ρ0_?MatrixQ, H_?MatrixQ, Ls_List, η : (None | _List) : None,
   V : (None | _List) : None, dt_?NumericQ, tf_?NumericQ] :=
  Module[
   {d = Length[ρ0], r0 = N[ρ0], Hn = N[H], Lsn = N[Ls], dtn = N[dt],
    ηv, Vsn, heff, lcorr, nsteps, dw, step},
   ηv = N[If[η === None, ConstantArray[1, Length[Ls]], η]];
   Vsn = If[V === None, {}, N[V]];
   heff = IdentityMatrix[d] - dtn (I Hn
       + (1/2) Total[ConjugateTranspose[#] . # & /@ Lsn]
       + (1/2) Total[ConjugateTranspose[#] . # & /@ Vsn]);
   lcorr = (dtn/2) Total[MapThread[#1 (#2 . #2) &, {ηv, Lsn}]];
   nsteps = Floor[tf/dt];
   dw = Transpose[RandomVariate[NormalDistribution[0, Sqrt[dtn]], {Length[Lsn], nsteps}]];
   step[ρ_, dwv_] := Module[{dy, s, m, num},
      dy = Sow[\[ScriptCapitalR][ρ, Lsn, dtn, dwv, ηv]];
      s = Total[MapThread[#1 Sqrt[#2] #3 &, {dy, ηv, Lsn}]];
      m = heff + s + (1/2) s . s - lcorr;
      num = m . ρ . ConjugateTranspose[m]
          + dtn Total[(# . ρ . ConjugateTranspose[#]) & /@ Vsn]
          + dtn Total[MapThread[(1 - #2) (#1 . ρ . ConjugateTranspose[#1]) &, {Lsn, ηv}]];
      num/Tr[num]];
   FoldList[step, r0, dw]];
```

A few design choices are worth reading off the code. The inputs are numericalized once up front (`N[ρ0]`, `N[H]`, `N[Ls]`), so passing exact matrices like $|+\rangle\langle+|$ costs nothing and the hot loop runs on packed numeric arrays. Nothing in the body refers to any external paclet: `ρ`, `H`, and the operators are ordinary dense complex matrices. The dimension `d` is read from the state, so the same function runs a qubit, a qutrit, or two qubits without change, which is the freedom Part VII spends. And the pattern guards (`_?MatrixQ`, `_List`, `(None | _List)`) let the two middle arguments be optional while keeping the call unambiguous, so `\[ScriptCapitalT][ρ0, H, Ls, dt, tf]` and `\[ScriptCapitalT][ρ0, H, Ls, η, V, dt, tf]` both parse correctly.

Next the primitives everything is built from. Define the Pauli matrices and the identity:

```wl
X = PauliMatrix[1]; Y = PauliMatrix[2]; Z = PauliMatrix[3]; Id = IdentityMatrix[2];
```

The qubit lowering operator is $J_- = (X - iY)/2$. Define it:

```wl
Jminus = (X - I Y)/2;
```

See what it does: display $J_-$ as a matrix and its action on the two basis states:

```wl
Column[{
   Row[{Subscript["J", "-"], " = ", MatrixForm[Jminus]}],
   Row[{Subscript["J", "-"], "|0\[RightAngleBracket] = ", Jminus . {1, 0}}],
   Row[{Subscript["J", "-"], "|1\[RightAngleBracket] = ", Jminus . {0, 1}}]}]
```

As one can see, $J_- = \left(\begin{smallmatrix}0&0\\1&0\end{smallmatrix}\right)$ sends $|0\rangle\to|1\rangle$ and annihilates $|1\rangle$. Note the convention this fixes: here $|0\rangle$ is the *upper* state that decays and $|1\rangle$ is the ground (dark) state, so spontaneous emission through $J_-$ drives the qubit toward $|1\rangle$. Some literature labels these the other way; we will keep $|0\rangle$ excited throughout.

A density matrix is the outer product $|\psi\rangle\langle\psi|$ of a state vector with its conjugate. Define that map:

```wl
dm[ψ_] := Outer[Times, ψ, Conjugate[ψ]];
```

Name the three basis states we will start from, the excited state, the ground state, and the equatorial superposition:

```wl
ket0 = {1, 0}; ket1 = {0, 1}; ketPlus = {1, 1}/Sqrt[2];
```

Show the density matrices of $|0\rangle$ and $|+\rangle$:

```wl
Column[{
   Row[{"|0\[RightAngleBracket]\[LeftAngleBracket]0| = ", MatrixForm[dm[ket0]]}],
   Row[{"|+\[RightAngleBracket]\[LeftAngleBracket]+| = ", MatrixForm[dm[ketPlus]]}]}]
```

Finally, the geometry. A qubit state $\rho$ is fully described by its Bloch vector $(\langle X\rangle, \langle Y\rangle, \langle Z\rangle)$ with $\langle A\rangle = \mathrm{Tr}[\rho A]$, real because $\rho$ and the Paulis are Hermitian. Define the Bloch-vector extractor:

```wl
blochVector[ρ_] := Re[Tr[ρ . #] & /@ {X, Y, Z}];
```

Check that $|0\rangle$ sits at the north pole $(0,0,1)$:

```wl
blochVector[dm[ket0]]
```

For the sphere itself we use the Function Repository resource `ResourceFunction["BlochSpherePlot"]`: called with no argument, `ResourceFunction["BlochSpherePlot"][]` draws an empty Bloch sphere, onto which we overlay each run's full Bloch trajectory as a curve with `Show`. That is the entire toolkit: the map $\mathcal{T}$, the readout $\mathcal{R}$, the primitives, the Bloch vector, and the sphere plotter. Everything below is these with a different Hamiltonian, jump operator, or efficiency.

## Part III: The Lindblad Limit: One Noisy Run and the Average of Many

Recall the central claim of Part I: one trajectory is stochastic, but the average of many is the deterministic master-equation state. This Part checks it directly, on the simplest interesting case, a driven qubit whose $Z$ is continuously measured, with the Hamiltonian a Rabi drive about $x$ and the monitored operator $Z$ at a rate large compared with the drive. Fix the drive, the monitored operator, the initial state, the step, and the horizon:

```wl
Ω = 1.0; H = Ω X; γ = 2.0; Ls = {Sqrt[γ] Z}; ρ0 = dm[ket0]; δt = 0.005; tf = 20.0;
```

To have something to compare the ensemble against, we need the unconditional (Lindblad) solution. The Function Repository's `LindbladSolve` integrates the master equation from the Hamiltonian, the jump operators, and the initial state, returning $\rho(t)$ as an interpolating function. Solve it for this scenario:

```wl
lind = ResourceFunction["LindbladSolve"][H, Ls, ρ0, {t, 0, tf}];
```

Set the time grid the trajectories will use:

```wl
tgrid = δt Range[0, Floor[tf/δt]];
```

Read the Lindblad Bloch trajectory off the interpolating solution:

```wl
lindBloch = blochVector[lind[#]] & /@ tgrid;
```

Now one conditioned run. Generate one trajectory on a fixed seed so the picture is reproducible:

```wl
oneRun = BlockRandom[SeedRandom[1]; \[ScriptCapitalT][ρ0, H, Ls, δt, tf]];
```

Extract its Bloch trajectory:

```wl
oneBloch = blochVector /@ oneRun;
```

Before plotting, cash in the promise of the scheme: every state in the run should be a legitimate density matrix. Check that all of them are positive semidefinite:

```wl
Tally[PositiveSemidefiniteMatrixQ /@ oneRun]
```

Every state passes, which is the scheme's positivity guarantee made concrete: the numerator is a sum of $A\rho A^\dagger$ terms, so positivity never fails, here or at any step size. Now overlay the single conditioned $\langle Z\rangle$ on the smooth Lindblad $\langle Z\rangle$:

```wl
ListLinePlot[{Transpose[{tgrid, oneBloch[[All, 3]]}], Transpose[{tgrid, lindBloch[[All, 3]]}]},
 PlotStyle -> {Directive[Opacity[0.6], ColorData[97, 1]], Directive[Thick, Red]},
 PlotLegends -> {"one trajectory ⟨Z⟩", "Lindblad ⟨Z⟩"}, Frame -> True, ImageSize -> 540,
 AspectRatio -> 1/2, GridLines -> Automatic, PlotRange -> {-1, 1.05},
 FrameLabel -> {"time t", "⟨Z⟩"},
 PlotLabel -> "One monitored run jitters around the smooth master-equation curve"]
```

As one can see, the single trajectory fluctuates continuously while the Lindblad curve slides smoothly downward. The strong $Z$ measurement pins the qubit near a $Z$ eigenstate and only slowly lets the drive tip it over, so both curves decay gently rather than oscillating (we will contrast this with the weak-measurement regime in Part IV). The trajectory is one realization of the detector's noise; the Lindblad curve is what remains when the noise is averaged away. To make that precise, average many trajectories. Run an ensemble of 200 conditioned trajectories on a fixed seed:

```wl
ensemble = BlockRandom[SeedRandom[42]; Table[\[ScriptCapitalT][ρ0, H, Ls, δt, tf], {200}]];
```

Average them into the mean state at each time:

```wl
meanStates = Mean[ensemble];
```

Read its Bloch trajectory:

```wl
meanBloch = blochVector /@ meanStates;
```

Compare the ensemble mean of $\langle Z\rangle$ against the Lindblad prediction:

```wl
ListLinePlot[{Transpose[{tgrid, meanBloch[[All, 3]]}], Transpose[{tgrid, lindBloch[[All, 3]]}]},
 PlotStyle -> {Directive[Thick, ColorData[97, 1]], Directive[Dashed, Red]},
 PlotLegends -> {"mean of 200 trajectories ⟨Z⟩", "Lindblad ⟨Z⟩"}, Frame -> True,
 ImageSize -> 540, AspectRatio -> 1/2, GridLines -> Automatic, PlotRange -> {-1, 1.05},
 FrameLabel -> {"time t", "⟨Z⟩"}, PlotLabel -> "200 trajectories average to the master equation"]
```

Report the largest discrepancy across the whole trajectory:

```wl
Max@Abs[meanBloch[[All, 3]] - lindBloch[[All, 3]]]
```

The mean tracks the Lindblad curve, and the residual sits at the Monte-Carlo scale $1/\sqrt{K}\approx 0.07$ for $K=200$ realizations; push $K$ higher and it shrinks as $1/\sqrt{K}$. This is the trajectory unravelling of the master equation made numerical: the smooth, deterministic Lindblad state is nothing but the average over these jittering conditioned runs, each of which is a state an experimenter could actually infer from a detector record. The next question is what the individual runs look like when the drive, not the measurement, dominates.

## Part IV: Fast versus Slow Monitoring: Measurement-Induced Dephasing

The comparison that makes continuous measurement vivid is strong versus weak monitoring of the same driven qubit. Recall the Bloch equations for $H=\Omega X$ with $Z$ measured: the drive rotates $\langle Y\rangle$ and $\langle Z\rangle$ into each other at frequency $2\Omega$, while the $Z$ dissipator damps $\langle X\rangle$ and $\langle Y\rangle$ at rate $2\gamma$. The ratio $\gamma/\Omega$ decides which wins. Keep the drive fixed and take two rates, one well above and one well below it:

```wl
Ωd = 1.0; Hd = Ωd X; ρ0d = dm[ket0]; δtd = 0.005; tfd = 20.0; γFast = 2.0; γSlow = 0.1;
```

Set the shared time grid:

```wl
tgridd = δtd Range[0, Floor[tfd/δtd]];
```

First the master-equation view, which is deterministic and needs no averaging. Solve Lindblad at the fast rate:

```wl
lindFast = ResourceFunction["LindbladSolve"][Hd, {Sqrt[γFast] Z}, ρ0d, {t, 0, tfd}];
```

Read its $\langle Z\rangle$ on the grid:

```wl
zFast = blochVector[lindFast[#]][[3]] & /@ tgridd;
```

Solve Lindblad at the slow rate:

```wl
lindSlow = ResourceFunction["LindbladSolve"][Hd, {Sqrt[γSlow] Z}, ρ0d, {t, 0, tfd}];
```

Read its $\langle Z\rangle$:

```wl
zSlow = blochVector[lindSlow[#]][[3]] & /@ tgridd;
```

Plot the two unconditional curves together:

```wl
ListLinePlot[{Transpose[{tgridd, zFast}], Transpose[{tgridd, zSlow}]},
 PlotStyle -> {Directive[Thick, ColorData[97, 1]], Directive[Thick, ColorData[97, 2]]},
 PlotLegends -> {"γ = 2 ≫ Ω (frozen)", "γ = 0.1 ≪ Ω (Rabi)"}, Frame -> True, ImageSize -> 540,
 AspectRatio -> 1/2, GridLines -> Automatic, PlotRange -> {-1.05, 1.05},
 FrameLabel -> {"time t", "⟨Z⟩"},
 PlotLabel -> "Strong measurement freezes the qubit; weak measurement lets it oscillate"]
```

As one can see, the two regimes could hardly look more different. At $\gamma\gg\Omega$ the qubit is nearly frozen: strong $Z$ monitoring keeps collapsing it back toward a $Z$ eigenstate faster than the drive can rotate it away, so $\langle Z\rangle$ decays slowly and monotonically at the Zeno rate $\sim\Omega^2/\gamma$. This is the quantum Zeno effect: watching a system in a basis suppresses the coherent evolution out of that basis. At $\gamma\ll\Omega$ the drive wins and $\langle Z\rangle$ executes Rabi oscillations at $2\Omega$, their envelope slowly eroded by the weak dephasing at rate $\sim\gamma$.

The trajectory view shows the same physics one record at a time. Take one conditioned run at the fast rate:

```wl
runFast = BlockRandom[SeedRandom[7]; \[ScriptCapitalT][ρ0d, Hd, {Sqrt[γFast] Z}, δtd, tfd]];
```

And one at the slow rate, on the same seed:

```wl
runSlow = BlockRandom[SeedRandom[7]; \[ScriptCapitalT][ρ0d, Hd, {Sqrt[γSlow] Z}, δtd, tfd]];
```

Overlay their $\langle Z\rangle$:

```wl
ListLinePlot[{Transpose[{tgridd, blochVector[#][[3]] & /@ runFast}],
   Transpose[{tgridd, blochVector[#][[3]] & /@ runSlow}]},
 PlotStyle -> {Directive[Opacity[0.7], ColorData[97, 1]], Directive[Opacity[0.7], ColorData[97, 2]]},
 PlotLegends -> {"γ = 2 trajectory", "γ = 0.1 trajectory"}, Frame -> True, ImageSize -> 540,
 AspectRatio -> 1/2, GridLines -> Automatic, PlotRange -> {-1.05, 1.05},
 FrameLabel -> {"time t", "⟨Z⟩"}, PlotLabel -> "One record apiece: frozen near ±1, or noisily oscillating"]
```

The strongly measured run clings near $\pm 1$, taking rare, quick excursions as the drive occasionally flips it: the record is telling the experimenter, most of the time, which $Z$ eigenstate the qubit is in. The weakly measured run instead follows the Rabi oscillation, dressed in noise, because each weak reading barely perturbs the coherent motion. Same drive, same detector, opposite behavior, set entirely by how hard we look.

## Part V: Spontaneous Emission: Diffusion on the Bloch Sphere

Now let the monitored channel be spontaneous emission rather than dephasing, so the jump operator is the lowering operator $J_-$. Homodyne detection of the emitted field gives a diffusive unravelling, and the conditioned state wanders continuously over the Bloch sphere instead of jumping. Take a driven, detuned qubit emitting at a modest rate, started on the equator:

```wl
Ωs = 2.0; Δs = 1.0; Hs = (Ωs X + Δs Z)/2; γs = 0.2; Lse = {Sqrt[γs] Jminus}; ρ0s = dm[ketPlus]; δts = 0.01; tfs = 10.0;
```

Set the time grid:

```wl
tgrids = δts Range[0, Floor[tfs/δts]];
```

The unconditional state decays toward the ground state $|1\rangle$ while the drive stirs it. Solve its master equation:

```wl
lindS = ResourceFunction["LindbladSolve"][Hs, Lse, ρ0s, {t, 0, tfs}];
```

Read its Bloch trajectory:

```wl
lindBlochS = blochVector[lindS[#]] & /@ tgrids;
```

Individual trajectories show the diffusion directly. Generate a handful of conditioned runs and take their full Bloch trajectories:

```wl
fewRuns = BlockRandom[SeedRandom[3]; Table[blochVector /@ \[ScriptCapitalT][ρ0s, Hs, Lse, δts, tfs], {5}]];
```

Draw the sphere with `BlochSpherePlot[]` and overlay each run's full trajectory as a curve:

```wl
Show[ResourceFunction["BlochSpherePlot"][],
 Graphics3D[{Thick, MapIndexed[{ColorData[97, First@#2], Line@#1} &, fewRuns]}],
 ImageSize -> 420]
```

As one can see, each conditioned run traces its own wandering curve over the sphere, drifting in from the equatorial start toward the ground-state pole as the qubit emits, and no two curves coincide because each is driven by its own detector noise. This is measurement-induced diffusion: the homodyne record continuously reshapes the state, and the spread of the curves is the spread of possible records. Now average many of them and confirm the mean returns the master equation. Run an ensemble of 200 trajectories:

```wl
ensembleS = BlockRandom[SeedRandom[11]; Table[\[ScriptCapitalT][ρ0s, Hs, Lse, δts, tfs], {200}]];
```

Read the mean Bloch trajectory:

```wl
meanBlochS = blochVector /@ Mean[ensembleS];
```

Plot the three ensemble-mean components against the Lindblad prediction:

```wl
Show[
 ListLinePlot[Transpose[{tgrids, #}] & /@ Transpose[meanBlochS],
  PlotStyle -> {ColorData[97, 1], ColorData[97, 2], ColorData[97, 3]}],
 ListLinePlot[Transpose[{tgrids, #}] & /@ Transpose[lindBlochS], PlotStyle -> Directive[Dashed, Black]],
 Frame -> True, ImageSize -> 540, AspectRatio -> 1/2, GridLines -> Automatic,
 PlotRange -> All, FrameLabel -> {"time t", "Bloch components"},
 PlotLegends -> SwatchLegend[{ColorData[97, 1], ColorData[97, 2], ColorData[97, 3], Black},
   {"⟨X⟩ mean", "⟨Y⟩ mean", "⟨Z⟩ mean", "Lindblad (all three)"}],
 PlotLabel -> "200 emission trajectories average to the master equation, all three components"]
```

The three ensemble means fall on the three dashed Lindblad curves to Monte-Carlo scatter: the diffusing pure trajectories, averaged, reconstruct the decaying mixed state exactly as in Part III, now for an amplitude-damping channel rather than dephasing.

### Reproducing Wiseman and Milburn, Figure 4.6

A classic picture of a single homodyne trajectory plots its angles on the sphere: the azimuth $\varphi = \arctan(\langle Y\rangle, \langle X\rangle)$ and $\cos\theta = \langle Z\rangle$. Take the resonant, undetuned drive of the textbook figure and its emission channel:

```wl
Ωw = 3.0; Hw = Ωw Y/2; γw = 1.0; Lw = {Sqrt[γw] Jminus}; ρ0w = dm[ketPlus]; δtw = 0.01; tfw = 10.0;
```

Generate one trajectory's Bloch record:

```wl
trajW = BlockRandom[SeedRandom[1]; blochVector /@ \[ScriptCapitalT][ρ0w, Hw, Lw, δtw, tfw]];
```

Plot its azimuth $\varphi$ and $\cos\theta$:

```wl
ListLinePlot[Transpose[{ArcTan[#1, #2], #3} & @@@ trajW],
 PlotLegends -> {"φ = arctan(⟨Y⟩,⟨X⟩)", "cos θ = ⟨Z⟩"}, Frame -> True, ImageSize -> 540,
 AspectRatio -> 1/2, GridLines -> Automatic, FrameLabel -> {"step", "angle"},
 PlotLabel -> "A single homodyne trajectory: azimuth φ and polar cos θ"]
```

The measurement phase is itself a knob: mixing the emitted field with a local oscillator shifted by $\pi/2$ multiplies the monitored operator by $e^{i\pi/2}=i$, so a different quadrature of the field is recorded. Rotate the measured quadrature:

```wl
Lw2 = {Exp[I π/2] Sqrt[γw] Jminus};
```

Rerun with the same seed:

```wl
trajW2 = BlockRandom[SeedRandom[1]; blochVector /@ \[ScriptCapitalT][ρ0w, Hw, Lw2, δtw, tfw]];
```

Plot its angles:

```wl
ListLinePlot[Transpose[{ArcTan[#1, #2], #3} & @@@ trajW2],
 PlotLegends -> {"φ", "cos θ"}, Frame -> True, ImageSize -> 540, AspectRatio -> 1/2,
 GridLines -> Automatic, FrameLabel -> {"step", "angle"},
 PlotLabel -> "Rotating the measured quadrature by π/2 changes the trajectory"]
```

Same qubit, same drive, same seed, but the recorded quadrature is different, and so is the path the state takes: the trajectory is a property of the measurement as much as of the system.

## Part VI: The Measurement Record

So far we have read the conditioned *state* off each run. But Part I was blunt about what an experiment directly records: only the *current*, the readout increment $dR = \sqrt{\eta}\,\mathrm{Tr}[(L+L^\dagger)\rho]\,dt + dW$, which $\mathcal{T}$ has been `Sow`-ing all along. The Bloch vectors and purities of Parts III through V were never measured; they were inferred from a record by the SME. This Part works with the record itself: first its two-point correlation, then its power spectrum, then the honest question the framing forces, whether the inferred single-run trajectories were ever real. That analysis needs a long record, so set a resonantly driven, emitting qubit and a long horizon:

```wl
Ωc = 3.5; Hc = Ωc Y/2; γc = 0.7; Lc = {Sqrt[γc] Jminus}; ρ0c = dm[ket0]; δtc = 0.01; tfc = 4000.0;
```

Generate the trajectory and reap the readout it sows:

```wl
{trajC, {records}} = BlockRandom[SeedRandom[1]; Reap@\[ScriptCapitalT][ρ0c, Hc, Lc, δtc, tfc]];
```

There is only one monitored operator, so there is one current. Pull out its time series:

```wl
current = Re@Transpose[records][[1]];
```

Check its length:

```wl
Length[current]
```

### A Two-Point Correlation of the Raw Record

One step of the record is signal plus noise, $dR = \sqrt{\eta}\,\mathrm{Tr}[(L+L^\dagger)\rho]\,dt + dW$: the signal is of order $dt$, while the shot noise $dW$ has standard deviation $\sqrt{dt}$, so at $dt = 0.01$ each sample is about ten parts noise to one part signal. Confirm the record is noise-dominated sample by sample:

```wl
{Mean[current], StandardDeviation[current], Sqrt[δtc]}
```

The standard deviation is essentially $\sqrt{dt} = 0.1$, so one reading tells you almost nothing. The two-point correlation asks a different question, whether the reading at time $t$ is statistically related to the reading a lag $\tau$ later, $C(\tau) = \langle dR(t)\,dR(t+\tau)\rangle - \langle dR\rangle^2$. It is worth computing for a modest reason: coherent dynamics leave a fingerprint in $C(\tau)$ that survives the shot noise, and $C(\tau)$ is the step toward the power spectrum a laboratory actually reports. It is one statistic among several, not the whole record: being second order and phase-insensitive, it reveals the frequency and decay of the coherent motion but not the state, and it discards everything in the record's higher moments. Define an estimator that, for each lag, correlates the record with a shifted copy of itself and subtracts the squared mean:

```wl
twoPointCorrelation[data_, hmax_, dt_ : None, steps_ : 1] :=
  Module[{n = Length[data], ave = Mean[data], vec1, mm, corr},
   vec1 = data[[1 ;; n - hmax]];
   mm = Table[data[[1 + i ;; n - hmax + i]], {i, 1, hmax, steps}];
   corr = (1/n) (mm . vec1 - ave^2);
   If[dt === None, corr, {dt Range[1, hmax, steps], corr}]];
```

One detail of the estimator governs everything below: its lags are `dt Range[1, hmax, steps]`, starting at one step, $\tau = dt$, and never reaching $\tau = 0$. The shot noise is *white*, so its autocovariance is concentrated at exactly $\tau = 0$ and vanishes at every other lag. Since $\tau = 0$ is never plotted, the shot noise contributes nothing to the estimate, and the raw correlation shows the signal's own structure with no filtering at all. Correlate the raw current out to a lag of eight time units:

```wl
corrRaw = twoPointCorrelation[current, Floor[8/δtc], δtc, 4];
```

Plot it against the lag $\tau$:

```wl
ListLinePlot[Transpose[corrRaw], Frame -> True, ImageSize -> 540, AspectRatio -> 1/2,
 GridLines -> Automatic, PlotRange -> All, FrameLabel -> {"lag τ", "two-point correlation C(τ)"},
 PlotLabel -> "The raw current's correlation oscillates at the Rabi frequency"]
```

The correlation falls from its largest plotted value at $\tau = dt$ and oscillates as $\cos(\Omega_c\tau)$: the first trough sits near $\tau \approx 0.9$ and the first peak near $\tau \approx 1.8$, which is the Rabi period $2\pi/\Omega_c \approx 1.8$. There is no spike at the origin, precisely because the estimator skips $\tau = 0$ where the white shot noise lives. The oscillation amplitude decays on the scale of the emission time $1/\gamma_c \approx 1.4$, and past $\tau \approx 2.5$ the later cycles sink into the estimator's statistical scatter, a few $\times 10^{-5}$ over a finite record. So the raw record already carries the coherent dynamics, but the time-domain correlation fixes the frequency only roughly. The clean read is the spectrum, below.

### Low-Pass Filtering a Shot-Noise Record

A natural next instinct is to denoise the current with a low-pass filter, and it is worth seeing what that does to a shot-noise-dominated record, because the result is easy to misread. Shot noise is *white*: its power is flat in frequency out to the sampling limit. A low-pass of cutoff $\omega_c$ keeps only the frequencies below $\omega_c$, discarding most of the noise power, but it also *colors* what survives, giving the filtered noise a correlation time of order $1/\omega_c$. The zero-lag spike is then no longer a spike: it is smeared into a bump a few $\times 1/\omega_c$ wide. Define a zero-phase Butterworth low-pass (forward then reversed, averaged, so features are not shifted in time):

```wl
butterworthFilter[data_, ωc_, dt_ : 1, order_ : 2] :=
  Module[{model, f, b},
   model = ToDiscreteTimeModel[ButterworthFilterModel[{"Lowpass", order, ωc}], dt];
   f = RecurrenceFilter[model, data];
   b = RecurrenceFilter[model, Reverse[data]];
   (f + Reverse[b])/2];
```

Correlate the filtered current the same way, at cutoff $\omega_c = 50$, so the filter's correlation time is $1/\omega_c = 0.02$, two steps:

```wl
corrSmooth = twoPointCorrelation[butterworthFilter[current, 50, δtc], Floor[8/δtc], δtc, 4];
```

Measure the small-lag values against the raw correlation's largest:

```wl
{corrRaw[[2, 1]], corrSmooth[[2, 1]], corrSmooth[[2, 1]]/Max[corrRaw[[2]]]}
```

The filtered correlation starts near $9 \times 10^{-4}$, about nine times the largest value the raw correlation ever reaches, and it falls back to the oscillation's amplitude within a few lags (below a tenth of its height by $\tau \approx 0.13$). That short-lag feature, a bump a few $\times 1/\omega_c$ wide, is why a plot of `corrSmooth` with `PlotRange -> All` is dominated by the origin and buries the oscillation. It is not the Rabi signal. The test is to run the identical filter on a pure Wiener record of the same length with no qubit at all:

```wl
pureNoise = BlockRandom[SeedRandom[2]; RandomVariate[NormalDistribution[0, Sqrt[δtc]], Length[current]]];
corrNoise = twoPointCorrelation[butterworthFilter[pureNoise, 50, δtc], Floor[8/δtc], δtc, 4];
ListLinePlot[Transpose /@ {corrSmooth, corrNoise}, Frame -> True, ImageSize -> 540,
 AspectRatio -> 1/2, GridLines -> Automatic, PlotRange -> All,
 PlotLegends -> {"filtered current", "filtered pure noise"},
 FrameLabel -> {"lag τ", "two-point correlation C(τ)"},
 PlotLabel -> "The small-lag bump is the filter acting on shot noise, not the qubit"]
```

The two bumps coincide: filtering pure noise reproduces almost all of the feature (about ninety percent of its height), and the pure-noise curve is flat past it, while only the current keeps an oscillation at larger $\tau$. The bump is the filter, not the qubit. To see what the filter does to the actual signal, plot both correlations again with the bump excluded, restricting the lag to $\tau \geq 0.2$:

```wl
pastBump = Select[Transpose[#], First[#] >= 0.2 &] & /@ {corrRaw, corrSmooth};
ListLinePlot[pastBump, PlotLegends -> {"raw current", "Butterworth-filtered"}, Frame -> True,
 ImageSize -> 540, AspectRatio -> 1/2, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"lag τ", "two-point correlation C(τ)"},
 PlotLabel -> "Past the bump the filter tracks the same oscillation, more smoothly"]
```

Past the bump the filtered correlation traces the same $\cos(\Omega_c\tau)$ as the raw one, with the high-frequency scatter removed, so the first Rabi cycle reads more cleanly. That is the honest accounting: the low-pass suppresses jitter at $\tau \gtrsim 0.2$ at the price of a spurious bump at $\tau \lesssim 0.13$, and it does *not* strip the noise and leave the signal clean. The raw correlation needed none of this, since the estimator already skips the lag where the shot noise sits. A low-pass earns its place when the *current itself* must be smoothed, for a display or a feedback loop, not as a way to read the correlation.

### The Power Spectrum: What the Experiment Reports

The frequency-domain view is the clean one. By the Wiener-Khinchin theorem the power spectrum of the current is the Fourier transform of its two-point correlation, and it separates signal from noise where the time-domain correlation cannot: white shot noise, flat in frequency, becomes a constant floor, while the coherent oscillation becomes a peak on top of it, at the Rabi frequency and with a width set by the emission rate. Estimate it by Welch's method, averaging the periodograms of many segments of the record so the floor comes out smooth:

```wl
powerSpectrum[data_, dt_, nseg_] := Module[{seglen, segs, half, psd, freqs},
   seglen = Floor[Length[data]/nseg]; segs = Partition[data, seglen]; half = Floor[seglen/2];
   psd = (dt/seglen) Mean[Abs[Fourier[# - Mean[#], FourierParameters -> {-1, 1}]]^2 & /@ segs];
   freqs = 2 Pi Range[0, seglen - 1]/(seglen dt);
   Transpose[{freqs[[2 ;; half]], psd[[2 ;; half]]}]];
```

Compute it and plot over a band that brackets the drive, marking $\Omega_c$ with a gridline:

```wl
spectrum = Select[powerSpectrum[current, δtc, 120], #[[1]] <= 8 &];
ListLinePlot[spectrum, Frame -> True, PlotRange -> All, GridLines -> {{Ωc}, None},
 ImageSize -> 540, AspectRatio -> 1/2, FrameLabel -> {"frequency ω", "power spectrum S(ω)"},
 PlotLabel -> "The homodyne spectrum peaks at the Rabi frequency"]
```

The spectrum is flat at the shot-noise floor away from the drive and rises to a peak right at the gridline, about a factor of two above the floor here, modest because the measurement is weak and sharper the harder one drives or the longer one averages. Three numbers come straight off it: the peak *position* is the Rabi frequency $\Omega_c$, the peak *width* is the emission rate, and the *floor* is the detector's shot noise. Read the peak position off the data, with no knowledge of the Hamiltonian that produced it:

```wl
MaximalBy[Select[spectrum, 2 <= #[[1]] <= 5 &], Last][[1, 1]]
```

which returns $\Omega_c = 3.5$ to the frequency resolution of the estimate, about $0.19$ here. Fitting the peak's width instead gives the emission rate, which is continuous-measurement spectroscopy, and the point is that the coherent motion the earlier Bloch plots showed is present here too, now read from the record alone, which is all an experiment ever has.

### Are the Earlier Trajectories Real, or Is Only the Record?

This is the question Part I set up and Part VI has to answer. Only the record is measured. The conditional state $\rho_c(t)$, and every $\langle Z\rangle$, purity, and Bloch vector plotted from a single run in Parts III through V, was *inferred* from a record, never read off a meter. So were those single-run plots a fiction?

No, and it is worth being precise about why. The conditional state is the observer's best estimate of the system given the record, and the SME is exactly the machine that computes it, the optimal quantum filter, the quantum analogue of a Kalman filter. An estimate is not a fiction, and this one is tied to measurable quantities three ways.

First, the inferred state predicts the record's own measurable statistics. The mean current is an observable, a number you get by averaging the readout; the inferred trajectory predicts it as $dt$ times the conditional signal $\mathrm{Tr}[(L+L^\dagger)\rho_c]$ averaged along the run. Compare the two:

```wl
Column[{
   Row[{"mean current, measured off the record:    ", Mean[current]}],
   Row[{"δt × mean inferred signal Tr[(L+L†)ρc]:   ",
      δtc Mean[Re[Tr[(Lc[[1]] + ConjugateTranspose[Lc[[1]]]) . #]] & /@ Most[trajC]]}]}]
```

They agree to the shot-noise floor: the current the detector produced is, on average, exactly what the inferred trajectory says it should be, and the spectral peak above is the same statement in frequency, the record's coherent frequency sitting where the inferred state oscillates. Part VI reproducing these record statistics *is* a test of the inferred trajectory, and it passes.

Second, averaging the inferred conditional quantities over many runs gives the unconditional master-equation state, and *that* is directly measurable, by repeating the experiment and averaging: $\langle Z\rangle(t)$ averaged over runs equals $\mathrm{Tr}[\rho(t)Z]$. Parts III and V computed exactly this, the mean of two hundred conditioned trajectories falling on the independent `LindbladSolve` curve to the Monte-Carlo scale $1/\sqrt{K}$. That agreement tests the inference against a quantity an experiment measures without ever reconstructing a single trajectory.

Third, laboratories reconstruct *and verify* the single conditional trajectory. Continuous-measurement optomechanics runs this same filter on the photocurrent in real time: a levitated nanoparticle (Magrini et al., Nature 595, 373 (2021)) and a soft-clamped membrane (Rossi et al., with the companion "Observing and verifying the quantum trajectory of a mechanical resonator", Phys. Rev. Lett. 123, 163601 (2019)) each reconstruct the conditional state from the record by a real-time Kalman filter and check its predictions against the measured data.

So the single-run conditional quantities are not fake. They are inferences, testable through the record statistics they predict and through the ensemble averages that are directly measurable, and they are reconstructed and verified in the laboratory. What they are not is a raw meter reading: they are the best estimate conditioned on one record, and when detection is inefficient (Part VII) even the best estimate is mixed.

### What Else a Record Gives You

The state, the spectrum, and the mean current are a few uses of a record, not all of them. With the whole record in hand, the *future* readings also constrain the state at an earlier time, and combining past and future gives a strictly sharper estimate, the quantum *smoother*. The record, or the state filtered from it, can drive the Hamiltonian in real time to stabilize, cool, or error-correct the qubit, which is quantum feedback, the setting the positivity-preserving scheme was built for, since a filter that returned a non-state would feed the controller nonsense. And the correlation and spectrum are only second order: the higher moments, the full distribution of a time-integrated current, and, for a photon-counting rather than homodyne record, waiting-time distributions and antibunching all come from the same stream. The filter itself is a free choice matched to the question, the Butterworth low-pass above merely denoising, a *matched* filter tuned to detect a specific feature, a *lock-in* pulling out one quadrature at the drive frequency, and a *Wiener* or *Kalman* filter giving the statistically optimal estimate. The record is the raw observable; the state, its spectrum, the system's parameters, and a feedback signal are all things one computes from it.

## Part VII: What the Dimension-General Map Unlocks

The map $\mathcal{T}$ takes plain matrices of any size and carries the detector efficiency honestly, and both facts open doors a qubit-only, unit-efficiency version cannot. Take efficiency first. Real detectors lose some of the signal, $\eta < 1$, so part of the emitted information escapes unrecorded. The conditioned state then no longer purifies: it is a genuine mixture, its Bloch vector shrinks inside the sphere, yet the *ensemble* average is unchanged, because the unconditional master equation never knew there was a detector. Set a driven, emitting qubit:

```wl
Ωe = 2.0; He = Ωe X/2; γe = 1.0; Le = {Sqrt[γe] Jminus}; ρ0e = dm[ketPlus]; δte = 0.01; tfe = 8.0;
```

Run one trajectory at unit efficiency:

```wl
runFull = BlockRandom[SeedRandom[5]; \[ScriptCapitalT][ρ0e, He, Le, δte, tfe]];
```

Run one at half efficiency, on the same seed:

```wl
runHalf = BlockRandom[SeedRandom[5]; \[ScriptCapitalT][ρ0e, He, Le, {0.5}, δte, tfe]];
```

Define the purity $\mathrm{Tr}(\rho^2)$:

```wl
purity[ρ_] := Re@Tr[ρ . ρ];
```

Compare the mean purity of the two runs:

```wl
Column[{
   Row[{"mean purity at η = 1:   ", Mean[purity /@ runFull]}],
   Row[{"mean purity at η = 0.5: ", Mean[purity /@ runHalf]}]}]
```

At unit efficiency the conditioned states stay essentially pure (purity near $1$): a perfect record leaves the observer certain of the state. At $\eta = 0.5$ the mean purity drops well below $1$, because half the information leaks away and the best inferred state is mixed. See it geometrically: draw the sphere with `BlochSpherePlot[]` and overlay both full trajectories, the unit-efficiency run and the half-efficiency run:

```wl
Show[ResourceFunction["BlochSpherePlot"][],
 Graphics3D[{Thick, ColorData[97, 1], Line[blochVector /@ runFull],
   ColorData[97, 2], Line[blochVector /@ runHalf]}], ImageSize -> 420]
```

The full-efficiency curve rides on the surface of the ball; the half-efficiency curve stays inside it, its shrunken radius the visible signature of unrecorded information. Yet both, averaged over many runs, return the same Lindblad state, since the ensemble average is independent of $\eta$. That $\eta$-independence of the average, together with the correct handling of the Ito correction $L_{\rm corr} = \tfrac12 dt\sum_r\eta_r L_r^2$, is exactly the content the efficiency factor carries.

Now dimension. Nothing in $\mathcal{T}$ mentions the number $2$, so a qutrit runs by the same call. Build the spin-1 lowering operator, a soft drive from it, and its emission channel, with the excited $|m=1\rangle$ as the initial state:

```wl
Jm3 = {{0, 0, 0}, {Sqrt[2], 0, 0}, {0, Sqrt[2], 0}}; Sx3 = (Jm3 + ConjugateTranspose[Jm3])/Sqrt[2]; H3 = 1.5 Sx3; L3 = {Sqrt[0.6] Jm3}; ρ03 = dm[{1, 0, 0}]; δt3 = 0.01; tf3 = 8.0;
```

Run one qutrit trajectory at $\eta = 0.4$:

```wl
run3 = BlockRandom[SeedRandom[2]; \[ScriptCapitalT][ρ03, H3, L3, {0.4}, δt3, tf3]];
```

Confirm the states are $3\times 3$ and every one is positive semidefinite:

```wl
Column[{
   Row[{"state dimension: ", Dimensions[First@run3]}],
   Row[{"all physical? ", Tally[PositiveSemidefiniteMatrixQ /@ run3]}]}]
```

The states are $3\times 3$ and physical: the same positivity guarantee holds in any dimension. Confirm the physics too, that the ensemble mean reproduces the $d=3$ master equation. Set the grid:

```wl
tgrid3 = δt3 Range[0, Floor[tf3/δt3]];
```

Solve the $d=3$ master equation with the same `LindbladSolve`:

```wl
lind3 = ResourceFunction["LindbladSolve"][H3, L3, ρ03, {t, 0, tf3}];
```

Run 200 qutrit trajectories:

```wl
ens3 = BlockRandom[SeedRandom[9]; Table[\[ScriptCapitalT][ρ03, H3, L3, {0.4}, δt3, tf3], {200}]];
```

Read the three level populations off a state:

```wl
pop3[ρ_] := Re@Diagonal[ρ];
```

Collect the ensemble-mean populations over time:

```wl
meanPops = Transpose[pop3 /@ Mean[ens3]];
```

Collect the master-equation populations over time:

```wl
lindPops = Transpose[pop3[lind3[#]] & /@ tgrid3];
```

Plot the ensemble means (solid) against the master equation (dashed):

```wl
ListLinePlot[
 Join[Transpose[{tgrid3, #}] & /@ meanPops, Transpose[{tgrid3, #}] & /@ lindPops],
 PlotStyle -> Join[{ColorData[97, 1], ColorData[97, 2], ColorData[97, 3]},
   ConstantArray[Directive[Dashed, Black], 3]],
 Frame -> True, ImageSize -> 540, AspectRatio -> 1/2, GridLines -> Automatic, PlotRange -> All,
 PlotLegends -> SwatchLegend[{ColorData[97, 1], ColorData[97, 2], ColorData[97, 3], Black},
   {"|1⟩ mean", "|0⟩ mean", "|-1⟩ mean", "Lindblad (all three)"}],
 FrameLabel -> {"time t", "level populations"},
 PlotLabel -> "A qutrit at η = 0.4: 200 trajectories average to the d = 3 master equation"]
```

The three level populations of the ensemble mean land on the three master-equation curves. This is the payoff of building the integrator on plain, dimension-general matrices with an explicit efficiency: inefficient detection and higher-dimensional systems are the same function call, not a rewrite, and both stay physical and correct in the mean.

Before we close, let us summarize the most important points we have computed and checked:

- The stochastic master equation is the Lindblad master equation plus one innovation term $\sqrt{\eta_r}\,\mathcal{H}[L_r]\rho\,dW_r$ per monitored channel, and the detector emits $dR_r = \sqrt{\eta_r}\,\mathrm{Tr}[(L_r+L_r^\dagger)\rho]\,dt + dW_r$.
- The positivity-preserving update writes $\rho_{n+1}$ as a normalized sum of $A\rho A^\dagger$ terms, so every conditioned state is positive semidefinite at any step size, which we confirmed with `PositiveSemidefiniteMatrixQ` in $d=2$ and $d=3$.
- Averaging many trajectories reproduces the master equation to the Monte-Carlo scale $1/\sqrt{K}$, for $Z$ dephasing, for $J_-$ emission, and for the qutrit.
- Strong $Z$ monitoring freezes the qubit (the quantum Zeno effect, decay at $\sim\Omega^2/\gamma$), while weak monitoring lets Rabi oscillations survive with a slow $\sim\gamma$ decay.
- Only the record is directly measured; the conditional state and every Bloch vector or purity plotted from a single run is inferred from it, the SME acting as the optimal filter.
- The power spectrum of the homodyne current is flat at the shot-noise floor and peaks at the Rabi frequency (its width the emission rate), which is how the coherent dynamics are read from the record; the small-lag bump in the Butterworth-filtered time-domain correlation is filtered shot noise, not signal, as filtering pure noise confirms.
- A record supports more than a spectrum: state filtering (what the Bloch plots are), parameter estimation, smoothing, and feedback, with the choice of filter set by the question.
- At $\eta<1$ the conditioned state is mixed and its Bloch vector sits inside the sphere, while the ensemble average stays $\eta$-independent and equal to the master-equation state.

### Where This Leaves Us (and What Comes Next)

You now have a complete, computation-first toolkit for continuously measured quantum systems: the stochastic master equation and its readout, a positivity-preserving map $\mathcal{T}$ over plain matrices of any dimension, the Function Repository's `LindbladSolve` as the master-equation reference and `BlochSpherePlot` for the sphere, exercised on dephasing, spontaneous emission, the homodyne current's spectrum, inefficient detection, and a qutrit. The single most consequential parameter is the ratio of measurement rate to internal dynamics: turn it up and the record pins the state and freezes coherent motion; turn it down and the state follows its Hamiltonian, barely disturbed. The natural continuations are quantum feedback, using the conditioned state or the record to drive the Hamiltonian in real time (the setting the scheme was built for), and larger monitored systems, two qubits and beyond, where the same dimension-general call runs unchanged and the trajectories carry entanglement conditioned on a shared record.

### Acknowledgment

The qubit examples follow the continuous-measurement lectures of Gabriel T. Landi and the homodyne-trajectory pictures of Wiseman and Milburn's *Quantum Measurement and Control*; the integrator implements the positivity-preserving scheme cited in Part II. The framing throughout, one conditioned run against the average of many, is the trajectory unravelling that makes an abstract master equation into something an experiment records.
