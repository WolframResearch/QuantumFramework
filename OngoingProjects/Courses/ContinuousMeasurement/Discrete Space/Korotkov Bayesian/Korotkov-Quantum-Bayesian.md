---
Template: Default
---

# Watching a Qubit: Korotkov's Quantum-Bayesian Formalism for Continuous Measurement

**A weak continuous measurement never collapses a qubit outright. It trickles information out through a noisy detector current, and the observer folds each noisy reading back into the state with Bayes' rule. This essay builds Korotkov's quantum-Bayesian update from primitives and runs it: we watch a frozen qubit localize by accumulating likelihood, watch a pure state stay exactly pure along a single record while the ensemble dephases, drive the qubit with a Rabi field and interleave the measurement as a symmetric split, cap the purity with a non-ideal detector, and finally recover the closed-form output spectrum from simulated records, where the coherent line can never rise more than four times above the detector noise floor. A last section closes the loop with the sibling Stratonovich stochastic differential equation and exposes the Ito trap that a careless integrator falls into.**

Mads Bahrami (last updated: August 18, 2026)

## Setting the Stage: How This Essay Flows

This essay is a computation-first tour of a single qubit under one symmetric broadband detector, in the Markov limit: how a Gaussian likelihood turns a noisy current into a state update, why the update is exactly Bayes' rule on the populations with the coherence riding along, how a Hamiltonian is interleaved with the measurement, what non-ideality costs, and why the oscillation line in the output spectrum is bounded by a factor of four. Every claim is computed on the spot from base Wolfram Language, so nothing here rests on a formula you cannot immediately test and modify.

In other words, I have tried to build a small laboratory for continuous measurement of a qubit. I strongly believe in a computation-first narrative for learning: in a sense, if I cannot compute it, I cannot claim to understand it. The rhythm throughout is concept, then computation, then interpretation.

The environment you see is a live Wolfram notebook. Evaluate the cells from top to bottom; the toolkit defined in the next section is used by every section that follows, so those dependencies matter. Some cells run Monte-Carlo ensembles and take a little while. My suggestion is to focus on the output and its meaning first, then unpack the input code. You are not locked into any of it: change the rates, the step, the seed, and rerun your own experiments.

Every equation below is anchored to the grounded physics sheet in this folder (`korotkov_physics_sheet.md`), which transcribes Korotkov's own papers with line references; I cite it as [sheet X] rather than re-derive from memory. This is the Bayesian-update door into the same physics that the sibling essays treat as a stochastic differential equation: `NDSolve vs Ito/` integrates the Stratonovich Bloch-Langevin equation for the state and readout, and `Manual StepLike Ito/` runs Rouchon trajectories. Section 7 closes the loop by matching that Stratonovich equation numerically.

Let's start!

## Notation and Conventions: The Detector, the Record, and the Rates

We work in units where the reduced Planck constant is one, $\hbar=1$. The qubit is written in the measured (localized) basis $\{|1\rangle,|0\rangle\}$ that the detector couples to, with $|1\rangle=\{1,0\}$ and $|0\rangle=\{0,1\}$, so the density matrix is a $2\times 2$ matrix $\rho$ and the single monitored coordinate is the population difference $z=\rho_{11}-\rho_{00}$ [sheet K]. A qubit density matrix carries three real numbers, the Bloch vector; the detector reads out only the projection $z$, and the other two ride along in the coherence.

The detector emits a current whose short-time average over a window $\tau$ is $\bar I$. Given a definite state $|1\rangle$ or $|0\rangle$ the reading is Gaussian about $I_1$ or $I_0$, and in general it is a linear readout of $z$ plus white noise,

$$ I(t)=I_c+\frac{\Delta I}{2}\,z(t)+\xi(t), \qquad \Delta I=I_1-I_0, \qquad I_c=\frac{I_1+I_0}{2}, $$

with $\xi$ a white noise of single-sided spectral density $S$ [sheet K]. Averaging that current over the window $\tau$ leaves a Gaussian of variance $D=S/(2\tau)$ about the state-conditioned mean [sheet E]. Two rates organize everything: the measurement (information-acquisition) rate and the total ensemble dephasing rate,

$$ \Gamma_m=\frac{(\Delta I)^2}{4S}, \qquad \Gamma_d \ge \Gamma_m, \qquad \gamma=\Gamma_d-\Gamma_m\ge 0, \qquad \eta=\frac{\Gamma_m}{\Gamma_d}\in(0,1] , $$

where $\gamma$ is any extra dephasing the detector dissipates internally before readout, and the ideality $\eta$ is one exactly for a quantum-limited detector such as a symmetric quantum point contact [sheet B, sheet J]. A symmetric detector has zero output-backaction correlation, so there is no unitary phase kick; the back-action is purely informational, along the meridians of the Bloch sphere [sheet K]. When a Hamiltonian is present it is the symmetric-qubit drive $H_{qb}=\tfrac{\Omega_R}{2}\sigma_x$ with Rabi frequency $\Omega_R$ [sheet A].

Fix the detector's two output levels and its noise once. We set $I_1=+1$, $I_0=-1$ (so $\Delta I=2$) and $I_c=0$; the noise density then follows from whatever measurement rate we choose, $S=(\Delta I)^2/(4\Gamma_m)$:

```wl
ClearAll[dI, i1, i0, ic];
dI = 2.; i1 = 1.; i0 = -1.; ic = 0.;
{i1, i0, dI}
```

As one can see, the two conditioned means sit symmetrically about zero, so a reading's sign is the only hint of which state produced it, and $\Delta I=2$ is the full separation the noise must beat.

The chance that a state $|j\rangle$ produced a given reading is a Gaussian likelihood; define it as a function of the window variance $D$:

```wl
ClearAll[gaussLike];
gaussLike[mu_, var_, ibar_] := Exp[-(ibar - mu)^2/(2 var)]/Sqrt[2 Pi var]
```

This is the single probabilistic object the whole formalism turns on: everything downstream is Bayes' rule applied to these two likelihoods.

Define the window variance $D=S/(2\tau)$ and the honest detector draw, sampling $\bar I$ from the Gaussian mixture $\rho_{11}\,\mathcal N(I_1,D)+\rho_{00}\,\mathcal N(I_0,D)$ that the current state predicts [sheet E]:

```wl
ClearAll[avgVar, drawCurrent];
avgVar[gm_, tau_] := (dI^2/(4 gm))/(2 tau);
drawCurrent[rho_, var_] :=
  RandomChoice[Clip[Re@{rho[[1, 1]], rho[[2, 2]]}, {0., 1.}] -> {i1, i0}] +
    RandomVariate[NormalDistribution[0, Sqrt[var]]]
```

In words, the detector first picks which bell to ring with the current populations as weights, then buries the level under Gaussian noise of variance $D$; the window $\tau$ is the only knob that sharpens it.

Now the heart of the formalism. The exact update over a step, in its cleanest form, keeps the populations Bayesian and lets each coherence carry the geometric mean of the two likelihoods, with an optional pure-dephasing factor for a non-ideal detector [sheet D, sheet J, sheet K]:

$$ \rho_{ij}\ \longrightarrow\ \rho_{ij}\,\frac{\sqrt{P_i\,P_j}}{\sum_k \rho_{kk}\,P_k}\ \times\ e^{-\gamma\tau}\ \text{(off-diagonal only)}, \qquad P_j=\text{gaussLike}(I_j,D,\bar I) . $$

Encode that update as one element-wise operation on the matrix:

```wl
ClearAll[bayesUpdate];
bayesUpdate[rho_, ibar_, var_, gt_] :=
  With[{p1 = gaussLike[i1, var, ibar], p0 = gaussLike[i0, var, ibar]},
   rho {{p1, Sqrt[p1 p0] Exp[-gt]}, {Sqrt[p1 p0] Exp[-gt], p0}}/
    (Re@rho[[1, 1]] p1 + Re@rho[[2, 2]] p0)]
```

The diagonal entries pick up $P_j$ and renormalize, which is textbook Gaussian Bayes on the populations; the off-diagonal entries pick up $\sqrt{P_1P_0}$, so the coherence is pinned to the populations rather than measured directly.

Confirm that this update is trace-preserving for any reading, symbolically, so that $\rho$ stays a density matrix:

```wl
With[{r = {{r11, r10}, {Conjugate[r10], 1 - r11}}},
 FullSimplify[Tr[bayesUpdate[r, ibar, var, gt]] == 1,
  Assumptions -> {0 < r11 < 1, var > 0, ibar \[Element] Reals, gt >= 0}]]
```

As expected, the two updated populations always sum to one: the normalizing denominator is exactly the predictive probability of the reading, so probability is conserved by construction.

A single step of the full dynamics interleaves a Hamiltonian half-rotation, the Bayesian kick, and another half-rotation, a symmetric (Strang) split that is second-order accurate in the step [sheet C, sheet E]; for a frozen qubit the half-rotation is just the identity. Define the step, together with the readouts we will plot:

```wl
ClearAll[measureStep, blochVec, purity, zComp];
measureStep[uHalf_, var_, gt_][rho_] :=
  With[{r1 = uHalf . rho . ConjugateTranspose[uHalf]},
   uHalf . bayesUpdate[r1, drawCurrent[r1, var], var, gt] . ConjugateTranspose[uHalf]];
blochVec[rho_] := Re@Table[Tr[rho . PauliMatrix[j]], {j, 3}];
purity[rho_] := Re@Tr[rho . rho];
zComp[rho_] := Re[rho[[1, 1]] - rho[[2, 2]]]
```

That is the entire integrator: a Gaussian draw, an element-wise Bayes kick, and two matrix products for the drive. For the output spectrum in Section 6 we will also want the record itself, so define a run that emits the stream of readings while it measures a driven qubit:

```wl
ClearAll[recordCurrents];
recordCurrents[gm_, tau_, gt_, om_, nst_, seed_] :=
  BlockRandom[SeedRandom[seed];
   Module[{uH = MatrixExp[-I (om/2) PauliMatrix[1] (tau/2)],
     v = avgVar[gm, tau], rho = {{0.5, 0.5}, {0.5, 0.5}}, r1, ib},
    Table[r1 = uH . rho . ConjugateTranspose[uH]; ib = drawCurrent[r1, v];
      rho = uH . bayesUpdate[r1, ib, v, gt] . ConjugateTranspose[uH]; ib, {nst}]]]
```

We will return to this record in Section 6; for now it is enough that each measurement produces a number, and the state and the record share the very same noise.

## 1. The Detector's Two Bells: Why One Look Cannot Separate Them

A weak measurement is weak because a single short look is ambiguous. For a state that is definitely $|1\rangle$ or definitely $|0\rangle$ the averaged current is a Gaussian centered at $I_1$ or $I_0$ with spread $\sqrt D=\sqrt{S/(2\tau)}$; the two bells overlap heavily until the window $\tau$ is long enough to pull them apart, and the crossover scale is the measurement time $\tau_m=2S/(\Delta I)^2$ [sheet B].

Plot the two likelihoods $P_1$ and $P_0$ against the reading $\bar I$ for windows short, comparable to, and long relative to $\tau_m$:

```wl
With[{gm = 0.5},
 With[{taum = 2 (dI^2/(4 gm))/dI^2},
  Row[Table[
    Plot[{gaussLike[i1, avgVar[gm, tau], ibar], gaussLike[i0, avgVar[gm, tau], ibar]},
     {ibar, -4, 4}, PlotRange -> All, Filling -> Axis, Frame -> True,
     GridLines -> Automatic, FrameLabel -> {"averaged current", "likelihood"},
     PlotLabel -> "window \[Tau] = " <> ToString[tau/taum] <> " \[Tau]m",
     ImageSize -> 240], {tau, taum {1/4, 1, 6}}]]]]
```

Notice that at a quarter of the measurement time the two bells sit almost on top of each other, so a single reading barely favors either state; only as $\tau$ grows past $\tau_m$ do the peaks separate and a reading starts to mean something. Information is nothing but this shrinking overlap.

## 2. Populations Follow Bayes: Localization Without a Hamiltonian

Switch off the Hamiltonian and the record can only report which state we are in, never rotate us between them. Start from the fully uncertain equal mixture and feed each reading through the update: the populations move by the likelihood ratio and drift toward one pole, a noisy random walk in $z$ rather than a jump [sheet E]. The pole that wins is random, with probability equal to the initial population, which is the Born rule emerging from accumulated evidence.

Run one frozen-qubit trajectory from $\rho=\mathrm{diag}(\tfrac12,\tfrac12)$ and watch $z(t)$:

```wl
With[{gm = 1., tau = 0.02, nst = 1000},
 With[{traj = BlockRandom[SeedRandom[17];
      NestList[measureStep[IdentityMatrix[2], avgVar[gm, tau], 0.],
       {{0.5, 0.5}, {0.5, 0.5}}, nst]]},
  ListLinePlot[Transpose[{Range[0, nst] tau, zComp /@ traj}],
   PlotRange -> {-1.05, 1.05}, Frame -> True, GridLines -> Automatic,
   FrameLabel -> {"time (units of 1/\[CapitalGamma]m)", "z = \[Rho]11 - \[Rho]00"},
   PlotLabel -> "one frozen-qubit record: localization by Bayes"]]]
```

The population difference wanders under the noisy evidence and then commits to a pole and stays there; the "collapse" is the accumulation of log-likelihood, not a discontinuous jump.

## 3. The Purity Ride-Along: A Pure State Stays Pure in a Single Run

Here is the fact that makes an ideal detector special. The coherence is not measured; it rides on the populations through the geometric mean $\sqrt{P_1P_0}$, so the ratio $|\rho_{10}|^2/(\rho_{11}\rho_{00})$ is left untouched by every update. A state that starts pure therefore stays exactly pure along the whole record, never leaving the surface of the Bloch sphere, even as its populations localize [sheet D, sheet C].

First confirm the invariance symbolically: check that the ride-along ratio $|\rho_{10}|^2/(\rho_{11}\rho_{00})$ is unchanged by an ideal update:

```wl
With[{r = {{r11, r10}, {Conjugate[r10], 1 - r11}}},
 With[{u = bayesUpdate[r, ibar, var, 0]},
  FullSimplify[
   (u[[1, 2]] u[[2, 1]])/(u[[1, 1]] u[[2, 2]]) == (r10 Conjugate[r10])/(r11 (1 - r11)),
   Assumptions -> {0 < r11 < 1, var > 0, ibar \[Element] Reals}]]]
```

The ratio is algebraically conserved, so purity is a constant of the ideal single-run motion. Watch it hold along a record that starts pure on the equator, $\rho=|{+}\rangle\langle{+}|$, while $z$ localizes:

```wl
With[{gm = 1., tau = 0.02, nst = 1000},
 With[{traj = BlockRandom[SeedRandom[8];
      NestList[measureStep[IdentityMatrix[2], avgVar[gm, tau], 0.],
       {{0.5, 0.5}, {0.5, 0.5}}, nst]]},
  ListLinePlot[{Transpose[{Range[0, nst] tau, purity /@ traj}],
    Transpose[{Range[0, nst] tau, zComp /@ traj}]}, PlotRange -> {-1.05, 1.05},
   Frame -> True, GridLines -> Automatic, FrameLabel -> {"time (units of 1/\[CapitalGamma]m)", "value"},
   PlotLegends -> {"purity Tr[\[Rho]^2]", "z"},
   PlotLabel -> "ideal detector: purity pinned while z localizes"]]]
```

As you may have noticed, the purity sits flat at one while $z$ wanders to a pole: information gain alone does not decohere a single run. The ensemble decoheres only because different runs localize to different poles, which is the subject of Section 7.

## 4. Rabi Under Measurement: A Bayes Kick Between Half-Rotations

Now switch on the drive $H_{qb}=\tfrac{\Omega_R}{2}\sigma_x$ and interleave it with the measurement by the symmetric split: a half rotation, the Bayesian kick, another half rotation [sheet C, sheet E]. The measurement pulls the state toward the $z$-poles while the field rotates it, so a single trajectory is a noisy Rabi oscillation whose phase slowly diffuses; a run that starts fully mixed at the center of the Bloch ball sharpens into that oscillation as the record accrues.

Run one noisy-Rabi trajectory from $|1\rangle$ and, beside it, a run that starts fully mixed and purifies:

```wl
With[{gm = 0.1, tau = 0.05, nst = 500, om = 1.},
 With[{uH = MatrixExp[-I (om/2) PauliMatrix[1] (tau/2)], tt = Range[0, nst] tau},
  With[{pure = blochVec /@ BlockRandom[SeedRandom[5];
       NestList[measureStep[uH, avgVar[gm, tau], 0.], {{1., 0.}, {0., 0.}}, nst]],
     mixed = blochVec /@ BlockRandom[SeedRandom[8];
       NestList[measureStep[uH, avgVar[gm, tau], 0.], {{0.5, 0.}, {0., 0.5}}, nst]]},
   Row[{ListLinePlot[Transpose[{tt, #}] & /@ Transpose[pure],
      PlotRange -> {-1.05, 1.05}, Frame -> True, GridLines -> Automatic,
      PlotLegends -> {"x", "y", "z"}, FrameLabel -> {"time (units of 1/\[CapitalGamma]m)", "Bloch vector"},
      PlotLabel -> "pure start from |1>: noisy Rabi", ImageSize -> 250],
     ListLinePlot[Transpose[{tt, #}] & /@ Transpose[mixed],
      PlotRange -> {-1.05, 1.05}, Frame -> True, GridLines -> Automatic,
      PlotLegends -> {"x", "y", "z"}, FrameLabel -> {"time (units of 1/\[CapitalGamma]m)", "Bloch vector"},
      PlotLabel -> "mixed start: purifies into oscillation", ImageSize -> 250]}]]]]
```

The pure run oscillates and wobbles but keeps its full Bloch length, while the mixed run starts at the origin and swells outward into the same phase-locked, slowly diffusing oscillation: the record both drives and sharpens the state at once.

## 5. A Non-Ideal Detector: Ideality Caps the Purity

A realistic detector loses some of the information it gathers before it reaches the output, and that lost information appears as an extra pure-dephasing factor $e^{-\gamma\tau}$ on the coherence, with $\gamma=\Gamma_d(1-\eta)$ [sheet J]. The measurement keeps localizing at the same informational rate, but the coherence is now bled faster than the populations sharpen, so the state settles short of the surface at a steady purity set entirely by the ideality $\eta$.

Overlay the ensemble-averaged purity for a driven qubit measured by detectors of ideality $\eta=1,\,0.7,\,0.4$, the extra dephasing $\gamma\tau=\Gamma_m(1/\eta-1)\tau$ per step:

```wl
With[{gm = 0.1, tau = 0.05, nst = 500, om = 1., ntr = 80},
 With[{uH = MatrixExp[-I (om/2) PauliMatrix[1] (tau/2)]},
  ListLinePlot[
   Table[Transpose[{Range[0, nst] tau,
      Mean[Table[purity /@ BlockRandom[SeedRandom[Round[1000 eta] + k];
          NestList[measureStep[uH, avgVar[gm, tau], gm (1/eta - 1) tau], {{1., 0.}, {0., 0.}}, nst]],
        {k, ntr}]]}], {eta, {1., 0.7, 0.4}}],
   PlotRange -> {0.4, 1.02}, Frame -> True, GridLines -> Automatic,
   FrameLabel -> {"time (units of 1/\[CapitalGamma]m)", "mean purity Tr[\[Rho]^2]"},
   PlotLegends -> {"\[Eta] = 1", "\[Eta] = 0.7", "\[Eta] = 0.4"},
   PlotLabel -> "ideality caps the steady purity"]]]
```

As expected, only the ideal detector holds the state pure; below $\eta=1$ the purity climbs to a ceiling that drops as $\eta$ drops. The extra dephasing $\gamma$ is potential information dissipated inside the detector, and only $\eta=1$ lets you monitor the wavefunction itself.

## 6. The Output Spectrum and the Factor of Four

Now the payoff. When the drive keeps the qubit oscillating, the Rabi motion modulates the very current the detector reports, so the record carries a coherent line at $\Omega_R$ riding on the white noise floor $S$. Korotkov and Averin derived the line in closed form, and its height above the floor is bounded: at most four times the noise, reached only for an ideal detector, and degrading as $4\eta$ [sheet I, sheet K].

$$ S_I(\omega)=S+4\eta\,S\,\frac{\Omega_R^2\,\Gamma^2}{(\omega^2-\Omega_R^2)^2+\Gamma^2\omega^2}, \qquad \Gamma=\Gamma_d=\frac{\Gamma_m}{\eta}. $$

Half of that four is the ordinary spectrum of $z(t)$; the other half is a purely nonclassical correlation between the detector noise $\xi$ and the back-action it exerts on $z$, which no classical harmonic signal can produce [sheet I]. The record we simulate carries both halves automatically, because the same reading that enters the spectrum also drives the Bayesian kick.

Simulate long records in the good-oscillation regime $\Omega_R\gg\Gamma_m$, average the periodogram over an ensemble for three idealities, and overlay the closed form (dashed):

```wl
With[{gm = 0.05, om = 1., tau = 0.05, nst = 14000, ntr = 200, etas = {1., 0.5, 0.25}},
 Module[{gd, sClosed, spec, specs},
  gd[eta_] := gm/eta;
  sClosed[eta_, w_] := 1 + 4 eta om^2 gd[eta]^2/((w^2 - om^2)^2 + gd[eta]^2 w^2);
  spec[recs_] := Module[{n = Length[recs[[1]]], half, ws, p, fl},
    half = Quotient[n, 2]; ws = Table[2. Pi (s - 1)/(n tau), {s, 2, half}];
    p = Mean[(Abs[Fourier[#]]^2)[[2 ;; half]] & /@ recs];
    fl = Median[Pick[p, UnitStep[ws - 2.5], 1]];
    Transpose[{ws[[2 ;; -2]], MovingAverage[p/fl, 3]}]];
  LaunchKernels[]; DistributeDefinitions["Global`"];
  specs = Table[
    spec[ParallelTable[recordCurrents[gm, tau, (gd[eta] - gm) tau, om, nst, 7000 + Round[1000 eta] + k],
      {k, ntr}]], {eta, etas}];
  Show[
   ListLinePlot[specs, PlotRange -> {{0, 2.2}, {0, 6}}, Frame -> True, GridLines -> Automatic,
    FrameLabel -> {"angular frequency \[Omega]", "spectral density / noise floor"},
    PlotLegends -> {"simulated \[Eta] = 1", "simulated \[Eta] = 0.5", "simulated \[Eta] = 0.25"},
    PlotLabel -> "output spectrum: the coherent line and the factor of four"],
   ListLinePlot[Table[Table[{w, sClosed[eta, w]}, {w, 0.2, 2.2, 0.002}], {eta, etas}],
    PlotStyle -> Directive[Gray, Thick, Dashed]], ImageSize -> 540]]]
```

Each simulated line sits on the closed-form Lorentzian, and the coherent peak climbs toward four times the floor for the ideal detector and falls off as $4\eta$: the line is bounded because half of its area is the noise-backaction correlation, which saturates only when every bit of dephasing is informational.

## 7. Cross-Check: The Stratonovich SDE and the Ito Trap

The same physics is a Bloch-Langevin stochastic differential equation in the sibling essay `NDSolve vs Ito/`. The Bayesian update and that Stratonovich equation must agree as $\tau\to 0$, and comparing them exposes the classic Ito trap. The ensemble coherence decays at the total rate $\Gamma_d$; in the Ito form of the coherence equation that decay is the term $+(\Delta I)^2/4S$, which is not physical decoherence at all but the Stratonovich-to-Ito drift correction [sheet C]. Hand the Stratonovich right-hand side to a plain Ito-Euler stepper without that term and the coherence decays at the wrong rate.

Overlay the ensemble coherence $\langle x\rangle(t)$ from three routes started on the equator: the Bayesian update, the correct Ito stepper with the drift term, and the naive stepper missing it, against the exact decay $e^{-\Gamma_m t}$ (dashed):

```wl
With[{gm = 1., tau = 0.005, nst = 300, ntr = 1500},
 Module[{v = avgVar[gm, tau], bayes, ito, tt = Range[0, nst] tau, run},
  bayes[rho_] := measureStep[IdentityMatrix[2], v, 0.][rho];
  ito[drift_][{x_, z_}] := With[{xb = RandomVariate[NormalDistribution[0, Sqrt[v]]]},
    {x - z (2 gm) x xb tau - drift gm x tau, z + (1 - z^2) (2 gm) xb tau}];
  run[stepper_, init_, read_, seed_] := Mean[Table[read /@
      BlockRandom[SeedRandom[seed + k]; NestList[stepper, init, nst]], {k, ntr}]];
  ListLinePlot[{
    Transpose[{tt, run[bayes, {{0.5, 0.5}, {0.5, 0.5}}, (2 Re@#[[1, 2]] &), 1]}],
    Transpose[{tt, run[ito[1], {1., 0.}, First, 2]}],
    Transpose[{tt, run[ito[0], {1., 0.}, First, 3]}]},
   PlotRange -> {0, 1.02}, Frame -> True, GridLines -> Automatic,
   FrameLabel -> {"time (units of 1/\[CapitalGamma]m)", "ensemble coherence \[LeftAngleBracket]x\[RightAngleBracket]"},
   PlotLegends -> {"Bayesian update", "Ito with +(\[CapitalDelta]I)^2/4S drift", "naive Ito, drift dropped"},
   PlotLabel -> "the Ito trap: the drift correction IS the dephasing",
   Epilog -> {Black, Dashed, Line[Transpose[{tt, Exp[-gm tt]}]]}]]]
```

This confirms that the Bayesian update and the correct Ito stepper both fall onto the exact $e^{-\Gamma_m t}$ curve, while the naive stepper never dephases at all: the missing $+(\Delta I)^2/4S$ term was the entire ensemble decoherence. The Bayesian form has no calculus ambiguity to begin with, because it is discrete Gaussian Bayes between unitary half-steps; the stochastic differential equation reproduces it only when the Stratonovich-to-Ito drift is kept.

## Where This Leaves Us (and What Comes Next)

You now have a complete, computation-first toolkit for one qubit under one symmetric broadband detector: a Gaussian likelihood, a Bayes kick on the populations with the coherence riding along, a symmetric split that adds a Hamiltonian, a pure-dephasing factor for non-ideality, and an ensemble of records whose spectrum recovers the Korotkov-Averin line and its factor-of-four bound. Before moving on, the points that were computed and verified along the way:

- The detector separates the two states only gradually; information is the shrinking overlap of two Gaussians, on the scale $\tau_m=2S/(\Delta I)^2$.
- With no Hamiltonian the populations follow Gaussian Bayes and localize to a pole with the initial-population probability, the Born rule as accumulated evidence.
- An ideal detector conserves $|\rho_{10}|^2/(\rho_{11}\rho_{00})$ exactly, so a pure state stays pure along a single record; the ensemble dephases only because runs localize to different poles.
- A symmetric split interleaves the Rabi drive with the Bayes kick, turning a mixed state into a phase-locked, slowly diffusing oscillation.
- Non-ideality is a pure-dephasing factor $e^{-\gamma\tau}$ with $\gamma=\Gamma_d(1-\eta)$; the steady purity is capped by $\eta$, and only $\eta=1$ monitors the wavefunction.
- The output spectrum carries a coherent line at $\Omega_R$ whose height above the floor is $4\eta$, at most four, half of it the nonclassical noise-backaction correlation.
- The Bayesian update and the Stratonovich SDE agree as $\tau\to 0$; a naive Ito stepper that drops the $(\Delta I)^2/4S$ drift loses the entire ensemble dephasing.

The natural continuations are the pieces this essay deliberately left out: a non-ideal detector with output-backaction correlation, a finite detector bandwidth where the qubit and the cavity field must be tracked together, measurement feedback to lock the oscillation, and two-qubit measurement-induced entanglement. Each is a small extension of the same update, and each is a door the sibling essays begin to open.
