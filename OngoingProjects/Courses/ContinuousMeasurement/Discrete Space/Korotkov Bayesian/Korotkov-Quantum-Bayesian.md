---
Template: Default
---

# Watching a Qubit: The Quantum-Bayesian Approach to Continuous Measurement, from One Qubit to Entanglement

**A weak continuous measurement never collapses a qubit outright. It trickles information out through a noisy detector current, and the observer folds each noisy reading back into the state with Bayes' rule. This essay builds the quantum-Bayesian update from primitives, runs it, and then generalizes it twice. Part I is one qubit under one detector: a frozen qubit localizing by accumulating likelihood, a pure state staying exactly pure along a single record while the ensemble dephases, a Rabi drive interleaved as a symmetric split, a non-ideal detector capping the purity, the output spectrum whose coherent line can never rise more than four times above the noise floor, and the phase back-action a general detector adds, ending on the Kraus operator both generalizations grow from. Part II points a second detector at a different axis, so the two Bayes kicks stop commuting and the state stops collapsing and starts diffusing on the Bloch sphere. Part III gives one detector two qubits and a blind spot, so a measurement alone drives them into a Bell state. Each generalization is carried out against a real superconducting-qubit experiment: the single-qubit trajectories of Murch and Siddiqi, the simultaneous non-commuting measurement of Hacohen-Gourgy and co-workers, and the remote measurement-induced entanglement of Roch and co-workers.**

Mads Bahrami (last updated: August 18, 2026)

## Notation and Conventions: The Detector, the Record, and the Rates

We work in units where the reduced Planck constant is one, $\hbar=1$. The qubit is written in the measured (localized) basis $\{|0\rangle,|1\rangle\}$ that the detector couples to, with $|0\rangle=\{1,0\}$ and $|1\rangle=\{0,1\}$, so the density matrix is a $2\times 2$ matrix $\rho$ and the single monitored coordinate is the population difference $z=\rho_{00}-\rho_{11}=\langle\sigma_z\rangle$. A qubit density matrix carries three real numbers, the Bloch vector; the detector reads out only the projection $z$, and the other two ride along in the coherence.

The detector emits a current whose short-time average over a window $\tau$ is $\bar I$. Its two calibrated outputs are $I_0$ for a qubit definitely in $|0\rangle$ and $I_1$ for one definitely in $|1\rangle$. For a general state, a linear detector therefore has the population-weighted mean

$$ \langle I(t)\rangle=\rho_{00}I_0+\rho_{11}I_1. $$

Because $\rho_{00}+\rho_{11}=1$ and $z=\rho_{00}-\rho_{11}$, the populations are $\rho_{00}=(1+z)/2$ and $\rho_{11}=(1-z)/2$. Substitution separates the mean into a state-independent midpoint and a qubit-dependent signal,

$$ \langle I(t)\rangle=\frac{I_0+I_1}{2}+\frac{I_0-I_1}{2}\,z(t). $$

The actual record fluctuates around this mean. Defining the midpoint $I_c$ and detector response $\Delta I$, and writing the zero-mean white noise as $\xi(t)$, gives

$$ I(t)=I_c+\frac{\Delta I}{2}\,z(t)+\xi(t), \qquad \Delta I=I_0-I_1, \qquad I_c=\frac{I_0+I_1}{2}. $$

Here $\xi(t)$ is zero-mean white detector noise. In the single-sided convention used here,

$$
\langle\xi(t)\rangle=0,
\qquad
\langle\xi(t)\xi(t')\rangle=\frac{S}{2}\,\delta(t-t').
$$

Ideal white noise has no finite pointwise variance; only a record averaged over a nonzero time window is an ordinary random variable. Over the same window $\tau$ used to define $\bar I$, define

$$
\bar\xi_\tau(t)=\frac{1}{\tau}\int_t^{t+\tau}\xi(t')\,dt'.
$$

The averaged noise is Gaussian, with

$$
\langle\bar\xi_\tau(t)\rangle=0,
\qquad
\mathrm{Var}(\bar\xi_\tau(t))=\frac{S}{2\tau}.
$$

Thus a longer averaging window suppresses the detector noise: the variance decreases as $S/(2\tau)$ and the standard deviation as $\sqrt{S/(2\tau)}$. The limits make the construction explicit: $z=+1$ returns the mean $I_0$, while $z=-1$ returns $I_1$.

### What Is Set, What Is Measured, and What Is Derived

At a fixed detector operating point, the quantities have different roles:

- **Calibrated detector quantities:** $I_0$, $I_1$, the noise spectral density $S$, and the total ensemble dephasing rate $\Gamma_d$ are obtained from separate experimental calibrations.
- **Chosen controls:** the averaging window $\tau$ is selected by the observer or numerical protocol, while the Rabi frequency $\Omega_R$ is set by the qubit drive.
- **Derived quantities:** $I_c$, $\Delta I$, the averaged-noise variance $D$, the measurement rate $\Gamma_m$, the measurement time $\tau_m$, the extra dephasing rate $\gamma$, and the efficiency $\eta$ follow from the calibrated quantities.
- **Random and dynamical quantities:** the reading $\bar I$ is a random outcome, while $z(t)$ and $\rho(t)$ are the conditioned qubit state inferred from the accumulated record.

The word *calibrated* is more precise here than *fundamental*. Experimental controls such as detector bias, measurement power, and amplifier gain establish an operating point. The effective model then takes the measured values $I_0$, $I_1$, $S$, and $\Gamma_d$ at that operating point as inputs.

**Step 1: Calibrate the signal contrast.** Prepare the qubit in $|0\rangle$ and $|1\rangle$ separately and measure the corresponding mean outputs $I_0$ and $I_1$. Their midpoint and separation are then fixed:

$$
I_c=\frac{I_0+I_1}{2},
\qquad
\Delta I=I_0-I_1.
$$

Once $I_0$ and $I_1$ have been calibrated, $I_c$ and $\Delta I$ are not additional independent parameters.

**Step 2: Calibrate the noise, then choose a window.** Measure the flat part of the detector's single-sided noise spectrum to obtain $S$. The observer then chooses a finite averaging window $\tau$. This choice produces the ordinary Gaussian variance

$$ D=\mathrm{Var}(\bar\xi_\tau)=\frac{S}{2\tau}. $$

The detector property is $S$; the analysis choice is $\tau$; and $D$ is derived from both. Making $\tau$ smaller makes each individual reading noisier, but it also produces more readings per unit time. It does not change the detector's intrinsic information-acquisition rate.

For the white-noise model, $\tau$ must be longer than the detector correlation time. For a discrete weak-measurement step with a driven qubit, it should also be short enough that $\Gamma_m\tau\ll1$ and $\Omega_R\tau\ll1$.

**Step 3: Derive the measurement rate.** Over a window $\tau$, the two possible readout distributions are Gaussians with the same variance $D$ and with means separated by $\Delta I$:

$$
P_j(\bar I\mid\tau)
=\frac{1}{\sqrt{2\pi D}}
\exp\!\left[-\frac{(\bar I-I_j)^2}{2D}\right],
\qquad j\in\{0,1\}.
$$

Their overlap is

$$
\mathcal O(\tau)
=\int_{-\infty}^{\infty}\sqrt{P_0(\bar I)P_1(\bar I)}\,d\bar I
=\exp\!\left[-\frac{(\Delta I)^2}{8D}\right]
=\exp\!\left[-\frac{(\Delta I)^2}{4S}\tau\right].
$$

Define the measurement rate by $\mathcal O(\tau)=e^{-\Gamma_m\tau}$. It follows that

$$ \Gamma_m=\frac{(\Delta I)^2}{4S}. $$

This formula says that information arrives faster when the two detector outputs are farther apart and slower when the detector is noisier. The units also work: $S$ has units of current squared times time, so $(\Delta I)^2/S$ has units of inverse time. Only two of the three quantities $\Delta I$, $S$, and $\Gamma_m$ are independent.

The conventional measurement time is

$$
\tau_m=\frac{2S}{(\Delta I)^2},
\qquad
\Gamma_m=\frac{1}{2\tau_m}.
$$

The factor of two is a convention worth keeping visible. At $\tau=\tau_m$, the two Gaussian means are separated by twice their common standard deviation. The likelihood overlap decays on the related time scale $1/\Gamma_m=2\tau_m$.

**Step 4: Calibrate total ensemble dephasing.** Repeat the experiment many times and ignore the individual detector records. The ensemble-averaged coherence decays as

$$ \left|\langle\rho_{01}(t)\rangle_{\mathrm{records}}\right|\propto e^{-\Gamma_d t}. $$

Thus $\Gamma_d$ is measured independently of $\Gamma_m$. Quantum mechanics requires the total dephasing associated with this measurement model to satisfy $\Gamma_d\geq\Gamma_m$: the qubit cannot reveal information in the accessible record without losing at least the corresponding ensemble coherence. The difference and ratio define

$$
\gamma=\Gamma_d-\Gamma_m\geq0,
\qquad
\eta=\frac{\Gamma_m}{\Gamma_d}\in(0,1].
$$

Here $\Gamma_m$ is the rate of information that reaches the observer, $\gamma$ is dephasing whose information is not present in the observed record, and $\eta$ is the fraction of the total dephasing that is accounted for by acquired information. Equivalently,

$$
\gamma=\Gamma_d(1-\eta)
=\Gamma_m\left(\frac{1}{\eta}-1\right).
$$

For a quantum-limited detector, $\eta=1$, so $\Gamma_d=\Gamma_m$ and $\gamma=0$. A non-ideal detector has $\eta<1$ and therefore destroys ensemble coherence faster than its observed record acquires information.

**Step 5: Keep detector symmetry separate from qubit symmetry.** For the symmetric detector used first in this essay, the output noise and phase back-action are uncorrelated. There is therefore no record-correlated unitary phase kick. In the ideal case, the conditioned back-action is purely informational and moves the state along Bloch-sphere meridians; if $\eta<1$, the unobserved part still contributes the extra dephasing $\gamma$. This detector assumption is distinct from the symmetric-qubit assumption of zero energy bias. With that qubit choice, the controlled Hamiltonian is

$$ H_{qb}=\frac{\Omega_R}{2}\sigma_x, $$

where $\Omega_R$ is an externally set Rabi frequency, not a quantity derived from $S$ or $\Gamma_m$.

The physical calibration chain is therefore

$$
\{I_0,I_1,S\}
\longrightarrow
\{I_c,\Delta I,\Gamma_m,\tau_m\},
\qquad
\{S,\tau\}\longrightarrow D,
\qquad
\{\Gamma_m,\Gamma_d\}\longrightarrow\{\gamma,\eta\}.
$$

For the numerical examples, it is convenient to reverse one part of this chain and use $\Gamma_m$ as the time unit. We normalize the detector outputs to $I_0=+1$ and $I_1=-1$, so $\Delta I=2$ and $I_c=0$, and then infer the equivalent noise density from $S=(\Delta I)^2/(4\Gamma_m)$:

```wl
ClearAll[dI, i0, i1, ic];
i0 = 1; i1 = -1; dI = i0 - i1; ic = (i0 + i1)/2;
{i0, i1, dI}
```

As one can see, the two conditioned means sit symmetrically about zero, so a reading's sign is the only hint of which state produced it, and $\Delta I=2$ is the full separation the noise must beat.

The chance that a state $|j\rangle$ produced a given reading is a Gaussian likelihood; define it as a function of the window variance $D$:

```wl
ClearAll[gaussLike];
gaussLike[mu_, var_, ibar_] := Exp[-(ibar - mu)^2/(2 var)]/Sqrt[2 Pi var]
```

This is the single probabilistic object the whole formalism turns on: everything downstream is Bayes' rule applied to these two likelihoods.

Define the window variance $D=S/(2\tau)$ and the honest detector draw, sampling $\bar I$ from the Gaussian mixture $\rho_{00}\,\mathcal N(I_0,D)+\rho_{11}\,\mathcal N(I_1,D)$ that the current state predicts:

```wl
ClearAll[avgVar, drawCurrent];
avgVar[gm_, tau_] := (dI^2/(4 gm))/(2 tau);
drawCurrent[rho_, var_] :=
  RandomChoice[Clip[Re@{rho[[1, 1]], rho[[2, 2]]}, {0., 1.}] -> {i0, i1}] +
    RandomVariate[NormalDistribution[0, Sqrt[var]]]
```

In words, the detector first picks which bell to ring with the current populations as weights, then buries the level under Gaussian noise of variance $D$. For the fixed detector used in these examples, increasing $\tau$ sharpens an individual averaged readout but does not alter $\Gamma_m$.

Now the heart of the formalism. The exact update over a step, in its cleanest form, keeps the populations Bayesian and lets each coherence carry the geometric mean of the two likelihoods, with an optional pure-dephasing factor for a non-ideal detector:

$$ \rho_{ij}\ \longrightarrow\ \rho_{ij}\,\frac{\sqrt{P_i\,P_j}}{\sum_k \rho_{kk}\,P_k}\ \times\ e^{-\gamma\tau}\ \text{(off-diagonal only)}, \qquad P_j=\text{gaussLike}(I_j,D,\bar I) . $$

Encode that update as one element-wise operation on the matrix:

```wl
ClearAll[bayesUpdate];
bayesUpdate[rho_, ibar_, var_, gt_] :=
  With[{p0 = gaussLike[i0, var, ibar], p1 = gaussLike[i1, var, ibar]},
   rho {{p0, Sqrt[p0 p1] Exp[-gt]}, {Sqrt[p0 p1] Exp[-gt], p1}}/
    (Re@rho[[1, 1]] p0 + Re@rho[[2, 2]] p1)]
```

The diagonal entries pick up $P_j$ and renormalize, which is textbook Gaussian Bayes on the populations; the off-diagonal entries pick up $\sqrt{P_0P_1}$, so the coherence is pinned to the populations rather than measured directly.

Confirm that this update is trace-preserving for any reading, symbolically, so that $\rho$ stays a density matrix:

```wl
With[{r = {{r00, r01}, {Conjugate[r01], 1 - r00}}},
 FullSimplify[Tr[bayesUpdate[r, ibar, var, gt]] == 1,
  Assumptions -> {0 < r00 < 1, var > 0, ibar \[Element] Reals, gt >= 0}]]
```

As expected, the two updated populations always sum to one: the normalizing denominator is exactly the predictive probability of the reading, so probability is conserved by construction.

A single step of the full dynamics interleaves a Hamiltonian half-rotation, the Bayesian kick, and another half-rotation, a symmetric (Strang) split that is second-order accurate in the step; for a frozen qubit the half-rotation is just the identity. Define the step, together with the readouts we will plot:

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
recordCurrents[gm_, tau_, gt_, om_, nst_, seed_] := BlockRandom[SeedRandom[seed];
  With[{uH = MatrixExp[-I (om/2) PauliMatrix[1] (tau/2)], v = avgVar[gm, tau]},
   FoldPairList[
    Function[{rho, k}, With[{r1 = uH . rho . ConjugateTranspose[uH]},
      With[{ib = drawCurrent[r1, v]},
       {ib, uH . bayesUpdate[r1, ib, v, gt] . ConjugateTranspose[uH]}]]],
    {{0.5, 0.5}, {0.5, 0.5}}, Range[nst]]]]
```

We will return to this record in Section 6; for now it is enough that each measurement produces a number, and the state and the record share the very same noise.

## 1. The Detector's Two Bells: Why One Look Cannot Separate Them

A weak measurement is weak because a single short look is ambiguous. For a state that is definitely $|0\rangle$ or definitely $|1\rangle$ the averaged current is a Gaussian centered at $I_0$ or $I_1$ with spread $\sqrt D=\sqrt{S/(2\tau)}$; the two bells overlap heavily until the window $\tau$ is long enough to pull them apart, and the crossover scale is the measurement time $\tau_m=2S/(\Delta I)^2$.

Plot the two likelihoods $P_0$ and $P_1$ against the reading $\bar I$ for windows short, comparable to, and long relative to $\tau_m$:

```wl
With[{gm = 0.5},
 With[{taum = 2 (dI^2/(4 gm))/dI^2},
  Row[Table[
    Plot[{gaussLike[i0, avgVar[gm, tau], ibar], gaussLike[i1, avgVar[gm, tau], ibar]},
     {ibar, -4, 4}, PlotRange -> All, Filling -> Axis, Frame -> True,
     GridLines -> Automatic, FrameLabel -> {"averaged current", "likelihood"},
     PlotLabel -> "window \[Tau] = " <> ToString[tau/taum] <> " \[Tau]m",
     ImageSize -> 240], {tau, taum {1/4, 1, 6}}]]]]
```

Notice that at a quarter of the measurement time the two bells sit almost on top of each other, so a single reading barely favors either state; only as $\tau$ grows past $\tau_m$ do the peaks separate and a reading starts to mean something. Information is nothing but this shrinking overlap.

## 2. Populations Follow Bayes: Localization Without a Hamiltonian

Switch off the Hamiltonian and the record can only report which state we are in, never rotate us between them. Start from the equatorial superposition $|{+}\rangle=(|0\rangle+|1\rangle)/\sqrt{2}$, whose two measured-basis populations are equal, and feed each reading through the update: the populations move by the likelihood ratio and drift toward one pole, a noisy random walk in $z$ rather than a jump. The pole that wins is random, with probability equal to the initial population, which is the Born rule emerging from accumulated evidence.

Run one frozen-qubit trajectory from $\rho=|{+}\rangle\langle{+}|$ and watch $z(t)$:

```wl
With[{gm = 1., tau = 0.02, nst = 1000},
 With[{traj = BlockRandom[SeedRandom[17];
      NestList[measureStep[IdentityMatrix[2], avgVar[gm, tau], 0.],
       {{0.5, 0.5}, {0.5, 0.5}}, nst]]},
  ListLinePlot[Transpose[{Range[0, nst] tau, zComp /@ traj}],
   PlotRange -> {-1.05, 1.05}, Frame -> True, GridLines -> Automatic,
   FrameLabel -> {"time (units of 1/\[CapitalGamma]m)", "z = \[Rho]00 - \[Rho]11"},
   PlotLabel -> "one frozen-qubit record: localization by Bayes"]]]
```

The population difference wanders under the noisy evidence and then commits to a pole and stays there; the "collapse" is the accumulation of log-likelihood, not a discontinuous jump.

## 3. The Purity Ride-Along: A Pure State Stays Pure in a Single Run

Here is the fact that makes an ideal detector special. The coherence is not measured; it rides on the populations through the geometric mean $\sqrt{P_0P_1}$, so the ratio $|\rho_{01}|^2/(\rho_{00}\rho_{11})$ is left untouched by every update. A state that starts pure therefore stays exactly pure along the whole record, never leaving the surface of the Bloch sphere, even as its populations localize.

First confirm the invariance symbolically, and get the lossy case for free: check that one update, for any reading, multiplies the ride-along ratio $|\rho_{01}|^2/(\rho_{00}\rho_{11})$ by exactly $e^{-2\gamma\tau}$, so an ideal update ($\gamma=0$) leaves it untouched:

```wl
Block[{r00, r01, ibar, var, gt},
 With[{r = {{r00, r01}, {Conjugate[r01], 1 - r00}}},
  With[{u = bayesUpdate[r, ibar, var, gt]},
   FullSimplify[
    (u[[1, 2]] u[[2, 1]])/(u[[1, 1]] u[[2, 2]]) == Exp[-2 gt] (r01 Conjugate[r01])/(r00 (1 - r00)),
    Assumptions -> {0 < r00 < 1, var > 0, ibar \[Element] Reals, gt >= 0}]]]]
```

The ratio is multiplied by the same factor whatever the record says, so for an ideal detector it is algebraically conserved and purity is a constant of the single-run motion, while a lossy detector bleeds it at the fixed rate $2\gamma$; that one factor is the mechanism behind Section 5's purity ceiling, and for an undriven record it will give the purity in closed form. Watch it hold along a record that starts pure on the equator, $\rho=|{+}\rangle\langle{+}|$, while $z$ localizes:

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

Now switch on the drive $H_{qb}=\tfrac{\Omega_R}{2}\sigma_x$ and interleave it with the measurement by the symmetric split: a half rotation, the Bayesian kick, another half rotation. The measurement pulls the state toward the $z$-poles while the field rotates it, so a single trajectory is a noisy Rabi oscillation whose phase slowly diffuses; a run that starts fully mixed at the center of the Bloch ball sharpens into that oscillation as the record accrues.

Run one noisy-Rabi trajectory from $|0\rangle$ and, beside it, a run that starts fully mixed and purifies:

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
      PlotLabel -> "pure start from |0>: noisy Rabi", ImageSize -> 250],
     ListLinePlot[Transpose[{tt, #}] & /@ Transpose[mixed],
      PlotRange -> {-1.05, 1.05}, Frame -> True, GridLines -> Automatic,
      PlotLegends -> {"x", "y", "z"}, FrameLabel -> {"time (units of 1/\[CapitalGamma]m)", "Bloch vector"},
      PlotLabel -> "mixed start: purifies into oscillation", ImageSize -> 250]}]]]]
```

The pure run oscillates and wobbles but keeps its full Bloch length, while the mixed run starts at the origin and swells outward into the same phase-locked, slowly diffusing oscillation: the record both drives and sharpens the state at once.

## 5. A Non-Ideal Detector: Ideality Caps the Purity

A realistic detector loses some of the information it gathers before it reaches the output, and that lost information appears as an extra pure-dephasing factor $e^{-\gamma\tau}$ on the coherence, with $\gamma=\Gamma_d(1-\eta)$. The measurement keeps localizing at the same informational rate, but the coherence is now bled faster than the populations sharpen, so the state settles short of the surface at a steady purity set entirely by the ideality $\eta$.

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

As expected, only the ideal detector holds the state pure; below $\eta=1$ the purity settles onto a ceiling that drops as $\eta$ drops. The mechanism is the per-step law Section 3 proved: every update multiplies the ride-along ratio by exactly $e^{-2\gamma\tau}$, so the ceiling sits where that steady bleed balances the localization that sharpens the populations. The extra dephasing $\gamma$ is potential information dissipated inside the detector, and only $\eta=1$ lets you monitor the wavefunction itself.

Switch the drive off and the same law hands over the whole answer. Along an undriven record from the equator the ratio starts at one, so at every instant the purity is pinned to the populations, $\mathrm{Tr}\,\rho^2=1-2\rho_{00}\rho_{11}\,(1-e^{-2\gamma t})$; confirm the identity along one lossy record:

```wl
With[{gm = 1., tau = 0.02, nst = 1000, eta = 0.4},
 With[{gt = gm (1/eta - 1) tau},
  With[{traj = BlockRandom[SeedRandom[8];
      NestList[measureStep[IdentityMatrix[2], avgVar[gm, tau], gt], {{0.5, 0.5}, {0.5, 0.5}}, nst]]},
   Max[Abs[MapIndexed[
      purity[#1] - (1 - 2 Re[#1[[1, 1]]] (1 - Re[#1[[1, 1]]]) (1 - Exp[-2 gt (First[#2] - 1)])) &,
      traj]]]]]]
```

The identity holds to machine precision along the whole record: with no drive, all the stochasticity lives in the populations, and the only thing the lossy detector adds is the deterministic $e^{-2\gamma t}$ bleed. The drive is what breaks this closed form, trading population and coherence between kicks, and that is why the driven ceiling above needed an ensemble.

## 6. The Output Spectrum and the Factor of Four

Now the payoff. When the drive keeps the qubit oscillating, the Rabi motion modulates the very current the detector reports, so the record carries a coherent line at $\Omega_R$ riding on the white noise floor $S$. Korotkov and Averin derived the line in closed form, and its height above the floor is bounded: at most four times the noise in the good-oscillation regime $\Omega_R\gg\Gamma$, reached only for an ideal detector, and degrading as $4\eta$. This spectrum is not just a calculation: Palacios-Laloy and co-workers measured exactly this Rabi line in a continuously monitored superconducting qubit, and read its bounded height as a Leggett-Garg (temporal Bell) test of macroscopic realism (Nature Physics 6, 442 (2010), arXiv:1005.3435).

$$ S_I(\omega)=S+4\eta\,S\,\frac{\Omega_R^2\,\Gamma^2}{(\omega^2-\Omega_R^2)^2+\Gamma^2\omega^2}, \qquad \Gamma=\Gamma_d=\frac{\Gamma_m}{\eta}. $$

Half of that four is the ordinary spectrum of $z(t)$; the other half is a purely nonclassical correlation between the detector noise $\xi$ and the back-action it exerts on $z$, which no classical harmonic signal can produce. The record we simulate carries both halves automatically, because the same reading that enters the spectrum also drives the Bayesian kick. My suggestion for the essay's longest-running cell: read the output first, the overlaid closed form second, and the input code last.

Before simulating, read the bound off the closed form itself. Solve it for its one stationary point:

```wl
Block[{eta, om, gd, w},
 With[{s = 1 + 4 eta om^2 gd^2/((w^2 - om^2)^2 + gd^2 w^2)},
  FullSimplify[Solve[D[s, w] == 0 && w > 0, w, Reals], {om > 0, gd > 0}]]]
```

The one stationary point, $\omega^2=\Omega_R^2-\Gamma^2/2$, exists only while $\Gamma<\sqrt2\,\Omega_R$, a weak-damping condition the good-oscillation regime $\Omega_R\gg\Gamma$ satisfies with room to spare. Confirm that point is a maximum:

```wl
Block[{eta, om, gd, w},
 With[{s = 1 + 4 eta om^2 gd^2/((w^2 - om^2)^2 + gd^2 w^2)},
  FullSimplify[(D[s, {w, 2}] /. w -> Sqrt[om^2 - gd^2/2]) < 0,
   Assumptions -> {om > 0, eta > 0, 0 < gd < Sqrt[2] om}]]]
```

It is, so the line has a single peak; take its height in the narrow-line limit $\Gamma\to 0$:

```wl
Block[{eta, om, gd, w},
 With[{s = 1 + 4 eta om^2 gd^2/((w^2 - om^2)^2 + gd^2 w^2)},
  Limit[FullSimplify[s /. w -> Sqrt[om^2 - gd^2/2]], gd -> 0, Assumptions -> {om > 0, eta > 0}]]]
```

The height settles at $1+4\eta$ only as the linewidth closes; at any finite $\Gamma/\Omega_R$ the closed form's peak sits marginally above it, so "at most four above the floor" is a statement about the good-oscillation limit, which is where we now simulate.

Simulate long records in the good-oscillation regime $\Omega_R\gg\Gamma_m$, average the periodogram over an ensemble for three idealities, and overlay the closed form (dashed):

```wl
With[{gm = 0.05, om = 1., tau = 0.05, nst = 14000, ntr = 200, etas = {1., 0.5, 0.25}},
 Block[{gd, sClosed, spec, specs},
  gd[eta_] := gm/eta;
  sClosed[eta_, w_] := 1 + 4 eta om^2 gd[eta]^2/((w^2 - om^2)^2 + gd[eta]^2 w^2);
  spec[recs_] := With[{half = Quotient[Length[recs[[1]]], 2]},
    With[{ws = Table[2. Pi (s - 1)/(Length[recs[[1]]] tau), {s, 2, half}],
      psd = Mean[(Abs[Fourier[#]]^2)[[2 ;; half]] & /@ recs]},
     Transpose[{ws[[2 ;; -2]],
       MovingAverage[psd/Median[Pick[psd, UnitStep[ws - 2.5], 1]], 3]}]]];
  If[Kernels[] === {}, LaunchKernels[]]; DistributeDefinitions["Global`"];
  specs = Table[
    With[{gt = (gd[eta] - gm) tau},
     spec[ParallelTable[recordCurrents[gm, tau, gt, om, nst, 7000 + Round[1000 eta] + k],
       {k, ntr}]]], {eta, etas}];
  Show[
   ListLinePlot[specs, PlotRange -> {{0, 2.2}, {0, 6}}, Frame -> True, GridLines -> Automatic,
    FrameLabel -> {"angular frequency \[Omega]", "spectral density / noise floor"},
    PlotLegends -> {"simulated \[Eta] = 1", "simulated \[Eta] = 0.5", "simulated \[Eta] = 0.25"},
    PlotLabel -> "output spectrum: the coherent line and the factor of four"],
   ListLinePlot[Table[Table[{w, sClosed[eta, w]}, {w, 0.2, 2.2, 0.002}], {eta, etas}],
    PlotRange -> All, PlotStyle -> Directive[Gray, Thick, Dashed]], ImageSize -> 540]]]
```

Each simulated line sits on the closed-form Lorentzian, and the coherent peak climbs toward four times the floor for the ideal detector and falls off as $4\eta$: the line is bounded because half of its area is the noise-backaction correlation, which saturates only when every bit of dephasing is informational. The dashed peak of the sharpest line sits slightly above its simulated partner; that gap is the finite frequency bins and the smoothing window shaving the top off a narrow peak, not physics.

## 7. Cross-Check: The Stratonovich SDE and the Ito Trap

The same physics is a Bloch-Langevin stochastic differential equation in the sibling essay `NDSolve vs Ito/`. The Bayesian update and that Stratonovich equation must agree as $\tau\to 0$, and comparing them exposes the classic Ito trap. The ensemble coherence decays at the total rate $\Gamma_d$; in the Ito form of the coherence equation that decay is the term $+(\Delta I)^2/4S$, which is not physical decoherence at all but the Stratonovich-to-Ito drift correction. In other words, the term that looks like a decoherence rate in the Ito equation is a bookkeeping consequence of the calculus, and dropping it does not remove a small correction, it removes the whole dephasing. Hand the Stratonovich right-hand side to a plain Ito-Euler stepper without that term and the coherence decays at the wrong rate.

The exact ensemble decay is one Gaussian integral away from `gaussLike`, so derive it before simulating it. Average one update over the readings the state predicts: each population picks up its own likelihood, which integrates to one, while the coherence picks up the geometric mean of the two, whose integral is the Gaussian overlap; evaluate both, the second at the detector's levels and window variance:

```wl
Block[{a0, a1, var, ibar, gm, tau},
 {Integrate[gaussLike[a0, var, ibar], {ibar, -Infinity, Infinity}, Assumptions -> var > 0],
  FullSimplify[
   Integrate[Sqrt[gaussLike[a0, var, ibar] gaussLike[a1, var, ibar]], {ibar, -Infinity, Infinity},
     Assumptions -> {var > 0, {a0, a1} \[Element] Reals}] /.
    {a0 -> i0, a1 -> i1, var -> avgVar[gm, tau]}, Assumptions -> {gm > 0, tau > 0}]}]
```

The populations are a martingale, untouched on average by measurement alone, and every step multiplies the ensemble-averaged coherence by exactly $e^{-\Gamma_m\tau}$ with no small-step expansion: the dashed reference in the plot below is derived, not fitted, and Section 3's claim that the ensemble dephases only because runs localize differently is now an identity about the overlap of the two bells.

Overlay the ensemble coherence $\langle x\rangle(t)$ from three routes started on the equator: the Bayesian update, the correct Ito stepper with the drift term, and the naive stepper missing it, against the exact decay $e^{-\Gamma_m t}$ (dashed):

```wl
With[{gm = 1., tau = 0.005, nst = 300, ntr = 1500},
 Block[{v = avgVar[gm, tau], bayes, ito, tt = Range[0, nst] tau, run},
  bayes[rho_] := measureStep[IdentityMatrix[2], v, 0.][rho];
  ito[drift_][{x_, z_}] := With[{xb = RandomVariate[NormalDistribution[0, Sqrt[v]]]},
    {x - z (2 gm) x xb tau - drift gm x tau, z + (1 - z^2) (2 gm) xb tau}];
  run[stepper_, init_, read_, seed_] := Mean[Table[read /@
      BlockRandom[SeedRandom[seed + k]; NestList[stepper, init, nst]], {k, ntr}]];
  ListLinePlot[{
    Transpose[{tt, run[bayes, {{0.5, 0.5}, {0.5, 0.5}}, (2 Re@#[[1, 2]] &), 1]}],
    Transpose[{tt, run[ito[1], {1., 0.}, First, 2]}],
    Transpose[{tt, run[ito[0], {1., 0.}, First, 3]}]},
   PlotRange -> {0, 1.06}, Frame -> True, GridLines -> Automatic,
   FrameLabel -> {"time (units of 1/\[CapitalGamma]m)", "ensemble coherence \[LeftAngleBracket]x\[RightAngleBracket]"},
   PlotLegends -> {"Bayesian update", "Ito with +(\[CapitalDelta]I)^2/4S drift", "naive Ito, drift dropped"},
   PlotLabel -> "the Ito trap: the drift correction IS the dephasing",
   Epilog -> {Black, Dashed, Line[Transpose[{tt, Exp[-gm tt]}]]}]]]
```

This confirms that the Bayesian update and the correct Ito stepper both fall onto the exact $e^{-\Gamma_m t}$ curve, while the naive stepper never dephases at all: the missing $+(\Delta I)^2/4S$ term was the entire ensemble decoherence. The Bayesian form has no calculus ambiguity to begin with, because it is discrete Gaussian Bayes between unitary half-steps; the stochastic differential equation reproduces it only when the Stratonovich-to-Ito drift is kept.

## 8. The Other Back-Action: The Phase Kick and the Second Quadrature

Every detector in this essay has been symmetric: a quantum point contact whose two output levels carry the same noise, so the correlation between the detector's output noise and the disturbance it feeds back to the qubit vanishes, $K=0$. That is exactly why the back-action has been purely informational. Each reading squeezes the populations and slides the state along a meridian of the Bloch sphere, changing $z$, and nothing turns the coherence. A general detector is not so gentle: alongside the informational disturbance it applies a second, physical one, an active rotation of the coherence about the $z$-axis, driven by the very same record.

In other words, there are two back-actions, and so far we have felt only one. Korotkov names them "spooky" and "realistic". The spooky (informational) back-action is the non-unitary Bayes factor we have used throughout: it carries no physical mechanism, it is pure information, and it moves the state along meridians, changing $z$. The realistic (phase) back-action is unitary: it is a genuine rotation of the coherence about $z$, along the parallels of constant latitude, and it needs a physical mechanism, such as the fluctuating level shift an asymmetric detector imprints on the qubit. For a symmetric detector the second one is switched off; for a general one it is not.

Concretely, the coherence picks up one extra factor. Where the ideal update multiplied $\rho_{01}$ by the geometric mean of the two likelihoods, a general detector multiplies it also by a record-driven phase,

$$ \rho_{01}\ \longrightarrow\ \rho_{01}\,\frac{\sqrt{P_0 P_1}}{\sum_k \rho_{kk} P_k}\,e^{-\,iK\bar I\tau}, $$

a rotation of the coherence through an angle $K\bar I\tau$ set by the correlation $K$ and the integrated current $\bar I\tau$ of that step. In circuit QED this second channel is not abstract: a phase-preserving amplifier reports two output quadratures at once, and while the in-phase quadrature $I(t)$ carries the informational back-action, the second quadrature $Q(t)$ is the one that carries the phase kick. Murch, Weber, Macklin, and Siddiqi made the choice between the two explicit (Nature 502, 211 (2013), arXiv:1305.7270): running a near-quantum-limited amplifier phase-sensitively, they amplified one quadrature and de-amplified the other, so the qubit felt only the back-action they selected. Reading the informational quadrature, single-shot trajectories rode a meridian of the Bloch sphere toward the poles; reading the phase quadrature, they rode the equator as the coherence turned, both tracked at a quantum efficiency of 0.49.

Because the phase kick is unitary, it never touches a single run's purity: the ride-along of Section 3 still holds exactly, and a pure state stays on the surface of the Bloch sphere. It shows up instead in the ensemble, because different records wind the coherence by different angles and the average of $e^{-iK\bar I\tau}$ shrinks toward zero. The total ensemble dephasing rate therefore splits into three pieces,

$$ \Gamma=\underbrace{\frac{(\Delta I)^2}{4S}}_{\text{informational}}+\underbrace{\frac{K^2 S}{4}}_{\text{phase}}+\underbrace{\gamma}_{\text{extra environment}}, $$

and the fraction of that rate which is informational is precisely the number that fixes how far the Section 6 line rises above its noise floor: the peak settles at $1+4\tilde\eta$ times the floor, with $\tilde\eta=(\Delta I)^2/(4S\Gamma)$. In Section 6 the rise degraded as $4\eta$ only because $K=0$ made $\tilde\eta=\eta$; in general the rise is $4\tilde\eta$, and the phase back-action pushes it down even when the detector wastes no information internally ($\eta=1$, $\gamma=0$). A phase-preserving amplifier splits its dephasing evenly between the two quadratures, so $\tilde\eta\le\tfrac12$ and its coherent line can rise at most twice above the noise floor, never four times.

The phase rate is one Gaussian integral away, exactly as the informational rate was. Average the phase factor $e^{-iK\bar I\tau}$ over the readings a state at level $\mu$ produces:

```wl
Block[{mu, kk, tau, ss, ibar},
 FullSimplify[Integrate[gaussLike[mu, ss/(2 tau), ibar] Exp[-I kk ibar tau],
   {ibar, -Infinity, Infinity}, Assumptions -> {ss > 0, tau > 0, {mu, kk} \[Element] Reals}],
  {ss > 0, tau > 0, {mu, kk} \[Element] Reals}]]
```

The argument advances with the mean level $\mu$; read off the modulus, the ensemble dephasing the phase kick leaves behind:

```wl
Block[{mu, kk, tau, ss, ibar},
 FullSimplify[
  Abs[Integrate[gaussLike[mu, ss/(2 tau), ibar] Exp[-I kk ibar tau],
     {ibar, -Infinity, Infinity}, Assumptions -> {ss > 0, tau > 0, {mu, kk} \[Element] Reals}]] ==
   Exp[-kk^2 ss tau/4], Assumptions -> {ss > 0, tau > 0, {mu, kk} \[Element] Reals}]]
```

The modulus decays at exactly $K^2S/4$ per unit time: the second piece of the dephasing split is derived from the same `gaussLike` as the first, and it is pure ensemble loss, invisible in any single run. The general weak-measurement (Kraus) operator for a reading $\bar I$ factorizes into a unitary times a positive part,

$$ M_{\bar I}=U_{\bar I}\,\sqrt{M_{\bar I}^\dagger M_{\bar I}}, \qquad U_{\bar I}=e^{-iK\bar I\tau\,\sigma_z/2}, \qquad \sqrt{M_{\bar I}^\dagger M_{\bar I}}=\sqrt{P_0}\,|0\rangle\langle 0|+\sqrt{P_1}\,|1\rangle\langle 1| . $$

The positive square-root factor is exactly the informational Bayes update we built in Section 2; the unitary $U_{\bar I}$ is the phase rotation, which we quietly set to the identity by taking $K=0$. Everything in this essay so far has been the special case $U_{\bar I}=\mathbb 1$.

Turn the phase kick on and watch it act. Run one record through both detectors from the same pure equatorial start $\rho=|{+}\rangle\langle{+}|$, and plot the coherence in the equatorial $(x,y)$ plane for $K\neq 0$ beside the symmetric $K=0$ run:

```wl
With[{gm = 0.3, tau = 0.05, nst = 80, kk = 1., seed = 3},
 With[{v = avgVar[gm, tau], plus = {{0.5, 0.5}, {0.5, 0.5}}},
  With[{phaseStep = Function[{kval}, Function[rho,
      With[{ib = drawCurrent[rho, v]},
       With[{uk = DiagonalMatrix[{Exp[-I kval ib tau/2], Exp[I kval ib tau/2]}]},
        uk . bayesUpdate[rho, ib, v, 0.] . ConjugateTranspose[uk]]]]]},
   With[{pathK = (Most @* blochVec) /@
        BlockRandom[SeedRandom[seed]; NestList[phaseStep[kk], plus, nst]],
      path0 = (Most @* blochVec) /@
        BlockRandom[SeedRandom[seed]; NestList[phaseStep[0.], plus, nst]]},
    ListLinePlot[{pathK, path0}, AspectRatio -> 1,
     PlotRange -> {{-1.05, 1.05}, {-1.05, 1.05}}, Frame -> True,
     GridLines -> Automatic, FrameLabel -> {"Bloch x", "Bloch y"},
     PlotLegends -> {"K \[NotEqual] 0: phase + informational", "K = 0: informational only"},
     PlotLabel -> "same record: the phase back-action winds the coherence about z",
     Epilog -> {GrayLevel[0.6], Dashed, Circle[{0, 0}, 1]}]]]]]
```

The same record that localizes the populations now also carries the coherence around the $z$-axis: the symmetric detector slides the state straight inward along a meridian with the phase pinned (the radial $K=0$ track), while the general detector adds the parallel motion, a unitary turn through an angle equal to $K$ times the integrated current (the winding $K\neq 0$ track), so every single run stays exactly as pure as before even though the ensemble of records now dephases faster. Turning $U_{\bar I}$ fully on, and letting two such detectors act at once, is where the essay goes next.

Everything up to here has been one detector on one axis. The Kraus operator of Section 8 is the hinge for the two generalizations that fill the rest of the essay: keep the single qubit but point a second detector at a different axis, so the two informational kicks no longer commute (Sections 9 to 15); then keep the single detector but let it read a joint operator on two qubits, so a measurement alone can entangle them (Sections 16 to 20). Both are the same Gaussian measurement operator, applied twice or applied to a larger space, and each one has now been carried out in a real superconducting-qubit laboratory.

## 9. A Second Detector: A Kraus Operator for Any Axis

Section 8 wrote one step of the measurement as $M_{\bar I}=U_{\bar I}\sqrt{M_{\bar I}^\dagger M_{\bar I}}$, the phase kick times the informational Bayes factor. Keep the phase kick off ($U_{\bar I}=\mathbb 1$, a phase-sensitive amplifier on the optimal quadrature) but stop requiring the measured axis to be $\sigma_z$: a linear detector can read any Pauli $A$, with the same Gaussian informational factor, now diagonal in the eigenbasis of $A$. Point two such detectors at one qubit at once, the first on $\sigma_z$ and the second on $\sigma_\varphi=\sigma_z\cos\varphi+\sigma_x\sin\varphi$, an axis in the $xz$-plane at angle $\varphi$ from $z$. They emit two noisy records,

$$ I_z(t)=\mathrm{Tr}[\sigma_z\rho] + \xi_z(t), \qquad I_\varphi(t)=\mathrm{Tr}[\sigma_\varphi\rho] + \xi_\varphi(t), $$

with independent white noises of single-sided density $S=1/\Gamma$; both channels are equally strong, with rate $\Gamma$, and ideal, $\eta=1$, so a step-averaged reading has window variance $D=1/(2\Gamma\Delta t)$ as in Part I.

The Bloch vector and purity are already defined in Part I; add the angled axis $\sigma_\varphi$ and the Bloch length $|\mathbf r|$:

```wl
ClearAll[sigmaPhi, rlen, blochRho];
sigmaPhi[phi_] := Cos[phi] PauliMatrix[3] + Sin[phi] PauliMatrix[1];
rlen[rho_] := Norm[blochVec[rho]];
blochRho[x_, y_, z_] := (IdentityMatrix[2] + x PauliMatrix[1] + y PauliMatrix[2] + z PauliMatrix[3])/2
```

The Bloch length is one for a pure state and shrinks for a mixed one, the single number that says whether the two records monitor a wavefunction or only its blurred image; `blochRho` writes the general qubit state from its Bloch coordinates, the state every symbolic check below runs on. For a qubit the Bloch length carries the same content as the purity, $\mathrm{Tr}\,\rho^2=(1+|\mathbf r|^2)/2$; confirm the two agree for any qubit:

```wl
Block[{x, y, z},
 FullSimplify[purity[blochRho[x, y, z]] == (1 + rlen[blochRho[x, y, z]]^2)/2, {x, y, z} \[Element] Reals]]
```

The identity holds, so a steady Bloch length and a steady purity are the same statement, which we use in Section 13. Now the operator the rest of Part II turns on: the informational factor for a reading $\bar I$ of a Pauli $A$, the matrix function that multiplies each eigen-amplitude by the square root of its likelihood:

```wl
ClearAll[kraus];
kraus[ibar_, A_, dvar_] := MatrixExp[(ibar/(2 dvar)) A]
```

Because $A^2=\mathbb 1$, this is $\cosh(\bar I/2D)\,\mathbb 1+\sinh(\bar I/2D)\,A$, applied as $\rho\to M\rho M/\mathrm{Tr}(M\rho M)$. For $A=\sigma_z$ it is diagonal and is exactly the Part I population update in operator form; look at it:

```wl
kraus[ibar, PauliMatrix[3], dvar] // MatrixForm
```

As one can see, the single-axis operator is diagonal, $\mathrm{diag}(e^{\bar I/2D},e^{-\bar I/2D})$: it scales the two amplitudes in the measured basis and leaves the coherence riding along. "Exactly the Part I update in operator form" is a checkable statement, so check it: conjugating any state by the $\sigma_z$ operator and normalizing must reproduce `bayesUpdate` for any reading:

```wl
Block[{r00, r01, ibar, var},
 With[{r = {{r00, r01}, {Conjugate[r01], 1 - r00}}, m = kraus[ibar, PauliMatrix[3], var]},
  FullSimplify[bayesUpdate[r, ibar, var, 0] == With[{s = m . r . ConjugateTranspose[m]}, s/Tr[s]],
   Assumptions -> {0 < r00 < 1, var > 0, ibar \[Element] Reals}]]]
```

The two forms agree identically, so everything Part I established, the Bayes localization, the ride-along, the ensemble decay, is inherited by the operator form. Everything below is what a second operator, diagonal in a different basis, does to it.

Define the honest two-outcome draw for a general axis, sampling $\pm 1$ with the Born weights $(1\pm\mathrm{Tr}[A\rho])/2$ and burying the level under window noise:

```wl
ClearAll[drawRead];
drawRead[rho_, A_, dvar_] :=
  RandomChoice[Clip[{(1 + Re@Tr[A . rho])/2, (1 - Re@Tr[A . rho])/2}, {0., 1.}] -> {1., -1.}] +
   RandomVariate[NormalDistribution[0, Sqrt[dvar]]]
```

Each detector picks which eigenvalue to ring with the Born weights of its own axis, then adds Gaussian noise of variance $D$; the two channels draw independently from the same state.

For the sampler and the operator to describe one and the same detector, two identities must hold. The first needs a measure on readings: take the Gaussian base measure $w(\bar I)=\mathcal N(0,D)\,e^{-1/(2D)}$, its constant factor chosen so that the squared kick operators resolve the identity and the continuum of readings is a legitimate POVM; confirm the resolution:

```wl
Block[{ib, dvar},
 With[{m = kraus[ib, PauliMatrix[3], dvar]},
  Integrate[gaussLike[0, dvar, ib] Exp[-1/(2 dvar)] m . m, {ib, -Infinity, Infinity},
   Assumptions -> dvar > 0]]]
```

The identity comes back, so the readings form a POVM over $w$. Second, the Born mixture `drawRead` samples must be exactly that POVM's predictive density $\mathrm{Tr}[M\rho M]\,w(\bar I)$; confirm the two densities agree for every state:

```wl
Block[{x, y, z, ib, dvar},
 With[{m = kraus[ib, PauliMatrix[3], dvar]},
  FullSimplify[Tr[m . blochRho[x, y, z] . m] gaussLike[0, dvar, ib] Exp[-1/(2 dvar)] ==
    (1 + z)/2 gaussLike[1, dvar, ib] + (1 - z)/2 gaussLike[-1, dvar, ib],
   Assumptions -> {dvar > 0, {x, y, z, ib} \[Element] Reals}]]]
```

They agree identically: the Monte-Carlo detector and the algebraic one are the same object, so what the ensembles below sample is exactly what the symbolic identities constrain.

## 10. Why the Closed Form Breaks: The Two Kicks Do Not Commute

A second detector on $\sigma_\varphi$ contributes its own Gaussian operator, diagonal in the $\sigma_\varphi$ basis rather than the $\sigma_z$ basis. To apply both over one step we multiply the two operators, and the order matters unless they commute. They do not, and the mismatch is the entire reason simultaneous non-commuting measurement is a distinct process. Compute the commutator of the two single-axis kicks for $\sigma_z$ and $\sigma_x$:

```wl
Block[{a, b, dvar},
 FullSimplify[kraus[a, PauliMatrix[3], dvar] . kraus[b, PauliMatrix[1], dvar] -
   kraus[b, PauliMatrix[1], dvar] . kraus[a, PauliMatrix[3], dvar]]]
```

The commutator is nonzero and off-diagonal, so there is no basis in which both kicks are diagonal at once, and therefore no closed-form population update. Confirm the mismatch is exactly a kick about the third axis, $2i\sinh(a/2D)\sinh(b/2D)\,\sigma_y$:

```wl
Block[{a, b, dvar},
 FullSimplify[
  kraus[a, PauliMatrix[3], dvar] . kraus[b, PauliMatrix[1], dvar] -
    kraus[b, PauliMatrix[1], dvar] . kraus[a, PauliMatrix[3], dvar] ==
   2 I Sinh[a/(2 dvar)] Sinh[b/(2 dvar)] PauliMatrix[2],
  Assumptions -> {a \[Element] Reals, b \[Element] Reals, dvar > 0}]]
```

The two kicks fail to commute by a rotation about $y$, inheriting the algebra $[\sigma_z,\sigma_x]=2i\sigma_y$ of the observables they measure. The mismatch is first order in each reading, so the two orderings agree as the step shrinks; but at any finite step the joint measurement mixes the two channels, and we must carry the full $2\times 2$ matrix. Split the step the same way Part I split `drawCurrent` from `bayesUpdate`: name the sandwiched step operator once, wrap it in a deterministic update for a given pair of readings, and let the step feed that update two honest draws:

```wl
ClearAll[stepOp, measureUpdateNC, measureStepNC];
stepOp[ibz_, ibp_, phi_, dvar_] :=
  kraus[ibz/2, PauliMatrix[3], dvar] . kraus[ibp, sigmaPhi[phi], dvar] . kraus[ibz/2, PauliMatrix[3], dvar];
measureUpdateNC[rho_, ibz_, ibp_, phi_, dvar_] :=
  With[{r1 = With[{mm = stepOp[ibz, ibp, phi, dvar]}, mm . rho . ConjugateTranspose[mm]]}, r1/Tr[r1]];
measureStepNC[phi_, dvar_][rho_] :=
  measureUpdateNC[rho, drawRead[rho, PauliMatrix[3], dvar], drawRead[rho, sigmaPhi[phi], dvar], phi, dvar]
```

Because the update is now a plain function of the readings, its physicality can be checked once and for all rather than at one seed. Confirm that it keeps any qubit state trace one, and that the step operator is Hermitian, so conjugating by it preserves Hermiticity, for any pair of readings at any angle:

```wl
Block[{x, y, z, ibz, ibp, phi, dvar},
 FullSimplify[{Tr[measureUpdateNC[blochRho[x, y, z], ibz, ibp, phi, dvar]] == 1,
   stepOp[ibz, ibp, phi, dvar] == ConjugateTranspose[stepOp[ibz, ibp, phi, dvar]]},
  Assumptions -> {x^2 + y^2 + z^2 <= 1, {x, y, z, ibz, ibp, phi} \[Element] Reals, dvar > 0}]]
```

Both hold identically: the denominator is the predictive probability of the reading pair, and a Hermitian sandwich keeps $\rho$ Hermitian, so the general update is a legitimate quantum operation even though it is not a closed-form Bayes rule. One more identity does what no plot can. Each single-axis factor is $e^{(\bar I/2D)A}$ with $A$ a traceless Pauli, so its determinant is $e^{(\bar I/2D)\,\mathrm{Tr}A}=1$; confirm the whole sandwiched step operator is unimodular:

```wl
Block[{ibz, ibp, phi, dvar},
 FullSimplify[Det[stepOp[ibz, ibp, phi, dvar]] == 1,
  Assumptions -> {{ibz, ibp, phi} \[Element] Reals, dvar > 0}]]
```

Since the determinant is multiplicative, the normalized update sends $\det\rho\to\det\rho/(\mathrm{Tr}\,M\rho M)^2$: a pure state, $\det\rho=0$, stays exactly pure for every reading, every angle, every step. The great-circle walk of the next section never leaves the sphere as a matter of algebra, not of numerical luck.

## 11. Orthogonal Axes: Collapse Becomes Diffusion

Point the second detector at right angles to the first, $\varphi=\pi/2$, so the two watched axes are $z$ and $x$, maximally incompatible. Now each detector is forever trying to collapse the state onto its own axis while the other pulls it away, and neither wins: instead of localizing to a pole, the state performs a persistent random walk. On ideal detectors a pure state never leaves the surface of the Bloch sphere (the determinant identity of Section 10), and a mixed start is driven onto it as the record accrues. The other half of the great-circle claim is that the walk never leaves the $xz$-plane, and that too is algebra, not luck: both measured axes lie in the $xz$-plane, so both Kraus factors are real, and a state with no $y$-component can never acquire one. Confirm it for any readings at any angle:

```wl
Block[{x, z, ibz, ibp, phi, dvar},
 FullSimplify[Tr[measureUpdateNC[blochRho[x, 0, z], ibz, ibp, phi, dvar] . PauliMatrix[2]] == 0,
  Assumptions -> {x^2 + z^2 <= 1, {x, z, ibz, ibp, phi} \[Element] Reals, dvar > 0}]]
```

With the sphere pinned by the determinant and the plane pinned by reality, the walk is confined to the $xz$ great circle exactly. Run one orthogonal-measurement trajectory from the fully mixed center, keep it for the next plot, and watch the Bloch coordinates:

```wl
With[{gm = 1., dt = 0.02, nst = 800},
 With[{dvar = 1./(2 gm dt), tt = Range[0, nst] dt},
  xzWalk = BlockRandom[SeedRandom[4];
    blochVec /@ NestList[measureStepNC[Pi/2, dvar], {{0.5, 0.}, {0., 0.5}}, nst]];
  ListLinePlot[Transpose[{tt, #}] & /@ Transpose[xzWalk],
   PlotRange -> {-1.05, 1.05}, Frame -> True, GridLines -> Automatic,
   PlotLegends -> {"x", "y", "z"}, FrameLabel -> {"time (units of 1/\[CapitalGamma])", "Bloch vector"},
   PlotLabel -> "orthogonal measurement: persistent diffusion, no collapse", ImageSize -> 460]]]
```

The $z$ and $x$ coordinates wander over the whole range and never settle, while $y$ sits at exactly zero, as just proved: the state diffuses on the $xz$ great circle, driven back and forth by the two competing records. There is no pole to fall into, so the collapse of the single-axis case is replaced by a walk that never ends. See the very same trajectory as a path on the sphere by plotting it in the $xz$-plane against the unit circle:

```wl
ListLinePlot[xzWalk[[All, {1, 3}]], AspectRatio -> 1, Frame -> True, GridLines -> Automatic,
 PlotRange -> {{-1.1, 1.1}, {-1.1, 1.1}}, FrameLabel -> {"x", "z"},
 PlotLabel -> "the same record as a walk on the xz great circle",
 Epilog -> {Gray, Dashed, Circle[]}, ImageSize -> 380]
```

The path hugs the unit circle, wandering around it without a preferred point: once the record has pulled the mixed start onto the surface, the two orthogonal detectors monitor the wavefunction perfectly and it stays pure, its position on the great circle left to diffuse freely. This is the isotropic persistent diffusion the single-basis picture has no room for.

## 12. Turning the Angle: From Collapse to Localized Diffusion to a Free Walk

The angle $\varphi$ between the two watched axes tunes the whole character of the motion. At $\varphi=0$ the two detectors watch the same axis, their kicks commute, and the state collapses to a $z$-pole exactly as one strong measurement would. As $\varphi$ opens up, the state can no longer reach a single point, but the axes still leave an imprint: it localizes to regions near the two eigen-directions and jitters there. At $\varphi=\pi/2$ the imprint vanishes and the walk is uniform.

Watching this takes an ensemble of trajectories, and an ensemble is where stepping one trajectory at a time gets slow. It does not have to: every operation inside `measureStepNC`, the Born weights, the draws, the $\cosh+\sinh$ form of Section 9's operator, the products, the normalization, is arithmetic on matrix entries, and arithmetic does not care whether each entry holds one trajectory's number or the vector of that number across every trajectory. So stack the ensemble into a single array, $n$ density matrices as one $n\times 2\times 2$ array, and advance them all in one step: entry $(i,j)$ of every matrix at once is the slice `rr[[All, i, j]]`, the matrix product is threaded over the stack, and the two random draws become two vector draws. One care keeps it fast: `Transpose` of a hand-assembled nested list comes back unpacked, so each helper repacks its result to keep the arithmetic on packed arrays. Encode the stacked step from those same ingredients, with the deterministic kick split out exactly as `measureUpdateNC` was:

```wl
ClearAll[ensembleDot, ensembleKraus, ensembleDraw, ensembleKick, ensembleStepNC];
ensembleDot[x_, y_] := Developer`ToPackedArray @ MapThread[Dot, {x, y}];
ensembleKraus[ib_, A_, dvar_] := With[{c = Cosh[ib/(2 dvar)], s = Sinh[ib/(2 dvar)]},
   Developer`ToPackedArray @ Transpose[
    {{c + s A[[1, 1]], s A[[1, 2]]}, {s A[[2, 1]], c + s A[[2, 2]]}}, {2, 3, 1}]];
ensembleDraw[rr_, A_, dvar_] := With[
   {pUp = (1. + Clip[Re[A[[1, 1]] rr[[All, 1, 1]] + A[[1, 2]] rr[[All, 2, 1]] +
          A[[2, 1]] rr[[All, 1, 2]] + A[[2, 2]] rr[[All, 2, 2]]], {-1., 1.}])/2},
   2. UnitStep[pUp - RandomReal[1., Length[rr]]] - 1. +
    RandomVariate[NormalDistribution[0, Sqrt[dvar]], Length[rr]]];
ensembleKick[rr_, ibz_, ibp_, phi_, dvar_] := With[{sz = N@PauliMatrix[3], sp = N@sigmaPhi[phi]},
   With[{mz = ensembleKraus[ibz/2, sz, dvar]},
    With[{r2 = With[{mm = ensembleDot[mz, ensembleDot[ensembleKraus[ibp, sp, dvar], mz]]},
        ensembleDot[ensembleDot[mm, rr], mm]]},
     r2/Re[r2[[All, 1, 1]] + r2[[All, 2, 2]]]]]];
ensembleStepNC[phi_, dvar_] := With[{sz = N@PauliMatrix[3], sp = N@sigmaPhi[phi]},
   Function[rr, ensembleKick[rr, ensembleDraw[rr, sz, dvar], ensembleDraw[rr, sp, dvar], phi, dvar]]]
```

The step is the same symmetric split as before, the half $z$-kick reused on both sides of the $\varphi$-kick with a single draw behind it, and `UnitStep` against a vector of uniform draws replaces the one-at-a-time eigenvalue choice; $M_zM_\varphi M_z$ is symmetric, so the conjugation needs no transpose. Before trusting the stacked code with any physics, pin it to the code the physics was proved about: a stack of one state against `measureUpdateNC`, over a grid of readings, angles, and variances:

```wl
With[{rho = {{0.6, 0.2 + 0.1 I}, {0.2 - 0.1 I, 0.4}}},
 Max @ Table[Max[Abs[First[ensembleKick[{rho}, {ibz}, {ibp}, phi, dvar]] -
      measureUpdateNC[rho, ibz, ibp, phi, dvar]]],
   {ibz, {-1.3, 0.7}}, {ibp, {-0.4, 2.1}}, {phi, {0.3, Pi/2}}, {dvar, {1.5, 25.}}]]
```

The two implementations agree at machine precision, so every identity Section 10 proved about the single update, trace, Hermiticity, determinant, the pinned $y$, transfers to the stacked run wholesale. Now run the stacked ensemble from the center for three angles and read, at the end of each run, how far the state sits from the equator, the mean of $|z|$, keyed by the angle:

```wl
With[{gm = 1., dt = 0.03, nst = 250, ntr = 200},
 With[{dvar = 1./(2 gm dt)},
  AssociationMap[
   Function[phi, BlockRandom[SeedRandom[1 + Round[100 phi]];
     With[{fin = Nest[ensembleStepNC[phi, dvar], ConstantArray[{{0.5, 0.}, {0., 0.5}}, ntr], nst]},
      Mean[Abs[Re[fin[[All, 1, 1]] - fin[[All, 2, 2]]]]]]]], {0., Pi/4, Pi/2}]]]
```

At $\varphi=0$ the mean of $|z|$ is one (the state has collapsed to a pole), at $\varphi=\pi/4$ it drops to an intermediate value (localized but pulled toward the poles), and at $\varphi=\pi/2$ it settles near the value a uniform point on the great circle would give, $2/\pi\approx 0.64$: the axes have stopped imprinting and the walk has gone isotropic. This is exactly the transition Hacohen-Gourgy and co-workers watched on a transmon coupled to two cavity modes, one amplifier per mode (Nature 538, 491 (2016), arXiv:1608.06652): as they rotated the two measured axes from parallel to orthogonal, the reconstructed motion passed from collapse through localized diffusion to a uniform walk on the sphere, precisely this smooth handover in $\varphi$.

This stepper carries two approximations: the interleaved split of the kicks, and the sampling, which draws both readings from the state at the top of the step instead of conditioning each on the other's kick. So give Part II the same cross-check Part I got in Section 7, and derive the reference first. Averaged over its readings, one $x$-kick multiplies $\langle z\rangle$ by exactly $e^{-\Gamma\Delta t}$: weight each reading by the Section 9 base measure and integrate:

```wl
Block[{x, y, z, ib, dvar},
 With[{mx = kraus[ib, PauliMatrix[1], dvar]},
  FullSimplify[
   Integrate[Tr[PauliMatrix[3] . (mx . blochRho[x, y, z] . mx)] gaussLike[0, dvar, ib] Exp[-1/(2 dvar)],
     {ib, -Infinity, Infinity}, Assumptions -> {dvar > 0, {x, y, z} \[Element] Reals}] ==
    z Exp[-1/(2 dvar)], Assumptions -> {dvar > 0}]]]
```

This is Section 7's Gaussian-overlap integral taken in the $x$ basis, closing with no small-step expansion, so the unconditional mean from $|0\rangle$ at $\varphi=\pi/2$ must decay as exactly $e^{-\Gamma t}$. Measure the stacked ensemble's deviation from that exact value at a fixed time, sweeping the step from reckless to careful, keyed by the step:

```wl
With[{gm = 1., tf = 1.2, ntr = 20000},
 AssociationMap[
  Function[dt, BlockRandom[SeedRandom[11];
    With[{fin = Nest[ensembleStepNC[Pi/2, 1./(2 gm dt)],
        ConstantArray[{{1., 0.}, {0., 0.}}, ntr], Round[tf/dt]]},
     Mean[Re[fin[[All, 1, 1]] - fin[[All, 2, 2]]]] - Exp[-gm tf]]]],
  {0.6, 0.3, 0.1, 0.04, 0.02, 0.01}]]
```

At $\Gamma\Delta t$ of order one the deviation is a large fraction of the exact value itself: each "weak" look is nearly projective, and interleaved near-projections are not simultaneous monitoring. Deep inside the window the deviation shrinks as the step does, consistent with a bias first order in $\Gamma\Delta t$, and the leading source is the interleaving itself. To see it, integrate the sandwiched step operator against the base measure, once per reading, and look at what comes back:

```wl
Block[{ibz, ibp, dvar},
 With[{m = stepOp[ibz, ibp, Pi/2, dvar]},
  FullSimplify[Integrate[gaussLike[0, dvar, ibz] gaussLike[0, dvar, ibp] Exp[-1/dvar] m . m,
    {ibz, -Infinity, Infinity}, {ibp, -Infinity, Infinity}, Assumptions -> dvar > 0]]]]
```

The sandwiched operators sum to a multiple of the identity; expand that prefactor in the small parameter $1/D=2\Gamma\Delta t$:

```wl
Block[{dvar}, Series[1/2 - Exp[-1/dvar]/2 + Exp[-1/(2 dvar)], {dvar, Infinity, 2}]]
```

It falls short of one at second order per step, so over a fixed time the miss accumulates to first order in $\Gamma\Delta t$: the interleaved readings are not quite a POVM, and even a perfectly sampled interleaved step lets $\langle z\rangle$ decay slightly more slowly than $e^{-\Gamma\Delta t}$. The split owns the miss; confirm the same integral with the two kicks applied one after the other instead:

```wl
Block[{ibz, ibp, dvar},
 Integrate[gaussLike[0, dvar, ibz] gaussLike[0, dvar, ibp] Exp[-1/dvar]
   kraus[ibz, PauliMatrix[3], dvar] . kraus[ibp, sigmaPhi[Pi/2], dvar] .
    kraus[ibp, sigmaPhi[Pi/2], dvar] . kraus[ibz, PauliMatrix[3], dvar],
  {ibz, -Infinity, Infinity}, {ibp, -Infinity, Infinity}, Assumptions -> dvar > 0]]
```

Exactly the identity: taken one after the other the two kicks form a legitimate POVM at any step size, and the interleaving is what buys the single shared $z$-draw at the price of a second-order miss per step. The remaining piece of the bias is the sampling, both readings drawn from the state at the top of the step rather than from the joint density of the pair. That sets the validity window for everything in this Part: keep $\Gamma\Delta t$ well below one; at the steps used here the offset sits at the percent level, comparable to the ensemble scatter of the plots. You are not locked into these choices: change the trajectory count or the step list in the sweep above and watch the deviation and its scatter move together.

## 13. The Uncertainty Floor and the Purity Ceiling

There is a reason the state can never settle once the axes differ. For it to stop, both detectors would have to fall silent at the same instant: the state would have to be a still point for both back-actions at once. A single measurement has two such points, its two eigenstates; two non-commuting measurements share none, and the uncertainty principle is exactly the statement that the total measurement disturbance is bounded below by their commutator. For a qubit the disturbance the two channels deliver per step is proportional to the variance sum $\Delta\sigma_z^2+\Delta\sigma_\varphi^2$, and because $A^2=\mathbb 1$ for any Pauli axis each variance is $\Delta A^2=1-\langle A\rangle^2$, so at $\varphi=\pi/2$ the sum is $(1-z^2)+(1-x^2)$. In short, two detectors watching incompatible axes have no state they can both leave alone, so something is always being kicked. Confirm, for any pure state, that this sum can never fall below the commutator magnitude $|\langle[\sigma_z,\sigma_x]\rangle|=2|y|$:

```wl
Block[{x, y, z},
 FullSimplify[(2 - x^2 - z^2) - 2 Abs[y] >= 0,
  Assumptions -> {x^2 + y^2 + z^2 == 1, x \[Element] Reals, y \[Element] Reals, z \[Element] Reals}]]
```

The inequality holds everywhere and, written as $(|y|-1)^2\ge 0$, saturates only at $|y|=1$. More telling is the left side by itself: on the sphere the variance sum equals $1+y^2$, never below one, because $\sigma_z$ and $\sigma_x$ share no eigenstate and so no state can still both detectors at once. The walk of Section 11 lives precisely on the great circle $y=0$, where the total disturbance sits at that minimum of one, the commutator bound gone slack but the kicking as inescapable as ever: the diffusion never stops because the disturbance cannot vanish, a consequence of the uncertainty principle rather than a numerical artifact.

A real detector loses some of the information it gathers, and that loss caps how pure the state can be kept. With efficiency $\eta<1$ per channel the coherence is bled faster than the two records sharpen it, so the isotropic walk settles onto a mixed state whose purity sits below one, set entirely by $\eta$. The loss dephases the state about a measured axis $A$: the part of $\rho$ diagonal in the $A$-basis stays, the off-diagonal part shrinks by $e^{-\gamma\Delta t}$. Because $A^2=\mathbb 1$, the diagonal-in-$A$ part of $\rho$ is $(\rho+A\rho A)/2$, so the whole channel is one product with the axis; confirm the two-projector form and this one-product form agree for the angled axis and any qubit:

```wl
Block[{x, y, z, phi},
 With[{A = sigmaPhi[phi], rho = blochRho[x, y, z]},
  With[{p = (IdentityMatrix[2] + A)/2, q = (IdentityMatrix[2] - A)/2},
   FullSimplify[p . rho . p + q . rho . q == (rho + A . rho . A)/2]]]]
```

The two projected pieces recombine into half the state plus half its conjugation by the axis, so the dephasing channel is $\rho\to e^{-\gamma\Delta t}\rho+(1-e^{-\gamma\Delta t})(\rho+A\rho A)/2$: the $A$-populations ride through untouched, and what shrinks is exactly the information the detector gathered but failed to deliver. Add this loss to each channel after the stacked Bayesian kick and read the steady purity for three efficiencies, keyed by $\eta$; the purity of every trajectory at once is the elementwise product of the stack with its own transpose:

```wl
With[{gm = 1., dt = 0.02, nst = 600, ntr = 150},
 With[{step = ensembleStepNC[Pi/2, 1./(2 gm dt)], zB = ConstantArray[N@PauliMatrix[3], ntr],
   xB = ConstantArray[N@sigmaPhi[Pi/2], ntr]},
  AssociationMap[
   Function[eta,
    With[{e = Exp[-gm (1./eta - 1.) dt]},
     With[{deph = Function[{rr, aB}, ((1 + e) rr + (1 - e) ensembleDot[aB, ensembleDot[rr, aB]])/2]},
      BlockRandom[SeedRandom[Round[1000 eta]];
       With[{fin = Nest[deph[deph[step[#], zB], xB] &,
           ConstantArray[{{0.5, 0.}, {0., 0.5}}, ntr], nst]},
        Mean[Re[Total[fin Transpose[fin, {1, 3, 2}], {2, 3}]]]]]]]], {1., 0.64, 0.36}]]]
```

The steady purity is the ceiling the efficiency allows: only the ideal detector holds a fully pure state, and below $\eta=1$ the ceiling falls, the two records no longer able to keep the state pure. Through the identity $\mathrm{Tr}\,\rho^2=(1+|\mathbf r|^2)/2$, this is a statement about the Bloch length: the two records feed purity in at a rate proportional to $(1-r^2)(2-r^2)$, while the information the detector fails to deliver drains it at a rate proportional to $(1/\eta-1)\,r^2$, and the length settles where the two rates balance; we read both rates off the ensemble-averaged step without deriving them here, and let the simulated ceiling be the check. A steady, efficiency-capped Bloch length of exactly this kind is what the two-mode experiment reconstructed at its own efficiencies. The lost information is the gap between monitoring an actual wavefunction and monitoring only its blurred image.

The balance is quadratic in $r^2$, so let `Solve` settle it on the physical interval: it returns a single root $u=r^2$ between zero and one, the steady purity $(1+u)/2$; encode that root and confirm it is `Solve`'s root:

```wl
ClearAll[puritySS];
puritySS[eta_] := (1 + (1 + 2 eta - Sqrt[1 - 4 (eta - 1) eta])/(2 eta))/2;
Block[{u, eta},
 FullSimplify[
  puritySS[eta] ==
   (1 + (u /. First[Solve[{(1 - u) (2 - u) == (1/eta - 1) u, 0 <= u <= 1}, u, Reals]]))/2,
  0 < eta < 1]]
```

`Solve` finds exactly one root on the interval and `puritySS` is identically that root (the discarded root has $r^2>1$, a Bloch vector outside the sphere). The root carries the condition $\eta<1$; at $\eta=1$ the drain vanishes, the balance collapses to $(1-r^2)(2-r^2)=0$, and the closed form continues through the boundary with $r=1$, the ideal detector holding the state on the sphere. Evaluate the root at the three efficiencies just simulated, keyed by $\eta$:

```wl
AssociationThread[{1., 0.64, 0.36}, puritySS[{1., 0.64, 0.36}]]
```

The simulated purities of the previous cell land within about a percent of these roots, and the residue at this ensemble size is scatter: the balance owns the mechanism and the whole $\eta$-dependence.

## 14. The Output Correlators and the Cosine at Zero Lag

Now the payoff. The state itself is hard to see, but the two records $I_z(t)$ and $I_\varphi(t)$ are right there, and their temporal correlators carry the physics in closed form. Atalaya, Hacohen-Gourgy, Martin, Siddiqi, and Korotkov derived the self- and cross-correlators of the two output signals directly from the Bayesian equations and measured them on the two-mode apparatus (arXiv:1702.08077); for equally strong ideal channels with no extra evolution they are sums of two decaying exponentials with rates $\Gamma(1\mp\cos\varphi)$. Encode the closed forms:

```wl
ClearAll[kzzCF, kzpCF];
kzzCF[phi_, gm_, tau_] := (1/2)(1 + Cos[phi]) Exp[-gm (1 - Cos[phi]) tau] +
   (1/2)(1 - Cos[phi]) Exp[-gm (1 + Cos[phi]) tau];
kzpCF[phi_, gm_, tau_] := (1/2)(Exp[-gm (1 - Cos[phi]) tau] - Exp[-gm (1 + Cos[phi]) tau]) +
   (Cos[phi]/2)(Exp[-gm (1 - Cos[phi]) tau] + Exp[-gm (1 + Cos[phi]) tau])
```

Read the two limits these forms make vivid: the self-correlator starts at one and the cross-correlator starts at the cosine of the angle between the axes:

```wl
{FullSimplify[kzzCF[phi, gm, 0]], FullSimplify[kzpCF[phi, gm, 0]]}
```

The self-correlator at zero lag is one and the cross-correlator at zero lag is $\cos\varphi$: aligned detectors ($\varphi=0$) see perfectly correlated records, orthogonal detectors ($\varphi=\pi/2$) see uncorrelated ones, and the cosine interpolates. In the experiment this $\cos\varphi$ was traced across eleven angles from the raw records, a purely dynamical readout of the angle between two incompatible measurements, obtained without ever reconstructing the state. Recover the closed forms from the records themselves: run an ensemble of angled-measurement records, form the ensemble-and-time-averaged two-time correlators, and overlay the closed forms. First a run that advances the stacked ensemble while emitting, at every step, the two vectors of readings that drive it, carrying a slow-rotation slot as the same symmetric split (held at zero until the next section), together with the two-time average, a pair of shifted slices multiplied elementwise:

```wl
ClearAll[ensembleRecordsNC, corr];
ensembleRecordsNC[phi_, omR_, dt_, dvar_, nst_, ntr_, seed_] := BlockRandom[SeedRandom[seed];
  With[{sz = N@PauliMatrix[3], sp = N@sigmaPhi[phi],
    u = MatrixExp[-I (omR/2) PauliMatrix[2] (dt/2)]},
   With[{uB = Developer`ToPackedArray@ConstantArray[u, ntr],
     ubB = Developer`ToPackedArray@ConstantArray[ConjugateTranspose[u], ntr]},
    FoldPairList[
     Function[{rr, k}, With[{r1 = ensembleDot[uB, ensembleDot[rr, ubB]]},
       With[{ibz = ensembleDraw[r1, sz, dvar], ibp = ensembleDraw[r1, sp, dvar]},
        {{ibz, ibp}, ensembleDot[uB, ensembleDot[ensembleKick[r1, ibz, ibp, phi, dvar], ubB]]}]]],
     ConstantArray[{{0.5, 0.}, {0., 0.5}}, ntr], Range[nst]]]]];
corr[a_, b_, j_] := Total[a[[All, ;; -1 - j]] b[[All, 1 + j ;;]], 2]/(Length[a] (Dimensions[a][[2]] - j))
```

`FoldPairList` threads the state stack forward while peeling the two reading vectors off at each step, so the records and the states share the very same noise, which is what makes the cross-correlator carry the back-action; `corr` averages the lag-$j$ product over every trajectory and every start time at once. Now simulate two angles and compare with the closed form:

```wl
With[{gm = 1., dt = 0.02, nst = 250, ntr = 1000, nlag = 10, phis = {Pi/4, Pi/2}},
 With[{taus = Range[nlag] dt},
  With[{simCorr = Table[
      With[{recs = ensembleRecordsNC[phi, 0, dt, 1./(2 gm dt), nst, ntr, 700 + Round[100 phi]]},
       With[{rz = Transpose[recs[[All, 1]]], rp = Transpose[recs[[All, 2]]]},
        {Table[corr[rz, rz, j], {j, nlag}], Table[corr[rz, rp, j], {j, nlag}]}]],
      {phi, phis}]},
   Show[
    ListPlot[Catenate[Table[{Transpose[{taus, simCorr[[i, 1]]}], Transpose[{taus, simCorr[[i, 2]]}]}, {i, Length[phis]}]],
     Frame -> True, GridLines -> Automatic, PlotRange -> {{0, nlag dt}, {-0.25, 1.25}},
     FrameLabel -> {"lag \[Tau] (units of 1/\[CapitalGamma])", "output correlator"},
     PlotLegends -> {"K_zz, \[Phi]=\[Pi]/4", "K_z\[Phi], \[Phi]=\[Pi]/4", "K_zz, \[Phi]=\[Pi]/2", "K_z\[Phi], \[Phi]=\[Pi]/2"},
     PlotLabel -> "simulated records recover the closed-form correlators", ImageSize -> 520],
    Plot[Evaluate@Flatten@Table[{kzzCF[phi, gm, tau], kzpCF[phi, gm, tau]}, {phi, phis}], {tau, 0, nlag dt},
     PlotStyle -> Directive[Gray, Dashed]]]]]]
```

Each simulated correlator tracks its closed-form curve, with the points scattering about it, correlated across neighbouring lags, at this ensemble size: the self-correlators decay from one, the $\varphi=\pi/4$ cross-correlator starts near $\cos(\pi/4)\approx 0.71$ and the $\varphi=\pi/2$ cross-correlator starts near zero, the cosine law read straight off the records. The angle between the two incompatible measurements is written directly in how their noisy outputs correlate.

## 15. A Hidden Rotation, and the Bacon-Shor Connection

One last use of the correlators, and one last connection. If a small unwanted rotation $\tilde\Omega_R$ about $y$ is present, it breaks the symmetry between the two orderings and shows up only in the antisymmetrized cross-correlator $K_{z\varphi}(\tau)-K_{\varphi z}(\tau)$, which is otherwise zero. For the equally strong channels of the experiment the two decay rates coincide at $\varphi=\pi/2$, and that antisymmetrized combination reduces to $2\sin\varphi\,e^{-\Gamma\tau}\sin(\tilde\Omega_R\tau)$, linear in $\tilde\Omega_R$ for a slow rotation. The rotation slot in `ensembleRecordsNC` was built for this moment: simulate at $\varphi=\pi/2$ with a rotation below the measurement rate, form the antisymmetrized cross-correlator, and overlay the closed form:

```wl
With[{gm = 1., dt = 0.02, nst = 250, ntr = 2000, nlag = 100, omR = 0.5, phi = Pi/2},
 With[{taus = Range[nlag] dt},
  With[{recs = ensembleRecordsNC[phi, omR, dt, 1./(2 gm dt), nst, ntr, 400]},
   With[{rz = Transpose[recs[[All, 1]]], rp = Transpose[recs[[All, 2]]]},
    Show[
     ListPlot[Transpose[{taus, Table[corr[rz, rp, j] - corr[rp, rz, j], {j, nlag}]}],
      Frame -> True, GridLines -> Automatic, PlotRange -> {{0, nlag dt}, {-0.7, 0.9}},
      FrameLabel -> {"lag \[Tau] (units of 1/\[CapitalGamma])", "K_z\[Phi] - K_\[Phi]z"},
      PlotLabel -> "the antisymmetrized cross-correlator reads the hidden rotation", ImageSize -> 500],
     Plot[2 Sin[phi] Exp[-gm tau] Sin[omR tau], {tau, 0, nlag dt},
      PlotStyle -> Directive[Gray, Dashed]]]]]]]
```

The simulated difference rises from zero, follows the dashed $2\sin\varphi\,e^{-\Gamma\tau}\sin(\tilde\Omega_R\tau)$ through its maximum, and decays with the correlation time; the points scatter about the curve at this ensemble size, and the scatter is correlated across neighbouring lags, so a stretch of them can wander off together while the residuals average to zero. A rotation with no visible trace in either record alone stands out cleanly in the ordering asymmetry of the two records. This is how Atalaya and co-workers read a residual Rabi frequency of 12 kHz, some three thousand times below the 40 MHz drive, out of the records of the two-mode experiment.

The whole construction has a direct afterlife in error correction, and it is the subject of a paper of my own (Atalaya, Bahrami, Pryadko, Korotkov, arXiv:1612.02096). A four-qubit Bacon-Shor code measures four two-qubit gauge operators continuously; in the code space the four-qubit state stays $a(t)|z{+}\rangle+b(t)|z{-}\rangle$, so the measurement drives a single spare degree of freedom, a gauge qubit, while leaving the logical qubit untouched. In that two-dimensional subspace the two $Z$-type gauge operators are each a $z$-measurement of the gauge qubit and the two $X$-type ones are each an $x$-measurement, so the gauge qubit is watched along two orthogonal axes at once, exactly the doubly-watched qubit of this Part: it diffuses on a great circle while the logical qubit rides along. The error syndrome is then read not from a projective parity but from the time-averaged cross-correlators of Section 14, which for the ideal gauge qubit take the two-time form $e^{-2\Gamma_m|t_1-t_2|}$ with same-time value one; a single-qubit error moves the state into a subspace where a cross-correlator flips sign, and continuous monitoring performs about as well as projective syndrome extraction once the collapse time is roughly an order of magnitude shorter than the projective measurement period. The one-qubit cross-correlator we pulled out of simulated noise is, in that setting, a continuously monitored error check.

The other generalization keeps a single detector but gives it two qubits to watch at once. If the detector is built so it cannot tell $|01\rangle$ from $|10\rangle$, the same Bayes kick that localized one qubit now entangles two, with an odd-parity coherence riding along exactly as the single-qubit coherence did in Part I.

## 16. From One Qubit to Two: The Half-Parity Meter

Two qubits, remote and never coupled, are read by a single probe that reflects off both, so its output records the pair. In the computational basis $\{|00\rangle,|01\rangle,|10\rangle,|11\rangle\}$ the probe is engineered so that $|01\rangle$ and $|10\rangle$ shift it by the same amount: a half-parity meter, sorting the four states into three bins, $|00\rangle$, $|11\rangle$, and the odd pair $\{|01\rangle,|10\rangle\}$, never resolving which of the two odd states it holds. The window-averaged reading $V$ is Gaussian about a center that takes three values, $-\delta v$ for $|00\rangle$, $+\delta v$ for $|11\rangle$, and zero for both odd states. The reading carries white noise of single-sided density $s$, so a window $\Delta t$ leaves a Gaussian of variance $D=s/(2\Delta t)$ exactly as in Part I, and the measurement strength is $\delta v^2/s$, the two-qubit counterpart of the rate $\Gamma_m$; one extra rate $\gamma$ carries whatever dephasing survives inside the odd subspace. With the fast-damping off-diagonal elements dropped, the state is the X-state, five appreciable numbers $\{x_1,\dots,x_5\}=\{\rho_{00,00},\rho_{01,01},\rho_{10,10},\rho_{11,11},|\rho_{01,10}|\}$.

Encode the three-center Gaussian likelihoods for a reading $V$ and the concurrence of an X-state:

```wl
ClearAll[likelis, concurrence];
likelis[v_, dvv_, dvar_] := Exp[-({-dvv, 0, 0, dvv} - v)^2/(2 dvar)];
concurrence[x_] := 2 Max[0., x[[5]] - Sqrt[x[[1]] x[[4]]]]
```

The concurrence runs from zero for a separable or classically mixed state to one for a Bell state, set by a race: the odd-parity coherence $x_5=|\rho_{01,10}|$ pulling it up against the geometric mean $\sqrt{x_1 x_4}$ of the two even populations pulling it down. This is the very quantity the remote-entanglement experiment of Section 20 reports.

That two-line formula is the Wootters concurrence specialized to the X-state, and it is worth seeing it come out. Wootters' construction takes the square roots of the eigenvalues of $\rho\tilde\rho$, with $\tilde\rho=(\sigma_y\otimes\sigma_y)\rho^*(\sigma_y\otimes\sigma_y)$, sorts them, and subtracts the rest from the largest; compute those roots for a general X-state:

```wl
Block[{x1, x2, x3, x4, x5},
 With[{xr = {{x1, 0, 0, 0}, {0, x2, x5, 0}, {0, x5, x3, 0}, {0, 0, 0, x4}},
   sy = KroneckerProduct[PauliMatrix[2], PauliMatrix[2]]},
  FullSimplify[Sqrt[Eigenvalues[xr . sy . Conjugate[xr] . sy]],
   Assumptions -> {x1 > 0, x2 > 0, x3 > 0, x4 > 0, x5 > 0}]]]
```

The four roots pair up: $\sqrt{x_1x_4}$ twice from the even block, and $\sqrt{x_2x_3}\pm x_5$ from the odd block (the last two entries are those squares under the root). Two orderings put the odd block on top. First, $x_5\le\sqrt{x_2x_3}$, the positivity of the odd block, preserved because the ride-along ratio never grows. Second, $\sqrt{x_1x_4}\le\sqrt{x_2x_3}$: the even centers sit a full $\delta v$ from the odd ones, so every reading multiplies $x_1x_4/(x_2x_3)$ by $e^{-\delta v^2/D}<1$, and the meter starts the two products equal; Section 17 proves both step laws. Under those two orderings, Wootters' combination, twice the largest root minus their total, must collapse to the two-line formula; confirm the telescoping:

```wl
Block[{x1, x2, x3, x4, x5},
 With[{lam = {Sqrt[x1 x4], Sqrt[x1 x4], Sqrt[x2 x3] + x5, Sqrt[x2 x3] - x5}},
  FullSimplify[Max[0, 2 Max[lam] - Total[lam]] == 2 Max[0, x5 - Sqrt[x1 x4]],
   Assumptions -> {x1 > 0, x2 > 0, x3 > 0, x4 > 0, x5 > 0, x5 <= Sqrt[x2 x3], x1 x4 <= x2 x3}]]]
```

The identity holds on the whole ordered region: the race between the odd coherence and the even populations is the full Wootters concurrence on every state this meter can produce.

Now the update, the single-qubit Bayes kick of Part I one level up: each of the four populations picks up its likelihood and renormalizes, and the odd-parity coherence rides on $\sqrt{P_2 P_3}$, the geometric mean of the two likelihoods of the states it links, with the extra dephasing factor. Split it as before, the deterministic update for a given reading and the honest step that draws the reading from the mixture the state predicts:

```wl
ClearAll[twoQubitUpdate, twoQubitStep];
twoQubitUpdate[x_, v_, dvv_, dvar_, gt_] := With[{p = likelis[v, dvv, dvar]},
   {x[[1]] p[[1]], x[[2]] p[[2]], x[[3]] p[[3]], x[[4]] p[[4]],
     x[[5]] Sqrt[p[[2]] p[[3]]] Exp[-gt]}/(x[[1 ;; 4]] . p)];
twoQubitStep[dvv_, ss_, gt_, dt_][x_] := With[{dvar = ss/(2 dt)},
   twoQubitUpdate[x,
    RandomChoice[Clip[x[[1 ;; 4]], {0., 1.}] -> {-dvv, 0., 0., dvv}] +
     RandomVariate[NormalDistribution[0, Sqrt[dvar]]], dvv, dvar, gt]]
```

That is the entire integrator: draw which bin the record leans toward with the current populations as weights, bury it under window noise, and fold the reading back into the five numbers, the four populations by classical Bayes and the odd coherence by the geometric-mean ride-along.

## 17. The Odd-Parity Coherence Rides Along

The whole entangling mechanism is the purity ride-along of Section 3 promoted to the odd-parity subspace. There, because the coherence multiplied by the geometric mean of the two population likelihoods, the ratio $|\rho_{01}|^2/(\rho_{00}\rho_{11})$ was untouched by every ideal update, so a pure state stayed pure. Here the same ratio $|\rho_{01,10}|^2/(\rho_{01,01}\rho_{10,10})$ is untouched, so a coherent superposition inside the odd subspace stays coherent while the meter localizes the pair into it. Confirm three step laws on the update itself, for any reading and any $\gamma$: the four populations keep summing to one, the ride-along ratio picks up exactly $e^{-2\gamma\Delta t}$ (Section 3's law, one level up), and the even-to-odd population product falls by exactly $e^{-\delta v^2/D}$, the ordering Section 16's telescoping rests on:

```wl
Block[{x1, x2, x3, x4, x5, v, dvv, dvar, gt},
 With[{u = twoQubitUpdate[{x1, x2, x3, x4, x5}, v, dvv, dvar, gt]},
  FullSimplify[{Total[u[[1 ;; 4]]] == 1,
    u[[5]]^2/(u[[2]] u[[3]]) == Exp[-2 gt] x5^2/(x2 x3),
    (u[[1]] u[[4]])/(u[[2]] u[[3]]) == Exp[-dvv^2/dvar] (x1 x4)/(x2 x3)},
   Assumptions -> {x1 > 0, x2 > 0, x3 > 0, x4 > 0, x5 > 0, x1 + x2 + x3 + x4 == 1,
     dvar > 0, {v, dvv} \[Element] Reals, gt >= 0}]]]
```

All three hold identically: the update is trace-preserving, the ratio is conserved for an ideal meter and bled at the fixed rate $2\gamma$ for a lossy one, whatever the record says, and the even product drains steadily against the odd one. So the odd-parity coherence rides its two populations exactly as a single qubit's coherence rides its own, which is why a measurement that only sorts the pair into the odd bin leaves a fully coherent Bell state rather than a classical mixture of $|01\rangle$ and $|10\rangle$. Entanglement here is not created by any interaction; it is the coherence already present surviving the localization, because the meter cannot look inside the subspace it is projecting onto.

## 18. One Record Steered into a Bell State

Start the pair in a product of equal superpositions, so all four populations are one quarter and the odd coherence is one quarter. The meter sorts each record into one of three bins: if the record leans toward $-\delta v$ or $+\delta v$ the pair collapses toward the product state $|00\rangle$ or $|11\rangle$ and the odd populations drain away; if it stays near zero the pair is steered into the odd subspace, where the surviving coherence makes it the Bell state $(|01\rangle+|10\rangle)/\sqrt2$. Watch one record that lands in the odd bin, tracking its concurrence, the two even populations, and the total odd population:

```wl
With[{dvv = 1., ss = 2., gam = 0.15, dt = 0.02, nst = 300},
 With[{tt = Range[0, nst] dt,
   traj = BlockRandom[SeedRandom[4];
     NestList[twoQubitStep[dvv, ss, gam dt, dt], {0.25, 0.25, 0.25, 0.25, 0.25}, nst]]},
  ListLinePlot[{Transpose[{tt, concurrence /@ traj}],
    Transpose[{tt, traj[[All, 1]]}], Transpose[{tt, traj[[All, 4]]}],
    Transpose[{tt, traj[[All, 2]] + traj[[All, 3]]}]},
   PlotRange -> {0, 1.02}, Frame -> True, GridLines -> Automatic,
   PlotLegends -> {"concurrence", "\[Rho]00,00", "\[Rho]11,11", "odd population"},
   FrameLabel -> {"time (units of measurement time)", "value"},
   PlotLabel -> "one record steered into the entangled subspace", ImageSize -> 480]]]
```

The two even populations drain toward zero while the odd population fills toward one, and the concurrence climbs and then eases back: the record has projected the pair into the odd subspace and the surviving coherence has become genuine entanglement. A different seed would drain the odd population instead and collapse to a product state with zero concurrence; the outcome is stochastic, set by which bin the noisy record favors, with the bin probabilities equal to the initial populations, the Born rule again.

## 19. The Closed-Form Concurrence Ceiling

How entangled can a record make the pair, and how does that depend on when we look? The key premise is that the whole X-state is fixed by the running readout, and it follows from how Gaussian likelihoods compose. Before its normalization the update just multiplies the five numbers by the five per-reading factors $\{P_1,P_2,P_3,P_4,\sqrt{P_2P_3}\}$, so two steps multiply by the product of two factor sets; confirm that product equals the factor set at the averaged reading with half the variance, up to one common constant, by showing the elementwise ratio of the two is the same in every slot:

```wl
Block[{v1, v2, dvv, dvar},
 With[{fac = Function[{v, dv}, With[{p = likelis[v, dvv, dv]},
      {p[[1]], p[[2]], p[[3]], p[[4]], Sqrt[p[[2]] p[[3]]]}]]},
  FullSimplify[Equal @@ (fac[v1, dvar] fac[v2, dvar]/fac[(v1 + v2)/2, dvar/2]),
   Assumptions -> {dvar > 0, {v1, v2, dvv} \[Element] Reals}]]]
```

The ratio is slot-independent, and a common constant dies in the normalization: two steps at readings $V_1,V_2$ are one step at their average with half the variance, so by induction the state after any record depends on it only through the running mean $V$, with the dephasing factors multiplying up on their own. The concurrence at any time is therefore a definite function of $V$, and it is largest for the record that sits exactly at $V=0$, the most balanced odd-bin outcome. That most-fortunate concurrence is a closed form, given by the same quantum-Bayesian analysis the experiment of Section 20 is analyzed with: a race between the extra dephasing $\gamma$ that erodes the coherence and the measurement rate $\delta v^2/s$ that first builds and then over-localizes it. Encode the maximum concurrence and the concurrence-readout relation:

```wl
ClearAll[cmaxCF, cOfVCF];
cmaxCF[gam_, dvv_, ss_, t_] := (Exp[-gam t] - Exp[-dvv^2 t/ss])/(1 + Exp[-dvv^2 t/ss]);
cOfVCF[gam_, dvv_, ss_, vt_, t_] := (Exp[-gam t] - Exp[-dvv^2 t/ss])/(1 + Cosh[2 vt dvv t/ss] Exp[-dvv^2 t/ss])
```

Show that the concurrence-readout relation is largest at $V=0$ and equals the closed-form ceiling there:

```wl
{FullSimplify[cOfVCF[gam, dvv, ss, 0, t] == cmaxCF[gam, dvv, ss, t]],
 FullSimplify[(D[cOfVCF[gam, dvv, ss, vt, t], vt] /. vt -> 0) == 0]}
```

Both are true: the concurrence is stationary in the readout at $V=0$ and the value there is the ceiling, because $\cosh$ is smallest at zero and sits in the denominator. No record, however lucky, can beat the $V=0$ outcome.

Now the payoff, and the five-number update stacks over an ensemble exactly as the two-axis step did: the state becomes an $n\times 5$ array, the bin draw becomes one vector of uniform draws sorted by two `UnitStep` thresholds, and the two odd centers coincide, so the geometric mean $\sqrt{P_2P_3}$ is simply the shared zero-center likelihood. Encode the stacked run:

```wl
ClearAll[twoQubitKick, twoQubitEnsemble];
twoQubitKick[xx_, v_, dvv_, dvar_, e_] := With[
   {pm = Exp[-(v + dvv)^2/(2 dvar)], pz = Exp[-v^2/(2 dvar)], pp = Exp[-(v - dvv)^2/(2 dvar)]},
   Transpose[{xx[[All, 1]] pm, xx[[All, 2]] pz, xx[[All, 3]] pz, xx[[All, 4]] pp, xx[[All, 5]] pz e}]/
    (xx[[All, 1]] pm + (xx[[All, 2]] + xx[[All, 3]]) pz + xx[[All, 4]] pp)];
twoQubitEnsemble[dvv_, ss_, gt_, dt_, ntr_, nst_, seed_] := BlockRandom[SeedRandom[seed];
  With[{dvar = ss/(2 dt), e = Exp[-gt]},
   NestList[
    Function[xx, With[
      {w = Clip[xx[[All, 1 ;; 4]], {0., 1.}], u = RandomReal[1., Length[xx]]},
      twoQubitKick[xx,
       dvv (UnitStep[u - w[[All, 1]] - w[[All, 2]] - w[[All, 3]]] - UnitStep[w[[All, 1]] - u]) +
        RandomVariate[NormalDistribution[0, Sqrt[dvar]], Length[xx]], dvv, dvar, e]]],
    ConstantArray[{0.25, 0.25, 0.25, 0.25, 0.25}, ntr], nst]]]
```

One uniform draw per trajectory picks the bin, falling below the $|00\rangle$ weight or clearing everything but the $|11\rangle$ weight, the drawn center is buried under the window noise, and the five columns then update elementwise, the same Bayes-and-ride-along arithmetic acting on the whole ensemble at once. Pin the stacked kick to the update the symbolic checks ran on, a stack of one state over a grid of readings, strengths, and variances:

```wl
With[{x = {0.21, 0.34, 0.19, 0.26, 0.24}},
 Max @ Table[Max[Abs[First[twoQubitKick[{x}, {v}, dvv, dvar, Exp[-0.003]]] -
      twoQubitUpdate[x, v, dvv, dvar, 0.003]]],
   {v, {-1.2, 0.1, 0.9}}, {dvv, {0.7, 1.}}, {dvar, {2., 50.}}]]
```

Agreement at machine precision again: the trace, ride-along, product-drain, and composition laws proved for `twoQubitUpdate` are properties of the stacked run as well. Run it from the product start and plot every trajectory's concurrence, one `Ramp` over the whole stack, against the closed-form ceiling:

```wl
With[{dvv = 1., ss = 2., gam = 0.15, dt = 0.02, nst = 300, ntr = 60},
 With[{tt = Range[0, nst] dt, traj = twoQubitEnsemble[dvv, ss, gam dt, dt, ntr, nst, 900]},
  With[{cmat = 2 Ramp[traj[[All, All, 5]] - Sqrt[traj[[All, All, 1]] traj[[All, All, 4]]]]},
   Show[
    ListLinePlot[Transpose[{tt, #}] & /@ Transpose[cmat], PlotStyle -> Directive[Opacity[0.35], Gray],
     PlotRange -> {0, 0.45}, Frame -> True, GridLines -> Automatic,
     FrameLabel -> {"time (units of measurement time)", "concurrence"},
     PlotLabel -> "records press against the closed-form ceiling", ImageSize -> 500],
    Plot[cmaxCF[gam, dvv, ss, t], {t, 0, nst dt}, PlotStyle -> Directive[Black, Thick, Dashed]]]]]]
```

Every trajectory stays below the dashed ceiling, and the luckiest ones ride right along it: the closed-form maximum concurrence is exactly the sharp cutoff of the simulated distribution, recovered from records without ever imposing it. The ceiling rises from zero as the measurement builds coherence into concurrence, peaks, and then decays as continued localization and the extra dephasing win, so there is an optimal moment to stop measuring. Look at where the trajectories end up, at a fixed late time, across a large ensemble:

```wl
With[{dvv = 1., ss = 2., gam = 0.15, dt = 0.02, nst = 300, ntr = 2000},
 With[{fin = Last[twoQubitEnsemble[dvv, ss, gam dt, dt, ntr, nst, 7]]},
  Histogram[2 Ramp[fin[[All, 5]] - Sqrt[fin[[All, 1]] fin[[All, 4]]]], 40, Frame -> True, GridLines -> Automatic,
   FrameLabel -> {"final concurrence", "counts"},
   PlotLabel -> "bimodal: a spike at zero and a peak beneath the ceiling", ImageSize -> 460]]]
```

The distribution is bimodal: a tall spike at zero concurrence from the records that collapsed to $|00\rangle$ or $|11\rangle$, and a separate peak of entangled outcomes that pushes up against, but never past, the ceiling. The two peaks are the two fates of the pair, a product state or the Bell subspace, and the gap between them is why the measurement must be conditioned, kept only when the record lands in the odd bin.

## 20. What Sets the Ceiling, and the Remote-Qubit Experiment

The height of the ceiling is a trade-off between the extra dephasing $\gamma$ and the measurement rate. A slow, clean measurement builds coherence into concurrence before losing it; a noisy or lossy one erodes the coherence faster than it can localize the pair, and the peak concurrence drops. The trade-off is already written in the ceiling's three limits: a perfectly clean meter, long times, and short times:

```wl
Block[{gam, dvv, ss, t},
 {FullSimplify[cmaxCF[0, dvv, ss, t], {dvv > 0, ss > 0, t > 0}],
  Limit[cmaxCF[gam, dvv, ss, t], t -> Infinity, Assumptions -> {gam > 0, dvv > 0, ss > 0}],
  Normal @ Series[cmaxCF[gam, dvv, ss, t], {t, 0, 1}]}]
```

A clean meter's ceiling is $\tanh(\delta v^2 t/2s)$, climbing monotonically to a full Bell state; any $\gamma>0$ sends the ceiling to zero at long times; and the short-time slope is $(\delta v^2-\gamma s)/2s$, so the concurrence rises at all only when the measurement outruns the loss, $\delta v^2>\gamma s$. Between the rising start and the dying tail sits a peak; read it for three values of the extra dephasing, keyed by $\gamma$:

```wl
With[{dvv = 1., ss = 2., dt = 0.01, nst = 800},
 AssociationMap[Function[gam, Max[cmaxCF[gam, dvv, ss, Range[nst] dt]]], {0.05, 0.15, 0.4}]]
```

The peak concurrence falls as the extra dephasing grows: a nearly ideal channel reaches deep into the entangled regime, while a lossy one barely leaves the separable one. This is exactly the trade-off that capped the first remote-entanglement experiment (Roch, Schwartz, Motzoi, Macklin, Vijay, Eddins, Korotkov, Whaley, Sarovar, Siddiqi, Phys. Rev. Lett. 112, 170501 (2014), arXiv:1402.1868), where two transmons joined by 1.3 meters of coaxial cable and a shared probe reached a peak concurrence of 0.35 at a measurement efficiency near 0.4, before intrinsic dephasing and the finite efficiency of the amplifier stopped the climb; with better hardware the same analysis projects about 0.7.

The update we ran is the same one that experiment is analyzed with, and it says so explicitly, confirming the validity of the quantum-Bayesian formalism for a cascaded system. Two qubits in separate cavities, probed in sequence so the signal reflects off one and then the other, are a cascaded system whose joint measurement is exactly this half-parity meter; the four populations follow classical Bayes and the odd-parity coherence rides the geometric mean, which is the whole content of `twoQubitStep`, and the concurrence they report is the very $C=2\max(0,|\rho_{01,10}|-\sqrt{\rho_{00,00}\rho_{11,11}})$ we built in Section 16. The most-likely-path analysis of their ensemble splits into three branches, the high-concurrence odd branch hugging the ceiling and the two low-concurrence branches collapsing to the product states $|00\rangle$ and $|11\rangle$, matching the bimodal histogram. And Part II is the same physics seen the other way: there one qubit is watched along two incompatible axes and diffuses on a great circle; here two qubits are watched along one parity axis and are steered into a Bell state. In both, the state is moved not by any Hamiltonian but by the information a noisy record carries, folded back by Bayes' rule.

## Where This Leaves Us (and What Comes Next)

You now have a computation-first toolkit for continuous quantum measurement built from one Gaussian measurement operator and run in three settings. One qubit and one detector give the Bayes kick on the populations with the coherence riding along, the Rabi split, the ideality-capped purity, the Korotkov-Averin spectrum and its factor of four, and the phase back-action that completes the Kraus operator $M_{\bar I}=U_{\bar I}\sqrt{M_{\bar I}^\dagger M_{\bar I}}$. Two detectors on non-commuting axes turn that operator on twice and watch collapse give way to great-circle diffusion, with an uncertainty floor that forbids the walk from stopping and output correlators whose zero-lag cross-value reads the angle between the axes. One detector on two qubits with a blind spot turns the same operator into a half-parity meter that steers the pair into a Bell state, with a closed-form concurrence ceiling recovered from the records. Each generalization is the same update applied twice or applied to a larger space, and each has been carried out in a superconducting-qubit laboratory. Before closing, the load-bearing results computed along the way:

- Information is the shrinking overlap of two Gaussians; with no Hamiltonian the populations follow Gaussian Bayes and localize to a pole with the Born-rule probability.
- An ideal detector conserves $|\rho_{01}|^2/(\rho_{00}\rho_{11})$, so a pure state stays pure along a single record; the ensemble dephases only because runs localize differently.
- The output spectrum's coherent line is $4\tilde\eta$ above the floor, at most four and at most two for a phase-preserving amplifier, half of it the nonclassical noise-backaction correlation; a naive Ito stepper that drops the $(\Delta I)^2/4S$ drift loses the entire ensemble dephasing.
- Every measurement is one Kraus operator $M_{\bar I}=U_{\bar I}\sqrt{M_{\bar I}^\dagger M_{\bar I}}$, the phase kick times the informational Bayes factor.
- Two Kraus operators on non-commuting axes fail to commute by $2i\sinh(a/2D)\sinh(b/2D)\,\sigma_y$; orthogonal detectors keep a pure state pure while forcing great-circle diffusion, bounded below by the commutator, with a steady Bloch length set by the balance between record-driven purification and the undelivered information.
- The two output correlators are sums of exponentials with rates $\Gamma(1\mp\cos\varphi)$; at zero lag the cross-correlator is $\cos\varphi$, and its antisymmetrized part reads a hidden slow rotation, the same cross-correlators that serve as a Bacon-Shor code's error syndromes.
- A half-parity meter's odd-parity coherence rides the geometric mean of its two populations, the single-qubit ride-along one level up, so a measurement alone leaves a Bell state; the concurrence it reaches has a sharp closed-form ceiling and a bimodal final distribution.

The natural continuations are the pieces the essay still sets aside: a finite detector bandwidth, where the qubit and the cavity field must be tracked together; real-time feedback to lock the Rabi oscillation or stabilize the entanglement instead of only heralding it; and the many-qubit parity measurements that turn the single joint check of Part III into the continuously monitored error syndrome of a stabilizer code, where it meets the Bacon-Shor gauge qubit of Part II. This Bayesian-update view is also the sibling of the stochastic-differential-equation essays in this folder, `NDSolve vs Ito/` and `Manual StepLike Ito/`, which integrate the same physics as a Bloch-Langevin equation.
