---
Template: Default
---

# Watching Two Axes at Once: Non-Commuting Continuous Measurement in the Quantum-Bayesian Picture

**A single detector watching one axis of a qubit collapses it: the populations follow Bayes' rule and the coherence rides along in closed form, exactly the sibling essay's story. Point a second detector at a different axis and that closed form breaks, because the two Bayesian kicks are Gaussian operators that do not commute. This essay builds the general update from primitives and runs it: we watch the single-axis update stay diagonal, watch the two kicks fail to commute by a rotation about the third axis, and then let the general update take over. We watch the state stop collapsing and start diffusing, walk the measurement angle from collapse through localized diffusion to an isotropic random walk on the Bloch sphere that never lands, read a purity ceiling set by detector efficiency, and finally pull the two closed-form output correlators back out of simulated records, where the cross-correlator at zero lag reads the cosine of the angle between the two watched axes. A last step connects this one-qubit toy to the gauge qubit of a Bacon-Shor code, whose error syndromes are exactly these cross-correlators.**

Mads Bahrami (last updated: August 18, 2026)

## Setting the Stage: How This Essay Flows

This essay is a computation-first tour of one qubit watched by two linear detectors at the same time, along two different Bloch-sphere axes: how the single-axis Bayesian update is a closed form on the populations, why a second, non-commuting axis destroys that closed form, what the general update does instead, how the collapse gives way to persistent diffusion as the two axes separate, why the uncertainty principle forbids the diffusion from ever stopping, and how the two noisy output records carry the whole story in their temporal correlators. Every claim is computed on the spot from base Wolfram Language, so nothing here rests on a formula you cannot immediately test and change.

In other words, I have tried to build a small laboratory for simultaneous measurement of incompatible observables. I strongly believe in a computation-first narrative for learning: in a sense, if I cannot compute it, I cannot claim to understand it. The rhythm throughout is concept, then computation, then interpretation.

The environment you see is a live Wolfram notebook. Evaluate the cells from top to bottom; the toolkit defined in the next section is used by every section that follows, so those dependencies matter. Two of the cells run Monte-Carlo ensembles and take a minute or two; I flag them where they appear. My suggestion is to focus on the output and its meaning first, then unpack the input code. You are not locked into any of it: change the angle, the rates, the seed, and rerun your own experiments.

This is the non-commuting door into the same continuous-measurement physics that the core essay, `Korotkov-Quantum-Bayesian.md`, treats for a single symmetric detector. That essay is the prerequisite: it builds the Gaussian likelihood, the Bayes kick on the populations with the coherence riding along, and the purity ride-along that we here watch break. We recall its single-axis update in Section 1 and then generalize.

Let's start!

## Notation and Conventions: Two Detectors, Two Records, One Qubit

We work in units where the reduced Planck constant is one, $\hbar=1$. The qubit is a $2\times 2$ density matrix $\rho=(\openone + x\,\sigma_x + y\,\sigma_y + z\,\sigma_z)/2$, with real Bloch coordinates $(x,y,z)$. Two linear detectors watch two directions on the Bloch sphere at the same time: the first watches $\sigma_z$, the second watches $\sigma_\varphi=\sigma_z\cos\varphi + \sigma_x\sin\varphi$, an axis lying in the $xz$-plane at angle $\varphi$ from $z$. The two produce noisy output records,

$$ I_z(t)=\mathrm{Tr}[\sigma_z\rho(t)] + \sqrt{\tau}\,\xi_z(t), \qquad I_\varphi(t)=\mathrm{Tr}[\sigma_\varphi\rho(t)] + \sqrt{\tau}\,\xi_\varphi(t), $$

with $\xi_z,\xi_\varphi$ independent white noises, $\langle\xi_i(t)\,\xi_j(t')\rangle=\delta_{ij}\,\delta(t-t')$. We assume phase-sensitive amplifiers on the optimal quadrature, so there is no phase back-action, only the informational (quantum-Bayesian) kind. Two numbers per channel organize everything: the ensemble dephasing rate $\Gamma$ that measurement imposes, and the detector efficiency $\eta\in(0,1]$, related to the collapse time $\tau$ by $\eta=1/(2\tau\Gamma)$; an ideal (quantum-limited) detector has $\eta=1$, so $\tau=1/(2\Gamma)$. When we average a record over a step $\Delta t$, the reading has Gaussian window variance $D=\tau/\Delta t=1/(2\Gamma\Delta t)$ for an ideal detector. We take the two channels equally strong, $\Gamma_z=\Gamma_\varphi=\Gamma$, unless we say otherwise.

Fix the two measured directions and the readouts we will plot: the angled observable $\sigma_\varphi$, the Bloch vector, the purity, and the Bloch length:

```wl
ClearAll[sigmaPhi, blochVec, purity, rlen];
sigmaPhi[phi_] := Cos[phi] PauliMatrix[3] + Sin[phi] PauliMatrix[1];
blochVec[rho_] := Re@Table[Tr[rho . PauliMatrix[j]], {j, 3}];
purity[rho_] := Re@Tr[rho . rho];
rlen[rho_] := Norm[blochVec[rho]]
```

The Bloch length $|\mathbf r|=\sqrt{x^2+y^2+z^2}$ is one for a pure state and shrinks toward zero for a mixed one, so it is the single number that tells us whether the record is monitoring an actual wavefunction or only an ensemble.

Now the object the whole essay turns on. Recall from the core essay that measuring a $\pm 1$-valued observable and reading $\bar I$ over a step multiplies the state, in the eigenbasis of that observable, by the square roots of the two Gaussian likelihoods. That element-wise multiplication is exactly conjugation by a single operator: a Gaussian Kraus operator that is a matrix function of the measured Pauli. Encode it, for a reading $\bar I$ of a Pauli $A$ with window variance $D$:

```wl
ClearAll[kraus];
kraus[ibar_, A_, dvar_] := MatrixExp[(ibar/(2 dvar)) A]
```

Because $A^2=\openone$, this exponential is just $\cosh(\bar I/2D)\,\openone + \sinh(\bar I/2D)\,A$: a single Hermitian operator built from the measured axis, applied as $\rho\to M\rho M/\mathrm{Tr}(M\rho M)$. This one object is the whole update, and everything downstream is what happens when we build it out of two non-commuting axes.

Define the honest detector draw, sampling a reading from the two-outcome Gaussian mixture that the current state predicts for the measured observable:

```wl
ClearAll[drawRead];
drawRead[rho_, A_, dvar_] :=
  RandomChoice[Clip[{(1 + Re@Tr[A . rho])/2, (1 - Re@Tr[A . rho])/2}, {0., 1.}] -> {1., -1.}] +
   RandomVariate[NormalDistribution[0, Sqrt[dvar]]]
```

In words, the detector first picks which eigenvalue to ring, $+1$ or $-1$, with the Born weights $(1\pm\mathrm{Tr}[A\rho])/2$, then buries the level under Gaussian window noise of variance $D$; the two channels draw independently from the same state.

## 1. One Axis: The Closed Form We Are About to Break

Start with a single detector on $\sigma_z$, which is the core essay's whole story. Because $\sigma_z$ is diagonal, its Kraus operator is diagonal too, so the update touches the two populations and the two coherences separately: the populations pick up the Gaussian likelihood ratio and renormalize (classical Bayes), and the coherence rides along on the geometric mean. There is nothing to solve. Look at the single-axis Kraus operator:

```wl
kraus[ibar, PauliMatrix[3], dvar] // MatrixForm
```

As one can see, it is diagonal, $\mathrm{diag}(e^{\bar I/2D}, e^{-\bar I/2D})$: measuring one axis is a scaling of the two amplitudes in that axis's own basis, and the population ratio $\rho_{11}/\rho_{00}$ simply multiplies by $e^{\bar I/D}$. This is why the single-axis update has a closed form, and why a pure state stays pure along a single record. The next question is what a second, differently pointed detector does to this tidy diagonal picture.

## 2. Two Axes: Why the Closed Form Breaks

A second detector on $\sigma_\varphi$ contributes its own Gaussian Kraus operator, diagonal in the $\sigma_\varphi$ basis rather than the $\sigma_z$ basis. To apply both over one step we must multiply the two operators, and the order matters unless they commute. They do not, and the failure is not a small correction: it is the entire reason simultaneous non-commuting measurement is a distinct physical process. Compute the commutator of the two single-axis kicks for $\sigma_z$ and $\sigma_x$:

```wl
FullSimplify[kraus[a, PauliMatrix[3], dvar] . kraus[b, PauliMatrix[1], dvar] -
   kraus[b, PauliMatrix[1], dvar] . kraus[a, PauliMatrix[3], dvar]]
```

The commutator is nonzero and off-diagonal, so there is no basis in which both kicks are diagonal at once, and therefore no closed-form population update: the two detectors cannot be reduced to two independent Bayesian scalings. Confirm that the mismatch is exactly a kick about the third axis, $2i\sinh(a/2D)\sinh(b/2D)\,\sigma_y$:

```wl
FullSimplify[
 kraus[a, PauliMatrix[3], dvar] . kraus[b, PauliMatrix[1], dvar] -
   kraus[b, PauliMatrix[1], dvar] . kraus[a, PauliMatrix[3], dvar] ==
  2 I Sinh[a/(2 dvar)] Sinh[b/(2 dvar)] PauliMatrix[2],
  Assumptions -> {a \[Element] Reals, b \[Element] Reals, dvar > 0}]
```

The two kicks fail to commute by a rotation about $y$, inheriting the algebra $[\sigma_z,\sigma_x]=2i\sigma_y$ of the observables they measure. Because the mismatch is first order in each reading, it scales like the step, so the two orderings agree as the step shrinks; but at any finite step the joint measurement genuinely mixes the two channels, and we must carry the full $2\times 2$ matrix.

The general update takes over. It is the entire step: draw both readings from the current state, apply the two Bayesian kicks as a symmetric split (half the $z$-kick, the full $\varphi$-kick, half the $z$-kick), and renormalize. Define it:

```wl
ClearAll[measureStepNC];
measureStepNC[phi_, dvar_][rho_] := Module[
   {szA = PauliMatrix[3], spA = sigmaPhi[phi], ibz, ibp, mm, r1},
   ibz = drawRead[rho, szA, dvar]; ibp = drawRead[rho, spA, dvar];
   mm = kraus[ibz/2, szA, dvar] . kraus[ibp, spA, dvar] . kraus[ibz/2, szA, dvar];
   r1 = mm . rho . ConjugateTranspose[mm]; r1/Tr[r1]]
```

That is the whole integrator: two Gaussian draws, three matrix products for the symmetric split, and one renormalization, with no basis in which it factors into scalar Bayes rules. Confirm it keeps $\rho$ a physical state, trace one and Hermitian, after a step from a generic mixed state:

```wl
With[{r0 = {{0.6, 0.2 + 0.1 I}, {0.2 - 0.1 I, 0.4}}},
 With[{r1 = BlockRandom[SeedRandom[1]; measureStepNC[Pi/2, 5.][r0]]},
  {Chop@Tr[r1], HermitianMatrixQ[Chop@r1]}]]
```

As expected, the step preserves the trace and Hermiticity: the general update is a legitimate quantum operation even though it is not a closed-form Bayes rule.

## 3. Orthogonal Axes: The State Stops Collapsing and Starts Diffusing

Point the second detector at right angles to the first, $\varphi=\pi/2$, so the two watched axes are $z$ and $x$, maximally incompatible. Now each detector is forever trying to collapse the state onto its own axis while the other detector pulls it away, and neither wins: instead of localizing to a pole, the pure state performs a persistent random walk. On the ideal detectors it never leaves the surface of the Bloch sphere, because information gain alone does not decohere a single record, so the walk is along a great circle. Run one orthogonal-measurement trajectory from the fully mixed center and watch the Bloch coordinates:

```wl
With[{gm = 1., dt = 0.02, nst = 800},
 With[{dvar = 1./(2 gm dt), tt = Range[0, 800] 0.02},
  With[{traj = BlockRandom[SeedRandom[4];
       blochVec /@ NestList[measureStepNC[Pi/2, dvar], {{0.5, 0.}, {0., 0.5}}, nst]]},
   ListLinePlot[Transpose[{tt, #}] & /@ Transpose[traj],
    PlotRange -> {-1.05, 1.05}, Frame -> True, GridLines -> Automatic,
    PlotLegends -> {"x", "y", "z"}, FrameLabel -> {"time (units of 1/\[CapitalGamma])", "Bloch vector"},
    PlotLabel -> "orthogonal measurement: persistent diffusion, no collapse", ImageSize -> 460]]]]
```

The $z$ and $x$ coordinates wander over the whole range and never settle, while $y$ stays near zero: the state diffuses on the $xz$ great circle, driven back and forth by the two competing records. There is no pole to fall into, so the "collapse" of the single-axis case is replaced by a walk that never ends. See it as a path on the sphere by plotting the same trajectory in the $xz$-plane against the unit circle:

```wl
With[{gm = 1., dt = 0.02, nst = 800},
 With[{dvar = 1./(2 gm dt)},
  With[{traj = BlockRandom[SeedRandom[4];
       blochVec /@ NestList[measureStepNC[Pi/2, dvar], {{0.5, 0.}, {0., 0.5}}, nst]]},
   ListLinePlot[traj[[All, {1, 3}]], AspectRatio -> 1, Frame -> True, GridLines -> Automatic,
    PlotRange -> {{-1.1, 1.1}, {-1.1, 1.1}}, FrameLabel -> {"x", "z"},
    PlotLabel -> "the same record as a walk on the xz great circle",
    Epilog -> {Gray, Dashed, Circle[]}, ImageSize -> 380]]]]
```

The path hugs the unit circle, wandering around it without a preferred point: the two orthogonal detectors together monitor the wavefunction perfectly (it stays pure) while leaving its position on the great circle to diffuse freely. This is the isotropic persistent diffusion that the single-basis picture has no room for.

## 4. Turning the Angle: From Collapse to Localized Diffusion to a Free Walk

The angle $\varphi$ between the two watched axes tunes the whole character of the motion. At $\varphi=0$ the two detectors watch the same axis, their kicks commute (both diagonal), and the state collapses to a $z$-pole exactly as one strong measurement would. As $\varphi$ opens up, the state can no longer reach a single point, but the axes still leave an imprint: it localizes to regions near the two eigen-directions and jitters there. At $\varphi=\pi/2$ the imprint vanishes and the walk is uniform. Watch this by running an ensemble from the center for three angles and reading, at the end of each run, how far the state sits from the equator, the mean of $|z|$, keyed by the angle (this cell runs a small ensemble and takes a few seconds):

```wl
With[{gm = 1., dt = 0.03, nst = 250, ntr = 200},
 With[{dvar = 1./(2 gm dt)},
  AssociationThread[{0, Pi/4, Pi/2},
   Table[
    Mean[Abs@Table[BlockRandom[SeedRandom[Round[100 phi] + k];
        blochVec[Last@NestList[measureStepNC[phi, dvar], {{0.5, 0.}, {0., 0.5}}, nst]][[3]]], {k, ntr}]],
    {phi, {0., Pi/4, Pi/2}}]]]]
```

As one can see, at $\varphi=0$ the mean of $|z|$ is one (the state has collapsed to a pole), at $\varphi=\pi/4$ it drops to an intermediate value (localized but pulled toward the poles), and at $\varphi=\pi/2$ it settles near the value a uniform point on the great circle would give, $2/\pi\approx 0.64$: the axes have stopped imprinting and the walk has gone isotropic. The transition is smooth in $\varphi$, from collapse through localized diffusion to a free walk, and, as the great-circle path of Section 3 already showed, the Bloch length holds at one throughout because purity is untouched.

## 5. The Uncertainty Floor: Why the Diffusion Cannot Stop

There is a reason the state can never settle once the axes differ. For it to stop, both detectors would have to fall silent at the same instant, that is, the state would have to be a still point for both measurement back-actions at once. A single measurement has two such points, its two eigenstates. Two non-commuting measurements share none, and the uncertainty principle is exactly the statement that the total measurement disturbance is bounded below by their commutator. For a qubit the disturbance per step is $(\Delta\sigma_z^2 + \Delta\sigma_\varphi^2)\,\Gamma\eta\,\Delta t$, the sum of the two variances, and the bound is the magnitude of the commutator. Confirm, for any pure state, that the sum of variances can never fall below the commutator magnitude $|\langle[\sigma_z,\sigma_x]\rangle|=2|y|$:

```wl
FullSimplify[(2 - x^2 - z^2) - 2 Abs[y] >= 0,
 Assumptions -> {x^2 + y^2 + z^2 == 1, x \[Element] Reals, y \[Element] Reals, z \[Element] Reals}]
```

The inequality holds everywhere and, written as $(|y|-1)^2\ge 0$, saturates only at $|y|=1$: the disturbance vanishes only if the state points fully along $y$, which is a still point for neither the $z$ nor the $x$ detector but the direction that carries zero signal in both. Away from that single unreachable direction the disturbance is strictly positive, so the state is always being kicked, and the diffusion of Section 3 is a direct consequence of the uncertainty principle rather than a numerical artifact.

A real detector loses some of the information it gathers, and that loss caps how pure the state can be kept. With efficiency $\eta<1$ per channel the coherence is bled faster than the two records sharpen it, so the isotropic walk settles onto a sphere of radius less than one. The steady radius is set entirely by the efficiency. First encode the non-ideal loss as a pure-dephasing factor $e^{-\gamma\Delta t}$ on the coherence about a measured axis $A$:

```wl
ClearAll[dephaseAxis];
dephaseAxis[rho_, A_, gt_] := With[{p = (IdentityMatrix[2] + A)/2, q = (IdentityMatrix[2] - A)/2},
   p . rho . p + q . rho . q + Exp[-gt] (rho - p . rho . p - q . rho . q)]
```

This shrinks only the part of $\rho$ that is off-diagonal in the $A$-basis, leaving the $A$-populations alone, which is exactly the information the detector gathered but failed to deliver. Now add this loss to each of the two channels after the Bayesian kick, and read the steady Bloch length for three efficiencies, keyed by $\eta$ (this cell runs a small ensemble and takes a few seconds):

```wl
With[{gm = 1., dt = 0.02, nst = 600, ntr = 150},
 With[{dvar = 1./(2 gm dt)},
  AssociationThread[{1., 0.64, 0.36},
   Table[
    With[{gt = gm (1./eta - 1.) dt},
     With[{step = Function[rho, dephaseAxis[dephaseAxis[measureStepNC[Pi/2, dvar][rho],
          PauliMatrix[3], gt], sigmaPhi[Pi/2], gt]]},
      Mean[Table[BlockRandom[SeedRandom[Round[1000 eta] + k];
         rlen@Last@NestList[step, {{0.5, 0.}, {0., 0.5}}, nst]], {k, ntr}]]]],
    {eta, {1., 0.64, 0.36}}]]]]
```

As expected, the steady radius tracks $\sqrt{\eta}$: only the ideal detector holds the state on the surface, and below $\eta=1$ the isotropic walk settles onto a smaller sphere whose radius drops as the square root of the efficiency. The lost information is the gap between monitoring an actual wavefunction and monitoring only its blurred image.

## 6. The Output Correlators and the Cosine at Zero Lag

Now the payoff. The state itself is hard to see, but the two records $I_z(t)$ and $I_\varphi(t)$ are right there, and their temporal correlators carry the physics in closed form. Korotkov and Atalaya derived the self- and cross-correlators of the two output signals directly from the Bayesian equations, and for equally strong ideal channels with no extra evolution they are sums of two decaying exponentials with rates $\Gamma(1\mp\cos\varphi)$. Encode the closed forms:

```wl
ClearAll[kzzCF, kzpCF];
kzzCF[phi_, gm_, tau_] := (1/2)(1 + Cos[phi]) Exp[-gm (1 - Cos[phi]) tau] +
   (1/2)(1 - Cos[phi]) Exp[-gm (1 + Cos[phi]) tau];
kzpCF[phi_, gm_, tau_] := (1/2)(Exp[-gm (1 - Cos[phi]) tau] - Exp[-gm (1 + Cos[phi]) tau]) +
   (Cos[phi]/2)(Exp[-gm (1 - Cos[phi]) tau] + Exp[-gm (1 + Cos[phi]) tau])
```

Read the two limits that these forms make vivid: the self-correlator starts at one and the cross-correlator starts at the cosine of the angle between the axes:

```wl
{FullSimplify[kzzCF[phi, gm, 0]], FullSimplify[kzpCF[phi, gm, 0]]}
```

The self-correlator at zero lag is one and the cross-correlator at zero lag is $\cos\varphi$: aligned detectors ($\varphi=0$) see perfectly correlated records, orthogonal detectors ($\varphi=\pi/2$) see uncorrelated ones, and the cosine interpolates. This is a purely dynamical readout of the angle between two incompatible measurements, available without ever reconstructing the state.

Now recover these closed forms from the records themselves. Run an ensemble of orthogonal-and-angled measurement records, form the ensemble-and-time-averaged two-time correlators, and overlay the closed forms. First define a run that emits both readings while it measures:

```wl
ClearAll[recordNC];
recordNC[phi_, dvar_, nst_, seed_] := BlockRandom[SeedRandom[seed];
   Module[{szA = PauliMatrix[3], spA = sigmaPhi[phi], rho = {{0.5, 0.}, {0., 0.5}}, ibz, ibp, mm},
    Transpose@Table[ibz = drawRead[rho, szA, dvar]; ibp = drawRead[rho, spA, dvar];
      mm = kraus[ibz/2, szA, dvar] . kraus[ibp, spA, dvar] . kraus[ibz/2, szA, dvar];
      rho = mm . rho . ConjugateTranspose[mm]; rho = rho/Tr[rho]; {ibz, ibp}, {nst}]]]
```

Each step appends the two readings while the same readings drive the Bayesian kick, so the records and the state share the very same noise, which is what makes the cross-correlator carry the back-action. Now simulate two angles and compare with the closed form (this is the heavy cell, an ensemble of records for each angle; it takes a couple of minutes):

```wl
With[{gm = 1., dt = 0.02, nst = 250, ntr = 1000, nlag = 10, phis = {Pi/4, Pi/2}},
 Module[{dvar = 1./(2 gm dt), taus, mc, cf},
  taus = Range[nlag] dt;
  mc[phi_] := Module[{recs = Table[recordNC[phi, dvar, nst, 700 + Round[100 phi] + k], {k, ntr}]},
    {Table[Mean[Flatten[Table[#[[1, t]] #[[1, t + j]], {t, nst - j}] & /@ recs]], {j, nlag}],
     Table[Mean[Flatten[Table[#[[2, t + j]] #[[1, t]], {t, nst - j}] & /@ recs]], {j, nlag}]}];
  Show[
   ListPlot[Flatten[Table[{Transpose[{taus, mc[phi][[1]]}], Transpose[{taus, mc[phi][[2]]}]}, {phi, phis}], 1],
    Frame -> True, GridLines -> Automatic, PlotRange -> {{0, nlag dt}, {-0.1, 1.05}},
    FrameLabel -> {"lag \[Tau] (units of 1/\[CapitalGamma])", "output correlator"},
    PlotLegends -> {"K_zz, \[Phi]=\[Pi]/4", "K_z\[Phi], \[Phi]=\[Pi]/4", "K_zz, \[Phi]=\[Pi]/2", "K_z\[Phi], \[Phi]=\[Pi]/2"},
    PlotLabel -> "simulated records recover the closed-form correlators", ImageSize -> 520],
   Plot[Evaluate@Flatten@Table[{kzzCF[phi, gm, tau], kzpCF[phi, gm, tau]}, {phi, phis}], {tau, 0, nlag dt},
    PlotStyle -> Directive[Gray, Dashed]]]]]
```

Each simulated correlator sits on its closed-form curve: the self-correlators decay from one, the $\varphi=\pi/4$ cross-correlator starts near $\cos(\pi/4)\approx 0.71$ and the $\varphi=\pi/2$ cross-correlator starts near zero, exactly the cosine law read off from the records. The state was never reconstructed; the angle between the two incompatible measurements is written directly in how their noisy outputs correlate.

## 7. Reading a Hidden Rotation, and the Bacon-Shor Connection

One last use of the correlators, and one last connection. If a small unwanted rotation $\tilde\Omega_R$ about $y$ is present, it breaks the symmetry between the two orderings and shows up only in the antisymmetrized cross-correlator $K_{z\varphi}(\tau)-K_{\varphi z}(\tau)$, which is otherwise zero. That difference is a sensitive probe of a rotation far too slow to see directly. Show that the difference is proportional to $\tilde\Omega_R$, vanishing when the rotation does:

```wl
With[{Gp = 1.6, Gm = 0.4, tau = 0.5, phi = Pi/2},
 Table[(2 om Sin[phi]/(Gp - Gm)) (Exp[-Gm tau] - Exp[-Gp tau]), {om, {0., 0.05}}]]
```

With no rotation the antisymmetrized cross-correlator is zero, and with a small rotation it is nonzero and grows linearly with $\tilde\Omega_R$: averaging many records turns a small, otherwise invisible frequency mismatch into a clearly measurable difference signal. In the experiment that grounds this theory, a residual Rabi frequency of a few kilohertz riding on a forty-megahertz drive was read out this way.

The whole toy has a direct afterlife. A four-qubit Bacon-Shor code measures four two-qubit gauge operators continuously, and because two pairs of them do not commute, the encoded state carries a spare degree of freedom, a gauge qubit, that is measured along two orthogonal axes at once, exactly as here: it diffuses on a great circle while the logical qubit rides untouched. The code's error syndromes are then nothing but the cross-correlators of Section 6, positive when the code is clean and flipping sign when an error moves the state to a subspace where the two watched axes have swapped orientation. The nine-qubit error-correcting version watches twelve gauge operators and reads its syndromes from triple correlators of three records at a time. The one-qubit cross-correlator we pulled out of simulated noise is, in that setting, a continuously monitored error check.

## Where This Leaves Us (and What Comes Next)

You now have a complete, computation-first toolkit for one qubit under two simultaneous non-commuting detectors: a Gaussian Kraus operator built from any measured axis, a symbolic proof that two such operators fail to commute by a rotation about the third axis, a general update that carries the full matrix where the single-axis closed form cannot, the transition from collapse through localized diffusion to an isotropic random walk as the axes separate, an uncertainty floor that forbids the walk from ever stopping, a purity ceiling set by efficiency, and two output correlators whose zero-lag cross-value reads the cosine of the angle between the watched axes. Before moving on, the points that were computed and verified along the way:

- Measuring one axis is conjugation by a diagonal Gaussian Kraus operator, so the single-axis update stays a closed-form Bayes rule on the populations with the coherence riding along.
- Two axes' Kraus operators fail to commute by $2i\sinh(a/2D)\sinh(b/2D)\,\sigma_y$, inheriting $[\sigma_z,\sigma_x]=2i\sigma_y$; there is no basis diagonalizing both, hence no closed-form population update.
- Orthogonal detectors keep a pure state pure while forcing it to diffuse on a great circle, replacing collapse with a persistent walk that never lands.
- The angle $\varphi$ tunes the motion from collapse ($\varphi=0$) through localized diffusion to an isotropic walk ($\varphi=\pi/2$), with the Bloch length untouched throughout.
- The measurement disturbance is bounded below by the commutator magnitude, saturating only along the unreadable $y$-direction, so the diffusion is a consequence of the uncertainty principle; efficiency $\eta<1$ caps the steady radius at $\sqrt\eta$.
- The output self- and cross-correlators are sums of exponentials with rates $\Gamma(1\mp\cos\varphi)$; at zero lag $K_{zz}=1$ and $K_{z\varphi}=\cos\varphi$, recoverable directly from simulated records.
- The antisymmetrized cross-correlator reads a hidden slow rotation; the same cross-correlators are the error syndromes of a Bacon-Shor code whose gauge qubit is exactly this doubly-watched qubit.

The natural continuations are the pieces this essay set aside: a phase back-action term that tilts the diffusion off the great circle, the finite detector bandwidth that couples the qubit to a cavity field, real-time feedback that locks the walk, and, in place of one qubit's two axes, two qubits watched jointly, whose measurement-induced entanglement is the subject of the sibling essay `Watching-Two-Qubits.md`.
