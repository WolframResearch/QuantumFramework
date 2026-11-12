#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

# What the Wolfram Language and QuantumFramework Do for Physicists

#| style: Subtitle

Sixteen topics, each with a computation you can run live

#| style: Subsubtitle

Mohammad Bahrami, Wolfram Research

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 1. Continuous Quantum Error Correction: Syndrome Estimation and Feedback

The three-qubit bit-flip code protected in continuous time: two stabilizers measured weakly, a classical controller turning the noisy records into an estimate, and the estimate driving corrective pulses.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **The syndrome is readable even though each sample is useless.** The stabilizer signal sits well below the shot noise per step, yet four syndrome probabilities track every flip, through an exact Bayes update checked against the full stochastic master equation.
- **Feedback breaks the four-number description.** With $H_{\rm fb}=\lambda\sum_i p_iX_i$ on, populations couple to coherences, so the controller must carry the full $8\times8$ estimate.
- **The controller does not need the logical state.** Correction acts on syndromes, not on the logical amplitudes.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **A wrong detector-efficiency model costs protection,** and paired runs on identical noise separate the estimation gap from the physical loss.
- **Staying in the code space is not keeping the qubit.** The residual loss splits exactly into leakage plus population in $|1_L\rangle$, which no stabilizer can see.
- **Feedback that protects the code space exactly.** With $F_i=i\sqrt\gamma\,[\Pi_{\rm code},X_i]$, every encoded state is a fixed point at unit efficiency; below it, a $\frac{1-\eta}{\eta}\mathcal D[F_i]$ term appears.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 1.1: tracking bit flips from two noisy records

Pauli matrices, the stabilizers $Z_1Z_2$ and $Z_2Z_3$, and the four sector projectors.

```wl
{id2, X, Y, Z} = Table[PauliMatrix[j], {j, 0, 3}];
op3[o_, i_] := KroneckerProduct @@ ReplacePart[{id2, id2, id2}, i -> o];
ket[bits_] := Normal@SparseArray[FromDigits[bits, 2] + 1 -> 1, 8];
S1 = op3[Z, 1] . op3[Z, 2]; S2 = op3[Z, 2] . op3[Z, 3];
signs = {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
proj = (IdentityMatrix[8] + #[[1]] S1) . (IdentityMatrix[8] + #[[2]] S2)/4 & /@ signs;
```

Everything below is built from plain $8\times8$ matrices.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 1.1 The measurement step and the trajectory

The positivity-preserving step of Rouchon and Ralph (2015), plus a trajectory generator that returns both the states and the record.

```wl
measurementStep[H_, watched_, effs_, unwatched_, dt_] :=
  With[{id = IdentityMatrix[Length[H]],
    drift = I H + Total[ConjugateTranspose[#] . # & /@ Join[watched, unwatched]]/2,
    corr = Total[MapThread[#1 (#2 . #2) &, {effs, watched}]]/2},
   Function[{rho, dJ}, With[{sig = Total[MapThread[Sqrt[#2] #1 #3 &, {dJ, effs, watched}]]},
     With[{M = id - drift dt + sig + sig . sig/2 - corr dt},
      With[{top = M . rho . ConjugateTranspose[M] +
          dt Total[# . rho . ConjugateTranspose[#] & /@ unwatched] +
          dt Total[MapThread[(1 - #2) #1 . rho . ConjugateTranspose[#1] &, {watched, effs}]]},
       top/Re@Tr[top]]]]]];
measurementRecord[rho_, watched_, effs_, dt_, kick_] :=
  MapThread[Sqrt[#3] Re@Tr[(#1 + ConjugateTranspose[#1]) . rho] dt + #2 &, {watched, kick, effs}];
trajectory[rho0_, H_, watched_, effs_, unwatched_, dt_, tf_, seed_] :=
  BlockRandom[SeedRandom[seed];
   With[{step = measurementStep[H, watched, effs, unwatched, dt],
     kicks = RandomVariate[NormalDistribution[0, Sqrt[dt]], {Round[tf/dt], Length[watched]}]},
    With[{states = FoldList[step[#1, measurementRecord[#1, watched, effs, dt, #2]] &, rho0, kicks]},
     <|"states" -> states,
      "record" -> MapThread[measurementRecord[#1, watched, effs, dt, #2] &, {Most[states], kicks}]|>]]];
```

Each step conditions the state on one record increment without ever leaving the set of valid states.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 1.1 One monitored run

Start in $|000\rangle$, let each qubit flip at rate $\gamma$, watch both stabilizers, and keep only the record.

```wl
kappa = 1.; eta = 1.; gamma = 0.05; dt = 0.004;
run = trajectory[KroneckerProduct[#, #] &@ket[{0, 0, 0}], 0 IdentityMatrix[8],
   Sqrt[kappa] {S1, S2}, {eta, eta}, Sqrt[gamma] (op3[X, #] & /@ Range[3]), dt, 40., 3];
```

The run holds the true conditional state and the noisy record it produced.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 1.1 The four-number filter

Relax each sector toward $1/4$ under the flips, multiply by the likelihood of the new pair of record increments, renormalize.

```wl
syndromeFilter[p0_, record_] := With[{c = 2 Sqrt[eta kappa]},
   FoldList[With[{w = (1/4 + Exp[-4 gamma dt] (#1 - 1/4)) Exp[c signs . #2]}, w/Total[w]] &,
    p0, record]];
pFilter = syndromeFilter[{1., 0., 0., 0.}, run["record"]];
```

The filter sees nothing but the record.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 1.1 Syndrome probabilities in time

The filter's four probabilities over the whole run.

```wl
ListLinePlot[Transpose[pFilter], DataRange -> {0, 40}, PlotRange -> {-0.05, 1.05},
 PlotLegends -> {"code", "X1 flipped", "X2 flipped", "X3 flipped"}, Frame -> True,
 FrameLabel -> {"time", "syndrome probability"}]
```

Certainty passes from sector to sector as flips accumulate.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 1.1 Four numbers against eight dimensions

Compare the filter with the sector populations of the full conditional state.

```wl
Max@Abs[pFilter - Map[Function[r, Re@Tr[# . r] & /@ proj], run["states"]]]
```

They agree up to the finite-step difference between two discretizations of one continuous filter.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 1.2: exact protection by Markovian feedback

The averaged dynamics needs only the Liouvillian.

```wl
liouvillian[H_, leaks_List, d_] := With[{id = IdentityMatrix[d]},
   -I (KroneckerProduct[H, id] - KroneckerProduct[id, Transpose[H]]) +
    Total[Function[c, KroneckerProduct[c, Conjugate[c]] -
        (KroneckerProduct[ConjugateTranspose[c] . c, id] +
           KroneckerProduct[id, Transpose[ConjugateTranspose[c] . c]])/2] /@ leaks]];
evolve[H_, leaks_, rho0_, t_] := With[{d = Length[rho0]},
   ArrayReshape[MatrixExp[liouvillian[H, leaks, d] t] . Flatten[rho0], {d, d}]];
```

`evolve` exponentiates it acting on the vectorized density matrix.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 1.2 The code space is invariant, for every rate

Keep the error rate symbolic, build the feedback operators from the code projector, and check the two identities.

```wl
codeProj = First[proj];
cs = Sqrt[\[Gamma]] (op3[X, #] & /@ Range[3]);
F = I (codeProj . # - # . codeProj) & /@ cs;
{Union@Flatten@Simplify@MapThread[#1 . #2 + #2 . #1 &, {cs, F}],
 Union@Flatten@Simplify[(# . codeProj) & /@ (cs - I F)]}
```

Both lists are `{0}`: $c_iF_i+F_ic_i=0$ and $(c_i-iF_i)\Pi_{\rm code}=0$ for every $\gamma$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 1.2 The price of an imperfect detector

Code-space population at $t=15$ with and without feedback, across detector efficiencies.

```wl
codePop[e_] := With[{c = cs /. \[Gamma] -> gamma, f = F /. \[Gamma] -> gamma},
   Re@Tr[codeProj . evolve[0 IdentityMatrix[8], Join[c - I f, Sqrt[(1 - e)/e] f],
       KroneckerProduct[#, #] &@ket[{0, 0, 0}], 15.]]];
bare = Re@Tr[codeProj . evolve[0 IdentityMatrix[8], cs /. \[Gamma] -> gamma,
     KroneckerProduct[#, #] &@ket[{0, 0, 0}], 15.]];
etas = {0.01, 0.03, 0.1, 0.3, 0.6, 1.};
ListLogLinearPlot[{Table[{e, codePop[e]}, {e, etas}], {#, bare} & /@ etas}, Joined -> True,
 PlotMarkers -> Automatic, PlotRange -> {0, 1.05}, Frame -> True,
 PlotLegends -> {"Markovian feedback", "no feedback"},
 FrameLabel -> {"detector efficiency \[Eta]", "code-space population at t = 15"}]
```

At $\eta=1$ the code space is held exactly; at low efficiency the loop feeds amplified detector noise back onto the qubits and does worse than no loop at all.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 2. A Computational Guide to Open Quantum Systems

Twenty open-system problems built from about ten shared functions on plain matrices, each stating its equation, computing, and checking itself against a closed form or a second method.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **Decoherence is the price of discarding the record.** Conditioned trajectories stay pure while the ensemble average mixes.
- **The naive stochastic master equation fails.** One large noise kick gives a negative eigenvalue; the Kraus-form step lands on the boundary of state space instead.
- **Cats.** Coherence dies at a rate set by the square of the separation, checked against the exact Kraus solution of photon loss.
- **Quantum Brownian motion.** High-temperature Caldeira-Leggett is not of Lindblad form, so the uncertainty bound is watched directly.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **One leak, three unravellings.** Counting, homodyne and heterodyne give different histories with the same average.
- **Two thresholds where one is expected.** A watched, driven qubit has an exceptional point at $k=\Omega$ and a Zeno freeze rate $\Omega^2/2k$.
- **Born's rule as a martingale.** Every run tends to a pole while the mean stays fixed, which forces $p_\pm=(1\pm z_0)/2$.
- **Quantum Kalman filter.** The steady conditional state sits at the Heisenberg floor for every measurement strength.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **Feedback in three guises.** Wiseman-Milburn purification onto a dark state, Jacobs rapid purification, and feedback cooling with a floor no gain can cross.
- **Thermalization derives the canonical ensemble.** Detailed balance forces the whole thermal distribution, not only its mean.
- **Mollow triplet and squeezed fluorescence,** from one Liouvillian by quantum regression.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 2.1: the quantum Zeno crossover as an exceptional point

Pauli matrices and the Liouvillian of a Lindblad equation.

```wl
{id2, X, Y, Z} = Table[PauliMatrix[j], {j, 0, 3}];
liouvillian[H_, leaks_List, d_] := With[{id = IdentityMatrix[d]},
   -I (KroneckerProduct[H, id] - KroneckerProduct[id, Transpose[H]]) +
    Total[Function[c, KroneckerProduct[c, Conjugate[c]] -
        (KroneckerProduct[ConjugateTranspose[c] . c, id] +
           KroneckerProduct[id, Transpose[ConjugateTranspose[c] . c]])/2] /@ leaks]];
```

One function covers every example in this topic.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 2.1 A watched, driven qubit, fully symbolic

A qubit Rabi-driven at $\Omega$ while $\hat\sigma_z$ is watched at strength $k$.

```wl
Eigenvalues@FullSimplify[liouvillian[(\[CapitalOmega]/2) X, {Sqrt[k] Z}, 2], \[CapitalOmega] > 0 && k > 0]
```

The Rabi pair $-k\pm\sqrt{k^2-\Omega^2}$ oscillates for $k<\Omega$, coalesces at $k=\Omega$, and splits into two real rates beyond it.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 2.1 The Zeno freeze

The slow rate at strong measurement.

```wl
Normal@Series[-k + Sqrt[k^2 - \[CapitalOmega]^2], {k, Infinity, 1}]
```

Watching harder slows the qubit as $\Omega^2/2k$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 2.1 The crossover in one picture

Rabi frequency and slowest decay rate against measurement strength.

```wl
Plot[{Re@Sqrt[1 - r^2], r - Re@Sqrt[r^2 - 1]}, {r, 0, 3}, Frame -> True, GridLines -> {{1}, None},
 PlotLegends -> {"Rabi frequency / \[CapitalOmega]", "slowest decay rate / \[CapitalOmega]"},
 FrameLabel -> {"k/\[CapitalOmega]", "Liouvillian eigenvalue / \[CapitalOmega]"}]
```

Watching harder speeds relaxation only up to $k=\Omega$; past it the qubit freezes.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 2.2: the Mollow triplet and squeezed fluorescence

Two-time correlations run on the same Liouvillian as the state, so a spectrum is one linear solve per frequency; at $\omega=0$ the trace condition picks the physical branch.

```wl
lower = (X - I Y)/2;
steady[H_, leaks_] := With[{v = ArrayReshape[First@NullSpace[liouvillian[H, leaks, 2]], {2, 2}]}, v/Tr[v]];
resolventSpectrum[big_, seedMat_, readout_, tones_] := With[{n = Length[big], d = Length[readout]},
   Table[2 Re@Tr[readout . ArrayReshape[If[Abs[w] < 10^-9,
       LeastSquares[Join[N[-big], {1. Flatten[IdentityMatrix[d]]}], Join[Flatten[seedMat], {0.}]],
       LinearSolve[I w IdentityMatrix[n] - big, Flatten[seedMat]]], {d, d}]], {w, tones}]];
```

Steady state and spectrum come from the same generator.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 2.2 Weak against strong drive

The inelastic fluorescence spectrum of a resonantly driven atom.

```wl
mollow[om_, tones_] := With[{big = liouvillian[(om/2) X, {lower}, 2], rss = steady[(om/2) X, {lower}]},
   With[{dm = lower - Tr[lower . rss] id2}, resolventSpectrum[big, dm . rss, ConjugateTranspose[dm], tones]]];
tones = Range[-9, 9, 0.05];
{weak, strong} = mollow[#, tones] & /@ {0.6, 6.};
ListLinePlot[{Transpose[{tones, weak/Max[weak]}], Transpose[{tones, strong/Max[strong]}]},
 Frame -> True, PlotRange -> All, PlotLegends -> {"\[CapitalOmega] = 0.6\[Gamma]", "\[CapitalOmega] = 6\[Gamma]"},
 FrameLabel -> {"frequency - drive frequency", "inelastic spectrum (scaled to peak)"}]
```

Strong driving dresses the atom and splits one line into three, with sidebands at $\pm\Omega$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 2.2 Homodyne noise at zero frequency

Compare the numerical homodyne spectrum with the closed form $S(0)=1+8\Omega^2(4\Omega^2-1)/(1+2\Omega^2)^3$.

```wl
homodyneS0[om_] := Module[{big = liouvillian[(om/2) X, {lower}, 2], rss = steady[(om/2) X, {lower}], cth, mq},
   cth = -I lower; mq = cth + ConjugateTranspose[cth];
   1 + First@resolventSpectrum[big, cth . rss + rss . ConjugateTranspose[cth] - Re@Tr[mq . rss] rss,
      mq - Re@Tr[mq . rss] id2, {0.}]];
Max@Table[Abs[homodyneS0[om] - (1 + 8 om^2 (4 om^2 - 1)/(1 + 2 om^2)^3)], {om, 0.05, 1.2, 0.05}]
```

The two agree to machine precision.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 2.2 The deepest squeezing, exactly

Minimize the closed form over the drive.

```wl
Minimize[{1 + 8 \[CapitalOmega]^2 (4 \[CapitalOmega]^2 - 1)/(1 + 2 \[CapitalOmega]^2)^3, \[CapitalOmega] > 0}, \[CapitalOmega]]
```

Fluorescence noise falls below the vacuum level for every $\Omega<\gamma/2$: nonclassical light from a single atom, visible only to phase-sensitive detection.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 3. Wong-Zakai Regularization for Continuous Measurement with NDSolve

A monitored, driven qubit and its readout as one stochastic equation, converted to Stratonovich form and integrated by `NDSolve` on a smooth interpolant of the noise.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **The exact algebra is checked symbolically, with all six rates symbolic.** `ItoProcess[StratonovichProcess[...]]` lands exactly on the Lindblad drift.
- **A positivity condition derived, not assumed.** The trajectory stays a density matrix if and only if $\Gamma_{CI}+\Gamma_{BA}\le2(\Gamma_d+\gamma_\phi)$.
- **The record is the noise.** Subtracting the integrated signal from the record returns the Wiener path that kicked the state.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **A closed law tying state to record.** With drive and relaxation off, $z(t)=\tanh(\operatorname{arctanh}z_0+\sqrt{\Gamma_{CI}}\,Q(t))$ on every path.
- **Honest numerics.** Means against the closed-form Lindblad solution, spreads against the built-in order-3/2 integrator, and the noise grid shown to be the step that matters.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 3.1: Itô-Stratonovich, Lindblad and positivity, all symbolic

The Stratonovich drift, the diffusion and the Itô drift on a symbolic rate vector $\{\Gamma_{CI},\Gamma_{BA},\Gamma_d,\gamma_\phi,\gamma_1,\Omega_x\}$.

```wl
driftStrat[{x_, y_, z_, q_}, {gCI_, gBA_, gd_, gphi_, g1_, om_}] := With[{g2 = g1/2 + gphi + gd},
   {((gBA + gCI)/2 - g2 - gCI z^2) x - Sqrt[gBA gCI] y z,
    ((gBA + gCI)/2 - g2 - gCI z^2) y + Sqrt[gBA gCI] x z - om z,
    g1 (1 - z) + gCI z (1 - z^2) + om y,
    Sqrt[gCI] z}];
diffStrat[{x_, y_, z_, q_}, {gCI_, gBA_, __}] :=
  {-Sqrt[gCI] x z - Sqrt[gBA] y, -Sqrt[gCI] y z + Sqrt[gBA] x, Sqrt[gCI] (1 - z^2), 1};
driftIto[{x_, y_, z_}, {_, _, gd_, gphi_, g1_, om_}] := With[{g2 = g1/2 + gphi + gd},
   {-g2 x, -g2 y - om z, g1 (1 - z) + om y}];
rs = {Subscript[\[CapitalGamma], "CI"], Subscript[\[CapitalGamma], "BA"], Subscript[\[CapitalGamma], "d"],
   Subscript[\[Gamma], "\[Phi]"], Subscript[\[Gamma], "1"], Subscript[\[CapitalOmega], "x"]};
```

Four components: three Bloch coordinates and the detector's accumulated readout $Q$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 3.1 Back to Itô form

Convert the Stratonovich process to Itô form and compare with the Lindblad drift plus the readout row.

```wl
Assuming[{rs[[1]] >= 0, rs[[2]] >= 0}, FullSimplify[
  (ItoProcess[StratonovichProcess[{driftStrat[{x, y, z, q}, rs], List /@ diffStrat[{x, y, z, q}, rs], {x, y, z, q}},
       {{x, y, z, q}, {x0, y0, z0, q0}}, t]]["Drift"] /. {x[t] -> x, y[t] -> y, z[t] -> z, q[t] -> q})
   == Append[driftIto[{x, y, z}, rs], Sqrt[rs[[1]]] z]]]
```

`True`, for all rates.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 3.1 The drift is the master equation

Build the Lindbladian with drive, relaxation and dephasing, project on the Pauli matrices, subtract.

```wl
With[{sx = PauliMatrix[1], sy = PauliMatrix[2], sz = PauliMatrix[3], sm = {{0, 1}, {0, 0}}},
 With[{rho = (IdentityMatrix[2] + x sx + y sy + z sz)/2},
  With[{diss = Function[L, L . rho . ConjugateTranspose[L] -
       (ConjugateTranspose[L] . L . rho + rho . ConjugateTranspose[L] . L)/2]},
   Simplify[(Tr[(-I (rs[[6]]/2) (sx . rho - rho . sx) + rs[[5]] diss[sm] +
          ((rs[[4]] + rs[[3]])/2) diss[sz]) . #] & /@ {sx, sy, sz}) - driftIto[{x, y, z}, rs]]]]]
```

`{0, 0, 0}`: the affine drift is the master equation.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 3.1 The flux of $r^2$ through the Bloch sphere

The quantity whose sign decides whether a trajectory can leave the set of density matrices.

```wl
With[{b = Most@diffStrat[{x, y, z, q}, rs]},
 FullSimplify[(2 {x, y, z} . driftIto[{x, y, z}, rs] + b . b) /.
   {x -> Sqrt[1 - z^2] Cos[\[Phi]], y -> Sqrt[1 - z^2] Sin[\[Phi]]}]]
```

The flux points inward everywhere exactly when $\Gamma_{CI}+\Gamma_{BA}\le2(\Gamma_d+\gamma_\phi)$: a detector cannot extract more information than the dephasing it causes pays for.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 3.2: NDSolve integrating a monitored qubit

Interpolate a sampled Wiener path with a cubic and give its derivative to `NDSolve` as a classical forcing.

```wl
wienerFun[tf_, h_, seed_] := BlockRandom[SeedRandom[seed];
   With[{grid = Subdivide[0., tf, Ceiling[tf/h]]},
    Interpolation[Transpose[{grid, Prepend[Accumulate[
         RandomVariate[NormalDistribution[0, Sqrt[grid[[2]]]], Length[grid] - 1]], 0.]}],
     InterpolationOrder -> 3]]];
oneTraj[r_, init_, tf_, h_, seed_] := Module[{Wf = wienerFun[tf, h, seed], v = {xx, yy, zz, QQ}, t},
   NDSolveValue[Join[
     Thread[Through[(Derivative[1] /@ v)[t]] ==
       driftStrat[Through[v[t]], r] + diffStrat[Through[v[t]], r] Wf'[t]],
     Thread[Through[v[0]] == init]], v, {t, 0, tf}, MaxStepSize -> h/4]];
```

The stochastic equation becomes an ordinary one, path by path.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 3.2 One trajectory and its readout

A driven, relaxing, quantum-limited readout starting in the ground state.

```wl
tr = oneTraj[{0.252, 0.1, 0.176, 1., 1., 3.}, {0., 0., 1., 0.}, 3., 1/120., 7];
Plot[{tr[[3]][t], tr[[4]][t]}, {t, 0, 3}, Frame -> True, PlotLegends -> {"z(t)", "Q(t)"}]
```

The monitored component and the readout move together because one noise drives both.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 3.2 State as a function of the record

Drive and relaxation off: test $z=\tanh(\operatorname{arctanh}z_0+\sqrt{\Gamma_{CI}}\,Q)$ along the path.

```wl
With[{rG = {0.352, 0., 0.176, 0., 0., 0.}},
 With[{trG = oneTraj[rG, {0., 0., 0.6, 0.}, 3., 1/(40 0.352), 11], s = Subdivide[0., 3., 1000]},
  Max@Abs[trG[[3]][s] - Tanh[ArcTanh[0.6] + Sqrt[rG[[1]]] trG[[4]][s]]]]]
```

The residual is at the solver's tolerance: the population is a function of what the detector has read so far.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 4. Continuous Quantum Measurement, from Bayes' Rule to the Quantum Kalman Filter

A particle whose position is weakly and continuously measured, in any potential: each step takes one noisy reading, updates the wavefunction by Bayes' rule, then evolves by a split-operator step.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **Measurement is literally Bayes' rule.** The amplitude update and "likelihood times prior" agree to machine precision.
- **Free particle: localization and heating.** The packet narrows to a steady width while the averaged momentum variance grows as $\lambda\hbar^2t$.
- **A cat collapses, and the readings alone tell which way.** A running log-odds from the outcomes tracks the hidden lobe mass.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **Squeezing below the zero-point at any strength.** The conditional covariances solve a Riccati equation in closed form, with $\lambda,\omega,m,\hbar$ symbolic.
- **The quantum Kalman filter.** Five numbers fed by the readings reproduce the conditional mean of the full wavefunction.
- **The price of an imperfect detector, in closed form.** The steady purity is exactly $\sqrt\eta$, independent of trap, mass and strength.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 4.1: squeezing below the zero-point, in closed form

The Riccati flow for the conditional covariances in a harmonic trap, its physical fixed point, and its determinant.

```wl
riccati = {2 Cc/m - 4 \[Lambda] Sx^2, Sp/m - m \[Omega]^2 Sx - 4 \[Lambda] Sx Cc,
   -2 m \[Omega]^2 Cc + \[Lambda] \[HBar]^2 - 4 \[Lambda] Cc^2};
steady = SolveValues[Thread[riccati == 0], {Sx, Cc, Sp}] // FullSimplify;
positiveQ[sol_] := TrueQ@Resolve[ForAll[{m, \[Lambda], \[Omega], \[HBar]},
     m > 0 && \[Lambda] > 0 && \[Omega] >= 0 && \[HBar] > 0,
     (sol[[1]] | sol[[3]]) \[Element] Reals && sol[[1]] > 0 && sol[[3]] > 0], Reals];
{SxSS, CcSS, SpSS} = SelectFirst[steady, positiveQ];
as = {\[Lambda] > 0, \[Omega] > 0, m > 0, \[HBar] > 0};
Simplify[SxSS SpSS - CcSS^2, as]
```

$\hbar^2/4$: the steady conditional state is pure.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 4.1 One function of one strength

In ground-state units the position variance depends only on $\mu=2\hbar\lambda/m\omega^2$.

```wl
sqRatio[u_] := Sqrt[2 (Sqrt[1 + u^2] - 1)]/u;
{Simplify[SxSS/(\[HBar]/(2 m \[Omega])) == sqRatio[2 \[HBar] \[Lambda]/(m \[Omega]^2)], as],
 Resolve[ForAll[u, u > 0, sqRatio[u] < 1], Reals],
 Asymptotic[sqRatio[u], {u, 0, 2}, Assumptions -> u >= 0],
 Asymptotic[sqRatio[u], {u, Infinity, 1}, Assumptions -> u >= 0]}
```

Squeezed below the zero-point for every strength, weakly as $1-\mu^2/8$ and strongly as $\sqrt{2/\mu}$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 4.1 Detector efficiency

A detector catching a fraction $\eta$ rescales $\lambda\to\eta\lambda$ and $\hbar\to\hbar/\sqrt\eta$ in the information terms only.

```wl
{SxE, CcE, SpE} = {SxSS, CcSS, SpSS} /. {\[Lambda] -> \[Eta] \[Lambda], \[HBar] -> \[HBar]/Sqrt[\[Eta]]};
{Simplify[\[HBar]/(2 Sqrt[SxE SpE - CcE^2]), {\[HBar] > 0, 0 < \[Eta] <= 1, \[Lambda] > 0, \[Omega] > 0, m > 0}],
 Simplify@Reduce[sqRatio[Sqrt[\[Eta]] u]/Sqrt[\[Eta]] <= 1 && u > 0 && 0 < \[Eta] < 1, u, Reals]}
```

The steady purity is exactly $\sqrt\eta$, and squeezing remains possible past a threshold $\mu^*=2\sqrt{1-\eta}/\eta^{3/2}$ that diverges as the detector fails.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 4.2: five numbers recover what 256 amplitudes compute

A periodic grid, a Gaussian packet, the Bayesian measurement step and a Strang split-operator step ($\hbar=m=1$).

```wl
hbar = 1; mass = 1;
grid[n_, l_] := With[{dx = l/n}, <|"n" -> n, "x" -> dx Range[-n/2, n/2 - 1],
    "p" -> (2 Pi hbar/l) RotateLeft[Range[-n/2, n/2 - 1], n/2]|>];
ft[psi_] := Fourier[psi, FourierParameters -> {0, -1}];
ift[c_] := InverseFourier[c, FourierParameters -> {0, -1}];
gaussian[g_, x0_, p0_, s_] := Normalize[Exp[-(g["x"] - x0)^2/(4 s^2) + I p0 g["x"]/hbar]];
xMean[g_, psi_] := g["x"] . Abs[psi]^2;
outcome[g_, lam_, dt_][psi_] := RandomChoice[Abs[psi]^2 -> g["x"]] +
   RandomVariate[NormalDistribution[0, 1/Sqrt[4 lam dt]]];
measUpdate[g_, lam_, dt_, xb_][psi_] := Normalize[Exp[-lam dt (g["x"] - xb)^2] psi];
potentialStep[g_, V_, dt_][psi_] := Exp[-I V[g["x"]] dt/(2 hbar)] psi;
kineticStep[g_, dt_][psi_] := ift[Exp[-I g["p"]^2 dt/(2 mass hbar)] ft[psi]];
evolveStep[g_, V_, dt_] := potentialStep[g, V, dt] @* kineticStep[g, dt] @* potentialStep[g, V, dt];
trajectory[g_, V_, lam_, dt_, nt_][psi0_] := NestList[
   evolveStep[g, V, dt][measUpdate[g, lam, dt, Sow@outcome[g, lam, dt][#]][#]] &, psi0, nt];
```

The measurement half never sees the Hamiltonian, so any $V(x)$ costs one pointwise phase.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 4.2 A watched particle in a trap

Keep both the conditional wavefunction and the stream of readings.

```wl
wR = 1.; lamH = 2.; SxZP = hbar/(2 mass wR); SpZP = hbar mass wR/2;
gH = grid[256, 24.]; dtH = 0.002; ntH = 1500;
{states, {outs}} = BlockRandom[SeedRandom[21];
   Reap[trajectory[gH, 0.5 mass wR^2 #^2 &, lamH, dtH, ntH][gaussian[gH, 0., 0., Sqrt[SxZP]]]]];
xGrid = xMean[gH, #] & /@ states;
```

The readings are the only thing an experimenter would have.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 4.2 The quantum Kalman filter

Two means nudged by the surprise in each reading, three covariances following the Riccati flow of Example 4.1.

```wl
num = {\[HBar] -> hbar, m -> mass, \[Omega] -> wR, \[Lambda] -> lamH};
covDrift = Function[{sx, cc, sp}, Evaluate[riccati /. num /. {Sx -> sx, Cc -> cc, Sp -> sp}]];
filterStep[dt_][{mx_, mp_, sx_, cc_, sp_}, dY_] := With[{dWf = dY - 2 Sqrt[lamH] mx dt},
   Join[{mx + mp/mass dt + 2 Sqrt[lamH] sx dWf, mp - mass wR^2 mx dt + 2 Sqrt[lamH] cc dWf},
    {sx, cc, sp} + dt covDrift[sx, cc, sp]]];
mFilter = FoldList[filterStep[dtH], {0., 0., SxZP, 0., SpZP}, 2 Sqrt[lamH] dtH outs][[All, 1]];
ListLinePlot[{xGrid, mFilter}, DataRange -> {0, ntH dtH},
 PlotStyle -> {Directive[Thick, Gray], Directive[Red, Dashed]},
 PlotLegends -> {"grid (256 amplitudes)", "Kalman filter (5 numbers)"}, Frame -> True,
 FrameLabel -> {"t", "\[LeftAngleBracket]x\[RightAngleBracket]"}]
```

The filter, fed nothing but the readings, sits on the conditional mean of the full wavefunction.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 4.2 How close

The largest gap, relative to the range the mean covers.

```wl
Max@Abs[mFilter - xGrid]/(Max@xGrid - Min@xGrid)
```

In a quadratic potential the conditional state is a Gaussian whose center follows the record, so the wavefunction is never needed.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 5. Transmon cQED: Wolfram and Qolab

The demo shown at APS 2026 with Qolab: from the transmon circuit Hamiltonian to a prepared qubit state, down to single-shot (I, Q) points at the digitizer, and the same analysis on real hardware data.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **From circuit to qubit.** The Cooper-pair-box Hamiltonian $4E_C(\hat n-n_g)^2-E_J\cos\hat\phi$, and the spectrum flattening against offset charge as $E_J/E_C$ grows.
- **Any state with two hardware operations.** A virtual Z followed by one calibrated pulse reaches any $(\theta,\varphi)$.
- **Measurement is a cavity, not a projector.** The readout resonator rings up along two qubit-dependent spirals whose separation sets the signal-to-noise.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **The whole chain, both directions.** Simulated single-shot IQ records, then Gaussian-mixture fitting, assignment errors and maximum-likelihood tomography.
- **Real hardware in the same pipeline.** `QolabExperiment` returns IQ data from a Qolab qubit, analyzed by the identical code.
- **Decoherence and DRAG,** with leakage to $|2\rangle$ in a three-level simulation.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 5.1: arbitrary state preparation from the pulse area

Load TransmonLab, define the device, a Gaussian drive pulse and a flat-top readout pulse.

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/mbahram/TransmonLab/TransmonLab/"];
Needs["TransmonLab`"];
sys = TransmonSystem["QubitFrequency" -> 5.10*10^9, "ResonatorFrequency" -> 6.00*10^9,
   "Anharmonicity" -> -220*10^6, "Coupling" -> 100*10^6, "CavityLinewidth" -> 1*10^6,
   "T1" -> 50*10^-6, "T2" -> 30*10^-6, "OmegaMax" -> 2 Pi*42*10^6, "EpsilonMax" -> 2 Pi*5*10^6, "NMax" -> 5];
gaussPulse = Pulse["Gaussian", 50]; readoutPulse = Pulse["FlatTop", 500, "RiseTime" -> 10];
```

The device is a list of named physical parameters.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 5.1 Calibration from the pulse area

On resonance only the pulse area matters: rotation angle per unit amplitude is $\Omega_{\max}\int{\rm env}\,dt$.

```wl
thetaPerA = sys["OmegaMax"] gaussPulse["EnvelopeArea"] 10^-9;
```

One number calibrates every rotation.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 5.1 Prepare $\theta=\pi/3$, $\varphi=\pi/6$

One frame rotation, one scaled pulse, then readout; compare the simulated Bloch vector with the target.

```wl
{thT, phT} = {Pi/3, Pi/6}; aT = N[thT/thetaPerA];
seq = PulseSequence[FrameRotation["Qubit", phT + Pi/2], Drive[Amplitude[aT]@gaussPulse, "Qubit"],
   Sync["Qubit", "Resonator"], Readout[readoutPulse, "Resonator", 500]];
res = TransmonPhysics[sys, seq, "Model" -> "QubitOnly"];
bz = #[gaussPulse["Duration"] 10^-9] & /@ Values[res["BlochTrajectories"][[1]]];
{bz, N@{Sin[thT] Cos[phT], Sin[thT] Sin[phT], Cos[thT]}}
```

The simulated state lands on the analytic target to the accuracy the truncated Gaussian envelope allows.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 5.1 The same constant as a Power Rabi sweep

Sweep the drive amplitude.

```wl
TransmonPhysics[sys, PulseSequence[Drive[gaussPulse, "Qubit"], Sync["Qubit", "Resonator"],
    Readout[readoutPulse, "Resonator", 500]], Sweep["Amplitude", Range[0, 1, 0.02]],
  "Model" -> "QubitOnly"]["Plot"]
```

The calibration constant appears experimentally as the period of the Rabi oscillation.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 5.2: from the digitizer back to the Bloch vector

The cavity pointer states, one spiral per qubit state.

```wl
dr = DispersiveReadout[sys, readoutPulse];
IQPlot[dr]
```

Their separation is what the readout can distinguish.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 5.2 Single shots in three bases

Prepare the target rotated into each tomography basis and simulate 200 noisy shots per basis. `IQAnalysis` warns that the linear-inversion estimate is unphysical.

```wl
bv = BlochVector[{Cos[thT/2], Exp[I phT] Sin[thT/2]}];
{thX, phX} = N@Rest@ToSphericalCoordinates[RotationMatrix[-Pi/2, {0, 1, 0}] . bv];
{thY, phY} = N@Rest@ToSphericalCoordinates[RotationMatrix[Pi/2, {1, 0, 0}] . bv];
iqt = IQTrajectory[sys, dr, "Bases" -> {"x", "y", "z"}, "NShots" -> 200,
   "InitialStates" -> <|"z" -> {1.0, thT, phT}, "x" -> {1.0, thX, phX}, "y" -> {1.0, thY, phY}|>,
   "Seed" -> 123];
iqa = IQAnalysis[iqt["MeasurementData"]];
```

The analysis sees only clouds of (I, Q) points.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 5.2 Maximum-likelihood reconstruction

Reconstruct inside the Bloch ball and compare with the target.

```wl
target = N@{Sin[thT] Cos[phT], Sin[thT] Sin[phT], Cos[thT]};
mleB = Re[Tr[iqa["DensityMatrixMLE"] . #] & /@ (PauliMatrix /@ {1, 2, 3})];
{mleB, Norm[mleB - target], iqa["BlochVectorAround"]}
```

The reconstructed Bloch vector lands near the prepared one; swapping in the table from `QolabExperiment` runs the identical analysis on hardware.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 6. The Gravitational Soliton

Can a quantum particle be held together by its own gravity? In the Schrodinger-Newton equation the wavefunction moves in the potential sourced by its own mass density, and quantum spreading balances self-gravity in a stationary soliton.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **The scale in one line.** Minimizing $\hbar^2/m\ell^2-Gm^2/\ell$ gives $\ell_\star=\hbar^2/Gm^3$: astronomical for an electron, laboratory-sized for a heavy nanoparticle.
- **A self-consistent loop.** Newton's equation and a radial Schrodinger eigenproblem on a finite-element mesh alternate until the energy stops changing.
- **The box is removed, not ignored.** $\Phi(R)=0$ shifts the energy by $+1/R$; a fit over six boxes recovers the Moroz-Penrose-Tod energy.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **A check the solver never saw.** The virial identity $2T+W=0$ forces $T/E=-1/3$; the raw solution misses it by exactly the boundary shift.
- **Excited solitons.** The first three bound states have zero, one and two radial nodes, each in its own $-1/r$-tailed potential.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 6.1: the natural length scale

Balance kinetic against gravitational energy.

```wl
Solve[D[\[HBar]^2/(m l^2) - G m^2/l, l] == 0, l]
```

$\ell=\hbar^2/Gm^3$: the only length the problem has.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 6.1 One solver for each half of the loop

The radial Newton operator, the Schrodinger component, and solvers for each.

```wl
Needs["NDSolve`FEM`"];
newtonOp = Laplacian[\[CapitalPhi][r], {r, \[Theta], \[Phi]}, "Spherical"];
hamiltonianOp = SchrodingerPDEComponent[{u[r], {r}},
   <|"ReducedPlanckConstant" -> 1, "SchrodingerPotential" -> \[CapitalPhi][r]|>];
sortStates[s_] := SortBy[Transpose[s], Re] // Transpose;
solveNewton[uFunc_, mesh_] := NDSolveValue[{newtonOp == uFunc^2/r^2,
    DirichletCondition[\[CapitalPhi][r] == 0, r == mesh["Bounds"][[1, 2]]]}, \[CapitalPhi], {r} \[Element] mesh];
solveSchrodinger[PhiIn_, mesh_] := sortStates@NDEigensystem[
   {hamiltonianOp /. \[CapitalPhi] -> PhiIn, DirichletCondition[u[r] == 0, True]}, u[r], {r} \[Element] mesh, 3,
   Method -> {"Eigensystem" -> {"Arnoldi", "MaxIterations" -> 3000, "BasisSize" -> 50, "Shift" -> -10},
     "Interpolation" -> {"ExtrapolationHandler" -> {Automatic, "WarningMessage" -> False}}}];
```

Each half is an ordinary finite-element problem.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 6.1 Iterate to self-consistency, in six boxes

Start from a normalized Gaussian and iterate until the eigenvalue settles, then fit $E(R)=E_\infty+C/R$. This cell takes several minutes.

```wl
solveSN[Rval_, nEl_] := Module[{mesh, s = 3, uk, phik, ek, eOld = Infinity, iter = 0},
   mesh = ToGradedMesh[{Line[{{0}, {Rval}}], <|"Alignment" -> "Left", "ElementCount" -> nEl|>}];
   uk = 2^(7/4)/(s^(3/2) Pi^(1/4)) r Exp[-r^2/s^2];
   phik = solveNewton[uk, mesh]; {ek, uk} = {#[[1, 1]], #[[2, 1]]} &@solveSchrodinger[phik, mesh];
   While[Abs[ek - eOld] > 10^-8 && iter < 120, eOld = ek; phik = solveNewton[uk, mesh];
    {ek, uk} = {#[[1, 1]], #[[2, 1]]} &@solveSchrodinger[phik, mesh]; iter++];
   <|"R" -> Rval, "Energy" -> ek, "Iterations" -> iter, "u" -> Head[uk], "Phi" -> phik|>];
results = Table[solveSN[Rv, 8000], {Rv, {20, 40, 60, 100, 200, 300}}];
fit = FindFit[{#R, #Energy} & /@ results, e + c/rr, {e, c}, rr]
```

The fit recovers the literature ground-state energy and $C$ at the predicted value one: the whole finite-box error is the missing $-1/r$ tail.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 6.1 Every box on one line

Raw energies and energies with $1/R$ removed.

```wl
ListLinePlot[{{#R, #Energy} & /@ results, {#R, #Energy - 1/#R} & /@ results}, PlotMarkers -> Automatic,
 PlotLegends -> {"E(R)", "E(R) - 1/R"}, Frame -> True, FrameLabel -> {"R", "E"},
 Epilog -> {Red, Dashed, InfiniteLine[{0, e /. fit}, {1, 0}]}]
```

Subtracting the boundary shift collapses every box onto the infinite-box value.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 6.2: the virial identity as an independent check

Kinetic and potential energy of the largest-box ground state.

```wl
sol = Last[results]; u0 = sol["u"]; Phi0 = sol["Phi"]; e0 = sol["Energy"];
t = 1/2 NIntegrate[D[u0[r], r]^2, {r, 0, 300}];
v = NIntegrate[Phi0[r] u0[r]^2, {r, 0, 300}];
{t + v - e0, t/e0, t/(e0 - 1/300), {-3 t - e0, -1/300.}}
```

$T+V$ reproduces the eigenvalue; the raw $T/E$ misses $-1/3$ by exactly the boundary shift, and removing it restores the virial value.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 7. The Stabilizer Formalism by Computation

A state named by the Pauli operators that fix it: $n$ generators replace $2^n$ amplitudes, Clifford gates become bit rewrites of a tableau, and one object carries the reader from one qubit to a thousand, then to codes, magic and noise.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **Why Paulis, precisely.** Every state is fixed by some operator; four properties of the Pauli group make the description pay.
- **Schrodinger and Heisenberg meet.** $UZ_1U^\dagger$ stays one Pauli string, while $TXT^\dagger=(X+Y)/\sqrt2$ shows where the method breaks.
- **Stabilizer states are flat.** A uniform superposition over an affine subspace of $\mathbb F_2^n$.
- **Entanglement is a rank.** Entropy across a cut is an $\mathbb F_2$ rank, fast even for a thousand qubits.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **Measurement without branching.** Each random outcome is carried as a formal bit, so correlations come out as Boolean relations.
- **Codes.** The five-qubit code's fifteen single-qubit errors map one-to-one onto the fifteen nonzero syndromes.
- **Past the Clifford edge.** Stabilizer rank grows as $2^t$ for $t$ T gates, independent of $n$.
- **An honest benchmark** against Stim and QuantumClifford.jl, stating where QF is behind and what it offers instead.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 7.1: the perfect code's syndrome table

The five-qubit code, its four checks, and the syndrome of each single-qubit error.

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
code = PauliStabilizer["5QubitCode"];
checks = Most[code["Stabilizers"]];
syndrome[err_] := (1 - (code[err]["Expectation", #] & /@ checks))/2;
errors = Flatten[Table[{p -> q}, {p, {"X", "Y", "Z"}}, {q, 1, 5}], 1];
AssociationMap[syndrome, errors]
```

Each of the fifteen errors flips a different pattern of the four checks.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 7.1 What "perfect" means

Are the patterns exactly the fifteen nonzero four-bit strings?

```wl
Sort[syndrome /@ errors] === Sort@DeleteCases[Tuples[{0, 1}, 4], {0, 0, 0, 0}]
```

None missing and none repeated.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 7.1 A thousand-qubit GHZ state

Entanglement entropy across one qubit and across half the chain, without a state vector.

```wl
ghz = PauliStabilizer[1000][Prepend[Table["CNOT" -> {q - 1, q}, {q, 2, 1000}], "H"]];
{AbsoluteTiming[ghz["Entropy", {1}]], ghz["Entropy", Range[500]]}
```

One bit either way: GHZ entanglement is a single shared parity, however the chain is cut.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 7.2: measurement outcomes as algebra

Copy the parity of two random qubits onto a third, then measure all three symbolically.

```wl
parityCheck = PauliStabilizer[3][{"H" -> 1, "H" -> 2, "CNOT" -> {1, 3}, "CNOT" -> {2, 3}}];
smPar = parityCheck["SymbolicMeasure", {1, 2, 3}];
smPar["Stabilizers"]
```

Only two fresh outcome bits appear, because the conserved $ZZZ$ fixes the third.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 7.2 Reading outcomes by substitution

Substitute the free bits over all assignments, then draw a thousand samples by substitution alone.

```wl
readBits[ps_] := First@Keys@Select[AssociationThread[Tuples[{0, 1}, ps["Qubits"]],
      Normal[ps["State"]["StateVector"]]], # != 0 &];
free = DeleteDuplicates@Cases[smPar["Phase"], \[FormalS][_], Infinity];
{readBits[smPar["SubstituteOutcomes", Thread[free -> #]]] & /@ Tuples[{0, 1}, Length[free]],
 SeedRandom[7]; Counts[readBits /@ smPar["SampleOutcomes", 1000]]}
```

The third outcome is the XOR of the other two, and no sample leaves the even-parity sector.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 7.2 The price of magic

Three T gates on a GHZ backbone at growing width, and a single T on $|+\rangle$.

```wl
magicState[n_] := PauliStabilizer[n][Join[{"H" -> 1}, Table["CNOT" -> {i, i + 1}, {i, 1, n - 1}],
     {"T" -> 1, "T" -> 2, "T" -> 3}]];
{Table[{n, magicState[n]["Length"], 2^n}, {n, {4, 10, 20}}],
 FullSimplify[(1 + PauliStabilizer[1][{"H", "T"}]["Expectation", "X"])/2]}
```

The frame size is set by the number of T gates, not qubits, and the Born probability $(2+\sqrt2)/4$ comes straight from a two-term frame.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 8. Dirac Matrices in a Curved Spacetime

Gamma matrices satisfying $\{\gamma^\mu,\gamma^\nu\}=2g^{\mu\nu}\mathbb 1$ for any real symmetric metric, any signature and any dimension, built from the Brauer-Weyl construction and a vielbein $g=e^T\eta e$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **One function, every signature and dimension.** `GammaMatrices[g]` for Euclidean, Minkowski, split signature, Schwarzschild at a point, or a random non-diagonal metric.
- **The textbook identities as tests.** Traces, $\gamma^\mu\gamma_\mu=n$, $\gamma^\mu\gamma_\nu\gamma_\mu=(2-n)\gamma_\nu$ and chirality, each one `Table` over all signatures.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **Representations on demand.** Brauer-Weyl, Dirac and Weyl bases by explicit similarity transforms, odd dimensions included.
- **The graded Clifford basis.** All $2^n$ antisymmetrized products, the operator basis behind Fierz identities.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 8.1: gamma matrices for Schwarzschild

Load the DiracMatrix paclet.

```wl
PacletDirectoryLoad["/Users/mohammadb/Nextcloud/Mads/Quantum/DiracMatrix/DiracMatrixPackage"];
Needs["DiracMatrix`"];
```

Every function below comes from it.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 8.1 The equator at $r=3M$

The Schwarzschild metric there is $\mathrm{diag}(-1/3,3,9,9)$.

```wl
g = DiagonalMatrix[{-1/3, 3, 9, 9}];
Gs = GammaMatrices[g];
{Outer[#1 . #2 + #2 . #1 &, Gs, Gs, 1] == 2 TensorProduct[g, IdentityMatrix[4]], MatrixForm[Gs[[1]]]}
```

Its gamma matrices obey the Clifford relation exactly.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 8.1 A random curved metric and its vielbein

Draw a non-diagonal Lorentzian metric and rebuild it from its vielbein.

```wl
SeedRandom[42]; gr = RandomCurvedMetric[1, 3, 1, 1];
{sig, e} = MetricVielbein[gr];
{Diagonal[sig], NumericZeroQ[Transpose[e] . sig . e - gr]}
```

The signature is read off the metric itself, and the vielbein reproduces the metric to numerical tolerance.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 8.2: the identities on every signature

For $n=2,\ldots,6$ and every signature, a random non-diagonal metric: the Clifford relation, $\gamma^\mu\gamma_\mu=n$ and the four-gamma trace.

```wl
SeedRandom[12345];
Table[Module[{eta = RandomCurvedMetric[p, n - p, 0.5, 1], G, A},
   G = GammaMatrices[eta]; A = TensorProduct[eta, eta];
   {NumericZeroQ[Outer[#1 . #2 + #2 . #1 &, G, G, 1] - 2 TensorProduct[eta, IdentityMatrix[2^Floor[n/2]]]],
    NumericZeroQ[Total@MapThread[Dot, {G, Inverse[eta] . G}] - n IdentityMatrix[2^Floor[n/2]]],
    NumericZeroQ[Outer[Tr@*Dot, G, G, G, G, 1] -
      2^Floor[n/2] (A - Transpose[A, {1, 3, 2, 4}] + Transpose[A, {1, 4, 2, 3}])]}],
  {n, 2, 6}, {p, n - 1}] // Flatten // Apply[And]
```

`True`: the identities a QFT course derives by hand for Minkowski space hold for every metric in the sweep.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 8.2 Representations by one call

The Weyl basis in Minkowski space.

```wl
MatrixForm /@ ToWeylBasis[GammaMatrices[FlatMetric[1, 3]]]
```

The chiral blocks appear directly.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 9. Snyder's Lorentz-Covariant Quantum Spacetime

Snyder (1947) made spacetime coordinates noncommuting without breaking Lorentz symmetry, $[x^\mu,x^\nu]=ia^2M^{\mu\nu}$; the algebra is declared once and the noncommutative-algebra engine does what is error-prone by hand.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **Consistency first.** All Jacobi identities vanish, and Lorentz covariance of $[x,x]$ comes out as a theorem.
- **Orderings no one does by hand.** A quartic word in the coordinates reduces to the ordered basis automatically.
- **The deformed center.** $x^2-\tfrac{a^2}{2}M^{\mu\nu}M_{\mu\nu}$ commutes with everything, while $x^2$ alone does not.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **Noncommutativity is orbital angular momentum.** Jacobi closure forces $M^{\mu\nu}=x^\mu p^\nu-x^\nu p^\mu$.
- **Physics in one stroke.** Deformed momentum addition from the BCH series, Snyder against κ-Minkowski, and a minimal length $\Delta x_{\min}=a$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 9.1: Snyder's algebra

Ten generators, the three families of commutation relations, and the algebra.

```wl
eta[m_, n_] := DiagonalMatrix[{1, -1, -1, -1}][[m + 1, n + 1]];
MGen[m_, n_] := Which[m == n, 0, m < n, M[m, n], True, -M[n, m]];
gen = Join[Table[x[m], {m, 0, 3}], Flatten@Table[M[m, n], {m, 0, 3}, {n, m + 1, 3}]];
xxRel = Flatten@Table[Commutator[x[n], x[m]] + I a^2 M[m, n], {m, 0, 3}, {n, m + 1, 3}];
mxRel = Flatten@Table[Commutator[M[m, n], x[r]] - I (eta[n, r] x[m] - eta[m, r] x[n]),
    {m, 0, 3}, {n, m + 1, 3}, {r, 0, 3}];
mmRel = Replace[Subsets[Flatten[Table[{m, n}, {m, 0, 3}, {n, m + 1, 3}], 1], {2}],
   {{m_, n_}, {r_, s_}} :> Commutator[M[r, s], M[m, n]] +
     I (eta[n, r] MGen[m, s] - eta[n, s] MGen[m, r] - eta[m, r] MGen[n, s] + eta[m, s] MGen[n, r]), {1}];
alg = NonCommutativeAlgebra[<|"Generators" -> {gen},
    "CommutationRelations" -> Join[xxRel, mxRel, mmRel], "ScalarVariables" -> {a}|>];
```

The whole model is this one declaration.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 9.1 Jacobi identities

Three coordinates, and two coordinates with a Lorentz generator, over every index choice.

```wl
jac[X_, Y_, Z_] := Commutator[X, Commutator[Y, Z]] + Commutator[Y, Commutator[Z, X]] +
   Commutator[Z, Commutator[X, Y]];
{Union@Flatten@Table[NonCommutativeExpand[jac[x[m], x[n], x[r]], alg], {m, 0, 3}, {n, 0, 3}, {r, 0, 3}],
 Union@Flatten@Table[NonCommutativeExpand[jac[x[m], x[n], MGen[r, s]], alg],
   {m, 0, 3}, {n, 0, 3}, {r, 0, 3}, {s, 0, 3}]}
```

Every Jacobiator vanishes: the algebra is consistent.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 9.1 Lorentz covariance, derived

$[M,[x,x]]$ through the $[M,x]$ relations against $ia^2[M,M]$ through the Lorentz algebra.

```wl
Union@Flatten@Table[NonCommutativeExpand[Commutator[MGen[r, s], Commutator[x[m], x[n]]], alg] -
    I a^2 I (eta[s, m] MGen[r, n] - eta[s, n] MGen[r, m] - eta[r, m] MGen[s, n] + eta[r, n] MGen[s, m]),
  {r, 0, 3}, {s, r + 1, 3}, {m, 0, 3}, {n, m + 1, 3}]
```

The difference is zero, so covariance is a consequence, not an assumption.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 9.1 The Casimir

Commute $x^2-a^2C_1$ with every coordinate and generator, and $x^2$ alone with $x^0$.

```wl
C1 = Sum[eta[m, m] eta[n, n] MGen[m, n] ** MGen[m, n], {m, 0, 3}, {n, m + 1, 3}];
xSq = Sum[eta[m, m] x[m] ** x[m], {m, 0, 3}];
{Union@Table[NonCommutativeExpand[Commutator[xSq - a^2 C1, x[r]], alg], {r, 0, 3}],
 Union@Flatten@Table[NonCommutativeExpand[Commutator[xSq - a^2 C1, MGen[r, s]], alg], {r, 0, 3}, {s, r + 1, 3}],
 NonCommutativeExpand[Commutator[xSq, x[0]], alg]}
```

The Casimir is central; the nonzero last entry is exactly the residue the $-a^2C_1$ term cancels.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 9.2: deformed momentum addition

Expand $\log(e^{i\alpha\cdot x}e^{i\beta\cdot x})$ to fourth order in the BCH series.

```wl
as = Table[al[m], {m, 0, 3}]; bs = Table[be[m], {m, 0, 3}];
algE = NonCommutativeAlgebra[<|"Generators" -> {gen}, "CommutationRelations" -> Join[xxRel, mxRel, mmRel],
    "ScalarVariables" -> Join[{a}, as, bs]|>];
bch4[X_, Y_] := X + Y + Commutator[X, Y]/2 +
   (Commutator[X, Commutator[X, Y]] - Commutator[Y, Commutator[X, Y]])/12 -
   Commutator[Y, Commutator[X, Commutator[X, Y]]]/24;
ax = as . Table[x[m], {m, 0, 3}]; bx = bs . Table[x[m], {m, 0, 3}];
oplus = Simplify@Table[Coefficient[NonCommutativeExpand[bch4[I ax, I bx], algE], x[m]]/I, {m, 0, 3}];
{oplus[[1]], Simplify[oplus - (oplus /. {al -> be, be -> al})]}
```

Snyder's $\alpha\oplus\beta$ picks up $a^2$ corrections but stays symmetric under exchange.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 9.2 The same driver on κ-Minkowski

$[x^0,x^i]=ix^i/\kappa$.

```wl
algK = NonCommutativeAlgebra[<|"Generators" -> {Table[x[m], {m, 0, 3}]},
    "CommutationRelations" -> Join[Table[Commutator[x[i], x[0]] - I/k x[i], {i, 1, 3}],
      Flatten@Table[Commutator[x[j], x[i]], {i, 1, 3}, {j, i + 1, 3}]],
    "ScalarVariables" -> Join[{k}, as, bs]|>];
Simplify@Table[Coefficient[NonCommutativeExpand[bch4[I ax, I bx], algK], x[m]]/I, {m, 0, 1}]
```

Spatial momenta pick up a term that flips sign under exchange: two Lorentz-compatible deformations with different multi-particle physics.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 9.2 The minimal length

The deformed uncertainty relation $\Delta x\ge(1+a^2\Delta p^2)/2\Delta p$.

```wl
Assuming[a > 0, Minimize[{(1 + a^2 dp^2)/(2 dp), dp > 0}, dp] // FullSimplify]
```

$\Delta x_{\min}=a$: no event can be localized below the algebra's length.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 10. Breaking Elliptic-Curve Cryptography: Shor's Algorithm for secp256k1

How close is a quantum computer to stealing a Bitcoin key? Shor's discrete-log algorithm on a toy curve, then the resource arithmetic behind Google's 2026 estimates.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **Point addition is a permutation unitary,** exactly like modular multiplication in factoring; the discrete log needs two registers, $f(a,b)=aG+bQ$.
- **Real cryptography in the same language.** `GenerateAsymmetricKeyPair["EllipticCurve"]` uses secp256k1 by default.
- **Precision is a theorem.** Continued fractions recover $r$ only when $2^t\ge2r^2$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **What Google optimized** is the reversible point-addition circuit, not the algorithm: windowed arithmetic and measurement-based uncomputation.
- **Two hardware roads.** Superconducting qubits for minutes-long attacks on transactions in flight, or neutral atoms with qLDPC codes for days-long attacks on dormant wallets.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 10.1: Shor's discrete log on a toy curve

Elliptic-curve arithmetic over GF(p), with `O` the point at infinity, by pattern-matched dispatch.

```wl
EllipticCurve[a_, b_, p_]["Add"][O, q_] := q;
EllipticCurve[a_, b_, p_]["Add"][q_, O] := q;
EllipticCurve[a_, b_, p_]["Add"][{x_, y1_}, {x_, y2_}] /; Mod[y1 + y2, p] == 0 := O;
EllipticCurve[a_, b_, p_]["Add"][{x_, y_}, {x_, y_}] := With[{s = Mod[(3 x^2 + a) PowerMod[2 y, -1, p], p]},
   With[{xr = Mod[s^2 - 2 x, p]}, {xr, Mod[s (x - xr) - y, p]}]];
EllipticCurve[a_, b_, p_]["Add"][{x1_, y1_}, {x2_, y2_}] := With[{s = Mod[(y2 - y1) PowerMod[x2 - x1, -1, p], p]},
   With[{xr = Mod[s^2 - x1 - x2, p]}, {xr, Mod[s (x1 - xr) - y1, p]}]];
EllipticCurve[a_, b_, p_]["Mult"][_, 0] := O;
EllipticCurve[a_, b_, p_]["Mult"][pt_, n_Integer] := With[{add = EllipticCurve[a, b, p]["Add"]},
   Fold[With[{d = add[#1, #1]}, If[#2 == 1, add[d, pt], d]] &, O, IntegerDigits[n, 2]]];
EllipticCurve[a_, b_, p_]["Points"] := Catenate@Table[With[{rhs = Mod[x^3 + a x + b, p]},
     Thread[{x, Select[Range[0, p - 1], Mod[#^2, p] == rhs &]}]], {x, 0, p - 1}];
```

Addition, doubling and scalar multiplication in a few rules.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 10.1 The curve and the public key

$y^2=x^3+x+1$ over GF(23), a generator $G$, its order $r$, and a public key $Q=7G$.

```wl
ec = EllipticCurve[1, 1, 23];
pts = ec["Points"]; r = Length[pts] + 1;
G = SelectFirst[pts, Length@DeleteDuplicates@Table[ec["Mult"][#, k], {k, 0, r - 1}] == r &];
Q = ec["Mult"][G, 7];
{r, G, Q}
```

The attacker knows $G$, $r$ and $Q$, not the 7.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 10.1 The quantum part, simulated exactly

Query $f(a,b)=aG+bQ$ on two registers, measure the point register, Fourier-transform both inputs.

```wl
elts = Prepend[pts, O]; idx = AssociationThread[elts -> Range[0, r - 1]];
t = 7; nn = 2^t;
f = Table[idx[ec["Add"][ec["Mult"][G, a], ec["Mult"][Q, b]]], {a, 0, nn - 1}, {b, 0, nn - 1}];
probs = Total[Table[Abs[Fourier[Boole[Map[# == w &, f, {2}]], FourierParameters -> {0, -1}]]^2,
      {w, 0, r - 1}]]/nn^2;
```

The output distribution peaks near $(j,jk)\cdot2^t/r$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 10.1 Twenty shots

Whenever $j$ is invertible mod $r$, the key is $l\,j^{-1}\bmod r$.

```wl
SeedRandom[1];
samples = RandomChoice[Flatten[probs] -> Tuples[Range[0, nn - 1], 2], 20];
keyFrom[{c_, d_}] := With[{j = Mod[Round[r c/nn], r], l = Mod[Round[r d/nn], r]},
   If[CoprimeQ[j, r], Mod[l PowerMod[j, -1, r], r], Missing[]]];
cands = keyFrom /@ samples;
{cands, Select[DeleteMissing@DeleteDuplicates[cands], ec["Mult"][G, #] === Q &]}
```

Every usable shot returns the key: it comes out of the interference pattern, never from scanning multiples of $G$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 10.2: the resource arithmetic

Google's Toffoli count at window size $w$, and the window that minimizes it.

```wl
T[w_, n_, tpa_] := (tpa + 3 2^w) (2 n/w - 4);
{T[16, 256, 2700000], First@MinimalBy[Range[4, 24], T[#, 256, 2700000] &]}
```

The formula reproduces the reported bound, and the optimum sits where the lookup-table cost starts to rival the point addition.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 10.2 Racing the next block

Blocks arrive as a Poisson process, so a 9-minute attack beats the next block with probability $e^{-9/T_{\rm block}}$, for Bitcoin, Litecoin, Zcash and Dogecoin.

```wl
N[Exp[-9/#]] & /@ {10, 2.5, 1.25, 1}
```

A real threat to Bitcoin transactions in flight, and negligible for the fast chains.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 11. What Makes QuantumFramework Unique

Seven claims, each argued by executed code: symbolic first, one object model for channels, POVMs, Fock space and phase space, measurement records as wires, three engines, ordinary Wolfram Language expressions, and a path to hardware.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **Time-ordered propagators in closed form,** including a Landau-Zener sweep in parabolic-cylinder functions.
- **Exact variational landscapes.** QAOA on $K_4$ gives a closed-form energy surface with a certified maximum.
- **One master-equation solution for every qubit,** with $T_2\le2T_1$ as a theorem.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **Breadth in one object model.** Channels, SIC-POVMs, coherent states, $g^{(2)}$, Wigner against Husimi, multiport interference, tomography.
- **Measurement without collapse.** A coherent error digitized into syndromes and corrected at fidelity exactly 1.
- **Three engines, exact interconversion.** State vector, tensor network and stabilizer tableau agreeing under `SameQ`.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **Language-level operators.** `Exp`, `D`, `Commutator` act on quantum objects; `Reduce` gives the Werner threshold $p<1/2$.
- **Out to hardware.** OpenQASM 3 with exact angles, transpilation, and a GHZ run on IBM hardware.
- **An honest boundary.** Where spectra are not solvable in radicals, QF returns `Root` objects.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 11.1: the master equation solved once

A symbolic density matrix, symbolic relaxation and dephasing rates, and a symbolic drive.

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
$Assumptions = \[Gamma]1 > 0 && \[Gamma]\[Phi] >= 0 && \[CapitalOmega] \[Element] Reals && t >= 0;
\[Sigma]m = QuantumOperator[{{0, 1}, {0, 0}}];
\[Rho]0 = QuantumState[Array[\[Rho], {2, 2}]];
\[Rho]t = QuantumEvolve[\[CapitalOmega]/2 QuantumOperator["Z"],
   {Sqrt[\[Gamma]1] \[Sigma]m, Sqrt[\[Gamma]\[Phi]] QuantumOperator["Z"]}, \[Rho]0, t];
FullSimplify /@ Normal[\[Rho]t["Computational"]["DensityMatrix"]] // MatrixForm
```

Populations relax at $\gamma_1$; coherences precess at $\Omega$ and decay at $\gamma_1/2+2\gamma_\phi$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 11.1 $T_2\le2T_1$ as a theorem

A statement about all admissible rates.

```wl
Reduce[ForAll[{\[Gamma]1, \[Gamma]\[Phi]}, \[Gamma]1 > 0 && \[Gamma]\[Phi] >= 0,
  1/(\[Gamma]1/2 + 2 \[Gamma]\[Phi]) <= 2/\[Gamma]1]]
```

`True`, with equality exactly when pure dephasing vanishes.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 11.2: a coherent error, digitized and corrected

Parity measurements whose outcomes land on record wires 0 and -1, and corrections controlled on those records.

```wl
$Assumptions = True;
id2 = IdentityMatrix[2]; zz = KroneckerProduct[PauliMatrix[3], PauliMatrix[3]];
syn[q1_, q2_] := QuantumMeasurementOperator[
   {QuantumOperator[(KroneckerProduct[id2, id2] + zz)/2, {q1, q2}],
    QuantumOperator[(KroneckerProduct[id2, id2] - zz)/2, {q1, q2}]}, {q1, q2}];
enc = QuantumCircuitOperator[{"CNOT", "CNOT" -> {1, 3}}, "Label" -> "Encode"];
qec = QuantumCircuitOperator[{enc, "RX"[2 \[Theta]] -> 2, syn[1, 2], syn[2, 3],
    "C"["X", {0}, {-1}], "C"["X" -> 2, {0, -1}, {}], "C"["X" -> 3, {-1}, {0}]}];
coded = qec[QuantumTensorProduct[QuantumState[{\[Alpha], \[Beta]}], QuantumState["00"]]];
```

The classical record is part of the circuit.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 11.2 The syndrome distribution

For a symbolic data qubit and a symbolic error angle.

```wl
as = {\[Theta] \[Element] Reals, \[Alpha] \[Element] Reals, \[Beta] \[Element] Reals, \[Alpha]^2 + \[Beta]^2 == 1};
Assuming[as, FullSimplify[ComplexExpand[
   Diagonal[Normal[QuantumPartialTrace[coded, {3, 4, 5}]["DensityMatrix"]]]]]]
```

$\{\cos^2\theta,0,0,\sin^2\theta\}$: the continuous rotation becomes "no error" or "qubit 2 flipped", independent of the data.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 11.2 After the feedforward

Fidelity of the data qubits with the ideal codeword.

```wl
Assuming[as, FullSimplify[ComplexExpand[QuantumSimilarity[QuantumPartialTrace[coded, {1, 2}],
    enc[QuantumTensorProduct[QuantumState[{\[Alpha], \[Beta]}], QuantumState["00"]]], "Fidelity"]]]]
```

Exactly 1 for every input and every angle: the digitization of errors from one symbolic run.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 12. DRAG: Fast Gates on a Transmon

How short can a π pulse on a transmon be? A short, strong pulse leaks into $|2\rangle$ and shifts the transition; a long one suffers $T_1$ and $T_2$. Five levels, a DRAG term, and calibration under the Lindblad equation.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **Leakage seen directly.** A short uncalibrated π pulse leaves a sizable population above the qubit; a long one leaves the flip incomplete.
- **DRAG alone is not the answer.** Leakage is smallest near $\lambda\approx1$, but transfer error is smallest well below it, so amplitude, detuning and $\lambda$ must be calibrated together.
- **The speed limit is set by decoherence, not leakage.** DRAG wins by allowing a shorter pulse.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **The model checks itself.** The averaged gate error tracks the transfer error, and three levels underestimate leakage while five and six agree.
- **All in one object.** `AnnihilationOperator`, a parametrized `QuantumOperator`, and `QuantumEvolve` for Schrodinger and Lindblad dynamics.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 12.1: leakage, and what DRAG does to it

The driven anharmonic oscillator with a DRAG term, and final populations after a pulse.

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
envelope = A (Pi/\[Tau]) (1 - Cos[2 Pi t/\[Tau]]);
hamiltonian[c_] := (\[Alpha]/2) (c["Dagger"] @ c["Dagger"] @ c @ c) + \[Delta] (c["Dagger"] @ c) +
   (envelope/2) (c + c["Dagger"]) - (\[Lambda] D[envelope, t]/(2 \[Alpha])) (I (c["Dagger"] - c));
anharmonicity = -2 Pi 0.2;
finalPopulations[n_, {amplitude_, drag_, detuning_, length_}] := QuantumEvolve[
   QuantumOperator[hamiltonian[AnnihilationOperator[n]], "Parameters" -> {\[Alpha], A, \[Lambda], \[Delta], \[Tau]}][
    anharmonicity, amplitude, drag, detuning, length],
   FockState[0, n], {t, 0, length}, Method -> "ExplicitRungeKutta", PrecisionGoal -> 12, AccuracyGoal -> 14][length][
  "ProbabilitiesList"];
```

One function returns the populations for any amplitude, DRAG, detuning and length.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 12.1 Uncalibrated π pulses

Populations of $|0\rangle,|1\rangle,|2\rangle$ after 8 ns and 35 ns.

```wl
{finalPopulations[5, {1, 0, 0, 8}][[;; 3]], finalPopulations[5, {1, 0, 0, 35}][[;; 3]]}
```

The short pulse loses population out of the qubit; the long one keeps it but leaves part of it in $|0\rangle$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 12.1 Scanning the DRAG coefficient

Leakage and transfer error at 8 ns.

```wl
ListLogPlot[Transpose@Table[With[{p = finalPopulations[5, {1, drag, 0, 8}]},
    {{drag, Total[p[[3 ;;]]]}, {drag, 1 - p[[2]]}}], {drag, 0, 2, 0.1}],
 Joined -> True, PlotMarkers -> Automatic, PlotLegends -> {"leakage", "transfer error"},
 Frame -> True, FrameLabel -> {"\[Lambda]", "population"}]
```

Suppressing leakage and completing the flip are different targets.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 12.1 Is the truncation safe?

Leakage with three to six levels.

```wl
Table[{n, Total[finalPopulations[n, {1, 0, 0, 8}][[3 ;;]]]}, {n, 3, 6}]
```

Three levels miss a large part of the leakage; five and six agree.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 12.2: calibrating the 8 ns pulse

Minimize the transfer error over amplitude and detuning, then over the DRAG coefficient as well. This search takes about five minutes.

```wl
err[a_?NumericQ, l_?NumericQ, d_?NumericQ] := 1 - finalPopulations[5, {a, l, d, 8}][[2]];
{NMinimize[{err[a, 0, d], 0 < a < 2 && -0.5 <= d <= 0.5}, {{a, 0.8, 1.2}, {d, -0.1, 0.1}},
   Method -> {"DifferentialEvolution", "SearchPoints" -> 10, "RandomSeed" -> 0}],
 NMinimize[{err[a, l, d], 0 < a < 2 && 0 <= l <= 3 && -0.5 <= d <= 0.5}, {{a, 0.8, 1.2}, {l, 0, 2}, {d, -0.1, 0.1}},
   Method -> {"DifferentialEvolution", "SearchPoints" -> 10, "RandomSeed" -> 0}]}
```

Without DRAG the short pulse stalls at a percent-level error; with it the error drops several-fold, with a nonzero detuning compensating the drive-induced shift.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 13. QuantumInterferometer: Symbolic Behaviour

QuantumFramework's linear-optics object, with beam splitters, phase shifters, squeezers and loss, when the inputs are symbols: angles, amplitudes, squeezing, loss rates, and above all photon numbers.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **Symbolic photon numbers.** A beam splitter's amplitude comes back as a Wigner $D$-function, the Schwinger realization of SU(2).
- **Hong-Ou-Mandel for every $n$.** On a balanced splitter $|n,n\rangle\to|k,2n-k\rangle$ is $d^n_{0,k-n}(\pi/2)$, and every odd split vanishes.
- **Coherent, squeezed and lossy light, exactly,** including the binomial loss law for symbolic $n$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **A physical-domain contract.** Photon numbers are nonnegative integers, angles real, $0\le\eta\le1$; results carry no conditions restating that domain.
- **Stated limits.** What is not supported symbolically is listed explicitly.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 13.1: a beam splitter with symbolic photon numbers

The amplitude from $|n_1,n_2\rangle$ to $|t_1,n_1+n_2-t_1\rangle$.

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
bs = QuantumInterferometer[{"BS"[th, ph]}];
bs["Amplitude", {n1, n2} -> {t1, n1 + n2 - t1}]
```

A two-mode beam splitter is a rotation of a spin of size $(n_1+n_2)/2$, and the photon numbers are its magnetic quantum numbers.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 13.1 Hong-Ou-Mandel for every $n$

A balanced splitter with $n$ photons in each port.

```wl
hom = QuantumInterferometer[{"BS"[Pi/4, 0]}]["Amplitude", {n, n} -> {k, 2 n - k}];
{hom, Table[hom /. {n -> m, k -> j}, {m, 1, 3}, {j, 1, 2 m - 1, 2}], Table[Abs[hom /. {n -> 2, k -> j}]^2, {j, 0, 4}]}
```

One formula for all $n$; every odd split is exactly zero, and $n=1$ is the HOM dip.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 13.1 Two splitters in a chain

A held sum over the internal photon number, checked against the concrete query.

```wl
chain = QuantumInterferometer[{"BS"[th1, ph1], "BS"[th2, ph2]}];
FullSimplify[Activate[chain["Amplitude", {n, 0} -> {t, n - t}] /. {n -> 2, t -> 1}] -
  chain["Amplitude", {2, 0} -> {1, 1}]]
```

Activated, the symbolic sum agrees exactly.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 13.2: loss and coherent light

Loss on $n$ photons with $n$ symbolic, and a coherent pair through a general splitter.

```wl
{QuantumInterferometer[{"Loss"[eta]}]["Probability", {n} -> {k}], bs["Coherent"[{a, b}]]}
```

The binomial loss law, and a coherent state that stays coherent with amplitudes rotated by the transfer matrix.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 13.2 Squeezed light

Two squeezed vacua on a balanced splitter, and the photon statistics of one squeezed vacuum.

```wl
{QuantumInterferometer[{"S"[x1] -> 1, "S"[x2] -> 2, "BS"[Pi/4, 0]}]["Amplitude", {0, 0} -> {1, 1}],
 QuantumInterferometer[{"S"[xi]}]["Probabilities", {0}, "Cutoff" -> 4, Assumptions -> xi > 0]}
```

Coincidences vanish for equal squeezing and peak for opposite squeezing; a squeezed vacuum populates only even photon numbers.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 14. Quantum Optics with QuantumInterferometer

Fourteen linear-optics experiments, each ending in a checked prediction. Fock amplitudes are permanents, squeezed amplitudes hafnians, and loss acts on probabilities directly.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **Exact zeros are the signal.** A Fourier tritter forbids whole families of patterns, which floating point would blur.
- **Laws for every photon number.** Twin-Fock inputs never give odd outputs, with the distribution closed in binomials for all $n$.
- **Tunable suppression.** $|2,2\rangle\to|2,2\rangle$ vanishes at two unbalanced reflectivities.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **Squeezed light, symbolically.** Opposite squeezers on a splitter give a two-mode squeezed vacuum.
- **Heralding and gates,** up to the KLM nonlinear sign gate.
- **Cross-checks everywhere.** Permanents against Fock-space circuits, native loss against a dilation, samples against exact distributions.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 14.1: the KLM nonlinear sign gate, exactly

A three-mode real unitary with one auxiliary photon; success is heralded by the auxiliary pattern $\{1,0\}$.

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
nsMatrix = {{1 - Sqrt[2], 2^(-1/4), Sqrt[3/Sqrt[2] - 2]}, {2^(-1/4), 1/2, (1 - Sqrt[2])/2},
   {Sqrt[3/Sqrt[2] - 2], (1 - Sqrt[2])/2, Sqrt[2] - 1/2}};
nsGate = QuantumInterferometer[nsMatrix];
{RootReduce[nsMatrix . Transpose[nsMatrix]] == IdentityMatrix[3],
 RootReduce[Table[nsGate["Amplitude", {n, 1, 0} -> {n, 1, 0}], {n, 0, 3}]]}
```

Amplitudes $\tfrac12,\tfrac12,-\tfrac12$ for 0, 1, 2 signal photons: a sign flip on $|2\rangle$ alone, which no phase shifter can do, and the three-photon amplitude shows where the gate stops working.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 14.2: squeezed pairs

Opposite squeezers and a balanced splitter, with $r$ symbolic.

```wl
symbolicPair = QuantumInterferometer[{"S"[r], "S"[-r] -> 2, "BS"[Pi/4]}];
{FullSimplify[Table[symbolicPair["Probability", {0, 0} -> {n, n}], {n, 0, 2}], r > 0],
 FullSimplify[symbolicPair["Probability", {0, 0} -> {2, 0}], r > 0]}
```

Photons come only in pairs, with the geometric law $\operatorname{sech}^2r\tanh^{2n}r$: a two-mode squeezed vacuum.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 14.2 Tunable four-photon suppression

The return amplitude of $|2,2\rangle$ on a splitter of angle $\theta$, and the reflectivities that kill it.

```wl
fourBS = QuantumInterferometer[{"BS"[\[Theta]]}];
fa = FullSimplify[fourBS["Amplitude", {2, 2} -> {2, 2}], \[Theta] \[Element] Reals];
{FullSimplify[fa == 1 - 6 Sin[\[Theta]]^2 Cos[\[Theta]]^2, \[Theta] \[Element] Reals],
 R /. Solve[1 - 6 R (1 - R) == 0 && 0 < R < 1, R, Reals]}
```

Zero at $R=(3\pm\sqrt3)/6$: unlike HOM, four-photon cancellation needs an unbalanced splitter.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 15. A Formal Semantics for QuantumFramework

What does a QF object mean? A coefficient array means nothing until its subsystem dimensions, its output/input split and its bases are fixed, and QF stores exactly those.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **Numbers versus maps.** Six numbers can be a qudit, a qubit and a qutrit, or a $2\times3$ map; stored $C$ and represented $A$ satisfy $A=F_{\rm out}CF_{\rm in}^{-1}$.
- **Bases are reconciled before composing.** Contracting stored numbers can give the wrong operator.
- **One operator-space model** for channels, evolution and phase space.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **One isometry, three readings.** The extra index of a Stinespring isometry can be kept, traced, or read.
- **Honest edges.** Constructors check shapes, not physics, and the stated limits are explicit.
- **Placement in the literature,** against categorical quantum mechanics, circuit languages and operator-algebraic semantics.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 15.1: stored numbers against the represented map

The operator $R$ whose stored coefficients are $\mathrm{diag}(1,-1)$ in the $X$ basis.

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
R = QuantumOperator[{{1, 0}, {0, -1}}, "PauliX"];
{Normal@R["Matrix"], Normal@R["MatrixRepresentation"],
 Normal@QuantumOperator[QuantumOperator["X"], "PauliX"]["Matrix"]}
```

Declaring a basis turns the numbers $\mathrm{diag}(1,-1)$ into the operator $X$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 15.1 Applying $R$ to $|0\rangle$

Stored vector and computational vector.

```wl
{Normal@R[QuantumState["0"]]["StateVector"], Normal@R[QuantumState["0"]]["Computational"]["StateVector"]}
```

$(1,-1)/\sqrt2$ in the $X$ basis is $|1\rangle$ computationally.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 15.1 Why composition needs the basis

$Z$ after $R$, against contracting the stored numbers naively.

```wl
{Normal@QuantumOperator["Z"][R]["MatrixRepresentation"], PauliMatrix[3] . Normal[R["Matrix"]]}
```

QF returns $ZX$; the naive contraction gives the identity, the wrong operator.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 15.2: maps on density matrices

The Liouvillian of $X$ acting on $|0\rangle\langle0|$.

```wl
Normal@QuantumOperator["Liouvillian"["X"]][QuantumState["0"]]["DensityMatrix"]
```

A map on density matrices, stored as an operator, returns the trace-zero commutator.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 15.2 A channel and a measurement

A list read as Kraus operators, and a list read as POVM effects.

```wl
ch = QuantumChannel[{Sqrt[3/4] IdentityMatrix[2], Sqrt[1/4] PauliMatrix[1]}];
m = QuantumMeasurementOperator[{{{3/4, 0}, {0, 1/4}}, {{1/4, 0}, {0, 3/4}}}];
{Normal@ch[QuantumState["0"]]["DensityMatrix"], ch["TracePreservingQ"], m[QuantumState["0"]]["ProbabilitiesList"]}
```

The same kind of data becomes a trace-preserving channel or a set of Born probabilities, depending on how it is read.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### 15.2 Branches are unnormalized

Measure $X$ on $R|0\rangle$.

```wl
mx = QuantumMeasurementOperator["X"][R[QuantumState["0"]]];
{mx["ProbabilitiesList"], Norm[Normal[#["StateVector"]]] & /@ mx["States"]}
```

Each branch has norm $\sqrt{p_m}$: branch, probability and conditional state are three distinct objects.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## 16. Why Boson Sampling Is Hard

Three tasks kept apart, computing one amplitude, collecting every probability, and drawing one sample, because the hardness argument lives in the gap between them.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice

- **Three particle statistics, one table.** Bosons (permanent), distinguishable particles (permanent of $|W|^2$) and fermions (determinant) on one splitter.
- **The permanent derived, not postulated.** Output amplitudes are coefficients of $\prod_k(\sum_jW_{jk}x_j)^{s_k}$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What to notice (continued)

- **Cost depends on occupation structure.** Collision-free and paired amplitudes grow exponentially with $n$; a fully bunched one stays essentially free.
- **Hard to compute is not hard to sample.** Distinguishable photons have permanent probabilities yet sample trivially; the Aaronson-Arkhipov argument is what turns amplitude cost into sampling hardness.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 16.1: bosons, classical particles and fermions

Three probability rules on one balanced beam splitter.

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
modeList[occ_List] := Catenate@MapIndexed[ConstantArray[First[#2], #1] &, occ];
bosonProb[w_, s_, t_] := Abs[Permanent[w[[modeList[t], modeList[s]]]]]^2/Times @@ Factorial[Join[s, t]];
distProb[w_, s_, t_] := Permanent[(Abs[w]^2)[[modeList[t], modeList[s]]]]/Times @@ Factorial[t];
fermionProb[w_, s_, t_] := Abs[Det[w[[modeList[t], modeList[s]]]]]^2;
Wbs = QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2]["TransferMatrix"];
Table[Through[{bosonProb, distProb, fermionProb}[Wbs, {1, 1}, t]], {t, {{2, 0}, {1, 1}, {0, 2}}}]
```

Bosons bunch, classical particles split, fermions antibunch: the only difference is the sign attached to exchanging two histories.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Example 16.2: what one amplitude costs

Collision-free, paired and fully bunched inputs in random networks with $m=2n$ modes, one amplitude each.

```wl
occupations[n_] := PadRight[#, 2 n] & /@ <|"collision-free" -> ConstantArray[1, n],
    "pairs" -> Append[ConstantArray[2, Quotient[n, 2]], Mod[n, 2]], "bunched" -> {n}|>;
randomInterferometer[n_] := BlockRandom[SeedRandom[1107 + 1000 n]; QuantumInterferometer["Random"[2 n]]];
timeAmplitude[q_, s_List] := With[{t = RotateRight[s, Total[s]]},
   ClearSystemCache[]; First@AbsoluteTiming[q["Amplitude", s -> t]]];
AssociationMap[Function[n, With[{q = randomInterferometer[n]}, timeAmplitude[q, #] & /@ occupations[n]]],
 {6, 12, 18, 20}]
```

Generic amplitudes pay the exponential cost of Ryser and Glynn while structured ones do not, which is why the hardness argument needs random networks.
