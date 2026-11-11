# Talk Topic Summaries (review draft, all 17 topics)

Each topic gets a summary of what should catch a physicist's attention, then one or two showcase examples. Every example below is a condensed, self-contained version of cells in the source essay and was run in a fresh kernel with `wolframscript` before being written here. The examples use plain matrices only; neither essay calls QuantumFramework.

## 1. Continuous Quantum Error Correction: Syndrome Estimation and Feedback

Source: `/Users/mohammadb/Documents/GitHub/mbahram/open-quantum-systems/Examples on Open Quantum Systems/three-qubit-bitflip-continuous-qec.md`

### Summary

The three-qubit bit-flip code is protected in continuous time: the two stabilizers $Z_1Z_2$ and $Z_2Z_3$ are measured weakly and continuously, a classical controller turns the two noisy records into an estimate, and the estimate drives corrective $X_i$ pulses. The essay simulates the physical qubits and the controller's belief as two separate density matrices fed by the same record, which is what an experiment such as Livingston et al. (2022) actually runs.

What should catch a physicist's attention:

- **The syndrome is readable even though each sample is useless.** In one time step the stabilizer signal sits about eight times below the shot noise, yet the four syndrome probabilities track every bit flip. With feedback off they obey a closed, four-number stochastic equation (a Wonham filter), $dp_j=\gamma(1-4p_j)\,dt+2\sqrt{\eta\kappa}\,p_j\sum_k(s_{jk}-\bar s_k)\,dW_k$, implemented as an exact Bayes update in three lines and checked against the full eight-dimensional stochastic master equation.
- **Feedback breaks the four-number description.** Once the drive $H_{\rm fb}=\lambda\sum_i p_iX_i$ is on, the sector populations couple to coherences between sectors, so the controller must carry the full $8\times8$ estimate, not four probabilities.
- **The controller does not need to know the logical state.** Starting the estimate at the maximally mixed code state while the qubits hold $|0_L\rangle$, the physical fidelity stays high even though the estimate gives $|0_L\rangle$ and $|1_L\rangle$ equal weight: correction acts on syndromes, not on the logical amplitudes.
- **A wrong detector-efficiency model costs protection, and the essay separates two effects.** It measures the gap between estimated and true fidelity, and separately the loss of physical fidelity caused by the wrong model, by pairing runs on identical noise.
- **Staying in the code space is not keeping the qubit.** The residual fidelity loss splits exactly into leakage out of the code space plus population in the orthogonal logical state $|1_L\rangle$, which no stabilizer can see.
- **A feedback that protects the code space exactly.** Monitoring the error channels themselves and feeding back $F_i=i\sqrt\gamma\,[\Pi_{\rm code},X_i]$ makes $(c_i-iF_i)\Pi_{\rm code}=0$, so at unit efficiency every encoded state is a fixed point of the averaged dynamics (Wiseman-Mancini-Wang feedback master equation). Below unit efficiency an extra $\frac{1-\eta}{\eta}\mathcal D[F_i]$ term appears, and at very poor efficiency this scheme does worse than no feedback at all.

### Example 1.1: tracking bit flips from two noisy records

The toolkit: Pauli matrices, the stabilizers, the four sector projectors, and the positivity-preserving measurement step of Rouchon and Ralph (2015), which conditions the state on a record increment without ever leaving the set of valid states.

```wl
{id2, X, Y, Z} = Table[PauliMatrix[j], {j, 0, 3}];
op3[o_, i_] := KroneckerProduct @@ ReplacePart[{id2, id2, id2}, i -> o];
ket[bits_] := Normal@SparseArray[FromDigits[bits, 2] + 1 -> 1, 8];
S1 = op3[Z, 1] . op3[Z, 2]; S2 = op3[Z, 2] . op3[Z, 3];
signs = {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
proj = (IdentityMatrix[8] + #[[1]] S1) . (IdentityMatrix[8] + #[[2]] S2)/4 & /@ signs;
```

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

Start in $|000\rangle$, let each qubit flip at rate $\gamma$, watch both stabilizers, and keep only the record:

```wl
kappa = 1.; eta = 1.; gamma = 0.05; dt = 0.004;
run = trajectory[KroneckerProduct[#, #] &@ket[{0, 0, 0}], 0 IdentityMatrix[8],
   Sqrt[kappa] {S1, S2}, {eta, eta}, Sqrt[gamma] (op3[X, #] & /@ Range[3]), dt, 40., 3];
```

The filter: relax each sector toward $1/4$ under the flips, then multiply by the likelihood of the new pair of record increments and renormalize:

```wl
syndromeFilter[p0_, record_] := With[{c = 2 Sqrt[eta kappa]},
   FoldList[With[{w = (1/4 + Exp[-4 gamma dt] (#1 - 1/4)) Exp[c signs . #2]}, w/Total[w]] &,
    p0, record]];
pFilter = syndromeFilter[{1., 0., 0., 0.}, run["record"]];
```

```wl
ListLinePlot[Transpose[pFilter], DataRange -> {0, 40}, PlotRange -> {-0.05, 1.05},
 PlotLegends -> {"code", "X1 flipped", "X2 flipped", "X3 flipped"}, Frame -> True,
 FrameLabel -> {"time", "syndrome probability"}]
```

The four probabilities hand certainty from sector to sector as flips accumulate. Four numbers, updated by Bayes' rule, follow the same syndrome the full eight-dimensional conditional state carries:

```wl
Max@Abs[pFilter - Map[Function[r, Re@Tr[# . r] & /@ proj], run["states"]]]
```

The two agree up to the finite-step difference between two discretizations of one continuous filter.

### Example 1.2: exact protection by Markovian feedback, and its price in detector efficiency

The averaged dynamics needs only the Liouvillian:

```wl
liouvillian[H_, leaks_List, d_] := With[{id = IdentityMatrix[d]},
   -I (KroneckerProduct[H, id] - KroneckerProduct[id, Transpose[H]]) +
    Total[Function[c, KroneckerProduct[c, Conjugate[c]] -
        (KroneckerProduct[ConjugateTranspose[c] . c, id] +
           KroneckerProduct[id, Transpose[ConjugateTranspose[c] . c]])/2] /@ leaks]];
evolve[H_, leaks_, rho0_, t_] := With[{d = Length[rho0]},
   ArrayReshape[MatrixExp[liouvillian[H, leaks, d] t] . Flatten[rho0], {d, d}]];
```

Keep the error rate symbolic, build the feedback operators from the code projector, and check the two identities that make the code space invariant, for every $\gamma$:

```wl
codeProj = First[proj];
cs = Sqrt[\[Gamma]] (op3[X, #] & /@ Range[3]);
F = I (codeProj . # - # . codeProj) & /@ cs;
{Union@Flatten@Simplify@MapThread[#1 . #2 + #2 . #1 &, {cs, F}],
 Union@Flatten@Simplify[(# . codeProj) & /@ (cs - I F)]}
```

Both lists are `{0}`: $c_iF_i+F_ic_i=0$ and $(c_i-iF_i)\Pi_{\rm code}=0$. Now the price of an imperfect detector, comparing the code-space population after $t=15$ with and without feedback:

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

At $\eta=1$ the code space is held exactly. As the efficiency drops, the $\frac{1-\eta}{\eta}\mathcal D[F_i]$ term feeds amplified detector noise back onto the qubits, and at a few percent efficiency the feedback leaves the code population below the no-feedback line: an inefficient loop is worse than no loop.

## 2. A Computational Guide to Open Quantum Systems: Lindblad Dynamics, Quantum Trajectories, and Feedback Control

Source: `/Users/mohammadb/Documents/GitHub/mbahram/open-quantum-systems/Examples on Open Quantum Systems/open-system-simulation-catalog.md`

### Summary

Twenty open-system problems, each built from about ten shared functions on plain matrices: the Lindbladian solved two independent ways (Liouvillian exponential and ODE), a steady-state finder that handles degenerate null spaces, a positivity-preserving measurement step, a trajectory generator, and a diagnostic that splits the gap between a trajectory ensemble and the master equation into Monte-Carlo scatter and time-step bias. Every example states its equation, computes, and checks itself against a closed form or a second method.

What should catch a physicist's attention:

- **Decoherence is the price of discarding the record.** Under ideal monitoring each conditioned trajectory stays pure while the ensemble average mixes; the essay shows this for homodyne, heterodyne and for a three-channel atom where only the complete record keeps the run pure.
- **Why the naive stochastic master equation fails.** One large noise kick in the textbook update gives a negative eigenvalue; the Kraus-form step lands exactly on the boundary of state space instead.
- **Cats.** A cavity cat's coherence dies at a rate set by the square of its separation (checked against the exact Kraus solution of photon loss); in the Haroche experiment the cat exists only once the atom is read in the right basis, and the two-atom correlation tracks the coherence; in dispersive readout the conditional cavity is a cat that turns into a mixture as the pointers separate, with a closed form checked against simulation.
- **Quantum Brownian motion.** The high-temperature Caldeira-Leggett equation is not of Lindblad form, so the essay watches the uncertainty bound directly; a warm bath kills a cat faster, and at fixed delay the surviving coherence peaks at a modest cat size.
- **One leak, three unravellings.** Photon counting, homodyne and heterodyne give different histories with the same average; the click record shows antibunching in its waiting-time distribution, and the homodyne section keeps a strict ledger of what is observed (the record) and what is inferred (state, signal, innovation).
- **Two thresholds where one is expected.** For a driven, watched qubit the Liouvillian eigenvalues are closed-form, with an exceptional point at $k=\Omega$ and a Zeno freeze rate $\Omega^2/2k$. In the charge qubit read by a quantum point contact, the spectral peaks merge already at $k=\Omega/\sqrt2$, before the eigenvalues coalesce at $k=\Omega$.
- **Collapse and Born's rule as a martingale.** Under pure measurement every run tends to a pole while the mean stays fixed, which forces $p_\pm=(1\pm z_0)/2$.
- **Quantum Kalman filter.** Watching position squeezes the conditional variance below the zero-point value, and the steady state sits exactly at the Heisenberg floor for every measurement strength, in closed form.
- **Feedback in three guises.** Wiseman-Milburn feedback purifies a maximally mixed qubit onto a dark state; Jacobs rapid purification beats the unsteered measurement on the mean at fixed time yet loses on mean first-passage time (the Wiseman-Ralph point, read off distributions); feedback cooling shows watching alone heats, an optimal gain exists, and the conditional covariance sets a floor no gain can cross.
- **Thermalization derives the canonical ensemble.** Detailed balance between $\gamma(n_T+1)$ and $\gamma n_T$ forces the whole thermal distribution, not only its mean.
- **Linear versus nonlinear trajectories.** The norm of the linear stochastic Schrodinger state is the likelihood ratio $dP/dQ$; equal-weight averaging gives a bias no ensemble size removes, and the effective sample size decays with time.
- **Mollow triplet and squeezed fluorescence.** From one Liouvillian, by quantum regression: a single line at weak drive, three lines at strong drive, and a homodyne spectrum that dips below shot noise.

### Example 2.1: the quantum Zeno crossover as an exceptional point

```wl
{id2, X, Y, Z} = Table[PauliMatrix[j], {j, 0, 3}];
liouvillian[H_, leaks_List, d_] := With[{id = IdentityMatrix[d]},
   -I (KroneckerProduct[H, id] - KroneckerProduct[id, Transpose[H]]) +
    Total[Function[c, KroneckerProduct[c, Conjugate[c]] -
        (KroneckerProduct[ConjugateTranspose[c] . c, id] +
           KroneckerProduct[id, Transpose[ConjugateTranspose[c] . c]])/2] /@ leaks]];
```

A qubit Rabi-driven at $\Omega$ while $\hat\sigma_z$ is watched at strength $k$, fully symbolic:

```wl
Eigenvalues@FullSimplify[liouvillian[(\[CapitalOmega]/2) X, {Sqrt[k] Z}, 2], \[CapitalOmega] > 0 && k > 0]
```

The spectrum is $\{0,\,-2k,\,-k\pm\sqrt{k^2-\Omega^2}\}$: the Rabi pair oscillates for $k<\Omega$, coalesces at the exceptional point $k=\Omega$, and splits into two real rates beyond it. The slow one is the Zeno freeze:

```wl
Normal@Series[-k + Sqrt[k^2 - \[CapitalOmega]^2], {k, Infinity, 1}]
```

```wl
Plot[{Re@Sqrt[1 - r^2], r - Re@Sqrt[r^2 - 1]}, {r, 0, 3}, Frame -> True, GridLines -> {{1}, None},
 PlotLegends -> {"Rabi frequency / \[CapitalOmega]", "slowest decay rate / \[CapitalOmega]"},
 FrameLabel -> {"k/\[CapitalOmega]", "Liouvillian eigenvalue / \[CapitalOmega]"}]
```

Watching harder speeds relaxation only up to $k=\Omega$; past it the slowest rate falls as $\Omega^2/2k$ and the qubit freezes.

### Example 2.2: the Mollow triplet and fluorescence squeezed below shot noise

Two-time correlations run on the same Liouvillian as the state (quantum regression), so a spectrum is one linear solve per frequency. At $\omega=0$ the system is singular, and the trace condition picks the physical branch:

```wl
lower = (X - I Y)/2;
steady[H_, leaks_] := With[{v = ArrayReshape[First@NullSpace[liouvillian[H, leaks, 2]], {2, 2}]}, v/Tr[v]];
resolventSpectrum[big_, seedMat_, readout_, tones_] := With[{n = Length[big], d = Length[readout]},
   Table[2 Re@Tr[readout . ArrayReshape[If[Abs[w] < 10^-9,
       LeastSquares[Join[N[-big], {1. Flatten[IdentityMatrix[d]]}], Join[Flatten[seedMat], {0.}]],
       LinearSolve[I w IdentityMatrix[n] - big, Flatten[seedMat]]], {d, d}]], {w, tones}]];
```

The inelastic fluorescence spectrum of a resonantly driven atom ($\gamma=1$), weak against strong drive:

```wl
mollow[om_, tones_] := With[{big = liouvillian[(om/2) X, {lower}, 2], rss = steady[(om/2) X, {lower}]},
   With[{dm = lower - Tr[lower . rss] id2}, resolventSpectrum[big, dm . rss, ConjugateTranspose[dm], tones]]];
tones = Range[-9, 9, 0.05];
{weak, strong} = mollow[#, tones] & /@ {0.6, 6.};
ListLinePlot[{Transpose[{tones, weak/Max[weak]}], Transpose[{tones, strong/Max[strong]}]},
 Frame -> True, PlotRange -> All, PlotLegends -> {"\[CapitalOmega] = 0.6\[Gamma]", "\[CapitalOmega] = 6\[Gamma]"},
 FrameLabel -> {"frequency - drive frequency", "inelastic spectrum (scaled to peak)"}]
```

Strong driving dresses the atom and splits one line into three, with sidebands at $\pm\Omega$. A homodyne detector reading the in-phase quadrature sees something a spectrometer cannot; at zero frequency:

```wl
homodyneS0[om_] := Module[{big = liouvillian[(om/2) X, {lower}, 2], rss = steady[(om/2) X, {lower}], cth, mq},
   cth = -I lower; mq = cth + ConjugateTranspose[cth];
   1 + First@resolventSpectrum[big, cth . rss + rss . ConjugateTranspose[cth] - Re@Tr[mq . rss] rss,
      mq - Re@Tr[mq . rss] id2, {0.}]];
Max@Table[Abs[homodyneS0[om] - (1 + 8 om^2 (4 om^2 - 1)/(1 + 2 om^2)^3)], {om, 0.05, 1.2, 0.05}]
```

The numerical spectrum is the closed form $S(0)=1+8\Omega^2(4\Omega^2-1)/(1+2\Omega^2)^3$ to machine precision, so the deepest squeezing follows exactly:

```wl
Minimize[{1 + 8 \[CapitalOmega]^2 (4 \[CapitalOmega]^2 - 1)/(1 + 2 \[CapitalOmega]^2)^3, \[CapitalOmega] > 0}, \[CapitalOmega]]
```

The fluorescence noise falls below the vacuum level for every drive $\Omega<\gamma/2$, deepest at $\Omega=\tfrac12\sqrt{3-\sqrt7}\,\gamma$, where $S(0)=(323-56\sqrt7)/243$: nonclassical light from a single driven atom, visible only to phase-sensitive detection.

## 3. Wong-Zakai Regularization for Continuous Quantum Measurement with NDSolve

Source: `/Users/mohammadb/Documents/GitHub/mbahram/continuous-measurement/Discrete Space/NDSolve vs Ito/Regularized-SDE-NDSolve-Simulation.md`

### Summary

A continuously monitored, driven qubit (homodyne readout of $\hat\sigma_z$, with relaxation, dephasing and phase backaction) is written as one four-component stochastic equation: three Bloch components and the detector's accumulated readout $Q$, all kicked by the same scalar noise. `NDSolve` cannot integrate white noise, so the essay converts the Itô equation to its Stratonovich form, replaces the Wiener path by a continuous piecewise-cubic interpolant, and hands `NDSolve` the resulting random ODE. By the Wong-Zakai theorem this converges to the Stratonovich (hence the original Itô) solution as the noise grid is refined, which works here because the noise is scalar.

What should catch a physicist's attention:

- **The exact algebra is checked symbolically, with all six rates symbolic.** `ItoProcess[StratonovichProcess[...]]` converts the Stratonovich field back and lands exactly on the affine Lindblad drift; the Lindblad generator projected on Pauli components reproduces the same drift.
- **A positivity condition derived, not assumed.** The Itô generator of $r^2$ on the Bloch sphere is $(1-z^2)(\Gamma_{CI}+\Gamma_{BA}-2(\Gamma_d+\gamma_\phi))-\gamma_1(1-z)^2$, so the trajectory stays a density matrix if and only if $\Gamma_{CI}+\Gamma_{BA}\le2(\Gamma_d+\gamma_\phi)$, i.e. detector efficiency $\eta\le1+\gamma_\phi/\Gamma_d$. Violate it and a path leaves the ball while the equation still integrates; saturate it and a pure state stays exactly on the sphere along every single path.
- **The record is the noise.** Subtracting the integrated signal $\sqrt{\Gamma_{CI}}\int z\,dt$ from one record returns the very Wiener path that kicked the state.
- **A closed law tying state to record.** With drive and relaxation off, $z(t)=\tanh(\operatorname{arctanh}z_0+\sqrt{\Gamma_{CI}}\,Q(t))$ on every path: the monitored coordinate is the readout read through a $\tanh$.
- **Honest numerics.** Means are checked against the closed-form Lindblad solution, spreads against the built-in order-3/2 `StratonovichProcess` integrator; the essay separates the two time steps (noise grid versus solver step) and shows the noise grid is the one that matters, including in a strong-readout regime where the Itô-Stratonovich correction is large.

### Example 3.1: the Itô-Stratonovich conversion, the Lindblad drift and the positivity boundary, all symbolic

Encode the Stratonovich drift, the diffusion and the Itô (Lindblad) drift on a rate vector $\{\Gamma_{CI},\Gamma_{BA},\Gamma_d,\gamma_\phi,\gamma_1,\Omega_x\}$:

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

Convert the Stratonovich process back to Itô form and compare with the Lindblad drift plus the readout row:

```wl
Assuming[{rs[[1]] >= 0, rs[[2]] >= 0}, FullSimplify[
  (ItoProcess[StratonovichProcess[{driftStrat[{x, y, z, q}, rs], List /@ diffStrat[{x, y, z, q}, rs], {x, y, z, q}},
       {{x, y, z, q}, {x0, y0, z0, q0}}, t]]["Drift"] /. {x[t] -> x, y[t] -> y, z[t] -> z, q[t] -> q})
   == Append[driftIto[{x, y, z}, rs], Sqrt[rs[[1]]] z]]]
```

`True`. That affine drift really is the master equation: build the Lindbladian with drive, relaxation and total dephasing, project on the Pauli matrices, subtract:

```wl
With[{sx = PauliMatrix[1], sy = PauliMatrix[2], sz = PauliMatrix[3], sm = {{0, 1}, {0, 0}}},
 With[{rho = (IdentityMatrix[2] + x sx + y sy + z sz)/2},
  With[{diss = Function[L, L . rho . ConjugateTranspose[L] -
       (ConjugateTranspose[L] . L . rho + rho . ConjugateTranspose[L] . L)/2]},
   Simplify[(Tr[(-I (rs[[6]]/2) (sx . rho - rho . sx) + rs[[5]] diss[sm] +
          ((rs[[4]] + rs[[3]])/2) diss[sz]) . #] & /@ {sx, sy, sz}) - driftIto[{x, y, z}, rs]]]]]
```

`{0, 0, 0}`. Last, the flux of $r^2$ through the Bloch sphere, the quantity whose sign decides whether a trajectory can leave the set of density matrices:

```wl
With[{b = Most@diffStrat[{x, y, z, q}, rs]},
 FullSimplify[(2 {x, y, z} . driftIto[{x, y, z}, rs] + b . b) /.
   {x -> Sqrt[1 - z^2] Cos[\[Phi]], y -> Sqrt[1 - z^2] Sin[\[Phi]]}]]
```

The result is $(1-z^2)(\Gamma_{CI}+\Gamma_{BA}-2\Gamma_d-2\gamma_\phi)-\gamma_1(1-z)^2$, independent of the azimuth. Near the ground pole the relaxation term vanishes quadratically while the measurement term is linear, so the flux points inward everywhere exactly when $\Gamma_{CI}+\Gamma_{BA}\le2(\Gamma_d+\gamma_\phi)$: a detector cannot extract more information than the dephasing it causes pays for.

### Example 3.2: NDSolve integrating a monitored qubit and its readout

Sample a Wiener path, interpolate it with a cubic, and give its derivative to `NDSolve` as a classical forcing of the Stratonovich field:

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

One trajectory of a driven, relaxing, quantum-limited readout, starting in the ground state with the record at zero:

```wl
tr = oneTraj[{0.252, 0.1, 0.176, 1., 1., 3.}, {0., 0., 1., 0.}, 3., 1/120., 7];
Plot[{tr[[3]][t], tr[[4]][t]}, {t, 0, 3}, Frame -> True, PlotLegends -> {"z(t)", "Q(t)"}]
```

The monitored component and the accumulated readout move together because one noise drives both. Switch off drive and relaxation and the two are locked by a closed law, $z=\tanh(\operatorname{arctanh}z_0+\sqrt{\Gamma_{CI}}\,Q)$, which holds on the regularized path exactly by the ordinary chain rule:

```wl
With[{rG = {0.352, 0., 0.176, 0., 0., 0.}},
 With[{trG = oneTraj[rG, {0., 0., 0.6, 0.}, 3., 1/(40 0.352), 11], s = Subdivide[0., 3., 1000]},
  Max@Abs[trG[[3]][s] - Tanh[ArcTanh[0.6] + Sqrt[rG[[1]]] trG[[4]][s]]]]]
```

The residual is at the solver's tolerance, not the noise scale: the qubit's population is literally a function of what the detector has read so far.

## 4. Continuous Quantum Measurement, from Bayes' Rule to the Quantum Kalman Filter

Source: `/Users/mohammadb/Documents/GitHub/mbahram/continuous-measurement/Continuous Space/measurement-in-a-potential.md`

### Summary

A particle whose position is weakly and continuously measured, in any potential, simulated on a periodic grid. The stochastic Schrodinger equation is run as what it is in a laboratory: each step receives one noisy reading $\bar x$, updates the wavefunction by Bayes' rule (multiply by the square root of a Gaussian likelihood, renormalize), then evolves by a Strang split-operator step. The measurement half never sees the Hamiltonian, so any $V(x)$ costs one pointwise phase.

What should catch a physicist's attention:

- **Measurement is literally Bayes' rule.** The amplitude update and "likelihood times prior" agree to machine precision, and the essay keeps a strict ledger: only the record is data; the conditional mean, width and noise are inferences.
- **Free particle: localization and heating.** A broad packet narrows to the steady width $(\hbar/4\lambda m)^{1/4}$; the averaged state's momentum variance grows exactly as $\lambda\hbar^2t$, the uncertainty-principle price of position information; the leftover offset is shown to be first-order splitting bias by halving the step.
- **A cat collapses, and the observer can tell which way from the readings alone.** A running log-odds built from the noisy outcomes tracks the hidden lobe mass on the grid; different noise histories pick different lobes.
- **Harmonic trap: squeezing below the zero-point at any strength.** The conditional covariances obey a deterministic Riccati equation solved in closed form with $\lambda,\omega,m,\hbar$ symbolic. With $\mu=2\hbar\lambda/m\omega^2$, the steady variance in ground-state units is $\sqrt{2(\sqrt{1+\mu^2}-1)}/\mu$, proved below one and decreasing for every $\mu>0$ (`Resolve` over all $\mu$), with limits $1-\mu^2/8$ and $\sqrt{2/\mu}$. The grid, which never assumes a Gaussian, lands on the predicted pure squeezed Gaussian.
- **The quantum Kalman filter.** Five numbers updated from the readings alone reproduce the conditional mean that the 256-amplitude wavefunction carries: the filter levitated-nanoparticle and membrane experiments run.
- **The price of an imperfect detector, in closed form.** At efficiency $\eta$ the steady purity is exactly $\sqrt\eta$, independent of trap, mass and strength; squeezing below the zero-point needs $\mu\ge2\sqrt{1-\eta}/\eta^{3/2}$, so at $\eta=0.34$ and the essay's strength the state sits above the vacuum because the threshold is not reached, not because inefficiency forbids squeezing (cf. a levitated nanoparticle tracked to about 1.3 times its zero-point motion).
- **The averaged state never settles.** Energy rises exactly linearly at $\lambda\hbar^2/2m$ in any potential (the double commutator $[\hat x,[\hat x,\hat H]]$ is a constant); the heating lives entirely in the random wandering of the packet's center, whose Itô power equals the Lindblad rate exactly.

### Example 4.1: squeezing below the zero-point, and what detector efficiency costs, in closed form

The conditional covariances in a harmonic trap obey a Riccati flow that never sees the record. Keep everything symbolic and select the physical fixed point:

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

$\hbar^2/4$: the steady state is pure. In ground-state units the position variance is one function of one dimensionless strength $\mu$:

```wl
sqRatio[u_] := Sqrt[2 (Sqrt[1 + u^2] - 1)]/u;
{Simplify[SxSS/(\[HBar]/(2 m \[Omega])) == sqRatio[2 \[HBar] \[Lambda]/(m \[Omega]^2)], as],
 Resolve[ForAll[u, u > 0, sqRatio[u] < 1], Reals],
 Asymptotic[sqRatio[u], {u, 0, 2}, Assumptions -> u >= 0],
 Asymptotic[sqRatio[u], {u, Infinity, 1}, Assumptions -> u >= 0]}
```

The position is squeezed below the zero-point variance for every measurement strength, weakly as $1-\mu^2/8$ and strongly as $\sqrt{2/\mu}$. A detector catching a fraction $\eta$ of the light rescales $\lambda\to\eta\lambda$ and $\hbar\to\hbar/\sqrt\eta$ in the information terms only:

```wl
{SxE, CcE, SpE} = {SxSS, CcSS, SpSS} /. {\[Lambda] -> \[Eta] \[Lambda], \[HBar] -> \[HBar]/Sqrt[\[Eta]]};
{Simplify[\[HBar]/(2 Sqrt[SxE SpE - CcE^2]), {\[HBar] > 0, 0 < \[Eta] <= 1, \[Lambda] > 0, \[Omega] > 0, m > 0}],
 Simplify@Reduce[sqRatio[Sqrt[\[Eta]] u]/Sqrt[\[Eta]] <= 1 && u > 0 && 0 < \[Eta] < 1, u, Reals]}
```

The steady purity is exactly $\sqrt\eta$, whatever the trap, the particle or the strength; and going below the vacuum is still possible at any nonzero efficiency, but only past $\mu^*=2\sqrt{1-\eta}/\eta^{3/2}$, a threshold that diverges as the detector fails.

### Example 4.2: five numbers recover what 256 amplitudes compute

A grid, a Gaussian packet, the Bayesian measurement step and a Strang split-operator step ($\hbar=m=1$):

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

Watch a particle in a trap starting from the ground state, keeping both the full conditional wavefunction and the stream of readings:

```wl
wR = 1.; lamH = 2.; SxZP = hbar/(2 mass wR); SpZP = hbar mass wR/2;
gH = grid[256, 24.]; dtH = 0.002; ntH = 1500;
{states, {outs}} = BlockRandom[SeedRandom[21];
   Reap[trajectory[gH, 0.5 mass wR^2 #^2 &, lamH, dtH, ntH][gaussian[gH, 0., 0., Sqrt[SxZP]]]]];
xGrid = xMean[gH, #] & /@ states;
```

The quantum Kalman filter: two means nudged by the surprise in each reading, three covariances following the record-blind Riccati flow from Example 4.1:

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

The dashed filter, fed nothing but the noisy readings, sits on the conditional mean of the full wavefunction; their largest gap is a small fraction of the range the mean covers, the cost of a finite time step:

```wl
Max@Abs[mFilter - xGrid]/(Max@xGrid - Min@xGrid)
```

The wavefunction is never needed to track the particle: in a quadratic potential the conditional state is a Gaussian whose shape is fixed in advance and whose center follows the record.

## 5. Transmon cQED: Wolfram x Qolab Collaboration

Source: no authored `.md` exists. The text below was read from the machine conversion `/Users/mohammadb/Documents/GitHub/mbahram/QIS book/source/hardware-second-half/00-official-tutorial-notebooks/converted-md/TransmonLab-APS-Wolfram-Qolab.md`; the live source is `/Users/mohammadb/Documents/GitHub/mbahram/TransmonLab/Examples/TransmonLab-APS-Wolfram-Qolab.nb`. The examples load the TransmonLab paclet from `/Users/mohammadb/Documents/GitHub/mbahram/TransmonLab/TransmonLab/`.

Notes for review: (1) the converted `.md`, committed in the QIS-book repository, contains the Qolab access token in plain text in its `SystemCredential` line; it should be removed and the token rotated. (2) The tomography cell reports the linear-inversion Bloch vector, which here has length about 1.12 (TransmonLab itself warns it is unphysical); the maximum-likelihood estimate is physical and closer to the target, so Example 5.2 uses it.

### Summary

The demo shown at APS 2026 with Qolab: take a transmon from its circuit Hamiltonian to a prepared qubit state, and then all the way down to the digitizer, where the state exists only as clouds of single-shot (I, Q) points. The same analysis then runs on real IQ data from Qolab's hardware.

What should catch a physicist's attention:

- **From circuit to qubit.** The Cooper-pair-box Hamiltonian $4E_C(\hat n-n_g)^2-E_J\cos\hat\phi$ in the charge basis, with the spectrum flattening against offset charge as $E_J/E_C$ grows: why the transmon works.
- **Any state with two hardware operations.** A virtual Z (a zero-duration phase update of the drive oscillator, $\varphi_{\rm rot}=\varphi+\pi/2$) followed by one calibrated pulse reaches any $(\theta,\varphi)$; the pulse-area theorem gives the amplitude, a simulated Power Rabi sweep calibrates it, and the simulated Bloch vector lands on the analytic target.
- **Measurement is a cavity, not a projector.** The readout resonator rings up along two qubit-dependent spirals $\dot\alpha_{0,1}=-(\kappa/2\pm i\chi)\alpha_{0,1}-i\varepsilon f(t)/2$ whose separation sets the signal-to-noise; with these parameters the library flags that Purcell decay through the cavity dominates $T_1$.
- **The whole chain, both directions.** Stochastic single-shot IQ records (raw and integrated, time-resolved) for three tomography bases, then the inverse: Gaussian-mixture fitting, assignment errors, linear-inversion and maximum-likelihood tomography, posterior distributions of the Bloch components.
- **Real hardware in the same pipeline.** One function call (`QolabExperiment`) returns IQ data from a Qolab qubit, analyzed by the identical code as the simulation.
- **Decoherence and DRAG.** Purity and fidelity decay on the $T_2$ scale after preparation, and a DRAG pulse reduces leakage to $|2\rangle$ in a three-level (Duffing) simulation.

### Example 5.1: arbitrary state preparation, calibrated from the pulse area

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/mbahram/TransmonLab/TransmonLab/"];
Needs["TransmonLab`"];
sys = TransmonSystem["QubitFrequency" -> 5.10*10^9, "ResonatorFrequency" -> 6.00*10^9,
   "Anharmonicity" -> -220*10^6, "Coupling" -> 100*10^6, "CavityLinewidth" -> 1*10^6,
   "T1" -> 50*10^-6, "T2" -> 30*10^-6, "OmegaMax" -> 2 Pi*42*10^6, "EpsilonMax" -> 2 Pi*5*10^6, "NMax" -> 5];
gaussPulse = Pulse["Gaussian", 50]; readoutPulse = Pulse["FlatTop", 500, "RiseTime" -> 10];
```

On resonance only the pulse area matters, so the rotation angle per unit amplitude is $\Omega_{\max}\int{\rm env}\,dt$:

```wl
thetaPerA = sys["OmegaMax"] gaussPulse["EnvelopeArea"] 10^-9;
```

Target $\theta=\pi/3$, $\varphi=\pi/6$: one frame rotation, one scaled pulse, then readout:

```wl
{thT, phT} = {Pi/3, Pi/6}; aT = N[thT/thetaPerA];
seq = PulseSequence[FrameRotation["Qubit", phT + Pi/2], Drive[Amplitude[aT]@gaussPulse, "Qubit"],
   Sync["Qubit", "Resonator"], Readout[readoutPulse, "Resonator", 500]];
res = TransmonPhysics[sys, seq, "Model" -> "QubitOnly"];
bz = #[gaussPulse["Duration"] 10^-9] & /@ Values[res["BlochTrajectories"][[1]]];
{bz, N@{Sin[thT] Cos[phT], Sin[thT] Sin[phT], Cos[thT]}}
```

The simulated Bloch vector at the end of the pulse sits on the analytic target to the accuracy allowed by the truncated Gaussian envelope. The same calibration constant appears experimentally as the period of a Power Rabi sweep:

```wl
TransmonPhysics[sys, PulseSequence[Drive[gaussPulse, "Qubit"], Sync["Qubit", "Resonator"],
    Readout[readoutPulse, "Resonator", 500]], Sweep["Amplitude", Range[0, 1, 0.02]],
  "Model" -> "QubitOnly"]["Plot"]
```

### Example 5.2: from the digitizer back to the Bloch vector

The cavity pointer states, one spiral per qubit state:

```wl
dr = DispersiveReadout[sys, readoutPulse];
IQPlot[dr]
```

Prepare the same target, rotated into each tomography basis, and simulate 200 noisy single shots per basis:

```wl
bv = BlochVector[{Cos[thT/2], Exp[I phT] Sin[thT/2]}];
{thX, phX} = N@Rest@ToSphericalCoordinates[RotationMatrix[-Pi/2, {0, 1, 0}] . bv];
{thY, phY} = N@Rest@ToSphericalCoordinates[RotationMatrix[Pi/2, {1, 0, 0}] . bv];
iqt = IQTrajectory[sys, dr, "Bases" -> {"x", "y", "z"}, "NShots" -> 200,
   "InitialStates" -> <|"z" -> {1.0, thT, phT}, "x" -> {1.0, thX, phX}, "y" -> {1.0, thY, phY}|>,
   "Seed" -> 123];
iqa = IQAnalysis[iqt["MeasurementData"]];
```

Reconstruct the state by maximum likelihood, which keeps it inside the Bloch ball, and compare with the target (the `IQAnalysis` call above prints two messages: `IQAnalysis::nonphysical` for the linear-inversion estimate, which is why it is not used, and a `FindMinimum::lstol` from its internal fit, which does not affect the result):

```wl
target = N@{Sin[thT] Cos[phT], Sin[thT] Sin[phT], Cos[thT]};
mleB = Re[Tr[iqa["DensityMatrixMLE"] . #] & /@ (PauliMatrix /@ {1, 2, 3})];
{mleB, Norm[mleB - target], iqa["BlochVectorAround"]}
```

From nothing but clouds of (I, Q) points, the reconstructed Bloch vector lands near the prepared one, with posterior error bars of a few hundredths per component from 200 shots. Swapping `iqt["MeasurementData"]` for the table returned by `QolabExperiment` runs the identical analysis on real hardware.

## 6. The Gravitational Soliton: Ground State, Finite-Element Extrapolation and Virial Verification

Source: `/Users/mohammadb/Nextcloud/Mads/Quantum/SchrdingerNewton/SN1_GroundStateSoliton_v3_reviewed-MB.md`

Notes for review: the file still carries reviewer dialogue (the line beginning "what about $\Phi(0)$? Mads: ... Oliver: ..." and the sentence "It looks as if at ϕ(0) we have a Neumann zero BC" after the potential plot), and the closing section says "the small virial mismatch at $R=400$" while the computation uses $R=300$.

### Summary

Can a quantum particle be held together by its own gravity? In the Schrodinger-Newton equation the wavefunction moves in the Newtonian potential sourced by its own mass density $m|\psi|^2$, a mean-field (not quantum-gravity, not collapse) model. Quantum spreading pushes out, self-gravity pulls in, and the balance is a stationary soliton.

What should catch a physicist's attention:

- **The scale in one line.** Minimizing $\hbar^2/m\ell^2-Gm^2/\ell$ gives $\ell_\star=\hbar^2/Gm^3$: astronomically large for an electron, near $10^{-7}$ m for a $10^{10}$ amu nanoparticle, which is why this matters for macroscopic-superposition experiments.
- **Three dimensions to one, then a self-consistent loop.** With $u=\sqrt{4\pi}\,rf$ the problem becomes a radial Schrodinger eigenproblem in a potential that the eigenfunction itself sources; Newton's equation (`NDSolveValue`) and the Schrodinger eigenproblem (`NDEigensystem` with a shifted Arnoldi solver) on a graded finite-element mesh alternate until the energy stops changing.
- **The box is removed, not ignored.** Because $\Phi\sim-1/r$, imposing $\Phi(R)=0$ shifts the energy by $+1/R$. Fitting $E(R)=E_\infty+C/R$ over six boxes reproduces the Moroz-Penrose-Tod value $E_\infty\approx-0.163$ and the predicted coefficient $C\approx1$; subtracting $1/R$ collapses every box onto one line.
- **A check the solver never saw.** The virial identity $2T+W=0$ forces $T/E=-1/3$. The raw solution misses it by exactly the $1/R$ boundary shift, and removing that shift restores $-1/3$ to five digits.
- **Excited solitons.** The first three bound states have zero, one and two radial nodes, all sitting in their own $-1/r$-tailed potentials.

### Example 6.1: the self-consistent soliton and the infinite-box energy

The natural length scale first:

```wl
Solve[D[\[HBar]^2/(m l^2) - G m^2/l, l] == 0, l]
```

The radial Newton and Schrodinger operators, and one solver for each half of the loop:

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

Iterate from a normalized Gaussian seed until the eigenvalue settles:

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

The fit returns $E_\infty$ at the literature ground-state energy and $C$ at the predicted value one: the whole finite-box error is the missing $-1/r$ tail.

```wl
ListLinePlot[{{#R, #Energy} & /@ results, {#R, #Energy - 1/#R} & /@ results}, PlotMarkers -> Automatic,
 PlotLegends -> {"E(R)", "E(R) - 1/R"}, Frame -> True, FrameLabel -> {"R", "E"},
 Epilog -> {Red, Dashed, InfiniteLine[{0, e /. fit}, {1, 0}]}]
```

The six boxes take several minutes in total (about seven in one full-kernel run of this file).

### Example 6.2: the virial identity as an independent check

Kinetic and potential energy of the $R=300$ ground state:

```wl
sol = Last[results]; u0 = sol["u"]; Phi0 = sol["Phi"]; e0 = sol["Energy"];
t = 1/2 NIntegrate[D[u0[r], r]^2, {r, 0, 300}];
v = NIntegrate[Phi0[r] u0[r]^2, {r, 0, 300}];
{t + v - e0, t/e0, t/(e0 - 1/300), {-3 t - e0, -1/300.}}
```

$T+V$ reproduces the eigenvalue, so the energy bookkeeping is right. The raw ratio $T/E$ misses the virial value $-1/3$, and the miss $-3T-E$ equals the boundary shift $-1/R$; removing that shift returns $T/E=-1/3$. A relation the solver never imposed confirms the solution and locates its only error in the far-field boundary condition.

## 7. The Stabilizer Formalism by Computation: Aaronson-Gottesman Tableaux, Clifford Circuits and QEC Codes

Source: `/Users/mohammadb/Documents/GitHub/mbahram/stabilizer-symbolic/tutorial/Stabilizer-Formalism-By-Computation.md` (uses QuantumFramework's `PauliStabilizer`, loaded from the local checkout).

### Summary

A chalkboard course on the stabilizer formalism in which every claim is computed. A state is named by the Pauli operators that fix it; $n$ commuting generators replace $2^n$ amplitudes; Clifford gates become $O(n)$ bit rewrites of an Aaronson-Gottesman tableau; and the same object carries the reader from one qubit to a thousand, then on to codes, graph states, magic and noise.

What should catch a physicist's attention:

- **Why Paulis, precisely.** Every state is fixed by *some* operator (the reflection $2|\psi\rangle\langle\psi|-1$), so "a state is what fixes it" is empty on its own; the essay isolates the four properties of the Pauli group (cheap to name, $\pm1$ spectrum, commute-or-anticommute, closed under Clifford conjugation) that make it pay.
- **Schrodinger and Heisenberg meet.** The same circuit run by state-vector evolution and by conjugating generators gives the same state up to a global phase; a `PauliDecompose` shows $UZ_1U^\dagger$ stays one Pauli string, while $TXT^\dagger=(X+Y)/\sqrt2$ shows exactly where the method breaks.
- **Stabilizer states are flat.** Every stabilizer state is a uniform superposition over an affine subspace of $\mathbb F_2^n$, its support cut out by the $Z$-type checks.
- **Entanglement is a rank.** Entropy across a cut is an $\mathbb F_2$ rank (Fattal et al.), computed for a 1000-qubit GHZ state in a fraction of a second.
- **Measurement without branching.** Measurement histories grow as trees (`NestTree`); the SymPhase method of Fang and Ying instead carries each random outcome as a formal bit, so correlations come out as explicit Boolean relations and one symbolic traversal serves any number of shots.
- **Codes.** The five-qubit code's fifteen single-qubit errors map one-to-one onto the fifteen nonzero syndromes, which is what "perfect" means; Steane's code shows the CSS split.
- **Past the Clifford edge.** Stabilizer frames make the stabilizer rank visible ($2^t$ for $t$ T gates, independent of $n$), compute the irrational Born probability of the magic state exactly, compress redundant terms, and on diagonal circuits collapse to a rank-one phase polynomial whose amplitudes cost one linear solve at any $n$.
- **An honest benchmark.** Phase-resolved timings against Stim and QuantumClifford.jl, stating where QF is 4 to 7 times behind (kernel) and what it offers in exchange (symbolic signs, exact arithmetic, one object model).

### Example 7.1: the perfect code's syndrome table, and a thousand-qubit GHZ state

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
code = PauliStabilizer["5QubitCode"];
checks = Most[code["Stabilizers"]];
syndrome[err_] := (1 - (code[err]["Expectation", #] & /@ checks))/2;
errors = Flatten[Table[{p -> q}, {p, {"X", "Y", "Z"}}, {q, 1, 5}], 1];
AssociationMap[syndrome, errors]
```

Each of the fifteen single-qubit errors flips a different pattern of the four checks. The patterns are exactly the fifteen nonzero four-bit strings, none missing and none repeated:

```wl
Sort[syndrome /@ errors] === Sort@DeleteCases[Tuples[{0, 1}, 4], {0, 0, 0, 0}]
```

The same object scales without forming a state vector. A 1000-qubit GHZ state and its entanglement entropy across one qubit and across half the chain:

```wl
ghz = PauliStabilizer[1000][Prepend[Table["CNOT" -> {q - 1, q}, {q, 2, 1000}], "H"]];
{AbsoluteTiming[ghz["Entropy", {1}]], ghz["Entropy", Range[500]]}
```

One bit either way: GHZ entanglement is a single shared parity, however the chain is cut.

### Example 7.2: measurement outcomes as algebra, and the price of magic

Copy the parity of two random qubits onto a third, then measure all three symbolically. Only two fresh outcome bits appear, because the conserved $ZZZ$ fixes the third:

```wl
parityCheck = PauliStabilizer[3][{"H" -> 1, "H" -> 2, "CNOT" -> {1, 3}, "CNOT" -> {2, 3}}];
smPar = parityCheck["SymbolicMeasure", {1, 2, 3}];
smPar["Stabilizers"]
```

Substituting the free bits over all assignments reads off the third outcome as their XOR, and a thousand samples drawn by substitution alone never leave the even-parity sector:

```wl
readBits[ps_] := First@Keys@Select[AssociationThread[Tuples[{0, 1}, ps["Qubits"]],
      Normal[ps["State"]["StateVector"]]], # != 0 &];
free = DeleteDuplicates@Cases[smPar["Phase"], \[FormalS][_], Infinity];
{readBits[smPar["SubstituteOutcomes", Thread[free -> #]]] & /@ Tuples[{0, 1}, Length[free]],
 SeedRandom[7]; Counts[readBits /@ smPar["SampleOutcomes", 1000]]}
```

Cross the Clifford boundary with T gates. Three T gates on a GHZ backbone give a stabilizer frame of eight terms whatever the width, while the state vector it stands for grows as $2^n$; and a single T on $|+\rangle$ gives a Born probability no stabilizer state can:

```wl
magicState[n_] := PauliStabilizer[n][Join[{"H" -> 1}, Table["CNOT" -> {i, i + 1}, {i, 1, n - 1}],
     {"T" -> 1, "T" -> 2, "T" -> 3}]];
{Table[{n, magicState[n]["Length"], 2^n}, {n, {4, 10, 20}}],
 FullSimplify[(1 + PauliStabilizer[1][{"H", "T"}]["Expectation", "X"])/2]}
```

The probability $(2+\sqrt2)/4$ comes straight from the two-term frame, never from a state vector: the cost of magic is set by the number of T gates, not by the number of qubits.

## 8. Dirac Matrices in a Curved Spacetime

Source: no talk essay exists. The original is the notebook `/Users/mohammadb/Nextcloud/Mads/Quantum/DiracMatrix/gamma-dirac-matrix-v3.nb`, read here through its machine conversion `/Users/mohammadb/Nextcloud/Mads/Quantum/Roadmap and planning/nb-converted/3432d4a045f86e4f6487314c8e8aed64.md`; the same functions now live in the DiracMatrix paclet, documented in `/Users/mohammadb/Nextcloud/Mads/Quantum/DiracMatrix/DiracMatrixPackage/DiracMatrix-Documentation.md`, which the examples load.

Notes for review: (1) **`CliffordBasis` is wrong on non-flat metrics, and the notebook and documentation blame the wrong function.** Both say the explicit antisymmetrized sum (`CliffordCanonicalBasis`) drifts at high grade "due to numerical issues" while the recursion (`CliffordBasis`) is stable. With a moderate non-diagonal Lorentzian metric in machine precision, the explicit sum reproduces the antisymmetrized products to machine precision while `CliffordBasis` is off by order ten already at grade 3. The cause is the sign in the recursion at `/Users/mohammadb/Nextcloud/Mads/Quantum/DiracMatrix/DiracMatrixPackage/Kernel/DiracMatrix.wl:186`: `(-1)^(j - 1)` should be `(-1)^(k - j)` with `k = Length[is]`. With that sign the recursion matches to machine precision. Flat metrics hide the bug because every off-diagonal $\eta_{\mu i}$ vanishes. (2) "Curved spacetime" here means the gamma matrices of an arbitrary metric at a point, through a vielbein; spin connection and the curved-space Dirac operator are not part of the material. (3) Loading with `Needs["DiracMatrix`"]`, as the documentation says, emits `Needs::nocont` (the paclet's primary context is ``WolframQuantumComputation`DiracMatrix` ``) even though the functions load.

### Summary

Gamma matrices satisfying $\{\gamma^\mu,\gamma^\nu\}=2g^{\mu\nu}\mathbb 1$ for any real symmetric metric $g$, of any signature $(p,q)$ and any dimension. Flat-space matrices come from the Brauer-Weyl (Pauli-Kronecker) construction, $\Gamma^{2k-1}=\mathbb 1^{\otimes(k-1)}\otimes\sigma_1\otimes\sigma_3^{\otimes(m-k)}$, $\Gamma^{2k}=\mathbb 1^{\otimes(k-1)}\otimes\sigma_2\otimes\sigma_3^{\otimes(m-k)}$; time-like directions pick up a factor $i$; a general metric is handled by a vielbein, $g=e^T\eta e$, $\gamma^\mu=e^\mu{}_a\Gamma^a$.

What should catch a physicist's attention:

- **One function, every signature and dimension.** `GammaMatrices[g]` returns $n$ matrices of size $2^{\lfloor n/2\rfloor}$ for Euclidean, Minkowski, split $(2,2)$, Schwarzschild at a point, an FLRW slice, or a random non-diagonal metric.
- **The textbook identities as tests, swept over all signatures.** Tracelessness, $\operatorname{Tr}\gamma^\mu\gamma^\nu=2^{\lfloor n/2\rfloor}g^{\mu\nu}$, the four-gamma trace, $\gamma^\mu\gamma_\mu=n$, $\gamma^\mu\gamma_\nu\gamma_\mu=(2-n)\gamma_\nu$, the Lorentz-generator commutator, chirality in even dimension: each is one `Table` over $n\le10$ and all $(p,q)$, flat and curved.
- **Representations on demand.** Brauer-Weyl, Dirac ($\gamma^0$ diagonal) and Weyl/chiral bases by explicit similarity transforms, with the odd-dimensional case (no chiral basis) handled.
- **The graded Clifford basis.** All $2^n$ antisymmetrized products, the operator basis behind Fierz identities, by explicit sum and by a grade-raising recursion (see the review note on the recursion).

### Example 8.1: gamma matrices for Schwarzschild and for a random curved metric

```wl
PacletDirectoryLoad["/Users/mohammadb/Nextcloud/Mads/Quantum/DiracMatrix/DiracMatrixPackage"];
Needs["DiracMatrix`"];
```

The Schwarzschild metric on the equator at $r=3M$ is $\mathrm{diag}(-1/3,3,9,9)$. Its gamma matrices obey the Clifford relation exactly:

```wl
g = DiagonalMatrix[{-1/3, 3, 9, 9}];
Gs = GammaMatrices[g];
{Outer[#1 . #2 + #2 . #1 &, Gs, Gs, 1] == 2 TensorProduct[g, IdentityMatrix[4]], MatrixForm[Gs[[1]]]}
```

A random non-diagonal metric of Lorentzian signature, its vielbein, and the check that the vielbein rebuilds the metric:

```wl
SeedRandom[42]; gr = RandomCurvedMetric[1, 3, 1, 1];
{sig, e} = MetricVielbein[gr];
{Diagonal[sig], NumericZeroQ[Transpose[e] . sig . e - gr]}
```

The signature $(+,-,-,-)$ is read off the metric itself, and the vielbein is exact to numerical tolerance.

### Example 8.2: the textbook identities, on random curved metrics of every signature

For every dimension $n=2,\ldots,6$ and every signature $(p,n-p)$, draw a random non-diagonal metric and test the Clifford relation, the contraction $\gamma^\mu\gamma_\mu=n$ and the four-gamma trace:

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

`True`: the identities that every QFT course derives by hand for Minkowski space hold, as computed facts, for every signature and every metric in the sweep. Representations follow by one call, for example the Weyl basis in Minkowski space:

```wl
MatrixForm /@ ToWeylBasis[GammaMatrices[FlatMetric[1, 3]]]
```

## 9. Snyder's Lorentz-Covariant Quantum Spacetime

Source: no `.md` exists. The notebook `/Users/mohammadb/Nextcloud/Mads/Quantum/2nd quantization and non-commutative/NCA examples/NCA talks/Snyder's Lorentz-covariant quantum spacetime.nb` was read through an `nb-reader` conversion; a shorter version (triple-commutator check only) is section P17 of `/Users/mohammadb/Nextcloud/Mads/Quantum/2nd quantization and non-commutative/NCA examples/NCA talks/Overll nb and md/NCA-talks-4-17.md`. The examples use the built-in `NonCommutativeAlgebra` and `NonCommutativeExpand` (Wolfram Language 15).

### Summary

Snyder (1947) quantized spacetime without breaking Lorentz symmetry: the coordinates stop commuting, and their commutator is a Lorentz generator times a fundamental length squared, $[x^\mu,x^\nu]=ia^2M^{\mu\nu}$. The notebook declares this ten-generator algebra once and lets the noncommutative-algebra engine do what is error-prone by hand.

What should catch a physicist's attention:

- **Consistency first.** All four families of Jacobi identities vanish, and Lorentz covariance of $[x,x]$ comes out as a theorem: $[M^{\rho\sigma},[x^\mu,x^\nu]]$ computed through the $[M,x]$ relations equals $ia^2[M^{\rho\sigma},M^{\mu\nu}]$ computed through the Lorentz algebra.
- **Orderings no one does by hand.** The quartic word $x^3x^2x^1x^0$ reduces to the ordered basis with 136 leaves of $x$ and $a^2xM$ terms; the cubic Weyl symbol follows the same way.
- **The deformed center.** $x^2$ is a Lorentz scalar but not central, and the Lorentz Casimir is not central either; the SO(4,1) combination $x^2-\tfrac{a^2}{2}M^{\mu\nu}M_{\mu\nu}$ commutes with everything (Maggiore's conjecture, checked).
- **Noncommutativity is orbital angular momentum.** With momenta and $[x^\mu,p^\nu]=i(\eta^{\mu\nu}+a^2p^\mu p^\nu)$, Jacobi closure forces $M^{\mu\nu}=x^\mu p^\nu-x^\nu p^\mu$: the generator in $[x,x]$ is the orbital one already present in phase space.
- **Physics in one stroke.** The BCH expansion of a product of plane waves gives the deformed momentum addition (symmetric for Snyder at this order, non-symmetric already at order $1/\kappa$ for κ-Minkowski, run through the same code); the deformed uncertainty relation has a minimum $\Delta x_{\min}=a$; and the Fermi-LAT GRB 090510 bound translates to $a<1.5\times10^{-27}$ m, about eight decades above the Planck length.

### Example 9.1: Snyder's algebra, its consistency, and its Casimir

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

Jacobi identities for three coordinates and for two coordinates with a Lorentz generator, over every index choice:

```wl
jac[X_, Y_, Z_] := Commutator[X, Commutator[Y, Z]] + Commutator[Y, Commutator[Z, X]] +
   Commutator[Z, Commutator[X, Y]];
{Union@Flatten@Table[NonCommutativeExpand[jac[x[m], x[n], x[r]], alg], {m, 0, 3}, {n, 0, 3}, {r, 0, 3}],
 Union@Flatten@Table[NonCommutativeExpand[jac[x[m], x[n], MGen[r, s]], alg],
   {m, 0, 3}, {n, 0, 3}, {r, 0, 3}, {s, 0, 3}]}
```

Lorentz covariance of $[x,x]$, derived rather than assumed:

```wl
Union@Flatten@Table[NonCommutativeExpand[Commutator[MGen[r, s], Commutator[x[m], x[n]]], alg] -
    I a^2 I (eta[s, m] MGen[r, n] - eta[s, n] MGen[r, m] - eta[r, m] MGen[s, n] + eta[r, n] MGen[s, m]),
  {r, 0, 3}, {s, r + 1, 3}, {m, 0, 3}, {n, m + 1, 3}]
```

The Snyder Casimir $x^2-\tfrac{a^2}{2}M^{\mu\nu}M_{\mu\nu}$ commutes with every coordinate and every Lorentz generator, while $x^2$ alone does not commute with $x$:

```wl
C1 = Sum[eta[m, m] eta[n, n] MGen[m, n] ** MGen[m, n], {m, 0, 3}, {n, m + 1, 3}];
xSq = Sum[eta[m, m] x[m] ** x[m], {m, 0, 3}];
{Union@Table[NonCommutativeExpand[Commutator[xSq - a^2 C1, x[r]], alg], {r, 0, 3}],
 Union@Flatten@Table[NonCommutativeExpand[Commutator[xSq - a^2 C1, MGen[r, s]], alg], {r, 0, 3}, {s, r + 1, 3}],
 NonCommutativeExpand[Commutator[xSq, x[0]], alg]}
```

Every list collapses to `{0}`; the last entry is the nonzero $a^2(xM+Mx)$ residue that the $-a^2C_1$ term exists to cancel.

### Example 9.2: deformed momentum addition, Snyder against κ-Minkowski, and the minimal length

Treat the momenta $\alpha,\beta$ as commuting parameters and expand $\log(e^{i\alpha\cdot x}e^{i\beta\cdot x})$ to fourth order in the BCH series:

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

Snyder's sum $\alpha\oplus\beta$ picks up $a^2$ corrections but stays symmetric under $\alpha\leftrightarrow\beta$. The same driver on κ-Minkowski, $[x^0,x^i]=ix^i/\kappa$:

```wl
algK = NonCommutativeAlgebra[<|"Generators" -> {Table[x[m], {m, 0, 3}]},
    "CommutationRelations" -> Join[Table[Commutator[x[i], x[0]] - I/k x[i], {i, 1, 3}],
      Flatten@Table[Commutator[x[j], x[i]], {i, 1, 3}, {j, i + 1, 3}]],
    "ScalarVariables" -> Join[{k}, as, bs]|>];
Simplify@Table[Coefficient[NonCommutativeExpand[bch4[I ax, I bx], algK], x[m]]/I, {m, 0, 1}]
```

Here time momenta add undeformed but spatial momenta pick up $(\alpha^0\beta^i-\alpha^i\beta^0)/2\kappa$, which flips sign under exchange: two Lorentz-compatible deformations with different multi-particle physics. Finally the deformed uncertainty relation $\Delta x\ge(1+a^2\Delta p^2)/2\Delta p$ has a floor:

```wl
Assuming[a > 0, Minimize[{(1 + a^2 dp^2)/(2 dp), dp > 0}, dp] // FullSimplify]
```

$\Delta x_{\min}=a$ at $\Delta p=1/a$: the algebra's parameter is a length below which no event can be localized.

## 10. Breaking Elliptic Curve Crypto: Quantum Oracles and Phase Estimation for secp256k1

Source: no `.md` for the computation exists. The notebook `/Users/mohammadb/Nextcloud/Mads/Quantum/Shor-Google-Security/ShorECDLP-Google2026-v4.nb` (2026-04-30, the newest of v1 to v4) was read through an `nb-reader` conversion. The markdown with the matching title, `/Users/mohammadb/Nextcloud/Mads/Quantum/Shor-Google-Security/google-quantum-crypto-report.md`, is a guide to Google's paper rather than this computation and was not used for the examples.

Notes for review: **the notebook's quantum step does not solve the discrete logarithm.** Its phase-estimation circuit uses only the oracle $|i\rangle\to|i+G\rangle$, so it recovers the order of $G$; the secret $k$ is then found by a classical scan over $kG$. But the order of the secp256k1 generator is published in the standard, so this step gives an attacker nothing, and the scan is exactly the hard problem. The summary's "Phase estimation on 15 qubits recovered the secret key $k=7$" is therefore not what the code does, and the "Pattern That Carries Over" section's claim that only the group changes between factoring and ECDLP needs the two-register correction. Shor's ECDLP algorithm queries $f(a,b)=aG+bQ$ on two registers; Fourier sampling returns pairs $(c,d)\approx(j,jk)\cdot2^t/r$, and $k=l\,j^{-1}\bmod r$ follows without ever scanning. Example 10.1 implements that version and recovers $k$ directly. The resource-estimate half of the notebook (Google 2026, Oratomic/Caltech, on-spend race) is unaffected; the numbers checked below match.

### Summary

How close is a quantum computer to stealing a Bitcoin key? The notebook walks Shor's algorithm from factoring 15 to the elliptic-curve discrete logarithm on a toy curve over GF(23), then reads Google's March 2026 resource estimates for secp256k1 and the Oratomic/Caltech neutral-atom compilation.

What should catch a physicist's attention:

- **The group is the only thing that changes.** Elliptic-curve point addition, defined from scratch with pattern-matched dispatch, gives a permutation unitary exactly like modular multiplication in factoring; its eigenphases are the 28th roots of unity of the cyclic group of order 28.
- **Real cryptography in the same language.** `GenerateAsymmetricKeyPair["EllipticCurve"]` uses secp256k1 by default, and a signature over the Bitcoin genesis-block headline verifies, so the object under attack is the real one.
- **Precision is a theorem, not a tuning knob.** Continued fractions recover $r$ only when $2^t\ge2r^2$ (Legendre); in the notebook's exercise, 6 control qubits fail in all 50 samples while 10 succeed because favored shots land on the right convergents.
- **What Google actually optimized.** Not the algorithm but the reversible point-addition circuit: windowed arithmetic and measurement-based uncomputation. The appendix formula $T=(T_{PA}+3\cdot2^w)(2n/w-4)$ reproduces about 81 million Toffolis against the reported 90 million at most.
- **Two hardware roads.** Under half a million superconducting qubits for minutes-long attacks, or about 10,000 to 26,000 neutral atoms with high-rate qLDPC codes for days-long ones; the first threatens transactions in flight (about 41% success against a Bitcoin block at a 9-minute attack), the second dormant wallets.

### Example 10.1: Shor's discrete-log algorithm on a toy curve, recovering the key from the public point alone

Elliptic-curve arithmetic over GF(p) from the notebook, with `O` the point at infinity:

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

The curve $y^2=x^3+x+1$ over GF(23), a generator $G$, its public order $r$, and a public key $Q=7G$:

```wl
ec = EllipticCurve[1, 1, 23];
pts = ec["Points"]; r = Length[pts] + 1;
G = SelectFirst[pts, Length@DeleteDuplicates@Table[ec["Mult"][#, k], {k, 0, r - 1}] == r &];
Q = ec["Mult"][G, 7];
{r, G, Q}
```

The quantum part queries $f(a,b)=aG+bQ$ on two $t$-qubit registers, measures the point register, and Fourier-transforms both input registers. Simulate it exactly: one indicator array per curve point, a 2D Fourier transform of each, probabilities summed over points:

```wl
elts = Prepend[pts, O]; idx = AssociationThread[elts -> Range[0, r - 1]];
t = 7; nn = 2^t;
f = Table[idx[ec["Add"][ec["Mult"][G, a], ec["Mult"][Q, b]]], {a, 0, nn - 1}, {b, 0, nn - 1}];
probs = Total[Table[Abs[Fourier[Boole[Map[# == w &, f, {2}]], FourierParameters -> {0, -1}]]^2,
      {w, 0, r - 1}]]/nn^2;
```

Each measured pair $(c,d)$ sits near $(j,jk)\cdot2^t/r$, so whenever $j$ is invertible mod $r$ the key is $l\,j^{-1}\bmod r$. Twenty shots:

```wl
SeedRandom[1];
samples = RandomChoice[Flatten[probs] -> Tuples[Range[0, nn - 1], 2], 20];
keyFrom[{c_, d_}] := With[{j = Mod[Round[r c/nn], r], l = Mod[Round[r d/nn], r]},
   If[CoprimeQ[j, r], Mod[l PowerMod[j, -1, r], r], Missing[]]];
cands = keyFrom /@ samples;
{cands, Select[DeleteMissing@DeleteDuplicates[cands], ec["Mult"][G, #] === Q &]}
```

Every usable shot returns 7 and none returns a wrong key: the secret comes out of the interference pattern, using $Q$ only as a point and never scanning multiples of $G$. That is the step whose cost Google's circuits bring down to about a thousand logical qubits at 256 bits.

### Example 10.2: the resource arithmetic behind the headline numbers

Google's appendix formula for the total Toffoli count at window size $w$, with the low-qubit point-addition cost, and the window that minimizes it at fixed $T_{PA}$:

```wl
T[w_, n_, tpa_] := (tpa + 3 2^w) (2 n/w - 4);
{T[16, 256, 2700000], First@MinimalBy[Range[4, 24], T[#, 256, 2700000] &]}
```

About 81 million Toffolis at $w=16$, against the reported bound of 90 million, and the optimum window sits at 17, where the lookup-table cost $3\cdot2^w$ starts to rival the point addition. Blocks arrive as a Poisson process, so the chance an on-spend attack of 9 minutes beats the next block is $e^{-9/T_{\rm block}}$ for Bitcoin, Litecoin, Zcash and Dogecoin:

```wl
N[Exp[-9/#]] & /@ {10, 2.5, 1.25, 1}
```

About 41% against Bitcoin and negligible against the fast chains.

## 11. What Makes Wolfram QuantumFramework Unique: A Verified Showcase

Source: `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Why QF/showcase/QF-Showcase.md` (examples run against the local QuantumFramework checkout).

### Summary

Seven claims, each argued by executed code: QF is symbolic-first; one object model spans channels, POVMs, Fock space and phase space; measurement records are quantum wires; three engines (state vector, tensor network, stabilizer tableau) compute the same objects; a QF object is an ordinary Wolfram Language expression; every object draws itself; and the same objects leave for OpenQASM, transpilers and IBM hardware. Each claim names its closest competitor and what remains distinctive.

What should catch a physicist's attention:

- **Time-ordered propagators in closed form.** A rotating-field qubit returns the lab-frame Rabi propagator; a Landau-Zener sweep returns a propagator in Kummer and parabolic-cylinder functions whose plot shows the Landau-Zener-Stückelberg oscillations the asymptotic formula $e^{-\pi\Delta^2/2v}$ misses.
- **Exact variational landscapes.** QAOA on $K_4$, built entirely from the graph, gives a closed-form energy surface whose maximum is a certified cubic root, below the true cut of 4.
- **One master-equation solution for every qubit.** A symbolic density matrix and symbolic rates give the complete relaxation law at once, and $T_2\le2T_1$ comes out as a theorem via `Reduce`.
- **Breadth in one object model.** Amplitude-damping purity loss $2\gamma(1-\gamma)\sin^4(\theta/2)$, a tetrahedral SIC-POVM on a symbolic state, coherent-state dynamics as an identity in $\alpha,\omega,t$, $g^{(2)}$ for Fock, coherent and thermal light, Wigner versus Husimi for a cat, exact zeros of multiphoton interference in Fourier multiports, a spin-1 NV Hamiltonian and a five-site ring diagonalized by named bases, and tomography from counts.
- **Measurement without collapse.** Teleportation with its classical channel as record wires; a coherent $R_X(2\theta)$ error digitized into $\cos^2\theta$/$\sin^2\theta$ syndromes and corrected at fidelity exactly 1 for every input and angle; a measurement undone by its dilation's adjoint.
- **Three engines, exact interconversion.** A Grover search contracted along a steered tensor-network path to an exact rational success probability; a 1000-qubit random Clifford state with an exact half-chain entropy; tableau, state vector and synthesized circuit agreeing under `SameQ`; and a code layer returning a code's logical error rate as an exact polynomial in $p$.
- **Language-level operators.** `Exp`, `D`, `Plus`, `Commutator` act on quantum objects; $H(t)=i\dot UU^\dagger$ recovers the drive of a target trajectory; a Werner-state criterion passed to `Reduce` gives the exact entanglement threshold $p<1/2$.
- **Out to hardware.** OpenQASM 3 export and import with exact angles (an imported CHSH circuit gives exactly the Tsirelson value), device transpilation, and a 4096-shot GHZ run on `ibm_fez` drawn next to the exact distribution.
- **An honest boundary.** Closed forms stop where spectra stop being solvable in radicals (Abel-Ruffini), and QF signals it by returning `Root` objects.

### Example 11.1: the master equation solved once, for every qubit and every rate

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

Every entry of the initial state stays a symbol: populations relax at $\gamma_1$, coherences precess at $\Omega$ and decay at $\Gamma_2=\gamma_1/2+2\gamma_\phi$. The textbook bound between the two times is then a statement about all admissible rates:

```wl
Reduce[ForAll[{\[Gamma]1, \[Gamma]\[Phi]}, \[Gamma]1 > 0 && \[Gamma]\[Phi] >= 0,
  1/(\[Gamma]1/2 + 2 \[Gamma]\[Phi]) <= 2/\[Gamma]1]]
```

`True`: $T_2\le2T_1$ is derived from the solution, not quoted, with equality exactly when pure dephasing vanishes.

### Example 11.2: a coherent error, digitized and corrected, with measurement records as wires

The parity measurements are degenerate projective measurements with two rank-two projectors; their outcomes land on record wires 0 and -1, and the corrections are gates controlled (and anti-controlled) on those records:

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

The joint distribution of the two syndrome records, for a symbolic data qubit and a symbolic error angle:

```wl
as = {\[Theta] \[Element] Reals, \[Alpha] \[Element] Reals, \[Beta] \[Element] Reals, \[Alpha]^2 + \[Beta]^2 == 1};
Assuming[as, FullSimplify[ComplexExpand[
   Diagonal[Normal[QuantumPartialTrace[coded, {3, 4, 5}]["DensityMatrix"]]]]]]
```

$\{\cos^2\theta,0,0,\sin^2\theta\}$: the continuous rotation has become "no error" or "qubit 2 flipped", independent of $\alpha,\beta$. After the feedforward, the data qubits match the ideal codeword:

```wl
Assuming[as, FullSimplify[ComplexExpand[QuantumSimilarity[QuantumPartialTrace[coded, {1, 2}],
    enc[QuantumTensorProduct[QuantumState[{\[Alpha], \[Beta]}], QuantumState["00"]]], "Fidelity"]]]]
```

Fidelity exactly 1 for every input and every angle: the digitization of errors, the founding theorem of quantum error correction, obtained from one symbolic run.

## 12. DRAG: Fast Gates on a Transmon

Source: `/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/Documentation/English/Tutorials/TransmonGates.md` (the paclet tutorial, byte-identical to `TransmonGatesV3.md`).

### Summary

How short can a π pulse on a transmon be? A transmon is a weakly anharmonic oscillator, so a short, strong pulse leaks population into $|2\rangle$ and shifts the qubit transition away from the drive, while a long pulse suffers $T_1$ and $T_2$. The tutorial writes the driven anharmonic oscillator with a DRAG term in the rotating frame, $H=\tfrac\alpha2b^{\dagger2}b^2+\delta b^\dagger b+\tfrac{\Omega}{2}(b+b^\dagger)-\tfrac{\lambda\dot\Omega}{2\alpha}i(b^\dagger-b)$, keeps five levels explicitly, and calibrates amplitude, detuning and DRAG coefficient together at every length under the Lindblad equation.

What should catch a physicist's attention:

- **Leakage seen directly.** An uncalibrated 8 ns π pulse at $\alpha/2\pi=-200$ MHz leaves about a tenth of the population in $|2\rangle$ and above; at 35 ns leakage is negligible but the flip is still incomplete because the drive shifts the qubit transition.
- **DRAG alone is not the answer.** Scanning $\lambda$ at fixed amplitude, leakage bottoms out near the weak-drive value $\lambda\approx1$, but the transfer error is smallest at a much smaller $\lambda$: the leftover error is population still in $|0\rangle$. Amplitude, detuning and $\lambda$ must be calibrated together.
- **The speed limit is set by decoherence, not leakage.** Calibrated under $T_1=T_2=100\,\mu$s, the DRAG curve turns up at a much shorter length and reaches a lower minimum than the plain pulse; at the optimum both errors are dominated by decay and dephasing, and DRAG wins by being shorter. A narrow dip near 10 ns echoes Da Costa Jesus et al. (2025).
- **The model checks itself.** The gate error averaged over six input states stays close to the transfer error (time-reversal symmetry of the pulse), and three levels underestimate leakage while five and six agree.
- **All in one object.** `AnnihilationOperator[5]`, a parametrized `QuantumOperator`, `QuantumEvolve` for both Schrödinger and Lindblad dynamics, and `"ReturnEquations"` to hand the 25 density-matrix equations to `NDSolveValue` and `NMinimize`.

### Example 12.1: leakage, and what DRAG does to it

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

The populations of $|0\rangle,|1\rangle,|2\rangle$ after an uncalibrated π pulse of 8 ns and of 35 ns:

```wl
{finalPopulations[5, {1, 0, 0, 8}][[;; 3]], finalPopulations[5, {1, 0, 0, 35}][[;; 3]]}
```

At 8 ns roughly a tenth of the population has left the qubit; at 35 ns leakage is gone but a fraction of a percent stays in $|0\rangle$. Now scan the DRAG coefficient at 8 ns and compare leakage with transfer error:

```wl
ListLogPlot[Transpose@Table[With[{p = finalPopulations[5, {1, drag, 0, 8}]},
    {{drag, Total[p[[3 ;;]]]}, {drag, 1 - p[[2]]}}], {drag, 0, 2, 0.1}],
 Joined -> True, PlotMarkers -> Automatic, PlotLegends -> {"leakage", "transfer error"},
 Frame -> True, FrameLabel -> {"\[Lambda]", "population"}]
```

Leakage is minimized near $\lambda\approx1$, the weak-drive value, but the transfer error is minimized well below it: suppressing leakage and completing the flip are different targets. And the truncation is safe:

```wl
Table[{n, Total[finalPopulations[n, {1, 0, 0, 8}][[3 ;;]]]}, {n, 3, 6}]
```

Three levels miss about half the leakage; five and six agree.

### Example 12.2: calibrating the 8 ns pulse, with and without DRAG

Minimize the transfer error over amplitude and detuning, and over the DRAG coefficient as well (closed system, a quick differential-evolution search):

```wl
err[a_?NumericQ, l_?NumericQ, d_?NumericQ] := 1 - finalPopulations[5, {a, l, d, 8}][[2]];
{NMinimize[{err[a, 0, d], 0 < a < 2 && -0.5 <= d <= 0.5}, {{a, 0.8, 1.2}, {d, -0.1, 0.1}},
   Method -> {"DifferentialEvolution", "SearchPoints" -> 10, "RandomSeed" -> 0}],
 NMinimize[{err[a, l, d], 0 < a < 2 && 0 <= l <= 3 && -0.5 <= d <= 0.5}, {{a, 0.8, 1.2}, {l, 0, 2}, {d, -0.1, 0.1}},
   Method -> {"DifferentialEvolution", "SearchPoints" -> 10, "RandomSeed" -> 0}]}
```

This search takes about five minutes, and even so it is coarse; it still shows the mechanism: with only amplitude and detuning the 8 ns pulse cannot get below the percent level, while adding the DRAG coefficient cuts the error several-fold, with a nonzero detuning compensating the drive-induced shift. The tutorial's full scan repeats this under the Lindblad equation from 6 to 40 ns (the longest computation in the tutorial, run in parallel) to find where decoherence takes over.

## 13. QuantumInterferometer: Symbolic Behaviour

Source: `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Quantum Interferometer/QI-symbolic-showcase.md` (a reference page of input/output pairs rather than an essay; `QuantumInterferometer` lives in the ``Wolfram`QuantumFramework`SecondQuantization` `` context, which must be loaded).

### Summary

`QuantumInterferometer` is QF's linear-optics object: beam splitters, phase shifters, squeezers and loss acting on Fock, coherent and squeezed light. This page catalogues what it does when the inputs are symbols, angles, amplitudes, squeezing parameters, loss rates, and above all photon numbers.

What should catch a physicist's attention:

- **Symbolic photon numbers.** A beam splitter's amplitude from $|n_1,n_2\rangle$ to $|t_1,t_2\rangle$ comes back as a Wigner $D$-function times a Kronecker delta for photon conservation (the Schwinger two-mode realization of SU(2)), for symbolic $n_1,n_2,t_1,t_2$. Chains of splitters return held sums over internal photon numbers with the right bounds, and path conditions appear as `Piecewise` only when photon flow does not imply them.
- **Hong-Ou-Mandel for every $n$.** On a balanced splitter $|n,n\rangle\to|k,2n-k\rangle$ is $d^n_{0,k-n}(\pi/2)$, and every odd split vanishes.
- **Coherent, squeezed and lossy light, exactly.** Coherent states map to coherent states $|\alpha\rangle\to|W\alpha\rangle$; squeezed-vacuum photon statistics come out as $\operatorname{sech}r\,\tanh^{2m}r$ series; two squeezers on a splitter cancel coincidences for equal squeezing and maximize them for opposite squeezing (a two-mode squeezed vacuum); loss gives the binomial law $\binom nk\eta^k(1-\eta)^{n-k}$ for symbolic $n$, multiplicative rates on one mode, and coherences scaled by $\sqrt\eta$.
- **A physical-domain contract.** Photon numbers are nonnegative integers, angles real, $0\le\eta\le1$; results carry no conditions restating that domain, user assumptions refine them, and contradictory assumptions fail rather than return nonsense.
- **Stated limits.** Section 8 lists precisely what is not supported symbolically, and the final section records the gaps found while running the list (substitute photon numbers before machine-number angles, which `Refine` simplifications are left to the caller).

### Example 13.1: a beam splitter with symbolic photon numbers, and Hong-Ou-Mandel for every n

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
bs = QuantumInterferometer[{"BS"[th, ph]}];
bs["Amplitude", {n1, n2} -> {t1, n1 + n2 - t1}]
```

The amplitude is $D^{(n_1+n_2)/2}_{(n_1-n_2)/2,\;t_1-(n_1+n_2)/2}(\phi,2\theta,-\phi)$: a two-mode beam splitter is a rotation of a spin of size $(n_1+n_2)/2$, and the photon numbers are its magnetic quantum numbers. On a balanced splitter with $n$ photons in each port:

```wl
hom = QuantumInterferometer[{"BS"[Pi/4, 0]}]["Amplitude", {n, n} -> {k, 2 n - k}];
{hom, Table[hom /. {n -> m, k -> j}, {m, 1, 3}, {j, 1, 2 m - 1, 2}], Table[Abs[hom /. {n -> 2, k -> j}]^2, {j, 0, 4}]}
```

One formula for all $n$; every odd split is exactly zero (the $n=1$ case is the Hong-Ou-Mandel dip), and for two photons per port the outputs are $\{3/8,0,1/4,0,3/8\}$. Chaining two splitters gives a held sum over the internal photon number that, once activated, agrees with the concrete query:

```wl
chain = QuantumInterferometer[{"BS"[th1, ph1], "BS"[th2, ph2]}];
FullSimplify[Activate[chain["Amplitude", {n, 0} -> {t, n - t}] /. {n -> 2, t -> 1}] -
  chain["Amplitude", {2, 0} -> {1, 1}]]
```

### Example 13.2: loss, squeezing and coherent light in closed form

Loss on $n$ photons, with $n$ symbolic, and a coherent pair through a general splitter:

```wl
{QuantumInterferometer[{"Loss"[eta]}]["Probability", {n} -> {k}], bs["Coherent"[{a, b}]]}
```

The binomial loss law, and a coherent state that stays coherent with amplitudes rotated by the transfer matrix. Two squeezed vacua meeting on a balanced splitter, and the photon statistics of one squeezed vacuum:

```wl
{QuantumInterferometer[{"S"[x1] -> 1, "S"[x2] -> 2, "BS"[Pi/4, 0]}]["Amplitude", {0, 0} -> {1, 1}],
 QuantumInterferometer[{"S"[xi]}]["Probabilities", {0}, "Cutoff" -> 4, Assumptions -> xi > 0]}
```

The coincidence amplitude is proportional to $e^{i\arg x_2}\tanh|x_2|-e^{i\arg x_1}\tanh|x_1|$: it vanishes for equal squeezing, where the output is two independent squeezed vacua, and is largest for $x_2=-x_1$, where the splitter makes a two-mode squeezed vacuum. The squeezed vacuum populates only even photon numbers, $\operatorname{sech}r$, $\tfrac12\operatorname{sech}r\tanh^2r$, $\tfrac38\operatorname{sech}r\tanh^4r$.

## 14. Quantum Optics with QuantumInterferometer

Source: `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Quantum Interferometer/interferometer-showcase.md`.

Note for review: the closing "Scope" section and section 11 say lossy sampling is unsupported and that `"Samples"`/`"SampleCounts"` support only lossless passive networks. On the current kernel (after commits c91ca9ad and 24f5679f, "lossy and coherent samples") a lossy network samples without complaint, for example `QuantumInterferometer[{"BS"[Pi/4], "Loss"[0.5] -> 1}]["SampleCounts", {1, 1}, 1000]` returns counts; those two sentences look stale.

### Summary

Fourteen linear-optics experiments, each ending in a checked prediction: Hong-Ou-Mandel, Reck/Clements/cosine-sine meshes, Fourier multiports, squeezed light and Gaussian boson sampling, Mach-Zehnder fringes, coherent states and dark ports, twin-Fock interference, tunable four-photon suppression, heralded state preparation, multiphoton fringes, disentangling squeezers and Bloch-Messiah forms, pure loss, waveguide arrays and quantum walks, and the KLM nonlinear sign gate. Fock amplitudes are permanents, squeezed amplitudes hafnians, and loss acts on probabilities and finite states directly.

What should catch a physicist's attention:

- **Exact zeros are the signal.** HOM's coincidence is $\cos^2 2\theta$; a Fourier tritter forbids every $\{2,1,0\}$ pattern; with four modes 25 of 35 patterns vanish exactly, which floating point would blur.
- **Laws for every photon number.** Twin-Fock inputs $|n,n\rangle$ never give odd outputs, explained by $(z_2^2-z_1^2)^n$ and closed as $\binom{2k}{k}\binom{2n-2k}{n-k}/4^n$, normalized by `Sum` for all $n$; pure loss is binomial for symbolic $n$ with mean $n\eta$ proved by summation.
- **Tunable suppression.** $|2,2\rangle\to|2,2\rangle$ has amplitude $1-6R(1-R)$, so it vanishes at the two unbalanced reflectivities $R=(3\pm\sqrt3)/6$, and the zeros move with photon number.
- **Squeezed light, symbolically.** Opposite squeezers on a balanced splitter give a two-mode squeezed vacuum with geometric pair law $\tanh^{2n}r/\cosh^2r$; undoing the splitter returns two independent squeezed vacua exactly; a four-mode hafnian amplitude matches a six-level Fock-space circuit to machine precision.
- **Heralding and gates.** A photon plus weak coherent light, conditioned on one detected photon, approximates $(|0\rangle-|1\rangle)/\sqrt2$; the KLM nonlinear sign gate built from a three-mode real unitary gives amplitudes $\{1,1,-1\}/2$ on signal photon numbers 0, 1, 2, success probability exactly $1/4$ for every input, and a different amplitude at three photons, where the gate is undefined.
- **Cross-checks everywhere.** Permanents against Fock-space circuits, native loss against a vacuum-mode dilation, waveguide propagators $e^{-izH}$ against beam-splitter meshes, sample frequencies against exact distributions, a discrete-time quantum walk against the classical binomial walk.

### Example 14.1: the KLM nonlinear sign gate, exactly

A three-mode real unitary with one auxiliary photon; success is heralded by the auxiliary pattern $\{1,0\}$:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
nsMatrix = {{1 - Sqrt[2], 2^(-1/4), Sqrt[3/Sqrt[2] - 2]}, {2^(-1/4), 1/2, (1 - Sqrt[2])/2},
   {Sqrt[3/Sqrt[2] - 2], (1 - Sqrt[2])/2, Sqrt[2] - 1/2}};
nsGate = QuantumInterferometer[nsMatrix];
{RootReduce[nsMatrix . Transpose[nsMatrix]] == IdentityMatrix[3],
 RootReduce[Table[nsGate["Amplitude", {n, 1, 0} -> {n, 1, 0}], {n, 0, 3}]]}
```

The heralded amplitudes for 0, 1 and 2 signal photons are $\tfrac12,\tfrac12,-\tfrac12$: a sign flip on $|2\rangle$ alone, which no phase shifter can do, with success probability $1/4$ independent of the input. The three-photon amplitude, $(4\sqrt2-5)/2$, shows where the gate's domain ends.

### Example 14.2: squeezed pairs and tunable four-photon suppression

Opposite single-mode squeezers and a balanced splitter, with $r$ symbolic:

```wl
symbolicPair = QuantumInterferometer[{"S"[r], "S"[-r] -> 2, "BS"[Pi/4]}];
{FullSimplify[Table[symbolicPair["Probability", {0, 0} -> {n, n}], {n, 0, 2}], r > 0],
 FullSimplify[symbolicPair["Probability", {0, 0} -> {2, 0}], r > 0]}
```

Photons come only in pairs across the two modes, with the geometric law $\operatorname{sech}^2r\tanh^{2n}r$, and unequal patterns are exactly zero: a two-mode squeezed vacuum. Now the return amplitude of $|2,2\rangle$ on a beam splitter of angle $\theta$, and the reflectivities that kill it:

```wl
fourBS = QuantumInterferometer[{"BS"[\[Theta]]}];
fa = FullSimplify[fourBS["Amplitude", {2, 2} -> {2, 2}], \[Theta] \[Element] Reals];
{FullSimplify[fa == 1 - 6 Sin[\[Theta]]^2 Cos[\[Theta]]^2, \[Theta] \[Element] Reals],
 R /. Solve[1 - 6 R (1 - R) == 0 && 0 < R < 1, R, Reals]}
```

The amplitude is $1-6R(1-R)$, zero at $R=(3\pm\sqrt3)/6$: unlike HOM, four-photon cancellation needs an unbalanced splitter, and the two settings are mirror images under $R\leftrightarrow1-R$.

## 15. A Formal Semantics for QuantumFramework

Source: `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Formal Semantics/QF-formal-semantics.md` (written against QF 2.1.1, kernel `b13a9077`; its `.nb` predates the latest text edits).

### Summary

What does a QF object *mean*? The document answers object by object: what QF stores, and what it does with what it stores when it applies, composes, tensors, adjoints, discards, measures and evolves. The single thread is that a coefficient array means nothing until its subsystem dimensions, its output/input split and its bases are fixed, and QF stores exactly those. One running example, the operator $R$ whose stored coefficients are $\mathrm{diag}(1,-1)$ in the $X$ basis (so $R$ *is* $X$), carries every section.

What should catch a physicist's attention:

- **Numbers versus maps.** Six numbers can be a six-level qudit, a qubit and a qutrit, or a $2\times3$ map such as photon subtraction; the stored matrix $C$ and the represented map $A$ are related by $A=F_{\rm out}CF_{\rm in}^{-1}$ for any invertible (not necessarily orthonormal) basis.
- **Why bases must be reconciled before composing.** Composing $Z$ after $R$ by contracting the stored numbers would give the identity; QF changes coordinates first and returns $ZX$. The correctness requirement for composition is stated as a theorem about computational matrices.
- **One operator-space model for channels, evolution and phase space.** Vectorization turns maps on density matrices into `QuantumOperator`s on the operator space (`"Left"`, `"Right"`, `"Liouvillian"`); a channel is stored as a Stinespring isometry $V=\sum_m|m\rangle\otimes M_m$; the phase-space picture is a change of coordinates of the same matrices, so it preserves composition.
- **One isometry, three readings.** The extra index of $V$ can be kept as a live wire (coherent branches), traced (the channel) or read (the measurement); and the same list of matrices means Kraus operators to `QuantumChannel` but POVM effects to `QuantumMeasurementOperator` (stored as $\sqrt{E_m}$).
- **Honest edges.** Constructors check shapes, not physics (normalization, positivity, trace preservation are readable but not enforced); the adjoint is the physical one only in orthonormal bases (Gram factors are omitted); a partial trace is physical only when both sides of a wire carry the same basis; `QuantumDistance`'s default "Fidelity" is a distance, $1-\operatorname{Re}\operatorname{Tr}\sqrt{\rho_a\rho_b}$.
- **Placement in the literature.** QF is read against categorical quantum mechanics (CPM construction, dagger-compact structure, spiders), Proto-Quipper-style circuit languages and operator-algebraic semantics; the claimed combination no surveyed system has is a basis stored with every object, three pictures connected by transforms (phase space unrepresented elsewhere), Stinespring channels, and a phase-carrying stabilizer round trip.

### Example 15.1: the running operator, stored numbers against the represented map

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
R = QuantumOperator[{{1, 0}, {0, -1}}, "PauliX"];
{Normal@R["Matrix"], Normal@R["MatrixRepresentation"],
 Normal@QuantumOperator[QuantumOperator["X"], "PauliX"]["Matrix"]}
```

Declaring a basis turns the numbers $\mathrm{diag}(1,-1)$ into the operator $X$; converting the already-built $X$ into its own eigenbasis gives the same stored numbers back. Applying $R$ to $|0\rangle$ stores $(1,-1)/\sqrt2$ in the $X$ basis, which is $|1\rangle$ computationally:

```wl
{Normal@R[QuantumState["0"]]["StateVector"], Normal@R[QuantumState["0"]]["Computational"]["StateVector"]}
```

The case that makes the semantics necessary: compose $Z$ after $R$, and compare with contracting the stored numbers naively:

```wl
{Normal@QuantumOperator["Z"][R]["MatrixRepresentation"], PauliMatrix[3] . Normal[R["Matrix"]]}
```

QF returns $ZX=\begin{pmatrix}0&1\\-1&0\end{pmatrix}$; the naive contraction gives the identity, the wrong operator. The basis is part of the object, and composition reconciles it before summing.

### Example 15.2: one family of matrices, read as a channel and as a measurement

The Liouvillian of $X$ acting on $|0\rangle\langle0|$ is a map on density matrices stored as an operator; its output is the trace-zero commutator:

```wl
Normal@QuantumOperator["Liouvillian"["X"]][QuantumState["0"]]["DensityMatrix"]
```

A list read as Kraus operators is a channel (stored as a Stinespring isometry); a list read as effects is a POVM:

```wl
ch = QuantumChannel[{Sqrt[3/4] IdentityMatrix[2], Sqrt[1/4] PauliMatrix[1]}];
m = QuantumMeasurementOperator[{{{3/4, 0}, {0, 1/4}}, {{1/4, 0}, {0, 3/4}}}];
{Normal@ch[QuantumState["0"]]["DensityMatrix"], ch["TracePreservingQ"], m[QuantumState["0"]]["ProbabilitiesList"]}
```

The bit-flip channel sends $|0\rangle\langle0|$ to $\mathrm{diag}(3/4,1/4)$ and is trace preserving; the POVM with effects $\mathrm{diag}(3/4,1/4)$ and $\mathrm{diag}(1/4,3/4)$ gives Born probabilities $\{3/4,1/4\}$. Finally, measuring the $X$ observable on $R|0\rangle$ returns branch states that are *unnormalized*, each of norm $\sqrt{p_m}$:

```wl
mx = QuantumMeasurementOperator["X"][R[QuantumState["0"]]];
{mx["ProbabilitiesList"], Norm[Normal[#["StateVector"]]] & /@ mx["States"]}
```

Probabilities $\{1/2,1/2\}$ and two branch kets of norm $1/\sqrt2$: the branch, its probability and the conditional state are three distinct objects, and QF keeps them distinct.

## 16. Boson Sampling Explained

Source: `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Boson Sampling/boson-sampling.md`.

Note for review: the essay builds every network from `BeamSplitterOperator` and `PhaseShiftOperator` in truncated Fock space, and Part IV concludes that QuantumFramework loses to plain-WL photon-sector methods by orders of magnitude with few photons in many modes. `QuantumInterferometer`, which now computes the same probabilities as permanents with no Fock cutoff (Topics 13, 14), does not appear; Part IV's comparison may need a QI column before it is shown as a statement about QF.

### Summary

A computation-first tour from two photons on a beam splitter to the classical hardness of boson sampling, with every result cross-checked three ways: Fock-space simulation in QuantumFramework, the permanent formula, and plain Wolfram Language.

What should catch a physicist's attention:

- **HOM and its dip, symbolically.** The coincidence amplitude of $|1,1\rangle$ is $\cos2\theta$; with photons distinguished by a polarization angle $\chi$ the coincidence probability is $\cos^4\theta+\sin^4\theta-\tfrac12\sin^22\theta\cos^2\chi$, reducing at the balanced splitter to the HOM dip $\tfrac12\sin^2\chi$. Photon-number conservation is checked as an exact commutator in the untruncated ladder-operator algebra (`BosonicNormalOrder`).
- **Truncation is physics, not bookkeeping.** $n$ photons need $n+1$ levels per mode; one level too few either leaks probability or, with the matrix-exponential construction, freezes the interference entirely (or invents outputs a suppression law forbids).
- **Any interferometer from beam splitters.** A Reck decomposition written as a Givens-rotation fold builds any $U(m)$ from $m(m-1)/2$ splitters and $m$ phase shifters, symbolic entries included.
- **Every probability is a permanent; fermions get the determinant.** Checked symbolically for general complex matrices, including multiply-occupied inputs, and for a symbolic three-splitter mesh.
- **Suppression laws and their limits.** The Fourier tritter darkens exactly the patterns whose mode-label sum is not a multiple of 3 (Tichy et al.); at six photons the permanent finds dark patterns the law does not predict.
- **Bunching as a theorem.** Fully bunched outputs are exactly $n!$ times likelier for identical photons in every interferometer; the Haar-average bunching probability follows from Schur's lemma and matches a 1000-sample Monte Carlo.
- **Where the cost lives.** The truncated Fock space $(n+1)^m$ against the photon sector $\binom{n+m-1}{n}$; two photon-sector representations (a polynomial in creation operators and a `SymmetrizedArray` of permanents); Glynn's permanent doubling with every photon; a $1000\times1000$ determinant (fermions) cheaper than a 24-photon permanent.

### Example 16.1: Hong-Ou-Mandel with partially distinguishable photons

Path and polarization modes $1_H,2_H,1_V,2_V$; the second photon is polarized at angle $\chi$; a non-polarizing beam splitter acts on both polarizations:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
inPol = Cos[\[Chi]] FockState[{1, 1, 0, 0}, 3] + Sin[\[Chi]] FockState[{1, 0, 0, 1}, 3];
bsPol = QuantumCircuitOperator[{
    BeamSplitterOperator[{\[Theta], 0}, 3, {1, 2}, Method -> "Recurrence"],
    BeamSplitterOperator[{\[Theta], 0}, 3, {3, 4}, Method -> "Recurrence"]}];
coincidence = FullSimplify@ComplexExpand@Total@KeySelect[KeyMap[First, bsPol[inPol]["Probability"]],
     #[[1]] + #[[3]] == 1 && #[[2]] + #[[4]] == 1 &];
{FullSimplify[coincidence == Cos[\[Theta]]^4 + Sin[\[Theta]]^4 - Sin[2 \[Theta]]^2 Cos[\[Chi]]^2/2,
   (\[Theta] | \[Chi]) \[Element] Reals], FullSimplify[coincidence /. \[Theta] -> Pi/4, \[Chi] \[Element] Reals]}
```

The coincidence probability interpolates between identical photons ($\chi=0$, $\cos^22\theta$) and classical particles ($\chi=\pi/2$, $\cos^4\theta+\sin^4\theta$); at the balanced splitter it is $\tfrac12\sin^2\chi$, the HOM dip, with the polarization overlap $\cos^2\chi$ setting its depth.

### Example 16.2: permanents, suppression laws, and the n! of bunching

The permanent formula for identical and distinguishable photons:

```wl
modeList[occ_List] := Catenate@MapIndexed[ConstantArray[First[#2], #1] &, occ];
bosonProbability[w_, s_, t_] := Abs[Permanent[w[[modeList[t], modeList[s]]]]]^2/(Times @@ (s!) Times @@ (t!));
distinguishableProbability[w_, s_, t_] := Permanent[Abs[w[[modeList[t], modeList[s]]]]^2]/Times @@ (t!);
photonPatterns[n_, m_] := ReverseSort@FrobeniusSolve[ConstantArray[1, m], n];
```

The tritter's dark outputs are exactly those the suppression law predicts; at six photons the law misses some:

```wl
{AllTrue[photonPatterns[3, 3], PossibleZeroQ[bosonProbability[FourierMatrix[3], {1, 1, 1}, #]] ===
    ! Divisible[Range[0, 2] . #, 3] &],
 Length@Select[photonPatterns[6, 6], Divisible[Range[0, 5] . #, 6] &&
     PossibleZeroQ[bosonProbability[FourierMatrix[6], ConstantArray[1, 6], #]] &]}
```

`True`, and twelve six-photon patterns that the law allows are dark anyway: the law is sufficient, not necessary. For a fully general complex $6\times6$ transfer matrix, every all-in-one-mode output is exactly $n!$ times likelier for identical photons:

```wl
With[{w = Array[\[FormalW], {6, 6}]}, AllTrue[Range[2, 4], Function[n, AllTrue[n IdentityMatrix[6],
    FullSimplify[bosonProbability[w, PadRight[ConstantArray[1, n], 6], #] ==
      n! distinguishableProbability[w, PadRight[ConstantArray[1, n], 6], #]] &]]]]
```

The $n!$ is the number of identical-row permutations in the permanent, which the determinant would cancel instead: the exchange sign decides both bunching and computational cost.

## 17. Why Boson Sampling Is Hard

Source: `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Boson Sampling/why-boson-sampling-is-hard.md` (has uncommitted edits; written on `QuantumInterferometer`).

### Summary

A compact companion to Topic 16 that keeps three tasks apart, computing one amplitude, collecting every probability, and drawing one sample, because the hardness argument lives in the gap between them. Built entirely on `QuantumInterferometer`.

What should catch a physicist's attention:

- **Three particle statistics, one table.** On a balanced splitter, identical bosons never give a coincidence (permanent), distinguishable particles give $\tfrac14,\tfrac12,\tfrac14$ (permanent of $|W|^2$), and identical fermions always give one (determinant); the minus sign in $W$ is what cancels the bosonic coincidence.
- **The permanent derived, not postulated.** Each input creation operator maps to $\sum_jW_{jk}a_j^\dagger$; the output amplitudes are coefficients of the polynomial $\prod_k(\sum_jW_{jk}x_j)^{s_k}$, which is how QI's `"Probabilities"` computes the whole distribution, checked complex amplitude by complex amplitude against independent permanents.
- **Cost depends on occupation structure.** In random $2n$-mode networks, a collision-free or paired amplitude grows exponentially with $n$, while a fully bunched one stays essentially free because it reduces to $W_{n+1,1}^n$; collision-free amplitudes reach minutes near 28-29 photons.
- **Hard to compute is not hard to sample.** Distinguishable photons have permanent probabilities yet are sampled trivially column by column; what turns amplitude cost into sampling hardness is the Aaronson-Arkhipov argument (exact sampling would collapse the polynomial hierarchy; approximate sampling under Gaussian-permanent conjectures), alongside Clifford-Clifford's $O(n2^n)$ exact sampler.

### Example 17.1: bosons, classical particles and fermions on one beam splitter

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

Rows are the outputs $|2,0\rangle,|1,1\rangle,|0,2\rangle$; columns are identical bosons, distinguishable particles and identical fermions. Bosons bunch ($\tfrac12,0,\tfrac12$), classical particles split ($\tfrac14,\tfrac12,\tfrac14$), fermions antibunch ($0,1,0$): the only difference is the sign attached to exchanging two histories.

### Example 17.2: what one amplitude costs, and why structure matters

Collision-free, paired and fully bunched inputs in random networks with $m=2n$ modes, one QI amplitude each:

```wl
occupations[n_] := PadRight[#, 2 n] & /@ <|"collision-free" -> ConstantArray[1, n],
    "pairs" -> Append[ConstantArray[2, Quotient[n, 2]], Mod[n, 2]], "bunched" -> {n}|>;
randomInterferometer[n_] := BlockRandom[SeedRandom[1107 + 1000 n]; QuantumInterferometer["Random"[2 n]]];
timeAmplitude[q_, s_List] := With[{t = RotateRight[s, Total[s]]},
   ClearSystemCache[]; First@AbsoluteTiming[q["Amplitude", s -> t]]];
AssociationMap[Function[n, With[{q = randomInterferometer[n]}, timeAmplitude[q, #] & /@ occupations[n]]],
 {6, 12, 18, 20}]
```

From 12 to 20 photons the collision-free and paired amplitudes grow by more than two orders of magnitude, while the fully bunched one barely moves: a generic permanent pays the $2^n$ of Ryser and Glynn, a structured one does not. That is why worst-case #P-hardness is a statement about generic instances, and why the sampling-hardness argument needs random networks.
