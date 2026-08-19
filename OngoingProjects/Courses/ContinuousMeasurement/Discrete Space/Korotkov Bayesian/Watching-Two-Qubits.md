---
Template: Default
---

# Watching Two Qubits at Once: Measurement-Induced Entanglement in the Quantum-Bayesian Picture

**Two qubits that never interact can be entangled by a measurement alone, if the measurement is built so it cannot tell $|01\rangle$ from $|10\rangle$. A single probe that reflects off both qubits reads out only the parity group, sorting the four computational states into three bins, $|00\rangle$, $|11\rangle$, and the odd pair, and each noisy record steers the pair stochastically into one of the three. This essay builds the two-qubit quantum-Bayesian update from the single-qubit one and runs it: we watch the odd-parity coherence ride along on the geometric mean of its two populations, exactly the single-qubit purity ride-along one level up, watch one record climb into a Bell state while its siblings collapse to a product pole, and pull the concurrence distribution out of an ensemble of records. It is bimodal, a spike at zero and a peak that presses against a sharp ceiling, and that ceiling is the closed-form maximum concurrence, the entanglement of the most fortunate record, which we recover from the simulated envelope. A final step traces the same update to the remote-qubit experiments it models and to the entanglement essay's place beside the non-commuting one.**

Mads Bahrami (last updated: August 18, 2026)

## Setting the Stage: How This Essay Flows

This essay is a computation-first tour of two qubits entangled by a joint continuous measurement and nothing else: how the two-qubit Bayesian update is the single-qubit one repeated on a larger state, why a measurement that hides the odd-parity subspace projects into a Bell state, how a single noisy record steers the pair, why the concurrence it can reach is capped in closed form, and what sets that cap. Every claim is computed on the spot from base Wolfram Language, so nothing here rests on a formula you cannot immediately test and change.

In other words, I have tried to build a small laboratory for entanglement by measurement. I strongly believe in a computation-first narrative for learning: in a sense, if I cannot compute it, I cannot claim to understand it. The rhythm throughout is concept, then computation, then interpretation.

The environment you see is a live Wolfram notebook. Evaluate the cells from top to bottom; the toolkit defined in the next section is used by every section that follows, so those dependencies matter. A couple of cells run small Monte-Carlo ensembles and take a few seconds. My suggestion is to focus on the output and its meaning first, then unpack the input code. You are not locked into any of it: change the rates, the step, the seed, and rerun your own experiments.

This is the two-qubit door into the continuous-measurement physics of the sibling essays. The core essay `Korotkov-Quantum-Bayesian.md` builds the single-qubit Bayes kick with the coherence riding along, which we here promote to two qubits; the essay `Watching-Two-Axes.md` watches one qubit along two incompatible axes, where this essay watches two qubits along one joint parity axis. The single-qubit purity ride-along is the prerequisite: the whole entangling mechanism is that same ride-along, one level up.

Let's start!

## Notation and Conventions: One Probe, Two Qubits, Three Bins

We work in units where the reduced Planck constant is one, $\hbar=1$. Two qubits, remote and never coupled, are jointly measured by a single probe that reflects off both, so its output records the qubit pair. Written in the computational basis $\{|00\rangle,|01\rangle,|10\rangle,|11\rangle\}$, the joint measurement is engineered so that $|01\rangle$ and $|10\rangle$ shift the probe by the same amount: it is a half-parity meter, distinguishing three bins, the state $|00\rangle$, the state $|11\rangle$, and the odd-parity pair $\{|01\rangle,|10\rangle\}$, but never resolving which of the two odd states it holds. The averaged record $V$ over a window is Gaussian about a center $\delta v_i$ that takes three values: $\delta v_1=-\delta v$ for $|00\rangle$, $\delta v_4=+\delta v$ for $|11\rangle$, and $\delta v_2=\delta v_3=0$ for the odd pair. One rate organizes the measurement, the strength $\delta v^2/s$ with $s=1/(2\eta_m)$ set by the measurement efficiency $\eta_m$; and one extra rate $\gamma$ carries whatever dephasing survives inside the protected odd subspace (imperfect matching of the two centers, photon loss between the cavities). We work with the so-called X-state, in which only the two-qubit density-matrix elements $\rho_{00,00}$, $\rho_{11,11}$, $\rho_{01,01}$, $\rho_{10,10}$, and the odd-parity coherence $\rho_{01,10}$ are appreciable, because the meter's distinguishability damps every other off-diagonal element fast.

Fix the five surviving X-state numbers as $\{x_1,x_2,x_3,x_4,x_5\}=\{\rho_{00,00},\rho_{01,01},\rho_{10,10},\rho_{11,11},|\rho_{01,10}|\}$. Encode the three-center Gaussian likelihoods for a reading $V$, and the concurrence that measures the entanglement of an X-state:

```wl
ClearAll[likelis, concurrence];
likelis[v_, dvv_, dvar_] := Exp[-({-dvv, 0., 0., dvv} - v)^2/(2 dvar)];
concurrence[x_] := 2 Max[0., x[[5]] - Sqrt[x[[1]] x[[4]]]]
```

The concurrence runs from zero for a separable or classically mixed state to one for a Bell state, and for this X-state it is set by a race: the odd-parity coherence $x_5=|\rho_{01,10}|$ pulling it up against the geometric mean $\sqrt{x_1 x_4}$ of the two even populations pulling it down.

Now the update the essay turns on. Recall from the single-qubit case that a Bayesian kick multiplies each population by its Gaussian likelihood and renormalizes, while the coherence rides along on the geometric mean of the two likelihoods it connects. The two-qubit version is identical, one level up: the four populations each pick up their likelihood, and the odd-parity coherence rides on $\sqrt{P_2 P_3}$, the geometric mean of the two likelihoods of the states it links, with an extra dephasing factor. Encode one honest Bayesian step: draw the reading from the mixture the state predicts, then update:

```wl
ClearAll[twoQubitStep];
twoQubitStep[dvv_, ss_, gt_, dt_][x_] := Module[{dvar = ss/(2 dt), v, p, nrm},
   v = RandomChoice[Clip[x[[1 ;; 4]], {0., 1.}] -> {-dvv, 0., 0., dvv}] +
     RandomVariate[NormalDistribution[0, Sqrt[dvar]]];
   p = likelis[v, dvv, dvar]; nrm = x[[1 ;; 4]] . p;
   {x[[1]] p[[1]], x[[2]] p[[2]], x[[3]] p[[3]], x[[4]] p[[4]],
     x[[5]] Sqrt[p[[2]] p[[3]]] Exp[-gt]}/nrm]
```

That is the entire integrator: draw which bin the record leans toward with the current populations as weights, bury it under window noise, and fold the reading back into the five numbers, the four populations by classical Bayes and the odd coherence by the geometric-mean ride-along.

## 1. The Odd-Parity Coherence Is the Purity Ride-Along, One Level Up

The whole entangling mechanism is the single-qubit purity ride-along promoted to the odd-parity subspace. In the single-qubit case, because the coherence multiplies by the geometric mean of the two population likelihoods, the ratio $|\rho_{10}|^2/(\rho_{11}\rho_{00})$ is left untouched by every ideal update, so a pure state stays pure. Here the same ratio $|\rho_{01,10}|^2/(\rho_{01,01}\rho_{10,10})$ is left untouched, so a coherent superposition inside the odd subspace stays coherent while the meter localizes the pair into that subspace. Confirm that an ideal update ($\gamma=0$) leaves the odd-parity ride-along ratio unchanged:

```wl
FullSimplify[
 With[{p2 = P2, p3 = P3}, (x5 Sqrt[p2 p3])^2/((x2 p2) (x3 p3))] == x5^2/(x2 x3),
 Assumptions -> {P2 > 0, P3 > 0, x2 > 0, x3 > 0}]
```

The ratio is algebraically conserved, so the odd-parity coherence rides its two populations exactly as a single qubit's coherence rides its own, which is why a measurement that only sorts the pair into the odd bin leaves a fully coherent Bell state rather than a classical mixture of $|01\rangle$ and $|10\rangle$. Entanglement here is not created by any interaction; it is the coherence that was already present surviving the localization, because the meter cannot look inside the subspace it is projecting onto.

## 2. The Half-Parity Meter: Three Bins, One Protected Subspace

Start the pair in a product of equal superpositions, so all four populations are $1/4$ and the odd coherence is $1/4$ as well. The meter now sorts each record into one of three bins. If the record leans to $\delta v=-\delta v$ or $+\delta v$, the pair collapses toward the product pole $|00\rangle$ or $|11\rangle$ and the odd populations drain away; if the record stays near zero, the pair is steered into the odd subspace, where the surviving coherence makes it the Bell state $(|01\rangle+|10\rangle)/\sqrt2$. Watch one record that lands in the odd bin: track its concurrence, the two even populations, and the total odd population:

```wl
With[{dvv = 1., ss = 2., gam = 0.15, dt = 0.02, nst = 300},
 With[{tt = Range[0, 300] 0.02,
   traj = BlockRandom[SeedRandom[7];
     NestList[twoQubitStep[dvv, ss, gam dt, dt], {0.25, 0.25, 0.25, 0.25, 0.25}, nst]]},
  ListLinePlot[{Transpose[{tt, concurrence /@ traj}],
    Transpose[{tt, traj[[All, 1]]}], Transpose[{tt, traj[[All, 4]]}],
    Transpose[{tt, traj[[All, 2]] + traj[[All, 3]]}]},
   PlotRange -> {0, 1.02}, Frame -> True, GridLines -> Automatic,
   PlotLegends -> {"concurrence", "\[Rho]00,00", "\[Rho]11,11", "odd population"},
   FrameLabel -> {"time (units of measurement time)", "value"},
   PlotLabel -> "one record steered into the entangled subspace", ImageSize -> 480]]]
```

As one can see, the two even populations drain toward zero while the odd population fills toward one, and the concurrence climbs and then eases back: the record has projected the pair into the odd subspace and the surviving coherence has become genuine entanglement. A different seed would drain the odd population instead and collapse to a product pole with zero concurrence; the outcome is stochastic, set by which bin the noisy record happens to favor, with the bin probabilities equal to the initial populations, the Born rule again.

## 3. The Sharp Ceiling: A Closed-Form Maximum Concurrence

How entangled can a record make the pair, and how does that depend on when we look? Because the whole X-state is fixed by the running readout $V$ together with the initial state and the two rates, the concurrence at any time is a definite function of $V$, and it is largest for the record that sits exactly at $V=0$, the most balanced odd-bin outcome. That most-fortunate concurrence is a closed form, a race between the extra dephasing $\gamma$ that erodes the coherence and the measurement rate $\delta v^2/s$ that first builds and then over-localizes it. Encode the maximum concurrence and confirm it is the concurrence-readout relation evaluated at $V=0$:

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

Both are true: the concurrence is stationary in the readout at $V=0$ and the value there is the ceiling, because $\cosh$ is smallest at zero and it sits in the denominator. No record, however lucky, can beat the $V=0$ outcome, so the concurrence a continuous measurement can produce has a sharp upper limit at every instant, set entirely by the competition of the two rates.

## 4. Pulling the Ceiling Out of an Ensemble of Records

Now the payoff. The ceiling was derived from the update; the records should press against it and never cross it. Run an ensemble from the product start and plot every trajectory's concurrence against the closed-form ceiling (this cell runs a small ensemble and takes a few seconds):

```wl
With[{dvv = 1., ss = 2., gam = 0.15, dt = 0.02, nst = 300, ntr = 60},
 With[{tt = Range[0, 300] 0.02},
  With[{cmat = Table[BlockRandom[SeedRandom[900 + k];
       concurrence /@ NestList[twoQubitStep[dvv, ss, gam dt, dt], {0.25, 0.25, 0.25, 0.25, 0.25}, nst]], {k, ntr}]},
   Show[
    ListLinePlot[Transpose[{tt, #}] & /@ cmat, PlotStyle -> Directive[Opacity[0.35], Gray],
     PlotRange -> {0, 0.45}, Frame -> True, GridLines -> Automatic,
     FrameLabel -> {"time (units of measurement time)", "concurrence"},
     PlotLabel -> "records press against the closed-form ceiling", ImageSize -> 500],
    Plot[cmaxCF[gam, dvv, ss, t], {t, 0, 300 0.02}, PlotStyle -> Directive[Black, Thick, Dashed]]]]]]
```

Every trajectory stays below the dashed ceiling, and the luckiest ones ride right along it: the closed-form maximum concurrence is exactly the sharp cutoff of the simulated distribution, recovered from records without ever imposing it. The ceiling rises from zero as the measurement builds coherence into concurrence, peaks, and then decays as continued localization and the extra dephasing win, so there is an optimal moment to stop measuring.

Look at where the trajectories actually end up, at a fixed late time, across a large ensemble:

```wl
With[{dvv = 1., ss = 2., gam = 0.15, dt = 0.02, nst = 300, ntr = 2000},
 With[{cend = Table[BlockRandom[SeedRandom[k];
      concurrence@Last@NestList[twoQubitStep[dvv, ss, gam dt, dt], {0.25, 0.25, 0.25, 0.25, 0.25}, nst]], {k, ntr}]},
  Histogram[cend, 40, Frame -> True, GridLines -> Automatic,
   FrameLabel -> {"final concurrence", "counts"},
   PlotLabel -> "bimodal: a spike at zero and a peak beneath the ceiling", ImageSize -> 460]]]
```

The distribution is bimodal: a tall spike at zero concurrence from the records that collapsed to $|00\rangle$ or $|11\rangle$, and a separate peak of entangled outcomes that pushes up against, but never past, the ceiling. The two peaks are the two fates of the pair, a product pole or the Bell subspace, and the gap between them is why the measurement must be conditioned, kept only when the record lands in the odd bin.

## 5. What Sets the Ceiling, and Where This Update Lives

The height of the ceiling is a trade-off between the extra dephasing $\gamma$ and the measurement rate. A slow, clean measurement builds coherence into concurrence before losing it; a noisy or lossy one erodes the coherence faster than it can localize the pair, and the peak concurrence drops. Read the peak of the ceiling for three values of the extra dephasing, keyed by $\gamma$:

```wl
With[{dvv = 1., ss = 2., dt = 0.01, nst = 800},
 AssociationThread[{0.05, 0.15, 0.4},
  Table[Max[cmaxCF[gam, dvv, ss, #] & /@ (Range[nst] dt)], {gam, {0.05, 0.15, 0.4}}]]]
```

As expected, the peak concurrence falls as the extra dephasing grows: a nearly ideal channel reaches deep into the entangled regime, while a lossy one barely leaves the separable one. This is exactly the trade-off that limited the first remote-qubit experiments, where two transmons a meter apart, joined only by a coaxial cable and a shared probe, reached a peak concurrence around a third before intrinsic dephasing and the finite efficiency of the amplifier stopped the climb.

The update we ran is the same one those experiments are analyzed with. Two qubits in separate cavities, probed in sequence so the signal reflects off one and then the other, are a cascaded system whose joint measurement is exactly this half-parity meter; the diagonal elements follow classical Bayes and the odd-parity coherence rides the geometric mean, which is the whole content of `twoQubitStep`. The most-likely-path analysis of the ensemble splits into three branches, the high-concurrence odd branch that hugs the ceiling and the two low-concurrence branches that collapse to the poles, matching the bimodal histogram. And the essay beside this one, `Watching-Two-Axes.md`, is the same physics seen the other way: there one qubit is watched along two incompatible axes and diffuses on a great circle; here two qubits are watched along one parity axis and are steered into a Bell state. In both, the state is moved not by any Hamiltonian but by the information a noisy record carries, folded back by Bayes' rule.

## Where This Leaves Us (and What Comes Next)

You now have a complete, computation-first toolkit for two qubits entangled by a joint half-parity measurement: a three-center likelihood, a Bayesian step that updates four populations classically and rides the odd-parity coherence on their geometric mean, a demonstration that this ride-along is the single-qubit purity ride-along one level up, a single record steered into a Bell state, a closed-form ceiling on the concurrence recovered from an ensemble of records, and the trade-off between measurement rate and dephasing that sets its height. Before moving on, the points that were computed and verified along the way:

- The two-qubit Bayesian update is the single-qubit one on a larger state: four populations update by classical Bayes, and the odd-parity coherence rides on the geometric mean of its two likelihoods.
- The odd-parity ride-along conserves $|\rho_{01,10}|^2/(\rho_{01,01}\rho_{10,10})$ exactly, so a measurement that only sorts the pair into the odd bin leaves a coherent Bell state, not a classical mixture; entanglement is the coherence surviving localization.
- A half-parity meter sorts each record into $|00\rangle$, $|11\rangle$, or the odd subspace, with bin probabilities equal to the initial populations; a record that lands in the odd bin is stochastically steered into a Bell state.
- The concurrence at any instant is a definite function of the readout, largest at $V=0$, so it has a sharp closed-form ceiling that rises from zero, peaks, and decays.
- An ensemble of records presses against that ceiling and never crosses it, and the final-concurrence distribution is bimodal: a spike at zero from the collapsed records and a peak beneath the ceiling from the entangled ones.
- The ceiling's height is a trade-off between the extra dephasing and the measurement rate, the same trade-off that capped the first remote-qubit experiments near a third.

The natural continuations are the pieces this essay set aside: the full four-by-four state with the fast-damping off-diagonal elements restored, the cavity and probe fields the half-parity meter is built from, real-time feedback that stabilizes the entanglement instead of only heralding it, and the many-qubit parity measurements that turn this single joint check into the continuously monitored error syndrome of a stabilizer code, which is where the non-commuting essay's Bacon-Shor connection and this essay's parity meter meet.
