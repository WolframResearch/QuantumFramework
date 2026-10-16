---
Template: TechNote
Name: TransmonGates
Title: Fast Gates on a Transmon
Context: Wolfram`QuantumFramework`
ContextPath: [Wolfram`QuantumFramework`SecondQuantization`]
CellContext: Global`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/TransmonGates
Keywords: [transmon, superconducting qubit, leakage, DRAG, pi pulse, anharmonicity, transfer error, gate error, pulse length, Lindblad, decoherence, T1, T2, optimization]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: [TimeEvolution, SecondQuantization]
---

A transmon is a superconducting circuit that behaves as a weakly anharmonic oscillator. A qubit is stored in its two lowest levels, and a microwave pulse at the qubit frequency rotates it between them. The anharmonicity is small, a few hundred megahertz against a qubit frequency of several gigahertz, so the transition from $|1\rangle$ to $|2\rangle$ lies close to the qubit transition from $|0\rangle$ to $|1\rangle$. A short pulse has a broad spectrum and drives both transitions, carrying population out of the qubit. It must also be strong to complete the rotation in time. A strong drive couples $|1\rangle$ to $|2\rangle$ and pushes the energy of $|1\rangle$, which moves the qubit transition away from the drive frequency, so the rotation misses its target unless the pulse's amplitude and frequency are calibrated. A long pulse avoids both problems but is exposed for longer to energy decay and dephasing.

The gate studied here flips the qubit, taking $|0\rangle$ to $|1\rangle$ and $|1\rangle$ to $|0\rangle$. It rotates the qubit by an angle π, so it is called a π pulse. The population a pulse carries out of the qubit, mostly into $|2\rangle$, is called leakage. A correction known as DRAG (derivative removal by adiabatic gate) reduces it by adding a second component to the microwave pulse.

For a given pulse shape, how long the best π pulse is, and how much DRAG shortens it, has no exact closed-form answer. Both follow from the Schrödinger and Lindblad equations solved with the levels the qubit leaks into kept explicitly, and with the amplitude, frequency and DRAG coefficient of the pulse calibrated together with its length. For an anharmonicity of −200 MHz and coherence times $T_1$ and $T_2$ of 100 μs, DRAG lets the pulse shrink until energy decay and dephasing dominate its error. Without DRAG the pulse must stay long enough for leakage to die out, and its best error is larger.

## The Driven Transmon

The transmon is well approximated by an anharmonic oscillator. It is driven by a microwave pulse near its qubit frequency, with a DRAG correction added to the pulse:

$$H_{\text{lab}}(t) = \omega_q b^\dagger b + \frac{\alpha}{2} b^\dagger b^\dagger b b + \left[\Omega(t)\cos(\omega_d t) - \frac{\lambda\,\dot\Omega(t)}{\alpha}\sin(\omega_d t)\right]\left(b + b^\dagger\right)$$

The first two terms describe the transmon, and the last term is the drive. Here $b$ lowers the number of excitations by one, $\omega_q$ is the qubit frequency and $\alpha$ is the anharmonicity. For a transmon $\alpha$ is negative: the transition from $|1\rangle$ to $|2\rangle$ lies $|\alpha|$ below the qubit transition.

The drive oscillates at $\omega_d = \omega_q - \delta$, detuned from the qubit by $\delta$, and its envelope $\Omega(t)$ sets its strength. A constant detuning compensates for the shift of the qubit transition that the drive causes. The second term in the bracket is the DRAG correction. It oscillates as $\sin(\omega_d t)$ instead of $\cos(\omega_d t)$, and its strength follows the slope of the envelope, scaled by the DRAG coefficient $\lambda$. Dividing the slope by $\alpha$ gives it the units of $\Omega$, so $\lambda$ is a pure number. For a weak drive and no detuning, $\lambda = 1$ removes the component of the drive at the leakage frequency. At $\lambda = 0$ the correction is off and the pulse is a plain microwave drive. The amplitude of the envelope, $\lambda$ and $\delta$ are the three controls of the pulse.

The drive is simplest to describe in a frame rotating with it, which subtracts $\omega_d$ for every excitation. In the rotating frame, some drive terms vary only with the pulse envelope. Others oscillate at $2\omega_d$, much faster than the envelope and populations change. These fast terms largely average out when $\Omega\ll\omega_q$; dropping them is the rotating-wave approximation. What remains is the Hamiltonian used throughout:

$$H(t) = \frac{\alpha}{2} b^\dagger b^\dagger b b + \delta\, b^\dagger b + \frac{\Omega(t)}{2}\left(b + b^\dagger\right) - \frac{\lambda\,\dot{\Omega}(t)}{2\alpha}\, i\left(b^\dagger - b\right)$$

In this frame the first two terms set the energies of the levels, and the last two couple neighboring levels. To write these terms as matrices, the oscillator is truncated to its five lowest levels, $|0\rangle$ to $|4\rangle$, which hold the qubit and the levels it leaks into; the last section checks how many levels are needed.

Define the lowering operator of a five-level oscillator:

```wl
b = AnnihilationOperator[5]
```

---

The drive couples neighboring levels with strengths 1, $\sqrt{2}$, $\sqrt{3}$ and 2, so it drives the transition from $|1\rangle$ to $|2\rangle$ more strongly than the qubit transition:

```wl
MatrixForm[(b + b["Dagger"])["Matrix"]]
```

---

The envelope is a raised cosine of length $\tau$ and amplitude $A$. It and its slope vanish at both ends, so the DRAG term also starts and ends at zero, and at $A=1$ its area is $\pi$:

```wl
envelope = A (Pi/\[Tau]) (1 - Cos[2 Pi t/\[Tau]])
```

---

Write $H(t)$ for any lowering operator `c`, so that the number of kept levels can be changed later:

```wl
hamiltonian[c_] := (\[Alpha]/2) (c["Dagger"] @ c["Dagger"] @ c @ c) + \[Delta] (c["Dagger"] @ c) +
   (envelope/2) (c + c["Dagger"]) - (\[Lambda] D[envelope, t]/(2 \[Alpha])) (I (c["Dagger"] - c));
```

---

Make $H(t)$ for five levels one operator. Its `"Parameters"` list fixes the order in which numbers are supplied: $\alpha$, $A$, $\lambda$, $\delta$ and $\tau$:

```wl
model = QuantumOperator[hamiltonian[b], "Parameters" -> {\[Alpha], A, \[Lambda], \[Delta], \[Tau]}]
```

---

Times are in nanoseconds and frequencies in radians per nanosecond, so an anharmonicity of $\alpha/2\pi=-200$ MHz is:

```wl
anharmonicity = -2 Pi 0.2
```

## Leakage

An uncalibrated 8 ns π pulse, with $A=1$, no DRAG and no detuning, starting from $|0\rangle$:

```wl
psi = QuantumEvolve[model[anharmonicity, 1, 0, 0, 8], FockState[0, 5], {t, 0, 8}]
```

---

The populations of the three lowest levels during the pulse, together with the shape of the pulse scaled to a peak of 1:

```wl
Plot[Evaluate[Append[psi[t]["ProbabilitiesList"][[;; 3]], envelope/(envelope /. t -> 4) /. {A -> 1, \[Tau] -> 8}]], {t, 0, 8},
  PlotStyle -> {Automatic, Automatic, Automatic, Directive[Gray, Dashed]},
  PlotLegends -> {"|0\[RightAngleBracket]", "|1\[RightAngleBracket]", "|2\[RightAngleBracket]", "pulse shape (scaled)"},
  PlotLabel -> "population change during a \[Pi] pulse",
  Frame -> True, FrameLabel -> {"t (ns)", "population"}]
```

Population reaches $|2\rangle$ through $|1\rangle$, since the drive couples only neighboring levels, and part of it stays there.

---

The populations at the end of the pulse:

```wl
psi[8]["Probabilities"]
```

---

Single-qubit gates on current devices last about 35 ns. The same pulse at that length:

```wl
QuantumEvolve[model[anharmonicity, 1, 0, 0, 35], FockState[0, 5], {t, 0, 35}][35]["Probabilities"]
```

Leakage is now small, but the flip is still incomplete, because the drive shifts the qubit transition.

## What DRAG Does

Define the final populations of a transmon truncated to $n$ levels, without decay and dephasing. With strict precision goals the default integration method can stop partway through some pulses, so an explicit Runge-Kutta method is used:

```wl
finalPopulations[n_, {amplitude_, drag_, detuning_, length_}] := QuantumEvolve[
   QuantumOperator[hamiltonian[AnnihilationOperator[n]], "Parameters" -> {\[Alpha], A, \[Lambda], \[Delta], \[Tau]}][
    anharmonicity, amplitude, drag, detuning, length],
   FockState[0, n], {t, 0, length}, Method -> "ExplicitRungeKutta", PrecisionGoal -> 12, AccuracyGoal -> 14][length]["ProbabilitiesList"]
```

---

At 8 ns, with $A=1$ and no detuning, vary $\lambda$ and compare the leakage, the final population outside the qubit levels, with the transfer error, the population missing from $|1\rangle$:

```wl
ListLogPlot[Transpose @ Table[With[{p = finalPopulations[5, {1, drag, 0, 8}]}, {{drag, Total[p[[3 ;;]]]}, {drag, 1 - p[[2]]}}], {drag, 0, 2, 0.1}],
  Joined -> True, PlotMarkers -> Automatic, PlotLegends -> {"leakage", "transfer error"},
  PlotLabel -> "8 ns \[Pi] pulse", Frame -> True, FrameLabel -> {"\[Lambda]", "population"}]
```

The transfer error is the population left in $|0\rangle$ plus the leakage, since the populations sum to one. Leakage is smallest just above the weak-drive value $\lambda=1$, but there almost all of the transfer error is population left in $|0\rangle$, and the transfer error is smallest at a much smaller $\lambda$: DRAG alone does not complete the flip. The amplitude and detuning are therefore calibrated together with $\lambda$.

## The Best Pulse Length

To find the best length, include energy decay and dephasing and, as on a real device, calibrate the pulse controls with them present.

We describe the state by a density matrix $\rho(t)$ in our truncated five-level model. Keeping the same Hamiltonian $H(t)$, it evolves under the Lindblad master equation:

$$\frac{d\rho}{dt} = -i[H(t),\rho] + \sum_{\mu\in\{1,\phi\}}\left(L_\mu\rho L_\mu^\dagger - \frac12 L_\mu^\dagger L_\mu\rho - \frac12\rho L_\mu^\dagger L_\mu\right).$$

The two Lindblad operators, with their rates included, are

$$L_1=\frac{b}{\sqrt{T_1}},\qquad L_\phi=\sqrt{\frac{2}{T_\phi}}\,b^\dagger b.$$

The first describes energy decay from $|n\rangle$ to $|n-1\rangle$ at rate $n/T_1$. The second describes pure dephasing: on its own, it damps the off-diagonal entries of $\rho$ while leaving the populations unchanged. The pure dephasing time $T_\phi$ is related to the qubit’s energy decay time $T_1$ and coherence time $T_2$ by

$$\frac{1}{T_\phi}=\frac{1}{T_2}-\frac{1}{2T_1}.$$

For this example, take $T_1=T_2=100\,\mu\mathrm{s}$. Express both times in nanoseconds:

```wl
T1 = T2 = 100000
```

---

Compute the corresponding pure dephasing time $T_\phi$, also in nanoseconds:

```wl
dephasingTime = 1/(1/T2 - 1/(2 T1))
```

---

Construct the Lindblad equation once, with the controls and the length as parameters: `"ReturnEquations" -> True` returns it without solving, and <code>[Most]()</code> drops the placeholder time range:

```wl
{openEquations, openState} = Most @ QuantumEvolve[model[anharmonicity, A, \[Lambda], \[Delta], \[Tau]],
   {b, b["Dagger"] @ b} -> {1/T1, 2/dephasingTime}, FockState[0, 5], {t, 0, 1}, "ReturnEquations" -> True];
```

---

Write it for the 25 entries of the density matrix, with <code>[Normal]()</code> turning sparse arrays into lists:

```wl
densityEntries = Array[\[Rho], {5, 5}];
densityEquations = Normal[openEquations] /. openState -> Function[s, Evaluate[Map[#[s] &, densityEntries, {2}]]];
```

---

Define `openSolution`, which solves these equations with `NDSolveValue` for numerical controls and length and returns the $5\times5$ density matrix as functions of time. The precision and accuracy goals are raised because the errors computed below are small:

```wl
openSolution[amplitude_?NumericQ, drag_?NumericQ, detuning_?NumericQ, length_?NumericQ] :=
 NDSolveValue[densityEquations /. {A -> amplitude, \[Lambda] -> drag, \[Delta] -> detuning, \[Tau] -> length},
  densityEntries, {t, 0, length}, PrecisionGoal -> 10, AccuracyGoal -> 12];
```

---

With the density matrix, the transfer error is one minus the second diagonal entry at the end of the pulse:

```wl
openError[amplitude_?NumericQ, drag_?NumericQ, detuning_?NumericQ, length_?NumericQ] :=
 1 - Re[openSolution[amplitude, drag, detuning, length][[2, 2]][length]]
```

---

Calibrate a pulse by minimizing its transfer error with <code>[NMinimize]()</code>: over $A$ and $\delta$ without DRAG, and also over $\lambda$ with DRAG. Both searches use the same bounds: $0<A<2$ keeps the pulse near area $\pi$, $0\le\lambda\le3$ is about twice the largest first-order leakage-cancelling value in the detuning range searched, and $|\delta|\le0.5$ rad/ns. The variables, with the ranges the search starts from, and the bounds:

```wl
controls[True] = {{amplitude, 0, 2}, {drag, 0, 3}, {detuning, -0.5, 0.5}};
controls[False] = {{amplitude, 0, 2}, {detuning, -0.5, 0.5}};
bounds = 0 < amplitude < 2 && 0 <= drag <= 3 && -0.5 <= detuning <= 0.5;
```

---

Define `calibrate`. The search uses differential evolution with 10 trial pulses, starting from a given random seed:

```wl
calibrate[withDRAG_, length_, seed_] := NMinimize[
  {openError[amplitude, drag, detuning, length], bounds} /. If[withDRAG, {}, drag -> 0], controls[withDRAG],
  Method -> {"DifferentialEvolution", "SearchPoints" -> 10, "RandomSeed" -> seed}]
```

---

Calibrate pulses from 6 to 40 ns, with and without DRAG. Each calibration runs with two random seeds and keeps the better result, so that a search stopping early does not show up as a false spike. The lengths are calibrated in parallel with <code>[ParallelMap]()</code>. This is the longest computation in the tutorial:

```wl
lengths = Join[Range[6, 24], Range[26, 40, 2]];
scan = AssociationMap[Function[withDRAG, ParallelMap[
     Function[length, First[MinimalBy[calibrate[withDRAG, length, #] & /@ {0, 1}, First]]], lengths,
     Method -> "FinestGrained"]], {True, False}];
```

---

The transfer error of the calibrated pulse at each length:

```wl
ListLogPlot[Transpose[{lengths, First /@ #}] & /@ Values[scan], Joined -> True, PlotMarkers -> Automatic,
  PlotLegends -> {"with DRAG", "without DRAG"}, PlotLabel -> "transfer error of the calibrated pulse",
  Frame -> True, FrameLabel -> {"pulse length (ns)", "transfer error"}]
```

Leakage dominates at short lengths, and decay and dephasing at long ones, where both curves rise together. DRAG removes most of the leakage, so its curve turns up at a much shorter length and reaches a lower minimum. Near 10 ns the DRAG curve also has a dip narrower than the 1 ns spacing of the scan, whose bottom the scan cannot resolve; Da Costa Jesus et al. found such a dip for this pulse shape.

---

The lowest point of each curve:

```wl
{bestWithDRAG, bestWithoutDRAG} = Function[withDRAG, With[{k = First[Ordering[First /@ scan[withDRAG], 1]]},
     {scan[withDRAG][[k, 1]], Append[scan[withDRAG][[k, 2]], length -> lengths[[k]]]}]] /@ {True, False};
```

---

The best pulses found in the scan: their lengths, errors with and without decay and dephasing, and controls:

```wl
summary[{error_, rules_}] := With[{pulse = {amplitude, drag, detuning, length} /. Append[rules, drag -> 0]},
   AssociationThread[{"length (ns)", "error", "error without decay\n and dephasing", "amplitude", "DRAG coefficient", "detuning (MHz)"},
    {pulse[[4]], ScientificForm@error, ScientificForm@(1 - finalPopulations[5, pulse][[2]]), pulse[[1]], pulse[[2]], 1000 pulse[[3]]/(2 Pi)}]];
Dataset[<|"with DRAG" -> summary[bestWithDRAG], "without DRAG" -> summary[bestWithoutDRAG]|>]
```

Without decay and dephasing both pulses are far better, so decay and dephasing dominate their errors. The DRAG pulse wins by being shorter.

## The Gate Error

The transfer error follows only $|0\rangle$. The gate error averages, over the six input states $|0\rangle$, $|1\rangle$, $(|0\rangle\pm|1\rangle)/\sqrt2$ and $(|0\rangle\pm i|1\rangle)/\sqrt2$, one minus the overlap of the final density matrix with the X image of the input. The Hamiltonian at time $\tau-t$ is the transpose of the Hamiltonian at time $t$, so without decay and dephasing the pulse flips $|1\rangle$ exactly as it flips $|0\rangle$. With them, the two errors stay close, as the gate errors below show.

Define the final density matrix of a pulse for any initial density matrix `rho0`. The first part of `densityEquations` is the Lindblad equation without its initial condition:

```wl
finalDensityMatrix[rho0_, {amplitude_, drag_, detuning_, length_}] := Map[#[length] &,
  NDSolveValue[{densityEquations[[1]], Map[#[0] &, densityEntries, {2}] == rho0} /.
     {A -> amplitude, \[Lambda] -> drag, \[Delta] -> detuning, \[Tau] -> length},
   densityEntries, {t, 0, length}, PrecisionGoal -> 10, AccuracyGoal -> 12], {2}];
```

---

Define the gate error:

```wl
gateError[amplitude_?NumericQ, drag_?NumericQ, detuning_?NumericQ, length_?NumericQ] := 1 - Mean[Map[
    With[{in = PadRight[#, 5], out = PadRight[Reverse[#], 5]},
      Re[Conjugate[out] . finalDensityMatrix[KroneckerProduct[in, Conjugate[in]], {amplitude, drag, detuning, length}] . out]] &,
    Normalize /@ {{1, 0}, {0, 1}, {1, 1}, {1, -1}, {1, I}, {1, -I}}]];
```

---

The gate error of the best DRAG pulse:

```wl
gateError[amplitude, drag, detuning, length] /. Last[bestWithDRAG]
```

---

The gate error of the best pulse without DRAG:

```wl
gateError[amplitude, 0, detuning, length] /. Last[bestWithoutDRAG]
```

## Five Levels Are Enough

Every result above keeps five levels. The leakage of the uncalibrated 8 ns pulse for a transmon truncated to three to six levels:

```wl
Table[{n, Total[finalPopulations[n, {1, 0, 0, 8}][[3 ;;]]]}, {n, 3, 6}]
```

---

The best DRAG pulse has the smallest error without decay and dephasing, so a change from the truncation would show most there. That error with five and with six levels:

```wl
AssociationMap[1 - finalPopulations[#, {amplitude, drag, detuning, length} /. Last[bestWithDRAG]][[2]] &, {5, 6}]
```

Three levels underestimate the leakage, while five and six agree for both pulses, so five levels are enough.

## References

- J. Koch et al., "Charge-insensitive qubit design derived from the Cooper pair box," [Phys. Rev. A 76, 042319 (2007)](https://doi.org/10.1103/PhysRevA.76.042319).
- F. Motzoi, J. M. Gambetta, P. Rebentrost and F. K. Wilhelm, "Simple pulses for elimination of leakage in weakly nonlinear qubits," [Phys. Rev. Lett. 103, 110501 (2009)](https://doi.org/10.1103/PhysRevLett.103.110501).
- J. M. Gambetta, F. Motzoi, S. T. Merkel and F. K. Wilhelm, "Analytic control methods for high-fidelity unitary operations in a weakly nonlinear oscillator," [Phys. Rev. A 83, 012308 (2011)](https://doi.org/10.1103/PhysRevA.83.012308).
- P. Krantz et al., "A quantum engineer's guide to superconducting qubits," [Appl. Phys. Rev. 6, 021318 (2019)](https://doi.org/10.1063/1.5089550).
- E. Rosenfeld, C. Gidney, G. Roberts et al., "Magic state cultivation on a superconducting quantum processor," [arXiv:2512.13908 (2025)](https://arxiv.org/abs/2512.13908).
- P. Stadler, M. Lodi, A. Khedri, R. Reiner, K. Bark, N. Vogt, M. Marthaler and J. Leppäkangas, "Demonstration of system-bath physics on a gate-based quantum computer," [arXiv:2404.18828 (2024)](https://arxiv.org/abs/2404.18828).
- J. D. Da Costa Jesus, B. Li, Y. Gao, R. Barends, F. A. Cárdenas-López and F. Motzoi, "Analytical blueprint for 99.999% fidelity X-gates on present superconducting hardware under strong driving," [arXiv:2512.19919 (2025)](https://arxiv.org/abs/2512.19919).
