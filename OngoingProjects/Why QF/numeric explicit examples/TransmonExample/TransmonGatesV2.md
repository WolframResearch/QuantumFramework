---
Template: TechNote
Name: TransmonGates
Title: Fast Gates on a Transmon
Context: Wolfram`QuantumFramework`
ContextPath: [Wolfram`QuantumFramework`SecondQuantization`]
CellContext: Global`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/TransmonGates
Keywords: [transmon, superconducting qubit, leakage, DRAG, pi pulse, anharmonicity, transfer error, pulse length, Lindblad, decoherence, T1, T2, optimization]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: [TimeEvolution, SecondQuantization]
---

A transmon is a superconducting circuit that behaves as a weakly anharmonic oscillator. A qubit is stored in its two lowest levels, and a microwave pulse at the qubit frequency rotates it between them. The anharmonicity is small, a few hundred megahertz against a qubit frequency of several gigahertz, so the transition from $|1\rangle$ to $|2\rangle$ lies close to the qubit transition from $|0\rangle$ to $|1\rangle$. A short pulse has a broad spectrum and drives both transitions, carrying population out of the qubit. It must also be strong to complete the rotation in time. A strong drive couples $|1\rangle$ to $|2\rangle$ and pushes the energy of $|1\rangle$, which moves the qubit transition away from the drive frequency, so the rotation misses its target unless the pulse's amplitude and frequency are calibrated. A long pulse avoids both problems but is exposed for longer to energy decay and dephasing.

The gate studied here flips the qubit, taking $|0\rangle$ to $|1\rangle$ and $|1\rangle$ to $|0\rangle$. It rotates the qubit by an angle π, so it is called a π pulse. The population a pulse carries out of the qubit, mostly into $|2\rangle$, is called leakage. A correction known as DRAG (derivative removal by adiabatic gate) reduces it by adding a second drive to the pulse, shaped so that, for a weak drive, the combined pulse has no component at the frequency of the transition from $|1\rangle$ to $|2\rangle$.

For a given pulse shape, how long the best π pulse is, and how much DRAG shortens it, has no exact closed-form answer. Both follow from the Schrödinger and Lindblad equations solved with the levels the qubit leaks into kept explicitly, and with the amplitude, frequency and DRAG coefficient of the pulse calibrated together with its length. For an anharmonicity of −200 MHz and coherence times $T_1$ and $T_2$ of 100 μs, DRAG lets the pulse shrink until energy decay and dephasing alone set its error. Without DRAG the pulse must stay long enough for leakage to die out, and its best error is larger.

## The Driven Transmon

The transmon is well approximated by an anharmonic oscillator. It is driven by a microwave pulse near its qubit frequency, with a DRAG correction added to the pulse:

$$H_{\text{lab}}(t) = \omega_q b^\dagger b + \frac{\alpha}{2} b^\dagger b^\dagger b b + \left[\Omega(t)\cos(\omega_d t) - \frac{\lambda\,\dot\Omega(t)}{\alpha}\sin(\omega_d t)\right]\left(b + b^\dagger\right)$$

The first two terms describe the transmon, and the last term is the drive. Here $b$ lowers the number of excitations by one, $\omega_q$ is the qubit frequency and $\alpha$ is the anharmonicity. For a transmon $\alpha$ is negative: the transition from $|1\rangle$ to $|2\rangle$ lies $|\alpha|$ below the qubit transition.

The drive oscillates at $\omega_d = \omega_q - \delta$, detuned from the qubit by $\delta$, and its envelope $\Omega(t)$ sets its strength. The detuning lets the drive follow the qubit transition, which a strong drive moves by pushing $|1\rangle$ away from $|2\rangle$. The second term in the bracket is the DRAG correction. It oscillates as $\sin(\omega_d t)$ instead of $\cos(\omega_d t)$, and its strength follows the slope of the envelope, scaled by the DRAG coefficient $\lambda$. Dividing the slope by $\alpha$ gives it the units of $\Omega$, so $\lambda$ is a pure number. For a weak drive and no detuning, $\lambda = 1$ removes the component of the drive at the leakage frequency. At $\lambda = 0$ the correction is off and the pulse is a plain microwave drive. The amplitude of the envelope, $\lambda$ and $\delta$ are the three controls of the pulse.

The drive is simplest to describe in a frame rotating with it, which subtracts $\omega_d$ for every excitation. In the rotating frame, some drive terms vary only with the pulse envelope. Others oscillate at $2\omega_d$, much faster than the envelope and populations change. These fast terms largely average out when $\Omega\ll\omega_q$; dropping them is the rotating-wave approximation. What remains is the Hamiltonian used throughout:

$$H(t) = \frac{\alpha}{2} b^\dagger b^\dagger b b + \delta\, b^\dagger b + \frac{\Omega(t)}{2}\left(b + b^\dagger\right) - \frac{\lambda\,\dot{\Omega}(t)}{2\alpha}\, i\left(b^\dagger - b\right)$$

In this frame the first two terms set the energies of the levels, and the last two couple neighboring levels. To write these terms as matrices, the oscillator is truncated to its five lowest levels, $|0\rangle$ to $|4\rangle$, which hold the qubit and the levels it leaks into; the last section checks that a sixth level does not change the results.

Define the lowering operator of a five-level oscillator:

```wl
b = AnnihilationOperator[5]
```

---

Without detuning, the qubit levels $|0\rangle$ and $|1\rangle$ have the same energy in this frame, and the higher levels sit at $\alpha$, $3\alpha$ and $6\alpha$:

```wl
MatrixForm[((\[Alpha]/2) (b["Dagger"] @ b["Dagger"] @ b @ b))["Matrix"]]
```

<!-- => MatrixForm[{{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, \[Alpha], 0, 0}, {0, 0, 0, 3*\[Alpha], 0}, {0, 0, 0, 0, 6*\[Alpha]}}] -->

---

The drive couples neighboring levels with strengths 1, $\sqrt{2}$, $\sqrt{3}$ and 2:

```wl
MatrixForm[(b + b["Dagger"])["Matrix"]]
```

<!-- => MatrixForm[{{0, 1, 0, 0, 0}, {1, 0, Sqrt[2], 0, 0}, {0, Sqrt[2], 0, Sqrt[3], 0}, {0, 0, Sqrt[3], 0, 2}, {0, 0, 0, 2, 0}}] -->

The coupling between $|1\rangle$ and $|2\rangle$ is larger than the coupling between $|0\rangle$ and $|1\rangle$ by a factor of $\sqrt{2}$, so the drive acts more strongly on the leakage transition than on the qubit transition. The leakage transition is also off resonance by only $|\alpha|$. Both make leakage easy.

## The Pulse

The envelope of a pulse of length $\tau$ must vanish at both ends, so the drive switches on and off smoothly. Its slope must vanish there too, so the DRAG correction, which follows the slope, also starts and ends at zero. A constant and one cosine harmonic are the fewest terms that meet both conditions. They give the raised cosine used here:

$$\Omega(t) = A\,\frac{\pi}{\tau}\left(1 - \cos\frac{2\pi t}{\tau}\right)$$

It has no shape parameter beyond its amplitude $A$ and its length $\tau$.

Define the envelope:

```wl
envelope = A (Pi/\[Tau]) (1 - Cos[2 Pi t/\[Tau]])
```

<!-- => (A*Pi*(1 - Cos[(2*Pi*t)/\[Tau]]))/\[Tau] -->

---

Without DRAG and detuning, and restricted to $|0\rangle$ and $|1\rangle$, the drive rotates the qubit at the rate $\Omega(t)$, so the rotation angle is the area under the envelope:

```wl
Integrate[envelope, {t, 0, \[Tau]}]
```

<!-- => A*Pi -->

At $A=1$, the area is $\pi$, giving a complete flip for the resonant two-level model above. In our truncated five-level model, the drive also transfers amplitude into higher levels and back. Amplitude returning to $|1\rangle$ has accumulated phase according to those levels’ energies and interferes with the amplitude already there. This affects the transfer between $|0\rangle$ and $|1\rangle$, while some population may also remain in higher levels. An area of $\pi$ therefore no longer guarantees a complete flip. We calibrate $A$, together with the other pulse controls, rather than fix it at 1.

---

With the envelope defined, write the Hamiltonian $H(t)$ for any lowering operator `c`, so that $H(t)$ can be built for any number of kept levels:

```wl
hamiltonian[c_] := (\[Alpha]/2) (c["Dagger"] @ c["Dagger"] @ c @ c) + \[Delta] (c["Dagger"] @ c) +
   (envelope/2) (c + c["Dagger"]) - (\[Lambda] D[envelope, t]/(2 \[Alpha])) (I (c["Dagger"] - c));
```

---

Make $H(t)$ for the five-level oscillator one operator. Its `"Parameters"` list fixes the order in which numbers are supplied later: the anharmonicity $\alpha$, the amplitude $A$, the DRAG coefficient $\lambda$, the detuning $\delta$ and the pulse length $\tau$:

```wl
model = QuantumOperator[hamiltonian[b], "Parameters" -> {\[Alpha], A, \[Lambda], \[Delta], \[Tau]}]
```

---

We use $\hbar=1$ and measure time in nanoseconds, so $\alpha$, $\delta$ and $\Omega$ are expressed in radians per nanosecond. The amplitude factor $A$ and DRAG coefficient $\lambda$ are dimensionless. An anharmonicity of $\alpha/(2\pi)=-200\,\mathrm{MHz}$ corresponds to $\alpha=-2\pi\times0.2\,\mathrm{rad/ns}$:

```wl
anharmonicity = -2 Pi 0.2
```

<!-- => -1.2566370614359172 -->

## Leakage

We now examine how coupling to higher levels affects the transfer of population from $|0\rangle$ to $|1\rangle$. Before calibrating the controls, we follow the population during an 8 ns pulse of area $\pi$. Starting from $|0\rangle$, evolve the transmon with $A=1$, no DRAG ($\lambda=0$) and no detuning ($\delta=0$):

```wl
psi = QuantumEvolve[model[anharmonicity, 1, 0, 0, 8], FockState[0, 5], {t, 0, 8}]
```

---

The drive is resonant with the $|0\rangle\to|1\rangle$ transition, while $|1\rangle\to|2\rangle$ is off resonance by $|\alpha|$. When $\Omega\ll|\alpha|$, excitation of $|2\rangle$ is weak. When $\Omega$ becomes comparable to $|\alpha|$, the drive can transfer substantial population into $|2\rangle$. We chose a short 8 ns pulse, whose drive is strong, to make this excitation visible. For the same pulse shape and area, a longer pulse has a weaker drive.

---

Plot the populations of $|0\rangle$, $|1\rangle$ and $|2\rangle$ during the pulse, together with the shape of the pulse scaled to a peak of 1:

```wl
Plot[Evaluate[Append[psi[t]["ProbabilitiesList"][[;; 3]], envelope/(envelope /. t -> 4) /. {A -> 1, \[Tau] -> 8}]], {t, 0, 8},
  PlotStyle -> {Automatic, Automatic, Automatic, Directive[Gray, Dashed]},
  PlotLegends -> {"|0\[RightAngleBracket]", "|1\[RightAngleBracket]", "|2\[RightAngleBracket]", "pulse shape (scaled)"},
  PlotLabel -> "population change during a \[Pi] pulse",
  Frame -> True, FrameLabel -> {"t (ns)", "population"}]
```

---

Population reaches $|2\rangle$ through $|1\rangle$, since the drive couples only neighboring levels. Its population grows as $|1\rangle$ becomes occupied and the drive strengthens. As the drive weakens, some population returns to the qubit, but $|2\rangle$ does not empty completely.

Compute the populations at the end of the pulse:

```wl
psi[8]["Probabilities"]
```

The final population outside $|1\rangle$ has two contributions: population in $|0\rangle$, indicating an incomplete flip, and population in higher levels, mainly $|2\rangle$, which is leakage.

---

The short 8 ns pulse makes leakage easy to see. Single-qubit gates on current devices last about 25 to 35 ns. The final populations of the same uncalibrated π pulse lasting 30 ns:

```wl
QuantumEvolve[model[anharmonicity, 1, 0, 0, 30], FockState[0, 5], {t, 0, 30}][30]["Probabilities"]
```

Leakage into $|2\rangle$ is now small, and most of the error is population left in $|0\rangle$. While the drive is on, its coupling to $|2\rangle$ pushes the energy of $|1\rangle$ up, so the drive falls slightly out of resonance with the qubit and the flip stays incomplete. We next calibrate the pulse controls to maximize the final population in $|1\rangle$, first without DRAG and then with it.

## DRAG and Calibration

We first ask how much calibration and DRAG can improve the transfer at the same pulse duration. Keeping the duration at 8 ns, we adjust the amplitude $A$ and detuning $\delta$ without DRAG, then adjust them together with the DRAG coefficient $\lambda$. We postpone optimizing the duration until we include energy decay and dephasing: a longer pulse can reduce leakage, but it also leaves more time for these processes to spoil the transfer. The next section includes both effects and optimizes the duration along with the controls.

---

Calibration requires simulating the pulse for many trial settings. We therefore construct the Schrödinger equation once, keeping the controls and duration as parameters, and solve it numerically for each set of values. Leaving the duration as a parameter also lets us evaluate pulses of different lengths later.

---

Set `"ReturnEquations" -> True` to obtain the Schrödinger equation and initial condition without solving them. The upper limit `1` in `{t, 0, 1}` is a placeholder; the solver below integrates from 0 to the actual pulse length. <code>[Most]()</code> drops the time range, the last of the three results:

```wl
{equations, state} = Most @ QuantumEvolve[model[anharmonicity, A, \[Lambda], \[Delta], \[Tau]],
   FockState[0, 5], {t, 0, 1}, "ReturnEquations" -> True];
```

The equations keep the amplitude $A$, the DRAG coefficient $\lambda$, the detuning $\delta$ and the length $\tau$ as symbols. For each trial pulse, we will supply their numerical values as arguments to the solution function constructed below.

---

Write the state as five amplitudes, each a function of time, corresponding in order to $|0\rangle$ through $|4\rangle$. Substituting this list into the equations expresses them in terms of the individual amplitudes. <code>[Normal]()</code> converts the sparse arrays in the equations to ordinary lists:

```wl
stateEntries = Array[\[Psi], 5];
stateEquations = Normal[equations] /. state -> Function[s, Evaluate[Through[stateEntries[s]]]];
```

---

Define `solution`, which substitutes numerical values for the controls and the duration into these equations and solves them with `NDSolveValue`. For each choice of values, `solution[amplitude, drag, detuning, length]` returns the five amplitudes as functions of time. We raise the precision and accuracy goals because the transfer errors calculated later can be very small. With goals this strict, the default integration method can stop partway through some pulses, so we use an explicit Runge-Kutta method instead:

```wl
solution[amplitude_?NumericQ, drag_?NumericQ, detuning_?NumericQ, length_?NumericQ] :=
 NDSolveValue[stateEquations /. {A -> amplitude, \[Lambda] -> drag, \[Delta] -> detuning, \[Tau] -> length},
  stateEntries, {t, 0, length}, Method -> "ExplicitRungeKutta", PrecisionGoal -> 12, AccuracyGoal -> 14]
```

---

Define the transfer error as one minus the final population in $|1\rangle$, starting from $|0\rangle$. The second amplitude returned by `solution` corresponds to $|1\rangle$, and its squared magnitude gives that population. This error measures the transfer for this initial state only, not the full qubit gate.

```wl
transferError[amplitude_?NumericQ, drag_?NumericQ, detuning_?NumericQ, length_?NumericQ] :=
 1 - Abs[solution[amplitude, drag, detuning, length][[2]][length]]^2
```

---

First calibrate the 8 ns pulse without DRAG ($\lambda=0$) by minimizing the transfer error over $A$ and $\delta$. We restrict the search to $A>0$, since reversing the drive’s sign leaves the populations unchanged when starting from $|0\rangle$. We also keep $A<2$, so that the pulse stays near area $\pi$, the π pulse studied here. Every later search keeps $0 < A < 2$. The detuning may have either sign:

```wl
calibratedNoDRAG = NMinimize[{transferError[amplitude, 0, detuning, 8], 0 < amplitude < 2}, {amplitude, detuning}]
```

<!-- => {0.023499555881839496, {amplitude -> 1.074674453741852, detuning -> -0.2029661727584701}} -->

---

Now calibrate the 8 ns pulse with DRAG, minimizing the transfer error over $A$, $\lambda$ and $\delta$ together. We keep the same amplitude bounds and search over non-negative DRAG coefficients ($\lambda\ge0$), since the cancellation happens at positive $\lambda$:

```wl
calibratedDRAG = NMinimize[{transferError[amplitude, drag, detuning, 8], 0 < amplitude < 2, drag >= 0}, {amplitude, drag, detuning}]
```

<!-- => {0.005068842633449533, {amplitude -> 0.8576417252660532, drag -> 1.414594274763659, detuning -> 0.23945543584404172}} -->

The calibrated $\lambda$ is above 1: the detuning moves the leakage frequency to $\alpha + \delta$, and the drive is not weak.

---

`solution` evaluates only for numbers, so the calibrated values can be substituted after the call. The five amplitudes as functions of time for each pulse:

```wl
psiNoDRAG = solution[amplitude, 0, detuning, 8] /. Last[calibratedNoDRAG];
psiDRAG = solution[amplitude, drag, detuning, 8] /. Last[calibratedDRAG];
```

---

Compare the population of $|2\rangle$ during the two calibrated pulses. A logarithmic vertical scale makes the small populations remaining at the end easier to compare:

```wl
LogPlot[{Abs[psiNoDRAG[[3]][t]]^2, Abs[psiDRAG[[3]][t]]^2}, {t, 0, 8},
  PlotStyle -> {Directive[ColorData[97][3], Dashed], ColorData[97][3]},
  PlotLegends -> {"|2\[RightAngleBracket] no DRAG", "|2\[RightAngleBracket] DRAG"},
  PlotRange -> {10^-4, 1}, PlotLabel -> "leakage during the calibrated 8 ns \[Pi] pulses", Frame -> True,
  FrameLabel -> {"t (ns)", "population of |2\[RightAngleBracket]"}]
```

Both pulses appreciably populate $|2\rangle$ during the pulse. Toward the end of the pulse, this population decreases. The DRAG pulse (solid curve) has a somewhat lower peak and leaves less population in $|2\rangle$ than the pulse without DRAG (dashed curve). DRAG therefore reduces the final leakage, although neither pulse empties this level completely.

## Decoherence and the Best Pulse Length

The previous section compared calibrated pulses at one fixed duration. We now search for the best duration, with and without DRAG. Longer pulses can use weaker drives and reduce leakage, but they also allow more energy decay and dephasing. We therefore include these processes and, as on a real device, optimize the duration and pulse controls together with them present.

We describe the state by a density matrix $\rho(t)$ in our truncated five-level model. Keeping the same Hamiltonian $H(t)$, we include energy decay and pure dephasing through the Lindblad master equation:

$$\frac{d\rho}{dt} = -i[H(t),\rho] + \sum_{\mu\in\{1,\phi\}}\left(L_\mu\rho L_\mu^\dagger - \frac12 L_\mu^\dagger L_\mu\rho - \frac12\rho L_\mu^\dagger L_\mu\right).$$

The two Lindblad operators, with their rates included, are

$$L_1=\frac{b}{\sqrt{T_1}},\qquad L_\phi=\sqrt{\frac{2}{T_\phi}}\,b^\dagger b.$$

The first describes energy decay from $|n\rangle$ to $|n-1\rangle$ at rate $n/T_1$. The second describes pure dephasing: on its own, it damps the off-diagonal entries of $\rho$ while leaving the populations unchanged. The pure dephasing time $T_\phi$ is related to the qubit’s energy decay time $T_1$ and coherence time $T_2$ by

$$\frac{1}{T_\phi}=\frac{1}{T_2}-\frac{1}{2T_1}.$$

For this example, take $T_1=T_2=100\,\mu\mathrm{s}$. Express both times in nanoseconds:

```wl
T1 = T2 = 100000
```

<!-- => 100000 -->

---

Compute the corresponding pure dephasing time $T_\phi$, also in nanoseconds:

```wl
dephasingTime = 1/(1/T2 - 1/(2 T1))
```

<!-- => 200000 -->

---

Return the Lindblad equation and initial condition for a transmon starting in $|0\rangle$, keeping the pulse controls and duration symbolic. Here `QuantumEvolve` takes the operators $b$ and $b^\dagger b$ with their rates $1/T_1$ and $2/T_\phi$ supplied separately:

```wl
{openEquations, openState} = Most @ QuantumEvolve[model[anharmonicity, A, \[Lambda], \[Delta], \[Tau]],
   {b, b["Dagger"] @ b} -> {1/T1, 2/dephasingTime}, FockState[0, 5], {t, 0, 1}, "ReturnEquations" -> True];
```

---

Represent each entry of the $5\times5$ density matrix by a function of time. Substituting this matrix into the Lindblad equation gives the coupled equations for its 25 entries:

```wl
densityEntries = Array[\[Rho], {5, 5}];
densityEquations = Normal[openEquations] /. openState -> Function[s, Evaluate[Map[#[s] &, densityEntries, {2}]]];
```

---

Define `openSolution` in the same way as `solution`: it substitutes numerical values for the controls and the duration into these equations and solves them with `NDSolveValue`. For each choice of values, `openSolution[amplitude, drag, detuning, length]` returns a $5\times5$ matrix whose entries are functions of time. We use lower precision and accuracy goals here because the transfer errors with decay and dephasing are larger than the smallest errors without them:

```wl
openSolution[amplitude_?NumericQ, drag_?NumericQ, detuning_?NumericQ, length_?NumericQ] :=
 NDSolveValue[densityEquations /. {A -> amplitude, \[Lambda] -> drag, \[Delta] -> detuning, \[Tau] -> length},
  densityEntries, {t, 0, length}, PrecisionGoal -> 10, AccuracyGoal -> 12]
```

---

Keep the same definition of transfer error: one minus the final population in $|1\rangle$. With the density matrix, this population is the second diagonal entry evaluated at the end of the pulse.

```wl
openError[amplitude_?NumericQ, drag_?NumericQ, detuning_?NumericQ, length_?NumericQ] :=
 1 - Re[openSolution[amplitude, drag, detuning, length][[2, 2]][length]]
```

---

With DRAG, start from the calibrated 8 ns pulse and let the length vary with the controls. A global search over 6 to 40 ns finds a worse DRAG pulse, so this search is local: it follows that pulse to a length where the controls remove leakage almost completely, so decay and dephasing set what remains:

```wl
bestWithDRAG = FindMinimum[{openError[amplitude, drag, detuning, length], 0 < amplitude < 2, drag >= 0},
  Append[List @@@ Last[calibratedDRAG], {length, 8}]]
```

<!-- => {0.000054324240369885146, {amplitude -> 0.8869112753575197, drag -> 1.4933734845453894, detuning -> 0.1750676673532984, length -> 10.22056822192202}} -->

---

Without DRAG, a local search from the calibrated 8 ns pulse stops at a nearby dip in the leakage, so the length is searched globally from 6 to 40 ns:

```wl
bestWithoutDRAG = NMinimize[{openError[amplitude, 0, detuning, length], 0 < amplitude < 2, 6 <= length <= 40}, {amplitude, detuning, length}]
```

<!-- => {0.000128569073187923, {amplitude -> 1.0081652308238274, detuning -> -0.022324783128465024, length -> 22.35262249497928}} -->

---

Compare the two best pulses. For each, the table gives its length, its transfer error, the error of the same pulse without decay and dephasing, and its controls:

```wl
summary[{error_, rules_}] := Association[{"length (ns)" -> length, "error" -> error,
    "error without decay and dephasing" -> transferError[amplitude, drag, detuning, length],
    "amplitude" -> amplitude, "DRAG coefficient" -> drag, "detuning (MHz)" -> 1000 detuning/(2 Pi)} /. Append[rules, drag -> 0]];
Dataset[<|"with DRAG" -> summary[bestWithDRAG], "without DRAG" -> summary[bestWithoutDRAG]|>]
```

Without decay and dephasing, the error of each best pulse is far smaller, so decay and dephasing set most of it for both pulses. The DRAG pulse reaches this limit at less than half the length of the pulse without DRAG, so it is exposed to decay and dephasing for less time and ends with the smaller error. The price is a stronger drive.

---

The transfer error follows only $|0\rangle$, while the gate must also flip $|1\rangle$ and superpositions. The Hamiltonian at time $\tau - t$ is the transpose of the Hamiltonian at time $t$: the DRAG term keeps its form because both $\dot\Omega$ and $i(b^\dagger - b)$ change sign. Without decay and dephasing, the pulse therefore carries $|1\rangle$ to $|0\rangle$ with exactly the same amplitude as $|0\rangle$ to $|1\rangle$, so the flip needs no phase correction to act as an X gate. To check this with decay and dephasing included, define the gate error: start from each of the six states $|0\rangle$, $|1\rangle$, $(|0\rangle \pm |1\rangle)/\sqrt2$ and $(|0\rangle \pm i|1\rangle)/\sqrt2$, and average one minus the overlap of the final density matrix with the X image of the starting state. The first part of `densityEquations` is the Lindblad equation without its initial condition:

```wl
finalDensityMatrix[rho0_, {amplitude_, drag_, detuning_, length_}] := Map[#[length] &,
  NDSolveValue[{densityEquations[[1]], Map[#[0] &, densityEntries, {2}] == rho0} /.
     {A -> amplitude, \[Lambda] -> drag, \[Delta] -> detuning, \[Tau] -> length},
   densityEntries, {t, 0, length}, PrecisionGoal -> 10, AccuracyGoal -> 12], {2}];
gateError[amplitude_?NumericQ, drag_?NumericQ, detuning_?NumericQ, length_?NumericQ] := 1 - Mean[Map[
    With[{in = PadRight[#, 5], out = PadRight[Reverse[#], 5]},
      Re[Conjugate[out] . finalDensityMatrix[KroneckerProduct[in, Conjugate[in]], {amplitude, drag, detuning, length}] . out]] &,
    Normalize /@ {{1, 0}, {0, 1}, {1, 1}, {1, -1}, {1, I}, {1, -I}}]]
```

---

The gate errors of the two best pulses:

```wl
<|"with DRAG" -> gateError[amplitude, drag, detuning, length] /. Last[bestWithDRAG],
  "without DRAG" -> gateError[amplitude, 0, detuning, length] /. Last[bestWithoutDRAG]|>
```

<!-- => <|"with DRAG" -> 0.000055152679705350316, "without DRAG" -> 0.0001288515869743989|> -->

For both pulses the gate error is within a few percent of the transfer error, and the DRAG pulse still has the smaller error. Calibrating on the transfer from $|0\rangle$ therefore gives pulses that are equally good X gates.

---

The pulse without DRAG needs its greater length for leakage to die down. Calibrate it instead at the length of the best DRAG pulse:

```wl
dragLength = length /. Last[bestWithDRAG];
sameLengthWithoutDRAG = NMinimize[{openError[amplitude, 0, detuning, dragLength], 0 < amplitude < 2}, {amplitude, detuning}]
```

<!-- => {0.0021068825909453803, {amplitude -> 1.0435618584579909, detuning -> -0.1093413354401814}} -->

---

The population it leaves in $|2\rangle$:

```wl
Re[(openSolution[amplitude, 0, detuning, dragLength] /. Last[sameLengthWithoutDRAG])[[3, 3]][dragLength]]
```

<!-- => 0.00197070830613074 -->

At the same length, the pulse without DRAG has a much larger error than the DRAG pulse, and almost all the population that misses $|1\rangle$ ends up in $|2\rangle$: its error is leakage. DRAG removes this leakage, which is why the DRAG pulse can be this short.

## Five Levels Are Enough

Every result above keeps five levels. Define the final populations of a transmon truncated to $n$ levels, with the same Hamiltonian and the same solver settings as `solution`:

```wl
finalPopulations[n_, {amplitude_, drag_, detuning_, length_}] := QuantumEvolve[
   QuantumOperator[hamiltonian[AnnihilationOperator[n]], "Parameters" -> {\[Alpha], A, \[Lambda], \[Delta], \[Tau]}][
    anharmonicity, amplitude, drag, detuning, length],
   FockState[0, n], {t, 0, length}, Method -> "ExplicitRungeKutta", PrecisionGoal -> 12, AccuracyGoal -> 14][length]["ProbabilitiesList"]
```

---

The smallest error in this tutorial is that of the best DRAG pulse without decay and dephasing, where a change from the truncation would matter most. Its transfer error with five and with six levels:

```wl
AssociationMap[1 - finalPopulations[#, {amplitude, drag, detuning, length} /. Last[bestWithDRAG]][[2]] &, {5, 6}]
```

<!-- => <|5 -> 3.589819741467437*^-8, 6 -> 3.578671192716598*^-8|> -->

A sixth level changes this error by well under a percent, so the levels above $|4\rangle$ do not affect the results.

## References

- J. Koch et al., "Charge-insensitive qubit design derived from the Cooper pair box," [Phys. Rev. A 76, 042319 (2007)](https://doi.org/10.1103/PhysRevA.76.042319).
- F. Motzoi, J. M. Gambetta, P. Rebentrost and F. K. Wilhelm, "Simple pulses for elimination of leakage in weakly nonlinear qubits," [Phys. Rev. Lett. 103, 110501 (2009)](https://doi.org/10.1103/PhysRevLett.103.110501).
- J. M. Gambetta, F. Motzoi, S. T. Merkel and F. K. Wilhelm, "Analytic control methods for high-fidelity unitary operations in a weakly nonlinear oscillator," [Phys. Rev. A 83, 012308 (2011)](https://doi.org/10.1103/PhysRevA.83.012308).
- P. Krantz et al., "A quantum engineer's guide to superconducting qubits," [Appl. Phys. Rev. 6, 021318 (2019)](https://doi.org/10.1063/1.5089550).
