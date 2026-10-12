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

A transmon is a superconducting circuit that behaves as a weakly anharmonic oscillator. A qubit is stored in its two lowest levels, and a microwave pulse at the qubit frequency rotates it between them. The anharmonicity is small, a few hundred megahertz against a qubit frequency of several gigahertz, so the transition from $|1\rangle$ to $|2\rangle$ lies close to the qubit transition from $|0\rangle$ to $|1\rangle$. A short pulse has a broad spectrum and drives both transitions, carrying population out of the qubit. It must also be strong to complete the rotation in time, and a strong drive shifts the qubit frequency, so the rotation misses its target unless the pulse's amplitude and frequency are calibrated. A long pulse avoids both problems but is exposed for longer to energy decay and dephasing.

The gate studied here flips the qubit, taking $|0\rangle$ to $|1\rangle$ and $|1\rangle$ to $|0\rangle$. It rotates the qubit by an angle π, so it is called a π pulse. The population a pulse carries out of the qubit, mostly into $|2\rangle$, is called leakage. A correction known as DRAG (derivative removal by adiabatic gate) reduces it by adding a second drive to the pulse, shaped so that, for a weak drive, the combined pulse has no component at the frequency of the transition from $|1\rangle$ to $|2\rangle$.

For a given pulse shape, how long the best π pulse is, and how much DRAG shortens it, has no closed-form answer. Both follow from the Schrödinger and Lindblad equations solved with the levels the qubit leaks into kept explicitly, and with the amplitude, frequency and DRAG coefficient of the pulse calibrated together with its length. For an anharmonicity of −200 MHz and coherence times $T_1$ and $T_2$ of 100 μs, DRAG lets the pulse shrink until energy decay and dephasing alone set its error. Without DRAG the pulse must stay long enough for leakage to die out, and its best error is larger.

## The Driven Transmon

The transmon is well approximated by an anharmonic oscillator. It is driven by a microwave pulse near its qubit frequency, with a DRAG correction added to the pulse:

$$H_{\text{lab}}(t) = \omega_q b^\dagger b + \frac{\alpha}{2} b^\dagger b^\dagger b b + \left[\Omega(t)\cos(\omega_d t) - \frac{\lambda\,\dot\Omega(t)}{\alpha}\sin(\omega_d t)\right]\left(b + b^\dagger\right)$$

The first two terms describe the transmon, and the last term is the drive. Here $b$ lowers the number of excitations by one, $\omega_q$ is the qubit frequency and $\alpha$ is the anharmonicity. For a transmon $\alpha$ is negative: the transition from $|1\rangle$ to $|2\rangle$ lies $|\alpha|$ below the qubit transition.

The drive oscillates at $\omega_d = \omega_q - \delta$, detuned from the qubit by $\delta$ so that it can follow the shift of the qubit frequency that a strong drive causes, and its envelope $\Omega(t)$ sets its strength. The second term in the bracket is the DRAG correction. It oscillates as $\sin(\omega_d t)$ instead of $\cos(\omega_d t)$, and its strength follows the slope of the envelope, scaled by the DRAG coefficient $\lambda$. Dividing the slope by $\alpha$ gives it the units of $\Omega$, so $\lambda$ is a pure number. At $\lambda = 0$ the correction is off and the pulse is a plain microwave drive. The amplitude of the envelope, $\lambda$ and $\delta$ are the three controls of the pulse.

The drive is simplest to describe in a frame rotating with it, which subtracts $\omega_d$ for every excitation. The drive terms then split into slow terms and terms that oscillate at twice the drive frequency. The fast terms average out because $\Omega \ll \omega_q$, and dropping them is the rotating-wave approximation. What remains is the Hamiltonian used throughout:

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

With the envelope defined, the Hamiltonian $H(t)$ becomes one operator, with five parameters in this order: the anharmonicity $\alpha$, the amplitude $A$, the DRAG coefficient $\lambda$, the detuning $\delta$ and the pulse length $\tau$. The `"Parameters"` list specifies the order in which we will supply their values:

```wl
model = QuantumOperator[
  (\[Alpha]/2) (b["Dagger"] @ b["Dagger"] @ b @ b) + \[Delta] (b["Dagger"] @ b) +
    (envelope/2) (b + b["Dagger"]) - (\[Lambda] D[envelope, t]/(2 \[Alpha])) (I (b["Dagger"] - b)),
  "Parameters" -> {\[Alpha], A, \[Lambda], \[Delta], \[Tau]}]
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

The drive is resonant with the $|0\rangle\to|1\rangle$ transition, while $|1\rangle\to|2\rangle$ is off resonance by $|\alpha|$. When $\Omega\ll|\alpha|$, excitation of $|2\rangle$ is weak. When $\Omega$ becomes comparable to $|\alpha|$, the drive can transfer substantial population into $|2\rangle$. For this pulse, the envelope peaks at $t=\tau/2=4\,\mathrm{ns}$. Compare its peak value with $|\alpha|$:

```wl
(envelope /. {A -> 1, \[Tau] -> 8, t -> 4})/Abs[anharmonicity]
```

<!-- => 0.625 -->

The peak drive is $0.625|\alpha|$, large enough to produce appreciable excitation of $|2\rangle$. We chose the short 8 ns pulse to make this excitation visible. For the same pulse shape and area, increasing the duration lowers the peak drive strength.

---

Plot the populations of $|0\rangle$, $|1\rangle$ and $|2\rangle$ during the pulse:

```wl
Plot[Evaluate[psi[t]["ProbabilitiesList"][[;; 3]]], {t, 0, 8},
   PlotLegends -> 
  LineLegend[
   Table[ColorData[97][j], {j, 4}], {"|0\[RightAngleBracket]", "|1\[RightAngleBracket]", 
    "|2\[RightAngleBracket]", "Pulse"}] ,
   PlotLabel -> "population change during a \[Pi] pulse",
   Frame -> True, 
 Prolog -> 
  With[{pnts = 
     Table[{t, envelope /. {A -> 1, \[Tau] -> 8}}, {t, 0, 
       8, .1}]}, {ColorData[97][4], Dashed, Line[pnts]}],
   FrameLabel -> {"t (ns)", "population"}]
```

---

Population reaches $|2\rangle$ through $|1\rangle$, since the drive couples only neighboring levels. Its population grows as $|1\rangle$ becomes occupied and the drive strengthens. As the drive weakens, some population returns to the qubit, but $|2\rangle$ does not empty completely.

Compute the populations at the end of the pulse:

```wl
psi[8]["Probabilities"]
```

The final population outside $|1\rangle$ has two contributions: population in $|0\rangle$, indicating an incomplete flip, and population in higher levels, mainly $|2\rangle$, which is leakage. We next calibrate the pulse controls to maximize the final population in $|1\rangle$, first without DRAG and then with it.

## DRAG and Calibration

We first ask how much calibration and DRAG can improve the transfer at the same pulse duration. Keeping the duration at 8 ns, we adjust the amplitude $A$ and detuning $\delta$ without DRAG, then adjust them together with the DRAG coefficient $\lambda$. We postpone optimizing the duration until we include energy decay and dephasing: a longer pulse can reduce leakage, but it also leaves more time for these processes to spoil the transfer. The next section includes both effects and optimizes the duration along with the controls.

---

Calibration requires simulating the pulse for many trial settings. We therefore construct the Schrödinger equation once, keeping the controls and duration as parameters, and solve it numerically for each set of values. Leaving the duration as a parameter also lets us evaluate pulses of different lengths later.

---

Set `"ReturnEquations" -> True` to obtain the Schrödinger equation and initial condition without solving them. The upper limit `1` in `{t, 0, 1}` is a placeholder; we will supply the actual interval `{t, 0, length}` to the numerical solver below:

```wl
{equations, state, range} =
  QuantumEvolve[
   model[anharmonicity, amplitude, drag, detuning, length],
   FockState[0, 5], {t, 0, 1}, "ReturnEquations" -> True];
```

Leave `amplitude`, `drag`, `detuning` and `length` unassigned so they remain parameters in the equations. For each trial pulse, we will supply their numerical values as arguments to the solution function constructed below.

---

Write the state as five amplitudes, each a function of time, corresponding in order to $|0\rangle$ through $|4\rangle$. Substituting this list into the equations expresses them in terms of the individual amplitudes. <code>[Normal]()</code> converts the sparse arrays in the equations to ordinary lists:

```wl
stateEntries = Array[\[Psi], 5];
stateEquations = Normal[equations] /. state -> Function[s, Evaluate[Through[stateEntries[s]]]];
```

---

Use `ParametricNDSolveValue` to construct `solution`. For each numerical choice of controls and duration, `solution[amplitude, drag, detuning, length]` returns the five amplitudes as functions of time. We raise the precision and accuracy goals because the transfer errors calculated later can be very small:

```wl
solution = ParametricNDSolveValue[stateEquations, stateEntries, {t, 0, length},
  {amplitude, drag, detuning, length}, PrecisionGoal -> 12, AccuracyGoal -> 14]
```

---

Define the transfer error as one minus the final population in $|1\rangle$, starting from $|0\rangle$. The second amplitude returned by `solution` corresponds to $|1\rangle$, and its squared magnitude gives that population. This error measures the transfer for this initial state only, not the full qubit gate.

```wl
transferError[amplitude_?NumericQ, drag_?NumericQ, detuning_?NumericQ, 
  length_?NumericQ] := 
 1 - Abs[solution[amplitude, drag, detuning, length][[2]][length]]^2
```

---

First calibrate the 8 ns pulse without DRAG ($\lambda=0$) by minimizing the transfer error over $A$ and $\delta$. We restrict the search to $A>0$, since reversing the drive’s sign leaves the populations unchanged when starting from $|0\rangle$. The detuning may have either sign:

```wl
calibratedNoDRAG = NMinimize[{transferError[a, 0, d, 8], a > 0}, {a, d}]
```

<!-- => {0.023499555881837497, {a -> 1.0746743862497976, d -> -0.202966140761959}} -->

---

Now calibrate the 8 ns pulse with DRAG, minimizing the transfer error over $A$, $\lambda$ and $\delta$ together. We keep the same amplitude and detuning constraints and search over non-negative DRAG coefficients ($\lambda\ge0$):

```wl
calibratedDRAG = NMinimize[{transferError[a, l, d, 8], a > 0, l >= 0}, {a, l, d}]
```

<!-- => {0.00506884263345686, {a -> 0.857641801469691, l -> 1.414593497582007, d -> 0.2394552525178285}} -->

---

Pass the calibrated parameter values to `solution` to obtain the five amplitudes as functions of time for each pulse:

```wl
psiNoDRAG = solution[Sequence @@ ({a, 0, d} /. Last[calibratedNoDRAG]), 8];
psiDRAG = solution[Sequence @@ ({a, l, d} /. Last[calibratedDRAG]), 8];
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

The calibrated pulses still leave population in higher levels. Increasing the duration allows a weaker drive and can reduce this leakage, but it also gives energy decay and dephasing more time to act. We now include these processes and optimize the duration together with the pulse controls to find the smallest final transfer error.

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
{openEquations, openState, openRange} = QuantumEvolve[model[anharmonicity, amplitude, drag, detuning, length],
   {b, b["Dagger"] @ b} -> {1/T1, 2/dephasingTime}, FockState[0, 5], {t, 0, 1}, "ReturnEquations" -> True];
```

---

Represent each entry of the $5\times5$ density matrix by a function of time. Substituting this matrix into the Lindblad equation gives the coupled equations for its 25 entries:

```wl
densityEntries = Array[\[Rho], {5, 5}];
densityEquations = Normal[openEquations] /. openState -> Function[s, Evaluate[Map[#[s] &, densityEntries, {2}]]];
```

---

Use `ParametricNDSolveValue` to construct `openSolution`. For each numerical choice of controls and duration, `openSolution[amplitude, drag, detuning, length]` returns a $5\times5$ matrix whose entries are functions of time. We use lower precision and accuracy goals here because the transfer errors with decay and dephasing are larger than the smallest errors without them:

```wl
openSolution = ParametricNDSolveValue[densityEquations, densityEntries, {t, 0, length},
  {amplitude, drag, detuning, length}, PrecisionGoal -> 10, AccuracyGoal -> 12]
```

---

Keep the same definition of transfer error: one minus the final population in $|1\rangle$. With the density matrix, this population is the second diagonal entry evaluated at the end of the pulse.

```wl
openError[a_?NumericQ, l_?NumericQ, d_?NumericQ, T_?NumericQ] := 1 - Re[openSolution[a, l, d, T][[2, 2]][T]]
```

---

A device is calibrated with its decoherence present, so the controls are now calibrated directly on this error, and the length is optimized together with them. With DRAG, start from the calibrated 8 ns pulse and let the length vary with the controls. The search is local: it follows that pulse to a length where the controls remove leakage almost completely, so decay and dephasing set what remains:

```wl
bestWithDRAG = With[{start = Last[calibratedDRAG]},
  FindMinimum[{openError[a, l, d, T], a > 0, l >= 0}, {{a, a /. start}, {l, l /. start}, {d, d /. start}, {T, 8}}]]
```

<!-- => {0.000054324240093550635, {a -> 0.8869111274290689, l -> 1.493374521005273, d -> 0.17506782738529808, T -> 10.22056757371642}} -->

---

The same pulse without decay and dephasing:

```wl
transferError[a, l, d, T] /. Last[bestWithDRAG]
```

<!-- => 3.590268227160465*^-8 -->

Without decay and dephasing the error of this pulse is far smaller, so they set nearly all of it.

---

Without DRAG, a local search from the calibrated 8 ns pulse stops at a nearby dip in the leakage, so the length is searched globally from 6 to 40 ns:

```wl
bestWithoutDRAG = NMinimize[{openError[a, 0, d, T], a > 0, 6 <= T <= 40}, {a, d, T}]
```

<!-- => {0.00012856907384595218, {a -> 1.0081657037356688, d -> -0.022326078002319033, T -> 22.3519815815493}} -->

---

The same pulse without decay and dephasing:

```wl
transferError[a, 0, d, T] /. Last[bestWithoutDRAG]
```

<!-- => 0.00001529904542985072 -->

Decay and dephasing also set most of the error of this pulse, which is long enough for leakage to have died down.

---

The density matrices during the two best pulses:

```wl
rhoNoDRAG = openSolution[Sequence @@ ({a, 0, d, T} /. Last[bestWithoutDRAG])];
rhoDRAG = openSolution[Sequence @@ ({a, l, d, T} /. Last[bestWithDRAG])];
```

---

The population of $|2\rangle$ during the two best pulses, on a logarithmic scale. The DRAG pulse is much shorter, so each pulse is drawn against the fraction of its own length:

```wl
LogPlot[{Re[rhoNoDRAG[[3, 3]][s (T /. Last[bestWithoutDRAG])]], Re[rhoDRAG[[3, 3]][s (T /. Last[bestWithDRAG])]]}, {s, 0, 1},
  PlotStyle -> {Directive[ColorData[97][3], Dashed], ColorData[97][3]},
  PlotLegends -> {"|2\[RightAngleBracket] no DRAG", "|2\[RightAngleBracket] DRAG"},
  PlotLabel -> "leakage during the best \[Pi] pulses", Frame -> True,
  FrameLabel -> {"fraction of the pulse", "population of |2\[RightAngleBracket]"}]
```

The DRAG pulse fills $|2\rangle$ more in the middle, because it is shorter and drives harder, and then empties it steeply at the end, finishing below the slower pulse without DRAG. It cannot empty $|2\rangle$ completely: dephasing scrambles the phase relation that brings the population back, so a small part stays behind, far more than the same pulse leaves without decay and dephasing.

---

The pulse without DRAG sends less population through $|2\rangle$ during the pulse only because it is slow. The length of the best DRAG pulse:

```wl
dragLength = T /. Last[bestWithDRAG]
```

<!-- => 10.22056757371642 -->

---

Calibrate the pulse without DRAG at that length. The search needs more than the default 100 iterations to converge:

```wl
sameLengthWithoutDRAG = NMinimize[{openError[a, 0, d, dragLength], a > 0}, {a, d}, MaxIterations -> 200]
```

<!-- => {0.002106884572926626, {a -> 1.0435618008821956, d -> -0.10934132396742124}} -->

---

Its error without decay and dephasing, which comes only from leakage and the rotation error:

```wl
transferError[a, 0, d, dragLength] /. Last[sameLengthWithoutDRAG]
```

<!-- => 0.0020523602727415824 -->

---

The population it leaves in $|2\rangle$:

```wl
Abs[solution[Sequence @@ ({a, 0, d} /. Last[sameLengthWithoutDRAG]), dragLength][[3]][dragLength]]^2
```

<!-- => 0.0019674447286612716 -->

At the same length, the pulse without DRAG has a much larger error than the DRAG pulse. Removing decay and dephasing barely changes this error. Almost all the population that misses $|1\rangle$ ends up in $|2\rangle$, so the error is leakage. DRAG removes this leakage, which is why the DRAG pulse can be this short.

---

The best pulses the two searches found, side by side, with the peak drive in units of $|\alpha|$:

```wl
summary[{e_, rules_}] := With[{r = Join[rules, {l -> 0}]},
  <|"length (ns)" -> T /. r, "error" -> e, "amplitude" -> a /. r, "DRAG coefficient" -> l /. r,
    "detuning (MHz)" -> 1000 d/(2 Pi) /. r, "peak drive / |\[Alpha]|" -> (2 Pi a/T)/Abs[anharmonicity] /. r|>];
Dataset[<|"with DRAG" -> summary[bestWithDRAG], "without DRAG" -> summary[bestWithoutDRAG]|>]
```

With DRAG the best pulse is short enough that its error comes almost entirely from decay and dephasing during the pulse: the calibrated controls remove leakage and the rotation error to far below that level. Without DRAG the pulse has to be much longer before leakage falls below the decoherence error, and the longer pulse pays for it in decay and dephasing. The price of the short DRAG pulse is a strong drive, a sizable fraction of the anharmonicity at its peak.

## Five Levels Are Enough

Every result above keeps five levels. Define the final populations of a transmon truncated to $n$ levels, with the same Hamiltonian:

```wl
finalPopulations[n_, {a_, l_, d_, T_}] := With[{c = AnnihilationOperator[n]},
  QuantumEvolve[QuantumOperator[
      (\[Alpha]/2) (c["Dagger"] @ c["Dagger"] @ c @ c) + \[Delta] (c["Dagger"] @ c) +
        (envelope/2) (c + c["Dagger"]) - (\[Lambda] D[envelope, t]/(2 \[Alpha])) (I (c["Dagger"] - c)),
      "Parameters" -> {\[Alpha], A, \[Lambda], \[Delta], \[Tau]}][anharmonicity, a, l, d, T],
    FockState[0, n], {t, 0, T}, PrecisionGoal -> 12, AccuracyGoal -> 14][T]["ProbabilitiesList"]]
```

---

The smallest error in this tutorial is that of the best DRAG pulse without decay and dephasing, where a change from the truncation would matter most. Its transfer error with five levels:

```wl
1 - finalPopulations[5, {a, l, d, T} /. Last[bestWithDRAG]][[2]]
```

<!-- => 3.590292496635783*^-8 -->

---

The same pulse with six levels:

```wl
1 - finalPopulations[6, {a, l, d, T} /. Last[bestWithDRAG]][[2]]
```

<!-- => 3.579106111484265*^-8 -->

A sixth level changes this error by well under a percent, so the levels above $|4\rangle$ do not affect the results.

## References

- J. Koch et al., "Charge-insensitive qubit design derived from the Cooper pair box," [Phys. Rev. A 76, 042319 (2007)](https://doi.org/10.1103/PhysRevA.76.042319).
- F. Motzoi, J. M. Gambetta, P. Rebentrost and F. K. Wilhelm, "Simple pulses for elimination of leakage in weakly nonlinear qubits," [Phys. Rev. Lett. 103, 110501 (2009)](https://doi.org/10.1103/PhysRevLett.103.110501).
- J. M. Gambetta, F. Motzoi, S. T. Merkel and F. K. Wilhelm, "Analytic control methods for high-fidelity unitary operations in a weakly nonlinear oscillator," [Phys. Rev. A 83, 012308 (2011)](https://doi.org/10.1103/PhysRevA.83.012308).
- P. Krantz et al., "A quantum engineer's guide to superconducting qubits," [Appl. Phys. Rev. 6, 021318 (2019)](https://doi.org/10.1063/1.5089550).
