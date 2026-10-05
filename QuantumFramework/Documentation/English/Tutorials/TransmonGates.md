---
Template: TechNote
Name: TransmonGates
Title: Fast Gates on a Transmon
Context: Wolfram`QuantumFramework`
ContextPath: [Wolfram`QuantumFramework`SecondQuantization`]
CellContext: Global`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/TransmonGates
Keywords: [transmon, superconducting qubit, leakage, DRAG, pi pulse, anharmonicity, gate error, pulse length, Lindblad, decoherence, T1, T2, qudit, optimization]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: [TimeEvolution, SecondQuantization]
Typeset: _SuperDagger -> StandardForm
---

A transmon is a superconducting circuit that behaves as a weakly anharmonic oscillator. A qubit is stored in its two lowest levels, and a microwave pulse at the qubit frequency rotates it between them. The anharmonicity is small, a few hundred megahertz against a qubit frequency of several gigahertz, so the transition from |1⟩ to |2⟩ lies close to the qubit transition. A short pulse has a broad spectrum and drives both, carrying population out of the qubit. A long pulse avoids that but is exposed for longer to energy decay and dephasing.

How long the best pulse is, and how much the derivative correction known as DRAG shortens it, has no closed-form answer. Both follow from the Schrödinger and Lindblad equations solved with the levels the qubit leaks into kept explicitly. For an anharmonicity of −200 MHz and coherence times $T_1$ and $T_2$ of 100 μs, the best π pulse below lasts 17.7 ns with an error of 9.8 × 10⁻⁵ when it uses DRAG, against 2.0 × 10⁻⁴ near 32 ns without it.

## Five Levels of a Transmon

In the frame rotating at the qubit frequency, the transmon is the anharmonic oscillator $H_0 = \frac{\alpha}{2} b^\dagger b^\dagger b b$, where $b$ lowers the excitation number by one and $\alpha$ is the anharmonicity. Five levels are kept here; the last section shows that fewer change the answer and more do not.

The lowering operator of a five-level oscillator:

```wl
b = AnnihilationOperator[5]
```

---

In the rotating frame the qubit levels |0⟩ and |1⟩ are degenerate, and the higher levels are shifted by $\alpha$, $3\alpha$ and $6\alpha$:

```wl
MatrixForm[Normal[((\[Alpha]/2) (SuperDagger[b] @ SuperDagger[b] @ b @ b))["Matrix"]]]
```

<!-- => MatrixForm[{{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, \[Alpha], 0, 0}, {0, 0, 0, 3*\[Alpha], 0}, {0, 0, 0, 0, 6*\[Alpha]}}] -->

---

A drive couples neighboring levels with strengths 1, $\sqrt{2}$, $\sqrt{3}$ and 2, so it drives the transition from |1⟩ to |2⟩ more strongly than the qubit transition:

```wl
MatrixForm[Normal[(b + SuperDagger[b])["Matrix"]]]
```

<!-- => MatrixForm[{{0, 1, 0, 0, 0}, {1, 0, Sqrt[2], 0, 0}, {0, Sqrt[2], 0, Sqrt[3], 0}, {0, 0, Sqrt[3], 0, 2}, {0, 0, 0, 2, 0}}] -->

## The Pulse

The drive has the cosine envelope $\Omega(t) = \frac{\pi}{\tau}\left(1 - \cos\frac{2\pi t}{\tau}\right)$, which starts and ends at zero. DRAG (Motzoi, Gambetta, Rebentrost and Wilhelm, 2009) adds an out-of-phase quadrature proportional to the derivative of the envelope, and a detuning $\delta$ moves the drive off the qubit frequency:

$$H(t) = \frac{\alpha}{2} b^\dagger b^\dagger b b + \delta b^\dagger b + \frac{\Omega(t)}{2}\left(b + b^\dagger\right) - \frac{\lambda \dot{\Omega}(t)}{2\alpha} i\left(b^\dagger - b\right)$$

The envelope, with $\tau$ the pulse length:

```wl
envelope = (Pi/\[Tau]) (1 - Cos[2 Pi t/\[Tau]])
```

<!-- => (Pi*(1 - Cos[(2*Pi*t)/\[Tau]]))/\[Tau] -->

---

Its area is $\pi$ for every length, which makes it a $\pi$ pulse:

```wl
Integrate[envelope, {t, 0, \[Tau]}]
```

<!-- => Pi -->

---

The whole Hamiltonian as one operator, with the anharmonicity, the DRAG coefficient, the detuning and the pulse length kept as parameters:

```wl
model = QuantumOperator[
  (\[Alpha]/2) (SuperDagger[b] @ SuperDagger[b] @ b @ b) + \[Delta] (SuperDagger[b] @ b) +
    (envelope/2) (b + SuperDagger[b]) - (\[Lambda] D[envelope, t]/(2 \[Alpha])) (I (SuperDagger[b] - b)),
  "Parameters" -> {\[Alpha], \[Lambda], \[Delta], \[Tau]}]
```

---

Numbers enter by supplying the parameters in that order. Times are in nanoseconds and frequencies in radians per nanosecond, so an anharmonicity of −200 MHz is:

```wl
anharmonicity = -2 Pi 0.2
```

<!-- => -1.2566370614359172 -->

## Leakage

An 8 ns pulse without DRAG or detuning, starting from |0⟩:

```wl
psi = QuantumEvolve[model[anharmonicity, 0, 0, 8], FockState[0, 5], {t, 0, 8}]
```

---

The populations of the three lowest levels during the pulse:

```wl
Plot[Evaluate[psi[t]["ProbabilitiesList"][[;; 3]]], {t, 0, 8},
  PlotStyle -> {StandardBlue, StandardOrange, StandardRed},
  PlotLegends -> {"|0\[RightAngleBracket]", "|1\[RightAngleBracket]", "|2\[RightAngleBracket]"},
  Frame -> True, FrameLabel -> {"t (ns)", "population"}]
```

---

Population passes through |2⟩ on its way to |1⟩, and part of it stays there. At the end of the pulse 9.6% of the population is in |2⟩:

```wl
psi[8]["ProbabilitiesList"]
```

<!-- => {0.06576708469399378, 0.8385151962984492, 0.0956994002347821, 0.000018318582408697576, 1.903662292943524*^-10} -->

---

Pulses are compared through a function of the three controls. Numeric code here names them `drag`, `detuning` and `length` rather than $\lambda$, $\delta$ and $\tau$: <code>[Table]()</code> and <code>[FindMinimum]()</code> localize their variables dynamically, so iterating over `\[Tau]` would also set it inside `model`. The final populations for a given DRAG coefficient, detuning and length:

```wl
populations[drag_?NumericQ, detuning_?NumericQ, length_?NumericQ] :=
  QuantumEvolve[model[anharmonicity, drag, detuning, length], FockState[0, 5], {t, 0, length}][length]["ProbabilitiesList"]
```

---

Leakage into |2⟩ for pulses from 4 to 24 ns:

```wl
leakageTable = Table[{length, populations[0, 0, length][[3]]}, {length, {4, 6, 8, 12, 16, 24}}]
```

<!-- => {{4, 0.5832544342748343}, {6, 0.3099802101994099}, {8, 0.0956994002347821}, {12, 0.0006402597186450065}, {16, 0.00017968211106241127}, {24, 0.000010554578361778353}} -->

---

On a logarithmic scale:

```wl
ListLogPlot[leakageTable, Joined -> True, Mesh -> All, PlotStyle -> StandardRed,
  Frame -> True, FrameLabel -> {"pulse length (ns)", "leakage"}]
```

Leakage falls by more than four orders of magnitude between 4 and 24 ns, but even at 24 ns a pulse without DRAG leaves 1.1 × 10⁻⁵ of the population in |2⟩.

## DRAG

To first order in $1/\alpha$, the DRAG quadrature removes the leakage at $\lambda$ = 1. The transfer error counts the population that does not reach |1⟩, leakage included. Leakage and transfer error at 8 ns for $\lambda$ = 0, 1/2 and 1:

```wl
Table[With[{p = populations[drag, 0, 8]}, {drag, p[[3]], 1 - p[[2]]}], {drag, {0, 1/2, 1}}]
```

<!-- => {{0, 0.0956994002347821, 0.16148480370155083}, {1/2, 0.02629621013523719, 0.029374548754029783}, {1, 0.001160165362067272, 0.08109951084166467}} -->

---

$\lambda$ = 1 cuts the leakage by a factor of 80 but leaves 8% of the population short of |1⟩. A scan over $\lambda$ shows the two errors reaching their minima at different values:

```wl
ListLogPlot[
  Transpose @ Table[With[{p = populations[drag, 0, 8]}, {{drag, p[[3]]}, {drag, 1 - p[[2]]}}], {drag, 0, 1.5, 0.125}],
  Joined -> True, Mesh -> All, PlotStyle -> {StandardRed, StandardBlue},
  PlotLegends -> {"leakage", "transfer error"}, Frame -> True, FrameLabel -> {"\[Lambda]", "error"}]
```

---

The leakage is smallest just above $\lambda$ = 1, and the transfer error near $\lambda$ = 0.6, so no single $\lambda$ serves both. A detuning $\delta$ added to the search resolves this. The transfer error as a function of the three controls:

```wl
transferError[drag_?NumericQ, detuning_?NumericQ, length_?NumericQ] := 1 - populations[drag, detuning, length][[2]]
```

---

The best DRAG coefficient and detuning for a 12 ns pulse:

```wl
optimum = FindMinimum[transferError[drag, detuning, 12], {{drag, 0.75}, {detuning, 0.03}}, Method -> "PrincipalAxis"]
```

<!-- => {0.00021363351241587836, {drag -> 0.6958278336084507, detuning -> 0.029704469880998167}} -->

---

The transfer error falls to 2.1 × 10⁻⁴ at $\lambda$ = 0.70 and a detuning of 0.030 rad/ns, 4.7 MHz. Without DRAG, the best detuning at 12 ns leaves ten times the error:

```wl
noDRAG = FindMinimum[transferError[0, detuning, 12], {detuning, 0}, Method -> "PrincipalAxis"]
```

<!-- => {0.002054588385573597, {detuning -> -0.07312769855059591}} -->

## Decoherence and the Best Pulse Length

Leakage falls with the pulse length, and decay and dephasing grow with it. Energy decay enters as the jump operator $b$ at rate $1/T_1$, which empties level $n$ at rate $n/T_1$. Pure dephasing enters as $b^\dagger b$ at rate $2/T_\phi$, with $1/T_\phi = 1/T_2 - 1/(2T_1)$. Both coherence times, 100 μs, in nanoseconds:

```wl
T1 = T2 = 100000
```

<!-- => 100000 -->

---

The pure dephasing time:

```wl
dephasingTime = 1/(1/T2 - 1/(2 T1))
```

<!-- => 200000 -->

---

The transfer error with decay and dephasing, from the Lindblad master equation:

```wl
openError[drag_?NumericQ, detuning_?NumericQ, length_?NumericQ] :=
  1 - QuantumEvolve[model[anharmonicity, drag, detuning, length],
      {b, SuperDagger[b] @ b} -> {1/T1, 2/dephasingTime}, FockState[0, 5], {t, 0, length}][length]["ProbabilitiesList"][[2]]
```

---

At the 12 ns optimum, decoherence adds 6 × 10⁻⁵:

```wl
openError[drag /. Last[optimum], detuning /. Last[optimum], 12]
```

<!-- => 0.0002768221344375821 -->

---

Across pulse lengths, $\lambda$ is held at its 12 ns optimum: it is dimensionless, and the derivative term scales with the pulse by itself. The detuning compensates a shift that grows with the drive power, which goes as $1/\tau^2$, so it is optimized again at each length, starting from the 12 ns value scaled by $(12/\tau)^2$. The error with decoherence at the best detuning:

```wl
bestError[drag_?NumericQ, reference_?NumericQ, length_?NumericQ] :=
  With[{opt = FindMinimum[transferError[drag, detuning, length], {detuning, reference (12/length)^2}, Method -> "PrincipalAxis"]},
    openError[drag, detuning /. Last[opt], length]]
```

---

With DRAG:

```wl
withDRAG = Table[{length, bestError[drag /. Last[optimum], detuning /. Last[optimum], length]}, {length, {10, 12, 16, 24, 40}}]
```

<!-- => {{10, 0.00395190079396146}, {12, 0.00027682220969571514}, {16, 0.00010701705683557883}, {24, 0.0001225838678494684}, {40, 0.00020070192706100887}} -->

---

Without DRAG:

```wl
withoutDRAG = Table[{length, bestError[0, detuning /. Last[noDRAG], length]}, {length, {12, 16, 24, 40, 60}}]
```

<!-- => {{12, 0.0021181989213668873}, {16, 0.0008053919564293688}, {24, 0.0002497622922331688}, {40, 0.00021702061423867214}, {60, 0.00030339291315584216}} -->

---

The error against pulse length, with and without DRAG:

```wl
ListLogLogPlot[{withDRAG, withoutDRAG}, Joined -> True, Mesh -> All,
  PlotStyle -> {StandardBlue, StandardOrange}, PlotLegends -> {"DRAG", "no DRAG"},
  Frame -> True, FrameLabel -> {"pulse length (ns)", "error"}]
```

---

Both curves have a minimum: to its left leakage and the residual shift dominate, to its right decay and dephasing. <code>[FindMinimum]()</code> locates it. With DRAG:

```wl
FindMinimum[bestError[drag /. Last[optimum], detuning /. Last[optimum], length], {length, 16, 20}]
```

<!-- => {0.00009843280913635066, {length -> 17.73811583486808}} -->

---

Without DRAG:

```wl
FindMinimum[bestError[0, detuning /. Last[noDRAG], length], {length, 30, 40}]
```

<!-- => {0.00020011547475240477, {length -> 31.713776892897897}} -->

DRAG roughly halves both the pulse length and its error, from 32 ns to 17.7 ns and from 2.0 × 10⁻⁴ to 9.8 × 10⁻⁵. The minimum without DRAG is flat. Started from 35 ns instead, the search stops at another length with almost the same error:

```wl
FindMinimum[bestError[0, detuning /. Last[noDRAG], length], {length, 35, 45}]
```

<!-- => {0.0002028144999299819, {length -> 34.88982838867432}} -->

Its location is soft and its value is not.

## Five Levels Are Enough

Every result above keeps five levels. The leakage of the 8 ns pulse without DRAG, for a transmon truncated to $d$ levels:

```wl
leakageWithLevels[d_] := With[{a = AnnihilationOperator[d]},
  QuantumEvolve[(anharmonicity/2) (SuperDagger[a] @ SuperDagger[a] @ a @ a) + ((envelope /. \[Tau] -> 8)/2) (a + SuperDagger[a]),
    FockState[0, d], {t, 0, 8}][8]["ProbabilitiesList"][[3]]]
```

---

From three to six levels:

```wl
Table[{d, leakageWithLevels[d]}, {d, 3, 6}]
```

<!-- => {{3, 0.0513114942431535}, {4, 0.09414255358445456}, {5, 0.0956994002347821}, {6, 0.09571468862246171}} -->

Three levels underestimate the leakage by almost half. Five and six agree to within 0.02%.