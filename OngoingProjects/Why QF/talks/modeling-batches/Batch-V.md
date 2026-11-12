```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Maximum-likelihood decoding: exact logical error rates of three codes

A code stores one qubit in $n$ physical qubits. Each qubit independently suffers $X$, $Y$ or $Z$ with probabilities $p_X$, $p_Y$, $p_Z$, the checks are read without error, and a decoder turns the syndrome into a correction. How much lower is the probability $P_L$ of a logical error when the decoder picks the most probable class of errors instead of the lightest error?

Each qubit passes through the Pauli channel

$$
\mathcal E(\rho)=(1-p_X-p_Y-p_Z)\,\rho+p_X\,X\rho X+p_Y\,Y\rho Y+p_Z\,Z\rho Z,
$$

depolarizing when $p_X=p_Y=p_Z=p/3$ and $Z$-biased here as $p_X=p_Y=p/10$, $p_Z=p$. An error $E$ on the $n$ qubits has a syndrome, the checks it anticommutes with, and a class, the logical operators it anticommutes with. The lookup decoder answers each syndrome with the class of the lightest error that produces it; the maximum-likelihood decoder answers with the class of largest total probability. A round fails when the answer is not the class of the error that occurred:

$$
P_L=\sum_{E}\Pr(E)\,\big[\operatorname{class}(E)\neq\operatorname{guess}(\operatorname{syndrome}(E))\big].
$$

QECCore, a development package beside QF and not part of the shipped paclet, sums over all $4^n$ errors with $p$ symbolic, choosing each syndrome's most probable class at small $p$, and returns $P_L$ as a polynomial for each decoder; the panels plot $P_L/p^2$, whose value at $p=0$ is the coefficient of $p^2$:

```wl
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/QEC/QECCore/QECCore.wl"]; Needs["Wolfram`QuantumFramework`QEC`"];
mlDecCases = {{"5QubitCode", "five-qubit code, depolarizing", QECNoiseModel["Depolarizing", p]}, {"ShorCode", "Shor code, depolarizing", QECNoiseModel["Depolarizing", p]},
   {"SteaneCode", "Steane code, depolarizing", QECNoiseModel["Depolarizing", p]}, {"SteaneCode", "Steane code, Z-biased", QECNoiseModel[<|"X" -> p/10, "Y" -> p/10, "Z" -> p|>]}};
mlDecRates[code_, noise_] := Expand[QECLogicalErrorRate[QECCode[code], noise, "Decoder" -> #]] & /@ {"MaximumLikelihood", "MinimumWeight"};
mlDecPanel[{code_, title_, noise_}] := With[{r = mlDecRates[code, noise]},
   Plot[Evaluate[Expand[r/p^2]], {p, 0, 0.2}, PlotStyle -> {Automatic, Dashed}, PlotRange -> All, PlotLegends -> {"maximum likelihood", "lookup table"},
    Frame -> True, FrameLabel -> {"p", Row[{Subscript["P", "L"], " / ", Superscript["p", 2]}]}, ImageSize -> 360,
    PlotLabel -> Style[Column[{title, If[SameQ @@ r, Row[{Subscript["P", "L"], " = ", TraditionalForm[First[r]]}],
         Row[{"order ", Superscript["p", 2], ":  ", TraditionalForm[Coefficient[r[[1]], p, 2] p^2], "  against  ", TraditionalForm[Coefficient[r[[2]], p, 2] p^2]}]]}, Alignment -> Center], 12, FontFamily -> "Helvetica"]]];
mlDecPanels = mlDecPanel /@ mlDecCases;
Legended[Grid[Partition[First /@ mlDecPanels, 2], Spacings -> {2, 2}], Last[First[mlDecPanels]]]
```

The five-qubit code gives one polynomial for both decoders: its fifteen nonzero syndromes come from the fifteen one-qubit errors, so the lightest error already lies in the most probable class. The Shor code is degenerate, with many light errors that differ by a check and so share a class, and adding up their probabilities is where maximum likelihood lowers the $p^2$ coefficient most. The Steane code gains less, under depolarizing and $Z$-biased noise alike, and in every panel the order of the two curves holds over the whole range of $p$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Noise inside the circuit that measures the checks

A distance-three code corrects any one error on its qubits when its checks are read perfectly. Here a circuit reads them, and every gate, reset and readout in it fails with probability $p$. What does noise inside the measuring circuit do to the protection the code gives?

Each check $g$ is measured with a reset ancilla: every data qubit in the support of $g$ is rotated so that its letter becomes $Z$ ($H$ for $X$), CNOTs copy those qubits' parity into the ancilla, the rotations are undone, and the ancilla is read. After each one-qubit gate $X$, $Y$ or $Z$ strikes with probability $p/3$ each, after each CNOT one of the 15 two-qubit Paulis with probability $p/15$ each, and each reset and each readout is flipped with probability $p$; waiting qubits do not err. One noisy round is followed by a perfect readout of the data, and the decoder reads the round's syndrome and its change at the final readout. The transversal extraction instead measures each check through a verified cat state with one qubit per data qubit of $g$, each coupled to its data qubit by one two-qubit gate, and keeps only shots that pass the verification. As a reference, code capacity puts depolarizing noise of rate $p$ on the data and reads the checks perfectly.

QECCore, a development package beside QF and not part of the shipped paclet, returns $P_L$ for the five-qubit code (checks the cyclic shifts of $XZZXI$) as an exact polynomial in $p$, and the transversal rate at six values of $p$; the second panel plots $P_L/p$ for the bit-flip code (checks $Z_1Z_2$, $Z_2Z_3$) and the phase-flip code (checks $X_1X_2$, $X_2X_3$):

```wl
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/QEC/QECCore/QECCore.wl"]; Needs["Wolfram`QuantumFramework`QEC`"];
faultCircNoise = <|"code capacity" -> QECNoiseModel["Depolarizing", p], "circuit, all locations" -> QECNoiseModel["Circuit", p],
   "circuit, one-qubit gates noiseless" -> QECNoiseModel["Circuit", <|"TwoQubit" -> p, "Measurement" -> p, "Reset" -> p|>]|>;
faultCircRate[code_, key_] := Expand @ If[key === "code capacity", QECLogicalErrorRate[QECCode[code], faultCircNoise[key]], QECLogicalErrorRate[QECCode[code], faultCircNoise[key], "Rounds" -> 1]];
faultCircFive = faultCircRate["5QubitCode", #] & /@ {"code capacity", "circuit, all locations"};
faultCircTransversal = Table[{q, QECLogicalErrorRate[QECCode["5QubitCode"], QECNoiseModel["Circuit", q], "Rounds" -> 1, "Extraction" -> "Transversal"]}, {q, 10.^Range[-4, -1.5, 0.5]}];
faultCircDual = Delete[Flatten[Table[faultCircRate[code, key], {key, Keys[faultCircNoise]}, {code, {"BitFlipCode", "PhaseFlipCode"}}]], 5];
faultCircLabel[s_] := Style[s, 12, FontFamily -> "Helvetica"];
Grid[{{
   ListLogLogPlot[Append[Table[{q, # /. p -> q}, {q, 10^Range[-4, -1.5, 1/20]}] & /@ faultCircFive, faultCircTransversal], Joined -> {True, True, False}, PlotMarkers -> {None, None, Automatic},
    PlotLegends -> Placed[{"code capacity", "bare-ancilla circuit", "transversal circuit, kept shots"}, Below], Frame -> True, FrameLabel -> {"p", Subscript["P", "L"]}, ImageSize -> 380,
    PlotLabel -> faultCircLabel[Column[{"five-qubit code", Row[{"bare ancilla:  ", Subscript["P", "L"], " = ", TraditionalForm[Coefficient[Last[faultCircFive], p, 1] p], " + O(", Superscript["p", 2], ")"}]}, Alignment -> Center]]],
   Plot[Evaluate[Expand[Cancel[faultCircDual/p]]], {p, 0, 0.05}, PlotStyle -> {Automatic, Dashed, Automatic, Dashed, DotDashed}, PlotRange -> All, Frame -> True,
    FrameLabel -> {"p", Row[{Subscript["P", "L"], " / p"}]}, ImageSize -> 380,
    PlotLegends -> Placed[LineLegend[Automatic, Delete[Flatten[Table[Row[{code, ", ", key}], {key, Keys[faultCircNoise]}, {code, {"bit flip", "phase flip"}}]], 5], LegendLayout -> "Column"], Below],
    PlotLabel -> faultCircLabel[Column[{"bit-flip and phase-flip codes", Row[{"all locations:  ", TraditionalForm[Coefficient[faultCircDual[[3]], p, 1] p], "  against  ", TraditionalForm[Coefficient[faultCircDual[[4]], p, 1] p]}]}, Alignment -> Center]]]
   }}, Alignment -> Top, Spacings -> 2]
```

With the bare ancilla the five-qubit code fails at first order in $p$: one fault on an ancilla partway through its CNOTs spreads to several data qubits. The transversal extraction restores the second-order slope of code capacity on the kept shots, although at the largest $p$ it fails more often than the bare ancilla. The bit-flip and phase-flip codes fail identically at code capacity but not once their checks are measured by a circuit, because only the $X$ checks need $H$ gates on the data: with one-qubit-gate noise switched off, the phase-flip code falls back onto the bit-flip curve.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A diagonal magic circuit as one phase polynomial

Start $n$ qubits in $|+\rangle^{\otimes n}$, apply $T$ to every qubit, $CZ$ to every neighboring pair and one $CCZ$ to qubits 1, 2, 3. Each $T$ doubles a sum of stabilizer states, so kept that way the state would need $2^n$ terms. What does one amplitude of the result cost as $n$ grows?

Every gate is diagonal in the computational basis, so it only multiplies $|x\rangle$ by a power of $\omega=e^{i\pi/4}$: $T_q$ by $\omega^{x_q}$, $CZ$ by $\omega^{4x_qx_{q+1}}$, $CCZ$ by $\omega^{4x_1x_2x_3}$. The whole circuit is one phase polynomial over $\mathbb Z_8$,

$$
U|+\rangle^{\otimes n}=2^{-n/2}\sum_{x\in\{0,1\}^n}\omega^{\varphi(x)}|x\rangle,\qquad \varphi(x)=\sum_{q=1}^{n}x_q+4\sum_{q=1}^{n-1}x_qx_{q+1}+4x_1x_2x_3 \bmod 8,
$$

so an amplitude is the single term $\langle y|U|+\rangle^{\otimes n}=2^{-n/2}\,\omega^{\varphi(y)}$.

QF builds $\varphi$ with `"Compress" -> "PhasePolynomial"` and reads the amplitude of $|10\cdots0\rangle$; the cell times the build, one amplitude (fastest of five), and the same amplitude read from the dense state vector of the circuit with an $H$ layer in front, in seconds:

```wl
phasePolyGates[n_] := Join[Thread["T" -> Range[n]], Table["CZ" -> {q, q + 1}, {q, n - 1}], {"C"["Z" -> 3, {1, 2}]}];
phasePolyState[n_] := QuantumCircuitOperator[phasePolyGates[n]][Method -> {"Stabilizer", "Compress" -> "PhasePolynomial"}];
phasePolyY[n_] := PadRight[{1}, n];
phasePolyTimes = {
   Table[{n, First @ AbsoluteTiming[phasePolyState[n]]}, {n, 10, 100, 10}],
   Table[With[{s = phasePolyState[n]}, {n, Min @ Table[First @ AbsoluteTiming[s["Amplitude", phasePolyY[n]]], 5]}], {n, 10, 100, 10}],
   Table[{n, First @ AbsoluteTiming[QuantumCircuitOperator[Join[Thread["H" -> Range[n]], phasePolyGates[n]]][]["StateVector"][[2^(n - 1) + 1]]]}, {n, 4, 18, 2}]};
ListLogPlot[phasePolyTimes, Joined -> True, PlotMarkers -> Automatic, Frame -> True, FrameLabel -> {"qubits n", "seconds"},
  PlotLegends -> {"build the phase polynomial", "one amplitude from it", "dense state vector"},
  PlotLabel -> Style[Row[{"\[LeftAngleBracket]10\[CenterEllipsis]0|\[Psi]\[RightAngleBracket] at n = 20:  ", TraditionalForm[phasePolyState[20]["Amplitude", phasePolyY[20]]]}], 13, FontFamily -> "Helvetica"],
  ImageSize -> 480]
```

The dense state vector grows exponentially, a straight line on the log scale once its fixed overhead is passed. Building $\varphi$ and reading one amplitude from it grow only as a low power of $n$, bending away from a straight line, because an amplitude is $\varphi$ evaluated at one bit string rather than a sum over $2^n$ terms. The label is the single term at $n=20$: the $2^{-n/2}$ of a flat superposition times the phase $\omega$ that the one $T$ acting on $x_1=1$ wrote into it.
