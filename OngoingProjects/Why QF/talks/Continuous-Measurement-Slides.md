<!-- Build: MarkdownToNotebook["Continuous-Measurement-Slides.md", "Continuous-Measurement-Slides.nb", "EvaluateSeparator" -> None]. The loading cell is on the last slide; mark it as an initialization cell in the built notebook. -->

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

# Continuous Quantum Measurement in the Wolfram Language

#| style: Subtitle

A watched quantum system, one problem per slide

#| style: Subsubtitle

Mads Bahrami, Wolfram Research, Nov 2026

#| style: Subsubtitle

quantum@wolfram.com

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### How the Wolfram Language handles a watched quantum system

Every slide poses one physical question, writes the model as an equation, and lets the language either solve it exactly or sample it, with the closed form and the simulation checked against each other. A few ideas recur:

- **Write the stochastic master equation once.** The conditioned Bloch-vector equations become one `ItoProcess`, and every step reuses it. `RandomFunction` samples quantum trajectories. The `"Derivative"` property applies Itô's lemma: for the purity $r^2$ it shows whether watching purifies the state and whether a trajectory can leave the Bloch ball, and for the log-likelihood ratio $\operatorname{arctanh}z$ it shows when the integrated readout alone fixes the population. `StratonovichProcess` gives the Stratonovich form, which `NDSolve` integrates along a smooth noise path.
- **The average over records is a Lindblad equation in QuantumFramework.** `QuantumOperator["Liouvillian"[H, jumps, rates]]` gives the superoperator as a matrix, and `QuantumEvolve` with jump operators and rates evolves the averaged state. This average is the reference line the single runs are drawn against.
- **Exact statements come from the Liouvillian matrix.** `Discriminant` of its `CharacteristicPolynomial` locates exceptional points. `Reduce` with `Exists` turns "is there a measurement strength where..." into a condition on the parameters. `NullSpace` gives the steady state, and `LinearSolve` of $i\omega-\mathcal L$ gives the detector's noise spectrum in closed form.
- **A filter is a fold over the record.** Each reading is drawn with `RandomVariate` from Born's probabilities blurred by detector noise, and `FoldList` or `NestList` applies Bayes' rule to whatever the filter keeps: two amplitudes, a density matrix, five Gaussian moments or four syndrome probabilities. The same readings drive a full-state simulation, so the filter can be checked against it.
- **Operator algebra replaces matrices where it can.** `NonCommutativeAlgebra` with `NonCommutativeExpand` derives the moment equations of a watched Gaussian from $[x,p]=i\hbar$, and the syndrome equation of a watched code from the Pauli relations, with `SolveAlways` matching coefficients for a general state. The Riccati equations for the widths are then solved symbolically, with `Resolve` and `ForAll` picking the physical root.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A driven qubit under continuous measurement

A qubit is driven with Rabi frequency $\Omega$, detuned by $\Delta$ from its transition, while a detector continuously measures $\hat\sigma_z$ at strength $k$ and catches a fraction $\eta$ of the signal. Does watching stop the oscillation?

With $\hbar=1$ and $\hat\sigma_z=|0\rangle\langle0|-|1\rangle\langle1|$, the driven qubit in the frame rotating at the drive frequency is

$$
\hat H=\frac{\omega_q}{2}\hat\sigma_z+\Omega\cos(\omega_d t)\,\hat\sigma_x\;\longrightarrow\;\frac{\Delta}{2}\hat\sigma_z+\frac{\Omega}{2}\hat\sigma_x,\qquad \Delta=\omega_q-\omega_d,
$$

and the state conditioned on the detector record obeys the Itô stochastic master equation

$$
d\rho=-i[\hat H,\rho]\,dt+k\,\mathcal D[\hat\sigma_z]\rho\,dt+\sqrt{\eta k}\left(\hat\sigma_z\rho+\rho\hat\sigma_z-2\langle\hat\sigma_z\rangle\rho\right)dW .
$$

Its average over records, without the $dW$ term, is the Lindblad equation. For the Bloch vector $(x,y,z)=(\langle\hat\sigma_x\rangle,\langle\hat\sigma_y\rangle,\langle\hat\sigma_z\rangle)$ it reads

$$
\begin{aligned}
dx&=(-2kx-\Delta y)\,dt-2\sqrt{\eta k}\,xz\,dW,\\
dy&=(\Delta x-2ky-\Omega z)\,dt-2\sqrt{\eta k}\,yz\,dW,\\
dz&=\Omega y\,dt+2\sqrt{\eta k}\,(1-z^2)\,dW.
\end{aligned}
$$

Three runs start from $|0\rangle$ ($z=1$), sampled with a step that shrinks as $k$ grows, against the Lindblad average (red, dashed); each run plots $\langle\hat\sigma_z\rangle$ conditioned on its own record, and rates are in units of $\Omega$:

```wl
measuredQubitSME[k_, η_, Δ_] := ItoProcess[{{-2 k x[t] - Δ y[t], Δ x[t] - 2 k y[t] - z[t], y[t]},
    {{-2 Sqrt[η k] x[t] z[t]}, {-2 Sqrt[η k] y[t] z[t]}, {2 Sqrt[η k] (1 - z[t]^2)}}}, {{x, y, z}, {0, 0, 1}}, {t, 0}];
measuredQubitAverage[k_, Δ_] := QuantumEvolve[Δ/2 QuantumOperator["Z"] + QuantumOperator["X"]/2,
    {QuantumOperator["Z"]} -> {k}, QuantumState["0"], {t, 0, 30}]["BlochVector"][[3]];
measuredQubitPanel[k_, {η_, Δ_}] := Show[
    ListLinePlot[#[[;; ;; Round[0.01/Min[0.01, 0.005/k]]]] & /@ RandomFunction[measuredQubitSME[k, η, Δ], {0., 30., Min[0.01, 0.005/k]}, 3, Method -> "StochasticRungeKuttaScalarNoise"]["PathComponent", 3]["Paths"], PlotRange -> {-1, 1}],
    Plot[Evaluate[measuredQubitAverage[k, Δ]], {t, 0, 30}, PlotRange -> {-1, 1}, PlotStyle -> Directive[Red, Dashed]],
    PlotRange -> {-1, 1}, Frame -> True, FrameLabel -> {"\[CapitalOmega] t", Row[{"\[LeftAngleBracket]", Subscript["\[Sigma]", "z"], "\[RightAngleBracket]"}]},
    Epilog -> Inset[Row[{"\[Eta] = ", η}], Scaled[{0.96, 0.96}], {Right, Top}], ImageSize -> 300];
measuredQubitPanels = BlockRandom[SeedRandom[1]; AssociationMap[Function[kd, measuredQubitPanel[First[kd], {#, Last[kd]}] & /@ {0.3, 1}], Tuples[{{0.1, 0.5, 1, 5, 10}, {0, 0.5, 1, 2}}]]];
measuredQubitLesson = <|0.1 -> "Runs oscillate; the average decays only because they fall out of phase.",
   0.5 -> "Runs still oscillate, kicked irregularly by the record.",
   1 -> "Crossover: oscillation stops and runs linger near the poles.",
   5 -> "Zeno: runs stick at a pole and jump; the average only decays.",
   10 -> "Stronger watching, rarer jumps."|>;
Manipulate[Legended[PlotGrid[{measuredQubitPanels[{k, d}]}, ImageSize -> 800],
   Placed[LineLegend[Append[ColorData[97] /@ Range[3], Directive[Red, Dashed]], {"run 1", "run 2", "run 3", "Lindblad average"}, LegendLayout -> "Row"], Below]],
  Row[{Column[{Control[{{k, 0.1, "measurement strength k/\[CapitalOmega]"}, {0.1, 0.5, 1, 5, 10}, SetterBar}], Control[{{d, 0, "detuning \[CapitalDelta]/\[CapitalOmega]"}, {0, 0.5, 1, 2}, SetterBar}]}],
    Spacer[20], Dynamic[Grid[{
      {Style[Row[{"k/\[CapitalOmega] = ", k}], Bold, 12], measuredQubitLesson[k]},
      {Style["\[Eta] < 1", Bold, 12], "Runs stay mixed and miss the poles."},
      {Style["\[CapitalDelta] \[NotEqual] 0", Bold, 12], "Weak watching only: smaller swing, slow drift between poles."}},
     Alignment -> {Left, Top}, ItemSize -> {{6, 44}}, Spacings -> {1, 0.6}, BaseStyle -> {12, LineIndent -> 0}]]}],
  SaveDefinitions -> True]
```

At $k/\Omega=0.1$ the runs keep Rabi-oscillating. At $k/\Omega=10$ they stay pinned near $z=1$, with brief spikes that fall back and now and then a full jump to the other pole: that is the quantum Zeno effect, which the average shows only as a slow decay. A detector that catches less of the signal leaves each run mixed, so under weak measurement its oscillation no longer reaches the poles. Detuning shrinks the Rabi swing to a fast wiggle near the starting pole, and under weak watching the record then carries the runs slowly from pole to pole; once the measurement is much stronger than the detuning, the runs look as they do on resonance.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The quantum Zeno crossover as an exceptional point

A qubit is driven with Rabi frequency $\Omega$, detuned by $\Delta$ from the drive, while a detector measures $\hat\sigma_z$ at strength $k$. Weak watching damps the Rabi oscillation and strong watching freezes it; what happens to the dynamics at the crossover?

With $\hbar=1$, $|e\rangle=|0\rangle$, $|g\rangle=|1\rangle$, $\hat\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, qubit frequency $\omega$ and drive frequency $\omega-\Delta$, the Hamiltonian in the frame rotating with the drive is $\hat H=\frac{\Delta}{2}\hat\sigma_z+\frac{\Omega}{2}\hat\sigma_x$, and averaging over the detector record gives

$$
\dot\rho=\mathcal L\rho=-i[\hat H,\rho]+k\,\mathcal D[\hat\sigma_z]\rho .
$$

Away from exceptional points, every solution is a sum of $e^{\lambda t}$ over the eigenvalues $\lambda$ of the superoperator $\mathcal L$. At an exceptional point two eigenvalues and their eigenvectors coincide, so $\mathcal L$ cannot be diagonalized; there the discriminant of its characteristic polynomial, with the steady state's root $\lambda=0$ divided out, vanishes.

On resonance the Rabi pair solves $\lambda^2+2k\lambda+\Omega^2=0$, so $\lambda=-k\pm\sqrt{k^2-\Omega^2}$: a damped oscillation for $k<\Omega$, and for $k>\Omega$ two real rates, one of which, about $\Omega^2/2k$, slows as the watching grows stronger, the Zeno effect. The two meet at $k=\Omega$. The discriminant of $\mathcal L$'s characteristic polynomial, with $\lambda=0$ divided out, has a zero at some $k>0$, an exceptional point, only for:

```wl
liouA5 = Normal[QuantumOperator["Liouvillian"[Δ/2 QuantumOperator["Z"] + Ω/2 QuantumOperator["X"], {QuantumOperator["Z"]}, {k}]]["Matrix"]];
discA5 = Discriminant[Cancel[CharacteristicPolynomial[liouA5 /. Ω -> 1, λ]/λ], λ];
TraditionalForm[Reduce[Exists[k, k > 0 && discA5 == 0], Δ, Reals] /. Δ -> HoldForm[Δ/Ω]]
```

The top row lets $k$ be complex, so each exceptional point becomes an isolated zero of the discriminant $D(k)$, marked by a dot at the center of its rings of constant $\log|D|$. The bottom row plots, on the same axis for real $k$, the qubit's relaxation rate, the smallest nonzero $-\mathrm{Re}\,\lambda$ of $\mathcal L$, and its oscillation frequency, the largest $|\mathrm{Im}\,\lambda|$; dashed lines mark the real parts of the exceptional points, and $k$ and rates are in units of $\Omega$:

```wl
ΔsA5 = {0, 1/4, 1/(2 Sqrt[2]), 1};
epA5[Δv_] := Select[ReIm[k /. NSolve[(discA5 /. Δ -> Δv) == 0, k]], First[#] > 0 &];
markA5[Δv_, {y0_, y1_}] := {GrayLevel[0.4], AbsoluteThickness[1], Dashed, Line[{{#, y0}, {#, y1}}] & /@ Union[Round[epA5[Δv][[All, 1]], 10^-6]]};
labelA5[Δv_] := Inset[Framed[Style[Row[{"Δ/Ω = ", If[Δv === 1/(2 Sqrt[2]), Row[{"1/(2", TraditionalForm[Sqrt[2]], ")"}], Style[TraditionalForm[Δv], FractionBoxOptions -> {Beveled -> True}]]}], Black],
    Background -> White, FrameStyle -> None], Scaled[{0.97, 0.97}], {Right, Top}];
discPanelA5[Δv_] := ComplexContourPlot[Evaluate[Log[Abs[discA5 /. {Δ -> Δv, k -> z}]]], {z, 0.1 - 0.7 I, 1.7 + 0.7 I}, Contours -> 20, PlotPoints -> 60,
   Epilog -> {{GrayLevel[0.3], AbsoluteThickness[1], Line[{{0.1, 0}, {1.7, 0}}]}, markA5[Δv, {-0.7, 0.7}], {White, PointSize[0.045], Point[epA5[Δv]], Black, PointSize[0.028], Point[epA5[Δv]]}, labelA5[Δv]},
   FrameLabel -> {"Re k/Ω", "Im k/Ω"}, ImageSize -> 230];
rateA5[Δv_?NumericQ, kv_?NumericQ] := Min[Select[-Re[Eigenvalues[N[liouA5 /. {Ω -> 1, Δ -> Δv, k -> kv}]]], # > 10^-9 &]];
freqA5[Δv_?NumericQ, kv_?NumericQ] := Max[Abs[Im[Eigenvalues[N[liouA5 /. {Ω -> 1, Δ -> Δv, k -> kv}]]]]];
ratePanelA5[Δv_] := Plot[{rateA5[Δv, kk], freqA5[Δv, kk]}, {kk, 0.1, 1.7}, PlotPoints -> 200, PlotRange -> {-0.06, 1.5}, 
   Epilog -> markA5[Δv, {-0.06, 1.5}], GridLines -> {None, Prepend[Table[{y, Directive[GrayLevel[0.8], Dashed]}, {y, 0.2, 1.4, 0.2}], {0, Directive[GrayLevel[0.4]]}]}, AspectRatio -> 1/2, Frame -> True, FrameLabel -> {"k/Ω", "rate / Ω"}, ImageSize -> 230];
Labeled[Legended[Column[{PlotGrid[{discPanelA5 /@ ΔsA5}, ImageSize -> 900], PlotGrid[{ratePanelA5 /@ ΔsA5}, ImageSize -> 900]}, Spacings -> 0], Placed[LineLegend[ColorData[97] /@ {1, 2}, {"relaxation rate (slowest decay)", "oscillation frequency"}, LegendLayout -> {"Row", 1}], Below]],
  "top: log|D(k)| for complex k, with D the discriminant of the Liouvillian's characteristic polynomial;  bottom: relaxation rate and oscillation frequency for real k", Top]
```

On resonance the slowest decay rate rises to a sharp peak at $k=\Omega$, directly below the exceptional point: the qubit relaxes fastest there, critically damped between Rabi oscillation and the Zeno slowdown, whose rate falls as $\Omega^2/2k$. A small detuning splits the exceptional point into two zeros on the real axis, one where the Rabi pair turns real and one where two real rates meet and turn complex again; at $|\Delta/\Omega|=1/(2\sqrt2)$ the two zeros merge, and for $|\Delta/\Omega|>1/(2\sqrt2)$ they leave the axis as a complex-conjugate pair, so no measurement strength reaches them, the Rabi pair never turns real, and the sharp peak of the decay rate becomes a smooth hump.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A charge qubit read by a point contact: when the Rabi peaks merge

One electron tunnels between two quantum dots, left $|L\rangle=|0\rangle$ and right $|R\rangle=|1\rangle$, with tunnelling $\Omega$ and level offset $\varepsilon$, while a nearby quantum point contact passes a current that depends on which dot is occupied. The current noise carries the fluctuations of $\hat\sigma_z=|L\rangle\langle L|-|R\rangle\langle R|$. As the contact reads more strongly, do the Rabi peaks of that noise spectrum merge where the Liouvillian's eigenvectors coalesce?

With $\hbar=1$, $\hat n_L=|L\rangle\langle L|$ and $\mathcal D[\hat c]\rho=\hat c\rho \hat c^\dagger-\frac12\{\hat c^\dagger \hat c,\rho\}$,

$$
\dot\rho=\mathcal L\rho=-i\left[\frac\varepsilon2\hat\sigma_z+\frac\Omega2\hat\sigma_x,\rho\right]+\kappa\,\mathcal D[\hat n_L]\rho ,\qquad k=\kappa/4 ,
$$

since $\mathcal D[\hat n_L]=\frac14\mathcal D[\hat\sigma_z]$, the watched-qubit Liouvillian of the Zeno crossover with $\varepsilon$ for $\Delta$. By quantum regression, the symmetrized spectrum is $S(\omega)=2\,\mathrm{Re}\,\mathrm{Tr}\left[\hat\sigma_z\,(i\omega-\mathcal L)^{-1}\frac12\{\hat\sigma_z-\langle\hat\sigma_z\rangle_{\rm ss},\rho_{\rm ss}\}\right]$ with steady state $\rho_{\rm ss}$. A side peak exists when $dS/d\omega=0$ has a root at $\omega>0$; an exceptional point is a $k$ where two nonzero eigenvalues of $\mathcal L$ coincide.

Inverting $\mathcal L$ for general $\varepsilon$, $\Omega$ and $k$ gives $S(\omega)$ in closed form, with the largest $k$ that keeps side peaks and the $k$ where the discriminant of $\mathcal L$'s characteristic polynomial vanishes; on resonance these are:

```wl
qpcL = Normal[QuantumOperator["Liouvillian"[ε/2 QuantumOperator["Z"] + Ω/2 QuantumOperator["X"], {QuantumOperator[{{1, 0}, {0, 0}}]}, {4 k}]]["Matrix"]];
qpcSS = With[{v = ArrayReshape[First[NullSpace[qpcL]], {2, 2}]}, v/Tr[v]];
qpcSpec = With[{dz = PauliMatrix[3] - Tr[PauliMatrix[3] . qpcSS] IdentityMatrix[2]},
    Simplify[2 ComplexExpand[Re[Tr[PauliMatrix[3] . ArrayReshape[LinearSolve[I ω IdentityMatrix[4] - qpcL, Flatten[(dz . qpcSS + qpcSS . dz)/2]], {2, 2}]]]]]];
qpcSideEdge[εv_] := MaxValue[{k, Reduce[Exists[u, u > 0 && Numerator[Together[D[qpcSpec /. {ω -> Sqrt[u], Ω -> 1, ε -> εv}, u]]] == 0] && k > 0, k]}, k];
qpcDisc = Discriminant[Cancel[CharacteristicPolynomial[qpcL, λ]/λ], λ];
qpcEP[εv_] := Cases[Solve[(qpcDisc /. {Ω -> 1, ε -> εv}) == 0 && k > 0, k, Reals], {k -> v_} :> v];
<|"S(ω)" -> Simplify[qpcSpec /. ε -> 0],
  "side peaks" -> Simplify[Reduce[SeriesCoefficient[qpcSpec /. ε -> 0, {ω, 0, 2}] > 0 && k > 0 && Ω > 0, k], Ω > 0],
  "exceptional point" -> k == Normal[k /. First[Solve[(qpcDisc /. ε -> 0) == 0 && k > 0 && Ω > 0, k]]]|>
```

Each panel plots $S(\omega)$ for three measurement strengths, then for the largest $k$ that keeps side peaks and, where one exists, for the $k$ of each exceptional point; the gray lines mark the bare Rabi frequencies $\pm\sqrt{\Omega^2+\varepsilon^2}$, and frequencies and $k$ are in units of $\Omega$:

```wl
qpcKs = {0.25, 0.5, 1.5};
qpcPanel[εv_] := With[{ks = Join[qpcKs, {N[qpcSideEdge[Rationalize[εv]]]}, N[qpcEP[Rationalize[εv]]]], w = Sqrt[1 + εv^2]},
   Plot[Evaluate[qpcSpec /. {Ω -> 1, ε -> εv, k -> #} & /@ ks], {ω, -3, 3}, PlotRange -> All,
    GridLines -> {{-w, w}, None}, GridLinesStyle -> Directive[Gray, Dashed], Frame -> True, FrameLabel -> {"ω/Ω", "Ω S(ω)"},
    PlotLegends -> Placed[LineLegend[Row[{"k/Ω = ", Round[#, 0.01] /. x_Real /; x == Round[x] :> Round[x]}] & /@ ks, LegendLayout -> {"Row", 2}], Below],
    PlotLabel -> Row[{"ε/Ω = ", εv}], ImageSize -> 300]];
PlotGrid[{qpcPanel /@ {0, 0.3, 1}}, ImageSize -> 900]
```

On resonance the two peaks slide inward and meet at zero frequency at $k=\Omega/\sqrt2$, while the eigenvectors of $\mathcal L$ coalesce only at the exceptional point $k=\Omega$, where the spectrum has long shown one central peak. An offset raises the center of the spectrum and removes the side peaks at a smaller $k$; at $\varepsilon/\Omega=1$ they vanish although no exceptional point exists at all, so the merging of the peaks is a property of the line shape, not of eigenvector coalescence.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Dispersive readout pointer states

A readout cavity of linewidth $\kappa$, coupled to a transmon, has its frequency pulled by $+\chi$ or $-\chi$ depending on the qubit state, and is driven by a readout pulse $\varepsilon_r(t)$ of peak amplitude $\epsilon_{\max}$. How does the cavity field separate for the two qubit states as $\kappa/2|\chi|$ changes?

The circuit and its Hamiltonian at three levels of approximation; hovering over a part of the circuit highlights its term, and Level 2 is the dispersive Hamiltonian used below:

```wl
TransmonCircuitPlot[]
```

With $\hbar=1$, in the frame of the readout drive, the cavity holds a coherent state $|\alpha_j\rangle$ for qubit state $j$, with

$$
\dot\alpha_{0}=-\left(\frac{\kappa}{2}+i\chi\right)\alpha_0-\frac{i}{2}\epsilon_{\max}\varepsilon_r(t),\qquad \dot\alpha_{1}=-\left(\frac{\kappa}{2}-i\chi\right)\alpha_1-\frac{i}{2}\epsilon_{\max}\varepsilon_r(t),
$$

so under a constant drive $\alpha_j$ settles at $-i\epsilon_{\max}/(\kappa\pm2i\chi)$. The field leaking out of the cavity, which the measurement records, is $\sqrt\kappa\,\alpha_j$.

For a 4 μs flat-top readout the figure plots, up to the fall of the drive, the intracavity field $\alpha_j(t)$ with its steady states as points, and the output separation $\sqrt\kappa\,|\alpha_0-\alpha_1|$ in units of $\sqrt{2|\chi|}$. The transmon has qubit-cavity coupling 50 MHz and $T_1=T_2=10$ μs, so $\chi$ comes from the dispersive formula and the drive is $\epsilon_{\max}\approx9|\chi|$:

```wl
pointerSystem[κHz_] := TransmonSystem["Coupling" -> 50*^6, "CavityLinewidth" -> κHz, "T1" -> 10.*^-6, "T2" -> 10.*^-6];
pointerChiHz = Abs[TransmonSystem["Coupling" -> 50*^6]["Chi"]];
pointerRs = {0.3, 0.6, 2};
pointerDR = AssociationMap[DispersiveReadout[pointerSystem[2 # pointerChiHz], Pulse["FlatTop", 4000, "RiseTime" -> 50]] &, pointerRs];
pointerEnd = 4000*^-9 - 50*^-9;
pointerSep[r_, x_] := Sqrt[r] Abs[pointerDR[r]["Alpha0"][x/(2 Pi pointerChiHz)] - pointerDR[r]["Alpha1"][x/(2 Pi pointerChiHz)]];
pointerXEnd = 2 Pi pointerChiHz pointerEnd;
pointerIQ = Flatten[Table[ReIm[{pointerDR[r]["Alpha0"][s], pointerDR[r]["Alpha1"][s]}], {r, pointerRs}, {s, 0, pointerEnd, pointerEnd/400}], 2];
pointerRange = {{Min[#], Max[#]} &[pointerIQ[[All, 1]]], {Min[#], Max[#]} &[pointerIQ[[All, 2]]]};
pointerSepMax = 1.05 Max[Table[pointerSep[r, x], {r, pointerRs}, {x, 0, pointerXEnd, pointerXEnd/400}]];
pointerPanel[r_] := With[{dr = pointerDR[r]},
   Show[ParametricPlot[Evaluate[ReIm[{dr["Alpha0"][s], dr["Alpha1"][s]}]], {s, 0, pointerEnd},
     PlotLegends -> Placed[LineLegend[{"qubit in |0⟩", "qubit in |1⟩"}, LegendMarkerSize -> 15, LegendLayout -> "Row"], {{0.5, 0.02}, {0.5, 0}}]], ListPlot[List /@ ReIm[Values[dr["SteadyState"]]], PlotStyle -> PointSize[Large]],
    GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.8], Dashed], PlotRange -> pointerRange, PlotRangePadding -> {Scaled[0.04], {Scaled[0.14], Scaled[0.04]}}, Axes -> False, AspectRatio -> 1, Frame -> True, FrameLabel -> {"Re α", "Im α"}, ImagePadding -> {{60, 10}, {45, 10}}, ImageSize -> 300]];
pointerOutput[r_] := Plot[pointerSep[r, x], {x, 0, pointerXEnd}, PlotStyle -> ColorData[97][3], PlotRange -> {0, pointerSepMax}, GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.8], Dashed], AspectRatio -> 1, Frame -> True,
   FrameLabel -> {"|χ| t", Row[{"output separation / ", Sqrt[Row[{2, "|χ|"}]]}]}, ImagePadding -> {{60, 10}, {45, 10}}, ImageSize -> 300];
pointerPanels = AssociationMap[{pointerPanel[#], pointerOutput[#]} &, pointerRs];
Manipulate[Grid[{pointerPanels[r]}], Row[{Control[{{r, 0.6, "cavity linewidth over dispersive shift κ/2|χ|"}, pointerRs, SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"κ/2|χ| = ", r}], <|0.3 -> "Narrow cavity: long spirals, steady states far apart inside.",
        0.6 -> "The output separation ends largest.",
        2 -> "Wide cavity: the output separation settles soonest."|>[r]}}]]}], SaveDefinitions -> True]
```

The cavity rings up along a different spiral for each qubit state and settles on its steady state; a narrower cavity spirals more and ends with the two states farther apart inside. The output separation is not ordered the same way: it settles soonest in the widest cavity, but ends largest at the intermediate linewidth.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Born's rule as a martingale: runs collapse to the poles with probability $(1\pm z_0)/2$

A qubit, Rabi-driven at frequency $\Omega$, starts in the pure state with $z=z_0$, and a perfect detector measures $\hat\sigma_z$ continuously at strength $k$. Without a drive, each run ends at a pole; with what probability does it end at $z=+1$?

With $\hbar=1$ and $\hat H=\frac{\Omega}{2}\hat\sigma_x$, the $z$ component of the conditioned state obeys the Itô equation $dz=\Omega y\,dt+2\sqrt{k}\,(1-z^2)\,dW$, and the detector records $dY=2\sqrt k\,z\,dt+dW$. In each step the record multiplies the amplitudes of $|0\rangle$ and $|1\rangle$ by $e^{\pm\sqrt k\,dY}$ and the drive rotates them.

The record is sampled as 400 runs for each start. The left panel takes seven starting values $z_0$ without a drive and plots the fraction of runs that end above $z=0$, with its standard error, against $(1+z_0)/2$; the right panel follows four runs and the mean of all 400 for the chosen start and drive. Time is in units of $1/k$:

```wl
bornRunsZ[Ω_, z0_, dt_ : 0.005] := With[{u = MatrixExp[-I Ω dt/2 PauliMatrix[1]]},
   Transpose[(Abs[#]^2 . {1, -1}) & /@ FoldList[
     Function[{amp, dw}, Normalize /@ ((amp Exp[Outer[Times, 2 (Abs[amp]^2 . {1, -1}) dt + dw, {1, -1}]]) . u)],
     ConstantArray[N[Sqrt[{1 + z0, 1 - z0}/2]], 400], RandomVariate[NormalDistribution[0, Sqrt[dt]], {800, 400}]]]];
bornStarts = {-0.5, 0.2, 0.8}; bornDrives = {0, 1};
bornFractions = BlockRandom[SeedRandom[7]; Function[z0, With[{f = Mean[UnitStep[bornRunsZ[0, z0][[All, -1]]]]}, {z0, Around[f, Sqrt[f (1 - f)/400]]}]] /@ {-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8}];
bornLawPanel[z0_] := Legended[Show[
     Plot[(1 + s)/2, {s, -1, 1}, PlotStyle -> Directive[Red, AbsoluteThickness[2]]],
     ListPlot[bornFractions, PlotMarkers -> {"\[FilledCircle]", 8}, PlotStyle -> ColorData[97][1]],
     PlotRange -> {{-1, 1}, {0, 1}}, GridLines -> {{z0}, None}, GridLinesStyle -> Directive[Gray, Dashed], Frame -> True,
     FrameLabel -> {Subscript["z", 0], "fraction ending at +1"}, PlotLabel -> "no drive", AspectRatio -> 1, ImageSize -> 250],
    Placed[Column[{PointLegend[{ColorData[97][1]}, {"fraction of 400 runs"}, LegendMarkers -> {"\[FilledCircle]", 8}],
       LineLegend[{Directive[Red, AbsoluteThickness[2]]}, {Row[{"(1 + ", Subscript["z", 0], ")/2"}]}]}], Below]];
bornRunStyle = Directive[GrayLevel[0.6], Opacity[0.7], AbsoluteThickness[1]];
bornRunsPanel[{Ω_, z0_}] := With[{zs = bornRunsZ[Ω, z0]}, Legended[Show[
     ListLinePlot[Take[zs, 4], DataRange -> {0, 4}, PlotStyle -> bornRunStyle],
     ListPlot[Mean[zs][[;; ;; 20]], DataRange -> {0, 4}, PlotMarkers -> {"\[FilledCircle]", 6}, PlotStyle -> ColorData[97][1]],
     PlotRange -> {-1.05, 1.05}, GridLines -> {None, {z0}}, GridLinesStyle -> Directive[Red, Dashed], Frame -> True,
     FrameLabel -> {"k t", Row[{"\[LeftAngleBracket]", Subscript["\[Sigma]", "z"], "\[RightAngleBracket]"}]}, AspectRatio -> 1, ImageSize -> 250],
    Placed[Column[{LineLegend[{bornRunStyle}, {"single runs"}], PointLegend[{ColorData[97][1]}, {"mean of 400 runs"}, LegendMarkers -> {"\[FilledCircle]", 6}],
       LineLegend[{Directive[Red, Dashed]}, {Subscript["z", 0]}]}], Below]]];
bornRuns = BlockRandom[SeedRandom[8]; AssociationMap[bornRunsPanel, Tuples[{bornDrives, bornStarts}]]];
bornLaw = AssociationMap[bornLawPanel, bornStarts];
Manipulate[Grid[{{bornLaw[bornZ0], bornRuns[{bornOmega, bornZ0}]}}, Alignment -> Top, Spacings -> 3],
  Row[{Column[{Control[{{bornZ0, 0.2, Subscript["z", 0]}, bornStarts, SetterBar}], Control[{{bornOmega, 0, "\[CapitalOmega]/k"}, bornDrives, SetterBar}]}], Spacer[20], Dynamic[controlNotes[{
      {Row[{Subscript["z", 0], " = ", bornZ0}], <|-0.5 -> "Undriven, a quarter of the runs end at z = +1.",
        0.2 -> "Undriven, 60% of the runs end at z = +1.",
        0.8 -> "Undriven, 90% of the runs end at z = +1."|>[bornZ0]},
      {Row[{"Ω/k = ", bornOmega}], <|0 -> "No drive: the mean stays at z₀ and every run settles on a pole.",
        1 -> "Drive: the mean leaves z₀ and runs keep jumping."|>[bornOmega]}}]]}], SaveDefinitions -> True]
```

Without a drive the drift of $z$ is zero: every run heads for a pole while the mean stays on $z_0$. A bounded process with constant mean that ends at $\pm1$ must end at $+1$ with probability $(1+z_0)/2$, and the fraction of runs ending above zero matches it within its standard error. A drive adds the drift $\Omega y$: the mean leaves $z_0$ and the runs keep jumping instead of settling.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Quantum Bayesian update: each reading divides the mixedness by $w^2$

A detector reads a Pauli observable $\hat A$ of a qubit over a short window and returns a noisy reading $I$. Two detectors take turns, one reading $\hat\sigma_z$ and one reading $\hat\sigma_\vartheta=\cos\vartheta\,\hat\sigma_z+\sin\vartheta\,\hat\sigma_x$. How does the mixedness $1-r^2$ of the state change with each reading?

Each reading updates the state by Bayes' rule, written as a Kraus operator (Korotkov):

$$
\rho\;\longrightarrow\;\frac{\hat M\rho\hat M}{w},\qquad \hat M(I)=e^{I\hat A/2D},\qquad w=\mathrm{Tr}\,\hat M\rho\hat M ,
$$

and $I$ is drawn from a Gaussian of variance $D$ centered on $+1$ or $-1$ with the Born probabilities of $\hat A$. $\det\hat M$, for every reading and every $\vartheta$, is:

```wl
krausGauss[i_, A_, d_] := MatrixExp[i/(2 d) A];
detKraus = FullSimplify[Det[krausGauss[iRead, Cos[ϑ] PauliMatrix[3] + Sin[ϑ] PauliMatrix[1], dVar]], {iRead, ϑ} ∈ Reals && dVar > 0]
```

A run of 20 readings starts from $(x,y,z)=(0.3,0,0.4)$, alternating the two detectors. The left panel takes every reading of ten runs at each $D$ and $\vartheta$ and plots how much it changed the mixedness $1-r^2=4\det\rho$ against $1/w^2$; the right panel shows five runs at each detector noise $D$:

```wl
bayesStep[{ρ_, _}, {A_, d_}] := With[{p = Re[Tr[ρ . (IdentityMatrix[2] + A)/2]]},
   {mk = krausGauss[RandomVariate[NormalDistribution[RandomChoice[{p, 1 - p} -> {1, -1}], Sqrt[d]]], A, d]}, {w = Re[Tr[mk . ρ . mk]]}, {mk . ρ . mk/w, w}];
bayesMixedness[ρ_] := 4 Re[ρ[[1, 1]] ρ[[2, 2]] - ρ[[1, 2]] ρ[[2, 1]]];
bayesRho0 = (IdentityMatrix[2] + 0.3 PauliMatrix[1] + 0.4 PauliMatrix[3])/2;
bayesRun[d_, ϑ_] := FoldList[bayesStep[#1, {#2, d}] &, {bayesRho0, 1.}, Flatten[Table[{PauliMatrix[3], Cos[ϑ] PauliMatrix[3] + Sin[ϑ] PauliMatrix[1]}, 10], 1]];
bayesSteps[run_] := MapThread[{1/Last[#2]^2, bayesMixedness[First[#2]]/bayesMixedness[First[#1]]} &, {Most[run], Rest[run]}];
bayesThetas = {0, π/4, π/2};
bayesLawPanel = With[{pts = BlockRandom[SeedRandom[6]; Table[Catenate[Catenate[Table[bayesSteps[bayesRun[d, ϑ]], {d, {2, 8}}, {10}]]], {ϑ, bayesThetas}]]},
   Legended[Show[
     ListLogLogPlot[pts, PlotMarkers -> {{"\[FilledCircle]", 8}, {"\[EmptySquare]", 9}, {"\[Times]", 9}}, PlotStyle -> ColorData[97] /@ Range[3]],
     LogLogPlot[x, {x, 0.05, 50}, PlotStyle -> Directive[Gray, Dashed]],
     PlotRange -> All, Frame -> True, FrameLabel -> {Row[{"1/", Superscript["w", 2]}], Row[{"det ", Subscript["\[Rho]", "after"], " / det ", Subscript["\[Rho]", "before"]}]},
     PlotLabel -> "every reading", GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.8], Dashed], AspectRatio -> 1, ImageSize -> 250],
    Placed[PointLegend[ColorData[97] /@ Range[3], Row[{"ϑ = ", #}] & /@ {0, "π/4", "π/2"}, LegendMarkers -> {{"\[FilledCircle]", 8}, {"\[EmptySquare]", 9}, {"\[Times]", 9}}, LegendLayout -> "Row"], Below]]];
bayesNoisePanel = With[{runs = BlockRandom[SeedRandom[7]; Table[bayesMixedness[First[#]] & /@ bayesRun[d, π/4], {d, {2, 8}}, {5}]]}, Legended[Show[
     MapThread[ListLogPlot[#1, Joined -> True, DataRange -> {0, 20}, PlotStyle -> Directive[#2, AbsoluteThickness[1.5]]] &, {runs, ColorData[97] /@ {4, 5}}],
     PlotRange -> All, Frame -> True, FrameLabel -> {"reading number", Row[{"1 \[Minus] ", Superscript["r", 2]}]}, PlotLabel -> "five runs each, ϑ = π/4", GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.8], Dashed], AspectRatio -> 1, ImageSize -> 250],
    Placed[LineLegend[ColorData[97] /@ {4, 5}, {"D = 2", "D = 8"}, LegendLayout -> "Row"], Below]]];
Grid[{{bayesLawPanel, bayesNoisePanel}}, Alignment -> Top, Spacings -> 3]
```

Every Kraus operator has unit determinant, for every reading and every $\vartheta$, so each update multiplies $\det\rho$, and with it the mixedness, by $1/w^2$: every reading lands on the diagonal, the unlikely ones that raise the mixedness (above 1) as well as the likely ones that lower it. At $D=2$ the mixedness falls by orders of magnitude along every run. At $D=8$ one reading carries little information, and twenty readings purify the state only partly.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Measurement purifies a qubit: Itô's lemma in spherical coordinates

A qubit with Hamiltonian $\hat H=\frac{\Delta}{2}\hat\sigma_z+\frac{\Omega}{2}\hat\sigma_x$ is measured continuously along $\hat\sigma_z$ at strength $k$ by a detector of efficiency $\eta$. Along one run, does watching it purify its state, as measured by the length $r$ of its Bloch vector, and what does an imperfect detector take away?

With $\hbar=1$, the state conditioned on the record follows the Itô stochastic master equation of a $\hat\sigma_z$ measurement: $\hat H$ rotates it, the term $k\,\mathcal D[\hat\sigma_z]\rho$ dephases it, and the record kicks it through $\sqrt{\eta k}\left(\hat\sigma_z\rho+\rho\hat\sigma_z-2\langle\hat\sigma_z\rangle\rho\right)dW$.

As an Itô process for the Bloch vector, Itô's lemma for $r$ in spherical coordinates $(r,\theta,\varphi)$, with $z=r\cos\theta$, gives the stochastic equation for $r$:

```wl
blochOfPurify[m_] := Simplify[Tr[m . #] & /@ PauliMatrix[{1, 2, 3}]];
blochSME[Δ_, Ω_, k_, η_, v0_] := With[{ρv = (IdentityMatrix[2] + {x[t], y[t], z[t]} . PauliMatrix[{1, 2, 3}])/2, hv = Δ/2 PauliMatrix[3] + Ω/2 PauliMatrix[1], sz = PauliMatrix[3]},
   ItoProcess[{blochOfPurify[-I (hv . ρv - ρv . hv) + k (sz . ρv . sz - ρv)],
     List /@ blochOfPurify[Sqrt[η k] (sz . ρv + ρv . sz - 2 z[t] ρv)]}, {{x, y, z}, v0}, t]];
drSph = FullSimplify[blochSME[Δ, Ω, k, η, {x0, y0, z0}]["Derivative", Sqrt[x[t]^2 + y[t]^2 + z[t]^2]] /.
    Thread[{x[t], y[t], z[t]} -> r {Sin[θ] Cos[φ], Sin[θ] Sin[φ], Cos[θ]}], r > 0 && 0 < θ < π && k > 0 && 0 <= η <= 1 && {Δ, Ω} ∈ Reals];
Row[{"dr = ", TraditionalForm[drSph[[1, 1]]], " dt + ", TraditionalForm[drSph[[2, 1, 1]] /. -1 + r^2 -> -HoldForm[1 - r^2]], " dW"}]
```

For each detector efficiency $\eta$, three runs start from $(x,y,z)=(0.3,0,0.4)$ under a resonant drive of strength $\Omega/k=0$, $1$ and $5$; the dashed line marks $r=\sqrt\eta$, where the drift changes sign, and time is in units of $1/k$:

```wl
panelPurify[ηv_, {Ωv_, Δv_}, label_] := ListLinePlot[Map[{#[[1]], Norm[#[[2]]]} &,
     RandomFunction[blochSME[Δv, Ωv, 1, ηv, {0.3, 0, 0.4}], {0., 4., 0.001}, 3, Method -> "StochasticRungeKuttaScalarNoise"]["Paths"], {2}],
    PlotRange -> {0, 1.05}, GridLines -> {Table[{x, Directive[GrayLevel[0.85], Dashed]}, {x, 1, 4}], Append[Table[{y, Directive[GrayLevel[0.85], Dashed]}, {y, 0.2, 1, 0.2}], {Sqrt[ηv], Directive[Gray, Dashed, AbsoluteThickness[1.5]]}]},
    Frame -> True, FrameLabel -> {"k t", "r"}, PlotLabel -> label, ImageSize -> 300];
purifyDrives = {{{0, 0}, "\[CapitalOmega]/k = 0"}, {{1, 0}, "\[CapitalOmega]/k = 1"}, {{5, 0}, "\[CapitalOmega]/k = 5"}};
purifyPanels = BlockRandom[SeedRandom[2]; AssociationMap[Function[ηv, panelPurify[ηv, ##] & @@@ purifyDrives], {0.5, 0.8, 1}]];
Manipulate[Legended[PlotGrid[{purifyPanels[purifyEta]}, ImageSize -> 800],
   Placed[LineLegend[Append[ColorData[97] /@ Range[3], Directive[Gray, Dashed]], {"run 1", "run 2", "run 3", TraditionalForm[HoldForm[r == Sqrt[η]]]}, LegendLayout -> "Row"], Below]],
  Row[{Control[{{purifyEta, 1, "\[Eta]"}, Keys[purifyPanels], SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"η = ", purifyEta}], <|0.5 -> "Drift turns inward above r = √η ≈ 0.71; runs purify only at a pole.",
        0.8 -> "Drift turns inward above r = √η ≈ 0.89; a strong drive keeps runs mixed.",
        1 -> "Every run climbs to r = 1 and stays pure."|>[purifyEta]}}]]}], SaveDefinitions -> True]
```

Neither term contains $\Delta$ or $\Omega$: the Hamiltonian only rotates the Bloch vector. The drift is positive inside $r=\sqrt\eta$ and the noise vanishes at $r=1$, so at $\eta=1$ every run climbs to the surface and stays there: a mixed state purifies and a pure state stays pure. At $\eta<1$ the drift at $r=1$ is negative except at the poles, so without a drive the runs purify by reaching a pole, a drive comparable to $k$ lets them reach the poles only part of the time, and a strong drive keeps them off the poles and below $r=1$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A monitored qubit stays in the Bloch ball: Itô's lemma on the sphere

A qubit is Rabi-driven at frequency $\Omega$, relaxes toward $z=1$ at rate $\gamma_1$, dephases at rate $\gamma_\phi$, and is read out along $\hat\sigma_z$ by a detector that dephases it at rate $\Gamma_d$ and whose noise enters $z$ at the informational rate $\Gamma_{CI}$ and rotates the phase at the back-action rate $\Gamma_{BA}$; its efficiency is $\eta=(\Gamma_{CI}+\Gamma_{BA})/(2\Gamma_d)$. Does every run stay a density matrix, $r=|\vec r|\le1$, for every choice of rates?

With $\hbar=1$ and one Wiener noise $W$, the Bloch vector obeys the Itô equations

$$
\begin{aligned}
dx&=-\left(\frac{\gamma_1}{2}+\gamma_\phi+\Gamma_d\right)x\,dt-\left(\sqrt{\Gamma_{CI}}\,xz+\sqrt{\Gamma_{BA}}\,y\right)dW,\\
dy&=\left[-\left(\frac{\gamma_1}{2}+\gamma_\phi+\Gamma_d\right)y-\Omega z\right]dt-\left(\sqrt{\Gamma_{CI}}\,yz-\sqrt{\Gamma_{BA}}\,x\right)dW,\\
dz&=\left[\gamma_1(1-z)+\Omega y\right]dt+\sqrt{\Gamma_{CI}}\,(1-z^2)\,dW .
\end{aligned}
$$

Itô's lemma for $r^2$ on the sphere $r=1$, where $(x,y)=\sqrt{1-z^2}(\cos\varphi,\sin\varphi)$, gives the noise of $r^2$ there and the efficiencies for which its drift points outward somewhere:

```wl
monitoredQubit = ItoProcess[{
    {-(γ1/2 + γφ + Γd) x[t], -(γ1/2 + γφ + Γd) y[t] - Ω z[t], γ1 (1 - z[t]) + Ω y[t]},
    {{-Sqrt[ΓCI] x[t] z[t] - Sqrt[ΓBA] y[t]}, {-Sqrt[ΓCI] y[t] z[t] + Sqrt[ΓBA] x[t]}, {Sqrt[ΓCI] (1 - z[t]^2)}}}, {{x, y, z}, {x0, y0, z0}}, t];
{r2Drift, r2Noise} = Factor[Simplify[Flatten[monitoredQubit["Derivative", x[t]^2 + y[t]^2 + z[t]^2]] /.
     {x[t] -> Sqrt[1 - z^2] Cos[φ], y[t] -> Sqrt[1 - z^2] Sin[φ], z[t] -> z, ΓBA -> 2 η Γd - ΓCI}]];
ballOutwardCondition = Simplify[Reduce[Exists[z, -1 < z < 1, r2Drift > 0] && Γd > 0 && γφ >= 0 && γ1 >= 0, η, Reals], Γd > 0 && γφ >= 0 && γ1 >= 0];
TraditionalForm[<|"noise of r² on the sphere" -> r2Noise, "drift outward somewhere" -> ballOutwardCondition|> /.
   {Γd -> Subscript[Γ, "d"], γφ -> Subscript[γ, "φ"], γ1 -> Subscript[γ, 1]}]
```

Three runs start from the mixed state $(x,y,z)=(0.3,0,0.3)$, with $\Omega=2\Gamma_d$, $\gamma_1=0.2\,\Gamma_d$, no extra dephasing, $\gamma_\phi=0$, and the detector noise split equally, $\Gamma_{CI}=\Gamma_{BA}=\eta\Gamma_d$. Each panel plots the Bloch vector length $r$ for one efficiency, with the region $r>1$, where $\rho$ is not a density matrix, shaded; time is in units of $1/\Gamma_d$:

```wl
ballEtas = {0.5, 1, 3};
ballRun[ηv_] := ListLinePlot[Map[{#[[1]], Norm[#[[2]]]} &,
     RandomFunction[monitoredQubit /. {γ1 -> 0.2, γφ -> 0, Γd -> 1, Ω -> 2, ΓCI -> ηv, ΓBA -> ηv, x0 -> 0.3, y0 -> 0, z0 -> 0.3}, {0., 5., 0.001}, 3,
       Method -> "StochasticRungeKuttaScalarNoise"]["Paths"][[All, ;; ;; 5]], {2}],
    PlotRange -> {0, 1.6}, Prolog -> {GrayLevel[0.88], Rectangle[{0, 1}, {5, 1.6}]}, GridLines -> {None, {1}},
    Epilog -> {Inset[Style["not a density matrix", Gray, 11], Scaled[{0.97, 0.95}], {Right, Top}], Inset[Row[{"\[Eta] = ", ηv}], Scaled[{0.97, 0.06}], {Right, Bottom}]},
    Frame -> True, FrameLabel -> {Row[{Subscript["\[CapitalGamma]", "d"], " t"}], "Bloch vector length r"}, ImageSize -> 280];
BlockRandom[SeedRandom[3]; PlotGrid[{ballRun /@ ballEtas}, ImageSize -> 800]]
```

The measurement does two things to the Bloch vector: it dephases the qubit at rate $\Gamma_d$, which pulls the state inward, and each update from the record pushes it toward purity, at a rate set by $\eta\Gamma_d$. On the sphere the kicks only slide the state along it, so the drift decides, and it points outward somewhere exactly when the information outruns the disturbance, $\eta\Gamma_d>\Gamma_d+\gamma_\phi$, whatever the drive and the relaxation. A real detector has $\eta\le1$, so the runs at $\eta=0.5$ and $1$ stay in the ball, while $\eta=3$, a detector that would learn more than it disturbs, drives them out. Extra dephasing $\gamma_\phi$ disturbs without informing, so it only widens the margin.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### When the readout record fixes the population

A qubit, Rabi-driven at frequency $\Omega$ and relaxing toward $z=1$ at rate $\gamma_1$, is read out along $\hat\sigma_z$; the detector dephases it at rate $\Gamma_d$, its noise moves $z$ at rate $\Gamma_{CI}$ and the phase at rate $\Gamma_{BA}$, and $Q$ is the integrated readout. Is $z(t)=\tanh\left(\operatorname{arctanh}z_0+\sqrt{\Gamma_{CI}}\,Q(t)\right)$, with $z_0=z(0)$, on every run?

With $\hbar=1$, the population and the readout obey the Itô equations, driven by one Wiener noise $W$,

$$
\begin{aligned}
dz&=\left[\gamma_1(1-z)+\Omega y\right]dt+\sqrt{\Gamma_{CI}}\,(1-z^2)\,dW,\\
dQ&=\sqrt{\Gamma_{CI}}\,z\,dt+dW ,
\end{aligned}
$$

and $x$, $y$ follow the Bloch equations. The Itô derivative of $\operatorname{arctanh}z-\sqrt{\Gamma_{CI}}\,Q$, for general rates, is:

```wl
readoutV0 = {0.8, 0, -0.6, 0};
readoutIto = ItoProcess[{
    {-(γ1/2 + Γd) x[t], -(γ1/2 + Γd) y[t] - Ω z[t], γ1 (1 - z[t]) + Ω y[t], Sqrt[ΓCI] z[t]},
    {{-Sqrt[ΓCI] x[t] z[t] - Sqrt[ΓBA] y[t]}, {-Sqrt[ΓCI] y[t] z[t] + Sqrt[ΓBA] x[t]}, {Sqrt[ΓCI] (1 - z[t]^2)}, {1}}},
   {{x, y, z, q}, readoutV0}, t];
readoutStrat = StratonovichProcess[readoutIto];
readoutGap = Simplify[readoutIto["Derivative", ArcTanh[z[t]] - Sqrt[ΓCI] q[t]] /. {x[t] -> x, y[t] -> y, z[t] -> z}];
Row[{"d(arctanh z \[Minus] ", TraditionalForm[Sqrt[Subscript["\[CapitalGamma]", "CI"]]], " Q) = (",
   TraditionalForm[Collect[readoutGap[[1, 1]], {γ1, Ω}, Simplify] /. {γ1 -> Subscript["\[Gamma]", 1], Γd -> Subscript["\[CapitalGamma]", "d"]}], ") dt + ",
   TraditionalForm[readoutGap[[2, 1, 1]]], " dW"}]
```

Integration along a smooth path converges to the Stratonovich solution, so all four equations are put in Stratonovich form and integrated on one cubic interpolant of a sampled Wiener path, the same in every panel, from $(x,y,z)=(0.8,0,-0.6)$ with $\Gamma_{CI}/\Gamma_d=\Gamma_{BA}/\Gamma_d=0.5$. Each panel compares $z(t)$ with the formula built from the record; time is in units of $1/\Gamma_d$:

```wl
smoothWiener = BlockRandom[SeedRandom[5]; Interpolation[Transpose[{Subdivide[0., 4., 400],
      Prepend[Accumulate[RandomVariate[NormalDistribution[0, Sqrt[0.01]], 400]], 0.]}], InterpolationOrder -> 3]];
panelTanh[Ωv_, g1_] := With[{rates = {Γd -> 1, ΓCI -> 0.5, ΓBA -> 0.5, Ω -> Ωv, γ1 -> g1}},
   With[{sol = NDSolveValue[{Thread[D[{x[t], y[t], z[t], q[t]}, t] == (readoutStrat["Drift"] + Flatten[readoutStrat["Diffusion"]] smoothWiener'[t] /. rates)],
        Thread[{x[0], y[0], z[0], q[0]} == readoutV0]}, {z, q}, {t, 0, 4}, MaxStepSize -> 0.0025]},
    Plot[Evaluate[{sol[[1]][s], Tanh[ArcTanh[readoutV0[[3]]] + Sqrt[ΓCI] sol[[2]][s]]} /. rates], {s, 0, 4}, PlotStyle -> {Automatic, Dashed}, PlotRange -> {-1, 1},
     Frame -> True, FrameLabel -> {Row[{Subscript["\[CapitalGamma]", "d"], " t"}], "z"},
     PlotLabel -> Row[{"\[CapitalOmega]/", Subscript["\[CapitalGamma]", "d"], " = ", Ωv, ",   ", Subscript["\[Gamma]", 1], "/", Subscript["\[CapitalGamma]", "d"], " = ", g1}], ImageSize -> 300]]];
Legended[PlotGrid[Outer[panelTanh, {0, 1}, {0, 0.5}], ImageSize -> 500],
  Placed[LineLegend[{ColorData[97][1], Directive[ColorData[97][2], Dashed]},
    {"z(t)", Row[{"tanh(arctanh ", Subscript["z", "0"], " + ", Sqrt[Subscript["\[CapitalGamma]", "CI"]], " Q(t))"}]}, LegendLayout -> "Column"], Right]]
```

The noise term is zero for every rate, and the drift vanishes when $\Omega=\gamma_1=0$: then the record alone fixes $z$ on every run, and the two curves coincide. A drive or relaxation moves $z$ without writing it into the record, and the curves part.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A trapped particle watched continuously: squeezing below the zero-point

A particle of mass $m$ in a harmonic trap of frequency $\omega$ has its position measured continuously at rate $\lambda$ by a detector that catches a fraction $\eta$ of the signal. Can watching make it narrower in position than its ground state?

With $\hat H=\frac{\hat p^2}{2m}+\frac{m\omega^2\hat x^2}{2}$ and $\mathcal D[\hat x]\rho=\hat x\rho\hat x-\frac12\{\hat x^2,\rho\}$, the state conditioned on the record obeys

$$
d\rho=-\frac{i}{\hbar}[\hat H,\rho]\,dt+\lambda\,\mathcal D[\hat x]\rho\,dt+\sqrt{\eta\lambda}\left(\hat x\rho+\rho\hat x-2\langle\hat x\rangle\rho\right)dW .
$$

A Gaussian state stays Gaussian, and its covariances $\Sigma_x$, $C$, $\Sigma_p$ obey a Riccati equation with no noise in it:

$$
\dot\Sigma_x=\frac{2C}{m}-4\eta\lambda\Sigma_x^2,\qquad \dot C=\frac{\Sigma_p}{m}-m\omega^2\Sigma_x-4\eta\lambda\Sigma_xC,\qquad \dot\Sigma_p=-2m\omega^2C+\lambda\hbar^2-4\eta\lambda C^2 .
$$

The steady state is in closed form, keeping the root proved positive, and so is the squeezing condition; the curves integrate the flow from the ground state, dashed lines are the closed form, $\Sigma_x$ is in units of the zero-point variance $\hbar/(2m\omega)$, the strength is $\mu=2\hbar\lambda/(m\omega^2)$ and time is in units of $1/\omega$:

```wl
trapRiccati = {2 cxp/m - 4 η λ vx^2, vp/m - m ω^2 vx - 4 η λ vx cxp, -2 m ω^2 cxp + λ ℏ^2 - 4 η λ cxp^2};
trapDomain = m > 0 && λ > 0 && ω > 0 && ℏ > 0 && 0 < η <= 1;
trapSteady = SelectFirst[SolveValues[Thread[trapRiccati == 0], {vx, cxp, vp}],
    Resolve[ForAll[{m, λ, ω, ℏ, η}, trapDomain, (#[[1]] | #[[3]]) ∈ Reals && #[[1]] > 0 && #[[3]] > 0], Reals] &];
trapRatio = FullSimplify[trapSteady[[1]]/(ℏ/(2 m ω)) /. λ -> μ m ω^2/(2 ℏ), m > 0 && μ > 0 && ω > 0 && ℏ > 0 && 0 < η <= 1];
trapFlow[η0_, μ0_] := NDSolveValue[{Thread[{vxT'[t], cT'[t], vpT'[t]} == (trapRiccati /.
       {vx -> vxT[t], cxp -> cT[t], vp -> vpT[t], m -> 1, ω -> 1, ℏ -> 1, λ -> μ0/2, η -> η0})],
     vxT[0] == 1/2, cT[0] == 0, vpT[0] == 1/2}, vxT, {t, 0, 20}];
trapPanel[η0_] := Show[
    Plot[Evaluate[2 trapFlow[η0, #][t] & /@ {1/2, 5, 50}], {t, 0, 20}, PlotRange -> All, PlotLegends -> (Row[{"μ = ", #}] & /@ {1/2, 5, 50})],
    Plot[Evaluate[trapRatio /. {η -> η0, μ -> #} & /@ {1/2, 5, 50}], {t, 0, 20}, PlotStyle -> Dashed],
    PlotRange -> {0, 2.5}, Axes -> False, Frame -> True, GridLines -> {None, {1}}, FrameLabel -> {"ω t", "Σx / zero-point variance"},
    PlotLabel -> Row[{"η = ", η0, ",   squeezed for ", Reduce[trapRatio < 1 && μ > 0 /. η -> η0, μ, Reals]}], ImageSize -> 300];
trapPanels = trapPanel /@ {1, 1/2, 1/5};
Legended[PlotGrid[{First /@ trapPanels}, ImageSize -> 900], legendBelow[Last[First[trapPanels]]]]
```

A perfect detector squeezes the position below the zero-point at every strength. An imperfect one squeezes only above a threshold strength that grows as $\eta$ falls; below it, watching leaves the particle wider than its ground state.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A watched Gaussian: the center wanders, the widths do not

A particle of mass $m$ in a harmonic trap of frequency $\omega$ is watched in position at rate $\lambda$ by a detector of efficiency $\eta$, and its conditional state is a Gaussian with means $\langle\hat x\rangle,\langle\hat p\rangle$, variances $\Sigma_x,\Sigma_p$ and covariance $C$. Every run moves the packet's center at random. Do the widths of a watched Gaussian wander too?

With $\hat H=\frac{\hat p^2}{2m}+\frac{m\omega^2\hat x^2}{2}$ and $[\hat x,\hat p]=i\hbar$, the conditional mean of any observable $\hat O$ obeys

$$
d\langle\hat O\rangle=\left\langle \frac{i}{\hbar}[\hat H,\hat O]-\frac{\lambda}{2}\left[\hat x,[\hat x,\hat O]\right]\right\rangle dt+\sqrt{\eta\lambda}\,\left\langle\{\hat x-\langle\hat x\rangle,\hat O\}\right\rangle\,dW,
$$

and the second moments pick up the Itô correction, e.g. $d\Sigma_x=d\langle\hat x^2\rangle-2\langle\hat x\rangle d\langle\hat x\rangle-(d\langle\hat x\rangle)^2$. The Gaussian has the ordered generating function $\langle e^{\alpha\hat x}e^{\beta\hat p}\rangle=\exp\!\left(\alpha\langle\hat x\rangle+\beta\langle\hat p\rangle+\frac12(\alpha^2\Sigma_x+2\alpha\beta C+\beta^2\Sigma_p)+\frac{i\hbar}{2}\alpha\beta\right)$, and its purity is $\hbar/(2\sqrt{\Sigma_x\Sigma_p-C^2})$.

Ordering every product with $x$ to the left of $p$ and reading each ordered monomial off the generating function gives the drift and noise of the five numbers for general parameters, and the purity of every steady state of the drift:

```wl
xpAlgebra = NonCommutativeAlgebra[<|"Generators" -> {{xo, po}}, "CommutationRelations" -> {po ** xo - xo ** po + I ℏ}|>];
gaussMean[o_] := Block[{α, β}, With[{mgf = Exp[α mx + β mp + (α^2 Σx + 2 α β Cxp + β^2 Σp)/2 + I ℏ α β/2],
     q = NonCommutativeExpand[o, xpAlgebra] /. {GeneralizedPower[NonCommutativeMultiply, g_, k_] :> g^k, NonCommutativeMultiply -> Times}},
    Total[(#2 D[mgf, {α, #1[[1]]}, {β, #1[[2]]}] /. {α -> 0, β -> 0}) & @@@ CoefficientRules[q, {xo, po}]]]];
gaussDrift[o_] := gaussMean[I/ℏ Commutator[po ** po/(2 m) + m ω^2 xo ** xo/2, o] - λ/2 Commutator[xo, Commutator[xo, o]]];
gaussNoise[o_] := Sqrt[η λ] (gaussMean[Anticommutator[xo, o]] - 2 mx gaussMean[o]);
gaussCov[{a_, ma_}, {b_, mb_}] := With[{ab = (a ** b + b ** a)/2},
   {gaussDrift[ab] - ma gaussDrift[b] - mb gaussDrift[a] - gaussNoise[a] gaussNoise[b], gaussNoise[ab] - ma gaussNoise[b] - mb gaussNoise[a]}];
gaussEqs = Simplify[Join[{gaussDrift[#], gaussNoise[#]} & /@ {xo, po}, gaussCov @@@ {{{xo, mx}, {xo, mx}}, {{xo, mx}, {po, mp}}, {{po, mp}, {po, mp}}}]];
gaussShow = {mx -> AngleBracket["x"], mp -> AngleBracket["p"], Σx -> Subscript["Σ", "x"], Σp -> Subscript["Σ", "p"], Cxp -> "C"};
gaussSteadyPurity = Union[FullSimplify[ℏ/(2 Sqrt[Σx Σp - Cxp^2]) /. Thread[{Σx, Cxp, Σp} -> #], m > 0 && ω > 0 && λ > 0 && ℏ > 0 && 0 < η <= 1] & /@
    SolveValues[Thread[gaussEqs[[3 ;;, 1]] == 0], {Σx, Cxp, Σp}]];
Style[Column[{Grid[Prepend[MapThread[{TraditionalForm[#1], TraditionalForm[#2 /. gaussShow], TraditionalForm[#3 /. gaussShow]} &,
      {{"⟨x⟩", "⟨p⟩", Subscript["Σ", "x"], "C", Subscript["Σ", "p"]}, gaussEqs[[All, 1]], gaussEqs[[All, 2]]}], {"", "drift (dt)", "noise (dW)"}],
    Frame -> All, Alignment -> Left, Spacings -> {2, 0.8}], Row[{"every steady state of the drift column: purity = ", Row[Riffle[TraditionalForm /@ gaussSteadyPurity, ", "]]}]}, Center],
   14, FontFamily -> "Helvetica"]
```

The two columns define an Itô process, sampled four times per panel from $\langle\hat x\rangle=2$ and three times the zero-point variances, with $\hbar=m=1$, $\mu=2\hbar\lambda/(m\omega^2)=4$, time in units of $1/\omega$ and $\sqrt\eta$ dashed:

```wl
gaussRuns[η0_] := BlockRandom[SeedRandom[4]; RandomFunction[ItoProcess[{gaussEqs[[All, 1]], List /@ gaussEqs[[All, 2]]} /. Thread[{mx, mp, Σx, Cxp, Σp} -> Through[{mx, mp, Σx, Cxp, Σp}[t]]] /.
        {m -> 1, ω -> 1, ℏ -> 1, λ -> 2, η -> η0}, {{mx, mp, Σx, Cxp, Σp}, {2, 0, 3/2, 0, 3/2}}, t], {0., 6., 0.001}, 4]["Paths"]];
gaussPanel[η0_, runs_, {label_, f_, range_, grid_, legend_}] := ListLinePlot[Map[{#[[1]], f[#[[2]]]} &, runs, {2}], PlotRange -> range, PlotStyle -> If[label === "purity", ColorData[97][5], Automatic], Frame -> True, Axes -> False,
    GridLines -> {Automatic, If[grid === None, Automatic, Append[Table[{y, Directive[GrayLevel[0.8], Dashed]}, {y, 0.2, 1, 0.2}], {First[grid], Directive[Gray, Dashed, AbsoluteThickness[1.5]]}]]},
    GridLinesStyle -> Directive[GrayLevel[0.8], Dashed], FrameLabel -> {"ω t", label}, PlotLabel -> Row[{"η = ", η0, ",   μ = 4"}],
    Epilog -> If[legend, Inset[LineLegend[ColorData[97] /@ Range[4], Row[{"run ", #}] & /@ Range[4], LegendMarkerSize -> {12, 8}, LabelStyle -> 9, Spacings -> {0.3, 0.1}, LegendFunction -> (Framed[#, Background -> White, FrameStyle -> GrayLevel[0.8], FrameMargins -> 1] &)], Scaled[{0.02, 0.02}], {Left, Bottom}], {}], ImageSize -> 300];
With[{rows = Transpose[Function[η0, With[{runs = gaussRuns[η0]}, gaussPanel[η0, runs, #] & /@
       {{"⟨x⟩", First, All, None, η0 == 1/5}, {"purity", 1/(2 Sqrt[#[[3]] #[[5]] - #[[4]]^2]) &, {0, 1.05}, {Sqrt[η0]}, False}}]] /@ {1/5, 1}]},
  Column[{PlotGrid[{First[rows]}, ImageSize -> 650], PlotGrid[{Show[#, AspectRatio -> 1/2] & /@ Last[rows]}, ImageSize -> 650]}, Spacings -> 0]]
```

The noise of all three widths is zero for every parameter: the four runs put the packet in four places, so the means fan out, but they give it the same widths at the same time, and the purity curves fall on one line. The widths obey the Riccati equation in the drift column, where the efficiency scales the information terms $\eta\lambda$ but not the heating $\lambda\hbar^2$, and every steady state of that equation has purity $\sqrt\eta$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The quantum Kalman filter: five numbers against a wavefunction

A particle of mass $m$ sits in the potential $V(x)=\frac12m\omega^2x^2+\epsilon x^4$, and a perfect detector reads its position continuously at rate $\lambda$. Can five numbers, updated from the readings alone, follow the mean position that the full conditional wavefunction carries?

In a step $dt$ the detector returns a reading $\bar x$, distributed as $|\psi(x)|^2$ blurred by a Gaussian of variance $\frac{1}{4\lambda\,dt}$. Bayes' rule updates the wavefunction, and the Hamiltonian $\hat H=\hat p^2/2m+V(\hat x)$ moves it ($\hbar=1$):

$$
\psi\;\to\;e^{-i\hat H\,dt}\,\frac{e^{-\lambda dt\,(\hat x-\bar x)^2}\psi}{\lVert e^{-\lambda dt\,(\hat x-\bar x)^2}\psi\rVert}.
$$

The Kalman filter assumes a Gaussian state in the harmonic trap, $\epsilon=0$, and carries only its means and covariances, corrected by the innovation $\bar x-\langle\hat x\rangle$, the reading minus the predicted mean:

$$
d\langle\hat x\rangle=\frac{\langle\hat p\rangle}{m}dt+4\lambda\Sigma_x(\bar x-\langle\hat x\rangle)dt,\qquad d\langle\hat p\rangle=-m\omega^2\langle\hat x\rangle dt+4\lambda C(\bar x-\langle\hat x\rangle)dt,
$$
and its covariances $\Sigma_x$, $C$, $\Sigma_p$ obey the Riccati equations of the trapped-particle slide at $\eta=1$, with $\hbar=1$.

The wavefunction runs on 256 grid points with a split-operator step, from a coherent state displaced to $x=3$, and the filter receives the same readings. Lengths are in units of $\sqrt{\hbar/(m\omega)}$, time in $1/\omega$, the strength is $\mu=2\hbar\lambda/(m\omega^2)$ and $\epsilon$ is in units of $m^2\omega^3/\hbar$:

```wl
kalmanX = 24. Range[-128, 127]/256; kalmanP = 2 Pi RotateLeft[Range[-128, 127], 128]/24.; kalmanDt = 0.002;
kalmanGrid[μ_, ε_] := With[{λ0 = μ/2, half = Exp[-I (kalmanX^2/2 + ε kalmanX^4) kalmanDt/2], kin = Exp[-I kalmanP^2 kalmanDt/2]},
   Reap[NestList[With[{xb = RandomChoice[Abs[#]^2 -> kalmanX] + RandomVariate[NormalDistribution[0, 1/Sqrt[4 λ0 kalmanDt]]]},
       half InverseFourier[kin Fourier[half Normalize[Exp[-λ0 kalmanDt (kalmanX - Sow[xb])^2] #]]]] &, Normalize[Exp[-(kalmanX - 3)^2/2]], 5000]]];
kalmanFilter[μ_, readings_] := With[{λ0 = μ/2}, FoldList[Function[{v, xb}, v /. {mx_, mp_, vx_, c_, vp_} :> With[{g = 4 λ0 (xb - mx) kalmanDt},
       {mx, mp, vx, c, vp} + {mp kalmanDt + g vx, -mx kalmanDt + g c, (2 c - 4 λ0 vx^2) kalmanDt, (vp - vx - 4 λ0 vx c) kalmanDt, (-2 c + λ0 - 4 λ0 c^2) kalmanDt}]],
     {3., 0., 0.5, 0., 0.5}, readings][[All, 1]]];
kalmanPanel[μ_, ε_] := With[{run = BlockRandom[SeedRandom[7]; kalmanGrid[μ, ε]]},
    ListLinePlot[{Abs[#]^2 . kalmanX & /@ First[run], kalmanFilter[μ, run[[2, 1]]]}, DataRange -> {0, 10}, PlotStyle -> {Automatic, Dashed},
     PlotRange -> All, Frame -> True, FrameLabel -> {"ω t", "⟨x⟩"}, PlotLabel -> Row[{"μ = ", μ, ",   ε = ", ε}], ImageSize -> 300,
     PlotLegends -> {"wavefunction, 256 amplitudes", "Kalman filter, 5 numbers"}]];
kalmanPanels = Outer[kalmanPanel, {1, 10}, {0, 0.01, 0.05}];
Legended[PlotGrid[Map[First, kalmanPanels, {2}], ImageSize -> 800], legendBelow[Last[kalmanPanels[[1, 1]]]]]
```

In the harmonic trap the filter lies on the wavefunction's mean at both strengths: the conditional state stays Gaussian, so the 256 amplitudes carry nothing the five numbers miss. A quartic term makes the state non-Gaussian and the filter's harmonic model wrong. It drifts away under gentle watching. Under strong watching each reading pulls it back, and it falls visibly behind the faster true oscillation only at the larger $\epsilon$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### How much cavity Wigner negativity the atomic record certifies

An excited atom meets a cavity field, and after a time $t$ only the atom's inversion $z=P_e-P_g$ is recorded. The cavity, coupled to the atom at strength $g$, detuned by $\Delta$ and leaking photons at rate $\kappa$, can develop a Wigner function $W$ with negative regions. Fields with the same photon-number distribution $p_n$ give the atom the same record; how much negative volume $\int\max(0,-W)\,dx\,dp$ do all of them share?

With $\hbar=1$, $|e\rangle=|0\rangle$, $|g\rangle=|1\rangle$ and $\hat\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, the atom and cavity, frequencies $\omega_a$ and $\omega_c$, follow the Jaynes-Cummings Hamiltonian in the frame rotating at $\omega_c$, $\hat H=\frac{\Delta}{2}\hat\sigma_z+g(\hat\sigma_-\hat a^\dagger+\hat\sigma_+\hat a)$ with $\Delta=\omega_a-\omega_c$, and the loss enters through $\dot\rho=-i[\hat H,\rho]+\kappa\,\mathcal D[\hat a]\rho$. The coupling conserves the excitation number $\hat\sigma_+\hat\sigma_-+\hat a^\dagger \hat a$ and loss only lowers it, so an atom starting in $|e\rangle$ writes into $z(t)$ only the photon-number distribution $p_n$ of the field. Averaging the field over its phase keeps $p_n$ and averages $W$ over the angle, which can cancel negative regions but never create them; the averaged $W$ depends only on the radius, so its negative volume is a radial integral.

For the phase-randomized field, with $r^2=x^2+p^2$ and $L_n$ the Laguerre polynomial, the negative volume follows from $p_n(t)$ alone:

$$
\bar N(t)=2\pi\int_0^\infty r\,\max\left(0,\,-\sum_n p_n(t)\,W_n(r)\right)\,dr,\qquad W_n(r)=\frac{(-1)^n}{\pi}\,L_n(2r^2)\,e^{-r^2}.
$$

The atom starts in $|e\rangle$ with a coherent field of mean photon number 4, and with its phase-randomized version $\sum_np_n|n\rangle\langle n|$. The left panel plots the record $z(t)=\langle\hat\sigma_z\rangle$; the right panel plots the negative volume $N(t)=\int\max(0,-W)\,dx\,dp$ of the coherent field (solid) and $\bar N(t)$ (dashed); time is in units of $1/g$:

```wl
negA = AnnihilationOperator[20, {2}]; negField = CoherentState[20][2.];
negH[Δ_] := Δ/2 QuantumOperator["Z"] + QuantumOperator["-"] @ SuperDagger[negA] + SuperDagger[QuantumOperator["-"]] @ negA;
negVolume[ρ_] := Total[Ramp[-Flatten[WignerRepresentation[ρ, {-5., 5.}, {-5., 5.}, "GridSize" -> 81]["ValuesOnGrid"]]]]/64;
negRadii = Range[0.005, 6, 0.01];
negFock = Table[(-1)^m LaguerreL[m, 2 negRadii^2] Exp[-negRadii^2]/π, {m, 0, 19}];
negRadial[ρ_] := 2 π 0.01 negRadii . Ramp[-Re[Diagonal[Normal[ρ["DensityMatrix"]]]] . negFock];
negRun[{κ_, Δ_}, field_, volume_] := With[{evo = QuantumEvolve[negH[Δ], {negA} -> {κ}, QuantumTensorProduct[QuantumState["0"], field], {t, 0, 16}]},
    AssociationThread[{"record", "negativity"} -> Transpose @ Table[{{s, QuantumPartialTrace[evo[s], {2}]["BlochVector"][[3]]},
       {s, volume[QuantumPartialTrace[evo[s], {1}]]}}, {s, 0, 16, 0.1}]]];
negCases = {{0, 0}, {0.1, 0}, {0, 2}, {0.1, 2}};
negData = AssociationMap[Merge[{negRun[#, negField, negVolume], negRun[#, QuantumState[DiagonalMatrix[negField["ProbabilitiesList"]], 20], negRadial]}, Identity] &, negCases];
negAxes = <|"record" -> {Row[{"\[LeftAngleBracket]", Subscript["\[Sigma]", "z"], "\[RightAngleBracket]"}], {-1, 1}}, "negativity" -> {"negative volume", {0, 0.4}}|>;
negPanel[key_, {κ_, Δ_}] := ListLinePlot[negData[{κ, Δ}, key], PlotStyle -> {Automatic, Dashed}, PlotRange -> negAxes[key][[2]],
    PlotLegends -> {"coherent field", "phase-randomized field"}, Frame -> True, FrameLabel -> {"g t", negAxes[key][[1]]},
    ImageSize -> 260];
negPanels = AssociationMap[Function[c, negPanel[#, c] & /@ Keys[negAxes]], negCases];
Manipulate[Legended[Grid[{First /@ negPanels[{k, d}]}], legendBelow[Last[First[negPanels[{k, d}]]]]],
  Row[{Column[{Control[{{k, 0, "\[Kappa]/g"}, {0, 0.1}, SetterBar}], Control[{{d, 0, "\[CapitalDelta]/g"}, {0, 2}, SetterBar}]}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"κ/g = ", k}], <|0 -> "No loss: negativity builds after the collapse and returns.",
        0.1 -> "Loss: the guaranteed part (dashed) vanishes first, then all of it."|>[k]},
      {Row[{"Δ/g = ", d}], <|0 -> "The record swings about zero.",
        2 -> "The record sits above zero, and less negativity builds."|>[d]}}]]}], SaveDefinitions -> True]
```

The cavity's Wigner function shows what the two curves measure; red is positive and blue negative, on resonance and without loss, at four times in units of $1/g$:

```wl
negTimes = {0, 2, 6, 12};
negColor = Function[w, Blend[{RGBColor[0.2, 0.3, 0.8], White, RGBColor[0.85, 0.2, 0.15]}, 0.5 + 0.5 Sign[w] Sqrt[Min[Abs[w]/0.32, 1]]]];
negEvo[field_] := QuantumEvolve[negH[0], {negA} -> {0}, QuantumTensorProduct[QuantumState["0"], field], {t, 0, 16}];
negTile[evo_, s_] := With[{w = Re[WignerRepresentation[QuantumPartialTrace[evo[s], {1}], {-4., 4.}, {-4., 4.}, "GridSize" -> 81]["ValuesOnGrid"]]},
   ListDensityPlot[w, DataRange -> {{-4, 4}, {-4, 4}}, PlotRange -> All, ColorFunctionScaling -> False,
    ColorFunction -> negColor, Frame -> True, FrameTicks -> None, FrameStyle -> GrayLevel[0.7], ImageSize -> 110]];
Legended[Grid[Prepend[Function[f, Prepend[negTile[negEvo[Last[f]], #] & /@ negTimes, First[f]]] /@
    {"coherent field" -> negField, "phase-randomized field" -> QuantumState[DiagonalMatrix[negField["ProbabilitiesList"]], 20]},
   Prepend[Row[{"g t = ", #}] & /@ negTimes, ""]], Spacings -> {0.5, 0.5}], Placed[BarLegend[{negColor, {-0.32, 0.32}}, LegendLabel -> "W(x, p)"], Right]]
```

The two records coincide in every panel, so the atom cannot tell the coherent field from its phase-randomized version; in phase space the coherent field carries strong negative stripes, while its phase-randomized twin keeps only a thin negative ring. Any field with the same $p_n$ has at least the dashed negative volume, and the phase-randomized field reaches it, so the dashed curve is what they all share; it wiggles because the atom keeps trading photons with the field and reshapes $p_n$ on the Rabi time scale. The rest, carried by coherences between photon numbers, is most of the negativity between collapse and revival, and with loss the guaranteed part is the first to vanish.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Two qubits entangled by a measurement alone

Two qubits that never interact are read by one probe, a half-parity meter: its output $V$, averaged over a window $\Delta t$, is Gaussian about $-\delta v$ for $|00\rangle$, $+\delta v$ for $|11\rangle$ and $0$ for both $|01\rangle$ and $|10\rangle$, with white noise of spectral density $s$, and the meter also dephases the odd pair at rate $\gamma$. Starting from $|+\rangle|+\rangle$, can the record alone entangle the pair?

A window leaves a reading of variance $D=s/(2\Delta t)$. Coherences between states the meter tells apart are assumed to dephase fast, so the state keeps five numbers $x=(\rho_{00,00},\rho_{01,01},\rho_{10,10},\rho_{11,11},|\rho_{01,10}|)$, and Korotkov's quantum Bayes rule updates them for each reading:

$$
x_i\to\frac{x_iP_i(V)}{\sum_{j=1}^4x_jP_j(V)}\ (i\le4),\qquad x_5\to\frac{x_5P_2(V)\,e^{-\gamma\Delta t}}{\sum_{j=1}^4x_jP_j(V)},\qquad P_i(V)=e^{-(V-v_i)^2/2D},\quad v=(-\delta v,0,0,\delta v).
$$

On the states this rule produces, the Wootters concurrence is $C=2\max\left(0,\,x_5-\sqrt{x_1x_4}\right)$. Gaussian likelihoods multiply, so from $|+\rangle|+\rangle$ a record enters only through its running mean $\bar V$, and after a time $t$

$$
C(\bar V,t)=\max\!\left(0,\ \frac{e^{-\gamma t}-e^{-\delta v^2t/s}}{1+e^{-\delta v^2t/s}\cosh(2\delta v\,\bar V t/s)}\right).
$$

Five records per panel follow the update rule (lines), with $C(\bar V,t)$ marked at each record's own $\bar V$ (markers) and its ceiling $C(0,t)$ dashed; time is in units of $s/\delta v^2$:

```wl
hpmUpdate[x_, v_, dv_, var_, gdt_] := With[{p = Exp[-({-dv, 0, 0, dv} - v)^2/(2 var)]},
   Append[Take[x, 4] p, x[[5]] p[[2]] Exp[-gdt]]/(Take[x, 4] . p)];
hpmConc[vb_, τ_] := (Exp[-γ τ] - Exp[-δv^2 τ/s])/(1 + Exp[-δv^2 τ/s] Cosh[2 δv vb τ/s]);
hpmC[x_] := 2 Ramp[x[[5]] - Sqrt[x[[1]] x[[4]]]];
hpmDt = 0.01; hpmVar = 1/(2 hpmDt);
hpmRun[gs_, seed_] := BlockRandom[SeedRandom[seed]; NestList[With[{v = RandomChoice[Take[First[#], 4] -> {-1., 0., 0., 1.}] + RandomVariate[NormalDistribution[0, Sqrt[hpmVar]]]},
      {hpmUpdate[First[#], v, 1., hpmVar, gs hpmDt], Last[#] + v}] &, {ConstantArray[0.25, 5], 0.}, Round[5/hpmDt]]];
hpmPanel[gs_] := With[{runs = hpmRun[gs, #] & /@ Range[5], conc = Function[{vb, τ}, Ramp[hpmConc[vb, τ] /. {δv -> 1, s -> 1, γ -> gs}]]},
   Show[ListLinePlot[hpmC /@ #[[All, 1]] & /@ runs, DataRange -> {0, 5}, PlotRange -> All],
    ListPlot[Table[{k hpmDt, conc[#[[k + 1, 2]]/k, k hpmDt]}, {k, 25, Round[5/hpmDt], 25}] & /@ runs, PlotMarkers -> {Automatic, 9}, PlotRange -> All],
    Plot[conc[0, τ], {τ, 0, 5}, PlotStyle -> Directive[Gray, Dashed, Thick], PlotRange -> All],
    PlotRange -> {{0, 5}, {0, 1}}, Frame -> True, FrameLabel -> {Row[{Superscript["δv", 2], "t/s"}], "concurrence C"},
    PlotLabel -> Row[{"γ s/", Superscript["δv", 2], " = ", gs}], ImageSize -> 300]];
Labeled[Legended[PlotGrid[{hpmPanel /@ {0, 0.3, 1.5}}, ImageSize -> 900], legendBelow[Column[{"each color: one record",
    LineLegend[{Gray}, {"concurrence from the update rule"}], PointLegend[{Gray}, {Row[{"C(", OverBar["V"], ", t) at the record's running mean"}]}, LegendMarkers -> {"\[FilledCircle]", 9}],
    LineLegend[{Directive[Gray, Dashed, Thick]}, {"ceiling C(0, t)"}]}]]],
  Style[Row[{"ceiling:  ", TraditionalForm[hpmConc[0, t] /. Exp[e_] :> Superscript["e", Row[{"-", -e}]]], "  =  ", TraditionalForm[Normal[Series[hpmConc[0, t], {t, 0, 1}]]], " + O(", Superscript["t", 2], ")"}],
   13, FontFamily -> "Helvetica"], Top]
```

Each record's concurrence lands on the closed form at its own running mean: records that settle near $\bar V=0$ press against the dashed ceiling, while records pushed toward $\pm\delta v$ lose their concurrence. The ceiling starts with slope $(\delta v^2-\gamma s)/2s$: without dephasing it climbs toward $C=1$, with weak dephasing it peaks and falls, and once $\gamma s>\delta v^2$ its numerator stays negative and no record entangles the pair.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Continuous syndrome tracking: four numbers follow the records

Three qubits carry the bit-flip code, and each qubit flips at rate $\gamma$. The stabilizers $\hat S_1=Z_1Z_2$ and $\hat S_2=Z_2Z_3$ are measured weakly and continuously at rate $\kappa$ by detectors of efficiency $\eta$. Can four numbers, fed only these two noisy records, track the syndrome?

With $\mathcal D[\hat c]\rho=\hat c\rho\hat c^\dagger-\frac12\{\hat c^\dagger\hat c,\rho\}$, the three-qubit state conditioned on the records $dJ_k=2\sqrt{\eta\kappa}\,\langle\hat S_k\rangle dt+dW_k$ obeys

$$
d\rho=\gamma\sum_{i=1}^3\mathcal D[X_i]\rho\,dt+\kappa\sum_{k=1}^2\mathcal D[\hat S_k]\rho\,dt+\sqrt{\eta\kappa}\sum_{k=1}^2\left(\hat S_k\rho+\rho\hat S_k-2\langle\hat S_k\rangle\rho\right)dW_k .
$$

The four syndrome sectors, with projectors $\hat\Pi_{s}=\frac14(1+s_1\hat S_1)(1+s_2\hat S_2)$, have probabilities $p_s$. The filter relaxes them toward $1/4$ under the flips and multiplies each by the likelihood of the new pair of readings:

$$
p_s\;\to\;\frac{\left(\frac14+e^{-4\gamma dt}(p_s-\frac14)\right)\,e^{2\sqrt{\eta\kappa}\,(s_1dJ_1+s_2dJ_2)}}{\sum_{s'}(\cdots)} .
$$

The $8\times8$ equation is integrated from $(|000\rangle+i|100\rangle)/\sqrt2$, a superposition of two sectors, with the positivity-preserving step of Rouchon and Ralph, at two time steps on the same noise, and the filter runs on each record and is compared with the sector populations $\mathrm{Tr}[\hat\Pi_s\rho]$ of the full state through their largest gap; time is in units of $1/\kappa$:

```wl
qecStabs = Normal[QuantumOperator[#]["Matrix"]] & /@ {"ZZI", "IZZ"};
qecFlips = Normal[QuantumOperator[#]["Matrix"]] & /@ {"XII", "IXI", "IIX"};
qecSigns = {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
qecProj = (IdentityMatrix[8] + #[[1]] qecStabs[[1]]) . (IdentityMatrix[8] + #[[2]] qecStabs[[2]])/4 & /@ qecSigns;
qecStep[γ_, η_, dt_][ρ_, dJ_] := With[{σ = Sqrt[η] dJ . qecStabs},
   With[{M = (1 - (2 + 3 γ) dt/2 - η dt) IdentityMatrix[8] + σ + σ . σ/2},
    #/Tr[#] &[M . ρ . M + γ dt Total[# . ρ . # & /@ qecFlips] + (1 - η) dt Total[# . ρ . # & /@ qecStabs]]]];
qecRun[γ_, η_, dt_, kicks_] := Rest @ FoldList[Function[{state, kick}, With[{dJ = 2 Sqrt[η] Re[Tr /@ (qecStabs . First[state])] dt + kick},
      {qecStep[γ, η, dt][First[state], dJ], dJ}]], {Outer[Times, #, Conjugate[#]] &[(UnitVector[8, 1] + I UnitVector[8, 5])/Sqrt[2.]], {0., 0.}}, kicks];
qecFilter[γ_, η_, dt_, record_] := FoldList[#/Total[#] &[(1/4 + Exp[-4 γ dt] (#1 - 1/4)) Exp[2 Sqrt[η] qecSigns . #2]] &, {0.5, 0.5, 0., 0.}, record];
qecTrack[γ_, η_, dt_, kicks_] := With[{run = qecRun[γ, η, dt, kicks]}, With[{filter = qecFilter[γ, η, dt, run[[All, 2]]]},
    {filter, Max[Abs[Rest[filter] - (Re[Tr /@ (qecProj . First[#])] & /@ run)]]}]];
qecPanel[η_, γ_] := With[{kicks = RandomVariate[NormalDistribution[0, Sqrt[0.005]], {8000, 2}]},
   With[{fine = qecTrack[γ, η, 0.005, kicks], coarse = qecTrack[γ, η, 0.01, Total /@ Partition[kicks, 2]]},
    ListLinePlot[Transpose[First[fine]], DataRange -> {0, 40}, PlotRange -> {-0.02, 1.02}, Frame -> True, FrameLabel -> {"κ t", "syndrome probability"},
     PlotLabel -> Column[{Row[{"η = ", η, ",   γ/κ = ", γ}], Row[{"largest gap at κ dt = 0.01, 0.005: ", NumberForm[Last[coarse], 2], ", ", NumberForm[Last[fine], 2]}]}, Center],
     PlotLegends -> (Row[{"(", Subscript["S", "1"], ", ", Subscript["S", "2"], ") = (", #[[1]], ", ", #[[2]], ")"}] & /@ qecSigns), ImageSize -> 300]]];
qecPanels = BlockRandom[SeedRandom[5]; Outer[qecPanel, {1, 0.2}, {0.01, 0.05, 0.2}]];
Legended[PlotGrid[Map[First, qecPanels, {2}], ImageSize -> 900], legendBelow[Last[qecPanels[[1, 1]]]]]
```

The gap roughly halves with the time step, so it is the integration's error: the four probabilities follow the full state, coherences included, and need nothing else from it. With a perfect detector and rare flips the filter jumps cleanly from sector to sector as each flip happens. With a poor detector and frequent flips it rarely settles before the next flip.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The noise term in the syndrome equation of a watched code

Three qubits carry the bit-flip code, each qubit flips at rate $\gamma$, and the stabilizers $\hat S_1=Z_1Z_2$, $\hat S_2=Z_2Z_3$ are measured continuously at rate $\kappa$ by detectors of efficiency $\eta$. The four syndrome sectors $j$, with signs $s_{jk}=\pm1$ and projectors $\hat\Pi_j=\frac14(1+s_{j1}\hat S_1)(1+s_{j2}\hat S_2)$, have probabilities $p_j=\mathrm{Tr}[\hat\Pi_j\rho]$. Which term does the measurement record add to the equation for the $p_j$?

With $\mathcal D[\hat c]\rho=\hat c\rho \hat c^\dagger-\frac12\{\hat c^\dagger \hat c,\rho\}$, records $dJ_k=2\sqrt{\eta\kappa}\,\langle\hat S_k\rangle dt+dW_k$ and Wiener increments $dW_k$, the state conditioned on the records obeys

$$
d\rho=\gamma\sum_{i=1}^3\mathcal D[X_i]\rho\,dt+\kappa\sum_{k=1}^2\mathcal D[\hat S_k]\rho\,dt+\sqrt{\eta\kappa}\sum_{k=1}^2\left(\hat S_k\rho+\rho \hat S_k-2\langle\hat S_k\rangle\rho\right)dW_k ,
$$

and $dp_j=\mathrm{Tr}[\hat\Pi_j\,d\rho]$.

No matrix is ever built: $\rho$ is the most general element of the Pauli algebra, with 64 free coefficients, and requiring the identity for every coefficient turns the trace of $\hat\Pi_j$ against each term into a combination of the $p_j$. The flip term, the averaged stabilizer back-action and the record's noise term of the four-number equation are:

```wl
wsyGens = {wsyX1, wsyZ1, wsyX2, wsyZ2, wsyX3, wsyZ3};
wsyAlg = NonCommutativeAlgebra[<|"Generators" -> {wsyGens}, "CommutationRelations" -> (#[[2]] ** #[[1]] + #[[1]] ** #[[2]] & /@ Partition[wsyGens, 2]),
    "Relations" -> (# ** # - 1 & /@ wsyGens)|>];
wsyStab = {wsyZ1 ** wsyZ2, wsyZ2 ** wsyZ3};
wsySigns = {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
wsyProj[{a_, b_}] := (1 + a wsyStab[[1]]) ** (1 + b wsyStab[[2]])/4;
wsyTr[w_] := 8 ((NonCommutativeExpand[w, wsyAlg] /. _NonCommutativeMultiply -> 0) /. Thread[wsyGens -> 0]);
wsyRho = Total[MapIndexed[wsyA[First[#2]] If[Length[#1] < 2, Times @@ #1, NonCommutativeMultiply @@ #1] &, Subsets[wsyGens]]];
wsyP = Array[Subscript[p, #] &, 4];
wsyInPops[tr_] := Array[wsyC, 4] . wsyP /. First[SolveAlways[tr == Array[wsyC, 4] . (wsyTr[wsyProj[#] ** wsyRho] & /@ wsySigns), Array[wsyA, 64]]];
wsyDrift = Function[j, Simplify[γ wsyInPops[Total[wsyTr[wsyProj[j] ** (# ** wsyRho ** # - wsyRho)] & /@ {wsyX1, wsyX2, wsyX3}]]]] /@ wsySigns;
wsyFlipM = D[wsyDrift, {wsyP}] /. γ -> 1;
wsyBack = Function[j, κ wsyInPops[Total[wsyTr[wsyProj[j] ** (# ** wsyRho ** # - wsyRho)] & /@ wsyStab]]] /@ wsySigns;
wsyMean = wsyInPops[wsyTr[# ** wsyRho]] & /@ wsyStab;
wsyRecord = Outer[Sqrt[η κ] (wsyInPops[wsyTr[wsyProj[#1] ** (#2 ** wsyRho + wsyRho ** #2)]] - 2 wsyInPops[wsyTr[#2 ** wsyRho]] wsyInPops[wsyTr[wsyProj[#1] ** wsyRho]]) &, wsySigns, wsyStab, 1];
wsyCheck = Simplify[wsyRecord == 2 Sqrt[η κ] wsyP (Transpose[Transpose[wsySigns] - wsyMean])];
wsyBar = Subscript[OverBar["s"], #] & /@ {1, 2};
wsyNoise = MapIndexed[2 Sqrt[η κ] Subscript[p, First[#2]] HoldForm[Evaluate[Simplify[#1/(2 Sqrt[η κ] Subscript[p, First[#2]])]]] &, wsyRecord /. Thread[wsyMean -> wsyBar], {2}];
Style[Column[{Grid[Prepend[MapThread[{#1, Row[{"(", #2[[1]], ", ", #2[[2]], ")"}], TraditionalForm[#3], TraditionalForm[#4], Sequence @@ (TraditionalForm /@ #5)} &,
      {Range[4], wsySigns, wsyDrift, wsyBack, wsyNoise}],
     {"sector j", Row[{"(", Subscript["s", "j1"], ", ", Subscript["s", "j2"], ")"}], "flips (dt)", "stabilizer back-action (dt)",
      Row[{"record (d", Subscript["W", 1], ")"}], Row[{"record (d", Subscript["W", 2], ")"}]}],
    Frame -> All, FrameStyle -> GrayLevel[0.75], Alignment -> Center, Spacings -> {1.5, 0.9}, Background -> {None, {GrayLevel[0.93]}}, ItemStyle -> {None, {Bold}}],
   Grid[{{Row[{Subscript[OverBar["s"], "k"], " = Tr[", Subscript["S", "k"], "\[Rho]]"}], Row[Riffle[MapThread[Row[{#1, " = ", TraditionalForm[#2]}] &, {wsyBar, wsyMean}], ",   "]]},
     {Row[{"record term = 2", Sqrt["\[Eta] \[Kappa]"], " ", Subscript["p", "j"], "(", Subscript["s", "jk"], " \[Minus] ", Subscript[OverBar["s"], "k"], ")  for every j, k"}], wsyCheck}},
    Alignment -> Left, Spacings -> {2, 0.8}]}, Center, Spacings -> 1.2], 13, FontFamily -> "Helvetica"]
```

The four-number equation is integrated from the code sector, three runs and the mean of 1000 against the mean from the flip term alone; time is in units of $1/\kappa$:

```wl
wsyRun[e_, g_, m_] := BlockRandom[SeedRandom[3]; NestList[Nest[With[{dW = RandomVariate[NormalDistribution[0, Sqrt[0.001]], {m, 2}]},
        # + g # . Transpose[wsyFlipM] 0.001 + 2 Sqrt[e] # (dW . Transpose[wsySigns] - Total[(# . wsySigns) dW, {2}])] &, #, 50] &,
     ConstantArray[{1., 0., 0., 0.}, m], 400]];
wsyPanel[e_, g_] := With[{run = wsyRun[e, g, 1000][[All, All, 1]],
    flipsOnly = MatrixExp[g wsyFlipM #][[1, 1]] & /@ Subdivide[0, 20, 400]},
   ListLinePlot[Join[Transpose[run[[All, 1 ;; 3]]], {Mean /@ run, flipsOnly}], DataRange -> {0, 20},
    PlotStyle -> {Thin, Thin, Thin, Thick, Dashed}, PlotLegends -> {"run 1", "run 2", "run 3", "mean of 1000 runs", "mean from the flip term alone"},
    PlotRange -> {{0, 20}, {-0.02, 1.02}}, Frame -> True, FrameLabel -> {"κ t", Subscript["p", "code"]}, PlotLabel -> Row[{"η = ", e, ",   γ/κ = ", g}], ImageSize -> 300]];
wsyPanels = Outer[wsyPanel, {1, 0.2}, {0.02, 0.1}];
Legended[PlotGrid[Map[First, wsyPanels, {2}], ImageSize -> 700], legendBelow[Last[wsyPanels[[1, 1]]]]]
```

Each sector is pushed by its own probability times how far its signs sit from the running means, so the four probabilities form a closed equation, and the averaged back-action of the stabilizer measurement moves no probability. Conditioned on its record, each run jumps between near 0 and near 1, more sharply for the efficient detector, while the mean over runs follows the relaxation set by the flips alone, because the noise term has zero mean.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Markovian feedback that protects the code space, and what the detector costs

Three qubits carry the bit-flip code, $|0_L\rangle=|000\rangle$, $|1_L\rangle=|111\rangle$, and each qubit flips at rate $\gamma$ through the channel $\hat c_i=\sqrt\gamma X_i$. Each channel is monitored by homodyne detection of efficiency $\eta$, and its current $I_i$ drives the qubits at once through the Hamiltonian $I_i\hat F_i$, with $\hat F_i=i\sqrt\gamma\,{}[\hat\Pi_{\rm code},X_i]$ and $\hat\Pi_{\rm code}=|000\rangle\langle000|+|111\rangle\langle111|$. How well does this loop hold an encoded state?

With $\hbar=1$ and $\mathcal D[\hat c]\rho=\hat c\rho\hat c^\dagger-\frac12\{\hat c^\dagger\hat c,\rho\}$, the Wiseman-Milburn feedback master equation is

$$
\dot\rho=-\frac{i}{2}\sum_{i=1}^3\left[\hat c_i^\dagger\hat F_i+\hat F_i\hat c_i,\rho\right]+\sum_{i=1}^3\mathcal D[\hat c_i-i\hat F_i]\rho+\frac{1-\eta}{\eta}\sum_{i=1}^3\mathcal D[\hat F_i]\rho .
$$

The two identities $\hat c_i\hat F_i+\hat F_i\hat c_i=0$, which removes the Hamiltonian term, and $(\hat c_i-i\hat F_i)\hat\Pi_{\rm code}=0$ both hold exactly. The master equation evolves from $|0_L\rangle$ and from $|+_L\rangle=(|000\rangle+|111\rangle)/\sqrt2$, with $\mathcal D[\hat c_i]$ alone as the case without feedback; time is in units of $1/\gamma$:

```wl
fbFlips = Normal[QuantumOperator[#]["Matrix"]] & /@ {"XII", "IXI", "IIX"};
fbCode = DiagonalMatrix[UnitVector[8, 1] + UnitVector[8, 8]];
fbF = I (fbCode . # - # . fbCode) & /@ fbFlips;
fbChecks = {Max[Abs[MapThread[#1 . #2 + #2 . #1 &, {fbFlips, fbF}]]], Max[Abs[(# . fbCode) & /@ (fbFlips - I fbF)]]};
fbRho[ops_, rates_, ψ_] := Normal[QuantumEvolve[QuantumOperator[ConstantArray[0, {8, 8}]], QuantumOperator /@ ops -> rates, QuantumState[ψ], {t, 0, 3}]["DensityMatrix"]];
fbEtas = {1, 1/2, 1/10, 1/50};
fbStarts = {{Row[{"|", Subscript["0", "L"], "⟩"}], UnitVector[8, 1]}, {Row[{"|", Subscript["+", "L"], "⟩"}], (UnitVector[8, 1] + UnitVector[8, 8])/Sqrt[2]}};
fbCurves[ψ_] := Append[Function[η, fbRho[Join[fbFlips - I fbF, fbF], Join[{1, 1, 1}, ConstantArray[(1 - η)/η, 3]], ψ]] /@ fbEtas, fbRho[fbFlips, {1, 1, 1}, ψ]];
fbObs = {{Re[Tr[fbCode . #1]] &, "code-space population"}, {Re[Conjugate[#2] . #1 . #2] &, "fidelity with the start"}};
fbPanel[{name_, ψ_}, ρs_, {obs_, label_}] := Plot[Evaluate[obs[#, ψ] & /@ ρs], {t, 0, 3}, PlotRange -> {0, 1.02}, Frame -> True, FrameLabel -> {"γ t", label}, GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.8], Dashed],
   PlotLabel -> Row[{"start ", name}], PlotStyle -> {Automatic, Automatic, Automatic, Automatic, Dashed},
   PlotLegends -> Append[Row[{"η = ", #}] & /@ fbEtas, "no feedback"], ImageSize -> 300];
fbPanels = Function[start, With[{ρs = fbCurves[Last[start]]}, fbPanel[start, ρs, #] & /@ fbObs]] /@ fbStarts;
Labeled[Legended[PlotGrid[Transpose[Map[First, fbPanels, {2}]], ImageSize -> 700], legendBelow[Last[fbPanels[[1, 1]]]]],
  Style[Row[{Subscript["c", "i"], Subscript["F", "i"], " + ", Subscript["F", "i"], Subscript["c", "i"], " = ", fbChecks[[1]], "   and   (",
     Subscript["c", "i"], " - i ", Subscript["F", "i"], ")", Subscript["Π", "code"], " = ", fbChecks[[2]], "   for every i"}], 14, FontFamily -> "Helvetica"], Top]
```

At $\eta=1$ the two identities make every encoded state a fixed point, for $|0_L\rangle$ and $|+_L\rangle$ alike. Below unit efficiency the $\mathcal D[\hat F_i]$ term feeds detector noise back onto the qubits: the plateau falls with $\eta$, and at low efficiency the loop empties the code space faster than no feedback. From $|0_L\rangle$ the fidelity falls below the code-space population while from $|+_L\rangle$ the two agree: the noisy loop also flips $|0_L\rangle$ into $|1_L\rangle$, which no stabilizer sees.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Concluding remarks

- **A single record shows what the average hides.** The Lindblad equation is the mean over all detector records. Individual runs freeze and jump under strong watching, and collapse to a pole with Born's probabilities, while the average only decays smoothly.
- **Measurement informs and disturbs, and the efficiency sets the balance.** A perfect detector purifies the state; a lossy one leaves it mixed. The same trade-off fixes the bound that keeps every run inside the Bloch ball, the purity $\sqrt\eta$ of a watched Gaussian, the threshold for squeezing below the zero-point, and the cost of feedback.
- **The record can be enough.** When the conditioned state stays in a small family, a few numbers updated from the readings follow it exactly: the population from the integrated readout, five Gaussian moments in the Kalman filter, four syndrome probabilities for a watched code.
- **Watching changes the dynamics sharply.** The crossover from Rabi oscillation to Zeno freezing is an exceptional point. There the qubit is critically damped and relaxes fastest, and detuning splits that point and then removes it from the reachable measurement strengths.
- **The record is a resource.** It entangles qubits that never interact, and fed back as a drive it holds an encoded state in its code space.

Each of these results was reached the same way: the model written once, solved exactly where possible, sampled where not, and the two checked against each other.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Initialization

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/QEC/QECCore/QECCore.wl"]; Needs["Wolfram`QuantumFramework`QEC`"];
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/mbahram/TransmonLab/TransmonLab"]; Needs["TransmonLab`"];
legendBelow[\[FormalL]_] := Placed[\[FormalL] /. Placed[\[FormalX]_, ___] :> \[FormalX] /. Column[\[FormalY]_List, ___] :> Row[\[FormalY], Spacer[20]] /.
    (\[FormalH] : LineLegend | PointLegend | SwatchLegend)[\[FormalA]__] :> \[FormalH][Sequence @@ DeleteCases[{\[FormalA]}, LegendLayout -> _], LegendLayout -> {"Row", 1}], Below];
controlNotes[\[FormalR]_List] := Grid[MapAt[Style[#, Bold] &, \[FormalR], {All, 1}], Alignment -> {Left, Top}, ItemSize -> {{Automatic, 40}},
    Spacings -> {1, 0.6}, BaseStyle -> {12, LineIndent -> 0}];
```
