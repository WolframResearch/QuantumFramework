<!-- Build: MarkdownToNotebook["EU-Tour-2026-Core-Selection.md", "EU-Tour-2026-Core-Selection.nb", "EvaluateSeparator" -> None]. All cells share one evaluation context; the loading cell runs first. Example bodies are taken verbatim from Modeling-Slides.md, Catalogue sources and the TensorNetworks documentation. -->

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

# Quantum in the Wolfram Language

#| style: Subtitle

What QuantumFramework, TensorNetworks and the Wolfram Language compute, one problem at a time

#| style: Subsubtitle

Mads Bahrami, Wolfram Research. Europe tour, November 2026. quantum@wolfram.com

The loading cell uses the working copies of QuantumFramework and TensorNetworks; the examples that need QECCore or TransmonLab load them where they are used.

Every example below is a physics question, the model, one runnable cell and what the result says. The chapters follow capabilities.

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"];
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
Needs["Wolfram`TensorNetworks`"]; Needs["Wolfram`TensorNetworks`Symmetry`"];
legendBelow = Function[\[FormalL], Placed[\[FormalL] /. Placed[\[FormalX]_, ___] :> \[FormalX] /. Column[\[FormalY]_List, ___] :> Row[\[FormalY], Spacer[20]] /.
    (\[FormalH] : LineLegend | PointLegend | SwatchLegend)[\[FormalA]__] :> \[FormalH][Sequence @@ DeleteCases[{\[FormalA]}, LegendLayout -> _], LegendLayout -> {"Row", 1}], Below]];
controlNotes = Function[\[FormalR], Grid[MapAt[Style[#, Bold] &, \[FormalR], {All, 1}], Alignment -> {Left, Top}, ItemSize -> {{Automatic, 40}},
    Spacings -> {1, 0.6}, BaseStyle -> {12, LineIndent -> 0}]];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Open systems and dynamics

The examples in this section are systems that lose energy and phase to their surroundings. Each is a Lindblad equation that QF builds from a Hamiltonian and jump operators with symbolic rates and evolves exactly.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Correlated decay: bright, dim and dark states

Two identical atoms, frequency $\omega$, sit close enough to radiate into the same field modes, so each decays at rate $\gamma$ and the field correlates their emission with strength $0\le r\le1$. Which states radiate faster or slower than a lone atom, and can one stop radiating altogether?

With $\hbar=1$, $|e\rangle=|0\rangle$ and $|g\rangle=|1\rangle$, the Lindblad equation in Kossakowski form, with rate matrix $\gamma\begin{pmatrix}1&r\\ r&1\end{pmatrix}$, is

$$
\dot\rho=-i\left[\frac{\omega}{2}\left(\hat\sigma_z^1+\hat\sigma_z^2\right),\rho\right]+\gamma\sum_{i,j=1}^{2}\begin{pmatrix}1&r\\ r&1\end{pmatrix}_{ij}\left(\hat\sigma_-^i\rho\,\hat\sigma_+^j-\frac12\left\{\hat\sigma_+^j\hat\sigma_-^i,\rho\right\}\right),
$$

and in the frame rotating at $\omega$ the Hamiltonian drops out; the field's coherent exchange between the atoms only shifts the energies of $|T\rangle$ and $|S\rangle$ below, without mixing them, so it leaves $\langle\hat N\rangle$ unchanged and is left out.

The rate matrix enters as it stands, and the solution for general $r$ gives the mean number of excited atoms, $\langle\hat N\rangle=\mathrm{Tr}\left[\rho\,(\hat n_1+\hat n_2)\right]=2\rho_{ee,ee}+\rho_{eg,eg}+\rho_{ge,ge}$ with $\hat n_i=|e\rangle\langle e|_i$, from four initial states, with $|T\rangle,|S\rangle=(|eg\rangle\pm|ge\rangle)/\sqrt2$. Each curve is $\langle\hat N\rangle$ divided by its starting value $N_0$, 2 for $|ee\rangle$ and 1 for the others, and the shading marks where a state decays faster (orange) or slower (blue) than a lone atom, whose $\langle\hat N\rangle=e^{-\gamma t}$; time is in units of $1/\gamma$:

```wl
statesA3 = <|"|ee\[RightAngleBracket]" -> QuantumState["00"], "|T\[RightAngleBracket]" -> QuantumState[{0, 1, 1, 0}/Sqrt[2]],
    "|S\[RightAngleBracket]" -> QuantumState[{0, 1, -1, 0}/Sqrt[2]], "|eg\[RightAngleBracket]" -> QuantumState["01"]|>;
excitationA3 = Simplify[ComplexExpand @ Tr[Normal[QuantumEvolve[QuantumOperator[0 IdentityMatrix[4]], {QuantumOperator["-", {1}], QuantumOperator["-", {2}]} -> {{1, r}, {r, 1}},
       #, t]["DensityMatrix"]] . DiagonalMatrix[{2, 1, 1, 0}]], 0 <= r <= 1 && t >= 0] & /@ statesA3;
perExcitationA3 = Values[excitationA3]/{2, 1, 1, 1};
{superColorA3, subColorA3} = {RGBColor[1, 0.55, 0.3], RGBColor[0.55, 0.75, 1]};
panelA3[rv_] := Show[
    Plot[Exp[-t], {t, 0, 5}, PlotStyle -> None, Filling -> {1 -> {Bottom, Directive[Opacity[0.22], superColorA3]}}],
    Plot[Exp[-t], {t, 0, 5}, PlotStyle -> None, Filling -> {1 -> {1.02, Directive[Opacity[0.22], subColorA3]}}],
    Plot[Evaluate[Limit[#, r -> rv] & /@ perExcitationA3], {t, 0, 5}, PlotStyle -> ColorData[97] /@ {3, 4, 5, 1}],
    PlotRange -> {{0, 5}, {0, 1.02}}, Frame -> True, Axes -> False, GridLines -> {Automatic, {0, 0.25, 0.5, 0.75, 1}}, GridLinesStyle -> Dashed,
    FrameLabel -> {"γ t", Row[{"\[LeftAngleBracket]N\[RightAngleBracket] / ", Subscript["N", 0]}]},
    PlotLabel -> Row[{"correlation r = ", Style[TraditionalForm[rv], FractionBoxOptions -> {Beveled -> True}]}], ImageSize -> 260];
Legended[PlotGrid[{panelA3 /@ {0, 1/2, 1}}, ImageSize -> 760],
  Placed[Row[{LineLegend[ColorData[97] /@ {3, 4, 5, 1}, Keys[statesA3], LegendLayout -> "Row"], Spacer[30],
     SwatchLegend[{Directive[Opacity[0.22], superColorA3], Directive[Opacity[0.22], subColorA3]}, {"faster than a lone atom: superradiant", "slower than a lone atom: subradiant"}, LegendLayout -> "Row"]}], Below]]
```

Without correlation, at $r=0$, the atoms decay independently and all four curves, $|ee\rangle$ per excitation included, lie on that of a lone atom, so they show as one curve. Correlation makes the symmetric state superradiant, decaying faster than a lone atom, and the antisymmetric one subradiant, decaying slower: $\langle\hat N\rangle_T=e^{-(1+r)\gamma t}$ and $\langle\hat N\rangle_S=e^{-(1-r)\gamma t}$. At $r=1$ the singlet is dark, so $|eg\rangle$, half singlet, keeps half an excitation, while $|ee\rangle$ empties completely and, per excitation, faster than a lone atom.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### When a noisy circuit matches the master equation

A Ramsey experiment applies a Hadamard to $|g\rangle$, waits a time $t$, applies a second Hadamard and reads $P(e)$. During the wait the qubit, detuned by $\delta$ from the drive, decays at rate $\gamma_1$ and dephases at rate $\gamma_\phi$, and a residual drive of Rabi frequency $\Omega$ may act. Device simulators model such a wait not with a master equation but as a circuit: a gate followed by textbook noise channels. When do the two give the same fringe?

With $\hbar=1$, $|g\rangle=|0\rangle$ (the state the amplitude-damping channel relaxes to), $|e\rangle=|1\rangle$, $\hat\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, qubit frequency $\omega$ and drive frequency $\omega_d=\omega-\delta$, the wait obeys, in the frame rotating at $\omega_d$,

$$
\dot\rho=-i\left[\frac{\delta}{2}\hat\sigma_z+\frac{\Omega}{2}\hat\sigma_x,\rho\right]+\gamma_1\,\mathcal D[\hat\sigma_-]\rho+\gamma_\phi\,\mathcal D[\hat\sigma_z]\rho .
$$

The circuit cuts the wait into $m$ equal slices. Each slice is the unitary $\mathcal U_{t/m}\rho=e^{-i\hat Ht/m}\rho\,e^{i\hat Ht/m}$, then amplitude damping $\mathcal A_\lambda$ with decay probability $\lambda=1-e^{-\gamma_1t/m}$, then phase damping $\mathcal P_\mu$ with $\mu=1-e^{-4\gamma_\phi t/m}$, which shrinks the coherence by $e^{-2\gamma_\phi t/m}$, the same factor as the dephasing term. With $\mathcal L$ the right-hand side of the master equation, the circuit reproduces the wait in the limit

$$
e^{t\mathcal L}\rho=\lim_{m\to\infty}\left(\mathcal P_\mu\circ\mathcal A_\lambda\circ\mathcal U_{t/m}\right)^m\rho ,
$$

and with $m=1$ the whole wait is one gate and two channels, the usual circuit noise model.

Lines are the master equation, dots the circuit with $m=1$ and, when the drive is on, squares the circuit with $m=8$; the dephasing rate and the drive are set separately, and time is in units of $1/\gamma_1$:

```wl
hamA4[δ_, Ω_] := -δ/2 QuantumOperator["Z"] + Ω/2 QuantumOperator["X"]; jumpsA4 = {QuantumOperator[{{0, 1}, {0, 0}}], QuantumOperator["Z"]};
sliceA4[h_, {γ1v_, γφv_}, τ_, {lu_, lad_, lpd_}] := {QuantumOperator[Exp[-I τ h], "Label" -> lu],
    QuantumChannel["AmplitudeDamping"[1 - Exp[-γ1v τ]], "Label" -> lad], QuantumChannel["PhaseDamping"[1 - Exp[-4 γφv τ]], "Label" -> lpd]};
pEA4[ρ_] := ρ["ProbabilitiesList"][[2]];
viaMEA4[h_, γs_, T_] := pEA4[QuantumOperator["H"][QuantumEvolve[h, jumpsA4 -> γs, QuantumOperator["H"][QuantumState["0"]], T]]];
circuitA4[h_, γs_, m_, s_, labels_ : {"U", "AD", "PD"}, group_ : Catenate] := QuantumCircuitOperator[Join[{"H"}, group[ConstantArray[sliceA4[h, γs, s/m, labels], m]], {"H"}]];
labelsA4[τ_] := {Subscript["\[ScriptCapitalU]", τ], Subscript["\[ScriptCapitalA]", "λ"], Subscript["\[ScriptCapitalP]", "μ"]};
diagramsA4 = {circuitA4[hamA4[δ, Ω], {γ1, γφ}, 1, t, labelsA4["t"], Map[QuantumCircuitOperator[#, "slice"] &]]["Diagram", ImageSize -> 330],
    circuitA4[hamA4[δ, Ω], {γ1, γφ}, 8, t, labelsA4["t/8"], Map[QuantumCircuitOperator[#, "slice"] &]]["Diagram", ImageSize -> 760]};
transferA4[h_, γs_] := With[{ρ = Array[\[FormalR], {2, 2}]}, Function[\[FormalS], Evaluate[D[Flatten[Normal[
      QuantumCircuitOperator[sliceA4[h, γs, \[FormalT], {"U", "AD", "PD"}]][QuantumState[ρ]]["DensityMatrix"]]], {Flatten[ρ]}] /. \[FormalT] -> \[FormalS]]]];
pointsA4[T_, m_] := With[{hm = Normal[QuantumOperator["H"]["Matrix"]]},
   Table[{x, Re[(hm . Partition[MatrixPower[T[x/m], m] . Flatten[hm . {{1, 0}, {0, 0}} . hm], 2] . hm)[[2, 2]]]}, {x, 0, 4, If[m == 1, 0.2, 0.4]}]];
panelA4[{γφv_, Ωv_}, δv_] := With[{me = viaMEA4[hamA4[δv, Ωv], {1, γφv}, {t, 0, 4}], T = transferA4[N[hamA4[δv, Ωv]], {1, γφv}]},
   ListPlot[Prepend[pointsA4[T, #] & /@ If[Ωv == 0, {1}, {1, 8}], Table[{x, me /. t -> x}, {x, 0, 4, 0.04}]],
    Joined -> {True, False, False}, PlotMarkers -> {None, {"\[FilledCircle]", 7}, {"\[EmptySquare]", 11}},
    PlotLegends -> Take[{"master equation", "circuit, m = 1", "circuit, m = 8"}, If[Ωv == 0, 2, 3]], PlotRange -> {0, 1}, Frame -> True, GridLines -> {None, {0.25, 0.5, 0.75}}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.35], Thin], FrameLabel -> {Row[{Subscript["γ", 1], " t"}], "P(e)"},
    Epilog -> Inset[Row[{"δ/", Subscript["γ", 1], " = ", δv}], Scaled[{0.96, 0.96}], {Right, Top}], ImageSize -> Small]];
panelsA4 = AssociationMap[Function[c, panelA4[c, #] & /@ {0, 3, 8}], Tuples[{{0, 1/2, 1}, {0, 1, 2, 4}}]];
Column[{Column[diagramsA4, Alignment -> Center, Spacings -> 1],
  With[{legendBelow = legendBelow, controlNotes = controlNotes, panelsA4 = panelsA4},
 Manipulate[Legended[PlotGrid[{First /@ panelsA4[{gφ, gΩ}]}, ImageSize -> 660], legendBelow[Last[panelsA4[{gφ, gΩ}][[1]]]]],
   Grid[{{Grid[{{Style["dephasing", 12], Control[{{gφ, 1/2, Row[{Subscript["γ", "φ"], "/", Subscript["γ", 1]}]}, {0 -> "0", 1/2 -> "0.5", 1 -> "1"}, SetterBar}]},
      {Style["drive", 12], Control[{{gΩ, 2, Row[{"Ω/", Subscript["γ", 1]}]}, {0 -> "0", 1 -> "1", 2 -> "2", 4 -> "4"}, SetterBar}]}}, Alignment -> {{Right, Left}, Center}, Spacings -> {1, 0.8}],
    Dynamic[Pane[Style[<|{0, 0} -> "No dephasing, no drive: the routes agree; the fringe decays at γ₁/2.",
        {1/2, 0} -> Row[{"Dephasing, no drive: the routes still agree; decay at 1/T₂ = γ₁/2 + 2", Subscript["γ", "φ"], "."}],
        {1, 0} -> Row[{"Dephasing, no drive: the routes still agree; decay at 1/T₂ = γ₁/2 + 2", Subscript["γ", "φ"], "."}],
        {0, 1} -> "Drive, no dephasing: the one‐step circuit misses; with 8 short steps it stays close.",
        {0, 2} -> "Stronger drive: the one‐step circuit misses further; 8 short steps still come close.",
        {0, 4} -> "Strong drive: the one‐step circuit is far off; 8 steps come close but not all the way.",
        {1/2, 1} -> "Drive and dephasing: a smaller gap for one step; 8 steps close most of it.",
        {1/2, 2} -> "The circuit matches only when the time between gates is very short.",
        {1/2, 4} -> "Strong drive and dephasing: even 8 steps leave a visible gap.",
        {1, 1} -> "Strong dephasing: the fringe swings less about 1/2; the one‐step gap shrinks, the 8‐step gap grows.",
        {1, 2} -> "Strong dephasing and drive: 8 steps help less than without dephasing.",
        {1, 4} -> "Strong drive and dephasing: the largest gap at m = 8; the steps must be shorter still."|>[{gφ, gΩ}], 12, LineIndent -> 0], ImageSize -> 520]]}}, Spacings -> 2, Alignment -> {Left, Center}]]]}, Alignment -> Center, Spacings -> 2]
```

Without a drive the two routes agree for every value of every parameter, because the rotation, the decay and the dephasing generators commute: the wait is exactly one rotation followed by two named channels. The fringe oscillates about $P(e)=1/2$ at the detuning, and $\gamma_1$ enters its decay only through $1/T_2$. A drive during the wait stops the generators commuting: off resonance the circuit then approaches the master equation only as the time between gates gets shorter, while at $\delta=0$ they still give the same fringe. In these panels dephasing narrows the swing of the fringe about $1/2$ and shrinks the gap of the one-step circuit, but turning it on widens the gap left by eight steps, so with dephasing on the steps must be shorter.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Fastest relaxation at the Zeno exceptional point

A qubit is driven with Rabi frequency $\Omega$, detuned by $\Delta$ from the drive, while a detector measures $\hat\sigma_z$ at strength $k$. Weak watching damps the Rabi oscillation and strong watching freezes it; what happens to the dynamics at the crossover?

With $\hbar=1$, $|e\rangle=|0\rangle$, $|g\rangle=|1\rangle$, $\hat\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, qubit frequency $\omega$ and drive frequency $\omega-\Delta$, the Hamiltonian in the frame rotating with the drive is $\hat H=\frac{\Delta}{2}\hat\sigma_z+\frac{\Omega}{2}\hat\sigma_x$, and averaging over the detector record gives

$$
\dot\rho=\mathcal L\rho=-i[\hat H,\rho]+k\,\mathcal D[\hat\sigma_z]\rho .
$$

Away from exceptional points, every solution is a sum of $e^{\lambda t}$ over the eigenvalues $\lambda$ of the superoperator $\mathcal L$. At an exceptional point two eigenvalues and their eigenvectors coincide, so $\mathcal L$ cannot be diagonalized; there the discriminant of its characteristic polynomial, with the steady state's root $\lambda=0$ divided out, vanishes.

On resonance the Rabi pair solves $\lambda^2+2k\lambda+\Omega^2=0$, so $\lambda=-k\pm\sqrt{k^2-\Omega^2}$: a damped oscillation for $k<\Omega$, and for $k>\Omega$ two real rates, one of which, about $\Omega^2/2k$, slows as the watching grows stronger, the Zeno effect. The two meet at $k=\Omega$. That discriminant has a zero at some $k>0$, an exceptional point, only for:

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
Column[{Labeled[PlotGrid[{discPanelA5 /@ ΔsA5}, ImageSize -> 800], Style["log|D(k)| against complex measurement strength k, with D(k) the discriminant of the Liouvillian's characteristic polynomial", 13, FontFamily -> "Helvetica"], Top],
  Labeled[Legended[PlotGrid[{ratePanelA5 /@ ΔsA5}, ImageSize -> 800], Placed[LineLegend[ColorData[97] /@ {1, 2}, {"relaxation rate (slowest decay)", "oscillation frequency"}, LegendLayout -> {"Row", 1}], Below]],
   Style["the qubit's relaxation rate and oscillation frequency for real k", 13, FontFamily -> "Helvetica"], Top]}, Spacings -> 2]
```

On resonance the slowest decay rate rises to a sharp peak at $k=\Omega$, directly below the exceptional point: the qubit relaxes fastest there, critically damped between Rabi oscillation and the Zeno slowdown, whose rate falls as $\Omega^2/2k$. A small detuning splits the exceptional point into two zeros on the real axis, one where the Rabi pair turns real and one where two real rates meet and turn complex again; at $|\Delta/\Omega|=1/(2\sqrt2)$ the two zeros merge, and for $|\Delta/\Omega|>1/(2\sqrt2)$ they leave the axis as a complex-conjugate pair, so no measurement strength reaches them, the Rabi pair never turns real, and the sharp peak of the decay rate becomes a smooth hump.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Entanglement sudden death

Two qubits A and B start in $\cos\theta\,|gg\rangle+\sin\theta\,|ee\rangle$, and each decays on its own, A at rate $\gamma_A$ and B at rate $\gamma_B$. Each qubit's excitation decays exponentially; does the entanglement also fade only gradually?

With $\hbar=1$, $|e\rangle=|0\rangle$, $|g\rangle=|1\rangle$ and qubit frequencies $\omega_A$, $\omega_B$, the Lindblad equation is

$$
\dot\rho=-i\left[\frac{\omega_A}{2}\hat\sigma_z^A+\frac{\omega_B}{2}\hat\sigma_z^B,\rho\right]+\gamma_A\,\mathcal D[\hat\sigma_-^A]\rho+\gamma_B\,\mathcal D[\hat\sigma_-^B]\rho,
$$

and in the frame rotating with each qubit the Hamiltonian drops out while the two dissipators stay as they are.

The concurrence of the evolving state, computed for $0\le\theta\le\pi/2$ and general rates and $t$, is

$$
C(t)=2\,e^{-(\gamma_A+\gamma_B)t/2}\sin\theta\,\max\left(0,\ \cos\theta-\sin\theta\sqrt{(1-e^{-\gamma_At})(1-e^{-\gamma_Bt})}\right),
$$

and the curves plot it, with a dot where a curve reaches zero; time is in units of $1/\gamma_A$:

```wl
concA2 = QuantumEntanglementMonotone[QuantumEvolve[QuantumOperator[0 IdentityMatrix[4]], {QuantumOperator["-", {1}], QuantumOperator["-", {2}]} -> {γA, γB},
     QuantumState[{Sin[θ], 0, 0, Cos[θ]}], t], "Concurrence"];
{θsA2, γBsA2} = {{π/8, π/4, π/3, 5 π/12}, {0, 3/10, 1}};
deathA2[γBv_, θv_] := If[θv > π/4 && γBv > 0, First[NSolveValues[(1 - Exp[-t]) (1 - Exp[-γBv t]) == Cot[θv]^2 && t > 0, t, Reals]], None];
panelA2[γBv_] := Plot[Evaluate[concA2 /. {γA -> 1, γB -> γBv, θ -> #} & /@ θsA2], {t, 0, 4}, PlotRange -> {-0.05, 1.02}, WorkingPrecision -> 20, Exclusions -> None,
    GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.5, 0.35], Thin], Epilog -> {PointSize[0.035], MapIndexed[If[NumericQ[#1], {ColorData[97][First[#2]], Point[{#1, 0}]}, Nothing] &, deathA2[γBv, #] & /@ θsA2]},
    PlotStyle -> {Automatic, Dashed, Automatic, Automatic}, PlotLegends -> (Row[{"θ = ", #}] & /@ θsA2),
    Frame -> True, FrameLabel -> {Row[{Subscript["γ", "A"], " t"}], "concurrence"}, PlotLabel -> Row[{Subscript["γ", "B"], "/", Subscript["γ", "A"], " = ", γBv}], ImageSize -> Small];
panelsA2 = panelA2 /@ γBsA2;
Legended[PlotGrid[{First /@ panelsA2}, ImageSize -> 660], legendBelow[Column[{Last[First[panelsA2]], PointLegend[{Gray}, {"entanglement dies"}, LegendMarkers -> {"\[FilledCircle]", 9}]}]]]
```

It vanishes once $(1-e^{-\gamma_At})(1-e^{-\gamma_Bt})$ reaches $\cot^2\theta$. That needs both qubits to decay and $\pi/4<\theta<\pi/2$, where $|ee\rangle$ outweighs $|gg\rangle$: those states lose all their entanglement at a finite time and never regain it, and the more slowly B decays, the later that time comes. At $\theta\le\pi/4$ (the boundary case $\theta=\pi/4$ is dashed), or if B does not decay, the entanglement only fades.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Measurement and estimation

The examples in this section follow a system while it is measured: single runs under continuous measurement, the limit a detector's efficiency sets, what randomized benchmarking can and cannot see, and how much entanglement helps a sensor that dephases. The two trajectory slides use stochastic master equations written as WL Itô processes; QF supplies the operators, the channels and the quantum Fisher information throughout.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Watched qubits pin and jump

A qubit is driven with Rabi frequency $\Omega$, detuned by $\Delta$ from its transition, while a detector continuously measures $\hat\sigma_z$ at strength $k$ and catches a fraction $\eta$ of the signal; averaged over the detector record, this is the qubit of the exceptional-point example, now followed run by run. Does watching stop the oscillation?

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
measuredQubitLesson = <|0.1 -> "Runs oscillate; at η = 1 the average decays only because they fall out of phase.",
   0.5 -> "Runs still oscillate, kicked irregularly by the record.",
   1 -> "Crossover (on resonance): oscillation stops and runs linger near the poles.",
   5 -> "Zeno: runs stick at a pole and jump; the average only decays.",
   10 -> "Stronger watching, rarer jumps."|>;
With[{measuredQubitPanels = measuredQubitPanels, measuredQubitLesson = measuredQubitLesson},
 Manipulate[Legended[PlotGrid[{measuredQubitPanels[{k, d}]}, ImageSize -> 800],
   Placed[LineLegend[Append[ColorData[97] /@ Range[3], Directive[Red, Dashed]], {"run 1", "run 2", "run 3", "Lindblad average"}, LegendLayout -> {"Row", 1}], Below]],
  Row[{Column[{Control[{{k, 0.1, "measurement strength k/\[CapitalOmega]"}, {0.1, 0.5, 1, 5, 10}, SetterBar}], Control[{{d, 0, "detuning \[CapitalDelta]/\[CapitalOmega]"}, {0, 0.5, 1, 2}, SetterBar}]}],
    Spacer[20], Dynamic[Grid[{
      {Style[Row[{"k/\[CapitalOmega] = ", k}], Bold, 12], measuredQubitLesson[k]},
      {Style["\[Eta] < 1", Bold, 12], "Runs stay mixed and miss the poles."},
      {Style["\[CapitalDelta] \[NotEqual] 0", Bold, 12], "Weak watching only: smaller swing, slow drift between poles."}},
     Alignment -> {Left, Top}, ItemSize -> {{6, 44}}, Spacings -> {1, 0.6}, BaseStyle -> {12, LineIndent -> 0}]]}]]]
```

At $k/\Omega=0.1$ the runs keep Rabi-oscillating. At $k/\Omega=10$ they stay pinned near $z=1$, with brief spikes that fall back and now and then a full jump to the other pole: that is the quantum Zeno effect, which the average shows only as a slow decay. A detector that catches less of the signal leaves each run mixed, so under weak measurement its oscillation no longer reaches the poles. Detuning shrinks the Rabi swing to a fast wiggle near the starting pole, and under weak watching the record then carries the runs slowly from pole to pole; once the measurement is much stronger than the detuning, the runs look as they do on resonance.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A monitored qubit stays in the Bloch ball: Itô's lemma on the sphere

A qubit is Rabi-driven at frequency $\Omega$, relaxes toward $z=1$ at rate $\gamma_1$ and dephases at rate $\gamma_\phi$. It is read out along $\hat\sigma_z$ by a detector that dephases it at rate $\Gamma_d$; the detector's noise enters $z$ at the informational rate $\Gamma_{CI}$ and rotates the phase at the back-action rate $\Gamma_{BA}$, and its efficiency is $\eta=(\Gamma_{CI}+\Gamma_{BA})/(2\Gamma_d)$. Does every run stay a density matrix, $r=|\vec r|\le1$, for every choice of rates?

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

### Randomized benchmarking: same decay, different floors

A qubit starts in $|0\rangle$, undergoes $m$ Clifford gates drawn at random and then the Clifford that undoes their product, and is measured; every gate is followed by the same noise. The noise is amplitude damping, which drops $|1\rangle$ to $|0\rangle$ with probability $\gamma$, or a coherent over-rotation $e^{-i\varepsilon X/2}$ with $\varepsilon$ chosen so that both noises have the same average gate fidelity $F$. How does the survival probability, the chance of finding $|0\rangle$ at the end, decay with $m$?

Any qubit channel moves the Bloch vector as $\vec r\mapsto M\vec r+\vec c$. The 24 Cliffords act on $\vec r$ as the rotations of a cube, and averaging the noise over them, which the random sequence does, leaves a depolarizing map that shrinks $\vec r$ by $p=\operatorname{Tr}M/3=2F-1$. With the final noise left over, the mean survival is exactly

$$
P(m)=\frac12\left(1+c_z+M_{zz}\,p^m\right)=A\,p^m+B .
$$

For damping $M=\mathrm{diag}(\sqrt{1-\gamma},\sqrt{1-\gamma},1-\gamma)$ and $\vec c=(0,0,\gamma)$; for the rotation $\operatorname{Tr}M=1+2\cos\varepsilon$, so equal $F$ means $\cos\varepsilon=\sqrt{1-\gamma}-\gamma/2$.

$M$ and $\vec c$ of each noise come from its action on Bloch states. For each $\gamma$, both noises run 200 random sequences at each length $m$; the plot shows their mean survival with the spread of single sequences as bars, against $P(m)$ (dashed) and each noise's floor. The sphere shows, for one noise, where the qubit ends after each sequence, colored by $m$, the mean as a path from $|0\rangle$, and one sequence of ten gates drawn gate by gate, Clifford rotations in gray and the noise after each gate in red. Below it, a side cut through the ball counts where 2000 sequences of four lengths end:

```wl
rbCliffords = Select[Flatten[Outer[DiagonalMatrix[#1] . #2 &, Tuples[{1, -1}, 3], Permutations[IdentityMatrix[3]], 1], 1], Det[#] == 1 &];
rbBlochMap[f_] := With[{c = f[QuantumState["BlochVector"[{0, 0, 0}]]]["BlochVector"]}, {Transpose[f[QuantumState["BlochVector"[#]]]["BlochVector"] - c & /@ IdentityMatrix[3]], c}];
rbNoises[γ_] := Chop[Re[N[rbBlochMap[#]]]] & /@ <|"damping" -> QuantumChannel["AmplitudeDamping"[γ]], "over-rotation" -> QuantumOperator["RX"[ArcCos[Sqrt[1 - γ] - γ/2]]]|>;
rbUndo[seq_] := Transpose[Dot @@ Reverse[seq]];
rbEnd[{mm_, c_}, seq_] := Fold[mm . (#2 . #1) + c &, {0., 0., 1.}, Append[seq, rbUndo[seq]]];
rbExact[{mm_, c_}][m_] := (1 + c[[3]] + mm[[3, 3]] (Tr[mm]/3)^m)/2;
rbMs = Range[1, 61, 4]; rbCutMs = {1, 5, 15, 40, 60};
rbColor[m_] := ColorData["SunsetColors"][0.15 + 0.8 (m - 1)/60];
rbNoiseColor = <|"damping" -> RGBColor[0, 0.75, 0.85], "over-rotation" -> RGBColor[0.9, 0.45, 0.8]|>;
(* one gate drawn as the arc it rotates through *)
rbAxisAngle[r_] := With[{θ = ArcCos[Clip[(Tr[r] - 1)/2, {-1, 1}]]}, Which[θ < 10^-9, {{0, 0, 1}, 0.},
    Abs[θ - Pi] < 10^-9, {Normalize[First[NullSpace[N[r] - IdentityMatrix[3]]]], N[Pi]},
    True, {Normalize[{r[[3, 2]] - r[[2, 3]], r[[1, 3]] - r[[3, 1]], r[[2, 1]] - r[[1, 2]]}], N[θ]}]];
rbArc[r_, v_] := With[{aa = rbAxisAngle[r]}, Table[RotationMatrix[s Last[aa], First[aa]] . v, {s, 0, 1, 1/24}]];
rbNoiseStep[{mm_, c_}, v_] := If[Norm[c] < 10^-12 && Max[Abs[mm . Transpose[mm] - IdentityMatrix[3]]] < 10^-9, rbArc[mm, v], {v, mm . v + c}];
rbRun[noise_, seq_] := Rest[FoldList[With[{g = rbArc[#2, Last[#1[[2]]]]}, {g, rbNoiseStep[noise, Last[g]]}] &, {{}, {{0., 0., 1.}}}, Append[seq, rbUndo[seq]]]];
rbSeq = BlockRandom[SeedRandom[4]; RandomChoice[rbCliffords, 10]];
rbData[γ_] := BlockRandom[SeedRandom[1]; With[{maps = rbNoises[γ]}, <|"maps" -> maps,
     "ends" -> Map[Function[map, AssociationMap[Table[rbEnd[map, RandomChoice[rbCliffords, #]], 200] &, rbMs]], maps],
     "cuts" -> Map[Function[map, AssociationMap[Table[rbEnd[map, RandomChoice[rbCliffords, #]], 2000] &, rbCutMs]], maps]|>]];
rbSphere[d_, k_] := Show[ResourceFunction["BlochSpherePlot"][], Graphics3D[{
     KeyValueMap[{rbColor[#1], Opacity[0.4], AbsolutePointSize[2.5], Point[#2]} &, d["ends", k]],
     {rbNoiseColor[k], AbsoluteThickness[3.5], Line[Prepend[Mean /@ Values[d["ends", k]], {0., 0., 1.}]]},
     KeyValueMap[{rbColor[#1], AbsolutePointSize[8], Point[Mean[#2]]} &, d["ends", k]],
     {RGBColor[0.45, 0.52, 0.65], Opacity[0.55], AbsoluteThickness[1.5], Line[First[#]] & /@ rbRun[d["maps", k], rbSeq]},
     {RGBColor[0.8, 0.32, 0.32], Opacity[0.8], AbsoluteThickness[2.5], Line[Last[#]] & /@ rbRun[d["maps", k], rbSeq]}}],
   ViewPoint -> {2.2, -2.4, 1.2}, SphericalRegion -> True, ImageSize -> 175];
rbSurvival[d_] := Show[
   Plot[Evaluate[rbExact[d["maps", #]][m] & /@ {"damping", "over-rotation"}], {m, 0, 62}, PlotStyle -> Directive[GrayLevel[0.55], Dashed]],
   Graphics[KeyValueMap[Function[{k, shift}, KeyValueMap[Function[{m, pts}, With[{s = (1 + pts[[All, 3]])/2},
        {rbNoiseColor[k], AbsoluteThickness[1.5], Line[{{m + shift, Mean[s] - StandardDeviation[s]}, {m + shift, Mean[s] + StandardDeviation[s]}}]}]], d["ends", k]]],
     <|"damping" -> -0.6, "over-rotation" -> 0.6|>]],
   ListLinePlot[KeyValueMap[{#1, Mean[(1 + #2[[All, 3]])/2]} &, d["ends", #]] & /@ {"damping", "over-rotation"},
    PlotStyle -> (Directive[rbNoiseColor[#], AbsoluteThickness[3]] & /@ {"damping", "over-rotation"}), PlotMarkers -> {"\[FilledCircle]", 8}],
   PlotRange -> {{0, 63}, {0, 1.02}}, Frame -> True, Axes -> False, FrameLabel -> {"sequence length m", "survival P(0)"}, ImageSize -> 250, AspectRatio -> 0.8,
   GridLines -> {None, {{rbExact[d["maps", "damping"]][Infinity], Directive[rbNoiseColor["damping"], Dashed]}, {1/2, Directive[rbNoiseColor["over-rotation"], Dashed]}}}];
rbCut[d_, k_, m_] := With[{pts = d["cuts", k][m]}, With[{counts = BinCounts[BlockRandom[SeedRandom[9]; MapThread[{#2 Norm[#1[[;; 2]]], #1[[3]]} &, {pts, RandomChoice[{-1, 1}, Length[pts]]}]], {-1.02, 1.02, 0.04}, {-1.02, 1.02, 0.04}]},
    Show[ArrayPlot[Reverse[Transpose[Log10[1. + counts]]], DataRange -> {{-1.02, 1.02}, {-1.02, 1.02}}, ColorFunction -> (ColorData["SunsetColors"][Clip[#/3., {0, 1}]] &), ColorFunctionScaling -> False, Frame -> False],
     Graphics[{{GrayLevel[0.7], Circle[]}, {GrayLevel[0.5], Dashed, Line[{{0, -1.1}, {0, 1.1}}]}, {rbNoiseColor[k], PointSize[0.06], Point[{Norm[Mean[pts][[;; 2]]], Mean[pts][[3]]}]},
       Text[Style["|0⟩", 11, White], {0, 1.16}], Text[Style["|1⟩", 11, White], {0, -1.16}]}],
     PlotRange -> {{-1.2, 1.2}, {-1.25, 1.25}}, ImageSize -> 100, Background -> Black, PlotLabel -> Style[Row[{"m = ", m}], 12, FontFamily -> "Helvetica"]]]];
rbView[d_, k_] := Grid[{
    {rbSphere[d, k], Legended[rbSurvival[d], Placed[LineLegend[Join[rbNoiseColor /@ {"damping", "over-rotation"}, {Directive[GrayLevel[0.55], Dashed]}, Directive[#, Dashed] & /@ (rbNoiseColor /@ {"damping", "over-rotation"})],
        {"damping: mean of 200 sequences, bars: their spread", "over-rotation: the same", "exact P(m) = A p^m + B", "damping floor (1 + γ)/2", "over-rotation floor 1/2"}], Right]]},
    {Column[{Style[Row[{k, ": side cut through the ball, where 2000 runs of length m end (log color)"}], 12, FontFamily -> "Helvetica"], Row[rbCut[d, k, #] & /@ rbCutMs, Spacer[6]]}, Center], SpanFromLeft}},
   Spacings -> {2, 1.5}, Alignment -> Center];
rbViews = Association[Flatten[Table[With[{d = rbData[γ]}, ({γ, #} -> rbView[d, #]) & /@ {"damping", "over-rotation"}], {γ, {1/50, 1/20, 1/5}}]]];
With[{rbViews = rbViews},
 Manipulate[rbViews[{rbγ, rbK}], Grid[{{Control[{{rbγ, 1/20, "noise strength γ"}, {1/50, 1/20, 1/5}, SetterBar}],
     Control[{{rbK, "damping", "on the sphere"}, {"damping", "over-rotation"}, SetterBar}]}}, Spacings -> 4]]]
```

The two noises lose the qubit at the same rate: averaging over random Cliffords keeps only $\operatorname{Tr}M$ of the noise, so both means follow the same $p^m$. Everything else differs. Damping shrinks every run into the ball and pulls it toward $|0\rangle$, so the runs stay together near the axis and settle at the floor $(1+\gamma)/2$, above one half. The over-rotation keeps every run pure, on the surface; the random gates turn its small rotations to random axes, so the runs spread over the sphere and only their average sinks to the center, at the floor $\frac12$, with a far larger spread from one sequence to the next. The decay rate that randomized benchmarking reports cannot tell a coherent error from a random one; the spread between sequences can.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Dephasing: GHZ gains nothing, squeezing a bounded factor

A frequency $\omega$ is written on $n$ qubits, each dephasing independently at rate $\gamma$, by Ramsey runs of length $t$ repeated over a total time $T$. Without noise a GHZ probe has quantum Fisher information $F_Q=n^2t^2$, against $F_Q=nt^2$ for the product probe $|+\rangle^{\otimes n}$; $F_Q$ sets the best precision any measurement can reach, $\delta\omega\ge1/\sqrt{\nu F_Q}$ for $\nu$ runs. Does the GHZ advantage survive once $t$ is chosen optimally?

With $\hbar=1$ and $\mathcal D[\hat c]\rho=\hat c\rho\hat c^\dagger-\frac12\{\hat c^\dagger\hat c,\rho\}$, each run evolves under

$$
\dot\rho=-i\left[\frac\omega2\sum_{k=1}^nZ_k,\rho\right]+\gamma\sum_{k=1}^n\mathcal D[Z_k]\rho ,
$$

which damps each qubit's coherence as $e^{-2\gamma t}$. The $T/t$ runs carry information $(T/t)\,F_Q(t)$, with $F_Q=2\sum_{\lambda_i+\lambda_j>0}|\langle i|\partial_\omega\rho|j\rangle|^2/(\lambda_i+\lambda_j)$ over the eigenpairs $(\lambda_i,|i\rangle)$ of $\rho(t)$; the noise commutes with the phase, so $F_Q$ does not depend on $\omega$ and is taken at $\omega=0$.

The GHZ coherence decays as $e^{-2n\gamma t}$, so $F_Q^\mathrm{GHZ}=n^2t^2e^{-4n\gamma t}$ and $F_Q^\mathrm{product}=nt^2e^{-4\gamma t}$, and the best run lengths $t=1/(4n\gamma)$ and $t=1/(4\gamma)$ give

$$
\max_t\frac{F_Q^\mathrm{GHZ}}{t}=\max_t\frac{F_Q^\mathrm{product}}{t}=\frac{n}{4e\gamma}.
$$

Under the phase and the dephasing a GHZ state stays within $|0\cdots0\rangle$ and $|1\cdots1\rangle$, where it behaves as one qubit whose $\omega$ and $\gamma$ are both multiplied by $n$, while the product probe is $n$ independent qubits. So $F_Q$ is computed for symbolic $n$ from one evolved qubit. For comparison, a spin‑squeezed probe is added: the one‑axis‑twisted state with the least spin noise, read out through $J_y$, whose dephased spin averages give a lower bound on its $F_Q$. The curves show $F_Q/(nt)$, the information per qubit per unit time, for $n=2,10,100$, with markers at the best run lengths of GHZ and the product probe. The dotted curve is the large‑$n$ limit for any state, $F_Q\le n\,t^2e^{-4\gamma t}/(1-e^{-4\gamma t})$ (Escher et al. 2011; Demkowicz‑Dobrzański et al. 2012); rates are in units of $\gamma$:

```wl
dgzRho1[w_, g_] := Simplify[Normal[QuantumEvolve[w/2 QuantumOperator["Z"], {QuantumOperator["Z"]} -> {g}, QuantumState["+"], t]["DensityMatrix"]], g > 0 && n > 0];
dgzQFI[ρ_] := With[{es = Eigensystem[ρ /. ω -> 0]}, {λ = First[es], V = Orthogonalize[Last[es]]}, {M = Conjugate[V] . (D[ρ, ω] /. ω -> 0) . Transpose[V]},
   Simplify[2 Total[MapThread[If[PossibleZeroQ[#2], 0, Abs[#1]^2/#2] &, {M, Outer[Plus, λ, λ]}, 2], 2], γ > 0 && t > 0 && n > 0]];
dgzAny = {dgzQFI[dgzRho1[n ω, n γ]], n dgzQFI[dgzRho1[ω, γ]]};
dgzBest = Function[f, With[{t0 = t /. First[Solve[D[f/t, t] == 0, t]]}, {Simplify[f/t /. t -> t0], t0}]] /@ (dgzAny /. γ -> 1);
dgzSpin[k_] := With[{m = Range[k/2, -k/2, -1]}, With[{jp = SparseArray[Band[{1, 2}] -> (Sqrt[k/2 (k/2 + 1) - # (# + 1)] & /@ Rest[m]), {k + 1, k + 1}]},
   N[{(jp + ConjugateTranspose[jp])/2, (jp - ConjugateTranspose[jp])/(2 I), SparseArray[DiagonalMatrix[m]]}]]];
dgzTwist[J_, x0_, mu_?NumericQ] := With[{v = Exp[-I mu/2 Diagonal[Normal[J[[3]]]]^2] x0}, With[{ev = Function[op, Conjugate[v] . op . v]},
   With[{y = ev[J[[2]]], z = ev[J[[3]]], c = ev[(J[[2]] . J[[3]] + J[[3]] . J[[2]])/2]},
    {Re[ev[J[[1]]]], Min[Eigenvalues[Re[{{ev[J[[2]] . J[[2]]] - y^2, c - y z}, {c - y z, ev[J[[3]] . J[[3]]] - z^2}}]]]}]]];
dgzSqueezed[k_] := With[{J = dgzSpin[k], x0 = N[Sqrt[Binomial[k, Range[0, k]]]/2^(k/2)]},
   Block[{xi2}, xi2[mu_?NumericQ] := k #[[2]]/#[[1]]^2 &[dgzTwist[J, x0, mu]];
    dgzTwist[J, x0, mu /. Last[NMinimize[{xi2[mu], 10^-4 <= mu <= 1}, mu]]]]];
dgzSqF[{jx_, v_}, k_] := With[{η = Exp[-2 γ t]}, η^2 t^2 jx^2/(η^2 v + (1 - η^2) k/4)];
dgzNs = {2, 10, 100}; dgzSqNs = {10, 100};
dgzSq = dgzSqF[dgzSqueezed[#], #] & /@ dgzSqNs;
dgzLimit = t Exp[-4 γ t]/(1 - Exp[-4 γ t]);
Labeled[LogLinearPlot[Evaluate[Join[(First[dgzAny]/(n t) /. {γ -> 1, n -> #}) & /@ dgzNs, {Last[dgzAny]/(n t) /. γ -> 1},
     MapThread[#1/(#2 t) /. γ -> 1 &, {dgzSq, dgzSqNs}], {dgzLimit /. γ -> 1}]], {t, 10^-4, 3},
   PlotRange -> {0, 0.27}, PlotStyle -> {ColorData[97][1], ColorData[97][2], ColorData[97][3], Gray,
     Directive[ColorData[97][2], Dashed], Directive[ColorData[97][3], Dashed], Directive[Gray, Dotted]}, Frame -> True,
   FrameLabel -> {"γ t", Row[{"γ ", Subscript["F", "Q"], "/(n t)"}]}, GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.5, 0.35], Thin],
   Epilog -> {PointSize[0.014], MapIndexed[{ColorData[97][First[#2]], Point[{Log[1/(4 #1)], 1/(4 E)}]} &, dgzNs], Gray, Point[{Log[1/4], 1/(4 E)}]},
   PlotLegends -> Join[Row[{"GHZ, n = ", #}] & /@ dgzNs, {"product, any n"}, Row[{"spin‐squeezed, n = ", #}] & /@ dgzSqNs, {"limit for any state, large n"}], ImageSize -> 620],
  Style[Row[{"any n:   GHZ largest γ ", Subscript["F", "Q"], "/t = ", TraditionalForm[dgzBest[[1, 1]]], " at γ t = ", TraditionalForm[dgzBest[[1, 2]]],
    ",   product largest γ ", Subscript["F", "Q"], "/t = ", TraditionalForm[dgzBest[[2, 1]]], " at γ t = ", TraditionalForm[dgzBest[[2, 2]]]}], 13, FontFamily -> "Helvetica"], Top]
```

Divided by $n$, the product curve is the same for every $n$, and each GHZ curve reaches exactly the same peak, only at a run $n$ times shorter: GHZ picks up the phase $n$ times faster but loses its coherence $n$ times faster, and at the best run length the two cancel. The squeezed probe does beat the product probe, and by more as $n$ grows, but it stays under the dotted limit, which lies a factor $e$ above the product peak. Under independent dephasing, entanglement buys at most a fixed factor, never the $n$ of the noiseless case.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Circuits and algorithms

The examples in this section are quantum circuits that QF builds from named gates and runs exactly: teleportation with classical records, search and factoring, a variational algorithm, a Trotter circuit and a fixed-point arithmetic circuit. Where a claim holds for every angle or every input, WL proves it by exact maximization or quantifier elimination.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Teleportation with the classical channel as record wires

Alice holds an unknown qubit and one half of a Bell pair; Bob holds the other half. She measures her two qubits in the Bell basis, and Bob corrects his qubit according to her two outcomes. The outcomes are kept as record qubits, so Bob's corrections are gates controlled by the records. Does Bob end up with Alice's state, for every input?

With the input $|\psi\rangle=\cos\frac\beta2|0\rangle+e^{i\phi}\sin\frac\beta2|1\rangle$ on qubit 1, whose Bloch vector is $(\sin\beta\cos\phi,\sin\beta\sin\phi,\cos\beta)$, and $|00\rangle$ on qubits 2 and 3, the protocol is

$$
\underbrace{\mathrm{CNOT}_{23}\,H_2}_{\text{Bell pair}},\qquad \underbrace{H_1\,\mathrm{CNOT}_{12}}_{\text{Bell basis}},\qquad \text{measure qubits 1, 2}\to m_1,m_2,\qquad \text{Bob applies } Z^{m_1}X^{m_2}.
$$

The protocol runs for general $\beta$ and $\phi$ with both record-controlled corrections, with one of them, and with none. Bob's Bloch vector comes out in closed form, and the sphere evaluates it at twelve random inputs: a line joins each of Alice's inputs to the state Bob ends with, and the control chooses which corrections Bob applies:

```wl
teleportLabel[s_] := Style[s, Black, 12, FontFamily -> "Helvetica"];
teleportCircuit[corrections_] := QuantumCircuitOperator[Join[{
    QuantumCircuitOperator[{"H" -> 2, "CNOT" -> {2, 3}}, teleportLabel["Bell pair"]], "Barrier",
    QuantumCircuitOperator[{"CNOT" -> {1, 2}, "H" -> 1}, teleportLabel["Alice"]], {1}, {2}, "Barrier"},
   If[corrections === {}, {}, {QuantumCircuitOperator[corrections, teleportLabel["Bob"]]}]]];
bobBloch[corrections_] := FullSimplify[QuantumPartialTrace[Replace[
      teleportCircuit[corrections][
       QuantumTensorProduct[QuantumState[{Cos[β/2], Exp[I ϕ] Sin[β/2]}], QuantumState["00"]]], m_QuantumMeasurement :> m["State"]],
     {1, 2, 3, 4}]["BlochVector"], {β, ϕ} ∈ Reals];
teleportWires = {-1 -> Placed["m₂ to Bob", Right], 0 -> Placed["m₁ to Bob", Right], 1 -> Placed["Alice  |ψ⟩", Left], 2 -> Placed["Alice", Left], 3 -> Placed["Bob", Right]};
teleportCases = <|"X and Z" -> {"CNOT" -> {-1, 3}, "CZ" -> {0, 3}}, "X only" -> {"CNOT" -> {-1, 3}}, "Z only" -> {"CZ" -> {0, 3}}, "none" -> {}|>;
teleportInputs = BlockRandom[SeedRandom[11]; Table[{ArcCos[RandomReal[{-1, 1}]], RandomReal[{0, 2 π}]}, 12]];
teleportAlice = ({Sin[#1] Cos[#2], Sin[#1] Sin[#2], Cos[#1]} &) @@@ teleportInputs;
teleportSphere[bob_] := With[{bobs = (bob /. {β -> #1, ϕ -> #2} &) @@@ teleportInputs},
   Show[ResourceFunction["BlochSpherePlot"][], Graphics3D[{
      {GrayLevel[0.45], AbsoluteThickness[1.5], Line[Transpose[{teleportAlice, bobs}]]},
      {ColorData[97][2], Opacity[0.55], AbsolutePointSize[15], Point[teleportAlice]},
      {ColorData[97][1], AbsolutePointSize[8], Point[bobs]}}], ImageSize -> 260]];
teleportViews = AssociationThread[Keys[teleportCases], MapThread[Function[{corr, bob}, Column[{
       Grid[{{(teleportCircuit[corr]["Diagram", "WireLabels" -> teleportWires, ImageSize -> 520] //. EventHandler[\[FormalG]_, ___] :> \[FormalG]), teleportSphere[bob]}}, Spacings -> 2],
       Nothing}, Center]],
    {Values[teleportCases], bobBloch /@ Values[teleportCases]}]];
teleportLegend = legendBelow[Column[{PointLegend[{Directive[ColorData[97][2], Opacity[0.55]], ColorData[97][1]}, {"Alice's input", "Bob's state"}, LegendMarkers -> {{"\[FilledCircle]", 15}, {"\[FilledCircle]", 8}}],
    LineLegend[{Directive[GrayLevel[0.45], AbsoluteThickness[1.5]]}, {"Alice to Bob"}]}]];
With[{teleportViews = teleportViews, teleportLegend = teleportLegend, teleportCases = teleportCases},
 Manipulate[Legended[teleportViews[tpCase], teleportLegend],
 {{tpCase, "X and Z", "corrections Bob applies"}, Keys[teleportCases], SetterBar}]]
```

With both corrections Bob's Bloch vector equals the input's for every $\beta$ and $\phi$. With only one correction, the Pauli Bob did not undo flips two Bloch components at random, so only the component it leaves unflipped survives and his dots collapse onto one axis. With neither, Bob holds the maximally mixed state at the center of the sphere.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Noise caps Grover's speedup

A Boolean oracle accepts one of the $N=2^n$ bit strings, and Grover's search finds it in about $\frac\pi4\sqrt N$ iterations, where guessing needs $N$ tries. If every qubit is phase-flipped with probability $p$ before each iteration, how much noise can the search take before the speedup is gone?

With $|s\rangle=|+\rangle^{\otimes n}$, $x_*$ the accepted string and $\mathcal E_p(\rho)=(1-p)\rho+pZ\rho Z$ on each qubit, the iteration, the state after $k$ iterations and the success probability are

$$
G=\left(2|s\rangle\langle s|-1\right)O,\qquad \rho_{k+1}=G\,\mathcal E_p^{\otimes n}(\rho_k)\,G^\dagger,\quad \rho_0=|s\rangle\langle s|,\qquad P_k=\langle x_*|\rho_k|x_*\rangle .
$$

A run of $k$ iterations, checked with one more oracle call and repeated until it succeeds, costs on average $(k+1)/P_k$ oracle calls; the best $k$ gives the cost of the noisy search, against $N$ for guessing.

$G$ is QF's Grover iteration for the oracle accepting $x_*=5$. Phase flips on all qubits multiply each element $\rho_{xy}$ by $(1-2p)^{d(x,y)}$, with $d(x,y)$ the number of bits in which $x$ and $y$ differ, the same map as QF's `"PhaseFlip"` channel on every qubit. Left: $P_k$ for $N=256$ at four noise levels; right: the expected oracle calls against $p$ for $N=64$ and $256$, with guessing dashed:

```wl
grnIter[n_] := N[Normal[QuantumCircuitOperator["GroverPhase"[BooleanFunction[2^5, n]]]["Matrix"]]];
grnMask[n_, p_] := Outer[(1 - 2 p)^DigitCount[BitXor[#1, #2], 2, 1] &, Range[0, 2^n - 1], Range[0, 2^n - 1]];
(* success probability after k = 0, 1, ... iterations, every qubit phase-flipped with probability p before each iteration *)
grnSuccess[n_, p_] := With[{g = grnIter[n], mask = grnMask[n, p]},
   Re[#[[6, 6]]] & /@ NestList[g . (mask #) . ConjugateTranspose[g] &, ConstantArray[2.^-n, {2^n, 2^n}], Ceiling[π Sqrt[2^n]]]];
(* expected oracle calls: k iterations plus one call to check, repeated until success *)
grnCost[n_, p_] := Min[Rest[MapIndexed[First[#2]/#1 &, grnSuccess[n, p]]]];
grnPs = {0, 0.003, 0.01, 0.03};
grnLeft = ListPlot[grnSuccess[8, #] & /@ grnPs, DataRange -> {0, Ceiling[π Sqrt[256]]}, Joined -> True, PlotMarkers -> {Automatic, 5}, PlotRange -> {0, 1.02},
   Frame -> True, FrameLabel -> {"iterations k", "success probability"}, GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.5, 0.35], Thin],
   PlotLegends -> Placed[Row[{"p = ", #}] & /@ grnPs, Below], PlotLabel -> "N = 256, one accepted string", ImageSize -> 360];
grnPgrid = Append[10^Range[-4., -0.375, 0.125], 0.49];
grnCosts = Table[{p, grnCost[n, p]}, {n, {6, 8}}, {p, grnPgrid}];
grnRight = ListLogLogPlot[grnCosts, Joined -> True, PlotMarkers -> {Automatic, 5}, PlotRange -> {{10^-4, 0.49}, {3, 2000}},
   Frame -> True, FrameLabel -> {"phase‐flip probability p per qubit per iteration", "expected oracle calls"},
   GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.5, 0.35], Thin],
   Epilog -> {Dashed, ColorData[97][1], Line[{{Log[10^-4], Log[64]}, {Log[0.49], Log[64]}}], ColorData[97][2], Line[{{Log[10^-4], Log[256]}, {Log[0.49], Log[256]}}]},
   FrameTicks -> {{Automatic, None}, {{#, ToString[#]} & /@ {0.0001, 0.001, 0.01, 0.1}, None}}, PlotLegends -> Placed[{"Grover, N = 64", "Grover, N = 256"}, Below], PlotLabel -> "dashed: guessing, N calls", ImageSize -> 360];
Grid[{{grnLeft, grnRight}}, Spacings -> 3, Alignment -> Top]
```

Noise lowers each peak and pulls it to fewer iterations, so the best run stops earlier and has to be repeated more often. While a whole search, $n$ qubits over about $\frac\pi4\sqrt N$ iterations, collects well under one phase flip, the cost stays near $\frac\pi4\sqrt N$; once it collects about one, the cost climbs toward the classical $N$. So the noise a search tolerates per qubit and iteration while keeping its noiseless cost is roughly one over the number of chances it gives a phase flip, $1/(n\sqrt N)$: the larger the search, the better the qubits must be. A search that must run at least one iteration loses to guessing outright only at heavy noise, $p$ of order $0.2$ to $0.3$; running none is guessing itself.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Shor's algorithm factors 15 and 21

To factor $N$, pick $a$ coprime to $N$ and find its order $r$, the smallest $r$ with $a^r\equiv1 \pmod N$. If $r$ is even and $a^{r/2}\not\equiv-1$, then $\gcd(a^{r/2}\pm1,N)$ are factors. Can phase estimation find $r$ from a single distribution?

The work register holds $y\in\{0,\dots,N-1\}$ in as many qubits as $N$ has bits, and the multiplier $U_a|y\rangle=|ay\bmod N\rangle$ (identity on $y\ge N$) is unitary because $a$ is invertible mod $N$. Starting from $|1\rangle=\frac1{\sqrt r}\sum_{k=0}^{r-1}|u_k\rangle$, an equal superposition of eigenvectors with $U_a|u_k\rangle=e^{2\pi ik/r}|u_k\rangle$, phase estimation with $t$ control qubits returns $x$ with $x/2^t\approx k/r$, each $k$ with probability $1/r$. When $r$ divides $2^t$ the peaks are exact; otherwise each spreads over its neighbors. The outcome nearest $k/r$ lies within $2^{-t-1}$ of it, and since $2^t>N^2$ its continued-fraction convergent with denominator below $N$ is $k/r$ in lowest terms, so $r$ is the least common multiple of those denominators; outcomes further out carry no such guarantee.

Phase estimation runs with its initial $X$ gates removed and the work register in $|1\rangle$, with $t=8$ for $N=15$ and $t=9$ for $N=21$; each outcome of probability above $0.05$ is labeled with its continued-fraction convergent $k/r$, and each title gives the order read from them and the two greatest common divisors:

```wl
shorU[nN_, a_] := With[{d = 2^BitLength[nN]}, QuantumOperator[SparseArray[Thread[Transpose[{Join[Mod[a Range[0, nN - 1], nN], Range[nN, d - 1]], Range[0, d - 1]}] + 1 -> 1], {d, d}]]];
shorProbs[nN_, a_, t_] := N[QuantumCircuitOperator["PhaseEstimation"[shorU[nN, a], t]][[BitLength[nN] + 1 ;;]]][
    QuantumState[StringRepeat["0", t] <> IntegerString[1, 2, BitLength[nN]]]]["ProbabilitiesList"];
shorPanel[nN_, a_, t_] := With[{p = shorProbs[nN, a, t]}, {peaks = Pick[Range[0, 2^t - 1], UnitStep[p - 0.05], 1]},
   {fr = Last[Select[Convergents[#/2^t], Denominator[#] < nN &]] & /@ peaks}, {r = LCM @@ Denominator[fr], pts = Transpose[{peaks/2^t, p[[peaks + 1]], fr}]},
   ListPlot[{Transpose[{Range[0, 2^t - 1]/2^t, p}], Select[pts, Denominator[Last[#]] == r &][[All, ;; 2]], Select[pts, Denominator[Last[#]] < r &][[All, ;; 2]]},
    Filling -> {1 -> Axis}, PlotMarkers -> {{"\[FilledCircle]", 2}, {"\[FilledCircle]", 10}, {"\[EmptyCircle]", 10}}, PlotRange -> {{-0.03, 1.03}, {0, 0.62}}, Frame -> True,
    Epilog -> (Text[Style[TraditionalForm[#3], 11], {#1, #2 + 0.07}] & @@@ pts),
    FrameLabel -> {Row[{"x / ", Superscript[2, t]}], "P(x)"}, FrameTicks -> {{Automatic, None}, {{#, Style[TraditionalForm[#], FractionBoxOptions -> {Beveled -> True}]} & /@ Range[0, 1, 1/4], None}},
    PlotLabel -> Row[{"N = ", nN, ",  a = ", a, ",  r = ", r, If[EvenQ[r], Row[{",  gcd = ", GCD[a^(r/2) - 1, nN], ", ", GCD[a^(r/2) + 1, nN]}], " (odd)"]}], ImageSize -> 300]];
Column[{Labeled[QuantumCircuitOperator["PhaseEstimation"[shorU[15, 2], 4]][[BitLength[15] + 1 ;;]]["Diagram", ImageSize -> 560], Style["order finding for N = 15, a = 2, drawn with 4 of the 8 control qubits, and the four-qubit multiplier", 13, FontFamily -> "Helvetica"], Top],
  Legended[PlotGrid[{shorPanel[15, #, 8] & /@ {2, 11, 14}, shorPanel[21, #, 9] & /@ {2, 4, 13}}, ImageSize -> 900],
  legendBelow[PointLegend[ColorData[97] /@ {2, 3}, {"continued fraction gives r", "gives only a divisor of r"}, LegendMarkers -> {{"\[FilledCircle]", 10}, {"\[EmptyCircle]", 10}}]]]}, Center, Spacings -> 2]
```

For $N=15$ every distribution is $r$ equal peaks at the multiples of $1/r$; for $N=21$ a peak at a $k/r$ that falls between grid points spreads over neighboring $x$, and the continued fractions of the outcomes above the cutoff still return the order. The factors follow whenever $r$ is even and $a^{r/2}\not\equiv-1$: $a=14\equiv-1 \pmod{15}$ gives only trivial divisors and $a=4$ has odd order mod 21, so in both cases the algorithm restarts with another $a$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### QAOA on MaxCut: the exact landscape and a certified optimum

Color the $n$ vertices of a graph with two colors so that as many edges as possible join different colors; take the triangle, the complete graph $K_4$ and rings of 5 and 6 vertices. Each layer of the quantum approximate optimization algorithm (QAOA) applies a cost rotation by an angle $\gamma$ and a mixing rotation by an angle $\theta$, so $p$ layers prepare a state that depends on $2p$ angles, scored by the expected number of cut edges. An optimizer searching over the angles finds good values but cannot prove that no better ones exist. Can the best score be certified exactly, and how many layers reach the true maximum cut?

With $Z_i$ and $X_i$ the Pauli operators on vertex $i$, the cost Hamiltonian, the $p$-layer state and its score are

$$
\hat H_C=\sum_{(i,j)\in\mathcal E}\frac{1-Z_iZ_j}{2},\qquad |\psi_p\rangle=U_X(\theta_p)U_C(\gamma_p)\cdots U_X(\theta_1)U_C(\gamma_1)\,|+\rangle^{\otimes n},\qquad E=\langle\psi_p|\hat H_C|\psi_p\rangle,
$$

with $\mathcal E$ the edge set, $U_C(\gamma)=\prod_{(i,j)\in\mathcal E}e^{-i\gamma Z_iZ_j/2}$ and $U_X(\theta)=\prod_i e^{-i\theta X_i/2}$.

$E$ is in closed form from each graph's edge list and is drawn as a fraction of the maximum cut on one color scale. At one layer the figure is the whole landscape over $(\gamma_1,\theta_1)$ with its exact maximum as a dot, and the title gives that maximum, with its minimal polynomial when it is irrational; with more layers it is a slice over the last layer's angles, the earlier ones held at their best values, and the title gives the best value found:

```wl
qaoaCircuit[graph_, p_] := QuantumCircuitOperator[Join[{"H" -> VertexList[graph]},
    Catenate[Table[Join[Thread["R"[Subscript[γ, j], "ZZ"] -> List @@@ EdgeList[graph]], {"RX"[Subscript[θ, j]] -> VertexList[graph]}], {j, p}]]]];
qaoaLandscape[graph_, p_] := With[{hc = Total[(1 - QuantumOperator["ZZ", #])/2 & /@ List @@@ EdgeList[graph]], ψ = qaoaCircuit[graph, p][]},
    ψ["Dagger"][hc[ψ]]["Scalar"]];
qaoaAngles[p_] := Flatten[Table[{Subscript[γ, j], Subscript[θ, j]}, {j, p}]];
qaoaBest[l_, 1] := Maximize[FullSimplify[l, Subscript[γ, 1] ∈ Reals && Subscript[θ, 1] ∈ Reals], qaoaAngles[1]];
qaoaBest[l_, p_] := BlockRandom[SeedRandom[1]; MaximalBy[Table[FindMaximum[l, Transpose[{qaoaAngles[p], RandomReal[{0, 2 π}, 2 p]}], AccuracyGoal -> 6, PrecisionGoal -> 6], 6], First][[1]]];
qaoaColor = ColorData["SunsetColors"];
qaoaTitle[best_, cut_, 1] := If[First[best] ∈ Rationals, {Row[{"max E = ", First[best]}]},
    {Row[{"max E = ", N[First[best], 4]}], Row[{"root of ", TraditionalForm[MinimalPolynomial[First[best], \[FormalX]] == 0]}]}];
qaoaTitle[best_, cut_, p_] := {Row[{"best E found ≈ ", NumberForm[First[best], {6, 4}], ", maximum cut ", cut}]};
qaoaView[graph_, p_] := With[{l = qaoaLandscape[graph, p], cut = First[FindMaximumCut[graph]]}, With[{best = qaoaBest[l, p]},
    Grid[{{(qaoaCircuit[graph, p]["Diagram", ImageSize -> 220 + 160 p] //. EventHandler[\[FormalG]_, ___] :> \[FormalG]),
       DensityPlot[Evaluate[(l /. Drop[Last[best], -2] /. {Subscript[γ, p] -> u, Subscript[θ, p] -> v})/cut], {u, 0, 2 π}, {v, 0, π}, PlotPoints -> 40,
        ColorFunction -> qaoaColor, ColorFunctionScaling -> False, PlotRange -> {0, 1}, FrameLabel -> {Subscript["γ", p], Subscript["θ", p]},
        Epilog -> {ColorData[97][1], PointSize[Large], Point[Mod[N[{Subscript[γ, p], Subscript[θ, p]} /. Last[best]], {2 π, π}]]},
        PlotLabel -> Column[Prepend[qaoaTitle[best, cut, p], Row[{Graph[graph, ImageSize -> 40], "   max cut = ", cut, ",   p = ", p}]], Alignment -> Center],
        ImageSize -> 280]}}, Spacings -> 3, Alignment -> Center]]];
qaoaGraphs = <|"triangle" -> {CycleGraph[3], 1}, Subscript["K", 4] -> {CompleteGraph[4], 2}, "ring of 5" -> {CycleGraph[5], 2}, "ring of 6" -> {CycleGraph[6], 3}|>;
qaoaViews = Association[KeyValueMap[Function[{name, gp}, Table[{name, p} -> qaoaView[First[gp], p], {p, Last[gp]}]], qaoaGraphs]];
With[{qaoaViews = qaoaViews, qaoaGraphs = qaoaGraphs, qaoaColor = qaoaColor},
 Manipulate[Legended[qaoaViews[{qaG, Min[qaP, Last[qaoaGraphs[qaG]]]}], Placed[BarLegend[{qaoaColor, {0, 1}}, LegendLabel -> "E / max cut"], Right]],
 Row[{Control[{{qaG, "triangle", "graph"}, Keys[qaoaGraphs], SetterBar}], Spacer[30], "layers p  ", Dynamic[SetterBar[Dynamic[qaP], Range[Last[qaoaGraphs[qaG]]]]]}],
 {{qaP, 1}, ControlType -> None}]]
```

Because the expected cut is an exact trigonometric formula in the angles, WL's `Maximize` returns its one-layer global maximum as an exact algebraic number: no pair of angles anywhere does better. With more layers the best value comes from local searches started at random angles. One layer reaches the maximum cut only on the triangle; on $K_4$ it stops at an algebraic number below it, and on both rings at $3/4$ of the edges. A $p$-layer circuit sees only $2p+1$ consecutive edges of a ring, so, as Farhi, Goldstone and Gutmann found, it cuts $(2p+1)/(2p+2)$ of the edges until it sees the whole ring: $5/6$ for the ring of 6 at two layers, and the full cut once $2p+1$ reaches the ring's length, at two layers for the ring of 5 and three for the ring of 6; $K_4$ needs two. A search that reaches the maximum cut to numerical precision has found the best score, since no state cuts more edges than the maximum.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### An infinite Ising chain from 2d − 1 sites: the light cone of a Trotter circuit

A chain of spins with coupling $J$ in a transverse field $h$ starts with every spin up, $|0\cdots0\rangle$, and evolves under the transverse-field Ising Hamiltonian, Trotterized with step $\delta t$. What is the magnetization $\langle Z_c\rangle$ of a site deep inside the chain after $d$ steps?

With $\hbar=1$, $\hat H=J\sum_iZ_iZ_{i+1}+h\sum_iX_i$, and one step is $\prod_iR_X^{(i)}(2h\,\delta t)\prod_iR_{ZZ}^{(i,i+1)}(2J\,\delta t)$ with $R_{ZZ}(\vartheta)=e^{-i\vartheta Z\otimes Z/2}$. After $d$ steps $\langle Z_c\rangle=\langle0|U^\dagger Z_cU|0\rangle$, unchanged when the sign of $\hat H$ flips, since $\hat H$ and the initial state are real, so the flip only complex-conjugates $U$; so the sign convention, opposite to that of the Ising-ring example, does not matter here. Every gate acts on one site or on two neighbors, so each step lets $Z_c$ reach one more site to each side; gates beyond that reach cancel against their conjugates. This light cone would span the $2d+1$ sites from $c-d$ to $c+d$, but the first $ZZ$ layer only multiplies $|0\cdots0\rangle$ by a phase, so it spans only the $2d-1$ sites from $c-d+1$ to $c+d-1$, and a window of that size gives the value for a chain of any length.

The circuit runs from $|0\cdots0\rangle$ on a window of $m$ sites centered on $c$ and gives $\langle Z_c\rangle$ for every odd $m$ up to $2d+3$, with $J\delta t=0.2$; the dashed lines mark $m=2d-1$ and the value of the largest window. The drawing shows the $2d+3$-site circuit, one box per Trotter step, with each gate colored by whether it reaches $Z_c$ (key above the drawing); the third control picks a window, framed on the circuit and circled on the plot:

```wl
lcJdt = 0.2;
lcTerms[m_, hJ_] := Join[QuantumOperator["ZZ", #] & /@ Partition[Range[m], 2, 1], hJ QuantumOperator["X", {#}] & /@ Range[m]];
lcTrotter[m_, d_, hJ_] := QuantumCircuitOperator["Trotterization"[lcTerms[m, hJ], 1, d, lcJdt d]];
lcZ[m_, d_, hJ_] := With[{st = QuantumCircuitOperator[{StringRepeat["0", m], lcTrotter[m, d, hJ]}][]}, Re[st["Dagger"][QuantumOperator["Z", {(m + 1)/2}][st]]["Scalar"]]];
lcStep[d_, j_, hJ_] := With[{m = 2 d + 3, c = d + 2, k = d - j + 1}, QuantumCircuitOperator[Join[
     QuantumOperator["R"[2 lcJdt, "ZZ"], {#, # + 1}, "Label" -> Which[j == 1, "phase", Max[Abs[# - c], Abs[# + 1 - c]] <= k, "in", True, "out"]] & /@ Join[Range[1, m - 1, 2], Range[2, m - 1, 2]],
     {"Barrier"}, QuantumOperator["RX"[2 hJ lcJdt], {#}, "Label" -> If[Abs[# - c] <= k - 1, "in", "out"]] & /@ Range[m]], "step " <> ToString[j]]];
lcCone[{hJ_, d_}] := With[{m = 2 d + 3, c = d + 2}, (QuantumCircuitOperator[lcStep[d, #, hJ] & /@ Range[d]]["Diagram", "ShowGateLabels" -> False, "BarrierStyle" -> Opacity[0],
    "GateBackgroundStyle" -> {"in" -> RGBColor[0.55, 0.45, 0.85], "out" -> GrayLevel[0.93], "phase" -> RGBColor[1, 0.87, 0.65]},
    "GateBoundaryStyle" -> {"in" -> RGBColor[0.4, 0.3, 0.7], "out" -> GrayLevel[0.8], "phase" -> RGBColor[0.9, 0.6, 0.2]},
    "WireLabels" -> Table[i -> Placed[Which[i == c, Style["c", Bold], i < c, "c \[Minus] " <> ToString[c - i], True, "c + " <> ToString[i - c]], Left], {i, m}], ImageSize -> {UpTo[520], UpTo[360]}] //. EventHandler[\[FormalG]_, ___] :> \[FormalG])];
lcWin = AssociationMap[Function[hd, Table[{m, lcZ[m, hd[[2]], hd[[1]]]}, {m, 1, 2 hd[[2]] + 3, 2}]], Tuples[{{0.5, 1, 2}, Range[2, 5]}]];
lcCones = AssociationMap[lcCone, Tuples[{{0.5, 1, 2}, Range[2, 5]}]];
lcColumn[i_] := If[OddQ[i], {0.625, 1.375}, {1.625, 2.375}];
lcWindow[hJ_, d_, m_] := With[{lo = d + 2 - (m - 1)/2, hi = d + 2 + (m - 1)/2, x = First[PlotRange /. AbsoluteOptions[lcCones[{hJ, d}], PlotRange]], p = 0.08},
   With[{xl = x[[1]] + 0.1, xr = x[[2]] - 0.05, yt = -lo + 0.45, yb = -hi - 0.45,
     top = If[lo > 1, Table[4 (j - 1) + lcColumn[lo - 1] + {-p, p}, {j, d}], {}], bot = If[hi < 2 d + 3, Table[4 (j - 1) + lcColumn[hi] + {-p, p}, {j, d}], {}]},
    Show[lcCones[{hJ, d}], Graphics[{ColorData[97][3], AbsoluteThickness[2.5], Dashed, Line[Join[
         {{xl, yt}}, Catenate[{{#[[1]], yt}, {#[[1]], -lo - 0.375 - p}, {#[[2]], -lo - 0.375 - p}, {#[[2]], yt}} & /@ top], {{xr, yt}, {xr, yb}},
         Catenate[{{#[[2]], yb}, {#[[2]], -hi + 0.375 + p}, {#[[1]], -hi + 0.375 + p}, {#[[1]], yb}} & /@ Reverse[bot]], {{xl, yb}, {xl, yt}}]]}]]]];
lcPurpleZ[hJ_, d_] := With[{m = 2 d + 3, c = d + 2}, With[{st = QuantumCircuitOperator[Join[{StringRepeat["0", m]}, Catenate[Table[With[{k = d - j + 1}, Join[
         QuantumOperator["R"[2 lcJdt, "ZZ"], {#, # + 1}] & /@ Select[Range[m - 1], j > 1 && Max[Abs[# - c], Abs[# + 1 - c]] <= k &],
         QuantumOperator["RX"[2 hJ lcJdt], {#}] & /@ Select[Range[m], Abs[# - c] <= k - 1 &]]], {j, d}]]]][]},
    Re[st["Dagger"][QuantumOperator["Z", {c}][st]]["Scalar"]]]];
lcPurple = AssociationMap[lcPurpleZ @@ # &, Tuples[{{0.5, 1, 2}, Range[2, 5]}]];
lcView3[hJ_, d_, m_] := With[{pts = lcWin[{hJ, d}], exact = lcWin[{hJ, d}][[-1, 2]]}, Grid[{{
     Labeled[lcWindow[hJ, d, m], Column[{Row[{"purple gates reach ", Subscript["Z", "c"], "; grey gates cancel; orange gates only add a phase to ", Ket[{"0\[CenterEllipsis]0"}]}],
        Row[{Style["dashed frame", ColorData[97][3]], ": the ", m, "-site window, stepping around the two-qubit gates it drops"}]}, Center] // Style[#, 12, FontFamily -> "Helvetica"] &, Top],
     ListPlot[{Select[pts, #[[1]] < 2 d - 1 &], Select[pts, #[[1]] >= 2 d - 1 &], {{2 d + 5, lcPurple[{hJ, d}]}}}, PlotStyle -> {Red, ColorData[97][1], RGBColor[0.55, 0.45, 0.85]},
      PlotMarkers -> {{"\[FilledCircle]", 12}, {"\[FilledCircle]", 12}, {"\[FivePointedStar]", 20}},
      PlotLegends -> Placed[PointLegend[{Red, ColorData[97][1], RGBColor[0.55, 0.45, 0.85]}, {"window cuts the light cone", "window contains the light cone", "purple gates only"},
        LegendMarkers -> {{"\[FilledCircle]", 12}, {"\[FilledCircle]", 12}, {"\[FivePointedStar]", 16}}], Below],
      GridLines -> {{2 d - 1}, {exact}}, GridLinesStyle -> Directive[Gray, Dashed], PlotRange -> {{0, 2 d + 6}, All}, PlotRangePadding -> {Automatic, {Scaled[0.05], Scaled[0.16]}},
      Epilog -> {Text[Style["light-cone width 2d \[Minus] 1", 10, Gray], Scaled[{(2 d - 1)/(2 d + 6), 0.5}], {0, -1.4}, {0, 1}], Text[Style["infinite chain", 10, Gray], {0.2, exact}, {-1, 1.4}],
        Inset[Style["\[EmptyCircle]", 30, ColorData[97][3]], {m, Cases[pts, {m, v_} :> v][[1]]}],
        Text[Style[d (2 #1 - 1), 9, GrayLevel[0.35]], {#1, #2}, {0, -2.8}] & @@@ pts, Text[Style[2 d^2 - d, 9, GrayLevel[0.35]], {2 d + 5, lcPurple[{hJ, d}]}, {0, -2.8}]},
      PlotLabel -> Style["numbers: gates in each circuit", 10, Gray],
      Frame -> True, FrameLabel -> {"window size m (sites)", Row[{"\[LeftAngleBracket]", Subscript["Z", "c"], "\[RightAngleBracket]", " at J t = ", Rationalize[lcJdt d]}]},
      FrameTicks -> {{Automatic, None}, {Append[Range[1, 2 d + 3, 2], {2 d + 5, "purple\nonly"}], None}}, ImageSize -> 320]}}, Spacings -> 3, Alignment -> Center]];
lcViews = Association[Flatten[Table[{hJ, d, m} -> lcView3[hJ, d, m], {hJ, {0.5, 1, 2}}, {d, 2, 5}, {m, 1, 2 d + 3, 2}]]];
With[{lcViews = lcViews},
 Manipulate[lcViews[{lcH, lcD, Min[lcM, 2 lcD + 3]}],
 Row[{Control[{{lcH, 1, "h/J"}, {0.5, 1, 2}, SetterBar}], Spacer[20], Control[{{lcD, 2, "Trotter steps d"}, Range[2, 5], SetterBar}], Spacer[20],
   "window m  ", Dynamic[SetterBar[Dynamic[lcM], Range[1, 2 lcD + 3, 2]]]}], {{lcM, 3}, ControlType -> None}]]
```

Every window that contains the light cone, $m\ge2d-1$, gives the same $\langle Z_c\rangle$, so the $2d-1$-site circuit already gives the value for a bulk site of an arbitrarily long chain, while every smaller window cuts the cone and misses it. Stepping through $d$, the magnetization falls faster as $h/J$ grows, and at $h/J=2$ it changes sign within the times shown.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A $1/\sqrt{x}$ phase circuit: its error, rounding included, proved for every input

A circuit that applies $|x\rangle\mapsto e^{2\pi i/\sqrt x}|x\rangle$ needs $1/\sqrt x$ in fixed-point arithmetic. Writing $x=4^km$ with $1\le m<4$ reduces it to $1/\sqrt m$, which the circuit replaces by a table: on each of $3\cdot2^r$ intervals $m=a+t$, $0\le t\le h=2^{-r}$, a quadratic with coefficients $c_0,c_1,c_2$ in $2^{-w}\mathbb Z$, evaluated as $(c_2t+c_1)t+c_0$ in $w$-bit registers that round after each step. Does the circuit's result stay within a relative error $\eta$ of $1/\sqrt m$, with $\eta$ proposed from finitely many samples, for every $m$ and every rounding?

For $r=1$ and $w=8$ the six quadratics lie on $1/\sqrt m$ (top left); on the first interval the relative error is known only at the samples that propose $\eta$ (bottom left):

```wl
intervalsCert[r_, w_] := With[{h = 2^-r, e = 2^-(w + 1), nodes = 2^-r (1 - Cos[Pi (2 Range[3] - 1)/6])/2, ts = Subdivide[0., 2.^-r, 200]}, Table[
     With[{cs = Round[2^w CoefficientList[InterpolatingPolynomial[Transpose[{nodes, N[1/Sqrt[a + nodes], 30]}], \[FormalT]], \[FormalT]]]/2^w},
      With[{horner = Function[{t, d1, d2}, (cs[[3]] t + cs[[2]] + d1) t + cs[[1]] + d2], corners = Tuples[{-e, e}, 2]},
       With[{eta = Ceiling[1.02 2^40 Max[Abs[Sqrt[a + ts] # - 1] & /@ (horner[ts, ##] & @@@ corners)]]/2^40},
        <|"a" -> a, "h" -> h, "poly" -> cs . \[FormalT]^Range[0, 2], "horner" -> horner, "corners" -> corners, "eta" -> eta,
         "proof" -> Resolve[ForAll[{\[FormalT], \[FormalD], \[FormalE]}, 0 <= \[FormalT] <= h && Abs[\[FormalD]] <= e && Abs[\[FormalE]] <= e,
            horner[\[FormalT], \[FormalD], \[FormalE]] > 0 && (1 - eta)^2 <= (a + \[FormalT]) horner[\[FormalT], \[FormalD], \[FormalE]]^2 <= (1 + eta)^2], Reals]|>]]], {a, 1, 4 - h, h}]];
worstCase[iv_, m_] := Max[Function[c, Abs[Sqrt[m] iv["horner"][m - iv["a"], c[[1]], c[[2]]] - 1]] /@ iv["corners"]];
lookupLabel = Style[Column[{Row[{"j \[RightArrow] ", Subscript["c", 0], ", ", Subscript["c", 1], ", ", Subscript["c", 2]}], Row[{"j: 1 of 3\[CenterDot]", Superscript[2, "r"], " pieces"}]}, Alignment -> Center], Black];
lookupStripes = LinearGradientFilling[Flatten[Table[{(j - 1)/6 -> Lighter[ColorData[97][j], 0.25], j/6 -> Lighter[ColorData[97][j], 0.25]}, {j, 6}]], Bottom];
blk[label_, wires_] := QuantumOperator[IdentityMatrix[2^Length[wires]], wires, "Label" -> label];
phaseSteps = QuantumCircuitOperator[{
    blk[Row[{"x = ", Superscript[4, "k"], " m"}], {1, 2}],
    blk[lookupLabel, {1, 3}],
    blk[Row[{"y = ", Subscript["c", 0], " + ", Subscript["c", 1], " t + ", Subscript["c", 2], " ", Superscript["t", 2]}], {1, 3, 4}],
    blk[Row[{Superscript[2, "\[Minus]k"], " y"}], {2, 4}],
    blk[Superscript["e", Row[{"2\[Pi]i y"}]], {4}],
    blk["undo", {1, 2, 3, 4}]}];
phaseKey = Style[Grid[{
     {Row[{"x = ", Superscript[4, "k"], " m"}], Row[{"normalization: leading‐one detection gives k, a barrel shift gives m in [1, 4)"}]},
     {Row[{"j \[RightArrow] ", Subscript["c", 0], ", ", Subscript["c", 1], ", ", Subscript["c", 2]}], Row[{"QROM table lookup of piece j's coefficients"}]},
     {Row[{"y = ", Subscript["c", 0], " + ", Subscript["c", 1], " t + ", Subscript["c", 2], " ", Superscript["t", 2]}], Row[{"Horner evaluation with fixed‐point multiply‐adds, rounded to w bits"}]},
     {Row[{Superscript[2, "\[Minus]k"], " y"}], Row[{"barrel shift by k, restoring the scale: y \[TildeTilde] 1/", Sqrt["x"]}]},
     {Superscript["e", Row[{"2\[Pi]i y"}]], Row[{"phase gradient: one controlled phase rotation per bit of y"}]},
     {"undo", Row[{"uncomputation of every register except the phase"}]}}, Alignment -> {Left, Top}, Spacings -> {1.5, 0.6}, ItemSize -> {{Automatic, 30}}], 12, FontFamily -> "Helvetica", LineIndent -> 0];
phaseCircuit = Column[{phaseSteps["Diagram", "ShowConnectors" -> True, "GateBackgroundStyle" -> {lookupLabel -> lookupStripes}, "WireLabels" -> {"|x\[RightAngleBracket]", "k", Column[{Row[{Subscript["c", 0], ", ", Subscript["c", 1], ", ", Subscript["c", 2]}], "w bits each"}, Alignment -> Right], "y"}, ImageSize -> 440], phaseKey}, Spacings -> 1.5, Alignment -> Center];
problemIvs = intervalsCert[1, 8];
problemPanels = With[{iv = First[problemIvs], ms = 1 + Subdivide[0, 1/2, 20]}, {
    Show[Plot[1/Sqrt[m], {m, 1, 4}, PlotStyle -> Directive[GrayLevel[0.6], AbsoluteThickness[7]]],
     Plot[Evaluate[Piecewise[{{#poly /. \[FormalT] -> m - #a, #a <= m < #a + #h}}, Indeterminate] & /@ problemIvs], {m, 1, 4}, PlotStyle -> (Directive[ColorData[97][#], AbsoluteThickness[2]] & /@ Range[6]), Exclusions -> None],
     Graphics[{GrayLevel[0.5], Dashed, Line[{{#, 0.45}, {#, 1.02}}] & /@ Range[1, 4, 1/2]}],
     Prolog -> {Opacity[0.18], ColorData[97][1], Rectangle[{1, 0.45}, {1.5, 1.02}]},
     Epilog -> {Text[Style[Column[{"zoomed", "below"}, Alignment -> Center], 11, FontFamily -> "Helvetica"], {1.25, 0.52}],
       Table[Text[Style[Row[{"j = ", j - 1}], 11, Bold, ColorData[97][j], FontFamily -> "Helvetica"], {1 + (j - 1/2)/2, 0.98}], {j, 6}]},
     PlotRange -> {{1, 4}, {0.45, 1.02}}, Frame -> True, Axes -> False, FrameLabel -> {"m", None}, ImageSize -> 340,
     PlotLabel -> Style[Column[{Row[{"gray: exact ", 1/Sqrt["m"]}], "colors: the circuit's six quadratics, one per piece"}, Alignment -> Center], 12, FontFamily -> "Helvetica"]],
    Show[ListPlot[Transpose[{ms, N[worstCase[iv, #] & /@ ms]}], PlotStyle -> Directive[ColorData[97][1], AbsolutePointSize[6]]],
     Plot[N[iv["eta"]], {m, 1, 3/2}, PlotStyle -> Directive[ColorData[97][2], Dashed]],
     PlotRange -> {{1, 3/2}, {0, 1.08 N[iv["eta"]]}}, PlotRangePadding -> Scaled[0.03], Prolog -> {Opacity[0.18], ColorData[97][1], Rectangle[Scaled[{0, 0}], Scaled[{1, 1}]]}, Frame -> True, Axes -> False, FrameLabel -> {"m", "|relative error|"}, ImageSize -> 320,
     Epilog -> Inset[Style[Grid[{{Graphics[{ColorData[97][1], Disk[]}, ImageSize -> 7], "samples"},
          {Graphics[{ColorData[97][2], AbsoluteDashing[{4, 3}], AbsoluteThickness[1.5], Line[{{0, 0}, {1, 0}}]}, ImageSize -> {26, 6}, AspectRatio -> 6/26, PlotRange -> {{0, 1}, {-0.1, 0.1}}], Row[{"bound ", "\[Eta]", " suggested by samples"}]}}, Alignment -> {{Center, Left}, Center}, Spacings -> {0.6, 0.4}], 11, FontFamily -> "Helvetica"], Scaled[{0.97, 0.04}], {Right, Bottom}]]}];
Row[{Column[problemPanels, Spacings -> 2, Alignment -> Center], phaseCircuit}, Spacer[30], Alignment -> Center]
```

On each interval the claim is a polynomial inequality in three real variables,

$$
\forall\,t\in{}[0,h],\ |\delta_1|,|\delta_2|\le2^{-w-1}:\quad P>0\ \wedge\ (1-\eta)^2\le(a+t)\,P^2\le(1+\eta)^2,\qquad P=(c_2t+c_1+\delta_1)\,t+c_0+\delta_2,
$$

for $0\le\eta<1$, decided exactly over the reals, with rational $a$, $h$, $\eta$ and coefficients and no rounding in the proof itself. The quadratic interpolates $1/\sqrt{a+t}$ at Chebyshev nodes before its coefficients are rounded, and $\eta$ is proposed from 201 samples at the four extreme roundings, with a 2% margin.

For each $r$ and $w$ the figure plots the worst case of $|\sqrt m\,P-1|$ over the rounding (solid) against $m$, with each interval's $\eta$ (dashed), and reports whether the claim is proved on all intervals:

```wl
certs = AssociationMap[intervalsCert @@ # &, Tuples[{{1, 2}, {8, 12, 16}}]];
panelCert[{r_, w_}] := With[{ivs = certs[{r, w}]},
   Show[Plot[Evaluate[Piecewise[Function[iv, {worstCase[iv, m], iv["a"] <= m < iv["a"] + iv["h"]}] /@ ivs]], {m, 1, 4}, PlotRange -> All, PlotPoints -> 100, Exclusions -> None],
    Plot[Evaluate[Piecewise[N[{#eta, #a <= m < #a + #h}] & /@ ivs]], {m, 1, 4}, PlotStyle -> Directive[ColorData[97][2], Dashed], PlotRange -> All, Exclusions -> None],
    PlotRange -> {0, All}, Axes -> False, Frame -> True, FrameLabel -> {"m", "|relative error|"}, ImageSize -> 280,
    PlotLabel -> Style[Column[{Row[{"r = ", r, ":  ", Length[ivs], " intervals of width ", 2^-r}], Row[{"w = ", w, ":  every register in steps of ", Superscript[2, -w]}]}, Alignment -> Center], 13, FontFamily -> "Helvetica"]]];
Legended[Grid[Map[panelCert, Outer[List, {1, 2}, {8, 12, 16}], {2}]], legendBelow[LineLegend[{ColorData[97][1], Directive[ColorData[97][2], Dashed]},
  {"worst case over the circuit's rounding", If[And @@ Flatten[Lookup[#, "proof"] & /@ Values[certs]], "bound \[Eta] on each interval, proved by Resolve on all of them", "bound \[Eta] on each interval: NOT proved on every interval"]}]]]
```

The bound is proved on every interval in every panel, so the circuit's error stays under the dashed line for every real $m$ and every rounding, not only at the samples that proposed $\eta$. With 8-bit registers the rounding of the coefficients and of each step sets the error, which jumps from interval to interval and barely shrinks when the intervals halve. At 12 bits rounding and the interpolation's error are comparable on the wider intervals, while on the narrower ones rounding takes over again. At 16 bits the interpolation's ripple sets the error where it is largest, on the wider intervals and at small $m$; where the ripple has shrunk, on the narrower intervals and at larger $m$, rounding sets the error again. The model takes $t$ as exact and the shift by $2^{-k}$ as error-free.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Error correction

The examples in this section are error correction: entanglement read from a thousand-qubit stabilizer state; codes that fail through bit flips on stored data, through Pauli errors on stored data and through faults in the circuit that measures them; continuous monitoring and feedback; and the $T$-gate cost of fault tolerance. QF runs the codes as stabilizer tableaux and returns most entropies and error rates as exact numbers and polynomials.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A thousand-qubit GHZ state: entanglement entropy as a rank over F₂

Take $n$ qubits in the GHZ state $(|0\cdots0\rangle+|1\cdots1\rangle)/\sqrt2$, and for comparison in the one-dimensional cluster state, made by putting every qubit in $|+\rangle$ and applying $CZ$ to each neighboring pair. How much entanglement separates a set $A$ of $k$ qubits from the rest, when $A$ is a block at the end of the chain and when $A$ is scattered at random, and does the answer hold at $n=1000$?

A stabilizer state is fixed by $n$ commuting Pauli generators: $X^{\otimes n}$ and $Z_iZ_{i+1}$ for GHZ, $Z_{i-1}X_iZ_{i+1}$ for the cluster. Write each generator as a binary vector $(x|z)\in\mathbb F_2^{2n}$, stack them as the rows of an $n\times2n$ matrix $G$, and let $G_A$ keep only the $x$ and $z$ columns of the qubits in $A$. Then the entanglement entropy, in bits, is (Fattal et al. 2004)

$$
S(A)=\operatorname{rank}_{\mathbb F_2}G_A-|A|.
$$

The same rank has a picture. Applying $H$ to qubits $2,\dots,n$ turns the GHZ state into $(|0\rangle|{+}\cdots{+}\rangle+|1\rangle|{-}\cdots{-}\rangle)/\sqrt2$, the graph state of a star centered on qubit 1, and leaves every $S(A)$ unchanged; the cluster state is the graph state of the chain. In a graph state each edge is a $CZ$, and $S(A)$ is the rank over $\mathbb F_2$ of the block of the adjacency matrix that links $A$ to the rest, whose ones are the edges that cross the cut (Hein, Eisert and Briegel 2004). For $n=20$, with $A$ shaded, the crossing edges in red and $S(A)$ read from each tableau:

```wl
ghzGN = 20;
ghzGStar = UndirectedEdge[1, #] & /@ Range[2, ghzGN];
ghzGChain = UndirectedEdge @@@ Partition[Range[ghzGN], 2, 1];
ghzGStates = {PauliStabilizer[ghzGN][Prepend[Table["CNOT" -> {q, q + 1}, {q, ghzGN - 1}], "H"]],
   PauliStabilizer[ghzGN][Join[Thread["H" -> Range[ghzGN]], Table["CZ" -> {q, q + 1}, {q, ghzGN - 1}]]]};
ghzGLayouts = {{Prepend[CirclePoints[{1, Pi/2}, ghzGN - 1], {0, 0}], 0.065, 230}, {Table[{q, 0}, {q, ghzGN}], 0.28, 560}};
ghzGPanel[name_, edges_, {pos_, r_, size_}, ps_, a_] := With[{cut = Select[edges, Xor @@ (MemberQ[a, #] & /@ List @@ #) &]},
   Graphics[{{Opacity[0.25], ColorData[97][2], Disk[pos[[#]], 1.7 r] & /@ a},
     {AbsoluteThickness[1.2], GrayLevel[0.8], Line[pos[[List @@ #]]] & /@ Complement[edges, cut]},
     {AbsoluteThickness[2.2], ColorData[97][4], Line[pos[[List @@ #]]] & /@ cut},
     {EdgeForm[GrayLevel[0.3]], If[MemberQ[a, #], ColorData[97][2], White], Disk[pos[[#]], r]} & /@ Range[ghzGN]},
    ImageSize -> size, PlotLabel -> Column[{name, Row[{Length[cut], If[Length[cut] == 1, " edge", " edges"], " cut,   S(A) = ", ps["Entropy", a]}]}, Alignment -> Center]]];
ghzGSet[k_, "block"] := Range[k];
ghzGSet[k_, "random"] := BlockRandom[SeedRandom[k]; Sort[RandomSample[Range[ghzGN], k]]];
ghzGViews = Association @ Flatten @ Table[{k, cut} -> Grid[{MapThread[ghzGPanel[#1, #2, #3, #4, ghzGSet[k, cut]] &,
       {{"GHZ (star)", "cluster (chain)"}, {ghzGStar, ghzGChain}, ghzGLayouts, ghzGStates}]}, Spacings -> 4, Alignment -> Center], {k, ghzGN - 1}, {cut, {"block", "random"}}];
With[{ghzGViews = ghzGViews, ghzGN = ghzGN}, Manipulate[
  Column[{ghzGViews[{ghzK, ghzCut}],
    Style[Row[{"A = ", If[ghzCut == "block", Row[{"qubits 1 to ", ghzK}], Row[{ghzK, " qubits chosen at random"}]], " (shaded);  S(A) = entanglement between A and the other ", ghzGN - ghzK, " qubits"}], 12, FontFamily -> "Helvetica"]}, Center, Spacings -> 1.5],
  {{ghzK, 6, "k = number of qubits in A"}, 1, 19, 1, Appearance -> "Labeled"},
  {{ghzCut, "random", "A is"}, {"block" -> "a block 1 to k", "random" -> "k qubits at random"}}]]
```

In the star every crossing edge ends on qubit 1, so the crossing edges have rank 1 however many there are. In the chain a block is crossed by one edge, while a scattered $A$ is crossed by many separate edges, and most of them add to $S(A)$.

Both states at $n=1000$ are held as tableaux, and the rank is evaluated for $k=50,100,\dots,950$, with one random set $A$ per $k$; no $2^{1000}$-entry vector is ever formed:

```wl
ghzChain = PauliStabilizer[1000][Prepend[Table["CNOT" -> {q, q + 1}, {q, 999}], "H"]];
clusterChain = PauliStabilizer[1000][Join[Thread["H" -> Range[1000]], Table["CZ" -> {q, q + 1}, {q, 999}]]];
entropyPanel[label_, ps_] := ListLogPlot[Table[{k, ps["Entropy", sel[k]]}, {sel, {Range, RandomSample[Range[1000], #] &}}, {k, 50, 950, 50}],
   Joined -> True, PlotMarkers -> {{"\[EmptyCircle]", 16}, {"\[EmptyUpTriangle]", 8}}, PlotStyle -> {Automatic, Dashed}, PlotRange -> {0.7, 500}, GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.8], Dashed],
   PlotLegends -> {"block {1, \[Ellipsis], k}", "k qubits at random"}, Frame -> True, FrameLabel -> {"k", "S(A)  (bits)"}, PlotLabel -> label, ImageSize -> 340];
entropyPanels = BlockRandom[SeedRandom[1]; {entropyPanel["GHZ", ghzChain], entropyPanel["cluster", clusterChain]}];
Legended[PlotGrid[{First /@ entropyPanels}, ImageSize -> 680], legendBelow[Last[First[entropyPanels]]]]
```

The GHZ state carries exactly one bit across every cut, a block or a scattered set of any size: its entanglement is a single shared parity. The cluster state also carries one bit across a block boundary, but far more when the same number of qubits is scattered, so its entanglement depends on the shape of the cut, not only on its size.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Repetition code: the error rate of the decoded bit is a polynomial in the flip probability

One bit is stored in an odd number $n$ of qubits as $|0\cdots0\rangle$. Each qubit is flipped independently with probability $p$, and the bit is decoded by majority vote. For which $p$ does a longer code lower the probability $P_{\mathrm{bit}}$ that the decoded bit is wrong?

Each qubit passes through the bit-flip channel

$$
\mathcal E_p(\rho)=(1-p)\,\rho+p\,X\rho X,
$$

and majority vote fails when more than $n/2$ qubits flipped, so

$$
P_{\mathrm{bit}}(n)=\sum_{k>n/2}\binom nk p^k(1-p)^{n-k}=I_p\left(\frac{n+1}2,\frac{n+1}2\right).
$$

Whether majority vote fails depends only on how many qubits flipped, not on which ones. So for each $k=0,\dots,n$ the stabilizer state $|0\cdots0\rangle$ gets its first $k$ qubits flipped by $X$ gates, qubit $j$ reads 1 when its generator $Z_j$ has sign $-1$, and a tableau that decodes wrong contributes $\binom nk p^k(1-p)^{n-k}$, the weight of all patterns with $k$ flips. This gives $P_{\mathrm{bit}}(p)$ exactly from $n+1$ tableaux instead of $2^n$ branches, which reaches $n=101$, and it equals $I_p$ for every $n$ shown:

```wl
repFlipped[n_, k_] := QuantumCircuitOperator[Table["X" -> q, {q, k}]][PauliStabilizer[n]];
repRate[n_] := Total[Table[Boole[Count[repFlipped[n, k]["StabilizerSigns"], -1] > n/2] Binomial[n, k] p^k (1 - p)^(n - k), {k, 0, n}]];
repNs = {1, 3, 5, 7, 21, 101}; repRates = repRate /@ repNs;
repPanels = {Plot[Evaluate[repRates], {p, 0, 1}, WorkingPrecision -> 20, PlotLegends -> (Row[{"n = ", #}] & /@ repNs), Frame -> True, FrameLabel -> {"flip probability p of each qubit", Row[{Subscript["P", "bit"], ":  probability the majority vote is wrong"}]},
    Prolog -> {{Opacity[0.18], RGBColor[0.3, 0.75, 0.45], Rectangle[{0, 0}, {1/2, 1}]}, {Opacity[0.18], RGBColor[0.9, 0.4, 0.35], Rectangle[{1/2, 0}, {1, 1}]},
      Text[Style["longer code helps", 11, FontFamily -> "Helvetica"], {1/4, 0.94}], Text[Style["longer code hurts", 11, FontFamily -> "Helvetica"], {3/4, 0.06}]},
    PlotLabel -> Row[{Subscript["P", "bit"], "(n+2) < ", Subscript["P", "bit"], "(n)  \[DoubleLongLeftRightArrow]  ",
       First[Union[MapThread[Reduce[#2 < #1 && 0 < p < 1, p] &, {Most[Take[repRates, 4]], Rest[Take[repRates, 4]]}]]]}], ImageSize -> 320],
   LogLogPlot[Evaluate[repRates], {p, 10^-3, 1/2},
    WorkingPrecision -> 20, PlotRange -> {{10^-3, 1/2}, {10^-16, 1}},
    Prolog -> {{Opacity[0.18], RGBColor[0.3, 0.75, 0.45], Rectangle[Scaled[{0, 0}], Scaled[{1, 1}]]}, Text[Style["longer code helps", 11, FontFamily -> "Helvetica"], Scaled[{0.03, 0.06}], {-1, 0}]}, Frame -> True, FrameLabel -> {"flip probability p of each qubit", Row[{Subscript["P", "bit"], ":  probability the majority vote is wrong"}]},
    PlotLabel -> Row[{Subscript["P", "bit"], " \[TildeTilde] ", TraditionalForm[HoldForm[Binomial[n, (n + 1)/2] p^((n + 1)/2)]], "  at small p"}], ImageSize -> 320]};
Legended[Row[{First[First[repPanels]], Last[repPanels]}, Spacer[10]], legendBelow[Last[First[repPanels]]]]
```

For $n=3$ the rate has no linear term, so at small $p$ the code fails far less often than a bare qubit; on the log-log axes each code becomes a straight line at small $p$, with slope equal to the smallest number of flips that defeats the majority. All curves meet at $p=1/2$ and tend to a step there as $n$ grows, and the exact threshold is the same for each pair of neighboring codes drawn: below one half the longer code is better, above it worse. For every odd $n=2m+1$ the difference factors as $P_{\mathrm{bit}}(n+2)-P_{\mathrm{bit}}(n)=\binom{2m+1}{m}\,{}[p(1-p)]^{m+1}(2p-1)$, which has the sign of $2p-1$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Maximum-likelihood decoding: exact logical error rates of three codes

A code stores one qubit in $n$ physical qubits. Each qubit independently suffers $X$, $Y$ or $Z$ with probabilities $p_X$, $p_Y$, $p_Z$, the checks are read without error, and a decoder turns the syndrome into a correction. How much lower is the probability $P_L$ of a logical error when the decoder picks the most probable class of errors instead of the lightest error?

Each qubit passes through the Pauli channel

$$
\mathcal E(\rho)=(1-p_X-p_Y-p_Z)\,\rho+p_X\,X\rho X+p_Y\,Y\rho Y+p_Z\,Z\rho Z,
$$

depolarizing when $p_X=p_Y=p_Z=p/3$ and $Z$-biased here as $p_X=p_Y=p/10$, $p_Z=p$. An error $E$ on the $n$ qubits has a syndrome, the checks it anticommutes with, and a class, the logical operators it anticommutes with. The minimum-weight decoder answers each syndrome with the class of the lightest error that produces it, taking the first in a fixed order when several are equally light; the maximum-likelihood decoder answers with the class of largest total probability. A round fails when the answer is not the class of the error that occurred:

$$
P_L=\sum_{E}\Pr(E)\,\left[\operatorname{class}(E)\neq\operatorname{guess}(\operatorname{syndrome}(E))\right].
$$

Each code's check matrix, one row per stabilizer with its $X$ part in blue and its $Z$ part in orange, shows where degeneracy comes from: the Shor code's six $Z$ checks each touch only two neighboring qubits, so many light errors differ by a check and share a class. Summing over all $4^n$ errors gives $P_L$ as an exact polynomial in $p$ for each decoder. The maximum-likelihood decoder picks, for each syndrome, the class that is most probable at small $p$, and keeps that choice for every $p$; for $p$ up to $0.2$, the range of the top row, that choice stays the most probable:

```wl
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/QEC/QECCore/QECCore.wl"]; Needs["Wolfram`QuantumFramework`QEC`"];
mlDecCases = Append[MapThread[{#1, #2 <> ", depolarizing", QECNoiseModel["Depolarizing", p]} &, {{"5QubitCode", "ShorCode", "SteaneCode"}, {"five-qubit code", "Shor code", "Steane code"}}],
   {"SteaneCode", "Steane code, Z-biased", QECNoiseModel[<|"X" -> p/10, "Y" -> p/10, "Z" -> p|>]}];
mlDecRates[code_, noise_] := Expand[QECLogicalErrorRate[QECCode[code], noise, "Decoder" -> #, "DecoderReach" -> 3]] & /@ {"MaximumLikelihood", "MinimumWeight"};
mlDecPL = Association[(#[[2]] -> mlDecRates[#[[1]], #[[3]]]) & /@ mlDecCases];
mlCheckPlot[code_, title_] := With[{c = QECCode[code]}, With[{h = c["CheckMatrix"], n = c["Qubits"]},
   Labeled[ArrayPlot[MapThread[Join[#1, 2 #2] &, {h[[All, ;; n]], h[[All, n + 1 ;;]]}], ColorRules -> {0 -> White, 1 -> ColorData[97][1], 2 -> ColorData[97][2]},
     Mesh -> All, MeshStyle -> GrayLevel[0.85], Frame -> True, FrameTicks -> None, Epilog -> {GrayLevel[0.3], AbsoluteThickness[2], Line[{{n, 0}, {n, Length[h]}}]}, ImageSize -> 20 n],
    Column[{Style[title, 12, FontFamily -> "Helvetica"], Style[Row[{"[[", Row[c["Parameters"], ","], "]]"}], 11, Gray, FontFamily -> "Helvetica"]}, Center], Top]]];
mlChecks = Labeled[Row[MapThread[mlCheckPlot, {{"5QubitCode", "SteaneCode", "ShorCode"}, {"five‐qubit code", "Steane code", "Shor code"}}], Spacer[30]],
   Style[Column[{"check matrices:", "one row per stabilizer,", "one column per qubit", Row[{Style["X part", ColorData[97][1]], " | ", Style["Z part", ColorData[97][2]]}]}, Alignment -> Left, Spacings -> 0.3], 11, FontFamily -> "Helvetica"], Right];
Column[{mlChecks, Grid[Prepend[KeyValueMap[{#1, TraditionalForm[#2[[1]]], TraditionalForm[#2[[2]]]} &, mlDecPL], Style[#, Bold] & /@ {Row[{Subscript["P", "L"], ": total probability\nof a harmful correction"}], "maximum‐likelihood decoder", "minimum‐weight decoder"}],
 Alignment -> Left, Dividers -> {None, {2 -> True}}, Spacings -> {2, 1.2}, BaseStyle -> {12, FontFamily -> "Helvetica"}, FrameMargins -> 8]}, Spacings -> 2, Alignment -> Center]
```

The top row plots $P_L/p^2$, whose value at $p=0$ is the coefficient of $p^2$, the term that sets $P_L$ at small $p$; the bottom row sets $P_L$ against $p_X+p_Y+p_Z$, the probability that an unencoded qubit suffers an error, and marks the pseudothreshold, the $p$ below which the encoded qubit fails less often:

```wl
mlDecPanel[{code_, title_, noise_}] := With[{r = mlDecPL[title]},
   With[{label = Style[title, 11, FontFamily -> "Helvetica"],
     order2 = Style[Row[{"order ", Superscript["p", 2], ":  ", TraditionalForm[Coefficient[r[[1]], p, 2] p^2], "  against  ", TraditionalForm[Coefficient[r[[2]], p, 2] p^2]}], 11, FontFamily -> "Helvetica"]},
    Plot[Evaluate[Expand[r/p^2]], {p, 0, 0.2}, PlotStyle -> {Automatic, Dashed}, PlotRange -> All, PlotLegends -> {"maximum‐likelihood decoder:\nundoes the most probable class of errors", "minimum‐weight decoder:\nundoes the lightest error"}, Axes -> False, GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.8], Dashed],
     Frame -> True, FrameLabel -> {"p", Row[{Subscript["P", "L"], " / ", Superscript["p", 2]}]}, AspectRatio -> 1/2, ImageSize -> 290, PlotLabel -> label,
     Epilog -> Inset[order2, Scaled[{0.97, 0.94}], {Right, Top}]]]];
mlDecPanels = mlDecPanel /@ mlDecCases;
mlPTPanel[{code_, title_, noise_}] := With[{r = mlDecPL[title], bare = noise["MeanWeight"]},
  With[{pts = First[SortBy[p /. Solve[# == bare && 0 < p < 1, p, Reals], N]] & /@ r},
   LogLogPlot[Evaluate[Append[r, bare]], {p, 10^-3, 1/2}, PlotStyle -> {ColorData[97][1], Directive[ColorData[97][2], Dashed], Directive[ColorData[97][3], Dashed]},
    Frame -> True, FrameLabel -> {"p", Subscript["P", "L"]}, AspectRatio -> 1/2, GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.5, 0.35], Thin],
    PlotLabel -> Style[Row[{"pseudothreshold  ", Row[NumberForm[N[#], 2] & /@ pts, "  against  "]}], 11, FontFamily -> "Helvetica"], FrameTicks -> {{Automatic, None}, {{0.001, 0.01, 0.1}, None}}, ImageSize -> 290,
    Prolog -> With[{a = Log[10.^-4], b = Log[2.], ya = Log[N[bare /. p -> 10^-4]], yb = Log[N[bare /. p -> 2]], lo = Log[10.^-8], hi = Log[3.]},
      {{Opacity[0.12], ColorData[97][3], Polygon[{{a, ya}, {b, yb}, {b, lo}, {a, lo}}]}, {Opacity[0.12], ColorData[97][4], Polygon[{{a, ya}, {b, yb}, {b, hi}, {a, hi}}]},
       Text[Style["error correction helps", 10, ColorData[97][3], FontFamily -> "Helvetica"], Scaled[{0.72, 0.1}]], Text[Style["error correction hurts", 10, ColorData[97][4], FontFamily -> "Helvetica"], Scaled[{0.28, 0.9}]]}],
    Epilog -> {PointSize[0.03], MapThread[{#2, Point[Log[N[{#1, bare /. p -> #1}]]]} &, {pts, {ColorData[97][1], ColorData[97][2]}}]}]]];
mlPTPanels = mlPTPanel /@ mlDecCases;
Legended[Column[{PlotGrid[{First /@ mlDecPanels}, ImageSize -> 1200], PlotGrid[{mlPTPanels}, ImageSize -> 1200]}, Spacings -> 1.5],
 Placed[LineLegend[{ColorData[97][1], Directive[ColorData[97][2], Dashed], Directive[ColorData[97][3], Dashed]},
   {"maximum\[Hyphen]likelihood decoder:\nchoices fixed at small p, exact up to p = 0.2", "minimum\[Hyphen]weight decoder:\nundoes the lightest error", Row[{"unencoded qubit: error probability ", Subscript["p", "X"], " + ", Subscript["p", "Y"], " + ", Subscript["p", "Z"]}]}, LegendLayout -> {"Row", 1}], Below]]
```

Under depolarizing noise the two decoders give the same polynomial for each of the three codes, the same in both columns of the table. $X$, $Y$ and $Z$ are equally likely, so the lightest errors of a syndrome are also its most probable ones, and for these three codes the class the lightest error points to also carries the most probability; for the five-qubit code this is plain, since its fifteen nonzero syndromes come from the fifteen one-qubit errors. Under $Z$-biased noise the two part. A Steane syndrome left by two errors can also be left by two other errors of the same weight, one pair an $X$ and a $Z$, the other an $X$ and a $Y$, ten times less likely, and in different classes; counting weight alone, the minimum-weight decoder sometimes undoes the unlikely pair, while maximum likelihood weighs the classes and lowers the $p^2$ coefficient. Over the plotted range the maximum-likelihood curve never lies above the minimum-weight one.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A surface code measured by noisy gates: the exact logical error rate

Nine data qubits on a $3\times3$ grid carry one logical qubit in the distance-3 surface code. Its eight checks, four of $X$ type and four of $Z$ type, each act on two or four neighboring qubits. One round of check measurement is made by a circuit in which every CNOT fails with probability $p_C$ and every readout with probability $p_M$. How likely is a decoder to lose the logical qubit, and is every single fault still corrected?

Each check has its own ancilla (qubits 10 to 17), starting in $|0\rangle$. For an $X$ check its data qubits are rotated by $H$, CNOTs copy their parity into the ancilla in a fixed order ($X$ checks row by row, $Z$ checks column by column), the rotations are undone and the ancilla is read:

```wl
surfChecks = {"IXXIIIIII", "XXIXXIIII", "IIIIXXIXX", "IIIIIIXXI", "ZIIZIIIII", "IZZIZZIII", "IIIZZIZZI", "IIIIIZIIZ"};
surfOrder = {{2, 3}, {1, 2, 4, 5}, {5, 6, 8, 9}, {7, 8}, {1, 4}, {2, 5, 3, 6}, {4, 7, 5, 8}, {6, 9}};
surfXCheckQ = StringContainsQ["X"];
surfCircuit = Catenate[MapThread[Function[{qs, a, check}, With[{rot = If[surfXCheckQ[check], Thread["H" -> qs], {}]},
      Join[rot, Thread["CNOT" -> Thread[{qs, a}]], rot, {"M" -> a}]]], {surfOrder, Range[10, 17], surfChecks}]];
surfPos[q_] := {Mod[q - 1, 3], -Quotient[q - 1, 3]};
surfOut[qs_] := With[{v = Cross[Subtract @@ (surfPos /@ qs)]}, If[v . (Mean[surfPos /@ qs] - {1, -1}) > 0, v, -v]];
surfFace[qs_] := With[{c = Mean[surfPos /@ qs]}, If[Length[qs] == 4, Polygon[SortBy[surfPos /@ qs, ArcTan @@ (# - c) &]],
    Polygon[Table[c + 0.5 RotationMatrix[a] . surfOut[qs], {a, -Pi/2, Pi/2, Pi/24}]]]];
surfLattice = Graphics[{
    MapThread[Function[{qs, check}, {EdgeForm[GrayLevel[0.4]], FaceForm[Opacity[0.35, If[surfXCheckQ[check], ColorData[97][1], ColorData[97][2]]]], surfFace[qs]}], {surfOrder, surfChecks}],
    {Opacity[0.75], AbsoluteThickness[7], CapForm["Round"], ColorData[97][1], Line[surfPos /@ {3, 6, 9}], ColorData[97][2], Line[surfPos /@ {1, 5, 9}]},
    MapThread[Function[{qs, check}, Text[Style[Row[Subscript[If[surfXCheckQ[check], "X", "Z"], #] & /@ Sort[qs]], 10, Bold, Black], Mean[surfPos /@ qs] + If[Length[qs] == 2, 0.27 surfOut[qs], {0, 0}]]], {surfOrder, surfChecks}],
    Text[Style[Row[{Overscript["X", "_"], " = ", Subscript["X", 3], Subscript["X", 6], Subscript["X", 9]}], 12, Bold, ColorData[97][1]], {2.1, -2.85}],
    Text[Style[Row[{Overscript["Z", "_"], " = ", Subscript["Z", 1], Subscript["Z", 5], Subscript["Z", 9]}], 12, Bold, ColorData[97][2]], {-0.1, -2.85}],
    {EdgeForm[Black], White, Disk[surfPos[#], 0.13], Text[Style[#, 11, Black], surfPos[#]]} & /@ Range[9]},
   PlotLabel -> Style["distance\[Hyphen]3 surface code, [[9, 1, 3]]", 13, FontFamily -> "Helvetica"], ImageSize -> 270, PlotRange -> {{-0.85, 2.85}, {-3.05, 0.85}}];
surfCheckPlot = With[{bits = Function[c, {Boole[MemberQ[{"X", "Y"}, #]] & /@ Characters[c], Boole[MemberQ[{"Z", "Y"}, #]] & /@ Characters[c]}] /@ surfChecks},
   Labeled[ArrayPlot[Join[#[[1]], 2 #[[2]]] & /@ bits, ColorRules -> {0 -> White, 1 -> ColorData[97][1], 2 -> ColorData[97][2]}, Mesh -> All, MeshStyle -> GrayLevel[0.85],
     Frame -> True, FrameTicks -> None, Epilog -> {GrayLevel[0.3], AbsoluteThickness[2], Line[{{9, 0}, {9, 8}}]}, ImageSize -> 230],
    {Style["check matrix, [[9, 1, 3]]", 12, FontFamily -> "Helvetica"],
     Style[Column[{"one row per check, one column per qubit", Row[{Style["X part", ColorData[97][1]], " | ", Style["Z part", ColorData[97][2]]}]}, Alignment -> Center, Spacings -> 0.3], 11, FontFamily -> "Helvetica"]}, {Top, Bottom}]];
Grid[{{surfLattice, surfCheckPlot, Labeled[QuantumCircuitOperator[surfCircuit /. ("M" -> a_) :> {a}]["Diagram", FontSize -> 7, ImageSize -> 430], Style[Column[{"qubits 1\[Dash]9: data;  10\[Dash]17: one ancilla per check", Row[{"faults strike after every CNOT (", Subscript["p", "C"], ") and on every readout (", Subscript["p", "M"], ")"}]}, Alignment -> Center], 11, FontFamily -> "Helvetica"], Bottom]}}, Alignment -> Center, Spacings -> 3]
```

After each CNOT one of the 15 two-qubit Paulis other than the identity strikes, each with probability $p_C/15$, and each ancilla readout is flipped with probability $p_M$; one-qubit gates and waiting qubits do not err. The round is followed by a perfect, noiseless measurement of all eight checks and of the logical operator. A detector is a check's reading in the round, or its change at that final measurement: both are $0$ when nothing goes wrong. A stabilizer simulation of the 17 qubits, started in the code state with $\bar Z$ fixed and once with $\bar X$ fixed, with an $X$ or a $Z$ inserted after a CNOT, gives the detectors that fault fires and whether it flips $\bar X$ or $\bar Z$; every other Pauli follows by adding these bit patterns, since the circuit is Clifford. The faults, counted by how many detectors they fire and by how many of them flip a logical operator:

```wl
surfLogicals = {"IIXIIXIIX", "ZIIIZIIIZ"};
surfPad[s_] := s <> StringRepeat["I", 8];
surfStart[logical_] := PauliStabilizer[Join[surfPad /@ surfChecks, {surfPad[logical]}, Table[StringReplacePart[StringRepeat["I", 17], "Z", {a, a}], {a, 10, 17}]]];
surfBit[ps_, op_] := Enclose[First[Normal[ConfirmBy[ps["M", op], Length[#] == 1 &]]], Abort[] &];
surfStep[{ps_, record_}, "M" -> a_] := With[{r = surfBit[ps, a]}, {Last[r], Append[record, First[r]]}];
surfStep[{ps_, record_}, gate_] := {ps[{gate}], record};
surfRun[logical_, circuit_] := With[{end = Fold[surfStep, {surfStart[logical], {}}, circuit]},
   Join[Last[end], BitXor[First[surfBit[First[end], surfPad[#]]] & /@ surfChecks, Last[end]], {First[surfBit[First[end], surfPad[logical]]]}]];
Enclose[ConfirmAssert[FreeQ[surfRun[#, surfCircuit] & /@ surfLogicals, 1]], Abort[] &];
surfEffect[i_, fault_] := With[{d = surfRun[#, Insert[surfCircuit, fault, i + 1]] & /@ surfLogicals}, Join[Most[Last[d]], Last /@ d]];
surfCNOTs = Flatten[Position[surfCircuit, "CNOT" -> _]];
surfBasis = Association[Flatten[Outer[Function[{i, l, k}, {i, l, k} -> surfEffect[i, l -> surfCircuit[[i, 2, k]]]], surfCNOTs, {"X", "Z"}, {1, 2}], 2]];
surfLetter[i_, k_, "I"] := ConstantArray[0, 18];
surfLetter[i_, k_, l : "X" | "Z"] := surfBasis[{i, l, k}];
surfLetter[i_, k_, "Y"] := BitXor[surfBasis[{i, "X", k}], surfBasis[{i, "Z", k}]];
surfFaults = Join[
   Catenate[Function[i, {i, BitXor[surfLetter[i, 1, #[[1]]], surfLetter[i, 2, #[[2]]]], pC/15} & /@ Rest[Tuples[{"I", "X", "Y", "Z"}, 2]]] /@ surfCNOTs],
   Table[{-j, Join[UnitVector[16, j] + UnitVector[16, j + 8], {0, 0}], pM}, {j, 8}]];
surfFaults = Enclose[ConfirmBy[surfFaults, MatrixQ[#[[All, 2]], MatchQ[0 | 1]] &], Abort[] &];
Labeled[Grid[MapThread[Prepend, {Transpose[KeyValueMap[Prepend[#2, #1] &, KeySort[GroupBy[surfFaults[[All, 2]], Total[#[[;; 16]]] &, {Length[#], Count[#[[All, 17 ;;]], {___, 1, ___}]} &]]]],
   {"detectors the fault fires", "single faults", "of these, faults that flip the logical qubit"}}],
  Alignment -> {{Right, {Center}}, Center}, Dividers -> {{2 -> True}, {2 -> True}}, Spacings -> {1.5, 0.8}, BaseStyle -> {12, FontFamily -> "Helvetica"}],
 Style["every single fault: a Pauli error after one CNOT, or one flipped readout", 12, Bold, FontFamily -> "Helvetica"], Top]
```

Two decoders read the detectors. The first knows the pattern of every single fault and undoes that fault, and counts any other pattern as a failure. The second also knows every pair of faults at different locations: for a pattern no single fault makes, it applies the logical correction whose pairs that make it are together most probable when $p_C=p_M$ (ties go to the first class found), and counts only the remaining patterns as failures. The logical error per round, $\varepsilon_L$, is the chance that this one round, run with noisy gates and readouts, fails: the decoder leaves the logical qubit flipped or refuses the detector pattern. Faults at different locations are independent, and their effects $s\in\{0,1\}^{18}$ (16 detectors, $\bar X$, $\bar Z$) add bit by bit, so in the Walsh-Hadamard transform each location $\ell$ contributes one factor:

$$
\varepsilon_L=2^{-18}\sum_{k\in\{0,1\}^{18}}\hat W(k)\prod_{\ell}\Bigl(1-2\sum_{j\in\ell,\ k\cdot s_j\ \mathrm{odd}}p_j\Bigr),
$$

with $\hat W$ the transform of a decoder's failure indicator and $p_j$ the probability of outcome $j$ at location $\ell$. At a CNOT either none or exactly eight of the fifteen Paulis are odd under $k$, so its factor is $1$ or $1-16p_C/15$; a readout contributes $1$ or $1-2p_M$. For general $p_C$ and $p_M$, plotted at $p_C=p_M=p$:

```wl
surfIdx = Range[0, 2^18 - 1];
surfSig = FromDigits[#, 2] & /@ surfFaults[[All, 2]];
surfSyndrome = BitShiftRight[#, 2] &;
surfOne = Association[0 -> 0, Thread[surfSyndrome[surfSig] -> BitAnd[surfSig, 3]]];
surfPairs = Select[Subsets[Range[Length[surfSig]], {2}], surfFaults[[#[[1]], 1]] =!= surfFaults[[#[[2]], 1]] &];
surfPairSig = BitXor[surfSig[[surfPairs[[All, 1]]]], surfSig[[surfPairs[[All, 2]]]]];
surfPairWeight = Times @@ (surfFaults[[#, 3]] & /@ Transpose[surfPairs]) /. {pC -> 1, pM -> 1};
surfTwo = Join[Map[First[Keys[TakeLargest[#, 1]]] &, GroupBy[Transpose[{surfSyndrome[surfPairSig], BitAnd[surfPairSig, 3], surfPairWeight}], First -> Rest, GroupBy[#, First -> Last, Total] &]], surfOne];
surfBits = IntegerDigits[surfIdx, 2, 18];
surfLocations = Values[PositionIndex[surfFaults[[All, 1]]]];
surfOdd = Total[Mod[surfBits . Transpose[surfFaults[[#, 2]]], 2], {2}] & /@ surfLocations;
surfCNOTQ = Positive[surfFaults[[First[#], 1]]] & /@ surfLocations;
surfCounts = Transpose[{Total[Unitize[Pick[surfOdd, surfCNOTQ]]], Total[Pick[surfOdd, surfCNOTQ, False]]}];
surfRateOf[table_] := Expand[2^-18 Total[KeyValueMap[#2 (1 - 16 pC/15)^#1[[1]] (1 - 2 pM)^#1[[2]] &,
     GroupBy[Transpose[{surfCounts, 2^9 DiscreteHadamardTransform[Unitize[Lookup[table, surfSyndrome[surfIdx], -1] - BitAnd[surfIdx, 3]], Method -> "BitComplement"]}], First -> Last, Total]]]];
surfRates = surfRateOf /@ {surfOne, surfTwo};
surfLead = Expand[Normal[Series[# /. {pC -> t pC, pM -> t pM}, {t, 0, 2}]] /. t -> 1] & /@ surfRates;
surfCross = First[SolveValues[(# /. {pC -> p, pM -> p}) == p && 0 < p < 1/10, p, Reals]] & /@ surfRates;
surfPlot = LogLogPlot[Evaluate[Append[surfRates /. {pC -> p, pM -> p}, p]], {p, 10^-4, 1/2}, PlotStyle -> {Automatic, Automatic, Dashed},
   Prolog -> With[{a = Log[10.^-4], b = Log[0.5], lo = Log[10.^-7], hi = Log[3.]},
     {{Opacity[0.12], ColorData[97][3], Polygon[{{a, a}, {b, b}, {b, lo}, {a, lo}}]}, {Opacity[0.12], ColorData[97][4], Polygon[{{a, a}, {b, b}, {b, hi}, {a, hi}}]},
      Text[Style["error correction helps", 11, ColorData[97][3], FontFamily -> "Helvetica"], {Log[0.05], Log[2. 10^-6]}],
      Text[Style["error correction hurts", 11, ColorData[97][4], FontFamily -> "Helvetica"], {Log[3. 10^-4], Log[0.3]}]}], Frame -> True, GridLines -> Automatic,
   GridLinesStyle -> Directive[GrayLevel[0.8], Dashed], FrameLabel -> {Row[{"error probability p of each CNOT and each readout  (", Subscript["p", "C"], " = ", Subscript["p", "M"], " = p)"}], Row[{Subscript["\[CurlyEpsilon]", "L"], ", logical error per round"}]}, ImageSize -> 480,
   PlotLegends -> {"decoder knows single faults", "decoder knows single faults and pairs", Row[{Subscript["\[CurlyEpsilon]", "L"], " = p"}]},
   PlotLabel -> Column[MapThread[Row[{#1, ":  ", Subscript["\[CurlyEpsilon]", "L"], " = ", TraditionalForm[#2 /. {pC -> Subscript["p", "C"], pM -> Subscript["p", "M"]}], " + \[Ellipsis],   ", Subscript["\[CurlyEpsilon]", "L"], " < p  for  p < ", NumberForm[N[#3], 3]}] &,
      {{"single faults", "and pairs"}, surfLead, surfCross}], Alignment -> Center]];
Legended[First[surfPlot], Placed[Last[surfPlot] /. Placed[\[FormalX]_, ___] :> \[FormalX], Right]]
```

Neither decoder's $\varepsilon_L$ has a linear term: every single fault, in a CNOT or in a readout, is corrected. Most of the first decoder's $p^2$ term comes from patterns it refuses, many of which never touched the logical qubit, and with readout errors alone it fails exactly when two or more of the eight readings are wrong, $1-(1-p_M)^8-8p_M(1-p_M)^7$. The second decoder has no $p_M^2$ term, since two wrong readings are always recognized, and with readout errors alone it first fails at three; its rate is several times lower, and the crossing with $\varepsilon_L=p$, below which a round loses the logical qubit less often than a single CNOT fails, moves to the larger value printed above the plot.

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

The flips are drawn as events at random times, so a true syndrome exists at every moment, and each record is its stabilizer value plus detector noise. The filter starts from the two sectors of $(|000\rangle+i|100\rangle)/\sqrt2$ and sees only the records. The same records also drive the full $8\times8$ equation, integrated with the positivity-preserving step of Rouchon and Ralph at four time steps. For each detector efficiency and flip rate, the left plot shades the true syndrome and draws the probability the filter gives it, with how often the filter's most likely syndrome is the true one; the right plot draws the largest gap between the filter and the sector populations $\mathrm{Tr}[\hat\Pi_s\rho]$ of the full state against the time step. Time is in units of $1/\kappa$:

```wl
qecStabs = Normal[QuantumOperator[#]["Matrix"]] & /@ {"ZZI", "IZZ"};
qecFlips = Normal[QuantumOperator[#]["Matrix"]] & /@ {"XII", "IXI", "IIX"};
qecSigns = {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
qecFlipSigns = {{-1, 1}, {-1, -1}, {1, -1}};
qecTruth[γ_, dt_, n_] := FoldList[#1 (Times @@ Pick[qecFlipSigns, Thread[#2 < γ dt]]) &, RandomChoice[{{1, 1}, {-1, 1}}], RandomReal[1, {n, 3}]];
qecRecord[η_, dt_, truth_] := 2 Sqrt[η] Most[truth] dt + RandomVariate[NormalDistribution[0, Sqrt[dt]], {Length[truth] - 1, 2}];
qecFilter[γ_, η_, dt_, record_] := FoldList[#/Total[#] &[(1/4 + Exp[-4 γ dt] (#1 - 1/4)) Exp[2 Sqrt[η] qecSigns . #2]] &, {0.5, 0.5, 0., 0.}, record];
qecProj = (IdentityMatrix[8] + #[[1]] qecStabs[[1]]) . (IdentityMatrix[8] + #[[2]] qecStabs[[2]])/4 & /@ qecSigns;
qecStep[γ_, η_, dt_][ρ_, dJ_] := With[{σ = Sqrt[η] dJ . qecStabs},
   With[{M = (1 - (2 + 3 γ) dt/2 - η dt) IdentityMatrix[8] + σ + σ . σ/2},
    #/Tr[#] &[M . ρ . M + γ dt Total[# . ρ . # & /@ qecFlips] + (1 - η) dt Total[# . ρ . # & /@ qecStabs]]]];
qecRho0 = Outer[Times, #, Conjugate[#]] &[(UnitVector[8, 1] + I UnitVector[8, 5])/Sqrt[2.]];
qecPopulations[γ_, η_, dt_, record_] := Re[Tr /@ (qecProj . #)] & /@ FoldList[qecStep[γ, η, dt], qecRho0, record];
qecGap[γ_, η_, dt_, record_] := Max[Abs[qecFilter[γ, η, dt, record] - qecPopulations[γ, η, dt, record]]];
qecColors = ColorData[97] /@ Range[4];
qecRow[η_, γ_] := With[{dt = 0.005, n = 8000}, With[{truth = qecTruth[γ, dt, n]}, With[{record = qecRecord[η, dt, truth]},
    With[{filter = qecFilter[γ, η, dt, record], idx = First[FirstPosition[qecSigns, #]] & /@ truth,
      gaps = Table[{m dt, qecGap[γ, η, m dt, Total /@ Partition[record, m]]}, {m, {1, 2, 4, 8}}]},
     With[{bands = {qecColors[[First[#]]], Opacity[0.5], Rectangle[{dt (First[Last[#]] - 1), 0}, {dt Last[Last[#]], 1}]} & /@
         ({#[[1, 1]], #[[All, 2]]} & /@ SplitBy[Transpose[{idx, Range[Length[idx]]}], First]),
       hit = Mean[Boole[MapThread[First[Ordering[#1, -1]] == #2 &, {filter, idx}]]]},
      Grid[{{ListLinePlot[MapThread[Part, {filter, idx}], DataRange -> {0, n dt}, PlotRange -> {{0, n dt}, {0, 1}}, PlotStyle -> Directive[Black, AbsoluteThickness[1]],
          Prolog -> bands, Frame -> True, FrameLabel -> {"κ t", "P(true syndrome)"},
          PlotLabel -> Row[{"most likely syndrome correct ", Round[100 hit], "% of the time"}], ImageSize -> 400, AspectRatio -> 0.5],
         ListLogLogPlot[gaps, Joined -> True, PlotMarkers -> Automatic, Frame -> True, FrameLabel -> {"κ dt", "largest gap"},
          GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.8], Dashed], PlotLabel -> "filter vs 8×8 state, same record", ImageSize -> 300]}},
       Spacings -> 3, Alignment -> Center]]]]]];
qecEtas = {1, 0.2};
qecRates = {0.01, 0.05, 0.2};
qecRows = BlockRandom[SeedRandom[5]; Association[Flatten[Outer[{#1, #2} -> qecRow[#1, #2] &, qecEtas, qecRates]]]];
qecLegend = SwatchLegend[Directive[#, Opacity[0.5]] & /@ qecColors, MapThread[Row[{"(", #1[[1]], ", ", #1[[2]], ")  ", #2}] &, {qecSigns, {"no flip", "qubit 1 flipped", "qubit 2 flipped", "qubit 3 flipped"}}], LegendLabel -> "S₁ = Z₁Z₂,   S₂ = Z₂Z₃", LegendLayout -> {"Row", 1}];
With[{qecRows = qecRows, qecLegend = qecLegend, qecEtas = qecEtas, qecRates = qecRates},
 Manipulate[Legended[qecRows[{qecEta, qecRate}], Placed[qecLegend, Below]],
  Grid[{{Control[{{qecEta, 1, "detector efficiency η"}, qecEtas, SetterBar}], Control[{{qecRate, 0.05, "flip rate γ/κ"}, qecRates, SetterBar}]}}, Spacings -> 4]]]
```

With a perfect detector and rare flips the filter moves to the new syndrome soon after each flip and its best guess is almost always right. With a poor detector and frequent flips it rarely settles before the next flip arrives. In every case the gap falls in proportion to the time step, so it is the integration's error: the four probabilities follow the full state, coherences included, and need nothing else from it. They close on themselves because both stabilizers are diagonal in the sectors and each flip only moves probability from one sector to another. For the same reason, records drawn from a true flip history have the same statistics as records from the quantum model.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Markovian feedback that protects the code space, and what the detector costs

Three qubits carry the bit-flip code, $|0_L\rangle=|000\rangle$, $|1_L\rangle=|111\rangle$, and each qubit flips at rate $\gamma$ through the channel $\hat c_i=\sqrt\gamma X_i$. Each channel is monitored by homodyne detection of efficiency $\eta$, and its current $I_i$ drives the qubits at once through the Hamiltonian $I_i\hat F_i$, with $\hat F_i=i\sqrt\gamma\,{}[\hat\Pi_\mathrm{code},X_i]$ and $\hat\Pi_\mathrm{code}=|000\rangle\langle000|+|111\rangle\langle111|$. How well does this loop hold an encoded state?

With $\hbar=1$ and $\mathcal D[\hat c]\rho=\hat c\rho\hat c^\dagger-\frac12\{\hat c^\dagger\hat c,\rho\}$, the Wiseman-Milburn feedback master equation is

$$
\dot\rho=-\frac{i}{2}\sum_{i=1}^3\left[\hat c_i^\dagger\hat F_i+\hat F_i\hat c_i,\rho\right]+\sum_{i=1}^3\mathcal D[\hat c_i-i\hat F_i]\rho+\frac{1-\eta}{\eta}\sum_{i=1}^3\mathcal D[\hat F_i]\rho .
$$

The raw current is fed back as it arrives, with no filtering, no estimate of the state and no delay, which is what keeps this equation closed in $\rho$; the part of the current that is detector noise reaches the qubits too, as the last term. The two identities $\hat c_i\hat F_i+\hat F_i\hat c_i=0$, which removes the Hamiltonian term, and $(\hat c_i-i\hat F_i)\hat\Pi_\mathrm{code}=0$ both hold exactly. The master equation evolves from $|0_L\rangle$ and from $|+_L\rangle=(|000\rangle+|111\rangle)/\sqrt2$, with $\mathcal D[\hat c_i]$ alone as the case without feedback; time is in units of $1/\gamma$:

```wl
fbFlips = Normal[QuantumOperator[#]["Matrix"]] & /@ {"XII", "IXI", "IIX"};
fbCode = DiagonalMatrix[UnitVector[8, 1] + UnitVector[8, 8]];
fbF = I (fbCode . # - # . fbCode) & /@ fbFlips;
fbChecks = {Max[Abs[MapThread[#1 . #2 + #2 . #1 &, {fbFlips, fbF}]]], Max[Abs[(# . fbCode) & /@ (fbFlips - I fbF)]]};
fbRho[ops_, rates_, ψ_] := Normal[QuantumEvolve[QuantumOperator[ConstantArray[0, {8, 8}]], QuantumOperator /@ ops -> rates, QuantumState[ψ], {t, 0, 2}]["DensityMatrix"]];
fbEtas = {1, 1/2, 1/10, 1/50};
fbStarts = {{Row[{"|", Subscript["0", "L"], "⟩"}], UnitVector[8, 1]}, {Row[{"|", Subscript["+", "L"], "⟩"}], (UnitVector[8, 1] + UnitVector[8, 8])/Sqrt[2]}};
fbCurves[ψ_] := Append[Function[η, fbRho[Join[fbFlips - I fbF, fbF], Join[{1, 1, 1}, ConstantArray[(1 - η)/η, 3]], ψ]] /@ fbEtas, fbRho[fbFlips, {1, 1, 1}, ψ]];
fbObs = {{Re[Tr[fbCode . #1]] &, "code-space population"}, {Re[Conjugate[#2] . #1 . #2] &, "fidelity with the start"}};
fbPanel[{name_, ψ_}, ρs_, {obs_, label_}] := With[{tag = Style[Row[{"start ", name}], 13, FontFamily -> "Helvetica"]},
   Plot[Evaluate[obs[#, ψ] & /@ ρs], {t, 0, 2}, PlotRange -> {0, 1.2}, AspectRatio -> 1/2, Frame -> True, FrameLabel -> {"γ t", label},
    FrameTicks -> {{Range[0, 1, 0.2], None}, {Automatic, None}}, GridLines -> {Automatic, Range[0, 1, 0.2]}, GridLinesStyle -> Directive[GrayLevel[0.8], Dashed],
    Epilog -> Inset[tag, Scaled[{0.97, 0.96}], {Right, Top}], PlotStyle -> {Automatic, Automatic, Automatic, Automatic, AbsoluteThickness[3.5]},
    PlotLegends -> Append[Row[{"η = ", #}] & /@ fbEtas, "no feedback"], ImageSize -> 340]];
fbPanels = Function[start, With[{ρs = fbCurves[Last[start]]}, fbPanel[start, ρs, #] & /@ fbObs]] /@ fbStarts;
Labeled[Legended[PlotGrid[Transpose[Map[First, fbPanels, {2}]], ImageSize -> 760], legendBelow[Last[fbPanels[[1, 1]]]]],
  Style[Row[{Subscript["c", "i"], Subscript["F", "i"], " + ", Subscript["F", "i"], Subscript["c", "i"], " = ", fbChecks[[1]], "   and   (",
     Subscript["c", "i"], " - i ", Subscript["F", "i"], ")", Subscript["Π", "code"], " = ", fbChecks[[2]], "   for each qubit i = 1, 2, 3"}], 14, FontFamily -> "Helvetica"], Top]
```

At $\eta=1$ the two identities make every encoded state a fixed point, for $|0_L\rangle$ and $|+_L\rangle$ alike. Below unit efficiency the $\mathcal D[\hat F_i]$ term feeds detector noise back onto the qubits: the plateau falls with $\eta$, and at low efficiency, with this feedback strength, the loop empties the code space faster than no feedback. From $|0_L\rangle$ the fidelity falls below the code-space population while from $|+_L\rangle$ the two agree: the noisy loop also flips $|0_L\rangle$ into $|1_L\rangle$, which no stabilizer sees.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Compiling a rotation from H and T gates: each T gate buys a fixed factor in precision

A fault-tolerant qubit offers only a finite gate set: the Clifford gates generated by $H$ and $S=T^2$, which are cheap, and the $T$ gate, which is expensive. Any target rotation $U$ can only be approximated by a product of them. How many $T$ gates does a given precision cost?

With $T=\mathrm{diag}(1,e^{i\pi/4})$ and $S=T^2$, every product of these gates equals, up to a global phase, a normal form (Matsumoto and Amano)

$$
W=\left(T\ \text{or}\ 1\right)\,s_1s_2\cdots s_k\,C,\qquad s_i\in\{HT,\ SHT\},\qquad C\in\text{the 24 Cliffords},
$$

whose $T$ count is $k$ or $k+1$, so from $T$ count 1 on, each extra $T$ gate doubles the number of candidates. The error of $W$ is $\epsilon=1-\frac12|\mathrm{Tr}\,W^\dagger U|$, which ignores the global phase and grows as the square of the rotation angle of $W^\dagger U$. Rotations up to phase fill a three-dimensional space, so $2^T$ circuits cannot come closer to a typical target than a distance $\propto2^{-T/3}$, which sets the floor $\epsilon\propto2^{-2T/3}$; whether the circuits spread evenly enough to reach it is what the enumeration tests.

Every normal form up to $T$ count 16 is enumerated for 200 rotations drawn uniformly at random, and the smallest error so far is kept at each $T$ count; the orange line has the slope $2^{-2T/3}$ and passes through the median at $T$ count 10:

```wl
{ctH, ctS, ctT} = N[{HadamardMatrix[2], DiagonalMatrix[{1, I}], DiagonalMatrix[{1, Exp[I π/4]}]}];
ctMatrix[w_] := Fold[#2 . #1 &, IdentityMatrix[2], w /. {"H" -> ctH, "S" -> ctS}];
ctCliffords = DeleteDuplicatesBy[Catenate[Table[Tuples[{"H", "S"}, l], {l, 0, 6}]],
   With[{f = Flatten[ctMatrix[#]]}, Round[f Exp[-I Arg[SelectFirst[f, Abs[#] > 0.1 &]]], 10.^-6]] &];
ctLevels = NestList[Join[{#[[1]] . ctH . ctT, Join[{"T", "H"}, #[[2]]]} & /@ #, {#[[1]] . ctS . ctH . ctT, Join[{"T", "H", "S"}, #[[2]]]} & /@ #] &, {{IdentityMatrix[2], {}}}, 16];
ctForms[t_] := If[t == 0, ctLevels[[1]], Join[ctLevels[[t + 1]], {ctT . #[[1]], Append[#[[2]], "T"]} & /@ ctLevels[[t]]]];
ctOverlaps[t_, u_] := Abs[Flatten[Conjugate[ctForms[t][[All, 1]]], {{1}, {2, 3}}] . Transpose[Flatten[u . ConjugateTranspose[ctMatrix[#]]] & /@ ctCliffords]]/2;
ctRandomU[] := With[{q = Normalize[RandomVariate[NormalDistribution[], 4]]}, {{q[[1]] + I q[[2]], q[[3]] + I q[[4]]}, {-q[[3]] + I q[[4]], q[[1]] - I q[[2]]}}];
ctTargets = BlockRandom[SeedRandom[7]; Table[ctRandomU[], 200]];
ctBest = Transpose[Table[Max[ctOverlaps[t, u]], {t, 0, 16}, {u, ctTargets}]];
ctRun = 1 - FoldList[Max, #] & /@ ctBest;
ctQuant = Transpose[Quantile[#, {0.1, 0.5, 0.9}] & /@ Transpose[ctRun]];
ctLine = ctQuant[[2, 11]] 2^(-2 (Range[0, 16] - 10)/3);
Legended[Show[ListLogPlot[ctRun, DataRange -> {0, 16}, Joined -> True, PlotStyle -> Directive[GrayLevel[0.6, 0.25], Thin]],
  ListLogPlot[{ctQuant[[1]], ctQuant[[3]]}, DataRange -> {0, 16}, Joined -> True, PlotStyle -> Directive[ColorData[97][1], Dashed], Filling -> {1 -> {2}}, FillingStyle -> Directive[ColorData[97][1], Opacity[0.12]]],
  ListLogPlot[{ctQuant[[2]], ctLine}, DataRange -> {0, 16}, Joined -> True, PlotStyle -> {Directive[ColorData[97][1], Thick], Directive[ColorData[97][4], Opacity[.7], Thick]},
   PlotMarkers -> {{Automatic, 8}, None}],
  PlotRange -> {{0, 16}, {Log[10^-6], Log[0.3]}}, Frame -> True, FrameLabel -> {"T count", "best error ε"}, GridLines -> {Automatic, Log[10.^Range[-6, -1]]}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.35], Thin],
  PlotLabel -> Row[{Length[ctTargets], " random targets, all normal forms over ", Length[ctCliffords], " Cliffords"}], AspectRatio -> 1/2, ImageSize -> 520],
 Placed[LineLegend[{Directive[GrayLevel[0.6], Thin], Directive[ColorData[97][1], Dashed], Directive[ColorData[97][1], Thick], Directive[ColorData[97][4], Opacity[.7], Thick]},
   {"one target", "10% and 90% of targets", "median", Row[{"slope ", Superscript[2, "−2T/3"], " from counting circuits"}]}], Right]]
```

Each target's error falls in irregular steps: a target can sit far from every circuit of one $T$ count and close to one of the next. Across targets the steps average out, and the median follows the counting slope over the range shown, gaining about one decimal digit in $\epsilon$ for every five $T$ gates. The band shows how far luck moves a single target: a few $T$ gates either way.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Light

The examples in this section are photons in linear optical networks and light in phase space. QF gives output probabilities from the network's transfer matrix with photon numbers, angles and phases left symbolic, and the Wigner functions of the states.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A beam splitter with symbolic photon numbers

$n_1$ and $n_2$ photons enter the two ports of a beam splitter with mixing angle $\theta$ and phase $\varphi$. What is the probability of finding $m$ photons in the first output, for every photon number at once?

The splitter rotates the creation operators,

$$
\hat a_1^\dagger\to\cos\theta\,\hat a_1^\dagger+e^{i\varphi}\sin\theta\,\hat a_2^\dagger,\qquad \hat a_2^\dagger\to-e^{-i\varphi}\sin\theta\,\hat a_1^\dagger+\cos\theta\,\hat a_2^\dagger,
$$

so it conserves the photon number, and the second output holds $n_1+n_2-m$ photons.

QF returns the amplitude for general $n_1,n_2,m,\theta,\varphi$ as a Wigner D-function, and the probability as its square:

```wl
splitterAmp = Simplify[QuantumInterferometer[{"BS"[θ, φ]}]["Amplitude", {n1, n2} -> {m, n1 + n2 - m}] /.
    WignerD[{j_, a_, b_}, ψ_, β_, χ_] :> Exp[I (a ψ + b χ)] WignerD[{j, a, b}, β], Element[{θ, φ, n1, n2, m}, Reals]]
splitterProb = Simplify[Abs[splitterAmp]^2, Element[{θ, φ, n1, n2, m}, Reals]]
```

Each panel fixes the input photon numbers:

```wl
twinPanel[{p_, q_}] := Plot[Evaluate[FullSimplify[Table[splitterProb /. {n1 -> p, n2 -> q, m -> j}, {j, 0, p + q}], Element[{θ, φ}, Reals]]],
    {θ, 0, π/2}, PlotRange -> {-0.04, 1.04}, PlotStyle -> Table[Directive[If[OddQ[j], Dashed, {}], AbsoluteThickness[If[j < 2, 3, 1.6]]], {j, 0, p + q}],
    PlotLegends -> Table[Row[{"m = ", j}], {j, 0, p + q}], GridLines -> {{π/4}, None}, Frame -> True,
    FrameTicks -> {{Automatic, None}, {{#, Style[TraditionalForm[#], FractionBoxOptions -> {Beveled -> True}]} & /@ Range[0, π/2, π/8], None}}, FrameLabel -> {"θ", "P(m)"},
    PlotLabel -> Row[{"|", p, ", ", q, "⟩  →  |m, ", p + q, " − m⟩"}], ImageSize -> 270];
twinPanels = twinPanel /@ {{1, 1}, {2, 2}, {1, 3}};
Legended[PlotGrid[{First /@ twinPanels}, ImageSize -> 810], legendBelow[Last[Last[twinPanels]]]]
```

The probability is the squared Wigner D-function of a spin $j=(n_1+n_2)/2$ rotated by $2\theta$: two modes holding $N$ photons carry the spin-$N/2$ representation of SU(2), and the splitter rotates it. The phase $\varphi$ enters the D-function only as a phase factor, so it drops out of the probability. At the balanced splitter (vertical line) equal inputs never leave an odd $m$ (dashed), with $n_1=n_2=1$ the Hong-Ou-Mandel zero, while one and three photons never leave two in each output.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Hong-Ou-Mandel with partially distinguishable photons

Two photons meet at a beam splitter with angle $\theta$ and phase $\varphi$. The first is horizontally polarized, the second is polarized at angle $\chi$, $\cos\chi\,|H\rangle+\sin\chi\,|V\rangle$, and the detectors do not resolve polarization. How does the coincidence probability depend on $\theta$ and $\chi$?

Each port carries two modes, one per polarization, so the network has four wires, ordered $1_H,2_H,2_V,1_V$. The one physical splitter, the transformation of the beam-splitter example, treats both polarizations alike, so it appears twice: once on the $H$ pair and once on the $V$ pair. Both photons enter horizontal, $|1,1,0,0\rangle$, and a waveplate on port 2 turns the second photon's polarization by $\chi$: it mixes $2_H$ and $2_V$ as a beam splitter of angle $\chi$ would. The coincidence sums over the four outputs $T$ with one photon in each port:

$$
P_\mathrm{coinc}=\sum_T\left|A(1100\to T)\right|^2 .
$$

The network, each wire labeled by its port and polarization, and $P_\mathrm{coinc}(\theta,\chi)$ for general $\theta,\varphi,\chi$ in closed form; it depends on $\chi$ only through $\cos2\chi$, so it is written as its value at $\chi=0$ times $\cos^2\chi$ plus its value at $\chi=\pi/2$ times $\sin^2\chi$:

```wl
polarizedHOM = QuantumInterferometer[{"BS"[χ, 0] -> {2, 3}, "BS"[θ, φ] -> {1, 2}, "BS"[θ, φ] -> {4, 3}}, 4];
homRaw = FullSimplify[ComplexExpand[Total[Abs[polarizedHOM["Amplitude", {1, 1, 0, 0} -> #]]^2 & /@ {{1, 1, 0, 0}, {1, 0, 1, 0}, {0, 1, 0, 1}, {0, 0, 1, 1}}]], Element[{θ, φ, χ}, Reals]];
partialCoincidence = FullSimplify[homRaw /. χ -> 0] Cos[χ]^2 + FullSimplify[homRaw /. χ -> π/2] Sin[χ]^2;
homDiagramPol = Show[polarizedHOM["Diagram", "ShowParameters" -> True, "ModeLabels" -> False],
   Graphics[{MapThread[Text[Column[{Style[#1, 11, FontFamily -> "Helvetica"], Style[#3, 10, Italic, GrayLevel[0.5]]}, Alignment -> Right, Spacings -> 0], {0.12, -#2}, {1, 0}] &,
       {{"port 1, H", "port 2, H", "port 2, V", "port 1, V"}, Range[4], {"photon a", "photon b", "empty", "empty"}}],
     Text[Style["waveplate on port 2", 11, GrayLevel[0.4]], {1, -3.3}, {0, 1}],
     Text[Style["H pair", 11, GrayLevel[0.4]], {2.85, -1.5}, {-1, 0}], Text[Style["V pair", 11, GrayLevel[0.4]], {2.85, -3.5}, {-1, 0}]}],
   PlotRange -> {{-1.05, 3.4}, {-4.75, -0.3}}, ImageSize -> 300];
homPlotPol = Plot[Evaluate[Table[partialCoincidence /. χ -> c, {c, {0, π/8, π/4, 3 π/8, π/2}}]], {θ, 0, π/2}, PlotRange -> {0, 1},
  PlotLegends -> {"χ = 0 (identical)", "χ = π/8", "χ = π/4", "χ = 3π/8", "χ = π/2 (orthogonal)"}, GridLines -> {{π/4}, None},
  Epilog -> Text[Style["balanced", 10, GrayLevel[0.5]], {π/4, 0.97}, {-1.1, 0}],
  Frame -> True, FrameTicks -> {{Automatic, None}, {Range[0, π/2, π/8], None}}, FrameLabel -> {"splitter angle θ", "P(one photon in each output port)"},
  PlotLabel -> Row[{"P = ", TraditionalForm[partialCoincidence]}], ImageSize -> 400];
Grid[{{homDiagramPol, homPlotPol}}, Spacings -> 3, Alignment -> Center]
```

The phase $\varphi$ drops out, and $P_\mathrm{coinc}$ is the Hong-Ou-Mandel curve $\cos^22\theta$ of identical photons weighted by the squared polarization overlap $\cos^2\chi$, plus the distinguishable-particle curve $(3+\cos4\theta)/4=\cos^4\theta+\sin^4\theta$ weighted by $\sin^2\chi$. At the balanced splitter only the second term survives, so the depth of the dip measures how indistinguishable the photons are.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The nonlinear sign gate of Knill, Laflamme and Milburn

A signal mode holds $c_0|0\rangle+c_1|1\rangle+c_2|2\rangle$; two auxiliary modes hold one photon and the vacuum. A three-mode network mixes them, and a run is kept only when the auxiliary detectors count $\{1,0\}$. A phase shifter cannot flip the sign of $c_2$ alone, since it multiplies $|n\rangle$ by $e^{in\alpha}$. What does the kept signal look like?

The network maps $\hat a_k^\dagger\to\sum_j W_{jk}\,\hat a_j^\dagger$ with the real transfer matrix (Knill, Laflamme and Milburn, Nature 409, 46 (2001))

$$
W=\begin{pmatrix}1-\sqrt2 & 2^{-1/4} & \sqrt{3/\sqrt2-2}\\ 2^{-1/4} & 1/2 & (1-\sqrt2)/2\\ \sqrt{3/\sqrt2-2} & (1-\sqrt2)/2 & \sqrt2-1/2\end{pmatrix},
$$

so the kept signal is $\sum_m c_m\langle m,1,0|\hat U|m,1,0\rangle\,|m\rangle$.

$W$ factors into beam splitters with angle and phase $(\theta,\varphi)$ and a phase shifter $\phi$; the network gives these amplitudes for $m=0,\dots,5$ (bars), and the first three give the kept signal for general $c_0,c_1,c_2$:

```wl
nsGate = QuantumInterferometer[{{1 - Sqrt[2], 2^(-1/4), Sqrt[3/Sqrt[2] - 2]}, {2^(-1/4), 1/2, (1 - Sqrt[2])/2},
    {Sqrt[3/Sqrt[2] - 2], (1 - Sqrt[2])/2, Sqrt[2] - 1/2}}];
nsAmps = RootReduce[Table[nsGate["Amplitude", {m, 1, 0} -> {m, 1, 0}], {m, 0, 5}]];
keptSignal = {c0, c1, c2} Take[nsAmps, 3];
Grid[{{nsGate["Diagram", "ShowParameters" -> True], Labeled[BarChart[nsAmps, PlotLabels -> Range[0, 5],
    LabelingFunction -> (Placed[TraditionalForm[#1], If[#1 >= 0, Above, Below]] &), PlotRangePadding -> {Automatic, Scaled[0.12]},
    Frame -> True, FrameLabel -> {"signal photons m", "⟨m,1,0|U|m,1,0⟩"}, ImageSize -> 400],
  Style[Row[{"kept signal  ", TraditionalForm[keptSignal]}], 13, FontFamily -> "Helvetica"], Top]}}, Spacings -> 3, Alignment -> Center]
```

For zero, one and two signal photons the amplitudes have the same size and only the two-photon one is negative, so every input comes out with $c_2$ flipped and the same overall factor: the gate works with a probability that does not depend on the input. From three signal photons on the amplitudes differ, so the gate works only on inputs with at most two signal photons.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A Schrödinger-cat state in phase space

A field mode holds a superposition of two coherent states of real amplitude $\pm\alpha$ with relative phase $\varphi$. What sets the sign of its Wigner function at the origin?

With $\hat x=(\hat a+\hat a^\dagger)/\sqrt2$ and $\hat p=(\hat a-\hat a^\dagger)/i\sqrt2$, the state and its Wigner function are

$$
|\psi\rangle=\mathcal N\left(|\alpha\rangle+e^{i\varphi}|{-\alpha}\rangle\right),\qquad W(x,p)=\frac1\pi\int\langle x-y|\rho|x+y\rangle\,e^{2ipy}\,dy,
$$

and the photon-number parity $\hat\Pi=(-1)^{\hat a^\dagger\hat a}$ maps $|\alpha\rangle$ to $|{-\alpha}\rangle$. Since $\pi W(0,0)=\langle\hat\Pi\rangle$ for every state, the overlap $\langle\alpha|{-\alpha}\rangle=e^{-2\alpha^2}$ fixes the value at the origin in closed form.

$\pi W$ comes from the state on a Fock space cut at 40 levels. Each panel draws $\pi W$ as a surface over phase space, prints the computed $\pi W(0,0)$, marks $W=0$ with a white line, and all panels share one color scale:

```wl
catColors = ColorData["Rainbow"][(1 + #)/2] &;
catParity = (Exp[-2 α^2] + Cos[φ])/(1 + Exp[-2 α^2] Cos[φ]);
catPanel[amp_, phase_] := With[{w = WignerRepresentation[CatState[40][N[amp], N[phase]], {-6, 6}, {-3, 3}, "GridSize" -> 61]},
  Plot3D[π w[x, p], {x, -6, 6}, {p, -3, 3}, PlotRange -> {-1, 1}, ColorFunction -> (catColors[#3] &), ColorFunctionScaling -> False, PlotPoints -> 40, MaxRecursion -> 0,
   MeshFunctions -> {#3 &}, Mesh -> {{0}}, MeshStyle -> Directive[GrayLevel[1, 0.4], Thin], BoxRatios -> {2, 1, 0.6}, AxesLabel -> {"x", "p", "πW"}, ImageSize -> 250, PlotStyle -> Opacity[0.6], ViewPoint -> {1.2, -2.8, 0.55},
   PlotLabel -> Column[{Row[{"α = ", amp, ",   φ = ", StringReplace[ToString[phase, InputForm], "Pi" -> "π"]}], Row[{"computed πW(0,0) = ", NumberForm[π w[0, 0], {4, 3}, ExponentFunction -> (Null &)]}]}, Alignment -> Center]]];
catPanels = AssociationThread[{0.5, 1, 1.5, 2, 3}, Outer[catPanel, {0.5, 1, 1.5, 2, 3}, {0, π/2, π}]];
With[{controlNotes = controlNotes, catPanels = catPanels, catColors = catColors, catParity = catParity},
 Manipulate[Labeled[Legended[Grid[{catPanels[amp]}], Placed[BarLegend[{catColors, {-1, 1}}, LegendLabel -> "πW", LegendLayout -> "Row"], Below]],
  Style[Row[{"πW(0,0) = ⟨Π⟩ = ", TraditionalForm[catParity]}], 15, FontFamily -> "Helvetica"], Top], Row[{Control[{{amp, 2, "α"}, Keys[catPanels], SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"α = ", amp}], <|0.5 -> "The lobes merge: the even cat looks like the vacuum, the odd cat like one photon, with πW = −1 at the origin.",
        1 -> "The lobes still overlap, yet W already goes negative.",
        1.5 -> "The lobes separate; fringes form between them.",
        2 -> "The lobes pull apart; fringes between them dip below zero.",
        3 -> "The fringes crowd together; at φ = π/2 the origin sits near zero."|>[amp]}}]]}]]]
```

As $\alpha$ grows the two coherent lobes, at $x=\pm\sqrt2\,\alpha$, pull apart, and the fringes along $p$ between them dip below zero and crowd together; at $\alpha=1$ the lobes still overlap, yet $W$ already goes negative. In every panel the origin takes the closed-form parity: $\varphi=0$ and $\varphi=\pi$ give parity eigenstates at $\pm1$ for every $\alpha$, while at $\varphi=\pi/2$ the parity is $e^{-2\alpha^2}$, so the origin falls toward zero, between a positive and a negative fringe.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Entanglement and many-body states

The examples in this section treat entanglement as a resource and as a property of many-body states. QF builds the states, and TensorNetworks contracts the matrix-product states without writing out their amplitudes.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### GHZ and W: arbitrarily close, never connected

Three qubits share $|\mathrm{GHZ}\rangle=(|000\rangle+|111\rangle)/\sqrt2$ or $|W\rangle=(|001\rangle+|010\rangle+|100\rangle)/\sqrt3$. Each party may apply any invertible local map $A$, $B$, $C$, even one that succeeds only with some probability. Can such maps turn GHZ into W?

The maps act as $|\psi\rangle\mapsto(A\otimes B\otimes C)|\psi\rangle$, followed by normalization. Cut the amplitudes $\psi_{ijk}$ into the two $2\times2$ slices $M_i=(\psi_{ijk})_{jk}$; the Cayley hyperdeterminant $\operatorname{Det}\psi$ is the discriminant of the quadratic $\det(M_0+xM_1)$, and the three-tangle of a normalized state is $\tau_3=4|\operatorname{Det}\psi|$. It is the part of qubit 1's entanglement that no single pair carries, $\tau_3=\mathcal C_{1|23}^2-\mathcal C_{12}^2-\mathcal C_{13}^2$ (Coffman, Kundu and Wootters), with $\mathcal C$ the concurrence: $1$ for GHZ and $\frac89-\frac49-\frac49=0$ for W. The family

$$
|\psi_\varepsilon\rangle\propto A_\varepsilon^{\otimes3}|\mathrm{GHZ}\rangle\propto(|0\rangle+\varepsilon|1\rangle)^{\otimes3}-|000\rangle,\qquad A_\varepsilon=\begin{pmatrix}1&-1\\ \varepsilon&0\end{pmatrix},\quad\det A_\varepsilon=\varepsilon,
$$

is GHZ moved by an invertible map for every $\varepsilon>0$.

With all twelve entries of $A$, $B$, $C$ left as symbols, both states are transformed and their hyperdeterminants compared with $(\det A\det B\det C)^2$; the family $\psi_\varepsilon$ gives $\tau_3$ and the infidelity with W in closed form:

```wl
hyperdetSLOCC[v_] := With[{m = ArrayReshape[Normal[v], {2, 2, 2}]}, With[{c = CoefficientList[Det[m[[1]] + \[FormalX] m[[2]]], \[FormalX], 3]}, c[[2]]^2 - 4 c[[1]] c[[3]]]];
{ghzSLOCC, wSLOCC} = QuantumState[#[3]] & /@ {"GHZ", "W"};
mapSLOCC = QuantumTensorProduct @@ (QuantumOperator[Array[#, {2, 2}]] & /@ {slA, slB, slC});
detsSLOCC = Times @@ (Det[Array[#, {2, 2}]] & /@ {slA, slB, slC});
{ghzMapSLOCC, wMapSLOCC} = Normal[mapSLOCC[#]["StateVector"]] & /@ {ghzSLOCC, wSLOCC};
{ghzRuleQ, wRuleQ} = {Expand[hyperdetSLOCC[ghzMapSLOCC] - detsSLOCC^2/4] === 0, Expand[hyperdetSLOCC[wMapSLOCC]] === 0};
vecSLOCC = ghzMapSLOCC /. Flatten[Thread[Flatten[Array[#, {2, 2}]] -> {1, -1, \[Epsilon], 0}] & /@ {slA, slB, slC}];
{tau3SLOCC, infidSLOCC} = Simplify[{4 Abs[hyperdetSLOCC[vecSLOCC]]/(vecSLOCC . vecSLOCC)^2, 1 - (Normal[wSLOCC["StateVector"]] . vecSLOCC)^2/(vecSLOCC . vecSLOCC)}, \[Epsilon] > 0];
Grid[{{Subscript["\[Tau]", 3], "=", TraditionalForm[tau3SLOCC]}, {Row[{"1 - |\[LeftAngleBracket]W|", Subscript["A", "\[Epsilon]"], Subscript["A", "\[Epsilon]"], Subscript["A", "\[Epsilon]"], "|GHZ\[RightAngleBracket]|²"}], "=", TraditionalForm[infidSLOCC]},
  {Row[{"Det((A\[CircleTimes]B\[CircleTimes]C)GHZ) = (det A det B det C", Superscript[")", 2], "/4"}], "", ghzRuleQ}, {"Det((A\[CircleTimes]B\[CircleTimes]C)W) = 0", "", wRuleQ}},
  Alignment -> {{Right, Center, Left}}, Spacings -> {1, 1.2}, BaseStyle -> {13, FontFamily -> "Helvetica"}]
```

The same rule for 300 invertible maps with random complex integer entries (left), and the family $\psi_\varepsilon$ as $\varepsilon$ shrinks (right):

```wl
(* 300 invertible maps with random Gaussian-integer entries, exact arithmetic *)
entriesSLOCC = Flatten[Array[#, {2, 2}] & /@ {slA, slB, slC}];
samplesSLOCC = BlockRandom[SeedRandom[3]; Select[Table[Thread[entriesSLOCC -> RandomInteger[{-9, 9}, 12] + I RandomInteger[{-9, 9}, 12]], 400], (detsSLOCC /. #) =!= 0 &][[;; 300]]];
pointsSLOCC = {Abs[detsSLOCC^2/4], Abs[hyperdetSLOCC[ghzMapSLOCC]], Abs[hyperdetSLOCC[wMapSLOCC]]} /. # & /@ samplesSLOCC;
sampledSLOCC = ListLogLogPlot[{pointsSLOCC[[All, {1, 2}]], {#, 10^3} & /@ pointsSLOCC[[All, 1]]},
 PlotRange -> {{10^5, 10^13}, {10^2.5, 10^13}}, Frame -> True,
 FrameLabel -> {Row[{"(det A det B det C", Superscript[")", 2], "/4"}], Row[{"hyperdeterminant of A\[CircleTimes]B\[CircleTimes]C|\[Psi]\[RightAngleBracket]"}]},
 FrameTicks -> {{{{10^3, 0}, {10^6, Superscript[10, 6]}, {10^9, Superscript[10, 9]}, {10^12, Superscript[10, 12]}}, None}, {Automatic, None}},
 Epilog -> {Directive[GrayLevel[0.5], Thin], Line[{{Log[10^5], Log[10^5]}, {Log[10^13], Log[10^13]}}]},
 PlotLegends -> Placed[{"|\[Psi]\[RightAngleBracket] = |GHZ\[RightAngleBracket]", "|\[Psi]\[RightAngleBracket] = |W\[RightAngleBracket]"}, {Left, Top}], PlotLabel -> "300 random invertible local maps", AspectRatio -> 0.75, ImageSize -> 400];
labelSLOCC = Style[Column[{Row[{Subscript["\[Tau]", 3], " = ", RawBoxes[ToBoxes[TraditionalForm[tau3SLOCC]]]}], Row[{"1 - |\[LeftAngleBracket]W|", Subscript["A", "\[Epsilon]"], Subscript["A", "\[Epsilon]"], Subscript["A", "\[Epsilon]"], "|GHZ\[RightAngleBracket]|²", " = ", RawBoxes[ToBoxes[TraditionalForm[infidSLOCC]]]}]}, Spacings -> 1], 12];
limitSLOCC = LogLogPlot[{tau3SLOCC, infidSLOCC}, {\[Epsilon], 10^-3, 1}, PlotRange -> All, PlotStyle -> {Automatic, Dashed},
 PlotLegends -> Placed[{Subscript["\[Tau]", 3], Row[{"1 - |\[LeftAngleBracket]W|", Subscript["A", "\[Epsilon]"], Subscript["A", "\[Epsilon]"], Subscript["A", "\[Epsilon]"], "|GHZ\[RightAngleBracket]|²"}]}, {Right, Bottom}], Frame -> True, FrameLabel -> {Row[{"\[Epsilon],  the determinant of each local map ", Subscript["A", "\[Epsilon]"]}], None}, AspectRatio -> 0.75, ImageSize -> 400, PlotLabel -> Row[{"GHZ moved toward W by ", Subscript["A", "\[Epsilon]"], " on every qubit"}],
 Epilog -> Inset[labelSLOCC, Scaled[{0.04, 0.96}], {Left, Top}]];
Grid[{{sampledSLOCC, limitSLOCC}}, Spacings -> 4, Alignment -> Top]
```

Every map puts GHZ on the diagonal: its hyperdeterminant after the map is exactly $(\det A\det B\det C)^2/4$, never zero, while W stays at zero. An invertible local map can only rescale the hyperdeterminant, so no such map turns GHZ into W or W into GHZ. Yet $\psi_\varepsilon$ stays in the GHZ class for every $\varepsilon>0$ while its $\tau_3$ and its infidelity with W both vanish as $\varepsilon^2$: W sits on the boundary of the GHZ class, reached only in the limit $\varepsilon\to0$, where the map stops being invertible.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The AKLT ring: from a matrix-product state to an exact ground state

A ring of $n$ spin-1 sites is given as a matrix-product state whose bonds carry only two states: one $2\times2\times3$ tensor, repeated on every site, fixes all $3^n$ amplitudes. Without ever writing those amplitudes out, is it an exact ground state of the AKLT Hamiltonian, and how fast do the spin correlations $C(r)=\langle S^z_1S^z_{1+r}\rangle$ decay with the distance $r$?

With $\mathbf S_i$ the spin-1 operators of site $i$ and $P_2^{(i,i+1)}$ the projector of two neighbors onto total spin 2, the Hamiltonian on the ring ($n+1\equiv1$) and the state, with $m_i\in\{1,0,-1\}$ and $\sigma^\pm$, $\sigma^z$ the $2\times2$ Pauli matrices, are

$$
\hat H=\sum_{i=1}^{n}P_2^{(i,i+1)}=\sum_{i=1}^{n}\Big[\frac12\,\mathbf S_i\cdot\mathbf S_{i+1}+\frac16\,(\mathbf S_i\cdot\mathbf S_{i+1})^2+\frac13\Big],\qquad
|\psi\rangle=\sum_{m_1\dots m_n}\operatorname{Tr}\big(A^{m_1}\cdots A^{m_n}\big)\,|m_1\cdots m_n\rangle,
$$

$$
A^{1}=\sqrt{\frac23}\,\sigma^+,\qquad A^{0}=-\sqrt{\frac13}\,\sigma^z,\qquad A^{-1}=-\sqrt{\frac23}\,\sigma^- .
$$

$\hat H$ is a sum of projectors, so $\langle\hat H\rangle\ge0$ and any state with $\langle\hat H\rangle=0$ is a ground state. The transfer matrix $\sum_mA^m\otimes\overline{A^m}$ has eigenvalues $1$ and $-\frac13$ (three times), which gives on the ring

$$
C(r)=\frac43\,\frac{(-1/3)^r+(-1/3)^{n-r}}{1+3\,(-1/3)^n}.
$$

Every quantity is a closed network: the ring of $A$ joined to its mirror ring of $\overline{A}$ at every spin-1 leg gives $\langle\psi|\psi\rangle$, and inserting $P_2$ on one bond, or $S^z$ on two sites, gives the energy and $C(r)$. Each network contracts exactly to one number at a cost that grows only linearly with $n$, so the control offers rings of $n=6$, $12$ and $24$ sites, and the same contraction stays exact for rings of hundreds of sites. The drawings show the energy network, with $P_2$ joining sites 1 and 2, and the correlation network for the chosen $r$, with $S^z$ on sites 1 and $1+r$; that $r$ is marked in red on the plot. In the plot the markers are $|C(r)|$ for $r\le n/2$, the dashed line is the ring formula and the dotted line the infinite-chain law $\frac43\,3^{-r}$, and the energy is in the plot label:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`"];
akltA = Transpose[{Sqrt[2/3] {{0, 1}, {0, 0}}, -Sqrt[1/3] {{1, 0}, {0, -1}}, -Sqrt[2/3] {{0, 0}, {1, 0}}}, {3, 1, 2}];
akltSS = Total[KroneckerProduct[#, #] & /@ {{{0, 1, 0}, {1, 0, 1}, {0, 1, 0}}/Sqrt[2], {{0, -I, 0}, {I, 0, -I}, {0, I, 0}}/Sqrt[2], DiagonalMatrix[{1, 0, -1}]}];
akltP2 = (akltSS + IdentityMatrix[9]) . (akltSS + 2 IdentityMatrix[9])/6;
akltClosed[n_, ops_] := With[{opSites = Flatten[Last /@ ops]}, TensorNetwork[Join[ConstantArray[akltA, n], ConstantArray[Conjugate[akltA], n], First /@ ops],
    Join[Table[{i, Mod[i, n] + 1, 2 n + i}, {i, n}], Table[{n + i, n + Mod[i, n] + 1, If[MemberQ[opSites, i], 3 n + i, 2 n + i]}, {i, n}], Join[3 n + #, 2 n + #] & /@ (Last /@ ops)]]];
akltData[n_] := With[{norm = TensorNetworkContract[akltClosed[n, {}]]},
   {Simplify[n TensorNetworkContract[akltClosed[n, {{ArrayReshape[akltP2, {3, 3, 3, 3}], {1, 2}}}]]/norm],
    Table[{r, Abs[Simplify[TensorNetworkContract[akltClosed[n, {{DiagonalMatrix[{1, 0, -1}], {1}}, {DiagonalMatrix[{1, 0, -1}], {1 + r}}}]]/norm]]}, {r, n/2}]}];
akltAngle[n_, i_] := 2 Pi i/n;
akltDraw[n_, ops_, size_] := With[{opSites = Flatten[Last /@ ops], net = akltClosed[n, ops]},
   With[{opEdges = Join[3 n + #, 2 n + #] & /@ (Last /@ ops)},
    RawBoxes[ToBoxes[net["Hypergraph", VertexLabels -> Join[Table[2 n + i -> i, {i, opSites}], Table[k -> None, {k, Complement[Join[Range[3 n], 3 n + opSites], 2 n + opSites]}]], EdgeLabels -> None, ImageSize -> size,
      EdgeStyle -> Join[Table[{i, Mod[i, n] + 1, 2 n + i} -> Directive[ColorData[97][1], Opacity[0.6]], {i, n}],
        Table[{n + i, n + Mod[i, n] + 1, If[MemberQ[opSites, i], 3 n + i, 2 n + i]} -> Directive[ColorData[97][3], Opacity[0.6]], {i, n}],
        Thread[opEdges -> Directive[Red, Opacity[0.9], AbsoluteThickness[5]]]],
      VertexCoordinates -> Join[Table[i -> 1.3 {Cos[2 Pi (i - 1/2)/n], Sin[2 Pi (i - 1/2)/n]}, {i, n}], Table[n + i -> 0.7 {Cos[2 Pi (i - 1/2)/n], Sin[2 Pi (i - 1/2)/n]}, {i, n}],
        Table[2 n + i -> If[MemberQ[opSites, i], 1.15, 1] {Cos[akltAngle[n, i]], Sin[akltAngle[n, i]]}, {i, n}],
        Table[3 n + i -> 0.85 {Cos[akltAngle[n, i]], Sin[akltAngle[n, i]]}, {i, opSites}]]]]]]];
akltP2Op = {{ArrayReshape[akltP2, {3, 3, 3, 3}], {1, 2}}};
akltSzOps[r_] := {{DiagonalMatrix[{1, 0, -1}], {1}}, {DiagonalMatrix[{1, 0, -1}], {1 + r}}};
akltView[n_, akltR_] := With[{res = akltResults[n]}, Grid[{{
    Labeled[akltDraw[n, akltP2Op, 200], Style[Column[{Row[{"\[LeftAngleBracket]\[Psi]|", OverHat["H"], "|\[Psi]\[RightAngleBracket] = n \[LeftAngleBracket]\[Psi]|", Subsuperscript["P", 2, "(1,2)"], "|\[Psi]\[RightAngleBracket]"}],
        "red: the projector on one bond"}, Alignment -> Center], 12, FontFamily -> "Helvetica"], Top],
    Labeled[akltDraw[n, akltSzOps[akltR], 200], Style[Column[{Row[{"\[LeftAngleBracket]\[Psi]|", Subsuperscript["S", 1, "z"], Subsuperscript["S", 1 + akltR, "z"], "|\[Psi]\[RightAngleBracket],  r = ", akltR}],
        Row[{"red: ", Superscript["S", "z"], " on two sites"}]}, Alignment -> Center], 12, FontFamily -> "Helvetica"], Top],
    ListLogPlot[{res[[2]], Table[{r, 4/3 (3^-r + 3^(r - n))/(1 + 3^(1 - n))}, {r, 1, n/2, 1/20}], Table[{r, 4/3 3^-r}, {r, 1, 12, 1/20}]}, Joined -> {False, True, True},
     PlotStyle -> {ColorData[97][1], Directive[ColorData[97][1], Dashed], Directive[Gray, Dotted]}, PlotMarkers -> {Automatic, None, None},
     PlotLegends -> Placed[{Row[{"contracted, ring of n = ", n}], "ring formula |C(r)|", Row[{"infinite chain  ", TraditionalForm[4/3], " ", TraditionalForm[HoldForm[3^-r]]}]}, {Right, Top}],
     PlotRange -> {{0.5, 12.5}, {10^-6, 1}}, Frame -> True, Epilog -> {Red, PointSize[0.025], Point[{akltR, Log[res[[2, akltR, 2]]]}]},
     FrameLabel -> {"distance r between the two sites", Row[{"spin correlation  |\[LeftAngleBracket]", Subsuperscript["S", 1, "z"], Subsuperscript["S", "1+r", "z"], "\[RightAngleBracket]|"}]},
     PlotLabel -> Row[{"\[LeftAngleBracket]", OverHat["H"], "\[RightAngleBracket] = ", res[[1]]}], ImageSize -> 400]}}, Spacings -> 2, Alignment -> Center]];
akltResults = AssociationMap[akltData, {6, 12, 24}];
akltViews = Association[Flatten[Table[{n, r} -> akltView[n, r], {n, {6, 12, 24}}, {r, n/2}]]];
With[{akltViews = akltViews}, Manipulate[akltViews[{akltN, Min[akltRr, akltN/2]}],
  Grid[{{"ring size n", Control[{{akltN, 6, ""}, {6, 12, 24}, SetterBar}]}, {"distance r", Dynamic[SetterBar[Dynamic[akltRr], Range[akltN/2]]]}}, Alignment -> {{Right, Left}, Center}],
  {{akltRr, 2}, ControlType -> None}]]
```

Contracting the ring with its mirror image gives an energy of exactly zero for each ring offered, so the two-state bonds hold an exact ground state of the AKLT chain, and the state's $3^n$ amplitudes are never written out. The correlations fall by a factor of three per site with alternating sign, a correlation length of $1/\ln3$ sites. The contracted values sit on the ring formula, which follows the infinite-chain line until near $r=n/2$. There the correlation reaching the other way round the ring, the term $(-1/3)^{n-r}$, has the same sign on an even ring, adds to it, and bends the curve up; on a larger ring the bend moves out to larger $r$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The entanglement profile of an Ising ring, and the bond dimension it needs

A ring of $n$ spins in its ground state is cut into a block of $\ell$ neighboring sites and the rest. How does the entanglement entropy $S(\ell)$ of the block depend on $\ell$, in the ordered phase, at the critical point and in the paramagnet, and what bond dimension $\chi$ does a matrix-product state of the ground state need? Keeping the $\chi$ largest Schmidt values at a cut discards the weight of the rest, and $\chi$ needed is the smallest number that keeps the discarded weight below $10^{-6}$.

With $Z_i$, $X_i$ the Pauli operators of site $i$ and $Z_{n+1}\equiv Z_1$, the transverse-field Ising ring is

$$
\hat H=-J\sum_{i=1}^{n}Z_iZ_{i+1}-h\sum_{i=1}^{n}X_i ,\qquad S(\ell)=-\operatorname{Tr}\rho_\ell\ln\rho_\ell ,
$$

with $\rho_\ell$ the reduced state of sites $1,\dots,\ell$. At the critical point $h=J$ conformal field theory predicts (Calabrese and Cardy 2004)

$$
S(\ell)=\frac c3\,\ln\!\Big[\frac n\pi\sin\frac{\pi\ell}{n}\Big]+s_1 ,
$$

with central charge $c$ and a constant $s_1$.

The ground state is found by exact diagonalization at $h/J=\frac12$, $1$ and $2$ and written as a matrix-product state with no truncation, and the Schmidt values at each bond give $S(\ell)$ and $\chi$. Top: the state for $n=16$, each bond drawn as thick as the $\chi$ it needs; bottom left: $S(\ell)$, with the critical profile fitted to the formula above; bottom right: $\chi$ at the middle cut for rings of $n=8$ to $18$:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`"];
epSite[n_, ops_Association] := KroneckerProduct @@ Table[Lookup[ops, k, IdentityMatrix[2, SparseArray]], {k, n}];
epH[n_, hJ_] := -Sum[epSite[n, <|i -> SparseArray[{{1, 0}, {0, -1}}], Mod[i, n] + 1 -> SparseArray[{{1, 0}, {0, -1}}]|>], {i, n}] - hJ Sum[epSite[n, <|i -> SparseArray[{{0, 1}, {1, 0}}]|>], {i, n}];
epGround[n_, hJ_] := First[Eigenvectors[-N[epH[n, hJ]], 1, Method -> {"Arnoldi", "Criteria" -> "RealPart"}]];
epStep[{_, r_}, _] := With[{m = ArrayReshape[r, {2 Length[r], Length[First[r]]/2}]},
   With[{svd = SingularValueDecomposition[m, Min[Dimensions[m]]]}, {svd[[1]], svd[[2]] . ConjugateTranspose[svd[[3]]]}]];
epMPS[n_, psi_] := With[{steps = Rest[FoldList[epStep, {None, {psi}}, Range[n - 1]]]},
   TensorNetwork[Join[{Transpose[steps[[1, 1]]]}, Transpose[ArrayReshape[#, {Length[#]/2, 2, Last[Dimensions[#]]}], {1, 3, 2}] & /@ steps[[2 ;;, 1]], {steps[[-1, 2]]}],
    Join[{{1, n + 1}}, Table[{i - 1, i, n + i}, {i, 2, n - 1}], {{n - 1, 2 n}}]]];
chiNeeded[s_, tol_] := Count[Reverse[Accumulate[Reverse[s^2]]], _?(# > tol &)];
epHs = {1/2, 1, 2}; epNames = {"ordered", "critical", "paramagnet"};
epMPS16 = epMPS[16, epGround[16, #]] & /@ epHs;
epChi = Table[chiNeeded[MPSSchmidtValues[mps, l], 10^-6], {mps, epMPS16}, {l, 15}];
epS = Table[MPSEntanglementEntropy[mps, l], {mps, epMPS16}, {l, 15}];
epFit = Fit[Transpose[{Log[(16/Pi) Sin[Pi Range[15]/16]], epS[[2]]}], {1, x}, x];
epDraw[k_] := Graphics[{
    MapIndexed[{ColorData[97][k], AbsoluteThickness[1 + 1.2 #1], Line[{{First[#2], 0}, {First[#2] + 1, 0}}], Black, Text[Style[#1, 9], {First[#2] + 0.5, 0.45}]} &, epChi[[k]]],
    Table[{GrayLevel[0.5], Line[{{i, 0}, {i, -0.6}}], EdgeForm[GrayLevel[0.3]], White, Disk[{i, 0}, 0.22]}, {i, 16}]},
   PlotLabel -> Row[{epNames[[k]], ",  h/J = ", epHs[[k]]}], ImageSize -> 560, PlotRange -> {{0.4, 16.6}, {-0.8, 0.8}}];
epChain = Labeled[Column[epDraw /@ Range[3], Spacings -> 1], Style[Row[{"n = 16: each bond drawn as thick as the bond dimension that keeps the discarded weight below ", Superscript[10, -6]}], 12, FontFamily -> "Helvetica"], Top];
epProfile = Legended[ListPlot[Append[Transpose[{Range[15], #}] & /@ epS, Table[{l, epFit /. x -> Log[(16/Pi) Sin[Pi l/16]]}, {l, 1, 15, 0.1}]],
   Joined -> {False, False, False, True}, PlotMarkers -> {Automatic, Automatic, Automatic, None}, PlotStyle -> {Automatic, Automatic, Automatic, Dashed},
   
   GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.5, 0.35], Thin], Frame -> True, FrameLabel -> {Row[{"block size ", Style["\[ScriptL]", FontFamily -> "Times"]}], Row[{"entanglement entropy S(", Style["\[ScriptL]", FontFamily -> "Times"], ")"}]}, PlotRange -> {0, 0.85}, PlotLabel -> "n = 16", ImageSize -> 330], Placed[LineLegend[{Directive[ColorData[97][4], Dashed]}, {Row[{"(c/3) ln[(n/\[Pi]) sin(\[Pi]", Style["\[ScriptL]", FontFamily -> "Times"], "/n)] + ", Subscript["s", 1], ", c = ", NumberForm[3 Coefficient[epFit, x], 2]}]}], Below]];
epNs = {8, 10, 12, 14, 16, 18};
epGrowth = Table[chiNeeded[MPSSchmidtValues[epMPS[n, epGround[n, hJ]], n/2], 10^-6], {hJ, epHs}, {n, epNs}];
epScale = ListPlot[Transpose[{epNs, #}] & /@ epGrowth, Joined -> True, PlotMarkers -> Automatic, PlotRange -> {0, 15},
   GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.5, 0.35], Thin], Frame -> True, FrameLabel -> {"ring size n", "bond dimension at the middle cut"}, PlotLabel -> Row[{"bond dimension needed (discarded weight < ", Superscript[10, -6], ")"}], ImageSize -> 330];
Column[{epChain, Row[{epProfile, Spacer[30], epScale}]}, Alignment -> Center, Spacings -> 2]
```

In the ordered phase and in the paramagnet the entropy hardly depends on $\ell$, so the bonds stop thickening a few sites in from the ends and $\chi$ at the middle cut stays the same as the ring grows. At $h/J=\frac12$ the entropy sits near $\ln2$, the one bit of a state close to the cat state $\frac1{\sqrt2}(|{\uparrow\cdots\uparrow}\rangle+|{\downarrow\cdots\downarrow}\rangle)$ that the symmetric ground state of a finite ring forms; at $h/J=2$ it stays small. At the critical point the entropy follows the conformal profile with $c$ close to $\frac12$, the central charge of the Ising critical point, and the bonds keep thickening toward the middle: over the sizes shown, the bond dimension a critical chain needs grows with the system, while that of a gapped chain does not.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Hardware and beyond

The examples in this section reach hardware and geometry: a transmon driven by a fast pulse, an entanglement pump run on an IBM processor, and the symmetry of the Riemann tensor. QF evolves the pulse Hamiltonian and models the processor gate by gate with noise channels, TransmonLab gives the transmon's circuit and spectrum, and the TensorNetworks Young projectors carry the tensor symmetry.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Leakage of a fast π pulse on a transmon

A transmon is an anharmonic oscillator, frequency $\omega_q$ and anharmonicity $\alpha<0$, so its $|1\rangle\to|2\rangle$ transition sits only $|\alpha|$ below the qubit transition. A π pulse of length $\tau$, detuned by $\delta$ from the qubit, is meant to flip $|0\rangle$ to $|1\rangle$, and a DRAG term of strength $\lambda$ adds a second quadrature that follows the slope of the envelope. Where does the population end up after the pulse: outside the qubit, or left in $|0\rangle$?

For a typical transmon, $E_J/E_C=50$ and $E_C/h=250$ MHz, TransmonLab draws the circuit and gives its lowest levels, drawn in the cosine potential of the junction:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/mbahram/TransmonLab/TransmonLab"]; Needs["TransmonLab`"];
(* E_J/E_C = 50 and E_C/h = 250 MHz, a typical transmon; energies in units of E_C *)
{tmRatio, tmEC} = {50, 0.25};
tmSpec = TransmonSpectrum[tmRatio, 4];
tmLevels = Min[Eigenvalues[N[CPBHamiltonian[tmRatio, 0, 20]]]] + tmSpec["Eigenvalues"];
tmCross[e_] := ArcCos[-e/tmRatio];
tmGHz[x_] := Row[{NumberForm[tmEC x, 3], " GHz"}];
tmLevelPlot = Show[Plot[-tmRatio Cos[φ], {φ, -2.2, 2.2}, PlotStyle -> GrayLevel[0.6]],
   Graphics[{MapIndexed[{ColorData[97][1], AbsoluteThickness[2], Line[{{-tmCross[#1], #1}, {tmCross[#1], #1}}],
        GrayLevel[0.5], Text[Style[Row[{"|", First[#2] - 1, "\[RightAngleBracket]"}], 12], {tmCross[#1] + 0.06, #1}, {-1, 0}]} &, tmLevels],
     Arrowheads[0.04], ColorData[97][2], Arrow[{{-0.12, tmLevels[[1]]}, {-0.12, tmLevels[[2]]}}], ColorData[97][4], Arrow[{{0.12, tmLevels[[2]]}, {0.12, tmLevels[[3]]}}],
     GrayLevel[0.5], Text[Style[Row[{Subscript["\[Omega]", "01"], "/2\[Pi] = ", tmGHz[tmSpec["Frequency"]]}], 11], {-0.3, Mean[tmLevels[[;; 2]]]}, {1, 0}],
     Text[Style[Row[{Subscript["\[Omega]", "12"], "/2\[Pi] = ", tmGHz[tmSpec["Frequency"] + tmSpec["Anharmonicity"]]}], 11], {0.3, Mean[tmLevels[[2 ;; 3]]]}, {-1, 0}]}],
   PlotRange -> {{-2.2, 2.2}, {-52, tmLevels[[4]] + 6}}, Axes -> False, Frame -> True, FrameLabel -> {"phase \[CurlyPhi] across the junction", Row[{"energy / ", Subscript["E", "C"]}]},
   PlotLabel -> Row[{Subscript["E", "J"], "/", Subscript["E", "C"], " = ", tmRatio, ",  \[Alpha]/2\[Pi] = ", tmGHz[tmSpec["Anharmonicity"]]}], ImageSize -> 250];
Grid[{{TransmonCircuitPlot["Circuit", Background -> White, ImageSize -> 550], tmLevelPlot}}, Spacings -> 2, Alignment -> Center]
```

The model below keeps only the lowest levels of this spectrum, through the qubit frequency $\omega_q$ and the anharmonicity $\alpha$, which is close to $-E_C$. With $\hbar=1$, the driven transmon and its form in the frame rotating at the drive frequency $\omega_d=\omega_q-\delta$, after the rotating-wave approximation, are

$$
\hat H=\omega_q \hat b^\dagger \hat b+\frac{\alpha}{2}\hat b^{\dagger2}\hat b^2+\left[\Omega(t)\cos\omega_d t-\frac{\lambda\dot\Omega(t)}{\alpha}\sin\omega_d t\right](\hat b+\hat b^\dagger)\;\longrightarrow\;\frac{\alpha}{2}\hat b^{\dagger2}\hat b^2+\delta\,\hat b^\dagger \hat b+\frac{\Omega(t)}{2}(\hat b+\hat b^\dagger)-\frac{\lambda\dot\Omega(t)}{2\alpha}\,i(\hat b^\dagger-\hat b),
$$

with the raised-cosine envelope $\Omega(t)=\frac{\pi}{\tau}\left(1-\cos\frac{2\pi t}{\tau}\right)$, whose area is $\pi$, and the oscillator kept to its five lowest levels.

Starting from $|0\rangle$, the figure plots the leakage $P(n\ge2)$ (solid) and the population left in $|0\rangle$, $P(0)$ (dashed), during the pulse, whose envelope is drawn in gray, and lists the leakage left at the end; times and frequencies are in units of $|\alpha|$:

```wl
dragB = AnnihilationOperator[5]; dragN = SuperDagger[dragB] @ dragB;
dragEnvelope[τ_] := Pi/τ (1 - Cos[2 Pi t/τ]);
dragH[α_, τ_, λ_, δ_] := α/2 (dragN @ dragN - dragN) + δ dragN + dragEnvelope[τ]/2 (dragB + SuperDagger[dragB]) -
    λ D[dragEnvelope[τ], t]/(2 α) I (SuperDagger[dragB] - dragB);
dragB2 = AnnihilationOperator[2];
dragPops[τ_, λ_, δ_] := QuantumEvolve[dragH[-1, τ, λ, δ], FockState[0, 5], {t, 0, τ}][t]["ProbabilitiesList"];
dragλs = {0, 0.5, 1};
dragPercent[x_] := If[x >= 10^-4, Row[{NumberForm[100 x, 2], "%"}], Row[{NumberForm[ScientificForm[100 x, 1]], "%"}]];
dragPanel[τ_, δ_] := With[{pops = dragPops[τ, #, δ] & /@ dragλs},
   Show[LogPlot[Evaluate[Total[#[[3 ;;]]] & /@ pops], {t, 0, τ}, PlotRange -> All, PlotLegends -> (Row[{"λ = ", #}] & /@ dragλs)],
    LogPlot[Evaluate[First /@ pops], {t, 0, τ}, PlotRange -> All, PlotStyle -> Dashed],
    PlotRange -> {Log[10^-9], Log[2]}, Frame -> True, GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.5, 0.35], Thin],
    Prolog -> With[{pts = Table[{s, Log[10^-9] - Log[10^-9] (1 - Cos[2 Pi s/τ])/2}, {s, 0, τ, τ/200}]},
      {{Opacity[0.12], GrayLevel[0.5], Polygon[Join[{{0, Log[10^-9]}}, pts, {{τ, Log[10^-9]}}]]}, {GrayLevel[0.5, 0.6], AbsoluteThickness[0.8], Line[pts]}}],
    FrameLabel -> {"|α| t", "population"},
    Epilog -> Inset[Style[Column[Prepend[MapThread[Style[Row[{"λ = ", #1, ":  ", dragPercent[Total[#2[[3 ;;]]] /. t -> τ]}], #3] &,
         {dragλs, pops, ColorData[97] /@ Range[3]}], "leaked at the end"], Alignment -> Center, Spacings -> 0.2], 10, FontFamily -> "Helvetica"],
      Scaled[{0.5, 0.04}], {Center, Bottom}], PlotLabel -> Row[{"δ/|α| = ", δ}], ImageSize -> 300]];
dragPanels = AssociationThread[{5, 10, 20}, Outer[dragPanel, {5, 10, 20}, {-0.05, 0, 0.05}]];
With[{legendBelow = legendBelow, controlNotes = controlNotes, dragPanels = dragPanels},
 Manipulate[Legended[PlotGrid[{First /@ dragPanels[τv]}, ImageSize -> 700], Placed[Row[{First[legendBelow[Last[dragPanels[5][[1]]]]], Spacer[30],
     LineLegend[{Directive[Gray], Directive[Gray, Dashed], Directive[GrayLevel[0.5, 0.6], AbsoluteThickness[0.8]]}, {"leakage P(n \[GreaterEqual] 2)", "left in |0\[RightAngleBracket]", "pulse envelope \[CapitalOmega](t) (shape)"}, LegendLayout -> {"Row", 1}]}], Below]], Row[{Control[{{τv, 10, Column[{"pulse length |α|τ", Style["α: anharmonicity", 10]}, Alignment -> Right, Spacings -> 0]}, Keys[dragPanels], SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {"", <|5 -> "Shortest pulse: without DRAG, leakage grows to the end.",
        10 -> "Leakage peaks and partly returns; the two shifts nearly cancel at λ = 1/2.",
        20 -> "The lowest leakage; the shifts nearly cancel at λ = 1/2."|>[τv]}}]]}]]]
```

The drive couples only neighboring levels, so leakage builds up through $|1\rangle$ during the pulse; only the shortest pulse without DRAG keeps leaking to the end, and otherwise it peaks and partly returns. At $\lambda=1$ the population left outside the qubit drops by orders of magnitude, and a longer pulse lowers it further. Most of the remaining error is population left in $|0\rangle$: a resonant π pulse on a two-level qubit would empty it, so what stays there mostly reflects how far the qubit transition has moved. The qubit transition is shifted by the $|2\rangle$ level and, in the opposite direction, by the DRAG quadrature itself; for the longer pulses the two shifts nearly cancel at $\lambda=1/2$, and a detuning partly undoes what is left, a negative one without DRAG and a positive one at $\lambda=1$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Entanglement held by a pump on an IBM processor: predicted, then measured

Two qubits are in the cluster state when both stabilizers $K_1=X_1Z_2$ and $K_2=Z_1X_2$ equal $+1$. A third qubit, the ancilla, can drive them there without anyone reading it: it picks up the value of one stabilizer, flips the pair back into the $+1$ branch when the value is $-1$, and is reset unread. Each reset carries entropy out of the pair, which works against the noise. How much entanglement does it hold on a real processor, and can a model built from the processor's own calibration say so before the run?

$$
K_1K_2=Y_1Y_2,\qquad F=\langle\psi|\rho|\psi\rangle=\frac{1}{4}\left(1+\langle X_1Z_2\rangle+\langle Z_1X_2\rangle+\langle Y_1Y_2\rangle\right),
$$

and every separable state has $F\le\frac{1}{2}$, so a fidelity above $\frac{1}{2}$ certifies entanglement. Each expectation value is a parity read from counts after rotating the pair into the matching basis.

One pump cycle (left), and the three jobs sent to the processor (right): the same three cycles, each half drawn as one block, then local rotations $R_{n_1}(\theta_1)$ and $R_{n_2}(\theta_2)$ into a basis $B$ that turns one stabilizer into a plain two-bit readout; the key gives each job's rotations.

```wl
pumpReset = QuantumOperator[QuantumState["0"], {{2}, {}}];
pumpHalves = {{"H" -> 2, "CNOT" -> {2, 1}, "CZ" -> {2, 3}, "H" -> 2, "CZ" -> {2, 1}, pumpReset}, {"H" -> 2, "CZ" -> {2, 1}, "CNOT" -> {2, 3}, "H" -> 2, "CZ" -> {2, 3}, pumpReset}};
pumpCycle = Catenate[pumpHalves];
pumpBasis = <|"X1Z2" -> {"H" -> 1}, "Z1X2" -> {"H" -> 3}, "Y1Y2" -> {"RZ"[-Pi/2] -> 1, "H" -> 1, "RZ"[-Pi/2] -> 3, "H" -> 3}|>;
pumpLabels = {Style["K\:2081 = X\:2081Z\:2082", 13], Style["K\:2082 = Z\:2081X\:2082", 13]};
pumpBox[label_, wires_] := QuantumOperator[IdentityMatrix[2^Length[wires]], wires, "Label" -> label];
pumpBasisBox[k_] := Row[{Subscript["R", Subscript["n", k]], "(", Subscript["\[Theta]", k], ")"}];
pumpJobLabels = {Style["K\:2081", 13], Style["K\:2082", 13]};
pumpJobsShown = QuantumCircuitOperator[Join[Flatten[Table[MapThread[QuantumCircuitOperator, {If[c == 3, {pumpHalves[[1]], Most[pumpHalves[[2]]]}, pumpHalves], pumpJobLabels}], {c, 3}]],
    {pumpBox[pumpBasisBox[1], {1}], pumpBox[pumpBasisBox[2], {3}], "M" -> 1, "M" -> 3}]];
pumpJobKey = Style[Row[{"job 1:  B = X\:2081Z\:2082  (H on qubit 1)", "job 2:  B = Z\:2081X\:2082  (H on qubit 2)",
     Row[{"job 3:  B = Y\:2081Y\:2082  (", Subscript["R", "Z"], "(\[Minus]\[Pi]/2) then H on both)"}]}, Spacer[30]], 13, FontFamily -> "Helvetica"];
Column[{Grid[{{Show[QuantumCircuitOperator[MapThread[QuantumCircuitOperator, {pumpHalves, pumpLabels}]]["Diagram", "WireLabels" -> {"1", "ancilla", "2"}],
     PlotRangeClipping -> False, ImagePadding -> {{50, 10}, {10, 40}}, ImageSize -> 350],
    pumpJobsShown["Diagram", "WireLabels" -> {"1", "ancilla", "2"}, "Expand" -> 0, "RotateGateLabel" -> False, ImageSize -> 450]}}, Spacings -> 3, Alignment -> Center], pumpJobKey}, Spacings -> 2, Alignment -> Center]
```

The pump runs three cycles on ibm_fez qubits 123, 122 and 121, with the ancilla 122 in the middle; among the ten best-calibrated lines of three neighbors, this one had the highest predicted fidelity in the calibration of 11 November 2026. The three circuits go to ibm_fez for 4096 shots each, on the chosen qubits; this needs an IBM Quantum account, with its API key in place of "xxx":

```wl
#| eval: false
ServiceConnect["IBMQuantumPlatform", "New", Authentication -> {"apikey" -> "xxx"}];
pumpJobs = IBMJobSubmit[QuantumCircuitOperator[Join[Most[Catenate[ConstantArray[pumpCycle, 3]]], #, {"M" -> 1, "M" -> 3}]], "ibm_fez",
     "Shots" -> 4096, "InitialLayout" -> {123, 122, 121}, "OptimizationLevel" -> 1] & /@ pumpBasis
```

The prediction runs the same circuits in QF as IBM compiles them, gate by gate, with QF channels for the noise the calibration lists: each gate's error as depolarizing noise, amplitude and phase damping from $T_1$ and $T_2$ while a qubit waits (the pair waits through every reset of the ancilla), a reset that leaves the ancilla excited with probability equal to its readout error, and each qubit's readout error on the counts. It was recorded before the jobs were submitted. One cycle as the model simulates it, with the channels named and their calibrated strengths left out:

```wl
pumpIdle[r_, t_, w_] := If[t <= 0, {}, With[{g = 1 - Exp[-t/r["T1"][w]], c = Min[1, Exp[-t/r["T2"][w] + t/(2 r["T1"][w])]]},
    {QuantumChannel["AmplitudeDamping"[g] -> {w}], QuantumChannel["PhaseDamping"[1 - c^2] -> {w}]}]];
pumpPaulis2 = Catenate[Outer[KroneckerProduct, #, #, 1]] &[{IdentityMatrix[2], PauliMatrix[1], PauliMatrix[2], PauliMatrix[3]}];
pumpDep2[lam_, ws_] := QuantumChannel[Prepend[Sqrt[lam/16] Rest[pumpPaulis2], Sqrt[1 - 15 lam/16] IdentityMatrix[4]], ws];
pumpGate[r_, {"rz", {t_}, {w_}}] := {QuantumOperator["RZ"[t], {w}]};
pumpGate[r_, {"sx", _, {w_}}] := {QuantumOperator["V", {w}], QuantumChannel["Depolarizing"[2 r["SX"][w]] -> {w}]};
pumpGate[r_, {"x", _, {w_}}] := {QuantumOperator["X", {w}], QuantumChannel["Depolarizing"[2 r["SX"][w]] -> {w}]};
pumpGate[r_, {"cz", _, ws_}] := {QuantumOperator["CZ", ws], pumpDep2[4/3 r["CZ"][Sort[ws]], ws]};
pumpGate[r_, {"reset", _, {w_}}] := {QuantumChannel["ResetError"[1 - r["Readout"][w], r["Readout"][w]] -> {w}]};
pumpDuration[r_, {"cz", _, ws_}] := r["CZLength"][Sort[ws]];
pumpDuration[r_, {"sx" | "x", _, {w_}}] := r["SXLength"][w];
pumpDuration[r_, {"reset", _, _}] := r["ResetLength"];
pumpDuration[_, _] := 0.;
pumpStep[r_][{els_, clock_}, i : {name_, _, ws_}] := If[MemberQ[{"measure", "barrier"}, name], {els, clock},
   With[{start = Max[clock /@ ws]}, {Join[els, Catenate[pumpIdle[r, start - clock[#], #] & /@ ws], pumpGate[r, i]],
     Join[clock, AssociationThread[ws, start + pumpDuration[r, i]]]}]];
pumpElements[r_, ins_] := With[{s = Fold[pumpStep[r], {{}, <|1 -> 0., 2 -> 0., 3 -> 0.|>}, ins]},
   With[{tEnd = Max[s[[2]] /@ {1, 3}]}, Join[s[[1]], Catenate[pumpIdle[r, tEnd - s[[2]][#], #] & /@ {1, 3}]]]];
pumpH[w_] := {{"rz", {Pi/2}, {w}}, {"sx", {}, {w}}, {"rz", {Pi/2}, {w}}};
pumpCNOT[c_, t_] := Join[pumpH[t], {{"cz", {}, {c, t}}}, pumpH[t]];
pumpNative[b_] := Join[Most[Catenate[ConstantArray[Join[pumpH[2], pumpCNOT[2, 1], {{"cz", {}, {2, 3}}}, pumpH[2], {{"cz", {}, {2, 1}}}, {{"reset", {}, {2}}},
       pumpH[2], {{"cz", {}, {2, 1}}}, pumpCNOT[2, 3], pumpH[2], {{"cz", {}, {2, 3}}}, {{"reset", {}, {2}}}], 3]]],
   <|"X1Z2" -> pumpH[1], "Z1X2" -> pumpH[3], "Y1Y2" -> Join[{{"rz", {-Pi/2}, {1}}}, pumpH[1], {{"rz", {-Pi/2}, {3}}}, pumpH[3]]|>[b]];
pumpShapeOnly = <|"T1" -> <|1 -> 100., 2 -> 100., 3 -> 100.|>, "T2" -> <|1 -> 80., 2 -> 80., 3 -> 80.|>, "SX" -> <|1 -> 0.001, 2 -> 0.001, 3 -> 0.001|>,
   "CZ" -> <|{1, 2} -> 0.005, {2, 3} -> 0.005|>, "Readout" -> <|1 -> 0.01, 2 -> 0.01, 3 -> 0.01|>,
   "SXLength" -> <|1 -> 0.03, 2 -> 0.03, 3 -> 0.03|>, "CZLength" -> <|{1, 2} -> 0.07, {2, 3} -> 0.07|>, "ResetLength" -> 1.5|>;
pumpCycleNative = Join[pumpH[2], pumpCNOT[2, 1], {{"cz", {}, {2, 3}}}, pumpH[2], {{"cz", {}, {2, 1}}}, {{"reset", {}, {2}}},
   pumpH[2], {{"cz", {}, {2, 1}}}, pumpCNOT[2, 3], pumpH[2], {{"cz", {}, {2, 3}}}, {{"reset", {}, {2}}}];
pumpModelEls = pumpElements[pumpShapeOnly, pumpCycleNative];
pumpTwoLine[a_, b_] := Column[{a, b}, Alignment -> Center, Spacings -> 0];
pumpChannelName[e_] := Which[MatchQ[e["Label"], "AmplitudeDamping"[___]], pumpTwoLine["Amplitude", "Damping"],
   MatchQ[e["Label"], "PhaseDamping"[___]], pumpTwoLine["Phase", "Damping"], MatchQ[e["Label"], "ResetError"[___]], pumpTwoLine["Reset", "Error"], True, "Depolarizing"];
pumpNamed = Which[MatchQ[#, _QuantumChannel], QuantumChannel[#, "Label" -> pumpChannelName[#]],
    MatchQ[#, _QuantumOperator] && #["Label"] === "V", QuantumOperator[#, "Label" -> Sqrt["X"]], True, #] & /@ pumpModelEls;
pumpHalfEnd = First[FirstPosition[pumpModelEls, e_ /; MatchQ[e, _QuantumChannel] && MatchQ[e["Label"], "ResetError"[___]], Missing[], {1}, Heads -> False]];
Show[QuantumCircuitOperator[MapThread[QuantumCircuitOperator, {{Take[pumpNamed, pumpHalfEnd], Drop[pumpNamed, pumpHalfEnd]}, pumpLabels}]]["Diagram",
   "WireLabels" -> {"1", "ancilla", "2"}, "RotateGateLabel" -> False], PlotRangeClipping -> False, ImagePadding -> {{50, 10}, {10, 40}}, ImageSize -> 1700]
```

Every step from the calibration to the counts, with the IBM queries written out, is in [pinned-pump-chain.md](/Users/mohammadb/Documents/GitHub/mbahram/qf-paclet-studies/IBM%20Heron%20Device/backbone/nov-pinned/pinned-pump-chain.md). The results saved from one such run (jobs dbqfceujjaos73ah7meg, dbqfcf8k09pc73at6h20 and dbqfcfok09pc73at6h30, 11 November 2026) are decoded by `IBMJob` and set against the exact cluster state and the recorded prediction:

```wl
pumpDir = "/Users/mohammadb/Documents/GitHub/mbahram/qf-paclet-studies/IBM Heron Device/backbone/nov-pinned";
pumpCounts = KeySort[IBMJob[<|"ID" -> "saved", "Status" -> "Completed", "Primitive" -> "sampler", "MeasuredQubits" -> {0, 1},
       "Raw" -> <|"Results" -> #|>|>]["Counts"]] & /@ Get[FileNameJoin[{pumpDir, "pinned_raw_results.m"}]]["raw"];
pumpPredicted = Get[FileNameJoin[{pumpDir, "prediction_recorded.m"}]]["Prediction"]["Parities"];
pumpCluster = QuantumState["Graph"[Graph[{1 <-> 2}]]];
pumpExact = AssociationThread[Keys[pumpCounts], Re[pumpCluster["Dagger"][QuantumOperator[#][pumpCluster]]["Scalar"]] & /@ {"XZ", "ZX", "YY"}];
pumpParity = Function[c, N[Total[KeyValueMap[(-1)^Total[#1] #2 &, c]]/Total[c]]];
pumpMeasured = pumpParity /@ pumpCounts;
pumpF = Function[p, (1 + Total[p])/4];
pumpFError = Sqrt[Total[(1 - Values[pumpMeasured]^2)/(Total /@ Values[pumpCounts])]]/4;
pumpColors = ColorData[97] /@ {1, 2, 3};
pumpLabelStyle = Sequence[13, FontFamily -> "Helvetica"];
Column[{Grid[{{
    BarChart[Transpose[Values[{pumpExact, pumpPredicted, pumpMeasured}]], ChartStyle -> pumpColors, ChartLabels -> {{"X\:2081Z\:2082", "Z\:2081X\:2082", "Y\:2081Y\:2082"}, None}, ImagePadding -> {{55, 10}, {30, 30}},
     Frame -> True, FrameLabel -> {None, "parity"}, PlotRange -> {0.7, 1.02}, PlotRangeClipping -> True, AspectRatio -> Full, ImageSize -> {440, 280},
     PlotLabel -> Style["the three stabilizers", pumpLabelStyle]],
    BarChart[{pumpF[pumpExact], pumpF[pumpPredicted], pumpF[pumpMeasured]}, ChartStyle -> pumpColors, Epilog -> {Directive[GrayLevel[0.5], Dashed, AbsoluteThickness[1.5]], Line[{{0.4, 1/2}, {3.6, 1/2}}]}, ImagePadding -> {{55, 10}, {30, 30}}, Frame -> True, FrameLabel -> {None, "fidelity F"}, FrameTicks -> {{Automatic, None}, {None, None}},
     PlotRange -> {{0.4, 3.6}, {0, 1.05}}, AspectRatio -> Full, ImageSize -> {240, 280},
     PlotLabel -> Style[Row[{"ibm_fez: F = ", NumberForm[pumpF[pumpMeasured], {3, 3}], " \[PlusMinus] ", NumberForm[pumpFError, 1]}], pumpLabelStyle]]}},
   Spacings -> 4, Alignment -> {Center, Top}],
  Row[{SwatchLegend[pumpColors, {"exact cluster state", "predicted before the run", "ibm_fez"}, LegendLayout -> "Row", LabelStyle -> {pumpLabelStyle}],
    Spacer[30], LineLegend[{Directive[GrayLevel[0.5], Dashed, AbsoluteThickness[1.5]]}, {"separable states: F \[LessEqual] 1/2"}, LabelStyle -> {pumpLabelStyle}]}]}, Spacings -> 1.5, Alignment -> Center]
```

Prediction and measurement both sit far above the line $F=\frac{1}{2}$, so the pair ends entangled after three pumped rounds of resets that nobody reads. The model gets the order of the three parities right: $Z_1X_2$, pumped last in each cycle, holds best, while $X_1Z_2$ and $Y_1Y_2$ lose more. It is optimistic by far more than the shot noise, so the run loses something the calibration-based model does not capture.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The first Bianchi identity from the Riemann tensor's Young diagram

The Riemann tensor $R_{abcd}$ in $d$ dimensions has the pair symmetries $R_{abcd}=-R_{bacd}=-R_{abdc}=R_{cdab}$ and the first Bianchi identity $R_{abcd}+R_{acdb}+R_{adbc}=0$. Suppose only that $R$ lies in the image of the Young projector of shape $\{2,2\}$, which is one irreducible representation of $GL(d)$ inside the $\{2,2\}$ sector of $S_4$ acting on the four slots. Does that projector alone impose the first Bianchi identity?

For the tableau with rows $\{1,3\},\{2,4\}$ the projector is

$$
P_{\{2,2\}}\,T=\frac{2}{4!}\;\mathcal A_{\{1,2\}}\mathcal A_{\{3,4\}}\,\mathcal S_{\{1,3\}}\mathcal S_{\{2,4\}}\,T,
$$

where $\mathcal S$ sums over the permutations of a row of slots, $\mathcal A$ sums over the permutations of a column with their signs, and $2$ is the dimension of the $S_4$ representation $\{2,2\}$.

Each relation is a linear map on four-index tensors, and every tensor in the projector's image obeys it exactly when the map times the projector is zero. On four-index tensors in $d=4$ the 24 slot permutations act faithfully (Schur–Weyl duality), so an identity between such maps that holds there holds in the group algebra of the permutations, and so in every dimension. The cell builds the projector and each relation as exact $256\times256$ matrices and multiplies them, with a control relation that is not a Riemann symmetry. It then compares, for symbolic $d$, the projector's rank, the dimension $s_{\{2,2\}}(d)$ of its $GL(d)$ representation from TensorNetworks' `SchurDimension`, with the number of components the pair symmetries and Bianchi leave: the pair symmetries make $R$ a symmetric matrix on the $N=d(d-1)/2$ antisymmetric pairs, and Bianchi then only sets the totally antisymmetric part $R_{[abcd]}$ to zero, which removes $\binom d4$ more:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`Symmetry`"];
riemannTableau = YoungTableau[{{1, 3}, {2, 4}}];
(* every operator as an exact 256 x 256 matrix on 4-index tensors in d = 4 *)
slotOperator[f_] := Transpose[Flatten[f[ArrayReshape[#, {4, 4, 4, 4}]]] & /@ IdentityMatrix[256]];
projectorM = slotOperator[YoungProject[#, riemannTableau] &];
relationM = <|
   Row[{Subscript["R", "abcd"], " + ", Subscript["R", "acdb"], " + ", Subscript["R", "adbc"]}] -> slotOperator[Function[T, Table[T[[a, b, c, e]] + T[[a, c, e, b]] + T[[a, e, b, c]], {a, 4}, {b, 4}, {c, 4}, {e, 4}]]],
   Row[{Subscript["R", "abcd"], " + ", Subscript["R", "bacd"]}] -> slotOperator[Function[T, T + Transpose[T, {2, 1, 3, 4}]]],
   Row[{Subscript["R", "abcd"], " + ", Subscript["R", "abdc"]}] -> slotOperator[Function[T, T + Transpose[T, {1, 2, 4, 3}]]],
   Row[{Subscript["R", "abcd"], " \[Minus] ", Subscript["R", "cdab"]}] -> slotOperator[Function[T, T - Transpose[T, {3, 4, 1, 2}]]]|>;
relationHolds = (Max[Abs[# . projectorM]] == 0) & /@ relationM;
controlM = slotOperator[Function[T, T + Transpose[T, {1, 3, 2, 4}]]];
controlNonzero = Count[Flatten[controlM . projectorM], x_ /; x != 0];
dimRiemann = With[{np = d (d - 1)/2}, np (np + 1)/2] - Binomial[d, 4];
dimImage = SchurDimension[{2, 2}, d];
Grid[Join[
   {{Item[Style["on every tensor in the image of the projector, exactly:", Bold], Alignment -> Left], SpanFromLeft}},
   KeyValueMap[{Row[{#1, " = 0"}], #2} &, relationHolds],
   {{Row[{"control:  ", Subscript["R", "abcd"], " + ", Subscript["R", "acbd"], " = 0"}], Row[{"False  (", controlNonzero, " of the 65536 entries are nonzero)"}]}},
   {{Item[Style["dimensions, for symbolic d:", Bold], Alignment -> Left], SpanFromLeft},
    {Row[{"rank of ", Subscript["P", "{2,2}"], " = ", Subscript["s", "{2,2}"], "(d)"}], TraditionalForm[Factor[dimImage]]},
    {Row[{"pair symmetries and Bianchi:  N(N+1)/2 \[Minus] ", TraditionalForm[HoldForm[Binomial[d, 4]]]}], TraditionalForm[Factor[FunctionExpand[dimRiemann]]]},
    {"equal for every d", Simplify[FunctionExpand[dimRiemann - dimImage] == 0]}}],
  Alignment -> {{Right, Left}, Center}, Spacings -> {2, 0.9}, BaseStyle -> {13, FontFamily -> "Helvetica"}]
```

The two counts and their difference against $d$, the closed forms as lines and, as dots, the counts for $d=2,\dots,6$ taken directly from `SymmetrizedIndependentComponents` and the projector's rank:

```wl
pairCount[dd_] := Length[SymmetrizedIndependentComponents[{dd, dd, dd, dd}, {{Cycles[{{1, 2}}], -1}, {Cycles[{{3, 4}}], -1}, {Cycles[{{1, 3}, {2, 4}}], 1}}]];
projectorRank[dd_] := MatrixRank[Flatten[YoungProject[ArrayReshape[#, {dd, dd, dd, dd}], riemannTableau]] & /@ IdentityMatrix[dd^4]];
countTable = Table[{dd, pairCount[dd], projectorRank[dd]}, {dd, 2, 6}];
Show[Plot[Evaluate[{dimRiemann + Binomial[d, 4], dimImage, Binomial[d, 4]} /. d -> x], {x, 2, 6.2}, PlotStyle -> (Directive[ColorData[97][#], Opacity[0.6]] & /@ {1, 2, 3}),
   PlotLegends -> Placed[{"pair symmetries only, N(N+1)/2", Row[{"rank of ", Subscript["P", "{2,2}"], ", the Riemann count"}], Row[{"difference, ", TraditionalForm[HoldForm[Binomial[d, 4]]]}]}, {Left, Top}]],
 ListPlot[{countTable[[All, {1, 2}]], countTable[[All, {1, 3}]], {#1, #2 - #3} & @@@ countTable}, PlotMarkers -> {Automatic, 9}, PlotStyle -> (ColorData[97] /@ {1, 2, 3})],
 PlotRange -> {{1.8, 6.3}, {-5, 125}}, Frame -> True, FrameTicks -> {{Automatic, None}, {Range[2, 6], None}}, Axes -> False,
 FrameLabel -> {"dimension d", "independent components"}, GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.5, 0.35], Thin],
 PlotLabel -> "dots: counted directly;  lines: closed forms", ImageSize -> 420]
```

The projector alone makes every tensor in its image antisymmetric in each pair, symmetric under exchanging the pairs, and obey the first Bianchi identity. Its rank equals the number of components those relations leave, for every $d$, so it imposes nothing else: the Riemann symmetries are exactly the statement that $R$ lies in the image of this tableau's Young projector of shape $\{2,2\}$, with Bianchi included rather than added. The identity removes no component for $d\le3$, one in four dimensions and $\binom d4$ in general.

