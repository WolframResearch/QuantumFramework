<!-- Build: MarkdownToNotebook["EU-Tour-2026-Showcase.md", "EU-Tour-2026-Showcase.nb", "EvaluateSeparator" -> None]. All cells share one evaluation context; the loading cell runs first. Example bodies are taken verbatim from Modeling-Slides.md, Catalogue sources and the TensorNetworks documentation. -->

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

# Quantum in the Wolfram Language

#| style: Subtitle

What QuantumFramework, TensorNetworks and the Wolfram Language compute, one problem at a time

#| style: Subsubtitle

Mads Bahrami, Wolfram Research. Europe tour, November 2026. quantum@wolfram.com

The loading cell uses the working copies of QuantumFramework and TensorNetworks; examples that need QECCore, TransmonLab, QUALink or the DiracMatrix package load them where they are used.

Every example below is a physics question, the model, one runnable cell and what the result says. The chapters follow capabilities.

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"];
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
legendBelow[\[FormalL]_] := Placed[\[FormalL] /. Placed[\[FormalX]_, ___] :> \[FormalX] /. Column[\[FormalY]_List, ___] :> Row[\[FormalY], Spacer[20]] /.
    (\[FormalH] : LineLegend | PointLegend | SwatchLegend)[\[FormalA]__] :> \[FormalH][Sequence @@ DeleteCases[{\[FormalA]}, LegendLayout -> _], LegendLayout -> {"Row", 1}], Below];
controlNotes[\[FormalR]_List] := Grid[MapAt[Style[#, Bold] &, \[FormalR], {All, 1}], Alignment -> {Left, Top}, ItemSize -> {{Automatic, 40}},
    Spacings -> {1, 0.6}, BaseStyle -> {12, LineIndent -> 0}];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Exact dynamics with symbolic Hamiltonians

A `QuantumOperator` can hold symbols, so `QuantumEvolve` returns propagators as formulas in the parameters rather than tables of numbers. These examples solve time-dependent problems exactly and read the physics off the closed form.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A qubit driven by a rotating field


A qubit with splitting $\omega_0$ is driven by a field of strength $\Omega$ rotating at frequency $\omega$, detuned by $\Delta=\omega_0-\omega$. How much population does the drive move out of $|0\rangle$, and how fast?

With $\hbar=1$, the lab-frame Hamiltonian is

$$
\hat H(t)=\frac{\omega_0}{2}\,\hat\sigma_z+\frac{\Omega}{2}\left(\cos\omega t\,\hat\sigma_x+\sin\omega t\,\hat\sigma_y\right),
$$

which does not commute with itself at different times; the propagator solves $i\,\dot U=\hat H(t)\,U$ with $U(0)=\mathbb 1$.

The propagator for general parameters gives $P(1)=|\langle1|U(t)|0\rangle|^2$:

```wl
rabiP = FullSimplify[ComplexExpand[Abs[Normal[QuantumEvolve[\[Omega]0/2 QuantumOperator["Z"] + \[CapitalOmega]/2 (Cos[\[Omega] t] QuantumOperator["X"] + Sin[\[Omega] t] QuantumOperator["Y"]), None, t]["Matrix"]][[2, 1]]]^2],
    {\[Omega]0, \[Omega], t} \[Element] Reals && \[CapitalOmega] > 0]
```

The plot takes three detunings, with time in units of $1/\Omega$:

```wl
rabiDetunings = {0, 1, 2};
Plot[Evaluate[rabiP /. {\[CapitalOmega] -> 1, \[Omega] -> \[Omega]0 - #} & /@ rabiDetunings], {t, 0, 4 Pi}, PlotRange -> {0, 1}, Frame -> True,
  FrameLabel -> {"\[CapitalOmega] t", "P(1)"}, ImageSize -> 420, PlotLegends -> (Row[{"\[CapitalDelta]/\[CapitalOmega] = ", #}] & /@ rabiDetunings)]
```

The drive moves at most $\Omega^2/(\Omega^2+\Delta^2)$ of the population, and it oscillates at the generalized Rabi frequency $\sqrt{\Omega^2+\Delta^2}$, faster than on resonance. Only the detuning enters: $\omega_0$ and $\omega$ separately drop out.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A Landau-Zener sweep


A qubit has its bias swept at rate $v$ through an avoided crossing with gap $\Delta$. It starts in $|0\rangle$ long before the crossing. How does it reach its final population?

With $\hbar=1$, the Hamiltonian is

$$
\hat H(t)=\frac{vt}{2}\,\hat\sigma_z+\frac{\Delta}{2}\,\hat\sigma_x,
$$

which does not commute with itself at different times; the propagator solves $i\,\dot U=\hat H(t)\,U$ with $U(0)=\mathbb 1$, and for an infinitely long sweep the probability to stay in $|0\rangle$ is the Landau-Zener value $e^{-\pi\Delta^2/(2v)}$.

The qubit starts in $|0\rangle$ at $\sqrt v\,t=-20$ under the propagator for general $v$, $\Delta$ and $t$; the figure plots $P(0)$ from $\sqrt v\,t=-10$, with the Landau-Zener value dashed:

```wl
lzU = Normal[QuantumEvolve[v t/2 QuantumOperator["Z"] + \[CapitalDelta]/2 QuantumOperator["X"], None, t]["Matrix"]];
lzStay[r_] := With[{u = N[First[lzU /. {v -> 1, \[CapitalDelta] -> Sqrt[r]}]]}, Abs[u . Conjugate[u /. t -> -20]]^2];
{lzRatios, lzTimes} = {{0.1, 0.6, 1.6}, Range[-10, 10, 0.05]};
Show[ListLinePlot[Table[Transpose[{lzTimes, lzStay[r] /. t -> lzTimes}], {r, lzRatios}], PlotRange -> All,
   PlotLegends -> (Row[{"\[CapitalDelta]\.b2/v = ", #}] & /@ lzRatios)],
  Plot[Evaluate[Exp[-Pi lzRatios/2]], {t, -10, 10}, PlotStyle -> Dashed],
  PlotRange -> {0, 1}, Frame -> True, Axes -> False, FrameLabel -> {Row[{Sqrt["v"], " t"}], "P(0)"}, ImageSize -> 420]
```

The qubit drops through the crossing at $t=0$ and then settles toward the Landau-Zener value through oscillations that fade as the levels separate; the asymptotic formula does not contain them. A slower sweep, larger $\Delta^2/v$, leaves less population in $|0\rangle$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A Gaussian pulse on a qubit


A qubit with splitting $\omega_0$ is driven at frequency $\omega$, detuned by $\Delta=\omega_0-\omega$, by a Gaussian pulse of peak Rabi frequency $\Omega_0$, width $\tau$, centered at $t_0$. How much population has the pulse moved into $|1\rangle$ by time $t$, on and off resonance?

With $\hbar=1$, the lab-frame Hamiltonian and its form in the frame rotating at $\omega$, after the rotating-wave approximation, are

$$
\hat H=\frac{\omega_0}{2}\,\hat\sigma_z+\Omega(t)\cos\omega t\,\hat\sigma_x\;\longrightarrow\;\frac{\Delta}{2}\,\hat\sigma_z+\frac{\Omega(t)}{2}\,\hat\sigma_x,\qquad \Omega(t)=\Omega_0\,e^{-(t-t_0)^2/2\tau^2},
$$

and the propagator solves $i\,\dot U=\hat H\,U$ with $U(0)=\mathbb 1$.

On resonance the whole $2\times2$ propagator is solved exactly for general parameters; off resonance it is solved numerically. Each panel plots $P(1)$ from $|0\rangle$ together with the pulse $\Omega(t)$, for pulses at $t_0=4\tau$ of total area $A=\sqrt{2\pi}\,\Omega_0\tau$; every $\Omega(t)$ is drawn in units of the $2\pi$ pulse's peak, and time is in units of $\tau$:

```wl
pulseEq[d_, w_] := {upul'[t] == -I (d/2 PauliMatrix[3] + w/2 Exp[-(t - t0)^2/(2 \[Tau]^2)] PauliMatrix[1]) . upul[t], upul[0] == IdentityMatrix[2]};
upulse = Normal[QuantumEvolve[\[CapitalOmega]0/2 Exp[-(t - t0)^2/(2 \[Tau]^2)] QuantumOperator["X"], None, t]["Matrix"]];
pulseNum[a_] := {\[CapitalOmega]0 -> a/Sqrt[2 Pi], \[Tau] -> 1, t0 -> 4};
pulseP1[0, a_] := Abs[upulse[[2, 1]]]^2 /. pulseNum[a];
pulseP1[d_, a_] := NDSolveValue[pulseEq[d, \[CapitalOmega]0] /. pulseNum[a], Abs[Indexed[upul[t], {2, 1}]]^2, {t, 0, 8}];
{pulseDetunings, pulseAreas} = {{0, 0.5, 2}, {{Pi/2, "\[Pi]/2"}, {Pi, "\[Pi]"}, {2 Pi, "2\[Pi]"}}};
pulsePanel[d_, {a_, label_}] := With[{p = pulseP1[d, a]},
    Plot[{p, a/(2 Pi) Exp[-(t - 4)^2/2]}, {t, 0, 8}, PlotRange -> {0, 1}, Frame -> True, FrameLabel -> {"t/\[Tau]", None}, ImageSize -> 300,
     PlotLabel -> Style[Row[{"A = ", label}], 13, FontFamily -> "Helvetica"]]];
panelsP = Outer[pulsePanel, pulseDetunings, pulseAreas, 1];
Manipulate[Legended[PlotGrid[{panelsP[[pulseK]]}, ImageSize -> 900], Placed[LineLegend[{ColorData[97][1], ColorData[97][2]}, {"P(1)", Row[{"\[CapitalOmega](t), in units of the 2\[Pi] pulse's peak"}]}, LegendLayout -> "Row"], Below]],
 Row[{Control[{{pulseK, 2, "\[Tau]\[CapitalDelta]"}, MapIndexed[First[#2] -> #1 &, pulseDetunings], SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"τΔ = ", pulseDetunings[[pulseK]]}], <|1 -> "Resonance: area π flips the qubit fully, 2π brings it back.",
        2 -> "The π pulse stops short of a full flip.",
        3 -> "Every pulse leaves almost nothing in |1⟩ at the end."|>[pulseK]}}]]}], SaveDefinitions -> True]
```

On resonance $\hat H$ points along $\hat\sigma_x$ at every instant, so the propagator is the rotation by the area swept so far and a pulse with $A=\pi$ is a full flip. Detuning tilts the rotation axis by an angle that changes along the pulse, so the same pulse moves less population: at $\tau\Delta=0.5$ the $\pi$ pulse stops short of a full flip, and at $\tau\Delta=2$ the pulses leave almost nothing in $|1\rangle$ at the end.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The Mandelstam-Tamm bound for any drive


A quantum system is driven by an arbitrary time-dependent Hamiltonian. Can its state ever move farther from where it started than its accumulated energy spread allows?

With $\hbar=1$ and a Hamiltonian $\hat H(t)$ of undetermined entries,

$$
\frac{d}{dt}|\psi\rangle=-i\hat H|\psi\rangle,\qquad \frac{d}{dt}\langle\psi|=i\langle\psi|\hat H,\qquad \langle\dot\psi|\dot\psi\rangle-|\langle\psi|\dot\psi\rangle|^2=\Delta E^2 .
$$

The third relation equates the squared Fubini-Study speed with the energy variance, which bounds the distance $\arccos|\langle\psi_0|\psi_t\rangle|$ by the path length $\int_0^t\Delta E\,dt'$, the Mandelstam-Tamm bound in the geometric form of Anandan and Aharonov. For the runs, $\hat H(t)=h_0\mathbb 1+h_x\hat\sigma_x+h_y\hat\sigma_y+h_z\hat\sigma_z$ on a qubit with random amplitudes and frequencies of order $\epsilon$, and the constant drive $\hat H=\epsilon\,\hat\sigma_x$.

The identity is checked on a three-level system with every entry of $\hat H$ and the state general, using the Schrödinger equation only as a rule for time derivatives; three random drives are then integrated from $|0\rangle$ in units of $1/\epsilon$, with the path length dashed:

```wl
ehH = Array[ehh[##][t] &, {3, 3}];
{ehKet, ehBra} = {Array[ehc[#][t] &, 3], Array[ehb[#][t] &, 3]};
ehRules = Join[Thread[D[ehKet, t] -> -I ehH . ehKet], Thread[D[ehBra, t] -> I ehBra . ehH]];
speedResidual = Simplify[((D[ehBra, t] . D[ehKet, t] - (ehBra . D[ehKet, t]) (D[ehBra, t] . ehKet)) /. ehRules) - (ehBra . ehH . ehH . ehKet - (ehBra . ehH . ehKet)^2)];
mtDrive[seed_] := Total[MapThread[(#1[[1]] + #1[[2]] Cos[3 #1[[3]] s + 2 Pi #1[[4]]]) #2 &,
    {BlockRandom[SeedRandom[seed]; RandomReal[{-1, 1}, {4, 4}]], PauliMatrix[Range[0, 3]]}]];
mtRun[hn_] := With[{psi = {ehC1[s], ehC2[s]}},
    NDSolveValue[{D[psi, s] == -I hn . psi, (psi /. s -> 0) == {1, 0},
      ehLen'[s] == Sqrt[Abs[Conjugate[psi] . hn . hn . psi - (Conjugate[psi] . hn . psi)^2]], ehLen[0] == 0},
      {ArcCos[Min[1, Abs[ehC1[s]]]], ehLen[s]}, {s, 0, 6}]];
mtRuns = mtRun /@ Append[mtDrive /@ {1, 2, 3}, PauliMatrix[1]];
Legended[Show[Plot[Evaluate[First /@ mtRuns], {s, 0, 6}, PlotRange -> All],
  Plot[Evaluate[Last /@ mtRuns], {s, 0, 6}, PlotStyle -> Dashed, PlotRange -> All],
  PlotRange -> {0, All}, Frame -> True, FrameLabel -> {"\[Epsilon] t", Row[{"arccos|\[LeftAngleBracket]", Subscript["\[Psi]", 0], "|", Subscript["\[Psi]", "t"], "\[RightAngleBracket]|"}]}, ImageSize -> 440,
  PlotLabel -> Style[Row[{"three levels:  speed\.b2 \[Minus] \[CapitalDelta]E\.b2 = ", speedResidual}], 13, FontFamily -> "Helvetica"]],
 Column[{LineLegend[ColorData[97] /@ Range[4], {"drive 1", "drive 2", "drive 3", Row[{"constant ", Subscript["\[Sigma]", "x"]}]}],
   LineLegend[{Directive[Gray, Dashed]}, {"path length \[Integral]\[CapitalDelta]E dt"}]}]]
```

The residual vanishes for every drive and state, and nothing in the algebra uses the dimension, so no drive can move the state farther than its path length. Along each random drive the distance stays below its dashed path length; it can turn back toward $|0\rangle$ while the path length only grows. A constant drive between two levels runs along the bound, distance equal to path length, until the state is orthogonal to where it started.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The Hamiltonian behind a chirped rotation


A qubit should follow the trajectory $U(t)=R_Z(\beta t^2)\,R_X(\Omega t)$: a rotation about $x$ at rate $\Omega$, followed by a $z$ rotation whose angle grows quadratically at chirp rate $\beta$. What control Hamiltonian drives it?

With $\hbar=1$, Pauli matrices $X,Y,Z$ and $R_Z(\vartheta)=e^{-i\vartheta Z/2}$, $R_X(\vartheta)=e^{-i\vartheta X/2}$, a propagator obeys $i\dot U=\hat HU$, so the Hamiltonian is $\hat H(t)=i\,\dot U(t)\,U^\dagger(t)$, written as $\hat H=h_x X+h_y Y+h_z Z$.

Differentiating $U(t)$ gives $i\dot UU^\dagger$, whose three Pauli coefficients are read off, and evolving under them is checked to return $U(t)$; units of $\Omega$, time in units of $1/\Omega$:

```wl
hfU = QuantumOperator["RZ"[β t^2]] @ QuantumOperator["RX"[Ω t]];
hfField = FullSimplify[ComplexExpand[#]] & /@ (I D[hfU, t] @ hfU["Dagger"])["PauliDecompose"];
hfBackQ = FullSimplify[Normal[QuantumEvolve[Total[KeyValueMap[#2 QuantumOperator[#1] &, hfField]], None, t]["Matrix"]] == Normal[hfU["Matrix"]], {β, Ω, t} ∈ Reals];
hfPanel[b_] := Plot[Evaluate[Lookup[hfField, {"X", "Y", "Z"}] /. {Ω -> 1, β -> b}], {t, 0, 6}, PlotRange -> All,
   PlotLegends -> {Subscript["h", "x"], Subscript["h", "y"], Subscript["h", "z"]}, Frame -> True, FrameLabel -> {"Ω t", "coefficient / Ω"},
   PlotLabel -> Row[{"β/", Superscript["Ω", 2], " = ", b}], ImageSize -> 300];
hfPanels = hfPanel /@ {0.05, 0.2, 0.5};
Labeled[Legended[PlotGrid[{First /@ hfPanels}, ImageSize -> 900], legendBelow[Last[First[hfPanels]]]],
 Style[Row[Append[Riffle[MapThread[Row[{Subscript["h", #1], " = ", TraditionalForm[#2]}] &, {{"x", "y", "z"}, Lookup[hfField, {"X", "Y", "Z"}]}], ",    "],
   Row[{";    evolving under H returns U(t): ", hfBackQ}]]], 13, FontFamily -> "Helvetica"], Top]
```

The $z$ field splits the levels by $2\beta t$, exactly the instantaneous frequency of the transverse field turning with phase $\beta t^2$, so the chirped drive stays resonant at every instant. In the frame that turns with it the qubit sees a constant field $\Omega/2$ along $x$, which is the rotation $R_X(\Omega t)$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A spin coherent state


A spin $j$ starts in $|j,j\rangle$, pointing along $+z$, and is rotated to point along the direction $\hat n$ with polar angle $\theta$ and azimuth $\varphi$. How is the rotated spin spread over the levels $m$ of $\hat J_z$?

With $\hat J_z|j,m\rangle=m|j,m\rangle$, the spin coherent state and its level populations are

$$
|j,\hat n\rangle=e^{-i\varphi\hat J_z}e^{-i\theta\hat J_y}|j,j\rangle,\qquad P(m)=|\langle j,m|j,\hat n\rangle|^2=\left|d^{\,j}_{mj}(\theta)\right|^2,
$$

where $d^{\,j}_{m'm}(\theta)=\langle j,m'|e^{-i\theta\hat J_y}|j,m\rangle$ is the Wigner small-$d$ matrix; $\varphi$ only multiplies each amplitude by a phase.

Applying $e^{-i\theta\hat J_y}$ to $|j,j\rangle$ for general $\theta$ gives $P(m)$, checked against the binomial law; the figure plots $P(m)$ against $m/j$, with the Gaussian of the same mean $j\cos\theta$ and variance $\frac j2\sin^2\theta$ dashed:

```wl
spinCohP[j_] := ComplexExpand[Abs[Normal[(QuantumOperator["WignerD"[j, -\[Theta]]] @ QuantumState[UnitVector[2 j + 1, 2 j + 1], 2 j + 1])["StateVector"]]]^2];
{spins, polars} = {{2, 10, 40}, {Pi/6, Pi/2, 5 Pi/6}};
spinCoh = AssociationMap[spinCohP, spins];
binomialQ = And @@ (Simplify[spinCoh[#] == Table[Binomial[2 #, # + m] Cos[\[Theta]/2]^(2 (# + m)) Sin[\[Theta]/2]^(2 (# - m)), {m, -#, #}], 0 < \[Theta] < Pi] & /@ spins);
spinPanel[j_] := Show[ListLinePlot[Table[Transpose[{Range[-j, j]/j, spinCoh[j] /. \[Theta] -> th}], {th, polars}], PlotMarkers -> Automatic, PlotRange -> All,
     PlotLegends -> (Row[{"\[Theta] = ", #}] & /@ polars)],
    Plot[Evaluate[PDF[NormalDistribution[j Cos[#], Sqrt[j/2] Sin[#]], j mj] & /@ polars], {mj, -1, 1}, PlotStyle -> Dashed, PlotRange -> All],
    PlotRange -> All, Axes -> False, Frame -> True, FrameLabel -> {"m/j", "P(m)"}, PlotLabel -> Row[{"j = ", j}], ImageSize -> 300];
spinPanels = spinPanel /@ spins;
Labeled[Legended[Grid[{First /@ spinPanels}], legendBelow[Last[First[spinPanels]]]],
  Style[Row[{"P(m) = ", TraditionalForm[HoldForm[Binomial[2 j, j + m] Cos[\[Theta]/2]^(2 (j + m)) Sin[\[Theta]/2]^(2 (j - m))]], " for 0 < \[Theta] < \[Pi] at j = 2, 10, 40: ", binomialQ}],
    13, FontFamily -> "Helvetica"], Top]
```

The populations are binomial: the spin-$j$ coherent state is $2j$ spin-$\frac12$ particles all pointing along $\hat n$, each found up with probability $\cos^2(\theta/2)$. As $j$ grows the binomial approaches the dashed Gaussian centered at $m=j\cos\theta$, and its width relative to $j$ shrinks as $1/\sqrt j$, so a large spin points like a classical arrow.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The spin-1 NV center in a tilted magnetic field


The ground state of a nitrogen-vacancy center in diamond is a spin $S=1$. The crystal field lifts $m_s=\pm1$ above $m_s=0$ by the zero-field splitting $D$, strain adds a transverse anisotropy $E$, and a magnetic field $B$, tilted by $\vartheta$ from the defect axis $z$, couples through the electron's gyromagnetic ratio $\gamma$. Which level crossings survive strain and a tilted field?

With $\hbar=1$ and $\hat S_x,\hat S_y,\hat S_z$ the spin-1 matrices,

$$
\hat H=D\,\hat S_z^2+E\left(\hat S_x^2-\hat S_y^2\right)+\gamma B\left(\cos\vartheta\,\hat S_z+\sin\vartheta\,\hat S_x\right).
$$

The untilted spectrum is in closed form for general $D$, $E$ and $\gamma B$; the figure plots the three energy levels against the field for $E/D=0,\,0.1$ and $\vartheta=0,\,5^\circ,\,30^\circ$, in units of $D$:

```wl
{nvSx, nvSy, nvSz} = QuantumOperator /@ {"JX"[1], "JY"[1], "JZ"[1]};
nvMatrix = Normal[(nvD nvSz @ nvSz + nvE (nvSx @ nvSx - nvSy @ nvSy) + nvB (Cos[nvTheta] nvSz + Sin[nvTheta] nvSx))["MatrixRepresentation"]];
nvAxial = Simplify[Eigenvalues[nvMatrix /. nvTheta -> 0]];
nvCross = nvB /. First[Solve[Equal @@ Take[nvAxial, 2] && nvB > 0 && nvD > nvE > 0, nvB]];
nvLevels[b_?NumericQ, e_, th_] := Sort[Eigenvalues[N[nvMatrix /. {nvD -> 1, nvB -> b, nvE -> e, nvTheta -> th}]]];
nvLevel[k_, b_?NumericQ, e_, th_] := nvLevels[b, e, th][[k]];
nvPanel[e_, th_] := Plot[Evaluate[Table[nvLevel[k, b, e, th], {k, 3}]], {b, 0, 2}, PlotRange -> All, Frame -> True, Axes -> False, FrameLabel -> {"\[Gamma]B / D", "energy / D"},
   PlotLabel -> Row[{"E/D = ", e, ",   \[CurlyTheta] = ", th}], ImageSize -> 260];
Labeled[Grid[Outer[nvPanel, {0, 0.1}, {0, 5 Degree, 30 Degree}]],
  Style[Row[{"\[CurlyTheta] = 0:  levels ", TraditionalForm[nvAxial /. {nvD -> "D", nvE -> "E", nvB -> HoldForm[\[Gamma] B]}],
    ",   crossing at \[Gamma]B = ", TraditionalForm[Normal[nvCross] /. {nvD -> "D", nvE -> "E"}]}], 13, FontFamily -> "Helvetica"], Top]
```

Strain couples only $m_s=+1$ to $m_s=-1$, so it opens their crossing at zero field into a gap $2E$ but leaves the crossing of the lower level with $m_s=0$, moved to a slightly lower field, intact. A tilted field mixes $m_s=0$ with $m_s=\pm1$ and turns the crossing into an avoided one, with a gap that widens with the tilt.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A five-site Heisenberg ring split into symmetry sectors


Five spins $\frac12$ sit on a ring with antiferromagnetic exchange $J>0$ between neighbors. The Hamiltonian commutes with the total magnetization $\hat M=\sum_j\hat S^z_j$ and with the translation $\hat T$ that moves every spin one site along the ring, whose eigenvalues are $e^{ik}$ with $k=2\pi m/5$, $m=-2,\dots,2$. Which magnetization and momentum does the ground state carry?

With $\hat{\mathbf S}_j=\hat{\boldsymbol\sigma}_j/2$, $|0\rangle$ as spin up, and site 6 identified with site 1,

$$
\hat H=J\sum_{j=1}^{5}\hat{\mathbf S}_j\cdot\hat{\mathbf S}_{j+1},\qquad [\hat H,\hat M]=[\hat H,\hat T]=[\hat M,\hat T]=0 ;
$$

flipping every spin commutes with $\hat H$ and sends $\hat M\to-\hat M$, so the sectors with $M<0$ repeat those with $M>0$.

On the 32-dimensional space each joint eigenspace of $\hat M$ and $\hat T$ gets an exact orthonormal basis, and $\hat H$ is diagonalized in each in closed form; the figure plots $|\hat H|$ in the product basis, $|\hat H|$ in the basis sorted by $(M,k)$, and the energies of each sector, in units of $J$:

```wl
ringPairs = Table[{j, Mod[j, 5] + 1}, {j, 5}];
ringH = Normal[(Total[QuantumOperator @@@ Tuples[{{"XX", "YY", "ZZ"}, ringPairs}]]/4)["Matrix"]];
ringT = Normal[QuantumOperator["Permutation"[ConstantArray[2, 5], Cycles[{Range[5]}]]]["Matrix"]];
ringM = Diagonal[Normal[(Total[QuantumOperator["Z", {#}] & /@ Range[5]]/2)["Matrix"]]];
ringBasis[mm_, m_] := With[{p = Simplify[ComplexExpand[Total[Exp[-2 Pi I m #/5] MatrixPower[ringT, #] & /@ Range[0, 4]]/5]]},
   Normalize /@ DeleteDuplicatesBy[DeleteCases[Transpose[p][[Flatten[Position[ringM, mm]]]], {0 ..}], Unitize]];
ringSectors = DeleteCases[AssociationMap[ringBasis @@ # &, Tuples[{Range[5/2, -5/2, -1], Range[-2, 2]}]], {}];
ringV = Join @@ Values[ringSectors];
ringE = FullSimplify[Eigenvalues[Conjugate[#] . ringH . Transpose[#]]] & /@ ringSectors;
ringMs = {1/2, 3/2, 5/2};
Labeled[Row[{ArrayPlot[Abs[ringH], Frame -> True, Axes -> False, FrameTicks -> None, PlotLabel -> "product basis", ImageSize -> 220],
   ArrayPlot[Abs[N[Conjugate[ringV] . ringH . Transpose[ringV]]], Frame -> True, Axes -> False, FrameTicks -> None, PlotLabel -> "sorted by (M, k)", ImageSize -> 220],
   ListPlot[Table[Catenate[KeyValueMap[Thread[{Last[#1], N[#2]}] &, KeySelect[ringE, First[#] == mm &]]], {mm, ringMs}],
    PlotMarkers -> {{"\[EmptyCircle]", 22}, {"\[EmptySquare]", 14}, {"\[FilledCircle]", 7}}, Frame -> True, PlotRange -> {{-2.5, 2.5}, All}, Axes -> False,
    FrameLabel -> {TraditionalForm[HoldForm[5 k/(2 π)]], "E / J"}, ImageSize -> 300, PlotLegends -> Placed[Row[{"M = ", #}] & /@ ringMs, Below]]}, "   "],
  Style[Row[{Length[ringSectors], " nonempty sectors;  ground-state energy ", TraditionalForm[First[MinimalBy[Catenate[Values[ringE]], N]] J]}], 13, FontFamily -> "Helvetica"], Top]
```

The 32 levels fall into blocks no larger than $2\times2$, one for each nonempty pair $(M,k)$, so every energy is a root of a quadratic. The ground state sits at $M=\pm\frac12$ and $k=\pm2\pi/5$, four degenerate states. Every level of a larger $M$ reappears at the same $k$ in each smaller $M$, as the total-spin multiplets of a rotation-invariant $\hat H$ require.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Open quantum systems

Jump operators and rates, including symbolic rates and full rate matrices, go straight into `QuantumEvolve` or `QuantumOperator["Liouvillian"[...]]`. The same objects also compose with Kraus channels in circuits, so a master equation and a noisy circuit can be compared exactly.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Why a qubit's coherence time T2 is at most 2T1


A qubit with level splitting $\omega$ loses its energy, decaying from $|e\rangle$ to $|g\rangle$ at rate $\gamma_1=1/T_1$, and loses its phase through the jump operator $\sqrt{\gamma_\phi}\,\hat\sigma_z$. How long can a superposition of $|e\rangle$ and $|g\rangle$ keep its phase?

With $\hbar=1$, $|e\rangle=|0\rangle$, $|g\rangle=|1\rangle$ and $\hat\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, the Lindblad equation is

$$
\dot\rho=-i\left[\frac{\omega}{2}\hat\sigma_z,\rho\right]+\gamma_1\,\mathcal D[\hat\sigma_-]\rho+\gamma_\phi\,\mathcal D[\hat\sigma_z]\rho,\qquad \mathcal D[\hat L]\rho=\hat L\rho \hat L^\dagger-\frac12\{\hat L^\dagger \hat L,\rho\}.
$$

Starting from $|+\rangle=(|e\rangle+|g\rangle)/\sqrt2$, the exact solution for general $\omega$, $\gamma_1$ and $\gamma_\phi$ is

$$
\langle\hat\sigma_x\rangle(t)=e^{-t/T_2}\cos\omega t,\qquad \frac{1}{T_2}=\frac{\gamma_1}{2}+2\gamma_\phi.
$$

The figure plots it with time in units of $T_1$ and $\omega T_1=6$. The dashed curves are $\pm e^{-t/T_2}$, and the gray curves are the same envelope at $T_2=2T_1$, the value with no dephasing:

```wl
sxA1 = FullSimplify[First[QuantumEvolve[ω/2 QuantumOperator["Z"], {QuantumOperator["-"], QuantumOperator["Z"]} -> {γ1, γφ},
      QuantumState["+"], t]["BlochVector"]], γ1 > 0 && γφ >= 0 && t >= 0 && ω ∈ Reals];
panelA1[γφv_] := Plot[Evaluate[{Exp[-t/2], -Exp[-t/2], sxA1 /. {γ1 -> 1, γφ -> γφv, ω -> 6}, Exp[-(1/2 + 2 γφv) t], -Exp[-(1/2 + 2 γφv) t]}], {t, 0, 6},
    PlotRange -> {-1.05, 1.05}, PlotPoints -> 200,
    PlotStyle -> {Gray, Gray, ColorData[97][1], Directive[ColorData[97][2], Dashed], Directive[ColorData[97][2], Dashed]},
    Frame -> True, FrameLabel -> {Row[{"t/", Subscript["T", 1]}], Row[{"\[LeftAngleBracket]", Subscript["σ", "x"], "\[RightAngleBracket]"}]},
    PlotLabel -> Row[{Subscript["T", 2], "/", Subscript["T", 1], " = ", ToString[1/(1/2 + 2 Rationalize[γφv]), InputForm]}], ImageSize -> Small];
Legended[PlotGrid[{panelA1 /@ {0, 0.1, 0.25}}, ImageSize -> 660],
  Placed[LineLegend[{ColorData[97][1], Directive[ColorData[97][2], Dashed], Gray},
    {Row[{"\[LeftAngleBracket]", Subscript["σ", "x"], "\[RightAngleBracket]"}],
     Row[{"\[PlusMinus]", Superscript["e", Row[{"\[Minus]t/", Subscript["T", 2]}]]}],
     Row[{"\[PlusMinus]", Superscript["e", Row[{"\[Minus]t/", Subscript["T", 2]}]], " at ", Subscript["T", 2], " = 2", Subscript["T", 1], " (no dephasing)"}]},
    LegendLayout -> {"Row", 1}], Below]]
```

Energy loss alone already removes the phase at rate $\gamma_1/2$, and dephasing can only add to that rate, so $T_2\le 2T_1$, with equality only when $\gamma_\phi=0$. In the figure the oscillation touches the dashed envelope at every peak; at $\gamma_\phi=0$ that envelope lies on the gray limit, and as $\gamma_\phi$ grows it pulls inside it.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Master equation, Kraus channels and a circuit: one Ramsey fringe


A Ramsey experiment applies a Hadamard to $|g\rangle$, waits a time $t$, applies a second Hadamard and reads $P(e)$. During the wait the qubit, detuned by $\delta$ from the drive, decays at rate $\gamma_1$ and dephases at rate $\gamma_\phi$, and a residual drive of Rabi frequency $\Omega$ may act. Device simulators model such a wait not with a master equation but as a circuit: a gate followed by textbook noise channels. When do the two give the same fringe?

With $\hbar=1$, $|g\rangle=|0\rangle$ (the state the amplitude-damping channel relaxes to), $|e\rangle=|1\rangle$, $\hat\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, qubit frequency $\omega$ and drive frequency $\omega_d=\omega-\delta$, the wait obeys, in the frame rotating at $\omega_d$,

$$
\dot\rho=-i\left[\frac{\delta}{2}\hat\sigma_z+\frac{\Omega}{2}\hat\sigma_x,\rho\right]+\gamma_1\,\mathcal D[\hat\sigma_-]\rho+\gamma_\phi\,\mathcal D[\hat\sigma_z]\rho .
$$

The circuit cuts the wait into $m$ equal slices. Each slice is the unitary $e^{-i\hat Ht/m}$, then amplitude damping with decay probability $\lambda=1-e^{-\gamma_1t/m}$, then phase damping with $\mu=1-e^{-4\gamma_\phi t/m}$, which shrinks the coherence by $e^{-2\gamma_\phi t/m}$, the same factor as the dephasing term. With $m=1$ the whole wait is one gate and two channels, the usual circuit noise model. Larger $m$ follows the continuous evolution more closely.

Lines are the master equation, dots the circuit with $m=1$ and squares with $m=8$; time is in units of $1/\gamma_1$:

```wl
hamA4[δ_, Ω_] := -δ/2 QuantumOperator["Z"] + Ω/2 QuantumOperator["X"]; jumpsA4 = {QuantumOperator[{{0, 1}, {0, 0}}], QuantumOperator["Z"]};
sliceA4[h_, {γ1v_, γφv_}, τ_, {lu_, lad_, lpd_}] := {QuantumOperator[Exp[-I τ h], "Label" -> lu],
    QuantumChannel["AmplitudeDamping"[1 - Exp[-γ1v τ]], "Label" -> lad], QuantumChannel["PhaseDamping"[1 - Exp[-4 γφv τ]], "Label" -> lpd]};
pEA4[ρ_] := ρ["ProbabilitiesList"][[2]];
viaMEA4[h_, γs_, T_] := pEA4[QuantumOperator["H"][QuantumEvolve[h, jumpsA4 -> γs, QuantumOperator["H"][QuantumState["0"]], T]]];
circuitA4[h_, γs_, m_, s_, labels_ : {"U", "AD", "PD"}, group_ : Catenate] := QuantumCircuitOperator[Join[{"H"}, group[ConstantArray[sliceA4[h, γs, s/m, labels], m]], {"H"}]];
diagramsA4 = {circuitA4[hamA4[δ, Ω], {γ1, γφ}, 1, t, {Superscript["e", "\[Minus]i H t"], Column[{"Amplitude", "damping"}, Alignment -> Center], Column[{"Phase", "damping"}, Alignment -> Center]}, Map[QuantumCircuitOperator[#, "slice"] &]]["Diagram", ImageSize -> 330],
    circuitA4[hamA4[δ, Ω], {γ1, γφ}, 8, t, {"U", "AD", "PD"}, Map[QuantumCircuitOperator[#, "slice"] &]]["Diagram", ImageSize -> 760]};
transferA4[h_, γs_] := With[{ρ = Array[\[FormalR], {2, 2}]}, Function[\[FormalS], Evaluate[D[Flatten[Normal[
      QuantumCircuitOperator[sliceA4[h, γs, \[FormalT], {"U", "AD", "PD"}]][QuantumState[ρ]]["DensityMatrix"]]], {Flatten[ρ]}] /. \[FormalT] -> \[FormalS]]]];
pointsA4[T_, m_] := With[{hm = Normal[QuantumOperator["H"]["Matrix"]]},
   Table[{x, Re[(hm . Partition[MatrixPower[T[x/m], m] . Flatten[hm . {{1, 0}, {0, 0}} . hm], 2] . hm)[[2, 2]]]}, {x, 0, 4, If[m == 1, 0.2, 0.4]}]];
panelA4[{γφv_, Ωv_}, δv_] := With[{me = viaMEA4[hamA4[δv, Ωv], {1, γφv}, {t, 0, 4}], T = transferA4[N[hamA4[δv, Ωv]], {1, γφv}]},
   ListPlot[Prepend[pointsA4[T, #] & /@ If[Ωv == 0, {1}, {1, 8}], Table[{x, me /. t -> x}, {x, 0, 4, 0.04}]],
    Joined -> {True, False, False}, PlotMarkers -> {None, {"\[FilledCircle]", 7}, {"\[EmptySquare]", 11}},
    PlotLegends -> Take[{"master equation", "circuit, m = 1", "circuit, m = 8"}, If[Ωv == 0, 2, 3]], PlotRange -> {0, 1}, Frame -> True, FrameLabel -> {Row[{Subscript["γ", 1], " t"}], "P(e)"},
    Epilog -> Inset[Row[{"δ/", Subscript["γ", 1], " = ", δv}], Scaled[{0.96, 0.96}], {Right, Top}], ImageSize -> Small]];
casesA4 = {{0, 0}, {1/2, 0}, {1/2, 2}};
panelsA4 = Outer[panelA4, casesA4, {0, 3, 8}, 1];
Column[{Column[diagramsA4, Alignment -> Center, Spacings -> 1],
  Manipulate[Legended[PlotGrid[{First /@ panelsA4[[c]]}, ImageSize -> 660], legendBelow[Last[panelsA4[[c, 1]]]]],
   Row[{Control[{{c, 3, "case"}, MapIndexed[First[#2] -> Row[{Subscript["γ", "φ"], "/", Subscript["γ", 1], " = ", #1[[1]], ",  Ω/", Subscript["γ", 1], " = ", #1[[2]]}] &, casesA4], SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"case ", c}], <|1 -> "No dephasing, no drive: the routes agree; the fringe decays at γ₁/2.",
        2 -> "Dephasing added: the routes still agree; decay at 1/T₂ = γ₁/2 + 2γ_φ.",
        3 -> "Drive on: the circuit matches only when the wait is sliced finely."|>[c]}}]]}], SaveDefinitions -> True]}, Alignment -> Center, Spacings -> 2]
```

Without a drive the two routes agree for every value of every parameter, because the rotation, the decay and the dephasing generators commute: the wait is exactly one rotation followed by two named channels. The fringe oscillates about $P(e)=1/2$ at the detuning, and $\gamma_1$ enters its decay only through $1/T_2$. A drive during the wait stops the generators commuting, and the circuit approaches the master equation only as the wait is sliced finer.

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

The concurrence of the evolving state, computed for general $\theta$, rates and $t$, is

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
    GridLines -> {None, {0}}, Epilog -> {PointSize[0.035], MapIndexed[If[NumericQ[#1], {ColorData[97][First[#2]], Point[{#1, 0}]}, Nothing] &, deathA2[γBv, #] & /@ θsA2]},
    PlotStyle -> {Automatic, Dashed, Automatic, Automatic}, PlotLegends -> (Row[{"θ = ", #}] & /@ θsA2),
    Frame -> True, FrameLabel -> {Row[{Subscript["γ", "A"], " t"}], "concurrence"}, PlotLabel -> Row[{Subscript["γ", "B"], "/", Subscript["γ", "A"], " = ", γBv}], ImageSize -> Small];
panelsA2 = panelA2 /@ γBsA2;
Legended[PlotGrid[{First /@ panelsA2}, ImageSize -> 660], legendBelow[Column[{Last[First[panelsA2]], PointLegend[{Gray}, {"entanglement dies"}, LegendMarkers -> {"\[FilledCircle]", 9}]}]]]
```

It vanishes once $(1-e^{-\gamma_At})(1-e^{-\gamma_Bt})$ reaches $\cot^2\theta$. That needs both qubits to decay and $\theta>\pi/4$, where $|ee\rangle$ outweighs $|gg\rangle$: those states lose all their entanglement at a finite time and never regain it, and the more slowly B decays, the later that time comes. At $\theta\le\pi/4$ (the boundary case $\theta=\pi/4$ is dashed), or if B does not decay, the entanglement only fades.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Two atoms in a shared field: superradiance, subradiance and a dark state


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
    PlotRange -> {{0, 5}, {0, 1.02}}, Frame -> True, Axes -> False,
    FrameLabel -> {"γ t", Row[{"\[LeftAngleBracket]N\[RightAngleBracket] / ", Subscript["N", 0]}]},
    PlotLabel -> Row[{"r = ", Style[TraditionalForm[rv], FractionBoxOptions -> {Beveled -> True}]}], ImageSize -> 260];
Legended[PlotGrid[{panelA3 /@ {0, 1/2, 1}}, ImageSize -> 760],
  Placed[Row[{LineLegend[ColorData[97] /@ {3, 4, 5, 1}, Keys[statesA3], LegendLayout -> "Row"], Spacer[30],
     SwatchLegend[{Directive[Opacity[0.22], superColorA3], Directive[Opacity[0.22], subColorA3]}, {"faster than a lone atom: superradiant", "slower than a lone atom: subradiant"}, LegendLayout -> "Row"]}], Below]]
```

Without correlation, at $r=0$, the atoms decay independently and all four curves, $|ee\rangle$ per excitation included, lie on that of a lone atom, so they show as one curve. Correlation makes the symmetric state superradiant, decaying faster than a lone atom, and the antisymmetric one subradiant, decaying slower: $\langle\hat N\rangle_T=e^{-(1+r)\gamma t}$ and $\langle\hat N\rangle_S=e^{-(1-r)\gamma t}$. At $r=1$ the singlet is dark, so $|eg\rangle$, half singlet, keeps half an excitation, while $|ee\rangle$ empties completely and, per excitation, faster than a lone atom.

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

### Mollow triplet and squeezed resonance fluorescence


A two-level atom, frequency $\omega_a$, decays at rate $\gamma$ and is driven by a laser at frequency $\omega_L$ with Rabi frequency $\Omega$. Does the scattered light carry only the Mollow triplet a spectrometer sees, or is one of its quadratures squeezed below shot noise?

With $\hbar=1$, $|e\rangle=|0\rangle$, $|g\rangle=|1\rangle$, $\hat\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$ and $\Delta=\omega_a-\omega_L$, the Hamiltonian in the frame rotating at $\omega_L$ (rotating-wave approximation) is $\hat H=\frac{\Delta}{2}\hat\sigma_z+\frac{\Omega}{2}\hat\sigma_x$, and

$$
\dot\rho=\mathcal L\rho=-i[\hat H,\rho]+\gamma\,\mathcal D[\hat\sigma_-]\rho .
$$

By the quantum regression theorem, two-time correlations evolve with the same $\mathcal L$ from the steady state $\rho_{ss}$, so at frequency $\nu$ from the laser the incoherent fluorescence spectrum and the homodyne spectrum of the quadrature $\hat x=\hat c+\hat c^\dagger$ of the emitted field $\hat c=-i\sqrt\gamma\,\hat\sigma_-$ (out of phase with the drive), normalized to shot noise, are

$$
S(\nu)=2\,\mathrm{Re}\,\mathrm{Tr}\!\left[\hat\sigma_+(i\nu-\mathcal L)^{-1}\delta\hat\sigma_-\rho_{ss}\right],\qquad
S_x(\nu)=1+2\,\mathrm{Re}\,\mathrm{Tr}\!\left[\hat x\,(i\nu-\mathcal L)^{-1}\left(\hat c\rho_{ss}+\rho_{ss}\hat c^\dagger-\langle\hat x\rangle\rho_{ss}\right)\right],
$$

with $\delta\hat\sigma_-=\hat\sigma_--\langle\hat\sigma_-\rangle$.

The null vector of $\mathcal L$, for general $\Omega$ and $\Delta$, gives $\rho_{ss}$, and inverting $i\nu-\mathcal L$ in closed form gives both spectra; each fluorescence spectrum is scaled to its own maximum, and rates are in units of $\gamma$:

```wl
liouA6 = Normal[QuantumOperator["Liouvillian"[Δ/2 QuantumOperator["Z"] + Ω/2 QuantumOperator["X"], {QuantumOperator["-"]}, {1}]]["Matrix"]];
smA6 = Normal[QuantumOperator["-"]["Matrix"]]; cA6 = -I smA6; xA6 = cA6 + ConjugateTranspose[cA6];
ssA6 = #/Tr[#] &[ArrayReshape[First[NullSpace[liouA6]], {2, 2}]];
spectrumA6[seed_, readout_] := 2 Re @ Tr[readout . ArrayReshape[LinearSolve[I ν IdentityMatrix[4] - liouA6, Flatten[seed]], {2, 2}]];
fluorA6 = spectrumA6[(smA6 - Tr[smA6 . ssA6] IdentityMatrix[2]) . ssA6, ConjugateTranspose[smA6]];
quadA6 = 1 + spectrumA6[cA6 . ssA6 + ssA6 . ConjugateTranspose[cA6] - Tr[xA6 . ssA6] ssA6, xA6];
{ΩsA6, ΔsA6} = {{0.3, 1, 5}, {0, 1, 3}};
curveA6[1, Ωv_, Δv_] := With[{f = fluorA6 /. {Ω -> Ωv, Δ -> Δv}}, f/Max[Table[f, {ν, -8., 8., 0.02}]]];
curveA6[2, Ωv_, Δv_] := quadA6 /. {Ω -> Ωv, Δ -> Δv};
panelA6[row_, Ωv_] := Plot[Evaluate[curveA6[row, Ωv, #] & /@ ΔsA6], {ν, -8, 8}, PlotRange -> All, PlotLegends -> (Row[{"Δ/γ = ", #}] & /@ ΔsA6),
    Frame -> True, Axes -> False, GridLines -> {None, If[row == 2, {1}, None]}, FrameLabel -> {"ν/γ", {"S(ν)/max", Row[{Subscript["S", "x"], "(ν)"}]}[[row]]},
    ImageSize -> 300];
panelsA6 = AssociationMap[Function[Ωv, panelA6[#, Ωv] & /@ {1, 2}], ΩsA6];
Manipulate[Legended[Grid[{First /@ panelsA6[w]}], legendBelow[Last[First[panelsA6[w]]]]], Row[{Control[{{w, 0.3, "Ω/γ"}, ΩsA6, SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"Ω/γ = ", w}], <|0.3 -> "One line at the laser; on resonance, squeezed below shot noise at ν = 0.",
        1 -> "One broadened line; the squeezing at ν = 0 is gone.",
        5 -> "Mollow triplet, sidebands at ±Ω; no squeezing."|>[w]},
      {"Δ ≠ 0", "Light moves into sidebands near ±Δ, pushed further out at strong drive."}}]]}], SaveDefinitions -> True]
```

At weak resonant driving the incoherent light is a single line at the laser frequency; at $\Omega/\gamma=5$ it splits into the Mollow triplet, whose sidebands detuning pushes further out, and with detuning even weak driving moves the light into sidebands near $\pm\Delta$. At weak resonant driving the out-of-phase quadrature dips below shot noise around $\nu=0$: a single atom emits squeezed light. Stronger driving or detuning lifts it back above shot noise there.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A driven Kerr cavity: sub- or super-Poissonian steady light


A cavity mode of frequency $\omega_c$ with a Kerr nonlinearity $K$ is driven at amplitude $F$ and frequency $\omega_d$, and leaks photons at rate $\kappa$. Drive and loss balance in a steady state. For which detuning, Kerr strength and drive is the steady light sub-Poissonian?

With $\hbar=1$, the Hamiltonian, and its form in the frame rotating at $\omega_d$, are

$$
\hat H=\omega_c\,\hat a^\dagger \hat a+\frac K2\,\hat a^{\dagger2}\hat a^2+F\left(\hat a\,e^{i\omega_d t}+\hat a^\dagger e^{-i\omega_d t}\right)\;\longrightarrow\;-\Delta\,\hat a^\dagger \hat a+\frac K2\,\hat a^{\dagger2}\hat a^2+F\left(\hat a+\hat a^\dagger\right),\qquad \Delta=\omega_d-\omega_c,
$$

and the loss enters through $\dot\rho=-i[\hat H,\rho]+\kappa\,\mathcal D[\hat a]\rho=\mathcal L\rho$. Light is sub-Poissonian, which no classical field is, when the Mandel parameter $Q=(\langle\hat n^2\rangle-\langle\hat n\rangle^2-\langle\hat n\rangle)/\langle\hat n\rangle$ is negative. The steady state is known in closed form (Drummond and Walls): with $c=(\kappa-2i\Delta)/(iK)$ and $x=2F/K$, $\langle\hat a^{\dagger j}\hat a^j\rangle=|x|^{2j}\frac{\Gamma(c)\Gamma(c^*)}{\Gamma(c+j)\Gamma(c^*+j)}\,\frac{{}_0F_2(;c+j,c^*+j;2|x|^2)}{{}_0F_2(;c,c^*;2|x|^2)}$. The smallest nonzero decay rate of $\mathcal L$, the gap $\lambda$, sets how long the slowest transient lasts.

The vacuum evolves under the exponential of $\mathcal L$, and the gap (legends) comes from its spectrum; the closed-form steady $Q$ is dashed; rates are in units of $\kappa$, with 24 Fock levels:

```wl
{kerrLevels, kerrTimes} = {24, Range[0.1, 20, 0.1]};
kerrA = AnnihilationOperator[kerrLevels];
kerrH[Δ_, K_, F_] := -Δ SuperDagger[kerrA] @ kerrA + K/2 SuperDagger[kerrA] @ SuperDagger[kerrA] @ kerrA @ kerrA + F (kerrA + SuperDagger[kerrA]);
kerrQ[ρ_] := With[{p = Re[Diagonal[ρ]/Tr[ρ]], m = Range[0, kerrLevels - 1]}, (m^2 . p - (m . p)^2)/(m . p) - 1];
kerrMoment[Δ_, K_, F_, j_] := With[{c = N[(1 - 2 I Δ)/(I K)], x2 = N[(2 F/K)^2]},
    x2^j Gamma[c] Gamma[Conjugate[c]]/(Gamma[c + j] Gamma[Conjugate[c] + j]) HypergeometricPFQ[{}, {c + j, Conjugate[c] + j}, 2 x2]/HypergeometricPFQ[{}, {c, Conjugate[c]}, 2 x2]];
kerrQss[Δ_, K_, F_] := With[{n1 = kerrMoment[Δ, K, F, 1]}, Re[(kerrMoment[Δ, K, F, 2] - n1^2)/n1]];
kerrSolve[Δ_, K_, F_] := With[{L = N @ Normal @ QuantumOperator["Liouvillian"[kerrH[Δ, K, F], {kerrA}, {1}]]["Matrix"]}, {step = MatrixExp[First[kerrTimes] L]},
    {Transpose[{kerrTimes, kerrQ[Partition[#, kerrLevels]] & /@ Rest @ NestList[step . # &, UnitVector[kerrLevels^2, 1], Length[kerrTimes]]}],
     kerrQss[Δ, K, F], Min @ Select[-Re @ Eigenvalues[L], # > 10^-8 &]}];
kerrPanel[Δ_, K_] := With[{Fs = {0.5, 1.5, 3}}, {ss = kerrSolve[Δ, K, #] & /@ Fs},
    Show[
      ListLinePlot[ss[[All, 1]], PlotRange -> All,
        PlotLegends -> MapThread[Row[{"F/\[Kappa] = ", #1, ",   \[Lambda]/\[Kappa] = ", NumberForm[#2, 2]}] &, {Fs, ss[[All, 3]]}]],
      Plot[Evaluate[ss[[All, 2]]], {t, 0, Last[kerrTimes]}, PlotStyle -> Directive[Thin, Dashed]],
      PlotRange -> All, Frame -> True, FrameLabel -> {"\[Kappa] t", "Mandel Q"},
      PlotLabel -> Row[{"\[CapitalDelta]/\[Kappa] = ", Δ, ",   K/\[Kappa] = ", K}], ImageSize -> 250]];
kerrPanels = AssociationMap[kerrPanel @@ # &, Tuples[{{0, 2.5}, {0.5, 2.5, 10}}]];
Manipulate[kerrPanels[{d, k}], Row[{Column[{Control[{{d, 0, "\[CapitalDelta]/\[Kappa]"}, {0, 2.5}, SetterBar}], Control[{{k, 0.5, "K/\[Kappa]"}, {0.5, 2.5, 10}, SetterBar}]}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"Δ/κ = ", d, ",  K/κ = ", k}], <|{0, 0.5} -> "Sub-Poissonian at every drive; settled within about ten lifetimes.",
        {0, 2.5} -> "Sub-Poissonian at every drive; settled within a few lifetimes.",
        {0, 10} -> "Sub-Poissonian at every drive, deepest at F/κ = 1.5.",
        {2.5, 0.5} -> "F/κ = 1.5 turns super-Poissonian; its gap nearly closes, unsettled at κt = 20.",
        {2.5, 2.5} -> "The weakest drive turns super-Poissonian; stronger drives stay sub-Poissonian.",
        {2.5, 10} -> "The weakest drive is near Poissonian; stronger drives stay sub-Poissonian."|>[{d, k}]}}]]}], SaveDefinitions -> True]
```

Every transient but one settles on the closed-form steady value within twenty photon lifetimes. On resonance every steady state is sub-Poissonian, and the cavity reaches it within a few photon lifetimes, with gaps near the empty-cavity value $\kappa/2$. Detuning makes the sign of $Q$ depend on the drive and the Kerr strength: some steady states become super-Poissonian, and in one of them the gap nearly closes, so $Q$ is still far from its steady value after twenty photon lifetimes.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A cat in a warm bath: Caldeira-Leggett decoherence


A particle of mass $M$ in a trap of frequency $\Omega$ couples to a bath at temperature $T$ that damps its momentum at rate $2\gamma$. It starts in a cat state, the superposition $|\alpha\rangle+|-\alpha\rangle$ of two coherent states on either side of the trap. How fast does the bath turn it into a mixture of the two lobes?

With $\hbar=1$, the high-temperature Caldeira-Leggett master equation is

$$
\dot\rho=-i\left[\frac{\hat p^2}{2M}+\frac12M\Omega^2\hat x^2,\rho\right]-i\gamma\left[\hat x,\{\hat p,\rho\}\right]-2\gamma Mk_BT\,\left[\hat x,[\hat x,\rho]\right]-\frac{\gamma}{8Mk_BT}\left[\hat p,[\hat p,\rho]\right],
$$

where the last term, smaller than the others by a factor of order $(\Omega/k_BT)^2$, makes the generator a Lindblad generator. With jump operators $\hat x$ and $\hat p$ it takes the Hamiltonian $\hat H+\frac\gamma2\{\hat x,\hat p\}$ and the rate matrix $\begin{pmatrix}4\gamma Mk_BT & -i\gamma\\ i\gamma & \gamma/4Mk_BT\end{pmatrix}$.

The equation is linear, so the difference between the cat and the mixture of its two lobes, $\Delta\rho=\rho_{\rm cat}-\rho_{\rm mix}\propto C-s\left(|\alpha\rangle\langle\alpha|+|-\alpha\rangle\langle-\alpha|\right)$, with the cross term $C=|\alpha\rangle\langle-\alpha|+|-\alpha\rangle\langle\alpha|$ and the overlap $s=\langle-\alpha|\alpha\rangle=e^{-2\alpha^2}$, evolves on its own. The coherence left at time $t$ is its trace norm relative to the start,

$$
D(t)=\frac{\|\mathcal E_t(\Delta\rho)\|_1}{\|\Delta\rho\|_1},
$$

which reaches zero once the cat has become the mixture; the cross term alone would not, because its trace $2s$ is conserved.

The $[\hat x,[\hat x,\rho]]$ term damps $\langle x|C|x'\rangle$ at rate $2\gamma Mk_BT(x-x')^2$; the lobes sit $2\sqrt2\,\alpha\cos\Omega t/\sqrt{M\Omega}$ apart for $\hat x=(\hat a+\hat a^\dagger)/\sqrt{2M\Omega}$, which, for lobes of no width, predicts $D(t)=e^{-8\gamma k_BT\alpha^2(\Omega t+\sin2\Omega t/2)/\Omega^2}$.

$\Delta\rho$ evolves under the exponential of $\mathcal L$ on 40 Fock levels of the trap; the prediction is dashed, and units are $M=\Omega=1$:

```wl
qbmLevels = 40; qbmA = AnnihilationOperator[qbmLevels]; qbmX = (qbmA + SuperDagger[qbmA])/Sqrt[2]; qbmP = I (SuperDagger[qbmA] - qbmA)/Sqrt[2];
qbmStep[γ_, kT_] := qbmStep[γ, kT] = MatrixExp[π/32 N @ Normal @ QuantumOperator["Liouvillian"[(qbmP @ qbmP + qbmX @ qbmX)/2 + γ/2 (qbmX @ qbmP + qbmP @ qbmX),
      {qbmX, qbmP}, {{4 γ kT, -I γ}, {I γ, γ/(4 kT)}}]]["Matrix"]];
qbmDiff[α_] := With[{u = Normal[CoherentState[qbmLevels][N[α]]["StateVector"]], v = Normal[CoherentState[qbmLevels][-N[α]]["StateVector"]]},
    Flatten[Outer[Times, u, Conjugate[v]] + Outer[Times, v, Conjugate[u]] - Conjugate[v] . u (Outer[Times, u, Conjugate[u]] + Outer[Times, v, Conjugate[v]])]];
qbmCoherence[γ_, kT_, α_] := With[{c = Total[SingularValueList[Partition[#, qbmLevels]]] & /@ NestList[qbmStep[γ, kT] . # &, qbmDiff[α], 64]},
    Transpose[{Range[0, 2 π, π/32], c/First[c]}]];
qbmγs = {0.00005, 0.0001, 0.0002};
qbmPanel[α_, kT_] := Show[
    ListLogPlot[qbmCoherence[#, kT, α] & /@ qbmγs, Joined -> True, PlotLegends -> (Row[{"\[Gamma]/\[CapitalOmega] = ", #}] & /@ qbmγs)],
    LogPlot[Evaluate[Exp[-8 # kT α^2 (t + Sin[2 t]/2)] & /@ qbmγs], {t, 0, 2 π}, PlotStyle -> Directive[Gray, Dashed]],
    PlotRange -> {Log[0.18], Log[1.05]}, Axes -> False, Frame -> True, FrameLabel -> {"\[CapitalOmega] t", "D"},
    FrameTicks -> {{{Log[#], #} & /@ {0.2, 0.3, 0.5, 0.7, 1}, None}, {Table[{k π/2, k π/2}, {k, 0, 4}], None}}, GridLines -> {Table[k π/2, {k, 4}], None}, AspectRatio -> 1/2,
    PlotLabel -> Row[{"\[Alpha] = ", α, ",   ", Subscript["k", "B"], "T/\[CapitalOmega] = ", kT}], ImageSize -> 340];
qbmPanels = Outer[qbmPanel, {3, 4}, {5, 10}];
Legended[PlotGrid[Map[First, qbmPanels, {2}], ImageSize -> 720], legendBelow[Column[{Last[qbmPanels[[1, 1]]], LineLegend[{Directive[Gray, Dashed]}, {"closed form"}]}]]]
```

The coherence falls in steps: it drops while the two lobes sit apart in position, near $\Omega t=0,\pi,2\pi$, and stalls near $\pi/2$ and $3\pi/2$, where the trap has turned their separation into momentum, which the $[\hat x,[\hat x,\rho]]$ term does not see. The rate grows as $\gamma k_BT\alpha^2$, so a larger cat, a hotter bath or stronger friction makes every step deeper. The friction here is weak enough that the bath adds little noise to each lobe within a period, so the lobes stay sharp and the solid curves follow the prediction, which treats each lobe as a point. A gap opens only for the large cat in the hotter bath at the strongest friction, where the bath begins to heat the lobes themselves and the coherence decays more slowly than the prediction.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A qutrit Wigner function under dephasing: does the negativity die at an exact time?


A qutrit state can be written as a discrete Wigner function, $W(q,p)=\frac13\mathrm{Tr}\left[\rho\,\hat A_{q,p}\right]$ with $\hat A_{q,p}=\hat D_{q,p}\hat P\hat D_{q,p}^\dagger$, $\hat D_{q,p}=X^qZ^p$ and $\hat P|j\rangle=|-j\rangle$: nine real numbers on a $3\times3$ phase-space grid that sum to one, whose negative entries make it a quasi-probability rather than a probability. A qutrit with degenerate levels, so no Hamiltonian, is dephased at rate $\gamma$. Does every negative entry reach zero at a finite time?

With $X|j\rangle=|j+1\rangle$ and $Z|j\rangle=e^{-2\pi ij/3}|j\rangle$ (indices mod 3), the dephasing has jump operators $Z$ and $Z^2$,

$$
\dot\rho=\frac{\gamma}{3}\left(\mathcal D[Z]\rho+\mathcal D[Z^2]\rho\right),
$$

which keeps the populations $\rho_{jj}$ and damps every coherence $\rho_{jk}$ as $e^{-\gamma t}$, so each entry of $W$ relaxes along one exponential toward $W_\infty$, the Wigner function of the populations alone:

$$
W_t(q,p)=W_\infty(q,p)+e^{-\gamma t}\left(W_0(q,p)-W_\infty(q,p)\right).
$$

The solution for general $\gamma$ and $t$ from three initial states gives $W(t)$ on the $3\times3$ grid, shown at $\gamma t=0$, $\gamma t^*/2$, $\gamma t^*$ and $\infty$, where $t^*$ is the time at which the last negative entry reaches zero (at $\gamma t=0,1,2,\infty$ when no such time exists); every grid shares one color scale, and time is in units of $1/\gamma$:

```wl
dephA7 = {QuantumOperator["Z"[3]], QuantumOperator["Z"[3]] @ QuantumOperator["Z"[3]]};
wignerA7[ψ_] := Simplify[ComplexExpand[Normal[QuantumEvolve[QuantumOperator[0 IdentityMatrix[3]], dephA7 -> {γ/3, γ/3},
       QuantumState[ψ], t]["PhaseSpace"]]] /. γ -> 1, t >= 0];
deathA7[w_] := Max[Replace[Solve[# == 0 && t > 0, t, Reals], {{} -> ∞, sol_ :> (t /. First[sol])}] & /@ Select[Flatten[w], (# /. t -> 0) < 0 &]];
statesA7 = KeyMap[TraditionalForm, <|HoldForm[Divide[Ket[{1}] - Ket[{2}], Sqrt[2]]] -> {0, 1, -1}/Sqrt[2],
    HoldForm[Divide[Ket[{0}] + I Ket[{1}] - Ket[{2}], Sqrt[3]]] -> {1, I, -1}/Sqrt[3],
    HoldForm[Divide[Sqrt[2] Ket[{0}] + Ket[{1}] - Ket[{2}], 2]] -> {Sqrt[2], 1, -1}/2|>];
colorA7 = ColorData["ThermometerColors"][Rescale[#, {-1/3, 1/3}]] &;
gridA7[w_, tv_] := With[{m = Limit[w, t -> tv]},
   MatrixPlot[Transpose[m], ColorFunction -> colorA7, ColorFunctionScaling -> False, DataReversed -> True, Mesh -> All,
    FrameTicks -> {{{{1, 0}, {2, 1}, {3, 2}}, None}, {{{1, 0}, {2, 1}, {3, 2}}, None}}, FrameLabel -> {"p", "q"},
    Epilog -> MapIndexed[Text[Style[If[#1 == 0, "0", NumberForm[N[#1], 2]], Black], #2 - 1/2] &, m, {2}],
    PlotLabel -> Row[{"γ t = ", TraditionalForm[tv]}], ImageSize -> Small]];
rowsA7 = Function[ψ, With[{w = wignerA7[ψ]}, {ts = deathA7[w]}, gridA7[w, #] & /@ If[ts === ∞, {0, 1, 2, ∞}, {0, ts/2, ts, ∞}]]] /@ statesA7;
Manipulate[Legended[Grid[{Values[rowsA7][[s]]}], Placed[BarLegend[{colorA7, {-1/3, 1/3}}, LegendLabel -> "W"], Right]], Row[{Control[{{s, Length[rowsA7], "state"}, MapIndexed[First[#2] -> #1 &, Keys[rowsA7]], SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Keys[rowsA7][[s]], <|1 -> "Its negative entry sits where W_∞ = 0: it only decays, never dies.",
        2 -> "Negative entries cross zero at different times, the last at γt = ln 2.",
        3 -> "Negative entries cross zero at different times, the last at γt = 3 ln 2/2."|>[s]}}]]}], SaveDefinitions -> True]
```

An entry that starts negative and whose limit $W_\infty$ is positive crosses zero at the exact time $\gamma t^*=\ln\left(1+|W_0|/W_\infty\right)$, and different phase-space points cross at different times. With no population in $|0\rangle$, the single negative entry of $(|1\rangle-|2\rangle)/\sqrt2$ sits where $W_\infty$ is zero, so it only decays and the negativity never dies.

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

### One photon between an atom and a leaky cavity


An excited two-level atom, frequency $\omega_a$, sits in an empty cavity, frequency $\omega_c$. It hands its excitation to the field at coupling $g$, decays on its own at rate $\gamma$, and the mirror leaks photons at rate $\kappa$. Is the photon ever reabsorbed?

With $\hbar=1$ and $\hat\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, the Jaynes-Cummings Hamiltonian, and its form in the frame rotating at $\omega_c$, are

$$
\hat H=\omega_c\,\hat a^\dagger \hat a+\frac{\omega_a}{2}\,\hat\sigma_z+g\left(\hat\sigma_-\hat a^\dagger+\hat\sigma_+\hat a\right)\;\longrightarrow\;\frac{\Delta}{2}\,\hat\sigma_z+g\left(\hat\sigma_-\hat a^\dagger+\hat\sigma_+\hat a\right),\qquad \Delta=\omega_a-\omega_c,
$$

and the losses enter through the Lindblad equation $\dot\rho=-i[\hat H,\rho]+\gamma\,\mathcal D[\hat\sigma_-]\rho+\kappa\,\mathcal D[\hat a]\rho$.

Starting from $|e,0\rangle$, the figure follows the excited-state population $P_e$; with a single excitation, two field levels suffice. Rates are in units of $g$:

```wl
a = AnnihilationOperator[2, {2}]; \[Sigma]m = QuantumOperator["-"];
h[\[CapitalDelta]_] := \[CapitalDelta]/2 QuantumOperator["Z"] + \[Sigma]m @ SuperDagger[a] + SuperDagger[\[Sigma]m] @ a;
{\[Kappa]s, \[CapitalDelta]s, \[Gamma]s} = {{0, 0.1, 1, 5}, {0, 2, 6}, {0, 0.1, 0.5}};
qf = AssociationMap[QuantumPartialTrace[QuantumEvolve[h[#[[2]]], {\[Sigma]m, a} -> {#[[3]], #[[1]]},
       QuantumTensorProduct[QuantumState["0"], FockState[0, 2]], {t, 0, 10}], {2}]["ProbabilitiesList"][[1]] &, Tuples[{\[Kappa]s, \[CapitalDelta]s, \[Gamma]s}]];
panel[\[Kappa]_, \[CapitalDelta]_] := Show[
    Plot[Evaluate[qf[{\[Kappa], \[CapitalDelta], #}] & /@ \[Gamma]s], {t, 0, 10}, PlotRange -> All, PlotLegends -> (Row[{"\[Gamma]/g = ", #}] & /@ \[Gamma]s)],
    If[\[Kappa] == Last[\[Kappa]s], Plot[Evaluate[Exp[-(# + 4 \[Kappa]/(\[Kappa]^2 + 4 \[CapitalDelta]^2)) t] & /@ \[Gamma]s], {t, 0, 10}, PlotStyle -> Dashed, PlotRange -> All], Graphics[]],
    PlotRange -> {0, 1}, Frame -> True, FrameLabel -> {"g t", "P(e)"}, ImageSize -> Small];
cavityPanels = AssociationMap[Function[\[Kappa]v, Show[First[panel[\[Kappa]v, #]], Graphics[Inset[Style[Row[{"\[CapitalDelta]/g = ", #}], 12], Scaled[{0.95, 0.95}], {Right, Top}]]] & /@ \[CapitalDelta]s], \[Kappa]s];
Manipulate[Legended[PlotGrid[{cavityPanels[\[Kappa]]}, ImageSize -> 660],
    Placed[LineLegend[Append[ColorData[97] /@ Range[Length[\[Gamma]s]], If[\[Kappa] == Last[\[Kappa]s], Directive[Gray, Dashed], Nothing]],
      Append[Row[{"\[Gamma]/g = ", #}] & /@ \[Gamma]s, If[\[Kappa] == Last[\[Kappa]s], "rate-equation decay", Nothing]], LegendLayout -> {"Row", 1}], Below]],
  Row[{Control[{{\[Kappa], 1, "\[Kappa]/g"}, \[Kappa]s, SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"κ/g = ", \[Kappa]}], <|0 -> "No leak: the photon swaps back and forth; only γ damps it.",
        0.1 -> "Weak leak: the photon returns many times before escaping.",
        1 -> "Damped exchange: one or two weak returns, then decay.",
        5 -> "Bad cavity: plain decay near the Purcell rate (dashed)."|>[\[Kappa]]},
      {"Δ ≠ 0", "Faster but incomplete exchange, slower leak."}}]]}], SaveDefinitions -> True]
```

At small $\kappa$ the photon returns and $P_e$ oscillates; by $\kappa/g=5$ the atom simply decays, near the Purcell rate $\gamma+4g^2\kappa/(\kappa^2+4\Delta^2)$ (dashed), which becomes exact only for $\kappa\gg g$. Detuning makes the exchange faster but incomplete and slows the leak through the cavity.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Event-triggered control of a decaying qubit


A qubit decays from $|e\rangle$ at rate $\gamma$ under a drive of Rabi frequency $\Omega$, detuned by $\Delta$. A controller that tracks the excited population $P_e$ from the model, without measuring the qubit, switches the drive off when $P_e$ rises above $p+0.02$ and back on when it falls below $p-0.02$. A constant drive lets $P_e$ settle at or below $1/2$; can the switching hold it at a target $p$ above that?

With $\hbar=1$, $|e\rangle=|0\rangle$, $|g\rangle=|1\rangle$ and $\hat\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, the driven qubit in the frame rotating at the drive frequency is

$$
\hat H=\frac{\omega_q}{2}\hat\sigma_z+u(t)\cos(\omega_d t)\,\hat\sigma_x\;\longrightarrow\;\frac{\Delta}{2}\hat\sigma_z+\frac{u(t)}{2}\hat\sigma_x,\qquad \Delta=\omega_q-\omega_d,\quad u(t)\in\{0,\Omega\},
$$

and $\dot\rho=-i[\hat H,\rho]+\gamma\,\mathcal D[\hat\sigma_-]\rho$ is, for the Bloch vector $(x,y,z)$ with $P_e=(1+z)/2$,

$$
\dot x=-\Delta y-\frac{\gamma}{2}x,\qquad \dot y=\Delta x-\frac{\gamma}{2}y-u\,z,\qquad \dot z=u\,y-\gamma(1+z).
$$

These are integrated from $|g\rangle$ with $u$ switched between its two values whenever the switching condition is met, against the constant drive $u=\Omega$ (dashed); time is in units of $1/\gamma$:

```wl
controlRun[Ω_, Δ_, p_] := NDSolveValue[{controlX'[t] == -Δ controlY[t] - controlX[t]/2,
    controlY'[t] == Δ controlX[t] - controlY[t]/2 - controlU[t] controlZ[t], controlZ'[t] == controlU[t] controlY[t] - 1 - controlZ[t],
    controlX[0] == 0, controlY[0] == 0, controlZ[0] == -1, controlU[0] == Ω,
    WhenEvent[(1 + controlZ[t])/2 > p + 0.02, controlU[t] -> 0], WhenEvent[(1 + controlZ[t])/2 < p - 0.02, controlU[t] -> Ω]},
    controlZ, {t, 0, 3}, DiscreteVariables -> {controlU}];
controlConstant[Ω_, Δ_] := First @ QuantumEvolve[Δ/2 QuantumOperator["Z"] + Ω/2 QuantumOperator["X"], {QuantumOperator["-"]} -> {1},
    QuantumState["1"], {t, 0, 3}]["ProbabilitiesList"];
controlPanel[Ω_, Δr_] := With[{ps = {0.6, 0.75, 0.9}},
    Plot[Evaluate[Append[(1 + controlRun[Ω, Δr Ω, #][t])/2 & /@ ps, controlConstant[Ω, Δr Ω]]], {t, 0, 3},
      PlotStyle -> {Automatic, Automatic, Automatic, Directive[Thick, Dashed]}, PlotPoints -> 400, MaxRecursion -> 0, PlotRange -> {0, 1},
      PlotLegends -> Append[Row[{"p = ", #}] & /@ ps, "constant drive"], Frame -> True, FrameLabel -> {"\[Gamma] t", Subscript["P", "e"]},
      PlotLabel -> Row[{"\[CapitalDelta]/\[CapitalOmega] = ", Δr}], ImageSize -> 300]];
controlPanels = AssociationThread[{5, 15, 40}, Outer[controlPanel, {5, 15, 40}, {0, 0.5}]];
Manipulate[Legended[PlotGrid[{First /@ controlPanels[Ω]}, ImageSize -> 600], legendBelow[Last[controlPanels[5][[1]]]]], Row[{Control[{{Ω, 15, "\[CapitalOmega]/\[Gamma]"}, Keys[controlPanels], SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"Ω/γ = ", Ω}], <|5 -> "Only p = 0.6 is held; p = 0.75 is clipped once; p = 0.9 is never reached.",
        15 -> "p = 0.6 and 0.75 are held; p = 0.9 is clipped once.",
        40 -> "Every band is held, a higher one for a shorter time."|>[Ω]},
      {"Δ/Ω = 0.5", "P_e stays below Ω²/(Ω²+Δ²); only the first peaks are clipped."}}]]}], SaveDefinitions -> True]
```

Without detuning, a drive that overshoots the band is held inside it, above $1/2$, for a time that shortens as the target rises, after which $P_e$ oscillates about $1/2$ as under constant drive; a drive whose peak barely reaches the band is only clipped once, and one that never reaches it gives the constant-drive curve. With $\Delta/\Omega=0.5$ the tilted axis keeps $P_e$ below $\Omega^2/(\Omega^2+\Delta^2)$; the switching clips only the first peaks, then the oscillation decays below every band.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Continuous measurement and stochastic calculus

Wolfram Language has Itô and Stratonovich processes built in (`ItoProcess`, `StratonovichProcess`, `RandomFunction`), so stochastic master equations are written as the equations themselves, and Itô's lemma and `Reduce` prove statements about every trajectory.

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

### A Hamiltonian with white noise: Stratonovich to Itô gives the dephasing


A qubit with Hamiltonian $\frac{\omega}{2}\hat\sigma_x$ also feels a white-noise field along $z$ with power $\gamma$: $\hat H(t)=\frac{\omega}{2}\hat\sigma_x+\sqrt\gamma\,\xi(t)\,\hat\sigma_z$, with $\xi$ white noise. Every run is a unitary rotation; what does the average over runs do?

With $\hbar=1$, white noise as the limit of smooth noise makes the equation of motion a Stratonovich equation,

$$
d\rho=-i\left[\frac{\omega}{2}\hat\sigma_x,\rho\right]dt-i\sqrt\gamma\,{}[\hat\sigma_z,\rho]\circ dW ,
$$

and only its Itô form, whose drift gains $\frac12(b\cdot\nabla)b$ from the noise $b=2\sqrt\gamma\,(-y,x,0)$ of the Bloch vector, gives the average. The Itô drift of $(x,y,z)$ is:

```wl
noisyBlochOf[m_] := Simplify[Tr[m . #] & /@ PauliMatrix[{1, 2, 3}]];
noisyH = With[{ρv = (IdentityMatrix[2] + {x[t], y[t], z[t]} . PauliMatrix[{1, 2, 3}])/2},
   StratonovichProcess[{noisyBlochOf[-I ω/2 (PauliMatrix[1] . ρv - ρv . PauliMatrix[1])],
     List /@ noisyBlochOf[-I Sqrt[γ] (PauliMatrix[3] . ρv - ρv . PauliMatrix[3])]}, {{x, y, z}, {0, 0, 1}}, t]];
TraditionalForm[ItoProcess[noisyH]["Drift"] /. {x[t] -> x, y[t] -> y, z[t] -> z}]
```

200 runs start from $z=1$ as alternating rotations, a random angle about $z$ and $\omega\,dt$ about $x$. Each panel draws three single runs (gray), the average of all 200 (dots) and the prediction of the Lindblad equation $\dot\rho=-i[\frac{\omega}{2}\hat\sigma_x,\rho]+\gamma\,\mathcal D[\hat\sigma_z]\rho$ (red); time is in units of $1/\omega$:

```wl
noisyRunsZ[γ_] := With[{dt = 0.01}, Transpose[FoldList[Function[{v, a}, With[{xr = v[[1]] Cos[a] - v[[2]] Sin[a], yr = v[[1]] Sin[a] + v[[2]] Cos[a]},
       {xr, Cos[dt] yr - Sin[dt] v[[3]], Sin[dt] yr + Cos[dt] v[[3]]}]], {ConstantArray[0., 200], ConstantArray[0., 200], ConstantArray[1., 200]},
     RandomVariate[NormalDistribution[0, 2 Sqrt[γ dt]], {2000, 200}]][[All, 3]]]];
lindbladNoisyH[γ_] := QuantumEvolve[QuantumOperator["X"]/2, {QuantumOperator["Z"]} -> {γ}, QuantumState["0"], {t, 0, 20}]["BlochVector"][[3]];
noisyRunStyle = Directive[GrayLevel[0.6], Opacity[0.6], AbsoluteThickness[1]];
panelNoisyH[γ_] := With[{zs = noisyRunsZ[γ]}, Show[
    ListLinePlot[Take[zs, 3], DataRange -> {0, 20}, PlotRange -> {-1, 1}, PlotStyle -> noisyRunStyle],
    Plot[Evaluate[lindbladNoisyH[γ]], {t, 0, 20}, PlotRange -> {-1, 1}, PlotStyle -> Directive[Red, AbsoluteThickness[2]]],
    ListPlot[Mean[zs][[;; ;; 40]], DataRange -> {0, 20}, PlotMarkers -> {"\[FilledCircle]", 7}, PlotStyle -> ColorData[97][1]],
    PlotRange -> {-1, 1}, Frame -> True, FrameLabel -> {"ω t", Row[{"\[LeftAngleBracket]", Subscript["\[Sigma]", "z"], "\[RightAngleBracket]"}]},
    PlotLabel -> Row[{"γ/ω = ", γ}], ImageSize -> 300]];
Legended[BlockRandom[SeedRandom[3]; PlotGrid[{panelNoisyH /@ {0.1, 0.5, 2}}, ImageSize -> 900]],
  legendBelow[Column[{LineLegend[{noisyRunStyle}, {"single runs"}], PointLegend[{ColorData[97][1]}, {"average of 200 runs"}, LegendMarkers -> {"\[FilledCircle]", 7}],
    LineLegend[{Directive[Red, AbsoluteThickness[2]]}, {"Lindblad"}]}]]]
```

The conversion adds $-2\gamma x$ and $-2\gamma y$ to the drift: dephasing at rate $2\gamma$ about the noise axis. The average over runs follows the Lindblad curve, which relaxes toward $z=0$. The decay is fastest near $\gamma=\omega$; stronger noise slows it again, because the fast phase noise freezes the rotation.

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

## Measurement, tomography and benchmarking

Measurements, POVMs and estimators are objects too: simulated counts go into `QuantumStateEstimate`, measurement dilations go into circuits, and channel averages give benchmarking curves in closed form.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### State tomography from finite counts


An unknown qubit state $\rho=\frac12(I+\vec r\cdot\hat{\vec\sigma})$ is measured on $N$ copies in each of the three Pauli bases, and the state is estimated from the counts. How does the error of the estimate fall with $N$ for a maximally mixed and for a pure state?

In setting $k$ the counts $n_\pm$ give $\langle\hat\sigma_k\rangle\approx(n_+-n_-)/N$, each an average of $N$ outcomes $\pm1$ with variance $1-\langle\hat\sigma_k\rangle^2$, so the linear-inversion estimate has mean squared error $(3-|\vec r|^2)/N$. For a pure state the estimate is held inside the Bloch ball, which removes the outward half of the error along $\vec r$: $\left(2-\frac12(1-\sum_kr_k^4)\right)/N$. The estimator used here inverts the Born rule linearly, makes the result a valid state, and then climbs the likelihood of the counts by random steps, so the estimate is always a state.

The counts come from a true $\vec r$ off every measurement axis, along $(4,3,8.5)$; the figure plots $\sqrt N$ times the root-mean-square error $|\vec r_{\rm est}-\vec r|$ over 100 runs, against these two references (dashed); a flat curve means the error falls as $N^{-1/2}$:

```wl
tomoOps = QuantumMeasurementOperator /@ {"X", "Y", "Z"};
tomoNs = {30, 300, 3000};
tomoErr[len_] := BlockRandom[SeedRandom[7]; With[{r = len Normalize[{4, 3, 8.5}]},
    Table[{n, Sqrt[n Mean @ Table[Norm[QuantumStateEstimate[QuantumMeasurementSimulation[QuantumState["BlochVector"[r]], tomoOps, n]][
          "MaximumLikelihoodState"]["BlochVector"] - r]^2, 100]]}, {n, tomoNs}]]];
tomoLens = {0, 1};
ListPlot[Join[tomoErr /@ tomoLens, Function[len, {#, Sqrt[3 - len^2 - len^2 (1 - Total[Normalize[{4, 3, 8.5}]^4])/2]} & /@ tomoNs] /@ tomoLens], ScalingFunctions -> {"Log", None}, Joined -> True,
 PlotMarkers -> {Automatic, Automatic, None, None}, PlotStyle -> {Automatic, Automatic, Dashed, Dashed}, PlotRange -> {0, 2},
 PlotLegends -> Join[Row[{"|r| = ", #}] & /@ tomoLens, Row[{"projected linear-inversion reference, |r| = ", #}] & /@ tomoLens], Frame -> True,
 FrameLabel -> {"shots per setting N", Row[{Sqrt["N"], " \[Times] rms Bloch-vector error"}]}, ImageSize -> 460]
```

Both curves are flat, so the error falls as $N^{-1/2}$, and each scatters about its reference. The pure state is estimated more accurately at every $N$, partly because its outcomes are less random and partly because the estimate cannot leave the Bloch ball.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Process tomography of a decaying qubit


A qubit relaxes from $|1\rangle$ to $|0\rangle$ with probability $\gamma$. The channel is treated as unknown: four probe states go through it, and each output is reconstructed from $N$ shots per Pauli setting. How does the reconstructed channel approach the true one as $N$ grows?

With $\hat\sigma_z=|0\rangle\langle0|-|1\rangle\langle1|$, the amplitude-damping channel has Kraus operators $\hat K_0=|0\rangle\langle0|+\sqrt{1-\gamma}\,|1\rangle\langle1|$ and $\hat K_1=\sqrt\gamma\,|0\rangle\langle1|$. Any qubit channel acts on Bloch vectors as $\vec r\mapsto T\vec r+\vec c$. The probes $|0\rangle,|1\rangle,|{+}\rangle,|{+}i\rangle$ have Bloch vectors $\hat z,-\hat z,\hat x,\hat y$, so their outputs give

$$
\vec c=\frac12(\vec r_0+\vec r_1),\qquad T=\left(\vec r_+-\vec c,\ \ \vec r_{+i}-\vec c,\ \ \frac12(\vec r_0-\vec r_1)\right)\ \text{(columns)},
$$

with each output $\vec r$ estimated by maximum likelihood from the counts.

The figure draws the image of the input great circle in the $xz$ plane under the reconstructed $(T,\vec c)$ and under the exact channel:

```wl
ptProbes = QuantumState /@ {"0", "1", "Plus", "Right"};
ptOps = QuantumMeasurementOperator /@ {"X", "Y", "Z"};
ptMap[{r0_, r1_, rp_, ri_}] := With[{c = (r0 + r1)/2}, {Transpose[{rp - c, ri - c, (r0 - r1)/2}], c}];
ptCurve[{tm_, c_}] := (tm . {Sin[α], 0, Cos[α]} + c)[[{1, 3}]];
ptOut[γ_, n_] := BlockRandom[SeedRandom[n];
   QuantumStateEstimate[QuantumMeasurementSimulation[QuantumChannel["AmplitudeDamping"[γ]][#], ptOps, n]]["MaximumLikelihoodState"]["BlochVector"] & /@ ptProbes];
ptNs = {30, 300, 3000};
ptPanel[γ_] := With[{exact = ptMap[QuantumChannel["AmplitudeDamping"[γ]][#]["BlochVector"] & /@ ptProbes], rec = ptMap[Re[ptOut[γ, #]]] & /@ ptNs},
   ParametricPlot[Evaluate[Join[ptCurve /@ rec, {ptCurve[exact], {Sin[α], Cos[α]}}]], {α, 0, 2 Pi},
    PlotStyle -> {Automatic, Automatic, Automatic, Dashed, Dotted}, PlotLegends -> Join[Row[{"N = ", #}] & /@ ptNs, {"exact channel", "input great circle"}],
    PlotRange -> {{-1.05, 1.05}, {-1.05, 1.05}}, AspectRatio -> 1, Axes -> False, Frame -> True, FrameLabel -> {"x", "z"}, ImageSize -> 300,
    Epilog -> Inset[Framed[Style[Row[{"γ = ", γ}], 10], Background -> White, FrameStyle -> GrayLevel[0.8], FrameMargins -> 2], Scaled[{0.97, 0.97}], {Right, Top}]]];
ptPanels = ptPanel /@ {0.2, 0.5, 0.8};
Legended[PlotGrid[{First /@ ptPanels}, ImageSize -> 900], legendBelow[Last[First[ptPanels]]]]
```

The exact channel shrinks the circle to an ellipse pulled toward $|0\rangle$ at the top, flatter along $z$ than along $x$, and smaller as $\gamma$ grows. The reconstruction error falls as $N$ grows. The exact image touches the sphere at $|0\rangle$ with the sphere's own curvature, so the reconstruction errors carry part of the image of the sphere outside the Bloch ball, by less as $N$ grows: although every output it was built from is a state, the reconstructed map is not a channel.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A tetrahedron measurement on an unknown qubit


A qubit is in an unknown state $\rho=\frac12(I+\vec r\cdot\hat{\vec\sigma})$, pure when $\vec r=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$, and a single detector with four outcomes measures it. What do the four outcome probabilities reveal about the state?

The detector is the symmetric informationally complete POVM of a qubit: four effects $\hat E_i=\frac12|\psi_i\rangle\langle\psi_i|$ with $\sum_i\hat E_i=I$ and equal overlaps $|\langle\psi_i|\psi_j\rangle|^2=\frac13$ for $i\neq j$, whose Bloch vectors are the vertices of a regular tetrahedron. Here the first vertex is $|\psi_1\rangle=|0\rangle$; the other three are $|\psi_k\rangle=(|0\rangle+\sqrt2\,e^{i\chi_k}|1\rangle)/\sqrt3$ with $(\chi_2,\chi_3,\chi_4)=(\frac{4\pi}3,\frac{2\pi}3,0)$. The Born rule gives $p_i=\operatorname{Tr}[\rho\hat E_i]$.

Applied to a general $\vec r$, the Born rule recovers $\vec r$ from the $p_i$; applied to general $\theta$ and $\phi$, it maps each probability over the sphere of pure states:

```wl
sicE = Normal /@ QuantumMeasurementOperator["TetrahedronSICPOVM"]["POVMElements"];
sicN = QuantumState[#/Tr[#]]["BlochVector"] & /@ sicE;
sicBorn[r_] := FullSimplify[ComplexExpand[Tr[Normal[QuantumState["BlochVector"[r]]["DensityMatrix"]] . #] & /@ sicE]];
sicP = sicBorn[{Sin[θ] Cos[φ], Sin[θ] Sin[φ], Cos[θ]}];
sicPanels = GeoDensityPlot[Rule @@ Transpose[Flatten[Table[{GeoPosition[{90 - θv 180/Pi, φv 180/Pi - 180}], # /. {θ -> θv, φ -> φv}}, {θv, 0, Pi, Pi/40}, {φv, 0, 2 Pi, Pi/40}], 1]],
     GeoProjection -> "Hammer", GeoBackground -> None, OpacityFunction -> None, GeoGridLines -> Automatic, GeoGridLinesStyle -> Gray,
     PlotRange -> {0, 1/2}, PlotLegends -> Automatic, ImageSize -> 300] & /@ sicP;
sicLabel[i_] := Style[Row[{Subscript["p", i], " = (1 + ", OverVector["r"], "\[CenterDot]", Subscript[OverVector["n"], i], ")/4"}], 13, FontFamily -> "Helvetica"];
Labeled[Legended[Grid[Partition[MapIndexed[Labeled[First[#1], sicLabel[First[#2]], Top] &, sicPanels], 2], Spacings -> {2, 2}], Placed[Last[First[sicPanels]], Right]],
 Style[Row[{"for ", OverVector["r"], " = (x, y, z):   ", Subscript["p", 1], " + ", Subscript["p", 2], " + ", Subscript["p", 3], " + ", Subscript["p", 4], " = ", Simplify[Total[sicBorn[{x, y, z}]]],
    ",   3 Σ ", Subscript["p", "i"], Subscript[OverVector["n"], "i"], " = ", FullSimplify[3 sicBorn[{x, y, z}] . sicN]}], 13, FontFamily -> "Helvetica"], Top]
```

Each probability is $\frac14$ plus a linear function of $\vec r$, largest, $\frac12$, at its own vertex and zero at the opposite point, and the four sum to one for every state. Inverting that linear map, $\vec r=3\sum_ip_i\vec n_i$ returns the Bloch vector of any state, mixed or pure: one measurement setting determines the state.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Randomized benchmarking: the decay of sequence fidelity


A qubit starts in $|0\rangle$, undergoes $m$ Clifford gates drawn at random and then the Clifford that undoes their product, and is measured; every gate is followed by the same noise. The noise is amplitude damping, which drops $|1\rangle$ to $|0\rangle$ with probability $\gamma$, or a coherent over-rotation $e^{-i\varepsilon X/2}$ with $\varepsilon$ chosen so that both noises have the same average gate fidelity $F$. How does the survival probability, the chance of finding $|0\rangle$ at the end, decay with $m$?

Any qubit channel moves the Bloch vector as $\vec r\mapsto M\vec r+\vec c$. The 24 Cliffords act on $\vec r$ as the rotations of a cube, and averaging the noise over them, which the random sequence does, leaves a depolarizing map that shrinks $\vec r$ by $p=\operatorname{Tr}M/3=2F-1$. With the final noise left over, the mean survival is exactly

$$
P(m)=\frac12\left(1+c_z+M_{zz}\,p^m\right)=A\,p^m+B .
$$

For damping $M=\mathrm{diag}(\sqrt{1-\gamma},\sqrt{1-\gamma},1-\gamma)$ and $\vec c=(0,0,\gamma)$; for the rotation $\operatorname{Tr}M=1+2\cos\varepsilon$, so equal $F$ means $\cos\varepsilon=\sqrt{1-\gamma}-\gamma/2$.

$M$ and $\vec c$ of each noise come from its action on Bloch states; the survival is averaged over 50 random sequences at each $m$, in three independent runs (markers), against $P(m)$ (dashed):

```wl
cliffordRB = Select[Flatten[Outer[DiagonalMatrix[#1] . #2 &, Tuples[{1, -1}, 3], Permutations[IdentityMatrix[3]], 1], 1], Det[#] == 1 &];
blochMapRB[f_] := With[{c = f[QuantumState["BlochVector"[{0, 0, 0}]]]["BlochVector"]}, {Transpose[f[QuantumState["BlochVector"[#]]]["BlochVector"] - c & /@ IdentityMatrix[3]], c}];
survivalRB[{mm_, c_}, m_] := With[{seq = RandomChoice[cliffordRB, m], step = mm . (#2 . #1) + c &},
   (1 + Last[step[Fold[step, {0, 0, 1}, seq], Transpose[Dot @@ Reverse[seq]]]])/2];
decayRB[{mm_, c_}][m_] := (1 + c[[3]] + mm[[3, 3]] (Tr[mm]/3)^m)/2;
panelRB[noise_, title_] := With[{map = Chop[Re[N[blochMapRB[noise]]]], ms = Range[1, 61, 4]},
   Show[ListPlot[Table[BlockRandom[SeedRandom[run]; {#, Mean[Table[survivalRB[map, #], 50]]} & /@ ms], {run, 3}], PlotMarkers -> {Automatic, 8}, PlotRange -> All],
    Plot[decayRB[map][m], {m, 0, 61}, PlotRange -> All, PlotStyle -> Dashed],
    PlotRange -> {All, {0.4, 1.01}}, Axes -> False, Frame -> True, FrameLabel -> {"m", "survival"}, ImageSize -> 300,
    Epilog -> Inset[Framed[Style[Column[{title, Row[{"p = Tr M/3 = ", NumberForm[Tr[First[map]]/3, 3]}]}, Alignment -> Right], 10], Background -> White, FrameStyle -> GrayLevel[0.8], FrameMargins -> 2],
      Scaled[{0.97, 0.97}], {Right, Top}]]];
\[Gamma]sRB = {1/20, 1/5};
panelsRB = {panelRB[QuantumChannel["AmplitudeDamping"[#]], Row[{"damping, \[Gamma] = ", #}]] & /@ \[Gamma]sRB,
   panelRB[QuantumOperator["RX"[ArcCos[Sqrt[1 - #] - #/2]]], Row[{"over-rotation matched to \[Gamma] = ", #}]] & /@ \[Gamma]sRB};
Legended[PlotGrid[panelsRB, ImageSize -> 700], legendBelow[Column[{PointLegend[ColorData[97] /@ Range[3], {"run 1", "run 2", "run 3"}, LegendMarkers -> Automatic],
   LineLegend[{Directive[ColorData[97][1], Dashed]}, {Row[{"P(m) = A ", Superscript["p", "m"], " + B"}]}]}]]]
```

In each column the mean survival of the two noises falls at the same rate: averaging over random Cliffords keeps only $\operatorname{Tr}M$ of the noise. The floors differ: damping pulls toward $|0\rangle$ and holds long sequences above $1/2$, while the rotation, which has no shift $\vec c$, settles at $1/2$. The rotation's runs scatter more, because within one sequence coherent errors can add up or cancel, and only the average over sequences follows the exponential.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Undoing a measurement while its record survives


A qubit in $|\psi\rangle=\cos\frac\theta2|0\rangle+e^{i\phi}\sin\frac\theta2|1\rangle$ is measured in the $Z$ basis, and the outcome is written on a pointer qubit instead of being read. Before the measurement is reversed, the environment partly decoheres the pointer, at strength $\lambda$. How much of $|\psi\rangle$ does the reversal return?

The measurement is the isometry $\hat V|\psi\rangle=\sum_k\hat P_k|\psi\rangle\otimes|k\rangle_{\rm P}$, $\hat P_k=|k\rangle\langle k|$, with $\hat V^\dagger\hat V=I$, so $\hat V^\dagger$ is its exact left inverse. Between the two, phase damping $\mathcal D_\lambda$ acts on the pointer with Kraus operators $|0\rangle\langle0|+\sqrt{1-\lambda}\,|1\rangle\langle1|$ and $\sqrt\lambda\,|1\rangle\langle1|$, which multiply the record's coherence by $\sqrt{1-\lambda}$: $\lambda=0$ keeps the record intact and $\lambda=1$ erases its coherence. From $\rho=|\psi\rangle\langle\psi|$ the returned state is $\rho_{\rm out}=\hat V^\dagger\,\mathcal D_\lambda(\hat V\rho\hat V^\dagger)\,\hat V$, with fidelity $F=\langle\psi|\rho_{\rm out}|\psi\rangle$.

The measurement's dilation $\hat V$ puts its pointer on wire 2, and measure, damp, un-measure run as one circuit for general $\theta$, $\phi$, $\lambda$; $\rho_{\rm out}$ and $F$ come out in closed form:

```wl
umPsi = QuantumState[{Cos[θ/2], Exp[I φ] Sin[θ/2]}];
umV = QuantumOperator[QuantumMeasurementOperator[{1}]["SuperOperator"], {{2, 1}, {1}}];
umRho = QuantumCircuitOperator[{umV, QuantumChannel["PhaseDamping"[λ], {2}], umV["Dagger"]}][umPsi]["DensityMatrix"];
umRhoOut = FullSimplify[Normal[umRho], {θ, φ} ∈ Reals && 0 <= λ <= 1];
umF = Collect[Expand[FullSimplify[ComplexExpand[Conjugate[Normal[umPsi["StateVector"]]] . umRhoOut . Normal[umPsi["StateVector"]]], {θ, φ} ∈ Reals && 0 <= λ <= 1] /.
     Cos[2 θ] -> 1 - 2 Sin[θ]^2], Sin[θ]^2, Simplify];
Style[Column[{Row[{Subscript["ρ", "out"], " = ", TraditionalForm[MatrixForm[umRhoOut]]}], Row[{"F = ", TraditionalForm[umF]}]}, Center, Spacings -> 1], 14, FontFamily -> "Helvetica"]
```

The fidelity is plotted against the input's polar angle:

```wl
umλs = {0, 0.36, 0.75, 1};
Plot[Evaluate[umF /. λ -> umλs], {θ, 0, Pi}, PlotRange -> {0.45, 1.02}, Frame -> True, FrameLabel -> {"θ", "fidelity F"},
 FrameTicks -> {{Automatic, None}, {{0, Pi/4, Pi/2, 3 Pi/4, Pi}, None}}, PlotLegends -> (Row[{"λ = ", #}] & /@ umλs), ImageSize -> 420]
```

The reversal returns $|\psi\rangle$ with its populations intact and its coherence multiplied by $\sqrt{1-\lambda}$: decohering the record dephases the qubit, so $F$ does not depend on $\phi$ and with the record intact every input returns exactly. As the record decoheres, the poles, which the $Z$ measurement does not disturb, still return intact while the equator loses most; at $\lambda=1$ the equatorial fidelity is $\frac12$, the value left by a completed measurement.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Quantum Fisher information: n for product probes, n² for GHZ


A phase $\phi$ is written on $n$ qubits by $R_Z(\phi)$ on each, and any unbiased estimate of $\phi$ from $\nu$ repetitions has variance at least $1/(\nu F_Q)$. How does the quantum Fisher information $F_Q$ grow with $n$ for the product probe $|+\rangle^{\otimes n}$ and for the GHZ probe $(|0\cdots0\rangle+|1\cdots1\rangle)/\sqrt2$?

For a pure probe the phase family and its Fisher information, four times the Fubini-Study metric, are

$$
|\psi_\phi\rangle=\bigotimes_{k=1}^n e^{-i\phi Z_k/2}\,|\psi_0\rangle,\qquad F_Q=4\left(\langle\partial_\phi\psi_\phi|\partial_\phi\psi_\phi\rangle-\left|\langle\psi_\phi|\partial_\phi\psi_\phi\rangle\right|^2\right)=4\,\mathrm{Var}_{\psi_0}\left(\sum_k\frac{Z_k}{2}\right).
$$

The metric comes from differentiating the state for general $\phi$, for $n=1$ to $6$:

```wl
Needs["Wolfram`QuantumFramework`QuantumOptimization`"];
qfiPhase[probe_] := 4 FubiniStudyMetricTensor[
     QuantumCircuitOperator[{"RZ"[ϕ] -> Range[probe["Qudits"]]}][probe], "Matrix", "Parameters" -> {ϕ}][[1, 1]];
ListPlot[Table[{m, qfiPhase[probe[m]]}, {probe, {QuantumState[StringRepeat["+", #]] &, QuantumState["GHZ"[#]] &}}, {m, 6}],
 PlotLegends -> {"product", "GHZ"}, PlotMarkers -> Automatic, Joined -> True, Frame -> True, FrameLabel -> {"n", Subscript["F", "Q"]},
 PlotRange -> {0, 40}, ImageSize -> 400]
```

The product probe gives $F_Q=n$, the standard quantum limit, since the variances of independent qubits add; the GHZ probe gives $n^2$, the Heisenberg limit, since $\sum_kZ_k/2$ takes only the values $\pm n/2$ on it. At $n=1$ the two probes are the same state.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### GHZ against product probes when every qubit dephases


A frequency $\omega$ is written on $n$ qubits, each dephasing independently at rate $\gamma$, by Ramsey runs of length $t$ repeated over a total time $T$. Without noise a GHZ probe has quantum Fisher information $n^2t^2$ against $nt^2$ for the product probe $|+\rangle^{\otimes n}$; does the GHZ advantage survive once $t$ is chosen optimally?

With $\hbar=1$ and $\mathcal D[\hat c]\rho=\hat c\rho\hat c^\dagger-\frac12\{\hat c^\dagger\hat c,\rho\}$, each run evolves under

$$
\dot\rho=-i\left[\frac\omega2\sum_{k=1}^nZ_k,\rho\right]+\gamma\sum_{k=1}^n\mathcal D[Z_k]\rho ,
$$

which damps each qubit's coherence as $e^{-2\gamma t}$. The $T/t$ runs carry information $(T/t)\,F_Q(t)$, with $F_Q=2\sum_{\lambda_i+\lambda_j>0}|\langle i|\partial_\omega\rho|j\rangle|^2/(\lambda_i+\lambda_j)$ over the eigenpairs $(\lambda_i,|i\rangle)$ of $\rho(t)$; the noise commutes with the phase, so $F_Q$ does not depend on $\omega$ and is taken at $\omega=0$.

The GHZ coherence decays as $e^{-2n\gamma t}$, so $F_Q^{\rm GHZ}=n^2t^2e^{-4n\gamma t}$ and $F_Q^{\rm product}=nt^2e^{-4\gamma t}$, and the best run lengths $t=1/(4n\gamma)$ and $t=1/(4\gamma)$ give

$$
\max_t\frac{F_Q^{\rm GHZ}}{t}=\max_t\frac{F_Q^{\rm product}}{t}=\frac{n}{4e\gamma}.
$$

For both probes and $n=2,3,4$, $F_Q$ is computed from $\rho(t)$ for general $\omega$, $\gamma$ and $t$, and $F_Q/t$ is maximized over $t$ (markers); rates are in units of $\gamma$:

```wl
dgzZ[n_, k_] := QuantumOperator["Z", {k}, n];
dgzRho[ψ_] := With[{n = ψ["Qudits"]}, Simplify[Normal[QuantumEvolve[Total[ω/2 dgzZ[n, #] & /@ Range[n]], dgzZ[n, #] & /@ Range[n] -> ConstantArray[γ, n], ψ, t]["DensityMatrix"]], γ > 0]];
dgzQFI[ρ_] := With[{es = Eigensystem[ρ /. ω -> 0]}, {λ = First[es], V = Orthogonalize[Last[es]]}, {M = Conjugate[V] . (D[ρ, ω] /. ω -> 0) . Transpose[V]},
   Simplify[2 Total[MapThread[If[PossibleZeroQ[#2], 0, Abs[#1]^2/#2] &, {M, Outer[Plus, λ, λ]}, 2], 2], γ > 0 && t > 0]];
dgzF = Table[{dgzQFI[dgzRho[QuantumState["GHZ"[n]]]], dgzQFI[dgzRho[QuantumState[StringRepeat["+", n]]]]}, {n, 2, 4}];
dgzBest = Map[Maximize[{#/t /. γ -> 1, t > 0}, t] &, dgzF, {2}];
dgzForms[k_] := {k^2 t^2 Exp[-4 k γ t], k t^2 Exp[-4 γ t]};
dgzMatch = Simplify[dgzF - Table[dgzForms[k], {k, 2, 4}]] === ConstantArray[0, {3, 2}];
dgzAnyN = Function[f, With[{t0 = t /. First[Solve[D[f/t, t] == 0, t]]}, {Simplify[f/t /. t -> t0], t0}]] /@ (dgzForms[n] /. γ -> 1);
dgzPanel[i_] := Show[Plot[Evaluate[dgzF[[i]]/t /. γ -> 1], {t, 0, 0.6}, PlotRange -> All, PlotStyle -> {Automatic, Dashed}, PlotLegends -> {"GHZ", "product"}],
   ListPlot[List /@ ({t /. Last[#], First[#]} & /@ dgzBest[[i]]), PlotMarkers -> {Automatic, 10}],
   PlotRange -> {{0, 0.6}, {0, 0.5}}, Frame -> True, FrameLabel -> {"γ t", Row[{"γ ", Subscript["F", "Q"], "/t"}]}, Epilog -> Inset[Framed[Style[Row[{"n = ", i + 1}], 10], Background -> White, FrameStyle -> GrayLevel[0.8], FrameMargins -> 2], Scaled[{0.97, 0.97}], {Right, Top}], ImageSize -> 280];
dgzPanels = dgzPanel /@ Range[3];
Labeled[Legended[PlotGrid[{First /@ dgzPanels}, ImageSize -> 840], legendBelow[Last[First[dgzPanels]]]],
  Style[Column[{Grid[Prepend[MapThread[{#1, TraditionalForm[#2[[1]]], TraditionalForm[#3[[1, 1]]], TraditionalForm[t /. #3[[1, 2]]], TraditionalForm[#2[[2]]], TraditionalForm[#3[[2, 1]]], TraditionalForm[t /. #3[[2, 2]]]} &, {Range[2, 4], dgzF, dgzBest}],
       {"n", Row[{"GHZ  ", Subscript["F", "Q"]}], Row[{"largest γ ", Subscript["F", "Q"], "/t"}], "at γ t", Row[{"product  ", Subscript["F", "Q"]}], Row[{"largest γ ", Subscript["F", "Q"], "/t"}], "at γ t"}], Spacings -> {2, 0.6}],
     Row[{"table = ", TraditionalForm[dgzForms[n]], " for n = 2, 3, 4: ", dgzMatch}],
     Row[{"any n:   GHZ largest γ ", Subscript["F", "Q"], "/t = ", TraditionalForm[dgzAnyN[[1, 1]]], " at γ t = ", TraditionalForm[dgzAnyN[[1, 2]]],
       ",   product largest γ ", Subscript["F", "Q"], "/t = ", TraditionalForm[dgzAnyN[[2, 1]]], " at γ t = ", TraditionalForm[dgzAnyN[[2, 2]]]}]}, Center, Spacings -> 1], 13, FontFamily -> "Helvetica"], Bottom]
```

The computed $F_Q$ equals both closed forms for $n=2,3,4$, and each marker sits at $n/4e\gamma$. For short runs the GHZ probe gathers information faster, but the $n^2$ in its $F_Q$ comes with a decay $n$ times faster, so its best run is $n$ times shorter and reaches exactly the same rate as the product probe: under independent dephasing the Heisenberg gain cancels at the optimum.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Circuits and algorithms

Circuits are built from named gates, Boolean functions and named algorithms, run symbolically, and their outputs feed ordinary Wolfram Language (`Maximize`, `MinimalPolynomial`, `ContinuedFraction`).

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### QAOA on MaxCut: the exact landscape and a certified optimum


Color the $n$ vertices of a graph with two colors so that as many edges as possible join different colors; take the triangle, the square and the complete graph $K_4$. One layer of the quantum approximate optimization algorithm (QAOA) prepares a state that depends on two angles $\gamma$ and $\theta$ and scores it by the expected number of cut edges. How close does the best pair of angles come to the true maximum cut?

With $Z_i$ and $X_i$ the Pauli operators on vertex $i$, the cost Hamiltonian, the one-layer state and its landscape are

$$
\hat H_C=\sum_{(i,j)\in\mathcal E}\frac{1-Z_iZ_j}{2},\qquad |\psi(\gamma,\theta)\rangle=\prod_i e^{-i\theta X_i/2}\prod_{(i,j)\in\mathcal E}e^{-i\gamma Z_iZ_j/2}\,|+\rangle^{\otimes n},\qquad E(\gamma,\theta)=\langle\psi|\hat H_C|\psi\rangle,
$$

with $\mathcal E$ the edge set.

$E(\gamma,\theta)$ is in closed form from each graph's edge list, with its exact maximum as a dot; each title compares that maximum, with its minimal polynomial when it is irrational, with the graph's maximum cut:

```wl
qaoaCircuit[graph_] := QuantumCircuitOperator[{"H" -> VertexList[graph], Splice[Thread["R"[γ, "ZZ"] -> List @@@ EdgeList[graph]]], "RX"[θ] -> VertexList[graph]}];
qaoaLandscape[graph_] := With[{hc = Total[(1 - QuantumOperator["ZZ", #])/2 & /@ List @@@ EdgeList[graph]], ψ = qaoaCircuit[graph][]},
    FullSimplify[ψ["Dagger"][hc[ψ]]["Scalar"], {γ, θ} ∈ Reals]];
qaoaPanel[graph_] := With[{l = qaoaLandscape[graph]}, {opt = Maximize[l, {γ, θ}]},
    DensityPlot[l, {γ, 0, 2 π}, {θ, 0, π}, PlotPoints -> 60, PlotLegends -> Automatic, FrameLabel -> {"γ", "θ"},
     Epilog -> {PointSize[Large], Point[Mod[N[{γ, θ} /. Last[opt]], {2 π, π}]]},
     PlotLabel -> Column[{Row[{Graph[graph, ImageSize -> 40], "   max cut = ", First[FindMaximumCut[graph]]}],
        Row[{"max E = ", If[First[opt] ∈ Rationals, First[opt],
           Row[{N[First[opt], 4], ", root of ", TraditionalForm[MinimalPolynomial[First[opt], x] == 0]}]]}]}, Alignment -> Center],
     ImageSize -> 300]];
Column[{Labeled[qaoaCircuit[CycleGraph[3]]["Diagram", ImageSize -> 360], Style["one QAOA layer on the triangle", 13, FontFamily -> "Helvetica"], Top],
  Grid[{qaoaPanel /@ {CycleGraph[3], CycleGraph[4], CompleteGraph[4]}}]}, Center, Spacings -> 2]
```

For the triangle the best angles reach the maximum cut. For the square one layer reaches exactly three quarters of the edges, the ratio Farhi, Goldstone and Gutmann found for long even rings, and for $K_4$ it also falls short. The maximum is exact, not the best point a search visited, so the shortfall is certified.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Teleportation with the classical channel as record wires


Alice holds an unknown qubit and one half of a Bell pair; Bob holds the other half. She measures her two qubits in the Bell basis, and Bob corrects his qubit according to her two outcomes. The outcomes are kept as record qubits, so Bob's corrections are gates controlled by the records. Does Bob end up with Alice's state, for every input?

With the input $|\psi\rangle=\cos\frac\beta2|0\rangle+e^{i\phi}\sin\frac\beta2|1\rangle$ on qubit 1, whose Bloch vector is $(\sin\beta\cos\phi,\sin\beta\sin\phi,\cos\beta)$, and $|00\rangle$ on qubits 2 and 3, the protocol is

$$
\underbrace{\mathrm{CNOT}_{23}\,H_2}_{\text{Bell pair}},\qquad \underbrace{H_1\,\mathrm{CNOT}_{12}}_{\text{Bell basis}},\qquad \text{measure qubits 1, 2}\to m_1,m_2,\qquad \text{Bob applies } Z^{m_1}X^{m_2}.
$$

The protocol runs for general $\beta$ and $\phi$ with both record-controlled corrections, with one of them, and with none. Bob's Bloch vector comes out in closed form; each sphere draws it (dots) and the input's (dashed) along the inputs $\phi=4\beta$, $0\le\beta\le\pi$:

```wl
teleportCircuit[corrections_] := QuantumCircuitOperator[Join[{"H" -> 2, "CNOT" -> {2, 3}, "CNOT" -> {1, 2}, "H" -> 1, {1}, {2}}, corrections]];
bobBloch[corrections_] := FullSimplify[QuantumPartialTrace[Replace[
      teleportCircuit[corrections][
       QuantumTensorProduct[QuantumState[{Cos[β/2], Exp[I ϕ] Sin[β/2]}], QuantumState["00"]]], m_QuantumMeasurement :> m["State"]],
     {1, 2, 3, 4}]["BlochVector"], {β, ϕ} ∈ Reals];
teleportPanel[label_, bob_] := Labeled[Show[ResourceFunction["BlochSpherePlot"][],
   ParametricPlot3D[{Sin[β] Cos[4 β], Sin[β] Sin[4 β], Cos[β]}, {β, 0, π}, PlotStyle -> Dashed],
   ListPointPlot3D[Table[bob /. ϕ -> 4 β, {β, 0, π, π/24}], PlotStyle -> PointSize[0.03]], ImageSize -> 240],
  Style[Column[{label, TraditionalForm[bob]}, Alignment -> Center], 13, FontFamily -> "Helvetica"], Top];
teleportPanels = MapThread[teleportPanel[#1, bobBloch[#2]] &, {{"X and Z", "X only", "Z only", "no correction"},
    {{"CNOT" -> {-1, 3}, "CZ" -> {0, 3}}, {"CNOT" -> {-1, 3}}, {"CZ" -> {0, 3}}, {}}}];
Column[{Labeled[teleportCircuit[{"CNOT" -> {-1, 3}, "CZ" -> {0, 3}}]["Diagram", ImageSize -> 520], Style["the protocol with both record-controlled corrections", 13, FontFamily -> "Helvetica"], Top],
  Legended[Grid[{teleportPanels}], legendBelow[Column[{PointLegend[{ColorData[97][1]}, {"Bob's Bloch vector"}, LegendMarkers -> {"\[FilledCircle]", 10}],
   LineLegend[{Directive[ColorData[97][1], Dashed]}, {"input Bloch vector"}]}]]]}, Center, Spacings -> 2]
```

With both corrections Bob's Bloch vector equals the input's for every $\beta$ and $\phi$. With only one correction, the Pauli Bob did not undo flips two Bloch components at random, so only the component it leaves unflipped survives and his dots collapse onto one axis. With neither, Bob holds the maximally mixed state at the center of the sphere.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Grover search: the success probability, iteration by iteration


A Boolean oracle accepts $M$ of the $N=2^n$ bit strings. Grover's search starts from the uniform superposition and repeats one iteration, a phase flip of the accepted strings followed by a reflection about the uniform state. How does the probability of finding an accepted string change with the number of iterations $k$?

With $f$ the oracle and $|s\rangle=|+\rangle^{\otimes n}$, one iteration, the state after $k$ of them and the success probability are

$$
G=\left(2|s\rangle\langle s|-1\right)O,\qquad O|x\rangle=(-1)^{f(x)}|x\rangle,\qquad |\psi_k\rangle=G^k|s\rangle,\qquad P_k=\sum_{f(x)=1}|\langle x|\psi_k\rangle|^2 .
$$

In the plane of the accepted and rejected uniform states, $|s\rangle$ sits at the angle $a$ from the rejected one, with $\sin a=\sqrt{M/N}$, and $G$ is a rotation by $2a$, so $P_k=\sin^2\left((2k+1)a\right)$.

$G$, built from the Boolean function, acts on the state again and again (markers), up to about $\pi\sqrt N$ iterations, against the closed form (dashed); the accepted strings are a few fixed integers below 16:

```wl
groverSuccess[n_, accepted_] := With[{g = Normal[QuantumCircuitOperator["GroverPhase"[BooleanFunction[Total[2^accepted], n]]]["Matrix"]]},
   Total[Abs[#[[accepted + 1]]]^2] & /@ NestList[g . # &, ConstantArray[2^(-n/2), 2^n], Ceiling[π Sqrt[2^n]]]];
groverAccepted = {{5}, {5, 9}, {5, 9, 12, 14}};
groverPanels = Table[Show[
    ListPlot[groverSuccess[n, #] & /@ groverAccepted, DataRange -> {0, Ceiling[π Sqrt[2^n]]}, PlotMarkers -> Automatic, PlotLegends -> (Row[{"M = ", Length[#]}] & /@ groverAccepted)],
    Plot[Evaluate[With[{a = ArcSin[Sqrt[Length[#]/2^n]]}, Sin[(2 k + 1) a]^2] & /@ groverAccepted], {k, 0, Ceiling[π Sqrt[2^n]]}, PlotStyle -> Dashed, PlotRange -> All],
    PlotRange -> {0, 1.02}, Frame -> True, FrameLabel -> {"iterations k", "success probability"}, PlotLabel -> Row[{"N = ", 2^n}], ImageSize -> 300], {n, {4, 6, 8}}];
Column[{Labeled[QuantumCircuitOperator["GroverPhase"[BooleanFunction[2^5, 4]]]["Diagram", ImageSize -> 520], Style["one Grover iteration, N = 16, marked item 5", 13, FontFamily -> "Helvetica"], Top],
  Legended[PlotGrid[{First /@ groverPanels}, ImageSize -> 900], legendBelow[Last[First[groverPanels]]]]}, Center, Spacings -> 2]
```

The simulated probabilities sit on the closed form, which depends on the accepted strings only through $M/N$. Each rotation by $2a$ carries the state toward the accepted strings and then past them, so the first peak comes near $k\approx\frac\pi4\sqrt{N/M}$: later for a larger search space, earlier for more accepted strings. For $N=16$ and $M=4$ one iteration finds an accepted string with certainty.

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

### Quantum linear solver: HHL on four qubits


Solve $A\boldsymbol x=\boldsymbol b$ for $A=\begin{pmatrix}1&-1/3\\-1/3&1\end{pmatrix}$ and the unit vector $\boldsymbol b=(\cos\beta,\sin\beta)$. The Harrow-Hassidim-Lloyd (HHL) circuit loads $\boldsymbol b$ into one qubit, writes the eigenvalues of $A$ into a two-qubit clock by phase estimation, rotates an ancilla by an angle set by each eigenvalue, clears the clock, and keeps only the runs where the ancilla reads 1. Does the kept qubit give the probabilities $|x_j|^2/|\boldsymbol x|^2$ of the solution?

$A$ has eigenvalues $\lambda_j=2/3,\,4/3$ with eigenvectors $|u_j\rangle$, and $|b\rangle=\sum_jb_j|u_j\rangle$. Phase estimation of $U=e^{iAt}$ with $t=3\pi/4$ writes $\tilde\lambda_j=2^2\lambda_jt/2\pi=1,2$ into the clock, exactly, because the eigenvalue ratio is 2. The controlled rotation $R_Y(\theta_j)$ with $\sin(\theta_j/2)=1/\tilde\lambda_j$, followed by clearing the clock, gives

$$
\sum_jb_j\,|u_j\rangle\,|0_a\rangle\;\longrightarrow\;\sum_jb_j\,|u_j\rangle\left(\sqrt{1-\tilde\lambda_j^{-2}}\,|0_a\rangle+\tilde\lambda_j^{-1}\,|1_a\rangle\right),
$$

and the branch with the ancilla in $|1\rangle$ is proportional to $A^{-1}|b\rangle$.

Each clock value sets one bit, so each rotation needs a single control. The circuit is phase estimation without its initial $X$ and final measurements; the ancilla and solution qubit are measured for general $\beta$ and compared with the direct solution of the linear system:

```wl
hhlA = {{1, -1/3}, {-1/3, 1}};
hhlQPE = QuantumCircuitOperator["PhaseEstimation"[QuantumOperator[MatrixExp[I 3 \[Pi]/4 hhlA]], 2]][[2 ;; -3]];
hhlCircuit = QuantumCircuitOperator[{QuantumState[{Cos[\[Beta]], Sin[\[Beta]]}] -> {4}, hhlQPE -> {2, 3, 4},
    Splice[MapThread["C"["RY"[#1]] -> {#2, 1} &, {2 ArcSin[1/{1, 2}], {3, 2}}]], hhlQPE["Dagger"] -> {2, 3, 4}, {1}, {4}}];
hhlP = FullSimplify[ComplexExpand[hhlCircuit[]["ProbabilitiesList"]], \[Beta] \[Element] Reals];
hhlExact = Last[Normalize[LinearSolve[hhlA, {Cos[\[Beta]], Sin[\[Beta]]}]]]^2;
{hhlKept, hhlCond} = {#[[3]] + #[[4]], #[[4]]/(#[[3]] + #[[4]])} &[hhlP];
Column[{Labeled[hhlCircuit["Diagram", ImageSize -> 600], Style["HHL: phase estimation, controlled rotation of the ancilla, uncompute", 13, FontFamily -> "Helvetica"], Top],
  Plot[Evaluate[{hhlExact, hhlCond, hhlKept}], {\[Beta], 0, \[Pi]}, PlotStyle -> {Thick, Dashed, Automatic}, PlotRange -> {0, 1.05},
 Frame -> True, FrameLabel -> {"\[Beta]", "probability"}, FrameTicks -> {{Automatic, None}, {Range[0, \[Pi], \[Pi]/4], None}},
 PlotLegends -> {Row[{"LinearSolve:  ", Subsuperscript["x", 2, 2], "/|x|\.b2"}], "HHL:  P(solution qubit 1 | ancilla 1)", "HHL:  P(ancilla 1)"},
 PlotLabel -> Style[Row[{"HHL - LinearSolve = ", FullSimplify[hhlCond - hhlExact, \[Beta] \[Element] Reals],
     ",    P(ancilla 1) = ", TraditionalForm[TrigReduce[hhlKept]]}], 13, FontFamily -> "Helvetica"], ImageSize -> 480]}, Center, Spacings -> 2]
```

The HHL probability equals the direct solution identically in $\beta$: the kept qubit carries $|x_2|^2/|\boldsymbol x|^2$. The price is the post-selection: a run is always kept when $\boldsymbol b$ lies along the eigenvector of the smaller eigenvalue, whose rotation turns the ancilla fully to $|1\rangle$, and only a quarter of the time along the other. The identity holds because the two-qubit clock represents both eigenvalues exactly.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Variational ground state: plain against natural gradient descent


A variational quantum eigensolver looks for the ground state of $\hat H=X$ by adjusting the two angles of a one-qubit circuit to lower $\langle\hat H\rangle$. Plain gradient descent steps against the gradient; quantum natural gradient descent first rescales it by the geometry of the state space. From different starting angles and learning rates $\eta$, which optimizer reaches the ground energy in fewer steps?

The circuit, its energy and the update rules are

$$
|\psi(\theta,\phi)\rangle=P(2\phi)\,R_Y(2\theta)|0\rangle=\cos\theta\,|0\rangle+e^{2i\phi}\sin\theta\,|1\rangle,\qquad E(\theta,\phi)=\langle\psi|X|\psi\rangle=\sin2\theta\cos2\phi,
$$

$$
\boldsymbol\theta_{k+1}=\boldsymbol\theta_k-\eta\,\nabla E\quad\text{(plain)},\qquad \boldsymbol\theta_{k+1}=\boldsymbol\theta_k-\eta\,g^{-1}\nabla E\quad\text{(natural)},
$$

with the Fubini-Study metric $g_{ij}=\mathrm{Re}\left(\langle\partial_i\psi|\partial_j\psi\rangle-\langle\partial_i\psi|\psi\rangle\langle\psi|\partial_j\psi\rangle\right)=\mathrm{diag}(1,\sin^22\theta)$. It vanishes in the $\phi$ direction at the poles $\theta=0,\pi/2$, where $\phi$ no longer changes the state.

The metric $g$ is that of the circuit's state; both optimizers run 40 steps on the exact gradient:

```wl
Needs["Wolfram`QuantumFramework`QuantumOptimization`"];
vqeState = QuantumCircuitOperator[{"RY"[2 \[FormalTheta]], "P"[2 \[FormalPhi]]}][];
vqeEnergy[th_, ph_] = ComplexExpand[vqeState["Dagger"][QuantumOperator["X"][vqeState]]["Scalar"]] /. {\[FormalTheta] -> th, \[FormalPhi] -> ph};
vqeGrad = Function[{th, ph}, Evaluate[D[vqeEnergy[th, ph], {{th, ph}}]]];
vqeMetric = FubiniStudyMetricTensor[vqeState, "Parameters" -> {\[FormalTheta], \[FormalPhi]}];
vqePanel[\[Eta]_, start_] := ListLinePlot[{
    vqeEnergy @@@ Re[QuantumNaturalGradientDescent[vqeEnergy, vqeMetric, "InitialPoint" -> start, "LearningRate" -> \[Eta], "MaxIterations" -> 40, "Gradient" -> vqeGrad]],
    vqeEnergy @@@ GradientDescent[vqeEnergy, start, "LearningRate" -> \[Eta], "MaxIterations" -> 40, "Gradient" -> vqeGrad]}, DataRange -> {0, 40},
   PlotRange -> {-1.1, 1.1}, PlotLegends -> {"natural gradient", "plain gradient"}, Axes -> False, Frame -> True, FrameLabel -> {"step", "E"},
   PlotLabel -> Row[{"\[Eta] = ", \[Eta], ",   (\[Theta], \[Phi]) = (", Row[start, ", "], ")"}], ImageSize -> 280];
vqePanels = Outer[vqePanel, {0.05, 0.2}, {{\[Pi]/12, \[Pi]/12}, {\[Pi]/3, \[Pi]/12}, {0.01, 1}}, 1];
Legended[Grid[Map[First, vqePanels, {2}]], legendBelow[Last[vqePanels[[1, 1]]]]]
```

At the smaller rate the natural gradient reaches the ground energy in fewer steps from every start. At the larger rate, from $(\pi/12,\pi/12)$ and $(0.01,1)$, its rescaled steps overshoot and the energy climbs back up before it falls, so plain descent gets there sooner. From $(\pi/3,\pi/12)$ plain descent, at either rate, stalls at $E=0$ at $(\theta,\phi)=(\pi/2,\pi/4)$, a point on the pole $\theta=\pi/2$ where $\partial_\theta E=-2\cos2\phi$ also vanishes, a saddle although $|1\rangle$ is not an eigenstate of $X$; the natural gradient follows the steepest descent on the states themselves and reaches the ground energy.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Phase estimation: how many qubits resolve two eigenphases


A one-qubit unitary has the eigenphases $\varphi_1$ and $\varphi_2=\varphi_1+\delta$, and phase estimation is run with $t$ control qubits on a target that is an equal superposition of the two eigenvectors. How many control qubits separate the two phases in the measured distribution?

With $U|u_k\rangle=e^{2\pi i\varphi_k}|u_k\rangle$, the controls in $|+\rangle^{\otimes t}$ pick up the phases of $U^{2^{t-1}},\dots,U^2,U$ and an inverse Fourier transform reads them out. Each eigenvector contributes with its weight, so the outcome $x$ appears with probability

$$
P(x)=\sum_{k=1,2}\frac{\left|c_k\right|^2}{4^t}\left|\sum_{y=0}^{2^t-1}\exp\left(2\pi i\,y\left(\varphi_k-\frac{x}{2^t}\right)\right)\right|^2,\qquad \left|c_1\right|^2=\left|c_2\right|^2=\frac12 ,
$$

a peak of width about $2^{-t}$ around each $\varphi_k$.

Phase estimation runs on $U=\mathrm{diag}(e^{2\pi i\varphi_1},e^{2\pi i\varphi_2})$ with $\varphi_1=0.17$ and the target in $|+\rangle$; each row is the distribution for one number of control qubits, $t=3,\dots,9$, with bars of width $2^{-t}$, and dashed lines mark the two eigenphases:

```wl
qpePhi1 = 0.17;
qpeProbs[t_, δ_] := N[QuantumCircuitOperator["PhaseEstimation"[QuantumOperator[DiagonalMatrix[Exp[2 π I {qpePhi1, qpePhi1 + δ}]]], t]]][
    QuantumState[StringRepeat["0", t] <> "+"]]["ProbabilitiesList"];
qpePanel[δ_] := Graphics[{Table[With[{p = qpeProbs[t, δ]}, {ColorData[97][1], MapIndexed[Rectangle[{(First[#2] - 1)/2^t, t}, {First[#2]/2^t, t + 0.9 #1}] &, p]}], {t, 3, 9}],
    {Dashed, Gray, Line[{{#, 2.8}, {#, 10}}] & /@ {qpePhi1, qpePhi1 + δ}}}, PlotRange -> {{0, 0.5}, {2.8, 10}}, AspectRatio -> 0.9, Frame -> True,
   FrameLabel -> {Row[{"measured phase x / ", Superscript[2, "t"]}], "control qubits t"}, FrameTicks -> {{Range[3, 9], None}, {Automatic, None}},
   PlotLabel -> Row[{"\[Delta] = ", δ}], ImageSize -> 330];
Column[{Labeled[QuantumCircuitOperator["PhaseEstimation"[QuantumOperator[DiagonalMatrix[Exp[2 π I {qpePhi1, qpePhi1 + 0.1}]]], 3]]["Diagram", ImageSize -> 480], Style["phase estimation with t = 3 control qubits", 13, FontFamily -> "Helvetica"], Top],
  Legended[Grid[{qpePanel /@ {0.1, 0.02}}], legendBelow[Column[{SwatchLegend[{ColorData[97][1]}, {"P(x): bar height in each row"}], LineLegend[{Directive[Gray, Dashed]}, {Row[{"eigenphases ", Subscript["\[CurlyPhi]", 1], ", ", Subscript["\[CurlyPhi]", 2]}]}]}]]]}, Center, Spacings -> 2]
```

While $2^{-t}$ exceeds $\delta$ both phases fall on the same few bars; once $2^{-t}$ is a few times smaller than $\delta$ each phase has its own peak, so a smaller separation needs more control qubits. A phase that falls between grid points $x/2^t$ spreads its probability over the nearest bars.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### One oracle query: Deutsch-Jozsa and Bernstein-Vazirani


An oracle evaluates a Boolean function $f$ of $n$ bits. Deutsch-Jozsa is promised that $f$ is constant or balanced (zero on exactly half of the inputs) and must say which; Bernstein-Vazirani is promised $f(x)=s\cdot x \bmod 2$ and must find the string $s$. Does a single query settle each problem?

Both algorithms query the phase oracle $O_f|x\rangle=(-1)^{f(x)}|x\rangle$ once on $|+\rangle^{\otimes n}$ and measure each qubit in the $X$ basis, which is $H^{\otimes n}$ followed by a computational-basis measurement. The Bernstein-Vazirani circuit here realizes the oracle with CNOTs onto an extra qubit in $|-\rangle$. The outcome $z$ appears with probability

$$
P(z)=\left|\frac{1}{2^n}\sum_{x}(-1)^{f(x)+z\cdot x}\right|^2,\qquad P(0)=\left(1-\frac{2w}{2^n}\right)^2,
$$

where $w$ is the number of inputs $x$ with $f(x)=1$, and for $f(x)=s\cdot x$ the sum is $2^n\delta_{z,s}$.

Deutsch-Jozsa runs on up to four random three-bit functions of each weight $w$ (points, offset sideways; the dashed curve is $P(0)$ above), and Bernstein-Vazirani on all sixteen four-bit secrets:

```wl
djN = 3;
djZeroProb[k_] := First[QuantumCircuitOperator["DeutschJozsaPhase"[BooleanFunction[k, djN]]][]["ProbabilitiesList"]];
djPoints = BlockRandom[SeedRandom[3]; Catenate[KeyValueMap[Function[{w, ks}, With[{sample = RandomSample[ks, UpTo[4]]},
       MapIndexed[{w + 0.12 (First[#2] - (Length[sample] + 1)/2), djZeroProb[#1]} &, sample]]],
     KeySort[GroupBy[Range[0, 2^2^djN - 1], DigitCount[#, 2, 1] &]]]]];
bvTable = Table[QuantumCircuitOperator["BernsteinVazirani"[IntegerString[s, 2, 4]]][]["ProbabilitiesList"], {s, 0, 15}];
Column[{Grid[{{Labeled[QuantumCircuitOperator["DeutschJozsaPhase"[BooleanFunction[15, djN]]]["Diagram", ImageSize -> 300], Style["Deutsch-Jozsa, one balanced f", 13, FontFamily -> "Helvetica"], Top],
    Labeled[QuantumCircuitOperator["BernsteinVazirani"["1011"]]["Diagram", ImageSize -> 360], Style["Bernstein-Vazirani, s = 1011", 13, FontFamily -> "Helvetica"], Top]}}, Spacings -> 4],
  Grid[{{Show[ListPlot[djPoints, PlotMarkers -> Automatic], Plot[(1 - 2 w/2^djN)^2, {w, 0, 2^djN}, PlotStyle -> Dashed], PlotRange -> {{-0.3, 2^djN + 0.3}, {-0.05, 1.05}},
    Frame -> True, FrameTicks -> {{Automatic, None}, {Range[0, 2^djN], None}}, FrameLabel -> {"weight w", Row[{"P(z = ", StringRepeat["0", djN], ")"}]}, PlotLabel -> Row[{"Deutsch-Jozsa, n = ", djN}], ImageSize -> 360],
   ArrayPlot[bvTable, Mesh -> All, Frame -> True, FrameLabel -> {"secret s", "outcome z"},
    FrameTicks -> {{Thread[{Range[1, 16, 3], Range[0, 15, 3]}], None}, {Thread[{Range[1, 16, 3], Range[0, 15, 3]}], None}},
    PlotLabel -> "Bernstein-Vazirani, n = 4: P(z | s), black = 1, white = 0", ImageSize -> 360]}}, Alignment -> Center]}, Center, Spacings -> 2]
```

The constant functions return all zeros with certainty and the balanced ones never do, so under the promise one query decides. Every other function lands on the curve in between, where one query cannot decide. For every secret, the Bernstein-Vazirani outcome is $s$ itself with probability 1.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Search on one qudit with no entangling gate


A single $d$-level system, basis $|0\rangle,\dots,|d-1\rangle$, hides one marked level $|w\rangle$. An oracle flips the sign of $|w\rangle$, and a diffusion step reflects the state about the uniform superposition $|s\rangle$. How many oracle calls does it take to find $w$?

Starting from $|s\rangle=F|0\rangle$, with $F$ the $d$-point Fourier transform, each call applies

$$
G=D\,O,\qquad O=\mathbb 1-2|w\rangle\langle w|,\qquad D=2|s\rangle\langle s|-\mathbb 1,\qquad P_k=\left|\langle w|G^k|s\rangle\right|^2 ,
$$

and $G$ rotates the state in the plane of $|s\rangle$ and $|w\rangle$ by $2\arcsin(1/\sqrt d)\approx2/\sqrt d$, so $P_k=\sin^2\!\left((2k+1)\arcsin(1/\sqrt d)\right)$, which approaches $\sin^2(2k/\sqrt d)$ (dashed) as $d$ grows, and reaching $|w\rangle$ takes about $(\pi/4)\sqrt d$ calls.

$G$ acts repeatedly with $w=2$ for four qudit sizes, and $P_k$ is plotted against $k/\sqrt d$:

```wl
quditSearch[d_, w_ : 2] := With[{basis = QuantumBasis[{d}]}, {step = QuantumOperator[2 ConstantArray[1./d, {d, d}] - IdentityMatrix[d], {1}, basis] @ QuantumOperator[DiagonalMatrix[1. - 2. UnitVector[d, w + 1]], {1}, basis]},
   MapIndexed[{(First[#2] - 1)/Sqrt[d], #1["ProbabilitiesList"][[w + 1]]} &,
    NestList[step, QuantumOperator["Fourier"[d]][QuantumState[UnitVector[d, 1], basis]], Floor[1.6 Sqrt[d]]]]];
quditSizes = {4, 16, 64, 256};
Show[ListLinePlot[quditSearch /@ quditSizes, PlotMarkers -> Automatic, PlotLegends -> (Row[{"d = ", #}] & /@ quditSizes)],
 Plot[Sin[2 u]^2, {u, 0, 1.6}, PlotStyle -> Dashed], PlotRange -> {{-0.05, 1.65}, {-0.03, 1.04}}, Frame -> True, FrameLabel -> {Row[{"k / ", Sqrt["d"]}], "P(w)"}, ImageSize -> 420]
```

Every gate here acts on one $d$-level particle, so there is no second subsystem and no entanglement at any step, yet the marked level is found after a number of calls that grows as $\sqrt d$: in the rescaled time the curves close in on the dashed limit as $d$ grows. The smallest qudit sits off the common curve because it starts a third of the way to $|w\rangle$ and turns by $\pi/3$ per call, so one call lands it on $|w\rangle$ with certainty.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### SU(2) double cover from a controlled rotation


A rotation by $2\pi$ returns every Bloch vector to itself but multiplies the spinor by $-1$, so the rotated qubit alone cannot reveal it. An ancilla in superposition decides whether a target qubit, in an arbitrary pure state, is rotated by the angle $\alpha$ about an arbitrary axis $\mathbf n$. How often does the ancilla read 1 as a function of $\alpha$?

With the target $|\psi\rangle=\cos\frac\beta2|0\rangle+e^{i\phi}\sin\frac\beta2|1\rangle$, the axis $\mathbf n=(\sin\chi\cos\eta,\sin\chi\sin\eta,\cos\chi)$ and $R(\alpha)=e^{-i\alpha\,\mathbf n\cdot\hat{\boldsymbol\sigma}/2}$, the ancilla starts in $|0\rangle$ and the circuit is

$$
(H\otimes 1)\left(|0\rangle\langle0|\otimes 1+|1\rangle\langle1|\otimes R(\alpha)\right)(H\otimes 1)\;|0\rangle|\psi\rangle,
$$

followed by a measurement of the ancilla.

For general $\alpha$, $\beta$, $\phi$, $\chi$ and $\eta$, the probability of reading 1 comes out in closed form:

```wl
su2Circuit = QuantumCircuitOperator[{"H" -> 1, "C"["R"[α, {Sin[χ] Cos[η], Sin[χ] Sin[η], Cos[χ]} . (QuantumOperator /@ {"X", "Y", "Z"})] -> 2, {1}], "H" -> 1, {1}}];
ancillaProbs = Simplify[ComplexExpand[su2Circuit[QuantumTensorProduct[QuantumState["0"], QuantumState[{Cos[β/2], Exp[I ϕ] Sin[β/2]}]]]["ProbabilitiesList"]], 0 < α < 4 π];
Column[{Labeled[su2Circuit["Diagram", ImageSize -> 260], Style["the interference circuit", 13, FontFamily -> "Helvetica"], Top],
  Plot[Last[ancillaProbs], {α, 0, 4 π}, Frame -> True,
 FrameLabel -> {"α", "P(ancilla reads 1)"}, FrameTicks -> {{Automatic, None}, {Range[0, 4 π, π], None}},
 PlotLabel -> Row[{"P(1) = ", TraditionalForm[Last[ancillaProbs]]}], ImageSize -> 450]}, Center, Spacings -> 2]
```

This probability has no trace of the state or the axis: the ancilla sees only the real part of $\langle\psi|R(\alpha)|\psi\rangle=\cos\frac\alpha2-i\sin\frac\alpha2\langle\psi|\mathbf n\cdot\hat{\boldsymbol\sigma}|\psi\rangle$, and the state-dependent part is imaginary. At $\alpha=2\pi$, where every Bloch vector is back where it started, the $-1$ of the spinor makes the ancilla read 1 with certainty; it returns to 0 only at $4\pi$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Compiling a rotation from H and T gates


A fault-tolerant qubit offers only a finite gate set: the Clifford gates generated by $H$ and $S=T^2$, which are cheap, and the $T$ gate, which is expensive. Any target rotation $U$ can only be approximated by a product of them. How does the best approximation improve as the number of $T$ gates grows?

With $T=\mathrm{diag}(1,e^{i\pi/4})$ and $S=T^2$, every product of these gates equals, up to a global phase, a normal form (Matsumoto and Amano)

$$
W=\left(T\ \text{or}\ 1\right)\,s_1s_2\cdots s_k\,C,\qquad s_i\in\{HT,\ SHT\},\qquad C\in\text{the 24 Cliffords},
$$

whose $T$ count is $k$ or $k+1$, so from $T$ count 1 on, each extra $T$ gate doubles the number of candidates. The error of $W$ is $\epsilon=1-\frac12|\mathrm{Tr}\,W^\dagger U|$, which ignores the global phase.

Every normal form up to $T$ count 16 is enumerated and the smallest error kept, for three rotations $R_{\hat n}(\alpha)=e^{-i\alpha\,\hat n\cdot\boldsymbol\sigma/2}$ ($R_{x+z}$ rotates about $(\hat x+\hat z)/\sqrt2$); the circuit that gives the $R_z(0.3)$ point at $T$ count 7 is drawn in time order:

```wl
{ctH, ctS, ctT} = N[{HadamardMatrix[2], DiagonalMatrix[{1, I}], DiagonalMatrix[{1, Exp[I π/4]}]}];
ctMatrix[w_] := Fold[#2 . #1 &, IdentityMatrix[2], w /. {"H" -> ctH, "S" -> ctS}];
ctCliffords = DeleteDuplicatesBy[Catenate[Table[Tuples[{"H", "S"}, l], {l, 0, 6}]],
   With[{f = Flatten[ctMatrix[#]]}, Round[f Exp[-I Arg[SelectFirst[f, Abs[#] > 0.1 &]]], 10.^-6]] &];
ctLevels = NestList[Join[{#[[1]] . ctH . ctT, Join[{"T", "H"}, #[[2]]]} & /@ #, {#[[1]] . ctS . ctH . ctT, Join[{"T", "H", "S"}, #[[2]]]} & /@ #] &, {{IdentityMatrix[2], {}}}, 16];
ctForms[t_] := If[t == 0, ctLevels[[1]], Join[ctLevels[[t + 1]], {ctT . #[[1]], Append[#[[2]], "T"]} & /@ ctLevels[[t]]]];
ctOverlaps[t_, u_] := Abs[Flatten[Conjugate[ctForms[t][[All, 1]]], {{1}, {2, 3}}] . Transpose[Flatten[u . ConjugateTranspose[ctMatrix[#]]] & /@ ctCliffords]]/2;
ctTargets = <|Row[{Subscript["R", "x+z"], "(0.8)"}] -> MatrixExp[-0.4 I (PauliMatrix[1] + PauliMatrix[3])/Sqrt[2]],
   Row[{Subscript["R", "z"], "(0.3)"}] -> MatrixExp[-0.15 I PauliMatrix[3]], Row[{Subscript["R", "x"], "(1)"}] -> MatrixExp[-0.5 I PauliMatrix[1]]|>;
ctBestWord[t_, u_] := With[{pos = First[Position[#, Max[#]]] &[ctOverlaps[t, u]]}, Join[ctCliffords[[pos[[2]]]], ctForms[t][[pos[[1]], 2]]]];
Column[{ListLogPlot[1 - FoldList[Max, Table[Max[ctOverlaps[t, #]], {t, 0, 16}]] & /@ Values[ctTargets], DataRange -> {0, 16}, Joined -> True,
    PlotMarkers -> Automatic, PlotLegends -> Keys[ctTargets], Frame -> True, FrameLabel -> {"T count", "best error ε"}, PlotLabel -> Row[{"normal forms over ", Length[ctCliffords], " Cliffords"}], ImageSize -> 450],
  QuantumCircuitOperator[ctBestWord[7, ctTargets[[2]]]]["Diagram", ImageSize -> 650]}, Alignment -> Center]
```

The best error falls in steps: a target can sit far from every circuit of one $T$ count and close to one of the next. $R_{x+z}(0.8)$ happens to lie close to a circuit with three $T$ gates that no circuit up to $T$ count 16 improves on.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Stabilizers and error correction

`PauliStabilizer` runs Clifford circuits at a thousand qubits, keeps random outcomes as formal bits, and compresses non-Clifford circuits. The QEC examples add exact logical error rates; two of them use the QECCore package that sits beside QF.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Five-qubit code: fifteen single-qubit errors, fifteen distinct syndromes


Five qubits store one logical qubit in the joint $+1$ eigenspace of four commuting checks. One or two qubits suffer an $X$, $Y$ or $Z$ error, and the four checks are measured. Which errors do the four outcomes tell apart?

The checks are the cyclic shifts of $XZZXI$, and the logical operators are

$$
g_1=XZZXI,\quad g_2=IXZZX,\quad g_3=XIXZZ,\quad g_4=ZXIXZ,\qquad \bar Z=ZZZZZ,\quad \bar X=XXXXX .
$$

After an error $E$ on a codeword $|\psi\rangle$, check $g_k$ reads the bit $s_k=\left(1-\langle\psi|E^\dagger g_kE|\psi\rangle\right)/2$, which is 1 exactly when $E$ anticommutes with $g_k$, whatever the codeword.

With $|0_L\rangle$ held as the tableau of $g_1,\dots,g_4,\bar Z$, each of the 15 single-qubit errors and each of the 90 two-qubit errors acts on it; the figure shows the four bits (black is 1) and counts the two-qubit errors that share each single-qubit syndrome:

```wl
fiveCode = PauliStabilizer["5QubitCode"]; fiveChecks = Most[fiveCode["Stabilizers"]];
fiveSyndrome[errors_] := (1 - fiveCode[errors]["Expectation", #] & /@ fiveChecks)/2;
fiveErrors = Flatten[Outer[Rule, {"X", "Y", "Z"}, Range[5]]];
fiveSingles = fiveSyndrome[{#}] & /@ fiveErrors;
fivePairs = Counts[fiveSyndrome /@ Flatten[Outer[Thread[#2 -> #1] &, Subsets[Range[5], {2}], Tuples[{"X", "Y", "Z"}, 2], 1], 1]];
fiveLabels = First[#] <> ToString[Last[#]] & /@ fiveErrors;
Column[{Labeled[ArrayPlot[Transpose[fiveSingles], Mesh -> All, Frame -> True, ImageSize -> 520,
    FrameTicks -> {{Transpose[{Range[4], fiveChecks}], None}, {Transpose[{Range[15], fiveLabels}], None}}],
   Style[Row[{"single-qubit errors: ", Length[DeleteCases[Union[fiveSingles], {0, 0, 0, 0}]], " distinct nonzero syndromes"}], 13, FontFamily -> "Helvetica"], Top],
  BarChart[Lookup[fivePairs, fiveSingles, 0], ChartLabels -> fiveLabels, Frame -> True, FrameLabel -> {None, "two-qubit errors"}, AspectRatio -> 1/4, ImageSize -> 520,
   PlotLabel -> Style[Row[{"two-qubit errors with the same syndrome; with zero syndrome: ", Lookup[fivePairs, Key[{0, 0, 0, 0}], 0]}], 13, FontFamily -> "Helvetica"]]}, Spacings -> 2]
```

Every single-qubit error gives a different nonzero four-bit word, and four bits have only fifteen nonzero values: the checks have no room to spare, which is why the code is called perfect. Every two-qubit error also trips the checks, but it shares its syndrome with a single-qubit error, so a decoder that inverts the table leaves a nontrivial error of weight at most three with zero syndrome, which can only be a logical operator. The code detects all two-qubit errors and corrects none of them.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A thousand-qubit GHZ state: entanglement entropy as a rank over F₂


Take $n=1000$ qubits in the GHZ state $(|0\cdots0\rangle+|1\cdots1\rangle)/\sqrt2$, and for comparison in the one-dimensional cluster state, made by putting every qubit in $|+\rangle$ and applying $CZ$ to each neighboring pair. How much entanglement separates a set $A$ of $k$ qubits from the rest, when $A$ is a block at the end of the chain and when $A$ is scattered at random?

A stabilizer state is fixed by $n$ commuting Pauli generators: $X^{\otimes n}$ and $Z_iZ_{i+1}$ for GHZ, $Z_{i-1}X_iZ_{i+1}$ for the cluster. Write each generator as a binary vector $(x|z)\in\mathbb F_2^{2n}$, stack them as the rows of an $n\times2n$ matrix $G$, and let $G_A$ keep only the $x$ and $z$ columns of the qubits in $A$. Then the entanglement entropy, in bits, is (Fattal et al. 2004)

$$
S(A)=\operatorname{rank}_{\mathbb F_2}G_A-|A|.
$$

Both states are held as tableaux, and the rank is evaluated for $k=50,100,\dots,950$, with one random set $A$ per $k$; no $2^{1000}$-entry vector is ever formed:

```wl
ghzChain = PauliStabilizer[1000][Prepend[Table["CNOT" -> {q, q + 1}, {q, 999}], "H"]];
clusterChain = PauliStabilizer[1000][Join[Thread["H" -> Range[1000]], Table["CZ" -> {q, q + 1}, {q, 999}]]];
entropyPanel[label_, ps_] := ListLogPlot[Table[{k, ps["Entropy", sel[k]]}, {sel, {Range, RandomSample[Range[1000], #] &}}, {k, 50, 950, 50}],
   Joined -> True, PlotMarkers -> {{"\[EmptyCircle]", 16}, {"\[EmptyUpTriangle]", 8}}, PlotStyle -> {Automatic, Dashed}, PlotRange -> {0.7, 500},
   PlotLegends -> {"block {1, \[Ellipsis], k}", "k qubits at random"}, Frame -> True, FrameLabel -> {"k", "S(A)  (bits)"}, PlotLabel -> label, ImageSize -> 340];
entropyPanels = BlockRandom[SeedRandom[1]; {entropyPanel["GHZ", ghzChain], entropyPanel["cluster", clusterChain]}];
Legended[PlotGrid[{First /@ entropyPanels}, ImageSize -> 680], legendBelow[Last[First[entropyPanels]]]]
```

The GHZ state carries exactly one bit across every cut, a block or a scattered set of any size: its entanglement is a single shared parity. The cluster state also carries one bit across a block boundary, but far more when the same number of qubits is scattered, so its entanglement depends on the shape of the cut, not only on its size.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Stabilizer measurement kept symbolic: random outcomes as formal bits


Two qubits start in $|+\rangle$, and two CNOTs copy their parity onto a third qubit. Then an $X$ on qubit 1, an $X$ on qubit 3, or a $Z$ on qubit 1 may strike, and all three qubits are measured in the $Z$ basis. How is the third outcome fixed by the first two?

Before the error the state and its stabilizer generators are

$$
\mathrm{CNOT}_{23}\,\mathrm{CNOT}_{13}\,H_1H_2|000\rangle=\frac12\sum_{m_1,m_2}|m_1,m_2,m_1\veebar m_2\rangle,\qquad \langle X_1X_3,\;X_2X_3,\;Z_1Z_2Z_3\rangle .
$$

A measurement whose outcome is random does not branch: it allocates a formal bit $s_k$, which enters its generator as the sign $(-1)^{s_k}=1-2s_k$. A measurement whose outcome is already determined is stored as a Boolean function of the earlier bits.

Each state is measured once with formal bits; each run reports how many formal bits were allocated and the stored third outcome $m_3$, and the bars count 1000 shots drawn by substituting random values of $s_1,s_2$:

```wl
parityState = PauliStabilizer[3][{"H" -> 1, "H" -> 2, "CNOT" -> {1, 3}, "CNOT" -> {2, 3}}];
parityPanel[label_, ps_] := With[{ms = ps["SymbolicMeasure", {1, 2, 3}]},
   {s = Union[Cases[ms["Phase"], _\[FormalS], Infinity]], m3 = First[Values[ms["Outcomes"]]]},
    BarChart[Lookup[Counts[First[Keys[#["M", {1, 2, 3}]]] & /@ ms["SampleOutcomes", 1000]], Tuples[{0, 1}, 3], 0],
     ChartLabels -> (StringJoin[ToString /@ #] & /@ Tuples[{0, 1}, 3]), Frame -> True, FrameLabel -> {Row[Subscript["m", #] & /@ Range[3]], "shots"},
     PlotLabel -> Row[{label, ":  ", Length[s], " formal bits,  ", Subscript["m", 3], " = ", TraditionalForm[BooleanFunction[
          Thread[Tuples[{0, 1}, 2] -> (m3 /. Thread[s -> #] & /@ Tuples[{0, 1}, 2])], {Subscript["s", 1], Subscript["s", 2]}, "ANF"]]}], ImageSize -> 260]];
BlockRandom[SeedRandom[1]; Grid[{{parityPanel["no error", parityState], parityPanel["X on qubit 1", parityState[{"X" -> 1}]]},
  {parityPanel["X on qubit 3", parityState[{"X" -> 3}]], parityPanel["Z on qubit 1", parityState[{"Z" -> 1}]]}}]]
```

Each state allocates two formal bits, and the shots split evenly over four strings: two fair coins and one dependent bit. Without an error $m_3=s_1\veebar s_2$ and every shot has even parity. An $X$ on qubit 1 or on qubit 3 negates the stored function and moves every shot to odd parity, which is how a parity check detects a bit flip; a $Z$ error changes nothing that a $Z$ measurement sees.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### T gates as a sum of stabilizer states: how many terms survive compression


The gate $T=\mathrm{diag}(1,e^{i\pi/4})$ is not a Clifford gate, so a state made with $t$ T gates is stored as a sum of stabilizer states. On $n$ qubits, two circuits use $t$ T gates: in the first, qubits $1,\dots,t$ each get $H$ then $T$, and a CNOT ladder then entangles the chain; in the second, the ladder first makes a GHZ state and the T gates follow. How many terms does the frame keep after compression as $n$ grows?

As an operator, $T$ on qubit $q$ is a sum of two Pauli operators,

$$
T_q=\frac{1+e^{i\pi/4}}{2}\,\mathbb 1+\frac{1-e^{i\pi/4}}{2}\,Z_q,\qquad |\psi\rangle=\sum_{i=1}^{2^t}c_i\,|s_i\rangle ,
$$

so each T gate doubles the list of stabilizer states $|s_i\rangle$ in the frame. The overlaps $\langle s_i|s_j\rangle$ follow from the tableaux, and from them the reduction keeps a largest linearly independent subset of the $|s_i\rangle$ and recomputes the $c_i$ exactly, so $|\psi\rangle$ is unchanged.

The figure counts each state's terms after this reduction against the number of T gates, $t=1,\dots,5$, on a log scale, computed at $n=5,10,\dots,30$:

```wl
magicLadder[n_] := Table["CNOT" -> {i, i + 1}, {i, n - 1}];
magicProduct[n_, t_] := Join[Thread["H" -> Range[t]], Thread["T" -> Range[t]], magicLadder[n]];
magicGHZ[n_, t_] := Join[{"H" -> 1}, magicLadder[n], Thread["T" -> Range[t]]];
magicTerms[n_, gates_] := PauliStabilizer[n][gates]["Compress"]["Length"];
magicCounts = Table[magicTerms[n, #[n, t]], {t, 5}, {n, 5, 30, 5}] & /@ {magicProduct, magicGHZ};
ListLogPlot[Map[Transpose[{Range[5], First /@ #}] &, magicCounts], Joined -> True, PlotMarkers -> {Automatic, 12},
  PlotLegends -> {"H, T on qubits 1..t, then ladder", "ladder (GHZ), then T on qubits 1..t"}, Frame -> True, AspectRatio -> 0.7,
  FrameLabel -> {"T gates t", "stabilizer terms after compression"}, FrameTicks -> {{2^Range[0, 5], None}, {Range[5], None}}, PlotRange -> {All, {0.8, 45}},
  PlotLabel -> If[AllTrue[Flatten[magicCounts, 1], Equal @@ # &], "the same counts at every n = 5, 10, \[Ellipsis], 30", "the counts depend on n"], ImageSize -> 460]
```

When each T gate acts on its own $|+\rangle$, the $2^t$ terms are linearly independent and compression removes none, at every $n$, while the state vector has $2^n$ amplitudes. On the GHZ state the T gates only change the relative phase of $|0\cdots0\rangle$ and $|1\cdots1\rangle$, and two terms suffice for every $t$. For even $t$ that phase is a power of $i$, so the state is a single stabilizer state, yet the frame still keeps two terms: the count bounds the stabilizer rank from above.

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

$\varphi$ is held as a phase polynomial, with the register starting in $|+\rangle^{\otimes n}$; the figure times, fastest of three, building it, reading the amplitude of $|10\cdots0\rangle$, and reading the same amplitude from the dense state vector of the circuit with an $H$ layer in front:

```wl
phasePolyGates[n_] := Join[Thread["T" -> Range[n]], Table["CZ" -> {q, q + 1}, {q, n - 1}], {"C"["Z" -> 3, {1, 2}]}];
phasePolyState[n_] := QuantumCircuitOperator[phasePolyGates[n]][Method -> {"Stabilizer", "Compress" -> "PhasePolynomial"}];
phasePolyY[n_] := PadRight[{1}, n];
SetAttributes[phasePolyTime, HoldFirst];
phasePolyTime[expr_] := Min @ Table[ClearSystemCache[]; First @ AbsoluteTiming[expr], 3];
phasePolyTimes = {
   Table[{n, phasePolyTime[phasePolyState[n]]}, {n, 10, 100, 10}],
   Table[With[{s = phasePolyState[n]}, {n, phasePolyTime[s["Amplitude", phasePolyY[n]]]}], {n, 10, 100, 10}],
   Table[{n, phasePolyTime[QuantumCircuitOperator[Join[Thread["H" -> Range[n]], phasePolyGates[n]]][]["StateVector"][[2^(n - 1) + 1]]]}, {n, 4, 18, 2}]};
ListLogPlot[phasePolyTimes, Joined -> True, PlotMarkers -> Automatic, Frame -> True, FrameLabel -> {"qubits n", "seconds"},
  PlotLegends -> {"build the phase polynomial", "one amplitude from it", "dense state vector"},
  PlotLabel -> Style[Row[{"\[LeftAngleBracket]10\[CenterEllipsis]0|\[Psi]\[RightAngleBracket] at n = 20:  ", TraditionalForm[phasePolyState[20]["Amplitude", phasePolyY[20]]]}], 13, FontFamily -> "Helvetica"],
  ImageSize -> 480]
```

The dense state vector rises from a fixed floor into exponential growth, a straight line on the log scale. Building $\varphi$ and reading one amplitude from it grow slowly with $n$ and continue to $n=100$, where a dense vector of $2^{100}$ entries cannot be stored, because an amplitude is $\varphi$ evaluated at one bit string rather than a sum over $2^n$ terms. At $n=20$ that amplitude is a single term: the $2^{-n/2}$ of a flat superposition times the phase $\omega$ that the one $T$ acting on $x_1=1$ wrote into it.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Repetition code: the logical error rate is a polynomial in the flip probability


One bit is stored in an odd number $n$ of qubits as $|0\cdots0\rangle$. Each qubit is flipped independently with probability $p$, and the bit is decoded by majority vote. For which $p$ does a longer code lower the probability $P_L$ that the decoded bit is wrong?

Each qubit passes through the bit-flip channel

$$
\mathcal E_p(\rho)=(1-p)\,\rho+p\,X\rho X,
$$

and majority vote fails when more than $n/2$ qubits flipped, so

$$
P_L(n)=\sum_{k>n/2}\binom nk p^k(1-p)^{n-k}=I_p\left(\frac{n+1}2,\frac{n+1}2\right).
$$

The bit-flip channel acts on each qubit of the stabilizer state $|0\cdots0\rangle$; every application splits each branch into an unflipped and a flipped tableau with weights $1-p$ and $p$, and qubit $k$ reads 1 when its generator $Z_k$ has sign $-1$. Summing the weights of the failing branches gives $P_L(p)$ exactly; the regularized beta function $I_p$ (dashed) extends it to long codes:

```wl
repFlip[branches_, q_] := Join @@ (Function[{w, ps}, {w #1, #2} & @@@ QuantumChannel["BitFlip"[p], {q}][ps]] @@@ branches);
repRate[n_] := Expand[Total[Cases[Fold[repFlip, {{1, PauliStabilizer[n]}}, Range[n]], {w_, ps_} /; Count[ps["StabilizerSigns"], -1] > n/2 :> w]]];
repRates = repRate /@ {1, 3, 5, 7};
repPanels = {Plot[Evaluate[Join[repRates, BetaRegularized[p, (# + 1)/2, (# + 1)/2] & /@ {21, 101}]], {p, 0, 1},
    PlotStyle -> {Automatic, Automatic, Automatic, Automatic, Dashed, Dashed}, PlotLegends -> (Row[{"n = ", #}] & /@ {1, 3, 5, 7, 21, 101}), Frame -> True, FrameLabel -> {"p", Subscript["P", "L"]},
    PlotLabel -> Row[{Subscript["P", "L"], "(n+2) < ", Subscript["P", "L"], "(n)  \[DoubleLongLeftRightArrow]  ",
       First[Union[MapThread[Reduce[#2 < #1 && 0 < p < 1, p] &, {Most[repRates], Rest[repRates]}]]]}], ImageSize -> 320],
   LogLogPlot[Evaluate[repRates], {p, 10^-3, 1/2}, Frame -> True, FrameLabel -> {"p", Subscript["P", "L"]},
    PlotLabel -> Row[{Subscript["P", "L"], "(3) = ", repRates[[2]]}], ImageSize -> 320]};
Legended[Row[{First[First[repPanels]], Last[repPanels]}, "  "], legendBelow[Last[First[repPanels]]]]
```

For $n=3$ the rate has no linear term, so at small $p$ the code fails far less often than a bare qubit; on the log-log axes each code becomes a straight line at small $p$, with slope equal to the smallest number of flips that defeats the majority. All curves meet at $p=1/2$ and tend to a step there as $n$ grows, and the exact threshold is the same for each pair of neighboring codes drawn: below one half the longer code is better, above it worse. For every odd $n=2m+1$ the difference factors as $P_L(n+2)-P_L(n)=\binom{2m+1}{m}\,{}[p(1-p)]^{m+1}(2p-1)$, which has the sign of $2p-1$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A coherent error on the three-qubit code, undone by measurement records


Encode a qubit in three, rotate the middle one by the angle $2\theta$ about $x$ (a coherent error, not a Pauli flip), measure the two parities, and flip back the qubit the parities point to; for comparison, rotate all three qubits instead. Each measurement outcome is kept as a record qubit, so each correction is a gate controlled by the records. Does the cycle restore the encoded state for every $\theta$?

With the input $|\psi\rangle=\cos\frac\beta2|0\rangle+e^{i\phi}\sin\frac\beta2|1\rangle$, the code word, the two errors and the measured parity bits are

$$
|\bar\psi\rangle=\cos\frac\beta2|000\rangle+e^{i\phi}\sin\frac\beta2|111\rangle,\qquad e^{-i\theta X_2}\ \text{ or }\ \prod_{k=1}^3e^{-i\theta X_k},\qquad m_{12}=\frac{1-Z_1Z_2}{2},\quad m_{23}=\frac{1-Z_2Z_3}{2},
$$

each measured by its two eigenspace projectors, which leave superpositions inside an eigenspace intact. The records $(m_{12},m_{23})=(1,0),(1,1),(0,1)$ trigger $X_1$, $X_2$, $X_3$.

For general $\theta$, $\beta$ and $\phi$, each panel shows the record probabilities and the fidelity $\sqrt{\langle\bar\psi|\rho|\bar\psi\rangle}$ of the three code qubits' state $\rho$ with $|\bar\psi\rangle$ for the inputs $|0\rangle$ and $|+\rangle$:

```wl
parityBit = (1 - QuantumOperator["ZZ"])/2;
qecIn = QuantumTensorProduct[QuantumState[{Cos[β/2], Exp[I ϕ] Sin[β/2]}], QuantumState["00"]];
qecEncoder = {"CNOT" -> {1, 2}, "CNOT" -> {1, 3}};
qecCode = QuantumCircuitOperator[qecEncoder][qecIn];
qecCircuit[err_] := QuantumCircuitOperator[Join[qecEncoder, {"RX"[2 θ] -> err,
    "Lueders"[parityBit] -> {1, 2}, "Lueders"[parityBit] -> {2, 3},
    "C"["X" -> 1, {0}, {-1}], "C"["X" -> 2, {0, -1}, {}], "C"["X" -> 3, {-1}, {0}]}]];
qecPanel[label_, err_] := With[{out = qecCircuit[err][qecIn]},
   {rec = FullSimplify[Diagonal[Normal[QuantumPartialTrace[out, {3, 4, 5}]["DensityMatrix"]]], {θ, β, ϕ} ∈ Reals],
     fid = FullSimplify[QuantumSimilarity[QuantumPartialTrace[out, {1, 2}], qecCode, "Fidelity"], {θ, β, ϕ} ∈ Reals]},
    Plot[Evaluate[{rec[[1]], rec[[2]] + rec[[3]], rec[[4]], fid /. {β -> 0, ϕ -> 0}, fid /. {β -> π/2, ϕ -> 0}}], {θ, 0, π},
     PlotRange -> {-0.05, 1.05}, PlotStyle -> {Automatic, Automatic, Automatic, Dashed, DotDashed}, Frame -> True,
     FrameLabel -> {"θ", "probability, fidelity"}, FrameTicks -> {{Automatic, None}, {{#, Style[TraditionalForm[#], FractionBoxOptions -> {Beveled -> True}]} & /@ Range[0, π, π/4], None}},
     PlotLegends -> {"records (0,0)", "one odd parity", "records (1,1)", "fidelity, input |0⟩", "fidelity, input |+⟩"},
     PlotLabel -> Column[{label, If[fid === 1, "fidelity 1 for all β, ϕ", "fidelity depends on β, ϕ"]}, Alignment -> Center],
     ImageSize -> 340]];
qecPanels = {qecPanel["RX(2θ) on qubit 2", {2}], qecPanel["RX(2θ) on every qubit", {1, 2, 3}]};
Column[{Labeled[qecCircuit[{2}]["Diagram", ImageSize -> 560], Style["encode, error on qubit 2, two parity measurements, corrections controlled by the records", 13, FontFamily -> "Helvetica"], Top],
  Legended[PlotGrid[{First /@ qecPanels}, ImageSize -> 680], legendBelow[Last[First[qecPanels]]]]}, Center, Spacings -> 2]
```

With the error on qubit 2 a single odd parity never occurs, so the measurement turns the rotation into either no error or a full flip of the middle qubit, and the record-controlled gates undo both branches for every $\theta$ and every input. When all three qubits rotate, two flips are mistaken for one and three go unseen, so each branch leaves a superposition of the code word and its logical flip $\bar X=XXX$: a logical rotation about $\bar x$. It lowers the fidelity of $|0\rangle$ but leaves $|+\rangle$, an eigenstate of $\bar X$, untouched, so the code fails and the damage depends on the input.

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
P_L=\sum_{E}\Pr(E)\,\left[\operatorname{class}(E)\neq\operatorname{guess}(\operatorname{syndrome}(E))\right].
$$

Summing over all $4^n$ errors for general $p$, with each syndrome's most probable class chosen at small $p$, gives $P_L$ as a polynomial for each decoder; the panels plot $P_L/p^2$, whose value at $p=0$ is the coefficient of $p^2$:

```wl
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/QEC/QECCore/QECCore.wl"]; Needs["Wolfram`QuantumFramework`QEC`"];
mlDecCases = Append[MapThread[{#1, #2 <> ", depolarizing", QECNoiseModel["Depolarizing", p]} &, {{"5QubitCode", "ShorCode", "SteaneCode"}, {"five-qubit code", "Shor code", "Steane code"}}],
   {"SteaneCode", "Steane code, Z-biased", QECNoiseModel[<|"X" -> p/10, "Y" -> p/10, "Z" -> p|>]}];
mlDecRates[code_, noise_] := Expand[QECLogicalErrorRate[QECCode[code], noise, "Decoder" -> #]] & /@ {"MaximumLikelihood", "MinimumWeight"};
mlDecPanel[{code_, title_, noise_}] := With[{r = mlDecRates[code, noise]},
   Plot[Evaluate[Expand[r/p^2]], {p, 0, 0.2}, PlotStyle -> {Automatic, Dashed}, PlotRange -> All, PlotLegends -> {"maximum likelihood", "lookup table"},
    Frame -> True, FrameLabel -> {"p", Row[{Subscript["P", "L"], " / ", Superscript["p", 2]}]}, ImageSize -> 360,
    PlotLabel -> Style[Column[{title, If[SameQ @@ r, Row[{Subscript["P", "L"], " = ", TraditionalForm[First[r]]}],
         Row[{"order ", Superscript["p", 2], ":  ", TraditionalForm[Coefficient[r[[1]], p, 2] p^2], "  against  ", TraditionalForm[Coefficient[r[[2]], p, 2] p^2]}]]}, Alignment -> Center], 12, FontFamily -> "Helvetica"]]];
mlDecPanels = mlDecPanel /@ mlDecCases;
Legended[Grid[Partition[First /@ mlDecPanels, 2], Spacings -> {2, 2}], legendBelow[Last[First[mlDecPanels]]]]
```

The five-qubit code gives one polynomial for both decoders: its fifteen nonzero syndromes come from the fifteen one-qubit errors, so the lightest error already lies in the most probable class. The Shor code is degenerate, with many light errors that differ by a check and so share a class, and adding up their probabilities is where maximum likelihood lowers the $p^2$ coefficient most. The Steane code gains less, under depolarizing and $Z$-biased noise alike, and in every panel the order of the two curves holds over the whole range of $p$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Noise inside the circuit that measures the checks


The five-qubit code corrects any one error when its checks are read perfectly; the bit-flip and phase-flip codes correct any one error of a single type. Here a circuit reads the checks, and every gate, reset and readout in it fails with probability $p$. What does noise inside the measuring circuit do to the protection the code gives?

Each check $g$ is measured with a reset ancilla: every data qubit in the support of $g$ is rotated so that its letter becomes $Z$ ($H$ for $X$), CNOTs copy those qubits' parity into the ancilla, the rotations are undone, and the ancilla is read. After each one-qubit gate $X$, $Y$ or $Z$ strikes with probability $p/3$ each, after each CNOT one of the 15 two-qubit Paulis with probability $p/15$ each, and each reset and each readout is flipped with probability $p$; waiting qubits do not err. One noisy round is followed by a perfect readout of the data, and the decoder reads the round's syndrome and its change at the final readout. The transversal extraction instead measures each check through a verified cat state with one qubit per data qubit of $g$, each coupled to its data qubit by one two-qubit gate, and keeps only shots that pass the verification. As a reference, code capacity puts depolarizing noise of rate $p$ on the data and reads the checks perfectly.

$P_L$ for the five-qubit code (checks the cyclic shifts of $XZZXI$) is an exact polynomial in $p$, shown with the transversal rate at six values of $p$; the second panel plots $P_L/p$ for the bit-flip code (checks $Z_1Z_2$, $Z_2Z_3$) and the phase-flip code (checks $X_1X_2$, $X_2X_3$):

```wl
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/QEC/QECCore/QECCore.wl"]; Needs["Wolfram`QuantumFramework`QEC`"];
faultCircNoise = <|"code capacity" -> QECNoiseModel["Depolarizing", p], "circuit, all locations" -> QECNoiseModel["Circuit", p],
   "circuit, one-qubit gates noiseless" -> QECNoiseModel["Circuit", <|"TwoQubit" -> p, "Measurement" -> p, "Reset" -> p|>]|>;
faultCircRate[code_, "code capacity"] := Expand[QECLogicalErrorRate[QECCode[code], faultCircNoise["code capacity"]]];
faultCircRate[code_, key_] := Expand[QECLogicalErrorRate[QECCode[code], faultCircNoise[key], "Rounds" -> 1]];
faultCircFive = faultCircRate["5QubitCode", #] & /@ {"code capacity", "circuit, all locations"};
faultCircTransversal = Table[{q, QECLogicalErrorRate[QECCode["5QubitCode"], QECNoiseModel["Circuit", q], "Rounds" -> 1, "Extraction" -> "Transversal"]}, {q, 10.^Range[-4, -1.5, 0.5]}];
faultCircPairs = Append[Tuples[{{"code capacity", "circuit, all locations"}, {"BitFlipCode", "PhaseFlipCode"}}], {"circuit, one-qubit gates noiseless", "PhaseFlipCode"}];
faultCircDual = faultCircRate[#2, #1] & @@@ faultCircPairs;
faultCircLabel[s_] := Style[s, 12, FontFamily -> "Helvetica"];
Grid[{{
   ListLogLogPlot[Append[Table[{q, # /. p -> q}, {q, 10.^Range[-4, -1.5, 0.05]}] & /@ faultCircFive, faultCircTransversal], Joined -> {True, True, False}, PlotMarkers -> {None, None, Automatic},
    PlotLegends -> Placed[{"code capacity", "bare-ancilla circuit", "transversal circuit, kept shots"}, Below], Frame -> True, FrameLabel -> {"p", Subscript["P", "L"]}, ImageSize -> 380,
    PlotLabel -> faultCircLabel[Column[{"five-qubit code", Row[{"bare ancilla:  ", Subscript["P", "L"], " = ", TraditionalForm[Coefficient[Last[faultCircFive], p, 1] p], " + O(", Superscript["p", 2], ")"}]}, Alignment -> Center]]],
   Plot[Evaluate[Expand[Cancel[faultCircDual/p]]], {p, 0, 0.05}, PlotStyle -> {Automatic, Dashed, Automatic, Dashed, DotDashed}, PlotRange -> All, Frame -> True,
    FrameLabel -> {"p", Row[{Subscript["P", "L"], " / p"}]}, ImageSize -> 380,
    PlotLegends -> Placed[LineLegend[Automatic, Row[{#2 /. {"BitFlipCode" -> "bit flip", "PhaseFlipCode" -> "phase flip"}, ", ", #1}] & @@@ faultCircPairs, LegendLayout -> "Column"], Below],
    PlotLabel -> faultCircLabel[Column[{"bit-flip and phase-flip codes", Row[{"all locations, bit flip ", TraditionalForm[Coefficient[faultCircDual[[3]], p, 1] p], ",  phase flip ", TraditionalForm[Coefficient[faultCircDual[[4]], p, 1] p]}]}, Alignment -> Center]]]
   }}, Alignment -> Top, Spacings -> 2]
```

With the bare ancilla the five-qubit code fails at first order in $p$: one fault on an ancilla partway through its CNOTs spreads to several data qubits. The transversal extraction restores the second-order slope of code capacity on the kept shots, although at the largest values of $p$ it fails more often than the bare ancilla. The bit-flip and phase-flip codes fail identically at code capacity but not once their checks are measured by a circuit, because only the $X$ checks need $H$ gates on the data: with one-qubit-gate noise switched off, the phase-flip code falls back onto the bit-flip curve.

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

## Photonics: interferometers with symbolic photon numbers

`QuantumInterferometer` computes linear-optics amplitudes for photon numbers that stay symbols, with loss, squeezers and coherent inputs, and converts an interferometer to a Fock-space circuit when a nonlinear element is needed.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A beam splitter with symbolic photon numbers


$n_1$ and $n_2$ photons enter the two ports of a beam splitter with mixing angle $\theta$ and phase $\varphi$. What is the probability of finding $t_1$ photons in the first output, for every photon number at once?

The splitter rotates the creation operators,

$$
\hat a_1^\dagger\to\cos\theta\,\hat a_1^\dagger+e^{i\varphi}\sin\theta\,\hat a_2^\dagger,\qquad \hat a_2^\dagger\to-e^{-i\varphi}\sin\theta\,\hat a_1^\dagger+\cos\theta\,\hat a_2^\dagger,
$$

so it conserves the photon number, and the second output holds $n_1+n_2-t_1$ photons.

The probability for general $n_1,n_2,t_1,\theta,\varphi$ comes out in closed form (shown at $\varphi=0$); each panel fixes the input photon numbers:

```wl
splitterProb = QuantumInterferometer[{"BS"[θ, φ]}]["Probability", {n1, n2} -> {t1, n1 + n2 - t1}];
twinPanel[{p_, q_}] := Plot[Evaluate[FullSimplify[Table[splitterProb /. {n1 -> p, n2 -> q, t1 -> j}, {j, 0, p + q}], Element[{θ, φ}, Reals]]],
    {θ, 0, π/2}, PlotRange -> {0, 1}, PlotStyle -> Table[Directive[If[OddQ[j], Dashed, {}], AbsoluteThickness[If[j < 2, 4, 1.6]]], {j, 0, p + q}],
    PlotLegends -> Table[Row[{Subscript["t", 1], " = ", j}], {j, 0, p + q}], GridLines -> {{π/4}, None}, Frame -> True,
    FrameTicks -> {{Automatic, None}, {{#, Style[TraditionalForm[#], FractionBoxOptions -> {Beveled -> True}]} & /@ Range[0, π/2, π/8], None}}, FrameLabel -> {"θ", Row[{"P(", Subscript["t", 1], ")"}]},
    PlotLabel -> Row[{Subscript["n", 1], " = ", p, ",   ", Subscript["n", 2], " = ", q}], ImageSize -> 300];
twinPanels = twinPanel /@ {{1, 1}, {2, 2}, {1, 3}};
Labeled[Legended[PlotGrid[{First /@ twinPanels}, ImageSize -> 900], legendBelow[Last[Last[twinPanels]]]],
  Style[TraditionalForm[splitterProb /. {φ -> 0, n1 -> Subscript[n, 1], n2 -> Subscript[n, 2], t1 -> Subscript[t, 1]}], 14, FontFamily -> "Times"], Top]
```

The probability is the squared Wigner D-function of a spin $j=(n_1+n_2)/2$ rotated by $2\theta$: two modes holding $N$ photons carry the spin-$N/2$ representation of SU(2), and the splitter rotates it. The phase $\varphi$ enters the D-function only as a phase factor, so it drops out of the probability. At the balanced splitter (vertical line) equal inputs never leave an odd $t_1$ (dashed), with $n_1=n_2=1$ the Hong-Ou-Mandel zero, while one and three photons never leave two in each output.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Photon loss, and a Hong-Ou-Mandel zero that survives it


One photon enters each input of a beam splitter with angle $\theta$ and phase $\varphi$, and attenuators of transmissivities $\eta_1$ and $\eta_2$ sit on the two outputs. Does the two-photon coincidence zero survive the loss?

Pure loss mixes the mode with a vacuum mode $\hat b$ and discards $\hat b$: $\hat a\to\sqrt\eta\,\hat a+\sqrt{1-\eta}\,\hat b$, so each photon survives independently with probability $\eta$. The splitter $\hat a_1^\dagger\to\cos\theta\,\hat a_1^\dagger+e^{i\varphi}\sin\theta\,\hat a_2^\dagger$, $\hat a_2^\dagger\to-e^{-i\varphi}\sin\theta\,\hat a_1^\dagger+\cos\theta\,\hat a_2^\dagger$ acts first, then the loss on each output, and the coincidence probability is

$$
P(1,1)=\eta_1\eta_2\cos^2 2\theta \quad\text{for every }\varphi .
$$

The coincidence probability for general $\theta,\varphi,\eta_1,\eta_2$ comes out in closed form:

```wl
lossyCoincidence = FullSimplify[QuantumInterferometer[{"BS"[θ, φ], "Loss"[η1] -> 1, "Loss"[η2] -> 2}]["Probability", {1, 1} -> {1, 1}]];
lossPairs = {{1, 1}, {0.9, 0.5}, {0.6, 0.3}};
Plot[Evaluate[lossyCoincidence /. Thread[{η1, η2} -> #] & /@ lossPairs], {θ, 0, π/2}, PlotRange -> {0, 1},
    PlotLegends -> (Row[{Subscript["η", 1], " = ", #1, ",  ", Subscript["η", 2], " = ", #2}] & @@@ lossPairs), GridLines -> {{π/4}, None}, Frame -> True,
    FrameTicks -> {{Automatic, None}, {Range[0, π/2, π/8], None}}, FrameLabel -> {"θ", "P(1, 1)"}, PlotLabel -> Style[TraditionalForm[lossyCoincidence /. {η1 -> Subscript[η, 1], η2 -> Subscript[η, 2]}], 14, FontFamily -> "Times"], ImageSize -> 420]
```

Loss only scales the coincidence curve by $\eta_1\eta_2$, the chance that both photons survive, so its zero at the balanced splitter stays for every pair of transmissivities, equal or not: a coincidence needs both photons to survive, and a surviving pair carries the lossless amplitude.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Hong-Ou-Mandel with partially distinguishable photons


Two photons meet at a beam splitter with angle $\theta$ and phase $\varphi$. The first is horizontally polarized, the second is polarized at angle $\chi$, $\cos\chi\,|H\rangle+\sin\chi\,|V\rangle$, and the detectors do not resolve polarization. How does the coincidence probability depend on $\theta$ and $\chi$?

Treat the polarizations as extra modes $1_H,2_H,1_V,2_V$, with the same splitter $\hat a_1^\dagger\to\cos\theta\,\hat a_1^\dagger+e^{i\varphi}\sin\theta\,\hat a_2^\dagger$, $\hat a_2^\dagger\to-e^{-i\varphi}\sin\theta\,\hat a_1^\dagger+\cos\theta\,\hat a_2^\dagger$ on the $H$ pair and on the $V$ pair. The input is $\cos\chi\,|1,1,0,0\rangle+\sin\chi\,|1,0,0,1\rangle$, and the coincidence sums over the four outputs $T$ with one photon in each path:

$$
P_{\rm coinc}=\sum_T\left|\cos\chi\,A(1100\to T)+\sin\chi\,A(1001\to T)\right|^2 .
$$

Summing the polarization-resolved probabilities gives $P_{\rm coinc}(\theta,\chi)$ for general $\theta,\varphi,\chi$ in closed form:

```wl
polarizedHOM = QuantumInterferometer[{"BS"[θ, φ] -> {1, 2}, "BS"[θ, φ] -> {3, 4}}, 4];
partialCoincidence = Collect[ComplexExpand @ FullSimplify[Total[Abs[Cos[χ] polarizedHOM["Amplitude", {1, 1, 0, 0} -> #] + Sin[χ] polarizedHOM["Amplitude", {1, 0, 0, 1} -> #]]^2 & /@
    {{1, 1, 0, 0}, {1, 0, 0, 1}, {0, 1, 1, 0}, {0, 0, 1, 1}}], Element[{θ, φ, χ}, Reals]], Sin[χ]^2, Simplify];
Plot[Evaluate[Table[partialCoincidence /. χ -> c, {c, {0, π/8, π/4, 3 π/8, π/2}}]], {θ, 0, π/2}, PlotRange -> {0, 1},
  PlotLegends -> {"χ = 0", "χ = π/8", "χ = π/4", "χ = 3π/8", "χ = π/2"}, GridLines -> {{π/4}, None},
  Frame -> True, FrameTicks -> {{Automatic, None}, {Range[0, π/2, π/8], None}}, FrameLabel -> {"θ", "P(coincidence)"}, PlotLabel -> TraditionalForm[partialCoincidence], ImageSize -> 420]
```

The phase $\varphi$ drops out, and $P_{\rm coinc}$ is the Hong-Ou-Mandel curve $\cos^22\theta$ of identical photons weighted by the squared polarization overlap $\cos^2\chi$, plus the distinguishable-particle curve $(3+\cos4\theta)/4=\cos^4\theta+\sin^4\theta$ weighted by $\sin^2\chi$. At the balanced splitter only the second term survives, so the depth of the dip measures how indistinguishable the photons are.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Two squeezed vacua on a balanced beam splitter


Each of two modes is squeezed from the vacuum, with complex parameters $\xi_1$ and $\xi_2$, and the modes meet on a balanced beam splitter. When does one photon appear in each output?

Each squeezer is $\hat S(\xi)=\exp\!\left[(\xi^*\hat a^2-\xi\,\hat a^{\dagger2})/2\right]$; the splitter then maps $\hat a_1^\dagger\to(\hat a_1^\dagger+\hat a_2^\dagger)/\sqrt2$, $\hat a_2^\dagger\to(\hat a_2^\dagger-\hat a_1^\dagger)/\sqrt2$.

The coincidence probability $P(1,1)$ for general $\xi_1,\xi_2$ is plotted for $\xi_1=r$, $\xi_2=q\,re^{i\phi}$ with $q=|\xi_2|/|\xi_1|$ together with the pair probabilities $P(m,m)$ for opposite squeezers $\pm r$ for general $r$, all in closed form:

```wl
pairProb = QuantumInterferometer[{"S"[ξ1] -> 1, "S"[ξ2] -> 2, "BS"[π/4]}]["Probability", {0, 0} -> {1, 1}];
oppositeSqueezers = QuantumInterferometer[{"S"[r] -> 1, "S"[-r] -> 2, "BS"[π/4]}];
squeezerPairLaw = FullSimplify[Table[oppositeSqueezers["Probability", {0, 0} -> {m, m}], {m, 0, 3}], r > 0];
coincidencePanel[q_] := Plot[Evaluate[Table[pairProb /. {ξ1 -> r, ξ2 -> q r Exp[I p]}, {p, {0, π/2, π}}]], {r, 0, 2},
    PlotRange -> {-0.005, 0.26}, PlotLegends -> {"ϕ = 0", "ϕ = π/2", "ϕ = π"}, Frame -> True, FrameLabel -> {"r", "P(1, 1)"},
    PlotLabel -> Row[{"|", Subscript["ξ", 2], "|/|", Subscript["ξ", 1], "| = ", ToString[q, InputForm]}], ImageSize -> 300];
coincidencePanels = coincidencePanel /@ {1/2, 1};
Grid[{{Labeled[Legended[Grid[{First /@ coincidencePanels}], legendBelow[Last[Last[coincidencePanels]]]],
    Style[TraditionalForm[pairProb /. {ξ1 -> Subscript[ξ, 1], ξ2 -> Subscript[ξ, 2]}], 13, FontFamily -> "Times"], Top],
  Plot[Evaluate[squeezerPairLaw], {r, 0, 2}, PlotRange -> {0, 1}, PlotLegends -> Table[Row[{"m = ", m - 1, ":  ", TraditionalForm[squeezerPairLaw[[m]]]}], {m, 4}],
    Frame -> True, FrameLabel -> {"r", "P(m, m)"}, PlotLabel -> Row[{Subscript["ξ", 1], " = r,  ", Subscript["ξ", 2], " = -r"}], ImageSize -> 300]}}]
```

The coincidence vanishes only for identical squeezers ($q=1$, $\phi=0$): the pair generator $\hat a_1^{\dagger2}+\hat a_2^{\dagger2}$ is unchanged by the splitter, so the output is again two independent squeezed vacua, each holding only even photon numbers. Any difference in strength or phase lets coincidences through, most for opposite squeezers, whose outputs carry $m$ photons each with the two-mode squeezed vacuum statistics $P(m,m)=\tanh^{2m}r/\cosh^2r$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Fourier interferometers with one photon per port


One photon enters each port of an $N$-mode interferometer whose transfer matrix is the discrete Fourier matrix. Which output patterns occur?

The network maps $\hat a_k^\dagger\to\sum_j W_{jk}\,\hat a_j^\dagger$ with $W_{jk}=\omega^{(j-1)(k-1)}/\sqrt N$ and $\omega=e^{2\pi i/N}$; each output amplitude is a permanent of an $N\times N$ matrix built from $W$.

The figure gives the exact probability of every output pattern $\{t_1,\dots,t_N\}$ from $|1,\dots,1\rangle$, for $N=3$ and $N=4$:

```wl
fourierProbabilities[d_] := QuantumInterferometer[FourierMatrix[d]]["Probabilities", ConstantArray[1, d]];
fourierChart[d_] := With[{pr = fourierProbabilities[d]}, BarChart[Values[pr], ChartLabels -> Placed[StringJoin[ToString /@ #["Name"]] & /@ Keys[pr], Axis, Rotate[#, Pi/2] &],
    PlotRange -> {0, 0.4}, Frame -> True, FrameLabel -> {None, "probability"}, PlotLabel -> Row[{"output patterns,  N = ", d}], AspectRatio -> Full, ImageSize -> {80 + 18 Length[pr], 280}]];
Row[fourierChart /@ {3, 4}, Spacer[20]]
```

Most patterns are exactly zero, among them every two-plus-one pattern at $N=3$ and the input pattern $\{1,1,1,1\}$ itself at $N=4$. Relabeling the input ports cyclically leaves the input unchanged but multiplies the amplitude of $\{t_1,\dots,t_N\}$ by $\omega^{\sum_j (j-1)t_j}$, so a pattern can occur only if that sum is a multiple of $N$, and in both networks every pattern that passes this test does occur.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Six photons that the suppression law does not darken


One photon enters each port of a six-mode lossless network whose transfer matrix is the discrete Fourier matrix. The suppression law of Tichy and coworkers says that an output pattern $\{t_1,\dots,t_6\}$ is dark whenever $\sum_j (j-1)\,t_j$ is not a multiple of 6. Is every pattern the law allows actually bright?

The network maps $\hat a_k^\dagger\to\sum_j W_{jk}\,\hat a_j^\dagger$ with $W_{jk}=\omega^{(j-1)(k-1)}/\sqrt6$ and $\omega=e^{2\pi i/6}$. With one photon per input port, a pattern $T$ has probability

$$
P(T)=\frac{\left|\operatorname{Per}W_{T}\right|^2}{\prod_j t_j!},
$$

where $W_T$ uses every column of $W$ once and repeats row $j$ $t_j$ times.

The probabilities of all $\binom{11}{6}$ output patterns from $|1,1,1,1,1,1\rangle$ are exact rational numbers, so a zero is exact; the patterns are sorted by $\sum_j (j-1)\,t_j \bmod 6$:

```wl
fourierProbs[m_] := KeyMap[#["Name"] &, QuantumInterferometer[FourierMatrix[m]]["Probabilities", ConstantArray[1, m]]];
sixPhotonProbs = fourierProbs[6];
sixPhotonResidue = Mod[Range[0, 5] . #, 6] &;
sixPhotonCounts = GroupBy[Keys[sixPhotonProbs], {sixPhotonResidue, sixPhotonProbs[#] === 0 &}, Length];
sixPhotonAllowedDark = Select[Keys[sixPhotonProbs], sixPhotonResidue[#] == 0 && sixPhotonProbs[#] === 0 &];
sixPhotonExtra[m_] := Count[Normal @ fourierProbs[m], (t_ -> p_) /; Divisible[Range[0, m - 1] . t, m] && Chop[N[p]] == 0];
Labeled[BarChart[Table[Lookup[sixPhotonCounts[r], {False, True}, 0], {r, 0, 5}], ChartLayout -> "Stacked",
    ChartLabels -> {Range[0, 5], None}, ChartLegends -> {"bright", "exactly dark"}, Frame -> True,
    FrameLabel -> {Row[{"\[Sum] (j \[Minus] 1) ", Subscript["t", "j"], "  mod 6"}], "output patterns"}, ImageSize -> 480],
  Style[Column[{Row[{"allowed by the law, yet exactly dark (", Subscript["t", "1"], "\[Ellipsis]", Subscript["t", "6"], "):"}], Grid[Partition[Row /@ sixPhotonAllowedDark, 6], Spacings -> {2, 0.5}],
      Row[{"allowed yet dark, m = 2, 3, 4, 5, 6 modes:  ", Row[Append[sixPhotonExtra /@ Range[2, 5], Length[sixPhotonAllowedDark]], ", "]}]}, Alignment -> Center], 13, FontFamily -> "Helvetica"], Top]
```

Every pattern whose sum is not a multiple of 6 is exactly dark, as the law requires. Most patterns with sum 0 are bright, but some are exactly dark as well, each with one or two ports holding two photons. The law is sufficient for darkness but not necessary, and six ports is the smallest Fourier network where it fails to be necessary.

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

The network gives these amplitudes for $m=0,\dots,5$ (bars); the first three give the kept signal for general $c_0,c_1,c_2$:

```wl
nsGate = QuantumInterferometer[{{1 - Sqrt[2], 2^(-1/4), Sqrt[3/Sqrt[2] - 2]}, {2^(-1/4), 1/2, (1 - Sqrt[2])/2},
    {Sqrt[3/Sqrt[2] - 2], (1 - Sqrt[2])/2, Sqrt[2] - 1/2}}];
nsAmps = RootReduce[Table[nsGate["Amplitude", {m, 1, 0} -> {m, 1, 0}], {m, 0, 5}]];
keptSignal = {c0, c1, c2} Take[nsAmps, 3];
Labeled[BarChart[nsAmps, PlotLabels -> Range[0, 5],
    LabelingFunction -> Above, Frame -> True, FrameLabel -> {"signal photons m", "⟨m,1,0|U|m,1,0⟩"}, ImageSize -> 480],
  Style[Row[{"kept signal  ", TraditionalForm[keptSignal]}], 13, FontFamily -> "Helvetica"], Top]
```

For zero, one and two signal photons the amplitudes have the same size and only the two-photon one is negative, so every input comes out with $c_2$ flipped and the same overall factor: the gate works with a probability that does not depend on the input. From three signal photons on the amplitudes differ, so the gate works only on inputs with at most two signal photons.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Bosons, distinguishable particles and fermions at one beam splitter


Two particles, one in each input of a beam splitter with angle $\theta$ and phase $\varphi$, can both be transmitted or both be reflected. How likely is it that they leave through different ports?

With transfer matrix $W=\begin{pmatrix}\cos\theta & -e^{-i\varphi}\sin\theta\\ e^{i\varphi}\sin\theta & \cos\theta\end{pmatrix}$, bosons add the two histories with a plus sign, fermions with a minus sign, and distinguishable particles add their probabilities:

$$
P_{\rm B}=|\operatorname{Per}W|^2,\qquad P_{\rm F}=|\det W|^2,\qquad P_{\rm D}=\operatorname{Per}\left(|W_{jk}|^2\right).
$$

For general $\theta,\varphi$ the three coincidence probabilities follow from $W$ in closed form (legend):

```wl
statisticsW = QuantumInterferometer[{"BS"[θ, φ]}]["Unitary"];
statisticsCoincidences = FullSimplify[{Abs[Permanent[statisticsW]]^2, Permanent[Abs[statisticsW]^2], Abs[Det[statisticsW]]^2}, Element[{θ, φ}, Reals]];
Plot[Evaluate[statisticsCoincidences], {θ, 0, π/2}, PlotRange -> {0, 1.05}, GridLines -> {{π/4}, None},
  PlotLegends -> MapThread[Row[{#1, ":  ", TraditionalForm[#2]}] &, {{"bosons", "distinguishable", "fermions"}, statisticsCoincidences}],
  Frame -> True, FrameTicks -> {{Automatic, None}, {Range[0, π/2, π/8], None}}, FrameLabel -> {"θ", "P(1, 1)"}, ImageSize -> 420]
```

The phase $\varphi$ drops out of all three. Bosons never separate at the balanced splitter, fermions separate at every $\theta$, and distinguishable particles lie between them. Bosons and fermions differ only in the sign of the exchange term; distinguishable particles have none, so their curve is the average of the other two.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Photon loss and the Heisenberg advantage of a NOON state


A NOON state $(|n,0\rangle+|0,n\rangle)/\sqrt2$ enters the two arms of an interferometer; arm 2 carries an unknown phase $\theta$, and arm $j$ transmits each photon with probability $\eta_j$. Any estimate of $\theta$ has error at least $1/\sqrt{F}$, where $F$ is the quantum Fisher information, and $n$ independent photons reach $F=n$, the shot-noise limit. Does the NOON advantage keep growing with $n$ once photons are lost?

The phase acts as $e^{-i\theta\hat n_2}$, each arm loses photons through the beam-splitter loss channel of transmission $\eta_j$, and for the lossy state $\rho=\sum_i\lambda_i|i\rangle\langle i|$

$$
F=2\sum_{\lambda_i+\lambda_j>0}\frac{(\lambda_i-\lambda_j)^2}{\lambda_i+\lambda_j}\,\left|\langle i|\hat n_2|j\rangle\right|^2 .
$$

For symbolic $n$ the interferometer gives the probability that all $n$ photons survive each arm, $\eta_1^n$ and $\eta_2^n$. The phase lives only in the coherence between $|n,0\rangle$ and $|0,n\rangle$, and only the branch in which no photon is lost keeps it; a single Kraus operator acts on that branch, so it stays pure and the coherence is $\frac12\sqrt{\eta_1^n\eta_2^n}$, which gives $F=2n^2\eta_1^n\eta_2^n/(\eta_1^n+\eta_2^n)$ for every $n$. Curves are this $F/n$, markers are $F/n$ from the eigenpairs of the full lossy state, and the gridline at 1 is shot noise:

```wl
noonState[n_] := (FockState[{n, 0}, n + 1] + FockState[{0, n}, n + 1])/Sqrt[2];
noonLoss[η1_, η2_] := QuantumInterferometer[{"Loss"[η1] -> 1, "Loss"[η2] -> 2}, 2];
noonSurvival = noonLoss[η1, η2]["Probability", # -> #, Assumptions -> n ∈ PositiveIntegers] & /@ {{n, 0}, {0, n}};
noonFisher = 2 n^2 (Times @@ noonSurvival)/Total[noonSurvival];
noonQFI[ρ_, G_] := With[{λ = Re[First[#]], v = Last[#]},
    2 Total[Outer[If[#1 + #2 > 10^-10, (#1 - #2)^2/(#1 + #2), 0] &, λ, λ] Abs[Conjugate[v] . G . Transpose[v]]^2, 2]] & @ Eigensystem[ρ];
noonGain[η1v_, η2v_, nv_] := noonQFI[Normal @ noonLoss[η1v, η2v][noonState[nv]]["DensityMatrix"], KroneckerProduct[IdentityMatrix[nv + 1], DiagonalMatrix[Range[0, nv]]]]/nv;
noonPanel[η1v_] := With[{η2s = {1, 0.95, 0.9, 0.8}}, Show[
    ListPlot[Table[{nv, noonGain[η1v, η2v, nv]}, {η2v, η2s}, {nv, 20}], PlotLegends -> (Row[{Subscript["\[Eta]", 2], " = ", #}] & /@ η2s), PlotMarkers -> {Automatic, Small}],
    Plot[Evaluate[noonFisher/n /. {η1 -> η1v, η2 -> #} & /@ η2s], {n, 1, 20}],
    PlotRange -> {0, 11}, Frame -> True, FrameLabel -> {"n", "F/n"}, GridLines -> {None, {1}},
    PlotLabel -> Row[{Subscript["\[Eta]", 1], " = ", η1v}], ImageSize -> Small]];
noonPanels = noonPanel /@ {0.95, 0.9, 0.8};
Labeled[Legended[PlotGrid[{First /@ noonPanels}, ImageSize -> 660], legendBelow[Last[First[noonPanels]]]],
  With[{sub = {η1 -> Subscript["\[Eta]", 1], η2 -> Subscript["\[Eta]", 2]}},
   Style[Row[{"all ", n, " photons survive:  ", TraditionalForm[noonSurvival /. sub], ",     F = ", TraditionalForm[noonFisher /. sub]}], 13, FontFamily -> "Helvetica"]], Top]
```

The markers sit on the curves, which for equal arms reduce to $F=n^2\eta^n$. Since $(\eta_1\eta_2)^n$ eventually beats $n^2$, the gain peaks even if only one arm loses photons, and the more loss, the smaller the $n$ at which it peaks; with both arms at $\eta=0.8$ it falls back below shot noise within the plotted range.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What one boson-sampling amplitude costs


$n$ identical photons enter a random lossless network of $m=2n$ modes with occupation list $S=\{s_1,\dots,s_m\}$, and we ask for the amplitude of one output list $T$. How does the cost of that one number depend on $n$ and on how the photons share the modes?

With $W$ the $m\times m$ transfer matrix and $W_{T,S}$ the $n\times n$ matrix that repeats column $k$ $s_k$ times and row $j$ $t_j$ times,

$$
A(T\leftarrow S)=\frac{\operatorname{Per}W_{T,S}}{\sqrt{\prod_k s_k!\,\prod_j t_j!}} .
$$

Ryser's and Glynn's formulas evaluate a permanent in about $n\,2^n$ operations. When modes hold several photons, the amplitude can instead be grouped into a sum of $\min\{\prod_k(s_k+1),\prod_j(t_j+1)\}$ terms, each costing about one multiply per pair of occupied input and output modes; the computation here takes this grouped sum only when its estimated cost is about 50 times below $n\,2^n$.

Each amplitude is timed as the fastest of five runs, each with cleared caches, with the input in the first $n$ modes and the output pattern moved to the last $n$:

```wl
bsOccupations[n_] := PadRight[#, 2 n] & /@ <|"collision-free" -> ConstantArray[1, n],
    "pairs" -> Append[ConstantArray[2, Quotient[n, 2]], Mod[n, 2]], "bunched" -> {n}|>;
bsNetwork[n_] := BlockRandom[SeedRandom[1107 + 1000 n]; QuantumInterferometer["Random"[2 n]]];
bsSeconds[q_, s_] := With[{t = RotateRight[s, Total[s]]}, Min @ Table[ClearSystemCache[]; First @ AbsoluteTiming[q["Amplitude", s -> t]], 5]];
bsTimes = AssociationMap[Function[n, With[{q = bsNetwork[n]}, bsSeconds[q, #] & /@ bsOccupations[n]]], {12, 14, 16, 18, 20}];
ListLogPlot[Transpose[bsTimes], Joined -> True, PlotMarkers -> Automatic,
  Frame -> True, FrameLabel -> {"photons n  (m = 2n modes)", "seconds per amplitude"}, ImageSize -> 450]
```

The collision-free amplitude grows exponentially, a straight line on the log scale. For pairs at these photon numbers the grouped sum is not cheap enough (it overtakes the permanent only near $n=44$), so the permanent is computed and the paired curve tracks the collision-free one within timing noise, while the fully bunched amplitude needs only $n+1$ terms and stays essentially free.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Exact beam-splitter settings for an exact interferometer


A lossless interferometer on $m$ modes is fixed by its unitary transfer matrix $W$, here the discrete Fourier matrices $F_3$ and $F_4$, with $(F_m)_{jk}=e^{2\pi i(j-1)(k-1)/m}/\sqrt m$. Any such $W$ can be built from $m(m-1)/2$ beam splitters and $m$ phase shifters. When $W$ is known exactly, are the splitter settings exact numbers too?

The network maps $\hat a_k^\dagger\to\sum_j W_{jk}\,\hat a_j^\dagger$. A splitter on modes $p,q$ acts there by the block

$$
B(\theta,\varphi)=\begin{pmatrix}\cos\theta & -e^{-i\varphi}\sin\theta\\ e^{i\varphi}\sin\theta & \cos\theta\end{pmatrix},
$$

and a phase shifter multiplies one mode by $e^{i\alpha}$. Each splitter of the mesh is chosen to null one entry of $W$, so its $\theta$ is an arctangent of a ratio of moduli of entries and its $\varphi$ a difference of their arguments.

$F_3$ and $F_4$ are decomposed into a Clements mesh, the settings reduced, and the rebuilt mesh compared with $F_m$:

```wl
meshSettings[m_] := FullSimplify[FullSimplify[QuantumInterferometer[FourierMatrix[m]]["Elements"]] /. Arg[r_Root] :> Arg[ToRadicals[RootReduce[r/Abs[r]]]]] /.
   {"BS"[t_, f_] :> "BS"[t, Mod[f, 2 π]], "PS"[f_] :> "PS"[Mod[f, 2 π]]};
meshRow["BS"[t_, f_] -> q_] := {"splitter", Row[q, ", "], TraditionalForm[t], TraditionalForm[f]};
meshRow["PS"[f_] -> q_] := {"phase", Row[q, ", "], "", TraditionalForm[f]};
meshPanel[m_] := With[{el = meshSettings[m]}, With[{mesh = QuantumInterferometer[el, m]}, Labeled[
    Column[{mesh["Diagram", ImageSize -> 320],
      Style[Grid[Prepend[meshRow /@ el, {"element", "modes", "θ", "φ or α"}], Frame -> All, Spacings -> {1.5, 0.6}], 13, FontFamily -> "Helvetica"]}, Center, Spacer[10]],
    Style[Row[{"m = ", m, ":   rebuilt mesh equals ", Subscript["F", m], ":  ", AllTrue[Flatten @ FullSimplify[mesh["Unitary"] - FourierMatrix[m]], # === 0 &]}], 13, FontFamily -> "Helvetica"], Top]]];
Grid[{{meshPanel[3], meshPanel[4]}}, Alignment -> {Center, Top}, Spacings -> 4]
```

Every setting of both meshes comes out in closed form, the phases as rational multiples of $\pi$, and the mesh rebuilt from them equals $F_m$ exactly, so the decomposition is exact algebra, not a numerical fit.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A beam splitter in the ladder-operator algebra, with no photon cutoff


Two modes with annihilation operators $\hat a_1,\hat a_2$ meet on a beam splitter with mixing angle $\theta$ and phase $\varphi$, and $n_1$, $n_2$ photons enter. Without truncating the Fock space, what does the splitter do to the creation operators?

The splitter is $\hat B=e^{\hat G}$ with $\hat G=\theta\left(e^{i\varphi}\hat a_1\hat a_2^\dagger-e^{-i\varphi}\hat a_1^\dagger \hat a_2\right)$ and $[\hat a_j,\hat a_k^\dagger]=\delta_{jk}$. By $e^{\hat G}\hat a_k^\dagger e^{-\hat G}=\hat a_k^\dagger+[\hat G,\hat a_k^\dagger]+\frac1{2!}[\hat G,[\hat G,\hat a_k^\dagger]]+\dots$, if $[\hat G,\cdot\,]$ maps $(\hat a_1^\dagger,\hat a_2^\dagger)$ linearly by a matrix $M$, then $\hat a_k^\dagger\to\sum_j W_{jk}\hat a_j^\dagger$ with $W=e^{M}$. The vacuum stays put, so $|n_1,n_2\rangle$ leaves as $(\sum_j W_{j1}x_j)^{n_1}(\sum_j W_{j2}x_j)^{n_2}/\sqrt{n_1!\,n_2!}$ with $x_j$ standing for $\hat a_j^\dagger$, and $x_1^{t_1}x_2^{t_2}$ is $\sqrt{t_1!\,t_2!}\,|t_1,t_2\rangle$.

The commutators for general $\theta,\varphi$ give $e^{M}$, compared with the network's transfer matrix; the curves are the probabilities of $t_1$ photons in output 1 read off the polynomial, and the markers the Fock-space beam-splitter operator at $\varphi=\pi/3$:

```wl
algBSVars = FieldVariables[{\[FormalA], \[FormalB]}];
algBSGen = θ (Exp[I φ] algBSVars[[2, 1]] ** algBSVars[[1, 2]] - Exp[-I φ] algBSVars[[1, 1]] ** algBSVars[[2, 2]]);
algBSComm = BosonicNormalOrder[Commutator[algBSGen, #], algBSVars] & /@ algBSVars[[1]];
algBSW = FullSimplify[MatrixExp[Transpose[Coefficient[#, algBSVars[[1]]] & /@ algBSComm]]];
algBSProb[{n1_, n2_}] := With[{x = Array[\[FormalX], 2]}, KeySort @ Association[(First[#1] -> FullSimplify[Abs[#2]^2 Times @@ (#1!)/(n1! n2!), Element[{θ, φ}, Reals]]) & @@@
     CoefficientRules[Expand[(x . algBSW[[All, 1]])^n1 (x . algBSW[[All, 2]])^n2], x]]];
algBSFock[{n1_, n2_}] := With[{d = n1 + n2 + 1}, Transpose @ Table[With[{pr = ArrayReshape[BeamSplitterOperator[{t, π/3}, d][FockState[{n1, n2}, d]]["ProbabilitiesList"], {d, d}]},
      Table[{t, pr[[j + 1, n1 + n2 - j + 1]]}, {j, 0, n1 + n2}]], {t, 0., N[π/2], N[π/16]}]];
algBSPanel[in_] := Show[Plot[Evaluate[Values[algBSProb[in]]], {θ, 0, π/2}, PlotRange -> All, PlotLegends -> Table[Row[{Subscript["t", 1], " = ", j}], {j, 0, Total[in]}]],
    ListPlot[algBSFock[in], PlotMarkers -> Automatic], PlotRange -> {-0.02, 1.02}, Frame -> True, FrameTicks -> {{Automatic, None}, {{#, Style[TraditionalForm[#], FractionBoxOptions -> {Beveled -> True}]} & /@ Range[0, π/2, π/8], None}},
    FrameLabel -> {"θ", Row[{"P(", Subscript["t", 1], ")"}]}, PlotLabel -> Row[{"input  |", Subscript["n", 1], ", ", Subscript["n", 2], "⟩ = |", in[[1]], ", ", in[[2]], "⟩"}], ImageSize -> 300];
algBSPanels = algBSPanel /@ {{2, 1}, {3, 1}, {3, 2}};
algBSCoef[c_] := If[MatchQ[c, Times[-1, __]], Row[{"\[Minus]", TraditionalForm[-c]}], TraditionalForm[c]];
Labeled[Legended[PlotGrid[{First /@ algBSPanels}, ImageSize -> 900], legendBelow[Last[Last[algBSPanels]]]],
  Style[Row[{"[G, ", Superscript[Subscript["a", 1], "\[Dagger]"], "] = ", algBSCoef[Coefficient[algBSComm[[1]], algBSVars[[1, 2]]]], " ", Superscript[Subscript["a", 2], "\[Dagger]"],
     ",   [G, ", Superscript[Subscript["a", 2], "\[Dagger]"], "] = ", algBSCoef[Coefficient[algBSComm[[2]], algBSVars[[1, 1]]]], " ", Superscript[Subscript["a", 1], "\[Dagger]"],
     ",     ", Superscript["e", "M"], " = W of QuantumInterferometer:  ", algBSW === QuantumInterferometer[{"BS"[θ, φ]}]["Unitary"]}], 13, FontFamily -> "Helvetica"], Top]
```

The commutator turns each creation operator into a multiple of the other, so the whole series sums to the $2\times2$ transfer matrix and the splitter simply rotates $\hat a_1^\dagger,\hat a_2^\dagger$. The output amplitudes are then coefficients of a polynomial of degree $n_1+n_2$, and their squares need no value of $\varphi$ and sit on the Fock-space markers for each input shown, although no photon-number cutoff entered the algebra.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Balanced homodyne detection of squeezed light


A signal mode $\hat a$ in the displaced squeezed state $\hat D(\alpha)\hat S(r)|0\rangle$, with real $\alpha$ and $r$, meets a local oscillator, a coherent state $|\beta e^{i\varphi}\rangle$ in mode $\hat b$, on a balanced beam splitter, and two detectors count the photons $\hat N_c$ and $\hat N_d$ in its outputs. How strong must the local oscillator be for the difference count to show the signal's quadrature noise?

With $\hat S(r)=e^{r(\hat a^2-\hat a^{\dagger2})/2}$ and $\hat D(\alpha)=e^{\alpha(\hat a^\dagger-\hat a)}$, the splitter and the quadrature that the oscillator's phase $\varphi$ selects are

$$
\hat B=e^{\frac\pi4\left(\hat a\,\hat b^\dagger-\hat a^\dagger \hat b\right)},\qquad \hat X(\varphi)=\frac12\left(\hat a\,e^{-i\varphi}+\hat a^\dagger e^{i\varphi}\right),
$$

$\hat B$ turns the difference count, up to its sign, into $\hat a^\dagger \hat b+\hat b^\dagger \hat a$, and its square is normal-ordered; with $\hat b$ coherent this gives $\operatorname{Var}(\hat N_c-\hat N_d)=4|\beta|^2\operatorname{Var}\hat X(\varphi)+\langle \hat a^\dagger \hat a\rangle$, the second term coming from the signal's own photons beating with the oscillator's vacuum fluctuations.

With the two modes sent through $\hat B$ on 40 Fock levels each, the markers are $\operatorname{Var}(\hat N_c-\hat N_d)/4|\beta|^2$ from the joint photon-count distribution and the curves the closed form $\operatorname{Var}\hat X(\varphi)+\langle \hat a^\dagger \hat a\rangle/4|\beta|^2$ with $\langle \hat a^\dagger \hat a\rangle=\alpha^2+\sinh^2r$, with $\operatorname{Var}\hat X(\varphi)$ from the signal's covariance matrix dashed:

```wl
hdVars = FieldVariables[{\[FormalA], \[FormalB]}];
hdCount = hdVars[[1, 1]] ** hdVars[[2, 2]] + hdVars[[1, 2]] ** hdVars[[2, 1]];
hdLevels = 40; hdK = Range[0, hdLevels - 1]; hdDiff = Outer[Subtract, hdK, hdK];
hdSplitter = BeamSplitterOperator[{N[Pi/4], 0}, hdLevels];
hdNoise[sig_, β_][φ_] := With[{pr = ArrayReshape[hdSplitter[QuantumTensorProduct[sig, PhaseShiftOperator[φ, hdLevels] @ CoherentState[hdLevels][β]]]["ProbabilitiesList"], {hdLevels, hdLevels}]},
   (Total[hdDiff^2 pr, 2] - Total[hdDiff pr, 2]^2)/(4 β^2)];
hdPanel[α_, r_] := With[{sig = (DisplacementOperator[N[α], hdLevels] @ SqueezeOperator[N[r], hdLevels] @ FockState[0, hdLevels])["Normalize"], φs = Subdivide[0., N[Pi], 12]},
   With[{varX = With[{σ = CovarianceMatrix[sig]}, {Cos[#], Sin[#]} . σ . {Cos[#], Sin[#]}/2 &]}, Show[
     LogPlot[Evaluate[Append[Table[varX[φ] + (α^2 + Sinh[r]^2)/(4 β^2), {β, {0.5, 1, 3}}], varX[φ]]], {φ, 0, Pi}, PlotStyle -> {Automatic, Automatic, Automatic, Dashed},
      PlotLegends -> {"|\[Beta]| = 0.5", "|\[Beta]| = 1", "|\[Beta]| = 3", "Var X(\[CurlyPhi])"}],
     ListLogPlot[Table[Transpose[{φs, hdNoise[sig, β] /@ φs}], {β, {0.5, 1., 3.}}], PlotMarkers -> Automatic],
    PlotRange -> {{0, Pi}, Log[{0.04, 6}]}, GridLines -> {None, {Log[1/4]}}, Axes -> False,
    Frame -> True, FrameTicks -> {{{Log[#], #} & /@ {0.05, 0.1, 0.25, 0.5, 1, 2, 5}, None}, {Range[0, Pi, Pi/4], None}},
    FrameLabel -> {"\[CurlyPhi]", Row[{"Var(", Subscript["N", "c"], " - ", Subscript["N", "d"], ") / 4|\[Beta]|\.b2"}]},
    PlotLabel -> Row[{"\[Alpha] = ", α, ",   r = ", r}], ImageSize -> 300]]];
hdPanels = Outer[hdPanel, {0.5, 1}, {0, 0.4, 0.8}];
Labeled[Legended[Grid[Map[First, hdPanels, {2}]], legendBelow[Last[hdPanels[[1, 1]]]]],
  Style[Row[{"normal-ordered  ", Superscript[Row[{"(", Superscript["a", "\[Dagger]"], "b + ", Superscript["b", "\[Dagger]"], "a)"}], 2], "  =  ",
     TraditionalForm[BosonicNormalOrder[hdCount ** hdCount, hdVars] /. GeneralizedPower[NonCommutativeMultiply, o_, k_] :> Superscript[o, k] /. NonCommutativeMultiply -> (Row[{##}] &) /.
       {SuperDagger[\[FormalA]] -> Superscript["a", "\[Dagger]"], SuperDagger[\[FormalB]] -> Superscript["b", "\[Dagger]"], \[FormalA] -> "a", \[FormalB] -> "b"}]}], 13, FontFamily -> "Helvetica"], Top]
```

The markers lie on the closed form and above the dashed curve, so a weak oscillator buries the quadrature noise under the signal's own photons, the more so for the brighter, more displaced signal. The strongest oscillator dips below the coherent state's flat $1/4$ (grid line) near $\varphi=0$ and $\pi$ in every squeezed panel, and $|\beta|=1$ does so only for the dimmer, less squeezed signal. Even then the markers stay above $\operatorname{Var}\hat X(\varphi)$, because the signal's photons must be small next to $|\beta|^2$ times the squeezed variance, not merely next to $|\beta|^2$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Heralding a state from one photon and coherent light


One photon enters port 1 of a beam splitter with mixing angle $\theta$ and phase $\varphi$, and a coherent state $|\beta\rangle$ with complex $\beta$ enters port 2. A photon counter on output 1 registers exactly one photon. What state does that leave in output 2?

The splitter maps $\hat a_1^\dagger\to\cos\theta\,\hat a_1^\dagger+e^{i\varphi}\sin\theta\,\hat a_2^\dagger$ and splits the coherent amplitude into $-e^{-i\varphi}\beta\sin\theta$ in output 1 and $\beta\cos\theta$ in output 2. Projecting output 1 onto $|1\rangle$ keeps two terms, the photon transmitted with the coherent part giving no photon, or the photon reflected with the coherent part giving one, so the unnormalized state of output 2 is

$$
|\psi\rangle=e^{-|\beta|^2\sin^2\theta/2}\left(\cos\theta-\beta\sin^2\theta\,\hat a^\dagger\right)|\beta\cos\theta\rangle ,
$$

with the phase $\varphi$ cancelling between the two terms. Since $\hat a^\dagger|\gamma\rangle=\hat D(\gamma)(\hat a^\dagger+\gamma^*)|0\rangle$, this is a displaced superposition of $|0\rangle$ and $|1\rangle$,

$$
|\psi\rangle\propto \hat D(\beta\cos\theta)\left[\cos\theta\,(1-|\beta|^2\sin^2\theta)\,|0\rangle-\beta\sin^2\theta\,|1\rangle\right].
$$

The amplitudes $\langle 1,k|\hat U|1,\beta\rangle$ for $k\le14$ and general $\theta,\varphi,\beta$ are compared with $\langle k|\psi\rangle$; the panels show $\pi W$ of the normalized state built from them, at $\arg\beta=\pi/4$, with quadratures $\hat x,\hat p$ defined by $\hat a=(\hat x+i\hat p)/\sqrt2$, every panel on one color scale from $-1$ to $1$, the bounds of $\pi W$; the white line encloses $\pi W<-0.01$:

```wl
heraldAmps = Table[QuantumInterferometer[{"BS"[θ, φ]}]["Amplitude", {1, "Coherent"[β]} -> {1, k}], {k, 0, 14}];
heraldModel = Table[Exp[-Abs[β]^2/2] (Cos[θ] (β Cos[θ])^k - k β Sin[θ]^2 (β Cos[θ])^(k - 1))/Sqrt[k!], {k, 0, 14}];
heraldWigner[t_, b_] := WignerRepresentation[QuantumState[N[heraldAmps /. {θ -> t, β -> b Exp[I π/4]}], 15]["Normalized"], {-4.5, 4.5}, {-4.5, 4.5}, "GridSize" -> 121];
heraldColors = ColorData["Rainbow"][(1 + #)/2] &;
heraldPanel[t_, b_] := With[{w = heraldWigner[t, b]},
   Plot3D[π w[x, p], {x, -4.5, 4.5}, {p, -4.5, 4.5}, PlotRange -> {-1, 1}, ColorFunction -> (heraldColors[#3] &), ColorFunctionScaling -> False, PlotPoints -> 40, MaxRecursion -> 0,
    MeshFunctions -> {#3 &}, Mesh -> {{-0.01}}, MeshStyle -> Directive[White, Thick], BoxRatios -> {1, 1, 0.5}, AxesLabel -> {"x", "p", "πW"},
    PlotLabel -> Row[{"θ = ", t, ",   |β| = ", b}], ImageSize -> 280]];
heraldPanels = AssociationThread[{π/8, π/4, 3 π/8}, Outer[heraldPanel, {π/8, π/4, 3 π/8}, {0.5, 1, 2}]];
heraldLabel = Style[Row[{"QF:  ", TraditionalForm[FullSimplify[heraldAmps[[;; 2]], Element[θ, Reals]]], " for k = 0, 1;    all k ≤ 14 equal ⟨k|ψ⟩:  ",
     FullSimplify[heraldAmps - heraldModel, Element[θ, Reals]] === ConstantArray[0, 15]}], 13, FontFamily -> "Helvetica"];
Manipulate[Labeled[Legended[Grid[{heraldPanels[θv]}], BarLegend[{heraldColors, {-1, 1}}, LegendLabel -> "πW"]], heraldLabel, Top],
  Row[{Control[{{θv, π/4, "θ"}, Keys[heraldPanels], SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"θ = ", θv}], <|π/8 -> "Weak light looks like the coherent bump |β cos θ⟩.",
        π/4 -> "The negative dip deepens as |β| grows.",
        3 π/8 -> "|β| = 1 heralds nearly a single photon; at |β| = 2 the dip is shallower again."|>[θv]}}]]}], SaveDefinitions -> True]
```

Output 2 is a displaced superposition of $|0\rangle$ and $|1\rangle$, so its Wigner function is negative somewhere whenever the reflected photon contributes; with weak light and a small $\theta$ that negative region is too shallow to reach the white contour at $\pi W=-0.01$, and output 2 looks like the coherent bump $|\beta\cos\theta\rangle$. A larger $\theta$ or a brighter $\beta$ deepens it, beside a crescent of positive $W$, until $|\beta|\sin\theta$ reaches about 1, where the click heralds a displaced single photon; past that the vacuum part returns and the dip grows shallower again, as at $\theta=3\pi/8$, $|\beta|=2$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A Kerr medium in a Mach-Zehnder interferometer


A beam splitter $\hat B(\theta,\varphi)$ opens a Mach-Zehnder interferometer, arm 1 carries a linear phase $\alpha$ and a Kerr medium with nonlinear phase $\chi$, and the inverse splitter $\hat B^\dagger$ closes it. Either one photon enters port 1, or one photon enters each port. What does the Kerr phase change that the linear phase cannot?

The network is

$$
\hat U=\hat B(\theta,\varphi)^\dagger\,e^{i\chi\,\hat n_1(\hat n_1-1)/2}\,e^{i\alpha\,\hat n_1}\,\hat B(\theta,\varphi),\qquad \hat B(\theta,\varphi)=e^{\theta\left(e^{i\varphi}\hat a_1\hat a_2^\dagger-e^{-i\varphi}\hat a_1^\dagger \hat a_2\right)},
$$

with $\hat n_1=\hat a_1^\dagger \hat a_1$; the Kerr factor is 1 on zero or one photon in the arm. The two linear pieces become Fock-space circuits on three levels per mode, enough for two photons, with $\theta,\varphi,\alpha$ still general; the Kerr medium is a diagonal operator between them.

The figure plots, for general $\theta,\varphi,\alpha,\chi$, the probability that the input pattern comes back; the markers are the permanent formula for the linear network ($\chi=0$) at $\varphi=\pi/3$:

```wl
kerrMZIn = QuantumInterferometer[{"BS"[θ, φ], "PS"[α] -> 1}];
kerrMZOut = QuantumInterferometer[{"BS"[θ, φ]}]["Dagger"];
kerrMZKerr = QuantumOperator[DiagonalMatrix[SparseArray[Exp[I χ Range[0, 2] (Range[0, 2] - 1)/2]]], {1}, 3];
kerrMZCircuit = QuantumCircuitOperator[{kerrMZIn["CircuitOperator", 3], kerrMZKerr, kerrMZOut["CircuitOperator", 3]}];
kerrMZReturn[in_] := FullSimplify[Abs[(SuperDagger[FockState[in, 3]] @ kerrMZCircuit[FockState[in, 3]])["Scalar"]]^2, Element[{θ, φ, α, χ}, Reals]];
{kerrMZ1, kerrMZ2} = kerrMZReturn /@ {{1, 0}, {1, 1}};
kerrMZLinear = kerrMZOut[kerrMZIn]["Probability", # -> #] & /@ {{1, 0}, {1, 1}};
kerrMZPanel[t_] := Show[
    Plot[Evaluate[Prepend[Table[kerrMZ2 /. {θ -> t, χ -> c}, {c, {0, π/2, π}}], kerrMZ1 /. θ -> t]], {α, 0, 2 π}, PlotRange -> All, PlotStyle -> {Dashed, Automatic, Automatic, Automatic},
     PlotLegends -> {"one photon |1, 0⟩", "pair |1, 1⟩,  χ = 0", "pair |1, 1⟩,  χ = π/2", "pair |1, 1⟩,  χ = π"}],
    ListPlot[Table[{a, # /. {θ -> t, φ -> π/3, α -> a}}, {a, 0, 2 π, π/8}] & /@ kerrMZLinear, PlotMarkers -> Automatic],
    PlotRange -> {-0.02, 1.02}, Frame -> True, FrameTicks -> {{Automatic, None}, {{#, Style[TraditionalForm[#], FractionBoxOptions -> {Beveled -> True}]} & /@ Range[0, 2 π, π/2], None}},
    FrameLabel -> {"α", "return probability"}, PlotLabel -> Row[{"θ = ", t}], ImageSize -> 320];
kerrMZPanels = kerrMZPanel /@ {π/4, π/8};
Labeled[Legended[PlotGrid[{First /@ kerrMZPanels}, ImageSize -> 640], legendBelow[Last[First[kerrMZPanels]]]],
  Style[Column[{Row[{"P(|1, 0⟩ → |1, 0⟩) = ", TraditionalForm[kerrMZ1]}], Row[{"P(|1, 1⟩ → |1, 1⟩) = ", TraditionalForm[kerrMZ2]}]}], 13, FontFamily -> "Helvetica"], Top]
```

Neither return probability contains $\varphi$, and the one-photon fringe contains no $\chi$ either: a single photon in the arm picks up no Kerr phase, so the Kerr medium acts on photon pairs alone. On the balanced splitter the pair swings fully between returning and never returning, and $\chi$ only shifts that swing in $\alpha$. At $\theta=\pi/8$ the one-photon fringe no longer empties, and $\chi$ changes the shape of the pair fringe as well as its position.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Phase space and quantum optics

Cat, squeezed, thermal and Fock states come with Wigner and Husimi functions, covariance matrices and photon statistics.

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
   MeshFunctions -> {#3 &}, Mesh -> {{0}}, MeshStyle -> Directive[White, Thick], BoxRatios -> {2, 1, 0.6}, AxesLabel -> {"x", "p", "πW"}, ImageSize -> 340,
   PlotLabel -> Column[{Row[{"α = ", amp, ",   φ = ", StringReplace[ToString[phase, InputForm], "Pi" -> "π"]}], Row[{"computed πW(0,0) = ", NumberForm[π w[0, 0], {4, 3}, ExponentFunction -> (Null &)]}]}, Alignment -> Center]]];
catPanels = AssociationThread[{1, 2, 3}, Outer[catPanel, {1, 2, 3}, {0, π/2, π}]];
Manipulate[Labeled[Legended[Grid[{catPanels[amp]}], BarLegend[{catColors, {-1, 1}}, LegendLabel -> "πW"]],
  Style[Row[{"πW(0,0) = ⟨Π⟩ = ", TraditionalForm[catParity]}], 15, FontFamily -> "Helvetica"], Top], Row[{Control[{{amp, 2, "α"}, Keys[catPanels], SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"α = ", amp}], <|1 -> "The lobes still overlap, yet W already goes negative.",
        2 -> "The lobes pull apart; fringes between them dip below zero.",
        3 -> "The fringes crowd together; at φ = π/2 the origin sits near zero."|>[amp]}}]]}], SaveDefinitions -> True]
```

As $\alpha$ grows the two coherent lobes, at $x=\pm\sqrt2\,\alpha$, pull apart, and the fringes along $p$ between them dip below zero and crowd together; at $\alpha=1$ the lobes still overlap, yet $W$ already goes negative. In every panel the origin takes the closed-form parity: $\varphi=0$ and $\varphi=\pi$ give parity eigenstates at $\pm1$ for every $\alpha$, while at $\varphi=\pi/2$ the parity is $e^{-2\alpha^2}$, so the origin falls toward zero, between a positive and a negative fringe.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Squeezed thermal light: Wigner against Husimi


A field mode in thermal equilibrium, with mean photon number $\bar n$, is squeezed by a degenerate parametric amplifier with squeeze parameter $r$. How do its Wigner and Husimi ellipses compare with the vacuum's?

With $\hat x=(\hat a+\hat a^\dagger)/\sqrt2$ and $\hat p=(\hat a-\hat a^\dagger)/i\sqrt2$, the state is

$$
\rho=\hat S(r)\,\rho_{\rm th}\,\hat S^\dagger(r),\qquad \hat S(r)=e^{r\left(\hat a^2-\hat a^{\dagger2}\right)/2},\qquad \rho_{\rm th}=\sum_n\frac{\bar n^{\,n}}{(1+\bar n)^{n+1}}\,|n\rangle\langle n| .
$$

The state is Gaussian, so its Wigner function is the Gaussian $W\propto e^{-v\cdot\sigma^{-1}v/2}$ in $v=(x,p)$, fixed by the quadrature covariance matrix $\sigma$ ($\sigma=\frac12$ for the vacuum). The Husimi function $Q(x,p)=\langle\alpha|\rho|\alpha\rangle/2\pi$, with $\alpha=(x+ip)/\sqrt2$, is $W$ smoothed by the vacuum's Wigner function, so it is the Gaussian with $\sigma+\frac12$. Half the peak of either lies on the ellipse $v\cdot\sigma^{-1}v=2\ln2$. Here $\sigma=(\bar n+\frac12)\,{\rm diag}(e^{-2r},e^{2r})$, so the $W$ ellipse enters the vacuum's circle once $r>\frac12\ln(2\bar n+1)$.

$\sigma$ comes from the state on a Fock space cut at 140 levels. Each panel draws the half-peak ellipses of $W$ and $Q$ and, dashed, those of the vacuum:

```wl
sqLevels = 140;
sqCovariance[nb_, r_] := CovarianceMatrix[SqueezeOperator[N[r], sqLevels][ThermalState[N[nb], sqLevels]]["Normalized"]];
sqPanel[nb_, r_] := With[{cov = sqCovariance[nb, r], v = {x, p}},
  ContourPlot[Evaluate[{v . Inverse[cov] . v == 2 Log[2], v . Inverse[cov + IdentityMatrix[2]/2] . v == 2 Log[2], v . v == Log[2], v . v == 2 Log[2]}],
   {x, -5, 5}, {p, -5, 5}, ContourStyle -> {Automatic, Automatic, Dashed, Dashed}, PlotLegends -> {"Wigner W", "Husimi Q", "vacuum W", "vacuum Q"},
   FrameLabel -> {"x", "p"}, PlotLabel -> Row[{OverBar["n"], " = ", nb, ",   r = ", r}], ImageSize -> 260]];
sqPanels = AssociationThread[{0, 0.5, 1}, Outer[sqPanel, {0, 0.5, 1}, {0, 0.5, 1}]];
Manipulate[Legended[PlotGrid[{First /@ sqPanels[nb]}, ImageSize -> 780], legendBelow[Last[sqPanels[0][[1]]]]], Row[{Control[{{nb, 0.5, OverBar["n"]}, Keys[sqPanels], SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{OverBar["n"], " = ", nb}], <|0 -> "Any squeezing pushes W inside the vacuum circle.",
        0.5 -> "r = 0.5 already pushes W inside the vacuum circle.",
        1 -> "r = 0.5 stays just outside; only r = 1 pushes W inside."|>[nb]}}]]}], SaveDefinitions -> True]
```

Squeezing narrows both ellipses along $x$ and stretches them along $p$, while thermal noise widens them in every direction, so at $\bar n=1$ the Wigner ellipse at $r=0.5$ still sits just outside the vacuum's circle and only $r=1$ pushes it inside, as the threshold predicts. $Q$ is $W$ smoothed by the vacuum, so every $Q$ ellipse is its $W$ ellipse widened by one vacuum's worth of noise: it never enters the vacuum's $W$ circle, yet it enters the vacuum's own $Q$ circle exactly when $W$ enters the vacuum's $W$ circle. Only the width of $W$ is the quadrature noise itself.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Photon statistics behind a beam splitter


Light from a source holding $\langle\hat n\rangle$ photons on average (a photon-number state $|n\rangle$, a laser in a coherent state, or a thermal lamp) falls on a beam splitter of mixing angle $\theta$, with vacuum in the other port, and a detector counts the photons in each output. Do the joint counts tell the three sources apart, whatever the splitting?

The splitter mixes the source mode $\hat a$ with the empty mode $\hat b$, and $\hat N_c,\hat N_d$ count the photons in its two outputs. The splitter, the source's second-order coherence and the detectors' normalized joint count are

$$
\hat B(\theta)=e^{\theta\left(\hat a\,\hat b^\dagger-\hat a^\dagger \hat b\right)},\qquad g^{(2)}(0)=\frac{\langle \hat a^{\dagger2}\hat a^2\rangle}{\langle \hat a^\dagger \hat a\rangle^2},\qquad \frac{\langle \hat N_c\hat N_d\rangle}{\langle \hat N_c\rangle\langle \hat N_d\rangle} .
$$

With vacuum in $\hat b$, $n$ source photons split binomially, $P(k\mid n)=\binom nk\cos^{2k}\theta\,\sin^{2(n-k)}\theta$, so $\langle \hat N_c\hat N_d\rangle=\cos^2\theta\sin^2\theta\,\langle \hat a^{\dagger2}\hat a^2\rangle$ and $\langle \hat N_c\rangle\langle \hat N_d\rangle=\cos^2\theta\sin^2\theta\,\langle \hat a^\dagger \hat a\rangle^2$: for every splitting the joint count is $g^{(2)}(0)$.

The markers are $g^{(2)}(0)$ of each source, and the lines the joint counts from the source's photon-number distribution and the splitter's exact law $P(k\mid n)$, on 70 levels:

```wl
hbtLevels = 70;
hbtSources[nb_] := hbtSources[nb] = <|"Fock |n⟩" -> FockState[nb, hbtLevels], "coherent" -> CoherentState[hbtLevels][N[Sqrt[nb]]],
    "thermal" -> ThermalState[N[nb], hbtLevels]|>;
hbtLaw = QuantumInterferometer[{"BS"[θ]}]["Probability", {n, 0} -> {k, n - k}];
hbtMoments[th_] := hbtMoments[th] = Transpose @ Table[With[{pk = N[Table[hbtLaw, {k, 0, nn}] /. {n -> nn, θ -> th}], ks = Range[0, nn]},
    {ks . pk, (nn - ks) . pk, (ks (nn - ks)) . pk}], {nn, 0, hbtLevels - 1}];
hbtRatio[th_][src_] := With[{p = src["ProbabilitiesList"], m = hbtMoments[th]}, p . m[[3]]/((p . m[[1]]) (p . m[[2]]))];
hbtData[f_] := Transpose[Table[{nb, f[#]} & /@ Values[hbtSources[nb]], {nb, 4}]];
hbtMarkers = ListPlot[hbtData[G2Coherence], PlotMarkers -> Automatic, PlotLegends -> Keys[hbtSources[1]]];
Legended[Show[First[hbtMarkers], ListLinePlot[hbtData[hbtRatio[π/4]], PlotStyle -> Dashed],
    ListLinePlot[hbtData[hbtRatio[π/8]], PlotStyle -> Dotted], PlotRange -> {{0.8, 4.2}, {-0.05, 2.2}}, Frame -> True,
    FrameTicks -> {{Automatic, None}, {Range[4], None}}, FrameLabel -> {"⟨n⟩", Row[{Superscript["g", "(2)"], "(0),   ⟨", Subscript["N", "c"], Subscript["N", "d"], "⟩/⟨",
      Subscript["N", "c"], "⟩⟨", Subscript["N", "d"], "⟩"}]}, ImageSize -> 420],
  {Last[hbtMarkers], LineLegend[Directive[Black, #] & /@ {Dashed, Dotted}, {"joint counts, θ = π/4", "joint counts, θ = π/8"}]}]
```

The joint counts at both splittings lie on the source's $g^{(2)}(0)$. The photon-number state sits below 1, at $1-1/n$, the laser at 1 and the thermal lamp at 2, so the joint counts sort the sources by their statistics while the mean photon number cannot.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Collapse and revival of an atom in a coherent field


An excited two-level atom, frequency $\omega_a$, enters a cavity of frequency $\omega_c$ that holds a coherent field with mean photon number $\bar n$, and couples to it at strength $g$. Does the atom's Rabi oscillation survive the spread of photon numbers in the field?

With $\hbar=1$ and $\hat\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, the Jaynes-Cummings Hamiltonian $\omega_c \hat a^\dagger \hat a+\frac{\omega_a}{2}\hat\sigma_z+g(\hat\sigma_-\hat a^\dagger+\hat\sigma_+\hat a)$ becomes $\hat H=\frac{\Delta}{2}\hat\sigma_z+g(\hat\sigma_-\hat a^\dagger+\hat\sigma_+\hat a)$, $\Delta=\omega_a-\omega_c$, in the frame rotating at $\omega_c$. $\hat H$ couples $|e,n\rangle$ only to $|g,n+1\rangle$, so the propagator is a $2\times2$ block for each $n$, with off-diagonal element $g\sqrt{n+1}$. Exponentiating the block for general $n$, $g$, $\Delta$ and $t$ gives the inversion $w_n(t)=P_e-P_g$ from $|e,n\rangle$ in closed form, weighted with the Poisson distribution of the coherent field; time is in units of $1/g$:

```wl
jcBlock = {{Δ/2, g Sqrt[n + 1]}, {g Sqrt[n + 1], -Δ/2}};
jcU = FullSimplify[MatrixExp[-I t jcBlock], g > 0 && n >= 0 && Δ ∈ Reals];
jcInversion = FullSimplify[ComplexExpand[Abs[jcU[[1, 1]]]^2 - Abs[jcU[[2, 1]]]^2], g > 0 && n >= 0 && t > 0 && Δ ∈ Reals];
jcRecord[nbar_, Δv_] := With[{ms = Range[0, Ceiling[nbar + 8 Sqrt[nbar] + 10]]}, N[PDF[PoissonDistribution[nbar], ms]] . (jcInversion /. {n -> ms, g -> 1, Δ -> Δv})];
jcPanel[nbar_, Δv_] := Plot[Evaluate[jcRecord[nbar, Δv]], {t, 0, 60}, PlotPoints -> 2500, MaxRecursion -> 0, PlotRange -> {-1, 1},
    Frame -> True, FrameLabel -> {"g t", Row[{Subscript["P", "e"], " - ", Subscript["P", "g"]}]},
    PlotLabel -> Row[{OverBar["n"], " = ", nbar, ",   \[CapitalDelta]/g = ", Δv}], ImageSize -> 330];
jcPanels = AssociationMap[jcPanel @@ # &, Tuples[{{5, 20}, {0, 2, 4, 6}}]];
With[{formula = Style[Row[{Subscript["w", "n"], "(t) = ", TraditionalForm[jcInversion]}], 13, FontFamily -> "Helvetica"]},
  Manipulate[jcPanels[{nb, d}],
    Row[{Column[{Control[{{nb, 20, OverBar["n"]}, {5, 20}, SetterBar}], Control[{{d, 0, "\[CapitalDelta]/g"}, {0, 2, 4, 6}, SetterBar}]}], Spacer[40], Column[{formula, Dynamic[controlNotes[{
      {Row[{OverBar["n"], " = ", nb}], <|5 -> "Revivals smear into one another.",
        20 -> "Revivals stand apart cleanly."|>[nb]},
      {Row[{"Δ/g = ", d}], <|0 -> "Collapses about zero, then revives.",
        2 -> "Oscillates about a raised level; the revival comes later.",
        4 -> "A higher level and a later revival.",
        6 -> "The inversion stays near 1: the atom emits little."|>[d]}}]]}, Spacings -> 1]}],
    SaveDefinitions -> True]]
```

The oscillation does not survive at first: the photon numbers rotate the atom at different frequencies $\sqrt{\Delta^2+4g^2(n+1)}$, so the inversion collapses within a few $1/g$, and it revives when their phases realign. On resonance the revivals stand apart cleanly only in the larger field. Detuning lifts the level around which the inversion oscillates, because each $w_n$ is pulled toward 1, and delays the revival; the larger the detuning, the closer the inversion stays to 1 and the less the atom emits.

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

### Where a two-mode squeezed vacuum is squeezed


Two field modes $\hat a$ and $\hat b$ start in the vacuum and are squeezed together, $|\psi\rangle=\hat S_2(r)|0,0\rangle$. Where in the four-dimensional phase space is the state squeezed below the vacuum?

The two-mode squeezer $\hat S_2(r)=\exp{}[r(\hat a\hat b-\hat a^\dagger \hat b^\dagger)]$ creates and destroys photons in pairs, one in each mode. With $\hat x_k=(\hat a_k+\hat a_k^\dagger)/\sqrt2$, $\hat p_k=(\hat a_k-\hat a_k^\dagger)/i\sqrt2$ and $\hat{\mathbf r}=(\hat x_1,\hat p_1,\hat x_2,\hat p_2)$, the covariance matrix is $\sigma_{kl}=\frac12\langle\{\Delta\hat r_k,\Delta\hat r_l\}\rangle$, with $\sigma=\frac12\mathbb 1$ for the vacuum; here
$$
\sigma=\frac12\begin{pmatrix}c&0&-s&0\\0&c&0&s\\-s&0&c&0\\0&s&0&c\end{pmatrix},\qquad c=\cosh2r,\quad s=\sinh2r,
$$
so $\hat x_1+\hat x_2$ and $\hat p_1-\hat p_2$ have variance $e^{-2r}$, against 1 for the vacuum. In each plane of two quadratures the Wigner marginal falls to $e^{-1/2}$ of its peak on the ellipse $\mathbf r^T\sigma^{-1}\mathbf r=1$ built from that plane's $2\times2$ block of $\sigma$.

$\hat S_2(r)$, a matrix exponential on 30 Fock levels per mode, acts on $|0,0\rangle$, and the covariance matrix $\sigma$ of the result gives the markers; each panel draws the closed-form ellipse in one plane (curves), with the vacuum's circle dashed:

```wl
tmsA = AnnihilationOperator[30, {1}]; tmsB = AnnihilationOperator[30, {2}];
tmsSigma[r_] := tmsSigma[r] = CovarianceMatrix[MatrixExp[r (tmsA @ tmsB - SuperDagger[tmsA] @ SuperDagger[tmsB]), FockState[{0, 0}, 30]]];
tmsClosed[r_] := {{Cosh[2 r], 0, -Sinh[2 r], 0}, {0, Cosh[2 r], 0, Sinh[2 r]}, {-Sinh[2 r], 0, Cosh[2 r], 0}, {0, Sinh[2 r], 0, Cosh[2 r]}}/2;
tmsEllipse[σ_] := Table[MatrixPower[σ, 1/2] . {Cos[u], Sin[u]}, {u, 0, 2 Pi, Pi/60}];
tmsPanel[{i_, j_}] := With[{labels = {Subscript["x", 1], Subscript["p", 1], Subscript["x", 2], Subscript["p", 2]}},
   Show[ListLinePlot[Append[tmsEllipse[tmsClosed[#][[{i, j}, {i, j}]]] & /@ {0.3, 0.6, 1.}, tmsEllipse[IdentityMatrix[2]/2]],
     PlotStyle -> {Automatic, Automatic, Automatic, Dashed}, PlotLegends -> {"r = 0.3", "r = 0.6", "r = 1", "vacuum"}],
    ListPlot[tmsEllipse[tmsSigma[#][[{i, j}, {i, j}]]][[;; ;; 10]] & /@ {0.3, 0.6, 1.}, PlotMarkers -> {Automatic, 7}],
    AspectRatio -> 1, PlotRange -> {{-2.5, 2.5}, {-2.5, 2.5}}, Frame -> True, FrameLabel -> labels[[{i, j}]],
    PlotLabel -> Row[{"(", labels[[i]], ", ", labels[[j]], ") plane"}], ImageSize -> 280]];
tmsPanels = tmsPanel /@ {{1, 2}, {1, 3}, {2, 4}};
Legended[Row[First /@ tmsPanels, Spacer[10]], legendBelow[Last[First[tmsPanels]]]]
```

On its own, mode 1 is a circle wider than the vacuum's, an isotropic Gaussian, so it is thermal light, noisier for larger $r$. The squeezing lives only in joint quadratures: $\hat x_1$ and $\hat x_2$ are anticorrelated and $\hat p_1$ and $\hat p_2$ correlated, so $\hat x_1+\hat x_2$ and $\hat p_1-\hat p_2$ fall inside the vacuum circle along the narrow axes of the tilted ellipses, further for larger $r$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Entanglement and information, proved

Quantifier elimination (`Reduce`, `Resolve`) turns statements about every state, every measurement or every parameter into proofs, and QF supplies the states, channels and partial transposes they are about.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### GHZ and W: two kinds of three-qubit entanglement


Three qubits share $|\mathrm{GHZ}\rangle=(|000\rangle+|111\rangle)/\sqrt2$ or $|W\rangle=(|001\rangle+|010\rangle+|100\rangle)/\sqrt3$. Each party may apply any invertible local map $A$, $B$, $C$, even one that succeeds only with some probability. Can such maps turn GHZ into W?

The maps act as $|\psi\rangle\mapsto(A\otimes B\otimes C)|\psi\rangle$, followed by normalization. Cut the amplitudes $\psi_{ijk}$ into the two $2\times2$ slices $M_i=(\psi_{ijk})_{jk}$; the Cayley hyperdeterminant $\operatorname{Det}\psi$ is the discriminant of the quadratic $\det(M_0+xM_1)$, and the three-tangle of a normalized state is $\tau_3=4|\operatorname{Det}\psi|$. The family

$$
|\psi_\varepsilon\rangle\propto A_\varepsilon^{\otimes3}|\mathrm{GHZ}\rangle\propto(|0\rangle+\varepsilon|1\rangle)^{\otimes3}-|000\rangle,\qquad A_\varepsilon=\begin{pmatrix}1&-1\\ \varepsilon&0\end{pmatrix},\quad\det A_\varepsilon=\varepsilon,
$$

is GHZ moved by an invertible map for every $\varepsilon>0$.

Transforming both states with all twelve entries of $A,B,C$ general gives $\tau_3$ and the infidelity $1-|\langle W|\psi_\varepsilon\rangle|^2$ in closed form:

```wl
hyperdetSLOCC[v_] := With[{m = ArrayReshape[Normal[v], {2, 2, 2}]}, With[{c = CoefficientList[Det[m[[1]] + \[FormalX] m[[2]]], \[FormalX], 3]}, c[[2]]^2 - 4 c[[1]] c[[3]]]];
{ghzSLOCC, wSLOCC} = QuantumState[#[3]] & /@ {"GHZ", "W"};
mapSLOCC = QuantumTensorProduct @@ (QuantumOperator[Array[#, {2, 2}]] & /@ {slA, slB, slC});
detsSLOCC = Times @@ (Det[Array[#, {2, 2}]] & /@ {slA, slB, slC});
ghzMapSLOCC = Normal[mapSLOCC[ghzSLOCC]["StateVector"]];
{ghzRuleQ, wRuleQ} = {Expand[hyperdetSLOCC[ghzMapSLOCC] - detsSLOCC^2/4] === 0, Expand[hyperdetSLOCC[mapSLOCC[wSLOCC]["StateVector"]]] === 0};
vecSLOCC = ghzMapSLOCC /. Flatten[Thread[Flatten[Array[#, {2, 2}]] -> {1, -1, \[Epsilon], 0}] & /@ {slA, slB, slC}];
{tau3SLOCC, infidSLOCC} = Simplify[{4 Abs[hyperdetSLOCC[vecSLOCC]]/(vecSLOCC . vecSLOCC)^2, 1 - (Normal[wSLOCC["StateVector"]] . vecSLOCC)^2/(vecSLOCC . vecSLOCC)}, \[Epsilon] > 0];
LogLogPlot[{tau3SLOCC, infidSLOCC}, {\[Epsilon], 10^-3, 1}, PlotRange -> All, PlotStyle -> {Automatic, Dashed},
 PlotLegends -> {Subscript["\[Tau]", 3], "1 - |\[LeftAngleBracket]W|\[Psi]\[RightAngleBracket]|²"}, Frame -> True, FrameLabel -> {"\[Epsilon]", None}, ImageSize -> 420,
 PlotLabel -> Style[Column[{Row[{"Det((A\[CircleTimes]B\[CircleTimes]C)GHZ) = (det A det B det C)²/4: ", ghzRuleQ}], Row[{"Det((A\[CircleTimes]B\[CircleTimes]C)W) = 0: ", wRuleQ}],
     Row[{Subscript["\[Tau]", 3], " = ", TraditionalForm[tau3SLOCC], ",   1 - |\[LeftAngleBracket]W|\[Psi]\[RightAngleBracket]|² = ", TraditionalForm[infidSLOCC]}]}], 13, FontFamily -> "Helvetica"]]
```

An invertible local map multiplies the hyperdeterminant by $(\det A\det B\det C)^2\neq0$, so $\tau_3$ of W stays zero and $\tau_3$ of GHZ stays positive, and no such map links the two in either direction. The family $\psi_\varepsilon$ is GHZ-class at every $\varepsilon>0$, yet its $\tau_3$ and its infidelity with W both vanish as $\varepsilon^2$: W sits on the boundary of the GHZ class, reachable only as a limit.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Monogamy of entanglement: the CKW inequality


Three qubits $A$, $B$, $C$ share $|\psi_{p,\phi}\rangle=\sqrt p\,|\mathrm{GHZ}\rangle+e^{i\phi}\sqrt{1-p}\,|W\rangle$, a superposition of the two orthogonal states $(|000\rangle+|111\rangle)/\sqrt2$ and $(|001\rangle+|010\rangle+|100\rangle)/\sqrt3$. How much of qubit $A$'s entanglement with the pair $BC$ can it hold with $B$ and with $C$ separately?

$A$'s entanglement with $BC$ is the tangle $\tau_{A|BC}=C_{A|BC}^2=4\det\rho_A$. What it shares with $B$ alone is the Wootters concurrence $C_{AB}$ of the two-qubit state $\rho_{AB}=\operatorname{Tr}_C|\psi\rangle\langle\psi|$, and likewise $C_{AC}$. The Coffman-Kundu-Wootters inequality says

$$
\tau_{A|BC}\;\ge\;C_{AB}^2+C_{AC}^2,\qquad \tau_3=\tau_{A|BC}-C_{AB}^2-C_{AC}^2\ge0 ,
$$

and for a pure state the residual is the three-tangle $\tau_3=4|\operatorname{Det}\psi|$, with $\operatorname{Det}$ the Cayley hyperdeterminant of the amplitudes: the residual tangle that the pairwise concurrences do not account for.

$\tau_{A|BC}$ comes from $\rho_A$ in closed form, the pairwise concurrences are the markers, and $\tau_3$ comes from the amplitudes:

```wl
{ghzCKW, wCKW} = Normal[QuantumState[#[3]]["StateVector"]] & /@ {"GHZ", "W"};
vecCKW[p_, \[Phi]_] := Sqrt[p] ghzCKW + Exp[I \[Phi]] Sqrt[1 - p] wCKW;
hyperdetCKW[v_] := With[{m = ArrayReshape[v, {2, 2, 2}]}, With[{c = CoefficientList[Det[m[[1]] + \[FormalX] m[[2]]], \[FormalX], 3]}, c[[2]]^2 - 4 c[[1]] c[[3]]]];
tangleCKW = Simplify[ComplexExpand[4 Det[Normal[QuantumPartialTrace[QuantumState[vecCKW[\[FormalP], \[FormalPhi]], {2, 2, 2}], {2, 3}]["DensityMatrix"]]]], 0 < \[FormalP] < 1];
tau3CKW = 4 Abs[hyperdetCKW[vecCKW[\[FormalP], \[FormalPhi]]]];
pairsCKW[p_, \[Phi]_] := Total[QuantumEntanglementMonotone[QuantumPartialTrace[QuantumState[vecCKW[p, \[Phi]], {2, 2, 2}], {#}], "Concurrence"]^2 & /@ {3, 2}];
panelCKW[\[Phi]_] := Show[
   Plot[Evaluate[{tangleCKW - tau3CKW, tangleCKW, tau3CKW} /. \[FormalPhi] -> \[Phi]], {\[FormalP], 0, 1}, PlotStyle -> {Dashed, Automatic, Dotted}, PlotRange -> All,
    PlotLegends -> {Row[{Subsuperscript["C", "AB", 2], " + ", Subsuperscript["C", "AC", 2], ": QF (markers), ", Subscript["\[Tau]", "A|BC"], " - ", Subscript["\[Tau]", 3], " (curve)"}], Subscript["\[Tau]", "A|BC"], Row[{Subscript["\[Tau]", 3], " = 4|Det \[Psi]|"}]}],
   ListPlot[Table[{p, pairsCKW[p, \[Phi]]}, {p, 0., 1., 0.1}], PlotMarkers -> {"\[EmptyCircle]", 12}],
   PlotRange -> {{0, 1}, {-0.02, 1.02}}, Frame -> True, FrameLabel -> {"p", None}, ImageSize -> 300, ImagePadding -> {{Automatic, Automatic}, {Automatic, 30}}, PlotLabel -> Style[Row[{"\[Phi] = ", \[Phi]}], 13, FontFamily -> "Helvetica"]];
panelsCKW = panelCKW /@ {0, Pi/2, Pi};
Legended[PlotGrid[{First /@ panelsCKW}, ImageSize -> 900], legendBelow[Last[First[panelsCKW]]]]
```

The pairwise sums fall on $\tau_{A|BC}-\tau_3$ at every sampled weight and phase, below the tangle: $A$ cannot hold more with $B$ and $C$ one at a time than with the pair, and the shortfall is exactly the three-tangle. W ($p=0$) spends all of it on pairs and GHZ ($p=1$) none. At $\phi=\pi$ the three-tangle touches zero at one intermediate weight, a superposition whose three-tangle vanishes as W's does, although, like W, it is entangled across every cut.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Which pure states LOCC can convert: Nielsen's majorization


Alice and Bob share a pure state and may use local operations and classical communication (LOCC). Two qutrit pairs have Schmidt weights $\lambda^\psi=(\frac12,\frac3{10},\frac15)$ and $\lambda^\phi=(\frac35,\frac15,\frac15)$, and two ququart pairs have $\lambda^\psi=(\frac25,\frac25,\frac1{10},\frac1{10})$ and $\lambda^\phi=(\frac12,\frac14,\frac14,0)$; a qubit pair $\chi$ with weights $(\frac35,\frac25)$ can be borrowed and must be returned. Which conversions succeed with certainty?

Write each state as $|\psi\rangle=\sum_k\sqrt{\lambda_k}\,|k\rangle_A|k\rangle_B$ with $\lambda_1\ge\lambda_2\ge\dots$. Nielsen's theorem: $|\psi\rangle\to|\phi\rangle$ by LOCC with certainty if and only if $\lambda^\psi$ is majorized by $\lambda^\phi$,

$$
\lambda^\psi\prec\lambda^\phi\iff\sum_{i\le k}\lambda^\psi_i\le\sum_{i\le k}\lambda^\phi_i\quad\text{for every }k .
$$

With $\chi$ borrowed, the question is $|\psi\rangle|\chi\rangle\to|\phi\rangle|\chi\rangle$, where Alice holds both first halves.

Each state is prepared with the discrete Fourier transform on both sides, so its Schmidt basis is hidden; the weights come back exactly as the eigenvalues of Alice's reduced state. The first panel shows, for every ordered qutrit target $(\lambda_1,\lambda_2,1-\lambda_1-\lambda_2)$, whether $\psi$ converts to it or it converts to $\psi$; the others plot the partial sums $\sum_{i\le k}\lambda_i$ against $k$, and $\psi\to\phi$ when the curve of $\phi$ is never below that of $\psi$:

```wl
schmidtNielsen[w_] := With[{f = FourierMatrix[Length[w]]}, QuantumState[KroneckerProduct[f, f] . Flatten[DiagonalMatrix[Sqrt[w]]], {1, 1} Length[w]]];
aliceNielsen[s_] := ReverseSort[Simplify[QuantumPartialTrace[s, Range[2, s["Qudits"], 2]]["Eigenvalues"]]];
majorizedNielsen[x_, y_] := And @@ Thread[Accumulate[x] <= Accumulate[y]];
{psiNielsen, phiNielsen} = aliceNielsen[schmidtNielsen[#]] & /@ {{1/2, 3/10, 1/5}, {3/5, 1/5, 1/5}};
regionNielsen = RegionPlot[Evaluate[(l1 >= l2 >= 1 - l1 - l2 >= 0 && #) & /@ {majorizedNielsen[psiNielsen, {l1, l2, 1 - l1 - l2}], majorizedNielsen[{l1, l2, 1 - l1 - l2}, psiNielsen]}],
   {l1, 1/3, 1}, {l2, 0, 1/2}, PlotPoints -> 60, PlotLegends -> {"\[Psi] \[RightArrow] target", "target \[RightArrow] \[Psi]"},
   Epilog -> {Black, {Dashed, Line[{{1/3, 1/3}, {1/2, 1/2}, {1, 0}, {1/3, 1/3}}]}, PointSize[Large], Point[{psiNielsen[[;; 2]], phiNielsen[[;; 2]]}], Text[Style["\[Psi]", 14], psiNielsen[[;; 2]], {-2, 0}], Text[Style["\[Phi]", 14], phiNielsen[[;; 2]], {-2, 0}]},
   Frame -> True, FrameLabel -> {Subscript["\[Lambda]", 1], Subscript["\[Lambda]", 2]}, ImageSize -> 300, PlotLabel -> Style["qutrit targets", 13, FontFamily -> "Helvetica"]];
catalystNielsen = schmidtNielsen[{3/5, 2/5}];
pairsNielsen = With[{w = {{2/5, 2/5, 1/10, 1/10}, {1/2, 1/4, 1/4, 0}}}, {schmidtNielsen /@ w, QuantumTensorProduct[schmidtNielsen[#], catalystNielsen] & /@ w}];
panelNielsen[{s1_, s2_}, title_] := With[{w1 = aliceNielsen[s1], w2 = aliceNielsen[s2]},
   ListLinePlot[{Accumulate[Prepend[w1, 0]], Accumulate[Prepend[w2, 0]]}, DataRange -> {0, Length[w1]}, PlotMarkers -> {Automatic, 9}, PlotStyle -> {Automatic, Dashed},
    PlotRange -> {0, 1.02}, Frame -> True, FrameTicks -> {{Automatic, None}, {Range[0, Length[w1]], None}}, FrameLabel -> {"k", Row[{"\[Sum]", Subscript["\[Lambda]", "i \[LessEqual] k"]}]}, ImageSize -> 300, PlotLegends -> {"\[Psi]", "\[Phi]"},
    PlotLabel -> Style[Column[{title, Row[{"\[Psi] \[RightArrow] \[Phi]: ", majorizedNielsen[w1, w2], ",   \[Phi] \[RightArrow] \[Psi]: ", majorizedNielsen[w2, w1]}]}, Alignment -> Center], 13, FontFamily -> "Helvetica"]]];
panelsNielsen = MapThread[panelNielsen, {pairsNielsen, {"ququart pairs", Row[{"ququart pairs with ", "\[Chi]"}]}}];
Grid[{{regionNielsen, Legended[Row[First /@ panelsNielsen], legendBelow[Last[First[panelsNielsen]]]]}}]
```

Inside the dashed triangle of ordered qutrit weights, every target in the right-hand region, $\phi$ among them, is reachable from $\psi$, the targets in the left-hand region reach $\psi$, and the rest are incomparable: LOCC can only make the weights steeper. The ququart curves cross, so neither state converts to the other, yet with $\chi$ borrowed the crossing disappears and $|\psi\rangle|\chi\rangle$ converts to $|\phi\rangle|\chi\rangle$, which hands $\chi$ back unchanged.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Bell-CHSH violation and the Tsirelson bound for a family of states


Two qubits share $\rho=p\,|\psi_\theta\rangle\langle\psi_\theta|+(1-p)\,I/4$, a partly entangled pure state $|\psi_\theta\rangle=\cos\theta|00\rangle+\sin\theta|11\rangle$ mixed with white noise. Alice measures $\vec a_1\cdot\hat{\vec\sigma}$ or $\vec a_2\cdot\hat{\vec\sigma}$, Bob $\vec b_1\cdot\hat{\vec\sigma}$ or $\vec b_2\cdot\hat{\vec\sigma}$, all with outcomes $\pm1$. For which $p$ and $\theta$ can the correlations violate the CHSH inequality?

Local models obey the CHSH inequality $|\langle\mathcal B\rangle|\le2$ for

$$
\mathcal B=\vec a_1\cdot\hat{\vec\sigma}\otimes(\vec b_1+\vec b_2)\cdot\hat{\vec\sigma}+\vec a_2\cdot\hat{\vec\sigma}\otimes(\vec b_1-\vec b_2)\cdot\hat{\vec\sigma} ,
$$

while quantum states reach at most $2\sqrt2$ (Tsirelson). With the correlation matrix $T_{ij}=\operatorname{Tr}\rho\,\hat\sigma_i\otimes\hat\sigma_j$, the Horodecki theorem gives the best value over all settings, $\max\langle\mathcal B\rangle=2\sqrt{t_1+t_2}$, where $t_1,t_2$ are the two largest eigenvalues of $T^{\mathsf T}T$. For two qubits $\rho$ is entangled exactly when its partial transpose has a negative eigenvalue (Peres-Horodecki).

For general $p$ and $\theta$ the optimal $\langle\mathcal B\rangle$ is in closed form; the figure maps, over every $(\theta,p)$, where the partial transpose has a negative eigenvalue, where the optimal settings violate CHSH, and where the settings optimal for a Bell state ($\vec a_{1,2}=\hat z,\hat x$, $\vec b_{1,2}=(\hat z\pm\hat x)/\sqrt2$) do:

```wl
rhoCHSH = \[FormalP] QuantumState[{Cos[\[FormalTheta]], 0, 0, Sin[\[FormalTheta]]}, 2]["MatrixState"] + (1 - \[FormalP]) QuantumState[IdentityMatrix[4]/4, {2, 2}];
corrCHSH = Simplify[Table[Tr[Normal[rhoCHSH["DensityMatrix"]] . KroneckerProduct[PauliMatrix[i], PauliMatrix[j]]], {i, 3}, {j, 3}], 0 < \[FormalTheta] < Pi/2];
bestCHSH = With[{ev = Eigenvalues[Transpose[corrCHSH] . corrCHSH]}, Simplify[2 Sqrt[Total[ev] - Min[ev]], 0 < \[FormalP] <= 1 && 0 < \[FormalTheta] < Pi/2]];
bellCHSH = Sqrt[2] (corrCHSH[[3, 3]] + corrCHSH[[1, 1]]);
ptCHSH = Simplify[QuantumPartialTranspose[rhoCHSH, {2}]["Eigenvalues"], 0 < \[FormalTheta] < Pi/2];
RegionPlot[Evaluate[{Min[ptCHSH] < 0, bestCHSH > 2, bellCHSH > 2}], {\[FormalTheta], 0, Pi/2}, {\[FormalP], 0, 1}, PlotPoints -> 50,
 PlotLegends -> {"entangled (negative partial transpose)", "violates CHSH, best settings", "violates CHSH, Bell-state settings"},
 Frame -> True, FrameLabel -> {"\[Theta]", "p"}, FrameTicks -> {{Automatic, None}, {Range[0, Pi/2, Pi/8], None}}, ImageSize -> 420,
 PlotLabel -> Style[Row[{"max\[LeftAngleBracket]\[ScriptCapitalB]\[RightAngleBracket] = ", TraditionalForm[bestCHSH /. {\[FormalP] -> p, \[FormalTheta] -> \[Theta]}]}], 13, FontFamily -> "Helvetica"]]
```

Every violating state is entangled, but below the violating region lies a band of entangled states that no CHSH settings detect, and as $\theta\to0$ it fills almost all of the entangled region. The optimal value reaches the Tsirelson value $2\sqrt2$ only for the Bell state at $p=1$. The Bell-state settings find only part of the violating region once $\theta\neq\pi/4$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A qutrit beats every noncontextual model


Five projectors $\hat\Pi_k=|v_k\rangle\langle v_k|$ on a qutrit form a pentagon: each can be measured jointly with its two neighbors, and no two neighbors both answer yes. If every test had a predetermined answer, whatever it is measured with, then $\sum_k\langle\hat\Pi_k\rangle\le2$ (Klyachko, Can, Binicioğlu and Shumovsky); the qutrit is in $\rho=p\,|\psi_\alpha\rangle\langle\psi_\alpha|+(1-p)\,\mathbb 1/3$, a pure state tilted by $\alpha$ from $|2\rangle$ mixed with white noise. Which of these states break the bound?

With $|0\rangle,|1\rangle,|2\rangle$ as the $x,y,z$ axes of a real space, the five rays sit on a cone around $|2\rangle$, one turned from the next by $4\pi/5$:

$$
v_k=R_z\!\left(\frac{4\pi k}{5}\right)\left(\sin\vartheta,\,0,\,\cos\vartheta\right),\qquad \cos^2\vartheta=\frac1{\sqrt5}\;\Longleftrightarrow\; v_k\cdot v_{k+1}=0,\qquad |\psi_\alpha\rangle=\sin\alpha\,|0\rangle+\cos\alpha\,|2\rangle .
$$

With $\hat K=\sum_k\hat\Pi_k$ and neighboring rays checked orthogonal, the figure maps where $\langle\hat K\rangle=\mathrm{Tr}\,\rho\hat K$ exceeds the noncontextual bound over every tilt and noise weight:

```wl
kcbsRays = Table[RotationMatrix[4 Pi k/5, {0, 0, 1}] . {Sqrt[1 - 1/Sqrt[5]], 0, 5^(-1/4)}, {k, 0, 4}];
kcbsExclusiveQ = FullSimplify[MapThread[Dot, {kcbsRays, RotateLeft[kcbsRays]}]] === ConstantArray[0, 5];
kcbsK = Total[KroneckerProduct[#, #] & /@ kcbsRays];
kcbsMean = FullSimplify[Tr[kcbsK . (\[FormalP] Outer[Times, #, #] &[{Sin[\[FormalAlpha]], 0, Cos[\[FormalAlpha]]}] + (1 - \[FormalP]) IdentityMatrix[3]/3)]];
RegionPlot[kcbsMean > 2, {\[FormalAlpha], 0, Pi/2}, {\[FormalP], 0, 1}, PlotPoints -> 50, Frame -> True, FrameLabel -> {"\[Alpha]", "p"},
 FrameTicks -> {{Automatic, None}, {Range[0, Pi/2, Pi/8], None}}, ImageSize -> 400, PlotLegends -> {Row[{"\[LeftAngleBracket]K\[RightAngleBracket] > 2"}]},
 PlotLabel -> Style[Column[{Row[{"neighbors orthogonal: ", kcbsExclusiveQ}], Row[{"\[LeftAngleBracket]K\[RightAngleBracket] = ", TraditionalForm[kcbsMean /. {\[FormalP] -> p, \[FormalAlpha] -> \[Alpha]}]}]}, Alignment -> Center], 13, FontFamily -> "Helvetica"]]
```

Only states close to the cone axis, with little noise, break the bound; the axis state tolerates the most noise, and tilting or mixing pulls $\langle\hat K\rangle$ down. A state tilted far enough toward the plane perpendicular to the axis stays below 2 even when pure.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Qutrit magic in the discrete Wigner function: Clifford gates only move it


A qutrit state maps to a quasi-probability $W(q,p)$ on the nine points of a $3\times3$ phase space. For odd dimension a pure state has $W\ge0$ exactly when it is a stabilizer state (Gross), so the sum-negativity $\mathrm{sn}=\sum_{W<0}|W|$ measures magic. Starting from the stabilizer state $|+\rangle$ and from $(|1\rangle-|2\rangle)/\sqrt2$, which of the Clifford gates $F$, $XZ$ and the non-Clifford gate $T$ changes the negativity?

With $\omega=e^{2\pi i/3}$, $\zeta=e^{2\pi i/9}$ and indices mod 3, the gates are

$$
X|j\rangle=|j+1\rangle,\quad Z|j\rangle=\omega^{-j}|j\rangle,\quad F|j\rangle=\frac1{\sqrt3}\sum_k\omega^{jk}|k\rangle,\quad T=\mathrm{diag}\!\left(1,\zeta,\zeta^{-1}\right),\qquad |+\rangle=F|0\rangle .
$$

After each gate $W$ is read from the state; each panel shows the nine values of $W$, negative ones in bold and framed, with its sum-negativity:

```wl
magicFmt[x_] := NumberForm[x, {3, 2}];
magicStates = <|TraditionalForm[Ket[{"+"}]] -> QuantumState[{1, 1, 1}/Sqrt[3], 3],
    TraditionalForm[HoldForm[Divide[Ket[{1}] - Ket[{2}], Sqrt[2]]]] -> QuantumState[{0, 1, -1}/Sqrt[2], 3]|>;
magicGates = <|"no gate" -> QuantumOperator["I"[3]], "F" -> QuantumOperator["H"[3]],
    "XZ" -> QuantumOperator["X"[3]] @ QuantumOperator["Z"[3]], "T" -> QuantumOperator[DiagonalMatrix[{1, Exp[2 Pi I/9], Exp[-2 Pi I/9]}], {1}, 3]|>;
qutritMagicPanel[gate_, state_] := With[{w = Re[N[Normal[magicGates[gate][magicStates[state]]["PhaseSpace"]]]]}, Style[
   Labeled[Grid[Append[Reverse[MapThread[Prepend, {Map[If[# < -10^-9, Item[Style[magicFmt[#], Bold], Frame -> True, FrameStyle -> Thick], magicFmt[#]] &, w, {2}], Range[0, 2]}]],
      {"p\\q", 0, 1, 2}], Dividers -> {{2 -> True}, {-2 -> True}}, ItemSize -> {3.5, 1.6}, Alignment -> Center],
    Column[{Row[{gate, ",  ", state}], Row[{"sn = ", magicFmt[Total[Select[Flatten[-w], Positive]]]}]}, Center], Top, LabelStyle -> {FontFamily -> "Helvetica", 13}],
   FontFamily -> "Helvetica", 13]];
Grid[Outer[qutritMagicPanel, Keys[magicGates], Keys[magicStates]]\[Transpose], Spacings -> {3, 2}]
```

The two Clifford gates only rearrange the nine numbers, so the sum-negativity never changes: $|+\rangle$ stays nonnegative, and the single negative entry of $(|1\rangle-|2\rangle)/\sqrt2$ moves to another phase-space point under $XZ$. $T$ changes the numbers themselves: it gives $|+\rangle$ three negative entries, creating magic from a stabilizer state, and it changes the negativity of the other state too.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Two qutrits entangled with a positive partial transpose


Five orthonormal product states of two qutrits, the Tiles basis, leave a four-dimensional subspace, and $\rho$ is the normalized projector onto it; $\rho_p=(1-p)\,\rho+p\,I/9$ adds white noise of weight $p$. Its partial transpose has no negative eigenvalue, so the Peres test is silent. Is $\rho$ nevertheless entangled?

With qutrit states $|0\rangle,|1\rangle,|2\rangle$, the Tiles states and the state are

$$
|0\rangle|0{-}1\rangle,\ \ |2\rangle|1{-}2\rangle,\ \ |0{-}1\rangle|2\rangle,\ \ |1{-}2\rangle|0\rangle,\ \ |0{+}1{+}2\rangle|0{+}1{+}2\rangle,\qquad \rho=\frac14\left(I-\sum_{k=1}^{5}|\psi_k\rangle\langle\psi_k|\right),
$$

with $|i{\pm}j\rangle$ normalized. A separable $\rho$ would be a mixture of product states in its range, that is, product states orthogonal to all five $|\psi_k\rangle$. The realignment test $\|R(\rho)\|_1>1$, with $R(\rho)_{ij,kl}=\rho_{ik,jl}$, is a second, independent proof of entanglement.

The search is for a product vector $|a\rangle|b\rangle$ with complex amplitudes, $|a\rangle,|b\rangle\neq0$, orthogonal to all five; the smallest eigenvalue of $\rho^{T_B}$ is exact, and the figure plots the realignment value $\|R(\rho_p)\|_1-1$ against $p$, marking where it crosses zero:

```wl
tilesStates = With[{e = IdentityMatrix[3], u = {1, 1, 1}}, Normalize[Flatten[KroneckerProduct @@ #]] & /@ {{e[[1]], e[[1]] - e[[2]]}, {e[[3]], e[[2]] - e[[3]]}, {e[[1]] - e[[2]], e[[3]]}, {e[[2]] - e[[3]], e[[1]]}, {u, u}}];
tilesRho = (IdentityMatrix[9] - Total[KroneckerProduct[#, #] & /@ tilesStates])/4;
tilesNoisy[p_] := (1 - p) tilesRho + p IdentityMatrix[9]/9;
tilesProductQ = With[{a = Array[tilesA, 3], b = Array[tilesB, 3]},
   Reduce[And @@ Thread[tilesStates . Flatten[KroneckerProduct[a, b]] == 0] && Or @@ Thread[a != 0] && Or @@ Thread[b != 0], Join[a, b]]];
tilesRegroup[r_, perm_] := ArrayReshape[Transpose[ArrayReshape[r, {3, 3, 3, 3}], perm], {9, 9}];
tilesMinPT = Min[QuantumPartialTranspose[QuantumState[tilesRho, {3, 3}], {2}]["Eigenvalues"]];
tilesRealign[p_?NumericQ] := Total[SingularValueList[N[tilesRegroup[tilesNoisy[p], {1, 3, 2, 4}]], Tolerance -> 0]] - 1;
tilesPStar = \[FormalP] /. FindRoot[tilesRealign[\[FormalP]] == 0, {\[FormalP], 0.1}];
Labeled[Plot[tilesRealign[p], {p, 0, 0.2}, PlotRange -> All, Frame -> True, GridLines -> {{tilesPStar}, {0}}, GridLinesStyle -> Dashed,
   FrameLabel -> {"noise weight p", Row[{"\[LeftDoubleBracketingBar]", "R(", Subscript["\[Rho]", "p"], ")", "\[RightDoubleBracketingBar]", Subscript["", 1], " - 1"}]}, ImageSize -> 420],
  Style[Column[{Row[{"product vector orthogonal to all five Tiles states exists:  ", tilesProductQ}], Row[{"smallest eigenvalue of ", Superscript["\[Rho]", Subscript["T", "B"]], ":  ", tilesMinPT}]}, Alignment -> Center], 13, FontFamily -> "Helvetica"], Top]
```

No product vector is orthogonal to the Tiles states, so the range of $\rho$ contains no product state and $\rho$ is entangled, although its partial transpose has no negative eigenvalue. Realignment sees the entanglement from $p=0$ up to the marked noise weight and is silent beyond it.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Coherent information of amplitude damping, maximized over the input


A qubit holds $\rho_p=(1-p)|0\rangle\langle0|+p\,|1\rangle\langle1|$, with $|1\rangle$ excited, and a reference qubit $R$ purifies it. Only the qubit crosses an amplitude-damping channel, which drops $|1\rangle$ to $|0\rangle$ with probability $\gamma$. Which input population $p$ carries the most quantum information across?

The channel and the coherent information, in bits, are

$$
\mathcal N_\gamma(\rho)=\hat K_0\rho \hat K_0^\dagger+\hat K_1\rho \hat K_1^\dagger,\quad \hat K_0=\begin{pmatrix}1&0\\0&\sqrt{1-\gamma}\end{pmatrix},\quad \hat K_1=\begin{pmatrix}0&\sqrt\gamma\\0&0\end{pmatrix},\qquad I_c=S(B)-S(RB),
$$

with $B$ the channel output and $S$ the von Neumann entropy. The channel commutes with phase rotations about $z$, and for $\gamma\le1/2$ it is degradable, so $I_c$ is concave in the input and a diagonal input is optimal: maximizing over $p$ maximizes over every input.

$I_c$ for general $p$ and $\gamma$ is compared with the binary-entropy form $h_2((1-\gamma)p)-h_2(\gamma p)$; the best input at each $\gamma<1/2$ is a dot:

```wl
icAD = QuantityMagnitude @ Simplify[With[{out = QuantumChannel["AmplitudeDamping"[γ], {1}][QuantumState[DiagonalMatrix[{1 - p, p}]]["Purify"]]},
     QuantumPartialTrace[out, {2}]["VonNeumannEntropy"] - out["VonNeumannEntropy"]], 0 < p < 1 && 0 < γ < 1];
h2AD[x_] := -x Log2[x] - (1 - x) Log2[1 - x];
icClosedQ = Simplify[icAD - (h2AD[(1 - γ) p] - h2AD[γ p]), 0 < p < 1 && 0 < γ < 1] === 0;
γsAD = {1/10, 1/4, 1/3, 1/2, 3/4};
belowHalfAD = Reduce[(D[h2AD[(1 - γ) p] - h2AD[γ p], p] /. p -> 1/2) >= 0 && 0 < γ < 1/2, γ, Reals];
bestAD = Through[{ArgMax, MaxValue}[{icAD /. γ -> #, 0 <= p <= 1}, p]] & /@ Select[γsAD, # < 1/2 &];
Show[Plot[Evaluate[icAD /. γ -> # & /@ γsAD], {p, 0, 1}, PlotRange -> All, PlotLegends -> (Row[{"γ = ", #}] & /@ γsAD)],
  ListPlot[List /@ N[bestAD], PlotMarkers -> {Automatic, 10}],
  PlotRange -> All, Frame -> True, FrameLabel -> {"p", Row[{Subscript["I", "c"], " (bits)"}]},
  PlotLabel -> Column[{Row[{Subscript["I", "c"], " = ", Subscript["h", 2], "((1-γ)p) - ", Subscript["h", 2], "(γp):  ", icClosedQ}],
    Row[{"slope at p = 1/2 ≥ 0 for some γ < 1/2:  ", belowHalfAD}]}, Alignment -> Center], ImageSize -> 420]
```

The output keeps the population $(1-\gamma)p$ and the environment receives $\gamma p$. Below $\gamma=1/2$ the slope at $p=1/2$ is negative for every $\gamma$, so the best input has $p<1/2$, fewer excitations than the maximally mixed state. At $\gamma=1/2$, $I_c$ vanishes for every input, and at $\gamma=3/4$ no input gives a positive $I_c$: the environment then gets the larger share.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Data-processing inequality under named channels


Two qubit states $\rho$ and $\sigma$ have Bloch vectors of length $3/5$ and $4/5$ at an angle $\theta$, so they do not commute. Both pass through the same channel at noise strength $p$. Can any $p$ make one of these pairs easier to tell apart?

Their distinguishability is the relative entropy $S(\rho\|\sigma)=\operatorname{Tr}\rho\,(\log_2\rho-\log_2\sigma)$, and the data-processing inequality says $S(\mathcal N_p(\rho)\|\mathcal N_p(\sigma))\le S(\rho\|\sigma)$ for every channel $\mathcal N_p(\rho)=\sum_k \hat K_k\rho \hat K_k^\dagger$. The four channels have Kraus operators

$$
\text{depolarizing: }\sqrt{1-\frac{3p}{4}}\,I,\ \frac{\sqrt p}{2}X,\ \frac{\sqrt p}{2}Y,\ \frac{\sqrt p}{2}Z;\qquad \text{bit flip: }\sqrt{1-p}\,I,\ \sqrt p\,X;\qquad \text{phase flip: }\sqrt{1-p}\,I,\ \sqrt p\,Z;\qquad \text{phase damping: }\begin{pmatrix}1&0\\0&\sqrt{1-p}\end{pmatrix},\ \begin{pmatrix}0&0\\0&\sqrt p\end{pmatrix}.
$$

Bit flip and phase flip mix the identity with a unitary $\hat U$, so joint convexity gives $S(\mathcal N_p(\rho)\|\mathcal N_p(\sigma))\le(1-p)\,S(\rho\|\sigma)+p\,S(\hat U\rho \hat U^\dagger\|\hat U\sigma \hat U^\dagger)=S(\rho\|\sigma)$ for every $p$.

The processed relative entropy is a closed form in $p$, with $\rho$ along $z$ and $\sigma$ in the $xz$ plane, and all of $0<p<1$ is searched exactly for a depolarizing or phase-damping value above $S(\rho\|\sigma)$ (dashed):

```wl
relEnt[a_, b_] := Simplify[QuantityMagnitude @ QuantumDistance[a, b, "RelativeEntropy"], 0 < p < 1];
channelsDPI = {"Depolarizing", "BitFlip", "PhaseFlip", "PhaseDamping"};
panelDPI[θ_] := With[{pair = QuantumState["BlochVector"[#]] & /@ {{0, 0, 3/5}, 4/5 {Sin[θ], 0, Cos[θ]}}},
   With[{before = relEnt @@ pair, after = relEnt @@ (QuantumChannel[#[p]] /@ pair) & /@ channelsDPI},
    Plot[Evaluate[after], {p, 0, 1}, PlotRange -> {0, All}, PlotLegends -> channelsDPI, GridLines -> {None, {before}}, GridLinesStyle -> Directive[Thick, Dashed],
     Frame -> True, FrameLabel -> {"p", "S(\[ScriptCapitalN](ρ)||\[ScriptCapitalN](σ)) (bits)"}, ImageSize -> 300,
     PlotLabel -> Row[{"θ = ", θ, ",   Reduce, p raising S: ", Or @@ (Reduce[# > before && 0 < p < 1, p, Reals] & /@ after[[{1, 4}]])}]]]];
panelsDPI = panelDPI /@ {Pi/3, Pi/2, 2 Pi/3};
Legended[Row[First /@ panelsDPI], legendBelow[Last[First[panelsDPI]]]]
```

No noise strength raises the relative entropy, for any $p$, not only the points drawn. Depolarizing noise only lowers it, down to zero at $p=1$ where both states become the maximally mixed state; phase damping only lowers it too. Bit flip and phase flip become the unitaries $X$ and $Z$ at $p=1$, so their curves climb back exactly to the dashed line there.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Two qubits have no bound entanglement


A bound-entangled state is entangled, yet its partial transpose $\rho^{T_B}$ is still positive, so no local protocol can distill a Bell pair from it. Can a two-qubit state be bound entangled?

Write every two-qubit state as $\rho=\sum_{i,j=0}^{3}c_{ij}\,\hat\sigma_i\otimes\hat\sigma_j$ with 16 real symbols $c_{ij}$. Wootters' concurrence $C=\max(0,\lambda_1-\lambda_2-\lambda_3-\lambda_4)$ vanishes exactly on separable states, where $\lambda_1\ge\dots\ge\lambda_4$ are the square roots of the eigenvalues of $R=\rho\,(\hat\sigma_y\otimes\hat\sigma_y)\,\rho^*(\hat\sigma_y\otimes\hat\sigma_y)$. Since $\sum_k\lambda_k^2=\operatorname{Tr}R$, $\sum_k\lambda_k^4=\operatorname{Tr}R^2$ and $\prod_k\lambda_k=\det\rho$, the claim is

$$
16\det\rho^{T_B}=(\operatorname{Tr}R)^2-2\operatorname{Tr}R^2+8\det\rho=\prod_{k=1}^{4}\left(2\lambda_k-\sum_j\lambda_j\right),
$$

and a nonnegative product forces $\lambda_1\le\lambda_2+\lambda_3+\lambda_4$, that is $C=0$.

The first equality is tested for general $c_{ij}$, all 16 of them, and the second for general $\lambda_k$; the second statement is decided over every ordered quadruple $\lambda_1\ge\lambda_2\ge\lambda_3\ge\lambda_4\ge0$, and the figure draws, in units of $\lambda_1$, where $C>0$ and where the product is nonnegative:

```wl
paulis2q = Flatten[Outer[KroneckerProduct, #, #, 1] &[PauliMatrix[{0, 1, 2, 3}]], 1];
yyPPT = KroneckerProduct[PauliMatrix[2], PauliMatrix[2]];
ρ2q = Array[c, 16] . paulis2q;
r2q = ρ2q . yyPPT . ComplexExpand[Conjugate[ρ2q]] . yyPPT;
detTB = Det[Normal @ QuantumPartialTranspose[QuantumState[ρ2q, 2], {2}]["DensityMatrix"]];
λsPPT = Array[l, 4];
identityQ = Expand[16 detTB - (Tr[r2q]^2 - 2 Tr[r2q . r2q] + 8 Det[ρ2q])] === 0 &&
   Expand[Times @@ (2 λsPPT - Total[λsPPT]) - (Total[λsPPT^2]^2 - 2 Total[λsPPT^4] + 8 Times @@ λsPPT)] === 0;
noBoundQ = With[{λ = λsPPT}, Resolve[ForAll[λ, GreaterEqual @@ Append[λ, 0] && Times @@ (2 λ - Total[λ]) >= 0, First[λ] - Total[Rest[λ]] <= 0], Reals]];
panelPPT[l4_] := RegionPlot[{v <= u && 1 - u - v - l4 > 0, v <= u && Times @@ (2 {1, u, v, l4} - (1 + u + v + l4)) >= 0}, {u, l4, 1}, {v, l4, 1}, PlotPoints -> 60, PlotLegends -> {"C > 0", Row[{"16 det ", Superscript["ρ", Subscript["T", "B"]], " ≥ 0"}]},
   Frame -> True, FrameLabel -> {Row[{Subscript["λ", 2], "/", Subscript["λ", 1]}], Row[{Subscript["λ", 3], "/", Subscript["λ", 1]}]},
   PlotLabel -> Row[{Subscript["λ", 4], "/", Subscript["λ", 1], " = ", l4}], ImageSize -> 260];
panelsPPT = panelPPT /@ {0, 1/10, 1/5};
Labeled[Legended[PlotGrid[{First /@ panelsPPT}, ImageSize -> 780], legendBelow[Last[First[panelsPPT]]]], Style[Row[{"both identities: ", identityQ, ",   Resolve: ", noBoundQ}], 13, FontFamily -> "Helvetica"], Top]
```

The two regions never overlap and together fill the ordered triangle $\lambda_3\le\lambda_2$ of every panel: wherever the determinant of the partial transpose is nonnegative, the concurrence is zero, and this is proved for every quadruple, not just the slices drawn. A positive partial transpose has a nonnegative determinant, so it forces $C=0$, which is separability. A two-qubit state is either separable or has a negative partial transpose, and every two-qubit state with a negative partial transpose can be distilled (Horodecki 1997), so nothing is left to be bound.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The Helstrom bound, minimized over every measurement


One of two qubit states arrives, $\rho_1$ with prior $p$ or $\rho_2$ with prior $1-p$; their Bloch vectors have equal length $r$ and meet at angle $\theta$. One two-outcome measurement $\{\hat E_1,\hat E_2\}$ must name it, with $\hat E_2$ meaning "it was $\rho_2$". What is the smallest possible error?

$$
P_{\rm err}=p\operatorname{Tr}\rho_1\hat E_2+(1-p)\operatorname{Tr}\rho_2\hat E_1=(1-p)+\operatorname{Tr}\Lambda \hat E_2,\qquad \Lambda=p\rho_1-(1-p)\rho_2,\qquad \hat E_1=I-\hat E_2,\quad 0\le \hat E_2\le I .
$$

A unitary maps the set of effects to itself, so in the eigenbasis of $\Lambda=\mathrm{diag}(l_1,l_2)$ the most general effect is $\hat E_2=\begin{pmatrix}a&b+ic\\b-ic&e\end{pmatrix}$ with $\hat E_2\ge0$ and $I-\hat E_2\ge0$.

Every $k$ with $\operatorname{Tr}\Lambda \hat E_2\ge k$ for all such effects is found exactly, giving the lowest error $(1-p)+\max k$ (dashed, thinner for larger $\theta$), against the projective $\hat E_2=\frac12(I+\hat m\cdot\hat{\vec\sigma})$ with $\hat m$ at angle $\phi$ from $\rho_1$'s Bloch vector in their plane (solid):

```wl
helstromK = Reduce[ForAll[{ea, eb, ec, ee}, 0 <= ea <= 1 && 0 <= ee <= 1 && ea ee >= eb^2 + ec^2 && (1 - ea) (1 - ee) >= eb^2 + ec^2,
    Tr[DiagonalMatrix[{l1, l2}] . {{ea, eb + I ec}, {eb - I ec, ee}}] >= k], k, Reals];
blochRho[v_] := (IdentityMatrix[2] + v . PauliMatrix[{1, 2, 3}])/2;
θsHel = {Pi/4, Pi/2, Pi};
panelHel[r_, pr_] := With[{ρ1 = blochRho[{0, 0, r}], ρ2s = blochRho[r {Sin[#], 0, Cos[#]}] & /@ θsHel, e2 = blochRho[{Sin[φ], 0, Cos[φ]}]},
   Show[Plot[Evaluate[pr Tr[ρ1 . e2] + (1 - pr) Tr[# . (IdentityMatrix[2] - e2)] & /@ ρ2s], {φ, 0, 2 Pi}, PlotRange -> All, PlotLegends -> (Row[{"θ = ", #}] & /@ θsHel)],
    Plot[Evaluate[(1 - pr) + MaxValue[{k, helstromK /. Thread[{l1, l2} -> Eigenvalues[pr ρ1 - (1 - pr) #]]}, k] & /@ ρ2s], {φ, 0, 2 Pi}, PlotStyle -> (Directive[Dashed, AbsoluteThickness[#]] & /@ {4, 2.5, 1})],
    PlotRange -> {0, 1}, Axes -> False, Frame -> True, FrameLabel -> {"φ", Subscript["P", "err"]}, PlotLabel -> Row[{"r = ", r, ",   p = ", pr}], ImageSize -> 280]];
panelsHel = Outer[panelHel, {1, 1/2}, {1/2, 3/4, 9/10}];
Labeled[Legended[Grid[Map[First, panelsHel, {2}]], legendBelow[Last[panelsHel[[1, 1]]]]], Style[Row[{"Reduce: ", TraditionalForm[helstromK /. {l1 -> Subscript["l", 1], l2 -> Subscript["l", 2]}]}], 13, FontFamily -> "Helvetica"], Top]
```

The largest $k$ is the sum of the negative eigenvalues of $\Lambda$, so with $\operatorname{Tr}\Lambda=2p-1$ the lowest error is the Helstrom bound $\frac12(1-\|\Lambda\|_1)$. When $\Lambda$ has one negative and one nonnegative eigenvalue, the best rank-one projector touches the bound; when both are negative, which needs $p<1/2$, the optimum is $\hat E_2=I$ (always guess $\rho_2$). For mixed states with a strong prior $\Lambda$ has no negative eigenvalue, the optimum is $\hat E_2=0$ (always guess $\rho_1$), and a rank-one projector reaches it only when $\Lambda$ has a zero eigenvalue.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The Werner-state entanglement threshold in every dimension


Two qudits of dimension $d$ share a Werner state: weight $p$ on the symmetric subspace and $1-p$ on the antisymmetric one, each spread evenly. At which $p$ does its partial transpose turn negative in each dimension $d$?

With $\hat F$ the swap of the two qudits,

$$
\rho_W(p)=\frac{p}{d_s}\,\hat\Pi_s+\frac{1-p}{d_a}\,\hat\Pi_a,\qquad \hat\Pi_{s,a}=\frac{I\pm \hat F}{2},\qquad d_{s,a}=\frac{d(d\pm1)}{2},
$$

so for $d=2$, $\hat\Pi_a$ is the singlet projector. A negative eigenvalue of the partial transpose $\rho_W^{T_B}$ proves the state entangled. Writing $\rho_W=\alpha I+\beta \hat F$ and using $\hat F^{T_B}=d\,|\Phi^+\rangle\langle\Phi^+|$, the eigenvalues of $\rho_W^{T_B}$ are $\alpha$ and $\alpha+d\beta$ in every dimension.

The region where $\alpha+d\beta<0$ is found for general $d\ge2$ and checked against the partial-transpose eigenvalues at $d=2,3,4$; the eigenvalue that can turn negative is drawn, with the entangled range shaded:

```wl
ptEigenW = AssociationMap[Union @ Simplify @ QuantumPartialTranspose[QuantumState["Werner"[p, #]], {2}]["Eigenvalues"] &, {2, 3, 4}];
{entangledW, tiedW, lowFormW} = With[{α = p/(d (d + 1)) + (1 - p)/(d (d - 1)), β = p/(d (d + 1)) - (1 - p)/(d (d - 1))},
   {Reduce[α + d β < 0 && 0 <= p <= 1 && d >= 2, p, Reals],
    AllTrue[Keys[ptEigenW], Simplify[SortBy[{α, α + d β} /. d -> #, # /. p -> 0 &] - SortBy[ptEigenW[#], # /. p -> 0 &]] === {0, 0} &], Factor[α + d β]}];
ptLowW = First[SortBy[#, # /. p -> 0 &]] & /@ Values[ptEigenW];
Plot[Evaluate[ptLowW], {p, 0, 1}, PlotLegends -> (Row[{"d = ", #}] & /@ Keys[ptEigenW]), PlotRange -> {{0, 1}, {-0.55, 0.55}},
  Prolog -> {Opacity[0.12], ColorData[97][4], Rectangle[{0, -1}, {1/2, 1}]}, Epilog -> Text[Style["entangled", 14, ColorData[97][4]], {0.25, 0.35}],
  Frame -> True, FrameLabel -> {"p", Row[{"eigenvalue α + dβ of ", Superscript["ρ", Subscript["T", "B"]]}]},
  PlotLabel -> Column[{Row[{"α + dβ = ", TraditionalForm[lowFormW]}],
     If[tiedW, Row[{"negative for ", TraditionalForm[entangledW]}], "QF eigenvalues differ from α + dβ"]}, Alignment -> Center], ImageSize -> 420]
```

In every dimension the same eigenvalue crosses zero at the same weight, so the threshold does not depend on $d$. Below it the state is entangled. For Werner states a positive partial transpose means separable in every $d$ (Werner 1989), so $p<1/2$ is exactly the entangled range for all $d$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Proving a circuit's phase approximation correct over the reals


A circuit that applies $|x\rangle\mapsto e^{2\pi i/\sqrt x}|x\rangle$ needs $1/\sqrt x$ in fixed-point arithmetic. Writing $x=4^km$ with $1\le m<4$ reduces it to $1/\sqrt m$, evaluated on $3\cdot2^r$ intervals $m=a+t$, $0\le t\le h=2^{-r}$, as a quadratic $p(t)$ with coefficients in $2^{-w}\mathbb Z$. Is the relative error below a bound $\eta$, proposed on each interval from finitely many samples, for every $m$ in it?

On each interval the claim is a polynomial inequality in one real variable,

$$
\forall\,t\in{}[0,h]:\quad p(t)>0\ \wedge\ (1-\eta)^2\le(a+t)\,p(t)^2\le(1+\eta)^2\iff\left|\sqrt{a+t}\,p(t)-1\right|\le\eta ,
$$

for $0\le\eta<1$, decided exactly over the reals, with rational $a$, $h$, $\eta$ and coefficients and no rounding. The quadratic interpolates $1/\sqrt{a+t}$ at Chebyshev nodes before rounding, and $\eta$ is proposed from 201 samples with a 2% margin.

For each $r$ and $w$ the figure plots $|\sqrt m\,p-1|$ against $m$ over every interval (solid) with each interval's $\eta$ (dashed), and reports whether the claim is proved on all intervals:

```wl
intervalsCert[r_, w_] := With[{h = 2^-r, nodes = 2^-r (1 - Cos[Pi (2 Range[3] - 1)/6])/2, ts = Subdivide[0., 2.^-r, 200]}, Table[
     With[{poly = Round[2^w CoefficientList[InterpolatingPolynomial[Transpose[{nodes, N[1/Sqrt[a + nodes], 30]}], \[FormalT]], \[FormalT]]] . \[FormalT]^Range[0, 2]/2^w},
      With[{eta = Ceiling[1.02 2^40 Max[Abs[Sqrt[a + ts] (poly /. \[FormalT] -> ts) - 1]]]/2^40},
       <|"a" -> a, "h" -> h, "poly" -> poly, "eta" -> eta,
        "proof" -> Resolve[ForAll[\[FormalT], 0 <= \[FormalT] <= h, poly > 0 && (1 - eta)^2 <= (a + \[FormalT]) poly^2 <= (1 + eta)^2], Reals]|>]], {a, 1, 4 - h, h}]];
panelCert[r_, w_] := With[{ivs = intervalsCert[r, w]},
   Show[Plot[Evaluate[Piecewise[N[{Abs[Sqrt[m] (#poly /. \[FormalT] -> m - #a) - 1], #a <= m < #a + #h}] & /@ ivs]], {m, 1, 4}, PlotRange -> All, PlotPoints -> 100, Exclusions -> None],
    Plot[Evaluate[Piecewise[N[{#eta, #a <= m < #a + #h}] & /@ ivs]], {m, 1, 4}, PlotStyle -> Directive[ColorData[97][2], Dashed], PlotRange -> All, Exclusions -> None],
    PlotRange -> {0, All}, Axes -> False, Frame -> True, FrameLabel -> {"m", "|relative error|"}, ImageSize -> 280,
    PlotLabel -> Style[Column[{Row[{"r = ", r, ",   w = ", w, " bits"}], If[And @@ Lookup[ivs, "proof"], Row[{"bound proved on all ", Length[ivs], " intervals"}], "bound not proved on every interval"]}, Alignment -> Center], 13, FontFamily -> "Helvetica"]]];
Legended[Grid[Outer[panelCert, {1, 2}, {8, 16, 24}]], legendBelow[LineLegend[{ColorData[97][1], Directive[ColorData[97][2], Dashed]},
  {Row[{"|", Sqrt["m"], " p \[Minus] 1|"}], Row[{"proved bound \[Eta] on each interval"}]}]]]
```

The bound is proved on every interval in every panel, so the error stays under the dashed line for every real $m$, not only at the samples that proposed $\eta$. With 8-bit coefficients rounding sets the error, which jumps from interval to interval and barely shrinks when the intervals halve. By 24 bits rounding is negligible and the error is the interpolation's ripple, which shrinks as the intervals halve and as $m$ grows.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Tensors, symmetry and tensor networks

The TensorNetworks paclet contracts networks and builds matrix-product states, and its Symmetry tools turn Young diagrams into projectors and dimensions; built-in `TensorReduce` proves tensor identities at symbolic dimension.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### One site of an arbitrarily long Ising chain after a quench


A chain of spins with coupling $J$ in a transverse field $h$ starts with every spin up, $|0\cdots0\rangle$, and evolves under the transverse-field Ising Hamiltonian, Trotterized with step $\delta t$. What is the magnetization $\langle Z_c\rangle$ of a site deep inside the chain after $d$ steps?

With $\hbar=1$, $\hat H=J\sum_iZ_iZ_{i+1}+h\sum_iX_i$, and one step is $\prod_iR_X^{(i)}(2h\,\delta t)\prod_iR_{ZZ}^{(i,i+1)}(2J\,\delta t)$ with $R_{ZZ}(\vartheta)=e^{-i\vartheta Z\otimes Z/2}$. After $d$ steps $\langle Z_c\rangle=\langle0|U^\dagger Z_cU|0\rangle$, and the first $ZZ$ layer only multiplies $|0\cdots0\rangle$ by a phase, so every gate outside the sites $c-d+1,\dots,c+d-1$ cancels against its conjugate and a window of $2d-1$ sites gives the value for a chain of any length.

The circuit runs from $|0\cdots0\rangle$ on a window of $m$ sites centered on $c$, giving $\langle Z_c\rangle$ for the $2d-1$-site window (filled), a larger $2d+3$-site window (open) and a smaller $2d-3$-site window (crosses), with $J\delta t=0.2$:

```wl
lcJdt = 0.2;
lcQuench[m_, d_, hJ_] := QuantumCircuitOperator @ Join[{StringRepeat["0", m]},
    Flatten @ ConstantArray[Join[Thread["R"[2 lcJdt, "ZZ"] -> Partition[Range[m], 2, 1]], Thread["RX"[2 hJ lcJdt] -> Range[m]]], d]];
lcZ[m_, d_, hJ_] := With[{st = lcQuench[m, d, hJ][]}, Re[st["Dagger"][QuantumOperator["Z", {(m + 1)/2}][st]]["Scalar"]]];
lcPanel[hJ_] := ListPlot[Table[{lcJdt d, lcZ[2 d + #, d, hJ]}, {d, 2, 5}] & /@ {-1, 3, -3}, PlotRange -> {{0, 1.1}, {-1, 1.05}},
    PlotMarkers -> {{"\[FilledCircle]", 9}, {"\[EmptyCircle]", 18}, {"\[Times]", 14}}, PlotLegends -> {Row[{2, "d - 1 sites"}], Row[{2, "d + 3 sites"}], Row[{2, "d - 3 sites"}]},
    Frame -> True, FrameLabel -> {"J t", Row[{"\[LeftAngleBracket]", Subscript["Z", "c"], "\[RightAngleBracket]"}]}, PlotLabel -> Row[{"h/J = ", hJ}], ImageSize -> 300];
lcPanels = lcPanel /@ {0.5, 1, 2};
Legended[PlotGrid[{First /@ lcPanels}, ImageSize -> 900], legendBelow[Last[First[lcPanels]]]]
```

At every depth the $2d-1$ and $2d+3$ windows give the same magnetization, so each filled point is exact for a bulk site of an arbitrarily long chain, while the $2d-3$ window cuts into the light cone and departs from them, visibly at $Jt=0.4$ for $h/J$ of 1 and above. The magnetization falls faster as $h/J$ grows, and at $h/J=2$ it changes sign within the times shown.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Tensor identities proved at symbolic dimension


Tensors in dimension $d$ are given only by their index symmetries, never by entries: an antisymmetric $A$, a symmetric $S$, an inverse metric $g$, and a tensor $R$ with the pair symmetries of the Riemann tensor; in three dimensions, the Levi-Civita symbol $\varepsilon$ and a Hessian $H_{jk}=\partial_j\partial_k f$. Which of these contractions vanish by symmetry alone?

$$
A_{ij}=-A_{ji},\quad S_{ij}=S_{ji},\quad g^{ij}=g^{ji},\quad R_{abcd}=-R_{bacd}=-R_{abdc}=R_{cdab},\quad \varepsilon_{ijk}\ \text{totally antisymmetric},\quad H_{jk}=H_{kj},
$$

and the Ricci tensor is $R_{bd}=R_{abcd}\,g^{ac}$ (repeated indices summed).

Each contraction reduces to a canonical form from these declarations alone:

```wl
tensorDeclarations = {aT ∈ Arrays[{d, d}, Reals, Antisymmetric[{1, 2}]], sT ∈ Arrays[{d, d}, Reals, Symmetric[{1, 2}]],
   gT ∈ Arrays[{d, d}, Reals, Symmetric[{1, 2}]], rT ∈ Arrays[{d, d, d, d}, Reals, {{Cycles[{{1, 2}}], -1}, {Cycles[{{3, 4}}], -1}, {Cycles[{{1, 3}, {2, 4}}], 1}}],
   epsT ∈ Arrays[{3, 3, 3}, Reals, Antisymmetric[{1, 2, 3}]], hT ∈ Arrays[{3, 3}, Reals, Symmetric[{1, 2}]]};
ricciT = TensorContract[TensorProduct[rT, gT], {{1, 5}, {3, 6}}];
reduceT = TensorReduce[#, Assumptions -> tensorDeclarations] &;
reducesToZero = MatchQ[reduceT[#], 0 | _SymbolicZerosArray] &;
Style[Grid[{{"contraction", "reduces to zero by symmetry"},
    {Row[{Subscript["A", "ij"], Subscript["S", "ij"]}], reducesToZero[TensorContract[TensorProduct[aT, sT], {{1, 3}, {2, 4}}]]},
    {Row[{Subscript["R", "bd"], " \[Minus] ", Subscript["R", "db"]}], reduceT[ricciT] === reduceT[TensorTranspose[ricciT, {2, 1}]]},
    {Row[{Subscript["\[CurlyEpsilon]", "ijk"], Subscript["H", "jk"]}], reducesToZero[TensorContract[TensorProduct[epsT, hT], {{2, 4}, {3, 5}}]]},
    {Row[{Subscript["A", "ij"], Subscript["A", "ij"]}], reducesToZero[TensorContract[TensorProduct[aT, aT], {{1, 3}, {2, 4}}]]}},
   Frame -> All, Alignment -> Left, Spacings -> {2, 1}], 14, FontFamily -> "Helvetica"]
```

Symmetry alone makes $A_{ij}S_{ij}$ vanish for every $d$, makes the Ricci tensor symmetric, and makes the curl of a gradient vanish in three dimensions. The sum of squares $A_{ij}A_{ij}$, which no symmetry forces to vanish, is correctly not reduced to zero.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Determinant and trace rules proved at symbolic dimension


$A$ and $B$ are arbitrary complex $n\times n$ matrices, and $K$ is an arbitrary $m\times m$ matrix acting on a second system of $m$ levels; no entry and no value of $n$ or $m$ is given. Which of the textbook determinant and trace rules can be proved for every $n$ and $m$ at once?

The candidates are

$$
\det(AB)=\det A\,\det B,\quad \det e^{A}=e^{\operatorname{Tr}A},\quad \operatorname{Tr}(BAB^{-1})=\operatorname{Tr}A,\quad (AB)^{-1}=B^{-1}A^{-1},\quad \det(A\otimes K)=(\det A)^m(\det K)^n,\quad \operatorname{Tr}(A\otimes K)=\operatorname{Tr}A\,\operatorname{Tr}K,
$$

and two rules that are false in general, $\det(A+B)=\det A+\det B$ and $\operatorname{Tr}e^{A}=e^{\operatorname{Tr}A}$.

Each column tries to reduce the difference of the two sides to zero, with index symmetries only, then with the matrix-algebra identities by expansion and by simplification, all at general $n$ and $m$; the last column substitutes random complex matrices for every $n=2,3,4$ and $m=2,3$:

```wl
matIdDims = matIdA ∈ Matrices[{n, n}] && matIdB ∈ Matrices[{n, n}] && matIdK ∈ Matrices[{m, m}];
matIdRules = {{"det(AB) = det A det B", Det[matIdA . matIdB] - Det[matIdA] Det[matIdB]},
   {"det exp(A) = exp(Tr A)", Det[MatrixExp[matIdA]] - Exp[Tr[matIdA]]},
   {Row[{"Tr(B A ", Superscript["B", -1], ") = Tr A"}], Tr[matIdB . matIdA . Inverse[matIdB]] - Tr[matIdA]},
   {Row[{Superscript["(AB)", -1], " = ", Superscript["B", -1], Superscript["A", -1]}], Inverse[matIdA . matIdB] - Inverse[matIdB] . Inverse[matIdA]},
   {Row[{"det(A ⊗ K) = ", Superscript["(det A)", "m"], Superscript["(det K)", "n"]}], Det[KroneckerProduct[matIdA, matIdK]] - Det[matIdA]^m Det[matIdK]^n},
   {"Tr(A ⊗ K) = Tr A Tr K", Tr[KroneckerProduct[matIdA, matIdK]] - Tr[matIdA] Tr[matIdK]},
   {"det(A + B) = det A + det B", Det[matIdA + matIdB] - Det[matIdA] - Det[matIdB]},
   {"Tr exp(A) = exp(Tr A)", Tr[MatrixExp[matIdA]] - Exp[Tr[matIdA]]}};
matIdZeroQ = If[MatchQ[#, 0 | _SymbolicZerosArray], "0", "not 0"] &;
matIdNumericQ[expr_] := BlockRandom[SeedRandom[3]; AllTrue[Tuples[{{2, 3, 4}, {2, 3}}],
     Max[Abs[expr /. {n -> #[[1]], m -> #[[2]], matIdA -> RandomComplex[{-1 - I, 1 + I}, {#[[1]], #[[1]]}], matIdB -> RandomComplex[{-1 - I, 1 + I}, {#[[1]], #[[1]]}],
          matIdK -> RandomComplex[{-1 - I, 1 + I}, {#[[2]], #[[2]]}]}]] < 10^-9 &]];
Grid[Prepend[{#[[1]], matIdZeroQ[TensorReduce[#[[2]], Assumptions -> matIdDims]], matIdZeroQ[ArrayExpand[#[[2]], matIdDims]],
      matIdZeroQ[ArraySimplify[#[[2]], matIdDims]], If[matIdNumericQ[#[[2]]], "holds", "fails"]} & /@ matIdRules,
   {"rule", "TensorReduce", "ArrayExpand", "ArraySimplify", "random matrices"}],
  Frame -> All, Alignment -> {{Left, Center, Center, Center, Center}}, Spacings -> {2, 1}, BaseStyle -> {14, FontFamily -> "Helvetica"}]
```

Index symmetries alone settle only the trace of a Kronecker product, which is a pure index contraction. The matrix-algebra identities prove all six true rules for every $n$ and $m$, and leave the two false rules standing, as the random matrices confirm.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The first Bianchi identity from the Riemann tensor's Young diagram


The Riemann tensor $R_{abcd}$ in $d$ dimensions has the pair symmetries $R_{abcd}=-R_{bacd}=-R_{abdc}=R_{cdab}$ and the first Bianchi identity $R_{abcd}+R_{acdb}+R_{adbc}=0$. Suppose only that $R$ lies in the image of the Young projector of shape $\{2,2\}$, which is one irreducible representation of $GL(d)$ inside the $\{2,2\}$ sector of $S_4$ acting on the four slots. Does that projector alone impose the first Bianchi identity?

For the tableau with rows $\{1,3\},\{2,4\}$ the projector is

$$
P_{\{2,2\}}\,T=\frac{2}{4!}\;\mathcal A_{\{1,2\}}\mathcal A_{\{3,4\}}\,\mathcal S_{\{1,3\}}\mathcal S_{\{2,4\}}\,T,
$$

where $\mathcal S$ sums over the permutations of a row of slots, $\mathcal A$ sums over the permutations of a column with their signs, and $2$ is the dimension of the $S_4$ representation $\{2,2\}$. The pair symmetries alone make $R$ a symmetric matrix on antisymmetric index pairs, with $N(N+1)/2$ components for $N=d(d-1)/2$; the image of $P_{\{2,2\}}$ has the dimension $s_{\{2,2\}}(d)$ of the $GL(d)$ representation.

The projector acts on a random integer tensor in $d=4$, and the slot-permuted copies of the result are searched for three-term sums that vanish, with slot 1 fixed; the components are then counted for $d=2,\dots,5$ and the two counts compared at general $d$:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`Symmetry`"];
riemannTableau = YoungTableau[{{1, 3}, {2, 4}}];
projectedT = BlockRandom[SeedRandom[7]; YoungProject[RandomInteger[{-9, 9}, {4, 4, 4, 4}], riemannTableau]];
threeTermRelations = Select[Subsets[Select[Permutations[Range[4]], First[#] == 1 &], {3}],
   Total[Transpose[projectedT, #] & /@ #] == ConstantArray[0, {4, 4, 4, 4}] &];
slotTerm[σ_] := Subscript["R", StringJoin[Characters["abcd"][[σ]]]];
pairCount[dd_] := Length[SymmetrizedIndependentComponents[{dd, dd, dd, dd}, {{Cycles[{{1, 2}}], -1}, {Cycles[{{3, 4}}], -1}, {Cycles[{{1, 3}, {2, 4}}], 1}}]];
projectorRank[dd_] := MatrixRank[Flatten[YoungProject[ArrayReshape[#, {dd, dd, dd, dd}], riemannTableau]] & /@ IdentityMatrix[dd^4]];
componentCounts = Table[{dd, pairCount[dd], projectorRank[dd]}, {dd, 2, 5}];
ListPlot[{componentCounts[[All, {1, 2}]], componentCounts[[All, {1, 3}]]}, Joined -> True, PlotMarkers -> Automatic,
  PlotLegends -> {"pair symmetries only", Row[{"rank of ", Subscript["P", "{2,2}"]}]}, Frame -> True,
  FrameLabel -> {"dimension d", "independent components"}, FrameTicks -> {{Automatic, None}, {Range[2, 5], None}}, ImageSize -> 440,
  PlotLabel -> Style[Column[{Row[{"vanishing three-term sums:  ", Column[Row[{Row[slotTerm /@ #, " + "], " = 0"}] & /@ threeTermRelations]}],
     Row[{"rank = ", Subscript["s", "{2,2}"], "(d) = ", TraditionalForm[Factor[SchurDimension[{2, 2}, d]]], " at d = 2, ..., 5:  ",
       componentCounts[[All, 3]] === (SchurDimension[{2, 2}, #] & /@ Range[2, 5])}],
     Row[{"pair count \[Minus] ", Subscript["s", "{2,2}"], "(d) = ", TraditionalForm[Factor[With[{np = d (d - 1)/2}, np (np + 1)/2] - SchurDimension[{2, 2}, d]]]}]}],
    13, FontFamily -> "Helvetica"]]
```

The search returns only the cyclic sum over slots 2, 3, 4 (once with each sign), which is the first Bianchi identity; nothing beyond the projector was supplied. The rank of the projector is the textbook count of Riemann components, and the gap to the pair-symmetry count vanishes for $d\le3$: only in four or more dimensions does Bianchi remove components.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Estimating a qubit's spectrum by Schur sampling


Take $n$ copies of a qubit in an unknown state $\rho=\frac12(\mathbb 1+\vec b\cdot\hat{\vec\sigma})$ and measure only which Young diagram $\lambda$, that is which irreducible representation of the group $S_n$ permuting the copies, the joint state occupies. What does the outcome reveal about $\rho$?

The outcome probabilities are

$$
\Pr(\lambda)=\operatorname{Tr}\!\left[\hat P_\lambda\,\rho^{\otimes n}\right]=f^\lambda\,s_\lambda(x,y),\qquad s_{\{\lambda_1,\lambda_2\}}(x,y)=(xy)^{\lambda_2}\sum_{k=0}^{\lambda_1-\lambda_2}x^ky^{\lambda_1-\lambda_2-k},\qquad x,y=\frac12(1\pm|\vec b|),
$$

where $\hat P_\lambda$ is the orthogonal projector onto the $\lambda$ block of $(\mathbb C^2)^{\otimes n}$, $f^\lambda$ counts the standard tableaux of shape $\lambda$, and $s_\lambda$ is the Schur polynomial at the eigenvalues of $\rho$ (Schur-Weyl duality); for qubits only diagrams with at most two rows occur. The projector itself is $\hat P_\lambda=\hat M\hat M^{+}$ with $\hat M=\sum_T\hat P_T$, the sum of the Young projectors of the standard tableaux $T$ and $\hat M^+$ its pseudoinverse.

The closed form is checked against $\hat P_\lambda$ built from the Young projectors at $n=4$ for general $\vec b$, and the figure plots $\Pr(\lambda)$ against $(\lambda_1-\lambda_2)/n$:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`Symmetry`"];
schurPr[{l1_, l2_ : 0}, b_] := With[{x = (1 + b)/2, y = (1 - b)/2},
   TableauDimension[DeleteCases[{l1, l2}, 0]] (x y)^l2 Total[x^Range[0, l1 - l2] y^Range[l1 - l2, 0, -1]]];
schurProjector[n_, λ_] := With[{tableaux = Select[YoungTableau[TakeList[#, λ]] & /@ Permutations[Range[n]], StandardTableauQ]},
   With[{m = Transpose[Flatten[Total[Function[t, YoungProject[ArrayReshape[#, ConstantArray[2, n]], t]] /@ tableaux]] & /@ IdentityMatrix[2^n]]},
    m . PseudoInverse[m]]];
schurRho = (IdentityMatrix[2] + {bx, by, bz} . PauliMatrix[{1, 2, 3}])/2;
schurCheck = AllTrue[IntegerPartitions[4, 2], Expand[Tr[schurProjector[4, #] . KroneckerProduct @@ ConstantArray[schurRho, 4]] -
      (Expand[schurPr[#, bLength]] /. bLength -> Sqrt[bx^2 + by^2 + bz^2])] === 0 &];
schurPanel[n_] := ListPlot[Table[{(n - 2 l2)/n, schurPr[{n - l2, l2}, b]}, {b, {0, 1/2, 0.9}}, {l2, 0, Floor[n/2]}], Joined -> True,
   PlotMarkers -> Automatic, PlotRange -> {{-0.03, 1.03}, All}, GridLines -> {{0, 1/2, 0.9}, None}, GridLinesStyle -> Dashed,
   PlotLegends -> (Row[{"|b| = ", #}] & /@ {0, 1/2, 0.9}), Frame -> True, FrameLabel -> {Row[{"(", Subscript["λ", 1], " \[Minus] ", Subscript["λ", 2], ")/n"}], "Pr(λ)"},
   PlotLabel -> Row[{"n = ", n}], ImageSize -> 300];
schurPanels = schurPanel /@ {6, 20, 60};
Labeled[Legended[Row[First /@ schurPanels], legendBelow[Last[First[schurPanels]]]],
  Style[Row[{"n = 4:  Tr[", Subscript["P", "λ"], Superscript["ρ", "⊗4"], "] = ", Superscript["f", "λ"], Subscript["s", "λ"], "(x, y) for every ", OverVector["b"], ":  ", schurCheck}], 13, FontFamily -> "Helvetica"], Top]
```

Every probability depends on $\vec b$ only through $|\vec b|$, so the diagram reveals the eigenvalues $(1\pm|\vec b|)/2$ and nothing about the eigenbasis. As $n$ grows the distribution of $(\lambda_1-\lambda_2)/n$ narrows around $|\vec b|$ (dashed), but it stays biased toward longer first rows, most visibly for the maximally mixed state.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Energies of all-to-all exchange from Young diagrams


$n$ identical sites, each with $d$ levels, are coupled pairwise by the exchange $\hat P_{ij}$, which swaps the states of sites $i$ and $j$, with one coupling $J$ for every pair; for $d=2$, $\hat P_{ij}=\frac12(1+\hat{\vec\sigma}_i\cdot\hat{\vec\sigma}_j)$ is the Heisenberg exchange. What are the energies and their degeneracies for every $n$ and $d$?

$$
\hat H=J\sum_{i<j}\hat P_{ij}\quad\text{on}\quad(\mathbb C^d)^{\otimes n}=\bigoplus_{\lambda\vdash n}V_\lambda\otimes W_\lambda(d),
$$

where the sum runs over Young diagrams $\lambda$ with $n$ boxes (Schur-Weyl duality): $V_\lambda$ carries the permutations of the sites, with $\dim V_\lambda=f^\lambda$ standard tableaux, and $W_\lambda(d)$ the common unitaries $U^{\otimes n}$, with $\dim W_\lambda(d)=s_\lambda(d)$, a polynomial in $d$ that vanishes when $\lambda$ has more than $d$ rows. The sum of all transpositions acts on block $\lambda$ as its content sum, so

$$
E_\lambda=J\,c(\lambda),\qquad c(\lambda)=\sum_{(r,k)\in\lambda}(k-r),\qquad \text{block dimension } f^\lambda s_\lambda(d),
$$

with $(r,k)$ the row and column of each box. The degeneracy of a level $E=Jc$ is the sum of $f^\lambda s_\lambda(d)$ over all diagrams with content sum $c$; up to $n=5$ every diagram has its own content sum, and the first coincidences come at $n=6$: $c(4,1,1)=c(3,3)=3$ and, for the conjugate diagrams, $c(3,1,1,1)=c(2,2,2)=-3$.

The figure plots each predicted level, with degeneracy $f^\lambda s_\lambda(d)$, against the eigenvalues of the $d^n\times d^n$ matrix of $\hat H$ built from qudit SWAPs, energies in units of $J$, with the $n=4$ levels for general $d$:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`Symmetry`"];
permSumContent[λ_] := Total[Flatten[Table[k - r, {r, Length[λ]}, {k, λ[[r]]}]]];
permSumLevels[nn_, dd_] := Select[{permSumContent[#], TableauDimension[#] SchurDimension[#, dd]} & /@ IntegerPartitions[nn], Last[#] > 0 &];
permSumEigen[nn_, dd_] := permSumEigen[nn, dd] = Eigenvalues[N[Normal[Total[QuantumOperator["SWAP"[dd], #] & /@ Subsets[Range[nn], {2}]]["Matrix"]]]];
permSumNumeric[nn_, dd_] := Tally[Round[permSumEigen[nn, dd]]];
permSumPanel[dd_, nn_] := ListPlot[{permSumLevels[nn, dd], permSumNumeric[nn, dd]}, Filling -> {1 -> Axis}, PlotMarkers -> {{"\[FilledCircle]", 9}, {"\[EmptyCircle]", 20}},
   PlotLegends -> {Row[{"J c(λ), ", Superscript["f", "λ"], Subscript["s", "λ"], "(d)"}], "eigenvalues of H"}, PlotRange -> {{-nn (nn - 1)/2 - 1, nn (nn - 1)/2 + 1}, {-0.06, 1.1} Max[permSumLevels[nn, dd][[All, 2]]]},
   Frame -> True, FrameLabel -> {"E / J", "degeneracy"}, PlotLabel -> Row[{"d = ", dd, ",   n = ", nn}], ImageSize -> 250];
permSumPanels = AssociationThread[{2, 3, 4}, Outer[permSumPanel, {2, 3, 4}, {3, 4, 5}]];
permSumLabel = Style[Column[Append[Row[{"λ = ", #, ":   E = ", permSumContent[#], " J,   degeneracy ", TraditionalForm[Factor[TableauDimension[#] SchurDimension[#, d]]]}] & /@ IntegerPartitions[4],
     Row[{"every eigenvalue within ", Superscript[10, -8], " J of an integer: ", Max[Abs[# - Round[#]] &[Flatten[Outer[permSumEigen, {3, 4, 5}, {2, 3, 4}]]]] < 10^-8}]]], 13, FontFamily -> "Helvetica"];
Manipulate[Labeled[Legended[Grid[{First /@ permSumPanels[dd]}], legendBelow[Last[permSumPanels[2][[1]]]]], permSumLabel, Top], Row[{Control[{{dd, 3, "d"}, Keys[permSumPanels], SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"d = ", dd}], <|2 -> "Two-row diagrams only: no antisymmetric level for n ≥ 3.",
        3 -> "n = 3 gains its antisymmetric level; n = 4, 5 still lack it.",
        4 -> "Only n = 5 still lacks its antisymmetric level."|>[dd]}}]]}], SaveDefinitions -> True]
```

Every eigenvalue of $\hat H$ sits on a predicted level with the predicted degeneracy, in every panel. A level exists only when its diagram fits in $d$ rows, so the most negative energies appear one at a time as $d$ grows, the fully antisymmetric level $-n(n-1)/2$ only once $d\ge n$. The degeneracy of a diagram with $r$ rows carries the factors $(d-1)\cdots(d-r+1)$, so it vanishes for every $d<r$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Exchange symmetry of two spins from Young projectors

Two particles of spin $j$, each with $d=2j+1$ states, share the space $\mathbb C^d\otimes\mathbb C^d$, which splits into a part symmetric under exchanging the particles and a part antisymmetric under it; for $j=\frac12$ these are the triplet and the singlet. Which total spins $S$ lie in each part, for every $j$?

With $\hat P_{12}$ the operator that swaps the two particles and $\hat{\vec s}$ the spin-$j$ matrices of one particle, the two parts are the images of the Young projectors of the one-row and the one-column diagram with two boxes, and the total spin is the sum of the two spins:

$$
\hat P_{\{2\}}=\frac12\left(\mathbb 1+\hat P_{12}\right),\qquad \hat P_{\{1,1\}}=\frac12\left(\mathbb 1-\hat P_{12}\right),\qquad \hat{\vec S}=\hat{\vec s}\otimes\mathbb 1+\mathbb 1\otimes\hat{\vec s},\qquad \hat{\vec S}^{\,2}=S(S+1),\quad S=0,1,\dots,2j .
$$

The figure gives every total spin $S$ found among the eigenvalues of $\hat{\vec S}^{\,2}$ on the image of each projector, for $j=\frac12,1,\dots,3$; filled circles mark the symmetric part and open squares the antisymmetric part:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`Symmetry`"];
twoSpinOps[j_] := With[{m = Range[j, -j, -1]}, With[{sp = SparseArray[{i_, k_} /; i == k - 1 :> Sqrt[j (j + 1) - m[[k]] (m[[k]] + 1)], {2 j + 1, 2 j + 1}]},
   {(sp + Transpose[sp])/2, (sp - Transpose[sp])/(2 I), DiagonalMatrix[m]}]];
twoSpinS2[j_] := Total[MatrixPower[KroneckerProduct[#, IdentityMatrix[2 j + 1]] + KroneckerProduct[IdentityMatrix[2 j + 1], #], 2] & /@ twoSpinOps[j]];
twoSpinBlock[j_, shape_] := With[{d = 2 j + 1}, With[{basis = Orthogonalize[DeleteCases[RowReduce[
        Flatten[YoungProject[ArrayReshape[#, {d, d}], YoungTableau[shape]]] & /@ IdentityMatrix[d^2]], {0 ..}]]},
   Union[(Sqrt[1 + 4 #] - 1)/2 & /@ Eigenvalues[Conjugate[basis] . twoSpinS2[j] . Transpose[basis]]]]];
twoSpinJs = Range[1/2, 3, 1/2];
ListPlot[Table[Flatten[Table[{j, s}, {j, twoSpinJs}, {s, twoSpinBlock[j, shape]}], 1], {shape, {{2}, {1, 1}}}],
  PlotMarkers -> {{"\[FilledCircle]", 14}, {"\[EmptySquare]", 16}}, PlotRange -> {{0.25, 3.25}, {-0.3, 6.3}},
  PlotLegends -> {Row[{"image of ", Subscript["P", "{2}"], ", symmetric"}], Row[{"image of ", Subscript["P", "{1,1}"], ", antisymmetric"}]},
  GridLines -> {None, Range[0, 6]}, Frame -> True, FrameLabel -> {"spin j of each particle", "total spin S"},
  FrameTicks -> {{Range[0, 6], None}, {{#, TraditionalForm[#]} & /@ twoSpinJs, None}}, ImageSize -> 440]
```

The two parts alternate from the top: the symmetric part holds $S=2j,2j-2,\dots$ and the antisymmetric part $S=2j-1,2j-3,\dots$, so a total spin is symmetric under exchange exactly when $2j-S$ is even. For $j=\frac12$ this is the triplet $S=1$ and the singlet $S=0$. Two identical spin-1 bosons in an exchange-symmetric orbital state, as in an s-wave collision, can therefore meet only with total spin 0 or 2.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Splitting n qudits into Young-diagram blocks

The space of $n$ sites with $d$ levels each carries two commuting actions, the permutations of the sites and the same unitary $U$ on every site, and it splits into blocks labelled by Young diagrams $\lambda$ with $n$ boxes. A standard tableau $T$ of shape $\lambda$ numbers the boxes $1,\dots,n$ increasing along every row and down every column, and its Young projector $\hat P_T$ picks out one copy of the unitary group's representation $\lambda$. Do the projections $\hat P_T\psi$ of a state $\psi$ over all standard tableaux add back to $\psi$, for every $n$?

$$
(\mathbb C^d)^{\otimes n}=\bigoplus_{\lambda\vdash n}V_\lambda\otimes W_\lambda(d),\qquad \dim\left(V_\lambda\otimes W_\lambda(d)\right)=f^\lambda s_\lambda(d),\qquad \hat P_T=\frac{f^\lambda}{n!}\,\mathcal A_T\,\mathcal S_T,\qquad \sum_T\hat P_T\overset{?}{=}\mathbb 1,
$$

where $f^\lambda=\dim V_\lambda$ is the number of standard tableaux of shape $\lambda$, $s_\lambda(d)=\dim W_\lambda(d)$ is a polynomial in $d$ that vanishes when $\lambda$ has more than $d$ rows, $\mathcal S_T$ sums over the permutations within each row of $T$, and $\mathcal A_T$ sums over the permutations within each column with their signs. Each $\hat P_T$ is idempotent, $\hat P_T^2=\hat P_T$.

The label writes out the block sizes $f^\lambda\times s_\lambda(d)$ for $n=3$ at general $d$, and compares the excess of the sum at $n=5$ with products of projectors of two tableaux of the same shape, on an integer tensor with $d=5$, where every diagram with five boxes fits; the figure gives the part of a random real state $\psi$ that the sum misses, for $n=2,\dots,7$:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`Symmetry`"];
youngStandard[n_] := Select[Flatten[Function[λ, YoungTableau[TakeList[#, λ]] & /@ Permutations[Range[n]]] /@ IntegerPartitions[n]], StandardTableauQ];
youngMiss[n_, dd_] := With[{ψ = BlockRandom[SeedRandom[n + 10 dd]; RandomReal[{-1, 1}, ConstantArray[dd, n]]]},
   Norm[Flatten[Total[YoungProject[ψ, #] & /@ youngStandard[n]] - ψ]]/Norm[Flatten[ψ]]];
youngPairs = Map[YoungTableau, {{{{1, 3, 5}, {2, 4}}, {{1, 2, 3}, {4, 5}}}, {{{1, 4}, {2, 5}, {3}}, {{1, 2}, {3, 4}, {5}}}}, {2}];
youngExcess = With[{ψ = BlockRandom[SeedRandom[5]; RandomInteger[{-9, 9}, ConstantArray[5, 5]]]},
   Norm[Flatten[Total[YoungProject[ψ, #] & /@ youngStandard[5]] - ψ - Total[YoungProject[YoungProject[ψ, #2], #1] & @@@ youngPairs]]]];
youngShow[t_] := Grid[TableauRows[t], Frame -> All, ItemSize -> 1.2];
youngP[i_] := Subscript["P", Subscript["T", i]];
Labeled[ListPlot[Table[{n, youngMiss[n, dd]}, {dd, {2, 3}}, {n, 2, 7}], Joined -> True, PlotMarkers -> Automatic, PlotRange -> {{1.7, 7.3}, {-0.03, 1}},
   PlotLegends -> {"qubits, d = 2", "qutrits, d = 3"}, Frame -> True, FrameTicks -> {{Automatic, None}, {Range[2, 7], None}},
   FrameLabel -> {"number of sites n", Row[{"\[LeftDoubleBracketingBar]", Subscript["\[Sum]", "T"], Subscript["P", "T"], "\[Psi] \[Minus] \[Psi]\[RightDoubleBracketingBar] / \[LeftDoubleBracketingBar]\[Psi]\[RightDoubleBracketingBar]"}]}, ImageSize -> 440],
  Style[Column[{Row[{"n = 3:   ", Row[Row[{TableauDimension[#], " \[Times] ", TraditionalForm[Factor[SchurDimension[#, d]]]}] & /@ IntegerPartitions[3], "  +  "], "  =  ",
       TraditionalForm[Expand[Total[TableauDimension[#] SchurDimension[#, d] & /@ IntegerPartitions[3]]]]}],
     Row[{"n = 5:   ", Row[MapThread[Row[{Subscript["T", #1], " = ", youngShow[#2]}] &, {Range[4], Flatten[youngPairs]}], ",   "]}],
     Row[{"integer \[Psi], d = 5:   \[LeftDoubleBracketingBar]", Subscript["\[Sum]", "T"], Subscript["P", "T"], "\[Psi] \[Minus] \[Psi] \[Minus] ", youngP[1], youngP[2], "\[Psi] \[Minus] ", youngP[3], youngP[4], "\[Psi]\[RightDoubleBracketingBar] = ", youngExcess}]},
    Spacings -> 1.2], 13, FontFamily -> "Helvetica"], Top]
```

At $n=3$ the blocks hold $d^3$ states in all; for three qubits they are a spin-$\frac32$ quartet, two spin-$\frac12$ doublets and no antisymmetric state. Through $n=4$ the projections add back to the state, for qubits and for qutrits; from $n=5$ to $n=7$ a sizable part is missing. At $n=5$ the excess is exactly $\hat P_{T_1}\hat P_{T_2}\psi+\hat P_{T_3}\hat P_{T_4}\psi$, so $\sum_T\hat P_T=\mathbb 1+\hat P_{T_1}\hat P_{T_2}+\hat P_{T_3}\hat P_{T_4}$: two pairs of standard tableaux of the same shape whose projectors do not annihilate each other. The pair of shape $\{2,2,1\}$ needs three levels, so for qubits only the first product remains.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Dicke and W states from the symmetric Young projector

$n$ qubits start in a computational basis state $|x\rangle$ with $k$ excitations: $k$ sites in $|1\rangle$ and the rest in $|0\rangle$. A measurement asks whether the state is symmetric under every permutation of the sites, that is, it projects with the Young projector of the one-row diagram. With what probability does the state pass, and what comes out?

$$
\hat P_{\{n\}}=\frac{1}{n!}\sum_{\pi\in S_n}\hat\pi,\qquad |D_n^k\rangle=\binom{n}{k}^{-1/2}\sum_{y:\;k\text{ ones}}|y\rangle,\qquad \hat S_z|D_n^k\rangle=\left(\frac n2-k\right)|D_n^k\rangle,\qquad \hat{\vec S}^{\,2}|D_n^k\rangle=\frac n2\left(\frac n2+1\right)|D_n^k\rangle,
$$

where $\hat\pi$ permutes the sites, the sum runs over all $\binom nk$ strings $y$ with $k$ ones, $\hat{\vec S}=\frac12\sum_i\hat{\vec\sigma}_i$ is the total spin with $\hat\sigma_z|0\rangle=|0\rangle$, and $|D_n^1\rangle$ is the W state. The image of $\hat P_{\{n\}}$ is the spin-$\frac n2$ multiplet, $n+1$ states out of $2^n$, spanned by the Dicke states $|D_n^k\rangle$.

The figure gives the passing probability $\|\hat P_{\{n\}}|x\rangle\|^2$ as markers for $n=4,6,8$ and every $k$, with $\binom nk^{-1}$ as lines, and the label gives the normalized output for $n=5$, $k=1$:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`Symmetry`"];
dickeString[n_, k_] := Normal[SparseArray[Join[ConstantArray[2, k], ConstantArray[1, n - k]] -> 1, ConstantArray[2, n]]];
dickeOut[n_, k_] := Flatten[YoungProject[dickeString[n, k], YoungTableau[{n}]]];
dickeNs = {4, 6, 8};
dickeW = With[{v = dickeOut[5, 1]}, With[{pos = Flatten[Position[v, _?Positive]]},
    Row[{TraditionalForm[v[[First[pos]]]/Norm[v]], " (", Row[Ket[{StringJoin[ToString /@ IntegerDigits[# - 1, 2, 5]]}] & /@ pos, " + "], ")"}]]];
Labeled[ListLogPlot[Join[Table[{k, 1/Binomial[n, k]}, {n, dickeNs}, {k, 0, n}], Table[{k, Norm[dickeOut[n, k]]^2}, {n, dickeNs}, {k, 0, n}]],
   Joined -> {True, True, True, False, False, False}, PlotStyle -> ColorData[97] /@ {1, 2, 3, 1, 2, 3}, PlotMarkers -> {None, None, None, {Automatic, 11}, {Automatic, 11}, {Automatic, 11}},
   PlotLegends -> {Row[{"n = ", 4}], Row[{"n = ", 6}], Row[{"n = ", 8}], None, None, None}, PlotRange -> {{-0.3, 8.3}, {0.01, 1.5}},
   Frame -> True, FrameLabel -> {"number of excitations k", "probability to pass"}, ImageSize -> 440],
  Style[Column[{Row[{"markers: \[LeftDoubleBracketingBar]", Subscript["P", "{n}"], Ket[{"1\[Ellipsis]10\[Ellipsis]0"}], Superscript["\[RightDoubleBracketingBar]", 2], ",   lines: ", TraditionalForm[HoldForm[1/Binomial[n, k]]]}],
     Row[{"n = 5, k = 1:   ", Subscript["P", "{5}"], Ket[{"10000"}], " / \[LeftDoubleBracketingBar]\[CenterEllipsis]\[RightDoubleBracketingBar]  =  ", dickeW}]}, Spacings -> 1.2], 13, FontFamily -> "Helvetica"], Top]
```

The projector averages $|x\rangle$ over all its rearrangements, and each of the $\binom nk$ distinct strings with $k$ ones appears equally often, so the output is the Dicke state $|D_n^k\rangle$, the W state with amplitudes $1/\sqrt5$ for $n=5$, $k=1$, and the state passes with probability $\binom nk^{-1}$, as the markers confirm at every $n$ and $k$ plotted. Only $k=0$ and $k=n$, already symmetric, pass for sure; the probability is smallest at half filling, $k=n/2$, and falls there as $n$ grows.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The AKLT ring: from a matrix-product state to an exact ground state

A ring of $n$ spin-1 sites is given as a matrix-product state whose bonds carry only two states: one $2\times2\times3$ tensor, repeated on every site, fixes all $3^n$ amplitudes. Contracted into those amplitudes, is it an exact ground state of the AKLT Hamiltonian, and how fast do the spin correlations $C(r)=\langle S^z_1S^z_{1+r}\rangle$ decay with the distance $r$?

With $\mathbf S_i$ the spin-1 operators of site $i$ and $P_2^{(i,i+1)}$ the projector of two neighbors onto total spin 2, the Hamiltonian on the ring ($n+1\equiv1$) and the state, with $m_i\in\{1,0,-1\}$ and $\sigma^\pm$, $\sigma^z$ the $2\times2$ Pauli matrices, are

$$
\hat H=\sum_{i=1}^{n}P_2^{(i,i+1)}=\sum_{i=1}^{n}\Big[\frac12\,\mathbf S_i\cdot\mathbf S_{i+1}+\frac16\,(\mathbf S_i\cdot\mathbf S_{i+1})^2+\frac13\Big],\qquad
|\psi\rangle=\sum_{m_1\dots m_n}\operatorname{Tr}\big(A^{m_1}\cdots A^{m_n}\big)\,|m_1\cdots m_n\rangle,
$$

$$
A^{1}=\sqrt{\frac23}\,\sigma^+,\qquad A^{0}=-\sqrt{\frac13}\,\sigma^z,\qquad A^{-1}=-\sqrt{\frac23}\,\sigma^- .
$$

$\hat H$ is a sum of projectors, so $\langle\hat H\rangle\ge0$ and any state with $\langle\hat H\rangle=0$ is a ground state. The transfer matrix $\sum_mA^m\otimes\bar A^m$ has eigenvalues $1$ and $-\frac13$ (three times), which gives on the ring

$$
C(r)=\frac43\,\frac{(-1/3)^r+(-1/3)^{n-r}}{1+3\,(-1/3)^n}.
$$

The network is contracted exactly for rings of $n=6$, $8$ and $10$ sites; the markers are $|C(r)|$ computed from the amplitudes for $r\le n/2$, the dashed line is the infinite-chain law $\frac43\,3^{-r}$, and the energy is in the plot label:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`"];
akltA = Transpose[{Sqrt[2/3] {{0, 1}, {0, 0}}, -Sqrt[1/3] {{1, 0}, {0, -1}}, -Sqrt[2/3] {{0, 0}, {1, 0}}}, {3, 1, 2}];
akltPsi[n_] := TensorNetworkContract[TensorNetwork[ConstantArray[akltA, n], Table[{i, Mod[i, n] + 1, n + i}, {i, n}]]];
akltSS = Total[KroneckerProduct[#, #] & /@ {{{0, 1, 0}, {1, 0, 1}, {0, 1, 0}}/Sqrt[2], {{0, -I, 0}, {I, 0, -I}, {0, I, 0}}/Sqrt[2], DiagonalMatrix[{1, 0, -1}]}];
akltP2 = (akltSS + IdentityMatrix[9]) . (akltSS + 2 IdentityMatrix[9])/6;
akltData[n_] := With[{psi = akltPsi[n]}, With[{m = ArrayReshape[psi, {9, 3^(n - 2)}], p = Abs[psi]^2/Total[Abs[psi]^2, Infinity]},
    {n Tr[akltP2 . m . ConjugateTranspose[m]]/Total[Abs[psi]^2, Infinity],
     Table[{r, Abs[Total[Total[Flatten[p, {{1}, {r + 1}, Complement[Range[n], {1, r + 1}]}], {3}] Outer[Times, {1, 0, -1}, {1, 0, -1}], 2]]}, {r, n/2}]}]];
akltRes = akltData /@ {6, 8, 10};
ListLogPlot[Append[akltRes[[All, 2]], Table[{r, 4/3 3^-r}, {r, 1, 5}]], Joined -> {False, False, False, True}, PlotStyle -> {Automatic, Automatic, Automatic, Dashed},
  PlotMarkers -> {Automatic, Automatic, Automatic, None}, PlotLegends -> Append[Row[{"n = ", #}] & /@ {6, 8, 10}, Row[{TraditionalForm[4/3], " ", TraditionalForm[HoldForm[3^-r]]}]],
  Frame -> True, FrameLabel -> {"r", Row[{"|\[LeftAngleBracket]", Subsuperscript["S", 1, "z"], Subsuperscript["S", "1+r", "z"], "\[RightAngleBracket]|"}]},
  PlotLabel -> Row[{"\[LeftAngleBracket]", OverHat["H"], "\[RightAngleBracket] = ", Row[akltRes[[All, 1]], ", "], "  for n = 6, 8, 10"}], ImageSize -> 460]
```

The contracted state has energy exactly zero on every ring, so the two-state bonds hold an exact ground state of the AKLT chain. The correlations fall by a factor of three per site with alternating sign, a correlation length of $1/\ln3$ sites, and they follow the infinite-chain line except near $r=n/2$, where the correlation reaching the other way round the ring, the term $(-1/3)^{n-r}$ of the ring formula, adds to it and lifts the points above the line.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The entanglement profile of an Ising ring, and the bond dimension it needs

A ring of $n$ spins in its ground state is cut into a block of $\ell$ neighboring sites and the rest. How does the entanglement entropy $S(\ell)$ of the block depend on $\ell$, in the ordered phase, at the critical point and in the paramagnet? The answer sets the bond dimension $\chi$ a matrix-product state needs, because the reduced state of the block has rank at least $e^{S(\ell)}$, so $\chi\ge e^{S(\ell)}$ at that cut.

With $Z_i$, $X_i$ the Pauli operators of site $i$ and $Z_{n+1}\equiv Z_1$, the transverse-field Ising ring is

$$
\hat H=-J\sum_{i=1}^{n}Z_iZ_{i+1}-h\sum_{i=1}^{n}X_i ,\qquad S(\ell)=-\operatorname{Tr}\rho_\ell\ln\rho_\ell ,
$$

with $\rho_\ell$ the reduced state of sites $1,\dots,\ell$. At the critical point $h=J$ conformal field theory predicts (Calabrese and Cardy 2004)

$$
S(\ell)=\frac c3\,\ln\!\Big[\frac n\pi\sin\frac{\pi\ell}{n}\Big]+s_1 ,
$$

with central charge $c$ and a constant $s_1$.

The ground state of a ring of $n=16$ spins is found exactly at $h/J=\frac12$, $1$ and $2$ and written as a matrix-product state with no truncation; $S(\ell)$ is read from its Schmidt values at each bond. The dashed lines are a straight-line fit to the critical points, whose slope gives $c$, and the level $\ln2$:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`"];
epN = 16;
epSite[ops_Association] := KroneckerProduct @@ Table[Lookup[ops, k, IdentityMatrix[2, SparseArray]], {k, epN}];
epH[hJ_] := -Sum[epSite[<|i -> SparseArray[{{1, 0}, {0, -1}}], Mod[i, epN] + 1 -> SparseArray[{{1, 0}, {0, -1}}]|>], {i, epN}] - hJ Sum[epSite[<|i -> SparseArray[{{0, 1}, {1, 0}}]|>], {i, epN}];
epGround[hJ_] := First[Eigenvectors[-N[epH[hJ]], 1, Method -> {"Arnoldi", "Criteria" -> "RealPart"}]];
epStep[{_, r_}, _] := With[{m = ArrayReshape[r, {2 Length[r], Length[First[r]]/2}]},
   With[{svd = SingularValueDecomposition[m, Min[Dimensions[m]]]}, {svd[[1]], svd[[2]] . ConjugateTranspose[svd[[3]]]}]];
epMPS[psi_] := With[{steps = Rest[FoldList[epStep, {None, {psi}}, Range[epN - 1]]]},
   TensorNetwork[Join[{Transpose[steps[[1, 1]]]}, Transpose[ArrayReshape[#, {Length[#]/2, 2, Last[Dimensions[#]]}], {1, 3, 2}] & /@ steps[[2 ;;, 1]], {steps[[-1, 2]]}],
    Join[{{1, epN + 1}}, Table[{i - 1, i, epN + i}, {i, 2, epN - 1}], {{epN - 1, 2 epN}}]]];
epChord = Log[(epN/Pi) Sin[Pi Range[epN - 1]/epN]];
epProfiles = Table[Transpose[{epChord, With[{mps = epMPS[epGround[hJ]]}, Table[MPSEntanglementEntropy[mps, l], {l, epN - 1}]]}], {hJ, {1/2, 1, 2}}];
epFit = Fit[epProfiles[[2]], {1, x}, x];
ListPlot[Join[epProfiles, {Table[{x, epFit}, {x, MinMax[epChord]}], Table[{x, Log[2]}, {x, MinMax[epChord]}]}], Joined -> {False, False, False, True, True},
  PlotStyle -> {Automatic, Automatic, Automatic, Dashed, Dashed}, PlotMarkers -> {Automatic, Automatic, Automatic, None, None}, PlotRange -> {0, 0.85},
  PlotLegends -> {Row[{"h/J = ", 1/2}], "h/J = 1", "h/J = 2", Row[{"fit at h/J = 1: c = ", NumberForm[3 Coefficient[epFit, x], 2]}], TraditionalForm[HoldForm[Log[2]]]},
  Frame -> True, FrameLabel -> {TraditionalForm[HoldForm[Log[n/\[Pi] Sin[\[Pi] \[ScriptL]/n]]]], Row[{"S(\[ScriptL])"}]}, PlotLabel -> "n = 16", ImageSize -> 480]
```

At the critical point the entropy grows along a straight line in $\ln[\frac n\pi\sin\frac{\pi\ell}n]$ whose slope gives $c$ close to $\frac12$, the central charge of the Ising critical point, so the bond dimension a critical chain needs grows as a power of the block size, $\chi\ge e^{S}\propto\ell^{1/6}$. In the paramagnet and in the ordered phase the entropy hardly depends on $\ell$: it stays small at $h/J=2$, and at $h/J=\frac12$ it sits near $\ln2$, the one bit of the cat state $\frac1{\sqrt2}(|{\uparrow\cdots\uparrow}\rangle+|{\downarrow\cdots\downarrow}\rangle)$ that the symmetric ground state of a finite ring forms. A gapped chain therefore needs a bond dimension that does not grow with the block.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Truncating a quenched Ising chain: does the error equal the discarded Schmidt weight?

An open chain of $n$ spins starts in $|0\cdots0\rangle$ and evolves under the critical Ising Hamiltonian. At time $t$ its matrix-product state is cut to bond dimension $\chi$ at every bond. How fast does the infidelity $1-F$ fall with $\chi$, and does it equal the Schmidt weight the cuts discard?

With $\hbar=1$, $Z_i$, $X_i$ the Pauli operators of site $i$, and $\lambda^{(k)}_1\ge\lambda^{(k)}_2\ge\dots$ the Schmidt values of $|\psi(t)\rangle$ at the bond between sites $k$ and $k+1$,

$$
\hat H=-J\sum_{i=1}^{n-1}Z_iZ_{i+1}-J\sum_{i=1}^{n}X_i ,\qquad |\psi(t)\rangle=e^{-i\hat Ht}|0\cdots0\rangle ,\qquad
F=|\langle\psi(t)|\psi_\chi\rangle|^2 ,\qquad \varepsilon_\chi=\sum_{k=1}^{n-1}\sum_{\alpha>\chi}\big(\lambda^{(k)}_\alpha\big)^2 ,
$$

where $|\psi_\chi\rangle$ is the normalized truncated state. The truncation sweeps from left to right and keeps the $\chi$ largest singular values at each bond. Those are the Schmidt values only when every tensor to the right of the bond is an isometry, the right-canonical form, so the state is brought to that form first.

The state of $n=14$ spins is computed exactly at $Jt=\frac12$, $1$ and $2$; filled circles are $1-F$ after truncation from the right-canonical form, the dashed line is $\varepsilon_\chi$, and open circles are $1-F$ when the same sweep starts from the left-canonical form:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`"];
qtN = 14;
qtSite[ops_Association] := KroneckerProduct @@ Table[Lookup[ops, k, IdentityMatrix[2, SparseArray]], {k, qtN}];
qtH = -Sum[qtSite[<|i -> SparseArray[{{1, 0}, {0, -1}}], i + 1 -> SparseArray[{{1, 0}, {0, -1}}]|>], {i, qtN - 1}] - Sum[qtSite[<|i -> SparseArray[{{0, 1}, {1, 0}}]|>], {i, qtN}];
qtStep[{_, r_}, _] := With[{m = ArrayReshape[r, {2 Length[r], Length[First[r]]/2}]},
   With[{svd = SingularValueDecomposition[m, Min[Dimensions[m]]]}, {svd[[1]], svd[[2]] . ConjugateTranspose[svd[[3]]]}]];
qtMPS[psi_] := With[{steps = Rest[FoldList[qtStep, {None, {psi}}, Range[qtN - 1]]]},
   TensorNetwork[Join[{Transpose[steps[[1, 1]]]}, Transpose[ArrayReshape[#, {Length[#]/2, 2, Last[Dimensions[#]]}], {1, 3, 2}] & /@ steps[[2 ;;, 1]], {steps[[-1, 2]]}],
    Join[{{1, qtN + 1}}, Table[{i - 1, i, qtN + i}, {i, 2, qtN - 1}], {{qtN - 1, 2 qtN}}]]];
qtPanel[Jt_] := With[{left = qtMPS[MatrixExp[-I N[Jt] qtH, N[UnitVector[2^qtN, 1]]]]},
   With[{right = MPSCanonicalForm[left, "Right"], schmidt = Table[MPSSchmidtValues[left, k], {k, qtN - 1}]},
    ListLogPlot[Select[#, Last[#] > 10^-13 &] & /@ {
       Table[{chi, 1 - Abs[MPSOverlap[MPSTruncate[right, chi], left]]^2}, {chi, 12}],
       Table[{chi, Total[Total[Drop[#, UpTo[chi]]^2] & /@ schmidt]}, {chi, 12}],
       Table[{chi, 1 - Abs[MPSOverlap[MPSTruncate[left, chi], left]]^2}, {chi, 12}]},
     Joined -> {False, True, False}, PlotStyle -> {Automatic, Dashed, Automatic}, PlotMarkers -> {{"\[FilledCircle]", 10}, None, {"\[EmptyCircle]", 14}},
     PlotLegends -> {"1 \[Minus] F, right-canonical", Row[{Subscript["\[CurlyEpsilon]", "\[Chi]"], ", discarded weight"}], "1 \[Minus] F, left-canonical"},
     PlotRange -> {10^-13, 2}, Frame -> True, FrameLabel -> {"\[Chi]", "error"}, PlotLabel -> Row[{"J t = ", Jt}], ImageSize -> 300]]];
qtPanels = qtPanel /@ {1/2, 1, 2};
Legended[PlotGrid[{First /@ qtPanels}, ImageSize -> 900], legendBelow[Last[First[qtPanels]]]]
```

From the right-canonical form the infidelity sits on, or just below, the summed discarded weight at every $\chi$ and every time, so the Schmidt values predict the error before any truncation is done. The error falls steeply with $\chi$, and the later the time, the more slowly: the entanglement the quench builds up raises the bond dimension a given accuracy needs. Started from the left-canonical form, where the singular values the sweep sees are not the Schmidt values, the same bond dimensions lose most of the state.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What it costs to compute the norm of a matrix-product state

A matrix-product state on $n$ sites with physical dimension $d$ and bond dimension $\chi$ has the norm $\langle\psi|\psi\rangle$, a closed network of $2n$ tensors. In which order should it be contracted, and how does the cost of the cheapest order grow with $n$ and $\chi$, compared with first building the $d^n$ amplitudes?

With $A^{(k)}$ the tensor of site $k$, bond indices $a_k$ of range $\chi$ and physical indices $m_k$ of range $d$,

$$
\langle\psi|\psi\rangle=\sum_{\{a\},\{a'\},\{m\}}A^{(1)}_{a_1m_1}A^{(2)}_{a_1a_2m_2}\cdots A^{(n)}_{a_{n-1}m_n}\;\bar A^{(1)}_{a'_1m_1}\bar A^{(2)}_{a'_1a'_2m_2}\cdots\bar A^{(n)}_{a'_{n-1}m_n} .
$$

Contracting two tensors costs the product of the ranges of all their indices, in multiply-adds; the cost of an order is the sum over its steps. Sweeping a $\chi\times\chi$ boundary matrix through the chain, one ket tensor and one bra tensor at a time, costs $2d\chi^3$ per site.

Only the shapes of the tensors decide the order, so the tensors are left as symbols. For $d=2$ and $\chi=4$ and $16$, filled circles are the cost of the cheapest order found by exhaustive search, open circles the cost of contracting the ket and the bra into $d^n$ amplitudes each and then taking their product, both divided by $2d\chi^3$; the dashed line is $n$:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`"];
ncNet[n_, chi_, d_] := TensorNetwork[Join @@ Table[Join[{ArraySymbol[s, {chi, d}]}, ConstantArray[ArraySymbol[s, {chi, chi, d}], n - 2], {ArraySymbol[s, {chi, d}]}], {s, {"A", "Abar"}}],
   Join @@ Table[Join[{{o + 1, 2 n - 1}}, Table[{o + i - 1, o + i, 2 n - 2 + i}, {i, 2, n - 1}], {{o + n - 1, 3 n - 2}}], {o, {0, n - 1}}], {}];
ncCost[net_, path_] := With[{dim = net["IndexDimensions"]},
   First[Fold[Function[{acc, step}, With[{all = Union @@ acc[[2]][[step]], rest = Delete[acc[[2]], List /@ step]},
       {acc[[1]] + Times @@ Lookup[dim, all], Append[rest, Select[all, MemberQ[Catenate[rest], #] &]]}]], {0, net["Hyperedges"]}, path]]];
ncVectorFirst[n_] := Join[{{1, 2}}, Table[{1, 2 n - k}, {k, n - 2}], {{1, 2}}, Table[{1, n - k + 1}, {k, n - 2}], {{1, 2}}];
ncPanel[chi_] := With[{ns = Range[4, 24, 2]},
   ListPlot[{Table[{n, ncCost[ncNet[n, chi, 2], OptimalContractionPath[ncNet[n, chi, 2], Method -> "flops"]]/(4 chi^3)}, {n, ns}],
     Table[{n, ncCost[ncNet[n, chi, 2], ncVectorFirst[n]]/(4 chi^3)}, {n, ns}], Table[{n, n}, {n, MinMax[ns]}]},
    Joined -> {False, False, True}, PlotStyle -> {Automatic, Automatic, Dashed}, PlotMarkers -> {{"\[FilledCircle]", 10}, {"\[EmptyCircle]", 14}, None},
    PlotLegends -> {"cheapest order", Row[{"amplitudes first"}], TraditionalForm[HoldForm[n]]}, PlotRange -> {0, 30}, Frame -> True,
    FrameLabel -> {"n", Row[{"cost / ", TraditionalForm[HoldForm[2 d \[Chi]^3]]}]}, PlotLabel -> Row[{"\[Chi] = ", chi, ",  d = 2"}], ImageSize -> 340]];
ncPanels = ncPanel /@ {4, 16};
Legended[PlotGrid[{First /@ ncPanels}, ImageSize -> 680], legendBelow[Last[First[ncPanels]]]]
```

Once the chain is long enough that $\chi\le d^{n/2}$, the cheapest order runs parallel to the dashed line: each added site costs exactly $2d\chi^3$, so the norm costs time linear in $n$ and cubic in $\chi$. Building the amplitudes first costs about as much only on the shortest chains; its cost, divided by $2d\chi^3$, grows as $d^n/\chi$ and leaves the plot within a few sites, later for the larger $\chi$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### How close GHZ and W states come to a product state

The geometric measure of entanglement asks how far a state is from the nearest product state. For the $n$-qubit GHZ and W states, what is the largest overlap $\Lambda^2$ with any product state, and how does it depend on $n$?

$$
|\mathrm{GHZ}_n\rangle=\frac{|0\cdots0\rangle+|1\cdots1\rangle}{\sqrt2},\qquad |W_n\rangle=\frac{|10\cdots0\rangle+|010\cdots0\rangle+\dots+|0\cdots01\rangle}{\sqrt n},\qquad
\Lambda^2(\psi)=\max_{\phi_1,\dots,\phi_n}\big|\langle\phi_1\cdots\phi_n|\psi\rangle\big|^2 ,
$$

where the overlap is the contraction $\langle\phi_1\cdots\phi_n|\psi\rangle=\sum_{i_1\dots i_n}\bar\phi^{(1)}_{i_1}\cdots\bar\phi^{(n)}_{i_n}\,\psi_{i_1\cdots i_n}$. Both states are symmetric under exchange of qubits, so a nearest product state can be taken symmetric, $|\phi\rangle^{\otimes n}$ (Hübener et al. 2009), and since their amplitudes are nonnegative, $|\phi\rangle=\cos\frac\theta2|0\rangle+\sin\frac\theta2|1\rangle$ with $0\le\theta\le\pi$.

The overlap is contracted exactly for general $\theta$ and maximized exactly for $n=2,\dots,8$; the labels give the exact maxima, and the dashed line is $1/e$:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`"];
geoGHZ[n_] := ArrayReshape[UnitVector[2^n, 1] + UnitVector[2^n, 2^n], ConstantArray[2, n]]/Sqrt[2];
geoW[n_] := ArrayReshape[Total[UnitVector[2^n, 2^# + 1] & /@ Range[0, n - 1]], ConstantArray[2, n]]/Sqrt[n];
geoOverlap[psi_] := With[{n = ArrayDepth[psi]},
   ActivateTensors[EinsteinSummation[Append[List /@ Range[n], Range[n]] -> {}, Append[ConstantArray[{Cos[\[Theta]/2], Sin[\[Theta]/2]}, n], psi]]]];
geoMax[psi_] := RootReduce[First[Maximize[{geoOverlap[psi]^2, 0 <= \[Theta] <= Pi}, \[Theta]]]];
geoTab = Table[{n, geoMax[geoGHZ[n]], geoMax[geoW[n]]}, {n, 2, 8}];
ListPlot[{geoTab[[All, {1, 2}]], Callout[{#1, #3}, TraditionalForm[#3], Below] & @@@ geoTab, {{2, 1/E}, {8, 1/E}}},
  Joined -> {False, False, True}, PlotStyle -> {Automatic, Automatic, Dashed}, PlotMarkers -> {{"\[FilledSquare]", 10}, {"\[FilledCircle]", 10}, None},
  PlotLegends -> {Row[{"|", Subscript["GHZ", "n"], "\[RightAngleBracket]"}], Row[{"|", Subscript["W", "n"], "\[RightAngleBracket]"}], TraditionalForm[HoldForm[1/E]]}, PlotRange -> {0.3, 0.55},
  Frame -> True, FrameLabel -> {"n", TraditionalForm[HoldForm[\[CapitalLambda]^2]]}, ImageSize -> 460]
```

The GHZ state reaches overlap exactly $\frac12$ for every $n$. The W state's maxima are exactly $\big(\frac{n-1}n\big)^{n-1}$, namely $\frac49,\frac{27}{64},\frac{256}{625},\dots$, falling toward $1/e$, so for every $n\ge3$ the W state is farther from all product states than the GHZ state, and at $n=2$ both are the same Bell state with $\Lambda^2=\frac12$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### How much entanglement each tensor-network family can carry

A tensor-network state of fixed bond dimension $\chi$ can hold only so much entanglement between a region $A$ of its sites and the rest. How does that maximum grow with the number of sites $|A|$ for a matrix-product state on a chain, a binary tree tensor network over a chain, and a projected entangled-pair state (PEPS) on a square grid?

Any cut through the network that separates the physical legs of $A$ from the others bounds the Schmidt rank by the product of the dimensions of the legs it crosses. With every leg, bond or physical, of dimension $\chi$,

$$
S(A)\le|\gamma_A|\,\ln\chi ,
$$

where $|\gamma_A|$ is the smallest number of legs any such cut crosses, the minimal cut.

The minimal cut is the maximum flow from the legs of $A$ to the rest, with one unit of capacity per leg. It is computed for a chain of 128 sites and a tree with 128 leaves, with $A$ a block of $|A|=1,2,4,\dots,64$ neighboring sites, and for a $16\times16$ PEPS with $A$ an $\ell\times\ell$ square, each averaged over every position of $A$ that does not touch the edge of the network; the dashed line is the perimeter of the square, $4\sqrt{|A|}$:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"]; Needs["Wolfram`TensorNetworks`"];
mcCut[edges_, region_] := With[{holders = GroupBy[Catenate[MapIndexed[Thread[#1 -> First[#2]] &, edges]], First -> Last]},
   With[{bonds = Values[Select[holders, Length[#] == 2 &]], phys = Keys[Select[holders, Length[#] == 1 &]]},
    FindMaximumFlow[Graph[Join[UndirectedEdge @@@ bonds, UndirectedEdge[First[holders[#]], {"leg", #}] & /@ phys,
       UndirectedEdge["A", {"leg", #}] & /@ region, UndirectedEdge[{"leg", #}, "rest"] & /@ Complement[phys, region]],
      EdgeCapacity -> Join[ConstantArray[1, Length[bonds] + Length[phys]], ConstantArray[Infinity, Length[phys]]]], "A", "rest"]]];
mcLegs[edges_] := Sort[Keys[Select[Counts[Catenate[edges]], # == 1 &]]];
mcBlocks[edges_] := With[{legs = mcLegs[edges]}, Table[{l, Mean[Table[mcCut[edges, legs[[a ;; a + l - 1]]], {a, 2, Length[legs] - l}]]}, {l, 2^Range[0, 6]}]];
mcMPS = mcBlocks[RandomTensorNetwork["MPS"[128, 2, 2]]["Hyperedges"]];
mcTree = mcBlocks[Select[RandomTensorNetwork["TTN"[8, 2, 2]]["Hyperedges"], Length[#] > 1 &]];
mcPEPS = With[{edges = RandomTensorNetwork["PEPS"[{16, 16}, 2, 2]]["Hyperedges"]}, With[{grid = Partition[mcLegs[edges], 16]},
    Table[{l^2, Mean[Catenate[Table[mcCut[edges, Catenate[grid[[r ;; r + l - 1, c ;; c + l - 1]]]], {r, 2, 16 - l}, {c, 2, 16 - l}]]]}, {l, 8}]]];
ListLogLinearPlot[{mcMPS, mcTree, mcPEPS, Table[{a, 4 Sqrt[a]}, {a, 1, 64}]}, Joined -> True, PlotStyle -> {Automatic, Automatic, Automatic, Dashed},
  PlotMarkers -> {Automatic, Automatic, Automatic, None}, PlotLegends -> {"matrix-product state", "binary tree", "PEPS, square region", TraditionalForm[HoldForm[4 Sqrt[Abs[A]]]]},
  Frame -> True, FrameLabel -> {"|A|", Row[{"|", Subscript["\[Gamma]", "A"], "|"}]}, ImageSize -> 460]
```

The three families separate by how their cuts grow. A block of two or more sites inside a matrix-product state is cut by exactly two bonds, whatever its size, so the state carries at most $2\ln\chi$: enough for a gapped chain, whose entropy saturates, but not for a critical chain, whose entropy grows as $\ln|A|$. The tree's cut grows by about the same amount each time the block doubles, the logarithmic growth of a critical chain. The PEPS cut is the perimeter $4\sqrt{|A|}$ for every square from $4\times4$ up, and smaller squares are cut more cheaply through their own physical legs: the area law of a gapped two-dimensional ground state.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Superconducting hardware: pulses, transmons and control

The same operators model a transmon's levels, pulses and couplings; TransmonLab models readout and QUALink compiles pulses to Quantum Machines control programs.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Leakage of a fast π pulse on a transmon


A transmon is an anharmonic oscillator, frequency $\omega_q$ and anharmonicity $\alpha<0$, so its $|1\rangle\to|2\rangle$ transition sits only $|\alpha|$ below the qubit transition. A π pulse of length $\tau$, detuned by $\delta$ from the qubit, is meant to flip $|0\rangle$ to $|1\rangle$, and a DRAG term of strength $\lambda$ adds a second quadrature that follows the slope of the envelope. Where does the population end up after the pulse: outside the qubit, or left in $|0\rangle$?

With $\hbar=1$, the driven transmon and its form in the frame rotating at the drive frequency $\omega_d=\omega_q+\delta$, after the rotating-wave approximation, are

$$
\hat H=\omega_q \hat b^\dagger \hat b+\frac{\alpha}{2}\hat b^{\dagger2}\hat b^2+\left[\Omega(t)\cos\omega_d t-\frac{\lambda\dot\Omega(t)}{\alpha}\sin\omega_d t\right](\hat b+\hat b^\dagger)\;\longrightarrow\;\frac{\alpha}{2}\hat b^{\dagger2}\hat b^2-\delta\,\hat b^\dagger \hat b+\frac{\Omega(t)}{2}(\hat b+\hat b^\dagger)-\frac{\lambda\dot\Omega(t)}{2\alpha}\,i(\hat b^\dagger-\hat b),
$$

with the raised-cosine envelope $\Omega(t)=\frac{\pi}{\tau}\left(1-\cos\frac{2\pi t}{\tau}\right)$, whose area is $\pi$, and the oscillator kept to its five lowest levels.

Starting from $|0\rangle$, the figure plots the leakage $P(n\ge2)$ (solid) and the population left in $|0\rangle$, $P(0)$ (dashed), during the pulse, and at $\delta=0$ also $P(0)$ for a two-level qubit with $\lambda=0$ (dotted); times and frequencies are in units of $|\alpha|$:

```wl
dragB = AnnihilationOperator[5]; dragN = SuperDagger[dragB] @ dragB;
dragEnvelope[τ_] := Pi/τ (1 - Cos[2 Pi t/τ]);
dragH[α_, τ_, λ_, δ_] := α/2 (dragN @ dragN - dragN) - δ dragN + dragEnvelope[τ]/2 (dragB + SuperDagger[dragB]) -
    λ D[dragEnvelope[τ], t]/(2 α) I (SuperDagger[dragB] - dragB);
dragB2 = AnnihilationOperator[2];
dragPops[τ_, λ_, δ_] := QuantumEvolve[dragH[-1, τ, λ, δ], FockState[0, 5], {t, 0, τ}][t]["ProbabilitiesList"];
dragλs = {0, 0.5, 1};
dragPanel[τ_, δ_] := With[{pops = dragPops[τ, #, δ] & /@ dragλs},
   Show[LogPlot[Evaluate[Total[#[[3 ;;]]] & /@ pops], {t, 0, τ}, PlotRange -> All, PlotLegends -> (Row[{"λ = ", #}] & /@ dragλs)],
    LogPlot[Evaluate[First /@ pops], {t, 0, τ}, PlotRange -> All, PlotStyle -> Dashed],
    If[δ == 0, LogPlot[Evaluate[First @ QuantumEvolve[dragEnvelope[τ]/2 (dragB2 + SuperDagger[dragB2]), FockState[0, 2], {t, 0, τ}][t]["ProbabilitiesList"]],
      {t, 0, τ}, PlotRange -> All, PlotStyle -> Dotted], Graphics[]], PlotRange -> {Log[10^-9], Log[2]}, Frame -> True,
    FrameLabel -> {"|α| t", "population"}, PlotLabel -> Row[{"|α|τ = ", τ, ",   δ/|α| = ", δ}], ImageSize -> 300]];
dragPanels = AssociationThread[{5, 10, 20}, Outer[dragPanel, {5, 10, 20}, {-0.05, 0, 0.05}]];
Manipulate[Legended[PlotGrid[{First /@ dragPanels[τv]}, ImageSize -> 900], legendBelow[Last[dragPanels[5][[1]]]]], Row[{Control[{{τv, 10, "|α|τ"}, Keys[dragPanels], SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"|α|τ = ", τv}], <|5 -> "Shortest pulse: without DRAG, leakage grows to the end.",
        10 -> "Leakage peaks and partly returns; the two shifts nearly cancel at λ = 1/2.",
        20 -> "The lowest leakage; the shifts nearly cancel at λ = 1/2."|>[τv]}}]]}], SaveDefinitions -> True]
```

The drive couples only neighboring levels, so leakage builds up through $|1\rangle$ during the pulse; only the shortest pulse without DRAG keeps leaking to the end, and otherwise it peaks and partly returns. At $\lambda=1$ the population left outside the qubit drops by orders of magnitude, and a longer pulse lowers it further. Most of the remaining error is population left in $|0\rangle$, which a resonant π pulse on a two-level qubit would empty (dotted). The qubit transition is shifted by the $|2\rangle$ level and, in the opposite direction, by the DRAG quadrature itself; for the longer pulses the two shifts nearly cancel at $\lambda=1/2$, and a detuning partly undoes what is left, a positive one without DRAG and a negative one at $\lambda=1$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The CZ gate from two coupled transmons


Two transmons, frequencies $\omega_j$ and anharmonicities $\alpha_j$ (equal, $\alpha_1=\alpha_2=\alpha$, in the figures), exchange excitations at coupling $g$. A flux pulse tunes transmon 1 until $|11\rangle$ is detuned by $\delta$ from $|02\rangle$, holds it there for a time $t$, and tunes it back, both tunings taken as instant. How does this hold become a controlled-Z, and what does it need?

With $\hbar=1$ and three levels per transmon, the pair Hamiltonian is

$$
\hat H=\sum_{j=1,2}\left[\omega_j\,\hat n_j+\frac{\alpha_j}{2}\,\hat n_j(\hat n_j-1)\right]+g\left(\hat b_1^\dagger \hat b_2+\hat b_2^\dagger \hat b_1\right),\qquad \hat n_j=\hat b_j^\dagger \hat b_j,\qquad \delta=\omega_1-\omega_2-\alpha_2 .
$$

It conserves the total excitation number, so $|01\rangle$ and $|10\rangle$ couple only to each other, and $|11\rangle$ only to $|02\rangle$ and $|20\rangle$. Each transmon alone gives its states phases; the gate is what remains once they are removed, which for $|11\rangle$ is the conditional amplitude $c=U_{11,11}\,e^{-i\arg(U_{01,01}U_{10,10}/U_{00,00})}$.

On resonance, $\delta=0$, the figure follows each computational state on all nine levels for one full swap, $t=\pi/(\sqrt2\,g)$, with $\omega_2=0$ and rates in units of $g$:

```wl
czLower = DiagonalMatrix[Sqrt[{1, 2}], 1];
{czB1, czB2} = {KroneckerProduct[czLower, IdentityMatrix[3]], KroneckerProduct[IdentityMatrix[3], czLower]};
czDuffing[bj_, w_, a_] := With[{nj = Transpose[bj] . bj}, w nj + a/2 nj . (nj - IdentityMatrix[9])];
czH = czDuffing[czB1, ω1, α1] + czDuffing[czB2, ω2, α2] + g (Transpose[czB1] . czB2 + Transpose[czB2] . czB1);
czIndex = AssociationMap[1 + FromDigits[#, 3] &, {{0, 0}, {0, 1}, {1, 0}, {1, 1}, {0, 2}, {2, 0}}];
czU[d_, a_, s_] := MatrixExp[-I s N[czH /. {ω2 -> 0, α1 -> a, α2 -> a, ω1 -> a + d, g -> 1}]];
czCond[u_] := u[[czIndex[{1, 1}], czIndex[{1, 1}]]] Exp[-I Arg[u[[czIndex[{0, 1}], czIndex[{0, 1}]]] u[[czIndex[{1, 0}], czIndex[{1, 0}]]]/u[[czIndex[{0, 0}], czIndex[{0, 0}]]]]];
czSwap = N[π/Sqrt[2]];
czTimes = Subdivide[0, czSwap, 200];
czKet[l_] := Ket[{Row[l]}];
czDynamics[a_] := With[{us = czU[0, a, #] & /@ czTimes},
   {Transpose[{czTimes, #}] & /@ Join[
      Table[Abs[#[[czIndex[f], czIndex[{1, 1}]]]]^2 & /@ us, {f, {{1, 1}, {0, 2}, {2, 0}}}],
      {Abs[#[[czIndex[{0, 1}], czIndex[{0, 1}]]]]^2 & /@ us}],
    Transpose[{czTimes, #}] & /@ Transpose[ReIm[czCond /@ us]]}];
czPopPanel[a_] := ListLinePlot[czDynamics[a][[1]], PlotStyle -> {Automatic, Automatic, Automatic, Dashed}, PlotRange -> {-0.02, 1.02},
   PlotLegends -> Placed[LineLegend[{"P(|11⟩ → |11⟩)", "P(|11⟩ → |02⟩)", "P(|11⟩ → |20⟩)", "P(|01⟩ → |01⟩) = P(|10⟩ → |10⟩)"}, LegendLayout -> {"Column", 2}], Below],
   Frame -> True, FrameLabel -> {"g t", "transition probability"}, ImagePadding -> {{50, 10}, {40, 10}}, ImageSize -> 330];
czAmpPanel[a_] := ListLinePlot[czDynamics[a][[2]], PlotRange -> {-1.05, 1.05},
   PlotLegends -> Placed[LineLegend[{"Re c", "Im c"}, LegendLayout -> "Row"], Below], GridLines -> {None, {0}}, Frame -> True, FrameLabel -> {"g t", "amplitude c of |11⟩"}, ImagePadding -> {{50, 10}, {40, 10}}, ImageSize -> 330];
czStep1 = AssociationMap[{czPopPanel[#], czAmpPanel[#]} &, {-3, -10, -40}];
Manipulate[Grid[{czStep1[czα]}, Alignment -> Top], Row[{Control[{{czα, -10, "α/g"}, Keys[czStep1], SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"α/g = ", czα}], <|-3 -> "|20⟩ joins in, |11⟩ does not fully return, |01⟩ and |10⟩ exchange.",
        -10 -> "|11⟩ swaps to |02⟩ and back, returning as c = -1.",
        -40 -> "A clean swap: |11⟩ returns close to c = -1."|>[czα]}}]]}], SaveDefinitions -> True]
```

From $|11\rangle$ the population moves to $|02\rangle$ and back, while $|01\rangle$ and $|10\rangle$ stay where they are. The amplitude $c$ stays nearly real, passes through zero when the state is entirely in $|02\rangle$, and returns close to $-1$ at $\alpha/g=-40$. When $|\alpha|$ is comparable to $g$ this blurs: $|20\rangle$ joins in, $|11\rangle$ does not fully return, and $|01\rangle$ and $|10\rangle$, now close in energy, exchange population.

Dropping $|20\rangle$ leaves a two-level system on $\{|11\rangle,|02\rangle\}$, coupled by $\sqrt2\,g$ and detuned by $\delta$. One full swap turns its Bloch vector once around, and a two-level system turned by $2\pi$ returns with a factor $-1$ for any $\delta$; off resonance the energy offset of $|11\rangle$ adds a phase $\delta t/2$, so after one swap

$$
\phi=\pi\left(1+\frac{\delta}{\sqrt{\delta^2+8g^2}}\right).
$$

The figure draws the turn on the Bloch sphere of the pair, and the phase after one swap against $\delta$, from this closed form (gray) and from all nine levels (points):

```wl
czBlock[d_] := {{d, Sqrt[2]}, {Sqrt[2], 0}};
czBlockPath[d_] := With[{T = 2 Pi/Sqrt[d^2 + 8]}, Table[QuantumState[MatrixExp[-I s N[czBlock[d]]] . {1, 0}]["BlochVector"], {s, Subdivide[0, T, 150]}]];
czδs = {0, 1, 3};
czBloch = Graphics3D[{{Opacity[0.12], Sphere[]}, {GrayLevel[0.6], Line[{{0, 0, -1}, {0, 0, 1}}], Line[Table[{Cos[u], Sin[u], 0}, {u, 0, 2 Pi, Pi/60}]]},
    MapIndexed[{ColorData[97][First[#2]], AbsoluteThickness[3], Line[czBlockPath[#1]]} &, czδs],
    Text[Style["|11⟩", 15], {0, 0, 1.22}], Text[Style["|02⟩", 15], {0, 0, -1.22}]},
   Boxed -> False, ViewPoint -> {2.6, -1.4, 0.8}, Lighting -> "Neutral", ImageSize -> 300];
czPhase[d_] := 1 + d/Sqrt[d^2 + 8];
czPhase9[d_, a_] := Mod[Arg[czCond[czU[d, a, 2 Pi/Sqrt[d^2 + 8]]]]/Pi, 2];
czPhasePanel = Show[Plot[czPhase[d], {d, -6, 6}, PlotStyle -> Gray, PlotRange -> {0, 2}, PlotLegends -> Placed[LineLegend[{"2×2 block"}], {{0.03, 0.97}, {0, 1}}]],
   ListPlot[Table[{d, czPhase9[d, a]}, {a, {-3, -10, -40}}, {d, -6, 6, 0.5}], PlotMarkers -> {Automatic, Small}, PlotLegends -> Placed[LineLegend[Row[{"α/g = ", #}] & /@ {-3, -10, -40}, LegendLayout -> "Row"], Below]],
   Frame -> True, FrameLabel -> {"δ/g", "conditional phase / π"}, GridLines -> {{0}, {1}}, AspectRatio -> 1, ImagePadding -> {{50, 10}, {40, 10}}, ImageSize -> 300];
Grid[{{Legended[czBloch, Placed[LineLegend[ColorData[97] /@ Range[3], Row[{"δ/g = ", #}] & /@ czδs, LegendLayout -> "Row"], Below]], czPhasePanel}}, Alignment -> Top]
```

On resonance the turn passes through $|02\rangle$ and the phase is exactly $\pi$; detuning tilts the turn away from $|02\rangle$ and adds to the phase, so the controlled-Z needs $\delta=0$. The nine-level phases follow the closed form closely at $\alpha/g=-10$ and $-40$ and fall away from it at $\alpha/g=-3$.

After one full swap, the gate on all nine levels, with each transmon's own phases removed, is shown by its real and imaginary parts; the white square marks the four computational states:

```wl
czOrder = {{0, 0}, {0, 1}, {1, 0}, {1, 1}, {0, 2}, {2, 0}, {1, 2}, {2, 1}, {2, 2}};
czGate[d_, a_] := With[{u = czU[d, a, 2 Pi/Sqrt[d^2 + 8]][[1 + FromDigits[#, 3] & /@ czOrder, 1 + FromDigits[#, 3] & /@ czOrder]]},
   With[{φ0 = Arg[u[[1, 1]]], φ2 = Arg[u[[2, 2]]] - Arg[u[[1, 1]]], φ1 = Arg[u[[3, 3]]] - Arg[u[[1, 1]]]},
    DiagonalMatrix[Exp[-I (φ0 + φ1 #[[1]] + φ2 #[[2]])] & /@ czOrder] . u]];
czColor = ColorData["ThermometerColors"][Rescale[#, {-1, 1}]] &;
czTicks = Transpose[{Range[9], "|" <> StringJoin[ToString /@ #] <> "\[RightAngleBracket]" & /@ czOrder}];
czMatrix[m_, label_] := MatrixPlot[m, ColorFunction -> czColor, ColorFunctionScaling -> False, Mesh -> All, MeshStyle -> GrayLevel[0.5],
   FrameTicks -> {{czTicks, None}, {czTicks, None}}, PlotLabel -> label,
   Epilog -> {White, AbsoluteThickness[3], Line[{{0, 5}, {4, 5}, {4, 9}, {0, 9}, {0, 5}}]}, ImageSize -> 260];
czStep3 = AssociationMap[With[{u = czGate @@ #}, {czMatrix[Re[u], "Re U"], czMatrix[Im[u], "Im U"]}] &, Tuples[{{0, 1, 3}, {-3, -10, -40}}]];
Manipulate[Legended[Grid[{czStep3[{czδ, czα}]}], Placed[BarLegend[{czColor, {-1, 1}}], Right]],
  Row[{Column[{Control[{{czδ, 0, "δ/g"}, {0, 1, 3}, SetterBar}], Control[{{czα, -40, "α/g"}, {-3, -10, -40}, SetterBar}]}], Spacer[20], Dynamic[controlNotes[{
      {Row[{"δ/g = ", czδ}], <|0 -> "Resonance: |11⟩ picks up close to -1, nearer as |α|/g grows.",
        1 -> "|11⟩ picks up a complex phase, not -1.",
        3 -> "|11⟩ picks up a complex phase, further from -1."|>[czδ]},
      {Row[{"α/g = ", czα}], <|-3 -> "|01⟩ and |10⟩ mix; |11⟩ leaks into |20⟩.",
        -10 -> "Nearly isolated from the other levels.",
        -40 -> "The corner is cut off from the other levels."|>[czα]}}]]}], SaveDefinitions -> True]
```

On resonance and with $|\alpha|\gg g$, the corner is $\mathrm{diag}(1,1,1,-1)$, a controlled-Z, and nothing connects it to the other levels. Detuning keeps the corner diagonal but gives $|11\rangle$ a complex phase; a small $|\alpha|$ mixes $|01\rangle$ with $|10\rangle$ inside the corner and leaks $|11\rangle$ into $|20\rangle$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Pulse area sets the rotation


A qubit is driven for a time $T$ by a pulse $\Omega(t)=\Omega_0\,\varepsilon(t)$ of shape $\varepsilon$ with peak 1, detuned by $\delta=\omega_d-\omega_q$. Does only the pulse area set the rotation, for any shape and detuning?

With $\hbar=1$, ground state $|0\rangle$ and $\hat\sigma_z|0\rangle=|0\rangle$, the qubit in the frame rotating at $\omega_d$ has

$$
\hat H=\frac{\delta}{2}\,\hat\sigma_z+\frac{\Omega(t)}{2}\,\hat\sigma_x,\qquad \theta=\int_0^T\Omega\,dt,
$$

where $\varepsilon(t)$ is a Gaussian of width $T/5$ centred at $T/2$ and cut off at the pulse edges, a raised cosine $\sin^2(\pi t/T)$, or a square. On resonance $\hat H$ commutes with itself at all times and rotates the qubit by the area $\theta$, so $P(1)=\sin^2(\theta/2)$; the square pulse is solvable at any detuning, $P(1)=\frac{\theta^2}{\theta^2+\delta^2T^2}\sin^2\frac{\sqrt{\theta^2+\delta^2T^2}}{2}$.

From $|0\rangle$, the figure plots $P(1)$ at the end of the pulse against $\theta$ (points), with time in units of $T$ and the exact square-pulse formula dashed:

```wl
rabiShapes = <|"Gaussian" -> Exp[-25 (t - 1/2)^2/2], "Cosine" -> Sin[Pi t]^2, "Square" -> 1|>;
rabiArea = Integrate[#, {t, 0, 1}] & /@ rabiShapes;
rabiP1[shape_, δT_, θ_] := QuantumEvolve[δT/2 QuantumOperator["Z"] + N[θ/(2 rabiArea[shape])] rabiShapes[shape] QuantumOperator["X"],
     QuantumState["0"], {t, 0, 1}][1]["ProbabilitiesList"][[2]];
rabiPanel[δT_] := ListPlot[Prepend[Table[{θ, rabiP1[s, δT, θ]}, {s, Keys[rabiShapes]}, {θ, Subdivide[0, 2 Pi, 24]}],
     Table[{θ, θ^2/4 Sinc[Sqrt[θ^2 + δT^2]/2]^2}, {θ, Subdivide[0, 2 Pi, 120]}]],
   Joined -> {True, False, False, False}, PlotStyle -> {Dashed, Automatic, Automatic, Automatic}, PlotMarkers -> {"", "\[FilledCircle]", "\[FilledSquare]", "\[FilledUpTriangle]"},
   PlotLegends -> Prepend[Keys[rabiShapes], "square, exact"], PlotRange -> {{0, 2 Pi}, {0, 1}}, Frame -> True, FrameTicks -> {{Automatic, None}, {{#, Style[TraditionalForm[#], FractionBoxOptions -> {Beveled -> True}]} & /@ Range[0, 2 Pi, Pi/2], None}}, FrameLabel -> {"θ", "P(1)"},
   PlotLabel -> Row[{"T δ = ", δT}], ImageSize -> 300];
rabiPanels = rabiPanel /@ {0, 1, 3};
Legended[PlotGrid[{First /@ rabiPanels}, ImageSize -> 900], legendBelow[Last[First[rabiPanels]]]]
```

On resonance all three shapes fall on $\sin^2(\theta/2)$: only the area matters. Detuning breaks this, and the square pulse, which follows the exact formula, departs most: it spreads the area over the whole window, so its spectrum is the narrowest and a given detuning falls furthest down it.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Dispersive readout pointer states


A readout cavity of linewidth $\kappa$, coupled to a transmon, has its frequency pulled by $+\chi$ or $-\chi$ depending on the qubit state, and is driven by a readout pulse $\varepsilon_r(t)$ of peak amplitude $\epsilon_{\max}$. How does the cavity field separate for the two qubit states as $\kappa/2|\chi|$ changes?

The circuit and its Hamiltonian at three levels of approximation; hovering over a part of the circuit highlights its term, and Level 2 is the dispersive Hamiltonian used below:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/mbahram/TransmonLab/TransmonLab"]; Needs["TransmonLab`"];
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

### Where a DRAG pulse puts no power


A DRAG pulse has a Gaussian envelope $I(t)=e^{-(t-t_0)^2/2\sigma^2}$ ($\sigma=5$ ns, $t_0=20$ ns, 40 ns long) and a quadrature $Q(t)=-\beta\dot I(t)$, with $\beta$ to be calibrated by a sweep; the controller plays 1 ns samples. Where in frequency does each compiled pulse put no power?

With $\omega$ measured from the carrier and $\tilde f(\omega)=\int f(t)\,e^{i\omega t}dt$, the transform of $\dot I$ is $-i\omega\tilde I$, so the complex envelope has

$$
\widetilde{I+iQ}(\omega)=(1-\beta\omega)\,\tilde I(\omega),
$$

a spectral zero at $\omega=1/\beta$. With $\beta=\lambda/\alpha$, the DRAG pulse of the leakage slide ($\lambda=1$) puts it at $|\omega|=|\alpha|$, on the $|1\rangle\to|2\rangle$ transition; the sweep here, $|\beta|=\sigma$, puts it at $|\omega|=1/\sigma$.

For the sweep $\beta/\sigma\in\{-1,0,1\}$ the figure plots the closed-form power spectrum $|\widetilde{I+iQ}|^2$ (line) and that of each compiled pulse's samples (points), with $\omega$ in rad/ns:

```wl
PacletDirectoryLoad[{"/Users/mohammadb/Documents/GitHub/mbahram/TransmonLab/TransmonLab", "/Users/mohammadb/Documents/GitHub/mbahram/QUALink/QUALink"}];
Needs["WolframQuantumComputation`QUALink`"];
quaσ = 5; quaI = Exp[-(t - 20)^2/(2 quaσ^2)];
quaProg = QUAProgram["Qubit" -> {5.105*^9, {quaI, -β D[quaI, t]}, {t, 0, 40}}, "Resonator" -> 6.012*^9, "β" -> quaσ {-1, 0, 1}];
quaSpectrum = FullSimplify[Abs[FourierTransform[quaI - I β D[quaI, t], t, ω, FourierParameters -> {1, 1}]]^2, {β, ω} ∈ Reals];
quaSamples[pulse_] := quaProg["Waveforms"] /@ Values[quaProg["Config"]["pulses"][pulse]["waveforms"]];
quaSweep = quaProg["SweepValues"];
quaPoints = MapThread[With[{iq = quaSamples[#1]}, Table[{w, Abs[(iq[[1]] + I iq[[2]]) . Exp[I w Range[0, 39]]]^2}, {w, -0.6, 0.6, 0.04}]] &,
   {Select[Keys[quaProg["Config"]["pulses"]], StringStartsQ["qubit_drive"]], quaSweep}];
Show[LogPlot[Evaluate[quaSpectrum /. β -> # & /@ quaSweep], {ω, -0.6, 0.6}, PlotRange -> All, PlotPoints -> 100,
   PlotLegends -> Placed[LineLegend[Automatic, Row[{"β/σ = ", #/quaσ}] & /@ quaSweep, LegendLayout -> {"Row", 1}], Below]],
  ListLogPlot[quaPoints],
  PlotRange -> {Log[10^-3], Log[10^3]}, Axes -> {False, True}, Frame -> True, FrameLabel -> {"ω (rad/ns)", "power (ns²)"}, ImageSize -> 400]
```

The played samples (points) follow the closed-form spectrum (line) for every sweep value. At $\beta=0$ the spectrum is the Gaussian's, symmetric about the carrier; a nonzero $\beta$ tilts it, and the closed form has a zero at $\omega=1/\beta$, on the side set by the sign of $\beta$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The transmon spectrum against offset charge


A Cooper-pair box is a superconducting island, charging energy $E_C$, joined to a reservoir by a Josephson junction of energy $E_J$; a gate sets the offset charge $n_g$, in units of $2e$. How do its lowest levels depend on $n_g$ as $E_J/E_C$ grows?

With $\hat n$ the number of Cooper pairs on the island and $\hat\varphi$ the conjugate phase,

$$
\hat H=4E_C\,(\hat n-n_g)^2-E_J\cos\hat\varphi .
$$

In the phase basis this is Mathieu's equation, so the levels are $E_C\,a_{\nu}(-E_J/2E_C)$ with $\nu=2(n_g+k)$ over integers $k$, where $a_\nu$ is replaced by $b_{|\nu|}$ for $\nu<0$; the two differ only at integer and half-integer $n_g$, where $b$ supplies the odd partner level.

The figure shows the three lowest levels $E_m(n_g)$, measured from $E_0$ at $n_g=1/2$ in units of the qubit splitting $E_{01}=E_1-E_0$ there:

```wl
cpbLevels[ej_, ng_] := Take[Sort @ Table[With[{ν = N[2 (ng + k)]}, If[ν >= 0, MathieuCharacteristicA[ν, -ej/2.], MathieuCharacteristicB[-ν, -ej/2.]]], {k, -3, 3}], 3];
cpbPanel[ej_] := With[{ref = cpbLevels[ej, 1/2]},
   ListLinePlot[Transpose[Table[(cpbLevels[ej, ng] - ref[[1]])/(ref[[2]] - ref[[1]]), {ng, -1, 1, 1/100}]], DataRange -> {-1, 1}, PlotRange -> All,
    GridLines -> {None, {2}}, Frame -> True, PlotLegends -> (Row[{"m = ", #}] & /@ Range[0, 2]), ImageSize -> 260,
    FrameLabel -> {Subscript["n", "g"], TraditionalForm[HoldForm[Divide[Subscript["E", "m"][Subscript["n", "g"]] - Subscript["E", "0"][1/2], Subscript["E", "01"][1/2]]]]},
    PlotLabel -> Row[{Subscript["E", "J"], "/", Subscript["E", "C"], " = ", ej}]]];
cpbPanels = AssociationMap[cpbPanel, {1, 5, 10, 50}];
Manipulate[Legended[First[cpbPanels[ej]], legendBelow[Last[cpbPanels[1]]]], Row[{Control[{{ej, 10, Row[{Subscript["E", "J"], "/", Subscript["E", "C"]}]}, Keys[cpbPanels], SetterBar}], Spacer[20], Dynamic[controlNotes[{
      {Row[{Subscript["E", "J"], "/", Subscript["E", "C"], " = ", ej}], <|1 -> "Charge parabolas split at n_g = 1/2; the qubit frequency swings strongly.",
        5 -> "The bands flatten, but the qubit frequency still moves with n_g.",
        10 -> "The qubit band is nearly flat; the m = 2 band still moves.",
        50 -> "Flat bands; E₁₂ - E₀₁ is small, negative and fixed."|>[ej]}}]]}], SaveDefinitions -> True]
```

At $E_J/E_C=1$ the levels are charge parabolas that the junction splits at $n_g=1/2$, and the qubit frequency swings strongly with offset charge. As $E_J/E_C$ grows the bands flatten. The price is the anharmonicity: at $E_J/E_C=50$ the $m=2$ band sits flat just below the gridline at 2, so $E_{12}-E_{01}$ is small, negative and fixed.

How flat the bands become is measured by the charge dispersion of the qubit frequency, its full range over offset charge,

$$
\epsilon_{01}=\left|E_{01}(0)-E_{01}(\frac12)\right|,
$$

which for large $E_J/E_C$ tends to $2^5\sqrt{2/\pi}\left[(E_J/2E_C)^{3/4}+2^4(E_J/2E_C)^{5/4}\right]E_C\,e^{-\sqrt{8E_J/E_C}}$ (Koch et al. 2007), dashed:

```wl
cpbE01[ej_, ng_] := Subtract @@ Reverse @ Take[cpbLevels[ej, ng], 2];
Show[ListLogPlot[Table[{Sqrt[8 ej], Abs[cpbE01[ej, 0] - cpbE01[ej, 1/2]]}, {ej, 1, 80}], Joined -> True],
   ListLogPlot[Table[{Sqrt[8 ej], 2^5 Sqrt[2/Pi] ((ej/2)^(3/4) + 2^4 (ej/2)^(5/4)) Exp[-Sqrt[8 ej]]}, {ej, 1, 80}], Joined -> True, PlotStyle -> Dashed],
   PlotRange -> {{2.5, All}, All}, Frame -> True, FrameLabel -> {Sqrt[Row[{8, Subscript["E", "J"], "/", Subscript["E", "C"]}]], Row[{Subscript["\[Epsilon]", "01"], " / ", Subscript["E", "C"]}]}, ImageSize -> 300]
```

Once the bands are nearly flat, $\epsilon_{01}$ falls on a straight line on this logarithmic scale: it decreases exponentially in $\sqrt{8E_J/E_C}$, the exponent of phase tunneling between neighboring wells of the cosine potential, which is the only way the offset charge reaches the levels. The asymptote overshoots where the bands are still charge parabolas and closes in on the exact dispersion from above in the transmon regime, where offset charge barely moves the qubit frequency; its next-order corrections still keep it visibly above the exact curve at the right edge of the plot.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Interoperability: OpenQASM and real hardware

Circuits leave and enter QF as OpenQASM in pure Wolfram Language, can be transpiled to a device target, and run on IBM hardware. The cells that need an IBM account are marked not to evaluate.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### An imported CHSH circuit reaches the Tsirelson value exactly


A referee sends random bits $x$ to Alice and $y$ to Bob, who share an entangled pair and answer with bits $a$ and $b$; they win when $x\wedge y=a\oplus b$, which no classical strategy does more than $\frac34$ of the time. The strategy arrives as an OpenQASM 3 program in which Bob's measurement is set by a rotation angle $\theta_B$. How does the winning probability depend on $\theta_B$?

The program prepares $(|00\rangle+|11\rangle)/\sqrt2$ on Alice's qubit $q_0$ and Bob's $q_3$, and puts $q_1,q_2$ in $|{+}\rangle$ so that measuring them supplies $x$ and $y$. Alice applies a Hadamard to $q_0$ when $x=1$; Bob applies $R_Y(\theta_B)=U(\theta_B,0,0)$ to $q_3$, then a Hadamard when $y=0$; all four qubits are measured in $Z$, giving $(a,x,y,b)$.

For $\theta_B=k\pi/8$, $k=0,\dots,16$, the winning probability comes exactly from the joint outcome distribution of the program's circuit:

```wl
chshProgram = "OPENQASM 3.0;
qubit[4] q;
bit[4] c;
reset q[0];
reset q[3];
U(pi/2, 0, pi) q[3];
U(pi/2, 0, pi) q[0];
swap q[0], q[3];
ctrl(0) @ negctrl(1) @ U(pi/2, pi, pi) q[0], q[3];
ctrl(1) @ negctrl(0) @ U(pi/2, 0, 0) q[0], q[3];
reset q[1];
U(pi/2, 0, pi) q[1];
reset q[2];
U(pi/2, 0, pi) q[2];
ctrl(1) @ negctrl(0) @ U(pi/2, 0, pi) q[1], q[0];
U(THETAB, 0, 0) q[3];
ctrl(0) @ negctrl(1) @ U(pi/2, 0, pi) q[2], q[3];
c[0] = measure q[0];
c[1] = measure q[1];
c[2] = measure q[2];
c[3] = measure q[3];";
chshWin = Table[{k/8, FullSimplify @ Probability[BitAnd[xBit, yBit] == BitXor[aBit, bBit], {aBit, xBit, yBit, bBit} \[Distributed]
       QuantumCircuitOperator[StringReplace[chshProgram, "THETAB" -> ToString[k] <> "*pi/8"]][]["MultivariateDistribution"]]}, {k, 0, 16}];
chshBest = First[MaximalBy[chshWin, N @* Last]];
Column[{Labeled[QuantumCircuitOperator[StringReplace[chshProgram, "THETAB" -> "pi/4"]]["Diagram", ImageSize -> 600], Style["the imported OpenQASM program, at θ_B = π/4", 13, FontFamily -> "Helvetica"], Top],
  Show[ListPlot[chshWin, PlotMarkers -> {Automatic, 9}, PlotLegends -> {"imported program"}],
 Plot[{3/4, 1/2 + Cos[π u - π/4]/(2 Sqrt[2])}, {u, 0, 2}, PlotStyle -> {Dotted, Dashed},
  PlotLegends -> {"classical bound 3/4", TraditionalForm[HoldForm[1/2 + Cos[Subscript[θ, "B"] - π/4]/(2 Sqrt[2])]]}],
 PlotRange -> {{-0.05, 2.05}, {0, 1}}, Frame -> True, FrameLabel -> {Row[{Subscript["θ", "B"], "/π"}], "winning probability"},
 PlotLabel -> Style[Row[{"largest: ", TraditionalForm[Last[chshBest]], "  at  ", Subscript["θ", "B"], "/π = ", First[chshBest]}], 13, FontFamily -> "Helvetica"], ImageSize -> 460]}, Center, Spacings -> 2]
```

Every angle stays an exact multiple of $\pi$ through the import, so each probability is an exact number. The points lie on $\frac12+\cos(\theta_B-\pi/4)/(2\sqrt2)$, which beats the classical bound only for $0<\theta_B<\pi/2$ and peaks at the Tsirelson value $\cos^2\frac\pi8$ at $\theta_B=\pi/4$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A Toffoli written out for hardware: six CNOTs, seven T gates, and the same unitary

The Toffoli gate flips a target qubit when both of its control qubits are 1, $|a,b,c\rangle\mapsto|a,b,c\oplus ab\rangle$. A device executes only one- and two-qubit gates, so a program sent to it in OpenQASM 3, the text format hardware stacks read, must spell the Toffoli out in those. How many CNOTs and $T$ gates does the written program need, and does reading the program back give the Toffoli unitary again?

Every one-qubit line of the program is a rotation $U(\theta,\phi,\lambda)$, and the gates that occur are the Hadamard $H$, the phase gate $T$ and its inverse $T^\dagger$:

$$
U(\theta,\phi,\lambda)=\begin{pmatrix}\cos\frac{\theta}{2}&-e^{i\lambda}\sin\frac{\theta}{2}\\ e^{i\phi}\sin\frac{\theta}{2}&e^{i(\phi+\lambda)}\cos\frac{\theta}{2}\end{pmatrix},\qquad H=U\left(\frac{\pi}{2},0,\pi\right),\qquad T=U\left(0,0,\frac{\pi}{4}\right),\qquad T^\dagger=U\left(0,0,\frac{7\pi}{4}\right).
$$

The Toffoli with controls $q_0,q_1$ and target $q_2$ is exported, each one-qubit line is named $H$, $T$ or $T^\dagger$ from its angles, and the program is read back and its unitary compared entry by entry with the Toffoli matrix:

```wl
toffoliProgram = QuantumQASM[QuantumCircuitOperator["Toffoli"]];
toffoliBack = QuantumQASM[toffoliProgram];
toffoliGateName = <|{1/2, 1} -> "H", {0, 1/4} -> "T", {0, 7/4} -> Superscript["T", "†"]|>;
toffoliCounts = KeyTake[Counts[StringCases[toffoliProgram, {"cx " -> "CNOT",
     "U(" ~~ θ : NumberString ~~ ", 0., " ~~ λ : NumberString ~~ ")" :> toffoliGateName[Rationalize[ToExpression[{θ, λ}]/π, 10^-12]]}]], {"CNOT", "T", Superscript["T", "†"], "H"}];
Grid[{{Framed[Style[toffoliProgram, 11, FontFamily -> "Courier"], FrameStyle -> GrayLevel[0.7], FrameMargins -> 10],
   Column[{Grid[{Keys[toffoliCounts], Values[toffoliCounts]}, Frame -> All, Spacings -> {2, 1}, BaseStyle -> {14, FontFamily -> "Helvetica"}],
     Style[Row[{"largest entry of |", Subscript["U", "read back"], " \[Minus] ", Subscript["U", "Toffoli"], "|:  ",
        ScientificForm[Max[Abs[Normal[toffoliBack["Matrix"]] - Normal[QuantumCircuitOperator["Toffoli"]["Matrix"]]]], 2]}], 13, FontFamily -> "Helvetica"]}, Center, Spacings -> 3]}},
  Spacings -> 4, Alignment -> Center]
```

The program uses six CNOTs and seven $T$-type gates, four $T$ and three $T^\dagger$, with two Hadamards on the target. Six CNOTs is the fewest any Toffoli circuit can use (Shende and Markov), and seven is the fewest $T$ gates without extra qubits (Gosset, Kliuchnikov, Mosca and Russo). The angles are written as 16-digit decimals, so the program read back reproduces the Toffoli matrix to rounding error, not exactly.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What a line of qubits costs the quantum Fourier transform

The quantum Fourier transform $F$ on $n$ qubits couples every pair of qubits through a controlled phase. A device offers a few native gates, here $\mathrm{CZ}$, $R_z$, $\sqrt X$ and $X$, and couples only some pairs; on a line, only neighbours. How many two-qubit gates does $F$ need on a line compared with full connectivity, and is the compiled circuit still $F$?

$F$ is a sequence of Hadamards and the $n(n-1)/2$ controlled phases $\mathrm{diag}(1,1,1,e^{2\pi i/2^{k}})$, $k\geq2$, ending in a reversal of the qubit order:

$$
F|x\rangle=\frac{1}{2^{n/2}}\sum_{y=0}^{2^n-1}e^{2\pi i\,xy/2^n}|y\rangle .
$$

Each controlled phase takes two $\mathrm{CZ}$ gates, and the reversal costs nothing if the outputs are relabeled, so full connectivity needs $n(n-1)$. On a line, following each controlled phase with a swap of the same pair carries every qubit past every other once and ends in the reversed order; a controlled phase together with a swap takes three $\mathrm{CZ}$ gates, so $\frac{3}{2}n(n-1)$ suffice.

For $n=2,\dots,8$ the circuit is compiled to the native gates on both couplings (best of five compiler seeds, input qubits held in place) and its $\mathrm{CZ}$ gates are counted; each compiled unitary $U$ must satisfy $UF^\dagger=P$ with $P$ a permutation of the qubits. The dashed and dotted curves are the two counts above:

```wl
qftCoupling = <|"all-to-all" -> (Select[Tuples[Range[0, # - 1], 2], Apply[Unequal]] &), "line" -> (Catenate[{{#, # + 1}, {# + 1, #}} & /@ Range[0, # - 2]] &)|>;
qftNative[n_, map_, seed_] := QuantumQASM[QuantumCircuitOperator["Fourier"[n]], "Target" -> <|"NumQubits" -> n, "BasisGates" -> {"cz", "rz", "sx", "x"},
     "CouplingMap" -> qftCoupling[map][n]|>, "OptimizationLevel" -> 3, "SeedTranspiler" -> seed, "InitialLayout" -> Range[0, n - 1]];
qftBest = Association @ Table[{map, n} -> First[MinimalBy[qftNative[n, map, #] & /@ Range[0, 4], StringCount[#, "cz "] &]], {map, Keys[qftCoupling]}, {n, 2, 8}];
qftOffPermutation[prog_, n_] := With[{m = Normal[QuantumQASM[prog]["Matrix"]] . ConjugateTranspose[Normal[QuantumCircuitOperator["Fourier"[n]]["Matrix"]]]},
   With[{p = Round[Abs[m]]}, Max[Abs[m - Tr[Transpose[p] . m]/2^n p]]]];
Show[ListPlot[Table[{n, StringCount[qftBest[{map, n}], "cz "]}, {map, Keys[qftCoupling]}, {n, 2, 8}], PlotMarkers -> {Automatic, 10},
   PlotLegends -> {"all-to-all", "line, nearest neighbours"}],
  Plot[{n (n - 1), 3 n (n - 1)/2}, {n, 2, 8}, PlotStyle -> {Dashed, Dotted},
   PlotLegends -> {TraditionalForm[HoldForm[n (n - 1)]], TraditionalForm[HoldForm[3 n (n - 1)/2]]}],
  Frame -> True, FrameLabel -> {"qubits n", "two-qubit gates (cz)"}, PlotRange -> {0, All}, ImageSize -> 480,
  PlotLabel -> Style[Row[{"largest deviation of ", Subscript["U", "transpiled"], Superscript["F", "\[Dagger]"], " from a qubit permutation: ",
      ScientificForm[Max[KeyValueMap[qftOffPermutation[#2, Last[#1]] &, qftBest]], 2]}], 13, FontFamily -> "Helvetica"]]
```

With every pair coupled the count is exactly $n(n-1)$: two $\mathrm{CZ}$ gates per controlled phase and none for the reversal. On the line the compiled count stays at or below $\frac{3}{2}n(n-1)$ up to $n=5$ and runs above it from $n=6$ on, so the compiler's heuristic routing misses the swap network. Every compiled circuit equals $F$ followed by a relabeling of qubits, to rounding error.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A GHZ state on an IBM processor: where the lost shots go

A Hadamard on qubit 1 and CNOTs from qubit 1 to 2 and from 2 to 3 prepare the three-qubit GHZ state. Measured in $Z$, it gives $000$ or $111$, each with probability $\frac{1}{2}$, and nothing else. On a real processor some shots land on the six other outcomes. How many, and is the loss spread evenly over the three qubits?

$$
|\mathrm{GHZ}\rangle=\mathrm{CNOT}_{2\to3}\,\mathrm{CNOT}_{1\to2}\,H_1|000\rangle=\frac{|000\rangle+|111\rangle}{\sqrt2},\qquad P(000)=P(111)=\frac{1}{2}.
$$

In each of the six other outcomes exactly one qubit disagrees with the other two, so the lost shots sort into three groups, one per qubit. The circuit, ending in a measurement of all three qubits, is sent to IBM's ibm_fez processor for 4096 shots; this needs an IBM Quantum account, with its API key and instance CRN in place of "xxx":

```wl
#| eval: false
LocalSymbol["ibm_crn"] = "xxx";
ServiceConnect["IBMQuantumPlatform", "New", Authentication -> {"apikey" -> "xxx"}];
job = IBMJobSubmit[QuantumCircuitOperator[{"GHZ"[3], {1, 2, 3}}], "ibm_fez", "Shots" -> 4096, "Wait" -> True];
job["Refresh"]["Counts"]
```

The counts recorded from one such run (job dah09eb9k43c73ae1a4g, submitted 10 September 2026) are set against the exact distribution, and the lost shots are grouped by the qubit that disagrees:

```wl
ghzCounts = <|{0, 0, 0} -> 1854, {0, 0, 1} -> 22, {0, 1, 0} -> 31, {0, 1, 1} -> 129, {1, 0, 0} -> 105, {1, 0, 1} -> 42, {1, 1, 0} -> 33, {1, 1, 1} -> 1880|>;
ghzShots = Total[ghzCounts];
ghzExact = ghzShots KeyMap[First, QuantumCircuitOperator[{"GHZ"[3], {1, 2, 3}}][]["Probabilities"]];
ghzOddQubit[b_] := If[Equal @@ b, 0, First[FirstPosition[b, If[Total[b] == 1, 1, 0]]]];
ghzByQubit[c_] := KeySort[GroupBy[Normal[c], ghzOddQubit @* First -> Last, Total]];
Grid[{{
   BarChart[Transpose[{Lookup[ghzExact, Keys[ghzCounts]], Values[ghzCounts]}], ChartLabels -> {Row /@ Keys[ghzCounts], None},
    ChartLegends -> {"exact", Row[{"ibm_fez, ", ghzShots, " shots"}]}, Frame -> True, FrameLabel -> {"outcome", "shots"}, ImageSize -> 470,
    PlotLabel -> Style[Row[{"shots on 000 or 111: ", ghzByQubit[ghzCounts][0], " of ", ghzShots}], 13, FontFamily -> "Helvetica"]],
   BarChart[Transpose[Values[KeyDrop[ghzByQubit[#], 0]] & /@ {ghzExact, ghzCounts}], ChartLabels -> {Row[{"qubit ", #}] & /@ Range[3], None},
    LabelingFunction -> Above, Frame -> True, FrameLabel -> {"the qubit that disagrees with the other two", "shots"}, ImageSize -> 330, PlotRange -> {0, All}]}},
 Spacings -> 3, Alignment -> Bottom]
```

About nine in ten shots stay on $000$ and $111$, split nearly evenly between the two, as in the exact state. The lost shots are not spread evenly: qubit 1 disagrees with the other two three to four times as often as qubit 2 or qubit 3.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

## Algebra beyond quantum states

The same language handles noncommutative algebras, Clifford algebras and finite-element eigenproblems, so the physics around the quantum model is computed in one place.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Snyder's Lorentz-covariant quantum spacetime


Snyder made the spacetime coordinates $x^\mu$ noncommuting without breaking Lorentz symmetry: their commutator is a Lorentz generator $M^{\mu\nu}$ times the square of a fundamental length $a$. Is there a quadratic combination of the ten generators that commutes with all of them?

With $\eta=\operatorname{diag}(1,-1,-1,-1)$, the algebra is

$$
[x^\mu,x^\nu]=ia^2M^{\mu\nu},\qquad [M^{\mu\nu},x^\rho]=i\left(\eta^{\nu\rho}x^\mu-\eta^{\mu\rho}x^\nu\right),\qquad [M^{\mu\nu},M^{\rho\sigma}]=i\left(\eta^{\nu\rho}M^{\mu\sigma}-\eta^{\nu\sigma}M^{\mu\rho}-\eta^{\mu\rho}M^{\nu\sigma}+\eta^{\mu\sigma}M^{\nu\rho}\right),
$$

and the ansatz is $x^2+c\,\frac12M_{\mu\nu}M^{\mu\nu}$, with $x^2=\eta_{\mu\nu}x^\mu x^\nu$ and the Lorentz Casimir $\frac12M_{\mu\nu}M^{\mu\nu}$.

With commutators reduced to an ordered basis by the relations, the commutator with $x^1$ fixes $c$, and each table entry says whether the reduced commutator with a generator is zero:

```wl
snyderEta[μ_, ν_] := DiagonalMatrix[{1, -1, -1, -1}][[μ + 1, ν + 1]];
lorentzM[μ_, ν_] := Which[μ == ν, 0, μ < ν, snyderM[μ, ν], True, -snyderM[ν, μ]];
lorentzPairs = Subsets[Range[0, 3], {2}];
snyderGenerators = Join[snyderX /@ Range[0, 3], snyderM @@@ lorentzPairs];
snyderAlgebra = NonCommutativeAlgebra[<|"Generators" -> {snyderGenerators}, "ScalarVariables" -> {snyderA, snyderC},
    "CommutationRelations" -> Join[
      Commutator[snyderX[#2], snyderX[#1]] + I snyderA^2 snyderM[#1, #2] & @@@ lorentzPairs,
      Flatten[Table[Commutator[snyderM[μ, ν], snyderX[ρ]] - I (snyderEta[ν, ρ] snyderX[μ] - snyderEta[μ, ρ] snyderX[ν]), {μ, 0, 3}, {ν, μ + 1, 3}, {ρ, 0, 3}]],
      Replace[Subsets[lorentzPairs, {2}], {{μ_, ν_}, {ρ_, σ_}} :> Commutator[snyderM[ρ, σ], snyderM[μ, ν]] +
         I (snyderEta[ν, ρ] lorentzM[μ, σ] - snyderEta[ν, σ] lorentzM[μ, ρ] - snyderEta[μ, ρ] lorentzM[ν, σ] + snyderEta[μ, σ] lorentzM[ν, ρ]), {1}]]|>];
xSquared = Sum[snyderEta[μ, μ] snyderX[μ] ** snyderX[μ], {μ, 0, 3}];
lorentzCasimir = Sum[snyderEta[μ, μ] snyderEta[ν, ν] snyderM[μ, ν] ** snyderM[μ, ν], {μ, 0, 3}, {ν, μ + 1, 3}];
snyderCoef = First[SolveValues[Coefficient[NonCommutativeExpand[Commutator[xSquared + snyderC lorentzCasimir, snyderX[1]], snyderAlgebra], snyderX[1]] == 0, snyderC]];
snyderCasimirLabel = Row[{"½ ", Subscript["M", "μν"], Superscript["M", "μν"]}];
Style[Labeled[Grid[Prepend[MapThread[Prepend, {Outer[If[NonCommutativeExpand[Commutator[#1, #2], snyderAlgebra] === 0, 0, "\[NotEqual] 0"] &,
       {xSquared, lorentzCasimir, xSquared + snyderCoef lorentzCasimir}, snyderGenerators],
      {Superscript["x", 2], snyderCasimirLabel, Row[{Superscript["x", 2], " + (", snyderCoef /. snyderA -> "a", ") ", snyderCasimirLabel}]}}],
    Prepend[Replace[snyderGenerators, {snyderX[μ_] :> Superscript["x", μ], snyderM[μ_, ν_] :> Superscript["M", Row[{μ, ν}]]}, {1}], "commutator with"]],
   Frame -> All, Spacings -> {1.5, 1}], Row[{"[", Superscript["x", 2], " + c ", snyderCasimirLabel, ", ", Superscript["x", 1], "] = 0  \[DoubleLongRightArrow]  c = ", snyderCoef /. snyderA -> "a"}], Top],
  16, FontFamily -> "Helvetica"]
```

$x^2$ and the Lorentz Casimir each commute with every $M^{\mu\nu}$, as Lorentz scalars must, but neither commutes with the coordinates. Commuting with one coordinate fixes $c=-a^2$, and that combination then commutes with all ten generators: Snyder's spacetime keeps a central quadratic element, the Minkowski $x^2$ corrected at order $a^2$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Adding two momenta on a noncommuting spacetime


A plane wave $e^{i\alpha_\mu x^\mu}$ carries momentum $\alpha_\mu$, and when the coordinates $x^\mu$ do not commute, the product of two plane waves defines a deformed sum $\alpha\oplus\beta$. Snyder's spacetime has a fundamental length $a$, and $\kappa$-Minkowski spacetime has a mass scale $\kappa$. Does the deformed sum depend on which momentum comes first?

With $\eta=\operatorname{diag}(1,-1,-1,-1)$ and Lorentz generators $M^{\mu\nu}$, the two algebras are

$$
\text{Snyder:}\ \ [x^\mu,x^\nu]=ia^2M^{\mu\nu},\quad [M^{\mu\nu},x^\rho]=i\left(\eta^{\nu\rho}x^\mu-\eta^{\mu\rho}x^\nu\right);\qquad \kappa\text{-Minkowski:}\ \ [x^0,x^j]=\frac{i}{\kappa}\,x^j,\quad [x^j,x^k]=0,
$$

with the Lorentz algebra for $[M,M]$. Writing $e^{X}e^{Y}=e^{Z}$, the Baker-Campbell-Hausdorff series to fourth order is

$$
Z=X+Y+\frac12[X,Y]+\frac1{12}\left([X,[X,Y]]+[Y,[Y,X]]\right)-\frac1{24}[Y,[X,[X,Y]]],\qquad X=i\alpha_\mu x^\mu,\ Y=i\beta_\mu x^\mu.
$$

In $\kappa$-Minkowski spacetime $Z$ is a combination of the $x^\mu$ alone, and $(\alpha\oplus\beta)_\mu$ is the coefficient of $ix^\mu$ in $Z$, the momentum of the symmetrically ordered plane wave. The matrices $x^0=\frac{i}{\kappa}E_{11}$, $x^1=E_{12}$ obey the relations, so $e^Xe^Y$ can be multiplied and its logarithm taken in closed form. Writing the same plane wave with $x^0$ to the right, $e^{ik_jx^j}e^{ik_0x^0}$, relabels the momentum as $k_0=\gamma_0$, $k_j=\gamma_j\,\kappa\left(1-e^{-\gamma_0/\kappa}\right)/\gamma_0$, and turns the sum into the bicrossproduct law $k\oplus q=\left(k_0+q_0,\ k_j+e^{-k_0/\kappa}q_j\right)$. Back in the symmetric labels the exact sum is $(\alpha\oplus\beta)_j=\left(\alpha_je^{\beta_0/\kappa}\varphi(\alpha_0/\kappa)+\beta_j\varphi(\beta_0/\kappa)\right)/\varphi\left((\alpha_0+\beta_0)/\kappa\right)$ with $\varphi(x)=(e^x-1)/x$: each spatial momentum is weighted by a factor that depends only on the energies.

In Snyder's algebra $Z$ also carries Lorentz generators, starting with $\frac12[X,Y]=-\frac{i}{2}a^2\alpha_\mu\beta_\nu M^{\mu\nu}$, and they cannot be dropped: $M^{\mu\nu}$ does not commute with $x^\rho$, so $e^{i\gamma_\mu x^\mu+i\omega_{\mu\nu}M^{\mu\nu}}$ is not a plane wave of momentum $\gamma$. Instead the product is split as $e^Xe^Y=e^{i(\alpha\oplus\beta)_\mu x^\mu}e^{i\omega_{\mu\nu}M^{\mu\nu}}$; the Lorentz factor leaves the Lorentz-invariant vacuum unchanged, so $\alpha\oplus\beta$ is the momentum of the product. Comparing the two Baker-Campbell-Hausdorff series, $(\alpha\oplus\beta)_\mu$ is the coefficient of $ix^\mu$ in

$$
X+Y-\frac16[X,[X,Y]]+\frac13[Y,[Y,X]],
$$

and since Snyder's relations are unchanged under $x\to-x$, $M\to M$, the next correction is fifth order. The $5\times5$ matrices of the de Sitter algebra, $x^\mu=aL^{\mu4}$ and $M^{\mu\nu}=L^{\mu\nu}$ with metric $\operatorname{diag}(1,-1,-1,-1,-1)$, obey Snyder's relations and give the exact sum: $e^Xe^Y$ carries the point $(0,0,0,0,1)$ of de Sitter space to a new point, and $\alpha\oplus\beta$ is the $\gamma$ for which $e^{i\gamma_\mu x^\mu}$ carries it to the same place. Each spatial $x^j$ rotates that point in a compact plane, so these coordinates label points uniquely only while $a\sqrt{-\gamma\cdot\gamma}<\pi$; beyond that, different $\gamma$ carry $(0,0,0,0,1)$ to the same point. The momenta below stay well inside that range.

Each set of relations reduces the series for general $\alpha$, $\beta$; the figure shows two momenta in the $(p_0,p_1)$ plane, in units of $1/a$ for Snyder and of $\kappa$ for $\kappa$-Minkowski, for $\alpha=s\,(1,\frac12,0,0)$ and $\beta=s\,(\frac35,-\frac45,0,0)$ at three scales $s$, with the exact sums from the $5\times5$ and $2\times2$ matrices:

```wl
composeEta = DiagonalMatrix[{1, -1, -1, -1}]; composeG = DiagonalMatrix[{1, -1, -1, -1, -1}];
composeL[A_, A_] = 0; composeL[A_, B_] := If[A < B, composeJ[A, B], -composeJ[B, A]];
composePairs = Subsets[Range[0, 4], {2}];
{composeAl, composeBe} = {composeα /@ Range[0, 3], composeβ /@ Range[0, 3]};
composeRelations["Snyder"] = Replace[Subsets[composePairs, {2}], {{A_, B_}, {C_, D_}} :> Commutator[composeJ[C, D], composeJ[A, B]] +
     I (composeG[[B + 1, C + 1]] composeL[A, D] - composeG[[B + 1, D + 1]] composeL[A, C] - composeG[[A + 1, C + 1]] composeL[B, D] + composeG[[A + 1, D + 1]] composeL[B, C]), {1}];
composeRelations["κ-Minkowski"] = Join[Commutator[composeT[#], composeT[0]] + I composeT[#]/composeK & /@ Range[3], Commutator[composeT[#2], composeT[#1]] & @@@ Subsets[Range[3], {2}]];
composeAlgebra["Snyder"] = NonCommutativeAlgebra[<|"Generators" -> {composeJ @@@ composePairs}, "ScalarVariables" -> Join[{composeA}, composeAl, composeBe], "CommutationRelations" -> composeRelations["Snyder"]|>];
composeAlgebra["κ-Minkowski"] = NonCommutativeAlgebra[<|"Generators" -> {composeT /@ Range[0, 3]}, "ScalarVariables" -> Join[{composeK}, composeAl, composeBe], "CommutationRelations" -> composeRelations["κ-Minkowski"]|>];
composeFrame["Snyder"] = {composeA, composeJ[#, 4] & /@ Range[0, 3]}; composeFrame["κ-Minkowski"] = {1, composeT /@ Range[0, 3]};
composeBCH["Snyder"][x_, y_] := x + y - Commutator[x, Commutator[x, y]]/6 + Commutator[y, Commutator[y, x]]/3;
composeBCH["κ-Minkowski"][x_, y_] := x + y + Commutator[x, y]/2 + (Commutator[x, Commutator[x, y]] + Commutator[y, Commutator[y, x]])/12 - Commutator[y, Commutator[x, Commutator[x, y]]]/24;
composeSeriesSum[name_] := composeSeriesSum[name] = With[{c = First[composeFrame[name]], e = Last[composeFrame[name]]},
    Simplify[Coefficient[NonCommutativeExpand[composeBCH[name][I c composeAl . e, I c composeBe . e], composeAlgebra[name]], #]/(I c) & /@ e]];
composeSwap = {composeα -> composeβ, composeβ -> composeα};
composeDot[u_, v_] := u . composeEta . v;
composeLmat[A_, B_] := I (composeG[[A + 1, A + 1]] SparseArray[{A + 1, B + 1} -> 1, {5, 5}] - composeG[[B + 1, B + 1]] SparseArray[{B + 1, A + 1} -> 1, {5, 5}]);
composeSinhOverTheta[z_] := Hypergeometric0F1[3/2, z/4]; composeCoshTerm[z_] := Hypergeometric0F1[3/2, z/16]^2/2;
composeBoost[α_] := With[{v = I α . (composeLmat[#, 4] & /@ Range[0, 3])}, IdentityMatrix[5] + composeSinhOverTheta[composeDot[α, α]] v + composeCoshTerm[composeDot[α, α]] v . v];
composeThetaOverSinh[c_] := Hypergeometric2F1[1, 1, 3/2, (1 - c)/2];
composeExactSum["Snyder"][α_, β_] := With[{p = composeBoost[composeA α] . composeBoost[composeA β] . UnitVector[5, 5]}, composeThetaOverSinh[p[[5]]] {-p[[1]], p[[2]], p[[3]], p[[4]]}/composeA];
composePhi[x_] := Exp[x/2] Hypergeometric0F1[3/2, x^2/16];
composeExactSum["κ-Minkowski"][α_, β_] := Prepend[(Rest[α] Exp[β[[1]]/composeK] composePhi[α[[1]]/composeK] + Rest[β] composePhi[β[[1]]/composeK])/composePhi[(α[[1]] + β[[1]])/composeK], α[[1]] + β[[1]]];
composeKappaRep = {{{I/composeK, 0}, {0, 0}}, {{0, 1}, {0, 0}}};
composeChecks = {
   "5×5 matrices obey Snyder's relations" -> AllTrue[Flatten[Normal[composeRelations["Snyder"] /. {p_ ** q_ :> (p /. composeJ -> composeLmat) . (q /. composeJ -> composeLmat), composeJ -> composeLmat}]], # === 0 &],
   "2×2 matrices: e^X e^Y = exp(i(α ⊕ β)·x), exact κ sum" -> (Simplify[MatrixExp[I composeAl[[;; 2]] . composeKappaRep] . MatrixExp[I composeBe[[;; 2]] . composeKappaRep] -
        MatrixExp[I composeExactSum["κ-Minkowski"][composeAl, composeBe][[;; 2]] . composeKappaRep]] === {{0, 0}, {0, 0}}),
   "exact − series = O(5th order), Snyder" -> (Simplify[Normal[Series[composeExactSum["Snyder"][t composeAl, t composeBe] - (composeSeriesSum["Snyder"] /. {composeα[μ_] :> t composeα[μ], composeβ[μ_] :> t composeβ[μ]}), {t, 0, 4}]]] === {0, 0, 0, 0}),
   "exact − series = O(1/κ⁴), κ-Minkowski" -> (Simplify[Normal[Series[(composeExactSum["κ-Minkowski"][composeAl, composeBe] - composeSeriesSum["κ-Minkowski"]) /. composeK -> 1/ε, {ε, 0, 3}]]] === {0, 0, 0, 0}),
   "Snyder, (0, 5/2, 0, 0) ⊕ (0, 3/2, 0, 0) = (0, 4 − 2π, 0, 0) in units 1/a" -> (Simplify[composeExactSum["Snyder"][{0, 5/2, 0, 0}, {0, 3/2, 0, 0}] /. composeA -> 1] === {0, 4 - 2 Pi, 0, 0})};
composeSnyderSplit = With[{u = Array[composeU, 3], v = Array[composeV, 3], s = {composeDot[composeAl, composeAl], composeDot[composeBe, composeBe], composeDot[composeAl, composeBe]}},
   {u, v} /. First[SolveAlways[Thread[composeSeriesSum["Snyder"] - (composeSeriesSum["Snyder"] /. composeSwap) == (u . s) composeAl + (v . s) composeBe], Join[composeAl, composeBe]]]];
composeAssoc[name_, unit_, a_, b_, c_] := With[{f = Chop[N[composeExactSum[name][#1, #2] /. unit]] &}, Chop[f[f[a, b], c] - f[a, f[b, c]]]];
composeTriple = 3/5 {{1, 1/2, 0, 0}, {3/5, -4/5, 0, 0}, {1/2, 0, 1/2, 0}};
composeShow = {composeα[μ_] :> Subscript["α", μ], composeβ[μ_] :> Subscript["β", μ], composeK -> "κ", composeA -> "a"};
composePanel[name_, unit_, s_] := With[{pts = Thread[Join[composeAl, composeBe] -> s Join[{1, 1/2, 0, 0}, {3/5, -4/5, 0, 0}]]},
   With[{α = composeAl /. pts, β = composeBe /. pts},
    ListPlot[List /@ Chop[N[{α + β, composeSeriesSum[name] /. pts, composeSeriesSum[name] /. composeSwap /. pts, composeExactSum[name][α, β], composeExactSum[name][β, α]} /. unit]][[All, ;; 2]],
     PlotMarkers -> {{"\[FilledSquare]", 14}, {"\[EmptyCircle]", 16}, {"\[EmptyCircle]", 30}, {"+", 14}, {"\[Times]", 14}},
     PlotLegends -> {"α + β", "α ⊕ β, fourth order", "β ⊕ α, fourth order", "α ⊕ β, exact", "β ⊕ α, exact"},
     Epilog -> {Gray, Arrow[{{0, 0}, α[[;; 2]]}], Text["α", 0.55 α[[;; 2]], {0, -1.5}], Dashed, Arrow[{{0, 0}, β[[;; 2]]}], Text["β", 0.55 β[[;; 2]], {-1.5, 0}]},
     PlotRange -> s {{-0.1, 1.9}, {-1.1, 0.7}}, AspectRatio -> Automatic, Frame -> True, FrameLabel -> {Subscript["p", 0], Subscript["p", 1]},
     PlotLabel -> Row[{name, ",   s = ", s}], ImageSize -> 260]]];
composePanels = Table[composePanel[Sequence @@ m, s], {m, {{"Snyder", composeA -> 1}, {"κ-Minkowski", composeK -> 1}}}, {s, {0.3, 0.6, 0.9}}];
Labeled[Legended[Grid[Map[First, composePanels, {2}]], legendBelow[Last[composePanels[[2, 1]]]]],
  Style[Column[{
     Row[{"Snyder, fourth order:  α ⊕ β − β ⊕ α = (", TraditionalForm[Factor[composeSnyderSplit[[1]] . {Superscript["α", 2], Superscript["β", 2], "α·β"}] /. composeShow], ") α + (",
       TraditionalForm[Factor[composeSnyderSplit[[2]] . {Superscript["α", 2], Superscript["β", 2], "α·β"}] /. composeShow], ") β"}],
     Row[{"κ-Minkowski, fourth order:  α ⊕ β − β ⊕ α = ", TraditionalForm[Factor[composeSeriesSum["κ-Minkowski"] - (composeSeriesSum["κ-Minkowski"] /. composeSwap)] /. composeShow]}],
     Row[{"grouping, exact, α, β, γ = ", Row[composeTriple, ",  "], " in units 1/a or κ:  (α ⊕ β) ⊕ γ − α ⊕ (β ⊕ γ) = ", NumberForm[composeAssoc["Snyder", composeA -> 1, Sequence @@ composeTriple], 2], " (Snyder),  ", composeAssoc["κ-Minkowski", composeK -> 1, Sequence @@ composeTriple], " (κ-Minkowski)"}],
     If[AllTrue[Values[composeChecks], TrueQ], "checks passed: the 5×5 matrices obey Snyder's relations, the 2×2 product gives the exact κ sum, both exact sums match their series, and a spatial Snyder momentum wraps around",
      Row[{"checks failed: ", Row[Keys[Select[Association[composeChecks], Not @* TrueQ]], ", "]}]]}], 13, FontFamily -> "Helvetica"], Top]
```

The deformed sum moves away from $\alpha+\beta$ in both spacetimes, more so at larger $s$, and in both it depends on which momentum comes first. In Snyder's spacetime the two orders split by $-\frac{a^2}{2}\left(\alpha^2\beta-\beta^2\alpha+(\alpha\cdot\beta)(\beta-\alpha)\right)$, cubic in the momenta and zero when they are parallel: the order is recorded in the Lorentz factor that the product leaves behind. In $\kappa$-Minkowski spacetime the split is already quadratic, a spatial term proportional to $\alpha_0\beta_j-\alpha_j\beta_0$. The exact $\kappa$-Minkowski sums sit on the fourth-order points; the exact Snyder sums begin to separate from theirs at $s=0.6$ and clearly leave them at the largest $s$, where the momenta approach $1/a$ and the fifth-order term is no longer small. Grouping separates the two: adding three momenta in Snyder's spacetime depends on which pair is added first, because its momentum space is de Sitter space, a coset rather than a group, while in $\kappa$-Minkowski spacetime the momenta form a group and the grouping does not matter.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Dirac matrices for a curved metric


At one point of a curved spacetime the metric $g_{\mu\nu}$ is a real symmetric matrix, in general not diagonal, with $p$ positive and $n-p$ negative eigenvalues. Does the vielbein construction give gamma matrices obeying the Clifford relation and the flat-space contraction identities for every such metric?

The gamma matrices must satisfy the Clifford relation, and a vielbein $e$ reduces $g$ to its signature $\eta=\operatorname{diag}(\pm1)$:

$$
\{\gamma_\mu,\gamma_\nu\}=2g_{\mu\nu}\mathbb 1,\qquad g=e^{T}\eta\,e,\qquad \gamma_\mu=e^{a}{}_{\mu}\,\Gamma_a ,
$$

where the flat $\Gamma_a$ are $2^{\lfloor n/2\rfloor}$-dimensional tensor products of Pauli matrices (the Brauer-Weyl construction), multiplied by $i$ where $\eta_{aa}=-1$. Indices are raised with $g^{\mu\nu}$. Any real symmetric nondegenerate $g$ has a real vielbein (Sylvester's law of inertia), and then

$$
\{\gamma_\mu,\gamma_\nu\}=e^a{}_\mu e^b{}_\nu\{\Gamma_a,\Gamma_b\}=2e^a{}_\mu e^b{}_\nu\eta_{ab}\mathbb 1=2g_{\mu\nu}\mathbb 1,
$$

while $\gamma^\mu\gamma_\nu\gamma_\mu=e^b{}_\nu\,\Gamma^a\Gamma_b\Gamma_a=(2-n)\gamma_\nu$ carries over from flat space.

The construction is tested on one random non-diagonal machine-number metric per $(n,p)$, with eigenvalue magnitudes between $0.1$ and $10$, for $n=2,\dots,8$ and $p=0,\dots,n$, and for one exact integer metric. On the left, the anticommutators for one random non-diagonal Lorentzian metric are drawn as a $4\times4$ grid of $4\times4$ blocks beside $2g_{\mu\nu}$, on one color scale; on the right, each cell gives the error of the Clifford relation, relative to the largest metric entry, for one signature $(n,p)$:

```wl
Needs["WolframQuantumComputation`DiracMatrix`", "/Users/mohammadb/Nextcloud/Mads/Quantum/DiracMatrix/DiracMatrixPackage/Kernel/DiracMatrix.wl"];
gammaErrors[n_, p_] := With[{g = RandomCurvedMetric[p, n - p, 1, 1]}, With[{γ = GammaMatrices[g]},
   {Max[Abs[Outer[#1 . #2 + #2 . #1 &, γ, γ, 1] - 2 TensorProduct[g, IdentityMatrix[Length[First[γ]]]]]]/Max[Abs[g]],
    Max[Max[Abs[Total[MapThread[Function[{up, dn}, up . # . dn], {Inverse[g] . γ, γ}]] - (2 - n) #]]/Max[Abs[#]] & /@ γ]}]];
gammaResults = BlockRandom[SeedRandom[1]; Table[gammaErrors[n, p], {n, 2, 8}, {p, 0, n}]];
gammaExactMetric = {{1, 2, 0}, {2, -3, 1}, {0, 1, 2}};
gammaLorentz = BlockRandom[SeedRandom[3]; RandomCurvedMetric[1, 3, 1, 1]];
gammaShade = With[{m = 2 Max[Abs[gammaLorentz]]}, ColorData["RedBlueTones"][Rescale[#, {m, -m}]] &];
gammaBlocks = Grid[{{ArrayPlot[Re[ArrayFlatten[Outer[#1 . #2 + #2 . #1 &, #, #, 1]]] &[GammaMatrices[gammaLorentz]], ColorFunction -> gammaShade, ColorFunctionScaling -> False,
     Mesh -> {{4, 8, 12}, {4, 8, 12}}, MeshStyle -> Gray, ImageSize -> 260, PlotLabel -> Row[{"{", Subscript["γ", "μ"], ", ", Subscript["γ", "ν"], "}: 4 × 4 blocks of 4 × 4"}]],
    ArrayPlot[2 gammaLorentz, ColorFunction -> gammaShade, ColorFunctionScaling -> False, Mesh -> All, MeshStyle -> Gray, ImageSize -> 150, PlotLabel -> Row[{"2", Subscript["g", "μν"]}]]}}, Spacings -> 2];
gammaMap = Legended[Graphics[Flatten[MapIndexed[With[{n = First[#2] + 1, p = Last[#2] - 1, e = Log10[First[#1]]},
       {ColorData["TemperatureMap"][Rescale[e, {-17, -14}]], Rectangle[{p - 0.5, n - 0.5}, {p + 0.5, n + 0.5}], Black, Text[Style[NumberForm[e, {3, 1}], 10], {p, n}]}] &, gammaResults, {2}], 2],
     Frame -> True, FrameLabel -> {"positive eigenvalues p", "dimension n"}, FrameTicks -> {{Range[2, 8], None}, {Range[0, 8], None}}, ImageSize -> 340,
     PlotLabel -> Row[{Subscript["log", 10], " of the Clifford error, every signature"}]], BarLegend[{"TemperatureMap", {-17, -14}}]];
Labeled[Grid[{{gammaBlocks, gammaMap}}, Spacings -> 4],
  Style[Column[{Row[{"random Lorentzian metric (left):  every block of {", Subscript["γ", "μ"], ", ", Subscript["γ", "ν"], "} is 2", Subscript["g", "μν"], " times the identity"}],
     Row[{Superscript["γ", "μ"], Subscript["γ", "ν"], Subscript["γ", "μ"], " = (2 \[Minus] n)", Subscript["γ", "ν"], ":  largest relative error over all signatures ", ScientificForm[Max[gammaResults[[All, All, 2]]], 1]}],
     Row[{"exact metric ", MatrixForm[gammaExactMetric], ":  Clifford residual = ",
       Max[Abs[RootReduce[Outer[#1 . #2 + #2 . #1 &, #, #, 1] &[GammaMatrices[gammaExactMetric]] - 2 TensorProduct[gammaExactMetric, IdentityMatrix[2]]]]]}]}, Center], 13, FontFamily -> "Helvetica"], Top]
```

Both relations hold to rounding error in every dimension and signature tried, and exactly for the integer metric. The identity above already guarantees them for every metric, so the figure tests that the vielbein and the flat $\Gamma_a$ are built correctly, including for indefinite, non-diagonal metrics.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A particle bound by its own gravity


A particle of mass $m$ moves in the Newtonian potential $\Phi$ sourced by its own mass density $m|\psi|^2$: quantum spreading pushes the wave out, self-gravity pulls it in. What is the energy of its stationary bound state once the finite computational box of radius $R$ is removed?

The Schrödinger-Newton equations are

$$
i\hbar\,\partial_t\psi=-\frac{\hbar^2}{2m}\nabla^2\psi+m\Phi\,\psi,\qquad \nabla^2\Phi=4\pi Gm\,|\psi|^2 .
$$

With $\hbar=m=G=1$, lengths in units of $\hbar^2/Gm^3$, a spherically symmetric stationary state $\psi=e^{-iEt}\,u(r)/(\sqrt{4\pi}\,r)$ obeys

$$
\Phi''+\frac{2}{r}\Phi'=\frac{u^2}{r^2},\qquad -\frac12u''+\Phi\,u=E\,u,\qquad u(0)=u(R)=\Phi(R)=0,\quad \int_0^R u^2\,dr=1 .
$$

Far from the particle $\Phi\approx-1/r$, so setting $\Phi(R)=0$ raises the energy by about $1/R$. Without a box, the ground state found by Moroz, Penrose and Tod (1998) follows by shooting: $S=u/r$ and $W=\Phi-E$ obey $S''+2S'/r=2WS$, $W''+2W'/r=S^2$ with $S(0)=1$, $W(0)$ is bisected until $S$ decays without a node, outside the particle $W=W(\infty)-\mathcal N/r$ with $\mathcal N=\int S^2r^2dr=r^2W'$ (Gauss's law), so $W(\infty)=W+rW'$, and rescaling to unit norm gives $E=-W(\infty)/\mathcal N^2$.

Finite elements on a mesh of 1000 elements graded toward $r=0$ solve Newton's equation for $\Phi$ given $u$, then the lowest eigenpair given $\Phi$, repeated from a normalized Gaussian until $E$ changes by less than $10^{-7}$, in four boxes; the shooting energy is the reference:

```wl
Needs["NDSolve`FEM`"];
selfPotential[uR_, mesh_] := NDSolveValue[{Laplacian[snPhi[snR], {snR, θ, φ}, "Spherical"] == uR^2/snR^2,
    DirichletCondition[snPhi[snR] == 0, snR == mesh["Bounds"][[1, 2]]]}, snPhi, {snR} ∈ mesh];
solitonGroundState[ϕ_, mesh_] := First /@ NDEigensystem[{SchrodingerPDEComponent[{snU[snR], {snR}}, <|"ReducedPlanckConstant" -> 1, "SchrodingerPotential" -> ϕ[snR]|>],
     DirichletCondition[snU[snR] == 0, True]}, snU[snR], {snR} ∈ mesh, 1, Method -> {"Eigensystem" -> {"Arnoldi", "Shift" -> -1}}];
solitonEnergy[boxR_] := With[{mesh = ToGradedMesh[{Line[{{0}, {boxR}}], <|"Alignment" -> "Left", "ElementCount" -> 1000|>}]},
   Re @ First @ NestWhile[solitonGroundState[selfPotential[Last[#], mesh], mesh] &,
      solitonGroundState[selfPotential[With[{snU0 = snR Exp[-snR^2/9]}, snU0/Sqrt[Integrate[snU0^2, {snR, 0, Infinity}]]], mesh], mesh], Abs[First[#1] - First[#2]] > 10^-7 &, 2, 60]];
boxEnergies = {#, solitonEnergy[#]} & /@ {20, 40, 80, 160};
snEquations[w_] := {snS''[snR] + 2 snS'[snR]/snR == 2 snW[snR] snS[snR], snW''[snR] + 2 snW'[snR]/snR == snS[snR]^2, snS[10^-6] == 1, snS'[10^-6] == 0, snW[10^-6] == w, snW'[10^-6] == 0};
snNodeQ[w_?NumericQ] := Reap[NDSolve[Join[snEquations[w], {WhenEvent[snS[snR] < 0, Sow[True]; "StopIntegration"], WhenEvent[snS'[snR] > 0, Sow[False]; "StopIntegration"]}],
      {snS, snW}, {snR, 10^-6, 60}, WorkingPrecision -> 30]][[2, 1, 1]];
snW0 = First[Nest[With[{mid = Mean[#]}, If[snNodeQ[mid], {mid, Last[#]}, {First[#], mid}]] &, N[{-2, 0}, 30], 45]];
snShooting = With[{sol = First[NDSolve[snEquations[snW0], {snS, snW}, {snR, 10^-6, 12}, WorkingPrecision -> 30]]}, -(snW[12] + 12 snW'[12])/(144 snW'[12])^2 /. sol];
ListPlot[{boxEnergies, {#1, #2 - 1/#1} & @@@ boxEnergies, {{0, snShooting}, {170, snShooting}}}, Joined -> True,
  PlotStyle -> {Automatic, Automatic, Directive[Thick, Dotted]}, PlotMarkers -> {{Automatic, 12}, {Automatic, 12}, None},
  PlotLegends -> {"E(R)", "E(R) \[Minus] 1/R", "free space, by shooting"}, PlotRange -> {{0, 170}, All}, Frame -> True, ImageSize -> 420,
  FrameLabel -> {Row[{"box radius R (", Superscript["ℏ", 2], "/G", Superscript["m", 3], ")"}], Row[{"energy E (", Superscript["G", 2], Superscript["m", 5], "/", Superscript["ℏ", 2], ")"}]}]
```

A bound state exists in every box, and its energy falls toward the free-space value as the box grows. Subtracting $1/R$ puts every box on the shooting energy (dotted): the box shifts the energy by the cut-off $-1/r$ tail, and what remains is below the resolution of the plot.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Stored numbers, declared bases and composition


An operator can be specified by the numbers it takes in a chosen basis: store the numbers $C=\operatorname{diag}(1,-1)$ together with that basis as the operator $R$. What does applying $Z$ after $R$ give?

If the columns of $F$ are the declared basis vectors, the stored matrix $C$ represents the map $A=F\,C\,F^{-1}$ in the computational basis; $F$ only needs to be invertible. For an orthonormal basis $A$ is the observable with outcomes $\pm1$ in that basis; for a basis that is not orthonormal, such as $\{|0\rangle,|+\rangle\}$, $A$ is not Hermitian, a linear map rather than an observable. Composition contracts the shared wire after both operators are written in one basis, so $Z\circ R$ is the product $ZA$. Contracting the stored numbers as they sit would give $ZC$ instead.

For $R$ in four declared bases, the table lists the map $A$ it represents, the composite $Z\circ R$, and whether that composite equals $ZA$, after $C$ and $ZC$:

```wl
semNames = {"{|0⟩, |1⟩}", "{|+⟩, |−⟩}", "{|+i⟩, |−i⟩}", "{|0⟩, |+⟩}"};
semBases = {"Computational", "PauliX", "PauliY", QuditBasis[{{1, 0}, {1, 1}/Sqrt[2]}]};
semC = PauliMatrix[3];
semRow[b_] := With[{r = QuantumOperator[semC, b]}, With[{am = Simplify[Normal[r["MatrixRepresentation"]]], zam = Simplify[Normal[QuantumOperator["Z"][r]["MatrixRepresentation"]]]},
    {TraditionalForm[MatrixForm[am]], TraditionalForm[MatrixForm[zam]], Simplify[zam - semC . am] === ConstantArray[0, {2, 2}]}]];
Style[Labeled[Grid[Prepend[MapThread[Prepend[semRow[#1], #2] &, {semBases, semNames}],
    {"declared basis", "represented A", "Z after R in QF", "= Z A"}], Frame -> All, Spacings -> {2, 1}],
  Row[{"C = ", TraditionalForm[MatrixForm[semC]], ",   Z C = ", TraditionalForm[MatrixForm[semC . semC]], " in every basis"}], Top], 13, FontFamily -> "Helvetica"]
```

The same numbers $C$ represent $Z$, $X$, $Y$, or, in the non-orthogonal basis, a map that is not Hermitian, and the composite is $Z$ times that map in every row. Multiplying the stored numbers gives the identity whatever the basis, which is the right answer only when the declared basis is the computational one.

