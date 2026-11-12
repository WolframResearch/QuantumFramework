```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A Schrödinger-cat state in phase space

A field mode holds a superposition of two coherent states of real amplitude $\pm\alpha$ with relative phase $\varphi$. What sets the sign of its Wigner function at the origin?

With $x=(a+a^\dagger)/\sqrt2$ and $p=(a-a^\dagger)/i\sqrt2$, the state and its Wigner function are

$$
|\psi\rangle=\mathcal N\left(|\alpha\rangle+e^{i\varphi}|{-\alpha}\rangle\right),\qquad W(x,p)=\frac1\pi\int\langle x-y|\rho|x+y\rangle\,e^{2ipy}\,dy,
$$

and the photon-number parity $\Pi=(-1)^{a^\dagger a}$ maps $|\alpha\rangle$ to $|{-\alpha}\rangle$.

QF computes $\pi W$ from the state on a Fock space cut at 40 levels. Each panel prints $\pi W(0,0)$ and $\langle\Pi\rangle$, its line marks $W=0$, and all panels share one color scale:

```wl
catColors = Blend["WL14DefaultDensityGradient", (1 + #)/2] &;
catPanel[α_, φ_] := With[{st = CatState[40][N[α], N[φ]]}, With[{w = WignerRepresentation[st, {-6, 6}, {-3, 3}, "GridSize" -> 101],
    parity = Re @ Total[(-1)^Range[0, 39] Diagonal[Normal[st["DensityMatrix"]]]]},
  DensityPlot[Pi w[x, p], {x, -6, 6}, {p, -3, 3}, PlotRange -> {-1, 1}, ColorFunction -> catColors, ColorFunctionScaling -> False, PlotPoints -> 70,
   MeshFunctions -> {#3 &}, Mesh -> {{0}}, FrameLabel -> {"x", "p"}, AspectRatio -> 1/2, ImageSize -> 320,
   PlotLabel -> Column[{Row[{"\[Alpha] = ", α, ",   \[CurlyPhi] = ", φ}],
     Row[{"\[Pi]W(0,0) = ", Round[Pi w[0, 0], 0.001], ",   \[LeftAngleBracket]\[CapitalPi]\[RightAngleBracket] = ", Round[parity, 0.001]}]}, Alignment -> Center]]]];
Legended[Grid[Outer[catPanel, {1, 2, 3}, {0, Pi/2, Pi}]], BarLegend[{catColors, {-1, 1}}, LegendLabel -> "\[Pi]W"]]
```

From $\alpha=2$ on, the two coherent lobes separate at $x=\pm\sqrt2\,\alpha$, and between them fringes along $p$ dip below zero, closer together as $\alpha$ grows. The value at the origin is the parity in every panel: $\varphi=0$ and $\varphi=\pi$ give parity eigenstates, so the origin sits at the extreme $\pm1$; at $\varphi=\pi/2$ the parity falls toward zero and the origin lands between a positive and a negative fringe.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Squeezed thermal light: Wigner against Husimi

A field mode in thermal equilibrium, with mean photon number $\bar n$, is squeezed by a degenerate parametric amplifier with squeeze parameter $r$. Which phase-space picture shows its quadrature noise directly, and how does each compare with the vacuum?

With $x=(a+a^\dagger)/\sqrt2$ and $p=(a-a^\dagger)/i\sqrt2$, the state is

$$
\rho=S(r)\,\rho_{\rm th}\,S^\dagger(r),\qquad S(r)=e^{r\left(a^2-a^{\dagger2}\right)/2},\qquad \rho_{\rm th}=\sum_n\frac{\bar n^{\,n}}{(1+\bar n)^{n+1}}\,|n\rangle\langle n| .
$$

Its Wigner function is $W(x,p)=\frac1\pi\int\langle x-y|\rho|x+y\rangle\,e^{2ipy}\,dy$, and its Husimi function $Q(x,p)=\langle\alpha|\rho|\alpha\rangle/2\pi$, with $\alpha=(x+ip)/\sqrt2$, is $W$ smoothed by the vacuum's Wigner function $e^{-x^2-p^2}/\pi$.

QF computes both on a Fock space cut at 70 levels. Each panel draws the contours at half the peak of $W$ and of $Q$, and, dashed, those of the vacuum's $W$ and $Q$:

```wl
sqThermal[nb_, r_] := SqueezeOperator[N[r], 70][ThermalState[N[nb], 70]];
sqPanel[nb_, r_] := With[{st = sqThermal[nb, r]}, With[{w = WignerRepresentation[st, {-5, 5}, {-5, 5}, "GridSize" -> 101],
    q = HusimiQRepresentation[st, {-5, 5}, {-5, 5}, "GridSize" -> 101]},
  ContourPlot[{w[x, p] == w[0, 0]/2, q[x, p] == q[0, 0]/2, x^2 + p^2 == Log[2], x^2 + p^2 == 2 Log[2]}, {x, -5, 5}, {p, -5, 5},
   ContourStyle -> {Automatic, Automatic, Dashed, Dashed}, PlotLegends -> {"Wigner W", "Husimi Q", "vacuum W", "vacuum Q"},
   FrameLabel -> {"x", "p"}, PlotLabel -> Row[{"\!\(\*OverscriptBox[\(n\), \(_\)]\) = ", nb, ",   r = ", r}], ImageSize -> 260]]];
sqPanels = Outer[sqPanel, {0, 0.5, 1}, {0, 0.5, 1}];
Legended[Grid[Map[First, sqPanels, {2}]], Last[sqPanels[[1, 1]]]]
```

Squeezing narrows both ellipses along $x$ and stretches them along $p$, while thermal noise widens them in every direction, so at $\bar n=1$ the Wigner ellipse at $r=0.5$ still sits just outside the vacuum's circle and only $r=1$ pushes it inside. $Q$ is $W$ smoothed by the vacuum, so every $Q$ ellipse is its $W$ ellipse widened by one vacuum's worth of noise: it never enters the vacuum's $W$ circle, yet it enters the vacuum's own $Q$ circle exactly when $W$ enters the vacuum's $W$ circle. Only the width of $W$ is the quadrature noise itself.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### One photon or a coherent pulse in a Mach-Zehnder interferometer

A single photon, or a coherent pulse of amplitude $\alpha$, enters port 2 of a Mach-Zehnder interferometer whose arms differ in phase by $\varphi$ and whose two beam splitters have mixing angle $\theta$. Detectors $D_1$ and $D_2$ count the photons leaving the two output ports. What does a single photon do that a classical pulse cannot?

The interferometer and the two inputs are

$$
U=B(\theta)\,e^{i\varphi\,a_1^\dagger a_1}\,B(\theta),\qquad B(\theta)=e^{\theta\left(a_1a_2^\dagger-a_1^\dagger a_2\right)},\qquad |0\rangle_1|1\rangle_2\ \ \text{or}\ \ |0\rangle_1|\alpha\rangle_2 .
$$

QF applies $U$ as a Fock-space operator with $\varphi$ symbolic, on 16 levels per mode, and plots the share of the photons that reaches $D_1$ and the probability that both detectors fire:

```wl
mzLevels = 16;
mzInputs = <|"|1\[RightAngleBracket]" -> {FockState[1, mzLevels], 1}, "|\[Alpha] = 1\[RightAngleBracket]" -> {CoherentState[mzLevels][1.], 1},
    "|\[Alpha] = 2\[RightAngleBracket]" -> {CoherentState[mzLevels][2.], 4}|>;
mzCounts[θ_, in_] := mzCounts[θ, in] = With[{bs = BeamSplitterOperator[{N[θ], 0}, mzLevels]},
    ArrayReshape[bs[PhaseShiftOperator[φ, mzLevels][bs[QuantumTensorProduct[FockState[0, mzLevels], in]]]]["ProbabilitiesList"], {mzLevels, mzLevels}]];
mzShare[θ_, {in_, n_}] := Range[0, mzLevels - 1] . Total[mzCounts[θ, in], {2}]/n;
mzBoth[θ_, {in_, _}] := With[{pr = mzCounts[θ, in]}, 1 - Total[pr[[1]]] - Total[pr[[All, 1]]] + pr[[1, 1]]];
mzPanel[obs_, label_][θ_] := Plot[Evaluate[obs[θ, #] & /@ Values[mzInputs]], {φ, 0, 2 Pi}, PlotRange -> {-0.02, 1.02}, FrameTicks -> {{Automatic, None}, {Range[0, 2 Pi, Pi/2], None}},
    PlotStyle -> {Thick, Dashed, Dashed}, PlotLegends -> Keys[mzInputs], Frame -> True, FrameLabel -> {"\[CurlyPhi]", label},
    PlotLabel -> Row[{"\[Theta] = ", θ}], ImageSize -> 300];
mzPanels = Outer[#1[#2] &, {mzPanel[mzShare, "share at \!\(\*SubscriptBox[\(D\), \(1\)]\)"], mzPanel[mzBoth, "P(both fire)"]}, {Pi/4, Pi/8}];
Legended[Grid[Map[First, mzPanels, {2}]], Last[mzPanels[[1, 1]]]]
```

The share reaching $D_1$ follows one fringe for the photon and for both pulses: the mean intensity does not see the statistics of the light, and unbalanced splitters shrink the fringe for all three alike. The single photon never fires both detectors, at any phase or splitting. A coherent pulse does, more often as $\alpha$ grows, except at the phases where the fringe sends all of its light to one port.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Photon statistics behind a beam splitter

Light from a source holding $\langle n\rangle$ photons on average (a photon-number state $|n\rangle$, a laser in a coherent state, or a thermal lamp) falls on a beam splitter of mixing angle $\theta$, with vacuum in the other port, and a detector counts the photons in each output. Do the joint counts tell the three sources apart, whatever the splitting?

The splitter mixes the source mode $a$ with the empty mode $b$, and $N_c,N_d$ count the photons in its two outputs. The splitter, the source's second-order coherence and the detectors' normalized joint count are

$$
B(\theta)=e^{\theta\left(a\,b^\dagger-a^\dagger b\right)},\qquad g^{(2)}(0)=\frac{\langle a^{\dagger2}a^2\rangle}{\langle a^\dagger a\rangle^2},\qquad \frac{\langle N_cN_d\rangle}{\langle N_c\rangle\langle N_d\rangle} .
$$

QF computes $g^{(2)}(0)$ of each source with `G2Coherence`, and the joint-count ratio from the two-mode photon-number distribution after the splitter, on 44 levels per mode:

```wl
hbtLevels = 44;
hbtSources[n_] := <|"Fock |n\[RightAngleBracket]" -> FockState[n, hbtLevels], "coherent" -> CoherentState[hbtLevels][N[Sqrt[n]]],
    "thermal" -> ThermalState[N[n], hbtLevels]|>;
hbtSplitter[θ_] := hbtSplitter[θ] = BeamSplitterOperator[{N[θ], 0}, hbtLevels];
hbtRatio[θ_][src_] := With[{k = Range[0, hbtLevels - 1],
    pr = ArrayReshape[hbtSplitter[θ][QuantumTensorProduct[src, FockState[0, hbtLevels]]]["ProbabilitiesList"], {hbtLevels, hbtLevels}]},
   (k . pr . k)/((k . Total[pr, {2}]) (Total[pr] . k))];
hbtPanel[f_, title_, ylabel_] := ListPlot[Transpose[Table[{n, f[#]} & /@ Values[hbtSources[n]], {n, 4}]], PlotMarkers -> Automatic,
    PlotRange -> {{0.8, 4.2}, {-0.05, 2.2}}, Frame -> True, FrameTicks -> {{Automatic, None}, {Range[4], None}},
    FrameLabel -> {"\[LeftAngleBracket]n\[RightAngleBracket]", ylabel}, PlotLabel -> title, PlotLegends -> Keys[hbtSources[1]], ImageSize -> 280];
hbtJoint = "\[LeftAngleBracket]\!\(\*SubscriptBox[\(N\), \(c\)]\)\!\(\*SubscriptBox[\(N\), \(d\)]\)\[RightAngleBracket]/\[LeftAngleBracket]\!\(\*SubscriptBox[\(N\), \(c\)]\)\[RightAngleBracket]\[LeftAngleBracket]\!\(\*SubscriptBox[\(N\), \(d\)]\)\[RightAngleBracket]";
hbtPanels = {hbtPanel[G2Coherence, "G2Coherence of the source", "\!\(\*SuperscriptBox[\(g\), \((2)\)]\)(0)"],
    hbtPanel[hbtRatio[Pi/4], "joint counts, \[Theta] = \[Pi]/4", hbtJoint], hbtPanel[hbtRatio[Pi/8], "joint counts, \[Theta] = \[Pi]/8", hbtJoint]};
Legended[Row[First /@ hbtPanels, Spacer[10]], Last[First[hbtPanels]]]
```

The three panels agree: the normalized joint count equals the source's $g^{(2)}(0)$ at either splitting. The photon-number state sits below 1, at $1-1/n$, the laser at 1 and the thermal lamp at 2, so the joint counts sort the sources by their statistics while the mean photon number cannot.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What one boson-sampling amplitude costs

$n$ identical photons enter a random lossless network of $m=2n$ modes with occupation list $S=\{s_1,\dots,s_m\}$, and we ask for the amplitude of one output list $T$. How does the cost of that one number depend on $n$ and on how the photons share the modes?

With $W$ the $m\times m$ transfer matrix and $W_{T,S}$ the $n\times n$ matrix that repeats column $k$ $s_k$ times and row $j$ $t_j$ times,

$$
A(T\leftarrow S)=\frac{\operatorname{Per}W_{T,S}}{\sqrt{\prod_k s_k!\,\prod_j t_j!}} .
$$

Ryser's and Glynn's formulas evaluate a permanent in about $n\,2^n$ operations. When modes hold several photons, QF can instead group them into a sum of $\min\{\prod_k(s_k+1),\prod_j(t_j+1)\}$ terms, and it takes whichever of the two routes its cost estimate favors.

The cell times each amplitude as the fastest of three runs, each after `ClearSystemCache[]`, with the input in the first $n$ modes and the output pattern moved to the last $n$:

```wl
bsOccupations[n_] := PadRight[#, 2 n] & /@ <|"collision-free" -> ConstantArray[1, n],
    "pairs" -> Append[ConstantArray[2, Quotient[n, 2]], Mod[n, 2]], "bunched" -> {n}|>;
bsNetwork[n_] := BlockRandom[SeedRandom[1107 + 1000 n]; QuantumInterferometer["Random"[2 n]]];
bsSeconds[q_, s_] := With[{t = RotateRight[s, Total[s]]}, Min @ Table[ClearSystemCache[]; First @ AbsoluteTiming[q["Amplitude", s -> t]], 3]];
bsTimes = AssociationMap[Function[n, With[{q = bsNetwork[n]}, bsSeconds[q, #] & /@ bsOccupations[n]]], {8, 12, 16, 20, 24}];
ListLogPlot[Values @ Merge[KeyValueMap[Function[{n, r}, {n, #} & /@ r], bsTimes], Identity], Joined -> True, PlotMarkers -> Automatic,
  Frame -> True, FrameLabel -> {"photons n  (m = 2n modes)", "seconds per amplitude"}, PlotLegends -> Keys[bsOccupations[4]], ImageSize -> 450]
```

The collision-free amplitude grows exponentially, a straight line on the log scale. For pairs QF's cost estimate still prefers the direct permanent, so the paired and collision-free curves coincide. The fully bunched amplitude needs only $n+1$ terms and stays essentially free. The cost of one amplitude is set by how the photons share the modes, not by the photon number alone.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Cascaded beam splitters with symbolic photon numbers

Three modes are chained by beam splitters: the first mixes modes 1 and 2, the second modes 2 and 3, the third modes 1 and 2 again, splitter $i$ with mixing angle $\theta_i$ and phase $\phi_i$. With $n_1$ and $n_2$ photons entering modes 1 and 2, and every photon number a symbol, how likely is the input pattern $|n,n,0\rangle$ to come back out, for every $n$ at once?

The network is

$$
U=B_{12}(\theta_3,\phi_3)\,B_{23}(\theta_2,\phi_2)\,B_{12}(\theta_1,\phi_1),\qquad B_{jk}(\theta,\phi)=e^{\theta\left(e^{i\phi}a_ja_k^\dagger-e^{-i\phi}a_j^\dagger a_k\right)} .
$$

For this return the phase $\phi_2$ cancels, and $\phi_1,\phi_3$ enter only through $\phi_1-\phi_3$. QF returns $\langle t_1,t_2,t_3|U|n_1,n_2,0\rangle$ with all photon numbers, angles and phases symbolic (top label, with summation index $K[1]$); the panels evaluate it at $\theta_1=\pi/4$ and $\phi_1=\phi_3$:

```wl
cascade = QuantumInterferometer[{"BS"[θ1, ϕ1] -> {1, 2}, "BS"[θ2, ϕ2] -> {2, 3}, "BS"[θ3, ϕ3] -> {1, 2}}, 3];
cascadeAmp = cascade["Amplitude", {n1, n2, 0} -> {t1, t2, t3}];
cascadeReturn[a2_, a3_][n_] := Abs[Activate[cascadeAmp /. {n1 -> n, n2 -> n, t1 -> n, t2 -> n, t3 -> 0,
     θ1 -> N[Pi/4], θ2 -> N[a2], θ3 -> N[a3], ϕ1 -> 0, ϕ2 -> 0, ϕ3 -> 0}]]^2;
cascadeAngles = {-Pi/4, Pi/8, Pi/3};
cascadePanel[a2_] := ListPlot[Table[{n, cascadeReturn[a2, a3][n]}, {a3, cascadeAngles}, {n, 12}], Joined -> True, PlotMarkers -> Automatic,
    PlotRange -> {-0.02, 1.02}, Frame -> True, FrameLabel -> {"n", "P(|n,n,0\[RightAngleBracket] \[Rule] |n,n,0\[RightAngleBracket])"},
    PlotLabel -> Row[{"\!\(\*SubscriptBox[\(\[Theta]\), \(2\)]\) = ", a2}], PlotLegends -> (Row[{"\!\(\*SubscriptBox[\(\[Theta]\), \(3\)]\) = ", #}] & /@ cascadeAngles), ImageSize -> 300];
cascadePanels = cascadePanel /@ {0, Pi/8, Pi/4};
Labeled[Legended[Row[First /@ cascadePanels, Spacer[10]], Last[First[cascadePanels]]], TraditionalForm[cascadeAmp], Top]
```

The amplitude is one finite sum over $K$, the photon number in mode 1 between the first and the last splitter; no detector sees $K$, so its values interfere, and $K\le\min(n_1+n_2,\,t_1+t_2)$ because those photons start and end in modes 1 and 2. One formula gives the return for every $n$: without a leak into mode 3, $\theta_3=-\pi/4$ undoes the first splitter and the pattern always returns, while any leak makes the return fall with $n$.
