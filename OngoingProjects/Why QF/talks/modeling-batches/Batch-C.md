```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Measurement purifies a qubit: Itô's lemma in spherical coordinates

A qubit with Hamiltonian $H=\frac{\Delta}{2}\sigma_z+\frac{\Omega}{2}\sigma_x$ is measured continuously along $\sigma_z$ at strength $k$ by a detector of efficiency $\eta$. What does the length $r$ of its Bloch vector do along one run?

With $\hbar=1$, the state conditioned on the record follows the Itô stochastic master equation of a $\sigma_z$ measurement: $H$ rotates it, the term $k\,\mathcal D[\sigma_z]\rho$ dephases it, and the record kicks it through $\sqrt{\eta k}\left(\sigma_z\rho+\rho\sigma_z-2\langle\sigma_z\rangle\rho\right)dW$.

WL writes this as an Itô process for the Bloch vector, applies Itô's lemma to $r$ in spherical coordinates $(r,\theta,\varphi)$, with $z=r\cos\theta$ (top), and samples three runs from $(x,y,z)=(0.3,0,0.4)$; time is in units of $1/k$:

```wl
blochSME[Δ_, Ω_, k_, η_, v0_] := With[{ρv = (IdentityMatrix[2] + {x[t], y[t], z[t]} . (PauliMatrix /@ {1, 2, 3}))/2},
   With[{bloch = Function[m, Simplify[Tr[m . PauliMatrix[#]] & /@ {1, 2, 3}]], h = Δ/2 PauliMatrix[3] + Ω/2 PauliMatrix[1], sz = PauliMatrix[3]},
    ItoProcess[{bloch[-I (h . ρv - ρv . h) + k (sz . ρv . sz - ρv)],
      List /@ bloch[Sqrt[η k] (sz . ρv + ρv . sz - 2 z[t] ρv)]}, {{x, y, z}, v0}, t]]];
drSph = FullSimplify[blochSME[Δ, Ω, k, η, {x0, y0, z0}]["Derivative", Sqrt[x[t]^2 + y[t]^2 + z[t]^2]] /.
    Thread[{x[t], y[t], z[t]} -> r {Sin[θ] Cos[φ], Sin[θ] Sin[φ], Cos[θ]}], r > 0 && 0 < θ < π && k > 0 && 0 <= η <= 1 && {Δ, Ω} ∈ Reals];
panelPurify[η_, {Ωv_, Δv_}] := ListLinePlot[Map[{#[[1]], Norm[#[[2]]]} &,
     RandomFunction[blochSME[Δv, Ωv, 1, η, {0.3, 0, 0.4}], {0., 4., 0.001}, 3, Method -> "StochasticRungeKuttaScalarNoise"]["Paths"], {2}],
    PlotRange -> {0, 1.05}, Frame -> True, FrameLabel -> {"k t", "r"}, PlotLabel -> Row[{"η = ", η, ",   Ω/k = ", Ωv, ",   Δ/k = ", Δv}], ImageSize -> 300];
SeedRandom[2]; Labeled[Grid[Outer[panelPurify, {0.5, 1}, {{0, 0}, {5, 0}, {5, 5}}, 1]],
  TraditionalForm[Row[{"dr = ", drSph[[1, 1]], " dt + ", drSph[[2, 1, 1]], " dW"}]], Top]
```

Neither term contains $\Delta$ or $\Omega$: the Hamiltonian only rotates the Bloch vector. The drift is positive inside $r=\sqrt\eta$ and the noise vanishes at $r=1$, so at $\eta=1$ every run climbs to the surface and stays there: a mixed state purifies and a pure state stays pure. At $\eta<1$ the drift is negative at $r=1$ except at the poles; without a drive the runs still reach $r=1$, which they can do only at a pole, where $\sin\theta=0$, while a drive keeps them off the poles and $r$ keeps wandering below 1.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A monitored qubit stays in the Bloch ball: Itô's lemma and one Reduce

A qubit is Rabi-driven at frequency $\Omega$, relaxes toward $z=1$ at rate $\gamma_1$, dephases at rate $\gamma_\phi$, and is read out along $\sigma_z$ by a detector that dephases it at rate $\Gamma_d$ and whose noise enters $z$ at the informational rate $\Gamma_{CI}$ and rotates the phase at the back-action rate $\Gamma_{BA}$; its efficiency is $\eta=(\Gamma_{CI}+\Gamma_{BA})/2\Gamma_d$. Does every run stay a density matrix, $r=|\vec r|\le1$, for every choice of rates?

With $\hbar=1$ and one Wiener noise $W$, the Bloch vector obeys the Itô equations

$$
\begin{aligned}
dx&=-\left(\tfrac{\gamma_1}{2}+\gamma_\phi+\Gamma_d\right)x\,dt-\left(\sqrt{\Gamma_{CI}}\,xz+\sqrt{\Gamma_{BA}}\,y\right)dW,\\
dy&=\left[-\left(\tfrac{\gamma_1}{2}+\gamma_\phi+\Gamma_d\right)y-\Omega z\right]dt-\left(\sqrt{\Gamma_{CI}}\,yz-\sqrt{\Gamma_{BA}}\,x\right)dW,\\
dz&=\left[\gamma_1(1-z)+\Omega y\right]dt+\sqrt{\Gamma_{CI}}\,(1-z^2)\,dW .
\end{aligned}
$$

WL applies Itô's lemma to $r^2$ on the sphere $r=1$, where $(x,y)=\sqrt{1-z^2}(\cos\varphi,\sin\varphi)$, and `Reduce` decides for which efficiencies the drift points outward somewhere (top). The panels plot that drift against $z$; rates are in units of $\Gamma_d$:

```wl
monitoredQubit = ItoProcess[{
    {-(γ1/2 + γφ + Γd) x[t], -(γ1/2 + γφ + Γd) y[t] - Ω z[t], γ1 (1 - z[t]) + Ω y[t]},
    {{-Sqrt[ΓCI] x[t] z[t] - Sqrt[ΓBA] y[t]}, {-Sqrt[ΓCI] y[t] z[t] + Sqrt[ΓBA] x[t]}, {Sqrt[ΓCI] (1 - z[t]^2)}}}, {{x, y, z}, {x0, y0, z0}}, t];
{r2Drift, r2Noise} = Factor[Simplify[Flatten[monitoredQubit["Derivative", x[t]^2 + y[t]^2 + z[t]^2]] /.
     {x[t] -> Sqrt[1 - z^2] Cos[φ], y[t] -> Sqrt[1 - z^2] Sin[φ], z[t] -> z, ΓBA -> 2 η Γd - ΓCI}]];
ballOutwardCondition = Reduce[Exists[z, -1 < z < 1, r2Drift > 0] && Γd > 0 && γφ >= 0 && γ1 >= 0, η, Reals];
panelBall[g1_, gφ_] := Plot[Evaluate[r2Drift /. {Γd -> 1, γ1 -> g1, γφ -> gφ, η -> #} & /@ {0.5, 1, 1.25, 2}], {z, -1, 1},
    PlotRange -> All, PlotLegends -> (Row[{"η = ", #}] & /@ {0.5, 1, 1.25, 2}), Frame -> True, FrameLabel -> {"z on the sphere", "drift of r² / Γd"}, Axes -> False,
    PlotLabel -> Row[{"γ1/Γd = ", g1, ",   γφ/Γd = ", gφ}], ImageSize -> 300];
panelsBall = Outer[panelBall, {0, 1}, {0, 0.5, 1.5}];
Labeled[Legended[Grid[Map[First, panelsBall, {2}]], Last[panelsBall[[1, 1]]]],
  TraditionalForm[Row[{"noise of r² on the sphere: ", r2Noise, ";    drift outward somewhere  \[DoubleLeftRightArrow]  ", ballOutwardCondition}]], Top]
```

`Reduce` finds an outward drift only for $\eta>1+\gamma_\phi/\Gamma_d$, whatever the drive and the relaxation, and the noise of $r^2$ vanishes on the sphere: a kick only slides the state along it. So for every real detector, $\eta\le1$, the drift on the sphere points inward or along it and no run leaves the ball; relaxation bends the curves but does not move the threshold.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A Hamiltonian with white noise: Stratonovich to Itô gives the dephasing

A qubit with Hamiltonian $\frac{\omega}{2}\sigma_x$ also feels a white-noise field along $z$ with power $\gamma$: $H(t)=\frac{\omega}{2}\sigma_x+\sqrt\gamma\,\xi(t)\,\sigma_z$, with $\xi$ white noise. Every run is a unitary rotation; what does the average over runs do?

With $\hbar=1$, white noise as the limit of smooth noise makes the equation of motion a Stratonovich equation,

$$
d\rho=-i\left[\tfrac{\omega}{2}\sigma_x,\rho\right]dt-i\sqrt\gamma\,[\sigma_z,\rho]\circ dW ,
$$

and only its Itô form, whose drift gains $\frac12(b\cdot\nabla)b$ from the noise $b=2\sqrt\gamma\,(-y,x,0)$ of the Bloch vector, gives the average. WL converts the Bloch form of the equation to Itô (top), samples 200 runs from $z=1$ and plots three of them with their average (thick); QF solves the Lindblad equation $\dot\rho=-i[\frac{\omega}{2}\sigma_x,\rho]+\gamma\,\mathcal D[\sigma_z]\rho$ (dashed). Time is in units of $1/\omega$:

```wl
ρNoisy = (IdentityMatrix[2] + {x[t], y[t], z[t]} . (PauliMatrix /@ {1, 2, 3}))/2;
blochOf[m_] := Simplify[Tr[m . PauliMatrix[#]] & /@ {1, 2, 3}];
noisyH[ω_, γ_] := StratonovichProcess[{blochOf[-I ω/2 (PauliMatrix[1] . ρNoisy - ρNoisy . PauliMatrix[1])],
    List /@ blochOf[-I Sqrt[γ] (PauliMatrix[3] . ρNoisy - ρNoisy . PauliMatrix[3])]}, {{x, y, z}, {0, 0, 1}}, t];
lindbladNoisyH[γ_] := Re @ Tr[Normal[QuantumEvolve[QuantumOperator["X"]/2, {QuantumOperator["Z"]} -> {γ}, QuantumState["0"], {t, 0, 20}]["DensityMatrix"]] . PauliMatrix[3]];
panelNoisyH[γ_] := With[{paths = RandomFunction[noisyH[1, γ], {0., 20., 0.01}, 200, Method -> "StochasticRungeKuttaScalarNoise"]["PathComponent", 3]["Paths"],
     avg = lindbladNoisyH[γ]},
    ListLinePlot[Join[Take[paths, 3], {Transpose[{paths[[1, All, 1]], Mean[paths[[All, All, 2]]]}], Table[{s, avg /. t -> s}, {s, 0., 20., 0.1}]}],
     PlotStyle -> {Automatic, Automatic, Automatic, Thick, Dashed}, PlotRange -> {-1, 1}, Frame -> True, FrameLabel -> {"ω t", "z"},
     PlotLabel -> Row[{"γ/ω = ", γ}], ImageSize -> 300]];
SeedRandom[3]; Labeled[Row[panelNoisyH /@ {0.05, 0.25, 1}],
  TraditionalForm[Row[{"Itô drift: ", ItoProcess[noisyH[ω, γ]]["Drift"] /. {x[t] -> x, y[t] -> y, z[t] -> z}}]], Top]
```

The conversion adds $-2\gamma x$ and $-2\gamma y$ to the drift: dephasing at rate $2\gamma$ about the noise axis. The runs rotate by different random angles about $z$ and scatter, and their average follows the Lindblad curve to $z=0$, faster as $\gamma/\omega$ grows toward 1.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A laser with a diffusing phase: the exact mean field and a non-stationary correlation

A laser field of frequency $\omega$ has a phase that diffuses at rate $\kappa$, starting from zero. What are the mean field and the two-time correlation $\langle E(t)E(t+\tau)\rangle$?

$$
E(t)=\cos\left(\omega t+\phi_t\right),\qquad d\phi=\sqrt{2\kappa}\,dW,\qquad \phi_0=0 .
$$

WL holds $E$ as a `TransformedProcess` of a Wiener process and computes both moments exactly, with $\omega$, $\kappa$, $t$ and $\tau$ as symbols (top). The upper panels show one sampled trace and the mean field (dashed); the lower panels show the correlation at three times $t$. Time is in units of $1/\omega$:

```wl
laserField = TransformedProcess[Cos[ω t + ϕ[t]], ϕ \[Distributed] WienerProcess[0, Sqrt[2 κ]], t];
laserMeanField = Mean[laserField[t]];
laserCorrelation = FullSimplify[ComplexExpand[CovarianceFunction[laserField, t, t + τ] + laserMeanField (laserMeanField /. t -> t + τ)],
    t > 0 && τ > 0 && ω ∈ Reals];
laserTracesPanel[k_] := ListLinePlot[Append[RandomFunction[laserField /. {ω -> 1, κ -> k}, {0., 30., 0.02}]["Paths"],
     Table[{s, laserMeanField /. {ω -> 1, κ -> k, t -> s}}, {s, 0., 30., 0.05}]], PlotStyle -> {Automatic, Dashed},
    PlotRange -> {-1, 1}, Frame -> True, FrameLabel -> {"ω t", "E(t)"}, PlotLabel -> Row[{"κ/ω = ", k}], ImageSize -> 300];
laserCorrelationPanel[k_] := Plot[Evaluate[laserCorrelation /. {ω -> 1, κ -> k, t -> #/k} & /@ {0, 0.25, 2}], {τ, 0, 20},
    PlotRange -> {-1, 1}, PlotLegends -> (Row[{"κ t = ", #}] & /@ {0, 0.25, 2}), Frame -> True,
    FrameLabel -> {"ω τ", "⟨E(t) E(t+τ)⟩"}, PlotLabel -> Row[{"κ/ω = ", k}], ImageSize -> 300];
SeedRandom[4]; With[{lower = laserCorrelationPanel /@ {0.02, 0.1, 0.5}},
  Labeled[Column[{Row[laserTracesPanel /@ {0.02, 0.1, 0.5}], Legended[Row[First /@ lower], Last[First[lower]]]}],
   Column[TraditionalForm /@ ({Row[{"⟨E(t)⟩ = ", laserMeanField}], Row[{"⟨E(t) E(t+τ)⟩ = ", laserCorrelation}]} /.
       Exp[Times[-1, e__]] :> Exp[-HoldForm[Times[e]]]), Center], Top]]
```

The mean field oscillates at $\omega$ and decays as $e^{-\kappa t}$, although no trace loses amplitude: the decay is the spread of phases among runs. The correlation has a part that depends only on the lag $\tau$ and a part that depends on $t$ itself, because every run starts from the same phase. That second part dies on the time scale $1/4\kappa$, and by $\kappa t=2$ only the stationary $\frac12e^{-\kappa\tau}\cos\omega\tau$ is left.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### NDSolve on a smooth noise path: the population read off the record

A qubit is Rabi-driven at frequency $\Omega$, relaxes toward $z=1$ at rate $\gamma_1$, and is read out along $\sigma_z$ by a detector that dephases it at rate $\Gamma_d$ and whose noise enters $z$ at the informational rate $\Gamma_{CI}$ and rotates the phase at the back-action rate $\Gamma_{BA}$. The detector's integrated readout is $Q$. Is $z(t)=\tanh\left(\operatorname{arctanh}z_0+\sqrt{\Gamma_{CI}}\,Q(t)\right)$, with $z_0=z(0)$, on every run?

With $\hbar=1$, the population and the readout obey the Itô equations, driven by one Wiener noise $W$,

$$
\begin{aligned}
dz&=\left[\gamma_1(1-z)+\Omega y\right]dt+\sqrt{\Gamma_{CI}}\,(1-z^2)\,dW,\\
dQ&=\sqrt{\Gamma_{CI}}\,z\,dt+dW ,
\end{aligned}
$$

and $x$, $y$ follow the Bloch equations in the cell. WL converts all four to Stratonovich form, where on a smooth path the chain rule is the ordinary one. `NDSolve` integrates them on one cubic interpolant of a sampled Wiener path, the same path in every panel, from the pure state $(x,y,z)=(0.8,0,-0.6)$ with $\Gamma_{CI}/\Gamma_d=\Gamma_{BA}/\Gamma_d=0.5$; time is in units of $1/\Gamma_d$:

```wl
readoutStrat = StratonovichProcess[ItoProcess[{
     {-(γ1/2 + Γd) x[t], -(γ1/2 + Γd) y[t] - Ω z[t], γ1 (1 - z[t]) + Ω y[t], Sqrt[ΓCI] z[t]},
     {{-Sqrt[ΓCI] x[t] z[t] - Sqrt[ΓBA] y[t]}, {-Sqrt[ΓCI] y[t] z[t] + Sqrt[ΓBA] x[t]}, {Sqrt[ΓCI] (1 - z[t]^2)}, {1}}},
    {{x, y, z, q}, {0.8, 0, -0.6, 0}}, t]];
smoothWiener = BlockRandom[SeedRandom[5]; Interpolation[Transpose[{Subdivide[0., 4., 400],
      Prepend[Accumulate[RandomVariate[NormalDistribution[0, Sqrt[0.01]], 400]], 0.]}], InterpolationOrder -> 3]];
panelTanh[Ωv_, g1_] := With[{rates = {Γd -> 1, ΓCI -> 0.5, ΓBA -> 0.5, Ω -> Ωv, γ1 -> g1}},
   With[{sol = NDSolveValue[{Thread[D[{x[t], y[t], z[t], q[t]}, t] == (readoutStrat["Drift"] + Flatten[readoutStrat["Diffusion"]] smoothWiener'[t] /. rates)],
        Thread[{x[0], y[0], z[0], q[0]} == {0.8, 0, -0.6, 0}]}, {z, q}, {t, 0, 4}, MaxStepSize -> 0.0025]},
    Plot[{sol[[1]][s], Tanh[ArcTanh[-0.6] + Sqrt[0.5] sol[[2]][s]]}, {s, 0, 4}, PlotStyle -> {Automatic, Dashed}, PlotRange -> {-1, 1},
     PlotLegends -> {"z(t)", "tanh(arctanh z0 + √ΓCI Q(t))"}, Frame -> True, FrameLabel -> {"Γd t", "z"},
     PlotLabel -> Row[{"Ω/Γd = ", Ωv, ",   γ1/Γd = ", g1}], ImageSize -> 300]]];
panelsTanh = Outer[panelTanh, {0, 1}, {0, 0.5}];
Legended[Grid[Map[First, panelsTanh, {2}]], Last[panelsTanh[[1, 1]]]]
```

Without drive or relaxation the population is the readout passed through $\tanh$ on this path, back-action included: in Stratonovich form the ordinary chain rule gives $d\operatorname{arctanh}z=\sqrt{\Gamma_{CI}}\,dQ$, so the record alone fixes the state. A drive or relaxation moves $z$ without writing it into the record, and the two curves part.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Quantum Bayesian update: each reading divides the mixedness by $w^2$

A detector reads a Pauli observable $A$ of a qubit over a short window and returns a noisy reading $I$. Two detectors take turns, one reading $\sigma_z$ and one reading $\sigma_\vartheta=\cos\vartheta\,\sigma_z+\sin\vartheta\,\sigma_x$. How does the mixedness $1-r^2$ of the state change with each reading?

Each reading updates the state by Bayes' rule, written as a Kraus operator (Korotkov):

$$
\rho\;\longrightarrow\;\frac{M\rho M}{w},\qquad M(I)=e^{IA/2D},\qquad w=\mathrm{Tr}\,M\rho M ,
$$

and $I$ is drawn from a Gaussian of variance $D$ centered on $+1$ or $-1$ with the Born probabilities of $A$. WL computes $\det M$ for every reading and every $\vartheta$ (top), then follows three runs of 20 readings from $(x,y,z)=(0.3,0,0.4)$, plotting $1-r^2$ of the updated state (lines) and $(1-r_0^2)/\prod_k w_k^2$ (markers):

```wl
krausGauss[i_, A_, d_] := MatrixExp[i/(2 d) A];
detKraus = FullSimplify[Det[krausGauss[iRead, Cos[ϑ] PauliMatrix[3] + Sin[ϑ] PauliMatrix[1], dVar]], {iRead, ϑ} ∈ Reals && dVar > 0];
bayesStep[{ρ_, m_}, A_, d_] := With[{p = Re[Tr[ρ . (IdentityMatrix[2] + A)/2]]},
   With[{mk = krausGauss[RandomVariate[NormalDistribution[RandomChoice[{p, 1 - p} -> {1, -1}], Sqrt[d]]], A, d]},
    With[{w = Re[Tr[mk . ρ . mk]]}, {mk . ρ . mk/w, m/w^2}]]];
mixedness[ρ_] := 1 - Total[(Re[Tr[ρ . PauliMatrix[#]]] & /@ {1, 2, 3})^2];
panelBayes[d_, ϑ_] := With[{runs = Table[FoldList[bayesStep[#1, #2, d] &, {(IdentityMatrix[2] + 0.3 PauliMatrix[1] + 0.4 PauliMatrix[3])/2, 0.75},
        Flatten[Table[{PauliMatrix[3], Cos[ϑ] PauliMatrix[3] + Sin[ϑ] PauliMatrix[1]}, 10], 1]], 3]},
    Show[ListLogPlot[Map[mixedness @* First, runs, {2}], Joined -> True, DataRange -> {0, 20}, PlotRange -> {10^-12, 1}], ListLogPlot[runs[[All, All, 2]], DataRange -> {0, 20}, PlotMarkers -> Automatic],
     Frame -> True, FrameLabel -> {"reading number", "1 - r²"}, PlotLabel -> Row[{"D = ", d, ",   ϑ = ", ϑ}], ImageSize -> 300]];
SeedRandom[6]; Labeled[Grid[Outer[panelBayes, {2, 8}, {0, π/4, π/2}]], TraditionalForm[Row[{"det M = ", detKraus}]], Top]
```

Every Kraus operator has unit determinant, for every reading and every $\vartheta$, and $\det\rho=(1-r^2)/4$, so each update divides the mixedness by exactly $w^2$: the markers sit on the lines. The mixedness falls whether the two detectors commute ($\vartheta=0$) or not, steeply at $D=2$ and slowly at $D=8$, where one reading carries little information.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Born's rule as a martingale: runs collapse to the poles with probability $(1\pm z_0)/2$

A qubit with no Hamiltonian starts at $z=z_0$ and a detector of efficiency $\eta$ measures $\sigma_z$ continuously at strength $k$. Each run ends at a pole; with what probability does it end at $z=+1$?

Without a Hamiltonian, the $z$ component of the conditioned state obeys the Itô equation $dz=2\sqrt{\eta k}\,(1-z^2)\,dW$.

WL samples 400 runs; each panel shows four of them, the mean over all 400 (thick), and the fraction of runs that end above $z=0$ next to $(1+z_0)/2$. Time is in units of $1/k$:

```wl
collapseZ[η_, z0_] := ItoProcess[{0, 2 Sqrt[η] (1 - z[t]^2)}, {z, z0}, t];
panelBorn[η_, z0_] := With[{paths = RandomFunction[collapseZ[η, z0], {0., 4., 0.005}, 400, Method -> "StochasticRungeKuttaScalarNoise"]["Paths"]},
   ListLinePlot[Append[Take[paths, 4], Transpose[{paths[[1, All, 1]], Mean[paths[[All, All, 2]]]}]],
    PlotStyle -> {Automatic, Automatic, Automatic, Automatic, Thick}, PlotRange -> {-1.05, 1.05}, Frame -> True, FrameLabel -> {"k t", "z"},
    PlotLabel -> Column[{Row[{"η = ", η, ",   z0 = ", z0}],
       Row[{"ending above 0: ", Round[N[Count[paths[[All, -1, 2]], _?Positive]/400], 0.01], ",   (1 + z0)/2 = ", (1 + z0)/2}]}, Center], ImageSize -> 300]];
SeedRandom[7]; Grid[Outer[panelBorn, {0.3, 1}, {-0.5, 0.2, 0.8}]]
```

The drift of $z$ is zero, so the mean over runs stays at $z_0$ while every run heads for a pole. A bounded process with constant mean that ends at $\pm1$ must end at $+1$ with probability $(1+z_0)/2$, and the fraction of runs ending above zero matches it within sampling error at both efficiencies; a weaker detector only slows the collapse.
