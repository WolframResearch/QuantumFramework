```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A driven Kerr cavity: sub- or super-Poissonian steady light

A cavity mode of frequency $\omega_c$ with a Kerr nonlinearity $K$ is driven at amplitude $F$ and frequency $\omega_d$, and leaks photons at rate $\kappa$. Drive and loss balance in a steady state. How do the light the cavity settles into, and the time it takes to get there, depend on detuning, Kerr strength and drive?

With $\hbar=1$, the Hamiltonian, and its form in the frame rotating at $\omega_d$, are

$$
H=\omega_c\,a^\dagger a+\frac K2\,a^{\dagger2}a^2+F\left(a\,e^{i\omega_d t}+a^\dagger e^{-i\omega_d t}\right)\;\longrightarrow\;-\Delta\,a^\dagger a+\frac K2\,a^{\dagger2}a^2+F\left(a+a^\dagger\right),\qquad \Delta=\omega_d-\omega_c,
$$

and the loss enters through $\dot\rho=-i[H,\rho]+\kappa\,\mathcal D[a]\rho=\mathcal L\rho$. Light is sub-Poissonian, which no classical field is, when the Mandel parameter $Q=(\langle n^2\rangle-\langle n\rangle^2-\langle n\rangle)/\langle n\rangle$ is negative.

QF builds $\mathcal L$ as a matrix. Its exponential carries the cavity from the vacuum and gives $Q(t)$, its null vector is the steady state (dashed), and its smallest nonzero decay rate, the gap $\lambda$ that sets the slowest transient, is in each legend. Rates are in units of $\kappa$, with 24 Fock levels:

```wl
kerrA = AnnihilationOperator[24];
kerrH[Δ_, K_, F_] := -Δ SuperDagger[kerrA] @ kerrA + K/2 SuperDagger[kerrA] @ SuperDagger[kerrA] @ kerrA @ kerrA + F (kerrA + SuperDagger[kerrA]);
kerrQ[ρ_] := With[{p = Re[Diagonal[ρ]/Tr[ρ]], m = Range[0, 23]}, (m^2 . p - (m . p)^2)/(m . p) - 1];
kerrSolve[Δ_, K_, F_] := With[{L = N @ Normal @ QuantumOperator["Liouvillian"[kerrH[Δ, K, F], {kerrA}, {1}]]["Matrix"]}, With[{step = MatrixExp[0.1 L]},
    {Transpose[{Range[0.1, 20, 0.1], kerrQ[Partition[#, 24]] & /@ Rest @ NestList[step . # &, UnitVector[576, 1], 200]}],
     kerrQ[Partition[First @ NullSpace[L], 24]], Min @ Select[-Re @ Eigenvalues[L], # > 10^-8 &]}]];
kerrPanel[Δ_, K_] := With[{Fs = {0.5, 1.5, 3}}, With[{ss = kerrSolve[Δ, K, #] & /@ Fs},
    Show[
      ListLinePlot[ss[[All, 1]], PlotRange -> All,
        PlotLegends -> MapThread[Row[{"F/\[Kappa] = ", #1, ",   \[Lambda]/\[Kappa] = ", NumberForm[#2, 2]}] &, {Fs, ss[[All, 3]]}]],
      Plot[Evaluate[ss[[All, 2]]], {t, 0, 20}, PlotStyle -> Dashed],
      PlotRange -> {-0.7, 2}, Frame -> True, FrameLabel -> {"\[Kappa] t", "Mandel Q"},
      PlotLabel -> Row[{"\[CapitalDelta]/\[Kappa] = ", Δ, ",   K/\[Kappa] = ", K}], ImageSize -> 330]]];
Grid[Outer[kerrPanel, {0, 2.5}, {0.5, 2.5, 10}]]
```

On resonance every steady state is sub-Poissonian, and the cavity reaches it within a few photon lifetimes, with gaps near the empty-cavity value $\kappa/2$. Detuning makes the sign of $Q$ depend on the drive and the Kerr strength: some steady states become super-Poissonian, and in one of them the gap nearly closes, so $Q$ is still far from its steady value after twenty photon lifetimes.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Photon loss and the Heisenberg advantage of a NOON state

A NOON state $(|n,0\rangle+|0,n\rangle)/\sqrt2$ enters the two arms of an interferometer; arm 2 carries an unknown phase $\theta$, and arm $j$ transmits each photon with probability $\eta_j$. Any estimate of $\theta$ has error at least $1/\sqrt{F}$, where $F$ is the quantum Fisher information, and $n$ independent photons reach $F=n$, the shot-noise limit. Does the NOON advantage keep growing with $n$ once photons are lost?

The phase acts as $e^{-i\theta\hat n_2}$, each arm loses photons through the beam-splitter loss channel of transmission $\eta_j$, and for the lossy state $\rho=\sum_i\lambda_i|i\rangle\langle i|$

$$
F=2\sum_{\lambda_i+\lambda_j>0}\frac{(\lambda_i-\lambda_j)^2}{\lambda_i+\lambda_j}\,\left|\langle i|\hat n_2|j\rangle\right|^2 .
$$

QF applies both losses with `QuantumInterferometer`, and the cell evaluates $F$ from the eigenpairs of the output; points are the gain $F/n$ over shot noise (the gridline at 1), dashed curves are $2n(\eta_1\eta_2)^n/(\eta_1^n+\eta_2^n)$:

```wl
noonState[n_] := (FockState[{n, 0}, n + 1] + FockState[{0, n}, n + 1])/Sqrt[2];
noonQFI[ρ_, G_] := With[{λ = Re[First[#]], v = Last[#]},
    2 Total[Outer[Subtract, λ, λ]^2/(Outer[Plus, λ, λ] + 10.^-12) Abs[Conjugate[v] . G . Transpose[v]]^2, 2]] & @ Eigensystem[ρ];
noonGain[η1_, η2_, n_] := noonQFI[Normal @ QuantumInterferometer[{"Loss"[η1] -> 1, "Loss"[η2] -> 2}, 2][noonState[n]]["DensityMatrix"],
    KroneckerProduct[IdentityMatrix[n + 1], DiagonalMatrix[Range[0, n]]]]/n;
noonPanel[η1_] := With[{η2s = {1, 0.95, 0.9, 0.8}}, Show[
    ListPlot[Table[{n, noonGain[η1, η2, n]}, {η2, η2s}, {n, 12}], PlotLegends -> (Row[{Subscript["\[Eta]", 2], " = ", #}] & /@ η2s)],
    Plot[Evaluate[2 n (η1 #)^n/(η1^n + #^n) & /@ η2s], {n, 1, 12}, PlotStyle -> Dashed],
    PlotRange -> {0, 9}, Frame -> True, FrameLabel -> {"n", "F/n"}, GridLines -> {None, {1}},
    PlotLabel -> Row[{Subscript["\[Eta]", 1], " = ", η1}], ImageSize -> 300]];
noonPanels = noonPanel /@ {0.95, 0.9, 0.8};
Legended[Row[First /@ noonPanels], Last[First[noonPanels]]]
```

The points sit on the dashed curves, which for equal arms reduce to $F=n^2\eta^n$. The phase is carried only by the coherence between $|n,0\rangle$ and $|0,n\rangle$, and loss multiplies it by $(\eta_1\eta_2)^{n/2}$, the amplitude for every photon of both branches to survive. Since $(\eta_1\eta_2)^n$ eventually beats $n^2$, the gain peaks even if only one arm loses photons, and the more loss, the smaller the $n$ at which it peaks; with both arms at $\eta=0.8$ it falls back below shot noise within the plotted range.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Collapse and revival of an atom in a coherent field

An excited two-level atom, frequency $\omega_a$, enters a cavity of frequency $\omega_c$ that holds a coherent field with mean photon number $\bar n$, and couples to it at strength $g$. Does the atom's Rabi oscillation survive the spread of photon numbers in the field?

With $\hbar=1$ and $\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, the Jaynes-Cummings Hamiltonian, and its form in the frame rotating at $\omega_c$, are

$$
H=\omega_c\,a^\dagger a+\frac{\omega_a}{2}\,\sigma_z+g\left(\sigma_-a^\dagger+\sigma_+a\right)\;\longrightarrow\;\frac{\Delta}{2}\,\sigma_z+g\left(\sigma_-a^\dagger+\sigma_+a\right),\qquad \Delta=\omega_a-\omega_c .
$$

$H$ couples $|e,n\rangle$ only to $|g,n+1\rangle$, so the propagator is a $2\times2$ block for each $n$, with off-diagonal element $g\sqrt{n+1}$. WL exponentiates the block with $n$, $g$, $\Delta$ and $t$ symbolic, takes the inversion $w_n(t)=P_e-P_g$ from $|e,n\rangle$ (shown above the figure), and weights it with the Poisson distribution of the coherent field; time is in units of $1/g$:

```wl
jcBlock = {{Δ/2, g Sqrt[n + 1]}, {g Sqrt[n + 1], -Δ/2}};
jcU = FullSimplify[MatrixExp[-I t jcBlock], g > 0 && n >= 0 && Δ ∈ Reals];
jcInversion = FullSimplify[ComplexExpand[Abs[jcU[[1, 1]]]^2 - Abs[jcU[[2, 1]]]^2], g > 0 && n >= 0 && t > 0 && Δ ∈ Reals];
jcRecord[nbar_, Δv_] := Total @ Table[N[PDF[PoissonDistribution[nbar], m]] (jcInversion /. {n -> m, g -> 1, Δ -> Δv}), {m, 0, Ceiling[nbar + 8 Sqrt[nbar] + 10]}];
jcPanel[nbar_, Δv_] := Plot[Evaluate[jcRecord[nbar, Δv]], {t, 0, 60}, PlotPoints -> 2500, MaxRecursion -> 0, PlotRange -> {-1, 1},
    Frame -> True, FrameLabel -> {"g t", Row[{Subscript["P", "e"], " - ", Subscript["P", "g"]}]},
    PlotLabel -> Row[{OverBar["n"], " = ", nbar, ",   \[CapitalDelta]/g = ", Δv}], ImageSize -> 330];
Labeled[Grid[Outer[jcPanel, {5, 20}, {0, 5, 20}]], Row[{Subscript["w", "n"], "(t) = ", TraditionalForm[jcInversion]}], Top]
```

The oscillation does not survive at first: the photon numbers rotate the atom at different frequencies $\sqrt{\Delta^2+4g^2(n+1)}$, so the inversion collapses within a few $1/g$, and it revives when their phases realign. The revivals stand apart cleanly only in the larger field. Detuning lifts the level around which the inversion oscillates, because each $w_n$ is pulled toward 1, and delays the revival; at $\Delta/g=20$ the atom hardly emits.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### How much cavity Wigner negativity the atomic record certifies

An excited atom meets a cavity field, and after a time $t$ only the atom's inversion $z=P_e-P_g$ is recorded. The cavity, coupled to the atom at strength $g$, detuned by $\Delta$ and leaking photons at rate $\kappa$, can develop a Wigner function $W$ with negative regions. Fields with the same photon-number distribution $p_n$ give the atom the same record; how much negative volume $\int\max(0,-W)\,dx\,dp$ do all of them share?

With $\hbar=1$ and $\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, the Jaynes-Cummings Hamiltonian, and its form in the frame rotating at $\omega_c$, are

$$
H=\omega_c\,a^\dagger a+\frac{\omega_a}{2}\,\sigma_z+g\left(\sigma_-a^\dagger+\sigma_+a\right)\;\longrightarrow\;\frac{\Delta}{2}\,\sigma_z+g\left(\sigma_-a^\dagger+\sigma_+a\right),\qquad \Delta=\omega_a-\omega_c,
$$

and the loss enters through $\dot\rho=-i[H,\rho]+\kappa\,\mathcal D[a]\rho$. The coupling conserves the excitation number $\sigma_+\sigma_-+a^\dagger a$ and loss only lowers it, so an atom starting in $|e\rangle$ writes into $z(t)$ only the photon-number distribution $p_n$ of the field. Averaging the field over its phase keeps $p_n$ and averages $W$ over the angle, which can cancel negative regions but never create them.

QF evolves the atom in $|e\rangle$ with two fields of mean photon number 4: a coherent state, and the same state with a uniformly random phase, $\sum_np_n|n\rangle\langle n|$. The top row is the record, the bottom row the negative volume of the cavity field; the phase-randomized field has a Wigner function that depends only on the radius, so its negative volume is a radial integral. Time is in units of $1/g$:

```wl
negA = AnnihilationOperator[20, {2}]; negField = CoherentState[20][2.];
negFields = {negField, QuantumState[DiagonalMatrix[negField["ProbabilitiesList"]], 20]};
negH[Δ_] := Δ/2 QuantumOperator["Z"] + QuantumOperator["-"] @ SuperDagger[negA] + SuperDagger[QuantumOperator["-"]] @ negA;
negVolume[ρ_] := Total[Ramp[-Flatten[WignerRepresentation[ρ, {-5., 5.}, {-5., 5.}, "GridSize" -> 81]["ValuesOnGrid"]]]]/64;
negRadii = Range[0.005, 6, 0.01];
negFock = Table[(-1)^m LaguerreL[m, 2 negRadii^2] Exp[-negRadii^2]/π, {m, 0, 19}];
negRadial[ρ_] := 2 π 0.01 negRadii . Ramp[-Re[Diagonal[Normal[ρ["DensityMatrix"]]]] . negFock];
negData[{κ_, Δ_}] := Transpose[MapThread[Function[{f, volume}, With[{evo = QuantumEvolve[negH[Δ], {negA} -> {κ}, QuantumTensorProduct[QuantumState["0"], f], {t, 0, 16}]},
    Table[{{s, Re @ Tr[PauliMatrix[3] . Normal @ QuantumPartialTrace[evo[s], {2}]["DensityMatrix"]]}, {s, volume[QuantumPartialTrace[evo[s], {1}]]}}, {s, 0, 16, 0.1}]]],
    {negFields, {negVolume, negRadial}}], {2, 3, 1}];
negPanel[{κ_, Δ_}, data_, {label_, range_}] := ListLinePlot[data, PlotStyle -> {Automatic, Dashed}, PlotRange -> range,
    PlotLegends -> {"coherent field", "phase-randomized field"}, Frame -> True, FrameLabel -> {"g t", label},
    PlotLabel -> Row[{"\[Kappa]/g = ", κ, ",   \[CapitalDelta]/g = ", Δ}], ImageSize -> 260];
negCases = {{0, 0}, {0.1, 0}, {0, 2}, {0.1, 2}};
negPanels = Transpose @ MapThread[Function[{c, d}, MapThread[negPanel[c, #1, #2] &, {d, {{"z", {-1, 1}}, {"negative volume", {0, 0.4}}}}]], {negCases, negData /@ negCases}];
Legended[Grid[Map[First, negPanels, {2}]], Last[negPanels[[1, 1]]]]
```

The two records coincide in every panel, so the atom cannot tell the coherent field from its phase-randomized version. Any field with the same $p_n$ has at least the dashed negative volume, and the phase-randomized field reaches it, so the dashed curve is what they all share; it wiggles because the atom keeps trading photons with the field and reshapes $p_n$ on the Rabi time scale. The rest, carried by coherences between photon numbers, is most of the negativity between collapse and revival, and with loss the guaranteed part is the first to vanish.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Event-triggered control of a decaying qubit

A qubit decays from $|e\rangle$ at rate $\gamma$ under a drive of Rabi frequency $\Omega$, detuned by $\Delta$. A controller that tracks the excited population $P_e$ from the model, without measuring the qubit, switches the drive off when $P_e$ rises above $p+0.02$ and back on when it falls below $p-0.02$. A constant drive lets $P_e$ settle at or below $1/2$; can the switching hold it at a target $p$ above that?

With $\hbar=1$ and $\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, the driven qubit in the frame rotating at the drive frequency is

$$
H=\frac{\omega_q}{2}\sigma_z+u(t)\cos(\omega_d t)\,\sigma_x\;\longrightarrow\;\frac{\Delta}{2}\sigma_z+\frac{u(t)}{2}\sigma_x,\qquad \Delta=\omega_q-\omega_d,\quad u(t)\in\{0,\Omega\},
$$

and $\dot\rho=-i[H,\rho]+\gamma\,\mathcal D[\sigma_-]\rho$ is, for the Bloch vector $(x,y,z)$ with $P_e=(1+z)/2$,

$$
\dot x=-\Delta y-\tfrac{\gamma}{2}x,\qquad \dot y=\Delta x-\tfrac{\gamma}{2}y-u\,z,\qquad \dot z=u\,y-\gamma(1+z).
$$

`NDSolveValue` integrates these from $|g\rangle$ with $u$ as a discrete variable that two `WhenEvent`s switch, and QF gives the constant drive $u=\Omega$ (dashed); time is in units of $1/\gamma$:

```wl
controlRun[Ω_, Δ_, p_] := NDSolveValue[{controlX'[t] == -Δ controlY[t] - controlX[t]/2,
    controlY'[t] == Δ controlX[t] - controlY[t]/2 - controlU[t] controlZ[t], controlZ'[t] == controlU[t] controlY[t] - 1 - controlZ[t],
    controlX[0] == 0, controlY[0] == 0, controlZ[0] == -1, controlU[0] == Ω,
    WhenEvent[controlZ[t] > 2 p - 0.96, controlU[t] -> 0], WhenEvent[controlZ[t] < 2 p - 1.04, controlU[t] -> Ω]},
    controlZ, {t, 0, 3}, DiscreteVariables -> {controlU}];
controlConstant[Ω_, Δ_] := Re @ Normal[QuantumEvolve[Δ/2 QuantumOperator["Z"] + Ω/2 QuantumOperator["X"], {QuantumOperator["-"]} -> {1},
    QuantumState["1"], {t, 0, 3}]["DensityMatrix"]][[1, 1]];
controlPanel[Ω_, Δr_] := With[{ps = {0.6, 0.75, 0.9}},
    Plot[Evaluate[Append[(1 + controlRun[Ω, Δr Ω, #][t])/2 & /@ ps, controlConstant[Ω, Δr Ω]]], {t, 0, 3},
      PlotStyle -> {Automatic, Automatic, Automatic, Directive[Thick, Dashed]}, PlotPoints -> 400, PlotRange -> {0, 1},
      PlotLegends -> Append[Row[{"p = ", #}] & /@ ps, "constant drive"], Frame -> True, FrameLabel -> {"\[Gamma] t", Subscript["P", "e"]},
      PlotLabel -> Row[{"\[CapitalOmega]/\[Gamma] = ", Ω, ",   \[CapitalDelta]/\[CapitalOmega] = ", Δr}], ImageSize -> 300]];
controlPanels = Outer[controlPanel, {5, 15, 40}, {0, 0.5}];
Legended[Grid[Map[First, controlPanels, {2}]], Last[controlPanels[[1, 1]]]]
```

Without detuning, whenever the drive can reach the band the switching holds $P_e$ inside it, above $1/2$, for a time that shortens as the target rises; then the drive stays on and $P_e$ oscillates about $1/2$ as under constant drive. A curve that never reaches its band coincides with the constant-drive curve. With $\Delta/\Omega=0.5$ the tilted axis keeps $P_e$ below $\Omega^2/(\Omega^2+\Delta^2)$; the switching clips only the first peaks, then the oscillation decays below every band.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A cavity mode thermalizes

A cavity mode of frequency $\omega$ exchanges photons with a bath at temperature $T$: it loses photons at rate $\kappa(n_{\rm th}+1)$ and gains them at rate $\kappa\,n_{\rm th}$, with $n_{\rm th}=1/(e^{\hbar\omega/k_BT}-1)$. Does the bath fix only the mean photon number, or the whole photon-number distribution?

With $\hbar=1$ the Hamiltonian is $H=\omega\,a^\dagger a$, and the Lindblad equation is

$$
\dot\rho=-i[\omega\,a^\dagger a,\rho]+\kappa\,(n_{\rm th}+1)\,\mathcal D[a]\rho+\kappa\,n_{\rm th}\,\mathcal D[a^\dagger]\rho .
$$

In the frame rotating at $\omega$ the Hamiltonian drops out and the two dissipators stay as they are. A steady state must balance the flow from $n$ down to $n-1$ against the flow back up, $\kappa(n_{\rm th}+1)\,n\,p_n=\kappa\,n_{\rm th}\,n\,p_{n-1}$, which fixes $p_n=n_{\rm th}^n/(n_{\rm th}+1)^{n+1}$.

QF evolves three initial states with 32 Fock levels and plots the photon-number distribution $p_n$ at three times (points); the dashed line is the thermal distribution:

```wl
thermA = AnnihilationOperator[32];
thermP[nth_, ψ0_] := With[{evo = QuantumEvolve[0 SuperDagger[thermA] @ thermA, {thermA, SuperDagger[thermA]} -> {nth + 1, nth}, ψ0, {t, 0, 8}]},
    Table[Re @ Diagonal @ Normal @ evo[s]["DensityMatrix"], {s, {0.25, 1, 8}}]];
thermPanel[nth_, {label_, ψ0_}] := Show[
    ListLogPlot[thermP[nth, ψ0][[All, ;; 21]], DataRange -> {0, 20}, PlotMarkers -> Automatic,
      PlotLegends -> (Row[{"\[Kappa] t = ", #}] & /@ {0.25, 1, 8})],
    ListLogPlot[Table[nth^m/(nth + 1)^(m + 1), {m, 0, 20}], DataRange -> {0, 20}, Joined -> True, PlotStyle -> Dashed],
    PlotRange -> {Log[10^-6], 0}, Axes -> False, Frame -> True, FrameLabel -> {"n", Subscript["p", "n"]},
    PlotLabel -> Row[{Subscript["n", "th"], " = ", nth, ",   start ", label}], ImageSize -> 300];
thermPanels = Outer[thermPanel, {0.5, 2}, {{"|0\[RightAngleBracket]", FockState[0, 32]}, {"|4\[RightAngleBracket]", FockState[4, 32]}, {"|\[Alpha] = 2\[RightAngleBracket]", CoherentState[32][2.]}}, 1];
Legended[Grid[Map[First, thermPanels, {2}]], Last[thermPanels[[1, 1]]]]
```

Whatever the start, by $\kappa t=8$ every distribution lies on the dashed thermal line, straight on this logarithmic scale with a slope set by $n_{\rm th}$ alone. The Fock state $|4\rangle$ and the coherent state $|\alpha=2\rangle$ start with the same mean photon number and different distributions, and the bath erases the difference: it fixes the whole distribution, not only its mean. From the vacuum the distribution fills from low $n$ upward.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A cat in a warm bath: Caldeira-Leggett decoherence

A particle of mass $M$ in a trap of frequency $\Omega$ couples to a bath at temperature $T$ that damps its momentum at rate $2\gamma$. It starts in a cat state, the superposition $|\alpha\rangle+|-\alpha\rangle$ of two coherent states on either side of the trap. How fast does the bath turn it into a mixture of the two lobes?

With $\hbar=1$, the high-temperature Caldeira-Leggett master equation is

$$
\dot\rho=-i\Big[\tfrac{p^2}{2M}+\tfrac12M\Omega^2x^2,\rho\Big]-i\gamma\big[x,\{p,\rho\}\big]-2\gamma Mk_BT\,\big[x,[x,\rho]\big]-\frac{\gamma}{8Mk_BT}\big[p,[p,\rho]\big],
$$

where the last term, smaller than the others by a factor of order $(\Omega/k_BT)^2$, makes the generator a Lindblad generator. With jump operators $x$ and $p$ it takes the Hamiltonian $H+\tfrac\gamma2\{x,p\}$ and the rate matrix $\begin{pmatrix}4\gamma Mk_BT & -i\gamma\\ i\gamma & \gamma/4Mk_BT\end{pmatrix}$, which QF accepts directly.

QF evolves the cat $N(|\alpha\rangle+|-\alpha\rangle)$, with $N^2=1/(2(1+e^{-2\alpha^2}))$, and the mixture of its two lobes through the same bath, with 26 Fock levels of the trap. The coherence is the trace norm of $\rho_{\rm cat}-N^2(|\alpha\rangle\langle\alpha|+|-\alpha\rangle\langle-\alpha|)$, the part made of the cross terms $|\alpha\rangle\langle-\alpha|$ and $|-\alpha\rangle\langle\alpha|$, relative to its initial value; units are $M=\Omega=1$:

```wl
qbmA = AnnihilationOperator[26]; qbmX = (qbmA + SuperDagger[qbmA])/Sqrt[2]; qbmP = I (SuperDagger[qbmA] - qbmA)/Sqrt[2];
qbmCoherence[γ_, kT_, α_] := With[{evs = QuantumEvolve[(qbmP @ qbmP + qbmX @ qbmX)/2 + γ/2 (qbmX @ qbmP + qbmP @ qbmX),
      {qbmX, qbmP} -> {{4 γ kT, -I γ}, {I γ, γ/(4 kT)}}, #, {t, 0, 2 π}] & /@
      {CatState[26][α, 0], QuantumState[(CoherentState[26][α]["DensityMatrix"] + CoherentState[26][-α]["DensityMatrix"])/2, 26]}},
    With[{c = Table[Total @ SingularValueList[Normal[evs[[1]][s]["DensityMatrix"] - evs[[2]][s]["DensityMatrix"]/(1 + Exp[-2 α^2])]], {s, 0, 2 π, π/32}]},
      Transpose[{Range[0, 2 π, π/32], c/First[c]}]]];
qbmPanel[α_, kT_] := ListLogPlot[qbmCoherence[#, kT, α] & /@ {0.001, 0.003, 0.01}, Joined -> True, PlotRange -> {10^-4, 1.1},
    PlotLegends -> (Row[{"\[Gamma]/\[CapitalOmega] = ", #}] & /@ {0.001, 0.003, 0.01}), Frame -> True, FrameLabel -> {"\[CapitalOmega] t", "coherence"},
    FrameTicks -> {{Automatic, None}, {Table[{k π/2, k π/2}, {k, 0, 4}], None}}, GridLines -> {Table[k π/2, {k, 4}], None},
    PlotLabel -> Row[{"\[Alpha] = ", α, ",   ", Subscript["k", "B"], "T/\[CapitalOmega] = ", kT}], ImageSize -> 300];
qbmPanels = Outer[qbmPanel, {1, 2}, {5, 20}];
Legended[Grid[Map[First, qbmPanels, {2}]], Last[qbmPanels[[1, 1]]]]
```

The coherence falls in steps. It drops while the two lobes sit apart in position, near $\Omega t=0,\pi,2\pi$, and stalls near $\pi/2$ and $3\pi/2$, where the trap has turned their separation into momentum, which the $[x,[x,\rho]]$ term does not see. A larger cat, a hotter bath or stronger friction makes each step deeper, and the coherence is lost over the plotted time, during which friction at rate $2\gamma\ll\Omega$ removes almost none of the motion.
