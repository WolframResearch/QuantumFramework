```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A driven, leaky cavity settles into a coherent state

A cavity mode of frequency $\omega_c$ is driven by a classical field of frequency $\omega_d$ and strength $G$, and its mirror leaks photons at rate $\kappa$. Does the field end in the same pure coherent state whatever its start?

With $\hbar=1$, the Hamiltonian, and its form in the frame rotating at $\omega_d$, are

$$
H=\omega_c\,a^\dagger a+G\left(a\,e^{i\omega_dt}+a^\dagger e^{-i\omega_dt}\right)\;\longrightarrow\;\Delta\,a^\dagger a+G\left(a+a^\dagger\right),\qquad \Delta=\omega_c-\omega_d,
$$

the loss enters as $\dot\rho=-i[H,\rho]+\kappa\,\mathcal D[a]\rho$ with $\mathcal D[L]\rho=L\rho L^\dagger-\tfrac12\{L^\dagger L,\rho\}$, and the mean field obeys $\langle\dot a\rangle=-(i\Delta+\kappa/2)\langle a\rangle-iG$, whose fixed point is $\beta=-iG/(\kappa/2+i\Delta)$.

QF evolves the vacuum, the number state $|3\rangle$ and a thermal state with $\bar n=2$ on 30 Fock levels and plots their overlap $\langle\beta|\rho(t)|\beta\rangle$ with the coherent state $|\beta\rangle$; rates are in units of $\kappa$:

```wl
balA = AnnihilationOperator[30];
balStarts = {FockState[0, 30], FockState[3, 30], ThermalState[2., 30]};
balOverlap[G_, Δ_] := With[{v = Normal[CoherentState[30][N[-I G/(1/2 + I Δ)]]["StateVector"]]},
   Re[Conjugate[v] . Normal[QuantumEvolve[Δ SuperDagger[balA] @ balA + G (balA + SuperDagger[balA]), {balA} -> {1}, #, {t, 0, 8}]["DensityMatrix"]] . v] & /@ balStarts];
balPanel[G_, Δ_] := Plot[Evaluate[balOverlap[G, Δ]], {t, 0, 8}, PlotRange -> {0, 1.02}, Frame -> True,
    PlotLegends -> {"|0\[RightAngleBracket]", "|3\[RightAngleBracket]", Row[{"thermal, ", OverBar["n"], " = 2"}]},
    FrameLabel -> {"\[Kappa] t", "\[LeftAngleBracket]\[Beta]|\[Rho]|\[Beta]\[RightAngleBracket]"},
    PlotLabel -> Row[{"G/\[Kappa] = ", G, ",   \[CapitalDelta]/\[Kappa] = ", Δ}], ImageSize -> 300];
balPanels = Outer[balPanel, {0.5, 1.5}, {0, 1, 3}];
Legended[Grid[Map[First, balPanels, {2}]], Last[balPanels[[1, 1]]]]
```

From every start, the mixed thermal state included, the overlap reaches 1 within a few $1/\kappa$: shifted by $\beta$, the master equation describes a detuned, lossy mode with no drive, which the loss empties into its vacuum, and that vacuum is $|\beta\rangle$. The drive and the detuning only set where $\beta$ lies, so the vacuum and thermal starts begin lowest where $|\beta|=G/\sqrt{\kappa^2/4+\Delta^2}$ is largest.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A parametric oscillator as a Langevin process in phase space

A cavity mode of frequency $\omega_c$ is pumped at $2\omega_p$ by a classical field that creates photons in pairs with strength $\varepsilon$, and leaks photons at rate $\kappa$. Can classical random trajectories in phase space reproduce the noise of its quantum steady state?

With $\hbar=1$, the Hamiltonian, and its form in the frame rotating at $\omega_p$, are

$$
H=\omega_c\,a^\dagger a+\frac{i\varepsilon}{2}\left(a^{\dagger2}e^{-2i\omega_pt}-a^2e^{2i\omega_pt}\right)\;\longrightarrow\;\Delta\,a^\dagger a+\frac{i\varepsilon}{2}\left(a^{\dagger2}-a^2\right),\qquad \Delta=\omega_c-\omega_p,
$$

with the loss $\kappa\,\mathcal D[a]\rho$ in the Lindblad equation. In the quadratures $x=(a+a^\dagger)/\sqrt2$, $p=(a-a^\dagger)/i\sqrt2$, a quadratic $H$ and a linear loss make the Wigner function obey an exact Fokker-Planck equation, which is the law of an Itô process $\mathbf r=(x,p)$ driven by two independent Wiener processes:

$$
\partial_tW=-\nabla\cdot\left(J\mathbf r\,W\right)+\frac\kappa4\nabla^2W
\quad\Longleftrightarrow\quad
d\mathbf r=J\mathbf r\,dt+\sqrt{\kappa/2}\;d\mathbf w,\qquad
J=\begin{pmatrix}\varepsilon-\kappa/2&\Delta\\-\Delta&-\varepsilon-\kappa/2\end{pmatrix}.
$$

WL samples 400 trajectories of the Itô process from the origin and QF evolves the master equation from the vacuum on 30 Fock levels, both to $\kappa t=60$; each panel shows the trajectory end points and the ellipses $\mathbf r^{T}\sigma^{-1}\mathbf r=4$ of their sample covariance and of the master equation's covariance matrix $\sigma$, with rates in units of $\kappa$:

```wl
opoA = AnnihilationOperator[30];
opoCloud[ε_, Δ_] := BlockRandom[SeedRandom[7]; RandomFunction[ItoProcess[{{{ε - 1/2, Δ}, {-Δ, -ε - 1/2}} . {x, p}, Sqrt[1/2] IdentityMatrix[2]}, {x, p}, t],
     {0, 60, 0.02}, 400]["SliceData", 60]];
opoSigma[ε_, Δ_] := CovarianceMatrix[QuantumEvolve[Δ SuperDagger[opoA] @ opoA + I ε/2 (SuperDagger[opoA] @ SuperDagger[opoA] - opoA @ opoA), {opoA} -> {1},
     FockState[0, 30], {t, 0, 60}][60]];
opoEllipse[σ_] := Table[2 MatrixPower[σ, 1/2] . {Cos[u], Sin[u]}, {u, 0, 2 Pi, Pi/60}];
opoPanel[ε_, Δ_] := With[{cloud = opoCloud[ε, Δ]},
   ListPlot[{cloud, opoEllipse[opoSigma[ε, Δ]], opoEllipse[Covariance[cloud]], opoEllipse[IdentityMatrix[2]/2]}, Joined -> {False, True, True, True},
    PlotStyle -> {PointSize[0.006], Thick, Directive[Thick, Dashed], Directive[Thick, Dotted]},
    PlotLegends -> {"Langevin end points", "master equation", "Langevin covariance", "vacuum"},
    AspectRatio -> 1, Axes -> False, PlotRange -> {{-4, 4}, {-4, 4}}, Frame -> True, FrameLabel -> {"x", "p"},
    PlotLabel -> Row[{"\[CurlyEpsilon]/\[Kappa] = ", ε, ",   \[CapitalDelta]/\[Kappa] = ", Δ}], ImageSize -> 300]];
opoPanels = Outer[opoPanel, {0.2, 0.4}, {0, 0.5}];
Legended[Grid[Transpose[Map[First, opoPanels, {2}]]], Last[opoPanels[[1, 1]]]]
```

In every panel the ellipse of the Langevin end points lies on the master equation's: the Fokker-Planck equation is exact here, so a cloud of classical random walkers carries the steady-state Wigner function, vacuum noise included. The pump stretches the cloud along one quadrature and pulls it inside the vacuum circle along the other, more strongly at the larger $\varepsilon$; detuning tilts the squeezed axis and weakens the squeezing.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The phase of a light field

A field mode holds $\bar n$ photons on average, either in the number state $|\bar n\rangle$, in the coherent state $|\alpha\rangle$ with $\alpha=\sqrt{\bar n}\,e^{i\pi/3}$, or in a thermal state. How sharply does each of them define a phase?

Following Susskind and Glogower, the phase exponential $E=(aa^\dagger)^{-1/2}a=\sum_n|n\rangle\langle n+1|$ has the unnormalizable eigenstates $|\varphi\rangle$, which are not orthogonal but resolve the identity, so a state $\rho$ has the phase distribution $P(\varphi)$:

$$
E|\varphi\rangle=e^{i\varphi}|\varphi\rangle,\qquad |\varphi\rangle=\sum_{n\ge0}e^{in\varphi}|n\rangle,\qquad \frac1{2\pi}\int_{-\pi}^{\pi}|\varphi\rangle\langle\varphi|\,d\varphi=1,\qquad P(\varphi)=\frac{\langle\varphi|\rho|\varphi\rangle}{2\pi}.
$$

QF builds the three states on 60 Fock levels, and the cell plots $2\pi P(\varphi)$, so a uniform distribution sits at 1:

```wl
phaseKet[φ_] := Exp[I φ Range[0, 59]];
phaseStates[nb_] := {FockState[nb, 60], CoherentState[60][N[Sqrt[nb] Exp[I Pi/3]]], ThermalState[N[nb], 60]};
phasePanel[nb_] := With[{ms = Normal[#["DensityMatrix"]] & /@ phaseStates[nb], φs = Subdivide[-N[Pi], N[Pi], 360]},
   ListLinePlot[Table[Re[Conjugate[phaseKet[φ]] . m . phaseKet[φ]], {m, ms}, {φ, φs}], DataRange -> {-Pi, Pi}, PlotRange -> All,
    PlotStyle -> {Thick, Automatic, Dashed}, PlotLegends -> {"number state", "coherent state", "thermal state"},
    Frame -> True, FrameTicks -> {{Automatic, None}, {Range[-Pi, Pi, Pi/3], None}}, FrameLabel -> {"\[CurlyPhi]", "2\[Pi] P(\[CurlyPhi])"},
    PlotLabel -> Row[{OverBar["n"], " = ", nb}], ImageSize -> 300]];
phasePanels = phasePanel /@ {1, 4, 16};
Legended[Row[First /@ phasePanels, Spacer[10]], Last[First[phasePanels]]]
```

The number and thermal states lie flat at 1 for every $\bar n$: their density matrices are diagonal in $n$, so $\langle\varphi|\rho|\varphi\rangle=\sum_n\rho_{nn}=1$ and no phase is preferred. The coherent state, a superposition of number states with locked relative phases, peaks at $\arg\alpha=\pi/3$, and the peak narrows as $\bar n$ grows.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Balanced homodyne detection of squeezed light

A signal mode $a$ in the displaced squeezed state $D(\alpha)S(r)|0\rangle$, with real $\alpha$ and $r$, meets a local oscillator, a coherent state $|\beta e^{i\varphi}\rangle$ in mode $b$, on a balanced beam splitter, and two detectors count the photons $N_c$ and $N_d$ in its outputs. How strong must the local oscillator be for the difference count to show the signal's quadrature noise?

With $S(r)=e^{r(a^2-a^{\dagger2})/2}$ and $D(\alpha)=e^{\alpha(a^\dagger-a)}$, the splitter and the quadrature that the oscillator's phase $\varphi$ selects are

$$
B=e^{\frac\pi4\left(a\,b^\dagger-a^\dagger b\right)},\qquad X(\varphi)=\tfrac12\left(a\,e^{-i\varphi}+a^\dagger e^{i\varphi}\right),
$$

and the noise of the difference count is exactly $\operatorname{Var}(N_c-N_d)=4|\beta|^2\operatorname{Var}X(\varphi)+\langle a^\dagger a\rangle$, the second term coming from the signal's own photons beating with the oscillator's vacuum fluctuations.

QF sends the two modes through $B$ on 40 Fock levels each, computes $\operatorname{Var}(N_c-N_d)/4|\beta|^2$ from the joint photon-count distribution (markers), and $\operatorname{Var}X(\varphi)$ from the signal's covariance matrix (dashed):

```wl
hdLevels = 40; hdK = Range[0, hdLevels - 1];
hdSplitter = BeamSplitterOperator[{N[Pi/4], 0}, hdLevels];
hdNoise[sig_, β_][φ_] := With[{lo = PhaseShiftOperator[φ, hdLevels] @ CoherentState[hdLevels][β], dn = Outer[Subtract, hdK, hdK]},
   With[{pr = ArrayReshape[hdSplitter[QuantumTensorProduct[sig, lo]]["ProbabilitiesList"], {hdLevels, hdLevels}]}, (Total[dn^2 pr, 2] - Total[dn pr, 2]^2)/(4 β^2)]];
hdPanel[α_, r_] := With[{sig = (DisplacementOperator[N[α], hdLevels] @ SqueezeOperator[N[r], hdLevels] @ FockState[0, hdLevels])["Normalize"], φs = Subdivide[0., N[Pi], 24]},
   With[{σ = CovarianceMatrix[sig]}, ListLogPlot[Append[Table[Transpose[{φs, hdNoise[sig, β] /@ φs}], {β, {0.5, 1., 3.}}],
      Table[{φ, {Cos[φ], Sin[φ]} . σ . {Cos[φ], Sin[φ]}/2}, {φ, Subdivide[0., N[Pi], 120]}]],
    Joined -> {False, False, False, True}, PlotStyle -> {Automatic, Automatic, Automatic, Dashed}, PlotMarkers -> {Automatic, Automatic, Automatic, None},
    PlotLegends -> {"|\[Beta]| = 0.5", "|\[Beta]| = 1", "|\[Beta]| = 3", "Var X(\[CurlyPhi])"}, PlotRange -> {0.04, 6},
    Frame -> True, FrameTicks -> {{Automatic, None}, {Range[0, Pi, Pi/4], None}},
    FrameLabel -> {"\[CurlyPhi]", Row[{"Var(", Subscript["N", "c"], " - ", Subscript["N", "d"], ") / 4|\[Beta]|\.b2"}]},
    PlotLabel -> Row[{"\[Alpha] = ", α, ",   r = ", r}], ImageSize -> 300]]];
hdPanels = Outer[hdPanel, {0.5, 1}, {0, 0.4, 0.8}];
Legended[Grid[Map[First, hdPanels, {2}]], Last[hdPanels[[1, 1]]]]
```

The markers sit above the dashed curve, so a weak oscillator buries the quadrature noise under the signal's own photons, the more so for the brighter, more displaced signal. Only the strongest oscillator dips below the coherent state's flat $1/4$ near $\varphi=0$ and $\pi$. Even it stays above $\operatorname{Var}X(\varphi)$ there, because the signal's photons must be small next to $|\beta|^2$ times the squeezed variance, not merely next to $|\beta|^2$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Dressed states and atom-field entanglement

A two-level atom with transition frequency $\omega_a$ enters a cavity of frequency $\omega_c$ that holds a coherent field $|\alpha\rangle$ with real $\alpha=\sqrt{\bar n}$, and couples to it at strength $g$. Which initial atomic states stay nearly unentangled with the field?

With $\hbar=1$, $|e\rangle=|0\rangle$ and $\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, the interaction $g(a\sigma_++a^\dagger\sigma_-)$ couples only $|e,n\rangle$ and $|g,n+1\rangle$. In the frame rotating at $\omega_c$, where $H=\tfrac\Delta2\sigma_z+g(a\sigma_++a^\dagger\sigma_-)$ with $\Delta=\omega_a-\omega_c$, each such pair splits into the dressed states

$$
|n,+\rangle=\cos\tfrac{\Phi_n}2|e,n\rangle+\sin\tfrac{\Phi_n}2|g,n+1\rangle,\qquad |n,-\rangle=-\sin\tfrac{\Phi_n}2|e,n\rangle+\cos\tfrac{\Phi_n}2|g,n+1\rangle,\qquad \cos\Phi_n=\frac{\Delta}{\Omega_n},\quad \sin\Phi_n=\frac{2g\sqrt{n+1}}{\Omega_n},
$$

with $\Omega_n=\sqrt{\Delta^2+4g^2(n+1)}$ and energies $\pm\Omega_n/2$ (plus $\omega_c(n+\tfrac12)$ in the lab frame). Replacing $\sqrt{n+1}$ by $\alpha$ gives the atomic states $|\pm_\alpha\rangle$, the eigenstates of the atom in the classical field $\alpha$.

QF evolves $|\psi_A\rangle|\alpha\rangle$ for $|\psi_A\rangle=|e\rangle,|+_\alpha\rangle,|-_\alpha\rangle$ on 40 Fock levels and plots the atom's von Neumann entropy $S=-\mathrm{Tr}\,\rho_A\log_2\rho_A$ with $\rho_A=\mathrm{Tr}_{\rm field}\,|\psi(t)\rangle\langle\psi(t)|$; time is in units of $1/g$:

```wl
jcdA = AnnihilationOperator[40, {2}]; jcdLower = QuantumOperator["-"];
jcdAtoms[nb_, Δ_] := With[{Φ = ArcTan[Δ, 2 Sqrt[nb]]}, QuantumState /@ N[{{1, 0}, {Cos[Φ/2], Sin[Φ/2]}, {-Sin[Φ/2], Cos[Φ/2]}}]];
jcdEntropy[nb_, Δ_][atom_] := With[{ρA = Normal[QuantumPartialTrace[QuantumEvolve[Δ/2 QuantumOperator["Z"] + SuperDagger[jcdLower] @ jcdA + jcdLower @ SuperDagger[jcdA],
       QuantumTensorProduct[atom, CoherentState[40][N[Sqrt[nb]]]], {t, 0, 30}], {2}]["DensityMatrix"]]},
   With[{λ = Clip[Re[Eigenvalues[ρA /. t -> #]], {10.^-16, 1}]}, -λ . Log2[λ]] &];
jcdPanel[nb_, Δ_] := With[{fs = jcdEntropy[nb, Δ] /@ jcdAtoms[nb, Δ]},
   Plot[{fs[[1]][s], fs[[2]][s], fs[[3]][s]}, {s, 0, 30}, PlotPoints -> 400, MaxRecursion -> 0, PlotRange -> {0, 1.02}, PlotStyle -> {Automatic, Thick, Dashed},
    PlotLegends -> {"|e\[RightAngleBracket]", Row[{"|", Subscript["+", "\[Alpha]"], "\[RightAngleBracket]"}], Row[{"|", Subscript["-", "\[Alpha]"], "\[RightAngleBracket]"}]},
    Frame -> True, FrameLabel -> {"g t", "S (bits)"}, PlotLabel -> Row[{OverBar["n"], " = ", nb, ",   \[CapitalDelta]/g = ", Δ}], ImageSize -> 300]];
jcdPanels = Outer[jcdPanel, {4, 16}, {0, 3}];
Legended[Grid[Map[First, jcdPanels, {2}]], Last[jcdPanels[[1, 1]]]]
```

An excited atom is strongly entangled within a few $1/g$; in the larger field its entropy dips near $gt=\pi\sqrt{\bar n}$, midway to the revival, almost to zero on resonance, before rising again. An atom started in a dressed state $|\pm_\alpha\rangle$ stays far less entangled throughout in the larger field, which acts more like the classical amplitude $\alpha$, but on resonance only briefly in the weaker one. Detuning lowers the excited atom's peaks but keeps its dip from reaching zero, and separates the two dressed starts, which coincide on resonance.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Where a two-mode squeezed vacuum is squeezed

Two field modes $a$ and $b$ start in the vacuum and are squeezed together, $|\psi\rangle=S_2(\xi)|0,0\rangle$ with $\xi=re^{i\theta}$. Where in the four-dimensional phase space is the state squeezed below the vacuum?

The two-mode squeezer creates and destroys photons in pairs, one in each mode,

$$
S_2(\xi)=\exp\left(\xi^*ab-\xi\,a^\dagger b^\dagger\right),
$$

and a phase $\theta$ only rotates mode $b$'s phase space, since $S_2(re^{i\theta})|0,0\rangle=e^{i\theta b^\dagger b}S_2(r)|0,0\rangle$, so $\xi=r$ here. With $a_1=a$, $a_2=b$, $x_k=(a_k+a_k^\dagger)/\sqrt2$ and $p_k=(a_k-a_k^\dagger)/i\sqrt2$ collected in $\mathbf r=(x_1,p_1,x_2,p_2)$, the covariance matrix is $\sigma_{kl}=\tfrac12\langle\{\Delta r_k,\Delta r_l\}\rangle$, equal to $\tfrac12\mathbb 1$ for the vacuum. $S_2$ is quadratic, so the state stays Gaussian: its Wigner function is $\propto e^{-\mathbf r^T\sigma^{-1}\mathbf r/2}$, and its marginal on any plane of two quadratures falls to $e^{-1/2}$ of its peak on the ellipse $\mathbf r^T\sigma^{-1}\mathbf r=1$ built from that plane's $2\times2$ block of $\sigma$.

QF builds $S_2(r)$ as a matrix exponential on 30 Fock levels per mode, applies it to $|0,0\rangle$, and takes $\sigma$ with `CovarianceMatrix`; each panel draws the ellipse in one plane, with the vacuum's circle dashed:

```wl
tmsA = AnnihilationOperator[30, {1}]; tmsB = AnnihilationOperator[30, {2}];
tmsSigma[r_] := tmsSigma[r] = CovarianceMatrix[MatrixExp[r (tmsA @ tmsB - SuperDagger[tmsA] @ SuperDagger[tmsB]), FockState[{0, 0}, 30]]];
tmsEllipse[σ_] := Table[MatrixPower[σ, 1/2] . {Cos[u], Sin[u]}, {u, 0, 2 Pi, Pi/60}];
tmsPanel[{i_, j_}] := With[{labels = {Subscript["x", 1], Subscript["p", 1], Subscript["x", 2], Subscript["p", 2]}},
   ListLinePlot[Append[tmsEllipse[tmsSigma[#][[{i, j}, {i, j}]]] & /@ {0.3, 0.6, 1.}, tmsEllipse[IdentityMatrix[2]/2]],
    PlotStyle -> {Automatic, Automatic, Automatic, Dashed}, PlotLegends -> {"r = 0.3", "r = 0.6", "r = 1", "vacuum"},
    AspectRatio -> 1, PlotRange -> {{-2.5, 2.5}, {-2.5, 2.5}}, Frame -> True, FrameLabel -> labels[[{i, j}]],
    PlotLabel -> Row[{"(", labels[[i]], ", ", labels[[j]], ") plane"}], ImageSize -> 280]];
tmsPanels = tmsPanel /@ {{1, 2}, {1, 3}, {2, 4}};
Legended[Row[First /@ tmsPanels, Spacer[10]], Last[First[tmsPanels]]]
```

On its own, mode 1 is a circle wider than the vacuum's, an isotropic Gaussian, so it is thermal light, noisier for larger $r$. The squeezing lives only in joint quadratures: $x_1$ and $x_2$ are anticorrelated and $p_1$ and $p_2$ correlated, so $x_1+x_2$ and $p_1-p_2$ fall inside the vacuum circle along the narrow axes of the tilted ellipses, further for larger $r$.
