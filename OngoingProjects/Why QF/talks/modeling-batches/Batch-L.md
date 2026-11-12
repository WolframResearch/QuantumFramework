```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Leakage of a fast π pulse on a transmon

A transmon is an anharmonic oscillator, frequency $\omega_q$ and anharmonicity $\alpha<0$, so its $|1\rangle\to|2\rangle$ transition sits only $|\alpha|$ below the qubit transition. A π pulse of length $\tau$, detuned by $\delta$ from the qubit, is meant to flip $|0\rangle$ to $|1\rangle$, and a DRAG term of strength $\lambda$ adds a second quadrature that follows the slope of the envelope. Where does the population end up after the pulse: outside the qubit, or left in $|0\rangle$?

With $\hbar=1$, the driven transmon and its form in the frame rotating at the drive frequency $\omega_d=\omega_q-\delta$, after the rotating-wave approximation, are

$$
H=\omega_q b^\dagger b+\frac{\alpha}{2}b^{\dagger2}b^2+\left[\Omega(t)\cos\omega_d t-\frac{\lambda\dot\Omega(t)}{\alpha}\sin\omega_d t\right](b+b^\dagger)\;\longrightarrow\;\frac{\alpha}{2}b^{\dagger2}b^2+\delta\,b^\dagger b+\frac{\Omega(t)}{2}(b+b^\dagger)-\frac{\lambda\dot\Omega(t)}{2\alpha}\,i(b^\dagger-b),
$$

with the raised-cosine envelope $\Omega(t)=\frac{\pi}{\tau}\left(1-\cos\frac{2\pi t}{\tau}\right)$, whose area is $\pi$, and the oscillator kept to its five lowest levels.

QF solves it from $|0\rangle$ and plots the leakage $P(n\ge2)$ (solid) and the population missing from $|1\rangle$, $1-P(1)$ (dashed), during the pulse. Times and frequencies are in units of $|\alpha|$:

```wl
dragB = AnnihilationOperator[5];
dragEnvelope[τ_] := Pi/τ (1 - Cos[2 Pi t/τ]);
dragH[α_, τ_, λ_, δ_] := α/2 SuperDagger[dragB] @ SuperDagger[dragB] @ dragB @ dragB + δ SuperDagger[dragB] @ dragB +
    dragEnvelope[τ]/2 (dragB + SuperDagger[dragB]) - λ D[dragEnvelope[τ], t]/(2 α) I (SuperDagger[dragB] - dragB);
dragPops[τ_, λ_, δ_] := QuantumEvolve[dragH[-1, τ, λ, δ], FockState[0, 5], {t, 0, τ}][t]["ProbabilitiesList"];
dragλs = {0, 0.5, 1};
dragPanel[τ_, δ_] := With[{pops = dragPops[τ, #, δ] & /@ dragλs},
   Show[LogPlot[Evaluate[Total[#[[3 ;;]]] & /@ pops], {t, 0, τ}, PlotRange -> {10^-9, 2}, PlotLegends -> (Row[{"λ = ", #}] & /@ dragλs)],
    LogPlot[Evaluate[1 - #[[2]] & /@ pops], {t, 0, τ}, PlotRange -> {10^-9, 2}, PlotStyle -> Dashed], Frame -> True, FrameLabel -> {"|α| t", "population"},
    PlotLabel -> Row[{"|α|τ = ", τ, ",   δ/|α| = ", δ}], ImageSize -> 300]];
dragPanels = Outer[dragPanel, {5, 10, 20}, {-0.05, 0, 0.05}];
Legended[Grid[Map[First, dragPanels, {2}]], Last[dragPanels[[1, 1]]]]
```

The drive couples only neighboring levels, so leakage builds up through $|1\rangle$ during the pulse; only the shortest pulse without DRAG keeps leaking to the end, and otherwise it peaks and partly returns. At $\lambda=1$ the population left outside the qubit drops by orders of magnitude, and a longer pulse lowers it further. The rest of the error is population left in $|0\rangle$: at $\lambda=1$ the dashed curve lies far above the solid one, depends strongly on $\delta$, and the best $\delta$ differs between $\lambda=0$ and $\lambda=1$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The CZ gate from two coupled transmons

Two transmons, frequencies $\omega_j$ and anharmonicities $\alpha_j$ (equal, $\alpha_1=\alpha_2=\alpha$, in the figure), exchange excitations at coupling $g$. A flux pulse tunes transmon 1 until $|11\rangle$ is detuned by $\delta$ from $|02\rangle$, holds it there for a time $t$, and tunes it back, both tunings taken as instant. What does $|11\rangle$ come back with?

With $\hbar=1$ and three levels per transmon, the pair Hamiltonian is

$$
H=\sum_{j=1,2}\left[\omega_j\,n_j+\frac{\alpha_j}{2}\,n_j(n_j-1)\right]+g\left(b_1^\dagger b_2+b_2^\dagger b_1\right),\qquad n_j=b_j^\dagger b_j,\qquad \delta=\omega_1-\omega_2-\alpha_2 .
$$

It conserves the total excitation number, so $|11\rangle$ couples only to $|02\rangle$ and $|20\rangle$; $|20\rangle$ lies about $|\alpha_1+\alpha_2|$ away, and dropping it leaves a $2\times2$ block on $\{|11\rangle,|02\rangle\}$. The gate acts on $|11\rangle$ through the conditional amplitude $c=U_{11,11}\,e^{-i\arg(U_{01,01}U_{10,10}/U_{00,00})}$, which removes the phases each transmon picks up alone; with block energies counted from $|02\rangle$, the block gives $c=e^{i\delta t}\,(e^{-iH_{\rm block}t})_{11,11}$.

WL derives the $\{|11\rangle,|02\rangle\}$ block symbolically and traces $c$ in the complex plane over one full swap, $0\le t\le 2\pi/\sqrt{\delta^2+8g^2}$, nine-level (solid) and block (dashed), with $\omega_2=0$; rates are in units of $g$:

```wl
czLower = DiagonalMatrix[Sqrt[{1, 2}], 1];
{czB1, czB2} = {KroneckerProduct[czLower, IdentityMatrix[3]], KroneckerProduct[IdentityMatrix[3], czLower]};
czDuffing[bj_, w_, a_] := With[{nj = Transpose[bj] . bj}, w nj + a/2 nj . (nj - IdentityMatrix[9])];
czH = czDuffing[czB1, ω1, α1] + czDuffing[czB2, ω2, α2] + g (Transpose[czB1] . czB2 + Transpose[czB2] . czB1);
czPos = 1 + FromDigits[#, 3] & /@ {{1, 1}, {0, 2}};
czBlock = Simplify[czH[[czPos, czPos]] - czH[[czPos[[2]], czPos[[2]]]] IdentityMatrix[2] /. ω1 -> ω2 + α2 + δ];
czBlockAmp = FullSimplify[Exp[I δ t] MatrixExp[-I czBlock t][[1, 1]]];
czFullAmp[d_, a_, s_] := With[{u = MatrixExp[-I N[czH /. {ω2 -> 0, α1 -> a, α2 -> a, ω1 -> a + d, g -> 1}] s]},
    u[[5, 5]] Exp[-I Arg[u[[2, 2]] u[[4, 4]]/u[[1, 1]]]]];
czαs = {-3, -10, -40};
czPanel[d_] := With[{T = N[2 Pi/Sqrt[d^2 + 8]]},
   With[{curves = Append[Table[ReIm[czFullAmp[d, a, s]], {a, czαs}, {s, 0, T, T/150}], Table[ReIm[czBlockAmp /. {δ -> d, g -> 1, t -> s}], {s, 0, T, T/150}]]},
    ListLinePlot[curves, PlotStyle -> {Automatic, Automatic, Automatic, Dashed}, AspectRatio -> 1, PlotRange -> {{-1.05, 1.05}, {-1.05, 1.05}},
     Frame -> True, FrameLabel -> {"Re c", "Im c"},
     PlotLegends -> Append[(Row[{"α/g = ", #}] & /@ czαs), "2×2 block"],
     PlotLabel -> Row[{"δ/g = ", d}], ImageSize -> 300]]];
czPanels = czPanel /@ {0, 1, 3};
Legended[Grid[{Prepend[First /@ czPanels, Labeled[MatrixForm[czBlock], "{|11⟩, |02⟩} block", Top]]}], Last[First[czPanels]]]
```

The derived block couples $|11\rangle$ to $|02\rangle$ at $\sqrt2g$, since $b^\dagger|1\rangle=\sqrt2|2\rangle$, and at the crossing its $c$ runs along the real axis from 1 through 0 to $-1$ at $t=\pi/(\sqrt2g)$: a controlled-Z. The nine-level evolution follows that path more closely as $|\alpha|/g$ grows and misses $-1$ at $\alpha/g=-3$. Away from the crossing the swap is incomplete, and the block brings $c$ back to the unit circle at a phase other than $\pi$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Pulse area sets the rotation

A qubit is driven for a time $T$ by a pulse of shape $\varepsilon(t)$, peak 1, amplitude $a$ and peak Rabi frequency $\Omega_{\max}$, at drive frequency $\omega_d$ detuned by $\delta=\omega_d-\omega_q$ from the qubit; it relaxes at $1/T_1$ and dephases at total rate $1/T_2$. Does only the pulse area set the rotation, for any shape and detuning?

With $\hbar=1$ and $\sigma_z|0\rangle=|0\rangle$, the qubit in the frame rotating at $\omega_d$ obeys

$$
\dot\rho=-i[H,\rho]+\frac{1}{T_1}\mathcal D[\sigma_-]\rho+\frac{\gamma_\phi}{2}\mathcal D[\sigma_z]\rho,\qquad H=\frac{\delta}{2}\,\sigma_z+\frac{\Omega_{\max}\,a\,\varepsilon(t)}{2}\,\sigma_x,\qquad \gamma_\phi=\frac1{T_2}-\frac1{2T_1}.
$$

On resonance $H$ commutes with itself at all times and rotates the qubit by the pulse area $\theta=\Omega_{\max}a\int_0^T\varepsilon\,dt$, so without decay $P(1)=\sin^2(\theta/2)$.

TransmonLab sweeps $a$ for three pulse shapes and plots $P(1)$ at the end of the pulse, where a readout would start, against $\theta$, with $\sin^2(\theta/2)$ dashed. The modeled qubit has $T=100$ ns, $\Omega_{\max}T=5\pi$ and $T_1=T_2=10$ μs:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/mbahram/TransmonLab/TransmonLab"]; Needs["TransmonLab`"];
rabiSystem = TransmonSystem["T1" -> 10.*^-6, "T2" -> 10.*^-6, "OmegaMax" -> 2 Pi 25*^6];
rabiShapes = {"Gaussian", "Cosine", "Square"};
rabiCurve[shape_, δT_] := With[{p = Pulse[shape, 100]}, With[{θa = 2 Pi 25*^6 p["EnvelopeArea"] 10^-9},
    With[{amps = Subdivide[0, 2 Pi/θa, 40]},
     Transpose[{θa amps, TransmonPhysics[rabiSystem, PulseSequence[Drive[p, "Qubit"], Sync["Qubit", "Resonator"], Readout[Pulse["FlatTop", 200, "RiseTime" -> 10], "Resonator"]],
         Sweep["Amplitude", amps], "DetuningQubit" -> δT/100*^-9, "ShowProgress" -> False]["Populations"]["1"]}]]]];
rabiPanel[δT_] := Show[ListLinePlot[rabiCurve[#, δT] & /@ rabiShapes, PlotLegends -> rabiShapes],
   Plot[Sin[θ/2]^2, {θ, 0, 2 Pi}, PlotStyle -> Dashed], PlotRange -> {{0, 2 Pi}, {0, 1}}, Frame -> True, FrameLabel -> {"θ", "P(1)"},
   PlotLabel -> Row[{"δ T = ", δT}], ImageSize -> 300];
rabiPanels = rabiPanel /@ {0, 1, 3};
Legended[Grid[{First /@ rabiPanels}], Last[First[rabiPanels]]]
```

On resonance the three shapes fall on one curve: only the area matters, so a single constant per shape turns any target angle into an amplitude. Detuning breaks this, and differently for each shape, since the state then depends on how the area is spread in time.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Dispersive readout pointer states

A readout cavity of linewidth $\kappa$, coupled to a transmon, has its frequency pulled by $+\chi$ or $-\chi$ depending on the qubit state, and is driven by a readout pulse $\varepsilon_r(t)$ of peak amplitude $\epsilon_{\max}$. How does the cavity field separate for the two qubit states as $\kappa/2|\chi|$ changes?

With $\hbar=1$, in the frame of the readout drive, the cavity holds a coherent state $|\alpha_j\rangle$ for qubit state $j$, with

$$
\dot\alpha_{0}=-\left(\frac{\kappa}{2}+i\chi\right)\alpha_0-\frac{i}{2}\epsilon_{\max}\varepsilon_r(t),\qquad \dot\alpha_{1}=-\left(\frac{\kappa}{2}-i\chi\right)\alpha_1-\frac{i}{2}\epsilon_{\max}\varepsilon_r(t),
$$

and the measured IQ point is the time integral of the output field, proportional to $\sqrt\kappa\,\alpha_j$.

TransmonLab solves both equations for a 2 μs flat-top readout and plots the intracavity field $\alpha_j(t)$ for three linewidths. The modeled system is TransmonLab's default transmon (qubit 5.1 GHz, anharmonicity $-220$ MHz, cavity 6 GHz, $\epsilon_{\max}/2\pi=5$ MHz) with qubit-cavity coupling 50 MHz, $T_1=T_2=10$ μs, and $\chi$ from the dispersive formula:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/mbahram/TransmonLab/TransmonLab"]; Needs["TransmonLab`"];
pointerSystem[κHz_] := TransmonSystem["Coupling" -> 50*^6, "CavityLinewidth" -> κHz, "T1" -> 10.*^-6, "T2" -> 10.*^-6];
pointerChiHz = Abs[pointerSystem[1*^6]["Chi"]];
pointerPanel[r_] := With[{dr = DispersiveReadout[pointerSystem[2 r pointerChiHz], Pulse["FlatTop", 2000, "RiseTime" -> 50]]},
   ParametricPlot[Evaluate[ReIm[{dr["Alpha0"][s], dr["Alpha1"][s]}]], {s, 0, dr["ReadoutDuration"]}, PlotRange -> All, AspectRatio -> 1,
    Frame -> True, FrameLabel -> {"Re α", "Im α"}, PlotLegends -> {"qubit in |0⟩", "qubit in |1⟩"},
    PlotLabel -> Row[{"κ/2|χ| = ", r}], ImageSize -> 300]];
pointerPanels = pointerPanel /@ {0.5, 1, 2};
Legended[Grid[{First /@ pointerPanels}], Last[First[pointerPanels]]]
```

The cavity rings up along a different spiral for each qubit state, toward two mirror-image points; a narrower cavity puts them farther apart and spirals more before settling.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A drive pulse as an expression, swept by substitution

A DRAG pulse has a Gaussian in-phase envelope $I(t)=e^{-(t-t_0)^2/2\sigma^2}$, here $\sigma=5$ ns and $t_0=20$ ns over 40 ns, and a quadrature $Q(t)=-\beta\,\dot I(t)$ whose coefficient $\beta$ is found by sweeping it. The program drives a qubit at 5.105 GHz and reads it through a resonator at 6.012 GHz; the controller plays sampled arrays at 1 ns. What does the compiled program play at each $\beta$?

QUALink keeps the pulse an expression in $t$ and turns the sweep $\beta/\sigma\in\{-1,0,1\}$ into a switch over three compiled pulses; the program itself is an association that `Cases` reads. Nothing here touches a device.

The cell compiles the sweep, finds the switch with `Cases`, and plots the samples each of its pulses plays (points) with the expression $\{I,Q\}$ at that $\beta$ (dashed):

```wl
PacletDirectoryLoad[{"/Users/mohammadb/Documents/GitHub/mbahram/TransmonLab/TransmonLab", "/Users/mohammadb/Documents/GitHub/mbahram/QUALink/QUALink"}];
Needs["WolframQuantumComputation`QUALink`"];
quaPulse = {Exp[-(t - 20)^2/50], -β D[Exp[-(t - 20)^2/50], t]};
quaProg = QUAProgram["Qubit" -> {5.105*^9, quaPulse, {t, 0, 40}}, "Resonator" -> 6.012*^9, "β" -> {-5, 0, 5}];
quaSamples[pulse_] := quaProg["Waveforms"] /@ Values[quaProg["Config"]["pulses"][pulse]["waveforms"]];
quaPanel[pulse_, b_] := Show[Plot[Evaluate[quaPulse /. β -> b], {t, 0, 40}, PlotStyle -> Directive[Dashed, Thick]],
   ListPlot[quaSamples[pulse], DataRange -> {0, 39}, PlotMarkers -> {Automatic, Small}, PlotLegends -> {"I samples", "Q samples"}], PlotRange -> {-1.05, 1.05}, Frame -> True,
   FrameLabel -> {"t (ns)", "amplitude"}, PlotLabel -> Row[{"β/σ = ", b/5}], ImageSize -> 300];
quaSwitch = First[Cases[quaProg["IR"], KeyValuePattern["Head" -> "SwitchPlay"], Infinity]];
quaPanels = MapThread[quaPanel, {quaSwitch["Pulses"], quaProg["SweepValues"]}];
Legended[Grid[{First /@ quaPanels}], Last[First[quaPanels]]]
```

Each sweep value is one substitution into the expression, and the samples the program plays lie on that expression; the quadrature changes sign with $\beta$ and vanishes at $\beta=0$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The transmon spectrum against offset charge

A Cooper-pair box is a superconducting island, charging energy $E_C$, joined to a reservoir by a Josephson junction of energy $E_J$; a gate sets the offset charge $n_g$, in units of $2e$. How do its lowest levels depend on $n_g$ as $E_J/E_C$ grows?

With $n$ the number of Cooper pairs on the island and $\varphi$ the conjugate phase,

$$
H=4E_C\,(n-n_g)^2-E_J\cos\varphi,\qquad \cos\varphi=\frac12\sum_n\big(|n\rangle\langle n+1|+|n+1\rangle\langle n|\big),
$$

so in the charge basis $H$ is tridiagonal: $4E_C(n-n_g)^2$ on the diagonal and $-E_J/2$ beside it.

WL diagonalizes it for $|n|\le15$ and plots the three lowest levels $E_m(n_g)$, measured from $E_0$ at $n_g=1/2$ and in units of the qubit splitting $E_{01}$ there; the matrix is built with $E_C=1$:

```wl
cpbLevels[ej_, ng_?NumericQ] := Take[Sort @ Eigenvalues[N[DiagonalMatrix[4 (Range[-15, 15] - ng)^2] - ej/2 (DiagonalMatrix[ConstantArray[1, 30], 1] + DiagonalMatrix[ConstantArray[1, 30], -1])]], 3];
cpbLevel[ej_, m_, ng_?NumericQ] := cpbLevels[ej, ng][[m + 1]];
cpbPanel[ej_] := With[{ref = cpbLevels[ej, 1/2]},
   Plot[Evaluate[Table[(cpbLevel[ej, m, ng] - ref[[1]])/(ref[[2]] - ref[[1]]), {m, 0, 2}]], {ng, -2, 2}, PlotRange -> All, Frame -> True,
    FrameLabel -> {Subscript["n", "g"], Row[{"(", Subscript["E", "m"], " - ", Subscript["E", "0"], ")/", Subscript["E", "01"]}]},
    PlotLegends -> (Row[{"m = ", #}] & /@ Range[0, 2]), PlotLabel -> Row[{Subscript["E", "J"], "/", Subscript["E", "C"], " = ", ej}], ImageSize -> 260]];
cpbPanels = cpbPanel /@ {1, 5, 10, 50};
Legended[Grid[{First /@ cpbPanels}], Last[First[cpbPanels]]]
```

At $E_J/E_C=1$ the levels are charge parabolas that the junction splits at $n_g=1/2$, and the qubit frequency swings strongly with offset charge. As $E_J/E_C$ grows the bands flatten, and by $E_J/E_C=50$ every level is nearly independent of $n_g$, which is the transmon. From $E_J/E_C=10$ on the second level lies below twice the first at every $n_g$, so the qubit transition can be addressed on its own.
