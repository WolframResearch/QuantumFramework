```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Relaxation of every qubit state, and why T2 is at most 2T1

A qubit with level splitting $\omega$ decays from $|e\rangle$ to $|g\rangle$ at rate $\gamma_1$ and is dephased by the jump operator $\sqrt{\gamma_\phi}\,\sigma_z$. How long can its coherence outlive its energy?

With $\hbar=1$, $|e\rangle=|0\rangle$, $|g\rangle=|1\rangle$ and $\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, the Lindblad equation is

$$
\dot\rho=-i\left[\tfrac{\omega}{2}\sigma_z,\rho\right]+\gamma_1\,\mathcal D[\sigma_-]\rho+\gamma_\phi\,\mathcal D[\sigma_z]\rho,\qquad \mathcal D[L]\rho=L\rho L^\dagger-\tfrac12\{L^\dagger L,\rho\}.
$$

QF solves it once, with a symbolic initial density matrix $\rho_0$ and symbolic rates, reads off $1/T_1$ and $1/T_2$, and plots $\langle\sigma_x\rangle$ and $\langle\sigma_z\rangle$ from Bloch states tilted by $\theta_0$ away from $|e\rangle$; time is in units of $1/\gamma_1$ and $\omega/\gamma_1=6$:

```wl
relaxA1 = FullSimplify[Normal[QuantumEvolve[ω/2 QuantumOperator["Z"], {QuantumOperator["-"], QuantumOperator["Z"]} -> {γ1, γφ},
     QuantumState[Array[ρ0, {2, 2}]], t]["DensityMatrix"]], γ1 > 0 && γφ >= 0 && t >= 0 && ω ∈ Reals];
{rateT1, rateT2} = Simplify[-Re[D[Log[#], t]], γ1 > 0 && γφ >= 0 && ω ∈ Reals] & /@ {relaxA1[[1, 1]], relaxA1[[1, 2]]};
tiltA1[θ0_] := Thread[Flatten[Array[ρ0, {2, 2}]] -> Flatten[{{1 + Cos[θ0], Sin[θ0]}, {Sin[θ0], 1 - Cos[θ0]}}/2]];
{γφsA1, θ0sA1} = {{0, 0.25, 1}, {π/4, π/2, 3 π/4}};
panelA1[γφv_, θ0_] := Plot[Evaluate[Re[{Tr[relaxA1 . PauliMatrix[1]], Tr[relaxA1 . PauliMatrix[3]]} /. tiltA1[θ0] /. {γ1 -> 1, γφ -> γφv, ω -> 6}]
     ~Append~ (Sin[θ0] Exp[-t/2])], {t, 0, 6}, PlotRange -> {-1.05, 1.05}, PlotStyle -> {Automatic, Automatic, Dashed},
    PlotLegends -> {Row[{"\[LeftAngleBracket]", Subscript["σ", "x"], "\[RightAngleBracket]"}], Row[{"\[LeftAngleBracket]", Subscript["σ", "z"], "\[RightAngleBracket]"}],
      Row[{"sin ", Subscript["θ", 0], " exp(-", Subscript["γ", 1], "t/2)"}]},
    Frame -> True, FrameLabel -> {Row[{Subscript["γ", 1], " t"}], None},
    PlotLabel -> Row[{Subscript["γ", "φ"], "/", Subscript["γ", 1], " = ", γφv, ",   ", Subscript["θ", 0], " = ", θ0}], ImageSize -> 300];
panelsA1 = Outer[panelA1, γφsA1, θ0sA1];
Labeled[Legended[Grid[Map[First, panelsA1, {2}]], Last[panelsA1[[1, 1]]]],
  Row[{"1/", Subscript["T", 1], " = ", TraditionalForm[rateT1], ",     1/", Subscript["T", 2], " = ", TraditionalForm[rateT2]} /. {γ1 -> Subscript["γ", 1], γφ -> Subscript["γ", "φ"]}], Top]
```

Both rates come out the same for every initial state, and energy decay alone already costs the coherence half of $\gamma_1$, so $T_2\le 2T_1$, with equality only when $\gamma_\phi=0$. In the figure the precessing $\langle\sigma_x\rangle$ fills the dashed envelope only in the top row, and dephasing pulls it inside, while $\langle\sigma_z\rangle$ relaxes toward $-1$ the same way in every row.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Entanglement sudden death

Two qubits A and B start in $\cos\theta\,|gg\rangle+\sin\theta\,|ee\rangle$, and each decays on its own, A at rate $\gamma_A$ and B at rate $\gamma_B$. Each qubit's excitation decays exponentially; does the entanglement also fade only gradually?

With $\hbar=1$, $|e\rangle=|0\rangle$, $|g\rangle=|1\rangle$ and qubit frequencies $\omega_A$, $\omega_B$, the Lindblad equation is

$$
\dot\rho=-i\left[\tfrac{\omega_A}{2}\sigma_z^A+\tfrac{\omega_B}{2}\sigma_z^B,\rho\right]+\gamma_A\,\mathcal D[\sigma_-^A]\rho+\gamma_B\,\mathcal D[\sigma_-^B]\rho,
$$

and in the frame rotating with each qubit the Hamiltonian drops out while the two dissipators stay as they are.

QF solves it once with $\theta$, both rates and $t$ symbolic, then evaluates QF's concurrence of the solution; time is in units of $1/\gamma_A$:

```wl
pairA2 = Normal[QuantumEvolve[QuantumOperator[0 IdentityMatrix[4]], {QuantumOperator["-", {1}], QuantumOperator["-", {2}]} -> {γA, γB},
     QuantumState[{Sin[θ], 0, 0, Cos[θ]}], t]["DensityMatrix"]];
{θsA2, γBsA2} = {{π/6, π/4, π/3, 5 π/12}, {1, 0.3, 0}};
concA2[θv_, γBv_] := Table[{s, QuantumEntanglementMonotone[QuantumState[pairA2 /. {θ -> θv, γA -> 1, γB -> γBv, t -> s}], {{1}, {2}}, "Concurrence"]}, {s, 0., 4., 0.025}];
panelA2[γBv_] := ListLinePlot[concA2[#, γBv] & /@ θsA2, PlotRange -> {0, 1.02}, PlotStyle -> {Automatic, Automatic, Dashed, Automatic}, PlotLegends -> (Row[{"θ = ", #}] & /@ θsA2),
    Frame -> True, FrameLabel -> {Row[{Subscript["γ", "A"], " t"}], "concurrence"}, PlotLabel -> Row[{Subscript["γ", "B"], "/", Subscript["γ", "A"], " = ", γBv}], ImageSize -> 300];
panelsA2 = panelA2 /@ γBsA2;
Legended[Row[First /@ panelsA2], Last[First[panelsA2]]]
```

When both qubits decay, the states with $\theta>\pi/4$, where $|ee\rangle$ outweighs $|gg\rangle$, lose all their entanglement at a finite time and never regain it, while at $\theta\le\pi/4$ it only fades. A slower qubit B delays the death, and if B does not decay at all the entanglement never dies; there $\theta$ and $\pi/2-\theta$ give the same curve.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Two atoms in a shared field: superradiance, subradiance and a dark state

Two identical atoms, frequency $\omega$, sit close enough to radiate into the same field modes, so each decays at rate $\gamma$ and the field correlates their emission with strength $0\le r\le1$. Which states radiate faster or slower than a lone atom, and can one stop radiating altogether?

With $\hbar=1$, $|e\rangle=|0\rangle$ and $|g\rangle=|1\rangle$, the Lindblad equation in Kossakowski form, with rate matrix $\gamma\begin{pmatrix}1&r\\ r&1\end{pmatrix}$, is

$$
\dot\rho=-i\left[\tfrac{\omega}{2}\left(\sigma_z^1+\sigma_z^2\right),\rho\right]+\gamma\sum_{i,j=1}^{2}\begin{pmatrix}1&r\\ r&1\end{pmatrix}_{ij}\left(\sigma_-^i\rho\,\sigma_+^j-\tfrac12\left\{\sigma_+^j\sigma_-^i,\rho\right\}\right),
$$

and in the frame rotating at $\omega$ the Hamiltonian drops out; the field's coherent exchange between the atoms commutes with $N$ and is left out.

QF takes the rate matrix as it stands and solves for the mean number of excited atoms $\langle N\rangle$ from four initial states, with $|T\rangle,|S\rangle=(|eg\rangle\pm|ge\rangle)/\sqrt2$; time is in units of $1/\gamma$:

```wl
statesA3 = <|"|ee\[RightAngleBracket]" -> QuantumState["00"], "|T\[RightAngleBracket]" -> QuantumState[{0, 1, 1, 0}/Sqrt[2]],
    "|S\[RightAngleBracket]" -> QuantumState[{0, 1, -1, 0}/Sqrt[2]], "|eg\[RightAngleBracket]" -> QuantumState["01"]|>;
excitationA3[r_, ψ_] := Re @ Tr[Normal[QuantumEvolve[QuantumOperator[0 IdentityMatrix[4]], {QuantumOperator["-", {1}], QuantumOperator["-", {2}]} -> {{1, r}, {r, 1}},
      ψ, t]["DensityMatrix"]] . DiagonalMatrix[{2, 1, 1, 0}]];
panelA3[r_] := Plot[Evaluate[excitationA3[r, #] & /@ Values[statesA3]], {t, 0, 5}, PlotRange -> {0, 2.05}, PlotStyle -> {Automatic, Automatic, Automatic, Dashed}, PlotLegends -> Keys[statesA3],
    Frame -> True, FrameLabel -> {"γ t", "\[LeftAngleBracket]N\[RightAngleBracket]"}, PlotLabel -> Row[{"r = ", r}], ImageSize -> 260];
panelsA3 = panelA3 /@ {0, 1/2, 9/10, 1};
Legended[Grid[Partition[First /@ panelsA3, 2]], Last[First[panelsA3]]]
```

Without correlation, at $r=0$, the atoms decay independently and $|T\rangle$, $|S\rangle$ and $|eg\rangle$ follow one curve. Correlation speeds up the symmetric $|T\rangle$ (superradiance) and slows down the antisymmetric $|S\rangle$ (subradiance); at $r=1$ the singlet stops radiating, so $|eg\rangle$, half singlet, keeps half an excitation, while $|ee\rangle$ still empties completely.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Master equation, Kraus channels and a circuit: one Ramsey fringe

A Ramsey experiment: a Hadamard gate, a wait of time $t$ during which the qubit, detuned by $\delta$ from the drive, decays at rate $\gamma_1$ and dephases through the jump operator $\sqrt{\gamma_\phi}\,\sigma_z$, and a second Hadamard. Is the master equation for the wait the same as a circuit of textbook channels?

With $\hbar=1$, $|g\rangle=|0\rangle$, $|e\rangle=|1\rangle$, $\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, qubit frequency $\omega$ and drive frequency $\omega_d=\omega-\delta$, the wait obeys, in the frame rotating at $\omega_d$,

$$
\dot\rho=-i\left[\tfrac{\delta}{2}\sigma_z,\rho\right]+\gamma_1\,\mathcal D[\sigma_-]\rho+\gamma_\phi\,\mathcal D[\sigma_z]\rho .
$$

The circuit replaces the wait by the rotation $e^{-i\delta t\sigma_z/2}$, the amplitude-damping channel with Kraus operators $\mathrm{diag}(1,\sqrt{1-\lambda})$ and $\sqrt\lambda\,|g\rangle\langle e|$ at $\lambda=1-e^{-\gamma_1t}$, and the phase-damping channel $\mathrm{diag}(1,\sqrt{1-\mu})$, $\sqrt\mu\,|e\rangle\langle e|$ at $\mu=1-e^{-4\gamma_\phi t}$.

QF computes the final $P(e)$ from $|g\rangle$ both ways with every parameter symbolic, tests them for equality, and plots the master equation as lines and the circuit as dots; time is in units of $1/\gamma_1$:

```wl
hadA4 = QuantumOperator["H"];
viaMEA4 = Normal[hadA4[QuantumEvolve[-δ/2 QuantumOperator["Z"], {QuantumOperator[{{0, 1}, {0, 0}}], QuantumOperator["Z"]} -> {γ1, γφ},
       hadA4[QuantumState["0"]], t]]["DensityMatrix"]][[2, 2]];
viaCircuitA4 = Normal[(hadA4 @ QuantumChannel["PhaseDamping"[1 - E^(-4 γφ t)]] @ QuantumChannel["AmplitudeDamping"[1 - E^(-γ1 t)]] @
       QuantumOperator["RZ"[-δ t]] @ hadA4)[QuantumState["0"]]["DensityMatrix"]][[2, 2]];
assumeA4 = γ1 > 0 && γφ >= 0 && δ ∈ Reals && t >= 0;
panelA4[γφv_, δv_] := With[{rules = {γ1 -> 1, γφ -> γφv, δ -> δv}}, Show[
    Plot[Evaluate[Re[viaMEA4 /. rules]], {t, 0, 4}, PlotRange -> All],
    ListPlot[Table[{s, Re[viaCircuitA4 /. rules /. t -> s]}, {s, 0, 4, 0.1}]],
    PlotRange -> {0, 1}, Frame -> True, FrameLabel -> {Row[{Subscript["γ", 1], " t"}], "P(e)"},
    PlotLabel -> Row[{Subscript["γ", "φ"], "/", Subscript["γ", 1], " = ", γφv, ",   δ/", Subscript["γ", 1], " = ", δv}], ImageSize -> 300]];
Labeled[Grid[Outer[panelA4, {0, 0.5}, {0, 3, 8}]],
  Row[{"P(e) = ", TraditionalForm[FullSimplify[ComplexExpand[viaMEA4], assumeA4] /. {γ1 -> Subscript["γ", 1], γφ -> Subscript["γ", "φ"]}], ",     equal for all parameters: ", FullSimplify[viaMEA4 == viaCircuitA4, assumeA4]}], Top]
```

The two routes agree for every value of every parameter: integrating the master equation for a time $t$ is the circuit of a rotation and two named channels. The fringe oscillates about $P(e)=1/2$ at the detuning, and $\gamma_1$ enters its decay only through $1/T_2$; at $\delta=0$ there is no fringe, only that relaxation.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The quantum Zeno crossover as an exceptional point

A qubit is driven with Rabi frequency $\Omega$, detuned by $\Delta$ from the drive, while a detector measures $\sigma_z$ at strength $k$. Weak watching damps the Rabi oscillation and strong watching freezes it; what happens to the dynamics at the crossover?

With $\hbar=1$, $|e\rangle=|0\rangle$, $|g\rangle=|1\rangle$, $\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$, qubit frequency $\omega$ and drive frequency $\omega-\Delta$, the Hamiltonian in the frame rotating with the drive is $H=\frac{\Delta}{2}\sigma_z+\frac{\Omega}{2}\sigma_x$, and averaging over the detector record gives

$$
\dot\rho=\mathcal L\rho=-i[H,\rho]+k\,\mathcal D[\sigma_z]\rho .
$$

Away from exceptional points, every solution is a sum of $e^{\lambda t}$ over the eigenvalues $\lambda$ of the superoperator $\mathcal L$. At an exceptional point two eigenvalues and their eigenvectors coincide, so $\mathcal L$ cannot be diagonalized, and the largest overlap $|\langle v_i|v_j\rangle|/(\|v_i\|\|v_j\|)$ between eigenvectors of the decaying eigenvalues (the steady state's $\lambda=0$ left out) reaches 1.

QF builds $\mathcal L$ with $\Delta$, $\Omega$ and $k$ symbolic, then plots its eigenvalues (top) and that overlap (bottom) as $k/\Omega$ runs from 0 to 3, in units of $\Omega$:

```wl
liouA5 = Normal[QuantumOperator["Liouvillian"[Δ/2 QuantumOperator["Z"] + Ω/2 QuantumOperator["X"], {QuantumOperator["Z"]}, {k}]]["Matrix"]];
specA5[Δv_, kv_] := Eigensystem[N[liouA5 /. {Ω -> 1, Δ -> Δv, k -> kv}]];
overlapA5[Δv_, kv_] := Max[Abs[Conjugate[#1] . #2]/(Norm[#1] Norm[#2]) & @@@ Subsets[Most[Last[specA5[Δv, kv]]], {2}]];
ksA5 = {0.5, 1, 2};
eigenPanelA5[Δv_] := ListPlot[Prepend[ReIm[Most[First[specA5[Δv, #]]]] & /@ ksA5, Catenate[ReIm[Most[First[specA5[Δv, #]]]] & /@ Range[0.02, 3, 0.02]]],
    PlotStyle -> {PointSize[0.006], PointSize[0.02], PointSize[0.02], PointSize[0.02]}, PlotLegends -> Prepend[Row[{"k/Ω = ", #}] & /@ ksA5, "0 < k/Ω < 3"],
    PlotRange -> All, Frame -> True, FrameLabel -> {"Re λ/Ω", "Im λ/Ω"}, PlotLabel -> Row[{"Δ/Ω = ", Δv}], AspectRatio -> 1, ImageSize -> 280];
overlapPanelA5[Δv_] := ListLinePlot[Table[{kv, overlapA5[Δv, kv]}, {kv, 0.02, 3, 0.01}], PlotRange -> {0, 1.02},
    Frame -> True, FrameLabel -> {"k/Ω", "eigenvector overlap"}, PlotLabel -> Row[{"Δ/Ω = ", Δv}], ImageSize -> 280];
eigenPanelsA5 = eigenPanelA5 /@ {0, 0.25, 1};
Legended[Grid[{First /@ eigenPanelsA5, overlapPanelA5 /@ {0, 0.25, 1}}], Last[First[eigenPanelsA5]]]
```

On resonance the damped Rabi pair falls onto the real axis at $k=\Omega$, where its two eigenvectors coincide: the crossover is an exceptional point, and past it one real rate moves toward zero, the slowing that is the Zeno effect. A small detuning splits it into two exceptional points, one where the Rabi pair turns real and one where two real rates meet and turn complex again. At $\Delta=\Omega$ no two eigenvectors coincide anywhere in this range.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Mollow triplet and squeezed resonance fluorescence

A two-level atom, frequency $\omega_a$, decays at rate $\gamma$ and is driven by a laser at frequency $\omega_L$ with Rabi frequency $\Omega$. What does the light it scatters look like, to a spectrometer and to a homodyne detector?

With $\hbar=1$, $|e\rangle=|0\rangle$, $|g\rangle=|1\rangle$, $\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$ and $\Delta=\omega_a-\omega_L$, the Hamiltonian in the frame rotating at $\omega_L$ (rotating-wave approximation) is $H=\frac{\Delta}{2}\sigma_z+\frac{\Omega}{2}\sigma_x$, and

$$
\dot\rho=\mathcal L\rho=-i[H,\rho]+\gamma\,\mathcal D[\sigma_-]\rho .
$$

By the quantum regression theorem, two-time correlations evolve with the same $\mathcal L$ from the steady state $\rho_{ss}$, so at frequency $\nu$ from the laser the incoherent fluorescence spectrum and the homodyne spectrum of the quadrature $x=c+c^\dagger$ of the emitted field $c=-i\sqrt\gamma\,\sigma_-$ (out of phase with the drive), normalized to shot noise, are

$$
S(\nu)=2\,\mathrm{Re}\,\mathrm{Tr}\!\left[\sigma_+(i\nu-\mathcal L)^{-1}\delta\sigma_-\rho_{ss}\right],\qquad
S_x(\nu)=1+2\,\mathrm{Re}\,\mathrm{Tr}\!\left[x\,(i\nu-\mathcal L)^{-1}\left(c\rho_{ss}+\rho_{ss}c^\dagger-\langle x\rangle\rho_{ss}\right)\right],
$$

with $\delta\sigma_-=\sigma_--\langle\sigma_-\rangle$.

QF builds $\mathcal L$, its null vector gives $\rho_{ss}$, and each frequency costs one linear solve; each fluorescence spectrum is scaled to its own maximum, and rates are in units of $\gamma$:

```wl
specA6[L_, seed_, readout_] := Table[{ν, 2 Re @ Tr[readout . ArrayReshape[LinearSolve[I ν IdentityMatrix[4] - L, Flatten[seed]], {2, 2}]]}, {ν, Subdivide[-8., 8., 401]}];
mollowA6[Ω_, Δ_] := With[{L = N @ Normal[QuantumOperator["Liouvillian"[Δ/2 QuantumOperator["Z"] + Ω/2 QuantumOperator["X"], {QuantumOperator["-"]}, {1}]]["Matrix"]],
     sm = N @ Normal[QuantumOperator["-"]["Matrix"]]},
    With[{ρss = #/Tr[#] &[ArrayReshape[First[NullSpace[L]], {2, 2}]], c = -I sm},
     {With[{f = specA6[L, (sm - Tr[sm . ρss] IdentityMatrix[2]) . ρss, ConjugateTranspose[sm]]}, MapAt[#/Max[f[[All, 2]]] &, f, {All, 2}]],
      MapAt[1 + # &, specA6[L, c . ρss + ρss . ConjugateTranspose[c] - Tr[(c + ConjugateTranspose[c]) . ρss] ρss, c + ConjugateTranspose[c]], {All, 2}]}]];
{ΩsA6, ΔsA6} = {{0.3, 1, 5}, {0, 1, 3}};
dataA6 = Outer[mollowA6, ΩsA6, ΔsA6];
panelA6[row_, i_] := ListLinePlot[dataA6[[i, All, row]], PlotRange -> All, PlotLegends -> (Row[{"Δ/γ = ", #}] & /@ ΔsA6),
    Frame -> True, Axes -> False, GridLines -> {None, If[row == 2, {1}, None]}, FrameLabel -> {"ν/γ", {"S(ν)/max", Row[{Subscript["S", "x"], "(ν)"}]}[[row]]},
    PlotLabel -> Row[{"Ω/γ = ", ΩsA6[[i]]}], ImageSize -> 300];
panelsA6 = Outer[panelA6, {1, 2}, {1, 2, 3}];
Legended[Grid[Map[First, panelsA6, {2}]], Last[panelsA6[[1, 1]]]]
```

At weak resonant driving the incoherent part of the spectrum is a single line at the laser frequency; at $\Omega/\gamma=5$ it splits into the Mollow triplet, whose sidebands detuning pushes further out. With detuning, even weak driving moves the incoherent light into sidebands near $\pm\Delta$. At weak resonant driving the out-of-phase quadrature fluctuates below shot noise around $\nu=0$: a single atom emits squeezed light. Stronger driving or detuning lifts this quadrature back above shot noise near $\nu=0$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A qutrit Wigner function under noise: the negativity dies at an exact time

A qutrit state can be written as a discrete Wigner function $W$: nine real numbers on a $3\times3$ phase-space grid that sum to one, whose negative entries make it a quasi-probability rather than a probability. A qutrit with degenerate levels, so no Hamiltonian, is subject only to depolarizing noise at rate $\gamma$, and its negative entries fill in. When, exactly, does the last one reach zero?

With $X|j\rangle=|j+1\rangle$ and $Z|j\rangle=e^{-2\pi ij/3}|j\rangle$ (indices mod 3), the noise is the Lindblad equation with the eight Weyl operators $X^aZ^b$, $(a,b)\neq(0,0)$, as jump operators at rate $\gamma/9$ each,

$$
\dot\rho=\frac{\gamma}{9}\sum_{(a,b)\neq(0,0)}\mathcal D\!\left[X^aZ^b\right]\rho .
$$

QF solves it with $\gamma$ and $t$ symbolic from three initial states, reads $W(t)$ off the solution, and solves for the time $t^*$ at which the last negative entry reaches zero (dashed); each curve is a set of phase-space points with equal $W$, and time is in units of $1/\gamma$:

```wl
{x3A7, z3A7} = Normal[QuantumOperator[#]["Matrix"]] & /@ {"X"[3], "Z"[3]};
weylA7 = QuantumOperator[MatrixPower[x3A7, #1] . MatrixPower[z3A7, #2]] & @@@ Rest[Tuples[Range[0, 2], 2]];
wignerA7[ψ_] := DeleteDuplicates[Flatten[Simplify[Normal[QuantumEvolve[QuantumOperator[0 IdentityMatrix[3]], weylA7 -> ConstantArray[γ/9, 8],
       QuantumState[ψ], t]["PhaseSpace"]] /. γ -> 1, t >= 0]]];
deathA7[w_] := Max[t /. Solve[# == 0 && t > 0, t, Reals] & /@ Select[w, (# /. t -> 0) < 0 &]];
statesA7 = <|"(|1\[RightAngleBracket] - |2\[RightAngleBracket])/\[Sqrt]2" -> {0, 1, -1}/Sqrt[2],
    "(|0\[RightAngleBracket] + i|1\[RightAngleBracket])/\[Sqrt]2" -> {1, I, 0}/Sqrt[2], "(|0\[RightAngleBracket] + |1\[RightAngleBracket])/\[Sqrt]2" -> {1, 1, 0}/Sqrt[2]|>;
panelA7[label_, ψ_] := With[{w = wignerA7[ψ]}, With[{tstar = deathA7[w]},
    Plot[Evaluate[w], {t, 0, 3}, PlotRange -> All, Frame -> True, GridLines -> {{tstar}, {0}}, GridLinesStyle -> Dashed,
     FrameLabel -> {"γ t", "W"}, PlotLabel -> Column[{label, Row[{"γ t* = ", tstar}]}, Center], ImageSize -> 280]]];
Row[KeyValueMap[panelA7, statesA7]]
```

Every entry relaxes toward the uniform value $1/9$ of the maximally mixed state, so each negative entry crosses zero at one exact time, later the more negative it starts. Of the three states, $(|1\rangle-|2\rangle)/\sqrt2$ starts most negative and holds out longest. After $t^*$ the Wigner function is a probability distribution.
