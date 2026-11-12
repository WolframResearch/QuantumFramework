```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A qubit driven by a rotating field

A qubit with splitting $\omega_0$ is driven by a field of strength $\Omega$ rotating at frequency $\omega$, detuned by $\delta=\omega-\omega_0$. How much population does the drive move out of $|0\rangle$, and how fast?

With $\hbar=1$, the lab-frame Hamiltonian is

$$
H(t)=\frac{\omega_0}{2}\,\sigma_z+\frac{\Omega}{2}\left(\cos\omega t\,\sigma_x+\sin\omega t\,\sigma_y\right),
$$

which does not commute with itself at different times; the propagator solves $i\,\dot U=H(t)\,U$ with $U(0)=\mathbb 1$.

QF returns the propagator with every parameter a symbol; the cell compares $P(1)=|\langle1|U(t)|0\rangle|^2$ with the Rabi formula and plots it, with time in units of $1/\Omega$:

```wl
rabiP = FullSimplify[ComplexExpand[Abs[Normal[QuantumEvolve[\[Omega]0/2 QuantumOperator["Z"] + \[CapitalOmega]/2 (Cos[\[Omega] t] QuantumOperator["X"] + Sin[\[Omega] t] QuantumOperator["Y"]), None, t]["Matrix"]][[2, 1]]]^2],
    {\[Omega]0, \[Omega], t} \[Element] Reals && \[CapitalOmega] > 0];
rabiQ = FullSimplify[(rabiP /. \[Omega] -> \[Omega]0 + \[Delta]) == \[CapitalOmega]^2/(\[CapitalOmega]^2 + \[Delta]^2) Sin[Sqrt[\[CapitalOmega]^2 + \[Delta]^2] t/2]^2, {\[Omega]0, \[Delta], t} \[Element] Reals && \[CapitalOmega] > 0];
rabiDeltas = {0, 1, 2};
Plot[Evaluate[rabiP /. {\[CapitalOmega] -> 1, \[Omega] -> \[Omega]0 + #} & /@ rabiDeltas], {t, 0, 4 Pi}, PlotRange -> {0, 1}, Frame -> True,
  FrameLabel -> {"\[CapitalOmega] t", "P(1)"}, ImageSize -> 420, PlotLegends -> (Row[{"\[Delta]/\[CapitalOmega] = ", #}] & /@ rabiDeltas),
  PlotLabel -> Row[{"P(1) = \!\(\*FractionBox[SuperscriptBox[\(\[CapitalOmega]\), \(2\)], \(\*SuperscriptBox[\(\[CapitalOmega]\), \(2\)] + \*SuperscriptBox[\(\[Delta]\), \(2\)]\)]\) \!\(\*SuperscriptBox[\(sin\), \(2\)]\)(\!\(\*SqrtBox[\(\*SuperscriptBox[\(\[CapitalOmega]\), \(2\)] + \*SuperscriptBox[\(\[Delta]\), \(2\)]\)]\) t/2): ", rabiQ}]]
```

The drive moves at most $\Omega^2/(\Omega^2+\delta^2)$ of the population, and it oscillates at the generalized Rabi frequency $\sqrt{\Omega^2+\delta^2}$, faster than on resonance. Only the detuning enters: $\omega_0$ and $\omega$ separately drop out.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A Landau-Zener sweep

A qubit has its bias swept at rate $v$ through an avoided crossing with gap $\Delta$. It starts in $|0\rangle$ long before the crossing. How does it reach its final population?

With $\hbar=1$, the Hamiltonian is

$$
H(t)=\frac{vt}{2}\,\sigma_z+\frac{\Delta}{2}\,\sigma_x,
$$

which does not commute with itself at different times; the propagator solves $i\,\dot U=H(t)\,U$ with $U(0)=\mathbb 1$, and for an infinitely long sweep the probability to stay in $|0\rangle$ is the Landau-Zener value $e^{-\pi\Delta^2/(2v)}$.

QF returns the propagator with $v$, $\Delta$ and $t$ symbolic. The cell starts the qubit in $|0\rangle$ at $\sqrt v\,t=-20$ and plots $P(0)$ from $\sqrt v\,t=-10$, with the Landau-Zener value dashed:

```wl
lzU = Normal[QuantumEvolve[v t/2 QuantumOperator["Z"] + \[CapitalDelta]/2 QuantumOperator["X"], None, t]["Matrix"]];
lzStay[r_] := With[{u = lzU /. {v -> 1, \[CapitalDelta] -> Sqrt[r]}}, Abs[(u . ConjugateTranspose[u /. t -> -20])[[1, 1]]]^2];
lzRatios = {0.1, 0.6, 1.6};
Show[Plot[Evaluate[lzStay /@ lzRatios], {t, -10, 10}, PlotRange -> All, PlotLegends -> (Row[{"\[CapitalDelta]\.b2/v = ", #}] & /@ lzRatios)],
  Plot[Evaluate[Exp[-Pi lzRatios/2]], {t, -10, 10}, PlotStyle -> Dashed],
  PlotRange -> {0, 1}, Frame -> True, Axes -> False, FrameLabel -> {Row[{Sqrt["v"], " t"}], "P(0)"}, ImageSize -> 420]
```

The qubit drops through the crossing at $t=0$ and then settles toward the Landau-Zener value through oscillations that fade as the levels separate; the asymptotic formula does not contain them. The small ripples before the crossing come from starting the sweep at a finite time. A slower sweep, larger $\Delta^2/v$, leaves less population in $|0\rangle$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A Gaussian pulse on a qubit

A qubit with splitting $\omega_0$ is driven at frequency $\omega$, detuned by $\Delta=\omega_0-\omega$, by a Gaussian pulse of peak Rabi frequency $\Omega_0$, width $\tau$, centered at $t_0$. Which gate has the pulse made by time $t$?

With $\hbar=1$, the lab-frame Hamiltonian and its form in the frame rotating at $\omega$, after the rotating-wave approximation, are

$$
H=\frac{\omega_0}{2}\,\sigma_z+\Omega(t)\cos\omega t\,\sigma_x\;\longrightarrow\;\frac{\Delta}{2}\,\sigma_z+\frac{\Omega(t)}{2}\,\sigma_x,\qquad \Omega(t)=\Omega_0\,e^{-(t-t_0)^2/2\tau^2},
$$

and the propagator solves $i\,\dot U=H\,U$ with $U(0)=\mathbb 1$.

`DSolveValue` solves it on resonance for the whole $2\times2$ matrix, with $\Omega_0$, $\tau$, $t_0$, $t$ symbolic, and compares it with the rotation by the pulse area so far, $e^{-i\theta(t)\sigma_x/2}$ with $\theta(t)=\int_0^t\Omega(s)\,ds$; off resonance `NDSolveValue` solves the same equation. The cell plots $P(1)$ from $|0\rangle$ for pulses centered at $t_0=4\tau$ with total area $A=\sqrt{2\pi}\,\Omega_0\tau$, with $\sin^2(\theta(t)/2)$ dashed and time in units of $\tau$:

```wl
pulseEq[d_, w_] := {upul'[t] == -I (d/2 PauliMatrix[3] + w/2 Exp[-(t - t0)^2/(2 \[Tau]^2)] PauliMatrix[1]) . upul[t], upul[0] == IdentityMatrix[2]};
upulse = DSolveValue[pulseEq[0, \[CapitalOmega]0], upul[t] \[Element] Matrices[{2, 2}], t];
pulseArea = Integrate[\[CapitalOmega]0 Exp[-(s - t0)^2/(2 \[Tau]^2)], {s, 0, t}, Assumptions -> \[Tau] > 0];
rotationQ = FullSimplify[upulse == MatrixExp[-I pulseArea/2 PauliMatrix[1]], \[Tau] > 0];
{pulseDetunings, pulseAreas} = {{0, 0.5, 2}, {Pi/2, Pi, 2 Pi}};
pulseP1[d_, a_] := With[{num = {\[CapitalOmega]0 -> a/Sqrt[2 Pi], \[Tau] -> 1, t0 -> 4}},
    If[d == 0, Abs[upulse[[2, 1]]]^2 /. num, NDSolveValue[pulseEq[d, \[CapitalOmega]0] /. num, Abs[Indexed[upul[t], {2, 1}]]^2, {t, 0, 8}]]];
pulsePanel[d_] := Show[
    Plot[Evaluate[pulseP1[d, #] & /@ pulseAreas], {t, 0, 8}, PlotRange -> All, PlotLegends -> (Row[{"A = ", #}] & /@ pulseAreas)],
    Plot[Evaluate[Sin[pulseArea/2]^2 /. {\[CapitalOmega]0 -> #/Sqrt[2 Pi], \[Tau] -> 1, t0 -> 4} & /@ pulseAreas], {t, 0, 8}, PlotStyle -> Dashed],
    PlotRange -> {0, 1}, Frame -> True, FrameLabel -> {"t/\[Tau]", "P(1)"}, ImageSize -> 300,
    PlotLabel -> If[d == 0, Row[{"\[CapitalDelta]\[Tau] = 0,  U = \!\(\*SuperscriptBox[\(e\), \(-i\[Theta](t)\*SubscriptBox[\(\[Sigma]\), \(x\)]/2\)]\): ", rotationQ}], Row[{"\[CapitalDelta]\[Tau] = ", d}]]];
panelsP = pulsePanel /@ pulseDetunings;
Legended[Grid[{First /@ panelsP}], Last[First[panelsP]]]
```

On resonance the solved propagator is exactly the rotation by the area swept so far, so the population follows $\sin^2(\theta(t)/2)$ and a pulse with $A=\pi$ is a full flip; this holds because the Hamiltonian then points along $\sigma_x$ at every instant and commutes with itself at different times. Detuning tilts the rotation axis out of the equator by an angle that changes along the pulse, so the population no longer follows the area, and at $\Delta\tau=2$ the pulses leave almost nothing in $|1\rangle$ at the end.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Two spins coupled by exchange

Two spins interact through an XXZ exchange with in-plane coupling $J$ and axial coupling $J_z$. The second spin starts in $|1\rangle$ and the first is tilted from $|0\rangle$ by a polar angle $\vartheta$. How entangled do they become?

With $\hbar=1$ and Pauli operators $X,Y,Z$ on each spin,

$$
H=J\,(XX+YY)+J_z\,ZZ,\qquad |\psi(0)\rangle=\left(\cos\tfrac{\vartheta}{2}|0\rangle+\sin\tfrac{\vartheta}{2}|1\rangle\right)\otimes|1\rangle .
$$

QF evolves the state with $J$, $J_z$, $\vartheta$ and $t$ symbolic and traces out the second spin; the cell plots the von Neumann entropy $S=-\mathrm{Tr}\,\rho_1\log_2\rho_1$ of the first spin, with time in units of $1/J$, and checks the reduced state from $|01\rangle$:

```wl
rhoXXZ = Normal[QuantumPartialTrace[QuantumEvolve[J (QuantumOperator["XX"] + QuantumOperator["YY"]) + Jz QuantumOperator["ZZ"],
    QuantumTensorProduct[QuantumState[{Cos[\[CurlyTheta]/2], Sin[\[CurlyTheta]/2]}], QuantumState["1"]], t], {2}]["DensityMatrix"]];
xxzEntropy[a_, jz_][s_?NumericQ] := QuantumState[N[rhoXXZ /. {J -> 1, Jz -> jz, \[CurlyTheta] -> a, t -> s}]]["VonNeumannEntropy", 2];
{tilts, axials} = {{0, Pi/3, 2 Pi/3}, {0, 2}};
xxzPanel[jz_] := Plot[Evaluate[xxzEntropy[#, jz][t] & /@ tilts], {t, 0, Pi/2}, PlotRange -> {0, 1.02}, Frame -> True,
    FrameLabel -> {"J t", "S (bits)"}, PlotLabel -> Row[{"\!\(\*SubscriptBox[\(J\), \(z\)]\)/J = ", jz}], ImageSize -> 330,
    PlotLegends -> (Row[{"\[CurlyTheta] = ", #}] & /@ tilts)];
xxzPanels = xxzPanel /@ axials;
Labeled[Legended[Grid[{First /@ xxzPanels}], Last[First[xxzPanels]]],
  Row[{"\!\(\*SubscriptBox[\(\[Rho]\), \(1\)]\) from |01\[RightAngleBracket] equals diag(cos\.b2 2Jt, sin\.b2 2Jt): ",
    FullSimplify[(rhoXXZ /. \[CurlyTheta] -> 0) == DiagonalMatrix[{Cos[2 J t]^2, Sin[2 J t]^2}], {J, Jz, t} \[Element] Reals]}], Top]
```

From $|01\rangle$ the reduced state is diagonal, so the excitation swaps between the spins at angular frequency $4J$ and the entropy is exactly the binary entropy $h_2(\cos^2 2Jt)$: one bit whenever the excitation is shared equally, zero whenever it sits on one spin. Tilting the first spin puts weight on $|11\rangle$, which the exchange cannot move, so the entanglement stays below one bit. The two panels are identical: spin 2 starts in $|1\rangle$ and the exchange conserves the number of excitations, so the state never acquires $|00\rangle$; on the other three states $J_zZZ$ equals $-J_z(1+Z_1+Z_2)$, a rotation of each spin about $z$ that commutes with $H$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A coherent state in a harmonic oscillator

A field mode of frequency $\omega$ starts in the coherent state $|\alpha\rangle$. Does it stay coherent?

With $\hbar=1$, $H=\omega\left(a^\dagger a+\tfrac12\right)$, and the expected answer is

$$
e^{-iHt}|\alpha\rangle=e^{-i\omega t/2}\,|\alpha e^{-i\omega t}\rangle .
$$

QF computes $e^{-iHt}$ with $\omega$ and $t$ symbolic in a Fock space of 16 levels, applies it to the unnormalized $|\alpha\rangle$ with $\alpha$ symbolic, and compares both sides amplitude by amplitude. The cell plots the Wigner function $W(x,p)$, with $x+ip=\sqrt2\,\alpha$, of the evolved state at $\alpha=2$:

```wl
uOsc = QuantumEvolve[\[Omega] (SuperDagger[AnnihilationOperator[]] @ AnnihilationOperator[] + 1/2), None, t];
evolvedAmps = (uOsc @ CoherentState["Normalized" -> False])["AmplitudesList"];
oscSameQ = Simplify[evolvedAmps == E^(-I \[Omega] t/2) CoherentState["Normalized" -> False][\[FormalAlpha] E^(-I \[Omega] t)]["AmplitudesList"]];
oscTimes = {0, Pi/3, 2 Pi/3, Pi};
oscPanel[wt_] := With[{w = WignerFunction[QuantumState[Normalize[evolvedAmps /. {\[FormalAlpha] -> 2, \[Omega] -> 1, t -> wt}], 16], {x, p}]},
    DensityPlot[w, {x, -4.5, 4.5}, {p, -4.5, 4.5}, PlotPoints -> 40, PlotRange -> All, FrameLabel -> {"x", "p"},
      PlotLabel -> "\[Omega] t = " <> ToString[wt, InputForm], ImageSize -> 220]];
Labeled[GraphicsRow[oscPanel /@ oscTimes, ImageSize -> 900],
  Row[{"\!\(\*SuperscriptBox[\(e\), \(-iHt\)]\)|\[Alpha]\[RightAngleBracket] = \!\(\*SuperscriptBox[\(e\), \(-i\[Omega]t/2\)]\)|\[Alpha]\!\(\*SuperscriptBox[\(e\), \(-i\[Omega]t\)]\)\[RightAngleBracket] for all 16 amplitudes: ", oscSameQ}], Top]
```

The identity holds for every $\alpha$, $\omega$ and $t$, phase $e^{-i\omega t/2}$ included: the propagator multiplies the $n$-th amplitude by $e^{-i\omega(n+1/2)t}$, so the cutoff never mixes levels. In phase space the round Gaussian of the vacuum moves clockwise on a circle of radius $\sqrt2\,|\alpha|$ without changing shape.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Ehrenfest and Mandelstam-Tamm for any drive

A qubit is driven by an arbitrary time-dependent Hamiltonian. Do Ehrenfest's theorem and the quantum speed limit hold without solving for the state?

With $\hbar=1$, four undetermined real drive functions and an arbitrary constant observable $A$,

$$
H(t)=h_0(t)\,\mathbb 1+h_x(t)\,\sigma_x+h_y(t)\,\sigma_y+h_z(t)\,\sigma_z,\qquad \tfrac{d}{dt}|\psi\rangle=-iH|\psi\rangle,\qquad \tfrac{d}{dt}\langle\psi|=i\langle\psi|H .
$$

Ehrenfest's theorem says $\tfrac{d}{dt}\langle A\rangle=i\langle[H,A]\rangle$; the squared Fubini-Study speed $\langle\dot\psi|\dot\psi\rangle-|\langle\psi|\dot\psi\rangle|^2$ should equal the energy variance $\Delta E^2$, so the distance $\arccos|\langle\psi_0|\psi_t\rangle|$ can never exceed the path length $\int_0^t\Delta E\,dt'$, the Mandelstam-Tamm bound. The cell computes both residuals with the Schrödinger equation used only as a rule for time derivatives, then integrates three random drives from $|0\rangle$, with drive amplitudes and frequencies of order $\epsilon$ and time in units of $1/\epsilon$; the path length is dashed:

```wl
ehH = h0[t] IdentityMatrix[2] + hx[t] PauliMatrix[1] + hy[t] PauliMatrix[2] + hz[t] PauliMatrix[3];
{ehKet, ehBra, ehA} = {{ehC1[t], ehC2[t]}, {ehB1[t], ehB2[t]}, Array[\[FormalCapitalA], {2, 2}]};
ehRules = Join[Thread[D[ehKet, t] -> -I ehH . ehKet], Thread[D[ehBra, t] -> I ehBra . ehH]];
ehrenfestResidual = Simplify[(D[ehBra . ehA . ehKet, t] /. ehRules) - I ehBra . (ehH . ehA - ehA . ehH) . ehKet];
speedResidual = Simplify[((D[ehBra, t] . D[ehKet, t] - (ehBra . D[ehKet, t]) (D[ehBra, t] . ehKet)) /. ehRules) - (ehBra . ehH . ehH . ehKet - (ehBra . ehH . ehKet)^2)];
mtRun[seed_] := With[{hn = Total[MapThread[(#1[[1]] + #1[[2]] Cos[3 #1[[3]] s + 2 Pi #1[[4]]]) #2 &,
      {BlockRandom[SeedRandom[seed]; RandomReal[{-1, 1}, {4, 4}]], Prepend[PauliMatrix /@ {1, 2, 3}, IdentityMatrix[2]]}]]},
    NDSolveValue[{{ehC1'[s], ehC2'[s]} == -I hn . {ehC1[s], ehC2[s]}, {ehC1[0], ehC2[0]} == {1, 0},
      ehLen'[s] == Sqrt[Abs[Conjugate[{ehC1[s], ehC2[s]}] . hn . hn . {ehC1[s], ehC2[s]} - (Conjugate[{ehC1[s], ehC2[s]}] . hn . {ehC1[s], ehC2[s]})^2]], ehLen[0] == 0},
      {ArcCos[Min[1, Abs[ehC1[s]]]], ehLen[s]}, {s, 0, 6}]];
mtRuns = mtRun /@ {1, 2, 3};
Show[Plot[Evaluate[First /@ mtRuns], {s, 0, 6}, PlotRange -> All, PlotLegends -> {"drive 1", "drive 2", "drive 3"}],
  Plot[Evaluate[Last /@ mtRuns], {s, 0, 6}, PlotStyle -> Dashed, PlotRange -> All],
  PlotRange -> {0, All}, Frame -> True, FrameLabel -> {"\[Epsilon] t", "arccos|\[LeftAngleBracket]\!\(\*SubscriptBox[\(\[Psi]\), \(0\)]\)|\!\(\*SubscriptBox[\(\[Psi]\), \(t\)]\)\[RightAngleBracket]|"}, ImageSize -> 440,
  PlotLabel -> Row[{"Ehrenfest residual: ", ehrenfestResidual, ",   speed\.b2 \[Minus] \[CapitalDelta]\!\(\*SuperscriptBox[\(E\), \(2\)]\): ", speedResidual}]]
```

Both residuals vanish for every drive, every observable and every state: the two identities follow from the equation of motion alone, so they hold even for drives whose propagator has no closed form. Along each random drive the distance meets its path length only at the start and then stays below it: the state can drift back toward $|0\rangle$, while the path length only grows.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A spin coherent state

A spin $j$ starts in $|j,j\rangle$, pointing along $+z$, and is rotated to point along the direction $\hat n$ with polar angle $\theta$ and azimuth $\varphi$. How is the rotated spin spread over the levels $m$ of $J_z$?

With $J_z|j,m\rangle=m|j,m\rangle$, the spin coherent state and its level populations are

$$
|j,\hat n\rangle=e^{-i\varphi J_z}e^{-i\theta J_y}|j,j\rangle,\qquad P(m)=|\langle j,m|j,\hat n\rangle|^2=\left|d^{\,j}_{mj}(\theta)\right|^2,
$$

where $d^{\,j}_{m'm}(\theta)=\langle j,m'|e^{-i\theta J_y}|j,m\rangle$ is the Wigner small-$d$ matrix; $\varphi$ only multiplies each amplitude by a phase.

QF builds the rotation as the named operator `"WignerD"` with $\theta$ symbolic, applies it to $|j,j\rangle$, and checks $P(m)$ against the binomial law; the cell plots $P(m)$ against $m/j$, with $\cos\theta$ dashed:

```wl
spinCohP[j_] := ComplexExpand[Abs[Normal[(QuantumOperator["WignerD"[j, \[Theta]]] @ QuantumState[UnitVector[2 j + 1, 2 j + 1], 2 j + 1])["StateVector"]]]^2];
{spins, polars} = {{2, 10, 40}, {Pi/6, Pi/2, 5 Pi/6}};
spinCoh = AssociationMap[spinCohP, spins];
binomialQ = And @@ (Simplify[spinCoh[#] == Table[Binomial[2 #, # + m] Cos[\[Theta]/2]^(2 (# + m)) Sin[\[Theta]/2]^(2 (# - m)), {m, -#, #}], 0 < \[Theta] < Pi] & /@ spins);
spinPanel[j_] := ListLinePlot[Table[Transpose[{Range[-j, j]/j, spinCoh[j] /. \[Theta] -> th}], {th, polars}], PlotMarkers -> Automatic, PlotRange -> All, Axes -> False,
    Frame -> True, FrameLabel -> {"m/j", "P(m)"}, PlotLabel -> Row[{"j = ", j}], GridLines -> {Cos[polars], None}, GridLinesStyle -> Dashed,
    ImageSize -> 300, PlotLegends -> (Row[{"\[Theta] = ", #}] & /@ polars)];
spinPanels = spinPanel /@ spins;
Labeled[Legended[Grid[{First /@ spinPanels}], Last[First[spinPanels]]],
  Row[{"P(m) = \!\(\*SubscriptBox[\(C\), \(2j, j+m\)]\) \!\(\*SuperscriptBox[\(cos\), \(2(j+m)\)]\)(\[Theta]/2) \!\(\*SuperscriptBox[\(sin\), \(2(j-m)\)]\)(\[Theta]/2) for 0 < \[Theta] < \[Pi] at j = 2, 10, 40: ", binomialQ}], Top]
```

The populations are binomial: the spin-$j$ coherent state is $2j$ spin-$\tfrac12$ particles all pointing along $\hat n$, each found up with probability $\cos^2(\theta/2)$. The distribution is centered at $m=j\cos\theta$, the projection of a classical vector of length $j$, and narrows relative to $j$ as $j$ grows, so a larger spin points ever more sharply, like a classical arrow.
