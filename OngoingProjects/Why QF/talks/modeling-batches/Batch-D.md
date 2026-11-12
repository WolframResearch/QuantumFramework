```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A trapped particle watched continuously: squeezing below the zero-point

A particle of mass $m$ in a harmonic trap of frequency $\omega$ has its position measured continuously at rate $\lambda$ by a detector that catches a fraction $\eta$ of the signal. Can watching make it narrower in position than its ground state?

With $H=\frac{p^2}{2m}+\frac{m\omega^2x^2}{2}$ and $\mathcal D[x]\rho=x\rho x-\frac12\{x^2,\rho\}$, the state conditioned on the record obeys

$$
d\rho=-\frac{i}{\hbar}[H,\rho]\,dt+\lambda\,\mathcal D[x]\rho\,dt+\sqrt{\eta\lambda}\left(x\rho+\rho x-2\langle x\rangle\rho\right)dW .
$$

A Gaussian state stays Gaussian, and its covariances $\Sigma_x$, $C$, $\Sigma_p$ obey a Riccati equation with no noise in it:

$$
\dot\Sigma_x=\frac{2C}{m}-4\eta\lambda\Sigma_x^2,\qquad \dot C=\frac{\Sigma_p}{m}-m\omega^2\Sigma_x-4\eta\lambda\Sigma_xC,\qquad \dot\Sigma_p=-2m\omega^2C+\lambda\hbar^2-4\eta\lambda C^2 .
$$

WL finds the steady state symbolically, keeps the root that `Resolve` proves positive, and reduces the squeezing condition (panel titles); the curves integrate the flow from the ground state, dashed lines are the closed form, $\Sigma_x$ is in units of the zero-point variance $\hbar/(2m\omega)$, the strength is $\mu=2\hbar\lambda/(m\omega^2)$ and time is in units of $1/\omega$:

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
Legended[Row[First /@ trapPanels], Last[First[trapPanels]]]
```

A perfect detector squeezes the position below the zero-point at every strength. An imperfect one squeezes only above a threshold strength that grows as $\eta$ falls; below it, watching leaves the particle wider than its ground state.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The quantum Kalman filter: five numbers against a wavefunction

A particle of mass $m$ sits in the potential $V(x)=\frac12m\omega^2x^2+\epsilon x^4$, and a perfect detector reads its position continuously at rate $\lambda$. Can five numbers, updated from the readings alone, follow the mean position that the full conditional wavefunction carries?

In a step $dt$ the detector returns a reading $\bar x$, distributed as $|\psi(x)|^2$ blurred by a Gaussian of variance $\frac{1}{4\lambda\,dt}$. Bayes' rule updates the wavefunction, and the Hamiltonian $H=p^2/2m+V(x)$ moves it ($\hbar=1$):

$$
\psi\;\to\;e^{-iH\,dt}\,\frac{e^{-\lambda dt\,(x-\bar x)^2}\psi}{\lVert e^{-\lambda dt\,(x-\bar x)^2}\psi\rVert}.
$$

The Kalman filter assumes a Gaussian state in the harmonic trap, $\epsilon=0$, and carries only its means and covariances, corrected by the innovation $\bar x-\langle x\rangle$, the reading minus the predicted mean:

$$
d\langle x\rangle=\frac{\langle p\rangle}{m}dt+4\lambda\Sigma_x(\bar x-\langle x\rangle)dt,\qquad d\langle p\rangle=-m\omega^2\langle x\rangle dt+4\lambda C(\bar x-\langle x\rangle)dt,
$$
$$
\dot\Sigma_x=\frac{2C}{m}-4\lambda\Sigma_x^2,\qquad \dot C=\frac{\Sigma_p}{m}-m\omega^2\Sigma_x-4\lambda\Sigma_xC,\qquad \dot\Sigma_p=-2m\omega^2C+\lambda-4\lambda C^2 .
$$

WL runs the wavefunction on 256 grid points with a split-operator step, from a coherent state displaced to $x=3$, and feeds the same readings to the filter. Lengths are in units of $\sqrt{\hbar/(m\omega)}$, time in $1/\omega$, the strength is $\mu=2\hbar\lambda/(m\omega^2)$ and $\epsilon$ is in units of $m^2\omega^3/\hbar$:

```wl
kalmanX = 24. Range[-128, 127]/256; kalmanP = 2 Pi RotateLeft[Range[-128, 127], 128]/24.;
kalmanGrid[μ_, ε_] := With[{λ0 = μ/2, dt = 0.002}, With[{half = Exp[-I (kalmanX^2/2 + ε kalmanX^4) dt/2], kin = Exp[-I kalmanP^2 dt/2]},
    Reap[NestList[With[{xb = RandomChoice[Abs[#]^2 -> kalmanX] + RandomVariate[NormalDistribution[0, 1/Sqrt[4 λ0 dt]]]},
        half InverseFourier[kin Fourier[half Normalize[Exp[-λ0 dt (kalmanX - Sow[xb])^2] #]]]] &, Normalize[Exp[-(kalmanX - 3)^2/2]], 5000]]]];
kalmanFilter[μ_, readings_] := With[{λ0 = μ/2, dt = 0.002}, FoldList[Function[{v, xb}, With[{k = 4 λ0 (xb - v[[1]]) dt},
      v + {v[[2]] dt + k v[[3]], -v[[1]] dt + k v[[4]], (2 v[[4]] - 4 λ0 v[[3]]^2) dt, (v[[5]] - v[[3]] - 4 λ0 v[[3]] v[[4]]) dt, (-2 v[[4]] + λ0 - 4 λ0 v[[4]]^2) dt}]],
    {3., 0., 0.5, 0., 0.5}, readings][[All, 1]]];
kalmanPanel[μ_, ε_] := With[{run = BlockRandom[SeedRandom[7]; kalmanGrid[μ, ε]]},
    ListLinePlot[{Abs[#]^2 . kalmanX & /@ First[run], kalmanFilter[μ, run[[2, 1]]]}, DataRange -> {0, 10}, PlotStyle -> {Automatic, Dashed},
     PlotRange -> All, Frame -> True, FrameLabel -> {"ω t", "⟨x⟩"}, PlotLabel -> Row[{"μ = ", μ, ",   ε = ", ε}], ImageSize -> 300,
     PlotLegends -> {"wavefunction, 256 amplitudes", "Kalman filter, 5 numbers"}]];
kalmanPanels = Outer[kalmanPanel, {1, 10}, {0, 0.01, 0.05}];
Legended[Grid[Map[First, kalmanPanels, {2}]], Last[kalmanPanels[[1, 1]]]]
```

In the harmonic trap the filter lies on the wavefunction's mean at both strengths: the conditional state stays Gaussian, so the 256 amplitudes carry nothing the five numbers miss. A quartic term makes the state non-Gaussian and the filter's harmonic model wrong. It drifts away under gentle watching. Under strong watching each reading pulls it back, and it falls visibly behind the faster true oscillation only at the larger $\epsilon$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A watched Gaussian: the center wanders, the widths do not

A particle of mass $m$ in a harmonic trap of frequency $\omega$ is watched in position at rate $\lambda$ by a detector of efficiency $\eta$, and its conditional state is a Gaussian with means $\langle x\rangle,\langle p\rangle$, variances $\Sigma_x,\Sigma_p$ and covariance $C$. Every run moves the packet's center at random. Do the widths of a watched Gaussian wander too?

With $H=\frac{p^2}{2m}+\frac{m\omega^2x^2}{2}$ and $[x,p]=i\hbar$, the conditional mean of any observable $O$ obeys

$$
d\langle O\rangle=\Big\langle \tfrac{i}{\hbar}[H,O]-\tfrac{\lambda}{2}\big[x,[x,O]\big]\Big\rangle dt+\sqrt{\eta\lambda}\,\big\langle\{x-\langle x\rangle,O\}\big\rangle\,dW,
$$

and the second moments pick up the Itô correction, e.g. $d\Sigma_x=d\langle x^2\rangle-2\langle x\rangle d\langle x\rangle-(d\langle x\rangle)^2$. The Gaussian has the ordered generating function $\langle e^{\alpha x}e^{\beta p}\rangle=\exp\!\big(\alpha\langle x\rangle+\beta\langle p\rangle+\tfrac12(\alpha^2\Sigma_x+2\alpha\beta C+\beta^2\Sigma_p)+\tfrac{i\hbar}{2}\alpha\beta\big)$, and its purity is $\hbar/(2\sqrt{\Sigma_x\Sigma_p-C^2})$.

WL hands the commutator to `NonCommutativeAlgebra`, which orders every product with $x$ to the left of $p$, reads each ordered monomial off the generating function, and returns the drift and noise of the five numbers with every parameter symbolic (top). The two columns then define an `ItoProcess`, sampled four times per panel from $\langle x\rangle=2$ and three times the zero-point variances, with $\hbar=m=1$, $\mu=2\hbar\lambda/(m\omega^2)=4$, time in units of $1/\omega$ and $\sqrt\eta$ dashed:

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
gaussRuns[η0_] := BlockRandom[SeedRandom[4]; RandomFunction[ItoProcess[{gaussEqs[[All, 1]], List /@ gaussEqs[[All, 2]]} /. Thread[{mx, mp, Σx, Cxp, Σp} -> Through[{mx, mp, Σx, Cxp, Σp}[t]]] /.
        {m -> 1, ω -> 1, ℏ -> 1, λ -> 2, η -> η0}, {{mx, mp, Σx, Cxp, Σp}, {2, 0, 3/2, 0, 3/2}}, t], {0., 6., 0.001}, 4]["Paths"]];
gaussPanel[η0_, runs_, {label_, f_, range_, grid_}] := ListLinePlot[Map[{#[[1]], f[#[[2]]]} &, runs, {2}], PlotRange -> range, Frame -> True, Axes -> False,
    GridLines -> {None, grid}, GridLinesStyle -> Dashed, FrameLabel -> {"ω t", label}, PlotLabel -> Row[{"η = ", η0, ",   μ = 4"}], ImageSize -> 300];
gaussShow = {mx -> AngleBracket["x"], mp -> AngleBracket["p"], Σx -> Subscript["Σ", "x"], Σp -> Subscript["Σ", "p"], Cxp -> "C"};
Labeled[Grid[Transpose[Function[η0, With[{runs = gaussRuns[η0]}, gaussPanel[η0, runs, #] & /@
       {{"⟨x⟩", First, All, None}, {"purity", 1/(2 Sqrt[#[[3]] #[[5]] - #[[4]]^2]) &, {0, 1.05}, {Sqrt[η0]}}}]] /@ {1, 1/5}]],
  Style[Grid[Prepend[MapThread[{TraditionalForm[#1], TraditionalForm[#2 /. gaussShow], TraditionalForm[#3 /. gaussShow]} &,
      {{"⟨x⟩", "⟨p⟩", Subscript["Σ", "x"], "C", Subscript["Σ", "p"]}, gaussEqs[[All, 1]], gaussEqs[[All, 2]]}], {"", "drift (dt)", "noise (dW)"}],
    Frame -> All, Alignment -> Left, Spacings -> {2, 0.8}], 14, FontFamily -> "Helvetica"], Top]
```

The noise of all three widths is zero for every parameter: the four runs put the packet in four places, so the means fan out, but they give it the same widths at the same time, so the purity curves fall on one line. Only the means feel the readings. The widths obey the Riccati equation in the drift column, where the efficiency scales the information terms $\eta\lambda$ while the heating $\lambda\hbar^2$ keeps the full rate, and the purity settles at $\sqrt\eta$ from any start.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Continuous error correction: four numbers track the syndrome

Three qubits carry the bit-flip code, and each qubit flips at rate $\gamma$. The stabilizers $S_1=Z_1Z_2$ and $S_2=Z_2Z_3$ are measured weakly and continuously at rate $\kappa$ by detectors of efficiency $\eta$. Can four numbers, fed only these two noisy records, track the syndrome?

With $\mathcal D[c]\rho=c\rho c^\dagger-\frac12\{c^\dagger c,\rho\}$, the three-qubit state conditioned on the records $dJ_k=2\sqrt{\eta\kappa}\,\langle S_k\rangle dt+dW_k$ obeys

$$
d\rho=\gamma\sum_{i=1}^3\mathcal D[X_i]\rho\,dt+\kappa\sum_{k=1}^2\mathcal D[S_k]\rho\,dt+\sqrt{\eta\kappa}\sum_{k=1}^2\big(S_k\rho+\rho S_k-2\langle S_k\rangle\rho\big)dW_k .
$$

The four syndrome sectors, with projectors $\Pi_{s}=\frac14(1+s_1S_1)(1+s_2S_2)$, have probabilities $p_s$. The filter relaxes them toward $1/4$ under the flips and multiplies each by the likelihood of the new pair of readings:

$$
p_s\;\to\;\frac{\big(\tfrac14+e^{-4\gamma dt}(p_s-\tfrac14)\big)\,e^{2\sqrt{\eta\kappa}\,(s_1dJ_1+s_2dJ_2)}}{\sum_{s'}(\cdots)} .
$$

WL integrates the $8\times8$ equation from $|000\rangle$ with the positivity-preserving step of Rouchon and Ralph, and runs the filter on the same records; panel titles give the largest gap between the filter and the sector populations $\mathrm{Tr}[\Pi_s\rho]$ of the full state, and time is in units of $1/\kappa$:

```wl
qecOp[o_, k_] := KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2], 3], k -> o];
qecStabs = {qecOp[PauliMatrix[3], 1] . qecOp[PauliMatrix[3], 2], qecOp[PauliMatrix[3], 2] . qecOp[PauliMatrix[3], 3]};
qecFlips = qecOp[PauliMatrix[1], #] & /@ Range[3];
qecSigns = {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
qecProj = (IdentityMatrix[8] + #[[1]] qecStabs[[1]]) . (IdentityMatrix[8] + #[[2]] qecStabs[[2]])/4 & /@ qecSigns;
qecStep[γ_, η_, dt_][ρ_, dJ_] := With[{σ = Sqrt[η] dJ . qecStabs},
   With[{M = (1 - (2 + 3 γ) dt/2 - η dt) IdentityMatrix[8] + σ + σ . σ/2},
    #/Tr[#] &[M . ρ . M + γ dt Total[# . ρ . # & /@ qecFlips] + (1 - η) dt Total[# . ρ . # & /@ qecStabs]]]];
qecRun[γ_, η_, dt_, n_] := Rest @ FoldList[Function[{state, kick}, With[{dJ = 2 Sqrt[η] (Tr[# . First[state]] & /@ qecStabs) dt + kick},
      {qecStep[γ, η, dt][First[state], dJ], dJ}]], {Outer[Times, #, #] &[UnitVector[8, 1]], {0., 0.}}, RandomVariate[NormalDistribution[0, Sqrt[dt]], {n, 2}]];
qecFilter[γ_, η_, dt_, record_] := FoldList[#/Total[#] &[(1/4 + Exp[-4 γ dt] (#1 - 1/4)) Exp[2 Sqrt[η] qecSigns . #2]] &, {1., 0., 0., 0.}, record];
qecPanel[η_, γ_] := With[{run = qecRun[γ, η, 0.01, 4000]}, With[{filter = qecFilter[γ, η, 0.01, run[[All, 2]]],
     full = Function[state, Re[Tr[# . First[state]] & /@ qecProj]] /@ run},
    ListLinePlot[Transpose[filter], DataRange -> {0, 40}, PlotRange -> {-0.02, 1.02}, Frame -> True, FrameLabel -> {"κ t", "syndrome probability"},
     PlotLabel -> Column[{Row[{"η = ", η, ",   γ/κ = ", γ}], Row[{"largest gap to the full state: ", NumberForm[Max[Abs[Rest[filter] - full]], 2]}]}, Center],
     PlotLegends -> (Row[{"(", Subscript["S", "1"], ", ", Subscript["S", "2"], ") = (", #[[1]], ", ", #[[2]], ")"}] & /@ qecSigns), ImageSize -> 300]]];
qecPanels = BlockRandom[SeedRandom[5]; Outer[qecPanel, {1, 0.2}, {0.01, 0.05, 0.2}]];
Legended[Grid[Map[First, qecPanels, {2}]], Last[qecPanels[[1, 1]]]]
```

In every panel the filter stays within the stated gap of the sector populations of the full eight-dimensional state: the four probabilities need nothing else from the state. With a perfect detector and rare flips the filter jumps cleanly from sector to sector as each flip happens. The records need a time of order $1/(\eta\kappa)$ to settle on a sector; once flips arrive faster than that, as with a poor detector and frequent flips, the filter never commits.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Markovian feedback that protects the code space, and what the detector costs

Three qubits carry the bit-flip code, $|0_L\rangle=|000\rangle$, $|1_L\rangle=|111\rangle$, and each qubit flips at rate $\gamma$ through the channel $c_i=\sqrt\gamma X_i$. Each channel is monitored by homodyne detection of efficiency $\eta$, and its current $I_i$ drives the qubits at once through the Hamiltonian $I_iF_i$, with $F_i=i\sqrt\gamma\,[\Pi_{\rm code},X_i]$ and $\Pi_{\rm code}=|000\rangle\langle000|+|111\rangle\langle111|$. How well does this loop hold an encoded state?

With $\hbar=1$ and $\mathcal D[c]\rho=c\rho c^\dagger-\frac12\{c^\dagger c,\rho\}$, the Wiseman-Milburn feedback master equation is

$$
\dot\rho=-\frac{i}{2}\sum_{i=1}^3\big[c_i^\dagger F_i+F_ic_i,\rho\big]+\sum_{i=1}^3\mathcal D[c_i-iF_i]\rho+\frac{1-\eta}{\eta}\sum_{i=1}^3\mathcal D[F_i]\rho .
$$

WL checks the two identities $c_iF_i+F_ic_i=0$, which removes the Hamiltonian term, and $(c_i-iF_i)\Pi_{\rm code}=0$ (top). QF then evolves the master equation from $|0_L\rangle$ and from $|+_L\rangle=(|000\rangle+|111\rangle)/\sqrt2$, with $\mathcal D[c_i]$ alone as the case without feedback; time is in units of $1/\gamma$:

```wl
fbOp[o_, k_] := KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2], 3], k -> o];
fbFlips = fbOp[PauliMatrix[1], #] & /@ Range[3];
fbCode = DiagonalMatrix[UnitVector[8, 1] + UnitVector[8, 8]];
fbF = I (fbCode . # - # . fbCode) & /@ fbFlips;
fbChecks = {Max[Abs[MapThread[#1 . #2 + #2 . #1 &, {fbFlips, fbF}]]], Max[Abs[(# . fbCode) & /@ (fbFlips - I fbF)]]};
fbRho[ops_, rates_, ψ_] := Normal[QuantumEvolve[QuantumOperator[ConstantArray[0, {8, 8}]], QuantumOperator /@ ops -> rates, QuantumState[ψ], {t, 0, 3}]["DensityMatrix"]];
fbEtas = {1, 1/2, 1/10, 1/50};
fbPanels = KeyValueMap[Function[{name, ψ}, With[{ρs = Append[fbRho[Join[fbFlips - I fbF, fbF], Join[{1, 1, 1}, ConstantArray[(1 - #)/#, 3]], ψ] & /@ fbEtas, fbRho[fbFlips, {1, 1, 1}, ψ]]},
     MapThread[Plot[Evaluate[#1 /@ ρs], {t, 0, 3}, PlotRange -> {0, 1.02}, Frame -> True, FrameLabel -> {"γ t", #2}, PlotLabel -> Row[{"start ", name}], ImageSize -> 300,
        PlotStyle -> {Automatic, Automatic, Automatic, Automatic, Dashed}, PlotLegends -> Append[Row[{"η = ", #}] & /@ fbEtas, "no feedback"]] &,
      {{Re[Tr[fbCode . #]] &, Re[Conjugate[ψ] . # . ψ] &}, {"code-space population", "fidelity with the start"}}]]],
   <|Row[{"|", Subscript["0", "L"], "⟩"}] -> UnitVector[8, 1], Row[{"|", Subscript["+", "L"], "⟩"}] -> (UnitVector[8, 1] + UnitVector[8, 8])/Sqrt[2]|>];
Labeled[Legended[Grid[Map[First, fbPanels, {2}]], Last[fbPanels[[1, 1]]]],
  Style[Row[{Subscript["c", "i"], Subscript["F", "i"], " + ", Subscript["F", "i"], Subscript["c", "i"], " = ", fbChecks[[1]], "   and   (",
     Subscript["c", "i"], " - i ", Subscript["F", "i"], ")", Subscript["Π", "code"], " = ", fbChecks[[2]], "   for every i"}], 14, FontFamily -> "Helvetica"], Top]
```

At $\eta=1$ the two identities make every encoded state a fixed point: the code-space population and the fidelity stay at one, for $|0_L\rangle$ and $|+_L\rangle$ alike. Below unit efficiency the $\mathcal D[F_i]$ term feeds detector noise back onto the qubits and the plateau falls with $\eta$. At low efficiency the loop empties the code space faster than no feedback, though it settles slightly above it. From $|0_L\rangle$ the fidelity falls below the code-space population, but from $|+_L\rangle$, which a logical flip leaves unchanged, the two agree: the noisy loop also carries $|0_L\rangle$ into $|1_L\rangle$, which no stabilizer can see.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Syndrome populations of a repetition code: a closed equation from a general state

Three qubits carry the bit-flip code with stabilizers $S_1=Z_1Z_2$ and $S_2=Z_2Z_3$, each qubit flips at rate $\gamma$, and a field may drive the first qubit through $H=\Omega X_1$. The four syndrome sectors, with projectors $\Pi_s=\frac14(1+s_1S_1)(1+s_2S_2)$, have populations $p_s=\mathrm{Tr}[\Pi_s\rho]$. Do these four numbers obey a closed equation of their own?

With $\hbar=1$, the three-qubit state obeys the Lindblad equation

$$
\dot\rho=-i\,\Omega\,[X_1,\rho]+\gamma\sum_{k=1}^3\big(X_k\rho X_k-\rho\big).
$$

WL puts the question to a general $8\times8$ matrix $\rho$ with 64 symbolic entries: `SolveAlways` looks for a $4\times4$ matrix $A$ with $\mathrm{Tr}[\Pi_s\dot\rho]=\sum_{s'}A_{ss'}p_{s'}$ for every entry, once for the flip term and once for the drive term (top). QF evolves the full equation from $|000\rangle$ and from a superposition of two sectors (solid), against $p(t)=e^{\gamma At}p(0)$ from the flip term alone (dashed); time is in units of $1/\gamma$:

```wl
synOp[o_, k_] := KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2], 3], k -> o];
synStabs = {synOp[PauliMatrix[3], 1] . synOp[PauliMatrix[3], 2], synOp[PauliMatrix[3], 2] . synOp[PauliMatrix[3], 3]};
synProj = (IdentityMatrix[8] + #[[1]] synStabs[[1]]) . (IdentityMatrix[8] + #[[2]] synStabs[[2]])/4 & /@ {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
synRho = Array[synEntry, {8, 8}];
synClosure[drho_] := SolveAlways[Thread[(Tr[# . drho] & /@ synProj) == Array[synA, {4, 4}] . (Tr[# . synRho] & /@ synProj)], Flatten[synRho]];
synFlipA = Array[synA, {4, 4}] /. First[synClosure[Sum[synOp[PauliMatrix[1], k] . synRho . synOp[PauliMatrix[1], k] - synRho, {k, 3}]]];
synDriveA = synClosure[-I (synOp[PauliMatrix[1], 1] . synRho - synRho . synOp[PauliMatrix[1], 1])];
synPanel[{name_, ψ0_}, Ω_] := With[{ρt = Normal[QuantumEvolve[QuantumOperator[Ω synOp[PauliMatrix[1], 1]],
       QuantumOperator[synOp[PauliMatrix[1], #]] & /@ Range[3] -> {1, 1, 1}, QuantumState[ψ0], {t, 0, 1.5}]["DensityMatrix"]]},
   Show[Plot[Evaluate[Re[Tr[# . ρt]] & /@ Take[synProj, 3]], {t, 0, 1.5}, PlotRange -> All,
     PlotLegends -> {"(S1, S2) = (1, 1)", "(S1, S2) = (-1, 1)", "(S1, S2) = (-1, -1) and (1, -1): equal by symmetry"}],
    Plot[Evaluate[Take[MatrixExp[synFlipA t] . (Re[Conjugate[ψ0] . # . ψ0] & /@ synProj), 3]], {t, 0, 1.5}, PlotStyle -> Dashed, PlotRange -> All],
    PlotRange -> {0, 1}, Axes -> False, Frame -> True, FrameLabel -> {"γ t", "sector population"},
    PlotLabel -> Row[{"start ", name, ",   Ω/γ = ", Ω}], ImageSize -> 300]];
synPanels = Outer[synPanel, {{"|000⟩", UnitVector[8, 1]}, {"(|000⟩ + i|100⟩)/√2", (UnitVector[8, 1] + I UnitVector[8, 5])/Sqrt[2]}}, {0, 1, 4}, 1];
Labeled[Legended[Grid[Map[First, synPanels, {2}]], Last[synPanels[[1, 1]]]],
  Style[Row[{"flip term:  A = ", TraditionalForm[MatrixForm[synFlipA]], "        drive term:  ", If[synDriveA === {}, "no closed equation", synDriveA]}], 14, FontFamily -> "Helvetica"], Top]
```

The flip term closes: all 64 entries go in and only the four populations come out, with $\dot p_s=\gamma\big(\sum_{s'}p_{s'}-4p_s\big)=\gamma(1-4p_s)$, so every sector relaxes to $1/4$ at rate $4\gamma$ whatever coherences the start carries. The drive term has no solution: $X_1$ turns coherences between sectors into populations, so once $\Omega\neq0$ the full evolution leaves the dashed curves, soonest from the superposition.
