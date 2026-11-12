```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Two qubits entangled by a measurement alone

Two qubits that never interact are read by one probe, a half-parity meter: its output $V$, averaged over a window $\Delta t$, is Gaussian about $-\delta v$ for $|00\rangle$, $+\delta v$ for $|11\rangle$ and $0$ for both $|01\rangle$ and $|10\rangle$, with white noise of spectral density $s$, and the meter also dephases the odd pair at rate $\gamma$. Starting from $|+\rangle|+\rangle$, can the record alone entangle the pair?

A window leaves a reading of variance $D=s/(2\Delta t)$. Coherences between states the meter tells apart are assumed to dephase fast, so the state keeps five numbers $x=(\rho_{00,00},\rho_{01,01},\rho_{10,10},\rho_{11,11},|\rho_{01,10}|)$, and Korotkov's quantum Bayes rule updates them for each reading:

$$
x_i\to\frac{x_iP_i(V)}{\sum_{j=1}^4x_jP_j(V)}\ (i\le4),\qquad x_5\to\frac{x_5P_2(V)\,e^{-\gamma\Delta t}}{\sum_{j=1}^4x_jP_j(V)},\qquad P_i(V)=e^{-(V-v_i)^2/2D},\quad v=(-\delta v,0,0,\delta v).
$$

On the states this rule produces, the Wootters concurrence is $C=2\max\big(0,\,x_5-\sqrt{x_1x_4}\big)$.

WL checks the update laws with every symbol free (top), derives $C$ after a record of running mean $\bar V$, and runs five records per panel as lines, with that closed form at each record's own $\bar V$ as markers and its ceiling $\bar V=0$ dashed; time is in units of $s/\delta v^2$:

```wl
hpmUpdate[x_, v_, dv_, var_, gdt_] := With[{p = Exp[-({-dv, 0, 0, dv} - v)^2/(2 var)]},
   Append[Take[x, 4] p, x[[5]] p[[2]] Exp[-gdt]]/(Take[x, 4] . p)];
hpmLaws = Block[{x1, x2, x3, x4, x5, v1, v2, dv, var, gdt}, With[{x = {x1, x2, x3, x4, x5}, u = hpmUpdate[{x1, x2, x3, x4, x5}, v1, dv, var, gdt]},
    FullSimplify[{Total[Take[u, 4]] == 1, u[[5]]^2/(u[[2]] u[[3]]) == Exp[-2 gdt] x5^2/(x2 x3),
      hpmUpdate[u, v2, dv, var, gdt] == hpmUpdate[x, (v1 + v2)/2, dv, var/2, 2 gdt]},
     Positive /@ {x1, x2, x3, x4, x5, var} && Element[{v1, v2, dv, gdt}, Reals] && x1 + x2 + x3 + x4 == 1]]];
hpmConc[vb_, t_] = Expand[Numerator[#] Exp[-δv^2 t/s]]/Expand[Denominator[#] Exp[-δv^2 t/s]] &[FullSimplify[
     2 (#[[5]] - Sqrt[#[[1]] #[[4]]]) &[hpmUpdate[{1, 1, 1, 1, 1}/4, vb, δv, s/(2 t), γ t]], s > 0 && t > 0 && δv > 0 && γ >= 0 && vb ∈ Reals]];
hpmRun[gs_, seed_] := BlockRandom[SeedRandom[seed]; NestList[With[{v = RandomChoice[Take[First[#], 4] -> {-1., 0., 0., 1.}] + RandomVariate[NormalDistribution[0, Sqrt[50.]]]},
      {hpmUpdate[First[#], v, 1., 50., gs/100], Last[#] + v}] &, {ConstantArray[0.25, 5], 0.}, 500]];
hpmPanel[gs_] := With[{runs = hpmRun[gs, #] & /@ Range[5], ceil = hpmConc[0, t] /. {δv -> 1, s -> 1, γ -> gs}},
   Show[ListLinePlot[Table[{k/100, 2 Ramp[#[[k + 1, 1, 5]] - Sqrt[#[[k + 1, 1, 1]] #[[k + 1, 1, 4]]]]}, {k, 0, 500}] & /@ runs, PlotRange -> All],
    ListPlot[Table[{k/100, Ramp[hpmConc[#[[k + 1, 2]]/k, k/100] /. {δv -> 1, s -> 1, γ -> gs}]}, {k, 25, 500, 25}] & /@ runs, PlotMarkers -> {Automatic, 9}, PlotRange -> All],
    Plot[Ramp[ceil], {t, 0, 5}, PlotStyle -> Directive[Dashed, Thick], PlotRange -> All],
    PlotRange -> {{0, 5}, {0, 1}}, Frame -> True, FrameLabel -> {Row[{Superscript["δv", 2], "t/s"}], "concurrence C"},
    PlotLabel -> Row[{"γ s/", Superscript["δv", 2], " = ", gs}], ImageSize -> 300]];
Labeled[Row[hpmPanel /@ {0, 0.3, 1.5}],
  Style[Column[{Row[{"trace kept, odd-pair purity ratio times ", Superscript["e", "-2γΔt"], ", two readings = one at their mean with variance D/2:   ", hpmLaws}],
     Row[{"C(", OverBar["V"], ", t) = ", TraditionalForm[hpmConc[OverBar["V"], t]], "        C(0, t) = ", TraditionalForm[Normal[Series[hpmConc[0, t], {t, 0, 1}]]], " + O(", Superscript["t", 2], ")"}]}, Center],
   13, FontFamily -> "Helvetica"], Top]
```

The update keeps the trace, scales the odd pair's purity ratio $x_5^2/(x_2x_3)$ by the dephasing alone, and turns two readings into one at their mean, so each record's concurrence lands on the closed form at its own running mean. Records that settle near $\bar V=0$ press against the dashed ceiling, while records pushed toward $\pm\delta v$ lose their concurrence. The ceiling starts with slope $(\delta v^2-\gamma s)/2s$: without dephasing it climbs toward $C=1$, with weak dephasing it peaks and falls, and once $\gamma s>\delta v^2$ the numerator $e^{-\gamma t}-e^{-\delta v^2t/s}$ stays negative and no record entangles the pair.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### GHZ against product probes when every qubit dephases

A frequency $\omega$ is written on $n$ qubits, each dephasing independently at rate $\gamma$, by Ramsey runs of length $t$ repeated over a total time $T$. Without noise a GHZ probe has quantum Fisher information $n^2t^2$ against $nt^2$ for the product probe $|+\rangle^{\otimes n}$; does the GHZ advantage survive once $t$ is chosen optimally?

With $\hbar=1$ and $\mathcal D[c]\rho=c\rho c^\dagger-\frac12\{c^\dagger c,\rho\}$, each run evolves under

$$
\dot\rho=-i\Big[\frac\omega2\sum_{k=1}^nZ_k,\rho\Big]+\gamma\sum_{k=1}^n\mathcal D[Z_k]\rho ,
$$

which damps each qubit's coherence as $e^{-2\gamma t}$. The $T/t$ runs carry information $(T/t)\,F_Q(t)$, with $F_Q=2\sum_{\lambda_i+\lambda_j>0}|\langle i|\partial_\omega\rho|j\rangle|^2/(\lambda_i+\lambda_j)$ over the eigenpairs $(\lambda_i,|i\rangle)$ of $\rho(t)$; the noise commutes with the phase, so $F_Q$ does not depend on $\omega$ and is taken at $\omega=0$.

QF evolves both probes with $\omega$, $\gamma$ and $t$ symbolic for $n=2,3,4$, WL takes $F_Q$ in closed form and maximizes $F_Q/t$ over $t$ (markers, and the table on top); rates are in units of $\gamma$:

```wl
dgzZ[n_, k_] := QuantumOperator["Z", {k}, n];
dgzRho[ψ_] := With[{n = ψ["Qudits"]}, Simplify[Normal[QuantumEvolve[Total[ω/2 dgzZ[n, #] & /@ Range[n]], dgzZ[n, #] & /@ Range[n] -> ConstantArray[γ, n], ψ, t]["DensityMatrix"]], γ > 0]];
dgzQFI[ρ_] := With[{λ = Eigenvalues[ρ /. ω -> 0], V = Orthogonalize[Eigenvectors[ρ /. ω -> 0], Dot]},
   With[{M = V . (D[ρ, ω] /. ω -> 0) . Transpose[V]},
    Simplify[2 Total[MapThread[If[PossibleZeroQ[#2], 0, Abs[#1]^2/#2] &, {M, Outer[Plus, λ, λ]}, 2], 2], γ > 0 && t > 0]]];
dgzF = Table[{dgzQFI[dgzRho[QuantumState["GHZ"[n]]]], dgzQFI[dgzRho[QuantumState[StringRepeat["+", n]]]]}, {n, 2, 4}];
dgzBest = Map[Maximize[{#/t /. γ -> 1, t > 0}, t] &, dgzF, {2}];
dgzPanel[i_] := Show[Plot[Evaluate[dgzF[[i]]/t /. γ -> 1], {t, 0, 0.6}, PlotRange -> All, PlotStyle -> {Automatic, Dashed}, PlotLegends -> {"GHZ", "product"}],
   ListPlot[List /@ ({t /. Last[#], First[#]} & /@ dgzBest[[i]]), PlotMarkers -> {Automatic, 10}],
   PlotRange -> {{0, 0.6}, {0, 0.5}}, Frame -> True, FrameLabel -> {"γ t", Row[{"γ ", Subscript["F", "Q"], "/t"}]}, PlotLabel -> Row[{"n = ", i + 1}], ImageSize -> 280];
dgzPanels = dgzPanel /@ Range[3];
Labeled[Legended[Row[First /@ dgzPanels], Last[First[dgzPanels]]],
  Style[Grid[Prepend[MapThread[{#1, TraditionalForm[#2[[1]]], TraditionalForm[#3[[1, 1]]], TraditionalForm[t /. #3[[1, 2]]], TraditionalForm[#2[[2]]], TraditionalForm[#3[[2, 1]]], TraditionalForm[t /. #3[[2, 2]]]} &, {Range[2, 4], dgzF, dgzBest}],
     {"n", Row[{"GHZ  ", Subscript["F", "Q"]}], Row[{"largest γ ", Subscript["F", "Q"], "/t"}], "at γ t", Row[{"product  ", Subscript["F", "Q"]}], Row[{"largest γ ", Subscript["F", "Q"], "/t"}], "at γ t"}], Spacings -> {2, 0.6}], 13, FontFamily -> "Helvetica"], Top]
```

The GHZ coherence decays $n$ times faster than one qubit's, so the $n^2$ in its $F_Q$ comes with $e^{-4n\gamma t}$, while the product probe keeps $e^{-4\gamma t}$. For short runs the GHZ probe gathers information faster, but its best run is $n$ times shorter and reaches exactly the same largest $F_Q/t$ as the product probe: under independent dephasing the Heisenberg gain cancels at the optimum.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The noise term in the syndrome equation of a watched code

Three qubits carry the bit-flip code, each qubit flips at rate $\gamma$, and the stabilizers $S_1=Z_1Z_2$, $S_2=Z_2Z_3$ are measured continuously at rate $\kappa$ by detectors of efficiency $\eta$. The four syndrome sectors $j$, with signs $s_{jk}=\pm1$ and projectors $\Pi_j=\frac14(1+s_{j1}S_1)(1+s_{j2}S_2)$, have probabilities $p_j=\mathrm{Tr}[\Pi_j\rho]$. Which term does the measurement record add to the equation for the $p_j$?

With $\mathcal D[c]\rho=c\rho c^\dagger-\frac12\{c^\dagger c,\rho\}$, records $dJ_k=2\sqrt{\eta\kappa}\,\langle S_k\rangle dt+dW_k$ and Wiener increments $dW_k$, the state conditioned on the records obeys

$$
d\rho=\gamma\sum_{i=1}^3\mathcal D[X_i]\rho\,dt+\kappa\sum_{k=1}^2\mathcal D[S_k]\rho\,dt+\sqrt{\eta\kappa}\sum_{k=1}^2\big(S_k\rho+\rho S_k-2\langle S_k\rangle\rho\big)dW_k ,
$$

and $dp_j=\mathrm{Tr}[\Pi_j\,d\rho]$.

WL never builds a matrix: it writes $\rho$ as the most general element of the Pauli algebra, with 64 free coefficients, and `SolveAlways` turns the trace of $\Pi_j$ against each term into a combination of the $p_j$ (top). It then integrates the resulting four-number equation from the code sector, three records and the mean of 1000 against the mean from the flip term alone; time is in units of $1/\kappa$:

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
wsyDrift = MapIndexed[Function[{j, i}, Simplify[γ wsyInPops[Total[wsyTr[wsyProj[j] ** (# ** wsyRho ** # - wsyRho)] & /@ {wsyX1, wsyX2, wsyX3}]] /.
      First[Solve[Total[wsyP] == 1, wsyP[[Mod[First[i], 4] + 1]]]]]], wsySigns];
wsyBack = Function[j, κ wsyInPops[Total[wsyTr[wsyProj[j] ** (# ** wsyRho ** # - wsyRho)] & /@ wsyStab]]] /@ wsySigns;
wsyMean = wsyInPops[wsyTr[# ** wsyRho]] & /@ wsyStab;
wsyRecord = Outer[Sqrt[η κ] (wsyInPops[wsyTr[wsyProj[#1] ** (#2 ** wsyRho + wsyRho ** #2)]] - 2 wsyInPops[wsyTr[#2 ** wsyRho]] wsyInPops[wsyTr[wsyProj[#1] ** wsyRho]]) &, wsySigns, wsyStab, 1];
wsyCheck = Simplify[wsyRecord == 2 Sqrt[η κ] wsyP (Transpose[Transpose[wsySigns] - wsyMean])];
wsyRun[e_, g_, m_] := BlockRandom[SeedRandom[3]; NestList[Nest[With[{dW = RandomVariate[NormalDistribution[0, Sqrt[0.001]], {m, 2}]},
        # + (wsyDrift /. {γ -> g, Subscript[p, j_] :> #[[All, j]]} // Transpose) 0.001 + 2 Sqrt[e] # (dW . Transpose[wsySigns] - Total[(# . wsySigns) dW, {2}])] &, #, 50] &,
     ConstantArray[{1., 0., 0., 0.}, m], 400]];
wsyPanel[e_, g_] := With[{run = wsyRun[e, g, 1000][[All, All, 1]],
    flipsOnly = DSolveValue[{y'[τ] == (First[wsyDrift] /. {γ -> g, Subscript[p, 1] -> y[τ]}), y[0] == 1}, y, τ]},
   ListLinePlot[Join[Transpose[run[[All, 1 ;; 3]]], {Mean /@ run, flipsOnly /@ Subdivide[0, 20, 400]}], DataRange -> {0, 20},
    PlotStyle -> {Thin, Thin, Thin, Thick, Dashed}, PlotLegends -> {"record 1", "record 2", "record 3", "mean of 1000 records", "mean from the flip term alone"},
    PlotRange -> {{0, 20}, {-0.02, 1.02}}, Frame -> True, FrameLabel -> {"κ t", Subscript["p", "code"]}, PlotLabel -> Row[{"η = ", e, ",   γ/κ = ", g}], ImageSize -> 300]];
wsyPanels = Outer[wsyPanel, {1, 0.2}, {0.02, 0.1}];
Labeled[Legended[Grid[Map[First, wsyPanels, {2}]], Last[wsyPanels[[1, 1]]]],
  Style[Column[{Row[{"flips:  ", wsyDrift, "        stabilizer back-action:  ", wsyBack}],
     Row[{"record, coefficient of d", Subscript["W", "k"], " equals 2", Sqrt["η κ"], " ", Subscript["p", "j"], "(", Subscript["s", "jk"], " - ", Subscript[OverBar["s"], "k"], "),  with  ",
       Subscript[OverBar["s"], "k"], " = Tr[", Subscript["S", "k"], "ρ] = ", wsyMean, ":   ", wsyCheck}]}, Center], 13, FontFamily -> "Helvetica"], Top]
```

The record adds $2\sqrt{\eta\kappa}\,p_j\sum_k(s_{jk}-\bar s_k)\,dW_k$: each sector is pushed by its own probability times how far its signs sit from the running means $\bar s_k$, so the four numbers still need nothing else from the state, and the back-action of the stabilizer measurement moves no probability at all. The records jump between near 0 and near 1, and the efficient detector commits more sharply than the weak one, while their mean follows the dashed relaxation set by the flips alone, because the noise term has zero mean.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A charge qubit read by a point contact: when the Rabi peaks merge

One electron tunnels between two quantum dots, left $|L\rangle=|0\rangle$ and right $|R\rangle=|1\rangle$, with tunnelling $\Omega$ and level offset $\varepsilon$, while a nearby quantum point contact passes a current that depends on which dot is occupied. The current noise carries the fluctuations of $\sigma_z=|L\rangle\langle L|-|R\rangle\langle R|$. As the contact reads more strongly, do the Rabi peaks of that noise spectrum merge where the Liouvillian's eigenvectors coalesce?

With $\hbar=1$, $n_L=|L\rangle\langle L|$ and $\mathcal D[c]\rho=c\rho c^\dagger-\frac12\{c^\dagger c,\rho\}$,

$$
\dot\rho=\mathcal L\rho=-i\Big[\frac\varepsilon2\sigma_z+\frac\Omega2\sigma_x,\rho\Big]+\kappa\,\mathcal D[n_L]\rho ,\qquad k=\kappa/4 ,
$$

since $\mathcal D[n_L]=\frac14\mathcal D[\sigma_z]$. By quantum regression, the symmetrized spectrum is $S(\omega)=2\,\mathrm{Re}\,\mathrm{Tr}\big[\sigma_z\,(i\omega-\mathcal L)^{-1}\tfrac12\{\sigma_z-\langle\sigma_z\rangle_{\rm ss},\rho_{\rm ss}\}\big]$ with steady state $\rho_{\rm ss}$. A side peak exists when $dS/d\omega=0$ has a root at $\omega>0$; an exceptional point is a $k$ where two nonzero eigenvalues of $\mathcal L$ coincide.

QF builds $\mathcal L$ with $\varepsilon$, $\Omega$ and $k$ symbolic, WL inverts it for $S(\omega)$ in closed form, finds with `Reduce` the largest $k$ that keeps side peaks and the $k$ where the discriminant of $\mathcal L$'s characteristic polynomial vanishes (panel titles, and the resonant results on top); frequencies and $k$ are in units of $\Omega$:

```wl
qpcL = Normal[QuantumOperator["Liouvillian"[ε/2 QuantumOperator["Z"] + Ω/2 QuantumOperator["X"], {QuantumOperator[{{1, 0}, {0, 0}}]}, {4 k}]]["Matrix"]];
qpcSS = With[{v = ArrayReshape[First[NullSpace[qpcL]], {2, 2}]}, v/Tr[v]];
qpcSpec = With[{dz = PauliMatrix[3] - Tr[PauliMatrix[3] . qpcSS] IdentityMatrix[2]},
    Simplify[2 ComplexExpand[Re[Tr[PauliMatrix[3] . ArrayReshape[LinearSolve[I ω IdentityMatrix[4] - qpcL, Flatten[(dz . qpcSS + qpcSS . dz)/2]], {2, 2}]]]]]];
qpcSideEdge[εv_] := MaxValue[{k, Reduce[Exists[u, u > 0 && Numerator[Together[D[qpcSpec /. {ω -> Sqrt[u], Ω -> 1, ε -> εv}, u]]] == 0] && k > 0, k]}, k];
qpcEP[εv_] := Cases[NSolve[Discriminant[Cancel[CharacteristicPolynomial[qpcL /. {Ω -> 1, ε -> εv}, λ]/λ], λ] == 0 && k > 0, k, Reals], {k -> v_} :> v];
qpcPanel[εv_] := Plot[Evaluate[qpcSpec /. {Ω -> 1, ε -> εv, k -> #} & /@ {0.25, 0.5, 0.75, 1, 1.5}], {ω, -3, 3}, PlotRange -> All,
   PlotLegends -> (Row[{"k/Ω = ", #}] & /@ {0.25, 0.5, 0.75, 1, 1.5}), Frame -> True, FrameLabel -> {"ω/Ω", "Ω S(ω)"},
   PlotLabel -> Column[{Row[{"ε/Ω = ", εv}], Row[{"side peaks for k/Ω < ", NumberForm[N[qpcSideEdge[Rationalize[εv]]], 3]}],
       Row[{"exceptional points at k/Ω = ", If[qpcEP[Rationalize[εv]] === {}, "none", NumberForm[#, 3] & /@ qpcEP[Rationalize[εv]]]}]}, Center], ImageSize -> 300];
qpcPanels = qpcPanel /@ {0, 0.3, 1};
Labeled[Legended[Row[First /@ qpcPanels], Last[First[qpcPanels]]],
  Style[Column[{Row[{"ε = 0:   S(ω) = ", TraditionalForm[Simplify[qpcSpec /. ε -> 0]]}],
     Row[{"side peaks while   ", Simplify[Reduce[SeriesCoefficient[qpcSpec /. ε -> 0, {ω, 0, 2}] > 0 && k > 0 && Ω > 0, k], Ω > 0],
       "        exceptional point at   ", Simplify[Reduce[Discriminant[Cancel[CharacteristicPolynomial[qpcL /. ε -> 0, λ]/λ], λ] == 0 && k > 0 && Ω > 0, k], Ω > 0]}]}, Center],
   13, FontFamily -> "Helvetica"], Top]
```

On resonance the two peaks slide inward and meet at zero frequency at $k=\Omega/\sqrt2$, while the eigenvectors of $\mathcal L$ coalesce only at the exceptional point $k=\Omega$: at $k/\Omega=0.75$ the spectrum already shows one central peak. An offset raises the center of the spectrum and removes the side peaks at a smaller $k$; at $\varepsilon/\Omega=1$ they vanish although no exceptional point exists at all, so the merging of the peaks is a property of the line shape, not of eigenvector coalescence.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### An atom in a thermal cavity field: where the populations settle

A two-level atom, frequency $\omega_a$, starts with excited-state population $c$ in a cavity of frequency $\omega_c$ whose field is thermal with mean photon number $\bar n$, and couples to it at strength $g$. On average over long times, what ratio of excited to ground population does the atom settle at?

With $\hbar=1$, $|e\rangle=|0\rangle$, $\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|$ and $\Delta=\omega_a-\omega_c$, the Jaynes-Cummings Hamiltonian and its form in the frame rotating at $\omega_c$ are

$$
H=\omega_c\,a^\dagger a+\frac{\omega_a}{2}\sigma_z+g\left(\sigma_-a^\dagger+\sigma_+a\right)\;\longrightarrow\;\frac\Delta2\sigma_z+g\left(\sigma_-a^\dagger+\sigma_+a\right),\qquad \rho_{\rm field}=\sum_{n\ge0}\frac{\bar n^{\,n}}{(\bar n+1)^{n+1}}|n\rangle\langle n| .
$$

$H$ couples $|e,n\rangle$ only to $|g,n+1\rangle$, and the long-time average keeps only the part of $\rho(0)$ diagonal in the eigenbasis of each such pair. An atomic coherence links $|e,n\rangle$ to $|g,n\rangle$, which sit in different pairs, so it drops out, and the atom starts as $c|e\rangle\langle e|+(1-c)|g\rangle\langle g|$.

WL averages each pair's populations from its eigenvectors and sums them over the thermal distribution in closed form (curves); QF time-averages $\rho(0)$ in the eigenbasis of $H$ on 100 Fock levels with $\bar n$ kept a symbol (markers); $\Delta$ is in units of $g$:

```wl
jctFromE[Δ_, m_] = Simplify[Total[#[[1]]^4/(# . #)^2 & /@ Eigenvectors[{{Δ/2, Sqrt[m + 1]}, {Sqrt[m + 1], -Δ/2}}]], m >= 0 && Δ ∈ Reals];
jctFromG[Δ_, m_] = Simplify[Total[#[[1]]^2 #[[2]]^2/(# . #)^2 & /@ Eigenvectors[{{Δ/2, Sqrt[m]}, {Sqrt[m], -Δ/2}}]], m > 0 && Δ ∈ Reals];
jctPe[nb_, c_, Δ_] := Sum[nb^m/(nb + 1)^(m + 1) jctFromE[Δ, m], {m, 0, Infinity}] c + Sum[nb^m/(nb + 1)^(m + 1) jctFromG[Δ, m], {m, 1, Infinity}] (1 - c);
jctRatio[nb_, c_, Δ_] := jctPe[nb, c, Δ]/(1 - jctPe[nb, c, Δ]);
jctA = AnnihilationOperator[100, {2}]; jctLower = QuantumOperator["-"];
jctQF[Δ_] := With[{es = Eigensystem[N[Normal[(Δ/2 QuantumOperator["Z"] + jctLower @ SuperDagger[jctA] + SuperDagger[jctLower] @ jctA)["Matrix"]]]]},
   With[{w = Total[Abs[Transpose[es[[2, #]]] . Conjugate[es[[2, #]]]]^2 & /@ GatherBy[Range[200], Round[es[[1, #]], 10^-8] &]],
     ρ0 = Diagonal[Normal[QuantumTensorProduct[QuantumState[{{c, 0}, {0, 1 - c}}], ThermalState[nb, 100]]["DensityMatrix"]]]},
    #/(1 - #) &[Total[w[[1 ;; 100]]] . ρ0]]];
jctMarks = AssociationMap[jctQF, {0, 2, 5}];
jctPanel[cv_] := Show[Plot[Evaluate[jctRatio[nb, cv, #] & /@ {0, 2, 5}], {nb, 0, 6}, PlotRange -> All, PlotLegends -> (Row[{"Δ/g = ", #}] & /@ {0, 2, 5})],
   ListPlot[Table[{x, jctMarks[#] /. {nb -> x, c -> cv}}, {x, {0.5, 1, 2, 3, 4.5, 6}}] & /@ {0, 2, 5}, PlotMarkers -> {Automatic, 9}],
   PlotRange -> All, Frame -> True, FrameLabel -> {OverBar["n"], Row[{OverBar[Subscript["P", "e"]], "/", OverBar[Subscript["P", "g"]]}]},
   PlotLabel -> Row[{"c = ", InputForm[cv]}], ImageSize -> 290];
jctPanels = jctPanel /@ {1, 1/2, 0};
Labeled[Legended[Row[First /@ jctPanels], Last[First[jctPanels]]],
  Style[Row[{"Δ = 0:   ", OverBar[Subscript["P", "e"]], "/", OverBar[Subscript["P", "g"]], " = ", TraditionalForm[FullSimplify[jctRatio[OverBar["n"], c, 0], OverBar["n"] > 0]]}], 14, FontFamily -> "Helvetica"], Top]
```

On resonance the ratio is $(\bar n+c)/(\bar n+2-c)$, which for $c=\tfrac12$ is $(2\bar n+1)/(2\bar n+3)$: the vacuum term $|g,0\rangle$, which never evolves, holds the ground state above the excited one, and its weight $1/(\bar n+1)$ fades as the field warms, so every start approaches equal populations. An excited atom on resonance settles at equal populations for every $\bar n$, while detuning keeps the atom nearer its starting state and slows the approach.
