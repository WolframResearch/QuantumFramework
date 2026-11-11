```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A tetrahedron measurement on an unknown qubit

A qubit is in an unknown pure state $|\psi\rangle=\cos\frac\theta2|0\rangle+e^{i\phi}\sin\frac\theta2|1\rangle$, and a single detector with four outcomes measures it. What do the four outcome probabilities reveal about $\theta$ and $\phi$?

The detector is the symmetric informationally complete POVM of a qubit: four effects $E_i=\tfrac12|\psi_i\rangle\langle\psi_i|$ with $\sum_iE_i=I$ and equal overlaps $|\langle\psi_i|\psi_j\rangle|^2=\tfrac13$ for $i\neq j$, whose Bloch vectors are the vertices of a regular tetrahedron. QF's first vertex is $|\psi_1\rangle=|0\rangle$; the other three are $|\psi_k\rangle=(|0\rangle+\sqrt2\,e^{i\chi_k}|1\rangle)/\sqrt3$ with $(\chi_2,\chi_3,\chi_4)=(\tfrac{4\pi}3,\tfrac{2\pi}3,0)$. The Born rule gives $p_i=\operatorname{Tr}[\rho E_i]$ with $\rho=|\psi\rangle\langle\psi|$.

QF supplies the four effects; the Born rule is applied to each with $\theta$ and $\phi$ symbolic, and each probability is mapped over the sphere of states:

```wl
sicPsi = QuantumState[{Cos[θ/2], Exp[I φ] Sin[θ/2]}];
sicP = FullSimplify[ComplexExpand[Tr[sicPsi["DensityMatrix"] . Normal[#]] & /@ QuantumMeasurementOperator["TetrahedronSICPOVM"]["POVMElements"]]];
sicPanels = DensityPlot[#, {φ, 0, 2 Pi}, {θ, 0, Pi}, PlotRange -> {0, 1/2}, PlotLegends -> Automatic, FrameLabel -> {"φ", "θ"},
     FrameTicks -> {{Range[0, Pi, Pi/4], None}, {Range[0, 2 Pi, Pi/2], None}}, ImageSize -> 300] & /@ sicP;
sicLabel[i_] := Style[Row[{Subscript["p", i], " = ", TraditionalForm[sicP[[i]]]}], 13, FontFamily -> "Helvetica", LineBreakWithin -> False];
Labeled[Legended[Grid[Partition[MapIndexed[Labeled[First[#1], sicLabel[First[#2]], Top] &, sicPanels], 2], Spacings -> {2, 2}], Last[First[sicPanels]]],
 Style[Row[{Subscript["p", 1], " + ", Subscript["p", 2], " + ", Subscript["p", 3], " + ", Subscript["p", 4], " = ", Simplify[Total[sicP]]}], 13, FontFamily -> "Helvetica"], Top]
```

Each probability is $\tfrac14$ plus a linear function of the Bloch components $(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$, largest, $\tfrac12$, at its own vertex and zero at the opposite point, and the four sum to one for every state. The four vertex directions span space, so the probabilities fix the Bloch vector and with it $\theta$ and $\phi$: one measurement setting determines the state.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### State tomography from finite counts

An unknown qubit state $\rho=\tfrac12(I+\vec r\cdot\vec\sigma)$ is measured on $N$ copies in each of the three Pauli bases, and the state is estimated from the counts. How does the error of the estimate fall with $N$ for a maximally mixed and for a pure state?

In setting $k$ the counts $n_\pm$ give $\langle\sigma_k\rangle\approx(n_+-n_-)/N$, each an average of $N$ outcomes $\pm1$ with variance $1-\langle\sigma_k\rangle^2$. QF's estimator inverts the Born rule linearly and then maximizes the likelihood of the counts over all states, so the estimate is always a state, even when noise would put the linear inversion outside the Bloch ball.

QF draws the counts for a true $\vec r$ off every measurement axis, along $(4,3,8.5)$, estimates the state, and plots $\sqrt N\,|\vec r_{\rm est}-\vec r|$ averaged over 100 runs; a flat curve means the error falls as $N^{-1/2}$:

```wl
tomoOps = QuantumMeasurementOperator /@ {"X", "Y", "Z"};
tomoNs = {30, 300, 3000, 30000};
tomoErr[len_] := BlockRandom[SeedRandom[7]; With[{r = len Normalize[{4, 3, 8.5}]},
    Table[{n, Sqrt[n] Mean @ Table[Norm[QuantumStateEstimate[QuantumMeasurementSimulation[QuantumState["BlochVector"[r]], tomoOps, n]][
          "MaximumLikelihoodState"]["BlochVector"] - r], 100]}, {n, tomoNs}]]];
tomoLens = {0, 1};
ListPlot[tomoErr /@ tomoLens, ScalingFunctions -> {"Log", None}, Joined -> True, PlotMarkers -> Automatic, PlotRange -> {0, 2},
 PlotLegends -> (Row[{"|r| = ", #}] & /@ tomoLens), Frame -> True,
 FrameLabel -> {"shots per setting N", Row[{Sqrt["N"], " \[Times] mean Bloch-vector error"}]}, ImageSize -> 460]
```

Both curves are flat, so for both states the error falls as $N^{-1/2}$. The pure state is estimated more accurately at every $N$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Process tomography of a decaying qubit

A qubit relaxes from $|1\rangle$ to $|0\rangle$ with probability $\gamma$. The channel is treated as unknown: four probe states go through it, and each output is reconstructed from $N$ shots per Pauli setting. How does the reconstructed channel approach the true one as $N$ grows?

With $\sigma_z=|0\rangle\langle0|-|1\rangle\langle1|$, the amplitude-damping channel has Kraus operators $K_0=|0\rangle\langle0|+\sqrt{1-\gamma}\,|1\rangle\langle1|$ and $K_1=\sqrt\gamma\,|0\rangle\langle1|$. Any qubit channel acts on Bloch vectors as $\vec r\mapsto T\vec r+\vec c$. The probes $|0\rangle,|1\rangle,|{+}\rangle,|{+}i\rangle$ have Bloch vectors $\hat z,-\hat z,\hat x,\hat y$, so their outputs give

$$
\vec c=\tfrac12(\vec r_0+\vec r_1),\qquad T=\big(\vec r_+-\vec c,\ \ \vec r_{+i}-\vec c,\ \ \tfrac12(\vec r_0-\vec r_1)\big)\ \text{(columns)},
$$

with each output $\vec r$ estimated by maximum likelihood from the counts.

QF reconstructs $(T,\vec c)$ and draws the image of the input great circle in the $xz$ plane under the reconstructed map and under the exact channel:

```wl
ptProbes = QuantumState /@ {"0", "1", "Plus", "Right"};
ptOps = QuantumMeasurementOperator /@ {"X", "Y", "Z"};
ptMap[{r0_, r1_, rp_, ri_}] := With[{c = (r0 + r1)/2}, {Transpose[{rp - c, ri - c, (r0 - r1)/2}], c}];
ptCurve[{tm_, c_}] := (tm . {Sin[α], 0, Cos[α]} + c)[[{1, 3}]];
ptOut[γ_, n_] := BlockRandom[SeedRandom[n];
   QuantumStateEstimate[QuantumMeasurementSimulation[QuantumChannel["AmplitudeDamping"[γ]][#], ptOps, n]]["MaximumLikelihoodState"]["BlochVector"] & /@ ptProbes];
ptNs = {30, 300, 3000};
ptPanel[γ_] := Show[
   ParametricPlot[Evaluate[ptCurve[ptMap[Re[ptOut[γ, #]]]] & /@ ptNs], {α, 0, 2 Pi}, PlotRange -> All, PlotLegends -> (Row[{"N = ", #}] & /@ ptNs)],
   ParametricPlot[Evaluate[ptCurve[ptMap[QuantumChannel["AmplitudeDamping"[γ]][#]["BlochVector"] & /@ ptProbes]]], {α, 0, 2 Pi}, PlotStyle -> Dashed, PlotRange -> All, PlotLegends -> {"exact channel"}],
   ParametricPlot[{Sin[α], Cos[α]}, {α, 0, 2 Pi}, PlotStyle -> Dotted, PlotLegends -> {"input great circle"}],
   PlotRange -> {{-1.05, 1.05}, {-1.05, 1.05}}, AspectRatio -> 1, Axes -> False, Frame -> True, FrameLabel -> {"x", "z"}, PlotLabel -> Row[{"γ = ", γ}], ImageSize -> 300];
ptPanels = ptPanel /@ {0.2, 0.5, 0.8};
Legended[Row[ptPanels //. Legended[g_, _] :> g], Flatten[Cases[First[ptPanels], Legended[_, l_] :> l, {0, Infinity}]]]
```

The exact channel shrinks the circle to an ellipse pulled toward $|0\rangle$ at the top, flatter along $z$ than along $x$, and smaller as $\gamma$ grows. At $N=30$ the reconstruction is visibly off and in every panel carries some states outside the Bloch ball, although every output it was built from is a state; by $N=3000$ it is close to the exact ellipse in every panel.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Undoing a measurement while its record survives

A qubit in $|\psi\rangle=\cos\frac\theta2|0\rangle+e^{i\phi}\sin\frac\theta2|1\rangle$ is measured in the $Z$ basis, and the outcome is written on a pointer qubit instead of being read. Before the measurement is reversed, the environment partly decoheres the pointer, at strength $\lambda$. How much of $|\psi\rangle$ does the reversal return?

The measurement is the isometry $V|\psi\rangle=\sum_kP_k|\psi\rangle\otimes|k\rangle_{\rm P}$, $P_k=|k\rangle\langle k|$, with $V^\dagger V=I$, so $V^\dagger$ is its exact left inverse. Between the two, phase damping $\mathcal D_\lambda$ acts on the pointer with Kraus operators $|0\rangle\langle0|+\sqrt{1-\lambda}\,|1\rangle\langle1|$ and $\sqrt\lambda\,|1\rangle\langle1|$, which multiply the record's coherence by $\sqrt{1-\lambda}$: $\lambda=0$ keeps the record intact and $\lambda=1$ erases its coherence. From $\rho=|\psi\rangle\langle\psi|$ the returned state is $\rho_{\rm out}=V^\dagger\,\mathcal D_\lambda(V\rho V^\dagger)\,V$, with fidelity $F=\langle\psi|\rho_{\rm out}|\psi\rangle$.

QF takes the dilation $V$ from the measurement object, places its pointer on wire 2, runs measure, damp, un-measure as one circuit with $\theta$, $\phi$, $\lambda$ symbolic, and plots $F$ against the input's polar angle:

```wl
umPsi = QuantumState[{Cos[θ/2], Exp[I φ] Sin[θ/2]}];
umV = QuantumOperator[QuantumMeasurementOperator[{1}]["SuperOperator"], {{2, 1}, {1}}];
umRho = QuantumCircuitOperator[{umV, QuantumChannel["PhaseDamping"[λ], {2}], umV["Dagger"]}][umPsi]["DensityMatrix"];
umF = FullSimplify[ComplexExpand[Conjugate[Normal[umPsi["StateVector"]]] . Normal[umRho] . Normal[umPsi["StateVector"]]], {θ, φ} ∈ Reals && 0 <= λ <= 1];
umλs = {0, 0.36, 0.75, 1};
Plot[Evaluate[umF /. λ -> umλs], {θ, 0, Pi}, PlotRange -> {0.45, 1.02}, Frame -> True, FrameLabel -> {"θ", "fidelity F"},
 FrameTicks -> {{Automatic, None}, {{0, Pi/4, Pi/2, 3 Pi/4, Pi}, None}}, PlotLegends -> (Row[{"λ = ", #}] & /@ umλs),
 PlotLabel -> Style[Row[{"F = ", TraditionalForm[umF]}], 13, FontFamily -> "Helvetica"], ImageSize -> 420]
```

$F$ does not depend on $\phi$, and with the record intact the measurement is undone exactly for every input. As the record decoheres, the poles, which the $Z$ measurement does not disturb, still return intact while the equator loses most; at $\lambda=1$ the equatorial fidelity is $\tfrac12$, the value left by a completed measurement.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The Hamiltonian behind a chirped rotation

A qubit should follow the trajectory $U(t)=R_Z(\beta t^2)\,R_X(\Omega t)$: a rotation about $x$ at rate $\Omega$, followed by a $z$ rotation whose angle grows quadratically at chirp rate $\beta$. What control Hamiltonian drives it?

With $\hbar=1$, Pauli matrices $X,Y,Z$ and $R_Z(\vartheta)=e^{-i\vartheta Z/2}$, $R_X(\vartheta)=e^{-i\vartheta X/2}$, a propagator obeys $i\dot U=HU$, so the Hamiltonian is $H(t)=i\,\dot U(t)\,U^\dagger(t)$, written as $H=h_x X+h_y Y+h_z Z$.

QF differentiates the operator $U(t)$ with `D`, forms $i\dot UU^\dagger$, and reads the three Pauli coefficients, in units of $\Omega$ and with time in units of $1/\Omega$:

```wl
hfU = QuantumOperator["RZ"[β t^2]] @ QuantumOperator["RX"[Ω t]];
hfField = FullSimplify[ComplexExpand[#]] & /@ (I D[hfU, t] @ hfU["Dagger"])["PauliDecompose"];
hfPanel[b_] := Plot[Evaluate[Lookup[hfField, {"X", "Y", "Z"}] /. {Ω -> 1, β -> b}], {t, 0, 6}, PlotRange -> All,
   PlotLegends -> {Subscript["h", "x"], Subscript["h", "y"], Subscript["h", "z"]}, Frame -> True, FrameLabel -> {"Ω t", "coefficient / Ω"},
   PlotLabel -> Row[{"β/", Superscript["Ω", 2], " = ", b}], ImageSize -> 300];
hfPanels = hfPanel /@ {0.05, 0.2, 0.5};
Labeled[Legended[Row[First /@ hfPanels], Last[First[hfPanels]]],
 Style[Row[Riffle[MapThread[Row[{Subscript["h", #1], " = ", TraditionalForm[#2]}] &, {{"x", "y", "z"}, Lookup[hfField, {"X", "Y", "Z"}]}], ",    "]], 13, FontFamily -> "Helvetica"], Top]
```

The drive is a transverse field of fixed strength $\Omega/2$ whose direction turns in the $xy$ plane with phase $\beta t^2$, a linear frequency chirp, plus a $z$ field that grows linearly as $\beta t$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### An imported CHSH circuit reaches the Tsirelson value exactly

A referee sends random bits $x$ to Alice and $y$ to Bob, who share an entangled pair and answer with bits $a$ and $b$; they win when $x\wedge y=a\oplus b$, which no classical strategy does more than $\tfrac34$ of the time. The strategy arrives as an OpenQASM 3 program in which Bob's measurement is set by a rotation angle $\theta_B$. How does the winning probability depend on $\theta_B$?

The program prepares $(|00\rangle+|11\rangle)/\sqrt2$ on Alice's qubit $q_0$ and Bob's $q_3$, and puts $q_1,q_2$ in $|{+}\rangle$ so that measuring them supplies $x$ and $y$. Alice applies a Hadamard to $q_0$ when $x=1$; Bob applies $R_Y(\theta_B)=U(\theta_B,0,0)$ to $q_3$, then a Hadamard when $y=0$; all four qubits are measured in $Z$, giving $(a,x,y,b)$.

QF imports the program for $\theta_B=k\pi/8$, $k=0,\dots,16$, and reads the winning probability exactly from the circuit's joint outcome distribution:

```wl
chshProgram = "OPENQASM 3.0;
qubit[4] q;
bit[4] c;
reset q[0];
reset q[3];
U(pi/2, 0, pi) q[3];
U(pi/2, 0, pi) q[0];
swap q[0] q[3];
ctrl(0) @ negctrl(1) @ U(pi/2, pi, pi) q[0] q[3];
ctrl(1) @ negctrl(0) @ U(pi/2, 0, 0) q[0] q[3];
reset q[1];
U(pi/2, 0, pi) q[1];
reset q[2];
U(pi/2, 0, pi) q[2];
ctrl(1) @ negctrl(0) @ U(pi/2, 0, pi) q[1] q[0];
U(THETAB, 0, 0) q[3];
ctrl(0) @ negctrl(1) @ U(pi/2, 0, pi) q[2] q[3];
c[0] = measure q[0];
c[1] = measure q[1];
c[2] = measure q[2];
c[3] = measure q[3];";
chshWin = Table[{k/8, FullSimplify @ Probability[BitAnd[xBit, yBit] == BitXor[aBit, bBit], {aBit, xBit, yBit, bBit} \[Distributed]
       QuantumCircuitOperator[StringReplace[chshProgram, "THETAB" -> ToString[k] <> "*pi/8"]][]["MultivariateDistribution"]]}, {k, 0, 16}];
chshBest = First[MaximalBy[chshWin, N @* Last]];
Show[ListPlot[chshWin, PlotMarkers -> {Automatic, 9}, PlotLegends -> {"imported program"}],
 Plot[3/4, {k, 0, 2}, PlotStyle -> Dashed, PlotLegends -> {"classical bound 3/4"}],
 PlotRange -> {{-0.05, 2.05}, {0, 1}}, Frame -> True, FrameLabel -> {Row[{Subscript["θ", "B"], "/π"}], "winning probability"},
 PlotLabel -> Style[Row[{"largest: ", TraditionalForm[Last[chshBest]], "  at  ", Subscript["θ", "B"], "/π = ", First[chshBest]}], 13, FontFamily -> "Helvetica"], ImageSize -> 460]
```

Every angle stays an exact multiple of $\pi$ through the import, so each probability is an exact number. Of the angles shown, only those strictly between $0$ and $\pi/2$ beat the classical bound, and the largest value is the Tsirelson bound $\cos^2\frac\pi8$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### One site of an arbitrarily long Ising chain after a quench

A chain of spins with coupling $J$ in a transverse field $h$ starts with every spin up, $|0\cdots0\rangle$, and evolves under the transverse-field Ising Hamiltonian, Trotterized with step $\delta t$. What is the magnetization $\langle Z_c\rangle$ of a site deep inside the chain after $d$ steps?

With $\hbar=1$, $H=J\sum_iZ_iZ_{i+1}+h\sum_iX_i$, and one step is $\prod_iR_X^{(i)}(2h\,\delta t)\prod_iR_{ZZ}^{(i,i+1)}(2J\,\delta t)$ with $R_{ZZ}(\vartheta)=e^{-i\vartheta Z\otimes Z/2}$. After $d$ steps $\langle Z_c\rangle=\langle0|U^\dagger Z_cU|0\rangle$, and every gate not wholly inside the sites $c-d,\dots,c+d$ commutes with the evolved operator it meets and cancels against its conjugate, so a window of $2d+1$ sites already gives the value for a chain of any length.

QF builds $U^\dagger Z_cU$ on a window of $m$ sites centered on $c$ as one circuit and contracts it as a tensor network, for the $2d+1$-site window (filled) and a $2d+5$-site window (open), with $J\delta t=0.2$:

```wl
lcQuench[m_, d_, θx_] := QuantumCircuitOperator @ Join[{StringJoin @ ConstantArray["0", m]},
    Flatten @ Table[Join[("R"[0.4, "ZZ"] -> # &) /@ Partition[Range[m], 2, 1], ("RX"[θx] -> # &) /@ Range[m]], {d}]];
lcZ[m_, d_, θx_] := With[{u = lcQuench[m, d, θx]},
    Re @ (u /* QuantumOperator["Z", {(m + 1)/2}] /* u["Dagger"])[Method -> {"TensorNetwork", "Path" -> "Greedy"}]["Scalar"]];
lcPanel[hJ_] := ListPlot[Table[{0.2 d, lcZ[2 d + #, d, 0.4 hJ]}, {d, 1, 5}] & /@ {1, 5}, PlotRange -> {{0, 1.1}, {-1, 1.05}},
    PlotMarkers -> {{"\[FilledCircle]", 9}, {"\[EmptyCircle]", 18}}, PlotLegends -> {Row[{2, "d + 1 sites"}], Row[{2, "d + 5 sites"}]},
    Frame -> True, FrameLabel -> {"J t", Row[{"\[LeftAngleBracket]", Subscript["Z", "c"], "\[RightAngleBracket]"}]}, PlotLabel -> Row[{"h/J = ", hJ}], ImageSize -> 300];
lcPanels = lcPanel /@ {0.5, 1, 2};
Legended[Row[First /@ lcPanels], Last[First[lcPanels]]]
```

At every depth the two windows give the same magnetization, so the sites beyond the light cone do not enter and each point is exact for a bulk site of an arbitrarily long chain. The magnetization falls faster as $h/J$ grows, and at $h/J=2$ it changes sign within the times shown.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Stored numbers, declared bases and composition

An observable is specified by its outcomes $\pm1$ in a chosen measurement basis: QF stores the numbers $C=\operatorname{diag}(1,-1)$ together with that basis as the operator $R$. What does applying $Z$ after $R$ give?

If the columns of $F$ are the declared basis vectors, the stored matrix $C$ represents the map $A=F\,C\,F^{-1}$ in the computational basis; $F$ only needs to be invertible. Composition contracts the shared wire after both operators are written in one basis, so $Z\circ R$ is the product $ZA$. Contracting the stored numbers as they sit would give $ZC$ instead.

QF builds $R$ in four declared bases and lists the map it represents and the composite $Z\circ R$, with $C$ and $ZC$ on top:

```wl
semNames = {"{|0⟩, |1⟩}", "{|+⟩, |−⟩}", "{|+i⟩, |−i⟩}", "{|0⟩, |+⟩}"};
semBases = {"Computational", "PauliX", "PauliY", QuditBasis[{{1, 0}, {1, 1}/Sqrt[2]}]};
semRow[b_] := With[{r = QuantumOperator[{{1, 0}, {0, -1}}, b]},
   TraditionalForm[MatrixForm[Simplify[Normal[#]]]] & /@ {r["MatrixRepresentation"], QuantumOperator["Z"][r]["MatrixRepresentation"]}];
Style[Labeled[Grid[Prepend[MapThread[Prepend[semRow[#1], #2] &, {semBases, semNames}],
    {"declared basis", "represented A", "Z after R in QF"}], Frame -> All, Spacings -> {2, 1}],
  Row[{"C = ", TraditionalForm[MatrixForm[{{1, 0}, {0, -1}}]], ",   Z C = ", TraditionalForm[MatrixForm[PauliMatrix[3] . {{1, 0}, {0, -1}}]], " in every basis"}], Top], 13, FontFamily -> "Helvetica"]
```

The same numbers $C$ represent $Z$, $X$, $Y$, or, in the non-orthogonal basis, a map that is not Hermitian, and QF's composite is $Z$ times that map in every row. Multiplying the stored numbers gives the identity whatever the basis, which is the right answer only when the declared basis is the computational one.
