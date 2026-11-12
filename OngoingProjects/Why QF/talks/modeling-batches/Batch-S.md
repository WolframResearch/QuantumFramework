```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Exact beam-splitter settings for an exact interferometer

A lossless interferometer on $m$ modes is fixed by its unitary transfer matrix $W$, here the discrete Fourier matrices $F_3$ and $F_4$, with $(F_m)_{jk}=e^{2\pi i(j-1)(k-1)/m}/\sqrt m$. Any such $W$ can be built from $m(m-1)/2$ beam splitters and $m$ phase shifters. When $W$ is known exactly, are the splitter settings exact numbers too?

The network maps $a_k^\dagger\to\sum_j W_{jk}\,a_j^\dagger$. A splitter on modes $p,q$ acts there by the block

$$
B(\theta,\varphi)=\begin{pmatrix}\cos\theta & -e^{-i\varphi}\sin\theta\\ e^{i\varphi}\sin\theta & \cos\theta\end{pmatrix},
$$

and a phase shifter multiplies one mode by $e^{i\alpha}$. Each splitter of the mesh is chosen to null one entry of $W$: for the two entries $a,b$ it mixes, $\tan\theta=|b|/|a|$ and $\varphi=\arg(b\,a^*)$.

QF decomposes $F_3$ and $F_4$ into its default (Clements) mesh, reduces the settings, rebuilds the mesh and compares it with $F_m$; the top label applies the nulling formula to a splitter with symbolic angles $\alpha,\beta$:

```wl
reckSettings[m_] := FullSimplify[FullSimplify[QuantumInterferometer[FourierMatrix[m]]["Elements"]] /. Arg[r_Root] :> Arg[ToRadicals[RootReduce[r/Abs[r]]]]] /.
   {"BS"[t_, f_] :> "BS"[t, Mod[f, 2 π]], "PS"[f_] :> "PS"[Mod[f, 2 π]]};
reckRow["BS"[t_, f_] -> q_] := {"splitter", Row[q, ", "], TraditionalForm[t], TraditionalForm[f]};
reckRow["PS"[f_] -> q_] := {"phase", Row[q, ", "], "", TraditionalForm[f]};
reckPanel[m_] := With[{el = reckSettings[m]}, Labeled[
    Column[{QuantumInterferometer[el, m]["Diagram", ImageSize -> 320],
      Style[Grid[Prepend[reckRow /@ el, {"element", "modes", "θ", "φ or α"}], Frame -> All, Spacings -> {1.5, 0.6}], 13, FontFamily -> "Helvetica"]}, Center, Spacer[10]],
    Style[Row[{"m = ", m, ":   rebuilt mesh equals ", Subscript["F", m], ":  ", FullSimplify[QuantumInterferometer[el, m]["Unitary"] - FourierMatrix[m]] == ConstantArray[0, {m, m}]}], 13, FontFamily -> "Helvetica"], Top]];
reckNull = FullSimplify[{ArcTan[Abs[Cos[α]], Abs[E^(I β) Sin[α]]], Arg[E^(I β) Sin[α] Cos[α]]}, 0 < α < π/2 && -π < β < π];
Labeled[Grid[{{reckPanel[3], reckPanel[4]}}, Alignment -> {Center, Top}, Spacings -> 4],
  Style[Row[{"nulling the splitter B(α, β) gives (θ, φ) = ", TraditionalForm[reckNull], "   for 0 < α < π/2, -π < β < π"}], 13, FontFamily -> "Helvetica"], Top]
```

Every setting of both meshes comes out in closed form, the phases as rational multiples of $\pi$, and the mesh rebuilt from them equals $F_m$ exactly. The nulling formula hands back a symbolic splitter's own angles, so the decomposition is exact algebra, not a numerical fit.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A beam splitter in the ladder-operator algebra, with no photon cutoff

Two modes with annihilation operators $a_1,a_2$ meet on a beam splitter with mixing angle $\theta$ and phase $\varphi$, and $n_1$, $n_2$ photons enter. Without truncating the Fock space, what does the splitter do to the creation operators?

The splitter is $B=e^{G}$ with $G=\theta\left(e^{i\varphi}a_1a_2^\dagger-e^{-i\varphi}a_1^\dagger a_2\right)$ and $[a_j,a_k^\dagger]=\delta_{jk}$. By $e^{G}a_k^\dagger e^{-G}=a_k^\dagger+[G,a_k^\dagger]+\tfrac1{2!}[G,[G,a_k^\dagger]]+\dots$, if $[G,\cdot\,]$ maps $(a_1^\dagger,a_2^\dagger)$ linearly by a matrix $M$, then $a_k^\dagger\to\sum_j W_{jk}a_j^\dagger$ with $W=e^{M}$. The vacuum stays put, so $|n_1,n_2\rangle$ leaves as $(\sum_j W_{j1}x_j)^{n_1}(\sum_j W_{j2}x_j)^{n_2}/\sqrt{n_1!\,n_2!}$ with $x_j$ standing for $a_j^\dagger$, and $x_1^{t_1}x_2^{t_2}$ is $\sqrt{t_1!\,t_2!}\,|t_1,t_2\rangle$.

QF computes the commutators in its bosonic algebra with $\theta,\varphi$ symbolic and compares $e^{M}$ with QuantumInterferometer's transfer matrix (top label); the curves are the probabilities of $t_1$ photons in output 1 read off the polynomial, and the markers are QF's Fock-space `BeamSplitterOperator` at $\varphi=\pi/3$:

```wl
algBSVars = Join[FieldVariables[\[FormalA]], FieldVariables[\[FormalB]]];
algBSGen = θ (Exp[I φ] algBSVars[[1]] ** algBSVars[[4]] - Exp[-I φ] algBSVars[[2]] ** algBSVars[[3]]);
algBSComm = BosonicNormalOrder[Commutator[algBSGen, #], algBSVars, "Scalars" -> {θ, φ}] & /@ algBSVars[[{2, 4}]];
algBSW = FullSimplify[MatrixExp[Transpose[Coefficient[#, algBSVars[[{2, 4}]]] & /@ algBSComm]]];
algBSProb[{n1_, n2_}] := With[{x = Array[\[FormalX], 2]}, KeySort @ Association[(First[#1] -> FullSimplify[Abs[#2]^2 Times @@ (#1!)/(n1! n2!), Element[{θ, φ}, Reals]]) & @@@
     CoefficientRules[Expand[(x . algBSW[[All, 1]])^n1 (x . algBSW[[All, 2]])^n2], x]]];
algBSFock[{n1_, n2_}] := With[{d = n1 + n2 + 1}, Transpose @ Table[With[{pr = ArrayReshape[BeamSplitterOperator[{t, π/3}, d][FockState[{n1, n2}, d]]["ProbabilitiesList"], {d, d}]},
      Table[{t, pr[[j + 1, n1 + n2 - j + 1]]}, {j, 0, n1 + n2}]], {t, 0., N[π/2], N[π/16]}]];
algBSPanel[in_] := Show[Plot[Evaluate[Values[algBSProb[in]]], {θ, 0, π/2}, PlotRange -> All, PlotLegends -> Table[Row[{Subscript["t", 1], " = ", j}], {j, 0, Total[in]}]],
    ListPlot[algBSFock[in], PlotMarkers -> Automatic], PlotRange -> {-0.02, 1.02}, Frame -> True, FrameTicks -> {{Automatic, None}, {Range[0, π/2, π/8], None}},
    FrameLabel -> {"θ", Row[{"P(", Subscript["t", 1], ")"}]}, PlotLabel -> Row[{"input  |", Subscript["n", 1], ", ", Subscript["n", 2], "⟩ = |", in[[1]], ", ", in[[2]], "⟩"}], ImageSize -> 300];
algBSPanels = algBSPanel /@ {{2, 1}, {3, 1}, {3, 2}};
Labeled[Legended[Row[First /@ algBSPanels, Spacer[10]], Last[Last[algBSPanels]]],
  Style[Row[{"[G, ", TraditionalForm[algBSVars[[2]]], "] = ", TraditionalForm[algBSComm[[1]]], ",   [G, ", TraditionalForm[algBSVars[[4]]], "] = ", TraditionalForm[algBSComm[[2]]],
     ",     ", Superscript["e", "M"], " = W of QuantumInterferometer:  ", algBSW === QuantumInterferometer[{"BS"[θ, φ]}]["Unitary"]} /.
     {\[FormalA] -> Subscript["a", 1], \[FormalB] -> Subscript["a", 2]}], 13, FontFamily -> "Helvetica"], Top]
```

The commutator turns each creation operator into a multiple of the other, so the whole series sums to the $2\times2$ transfer matrix and the splitter simply rotates $a_1^\dagger,a_2^\dagger$. The output amplitudes are then coefficients of a polynomial of degree $n_1+n_2$, and their squares need no value of $\varphi$ and sit on the Fock-space markers for every input, although no photon-number cutoff entered the algebra.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Heralding a state from one photon and coherent light

One photon enters port 1 of a beam splitter with mixing angle $\theta$ and phase $\varphi$, and a coherent state $|\beta\rangle$ with complex $\beta$ enters port 2. A photon counter on output 1 registers exactly one photon. What state does that leave in output 2?

The splitter maps $a_1^\dagger\to\cos\theta\,a_1^\dagger+e^{i\varphi}\sin\theta\,a_2^\dagger$ and splits the coherent amplitude into $-e^{-i\varphi}\beta\sin\theta$ in output 1 and $\beta\cos\theta$ in output 2. Projecting output 1 onto $|1\rangle$ keeps two terms, the photon transmitted with the coherent part giving no photon, or the photon reflected with the coherent part giving one, so the unnormalized state of output 2 is

$$
|\psi\rangle=e^{-|\beta|^2\sin^2\theta/2}\left(\cos\theta-\beta\sin^2\theta\,a^\dagger\right)|\beta\cos\theta\rangle ,
$$

with the phase $\varphi$ cancelling between the two terms.

QF computes the amplitudes $\langle 1,k|U|1,\beta\rangle$ for $k\le14$ with $\theta,\varphi,\beta$ symbolic and compares them with $\langle k|\psi\rangle$ (top label); the panels show $\pi$ times QF's Wigner function $W$ of the normalized state built from them, at $\arg\beta=\pi/4$, with quadratures $x,p$ defined by $a=(x+ip)/\sqrt2$; the dashed line encloses $\pi W<-0.01$:

```wl
heraldAmps = Table[QuantumInterferometer[{"BS"[θ, φ]}]["Amplitude", {1, "Coherent"[β]} -> {1, k}], {k, 0, 14}];
heraldModel = Table[Exp[-Abs[β]^2/2] (Cos[θ] (β Cos[θ])^k - k β Sin[θ]^2 (β Cos[θ])^(k - 1))/Sqrt[k!], {k, 0, 14}];
heraldGrid = Range[-4.5, 4.5, 0.075];
heraldWigner[t_, b_] := With[{w = WignerFunction[QuantumState[N[heraldAmps /. {θ -> t, β -> b Exp[I π/4]}], 15]["Normalized"], {x, p}]},
    Re[w /. {x -> ConstantArray[heraldGrid, Length[heraldGrid]], p -> Transpose[ConstantArray[heraldGrid, Length[heraldGrid]]]}]];
heraldColors = Blend["WL14DefaultDensityGradient", (1 + #)/2] &;
heraldPanel[t_, b_] := With[{data = π heraldWigner[t, b]}, Show[
    ListDensityPlot[data, DataRange -> {{-4.5, 4.5}, {-4.5, 4.5}}, PlotRange -> {-1, 1}, ColorFunction -> heraldColors, ColorFunctionScaling -> False],
    ListContourPlot[data, DataRange -> {{-4.5, 4.5}, {-4.5, 4.5}}, Contours -> {-0.01}, ContourShading -> None, ContourStyle -> Dashed],
    Frame -> True, FrameLabel -> {"x", "p"}, PlotLabel -> Row[{"θ = ", t, ",   |β| = ", b}], ImageSize -> 230]];
heraldPanels = Outer[heraldPanel, {π/8, π/4, 3 π/8}, {0.5, 1, 2}];
Labeled[Legended[Grid[heraldPanels], BarLegend[{heraldColors, {-1, 1}}, LegendLabel -> "πW"]],
  Style[Row[{"QF:  ", TraditionalForm[FullSimplify[heraldAmps[[;; 2]], Element[θ, Reals]]], " for k = 0, 1;    all k ≤ 14 equal ⟨k|ψ⟩:  ",
     FullSimplify[heraldAmps - heraldModel, Element[θ, Reals]] === ConstantArray[0, 15]}], 13, FontFamily -> "Helvetica"], Top]
```

With weak coherent light and a small $\theta$ the click almost always comes from the transmitted photon, and output 2 keeps a positive bump close to the coherent state $|\beta\cos\theta\rangle$. A larger $\theta$ or a brighter $\beta$ gives the reflected photon, the $a^\dagger$ term, more weight, and a negative region opens beside a crescent of positive $W$, which no classical field can produce. At $\theta=3\pi/8$ every $|\beta|$ shown leaves a negative center.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A Kerr medium in a Mach-Zehnder interferometer

A beam splitter $B(\theta,\varphi)$ opens a Mach-Zehnder interferometer, arm 1 carries a linear phase $\alpha$ and a Kerr medium with nonlinear phase $\chi$, and the inverse splitter $B^\dagger$ closes it. Either one photon enters port 1, or one photon enters each port. What does the Kerr phase change that the linear phase cannot?

The network is

$$
U=B(\theta,\varphi)^\dagger\,e^{i\chi\,n_1(n_1-1)/2}\,e^{i\alpha\,n_1}\,B(\theta,\varphi),\qquad B(\theta,\varphi)=e^{\theta\left(e^{i\varphi}a_1a_2^\dagger-e^{-i\varphi}a_1^\dagger a_2\right)},
$$

with $n_1=a_1^\dagger a_1$; the Kerr factor is 1 on zero or one photon in the arm. The two linear pieces are QuantumInterferometer networks, and their `"CircuitOperator"` turns them into Fock-space circuits on three levels per mode, enough for two photons, with $\theta,\varphi,\alpha$ still symbolic; the Kerr medium is a diagonal QF operator between them.

QF composes the circuit with $\theta,\varphi,\alpha,\chi$ symbolic and plots the probability that the input pattern comes back; the markers are QuantumInterferometer's permanent route for the linear network ($\chi=0$) at $\varphi=\pi/3$:

```wl
kerrMZIn = QuantumInterferometer[{"BS"[θ, φ], "PS"[α] -> 1}];
kerrMZOut = QuantumInterferometer[{"BS"[θ, φ]}]["Dagger"];
kerrMZKerr = QuantumOperator[DiagonalMatrix[SparseArray[Exp[I χ Range[0, 2] (Range[0, 2] - 1)/2]]], {1}, 3];
kerrMZCircuit = QuantumCircuitOperator[{kerrMZIn["CircuitOperator", 3], kerrMZKerr, kerrMZOut["CircuitOperator", 3]}];
kerrMZReturn[in_] := FullSimplify[Abs[(SuperDagger[FockState[in, 3]] @ kerrMZCircuit[FockState[in, 3]])["Scalar"]]^2, Element[{θ, φ, α, χ}, Reals]];
{kerrMZ1, kerrMZ2} = kerrMZReturn /@ {{1, 0}, {1, 1}};
kerrMZLinear = kerrMZOut[kerrMZIn]["Probability", # -> #] & /@ {{1, 0}, {1, 1}};
kerrMZPanel[t_] := Show[
    Plot[Evaluate[Prepend[Table[kerrMZ2 /. {θ -> t, χ -> c}, {c, {0, π/2, π}}], kerrMZ1 /. θ -> t]], {α, 0, 2 π}, PlotRange -> All, PlotStyle -> {Dashed, Automatic, Automatic, Automatic},
     PlotLegends -> {"one photon |1, 0⟩", "pair |1, 1⟩,  χ = 0", "pair |1, 1⟩,  χ = π/2", "pair |1, 1⟩,  χ = π"}],
    ListPlot[Table[{a, # /. {θ -> t, φ -> π/3, α -> a}}, {a, 0, 2 π, π/8}] & /@ kerrMZLinear, PlotMarkers -> Automatic],
    PlotRange -> {-0.02, 1.02}, Frame -> True, FrameTicks -> {{Automatic, None}, {Range[0, 2 π, π/2], None}},
    FrameLabel -> {"α", "return probability"}, PlotLabel -> Row[{"θ = ", t}], ImageSize -> 320];
kerrMZPanels = kerrMZPanel /@ {π/4, π/8};
Labeled[Legended[Row[First /@ kerrMZPanels, Spacer[10]], Last[First[kerrMZPanels]]],
  Style[Column[{Row[{"P(|1, 0⟩ → |1, 0⟩) = ", TraditionalForm[kerrMZ1]}], Row[{"P(|1, 1⟩ → |1, 1⟩) = ", TraditionalForm[kerrMZ2]}]}], 13, FontFamily -> "Helvetica"], Top]
```

Neither return probability contains $\varphi$, and the one-photon fringe contains no $\chi$ either: a single photon in the arm picks up no Kerr phase, so the Kerr medium acts on photon pairs alone. On the balanced splitter the pair swings fully between returning and never returning, and $\chi$ only shifts that swing in $\alpha$. At $\theta=\pi/8$ the one-photon fringe no longer empties, and $\chi$ changes the shape of the pair fringe as well as its position.
