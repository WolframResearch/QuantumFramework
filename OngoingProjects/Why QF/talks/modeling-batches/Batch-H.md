```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A beam splitter with symbolic photon numbers

$n_1$ and $n_2$ photons enter the two ports of a beam splitter with mixing angle $\theta$ and phase $\varphi$. What is the probability of finding $t_1$ photons in the first output, for every photon number at once?

The splitter rotates the creation operators,

$$
a_1^\dagger\to\cos\theta\,a_1^\dagger+e^{i\varphi}\sin\theta\,a_2^\dagger,\qquad a_2^\dagger\to-e^{-i\varphi}\sin\theta\,a_1^\dagger+\cos\theta\,a_2^\dagger,
$$

so it conserves the photon number, and the second output holds $n_1+n_2-t_1$ photons.

QF returns the probability with $n_1,n_2,t_1,\theta,\varphi$ all symbolic (top label); the panels substitute $n$ photons in each input and show $t_1\le n$; for equal inputs $|d^{\,n}_{0,m}|=|d^{\,n}_{0,-m}|$ with $m=t_1-n$, so $t_1$ and $2n-t_1$ are equally likely:

```wl
splitterProb = QuantumInterferometer[{"BS"[θ, φ]}]["Probability", {n1, n2} -> {t1, n1 + n2 - t1}];
twinPanel[m_] := Plot[Evaluate[FullSimplify[Table[splitterProb /. {n1 -> m, n2 -> m, t1 -> j}, {j, 0, m}], Element[{θ, φ}, Reals]]],
    {θ, 0, π/2}, PlotRange -> {0, 1}, PlotStyle -> Table[If[OddQ[j], Dashed, Automatic], {j, 0, m}],
    PlotLegends -> Table[Row[{"t1 = ", j}], {j, 0, m}], GridLines -> {{π/4}, None}, Frame -> True, FrameTicks -> {{Automatic, None}, {Range[0, π/2, π/8], None}},
    FrameLabel -> {"θ", "P(t1)"}, PlotLabel -> Row[{"n1 = n2 = ", m}], ImageSize -> 300];
twinPanels = twinPanel /@ {1, 2, 3};
Labeled[Legended[Row[First /@ twinPanels], Last[Last[twinPanels]]], TraditionalForm[splitterProb], Top]
```

The probability is the squared Wigner D-function of a spin $j=(n_1+n_2)/2$ rotated by $2\theta$: two modes holding $N$ photons carry the spin-$N/2$ representation of SU(2), and the splitter rotates it. The phase $\varphi$ enters the D-function only as a phase factor, so it drops out of the probability and the panels need no value for it. At the balanced splitter (vertical line) every odd $t_1$ (dashed) vanishes in each panel, so equal inputs leave with an even number of photons in each output; $n=1$ is the Hong-Ou-Mandel zero.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Photon loss, and a Hong-Ou-Mandel zero that survives it

One photon enters each input of a beam splitter with angle $\theta$ and phase $\varphi$, and an attenuator of transmissivity $\eta$ sits on each output. If both outputs lose light at the same $\eta$, does the two-photon coincidence zero survive, and how does $\eta$ enter?

Pure loss mixes the mode with a vacuum mode $b$ and discards $b$: $a\to\sqrt\eta\,a+\sqrt{1-\eta}\,b$. Each photon then survives independently, so $n$ photons leave $k$ survivors with the binomial probability $P(k\mid n)=\binom nk\eta^k(1-\eta)^{n-k}$. The splitter $a_1^\dagger\to\cos\theta\,a_1^\dagger+e^{i\varphi}\sin\theta\,a_2^\dagger$, $a_2^\dagger\to-e^{-i\varphi}\sin\theta\,a_1^\dagger+\cos\theta\,a_2^\dagger$ acts first, then the loss on each output.

QF gives $P(k\mid n)$ with $n,k,\eta$ symbolic (left label), shown at $n=4$, and the coincidence probability with $\theta,\varphi,\eta$ symbolic (right label):

```wl
lossLaw = QuantumInterferometer[{"Loss"[η]}]["Probability", {n} -> {k}];
lossyCoincidence = FullSimplify[QuantumInterferometer[{"BS"[θ, φ], "Loss"[η] -> 1, "Loss"[η] -> 2}]["Probability", {1, 1} -> {1, 1}]];
Row[{
  BarChart[Table[lossLaw /. {n -> 4, k -> j, η -> e}, {j, 0, 4}, {e, {0.3, 0.6, 0.9}}],
    PlotLabels -> <|"Groups" -> Table[Row[{"k = ", j}], {j, 0, 4}]|>, PlotLegends -> <|"Elements" -> {"η = 0.3", "η = 0.6", "η = 0.9"}|>,
    Frame -> True, FrameLabel -> {None, "P(k | n = 4)"}, PlotLabel -> TraditionalForm[lossLaw], ImageSize -> 320],
  Plot[Evaluate[Table[lossyCoincidence /. η -> e, {e, {1, 0.7, 0.4}}]], {θ, 0, π/2}, PlotRange -> {0, 1},
    PlotLegends -> Table[Row[{"η = ", e}], {e, {1, 0.7, 0.4}}], GridLines -> {{π/4}, None}, Frame -> True, FrameTicks -> {{Automatic, None}, {Range[0, π/2, π/8], None}},
    FrameLabel -> {"θ", "P(1, 1)"}, PlotLabel -> TraditionalForm[lossyCoincidence], ImageSize -> 320]}]
```

Loss only scales the coincidence curve by $\eta^2$, the chance that both photons survive. Its zero at the balanced splitter stays in place for every $\eta$: a coincidence needs both photons to survive, and a surviving pair carries the lossless amplitude.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Two squeezed vacua on a balanced beam splitter

Each of two modes is squeezed from the vacuum, with complex parameters $\xi_1$ and $\xi_2$, and the modes meet on a balanced beam splitter. When does one photon appear in each output?

Each squeezer is $S(\xi)=\exp\!\big[(\xi^*a^2-\xi\,a^{\dagger2})/2\big]$; the splitter then maps $a_1^\dagger\to(a_1^\dagger+a_2^\dagger)/\sqrt2$, $a_2^\dagger\to(a_2^\dagger-a_1^\dagger)/\sqrt2$.

QF computes the coincidence amplitude $\langle1,1|U|0,0\rangle$ with $\xi_1,\xi_2$ symbolic, plotted for $\xi_1=r$, $\xi_2=q\,re^{i\phi}$ with $q=|\xi_2|/|\xi_1|$ (first two panels), and the pair probabilities $P(m,m)$ for opposite squeezers $\pm r$ with $r$ symbolic (third panel, closed forms in the legend):

```wl
pairAmplitude = QuantumInterferometer[{"S"[ξ1] -> 1, "S"[ξ2] -> 2, "BS"[π/4]}]["Amplitude", {0, 0} -> {1, 1}];
oppositeSqueezers = QuantumInterferometer[{"S"[r] -> 1, "S"[-r] -> 2, "BS"[π/4]}];
squeezerPairLaw = Table[FullSimplify[oppositeSqueezers["Probability", {0, 0} -> {m, m}, Assumptions -> r > 0], r > 0], {m, 0, 3}];
coincidencePanel[q_] := Plot[Evaluate[Table[Abs[pairAmplitude /. {ξ1 -> r, ξ2 -> q r Exp[I p]}]^2, {p, {0, π/2, π}}]], {r, 0, 2},
    PlotRange -> {-0.005, 0.26}, PlotLegends -> {"ϕ = 0", "ϕ = π/2", "ϕ = π"}, Frame -> True, FrameLabel -> {"r", "P(1, 1)"},
    PlotLabel -> Row[{"|ξ2|/|ξ1| = ", q}], ImageSize -> 300];
Grid[{{Legended[Grid[{First /@ coincidencePanel /@ {1/2, 1}}], Last[coincidencePanel[1]]],
  Plot[Evaluate[squeezerPairLaw], {r, 0, 2}, PlotRange -> {0, 1}, PlotLegends -> Table[Row[{"m = ", m - 1, ":  ", TraditionalForm[squeezerPairLaw[[m]]]}], {m, 4}],
    Frame -> True, FrameLabel -> {"r", "P(m, m)"}, PlotLabel -> "ξ1 = r,  ξ2 = -r", ImageSize -> 300]}}]
```

The coincidence vanishes only for identical squeezers ($q=1$, $\phi=0$): the pair generator $a_1^{\dagger2}+a_2^{\dagger2}$ is unchanged by the splitter, so the output is again two independent squeezed vacua, each holding only even photon numbers. Any difference in strength or phase lets coincidences through, most for opposite squeezers, whose outputs carry $m$ photons each with probability $\operatorname{sech}^2r\tanh^{2m}r$, the photon statistics of a two-mode squeezed vacuum.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The three-mode Fourier interferometer

One photon enters each port of a three-mode interferometer whose transfer matrix is the discrete Fourier matrix. Which output patterns occur?

The network maps $a_k^\dagger\to\sum_j W_{jk}\,a_j^\dagger$ with $W_{jk}=\omega^{(j-1)(k-1)}/\sqrt3$ and $\omega=e^{2\pi i/3}$; each output amplitude is a permanent of a $3\times3$ matrix built from $W$.

QF computes the exact probability of every output pattern $\{t_1,t_2,t_3\}$ from $|1,1,1\rangle$:

```wl
tritterProbabilities = QuantumInterferometer[FourierMatrix[3]]["Probabilities", {1, 1, 1}];
BarChart[Values[tritterProbabilities], PlotLabels -> (Row[#["Name"], ","] & /@ Keys[tritterProbabilities]),
  LabelingFunction -> Above, PlotRange -> {0, 0.4}, Frame -> True, FrameLabel -> {"output pattern", "probability"}, ImageSize -> 520]
```

Every pattern with two photons in one port and one in another is exactly zero, while $\{1,1,1\}$ and the three fully bunched patterns share all the probability. Relabeling the input ports cyclically leaves $|1,1,1\rangle$ unchanged but multiplies the amplitude of pattern $\{t_1,t_2,t_3\}$ by $\omega^{\sum_j (j-1)t_j}$, so the amplitude must vanish unless that sum is a multiple of 3, which no 2+1 pattern satisfies.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The nonlinear sign gate of Knill, Laflamme and Milburn

A signal mode holds $c_0|0\rangle+c_1|1\rangle+c_2|2\rangle$; two auxiliary modes hold one photon and the vacuum. A three-mode network mixes them, and a run is kept only when the auxiliary detectors count $\{1,0\}$. A phase shifter cannot flip the sign of $c_2$ alone, since it multiplies $|n\rangle$ by $e^{in\alpha}$. What does the kept signal look like?

The network maps $a_k^\dagger\to\sum_j W_{jk}\,a_j^\dagger$ with the real transfer matrix (Knill, Laflamme and Milburn, Nature 409, 46 (2001))

$$
W=\begin{pmatrix}1-\sqrt2 & 2^{-1/4} & \sqrt{3/\sqrt2-2}\\ 2^{-1/4} & 1/2 & (1-\sqrt2)/2\\ \sqrt{3/\sqrt2-2} & (1-\sqrt2)/2 & \sqrt2-1/2\end{pmatrix},
$$

so the kept signal is $\sum_m c_m\langle m,1,0|U|m,1,0\rangle\,|m\rangle$.

QF computes these amplitudes for $m=0,\dots,5$ (bars); the first three give the kept signal for symbolic $c_0,c_1,c_2$ (top label):

```wl
nsGate = QuantumInterferometer[{{1 - Sqrt[2], 2^(-1/4), Sqrt[3/Sqrt[2] - 2]}, {2^(-1/4), 1/2, (1 - Sqrt[2])/2},
    {Sqrt[3/Sqrt[2] - 2], (1 - Sqrt[2])/2, Sqrt[2] - 1/2}}];
nsAmps = RootReduce[Table[nsGate["Amplitude", {m, 1, 0} -> {m, 1, 0}], {m, 0, 5}]];
keptSignal = {c0, c1, c2} Take[nsAmps, 3];
Labeled[BarChart[nsAmps, PlotLabels -> Range[0, 5],
    LabelingFunction -> Above, Frame -> True, FrameLabel -> {"signal photons m", "⟨m,1,0|U|m,1,0⟩"}, ImageSize -> 480],
  Row[{"kept signal: ", TraditionalForm[keptSignal]}], Top]
```

For zero, one and two signal photons the amplitudes have the same size and only the two-photon one is negative, so every input comes out with $c_2$ flipped and the same overall factor: the gate works with a probability that does not depend on the input. From three signal photons on the amplitudes differ, so the gate works only on inputs with at most two signal photons.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Bosons, distinguishable particles and fermions at one beam splitter

Two particles, one in each input of a beam splitter with angle $\theta$ and phase $\varphi$, can both be transmitted or both be reflected. How likely is it that they leave through different ports?

With transfer matrix $W=\begin{pmatrix}\cos\theta & -e^{-i\varphi}\sin\theta\\ e^{i\varphi}\sin\theta & \cos\theta\end{pmatrix}$, bosons add the two histories with a plus sign, fermions with a minus sign, and distinguishable particles add their probabilities:

$$
P_{\rm B}=|\operatorname{Per}W|^2,\qquad P_{\rm F}=|\det W|^2,\qquad P_{\rm D}=\operatorname{Per}\big(|W_{jk}|^2\big).
$$

QF gives $W$ with $\theta,\varphi$ symbolic, and the three coincidence probabilities follow in closed form (legend):

```wl
statisticsW = QuantumInterferometer[{"BS"[θ, φ]}]["Unitary"];
statisticsCoincidences = FullSimplify[{Abs[Permanent[statisticsW]]^2, Permanent[Abs[statisticsW]^2], Abs[Det[statisticsW]]^2}, Element[{θ, φ}, Reals]];
Plot[Evaluate[statisticsCoincidences], {θ, 0, π/2}, PlotRange -> {0, 1.05}, GridLines -> {{π/4}, None},
  PlotLegends -> MapThread[Row[{#1, ":  ", TraditionalForm[#2]}] &, {{"bosons", "distinguishable", "fermions"}, statisticsCoincidences}],
  Frame -> True, FrameTicks -> {{Automatic, None}, {Range[0, π/2, π/8], None}}, FrameLabel -> {"θ", "P(1, 1)"}, ImageSize -> 420]
```

The phase $\varphi$ drops out of all three. Bosons never separate at the balanced splitter, fermions separate at every $\theta$, and distinguishable particles lie between them. Bosons and fermions differ only in the sign of the exchange term; distinguishable particles have none, so their curve is the average of the other two.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Hong-Ou-Mandel with partially distinguishable photons

Two photons meet at a beam splitter with angle $\theta$ and phase $\varphi$. The first is horizontally polarized, the second is polarized at angle $\chi$, $\cos\chi\,|H\rangle+\sin\chi\,|V\rangle$, and the detectors do not resolve polarization. How does the coincidence probability depend on $\theta$ and $\chi$?

Treat the polarizations as extra modes $1_H,2_H,1_V,2_V$, with the same splitter $a_1^\dagger\to\cos\theta\,a_1^\dagger+e^{i\varphi}\sin\theta\,a_2^\dagger$, $a_2^\dagger\to-e^{-i\varphi}\sin\theta\,a_1^\dagger+\cos\theta\,a_2^\dagger$ on the $H$ pair and on the $V$ pair. The input is $\cos\chi\,|1,1,0,0\rangle+\sin\chi\,|1,0,0,1\rangle$, and the coincidence sums over the four outputs $T$ with one photon in each path:

$$
P_{\rm coinc}=\sum_T\big|\cos\chi\,A(1100\to T)+\sin\chi\,A(1001\to T)\big|^2 .
$$

QF computes the amplitudes $A$ with $\theta,\varphi,\chi$ symbolic; the closed form is the plot label:

```wl
polarizedHOM = QuantumInterferometer[{"BS"[θ, φ] -> {1, 2}, "BS"[θ, φ] -> {3, 4}}, 4];
partialCoincidence = FullSimplify[Collect[Simplify[ComplexExpand[Total[Abs[Cos[χ] polarizedHOM["Amplitude", {1, 1, 0, 0} -> #] + Sin[χ] polarizedHOM["Amplitude", {1, 0, 0, 1} -> #]]^2 & /@
    {{1, 1, 0, 0}, {1, 0, 0, 1}, {0, 1, 1, 0}, {0, 0, 1, 1}}]] /. Sin[χ]^2 -> 1 - Cos[χ]^2], Cos[χ]^2, FullSimplify]];
Plot[Evaluate[Table[partialCoincidence /. χ -> c, {c, {0, π/8, π/4, 3 π/8, π/2}}]], {θ, 0, π/2}, PlotRange -> {0, 1},
  PlotLegends -> {"χ = 0", "χ = π/8", "χ = π/4", "χ = 3π/8", "χ = π/2"}, GridLines -> {{π/4}, None},
  Frame -> True, FrameTicks -> {{Automatic, None}, {Range[0, π/2, π/8], None}}, FrameLabel -> {"θ", "P(coincidence)"}, PlotLabel -> TraditionalForm[partialCoincidence], ImageSize -> 420]
```

Identical photons ($\chi=0$) give the Hong-Ou-Mandel curve $\cos^22\theta$, and orthogonal ones ($\chi=\pi/2$) behave as distinguishable particles. The phase $\varphi$ drops out, and the label is the distinguishable-particle value $\cos^4\theta+\sin^4\theta$ minus an interference term weighted by the squared overlap $\cos^2\chi$ of the two polarizations, so at the balanced splitter the coincidence is $\tfrac12\sin^2\chi$: the depth of the dip measures how indistinguishable the photons are.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### When n photons per port return to their inputs

$n$ photons enter each port of a beam splitter with reflectivity $R=\sin^2\theta$ and phase $\varphi$. For which reflectivities does the probability that they leave as they came, $n$ in each output, vanish?

The splitter maps $a_1^\dagger\to\cos\theta\,a_1^\dagger+e^{i\varphi}\sin\theta\,a_2^\dagger$, $a_2^\dagger\to-e^{-i\varphi}\sin\theta\,a_1^\dagger+\cos\theta\,a_2^\dagger$.

QF returns the amplitude $\langle n,n|U|n,n\rangle$ with $n,\theta,\varphi$ symbolic (top label); the curves substitute $n=1,\dots,4$ and $\theta=\arcsin\sqrt R$:

```wl
returnAmplitude = QuantumInterferometer[{"BS"[θ, φ]}]["Amplitude", {n, n} -> {n, n}];
Labeled[Plot[Evaluate[Table[Abs[returnAmplitude /. {n -> j, θ -> ArcSin[Sqrt[R]]}]^2, {j, 4}]], {R, 0, 1}, PlotRange -> {0, 1},
    PlotLegends -> Table[Row[{"n = ", j}], {j, 4}], Frame -> True, FrameLabel -> {"R", "P(n, n → n, n)"}, ImageSize -> 420],
  TraditionalForm[returnAmplitude], Top]
```

The return amplitude is the rotation matrix element $d^{\,n}_{00}(2\theta)$, which is the Legendre polynomial $P_n(\cos2\theta)=P_n(1-2R)$, so $n$ photons per port never return at the $n$ zeros of $P_n$. For even $n$ the balanced splitter does not suppress the return, and for $n=2$ the two zeros sit at unbalanced reflectivities symmetric about $R=1/2$; every odd $n$ keeps a zero at the balanced splitter.
