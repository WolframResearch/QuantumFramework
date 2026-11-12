```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### GHZ and W: two kinds of three-qubit entanglement

Three qubits share $|\mathrm{GHZ}\rangle=(|000\rangle+|111\rangle)/\sqrt2$ or $|W\rangle=(|001\rangle+|010\rangle+|100\rangle)/\sqrt3$. Each party may apply any invertible local map $A$, $B$, $C$, even one that succeeds only with some probability. Can such maps turn GHZ into W?

The maps act as $|\psi\rangle\mapsto(A\otimes B\otimes C)|\psi\rangle$, followed by normalization. Cut the amplitudes $\psi_{ijk}$ into the two $2\times2$ slices $M_i=(\psi_{ijk})_{jk}$; the Cayley hyperdeterminant $\operatorname{Det}\psi$ is the discriminant of the quadratic $\det(M_0+xM_1)$, and the three-tangle of a normalized state is $\tau_3=4|\operatorname{Det}\psi|$. The family

$$
|\psi_\varepsilon\rangle\propto A_\varepsilon^{\otimes3}|\mathrm{GHZ}\rangle\propto(|0\rangle+\varepsilon|1\rangle)^{\otimes3}-|000\rangle,\qquad A_\varepsilon=\begin{pmatrix}1&-1\\ \varepsilon&0\end{pmatrix},\quad\det A_\varepsilon=\varepsilon,
$$

is GHZ moved by an invertible map for every $\varepsilon>0$.

The cell transforms both states with all twelve entries of $A,B,C$ symbolic, then gives $\tau_3$ and the infidelity $1-|\langle W|\psi_\varepsilon\rangle|^2$ in closed form:

```wl
hyperdetSLOCC[v_] := With[{c = CoefficientList[Det[#[[1]] + \[FormalX] #[[2]]], \[FormalX], 3]}, c[[2]]^2 - 4 c[[1]] c[[3]]] &[ArrayReshape[Normal[v], {2, 2, 2}]];
{ghzSLOCC, wSLOCC} = QuantumState[#[3]] & /@ {"GHZ", "W"};
mapSLOCC = QuantumTensorProduct @@ (QuantumOperator[Array[#, {2, 2}]] & /@ {ma, mb, mc});
detsSLOCC = Times @@ (Det[Array[#, {2, 2}]] & /@ {ma, mb, mc});
{ghzRuleQ, wRuleQ} = {Expand[hyperdetSLOCC[mapSLOCC[ghzSLOCC]["StateVector"]] - detsSLOCC^2/4] === 0, Expand[hyperdetSLOCC[mapSLOCC[wSLOCC]["StateVector"]]] === 0};
vecSLOCC = Normal[(QuantumTensorProduct @@ ConstantArray[QuantumOperator[{{1, -1}, {\[Epsilon], 0}}], 3])[ghzSLOCC]["StateVector"]];
{tau3SLOCC, infidSLOCC} = Simplify[{4 Abs[hyperdetSLOCC[vecSLOCC]]/(vecSLOCC . vecSLOCC)^2, 1 - (Normal[wSLOCC["StateVector"]] . vecSLOCC)^2/(vecSLOCC . vecSLOCC)}, \[Epsilon] > 0];
LogLogPlot[{tau3SLOCC, infidSLOCC}, {\[Epsilon], 10^-3, 1}, PlotRange -> All, PlotStyle -> {Automatic, Dashed},
 PlotLegends -> {Subscript["\[Tau]", 3], "1 - |\[LeftAngleBracket]W|\[Psi]\[RightAngleBracket]|²"}, Frame -> True, FrameLabel -> {"\[Epsilon]", None}, ImageSize -> 420,
 PlotLabel -> Style[Column[{Row[{"Det((A\[CircleTimes]B\[CircleTimes]C)GHZ) = (det A det B det C)²/4: ", ghzRuleQ}], Row[{"Det((A\[CircleTimes]B\[CircleTimes]C)W) = 0: ", wRuleQ}],
     Row[{Subscript["\[Tau]", 3], " = ", TraditionalForm[tau3SLOCC], ",   1 - |\[LeftAngleBracket]W|\[Psi]\[RightAngleBracket]|² = ", TraditionalForm[infidSLOCC]}]}], 13, FontFamily -> "Helvetica"]]
```

An invertible local map multiplies the hyperdeterminant by $(\det A\det B\det C)^2\neq0$, so $\tau_3$ of W stays zero and $\tau_3$ of GHZ stays positive, and no such map links the two in either direction. The family $\psi_\varepsilon$ is GHZ-class at every $\varepsilon>0$, yet its $\tau_3$ and its infidelity with W both vanish as $\varepsilon^2$: W sits on the boundary of the GHZ class, reachable only as a limit.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Monogamy of entanglement: the CKW inequality

Three qubits $A$, $B$, $C$ share $|\psi_{p,\phi}\rangle=\sqrt p\,|\mathrm{GHZ}\rangle+e^{i\phi}\sqrt{1-p}\,|W\rangle$, a superposition of the two orthogonal states $(|000\rangle+|111\rangle)/\sqrt2$ and $(|001\rangle+|010\rangle+|100\rangle)/\sqrt3$. How much of qubit $A$'s entanglement with the pair $BC$ can it hold with $B$ and with $C$ separately?

$A$'s entanglement with $BC$ is the tangle $\tau_{A|BC}=C_{A|BC}^2=4\det\rho_A$. What it shares with $B$ alone is the Wootters concurrence $C_{AB}$ of the two-qubit state $\rho_{AB}=\operatorname{Tr}_C|\psi\rangle\langle\psi|$, and likewise $C_{AC}$. The Coffman-Kundu-Wootters inequality says

$$
\tau_{A|BC}\;\ge\;C_{AB}^2+C_{AC}^2,\qquad \tau_3=\tau_{A|BC}-C_{AB}^2-C_{AC}^2\ge0 ,
$$

and the residual $\tau_3$ is the three-tangle, the part of the entanglement that no pair holds.

QF computes all three concurrences with `QuantumEntanglementMonotone`, from the pure state and from its two-qubit reductions, on a grid of $p$:

```wl
{ghzCKW, wCKW} = Normal[QuantumState[#[3]]["StateVector"]] & /@ {"GHZ", "W"};
ckwTerms[p_, \[Phi]_] := With[{s = QuantumState[Sqrt[p] ghzCKW + Exp[I \[Phi]] Sqrt[1 - p] wCKW, {2, 2, 2}]},
   Re @ With[{tangle = QuantumEntanglementMonotone[s, {{1}, {2, 3}}, "Concurrence"]^2, pairs = Total[QuantumEntanglementMonotone[QuantumPartialTrace[s, {#}], "Concurrence"]^2 & /@ {3, 2}]},
     {tangle, pairs, tangle - pairs}]];
panelCKW[\[Phi]_] := With[{data = Table[{p, #} & /@ ckwTerms[p, \[Phi]], {p, Union[Subdivide[0., 1., 50], Subdivide[0.6, 0.66, 60]]}]},
   ListLinePlot[Transpose[data], PlotRange -> {-0.02, 1.02}, PlotStyle -> {Automatic, Dashed, Dotted}, Frame -> True, FrameLabel -> {"p", None}, ImageSize -> 300,
    PlotLegends -> {Subscript["\[Tau]", "A|BC"], Row[{Subsuperscript["C", "AB", 2], " + ", Subsuperscript["C", "AC", 2]}], Row[{Subscript["\[Tau]", 3], " (residual)"}]},
    PlotLabel -> Style[Row[{"\[Phi] = ", \[Phi]}], 13, FontFamily -> "Helvetica"]]];
panelsCKW = panelCKW /@ {0, Pi/2, Pi};
Legended[Row[First /@ panelsCKW], Last[First[panelsCKW]]]
```

At every weight and at all three phases the pairwise sum never exceeds the tangle: $A$ cannot hold more with $B$ and $C$ one at a time than it holds with the pair. At $p=0$ the W state spends all of it on pairs, and at $p=1$ GHZ spends none, so all of it is three-way. The phase moves the split: at $\phi=\pi$ the residual touches zero at one intermediate weight, a superposition that, like W, has no three-way entanglement.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Which pure states LOCC can convert: Nielsen's majorization

Alice and Bob share a pure state and may use local operations and classical communication (LOCC). Write each state in Schmidt form, $|\psi\rangle=\sum_k\sqrt{\lambda_k}\,|k\rangle_A|k\rangle_B$, with weights $\lambda_1\ge\lambda_2\ge\dots$. Two qutrit pairs have weights $\lambda^\psi=(\tfrac12,\tfrac3{10},\tfrac15)$ and $\lambda^\phi=(\tfrac35,\tfrac15,\tfrac15)$; two ququart pairs, again called $\psi$ and $\phi$, have $\lambda^\psi=(\tfrac25,\tfrac25,\tfrac1{10},\tfrac1{10})$ and $\lambda^\phi=(\tfrac12,\tfrac14,\tfrac14,0)$; a third pair $\chi$, a qubit pair with weights $(\tfrac35,\tfrac25)$, can be borrowed and must be returned. Which conversions succeed with certainty?

Nielsen's theorem: $|\psi\rangle\to|\phi\rangle$ by LOCC with certainty if and only if $\lambda^\psi$ is majorized by $\lambda^\phi$,

$$
\lambda^\psi\prec\lambda^\phi\iff\sum_{i\le k}\lambda^\psi_i\le\sum_{i\le k}\lambda^\phi_i\quad\text{for every }k .
$$

With $\chi$ borrowed, the question is $|\psi\rangle|\chi\rangle\to|\phi\rangle|\chi\rangle$, where Alice holds both first halves.

QF reads the weights as the eigenvalues of Alice's reduced state, and the cell plots the partial sums $\sum_{i\le k}\lambda_i$ against $k$; $\psi\to\phi$ when the curve of $\phi$ is never below that of $\psi$:

```wl
schmidtNielsen[w_] := QuantumState[Flatten[DiagonalMatrix[Sqrt[w]]], {Length[w], Length[w]}];
aliceNielsen[s_] := ReverseSort[QuantumPartialTrace[s, Range[2, s["Qudits"], 2]]["Eigenvalues"]];
majorizedNielsen[x_, y_] := And @@ Thread[Accumulate[x] <= Accumulate[y]];
catalystNielsen = schmidtNielsen[{3/5, 2/5}];
pairsNielsen = {schmidtNielsen /@ {{1/2, 3/10, 1/5}, {3/5, 1/5, 1/5}}, schmidtNielsen /@ {{2/5, 2/5, 1/10, 1/10}, {1/2, 1/4, 1/4, 0}},
   QuantumTensorProduct[schmidtNielsen[#], catalystNielsen] & /@ {{2/5, 2/5, 1/10, 1/10}, {1/2, 1/4, 1/4, 0}}};
panelNielsen[{s1_, s2_}, title_] := With[{w1 = aliceNielsen[s1], w2 = aliceNielsen[s2]},
   ListLinePlot[{Accumulate[Prepend[w1, 0]], Accumulate[Prepend[w2, 0]]}, DataRange -> {0, Length[w1]}, PlotMarkers -> {Automatic, 9}, PlotStyle -> {Automatic, Dashed},
    PlotRange -> {0, 1.02}, Frame -> True, FrameTicks -> {{Automatic, None}, {Range[0, Length[w1]], None}}, FrameLabel -> {"k", Row[{"\[Sum]", Subscript["\[Lambda]", "i \[LessEqual] k"]}]}, ImageSize -> 300, PlotLegends -> {"\[Psi]", "\[Phi]"},
    PlotLabel -> Style[Column[{title, Row[{"\[Psi] \[RightArrow] \[Phi]: ", majorizedNielsen[w1, w2], ",   \[Phi] \[RightArrow] \[Psi]: ", majorizedNielsen[w2, w1]}]}, Alignment -> Center], 13, FontFamily -> "Helvetica"]]];
panelsNielsen = MapThread[panelNielsen, {pairsNielsen, {"qutrit pairs", "ququart pairs", Row[{"ququart pairs with ", "\[Chi]"}]}}];
Legended[Row[First /@ panelsNielsen], Last[First[panelsNielsen]]]
```

The qutrit pair converts one way only: $\psi$ has the flatter weights, and LOCC can only make weights steeper. The ququart curves cross, so neither state converts to the other, yet with $\chi$ borrowed the crossing disappears and $|\psi\rangle|\chi\rangle$ converts to $|\phi\rangle|\chi\rangle$, which hands $\chi$ back unchanged.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The Holevo bound on accessible information

A sender picks $X\in\{1,2\}$ with $P(X=1)=p$ and sends a qubit $\rho_X$ whose Bloch vector $\vec r_X$ has length $r$; $\vec r_1$ and $\vec r_2$ are $\theta$ apart. The receiver measures and gets an outcome $Y$. How much of $X$ can a measurement recover?

The recovered information is the mutual information $I(X{:}Y)=H(Y)-\sum_xP(x)H(Y|x)$, in bits. Holevo's theorem bounds it for every measurement by

$$
\chi=S\Big(\sum_x P(x)\,\rho_x\Big)-\sum_xP(x)\,S(\rho_x)\;\ge\;I(X{:}Y),
$$

with $S$ the von Neumann entropy. A projective measurement along a unit vector $\hat n$ returns $\pm1$ with probability $P(\pm|x)=\tfrac12(1\pm\vec r_x\cdot\hat n)$.

QF gives $\chi$ (dashed) from the entropies of the two states and of their mixture; the cell plots $I(X{:}Y)$ (solid) for $\hat n$ in the plane of the two Bloch vectors, at angle $\phi$ from $\vec r_1$:

```wl
blochHolevo[r_, \[Theta]_] := QuantumState["BlochVector"[#]] & /@ {{0, 0, r}, r {Sin[\[Theta]], 0, Cos[\[Theta]]}};
chiHolevo[pr_, qs_] := QuantityMagnitude[Total[{pr, 1 - pr} qs]["VonNeumannEntropy"] - {pr, 1 - pr} . Through[qs["VonNeumannEntropy"]]];
h2Holevo[x_] := If[0 < x < 1, -x Log2[x] - (1 - x) Log2[1 - x], 0];
infoHolevo[pr_, r_, \[Theta]_, \[CurlyPhi]_] := With[{plus = (1 + r Cos[{0, \[Theta]} - \[CurlyPhi]])/2}, h2Holevo[{pr, 1 - pr} . plus] - {pr, 1 - pr} . (h2Holevo /@ plus)];
\[Theta]sHolevo = {Pi/4, Pi/2, 3 Pi/4};
panelHolevo[r_, pr_] := Show[
   Plot[Evaluate[infoHolevo[pr, r, #, \[CurlyPhi]] & /@ \[Theta]sHolevo], {\[CurlyPhi], 0, Pi}, PlotRange -> All, PlotLegends -> (Row[{"\[Theta] = ", #}] & /@ \[Theta]sHolevo)],
   Plot[Evaluate[chiHolevo[pr, blochHolevo[r, #]] & /@ \[Theta]sHolevo], {\[CurlyPhi], 0, Pi}, PlotStyle -> Dashed, PlotRange -> All],
   PlotRange -> {0, All}, Frame -> True, FrameLabel -> {"\[Phi]", "bits"}, ImageSize -> 260, FrameTicks -> {{Automatic, None}, {{0, Pi/4, Pi/2, 3 Pi/4, Pi}, None}},
   PlotLabel -> Style[Row[{"r = ", r, ",   p = ", pr}], 13, FontFamily -> "Helvetica"]];
panelsHolevo = Outer[panelHolevo, {1, 1/2}, {1/2, 3/4, 9/10}];
Legended[Grid[Map[First, panelsHolevo, {2}]], Last[panelsHolevo[[1, 1]]]]
```

In every panel each solid curve stays below its dashed $\chi$ for every direction $\phi$. It falls to zero at $\phi=\theta/2$, where $\hat n$ is perpendicular to $\vec r_1-\vec r_2$ and both states give the same outcome probabilities, and peaks near the direction of $\vec r_1-\vec r_2$. For pure states ($r=1$) the peak stays well below $\chi$; mixing them to $r=1/2$ lowers $\chi$ itself, and the best direction then comes close to it.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Classical and quantum capacity of a noisy qubit channel

A qubit crosses a dephasing channel or a depolarizing channel at noise strength $p$. A classical sender encodes one bit with equal priors in the $z$ basis $\{|0\rangle,|1\rangle\}$ or the $x$ basis $\{|+\rangle,|-\rangle\}$; a quantum sender keeps a reference $R$ entangled with the input. How many bits per use does each sender get across?

With Pauli matrices $X,Y,Z$ the two channels are

$$
\mathcal N_{\rm deph}(\rho)=(1-p)\rho+pZ\rho Z,\qquad \mathcal N_{\rm dep}(\rho)=\big(1-\tfrac{3p}4\big)\rho+\tfrac p4\,(X\rho X+Y\rho Y+Z\rho Z).
$$

A pair of input states $\{\rho_1,\rho_2\}$ carries at most its Holevo quantity $\chi=S\big(\mathcal N(\tfrac{\rho_1+\rho_2}2)\big)-\tfrac12\sum_iS(\mathcal N(\rho_i))\le\log_22=1$ bit. The one-use quantum rate is $Q^{(1)}=\max_\rho I_c$, with coherent information $I_c=S(B)-S(RB)$ for the purified input. Every depolarizing input is a rotation of a diagonal one with the same $I_c$; the dephasing channel is degradable, so $I_c$ is concave, and averaging an input over phase rotations keeps its diagonal and cannot lower $I_c$. A diagonal input $\rho_q=\mathrm{diag}(1-q,q)$ therefore suffices, and for a degradable channel $Q^{(1)}$ is the quantum capacity.

QF applies the named channels, with $p$ symbolic, to each basis state and to the purified $\rho_q$; the entropies come from the output spectra in bits, and $I_c$ is maximized over $q$ on a grid:

```wl
codesCap = Map[QuantumState[#]["MatrixState"] &, {{"0", "1"}, {"+", "-"}}, {2}];
entCap[m_, rules_] := With[{ev = Select[Re @ Eigenvalues[N[m /. rules]], # > 10^-12 &]}, -ev . Log2[ev]];
outCap[name_, st_] := Normal[QuantumChannel[name[\[FormalP]], {1}][st]["DensityMatrix"]];
panelCap[name_, title_] := With[{chis = Function[code, outCap[name, #] & /@ Prepend[code, Total[code]/2]] /@ codesCap,
    rb = outCap[name, QuantumState[DiagonalMatrix[{1 - \[FormalQ], \[FormalQ]}]]["Purify"]], ps = Subdivide[0, 1, 40]},
   With[{b = Normal[QuantumPartialTrace[QuantumState[rb], {2}]["DensityMatrix"]]},
    ListLinePlot[Append[Table[{#, entCap[c[[1]], \[FormalP] -> #] - (entCap[c[[2]], \[FormalP] -> #] + entCap[c[[3]], \[FormalP] -> #])/2} & /@ ps, {c, chis}],
      {#, Max[Table[entCap[b, {\[FormalP] -> #, \[FormalQ] -> q}] - entCap[rb, {\[FormalP] -> #, \[FormalQ] -> q}], {q, 0, 1, 1/20}]]} & /@ ps],
     PlotStyle -> {Automatic, Dashed, Directive[Thick, Dotted]}, PlotRange -> {-0.02, 1.02}, Frame -> True, FrameLabel -> {"p", "bits per use"}, ImageSize -> 320,
     PlotLegends -> {"\[Chi], z basis", "\[Chi], x basis", Superscript["Q", "(1)"]}, PlotLabel -> Style[title, 13, FontFamily -> "Helvetica"]]]];
panelsCap = MapThread[panelCap, {{"PhaseFlip", "Depolarizing"}, {"dephasing", "depolarizing"}}];
Legended[Row[First /@ panelsCap], Last[First[panelsCap]]]
```

Dephasing never touches the $z$-basis states, which reach the one-bit ceiling at every $p$, while the $x$-basis states and the quantum rate fall together to zero at $p=1/2$: there a full classical bit crosses a channel that carries no quantum information. The depolarizing channel treats both bases alike, and its $Q^{(1)}$ hits zero at a noise level where classical bits still get through. That curve is a one-use rate, so for this channel it is only a lower bound on the quantum capacity.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Randomized benchmarking: the decay of sequence fidelity

A qubit starts in $|0\rangle$, undergoes $m$ Clifford gates drawn at random and then the Clifford that undoes their product, and is measured; every gate is followed by the same noise. The noise is amplitude damping, which drops $|1\rangle$ to $|0\rangle$ with probability $\gamma$, or a coherent over-rotation $e^{-i\varepsilon X/2}$ with $\varepsilon$ chosen so that both noises have the same average gate fidelity $F$. How does the survival probability, the chance of finding $|0\rangle$ at the end, decay with $m$?

Any qubit channel moves the Bloch vector as $\vec r\mapsto M\vec r+\vec c$. The 24 Cliffords act on $\vec r$ as the rotations of a cube, and averaging the noise over them, which the random sequence does, leaves a depolarizing map that shrinks $\vec r$ by $p=\operatorname{Tr}M/3=2F-1$. With the final noise left over, the mean survival is exactly

$$
P(m)=\tfrac12\big(1+c_z+M_{zz}\,p^m\big)=A\,p^m+B .
$$

For damping $M=\mathrm{diag}(\sqrt{1-\gamma},\sqrt{1-\gamma},1-\gamma)$ and $\vec c=(0,0,\gamma)$; for the rotation $\operatorname{Tr}M=1+2\cos\varepsilon$, so equal $F$ means $\cos\varepsilon=\sqrt{1-\gamma}-\gamma/2$.

QF reads $M$ and $\vec c$ off each noise by applying it to Bloch states; the cell averages the survival over 50 random sequences at each $m$, in three independent runs (markers), and draws $P(m)$ (line):

```wl
cliffordRB = Select[Flatten[Outer[DiagonalMatrix[#1] . #2 &, Tuples[{1, -1}, 3], Permutations[IdentityMatrix[3]], 1], 1], Det[#] == 1 &];
blochMapRB[f_] := With[{c = f[QuantumState["BlochVector"[{0, 0, 0}]]]["BlochVector"]}, {Transpose[f[QuantumState["BlochVector"[#]]]["BlochVector"] - c & /@ IdentityMatrix[3]], c}];
survivalRB[{mm_, c_}, m_] := With[{seq = RandomChoice[cliffordRB, m], step = mm . (#2 . #1) + c &},
   (1 + Last[step[Fold[step, {0, 0, 1}, seq], Transpose[Dot @@ Reverse[seq]]]])/2];
panelRB[noise_, title_] := With[{map = Chop[Re[N[blochMapRB[noise]]]], ms = Range[1, 61, 4]},
   Show[ListPlot[Table[BlockRandom[SeedRandom[run]; {#, Mean[Table[survivalRB[map, #], 50]]} & /@ ms], {run, 3}], PlotMarkers -> {Automatic, 8}, PlotRange -> All],
    Plot[(1 + map[[2, 3]] + map[[1, 3, 3]] (Tr[map[[1]]]/3)^m)/2, {m, 0, 61}, PlotRange -> All],
    PlotRange -> {All, {0.4, 1.01}}, Axes -> False, Frame -> True, FrameLabel -> {"m", "survival"}, ImageSize -> 300,
    PlotLabel -> Style[Column[{title, Row[{"p = Tr M/3 = ", NumberForm[Tr[map[[1]]]/3, 3]}]}, Alignment -> Center], 13, FontFamily -> "Helvetica"]]];
\[Gamma]sRB = {1/20, 1/5};
Grid[{panelRB[QuantumChannel["AmplitudeDamping"[#]], Row[{"damping, \[Gamma] = ", #}]] & /@ \[Gamma]sRB,
   panelRB[QuantumOperator["RX"[ArcCos[Sqrt[1 - #] - #/2]]], Row[{"over-rotation matched to \[Gamma] = ", #}]] & /@ \[Gamma]sRB}]
```

In each column the mean survival of the two noises falls at the same rate: averaging over random Cliffords keeps only $\operatorname{Tr}M$ of the noise. The floors differ: damping pulls toward $|0\rangle$ and holds long sequences above $1/2$, while the rotation, which has no shift $\vec c$, settles at $1/2$. The rotation's runs scatter more, because within one sequence coherent errors can add up or cancel, and only the average over sequences follows the exponential.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Bell-CHSH violation and the Tsirelson bound for a family of states

Two qubits share $\rho=p\,|\psi_\theta\rangle\langle\psi_\theta|+(1-p)\,I/4$, a partly entangled pure state $|\psi_\theta\rangle=\cos\theta|00\rangle+\sin\theta|11\rangle$ mixed with white noise. Alice measures $\vec a_1\cdot\vec\sigma$ or $\vec a_2\cdot\vec\sigma$, Bob $\vec b_1\cdot\vec\sigma$ or $\vec b_2\cdot\vec\sigma$, all with outcomes $\pm1$. For which $p$ and $\theta$ can the correlations violate the CHSH inequality?

Local models obey the CHSH inequality $|\langle\mathcal B\rangle|\le2$ for

$$
\mathcal B=\vec a_1\cdot\vec\sigma\otimes(\vec b_1+\vec b_2)\cdot\vec\sigma+\vec a_2\cdot\vec\sigma\otimes(\vec b_1-\vec b_2)\cdot\vec\sigma ,
$$

while quantum states reach at most $2\sqrt2$ (Tsirelson). With the correlation matrix $T_{ij}=\operatorname{Tr}\rho\,\sigma_i\otimes\sigma_j$, the Horodecki theorem gives the best value over all settings, $\max\langle\mathcal B\rangle=2\sqrt{t_1+t_2}$, where $t_1,t_2$ are the two largest eigenvalues of $T^{\mathsf T}T$. For two qubits $\rho$ is entangled exactly when its partial transpose has a negative eigenvalue (Peres-Horodecki).

QF builds $\rho$ with $p$ and $\theta$ symbolic; the cell gives the optimal $\langle\mathcal B\rangle$ in closed form, compares it with the settings that are optimal for a Bell state ($\vec a_{1,2}=\hat z,\hat x$, $\vec b_{1,2}=(\hat z\pm\hat x)/\sqrt2$, dashed), and asks `Reduce` where $\rho$ violates CHSH and where its partial transpose has a negative eigenvalue:

```wl
rhoCHSH = \[FormalP] QuantumState[{Cos[\[FormalTheta]], 0, 0, Sin[\[FormalTheta]]}, 2]["MatrixState"] + (1 - \[FormalP]) QuantumState[IdentityMatrix[4]/4, {2, 2}];
corrCHSH = Simplify[Table[Tr[Normal[rhoCHSH["DensityMatrix"]] . KroneckerProduct[PauliMatrix[i], PauliMatrix[j]]], {i, 3}, {j, 3}], 0 < \[FormalTheta] < Pi/2];
bestCHSH = With[{ev = Eigenvalues[Transpose[corrCHSH] . corrCHSH]}, Simplify[2 Sqrt[Total[ev] - Min[ev]], 0 < \[FormalP] <= 1 && 0 < \[FormalTheta] < Pi/2]];
bellCHSH = Sqrt[2] (corrCHSH[[3, 3]] + corrCHSH[[1, 1]]);
ptCHSH = Simplify[QuantumPartialTranspose[rhoCHSH, {2}]["Eigenvalues"], 0 < \[FormalTheta] < Pi/2];
panelCHSH[\[Theta]_] := With[{best = bestCHSH /. \[FormalTheta] -> \[Theta], bell = bellCHSH /. \[FormalTheta] -> \[Theta]},
   With[{violate = Reduce[best > 2 && 0 <= \[FormalP] <= 1, \[FormalP], Reals], entangled = Reduce[Min[ptCHSH /. \[FormalTheta] -> \[Theta]] < 0 && 0 <= \[FormalP] <= 1, \[FormalP], Reals]},
    Plot[{best, bell}, {\[FormalP], 0, 1}, PlotRange -> {0, 3}, PlotStyle -> {Automatic, Dashed}, PlotLegends -> {"best settings", "Bell-state settings"},
     GridLines -> {{entangled[[1]], violate[[1]]}, {2, 2 Sqrt[2]}}, GridLinesStyle -> Dashed, Frame -> True, FrameLabel -> {"p", "\[LeftAngleBracket]\[ScriptCapitalB]\[RightAngleBracket]"}, ImageSize -> 300,
     PlotLabel -> Style[Column[{Row[{"\[Theta] = ", \[Theta]}], Row[{"entangled for p > ", FullSimplify[entangled[[1]]]}], Row[{"violates CHSH for p > ", FullSimplify[violate[[1]]]}]}, Alignment -> Center], 13, FontFamily -> "Helvetica"]]]];
panelsCHSH = panelCHSH /@ {Pi/4, Pi/8, Pi/16};
Labeled[Legended[Row[First /@ panelsCHSH], Last[First[panelsCHSH]]], Style[Row[{"max\[LeftAngleBracket]\[ScriptCapitalB]\[RightAngleBracket] = ", TraditionalForm[bestCHSH /. {\[FormalP] -> p, \[FormalTheta] -> \[Theta]}]}], 13, FontFamily -> "Helvetica"], Top]
```

The optimal value reaches the Tsirelson value $2\sqrt2$ (upper horizontal line) only for the Bell state at $p=1$, and crosses the local bound $2$ (lower horizontal line) at the right-hand vertical line, always at a larger $p$ than the left-hand line where entanglement begins: between the two, the state is entangled yet no choice of CHSH settings detects it. A less entangled $|\psi_\theta\rangle$ needs less noise to stop violating. The Bell-state settings fall short once $\theta\neq\pi/4$, and at $\theta=\pi/16$ they miss a violation that other settings find.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Proving a circuit's phase approximation correct with Resolve

A phase circuit must apply $|x\rangle\mapsto e^{2\pi i/\sqrt x}|x\rangle$. Writing $x=4^km$ with $1\le m<4$ reduces the job to $1/\sqrt m$, and the circuit evaluates it on $3\cdot2^r$ intervals $m=a+t$, $0\le t\le h=2^{-r}$, as a quadratic $p(t)$ whose coefficients are integers over $2^w$, that is, fixed-point numbers with $w$ fractional bits. Is the relative error below a bound $\eta$, proposed on each interval from finitely many samples, for every $m$ in it?

On each interval the claim is a polynomial inequality in one real variable,

$$
\forall\,t\in[0,h]:\quad p(t)>0\ \wedge\ (1-\eta)^2\le(a+t)\,p(t)^2\le(1+\eta)^2\iff\big|\sqrt{a+t}\,p(t)-1\big|\le\eta ,
$$

for $0\le\eta<1$, which `Resolve` decides exactly over the reals, with rational $a$, $h$, $\eta$ and coefficients and no rounding. The quadratic interpolates $1/\sqrt{a+t}$ at the three Chebyshev nodes of the interval before its coefficients are rounded to $w$ bits, and $\eta$ is proposed from the error at 201 sample points plus 2%, rounded up to a multiple of $2^{-40}$.

The cell builds every interval for each $r$ and $w$, plots $|\sqrt m\,p-1|$ against $m$ (solid) with each interval's $\eta$ (dashed), and reports whether `Resolve` proves the claim on all intervals:

```wl
intervalsCert[r_, w_] := With[{h = 2^-r}, Table[With[{a = 1 + i h, nodes = h (1 - Cos[Pi (2 Range[3] - 1)/6])/2},
     With[{poly = Round[2^w CoefficientList[InterpolatingPolynomial[Transpose[{nodes, N[1/Sqrt[a + nodes], 30]}], \[FormalT]], \[FormalT]]] . \[FormalT]^Range[0, 2]/2^w},
      With[{eta = Ceiling[1.02 2^40 Max[Abs[Sqrt[a + #] (poly /. \[FormalT] -> #) - 1] & /@ Subdivide[0., N[h], 200]]]/2^40},
       <|"a" -> a, "h" -> h, "poly" -> poly, "eta" -> eta,
        "proof" -> Resolve[ForAll[\[FormalT], 0 <= \[FormalT] <= h, poly > 0 && (1 - eta)^2 <= (a + \[FormalT]) poly^2 <= (1 + eta)^2], Reals]|>]]], {i, 0, 3/h - 1}]];
panelCert[r_, w_] := With[{ivs = intervalsCert[r, w]},
   Show[Plot[Evaluate[Piecewise[N[{Abs[Sqrt[m] (#poly /. \[FormalT] -> m - #a) - 1], #a <= m < #a + #h}] & /@ ivs]], {m, 1, 4}, PlotRange -> All, PlotPoints -> 100, Exclusions -> None],
    Plot[Evaluate[Piecewise[N[{#eta, #a <= m < #a + #h}] & /@ ivs]], {m, 1, 4}, PlotStyle -> Dashed, PlotRange -> All, Exclusions -> None],
    PlotRange -> {0, All}, Axes -> False, Frame -> True, FrameLabel -> {"m", "|relative error|"}, ImageSize -> 280,
    PlotLabel -> Style[Column[{Row[{"r = ", r, ",   w = ", w, " bits"}], Row[{"Resolve, all ", Length[ivs], " intervals: ", And @@ Lookup[ivs, "proof"]}]}, Alignment -> Center], 13, FontFamily -> "Helvetica"]]];
Grid[Outer[panelCert, {1, 2}, {8, 12, 16}]]
```

`Resolve` proves the bound on every interval in every panel, so the error stays under the dashed line for every real $m$, not only at the samples that proposed $\eta$. With 8-bit coefficients rounding sets the error, which jumps from interval to interval and barely shrinks when the intervals halve. By 16 bits rounding is negligible and the error is the interpolation's ripple, which shrinks as the intervals halve and as $m$ grows.
