```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Coherent information of amplitude damping, maximized over the input

A qubit holds $\rho_p=(1-p)|0\rangle\langle0|+p\,|1\rangle\langle1|$, with $|1\rangle$ excited, and a reference qubit $R$ purifies it. Only the qubit crosses an amplitude-damping channel, which drops $|1\rangle$ to $|0\rangle$ with probability $\gamma$. Which input population $p$ carries the most quantum information across?

The channel and the coherent information, in bits, are

$$
\mathcal N_\gamma(\rho)=K_0\rho K_0^\dagger+K_1\rho K_1^\dagger,\quad K_0=\begin{pmatrix}1&0\\0&\sqrt{1-\gamma}\end{pmatrix},\quad K_1=\begin{pmatrix}0&\sqrt\gamma\\0&0\end{pmatrix},\qquad I_c=S(B)-S(RB),
$$

with $B$ the channel output and $S$ the von Neumann entropy.

QF computes $I_c$ with $p$ and $\gamma$ symbolic and tests it against the binary-entropy form $h_2((1-\gamma)p)-h_2(\gamma p)$; `ArgMax` and `MaxValue` then give the best input at each $\gamma<1/2$ (dots):

```wl
icAD = QuantityMagnitude @ FullSimplify[With[{out = QuantumChannel["AmplitudeDamping"[γ], {1}][QuantumState[DiagonalMatrix[{1 - p, p}]]["Purify"]]},
     QuantumPartialTrace[out, {2}]["VonNeumannEntropy"] - out["VonNeumannEntropy"]], 0 < p < 1 && 0 < γ < 1];
h2AD[x_] := -x Log2[x] - (1 - x) Log2[1 - x];
icClosedQ = FullSimplify[icAD == h2AD[(1 - γ) p] - h2AD[γ p], 0 < p < 1 && 0 < γ < 1];
γsAD = {1/10, 1/4, 1/3, 1/2, 3/4};
bestAD = Through[{ArgMax, MaxValue}[{h2AD[(1 - #) p] - h2AD[# p], 0 <= p <= 1}, p]] & /@ Select[γsAD, # < 1/2 &];
Show[Plot[Evaluate[h2AD[(1 - #) p] - h2AD[# p] & /@ γsAD], {p, 0, 1}, PlotRange -> All, PlotLegends -> (Row[{"γ = ", #}] & /@ γsAD)],
  ListPlot[List /@ N[bestAD], PlotMarkers -> {Automatic, 10}],
  PlotRange -> All, Frame -> True, FrameLabel -> {"p", Row[{Subscript["I", "c"], " (bits)"}]},
  PlotLabel -> Row[{Subscript["I", "c"], " = ", Subscript["h", 2], "((1-γ)p) - ", Subscript["h", 2], "(γp):  ", icClosedQ}], ImageSize -> 420]
```

The output keeps the population $(1-\gamma)p$ and the environment receives $\gamma p$. Below $\gamma=1/2$ the best input has $p<1/2$, fewer excitations than the maximally mixed state. At $\gamma=1/2$, $I_c$ vanishes for every input, and at $\gamma=3/4$ no input gives a positive $I_c$: the environment then gets the larger share.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Data-processing inequality under named channels

Two qubit states $\rho$ and $\sigma$ have Bloch vectors of length $3/5$ and $4/5$ at an angle $\theta$, so they do not commute. Both pass through the same channel at noise strength $p$. Can any $p$ make them easier to tell apart?

Their distinguishability is the relative entropy $S(\rho\|\sigma)=\operatorname{Tr}\rho\,(\log_2\rho-\log_2\sigma)$, and the data-processing inequality says $S(\mathcal N_p(\rho)\|\mathcal N_p(\sigma))\le S(\rho\|\sigma)$ for every channel $\mathcal N_p(\rho)=\sum_k K_k\rho K_k^\dagger$. The four channels have Kraus operators

$$
\text{depolarizing: }\sqrt{1-\tfrac{3p}{4}}\,I,\ \tfrac{\sqrt p}{2}X,\ \tfrac{\sqrt p}{2}Y,\ \tfrac{\sqrt p}{2}Z;\qquad \text{bit flip: }\sqrt{1-p}\,I,\ \sqrt p\,X;\qquad \text{phase flip: }\sqrt{1-p}\,I,\ \sqrt p\,Z;\qquad \text{phase damping: }\begin{pmatrix}1&0\\0&\sqrt{1-p}\end{pmatrix},\ \begin{pmatrix}0&0\\0&\sqrt p\end{pmatrix}.
$$

QF gives the processed relative entropy as a closed form in $p$, with $\rho$ along $z$ and $\sigma$ in the $xz$ plane, and `Reduce` searches all of $0<p<1$ for a value above $S(\rho\|\sigma)$ (dashed):

```wl
relEnt[a_, b_] := FullSimplify[QuantityMagnitude @ QuantumDistance[a, b, "RelativeEntropy"], 0 < p < 1];
channelsDPI = {"Depolarizing", "BitFlip", "PhaseFlip", "PhaseDamping"};
panelDPI[θ_] := With[{pair = QuantumState["BlochVector"[#]] & /@ {{0, 0, 3/5}, 4/5 {Sin[θ], 0, Cos[θ]}}},
   With[{before = relEnt @@ pair, after = relEnt @@ (QuantumChannel[#[p]] /@ pair) & /@ channelsDPI},
    Plot[Evaluate[after], {p, 0, 1}, PlotRange -> {0, All}, PlotLegends -> channelsDPI, GridLines -> {None, {before}}, GridLinesStyle -> Directive[Thick, Dashed],
     Frame -> True, FrameLabel -> {"p", "S(\[ScriptCapitalN](ρ)||\[ScriptCapitalN](σ)) (bits)"}, ImageSize -> 300,
     PlotLabel -> Row[{"θ = ", θ, ",   p raising S: ", Or @@ (Reduce[# > before && 0 < p < 1, p, Reals] & /@ after)}]]]];
panelsDPI = panelDPI /@ {Pi/3, Pi/2, 2 Pi/3};
Legended[Row[First /@ panelsDPI], Last[First[panelsDPI]]]
```

For every channel and each of the three pairs, `Reduce` finds no noise strength that raises the relative entropy, so each curve stays under the dashed line for all $p$, not just at the points drawn. Depolarizing noise only lowers it, down to zero at $p=1$ where both states become the maximally mixed state; phase damping only lowers it too. Bit flip and phase flip become the unitaries $X$ and $Z$ at $p=1$, so their curves climb back exactly to the dashed line there.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Two qubits have no bound entanglement

A bound-entangled state is entangled, yet its partial transpose $\rho^{T_B}$ is still positive, so no local protocol can distill a Bell pair from it. Can a two-qubit state be bound entangled?

Write every two-qubit state as $\rho=\sum_{i,j=0}^{3}c_{ij}\,\sigma_i\otimes\sigma_j$ with 16 real symbols $c_{ij}$. Wootters' concurrence $C=\max(0,\lambda_1-\lambda_2-\lambda_3-\lambda_4)$ vanishes exactly on separable states, where $\lambda_1\ge\dots\ge\lambda_4$ are the square roots of the eigenvalues of $R=\rho\,(\sigma_y\otimes\sigma_y)\,\rho^*(\sigma_y\otimes\sigma_y)$. Since $\sum_k\lambda_k^2=\operatorname{Tr}R$, $\sum_k\lambda_k^4=\operatorname{Tr}R^2$ and $\prod_k\lambda_k=\det\rho$, the claim is

$$
16\det\rho^{T_B}=(\operatorname{Tr}R)^2-2\operatorname{Tr}R^2+8\det\rho=\prod_{k=1}^{4}\Big(2\lambda_k-\sum_j\lambda_j\Big),
$$

and a nonnegative product forces $\lambda_1\le\lambda_2+\lambda_3+\lambda_4$, that is $C=0$.

The cell tests the first equality in all 16 symbols $c_{ij}$ and the second in four symbols $\lambda_k$, asks `Resolve` for the second statement over every ordered quadruple $\lambda_1\ge\lambda_2\ge\lambda_3\ge\lambda_4\ge0$, and draws, in units of $\lambda_1$, where $C>0$ and where the product is nonnegative:

```wl
paulis2q = Flatten[Outer[KroneckerProduct, #, #, 1] &[PauliMatrix[{0, 1, 2, 3}]], 1];
yyPPT = KroneckerProduct[PauliMatrix[2], PauliMatrix[2]];
ρ2q = Array[c, 16] . paulis2q;
r2q = ρ2q . yyPPT . (Array[c, 16] . Conjugate[paulis2q]) . yyPPT;
detTB = Det[ArrayReshape[Transpose[ArrayReshape[ρ2q, {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}]];
λsPPT = Array[l, 4];
identityQ = Expand[16 detTB - (Tr[r2q]^2 - 2 Tr[r2q . r2q] + 8 Det[ρ2q])] === 0 &&
   Expand[Times @@ (2 λsPPT - Total[λsPPT]) - (Total[λsPPT^2]^2 - 2 Total[λsPPT^4] + 8 Times @@ λsPPT)] === 0;
noBoundQ = With[{λ = λsPPT}, Resolve[ForAll[λ, GreaterEqual @@ Append[λ, 0] && Times @@ (2 λ - Total[λ]) >= 0, First[λ] - Total[Rest[λ]] <= 0], Reals]];
panelPPT[l4_] := RegionPlot[{v <= u && 1 - u - v - l4 > 0, v <= u && Times @@ (2 {1, u, v, l4} - (1 + u + v + l4)) >= 0}, {u, l4, 1}, {v, l4, 1}, PlotPoints -> 60, PlotLegends -> {"C > 0", Row[{"16 det ", Superscript["ρ", Subscript["T", "B"]], " ≥ 0"}]},
   Frame -> True, FrameLabel -> {Row[{Subscript["λ", 2], "/", Subscript["λ", 1]}], Row[{Subscript["λ", 3], "/", Subscript["λ", 1]}]},
   PlotLabel -> Row[{Subscript["λ", 4], "/", Subscript["λ", 1], " = ", l4}], ImageSize -> 260];
panelsPPT = panelPPT /@ {0, 1/10, 1/5};
Labeled[Legended[Row[First /@ panelsPPT], Last[First[panelsPPT]]], Style[Row[{"both identities: ", identityQ, ",   Resolve: ", noBoundQ}], 13, FontFamily -> "Helvetica"], Top]
```

The two regions never overlap and together fill every panel: wherever the determinant of the partial transpose is nonnegative, the concurrence is zero, and `Resolve` makes this hold for every quadruple, not just the slices drawn. A positive partial transpose has a nonnegative determinant, so it forces $C=0$, which is separability. A two-qubit state is either separable or has a negative partial transpose, and nothing is left to be bound.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Robertson-Schrödinger slack for every qubit state and observable pair

Two observables $A=\vec a\cdot\vec\sigma$ and $B=\vec b\cdot\vec\sigma$ act on a qubit $\rho=\tfrac12(I+\vec r\cdot\vec\sigma)$ with $|\vec r|\le1$; a multiple of the identity added to either drops out of every term below, so these are all qubit observables. How far above the bound does each state sit?

With $\langle O\rangle=\operatorname{Tr}O\rho$ and $\sigma_O^2=\langle O^2\rangle-\langle O\rangle^2$, the Robertson-Schrödinger relation is

$$
\sigma_A^2\,\sigma_B^2\;\ge\;\Big(\tfrac12\langle\{A,B\}\rangle-\langle A\rangle\langle B\rangle\Big)^2+\Big(\tfrac1{2i}\langle[A,B]\rangle\Big)^2 .
$$

The cell computes both sides in all nine real symbols, tests the slack (left minus right) against $|\vec a\times\vec b|^2(1-|\vec r|^2)$, and asks `Reduce` for any violation. The figure takes unit $\vec a$ along $z$, unit $\vec b$ at angle $\gamma$ from it, and $\vec r$ at angle $\beta$ from $z$ in the same plane; left side solid, right side dashed, and at $|\vec r|=1$ the dashed curve lies exactly under the solid one:

```wl
{lhsRS, rhsRS} = With[{ρq = (IdentityMatrix[2] + {r1, r2, r3} . PauliMatrix[{1, 2, 3}])/2, opA = {a1, a2, a3} . PauliMatrix[{1, 2, 3}], opB = {b1, b2, b3} . PauliMatrix[{1, 2, 3}]},
   With[{ev = Tr[# . ρq] &}, Expand[{(ev[opA . opA] - ev[opA]^2) (ev[opB . opB] - ev[opB]^2),
      (ev[(opA . opB + opB . opA)/2] - ev[opA] ev[opB])^2 + (ev[opA . opB - opB . opA]/(2 I))^2}]]];
factorsQ = Expand[lhsRS - rhsRS - Cross[{a1, a2, a3}, {b1, b2, b3}] . Cross[{a1, a2, a3}, {b1, b2, b3}] (1 - {r1, r2, r3} . {r1, r2, r3})] === 0;
violationRS = Reduce[lhsRS - rhsRS < 0 && r1^2 + r2^2 + r3^2 <= 1, {r1, r2, r3, a1, a2, a3, b1, b2, b3}, Reals];
rLengthsRS = {0, 1/2, 1};
panelRS[γ_] := With[{sides = {lhsRS, rhsRS} /. {a1 -> 0, a2 -> 0, a3 -> 1, b1 -> Sin[γ], b2 -> 0, b3 -> Cos[γ]}},
   With[{curves = Table[sides /. {r1 -> ρr Sin[β], r2 -> 0, r3 -> ρr Cos[β]}, {ρr, rLengthsRS}]},
    Show[Plot[Evaluate[curves[[All, 1]]], {β, 0, 2 Pi}, PlotRange -> All, PlotLegends -> (Row[{"|r| = ", #}] & /@ rLengthsRS)],
     Plot[Evaluate[curves[[All, 2]]], {β, 0, 2 Pi}, PlotRange -> All, PlotStyle -> Dashed],
     PlotRange -> {-0.02, 1.02}, Frame -> True, FrameLabel -> {"β", Row[{Superscript[Subscript["σ", "A"], 2], Superscript[Subscript["σ", "B"], 2], " (solid),  right side (dashed)"}]}, PlotLabel -> Row[{"γ = ", γ}], ImageSize -> 300]]];
panelsRS = panelRS /@ {Pi/6, Pi/3, Pi/2};
Labeled[Legended[Row[First /@ panelsRS], Last[First[panelsRS]]], Style[Row[{"slack = |a×b|²(1-|r|²): ", factorsQ, ",   violation found by Reduce: ", violationRS}], 13, FontFamily -> "Helvetica"], Top]
```

For $|\vec r|>0$ both sides swing as the state turns, but the gap between them does not move: it depends on the state only through its mixedness $1-|\vec r|^2$ and on the observables only through how far they are from commuting. The gap closes for every pure state, at every $\beta$, and `Reduce` finds no state and no pair of observables anywhere in the nine-dimensional space that violates the relation.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Helstrom bound by Reduce over every measurement

One of two qubit states arrives, $\rho_1$ with prior $p$ or $\rho_2$ with prior $1-p$; their Bloch vectors have equal length $r$ and meet at angle $\theta$. One two-outcome measurement $\{E_1,E_2\}$ must name it, with $E_2$ meaning "it was $\rho_2$". What is the smallest possible error?

$$
P_{\rm err}=p\operatorname{Tr}\rho_1E_2+(1-p)\operatorname{Tr}\rho_2E_1=(1-p)+\operatorname{Tr}\Lambda E_2,\qquad \Lambda=p\rho_1-(1-p)\rho_2,\qquad E_1=I-E_2,\quad 0\le E_2\le I .
$$

A unitary maps the set of effects to itself, so in the eigenbasis of $\Lambda=\mathrm{diag}(l_1,l_2)$ the most general effect is $E_2=\begin{pmatrix}a&b+ic\\b-ic&e\end{pmatrix}$ with $E_2\ge0$ and $I-E_2\ge0$.

`Reduce` finds every $k$ with $\operatorname{Tr}\Lambda E_2\ge k$ for all such effects, and $(1-p)$ plus the largest $k$ is the lowest error (dashed). Solid: the error of the projective measurement $E_2=\tfrac12(I+\hat m\cdot\vec\sigma)$, with $\hat m$ at angle $\phi$ from $\rho_1$'s Bloch vector in the plane of both:

```wl
helstromK = Reduce[ForAll[{ea, eb, ec, ee}, 0 <= ea <= 1 && 0 <= ee <= 1 && ea ee >= eb^2 + ec^2 && (1 - ea) (1 - ee) >= eb^2 + ec^2,
    Tr[DiagonalMatrix[{l1, l2}] . {{ea, eb + I ec}, {eb - I ec, ee}}] >= k], k, Reals];
blochRho[v_] := (IdentityMatrix[2] + v . PauliMatrix[{1, 2, 3}])/2;
θsHel = {Pi/4, Pi/2, Pi};
panelHel[r_, pr_] := With[{ρ1 = blochRho[{0, 0, r}], ρ2s = blochRho[r {Sin[#], 0, Cos[#]}] & /@ θsHel, e2 = blochRho[{Sin[φ], 0, Cos[φ]}]},
   Show[Plot[Evaluate[pr Tr[ρ1 . e2] + (1 - pr) Tr[# . (IdentityMatrix[2] - e2)] & /@ ρ2s], {φ, 0, 2 Pi}, PlotRange -> All, PlotLegends -> (Row[{"θ = ", #}] & /@ θsHel)],
    Plot[Evaluate[(1 - pr) + MaxValue[{k, helstromK /. Thread[{l1, l2} -> Eigenvalues[pr ρ1 - (1 - pr) #]]}, k] & /@ ρ2s], {φ, 0, 2 Pi}, PlotStyle -> Dashed],
    PlotRange -> {0, 1}, Axes -> False, Frame -> True, FrameLabel -> {"φ", Subscript["P", "err"]}, PlotLabel -> Row[{"r = ", r, ",   p = ", pr}], ImageSize -> 280]];
panelsHel = Outer[panelHel, {1, 1/2}, {1/2, 3/4, 9/10}];
Labeled[Legended[Grid[Map[First, panelsHel, {2}]], Last[panelsHel[[1, 1]]]], Style[Row[{"Reduce: ", TraditionalForm[helstromK /. {l1 -> Subscript["l", 1], l2 -> Subscript["l", 2]}]}], 13, FontFamily -> "Helvetica"], Top]
```

The largest $k$ is the sum of the negative eigenvalues of $\Lambda$, so with $\operatorname{Tr}\Lambda=2p-1$ the lowest error is the Helstrom bound $\tfrac12(1-\|\Lambda\|_1)$. Whenever $\Lambda$ has a negative eigenvalue, the best rank-one projector touches the bound. For mixed states with a strong prior $\Lambda$ has no negative eigenvalue, the optimum is $E_2=0$ (always guess $\rho_1$), and a rank-one projector reaches it only when $\Lambda$ has a zero eigenvalue.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Werner state entanglement threshold by Reduce

Two qudits of dimension $d$ share a Werner state: weight $p$ on the symmetric subspace and $1-p$ on the antisymmetric one, each spread evenly. At which $p$ does its partial transpose turn negative in each dimension $d$?

With $F$ the swap of the two qudits,

$$
\rho_W(p)=\frac{p}{d_s}\,\Pi_s+\frac{1-p}{d_a}\,\Pi_a,\qquad \Pi_{s,a}=\frac{I\pm F}{2},\qquad d_{s,a}=\frac{d(d\pm1)}{2},
$$

so for $d=2$, $\Pi_a$ is the singlet projector. A negative eigenvalue of the partial transpose $\rho_W^{T_B}$ proves the state entangled.

QF builds $\rho_W(p)$ with $p$ symbolic, takes its partial transpose on the second qudit, and `Reduce` finds where the smallest eigenvalue is negative. The figure shows the two distinct eigenvalues of $\rho_W^{T_B}$, the one that stays positive dashed:

```wl
ptEigenW = AssociationMap[Union @ Simplify @ QuantumPartialTranspose[QuantumState["Werner"[p, #]], {2}]["Eigenvalues"] &, {2, 3, 4}];
entangledW = Reduce[Min[#] < 0 && 0 <= p <= 1, p, Reals] & /@ ptEigenW;
{ptLowW, ptHighW} = Transpose[SortBy[#, # /. p -> 0 &] & /@ Values[ptEigenW]];
Show[Plot[Evaluate[ptLowW], {p, 0, 1}, PlotRange -> All, PlotLegends -> (Row[{"d = ", #}] & /@ Keys[ptEigenW])],
  Plot[Evaluate[ptHighW], {p, 0, 1}, PlotRange -> All, PlotStyle -> Dashed],
  PlotRange -> All, Frame -> True, FrameLabel -> {"p", Row[{"eigenvalues of ", Superscript["ρ", Subscript["T", "B"]]}]},
  PlotLabel -> Column[KeyValueMap[Row[{"d = ", #1, ":  negative for  ", #2}] &, entangledW], Spacings -> 1], ImageSize -> 420]
```

In every dimension the same eigenvalue crosses zero at the same weight, so the threshold does not depend on $d$. Below it the state is entangled. For two qubits a positive partial transpose also means separable (the no-bound-entanglement slide), so there $p<1/2$ is exactly the entangled range.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Entropic uncertainty for a qubit in two mutually unbiased bases

A qubit with Bloch vector $(x,y,z)$, $x^2+y^2+z^2\le1$, is measured in the $X$ basis or in the $Z$ basis. The two bases are mutually unbiased: every overlap $|\langle\pm|0\rangle|^2$, $|\langle\pm|1\rangle|^2$ is $1/2$. Can any state be sharp in both?

The Maassen-Uffink relation bounds the Shannon entropies of the two outcome distributions for every state by the largest overlap $c^2=1/2$:

$$
H(X)+H(Z)\;\ge\;-\log_2c^2=1\ \text{bit},\qquad H(X)=h(x),\quad H(Z)=h(z),\quad h(u)=-\tfrac{1+u}{2}\log_2\tfrac{1+u}{2}-\tfrac{1-u}{2}\log_2\tfrac{1-u}{2}.
$$

`Reduce` shows that $h(u)<1-u^2$ has no solution on $-1<u<1$, so $H(X)+H(Z)\ge2-x^2-z^2\ge1$, and finds where $h(u)=1-u^2$ there; at $u=\pm1$ both sides are $0$. The figure draws $H(X)+H(Z)$ over every $(x,z)$ a state can have, with the four states on the bound marked:

```wl
hMU[u_] := -(1 + u)/2 Log2[(1 + u)/2] - (1 - u)/2 Log2[(1 - u)/2];
stepMU = Reduce[hMU[u] < 1 - u^2 && -1 < u < 1, u, Reals];
eqMU = Reduce[hMU[u] == 1 - u^2 && -1 < u < 1, u, Reals];
Show[ContourPlot[hMU[u] + hMU[w], {u, -1, 1}, {w, -1, 1}, RegionFunction -> (#1^2 + #2^2 < 1 &), Contours -> Range[1.1, 1.9, 0.1],
   PlotPoints -> 60, PlotLegends -> Automatic],
  ListPlot[{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}, PlotStyle -> PointSize[0.025]],
  Frame -> True, FrameLabel -> {"x", "z"}, PlotRange -> {{-1.05, 1.05}, {-1.05, 1.05}},
  PlotLabel -> Column[{"H(X) + H(Z) in bits", Row[{"on -1 < u < 1:   h(u) < 1 - ", Superscript["u", 2], ": ", stepMU, ";   h(u) = 1 - ", Superscript["u", 2], ": ", eqMU}]}, Alignment -> Center], ImageSize -> 420]
```

The sum never drops below one bit. Equality needs $x^2+z^2=1$ and $x,z\in\{0,\pm1\}$, which are the four marked eigenstates of $X$ and $Z$: a state sharp in one basis is completely uncertain in the other. Every mixed state, and every state with $y\ne0$, has $x^2+z^2<1$ and sits strictly above the bound.
