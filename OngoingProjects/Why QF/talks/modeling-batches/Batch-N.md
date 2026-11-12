```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A qutrit beats every noncontextual model

Five projectors $\Pi_k=|v_k\rangle\langle v_k|$ on a qutrit form a pentagon: each can be measured jointly with its two neighbors, and no two neighbors both answer yes. If every test had a predetermined answer, whatever it is measured with, then $\sum_k\langle\Pi_k\rangle\le2$ (Klyachko, Can, Binicioğlu and Shumovsky); the qutrit is in $\rho=p\,|\psi_\alpha\rangle\langle\psi_\alpha|+(1-p)\,\mathbb 1/3$, a pure state tilted by $\alpha$ from $|2\rangle$ mixed with white noise. Which of these states break the bound?

With $|0\rangle,|1\rangle,|2\rangle$ as the $x,y,z$ axes of a real space, the five rays sit on a cone around $|2\rangle$, one turned from the next by $4\pi/5$:

$$
v_k=R_z\!\left(\tfrac{4\pi k}{5}\right)\left(\sin\vartheta,\,0,\,\cos\vartheta\right),\qquad \cos^2\vartheta=\tfrac1{\sqrt5}\;\Longleftrightarrow\; v_k\cdot v_{k+1}=0,\qquad |\psi_\alpha\rangle=\sin\alpha\,|0\rangle+\cos\alpha\,|2\rangle .
$$

The cell builds $K=\sum_k\Pi_k$ and plots $\langle K\rangle=\mathrm{Tr}\,\rho K$ against $p$ for four tilts, with the noncontextual bound dashed:

```wl
kcbsRays = Table[RotationMatrix[4 Pi k/5, {0, 0, 1}] . {Sqrt[1 - 1/Sqrt[5]], 0, 5^(-1/4)}, {k, 0, 4}];
kcbsK = Total[Outer[Times, #, #] & /@ kcbsRays];
kcbsMean[p_, α_] := Tr[kcbsK . (p Outer[Times, {Sin[α], 0, Cos[α]}, {Sin[α], 0, Cos[α]}] + (1 - p) IdentityMatrix[3]/3)];
kcbsαs = {0, Pi/8, Pi/4, Pi/2};
Plot[Evaluate[kcbsMean[p, #] & /@ kcbsαs], {p, 0, 1}, PlotRange -> All, Frame -> True, GridLines -> {None, {2}}, GridLinesStyle -> Dashed,
  FrameLabel -> {"p", "\[LeftAngleBracket]K\[RightAngleBracket]"}, ImageSize -> 400, PlotLegends -> (Row[{"\[Alpha] = ", #}] & /@ kcbsαs)]
```

Only states close to the cone axis, with little noise, rise above the dashed line, and of the four tilts the axis state does best. Tilting the state toward the plane perpendicular to the cone axis or mixing in noise pulls $\langle K\rangle$ down, and a state in that plane stays below 2 even when pure.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Qutrit magic in the discrete Wigner function: Clifford gates only move it

A qutrit state maps to a quasi-probability $W(q,p)$ on the nine points of a $3\times3$ phase space. For odd dimension a pure state has $W\ge0$ exactly when it is a stabilizer state (Gross), so the sum-negativity $\mathrm{sn}=\sum_{W<0}|W|$ measures magic. Starting from the stabilizer state $|+\rangle$ and from $(|1\rangle-|2\rangle)/\sqrt2$, which of the Clifford gates $F$, $XZ$ and the non-Clifford gate $T$ changes the negativity?

With $\omega=e^{2\pi i/3}$, $\zeta=e^{2\pi i/9}$ and indices mod 3, the gates are

$$
X|j\rangle=|j+1\rangle,\quad Z|j\rangle=\omega^{-j}|j\rangle,\quad F|j\rangle=\tfrac1{\sqrt3}\sum_k\omega^{jk}|k\rangle,\quad T=\mathrm{diag}\!\left(1,\zeta,\zeta^{-1}\right),\qquad |+\rangle=F|0\rangle .
$$

QF applies each gate and reads $W$ from the state's `"PhaseSpace"` property; each panel shows the nine values of $W$, negative ones in bold and framed, with the sum-negativity on top:

```wl
magicBasis = QuantumBasis[{3}]; magicFmt[x_] := ToString[NumberForm[x, {3, 2}]];
magicStates = <|"|+\[RightAngleBracket]" -> QuantumOperator["H"[3]][QuantumState[{1, 0, 0}, magicBasis]],
    "(|1\[RightAngleBracket] - |2\[RightAngleBracket])/\[Sqrt]2" -> QuantumState[{0, 1, -1}/Sqrt[2], magicBasis]|>;
magicGates = <|"no gate" -> QuantumOperator[IdentityMatrix[3], {1}, magicBasis], "F" -> QuantumOperator["H"[3]],
    "XZ" -> QuantumOperator["X"[3]] @ QuantumOperator["Z"[3]], "T" -> QuantumOperator[DiagonalMatrix[{1, Exp[2 Pi I/9], Exp[-2 Pi I/9]}], {1}, magicBasis]|>;
qutritMagicPanel[gate_, state_] := With[{w = Re[N[Normal[magicGates[gate][magicStates[state]]["PhaseSpace"]]]]}, Style[
   Labeled[Grid[Append[Reverse[MapThread[Prepend, {Map[If[# < -10^-9, Item[Style[magicFmt[#], Bold], Frame -> True, FrameStyle -> Thick], magicFmt[#]] &, w, {2}], Range[0, 2]}]],
      {"p\\q", 0, 1, 2}], Dividers -> {{2 -> True}, {-2 -> True}}, ItemSize -> {3.5, 1.6}, Alignment -> Center],
    Column[{Row[{gate, ",  ", state}], Row[{"sn = ", magicFmt[Total[Select[Flatten[-w], Positive]]]}]}, Center], Top, LabelStyle -> {FontFamily -> "Helvetica", 13}],
   BaseStyle -> {FontFamily -> "Helvetica", 13}]];
Grid[Outer[qutritMagicPanel, Keys[magicGates], Keys[magicStates]]\[Transpose], Spacings -> {3, 2}]
```

The two Clifford gates only rearrange the nine numbers, so the sum-negativity never changes: $|+\rangle$ stays nonnegative, and the single negative entry of $(|1\rangle-|2\rangle)/\sqrt2$ moves to another phase-space point under $XZ$. $T$ changes the numbers themselves: it gives $|+\rangle$ three negative entries, creating magic from a stabilizer state, and it changes the negativity of the other state too.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Search on one qudit with no entangling gate

A single $d$-level system, basis $|0\rangle,\dots,|d-1\rangle$, hides one marked level $|w\rangle$. An oracle flips the sign of $|w\rangle$, and a diffusion step reflects the state about the uniform superposition $|s\rangle$. How many oracle calls does it take to find $w$?

Starting from $|s\rangle=F|0\rangle$, with $F$ the $d$-point Fourier transform, each call applies

$$
G=D\,O,\qquad O=\mathbb 1-2|w\rangle\langle w|,\qquad D=2|s\rangle\langle s|-\mathbb 1,\qquad P_k=\left|\langle w|G^k|s\rangle\right|^2 ,
$$

and $G$ rotates the state in the plane of $|s\rangle$ and $|w\rangle$ by $2\arcsin(1/\sqrt d)\approx2/\sqrt d$, so reaching $|w\rangle$ takes about $(\pi/4)\sqrt d$ calls (dashed).

QF runs the circuit with $w=2$ for four qudit sizes and plots $P_k$ against $k/\sqrt d$:

```wl
quditSearch[d_] := With[{basis = QuantumBasis[{d}], uniform = ConstantArray[1, d]/Sqrt[d]},
   With[{step = QuantumOperator[2 KroneckerProduct[uniform, uniform] - IdentityMatrix[d], {1}, basis] @
       QuantumOperator[DiagonalMatrix[ReplacePart[ConstantArray[1, d], 3 -> -1]], {1}, basis]},
    MapIndexed[{(First[#2] - 1)/Sqrt[d], #1["ProbabilitiesList"][[3]]} &,
     NestList[step, QuantumOperator["Fourier"[d]][QuantumState[UnitVector[d, 1], basis]], Floor[1.6 Sqrt[d]]]]]];
quditSizes = {4, 16, 64, 256};
ListLinePlot[quditSearch /@ quditSizes, PlotMarkers -> Automatic, PlotRange -> {{-0.05, 1.65}, {-0.03, 1.04}}, Frame -> True, GridLines -> {{Pi/4}, None},
 GridLinesStyle -> Dashed, FrameLabel -> {Row[{"k / ", Sqrt["d"]}], "P(w)"}, ImageSize -> 420, PlotLegends -> (Row[{"d = ", #}] & /@ quditSizes)]
```

Every gate here acts on one $d$-level particle, so there is no second subsystem and no entanglement at any step, yet the marked level is found after a number of calls that grows as $\sqrt d$: in the rescaled time the curves close in on one as $d$ grows, peaking near the dashed line. The smallest qudit sits off the common curve because it starts a third of the way to $|w\rangle$ and turns by $\pi/3$ per call, so one call lands it on $|w\rangle$ with certainty.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The spin-1 NV center in a tilted magnetic field

The ground state of a nitrogen-vacancy center in diamond is a spin $S=1$. The crystal field lifts $m_s=\pm1$ above $m_s=0$ by the zero-field splitting $D$, strain adds a transverse anisotropy $E$, and a magnetic field $B$, tilted by $\vartheta$ from the defect axis $z$, couples through the electron's gyromagnetic ratio $\gamma$. Which level crossings survive strain and a tilted field?

With $\hbar=1$ and $S_x,S_y,S_z$ the spin-1 matrices,

$$
H=D\,S_z^2+E\left(S_x^2-S_y^2\right)+\gamma B\left(\cos\vartheta\,S_z+\sin\vartheta\,S_x\right).
$$

QF builds $H$ from its spin-1 operators with $E$, $B$ and $\vartheta$ symbolic and $D=1$, and the cell plots the three energy levels against the field, in units of $D$:

```wl
{nvSx, nvSy, nvSz} = QuantumOperator /@ {"JX"[1], "JY"[1], "JZ"[1]};
nvMatrix = Normal[QuantumOperator[nvSz @ nvSz + nvE (nvSx @ nvSx - nvSy @ nvSy) + nvB (Cos[nvTheta] nvSz + Sin[nvTheta] nvSx), QuantumBasis[3]]["Matrix"]];
nvLevels[b_?NumericQ, e_, th_] := Sort[Eigenvalues[N[nvMatrix /. {nvB -> b, nvE -> e, nvTheta -> th}]]];
nvPanel[e_, th_] := Plot[nvLevels[b, e, th], {b, 0, 2}, PlotRange -> All, Frame -> True, Axes -> False, FrameLabel -> {"\[Gamma]B / D", "energy / D"},
   PlotLabel -> Row[{"E/D = ", e, ",   \[CurlyTheta] = ", th}], ImageSize -> 260];
Grid[Outer[nvPanel, {0, 0.1}, {0, 5 Degree, 30 Degree}]]
```

Along the axis the field conserves $m_s$, so the $m_s=-1$ level falls straight through $m_s=0$ near $\gamma B=D$, and strain, which couples only $m_s=+1$ to $m_s=-1$, splits those two at zero field but leaves that crossing intact. A tilted field mixes $m_s=0$ with $m_s=\pm1$ and turns the crossing into an avoided one, with a gap that widens with the tilt.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A five-site Heisenberg ring split into symmetry sectors

Five spins $\tfrac12$ sit on a ring with antiferromagnetic exchange $J>0$ between neighbors. The Hamiltonian commutes with the total magnetization $M=\sum_jS^z_j$ and with the translation $T$ that moves every spin one site along the ring, whose eigenvalues are $e^{ik}$ with $k=2\pi m/5$, $m=-2,\dots,2$. Which magnetization and momentum does the ground state carry?

With $\mathbf S_j=\boldsymbol\sigma_j/2$, $|0\rangle$ as spin up, and site 6 identified with site 1,

$$
H=J\sum_{j=1}^{5}\mathbf S_j\cdot\mathbf S_{j+1},\qquad [H,M]=[H,T]=[M,T]=0 ;
$$

flipping every spin commutes with $H$ and sends $M\to-M$, so the sectors with $M<0$ repeat those with $M>0$.

QF builds $H$, $M$ and $T$ on the 32-dimensional space; the cell projects onto each joint eigenspace of $M$ and $T$ and plots $|H|$ in the product basis, $|H|$ in the basis sorted by $(M,k)$, and the energies of each sector, in units of $J$:

```wl
ringPairs = Table[{j, Mod[j, 5] + 1}, {j, 5}];
ringH = Normal[(Total[QuantumOperator @@@ Tuples[{{"XX", "YY", "ZZ"}, ringPairs}]]/4)["Matrix"]];
ringT = Normal[QuantumOperator["Permutation"[ConstantArray[2, 5], Cycles[{Range[5]}]]]["Matrix"]];
ringM = Diagonal[Normal[(Total[QuantumOperator["Z", {#}] & /@ Range[5]]/2)["Matrix"]]];
ringBasis[mm_, m_] := With[{p = N[(Total[Exp[-2 Pi I m #/5] MatrixPower[ringT, #] & /@ Range[0, 4]]/5) . DiagonalMatrix[Boole[Thread[ringM == mm]]]]},
   With[{es = Eigensystem[(p + ConjugateTranspose[p])/2]}, Pick[es[[2]], Thread[Re[es[[1]]] > 1/2]]]];
ringEnergies[mm_, m_] := With[{v = ringBasis[mm, m]}, If[v === {}, {}, Re[Eigenvalues[Conjugate[v] . ringH . Transpose[v]]]]];
ringV = Join @@ (ringBasis @@@ Tuples[{{5/2, 3/2, 1/2, -1/2, -3/2, -5/2}, Range[-2, 2]}]);
ringMs = {1/2, 3/2, 5/2};
Row[{ArrayPlot[Abs[ringH], Frame -> True, Axes -> False, FrameTicks -> None, PlotLabel -> "product basis", ImageSize -> 220],
  ArrayPlot[Abs[Chop[Conjugate[ringV] . ringH . Transpose[ringV]]], Frame -> True, Axes -> False, FrameTicks -> None, PlotLabel -> "sorted by (M, k)", ImageSize -> 220],
  ListPlot[Table[Flatten[Table[{m, en}, {m, -2, 2}, {en, ringEnergies[mm, m]}], 1], {mm, ringMs}],
   PlotMarkers -> {{"\[EmptyCircle]", 22}, {"\[EmptySquare]", 14}, {"\[FilledCircle]", 7}}, Frame -> True, PlotRange -> {{-2.5, 2.5}, All}, Axes -> False,
   FrameLabel -> {"5k / 2\[Pi]", "E / J"}, ImageSize -> 300, PlotLegends -> (Row[{"M = ", #}] & /@ ringMs)]}, "   "]
```

The 32 levels fall into blocks no larger than $2\times2$, one for each pair $(M,k)$. The ground state sits at $M=\pm\tfrac12$ and $k=\pm2\pi/5$, four degenerate states. Every level of a larger $M$ reappears at the same $k$ in each smaller $M$, as the total-spin multiplets of a rotation-invariant $H$ require, and the energies at $k$ and $-k$ agree.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A variational upper bound on the ground-state energy

Six spins on a ring have a ferromagnetic Ising coupling $J$ along $z$ and sit in a transverse field $h$ along $x$. The trial state tilts every spin by the same angle $\theta$ from $|0\rangle$ toward $|1\rangle$, a product state with no entanglement. How close does the best such state come to the true ground-state energy $E_0$?

With $Z_7\equiv Z_1$, the Hamiltonian, the trial state and the Rayleigh-Ritz bound are

$$
H=-J\sum_{j=1}^{6}Z_jZ_{j+1}-h\sum_{j=1}^{6}X_j,\qquad |\theta\rangle=\left(\cos\tfrac\theta2\,|0\rangle+\sin\tfrac\theta2\,|1\rangle\right)^{\otimes6},\qquad E(\theta)=\langle\theta|H|\theta\rangle\ge E_0 .
$$

QF evaluates $E(\theta)$ in closed form with $\theta$, $J$ and $h$ symbolic, finds $E_0$ by diagonalizing $H$ (dashed), and plots both per spin in units of $J$:

```wl
isingH = -isingJ Total[QuantumOperator["ZZ", {#, Mod[#, 6] + 1}] & /@ Range[6]] - hx Total[QuantumOperator["X", {#}] & /@ Range[6]];
isingTrial = QuantumTensorProduct @@ ConstantArray[QuantumState[{Cos[isingTheta/2], Sin[isingTheta/2]}], 6];
isingE = FullSimplify[(isingTrial["Dagger"] @ isingH @ isingTrial)["Scalar"]/6, isingTheta ∈ Reals];
isingE0[field_] := Min[Eigenvalues[N[Normal[isingH["Matrix"]] /. {isingJ -> 1, hx -> field}]]]/6;
isingFields = {0.5, 1, 2, 4};
Show[Plot[Evaluate[isingE /. {isingJ -> 1, hx -> #} & /@ isingFields], {isingTheta, 0, Pi}, PlotRange -> All, PlotLegends -> (Row[{"h/J = ", #}] & /@ isingFields)],
 Plot[Evaluate[isingE0 /@ isingFields], {isingTheta, 0, Pi}, PlotStyle -> Dashed, PlotRange -> All],
 PlotRange -> {{0, Pi}, {-4.2, -0.4}}, Frame -> True, Axes -> False, FrameLabel -> {"\[Theta]", "energy per spin / J"}, ImageSize -> 420]
```

Each trial curve stays above its dashed line, as the bound requires, and its minimum gives the best product-state estimate. The gap at the minimum, the energy a product state cannot reach, is tiny at $h/J=0.5$, grows up to $h/J=2$, where the best product state points fully along $x$, and shrinks again at $h/J=4$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Fermi's golden rule from a level coupled to a band

A discrete level $|i\rangle$ at energy zero couples with the same matrix element $V$ to each level $|k\rangle$ of a band, levels spaced by $\delta$ from $-W/2$ to $W/2$. First-order perturbation theory predicts that $|i\rangle$ empties at the golden-rule rate $\Gamma$. When does the exact population of $|i\rangle$ follow $e^{-\Gamma t}$?

With $\hbar=1$ and $E_k=k\delta$ for $k=-K,\dots,K$, $W=2K\delta$,

$$
H=\sum_kE_k|k\rangle\langle k|+V\sum_k\left(|i\rangle\langle k|+|k\rangle\langle i|\right),\qquad 1-P_i(t)\approx\sum_k\frac{4V^2\sin^2(E_kt/2)}{E_k^2}\;\longrightarrow\;\Gamma t,\qquad \Gamma=\frac{2\pi V^2}{\delta},
$$

where the limit treats the band as a continuum of density $1/\delta$.

QF solves the Schrödinger equation from $|i\rangle$ for $P_i(t)$, with $V$ set by $\Gamma$, and compares it with $e^{-\Gamma t}$ (dashed); time is in units of $1/\Gamma$:

```wl
fgrP[gd_, wg_] := With[{kk = Round[wg gd/2], v = Sqrt[1/(2 Pi gd)]}, With[{basis = QuantumBasis[{2 kk + 2}]},
    QuantumEvolve[QuantumOperator[SparseArray[{{1, j_} /; j > 1 -> v, {j_, 1} /; j > 1 -> v, {j_, j_} /; j > 1 -> (j - kk - 2)/gd}, {2 kk + 2, 2 kk + 2}], {1}, basis],
      QuantumState[UnitVector[2 kk + 2, 1], basis], {t, 0, 15}]["ProbabilitiesList"][[1]]]];
fgrPanel[wg_, gd_] := Show[Plot[Evaluate[fgrP[gd, wg]], {t, 0, 15}, PlotRange -> All], Plot[Exp[-t], {t, 0, 15}, PlotStyle -> Dashed, PlotRange -> All],
   PlotRange -> {0, 1.02}, Frame -> True, Axes -> False, FrameLabel -> {"\[CapitalGamma] t", Subscript["P", "i"]},
   PlotLabel -> Row[{"\[CapitalGamma]/\[Delta] = ", gd, ",   W/\[CapitalGamma] = ", wg}], ImageSize -> 280];
Grid[Outer[fgrPanel, {4, 40}, {0.5, 2, 8}]]
```

In the wide band the population follows $e^{-\Gamma t}$ until $t=2\pi/\delta$, when the levels come back into phase and return part of it; the smaller $\Gamma/\delta$, the earlier and larger this revival, and at $\Gamma/\delta=8$ it comes after the plotted time. The narrow band breaks the rule at short times: the population starts flat, stays above $e^{-\Gamma t}$, then falls below it.
