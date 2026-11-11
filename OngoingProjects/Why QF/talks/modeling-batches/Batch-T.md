```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Adding two momenta on a noncommuting spacetime

A plane wave $e^{i\alpha_\mu x^\mu}$ carries momentum $\alpha_\mu$, and when the coordinates $x^\mu$ do not commute, the product of two plane waves defines a deformed sum $\alpha\oplus\beta$. Snyder's spacetime has a fundamental length $a$, and $\kappa$-Minkowski spacetime has a mass scale $\kappa$. Does the deformed sum depend on which momentum comes first?

With $\eta=\operatorname{diag}(1,-1,-1,-1)$ and Lorentz generators $M^{\mu\nu}$, the two algebras are

$$
\text{Snyder:}\ \ [x^\mu,x^\nu]=ia^2M^{\mu\nu},\quad [M^{\mu\nu},x^\rho]=i\left(\eta^{\nu\rho}x^\mu-\eta^{\mu\rho}x^\nu\right);\qquad \kappa\text{-Minkowski:}\ \ [x^0,x^j]=\tfrac{i}{\kappa}\,x^j,\quad [x^j,x^k]=0,
$$

with the Lorentz algebra for $[M,M]$. Writing $e^{X}e^{Y}=e^{Z}$, the Baker-Campbell-Hausdorff series to fourth order is

$$
Z=X+Y+\tfrac12[X,Y]+\tfrac1{12}\big([X,[X,Y]]+[Y,[Y,X]]\big)-\tfrac1{24}[Y,[X,[X,Y]]],\qquad X=i\alpha_\mu x^\mu,\ Y=i\beta_\mu x^\mu,
$$

and $(\alpha\oplus\beta)_\mu$ is the coefficient of $ix^\mu$ in $Z$. In Snyder's algebra $Z$ also carries Lorentz generators, starting with $\tfrac12[X,Y]=-\tfrac{i}{2}a^2\alpha_\mu\beta_\nu M^{\mu\nu}$; they are not part of the momentum.

`NonCommutativeAlgebra` holds each set of relations and `NonCommutativeExpand` reduces $Z$ with $\alpha$, $\beta$ symbolic; the figure shows two momenta in the $(p_0,p_1)$ plane, in units of $1/a$ for Snyder and of $\kappa$ for $\kappa$-Minkowski, for $\alpha=s\,(1,\tfrac12,0,0)$ and $\beta=s\,(\tfrac35,-\tfrac45,0,0)$ at three scales $s$:

```wl
composeEta[μ_, ν_] := DiagonalMatrix[{1, -1, -1, -1}][[μ + 1, ν + 1]];
composeL[μ_, ν_] := Which[μ == ν, 0, μ < ν, composeM[μ, ν], True, -composeM[ν, μ]];
composePairs = Subsets[Range[0, 3], {2}];
composeX = composeQ /@ Range[0, 3];
composeSnyder = NonCommutativeAlgebra[<|"Generators" -> {Join[composeX, composeM @@@ composePairs]}, "CommutationRelations" -> Join[
     Commutator[composeQ[#2], composeQ[#1]] + I composeA^2 composeM[#1, #2] & @@@ composePairs,
     Flatten[Table[Commutator[composeM[μ, ν], composeQ[ρ]] - I (composeEta[ν, ρ] composeQ[μ] - composeEta[μ, ρ] composeQ[ν]), {μ, 0, 3}, {ν, μ + 1, 3}, {ρ, 0, 3}]],
     Replace[Subsets[composePairs, {2}], {{μ_, ν_}, {ρ_, σ_}} :> Commutator[composeM[ρ, σ], composeM[μ, ν]] +
        I (composeEta[ν, ρ] composeL[μ, σ] - composeEta[ν, σ] composeL[μ, ρ] - composeEta[μ, ρ] composeL[ν, σ] + composeEta[μ, σ] composeL[ν, ρ]), {1}]]|>];
composeKappa = NonCommutativeAlgebra[<|"Generators" -> {composeX}, "CommutationRelations" ->
     Join[Commutator[composeQ[#], composeQ[0]] + I composeQ[#]/composeK & /@ Range[3], Commutator[composeQ[#2], composeQ[#1]] & @@@ Subsets[Range[3], {2}]]|>];
composeBCH[x_, y_] := x + y + Commutator[x, y]/2 + (Commutator[x, Commutator[x, y]] + Commutator[y, Commutator[y, x]])/12 - Commutator[y, Commutator[x, Commutator[x, y]]]/24;
{composeAl, composeBe} = {composeα /@ Range[0, 3], composeβ /@ Range[0, 3]};
composeSum[alg_] := composeSum[alg] = Simplify[Coefficient[NonCommutativeExpand[composeBCH[I composeAl . composeX, I composeBe . composeX], alg], #]/I & /@ composeX];
composeSwap = {composeα -> composeβ, composeβ -> composeα};
composeShow = {composeα[μ_] :> Subscript["α", μ], composeβ[μ_] :> Subscript["β", μ], composeK -> "κ"};
composeDifference = Factor[composeSum[#] - (composeSum[#] /. composeSwap)] & /@ {composeSnyder, composeKappa};
composePanel[alg_, unit_, name_, s_] := With[{pts = Thread[Join[composeAl, composeBe] -> s Join[{1, 1/2, 0, 0}, {3/5, -4/5, 0, 0}]]},
   With[{α = composeAl[[;; 2]] /. pts, β = composeBe[[;; 2]] /. pts},
    ListPlot[{{{0, 0}, α}, {{0, 0}, β}, {α + β}, {composeSum[alg][[;; 2]] /. pts /. unit}, {(composeSum[alg] /. composeSwap)[[;; 2]] /. pts /. unit}},
     Joined -> {True, True, False, False, False}, PlotMarkers -> {None, None, {"\[FilledSquare]", 14}, {"\[FilledCircle]", 10}, {"\[EmptyCircle]", 24}},
     PlotStyle -> {Automatic, Dashed, Automatic, Automatic, Automatic}, PlotLegends -> {"α", "β", "α + β", "α ⊕ β", "β ⊕ α"},
     PlotRange -> s {{-0.1, 1.9}, {-1.1, 0.7}}, AspectRatio -> Automatic, Frame -> True, FrameLabel -> {Subscript["p", 0], Subscript["p", 1]},
     PlotLabel -> Row[{name, ",   s = ", s}], ImageSize -> 260]]];
composePanels = Table[composePanel[Sequence @@ m, s], {m, {{composeSnyder, composeA -> 1, "Snyder"}, {composeKappa, composeK -> 1, "κ-Minkowski"}}}, {s, {0.3, 0.6, 0.9}}];
Labeled[Legended[Grid[Map[First, composePanels, {2}]], Last[composePanels[[1, 1]]]],
  Style[Column[{Row[{"Snyder:  α ⊕ β − β ⊕ α = ", composeDifference[[1]]}],
     Row[{"κ-Minkowski:  α ⊕ β − β ⊕ α = ", TraditionalForm[composeDifference[[2]] /. composeShow]}]}], 13, FontFamily -> "Helvetica"], Top]
```

The deformed sum moves away from $\alpha+\beta$ in both spacetimes, more so at larger $s$. In Snyder's spacetime both orders land on the same point, and the exchange difference vanishes at every order computed. In $\kappa$-Minkowski spacetime the two orders split, and the difference is a spatial term proportional to $\alpha_0\beta_j-\alpha_j\beta_0$, so composing two momenta depends on their order.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Determinant and trace rules proved at symbolic dimension

$A$ and $B$ are arbitrary complex $n\times n$ matrices, for example generators $A=-iHt$ of an $n$-level system, and $K$ is an arbitrary $m\times m$ matrix acting on a second system of $m$ levels; no entry and no value of $n$ or $m$ is given. Which of the textbook determinant and trace rules can be proved for every $n$ and $m$ at once?

The candidates are

$$
\det(AB)=\det A\,\det B,\quad \det e^{A}=e^{\operatorname{Tr}A},\quad \operatorname{Tr}(BAB^{-1})=\operatorname{Tr}A,\quad (AB)^{-1}=B^{-1}A^{-1},\quad \det(A\otimes K)=(\det A)^m(\det K)^n,\quad \operatorname{Tr}(A\otimes K)=\operatorname{Tr}A\,\operatorname{Tr}K,
$$

and two rules that are false in general, $\det(A+B)=\det A+\det B$ and $\operatorname{Tr}e^{A}=e^{\operatorname{Tr}A}$.

Each column tries to reduce the difference of the two sides to zero: `TensorReduce` with index symmetries only, `ArrayExpand` and `ArraySimplify` with the matrix-algebra identities, all at symbolic $n$ and $m$; the last column substitutes random real matrices for every $n=2,3,4$ and $m=2,3$:

```wl
matIdDims = matIdA ∈ Matrices[{n, n}] && matIdB ∈ Matrices[{n, n}] && matIdK ∈ Matrices[{m, m}];
matIdRules = {{"det(AB) = det A det B", Det[matIdA . matIdB] - Det[matIdA] Det[matIdB]},
   {"det exp(A) = exp(Tr A)", Det[MatrixExp[matIdA]] - Exp[Tr[matIdA]]},
   {Row[{"Tr(B A ", Superscript["B", -1], ") = Tr A"}], Tr[matIdB . matIdA . Inverse[matIdB]] - Tr[matIdA]},
   {Row[{Superscript["(AB)", -1], " = ", Superscript["B", -1], Superscript["A", -1]}], Inverse[matIdA . matIdB] - Inverse[matIdB] . Inverse[matIdA]},
   {Row[{"det(A ⊗ K) = ", Superscript["(det A)", "m"], Superscript["(det K)", "n"]}], Det[KroneckerProduct[matIdA, matIdK]] - Det[matIdA]^m Det[matIdK]^n},
   {"Tr(A ⊗ K) = Tr A Tr K", Tr[KroneckerProduct[matIdA, matIdK]] - Tr[matIdA] Tr[matIdK]},
   {"det(A + B) = det A + det B", Det[matIdA + matIdB] - Det[matIdA] - Det[matIdB]},
   {"Tr exp(A) = exp(Tr A)", Tr[MatrixExp[matIdA]] - Exp[Tr[matIdA]]}};
matIdZeroQ = If[MatchQ[#, 0 | _SymbolicZerosArray], "0", "not 0"] &;
matIdNumericQ[expr_] := BlockRandom[SeedRandom[3]; AllTrue[Tuples[{{2, 3, 4}, {2, 3}}],
     Max[Abs[expr /. {n -> #[[1]], m -> #[[2]], matIdA -> RandomReal[{-1, 1}, {#[[1]], #[[1]]}], matIdB -> RandomReal[{-1, 1}, {#[[1]], #[[1]]}],
          matIdK -> RandomReal[{-1, 1}, {#[[2]], #[[2]]}]}]] < 10^-9 &]];
Grid[Prepend[{#[[1]], matIdZeroQ[TensorReduce[#[[2]], Assumptions -> matIdDims]], matIdZeroQ[ArrayExpand[#[[2]], matIdDims]],
      matIdZeroQ[ArraySimplify[#[[2]], matIdDims]], If[matIdNumericQ[#[[2]]], "holds", "fails"]} & /@ matIdRules,
   {"rule", "TensorReduce", "ArrayExpand", "ArraySimplify", "random matrices"}],
  Frame -> All, Alignment -> {{Left, Center, Center, Center, Center}}, Spacings -> {2, 1}, BaseStyle -> {14, FontFamily -> "Helvetica"}]
```

Index symmetries alone settle only the trace of a Kronecker product, which is a pure index contraction. The matrix-algebra identities inside `ArrayExpand` and `ArraySimplify` prove all six true rules for every $n$ and $m$, and leave the two false rules standing, as the random matrices confirm.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Energies of all-to-all exchange from Young diagrams

$n$ identical sites, each with $d$ levels, are coupled pairwise by the exchange $P_{ij}$, which swaps the states of sites $i$ and $j$, with one coupling $J$ for every pair; for $d=2$, $P_{ij}=\tfrac12(1+\vec\sigma_i\cdot\vec\sigma_j)$ is the Heisenberg exchange. What are the energies and their degeneracies for every $n$ and $d$?

$$
H=J\sum_{i<j}P_{ij}\quad\text{on}\quad(\mathbb C^d)^{\otimes n}=\bigoplus_{\lambda\vdash n}V_\lambda\otimes W_\lambda(d),
$$

where the sum runs over Young diagrams $\lambda$ with $n$ boxes (Schur-Weyl duality): $V_\lambda$ carries the permutations of the sites, with $\dim V_\lambda=f^\lambda$ standard tableaux, and $W_\lambda(d)$ the common unitaries $U^{\otimes n}$, with $\dim W_\lambda(d)=s_\lambda(d)$, a polynomial in $d$ that vanishes when $\lambda$ has more than $d$ rows. The sum of all transpositions acts on block $\lambda$ as its content sum, so

$$
E_\lambda=J\,c(\lambda),\qquad c(\lambda)=\sum_{(r,k)\in\lambda}(k-r),\qquad \text{degeneracy } f^\lambda s_\lambda(d),
$$

with $(r,k)$ the row and column of each box.

The TensorNetworks functions `TableauDimension` and `SchurDimension` give $f^\lambda$ and $s_\lambda(d)$; the figure plots each predicted level with its degeneracy against the eigenvalues of the $d^n\times d^n$ matrix of $H$, energies in units of $J$, and the top label gives the $n=4$ levels for symbolic $d$:

```wl
Needs["Wolfram`TensorNetworks`Symmetry`"];
permSumContent[λ_] := Total[Flatten[Table[k - r, {r, Length[λ]}, {k, λ[[r]]}]]];
permSumLevels[nn_, dd_] := Select[{permSumContent[#], TableauDimension[#] SchurDimension[#, dd]} & /@ IntegerPartitions[nn], Last[#] > 0 &];
permSumSwap[nn_, dd_, {i_, j_}] := SparseArray[MapIndexed[{First[#2], FromDigits[Permute[#1, Cycles[{{i, j}}]] - 1, dd] + 1} -> 1 &, Tuples[Range[dd], nn]], {dd^nn, dd^nn}];
permSumNumeric[nn_, dd_] := Tally[Round[Eigenvalues[N[Normal[Total[permSumSwap[nn, dd, #] & /@ Subsets[Range[nn], {2}]]]]]]];
permSumPanel[dd_, nn_] := ListPlot[{permSumLevels[nn, dd], permSumNumeric[nn, dd]}, Filling -> {1 -> Axis}, PlotMarkers -> {{"\[FilledCircle]", 9}, {"\[EmptyCircle]", 20}},
   PlotLegends -> {Row[{"J c(λ), ", Superscript["f", "λ"], Subscript["s", "λ"], "(d)"}], "eigenvalues of H"}, PlotRange -> {{-nn (nn - 1)/2 - 1, nn (nn - 1)/2 + 1}, {-0.06, 1.1} Max[permSumLevels[nn, dd][[All, 2]]]},
   Frame -> True, FrameLabel -> {"E / J", "degeneracy"}, PlotLabel -> Row[{"d = ", dd, ",   n = ", nn}], ImageSize -> 250];
permSumPanels = Outer[permSumPanel, {2, 3, 4}, {3, 4, 5}];
Labeled[Legended[Grid[Map[First, permSumPanels, {2}]], Last[permSumPanels[[1, 1]]]],
  Style[Column[Row[{"λ = ", #, ":   E = ", permSumContent[#], " J,   degeneracy ", TraditionalForm[Factor[TableauDimension[#] SchurDimension[#, d]]]}] & /@ IntegerPartitions[4]], 13, FontFamily -> "Helvetica"], Top]
```

Every eigenvalue of $H$ sits on a predicted level with the predicted degeneracy, in every panel. A level exists only when its diagram fits in $d$ rows, so the most negative energies appear one at a time as $d$ grows, the fully antisymmetric level $-n(n-1)/2$ only once $d\ge n$. In the top label the degeneracy of a diagram with $r$ rows carries the factors $(d-1)\cdots(d-r+1)$, so it vanishes for every $d<r$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Fidelity of two mixed qubits in closed form

Two qubit states $\rho_1=\tfrac12(I+\vec r_1\cdot\vec\sigma)$ and $\rho_2=\tfrac12(I+\vec r_2\cdot\vec\sigma)$ have Bloch vectors of lengths $r_1,r_2\le1$ at an angle $\theta$. Uhlmann's fidelity is defined through nested matrix square roots; does it reduce to a formula in the Bloch vectors alone?

$$
F(\rho_1,\rho_2)=\Big(\operatorname{Tr}\sqrt{\sqrt{\rho_1}\,\rho_2\sqrt{\rho_1}}\Big)^2\;\overset{?}{=}\;\tfrac12\Big(1+\vec r_1\cdot\vec r_2+\sqrt{(1-r_1^2)(1-r_2^2)}\Big).
$$

$F$ is unchanged when both states are rotated together, so $\vec r_1=r_1\hat z$ and $\vec r_2=r_2(\sin\theta,0,\cos\theta)$ cover every pair.

`FullSimplify` compares the definition, square roots included, with the formula for symbolic $r_1$, $r_2$, $\theta$; the figure draws the formula against $\theta$ and, as markers, the definition evaluated numerically on random pairs of states with the same lengths and angle, each pair turned by a random rotation:

```wl
uhlBloch[r_] := (IdentityMatrix[2] + r . PauliMatrix[{1, 2, 3}])/2;
uhlDefinition[r1_, r2_] := With[{s = MatrixPower[uhlBloch[r1], 1/2]}, Tr[MatrixPower[s . uhlBloch[r2] . s, 1/2]]^2];
uhlFormula[r1_, r2_] := (1 + r1 . r2 + Sqrt[(1 - r1 . r1) (1 - r2 . r2)])/2;
uhlPair[l1_, l2_, θ_] := {{0, 0, l1}, l2 {Sin[θ], 0, Cos[θ]}};
uhlIdentity = FullSimplify[uhlDefinition @@ uhlPair[l1, l2, θ] == uhlFormula @@ uhlPair[l1, l2, θ], 0 <= l1 <= 1 && 0 <= l2 <= 1 && 0 <= θ <= π];
uhlRandom[l1_, l2_, θ_] := With[{rot = RotationMatrix[RandomReal[{0, 2 π}], RandomPoint[Sphere[]]]}, Re[uhlDefinition @@ (rot . # & /@ uhlPair[l1, l2, θ])]];
uhlPanel[l1_] := Show[
   Plot[Evaluate[uhlFormula @@ uhlPair[l1, #, θ] & /@ {0, 0.5, 1}], {θ, 0, π}, PlotRange -> All,
    PlotLegends -> (Row[{Subscript["r", 2], " = ", #}] & /@ {0, 0.5, 1}), PlotStyle -> {Automatic, Automatic, Dashed}],
   BlockRandom[SeedRandom[5]; ListPlot[Table[{θ, uhlRandom[l1, #, θ]}, {θ, 0, π, π/8}] & /@ {0, 0.5, 1}, PlotMarkers -> {"\[EmptyCircle]", 14}]],
   PlotRange -> {{0, π}, {-0.02, 1.02}}, Axes -> False, Frame -> True, FrameLabel -> {"angle θ between the Bloch vectors", "F"},
   FrameTicks -> {{Automatic, None}, {Table[{x, x}, {x, 0, π, π/4}], None}}, PlotLabel -> Row[{Subscript["r", 1], " = ", l1}], ImageSize -> 300];
uhlPanels = Map[uhlPanel, {{0, 0.5}, {0.9, 1}}, {2}];
Labeled[Legended[Grid[Map[First, uhlPanels, {2}]], Last[uhlPanels[[1, 1]]]],
  Style[Row[{"definition = formula for all ", Subscript["r", 1], ", ", Subscript["r", 2], ", θ:  ", uhlIdentity}], 13, FontFamily -> "Helvetica"], Top]
```

The nested square roots collapse exactly: the definition equals the formula for every pair of lengths and every angle, and the randomly rotated pairs fall on the curves. $F$ reaches zero, at $\theta=\pi$, only when both states are pure; any mixedness keeps opposite Bloch vectors from being perfectly distinguishable. The maximally mixed state ($r_1=0$) has the same fidelity with every state of a given length $r_2$, whatever its direction.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Two qutrits entangled with a positive partial transpose

Five orthonormal product states of two qutrits, the Tiles basis, leave a four-dimensional subspace, and $\rho$ is the normalized projector onto it; $\rho_p=(1-p)\,\rho+p\,I/9$ adds white noise of weight $p$. Its partial transpose has no negative eigenvalue, so the Peres test is silent. Is $\rho$ nevertheless entangled?

With qutrit states $|0\rangle,|1\rangle,|2\rangle$, the Tiles states and the state are

$$
|0\rangle|0{-}1\rangle,\ \ |2\rangle|1{-}2\rangle,\ \ |0{-}1\rangle|2\rangle,\ \ |1{-}2\rangle|0\rangle,\ \ |0{+}1{+}2\rangle|0{+}1{+}2\rangle,\qquad \rho=\tfrac14\Big(I-\sum_{k=1}^{5}|\psi_k\rangle\langle\psi_k|\Big),
$$

with $|i{\pm}j\rangle$ normalized. A separable $\rho$ would be a mixture of product states in its range, that is, product states orthogonal to all five $|\psi_k\rangle$. The realignment test $\|R(\rho)\|_1>1$, with $R(\rho)_{ij,kl}=\rho_{ik,jl}$, is a second, independent proof of entanglement.

`Reduce` asks for a product vector $|a\rangle|b\rangle$ with complex amplitudes, $|a\rangle,|b\rangle\neq0$, orthogonal to all five; the figure plots against $p$ the smallest eigenvalue of $\rho_p^{T_B}$, exactly, and the realignment value $\|R(\rho_p)\|_1-1$:

```wl
tilesKet = UnitVector[3, # + 1] &;
tilesStates = Normalize[Flatten[KroneckerProduct @@ #]] & /@ {{tilesKet[0], tilesKet[0] - tilesKet[1]}, {tilesKet[2], tilesKet[1] - tilesKet[2]},
    {tilesKet[0] - tilesKet[1], tilesKet[2]}, {tilesKet[1] - tilesKet[2], tilesKet[0]}, {Total[tilesKet /@ Range[0, 2]], Total[tilesKet /@ Range[0, 2]]}};
tilesRho = (IdentityMatrix[9] - Total[KroneckerProduct[#, #] & /@ tilesStates])/4;
tilesNoisy[p_] := (1 - p) tilesRho + p IdentityMatrix[9]/9;
tilesProductQ = With[{a = Array[tilesA, 3], b = Array[tilesB, 3]},
   Reduce[And @@ Thread[tilesStates . Flatten[KroneckerProduct[a, b]] == 0] && Or @@ Thread[a != 0] && Or @@ Thread[b != 0], Join[a, b]]];
tilesRegroup[r_, perm_] := ArrayReshape[Transpose[ArrayReshape[r, {3, 3, 3, 3}], perm], {9, 9}];
tilesMinPT = Min[Eigenvalues[tilesRegroup[tilesNoisy[p], {1, 4, 3, 2}]]];
tilesRealign[p_?NumericQ] := Total[SingularValueList[N[tilesRegroup[tilesNoisy[p], {1, 3, 2, 4}]], Tolerance -> 0]] - 1;
Labeled[Plot[{tilesMinPT, tilesRealign[p]}, {p, 0, 0.2}, PlotRange -> All, Frame -> True, GridLines -> {None, {0}},
   PlotLegends -> {Row[{"smallest eigenvalue of ", Superscript[Subscript["ρ", "p"], Subscript["T", "B"]]}], Row[{"realignment ", "\[LeftDoubleBracketingBar]", "R(", Subscript["ρ", "p"], ")", "\[RightDoubleBracketingBar]", Subscript["", 1], " - 1"}]},
   PlotStyle -> {Automatic, Dashed}, FrameLabel -> {"noise weight p", None}, ImageSize -> 460],
  Style[Row[{"product vector orthogonal to all five Tiles states exists:  ", tilesProductQ}], 13, FontFamily -> "Helvetica"], Top]
```

`Reduce` finds no product vector orthogonal to the Tiles states, so the range of $\rho$ contains no product state and $\rho$ is entangled, although the smallest eigenvalue of its partial transpose stays nonnegative at every noise weight drawn. Realignment sees the entanglement at $p=0$ and keeps seeing it under weak noise, until the dashed curve crosses zero; beyond that point neither test decides.
