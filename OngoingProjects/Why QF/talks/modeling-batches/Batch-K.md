```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Tensor identities proved at symbolic dimension

Tensors in dimension $d$ are given only by their index symmetries, never by entries: an antisymmetric $A$, a symmetric $S$, an inverse metric $g$, and a tensor $R$ with the pair symmetries of the Riemann tensor; in three dimensions, the Levi-Civita symbol $\varepsilon$ and a Hessian $H_{jk}=\partial_j\partial_k f$. Which contractions vanish for every $d$ by symmetry alone?

$$
A_{ij}=-A_{ji},\quad S_{ij}=S_{ji},\quad g^{ij}=g^{ji},\quad R_{abcd}=-R_{bacd}=-R_{abdc}=R_{cdab},\quad \varepsilon_{ijk}\ \text{totally antisymmetric},\quad H_{jk}=H_{kj},
$$

and the Ricci tensor is $R_{bd}=R_{abcd}\,g^{ac}$ (repeated indices summed).

`TensorReduce` brings each contraction to a canonical form from these declarations alone:

```wl
tensorDeclarations = {aT ∈ Arrays[{d, d}, Reals, Antisymmetric[{1, 2}]], sT ∈ Arrays[{d, d}, Reals, Symmetric[{1, 2}]],
   gT ∈ Arrays[{d, d}, Reals, Symmetric[{1, 2}]], rT ∈ Arrays[{d, d, d, d}, Reals, {{Cycles[{{1, 2}}], -1}, {Cycles[{{3, 4}}], -1}, {Cycles[{{1, 3}, {2, 4}}], 1}}],
   epsT ∈ Arrays[{3, 3, 3}, Reals, Antisymmetric[{1, 2, 3}]], hT ∈ Arrays[{3, 3}, Reals, Symmetric[{1, 2}]]};
ricciT = TensorContract[TensorProduct[rT, gT], {{1, 5}, {3, 6}}];
reduceT = TensorReduce[#, Assumptions -> tensorDeclarations] &;
reducesToZero = MatchQ[reduceT[#], 0 | _SymbolicZerosArray] &;
Grid[{{"contraction", "reduces to zero for every d"},
   {"\!\(\*SubscriptBox[\(A\), \(ij\)]\)\!\(\*SubscriptBox[\(S\), \(ij\)]\)", reducesToZero[TensorContract[TensorProduct[aT, sT], {{1, 3}, {2, 4}}]]},
   {"\!\(\*SubscriptBox[\(R\), \(bd\)]\) - \!\(\*SubscriptBox[\(R\), \(db\)]\)", reduceT[ricciT] === reduceT[TensorTranspose[ricciT, {2, 1}]]},
   {"\!\(\*SubscriptBox[\(\[CurlyEpsilon]\), \(ijk\)]\)\!\(\*SubscriptBox[\(H\), \(jk\)]\)", reducesToZero[TensorContract[TensorProduct[epsT, hT], {{2, 4}, {3, 5}}]]},
   {"\!\(\*SubscriptBox[\(A\), \(ij\)]\)\!\(\*SubscriptBox[\(A\), \(ij\)]\)", reducesToZero[TensorContract[TensorProduct[aT, aT], {{1, 3}, {2, 4}}]]}},
  Frame -> All, Alignment -> Left, Spacings -> {2, 1}]
```

Symmetry alone makes $A_{ij}S_{ij}$ vanish for every $d$, makes the Ricci tensor symmetric, and makes the curl of a gradient vanish in three dimensions. The sum of squares $A_{ij}A_{ij}$, which no symmetry forces to vanish, is correctly not reduced to zero.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### The first Bianchi identity from the Riemann tensor's Young diagram

The Riemann tensor $R_{abcd}$ in $d$ dimensions has the pair symmetries $R_{abcd}=-R_{bacd}=-R_{abdc}=R_{cdab}$ and the first Bianchi identity $R_{abcd}+R_{acdb}+R_{adbc}=0$. Suppose only that $R$ lies in the image of the Young projector of shape $\{2,2\}$, which is one irreducible representation of $GL(d)$ inside the $\{2,2\}$ sector of $S_4$ acting on the four slots. Does that projector impose Bianchi, and so cut the count of independent components below what the pair symmetries leave?

For the tableau with rows $\{1,3\},\{2,4\}$ the projector is

$$
P_{\{2,2\}}\,T=\frac{2}{4!}\;\mathcal A_{\{1,2\}}\mathcal A_{\{3,4\}}\,\mathcal S_{\{1,3\}}\mathcal S_{\{2,4\}}\,T,
$$

where $\mathcal S$ sums over the permutations of a row of slots, $\mathcal A$ sums over the permutations of a column with their signs, and $2$ is the dimension of the $S_4$ representation $\{2,2\}$.

The TensorNetworks projector acts on a random integer tensor in $d=4$, and the code searches the slot-permuted copies of the result for three-term sums that vanish, with slot 1 fixed. It then counts independent components for $d=2,\dots,5$ two ways: from the pair symmetries alone, and as the rank of $P_{\{2,2\}}$ on all $d^4$ tensors:

```wl
Needs["Wolfram`TensorNetworks`Symmetry`"];
riemannTableau = YoungTableau[{{1, 3}, {2, 4}}];
SeedRandom[7]; projectedT = YoungProject[RandomInteger[{-9, 9}, {4, 4, 4, 4}], riemannTableau];
threeTermRelations = Select[Subsets[Select[Permutations[Range[4]], First[#] == 1 &], {3}],
   Total[Transpose[projectedT, #] & /@ #] == ConstantArray[0, {4, 4, 4, 4}] &];
slotTerm[σ_] := Subscript["R", StringJoin[Characters["abcd"][[σ]]]];
pairCount[dd_] := Length[SymmetrizedIndependentComponents[{dd, dd, dd, dd}, {{Cycles[{{1, 2}}], -1}, {Cycles[{{3, 4}}], -1}, {Cycles[{{1, 3}, {2, 4}}], 1}}]];
projectorRank[dd_] := MatrixRank[Flatten[YoungProject[ArrayReshape[#, {dd, dd, dd, dd}], riemannTableau]] & /@ IdentityMatrix[dd^4]];
componentCounts = Table[{dd, pairCount[dd], projectorRank[dd]}, {dd, 2, 5}];
ListPlot[{componentCounts[[All, {1, 2}]], componentCounts[[All, {1, 3}]]}, Joined -> True, PlotMarkers -> Automatic,
  PlotLegends -> {"pair symmetries only", "rank of \!\(\*SubscriptBox[\(P\), \({2, 2}\)]\)"}, Frame -> True,
  FrameLabel -> {"dimension d", "independent components"}, FrameTicks -> {{Automatic, None}, {Range[2, 5], None}},
  PlotLabel -> Column[{Row[{"vanishing three-term sums:  ", Column[Row[{Row[slotTerm /@ #, " + "], " = 0"}] & /@ threeTermRelations]}],
     Row[{"rank = \!\(\*SuperscriptBox[\(d\), \(2\)]\)(\!\(\*SuperscriptBox[\(d\), \(2\)]\)-1)/12 for every d:  ", componentCounts[[All, 3]] == (#^2 (#^2 - 1)/12 & /@ Range[2, 5])}]}],
  ImageSize -> 420]
```

The search returns only the cyclic sum over slots 2, 3, 4 (once with each sign), which is the first Bianchi identity; nothing beyond the projector was supplied. The rank of the projector is the textbook count $d^2(d^2-1)/12$. It equals the pair-symmetry count for $d\le3$ and falls below it from $d=4$ on: only in four or more dimensions does Bianchi remove components.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Estimating a qubit's spectrum by Schur sampling

Take $n$ copies of a qubit in an unknown state $\rho=\tfrac12(\mathbb 1+\vec b\cdot\vec\sigma)$ and measure only which Young diagram $\lambda$, that is which irreducible representation of the group $S_n$ permuting the copies, the joint state occupies. What does the outcome reveal about $\rho$?

The outcome probabilities are

$$
\Pr(\lambda)=\operatorname{Tr}\!\left[P_\lambda\,\rho^{\otimes n}\right],\qquad P_\lambda=MM^{+},\quad M=\sum_{T}P_T ,
$$

where $P_T$ is the Young projector of a standard tableau $T$ of shape $\lambda$, $M^+$ is the pseudoinverse, and $P_\lambda$ is the orthogonal projector onto the $\lambda$ block of $(\mathbb C^2)^{\otimes n}$. For qubits only diagrams with at most two rows occur.

The code builds $P_\lambda$ from the TensorNetworks Young projectors, computes $\Pr(\lambda)$ with $\vec b$ symbolic, then substitutes $b_z=\sqrt{|\vec b|^2-b_x^2-b_y^2}$:

```wl
Needs["Wolfram`TensorNetworks`Symmetry`"];
schurProjector[n_, λ_] := With[{tableaux = Select[YoungTableau[TakeList[#, λ]] & /@ Permutations[Range[n]], StandardTableauQ]},
   With[{m = Transpose[Flatten[Total[Function[t, YoungProject[ArrayReshape[#, ConstantArray[2, n]], t]] /@ tableaux]] & /@ IdentityMatrix[2^n]]},
    m . PseudoInverse[m]]];
schurRho = (IdentityMatrix[2] + bx PauliMatrix[1] + by PauliMatrix[2] + bz PauliMatrix[3])/2;
schurProbabilities[n_] := AssociationMap[Simplify[Total[schurProjector[n, #] Transpose[KroneckerProduct @@ ConstantArray[schurRho, n]], 2] /.
      bz -> Sqrt[bLength^2 - bx^2 - by^2]] &, IntegerPartitions[n, 2]];
schurPanel[n_] := With[{pr = schurProbabilities[n]},
   Plot[Evaluate[Values[pr]], {bLength, 0, 1}, PlotRange -> {0, 1}, PlotLegends -> (Row[{"λ = ", #}] & /@ Keys[pr]), Frame -> True,
    FrameLabel -> {"Bloch length |b|", "Pr(λ)"}, PlotLabel -> Row[{"n = ", n, ",   depends on |b| only: ", FreeQ[pr, bx | by]}], ImageSize -> 300]];
Grid[{schurPanel /@ {2, 4, 6}}]
```

Every probability depends on $\vec b$ only through $|\vec b|$, so the diagram reveals the eigenvalues $(1\pm|\vec b|)/2$ and nothing about the eigenbasis, and a pure state always gives the one-row diagram $\{n\}$. As $|\vec b|$ grows the most probable diagram moves to a longer first row, but at six copies the estimate is coarse and biased high: even the maximally mixed state most often gives $\{4,2\}$ rather than $\{3,3\}$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Snyder's Lorentz-covariant quantum spacetime

Snyder made the spacetime coordinates $x^\mu$ noncommuting without breaking Lorentz symmetry: their commutator is a Lorentz generator $M^{\mu\nu}$ times the square of a fundamental length $a$. Is there a quadratic combination of the ten generators that commutes with all of them?

With $\eta=\operatorname{diag}(1,-1,-1,-1)$, the algebra is

$$
[x^\mu,x^\nu]=ia^2M^{\mu\nu},\qquad [M^{\mu\nu},x^\rho]=i\left(\eta^{\nu\rho}x^\mu-\eta^{\mu\rho}x^\nu\right),\qquad [M^{\mu\nu},M^{\rho\sigma}]=i\left(\eta^{\nu\rho}M^{\mu\sigma}-\eta^{\nu\sigma}M^{\mu\rho}-\eta^{\mu\rho}M^{\nu\sigma}+\eta^{\mu\sigma}M^{\nu\rho}\right),
$$

and the candidates are $x^2=\eta_{\mu\nu}x^\mu x^\nu$, the Lorentz Casimir $\tfrac12M_{\mu\nu}M^{\mu\nu}$, and the combination $x^2-\tfrac{a^2}{2}M_{\mu\nu}M^{\mu\nu}$.

`NonCommutativeAlgebra` holds the relations, and `NonCommutativeExpand` reduces the commutator of each candidate with each generator to an ordered basis; each entry says whether the reduced commutator is zero:

```wl
snyderEta[μ_, ν_] := DiagonalMatrix[{1, -1, -1, -1}][[μ + 1, ν + 1]];
lorentzM[μ_, ν_] := Which[μ == ν, 0, μ < ν, snyderM[μ, ν], True, -snyderM[ν, μ]];
lorentzPairs = Subsets[Range[0, 3], {2}];
snyderGenerators = Join[snyderX /@ Range[0, 3], snyderM @@@ lorentzPairs];
snyderAlgebra = NonCommutativeAlgebra[<|"Generators" -> {snyderGenerators}, "ScalarVariables" -> {snyderA},
    "CommutationRelations" -> Join[
      Commutator[snyderX[#2], snyderX[#1]] + I snyderA^2 snyderM[#1, #2] & @@@ lorentzPairs,
      Flatten[Table[Commutator[snyderM[μ, ν], snyderX[ρ]] - I (snyderEta[ν, ρ] snyderX[μ] - snyderEta[μ, ρ] snyderX[ν]), {μ, 0, 3}, {ν, μ + 1, 3}, {ρ, 0, 3}]],
      Replace[Subsets[lorentzPairs, {2}], {{μ_, ν_}, {ρ_, σ_}} :> Commutator[snyderM[ρ, σ], snyderM[μ, ν]] +
         I (snyderEta[ν, ρ] lorentzM[μ, σ] - snyderEta[ν, σ] lorentzM[μ, ρ] - snyderEta[μ, ρ] lorentzM[ν, σ] + snyderEta[μ, σ] lorentzM[ν, ρ]), {1}]]|>];
xSquared = Sum[snyderEta[μ, μ] snyderX[μ] ** snyderX[μ], {μ, 0, 3}];
lorentzCasimir = Sum[snyderEta[μ, μ] snyderEta[ν, ν] snyderM[μ, ν] ** snyderM[μ, ν], {μ, 0, 3}, {ν, μ + 1, 3}];
Grid[Prepend[MapThread[Prepend, {Outer[If[NonCommutativeExpand[Commutator[#1, #2], snyderAlgebra] === 0, 0, "\[NotEqual] 0"] &,
      {xSquared, lorentzCasimir, xSquared - snyderA^2 lorentzCasimir}, snyderGenerators],
     {"\!\(\*SuperscriptBox[\(x\), \(2\)]\)", "\!\(\*FractionBox[\(1\), \(2\)]\)\!\(\*SubscriptBox[\(M\), \(\[Mu]\[Nu]\)]\)\!\(\*SuperscriptBox[\(M\), \(\[Mu]\[Nu]\)]\)",
      "\!\(\*SuperscriptBox[\(x\), \(2\)]\) - \!\(\*FractionBox[SuperscriptBox[\(a\), \(2\)], \(2\)]\)\!\(\*SubscriptBox[\(M\), \(\[Mu]\[Nu]\)]\)\!\(\*SuperscriptBox[\(M\), \(\[Mu]\[Nu]\)]\)"}}],
   Prepend[Replace[snyderGenerators, {snyderX[μ_] :> Superscript["x", μ], snyderM[μ_, ν_] :> Superscript["M", Row[{μ, ν}]]}, {1}], "commutator with"]],
  Frame -> All, Spacings -> {1.5, 1}, BaseStyle -> 16]
```

$x^2$ and the Lorentz Casimir each commute with every $M^{\mu\nu}$, as Lorentz scalars must, but neither commutes with the coordinates. Their combination $x^2-\tfrac{a^2}{2}M_{\mu\nu}M^{\mu\nu}$ commutes with all ten generators: Snyder's spacetime keeps a central quadratic element, the Minkowski $x^2$ corrected at order $a^2$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Dirac matrices for a curved metric

At one point of a curved spacetime the metric $g_{\mu\nu}$ is a real symmetric matrix, in general not diagonal, with $p$ positive and $n-p$ negative eigenvalues. Does the vielbein construction give gamma matrices obeying the Clifford relation and the flat-space contraction identities for every such metric?

The gamma matrices must satisfy the Clifford relation, and a vielbein $e$ reduces $g$ to its signature $\eta=\operatorname{diag}(\pm1)$:

$$
\{\gamma_\mu,\gamma_\nu\}=2g_{\mu\nu}\mathbb 1,\qquad g=e^{T}\eta\,e,\qquad \gamma_\mu=e^{a}{}_{\mu}\,\Gamma_a ,
$$

where the flat $\Gamma_a$ are $2^{\lfloor n/2\rfloor}$-dimensional tensor products of Pauli matrices (the Brauer-Weyl construction), multiplied by $i$ where $\eta_{aa}=-1$. Indices are raised with $g^{\mu\nu}$, and in flat space $\gamma^\mu\gamma_\nu\gamma_\mu=(2-n)\gamma_\nu$.

`GammaMatrices` from the DiracMatrix paclet builds $\gamma_\mu$ for one random non-diagonal metric per $(n,p)$, with eigenvalue magnitudes between $0.1$ and $10$, for $n=2,\dots,8$ and $p=0,\dots,n$. The code measures the error of the Clifford relation relative to the largest metric entry, and the error of $\gamma^\mu\gamma_\nu\gamma_\mu=(2-n)\gamma_\nu$ relative to the largest entry of $\gamma_\nu$:

```wl
Needs["WolframQuantumComputation`DiracMatrix`", "/Users/mohammadb/Nextcloud/Mads/Quantum/DiracMatrix/DiracMatrixPackage/Kernel/DiracMatrix.wl"];
gammaChecks[n_, p_] := With[{g = RandomCurvedMetric[p, n - p, 1, 1]},
   With[{γ = GammaMatrices[g]}, With[{γUp = Inverse[g] . γ},
     {{n, Max[Abs[Outer[#1 . #2 + #2 . #1 &, γ, γ, 1] - 2 TensorProduct[g, IdentityMatrix[2^Floor[n/2]]]]]/Max[Abs[g]]},
      {n, Max[Table[Max[Abs[Total[MapThread[Dot[#1, γ[[ν]], #2] &, {γUp, γ}]] - (2 - n) γ[[ν]]]]/Max[Abs[γ[[ν]]]], {ν, n}]]}}]]];
SeedRandom[1]; gammaResults = Flatten[Table[gammaChecks[n, p], {n, 2, 8}, {p, 0, n}], 1];
ListLogPlot[{gammaResults[[All, 1]], gammaResults[[All, 2]]}, PlotRange -> {10^-17, 10^-12}, Frame -> True, ImageSize -> 420,
  PlotMarkers -> Automatic, FrameLabel -> {"dimension n", "relative error"}, PlotLabel -> "random metrics of every signature",
  PlotLegends -> {"\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\) + \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\) = 2\!\(\*SubscriptBox[\(g\), \(\[Mu]\[Nu]\)]\)",
    "\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\) = (2 - n)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\)"}]
```

Both the Clifford relation and $\gamma^\mu\gamma_\nu\gamma_\mu=(2-n)\gamma_\nu$ hold to rounding error in every dimension and signature tried, so the vielbein construction gives valid gamma matrices for an arbitrary metric, and the flat-space contraction identity carries over unchanged.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A particle bound by its own gravity

A particle of mass $m$ moves in the Newtonian potential $\Phi$ sourced by its own mass density $m|\psi|^2$: quantum spreading pushes the wave out, self-gravity pulls it in. What is the energy of its stationary bound state once the finite computational box of radius $R$ is removed?

The Schrödinger-Newton equations are

$$
i\hbar\,\partial_t\psi=-\frac{\hbar^2}{2m}\nabla^2\psi+m\Phi\,\psi,\qquad \nabla^2\Phi=4\pi Gm\,|\psi|^2 .
$$

With $\hbar=m=G=1$, lengths in units of $\hbar^2/Gm^3$, a spherically symmetric stationary state $\psi=e^{-iEt}\,u(r)/(\sqrt{4\pi}\,r)$ obeys

$$
\Phi''+\frac{2}{r}\Phi'=\frac{u^2}{r^2},\qquad -\frac12u''+\Phi\,u=E\,u,\qquad u(0)=u(R)=\Phi(R)=0,\quad \int_0^R u^2\,dr=1 .
$$

Far from the particle $\Phi\approx-1/r$, so setting $\Phi(R)=0$ raises the energy by about $1/R$. Moroz, Penrose and Tod (1998) give the free-space ground-state energy $E=-0.163\,G^2m^5/\hbar^2$.

Finite elements on a mesh of 2000 elements graded toward $r=0$ solve Newton's equation for $\Phi$ given $u$, then the lowest eigenpair given $\Phi$, repeated from a normalized Gaussian until $E$ changes by less than $10^{-7}$. The code does this in four boxes:

```wl
Needs["NDSolve`FEM`"];
selfPotential[uR_, mesh_] := NDSolveValue[{Laplacian[Φ[r], {r, θ, φ}, "Spherical"] == uR^2/r^2,
    DirichletCondition[Φ[r] == 0, r == mesh["Bounds"][[1, 2]]]}, Φ, {r} ∈ mesh];
solitonGroundState[ϕ_, mesh_] := First @ SortBy[Transpose @ NDEigensystem[{SchrodingerPDEComponent[{u[r], {r}}, <|"ReducedPlanckConstant" -> 1, "SchrodingerPotential" -> ϕ[r]|>],
       DirichletCondition[u[r] == 0, True]}, u[r], {r} ∈ mesh, 3,
      Method -> {"Eigensystem" -> {"Arnoldi", "MaxIterations" -> 3000, "BasisSize" -> 50, "Shift" -> -10}}], Re @* First];
solitonEnergy[boxR_] := With[{mesh = ToGradedMesh[{Line[{{0}, {boxR}}], <|"Alignment" -> "Left", "ElementCount" -> 2000|>}]},
   Re @ First @ NestWhile[solitonGroundState[selfPotential[Last[#], mesh], mesh] &,
      solitonGroundState[selfPotential[2^(7/4)/(3^(3/2) Pi^(1/4)) r Exp[-r^2/9], mesh], mesh], Abs[First[#1] - First[#2]] > 10^-7 &, 2, 60]];
boxEnergies = {#, solitonEnergy[#]} & /@ {20, 40, 80, 160};
mptEnergy = -0.163;
ListPlot[{boxEnergies, {#1, #2 - 1/#1} & @@@ boxEnergies, {{0, mptEnergy}, {170, mptEnergy}}}, Joined -> {False, False, True},
  PlotStyle -> {Automatic, Automatic, Directive[Thick, Dotted]}, PlotMarkers -> {{Automatic, 12}, {Automatic, 12}, None},
  PlotLegends -> {"E(R)", "E(R) - 1/R", "Moroz-Penrose-Tod"}, PlotRange -> {{0, 170}, All}, Frame -> True, ImageSize -> 420,
  FrameLabel -> {"box radius R (\!\(\*SuperscriptBox[\(\[HBar]\), \(2\)]\)/(G \!\(\*SuperscriptBox[\(m\), \(3\)]\)))", "energy E (\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(m\), \(5\)]\)/\!\(\*SuperscriptBox[\(\[HBar]\), \(2\)]\))"}]
```

A bound state exists in every box, and its energy falls toward the free-space value as the box grows. Subtracting $1/R$ puts every box on the Moroz-Penrose-Tod energy (dotted): the whole finite-box error is the cut-off $-1/r$ tail.
