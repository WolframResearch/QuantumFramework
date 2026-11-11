<!-- Candidate tensor-network slides for Modeling-Slides.md, in the deck's slide format. Each slide's cell is self-contained once the loading cell below has run. Not yet merged into the deck. -->

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/TensorNetworks/TensorNetworks"];
Needs["Wolfram`TensorNetworks`"]; Needs["Wolfram`TensorNetworks`Symmetry`"];
rowLegend[l_] := l /. Placed[x_, ___] :> x /. (h : LineLegend | PointLegend | SwatchLegend)[x_, y_, o___] :> h[x, y, LegendLayout -> "Row", o];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### How much entanglement can one cut of a GHZ, W or Dicke state carry?

Split $n$ qubits after the first $l$. The GHZ state and the Dicke states with $k$ excitations,

$$
|\mathrm{GHZ}_n\rangle=\frac{|0\cdots0\rangle+|1\cdots1\rangle}{\sqrt2},\qquad |D_n^k\rangle=\binom{n}{k}^{-1/2}\sum_{|x|=k}|x\rangle,\qquad |W_n\rangle=|D_n^1\rangle,
$$

are all multipartite entangled. Does the entanglement entropy of a cut grow with $n$, and for which of them?

A Dicke state has Schmidt weights $p_j=\binom{l}{j}\binom{n-l}{k-j}\big/\binom{n}{k}$, the chance that $j$ of the $k$ excitations sit left of the cut, so $S=-\sum_j p_j\ln p_j$; GHZ has $S=\ln 2$ at every cut. Each state is written exactly as a matrix product state whose bond counts the excitations to its left (bond dimension $k+1$; 2 for GHZ), and the entropy is read off each bond. Entropies are in nats; dashed curves are the formula above and, on the right, the large-$n$ limit $\frac12\ln\frac{\pi e n}{8}$ of the middle cut of $|D_n^{n/2}\rangle$:

```wl
dickeA[k_] := N @ Table[Boole[m2 == m1 + s], {m1, 0, k}, {m2, 0, k}, {s, 0, 1}];
chain[first_, bulk_, last_, n_] := TensorNetwork[Join[{first}, ConstantArray[bulk, n - 2], {last}],
   Join[{{1, n}}, Table[{i - 1, i, n - 1 + i}, {i, 2, n - 1}], {{n - 1, 2 n - 1}}], Range[n, 2 n - 1]];
dicke[n_, k_] := With[{a = dickeA[k]}, chain[a[[1]], a, a[[All, -1]], n]];
ghz[n_] := chain[N @ IdentityMatrix[2], N @ {{{1, 0}, {0, 0}}, {{0, 0}, {0, 1}}}, N @ IdentityMatrix[2], n];
hyperS[n_, k_, l_] := -Total[# Log[#] &[DeleteCases[PDF[HypergeometricDistribution[l, k, n], Range[0, k]], 0]]];
n = 24; ks = {1, 2, n/2}; ns = Range[4, 64, 4];
labels = {"GHZ", "W", Row[{"Dicke k = ", 2}], Row[{"Dicke k = n/2"}]};
profile = ListPlot[Join[{Table[{l, MPSEntanglementEntropy[ghz[n], l]}, {l, n - 1}]},
     Table[{l, MPSEntanglementEntropy[dicke[n, k], l]}, {k, ks}, {l, n - 1}]], PlotMarkers -> Automatic, PlotLegends -> labels,
    Frame -> True, FrameLabel -> {"l", "S"}, PlotLabel -> Row[{"n = ", n}], PlotRange -> {0, All}, ImageSize -> 330];
profileF = ListLinePlot[Join[{Table[{l, Log[2.]}, {l, n - 1}]}, Table[{l, hyperS[n, k, l]}, {k, ks}, {l, n - 1}]], PlotStyle -> Directive[Gray, Dashed]];
middle = ListPlot[Join[{Table[{m, MPSEntanglementEntropy[ghz[m], m/2]}, {m, ns}]},
     Table[{m, MPSEntanglementEntropy[dicke[m, kf[m]], m/2]}, {kf, {1 &, 2 &, #/2 &}}, {m, ns}]], PlotMarkers -> Automatic,
    Frame -> True, FrameLabel -> {"n", "S at l = n/2"}, PlotLabel -> "middle cut", PlotRange -> {0, All}, ImageSize -> 330];
middleF = Plot[Log[Pi E m/8]/2, {m, 4, 64}, PlotStyle -> Directive[Gray, Dashed]];
Legended[Row[{Show[First[profile], profileF], Show[middle, middleF]}, Spacer[10]],
  Placed[Row[{rowLegend[Last[profile]], LineLegend[{Directive[Gray, Dashed]}, {"exact Schmidt formula and large-n limit"}]}, Spacer[20]], Below]]
```

GHZ carries exactly one ebit, $\ln 2$, at every cut and every $n$, and W carries less, $h(l/n)$ with $h$ the binary entropy, peaking at $\ln 2$ at the middle. A Dicke state with a fixed number $k$ of excitations has only $k+1$ Schmidt values, so its entropy stays below $\ln(k+1)$ for every $n$. The half-filled Dicke state is the one whose middle cut keeps growing, along $\frac12\ln\frac{\pi e n}{8}$, because the number of excitations on the left half fluctuates by $\sqrt n/4$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Is the AKLT spin-1 chain disordered?

A ring of $L$ spin-1 sites with the Affleck-Kennedy-Lieb-Tasaki Hamiltonian

$$
H=\sum_{i}\left[\mathbf S_i\cdot\mathbf S_{i+1}+\frac13\left(\mathbf S_i\cdot\mathbf S_{i+1}\right)^2\right]
$$

has a unique ground state with no magnetic order: spin correlations die exponentially. Is there an order that does not die?

The ground state is the matrix product state with $A^{+}=\sqrt{2/3}\,\sigma_+$, $A^{0}=-\sqrt{1/3}\,\sigma_z$, $A^{-}=-\sqrt{2/3}\,\sigma_-$, and on the infinite chain

$$
C(r)=\langle S^z_0S^z_r\rangle=\frac43\left(-\frac13\right)^{r},\qquad O(r)=\left\langle S^z_0\,e^{i\pi\sum_{0<k<r}S^z_k}\,S^z_r\right\rangle=-\frac49 .
$$

Both correlators are computed on a ring of $L=40$ sites by contracting the ket, the bra and the operator insertions into a number; dashed lines are the infinite-chain values:

```wl
aklt = Transpose[N @ {Sqrt[2/3] {{0, 1}, {0, 0}}, -Sqrt[1/3] {{1, 0}, {0, -1}}, -Sqrt[2/3] {{0, 0}, {1, 0}}}, {3, 1, 2}];
sz = {1, 0, -1}; one = {1, 1, 1}; L = 40; rs = Range[16];
ring[w_] := TensorNetwork[Join[(aklt . DiagonalMatrix[#] &) /@ w, ConstantArray[aklt, L]],
   Join[Table[{Mod[i - 2, L] + 1, i, 2 L + i}, {i, L}], Table[{L + Mod[i - 2, L] + 1, L + i, 2 L + i}, {i, L}]], {}];
norm = TensorNetworkContract[ring[ConstantArray[one, L]], "Greedy"];
c[r_] := TensorNetworkContract[ring[Join[{sz}, ConstantArray[one, r - 1], {sz}, ConstantArray[one, L - r - 1]]], "Greedy"]/norm;
o[r_] := TensorNetworkContract[ring[Join[{sz}, ConstantArray[(-1)^sz, r - 1], {sz}, ConstantArray[one, L - r - 1]]], "Greedy"]/norm;
data = ListLogPlot[{{#, Abs[c[#]]} & /@ rs, {#, Abs[o[#]]} & /@ rs}, PlotMarkers -> Automatic,
   PlotLegends -> {Row[{"|C(r)|, two-point"}], Row[{"|O(r)|, string"}]}, Frame -> True, FrameLabel -> {"r", "|correlator|"}, ImageSize -> 420];
Legended[Show[First[data], LogPlot[{4/3 3^-r, 4/9}, {r, 1, 16}, PlotStyle -> Directive[Gray, Dashed]], PlotRange -> All],
  Placed[Row[{rowLegend[Last[data]], LineLegend[{Directive[Gray, Dashed]}, {Row[{"(4/3) ", Superscript[3, "-r"], " and 4/9"}]}]}, Spacer[20]], Below]]
```

The two-point correlator falls by a factor of 3 per site, a correlation length of $1/\ln 3$ sites, so no spin knows the orientation of a distant one. The string correlator stays at $4/9$ at every distance: after a site with $S^z=+1$, the next nonzero spin along the chain is always $-1$, however many zeros sit between them. The disorder is only in where the zeros fall; the alternation of the nonzero spins is perfect.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Why an exact 2D partition function costs exponentially in the width, not the area

The Ising model on an $L\times L$ square lattice with open boundaries has

$$
Z=\sum_{\{s_i=\pm1\}}\exp\Big(K\sum_{\langle ij\rangle}s_is_j\Big)=\sum_{\text{bonds}}\prod_{\text{sites}}T^{(i)},\qquad T^{(i)}_{a_1\cdots a_m}=\sum_{s=\pm1}\prod_{\mu=1}^{m}\big(\sqrt W\big)_{s\,a_\mu},\quad W=\begin{pmatrix}e^{K}&e^{-K}\\e^{-K}&e^{K}\end{pmatrix},
$$

one tensor per site with one leg per bond. Computing $Z$ exactly is a sequence of pairwise contractions. How large must the biggest intermediate tensor be?

Every leg has dimension 2, so a tensor with $m$ open legs has $2^m$ entries, for any $K$. The figure compares three orders: sweeping the sites row by row, a greedy search for a cheap order, and sweeping the sites in a random order; dashed is $L+1$:

```wl
isingTN[L_, K_] := With[{sites = Tuples[Range[L], 2], q = MatrixPower[N @ {{Exp[K], Exp[-K]}, {Exp[-K], Exp[K]}}, 1/2]},
   With[{bonds = Join[{#, # + {0, 1}} & /@ Select[sites, Last[#] < L &], {#, # + {1, 0}} & /@ Select[sites, First[#] < L &]]},
    With[{legs = Position[bonds, #, {2}][[All, 1]] & /@ sites},
     TensorNetwork[Total[Function[row, Outer[Times, ##] & @@ ConstantArray[row, Length[#]]] /@ q] & /@ legs, legs, {}]]]];
peak[legs_, path_] := Max[Length @* Last /@ Rest[FoldList[Function[{st, p},
       Append[Delete[st, List /@ p], Complement[Union @@ st[[p]], Intersection @@ st[[p]]]]], legs, path]]];
sweep[n_] := Prepend[Table[{1, n - k + 1}, {k, 2, n - 1}], {1, 2}];
Ls = Range[2, 12];
cost = Table[With[{tn = isingTN[L, 0.4]}, With[{legs = TensorNetworkIndices[tn] /. Superscript[_, b_] :> b},
     {peak[legs, sweep[L^2]], peak[legs, GreedyContractionPath[tn]], BlockRandom[SeedRandom[L]; peak[legs[[RandomSample[Range[L^2]]]], sweep[L^2]]]}]], {L, Ls}];
data = ListLogPlot[Transpose[MapThread[Thread[{#1, #2}] &, {Ls, cost}]], PlotMarkers -> Automatic, Joined -> True,
   PlotLegends -> {"row-by-row sweep", "greedy search", "random sweep"}, Frame -> True,
   FrameLabel -> {"L", Row[{Subscript["log", 2], " entries of the largest intermediate"}]}, ImageSize -> 420];
Legended[Show[First[data], LogPlot[L + 1, {L, 2, 12}, PlotStyle -> Directive[Gray, Dashed]], PlotRange -> All],
  Placed[Row[{rowLegend[Last[data]], LineLegend[{Directive[Gray, Dashed]}, {"L + 1"}]}, Spacer[20]], Below]]
```

Sweeping row by row never holds more than $L+1$ open legs: the tensor carries one row of $L$ vertical bonds, plus the one horizontal bond of the site being absorbed. That tensor is the row-to-row transfer vector, a state of $L$ spins, so the exact $Z$ costs $2^{L+1}$ numbers however tall the lattice is. The greedy search stays within a few legs of this, while a random order keeps many cut bonds open at once and its largest tensor grows with the area.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### What the entanglement of a critical chain measures: the central charge

The transverse-field Ising chain of $L$ spins with open ends,

$$
H=-\sum_{i=1}^{L-1}\sigma^z_i\sigma^z_{i+1}-g\sum_{i=1}^{L}\sigma^x_i ,
$$

is critical at $g=1$. Conformal field theory predicts that the entanglement entropy of the first $l$ spins is

$$
S(l)=\frac c6\ln\!\left[\frac{2L}{\pi}\sin\frac{\pi l}{L}\right]+s_0 ,
$$

with central charge $c=\frac12$ for this chain. Does the ground state follow it, and what happens away from $g=1$?

For $L=18$ the ground state is found by exact diagonalization and written as a matrix product state by successive Schmidt decompositions; the entropy is read off each bond. Plotted against $x=\frac16\ln\big[\frac{2L}{\pi}\sin\frac{\pi l}{L}\big]$, the prediction is a line of slope $c$; the dashed line has slope $\frac12$, with $s_0$ chosen to match the $g=1$ data on average:

```wl
toMPS[psi_, n_] := With[{steps = FoldList[Function[{st, k}, With[{m = ArrayReshape[Last[st], {2 Length[Last[st]], Length[Flatten[Last[st]]]/(2 Length[Last[st]])}]},
        With[{svd = SingularValueDecomposition[m, Min[Dimensions[m]]]}, With[{r = Count[Diagonal[svd[[2]]], x_ /; x > 10^-12]},
          {ArrayReshape[svd[[1]][[All, ;; r]], {Length[m]/2, 2, r}], svd[[2]][[;; r, ;; r]] . ConjugateTranspose[svd[[3]][[All, ;; r]]]}]]]],
      {None, {psi}}, Range[n - 1]]},
   TensorNetwork[Join[{Transpose[steps[[2, 1, 1]]]}, Transpose[#[[1]], {1, 3, 2}] & /@ steps[[3 ;;]], {ArrayReshape[steps[[-1, 2]], {Length[steps[[-1, 2]]], 2}]}],
    Join[{{1, n}}, Table[{i - 1, i, n - 1 + i}, {i, 2, n - 1}], {{n - 1, 2 n - 1}}], Range[n, 2 n - 1]]];
op[m_, i_, L_] := KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2, SparseArray], L], i -> SparseArray[m]];
hTFIM[g_, L_] := -Sum[op[PauliMatrix[3], i, L] . op[PauliMatrix[3], i + 1, L], {i, L - 1}] - g Sum[op[PauliMatrix[1], i, L], {i, L}];
L = 18; gs = {1, 1.5, 3}; x = Table[Log[2 L/Pi Sin[Pi l/L]]/6, {l, L - 1}];
ent = Table[With[{mps = toMPS[First[Eigenvectors[N[hTFIM[g, L]] - (1 + g) L IdentityMatrix[2^L, SparseArray], 1]], L]},
     Table[MPSEntanglementEntropy[mps, l], {l, L - 1}]], {g, gs}];
s0 = Mean[First[ent] - x/2];
data = ListPlot[Thread[{x, #}] & /@ ent, PlotMarkers -> Automatic, PlotLegends -> (Row[{"g = ", #}] & /@ gs),
   Frame -> True, FrameLabel -> {Row[{"x = ", TraditionalForm[HoldForm[Log[2 L/Pi Sin[Pi l/L]]/6]]}], "S(l)"}, PlotLabel -> Row[{"L = ", L}], ImageSize -> 420];
Legended[Show[First[data], Plot[s/2 + s0, {s, Min[x], Max[x]}, PlotStyle -> Directive[Gray, Dashed]], PlotRange -> All],
  Placed[Row[{rowLegend[Last[data]], LineLegend[{Directive[Gray, Dashed]}, {"slope c = 1/2"}]}, Spacer[20]], Below]]
```

At $g=1$ the entropies fall on the line of slope $\frac12$ except at the cuts next to the ends, where the open boundary adds its own correction: the entanglement grows logarithmically with the size of the cut, and the rate is the central charge of the Ising field theory. Away from the critical point the chain has a gap, the entropy stops growing once the cut is longer than the correlation length, and the points bend over to a constant; the farther from $g=1$, the lower the plateau.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Can n copies of a qubit reveal its spectrum without choosing a basis?

Take $n$ copies of a qubit in the state $\rho$ with eigenvalues $p\ge 1-p$, and measure only which irreducible representation of the permutation group $S_n$ (equivalently, which total spin) the copies occupy, labelled by a two-row Young diagram $\lambda=(n-k,k)$. Schur-Weyl duality gives

$$
\Pr(\lambda)=f^{\lambda}\,s_\lambda(p,1-p),\qquad s_{(n-k,k)}(x,y)=(xy)^k\sum_{j=0}^{n-2k}x^{j}y^{\,n-2k-j},
$$

where $f^\lambda$ is the number of standard Young tableaux of shape $\lambda$ (the hook-length formula). Does $\lambda_1/n$ estimate $p$?

The measurement commutes with every rotation $U^{\otimes n}$, so no basis is chosen. The figure shows $\Pr(\lambda)$ against $\lambda_1/n$ for three values of $p$; dashed is $\lambda_1/n=p$:

```wl
prob[n_, p_] := Table[{(n - k)/n, TableauDimension[DeleteCases[{n - k, k}, 0]] (p (1 - p))^k Sum[p^j (1 - p)^(n - 2 k - j), {j, 0, n - 2 k}]}, {k, 0, Floor[n/2]}];
{ps, ns} = {{0.5, 0.7, 0.9}, {10, 40, 160}};
panel[p_] := ListLinePlot[prob[#, p] & /@ ns, PlotMarkers -> Automatic, PlotLegends -> (Row[{"n = ", #}] & /@ ns),
   Epilog -> {Directive[Gray, Dashed], Line[{{p, 0}, {p, 1}}]}, PlotRange -> {{0.45, 1.02}, {0, All}},
   Frame -> True, FrameLabel -> {Row[{Subscript["\[Lambda]", 1], "/n"}], "Pr(\[Lambda])"}, PlotLabel -> Row[{"p = ", p}], ImageSize -> 260];
panels = panel /@ ps;
Legended[Row[First /@ panels, Spacer[5]],
  Placed[Row[{rowLegend[Last[First[panels]]], LineLegend[{Directive[Gray, Dashed]}, {Row[{Subscript["\[Lambda]", 1], "/n = p"}]}]}, Spacer[20]], Below]]
```

As $n$ grows the distribution of $\lambda_1/n$ narrows onto the larger eigenvalue $p$, so the shape of the Young diagram alone estimates the spectrum. At finite $n$ it sits above $p$: the estimate is biased toward a purer state, and at $p=\frac12$ every outcome has $\lambda_1/n\ge\frac12$, so the bias shrinks only as the distribution narrows.
