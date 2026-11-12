# Audit: "Adding two momenta on a noncommuting spacetime"

## Verdict

The κ-Minkowski half of the slide is right. The Snyder half reaches the wrong answer to the slide's own question: in Snyder's spacetime the sum of two momenta does depend on their order, and the slide finds that it does not only because it throws away the Lorentz part of $\log(e^Xe^Y)$. That part cannot be dropped, because $M^{\mu\nu}$ does not commute with $x^\mu$, so $e^{i\gamma\cdot x+i\omega\cdot M}$ is not a plane wave of momentum $\gamma$. Once the Lorentz factor is split off properly, $e^Xe^Y=e^{i(\alpha\oplus\beta)\cdot x}e^{i\omega\cdot M}$, the sum reproduces the published Snyder addition law term by term and is neither commutative nor associative. This is not a difference of convention: Meljanac, Meljanac, Mignemi and Štrajn show the order dependence at first order in $a^2$ for every realization of a Snyder spacetime, and a relabeling of momenta cannot turn a noncommutative sum into a commutative one. The code computes exactly what the prose describes, and the figure shows what the closing paragraph says; the error sits in the definition of the Snyder momentum, so the figure and the closing paragraph faithfully show a wrong result. A corrected paragraph, closing paragraph and cell are at the end; the corrected cell runs with no messages and its figure is [Momentum-Composition-Audit-corrected.png](Momentum-Composition-Audit-corrected.png) (the current figure is [Momentum-Composition-Audit-current.png](Momentum-Composition-Audit-current.png)).

## Sources read in full

- M. V. Battisti and S. Meljanac, "Scalar field theory on non-commutative Snyder space-time", Phys. Rev. D 82, 024028 (2010), arXiv:1003.2108.
- S. Meljanac, D. Meljanac, S. Mignemi and R. Štrajn, "Snyder-type spacetimes, twisted Poincaré algebra and addition of momenta", Int. J. Mod. Phys. A 32, 1750172 (2017), arXiv:1608.06207.
- J. Kowalski-Glikman and S. Nowak, "Doubly special relativity theories as different bases of κ-Poincaré algebra", Phys. Lett. B 539, 126 (2002), arXiv:hep-th/0203040.

The PDFs and their text are in the side-session notes folder, `momcomp/papers/`. All kernel checks are the scripts in `momcomp/` there (`checks.wls`, `kappa.wls`, `fact.wls`, `snyder5d.wls`), run with `wolframscript -file` against the local QF paclet loaded as in the deck's title slide.

## Findings

### 1. The Snyder sum is order dependent; the slide's symmetry comes from discarding the Lorentz part (real error)

What the slide does: it defines $(\alpha\oplus\beta)_\mu$ as the coefficient of $ix^\mu$ in $Z=\log(e^Xe^Y)$ and says the Lorentz generators in $Z$ "are not part of the momentum". Its argument that this coefficient is symmetric is correct: $x\to-x$, $M\to M$ preserves Snyder's relations, so the $x$ part of $Z$ sits in odd orders, and $e^Ye^X=e^{-Z(-X,-Y)}$ makes odd orders symmetric. The kernel confirms this to all orders: with the exact matrix logarithm in a faithful $5\times5$ representation (below), the coefficient is symmetric to machine precision for six random pairs of 4-momenta of size up to $0.6/a$, and the $\{x,M\}$ decomposition of the logarithm has residual below $10^{-28}$, the noise floor of the working precision.

Why it is the wrong object: the product $e^Xe^Y=e^{i\gamma\cdot x+i\omega\cdot M}$ is not the plane wave $e^{i\gamma\cdot x}$ times something that can be ignored. The literature reads the momentum of a product off its action on the Lorentz-invariant vacuum: Battisti and Meljanac write $e^{ik\cdot\tilde x}e^{iq\cdot\tilde x}\triangleright 1=e^{iD(k,q)\cdot x}$ (their eqs. 21 to 23), and in the coset picture momentum space is de Sitter space $SO(4,1)/SO(3,1)$, which is not a group (their Sec. III.C). The Lorentz generators annihilate the vacuum, but $e^{i\gamma\cdot x+i\omega\cdot M}\triangleright1\neq e^{i\gamma\cdot x}\triangleright1$ because $M$ and $x$ do not commute. The momentum of the product is $\alpha\oplus\beta$ defined by splitting the Lorentz factor to the right, $e^Xe^Y=e^{i(\alpha\oplus\beta)\cdot x}e^{i\omega\cdot M}$.

Evidence:

- Faithful matrices. With $G=\operatorname{diag}(1,-1,-1,-1,-1)$ and $L^{AB}=i\,(G_{AA}E_{AB}-G_{BB}E_{BA})$, the choice $x^\mu=aL^{\mu4}$, $M^{\mu\nu}=L^{\mu\nu}$ satisfies all three of the slide's commutation relations exactly (residual 0; the other three sign choices fail). This is the de Sitter algebra, Snyder's original momentum space. The $M^{\mu\nu}$ annihilate $o=(0,0,0,0,1)$, so the momentum of $g$ is read from the point $g\,o$: $e^{i\gamma\cdot x}o=\cosh\theta\,o+\frac{\sinh\theta}{\theta}(-\gamma_0,\gamma_1,\gamma_2,\gamma_3,0)$ with $\theta^2=\gamma\cdot\gamma$, which inverts in closed form.
- Exact sums at the slide's own points ($a=1$, $\alpha=s(1,\frac12,0,0)$, $\beta=s(\frac35,-\frac45,0,0)$): at $s=0.9$, $\alpha\oplus\beta=(1.337,-0.050)$ and $\beta\oplus\alpha=(1.446,-0.553)$, while the slide's symmetric definition gives the single point $(1.360,-0.295)$. For six random pairs the two orders differ by $0.001$ to $0.25$ in units of $1/a$.
- Third order, symbolic. Expanding the coset sum, it equals Battisti and Meljanac's first-order coproduct in the Weyl realization (their eq. 27 with $c=1/6$, the realization in which $e^{ik\cdot\tilde x}\triangleright1=e^{ik\cdot x}$), term by term, once their metric $\operatorname{diag}(-,+,+,+)$ is translated to the slide's and their $s$ is identified with $a^2$. Its antisymmetric part is
  $$
  \alpha\oplus\beta-\beta\oplus\alpha=-\frac{a^2}{2}\left(\alpha^2\beta-\beta^2\alpha+(\alpha\cdot\beta)(\beta-\alpha)\right),
  $$
  the same structure as Meljanac et al. eq. 31, whose antisymmetric part is $\frac{\beta s}{2}\left[(k\cdot q)(q-k)+k^2q-q^2k\right]$ for every realization (their $s_1$, $s_2$), vanishing only for $s=0$, commutative spacetime (their eq. 33). It vanishes when $\beta\parallel\alpha$, as it must.
- The slide's third-order sum is exactly the symmetric part of the true one: $\frac12(\alpha\oplus\beta+\beta\oplus\alpha)$ (kernel: difference $\{0,0,0,0\}$). Dropping the Lorentz part removes precisely the order dependence the slide asks about.
- Through fourth order the correct sum is the coefficient of $ix^\mu$ in $X+Y-\frac16[X,[X,Y]]+\frac13[Y,[Y,X]]$ (kernel: equal to the coset sum, difference $\{0,0,0,0\}$). It follows from BCH: the $x$ part of $\log(e^{i\gamma\cdot x}e^{i\omega\cdot M})$ at third order contains $\frac12[i\gamma\cdot x,i\omega\cdot M]$ with $i\omega\cdot M=\frac12[X,Y]$, so $\gamma=Z_x-\frac14([X,[X,Y]]-[Y,[Y,X]])_x$.

Reduction modulo conventions: any change of momentum labels $k\to K(k)$ maps a commutative sum to a commutative sum, and Meljanac et al. eq. 31 covers the general first-order realization, so no relabeling turns the published Snyder law into the slide's symmetric one. The slide's sign of $a^2$ and its metric match Snyder's de Sitter case (checked through the $5\times5$ representation), so the disagreement is not a sign convention either.

### 2. Nonassociativity is not mentioned (omission)

The slide's question is only about order, so this is not an error, but it is the sharpest contrast between the two spacetimes. The true Snyder sum is not associative: $(\alpha\oplus\beta)\oplus\gamma-\alpha\oplus(\beta\oplus\gamma)$ is $0.007$ to $0.17$ for six random triples of size up to $0.6/a$, as Battisti and Meljanac (Sec. III.C) and Meljanac et al. (eqs. 40 to 42) state. The κ-Minkowski sum is associative to $10^{-15}$, since κ-Minkowski momentum space is a group. The replacement cell prints this associator for one triple with an out-of-plane third momentum, and the proposed closing paragraph states it.

### 3. κ-Minkowski: correct, with the basis left implicit (minor)

Every κ-Minkowski claim checks out:

- $x^0=\frac{i}{\kappa}E_{11}$, $x^1=E_{12}$ satisfy $[x^0,x^1]=\frac{i}{\kappa}x^1$.
- The slide's closed form for $(\alpha\oplus\beta)_1$ agrees with a numerical matrix logarithm to $3\times10^{-15}$ for eight random pairs, including negative energies, so the `PowerExpand` is safe for real momenta.
- The order splitting is exactly proportional to $\alpha_0\beta_j-\alpha_j\beta_0$ to all orders: the exact difference divided by $\alpha_1\beta_0-\alpha_0\beta_1$ is $\frac{(\alpha_0+\beta_0)(e^{\alpha_0/\kappa}-1)(e^{\beta_0/\kappa}-1)}{\alpha_0\beta_0(e^{(\alpha_0+\beta_0)/\kappa}-1)}$, free of the spatial components, with series $\frac1\kappa-\frac{\alpha_0\beta_0}{12\kappa^3}+\dots$ matching the fourth-order label.
- The fourth-order points lie within $2.4\times10^{-3}\kappa$ of the exact ones at $s=0.9$, invisible at the plot's scale, so "the exact sums sit on the fourth-order points" is fair.

Conventions: taking the coefficient of $ix^\mu$ in $Z$ is the symmetric (Weyl) ordering of plane waves. Writing the plane wave with $x^0$ to the right, $e^{ik_jx^j}e^{ik_0x^0}$, gives in the kernel $k\oplus q=(k_0+q_0,\ k_j+e^{-k_0/\kappa}q_j)$, which is the bicrossproduct coproduct $\Delta p_i=p_i\otimes1+e^{-p_0/\kappa}\otimes p_i$, $\Delta p_0=p_0\otimes1+1\otimes p_0$ of Kowalski-Glikman and Nowak (their eq. 5). Their commutator is $[x_0,x_i]=-\frac{i}{\kappa}x_i$ with metric $\operatorname{diag}(-1,1,1,1)$, the opposite sign to the slide's, but the two lead to the same composition law, so the sign difference is one of metric and plane-wave convention. The two momentum labels are related by $k_0=\gamma_0$, $k_j=\gamma_j\,\kappa(1-e^{-\gamma_0/\kappa})/\gamma_0$. Order dependence is basis independent here because the momentum group is non-abelian. The slide never says which basis it uses; one sentence fixes that (included below).

### 4. Units and parameters (no error)

$a$ (length) and $\kappa$ (mass) are used consistently, the plot units $1/a$ and $\kappa$ match the substitutions `composeA -> 1` and `composeK -> 1`, and $\frac12[X,Y]=-\frac{i}{2}a^2\alpha_\mu\beta_\nu M^{\mu\nu}$ is correct. Notation agrees with the neighboring Snyder Casimir slide (same metric, same $[x,x]=ia^2M$, same $[M,x]$).

### 5. The question (well posed once the momentum of a product is defined)

"Does the deformed sum depend on which momentum comes first?" is well posed and has a sharp answer once the momentum of a product of plane waves is defined by its action on the vacuum, which is the literature's definition. The answer is yes in both spacetimes. With the slide's definition the Snyder answer is no, which is the error in Finding 1.

### 6. Low-momentum limit (correct)

Both sums reduce to $\alpha+\beta$ as $s\to0$: the Snyder correction is cubic in momentum times $a^2$, the κ correction quadratic over $\kappa$. In the corrected figure the exact Snyder sums begin to separate from the fourth-order points at $s=0.6$ and clearly leave them at $s=0.9$, where the momenta approach $1/a$ and the fifth-order term matters.

## Proposed corrections

### Replacement for the paragraph after the BCH formula

Current paragraph:

> and $(\alpha\oplus\beta)_\mu$ is the coefficient of $ix^\mu$ in $Z$. In Snyder's algebra $Z$ also carries Lorentz generators, starting with $\frac12[X,Y]=-\frac{i}{2}a^2\alpha_\mu\beta_\nu M^{\mu\nu}$; they are not part of the momentum. Snyder's relations are unchanged under $x\to-x$, $M\to M$, so the momentum collects only the odd orders of $Z$, and since $e^Ye^X=e^{-Z(-X,-Y)}$ these are symmetric under exchange to all orders. In $\kappa$-Minkowski the matrices $x^0=\frac{i}{\kappa}E_{11}$, $x^1=E_{12}$ obey the relations, so $e^Xe^Y$ can be multiplied and its logarithm taken in closed form.

Proposed (three paragraphs and one display):

> In $\kappa$-Minkowski spacetime $Z$ is a combination of the $x^\mu$ alone, and $(\alpha\oplus\beta)_\mu$ is the coefficient of $ix^\mu$ in $Z$, the momentum of the symmetrically ordered plane wave. The matrices $x^0=\frac{i}{\kappa}E_{11}$, $x^1=E_{12}$ obey the relations, so $e^Xe^Y$ can be multiplied and its logarithm taken in closed form. Writing the same plane wave with $x^0$ to the right, $e^{ik_jx^j}e^{ik_0x^0}$, relabels the momentum as $k_0=\gamma_0$, $k_j=\gamma_j\,\kappa\left(1-e^{-\gamma_0/\kappa}\right)/\gamma_0$, and turns the sum into the bicrossproduct law $k\oplus q=\left(k_0+q_0,\ k_j+e^{-k_0/\kappa}q_j\right)$. Back in the symmetric labels the exact sum is $(\alpha\oplus\beta)_j=\left(\alpha_je^{\beta_0/\kappa}\varphi(\alpha_0/\kappa)+\beta_j\varphi(\beta_0/\kappa)\right)/\varphi\left((\alpha_0+\beta_0)/\kappa\right)$ with $\varphi(x)=(e^x-1)/x$: each spatial momentum is weighted by a factor that depends only on the energies.
>
> In Snyder's algebra $Z$ also carries Lorentz generators, starting with $\frac12[X,Y]=-\frac{i}{2}a^2\alpha_\mu\beta_\nu M^{\mu\nu}$, and they cannot be dropped: $M^{\mu\nu}$ does not commute with $x^\rho$, so $e^{i\gamma_\mu x^\mu+i\omega_{\mu\nu}M^{\mu\nu}}$ is not a plane wave of momentum $\gamma$. Instead the product is split as $e^Xe^Y=e^{i(\alpha\oplus\beta)_\mu x^\mu}e^{i\omega_{\mu\nu}M^{\mu\nu}}$; the Lorentz factor leaves the Lorentz-invariant vacuum unchanged, so $\alpha\oplus\beta$ is the momentum of the product. Comparing the two Baker-Campbell-Hausdorff series, $(\alpha\oplus\beta)_\mu$ is the coefficient of $ix^\mu$ in
>
> $$
> X+Y-\frac16[X,[X,Y]]+\frac13[Y,[Y,X]],
> $$
>
> and since Snyder's relations are unchanged under $x\to-x$, $M\to M$, the next correction is fifth order. The $5\times5$ matrices of the de Sitter algebra, $x^\mu=aL^{\mu4}$ and $M^{\mu\nu}=L^{\mu\nu}$ with metric $\operatorname{diag}(1,-1,-1,-1,-1)$, obey Snyder's relations and give the exact sum: $e^Xe^Y$ carries the point $(0,0,0,0,1)$ of de Sitter space to a new point, and $\alpha\oplus\beta$ is the $\gamma$ for which $e^{i\gamma_\mu x^\mu}$ carries it to the same place. Each spatial $x^j$ rotates that point in a compact plane, so these coordinates label points uniquely only while $a\sqrt{-\gamma\cdot\gamma}<\pi$; beyond that, different $\gamma$ carry $(0,0,0,0,1)$ to the same point. The momenta below stay well inside that range.

### Replacement for the sentence before the cell

Current:

> Each set of relations reduces $Z$ for general $\alpha$, $\beta$; the figure shows two momenta in the $(p_0,p_1)$ plane, in units of $1/a$ for Snyder and of $\kappa$ for $\kappa$-Minkowski, for $\alpha=s\,(1,\frac12,0,0)$ and $\beta=s\,(\frac35,-\frac45,0,0)$ at three scales $s$, with the exact $\kappa$-Minkowski sum from the $2\times2$ matrices:

Proposed:

> Each set of relations reduces the series for general $\alpha$, $\beta$; the figure shows two momenta in the $(p_0,p_1)$ plane, in units of $1/a$ for Snyder and of $\kappa$ for $\kappa$-Minkowski, for $\alpha=s\,(1,\frac12,0,0)$ and $\beta=s\,(\frac35,-\frac45,0,0)$ at three scales $s$, with the exact sums from the $5\times5$ and $2\times2$ matrices:

### Replacement cell

Runs with `wolframscript -file` against the local QF paclet in about 16 s with no messages (most of it in the symbolic checks). The κ-Minkowski computation gives the same result as the original; its series and exact sum are restructured so that both spacetimes share one set of definitions. What is new: Snyder's relations come from one family, the $[L,L]$ relations of the de Sitter algebra with $x^\mu=aL^{\mu4}$, so the series algebra and the $5\times5$ matrices are the same object; the Snyder momentum is read after the Lorentz factor is split off; both exact sums are closed forms in $\alpha$, $\beta$ (and $a$, $\kappa$), written with hypergeometric forms of $\sinh\theta/\theta$ and $(e^x-1)/x$ so that they need no branch choice and stay finite at zero energy; four checks are printed (the $5\times5$ matrices obey the relations, the $2\times2$ product reproduces the exact κ sum, and each exact sum agrees with its series to the order kept), plus a fifth showing a spatial Snyder momentum wrapping around; the Snyder splitting is displayed as computed coefficients of $\alpha$ and $\beta$; and one line shows that grouping three momenta matters for Snyder and not for κ-Minkowski. The cached series is named `composeSeriesSum` and keyed by the spacetime's name, so the cell gives the same result whether or not the original cell ran first in the same kernel.

```wl
composeEta = DiagonalMatrix[{1, -1, -1, -1}]; composeG = DiagonalMatrix[{1, -1, -1, -1, -1}];
composeL[A_, A_] = 0; composeL[A_, B_] := If[A < B, composeJ[A, B], -composeJ[B, A]];
composePairs = Subsets[Range[0, 4], {2}];
{composeAl, composeBe} = {composeα /@ Range[0, 3], composeβ /@ Range[0, 3]};
composeRelations["Snyder"] = Replace[Subsets[composePairs, {2}], {{A_, B_}, {C_, D_}} :> Commutator[composeJ[C, D], composeJ[A, B]] +
     I (composeG[[B + 1, C + 1]] composeL[A, D] - composeG[[B + 1, D + 1]] composeL[A, C] - composeG[[A + 1, C + 1]] composeL[B, D] + composeG[[A + 1, D + 1]] composeL[B, C]), {1}];
composeRelations["κ-Minkowski"] = Join[Commutator[composeT[#], composeT[0]] + I composeT[#]/composeK & /@ Range[3], Commutator[composeT[#2], composeT[#1]] & @@@ Subsets[Range[3], {2}]];
composeAlgebra["Snyder"] = NonCommutativeAlgebra[<|"Generators" -> {composeJ @@@ composePairs}, "ScalarVariables" -> Join[{composeA}, composeAl, composeBe], "CommutationRelations" -> composeRelations["Snyder"]|>];
composeAlgebra["κ-Minkowski"] = NonCommutativeAlgebra[<|"Generators" -> {composeT /@ Range[0, 3]}, "ScalarVariables" -> Join[{composeK}, composeAl, composeBe], "CommutationRelations" -> composeRelations["κ-Minkowski"]|>];
composeFrame["Snyder"] = {composeA, composeJ[#, 4] & /@ Range[0, 3]}; composeFrame["κ-Minkowski"] = {1, composeT /@ Range[0, 3]};
composeBCH["Snyder"][x_, y_] := x + y - Commutator[x, Commutator[x, y]]/6 + Commutator[y, Commutator[y, x]]/3;
composeBCH["κ-Minkowski"][x_, y_] := x + y + Commutator[x, y]/2 + (Commutator[x, Commutator[x, y]] + Commutator[y, Commutator[y, x]])/12 - Commutator[y, Commutator[x, Commutator[x, y]]]/24;
composeSeriesSum[name_] := composeSeriesSum[name] = With[{c = First[composeFrame[name]], e = Last[composeFrame[name]]},
    Simplify[Coefficient[NonCommutativeExpand[composeBCH[name][I c composeAl . e, I c composeBe . e], composeAlgebra[name]], #]/(I c) & /@ e]];
composeSwap = {composeα -> composeβ, composeβ -> composeα};
composeDot[u_, v_] := u . composeEta . v;
composeLmat[A_, B_] := I (composeG[[A + 1, A + 1]] SparseArray[{A + 1, B + 1} -> 1, {5, 5}] - composeG[[B + 1, B + 1]] SparseArray[{B + 1, A + 1} -> 1, {5, 5}]);
composeSinhOverTheta[z_] := Hypergeometric0F1[3/2, z/4]; composeCoshTerm[z_] := Hypergeometric0F1[3/2, z/16]^2/2;
composeBoost[α_] := With[{v = I α . (composeLmat[#, 4] & /@ Range[0, 3])}, IdentityMatrix[5] + composeSinhOverTheta[composeDot[α, α]] v + composeCoshTerm[composeDot[α, α]] v . v];
composeThetaOverSinh[c_] := Hypergeometric2F1[1, 1, 3/2, (1 - c)/2];
composeExactSum["Snyder"][α_, β_] := With[{p = composeBoost[composeA α] . composeBoost[composeA β] . UnitVector[5, 5]}, composeThetaOverSinh[p[[5]]] {-p[[1]], p[[2]], p[[3]], p[[4]]}/composeA];
composePhi[x_] := Exp[x/2] Hypergeometric0F1[3/2, x^2/16];
composeExactSum["κ-Minkowski"][α_, β_] := Prepend[(Rest[α] Exp[β[[1]]/composeK] composePhi[α[[1]]/composeK] + Rest[β] composePhi[β[[1]]/composeK])/composePhi[(α[[1]] + β[[1]])/composeK], α[[1]] + β[[1]]];
composeKappaRep = {{{I/composeK, 0}, {0, 0}}, {{0, 1}, {0, 0}}};
composeChecks = {
   "5×5 matrices obey Snyder's relations" -> AllTrue[Flatten[Normal[composeRelations["Snyder"] /. {p_ ** q_ :> (p /. composeJ -> composeLmat) . (q /. composeJ -> composeLmat), composeJ -> composeLmat}]], # === 0 &],
   "2×2 matrices: e^X e^Y = exp(i(α ⊕ β)·x), exact κ sum" -> (Simplify[MatrixExp[I composeAl[[;; 2]] . composeKappaRep] . MatrixExp[I composeBe[[;; 2]] . composeKappaRep] -
        MatrixExp[I composeExactSum["κ-Minkowski"][composeAl, composeBe][[;; 2]] . composeKappaRep]] === {{0, 0}, {0, 0}}),
   "exact − series = O(5th order), Snyder" -> (Simplify[Normal[Series[composeExactSum["Snyder"][t composeAl, t composeBe] - (composeSeriesSum["Snyder"] /. {composeα[μ_] :> t composeα[μ], composeβ[μ_] :> t composeβ[μ]}), {t, 0, 4}]]] === {0, 0, 0, 0}),
   "exact − series = O(1/κ⁴), κ-Minkowski" -> (Simplify[Normal[Series[(composeExactSum["κ-Minkowski"][composeAl, composeBe] - composeSeriesSum["κ-Minkowski"]) /. composeK -> 1/ε, {ε, 0, 3}]]] === {0, 0, 0, 0}),
   "Snyder, (0, 5/2, 0, 0) ⊕ (0, 3/2, 0, 0) = (0, 4 − 2π, 0, 0) in units 1/a" -> (Simplify[composeExactSum["Snyder"][{0, 5/2, 0, 0}, {0, 3/2, 0, 0}] /. composeA -> 1] === {0, 4 - 2 Pi, 0, 0})};
composeSnyderSplit = With[{u = Array[composeU, 3], v = Array[composeV, 3], s = {composeDot[composeAl, composeAl], composeDot[composeBe, composeBe], composeDot[composeAl, composeBe]}},
   {u, v} /. First[SolveAlways[Thread[composeSeriesSum["Snyder"] - (composeSeriesSum["Snyder"] /. composeSwap) == (u . s) composeAl + (v . s) composeBe], Join[composeAl, composeBe]]]];
composeAssoc[name_, unit_, a_, b_, c_] := With[{f = Chop[N[composeExactSum[name][#1, #2] /. unit]] &}, Chop[f[f[a, b], c] - f[a, f[b, c]]]];
composeTriple = 3/5 {{1, 1/2, 0, 0}, {3/5, -4/5, 0, 0}, {1/2, 0, 1/2, 0}};
composeShow = {composeα[μ_] :> Subscript["α", μ], composeβ[μ_] :> Subscript["β", μ], composeK -> "κ", composeA -> "a"};
composePanel[name_, unit_, s_] := With[{pts = Thread[Join[composeAl, composeBe] -> s Join[{1, 1/2, 0, 0}, {3/5, -4/5, 0, 0}]]},
   With[{α = composeAl /. pts, β = composeBe /. pts},
    ListPlot[List /@ Chop[N[{α + β, composeSeriesSum[name] /. pts, composeSeriesSum[name] /. composeSwap /. pts, composeExactSum[name][α, β], composeExactSum[name][β, α]} /. unit]][[All, ;; 2]],
     PlotMarkers -> {{"\[FilledSquare]", 14}, {"\[EmptyCircle]", 16}, {"\[EmptyCircle]", 30}, {"+", 14}, {"\[Times]", 14}},
     PlotLegends -> {"α + β", "α ⊕ β, fourth order", "β ⊕ α, fourth order", "α ⊕ β, exact", "β ⊕ α, exact"},
     Epilog -> {Gray, Arrow[{{0, 0}, α[[;; 2]]}], Text["α", 0.55 α[[;; 2]], {0, -1.5}], Dashed, Arrow[{{0, 0}, β[[;; 2]]}], Text["β", 0.55 β[[;; 2]], {-1.5, 0}]},
     PlotRange -> s {{-0.1, 1.9}, {-1.1, 0.7}}, AspectRatio -> Automatic, Frame -> True, FrameLabel -> {Subscript["p", 0], Subscript["p", 1]},
     PlotLabel -> Row[{name, ",   s = ", s}], ImageSize -> 260]]];
composePanels = Table[composePanel[Sequence @@ m, s], {m, {{"Snyder", composeA -> 1}, {"κ-Minkowski", composeK -> 1}}}, {s, {0.3, 0.6, 0.9}}];
Labeled[Legended[Grid[Map[First, composePanels, {2}]], Last[composePanels[[2, 1]]]],
  Style[Column[{
     Row[{"Snyder, fourth order:  α ⊕ β − β ⊕ α = (", TraditionalForm[Factor[composeSnyderSplit[[1]] . {Superscript["α", 2], Superscript["β", 2], "α·β"}] /. composeShow], ") α + (",
       TraditionalForm[Factor[composeSnyderSplit[[2]] . {Superscript["α", 2], Superscript["β", 2], "α·β"}] /. composeShow], ") β"}],
     Row[{"κ-Minkowski, fourth order:  α ⊕ β − β ⊕ α = ", TraditionalForm[Factor[composeSeriesSum["κ-Minkowski"] - (composeSeriesSum["κ-Minkowski"] /. composeSwap)] /. composeShow]}],
     Row[{"grouping, exact, α, β, γ = ", Row[composeTriple, ",  "], " in units 1/a or κ:  (α ⊕ β) ⊕ γ − α ⊕ (β ⊕ γ) = ", NumberForm[composeAssoc["Snyder", composeA -> 1, Sequence @@ composeTriple], 2], " (Snyder),  ", composeAssoc["κ-Minkowski", composeK -> 1, Sequence @@ composeTriple], " (κ-Minkowski)"}],
     Row[{"checks:  ", Row[KeyValueMap[Row[{#1, " ", #2}] &, Association[composeChecks]], ";   "]}]}], 13, FontFamily -> "Helvetica"], Top]
```

### Replacement for the closing paragraph

Current:

> The deformed sum moves away from $\alpha+\beta$ in both spacetimes, more so at larger $s$. In Snyder's spacetime both orders land on the same point, as the odd-order argument requires. In $\kappa$-Minkowski spacetime the two orders split by a spatial term proportional to $\alpha_0\beta_j-\alpha_j\beta_0$, and the exact sums sit on the fourth-order points, so composing two momenta depends on their order.

Proposed:

> The deformed sum moves away from $\alpha+\beta$ in both spacetimes, more so at larger $s$, and in both it depends on which momentum comes first. In Snyder's spacetime the two orders split by $-\frac{a^2}{2}\left(\alpha^2\beta-\beta^2\alpha+(\alpha\cdot\beta)(\beta-\alpha)\right)$, cubic in the momenta and zero when they are parallel: the order is recorded in the Lorentz factor that the product leaves behind. In $\kappa$-Minkowski spacetime the split is already quadratic, a spatial term proportional to $\alpha_0\beta_j-\alpha_j\beta_0$. The exact $\kappa$-Minkowski sums sit on the fourth-order points; the exact Snyder sums begin to separate from theirs at $s=0.6$ and clearly leave them at the largest $s$, where the momenta approach $1/a$ and the fifth-order term is no longer small. Grouping separates the two: adding three momenta in Snyder's spacetime depends on which pair is added first, because its momentum space is de Sitter space, a coset rather than a group, while in $\kappa$-Minkowski spacetime the momenta form a group and the grouping does not matter.

## Open items from the review rounds

Three rounds of correctness review and three rounds of code-quality review were run. Every finding from the correctness rounds was fixed, but the fixes after the third correctness round were not re-checked by a fresh correctness reviewer. The final quality round found the physics and figure sound and left five items open; none changes the verdict:

- No Snyder sum exists for some pairs. The exponential of a single momentum only reaches de Sitter points with fifth coordinate $p_5\ge-1$. For example, $\alpha=(2,0,0,0)$, $\beta=(0,3,0,0)$ at $a=1$ lands at $p_5\approx-3.7$, and the cell then returns complex numbers without a message. At exactly $a\sqrt{-\gamma\cdot\gamma}=\pi$ it returns `Indeterminate`. A guard on $p_5\ge-1$ and one sentence in the prose would cover this. The plotted momenta are far from this edge.
- The grouping test is numeric at one triple. Both cases close symbolically: κ-Minkowski is exactly associative, and the leading Snyder associator at order $a^2$ is an $(\alpha\wedge\beta)\cdot\gamma$ structure.
- The cell does not assert Lorentz covariance of the Snyder sum, $\Lambda(\alpha\oplus\beta)=(\Lambda\alpha)\oplus(\Lambda\beta)$. The reviewer checked it in the kernel and it holds; κ-Minkowski breaks it.
- The exact Snyder sum is checked against the series only through fourth order, plus one collinear exact point. An independent matrix-exponential check at a non-collinear pair agrees, for example at the $s=0.9$ pair, but the cell does not show that check.
- Two possible simplifications: assert the split formula directly instead of fitting it with `SolveAlways`, and build one legend instead of six.
