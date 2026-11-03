# Perelomov section (P-001 to P-049, P2-001 to P2-058), exact computation only

Every check below is exact: QF's symbolic second-quantization algebra, exact symbolic Wolfram computation, exact spin, 2x2 group, coefficient or Jordan-Wigner fermion matrices, or a hand derivation. No truncated Fock matrix is used anywhere. Each `wl` cell runs as written after the loading lines

```wl
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

and its Output is copied from that run (`drafts/exact/P-work/driver.wl` runs every cell in `cells/` in a cleared context; outputs in `out_cells.txt`). The exact checks without a QF route are in `drafts/exact/P-work/wlchecks.wl` (outputs in `out_wl.txt`).

## Errors found

No new error in the atlas transcription. Every source correction the atlas records holds under exact computation:
- P-029: the printed middle exponents of the SU(1,1) disentangling are reversed; the corrected form holds exactly in the faithful 2x2 matrices and in QF's single-mode realization, and the printed form fails in both.
- P-019: the factor 2 in the $P\to W$ and $W\to Q$ smoothing kernels is required (exact Gaussian integrals).
- P2-013: the corrected phase equation solves the Schrodinger equation exactly; the printed sign leaves a residual proportional to $k\,\mathrm{Im}(\bar\beta\zeta)$.
- P2-017: the printed normalization is wrong already for real parameters; the corrected one gives norm 1 exactly for complex $\zeta,\alpha$ (exact Bargmann integral).
- P2-023: the trace-log expansion needs a plus sign (exact series).
- P2-054: the extra $e^{|\gamma|^2}$ in the printed generating function is wrong (QF exact matrix elements).
- P2-038: the corrected sign in the Gaussian kernel holds (normalized Gaussian, exact integral).

Convention gaps the atlas leaves implicit (formulas are right only under a specific index reading):
- P-045: with $\hat X_{ij}=\hat c_i\hat c_j$, the displayed $A$ equals $\langle \hat c_j\hat c_i\rangle=-\langle\hat c_i\hat c_j\rangle$, and $B_{ij}=\langle\hat c_j^\dagger\hat c_i\rangle-\delta_{ij}/2$. Read naively as $A_{ij}=\langle\hat X_{ij}\rangle$, $A$ has the wrong sign.
- P-039: $A_{ij}=\langle\hat a_i\hat a_j\rangle$ holds, but $B_{ij}=\langle\hat a_j^\dagger\hat a_i\rangle+\delta_{ij}/2$, the transpose of $\langle\hat X^i{}_j\rangle$ under the atlas's definition of $\hat X^i{}_j$.
- P2-012: $\hat K_2=\tfrac14(\hat p\hat q+\hat q\hat p)$ equals $-(\hat K_+-\hat K_-)/(2i)$ for $\hat K_\pm=\tfrac12\hat a^{\dagger2},\tfrac12\hat a^2$. This is a consistent su(1,1) basis only together with the source's $\hat K_1$, which flips sign too (the relabeling $\hat a\to i\hat a$).

QF defect seen while writing the cells: `BosonicNormalOrder` given a List of expressions returns pure functions for some elements instead of mapping over the list; for example `BosonicNormalOrder[{a ** ad, ad ** a}, {a, ad}]` returns `{{0 & }, SuperDagger[\[FormalA]1]**\[FormalA]1}`, while `BosonicNormalOrder[a ** ad, {a, ad}]` correctly returns `1 + SuperDagger[\[FormalA]1]**\[FormalA]1`. The cells below map it over lists explicitly.

---

### Canonical commutators, ladder operators and the multimode displacement
Topic: Displacement
Appears in: P-001

The canonical pair fixes the ladder operators and their unit commutator, and the multimode displacement is the exponential of the linear form in the ladder operators.

$$ [\hat q,\hat p]=i\hbar \hat I,\quad \hat a=\frac{\hat q+i\hat p}{\sqrt{2\hbar}},\quad \hat a^\dagger=\frac{\hat q-i\hat p}{\sqrt{2\hbar}},\quad [\hat a,\hat a^\dagger]=\hat I,\quad \hat D(\boldsymbol\alpha)=\exp\sum_i(\alpha_i \hat a_i^\dagger-\bar\alpha_i \hat a_i) $$

Audit: Correct. Checked by: QF exact (unit commutator, and $[\hat q,\hat p]=i\hbar$ with $\hat q,\hat p$ written through $\hat a,\hat a^\dagger$); hand derivation for the multimode form: the modes commute, so the exponential splits into one-mode displacements.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
q = Sqrt[hb/2] (a + ad); p = -I Sqrt[hb/2] (a - ad);
{BosonicNormalOrder[a ** ad - ad ** a, {a, ad}], BosonicNormalOrder[q ** p - p ** q, {a, ad}], DisplacementOperator[al, Infinity]}
```
Output: `{1, I*hb, E^(-(\[FormalA]1*Conjugate[al]) + al*SuperDagger[\[FormalA]1])}`. Matches: the first two entries are $[\hat a,\hat a^\dagger]$ and $[\hat q,\hat p]$, the third is $\hat D(\alpha)$ for one mode.

### Central-commutator BCH
Topic: BCH and Hadamard
Appears in: P-002 (also used by P-010, P2-035)

When the commutator of two operators commutes with both, the product of their exponentials equals the exponential of the sum times the exponential of half the commutator.

$$ e^{\hat A}e^{\hat B}=e^{[\hat A,\hat B]/2}e^{\hat A+\hat B},\qquad e^{\hat A+\hat B}=e^{-[\hat A,\hat B]/2}e^{\hat A}e^{\hat B} $$

Conditions: $[\hat A,[\hat A,\hat B]]=[\hat B,[\hat A,\hat B]]=0$; common domain for unbounded generators.
Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
BosonicExpSimplify[Exp[x ad] ** Exp[y a]]
```
Output: `E^(\[FormalA]1*y - (x*y)/2 + x*SuperDagger[\[FormalA]1])`. Matches: $[x\hat a^\dagger,y\hat a]=-xy$, so the result is $e^{[A,B]/2}e^{A+B}$; the second form is the same identity solved for $e^{A+B}$.

### Weyl composition law and adjoint
Topic: Displacement
Appears in: P-003

Two displacements compose to the displacement by the summed amplitude times a phase set by the symplectic area, and the adjoint of a displacement is the displacement by the opposite amplitude.

$$ \hat D(\alpha)\hat D(\beta)=e^{(\alpha\bar\beta-\bar\alpha\beta)/2}\hat D(\alpha+\beta),\qquad \hat D(\alpha)^\dagger=\hat D(-\alpha) $$

Audit: Correct. Checked by: QF exact; the adjoint follows because the generator $\alpha\hat a^\dagger-\bar\alpha\hat a$ is anti-Hermitian, so $\hat D(\alpha)^\dagger=\hat D(-\alpha)$, which the second output confirms is the inverse.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
d[x_] := DisplacementOperator[x, Infinity] /. Conjugate[s_Plus] :> Map[Conjugate, s];
{Simplify[BosonicExpSimplify[d[al] ** d[be]]/(Exp[(al Conjugate[be] - Conjugate[al] be)/2] d[al + be])],
 BosonicExpSimplify[d[al] ** d[-al]]}
```
Output: `{1, 1}`. Matches: the ratio of $\hat D(\alpha)\hat D(\beta)$ to $e^{(\alpha\bar\beta-\bar\alpha\beta)/2}\hat D(\alpha+\beta)$ is 1, and $\hat D(\alpha)\hat D(-\alpha)=\hat I$.

### Displacement commutation phase, group commutator and polygon phase
Topic: Displacement
Appears in: P-003

Swapping two displacements costs twice the composition phase, the group commutator is that phase times the identity, and an ordered product of several displacements is the total displacement times the phase $\mathrm{Im}\sum_{j>k}\alpha_j\bar\alpha_k$.

$$ \hat D(\alpha)\hat D(\beta)=e^{\alpha\bar\beta-\bar\alpha\beta}\hat D(\beta)\hat D(\alpha),\quad \hat D(\alpha)\hat D(\beta)\hat D(-\alpha)\hat D(-\beta)=e^{\alpha\bar\beta-\bar\alpha\beta}\hat I,\quad \hat D(\alpha_n)\cdots\hat D(\alpha_1)=e^{i\,\mathrm{Im}\sum_{j>k}\alpha_j\bar\alpha_k}\hat D\big(\textstyle\sum_j\alpha_j\big) $$

Audit: Correct (the polygon phase means the factor $e^{i\,\mathrm{Im}\sum_{j>k}\alpha_j\bar\alpha_k}$). Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
d[x_] := DisplacementOperator[x, Infinity];
ex[e_] := e /. Power[E, x_] :> Power[E, Expand[x]];
s = be Conjugate[al] + g Conjugate[al] + g Conjugate[be];
{ex[BosonicExpSimplify[d[al] ** d[be]]/BosonicExpSimplify[d[be] ** d[al]]],
 ex[BosonicExpSimplify[d[al] ** d[be] ** d[-al] ** d[-be]]],
 Simplify[ex[BosonicExpSimplify[d[g] ** d[be] ** d[al] ** d[-(al + be + g)]]]/Exp[(s - Conjugate[s])/2]]}
```
Output: `{E^(-(be*Conjugate[al]) + al*Conjugate[be]), E^(-(be*Conjugate[al]) + al*Conjugate[be]), 1}`. Matches: swap ratio and group commutator are both $e^{\alpha\bar\beta-\bar\alpha\beta}$; the three-factor product times $\hat D(-(\alpha_1+\alpha_2+\alpha_3))$ divided by $e^{i\,\mathrm{Im}\,s}$ is 1, with $s=\alpha_2\bar\alpha_1+\alpha_3\bar\alpha_1+\alpha_3\bar\alpha_2$.

### Lattice commutation condition and critical coherent lattice
Topic: Displacement
Appears in: P-014, P2-050

Two displacements commute exactly when their symplectic area is a multiple of $\pi$; a coherent-state lattice is complete when its cell area is at most $\pi$ and incomplete above.

$$ [\hat D(\omega_1),\hat D(\omega_2)]=0 \iff \mathrm{Im}(\omega_1\bar\omega_2)\in\pi\mathbb Z;\qquad \text{lattice of cell area } S:\ \text{complete for } S\le\pi,\ \text{incomplete for } S>\pi $$

Audit: Correct for the commutation condition; the completeness statements are Unchecked: cited theorem, source proof via the Weierstrass sigma function. Checked by: QF exact plus hand derivation: the group commutator is $e^{\omega_1\bar\omega_2-\bar\omega_1\omega_2}=e^{2i\,\mathrm{Im}(\omega_1\bar\omega_2)}$, which is 1 exactly when $\mathrm{Im}(\omega_1\bar\omega_2)\in\pi\mathbb Z$. P-014 and P2-050 agree (overcomplete below $\pi$, complete at $\pi$ with one deletable state).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
d[x_] := DisplacementOperator[x, Infinity];
BosonicExpSimplify[d[w1] ** d[w2] ** d[-w1] ** d[-w2]]
```
Output: `E^(-(w2*Conjugate[w1]) + w1*Conjugate[w2])`. Matches.

### Fock basis and ladder action
Topic: Fock matrix elements and generating functions
Appears in: P-005

The vacuum is killed by the annihilator, number states are normalized powers of the creator on the vacuum, and the ladder operators move one step with the square-root factors.

$$ \hat a|0\rangle=0,\quad |n\rangle=\frac{(\hat a^\dagger)^n}{\sqrt{n!}}|0\rangle,\quad \hat a|n\rangle=\sqrt n|n-1\rangle,\quad \hat a^\dagger|n\rangle=\sqrt{n+1}|n+1\rangle $$

Audit: Correct. Checked by: QF exact, with symbolic Fock indices.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicMatrixElement[{m, n}, a], BosonicMatrixElement[{m, n}, ad],
 Table[BosonicMatrixElement[{k, 0}, ad ** ad ** ad]/Sqrt[3!], {k, 0, 4}]}
```
Output: `{Sqrt[n]*KroneckerDelta[m, -1 + n], Sqrt[m]*KroneckerDelta[-1 + m, n], {0, 0, 0, 1, 0}}`. Matches: $\langle m|\hat a|n\rangle=\sqrt n\,\delta_{m,n-1}$, $\langle m|\hat a^\dagger|n\rangle=\sqrt m\,\delta_{m-1,n}$, and $(\hat a^\dagger)^3|0\rangle/\sqrt{3!}=|3\rangle$.

### Displacement in position representation
Topic: Displacement
Appears in: P-006

The displacement factors into a momentum kick and a position translation with a half-area phase, so on wavefunctions it shifts the argument by $Q$ and multiplies by a plane-wave phase.

$$ \hat D(\alpha)=e^{-iPQ/(2\hbar)}e^{iP\hat q/\hbar}e^{-iQ\hat p/\hbar},\quad [\hat D(\alpha)\psi](q)=e^{iP(q-Q/2)/\hbar}\psi(q-Q),\quad \alpha=\frac{Q+iP}{\sqrt{2\hbar}} $$

Audit: Correct. Checked by: QF exact for the operator factorization ($\hbar=1$, $\hat q=(\hat a+\hat a^\dagger)/\sqrt2$, $\hat p=-i(\hat a-\hat a^\dagger)/\sqrt2$); hand derivation for the action: $e^{-iQ\hat p}$ translates by $Q$ and $e^{iP\hat q}$ multiplies by $e^{iPq}$, which with the phase $e^{-iPQ/2}$ gives $e^{iP(q-Q/2)}\psi(q-Q)$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2]; p = -I (a - ad)/Sqrt[2];
Simplify[Exp[-I P Q/2] BosonicExpSimplify[Exp[I P q] ** Exp[-I Q p]]/DisplacementOperator[(Q + I P)/Sqrt[2], Infinity],
  Element[{P, Q}, Reals]]
```
Output: `1`. Matches: $e^{-iPQ/2}e^{iP\hat q}e^{-iQ\hat p}\,\hat D(\alpha)^{-1}=1$ for real $P,Q$.

### Overlap of displaced fiducial states
Topic: Coherent states (general)
Appears in: P-007

The overlap of two displaced copies of one fiducial state is the Weyl phase times the expectation of the relative displacement in the fiducial.

$$ \langle\beta;\psi_0|\alpha;\psi_0\rangle=e^{(\alpha\bar\beta-\bar\alpha\beta)/2}\langle\psi_0|\hat D(\alpha-\beta)|\psi_0\rangle $$

Conditions: normalized fiducial; resolution of identity $\int d^2\alpha\,|\alpha;\psi_0\rangle\langle\alpha;\psi_0|/\pi=\hat I$ in the weak sense.
Audit: Correct for the overlap; the resolution of identity is Unchecked here (standard Schur-orthogonality result). Checked by: QF exact for the operator identity $\hat D(\beta)^\dagger\hat D(\alpha)=\hat D(-\beta)\hat D(\alpha)=e^{(\alpha\bar\beta-\bar\alpha\beta)/2}\hat D(\alpha-\beta)$, which gives the overlap for every fiducial.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
d[x_] := DisplacementOperator[x, Infinity];
Simplify[BosonicExpSimplify[d[-be] ** d[al]]/(Exp[(al Conjugate[be] - Conjugate[al] be)/2] d[al - be])]
```
Output: `1`. Matches.

### Displacement conjugation of a, q, p
Topic: Displacement
Appears in: P-008

Conjugating by a displacement shifts the annihilator by $\alpha$ and the quadratures by $Q$ and $P$.

$$ \hat D(\alpha)^\dagger \hat a\hat D(\alpha)=\hat a+\alpha\hat I,\quad \hat D^\dagger \hat q\hat D=\hat q+Q\hat I,\quad \hat D^\dagger \hat p\hat D=\hat p+P\hat I $$

Audit: Correct. Checked by: QF exact, with $\hat q,\hat p$ written through $\hat a,\hat a^\dagger$ and symbolic $\hbar$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
dag[x_] := BosonicExpSimplify[Exp[-(g ad - Conjugate[g] a)] ** x ** Exp[g ad - Conjugate[g] a]];
q = Sqrt[hb/2] (a + ad); p = -I Sqrt[hb/2] (a - ad);
Simplify[{dag[a] - a, dag[q] - q, dag[p] - p} /. g -> (Q + I P)/Sqrt[2 hb], Element[{P, Q}, Reals] && hb > 0]
```
Output: `{(I*P + Q)/(Sqrt[2]*Sqrt[hb]), Q, P}`. Matches: the shifts are $\alpha=(Q+iP)/\sqrt{2\hbar}$, $Q$ and $P$.

### Bogoliubov combination keeps the canonical commutator
Topic: Bosonic canonical transformations
Appears in: P-009

A linear combination of $\hat a$ and $\hat a^\dagger$ is again a canonical annihilator exactly when $|u|^2-|v|^2=1$.

$$ \hat b=u\hat a+v\hat a^\dagger,\qquad [\hat b,\hat b^\dagger]=(|u|^2-|v|^2)\hat I=\hat I \text{ when } |u|^2-|v|^2=1 $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
b = u a + v ad; bd = Conjugate[u] ad + Conjugate[v] a;
BosonicNormalOrder[b ** bd - bd ** b, {a, ad}]
```
Output: `u*Conjugate[u] - v*Conjugate[v]`. Matches.

### Displacement: normal and antinormal factorizations
Topic: Displacement
Appears in: P-010, P2-053

The displacement in normal order carries the scalar $e^{-|\alpha|^2/2}$ and in antinormal order carries $e^{+|\alpha|^2/2}$.

$$ \hat D(\alpha)=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}e^{-\bar\alpha \hat a}=e^{+|\alpha|^2/2}e^{-\bar\alpha \hat a}e^{\alpha \hat a^\dagger} $$

Audit: Correct. Checked by: QF exact (the built-in ordered forms and `BosonicExpOrder` of the generator agree).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{DisplacementOperator[al, Infinity, "Ordering" -> "Normal"], DisplacementOperator[al, Infinity, "Ordering" -> "Antinormal"],
 BosonicExpOrder[Exp[al ad - Conjugate[al] a]], BosonicExpOrder[Exp[al ad - Conjugate[al] a], "Ordering" -> "Antinormal"]}
```
Output: `{E^(al*SuperDagger[\[FormalA]1])**E^(-(\[FormalA]1*Conjugate[al]))/E^((al*Conjugate[al])/2), E^((al*Conjugate[al])/2)*E^(-(\[FormalA]1*Conjugate[al]))**E^(al*SuperDagger[\[FormalA]1]), <same two forms again>}`. Matches.

### Differential correspondences for multiplying a displacement
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: P-011

Multiplying a displacement by a ladder operator on either side equals a first-order Wirtinger derivative of the displacement plus a linear term.

$$ \hat a^\dagger \hat D=(\partial_\alpha+\bar\alpha/2)\hat D,\quad \hat D\hat a^\dagger=(\partial_\alpha-\bar\alpha/2)\hat D,\quad \hat a\hat D=(-\partial_{\bar\alpha}+\alpha/2)\hat D,\quad \hat D\hat a=(-\partial_{\bar\alpha}-\alpha/2)\hat D $$

Conditions: $\alpha,\bar\alpha$ independent Wirtinger variables.
Audit: Correct. Checked by: QF exact for the two commutation steps, and hand derivation from the normal form $\hat D=e^{-\alpha\bar\alpha/2}e^{\alpha\hat a^\dagger}e^{-\bar\alpha\hat a}$ with $\alpha,\bar\alpha$ independent: $\partial_\alpha\hat D=(\hat a^\dagger-\bar\alpha/2)\hat D$ and $\partial_{\bar\alpha}\hat D=-\hat D(\hat a+\alpha/2)$ give the first and fourth lines; $e^{-\bar\alpha\hat a}\hat a^\dagger=(\hat a^\dagger-\bar\alpha)e^{-\bar\alpha\hat a}$ and $\hat a e^{\alpha\hat a^\dagger}=e^{\alpha\hat a^\dagger}(\hat a+\alpha)$ turn them into the second and third.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicExpSimplify[Exp[-alb a] ** ad ** Exp[alb a]], BosonicExpSimplify[Exp[-al ad] ** a ** Exp[al ad]]}
```
Output: `{-alb + SuperDagger[\[FormalA]1], \[FormalA]1 + al}`. Matches: $e^{-\bar\alpha\hat a}\hat a^\dagger e^{\bar\alpha\hat a}=\hat a^\dagger-\bar\alpha$ and $e^{-\alpha\hat a^\dagger}\hat a e^{\alpha\hat a^\dagger}=\hat a+\alpha$.

### Coherent state expansion, eigenvalue equation, overlap and resolution
Topic: Coherent states (general)
Appears in: P-012

A coherent state is the normalized exponential of $\alpha\hat a^\dagger$ on the vacuum, it is an eigenvector of $\hat a$, two coherent states overlap by a Gaussian, and coherent projectors resolve the identity.

$$ |\alpha\rangle=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}|0\rangle=e^{-|\alpha|^2/2}\sum_{n}\frac{\alpha^n}{\sqrt{n!}}|n\rangle,\quad \hat a|\alpha\rangle=\alpha|\alpha\rangle,\quad \langle\alpha|\beta\rangle=e^{-(|\alpha|^2+|\beta|^2)/2+\bar\alpha\beta},\quad |\langle\alpha|\beta\rangle|^2=e^{-|\alpha-\beta|^2},\quad \int\frac{d^{2N}\alpha}{\pi^N}|\boldsymbol\alpha\rangle\langle\boldsymbol\alpha|=\hat I $$

Audit: Correct. Checked by: QF exact for the expansion (symbolic $n$), eigenvalue and overlap; exact WL computation for the resolution: the angular integral kills off-diagonal terms and the radial integral $\int_0^\infty 2r\,e^{-r^2}r^{2n}dr/n!$ returns 1 for every $n$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicMatrixElement[{n, 0}, Exp[al ad]], BosonicExpSimplify[a ** Exp[al ad]],
 BosonicMatrixElement[{0, 0}, Exp[Conjugate[al] a] ** Exp[be ad]]}
```
Output: `{Sqrt[n!]*Piecewise[{{al^n/n!, n >= 0}}, 0], E^(al*SuperDagger[\[FormalA]1])**(\[FormalA]1 + al), E^(be*Conjugate[al])}`. Matches: $\langle n|e^{\alpha\hat a^\dagger}|0\rangle=\alpha^n/\sqrt{n!}$; $\hat a e^{\alpha\hat a^\dagger}|0\rangle=\alpha e^{\alpha\hat a^\dagger}|0\rangle$; $\langle0|e^{\bar\alpha\hat a}e^{\beta\hat a^\dagger}|0\rangle=e^{\bar\alpha\beta}$, which with the two normalizations is the stated overlap.

### Bargmann realization and displacement as weighted translation
Topic: Coherent states (general)
Appears in: P-013

In the Bargmann space $\hat a^\dagger$ is multiplication by $z$ and $\hat a$ is $\partial_z$, and a displacement acts as a weighted translation of the entire function.

$$ \hat a\leftrightarrow\partial_z,\quad \hat a^\dagger\leftrightarrow z,\quad \|f\|^2=\int\frac{d^2z}{\pi}e^{-|z|^2}|f(z)|^2,\quad [\hat D(\alpha)f](z)=e^{-|\alpha|^2/2+\alpha z}f(z-\bar\alpha) $$

Audit: Correct, including the boundary example $e^{z^2/2}$. Checked by: QF exact plus hand derivation: $\hat D f(\hat a^\dagger)|0\rangle=e^{-|\alpha|^2/2}e^{\alpha\hat a^\dagger}e^{-\bar\alpha\hat a}f(\hat a^\dagger)|0\rangle=e^{-|\alpha|^2/2}e^{\alpha\hat a^\dagger}f(\hat a^\dagger-\bar\alpha)|0\rangle$, which reads $e^{-|\alpha|^2/2+\alpha z}f(z-\bar\alpha)$. Exact WL computation: the norm series of $e^{z^2/2}$, $\sum_k(2k)!/(4^k(k!)^2)$, diverges (`SumConvergence` returns False).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{DisplacementOperator[al, Infinity, "Ordering" -> "Normal"], BosonicExpSimplify[Exp[-Conjugate[al] a] ** ad ** Exp[Conjugate[al] a]]}
```
Output: `{E^(al*SuperDagger[\[FormalA]1])**E^(-(\[FormalA]1*Conjugate[al]))/E^((al*Conjugate[al])/2), -Conjugate[al] + SuperDagger[\[FormalA]1]}`. Matches: the normal form of $\hat D$, and $e^{-\bar\alpha\hat a}\hat a^\dagger e^{\bar\alpha\hat a}=\hat a^\dagger-\bar\alpha$.

### Coordinate to Bargmann kernel
Topic: Fock matrix elements and generating functions
Appears in: P-013

The kernel that takes a position wavefunction to its Bargmann function is a Gaussian in $z$ and $q$.

$$ K(z,q)=(\pi\hbar)^{-1/4}\exp\left(-\frac{z^2}{2}+\sqrt{\frac2\hbar}zq-\frac{q^2}{2\hbar}\right) $$

Audit: Correct. Checked by: exact WL computation: $\sum_n z^n\psi_n(q)/\sqrt{n!}$ with Hermite functions, summed in closed form through the Hermite generating function $\sum_n t^nH_n(x)/n!=e^{2tx-t^2}$ at $t=z/\sqrt2$, $x=q/\sqrt\hbar$, equals the stated kernel for symbolic $\hbar$.
QF code (exact): none. QF's second-quantization layer has no position representation of states.

### Covariant (Q) and contravariant (P) symbols from ordered expansions
Topic: Normal symbols and normal-order exponentials
Appears in: P-015

The covariant symbol is the coherent expectation value and is read from the normal-ordered expansion; the contravariant symbol is the weight in the coherent-projector expansion and is read from the antinormal-ordered expansion.

$$ Q_{\hat A}(\alpha)=\langle\alpha|\hat A|\alpha\rangle,\ \ \hat A=\sum A_{mn}(\hat a^\dagger)^m\hat a^n\Rightarrow Q_{\hat A}=\sum A_{mn}\bar\alpha^m\alpha^n;\quad \hat A=\int\frac{d^2\alpha}{\pi}P_{\hat A}|\alpha\rangle\langle\alpha|,\ \ \hat A=\sum \widetilde A_{mn}\hat a^n(\hat a^\dagger)^m\Rightarrow P_{\hat A}=\sum \widetilde A_{mn}\bar\alpha^m\alpha^n $$

Audit: Correct. Checked by: QF exact for the covariant rule; hand derivation for the contravariant rule ($\int\alpha^n\bar\alpha^m|\alpha\rangle\langle\alpha|d^2\alpha/\pi=\hat a^n(\int|\alpha\rangle\langle\alpha|d^2\alpha/\pi)\hat a^{\dagger m}$).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
BosonicMatrixElement[{al, al}, ad ** ad ** a, "Basis" -> "Coherent"] /. Abs[al]^2 -> al Conjugate[al]
```
Output: `al*Conjugate[al]^2`. Matches: $Q_{\hat a^{\dagger2}\hat a}=\bar\alpha^2\alpha$.

### Q symbol is the Gaussian smoothing of the P symbol
Topic: Normal symbols and normal-order exponentials
Appears in: P-015

The covariant symbol is the contravariant symbol smoothed with the unit Gaussian.

$$ Q_{\hat A}(\alpha)=\int\frac{d^2\beta}{\pi}e^{-|\alpha-\beta|^2}P_{\hat A}(\beta) $$

Audit: Correct. Checked by: QF exact plus exact WL integral, on $\hat a\hat a^\dagger$, whose contravariant symbol is $|\beta|^2$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
qsym = BosonicMatrixElement[{al, al}, a ** ad, "Basis" -> "Coherent"];
smooth = Integrate[Exp[-((x - x0)^2 + (y - y0)^2)] (x^2 + y^2), {x, -Infinity, Infinity}, {y, -Infinity, Infinity}]/Pi;
Simplify[ComplexExpand[qsym /. al -> x0 + I y0] - smooth]
```
Output: `0`. Matches: the covariant symbol from QF minus the Gaussian smoothing of $|\beta|^2$ is 0.

### Coherent-symbol bounds on exponential traces (Berezin-Lieb, Golden-Thompson)
Topic: Thermal states and traces
Appears in: P-016

The trace of a Gibbs exponential is bounded below by the integral of $e^{-tQ}$ and above by the integral of $e^{-tP}$; Golden-Thompson bounds the trace of an exponential of a sum by the trace of the product.

$$ \int d\mu(\alpha)e^{-tQ_{\hat H}(\alpha)}\le\operatorname{Tr}e^{-t\hat H}\le\int d\mu(\alpha)e^{-tP_{\hat H}(\alpha)},\qquad \operatorname{Tr}e^{-t(\hat A+\hat B)}\le\operatorname{Tr}(e^{-t\hat A}e^{-t\hat B}) $$

Conditions: self-adjoint, bounded below, trace class, convergent symbol integrals; $d\mu=d^2\alpha/\pi$.
Audit: Correct. Checked by: hand derivation for the lower bound (Jensen: $e^{-t\langle\alpha|\hat H|\alpha\rangle}\le\langle\alpha|e^{-t\hat H}|\alpha\rangle$, then the resolution of identity); upper bound and Golden-Thompson cited as standard theorems.
QF code (exact): none. QF has no operator-inequality or trace-bound function.

### Q and P symbols of a displacement
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: P-017

The covariant symbol of a displacement is a Gaussian-damped plane wave in phase space, and its contravariant symbol has the opposite sign of the $|\gamma|^2/2$ term.

$$ Q_{\hat D(\gamma)}(\alpha)=e^{-|\gamma|^2/2+\gamma\bar\alpha-\bar\gamma\alpha},\qquad P_{\hat D(\gamma)}(\alpha)=e^{+|\gamma|^2/2+\gamma\bar\alpha-\bar\gamma\alpha} $$

Audit: Correct. Checked by: QF exact for $Q$ (vacuum element of $e^{\bar\alpha\hat a}\hat D(\gamma)e^{\alpha\hat a^\dagger}$, times $e^{-|\alpha|^2}$); exact WL Gaussian integral showing that smoothing the stated $P$ gives the stated $Q$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
q = Exp[-al Conjugate[al]] (BosonicMatrixElement[{0, 0}, Exp[cal a] ** DisplacementOperator[g, Infinity] ** Exp[al ad]] /. cal -> Conjugate[al]);
Simplify[q/Exp[-g Conjugate[g]/2 + g Conjugate[al] - Conjugate[g] al]]
```
Output: `1`. Matches: the ratio to $e^{-|\gamma|^2/2+\gamma\bar\alpha-\bar\gamma\alpha}$ is 1.

### Distributional trace of displacement and operator reconstruction
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: P-017, P-020 (the $\hat\rho=\int\chi_W\hat D(-\gamma)$ line)

A displacement has a delta-function trace, displacements are orthogonal in the Hilbert-Schmidt sense, and any operator is rebuilt from its displacement transform.

$$ \operatorname{Tr}\hat D(\gamma)=\pi\delta^{(2)}(\gamma),\quad \operatorname{Tr}[\hat D(\alpha)\hat D(\beta)^\dagger]=\pi\delta^{(2)}(\alpha-\beta),\quad \hat A=\int\frac{d^2\gamma}{\pi}\operatorname{Tr}[\hat D(\gamma)\hat A]\hat D(-\gamma) $$

Conditions: Hilbert-Schmidt / distributional sense; $d^2\gamma=d\,\mathrm{Re}\gamma\,d\,\mathrm{Im}\gamma$.
Audit: Correct. Checked by: hand derivation: $\operatorname{Tr}\hat D(\gamma)=\sum_ne^{-|\gamma|^2/2}L_n(|\gamma|^2)$, and $\sum_nL_n(x)$ is $\pi\delta^{(2)}$ in this measure (Glauber's completeness of displacements).
QF code (exact): none. A distributional trace has no exact QF counterpart.

### Displacement matrix elements in associated Laguerre form
Topic: Fock matrix elements and generating functions
Appears in: P-018, P2-052, P2-011 (transition probability $W_{mn}=|\langle m|\hat D(\gamma)|n\rangle|^2$)

Fock matrix elements of a displacement are associated Laguerre polynomials times a Gaussian and a power of the amplitude.

$$ \langle m|\hat D(\gamma)|n\rangle=e^{-|\gamma|^2/2}\begin{cases}\sqrt{n!/m!}\,\gamma^{m-n}L_n^{(m-n)}(|\gamma|^2),&m\ge n,\\ \sqrt{m!/n!}\,(-\bar\gamma)^{n-m}L_m^{(n-m)}(|\gamma|^2),&n\ge m\end{cases} $$

Audit: Correct. Checked by: QF exact: the matrix elements of the normal-ordered product $e^{-|\gamma|^2/2}e^{\gamma\hat a^\dagger}e^{-\bar\gamma\hat a}$ (a route independent of QF's built-in Laguerre formula) equal the stated Laguerre form for all $m,n\le6$, with $\gamma,\bar\gamma$ symbolic.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
lag[m_, n_] := Exp[-g gb/2] If[m >= n, Sqrt[n!/m!] g^(m - n) LaguerreL[n, m - n, g gb], Sqrt[m!/n!] (-gb)^(n - m) LaguerreL[m, n - m, g gb]];
nrm[m_, n_] := BosonicMatrixElement[{m, n}, DisplacementOperator[g, Infinity, "Ordering" -> "Normal"]] /. Conjugate[g] -> gb;
Union[Flatten[Table[Simplify[nrm[m, n] - lag[m, n]], {m, 0, 6}, {n, 0, 6}]]]
```
Output: `{0}`. Matches: all 49 differences are 0.

### Coherent generating kernel of the displacement
Topic: Fock matrix elements and generating functions
Appears in: P-018, P2-051 (same identity: $\sum_{mn}D_{mn}\bar\alpha^m\beta^n/\sqrt{m!n!}$ is $e^{(|\alpha|^2+|\beta|^2)/2}\langle\alpha|\hat D|\beta\rangle$)

The unnormalized coherent matrix element of a displacement is a single Gaussian exponential in the two coherent labels and $\gamma$.

$$ e^{(|\alpha|^2+|\beta|^2)/2}\langle\alpha|\hat D(\gamma)|\beta\rangle=e^{-|\gamma|^2/2+\bar\alpha\beta+\bar\alpha\gamma-\beta\bar\gamma} $$

Audit: Correct. Checked by: QF exact, $e^{(|\alpha|^2+|\beta|^2)/2}\langle\alpha|\hat D(\gamma)|\beta\rangle=\langle0|e^{\bar\alpha\hat a}\hat D(\gamma)e^{\beta\hat a^\dagger}|0\rangle$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
k = BosonicMatrixElement[{0, 0}, Exp[cal a] ** DisplacementOperator[g, Infinity] ** Exp[be ad]];
Simplify[k/Exp[-g Conjugate[g]/2 + cal be + cal g - be Conjugate[g]]]
```
Output: `1`. Matches: the ratio to the stated kernel is 1 (`cal` stands for $\bar\alpha$).

### Displacement matrix elements as a finite sum
Topic: Fock matrix elements and generating functions
Appears in: P2-051

Each displacement matrix element is a finite sum over the number of contractions, with all phases explicit.

$$ D_{mn}=e^{-|\gamma|^2/2}\sqrt{m!n!}\sum_{r=0}^{\min(m,n)}\frac{\gamma^{m-r}(-\bar\gamma)^{n-r}}{r!(m-r)!(n-r)!} $$

Audit: Correct. Checked by: QF exact, against the matrix elements of the normal-ordered product for all $m,n\le6$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
fs[m_, n_] := Exp[-g gb/2] Sqrt[m! n!] Sum[g^(m - r) (-gb)^(n - r)/(r! (m - r)! (n - r)!), {r, 0, Min[m, n]}];
nrm[m_, n_] := BosonicMatrixElement[{m, n}, DisplacementOperator[g, Infinity, "Ordering" -> "Normal"]] /. Conjugate[g] -> gb;
Union[Flatten[Table[Simplify[nrm[m, n] - fs[m, n]], {m, 0, 6}, {n, 0, 6}]]]
```
Output: `{0}`. Matches.

### Weyl operator from displaced parity, Weyl symbol and reconstruction
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: P-019

The Weyl operator is twice a displaced parity, parity flips the sign of a displacement, the Weyl symbol is the trace against it, and these operators are delta-orthogonal.

$$ \hat T(\alpha)=2\hat D(\alpha)\hat\Pi \hat D(-\alpha)=2\hat D(2\alpha)\hat\Pi,\quad \hat\Pi \hat D(\alpha)\hat\Pi=\hat D(-\alpha),\quad W_{\hat A}(\alpha)=\operatorname{Tr}[\hat T(\alpha)\hat A],\quad \hat A=\int d\mu(\alpha)W_{\hat A}(\alpha)\hat T(\alpha),\quad \pi^{-1}\operatorname{Tr}[\hat T(\alpha)\hat T(\beta)]=\delta^{(2)}(\alpha-\beta) $$

Conditions: $d\mu=d^2\alpha/\pi$, $\hat\Pi=(-1)^{\hat a^\dagger\hat a}$.
Audit: Correct. Checked by: QF exact for $\hat\Pi\hat D(\alpha)\hat\Pi=\hat D(-\alpha)$ and $\hat D(\alpha)^2=\hat D(2\alpha)$, which force $2\hat D(\alpha)\hat\Pi\hat D(-\alpha)=2\hat D(\alpha)\hat D(\alpha)\hat\Pi=2\hat D(2\alpha)\hat\Pi$; hand derivation for the orthogonality ($4\operatorname{Tr}[\hat D(2\alpha-2\beta)]=4\pi\delta^{(2)}(2(\alpha-\beta))=\pi\delta^{(2)}(\alpha-\beta)$).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
d[x_] := DisplacementOperator[x, Infinity];
{BosonicExpSimplify[Exp[I Pi ad ** a] ** d[al] ** Exp[-I Pi ad ** a]], BosonicExpSimplify[d[al] ** d[al]],
 Simplify[ComplexExpand[-2 BosonicMatrixElement[{1, 1}, DisplacementOperator[2 (x + I y)]]]]}
```
Output: `{E^(\[FormalA]1*Conjugate[al] - al*SuperDagger[\[FormalA]1]), E^(-2*\[FormalA]1*Conjugate[al] + 2*al*SuperDagger[\[FormalA]1]), (-2 + 8*x^2 + 8*y^2)/E^(2*(x^2 + y^2))}`. Matches; the third output is the Weyl symbol of $|1\rangle\langle1|$, $-2\langle1|\hat D(2\alpha)|1\rangle=-2e^{-2|\alpha|^2}(1-4|\alpha|^2)$ at $\alpha=x+iy$.

### Gaussian smoothing kernels between P, W and Q (derived normalization correction)
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: P-019

With the measure $d^2\alpha/\pi$, the smoothing kernels from $P$ to $W$ and from $W$ to $Q$ are $2e^{-2|\alpha-\beta|^2}$, and from $P$ to $Q$ it is $e^{-|\alpha-\beta|^2}$.

$$ W=\int\frac{d^2\beta}{\pi}\,2e^{-2|\alpha-\beta|^2}P(\beta),\qquad Q=\int\frac{d^2\beta}{\pi}\,2e^{-2|\alpha-\beta|^2}W(\beta),\qquad Q=\int\frac{d^2\beta}{\pi}e^{-|\alpha-\beta|^2}P(\beta) $$

Audit: Correct (the atlas's corrected form; the printed kernels without the factor 2 give half the result). Checked by: QF exact plus exact WL integral for $W\to Q$ on $|1\rangle$: the QF Weyl symbol smoothed with $2e^{-2|\alpha-\beta|^2}$ is exactly the QF covariant symbol $|\alpha|^2e^{-|\alpha|^2}$. Exact WL computation for $P\to W$ on a thermal state: smoothing $P=e^{-|\beta|^2/\bar n}/\bar n$ gives $e^{-|\alpha|^2/(\bar n+1/2)}/(\bar n+1/2)$, equal to $\operatorname{Tr}[\hat T(\alpha)\hat\rho]=2e^{-2|\alpha|^2}\sum_np_n(-1)^nL_n(4|\alpha|^2)$ summed exactly.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
w1 = ComplexExpand[-2 BosonicMatrixElement[{1, 1}, DisplacementOperator[2 (x + I y)]]];
q1 = ComplexExpand[Abs[Exp[-al Conjugate[al]/2] BosonicMatrixElement[{1, 0}, Exp[al ad]]]^2 /. al -> x0 + I y0];
conv = Integrate[2 Exp[-2 ((x0 - x)^2 + (y0 - y)^2)] w1, {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, Assumptions -> Element[{x0, y0}, Reals]]/Pi;
Simplify[{q1, conv - q1}]
```
Output: `{E^(-x0^2 - y0^2)*(x0^2 + y0^2), 0}`. Matches: the first output is $Q_{|1\rangle\langle1|}$, the second is the smoothed Weyl symbol minus it.

### Ordered characteristic functions differ by Gaussian factors
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: P-020, P2-035

Normal, Weyl and antinormal characteristic functions differ only by Gaussian factors $e^{\pm|\gamma|^2/2}$, because the ordered exponentials differ from $\hat D(\gamma)$ by those scalars.

$$ \chi_W=\operatorname{Tr}[\hat\rho\hat D(\gamma)],\ \chi_N=\operatorname{Tr}[\hat\rho e^{\gamma\hat a^\dagger}e^{-\bar\gamma\hat a}],\ \chi_A=\operatorname{Tr}[\hat\rho e^{-\bar\gamma\hat a}e^{\gamma\hat a^\dagger}];\quad \chi_N=e^{|\gamma|^2/2}\chi_W,\ \chi_A=e^{-|\gamma|^2/2}\chi_W;\quad \chi_k=e^{-\sigma_k|\gamma|^2}\chi_N,\ \sigma_N=0,\ \sigma_O=\tfrac12,\ \sigma_A=1 $$

Audit: Correct. P2-035 is the same identity written relative to $\chi_N$ (source letter $O$ for Weyl, argument $\eta$). Checked by: QF exact, at the operator level, so it holds for every density operator.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
dd = DisplacementOperator[g, Infinity];
Simplify[{BosonicExpSimplify[Exp[g ad] ** Exp[-Conjugate[g] a]]/dd, BosonicExpSimplify[Exp[-Conjugate[g] a] ** Exp[g ad]]/dd}]
```
Output: `{E^((g*Conjugate[g])/2), E^(-1/2*(g*Conjugate[g]))}`. Matches: $e^{\gamma\hat a^\dagger}e^{-\bar\gamma\hat a}=e^{|\gamma|^2/2}\hat D(\gamma)$ and $e^{-\bar\gamma\hat a}e^{\gamma\hat a^\dagger}=e^{-|\gamma|^2/2}\hat D(\gamma)$.

### Normally ordered moments from the normal characteristic function
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: P-020

Normally ordered moments are mixed Wirtinger derivatives of the normal characteristic function at the origin.

$$ \langle(\hat a^\dagger)^m\hat a^n\rangle=(-1)^n\left.\partial_\gamma^m\partial_{\bar\gamma}^n\chi_N(\gamma)\right|_0 $$

Audit: Correct. Checked by: hand derivation; each $\partial_{\bar\gamma}$ brings down $-\hat a$ on the right of $e^{\gamma\hat a^\dagger}e^{-\bar\gamma\hat a}$ and each $\partial_\gamma$ brings down $\hat a^\dagger$ on the left. QF exact on a coherent state: $\chi_N=e^{\gamma\bar\beta-\bar\gamma\beta}$ there, and the derivative formula reproduces QF's $\langle\beta|\hat a^{\dagger2}\hat a|\beta\rangle$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
chiN = Exp[g Conjugate[be] - gb be];
{BosonicMatrixElement[{be, be}, ad ** ad ** a, "Basis" -> "Coherent"] /. Abs[be]^2 -> be Conjugate[be],
 -D[chiN, {g, 2}, {gb, 1}] /. {g -> 0, gb -> 0}}
```
Output: `{be*Conjugate[be]^2, be*Conjugate[be]^2}`. Matches: both entries are $\bar\beta^2\beta$.

### Weyl characteristic function of a number-diagonal state
Topic: Thermal states and traces
Appears in: P-020

For a state diagonal in the number basis, the Weyl characteristic function is a Gaussian times the population-weighted sum of Laguerre polynomials.

$$ \chi_W(\gamma)=e^{-|\gamma|^2/2}\sum_np_nL_n(|\gamma|^2) $$

Audit: Correct. Checked by: QF exact for $\langle n|\hat D(\gamma)|n\rangle=e^{-|\gamma|^2/2}L_n(|\gamma|^2)$ ($n\le8$, from the normal-ordered product; for general $n$ see the Laguerre entry), and exact WL sum for thermal $p_n=(1-x)x^n$, which gives $\chi_W=e^{-(\bar n+1/2)|\gamma|^2}$ with $\bar n=x/(1-x)$, $y=|\gamma|^2$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
diag[n_] := BosonicMatrixElement[{n, n}, DisplacementOperator[g, Infinity, "Ordering" -> "Normal"]] /. Conjugate[g] -> gb;
{Union[Table[Simplify[diag[n] - Exp[-g gb/2] LaguerreL[n, g gb]], {n, 0, 8}]],
 Simplify[Exp[-y/2] Sum[(1 - x) x^n LaguerreL[n, y], {n, 0, Infinity}], 0 < x < 1]}
```
Output: `{{0}, E^(((1 + x)*y)/(2*(-1 + x)))}`. Matches.

### su(2) commutators, ladder action and holomorphic realization
Topic: Lie-group decompositions and generalized coherent states
Appears in: P-024

The spin generators obey the su(2) commutators, act on $|j,m\rangle$ with the standard square-root ladder factors, and are realized on polynomials of degree at most $2j$ by first-order differential operators.

$$ [\hat J_0,\hat J_\pm]=\pm \hat J_\pm,\ [\hat J_+,\hat J_-]=2\hat J_0,\ \hat J_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}|j,m\pm1\rangle;\quad \hat J_-=\partial_z,\ \hat J_+=-z^2\partial_z+2jz,\ \hat J_0=z\partial_z-j $$

Audit: Correct. Checked by: exact WL computation: commutators of the differential operators on a symbolic $f(z)$ vanish identically, and the exact spin matrices for $j=1/2,1,3/2,2$ satisfy the algebra with zero residual.
QF code (exact): none. QF's exact algebra layer is bosonic; it has no abstract su(2) generators.

### SU(2) disentangling
Topic: Lie-group decompositions and generalized coherent states
Appears in: P-025

A compact rotation exponential factors into raising, diagonal and lowering exponentials with $\zeta=e^{i\phi}\tan r$, in either order, on the chart $\cos r>0$.

$$ e^{\xi \hat J_+-\bar\xi \hat J_-}=e^{\zeta \hat J_+}e^{\log(1+|\zeta|^2)\hat J_0}e^{-\bar\zeta \hat J_-}=e^{-\bar\zeta \hat J_-}e^{-\log(1+|\zeta|^2)\hat J_0}e^{\zeta \hat J_+},\quad \xi=re^{i\phi},\ \zeta=e^{i\phi}\tan r $$

Conditions: $\cos r>0$.
Audit: Correct on the stated chart. Checked by: exact WL computation in the faithful 2x2 matrices with symbolic $r,\phi$: both orders equal the left side for $0<r<\pi/2$, and for $\pi/2<r<\pi$ both equal minus the left side. That sign is the central element $-I$, which acts as $(-1)^{2j}$ in spin $j$, so the identity holds in every representation on the chart and fails for half-integer $j$ beyond it.
QF code (exact): none. QF's exact algebra layer has no su(2) generators.

### SU(2) coherent state expansion
Topic: Lie-group decompositions and generalized coherent states
Appears in: P-026

The spin coherent state is the normalized exponential of $\zeta\hat J_+$ on the lowest weight, with binomial square-root amplitudes.

$$ |\zeta;j\rangle=(1+|\zeta|^2)^{-j}e^{\zeta \hat J_+}|j,-j\rangle=(1+|\zeta|^2)^{-j}\sum_{n=0}^{2j}\sqrt{\binom{2j}{n}}\zeta^n|j,-j+n\rangle $$

Audit: Correct. Checked by: exact WL computation with exact spin matrices ($j=3/2$ and 2, symbolic $\zeta$), and hand derivation for general $j$: $\hat J_+^n|j,-j\rangle=\sqrt{n!\,(2j)!/(2j-n)!}\,|j,-j+n\rangle$.
QF code (exact): none. No su(2) layer in QF.

### SU(2) coherent overlap and annihilator
Topic: Lie-group decompositions and generalized coherent states
Appears in: P-026

Spin coherent states overlap by $(1+\bar\zeta\eta)^{2j}$ over the normalizations, and each is annihilated by $\hat J_-+2\zeta\hat J_0-\zeta^2\hat J_+$.

$$ \langle\zeta;j|\eta;j\rangle=\frac{(1+\bar\zeta\eta)^{2j}}{(1+|\zeta|^2)^j(1+|\eta|^2)^j},\qquad (\hat J_-+2\zeta \hat J_0-\zeta^2\hat J_+)|\zeta;j\rangle=0 $$

Audit: Correct. Checked by: exact WL computation ($j=3/2$ and 2, symbolic $\zeta,\eta$, zero residuals); hand derivation of the annihilator: $e^{\zeta\hat J_+}\hat J_-e^{-\zeta\hat J_+}=\hat J_-+2\zeta\hat J_0-\zeta^2\hat J_+$.
QF code (exact): none. No su(2) layer in QF.

### SU(2) resolution of identity
Topic: Lie-group decompositions and generalized coherent states
Appears in: P-026

Spin coherent projectors integrated with the invariant sphere measure resolve the identity.

$$ \frac{2j+1}{\pi}\int\frac{d^2\zeta}{(1+|\zeta|^2)^2}|\zeta;j\rangle\langle\zeta;j|=\hat I $$

Audit: Correct. Checked by: exact WL computation; the diagonal element $(2j+1)\binom{2j}{n}\int_0^\infty t^n(1+t)^{-2j-2}dt$ simplifies to 1 for symbolic $j$ and $n$.
QF code (exact): none. No su(2) layer in QF.

### su(1,1) algebra, Casimir, discrete-series ladder and holomorphic realization
Topic: Lie-group decompositions and generalized coherent states
Appears in: P-028

The su(1,1) generators close with $[\hat K_-,\hat K_+]=2\hat K_0$, the Casimir is $k(k-1)$ in the positive discrete series, and the series has explicit ladder factors and a holomorphic realization.

$$ [\hat K_0,\hat K_\pm]=\pm \hat K_\pm,\ [\hat K_-,\hat K_+]=2\hat K_0,\ \hat C=\hat K_0^2-\tfrac12(\hat K_+\hat K_-+\hat K_-\hat K_+)=k(k-1)\hat I,\ \hat K_+|k,n\rangle=\sqrt{(n+1)(2k+n)}|k,n+1\rangle,\ \hat K_-=\partial_z,\ \hat K_+=z^2\partial_z+2kz,\ \hat K_0=z\partial_z+k $$

Audit: Correct. Checked by: exact WL computation: the differential realization satisfies the commutators and the Casimir on a symbolic $f(z)$, and the ladder factors give $[\hat K_-,\hat K_+]|k,n\rangle=2(k+n)|k,n\rangle$ and Casimir $k(k-1)$ for symbolic $k,n$.
QF code (exact): none. QF has no abstract su(1,1) generators; the oscillator realizations are checked in the single-mode and two-mode realization entries.

### SU(1,1) disentangling (middle-exponent sign corrected)
Topic: Lie-group decompositions and generalized coherent states
Appears in: P-029 (also referenced in P2-058)

A noncompact su(1,1) exponential factors into raising, diagonal and lowering exponentials with $\zeta=e^{i\phi}\tanh r$; the diagonal exponent is $+\log(1-|\zeta|^2)\hat K_0$ in normal order and the opposite in antinormal order.

$$ e^{\xi \hat K_+-\bar\xi \hat K_-}=e^{\zeta \hat K_+}e^{\log(1-|\zeta|^2)\hat K_0}e^{-\bar\zeta \hat K_-}=e^{-\bar\zeta \hat K_-}e^{-\log(1-|\zeta|^2)\hat K_0}e^{\zeta \hat K_+},\quad \zeta=e^{i\phi}\tanh r $$

Audit: Correct as corrected in the atlas; the printed source signs are wrong. Checked by: QF exact in the single-mode realization ($\hat K_+=\hat a^{\dagger2}/2$, $\hat K_-=\hat a^2/2$, $\hat K_0=(\hat a^\dagger\hat a+1/2)/2$), where the corrected middle factor is $(1-\tanh^2r)^{\hat K_0}=(\mathrm{sech}\,r)^{\hat a^\dagger\hat a+1/2}$ and the printed one would be $(\cosh r)^{\hat a^\dagger\hat a+1/2}$; and exact WL computation in the faithful 2x2 matrices ($\hat K_+=E_{12}$, $\hat K_-=-E_{21}$, $\hat K_0=\mathrm{diag}(1/2,-1/2)$, commutators verified), symbolic $r>0$, $\phi$: corrected forms give zero residual, printed forms do not. Both sides agree at $r=0$ and vary continuously, so the 2x2 identity lifts to every positive discrete series.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
kp = ad ** ad/2; km = a ** a/2;
gen = rr Exp[I ph] kp - rr Exp[-I ph] km;
{Simplify[BosonicExpOrder[Exp[gen]], rr > 0], Simplify[BosonicExpOrder[Exp[gen], "Ordering" -> "Antinormal"], rr > 0]}
```
Output: `{E^((E^(I*ph)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]*Tanh[rr])/2)**Sech[rr]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^(-1/2*(\[FormalA]1**\[FormalA]1*Tanh[rr])/E^(I*ph)), E^(-1/2*(\[FormalA]1**\[FormalA]1*Tanh[rr])/E^(I*ph))**Sech[rr]^(1/2 - \[FormalA]1**SuperDagger[\[FormalA]1])**E^((E^(I*ph)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]*Tanh[rr])/2)}`. Matches the corrected signs: normal order has middle factor $\mathrm{sech}(r)^{1/2+\hat a^\dagger\hat a}$; antinormal order has $\mathrm{sech}(r)^{1/2-\hat a\hat a^\dagger}=\cosh(r)^{\hat a^\dagger\hat a+1/2}$.

### SU(1,1) coherent state expansion and normalization
Topic: Lie-group decompositions and generalized coherent states
Appears in: P-030

The su(1,1) coherent state is the normalized exponential of $\zeta\hat K_+$ on the lowest weight, with Pochhammer square-root amplitudes, and it is normalized for $|\zeta|<1$.

$$ |\zeta;k\rangle=(1-|\zeta|^2)^ke^{\zeta \hat K_+}|k,0\rangle=(1-|\zeta|^2)^k\sum_{n\ge0}\sqrt{\frac{\Gamma(n+2k)}{n!\Gamma(2k)}}\zeta^n|k,n\rangle $$

Audit: Correct. Checked by: exact WL computation: the product of ladder factors $\prod_{j<n}\sqrt{(j+1)(2k+j)}/n!$ equals $\sqrt{\Gamma(n+2k)/(n!\Gamma(2k))}$ for symbolic $k,n$, and $\sum_n\Gamma(n+2k)x^n/(n!\Gamma(2k))=(1-x)^{-2k}$ exactly.
QF code (exact): none. No su(1,1) layer in QF.

### SU(1,1) coherent overlap and annihilator
Topic: Lie-group decompositions and generalized coherent states
Appears in: P-030

Two su(1,1) coherent states overlap by $(1-\bar\zeta\eta)^{-2k}$ times the normalizations, and each is annihilated by $\hat K_--2\zeta\hat K_0+\zeta^2\hat K_+$.

$$ \langle\zeta;k|\eta;k\rangle=\frac{(1-|\zeta|^2)^k(1-|\eta|^2)^k}{(1-\bar\zeta\eta)^{2k}},\qquad (\hat K_--2\zeta \hat K_0+\zeta^2\hat K_+)|\zeta;k\rangle=0 $$

Audit: Correct. Checked by: exact WL computation: the overlap series sums to $(1-\bar\zeta\eta)^{-2k}$, and the coefficient of $|k,n\rangle$ in the annihilator equation vanishes for symbolic $k,n$ (neighbouring amplitudes differ by $\zeta\sqrt{(n+2k)/(n+1)}$, so the three terms give $\zeta(n+2k)-2\zeta(k+n)+\zeta n=0$).
QF code (exact): none. No su(1,1) layer in QF.

### SU(1,1) resolution of identity
Topic: Lie-group decompositions and generalized coherent states
Appears in: P-030

su(1,1) coherent projectors integrated with the invariant disk measure resolve the identity only for $k>1/2$.

$$ \frac{2k-1}{\pi}\int_{|\zeta|<1}\frac{d^2\zeta}{(1-|\zeta|^2)^2}|\zeta;k\rangle\langle\zeta;k|=\hat I,\quad k>1/2 $$

Audit: Correct, including the restriction $k>1/2$. Checked by: exact WL computation; the diagonal element reduces to $(2k-1)\Gamma(n+2k)/(n!\Gamma(2k))\int_0^1t^n(1-t)^{2k-2}dt$, which the kernel returns as 1 under the condition $k>1/2$ (the integral diverges otherwise).
QF code (exact): none. No su(1,1) layer in QF.

### Single-mode su(1,1) realization and Casimir
Topic: Single-mode squeezing
Appears in: P-031, P2-012 (parity sectors $k=1/4,3/4$)

Half the squares of the ladder operators and half the shifted number operator realize su(1,1) on one mode, with Casimir $-3/16$, so the even and odd sectors carry $k=1/4$ and $3/4$.

$$ \hat K_+=\tfrac12(\hat a^\dagger)^2,\quad \hat K_-=\tfrac12\hat a^2,\quad \hat K_0=\tfrac12(\hat a^\dagger \hat a+\tfrac12\hat I),\quad \hat C=-\tfrac{3}{16}\hat I $$

Audit: Correct; $k(k-1)=-3/16$ for $k=1/4$ and $3/4$. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
kp = ad ** ad/2; km = a ** a/2; k0 = (ad ** a + 1/2)/2;
BosonicNormalOrder[#, {a, ad}] & /@ {k0 ** kp - kp ** k0 - kp, km ** kp - kp ** km - 2 k0, k0 ** k0 - (kp ** km + km ** kp)/2}
```
Output: `{0, 0, -3/16}`. Matches: $[\hat K_0,\hat K_+]=\hat K_+$, $[\hat K_-,\hat K_+]=2\hat K_0$, Casimir $-3/16$.

### Single-mode squeezing: normal-ordered factorization
Topic: Single-mode squeezing
Appears in: P-031

The one-mode squeeze operator in normal order is a pair-creation exponential, a number-diagonal factor $(\mathrm{sech}\,r)^{\hat a^\dagger\hat a+1/2}$ and a pair-annihilation exponential, so its vacuum amplitude is $(\cosh r)^{-1/2}$.

$$ \hat S(\xi)=e^{\frac12(\bar\xi\hat a^2-\xi\hat a^{\dagger2})}=e^{-\frac12e^{i\theta}\tanh r\,(\hat a^\dagger)^2}e^{-\log\cosh r\,(\hat a^\dagger \hat a+1/2)}e^{\frac12e^{-i\theta}\tanh r\,\hat a^2},\quad \xi=re^{i\theta} $$

Audit: Correct. Checked by: QF exact, both the built-in ordered form and `BosonicExpOrder` of the generator.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
nf = SqueezeOperator[xi, Infinity, "Ordering" -> "Normal"];
{nf, Simplify[BosonicExpOrder[Exp[(rr Exp[-I th] a ** a - rr Exp[I th] ad ** ad)/2]], rr > 0],
 BosonicMatrixElement[{0, 0}, nf]}
```
Output: `{E^(-1/2*(E^(I*Arg[xi])*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]1]*Tanh[Abs[xi]]))**Sech[Abs[xi]]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^((\[FormalA]1**\[FormalA]1*Tanh[Abs[xi]])/(2*E^(I*Arg[xi]))), <the same three factors with xi -> rr E^(I th)>, Sqrt[Sech[Abs[xi]]]}`. Matches: with $\xi=re^{i\theta}$ the two ordered forms coincide and the vacuum element is $\mathrm{sech}(r)^{1/2}$.

### Two-mode su(1,1) realization and Casimir
Topic: Two-mode squeezing
Appears in: P-032, P2-019

The pair operators of two modes realize su(1,1), and the Casimir is fixed by the number difference, so a sector with difference $d$ has $k=(1+|d|)/2$.

$$ \hat K_+=\hat a^\dagger \hat b^\dagger,\quad \hat K_-=\hat a\hat b,\quad \hat K_0=\tfrac12(\hat N_a+\hat N_b+\hat I),\quad \hat C=\tfrac14[(\hat N_a-\hat N_b)^2-\hat I] $$

Audit: Correct; $k(k-1)=(d^2-1)/4$ gives $k=(1+|d|)/2$. Checked by: QF exact.
QF code (exact):
```wl
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
kp = ad1 ** ad2; km = a1 ** a2; n1 = ad1 ** a1; n2 = ad2 ** a2; k0 = (n1 + n2 + 1)/2;
BosonicNormalOrder[#, {a1, ad1, a2, ad2}] & /@ {k0 ** kp - kp ** k0 - kp, km ** kp - kp ** km - 2 k0,
  k0 ** k0 - (kp ** km + km ** kp)/2 - ((n1 - n2) ** (n1 - n2) - 1)/4}
```
Output: `{0, 0, 0}`. Matches: both commutators and the Casimir identity reduce to 0.

### Two-mode squeezing: normal-ordered factorization and pair state
Topic: Two-mode squeezing
Appears in: P-032, P2-019

The two-mode squeeze operator in normal order is a pair-creation exponential, $(\mathrm{sech}\,r)^{\hat N_a+\hat N_b+1}$ and a pair-annihilation exponential; on the vacuum it gives a pair state with geometric pair-number distribution.

$$ e^{r(e^{i\phi}\hat a^\dagger \hat b^\dagger-e^{-i\phi}\hat a\hat b)}=e^{e^{i\phi}\tanh r\,\hat a^\dagger \hat b^\dagger}e^{-\log\cosh r\,(\hat N_a+\hat N_b+\hat I)}e^{-e^{-i\phi}\tanh r\,\hat a\hat b};\qquad |\zeta\rangle=\sqrt{1-|\zeta|^2}\,e^{\zeta\hat a^\dagger\hat b^\dagger}|0,0\rangle,\ W_n=(1-\rho)\rho^n,\ \rho=|\zeta|^2 $$

Audit: Correct. P-032 writes the pair state with $\hat a^\dagger\hat b^\dagger$, P2-019 with $\hat b^\dagger\hat a^\dagger$ (same operator). Checked by: QF exact for the factorization; hand derivation for the pair state: $(\hat a^\dagger\hat b^\dagger)^n|0,0\rangle=n!\,|n,n\rangle$ (QF confirms $n=3$), so $e^{\zeta\hat a^\dagger\hat b^\dagger}|0,0\rangle=\sum_n\zeta^n|n,n\rangle$ with squared norm $(1-|\zeta|^2)^{-1}$ and $W_n=(1-\rho)\rho^n$.
QF code (exact):
```wl
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
{Simplify[BosonicExpOrder[Exp[rr Exp[I ph] ad1 ** ad2 - rr Exp[-I ph] a1 ** a2]], rr > 0],
 BosonicMatrixElement[{{3, 3}, {0, 0}}, (ad1 ** ad2) ** (ad1 ** ad2) ** (ad1 ** ad2)]/3!}
```
Output: `{E^(E^(I*ph)*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]2]*Tanh[rr])**Sech[rr]^(SuperDagger[\[FormalA]1]**\[FormalA]1 + SuperDagger[\[FormalA]2]**\[FormalA]2)**E^(-((\[FormalA]1**\[FormalA]2*Tanh[rr])/E^(I*ph)))*Sech[rr], 1}`. Matches for real $r>0$.

### Finite-mode bosonic Bogoliubov conditions
Topic: Bosonic canonical transformations
Appears in: P-034

A linear canonical map of $N$ bosonic modes preserves the indefinite form $K$, which fixes the block relations among $U$ and $V$, the inverse, and makes $Z=U^{-1}V$ symmetric and a strict contraction.

$$ \hat{\mathbf b}=U\hat{\mathbf a}+V\hat{\mathbf a}^\dagger:\ UU^\dagger-VV^\dagger=I_N,\ UV^T=VU^T,\ MKM^\dagger=K,\ M^{-1}=\begin{pmatrix}U^\dagger&-V^T\\-V^\dagger&U^T\end{pmatrix},\ U^\dagger U-V^T\bar V=I_N,\ U^\dagger V=V^T\bar U,\ Z=U^{-1}V=Z^T,\ I_N-ZZ^\dagger=(U^\dagger U)^{-1} $$

Audit: Correct. Checked by: exact WL computation, $N=2$, with an exact group element obtained as the Cayley transform $(I+L/2)(I-L/2)^{-1}$ of a rational generator $L$ of the P-035 form: all ten relations hold with zero residual. Forcing line for the last one: $ZZ^\dagger=U^{-1}(UU^\dagger-I)U^{-\dagger}=I-(U^\dagger U)^{-1}$.
QF code (exact): none. QF has no function that tests Bogoliubov coefficient matrices.

### Quadratic bosonic generators as sp(2N,R) matrices
Topic: Bosonic canonical transformations
Appears in: P-035

A quadratic bosonic generator corresponds to a block matrix with anti-Hermitian number block and symmetric pairing block.

$$ L=\begin{pmatrix}C&A\\\bar A&\bar C\end{pmatrix},\qquad C^\dagger=-C,\quad A^T=A $$

Audit: Correct. Checked by: hand derivation, $LK+KL^\dagger$ has blocks $C+C^\dagger$ and $A^T-A$, so it vanishes exactly when $C^\dagger=-C$ and $A^T=A$; exact WL computation confirms that block form for a general symbolic $L$.
QF code (exact): none. QF has no sp(2N) coefficient-matrix layer.

### Bosonic quadratic commutator table (E, F, G)
Topic: Bosonic canonical transformations
Appears in: P-035

The bilinears $\hat E,\hat F,\hat G$ close under commutators with the stated structure constants.

$$ [\hat E_{ij},\hat E_{kl}]=\delta_{jk}\hat E_{il}-\delta_{il}\hat E_{kj},\ [\hat E_{ij},\hat F_{kl}]=-\delta_{ik}\hat F_{jl}-\delta_{il}\hat F_{kj},\ [\hat E_{ij},\hat G_{kl}]=\delta_{jk}\hat G_{il}+\delta_{jl}\hat G_{ki},\ [\hat F_{ij},\hat G_{kl}]=\delta_{jk}\hat E_{li}+\delta_{jl}\hat E_{ki}+\delta_{ik}\hat E_{lj}+\delta_{il}\hat E_{kj} $$

Conditions: $\hat E_{ij}=\hat a_i^\dagger\hat a_j+\delta_{ij}/2$, $\hat F_{ij}=\hat a_i\hat a_j$, $\hat G_{ij}=\hat a_i^\dagger\hat a_j^\dagger$; $[\hat F,\hat F]=[\hat G,\hat G]=0$.
Audit: Correct. Checked by: QF exact, all 16 index choices over two modes, six relations each.
QF code (exact):
```wl
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
an = {a1, a2}; cr = {ad1, ad2}; d = KroneckerDelta; cm[x_, y_] := x ** y - y ** x;
e[i_, j_] := cr[[i]] ** an[[j]] + d[i, j]/2; f[i_, j_] := an[[i]] ** an[[j]]; g[i_, j_] := cr[[i]] ** cr[[j]];
res = Flatten[Table[BosonicNormalOrder[#, {a1, ad1, a2, ad2}] & /@ {
     cm[e[i, j], e[k, l]] - (d[j, k] e[i, l] - d[i, l] e[k, j]),
     cm[e[i, j], f[k, l]] - (-d[i, k] f[j, l] - d[i, l] f[k, j]),
     cm[e[i, j], g[k, l]] - (d[j, k] g[i, l] + d[j, l] g[k, i]),
     cm[f[i, j], g[k, l]] - (d[j, k] e[l, i] + d[j, l] e[k, i] + d[i, k] e[l, j] + d[i, l] e[k, j]),
     cm[f[i, j], f[k, l]], cm[g[i, j], g[k, l]]}, {i, 2}, {j, 2}, {k, 2}, {l, 2}]];
{Length[res], Union[Expand[res]]}
```
Output: `{96, {0}}`. Matches: 96 results, all 0.

### Bosonic squeezed vacuum chart and its annihilator
Topic: Gaussian vacua and overlaps
Appears in: P-036

A general bosonic group operator factors into pair-creation, number-conserving and pair-annihilation exponentials; on the vacuum this leaves a pair-creation Gaussian, which is annihilated by $\hat{\mathbf a}+Z\hat{\mathbf a}^\dagger$.

$$ \hat T(g)=\mathcal N\exp[-\tfrac12\xi_{ij}\hat X^{ij}]\exp[\alpha^k{}_l\hat X^l{}_k]\exp[-\tfrac12\eta^{ij}\hat X_{ij}];\qquad |Z\rangle=\mathcal N_B(Z)e^{-\frac12\hat{\mathbf a}^{\dagger T}Z\hat{\mathbf a}^\dagger}|0\rangle,\quad(\hat{\mathbf a}+Z\hat{\mathbf a}^\dagger)|Z\rangle=0 $$

Audit: Correct for the state and annihilator; the triangular factorization is structural (no coefficient functions given) and Unchecked beyond its single-mode instance. Checked by: QF exact, one mode and two modes with symmetric $Z$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
gen2 = -(z11 ad1 ** ad1 + 2 z12 ad1 ** ad2 + z22 ad2 ** ad2)/2;
{BosonicExpSimplify[a ** Exp[-z/2 ad ** ad]], BosonicExpSimplify[a1 ** Exp[gen2]], BosonicExpSimplify[a2 ** Exp[gen2]]}
```
Output: `{E^(-1/2*(z*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]1]))**(\[FormalA]1 - z*SuperDagger[\[FormalA]1]), E^(<gen2>)**(\[FormalA]1 - z11*SuperDagger[\[FormalA]1] - z12*SuperDagger[\[FormalA]2]), E^(<gen2>)**(\[FormalA]2 - z12*SuperDagger[\[FormalA]1] - z22*SuperDagger[\[FormalA]2])}`. Matches: $\hat a_i e^{-\frac12\hat{\mathbf a}^{\dagger T}Z\hat{\mathbf a}^\dagger}=e^{-\frac12\hat{\mathbf a}^{\dagger T}Z\hat{\mathbf a}^\dagger}(\hat a_i-\sum_jZ_{ij}\hat a_j^\dagger)$, so on the vacuum $(\hat a_i+\sum_jZ_{ij}\hat a_j^\dagger)|Z\rangle=0$.

### Bosonic Gaussian normalization and overlap
Topic: Gaussian vacua and overlaps
Appears in: P-037

A pair-creation Gaussian is normalized by $\det(I-ZZ^\dagger)^{1/4}$, and two of them overlap by the normalizations over $\det(I-WZ^\dagger)^{1/2}$.

$$ \mathcal N_B(Z)=\det(I_N-ZZ^\dagger)^{1/4},\qquad \langle Z|W\rangle=\frac{[\det(I_N-ZZ^\dagger)\det(I_N-WW^\dagger)]^{1/4}}{\det(I_N-WZ^\dagger)^{1/2}} $$

Conditions: $Z^T=Z$, $I_N-ZZ^\dagger>0$; square-root branch continued from $Z=W=0$.
Audit: Correct. Checked by: QF exact: the one-mode vacuum element in closed form, and for two modes the vacuum element $\langle0|e^{-\frac12\hat{\mathbf a}^T\bar Z\hat{\mathbf a}}e^{-\frac12\hat{\mathbf a}^{\dagger T}W\hat{\mathbf a}^\dagger}|0\rangle$ expanded in powers of $\bar Z,W$ matches the series of $\det(I-W\bar Z)^{-1/2}$ through third order with symbolic symmetric $Z,W$. All orders follow from the Bargmann Gaussian integral $\int\prod d^2w/\pi\,e^{-|w|^2-\frac12\bar w^T\bar Z\bar w-\frac12w^TWw}=\det(I-W\bar Z)^{-1/2}$ (standard).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
one = BosonicMatrixElement[{0, 0}, Exp[-Conjugate[z]/2 a ** a] ** Exp[-w/2 ad ** ad]];
zq = (zb11 a1 ** a1 + 2 zb12 a1 ** a2 + zb22 a2 ** a2)/2;
wq = (w11 ad1 ** ad1 + 2 w12 ad1 ** ad2 + w22 ad2 ** ad2)/2;
ov[k_] := BosonicMatrixElement[{{0, 0}, {0, 0}}, NonCommutativeMultiply @@ Join[Table[zq, k], Table[wq, k]]]/(k!)^2;
target = Series[Det[IdentityMatrix[2] - s {{w11, w12}, {w12, w22}} . {{zb11, zb12}, {zb12, zb22}}]^(-1/2), {s, 0, 3}];
{one, Expand[Table[ov[k] - SeriesCoefficient[target, k], {k, 1, 3}]]}
```
Output: `{1/Sqrt[1 - w*Conjugate[z]], {0, 0, 0}}`. Matches: one mode gives $(1-\bar zw)^{-1/2}$; the two-mode series differences are 0 at orders 1, 2, 3.

### Fractional-linear action on symmetric Gaussian parameters
Topic: Bosonic canonical transformations
Appears in: P-038

Canonical maps act on the symmetric Gaussian parameter by a matrix fractional-linear map, and the analytic representation carries a determinant multiplier.

$$ Z\mapsto (U+Z\bar V)^{-1}(V+Z\bar U),\qquad [\hat T^{(+)}(g)f](Z)=\det(U+Z\bar V)^{-1/2}f\left((U+Z\bar V)^{-1}(V+Z\bar U)\right) $$

Audit: Correct for the map; the determinant multiplier is Unchecked (its sign and branch depend on the source's function-space convention). Checked by: exact WL computation ($N=2$, two exact Cayley group elements, rational $Z$): the image is symmetric, $I-Z'Z'^\dagger$ has positive eigenvalues, and $g_1\cdot(g_2\cdot Z)=(g_2g_1)\cdot Z$ with zero residual, which is the composition order that makes $f\mapsto f(g\cdot Z)$ a representation.
QF code (exact): none. QF has no Gaussian-parameter action.

### Bosonic Gaussian bilinear symbols and purity constraints
Topic: Gaussian vacua and overlaps
Appears in: P-039

On a pure bosonic Gaussian the pairing matrix $A$ and the number matrix $B$ are rational functions of $Z$ and obey quadratic purity constraints.

$$ A=-\hbar Z(I_N-Z^\dagger Z)^{-1},\quad B=\tfrac\hbar2 I_N+\hbar Z(I_N-Z^\dagger Z)^{-1}Z^\dagger,\quad B^2-A\bar A=\tfrac{\hbar^2}{4}I_N,\quad BA=A\bar B,\quad B^\dagger=B,\quad A^T=A $$

Audit: Convention: on the two-mode state $|Z\rangle$ ($\hbar=1$), $A_{ij}=\langle\hat a_i\hat a_j\rangle$ and $B_{ij}=\langle\hat a_j^\dagger\hat a_i\rangle+\delta_{ij}/2$; reading $B_{ij}$ as $\langle\hat a_i^\dagger\hat a_j\rangle+\delta_{ij}/2$ (that is, $\langle\hat X^i{}_j\rangle$) gives the transpose, which differs. The four matrix constraints are Correct. Checked by: exact WL computation with an exact rational complex symmetric $Z$: the second moments are solved exactly from the annihilator relation $\hat a_i|Z\rangle=-\sum_kZ_{ik}\hat a_k^\dagger|Z\rangle$ (QF exact in the previous entry) and the CCR, then compared with the formulas; the constraints hold with zero residual.
QF code (exact): none. `BosonicMatrixElement` takes polynomials between Fock states; it does not give moments of a multimode Gaussian state, and `CovarianceMatrix` works only on cut-off states.

### Fermionic Bogoliubov constraints
Topic: Fermionic canonical transformations
Appears in: P-041

A fermionic canonical map keeps the anticommutators exactly when the doubled coefficient matrix is unitary, and $Z=U^{-1}V$ is then skew-symmetric.

$$ \{\hat c_i,\hat c_j^\dagger\}=\delta_{ij},\ \hat{\mathbf b}=U\hat{\mathbf c}+V\hat{\mathbf c}^\dagger:\ UU^\dagger+VV^\dagger=I_N,\ UV^T+VU^T=0,\ M^\dagger M=I_{2N},\ Z=U^{-1}V=-Z^T $$

Audit: Correct. Checked by: exact fermion matrices (3 modes), with an exact unitary group element from the Cayley transform of a rational anti-Hermitian generator of block form: the new operators satisfy the CAR exactly and all block relations hold with zero residual.
QF code (exact): none. no fermionic layer.

### Fermionic quadratic commutator table
Topic: Fermionic canonical transformations
Appears in: P-042

Fermionic bilinears close under commutators with structure constants that differ from the bosonic table in the $[\hat F,\hat G]$ signs.

$$ [\hat F_{ij},\hat G_{kl}]=-\delta_{jk}\hat E_{li}+\delta_{ik}\hat E_{lj}+\delta_{jl}\hat E_{ki}-\delta_{il}\hat E_{kj},\quad [\hat F,\hat F]=[\hat G,\hat G]=0,\ \text{first three bosonic relations unchanged} $$

Conditions: $\hat E_{ij}=\hat c_i^\dagger\hat c_j-\delta_{ij}/2$, $\hat F_{ij}=\hat c_i\hat c_j$, $\hat G_{ij}=\hat c_i^\dagger\hat c_j^\dagger$.
Audit: Correct. Checked by: exact fermion matrices (2 modes), all 16 index choices, six relations: every residual exactly 0.
QF code (exact): none. no fermionic layer.

### Fermionic pair-exponential state, normalization and overlap
Topic: Gaussian vacua and overlaps
Appears in: P-043

A fermionic pair-exponential state with skew $Z$ is normalized by $\det(I+ZZ^\dagger)^{-1/4}$, is annihilated by $\hat{\mathbf c}+Z\hat{\mathbf c}^\dagger$, and two such states overlap by $\det(I+WZ^\dagger)^{1/2}$ over the normalizations.

$$ |Z\rangle=\det(I_N+ZZ^\dagger)^{-1/4}e^{-\frac12\hat{\mathbf c}^{\dagger T}Z\hat{\mathbf c}^\dagger}|0\rangle,\ Z^T=-Z,\ (\hat{\mathbf c}+Z\hat{\mathbf c}^\dagger)|Z\rangle=0,\quad \langle Z|W\rangle=\frac{\det(I_N+WZ^\dagger)^{1/2}}{[\det(I_N+ZZ^\dagger)\det(I_N+WW^\dagger)]^{1/4}} $$

Audit: Correct. Checked by: exact fermion matrices (4 modes), exact rational complex skew $Z,W$: the squared norm squared equals $\det(I+ZZ^\dagger)$, the annihilator residual is exactly 0, and the unnormalized overlap equals the principal square root of $\det(I+WZ^\dagger)$ exactly (the atlas's branch caveat stands for other parameters).
QF code (exact): none. no fermionic layer.

### Fermionic analytic group action
Topic: Fermionic canonical transformations
Appears in: P-044

Fermionic canonical maps act on skew $Z$ by the same fractional-linear map as bosons, with determinant exponent $+1/2$.

$$ [\hat T^{(+)}(g)f](Z)=\det(U+Z\bar V)^{1/2}f\left((U+Z\bar V)^{-1}(V+Z\bar U)\right) $$

Audit: Unchecked: the determinant exponent depends on the source's function-space realization, which is not reconstructed here. The skew-symmetry of the image follows by the same algebra as the bosonic map with signs from $UV^T=-VU^T$.
QF code (exact): none. no fermionic layer.

### Fermionic bilinear symbols and purity constraints
Topic: Gaussian vacua and overlaps
Appears in: P-045

On a pure fermionic Gaussian the skew pairing matrix $A$ and the number matrix $B$ are rational in $Z$ and obey the same quadratic constraints as for bosons with the signs changed.

$$ A=-\hbar Z(I_N+Z^\dagger Z)^{-1},\quad B=\hbar Z(I_N+Z^\dagger Z)^{-1}Z^\dagger-\tfrac\hbar2I_N,\quad B^2-A\bar A=\tfrac{\hbar^2}{4}I_N,\quad BA=A\bar B,\quad B^\dagger=B,\quad A^T=-A $$

Audit: Convention: on the four-mode state of the pair-exponential entry ($\hbar=1$), $A_{ij}=\langle\hat c_j\hat c_i\rangle=-\langle\hat c_i\hat c_j\rangle$ and $B_{ij}=\langle\hat c_j^\dagger\hat c_i\rangle-\delta_{ij}/2$; read as $\langle\hat X_{ij}\rangle$ with $\hat X_{ij}=\hat c_i\hat c_j$, the sign of $A$ is wrong. The constraints are Correct. Checked by: exact fermion matrices (4 modes), exact rational $Z$, all residuals exactly 0.
QF code (exact): none. no fermionic layer.

### Root-space coherent states and displacement-like generator
Topic: Lie-group decompositions and generalized coherent states
Appears in: P-047

For a general semisimple group, coherent states are exponentials of positive-root generators on a lowest weight, and the unitary displacement uses positive and negative roots.

$$ |\boldsymbol\zeta\rangle=\mathcal N(\boldsymbol\zeta)\exp\left(\sum_{\alpha>0}\zeta_\alpha \hat E_\alpha\right)|-\lambda\rangle,\qquad \hat D(\boldsymbol\xi)=\exp\left(\sum_{\alpha>0}\xi_\alpha \hat E_\alpha-\bar\xi_\alpha \hat E_{-\alpha}\right) $$

Audit: Unchecked: definitions on a general semisimple group; the conversion from $\xi$ to $\zeta$ is not given. Rank-one instances are the SU(2) and SU(1,1) disentangling entries.
QF code (exact): none. QF has no root-space layer.

### Determinant and log-determinant series
Topic: Gaussian integrals and Pfaffians
Appears in: P-048

The determinant is the exponential of the trace of the logarithm, and $-\log\det(I-X)$ is the trace power series, convergent for spectral radius below 1.

$$ \det A=e^{\operatorname{Tr}\log A},\qquad -\log\det(I_N-X)=\sum_{m\ge1}\frac{\operatorname{Tr}X^m}{m},\ \rho(X)<1;\qquad \psi_Z(W)=B(Z,Z)^{-1/2}B(Z,W) $$

Audit: Correct (the kernel $\psi_Z$ is a definition). Checked by: exact WL computation: for a fully symbolic $3\times3$ matrix the series of $-\log\det(I-tX)$ agrees with $\sum_mt^m\operatorname{Tr}X^m/m$ through $t^5$.
QF code (exact): none. Not an operator identity.

### Siegel upper half-space and symplectic fractional-linear action
Topic: Bosonic canonical transformations
Appears in: P2-002

Real symplectic matrices act on the Siegel upper half-space by a fractional-linear map, and their blocks obey three quadratic relations.

$$ \mathfrak H_p=\{Z=X+iY:Z^T=Z,\ Y>0\},\quad Z\mapsto(AZ+B)(CZ+D)^{-1},\quad A^TC=C^TA,\ B^TD=D^TB,\ A^TD-C^TB=I_p $$

Audit: Correct. Checked by: exact WL computation: an exact rational real symplectic $4\times4$ matrix (Cayley transform of $JS$, $S$ symmetric) satisfies the block relations with zero residual, and the image of a point of the half-space is symmetric with positive-definite imaginary part.
QF code (exact): none. Not an operator identity.

### Off-diagonal covariant symbol by analytic continuation
Topic: Normal symbols and normal-order exponentials
Appears in: P2-005

The covariant symbol is the coherent expectation value divided by the overlap, and its off-diagonal version uses two different coherent labels.

$$ A(z,\bar z)=\frac{\langle z|\hat A|z\rangle}{\langle z|z\rangle},\qquad A(z,\bar w)=\frac{\langle z|\hat A|w\rangle}{\langle z|w\rangle} $$

Audit: Correct (definitions; the continuation statement holds because the off-diagonal symbol is holomorphic in one label and antiholomorphic in the other). Checked by: hand derivation; QF exact for a flat-space instance.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
Simplify[BosonicMatrixElement[{al, be}, ad ** a, "Basis" -> "Coherent"]/BosonicMatrixElement[{al, be}, 1, "Basis" -> "Coherent"]]
```
Output: `be*Conjugate[al]`. Matches: the flat off-diagonal symbol of $\hat a^\dagger\hat a$ is $\bar\alpha\beta$, the analytic continuation of $|\alpha|^2$; QF has no disk coherent states.

### Disk star product to first order in h
Topic: Normal symbols and normal-order exponentials
Appears in: P2-006

On the disk, the Berezin transform is the identity plus $h$ times the invariant Laplacian, and the symbol product differs from the pointwise product at order $h$ by a derivative pairing.

$$ \hat T_hf=f+h\hat\Delta f+o(h),\ \hat\Delta=(1-|z|^2)^2\partial_z\partial_{\bar z},\qquad A_1\star A_2=A_1A_2+h(1-|z|^2)^2(\partial_{\bar z}A_1)(\partial_zA_2)+o(h) $$

Audit: Unchecked: asymptotic statement not re-derived. The pairing $\partial_{\bar z}A_1\,\partial_zA_2$ is the one expected when $z$ labels the creation side.
QF code (exact): none. QF has no disk coherent states.

### Compact Kahler kernel
Topic: Lie-group decompositions and generalized coherent states
Appears in: P2-007

On the sphere the reproducing kernel is $(1+z\bar w)^{1/h}$ with $1/h$ an integer.

$$ K(z,\bar w)=(1+z\bar w)^{1/h},\quad h^{-1}\in\mathbb Z $$

Audit: Correct: with $1/h=2j$ this is the unnormalized SU(2) coherent overlap. Checked by: comparison with the SU(2) overlap check (exact).
QF code (exact): none. No su(2) layer in QF.

### Forced oscillator: displacement plus accumulated phase
Topic: Forced and parametric oscillators
Appears in: P2-010

A linearly forced oscillator stays coherent, the amplitude integrates the force, a phase accumulates as $\dot\phi=\mathrm{Im}(\bar\beta\alpha)$, and from the vacuum the number distribution is Poisson.

$$ |\psi(t)\rangle=e^{-i\phi(t)}|\alpha(t)\rangle,\quad \dot\alpha=\beta(t),\quad \dot\phi=\operatorname{Im}(\bar\beta\alpha),\quad W_n=e^{-|\gamma|^2}|\gamma|^{2n}/n! $$

Conditions: interaction picture, generator $\beta\hat a^\dagger-\bar\beta\hat a$ with no extra c-number term.
Audit: Correct. Checked by: QF exact for the step $\hat D(\beta\,dt)\hat D(\alpha)=e^{dt(\beta\bar\alpha-\bar\beta\alpha)/2}\hat D(\alpha+\beta\,dt)$, so $\dot\phi=-\mathrm{Im}(\beta\bar\alpha)=\mathrm{Im}(\bar\beta\alpha)$, and for the amplitude $\langle n|\hat D(\gamma)|0\rangle=e^{-|\gamma|^2/2}\gamma^n/\sqrt{n!}$ (symbolic $n$), whose square is the Poisson law.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
d[x_] := DisplacementOperator[x, Infinity] /. Conjugate[s_Plus] :> Map[Conjugate, s];
{Simplify[BosonicExpSimplify[d[b dt] ** d[al]]/d[al + b dt], Element[dt, Reals]],
 Exp[-g Conjugate[g]/2] BosonicMatrixElement[{n, 0}, Exp[g ad]]}
```
Output: `{E^((dt*(b*Conjugate[al] - al*Conjugate[b]))/2), (Sqrt[n!]*Piecewise[{{g^n/n!, n >= 0}}, 0])/E^((g*Conjugate[g])/2)}`. Matches.

### Parametric oscillator su(1,1) generators
Topic: Forced and parametric oscillators
Appears in: P2-012

The quadratic oscillator operators $\hat K_0$ and $\hat K_2$ realize su(1,1); $\hat K_0$ is half the shifted number operator.

$$ \hat K_0=\tfrac14(\hat p^2/\omega_0+\omega_0\hat q^2),\qquad \hat K_2=\tfrac14(\hat p\hat q+\hat q\hat p) $$

Audit: Convention: $\hat K_0=\tfrac12(\hat a^\dagger\hat a+\tfrac12)$ is Correct; $\hat K_2=\tfrac{i}{4}(\hat a^{\dagger2}-\hat a^2)=-(\hat K_+-\hat K_-)/(2i)$ for the single-mode $\hat K_\pm$, and the source's $\hat K_1=\tfrac14(\hat p^2/\omega_0-\omega_0\hat q^2)=-(\hat a^{\dagger2}+\hat a^2)/4$ flips sign too, so this basis matches the single-mode realization after $\hat a\to i\hat a$; the algebra closes either way. Checked by: QF exact, with $\hat q=(\hat a+\hat a^\dagger)/\sqrt{2\omega_0}$, $\hat p=-i\sqrt{\omega_0/2}(\hat a-\hat a^\dagger)$, $\hbar=1$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2 w0]; p = -I Sqrt[w0/2] (a - ad);
BosonicNormalOrder[#, {a, ad}] & /@ {(p ** p/w0 + w0 q ** q)/4, (p ** q + q ** p)/4, (p ** p/w0 - w0 q ** q)/4}
```
Output: `{1/4 + SuperDagger[\[FormalA]1]**\[FormalA]1/2, (-1/4*I)*GeneralizedPower[NonCommutativeMultiply, \[FormalA]1, 2] + (I/4)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2], -1/4*GeneralizedPower[NonCommutativeMultiply, \[FormalA]1, 2] - GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]/4}`. Matches: the outputs are $\hat K_0$, $\hat K_2$ and $\hat K_1$ in ladder form.

### Riccati equation and phase for squeezing dynamics (phase sign corrected)
Topic: Quadratic propagators and Riccati
Appears in: P2-013

Under a Hamiltonian linear in the su(1,1) generators, a coherent state stays coherent; its label obeys a Riccati equation and its phase obeys a first-order equation.

$$ \hat H=i(\beta \hat K_+-\bar\beta \hat K_--i\gamma \hat K_0),\ |\psi\rangle=e^{-i\phi}|\zeta;k\rangle:\quad \dot\zeta=\beta-\bar\beta\zeta^2-i\gamma\zeta,\qquad i\dot\phi=k(\bar\beta\zeta-\beta\bar\zeta+i\gamma) $$

Audit: Correct as corrected in the atlas; the printed phase sign is wrong. Checked by: exact WL computation in the holomorphic realization, $\psi=e^{-i\phi}(1-\zeta\bar\zeta)^k(1-\zeta z)^{-2k}$ with $\zeta,\bar\zeta,\phi$ functions of $t$: the Schrodinger residual divided by $\psi$ is exactly 0 with the corrected phase, and $2ik(\bar\beta\zeta-\beta\bar\zeta)=-4k\,\mathrm{Im}(\bar\beta\zeta)$ with the printed one.
QF code (exact): none. No su(1,1) layer in QF.

### Floquet classification of a quadratic oscillator
Topic: Quadratic propagators and Riccati
Appears in: P2-014

The sign of the Lorentz norm of the generator classifies a periodic quadratic oscillator: discrete spectrum, two-sided continuum, or one-sided continuum.

$$ \Omega_0^2-\Omega_1^2-\Omega_2^2>0:\ \text{spectrum } \pm\hbar\Omega(k+n);\quad <0:\ \text{two-sided continuous};\quad =0,\ \Omega\ne0:\ \text{one-sided continuous} $$

Audit: Correct as a statement about Lie-algebra generators (timelike elements are conjugate to $\Omega\hat K_0$, whose spectrum is $\Omega(k+n)$); the existence of a real logarithm of the period map is a separate assumption, as the atlas says. Checked by: hand derivation for the timelike case; other cases Unchecked.
QF code (exact): none. QF has no Floquet classification.

### Singular oscillator ground state and su(1,1) label
Topic: Forced and parametric oscillators
Appears in: P2-015

The half-line inverse-square oscillator has ground state $x^ae^{-x^2/2}$ with $a(a-1)=2g^2$, energy $a+1/2$, and an equally spaced ladder with su(1,1) label $k=E_0/2$.

$$ \psi_0\sim x^a e^{-x^2/2},\quad a(a-1)=2g^2,\quad E_0=a+\tfrac12,\quad k=\tfrac12(a+\tfrac12),\quad E_n=2n+E_0 $$

Conditions: $\hat H=\tfrac12(\hat p^2+x^2)+g^2/x^2$ on the half line, regular boundary condition.
Audit: Correct. Checked by: exact WL computation, $\hat H(x^ae^{-x^2/2})/(x^ae^{-x^2/2})$ simplifies to $a+1/2$ when the inverse-square coefficient is $a(a-1)/2$; $k=E_0/2$ since $\hat K_0=\hat H/2$.
QF code (exact): none. Not a Fock-space operator identity.

### Calogero degeneracy generating function
Topic: Infinite modes and domains
Appears in: P2-016

The degeneracies of the Calogero model after removing the center of mass are generated by $\prod_{j=2}^N(1-t^j)^{-1}$.

$$ \prod_{j=2}^N(1-t^j)^{-1} $$

Audit: Correct by counting: symmetric polynomials in $N$ variables are generated by power sums of degree $1..N$; removing center-of-mass motion removes degree 1. Checked by: hand derivation.
QF code (exact): none. Not an operator identity.

### Affine Gaussian state normalization (derived correction)
Topic: Gaussian vacua and overlaps
Appears in: P2-017

The state $e^{\zeta\hat a^{\dagger2}/2+\alpha\hat a^\dagger}|0\rangle$ is normalized by $(1-|\zeta|^2)^{1/4}$ times a Gaussian in $\alpha$ that includes the cross term $\mathrm{Re}(\bar\zeta\alpha^2)$.

$$ |\zeta,\alpha\rangle=Ne^{\zeta\hat a^{\dagger2}/2+\alpha\hat a^\dagger}|0\rangle,\qquad N=(1-|\zeta|^2)^{1/4}\exp\left[-\frac{|\alpha|^2+\operatorname{Re}(\bar\zeta\alpha^2)}{2(1-|\zeta|^2)}\right] $$

Audit: Correct as corrected; the printed $N=(1-|\zeta|^2)^{-1/2}e^{-|\alpha|^2/2}$ is wrong. Checked by: exact WL computation: the squared norm, the Bargmann Gaussian integral $\int d^2w/\pi\,e^{-|w|^2}|e^{\zeta w^2/2+\alpha w}|^2$ with symbolic complex $\zeta,\alpha$, times the corrected $N^2$ is exactly 1; with the printed $N$ it is $(1-\zeta^2)^{-3/2}e^{\alpha^2\zeta/(1-\zeta)}\ne1$ already for real $\zeta,\alpha$.
QF code (exact): none. `BosonicMatrixElement[{0, 0}, Exp[ca a] ** Exp[cz/2 a ** a] ** Exp[z/2 ad ** ad] ** Exp[al ad]]` returns unevaluated, and `BosonicExpOrder` of that product orders the quadratic factors but leaves the linear wings unordered.

### Forced parametric oscillator generating kernel and vacuum amplitude
Topic: Forced and parametric oscillators
Appears in: P2-018

For a forced parametric oscillator the vacuum amplitude has modulus $(1-\rho)^{1/4}e^{-\nu(1-\sqrt\rho\cos2\varphi)/2}$.

$$ |\langle0|\psi\rangle|=(1-\rho)^{1/4}e^{-\nu(1-\sqrt\rho\cos2\varphi)/2} $$

Audit: Unchecked: the map from the force and frequency history to $(\rho,\nu,\varphi)$ is a source lookup. At $\nu=0$ it reduces to the squeezed-vacuum amplitude $(1-|\zeta|^2)^{1/4}$ of the affine Gaussian entry, consistent with $\rho=|\zeta|^2$.
QF code (exact): none. The source parameters are not reconstructed.

### Multimode particle-antiparticle bosonic blocks
Topic: Bosonic canonical transformations
Appears in: P2-020

A canonical map mixing particles and antiparticles preserves the indefinite form $E$, which fixes the block relations and the inverse.

$$ MEM^\dagger=E,\quad AA^\dagger-BB^\dagger=I_N,\quad DD^\dagger-CC^\dagger=I_N,\quad AC^\dagger=BD^\dagger,\quad M^{-1}=EM^\dagger E $$

Audit: Correct. Checked by: hand derivation, blockwise expansion of $MEM^\dagger=E$; same structure as the finite-mode bosonic Bogoliubov entry with $(\hat a,\hat b^\dagger)$ in place of $(\hat a,\hat a^\dagger)$.
QF code (exact): none. QF has no coefficient-matrix layer.

### Multimode paired vacuum normalization
Topic: Two-mode squeezing
Appears in: P2-021

The multimode pair-creation exponential on the vacuum has squared norm $\det(I_N-Z^\dagger Z)^{-1}$, so the stated prefactor normalizes it.

$$ |Z\rangle=\det(I_N-Z^\dagger Z)^{1/2}\exp\left(\sum_{ij}Z_{ij}\hat a_i^\dagger \hat b_j^\dagger\right)|0,0\rangle,\quad I_N-Z^\dagger Z>0 $$

Audit: Correct (the relation $Z=(A^\dagger)^{-1}C^\dagger=BD^{-1}$ is convention-dependent and Unchecked). Checked by: QF exact, $N=2$ (four modes, symbolic complex $Z$): the vacuum element of $e^{\sum\bar Z_{ij}\hat a_i\hat b_j}e^{\sum Z_{ij}\hat a_i^\dagger\hat b_j^\dagger}$ expanded in powers of $Z$ matches $\det(I-Z^\dagger Z)^{-1}$ through third order; hand derivation for all orders: in Bargmann variables the $v$ integral of $e^{-|v|^2+u^TZv+\overline{u^TZv}}$ gives $e^{u^TZZ^\dagger\bar u}$, and the $u$ integral then gives $\det(I-Z^\dagger Z)^{-1}$ for $I-Z^\dagger Z>0$.
QF code (exact):
```wl
{a1, a2, b1, b2} = FieldVariables[{1, 2, 3, 4}][[{1, 3, 5, 7}]];
{ad1, ad2, bd1, bd2} = SuperDagger /@ {a1, a2, b1, b2};
zm = {{z11, z12}, {z21, z22}};
cre = Total[Flatten[zm {{ad1 ** bd1, ad1 ** bd2}, {ad2 ** bd1, ad2 ** bd2}}]];
ann = Total[Flatten[Conjugate[zm] {{a1 ** b1, a1 ** b2}, {a2 ** b1, a2 ** b2}}]];
nrm[k_] := BosonicMatrixElement[{{0, 0, 0, 0}, {0, 0, 0, 0}}, NonCommutativeMultiply @@ Join[Table[ann, k], Table[cre, k]],
    "Generators" -> {a1, a2, b1, b2}]/(k!)^2;
target = Series[1/Det[IdentityMatrix[2] - s ConjugateTranspose[zm] . zm], {s, 0, 3}];
Expand[Table[nrm[k] - SeriesCoefficient[target, k], {k, 1, 3}]]
```
Output: `{0, 0, 0}`. Matches: the differences at orders 1, 2, 3 are 0.

### Multimode paired vacuum overlap
Topic: Two-mode squeezing
Appears in: P2-022

Two multimode paired vacua overlap by the product of square-root normalizations over $\det(I_N-Z_1^\dagger Z_2)$.

$$ \langle Z_1|Z_2\rangle=\frac{\det(I_N-Z_1^\dagger Z_1)^{1/2}\det(I_N-Z_2^\dagger Z_2)^{1/2}}{\det(I_N-Z_1^\dagger Z_2)},\qquad Z'=(AZ-B)(-CZ+D)^{-1} $$

Audit: Correct for the overlap; the group action is Unchecked (convention-dependent orientation, as the atlas notes). Checked by: QF exact, $N=2$ with independent symbolic $Z_1,Z_2$, through third order; all orders by the same Bargmann integral with $u^TZ_2v$ and $\overline{u^TZ_1v}$, which gives $\det(I-Z_1^\dagger Z_2)^{-1}$.
QF code (exact):
```wl
{a1, a2, b1, b2} = FieldVariables[{1, 2, 3, 4}][[{1, 3, 5, 7}]];
{ad1, ad2, bd1, bd2} = SuperDagger /@ {a1, a2, b1, b2};
z1 = {{x11, x12}, {x21, x22}}; z2 = {{y11, y12}, {y21, y22}};
cre = Total[Flatten[z2 {{ad1 ** bd1, ad1 ** bd2}, {ad2 ** bd1, ad2 ** bd2}}]];
ann = Total[Flatten[Conjugate[z1] {{a1 ** b1, a1 ** b2}, {a2 ** b1, a2 ** b2}}]];
ov[k_] := BosonicMatrixElement[{{0, 0, 0, 0}, {0, 0, 0, 0}}, NonCommutativeMultiply @@ Join[Table[ann, k], Table[cre, k]],
    "Generators" -> {a1, a2, b1, b2}]/(k!)^2;
target = Series[1/Det[IdentityMatrix[2] - s ConjugateTranspose[z1] . z2], {s, 0, 3}];
Expand[Table[ov[k] - SeriesCoefficient[target, k], {k, 1, 3}]]
```
Output: `{0, 0, 0}`. Matches: the differences at orders 1, 2, 3 are 0.

### Bosonic pair-number generating determinant (trace-log sign corrected)
Topic: Two-mode squeezing
Appears in: P2-023

The pair-number distribution of a multimode paired vacuum has generating function $\det(I-Z^\dagger Z)/\det(I-tZ^\dagger Z)$, and its logarithm is the trace series with a plus sign.

$$ \mathcal G(t)=\sum_nW_nt^n=\frac{\det(I_N-Z^\dagger Z)}{\det(I_N-tZ^\dagger Z)},\qquad \log\mathcal G=\log\det(I_N-Z^\dagger Z)+\sum_k\frac{\operatorname{Tr}[(Z^\dagger Z)^k]t^k}{k} $$

Audit: Correct as corrected; the printed minus sign before the trace sum is wrong. Checked by: hand derivation from the overlap entry: $t^{\hat N_a}\hat a_i^\dagger t^{-\hat N_a}=t\hat a_i^\dagger$ (QF exact) gives $t^{\hat N_a}|Z\rangle\propto|tZ\rangle$, so $\mathcal G(t)=\det(I-Z^\dagger Z)\langle0|\ldots\rangle_{Z,tZ}=\det(I-Z^\dagger Z)/\det(I-tZ^\dagger Z)$; exact WL computation for a symbolic $2\times2$ matrix $M=Z^\dagger Z$: the plus-sign trace series matches $\log\mathcal G$ through $t^4$, and the minus sign already fails at order $t$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
mm = {{m11, m12}, {m21, m22}}; id = IdentityMatrix[2];
logG = Log[Det[id - mm]/Det[id - t mm]];
tr[sg_] := Log[Det[id - mm]] + sg Sum[Tr[MatrixPower[mm, k]] t^k/k, {k, 1, 4}];
{Simplify[TrigToExp[BosonicExpSimplify[Exp[lam ad ** a] ** ad ** Exp[-lam ad ** a]]], lam > 0],
 Simplify[Normal[Series[logG - tr[+1], {t, 0, 4}]]], Simplify[Normal[Series[logG - tr[-1], {t, 0, 1}]]]}
```
Output: `{E^lam*SuperDagger[\[FormalA]1], 0, 2*(m11 + m22)*t}`. Matches: $e^{\lambda\hat N}\hat a^\dagger e^{-\lambda\hat N}=e^\lambda\hat a^\dagger$ (simplified for $\lambda>0$; the raw output is even in $\sqrt{\lambda^2}$, so it holds for all $\lambda$); plus sign residual 0; minus sign residual $2\operatorname{Tr}M\,t$.

### Multimode pair Riccati equation and vacuum phase
Topic: Quadratic propagators and Riccati
Appears in: P2-024

The pair matrix of a multimode paired vacuum obeys a matrix Riccati equation, and its vacuum phase obeys a trace equation.

$$ i\dot Z=A+ZA^\dagger Z+CZ+ZD^T,\qquad \dot\phi=\tfrac12\operatorname{Tr}(A^\dagger Z+Z^\dagger A) $$

Conditions: $\hat H=A_{ij}\hat a_i^\dagger\hat b_j^\dagger+\bar A_{ij}\hat a_i\hat b_j+C_{ij}\hat a_i^\dagger\hat a_j+D_{ij}\hat b_i^\dagger\hat b_j$.
Audit: Correct. Checked by: hand derivation. The pair-creation coefficient gives $[C_{ij}\hat a_i^\dagger\hat a_j,\,Z_{kl}\hat a_k^\dagger\hat b_l^\dagger]\to(CZ)$, the $D$ term gives $ZD^T$, and the double commutator of $\bar A_{ij}\hat a_i\hat b_j$ gives $ZA^\dagger Z$. The vacuum coefficient gives $\dot\phi+i\dot N/N=\operatorname{Tr}(A^\dagger Z)$; its real part is the stated phase, and its imaginary part reproduces $\frac{d}{dt}\tfrac12\log\det(I-Z^\dagger Z)$.
QF code (exact): none. QF has no time-dependent Gaussian-state equations.

### Single fermionic pair as spin 1/2
Topic: Fermionic canonical transformations
Appears in: P2-025

One fermionic particle-antiparticle pair is a spin-1/2 problem: the pair operators act as $\hat J_\pm,\hat J_0$ on the vacuum and the occupied pair.

$$ \hat J_+=\hat a^\dagger\hat b^\dagger,\ \hat J_-=\hat b\hat a,\ \hat J_0=\tfrac12(\hat a^\dagger\hat a+\hat b^\dagger\hat b-\hat I),\qquad |\zeta\rangle=(1+|\zeta|^2)^{-1/2}(\hat I+\zeta\hat a^\dagger\hat b^\dagger)|0\rangle $$

Audit: Correct. Checked by: hand derivation; on $\{|0\rangle,\hat a^\dagger\hat b^\dagger|0\rangle\}$, $[\hat J_+,\hat J_-]=\mathrm{diag}(-1,1)=2\hat J_0$; exact fermion matrices (4 modes) confirm $[\hat J_+,\hat J_-]-2\hat J_0$ annihilates the vacuum. The Dirac-field reductions are Unchecked.
QF code (exact): none. no fermionic layer.

### Multimode fermionic blocks
Topic: Fermionic canonical transformations
Appears in: P2-026

The doubled fermionic coefficient matrix is unitary, which fixes the block relations.

$$ MM^\dagger=I_{2N},\quad AA^\dagger+BB^\dagger=I_N,\quad DD^\dagger+CC^\dagger=I_N,\quad AC^\dagger=-BD^\dagger $$

Audit: Correct. Checked by: hand derivation, blockwise expansion of $MM^\dagger=I$ (same algebra as the exact fermionic Bogoliubov check).
QF code (exact): none. no fermionic layer.

### Fermionic paired vacuum normalization
Topic: Gaussian vacua and overlaps
Appears in: P2-027

The fermionic paired vacuum $e^{\sum Z_{ij}\hat a_i^\dagger\hat b_j^\dagger}|0,0\rangle$ has squared norm $\det(I_N+Z^\dagger Z)$.

$$ |Z\rangle=\det(I_N+Z^\dagger Z)^{-1/2}e^{\sum_{ij}Z_{ij}\hat a_i^\dagger \hat b_j^\dagger}|0,0\rangle $$

Audit: Correct (the chart relation to $A,B,C,D$ is Unchecked). Checked by: exact fermion matrices (4 modes), exact rational complex $2\times2$ $Z$: squared norm minus $\det(I+Z^\dagger Z)$ is exactly 0.
QF code (exact): none. no fermionic layer.

### Fermionic paired vacuum overlap
Topic: Gaussian vacua and overlaps
Appears in: P2-028

Two fermionic paired vacua overlap by $\det(I_N+Z_1^\dagger Z_2)$ over the normalizations.

$$ \langle Z_1|Z_2\rangle=\frac{\det(I_N+Z_1^\dagger Z_2)}{\det(I_N+Z_1^\dagger Z_1)^{1/2}\det(I_N+Z_2^\dagger Z_2)^{1/2}},\qquad Z'=(AZ+B)(CZ+D)^{-1} $$

Audit: Correct for the overlap; the group action is Unchecked (convention). Checked by: exact fermion matrices (4 modes), exact rational $Z_1,Z_2$: residual exactly 0.
QF code (exact): none. no fermionic layer.

### Fermionic pair-number generating determinant
Topic: Gaussian vacua and overlaps
Appears in: P2-029

The fermionic pair-number generating function is $\det(I+tZ^\dagger Z)/\det(I+Z^\dagger Z)$, a product of Bernoulli laws.

$$ \mathcal F(t)=\sum_{n=0}^NW_nt^n=\frac{\det(I_N+tZ^\dagger Z)}{\det(I_N+Z^\dagger Z)},\qquad p_j=\frac{\lambda_j}{1+\lambda_j} $$

Audit: Correct. Checked by: exact fermion matrices (4 modes), $t^{\hat N_a}=\prod_i(1+(t-1)\hat n_{a_i})$ with symbolic $t$: residual exactly 0 as a polynomial in $t$. Bernoulli factorization follows from the eigenvalues of $Z^\dagger Z$.
QF code (exact): none. no fermionic layer.

### Quasiclassical partition-function bounds
Topic: Thermal states and traces
Appears in: P2-032

For lattice spin systems the quantum partition function is squeezed between classical ones at $\beta$ and at a slightly larger $\beta$.

$$ Z_{\rm cl}(\beta)\le Z_Q^{(N)}(\beta)\le Z_{\rm cl}(\beta(1+a/N)) $$

Audit: Unchecked: Simon's theorem cited from the source; the lower bound is the Berezin-Lieb lower bound applied to a product of spin coherent states.
QF code (exact): none. QF has no partition-function bounds.

### Thermal oscillator master equation
Topic: Thermal states and traces
Appears in: P2-034

The thermal damped-oscillator master equation makes $\langle\hat a\rangle$ decay at rate $\gamma/2$ and $\langle\hat n\rangle$ relax to $\nu$ at rate $\gamma$.

$$ \dot{\hat\rho}=\gamma(\nu+1)\left(\hat a\hat\rho \hat a^\dagger-\tfrac12\{\hat a^\dagger \hat a,\hat\rho\}\right)+\gamma\nu\left(\hat a^\dagger\hat\rho \hat a-\tfrac12\{\hat a\hat a^\dagger,\hat\rho\}\right),\quad \nu=(e^{\hbar\omega/kT}-1)^{-1} $$

Audit: Correct (standard). Checked by: QF exact, the adjoint generator $\mathcal L^\dagger(X)=\gamma(\nu+1)(\hat a^\dagger X\hat a-\tfrac12\{\hat a^\dagger\hat a,X\})+\gamma\nu(\hat aX\hat a^\dagger-\tfrac12\{\hat a\hat a^\dagger,X\})$ applied to $\hat a$ and $\hat a^\dagger\hat a$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
adj[x_] := g (nu + 1) (ad ** x ** a - (ad ** a ** x + x ** ad ** a)/2) + g nu (a ** x ** ad - (a ** ad ** x + x ** a ** ad)/2);
BosonicNormalOrder[#, {a, ad}] & /@ {adj[a], adj[ad ** a]}
```
Output: `{\[FormalA]1*((g*nu)/2 - (g*(1 + nu))/2), g*nu + (g*nu - g*(1 + nu))*SuperDagger[\[FormalA]1]**\[FormalA]1}`. Matches: $\mathcal L^\dagger(\hat a)=-\tfrac\gamma2\hat a$ and $\mathcal L^\dagger(\hat a^\dagger\hat a)=-\gamma(\hat a^\dagger\hat a-\nu)$.

### Exact attenuation solution for ordered characteristic functions
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: P2-036, P2-039 (thermal $\mu(t)=p\nu_0+q\nu$)

Under thermal attenuation each ordered characteristic function is the initial one at contracted argument times a Gaussian that grows to the noise level.

$$ \chi_k(\eta,t)=e^{-q(\nu+\sigma_k)|\eta|^2}\chi_k(\sqrt p\,\eta,0),\qquad \partial_t\chi_k=-\gamma[\tfrac12\eta_j\partial_{\eta_j}+(\nu+\sigma_k)|\eta|^2]\chi_k,\quad p=e^{-\gamma t},\ q=1-p $$

Audit: Correct. Checked by: exact WL computation: with $\eta,\bar\eta$ independent and an arbitrary initial function, the stated solution satisfies the PDE identically. For thermal input $\chi_N=e^{-\nu_0|\eta|^2}$ the solution is thermal with $\mu=p\nu_0+q\nu$ (hand).
QF code (exact): none. QF has no characteristic-function evolution.

### Fokker-Planck drift, diffusion and covariance evolution
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: P2-037

The ordered quasiprobabilities obey Fokker-Planck equations with drift $-\gamma\alpha/2$, and their covariance relaxes linearly to the noise level plus the ordering offset.

$$ W_N=P,\ W_O=W_{\rm Wigner},\ W_A=Q/\pi;\qquad D(t)=pD(0)+q(\nu+\sigma_k)I_2/2 $$

Audit: Correct. Checked by: hand derivation; a Gaussian factor $e^{-s|\eta|^2}$ in $\chi_k$ corresponds to variance $s/2$ per real quadrature, giving the additive $q(\nu+\sigma_k)/2$.
QF code (exact): none. QF has no quasiprobability evolution.

### Gaussian convolution kernel for attenuation (sign corrected)
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: P2-038

Attenuation acts on each ordered quasiprobability as convolution with a normalized Gaussian centered at $\sqrt p\,\alpha'$.

$$ W_k(\alpha,t)=\int d^2\alpha'\,W_k(\alpha',0)\frac{\exp[-|\alpha-\sqrt p\,\alpha'|^2/(q(\nu+\sigma_k))]}{\pi q(\nu+\sigma_k)} $$

Audit: Correct as corrected; the printed positive exponent is not integrable. Checked by: exact WL computation: the corrected kernel integrates to 1, and its Fourier transform is $e^{-s|k|^2/4}$ with $s=q(\nu+\sigma_k)$, the Gaussian factor of the characteristic-function solution.
QF code (exact): none. QF has no quasiprobability evolution.

### Analytic density kernel
Topic: Fock matrix elements and generating functions
Appears in: P2-040

The unnormalized coherent matrix element of a density operator is a generating function whose Taylor coefficients are the Fock matrix elements.

$$ R(\bar\alpha,\beta)=e^{(|\alpha|^2+|\beta|^2)/2}\langle\alpha|\hat\rho|\beta\rangle=\sum_{mn}\rho_{mn}\frac{\bar\alpha^m\beta^n}{\sqrt{m!n!}} $$

Audit: Correct (with $z_1\leftrightarrow\bar\alpha$, $\bar z_2\leftrightarrow\beta$). Checked by: hand derivation from the coherent expansion $e^{|\alpha|^2/2}|\alpha\rangle=\sum_n\alpha^n|n\rangle/\sqrt{n!}$.
QF code (exact): none. The displaced thermal entry below exercises this kernel with QF.

### Displaced thermal density kernel
Topic: Thermal states and traces
Appears in: P2-043

The coherent kernel of a displaced thermal state is a single Gaussian in the two labels, with mean occupation $\mu$ and displacement $\beta$.

$$ R(z_1,\bar z_2;t)=\frac1{1+\mu}\exp\left[\frac{\mu z_1\bar z_2+\beta z_1+\bar\beta\bar z_2-|\beta|^2}{1+\mu}\right] $$

Audit: Correct. Checked by: QF exact: $\hat\rho=\hat D(\beta)\,e^{\lambda\hat N}\hat D(\beta)^\dagger/(1+\mu)$ with $e^\lambda=\mu/(1+\mu)$; the displacements move the coherent labels with the Weyl phases (`BosonicExpSimplify`), and the coherent element of $e^{\lambda\hat N}$ is QF's closed form. The kernel with $z_1=\bar\alpha$, $\bar z_2=\alpha'$ divided by the stated $R$ is 1.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
d[x_] := DisplacementOperator[x, Infinity] /. Conjugate[s_Plus] :> Map[Conjugate, s];
ph[x_] := Simplify[BosonicExpSimplify[d[-b] ** d[x]]/d[x - b]];
inner = BosonicMatrixElement[{al - b, alp - b}, Exp[Log[mu/(1 + mu)] ad ** a], "Basis" -> "Coherent"]/(1 + mu);
kern = Exp[(al Conjugate[al] + alp Conjugate[alp])/2] Conjugate[ph[al]] ph[alp] inner;
target = Exp[(mu Conjugate[al] alp + b Conjugate[al] + Conjugate[b] alp - b Conjugate[b])/(1 + mu)]/(1 + mu);
re = {al -> x1 + I y1, alp -> x2 + I y2, b -> u + I v};
Simplify[ComplexExpand[kern/target /. re], mu > 0]
```
Output: `1`. Matches.

### Populations of a displaced thermal state
Topic: Thermal states and traces
Appears in: P2-044

The populations of a displaced thermal state are a geometric factor times a Gaussian times a Laguerre polynomial of negative argument.

$$ W_n=\frac{\mu^n}{(1+\mu)^{n+1}}e^{-|\beta|^2/(1+\mu)}L_n\left(-\frac{|\beta|^2}{\mu(1+\mu)}\right) $$

Audit: Correct, including the Poisson limit $\mu\to0$. Checked by: exact WL computation: $n!$ times the $z_1^n\bar z_2^n$ Taylor coefficient of the kernel $R$ (QF exact in the previous entry) equals the stated $W_n$ for $n\le6$ with symbolic $\mu,\beta$.
QF code (exact): none. `BosonicMatrixElement` handles $\hat D(\beta)$ or a function of $\hat N$ alone, not their product; the populations are read off the exact kernel of the previous entry.

### Drive equation and displacement channel factor
Topic: Forced and parametric oscillators
Appears in: P2-042

A linear drive adds $-i\epsilon$ to the decay of the mean field, and the drive enters the state as a displacement.

$$ \hat\rho\mapsto \hat D(v)\hat\rho \hat D(v)^\dagger,\quad \hat D(v)=e^{v\hat a^\dagger-\bar v \hat a};\qquad \hat H_d=\epsilon\hat a^\dagger+\bar\epsilon\hat a\ \Rightarrow\ \dot\beta=-(\gamma/2)\beta-i\epsilon(t) $$

Audit: Correct. Checked by: QF exact for the drive term, $-i[\hat a,\hat H_d]=-i\epsilon$; the damping term $-\tfrac\gamma2\langle\hat a\rangle$ is the master-equation entry.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
hd = eps ad + Conjugate[eps] a;
BosonicNormalOrder[-I (a ** hd - hd ** a), {a, ad}]
```
Output: `(-I)*eps`. Matches.

### Free oscillator propagator and thermal trace in coherent states
Topic: Thermal states and traces
Appears in: P2-046

The coherent matrix element of $e^{-it\omega\hat a^\dagger\hat a}$ is the overlap with a rotated label, and the thermal trace is the geometric series.

$$ \langle\alpha|e^{-it\omega\hat a^\dagger\hat a}|\beta\rangle=\exp\left[\bar\alpha\beta e^{-i\omega t}-\tfrac12(|\alpha|^2+|\beta|^2)\right],\qquad Z=(1-e^{-\omega/T})^{-1} $$

Audit: Correct (and restoring vacuum energy multiplies by $e^{-\omega/(2T)}$). Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicMatrixElement[{al, be}, Exp[lam ad ** a], "Basis" -> "Coherent"], BosonicExpOrder[Exp[lam ad ** a]],
 Sum[Exp[-w k/T], {k, 0, Infinity}, Assumptions -> w/T > 0]}
```
Output: `{E^((-Abs[al]^2 - Abs[be]^2)/2 + be*E^lam*Conjugate[al]), NormalOrdered[E^((-1 + E^lam)*SuperDagger[\[FormalA]1]**\[FormalA]1)], E^(w/T)/(-1 + E^(w/T))}`. Matches: with $\lambda=-i\omega t$ or $-\omega/T$ the element is $e^{\bar\alpha\beta e^\lambda-(|\alpha|^2+|\beta|^2)/2}$, and the trace is $(1-e^{-\omega/T})^{-1}$.

### Stable oscillator proper-time trace
Topic: Thermal states and traces
Appears in: P2-048

The trace of $e^{-s(2\hat N+1)}$ is $1/(2\sinh s)$ for $s>0$; the source's diagonalizing $\hat U$ is a complex continuation, not a unitary.

$$ \operatorname{Tr}e^{-s(2\hat N+\hat I)}=[2\sinh s]^{-1},\ s>0;\qquad \hat U=e^{i\pi(\hat b^{\dagger2}-\hat b^2)/8}\ \text{(not unitary)} $$

Audit: Correct, with the atlas's qualification that $\hat U$ has a Hermitian generator and is a complex canonical continuation. Checked by: exact WL computation (the sum minus $1/(2\sinh s)$ simplifies to 0); hand check that $i(\hat b^{\dagger2}-\hat b^2)$ is Hermitian.
QF code (exact): none. A trace of a function of $\hat N$ is an ordinary series; QF adds nothing.

### Laguerre generating function in the lower index
Topic: Fock matrix elements and generating functions
Appears in: P2-053

The generating function of $L_n^k$ in the lower index is $(1-t)^{-k-1}e^{-tx/(1-t)}$.

$$ \sum_{n\ge0}t^nL_n^k(x)=(1-t)^{-k-1}e^{-tx/(1-t)},\quad |t|<1 $$

Audit: Correct. Checked by: exact WL computation, series coefficients through $t^6$ with symbolic $k$ and $x$ all equal.
QF code (exact): none. Not an operator identity.

### Laguerre generating function along the diagonal m-n
Topic: Fock matrix elements and generating functions
Appears in: P2-053

Along the diagonal $m-n$ the Laguerre generating function is $(1+t)^me^{-tx}$.

$$ \sum_{n\ge0}t^nL_n^{m-n}(x)=(1+t)^m e^{-tx} $$

Audit: Correct. Checked by: exact WL computation, $m=0..6$ through order $t^{10}$, all differences 0.
QF code (exact): none. Not an operator identity.

### Generating function of displacement transition probabilities (prefactor corrected)
Topic: Fock matrix elements and generating functions
Appears in: P2-054

The double generating function of the displacement transition probabilities is $(1-uv)^{-1}\exp[-|\gamma|^2(1-u)(1-v)/(1-uv)]$.

$$ F(u,v;\gamma)=\sum_{mn}u^mv^n|D_{mn}(\gamma)|^2=\frac1{1-uv}\exp\left[-\frac{|\gamma|^2(1-u)(1-v)}{1-uv}\right] $$

Audit: Correct as corrected; the printed extra $e^{|\gamma|^2}$ is wrong. Checked by: QF exact: the exact matrix elements $D_{mn}(\gamma)$ with $\gamma=g_x+ig_y$ symbolic give $\sum_{m,n\le5}u^mv^n|D_{mn}|^2$, which equals the double Taylor polynomial of the corrected form through $u^5v^5$; at $u=v=0$ the sum is $e^{-|\gamma|^2}$, while the printed form gives 1.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
dm[m_, n_] := ComplexExpand[BosonicMatrixElement[{m, n}, DisplacementOperator[gx + I gy, Infinity, "Ordering" -> "Normal"]]];
f = Sum[u^m v^n dm[m, n] ComplexExpand[Conjugate[dm[m, n]]], {m, 0, 5}, {n, 0, 5}];
rhs = Normal[Series[Exp[-(gx^2 + gy^2) (1 - u) (1 - v)/(1 - u v)]/(1 - u v), {u, 0, 5}, {v, 0, 5}]];
{Simplify[Expand[f - rhs]], f /. {u -> 0, v -> 0}}
```
Output: `{0, E^(-gx^2 - gy^2)}`. Matches the corrected form: difference 0, and $F(0,0)=e^{-|\gamma|^2}$.

### Displacement matrix elements as Wirtinger derivatives
Topic: Fock matrix elements and generating functions
Appears in: P2-055

Each displacement matrix element is a mixed Wirtinger derivative of $e^{-|\alpha|^2}$, rescaled; transitions are classically allowed where the two number circles intersect.

$$ D_{mn}(\alpha)=\frac{e^{|\alpha|^2/2}}{\sqrt{m!n!}}(-\partial_{\bar\alpha})^m(\partial_\alpha)^n e^{-|\alpha|^2},\qquad |\sqrt n-|\gamma||<\sqrt m<\sqrt n+|\gamma|\ \text{(semiclassical allowed region)} $$

Audit: Correct for the exact formula; the semiclassical region is Unchecked (asymptotic). Checked by: QF exact, symbolic with $\alpha,\bar\alpha$ independent, against the normal-ordered matrix elements for $m,n\le5$: all differences 0.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
der[m_, n_] := Exp[al alb/2]/Sqrt[m! n!] (-1)^m D[Exp[-al alb], {alb, m}, {al, n}];
nrm[m_, n_] := BosonicMatrixElement[{m, n}, DisplacementOperator[al, Infinity, "Ordering" -> "Normal"]] /. Conjugate[al] -> alb;
Union[Flatten[Table[Simplify[nrm[m, n] - der[m, n]], {m, 0, 5}, {n, 0, 5}]]]
```
Output: `{0}`. Matches.

### Jacobians of classical-domain transformations
Topic: Lie-group decompositions and generalized coherent states
Appears in: P2-056

The complex Jacobian of the fractional-linear action on a type-I domain is $\det(B^Tz+D)^{-(p+q)}$.

$$ \det(B^Tz+D)^{-(p+q)} $$

Audit: Unchecked: depends on the source's transposed-block convention for type-I domains; not re-derived.
QF code (exact): none. Not an operator identity.

### Completeness criterion for discrete SU(1,1) subsystems
Topic: Lie-group decompositions and generalized coherent states
Appears in: P2-004

A discrete SU(1,1) coherent subsystem is complete when the fundamental area is below $\pi/(2k-1)$ and incomplete for large $k$.

$$ S_\Gamma<\frac{\pi}{2k-1}\ \Rightarrow\ \text{complete};\qquad k>m_0^++\tfrac12\ \Rightarrow\ \text{incomplete} $$

Audit: Unchecked: theorem cited from the source (automorphic-form argument), area convention source-specific.
QF code (exact): none. QF has no su(1,1) layer.

## Non-formula records

- P-004: Stone-von Neumann uniqueness and its domain qualification (prose).
- P-021: group-orbit coherent states, stabilizer phase, Schur resolution (prose).
- P-022: isotropy subalgebra test for analytic realizations (prose).
- P-023: Gauss decomposition on the big cell of SL(2,C), singular at a vanishing pivot (prose; the 2x2 instance is the checked P-029 rep).
- P-027: spin-to-oscillator contraction $\hat J_\pm\sim\sqrt{2j}(\hat a^\dagger,\hat a)$, $\zeta=\alpha/\sqrt{2j}$ (asymptotic prose; consistent with P-026, since $(1+\bar\zeta\eta)^{2j}\to e^{\bar\alpha\beta}$).
- P-033: principal and Lorentz series scope boundary (prose).
- P-040: Siegel disk geometry and the $k>2N$ integral restriction (prose, Unchecked).
- P-046: nilpotent groups, truncating BCH, Kirillov orbits (prose).
- P-049: principal-series and Iwasawa scope record (prose).
- P2-001: principal-series systems are not square integrable (prose).
- P2-003: horospherical kernels (prose, source lookup).
- P2-008: general Kahler quantization (prose, source lookup).
- P2-009: linear Lie-algebra Hamiltonians preserve coherent orbits (prose).
- P2-030: Casimir multiplicities and finite/infinite block restrictions (prose).
- P2-031: Clebsch-Gordan generating polynomial (source lookup, normalization not transcribed).
- P2-033: bilinear fermion algebras and large-N mechanics (prose).
- P2-041: driven relaxation Green kernel (source lookup, not transcribed).
- P2-045: number-state relaxation and spin relaxation (prose, source lookup).
- P2-047: Landau-level and Fermi-gas applications (prose).
- P2-049: synchrotron and Wehrl-entropy applications (prose, historical).
- P2-050: coherent-lattice completeness theorem (merged with the lattice entry above).
- P2-057: addendum and references coverage record.
- P2-058: source-audit cautions (cross-references to P2-017, P2-023, P2-038, P2-054, P2-048, P-029; all confirmed above).
