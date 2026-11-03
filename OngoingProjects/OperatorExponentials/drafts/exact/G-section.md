# Section G: Gerry and Knight, Introductory Quantum Optics (atlas lines 4361 to 5648)

## Errors found

No error in the atlas formulas of this section; every entry was rechecked with exact algebra only. The exact checks confirm the atlas's corrections of three printed equations: the antinormal characteristic function of Eq. (3.128c) (a creation-first product is the normal-ordered function, which differs from the antinormal one by e^{|lambda|^2}; G-028), the odd-number phase of Eq. (7.142) (it is -i, not +i; G-077), and the missing number shift in a e^{-2iKt n} (off by e^{2iKt}; G-079). One remark stays unchecked: the interferometer fringe sign in G-050 depends on the book's Eq. (11.2) layout, which the atlas does not reproduce.

Every `wl` cell below runs as written after these loading lines, with `{a, ad} = FieldVariables[{1}]` (or `FieldVariables[{1, 2}]`, `FieldVariables[{1, 2, 3}]` for several modes) at its top:

```wl
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

The cells are the files `drafts/exact/G-work/cells/NNN.wl` (NNN = entry position), run by `drafts/exact/G-work/driver.wl` through `runner.wls`; their raw output is `drafts/exact/G-work/out-cells.txt`, and every cell ran without a message. Exact WL checks for entries without a QF cell are in `drafts/exact/G-work/wl/NNN.wl`, output `out-wl.txt`. No truncated Fock matrix is used anywhere. Symbols such as `al`, `be`, `xi`, `r`, `th` are left symbolic; `r > 0` and real `th` are passed as assumptions where a cell simplifies.

### Canonical oscillator realization of the CCR
Topic: Coherent states (general)
Appears in: G-001

Writing the annihilation operator as a complex combination of position and momentum turns the canonical commutator into the unit bosonic commutator, and the oscillator Hamiltonian becomes the number operator plus a zero-point half.

$$ [\hat q,\hat p]=i\hbar \hat I,\qquad \hat a=\frac{\omega \hat q+i\hat p}{\sqrt{2\hbar\omega}},\qquad [\hat a,\hat a^\dagger]=\hat I,\qquad \hat H=\hbar\omega(\hat n+\tfrac12\hat I). $$

Audit: Correct. Checked by: QF exact (with q and p rebuilt from a and a^dagger, hbar and omega symbolic).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
q = Sqrt[hb/(2 w)] (a + ad); p = I Sqrt[hb w/2] (ad - a);
PowerExpand[{BosonicNormalOrder[a ** ad - ad ** a, {a, ad}], BosonicNormalOrder[q ** p - p ** q, {a, ad}], BosonicNormalOrder[(p ** p + w^2 q ** q)/2, {a, ad}]}]
```
Output: `{1, I*hb, (hb*w)/2 + hb*w*SuperDagger[\[FormalA]1]**\[FormalA]1}`. Matches: [a, a^dagger] = 1, [q, p] = i hbar, and H = hbar omega (n + 1/2).

### Hadamard conjugation series
Topic: BCH and Hadamard
Appears in: G-002 (applied in G-018, G-044, G-052, G-068)

Conjugating B by an exponential of A sums the nested commutators of A with B; the series closes in finitely many terms when the nested commutators vanish or repeat.

$$ e^{s\hat A}\hat Be^{-s\hat A}=\sum_{k\ge0}\frac{s^k}{k!}\widehat{\operatorname{ad}}_{\hat A}^k(\hat B),\quad\widehat{\operatorname{ad}}_{\hat A}(\hat B)=[\hat A,\hat B]. $$

Audit: Correct (standard). Checked by: hand derivation (both sides solve dX/ds = [A, X] with X(0) = B); QF exact on two instances.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
FullSimplify[{BosonicExpSimplify[Exp[s ad ** a] ** ad ** Exp[-s ad ** a]], BosonicExpSimplify[Exp[s a ** a] ** ad ** Exp[-s a ** a]]}, s \[Element] Reals]
```
Output: `{E^s*SuperDagger[\[FormalA]1], 2*\[FormalA]1*s + SuperDagger[\[FormalA]1]}`. Matches the series: e^{s n} a^dagger e^{-s n} = e^s a^dagger (rotation closure) and e^{s a^2} a^dagger e^{-s a^2} = a^dagger + 2 s a (series stops after one commutator).

### Central-commutator BCH
Topic: BCH and Hadamard
Appears in: G-003

When the commutator of A and B commutes with both, the exponential of the sum splits into the product of exponentials times a scalar exponential of half the commutator, with the sign set by the order.

$$ e^{\hat A+\hat B}=e^{\hat A}e^{\hat B}e^{-[\hat A,\hat B]/2}=e^{\hat B}e^{\hat A}e^{[\hat A,\hat B]/2} $$

Conditions: $[\hat A,[\hat A,\hat B]]=[\hat B,[\hat A,\hat B]]=0$.
Audit: Correct. Checked by: QF exact on A = x a^dagger, B = y a (so [A, B] = -x y, central).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicExpOrder[Exp[x ad + y a]], BosonicExpOrder[Exp[x ad + y a], "Ordering" -> "Antinormal"]}
```
Output: `{E^((x*y)/2)*E^(x*SuperDagger[\[FormalA]1])**E^(\[FormalA]1*y), E^(\[FormalA]1*y)**E^(x*SuperDagger[\[FormalA]1])/E^((x*y)/2)}`. Matches: normal order gives e^{xy/2} e^{x a^dagger} e^{y a} = e^A e^B e^{-[A,B]/2}, antinormal order gives e^{-xy/2} e^{y a} e^{x a^dagger} = e^B e^A e^{[A,B]/2}.

### Number commutators and free rotation of the field
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-004 (sign convention reused in G-101)

The number operator lowers and raises with a and a^dagger, so a free rotation multiplies a by a phase and carries a coherent state around a circle, the zero-point term adding only a global phase.

$$ [\hat n,\hat a]=-\hat a,\quad[\hat n,\hat a^\dagger]=\hat a^\dagger,\quad e^{i\omega t \hat n}\hat ae^{-i\omega t \hat n}=e^{-i\omega t}\hat a,\quad e^{-i\omega t(\hat n+\tfrac12\hat I)}|\alpha\rangle=e^{-i\omega t/2}|\alpha e^{-i\omega t}\rangle. $$

Audit: Correct. Checked by: QF exact (commutators, conjugation, and coherent-basis matrix elements; the coherent-state identity holds for every bra beta, so it holds as a vector identity).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicNormalOrder[ad ** a ** a - a ** ad ** a, {a, ad}], BosonicNormalOrder[ad ** a ** ad - ad ** ad ** a, {a, ad}],
 FullSimplify[BosonicExpSimplify[Exp[I wt ad ** a] ** a ** Exp[-I wt ad ** a]], wt > 0],
 FullSimplify[Exp[-I wt/2] BosonicMatrixElement[{be, al}, Exp[-I wt ad ** a], "Basis" -> "Coherent"] == Exp[-I wt/2] BosonicMatrixElement[{be, al Exp[-I wt]}, 1, "Basis" -> "Coherent"], wt > 0]}
```
Output: `{-\[FormalA]1, SuperDagger[\[FormalA]1], \[FormalA]1/E^(I*wt), True}`. Matches: [n, a] = -a, [n, a^dagger] = a^dagger, e^{i wt n} a e^{-i wt n} = e^{-i wt} a, and <beta|e^{-i wt (n + 1/2)}|alpha> = e^{-i wt/2} <beta|alpha e^{-i wt}> for all beta.

### Fock ladder action and resolution of identity
Topic: Fock matrix elements and generating functions
Appears in: G-005

The ladder operators move one step on the number basis with square-root weights, the n-photon state is the normalized n-th power of a^dagger on the vacuum, and the number states resolve the identity.

$$ \hat a|n\rangle=\sqrt n|n-1\rangle,\quad \hat a^\dagger|n\rangle=\sqrt{n+1}|n+1\rangle,\quad |n\rangle=\frac{(\hat a^\dagger)^n}{\sqrt{n!}}|0\rangle,\quad\sum_{n\ge0}|n\rangle\langle n|=\hat I. $$

Audit: Correct. Checked by: QF exact (symbolic Fock indices; norm of (a^dagger)^3|0>/sqrt(3!)); hand derivation for the completeness sum.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicMatrixElement[{m, n}, a], BosonicMatrixElement[{m, n}, ad], BosonicVEV[a ** a ** a ** ad ** ad ** ad, {a, ad}]/3!}
```
Output: `{Sqrt[n]*KroneckerDelta[m, -1 + n], Sqrt[m]*KroneckerDelta[-1 + m, n], 1}`. Matches.

### Number-shift rules for functions of the number operator
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-006 (used in G-079, G-080)

Moving a through a function of the number operator shifts the argument of the function by one.

$$ \hat a f(\hat n)=f(\hat n+\hat I)\hat a,\quad f(\hat n)\hat a=\hat a f(\hat n-\hat I) $$

Audit: Correct. Checked by: QF exact for f(n) = n^3, n^2 and e^{s n}; hand derivation for general f (both sides have the Fock element f(n) sqrt(n) between <n-1| and |n>).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}]; nn = ad ** a;
{BosonicNormalOrder[a ** nn ** nn ** nn - (nn + 1) ** (nn + 1) ** (nn + 1) ** a, {a, ad}],
 BosonicNormalOrder[nn ** nn ** a - a ** (nn - 1) ** (nn - 1), {a, ad}],
 FullSimplify[BosonicExpSimplify[Exp[-s nn] ** a ** Exp[s nn]], s > 0]}
```
Output: `{0, 0, \[FormalA]1*E^s}`. Matches: both polynomial differences normal-order to 0, and e^{-s n} a e^{s n} = e^s a is a f(n) = f(n+1) a with f = e^{s n}.

### Normal-ordered powers as falling factorials
Topic: Normal symbols and normal-order exponentials
Appears in: G-006, G-078 ($\hat a^{\dagger2}\hat a^2=\hat n(\hat n-\hat I)$), G-097 (eigenvalue of ${:\hat n^r:}$)

Normal-ordered powers a^dagger^r a^r equal the falling factorial of the number operator.

$$ (\hat a^\dagger)^r \hat a^r=\hat n(\hat n-\hat I)\cdots(\hat n-r\hat I+\hat I) $$

Audit: Correct. Checked by: QF exact (operator identity for r = 3, symbolic Fock element for r = 4); hand derivation for general r (the Fock element is n!/(n-r)!).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}]; nn = ad ** a;
{BosonicNormalOrder[nn ** (nn - 1) ** (nn - 2) - ad ** ad ** ad ** a ** a ** a, {a, ad}],
 Factor[BosonicMatrixElement[{n, n}, ad ** ad ** ad ** ad ** a ** a ** a ** a]]}
```
Output: `{0, FactorialPower[n, 4]}`. Matches: the r = 3 difference is 0 and <n|a^dagger^4 a^4|n> = n(n-1)(n-2)(n-3).

### Quadrature commutator and number-state variance
Topic: Coherent states (general)
Appears in: G-007

The two quadratures have commutator i/2, so their variances obey the 1/16 product bound, and a number state has variance (2n+1)/4 in both.

$$ [\hat X_1,\hat X_2]=\frac i2\hat I,\quad \operatorname{Var}\hat X_1\operatorname{Var}\hat X_2\ge\frac1{16},\quad\operatorname{Var}_{|n\rangle}\hat X_j=\frac{2n+1}{4}. $$

Audit: Correct. Checked by: QF exact (commutator and symbolic-n variances); the bound is Robertson's inequality with this commutator.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
X1 = (a + ad)/2; X2 = (a - ad)/(2 I);
{BosonicNormalOrder[X1 ** X2 - X2 ** X1, {a, ad}],
 Simplify[BosonicMatrixElement[{n, n}, X1 ** X1] - BosonicMatrixElement[{n, n}, X1]^2],
 Simplify[BosonicMatrixElement[{n, n}, X2 ** X2] - BosonicMatrixElement[{n, n}, X2]^2]}
```
Output: `{I/2, (1 + 2*n)/4, (1 + 2*n)/4}`. Matches.

### Independent-mode CCR and multimode Fock states
Topic: Fock matrix elements and generating functions
Appears in: G-008

Different modes commute, each mode has its own unit commutator, and multimode number states are products of normalized creation powers on the joint vacuum.

$$ [\hat a_j,\hat a_k^\dagger]=\delta_{jk}\hat I,\quad[\hat a_j,\hat a_k]=[\hat a_j^\dagger,\hat a_k^\dagger]=0,\quad |\{n_j\}\rangle=\prod_j\frac{(\hat a_j^\dagger)^{n_j}}{\sqrt{n_j!}}|\{0\}\rangle. $$

Audit: Correct. Checked by: QF exact (two modes, and the norm of a1^dagger^2 a2^dagger |0,0>/sqrt(2! 1!)).
QF code (exact):
```wl
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
{BosonicNormalOrder[a1 ** ad2 - ad2 ** a1, {a1, ad1, a2, ad2}], BosonicNormalOrder[a1 ** a2 - a2 ** a1, {a1, ad1, a2, ad2}],
 BosonicNormalOrder[a1 ** ad1 - ad1 ** a1, {a1, ad1, a2, ad2}], BosonicMatrixElement[{{2, 1}, {0, 0}}, ad1 ** ad1 ** ad2/Sqrt[2! 1!]]}
```
Output: `{0, 0, 1, 1}`. Matches.

### Thermal state as a number exponential
Topic: Thermal states and traces
Appears in: G-009

A thermal state is a normalized exponential of the number operator; its eigenvalues form a geometric distribution in the mean photon number, and the partition function carries the zero-point factor.

$$ Z=\frac{e^{-x/2}}{1-e^{-x}},\quad \hat\rho_{\rm th}=(1-e^{-x})e^{-x\hat n}=\frac1{1+\bar n}\left(\frac{\bar n}{1+\bar n}\right)^{\hat n},\quad\bar n=(e^x-1)^{-1}. $$

Conditions: $x=\hbar\omega/(k_BT)>0$.
Audit: Correct. Checked by: QF exact (Fock diagonal of e^{-x n}) and exact WL sums.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
pm = (1 - E^-x) BosonicMatrixElement[{m, m}, Exp[-x ad ** a]];
{pm, FullSimplify[pm == nb^m/(1 + nb)^(m + 1) /. nb -> 1/(E^x - 1), x > 0 && m >= 0],
 Simplify[Sum[pm, {m, 0, Infinity}, Assumptions -> x > 0]], Simplify[Sum[Exp[-x (m + 1/2)], {m, 0, Infinity}, Assumptions -> x > 0] - E^(-x/2)/(1 - E^-x)]}
```
Output: `{(1 - E^(-x))/E^(m*x), True, 1, 0}`. Matches: the eigenvalues equal nbar^m/(1+nbar)^{m+1} with nbar = 1/(e^x - 1), they sum to 1, and the partition sum equals e^{-x/2}/(1 - e^{-x}).

### Thermal factorial moments and number variance
Topic: Thermal states and traces
Appears in: G-010

In a thermal state the r-th factorial moment is r! nbar^r, which gives the Bose-Einstein number variance nbar(1 + nbar).

$$ \langle(\hat a^\dagger)^r \hat a^r\rangle_{\rm th}=r!\bar n^r,\quad\langle \hat n^2\rangle=\bar n+2\bar n^2,\quad\operatorname{Var}\hat n=\bar n(1+\bar n). $$

Audit: Correct. Checked by: QF exact (symbolic Fock elements of a^dagger^r a^r) with exact WL sums over the thermal distribution.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}]; nn = ad ** a;
pm = nb^m/(1 + nb)^(m + 1);
mom[op_] := FullSimplify[Sum[pm BosonicMatrixElement[{m, m}, op], {m, 0, Infinity}, Assumptions -> nb > 0], nb > 0];
fac[r_] := NonCommutativeMultiply @@ Join[ConstantArray[ad, r], ConstantArray[a, r]];
{Table[mom[fac[r]] - r! nb^r, {r, 1, 4}], mom[nn ** nn], FullSimplify[mom[nn ** nn] - mom[nn]^2]}
```
Output: `{{0, 0, 0, 0}, nb*(1 + 2*nb), nb*(1 + nb)}`. Matches for r = 1..4, <n^2> = nbar(1 + 2 nbar), and Var n = nbar(1 + nbar).

### No Hermitian phase operator canonically conjugate to number
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-011

A Hermitian phase operator canonically conjugate to the number operator cannot exist, because its diagonal number-basis elements would have to satisfy 0 = i.

$$ (m-n)\langle m|\hat\phi|n\rangle=i\delta_{mn} $$

Audit: Correct. Checked by: hand derivation (the m = n element gives 0 = i).
QF code (exact): none. an impossibility statement; there is nothing to compute.

### Susskind-Glogower one-sided shift
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-012

The Susskind-Glogower operator shifts the number down by one with unit weight; it is an isometry from one side only, missing exactly the vacuum projector on the other.

$$ \hat E=(\hat n+\hat I)^{-1/2}\hat a=\sum_{m\ge0}|m\rangle\langle m+1|,\quad \hat E\hat E^\dagger=\hat I,\quad \hat E^\dagger \hat E=\hat I-|0\rangle\langle0|. $$

Audit: Correct. Checked by: hand derivation (E^dagger E = sum_m |m+1><m+1|) with the QF Fock element of a.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicMatrixElement[{m, n}, a], BosonicMatrixElement[{m, m + 1}, a]/Sqrt[m + 1]}
```
Output: `{Sqrt[n]*KroneckerDelta[m, -1 + n], 1}`. Matches: the only nonzero element of (n+1)^{-1/2} a is <m|.|m+1> = 1, which is the one-sided shift.

### Coherent state as eigenvector and its Fock series
Topic: Coherent states (general)
Appears in: G-013, G-016 (vacuum displacement), G-021 (Poisson coefficients)

A coherent state is the right eigenvector of a, obtained by displacing the vacuum, with Poisson-weighted Fock amplitudes.

$$ \hat a|\alpha\rangle=\alpha|\alpha\rangle,\quad |\alpha\rangle=e^{-|\alpha|^2/2}\sum_{m\ge0}\frac{\alpha^m}{\sqrt{m!}}|m\rangle,\quad\langle\alpha|\hat a^\dagger=\alpha^*\langle\alpha|. $$

Audit: Correct. Checked by: QF exact (D^dagger a D = a + alpha gives a D|0> = alpha D|0>; Fock amplitudes for m = 0..8); hand derivation for general m from the eigenvalue relation.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicExpSimplify[Exp[-(al ad - Conjugate[al] a)] ** a ** Exp[al ad - Conjugate[al] a]],
 Table[Simplify[BosonicMatrixElement[{m, 0}, DisplacementOperator[al]] - Exp[-Abs[al]^2/2] al^m/Sqrt[m!]], {m, 0, 8}]}
```
Output: `{\[FormalA]1 + al, {0, 0, 0, 0, 0, 0, 0, 0, 0}}`. Matches.

### No normalizable right eigenvector of the creation operator
Topic: Coherent states (general)
Appears in: G-014

The creation operator has no normalizable right eigenvector: the eigenvalue equation forces every Fock component to vanish.

$$ \hat a^\dagger|\psi\rangle=\lambda|\psi\rangle\ \Rightarrow\ |\psi\rangle=0 $$

Audit: Correct. Checked by: hand derivation: the components read 0 = lambda c_0 and sqrt(m) c_{m-1} = lambda c_m; for lambda != 0 this gives c_0 = 0 and then every c_m = 0, for lambda = 0 it gives every c_{m-1} = 0.
QF code (exact): none. a nonexistence statement; there is nothing to compute.

### Displacement: normal and antinormal factorizations, adjoint
Topic: Displacement
Appears in: G-015 (also G-101 at $\alpha=-i\eta$)

The displacement exponential factors into normal or antinormal order with opposite scalar factors, and its adjoint is the displacement by minus the amplitude.

$$ \hat D(\alpha)=e^{\alpha \hat a^\dagger-\alpha^*\hat a}=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}e^{-\alpha^*\hat a}=e^{|\alpha|^2/2}e^{-\alpha^*\hat a}e^{\alpha \hat a^\dagger},\quad \hat D^\dagger(\alpha)=\hat D(-\alpha). $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{DisplacementOperator[al, Infinity, "Ordering" -> "Normal"], DisplacementOperator[al, Infinity, "Ordering" -> "Antinormal"],
 BosonicExpSimplify[DisplacementOperator[al, Infinity] ** DisplacementOperator[-al, Infinity]]}
```
Output: `{E^(al*SuperDagger[\[FormalA]1])**E^(-(\[FormalA]1*Conjugate[al]))/E^((al*Conjugate[al])/2), E^((al*Conjugate[al])/2)*E^(-(\[FormalA]1*Conjugate[al]))**E^(al*SuperDagger[\[FormalA]1]), 1}`. Matches both ordered forms; D(alpha) D(-alpha) = 1, and since the generator is anti-Hermitian D^dagger(alpha) = D(-alpha).

### Vacuum action of the ordered displacement
Topic: Displacement
Appears in: G-016

In the normal-ordered displacement the annihilation factor acts trivially on the vacuum, so the coherent state is a normalized exponential of a^dagger on the vacuum.

$$ e^{-\alpha^*\hat a}|0\rangle=|0\rangle,\quad \hat D(\alpha)|0\rangle=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}|0\rangle=|\alpha\rangle,\quad \|e^{\alpha\hat a^\dagger}|0\rangle\|=e^{|\alpha|^2/2}. $$

Audit: Correct. Checked by: QF exact (Fock amplitudes m = 0..8; antinormal-to-normal reordering and vacuum element).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{Table[Simplify[BosonicMatrixElement[{m, 0}, DisplacementOperator[al]] - Exp[-Abs[al]^2/2] BosonicMatrixElement[{m, 0}, Exp[al ad]]], {m, 0, 8}],
 BosonicExpOrder[Exp[Conjugate[al] a] ** Exp[al ad]], BosonicMatrixElement[{0, 0}, BosonicExpOrder[Exp[Conjugate[al] a] ** Exp[al ad]]]}
```
Output: `{{0, 0, 0, 0, 0, 0, 0, 0, 0}, E^(al*Conjugate[al])*E^(al*SuperDagger[\[FormalA]1])**E^(\[FormalA]1*Conjugate[al]), E^(al*Conjugate[al])}`. Matches: <m|D|0> = e^{-|alpha|^2/2} <m|e^{alpha a^dagger}|0>, and <0|e^{alpha^* a} e^{alpha a^dagger}|0> = e^{|alpha|^2}, so the norm is e^{|alpha|^2/2}.

### Weyl composition law and coherent displacement phase
Topic: Displacement
Appears in: G-017, G-022 (phase of composed drives)

Two displacements compose to the displacement by the summed amplitude times a phase, so displacements commute up to the phase e^{alpha beta^* - alpha^* beta}.

$$ \hat D(\alpha)\hat D(\beta)=e^{(\alpha\beta^*-\alpha^*\beta)/2}\hat D(\alpha+\beta),\quad \hat D(\alpha)|\beta\rangle=e^{i\operatorname{Im}(\alpha\beta^*)}|\alpha+\beta\rangle,\quad \hat D(\alpha)\hat D(\beta)=e^{\alpha\beta^*-\alpha^*\beta}\hat D(\beta)\hat D(\alpha). $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
expo[e_] := Expand[e /. Power[E, x_] :> x];
dd = BosonicExpSimplify[DisplacementOperator[al, Infinity] ** DisplacementOperator[be, Infinity]];
ddRev = BosonicExpSimplify[DisplacementOperator[be, Infinity] ** DisplacementOperator[al, Infinity]];
{Simplify[expo[dd] - ((al Conjugate[be] - Conjugate[al] be)/2 + (al + be) ad - (Conjugate[al] + Conjugate[be]) a)],
 Simplify[expo[dd] - expo[ddRev] - (al Conjugate[be] - Conjugate[al] be)],
 Simplify[ComplexExpand[(al Conjugate[be] - Conjugate[al] be)/2 - I Im[al Conjugate[be]] /. {al -> x1 + I y1, be -> x2 + I y2}]]}
```
Output: `{0, 0, 0}`. Matches all three: the merged exponent, the commutation phase, and the phase (alpha beta^* - alpha^* beta)/2 = i Im(alpha beta^*).

### Displacement as a translation of the field
Topic: Displacement
Appears in: G-018 (used in G-059, G-061, G-093)

Conjugation by a displacement shifts a and a^dagger by the amplitude and its conjugate; the inverse conjugation shifts the other way.

$$ \hat D^\dagger(\alpha)\hat a\hat D(\alpha)=\hat a+\alpha\hat I,\quad \hat D^\dagger(\alpha)\hat a^\dagger \hat D(\alpha)=\hat a^\dagger+\alpha^*\hat I,\quad \hat D\hat a\hat D^\dagger=\hat a-\alpha\hat I. $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
dg = al ad - Conjugate[al] a;
{BosonicExpSimplify[Exp[-dg] ** a ** Exp[dg]], BosonicExpSimplify[Exp[-dg] ** ad ** Exp[dg]], BosonicExpSimplify[Exp[dg] ** a ** Exp[-dg]]}
```
Output: `{\[FormalA]1 + al, Conjugate[al] + SuperDagger[\[FormalA]1], \[FormalA]1 - al}`. Matches.

### Coherent overlap and resolution of identity
Topic: Coherent states (general)
Appears in: G-019, G-104 (trace formula)

Coherent states overlap with a Gaussian in their separation and resolve the identity with measure d^2 alpha/pi.

$$ \langle\beta|\alpha\rangle=e^{-(|\alpha|^2+|\beta|^2)/2+\beta^*\alpha},\quad |\langle\beta|\alpha\rangle|^2=e^{-|\alpha-\beta|^2},\quad\int_{\mathbb C}|\alpha\rangle\langle\alpha|\frac{d^2\alpha}{\pi}=\hat I. $$

Audit: Correct. Checked by: QF exact for the overlap; exact WL integrals for the resolution (the angular integral vanishes for m != n, the radial one equals 1).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
ov = BosonicMatrixElement[{be, al}, 1, "Basis" -> "Coherent"];
{ov, Simplify[ComplexExpand[ov Conjugate[ov] - Exp[-Abs[al - be]^2] /. {al -> x1 + I y1, be -> x2 + I y2}]],
 Integrate[2 rho^(2 m + 1) Exp[-rho^2]/m!, {rho, 0, Infinity}, Assumptions -> m >= 0],
 Simplify[Integrate[Exp[I k phi], {phi, 0, 2 Pi}], k \[Element] Integers && k != 0]}
```
Output: `{E^((-Abs[al]^2 - Abs[be]^2)/2 + al*Conjugate[be]), 0, 1, 0}`. Matches.

### Normally ordered coherent matrix elements
Topic: Normal symbols and normal-order exponentials
Appears in: G-020

Between coherent states a normally ordered operator is replaced by its symbol with a^dagger -> beta^* and a -> alpha, times the overlap.

$$ \langle\beta|(\hat a^\dagger)^p \hat a^q|\alpha\rangle=(\beta^*)^p\alpha^q\langle\beta|\alpha\rangle,\quad\langle\beta|{:F(\hat a^\dagger,\hat a):}|\alpha\rangle=F(\beta^*,\alpha)\langle\beta|\alpha\rangle. $$

Audit: Correct. Checked by: QF exact for p = 2, q = 1; hand derivation for general normal-ordered F (a acts on the ket, a^dagger on the bra).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
Simplify[BosonicMatrixElement[{be, al}, ad ** ad ** a, "Basis" -> "Coherent"]/BosonicMatrixElement[{be, al}, 1, "Basis" -> "Coherent"]]
```
Output: `al*Conjugate[be]^2`. Matches.

### Coherent-state photon statistics
Topic: Coherent states (general)
Appears in: G-021

A coherent state has Poisson photon statistics, factorial moments |alpha|^{2r}, and vacuum quadrature noise.

$$ P_m=e^{-|\alpha|^2}\frac{|\alpha|^{2m}}{m!},\quad \langle \hat n\rangle=\operatorname{Var}\hat n=|\alpha|^2,\quad \langle(\hat a^\dagger)^r \hat a^r\rangle=|\alpha|^{2r},\quad\operatorname{Var}\hat X_j=\frac14. $$

Audit: Correct. Checked by: QF exact (coherent-basis expectation values and Fock amplitudes m = 0..5).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}]; nn = ad ** a; X1 = (a + ad)/2;
ev[op_] := BosonicMatrixElement[{al, al}, op, "Basis" -> "Coherent"];
{FullSimplify[{ev[nn], ev[nn ** nn] - ev[nn]^2, ev[ad ** ad ** ad ** a ** a ** a], ev[X1 ** X1] - ev[X1]^2}],
 Table[Simplify[Abs[BosonicMatrixElement[{m, 0}, DisplacementOperator[al]]]^2 - Exp[-Abs[al]^2] Abs[al]^(2 m)/m!], {m, 0, 5}]}
```
Output: `{{Abs[al]^2, Abs[al]^2, Abs[al]^6, 1/4}, {0, 0, 0, 0, 0, 0}}`. Matches.

### Linearly driven evolution equals a displacement up to phase
Topic: Forced and parametric oscillators
Appears in: G-022

A linearly driven oscillator, built as a product of short displacements, is one displacement by the accumulated amplitude times an accumulated phase.

$$ \hat U(T,0)=e^{i\Phi(T)}\hat D(\alpha(T)),\qquad\alpha(T)=\int_0^T u(t)\,dt. $$

Audit: Correct (structure). Checked by: QF exact on a product of three displacements, and hand derivation from the Weyl law; the time-dependent drive u is book-specific and not rederived.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
BosonicExpSimplify[DisplacementOperator[x1, Infinity] ** DisplacementOperator[x2, Infinity] ** DisplacementOperator[x3, Infinity]]
```
Output: `E^(-(\[FormalA]1*Conjugate[x1]) - \[FormalA]1*Conjugate[x2] + (-(x2*Conjugate[x1]) + x1*Conjugate[x2])/2 - \[FormalA]1*Conjugate[x3] + (-(x3*(Conjugate[x1] + Conjugate[x2])) + (x1 + x2)*Conjugate[x3])/2 + x1*SuperDagger[\[FormalA]1] + x2*SuperDagger[\[FormalA]1] + x3*SuperDagger[\[FormalA]1])`. Matches: the product is a single displacement by x1 + x2 + x3 times a scalar phase built from the pairwise Weyl phases.

### Bargmann (entire-function) representation
Topic: Coherent states (general)
Appears in: G-023

Projecting a state on a coherent bra gives a Gaussian times an entire function of alpha^*, on which a^dagger acts as multiplication by z and a as d/dz.

$$ \langle\alpha|\psi\rangle=e^{-|\alpha|^2/2}\psi(\alpha^*),\quad \psi(z)=\sum_{m\ge0}\psi_m\frac{z^m}{\sqrt{m!}},\quad \hat a^\dagger\mapsto z,\ \hat a\mapsto\partial_z. $$

Audit: Correct. Checked by: QF exact (<alpha|m> for m = 0..6) and exact WL derivatives of z^m/sqrt(m!) matching the QF ladder elements for symbolic m.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{Table[Simplify[Conjugate[BosonicMatrixElement[{m, 0}, DisplacementOperator[al]]] - Exp[-Abs[al]^2/2] Conjugate[al]^m/Sqrt[m!]], {m, 0, 6}],
 FullSimplify[{D[z^m/Sqrt[m!], z] - Sqrt[m] z^(m - 1)/Sqrt[(m - 1)!], z z^m/Sqrt[m!] - Sqrt[m + 1] z^(m + 1)/Sqrt[(m + 1)!]}, m >= 1],
 BosonicMatrixElement[{m - 1, m}, a], BosonicMatrixElement[{m + 1, m}, ad]}
```
Output: `{{0, 0, 0, 0, 0, 0, 0}, {0, 0}, Sqrt[m], Sqrt[1 + m]}`. Matches.

### Fock matrix elements from the coherent diagonal symbol
Topic: Fock matrix elements and generating functions
Appears in: G-024

The Fock matrix elements of an operator are Taylor coefficients of its coherent-state diagonal times e^{|alpha|^2}, with alpha and alpha^* as independent variables.

$$ \langle m|\hat F|n\rangle=\frac1{\sqrt{m!n!}}\left.\partial_{\alpha^*}^{m}\partial_\alpha^n\left(e^{|\alpha|^2}\langle\alpha|\hat F|\alpha\rangle\right)\right|_{\alpha=\alpha^*=0}. $$

Audit: Correct where the series converges. Checked by: QF exact (coherent diagonal and Fock elements for a polynomial and for e^{c n}, all m, n <= 4).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
recon[f_] := Module[{g = Exp[al Conjugate[al]] BosonicMatrixElement[{al, al}, f, "Basis" -> "Coherent"] /. Abs[al] -> Sqrt[al Conjugate[al]] /. Conjugate[al] -> u /. al -> v},
  Table[Simplify[(D[g, {u, m}, {v, n}] /. {u -> 0, v -> 0})/Sqrt[m! n!] - BosonicMatrixElement[{m, n}, f]], {m, 0, 4}, {n, 0, 4}]];
{recon[ad ** a ** a + 3 a ** ad ** ad ** ad], recon[Exp[c ad ** a]]}
```
Output: `{{{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}, {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}}`. Matches: every coefficient difference is 0.

### Coherent-projector differential correspondences
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: G-025

Multiplying a coherent projector by a^dagger on the left or a on the right acts on its matrix elements as a first-order differential operator in alpha and alpha^*.

$$ \hat a^\dagger|\alpha\rangle\langle\alpha|=(\alpha^*+\partial_\alpha)|\alpha\rangle\langle\alpha|,\quad |\alpha\rangle\langle\alpha|\hat a=(\alpha+\partial_{\alpha^*})|\alpha\rangle\langle\alpha|. $$

Conditions: $\alpha,\alpha^*$ independent Wirtinger variables.
Audit: Correct. Checked by: QF exact coherent overlaps with alpha and alpha^* made independent, then exact WL differentiation; the identity holds between all coherent bras and kets.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
ov[x_, y_] := BosonicMatrixElement[{x, y}, 1, "Basis" -> "Coherent"];
f = ov[be, al] ov[al, ga] /. Abs[al] -> Sqrt[al Conjugate[al]] /. Conjugate[al] -> u /. al -> v;
{Simplify[u f + D[f, v] - Conjugate[be] f], Simplify[v f + D[f, u] - ga f]}
```
Output: `{0, 0}`. Matches.

### P representation, optical equivalence, and Q function
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: G-026

A density operator is a P-weighted mixture of coherent projectors, normal-ordered expectations are P-averages of the symbol, Q is the coherent diagonal over pi, and normal ordering of n^2 drops a term n.

$$ \hat\rho=\int P(\alpha)|\alpha\rangle\langle\alpha|d^2\alpha,\quad\langle{:G(\hat a^\dagger,\hat a):}\rangle=\int P(\alpha)G(\alpha^*,\alpha)d^2\alpha,\quad Q(\alpha)=\pi^{-1}\langle\alpha|\hat\rho|\alpha\rangle,\quad {:\hat n^2:}=\hat a^{\dagger2}\hat a^2\ne\hat n^2. $$

Audit: Correct. Checked by: QF exact for the ordering remark; hand derivation for optical equivalence (the G-020 rule inside the P integral) and the Q definition.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}]; nn = ad ** a;
{BosonicNormalOrder[nn ** nn, {a, ad}], BosonicNormalOrder[nn ** nn - ad ** ad ** a ** a, {a, ad}]}
```
Output: `{GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]**GeneralizedPower[NonCommutativeMultiply, \[FormalA]1, 2] + SuperDagger[\[FormalA]1]**\[FormalA]1, SuperDagger[\[FormalA]1]**\[FormalA]1}`. Matches: n^2 = a^dagger^2 a^2 + n, so :n^2: differs from n^2 by n.

### Glauber P function of a number state
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: G-027

The P function of a coherent state is a delta at its amplitude; that of a number state is a derivative of the delta, defined only under integration against test functions.

$$ P_{|\beta\rangle}(\alpha)=\delta^{(2)}(\alpha-\beta),\qquad P_m(\alpha)=\frac{e^{|\alpha|^2}}{m!}\partial_\alpha^m\partial_{\alpha^*}^m\delta^{(2)}(\alpha). $$

Audit: Correct as distributions. Checked by: exact WL computation with QF coherent overlaps: integrating P_m against <beta|alpha><alpha|gamma> by parts gives (1/m!) d^m d^m [e^{alpha alpha^*} f] at 0, compared with <beta|m><m|gamma> for m = 0..6; hand derivation for the coherent case.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
ov[x_, y_] := BosonicMatrixElement[{x, y}, 1, "Basis" -> "Coherent"];
f = ov[be, al] ov[al, ga] /. Abs[al] -> Sqrt[al Conjugate[al]] /. Conjugate[al] -> u /. al -> v;
Table[Simplify[(D[Exp[u v] f, {v, m}, {u, m}] /. {u -> 0, v -> 0})/m! - Conjugate[BosonicMatrixElement[{m, 0}, DisplacementOperator[be]]] BosonicMatrixElement[{m, 0}, DisplacementOperator[ga]]], {m, 0, 6}]
```
Output: `{0, 0, 0, 0, 0, 0, 0}`. Matches: the weak action of P_m reproduces |m><m| between all coherent states.

### Three ordered characteristic functions and their Gaussian conversion
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: G-028, G-030 (s-ordered form $C(\lambda,s)=e^{s|\lambda|^2/2}C_W$)

The Weyl, normal and antinormal characteristic functions are traces of the displacement in three orders; they differ only by the Gaussian factors e^{-+|lambda|^2/2}, which is the s-ordered family.

$$ C_W(\lambda)=\operatorname{Tr}(\hat\rho \hat D(\lambda)),\quad C_N(\lambda)=\operatorname{Tr}(\hat\rho e^{\lambda \hat a^\dagger}e^{-\lambda^*\hat a}),\quad C_A(\lambda)=\operatorname{Tr}(\hat\rho e^{-\lambda^*\hat a}e^{\lambda \hat a^\dagger}),\quad C_W=e^{-|\lambda|^2/2}C_N=e^{|\lambda|^2/2}C_A. $$

Audit: Correct, including the atlas's correction of the printed Eq. (3.128c). Checked by: QF exact (the two ordered displacements hold for any rho under the trace); hand derivation that a creation-first product is the normal function, which differs from the antinormal one by e^{|lambda|^2}.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{DisplacementOperator[l, Infinity, "Ordering" -> "Normal"], DisplacementOperator[l, Infinity, "Ordering" -> "Antinormal"]}
```
Output: `{E^(l*SuperDagger[\[FormalA]1])**E^(-(\[FormalA]1*Conjugate[l]))/E^((l*Conjugate[l])/2), E^((l*Conjugate[l])/2)*E^(-(\[FormalA]1*Conjugate[l]))**E^(l*SuperDagger[\[FormalA]1])}`. Matches: D = e^{-|lambda|^2/2} e^{lambda a^dagger} e^{-lambda^* a} = e^{|lambda|^2/2} e^{-lambda^* a} e^{lambda a^dagger}, giving C_W = e^{-|lambda|^2/2} C_N = e^{|lambda|^2/2} C_A.

### Ordered moments as derivatives of characteristic functions
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: G-029

Ordered moments are read off as mixed derivatives of the matching ordered characteristic function at the origin.

$$ \langle(\hat a^\dagger)^m \hat a^n\rangle=\left.\partial_\lambda^m\partial_{(-\lambda^*)}^nC_N\right|_0,\quad\langle \hat a^m(\hat a^\dagger)^n\rangle=\left.\partial_\lambda^n\partial_{(-\lambda^*)}^mC_A\right|_0. $$

Audit: Correct. Checked by: hand derivation (Taylor expansion of the ordered exponentials under the trace).
QF code (exact): none. a Taylor-coefficient identity under a trace; no QF call produces it, and the coherent-basis BosonicMatrixElement leaves exponentials of a or a^dagger unevaluated.

### Complex Fourier transform between P and the normal characteristic function
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: G-031

P and the normal characteristic function form a complex Fourier pair with kernel e^{lambda alpha^* - lambda^* alpha} and normalization pi^{-2}.

$$ C_N(\lambda)=\int P(\alpha)e^{\lambda\alpha^*-\lambda^*\alpha}d^2\alpha,\quad P(\alpha)=\frac1{\pi^2}\int C_N(\lambda)e^{\lambda^*\alpha-\lambda\alpha^*}d^2\lambda. $$

Audit: Correct (sign and pi^{-2}). Checked by: exact WL computation (script wl/031.wl): the inverse transform of e^{-nbar|lambda|^2} gives e^{-|alpha|^2/nbar}/(pi nbar) and the forward transform returns e^{-nbar|lambda|^2}.
QF code (exact): none. QF has no P function or characteristic-function transform.

### Thermal characteristic functions and phase-space densities
Topic: Thermal states and traces
Appears in: G-032

A thermal state has Gaussian characteristic functions with widths nbar + 1, nbar, nbar + 1/2, a Gaussian P of width nbar, and a Gaussian Q of width nbar + 1.

$$ C_A=e^{-(\bar n+1)|\lambda|^2},\quad C_N=e^{-\bar n|\lambda|^2},\quad C_W=e^{-(\bar n+1/2)|\lambda|^2},\quad P(\alpha)=\frac{e^{-|\alpha|^2/\bar n}}{\pi\bar n},\quad Q(\alpha)=\frac{e^{-|\alpha|^2/(1+\bar n)}}{\pi(1+\bar n)}. $$

Audit: Correct. Checked by: QF exact (coherent diagonal of e^{-x n} for Q; Fock elements <m|e^{lambda a^dagger} e^{-lambda^* a}|m> = L_m(|lambda|^2) for m = 0..5) with the exact Laguerre generating sum; C_A, C_W by the conversion factors of the G-028 entry, P by the G-031 integral.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
q = (1 - E^-x)/Pi BosonicMatrixElement[{al, al}, Exp[-x ad ** a], "Basis" -> "Coherent"];
{FullSimplify[q /. x -> Log[1 + 1/nb], nb > 0],
 Table[Simplify[BosonicMatrixElement[{m, m}, Exp[l ad] ** Exp[-Conjugate[l] a]] - LaguerreL[m, l Conjugate[l]]], {m, 0, 5}],
 Simplify[Sum[(1 - t) t^m LaguerreL[m, z], {m, 0, Infinity}] /. t -> nb/(1 + nb), nb > 0]}
```
Output: `{1/(E^((al*Conjugate[al])/(1 + nb))*(Pi + nb*Pi)), {0, 0, 0, 0, 0, 0}, E^(-(nb*z))}`. Matches: Q = e^{-|alpha|^2/(1+nbar)}/(pi(1+nbar)) and C_N = e^{-nbar|lambda|^2}.

### Q and Wigner functions of coherent and number states
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: G-033

Coherent states have Gaussian Q and W centered at their amplitude; number states have a Poisson-shaped Q and a Laguerre-ringed W of sign (-1)^m at the origin.

$$ Q_{|\beta\rangle}(\alpha)=\pi^{-1}e^{-|\alpha-\beta|^2},\quad W_{|\beta\rangle}(\alpha)=\frac2\pi e^{-2|\alpha-\beta|^2},\quad Q_{|m\rangle}(\alpha)=\frac{e^{-|\alpha|^2}|\alpha|^{2m}}{\pi m!},\quad W_{|m\rangle}(\alpha)=\frac2\pi(-1)^m L_m(4|\alpha|^2)e^{-2|\alpha|^2}. $$

Audit: Correct. Checked by: QF exact (coherent overlaps, parity in the coherent basis, Fock elements of D(alpha) and D(2 alpha) for m = 0..5); W through the displaced-parity identity of the next entry.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
ov[x_, y_] := BosonicMatrixElement[{x, y}, 1, "Basis" -> "Coherent"];
{Simplify[ComplexExpand[ov[al, be] ov[be, al] - Exp[-Abs[al - be]^2] /. {al -> x1 + I y1, be -> x2 + I y2}]],
 Simplify[ComplexExpand[BosonicMatrixElement[{d, d}, Exp[I Pi ad ** a], "Basis" -> "Coherent"] - Exp[-2 Abs[d]^2] /. d -> x + I y]],
 Table[Simplify[Abs[BosonicMatrixElement[{m, 0}, DisplacementOperator[al]]]^2/Pi - Exp[-Abs[al]^2] Abs[al]^(2 m)/(Pi m!)], {m, 0, 5}],
 Table[Simplify[2/Pi (-1)^m BosonicMatrixElement[{m, m}, DisplacementOperator[2 al]] - 2/Pi (-1)^m LaguerreL[m, 4 Abs[al]^2] Exp[-2 Abs[al]^2]], {m, 0, 5}]}
```
Output: `{0, 0, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}}`. Matches all four.

### Wigner function as displaced parity
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: G-034

The Wigner function is 2/pi times the expectation of the parity operator displaced to alpha, and the displaced parity equals D(2 alpha) times the parity.

$$ W_{\hat\rho}(\alpha)=\frac2\pi\operatorname{Tr}[\hat\rho \hat D(\alpha)\hat\Pi \hat D^\dagger(\alpha)],\qquad\hat D(\alpha)\hat\Pi \hat D^\dagger(\alpha)=\hat D(2\alpha)\hat\Pi. $$

Audit: Correct (the Fock-sum form in the atlas is the same trace). Checked by: QF exact: Pi D(-alpha) Pi^{-1} = D(alpha) and D(alpha)^2 = D(2 alpha), so D(alpha) Pi D(-alpha) = D(alpha) D(alpha) Pi.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicExpSimplify[Exp[I Pi ad ** a] ** DisplacementOperator[-al, Infinity] ** Exp[-I Pi ad ** a]],
 BosonicExpSimplify[DisplacementOperator[al, Infinity] ** DisplacementOperator[al, Infinity]]}
```
Output: `{E^(-(\[FormalA]1*Conjugate[al]) + al*SuperDagger[\[FormalA]1]), E^(-2*\[FormalA]1*Conjugate[al] + 2*al*SuperDagger[\[FormalA]1])}`. Matches.

### Normally ordered correlation functions
Topic: Normal symbols and normal-order exponentials
Appears in: G-035

Ideal photodetection measures normally ordered field correlations; for one mode the normalized r-th order coherence is the r-th factorial moment over the mean to the r-th power.

$$ G^{(r)}=\langle \hat E^{(-)}(x_1)\cdots \hat E^{(-)}(x_r)\hat E^{(+)}(x_r)\cdots \hat E^{(+)}(x_1)\rangle,\qquad g^{(r)}=\langle(\hat a^\dagger)^r \hat a^r\rangle/\langle \hat n\rangle^r. $$

Conditions: $\langle\hat n\rangle\ne0$.
Audit: Correct (definition plus the ideal-absorber argument). Checked by: prior audit only for the detection model; the single-mode g^{(2)} case is checked exactly in the next entry.
QF code (exact): none. a definition of a measured quantity; QF's G1Correlation and G2Coherence act on finite Fock states only.

### Second-order coherence and Mandel Q
Topic: Normal symbols and normal-order exponentials
Appears in: G-036

g^{(2)} is the second factorial moment over the squared mean, equal to 1 plus the excess number variance; it is 1 for coherent light, 2 for thermal light and 1 - 1/m for a number state, and undefined for the vacuum.

$$ g^{(2)}(0)=\frac{\langle \hat n(\hat n-\hat I)\rangle}{\langle \hat n\rangle^2}=1+\frac{\operatorname{Var}\hat n-\langle \hat n\rangle}{\langle \hat n\rangle^2},\quad Q_M=\langle \hat n\rangle(g^{(2)}(0)-1);\quad g^{(2)}=1,\ 2,\ 1-1/m. $$

Audit: Correct. Checked by: QF exact (coherent-basis and symbolic-m Fock moments, exact thermal sums); hand derivation for the variance form and the vacuum 0/0.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}]; nn = ad ** a;
pm = nb^m/(1 + nb)^(m + 1);
th[op_] := Sum[pm BosonicMatrixElement[{m, m}, op], {m, 0, Infinity}, Assumptions -> nb > 0];
{FullSimplify[BosonicMatrixElement[{al, al}, ad ** ad ** a ** a, "Basis" -> "Coherent"]/BosonicMatrixElement[{al, al}, nn, "Basis" -> "Coherent"]^2],
 FullSimplify[th[ad ** ad ** a ** a]/th[nn]^2, nb > 0],
 Simplify[FunctionExpand[BosonicMatrixElement[{m, m}, ad ** ad ** a ** a]/BosonicMatrixElement[{m, m}, nn]^2]],
 BosonicNormalOrder[nn ** (nn - 1), {a, ad}]}
```
Output: `{1, 2, (-1 + m)/m, GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]**GeneralizedPower[NonCommutativeMultiply, \[FormalA]1, 2]}`. Matches: 1, 2, (m-1)/m, and n(n-1) = a^dagger^2 a^2.

### Atomic transition algebra and Jaynes-Cummings excitation number
Topic: SU(2) and passive mode mixing
Appears in: G-037

The atomic raising and lowering operators close on sigma_3, and the Jaynes-Cummings Hamiltonian conserves the total excitation number of field plus atom.

$$ [\hat\sigma_+,\hat\sigma_-]=\hat\sigma_3,\quad[\hat\sigma_3,\hat\sigma_\pm]=\pm2\hat\sigma_\pm,\quad[\hat H_{\rm JC},\hat n+|e\rangle\langle e|]=0. $$

Audit: Correct. Checked by: exact 2x2 atomic matrices with QF normal-ordered field entries.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}]; nn = ad ** a;
sp = {{0, 1}, {0, 0}}; sm = Transpose[sp]; s3 = {{1, 0}, {0, -1}}; id = IdentityMatrix[2]; pe = {{1, 0}, {0, 0}};
ncm[1, y_] := y; ncm[x_, 1] := x; ncm[x_, y_] := BosonicNormalOrder[x ** y, {a, ad}];
prod[x_, y_] := Total[Flatten[Outer[#1[[1]] #2[[1]] ncm[#1[[2]], #2[[2]]] (#1[[3]] . #2[[3]]) &, x, y, 1], 1]];
H = {{w, nn, id}, {w0/2, 1, s3}, {lam, a, sp}, {lam, ad, sm}};
M = {{1, nn, id}, {1, 1, pe}};
{sp . sm - sm . sp == s3, s3 . sp - sp . s3 == 2 sp, s3 . sm - sm . s3 == -2 sm, Simplify[Expand[prod[H, M] - prod[M, H]]]}
```
Output: `{True, True, True, {{0, 0}, {0, 0}}}`. Matches: the three Pauli relations hold and [H_JC, n + |e><e|] is the zero operator matrix.

### Exact resonant Jaynes-Cummings propagator
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-038

On resonance the Jaynes-Cummings propagator rotates each doublet {|e,m>, |g,m+1>} by the angle lambda t sqrt(m+1), which the operator-valued matrix with its square-root number functions encodes.

$$ e^{-i\lambda t(\hat a\hat\sigma_++\hat a^\dagger\hat\sigma_-)}=\begin{pmatrix}\cos(\lambda t\sqrt{\hat n+\hat I})&-i\hat a\,\dfrac{\sin(\lambda t\sqrt{\hat n})}{\sqrt{\hat n}}\\-i\hat a^\dagger\,\dfrac{\sin(\lambda t\sqrt{\hat n+\hat I})}{\sqrt{\hat n+\hat I}}&\cos(\lambda t\sqrt{\hat n})\end{pmatrix}. $$

Audit: Correct, including factor placement. Checked by: exact WL MatrixExp of the 2x2 generator on each invariant doublet (symbolic m), coupling from the QF Fock element; |g,0> is invariant (hand).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
g = BosonicMatrixElement[{m, m + 1}, a];
u = MatrixExp[-I lt {{0, g}, {g, 0}}];
claimed = {{Cos[lt Sqrt[m + 1]], -I BosonicMatrixElement[{m, m + 1}, a] Sin[lt Sqrt[m + 1]]/Sqrt[m + 1]},
   {-I BosonicMatrixElement[{m + 1, m}, ad] Sin[lt Sqrt[m + 1]]/Sqrt[m + 1], Cos[lt Sqrt[m + 1]]}};
{g, Simplify[ExpToTrig[u] - claimed, m >= 0]}
```
Output: `{Sqrt[1 + m], {{0, 0}, {0, 0}}}`. Matches: every entry of the claimed operator matrix equals the doublet rotation.

### Detuned Jaynes-Cummings doublet propagator
Topic: SU(2) and passive mode mixing
Appears in: G-039

Detuning turns each doublet into a spin-1/2 precession at the generalized Rabi frequency Omega_m, because the doubled doublet Hamiltonian squares to Omega_m^2.

$$ \hat h_m/\hbar=\frac12\begin{pmatrix}\Delta&2\lambda\sqrt{m+1}\\2\lambda\sqrt{m+1}&-\Delta\end{pmatrix},\quad\Omega_m=\sqrt{\Delta^2+4\lambda^2(m+1)},\quad e^{-i\hat h_mt/\hbar}=\hat I_2\cos(\Omega_mt/2)-i\frac{2\hat h_m}{\hbar\Omega_m}\sin(\Omega_mt/2). $$

Audit: Correct; dressed energies hbar omega(m+1/2) +- hbar Omega_m/2 (atlas's correction of printed Eq. 4.145 holds). Checked by: exact WL computation (script wl/039.wl): (2h)^2 = Omega_m^2 I, the closed-form propagator equals MatrixExp, and the eigenvalues are +-Omega_m/2.
QF code (exact): none. a 2x2 atomic block with no field operator left; nothing for QF's bosonic algebra to do.

### Dispersive Hamiltonian as conditional number rotation
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-040

In the dispersive limit the atom picks a number-dependent rotation of the field: e^{i chi t n} in the ground state and e^{-i chi t(n+1)} in the excited state.

$$ \hat H_{\rm eff}=\hbar\chi\bigl(|e\rangle\langle e|+\hat n\hat\sigma_3\bigr),\quad\hat U_{\rm eff}=|g\rangle\langle g|\,e^{i\chi t\hat n}+|e\rangle\langle e|\,e^{-i\chi t(\hat n+\hat I)}. $$

Audit: Correct (excited-state sign as in the atlas). Checked by: hand derivation (sigma_3 = -1 on g gives -chi n; +1 on e gives chi(n+1)).
QF code (exact): none. the propagator is diagonal in the atom-number basis; the field rotation it implies is the G-004 check.

### Raman number-controlled atomic rotation
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-041

The Raman coupling lambda n sigma_1 rotates the atom by an angle proportional to the photon number, equivalently it rotates the field in opposite senses on the two sigma_1 eigenstates.

$$ \hat U_I(t)=\cos(\lambda t\hat n)\hat I_A-i\sin(\lambda t\hat n)\hat\sigma_1=|+\rangle_x\langle+|e^{-i\lambda t\hat n}+|-\rangle_x\langle-|e^{i\lambda t\hat n}. $$

Audit: Correct. Checked by: exact WL computation (script wl/041.wl): on each number sector the exact 2x2 exponential equals both stated forms, symbolic m.
QF code (exact): none. no spin-boson layer; the operator acts as an exact 2x2 matrix in each number sector.

### Multiphoton Jaynes-Cummings couplings
Topic: Fock matrix elements and generating functions
Appears in: G-042

The two-photon and two-mode Jaynes-Cummings couplings connect |e,m> to |g,m+2> and |e,m,k> to |g,m+1,k+1> with the products of ladder square roots.

$$ \langle g,m+2|\hat a^{\dagger2}\hat\sigma_-|e,m\rangle=\sqrt{(m+1)(m+2)},\qquad\langle g,m+1,k+1|\hat a^\dagger\hat b^\dagger\hat\sigma_-|e,m,k\rangle=\sqrt{(m+1)(k+1)}. $$

Audit: Correct. Checked by: QF exact (symbolic Fock indices).
QF code (exact):
```wl
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
{BosonicMatrixElement[{m + 2, m}, ad1 ** ad1], BosonicMatrixElement[{{m + 1, k + 1}, {m, k}}, ad1 ** ad2]}
```
Output: `{Sqrt[FactorialPower[2 + m, 2]], Sqrt[1 + k]*Sqrt[1 + m]}`. Matches: sqrt((m+1)(m+2)) and sqrt((m+1)(k+1)).

### Passive two-port transform must be unitary
Topic: SU(2) and passive mode mixing
Appears in: G-043

A linear two-port keeps the output commutators bosonic only if its transfer matrix is unitary.

$$ [\hat a_{{\rm out},j},\hat a_{{\rm out},k}^\dagger]=(TT^\dagger)_{jk}\hat I,\quad\text{hence }TT^\dagger=I_2. $$

Audit: Correct. Checked by: QF exact (general symbolic 2x2 T).
QF code (exact):
```wl
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
t = {{t11, t12}, {t21, t22}};
out = t . {a1, a2}; outd = Conjugate[t] . {ad1, ad2};
Simplify[Table[BosonicNormalOrder[out[[j]] ** outd[[k]] - outd[[k]] ** out[[j]], {a1, ad1, a2, ad2}], {j, 2}, {k, 2}] - t . ConjugateTranspose[t]]
```
Output: `{{0, 0}, {0, 0}}`. Matches: the output commutator matrix equals T T^dagger entry by entry.

### Beam splitter with i-reflection: exponential and adjoint action
Topic: SU(2) and passive mode mixing
Appears in: G-044, G-086 (50:50 case $\vartheta=\pi/2$)

The symmetric beam splitter exp(i theta J_1) mixes the two modes with cos(theta/2) transmission and i sin(theta/2) reflection.

$$ \hat B_1(\vartheta)=e^{i\vartheta(\hat a^\dagger \hat b+\hat b^\dagger \hat a)/2},\quad\hat B_1^\dagger\begin{pmatrix}\hat a\\\hat b\end{pmatrix}\hat B_1=\begin{pmatrix}\cos(\vartheta/2)&i\sin(\vartheta/2)\\i\sin(\vartheta/2)&\cos(\vartheta/2)\end{pmatrix}\begin{pmatrix}\hat a\\\hat b\end{pmatrix}. $$

Audit: Correct. Checked by: QF exact (symbolic theta).
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}];
g = I th (ad ** b + bd ** a)/2;
FullSimplify[{BosonicExpSimplify[Exp[-g] ** a ** Exp[g]], BosonicExpSimplify[Exp[-g] ** b ** Exp[g]]}, th > 0]
```
Output: `{\[FormalA]1*Cos[th/2] + I*\[FormalA]2*Sin[th/2], \[FormalA]2*Cos[th/2] + I*\[FormalA]1*Sin[th/2]}`. Matches both rows.

### Real beam-splitter rotation
Topic: SU(2) and passive mode mixing
Appears in: G-045

The real beam splitter exp(i theta J_2) is a plain rotation of the mode pair by theta/2.

$$ \hat B_2(\vartheta)=e^{i\vartheta(\hat a^\dagger \hat b-\hat b^\dagger \hat a)/(2i)},\quad \hat B_2^\dagger\begin{pmatrix}\hat a\\\hat b\end{pmatrix}\hat B_2=\begin{pmatrix}c&s\\-s&c\end{pmatrix}\begin{pmatrix}\hat a\\\hat b\end{pmatrix},\quad c=\cos(\vartheta/2),\ s=\sin(\vartheta/2). $$

Audit: Correct. Checked by: QF exact (symbolic theta).
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}];
g = I th (ad ** b - bd ** a)/(2 I);
FullSimplify[{BosonicExpSimplify[Exp[-g] ** a ** Exp[g]], BosonicExpSimplify[Exp[-g] ** b ** Exp[g]]}, th > 0]
```
Output: `{\[FormalA]1*Cos[th/2] + \[FormalA]2*Sin[th/2], \[FormalA]2*Cos[th/2] - \[FormalA]1*Sin[th/2]}`. Matches both rows (c, s) and (-s, c).

### Schwinger SU(2) algebra, Casimir and J3 phase action
Topic: SU(2) and passive mode mixing
Appears in: G-046

Two modes realize SU(2) through the Schwinger generators, with Casimir J_0(J_0+1); a J_3 rotation gives the two modes opposite half-phases.

$$ \hat J_3=(\hat n_a-\hat n_b)/2,\ \hat J_0=(\hat n_a+\hat n_b)/2,\ [\hat J_i,\hat J_j]=i\epsilon_{ijk}\hat J_k,\ \hat{\mathbf J}^2=\hat J_0(\hat J_0+\hat I);\quad e^{-i\varphi\hat J_3}\hat a e^{i\varphi\hat J_3}=e^{i\varphi/2}\hat a,\ \hat b\mapsto e^{-i\varphi/2}\hat b. $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}]; v = {a, ad, b, bd};
j1 = (ad ** b + bd ** a)/2; j2 = (ad ** b - bd ** a)/(2 I); j3 = (ad ** a - bd ** b)/2; j0 = (ad ** a + bd ** b)/2;
cm[x_, y_] := BosonicNormalOrder[x ** y - y ** x, v];
{Simplify[{cm[j1, j2] - BosonicNormalOrder[I j3, v], cm[j2, j3] - BosonicNormalOrder[I j1, v], cm[j3, j1] - BosonicNormalOrder[I j2, v], cm[j0, j1], cm[j0, j2], cm[j0, j3],
  BosonicNormalOrder[j1 ** j1 + j2 ** j2 + j3 ** j3 - j0 ** j0 - j0, v]}],
 FullSimplify[{BosonicExpSimplify[Exp[-I ph j3] ** a ** Exp[I ph j3]], BosonicExpSimplify[Exp[-I ph j3] ** b ** Exp[I ph j3]]}, ph > 0]}
```
Output: `{{0, 0, 0, 0, 0, 0, 0}, {\[FormalA]1*E^((I/2)*ph), \[FormalA]2/E^((I/2)*ph)}}`. Matches: all three cyclic commutators, J_0 commutes with every J_i, J^2 = J_0(J_0+1), and e^{-i phi J_3} a e^{i phi J_3} = e^{i phi/2} a, b -> e^{-i phi/2} b.

### Beam-splitter output of a number state (binomial state)
Topic: SU(2) and passive mode mixing
Appears in: G-047

Sending N photons into one port of a beam splitter produces the binomial superposition over the two outputs.

$$ \hat B_1(\vartheta)|0,N\rangle=\sum_{k=0}^N\sqrt{\binom Nk}\,[i\sin(\vartheta/2)]^k[\cos(\vartheta/2)]^{N-k}|k,N-k\rangle. $$

Audit: Correct. Checked by: QF exact: B b^dagger B^dagger = i sin(theta/2) a^dagger + cos(theta/2) b^dagger and B|0,0> = |0,0> give B|0,N> = (i s a^dagger + c b^dagger)^N/sqrt(N!)|0,0>; Fock amplitudes for N = 3, symbolic theta; hand binomial expansion for general N.
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}];
g = I th (ad ** b + bd ** a)/2;
bdOut = FullSimplify[BosonicExpSimplify[Exp[g] ** bd ** Exp[-g]], th > 0];
nN = 3;
{bdOut, Table[FullSimplify[BosonicMatrixElement[{{k, nN - k}, {0, 0}}, (NonCommutativeMultiply @@ ConstantArray[bdOut, nN])/Sqrt[nN!]] - Sqrt[Binomial[nN, k]] (I Sin[th/2])^k Cos[th/2]^(nN - k)], {k, 0, nN}]}
```
Output: `{I*Sin[th/2]*SuperDagger[\[FormalA]1] + Cos[th/2]*SuperDagger[\[FormalA]2], {0, 0, 0, 0}}`. Matches.

### Coherent products stay products under passive mixing
Topic: SU(2) and passive mode mixing
Appears in: G-048

Passive mixing maps a product of coherent states to a product of coherent states with amplitudes transformed by the beam-splitter matrix.

$$ \hat B_1(\vartheta)|\alpha,\beta\rangle=|c\alpha+is\beta,\;is\alpha+c\beta\rangle. $$

Audit: Correct. Checked by: QF exact (conjugation of the two-mode displacement exponent, symbolic theta), with B|0,0> = |0,0> (hand).
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}];
g = I th (ad ** b + bd ** a)/2;
FullSimplify[BosonicExpSimplify[Exp[g] ** Exp[al ad - Conjugate[al] a + be bd - Conjugate[be] b] ** Exp[-g]], th > 0]
```
Output: `E^(I*Sin[th/2]*(\[FormalA]2*Conjugate[al] + \[FormalA]1*Conjugate[be] + be*SuperDagger[\[FormalA]1] + al*SuperDagger[\[FormalA]2]) + Cos[th/2]*(-(\[FormalA]1*Conjugate[al]) - \[FormalA]2*Conjugate[be] + al*SuperDagger[\[FormalA]1] + be*SuperDagger[\[FormalA]2]))`. Matches: the transformed displacement has amplitudes c alpha + i s beta on mode 1 and i s alpha + c beta on mode 2.

### Two-photon interference and twin-Fock output at 50:50
Topic: SU(2) and passive mode mixing
Appears in: G-049

At 50:50, two single photons leave together in one port (Hong-Ou-Mandel), and N photons in each port give only even occupations with the stated weights times i^N.

$$ \hat B_1|1,1\rangle=\frac i{\sqrt2}(|2,0\rangle+|0,2\rangle),\qquad \hat B_1|N,N\rangle=i^N\sum_{k=0}^N\frac{\sqrt{(2k)!(2N-2k)!}}{2^N k!(N-k)!}|2k,2N-2k\rangle. $$

Conditions: 50:50 ($\vartheta=\pi/2$).
Audit: Correct, with the overall phase i^N (the atlas writes "up to the overall phase i^N"). Checked by: QF exact for N = 1, 2, 3 (B a^dagger B^dagger and B b^dagger B^dagger at theta = pi/2, then two-mode Fock amplitudes); odd components vanish.
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}];
g = I Pi (ad ** b + bd ** a)/4;
adO = BosonicExpSimplify[Exp[g] ** ad ** Exp[-g]]; bdO = BosonicExpSimplify[Exp[g] ** bd ** Exp[-g]];
st[nN_] := (NonCommutativeMultiply @@ Join[ConstantArray[adO, nN], ConstantArray[bdO, nN]])/nN!;
out[nN_] := Table[Simplify[BosonicMatrixElement[{{2 k, 2 nN - 2 k}, {0, 0}}, st[nN]] - I^nN Sqrt[(2 k)! (2 nN - 2 k)!]/(2^nN k! (nN - k)!)], {k, 0, nN}];
{adO, bdO, out[1], out[2], out[3], Simplify[BosonicMatrixElement[{{1, 1}, {0, 0}}, st[1]]], Simplify[BosonicMatrixElement[{{3, 1}, {0, 0}}, st[2]]]}
```
Output: `{(SuperDagger[\[FormalA]1] + I*SuperDagger[\[FormalA]2])/Sqrt[2], (I*SuperDagger[\[FormalA]1] + SuperDagger[\[FormalA]2])/Sqrt[2], {0, 0}, {0, 0, 0}, {0, 0, 0, 0}, 0, 0}`. Matches for N = 1, 2, 3, and the |1,1> and |3,1> amplitudes are 0.

### Phase shifter and parity exponential
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-050, G-076 (Kerr at $t=\pi/K$ equals parity)

A phase shifter rotates a coherent amplitude; at angle pi it is the parity operator, which squares to the identity and maps |alpha> to |-alpha>.

$$ \hat R(\phi)=e^{i\phi \hat n},\quad \hat R(\phi)|\alpha\rangle=|\alpha e^{i\phi}\rangle,\quad\hat\Pi=e^{i\pi \hat n}=(-1)^{\hat n},\quad\hat\Pi^2=\hat I,\quad\hat\Pi|\alpha\rangle=|-\alpha\rangle. $$

Audit: Correct. The interferometer fringe sign remark (-cos 2 phi) is Unchecked: it depends on the book's Eq. (11.2) setup, not reproduced in the atlas; prior audit only. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}]; nn = ad ** a;
ov[x_, y_] := BosonicMatrixElement[{x, y}, 1, "Basis" -> "Coherent"];
{FullSimplify[BosonicMatrixElement[{be, al}, Exp[I ph nn], "Basis" -> "Coherent"] == ov[be, al Exp[I ph]], ph \[Element] Reals],
 Simplify[BosonicMatrixElement[{m, n}, Exp[I Pi nn]], n \[Element] Integers],
 Simplify[BosonicMatrixElement[{m, n}, Exp[2 I Pi nn]], n \[Element] Integers],
 Simplify[BosonicMatrixElement[{be, al}, Exp[I Pi nn], "Basis" -> "Coherent"] == ov[be, -al]]}
```
Output: `{True, (-1)^n*KroneckerDelta[m, n], KroneckerDelta[m, n], True}`. Matches: <beta|e^{i phi n}|alpha> = <beta|alpha e^{i phi}>, e^{i pi n} has Fock elements (-1)^n, e^{2 pi i n} is the identity on integer n, and <beta|Pi|alpha> = <beta|-alpha>.

### Single-mode squeeze operator
Topic: Single-mode squeezing
Appears in: G-051

The squeeze operator is the exponential of a quadratic anti-Hermitian generator, so its adjoint is the squeeze with opposite parameter, and it squeezes the vacuum.

$$ \hat S(\xi)=\exp\!\left[\tfrac12(\xi^*\hat a^2-\xi \hat a^{\dagger2})\right],\quad \hat S^\dagger(\xi)=\hat S(-\xi),\quad|\xi\rangle=\hat S(\xi)|0\rangle. $$

Audit: Correct; same sign convention as QF. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[xi, Infinity];
{s, BosonicExpSimplify[s ** SqueezeOperator[-xi, Infinity]] /. Power[E, x_] :> Power[E, Expand[x]]}
```
Output: `{E^((Conjugate[xi]*\[FormalA]1**\[FormalA]1 - xi*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]1])/2), 1}`. Matches: QF's symbolic S(xi) has the book's exponent and S(xi) S(-xi) = 1.

### Single-mode Bogoliubov transformation
Topic: Single-mode squeezing
Appears in: G-052

Squeezing transforms a into the Bogoliubov combination mu a - nu a^dagger with mu = cosh r and nu = e^{i theta} sinh r.

$$ \hat S^\dagger \hat a\hat S=\mu \hat a-\nu \hat a^\dagger,\quad\hat S^\dagger \hat a^\dagger \hat S=\mu \hat a^\dagger-\nu^*\hat a,\quad\mu=\cosh r,\ \nu=e^{i\theta}\sinh r. $$

Audit: Correct. Checked by: QF exact (symbolic r, theta).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
FullSimplify[{BosonicExpSimplify[sd ** a ** s], BosonicExpSimplify[sd ** ad ** s]}, r > 0 && th \[Element] Reals]
```
Output: `{\[FormalA]1*Cosh[r] - E^(I*th)*Sinh[r]*SuperDagger[\[FormalA]1], -((\[FormalA]1*Sinh[r])/E^(I*th)) + Cosh[r]*SuperDagger[\[FormalA]1]}`. Matches both lines.

### Quadratic SU(1,1) algebra of one mode
Topic: Single-mode squeezing
Appears in: G-053

Half the squares of the ladder operators, with half of n + 1/2, close into the su(1,1) algebra.

$$ \hat K_+=\tfrac12\hat a^{\dagger2},\ \hat K_-=\tfrac12\hat a^2,\ \hat K_0=\tfrac12(\hat n+\tfrac12\hat I),\quad[\hat K_0,\hat K_\pm]=\pm \hat K_\pm,\quad[\hat K_-,\hat K_+]=2\hat K_0. $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}]; v = {a, ad};
kp = ad ** ad/2; km = a ** a/2; k0 = (ad ** a + 1/2)/2;
cm[x_, y_] := BosonicNormalOrder[x ** y - y ** x, v];
{cm[k0, kp] - BosonicNormalOrder[kp, v], cm[k0, km] + BosonicNormalOrder[km, v], cm[km, kp] - BosonicNormalOrder[2 k0, v]}
```
Output: `{0, 0, 0}`. Matches all three commutators.

### Squeeze operator: normal-ordered Gauss factorization
Topic: Single-mode squeezing
Appears in: G-054

The squeeze operator factors into a normal-ordered product: an exponential of a^dagger^2, a power of sech r in the number operator, and an exponential of a^2.

$$ \hat S(\xi)=e^{-\tau \hat a^{\dagger2}/2}e^{-\ln(\cosh r)(\hat n+\tfrac12\hat I)}e^{\tau^*\hat a^2/2},\qquad\tau=e^{i\theta}\tanh r. $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{SqueezeOperator[xi, Infinity, "Ordering" -> "Normal"], BosonicExpOrder[Exp[c1 a ** a + c2 ad ** ad]]}
```
Output: `{E^(-1/2*(E^(I*Arg[xi])*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]1]*Tanh[Abs[xi]]))**Sech[Abs[xi]]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^((\[FormalA]1**\[FormalA]1*Tanh[Abs[xi]])/(2*E^(I*Arg[xi]))), E^((c2*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]*Tan[2*Sqrt[c1*c2]])/(2*Sqrt[c1*c2]))**Sec[2*Sqrt[c1*c2]]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^((c1*GeneralizedPower[NonCommutativeMultiply, \[FormalA]1, 2]*Tan[2*Sqrt[c1*c2]])/(2*Sqrt[c1*c2]))}`. Matches: Sech[r]^(1/2 + n) = e^{-ln(cosh r)(n + 1/2)} and tau = e^{i Arg xi} Tanh[Abs xi]; the general BosonicExpOrder form reduces to it with c1 = xi^*/2, c2 = -xi/2, sqrt(c1 c2) = i r/2, Tan[i r]/(i r) = Tanh[r]/r.

### Squeeze operator: reversed (antinormal) Gauss factorization
Topic: Single-mode squeezing
Appears in: G-055

The squeeze operator also factors in reversed order, with the middle number exponential changing sign.

$$ \hat S(\xi)=e^{\tau^*\hat a^2/2}e^{+\ln(\cosh r)(\hat n+\tfrac12\hat I)}e^{-\tau \hat a^{\dagger2}/2}. $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{SqueezeOperator[xi, Infinity, "Ordering" -> "Antinormal"], BosonicExpOrder[Exp[c1 a ** a + c2 ad ** ad], "Ordering" -> "Antinormal"]}
```
Output: `{E^((\[FormalA]1**\[FormalA]1*Tanh[Abs[xi]])/(2*E^(I*Arg[xi])))**Cosh[Abs[xi]]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^(-1/2*(E^(I*Arg[xi])*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]1]*Tanh[Abs[xi]])), E^((c1*GeneralizedPower[NonCommutativeMultiply, \[FormalA]1, 2]*Tan[2*Sqrt[c1*c2]])/(2*Sqrt[c1*c2]))**Sec[2*Sqrt[c1*c2]]^(1/2 - \[FormalA]1**SuperDagger[\[FormalA]1])**E^((c2*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]*Tan[2*Sqrt[c1*c2]])/(2*Sqrt[c1*c2]))}`. Matches: Cosh[r]^(1/2 + n) = e^{+ln(cosh r)(n + 1/2)}; the general antinormal form has Sec^(1/2 - a a^dagger), and a a^dagger = n + 1 turns it into the same Cosh^(1/2 + n) after c1 = xi^*/2, c2 = -xi/2.

### Squeezed vacuum: nullifier and Fock expansion
Topic: Single-mode squeezing
Appears in: G-056 (also G-062 at $\alpha=0$, G-063 overlap)

Squeezed vacuum is annihilated by the Bogoliubov combination mu a + nu a^dagger and contains only even photon numbers with the stated amplitudes.

$$ (\mu \hat a+\nu \hat a^\dagger)|\xi\rangle=0,\quad |\xi\rangle=\frac1{\sqrt{\cosh r}}\sum_{m\ge0}\frac{\sqrt{(2m)!}}{2^m m!}(-e^{i\theta}\tanh r)^m|2m\rangle,\quad P_{2m}=\frac{(2m)!\tanh^{2m}r}{2^{2m}(m!)^2\cosh r}. $$

Audit: Correct (including the exclusion of an odd solution). Checked by: QF exact (nullifier by conjugation; Fock amplitudes for 2m = 0..10 and odd ones zero) and exact WL normalization sum.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
mu = Cosh[r]; nu = E^(I th) Sinh[r];
null = FullSimplify[BosonicNormalOrder[mu BosonicExpSimplify[sd ** a ** s] + nu BosonicExpSimplify[sd ** ad ** s], {a, ad}], r > 0 && th \[Element] Reals];
coef[m_] := FullSimplify[BosonicMatrixElement[{2 m, 0}, SqueezeOperator[r E^(I th)]] - Sqrt[(2 m)!]/(2^m m! Sqrt[Cosh[r]]) (-E^(I th) Tanh[r])^m, r > 0 && th \[Element] Reals];
{null, Table[coef[m], {m, 0, 5}], Table[BosonicMatrixElement[{2 m + 1, 0}, SqueezeOperator[r E^(I th)]], {m, 0, 2}],
 FullSimplify[Sum[(2 m)! Tanh[r]^(2 m)/(2^(2 m) (m!)^2 Cosh[r]), {m, 0, Infinity}], r > 0]}
```
Output: `{\[FormalA]1, {0, 0, 0, 0, 0, 0}, {0, 0, 0}, 1}`. Matches: S^dagger (mu a + nu a^dagger) S = a, so the combination annihilates S|0>; the even amplitudes match, odd ones are 0, and the probabilities sum to 1.

### Rotated-quadrature variance of squeezed vacuum
Topic: Single-mode squeezing
Appears in: G-057, G-073 (collective mode)

The quadrature at angle vartheta of squeezed vacuum has variance (cosh 2r - sinh 2r cos(2 vartheta - theta))/4.

$$ \hat X(\vartheta)=\frac{\hat ae^{-i\vartheta}+\hat a^\dagger e^{i\vartheta}}2,\quad\operatorname{Var}_{\xi}\hat X(\vartheta)=\frac14[\cosh2r-\sinh2r\cos(2\vartheta-\theta)]. $$

Audit: Correct. Checked by: QF exact (Bogoliubov-transformed quadrature, vacuum expectation, symbolic r, theta, vartheta).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
sa = BosonicExpSimplify[sd ** a ** s]; sad = BosonicExpSimplify[sd ** ad ** s];
x = (E^(-I v) sa + E^(I v) sad)/2;
FullSimplify[BosonicVEV[x ** x, {a, ad}] - BosonicVEV[x, {a, ad}]^2 - (Cosh[2 r] - Sinh[2 r] Cos[2 v - th])/4, r > 0 && th \[Element] Reals && v \[Element] Reals]
```
Output: `0`. Matches.

### Quadrature covariance of squeezed vacuum
Topic: Single-mode squeezing
Appears in: G-058

The two Cartesian quadratures of squeezed vacuum are correlated by -sinh 2r sin theta/4, and the covariance determinant stays at the minimum 1/16.

$$ \operatorname{Cov}(\hat X_1,\hat X_2)=\tfrac12\langle\{\Delta \hat X_1,\Delta \hat X_2\}\rangle=-\frac14\sinh2r\sin\theta,\qquad\det V=\frac1{16}. $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
sa = BosonicExpSimplify[sd ** a ** s]; sad = BosonicExpSimplify[sd ** ad ** s];
x1 = (sa + sad)/2; x2 = (sa - sad)/(2 I);
ev[o_] := BosonicVEV[o, {a, ad}];
cov = ev[x1 ** x2 + x2 ** x1]/2 - ev[x1] ev[x2];
vm = {{ev[x1 ** x1] - ev[x1]^2, cov}, {cov, ev[x2 ** x2] - ev[x2]^2}};
FullSimplify[{cov + Sinh[2 r] Sin[th]/4, Det[vm] - 1/16}, r > 0 && th \[Element] Reals]
```
Output: `{0, 0}`. Matches both.

### Displaced squeezed state: low moments
Topic: Single-mode squeezing
Appears in: G-059

Displacing squeezed vacuum shifts the first moment by alpha, adds alpha^2 to <a^2>, and adds |alpha|^2 to the mean photon number.

$$ |\alpha,\xi\rangle=\hat D(\alpha)\hat S(\xi)|0\rangle,\quad\langle \hat a\rangle=\alpha,\quad\langle \hat a^2\rangle=\alpha^2-e^{i\theta}\sinh r\cosh r,\quad\langle \hat n\rangle=|\alpha|^2+\sinh^2r. $$

Audit: Correct. Checked by: QF exact (S^dagger D^dagger a D S = mu a - nu a^dagger + alpha, then vacuum expectations).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
tr = {a -> BosonicExpSimplify[sd ** a ** s] + al, ad -> BosonicExpSimplify[sd ** ad ** s] + Conjugate[al]};
ev[o_] := BosonicVEV[o /. tr, {a, ad}];
FullSimplify[{ev[a] - al, ev[a ** a] - (al^2 - E^(I th) Sinh[r] Cosh[r]), ev[ad ** a] - (al Conjugate[al] + Sinh[r]^2)}, r > 0 && th \[Element] Reals]
```
Output: `{0, 0, 0}`. Matches all three.

### Interchanging displacement and squeezing
Topic: Single-mode squeezing
Appears in: G-060

Squeezing and displacement commute up to a Bogoliubov change of the displacement amplitude, with no extra phase.

$$ \hat S(\xi)\hat D(\beta)=\hat D(\mu\beta-\nu\beta^*)\hat S(\xi),\quad \hat D(\alpha)\hat S(\xi)=\hat S(\xi)\hat D(\mu\alpha+\nu\alpha^*). $$

Audit: Correct, no extra phase. Checked by: QF exact (conjugation of a displacement by a squeeze, symbolic r, theta).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
mu = Cosh[r]; nu = E^(I th) Sinh[r];
expo[e_] := Expand[e /. Power[E, x_] :> x];
lhs1 = BosonicExpSimplify[s ** DisplacementOperator[be, Infinity] ** sd];
lhs2 = BosonicExpSimplify[sd ** DisplacementOperator[al, Infinity] ** s];
FullSimplify[{expo[lhs1] - expo[DisplacementOperator[mu be - nu Conjugate[be], Infinity]], expo[lhs2] - expo[DisplacementOperator[mu al + nu Conjugate[al], Infinity]]}, r > 0 && th \[Element] Reals]
```
Output: `{0, 0}`. Matches both exponents exactly.

### Displaced squeezed state as a Bogoliubov eigenstate
Topic: Single-mode squeezing
Appears in: G-061

A displaced squeezed state is the eigenstate of the Bogoliubov combination mu a + nu a^dagger with eigenvalue gamma = mu alpha + nu alpha^*, which fixes the rotated-quadrature eigenvalue.

$$ (\mu \hat a+\nu \hat a^\dagger)|\alpha,\xi\rangle=\gamma|\alpha,\xi\rangle,\quad\gamma=\mu\alpha+\nu\alpha^*,\qquad(\hat Y_1+ie^{-2r}\hat Y_2)|\alpha,\xi\rangle=e^{-r-i\theta/2}\gamma|\alpha,\xi\rangle. $$

Audit: Correct (the corrected rotated eigenvalue holds). Checked by: QF exact: (D S)^dagger (mu a + nu a^dagger) D S = a + gamma, and Y_1 + i e^{-2r} Y_2 = e^{-r - i theta/2}(mu a + nu a^dagger) as an operator identity.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
mu = Cosh[r]; nu = E^(I th) Sinh[r]; dg = al ad - Conjugate[al] a;
step1 = BosonicExpSimplify[Exp[-dg] ** (mu a + nu ad) ** Exp[dg]];
step2 = FullSimplify[BosonicNormalOrder[step1 /. {a -> BosonicExpSimplify[sd ** a ** s], ad -> BosonicExpSimplify[sd ** ad ** s]}, {a, ad}], r > 0 && th \[Element] Reals];
y1 = (a E^(-I th/2) + ad E^(I th/2))/2; y2 = (a E^(-I th/2) - ad E^(I th/2))/(2 I);
{step2, Simplify[BosonicNormalOrder[y1 + I E^(-2 r) y2 - E^(-r - I th/2) (mu a + nu ad), {a, ad}]]}
```
Output: `{\[FormalA]1 + al*Cosh[r] + E^(I*th)*Conjugate[al]*Sinh[r], 0}`. Matches.

### Fock coefficients of a displaced squeezed vacuum (Hermite form)
Topic: Fock matrix elements and generating functions
Appears in: G-062

The Fock amplitudes of a displaced squeezed vacuum are Hermite polynomials in gamma over sqrt(e^{i theta} sinh 2r), times a Gaussian prefactor.

$$ \langle m|\alpha,\xi\rangle=\frac{e^{-|\alpha|^2/2-\alpha^{*2}e^{i\theta}\tanh r/2}}{\sqrt{m!\cosh r}}\left(\frac{e^{i\theta}\tanh r}{2}\right)^{m/2}H_m\!\left(\frac{\gamma}{\sqrt{e^{i\theta}\sinh2r}}\right). $$

Audit: Correct (principal square roots used consistently). Checked by: QF exact normal forms of D and S and the braid e^{-alpha^* a} e^{-tau a^dagger^2/2} = e^{-tau (a^dagger - alpha^*)^2/2} e^{-alpha^* a}, so D S|0> = F(a^dagger)|0> with F(z) = e^{-|alpha|^2/2} cosh(r)^{-1/2} e^{alpha z - tau (z - alpha^*)^2/2}; the Taylor coefficients of F (m = 0..6) were compared with the Hermite form at three random exact points to 30 digits, because FullSimplify leaves a root-branch residual at m = 1 that is zero for r > 0.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
braid = BosonicExpBraid[Exp[-Conjugate[al] a] ** Exp[-tau/2 ad ** ad]];
f = Exp[-Abs[al]^2/2]/Sqrt[Cosh[r]] Exp[al z - tau/2 (z - Conjugate[al])^2];
coefs = Table[Sqrt[m!] SeriesCoefficient[f, {z, 0, m}], {m, 0, 6}] /. tau -> E^(I th) Tanh[r];
herm = Table[Exp[-Abs[al]^2/2 - Conjugate[al]^2 E^(I th) Tanh[r]/2]/Sqrt[m! Cosh[r]] (E^(I th) Tanh[r]/2)^(m/2) HermiteH[m, (Cosh[r] al + E^(I th) Sinh[r] Conjugate[al])/Sqrt[E^(I th) Sinh[2 r]]], {m, 0, 6}];
pts = {{r -> 3/7, th -> 6/5, al -> 2/5 - 7 I/10}, {r -> 6/5, th -> -2, al -> -1 + I/3}, {r -> 1/9, th -> 3, al -> 5/4 + 2 I}};
{DisplacementOperator[al, Infinity, "Ordering" -> "Normal"], SqueezeOperator[r E^(I th), Infinity, "Ordering" -> "Normal"], braid,
 Max[Abs[N[coefs /. pts, 30] - N[herm /. pts, 30]]]}
```
Output: ``` {E^(al*SuperDagger[\[FormalA]1])**E^(-(\[FormalA]1*Conjugate[al]))/E^((al*Conjugate[al])/2), E^(-1/2*(E^(I*Arg[E^(I*th)*r])*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]1]*Tanh[Abs[r]/E^Im[th]]))**Sech[Abs[r]/E^Im[th]]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^((\[FormalA]1**\[FormalA]1*Tanh[Abs[r]/E^Im[th]])/(2*E^(I*Arg[E^(I*th)*r]))), E^(-1/2*(tau*(-Conjugate[al] + SuperDagger[\[FormalA]1])**(-Conjugate[al] + SuperDagger[\[FormalA]1])))**E^(-(\[FormalA]1*Conjugate[al])), 0``29.803179072070474} ```. Matches: the largest difference over m = 0..6 and the three points is zero to 30 digits.

### Coherent overlap of squeezed vacuum and Q function of a displaced squeezed state
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: G-063

The coherent overlap of squeezed vacuum is a Gaussian in beta with the e^{i theta} tanh r quadratic term, and the Q function of the displaced state is that overlap squared, centered at alpha.

$$ \langle\beta|\hat S(\xi)|0\rangle=\frac1{\sqrt{\cosh r}}e^{-|\beta|^2/2-e^{i\theta}\beta^{*2}\tanh r/2},\quad Q(\beta)=\frac1{\pi\cosh r}\exp\{-|\delta|^2-\tanh r\,\operatorname{Re}(e^{-i\theta}\delta^2)\},\ \delta=\beta-\alpha. $$

Audit: Correct (the atlas's correction of printed Eq. 7.82). Checked by: exact WL computation (script wl/063.wl): the Fock sum of <beta|2m><2m|S|0> with the G-056 amplitudes gives the overlap, and |<beta - alpha|S|0>|^2/pi equals the stated Q; hand derivation of D^dagger(alpha)|beta> = phase |beta - alpha>.
QF code (exact): none. BosonicMatrixElement[{beta, 0}, Exp[c ad ** ad], "Basis" -> "Coherent"] returns unevaluated (exponentials in the coherent basis are a known missing feature).

### Wigner function of a displaced squeezed state
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: G-064

A displaced squeezed state has a Gaussian Wigner function centered at alpha, narrowed by e^{-r} along the squeezing axis at angle theta/2 and widened by e^{r} across it.

$$ W(\beta)=\frac2\pi\exp[-2e^{2r}x^2-2e^{-2r}y^2],\qquad x+iy=e^{-i\theta/2}(\beta-\alpha). $$

Audit: Correct (the atlas's correction of printed Eq. 7.83). Checked by: QF exact for C_W: S^dagger D(lambda) S = D(mu lambda + nu lambda^*) and <0|D(z)|0> = e^{-|z|^2/2}; exact WL Gaussian integral of the inverse transform (script wl/064.wl, standard formula for a positive quadratic form, determinant 1).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
{FullSimplify[BosonicExpSimplify[sd ** DisplacementOperator[l, Infinity] ** s], r > 0 && th \[Element] Reals], BosonicMatrixElement[{0, 0}, DisplacementOperator[z]]}
```
Output: `{E^(-(\[FormalA]1*Conjugate[l]*Cosh[r]) - (\[FormalA]1*l*Sinh[r])/E^(I*th) + l*Cosh[r]*SuperDagger[\[FormalA]1] + E^(I*th)*Conjugate[l]*Sinh[r]*SuperDagger[\[FormalA]1]), E^(-1/2*Abs[z]^2)}`. Matches: the inverse transform of C_W(lambda) = e^{lambda alpha^* - lambda^* alpha} e^{-|mu lambda + nu lambda^*|^2/2} equals the stated W exactly.

### Number variance of a displaced squeezed vacuum
Topic: Single-mode squeezing
Appears in: G-065

The number variance of a displaced squeezed vacuum is the coherent term |mu alpha - nu alpha^*|^2 plus the squeezed-vacuum term 2 sinh^2 r cosh^2 r.

$$ \operatorname{Var}\hat n=|\mu\alpha-\nu\alpha^*|^2+2\sinh^2r\cosh^2r. $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
tr = {a -> BosonicExpSimplify[sd ** a ** s] + al, ad -> BosonicExpSimplify[sd ** ad ** s] + Conjugate[al]};
ev[o_] := BosonicVEV[o /. tr, {a, ad}];
FullSimplify[ev[ad ** a ** ad ** a] - ev[ad ** a]^2 - ((Cosh[r] al - E^(I th) Sinh[r] Conjugate[al]) (Cosh[r] Conjugate[al] - E^(-I th) Sinh[r] al) + 2 Sinh[r]^2 Cosh[r]^2), r > 0 && th \[Element] Reals]
```
Output: `0`. Matches.

### Degenerate parametric amplifier gives a squeeze of parameter 2 eta t
Topic: Forced and parametric oscillators
Appears in: G-066, G-067 (four-wave mixing, same exponent with $\eta=\chi^{(3)}\beta^2$)

Under an undepleted classical pump at resonance the degenerate amplifier generates a single-mode squeeze with parameter 2 eta t.

$$ \hat H_I=i\hbar(\eta^*\hat a^2-\eta \hat a^{\dagger2}),\quad\hat U_I(t)=e^{t(\eta^*\hat a^2-\eta \hat a^{\dagger2})}=\hat S(2\eta t). $$

Conditions: undepleted classical pump, resonance.
Audit: Correct. Checked by: hand derivation, (xi^* a^2 - xi a^dagger^2)/2 at xi = 2 eta t is t(eta^* a^2 - eta a^dagger^2); QF exact on the exponent.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
Simplify[SqueezeOperator[2 eta t, Infinity] /. Power[E, x_] :> Power[E, Expand[x]], t > 0]
```
Output: `E^(-(eta*t*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]) + t*Conjugate[eta]*\[FormalA]1**\[FormalA]1)`. Matches: QF's S(2 eta t) exponent is t eta^* a^2 - t eta a^dagger^2 for real t. The quantized-pump Hamiltonians of G-066 and G-067 have no QF closed form.

### Two-mode squeeze operator and Bogoliubov action
Topic: Two-mode squeezing
Appears in: G-068

The two-mode squeeze mixes each mode with the conjugate of the other: a -> mu a - nu b^dagger and b -> mu b - nu a^dagger.

$$ \hat S_2(\xi)=e^{\xi^*\hat a\hat b-\xi \hat a^\dagger \hat b^\dagger},\quad\hat S_2^\dagger \hat a\hat S_2=\mu \hat a-\nu \hat b^\dagger,\quad\hat S_2^\dagger \hat b\hat S_2=\mu \hat b-\nu \hat a^\dagger. $$

Audit: Correct. Checked by: QF exact (symbolic r, theta).
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}];
xi = r E^(I th); g = Conjugate[xi] a ** b - xi ad ** bd;
FullSimplify[{BosonicExpSimplify[Exp[-g] ** a ** Exp[g]], BosonicExpSimplify[Exp[-g] ** b ** Exp[g]]}, r > 0 && th \[Element] Reals]
```
Output: `{\[FormalA]1*Cosh[r] - E^(I*th)*Sinh[r]*SuperDagger[\[FormalA]2], \[FormalA]2*Cosh[r] - E^(I*th)*Sinh[r]*SuperDagger[\[FormalA]1]}`. Matches both.

### Two-mode SU(1,1) algebra and conserved number difference
Topic: Two-mode squeezing
Appears in: G-069

The pair operators a^dagger b^dagger, a b and half the total number plus one close into su(1,1), and the number difference commutes with all of them.

$$ \hat K_+=\hat a^\dagger \hat b^\dagger,\ \hat K_-=\hat a\hat b,\ \hat K_0=(\hat n_a+\hat n_b+\hat I)/2,\quad[\hat K_0,\hat K_\pm]=\pm \hat K_\pm,\ [\hat K_-,\hat K_+]=2\hat K_0,\ [\hat n_a-\hat n_b,\hat K_\pm]=0. $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}]; v = {a, ad, b, bd};
kp = ad ** bd; km = a ** b; k0 = (ad ** a + bd ** b + 1)/2; dn = ad ** a - bd ** b;
cm[x_, y_] := BosonicNormalOrder[x ** y - y ** x, v];
{cm[k0, kp] - BosonicNormalOrder[kp, v], cm[k0, km] + BosonicNormalOrder[km, v], cm[km, kp] - BosonicNormalOrder[2 k0, v], cm[dn, kp], cm[dn, km]}
```
Output: `{0, 0, 0, 0, 0}`. Matches all five commutators.

### Two-mode squeeze: normal and reversed Gauss factorizations
Topic: Two-mode squeezing
Appears in: G-070, G-071 (reversed form)

The two-mode squeeze factors in normal order and in reversed order, the middle factor being a power of cosh r in the total number plus one.

$$ \hat S_2(\xi)=e^{-\tau \hat a^\dagger \hat b^\dagger}e^{-\ln(\cosh r)(\hat n_a+\hat n_b+\hat I)}e^{\tau^*\hat a\hat b}=e^{\tau^*\hat a\hat b}e^{+\ln(\cosh r)(\hat n_a+\hat n_b+\hat I)}e^{-\tau \hat a^\dagger \hat b^\dagger}. $$

Audit: Correct. Checked by: QF exact (normal and antinormal BosonicExpOrder).
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}];
{BosonicExpOrder[Exp[c1 a ** b + c2 ad ** bd]], BosonicExpOrder[Exp[c1 a ** b + c2 ad ** bd], "Ordering" -> "Antinormal"]}
```
Output: `{E^((c2*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]2]*Tan[Sqrt[c1*c2]])/Sqrt[c1*c2])**Sec[Sqrt[c1*c2]]^(SuperDagger[\[FormalA]1]**\[FormalA]1 + SuperDagger[\[FormalA]2]**\[FormalA]2)**E^((c1*\[FormalA]1**\[FormalA]2*Tan[Sqrt[c1*c2]])/Sqrt[c1*c2])*Sec[Sqrt[c1*c2]], E^((c1*\[FormalA]1**\[FormalA]2*Tan[Sqrt[c1*c2]])/Sqrt[c1*c2])**Sec[Sqrt[c1*c2]]^(-\[FormalA]1**SuperDagger[\[FormalA]1] - \[FormalA]2**SuperDagger[\[FormalA]2])**E^((c2*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]2]*Tan[Sqrt[c1*c2]])/Sqrt[c1*c2])*Sec[Sqrt[c1*c2]]}`. Matches after c1 = xi^*, c2 = -xi (sqrt(c1 c2) = i r, Tan[i r]/(i r) = Tanh[r]/r, Sec[i r] = Sech[r]); in the antinormal form a a^dagger + b b^dagger = n_a + n_b + 2 turns Sec^(-(...)) Sec into Cosh^(n_a + n_b + 1).

### Two-mode squeeze as two opposite single-mode squeezes of the sum and difference modes
Topic: Two-mode squeezing
Appears in: G-071

In the sum and difference modes the two-mode squeeze is a product of commuting single-mode squeezes with opposite parameters.

$$ \hat S_{2,ab}(\xi)=\hat S_c(\xi)\hat S_d(-\xi),\qquad\hat c=(\hat a+\hat b)/\sqrt2,\ \hat d=(\hat a-\hat b)/\sqrt2. $$

Audit: Correct. Checked by: QF exact (generator identity a b = (c^2 - d^2)/2 and commutation of the two single-mode generators).
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}]; v = {a, ad, b, bd};
c = (a + b)/Sqrt[2]; cd = (ad + bd)/Sqrt[2]; d = (a - b)/Sqrt[2]; dd = (ad - bd)/Sqrt[2];
xi = r E^(I th);
gc = (Conjugate[xi] c ** c - xi cd ** cd)/2; gd = (-Conjugate[xi] d ** d + xi dd ** dd)/2;
{BosonicNormalOrder[gc + gd - (Conjugate[xi] a ** b - xi ad ** bd), v], BosonicNormalOrder[gc ** gd - gd ** gc, v], BosonicNormalOrder[c ** dd - dd ** c, v], BosonicNormalOrder[c ** cd - cd ** c, v]}
```
Output: `{0, 0, 0, 1}`. Matches: the generators add up to the two-mode generator, they commute, and c, d^dagger commute while [c, c^dagger] = 1.

### Twin-beam vacuum: Fock expansion and nullifiers
Topic: Two-mode squeezing
Appears in: G-072

The twin-beam vacuum is a sum over equal photon pairs with weights (-e^{i theta} tanh r)^m/cosh r, annihilated by two Bogoliubov combinations and by the number difference.

$$ |\xi\rangle_2=\frac1{\cosh r}\sum_{m\ge0}(-e^{i\theta}\tanh r)^m|m,m\rangle,\quad(\mu \hat a+\nu \hat b^\dagger)|\xi\rangle_2=(\mu \hat b+\nu \hat a^\dagger)|\xi\rangle_2=0,\quad(\hat n_a-\hat n_b)|\xi\rangle_2=0. $$

Audit: Correct. Checked by: QF exact (nullifiers by conjugation, number difference commuting with the generator, and (a^dagger b^dagger)^m/m! |0,0> = |m,m>); the expansion follows from the normal form of the G-070 entry acting on the vacuum.
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}]; v = {a, ad, b, bd};
xi = r E^(I th); g = Conjugate[xi] a ** b - xi ad ** bd; mu = Cosh[r]; nu = E^(I th) Sinh[r]; tau = E^(I th) Tanh[r];
tf[o_] := BosonicExpSimplify[Exp[-g] ** o ** Exp[g]];
{FullSimplify[{BosonicNormalOrder[mu tf[a] + nu tf[bd], v], BosonicNormalOrder[mu tf[b] + nu tf[ad], v]}, r > 0 && th \[Element] Reals],
 BosonicNormalOrder[(ad ** a - bd ** b) ** g - g ** (ad ** a - bd ** b), v],
 Table[BosonicMatrixElement[{{m, m}, {0, 0}}, (NonCommutativeMultiply @@ Join[ConstantArray[ad, m], ConstantArray[bd, m]])/m!], {m, 1, 4}]}
```
Output: `{{\[FormalA]1, \[FormalA]2}, 0, {1, 1, 1, 1}}`. Matches.

### Twin-beam number statistics, collective squeezing, and the intermode Cauchy-Schwarz test
Topic: Two-mode squeezing
Appears in: G-073, G-091

Each twin-beam mode has mean sinh^2 r and Bose-Einstein variance, the two photon numbers are perfectly correlated, the normally ordered Cauchy-Schwarz inequality fails, and the collective quadratures are squeezed by e^{-2r}.

$$ \langle \hat n_a\rangle=N=\sinh^2r,\quad\operatorname{Var}\hat n_a=\operatorname{Cov}(\hat n_a,\hat n_b)=N(N+1),\quad\langle \hat a^{\dagger2}\hat a^2\rangle\langle \hat b^{\dagger2}\hat b^2\rangle=4N^4<(2N^2+N)^2=\langle \hat n_a\hat n_b\rangle^2. $$

Audit: Correct; the collective quadrature variances at theta = 0 are e^{-2r}/4 and e^{2r}/4. Checked by: QF exact (all moments and both collective variances through the Bogoliubov-transformed operators).
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}]; v = {a, ad, b, bd};
xi = r E^(I th); g = Conjugate[xi] a ** b - xi ad ** bd;
tf[o_] := BosonicExpSimplify[Exp[-g] ** o ** Exp[g]];
tr = {a -> tf[a], ad -> tf[ad], b -> tf[b], bd -> tf[bd]};
ev[o_] := BosonicVEV[o /. tr, v];
nN = Sinh[r]^2; na = ad ** a; nbb = bd ** b;
c = (a + b)/Sqrt[2]; cd = (ad + bd)/Sqrt[2]; x1 = (c + cd)/2; x2 = (c - cd)/(2 I);
FullSimplify[{ev[na] - nN, ev[nbb] - nN, ev[na ** na] - ev[na]^2 - nN (nN + 1), ev[na ** nbb] - ev[na] ev[nbb] - nN (nN + 1),
  ev[ad ** ad ** a ** a] ev[bd ** bd ** b ** b] - 4 nN^4, ev[na ** nbb] - (2 nN^2 + nN),
  (ev[x1 ** x1] - ev[x1]^2 - E^(-2 r)/4) /. th -> 0, (ev[x2 ** x2] - ev[x2]^2 - E^(2 r)/4) /. th -> 0}, r > 0 && th \[Element] Reals]
```
Output: `{0, 0, 0, 0, 0, 0, 0, 0}`. Matches all eight.

### Reduced twin-beam state is thermal; entanglement entropy
Topic: Thermal states and traces
Appears in: G-074, G-104 (Schmidt spectrum justification)

Tracing one twin-beam mode leaves a thermal state with tanh^2 r as the Boltzmann factor, whose entropy is the Bose-Einstein entropy of mean sinh^2 r.

$$ \hat\rho_a=\operatorname{sech}^2r\sum_{m\ge0}\tanh^{2m}r|m\rangle\langle m|,\quad e^{-\hbar\omega_a/(k_BT_{\rm eff})}=\tanh^2r,\quad S(\hat\rho_a)=(N+1)\ln(N+1)-N\ln N. $$

Audit: Correct. Checked by: QF exact amplitudes |c_m|^2 = sech^2 r tanh^{2m} r (Schmidt form, so the reduced state is diagonal with these weights) and exact WL entropy sum.
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}];
tau = E^(I th) Tanh[r];
cm = Table[BosonicMatrixElement[{{m, m}, {0, 0}}, (NonCommutativeMultiply @@ Join[ConstantArray[ad, m], ConstantArray[bd, m]])/m!] (-tau)^m/Cosh[r], {m, 1, 4}];
ent = -Sum[(1 - t) t^m (Log[1 - t] + m Log[t]), {m, 0, Infinity}, Assumptions -> 0 < t < 1];
{FullSimplify[Abs[cm]^2 - Table[Sech[r]^2 Tanh[r]^(2 m), {m, 1, 4}], r > 0 && th \[Element] Reals],
 FullSimplify[(ent /. t -> nN/(nN + 1)) - ((nN + 1) Log[nN + 1] - nN Log[nN]), nN > 0], FullSimplify[Tanh[r]^2 - Sinh[r]^2/(Sinh[r]^2 + 1)]}
```
Output: `{{0, 0, 0, 0}, 0, 0}`. Matches: the weights, S = (N+1) ln(N+1) - N ln N, and tanh^2 r = N/(N+1).

### Nondegenerate parametric amplifier gives a two-mode squeeze of parameter eta t
Topic: Forced and parametric oscillators
Appears in: G-075

Under an undepleted pump the nondegenerate amplifier generates a two-mode squeeze with parameter eta t, with no factor 2.

$$ \hat H_I=i\hbar(\eta^*\hat a\hat b-\eta \hat a^\dagger \hat b^\dagger),\quad \hat U_I=\hat S_2(\eta t). $$

Audit: Correct. Checked by: hand derivation (-i H_I t/hbar = t(eta^* a b - eta a^dagger b^dagger), no factor 1/2 in S_2).
QF code (exact): none. the statement is the one-line exponent identity; the exponential it names is the two-mode squeeze checked in the G-068 and G-070 entries.

### Self-Kerr evolution of a coherent state
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-076

Self-Kerr evolution multiplies each Fock amplitude by e^{-i K t m^2}; at t = pi/K this is the parity operator.

$$ \hat U_K(t)=e^{-iKt \hat n^2},\quad \hat U_K(t)|\alpha\rangle=e^{-|\alpha|^2/2}\sum_{m\ge0}\frac{\alpha^m}{\sqrt{m!}}e^{-iKtm^2}|m\rangle,\quad \hat U_K(\pi/K)=\hat\Pi. $$

Audit: Correct. Checked by: exact WL computation (script wl/075.wl): e^{-i pi m^2} = (-1)^m for integer m.
QF code (exact): none. BosonicMatrixElement[{m, n}, Exp[-I K t nn ** nn]] returns unevaluated (no rule for the square of a^dagger a in an exponent).

### Half-revival cat identity
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-077

At the half revival the Kerr phase is 1 on even and -i on odd photon numbers, so the operator is a combination of identity and parity and turns a coherent state into a cat.

$$ e^{-i\pi \hat n^2/2}=\frac{e^{-i\pi/4}}{\sqrt2}(\hat I+i\hat\Pi),\quad e^{-i\pi \hat n^2/2}|\alpha\rangle=\frac{e^{-i\pi/4}}{\sqrt2}(|\alpha\rangle+i|-\alpha\rangle). $$

Audit: Correct; confirms the atlas's correction of printed Eq. (7.142). Checked by: exact WL computation (script wl/076.wl): e^{-i pi m^2/2} is 1 for m = 2k and -i for m = 2k+1, the same values as e^{-i pi/4}(1 + i(-1)^m)/sqrt(2); the printed +i on odd m is wrong.
QF code (exact): none. same missing n^2 rule as the previous entry.

### Normally ordered Kerr differs by a rotation
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-078

The normally ordered Kerr generator n(n-1) differs from n^2 by a linear term, so its evolution differs by a free rotation and is periodic with period pi/K.

$$ \hat a^{\dagger2}\hat a^2=\hat n(\hat n-\hat I),\quad e^{-iKt \hat n(\hat n-\hat I)}=e^{-iKt \hat n^2}e^{iKt\hat n},\quad e^{-i\pi\hat n(\hat n-\hat I)}=\hat I. $$

Audit: Correct (period pi/K). Checked by: QF exact for the ordering; exact WL for the phases on integer m.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}]; nn = ad ** a;
{BosonicNormalOrder[nn ** (nn - 1) - ad ** ad ** a ** a, {a, ad}], FullSimplify[Exp[-I Pi m (m - 1)], m \[Element] Integers],
 Simplify[Exp[-I K t m (m - 1)] - Exp[-I K t m^2] Exp[I K t m]]}
```
Output: `{0, 1, 0}`. Matches: n(n-1) = a^dagger^2 a^2, e^{-i pi m(m-1)} = 1, and e^{-iKt m(m-1)} = e^{-iKt m^2} e^{iKt m}.

### Number-function conjugation and Kerr Heisenberg solution
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-079

Conjugating a by any function of n multiplies it by the phase of the difference f(n) - f(n-1); for the normally ordered Kerr this gives a(t) = a e^{-2iKt(n-1)} = e^{-2iKtn} a.

$$ e^{itf(\hat n)}\hat ae^{-itf(\hat n)}=\hat a\,e^{-it[f(\hat n)-f(\hat n-\hat I)]},\qquad\hat a(t)=\hat a\,e^{-2iKt(\hat n-\hat I)}=e^{-2iKt\hat n}\hat a\quad(f=K\hat n(\hat n-\hat I)). $$

Audit: Correct. Checked by: QF exact Fock elements with exact WL algebra for general f (both sides are compared on their only nonzero element <n-1|.|n>).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
amp = BosonicMatrixElement[{n - 1, n}, a];
lhs = Exp[I t f[n - 1]] amp Exp[-I t f[n]];
rhs = amp Exp[-I t (f[n] - f[n - 1])];
kerr = f -> Function[x, K x (x - 1)];
{amp, Simplify[lhs - rhs], Simplify[(lhs /. kerr) - amp Exp[-2 I K t (n - 1)]],
 Simplify[(lhs /. kerr) - BosonicMatrixElement[{n - 1, n - 1}, Exp[-2 I K t ad ** a]] amp], Simplify[(lhs /. kerr)/(amp Exp[-2 I K t n])]}
```
Output: `{Sqrt[n], 0, 0, 0, E^((2*I)*K*t)}`. Matches both ordered forms; the form a e^{-2iKt n} without the shift is off by the factor e^{2iKt}, as the atlas warns.

### Kerr coherent-state moments and quadrature variance
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-080

Under the normally ordered Kerr Hamiltonian a coherent state's moments <a^m> acquire the phase e^{-iKt m(m-1)} and the factor exp(N(e^{-2imKt} - 1)), and quadrature variances follow from the low moments.

$$ \langle \hat a^m\rangle_t=\alpha^m e^{-iKtm(m-1)}\exp\{N(e^{-2imKt}-1)\},\quad\operatorname{Var}\hat X(\vartheta)=\frac14\{1+2[\langle \hat n\rangle-|\langle \hat a\rangle|^2]+2\operatorname{Re}[e^{-2i\vartheta}(\langle \hat a^2\rangle-\langle \hat a\rangle^2)]\}. $$

Conditions: normally ordered Kerr generator $K\hat n(\hat n-\hat I)$, $N=|\alpha|^2$.
Audit: Correct. Checked by: QF exact Fock elements of a^m (m = 1, 2, 3) with exact WL sums over the coherent distribution; QF exact normal-ordered X(vartheta)^2 for the variance identity (hand: Var X = <X^2> - <X>^2 with that normal form).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
f[x_] := K x (x - 1);
cn[n_] := Exp[-rho^2/2] rho^n E^(I n p)/Sqrt[n!]; cnc[n_] := Exp[-rho^2/2] rho^n E^(-I n p)/Sqrt[n!];
term[m_] := FullSimplify[FunctionExpand[(cnc[n - m] cn[n] BosonicMatrixElement[{n - m, n}, NonCommutativeMultiply @@ ConstantArray[a, m]] Exp[I t (f[n - m] - f[n])]) /. n -> k + m], k >= 0 && rho > 0];
mom[m_] := Sum[Evaluate[term[m]], {k, 0, Infinity}];
xv = (a E^(-I v) + ad E^(I v))/2;
{Table[FullSimplify[mom[m] - (rho E^(I p))^m Exp[-I K t m (m - 1)] Exp[rho^2 (Exp[-2 I m K t] - 1)], rho > 0 && Element[{p, K, t}, Reals]], {m, 1, 3}],
 BosonicNormalOrder[xv ** xv, {a, ad}]}
```
Output: `{{0, 0, 0}, 1/4 + GeneralizedPower[NonCommutativeMultiply, \[FormalA]1, 2]/(4*E^((2*I)*v)) + (E^((2*I)*v)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2])/4 + SuperDagger[\[FormalA]1]**\[FormalA]1/2}`. Matches for m = 1, 2, 3; X(vartheta)^2 normal-orders to (e^{-2i vartheta} a^2 + e^{2i vartheta} a^dagger^2 + 2 a^dagger a + 1)/4, which gives the stated variance.

### Cross-Kerr conditional rotation
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-081

Cross-Kerr coupling rotates mode a by an angle proportional to the photon number in mode b.

$$ \hat U_{\rm CK}=e^{-iKt\hat n_a\hat n_b},\quad \hat U_{\rm CK}|\alpha\rangle_a|m\rangle_b=|\alpha e^{-iKtm}\rangle_a|m\rangle_b,\quad\hat U_{\rm CK}^\dagger \hat a\hat U_{\rm CK}=\hat a e^{-iKt\hat n_b}. $$

Audit: Correct. Checked by: QF exact on each number sector of b (n_b = m), hand derivation that the cross-Kerr operator acts there as e^{-iKtm n_a}.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{FullSimplify[BosonicExpSimplify[Exp[I K t m ad ** a] ** a ** Exp[-I K t m ad ** a]], K > 0 && t > 0 && m > 0],
 FullSimplify[BosonicMatrixElement[{be, al}, Exp[-I K t m ad ** a], "Basis" -> "Coherent"] == BosonicMatrixElement[{be, al Exp[-I K t m]}, 1, "Basis" -> "Coherent"], Element[{K, t, m}, Reals]]}
```
Output: `{\[FormalA]1/E^(I*K*m*t), True}`. Matches: e^{iKtm n} a e^{-iKtm n} = e^{-iKtm} a and <beta|e^{-iKtm n}|alpha> = <beta|alpha e^{-iKtm}>.

### Coherent-state cat: normalization and a^2 eigenvalue
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-082, G-085 (parity projection prepares the even/odd cat)

A superposition of |alpha> and |-alpha> with relative phase Phi has squared norm 2(1 + e^{-2|alpha|^2} cos Phi) and is an eigenstate of a^2 with eigenvalue alpha^2.

$$ |C_\Phi\rangle=\frac{|\alpha\rangle+e^{i\Phi}|-\alpha\rangle}{\sqrt{2(1+e^{-2|\alpha|^2}\cos\Phi)}},\qquad \hat a^2|C_\Phi\rangle=\alpha^2|C_\Phi\rangle. $$

Audit: Correct. Checked by: QF exact (coherent overlaps and coherent-basis a^2).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
ov[x_, y_] := BosonicMatrixElement[{x, y}, 1, "Basis" -> "Coherent"];
nrm = ov[al, al] + ov[-al, -al] + E^(I ph) ov[al, -al] + E^(-I ph) ov[-al, al];
{FullSimplify[nrm - 2 (1 + Exp[-2 Abs[al]^2] Cos[ph]), ph \[Element] Reals],
 Simplify[{BosonicMatrixElement[{be, al}, a ** a, "Basis" -> "Coherent"]/ov[be, al], BosonicMatrixElement[{be, -al}, a ** a, "Basis" -> "Coherent"]/ov[be, -al]}]}
```
Output: `{0, {al^2, al^2}}`. Matches.

### Cat photon-number distribution
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-083

The photon statistics of a cat are the Poisson weights times the interference factor |1 + e^{i Phi}(-1)^m|^2 over the normalization.

$$ |C_\pm\rangle=\frac{(\hat I\pm\hat\Pi)|\alpha\rangle}{\sqrt{2(1\pm e^{-2N})}},\quad P_m=\frac{e^{-N}N^m}{m!}\frac{|1+e^{i\Phi}(-1)^m|^2}{2(1+e^{-2N}\cos\Phi)}. $$

Audit: Correct. Checked by: QF exact Fock amplitudes of |alpha> and |-alpha> for m = 0..6, symbolic Phi; hand derivation for general m; (I +- Pi)|alpha> = |alpha> +- |-alpha> by the parity entry.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
amp[m_, x_] := BosonicMatrixElement[{m, 0}, DisplacementOperator[x]];
Table[FullSimplify[Abs[amp[m, al] + E^(I ph) amp[m, -al]]^2/(2 (1 + Exp[-2 Abs[al]^2] Cos[ph])) - Exp[-Abs[al]^2] Abs[al]^(2 m)/m! Abs[1 + E^(I ph) (-1)^m]^2/(2 (1 + Exp[-2 Abs[al]^2] Cos[ph])), ph \[Element] Reals], {m, 0, 6}]
```
Output: `{0, 0, 0, 0, 0, 0, 0}`. Matches.

### Cat Wigner function
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: G-084

The cat Wigner function is two Gaussians at +-alpha plus an interference fringe cos(4 alpha y + Phi) at the origin, over the cat normalization.

$$ W_{C_\Phi}(x,y)=\frac{e^{-2(x-\alpha)^2-2y^2}+e^{-2(x+\alpha)^2-2y^2}+2e^{-2x^2-2y^2}\cos(4\alpha y+\Phi)}{\pi(1+e^{-2\alpha^2}\cos\Phi)}. $$

Conditions: real $\alpha$, $\beta=x+iy$.
Audit: Correct (odd-cat denominator with a minus). Checked by: QF exact coherent overlaps, Weyl phase and parity action assembled into W = (2/pi)<C|D(2 beta) Pi|C>, then exact WL simplification for real alpha, symbolic Phi, x, y.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
ov[x_, y_] := BosonicMatrixElement[{x, y}, 1, "Basis" -> "Coherent"];
dp[x_, y_] := Exp[(2 be Conjugate[-y] - Conjugate[2 be] (-y))/2] ov[x, 2 be - y];
kets = {{1, al}, {E^(I ph), -al}};
w = 2/Pi Total[Flatten[Outer[Conjugate[#1[[1]]] #2[[1]] dp[#1[[2]], #2[[2]]] &, kets, kets, 1]]]/(2 (1 + Exp[-2 al^2] Cos[ph]));
target = (Exp[-2 (x - al)^2 - 2 y^2] + Exp[-2 (x + al)^2 - 2 y^2] + 2 Exp[-2 x^2 - 2 y^2] Cos[4 al y + ph])/(Pi (1 + Exp[-2 al^2] Cos[ph]));
FullSimplify[ComplexExpand[w /. be -> x + I y] - target, Element[{al, ph, x, y}, Reals]]
```
Output: `0`. Matches: the difference simplifies to 0.

### Cross-Kerr gate conjugated by a 50:50 beam splitter
Topic: SU(2) and passive mode mixing
Appears in: G-086

A cross-Kerr gate between b and c, conjugated by a 50:50 beam splitter on a and b, becomes a total-number phase times a c-controlled real rotation of a and b.

$$ \hat B^\dagger \hat U_K\hat B=e^{i\eta \hat n_c(\hat n_a+\hat n_b)/2}e^{(\eta \hat n_c/2)(\hat a^\dagger \hat b-\hat b^\dagger \hat a)},\qquad\hat B=e^{i\pi(\hat a^\dagger \hat b+\hat a\hat b^\dagger)/4},\ \hat U_K=e^{i\eta \hat n_b \hat n_c}. $$

Audit: Correct. Checked by: QF exact (three modes: B^dagger n_b B, the exponent identity, and the commutation of the two factors).
QF code (exact):
```wl
{a, ad, b, bd, c, cd} = FieldVariables[{1, 2, 3}]; v = {a, ad, b, bd, c, cd};
g = I Pi (ad ** b + a ** bd)/4;
nbT = BosonicNormalOrder[BosonicExpSimplify[Exp[-g] ** bd ** Exp[g]] ** BosonicExpSimplify[Exp[-g] ** b ** Exp[g]], v];
gen1 = I eta (ad ** a + bd ** b)/2; gen2 = eta (ad ** b - bd ** a)/2;
{Simplify[nbT], Simplify[BosonicNormalOrder[I eta nbT ** cd ** c - (gen1 + gen2) ** cd ** c, v]],
 BosonicNormalOrder[(gen1 ** cd ** c) ** (gen2 ** cd ** c) - (gen2 ** cd ** c) ** (gen1 ** cd ** c), v]}
```
Output: `{(SuperDagger[\[FormalA]1]**(\[FormalA]1 - I*\[FormalA]2) + SuperDagger[\[FormalA]2]**(I*\[FormalA]1 + \[FormalA]2))/2, 0, 0}`. Matches: B^dagger n_b B = (n_a + n_b)/2 - (i/2)(a^dagger b - b^dagger a), the exponent equals the stated sum, and the two generators commute.

### Linear quadrature exponential: normal-order conversion factor
Topic: Normal symbols and normal-order exponentials
Appears in: G-087

The exponential of a quadrature equals its normal-ordered form times e^{y^2/8}, the Gaussian factor of the central commutator.

$$ e^{y\Delta \hat X}={:e^{y\Delta \hat X}:}\,e^{y^2/8}. $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
BosonicExpOrder[Exp[y/2 (a + ad)]]
```
Output: `E^(y^2/8)*E^((y*SuperDagger[\[FormalA]1])/2)**E^((\[FormalA]1*y)/2)`. Matches (the shift by <X> is a scalar factor on both sides).

### Centered quadrature moments from normal moments
Topic: Normal symbols and normal-order exponentials
Appears in: G-088

Centered quadrature moments expand into normally ordered moments with Gaussian contraction weights N!/((N-2k)! k! 8^k); for a coherent state only the full contraction survives.

$$ \langle(\Delta \hat X)^N\rangle=\sum_{k=0}^{\lfloor N/2\rfloor}\frac{N!}{(N-2k)!k!8^k}\langle{: (\Delta \hat X)^{N-2k}:}\rangle,\qquad\langle(\Delta \hat X)^{2\ell}\rangle_{\rm coh}=(2\ell-1)!!/4^\ell. $$

Audit: Correct (fourth-order coefficients 3/2 and 3/16). Checked by: QF exact operator identity for N = 1..6 (the scalar shift by <X> leaves the commutator and hence the identity unchanged, hand); exact WL for the coherent values.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
x = (a + ad)/2;
nord[j_] := Total[Table[Binomial[j, i] NonCommutativeMultiply @@ Join[ConstantArray[ad, i], ConstantArray[a, j - i]], {i, 0, j}]]/2^j;
pw[j_] := NonCommutativeMultiply @@ ConstantArray[x, j];
rule[nN_] := BosonicNormalOrder[pw[nN] - Sum[nN!/((nN - 2 k)! k! 8^k) If[nN - 2 k == 0, 1, nord[nN - 2 k]], {k, 0, Floor[nN/2]}], {a, ad}];
{Table[rule[nN], {nN, 1, 6}], Table[(2 l)!/(l! 8^l) - (2 l - 1)!!/4^l, {l, 1, 5}]}
```
Output: `{{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}`. Matches: all six differences normal-order to 0 and (2l)!/(l! 8^l) = (2l-1)!!/4^l.

### Squared-field quadratures and their commutators
Topic: Single-mode squeezing
Appears in: G-089

The squared-field quadratures close with K_3 into su(1,1) with a state-dependent commutator, giving the bound Var Q_1 Var Q_2 >= 4 <K_3>^2.

$$ [\hat Q_1,\hat Q_2]=-4i\hat K_3,\quad[\hat K_3,\hat Q_1]=i\hat Q_2,\quad[\hat Q_2,\hat K_3]=i\hat Q_1,\qquad\operatorname{Var}\hat Q_1\operatorname{Var}\hat Q_2\ge4\langle \hat K_3\rangle^2. $$

Conditions: $\hat Q_1=(\hat a^{\dagger2}+\hat a^2)/2$, $\hat Q_2=(\hat a^{\dagger2}-\hat a^2)/(2i)$, $\hat K_3=(\hat n+\tfrac12\hat I)/2$.
Audit: Correct. Checked by: QF exact for the three commutators; the bound is Robertson's inequality with the first commutator.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}]; v = {a, ad};
q1 = (ad ** ad + a ** a)/2; q2 = (ad ** ad - a ** a)/(2 I); k3 = (ad ** a + 1/2)/2;
cm[x_, y_] := BosonicNormalOrder[x ** y - y ** x, v];
Simplify[{cm[q1, q2] - BosonicNormalOrder[-4 I k3, v], cm[k3, q1] - BosonicNormalOrder[I q2, v], cm[q2, k3] - BosonicNormalOrder[I q1, v]}]
```
Output: `{0, 0, 0}`. Matches all three.

### Pair coherent states
Topic: Two-mode squeezing
Appears in: G-090

Pair coherent states are eigenstates of a b with fixed number difference q; their amplitudes follow from a one-step recurrence and their normalization is a modified Bessel function.

$$ |\eta,q\rangle=\left[\sum_{m\ge0}\frac{|\eta|^{2m}}{m!(m+q)!}\right]^{-1/2}\sum_{m\ge0}\frac{\eta^m}{\sqrt{m!(m+q)!}}|m+q,m\rangle,\qquad\sum_{m\ge0}\frac{|\eta|^{2m}}{m!(m+q)!}=|\eta|^{-q}I_q(2|\eta|). $$

Audit: Correct. Checked by: QF exact (two-mode Fock elements with symbolic m, q) and exact WL Bessel sum.
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}];
el = BosonicMatrixElement[{{m + q, m}, {m + q + 1, m + 1}}, a ** b];
{el, FullSimplify[eta^(m + 1)/Sqrt[(m + 1)! (m + q + 1)!] el - eta eta^m/Sqrt[m! (m + q)!], m >= 0 && q >= 0],
 BosonicMatrixElement[{{m + q, m}, {m + q, m}}, ad ** a - bd ** b],
 FullSimplify[Sum[x^(2 m)/(m! (m + q)!), {m, 0, Infinity}] - x^(-q) BesselI[q, 2 x], x > 0 && q \[Element] Integers && q >= 0]}
```
Output: `{Sqrt[1 + m]*Sqrt[1 + m + q], 0, q, 0}`. Matches: the recurrence c_{m+1} sqrt((m+1)(m+q+1)) = eta c_m holds, n_a - n_b = q, and the sum equals |eta|^{-q} I_q(2|eta|).

### Frequency conversion as an SU(2) rotation
Topic: SU(2) and passive mode mixing
Appears in: G-092

Frequency conversion with U = e^{-tG} is a passive SU(2) rotation of modes b and c by angle g t with the pump phase on the cross terms.

$$ \hat U=e^{-t\hat G},\ \hat G=\eta\hat b^\dagger\hat c-\eta^*\hat b\hat c^\dagger,\quad\hat U^\dagger \hat b\hat U=\hat b\cos gt-e^{i\varphi}\hat c\sin gt,\quad\hat U^\dagger \hat c\hat U=\hat c\cos gt+e^{-i\varphi}\hat b\sin gt. $$

Conditions: $\hat H_I=-i\hbar\hat G$, $\eta=ge^{i\varphi}$.
Audit: Correct (the atlas's sign correction of the printed e^{+tG} holds). Checked by: QF exact (symbolic g, phi, t).
QF code (exact):
```wl
{b, bd, c, cd} = FieldVariables[{1, 2}];
eta = gg E^(I ph); gen = eta bd ** c - Conjugate[eta] b ** cd;
FullSimplify[{BosonicExpSimplify[Exp[t gen] ** b ** Exp[-t gen]], BosonicExpSimplify[Exp[t gen] ** c ** Exp[-t gen]]}, gg > 0 && t > 0 && ph \[Element] Reals]
```
Output: `{\[FormalA]1*Cos[gg*t] - \[FormalA]2*E^(I*ph)*Sin[gg*t], \[FormalA]2*Cos[gg*t] + (\[FormalA]1*Sin[gg*t])/E^(I*ph)}`. Matches both lines.

### Photon subtraction and addition
Topic: Coherent states (general)
Appears in: G-093

Subtracting a photon conditions the mean number to <n(n-1)>/<n>; adding a photon to a coherent state gives a displaced superposition of |1> and alpha^*|0>.

$$ \langle \hat n\rangle_- =\frac{\langle \hat n(\hat n-\hat I)\rangle}{\langle \hat n\rangle}=\langle \hat n\rangle-1+\frac{\operatorname{Var}\hat n}{\langle \hat n\rangle},\qquad|\alpha;+\rangle=\frac{\hat a^\dagger|\alpha\rangle}{\sqrt{1+|\alpha|^2}}=\frac{\hat D(\alpha)(|1\rangle+\alpha^*|0\rangle)}{\sqrt{1+|\alpha|^2}}. $$

Audit: Correct. Checked by: QF exact (a^dagger n a = n(n-1), D^dagger a^dagger D = a^dagger + alpha^*, norm <alpha|a a^dagger|alpha> = 1 + |alpha|^2); hand for the variance form.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}]; nn = ad ** a;
dg = al ad - Conjugate[al] a;
{BosonicNormalOrder[ad ** nn ** a - nn ** (nn - 1), {a, ad}], BosonicExpSimplify[Exp[-dg] ** ad ** Exp[dg]], FullSimplify[BosonicMatrixElement[{al, al}, a ** ad, "Basis" -> "Coherent"]]}
```
Output: `{0, Conjugate[al] + SuperDagger[\[FormalA]1], 1 + Abs[al]^2}`. Matches.

### Collective broadband commutator
Topic: Infinite modes and domains
Appears in: G-094

A collective mode summed over a band keeps a bosonic commutator scaled by C = sum |l_k|^2, and its quadratures have commutator iC/2.

$$ [\hat A(t),\hat A^\dagger(t)]=C\hat I,\quad C=\sum_{k\in\mathcal B}|l_k|^2,\quad[\hat X_1^{(c)},\hat X_2^{(c)}]=(iC/2)\hat I. $$

Conditions: $C$ finite.
Audit: Correct. Checked by: QF exact for three modes with symbolic weights and frequencies; hand derivation for a general finite band.
QF code (exact):
```wl
{a1, ad1, a2, ad2, a3, ad3} = FieldVariables[{1, 2, 3}]; v = {a1, ad1, a2, ad2, a3, ad3};
aA = l1 E^(-I w1 t) a1 + l2 E^(-I w2 t) a2 + l3 E^(-I w3 t) a3;
aAd = Conjugate[l1] E^(I w1 t) ad1 + Conjugate[l2] E^(I w2 t) ad2 + Conjugate[l3] E^(I w3 t) ad3;
x1 = (aA + aAd)/2; x2 = (aA - aAd)/(2 I);
Simplify[{BosonicNormalOrder[aA ** aAd - aAd ** aA, v], BosonicNormalOrder[x1 ** x2 - x2 ** x1, v]}]
```
Output: `{l1*Conjugate[l1] + l2*Conjugate[l2] + l3*Conjugate[l3], (I/2)*(l1*Conjugate[l1] + l2*Conjugate[l2] + l3*Conjugate[l3])}`. Matches.

### Low-gain down-conversion to second order
Topic: Two-mode squeezing
Appears in: G-095

At low gain the down-conversion propagator on vacuum is vacuum, minus i mu times one pair, minus mu^2 times two pairs, with the norm correction -|mu|^2/2 on vacuum.

$$ \hat U|0,0\rangle=(1-|\mu|^2/2)|0,0\rangle-i\mu|1,1\rangle-\mu^2|2,2\rangle+O(t^3),\qquad\mu=\eta t. $$

Conditions: $\hat H_I=\hbar(\eta\hat a_s^\dagger\hat a_i^\dagger+\eta^*\hat a_s\hat a_i)$.
Audit: Correct. Checked by: QF exact (two-mode Fock amplitudes of the first- and second-order Taylor terms of e^{-i(mu a^dagger b^dagger + mu^* a b)}; the third-order term is cubic in mu by construction).
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}];
h = mu ad ** bd + Conjugate[mu] a ** b;
term[k_] := (-I)^k (NonCommutativeMultiply @@ ConstantArray[h, k])/k!;
Table[Simplify[BosonicMatrixElement[{{j, j}, {0, 0}}, term[k]]], {k, 1, 2}, {j, 0, 2}]
```
Output: `{{0, (-I)*mu, 0}, {-1/2*(mu*Conjugate[mu]), 0, -mu^2}}`. Matches: first order gives -i mu on |1,1>, second order gives -|mu|^2/2 on |0,0> and -mu^2 on |2,2>.

### No-jump number exponential on a coherent state
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-096

A no-jump number exponential shrinks a coherent amplitude by e^{-s} and multiplies the state by e^{-|alpha|^2(1 - e^{-2s})/2}.

$$ e^{-s\hat n}|\alpha\rangle=\exp[-\tfrac12|\alpha|^2(1-e^{-2s})]|\alpha e^{-s}\rangle,\qquad\|e^{-s\hat n}|\alpha\rangle\|^2=e^{-|\alpha|^2(1-e^{-2s})}. $$

Audit: Correct. Checked by: QF exact (coherent-basis matrix elements).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{FullSimplify[BosonicMatrixElement[{be, al}, Exp[-s ad ** a], "Basis" -> "Coherent"] == Exp[-Abs[al]^2 (1 - Exp[-2 s])/2] BosonicMatrixElement[{be, al Exp[-s]}, 1, "Basis" -> "Coherent"], s > 0],
 FullSimplify[BosonicMatrixElement[{al, al}, Exp[-2 s ad ** a], "Basis" -> "Coherent"], s > 0]}
```
Output: `{True, E^(al*(-1 + E^(-2*s))*Conjugate[al])}`. Matches both.

### Normal ordering of a number exponential and the vacuum projector
Topic: Normal symbols and normal-order exponentials
Appears in: G-097

The exponential of lambda n is the normal-ordered exponential of (e^lambda - 1) n; at lambda -> -infinity this is the vacuum projector :e^{-n}:.

$$ e^{\lambda \hat n}={:\exp[(e^\lambda-1)\hat n]:},\qquad|0\rangle\langle0|={:e^{-\hat n}:}. $$

Audit: Correct. Checked by: QF exact; exact WL binomial sum for the Fock eigenvalues, and Fock elements of :e^{-n}: for m = 0..5.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicExpOrder[Exp[lam ad ** a]], Simplify[Sum[Binomial[m, k] (E^lam - 1)^k, {k, 0, m}]],
 Table[1 + Sum[(-1)^k/k! BosonicMatrixElement[{m, m}, NonCommutativeMultiply @@ Join[ConstantArray[ad, k], ConstantArray[a, k]]], {k, 1, m}], {m, 0, 5}]}
```
Output: `{NormalOrdered[E^((-1 + E^lam)*SuperDagger[\[FormalA]1]**\[FormalA]1)], E^(lam*m), {1, 0, 0, 0, 0, 0}}`. Matches: NormalOrdered[E^((E^lam - 1) n)], sum_k binom(m,k)(e^lam - 1)^k = e^{lam m}, and :e^{-n}: has eigenvalues 1, 0, 0, 0, 0, 0.

### Thermal state in normally ordered exponential form
Topic: Thermal states and traces
Appears in: G-097

Combining the previous entry with the thermal exponential writes the thermal state as a normal-ordered Gaussian with weight 1/(1 + nbar).

$$ \hat\rho_{\rm th}=\frac1{1+\bar n}{:\exp[-\hat n/(1+\bar n)]:}. $$

Audit: Correct. Checked by: QF exact (BosonicExpOrder and the substitution x = ln(1 + 1/nbar)).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicExpOrder[Exp[-x ad ** a]], Simplify[{1 - E^-x, E^-x - 1} /. x -> Log[1 + 1/nb], nb > 0]}
```
Output: `{NormalOrdered[E^((-1 + E^(-x))*SuperDagger[\[FormalA]1]**\[FormalA]1)], {(1 + nb)^(-1), -(1 + nb)^(-1)}}`. Matches: 1 - e^{-x} = 1/(1+nbar) and e^{-x} - 1 = -1/(1+nbar).

### Amplitude damping: Kraus form and matrix-element solution
Topic: Thermal states and traces
Appears in: G-098

Vacuum amplitude damping is a Kraus map with operators sqrt((1-T)^l/l!) T^{n/2} a^l, whose Fock components give the stated matrix-element sum and which solves the damping master equation.

$$ \rho_{mn}(t)=T^{(m+n)/2}\sum_{\ell\ge0}\frac{(1-T)^\ell}{\ell!}\sqrt{\frac{(m+\ell)!(n+\ell)!}{m!n!}}\,\rho_{m+\ell,n+\ell}(0),\qquad\hat M_\ell=\sqrt{(1-T)^\ell/\ell!}\,T^{\hat n/2}\hat a^\ell,\ \sum_\ell\hat M_\ell^\dagger\hat M_\ell=\hat I. $$

Conditions: $T=e^{-\gamma t}$, $\dot{\hat\rho}=\gamma(\hat a\hat\rho\hat a^\dagger-\{\hat n,\hat\rho\}/2)$.
Audit: Correct. Checked by: QF exact Fock elements (completeness for m = 0..5, <m|a^l|m+l>^2 = (m+l)!/m! for l = 1..3) and exact WL differentiation of the component coefficients against the Lindblad component equation (symbolic m, n, l).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
fac[l_] := NonCommutativeMultiply @@ Join[ConstantArray[ad, l], ConstantArray[a, l]];
low[l_] := NonCommutativeMultiply @@ ConstantArray[a, l];
cf[m_, n_, l_] := tT^((m + n)/2) (1 - tT)^l/l! Sqrt[(m + l)! (n + l)!/(m! n!)];
{Table[Simplify[tT^m + Sum[(1 - tT)^l/l! tT^(m - l) BosonicMatrixElement[{m, m}, fac[l]], {l, 1, m}]], {m, 0, 5}],
 Table[FullSimplify[FunctionExpand[BosonicMatrixElement[{m, m + l}, low[l]]^2 - (m + l)!/m!], m >= 0], {l, 1, 3}],
 FullSimplify[{D[cf[m, n, l] /. tT -> E^(-ga t), t] - ga (Sqrt[(m + 1) (n + 1)] cf[m + 1, n + 1, l - 1] - (m + n)/2 cf[m, n, l]) /. tT -> E^(-ga t),
   D[cf[m, n, 0] /. tT -> E^(-ga t), t] + ga (m + n)/2 cf[m, n, 0] /. tT -> E^(-ga t)}, m >= 0 && n >= 0 && l >= 1 && ga > 0 && t > 0]}
```
Output: `{{1, 1, 1, 1, 1, 1}, {0, 0, 0}, {0, 0}}`. Matches: completeness, the component amplitudes, and d/dt of each coefficient equals the master-equation right side.

### Damping of coherent dyads and of the normal characteristic function
Topic: Thermal states and traces
Appears in: G-099

Damping maps a coherent dyad to the dyad of shrunken amplitudes times a coherence factor that decays with the squared distance between alpha and beta.

$$ C_N(\lambda,t)=C_N(\sqrt T\lambda,0),\qquad\widehat\Lambda_T(|\alpha\rangle\langle\beta|)=\exp\!\left[-\frac{1-T}{2}(|\alpha|^2+|\beta|^2-2\beta^*\alpha)\right]|\sqrt T\alpha\rangle\langle\sqrt T\beta|. $$

Audit: Correct. Checked by: QF exact (T^{n/2}|alpha> = e^{-(1-T)|alpha|^2/2} |sqrt(T) alpha> in the coherent basis) with the exact WL sum over Kraus indices; the C_N statement follows by optical equivalence (hand).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
ov[x_, y_] := BosonicMatrixElement[{x, y}, 1, "Basis" -> "Coherent"];
k[x_] := FullSimplify[BosonicMatrixElement[{u, x}, Exp[Log[tT]/2 ad ** a], "Basis" -> "Coherent"]/ov[u, Sqrt[tT] x], 0 < tT < 1];
{k[al], FullSimplify[k[al] Conjugate[k[be]] Sum[(1 - tT)^l/l! (al Conjugate[be])^l, {l, 0, Infinity}] - Exp[-(1 - tT)/2 (Abs[al]^2 + Abs[be]^2 - 2 Conjugate[be] al)], 0 < tT < 1]}
```
Output: `{E^((al*(-1 + tT)*Conjugate[al])/2), 0}`. Matches.

### Driven damped oscillator: coherent amplitude and steady state
Topic: Forced and parametric oscillators
Appears in: G-100

A driven damped oscillator keeps coherent states coherent; the amplitude relaxes to -2iG/gamma.

$$ \dot\alpha=-iG-\gamma\alpha/2,\quad\alpha(t)=\alpha_0e^{-\gamma t/2}-\frac{2iG}{\gamma}(1-e^{-\gamma t/2}),\quad\alpha_{\rm ss}=-2iG/\gamma. $$

Audit: Correct (factor 2 versus the printed -iG/gamma). Checked by: QF exact coherent-basis matrix elements of the Lindblad generator applied to |beta><beta| at beta = -2iG/gamma (zero between all coherent bras and kets), and exact WL DSolve for the amplitude equation.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}]; nn = ad ** a;
me[x_, op_, y_] := BosonicMatrixElement[{x, y}, op, "Basis" -> "Coherent"];
bb = -2 I gG/ga;
lrho = -I gG (me[u, a + ad, bb] me[bb, 1, v] - me[u, 1, bb] me[bb, a + ad, v]) + ga (me[u, a, bb] me[bb, ad, v] - (me[u, nn, bb] me[bb, 1, v] + me[u, 1, bb] me[bb, nn, v])/2);
sol = alpha[t] /. First[DSolve[{alpha'[t] == -I gG - ga alpha[t]/2, alpha[0] == a0}, alpha[t], t]];
{FullSimplify[lrho, gG > 0 && ga > 0], FullSimplify[sol - (a0 E^(-ga t/2) - 2 I gG/ga (1 - E^(-ga t/2)))], Limit[sol, t -> Infinity, Assumptions -> ga > 0]}
```
Output: `{0, 0, ((-2*I)*gG)/ga}`. Matches: the generator vanishes on the steady state, and the ODE solution and its limit are as stated.

### Recoil exponential as a displacement and its interaction-picture form
Topic: Displacement
Appears in: G-101

The recoil exponential is a displacement by -i eta, and the free oscillator rotation turns its generator into the time-dependent quadrature.

$$ e^{-i\eta(\hat a+\hat a^\dagger)}=\hat D(-i\eta)=e^{-\eta^2/2}e^{-i\eta \hat a^\dagger}e^{-i\eta \hat a},\qquad\hat U_0^\dagger e^{-i\eta(\hat a+\hat a^\dagger)}\hat U_0=\exp[-i\eta(\hat ae^{-i\nu t}+\hat a^\dagger e^{i\nu t})]. $$

Conditions: real $\eta$.
Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}]; nn = ad ** a;
expo[e_] := Expand[e /. Power[E, x_] :> x];
{DisplacementOperator[-I eta, Infinity, "Ordering" -> "Normal"], BosonicExpOrder[Exp[-I eta (a + ad)]],
 FullSimplify[ExpToTrig[expo[BosonicExpSimplify[Exp[I nu t nn] ** Exp[-I eta (a + ad)] ** Exp[-I nu t nn]]] + I eta (a E^(-I nu t) + ad E^(I nu t))], nu > 0 && t > 0]}
```
Output: `{E^((-I)*eta*SuperDagger[\[FormalA]1])**E^((-I)*\[FormalA]1*Conjugate[eta])/E^((eta*Conjugate[eta])/2), E^((-I)*eta*SuperDagger[\[FormalA]1])**E^((-I)*\[FormalA]1*eta)/E^(eta^2/2), 0}`. Matches for real eta (and BosonicExpOrder gives the same normal form for any eta); the interaction-picture exponent difference is 0.

### Second-order Dyson elimination of a fast oscillating coupling
Topic: Wick, Dyson and auxiliary fields
Appears in: G-102

To second order a fast oscillating coupling leaves a secular term linear in t proportional to [A, A^dagger], which defines the effective Hamiltonian.

$$ \hat U_I\approx \hat I-\frac{ig^2t}{\Delta}[\hat A,\hat A^\dagger],\qquad\hat H_{\rm eff}=\frac{\hbar g^2}{\Delta}[\hat A,\hat A^\dagger]. $$

Audit: Correct (sign of the secular term). Checked by: exact WL computation (script wl/100.wl): the time-ordered double integrals of e^{+-i Delta(t1 - t2)} have secular parts +-i t/Delta and those of e^{+-i Delta(t1 + t2)} are bounded, giving -i g^2 t [A, A^dagger]/Delta.
QF code (exact): none. a time-ordered perturbative statement; QF has no Dyson-series tool.

### Effective Kerr Hamiltonian from an off-resonant parametric coupling
Topic: Nonlinear exponentials (number shift, Kerr, cats)
Appears in: G-103

For A = i a^dagger^2 b the commutator [A, A^dagger] is a polynomial in the two number operators, so the effective Hamiltonian is a sum of self-Kerr, cross-Kerr and pump-rotation terms.

$$ [\hat A,\hat A^\dagger]=\hat n_a(\hat n_a-\hat I)-(4\hat n_a+2\hat I)\hat n_b,\qquad\hat A=i\hat a^{\dagger2}\hat b. $$

Audit: Correct; the resulting H_eff follows by the G-102 formula. Checked by: QF exact.
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}]; v = {a, ad, b, bd};
aA = I ad ** ad ** b; aAd = -I bd ** a ** a;
na = ad ** a; nbb = bd ** b;
BosonicNormalOrder[aA ** aAd - aAd ** aA - (na ** (na - 1) - (4 na + 2) ** nbb), v]
```
Output: `0`. Matches: the difference normal-orders to 0.

### Coherent trace formula and Schmidt spectrum
Topic: Thermal states and traces
Appears in: G-104

The trace of an operator is the integral of its coherent diagonal over d^2 alpha/pi, and the entropy of a pure bipartite state is the Shannon entropy of its squared Schmidt coefficients.

$$ \operatorname{Tr}\hat F=\int\langle\alpha|\hat F|\alpha\rangle\frac{d^2\alpha}{\pi},\quad\hat\rho(t)=\hat U\hat\rho(0)\hat U^\dagger,\quad S=-\sum_j|g_j|^2\ln|g_j|^2. $$

Conditions: trace-class $\hat F$; orthonormal Schmidt families.
Audit: Correct with the stated Schmidt premise. Checked by: QF exact coherent diagonal of e^{-x n} with an exact WL integral against the Fock trace; hand derivation (G-019 resolution) for general trace-class F and for the Schmidt statement.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
d2 = FullSimplify[BosonicMatrixElement[{al, al}, Exp[-x ad ** a], "Basis" -> "Coherent"] /. al -> rho E^(I phi), rho > 0 && phi \[Element] Reals];
{d2, Integrate[2 rho d2, {rho, 0, Infinity}, Assumptions -> x > 0], Sum[E^(-x m), {m, 0, Infinity}]}
```
Output: `{E^((-1 + E^(-x))*rho^2), E^x/(-1 + E^x), E^x/(-1 + E^x)}`. Matches: both sides equal e^x/(e^x - 1).

### Continuum CCR and first-order pair-creation kernel
Topic: Infinite modes and domains
Appears in: G-105

Continuum modes have a delta-function commutator, and to first order a pump photon splits into a pair constrained by energy and momentum conservation.

$$ [\hat a(\mathbf k),\hat a^\dagger(\mathbf k')]=\delta^{(3)}(\mathbf k-\mathbf k')\hat I,\qquad\int d^3k_s\,d^3k_i\,\delta(\omega_P-\omega_s-\omega_i)\delta^{(3)}(\mathbf k_P-\mathbf k_s-\mathbf k_i)\hat a_s^\dagger(\mathbf k_s)\hat a_i^\dagger(\mathbf k_i)|0\rangle. $$

Audit: Correct as a distributional statement; the smeared form [a(f), a^dagger(g)] = <f, g> I is right. Checked by: hand derivation.
QF code (exact): none. continuum modes; QF works with finitely many modes.

## Non-formula records

- G-085: conditional measurement prepares cats; the only formula ($\propto\hat I\pm\hat\Pi$) is covered by the cat parity entry. Listed in "Coherent-state cat" Appears in.

## QF notes from this section

- No wrong value from an exact QF call. `BosonicNormalOrder`, `BosonicExpOrder` (one- and two-mode, normal and antinormal), `BosonicExpSimplify` (one-, two- and three-mode conjugations, products of displacements, conjugation of an exponential), `BosonicExpBraid`, `BosonicVEV`, `BosonicMatrixElement` (symbolic and integer Fock indices, one and several modes, coherent basis), and the symbolic `DisplacementOperator`/`SqueezeOperator` all agree with the atlas.
- Missing features met here, all returning the input unevaluated rather than a wrong value: `BosonicMatrixElement` with an exponential of `ad ** ad` in the coherent basis; `BosonicMatrixElement[{m, m}, DisplacementOperator[al]]` and `BosonicMatrixElement[{2 m, 0}, SqueezeOperator[xi]]` with symbolic m (integer m works); exponentials of `(ad ** a) ** (ad ** a)` (Kerr); `BosonicVEV` of a product of exponentials (route through `BosonicExpOrder` and `BosonicMatrixElement[{0, 0}, ...]`); `BosonicExpSimplify` of the Kerr conjugation, whose commutator series does not close.
- `BosonicExpSimplify[SqueezeOperator[xi, Infinity] ** SqueezeOperator[-xi, Infinity]]` returns `E^(X/2 + (-X)/2)` with the two halves not combined; expanding the exponent gives 1. Not a wrong value.
- `BosonicNormalOrder` needs the variable list as its second argument; with one argument it returns unevaluated.
