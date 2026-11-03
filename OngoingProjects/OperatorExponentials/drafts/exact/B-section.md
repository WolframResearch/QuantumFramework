# B: Berezin, The Method of Second Quantization (atlas lines 338 to 1397, IDs B-001 to B-108)

Every QF cell below runs after these loading lines:

```wl
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

Each cell is stored verbatim as `drafts/exact/B-work/cells/NN.wl` (NN is the entry's position) and was run by `drafts/exact/B-work/runcells.wls` in a fresh Global context; outputs are in `drafts/exact/B-work/cells-out-all.txt`. The exact computations outside QF (Gaussian integrals, exact Bogoliubov matrices, the Bargmann vacuum equation, the Riccati fixed point) are in `drafts/exact/B-work/checks.wl`, output `drafts/exact/B-work/checks-out.txt`. No truncated Fock matrix is used anywhere; fermionic checks use exact Jordan-Wigner matrices or an explicit Grassmann algebra.

## Errors found

1. **B-071, printed (6.23) lacks the 1/2.** Rechecked exactly: inserting the corrected c(t) into the vacuum column of e^{itH} in the Bargmann realization (a acts as d/dz, a* as z) leaves residual exactly 0 for one mode with symbolic t, and below 1e-57 for two modes at 60 digits; the printed coefficient leaves a nonzero residual.
2. **B-081, printed sign of the constant.** Rechecked with QF's exact normal ordering: H - (omega b*b - g^2|f|^2/omega) is exactly 0, while the printed + sign leaves 2 g^2 |f|^2/omega.
3. **B-099, wording that flips the fermionic sign.** The atlas says "the first derivative is the right derivative, the second the left derivative". Read by position this gives the wrong fermionic contraction sign and B-101, B-102 fail; read by variable (right derivative in a*, left derivative in a) the contraction is -Delta and B-101, B-102 hold exactly (explicit two-mode Grassmann algebra).
4. **Atlas corrections confirmed:** B-036 (normalized paired integral equals (-1)^{n(n-1)/2} Pf, explicit Grassmann algebra n = 1, 2, 3); B-075 (the integral needs Re z below the spectrum: WL's exact `Integrate` returns the condition Re[c n] > Re[z]); B-077 (the trace-class condition is on Phi e^{itC} - I; no formula to compute).
5. **B-045, the warning is now exercised.** For a two-mode transformation built from QF's exact conjugations, ||Z|| = ||W|| = 4/5 while ||Psi Phi^{-1}|| = 4/3, so the unbarred contraction does exceed 1, as the atlas says.

No QF exact call returned a wrong value. Two exact calls stop short: `BosonicExpOrder[Exp[-(w + mu) ad ** a]]` returns its input unchanged when the coefficient of a*a is a sum (the same call with a single symbol in place of the sum normal-orders correctly), and `BosonicExpSimplify` leaves e^{G} a e^{-G} unevaluated when G mixes quadratic and linear terms, and also e^{G} (a + f) e^{-G} for a pure squeeze G (the conjugation of a alone closes).

---

### Smeared CCR and CAR
Topic: Infinite modes and domains
Appears in: B-003

Smeared bosonic annihilators and creators commute to the bilinear pairing of their test functions, and smeared fermionic ones anticommute to the same pairing; like-type pairs commute (bosons) or anticommute (fermions) to zero.

$$ [\hat a(f),\hat a^*(g)]=(f,g^*)\hat I,\qquad \{\hat a(f),\hat a^*(g)\}=(f,g^*)\hat I $$

Conditions: bilinear pairing fg = (f, g*); like-type (anti)commutators vanish.
Audit: Correct. Checked by: QF exact (bosons, two modes) and exact fermion matrices (2 modes).
QF code (exact):
```wl
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
v = {a1, a1d, a2, a2d};
af = f1 a1 + f2 a2; adg = g1 a1d + g2 a2d;
BosonicNormalOrder[#, v] & /@ {af ** adg - adg ** af, a1 ** a2 - a2 ** a1, a1d ** a2d - a2d ** a1d}
```
Output: `{f1*g1 + f2*g2, 0, 0}`. Matches: the smeared commutator is f1 g1 + f2 g2 = (f, g*), the like-type commutators vanish.

### Fock vector from the vacuum
Topic: Fock matrix elements and generating functions
Appears in: B-004

Any Fock vector is a sum over particle number of creators applied to the vacuum, weighted by a symmetric (bosons) or antisymmetric (fermions) kernel; the 1/sqrt(n!) makes the squared norm the plain sum of integrals of |K_n|^2.

$$ \Phi=\sum_{n\ge0}\frac1{\sqrt{n!}}\int K_n(x_1,\ldots,x_n)\hat a^*(x_1)\cdots\hat a^*(x_n)\,d^nx\,\Omega $$

Audit: Correct. Checked by: QF exact (two-particle sector, two modes, symmetric K_2).
QF code (exact):
```wl
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
k = {{k11, k12}, {k12, k22}}; kb = {{kb11, kb12}, {kb12, kb22}};
ket = Sum[k[[i, j]] {a1d, a2d}[[i]] ** {a1d, a2d}[[j]], {i, 2}, {j, 2}]/Sqrt[2];
bra = Sum[kb[[i, j]] {a1, a2}[[j]] ** {a1, a2}[[i]], {i, 2}, {j, 2}]/Sqrt[2];
Expand[BosonicVEV[bra ** ket, {a1, a1d, a2, a2d}]]
```
Output: `k11*kb11 + 2*k12*kb12 + k22*kb22`. Matches: the squared norm of (1/sqrt 2!) Sum K_ij a*_i a*_j Omega is Sum_ij |K_ij|^2.

### Norm bound for fermionic fields
Topic: Fermionic canonical transformations
Appears in: B-005

A smeared fermionic field is a bounded operator whose norm is at most the norm of its test function; in fact the bound is attained.

$$ \|\hat a(f)\|,\ \|\hat a^*(f)\|\le\|f\| $$

Audit: Correct, with equality: a(f)^2 = 0 and {a(f), a(f)*} = ||f||^2 make a(f)* a(f) / ||f||^2 a projector. Checked by: hand derivation: (a(f)* a(f))^2 = a(f)* (||f||^2 - a(f)* a(f)) a(f) = ||f||^2 a(f)* a(f), and exact fermion matrices (2 modes, operator norm equal to ||f||).
QF code (exact): none. No fermionic layer.

### Bosonic quadratures
Topic: Infinite modes and domains
Appears in: B-006

The smeared position and momentum quadratures are the Hermitian combinations of a smeared field and its adjoint, normalized so that [Q, P] = i for a normalized real test function.

$$ \hat Q(f)=(\hat a(f)+\hat a^*(f^*))/\sqrt2,\qquad \hat P(f)=(\hat a(f)-\hat a^*(f^*))/(i\sqrt2) $$

Audit: Correct as a definition; the essential self-adjointness on finite-particle vectors is a source theorem, Unchecked (domain statement, no finite algebraic content). Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2]; p = (a - ad)/(I Sqrt[2]);
BosonicNormalOrder[q ** p - p ** q, {a, ad}]
```
Output: `I`. Matches ([Q, P] = i). QF's own `QuadratureOperators` uses (a + a*)/2, (a - a*)/(2i), a factor sqrt 2 different.

### Vacuum projector as a normal-ordered exponential
Topic: Normal symbols and normal-order exponentials
Appears in: B-013 (related: B-107 normal-order identity)

The projector onto the vacuum is the normal-ordered exponential of minus the number operator: the normal-ordered exponential :e^{c a*a}: multiplies |n> by (1 + c)^n, and c = -1 keeps only n = 0.

$$ \hat P_0=:\!e^{-\hat a^*\hat a}\!: $$

Audit: Correct. Checked by: QF exact (symbolic Fock index).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicMatrixElement[{m, m}, NormalOrdered[Exp[-ad ** a]]],
 Table[BosonicMatrixElement[{j, k}, NormalOrdered[Exp[-ad ** a]]], {j, 0, 3}, {k, 0, 3}]}
```
Output: `{KroneckerDelta[0, m], {{1, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}}`. Matches: the diagonal element is 1 at m = 0 and 0 otherwise, and every off-diagonal element vanishes.

### Matrix kernel versus normal symbol
Topic: Normal symbols and normal-order exponentials
Appears in: B-014

The matrix kernel of an operator (its matrix element between unnormalized coherent vectors) equals its normal symbol times the coherent overlap e^{a*a}.

$$ N_{\hat A}(a^*,a)=K_{\hat A}(a^*,a)e^{-a^*a},\qquad K_{\hat A}=N_{\hat A}e^{a^*a} $$

Audit: Correct. Checked by: QF exact (symbolic coherent amplitudes).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
toSym[e_] := e /. {ad -> zb, a -> z} /. GeneralizedPower[NonCommutativeMultiply, x_, n_] :> x^n /. NonCommutativeMultiply -> Times;
op = a ** ad ** a ** ad;
nsym = Expand[toSym[BosonicNormalOrder[op, {a, ad}]]];
kcoh = BosonicMatrixElement[{al, be}, op, "Basis" -> "Coherent"];
{nsym, Simplify[kcoh Exp[(Abs[al]^2 + Abs[be]^2)/2] - (nsym /. {zb -> Conjugate[al], z -> be}) Exp[Conjugate[al] be]]}
```
Output: `{1 + 3*z*zb + z^2*zb^2, 0}`. Matches: with a* = conj(al), a = be, the unnormalized coherent matrix element equals N times e^{conj(al) be} exactly.

### Symbols as coherent matrix elements
Topic: Normal symbols and normal-order exponentials
Appears in: B-019

The normal symbol on the diagonal a* = conj(a) is the expectation value of the operator in the normalized coherent state of amplitude a.

$$ N_{\hat A}(a^*,a)=\langle P_a,\hat A P_a\rangle $$

Audit: Correct. Checked by: QF exact (symbolic amplitude).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
toSym[e_] := e /. {ad -> zb, a -> z} /. GeneralizedPower[NonCommutativeMultiply, x_, n_] :> x^n /. NonCommutativeMultiply -> Times;
op = a ** ad ** a ** ad;
nsym = Expand[toSym[BosonicNormalOrder[op, {a, ad}]]];
diag = BosonicMatrixElement[{al, al}, op, "Basis" -> "Coherent"] /. Abs[al]^2 -> al Conjugate[al];
{diag, Simplify[diag - (nsym /. {zb -> Conjugate[al], z -> al})]}
```
Output: `{1 + 3*al*Conjugate[al] + al^2*Conjugate[al]^2, 0}`. Matches.

### Differential realization on normal symbols
Topic: Normal symbols and normal-order exponentials
Appears in: B-015

Multiplying an operator by a on the left, or by a* on the right, acts on its normal symbol as multiplication plus a derivative in the other variable, because only those two products need reordering.

$$ N_{\hat a\hat A}=(a+\partial_{a^*})N_{\hat A},\qquad N_{\hat A\hat a^*}=(a^*+\partial_{a})N_{\hat A} $$

Conditions: polynomial or analytic symbols, smearing implicit.
Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
toSym[e_] := e /. {ad -> zb, a -> z} /. GeneralizedPower[NonCommutativeMultiply, x_, n_] :> x^n /. NonCommutativeMultiply -> Times;
q = ad ** ad ** a + 2 ad + a ** a; s = toSym[q];
{Expand[toSym[BosonicNormalOrder[a ** q, {a, ad}]] - (z s + D[s, zb])],
 Expand[toSym[BosonicNormalOrder[q ** ad, {a, ad}]] - (zb s + D[s, z])]}
```
Output: `{0, 0}`. Matches (both differences vanish for q = a*^2 a + 2 a* + a^2).

### Fermionic symbols: which side the derivative acts on
Topic: Grassmann calculus and fermionic coherent states
Appears in: B-016

For fermions, multiplying by a on the right or by a* on the left keeps normal order, so the normal symbol is just multiplied on the same side by the Grassmann variable.

$$ N_{\hat A\hat a}=N_{\hat A}a,\qquad N_{\hat a^*\hat A}=a^*N_{\hat A} $$

Audit: Correct: right multiplication by a and left multiplication by a* keep normal order, so the symbol is multiplied on the same side; the derivative sides only enter for the other two products. Checked by: hand derivation: a* (a*^p a^q) and (a*^p a^q) a are already normal ordered.
QF code (exact): none. No fermionic layer.

### Poisson (coherent) vectors
Topic: Coherent states (general)
Appears in: B-017

The unnormalized coherent vector e^{a*f} has squared norm e^{|f|^2}; dividing by e^{|f|^2/2} gives the normalized coherent vector, whose photon-number distribution is Poisson with mean |f|^2.

$$ P_f(a^*)=e^{a^*f-\|f\|^2/2},\qquad Q_f(a^*)=e^{a^*f},\quad \|P_f\|=1,\quad\|Q_f\|^2=e^{\|f\|^2} $$

Audit: Correct, including the Poisson occupation law. Checked by: QF exact (symbolic Fock index, exact sum).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
qn = Refine[BosonicMatrixElement[{n, 0}, Exp[f ad]], n >= 0];
w = qn (qn /. f -> fb);
{DisplacementOperator[f, Infinity, "Ordering" -> "Normal"], qn, Sum[w, {n, 0, Infinity}], Exp[-f fb] w}
```
Output: `{E^(f*SuperDagger[\[FormalA]1])**E^(-(\[FormalA]1*Conjugate[f]))/E^((f*Conjugate[f])/2), f^n/Sqrt[n!], E^(f*fb), (f^n*fb^n)/(E^(f*fb)*n!)}`. Matches after writing fb for conj(f): D(f) = e^{-|f|^2/2} e^{f a*} e^{-conj(f) a} so P_f = e^{-|f|^2/2} Q_f; <n|Q_f> = f^n/sqrt(n!); ||Q_f||^2 = e^{|f|^2}; the occupation law is e^{-|f|^2}|f|^{2n}/n!, which sums to 1.

### Rank-one projector symbols
Topic: Normal symbols and normal-order exponentials
Appears in: B-020

The projector onto a normalized vector Phi has matrix kernel Phi(a*) Phi*(a) and normal symbol that times e^{-a*a}.

$$ K_{|\Phi\rangle\langle\Phi|}=\Phi(a^*)\Phi^*(a),\qquad N_{|\Phi\rangle\langle\Phi|}=\Phi(a^*)\Phi^*(a)e^{-a^*a} $$

Audit: Correct (bosons; fermions with the reversed factor order). Checked by: QF exact for Phi = Omega (Phi = 1, so N = e^{-a*a}) and hand derivation for Phi = |1>: <z|1><1|z> = |z|^2 e^{-|z|^2}, which is z* z e^{-z* z}.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
BosonicMatrixElement[{al, al}, NormalOrdered[Exp[-ad ** a]], "Basis" -> "Coherent"]
```
Output: `E^(-Abs[al]^2)`. Matches (Phi = Omega). For Phi = |1>, `BosonicMatrixElement[{al, al}, ad ** NormalOrdered[Exp[-ad ** a]] ** a, "Basis" -> "Coherent"]` returns unevaluated.

### Real Gaussian integral
Topic: Gaussian integrals and Pfaffians
Appears in: B-022

The normalized real Gaussian integral of e^{-(Af, f)} is det A^{-1/2}.

$$ \int e^{-(\hat Af,f)}\Pi df=(\det \hat A)^{-1/2},\qquad \hat A>0 $$

Conditions: normalized measure pi^{-n/2} d^n f; infinite dimension needs A - I trace class.
Audit: Correct. Checked by: exact WL computation (2x2 positive rational A, exact `Integrate`, difference 0).
QF code (exact): none. A c-number Gaussian integral; no second-quantization call is involved.

### Doubled complex Gaussian with sources
Topic: Gaussian integrals and Pfaffians
Appears in: B-023

A complex Gaussian with both holomorphic and antiholomorphic quadratic terms and linear sources integrates to the inverse square root of a row-swapped block determinant times the exponential of half the inverse form on the sources.

$$ \int e^{-\frac12(f\ f^*)\widehat{\mathcal A}\binom f{f^*}+(\varphi_1\ \varphi_2)\binom f{f^*}}\Pi df\,df^* =\det\!\begin{pmatrix}\hat A_{21}&\hat A_{22}\\\hat A_{11}&\hat A_{12}\end{pmatrix}^{-1/2} e^{\frac12(\varphi_1\ \varphi_2)\widehat{\mathcal A}^{-1}\binom{\varphi_1}{\varphi_2}} $$

Conditions: A transpose-symmetric, A12 - I and A21 - I trace class, A11, A22 Hilbert-Schmidt, real part of the form bounded below.
Audit: Correct, with the row-swapped determinant as printed (for one mode it equals -det A). Checked by: exact WL computation (one mode, complex A11, A22, A12 = A21 = 1 and complex sources, measure d^2f/pi: the exact `Integrate` and the formula agree; `Simplify` stalls on the nested radicals, and the difference evaluates to 0 at 99 digits).
QF code (exact): none. A c-number Gaussian integral; no second-quantization call is involved.

### Kernel acting on a vector
Topic: Normal symbols and normal-order exponentials
Appears in: B-026

An operator acts on a vector functional by integrating its matrix kernel against the vector with the Bargmann weight e^{-alpha* alpha}.

$$ (\hat A\Phi)(a^*)=\int K_{\hat A}(a^*,\alpha)\Phi(\alpha^*)e^{-\alpha^*\alpha}\Pi d\alpha^*d\alpha $$

Audit: Correct. Checked by: hand derivation: insert the coherent resolution of identity int |alpha><alpha| e^{-alpha* alpha} d^2alpha/pi = I between A and Phi.
QF code (exact): none. QF has no Bargmann-integral representation of operator action; `BosonicMatrixElement` gives matrix elements only.

### Product of operators: matrix kernels
Topic: Normal symbols and normal-order exponentials
Appears in: B-027

The matrix kernel of a product is the Bargmann-weighted integral of the two kernels over the intermediate coherent label.

$$ K_{\hat A\hat B}(a^*,a)=\int K_{\hat A}(a^*,\alpha)K_{\hat B}(\alpha^*,a)e^{-\alpha^*\alpha}\Pi d\alpha^*d\alpha $$

Audit: Correct. Checked by: hand derivation: the same resolution of identity as B-026, inserted between A and B.
QF code (exact): none. QF has no kernel-composition integral; the normal-symbol version is checked in the next entry.

### Product of operators: normal symbols
Topic: Normal symbols and normal-order exponentials
Appears in: B-027

The normal symbol of a product is a Gaussian convolution of the two normal symbols; for polynomials it equals the series N_A exp(left d/da, right d/da*) N_B.

$$ N_{\hat A\hat B}(a^*,a)=\int N_{\hat A}(a^*,\alpha)N_{\hat B}(\alpha^*,a)e^{-(\alpha^*-a^*)(\alpha-a)}\Pi d\alpha^*d\alpha $$

Audit: Correct. Checked by: QF exact (the series form) and exact WL computation (the Gaussian integral at a = 3/10 + i/5, exact `Integrate`, difference 0).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
toSym[e_] := e /. {ad -> zb, a -> z} /. GeneralizedPower[NonCommutativeMultiply, x_, n_] :> x^n /. NonCommutativeMultiply -> Times;
p1 = ad ** a ** a; p2 = ad ** ad ** a + a;
{s1, s2} = toSym /@ {p1, p2};
s12 = Expand[toSym[BosonicNormalOrder[p1 ** p2, {a, ad}]]];
{s12, Expand[s12 - Sum[D[s1, {z, k}] D[s2, {zb, k}]/k!, {k, 0, 4}]]}
```
Output: `{2*z*zb + z^3*zb + 4*z^2*zb^2 + z^3*zb^3, 0}`. Matches: the normal symbol of the product equals the derivative series exactly.

### Bosonic trace
Topic: Thermal states and traces
Appears in: B-028

The trace of a trace-class operator is the integral of its normal symbol over the normalized complex plane (equivalently, of its kernel times e^{-a*a}).

$$ \operatorname{Tr}\hat A=\int K_{\hat A}(a^*,a)e^{-a^*a}\Pi da^*da=\int N_{\hat A}(a^*,a)\Pi da^*da $$

Conditions: A trace class.
Audit: Correct. Checked by: QF exact plus exact WL computation, for A = e^{-bN}.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
tr = Integrate[Exp[(E^-b - 1) (x^2 + y^2)]/Pi, {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, Assumptions -> b > 0];
{BosonicExpOrder[Exp[-b ad ** a]], tr, Simplify[tr - Sum[E^(-b n), {n, 0, Infinity}], b > 0]}
```
Output: `{NormalOrdered[E^((-1 + E^(-b))*SuperDagger[\[FormalA]1]**\[FormalA]1)], E^b/(-1 + E^b), 0}`. Matches: the normal symbol e^{(e^{-b} - 1)|z|^2} integrates to 1/(1 - e^{-b}), which equals Sum_n e^{-bn}.

### Grassmann algebra relations
Topic: Grassmann calculus and fermionic coherent states
Appears in: B-029

Grassmann generators anticommute and square to zero, so n of them span an algebra of dimension 2^n.

$$ x_i x_j=-x_jx_i,\qquad x_i^2=0,\qquad \dim=2^n $$

Audit: Correct. Checked by: hand derivation: a basis is the 2^n ordered monomials over subsets of the generators.
QF code (exact): none. No fermionic layer.

### Left versus right derivatives
Topic: Grassmann calculus and fermionic coherent states
Appears in: B-030

On a monomial of degree p, the right and left Grassmann derivatives differ by the sign (-1)^{p-1}.

$$ \partial^{R}_{x}m=(-1)^{p-1}\,\partial^{L}_{x}m\quad\text{on a monomial }m\text{ of degree }p $$

Audit: Correct. Checked by: hand derivation: moving the target generator to the right end instead of the left end takes p - 1 transpositions more in total.
QF code (exact): none. No fermionic layer.

### Berezin integration rules
Topic: Grassmann calculus and fermionic coherent states
Appears in: B-031

Berezin integration kills constants and gives 1 on the integration variable itself.

$$ \int dx_i=0,\qquad \int x_i\,dx_i=1 $$

Audit: Correct as the defining rules. The partial-integration warning is a source statement, Unchecked. Checked by: hand derivation (definition).
QF code (exact): none. No fermionic layer.

### Grassmann change of variables
Topic: Grassmann calculus and fermionic coherent states
Appears in: B-032

Under a linear change of Grassmann variables the integral picks up the inverse determinant, the opposite of commuting variables.

$$ \int f(x)\,dx_n\cdots dx_1=(\det A)^{-1}\int f(Ay)\,dy_n\cdots dy_1 $$

Audit: Correct. Checked by: hand derivation: for f = x_1...x_n, f(Ay) = det A y_1...y_n, so the right side is (det A)^{-1} det A = 1, the left side.
QF code (exact): none. No fermionic layer.

### Fermionic Gaussian is a Pfaffian
Topic: Gaussian integrals and Pfaffians
Appears in: B-033

The Berezin integral of a Grassmann Gaussian with skew matrix M is the Pfaffian of M.

$$ \int e^{\frac12x^TMx}\,dx_n\cdots dx_1=\operatorname{Pf}(M),\quad M^T=-M $$

Audit: Correct. Checked by: exact WL computation (explicit Grassmann algebra, 4 generators, complex skew M: top coefficient minus Pf is 0).
QF code (exact): none. No fermionic layer.

### Clifford operators from multiplication and derivative
Topic: Grassmann calculus and fermionic coherent states
Appears in: B-034

Grassmann multiplication and left differentiation satisfy the CAR, so their sum and difference generate a Clifford algebra.

$$ \{\hat Q_i,\hat Q_j\}=\{\hat P_i,\hat P_j\}=2\delta_{ij}\hat I,\quad \{\hat P_i,\hat Q_j\}=0,\qquad \hat Q_j=\hat X_j+\hat D_j,\ \hat P_j=(\hat X_j-\hat D_j)/i $$

Audit: Correct. Checked by: exact fermion matrices (2 modes; X = c^dag, D = c; all residuals 0).
QF code (exact): none. No fermionic layer.

### Doubled fermionic Gaussian (corrected determinant)
Topic: Gaussian integrals and Pfaffians
Appears in: B-036

The paired Grassmann Gaussian integral equals the square root of the unpermuted block determinant, fixed in sign as (-1)^{n(n-1)/2} times the Pfaffian.

$$ \int e^{\frac12(\alpha\ \alpha^*)\widehat{\mathcal A}\binom\alpha{\alpha^*}}\Pi d\alpha^*d\alpha =\det\!\begin{pmatrix}\hat A_{11}&\hat A_{12}\\\hat A_{21}&\hat A_{22}\end{pmatrix}^{1/2} = (-1)^{n(n-1)/2}\operatorname{Pf}(\mathcal A) $$

Conditions: normalized paired integration, int e^{-alpha* alpha} = 1; block order (alpha_1..alpha_n, alpha*_1..alpha*_n).
Audit: Correct as corrected in the atlas. Checked by: exact WL computation (explicit Grassmann algebra, n = 1, 2, 3, complex skew matrices: value matches (-1)^{n(n-1)/2} Pf and its square matches det; the plain Pf fails for n = 2, 3).
QF code (exact): none. No fermionic layer.

### Fermionic trace
Topic: Thermal states and traces
Appears in: B-040

The fermionic trace is a Berezin integral of the kernel with weight e^{-aa*}, or of the normal symbol with the doubled weight e^{-2aa*}.

$$ \operatorname{Tr}\hat A=\int K_{\hat A}(a^*,a)e^{-aa^*}\Pi da\,da^* =\int N_{\hat A}(a^*,a)e^{-2aa^*}\Pi da\,da^* $$

Audit: Correct. Checked by: hand derivation (one mode): N = 1 gives int (1 - 2 a a*) = 2 = Tr I; N = a*a gives 1 = Tr N; odd symbols give 0; the formula is multiplicative over modes.
QF code (exact): none. No fermionic layer.

### Normalized fermionic trace
Topic: Thermal states and traces
Appears in: B-041

Rescaling the symbol's arguments by 1/sqrt 2 and using the single weight gives the trace divided by 2^n, the normalized trace that is 1 on the identity.

$$ \tau(\hat A)=\int N_{\hat A}(a^*/\sqrt2,a/\sqrt2)e^{-aa^*}\Pi da\,da^*=2^{-n}\operatorname{Tr}\hat A $$

Audit: Correct. Checked by: hand derivation (one mode): N = 1 gives 1, N = a*a gives 1/2.
QF code (exact): none. No fermionic layer.

### Bosonic affine canonical transformation
Topic: Bosonic canonical transformations
Appears in: B-042

A bosonic canonical transformation sends a to a linear combination of a and a* plus a constant shift; the transformed creator is the adjoint of that expression.

$$ \hat b=\hat\Phi\hat a+\hat\Psi\hat a^*+f\hat I,\qquad \hat b^*=\bar{\hat\Psi}\hat a+\bar{\hat\Phi}\hat a^*+f^*\hat I $$

Audit: Correct as a definition (the second line is the adjoint of the first in the bilinear convention). Checked by: QF exact (one mode: the displacement gives the shift, the squeezer gives Phi and Psi, and the creator transforms with the conjugate coefficients).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
g = (xb a ** a - x ad ** ad)/2;
{BosonicExpSimplify[Exp[fb a - f ad] ** a ** Exp[-(fb a - f ad)], {a, ad}],
 BosonicExpSimplify[Exp[g] ** a ** Exp[-g], {a, ad}],
 BosonicExpSimplify[Exp[g] ** ad ** Exp[-g], {a, ad}]}
```
Output: `{\[FormalA]1 + f, \[FormalA]1*Cosh[Sqrt[x*xb]] + (x*Sinh[Sqrt[x*xb]]*SuperDagger[\[FormalA]1])/Sqrt[x*xb], (\[FormalA]1*xb*Sinh[Sqrt[x*xb]])/Sqrt[x*xb] + Cosh[Sqrt[x*xb]]*SuperDagger[\[FormalA]1]}`. Matches: with U = S U_f, b = U a U^{-1} = cosh|x| a + (x/|x|) sinh|x| a* + f, and b* has the conjugate coefficients. The combined conjugation `BosonicExpSimplify[Exp[g] ** (a + f) ** Exp[-g], {a, ad}]` returns unevaluated, so linearity is applied by hand.

### Bosonic block relations and inverse
Topic: Bosonic canonical transformations
Appears in: B-043

The coefficient blocks of a bosonic canonical transformation satisfy Phi Phi* - Psi Psi* = I and Phi Psi' = Psi Phi', which make the doubled matrix invertible with the displayed inverse.

$$ \hat S=\begin{pmatrix}\hat\Phi&\hat\Psi\\\bar{\hat\Psi}&\bar{\hat\Phi}\end{pmatrix},\quad \hat S^{-1}=\begin{pmatrix}\hat\Phi^*&-\hat\Psi'\\-\hat\Psi^*&\hat\Phi'\end{pmatrix},\quad \hat\Phi\hat\Phi^*-\hat\Psi\hat\Psi^*=\hat I,\ \hat\Phi\hat\Psi'=\hat\Psi\hat\Phi' $$

Audit: Correct, including the inverse relations Phi* Phi - Psi' conj(Psi) = I and Phi* Psi = Psi' conj(Phi). Checked by: QF exact (two modes: doubled matrices of a squeezer on mode 1 and a beam splitter read off QF's exact conjugations, multiplied, then specialized to r = log 3, theta = pi/4).
QF code (exact):
```wl
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
v = {a1, a1d, a2, a2d};
dbl[g_] := Table[Coefficient[Expand[BosonicExpSimplify[Exp[g] ** x ** Exp[-g], v, "Scalars" -> {r, th}]], y], {x, {a1, a2, a1d, a2d}}, {y, {a1, a2, a1d, a2d}}];
s = FullSimplify[dbl[r (a1 ** a1 - a1d ** a1d)/2] . dbl[I th (a1d ** a2 + a2d ** a1)] /. {r -> Log[3], th -> Pi/4}];
{phi, psi} = {s[[1 ;; 2, 1 ;; 2]], s[[1 ;; 2, 3 ;; 4]]};
sinv = ArrayFlatten[{{ConjugateTranspose[phi], -Transpose[psi]}, {-ConjugateTranspose[psi], Transpose[phi]}}];
{s, Simplify[{s[[3 ;; 4]] - ArrayFlatten[{{Conjugate[psi], Conjugate[phi]}}],
  phi . ConjugateTranspose[phi] - psi . ConjugateTranspose[psi] - IdentityMatrix[2],
  phi . Transpose[psi] - psi . Transpose[phi], s . sinv - IdentityMatrix[4],
  ConjugateTranspose[phi] . phi - Transpose[psi] . Conjugate[psi] - IdentityMatrix[2],
  ConjugateTranspose[phi] . psi - Transpose[psi] . Conjugate[phi]}]}
```
Output: `{{{5/(3*Sqrt[2]), ((-5*I)/3)/Sqrt[2], (2*Sqrt[2])/3, ((2*I)/3)*Sqrt[2]}, {(-I)/Sqrt[2], 1/Sqrt[2], 0, 0}, <the two conjugate rows>}, {<six zero arrays>}}`. Matches: the lower blocks are the conjugates of the upper ones, and all five relations and S.S^{-1} = I hold exactly.

### Symplectic form and inverse in quadratures
Topic: Bosonic canonical transformations
Appears in: B-044

In real quadrature coordinates the transformation matrix is symplectic, and its inverse is obtained by transposing and rearranging the blocks.

$$ M^TJM=J,\qquad M^{-1}=\begin{pmatrix}D^T&-B^T\\-C^T&A^T\end{pmatrix} $$

Audit: Correct. Checked by: hand derivation: M^{-1} = -J M^T J, which multiplies out to the displayed blocks for J = ((0, I), (-I, 0)).
QF code (exact): none. A statement about real 2N x 2N matrices; no second-quantization call is involved.

### Squeezing contraction, both charts
Topic: Bosonic canonical transformations
Appears in: B-045

The contractions Z = Phi^{-1} Psi and W = conj(Psi) Phi^{-1} that parametrize the squeezed vacuum have norm below 1, with Z symmetric; the unbarred Psi Phi^{-1} need not.

$$ \hat Z=\hat\Phi^{-1}\hat\Psi,\ \hat Z\hat Z^*=\hat I-(\hat\Phi^*\hat\Phi)^{-1};\qquad \hat W=\bar{\hat\Psi}\hat\Phi^{-1},\ \hat W^*\hat W=\hat I-(\hat\Phi\hat\Phi^*)^{-1} $$

Audit: Correct; Z symmetric and both norms below 1. The atlas warning that ||Psi Phi^{-1}|| can exceed 1 is confirmed by the example below. Checked by: QF exact (same two-mode transformation as B-043) plus hand derivation: Psi Psi* = Phi Phi* - I gives Z Z* = I - Phi^{-1} Phi*^{-1}.
QF code (exact):
```wl
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
v = {a1, a1d, a2, a2d};
dbl[g_] := Table[Coefficient[Expand[BosonicExpSimplify[Exp[g] ** x ** Exp[-g], v, "Scalars" -> {r, th}]], y], {x, {a1, a2, a1d, a2d}}, {y, {a1, a2, a1d, a2d}}];
s = FullSimplify[dbl[r (a1 ** a1 - a1d ** a1d)/2] . dbl[I th (a1d ** a2 + a2d ** a1)] /. {r -> Log[3], th -> Pi/4}];
{phi, psi} = {s[[1 ;; 2, 1 ;; 2]], s[[1 ;; 2, 3 ;; 4]]};
z = Inverse[phi] . psi; w = Conjugate[psi] . Inverse[phi];
Simplify[{z - Transpose[z], z . ConjugateTranspose[z] - IdentityMatrix[2] + Inverse[ConjugateTranspose[phi] . phi],
  ConjugateTranspose[w] . w - IdentityMatrix[2] + Inverse[phi . ConjugateTranspose[phi]],
  Norm[z], Norm[w], Norm[psi . Inverse[phi]]}]
```
Output: `{{{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, 4/5, 4/5, 4/3}`. Matches: Z is symmetric, both identities hold exactly, ||Z|| = ||W|| = 4/5 < 1, while ||Psi Phi^{-1}|| = 4/3 > 1.

### Transformed vacuum (homogeneous part)
Topic: Gaussian vacua and overlaps
Appears in: B-048

The vacuum of the transformed annihilators is a Gaussian in a* with pairing matrix Phi^{-1} Psi, normalized by det(Phi Phi*)^{-1/4} up to a phase.

$$ F_0(a^*)=c\exp\!\left[-\tfrac12a^*\hat\Phi^{-1}\hat\Psi a^*-a^*\hat\Phi^{-1}f\right],\qquad c=\theta\det(\hat\Phi\hat\Phi^*)^{-1/4}\ (f=0) $$

Audit: Correct: (Phi d/dz + Psi z + f) F = 0 forces the exponent, and ||F_0|| = 1 with the Gaussian integral gives |c|^2 det(I - Z Z*)^{-1/2} = 1, i.e. |c| = det(Phi* Phi)^{-1/4} by B-045. Checked by: QF exact (two-mode squeezer, f = 0) and exact WL computation (two modes, exact non-commuting Bogoliubov blocks with f: the annihilation equation holds identically).
QF code (exact):
```wl
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
v = {a1, a1d, a2, a2d};
g = r (eb a1 ** a2 - e a1d ** a2d);
uN = BosonicExpOrder[Exp[g]];
b = Expand[BosonicExpSimplify[Exp[g] ** # ** Exp[-g], v, "Scalars" -> {r, e, eb}]] & /@ {a1, a2};
{phi, psi} = Table[Coefficient[bi, y], {bi, b}, {y, #}] & /@ {{a1, a2}, {a1d, a2d}} /. eb -> 1/e;
vac = uN /. {Exp[x_] /; ! FreeQ[x, a1] && FreeQ[x, a1d] :> 1, Power[s_, n_] /; s =!= E && ! FreeQ[n, a1d] :> 1} /. {a1d -> z1, a2d -> z2} /. NonCommutativeMultiply -> Times;
zz = Inverse[phi] . psi;
{uN, FullSimplify[(vac /. eb -> 1/e) - Det[phi . ConjugateTranspose[phi]]^(-1/4) Exp[-{z1, z2} . zz . {z1, z2}/2], r > 0]}
```
Output: `{E^(-((e*r*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]2]*Tan[Sqrt[-(e*eb*r^2)]])/Sqrt[-(e*eb*r^2)]))**Sec[Sqrt[-(e*eb*r^2)]]^(SuperDagger[\[FormalA]1]**\[FormalA]1 + SuperDagger[\[FormalA]2]**\[FormalA]2)**E^(<conjugate pair term>)*Sec[Sqrt[-(e*eb*r^2)]], 0}`. Matches: the vacuum column of QF's normal-ordered two-mode squeezer equals det(Phi Phi*)^{-1/4} exp(-a* Z a*/2) exactly (with the unit phase eb = 1/e).

### Transformed vacuum: affine normalization (4.25)
Topic: Gaussian vacua and overlaps
Appears in: B-048, B-049 ((4.26'))

With a shift f, the vacuum amplitude acquires a real Gaussian factor in f, so its modulus is det(Phi Phi*)^{-1/4} times the displayed exponential.

$$ c=\frac{\theta}{\det(\hat\Phi\hat\Phi^*)^{1/4}} \exp\!\left[\frac14(f\ f^*) \begin{pmatrix}\hat\Phi'^{-1}\hat\Psi^*&-\hat I\\-\hat I&\hat\Phi^{*-1}\hat\Psi'\end{pmatrix} \binom f{f^*}\right] $$

Audit: Correct (the exponent is real, as it must be for a modulus). Checked by: QF exact (one mode, U = S(xi) U_f, symbolic r, phi, f) and exact WL computation (two modes, exact Bogoliubov blocks and complex f: the Gaussian norm of F_0 gives the same |c|, `FullSimplify` of the difference is 0).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{sN, fN} = {SqueezeOperator[xi, Infinity, "Ordering" -> "Normal"], DisplacementOperator[-f, Infinity, "Ordering" -> "Normal"]};
mid = BosonicExpOrder[BosonicExpBraid[Exp[E^(-I phi) Tanh[r] a ** a/2] ** Exp[-f ad]]];
c = Sech[r]^(1/2) E^(-f fb/2) BosonicMatrixElement[{0, 0}, Evaluate[mid]];
{ph, ps} = {Cosh[r], E^(I phi) Sinh[r]};
formula = (ph ph)^(-1/4) Exp[(f^2 Conjugate[ps]/ph - 2 f fb + fb^2 ps/ph)/4];
cplx = {f -> fx + I fy, fb -> fx - I fy};
{sN, fN, mid, c, FullSimplify[ComplexExpand[Abs[c /. cplx]] - ComplexExpand[formula /. cplx], r > 0]}
```
Output: `{<S(xi) normal form, as in B-051>, E^(-(f*SuperDagger[\[FormalA]1]))**E^(\[FormalA]1*Conjugate[f])/E^((f*Conjugate[f])/2), E^((f^2*Tanh[r])/(2*E^(I*phi)))*E^(-(f*SuperDagger[\[FormalA]1]))**E^(((-2*\[FormalA]1*f + \[FormalA]1**\[FormalA]1)*Tanh[r])/(2*E^(I*phi))), E^(-1/2*(f*fb) + (f^2*Tanh[r])/(2*E^(I*phi)))*Sqrt[Sech[r]], 0}`. Matches after writing f = fx + i fy: |<0|S U_f|0>| equals the formula exactly. The middle step uses `BosonicExpBraid` and `BosonicExpOrder` because QF cannot order the whole product S ** U_f at once.

### Full bosonic Gaussian kernel
Topic: Bosonic canonical transformations
Appears in: B-049 (affine form also used by B-071)

The matrix kernel of the implementer is the vacuum normalization c times a Gaussian whose quadratic part is built from Phi^{-1} and Psi and whose linear part carries the shift.

$$ K_{\hat U}(a^*,a)=c\exp\!\left[ \tfrac12(a\ a^*) \begin{pmatrix}\bar{\hat\Psi}\hat\Phi^{-1}&\hat\Phi'^{-1}\\\hat\Phi^{-1}&-\hat\Phi^{-1}\hat\Psi\end{pmatrix} \binom a{a^*} +a(f^*-\bar{\hat\Psi}\hat\Phi^{-1}f)-a^*\hat\Phi^{-1}f\right], \quad N_{\hat U}=K_{\hat U}e^{-a^*a} $$

Audit: Correct. Checked by: QF exact (one mode, affine U = S(xi) U_f brought to normal order with QF's braid, ordering and conjugation rules; symbolic r, phi, f) and exact WL computation (two modes, exact non-commuting Bogoliubov blocks and complex f: K satisfies b U = U a and b* U = U a* as differential equations, residual 0).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
tau = E^(I phi) Tanh[r]; taub = E^(-I phi) Tanh[r];
mid = BosonicExpOrder[BosonicExpBraid[Exp[taub a ** a/2] ** Exp[-f ad]]];
mv = FullSimplify[BosonicExpSimplify[Exp[l ad ** a] ** Exp[-f ad] ** Exp[-l ad ** a], {a, ad}] /. l -> -Log[Cosh[r]], r > 0];
u = Sech[r]^(1/2) E^(-f fb/2) E^(f^2 taub/2) Exp[-tau ad ** ad/2] ** mv ** Sech[r]^(ad ** a) ** Exp[taub (a ** a - 2 f a)/2] ** Exp[fb a];
kern = u /. Sech[r]^(ad ** a) -> Exp[Sech[r] zs z] /. {ad -> zs, a -> z} /. NonCommutativeMultiply -> Times;
{ph, ps} = {Cosh[r], E^(I phi) Sinh[r]};
expo = (Conjugate[ps]/ph z^2 + 2 z zs/ph - ps/ph zs^2)/2 + z (fb - Conjugate[ps] f/ph) - zs f/ph;
{mid, mv, FullSimplify[PowerExpand[Log[kern/(kern /. {z -> 0, zs -> 0})]] - expo, r > 0 && Element[phi, Reals]]}
```
Output: `{E^((f^2*Tanh[r])/(2*E^(I*phi)))*E^(-(f*SuperDagger[\[FormalA]1]))**E^(((-2*\[FormalA]1*f + \[FormalA]1**\[FormalA]1)*Tanh[r])/(2*E^(I*phi))), E^(-(f*Sech[r]*SuperDagger[\[FormalA]1])), 0}`. Matches: log(K/c) minus the displayed exponent is exactly 0. The normal-ordered product u is assembled from the QF outputs `mid` (the squeeze factor moved past the displacement) and `mv` (Sech[r]^N moved past e^{-f a*}).

### Displacement in the source direction
Topic: Displacement
Appears in: B-050

Berezin's displacement U_f shifts a by +f under U a U^{-1}; it is D(-f) in the common convention, and its normal symbol is a Gaussian-normalized exponential of linear terms.

$$ \hat U_f=e^{\hat a(f^*)-\hat a^*(f)},\quad \hat U_f\hat a\hat U_f^{-1}=\hat a+f\hat I, \qquad N_{\hat U_f}=e^{-\|f\|^2/2+af^*-a^*f} $$

Audit: Convention: U_f = D(-f) = D(f)^dag in the common convention D(f) = exp(f a^dag - f* a). Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{DisplacementOperator[-f, Infinity, "Ordering" -> "Normal"],
 BosonicExpSimplify[Exp[fb a - f ad] ** a ** Exp[-(fb a - f ad)], {a, ad}]}
```
Output: `{E^(-(f*SuperDagger[\[FormalA]1]))**E^(\[FormalA]1*Conjugate[f])/E^((f*Conjugate[f])/2), \[FormalA]1 + f}`. Matches: the normal-ordered form has symbol e^{-|f|^2/2} e^{-f a*} e^{conj(f) a}, and the conjugation gives a + f.

### One-mode squeezing implementer and kernel
Topic: Single-mode squeezing
Appears in: B-051

The one-mode squeeze exp[(sigma/2)(conj(eps) a^2 - eps a*^2)] maps a to a cosh sigma + eps a* sinh sigma, and its kernel is a Gaussian with coefficients tanh sigma and sech sigma.

$$ \hat U_\sigma=\exp\!\left[\tfrac\sigma2(\bar\epsilon\hat a^2-\epsilon\hat a^{*2})\right], \quad K_{\hat U_\sigma}=(\cosh\sigma)^{-1/2} e^{\frac12\bar\epsilon\tanh\sigma\,a^2+\operatorname{sech}\sigma\,a^*a-\frac12\epsilon\tanh\sigma\,a^{*2}} $$

Conditions: b = a cosh sigma + epsilon a* sinh sigma, |epsilon| = 1, sigma real.
Audit: Correct; U_sigma is QF's S(xi) with xi = sigma epsilon. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
kern = Sech[r]^(1/2) Exp[E^(-I phi) Tanh[r] z^2/2 + Sech[r] zs z - E^(I phi) Tanh[r] zs^2/2];
{SqueezeOperator[xi, Infinity, "Ordering" -> "Normal"],
 BosonicExpSimplify[Exp[(xb a ** a - x ad ** ad)/2] ** a ** Exp[-(xb a ** a - x ad ** ad)/2], {a, ad}],
 FullSimplify[Table[BosonicMatrixElement[{m, n}, SqueezeOperator[r E^(I phi)]] - Sqrt[m! n!] SeriesCoefficient[kern, {zs, 0, m}, {z, 0, n}], {m, 0, 4}, {n, 0, 4}], r > 0 && -Pi < phi < Pi]}
```
Output: `{E^(-1/2*(E^(I*Arg[xi])*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]1]*Tanh[Abs[xi]]))**Sech[Abs[xi]]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^((\[FormalA]1**\[FormalA]1*Tanh[Abs[xi]])/(2*E^(I*Arg[xi]))), \[FormalA]1*Cosh[Sqrt[x*xb]] + (x*Sinh[Sqrt[x*xb]]*SuperDagger[\[FormalA]1])/Sqrt[x*xb], <5 x 5 array of 0>}`. Matches: the conjugation gives b with x = sigma epsilon, and every coefficient of the kernel (m, n up to 4) equals QF's closed-form <m|S(r e^{i phi})|n>.

### Passive transformation of normal symbols
Topic: SU(2) and passive mode mixing
Appears in: B-053

A passive (number-conserving) transformation acts on normal symbols by substituting conj(Phi) a* and Phi a, because it never mixes a with a* and so keeps normal order.

$$ N_{\hat U\hat A\hat U^{-1}}(a^*,a)=N_{\hat A}(\bar{\hat\Phi} a^*,\hat\Phi a) $$

Conditions: Psi = f = 0, Phi unitary.
Audit: Correct. Checked by: QF exact (two-mode beam splitter with symbolic angle and phase, A = a1*^2 a2 + 3 a2* a1 a2).
QF code (exact):
```wl
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
v = {a1, a1d, a2, a2d};
toSym[e_] := e /. {a1d -> zs1, a2d -> zs2, a1 -> z1, a2 -> z2} /. GeneralizedPower[NonCommutativeMultiply, x_, n_] :> x^n /. NonCommutativeMultiply -> Times;
g = th (E^(I chi) a1d ** a2 - E^(-I chi) a2d ** a1);
conj[x_] := BosonicExpSimplify[Exp[g] ** x ** Exp[-g], v];
phi = Table[Coefficient[Expand[conj[x]], y], {x, {a1, a2}}, {y, {a1, a2}}];
nA[zs_, z_] := zs[[1]]^2 z[[2]] + 3 zs[[2]] z[[1]] z[[2]];
lhs = toSym[BosonicNormalOrder[conj[a1d ** a1d ** a2] + 3 conj[a2d ** a1 ** a2], v]];
{FullSimplify[phi, Element[th, Reals]], FullSimplify[lhs - nA[Conjugate[phi] . {zs1, zs2}, phi . {z1, z2}], Element[th | chi, Reals]]}
```
Output: `{{{Cos[th], -(E^(I*chi)*Sin[th])}, {Sin[th]/E^(I*chi), Cos[th]}}, 0}`. Matches: Phi is the displayed unitary, and the normal symbol of U A U^{-1} equals N_A(conj(Phi) a*, Phi a) exactly. A sum inside one conjugation is not distributed by `BosonicExpSimplify`, so each monomial is conjugated separately.

### General bosonic transformation kernel for normal symbols
Topic: Bosonic canonical transformations
Appears in: B-054

For a mixing transformation, the conjugated normal symbol is a limit of Gaussian integrals of the original symbol against a kernel built from the displayed matrix L.

$$ N_{\hat U\hat A\hat U^{-1}}=\lim_n\int\mathcal K_n(a^*,a\mid\alpha^*,\alpha)N_{\hat A}(\alpha^*,\alpha)\,d^nt\,d^n\tau,\qquad \widehat{\mathcal L}=\begin{pmatrix}\bar{\hat\Phi}\bar{\hat\Psi}^{-1}&-\hat I\\-\hat I&\hat\Phi\hat\Psi^{-1}\end{pmatrix} $$

Audit: Unchecked: the Gaussian kernel (4.39) and its contour are not transcribed in the atlas, so there is nothing complete to test. Checked by: prior audit only.
QF code (exact): none. The kernel (4.39) is not transcribed, so there is no formula to compare a QF result against.

### Gaussian normal symbols stay Gaussian under conjugation
Topic: Bosonic canonical transformations
Appears in: B-055

Conjugating an operator whose normal symbol is Gaussian gives another Gaussian normal symbol, with prefactor det C^{-1/2} for the displayed block matrix C.

$$ \hat C=\begin{pmatrix}\bar{\hat\Phi}+\hat A_{11}\hat\Psi&(\hat A_{12}-\hat I)\bar{\hat\Psi}\\(\hat A_{21}-\hat I)\hat\Psi&\hat\Phi+\hat A_{22}\bar{\hat\Psi}\end{pmatrix},\qquad \text{prefactor }\det\hat C^{-1/2} $$

Audit: Correct for the part that can be tested (the full exponent is not transcribed). For N_A = e^{-lambda a* a}, i.e. A = (1 - lambda)^N, and U the squeeze of B-051, the prefactor is the vacuum value of the conjugated symbol, <0|U A U^{-1}|0> = det C^{-1/2}. Checked by: QF exact (symbolic r, phi, lambda).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
nA = BosonicExpOrder[Exp[l ad ** a]] /. l -> Log[1 - lam];
sdag = SqueezeOperator[-xi, Infinity, "Ordering" -> "Normal"];
moved = BosonicExpSimplify[Exp[l ad ** a] ** Exp[y ad ** ad] ** Exp[-l ad ** a], {a, ad}];
vev = BosonicMatrixElement[{0, 0}, Evaluate[BosonicExpOrder[Exp[x a ** a] ** Exp[yy ad ** ad]]]];
yy1 = moved /. Exp[q_] :> q /. GeneralizedPower[NonCommutativeMultiply, ad, 2] -> 1;
detC = Det[{{Cosh[r], (lam - 1) E^(-I phi) Sinh[r]}, {(lam - 1) E^(I phi) Sinh[r], Cosh[r]}}];
{nA, sdag, moved, vev, FullSimplify[TrigToExp[(Sech[r] vev /. {x -> E^(-I phi) Tanh[r]/2, yy -> yy1} /. {y -> E^(I phi) Tanh[r]/2, l -> Log[1 - lam]}) - detC^(-1/2)], r > 0 && 0 < lam < 1 && Element[phi, Reals]]}
```
Output: `{NormalOrdered[E^(-(lam*SuperDagger[\[FormalA]1]**\[FormalA]1))], <S(-xi) normal form>, E^(y*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]*(Cosh[2*Sqrt[l^2]] + (l*Sinh[2*Sqrt[l^2]])/Sqrt[l^2])), 1/Sqrt[1 - 4*x*yy], 0}`. Matches: <0|S A S^dag|0> = sech r (1 - (1 - lambda)^2 tanh^2 r)^{-1/2} equals det C^{-1/2} exactly. The steps are: S^dag|0> = sech^{1/2} r e^{tau a*^2/2}|0> (from `sdag`), A moves past the squeeze wing as (1 - lambda)^2 (from `moved`, with l = log(1 - lambda)), and the vacuum value of e^{x a^2} e^{y a*^2} is (1 - 4xy)^{-1/2} (from `vev`).

### Fermionic canonical relations
Topic: Fermionic canonical transformations
Appears in: B-056 (also B-057)

A fermionic canonical transformation has blocks with Phi Phi* + Psi Psi* = I and Phi Psi' + Psi Phi' = 0, the plus signs replacing the bosonic minus signs.

$$ \hat b=\hat\Phi\hat a+\hat\Psi\hat a^*,\quad \hat\Phi\hat\Phi^*+\hat\Psi\hat\Psi^*=\hat I,\quad\hat\Phi\hat\Psi'+\hat\Psi\hat\Phi'=0,\qquad \hat S^{-1}=\begin{pmatrix}\hat\Phi^*&\hat\Psi'\\\hat\Psi^*&\hat\Phi'\end{pmatrix} $$

Audit: Correct. Checked by: exact fermion matrices (2 modes; U = exp(itH) for complex pairing and complex Hermitian hopping; Phi, Psi read off U c U^dag by trace projection; all relations hold).
QF code (exact): none. No fermionic layer.

### Fermionic Gaussian implementer kernel
Topic: Fermionic canonical transformations
Appears in: B-060

The fermionic implementer's kernel is a Grassmann Gaussian with the same block structure as the bosonic one, but with prefactor det(Phi Phi*)^{+1/4}.

$$ K_{\hat U}=\theta\det(\hat\Phi\hat\Phi^*)^{1/4} \exp\!\left[-\tfrac12\left(a\bar{\hat\Psi}\hat\Phi^{-1}a-2a^*\hat\Phi^{-1}a+a^*\hat\Phi^{-1}\hat\Psi a^*\right)\right], \qquad N_{\hat U}=K_{\hat U}e^{-a^*a} $$

Conditions: two-sided Phi^{-1}; Grassmann order as displayed.
Audit: Correct. Checked by: exact fermion matrices (2 modes): the Grassmann generating function Sum_{m,n} <m|U|n> (a*)^m ((a*)^n)* of the Jordan-Wigner U equals the displayed exponential times one scalar (all 16 coefficients), and that scalar has modulus det(Phi Phi*)^{1/4}.
QF code (exact): none. No fermionic layer.

### Fermionic paired vacuum
Topic: Gaussian vacua and overlaps
Appears in: B-061

The fermionic transformed vacuum is a pair exponential with skew pairing matrix Phi^{-1} Psi, normalized by det(Phi Phi*)^{1/4}.

$$ F(a^*)=\theta\det(\hat\Phi\hat\Phi^*)^{1/4} \exp[-\tfrac12a^*\hat\Phi^{-1}\hat\Psi a^*] $$

Audit: Correct. Checked by: exact fermion matrices (2 modes: <11|U|00>/<00|U|00> = -(Phi^{-1} Psi)_12).
QF code (exact): none. No fermionic layer.

### Fermionic transformation of normal symbols with mixing
Topic: Fermionic canonical transformations
Appears in: B-063

For a mixing fermionic transformation the conjugated normal symbol is a Grassmann Gaussian integral of the original one, with the displayed matrix and determinant factor.

$$ \widehat{\mathcal L}=\begin{pmatrix}\bar{\hat\Phi}\bar{\hat\Psi}^{-1}&-\hat I\\\hat I&-\hat\Phi\hat\Psi^{-1}\end{pmatrix},\qquad \det(\hat P_n\hat\Psi\hat\Psi^*\hat P_n)^{1/2} $$

Audit: Unchecked: the integral kernel (5.19)-(5.20) is described but not transcribed; the passive case (substitution) is correct by the same argument as B-053. Checked by: prior audit only for the mixing part.
QF code (exact): none. No fermionic layer.

### General quadratic expression
Topic: Quadratic propagators and Riccati
Appears in: B-065

The general quadratic Hamiltonian has pair-creation, pair-annihilation and number-type kernels plus linear and constant terms; only the symmetric (bosons) or antisymmetric (fermions) part of the pair kernels contributes.

$$ \hat H=\tfrac12\int[A(x,y)\hat a^*(x)\hat a^*(y)+B(x,y)\hat a(x)\hat a(y)+2C(x,y)\hat a^*(x)\hat a(y)]\,dxdy +\hat a^*(f)+\hat a(g)+k\hat I $$

Audit: Correct as a definition (A, B symmetric for bosons, skew for fermions). Checked by: QF exact (bosons: a1 a2 = a2 a1 and a1* a2* = a2* a1*, so the antisymmetric part of A, B drops out) and hand derivation (fermions: the symmetric part drops out by anticommutation).
QF code (exact):
```wl
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
BosonicNormalOrder[#, {a1, a1d, a2, a2d}] & /@ {a1 ** a2 - a2 ** a1, a1d ** a2d - a2d ** a1d}
```
Output: `{0, 0}`. Matches.

### Bound for exponentials of bounded perturbations
Topic: Wick, Dyson and auxiliary fields
Appears in: B-068

If A is a self-adjoint operator plus a bounded perturbation T, then e^{itA} grows at most like e^{|t| ||T||}.

$$ \|e^{it\hat A}\|\le e^{|t|\|\hat T\|},\qquad \hat A=\hat S+\hat T,\ \hat S=\hat S^* $$

Audit: Correct. Checked by: hand derivation: in the interaction-picture series the n-th term is bounded by (|t| ||T||)^n / n!.
QF code (exact): none. An operator-norm bound in a general Hilbert space; no second-quantization call is involved.

### Classical generator and Heisenberg action of a quadratic Hamiltonian
Topic: Quadratic propagators and Riccati
Appears in: B-069, B-071, B-072

A quadratic Hamiltonian evolves the fields linearly: conjugation by e^{itH} acts on (a, a*) by the matrix exponential of the classical generator L, plus a constant shift g from the linear terms.

$$ \widehat{\mathcal L}=\begin{pmatrix}-\hat C&-\hat A\\\bar{\hat A}&\bar{\hat C}\end{pmatrix},\quad e^{it\widehat{\mathcal L}}=\begin{pmatrix}\hat\Phi&\hat\Psi\\\bar{\hat\Psi}&\bar{\hat\Phi}\end{pmatrix},\quad \hat U(t)\binom{\hat a}{\hat a^*}\hat U(t)^{-1} =e^{it\widehat{\mathcal L}}\binom{\hat a}{\hat a^*}+\binom g{g^*}\hat I,\quad \binom g{g^*}=i\int_0^t e^{i\tau\widehat{\mathcal L}}\,d\tau\binom{-f}{f^*} $$

Conditions: U(t) = e^{itH}; the Hilbert-Schmidt criterion on F(t) in B-069 is a source theorem, Unchecked.
Audit: Correct. Checked by: QF exact (two modes, symbolic A, conj(A), C, f) plus hand derivation: d/dt (U a U^{-1}) = i U [H, a] U^{-1} with [H, a] = -C a - A a* - f and [H, a*] = conj(A) a + C' a* + conj(f), a linear system whose solution is the displayed formula.
QF code (exact):
```wl
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
v = {a1, a1d, a2, a2d}; an = {a1, a2}; cr = {a1d, a2d};
amat = {{A11, A12}, {A12, A22}}; abar = {{Ab11, Ab12}, {Ab12, Ab22}}; cmat = {{c11, c12}, {c21, c22}};
h = Sum[amat[[i, j]] cr[[i]] ** cr[[j]] + abar[[i, j]] an[[i]] ** an[[j]] + 2 cmat[[i, j]] cr[[i]] ** an[[j]], {i, 2}, {j, 2}]/2 + {f1, f2} . cr + {fb1, fb2} . an;
{Expand[BosonicNormalOrder[h ** # - # ** h, v] & /@ an + cmat . an + amat . cr + {f1, f2}],
 Expand[BosonicNormalOrder[h ** # - # ** h, v] & /@ cr - abar . an - Transpose[cmat] . cr - {fb1, fb2}]}
```
Output: `{{0, 0}, {0, 0}}`. Matches: [H, a] + C a + A a* + f and [H, a*] - conj(A) a - C' a* - conj(f) vanish exactly. The full conjugation `BosonicExpSimplify[Exp[I t H] ** a ** Exp[-I t H]]` returns unevaluated for this H.

### Quadratic propagator vacuum coefficient (corrected (6.23))
Topic: Quadratic propagators and Riccati
Appears in: B-071 (B-077 cites the same scalar terms)

The vacuum amplitude of e^{itH} is det(Phi e^{itC})^{-1/2} times the exponential of a time integral of terms quadratic and linear in the shift g, the quadratic one with coefficient 1/2.

$$ c(t)=\det(\hat\Phi e^{it\hat C})^{-1/2} \exp\!\left[i\int_0^t\left(\tfrac12g\hat\Phi'^{-1}\bar{\hat A}\hat\Phi^{-1}g-f^*\hat\Phi^{-1}g\right)d\tau\right] $$

Conditions: square root continued from c(0) = 1; constant k adds e^{itk}.
Audit: Correct as corrected in the atlas (the printed source omits the 1/2). Checked by: exact WL computation: in the Bargmann realization (a acts as d/dz, a* as z), U(t)|0> = c(t) exp(z R z/2 + l z) with R = -Phi^{-1} Psi, l = -Phi^{-1} g from exact e^{itL}; the residual of d/dt - iH on this vector vanishes identically (one mode, A = 3/10 + 2i/5, C = 1, f = 1/2 + i/3, symbolic t: residual 0; two modes, complex A, Hermitian C, complex f at t = 7/10 with 60-digit arithmetic: residual 0 to 58 digits); the printed coefficient (no 1/2) leaves a constant residual of about 0.03 in both cases. Forcing line: the constant term of (1/2) conj(A) d^2/dz^2 acting on e^{lz} is (1/2) conj(A) l^2.
QF code (exact): none. `BosonicExpOrder[Exp[I t H]]` returns its input unchanged for a quadratic H with linear terms, and `BosonicExpSimplify[Exp[I t H] ** a ** Exp[-I t H]]` stays unevaluated, so QF gives no exact vacuum amplitude for e^{itH}.

### Logarithmic derivative of a determinant
Topic: Quadratic propagators and Riccati
Appears in: B-073

The time derivative of log det M is the trace of M^{-1} times the derivative of M (Jacobi's formula).

$$ \frac{d}{dt}\log\det \hat M=\operatorname{Tr}(\hat M^{-1}\dot{\hat M}) $$

Audit: Correct (Jacobi's formula; determinant-class assumptions as stated). Checked by: hand derivation: log det M = Tr log M, differentiated under the trace.
QF code (exact): none. A matrix identity; no second-quantization call is involved.

### Fermionic quadratic propagator
Topic: Fermionic canonical transformations
Appears in: B-074

For fermions the propagator kernel is det(Phi e^{itC})^{+1/2} times a Grassmann Gaussian built from the same classical generator.

$$ K_{e^{it\hat H}}=\det(\hat\Phi e^{it\hat C})^{1/2} \exp\!\left[-\tfrac12(a\ a^*) \begin{pmatrix}\bar{\hat\Psi}\hat\Phi^{-1}&\hat\Phi'^{-1}\\-\hat\Phi^{-1}&\hat\Phi^{-1}\hat\Psi\end{pmatrix}\binom a{a^*}\right] $$

Audit: Correct. The exponent equals that of B-060 because a Phi'^{-1} a* = -a* Phi^{-1} a for Grassmann variables; the same generator L holds for fermions. Checked by: exact fermion matrices (2 modes: Phi, Psi from U c U^dag equal the blocks of exp(itL); <00|U|00> = det(Phi e^{itC})^{1/2} at t = 0.7).
QF code (exact): none. No fermionic layer.

### Resolvent kernel of a number-type Hamiltonian
Topic: Normal symbols and normal-order exponentials
Appears in: B-075

The resolvent of H = a* C a is an integral over xi in (0, 1) of the kernel of xi^H weighted by xi^{-z-1}, valid only where that integral converges.

$$ K_{(\hat H-z\hat I)^{-1}}=\int_0^1\exp(a^*\xi^{\hat C}a)\,\xi^{-z-1}d\xi $$

Conditions: the atlas's added restriction: needs Re z below the spectrum of H (Re z < 0 for C >= 0), then analytic continuation.
Audit: Correct as qualified in the atlas. Checked by: QF exact (one mode: xi^{cN} has normal form :exp((xi^c - 1) a*a):, so its kernel is exp(a* xi^c a)) and exact WL computation (int_0^1 xi^{cn - z - 1} d xi = 1/(cn - z) only under Re(cn) > Re(z)).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicExpOrder[Exp[s ad ** a]] /. s -> c Log[xi], Integrate[xi^(c n - z - 1), {xi, 0, 1}]}
```
Output: `{NormalOrdered[E^((-1 + xi^c)*SuperDagger[\[FormalA]1]**\[FormalA]1)], ConditionalExpression[(c*n - z)^(-1), Re[c*n] > Re[z]]}`. Matches, with the convergence condition the atlas adds.

### Weighted sufficient conditions for a quadratic normal form
Topic: Infinite modes and domains
Appears in: B-078

When C is multiplication by c(x), finiteness of these two weighted integrals suffices for the formal quadratic expression to be reduced to a quadratic normal form.

$$ \int\frac{|f(x)|^2}{|c(x)-z|}dx<\infty,\quad \iint\frac{|A(x,y)|^2}{|c(x)+c(y)-z|}dxdy<\infty $$

Audit: Unchecked: infinite-mode sufficient condition from a source theorem; no finite-dimensional content. Checked by: prior audit only.
QF code (exact): none. A continuum domain condition; no finite second-quantization content.

### Forced oscillator: completing the square
Topic: Forced and parametric oscillators
Appears in: B-081

A linearly forced oscillator is a free oscillator in the displaced field b = a + gf/omega, lowered by the constant g^2 |f|^2/omega.

$$ \hat H=\int\omega\hat b^*\hat b\,dk- g^2\int\frac{|f|^2}{\omega}\,dk\,\hat I,\qquad \hat a=\hat b-gf/\omega\,\hat I $$

Audit: Correct as corrected in the atlas (the source prints +). Checked by: QF exact (one mode, symbolic omega, g, f) plus hand derivation: b = U a U^{-1} for a displacement, so b*b has spectrum 0, 1, 2, ... and the ground energy is -g^2|f|^2/omega.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
h = w ad ** a + g (f ad + fb a);
{BosonicNormalOrder[w (ad + g fb/w) ** (a + g f/w) - g^2 f fb/w - h, {a, ad}],
 BosonicNormalOrder[w (ad + g fb/w) ** (a + g f/w) + g^2 f fb/w - h, {a, ad}]}
```
Output: `{0, (2*f*fb*g^2)/w}`. Matches with the minus sign (first entry 0); the printed plus sign leaves 2 g^2 |f|^2/omega (second entry).

### Regularized forced Hamiltonian
Topic: Forced and parametric oscillators
Appears in: B-083

Writing the forced Hamiltonian as the displaced number operator adds back the constant g^2 int |f|^2/omega, the term that diverges when only int |f|^2/omega^2 is finite.

$$ \hat H_{\mathrm{reg}}=\int\omega(k)\left(\hat a^*+g\frac{f^*}{\omega}\hat I\right) \left(\hat a+g\frac f\omega\hat I\right)dk $$

Audit: Correct and consistent with the corrected B-081: expanding gives the original H plus g^2 int |f|^2/omega. Checked by: QF exact (one mode).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
BosonicNormalOrder[w (ad + g fb/w) ** (a + g f/w) - (w ad ** a + g (f ad + fb a)) - g^2 f fb/w, {a, ad}]
```
Output: `0`. Matches.

### Bosonic canonical diagonalization
Topic: Bosonic canonical transformations
Appears in: B-084

A stable quadratic Hamiltonian (S plus or minus T positive) is a free Hamiltonian in squeezed modes, with frequency sqrt(S^2 - T^2) in one mode and ground-energy shift (Omega - S)/2.

$$ \hat H=\tfrac12(\hat a^*\hat T\hat a^*+\hat a\bar{\hat T}\hat a+2\hat a^*\hat S\hat a)\ \to\ \hat a^*\hat\Omega\hat a+k\hat I,\qquad \hat S\pm \hat T>\mu \hat I>0 $$

Audit: Correct in the finite one-mode case: Omega = sqrt(S^2 - T^2), k = (Omega - S)/2. The infinite-mode theorem is Unchecked. Checked by: QF exact (symbolic s > t > 0) plus hand derivation: b = cosh r a + sinh r a* with cosh^2 r - sinh^2 r = 1 is a squeezed annihilator, so Omega b*b + k has spectrum Omega n + k.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
h = (t ad ** ad + t a ** a + 2 s ad ** a)/2;
om = Sqrt[s^2 - t^2]; ch = Sqrt[(s/om + 1)/2]; sh = Sqrt[(s/om - 1)/2];
{b, bd} = {ch a + sh ad, ch ad + sh a};
Simplify[{BosonicNormalOrder[b ** bd - bd ** b, {a, ad}], BosonicNormalOrder[om bd ** b + (om - s)/2 - h, {a, ad}]}, s > t > 0]
```
Output: `{1, 0}`. Matches: [b, b*] = 1 and H = Omega b* b + (Omega - S)/2 exactly.

### Riccati fixed point for the diagonalizing transformation
Topic: Quadratic propagators and Riccati
Appears in: B-085

The ground state of H = (q beta q + p gamma p)/2 is the Gaussian exp(i q Z q/2), with Z = iW and W the positive solution of W gamma W = beta, given in closed form by the displayed square roots.

$$ \hat Z=i\hat\gamma^{-1/2}(\hat\gamma^{1/2}\hat\beta\hat\gamma^{1/2})^{1/2}\hat\gamma^{-1/2},\qquad \hat\beta=\hat S+\hat T,\quad\hat\gamma=\hat S-\hat T $$

Audit: Correct. Checked by: hand derivation (W gamma W = gamma^{-1/2} X^{1/2} X^{1/2} gamma^{-1/2} = beta for X = gamma^{1/2} beta gamma^{1/2}) and exact WL computation (exact rational 2x2 positive beta, gamma: W gamma W - beta and W - W^T reduce to 0 with `RootReduce`, and W has positive eigenvalues).
QF code (exact): none. A real matrix equation; no second-quantization call is involved.

### Fermionic diagonalization conditions
Topic: Fermionic canonical transformations
Appears in: B-086

A fermionic quadratic Hamiltonian is diagonalized by a proper transformation when S and S^2 - [S, T] - T^2 are bounded below by a positive constant.

$$ \hat S>\mu \hat I>0,\qquad \hat S^2-[\hat S,\hat T]-\hat T^2>\mu \hat I>0 $$

Audit: Consistent in finite dimension: with T real skew, S^2 - [S,T] - T^2 = (S+T)(S+T)^T, so the second condition is a gap condition (no zero quasiparticle energy). The infinite-mode theorem is Unchecked. Checked by: hand derivation: (S + T)(S - T) = S^2 - ST + TS - T^2 and (S + T)^T = S - T.
QF code (exact): none. No fermionic layer.

### Square-root perturbation stays in the ideal
Topic: Infinite modes and domains
Appears in: B-087

The correction X = (C^2 + K)^{1/2} - C solves X^2 + CX + XC = K, which the source uses to show X is as small (trace class, Hilbert-Schmidt) as K.

$$ \hat X=(\hat C^2+\hat K)^{1/2}-\hat C,\qquad \hat X^2+\hat C\hat X+\hat X\hat C=\hat K $$

Audit: Correct for the identity. The ideal statement is a source theorem. Checked by: hand derivation: (X + C)^2 = C^2 + K.
QF code (exact): none. An operator identity; no second-quantization call is involved.

### Scattering unitary in the Thirring model
Topic: Infinite modes and domains
Appears in: B-092

The Thirring model's scattering operator is a phase times the exponential of the product of the two zero-momentum charges.

$$ \hat S=\theta e^{-\pi igF(0)\hat Q_1'(0)\hat Q_2'(0)} $$

Audit: Unchecked: continuum model result, not reconstructed. Checked by: prior audit only.
QF code (exact): none. A continuum fermionic model; no fermionic layer.

### Pair potential on the n-particle sector
Topic: Wick, Dyson and auxiliary fields
Appears in: B-093

A quartic normal-ordered interaction acts on the n-particle sector as a sum over particle pairs of the pair potential; the factor 2 corresponds to the quartic term written without the usual 1/2.

$$ \hat H_n=-\sum_i\Delta_i+2\sum_{i<j}u(x_i-x_j)\hat I $$

Audit: Convention: the 2 corresponds to a quartic term Sum u a*a*aa without the usual 1/2; the source's quartic coefficient is not transcribed, so the factor is Unchecked. Checked by: QF exact for the general rule (two sites: Sum_{ij} u_ij a_i* a_j* a_j a_i gives 2 u_12 on |1,1> and 2 u_11 on |2,0>).
QF code (exact):
```wl
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
cr = {a1d, a2d}; an = {a1, a2}; u = {{u11, u12}, {u12, u22}};
vq = Sum[u[[i, j]] cr[[i]] ** cr[[j]] ** an[[j]] ** an[[i]], {i, 2}, {j, 2}];
{BosonicMatrixElement[{{1, 1}, {1, 1}}, vq], BosonicMatrixElement[{{2, 0}, {2, 0}}, vq], BosonicMatrixElement[{{1, 1}, {2, 0}}, vq]}
```
Output: `{2*u12, 2*u11, 0}`. Matches the rule (a quartic term without 1/2 gives 2 u per pair).

### Dyson series with the e^{+itH} convention
Topic: Wick, Dyson and auxiliary fields
Appears in: B-094

With Berezin's sign convention e^{+itH}, the interaction-picture propagator is the time-ordered exponential with powers of +i.

$$ \hat U(t,\tau)=\sum_{n\ge0}\frac{i^n}{n!}\int_\tau^t\cdots\int_\tau^t T[\hat V(t_1)\cdots \hat V(t_n)]\,dt_1\cdots dt_n $$

Conditions: U(t, tau) = e^{-itH0} e^{i(t-tau)H} e^{i tau H0}, V(t) = e^{-itH0} V e^{itH0}; formal.
Audit: Correct. Checked by: hand derivation: dU/dt = i V(t) U with U(tau, tau) = I.
QF code (exact): none. A formal time-ordered series; no exact QF call evaluates it.

### Free time dependence and contraction kernel
Topic: Wick, Dyson and auxiliary fields
Appears in: B-095

In Berezin's convention the free annihilator acquires the phase e^{+it omega}, the creator e^{-it omega}, and their unequal-time commutator is the scalar D.

$$ \hat a(p,t)=e^{it\omega(p)}\hat a(p),\quad \hat a^*(p,t)=e^{-it\omega(p)}\hat a^*(p),\qquad D(p_1,t_1\mid p_2,t_2)=\delta(p_1-p_2)e^{i(t_1-t_2)\omega(p_1)} $$

Audit: Correct (consistent with V(t) = e^{-itH0} V e^{itH0} of B-094). Checked by: QF exact (one mode; the commutator D is checked in the next entry).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
FullSimplify[{BosonicExpSimplify[Exp[-I t w ad ** a] ** a ** Exp[I t w ad ** a], {a, ad}],
  BosonicExpSimplify[Exp[-I t w ad ** a] ** ad ** Exp[I t w ad ** a], {a, ad}]}, Element[t | w, Reals]]
```
Output: `{\[FormalA]1*E^(I*t*w), SuperDagger[\[FormalA]1]/E^(I*t*w)}`. Matches.

### Retarded contraction and the two-field time-ordered product
Topic: Wick, Dyson and auxiliary fields
Appears in: B-097

The retarded contraction is D when the annihilator's time is later and 0 otherwise, so a two-field time-ordered product is its normal product plus that contraction.

$$ \Delta(p_1,t_1\mid p_2,t_2)= \begin{cases}D(p_1,t_1\mid p_2,t_2),&t_1>t_2,\\0,&t_1\le t_2,\end{cases}\qquad T[\hat a^*(1)\hat a(2)]=\hat a^*(1)\hat a(2)+\Delta(2\mid1)\hat I $$

Audit: Correct for bosons (for t2 > t1, a(2)a*(1) = a*(1)a(2) + D(2|1)); for fermions the term is -Delta(2|1). Checked by: QF exact (bosons) and hand derivation (fermions: the anticommutator replaces the commutator and time ordering carries the sign).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{x1, x2} = {E^(-I t1 w) ad, E^(I t2 w) a};
BosonicNormalOrder[x2 ** x1 - x1 ** x2, {a, ad}]
```
Output: `E^((-I)*t1*w + I*t2*w)`. Matches: a(2) a*(1) - a*(1) a(2) = e^{i(t2 - t1) omega} = D(2|1).

### Wick's theorem as an exponential of the contraction operator
Topic: Wick, Dyson and auxiliary fields
Appears in: B-099

The normal symbol of a time-ordered product is obtained by applying the exponential of the second-order contraction operator to the plain product of symbols.

$$ \widehat{\mathcal D}=\int\frac{\delta}{\delta a(q,\tau)}\Delta(q,\tau\mid p,t) \frac{\delta}{\delta a^*(p,t)}\,dp\,dq\,dt\,d\tau,\qquad T\text{-symbol}=e^{\widehat{\mathcal D}}F $$

Audit: Correct for bosons. For fermions the derivative sides in the atlas wording are ambiguous and one reading gives the wrong sign (see Errors found): the right derivative must be the one in a* and the left derivative the one in a, giving D(a*_p a_q) = -Delta_qp. Checked by: QF exact (bosons, four fields a*(t=1), a(t=4), a*(t=2), a(t=5), including both double contractions) and exact WL computation (fermions, two-mode Grassmann algebra, three readings compared).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
toSym[e_] := e /. {ad -> zb, a -> z} /. GeneralizedPower[NonCommutativeMultiply, x_, n_] :> x^n /. NonCommutativeMultiply -> Times;
fld = {{ad, -1, 1}, {a, 1, 4}, {ad, -1, 2}, {a, 1, 5}};
tProd = toSym[BosonicNormalOrder[NonCommutativeMultiply @@ (E^(I #[[2]] #[[3]] w) #[[1]] & /@ SortBy[fld, -Last[#] &]), {a, ad}]];
{ann, cre} = Flatten[Position[fld[[All, 1]], #, {1}]] & /@ {a, ad};
dop[e_] := Sum[If[fld[[q, 3]] > fld[[p, 3]], E^(I (fld[[q, 3]] - fld[[p, 3]]) w), 0] D[e, x[q], x[p]], {q, ann}, {p, cre}];
eD = Sum[Nest[dop, Product[x[j], {j, 4}], k]/k!, {k, 0, 2}];
Expand[tProd - (eD /. x[j_] :> E^(I fld[[j, 2]] fld[[j, 3]] w) If[fld[[j, 1]] === a, z, zb])]
```
Output: `0`. Matches: the normal symbol of the time-ordered product minus e^{D} applied to the product of the four symbols, with the free phases inserted afterwards, is exactly 0.

### Normal symbol of interaction-picture evolution
Topic: Wick, Dyson and auxiliary fields
Appears in: B-100

Applying Wick's theorem term by term to the Dyson series gives the normal symbol of the interaction-picture propagator as e^{D} acting on the exponential of the integrated interaction symbol.

$$ N_{\hat U,1}(t,\tau\mid a^*,a)=e^{\widehat{\mathcal D}} \exp\!\left[i\int_\tau^t V(s\mid a^*,a)ds\right] $$

Audit: Correct as a formal identity: B-094 plus B-099 term by term. Checked by: hand derivation.
QF code (exact): none. A formal series over arbitrary interactions; no exact QF call evaluates it.

### Contractions acting on an exponential source
Topic: Wick, Dyson and auxiliary fields
Appears in: B-101

An exponential of linear Grassmann sources is an eigenfunction of the contraction operator, picking up the factor e^{f* Delta f}.

$$ e^{\widehat{\mathcal D}}e^{f^*a+a^*f}=e^{f^*\Delta f+a^*f+f^*a} $$

Audit: Correct under the reading of B-099 that gives the fermionic sign; fails under the other reading. Checked by: exact WL computation (explicit Grassmann algebra, two modes, complex Delta, Grassmann sources: exact match).
QF code (exact): none. No fermionic layer.

### Determinant identity for a quadratic symbol
Topic: Wick, Dyson and auxiliary fields
Appears in: B-102

Contracting a fermionic quadratic exponential gives det(I + K Delta) times a quadratic exponential with the dressed kernel (K^{-1} + Delta)^{-1}.

$$ e^{\widehat{\mathcal D}}e^{-a^*\hat Ka} =\det(\hat I+\hat K\hat\Delta)\exp[-a^*(\hat K^{-1}+\hat\Delta)^{-1}a] $$

Audit: Correct (fermionic). Checked by: exact WL computation (explicit Grassmann algebra, two modes, complex non-commuting K and Delta: exact match under the correct reading).
QF code (exact): none. No fermionic layer.

### Auxiliary-field decoupling of a quartic term
Topic: Wick, Dyson and auxiliary fields
Appears in: B-103

Factorizing the quartic kernel through an auxiliary label lets a commuting Gaussian integral over z, z* turn the quartic exponential into a quadratic one.

$$ \int V_1(p_1,q_1\mid s)V_2(s\mid p_2,q_2)ds=iV(p_1,p_2\mid q_1,q_2) $$

Audit: Correct in mechanism: int e^{-z*z + u z + z* v} d^2z/pi = e^{uv} turns the factorized kernel into a quadratic coupling. Signs and the factor i in the full Appendix formula are Unchecked (not transcribed). Checked by: hand derivation: complete the square, z -> z + v, z* -> z* + u.
QF code (exact): none. The full Appendix formula is not transcribed, so there is nothing for QF to reproduce.

### Retarded loops vanish
Topic: Wick, Dyson and auxiliary fields
Appears in: B-104

Closed loops of retarded contractions with a time-local kernel vanish, so the auxiliary determinant equals 1.

$$ \operatorname{Tr}(\hat R\hat\Delta)^n=0,\qquad \det(\hat I+\hat R\hat\Delta)=1 $$

Audit: Correct as a formal result. Checked by: hand derivation: a closed loop needs t_1 > t_2 > ... > t_n > t_1, impossible, and n = 1 vanishes by Delta(t, t) = 0.
QF code (exact): none. A formal statement about continuum time kernels; no exact QF call applies.

### Auxiliary-field representation of the propagator
Topic: Wick, Dyson and auxiliary fields
Appears in: B-105

The normal symbol of the full propagator is a Gaussian average over the auxiliary field of a quadratic normal-ordered exponential; with no interaction it reduces to the normal symbol of e^{it omega N}.

$$ N_{e^{it\hat H}}=\int\exp\!\left[a^*\bigl(e^{it\hat\omega}(\hat I-\hat T)-\hat I\bigr)a-zz^*\right]\Pi dz\,dz^* $$

Audit: Checked only in the free limit T = 0, where it reduces to the normal symbol of e^{it omega N}, which is correct; the interacting T is not transcribed, Unchecked. Checked by: QF exact (free limit, one mode).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
BosonicExpOrder[Exp[I t w ad ** a]]
```
Output: `NormalOrdered[E^((-1 + E^(I*t*w))*SuperDagger[\[FormalA]1]**\[FormalA]1)]`. Matches the free limit.

### Grand partition functions
Topic: Thermal states and traces
Appears in: B-106

The grand partition function is an auxiliary-field average of 1/det(I - M) for bosons and det(I + M) for fermions, with M the one-particle Boltzmann factor dressed by the interaction.

$$ \Xi_{B}=\int\det(\hat I-\hat M)^{-1}e^{-zz^*}\Pi dz\,dz^*,\qquad \Xi_{F}=\int\det(\hat I+\hat M)e^{-zz^*}\Pi dz\,dz^*,\quad \hat M=e^{-\beta(\hat\omega+\mu\hat I)}(\hat I-\hat T) $$

Audit: Correct in the free limit (T = 0 gives prod 1/(1 - e^{-beta(omega+mu)}) and prod (1 + e^{-beta(omega+mu)}), the traces of e^{-beta(H0 + mu N)} with the source's plus sign); interacting case Unchecked. Checked by: QF exact (bosonic free limit, one mode) and hand derivation (fermionic free limit: one mode has the two levels 1 and e^{-beta(omega+mu)}).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicExpOrder[Exp[-x ad ** a]] /. x -> beta (w + mu),
 Integrate[Exp[(E^(-beta (w + mu)) - 1) (p^2 + q^2)]/Pi, {p, -Infinity, Infinity}, {q, -Infinity, Infinity}, Assumptions -> beta (w + mu) > 0]}
```
Output: `{NormalOrdered[E^((-1 + E^(-(beta*(mu + w))))*SuperDagger[\[FormalA]1]**\[FormalA]1)], 1 + (-1 + E^(beta*(mu + w)))^(-1)}`. Matches: 1 + 1/(e^{beta(omega+mu)} - 1) = 1/(1 - e^{-beta(omega+mu)}). Written with -(w + mu) directly in the exponent, `BosonicExpOrder` returns its input unchanged, so a single symbol x stands for beta (w + mu) inside the call.

### Free Fock-space traces of a second-quantized one-particle map
Topic: Thermal states and traces
Appears in: B-107

The trace over bosonic Fock space of the second quantization of M is 1/det(I - M), and over fermionic Fock space det(I + M).

$$ \operatorname{Tr}_{\mathcal F_{B}}\widehat\Gamma(\hat M)=\det(\hat I-\hat M)^{-1},\qquad \operatorname{Tr}_{\mathcal F_{F}}\widehat\Gamma(\hat M)=\det(\hat I+\hat M) $$

Conditions: bosons need ||M|| < 1 (and trace class in infinite dimension); fermions need trace class only.
Audit: Correct, including the M = 2I counterexample (divergent occupation sum). Checked by: QF exact (one mode: the normal form of m^N and the geometric sum; 2^n is not summable), exact WL computation (two modes, M = ((1/2, 1/8), (1/8, 1/4)): the Gaussian integral of the normal symbol exp(a*(M - I)a) is exactly 1/det(I - M)), hand derivation (Gamma(V M V^{-1}) = Gamma(V) Gamma(M) Gamma(V)^{-1}, so the trace is prod 1/(1 - m_i) over the eigenvalues) and exact fermion matrices (2 modes, complex non-Hermitian h: Tr Gamma(e^h) = det(I + e^h)).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicExpOrder[Exp[h ad ** a]], Sum[m^n, {n, 0, Infinity}], SumConvergence[2^n, n]}
```
Output: `{NormalOrdered[E^((-1 + E^h)*SuperDagger[\[FormalA]1]**\[FormalA]1)], (1 - m)^(-1), False}`. Matches (with m = e^h the one-mode Boltzmann factor).

### Number-type exponential in normal order
Topic: Normal symbols and normal-order exponentials
Appears in: B-107 (special cases in B-013, B-028, B-105)

The exponential of a quadratic number-conserving form a* h a is the normal-ordered exponential of a*(e^h - I)a.

$$ e^{\hat a^*\hat h\hat a}=:\!\exp[\hat a^*(e^{\hat h}-\hat I)\hat a]\!: $$

Audit: Correct. Checked by: QF exact (one mode) and exact WL computation (two modes, non-normal complex h: in the Bargmann realization z.h.d/dz generates z -> e^{h^T} z, so e^{a* h a} sends exp(z.w) to exp(z.e^h w), residual 0; hence K = exp(a* e^h a) and N = K e^{-a*a}).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
BosonicExpOrder[Exp[h ad ** a]]
```
Output: `NormalOrdered[E^((-1 + E^h)*SuperDagger[\[FormalA]1]**\[FormalA]1)]`. Matches. For two modes with a general matrix h, `BosonicExpOrder` returns its input unchanged (QF has only the two-mode beam-splitter form).

---

## Non-formula records

- B-001: Fock spaces and dimensions (2^N fermionic, infinite bosonic); correct.
- B-002: involution, transpose, regularized determinants; definitions and a source theorem.
- B-007: number operator and the form domain of N^{1/2}; domain statement.
- B-008: commutant and irreducibility; source theorem.
- B-009: uniqueness of the vacuum representation; source theorem.
- B-010: meaning of normal form analytically; definition.
- B-011: bilinear forms weaker than operators; domain statement.
- B-012: bounded operators admit normal form; source theorem.
- B-018: entire analyticity of generating series; source theorem.
- B-021: continual integrals as projection-independent limits; definition.
- B-024: Gaussian moments (holomorphic and antiholomorphic degrees must match, sum of pairings); prose, consistent with int z*^k z^l e^{-|z|^2} d^2z/pi = delta_kl k!.
- B-025: analytic (Bargmann) Hilbert-space realization; source theorem.
- B-035: infinite Grassmann functionals; source construction.
- B-037: signed fermionic Gaussian moments; prose (permutation parity), consistent with B-101.
- B-038: fermionic Fock inner product as a Grassmann integral; prose.
- B-039: fermionic kernel composition has the B-026/B-027 form; prose.
- B-046: bosonic implementability criterion (Psi Hilbert-Schmidt, f in L); source theorem.
- B-047: vacuum criterion for an implementer; source lemma.
- B-052: projective phases and metaplectic double cover; prose (a 2 pi rotation e^{2 pi i (N + 1/2)} = -I illustrates it).
- B-057: fermionic boundedness and orthogonal geometry; prose (relation itself covered in B-056).
- B-058: fermionic implementability with property A; source theorem.
- B-059: singular fermionic charts and filled sectors; source lemmas.
- B-062: fermionic spin representation; prose.
- B-064: fermionic Gaussian symbol conjugation, C replaced by C1; the block C is not transcribed.
- B-066: vacuum-domain necessity for quadratic normal forms; domain statement.
- B-067: self-adjointness of quadratic normal forms; source theorem.
- B-070: determinant-class phase condition (Phi e^{itC} - I trace class); source lemma.
- B-076: necessary condition for a singular quadratic Hamiltonian; source theorem.
- B-077: sufficient condition with the missing -I; the correction is consistent with B-070 and with t = 0.
- B-079: weaker conditions and infinite additive constants; source statement.
- B-080: continuum translation-invariant pairing is not Hilbert-Schmidt; domain statement.
- B-082: weighted coupling regimes and domain dependence; source statement.
- B-088: inequivalent CAR representations in the continuum; source statement.
- B-089: formal CAR rearrangement giving a false k = 0; source footnote.
- B-090: approximation before singular limits; source statement.
- B-091: evolved fields without a Fock-space Hamiltonian; source statement.
- B-096: graded time ordering and equal-time normal-order convention; definition.
- B-098: Wick's theorem as a sum over contractions; prose (formula form in B-099).
- B-108: coverage limits of the book; coverage note.

