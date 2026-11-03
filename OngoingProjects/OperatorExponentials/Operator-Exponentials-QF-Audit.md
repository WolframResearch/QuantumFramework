# Operator Exponentials: Formula Audit and Exact QuantumFramework Code

Each entry below is one distinct formula from the atlas [Bosonic-Fermionic-Operator-Exponentials.md](Bosonic-Fermionic-Operator-Exponentials.md). It gives a short statement of what the formula says, the formula, every atlas ID where it appears, an audit verdict with how it was checked, and the QuantumFramework (QF) code that produces it with the kernel's output. When several books state the same identity, the entry shows one version first and each other book's version under "Also stated in". All 428 atlas IDs are covered: 326 formula entries, of which 173 have an exact QF cell, plus 74 records that carry no formula.

Everything here is exact. No truncated Fock matrix is used anywhere, neither as QF code nor as a check. A QF cell uses only the symbolic second-quantization algebra (`**`, `FieldVariables`, `BosonicNormalOrder`, `BosonicExpOrder`, `BosonicExpSimplify`, `BosonicExpBraid`, `BosonicMatrixElement`, `BosonicVEV`, `NormalOrdered`, and `DisplacementOperator` / `SqueezeOperator` with size `Infinity`). Where QF has no exact route, the entry says "QF code (exact): none" and names where QF stops; the formula is then checked by exact Wolfram computation (symbolic algebra, series to a stated order, exact integrals, Bargmann or position realizations acting on symbolic functions, exact 2×2 or 2N×2N group matrices, exact Jordan-Wigner matrices for a few fermion modes) or by a hand derivation. In a few places two exact closed forms are compared at random points to high precision because `FullSimplify` stalls; the entry says so.

## Loading the QF code

Every cell assumes these two lines have been run. Without the last `Needs` the second-quantization symbols stay unevaluated and raise no error.

```wl
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

## Conventions to fix before comparing with QF

- `SqueezeOperator[ξ, Infinity]` is exp((ξ̄ a² − ξ a†²)/2), the sign of the atlas's common convention (X-002) and of Gerry and Knight. Combescure and Robert's S(ξ) in C-025, C-027 and C-048 is QF's `SqueezeOperator[-ξ, Infinity]`.
- `DisplacementOperator[α, Infinity]` satisfies D(α)†aD(α) = a + α. Berezin's implementer, oriented as U a U⁻¹ = a + f, is QF's D(−f). Without `Infinity`, `DisplacementOperator[α]` is a truncated matrix; it stays exact only when held inside `BosonicMatrixElement` with integer indices.
- `DisplacementOperator[α, Infinity]` and `SqueezeOperator[ξ, Infinity]` are written in `\[FormalA]1` (from `FieldVariables[{1}]`), while `FieldVariables[]` returns `\[FormalA]`. QF treats every distinct formal symbol as its own mode, so mixing the two silently makes a two-mode expression. Build everything from `FieldVariables[{1}]`.

## QF problems found (exact layer)

- Fixed in commit fbebae87: `BosonicMatrixElement[{0, 0}, NormalOrdered[Exp[-a† ** a]]]` returned `Indeterminate` (a 0⁰ in `numberDiagonal`); it now returns 1 at the vacuum and 0 elsewhere.
- Open, minor: `BosonicMatrixElement[{0, 0}, 0^(-a† ** a)]` and `{1, 1}` of `0^(a† ** a - 1)` still give `Indeterminate`, because the fix covers a zero base only when the exponent is provably nonnegative.
- Open: `BosonicNormalOrder[{e1, e2}, vars]` with a list as first argument returns a wrong value with no message (for example `{{0 &}, a† ** a}` for `{a ** a†, a† ** a}`). The usage documents a single expression; mapping over the list gives the right answer.
- Missing rules (the call returns unevaluated or unchanged, never a wrong value):
  - `BosonicMatrixElement` in the coherent basis does not take exponentials of a or a† (`BosonicAlgebra.m:406`); with several modes it does not take exponentials in the Fock basis; with symbolic m it does not reduce {m, m} of D or {2m, 0} of S.
  - `BosonicExpOrder` leaves exp(x a² + y a†² + w a†a), exp(c a†a + linear), a general quadratic with linear terms, the Kerr exponential exp(λ (a†a)²), and exp(−(w + μ) a†a) with a sum as coefficient unchanged.
  - `BosonicExpSimplify` leaves conjugation of a sum of monomials, of a + f by a squeeze, by a mixed quadratic and linear generator, and two-mode conjugation by a general complex-symmetric squeeze unevaluated; it leaves `S(ξ) ** S(−ξ)` as `E^(X/2 − X/2)` uncombined.
  - `BosonicVEV` does not take a product of exponentials.

These gaps, and minimal ways to close them, are worked out in [Minimal-Bosonic-Design.md](Minimal-Bosonic-Design.md).

## Atlas issues found

No formula in the atlas has a wrong sign, factor or condition beyond the source errors it already flags. The corrections the atlas makes were rechecked exactly, and in each case the printed form fails while the corrected one holds: B-036, B-071, B-075, B-077, B-081, P-019, P-029, P2-013, P2-017, P2-023, P2-038, P2-054, C-050, C-085, G-028, G-077, G-079. These statements are ambiguous or loose and should be tightened:

- B-099 says "the first derivative is the right derivative". Read as the derivative written first (in a), it gives the wrong fermionic contraction sign and makes B-101 and B-102 fail; read as the right derivative in a* and the left one in a, both hold. The text should name the variable, not the position.
- X-006 refers to a matrix U that it never defines; it means the canonical block u of X-009.
- X-007 calls 0 ≤ M < I equivalent to ‖M‖ < 1. That holds for positive M; for a general M the trace converges when the spectral radius is below 1, and ‖M‖ < 1 is only sufficient.
- P-045, P-039 and P2-012 are right only under an index or basis reading the atlas leaves implicit: P-045's A is ⟨c_j c_i⟩, P-039's B is the transpose of the natural reading, and P2-012's K₂ pairs with a K₁ of flipped sign (the relabeling a → ia).
- G-050's fringe sign depends on an interferometer layout the atlas does not reproduce; it stays unchecked.

Every cell and check is stored in [drafts/exact](drafts/exact) (one `<section>-work` folder per book section, with the cells as run and their raw output).

## Contents

- [BCH and Hadamard](#bch-and-hadamard) (5)
- [Displacement](#displacement) (13)
- [Coherent states (general)](#coherent-states-general) (11)
- [Single-mode squeezing](#single-mode-squeezing) (17)
- [Two-mode squeezing](#two-mode-squeezing) (12)
- [SU(2) and passive mode mixing](#su-2-and-passive-mode-mixing) (14)
- [Bosonic canonical transformations](#bosonic-canonical-transformations) (20)
- [Fermionic canonical transformations](#fermionic-canonical-transformations) (17)
- [Gaussian vacua and overlaps](#gaussian-vacua-and-overlaps) (14)
- [Quadratic propagators and Riccati](#quadratic-propagators-and-riccati) (16)
- [Forced and parametric oscillators](#forced-and-parametric-oscillators) (13)
- [Normal symbols and normal-order exponentials](#normal-symbols-and-normal-order-exponentials) (19)
- [Weyl and antinormal symbols, characteristic functions](#weyl-and-antinormal-symbols-characteristic-functions) (32)
- [Grassmann calculus and fermionic coherent states](#grassmann-calculus-and-fermionic-coherent-states) (16)
- [Gaussian integrals and Pfaffians](#gaussian-integrals-and-pfaffians) (9)
- [Thermal states and traces](#thermal-states-and-traces) (24)
- [Wick, Dyson and auxiliary fields](#wick-dyson-and-auxiliary-fields) (13)
- [Nonlinear exponentials (number shift, Kerr, cats)](#nonlinear-exponentials-number-shift-kerr-cats) (18)
- [Fock matrix elements and generating functions](#fock-matrix-elements-and-generating-functions) (14)
- [Infinite modes and domains](#infinite-modes-and-domains) (12)
- [Lie-group decompositions and generalized coherent states](#lie-group-decompositions-and-generalized-coherent-states) (17)
- [Non-formula records](#non-formula-records)
- [Index of atlas IDs](#index-of-atlas-ids)

## Entries marked Unchecked

- General bosonic transformation kernel for normal symbols
- Fermionic transformation of normal symbols with mixing
- Fermionic analytic group action
- Fermionic Mehler kernel
- Forced parametric oscillator generating kernel and vacuum amplitude
- Periodic quadratic evolution and quasi-energies
- Disk star product to first order in h
- Fermionic Weyl translation expansion
- Quasiclassical partition-function bounds
- Weighted sufficient conditions for a quadratic normal form
- Scattering unitary in the Thirring model
- Periodic-orbit determinant in regularized traces
- Root-space coherent states and displacement-like generator
- Jacobians of classical-domain transformations
- Completeness criterion for discrete SU(1,1) subsystems

<a id="bch-and-hadamard"></a>
## BCH and Hadamard

<a id="hadamard-series-for-a-closed-linear-commutator-action"></a>
### Hadamard series for a closed linear commutator action
Appears in: X-005 (source records cited: P-034, P-035, P-041, P-042, B-071 to B-074, C-018 to C-024)

When commuting with $\hat Q$ maps a list of operators linearly into itself, conjugation by $e^{t\hat Q}$ is the matrix exponential of that linear map.

$$ [\hat Q,\hat R_i]=\sum_jL_{ij}\hat R_j\ \Longrightarrow\ e^{t\hat Q}\hat R e^{-t\hat Q}=e^{tL}\hat R $$

Conditions: $L$ scalar; time-dependent noncommuting generators need time ordering.
Audit: Correct. Checked by: hand derivation: $\operatorname{ad}_{\hat Q}^k\hat R=L^k\hat R$ because $L$ is scalar, so the Hadamard series sums to $e^{tL}\hat R$; and QF exact for the general single-mode quadratic $\hat Q=p\hat a^2+q\hat a^{\dagger2}+s\hat a^\dagger\hat a$, where $[\hat Q,(\hat a,\hat a^\dagger)]$ gives $L=\begin{pmatrix}-s&-2q\\2p&s\end{pmatrix}$. The atlas remark that the adjoint action loses the overall phase is right.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
q = p a ** a + qq ad ** ad + s ad ** a;
lhs = BosonicExpSimplify[Exp[t q] ** a ** Exp[-t q]];
rhs = First[MatrixExp[t {{-s, -2 qq}, {2 p, s}}] . {a, ad}];
{lhs, FullSimplify[lhs - rhs], Chop[Expand[(lhs - rhs) /. {p -> 0.3 + 0.2 I, qq -> -0.7 + 0.1 I, s -> 0.4 - 0.9 I, t -> 1.3}]]}
```
Output: `{\[FormalA]1*Cosh[Sqrt[(-4*p*qq + s^2)*t^2]] - (t*Sinh[Sqrt[(-4*p*qq + s^2)*t^2]]*(\[FormalA]1*s + 2*qq*SuperDagger[\[FormalA]1]))/Sqrt[(-4*p*qq + s^2)*t^2], 0, 0}`. Matches: the difference from the first row of $e^{tL}(\hat a,\hat a^\dagger)^T$ is zero symbolically and at a random complex point.

<a id="central-commutator-bch"></a>
### Central-commutator BCH
Appears in: P-002 (also used by P-010, P2-035); C-002 (also used in C-006, C-059, C-071); G-003

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

**Also stated in C-002 (Central-commutator BCH):**

When the commutator of two generators commutes with both, the exponential of their sum is the product of their exponentials times the scalar exponential of minus half the commutator.

$$ e^{\hat A+\hat B}=e^{-[\hat A,\hat B]/2}e^{\hat A}e^{\hat B} $$

Conditions: [A,[A,B]] = [B,[A,B]] = 0 on an invariant common core; anti-self-adjoint generators.
Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
BosonicExpOrder[Exp[u ad + v a]]
```
Output: `E^((u*v)/2)*E^(u*SuperDagger[\[FormalA]1])**E^(\[FormalA]1*v)`. Matches: with A = u a^dag, B = v a one has [A, B] = -u v, so e^{-[A,B]/2} = e^{u v/2}.

**Also stated in G-003 (Central-commutator BCH):**

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

<a id="duhamel-derivative-of-an-exponential"></a>
### Duhamel derivative of an exponential
Appears in: C-051

The derivative of the exponential of a moving generator is an average of the derivative conjugated along the exponential, not simply L' e^L.

$$ \frac d{dt}e^{\hat L(t)}=\int_0^1e^{s\hat L(t)}\dot{\hat L}(t)e^{(1-s)\hat L(t)}ds $$

Audit: Correct. Checked by: WL computation with finite 3-by-3 matrices (central difference against quadrature of the Duhamel integral).
QF code (exact): none. A general matrix identity; QF's second-quantization layer does not apply.

<a id="adjoint-exponential"></a>
### Adjoint exponential
Appears in: C-094

Conjugation by an exponential is the exponential of the adjoint map, which computes the conjugation whenever the commutator span closes.

$$ \operatorname{Ad}(e^X)=e^{\operatorname{ad}X},\qquad \operatorname{ad}X(Y)=[X,Y] $$

Audit: Correct; the atlas's rejection of a global diffeomorphism claim is right (on SU(2), exp(2 pi i sigma_3) = exp(0)). Checked by: QF exact for a closing case plus hand derivation.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
ada = BosonicExpSimplify[Exp[s ad ** a] ** a ** Exp[-s ad ** a]];
{ada, Simplify[ada/a - Exp[-s], Element[s, Reals]]}
```
Output: `{\[FormalA]1*(Cosh[Sqrt[s^2]] - (s*Sinh[Sqrt[s^2]])/Sqrt[s^2]), 0}`. Matches: e^{sN} a e^{-sN} = e^{-s} a = e^{ad(sN)} a, since [sN, a] = -s a.

<a id="hadamard-conjugation-series"></a>
### Hadamard conjugation series
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

<a id="displacement"></a>
## Displacement

<a id="displacement-normal-and-antinormal-factorization"></a>
### Displacement: normal and antinormal factorization
Appears in: X-001 (source records cited: P-002, P-003, P-010, G-015, G-017, B-050, C-006); P-010, P2-053; C-006; same identity in C-059 (infinite-mode form); G-015 (also G-101 at $\alpha=-i\eta$)

The displacement operator splits into a creation exponential times an annihilation exponential (or the reverse), and the only cost of reordering is a Gaussian scalar.

$$ \hat D(\alpha)=e^{\alpha \hat a^\dagger-\bar\alpha \hat a}=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}e^{-\bar\alpha \hat a}=e^{|\alpha|^2/2}e^{-\bar\alpha \hat a}e^{\alpha \hat a^\dagger},\qquad [\hat a,\hat a^\dagger]=\hat I $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{DisplacementOperator[al, Infinity, "Ordering" -> "Normal"], DisplacementOperator[al, Infinity, "Ordering" -> "Antinormal"], BosonicExpOrder[Exp[al ad - Conjugate[al] a]], BosonicExpOrder[Exp[al ad - Conjugate[al] a], "Ordering" -> "Antinormal"]}
```
Output: `{E^(al*SuperDagger[\[FormalA]1])**E^(-(\[FormalA]1*Conjugate[al]))/E^((al*Conjugate[al])/2), E^((al*Conjugate[al])/2)*E^(-(\[FormalA]1*Conjugate[al]))**E^(al*SuperDagger[\[FormalA]1]), <same two forms again from BosonicExpOrder>}`. Matches: both orderings are the atlas forms verbatim, from the named operator and from ordering the bare exponential.

**Also stated in P-010 (Displacement: normal and antinormal factorizations):**

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

**Also stated in C-006 (Normal and antinormal displacement factorization, displacement conjugation):**

The displacement operator equals its normal-ordered product times e^{-|alpha|^2/2}, or its antinormal product times e^{+|alpha|^2/2}, and it shifts the annihilator by -alpha.

$$ \hat D(\alpha)=e^{\alpha\cdot \hat a^\dagger-\bar\alpha\cdot \hat a}=e^{-|\alpha|^2/2}e^{\alpha\cdot \hat a^\dagger}e^{-\bar\alpha\cdot \hat a}=e^{|\alpha|^2/2}e^{-\bar\alpha\cdot\hat a}e^{\alpha\cdot\hat a^\dagger},\quad \hat D(\alpha)\hat a\hat D(\alpha)^{-1}=\hat a-\alpha\hat I $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{DisplacementOperator[al, Infinity, "Ordering" -> "Normal"],
 DisplacementOperator[al, Infinity, "Ordering" -> "Antinormal"],
 BosonicExpSimplify[Exp[al ad - Conjugate[al] a] ** a ** Exp[-al ad + Conjugate[al] a]]}
```
Output: `{E^(al*SuperDagger[\[FormalA]1])**E^(-(\[FormalA]1*Conjugate[al]))/E^((al*Conjugate[al])/2), E^((al*Conjugate[al])/2)*E^(-(\[FormalA]1*Conjugate[al]))**E^(al*SuperDagger[\[FormalA]1]), \[FormalA]1 - al}`. Matches all three statements.

**Also stated in G-015 (Displacement: normal and antinormal factorizations, adjoint):**

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

<a id="displacement-adjoint-action-on-a"></a>
### Displacement: adjoint action on a
Appears in: X-001; P-008; G-018 (used in G-059, G-061, G-093)

Conjugating the annihilator by the displacement shifts it by the c-number $\alpha$.

$$ \hat D(\alpha)^\dagger \hat a\hat D(\alpha)=\hat a+\alpha\hat I $$

Conditions: the atlas notes that Berezin's direction $\hat U\hat a\hat U^{-1}$ reverses the sign of the parameter; with $\hat U=\hat D(\alpha)$ that direction gives $\hat a-\alpha$.
Audit: Correct, including the Berezin-direction remark. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
g = al ad - Conjugate[al] a;
{BosonicExpSimplify[Exp[-g] ** a ** Exp[g]], BosonicExpSimplify[Exp[g] ** a ** Exp[-g]]}
```
Output: `{\[FormalA]1 + al, \[FormalA]1 - al}`. Matches: the first is $\hat D^\dagger\hat a\hat D$, the second the reversed direction $\hat D\hat a\hat D^\dagger$.

**Also stated in P-008 (Displacement conjugation of a, q, p):**

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

**Also stated in G-018 (Displacement as a translation of the field):**

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

<a id="displacement-composition-law"></a>
### Displacement: composition law
Appears in: X-001; P-003; C-003; same law as C-059 (one-mode re-expression in terms of Im<f,g>); G-017, G-022 (phase of composed drives)

Two displacements compose to the displacement by the sum, up to a phase set by the symplectic area of the two amplitudes.

$$ \hat D(\alpha)\hat D(\beta)=e^{(\alpha\bar\beta-\bar\alpha\beta)/2}\hat D(\alpha+\beta) $$

Audit: Correct (with the opposite phase sign the scalar would not match the QF result). Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
lhs = BosonicExpSimplify[Exp[al ad - Conjugate[al] a] ** Exp[be ad - Conjugate[be] a]];
{lhs, Simplify[lhs - Exp[(al Conjugate[be] - Conjugate[al] be)/2] Exp[(al + be) ad - (Conjugate[al] + Conjugate[be]) a]]}
```
Output: `{E^(-(\[FormalA]1*Conjugate[al]) - \[FormalA]1*Conjugate[be] + (-(be*Conjugate[al]) + al*Conjugate[be])/2 + al*SuperDagger[\[FormalA]1] + be*SuperDagger[\[FormalA]1]), 0}`. Matches: the merged exponent carries the scalar $(\alpha\bar\beta-\bar\alpha\beta)/2$ and the difference from the atlas right side is zero.

**Also stated in P-003 (Weyl composition law and adjoint):**

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

**Also stated in C-003 (Weyl translation product and phase):**

Two phase-space translations compose into the translation by the summed vector times a phase set by the symplectic area between them; reversing the order reverses that phase, and the inverse is the translation by -z.

$$ \hat T(z)\hat T(z')=e^{-i\sigma(z,z')/(2\hbar)}\hat T(z+z'),\quad \hat T(z)\hat T(z')=e^{-i\sigma(z,z')/\hbar}\hat T(z')\hat T(z),\quad \hat T(z)^{-1}=\hat T(-z) $$

Audit: Correct. With alpha = (q+ip)/sqrt(2 hbar), alpha conj(beta) - conj(alpha) beta = -i sigma/hbar, which forces the phase. Checked by: QF exact (hbar = 1).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
d[x_] := Exp[x ad - Conjugate[x] a];
{BosonicExpSimplify[d[al] ** d[be]], BosonicExpSimplify[d[be] ** d[al]], BosonicExpSimplify[d[al] ** d[-al]],
 ComplexExpand[(al Conjugate[be] - be Conjugate[al])/2 /. {al -> (q1 + I p1)/Sqrt[2], be -> (q2 + I p2)/Sqrt[2]}]}
```
Output: `{E^(-(\[FormalA]1*Conjugate[al]) - \[FormalA]1*Conjugate[be] + (-(be*Conjugate[al]) + al*Conjugate[be])/2 + al*SuperDagger[\[FormalA]1] + be*SuperDagger[\[FormalA]1]), E^(-(\[FormalA]1*Conjugate[al]) - \[FormalA]1*Conjugate[be] + (be*Conjugate[al] - al*Conjugate[be])/2 + al*SuperDagger[\[FormalA]1] + be*SuperDagger[\[FormalA]1]), 1, I*(-1/2*(p2*q1) + (p1*q2)/2)}`. Matches: the scalar in the first product is (al conj(be) - be conj(al))/2, which the last entry shows is -i sigma(z,z')/2; the reversed product carries the opposite scalar, so their ratio is e^{-i sigma}; D(al) D(-al) = 1.

**Also stated in G-017 (Weyl composition law and coherent displacement phase):**

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

<a id="displacement-in-the-source-direction"></a>
### Displacement in the source direction
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

<a id="canonical-commutators-ladder-operators-and-the-multimode-displacement"></a>
### Canonical commutators, ladder operators and the multimode displacement
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

<a id="displacement-commutation-phase-group-commutator-and-polygon-phase"></a>
### Displacement commutation phase, group commutator and polygon phase
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

<a id="lattice-commutation-condition-and-critical-coherent-lattice"></a>
### Lattice commutation condition and critical coherent lattice
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

<a id="displacement-in-position-representation"></a>
### Displacement in position representation
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

<a id="displacement-conjugation-and-coordinate-factorization"></a>
### Displacement conjugation and coordinate factorization
Appears in: C-004

Conjugating by a translation shifts position and momentum by the classical displacement; the translation factors into a position-space kick and a momentum-generated shift times a scalar phase, and in position space it shifts the argument and multiplies by a plane wave.

$$ \hat T(z)\hat Q\hat T(z)^{-1}=\hat Q-q\hat I,\quad \hat T(z)=e^{-iq\cdot p/(2\hbar)}e^{ip\cdot \hat Q/\hbar}e^{-iq\cdot \hat P/\hbar},\quad (\hat T(z)u)(x)=e^{ip\cdot(x-q/2)/\hbar}u(x-q) $$

Audit: Correct. The position action follows by hand: e^{-iqP/hbar} u(x) = u(x-q) and e^{ipQ/hbar} multiplies by e^{ipx/hbar}, so with the phase e^{-iqp/(2 hbar)} one gets e^{ip(x-q/2)/hbar} u(x-q). Checked by: QF exact (hbar = 1) for the conjugation and the factorization; hand derivation: e^{ipx} e^{-iqp/2} u(x-q) for the position action.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2]; p = -I (a - ad)/Sqrt[2];
{BosonicExpSimplify[Exp[I (p0 q - q0 p)] ** q ** Exp[-I (p0 q - q0 p)]],
 BosonicExpSimplify[Exp[I (p0 q - q0 p)] ** p ** Exp[-I (p0 q - q0 p)]],
 BosonicExpSimplify[Exp[-I q0 p0/2] Exp[I p0 q] ** Exp[-I q0 p]]}
```
Output: `{-q0 + (\[FormalA]1 + SuperDagger[\[FormalA]1])/Sqrt[2], -p0 - (I*(\[FormalA]1 - SuperDagger[\[FormalA]1]))/Sqrt[2], E^(-((q0*(\[FormalA]1 - SuperDagger[\[FormalA]1]))/Sqrt[2]) + (I*p0*(\[FormalA]1 + SuperDagger[\[FormalA]1]))/Sqrt[2])}`. Matches: T Q T^{-1} = Q - q, T P T^{-1} = P - p, and e^{-iqp/2} e^{ipQ} e^{-iqP} returns exp(i(p Q - q P)) = T(z).

<a id="displacement-conjugation-of-a-field-and-of-n"></a>
### Displacement conjugation of a field and of N
Appears in: C-060

Conjugating by a Weyl operator shifts the field a(g) by <g,f> and the number operator by the linear field terms plus ||f||^2.

$$ \hat T(f)^*\hat a(g)\hat T(f)=\hat a(g)+\langle g,f\rangle\hat I,\quad \hat T(f)^*\hat N\hat T(f)=\hat N+\hat a(f)+\hat a^\dagger(f)+\|f\|^2\hat I $$

Audit: Correct (one mode: a(g) = conj(g) a, so a(g) + <g,f> is conj(g)(a + f)). Checked by: QF exact (one mode; several modes add mode by mode).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicExpSimplify[Exp[-al ad + Conjugate[al] a] ** a ** Exp[al ad - Conjugate[al] a]],
 BosonicExpSimplify[Exp[-al ad + Conjugate[al] a] ** (ad ** a) ** Exp[al ad - Conjugate[al] a]]}
```
Output: `{\[FormalA]1 + al, (Conjugate[al] + SuperDagger[\[FormalA]1])**(\[FormalA]1 + al)}`. Matches: D(-al) a D(al) = a + al and D(-al) N D(al) = (a^dag + conj(al))(a + al) = N + al a^dag + conj(al) a + |al|^2.

<a id="derivative-of-a-moving-displacement"></a>
### Derivative of a moving displacement
Appears in: C-061

The time derivative of a moving Weyl operator is the operator times the field of the velocity plus the scalar i Im<f, f'>.

$$ \frac d{dt}\hat T(f_t)=\hat T(f_t)\big[\hat a^\dagger(\dot f_t)-\hat a(\dot f_t)+i\Im\langle f_t,\dot f_t\rangle\hat I\big] $$

Audit: Correct. Checked by: QF exact (one mode): the exponent of D(-alpha) D(alpha + e beta) is linear in e with coefficient beta a^dag - conj(beta) a + (conj(alpha) beta - alpha conj(beta))/2, and (conj(alpha) beta - alpha conj(beta))/2 = i Im(conj(alpha) beta).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
d[x_, xc_] := Exp[x ad - xc a];
gen = BosonicExpSimplify[d[-al, -alc] ** d[al + e be, alc + e bec]] /. E^x_ :> x;
{gen, Simplify[D[gen, e] - (be ad - bec a + (alc be - al bec)/2)]}
```
Output: `{\[FormalA]1*alc + ((alc*be - al*bec)*e)/2 - \[FormalA]1*(alc + bec*e) - al*SuperDagger[\[FormalA]1] + (al + be*e)*SuperDagger[\[FormalA]1], 0}`. Matches: the difference is 0 (conjugates entered as independent symbols alc, bec).

<a id="vacuum-action-of-the-ordered-displacement"></a>
### Vacuum action of the ordered displacement
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

<a id="recoil-exponential-as-a-displacement-and-its-interaction-picture-form"></a>
### Recoil exponential as a displacement and its interaction-picture form
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

<a id="coherent-states-general"></a>
## Coherent states (general)

<a id="poisson-coherent-vectors"></a>
### Poisson (coherent) vectors
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

<a id="overlap-of-displaced-fiducial-states"></a>
### Overlap of displaced fiducial states
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

<a id="coherent-state-expansion-eigenvalue-equation-overlap-and-resolution"></a>
### Coherent state expansion, eigenvalue equation, overlap and resolution
Appears in: P-012; C-007; infinite-mode form Psi(f) in C-059; G-013, G-016 (vacuum displacement), G-021 (Poisson coefficients)

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

**Also stated in C-007 (Coherent state: vacuum action, eigenvalue, number expansion):**

A coherent state is the vacuum displaced by alpha; it is an eigenvector of the annihilator with eigenvalue alpha and has amplitudes alpha^n/sqrt(n!) times e^{-|alpha|^2/2} in the number basis.

$$ |\alpha\rangle=e^{-|\alpha|^2/2}e^{\alpha\cdot \hat a^\dagger}|0\rangle,\quad \hat a_j|\alpha\rangle=\alpha_j|\alpha\rangle,\quad |\alpha\rangle=e^{-|\alpha|^2/2}\sum_{\nu}\frac{\alpha^\nu}{\sqrt{\nu!}}|\nu\rangle $$

Audit: Correct. Checked by: QF exact for n = 0 to 5 and for the eigenvalue; hand derivation for general n: e^{alpha a^dag}|0> = sum_n alpha^n (a^dag)^n/n! |0> = sum_n alpha^n/sqrt(n!) |n>.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{Simplify[Table[BosonicMatrixElement[{n, 0}, DisplacementOperator[al]] - Exp[-Abs[al]^2/2] al^n/Sqrt[n!], {n, 0, 5}]],
 BosonicExpSimplify[a ** Exp[al ad]]}
```
Output: `{{0, 0, 0, 0, 0, 0}, E^(al*SuperDagger[\[FormalA]1])**(\[FormalA]1 + al)}`. Matches: the amplitude differences vanish, and a e^{al a^dag} = e^{al a^dag}(a + al), so a acting on e^{al a^dag}|0> returns al times the same vector.

**Also stated in G-013 (Coherent state as eigenvector and its Fock series):**

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

<a id="bargmann-realization-and-displacement-as-weighted-translation"></a>
### Bargmann realization and displacement as weighted translation
Appears in: P-013; C-011; G-023

In the Bargmann space $\hat a^\dagger$ is multiplication by $z$ and $\hat a$ is $\partial_z$, and a displacement acts as a weighted translation of the entire function.

$$ \hat a\leftrightarrow\partial_z,\quad \hat a^\dagger\leftrightarrow z,\quad \|f\|^2=\int\frac{d^2z}{\pi}e^{-|z|^2}|f(z)|^2,\quad [\hat D(\alpha)f](z)=e^{-|\alpha|^2/2+\alpha z}f(z-\bar\alpha) $$

Audit: Correct, including the boundary example $e^{z^2/2}$. Checked by: QF exact plus hand derivation: $\hat D f(\hat a^\dagger)|0\rangle=e^{-|\alpha|^2/2}e^{\alpha\hat a^\dagger}e^{-\bar\alpha\hat a}f(\hat a^\dagger)|0\rangle=e^{-|\alpha|^2/2}e^{\alpha\hat a^\dagger}f(\hat a^\dagger-\bar\alpha)|0\rangle$, which reads $e^{-|\alpha|^2/2+\alpha z}f(z-\bar\alpha)$. Exact WL computation: the norm series of $e^{z^2/2}$, $\sum_k(2k)!/(4^k(k!)^2)$, diverges (`SumConvergence` returns False).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{DisplacementOperator[al, Infinity, "Ordering" -> "Normal"], BosonicExpSimplify[Exp[-Conjugate[al] a] ** ad ** Exp[Conjugate[al] a]]}
```
Output: `{E^(al*SuperDagger[\[FormalA]1])**E^(-(\[FormalA]1*Conjugate[al]))/E^((al*Conjugate[al])/2), -Conjugate[al] + SuperDagger[\[FormalA]1]}`. Matches: the normal form of $\hat D$, and $e^{-\bar\alpha\hat a}\hat a^\dagger e^{\bar\alpha\hat a}=\hat a^\dagger-\bar\alpha$.

**Also stated in C-011 (Bargmann realization):**

In the Bargmann space of entire functions, creation acts as multiplication by zeta and annihilation as the derivative, so the linear and number exponentials act as a multiplication, a shift, and a dilation.

$$ \hat a^\dagger_j\mapsto\zeta_j,\quad \hat a_j\mapsto\partial_{\zeta_j},\quad \hat N\mapsto\zeta\cdot\partial_\zeta;\quad e^{v\cdot \hat a}f(\zeta)=f(\zeta+v),\ e^{u\cdot \hat a^\dagger}f=e^{u\cdot\zeta}f,\ e^{s\hat N}f(\zeta)=f(e^s\zeta) $$

Conditions: hbar = 1, analytic vectors.
Audit: Correct (Taylor shift, multiplication, Euler scaling). Checked by: hand derivation: e^{v d/dzeta} f(zeta) = sum_k v^k f^{(k)}(zeta)/k! = f(zeta + v), and zeta d/dzeta zeta^k = k zeta^k gives e^{sN} zeta^k = (e^s zeta)^k.
QF code (exact): none. QF has no function-space (Bargmann) realization; the operator statement behind the dilation, e^{sN} a e^{-sN} = e^{-s} a, is the C-094 call.

**Also stated in G-023 (Bargmann (entire-function) representation):**

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

<a id="vacuum-displacement-amplitude-and-coherent-overlap"></a>
### Vacuum displacement amplitude and coherent overlap
Appears in: C-008; same overlap in C-059 (as exp(-|f-g|^2/2 + i Im<f,g>)); C-010; G-019, G-104 (trace formula)

The vacuum expectation of a translation is a Gaussian in the displacement, and two coherent states overlap as a Gaussian in their separation times a symplectic phase.

$$ \langle0|\hat T(z)|0\rangle=e^{-|z|^2/(4\hbar)},\quad \langle\phi_z,\phi_{z'}\rangle=e^{i\sigma(z,z')/(2\hbar)}e^{-|z-z'|^2/(4\hbar)},\quad \hat T(z)\phi_{z'}=e^{-i\sigma(z,z')/(2\hbar)}\phi_{z+z'},\quad \langle\alpha|\beta\rangle=e^{-(|\alpha|^2+|\beta|^2)/2+\bar\alpha\cdot\beta} $$

Audit: Correct. |alpha|^2 = |z|^2/(2 hbar) and i Im(conj(alpha) beta) = i sigma(z,z')/(2 hbar) convert one form into the other; T(z) phi_{z'} follows from C-003 applied to T(z') acting on the vacuum. Checked by: QF exact (hbar = 1) plus hand derivation.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
vac = BosonicMatrixElement[{0, 0}, DisplacementOperator[al]];
ov = BosonicMatrixElement[{al, be}, 1, "Basis" -> "Coherent"];
toZ = {al -> (q1 + I p1)/Sqrt[2], be -> (q2 + I p2)/Sqrt[2]};
{vac, ov, Simplify[ComplexExpand[{vac, ov} /. E^x_ :> x /. toZ] -
   {-(q1^2 + p1^2)/4, I (q1 p2 - p1 q2)/2 - ((q1 - q2)^2 + (p1 - p2)^2)/4}]}
```
Output: `{E^(-1/2*Abs[al]^2), E^((-Abs[al]^2 - Abs[be]^2)/2 + be*Conjugate[al]), {0, 0}}`. Matches: both exponents equal the phase-space forms (last entry zero).

**Also stated in C-010 (Coherent resolution of identity, trace and Hilbert-Schmidt norm):**

Coherent projectors integrated over phase space with weight (2 pi hbar)^{-n} give the identity, which turns traces and Hilbert-Schmidt norms into phase-space integrals of coherent expectation values.

$$ (2\pi\hbar)^{-n}\int|\phi_z\rangle\langle\phi_z|\,dz=\hat I,\quad \operatorname{Tr}\hat B=(2\pi\hbar)^{-n}\int\langle\phi_z,\hat B\phi_z\rangle dz,\quad \|\hat B\|_{\rm HS}^2=(2\pi\hbar)^{-n}\int\|\hat B\phi_z\|^2dz $$

Conditions: trace class, resp. Hilbert-Schmidt.
Audit: Correct. dq dp = 2 hbar d^2 alpha turns the prefactor into 1/pi per mode; the number-basis element <m|...|n> has an angular integral that vanishes for m != n and a radial part (2/n!) times the integral of r^{2n+1} e^{-r^2}, which exact integration gives as 1. Checked by: exact WL computation plus hand derivation.
QF code (exact): none. QF has no phase-space integration; `BosonicMatrixElement[{al, be}, 1, "Basis" -> "Coherent"]` gives only the overlaps in the integrand.

**Also stated in G-019 (Coherent overlap and resolution of identity):**

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

<a id="harmonic-rotation-of-displacements-and-coherent-states"></a>
### Harmonic rotation of displacements and coherent states
Appears in: C-009

The oscillator propagator rotates phase space rigidly: it conjugates a translation into the translation by the rotated vector and sends a coherent state to the coherent state labelled e^{-it} alpha, times the zero-point phase.

$$ \hat U_t\hat T(z)\hat U_t^*=\hat T(F_tz),\quad F_t=\begin{pmatrix}\cos t&\sin t\\-\sin t&\cos t\end{pmatrix},\quad \hat U_t|\alpha\rangle=e^{-int/2}|e^{-it}\alpha\rangle $$

Conditions: U_t = exp(-i t (N + n/2)).
Audit: Correct (e^{-it} alpha corresponds to F_t z). QF gives e^{-itN} a e^{itN} = e^{it} a, so U_t (alpha a^dag - conj(alpha) a) U_t^* is the generator of D(e^{-it} alpha); with U_t|0> = e^{-int/2}|0> this gives the state formula. Checked by: QF exact (one mode) plus hand derivation.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
rot = BosonicExpSimplify[Exp[-I t ad ** a] ** a ** Exp[I t ad ** a]];
{rot, Simplify[rot/a - Exp[I t], t > 0]}
```
Output: `{\[FormalA]1*(Cosh[Sqrt[-t^2]] + (I*t*Sinh[Sqrt[-t^2]])/Sqrt[-t^2]), 0}`. Matches: Cosh[Sqrt[-t^2]] + I t Sinh[Sqrt[-t^2]]/Sqrt[-t^2] is cos t + i sin t = e^{it} (second entry zero).

<a id="canonical-oscillator-realization-of-the-ccr"></a>
### Canonical oscillator realization of the CCR
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

<a id="quadrature-commutator-and-number-state-variance"></a>
### Quadrature commutator and number-state variance
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

<a id="no-normalizable-right-eigenvector-of-the-creation-operator"></a>
### No normalizable right eigenvector of the creation operator
Appears in: G-014

The creation operator has no normalizable right eigenvector: the eigenvalue equation forces every Fock component to vanish.

$$ \hat a^\dagger|\psi\rangle=\lambda|\psi\rangle\ \Rightarrow\ |\psi\rangle=0 $$

Audit: Correct. Checked by: hand derivation: the components read 0 = lambda c_0 and sqrt(m) c_{m-1} = lambda c_m; for lambda != 0 this gives c_0 = 0 and then every c_m = 0, for lambda = 0 it gives every c_{m-1} = 0.
QF code (exact): none. a nonexistence statement; there is nothing to compute.

<a id="coherent-state-photon-statistics"></a>
### Coherent-state photon statistics
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

<a id="photon-subtraction-and-addition"></a>
### Photon subtraction and addition
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

<a id="single-mode-squeezing"></a>
## Single-mode squeezing

<a id="single-mode-squeeze-normal-ordered-factorization"></a>
### Single-mode squeeze: normal-ordered factorization
Appears in: X-002 (source records cited: P-028 to P-031, G-051 to G-061, B-051, C-048, C-049); P-031; G-054; C-048 (two-mode version: Two-mode squeezing)

The squeeze operator factors into a pair-creation exponential, a power of the number operator, and a pair-annihilation exponential, with $\tanh r$ in the outer factors.

$$ \hat S(\xi)=\exp\!\left[\tfrac12\left(\bar\xi \hat a^2-\xi(\hat a^\dagger)^2\right)\right]=e^{-\frac12e^{i\theta}\tanh r\,(\hat a^\dagger)^2}e^{-\log\cosh r\,(\hat N+\tfrac12\hat I)}e^{\frac12e^{-i\theta}\tanh r\,\hat a^2} $$

Conditions: $\xi=re^{i\theta}$, $r\ge0$.
Audit: Correct, including the $\tfrac12$ in the number exponent. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{SqueezeOperator[xi, Infinity, "Ordering" -> "Normal"], FullSimplify[BosonicExpOrder[Exp[(r Exp[-I th] a ** a - r Exp[I th] ad ** ad)/2]], r > 0 && Element[th, Reals]]}
```
Output: `{E^(-1/2*(E^(I*Arg[xi])*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]1]*Tanh[Abs[xi]]))**Sech[Abs[xi]]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^((\[FormalA]1**\[FormalA]1*Tanh[Abs[xi]])/(2*E^(I*Arg[xi]))), E^(-1/2*(E^(I*th)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]*Tanh[r]))**Sech[r]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^((\[FormalA]1**\[FormalA]1*Tanh[r])/(2*E^(I*th)))}`. Matches: $\operatorname{sech}r^{\,\hat N+1/2}=e^{-\log\cosh r(\hat N+\frac12)}$, so both the named operator and the ordered bare exponential are the atlas form.

**Also stated in P-031 (Single-mode squeezing: normal-ordered factorization):**

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

**Also stated in G-054 (Squeeze operator: normal-ordered Gauss factorization):**

The squeeze operator factors into a normal-ordered product: an exponential of a^dagger^2, a power of sech r in the number operator, and an exponential of a^2.

$$ \hat S(\xi)=e^{-\tau \hat a^{\dagger2}/2}e^{-\ln(\cosh r)(\hat n+\tfrac12\hat I)}e^{\tau^*\hat a^2/2},\qquad\tau=e^{i\theta}\tanh r. $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{SqueezeOperator[xi, Infinity, "Ordering" -> "Normal"], BosonicExpOrder[Exp[c1 a ** a + c2 ad ** ad]]}
```
Output: `{E^(-1/2*(E^(I*Arg[xi])*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]1]*Tanh[Abs[xi]]))**Sech[Abs[xi]]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^((\[FormalA]1**\[FormalA]1*Tanh[Abs[xi]])/(2*E^(I*Arg[xi]))), E^((c2*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]*Tan[2*Sqrt[c1*c2]])/(2*Sqrt[c1*c2]))**Sec[2*Sqrt[c1*c2]]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^((c1*GeneralizedPower[NonCommutativeMultiply, \[FormalA]1, 2]*Tan[2*Sqrt[c1*c2]])/(2*Sqrt[c1*c2]))}`. Matches: Sech[r]^(1/2 + n) = e^{-ln(cosh r)(n + 1/2)} and tau = e^{i Arg xi} Tanh[Abs xi]; the general BosonicExpOrder form reduces to it with c1 = xi^*/2, c2 = -xi/2, sqrt(c1 c2) = i r/2, Tan[i r]/(i r) = Tanh[r]/r.

**Also stated in C-048 (Noncompact Gauss factorization of squeezers):**

A squeeze exponential factors into a pair-creation exponential, a number exponential with log cosh r, and a pair-annihilation exponential with tanh r; the same holds for two-mode pairs.

$$ e^{\xi\hat K_+-\bar\xi\hat K_-}=e^{\zeta\hat K_+}e^{-2\log(\cosh r)\hat K_0}e^{-\bar\zeta\hat K_-},\ \zeta=\frac\xi r\tanh r;\qquad \hat S(\xi)=e^{\zeta\hat a^{\dagger2}/2}e^{-\log(\cosh r)(\hat N+\frac12)}e^{-\bar\zeta\hat a^2/2} $$

and the two-mode form with a1^dag a2^dag, a1 a2, N1 + N2 + 1.

Audit: Correct. Checked by: exact WL computation (2-by-2 representation of su(1,1), symbolic r and theta) and QF exact (one mode with real r and complex xi, two modes; here the atlas's S(xi) is entered as an explicit exponential, so no sign flip).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
{BosonicExpOrder[Exp[(r ad ** ad - r a ** a)/2], r > 0],
 BosonicExpOrder[Exp[(xi ad ** ad - Conjugate[xi] a ** a)/2]],
 BosonicExpOrder[Exp[xi a1d ** a2d - Conjugate[xi] a1 ** a2]],
 Simplify[{Tan[Sqrt[-xi Conjugate[xi]]]/Sqrt[-xi Conjugate[xi]], Sec[Sqrt[-xi Conjugate[xi]]]} /. xi -> r E^(I th), r > 0 && Element[th, Reals]]}
```
Output: `{E^((GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]*Tanh[r])/2)**Sech[r]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^(-1/2*(GeneralizedPower[NonCommutativeMultiply, \[FormalA]1, 2]*Tanh[r])), <one-mode complex xi form with Tan[Sqrt[-(xi*Conjugate[xi])]]>, E^((xi*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]2]*Tan[Sqrt[-(xi*Conjugate[xi])]])/Sqrt[-(xi*Conjugate[xi])])**Sec[Sqrt[-(xi*Conjugate[xi])]]^(SuperDagger[\[FormalA]1]**\[FormalA]1 + SuperDagger[\[FormalA]2]**\[FormalA]2)**E^((Sqrt[-(xi*Conjugate[xi])]*\[FormalA]1**\[FormalA]2*Tan[Sqrt[-(xi*Conjugate[xi])]])/xi)*Sec[Sqrt[-(xi*Conjugate[xi])]], {Tanh[r]/r, Sech[r]}}`. Matches: with Tan[Sqrt[-|xi|^2]]/Sqrt[-|xi|^2] = Tanh[r]/r and Sec[Sqrt[-|xi|^2]] = Sech[r] (last entry), the one-mode result is e^{zeta a^dag^2/2} Sech[r]^(N+1/2) e^{-conj(zeta) a^2/2} and the two-mode result e^{zeta a1^dag a2^dag} Sech[r]^(N1+N2+1) e^{-conj(zeta) a1 a2}, zeta = xi tanh(r)/r.

<a id="single-mode-squeeze-adjoint-action-on-a"></a>
### Single-mode squeeze: adjoint action on a
Appears in: X-002; G-052

Squeezing mixes the annihilator with the creator, with $\cosh r$ and $\sinh r$ weights.

$$ \hat S^\dagger \hat a\hat S=\hat a\cosh r-e^{i\theta}\hat a^\dagger\sinh r $$

Conditions: the atlas notes that $\hat S^{-1}$ flips the sign of every squeeze parameter.
Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
g = (r Exp[-I th] a ** a - r Exp[I th] ad ** ad)/2;
Simplify[BosonicExpSimplify[Exp[-g] ** a ** Exp[g]], r > 0]
```
Output: `\[FormalA]1*Cosh[r] - E^(I*th)*Sinh[r]*SuperDagger[\[FormalA]1]`. Matches.

**Also stated in G-052 (Single-mode Bogoliubov transformation):**

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

<a id="single-mode-squeezed-vacuum"></a>
### Single-mode squeezed vacuum
Appears in: X-002; G-056 (also G-062 at $\alpha=0$, G-063 overlap)

The squeezed vacuum is a normalized Gaussian superposition of even photon numbers.

$$ \hat S|0\rangle=(\cosh r)^{-1/2}e^{-\frac12e^{i\theta}\tanh r\,(\hat a^\dagger)^2}|0\rangle $$

Audit: Correct; it is X-006 with $Z=e^{i\theta}\tanh r$, since $(1-\tanh^2r)^{1/4}=(\cosh r)^{-1/2}$. Checked by: QF exact, and hand derivation: in the QF normal form above the annihilation factor acts as $\hat I$ on $|0\rangle$ and $\operatorname{sech}r^{\,\hat N+1/2}|0\rangle=(\cosh r)^{-1/2}|0\rangle$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
tau = Exp[I th] Tanh[r];
formula[m_] := If[EvenQ[m], Cosh[r]^(-1/2) (-tau/2)^(m/2) Sqrt[m!]/(m/2)!, 0];
{Table[FullSimplify[BosonicMatrixElement[{m, 0}, SqueezeOperator[r Exp[I th], Infinity]] - formula[m], r > 0 && Element[th, Reals]], {m, 0, 10}], FullSimplify[Sum[Abs[formula[2 k]]^2, {k, 0, Infinity}], r > 0 && Element[th, Reals]]}
```
Output: `{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 0}`. Matches: QF's closed-form amplitudes $\langle m|\hat S|0\rangle$ equal the formula's for $m=0,\dots,10$ with $r,\theta$ symbolic, and the formula's state has norm exactly 1 (the second entry is norm minus 1, summed in closed form over all $m$).

**Also stated in G-056 (Squeezed vacuum: nullifier and Fock expansion):**

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

<a id="one-mode-squeezing-implementer-and-kernel"></a>
### One-mode squeezing implementer and kernel
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

<a id="single-mode-su-1-1-realization-and-casimir"></a>
### Single-mode su(1,1) realization and Casimir
Appears in: P-031, P2-012 (parity sectors $k=1/4,3/4$); G-053

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

**Also stated in G-053 (Quadratic SU(1,1) algebra of one mode):**

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

<a id="multimode-squeeze-operator-and-bounded-squeeze-coordinate"></a>
### Multimode squeeze operator and bounded squeeze coordinate
Appears in: C-025

The multimode squeeze operator is the exponential of a quadratic in creators minus its adjoint; its inverse flips the sign of B, and B and the bounded coordinate W are related by tanh and the arctanh series.

$$ \hat S(B)=\exp\tfrac12(\hat a^\dagger\cdot B\hat a^\dagger-\hat a\cdot B^*\hat a),\quad \hat S(B)^{-1}=\hat S(-B),\quad B=U\operatorname{arctanh}|W|=W\sum_k\frac{(W^*W)^k}{2k+1},\quad W=U\tanh|B| $$

Conditions: B complex symmetric, W^*W < I.
Audit: Correct (the scalar series sums to ArcTanh; for matrices apply it on the eigenbasis of W^*W, where W sum_k (W^*W)^k/(2k+1) = U arctanh|W|). Convention: QF's `SqueezeOperator[xi]` is this S(-xi) for one mode. Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
g = (xi ad ** ad - Conjugate[xi] a ** a)/2;
{SqueezeOperator[xi, Infinity], BosonicExpSimplify[Exp[g] ** Exp[-g]] /. E^x_ :> E^Expand[x],
 Sum[x^(2 k + 1)/(2 k + 1), {k, 0, Infinity}]}
```
Output: `{E^((Conjugate[xi]*\[FormalA]1**\[FormalA]1 - xi*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]1])/2), 1, ArcTanh[x]}`. Matches: QF's squeezer is exp((conj(xi) a^2 - xi a^dag^2)/2), the atlas S(-xi); S S^{-1} = 1; the series is ArcTanh.

<a id="squeeze-conjugation-in-block-form"></a>
### Squeeze conjugation in block form
Appears in: C-026; one-mode case also in C-049 (first display)

Squeezing acts on the column (a, a^dag) by a block matrix built from W; for one mode it gives a cosh r - e^{i theta} a^dag sinh r.

$$ \hat S(B)\binom{\hat a}{\hat a^\dagger}\hat S(-B)=\begin{pmatrix}(I-WW^*)^{-1/2}&-W(I-W^*W)^{-1/2}\\-(I-W^*W)^{-1/2}W^*&(I-W^*W)^{-1/2}\end{pmatrix}\binom{\hat a}{\hat a^\dagger},\qquad \hat S\hat a\hat S^{-1}=\hat a\cosh r-e^{i\theta}\hat a^\dagger\sinh r $$

Audit: Correct. Checked by: QF exact (one mode: symbolic r, theta, with the atlas S(B) entered directly; two modes: QF gives the commutators of the generator with (a1, a2, a1^dag, a2^dag) exactly, the conjugation is the exponential of that 4-by-4 matrix by C-094, and it equals the block formula for a complex symmetric B at 40 digits).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
one = BosonicExpSimplify[Exp[(r E^(I th) ad ** ad - r E^(-I th) a ** a)/2] ** a **
    Exp[-(r E^(I th) ad ** ad - r E^(-I th) a ** a)/2]];
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
b = {{1/4 + I/10, 3/20 - I/20}, {3/20 - I/20, -1/10 + I/5}};
cr = {a1d, a2d}; an = {a1, a2}; ops = Join[an, cr];
g = (Sum[b[[i, j]] cr[[i]] ** cr[[j]], {i, 2}, {j, 2}] - Sum[Conjugate[b[[i, j]]] an[[i]] ** an[[j]], {i, 2}, {j, 2}])/2;
adg = Table[Coefficient[BosonicNormalOrder[g ** x - x ** g, ops], y], {x, ops}, {y, ops}];
mb = MatrixPower[N[ConjugateTranspose[b] . b, 40], 1/2];
w = N[b, 40] . Inverse[mb] . MatrixFunction[Tanh, mb];
blk = ArrayFlatten[{{MatrixPower[IdentityMatrix[2] - w . ConjugateTranspose[w], -1/2], -w . MatrixPower[IdentityMatrix[2] - ConjugateTranspose[w] . w, -1/2]},
    {-MatrixPower[IdentityMatrix[2] - ConjugateTranspose[w] . w, -1/2] . ConjugateTranspose[w], MatrixPower[IdentityMatrix[2] - ConjugateTranspose[w] . w, -1/2]}}];
{one, adg, Chop[Max[Abs[MatrixExp[N[adg, 40]] - blk]], 10^-30]}
```
Output: `{\[FormalA]1*Cosh[Sqrt[r^2]] - (E^(I*th)*r*Sinh[Sqrt[r^2]]*SuperDagger[\[FormalA]1])/Sqrt[r^2], {{0, 0, -1/4 - I/10, -3/20 + I/20}, {0, 0, -3/20 + I/20, 1/10 - I/5}, {-1/4 + I/10, -3/20 - I/20, 0, 0}, {-3/20 - I/20, 1/10 + I/5, 0, 0}}, 0}`. Matches: one mode gives a cosh r - e^{i th} sinh r a^dag; the commutator matrix is {{0, -B}, {-conj(B), 0}}; the block residual is 0.

<a id="squeezed-number-operator"></a>
### Squeezed number operator
Appears in: C-049

Squeezing maps the annihilator to (a - omega a^dag)/sqrt(1 - |omega|^2) and rotates K0 into a hyperbolic combination of K0, K+ and K-.

$$ \hat S\hat a\hat S^{-1}=\frac{\hat a-\omega\hat a^\dagger}{\sqrt{1-|\omega|^2}},\quad \hat S\hat K_0\hat S^{-1}=\cosh(2r)\hat K_0-\tfrac12\sinh(2r)(e^{i\theta}\hat K_++e^{-i\theta}\hat K_-) $$

Audit: Correct. Checked by: QF exact (symbolic r > 0 and theta, atlas S(xi) entered directly).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
g = (r E^(I th) ad ** ad - r E^(-I th) a ** a)/2;
k0 = (ad ** a + 1/2)/2; kp = ad ** ad/2; km = a ** a/2;
sa = BosonicExpSimplify[Exp[g] ** a ** Exp[-g]];
sk = BosonicExpSimplify[Exp[g] ** k0 ** Exp[-g]];
{sa, sk, Simplify[BosonicNormalOrder[sk - (Cosh[2 r] k0 - Sinh[2 r]/2 (E^(I th) kp + E^(-I th) km)), {a, ad}], r > 0],
 Simplify[sa - (a - E^(I th) Tanh[r] ad)/Sqrt[1 - Tanh[r]^2], r > 0]}
```
Output: `{\[FormalA]1*Cosh[Sqrt[r^2]] - (E^(I*th)*r*Sinh[Sqrt[r^2]]*SuperDagger[\[FormalA]1])/Sqrt[r^2], (Cosh[2*Sqrt[r^2]]*(1 + 2*SuperDagger[\[FormalA]1]**\[FormalA]1) - (r*(E^((2*I)*th)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2] + \[FormalA]1**\[FormalA]1)*Sinh[2*Sqrt[r^2]])/(E^(I*th)*Sqrt[r^2]))/4, 0, 0}`. Matches: both differences are 0.

<a id="single-mode-squeeze-operator"></a>
### Single-mode squeeze operator
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

<a id="squeeze-operator-reversed-antinormal-gauss-factorization"></a>
### Squeeze operator: reversed (antinormal) Gauss factorization
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

<a id="rotated-quadrature-variance-of-squeezed-vacuum"></a>
### Rotated-quadrature variance of squeezed vacuum
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

<a id="quadrature-covariance-of-squeezed-vacuum"></a>
### Quadrature covariance of squeezed vacuum
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

<a id="displaced-squeezed-state-low-moments"></a>
### Displaced squeezed state: low moments
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

<a id="interchanging-displacement-and-squeezing"></a>
### Interchanging displacement and squeezing
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

<a id="displaced-squeezed-state-as-a-bogoliubov-eigenstate"></a>
### Displaced squeezed state as a Bogoliubov eigenstate
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

<a id="number-variance-of-a-displaced-squeezed-vacuum"></a>
### Number variance of a displaced squeezed vacuum
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

<a id="squared-field-quadratures-and-their-commutators"></a>
### Squared-field quadratures and their commutators
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

<a id="two-mode-squeezing"></a>
## Two-mode squeezing

<a id="two-mode-squeeze-normal-ordered-factorization"></a>
### Two-mode squeeze: normal-ordered factorization
Appears in: X-003 (source records cited: P-032, G-068 to G-075, P-029); P-032, P2-019; G-070, G-071 (reversed form)

The two-mode squeezer factors the same way as the single-mode one, with the pair $\hat a^\dagger\hat b^\dagger$ in place of $(\hat a^\dagger)^2/2$.

$$ \hat S_2(\xi)=e^{\bar\xi \hat a\hat b-\xi \hat a^\dagger \hat b^\dagger}=e^{-e^{i\theta}\tanh r\,\hat a^\dagger \hat b^\dagger}e^{-\log\cosh r\,(\hat N_a+\hat N_b+\hat I)}e^{e^{-i\theta}\tanh r\,\hat a\hat b} $$

Conditions: closed algebra $\hat K_+=\hat a^\dagger\hat b^\dagger$, $\hat K_-=\hat a\hat b$, $\hat K_0=(\hat N_a+\hat N_b+\hat I)/2$; $\hat N_a-\hat N_b$ conserved; P-032's convention is $\xi\mapsto-\xi$.
Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
FullSimplify[BosonicExpOrder[Exp[r Exp[-I th] a1 ** a2 - r Exp[I th] ad1 ** ad2]], r > 0 && Element[th, Reals]]
```
Output: `E^(-(E^(I*th)*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]2]*Tanh[r]))**Sech[r]^(SuperDagger[\[FormalA]1]**\[FormalA]1 + SuperDagger[\[FormalA]2]**\[FormalA]2)**E^((\[FormalA]1**\[FormalA]2*Tanh[r])/E^(I*th))*Sech[r]`. Matches: the trailing `Sech[r]` is the $\hat I$ in the middle exponent.

**Also stated in P-032 (Two-mode squeezing: normal-ordered factorization and pair state):**

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

**Also stated in G-070 (Two-mode squeeze: normal and reversed Gauss factorizations):**

The two-mode squeeze factors in normal order and in reversed order, the middle factor being a power of cosh r in the total number plus one.

$$ \hat S_2(\xi)=e^{-\tau \hat a^\dagger \hat b^\dagger}e^{-\ln(\cosh r)(\hat n_a+\hat n_b+\hat I)}e^{\tau^*\hat a\hat b}=e^{\tau^*\hat a\hat b}e^{+\ln(\cosh r)(\hat n_a+\hat n_b+\hat I)}e^{-\tau \hat a^\dagger \hat b^\dagger}. $$

Audit: Correct. Checked by: QF exact (normal and antinormal BosonicExpOrder).
QF code (exact):
```wl
{a, ad, b, bd} = FieldVariables[{1, 2}];
{BosonicExpOrder[Exp[c1 a ** b + c2 ad ** bd]], BosonicExpOrder[Exp[c1 a ** b + c2 ad ** bd], "Ordering" -> "Antinormal"]}
```
Output: `{E^((c2*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]2]*Tan[Sqrt[c1*c2]])/Sqrt[c1*c2])**Sec[Sqrt[c1*c2]]^(SuperDagger[\[FormalA]1]**\[FormalA]1 + SuperDagger[\[FormalA]2]**\[FormalA]2)**E^((c1*\[FormalA]1**\[FormalA]2*Tan[Sqrt[c1*c2]])/Sqrt[c1*c2])*Sec[Sqrt[c1*c2]], E^((c1*\[FormalA]1**\[FormalA]2*Tan[Sqrt[c1*c2]])/Sqrt[c1*c2])**Sec[Sqrt[c1*c2]]^(-\[FormalA]1**SuperDagger[\[FormalA]1] - \[FormalA]2**SuperDagger[\[FormalA]2])**E^((c2*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]2]*Tan[Sqrt[c1*c2]])/Sqrt[c1*c2])*Sec[Sqrt[c1*c2]]}`. Matches after c1 = xi^*, c2 = -xi (sqrt(c1 c2) = i r, Tan[i r]/(i r) = Tanh[r]/r, Sec[i r] = Sech[r]); in the antinormal form a a^dagger + b b^dagger = n_a + n_b + 2 turns Sec^(-(...)) Sec into Cosh^(n_a + n_b + 1).

<a id="two-mode-squeezed-vacuum"></a>
### Two-mode squeezed vacuum
Appears in: X-003; G-072

The two-mode squeezed vacuum is a sum of twin-photon-number states with geometric weights.

$$ \hat S_2|0,0\rangle=\operatorname{sech}r\sum_{n\ge0}(-e^{i\theta}\tanh r)^n|n,n\rangle $$

Audit: Correct. Checked by: hand derivation: in the QF normal form above the annihilation factor and $\operatorname{sech}r^{\,\hat N_a+\hat N_b}$ act as $\hat I$ on $|0,0\rangle$, and $(\hat a^\dagger\hat b^\dagger)^n|0,0\rangle=n!\,|n,n\rangle$; and exact WL computation: in the Bargmann realization ($\hat a_i\to\partial_{z_i}$, $\hat a_i^\dagger\to z_i$) the function $\psi(s)=\operatorname{sech}(sr)\exp(-e^{i\theta}\tanh(sr)\,z_1z_2)$ satisfies $\partial_s\psi=(\bar\xi\partial_1\partial_2-\xi z_1z_2)\psi$ with $\psi(0)=1$ exactly, so it is $\hat S_2(s\xi)|0,0\rangle$.
QF code (exact): none. `BosonicMatrixElement` takes multimode Fock indices only for polynomials; `BosonicMatrixElement[{{2, 2}, {0, 0}}, Exp[x a1 ** a2 + y ad1 ** ad2]]` returns unevaluated.

**Also stated in G-072 (Twin-beam vacuum: Fock expansion and nullifiers):**

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

<a id="two-mode-su-1-1-realization-and-casimir"></a>
### Two-mode su(1,1) realization and Casimir
Appears in: P-032, P2-019; G-069

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

**Also stated in G-069 (Two-mode SU(1,1) algebra and conserved number difference):**

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

<a id="multimode-paired-vacuum-normalization"></a>
### Multimode paired vacuum normalization
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

<a id="multimode-paired-vacuum-overlap"></a>
### Multimode paired vacuum overlap
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

<a id="bosonic-pair-number-generating-determinant-trace-log-sign-corrected"></a>
### Bosonic pair-number generating determinant (trace-log sign corrected)
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

<a id="bosonic-realizations-of-su-1-1"></a>
### Bosonic realizations of SU(1,1)
Appears in: C-046

Squares of one mode and pair operators of two modes realize su(1,1); the one-mode Casimir is -3/16 and the two-mode Casimir is ((N1 - N2)^2 - 1)/4.

$$ \hat K_+=\tfrac12\hat a^{\dagger2},\ \hat K_-=\tfrac12\hat a^2,\ \hat K_0=\tfrac12(\hat N+\tfrac12),\ \hat C=-\tfrac3{16};\qquad \hat K_+=\hat a_1^\dagger\hat a_2^\dagger,\ \hat K_-=\hat a_1\hat a_2,\ \hat K_0=\tfrac12(\hat N_1+\hat N_2+1),\ \hat C=\tfrac14[(\hat N_1-\hat N_2)^2-1] $$

Audit: Correct (k = 1/4, 3/4 and k = (1+|nu|)/2 reproduce k(k-1)). Checked by: QF exact.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
k0 = (ad ** a + 1/2)/2; kp = ad ** ad/2; km = a ** a/2;
one = {BosonicNormalOrder[k0 ** k0 - (kp ** km + km ** kp)/2, {a, ad}],
   BosonicNormalOrder[km ** kp - kp ** km - 2 k0, {a, ad}], BosonicNormalOrder[k0 ** kp - kp ** k0 - kp, {a, ad}]};
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
n1 = a1d ** a1; n2 = a2d ** a2;
kk0 = (n1 + n2 + 1)/2; kkp = a1d ** a2d; kkm = a1 ** a2;
two = BosonicNormalOrder[kk0 ** kk0 - (kkp ** kkm + kkm ** kkp)/2 - ((n1 - n2) ** (n1 - n2) - 1)/4, {a1, a1d, a2, a2d}];
{one, two}
```
Output: `{{-3/16, 0, 0}, 0}`. Matches: the one-mode Casimir is -3/16 and both one-mode commutation relations give 0; the two-mode Casimir minus ((N1-N2)^2 - 1)/4 normal-orders to 0.

<a id="two-mode-squeeze-operator-and-bogoliubov-action"></a>
### Two-mode squeeze operator and Bogoliubov action
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

<a id="two-mode-squeeze-as-two-opposite-single-mode-squeezes-of-the-sum-and-difference-modes"></a>
### Two-mode squeeze as two opposite single-mode squeezes of the sum and difference modes
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

<a id="twin-beam-number-statistics-collective-squeezing-and-the-intermode-cauchy-schwarz-test"></a>
### Twin-beam number statistics, collective squeezing, and the intermode Cauchy-Schwarz test
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

<a id="pair-coherent-states"></a>
### Pair coherent states
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

<a id="low-gain-down-conversion-to-second-order"></a>
### Low-gain down-conversion to second order
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

<a id="su-2-and-passive-mode-mixing"></a>
## SU(2) and passive mode mixing

<a id="passive-transformation-of-normal-symbols"></a>
### Passive transformation of normal symbols
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

<a id="exponential-of-an-su-2-matrix"></a>
### Exponential of an su(2) matrix
Appears in: C-039

A 2-by-2 su(2) matrix squares to minus its determinant, so its exponential is a cosine plus a sine term.

$$ e^{tA}=I_2\cos t+A\sin t\ (\det A=1),\qquad e^{it\mathbf v\cdot\boldsymbol\sigma}=\cos(t|\mathbf v|)I_2+i\sin(t|\mathbf v|)\mathbf n\cdot\boldsymbol\sigma $$

Audit: Correct (A^2 = -I). Checked by: exact WL computation (2-by-2 matrices).
QF code (exact): none. A 2-by-2 matrix identity; QF's second-quantization layer does not apply.

<a id="second-quantization-of-a-one-particle-unitary"></a>
### Second quantization of a one-particle unitary
Appears in: C-058

Second quantization turns the one-particle unitary e^{-ith} into the exponential of dGamma(h); a number-conserving quadratic exponential acts on creators by the one-particle matrix and maps Psi(f) to Psi(Uf).

$$ \widehat\Gamma(e^{-it\hat h})=e^{-it\,\widehat{d\Gamma}(\hat h)},\quad \widehat{d\Gamma}(\hat I)=\hat N,\quad \widehat\Gamma(\hat U)\Psi(f)=\Psi(\hat Uf) $$

Audit: Correct. Checked by: QF exact (two modes, beam-splitter generator, symbolic theta and phi): e^G a_j^dag e^{-G} = sum_k (e^g)_{kj} a_k^dag with g the one-particle matrix; since G annihilates the vacuum, e^G Psi(f) = Psi(e^g f) follows by hand.
QF code (exact):
```wl
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
g = th (Exp[I ph] a2d ** a1 - Exp[-I ph] a1d ** a2);
out = {BosonicExpSimplify[Exp[g] ** a1d ** Exp[-g]], BosonicExpSimplify[Exp[g] ** a2d ** Exp[-g]]};
u = MatrixExp[{{0, -th Exp[-I ph]}, {th Exp[I ph], 0}}];
{out, Simplify[Table[Coefficient[out[[j]], {a1d, a2d}[[k]]] - u[[k, j]], {j, 2}, {k, 2}], th > 0 && Element[ph, Reals]]}
```
Output: `{{Cosh[Sqrt[-th^2]]*SuperDagger[\[FormalA]1] + (E^(I*ph)*th*Sinh[Sqrt[-th^2]]*SuperDagger[\[FormalA]2])/Sqrt[-th^2], (Sqrt[-th^2]*Sinh[Sqrt[-th^2]]*SuperDagger[\[FormalA]1])/(E^(I*ph)*th) + Cosh[Sqrt[-th^2]]*SuperDagger[\[FormalA]2]}, {{0, 0}, {0, 0}}}`. Matches: all four coefficient differences are 0.

<a id="atomic-transition-algebra-and-jaynes-cummings-excitation-number"></a>
### Atomic transition algebra and Jaynes-Cummings excitation number
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

<a id="detuned-jaynes-cummings-doublet-propagator"></a>
### Detuned Jaynes-Cummings doublet propagator
Appears in: G-039

Detuning turns each doublet into a spin-1/2 precession at the generalized Rabi frequency Omega_m, because the doubled doublet Hamiltonian squares to Omega_m^2.

$$ \hat h_m/\hbar=\frac12\begin{pmatrix}\Delta&2\lambda\sqrt{m+1}\\2\lambda\sqrt{m+1}&-\Delta\end{pmatrix},\quad\Omega_m=\sqrt{\Delta^2+4\lambda^2(m+1)},\quad e^{-i\hat h_mt/\hbar}=\hat I_2\cos(\Omega_mt/2)-i\frac{2\hat h_m}{\hbar\Omega_m}\sin(\Omega_mt/2). $$

Audit: Correct; dressed energies hbar omega(m+1/2) +- hbar Omega_m/2 (atlas's correction of printed Eq. 4.145 holds). Checked by: exact WL computation (script wl/039.wl): (2h)^2 = Omega_m^2 I, the closed-form propagator equals MatrixExp, and the eigenvalues are +-Omega_m/2.
QF code (exact): none. a 2x2 atomic block with no field operator left; nothing for QF's bosonic algebra to do.

<a id="passive-two-port-transform-must-be-unitary"></a>
### Passive two-port transform must be unitary
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

<a id="beam-splitter-with-i-reflection-exponential-and-adjoint-action"></a>
### Beam splitter with i-reflection: exponential and adjoint action
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

<a id="real-beam-splitter-rotation"></a>
### Real beam-splitter rotation
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

<a id="schwinger-su-2-algebra-casimir-and-j3-phase-action"></a>
### Schwinger SU(2) algebra, Casimir and J3 phase action
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

<a id="beam-splitter-output-of-a-number-state-binomial-state"></a>
### Beam-splitter output of a number state (binomial state)
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

<a id="coherent-products-stay-products-under-passive-mixing"></a>
### Coherent products stay products under passive mixing
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

<a id="two-photon-interference-and-twin-fock-output-at-50-50"></a>
### Two-photon interference and twin-Fock output at 50:50
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

<a id="cross-kerr-gate-conjugated-by-a-50-50-beam-splitter"></a>
### Cross-Kerr gate conjugated by a 50:50 beam splitter
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

<a id="frequency-conversion-as-an-su-2-rotation"></a>
### Frequency conversion as an SU(2) rotation
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

<a id="bosonic-canonical-transformations"></a>
## Bosonic canonical transformations

<a id="multimode-canonical-action-and-its-constraints"></a>
### Multimode canonical action and its constraints
Appears in: X-009 (source records cited: B-049, B-060, B-061, B-062, P-036, P-042, P-043); P-034

A Gaussian unitary maps the annihilators to linear combinations of annihilators and creators, and the two matrix constraints say that the new operators obey the same commutation (or anticommutation) relations.

$$ \hat U\hat{\mathbf d}\hat U^{-1}=u\hat{\mathbf d}+v\hat{\mathbf d}^{\dagger},\qquad uu^\dagger-\eta vv^\dagger=I_N,\quad uv^T=\eta vu^T $$

Conditions: $\eta=+1$ bosons, $-1$ fermions; $\hat{\mathbf d}^\dagger$ a column of creators.
Audit: Correct (the two constraints are the CCR/CAR preserved by the map). Checked by: QF exact (bosons: with symbolic $u,v$, $[\hat b_i,\hat b_j^\dagger]=(uu^\dagger-vv^\dagger)_{ij}$ and $[\hat b_i,\hat b_j]=(uv^T-vu^T)_{ij}$ for $\hat b=u\hat a+v\hat a^\dagger$; and QF's exact Bogoliubov matrices of a network of two single-mode squeezers and a beam splitter satisfy both constraints exactly); exact fermion matrices (2 modes, symbolic $u,v$: $\{\hat b_i,\hat b_j^\dagger\}=(uu^\dagger+vv^\dagger)_{ij}$, $\{\hat b_i,\hat b_j\}=(uv^T+vu^T)_{ij}$; 3 and 4 modes, $u,v$ read off from $e^{-i\hat H}$ of a random complex even quadratic $\hat H$: both constraints hold to machine precision).
QF code (exact):
```wl
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
vars = {a1, ad1, a2, ad2};
b = Table[Sum[u[i, j] {a1, a2}[[j]] + v[i, j] {ad1, ad2}[[j]], {j, 2}], {i, 2}];
bd = Table[Sum[Conjugate[u[i, j]] {ad1, ad2}[[j]] + Conjugate[v[i, j]] {a1, a2}[[j]], {j, 2}], {i, 2}];
uu = Array[u, {2, 2}]; vv = Array[v, {2, 2}];
qi = QuantumInterferometer[{"S"[1/3 Exp[I/2]] -> 1, "S"[1/4 Exp[-I]] -> 2, "BS"[3/5, 9/10] -> {1, 2}}, 2];
{al, be, c} = First[qi]["Bogoliubov"];
{uq, vq} = ConjugateTranspose /@ {al, be};
{Simplify[Table[BosonicNormalOrder[b[[i]] ** bd[[j]] - bd[[j]] ** b[[i]], vars], {i, 2}, {j, 2}] - (uu . ConjugateTranspose[uu] - vv . ConjugateTranspose[vv])],
 Simplify[Table[BosonicNormalOrder[b[[i]] ** b[[j]] - b[[j]] ** b[[i]], vars], {i, 2}, {j, 2}] - (uu . Transpose[vv] - vv . Transpose[uu])],
 FullSimplify[uq . ConjugateTranspose[uq] - vq . ConjugateTranspose[vq] - IdentityMatrix[2]], FullSimplify[uq . Transpose[vq] - vq . Transpose[uq]]}
```
Output: `{{{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}}`. Matches: the first two are the commutators minus the constraint matrices (symbolic), the last two are the constraints on QF's exact $u=\alpha^\dagger$, $v=\beta^\dagger$ (QF stores the creator action $\hat U\hat a_k^\dagger\hat U^\dagger=\sum_j\alpha_{jk}\hat a_j^\dagger+\beta_{jk}\hat a_j$).

**Also stated in P-034 (Finite-mode bosonic Bogoliubov conditions):**

A linear canonical map of $N$ bosonic modes preserves the indefinite form $K$, which fixes the block relations among $U$ and $V$, the inverse, and makes $Z=U^{-1}V$ symmetric and a strict contraction.

$$ \hat{\mathbf b}=U\hat{\mathbf a}+V\hat{\mathbf a}^\dagger:\ UU^\dagger-VV^\dagger=I_N,\ UV^T=VU^T,\ MKM^\dagger=K,\ M^{-1}=\begin{pmatrix}U^\dagger&-V^T\\-V^\dagger&U^T\end{pmatrix},\ U^\dagger U-V^T\bar V=I_N,\ U^\dagger V=V^T\bar U,\ Z=U^{-1}V=Z^T,\ I_N-ZZ^\dagger=(U^\dagger U)^{-1} $$

Audit: Correct. Checked by: exact WL computation, $N=2$, with an exact group element obtained as the Cayley transform $(I+L/2)(I-L/2)^{-1}$ of a rational generator $L$ of the P-035 form: all ten relations hold with zero residual. Forcing line for the last one: $ZZ^\dagger=U^{-1}(UU^\dagger-I)U^{-\dagger}=I-(U^\dagger U)^{-1}$.
QF code (exact): none. QF has no function that tests Bogoliubov coefficient matrices.

<a id="multimode-creator-number-annihilator-factorization"></a>
### Multimode creator-number-annihilator factorization
Appears in: X-009

Any Gaussian unitary with invertible $u$ factors into a pair-creation exponential, a number-conserving exponential, and a pair-annihilation exponential, times its vacuum amplitude.

$$ \hat U=c_0\exp\!\left[-\tfrac12\hat{\mathbf d}^{\dagger T}Z\hat{\mathbf d}^{\dagger}\right]\exp\!\left[\hat{\mathbf d}^{\dagger T}h\hat{\mathbf d}\right]\exp\!\left[\tfrac\eta2\hat{\mathbf d}^{T}W\hat{\mathbf d}\right],\quad Z=u^{-1}v,\ W=\bar v u^{-1},\ e^h=u^{-1},\ c_0=\langle0|\hat U|0\rangle,\ |c_0|=\det(uu^\dagger)^{-\eta/4} $$

Conditions: $u$ invertible (always so for bosons, since $uu^\dagger=I+vv^\dagger$; fermionic $u$ can be singular, e.g. the pair rotation of X-004 at $r=\pi/2$). The phase of $c_0$ is not fixed by $u,v$.
Audit: Correct, including the $\eta/2$ sign of the pair-annihilation factor ($Z$ and $W$ come out symmetric for bosons, antisymmetric for fermions). Checked by: hand derivation for bosons: QF gives $e^{\hat X}\hat a e^{-\hat X}=\hat a+Z\hat a^\dagger$ and $e^{\hat Y}\hat a^\dagger e^{-\hat Y}=\hat a^\dagger+W\hat a$ for the outer exponents $\hat X,\hat Y$, and the middle factor sends $\hat a\to u\hat a$, $\hat a^\dagger\to(u^{-1})^T\hat a^\dagger$; so the product sends $\hat a\to u\hat a+v\hat a^\dagger$ and $\hat a^\dagger\to W u\,\hat a+((u^{-1})^T+WuZ)\hat a^\dagger$, which is $\bar v\hat a+\bar u\hat a^\dagger$ because $u^T\bar u-v^\dagger v=I$ and $u^T\bar v=v^\dagger u$; canonical action plus the vacuum amplitude fix $\hat U$, and $|c_0|$ follows from the X-006 normalization since $\det(I-ZZ^\dagger)=\det(uu^\dagger)^{-1}$. The needed matrix identity, symmetry of $Z,W$, and $|c_0|=\det(uu^\dagger)^{-1/4}$ with QF's $c$ were confirmed on QF's exact Bogoliubov matrices of a four-element two-mode network, evaluated at 40 digits. Exact fermion matrices (3 and 4 modes): the factorization equals $e^{-i\hat H}$ to machine precision and $|c_0|=\det(uu^\dagger)^{1/4}$; with $+\tfrac12$ in place of $\eta/2$ the mismatch is O(1), so the sign is load-bearing.
QF code (exact):
```wl
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
g = r Exp[-I th] a1 ** a2 - r Exp[I th] ad1 ** ad2;
act = Simplify[{BosonicExpSimplify[Exp[g] ** a1 ** Exp[-g]], BosonicExpSimplify[Exp[g] ** a2 ** Exp[-g]]}, r > 0];
u = Table[Coefficient[act[[i]], x], {i, 2}, {x, {a1, a2}}];
v = Table[Coefficient[act[[i]], x], {i, 2}, {x, {ad1, ad2}}];
z = Simplify[Inverse[u] . v]; w = Simplify[Conjugate[v] . Inverse[u], r > 0 && Element[th, Reals]];
{act, z, w, Inverse[u], Simplify[Det[u . ComplexExpand[ConjugateTranspose[u]]]^(-1/4), r > 0], FullSimplify[BosonicExpOrder[Exp[g]], r > 0 && Element[th, Reals]]}
```
Output: `{{\[FormalA]1*Cosh[r] + E^(I*th)*Sinh[r]*SuperDagger[\[FormalA]2], \[FormalA]2*Cosh[r] + E^(I*th)*Sinh[r]*SuperDagger[\[FormalA]1]}, {{0, E^(I*th)*Tanh[r]}, {E^(I*th)*Tanh[r], 0}}, {{0, Tanh[r]/E^(I*th)}, {Tanh[r]/E^(I*th), 0}}, {{Sech[r], 0}, {0, Sech[r]}}, Sech[r], E^(-(E^(I*th)*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]2]*Tanh[r]))**Sech[r]^(SuperDagger[\[FormalA]1]**\[FormalA]1 + SuperDagger[\[FormalA]2]**\[FormalA]2)**E^((\[FormalA]1**\[FormalA]2*Tanh[r])/E^(I*th))*Sech[r]}`. Matches for the two-mode squeezer $\hat U=\hat S_2(re^{i\theta})$: QF's canonical action gives $Z$, $W$, $u^{-1}=\operatorname{sech}r\,I$ and $\det(uu^\dagger)^{-1/4}=\operatorname{sech}r$, and the atlas product built from them, $c_0e^{-e^{i\theta}\tanh r\,\hat a_1^\dagger\hat a_2^\dagger}\operatorname{sech}r^{\,\hat N_1+\hat N_2}e^{e^{-i\theta}\tanh r\,\hat a_1\hat a_2}$, is QF's own ordering of $\hat U$.

<a id="bosonic-coherent-state-kernel-of-the-factorized-implementer"></a>
### Bosonic coherent-state kernel of the factorized implementer
Appears in: X-009

Between unnormalized coherent states the factorized implementer has a Gaussian kernel whose quadratic form is read off from $Z$, $u^{-1}$ and $W$.

$$ K_{\hat U}(x,y)=c_0\exp\!\left[-\tfrac12x^TZx+x^Tu^{-1}y+\tfrac12y^TWy\right] $$

Conditions: kernel taken between unnormalized coherent states, $x=\bar z$, $y=w$ in $\langle z|\hat U|w\rangle$ with $|w\rangle=e^{w\cdot\hat a^\dagger}|0\rangle$. The atlas does not state this convention; it is the one that makes the formula hold.
Audit: Correct. Checked by: hand derivation: $\langle z|e^{-\frac12\hat a^{\dagger T}Z\hat a^\dagger}=e^{-\frac12\bar z^TZ\bar z}\langle z|$, $e^{\hat a^{\dagger T}h\hat a}|w\rangle=|e^hw\rangle$, $e^{\frac12\hat a^TW\hat a}|w\rangle=e^{\frac12w^TWw}|w\rangle$, and $\langle z|w'\rangle=e^{\bar z^Tw'}$; and QF exact: QF's package-private `gaussianKernel` on exact Bogoliubov matrices returns exactly the atlas blocks, and the single-mode middle line is QF's coherent matrix element of $e^{h\hat a^\dagger\hat a}$.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
qi = QuantumInterferometer[{"S"[1/3 Exp[I/2]] -> 1, "S"[1/4 Exp[-I]] -> 2, "BS"[3/5, 9/10] -> {1, 2}}, 2];
{al, be, c} = First[qi]["Bogoliubov"];
{u, v} = ConjugateTranspose /@ {al, be};
gk = Symbol[First[Names["Wolfram`QuantumFramework`*`gaussianKernel"]]][qi];
{FullSimplify[First[gk] - ArrayFlatten[{{-Inverse[u] . v, Inverse[u]}, {Transpose[Inverse[u]], Conjugate[v] . Inverse[u]}}]], Last[gk] === c, BosonicMatrixElement[{z, w}, Exp[h ad ** a], "Basis" -> "Coherent"]}
```
Output: `{{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}, True, E^((-Abs[w]^2 - Abs[z]^2)/2 + E^h*w*Conjugate[z])}`. Matches: QF's kernel matrix minus $\begin{pmatrix}-Z&u^{-1}\\(u^{-1})^T&W\end{pmatrix}$ is exactly zero, its scalar is QF's $c$, and the normalized coherent element of $e^{h\hat a^\dagger\hat a}$ is $e^{e^h\bar zw}$ times the normalization factors.

---

<a id="bosonic-affine-canonical-transformation"></a>
### Bosonic affine canonical transformation
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

<a id="bosonic-block-relations-and-inverse"></a>
### Bosonic block relations and inverse
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

<a id="symplectic-form-and-inverse-in-quadratures"></a>
### Symplectic form and inverse in quadratures
Appears in: B-044

In real quadrature coordinates the transformation matrix is symplectic, and its inverse is obtained by transposing and rearranging the blocks.

$$ M^TJM=J,\qquad M^{-1}=\begin{pmatrix}D^T&-B^T\\-C^T&A^T\end{pmatrix} $$

Audit: Correct. Checked by: hand derivation: M^{-1} = -J M^T J, which multiplies out to the displayed blocks for J = ((0, I), (-I, 0)).
QF code (exact): none. A statement about real 2N x 2N matrices; no second-quantization call is involved.

<a id="squeezing-contraction-both-charts"></a>
### Squeezing contraction, both charts
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

<a id="full-bosonic-gaussian-kernel"></a>
### Full bosonic Gaussian kernel
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

<a id="general-bosonic-transformation-kernel-for-normal-symbols"></a>
### General bosonic transformation kernel for normal symbols
Appears in: B-054

For a mixing transformation, the conjugated normal symbol is a limit of Gaussian integrals of the original symbol against a kernel built from the displayed matrix L.

$$ N_{\hat U\hat A\hat U^{-1}}=\lim_n\int\mathcal K_n(a^*,a\mid\alpha^*,\alpha)N_{\hat A}(\alpha^*,\alpha)\,d^nt\,d^n\tau,\qquad \widehat{\mathcal L}=\begin{pmatrix}\bar{\hat\Phi}\bar{\hat\Psi}^{-1}&-\hat I\\-\hat I&\hat\Phi\hat\Psi^{-1}\end{pmatrix} $$

Audit: Unchecked: the Gaussian kernel (4.39) and its contour are not transcribed in the atlas, so there is nothing complete to test. Checked by: prior audit only.
QF code (exact): none. The kernel (4.39) is not transcribed, so there is no formula to compare a QF result against.

<a id="gaussian-normal-symbols-stay-gaussian-under-conjugation"></a>
### Gaussian normal symbols stay Gaussian under conjugation
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

<a id="bosonic-canonical-diagonalization"></a>
### Bosonic canonical diagonalization
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

<a id="bogoliubov-combination-keeps-the-canonical-commutator"></a>
### Bogoliubov combination keeps the canonical commutator
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

<a id="quadratic-bosonic-generators-as-sp-2n-r-matrices"></a>
### Quadratic bosonic generators as sp(2N,R) matrices
Appears in: P-035

A quadratic bosonic generator corresponds to a block matrix with anti-Hermitian number block and symmetric pairing block.

$$ L=\begin{pmatrix}C&A\\\bar A&\bar C\end{pmatrix},\qquad C^\dagger=-C,\quad A^T=A $$

Audit: Correct. Checked by: hand derivation, $LK+KL^\dagger$ has blocks $C+C^\dagger$ and $A^T-A$, so it vanishes exactly when $C^\dagger=-C$ and $A^T=A$; exact WL computation confirms that block form for a general symbolic $L$.
QF code (exact): none. QF has no sp(2N) coefficient-matrix layer.

<a id="bosonic-quadratic-commutator-table-e-f-g"></a>
### Bosonic quadratic commutator table (E, F, G)
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

<a id="fractional-linear-action-on-symmetric-gaussian-parameters"></a>
### Fractional-linear action on symmetric Gaussian parameters
Appears in: P-038

Canonical maps act on the symmetric Gaussian parameter by a matrix fractional-linear map, and the analytic representation carries a determinant multiplier.

$$ Z\mapsto (U+Z\bar V)^{-1}(V+Z\bar U),\qquad [\hat T^{(+)}(g)f](Z)=\det(U+Z\bar V)^{-1/2}f\left((U+Z\bar V)^{-1}(V+Z\bar U)\right) $$

Audit: Correct for the map; the determinant multiplier is Unchecked (its sign and branch depend on the source's function-space convention). Checked by: exact WL computation ($N=2$, two exact Cayley group elements, rational $Z$): the image is symmetric, $I-Z'Z'^\dagger$ has positive eigenvalues, and $g_1\cdot(g_2\cdot Z)=(g_2g_1)\cdot Z$ with zero residual, which is the composition order that makes $f\mapsto f(g\cdot Z)$ a representation.
QF code (exact): none. QF has no Gaussian-parameter action.

<a id="siegel-upper-half-space-and-symplectic-fractional-linear-action"></a>
### Siegel upper half-space and symplectic fractional-linear action
Appears in: P2-002

Real symplectic matrices act on the Siegel upper half-space by a fractional-linear map, and their blocks obey three quadratic relations.

$$ \mathfrak H_p=\{Z=X+iY:Z^T=Z,\ Y>0\},\quad Z\mapsto(AZ+B)(CZ+D)^{-1},\quad A^TC=C^TA,\ B^TD=D^TB,\ A^TD-C^TB=I_p $$

Audit: Correct. Checked by: exact WL computation: an exact rational real symplectic $4\times4$ matrix (Cayley transform of $JS$, $S$ symmetric) satisfies the block relations with zero residual, and the image of a point of the half-space is symmetric with positive-definite imaginary part.
QF code (exact): none. Not an operator identity.

<a id="multimode-particle-antiparticle-bosonic-blocks"></a>
### Multimode particle-antiparticle bosonic blocks
Appears in: P2-020

A canonical map mixing particles and antiparticles preserves the indefinite form $E$, which fixes the block relations and the inverse.

$$ MEM^\dagger=E,\quad AA^\dagger-BB^\dagger=I_N,\quad DD^\dagger-CC^\dagger=I_N,\quad AC^\dagger=BD^\dagger,\quad M^{-1}=EM^\dagger E $$

Audit: Correct. Checked by: hand derivation, blockwise expansion of $MEM^\dagger=E$; same structure as the finite-mode bosonic Bogoliubov entry with $(\hat a,\hat b^\dagger)$ in place of $(\hat a,\hat a^\dagger)$.
QF code (exact): none. QF has no coefficient-matrix layer.

<a id="symplectic-covariance-of-weyl-symbols"></a>
### Symplectic covariance of Weyl symbols
Appears in: C-015

A metaplectic operator transforms Weyl symbols and Wigner functions by composing them with the classical linear canonical map.

$$ \hat R(F)^*\hat A\hat R(F)=\widehat{A\circ F},\qquad W_{\hat R(F)\psi,\hat R(F)\eta}(z)=W_{\psi,\eta}(F^{-1}z) $$

Audit: Correct (the second follows from the first through C-014). Checked by: hand derivation.
QF code (exact): none. QF has no Weyl-symbol calculus.

<a id="metaplectic-generator-from-a-symplectic-path"></a>
### Metaplectic generator from a symplectic path
Appears in: C-023

The quadratic generator of a symplectic path is read off as -J F' F^{-1}, and implementers compose only up to a sign (the metaplectic double cover).

$$ S_t=-J\dot F_tF_t^{-1},\qquad \hat R(F_1)\hat R(F_2)=\pm\hat R(F_1F_2) $$

Audit: Correct (from F' = J S F and J^{-1} = -J; the sign is the metaplectic double cover). Checked by: hand derivation.
QF code (exact): none. Classical matrix statement; QF has no metaplectic-lift routine.

<a id="bogoliubov-blocks-of-a-symplectic-matrix"></a>
### Bogoliubov blocks of a symplectic matrix
Appears in: C-024

The Heisenberg action of a quadratic propagator on the annihilator is a Bogoliubov mixing u a + v a^dag whose blocks are fixed combinations of the classical flow blocks.

$$ \hat U^*\hat a\hat U=u\hat a+v\hat a^\dagger,\quad u=\tfrac12[A+D+i(C-B)],\quad v=\tfrac12[A-D+i(C+B)],\quad uu^\dagger-vv^\dagger=I_n,\ uv^T=vu^T $$

Audit: Correct. Checked by: QF exact (one mode, symbolic S, hbar = 1): u and v read off U^* a U equal (A+D+i(C-B))/2 and (A-D+i(C+B))/2 for the blocks of e^{tJS}, and |u|^2 - |v|^2 = 1; for several modes u v^T = v u^T is the symplectic condition written in these blocks (hand).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2]; p = -I (a - ad)/Sqrt[2];
h = (s11 q ** q + s12 (q ** p + p ** q) + s22 p ** p)/2;
heis = BosonicExpSimplify[Exp[I t h] ** a ** Exp[-I t h]];
{{fa, fb}, {fc, fd}} = MatrixExp[t {{0, 1}, {-1, 0}} . {{s11, s12}, {s12, s22}}];
{u, v} = {Coefficient[heis, a], Coefficient[heis, ad]};
Simplify[{u - (fa + fd + I (fc - fb))/2, v - (fa - fd + I (fc + fb))/2,
  ComplexExpand[u Conjugate[u] - v Conjugate[v]]}, t > 0 && s12^2 - s11 s22 > 0]
```
Output: `{0, 0, 1}`. Matches: both differences vanish and |u|^2 - |v|^2 = 1.

<a id="fermionic-canonical-transformations"></a>
## Fermionic canonical transformations

<a id="fermion-pair-realization-of-su-2"></a>
### Fermion-pair realization of su(2)
Appears in: X-004 (marked by the atlas as derived); P2-025; C-087

Creating or destroying a fermion pair in two modes, together with half the shifted number, closes on the compact su(2) algebra.

$$ \hat J_+=\hat c_1^\dagger \hat c_2^\dagger,\quad \hat J_-=\hat c_2\hat c_1=\hat J_+^\dagger,\quad \hat J_0=(\hat N_1+\hat N_2-\hat I)/2 $$

Audit: Correct; these satisfy the compact algebra. The order $\hat c_2\hat c_1$ is essential, as the atlas says: with $\hat c_1\hat c_2$ (which equals $-\hat J_+^\dagger$) the relation $[\hat J_+,\hat J_-]=2\hat J_0$ fails. Checked by: exact fermion matrices (2 modes).
QF code (exact): none. no fermionic layer.

**Also stated in P2-025 (Single fermionic pair as spin 1/2):**

One fermionic particle-antiparticle pair is a spin-1/2 problem: the pair operators act as $\hat J_\pm,\hat J_0$ on the vacuum and the occupied pair.

$$ \hat J_+=\hat a^\dagger\hat b^\dagger,\ \hat J_-=\hat b\hat a,\ \hat J_0=\tfrac12(\hat a^\dagger\hat a+\hat b^\dagger\hat b-\hat I),\qquad |\zeta\rangle=(1+|\zeta|^2)^{-1/2}(\hat I+\zeta\hat a^\dagger\hat b^\dagger)|0\rangle $$

Audit: Correct. Checked by: hand derivation; on $\{|0\rangle,\hat a^\dagger\hat b^\dagger|0\rangle\}$, $[\hat J_+,\hat J_-]=\mathrm{diag}(-1,1)=2\hat J_0$; exact fermion matrices (4 modes) confirm $[\hat J_+,\hat J_-]-2\hat J_0$ annihilates the vacuum. The Dirac-field reductions are Unchecked.
QF code (exact): none. no fermionic layer.

**Also stated in C-087 (Fermion pair rotation as SU(2)):**

Fermion pair creation and annihilation close an su(2), so a pair rotation sends the vacuum to cos r times the vacuum plus sin r times a filled pair.

$$ \hat J_+=\hat c_1^\dagger\hat c_2^\dagger,\ \hat J_-=\hat c_2\hat c_1,\ \hat J_3=\tfrac12(\hat N_1+\hat N_2-1),\qquad \hat U\Omega=\cos r\,\Omega+\frac\xi r\sin r\,\hat c_1^\dagger\hat c_2^\dagger\Omega $$

Audit: Correct. Checked by: exact fermion matrices (2 modes: su(2) relations and the vacuum image).
QF code (exact): none. no fermionic layer.

<a id="fermion-pair-rotation-of-the-vacuum"></a>
### Fermion-pair rotation of the vacuum
Appears in: X-004

The pair rotation turns the empty state into a superposition of empty and doubly occupied, with $\cos r$ and $\sin r$ amplitudes.

$$ e^{\xi \hat c_1^\dagger \hat c_2^\dagger-\bar\xi \hat c_2\hat c_1}|00\rangle=\cos r\,|00\rangle+e^{i\phi}\sin r\,|11\rangle $$

Conditions: $|11\rangle=\hat c_1^\dagger\hat c_2^\dagger|00\rangle$ (this ordering fixes the sign of the second term). It holds for all r; no chart condition is needed for the state.
Audit: Correct. Checked by: hand derivation: on span{|00>, |11>} the generator is the matrix ((0, -ξ̄), (ξ, 0)), whose exponential sends |00> to cos r |00> + (ξ/r) sin r |11>; and exact fermion matrices (2 modes, symbolic $r,\phi$, difference exactly zero).
QF code (exact): none. no fermionic layer.

<a id="norm-bound-for-fermionic-fields"></a>
### Norm bound for fermionic fields
Appears in: B-005

A smeared fermionic field is a bounded operator whose norm is at most the norm of its test function; in fact the bound is attained.

$$ \|\hat a(f)\|,\ \|\hat a^*(f)\|\le\|f\| $$

Audit: Correct, with equality: a(f)^2 = 0 and {a(f), a(f)*} = ||f||^2 make a(f)* a(f) / ||f||^2 a projector. Checked by: hand derivation: (a(f)* a(f))^2 = a(f)* (||f||^2 - a(f)* a(f)) a(f) = ||f||^2 a(f)* a(f), and exact fermion matrices (2 modes, operator norm equal to ||f||).
QF code (exact): none. No fermionic layer.

<a id="fermionic-gaussian-implementer-kernel"></a>
### Fermionic Gaussian implementer kernel
Appears in: B-060

The fermionic implementer's kernel is a Grassmann Gaussian with the same block structure as the bosonic one, but with prefactor det(Phi Phi*)^{+1/4}.

$$ K_{\hat U}=\theta\det(\hat\Phi\hat\Phi^*)^{1/4} \exp\!\left[-\tfrac12\left(a\bar{\hat\Psi}\hat\Phi^{-1}a-2a^*\hat\Phi^{-1}a+a^*\hat\Phi^{-1}\hat\Psi a^*\right)\right], \qquad N_{\hat U}=K_{\hat U}e^{-a^*a} $$

Conditions: two-sided Phi^{-1}; Grassmann order as displayed.
Audit: Correct. Checked by: exact fermion matrices (2 modes): the Grassmann generating function Sum_{m,n} <m|U|n> (a*)^m ((a*)^n)* of the Jordan-Wigner U equals the displayed exponential times one scalar (all 16 coefficients), and that scalar has modulus det(Phi Phi*)^{1/4}.
QF code (exact): none. No fermionic layer.

<a id="fermionic-transformation-of-normal-symbols-with-mixing"></a>
### Fermionic transformation of normal symbols with mixing
Appears in: B-063

For a mixing fermionic transformation the conjugated normal symbol is a Grassmann Gaussian integral of the original one, with the displayed matrix and determinant factor.

$$ \widehat{\mathcal L}=\begin{pmatrix}\bar{\hat\Phi}\bar{\hat\Psi}^{-1}&-\hat I\\\hat I&-\hat\Phi\hat\Psi^{-1}\end{pmatrix},\qquad \det(\hat P_n\hat\Psi\hat\Psi^*\hat P_n)^{1/2} $$

Audit: Unchecked: the integral kernel (5.19)-(5.20) is described but not transcribed; the passive case (substitution) is correct by the same argument as B-053. Checked by: prior audit only for the mixing part.
QF code (exact): none. No fermionic layer.

<a id="fermionic-quadratic-propagator"></a>
### Fermionic quadratic propagator
Appears in: B-074

For fermions the propagator kernel is det(Phi e^{itC})^{+1/2} times a Grassmann Gaussian built from the same classical generator.

$$ K_{e^{it\hat H}}=\det(\hat\Phi e^{it\hat C})^{1/2} \exp\!\left[-\tfrac12(a\ a^*) \begin{pmatrix}\bar{\hat\Psi}\hat\Phi^{-1}&\hat\Phi'^{-1}\\-\hat\Phi^{-1}&\hat\Phi^{-1}\hat\Psi\end{pmatrix}\binom a{a^*}\right] $$

Audit: Correct. The exponent equals that of B-060 because a Phi'^{-1} a* = -a* Phi^{-1} a for Grassmann variables; the same generator L holds for fermions. Checked by: exact fermion matrices (2 modes: Phi, Psi from U c U^dag equal the blocks of exp(itL); <00|U|00> = det(Phi e^{itC})^{1/2} at t = 0.7).
QF code (exact): none. No fermionic layer.

<a id="fermionic-diagonalization-conditions"></a>
### Fermionic diagonalization conditions
Appears in: B-086

A fermionic quadratic Hamiltonian is diagonalized by a proper transformation when S and S^2 - [S, T] - T^2 are bounded below by a positive constant.

$$ \hat S>\mu \hat I>0,\qquad \hat S^2-[\hat S,\hat T]-\hat T^2>\mu \hat I>0 $$

Audit: Consistent in finite dimension: with T real skew, S^2 - [S,T] - T^2 = (S+T)(S+T)^T, so the second condition is a gap condition (no zero quasiparticle energy). The infinite-mode theorem is Unchecked. Checked by: hand derivation: (S + T)(S - T) = S^2 - ST + TS - T^2 and (S + T)^T = S - T.
QF code (exact): none. No fermionic layer.

<a id="fermionic-bogoliubov-constraints"></a>
### Fermionic Bogoliubov constraints
Appears in: P-041; B-056 (also B-057)

A fermionic canonical map keeps the anticommutators exactly when the doubled coefficient matrix is unitary, and $Z=U^{-1}V$ is then skew-symmetric.

$$ \{\hat c_i,\hat c_j^\dagger\}=\delta_{ij},\ \hat{\mathbf b}=U\hat{\mathbf c}+V\hat{\mathbf c}^\dagger:\ UU^\dagger+VV^\dagger=I_N,\ UV^T+VU^T=0,\ M^\dagger M=I_{2N},\ Z=U^{-1}V=-Z^T $$

Audit: Correct. Checked by: exact fermion matrices (3 modes), with an exact unitary group element from the Cayley transform of a rational anti-Hermitian generator of block form: the new operators satisfy the CAR exactly and all block relations hold with zero residual.
QF code (exact): none. no fermionic layer.

**Also stated in B-056 (Fermionic canonical relations):**

A fermionic canonical transformation has blocks with Phi Phi* + Psi Psi* = I and Phi Psi' + Psi Phi' = 0, the plus signs replacing the bosonic minus signs.

$$ \hat b=\hat\Phi\hat a+\hat\Psi\hat a^*,\quad \hat\Phi\hat\Phi^*+\hat\Psi\hat\Psi^*=\hat I,\quad\hat\Phi\hat\Psi'+\hat\Psi\hat\Phi'=0,\qquad \hat S^{-1}=\begin{pmatrix}\hat\Phi^*&\hat\Psi'\\\hat\Psi^*&\hat\Phi'\end{pmatrix} $$

Audit: Correct. Checked by: exact fermion matrices (2 modes; U = exp(itH) for complex pairing and complex Hermitian hopping; Phi, Psi read off U c U^dag by trace projection; all relations hold).
QF code (exact): none. No fermionic layer.

<a id="fermionic-quadratic-commutator-table"></a>
### Fermionic quadratic commutator table
Appears in: P-042

Fermionic bilinears close under commutators with structure constants that differ from the bosonic table in the $[\hat F,\hat G]$ signs.

$$ [\hat F_{ij},\hat G_{kl}]=-\delta_{jk}\hat E_{li}+\delta_{ik}\hat E_{lj}+\delta_{jl}\hat E_{ki}-\delta_{il}\hat E_{kj},\quad [\hat F,\hat F]=[\hat G,\hat G]=0,\ \text{first three bosonic relations unchanged} $$

Conditions: $\hat E_{ij}=\hat c_i^\dagger\hat c_j-\delta_{ij}/2$, $\hat F_{ij}=\hat c_i\hat c_j$, $\hat G_{ij}=\hat c_i^\dagger\hat c_j^\dagger$.
Audit: Correct. Checked by: exact fermion matrices (2 modes), all 16 index choices, six relations: every residual exactly 0.
QF code (exact): none. no fermionic layer.

<a id="fermionic-analytic-group-action"></a>
### Fermionic analytic group action
Appears in: P-044

Fermionic canonical maps act on skew $Z$ by the same fractional-linear map as bosons, with determinant exponent $+1/2$.

$$ [\hat T^{(+)}(g)f](Z)=\det(U+Z\bar V)^{1/2}f\left((U+Z\bar V)^{-1}(V+Z\bar U)\right) $$

Audit: Unchecked: the determinant exponent depends on the source's function-space realization, which is not reconstructed here. The skew-symmetry of the image follows by the same algebra as the bosonic map with signs from $UV^T=-VU^T$.
QF code (exact): none. no fermionic layer.

<a id="multimode-fermionic-blocks"></a>
### Multimode fermionic blocks
Appears in: P2-026

The doubled fermionic coefficient matrix is unitary, which fixes the block relations.

$$ MM^\dagger=I_{2N},\quad AA^\dagger+BB^\dagger=I_N,\quad DD^\dagger+CC^\dagger=I_N,\quad AC^\dagger=-BD^\dagger $$

Audit: Correct. Checked by: hand derivation, blockwise expansion of $MM^\dagger=I$ (same algebra as the exact fermionic Bogoliubov check).
QF code (exact): none. no fermionic layer.

<a id="car-and-finite-fermionic-fock-space"></a>
### CAR and finite fermionic Fock space
Appears in: C-064

Fermionic modes satisfy canonical anticommutation relations, so n modes span a 2^n-dimensional Fock space and a smeared annihilator has operator norm ||f||.

$$ \{\hat c_j,\hat c_k^\dagger\}=\delta_{jk}\hat I,\ \{\hat c_j,\hat c_k\}=0,\ \|\hat c(f)\|=\|f\| $$

Audit: Correct. Checked by: exact fermion matrices (3 modes) for the CAR; hand derivation for the norm (c(f)^dag c(f) has eigenvalues 0 and ||f||^2 by the CAR).
QF code (exact): none. no fermionic layer.

<a id="linear-fermionic-exponential-with-complex-coefficients"></a>
### Linear fermionic exponential with complex coefficients
Appears in: C-065

A linear combination of fermionic creators and annihilators with complex coefficients squares to the scalar u.v, so its exponential is cosh plus a sinh-over-root multiple of the generator.

$$ e^{\hat L}=\cosh\sqrt{u\cdot v}\,\hat I+\frac{\sinh\sqrt{u\cdot v}}{\sqrt{u\cdot v}}(u\cdot\hat c^\dagger+v\cdot\hat c),\qquad \hat L^2=(u\cdot v)\hat I $$

Audit: Correct. Checked by: exact fermion matrices (2 modes, complex u, v).
QF code (exact): none. no fermionic layer.

<a id="bdg-hamiltonian-and-its-canonical-matrix"></a>
### BdG Hamiltonian and its canonical matrix
Appears in: C-080

A general quadratic fermion Hamiltonian is encoded by a 2n-by-2n BdG matrix Sigma, whose exponential gives the Heisenberg action on (c, c^dag).

$$ \hat H=\tfrac12(\hat c^\dagger M\hat c+\hat c^\dagger L\hat c^\dagger-\hat c\bar L\hat c)-\tfrac14(\operatorname{Tr}M)\hat I,\quad \Sigma=\begin{pmatrix}M/2&L\\-\bar L&-\bar M/2\end{pmatrix},\quad e^{-it\hat H}\binom{\hat c}{\hat c^\dagger}e^{it\hat H}=e^{it\Sigma}\binom{\hat c}{\hat c^\dagger} $$

Audit: Correct. Checked by: exact fermion matrices (2 modes, random Hermitian M, antisymmetric complex L; H Hermitian; Heisenberg action matches).
QF code (exact): none. no fermionic layer.

<a id="fermionic-pair-vacuum"></a>
### Fermionic pair vacuum
Appears in: C-082

A quadratic fermion propagator maps the vacuum to a pair Gaussian with antisymmetric Gamma_t and amplitude det(A_t)^{1/2}, and conjugates translations by the inverse canonical matrix.

$$ \hat U_t\Omega=\det(A_t)^{1/2}\exp\!\left[\tfrac12\hat c^\dagger\Gamma_t\hat c^\dagger\right]\Omega,\quad \Gamma_t=B_t^T(A_t^T)^{-1}=-\Gamma_t^T;\qquad \hat U_t\hat T(\alpha)\hat U_t^{-1}=\hat T(F_t^{-1}\alpha) $$

Audit: Correct (the atlas's inverse label). Checked by: exact fermion matrices for the pair vacuum (2 modes, t = 0.6); hand derivation for the label in one number mode (U c U^{-1} = e^{itm/2} c, so T(alpha) maps to T(e^{-itm/2} alpha) = T(F^{-1} alpha)).
QF code (exact): none. no fermionic layer.

<a id="fermionic-riccati-and-vacuum-amplitude"></a>
### Fermionic Riccati and vacuum amplitude
Appears in: C-083

The fermionic pairing matrix obeys a Riccati equation and the vacuum amplitude a linear equation with a trace, solved by det(A_t)^{1/2}.

$$ i\dot\Gamma=\tfrac12(M\Gamma+\Gamma\bar M)+L+\Gamma\bar L\Gamma,\quad i\dot s=-\tfrac12\operatorname{Tr}(\bar L\Gamma+M/2)s,\quad s=\det(A_t)^{1/2}=\langle\Omega|\hat U_t|\Omega\rangle $$

Audit: Correct. Checked by: exact fermion matrices (2 modes, the BdG data of C-080; derivatives by central differences).
QF code (exact): none. no fermionic layer.

<a id="fermionic-mehler-kernel"></a>
### Fermionic Mehler kernel
Appears in: C-084

The holomorphic kernel of a fermionic Gaussian unitary is the vacuum amplitude times a Grassmann Gaussian built from A^{-1} and B.

$$ K_{\hat U}(a^*,a)=a_0\exp\!\left[-\tfrac12\left(a^T\bar BA^{-1}a-2a^{*T}A^{-1}a+a^{*T}A^{-1}Ba^*\right)\right],\quad a_0=\langle0|\hat U|0\rangle $$

Audit: Unchecked here: needs the full Grassmann kernel convention; prior audit only (finite CAR checks).
QF code (exact): none. no fermionic layer.

<a id="gaussian-vacua-and-overlaps"></a>
## Gaussian vacua and overlaps

<a id="bosonic-multimode-gaussian-vacuum"></a>
### Bosonic multimode Gaussian vacuum
Appears in: X-006 (source records cited: P-036 to P-045, B-048, B-061); P-036; B-048

The multimode squeezed vacuum is a normalized pair-creation exponential on the vacuum, annihilated by the combinations $\hat a+Z\hat a^\dagger$.

$$ |Z\rangle_B=\det(I_N-ZZ^\dagger)^{1/4}e^{-\frac12\hat a^{\dagger T}Z\hat a^\dagger}|0\rangle,\qquad (\hat a+Z\hat a^\dagger)|Z\rangle=0 $$

Conditions: $Z^T=Z$, $I_N-ZZ^\dagger>0$.
Audit: Correct (normalization and annihilation condition). Consistent with X-009 ($Z=u^{-1}v$) and X-002. Checked by: QF exact for the annihilation condition: with $\hat X=-\frac12\hat a^{\dagger T}Z\hat a^\dagger$, QF gives $e^{\hat X}\hat a e^{-\hat X}=\hat a+Z\hat a^\dagger$, so $(\hat a+Z\hat a^\dagger)e^{\hat X}|0\rangle=e^{\hat X}\hat a|0\rangle=0$. Exact WL computation for the norm: in the Bargmann realization the squared norm is the real Gaussian integral of $\exp(-|z|^2-\operatorname{Re}z^TZz)$, whose real $2N\times2N$ form $A=\begin{pmatrix}I+X&-Y\\-Y&I-X\end{pmatrix}$ ($Z=X+iY$) has $\det A=\det(I_N-Z\bar Z)$ identically for $N=1,2,3$, so the squared norm before the prefactor is $\det(I_N-ZZ^\dagger)^{-1/2}$.
QF code (exact):
```wl
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
x = -(z11 ad1 ** ad1 + 2 z12 ad1 ** ad2 + z22 ad2 ** ad2)/2;
{BosonicExpSimplify[Exp[x] ** a1 ** Exp[-x]], BosonicExpSimplify[Exp[x] ** a2 ** Exp[-x]]}
```
Output: `{\[FormalA]1 + z11*SuperDagger[\[FormalA]1] + z12*SuperDagger[\[FormalA]2], \[FormalA]2 + z12*SuperDagger[\[FormalA]1] + z22*SuperDagger[\[FormalA]2]}`. Matches: these are the components of $\hat a+Z\hat a^\dagger$.

**Also stated in P-036 (Bosonic squeezed vacuum chart and its annihilator):**

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

**Also stated in B-048 (Transformed vacuum (homogeneous part)):**

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

<a id="fermionic-multimode-gaussian-vacuum"></a>
### Fermionic multimode Gaussian vacuum
Appears in: X-006; B-061

The fermionic analogue uses an antisymmetric $Z$ and a determinant with a plus sign, raised to the power $-1/4$.

$$ |Z\rangle_F=\det(I_N+ZZ^\dagger)^{-1/4}e^{-\frac12\hat c^{\dagger T}Z\hat c^\dagger}|0\rangle,\qquad (\hat c+Z\hat c^\dagger)|Z\rangle=0 $$

Conditions: $Z^T=-Z$. Exists on the chart where the canonical block $u$ is invertible (the atlas's "$U$" here is undefined in X-006; see Errors found).
Audit: Correct. Checked by: exact fermion matrices (4 modes, exact rational complex antisymmetric $Z$): the annihilation condition is exactly zero and the squared norm of the unnormalized state, squared, equals $\det(I_N+ZZ^\dagger)$ exactly.
QF code (exact): none. no fermionic layer.

**Also stated in B-061 (Fermionic paired vacuum):**

The fermionic transformed vacuum is a pair exponential with skew pairing matrix Phi^{-1} Psi, normalized by det(Phi Phi*)^{1/4}.

$$ F(a^*)=\theta\det(\hat\Phi\hat\Phi^*)^{1/4} \exp[-\tfrac12a^*\hat\Phi^{-1}\hat\Psi a^*] $$

Audit: Correct. Checked by: exact fermion matrices (2 modes: <11|U|00>/<00|U|00> = -(Phi^{-1} Psi)_12).
QF code (exact): none. No fermionic layer.

<a id="transformed-vacuum-affine-normalization-4-25"></a>
### Transformed vacuum: affine normalization (4.25)
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

<a id="bosonic-gaussian-normalization-and-overlap"></a>
### Bosonic Gaussian normalization and overlap
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

<a id="bosonic-gaussian-bilinear-symbols-and-purity-constraints"></a>
### Bosonic Gaussian bilinear symbols and purity constraints
Appears in: P-039

On a pure bosonic Gaussian the pairing matrix $A$ and the number matrix $B$ are rational functions of $Z$ and obey quadratic purity constraints.

$$ A=-\hbar Z(I_N-Z^\dagger Z)^{-1},\quad B=\tfrac\hbar2 I_N+\hbar Z(I_N-Z^\dagger Z)^{-1}Z^\dagger,\quad B^2-A\bar A=\tfrac{\hbar^2}{4}I_N,\quad BA=A\bar B,\quad B^\dagger=B,\quad A^T=A $$

Audit: Convention: on the two-mode state $|Z\rangle$ ($\hbar=1$), $A_{ij}=\langle\hat a_i\hat a_j\rangle$ and $B_{ij}=\langle\hat a_j^\dagger\hat a_i\rangle+\delta_{ij}/2$; reading $B_{ij}$ as $\langle\hat a_i^\dagger\hat a_j\rangle+\delta_{ij}/2$ (that is, $\langle\hat X^i{}_j\rangle$) gives the transpose, which differs. The four matrix constraints are Correct. Checked by: exact WL computation with an exact rational complex symmetric $Z$: the second moments are solved exactly from the annihilator relation $\hat a_i|Z\rangle=-\sum_kZ_{ik}\hat a_k^\dagger|Z\rangle$ (QF exact in the previous entry) and the CCR, then compared with the formulas; the constraints hold with zero residual.
QF code (exact): none. `BosonicMatrixElement` takes polynomials between Fock states; it does not give moments of a multimode Gaussian state, and `CovarianceMatrix` works only on cut-off states.

<a id="fermionic-pair-exponential-state-normalization-and-overlap"></a>
### Fermionic pair-exponential state, normalization and overlap
Appears in: P-043

A fermionic pair-exponential state with skew $Z$ is normalized by $\det(I+ZZ^\dagger)^{-1/4}$, is annihilated by $\hat{\mathbf c}+Z\hat{\mathbf c}^\dagger$, and two such states overlap by $\det(I+WZ^\dagger)^{1/2}$ over the normalizations.

$$ |Z\rangle=\det(I_N+ZZ^\dagger)^{-1/4}e^{-\frac12\hat{\mathbf c}^{\dagger T}Z\hat{\mathbf c}^\dagger}|0\rangle,\ Z^T=-Z,\ (\hat{\mathbf c}+Z\hat{\mathbf c}^\dagger)|Z\rangle=0,\quad \langle Z|W\rangle=\frac{\det(I_N+WZ^\dagger)^{1/2}}{[\det(I_N+ZZ^\dagger)\det(I_N+WW^\dagger)]^{1/4}} $$

Audit: Correct. Checked by: exact fermion matrices (4 modes), exact rational complex skew $Z,W$: the squared norm squared equals $\det(I+ZZ^\dagger)$, the annihilator residual is exactly 0, and the unnormalized overlap equals the principal square root of $\det(I+WZ^\dagger)$ exactly (the atlas's branch caveat stands for other parameters).
QF code (exact): none. no fermionic layer.

<a id="fermionic-bilinear-symbols-and-purity-constraints"></a>
### Fermionic bilinear symbols and purity constraints
Appears in: P-045

On a pure fermionic Gaussian the skew pairing matrix $A$ and the number matrix $B$ are rational in $Z$ and obey the same quadratic constraints as for bosons with the signs changed.

$$ A=-\hbar Z(I_N+Z^\dagger Z)^{-1},\quad B=\hbar Z(I_N+Z^\dagger Z)^{-1}Z^\dagger-\tfrac\hbar2I_N,\quad B^2-A\bar A=\tfrac{\hbar^2}{4}I_N,\quad BA=A\bar B,\quad B^\dagger=B,\quad A^T=-A $$

Audit: Convention: on the four-mode state of the pair-exponential entry ($\hbar=1$), $A_{ij}=\langle\hat c_j\hat c_i\rangle=-\langle\hat c_i\hat c_j\rangle$ and $B_{ij}=\langle\hat c_j^\dagger\hat c_i\rangle-\delta_{ij}/2$; read as $\langle\hat X_{ij}\rangle$ with $\hat X_{ij}=\hat c_i\hat c_j$, the sign of $A$ is wrong. The constraints are Correct. Checked by: exact fermion matrices (4 modes), exact rational $Z$, all residuals exactly 0.
QF code (exact): none. no fermionic layer.

<a id="affine-gaussian-state-normalization-derived-correction"></a>
### Affine Gaussian state normalization (derived correction)
Appears in: P2-017

The state $e^{\zeta\hat a^{\dagger2}/2+\alpha\hat a^\dagger}|0\rangle$ is normalized by $(1-|\zeta|^2)^{1/4}$ times a Gaussian in $\alpha$ that includes the cross term $\mathrm{Re}(\bar\zeta\alpha^2)$.

$$ |\zeta,\alpha\rangle=Ne^{\zeta\hat a^{\dagger2}/2+\alpha\hat a^\dagger}|0\rangle,\qquad N=(1-|\zeta|^2)^{1/4}\exp\left[-\frac{|\alpha|^2+\operatorname{Re}(\bar\zeta\alpha^2)}{2(1-|\zeta|^2)}\right] $$

Audit: Correct as corrected; the printed $N=(1-|\zeta|^2)^{-1/2}e^{-|\alpha|^2/2}$ is wrong. Checked by: exact WL computation: the squared norm, the Bargmann Gaussian integral $\int d^2w/\pi\,e^{-|w|^2}|e^{\zeta w^2/2+\alpha w}|^2$ with symbolic complex $\zeta,\alpha$, times the corrected $N^2$ is exactly 1; with the printed $N$ it is $(1-\zeta^2)^{-3/2}e^{\alpha^2\zeta/(1-\zeta)}\ne1$ already for real $\zeta,\alpha$.
QF code (exact): none. `BosonicMatrixElement[{0, 0}, Exp[ca a] ** Exp[cz/2 a ** a] ** Exp[z/2 ad ** ad] ** Exp[al ad]]` returns unevaluated, and `BosonicExpOrder` of that product orders the quadratic factors but leaves the linear wings unordered.

<a id="fermionic-paired-vacuum-normalization"></a>
### Fermionic paired vacuum normalization
Appears in: P2-027

The fermionic paired vacuum $e^{\sum Z_{ij}\hat a_i^\dagger\hat b_j^\dagger}|0,0\rangle$ has squared norm $\det(I_N+Z^\dagger Z)$.

$$ |Z\rangle=\det(I_N+Z^\dagger Z)^{-1/2}e^{\sum_{ij}Z_{ij}\hat a_i^\dagger \hat b_j^\dagger}|0,0\rangle $$

Audit: Correct (the chart relation to $A,B,C,D$ is Unchecked). Checked by: exact fermion matrices (4 modes), exact rational complex $2\times2$ $Z$: squared norm minus $\det(I+Z^\dagger Z)$ is exactly 0.
QF code (exact): none. no fermionic layer.

<a id="fermionic-paired-vacuum-overlap"></a>
### Fermionic paired vacuum overlap
Appears in: P2-028

Two fermionic paired vacua overlap by $\det(I_N+Z_1^\dagger Z_2)$ over the normalizations.

$$ \langle Z_1|Z_2\rangle=\frac{\det(I_N+Z_1^\dagger Z_2)}{\det(I_N+Z_1^\dagger Z_1)^{1/2}\det(I_N+Z_2^\dagger Z_2)^{1/2}},\qquad Z'=(AZ+B)(CZ+D)^{-1} $$

Audit: Correct for the overlap; the group action is Unchecked (convention). Checked by: exact fermion matrices (4 modes), exact rational $Z_1,Z_2$: residual exactly 0.
QF code (exact): none. no fermionic layer.

<a id="fermionic-pair-number-generating-determinant"></a>
### Fermionic pair-number generating determinant
Appears in: P2-029

The fermionic pair-number generating function is $\det(I+tZ^\dagger Z)/\det(I+Z^\dagger Z)$, a product of Bernoulli laws.

$$ \mathcal F(t)=\sum_{n=0}^NW_nt^n=\frac{\det(I_N+tZ^\dagger Z)}{\det(I_N+Z^\dagger Z)},\qquad p_j=\frac{\lambda_j}{1+\lambda_j} $$

Audit: Correct. Checked by: exact fermion matrices (4 modes), $t^{\hat N_a}=\prod_i(1+(t-1)\hat n_{a_i})$ with symbolic $t$: residual exactly 0 as a polynomial in $t$. Bernoulli factorization follows from the eigenvalues of $Z^\dagger Z$.
QF code (exact): none. no fermionic layer.

<a id="squeezed-vacuum-from-the-bargmann-riccati-equation"></a>
### Squeezed vacuum from the Bargmann Riccati equation
Appears in: C-027

The squeezed vacuum is the Gaussian e^{a^dag W a^dag/2} on the vacuum with normalization det(I - W^*W)^{1/4}, obtained from a matrix Riccati equation in Bargmann space; its vacuum amplitude is det(cosh|B|)^{-1/2}.

$$ \dot M=B-MB^*M,\ M_0=0,\quad \dot a=-\tfrac12a\operatorname{Tr}(B^*M),\quad \hat S(B)|0\rangle=\det(I_n-W^*W)^{1/4}\exp\!\left(\tfrac12\hat a^\dagger\cdot W\hat a^\dagger\right)|0\rangle,\quad \langle0|\hat S(B)|0\rangle=\det(\cosh|B|)^{-1/2} $$

Audit: Correct. Forcing line for the Riccati: the generator acting on a e^{zeta M zeta/2} gives (1/2)[zeta B zeta - Tr(B^* M) - zeta M B^* M zeta] times the Gaussian. Several modes: with the Takagi form B = V D V^T (V unitary, D >= 0), M_t = V tanh(tD) V^T solves the Riccati equation because V^T conj(V) = I reduces it to D sech^2(tD) = D - tanh(tD) D tanh(tD), and the amplitude integrates to det(cosh(tD))^{-1/2}. Checked by: QF exact (one mode) and hand derivation (n modes).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{SqueezeOperator[xi, Infinity, "Ordering" -> "Normal"],
 BosonicMatrixElement[{0, 0}, SqueezeOperator[xi]],
 Simplify[BosonicMatrixElement[{2, 0}, SqueezeOperator[xi]] - (1 - Abs[w]^2)^(1/4) w/Sqrt[2] /. w -> -Exp[I Arg[xi]] Tanh[Abs[xi]]]}
```
Output: `{E^(-1/2*(E^(I*Arg[xi])*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]1]*Tanh[Abs[xi]]))**Sech[Abs[xi]]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^((\[FormalA]1**\[FormalA]1*Tanh[Abs[xi]])/(2*E^(I*Arg[xi]))), Sqrt[Sech[Abs[xi]]], 0}`. Matches after xi -> -xi (QF convention): the normal-ordered squeezer is the atlas Gaussian, <0|S|0> = Sech[|xi|]^(1/2) = det(cosh|B|)^{-1/2}, and <2|S|0> equals the atlas value (difference 0).

<a id="coordinate-form-of-the-squeezed-vacuum"></a>
### Coordinate form of the squeezed vacuum
Appears in: C-028; one-mode case in C-049 (delta = (1-omega)/(1+omega))

In position space the squeezed vacuum is a normalized complex Gaussian whose covariance is the Cayley transform of W.

$$ \Gamma=i(I_n-W)(I_n+W)^{-1},\qquad \psi_B(x)=(\pi\hbar)^{-n/4}\det(I_n-W^*W)^{1/4}\det(I_n+W)^{-1/2}e^{ix^T\Gamma x/(2\hbar)} $$

Audit: Correct (the atlas's corrected amplitude). One mode, W = tanh r: the generator is -r(x d/dx + 1/2), so S|0> = pi^{-1/4} e^{-r/2} e^{-e^{-2r}x^2/2}, and (1 - tanh^2)^{1/4}(1 + tanh)^{-1/2} = e^{-r/2}; complex W: the norm integral gives 1. Checked by: hand derivation.
QF code (exact): none. A position-space wavefunction; QF has no exact position-representation call for bosonic states.

<a id="reduced-gaussian-channel-from-a-joint-quadratic-evolution"></a>
### Reduced Gaussian channel from a joint quadratic evolution
Appears in: C-056

A joint quadratic evolution with a Gaussian environment acts on the system by a linear map of the characteristic-function argument times a Gaussian noise factor.

$$ \chi_{S,\rm out}(k)=\chi_{S,\rm in}(X^Tk)\exp[-\tfrac12k^TYV_EY^Tk],\qquad V_{S,\rm out}=XV_SX^T+YV_EY^T $$

Audit: Correct. Forcing line: e^{ik.(Xz+Yu)} factors into system and environment parts, and the centered Gaussian environment has chi_E(m) = exp(-m^T V_E m/2). Checked by: hand derivation.
QF code (exact): none. QF has no characteristic-function calculus for bosonic channels.

<a id="quadratic-propagators-and-riccati"></a>
## Quadratic propagators and Riccati

<a id="general-quadratic-expression"></a>
### General quadratic expression
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

<a id="classical-generator-and-heisenberg-action-of-a-quadratic-hamiltonian"></a>
### Classical generator and Heisenberg action of a quadratic Hamiltonian
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

<a id="quadratic-propagator-vacuum-coefficient-corrected-6-23"></a>
### Quadratic propagator vacuum coefficient (corrected (6.23))
Appears in: B-071 (B-077 cites the same scalar terms)

The vacuum amplitude of e^{itH} is det(Phi e^{itC})^{-1/2} times the exponential of a time integral of terms quadratic and linear in the shift g, the quadratic one with coefficient 1/2.

$$ c(t)=\det(\hat\Phi e^{it\hat C})^{-1/2} \exp\!\left[i\int_0^t\left(\tfrac12g\hat\Phi'^{-1}\bar{\hat A}\hat\Phi^{-1}g-f^*\hat\Phi^{-1}g\right)d\tau\right] $$

Conditions: square root continued from c(0) = 1; constant k adds e^{itk}.
Audit: Correct as corrected in the atlas (the printed source omits the 1/2). Checked by: exact WL computation: in the Bargmann realization (a acts as d/dz, a* as z), U(t)|0> = c(t) exp(z R z/2 + l z) with R = -Phi^{-1} Psi, l = -Phi^{-1} g from exact e^{itL}; the residual of d/dt - iH on this vector vanishes identically (one mode, A = 3/10 + 2i/5, C = 1, f = 1/2 + i/3, symbolic t: residual 0; two modes, complex A, Hermitian C, complex f at t = 7/10 with 60-digit arithmetic: residual 0 to 58 digits); the printed coefficient (no 1/2) leaves a constant residual of about 0.03 in both cases. Forcing line: the constant term of (1/2) conj(A) d^2/dz^2 acting on e^{lz} is (1/2) conj(A) l^2.
QF code (exact): none. `BosonicExpOrder[Exp[I t H]]` returns its input unchanged for a quadratic H with linear terms, and `BosonicExpSimplify[Exp[I t H] ** a ** Exp[-I t H]]` stays unevaluated, so QF gives no exact vacuum amplitude for e^{itH}.

<a id="logarithmic-derivative-of-a-determinant"></a>
### Logarithmic derivative of a determinant
Appears in: B-073

The time derivative of log det M is the trace of M^{-1} times the derivative of M (Jacobi's formula).

$$ \frac{d}{dt}\log\det \hat M=\operatorname{Tr}(\hat M^{-1}\dot{\hat M}) $$

Audit: Correct (Jacobi's formula; determinant-class assumptions as stated). Checked by: hand derivation: log det M = Tr log M, differentiated under the trace.
QF code (exact): none. A matrix identity; no second-quantization call is involved.

<a id="riccati-fixed-point-for-the-diagonalizing-transformation"></a>
### Riccati fixed point for the diagonalizing transformation
Appears in: B-085

The ground state of H = (q beta q + p gamma p)/2 is the Gaussian exp(i q Z q/2), with Z = iW and W the positive solution of W gamma W = beta, given in closed form by the displayed square roots.

$$ \hat Z=i\hat\gamma^{-1/2}(\hat\gamma^{1/2}\hat\beta\hat\gamma^{1/2})^{1/2}\hat\gamma^{-1/2},\qquad \hat\beta=\hat S+\hat T,\quad\hat\gamma=\hat S-\hat T $$

Audit: Correct. Checked by: hand derivation (W gamma W = gamma^{-1/2} X^{1/2} X^{1/2} gamma^{-1/2} = beta for X = gamma^{1/2} beta gamma^{1/2}) and exact WL computation (exact rational 2x2 positive beta, gamma: W gamma W - beta and W - W^T reduce to 0 with `RootReduce`, and W has positive eigenvalues).
QF code (exact): none. A real matrix equation; no second-quantization call is involved.

<a id="riccati-equation-and-phase-for-squeezing-dynamics-phase-sign-corrected"></a>
### Riccati equation and phase for squeezing dynamics (phase sign corrected)
Appears in: P2-013

Under a Hamiltonian linear in the su(1,1) generators, a coherent state stays coherent; its label obeys a Riccati equation and its phase obeys a first-order equation.

$$ \hat H=i(\beta \hat K_+-\bar\beta \hat K_--i\gamma \hat K_0),\ |\psi\rangle=e^{-i\phi}|\zeta;k\rangle:\quad \dot\zeta=\beta-\bar\beta\zeta^2-i\gamma\zeta,\qquad i\dot\phi=k(\bar\beta\zeta-\beta\bar\zeta+i\gamma) $$

Audit: Correct as corrected in the atlas; the printed phase sign is wrong. Checked by: exact WL computation in the holomorphic realization, $\psi=e^{-i\phi}(1-\zeta\bar\zeta)^k(1-\zeta z)^{-2k}$ with $\zeta,\bar\zeta,\phi$ functions of $t$: the Schrodinger residual divided by $\psi$ is exactly 0 with the corrected phase, and $2ik(\bar\beta\zeta-\beta\bar\zeta)=-4k\,\mathrm{Im}(\bar\beta\zeta)$ with the printed one.
QF code (exact): none. No su(1,1) layer in QF.

<a id="floquet-classification-of-a-quadratic-oscillator"></a>
### Floquet classification of a quadratic oscillator
Appears in: P2-014

The sign of the Lorentz norm of the generator classifies a periodic quadratic oscillator: discrete spectrum, two-sided continuum, or one-sided continuum.

$$ \Omega_0^2-\Omega_1^2-\Omega_2^2>0:\ \text{spectrum } \pm\hbar\Omega(k+n);\quad <0:\ \text{two-sided continuous};\quad =0,\ \Omega\ne0:\ \text{one-sided continuous} $$

Audit: Correct as a statement about Lie-algebra generators (timelike elements are conjugate to $\Omega\hat K_0$, whose spectrum is $\Omega(k+n)$); the existence of a real logarithm of the period map is a separate assumption, as the atlas says. Checked by: hand derivation for the timelike case; other cases Unchecked.
QF code (exact): none. QF has no Floquet classification.

<a id="multimode-pair-riccati-equation-and-vacuum-phase"></a>
### Multimode pair Riccati equation and vacuum phase
Appears in: P2-024

The pair matrix of a multimode paired vacuum obeys a matrix Riccati equation, and its vacuum phase obeys a trace equation.

$$ i\dot Z=A+ZA^\dagger Z+CZ+ZD^T,\qquad \dot\phi=\tfrac12\operatorname{Tr}(A^\dagger Z+Z^\dagger A) $$

Conditions: $\hat H=A_{ij}\hat a_i^\dagger\hat b_j^\dagger+\bar A_{ij}\hat a_i\hat b_j+C_{ij}\hat a_i^\dagger\hat a_j+D_{ij}\hat b_i^\dagger\hat b_j$.
Audit: Correct. Checked by: hand derivation. The pair-creation coefficient gives $[C_{ij}\hat a_i^\dagger\hat a_j,\,Z_{kl}\hat a_k^\dagger\hat b_l^\dagger]\to(CZ)$, the $D$ term gives $ZD^T$, and the double commutator of $\bar A_{ij}\hat a_i\hat b_j$ gives $ZA^\dagger Z$. The vacuum coefficient gives $\dot\phi+i\dot N/N=\operatorname{Tr}(A^\dagger Z)$; its real part is the stated phase, and its imaginary part reproduces $\frac{d}{dt}\tfrac12\log\det(I-Z^\dagger Z)$.
QF code (exact): none. QF has no time-dependent Gaussian-state equations.

<a id="mehler-kernel-and-euclidean-oscillator-kernel"></a>
### Mehler kernel and Euclidean oscillator kernel
Appears in: C-012

The sum of products of Hermite functions weighted by w^k has a closed Gaussian form; with w = e^{-tau} and the zero-point factor it becomes the imaginary-time oscillator kernel.

$$ \sum_{k}w^{|k|}\phi_k(x)\phi_k(y)=\pi^{-n/2}(1-w^2)^{-n/2}\exp\!\left[-\frac{1+w^2}{2(1-w^2)}(|x|^2+|y|^2)+\frac{2w\,x\cdot y}{1-w^2}\right],\quad |w|<1 $$

and, n = 1, kernel of e^{-tau H_os} = (2 pi sinh tau)^{-1/2} exp[-(x^2+y^2)/(2 tanh tau) + xy/sinh tau].

Audit: Correct. Checked by: exact WL computation (n = 1, symbolic x, y: the Taylor coefficients of the closed form in w agree with the Hermite-function sum through order w^10; with w = e^{-tau} and the factor e^{-tau/2} the closed form divided by the tau kernel simplifies to 1); several dimensions factor into one-dimensional kernels.
QF code (exact): none. A position-space kernel; QF has no exact position-representation call for bosonic operators.

<a id="linear-quadratic-flow-and-heisenberg-action"></a>
### Linear quadratic flow and Heisenberg action
Appears in: C-018 (constant-S case also in C-023, C-030)

For a quadratic Hamiltonian the Heisenberg evolution of (Q, P) is linear and given by the classical flow matrix solving F' = J S F, which is e^{tJS} for constant S.

$$ \dot F_t=JS_tF_t,\ F_0=I_{2n},\quad \hat U_t^*(\hat Q,\hat P)^T\hat U_t=F_t(\hat Q,\hat P)^T,\quad F_t=e^{tJS}\ (S\ \text{constant}) $$

Conditions: H = z^T S z / 2 quantized symmetrically with the factor 1/2.
Audit: Correct. Checked by: QF exact (one mode, symbolic S, hbar = 1): i[H, (Q, P)] = J S (Q, P), and the full conjugation U^* (Q, P) U equals e^{tJS} (Q, P) in closed form.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2]; p = -I (a - ad)/Sqrt[2];
h = (s11 q ** q + s12 (q ** p + p ** q) + s22 p ** p)/2;
qp[x_] := Simplify[{Coefficient[x, a] + Coefficient[x, ad], I (Coefficient[x, a] - Coefficient[x, ad])}/Sqrt[2]];
gen = qp /@ {BosonicNormalOrder[I (h ** q - q ** h), {a, ad}], BosonicNormalOrder[I (h ** p - p ** h), {a, ad}]};
flow = qp /@ {BosonicExpSimplify[Exp[I t h] ** q ** Exp[-I t h]], BosonicExpSimplify[Exp[I t h] ** p ** Exp[-I t h]]};
{gen, gen - {{0, 1}, {-1, 0}} . {{s11, s12}, {s12, s22}},
 Simplify[flow - MatrixExp[t {{0, 1}, {-1, 0}} . {{s11, s12}, {s12, s22}}], t > 0 && s12^2 - s11 s22 > 0]}
```
Output: `{{{s12, s22}, {-s11, -s12}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}}`. Matches: the generator is J S and the flow difference vanishes (Simplify run in the case s12^2 > s11 s22; both sides are entire in the entries of S).

<a id="gaussian-propagation-by-a-linear-fractional-map"></a>
### Gaussian propagation by a linear-fractional map
Appears in: C-019

A Gaussian e^{ix^T Gamma x/(2 hbar)} stays Gaussian under a quadratic propagator: its covariance moves by the linear-fractional map of the flow blocks and its amplitude picks up det(A_t + B_t Gamma)^{-1/2}.

$$ \Gamma_t=(C_t+D_t\Gamma)(A_t+B_t\Gamma)^{-1},\qquad \hat U_t\big[a_\Gamma e^{ix^T\Gamma x/(2\hbar)}\big]=a_\Gamma\det(A_t+B_t\Gamma)^{-1/2}e^{ix^T\Gamma_tx/(2\hbar)} $$

Conditions: Gamma symmetric, Im Gamma > 0; root continued from t = 0.
Audit: Correct. Oscillator: Gamma = i stays i and det^{-1/2} = e^{-it/2} (zero-point phase); free particle: Gamma_t = i/(1+it), the known spreading packet. Checked by: hand derivation.
QF code (exact): none. A position-space Gaussian statement; QF has no exact position-representation call for bosonic states.

<a id="siegel-domain-invariance-and-determinant-cocycle"></a>
### Siegel-domain invariance and determinant cocycle
Appears in: C-020

The linear-fractional action of a symplectic matrix keeps a complex symmetric matrix with positive imaginary part inside that domain, composes as a group action, and its determinant factor is a cocycle.

$$ \Im\Gamma_F=(A+B\Gamma)^{-*}(\Im\Gamma)(A+B\Gamma)^{-1},\quad \alpha(F_1F_2)=\alpha(F_1)\alpha(F_2),\quad \delta(F_1F_2,\Gamma)=\delta(F_1,\alpha(F_2)\Gamma)\delta(F_2,\Gamma) $$

Audit: Correct. Checked by: WL computation with finite 4-by-4 symplectic matrices (two degrees of freedom, F1 and F2 as matrix exponentials, complex Gamma: all three identities hold to machine precision; Gamma_F is symmetric with positive imaginary part).
QF code (exact): none. Classical matrix algebra; QF's second-quantization layer does not apply.

<a id="matrix-riccati-and-amplitude-equations"></a>
### Matrix Riccati and amplitude equations
Appears in: C-021

The complex covariance of a propagated Gaussian obeys a matrix Riccati equation and its amplitude a linear equation with a trace.

$$ \dot\Gamma=-K-L\Gamma-\Gamma L^T-\Gamma G\Gamma,\qquad \dot a=-\tfrac12\operatorname{Tr}(L+G\Gamma)a,\qquad S=\begin{pmatrix}K&L\\L^T&G\end{pmatrix} $$

Audit: Correct. Forcing line: with M = A + B Gamma0, N = C + D Gamma0, d/dt(M;N) = JS(M;N) = (L^T M + G N; -K M - L N), and Gamma = N M^{-1} gives the Riccati equation; a ~ det(M)^{-1/2} gives the trace. Checked by: hand derivation, plus a WL computation with finite matrices (central-difference derivative of the linear-fractional map agrees to finite-difference accuracy).
QF code (exact): none. Classical matrix algebra; QF's second-quantization layer does not apply.

<a id="coherent-matrix-elements-of-a-metaplectic-operator"></a>
### Coherent matrix elements of a metaplectic operator
Appears in: C-022

The coherent-state matrix element of the propagator of a quadratic Hamiltonian is a closed Gaussian in the two labels, with det(Y)^{-1/2}, Y = A + D + i(B - C), as prefactor.

$$ \langle\phi_z,\hat U\phi_{z'}\rangle=2^{n/2}\det(Y)^{-1/2}e^{-i\sigma(Fz',z)/2}e^{-(|x|^2+|\xi|^2)/4}e^{-(\xi+ix)^TW(\xi+ix)/4},\quad Y=A+D+i(B-C),\ W=ZY^{-1},\ (x,\xi)=z-Fz' $$

Conditions: hbar = 1, Z = A - D + i(B+C), root along the path.
Audit: Correct. Checked by: exact WL computation (one mode, S = {{9/10, 3/10}, {3/10, 3/5}}, t = 4/5, hbar = 1: H written in su(1,1) generators, checked by QF normal ordering; U = e^{-itH} disentangled as e^{c K+} e^{d K0} e^{f K-} from the exact 2-by-2 representation; its coherent elements written in closed form and compared with the atlas formula at four label pairs at 40 digits: equal, while dropping the sigma phase leaves an order-one residual; the 2-by-2 entry whose square root is taken stays away from the negative axis along the path, so the principal root is the continuous one).
QF code (exact): none. `BosonicExpOrder[Exp[x1 a**a + x2 ad**ad + x3 ad**a]]` returns the exponential unchanged, and `BosonicMatrixElement[..., "Basis" -> "Coherent"]` leaves exponentials unevaluated.

<a id="squeeze-rotation-factorization-and-its-vacuum-phase"></a>
### Squeeze-rotation factorization and its vacuum phase
Appears in: C-029

A quadratic propagator factors into a squeeze and a number-preserving rotation; the rotation multiplies the vacuum by a unit-modulus phase fixed by det(Y_t) and det(u_t u_t^dag).

$$ \hat U_t=\hat S(B_t)\hat R_t,\quad \hat R_t|0\rangle=\vartheta_t|0\rangle,\quad \vartheta_t=\frac{2^{n/2}\det(Y_t)^{-1/2}}{\det(u_tu_t^\dagger)^{-1/4}},\quad |\vartheta_t|=1,\quad Y_t=2\bar u_t $$

Audit: Correct. |det Y|^{-1/2} 2^{n/2} = |det u|^{-1/2} = det(u u^dag)^{-1/4}, so |vartheta| = 1; Y = 2 conj(u) follows from C-024. Checked by: hand derivation (consistency of C-022, C-024, C-027).
QF code (exact): none. QF's exact layer has no metaplectic factorization; `QuantumInterferometer` works on cut-off Fock spaces and is not used.

<a id="herman-kluk-representation-exact-for-quadratics"></a>
### Herman-Kluk representation, exact for quadratics
Appears in: C-054

For a quadratic Hamiltonian, the phase-space integral of |phi_{F_t z}><phi_z| with the prefactor det(u_t)^{1/2} reproduces the propagator exactly.

$$ \hat U_t=(2\pi\hbar)^{-n}\int C_t\,|\phi_{F_tz}\rangle\langle\phi_z|\,dz,\qquad C_t=\det\!\left[\frac{A_t+D_t+i(C_t^{\rm block}-B_t)}2\right]^{1/2} $$

Audit: Correct (the atlas's normalized re-expression; C_t = det u^{1/2}). Checked by: exact WL computation (the C-022 data: the coherent element of the Herman-Kluk integral is a two-dimensional Gaussian integral, done in closed form with positive-definite real part, and it equals the exact disentangled matrix element at four label pairs at 40 digits; coherent elements determine the operator).
QF code (exact): none. `BosonicMatrixElement[..., "Basis" -> "Coherent"]` leaves exponentials unevaluated and QF has no phase-space integration.

<a id="forced-and-parametric-oscillators"></a>
## Forced and parametric oscillators

<a id="forced-oscillator-completing-the-square"></a>
### Forced oscillator: completing the square
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

<a id="regularized-forced-hamiltonian"></a>
### Regularized forced Hamiltonian
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

<a id="forced-oscillator-displacement-plus-accumulated-phase"></a>
### Forced oscillator: displacement plus accumulated phase
Appears in: P2-010; G-022

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

**Also stated in G-022 (Linearly driven evolution equals a displacement up to phase):**

A linearly driven oscillator, built as a product of short displacements, is one displacement by the accumulated amplitude times an accumulated phase.

$$ \hat U(T,0)=e^{i\Phi(T)}\hat D(\alpha(T)),\qquad\alpha(T)=\int_0^T u(t)\,dt. $$

Audit: Correct (structure). Checked by: QF exact on a product of three displacements, and hand derivation from the Weyl law; the time-dependent drive u is book-specific and not rederived.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
BosonicExpSimplify[DisplacementOperator[x1, Infinity] ** DisplacementOperator[x2, Infinity] ** DisplacementOperator[x3, Infinity]]
```
Output: `E^(-(\[FormalA]1*Conjugate[x1]) - \[FormalA]1*Conjugate[x2] + (-(x2*Conjugate[x1]) + x1*Conjugate[x2])/2 - \[FormalA]1*Conjugate[x3] + (-(x3*(Conjugate[x1] + Conjugate[x2])) + (x1 + x2)*Conjugate[x3])/2 + x1*SuperDagger[\[FormalA]1] + x2*SuperDagger[\[FormalA]1] + x3*SuperDagger[\[FormalA]1])`. Matches: the product is a single displacement by x1 + x2 + x3 times a scalar phase built from the pairwise Weyl phases.

<a id="parametric-oscillator-su-1-1-generators"></a>
### Parametric oscillator su(1,1) generators
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

<a id="singular-oscillator-ground-state-and-su-1-1-label"></a>
### Singular oscillator ground state and su(1,1) label
Appears in: P2-015

The half-line inverse-square oscillator has ground state $x^ae^{-x^2/2}$ with $a(a-1)=2g^2$, energy $a+1/2$, and an equally spaced ladder with su(1,1) label $k=E_0/2$.

$$ \psi_0\sim x^a e^{-x^2/2},\quad a(a-1)=2g^2,\quad E_0=a+\tfrac12,\quad k=\tfrac12(a+\tfrac12),\quad E_n=2n+E_0 $$

Conditions: $\hat H=\tfrac12(\hat p^2+x^2)+g^2/x^2$ on the half line, regular boundary condition.
Audit: Correct. Checked by: exact WL computation, $\hat H(x^ae^{-x^2/2})/(x^ae^{-x^2/2})$ simplifies to $a+1/2$ when the inverse-square coefficient is $a(a-1)/2$; $k=E_0/2$ since $\hat K_0=\hat H/2$.
QF code (exact): none. Not a Fock-space operator identity.

<a id="forced-parametric-oscillator-generating-kernel-and-vacuum-amplitude"></a>
### Forced parametric oscillator generating kernel and vacuum amplitude
Appears in: P2-018

For a forced parametric oscillator the vacuum amplitude has modulus $(1-\rho)^{1/4}e^{-\nu(1-\sqrt\rho\cos2\varphi)/2}$.

$$ |\langle0|\psi\rangle|=(1-\rho)^{1/4}e^{-\nu(1-\sqrt\rho\cos2\varphi)/2} $$

Audit: Unchecked: the map from the force and frequency history to $(\rho,\nu,\varphi)$ is a source lookup. At $\nu=0$ it reduces to the squeezed-vacuum amplitude $(1-|\zeta|^2)^{1/4}$ of the affine Gaussian entry, consistent with $\rho=|\zeta|^2$.
QF code (exact): none. The source parameters are not reconstructed.

<a id="drive-equation-and-displacement-channel-factor"></a>
### Drive equation and displacement channel factor
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

<a id="propagator-of-a-driven-quadratic-hamiltonian"></a>
### Propagator of a driven quadratic Hamiltonian
Appears in: C-032

For a Hamiltonian that is quadratic plus linear, the propagator is the classical action phase times a translation along the classical orbit, the quadratic propagator, and the inverse translation at the start time.

$$ \delta_t=\int_0^t[p_s\dot q_s-H(z_s,s)]ds-\tfrac12(q_tp_t-qp),\qquad \hat U_2(t,s)=e^{i(\delta_t-\delta_s)/\hbar}\hat T(z_t)\hat S(B_t)\hat R_t\hat R_s^{-1}\hat S(-B_s)\hat T(-z_s) $$

Audit: Correct for an affine quadratic Hamiltonian. Checked by: exact WL computation (H = N + 1/2 + f Q with symbolic f, q, p, t, hbar = 1, B = 0, R_t = e^{-it(N+1/2)}, s = 0: in the Bargmann realization a -> d/dz, a^dag -> z, the right-hand side applied to an arbitrary analytic function solves i d_t psi = H psi identically, FullSimplify giving 0, with the classical orbit and delta_t computed symbolically and the product equal to the identity at t = 0); general s follows from U(t,s) = U(t,0) U(s,0)^{-1}.
QF code (exact): none. `BosonicExpOrder[Exp[-I t (ad**a + 1/2 + f (a + ad)/Sqrt[2])]]` returns the exponential unfactored.

<a id="periodic-quadratic-evolution-and-quasi-energies"></a>
### Periodic quadratic evolution and quasi-energies
Appears in: C-035

For a periodic quadratic Hamiltonian in a stable region, the one-period propagator has eigenvalues of the form e^{-i rho pi (k+1/2)}.

$$ \text{Floquet eigenvalues } e^{-i\rho\pi(k+1/2)} $$

Audit: Unchecked: needs Mathieu-equation stability data; only the constant-frequency limit (rho = omega, period pi) was checked by hand.
QF code (exact): none. QF has no Floquet or classical-stability routine.

<a id="time-dependent-su-1-1-hamiltonian"></a>
### Time-dependent SU(1,1) Hamiltonian
Appears in: C-050

For a Hamiltonian linear in the su(1,1) generators with time-dependent coefficients, the propagator is a squeeze, a K0 phase rotation and an inverse squeeze, with the squeeze coordinate solving a scalar Riccati equation.

$$ i\dot\omega=\bar\lambda\omega^2+\mu\omega+\lambda,\ \dot\gamma=-\lambda\bar\omega-\bar\lambda\omega-\mu,\qquad \hat U(t,s)=\hat S(\beta_t)e^{i(\gamma_t-\gamma_s)\hat K_0}\hat S(-\beta_s),\ \omega_t=\frac{\beta_t}{|\beta_t|}\tanh|\beta_t| $$

Audit: Correct (the atlas's corrected endpoint). Checked by: WL computation in the exact 2-by-2 representation of su(1,1) (constant lambda, mu; Riccati integrated numerically: U(t,0) = S(beta_t) e^{i gamma_t K0} agrees at three times to the integration accuracy; the corrected two-time form agrees; the printed form with S(-beta_t) fails by an order-one residual).
QF code (exact): none. QF has no time-dependent disentangling routine.

<a id="degenerate-parametric-amplifier-gives-a-squeeze-of-parameter-2-eta-t"></a>
### Degenerate parametric amplifier gives a squeeze of parameter 2 eta t
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

<a id="nondegenerate-parametric-amplifier-gives-a-two-mode-squeeze-of-parameter-eta-t"></a>
### Nondegenerate parametric amplifier gives a two-mode squeeze of parameter eta t
Appears in: G-075

Under an undepleted pump the nondegenerate amplifier generates a two-mode squeeze with parameter eta t, with no factor 2.

$$ \hat H_I=i\hbar(\eta^*\hat a\hat b-\eta \hat a^\dagger \hat b^\dagger),\quad \hat U_I=\hat S_2(\eta t). $$

Audit: Correct. Checked by: hand derivation (-i H_I t/hbar = t(eta^* a b - eta a^dagger b^dagger), no factor 1/2 in S_2).
QF code (exact): none. the statement is the one-line exponent identity; the exponential it names is the two-mode squeeze checked in the G-068 and G-070 entries.

<a id="driven-damped-oscillator-coherent-amplitude-and-steady-state"></a>
### Driven damped oscillator: coherent amplitude and steady state
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

<a id="normal-symbols-and-normal-order-exponentials"></a>
## Normal symbols and normal-order exponentials

<a id="number-conserving-exponential-in-normal-order"></a>
### Number-conserving exponential in normal order
Appears in: X-007 (source records cited: B-013, B-028, B-040, B-107, G-009, G-097); B-107 (special cases in B-013, B-028, B-105); G-097

The exponential of a number-conserving quadratic equals the normal-ordered exponential of the same form with $h$ replaced by $e^h-I_N$.

$$ e^{\hat a^\dagger h \hat a}=:\!\exp[\hat a^\dagger(e^h-I_N)\hat a]\!: $$

Conditions: finite modes; the same identity for fermions with CAR normal ordering.
Audit: Correct for both statistics. Checked by: QF exact (single mode); exact WL computation for two bosonic modes in the Bargmann realization, where $\hat a^\dagger h\hat a\to\sum_{ij}z_ih_{ij}\partial_j$: the normal-ordered exponential acts on a generic cubic as $f(z)\mapsto f(z+B^Tz)$ exactly ($B=e^h-I$), and $f(e^{th^T}z)$ solves the flow equation of $e^{t\hat a^\dagger h\hat a}$ (residual zero at a random numeric point, symbolic $2\times2$ matrix exponential), so both sides send $f(z)$ to $f(e^{h^T}z)$; exact fermion matrices (3 modes, exact $h=P\,\mathrm{diag}(\log2,\log3,\log\frac15)P^{-1}$ with complex $P$, normal order $c^\dagger_{i_1}\cdots c^\dagger_{i_k}c_{j_k}\cdots c_{j_1}$): difference exactly zero.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
BosonicExpOrder[Exp[lam ad ** a]]
```
Output: `NormalOrdered[E^((-1 + E^lam)*SuperDagger[\[FormalA]1]**\[FormalA]1)]`. Matches (single mode).

**Also stated in B-107 (Number-type exponential in normal order):**

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

**Also stated in G-097 (Normal ordering of a number exponential and the vacuum projector):**

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

<a id="vacuum-projector-as-a-normal-ordered-exponential"></a>
### Vacuum projector as a normal-ordered exponential
Appears in: X-007; B-013 (related: B-107 normal-order identity)

Normal ordering $e^{-\hat N}$ gives the projector onto the vacuum; for fermions the projector is the product of the single-mode empty projectors.

$$ |0\rangle\langle0|=:\!e^{-\sum_i\hat N_i}\!:\quad\text{(bosons)},\qquad |0\rangle\langle0|=\prod_i(\hat I-\hat N_i)\quad\text{(fermions)} $$

Audit: Correct; for fermions $:\!e^{-\sum_i\hat N_i}\!:$ also equals $\prod_i(\hat I-\hat N_i)$. Checked by: QF exact (single boson mode, all Fock indices symbolic, and coherent-state matrix elements; several modes follow because the normal-ordered exponential factorizes over modes) and exact fermion matrices (3 modes, both forms exactly $|0\rangle\langle0|$).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicMatrixElement[{m, n}, NormalOrdered[Exp[-ad ** a]]], BosonicMatrixElement[{z, w}, NormalOrdered[Exp[-ad ** a]], "Basis" -> "Coherent"], NormalOrdered[Exp[-ad ** a]]["Series"]}
```
Output: `{KroneckerDelta[0, n]*KroneckerDelta[m, n], E^((-Abs[w]^2 - Abs[z]^2)/2), Inactive[Sum][(-1)^k a†^k ** a^k / k!, {k, 0, Infinity}]}` (series index shortened). Matches: $\langle m|\cdot|n\rangle=\delta_{m0}\delta_{n0}$, and between normalized coherent states it is $\langle z|0\rangle\langle0|w\rangle$.

**Also stated in B-013 (Vacuum projector as a normal-ordered exponential):**

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

<a id="matrix-kernel-versus-normal-symbol"></a>
### Matrix kernel versus normal symbol
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

<a id="symbols-as-coherent-matrix-elements"></a>
### Symbols as coherent matrix elements
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

<a id="differential-realization-on-normal-symbols"></a>
### Differential realization on normal symbols
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

<a id="rank-one-projector-symbols"></a>
### Rank-one projector symbols
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

<a id="kernel-acting-on-a-vector"></a>
### Kernel acting on a vector
Appears in: B-026

An operator acts on a vector functional by integrating its matrix kernel against the vector with the Bargmann weight e^{-alpha* alpha}.

$$ (\hat A\Phi)(a^*)=\int K_{\hat A}(a^*,\alpha)\Phi(\alpha^*)e^{-\alpha^*\alpha}\Pi d\alpha^*d\alpha $$

Audit: Correct. Checked by: hand derivation: insert the coherent resolution of identity int |alpha><alpha| e^{-alpha* alpha} d^2alpha/pi = I between A and Phi.
QF code (exact): none. QF has no Bargmann-integral representation of operator action; `BosonicMatrixElement` gives matrix elements only.

<a id="product-of-operators-matrix-kernels"></a>
### Product of operators: matrix kernels
Appears in: B-027

The matrix kernel of a product is the Bargmann-weighted integral of the two kernels over the intermediate coherent label.

$$ K_{\hat A\hat B}(a^*,a)=\int K_{\hat A}(a^*,\alpha)K_{\hat B}(\alpha^*,a)e^{-\alpha^*\alpha}\Pi d\alpha^*d\alpha $$

Audit: Correct. Checked by: hand derivation: the same resolution of identity as B-026, inserted between A and B.
QF code (exact): none. QF has no kernel-composition integral; the normal-symbol version is checked in the next entry.

<a id="product-of-operators-normal-symbols"></a>
### Product of operators: normal symbols
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

<a id="resolvent-kernel-of-a-number-type-hamiltonian"></a>
### Resolvent kernel of a number-type Hamiltonian
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

<a id="covariant-q-and-contravariant-p-symbols-from-ordered-expansions"></a>
### Covariant (Q) and contravariant (P) symbols from ordered expansions
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

<a id="off-diagonal-covariant-symbol-by-analytic-continuation"></a>
### Off-diagonal covariant symbol by analytic continuation
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

<a id="disk-star-product-to-first-order-in-h"></a>
### Disk star product to first order in h
Appears in: P2-006

On the disk, the Berezin transform is the identity plus $h$ times the invariant Laplacian, and the symbol product differs from the pointwise product at order $h$ by a derivative pairing.

$$ \hat T_hf=f+h\hat\Delta f+o(h),\ \hat\Delta=(1-|z|^2)^2\partial_z\partial_{\bar z},\qquad A_1\star A_2=A_1A_2+h(1-|z|^2)^2(\partial_{\bar z}A_1)(\partial_zA_2)+o(h) $$

Audit: Unchecked: asymptotic statement not re-derived. The pairing $\partial_{\bar z}A_1\,\partial_zA_2$ is the one expected when $z$ labels the creation side.
QF code (exact): none. QF has no disk coherent states.

<a id="normal-ordered-powers-as-falling-factorials"></a>
### Normal-ordered powers as falling factorials
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

<a id="normally-ordered-coherent-matrix-elements"></a>
### Normally ordered coherent matrix elements
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

<a id="normally-ordered-correlation-functions"></a>
### Normally ordered correlation functions
Appears in: G-035

Ideal photodetection measures normally ordered field correlations; for one mode the normalized r-th order coherence is the r-th factorial moment over the mean to the r-th power.

$$ G^{(r)}=\langle \hat E^{(-)}(x_1)\cdots \hat E^{(-)}(x_r)\hat E^{(+)}(x_r)\cdots \hat E^{(+)}(x_1)\rangle,\qquad g^{(r)}=\langle(\hat a^\dagger)^r \hat a^r\rangle/\langle \hat n\rangle^r. $$

Conditions: $\langle\hat n\rangle\ne0$.
Audit: Correct (definition plus the ideal-absorber argument). Checked by: prior audit only for the detection model; the single-mode g^{(2)} case is checked exactly in the next entry.
QF code (exact): none. a definition of a measured quantity; QF's G1Correlation and G2Coherence act on finite Fock states only.

<a id="second-order-coherence-and-mandel-q"></a>
### Second-order coherence and Mandel Q
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

<a id="linear-quadrature-exponential-normal-order-conversion-factor"></a>
### Linear quadrature exponential: normal-order conversion factor
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

<a id="centered-quadrature-moments-from-normal-moments"></a>
### Centered quadrature moments from normal moments
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

<a id="weyl-and-antinormal-symbols-characteristic-functions"></a>
## Weyl and antinormal symbols, characteristic functions

<a id="differential-correspondences-for-multiplying-a-displacement"></a>
### Differential correspondences for multiplying a displacement
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

<a id="q-and-p-symbols-of-a-displacement"></a>
### Q and P symbols of a displacement
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

<a id="distributional-trace-of-displacement-and-operator-reconstruction"></a>
### Distributional trace of displacement and operator reconstruction
Appears in: P-017, P-020 (the $\hat\rho=\int\chi_W\hat D(-\gamma)$ line)

A displacement has a delta-function trace, displacements are orthogonal in the Hilbert-Schmidt sense, and any operator is rebuilt from its displacement transform.

$$ \operatorname{Tr}\hat D(\gamma)=\pi\delta^{(2)}(\gamma),\quad \operatorname{Tr}[\hat D(\alpha)\hat D(\beta)^\dagger]=\pi\delta^{(2)}(\alpha-\beta),\quad \hat A=\int\frac{d^2\gamma}{\pi}\operatorname{Tr}[\hat D(\gamma)\hat A]\hat D(-\gamma) $$

Conditions: Hilbert-Schmidt / distributional sense; $d^2\gamma=d\,\mathrm{Re}\gamma\,d\,\mathrm{Im}\gamma$.
Audit: Correct. Checked by: hand derivation: $\operatorname{Tr}\hat D(\gamma)=\sum_ne^{-|\gamma|^2/2}L_n(|\gamma|^2)$, and $\sum_nL_n(x)$ is $\pi\delta^{(2)}$ in this measure (Glauber's completeness of displacements).
QF code (exact): none. A distributional trace has no exact QF counterpart.

<a id="weyl-operator-from-displaced-parity-weyl-symbol-and-reconstruction"></a>
### Weyl operator from displaced parity, Weyl symbol and reconstruction
Appears in: P-019; G-034

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

**Also stated in G-034 (Wigner function as displaced parity):**

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

<a id="gaussian-smoothing-kernels-between-p-w-and-q-derived-normalization-correction"></a>
### Gaussian smoothing kernels between P, W and Q (derived normalization correction)
Appears in: P-019; P-015; C-017

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

**Also stated in P-015 (Q symbol is the Gaussian smoothing of the P symbol):**

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

**Also stated in C-017 (Wick, anti-Wick and Weyl symbols related by Gaussian smoothing):**

The Weyl symbol is the anti-Wick symbol smoothed by a Gaussian of width set by hbar, and the Wick symbol is the Weyl symbol smoothed once more by the same Gaussian.

$$ A_w=G_\hbar*A^c=e^{\hbar\Delta/4}A^c,\qquad A_c=G_\hbar*A_w=e^{\hbar\Delta/4}A_w,\qquad G_\hbar(X)=(\pi\hbar)^{-n}e^{-|X|^2/\hbar} $$

Audit: Correct. The Weyl symbol of |phi_z><phi_z| is 2^n e^{-|X-z|^2/hbar}; with the (2 pi hbar)^{-n} measure this gives the kernel G; per coordinate variance hbar/2 gives e^{hbar Delta/4}. Checked by: hand derivation.
QF code (exact): none. QF's `WignerFunction`, `HusimiQFunction` and `SOrderedFunction` act on finite-size states, not on operator symbols.

<a id="ordered-characteristic-functions-differ-by-gaussian-factors"></a>
### Ordered characteristic functions differ by Gaussian factors
Appears in: P-020, P2-035; G-028, G-030 (s-ordered form $C(\lambda,s)=e^{s|\lambda|^2/2}C_W$)

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

**Also stated in G-028 (Three ordered characteristic functions and their Gaussian conversion):**

The Weyl, normal and antinormal characteristic functions are traces of the displacement in three orders; they differ only by the Gaussian factors e^{-+|lambda|^2/2}, which is the s-ordered family.

$$ C_W(\lambda)=\operatorname{Tr}(\hat\rho \hat D(\lambda)),\quad C_N(\lambda)=\operatorname{Tr}(\hat\rho e^{\lambda \hat a^\dagger}e^{-\lambda^*\hat a}),\quad C_A(\lambda)=\operatorname{Tr}(\hat\rho e^{-\lambda^*\hat a}e^{\lambda \hat a^\dagger}),\quad C_W=e^{-|\lambda|^2/2}C_N=e^{|\lambda|^2/2}C_A. $$

Audit: Correct, including the atlas's correction of the printed Eq. (3.128c). Checked by: QF exact (the two ordered displacements hold for any rho under the trace); hand derivation that a creation-first product is the normal function, which differs from the antinormal one by e^{|lambda|^2}.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{DisplacementOperator[l, Infinity, "Ordering" -> "Normal"], DisplacementOperator[l, Infinity, "Ordering" -> "Antinormal"]}
```
Output: `{E^(l*SuperDagger[\[FormalA]1])**E^(-(\[FormalA]1*Conjugate[l]))/E^((l*Conjugate[l])/2), E^((l*Conjugate[l])/2)*E^(-(\[FormalA]1*Conjugate[l]))**E^(l*SuperDagger[\[FormalA]1])}`. Matches: D = e^{-|lambda|^2/2} e^{lambda a^dagger} e^{-lambda^* a} = e^{|lambda|^2/2} e^{-lambda^* a} e^{lambda a^dagger}, giving C_W = e^{-|lambda|^2/2} C_N = e^{|lambda|^2/2} C_A.

<a id="normally-ordered-moments-from-the-normal-characteristic-function"></a>
### Normally ordered moments from the normal characteristic function
Appears in: P-020; G-029

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

**Also stated in G-029 (Ordered moments as derivatives of characteristic functions):**

Ordered moments are read off as mixed derivatives of the matching ordered characteristic function at the origin.

$$ \langle(\hat a^\dagger)^m \hat a^n\rangle=\left.\partial_\lambda^m\partial_{(-\lambda^*)}^nC_N\right|_0,\quad\langle \hat a^m(\hat a^\dagger)^n\rangle=\left.\partial_\lambda^n\partial_{(-\lambda^*)}^mC_A\right|_0. $$

Audit: Correct. Checked by: hand derivation (Taylor expansion of the ordered exponentials under the trace).
QF code (exact): none. a Taylor-coefficient identity under a trace; no QF call produces it, and the coherent-basis BosonicMatrixElement leaves exponentials of a or a^dagger unevaluated.

<a id="exact-attenuation-solution-for-ordered-characteristic-functions"></a>
### Exact attenuation solution for ordered characteristic functions
Appears in: P2-036, P2-039 (thermal $\mu(t)=p\nu_0+q\nu$)

Under thermal attenuation each ordered characteristic function is the initial one at contracted argument times a Gaussian that grows to the noise level.

$$ \chi_k(\eta,t)=e^{-q(\nu+\sigma_k)|\eta|^2}\chi_k(\sqrt p\,\eta,0),\qquad \partial_t\chi_k=-\gamma[\tfrac12\eta_j\partial_{\eta_j}+(\nu+\sigma_k)|\eta|^2]\chi_k,\quad p=e^{-\gamma t},\ q=1-p $$

Audit: Correct. Checked by: exact WL computation: with $\eta,\bar\eta$ independent and an arbitrary initial function, the stated solution satisfies the PDE identically. For thermal input $\chi_N=e^{-\nu_0|\eta|^2}$ the solution is thermal with $\mu=p\nu_0+q\nu$ (hand).
QF code (exact): none. QF has no characteristic-function evolution.

<a id="fokker-planck-drift-diffusion-and-covariance-evolution"></a>
### Fokker-Planck drift, diffusion and covariance evolution
Appears in: P2-037

The ordered quasiprobabilities obey Fokker-Planck equations with drift $-\gamma\alpha/2$, and their covariance relaxes linearly to the noise level plus the ordering offset.

$$ W_N=P,\ W_O=W_{\rm Wigner},\ W_A=Q/\pi;\qquad D(t)=pD(0)+q(\nu+\sigma_k)I_2/2 $$

Audit: Correct. Checked by: hand derivation; a Gaussian factor $e^{-s|\eta|^2}$ in $\chi_k$ corresponds to variance $s/2$ per real quadrature, giving the additive $q(\nu+\sigma_k)/2$.
QF code (exact): none. QF has no quasiprobability evolution.

<a id="gaussian-convolution-kernel-for-attenuation-sign-corrected"></a>
### Gaussian convolution kernel for attenuation (sign corrected)
Appears in: P2-038

Attenuation acts on each ordered quasiprobability as convolution with a normalized Gaussian centered at $\sqrt p\,\alpha'$.

$$ W_k(\alpha,t)=\int d^2\alpha'\,W_k(\alpha',0)\frac{\exp[-|\alpha-\sqrt p\,\alpha'|^2/(q(\nu+\sigma_k))]}{\pi q(\nu+\sigma_k)} $$

Audit: Correct as corrected; the printed positive exponent is not integrable. Checked by: exact WL computation: the corrected kernel integrates to 1, and its Fourier transform is $e^{-s|k|^2/4}$ with $s=q(\nu+\sigma_k)$, the Gaussian factor of the characteristic-function solution.
QF code (exact): none. QF has no quasiprobability evolution.

<a id="weyl-quantization-kernel-convention"></a>
### Weyl quantization kernel convention
Appears in: C-013

An operator acts through its Weyl symbol evaluated at the midpoint, and the symbol is recovered from the integral kernel by a Fourier transform in the difference variable.

$$ [\hat A\psi](x)=(2\pi\hbar)^{-n}\int A((x+y)/2,\xi)e^{i(x-y)\cdot\xi/\hbar}\psi(y)\,dy\,d\xi,\qquad A(x,\xi)=\int e^{-iu\cdot\xi/\hbar}K(x+u/2,x-u/2)\,du $$

Audit: Correct (the second is the Fourier inverse of the first in u = x - y). Checked by: hand derivation.
QF code (exact): none. QF has no Weyl-symbol map for operators.

<a id="symbol-traces-and-wigner-overlaps"></a>
### Symbol traces and Wigner overlaps
Appears in: C-014

With the kernel convention of C-013, Hilbert-Schmidt norms, traces of products and expectation values become phase-space integrals of symbols with the factor (2 pi hbar)^{-n}.

$$ \|\hat A\|_{\rm HS}^2=(2\pi\hbar)^{-n}\int|A|^2,\quad\operatorname{Tr}(\hat A\hat B)=(2\pi\hbar)^{-n}\int AB,\quad \langle\psi,\hat A\phi\rangle=(2\pi\hbar)^{-n}\int A\,W_{\psi,\phi},\quad |\langle\psi,\phi\rangle|^2=(2\pi\hbar)^{-n}\int W_\psi W_\phi $$

Conditions: W is the Weyl symbol of the rank-one operator (unnormalized, integral (2 pi hbar)^n).
Audit: Correct for the C-013 kernel convention. Checked by: hand derivation.
QF code (exact): none. QF has no Weyl-symbol calculus; its `WignerFunction` takes a finite-size state and returns a normalized density.

<a id="moyal-product-and-exact-quadratic-commutator"></a>
### Moyal product and exact quadratic commutator
Appears in: C-016

The Weyl symbol of a product is the Moyal star product, and when one factor is at most quadratic the commutator is exactly i hbar times the quantized Poisson bracket.

$$ A\star B=A\exp\!\left[\frac{i\hbar}{2}\left(\overleftarrow\partial_q\cdot\overrightarrow\partial_p-\overleftarrow\partial_p\cdot\overrightarrow\partial_q\right)\right]B,\qquad [\hat H,\hat B]=i\hbar\,\widehat{\{H,B\}}\ (\deg H\le2) $$

Audit: Correct. q*p - p*q = i hbar reproduces [Q,P] = i hbar, and H = q, B = p gives i hbar {q,p} = i hbar. Checked by: hand derivation.
QF code (exact): none. QF has no Moyal product.

<a id="weyl-symbol-of-a-metaplectic-operator-cayley-form"></a>
### Weyl symbol of a metaplectic operator (Cayley form)
Appears in: C-030; oscillator case C-053

Away from det(I + F) = 0, the Weyl symbol of a metaplectic operator is a Gaussian in the Cayley transform of F, with prefactor 2^n |det(I + F)|^{-1/2} and a path-dependent phase.

$$ R_w(F;X)=\frac{2^n e^{i\pi\nu}}{|\det(I_{2n}+F)|^{1/2}}\exp\!\left[-\frac{i}{\hbar}X^TJ(I_{2n}-F)(I_{2n}+F)^{-1}X\right] $$

Conditions: det(I+F) != 0; nu and phase from the path.
Audit: Correct (with the atlas's added 2^n). F = I gives 1. Oscillator F = e^{tJ}: (I-F)(I+F)^{-1} = -tan(t/2) J, so the exponent is -i tan(t/2)|X|^2/hbar and the prefactor 1/cos(t/2), the known symbol of e^{-itH_os}. Checked by: hand derivation.
QF code (exact): none. QF has no Weyl-symbol calculus.

<a id="translation-coefficient-of-a-metaplectic-operator"></a>
### Translation coefficient of a metaplectic operator
Appears in: C-031

The translation coefficient Tr[T(-z) R(F)] of a metaplectic operator is the Fresnel transform of its Weyl symbol: a Gaussian in the inverse Cayley transform with a signature phase.

$$ \widetilde R(F;z)=\frac{c_w\,e^{-i\pi\operatorname{sgn}(Q_F)/4}}{2^n|\det Q_F|^{1/2}}\exp\!\left[-\frac{i}{4\hbar}z^TJ(I+F)(I-F)^{-1}z\right],\quad Q_F=J(I-F)(I+F)^{-1};\qquad \widetilde R=\frac{1}{2i\sin(t/2)}e^{\frac{i}{4\hbar}\cot(t/2)|z|^2}\ \text{for}\ \hat R=e^{-it(\hat N+1/2)} $$

Conditions: I + F and I - F invertible; 0 < t < pi for the example.
Audit: Correct. Direct trace: Tr[D(-alpha) e^{-it(N+1/2)}] = e^{-it/2} e^{-x/2} sum L_n(x) w^n = 1/(2i sin(t/2)) exp(i x cot(t/2)/2), w = e^{-it}, x = |z|^2/(2 hbar), via the Laguerre generating function; the general chart formula with c_w = 1/cos(t/2), Q_F = tan(t/2) I, sgn = 2 gives the same. Checked by: hand derivation.
QF code (exact): none. QF has no trace of an operator on the full Fock space; its exact Laguerre matrix elements of `DisplacementOperator[al]` are one input of the hand derivation, not a trace.

<a id="wigner-function-of-a-one-mode-squeezed-gaussian"></a>
### Wigner function of a one-mode squeezed Gaussian
Appears in: C-049

The Wigner function of the one-mode squeezed Gaussian e^{-delta x^2/(2 hbar)} is a tilted Gaussian ellipse fixed by the real and imaginary parts of delta.

$$ W(q,p)=2\exp[-q^2\Re\delta/\hbar-(p+q\Im\delta)^2/(\hbar\Re\delta)],\qquad \psi\propto e^{-\delta x^2/(2\hbar)},\ \delta=\frac{1-\omega}{1+\omega} $$

Audit: Correct (unnormalized Weyl symbol of the projector; integral 2 pi hbar). Forcing line: psi(q+u/2) conj(psi(q-u/2)) = exp(-Re(delta) q^2 - Re(delta) u^2/4 - i Im(delta) q u), whose u-Fourier transform gives the (p + q Im delta)^2 term. Checked by: hand derivation.
QF code (exact): none. QF's `WignerFunction` takes a finite-size state; QF has no exact Weyl symbol of an operator.

<a id="weyl-symbol-of-the-oscillator-exponential"></a>
### Weyl symbol of the oscillator exponential
Appears in: C-053 (special case of C-030)

The Weyl symbol of the oscillator propagator is (cos(t/2))^{-n} e^{i tan(t/2)|X|^2}, and in imaginary time (cosh(tau/2))^{-n} e^{-tanh(tau/2)|X|^2}.

$$ [e^{it\hat H_{\rm os}}]_w=(\cos(t/2))^{-n}\exp[i\tan(t/2)|X|^2],\qquad [e^{-\tau\hat H_{\rm os}}]_w=(\cosh(\tau/2))^{-n}\exp[-\tanh(\tau/2)|X|^2] $$

Conditions: hbar = 1, |t| < pi.
Audit: Correct (power -n). Check: (2 pi)^{-1} times the integral of the tau symbol is 1/(2 sinh(tau/2)) = Tr e^{-tau H_os}. Checked by: hand derivation.
QF code (exact): none. QF has no Weyl-symbol calculus.

<a id="fermionic-weyl-translation-expansion"></a>
### Fermionic Weyl translation expansion
Appears in: C-074

An operator expands in fermionic translations with coefficient Str(A T(xi)).

$$ \hat A=\int A_w(\xi)\hat T(-\xi)d^2\xi,\qquad A_w(\xi)=\operatorname{Str}(\hat A\hat T(\xi)) $$

Audit: Unchecked here: depends on the source's d^2 xi normalization; prior audit only (verified finite).
QF code (exact): none. no fermionic layer.

<a id="normal-antinormal-and-s-ordered-fermionic-translations"></a>
### Normal, antinormal and s-ordered fermionic translations
Appears in: C-075

Normal and antinormal fermionic translations differ from the symmetric one by the scalar Grassmann factor e^{+gamma^* gamma/2} or e^{-gamma^* gamma/2}.

$$ \hat T_N(\gamma)=e^{\hat c^\dagger\gamma}e^{-\gamma^*\hat c}=\hat T(\gamma)e^{\gamma^*\gamma/2},\quad \hat T_A(\gamma)=e^{-\gamma^*\hat c}e^{\hat c^\dagger\gamma}=\hat T(\gamma)e^{-\gamma^*\gamma/2},\quad \hat T_s=\hat T e^{s\gamma^*\gamma/2} $$

Audit: Correct. Checked by: exact fermion matrices (extended Jordan-Wigner realization).
QF code (exact): none. no fermionic layer.

<a id="ordinary-trace-from-mixed-fermionic-symbols"></a>
### Ordinary trace from mixed fermionic symbols
Appears in: C-076

The ordinary trace of a product of fermionic operators is a paired Berezin integral of one contravariant and one covariant symbol at a reflected, doubled argument, and the trace of one operator is 2^n times its symbol at zero.

$$ \operatorname{Tr}(\hat H\hat G)=2^{-n}\widehat I_\zeta\!\left[H^w(\zeta)G_w(-2\zeta)\right],\qquad \operatorname{Tr}\hat A=2^nA^w(0) $$

Audit: Correct as far as checked: Tr A = 2^n A^w(0) holds by hand for I and c^dag c (symbol alpha^* alpha + 1/2 gives 1). The reflected argument -2 zeta: prior audit only.
QF code (exact): none. no fermionic layer.

<a id="fermionic-star-product"></a>
### Fermionic star product
Appears in: C-077

The fermionic star product contracts odd derivatives symmetrically, and the propagator's symbol solves the star Schrodinger equation.

$$ G\star H=G\exp\!\left[\tfrac12\left(\overleftarrow\partial_{\alpha^*}\cdot\overrightarrow\partial_\alpha+\overleftarrow\partial_\alpha\cdot\overrightarrow\partial_{\alpha^*}\right)\right]H,\qquad i\partial_tU^w=H\star U^w $$

Audit: Correct (one-mode consistency: alpha * alpha^* = alpha alpha^* + 1/2 and alpha^* * alpha = alpha^* alpha + 1/2 reproduce the symbols of c c^dag and c^dag c, which sum to 1). Checked by: hand derivation.
QF code (exact): none. no fermionic layer.

<a id="fermi-oscillator-exponential-and-its-weyl-symbol"></a>
### Fermi oscillator exponential and its Weyl symbol
Appears in: C-078

The exponential of one fermion number operator is linear in it, and its Weyl symbol is (e^z+1)/2 + (e^z-1) alpha^* alpha.

$$ e^{z\hat c^\dagger\hat c}=\hat I+(e^z-1)\hat c^\dagger\hat c,\quad (\hat c^\dagger\hat c)^w=\alpha^*\alpha+\tfrac12,\quad [e^{z\hat c^\dagger\hat c}]^w=\frac{e^z+1}2\exp[2\tanh(z/2)\alpha^*\alpha]=\frac{e^z+1}2+(e^z-1)\alpha^*\alpha $$

Audit: Correct ((e^z + 1) tanh(z/2) = e^z - 1). Checked by: hand derivation.
QF code (exact): none. no fermionic layer.

<a id="fermionic-cayley-weyl-symbol-corrected-sign"></a>
### Fermionic Cayley Weyl symbol (corrected sign)
Appears in: C-085

The Weyl symbol of a quadratic fermion propagator is a Grassmann Gaussian in the Cayley transform of F_t; the atlas corrects the sign of the exponent.

$$ U_t^w(\xi)=\det\!\left(\frac{F_t+I}{2}\right)^{1/2}\exp\!\left[-\binom{\xi^*}{\xi}\cdot(F_t-I)(F_t+I)^{-1}\binom\xi{\xi^*}\right]\ \text{(atlas correction: minus sign)} $$

Audit: Correct as corrected by the atlas. One mode, L = 0: F = diag(e^{2i theta}, e^{-2i theta}), theta = tm/4; the printed plus sign gives +2i tan(theta) xi^* xi, while C-078 gives cos(theta) exp(-2i tan(theta) xi^* xi). Checked by: hand derivation.
QF code (exact): none. no fermionic layer.

<a id="super-star-product"></a>
### Super star product
Appears in: C-089

The super star product adds the antisymmetric bosonic contraction to the symmetric fermionic one.

$$ G\star H=G\exp\tfrac12\left[\overleftarrow\partial_\zeta\overrightarrow\partial_{\bar\zeta}-\overleftarrow\partial_{\bar\zeta}\overrightarrow\partial_\zeta+\overleftarrow\partial_{\alpha^*}\overrightarrow\partial_\alpha+\overleftarrow\partial_\alpha\overrightarrow\partial_{\alpha^*}\right]H $$

Audit: Convention: correct when zeta is the symbol of a (then zeta * conj(zeta) - conj(zeta) * zeta = 1 = [a, a^dag]); with zeta as the symbol of a^dag the bosonic sign flips. Checked by: hand derivation.
QF code (exact): none. QF has no symbol calculus and no fermionic layer.

<a id="coherent-projector-differential-correspondences"></a>
### Coherent-projector differential correspondences
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

<a id="p-representation-optical-equivalence-and-q-function"></a>
### P representation, optical equivalence, and Q function
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

<a id="glauber-p-function-of-a-number-state"></a>
### Glauber P function of a number state
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

<a id="complex-fourier-transform-between-p-and-the-normal-characteristic-function"></a>
### Complex Fourier transform between P and the normal characteristic function
Appears in: G-031

P and the normal characteristic function form a complex Fourier pair with kernel e^{lambda alpha^* - lambda^* alpha} and normalization pi^{-2}.

$$ C_N(\lambda)=\int P(\alpha)e^{\lambda\alpha^*-\lambda^*\alpha}d^2\alpha,\quad P(\alpha)=\frac1{\pi^2}\int C_N(\lambda)e^{\lambda^*\alpha-\lambda\alpha^*}d^2\lambda. $$

Audit: Correct (sign and pi^{-2}). Checked by: exact WL computation (script wl/031.wl): the inverse transform of e^{-nbar|lambda|^2} gives e^{-|alpha|^2/nbar}/(pi nbar) and the forward transform returns e^{-nbar|lambda|^2}.
QF code (exact): none. QF has no P function or characteristic-function transform.

<a id="q-and-wigner-functions-of-coherent-and-number-states"></a>
### Q and Wigner functions of coherent and number states
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

<a id="coherent-overlap-of-squeezed-vacuum-and-q-function-of-a-displaced-squeezed-state"></a>
### Coherent overlap of squeezed vacuum and Q function of a displaced squeezed state
Appears in: G-063

The coherent overlap of squeezed vacuum is a Gaussian in beta with the e^{i theta} tanh r quadratic term, and the Q function of the displaced state is that overlap squared, centered at alpha.

$$ \langle\beta|\hat S(\xi)|0\rangle=\frac1{\sqrt{\cosh r}}e^{-|\beta|^2/2-e^{i\theta}\beta^{*2}\tanh r/2},\quad Q(\beta)=\frac1{\pi\cosh r}\exp\{-|\delta|^2-\tanh r\,\operatorname{Re}(e^{-i\theta}\delta^2)\},\ \delta=\beta-\alpha. $$

Audit: Correct (the atlas's correction of printed Eq. 7.82). Checked by: exact WL computation (script wl/063.wl): the Fock sum of <beta|2m><2m|S|0> with the G-056 amplitudes gives the overlap, and |<beta - alpha|S|0>|^2/pi equals the stated Q; hand derivation of D^dagger(alpha)|beta> = phase |beta - alpha>.
QF code (exact): none. BosonicMatrixElement[{beta, 0}, Exp[c ad ** ad], "Basis" -> "Coherent"] returns unevaluated (exponentials in the coherent basis are a known missing feature).

<a id="wigner-function-of-a-displaced-squeezed-state"></a>
### Wigner function of a displaced squeezed state
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

<a id="cat-wigner-function"></a>
### Cat Wigner function
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

<a id="grassmann-calculus-and-fermionic-coherent-states"></a>
## Grassmann calculus and fermionic coherent states

<a id="fermionic-symbols-which-side-the-derivative-acts-on"></a>
### Fermionic symbols: which side the derivative acts on
Appears in: B-016

For fermions, multiplying by a on the right or by a* on the left keeps normal order, so the normal symbol is just multiplied on the same side by the Grassmann variable.

$$ N_{\hat A\hat a}=N_{\hat A}a,\qquad N_{\hat a^*\hat A}=a^*N_{\hat A} $$

Audit: Correct: right multiplication by a and left multiplication by a* keep normal order, so the symbol is multiplied on the same side; the derivative sides only enter for the other two products. Checked by: hand derivation: a* (a*^p a^q) and (a*^p a^q) a are already normal ordered.
QF code (exact): none. No fermionic layer.

<a id="grassmann-algebra-relations"></a>
### Grassmann algebra relations
Appears in: B-029

Grassmann generators anticommute and square to zero, so n of them span an algebra of dimension 2^n.

$$ x_i x_j=-x_jx_i,\qquad x_i^2=0,\qquad \dim=2^n $$

Audit: Correct. Checked by: hand derivation: a basis is the 2^n ordered monomials over subsets of the generators.
QF code (exact): none. No fermionic layer.

<a id="left-versus-right-derivatives"></a>
### Left versus right derivatives
Appears in: B-030

On a monomial of degree p, the right and left Grassmann derivatives differ by the sign (-1)^{p-1}.

$$ \partial^{R}_{x}m=(-1)^{p-1}\,\partial^{L}_{x}m\quad\text{on a monomial }m\text{ of degree }p $$

Audit: Correct. Checked by: hand derivation: moving the target generator to the right end instead of the left end takes p - 1 transpositions more in total.
QF code (exact): none. No fermionic layer.

<a id="berezin-integration-rules"></a>
### Berezin integration rules
Appears in: B-031

Berezin integration kills constants and gives 1 on the integration variable itself.

$$ \int dx_i=0,\qquad \int x_i\,dx_i=1 $$

Audit: Correct as the defining rules. The partial-integration warning is a source statement, Unchecked. Checked by: hand derivation (definition).
QF code (exact): none. No fermionic layer.

<a id="grassmann-change-of-variables"></a>
### Grassmann change of variables
Appears in: B-032

Under a linear change of Grassmann variables the integral picks up the inverse determinant, the opposite of commuting variables.

$$ \int f(x)\,dx_n\cdots dx_1=(\det A)^{-1}\int f(Ay)\,dy_n\cdots dy_1 $$

Audit: Correct. Checked by: hand derivation: for f = x_1...x_n, f(Ay) = det A y_1...y_n, so the right side is (det A)^{-1} det A = 1, the left side.
QF code (exact): none. No fermionic layer.

<a id="clifford-operators-from-multiplication-and-derivative"></a>
### Clifford operators from multiplication and derivative
Appears in: B-034

Grassmann multiplication and left differentiation satisfy the CAR, so their sum and difference generate a Clifford algebra.

$$ \{\hat Q_i,\hat Q_j\}=\{\hat P_i,\hat P_j\}=2\delta_{ij}\hat I,\quad \{\hat P_i,\hat Q_j\}=0,\qquad \hat Q_j=\hat X_j+\hat D_j,\ \hat P_j=(\hat X_j-\hat D_j)/i $$

Audit: Correct. Checked by: exact fermion matrices (2 modes; X = c^dag, D = c; all residuals 0).
QF code (exact): none. No fermionic layer.

<a id="grassmann-derivatives-leibniz-rule-and-left-right-relation"></a>
### Grassmann derivatives: Leibniz rule and left/right relation
Appears in: C-066

The Grassmann derivative obeys a graded Leibniz rule, and the right and left derivatives differ by a parity sign.

$$ \partial_\theta(fg)=(\partial_\theta f)g+(-1)^{\pi(f)}f(\partial_\theta g),\qquad f\overleftarrow\partial_\theta=(-1)^{1+\pi(f)}\overrightarrow\partial_\theta f $$

Audit: Correct. Checked by: exact WL computation (exterior algebra, odd and even test elements).
QF code (exact): none. no fermionic layer.

<a id="berezin-integration-determinant-rule-and-paired-normalization"></a>
### Berezin integration: determinant rule and paired normalization
Appears in: C-067; the same paired functional is used in C-068, C-069, C-098, C-099

Berezin integration rescales by det A under a linear change of generators, and the paired convention normalizes e^{-theta^* theta} to 1.

$$ \int f(A\theta)\,d\theta=\det A\int f(\theta)\,d\theta,\qquad \widehat I_\theta\Big[\prod_j\theta_j\theta_j^*\Big]=1,\ \widehat I_\theta[e^{-\theta^*\theta}]=1,\ \widehat I_\theta=\partial_{\theta_1^*}\partial_{\theta_1}\cdots\partial_{\theta_n^*}\partial_{\theta_n} $$

Audit: Correct as a declared convention (the atlas's repair). Checked by: exact WL computation (exterior algebra: 3-by-3 determinant rule; paired Gaussian gives 1 for n = 2).
QF code (exact): none. no fermionic layer.

<a id="normal-ordered-monomial-expansion-of-fermionic-operators"></a>
### Normal-ordered monomial expansion of fermionic operators
Appears in: C-070

Every operator on n fermionic modes is a unique combination of normal-ordered monomials.

$$ \hat A=\sum_{\epsilon,\epsilon'\in\{0,1\}^n}A_{\epsilon\epsilon'}(\hat c^\dagger)^\epsilon\hat c^{\epsilon'} $$

Audit: Correct (4^n monomials span the 2^n-by-2^n matrices). Checked by: hand derivation (dimension count).
QF code (exact): none. no fermionic layer.

<a id="grassmann-weyl-translation"></a>
### Grassmann Weyl translation
Appears in: C-071

The fermionic translation with Grassmann labels has an even generator, so the central BCH formula factorizes it and composes it up to a scalar Grassmann phase, and it shifts c by gamma.

$$ \hat T(\gamma)=e^{\hat c^\dagger\cdot\gamma-\gamma^*\cdot\hat c}=e^{-\gamma^*\cdot\gamma/2}e^{\hat c^\dagger\cdot\gamma}e^{-\gamma^*\cdot\hat c},\quad \hat T(\alpha)\hat T(\gamma)=\hat T(\alpha+\gamma)e^{-(\alpha^*\cdot\gamma+\alpha\cdot\gamma^*)/2},\quad \hat T^{-1}\hat c\hat T=\hat c+\gamma,\ \hat T^{-1}\hat c^\dagger\hat T=\hat c^\dagger+\gamma^* $$

Audit: Correct. Checked by: exact fermion matrices (one mode plus four extra Jordan-Wigner modes standing for the odd labels; all identities exact; T(gamma)^{-1} = T(-gamma)). The adjoint statement T^* = T(-gamma) needs the Grassmann conjugation and was not tested.
QF code (exact): none. no fermionic layer.

<a id="fermionic-coherent-state-overlap-resolution"></a>
### Fermionic coherent state, overlap, resolution
Appears in: C-072

The fermionic coherent state is an exponential of theta.gamma, an eigenvector of each c_j with Grassmann eigenvalue gamma_j; overlaps are Gaussian in the labels and the states resolve the identity.

$$ \psi_\gamma=e^{-\gamma^*\gamma/2}e^{\theta\cdot\gamma},\ \hat c_j\psi_\gamma=\gamma_j\psi_\gamma,\quad \langle\psi_\gamma,\psi_\alpha\rangle=\exp[\gamma^*\alpha-(\gamma^*\gamma+\alpha^*\alpha)/2],\quad \psi=\int\langle\psi_\gamma,\psi\rangle\psi_\gamma\,d\gamma^*d\gamma $$

Audit: Correct for the eigenvalue. Checked by: exact fermion matrices (extended Jordan-Wigner realization). Overlap and resolution: standard results, prior audit only.
QF code (exact): none. no fermionic layer.

<a id="supertrace-and-berezinian"></a>
### Supertrace and Berezinian
Appears in: C-073

The supertrace weights the trace by the parity operator and is graded-cyclic; the Berezinian is the super analogue of the determinant and satisfies Ber e^H = e^{Str H}.

$$ \operatorname{Str}\hat A=\operatorname{Tr}(\hat P\hat A),\ \hat P=(-1)^{\hat N},\ \operatorname{Str}(\hat A\hat B)=(-1)^{\pi(A)\pi(B)}\operatorname{Str}(\hat B\hat A),\quad \operatorname{Ber}\begin{pmatrix}A&B\\C&D\end{pmatrix}=\frac{\det(A-BD^{-1}C)}{\det D},\ \operatorname{Ber}(e^H)=e^{\operatorname{Str}H} $$

Audit: Correct (standard). Checked by: hand derivation (odd A, B: Tr(P A B) = -Tr(P B A) because P anticommutes with odd matrices).
QF code (exact): none. no fermionic layer.

<a id="mixed-boson-fermion-displacement"></a>
### Mixed boson-fermion displacement
Appears in: C-088

A mixed boson-fermion displacement factors into its bosonic and fermionic normal-ordered parts, since the even generators of the two parts commute.

$$ \hat T(z,\gamma)=e^{-(|z|^2+\gamma^*\gamma)/2}e^{z\hat a_B^\dagger}e^{\hat c_F^\dagger\gamma}e^{-\bar z\hat a_B}e^{-\gamma^*\hat c_F} $$

Audit: Correct (each part factorizes by C-006 and C-071). Checked by: hand derivation.
QF code (exact): none. the bosonic factor is QF's `DisplacementOperator[al, Infinity, "Ordering" -> "Normal"]` (C-006), but there is no fermionic layer for the other factor.

<a id="supercommutator"></a>
### Supercommutator
Appears in: C-090

The supercommutator is the commutator with a sign that turns it into an anticommutator when both operators are odd.

$$ [\hat A,\hat B]_s=\hat A\hat B-(-1)^{\pi(\hat A)\pi(\hat B)}\hat B\hat A $$

Audit: Correct (definition). Checked by: hand derivation.
QF code (exact): none. QF has no graded (super) algebra; no fermionic layer.

<a id="grassmann-fourier-involution"></a>
### Grassmann Fourier involution
Appears in: C-098

The Grassmann Fourier transform is an involution in the paired convention; it maps 1 to the ordered top product of the alpha_j alpha_j^*.

$$ (\widehat{\mathcal F}f)(\alpha,\alpha^*)=\widehat I_\xi[e^{\xi^*\cdot\alpha-\alpha^*\cdot\xi}f],\quad \widehat{\mathcal F}^2=\hat I,\quad \widehat{\mathcal F}1=\alpha\alpha^*,\ \widehat{\mathcal F}\xi=\alpha,\ \widehat{\mathcal F}\xi^*=-\alpha^*,\ \widehat{\mathcal F}(\xi\xi^*)=1,\quad \widehat{\mathcal F}1=\prod_j\alpha_j\alpha_j^* $$

Audit: Correct (the atlas's corrected F1). Checked by: exact WL computation (exterior algebra: n = 1 table exact; inversion formula and F^2 = identity on all monomials for n = 1, 2; F1 for n = 2 is the ordered top product).
QF code (exact): none. no fermionic layer.

<a id="holomorphic-grassmann-reproducing-kernel"></a>
### Holomorphic Grassmann reproducing kernel
Appears in: C-099

The exponential Grassmann kernel e^{-(theta - gamma).gamma^*} reproduces every holomorphic Grassmann function under paired integration.

$$ \psi(\theta)=\widehat I_\gamma[\psi(\gamma)E(\theta,\gamma,\gamma^*)],\qquad E=e^{-(\theta-\gamma)\cdot\gamma^*}=\prod_j[1-(\theta_j-\gamma_j)\gamma_j^*] $$

Audit: Correct. Checked by: exact WL computation (exterior algebra, all four holomorphic monomials, n = 2).
QF code (exact): none. no fermionic layer.

<a id="gaussian-integrals-and-pfaffians"></a>
## Gaussian integrals and Pfaffians

<a id="real-gaussian-integral"></a>
### Real Gaussian integral
Appears in: B-022

The normalized real Gaussian integral of e^{-(Af, f)} is det A^{-1/2}.

$$ \int e^{-(\hat Af,f)}\Pi df=(\det \hat A)^{-1/2},\qquad \hat A>0 $$

Conditions: normalized measure pi^{-n/2} d^n f; infinite dimension needs A - I trace class.
Audit: Correct. Checked by: exact WL computation (2x2 positive rational A, exact `Integrate`, difference 0).
QF code (exact): none. A c-number Gaussian integral; no second-quantization call is involved.

<a id="doubled-complex-gaussian-with-sources"></a>
### Doubled complex Gaussian with sources
Appears in: B-023

A complex Gaussian with both holomorphic and antiholomorphic quadratic terms and linear sources integrates to the inverse square root of a row-swapped block determinant times the exponential of half the inverse form on the sources.

$$ \int e^{-\frac12(f\ f^*)\widehat{\mathcal A}\binom f{f^*}+(\varphi_1\ \varphi_2)\binom f{f^*}}\Pi df\,df^* =\det\!\begin{pmatrix}\hat A_{21}&\hat A_{22}\\\hat A_{11}&\hat A_{12}\end{pmatrix}^{-1/2} e^{\frac12(\varphi_1\ \varphi_2)\widehat{\mathcal A}^{-1}\binom{\varphi_1}{\varphi_2}} $$

Conditions: A transpose-symmetric, A12 - I and A21 - I trace class, A11, A22 Hilbert-Schmidt, real part of the form bounded below.
Audit: Correct, with the row-swapped determinant as printed (for one mode it equals -det A). Checked by: exact WL computation (one mode, complex A11, A22, A12 = A21 = 1 and complex sources, measure d^2f/pi: the exact `Integrate` and the formula agree; `Simplify` stalls on the nested radicals, and the difference evaluates to 0 at 99 digits).
QF code (exact): none. A c-number Gaussian integral; no second-quantization call is involved.

<a id="fermionic-gaussian-is-a-pfaffian"></a>
### Fermionic Gaussian is a Pfaffian
Appears in: B-033

The Berezin integral of a Grassmann Gaussian with skew matrix M is the Pfaffian of M.

$$ \int e^{\frac12x^TMx}\,dx_n\cdots dx_1=\operatorname{Pf}(M),\quad M^T=-M $$

Audit: Correct. Checked by: exact WL computation (explicit Grassmann algebra, 4 generators, complex skew M: top coefficient minus Pf is 0).
QF code (exact): none. No fermionic layer.

<a id="doubled-fermionic-gaussian-corrected-determinant"></a>
### Doubled fermionic Gaussian (corrected determinant)
Appears in: B-036

The paired Grassmann Gaussian integral equals the square root of the unpermuted block determinant, fixed in sign as (-1)^{n(n-1)/2} times the Pfaffian.

$$ \int e^{\frac12(\alpha\ \alpha^*)\widehat{\mathcal A}\binom\alpha{\alpha^*}}\Pi d\alpha^*d\alpha =\det\!\begin{pmatrix}\hat A_{11}&\hat A_{12}\\\hat A_{21}&\hat A_{22}\end{pmatrix}^{1/2} = (-1)^{n(n-1)/2}\operatorname{Pf}(\mathcal A) $$

Conditions: normalized paired integration, int e^{-alpha* alpha} = 1; block order (alpha_1..alpha_n, alpha*_1..alpha*_n).
Audit: Correct as corrected in the atlas. Checked by: exact WL computation (explicit Grassmann algebra, n = 1, 2, 3, complex skew matrices: value matches (-1)^{n(n-1)/2} Pf and its square matches det; the plain Pf fails for n = 2, 3).
QF code (exact): none. No fermionic layer.

<a id="determinant-and-log-determinant-series"></a>
### Determinant and log-determinant series
Appears in: P-048

The determinant is the exponential of the trace of the logarithm, and $-\log\det(I-X)$ is the trace power series, convergent for spectral radius below 1.

$$ \det A=e^{\operatorname{Tr}\log A},\qquad -\log\det(I_N-X)=\sum_{m\ge1}\frac{\operatorname{Tr}X^m}{m},\ \rho(X)<1;\qquad \psi_Z(W)=B(Z,Z)^{-1/2}B(Z,W) $$

Audit: Correct (the kernel $\psi_Z$ is a definition). Checked by: exact WL computation: for a fully symbolic $3\times3$ matrix the series of $-\log\det(I-tX)$ agrees with $\sum_mt^m\operatorname{Tr}X^m/m$ through $t^5$.
QF code (exact): none. Not an operator identity.

<a id="fermionic-gaussian-integral-with-sources"></a>
### Fermionic Gaussian integral with sources
Appears in: C-068

The paired Berezin integral of e^{-theta^* B theta} is det B, and linear sources complete the square to det B e^{gamma^* B^{-1} gamma}.

$$ \widehat I_\theta e^{-\theta^*B\theta}=\det B,\qquad \widehat I_\theta e^{-\theta^*B\theta+\gamma^*\theta+\theta^*\gamma}=\det B\,e^{\gamma^*B^{-1}\gamma} $$

Audit: Correct in the C-067 convention. Checked by: exact WL computation (exterior algebra, n = 2, complex non-Hermitian B, every Grassmann coefficient matches).
QF code (exact): none. no fermionic layer.

<a id="pfaffian-gaussian-integral"></a>
### Pfaffian Gaussian integral
Appears in: C-069

A general Grassmann Gaussian integrates to a Pfaffian with the block-to-paired sign (-1)^{n(n-1)/2}, and sources put the inverse matrix in the exponent.

$$ \widehat I_\gamma e^{\Phi+\gamma^*\xi-\xi^*\gamma}=(-1)^{n(n-1)/2}\operatorname{Pf}(\Lambda)\exp\!\left[\tfrac12\binom{\xi^*}{\xi}^{T}\Lambda^{-1}\binom{\xi^*}{\xi}\right],\quad \Phi=\tfrac12x^T\Lambda x $$

Audit: Correct (with the atlas's block-to-paired sign). Checked by: exact WL computation (exterior algebra, n = 2, complex K, L, M: all coefficients match; the example Phi = -gamma^* gamma gives integral 1 and block Pfaffian -1; Pf^2 = det).
QF code (exact): none. no fermionic layer.

<a id="fermionic-translation-coefficient-and-pfaffian-amplitude"></a>
### Fermionic translation coefficient and Pfaffian amplitude
Appears in: C-086

The translation coefficient of a quadratic fermion propagator is a Grassmann Gaussian in the inverse Cayley transform, with amplitude the supertrace, whose square is (-1)^n det(F_t - I) and which equals an ordered Pfaffian.

$$ U_{t,w}(\gamma)=a_t\exp\!\left[-\tfrac14\binom\gamma{\gamma^*}^{T}J_F(F_t+I)(F_t-I)^{-1}\binom\gamma{\gamma^*}\right],\ a_t=\operatorname{Str}\hat U_t,\ a_t^2=(-1)^n\det(F_t-I),\ a_t=(-1)^{n(n-1)/2}p_t\operatorname{Pf}(-2J_FC_t) $$

Audit: Correct. One mode: Str U = 2i sin(theta), the Grassmann coefficient of Str(U T(gamma)) is -cos(theta), the exponent gives (i/2) cot(theta) gamma^* gamma, det(F - I) = 4 sin^2(theta), and the Pfaffian form gives cos(theta) 2i tan(theta). Checked by: hand derivation.
QF code (exact): none. no fermionic layer.

<a id="commuting-gaussian-integral-and-fresnel-factor"></a>
### Commuting Gaussian integral and Fresnel factor
Appears in: C-091

The real Gaussian integral of e^{-x^T M x/2 + j^T x} is (2 pi)^{m/2} det(M)^{-1/2} e^{j^T M^{-1} j/2}, and its oscillatory limit carries the Fresnel phase e^{i pi sgn(A)/4}.

$$ \int_{\mathbb R^m}e^{-x^TMx/2+j^Tx}dx=(2\pi)^{m/2}\det(M)^{-1/2}e^{j^TM^{-1}j/2},\qquad \text{Fresnel factor } e^{i\pi\operatorname{sgn}(A)/4}|\det A|^{-1/2} $$

Conditions: M symmetric, Re M > 0; root positive on positive M.
Audit: Correct (standard). Checked by: hand derivation.
QF code (exact): none. A classical integral; QF's second-quantization layer does not apply.

<a id="thermal-states-and-traces"></a>
## Thermal states and traces

<a id="traces-of-second-quantized-one-particle-operators"></a>
### Traces of second-quantized one-particle operators
Appears in: X-007; B-107

The Fock-space trace of the second quantization of a one-particle matrix $M$ is the inverse determinant of $I-M$ for bosons and the determinant of $I+M$ for fermions.

$$ \operatorname{Tr}_{\mathcal F_B}\widehat\Gamma(M)=\det(I_N-M)^{-1},\qquad \operatorname{Tr}_{\mathcal F_F}\widehat\Gamma(M)=\det(I_N+M) $$

Conditions: $\widehat\Gamma(e^h)=e^{\hat a^\dagger h\hat a}$; bosons need convergence: spectral radius of $M$ below 1 (the atlas gives $0\le M<I_N$ for the positive case); infinite modes need trace class and $\|\hat M\|<1$ for bosons.
Audit: Correct. The atlas's "equivalently $\|M\|<1$" is only an equivalence for positive $M$ (see Errors found). Checked by: QF exact (single mode, with the convergence condition returned); exact WL computation for two bosonic modes: with $\langle z|\widehat\Gamma(M)|w\rangle=e^{\bar z^TMw}$ between unnormalized coherent states (the previous entry), the trace is the Gaussian integral of $e^{-\bar z^T(I-M)z}$ over $\mathbb C^2$ divided by $\pi^2$, evaluated exactly for a rational non-Hermitian $M$ with spectral radius below 1: it equals $1/\det(I_N-M)$; exact fermion matrices (3 modes, same exact $h$ as above): $\operatorname{Tr}e^{\hat c^\dagger h\hat c}-\det(I_N+e^h)$ is exactly zero.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
Sum[BosonicMatrixElement[{n, n}, Exp[lam ad ** a]], {n, 0, Infinity}, GenerateConditions -> True]
```
Output: `ConditionalExpression[-(-1 + E^lam)^(-1), E^Re[lam] < 1]`. Matches: $1/(1-e^\lambda)$ with the convergence condition $|e^\lambda|<1$.

**Also stated in B-107 (Free Fock-space traces of a second-quantized one-particle map):**

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

<a id="bosonic-trace"></a>
### Bosonic trace
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

<a id="fermionic-trace"></a>
### Fermionic trace
Appears in: B-040

The fermionic trace is a Berezin integral of the kernel with weight e^{-aa*}, or of the normal symbol with the doubled weight e^{-2aa*}.

$$ \operatorname{Tr}\hat A=\int K_{\hat A}(a^*,a)e^{-aa^*}\Pi da\,da^* =\int N_{\hat A}(a^*,a)e^{-2aa^*}\Pi da\,da^* $$

Audit: Correct. Checked by: hand derivation (one mode): N = 1 gives int (1 - 2 a a*) = 2 = Tr I; N = a*a gives 1 = Tr N; odd symbols give 0; the formula is multiplicative over modes.
QF code (exact): none. No fermionic layer.

<a id="normalized-fermionic-trace"></a>
### Normalized fermionic trace
Appears in: B-041

Rescaling the symbol's arguments by 1/sqrt 2 and using the single weight gives the trace divided by 2^n, the normalized trace that is 1 on the identity.

$$ \tau(\hat A)=\int N_{\hat A}(a^*/\sqrt2,a/\sqrt2)e^{-aa^*}\Pi da\,da^*=2^{-n}\operatorname{Tr}\hat A $$

Audit: Correct. Checked by: hand derivation (one mode): N = 1 gives 1, N = a*a gives 1/2.
QF code (exact): none. No fermionic layer.

<a id="grand-partition-functions"></a>
### Grand partition functions
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

<a id="coherent-symbol-bounds-on-exponential-traces-berezin-lieb-golden-thompson"></a>
### Coherent-symbol bounds on exponential traces (Berezin-Lieb, Golden-Thompson)
Appears in: P-016

The trace of a Gibbs exponential is bounded below by the integral of $e^{-tQ}$ and above by the integral of $e^{-tP}$; Golden-Thompson bounds the trace of an exponential of a sum by the trace of the product.

$$ \int d\mu(\alpha)e^{-tQ_{\hat H}(\alpha)}\le\operatorname{Tr}e^{-t\hat H}\le\int d\mu(\alpha)e^{-tP_{\hat H}(\alpha)},\qquad \operatorname{Tr}e^{-t(\hat A+\hat B)}\le\operatorname{Tr}(e^{-t\hat A}e^{-t\hat B}) $$

Conditions: self-adjoint, bounded below, trace class, convergent symbol integrals; $d\mu=d^2\alpha/\pi$.
Audit: Correct. Checked by: hand derivation for the lower bound (Jensen: $e^{-t\langle\alpha|\hat H|\alpha\rangle}\le\langle\alpha|e^{-t\hat H}|\alpha\rangle$, then the resolution of identity); upper bound and Golden-Thompson cited as standard theorems.
QF code (exact): none. QF has no operator-inequality or trace-bound function.

<a id="weyl-characteristic-function-of-a-number-diagonal-state"></a>
### Weyl characteristic function of a number-diagonal state
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

<a id="quasiclassical-partition-function-bounds"></a>
### Quasiclassical partition-function bounds
Appears in: P2-032

For lattice spin systems the quantum partition function is squeezed between classical ones at $\beta$ and at a slightly larger $\beta$.

$$ Z_{\rm cl}(\beta)\le Z_Q^{(N)}(\beta)\le Z_{\rm cl}(\beta(1+a/N)) $$

Audit: Unchecked: Simon's theorem cited from the source; the lower bound is the Berezin-Lieb lower bound applied to a product of spin coherent states.
QF code (exact): none. QF has no partition-function bounds.

<a id="thermal-oscillator-master-equation"></a>
### Thermal oscillator master equation
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

<a id="displaced-thermal-density-kernel"></a>
### Displaced thermal density kernel
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

<a id="populations-of-a-displaced-thermal-state"></a>
### Populations of a displaced thermal state
Appears in: P2-044

The populations of a displaced thermal state are a geometric factor times a Gaussian times a Laguerre polynomial of negative argument.

$$ W_n=\frac{\mu^n}{(1+\mu)^{n+1}}e^{-|\beta|^2/(1+\mu)}L_n\left(-\frac{|\beta|^2}{\mu(1+\mu)}\right) $$

Audit: Correct, including the Poisson limit $\mu\to0$. Checked by: exact WL computation: $n!$ times the $z_1^n\bar z_2^n$ Taylor coefficient of the kernel $R$ (QF exact in the previous entry) equals the stated $W_n$ for $n\le6$ with symbolic $\mu,\beta$.
QF code (exact): none. `BosonicMatrixElement` handles $\hat D(\beta)$ or a function of $\hat N$ alone, not their product; the populations are read off the exact kernel of the previous entry.

<a id="free-oscillator-propagator-and-thermal-trace-in-coherent-states"></a>
### Free oscillator propagator and thermal trace in coherent states
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

<a id="stable-oscillator-proper-time-trace"></a>
### Stable oscillator proper-time trace
Appears in: P2-048

The trace of $e^{-s(2\hat N+1)}$ is $1/(2\sinh s)$ for $s>0$; the source's diagonalizing $\hat U$ is a complex continuation, not a unitary.

$$ \operatorname{Tr}e^{-s(2\hat N+\hat I)}=[2\sinh s]^{-1},\ s>0;\qquad \hat U=e^{i\pi(\hat b^{\dagger2}-\hat b^2)/8}\ \text{(not unitary)} $$

Audit: Correct, with the atlas's qualification that $\hat U$ has a Hermitian generator and is a complex canonical continuation. Checked by: exact WL computation (the sum minus $1/(2\sinh s)$ simplifies to 0); hand check that $i(\hat b^{\dagger2}-\hat b^2)$ is Hermitian.
QF code (exact): none. A trace of a function of $\hat N$ is an ordinary series; QF adds nothing.

<a id="berezin-lieb-bounds-for-spin"></a>
### Berezin-Lieb bounds for spin
Appears in: C-043

The trace of a convex function of a spin Hamiltonian is bounded below by the integral over its covariant symbol and above by the integral over its contravariant symbol; for J_l these symbols are j n_l and (j+1) n_l.

$$ \int\chi(H_c)d\mu\le\operatorname{Tr}\chi(\hat H)\le\int\chi(H^c)d\mu,\qquad (\hat J_\ell)_c=jn_\ell,\ (\hat J_\ell)^c=(j+1)n_\ell $$

Audit: Correct (standard Berezin-Lieb inequality and spin symbols). Checked by: hand derivation for the covariant symbol (<j,-j| rotated J_l |j,-j> = j n_l); contravariant value and the inequality from the standard result, not recomputed.
QF code (exact): none. Spin symbols; QF's second-quantization layer does not apply.

<a id="gaussian-density-operators-purity-thermal-factors"></a>
### Gaussian density operators, purity, thermal factors
Appears in: C-055

A centered Gaussian Weyl symbol is a density operator when V + iJ/2 >= 0; its purity is 2^{-n} det(V)^{-1/2}, the oscillator partition function is a product of 1/(2 sinh) factors, and each thermal mode has covariance (hbar/2) coth(beta hbar omega/2).

$$ V+iJ/2\ge0,\quad \operatorname{Tr}\hat\rho^2=2^{-n}(\det V)^{-1/2},\quad \operatorname{Tr}e^{-\beta\hat H}=\prod_j[2\sinh(\beta\hbar\omega_j/2)]^{-1},\quad V_j=\frac\hbar2\coth(\beta\hbar\omega_j/2)I_2 $$

Audit: Correct (vacuum V = I/2 gives purity 1; thermal purity tanh(beta omega/2) = 1/(2 nbar + 1)). Checked by: QF exact (number-basis diagonals of Q^2, P^2 and (QP+PQ)/2 for symbolic n) plus exact WL sums over the thermal weights (hbar = 1): partition function, purity, and covariance; the positivity criterion by hand (standard).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2]; p = -I (a - ad)/Sqrt[2];
diag = Simplify[{BosonicMatrixElement[{n, n}, q ** q], BosonicMatrixElement[{n, n}, p ** p],
    BosonicMatrixElement[{n, n}, (q ** p + p ** q)/2]}, n >= 0 && Element[n, Integers]];
z[b_] := Sum[Exp[-b w (n + 1/2)], {n, 0, Infinity}];
cov = Sum[Exp[-b w (n + 1/2)] diag, {n, 0, Infinity}]/z[b];
{diag, Simplify[{z[b] - 1/(2 Sinh[b w/2]), z[2 b]/z[b]^2 - Tanh[b w/2],
   z[2 b]/z[b]^2 - 1/(2 (1/(Exp[b w] - 1)) + 1), cov - {Coth[b w/2]/2, Coth[b w/2]/2, 0}}, b > 0 && w > 0]}
```
Output: `{{1/2 + n, 1/2 + n, 0}, {0, 0, 0, {0, 0, 0}}}`. Matches: the diagonals are n + 1/2, n + 1/2, 0, and all four differences are 0.

<a id="fermionic-partition-determinant"></a>
### Fermionic partition determinant
Appears in: C-079

Traces of fermionic quadratic exponentials are determinants: det(I + e^h) for the trace and det(I - e^h) for the supertrace.

$$ \operatorname{Tr}e^{\sum_jz_j\hat c_j^\dagger\hat c_j}=\prod_j(1+e^{z_j}),\quad \operatorname{Tr}e^{\hat c^\dagger h\hat c}=\det(I_n+e^h),\quad \operatorname{Str}e^{\hat c^\dagger h\hat c}=\det(I_n-e^h) $$

Audit: Correct. Checked by: exact fermion matrices (3 modes, complex non-Hermitian h).
QF code (exact): none. no fermionic layer.

<a id="thermal-state-as-a-number-exponential"></a>
### Thermal state as a number exponential
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

<a id="thermal-factorial-moments-and-number-variance"></a>
### Thermal factorial moments and number variance
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

<a id="thermal-characteristic-functions-and-phase-space-densities"></a>
### Thermal characteristic functions and phase-space densities
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

<a id="reduced-twin-beam-state-is-thermal-entanglement-entropy"></a>
### Reduced twin-beam state is thermal; entanglement entropy
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

<a id="thermal-state-in-normally-ordered-exponential-form"></a>
### Thermal state in normally ordered exponential form
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

<a id="amplitude-damping-kraus-form-and-matrix-element-solution"></a>
### Amplitude damping: Kraus form and matrix-element solution
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

<a id="damping-of-coherent-dyads-and-of-the-normal-characteristic-function"></a>
### Damping of coherent dyads and of the normal characteristic function
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

<a id="coherent-trace-formula-and-schmidt-spectrum"></a>
### Coherent trace formula and Schmidt spectrum
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

<a id="wick-dyson-and-auxiliary-fields"></a>
## Wick, Dyson and auxiliary fields

<a id="bound-for-exponentials-of-bounded-perturbations"></a>
### Bound for exponentials of bounded perturbations
Appears in: B-068

If A is a self-adjoint operator plus a bounded perturbation T, then e^{itA} grows at most like e^{|t| ||T||}.

$$ \|e^{it\hat A}\|\le e^{|t|\|\hat T\|},\qquad \hat A=\hat S+\hat T,\ \hat S=\hat S^* $$

Audit: Correct. Checked by: hand derivation: in the interaction-picture series the n-th term is bounded by (|t| ||T||)^n / n!.
QF code (exact): none. An operator-norm bound in a general Hilbert space; no second-quantization call is involved.

<a id="pair-potential-on-the-n-particle-sector"></a>
### Pair potential on the n-particle sector
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

<a id="dyson-series-with-the-e-ith-convention"></a>
### Dyson series with the e^{+itH} convention
Appears in: B-094

With Berezin's sign convention e^{+itH}, the interaction-picture propagator is the time-ordered exponential with powers of +i.

$$ \hat U(t,\tau)=\sum_{n\ge0}\frac{i^n}{n!}\int_\tau^t\cdots\int_\tau^t T[\hat V(t_1)\cdots \hat V(t_n)]\,dt_1\cdots dt_n $$

Conditions: U(t, tau) = e^{-itH0} e^{i(t-tau)H} e^{i tau H0}, V(t) = e^{-itH0} V e^{itH0}; formal.
Audit: Correct. Checked by: hand derivation: dU/dt = i V(t) U with U(tau, tau) = I.
QF code (exact): none. A formal time-ordered series; no exact QF call evaluates it.

<a id="free-time-dependence-and-contraction-kernel"></a>
### Free time dependence and contraction kernel
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

<a id="retarded-contraction-and-the-two-field-time-ordered-product"></a>
### Retarded contraction and the two-field time-ordered product
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

<a id="wick-s-theorem-as-an-exponential-of-the-contraction-operator"></a>
### Wick's theorem as an exponential of the contraction operator
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

<a id="normal-symbol-of-interaction-picture-evolution"></a>
### Normal symbol of interaction-picture evolution
Appears in: B-100

Applying Wick's theorem term by term to the Dyson series gives the normal symbol of the interaction-picture propagator as e^{D} acting on the exponential of the integrated interaction symbol.

$$ N_{\hat U,1}(t,\tau\mid a^*,a)=e^{\widehat{\mathcal D}} \exp\!\left[i\int_\tau^t V(s\mid a^*,a)ds\right] $$

Audit: Correct as a formal identity: B-094 plus B-099 term by term. Checked by: hand derivation.
QF code (exact): none. A formal series over arbitrary interactions; no exact QF call evaluates it.

<a id="contractions-acting-on-an-exponential-source"></a>
### Contractions acting on an exponential source
Appears in: B-101

An exponential of linear Grassmann sources is an eigenfunction of the contraction operator, picking up the factor e^{f* Delta f}.

$$ e^{\widehat{\mathcal D}}e^{f^*a+a^*f}=e^{f^*\Delta f+a^*f+f^*a} $$

Audit: Correct under the reading of B-099 that gives the fermionic sign; fails under the other reading. Checked by: exact WL computation (explicit Grassmann algebra, two modes, complex Delta, Grassmann sources: exact match).
QF code (exact): none. No fermionic layer.

<a id="determinant-identity-for-a-quadratic-symbol"></a>
### Determinant identity for a quadratic symbol
Appears in: B-102

Contracting a fermionic quadratic exponential gives det(I + K Delta) times a quadratic exponential with the dressed kernel (K^{-1} + Delta)^{-1}.

$$ e^{\widehat{\mathcal D}}e^{-a^*\hat Ka} =\det(\hat I+\hat K\hat\Delta)\exp[-a^*(\hat K^{-1}+\hat\Delta)^{-1}a] $$

Audit: Correct (fermionic). Checked by: exact WL computation (explicit Grassmann algebra, two modes, complex non-commuting K and Delta: exact match under the correct reading).
QF code (exact): none. No fermionic layer.

<a id="auxiliary-field-decoupling-of-a-quartic-term"></a>
### Auxiliary-field decoupling of a quartic term
Appears in: B-103

Factorizing the quartic kernel through an auxiliary label lets a commuting Gaussian integral over z, z* turn the quartic exponential into a quadratic one.

$$ \int V_1(p_1,q_1\mid s)V_2(s\mid p_2,q_2)ds=iV(p_1,p_2\mid q_1,q_2) $$

Audit: Correct in mechanism: int e^{-z*z + u z + z* v} d^2z/pi = e^{uv} turns the factorized kernel into a quadratic coupling. Signs and the factor i in the full Appendix formula are Unchecked (not transcribed). Checked by: hand derivation: complete the square, z -> z + v, z* -> z* + u.
QF code (exact): none. The full Appendix formula is not transcribed, so there is nothing for QF to reproduce.

<a id="retarded-loops-vanish"></a>
### Retarded loops vanish
Appears in: B-104

Closed loops of retarded contractions with a time-local kernel vanish, so the auxiliary determinant equals 1.

$$ \operatorname{Tr}(\hat R\hat\Delta)^n=0,\qquad \det(\hat I+\hat R\hat\Delta)=1 $$

Audit: Correct as a formal result. Checked by: hand derivation: a closed loop needs t_1 > t_2 > ... > t_n > t_1, impossible, and n = 1 vanishes by Delta(t, t) = 0.
QF code (exact): none. A formal statement about continuum time kernels; no exact QF call applies.

<a id="auxiliary-field-representation-of-the-propagator"></a>
### Auxiliary-field representation of the propagator
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

<a id="second-order-dyson-elimination-of-a-fast-oscillating-coupling"></a>
### Second-order Dyson elimination of a fast oscillating coupling
Appears in: G-102

To second order a fast oscillating coupling leaves a secular term linear in t proportional to [A, A^dagger], which defines the effective Hamiltonian.

$$ \hat U_I\approx \hat I-\frac{ig^2t}{\Delta}[\hat A,\hat A^\dagger],\qquad\hat H_{\rm eff}=\frac{\hbar g^2}{\Delta}[\hat A,\hat A^\dagger]. $$

Audit: Correct (sign of the secular term). Checked by: exact WL computation (script wl/100.wl): the time-ordered double integrals of e^{+-i Delta(t1 - t2)} have secular parts +-i t/Delta and those of e^{+-i Delta(t1 + t2)} are bounded, giving -i g^2 t [A, A^dagger]/Delta.
QF code (exact): none. a time-ordered perturbative statement; QF has no Dyson-series tool.

<a id="nonlinear-exponentials-number-shift-kerr-cats"></a>
## Nonlinear exponentials (number shift, Kerr, cats)

<a id="number-commutators-and-free-rotation-of-the-field"></a>
### Number commutators and free rotation of the field
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

<a id="number-shift-rules-for-functions-of-the-number-operator"></a>
### Number-shift rules for functions of the number operator
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

<a id="no-hermitian-phase-operator-canonically-conjugate-to-number"></a>
### No Hermitian phase operator canonically conjugate to number
Appears in: G-011

A Hermitian phase operator canonically conjugate to the number operator cannot exist, because its diagonal number-basis elements would have to satisfy 0 = i.

$$ (m-n)\langle m|\hat\phi|n\rangle=i\delta_{mn} $$

Audit: Correct. Checked by: hand derivation (the m = n element gives 0 = i).
QF code (exact): none. an impossibility statement; there is nothing to compute.

<a id="susskind-glogower-one-sided-shift"></a>
### Susskind-Glogower one-sided shift
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

<a id="exact-resonant-jaynes-cummings-propagator"></a>
### Exact resonant Jaynes-Cummings propagator
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

<a id="dispersive-hamiltonian-as-conditional-number-rotation"></a>
### Dispersive Hamiltonian as conditional number rotation
Appears in: G-040

In the dispersive limit the atom picks a number-dependent rotation of the field: e^{i chi t n} in the ground state and e^{-i chi t(n+1)} in the excited state.

$$ \hat H_{\rm eff}=\hbar\chi\bigl(|e\rangle\langle e|+\hat n\hat\sigma_3\bigr),\quad\hat U_{\rm eff}=|g\rangle\langle g|\,e^{i\chi t\hat n}+|e\rangle\langle e|\,e^{-i\chi t(\hat n+\hat I)}. $$

Audit: Correct (excited-state sign as in the atlas). Checked by: hand derivation (sigma_3 = -1 on g gives -chi n; +1 on e gives chi(n+1)).
QF code (exact): none. the propagator is diagonal in the atom-number basis; the field rotation it implies is the G-004 check.

<a id="raman-number-controlled-atomic-rotation"></a>
### Raman number-controlled atomic rotation
Appears in: G-041

The Raman coupling lambda n sigma_1 rotates the atom by an angle proportional to the photon number, equivalently it rotates the field in opposite senses on the two sigma_1 eigenstates.

$$ \hat U_I(t)=\cos(\lambda t\hat n)\hat I_A-i\sin(\lambda t\hat n)\hat\sigma_1=|+\rangle_x\langle+|e^{-i\lambda t\hat n}+|-\rangle_x\langle-|e^{i\lambda t\hat n}. $$

Audit: Correct. Checked by: exact WL computation (script wl/041.wl): on each number sector the exact 2x2 exponential equals both stated forms, symbolic m.
QF code (exact): none. no spin-boson layer; the operator acts as an exact 2x2 matrix in each number sector.

<a id="phase-shifter-and-parity-exponential"></a>
### Phase shifter and parity exponential
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

<a id="self-kerr-evolution-of-a-coherent-state"></a>
### Self-Kerr evolution of a coherent state
Appears in: G-076

Self-Kerr evolution multiplies each Fock amplitude by e^{-i K t m^2}; at t = pi/K this is the parity operator.

$$ \hat U_K(t)=e^{-iKt \hat n^2},\quad \hat U_K(t)|\alpha\rangle=e^{-|\alpha|^2/2}\sum_{m\ge0}\frac{\alpha^m}{\sqrt{m!}}e^{-iKtm^2}|m\rangle,\quad \hat U_K(\pi/K)=\hat\Pi. $$

Audit: Correct. Checked by: exact WL computation (script wl/075.wl): e^{-i pi m^2} = (-1)^m for integer m.
QF code (exact): none. BosonicMatrixElement[{m, n}, Exp[-I K t nn ** nn]] returns unevaluated (no rule for the square of a^dagger a in an exponent).

<a id="half-revival-cat-identity"></a>
### Half-revival cat identity
Appears in: G-077

At the half revival the Kerr phase is 1 on even and -i on odd photon numbers, so the operator is a combination of identity and parity and turns a coherent state into a cat.

$$ e^{-i\pi \hat n^2/2}=\frac{e^{-i\pi/4}}{\sqrt2}(\hat I+i\hat\Pi),\quad e^{-i\pi \hat n^2/2}|\alpha\rangle=\frac{e^{-i\pi/4}}{\sqrt2}(|\alpha\rangle+i|-\alpha\rangle). $$

Audit: Correct; confirms the atlas's correction of printed Eq. (7.142). Checked by: exact WL computation (script wl/076.wl): e^{-i pi m^2/2} is 1 for m = 2k and -i for m = 2k+1, the same values as e^{-i pi/4}(1 + i(-1)^m)/sqrt(2); the printed +i on odd m is wrong.
QF code (exact): none. same missing n^2 rule as the previous entry.

<a id="normally-ordered-kerr-differs-by-a-rotation"></a>
### Normally ordered Kerr differs by a rotation
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

<a id="number-function-conjugation-and-kerr-heisenberg-solution"></a>
### Number-function conjugation and Kerr Heisenberg solution
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

<a id="kerr-coherent-state-moments-and-quadrature-variance"></a>
### Kerr coherent-state moments and quadrature variance
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

<a id="cross-kerr-conditional-rotation"></a>
### Cross-Kerr conditional rotation
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

<a id="coherent-state-cat-normalization-and-a-2-eigenvalue"></a>
### Coherent-state cat: normalization and a^2 eigenvalue
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

<a id="cat-photon-number-distribution"></a>
### Cat photon-number distribution
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

<a id="no-jump-number-exponential-on-a-coherent-state"></a>
### No-jump number exponential on a coherent state
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

<a id="effective-kerr-hamiltonian-from-an-off-resonant-parametric-coupling"></a>
### Effective Kerr Hamiltonian from an off-resonant parametric coupling
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

<a id="fock-matrix-elements-and-generating-functions"></a>
## Fock matrix elements and generating functions

<a id="fock-basis-and-ladder-action"></a>
### Fock basis and ladder action
Appears in: P-005; G-005; B-004

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

**Also stated in G-005 (Fock ladder action and resolution of identity):**

The ladder operators move one step on the number basis with square-root weights, the n-photon state is the normalized n-th power of a^dagger on the vacuum, and the number states resolve the identity.

$$ \hat a|n\rangle=\sqrt n|n-1\rangle,\quad \hat a^\dagger|n\rangle=\sqrt{n+1}|n+1\rangle,\quad |n\rangle=\frac{(\hat a^\dagger)^n}{\sqrt{n!}}|0\rangle,\quad\sum_{n\ge0}|n\rangle\langle n|=\hat I. $$

Audit: Correct. Checked by: QF exact (symbolic Fock indices; norm of (a^dagger)^3|0>/sqrt(3!)); hand derivation for the completeness sum.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
{BosonicMatrixElement[{m, n}, a], BosonicMatrixElement[{m, n}, ad], BosonicVEV[a ** a ** a ** ad ** ad ** ad, {a, ad}]/3!}
```
Output: `{Sqrt[n]*KroneckerDelta[m, -1 + n], Sqrt[m]*KroneckerDelta[-1 + m, n], 1}`. Matches.

**Also stated in B-004 (Fock vector from the vacuum):**

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

<a id="coordinate-to-bargmann-kernel"></a>
### Coordinate to Bargmann kernel
Appears in: P-013

The kernel that takes a position wavefunction to its Bargmann function is a Gaussian in $z$ and $q$.

$$ K(z,q)=(\pi\hbar)^{-1/4}\exp\left(-\frac{z^2}{2}+\sqrt{\frac2\hbar}zq-\frac{q^2}{2\hbar}\right) $$

Audit: Correct. Checked by: exact WL computation: $\sum_n z^n\psi_n(q)/\sqrt{n!}$ with Hermite functions, summed in closed form through the Hermite generating function $\sum_n t^nH_n(x)/n!=e^{2tx-t^2}$ at $t=z/\sqrt2$, $x=q/\sqrt\hbar$, equals the stated kernel for symbolic $\hbar$.
QF code (exact): none. QF's second-quantization layer has no position representation of states.

<a id="displacement-matrix-elements-in-associated-laguerre-form"></a>
### Displacement matrix elements in associated Laguerre form
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

<a id="coherent-generating-kernel-of-the-displacement"></a>
### Coherent generating kernel of the displacement
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

<a id="displacement-matrix-elements-as-a-finite-sum"></a>
### Displacement matrix elements as a finite sum
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

<a id="analytic-density-kernel"></a>
### Analytic density kernel
Appears in: P2-040

The unnormalized coherent matrix element of a density operator is a generating function whose Taylor coefficients are the Fock matrix elements.

$$ R(\bar\alpha,\beta)=e^{(|\alpha|^2+|\beta|^2)/2}\langle\alpha|\hat\rho|\beta\rangle=\sum_{mn}\rho_{mn}\frac{\bar\alpha^m\beta^n}{\sqrt{m!n!}} $$

Audit: Correct (with $z_1\leftrightarrow\bar\alpha$, $\bar z_2\leftrightarrow\beta$). Checked by: hand derivation from the coherent expansion $e^{|\alpha|^2/2}|\alpha\rangle=\sum_n\alpha^n|n\rangle/\sqrt{n!}$.
QF code (exact): none. The displaced thermal entry below exercises this kernel with QF.

<a id="laguerre-generating-function-in-the-lower-index"></a>
### Laguerre generating function in the lower index
Appears in: P2-053

The generating function of $L_n^k$ in the lower index is $(1-t)^{-k-1}e^{-tx/(1-t)}$.

$$ \sum_{n\ge0}t^nL_n^k(x)=(1-t)^{-k-1}e^{-tx/(1-t)},\quad |t|<1 $$

Audit: Correct. Checked by: exact WL computation, series coefficients through $t^6$ with symbolic $k$ and $x$ all equal.
QF code (exact): none. Not an operator identity.

<a id="laguerre-generating-function-along-the-diagonal-m-n"></a>
### Laguerre generating function along the diagonal m-n
Appears in: P2-053

Along the diagonal $m-n$ the Laguerre generating function is $(1+t)^me^{-tx}$.

$$ \sum_{n\ge0}t^nL_n^{m-n}(x)=(1+t)^m e^{-tx} $$

Audit: Correct. Checked by: exact WL computation, $m=0..6$ through order $t^{10}$, all differences 0.
QF code (exact): none. Not an operator identity.

<a id="generating-function-of-displacement-transition-probabilities-prefactor-corrected"></a>
### Generating function of displacement transition probabilities (prefactor corrected)
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

<a id="displacement-matrix-elements-as-wirtinger-derivatives"></a>
### Displacement matrix elements as Wirtinger derivatives
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

<a id="independent-mode-ccr-and-multimode-fock-states"></a>
### Independent-mode CCR and multimode Fock states
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

<a id="fock-matrix-elements-from-the-coherent-diagonal-symbol"></a>
### Fock matrix elements from the coherent diagonal symbol
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

<a id="multiphoton-jaynes-cummings-couplings"></a>
### Multiphoton Jaynes-Cummings couplings
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

<a id="fock-coefficients-of-a-displaced-squeezed-vacuum-hermite-form"></a>
### Fock coefficients of a displaced squeezed vacuum (Hermite form)
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

<a id="infinite-modes-and-domains"></a>
## Infinite modes and domains

<a id="smeared-ccr-and-car"></a>
### Smeared CCR and CAR
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

<a id="bosonic-quadratures"></a>
### Bosonic quadratures
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

<a id="weighted-sufficient-conditions-for-a-quadratic-normal-form"></a>
### Weighted sufficient conditions for a quadratic normal form
Appears in: B-078

When C is multiplication by c(x), finiteness of these two weighted integrals suffices for the formal quadratic expression to be reduced to a quadratic normal form.

$$ \int\frac{|f(x)|^2}{|c(x)-z|}dx<\infty,\quad \iint\frac{|A(x,y)|^2}{|c(x)+c(y)-z|}dxdy<\infty $$

Audit: Unchecked: infinite-mode sufficient condition from a source theorem; no finite-dimensional content. Checked by: prior audit only.
QF code (exact): none. A continuum domain condition; no finite second-quantization content.

<a id="square-root-perturbation-stays-in-the-ideal"></a>
### Square-root perturbation stays in the ideal
Appears in: B-087

The correction X = (C^2 + K)^{1/2} - C solves X^2 + CX + XC = K, which the source uses to show X is as small (trace class, Hilbert-Schmidt) as K.

$$ \hat X=(\hat C^2+\hat K)^{1/2}-\hat C,\qquad \hat X^2+\hat C\hat X+\hat X\hat C=\hat K $$

Audit: Correct for the identity. The ideal statement is a source theorem. Checked by: hand derivation: (X + C)^2 = C^2 + K.
QF code (exact): none. An operator identity; no second-quantization call is involved.

<a id="scattering-unitary-in-the-thirring-model"></a>
### Scattering unitary in the Thirring model
Appears in: B-092

The Thirring model's scattering operator is a phase times the exponential of the product of the two zero-momentum charges.

$$ \hat S=\theta e^{-\pi igF(0)\hat Q_1'(0)\hat Q_2'(0)} $$

Audit: Unchecked: continuum model result, not reconstructed. Checked by: prior audit only.
QF code (exact): none. A continuum fermionic model; no fermionic layer.

<a id="calogero-degeneracy-generating-function"></a>
### Calogero degeneracy generating function
Appears in: P2-016

The degeneracies of the Calogero model after removing the center of mass are generated by $\prod_{j=2}^N(1-t^j)^{-1}$.

$$ \prod_{j=2}^N(1-t^j)^{-1} $$

Audit: Correct by counting: symmetric polynomials in $N$ variables are generated by power sums of degree $1..N$; removing center-of-mass motion removes degree 1. Checked by: hand derivation.
QF code (exact): none. Not an operator identity.

<a id="canonical-commutator-of-position-and-momentum"></a>
### Canonical commutator of position and momentum
Appears in: C-001

Momentum and position along the same axis fail to commute by a fixed imaginary constant; along different axes they commute.

$$ [\hat P_j,\hat Q_k]=-i\hbar\delta_{jk}\hat I $$

Conditions: on the product domain; real phase-space labels.
Audit: Correct. Checked by: QF exact (one mode, hbar = 1).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2]; p = -I (a - ad)/Sqrt[2];
BosonicNormalOrder[p ** q - q ** p, {a, ad}]
```
Output: `-I`. Matches: [P, Q] = -i at hbar = 1; different modes commute because their field variables commute.

<a id="periodic-orbit-determinant-in-regularized-traces"></a>
### Periodic-orbit determinant in regularized traces
Appears in: C-036

Energy-localized traces of the propagator are asymptotically sums over periodic orbits, each weighted by |det(I - P_gamma)|^{-1/2} and a Maslov phase.

$$ \operatorname{Tr}[\chi(\hat H)\hat A\chi(\hat H)g((E\hat I-\hat H)/\hbar)]\ \sim\ \ldots\,|\det(I_{d_\gamma}-P_\gamma)|^{-1/2}e^{i\pi(\cdot)} $$

Audit: Unchecked: semiclassical stationary-phase asymptotics; no finite identity to test.
QF code (exact): none. QF has no semiclassical trace routine.

<a id="segal-fields"></a>
### Segal fields
Appears in: C-057

Smeared annihilation and creation operators satisfy the canonical commutator with the one-particle inner product, and the field operator is their real combination.

$$ [\hat a(f),\hat a^\dagger(g)]=\langle f,g\rangle\hat I,\qquad \hat\Phi(f)=\frac{\hat a(f)+\hat a^\dagger(f)}{\sqrt2} $$

Audit: Correct (standard; a(f) antilinear). Checked by: QF exact (two modes, symbolic f, g).
QF code (exact):
```wl
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
af[f_] := Conjugate[f[[1]]] a1 + Conjugate[f[[2]]] a2;
adf[f_] := f[[1]] a1d + f[[2]] a2d;
{BosonicNormalOrder[af[{f1, f2}] ** adf[{g1, g2}] - adf[{g1, g2}] ** af[{f1, f2}], {a1, a1d, a2, a2d}, "Scalars" -> {f1, f2, g1, g2}],
 Conjugate[{f1, f2}] . {g1, g2}}
```
Output: `{g1*Conjugate[f1] + g2*Conjugate[f2], g1*Conjugate[f1] + g2*Conjugate[f2]}`. Matches: [a(f), a^dag(g)] = conj(f).g = <f, g>.

<a id="infinite-mode-displacement-exponential-vector-product-overlap"></a>
### Infinite-mode displacement, exponential vector, product, overlap
Appears in: C-059 (one-mode content identical to C-003, C-006, C-007, C-008)

On Fock space over any one-particle space, the Weyl operator factors in normal order, composes with the phase e^{-i Im<f,g>}, and coherent vectors overlap as e^{-|f-g|^2/2 + i Im<f,g>}; occupation numbers are Poisson with mean ||f||^2.

$$ \hat T(f)=e^{-\|f\|^2/2}e^{\hat a^\dagger(f)}e^{-\hat a(f)},\quad \hat T(f)\hat T(g)=e^{-i\Im\langle f,g\rangle}\hat T(f+g),\quad \langle\Psi(f),\Psi(g)\rangle=e^{-\|f-g\|^2/2+i\Im\langle f,g\rangle} $$

Audit: Correct (alpha conj(beta) - conj(alpha) beta = -2i Im(conj(alpha) beta)). Checked by: QF exact (one mode, symbolic f, g; several modes give the product over an orthonormal basis since different modes commute) and hand derivation for the Poisson law (|<n|Psi(f)>|^2 = e^{-|f|^2}|f|^{2n}/n!, from C-007).
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
d[x_] := Exp[x ad - Conjugate[x] a];
cpx = {f -> f1 + I f2, g -> g1 + I g2};
prod = BosonicExpSimplify[d[f] ** d[g]];
ov = BosonicMatrixElement[{f, g}, 1, "Basis" -> "Coherent"];
{prod, ov, Simplify[ComplexExpand[{(f Conjugate[g] - g Conjugate[f])/2 + I Im[Conjugate[f] g],
    (ov /. E^x_ :> x) - (-Abs[f - g]^2/2 + I Im[Conjugate[f] g])} /. cpx]]}
```
Output: `{E^(-(\[FormalA]1*Conjugate[f]) - \[FormalA]1*Conjugate[g] + (-(g*Conjugate[f]) + f*Conjugate[g])/2 + f*SuperDagger[\[FormalA]1] + g*SuperDagger[\[FormalA]1]), E^((-Abs[f]^2 - Abs[g]^2)/2 + g*Conjugate[f]), {0, 0}}`. Matches: the product scalar equals -i Im<f,g> and the overlap exponent equals -|f-g|^2/2 + i Im<f,g> (both differences 0).

<a id="collective-broadband-commutator"></a>
### Collective broadband commutator
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

<a id="continuum-ccr-and-first-order-pair-creation-kernel"></a>
### Continuum CCR and first-order pair-creation kernel
Appears in: G-105

Continuum modes have a delta-function commutator, and to first order a pump photon splits into a pair constrained by energy and momentum conservation.

$$ [\hat a(\mathbf k),\hat a^\dagger(\mathbf k')]=\delta^{(3)}(\mathbf k-\mathbf k')\hat I,\qquad\int d^3k_s\,d^3k_i\,\delta(\omega_P-\omega_s-\omega_i)\delta^{(3)}(\mathbf k_P-\mathbf k_s-\mathbf k_i)\hat a_s^\dagger(\mathbf k_s)\hat a_i^\dagger(\mathbf k_i)|0\rangle. $$

Audit: Correct as a distributional statement; the smeared form [a(f), a^dagger(g)] = <f, g> I is right. Checked by: hand derivation.
QF code (exact): none. continuum modes; QF works with finitely many modes.

<a id="lie-group-decompositions-and-generalized-coherent-states"></a>
## Lie-group decompositions and generalized coherent states

<a id="su-2-disentangling-compact-case"></a>
### SU(2) disentangling (compact case)
Appears in: X-004 (source records cited: P-025, P-029, C-041, C-048, G-046, P-041 to P-043, P2-025); P-025; C-041; fermion-pair application in C-087

A compact group element generated by $\hat J_\pm$ factors into raising, Cartan, and lowering exponentials with $\zeta=e^{i\phi}\tan r$, on the chart where $\cos r>0$.

$$ e^{\xi \hat J_+-\bar\xi \hat J_-}=e^{\zeta \hat J_+}e^{\log(1+|\zeta|^2)\hat J_0}e^{-\bar\zeta \hat J_-},\qquad \zeta=e^{i\phi}\tan r $$

Conditions: $[\hat J_0,\hat J_\pm]=\pm\hat J_\pm$, $[\hat J_+,\hat J_-]=2\hat J_0$, $\xi=re^{i\phi}$, $\cos r>0$. The chart condition is needed: for $\cos r<0$ the middle factor loses the sign of $\cos r$, which matters in half-integer representations.
Audit: Correct, including the chart remark. Checked by: exact WL computation: symbolic $r,\phi$ in the spin-1/2 and spin-1 matrices; the difference is zero for $0<r<\pi/2$; for $\pi/2<r<\pi$ the spin-1/2 right side is exactly minus the left side, while spin 1 still matches.
QF code (exact): none. QF has no abstract SU(2) disentangling.

**Also stated in P-025 (SU(2) disentangling):**

A compact rotation exponential factors into raising, diagonal and lowering exponentials with $\zeta=e^{i\phi}\tan r$, in either order, on the chart $\cos r>0$.

$$ e^{\xi \hat J_+-\bar\xi \hat J_-}=e^{\zeta \hat J_+}e^{\log(1+|\zeta|^2)\hat J_0}e^{-\bar\zeta \hat J_-}=e^{-\bar\zeta \hat J_-}e^{-\log(1+|\zeta|^2)\hat J_0}e^{\zeta \hat J_+},\quad \xi=re^{i\phi},\ \zeta=e^{i\phi}\tan r $$

Conditions: $\cos r>0$.
Audit: Correct on the stated chart. Checked by: exact WL computation in the faithful 2x2 matrices with symbolic $r,\phi$: both orders equal the left side for $0<r<\pi/2$, and for $\pi/2<r<\pi$ both equal minus the left side. That sign is the central element $-I$, which acts as $(-1)^{2j}$ in spin $j$, so the identity holds in every representation on the chart and fails for half-integer $j$ beyond it.
QF code (exact): none. QF's exact algebra layer has no su(2) generators.

**Also stated in C-041 (SU(2) normal factorization and spin coherent state):**

The SU(2) displacement factors into a raising exponential, a J_3 exponential and a lowering exponential with tan r, and acting on the lowest weight it gives the normalized spin coherent state.

$$ e^{\xi\hat J_+-\bar\xi\hat J_-}=e^{\zeta\hat J_+}e^{-2\log(\cos r)\hat J_3}e^{-\bar\zeta\hat J_-},\ \zeta=\frac{\xi}{r}\tan r;\qquad \hat G(\xi)|j,-j\rangle=(1+|\zeta|^2)^{-j}e^{\zeta\hat J_+}|j,-j\rangle $$

Conditions: cos r != 0.
Audit: Correct. Checked by: exact WL computation (j = 3/2 complete spin matrices, xi = 0.7 e^{0.4 i}: both identities hold to machine precision).
QF code (exact): none. Finite spin matrices; QF's second-quantization layer does not apply.

<a id="su-1-1-disentangling-noncompact-case"></a>
### SU(1,1) disentangling (noncompact case)
Appears in: X-004; P-029 (also referenced in P2-058)

The noncompact analogue uses $\zeta=e^{i\phi}\tanh r$ and $\log(1-|\zeta|^2)$, and needs no chart restriction.

$$ e^{\xi \hat K_+-\bar\xi \hat K_-}=e^{\zeta \hat K_+}e^{\log(1-|\zeta|^2)\hat K_0}e^{-\bar\zeta \hat K_-},\qquad \zeta=e^{i\phi}\tanh r $$

Conditions: $[\hat K_0,\hat K_\pm]=\pm\hat K_\pm$, $[\hat K_-,\hat K_+]=2\hat K_0$; no chart restriction.
Audit: Correct. With $\hat K_+=\hat a^{\dagger2}/2$, $\hat K_-=\hat a^2/2$, $\hat K_0=(\hat N+\tfrac12)/2$ and $\xi\to-\xi$ it is the X-002 factorization; with the two-mode $\hat K$'s it is X-003. Checked by: exact WL computation (the 2x2 non-unitary representation $\hat K_0=\sigma_z/2$, $\hat K_+=\sigma_+$, $\hat K_-=-\sigma_-$ with symbolic $r,\phi$: algebra and identity both exact zero) and QF exact for the single-mode boson realization.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
FullSimplify[BosonicExpOrder[Exp[r Exp[I ph] ad ** ad/2 - r Exp[-I ph] a ** a/2]], r > 0 && Element[ph, Reals]]
```
Output: `E^((E^(I*ph)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]*Tanh[r])/2)**Sech[r]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^(-1/2*(\[FormalA]1**\[FormalA]1*Tanh[r])/E^(I*ph))`. Matches: this is $e^{\zeta\hat K_+}(1-|\zeta|^2)^{\hat K_0}e^{-\bar\zeta\hat K_-}$, since $(1-\tanh^2r)^{(\hat N+1/2)/2}=\operatorname{sech}r^{\,\hat N+1/2}$.

**Also stated in P-029 (SU(1,1) disentangling (middle-exponent sign corrected)):**

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

<a id="su-2-commutators-ladder-action-and-holomorphic-realization"></a>
### su(2) commutators, ladder action and holomorphic realization
Appears in: P-024; C-038

The spin generators obey the su(2) commutators, act on $|j,m\rangle$ with the standard square-root ladder factors, and are realized on polynomials of degree at most $2j$ by first-order differential operators.

$$ [\hat J_0,\hat J_\pm]=\pm \hat J_\pm,\ [\hat J_+,\hat J_-]=2\hat J_0,\ \hat J_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}|j,m\pm1\rangle;\quad \hat J_-=\partial_z,\ \hat J_+=-z^2\partial_z+2jz,\ \hat J_0=z\partial_z-j $$

Audit: Correct. Checked by: exact WL computation: commutators of the differential operators on a symbolic $f(z)$ vanish identically, and the exact spin matrices for $j=1/2,1,3/2,2$ satisfy the algebra with zero residual.
QF code (exact): none. QF's exact algebra layer is bosonic; it has no abstract su(2) generators.

**Also stated in C-038 (SU(2) commutators, ladder matrices, Casimir):**

The spin-j ladder operators satisfy the su(2) relations, raise and lower m with the standard square-root factors, and give the Casimir j(j+1).

$$ [\hat J_3,\hat J_\pm]=\pm\hat J_\pm,\ [\hat J_+,\hat J_-]=2\hat J_3,\ \hat J_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}|j,m\pm1\rangle,\ \hat J^2=j(j+1)\hat I $$

Audit: Correct. Checked by: exact WL computation (j = 3/2, complete 4-by-4 spin matrices).
QF code (exact): none. Finite spin matrices; QF's second-quantization layer does not apply.

<a id="su-2-coherent-state-expansion"></a>
### SU(2) coherent state expansion
Appears in: P-026

The spin coherent state is the normalized exponential of $\zeta\hat J_+$ on the lowest weight, with binomial square-root amplitudes.

$$ |\zeta;j\rangle=(1+|\zeta|^2)^{-j}e^{\zeta \hat J_+}|j,-j\rangle=(1+|\zeta|^2)^{-j}\sum_{n=0}^{2j}\sqrt{\binom{2j}{n}}\zeta^n|j,-j+n\rangle $$

Audit: Correct. Checked by: exact WL computation with exact spin matrices ($j=3/2$ and 2, symbolic $\zeta$), and hand derivation for general $j$: $\hat J_+^n|j,-j\rangle=\sqrt{n!\,(2j)!/(2j-n)!}\,|j,-j+n\rangle$.
QF code (exact): none. No su(2) layer in QF.

<a id="su-2-coherent-overlap-and-annihilator"></a>
### SU(2) coherent overlap and annihilator
Appears in: P-026; P-026; C-042

Spin coherent states overlap by $(1+\bar\zeta\eta)^{2j}$ over the normalizations, and each is annihilated by $\hat J_-+2\zeta\hat J_0-\zeta^2\hat J_+$.

$$ \langle\zeta;j|\eta;j\rangle=\frac{(1+\bar\zeta\eta)^{2j}}{(1+|\zeta|^2)^j(1+|\eta|^2)^j},\qquad (\hat J_-+2\zeta \hat J_0-\zeta^2\hat J_+)|\zeta;j\rangle=0 $$

Audit: Correct. Checked by: exact WL computation ($j=3/2$ and 2, symbolic $\zeta,\eta$, zero residuals); hand derivation of the annihilator: $e^{\zeta\hat J_+}\hat J_-e^{-\zeta\hat J_+}=\hat J_-+2\zeta\hat J_0-\zeta^2\hat J_+$.
QF code (exact): none. No su(2) layer in QF.

**Also stated in P-026 (SU(2) resolution of identity):**

Spin coherent projectors integrated with the invariant sphere measure resolve the identity.

$$ \frac{2j+1}{\pi}\int\frac{d^2\zeta}{(1+|\zeta|^2)^2}|\zeta;j\rangle\langle\zeta;j|=\hat I $$

Audit: Correct. Checked by: exact WL computation; the diagonal element $(2j+1)\binom{2j}{n}\int_0^\infty t^n(1+t)^{-2j-2}dt$ simplifies to 1 for symbolic $j$ and $n$.
QF code (exact): none. No su(2) layer in QF.

**Also stated in C-042 (Spin coherent overlap and resolution of identity):**

Spin coherent states overlap as (1 + conj(eta) zeta)^{2j} over the norms, with modulus set by the angle between the two directions, and resolve the identity with the measure (2j+1)/pi (1+|zeta|^2)^{-2}.

$$ \langle\eta|\zeta\rangle=\frac{(1+\bar\eta\zeta)^{2j}}{(1+|\eta|^2)^j(1+|\zeta|^2)^j},\quad |\langle n'|n\rangle|^2=\left(\frac{1+n\cdot n'}2\right)^{2j},\quad \frac{2j+1}{\pi}\int\frac{|\zeta\rangle\langle\zeta|}{(1+|\zeta|^2)^2}d^2\zeta=\hat I $$

Audit: Correct. Checked by: exact WL computation for both overlaps (j = 3/2 complete spin matrices); resolution by hand (each diagonal element is (2j+1) C(2j,k) B(k+1, 2j+1-k) = 1).
QF code (exact): none. Finite spin matrices; QF's second-quantization layer does not apply.

<a id="su-1-1-algebra-casimir-discrete-series-ladder-and-holomorphic-realization"></a>
### su(1,1) algebra, Casimir, discrete-series ladder and holomorphic realization
Appears in: P-028; C-045

The su(1,1) generators close with $[\hat K_-,\hat K_+]=2\hat K_0$, the Casimir is $k(k-1)$ in the positive discrete series, and the series has explicit ladder factors and a holomorphic realization.

$$ [\hat K_0,\hat K_\pm]=\pm \hat K_\pm,\ [\hat K_-,\hat K_+]=2\hat K_0,\ \hat C=\hat K_0^2-\tfrac12(\hat K_+\hat K_-+\hat K_-\hat K_+)=k(k-1)\hat I,\ \hat K_+|k,n\rangle=\sqrt{(n+1)(2k+n)}|k,n+1\rangle,\ \hat K_-=\partial_z,\ \hat K_+=z^2\partial_z+2kz,\ \hat K_0=z\partial_z+k $$

Audit: Correct. Checked by: exact WL computation: the differential realization satisfies the commutators and the Casimir on a symbolic $f(z)$, and the ladder factors give $[\hat K_-,\hat K_+]|k,n\rangle=2(k+n)|k,n\rangle$ and Casimir $k(k-1)$ for symbolic $k,n$.
QF code (exact): none. QF has no abstract su(1,1) generators; the oscillator realizations are checked in the single-mode and two-mode realization entries.

**Also stated in C-045 (SU(1,1) algebra and Bergman realization):**

The su(1,1) generators close with [K-, K+] = 2 K0 and Casimir k(k-1), and act on holomorphic functions on the disk as first-order differential operators.

$$ [\hat K_-,\hat K_+]=2\hat K_0,\ [\hat K_0,\hat K_\pm]=\pm\hat K_\pm,\ \hat C=\hat K_0^2-\tfrac12(\hat K_+\hat K_-+\hat K_-\hat K_+)=k(k-1)\hat I;\quad \hat K_-=\partial_z,\ \hat K_+=2kz+z^2\partial_z,\ \hat K_0=k+z\partial_z $$

Audit: Correct. Checked by: exact WL computation (symbolic action on an arbitrary f(z): all four relations vanish).
QF code (exact): none. An abstract differential realization; the bosonic realizations are the QF call of C-046.

<a id="su-1-1-coherent-state-expansion-and-normalization"></a>
### SU(1,1) coherent state expansion and normalization
Appears in: P-030; P-030; C-047; overlap modulus reused in C-097

The su(1,1) coherent state is the normalized exponential of $\zeta\hat K_+$ on the lowest weight, with Pochhammer square-root amplitudes, and it is normalized for $|\zeta|<1$.

$$ |\zeta;k\rangle=(1-|\zeta|^2)^ke^{\zeta \hat K_+}|k,0\rangle=(1-|\zeta|^2)^k\sum_{n\ge0}\sqrt{\frac{\Gamma(n+2k)}{n!\Gamma(2k)}}\zeta^n|k,n\rangle $$

Audit: Correct. Checked by: exact WL computation: the product of ladder factors $\prod_{j<n}\sqrt{(j+1)(2k+j)}/n!$ equals $\sqrt{\Gamma(n+2k)/(n!\Gamma(2k))}$ for symbolic $k,n$, and $\sum_n\Gamma(n+2k)x^n/(n!\Gamma(2k))=(1-x)^{-2k}$ exactly.
QF code (exact): none. No su(1,1) layer in QF.

**Also stated in P-030 (SU(1,1) coherent overlap and annihilator):**

Two su(1,1) coherent states overlap by $(1-\bar\zeta\eta)^{-2k}$ times the normalizations, and each is annihilated by $\hat K_--2\zeta\hat K_0+\zeta^2\hat K_+$.

$$ \langle\zeta;k|\eta;k\rangle=\frac{(1-|\zeta|^2)^k(1-|\eta|^2)^k}{(1-\bar\zeta\eta)^{2k}},\qquad (\hat K_--2\zeta \hat K_0+\zeta^2\hat K_+)|\zeta;k\rangle=0 $$

Audit: Correct. Checked by: exact WL computation: the overlap series sums to $(1-\bar\zeta\eta)^{-2k}$, and the coefficient of $|k,n\rangle$ in the annihilator equation vanishes for symbolic $k,n$ (neighbouring amplitudes differ by $\zeta\sqrt{(n+2k)/(n+1)}$, so the three terms give $\zeta(n+2k)-2\zeta(k+n)+\zeta n=0$).
QF code (exact): none. No su(1,1) layer in QF.

**Also stated in C-047 (SU(1,1) coherent series and overlap):**

The SU(1,1) coherent state is the normalized exponential of the raising operator on the lowest-weight vector; its coefficients are square roots of Gamma-function ratios and its overlap is a power of (1 - conj(eta) zeta).

$$ |\zeta;k\rangle=(1-|\zeta|^2)^k\sum_m\sqrt{\frac{\Gamma(2k+m)}{\Gamma(2k)m!}}\zeta^m|k,m\rangle,\qquad \langle\eta;k|\zeta;k\rangle=\frac{(1-|\eta|^2)^k(1-|\zeta|^2)^k}{(1-\bar\eta\zeta)^{2k}} $$

Conditions: |zeta| < 1; resolution measure (2k-1)/pi needs k > 1/2.
Audit: Correct. For general k the coefficients follow by hand from K+|k,m> = sqrt((m+1)(m+2k)) |k,m+1>. In the odd oscillator sector (k = 3/4, |k,m> = |2m+1>, K+ = a^dag^2/2) QF gives the coefficients for m = 0 to 4, an exact FullSimplify shows the Fock form (1/2)^m sqrt((2m+1)!)/m! equals the Gamma form for every m, and the exact sum of Gamma(2k+m)/(Gamma(2k) m!) x^m is (1-x)^{-2k}, which gives the overlap and the normalization. Checked by: QF exact plus exact WL computation.
QF code (exact):
```wl
{a, ad} = FieldVariables[{1}];
qf = Table[BosonicMatrixElement[{2 m + 1, 1}, Exp[ze ad ** ad/2]], {m, 0, 4}];
atlas = Table[Sqrt[Gamma[3/2 + m]/(Gamma[3/2] m!)] ze^m, {m, 0, 4}];
{qf, Simplify[qf - atlas],
 FullSimplify[(1/2)^m Sqrt[(2 m + 1)!]/m! - Sqrt[Gamma[3/2 + m]/(Gamma[3/2] m!)], m >= 0 && Element[m, Integers]],
 Sum[Gamma[2 k + m]/(Gamma[2 k] m!) x^m, {m, 0, Infinity}]}
```
Output: `{{1, Sqrt[3/2]*ze, (Sqrt[15/2]*ze^2)/2, (Sqrt[35]*ze^3)/4, (3*Sqrt[35/2]*ze^4)/8}, {0, 0, 0, 0, 0}, 0, (1 - x)^(-2*k)}`. Matches: all differences are 0 and the series sums to (1 - x)^(-2 k).

<a id="su-1-1-resolution-of-identity"></a>
### SU(1,1) resolution of identity
Appears in: P-030

su(1,1) coherent projectors integrated with the invariant disk measure resolve the identity only for $k>1/2$.

$$ \frac{2k-1}{\pi}\int_{|\zeta|<1}\frac{d^2\zeta}{(1-|\zeta|^2)^2}|\zeta;k\rangle\langle\zeta;k|=\hat I,\quad k>1/2 $$

Audit: Correct, including the restriction $k>1/2$. Checked by: exact WL computation; the diagonal element reduces to $(2k-1)\Gamma(n+2k)/(n!\Gamma(2k))\int_0^1t^n(1-t)^{2k-2}dt$, which the kernel returns as 1 under the condition $k>1/2$ (the integral diverges otherwise).
QF code (exact): none. No su(1,1) layer in QF.

<a id="root-space-coherent-states-and-displacement-like-generator"></a>
### Root-space coherent states and displacement-like generator
Appears in: P-047

For a general semisimple group, coherent states are exponentials of positive-root generators on a lowest weight, and the unitary displacement uses positive and negative roots.

$$ |\boldsymbol\zeta\rangle=\mathcal N(\boldsymbol\zeta)\exp\left(\sum_{\alpha>0}\zeta_\alpha \hat E_\alpha\right)|-\lambda\rangle,\qquad \hat D(\boldsymbol\xi)=\exp\left(\sum_{\alpha>0}\xi_\alpha \hat E_\alpha-\bar\xi_\alpha \hat E_{-\alpha}\right) $$

Audit: Unchecked: definitions on a general semisimple group; the conversion from $\xi$ to $\zeta$ is not given. Rank-one instances are the SU(2) and SU(1,1) disentangling entries.
QF code (exact): none. QF has no root-space layer.

<a id="compact-kahler-kernel"></a>
### Compact Kahler kernel
Appears in: P2-007

On the sphere the reproducing kernel is $(1+z\bar w)^{1/h}$ with $1/h$ an integer.

$$ K(z,\bar w)=(1+z\bar w)^{1/h},\quad h^{-1}\in\mathbb Z $$

Audit: Correct: with $1/h=2j$ this is the unnormalized SU(2) coherent overlap. Checked by: comparison with the SU(2) overlap check (exact).
QF code (exact): none. No su(2) layer in QF.

<a id="jacobians-of-classical-domain-transformations"></a>
### Jacobians of classical-domain transformations
Appears in: P2-056

The complex Jacobian of the fractional-linear action on a type-I domain is $\det(B^Tz+D)^{-(p+q)}$.

$$ \det(B^Tz+D)^{-(p+q)} $$

Audit: Unchecked: depends on the source's transposed-block convention for type-I domains; not re-derived.
QF code (exact): none. Not an operator identity.

<a id="completeness-criterion-for-discrete-su-1-1-subsystems"></a>
### Completeness criterion for discrete SU(1,1) subsystems
Appears in: P2-004

A discrete SU(1,1) coherent subsystem is complete when the fundamental area is below $\pi/(2k-1)$ and incomplete for large $k$.

$$ S_\Gamma<\frac{\pi}{2k-1}\ \Rightarrow\ \text{complete};\qquad k>m_0^++\tfrac12\ \Rightarrow\ \text{incomplete} $$

Audit: Unchecked: theorem cited from the source (automorphic-form argument), area convention source-specific.
QF code (exact): none. QF has no su(1,1) layer.

<a id="clock-and-shift-matrices-on-the-torus"></a>
### Clock and shift matrices on the torus
Appears in: C-037

On the quantum torus with 2 pi hbar N = 1, translations become N-by-N clock and shift matrices of order N that commute up to a root of unity.

$$ 2\pi\hbar N=1,\quad \hat X^N=\hat Z^N=\hat I,\quad \hat Z\hat X=\omega\hat X\hat Z,\ \omega=e^{2\pi i/N} $$

Audit: Correct for Z|j> = omega^j |j>, X|j> = |j+1>. Checked by: exact WL computation (N = 5; the torus space is finite, so these are the complete operators). The book's (6.18)-(6.19) matrices were not seen; prior audit only for the claim that they fix this orientation.
QF code (exact): none. Finite clock and shift matrices; QF's second-quantization layer does not apply.

<a id="gauss-decomposition-of-sl-2-c"></a>
### Gauss decomposition of SL(2,C)
Appears in: C-040

A unimodular 2-by-2 matrix with nonzero lower-right entry factors uniquely into upper-triangular, diagonal and lower-triangular pieces.

$$ g=\begin{pmatrix}1&\beta/\delta\\0&1\end{pmatrix}\begin{pmatrix}\delta^{-1}&0\\0&\delta\end{pmatrix}\begin{pmatrix}1&0\\\gamma/\delta&1\end{pmatrix} $$

Conditions: alpha delta - beta gamma = 1, delta != 0.
Audit: Correct. Checked by: exact WL computation (symbolic product returns g).
QF code (exact): none. A 2-by-2 matrix identity; QF's second-quantization layer does not apply.

<a id="spin-to-oscillator-contraction"></a>
### Spin-to-oscillator contraction
Appears in: C-044

As j grows, the spin coherent kernel at scale zeta/sqrt(2j) tends to the bosonic coherent kernel.

$$ \text{spin kernel at scale }\zeta/\sqrt{2j}\ \to\ e^{\bar\eta\zeta-|\eta|^2/2}\quad (j\to\infty) $$

Audit: Correct as a limit ((1 + conj(eta) zeta/(2j))^{2j} -> e^{conj(eta) zeta}). Checked by: hand derivation.
QF code (exact): none. A large-spin limit; QF's second-quantization layer does not apply.

<a id="singular-oscillator-lowest-weight"></a>
### Singular oscillator lowest weight
Appears in: C-052

For the oscillator with an inverse-square potential on the half line, the lowest-weight state is x^s e^{-x^2/2} with s fixed by g, and its Bargmann index is k = (2s+1)/4.

$$ \psi\propto x^se^{-x^2/2},\quad s=\tfrac12+\sqrt{\tfrac14+g^2},\quad k=\frac{2s+1}{4} $$

Conditions: H = (-d^2 + x^2 + g^2/x^2)/2 on L^2(R+), chosen self-adjoint realization.
Audit: Correct (energy (2s+1)/2 = 2k with K0 = H/2). Checked by: exact WL computation (symbolic eigenvalue check).
QF code (exact): none. A differential operator on L^2(R+), outside QF's Fock-space layer.

<a id="berezin-symbol-product-on-the-disk"></a>
### Berezin symbol product on the disk
Appears in: C-097

On the disk, the covariant symbol of a product is an integral of the two off-diagonal symbols against the squared overlap of SU(1,1) coherent states.

$$ (\hat A\hat B)_c(z,\bar z)=\int_DA_c(z,\bar w)B_c(w,\bar z)|\langle\psi_z,\psi_w\rangle|^2d\nu_n(w),\qquad |\langle\psi_z,\psi_w\rangle|^2=\left[\frac{(1-|z|^2)(1-|w|^2)}{|1-\bar zw|^2}\right]^n $$

Audit: Correct (insert the resolution of identity; overlap is C-047 with n = 2k). Checked by: hand derivation.
QF code (exact): none. QF has no SU(1,1) coherent-state symbol calculus.

<a id="non-formula-records"></a>
## Non-formula records

These atlas records state domains, theorems in prose, or coverage notes and carry no formula to compute.

- X-008: a five-step derivation route (commutator closure, ordered ansatz, Riccati ODEs, initial conditions and charts, limits for cubic and quartic generators) plus a statement of what was checked when the atlas was assembled. No formula; the claims are consistent with the checks above (the 2x2 compact and noncompact factorizations and the two-mode CAR realization all match).
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
- C-005: analytic core for complex labels (domain statement).
- C-033: quadratic replacement gives O(hbar^{1/2}) and Ehrenfest-time errors; asymptotic, Unchecked.
- C-034: branch prescription for Gaussian determinants (Maslov index), prose.
- C-062: Bogoliubov fluctuation dynamics around a Hartree solution; convergence statement, Unchecked.
- C-063: Dyson/Gevrey/Borel summability statement, Unchecked.
- C-081: Spin double cover and O(2n) versus SO(2n) implementers, prose (consistent with C-087 and with c <-> c^dag being odd).
- C-092: stationary phase is asymptotic, prose.
- C-093: Moyal remainder and functional calculus, prose.
- C-095: Schur uniqueness and square integrability for group coherent states, prose.
- C-096: Cartan polar decomposition, prose.
- G-085: conditional measurement prepares cats; the only formula ($\propto\hat I\pm\hat\Pi$) is covered by the cat parity entry. Listed in "Coherent-state cat" Appears in.

<a id="index-of-atlas-ids"></a>
## Index of atlas IDs

- X-001: [Displacement: normal and antinormal factorization](#displacement-normal-and-antinormal-factorization); [Displacement: adjoint action on a](#displacement-adjoint-action-on-a); [Displacement: composition law](#displacement-composition-law)
- X-002: [Single-mode squeeze: normal-ordered factorization](#single-mode-squeeze-normal-ordered-factorization); [Single-mode squeeze: adjoint action on a](#single-mode-squeeze-adjoint-action-on-a); [Single-mode squeezed vacuum](#single-mode-squeezed-vacuum)
- X-003: [Two-mode squeeze: normal-ordered factorization](#two-mode-squeeze-normal-ordered-factorization); [Two-mode squeezed vacuum](#two-mode-squeezed-vacuum)
- X-004: [SU(2) disentangling (compact case)](#su-2-disentangling-compact-case); [SU(1,1) disentangling (noncompact case)](#su-1-1-disentangling-noncompact-case); [Fermion-pair realization of su(2)](#fermion-pair-realization-of-su-2); [Fermion-pair rotation of the vacuum](#fermion-pair-rotation-of-the-vacuum)
- X-005: [Hadamard series for a closed linear commutator action](#hadamard-series-for-a-closed-linear-commutator-action)
- X-006: [Bosonic multimode Gaussian vacuum](#bosonic-multimode-gaussian-vacuum); [Fermionic multimode Gaussian vacuum](#fermionic-multimode-gaussian-vacuum)
- X-007: [Number-conserving exponential in normal order](#number-conserving-exponential-in-normal-order); [Vacuum projector as a normal-ordered exponential](#vacuum-projector-as-a-normal-ordered-exponential); [Traces of second-quantized one-particle operators](#traces-of-second-quantized-one-particle-operators)
- X-008: [Non-formula records](#non-formula-records)
- X-009: [Multimode canonical action and its constraints](#multimode-canonical-action-and-its-constraints); [Multimode creator-number-annihilator factorization](#multimode-creator-number-annihilator-factorization); [Bosonic coherent-state kernel of the factorized implementer](#bosonic-coherent-state-kernel-of-the-factorized-implementer)
- B-001: [Non-formula records](#non-formula-records)
- B-002: [Non-formula records](#non-formula-records)
- B-003: [Smeared CCR and CAR](#smeared-ccr-and-car)
- B-004: [Fock basis and ladder action](#fock-basis-and-ladder-action)
- B-005: [Norm bound for fermionic fields](#norm-bound-for-fermionic-fields)
- B-006: [Bosonic quadratures](#bosonic-quadratures)
- B-007: [Non-formula records](#non-formula-records)
- B-008: [Non-formula records](#non-formula-records)
- B-009: [Non-formula records](#non-formula-records)
- B-010: [Non-formula records](#non-formula-records)
- B-011: [Non-formula records](#non-formula-records)
- B-012: [Non-formula records](#non-formula-records)
- B-013: [Number-conserving exponential in normal order](#number-conserving-exponential-in-normal-order); [Vacuum projector as a normal-ordered exponential](#vacuum-projector-as-a-normal-ordered-exponential)
- B-014: [Matrix kernel versus normal symbol](#matrix-kernel-versus-normal-symbol)
- B-015: [Differential realization on normal symbols](#differential-realization-on-normal-symbols)
- B-016: [Fermionic symbols: which side the derivative acts on](#fermionic-symbols-which-side-the-derivative-acts-on)
- B-017: [Poisson (coherent) vectors](#poisson-coherent-vectors)
- B-018: [Non-formula records](#non-formula-records)
- B-019: [Symbols as coherent matrix elements](#symbols-as-coherent-matrix-elements)
- B-020: [Rank-one projector symbols](#rank-one-projector-symbols)
- B-021: [Non-formula records](#non-formula-records)
- B-022: [Real Gaussian integral](#real-gaussian-integral)
- B-023: [Doubled complex Gaussian with sources](#doubled-complex-gaussian-with-sources)
- B-024: [Non-formula records](#non-formula-records)
- B-025: [Non-formula records](#non-formula-records)
- B-026: [Kernel acting on a vector](#kernel-acting-on-a-vector)
- B-027: [Product of operators: matrix kernels](#product-of-operators-matrix-kernels); [Product of operators: normal symbols](#product-of-operators-normal-symbols)
- B-028: [Number-conserving exponential in normal order](#number-conserving-exponential-in-normal-order); [Bosonic trace](#bosonic-trace)
- B-029: [Grassmann algebra relations](#grassmann-algebra-relations)
- B-030: [Left versus right derivatives](#left-versus-right-derivatives)
- B-031: [Berezin integration rules](#berezin-integration-rules)
- B-032: [Grassmann change of variables](#grassmann-change-of-variables)
- B-033: [Fermionic Gaussian is a Pfaffian](#fermionic-gaussian-is-a-pfaffian)
- B-034: [Clifford operators from multiplication and derivative](#clifford-operators-from-multiplication-and-derivative)
- B-035: [Non-formula records](#non-formula-records)
- B-036: [Doubled fermionic Gaussian (corrected determinant)](#doubled-fermionic-gaussian-corrected-determinant)
- B-037: [Non-formula records](#non-formula-records)
- B-038: [Non-formula records](#non-formula-records)
- B-039: [Non-formula records](#non-formula-records)
- B-040: [Number-conserving exponential in normal order](#number-conserving-exponential-in-normal-order); [Fermionic trace](#fermionic-trace)
- B-041: [Normalized fermionic trace](#normalized-fermionic-trace)
- B-042: [Bosonic affine canonical transformation](#bosonic-affine-canonical-transformation)
- B-043: [Bosonic block relations and inverse](#bosonic-block-relations-and-inverse)
- B-044: [Symplectic form and inverse in quadratures](#symplectic-form-and-inverse-in-quadratures)
- B-045: [Squeezing contraction, both charts](#squeezing-contraction-both-charts)
- B-046: [Non-formula records](#non-formula-records)
- B-047: [Non-formula records](#non-formula-records)
- B-048: [Bosonic multimode Gaussian vacuum](#bosonic-multimode-gaussian-vacuum); [Transformed vacuum: affine normalization (4.25)](#transformed-vacuum-affine-normalization-4-25)
- B-049: [Multimode canonical action and its constraints](#multimode-canonical-action-and-its-constraints); [Transformed vacuum: affine normalization (4.25)](#transformed-vacuum-affine-normalization-4-25); [Full bosonic Gaussian kernel](#full-bosonic-gaussian-kernel)
- B-050: [Displacement: normal and antinormal factorization](#displacement-normal-and-antinormal-factorization); [Displacement in the source direction](#displacement-in-the-source-direction)
- B-051: [Single-mode squeeze: normal-ordered factorization](#single-mode-squeeze-normal-ordered-factorization); [One-mode squeezing implementer and kernel](#one-mode-squeezing-implementer-and-kernel)
- B-052: [Non-formula records](#non-formula-records)
- B-053: [Passive transformation of normal symbols](#passive-transformation-of-normal-symbols)
- B-054: [General bosonic transformation kernel for normal symbols](#general-bosonic-transformation-kernel-for-normal-symbols)
- B-055: [Gaussian normal symbols stay Gaussian under conjugation](#gaussian-normal-symbols-stay-gaussian-under-conjugation)
- B-056: [Fermionic Bogoliubov constraints](#fermionic-bogoliubov-constraints)
- B-057: [Fermionic Bogoliubov constraints](#fermionic-bogoliubov-constraints); [Non-formula records](#non-formula-records)
- B-058: [Non-formula records](#non-formula-records)
- B-059: [Non-formula records](#non-formula-records)
- B-060: [Multimode canonical action and its constraints](#multimode-canonical-action-and-its-constraints); [Fermionic Gaussian implementer kernel](#fermionic-gaussian-implementer-kernel)
- B-061: [Bosonic multimode Gaussian vacuum](#bosonic-multimode-gaussian-vacuum); [Fermionic multimode Gaussian vacuum](#fermionic-multimode-gaussian-vacuum); [Multimode canonical action and its constraints](#multimode-canonical-action-and-its-constraints)
- B-062: [Multimode canonical action and its constraints](#multimode-canonical-action-and-its-constraints); [Non-formula records](#non-formula-records)
- B-063: [Fermionic transformation of normal symbols with mixing](#fermionic-transformation-of-normal-symbols-with-mixing)
- B-064: [Non-formula records](#non-formula-records)
- B-065: [General quadratic expression](#general-quadratic-expression)
- B-066: [Non-formula records](#non-formula-records)
- B-067: [Non-formula records](#non-formula-records)
- B-068: [Bound for exponentials of bounded perturbations](#bound-for-exponentials-of-bounded-perturbations)
- B-069: [Classical generator and Heisenberg action of a quadratic Hamiltonian](#classical-generator-and-heisenberg-action-of-a-quadratic-hamiltonian)
- B-070: [Non-formula records](#non-formula-records)
- B-071: [Hadamard series for a closed linear commutator action](#hadamard-series-for-a-closed-linear-commutator-action); [Full bosonic Gaussian kernel](#full-bosonic-gaussian-kernel); [Classical generator and Heisenberg action of a quadratic Hamiltonian](#classical-generator-and-heisenberg-action-of-a-quadratic-hamiltonian); [Quadratic propagator vacuum coefficient (corrected (6.23))](#quadratic-propagator-vacuum-coefficient-corrected-6-23)
- B-072: [Classical generator and Heisenberg action of a quadratic Hamiltonian](#classical-generator-and-heisenberg-action-of-a-quadratic-hamiltonian)
- B-073: [Logarithmic derivative of a determinant](#logarithmic-derivative-of-a-determinant)
- B-074: [Hadamard series for a closed linear commutator action](#hadamard-series-for-a-closed-linear-commutator-action); [Fermionic quadratic propagator](#fermionic-quadratic-propagator)
- B-075: [Resolvent kernel of a number-type Hamiltonian](#resolvent-kernel-of-a-number-type-hamiltonian)
- B-076: [Non-formula records](#non-formula-records)
- B-077: [Quadratic propagator vacuum coefficient (corrected (6.23))](#quadratic-propagator-vacuum-coefficient-corrected-6-23); [Non-formula records](#non-formula-records)
- B-078: [Weighted sufficient conditions for a quadratic normal form](#weighted-sufficient-conditions-for-a-quadratic-normal-form)
- B-079: [Non-formula records](#non-formula-records)
- B-080: [Non-formula records](#non-formula-records)
- B-081: [Forced oscillator: completing the square](#forced-oscillator-completing-the-square)
- B-082: [Non-formula records](#non-formula-records)
- B-083: [Regularized forced Hamiltonian](#regularized-forced-hamiltonian)
- B-084: [Bosonic canonical diagonalization](#bosonic-canonical-diagonalization)
- B-085: [Riccati fixed point for the diagonalizing transformation](#riccati-fixed-point-for-the-diagonalizing-transformation)
- B-086: [Fermionic diagonalization conditions](#fermionic-diagonalization-conditions)
- B-087: [Square-root perturbation stays in the ideal](#square-root-perturbation-stays-in-the-ideal)
- B-088: [Non-formula records](#non-formula-records)
- B-089: [Non-formula records](#non-formula-records)
- B-090: [Non-formula records](#non-formula-records)
- B-091: [Non-formula records](#non-formula-records)
- B-092: [Scattering unitary in the Thirring model](#scattering-unitary-in-the-thirring-model)
- B-093: [Pair potential on the n-particle sector](#pair-potential-on-the-n-particle-sector)
- B-094: [Dyson series with the e^{+itH} convention](#dyson-series-with-the-e-ith-convention)
- B-095: [Free time dependence and contraction kernel](#free-time-dependence-and-contraction-kernel)
- B-096: [Non-formula records](#non-formula-records)
- B-097: [Retarded contraction and the two-field time-ordered product](#retarded-contraction-and-the-two-field-time-ordered-product)
- B-098: [Non-formula records](#non-formula-records)
- B-099: [Wick's theorem as an exponential of the contraction operator](#wick-s-theorem-as-an-exponential-of-the-contraction-operator)
- B-100: [Normal symbol of interaction-picture evolution](#normal-symbol-of-interaction-picture-evolution)
- B-101: [Contractions acting on an exponential source](#contractions-acting-on-an-exponential-source)
- B-102: [Determinant identity for a quadratic symbol](#determinant-identity-for-a-quadratic-symbol)
- B-103: [Auxiliary-field decoupling of a quartic term](#auxiliary-field-decoupling-of-a-quartic-term)
- B-104: [Retarded loops vanish](#retarded-loops-vanish)
- B-105: [Number-conserving exponential in normal order](#number-conserving-exponential-in-normal-order); [Auxiliary-field representation of the propagator](#auxiliary-field-representation-of-the-propagator)
- B-106: [Grand partition functions](#grand-partition-functions)
- B-107: [Number-conserving exponential in normal order](#number-conserving-exponential-in-normal-order); [Vacuum projector as a normal-ordered exponential](#vacuum-projector-as-a-normal-ordered-exponential); [Traces of second-quantized one-particle operators](#traces-of-second-quantized-one-particle-operators)
- B-108: [Non-formula records](#non-formula-records)
- P-001: [Canonical commutators, ladder operators and the multimode displacement](#canonical-commutators-ladder-operators-and-the-multimode-displacement)
- P-002: [Displacement: normal and antinormal factorization](#displacement-normal-and-antinormal-factorization); [Central-commutator BCH](#central-commutator-bch)
- P-003: [Displacement: normal and antinormal factorization](#displacement-normal-and-antinormal-factorization); [Displacement: composition law](#displacement-composition-law); [Displacement commutation phase, group commutator and polygon phase](#displacement-commutation-phase-group-commutator-and-polygon-phase)
- P-004: [Non-formula records](#non-formula-records)
- P-005: [Fock basis and ladder action](#fock-basis-and-ladder-action)
- P-006: [Displacement in position representation](#displacement-in-position-representation)
- P-007: [Overlap of displaced fiducial states](#overlap-of-displaced-fiducial-states)
- P-008: [Displacement: adjoint action on a](#displacement-adjoint-action-on-a)
- P-009: [Bogoliubov combination keeps the canonical commutator](#bogoliubov-combination-keeps-the-canonical-commutator)
- P-010: [Displacement: normal and antinormal factorization](#displacement-normal-and-antinormal-factorization); [Central-commutator BCH](#central-commutator-bch)
- P-011: [Differential correspondences for multiplying a displacement](#differential-correspondences-for-multiplying-a-displacement)
- P-012: [Coherent state expansion, eigenvalue equation, overlap and resolution](#coherent-state-expansion-eigenvalue-equation-overlap-and-resolution)
- P-013: [Bargmann realization and displacement as weighted translation](#bargmann-realization-and-displacement-as-weighted-translation); [Coordinate to Bargmann kernel](#coordinate-to-bargmann-kernel)
- P-014: [Lattice commutation condition and critical coherent lattice](#lattice-commutation-condition-and-critical-coherent-lattice)
- P-015: [Covariant (Q) and contravariant (P) symbols from ordered expansions](#covariant-q-and-contravariant-p-symbols-from-ordered-expansions); [Gaussian smoothing kernels between P, W and Q (derived normalization correction)](#gaussian-smoothing-kernels-between-p-w-and-q-derived-normalization-correction)
- P-016: [Coherent-symbol bounds on exponential traces (Berezin-Lieb, Golden-Thompson)](#coherent-symbol-bounds-on-exponential-traces-berezin-lieb-golden-thompson)
- P-017: [Q and P symbols of a displacement](#q-and-p-symbols-of-a-displacement); [Distributional trace of displacement and operator reconstruction](#distributional-trace-of-displacement-and-operator-reconstruction)
- P-018: [Displacement matrix elements in associated Laguerre form](#displacement-matrix-elements-in-associated-laguerre-form); [Coherent generating kernel of the displacement](#coherent-generating-kernel-of-the-displacement)
- P-019: [Weyl operator from displaced parity, Weyl symbol and reconstruction](#weyl-operator-from-displaced-parity-weyl-symbol-and-reconstruction); [Gaussian smoothing kernels between P, W and Q (derived normalization correction)](#gaussian-smoothing-kernels-between-p-w-and-q-derived-normalization-correction)
- P-020: [Distributional trace of displacement and operator reconstruction](#distributional-trace-of-displacement-and-operator-reconstruction); [Ordered characteristic functions differ by Gaussian factors](#ordered-characteristic-functions-differ-by-gaussian-factors); [Normally ordered moments from the normal characteristic function](#normally-ordered-moments-from-the-normal-characteristic-function); [Weyl characteristic function of a number-diagonal state](#weyl-characteristic-function-of-a-number-diagonal-state)
- P-021: [Non-formula records](#non-formula-records)
- P-022: [Non-formula records](#non-formula-records)
- P-023: [Non-formula records](#non-formula-records)
- P-024: [su(2) commutators, ladder action and holomorphic realization](#su-2-commutators-ladder-action-and-holomorphic-realization)
- P-025: [SU(2) disentangling (compact case)](#su-2-disentangling-compact-case)
- P-026: [SU(2) coherent state expansion](#su-2-coherent-state-expansion); [SU(2) coherent overlap and annihilator](#su-2-coherent-overlap-and-annihilator)
- P-027: [Non-formula records](#non-formula-records)
- P-028: [Single-mode squeeze: normal-ordered factorization](#single-mode-squeeze-normal-ordered-factorization); [su(1,1) algebra, Casimir, discrete-series ladder and holomorphic realization](#su-1-1-algebra-casimir-discrete-series-ladder-and-holomorphic-realization)
- P-029: [Two-mode squeeze: normal-ordered factorization](#two-mode-squeeze-normal-ordered-factorization); [SU(2) disentangling (compact case)](#su-2-disentangling-compact-case); [SU(1,1) disentangling (noncompact case)](#su-1-1-disentangling-noncompact-case)
- P-030: [SU(1,1) coherent state expansion and normalization](#su-1-1-coherent-state-expansion-and-normalization); [SU(1,1) resolution of identity](#su-1-1-resolution-of-identity)
- P-031: [Single-mode squeeze: normal-ordered factorization](#single-mode-squeeze-normal-ordered-factorization); [Single-mode su(1,1) realization and Casimir](#single-mode-su-1-1-realization-and-casimir)
- P-032: [Two-mode squeeze: normal-ordered factorization](#two-mode-squeeze-normal-ordered-factorization); [Two-mode su(1,1) realization and Casimir](#two-mode-su-1-1-realization-and-casimir)
- P-033: [Non-formula records](#non-formula-records)
- P-034: [Hadamard series for a closed linear commutator action](#hadamard-series-for-a-closed-linear-commutator-action); [Multimode canonical action and its constraints](#multimode-canonical-action-and-its-constraints)
- P-035: [Hadamard series for a closed linear commutator action](#hadamard-series-for-a-closed-linear-commutator-action); [Quadratic bosonic generators as sp(2N,R) matrices](#quadratic-bosonic-generators-as-sp-2n-r-matrices); [Bosonic quadratic commutator table (E, F, G)](#bosonic-quadratic-commutator-table-e-f-g)
- P-036: [Bosonic multimode Gaussian vacuum](#bosonic-multimode-gaussian-vacuum); [Multimode canonical action and its constraints](#multimode-canonical-action-and-its-constraints)
- P-037: [Bosonic Gaussian normalization and overlap](#bosonic-gaussian-normalization-and-overlap)
- P-038: [Fractional-linear action on symmetric Gaussian parameters](#fractional-linear-action-on-symmetric-gaussian-parameters)
- P-039: [Bosonic Gaussian bilinear symbols and purity constraints](#bosonic-gaussian-bilinear-symbols-and-purity-constraints)
- P-040: [Non-formula records](#non-formula-records)
- P-041: [SU(2) disentangling (compact case)](#su-2-disentangling-compact-case); [Hadamard series for a closed linear commutator action](#hadamard-series-for-a-closed-linear-commutator-action); [Fermionic Bogoliubov constraints](#fermionic-bogoliubov-constraints)
- P-042: [Hadamard series for a closed linear commutator action](#hadamard-series-for-a-closed-linear-commutator-action); [Multimode canonical action and its constraints](#multimode-canonical-action-and-its-constraints); [Fermionic quadratic commutator table](#fermionic-quadratic-commutator-table)
- P-043: [SU(2) disentangling (compact case)](#su-2-disentangling-compact-case); [Multimode canonical action and its constraints](#multimode-canonical-action-and-its-constraints); [Fermionic pair-exponential state, normalization and overlap](#fermionic-pair-exponential-state-normalization-and-overlap)
- P-044: [Fermionic analytic group action](#fermionic-analytic-group-action)
- P-045: [Bosonic multimode Gaussian vacuum](#bosonic-multimode-gaussian-vacuum); [Fermionic bilinear symbols and purity constraints](#fermionic-bilinear-symbols-and-purity-constraints)
- P-046: [Non-formula records](#non-formula-records)
- P-047: [Root-space coherent states and displacement-like generator](#root-space-coherent-states-and-displacement-like-generator)
- P-048: [Determinant and log-determinant series](#determinant-and-log-determinant-series)
- P-049: [Non-formula records](#non-formula-records)
- P2-001: [Non-formula records](#non-formula-records)
- P2-002: [Siegel upper half-space and symplectic fractional-linear action](#siegel-upper-half-space-and-symplectic-fractional-linear-action)
- P2-003: [Non-formula records](#non-formula-records)
- P2-004: [Completeness criterion for discrete SU(1,1) subsystems](#completeness-criterion-for-discrete-su-1-1-subsystems)
- P2-005: [Off-diagonal covariant symbol by analytic continuation](#off-diagonal-covariant-symbol-by-analytic-continuation)
- P2-006: [Disk star product to first order in h](#disk-star-product-to-first-order-in-h)
- P2-007: [Compact Kahler kernel](#compact-kahler-kernel)
- P2-008: [Non-formula records](#non-formula-records)
- P2-009: [Non-formula records](#non-formula-records)
- P2-010: [Forced oscillator: displacement plus accumulated phase](#forced-oscillator-displacement-plus-accumulated-phase)
- P2-011: [Displacement matrix elements in associated Laguerre form](#displacement-matrix-elements-in-associated-laguerre-form)
- P2-012: [Single-mode su(1,1) realization and Casimir](#single-mode-su-1-1-realization-and-casimir); [Parametric oscillator su(1,1) generators](#parametric-oscillator-su-1-1-generators)
- P2-013: [Riccati equation and phase for squeezing dynamics (phase sign corrected)](#riccati-equation-and-phase-for-squeezing-dynamics-phase-sign-corrected)
- P2-014: [Floquet classification of a quadratic oscillator](#floquet-classification-of-a-quadratic-oscillator)
- P2-015: [Singular oscillator ground state and su(1,1) label](#singular-oscillator-ground-state-and-su-1-1-label)
- P2-016: [Calogero degeneracy generating function](#calogero-degeneracy-generating-function)
- P2-017: [Affine Gaussian state normalization (derived correction)](#affine-gaussian-state-normalization-derived-correction)
- P2-018: [Forced parametric oscillator generating kernel and vacuum amplitude](#forced-parametric-oscillator-generating-kernel-and-vacuum-amplitude)
- P2-019: [Two-mode squeeze: normal-ordered factorization](#two-mode-squeeze-normal-ordered-factorization); [Two-mode su(1,1) realization and Casimir](#two-mode-su-1-1-realization-and-casimir)
- P2-020: [Multimode particle-antiparticle bosonic blocks](#multimode-particle-antiparticle-bosonic-blocks)
- P2-021: [Multimode paired vacuum normalization](#multimode-paired-vacuum-normalization)
- P2-022: [Multimode paired vacuum overlap](#multimode-paired-vacuum-overlap)
- P2-023: [Bosonic pair-number generating determinant (trace-log sign corrected)](#bosonic-pair-number-generating-determinant-trace-log-sign-corrected)
- P2-024: [Multimode pair Riccati equation and vacuum phase](#multimode-pair-riccati-equation-and-vacuum-phase)
- P2-025: [SU(2) disentangling (compact case)](#su-2-disentangling-compact-case); [Fermion-pair realization of su(2)](#fermion-pair-realization-of-su-2)
- P2-026: [Multimode fermionic blocks](#multimode-fermionic-blocks)
- P2-027: [Fermionic paired vacuum normalization](#fermionic-paired-vacuum-normalization)
- P2-028: [Fermionic paired vacuum overlap](#fermionic-paired-vacuum-overlap)
- P2-029: [Fermionic pair-number generating determinant](#fermionic-pair-number-generating-determinant)
- P2-030: [Non-formula records](#non-formula-records)
- P2-031: [Non-formula records](#non-formula-records)
- P2-032: [Quasiclassical partition-function bounds](#quasiclassical-partition-function-bounds)
- P2-033: [Non-formula records](#non-formula-records)
- P2-034: [Thermal oscillator master equation](#thermal-oscillator-master-equation)
- P2-035: [Central-commutator BCH](#central-commutator-bch); [Ordered characteristic functions differ by Gaussian factors](#ordered-characteristic-functions-differ-by-gaussian-factors)
- P2-036: [Exact attenuation solution for ordered characteristic functions](#exact-attenuation-solution-for-ordered-characteristic-functions)
- P2-037: [Fokker-Planck drift, diffusion and covariance evolution](#fokker-planck-drift-diffusion-and-covariance-evolution)
- P2-038: [Gaussian convolution kernel for attenuation (sign corrected)](#gaussian-convolution-kernel-for-attenuation-sign-corrected)
- P2-039: [Exact attenuation solution for ordered characteristic functions](#exact-attenuation-solution-for-ordered-characteristic-functions)
- P2-040: [Analytic density kernel](#analytic-density-kernel)
- P2-041: [Non-formula records](#non-formula-records)
- P2-042: [Drive equation and displacement channel factor](#drive-equation-and-displacement-channel-factor)
- P2-043: [Displaced thermal density kernel](#displaced-thermal-density-kernel)
- P2-044: [Populations of a displaced thermal state](#populations-of-a-displaced-thermal-state)
- P2-045: [Non-formula records](#non-formula-records)
- P2-046: [Free oscillator propagator and thermal trace in coherent states](#free-oscillator-propagator-and-thermal-trace-in-coherent-states)
- P2-047: [Non-formula records](#non-formula-records)
- P2-048: [Stable oscillator proper-time trace](#stable-oscillator-proper-time-trace)
- P2-049: [Non-formula records](#non-formula-records)
- P2-050: [Lattice commutation condition and critical coherent lattice](#lattice-commutation-condition-and-critical-coherent-lattice); [Non-formula records](#non-formula-records)
- P2-051: [Coherent generating kernel of the displacement](#coherent-generating-kernel-of-the-displacement); [Displacement matrix elements as a finite sum](#displacement-matrix-elements-as-a-finite-sum)
- P2-052: [Displacement matrix elements in associated Laguerre form](#displacement-matrix-elements-in-associated-laguerre-form)
- P2-053: [Displacement: normal and antinormal factorization](#displacement-normal-and-antinormal-factorization); [Laguerre generating function in the lower index](#laguerre-generating-function-in-the-lower-index); [Laguerre generating function along the diagonal m-n](#laguerre-generating-function-along-the-diagonal-m-n)
- P2-054: [Generating function of displacement transition probabilities (prefactor corrected)](#generating-function-of-displacement-transition-probabilities-prefactor-corrected)
- P2-055: [Displacement matrix elements as Wirtinger derivatives](#displacement-matrix-elements-as-wirtinger-derivatives)
- P2-056: [Jacobians of classical-domain transformations](#jacobians-of-classical-domain-transformations)
- P2-057: [Non-formula records](#non-formula-records)
- P2-058: [SU(1,1) disentangling (noncompact case)](#su-1-1-disentangling-noncompact-case); [Non-formula records](#non-formula-records)
- C-001: [Canonical commutator of position and momentum](#canonical-commutator-of-position-and-momentum)
- C-002: [Central-commutator BCH](#central-commutator-bch)
- C-003: [Displacement: composition law](#displacement-composition-law); [Infinite-mode displacement, exponential vector, product, overlap](#infinite-mode-displacement-exponential-vector-product-overlap)
- C-004: [Displacement conjugation and coordinate factorization](#displacement-conjugation-and-coordinate-factorization)
- C-005: [Non-formula records](#non-formula-records)
- C-006: [Displacement: normal and antinormal factorization](#displacement-normal-and-antinormal-factorization); [Central-commutator BCH](#central-commutator-bch); [Infinite-mode displacement, exponential vector, product, overlap](#infinite-mode-displacement-exponential-vector-product-overlap)
- C-007: [Coherent state expansion, eigenvalue equation, overlap and resolution](#coherent-state-expansion-eigenvalue-equation-overlap-and-resolution); [Infinite-mode displacement, exponential vector, product, overlap](#infinite-mode-displacement-exponential-vector-product-overlap)
- C-008: [Vacuum displacement amplitude and coherent overlap](#vacuum-displacement-amplitude-and-coherent-overlap); [Infinite-mode displacement, exponential vector, product, overlap](#infinite-mode-displacement-exponential-vector-product-overlap)
- C-009: [Harmonic rotation of displacements and coherent states](#harmonic-rotation-of-displacements-and-coherent-states)
- C-010: [Vacuum displacement amplitude and coherent overlap](#vacuum-displacement-amplitude-and-coherent-overlap)
- C-011: [Bargmann realization and displacement as weighted translation](#bargmann-realization-and-displacement-as-weighted-translation)
- C-012: [Mehler kernel and Euclidean oscillator kernel](#mehler-kernel-and-euclidean-oscillator-kernel)
- C-013: [Weyl quantization kernel convention](#weyl-quantization-kernel-convention)
- C-014: [Symbol traces and Wigner overlaps](#symbol-traces-and-wigner-overlaps)
- C-015: [Symplectic covariance of Weyl symbols](#symplectic-covariance-of-weyl-symbols)
- C-016: [Moyal product and exact quadratic commutator](#moyal-product-and-exact-quadratic-commutator)
- C-017: [Gaussian smoothing kernels between P, W and Q (derived normalization correction)](#gaussian-smoothing-kernels-between-p-w-and-q-derived-normalization-correction)
- C-018: [Hadamard series for a closed linear commutator action](#hadamard-series-for-a-closed-linear-commutator-action); [Linear quadratic flow and Heisenberg action](#linear-quadratic-flow-and-heisenberg-action)
- C-019: [Gaussian propagation by a linear-fractional map](#gaussian-propagation-by-a-linear-fractional-map)
- C-020: [Siegel-domain invariance and determinant cocycle](#siegel-domain-invariance-and-determinant-cocycle)
- C-021: [Matrix Riccati and amplitude equations](#matrix-riccati-and-amplitude-equations)
- C-022: [Coherent matrix elements of a metaplectic operator](#coherent-matrix-elements-of-a-metaplectic-operator)
- C-023: [Linear quadratic flow and Heisenberg action](#linear-quadratic-flow-and-heisenberg-action); [Metaplectic generator from a symplectic path](#metaplectic-generator-from-a-symplectic-path)
- C-024: [Hadamard series for a closed linear commutator action](#hadamard-series-for-a-closed-linear-commutator-action); [Bogoliubov blocks of a symplectic matrix](#bogoliubov-blocks-of-a-symplectic-matrix)
- C-025: [Multimode squeeze operator and bounded squeeze coordinate](#multimode-squeeze-operator-and-bounded-squeeze-coordinate)
- C-026: [Squeeze conjugation in block form](#squeeze-conjugation-in-block-form)
- C-027: [Squeezed vacuum from the Bargmann Riccati equation](#squeezed-vacuum-from-the-bargmann-riccati-equation)
- C-028: [Coordinate form of the squeezed vacuum](#coordinate-form-of-the-squeezed-vacuum)
- C-029: [Squeeze-rotation factorization and its vacuum phase](#squeeze-rotation-factorization-and-its-vacuum-phase)
- C-030: [Linear quadratic flow and Heisenberg action](#linear-quadratic-flow-and-heisenberg-action); [Weyl symbol of a metaplectic operator (Cayley form)](#weyl-symbol-of-a-metaplectic-operator-cayley-form); [Weyl symbol of the oscillator exponential](#weyl-symbol-of-the-oscillator-exponential)
- C-031: [Translation coefficient of a metaplectic operator](#translation-coefficient-of-a-metaplectic-operator)
- C-032: [Propagator of a driven quadratic Hamiltonian](#propagator-of-a-driven-quadratic-hamiltonian)
- C-033: [Non-formula records](#non-formula-records)
- C-034: [Non-formula records](#non-formula-records)
- C-035: [Periodic quadratic evolution and quasi-energies](#periodic-quadratic-evolution-and-quasi-energies)
- C-036: [Periodic-orbit determinant in regularized traces](#periodic-orbit-determinant-in-regularized-traces)
- C-037: [Clock and shift matrices on the torus](#clock-and-shift-matrices-on-the-torus)
- C-038: [su(2) commutators, ladder action and holomorphic realization](#su-2-commutators-ladder-action-and-holomorphic-realization)
- C-039: [Exponential of an su(2) matrix](#exponential-of-an-su-2-matrix)
- C-040: [Gauss decomposition of SL(2,C)](#gauss-decomposition-of-sl-2-c)
- C-041: [SU(2) disentangling (compact case)](#su-2-disentangling-compact-case)
- C-042: [SU(2) coherent overlap and annihilator](#su-2-coherent-overlap-and-annihilator)
- C-043: [Berezin-Lieb bounds for spin](#berezin-lieb-bounds-for-spin)
- C-044: [Spin-to-oscillator contraction](#spin-to-oscillator-contraction)
- C-045: [su(1,1) algebra, Casimir, discrete-series ladder and holomorphic realization](#su-1-1-algebra-casimir-discrete-series-ladder-and-holomorphic-realization)
- C-046: [Bosonic realizations of SU(1,1)](#bosonic-realizations-of-su-1-1)
- C-047: [SU(1,1) coherent state expansion and normalization](#su-1-1-coherent-state-expansion-and-normalization)
- C-048: [Single-mode squeeze: normal-ordered factorization](#single-mode-squeeze-normal-ordered-factorization); [SU(2) disentangling (compact case)](#su-2-disentangling-compact-case)
- C-049: [Single-mode squeeze: normal-ordered factorization](#single-mode-squeeze-normal-ordered-factorization); [Squeeze conjugation in block form](#squeeze-conjugation-in-block-form); [Coordinate form of the squeezed vacuum](#coordinate-form-of-the-squeezed-vacuum); [Squeezed number operator](#squeezed-number-operator); [Wigner function of a one-mode squeezed Gaussian](#wigner-function-of-a-one-mode-squeezed-gaussian)
- C-050: [Time-dependent SU(1,1) Hamiltonian](#time-dependent-su-1-1-hamiltonian)
- C-051: [Duhamel derivative of an exponential](#duhamel-derivative-of-an-exponential)
- C-052: [Singular oscillator lowest weight](#singular-oscillator-lowest-weight)
- C-053: [Weyl symbol of a metaplectic operator (Cayley form)](#weyl-symbol-of-a-metaplectic-operator-cayley-form); [Weyl symbol of the oscillator exponential](#weyl-symbol-of-the-oscillator-exponential)
- C-054: [Herman-Kluk representation, exact for quadratics](#herman-kluk-representation-exact-for-quadratics)
- C-055: [Gaussian density operators, purity, thermal factors](#gaussian-density-operators-purity-thermal-factors)
- C-056: [Reduced Gaussian channel from a joint quadratic evolution](#reduced-gaussian-channel-from-a-joint-quadratic-evolution)
- C-057: [Segal fields](#segal-fields)
- C-058: [Second quantization of a one-particle unitary](#second-quantization-of-a-one-particle-unitary)
- C-059: [Displacement: normal and antinormal factorization](#displacement-normal-and-antinormal-factorization); [Displacement: composition law](#displacement-composition-law); [Central-commutator BCH](#central-commutator-bch); [Coherent state expansion, eigenvalue equation, overlap and resolution](#coherent-state-expansion-eigenvalue-equation-overlap-and-resolution); [Vacuum displacement amplitude and coherent overlap](#vacuum-displacement-amplitude-and-coherent-overlap); [Infinite-mode displacement, exponential vector, product, overlap](#infinite-mode-displacement-exponential-vector-product-overlap)
- C-060: [Displacement conjugation of a field and of N](#displacement-conjugation-of-a-field-and-of-n)
- C-061: [Derivative of a moving displacement](#derivative-of-a-moving-displacement)
- C-062: [Non-formula records](#non-formula-records)
- C-063: [Non-formula records](#non-formula-records)
- C-064: [CAR and finite fermionic Fock space](#car-and-finite-fermionic-fock-space)
- C-065: [Linear fermionic exponential with complex coefficients](#linear-fermionic-exponential-with-complex-coefficients)
- C-066: [Grassmann derivatives: Leibniz rule and left/right relation](#grassmann-derivatives-leibniz-rule-and-left-right-relation)
- C-067: [Berezin integration: determinant rule and paired normalization](#berezin-integration-determinant-rule-and-paired-normalization)
- C-068: [Berezin integration: determinant rule and paired normalization](#berezin-integration-determinant-rule-and-paired-normalization); [Fermionic Gaussian integral with sources](#fermionic-gaussian-integral-with-sources)
- C-069: [Berezin integration: determinant rule and paired normalization](#berezin-integration-determinant-rule-and-paired-normalization); [Pfaffian Gaussian integral](#pfaffian-gaussian-integral)
- C-070: [Normal-ordered monomial expansion of fermionic operators](#normal-ordered-monomial-expansion-of-fermionic-operators)
- C-071: [Central-commutator BCH](#central-commutator-bch); [Grassmann Weyl translation](#grassmann-weyl-translation)
- C-072: [Fermionic coherent state, overlap, resolution](#fermionic-coherent-state-overlap-resolution)
- C-073: [Supertrace and Berezinian](#supertrace-and-berezinian)
- C-074: [Fermionic Weyl translation expansion](#fermionic-weyl-translation-expansion)
- C-075: [Normal, antinormal and s-ordered fermionic translations](#normal-antinormal-and-s-ordered-fermionic-translations)
- C-076: [Ordinary trace from mixed fermionic symbols](#ordinary-trace-from-mixed-fermionic-symbols)
- C-077: [Fermionic star product](#fermionic-star-product)
- C-078: [Fermi oscillator exponential and its Weyl symbol](#fermi-oscillator-exponential-and-its-weyl-symbol)
- C-079: [Fermionic partition determinant](#fermionic-partition-determinant)
- C-080: [BdG Hamiltonian and its canonical matrix](#bdg-hamiltonian-and-its-canonical-matrix)
- C-081: [Non-formula records](#non-formula-records)
- C-082: [Fermionic pair vacuum](#fermionic-pair-vacuum)
- C-083: [Fermionic Riccati and vacuum amplitude](#fermionic-riccati-and-vacuum-amplitude)
- C-084: [Fermionic Mehler kernel](#fermionic-mehler-kernel)
- C-085: [Fermionic Cayley Weyl symbol (corrected sign)](#fermionic-cayley-weyl-symbol-corrected-sign)
- C-086: [Fermionic translation coefficient and Pfaffian amplitude](#fermionic-translation-coefficient-and-pfaffian-amplitude)
- C-087: [SU(2) disentangling (compact case)](#su-2-disentangling-compact-case); [Fermion-pair realization of su(2)](#fermion-pair-realization-of-su-2)
- C-088: [Mixed boson-fermion displacement](#mixed-boson-fermion-displacement)
- C-089: [Super star product](#super-star-product)
- C-090: [Supercommutator](#supercommutator)
- C-091: [Commuting Gaussian integral and Fresnel factor](#commuting-gaussian-integral-and-fresnel-factor)
- C-092: [Non-formula records](#non-formula-records)
- C-093: [Non-formula records](#non-formula-records)
- C-094: [Adjoint exponential](#adjoint-exponential)
- C-095: [Non-formula records](#non-formula-records)
- C-096: [Non-formula records](#non-formula-records)
- C-097: [SU(1,1) coherent state expansion and normalization](#su-1-1-coherent-state-expansion-and-normalization); [Berezin symbol product on the disk](#berezin-symbol-product-on-the-disk)
- C-098: [Berezin integration: determinant rule and paired normalization](#berezin-integration-determinant-rule-and-paired-normalization); [Grassmann Fourier involution](#grassmann-fourier-involution)
- C-099: [Berezin integration: determinant rule and paired normalization](#berezin-integration-determinant-rule-and-paired-normalization); [Holomorphic Grassmann reproducing kernel](#holomorphic-grassmann-reproducing-kernel)
- G-001: [Canonical oscillator realization of the CCR](#canonical-oscillator-realization-of-the-ccr)
- G-002: [Hadamard conjugation series](#hadamard-conjugation-series)
- G-003: [Central-commutator BCH](#central-commutator-bch)
- G-004: [Number commutators and free rotation of the field](#number-commutators-and-free-rotation-of-the-field)
- G-005: [Fock basis and ladder action](#fock-basis-and-ladder-action)
- G-006: [Number-shift rules for functions of the number operator](#number-shift-rules-for-functions-of-the-number-operator); [Normal-ordered powers as falling factorials](#normal-ordered-powers-as-falling-factorials)
- G-007: [Quadrature commutator and number-state variance](#quadrature-commutator-and-number-state-variance)
- G-008: [Independent-mode CCR and multimode Fock states](#independent-mode-ccr-and-multimode-fock-states)
- G-009: [Number-conserving exponential in normal order](#number-conserving-exponential-in-normal-order); [Thermal state as a number exponential](#thermal-state-as-a-number-exponential)
- G-010: [Thermal factorial moments and number variance](#thermal-factorial-moments-and-number-variance)
- G-011: [No Hermitian phase operator canonically conjugate to number](#no-hermitian-phase-operator-canonically-conjugate-to-number)
- G-012: [Susskind-Glogower one-sided shift](#susskind-glogower-one-sided-shift)
- G-013: [Coherent state expansion, eigenvalue equation, overlap and resolution](#coherent-state-expansion-eigenvalue-equation-overlap-and-resolution)
- G-014: [No normalizable right eigenvector of the creation operator](#no-normalizable-right-eigenvector-of-the-creation-operator)
- G-015: [Displacement: normal and antinormal factorization](#displacement-normal-and-antinormal-factorization)
- G-016: [Coherent state expansion, eigenvalue equation, overlap and resolution](#coherent-state-expansion-eigenvalue-equation-overlap-and-resolution); [Vacuum action of the ordered displacement](#vacuum-action-of-the-ordered-displacement)
- G-017: [Displacement: normal and antinormal factorization](#displacement-normal-and-antinormal-factorization); [Displacement: composition law](#displacement-composition-law)
- G-018: [Displacement: adjoint action on a](#displacement-adjoint-action-on-a); [Hadamard conjugation series](#hadamard-conjugation-series)
- G-019: [Vacuum displacement amplitude and coherent overlap](#vacuum-displacement-amplitude-and-coherent-overlap)
- G-020: [Normally ordered coherent matrix elements](#normally-ordered-coherent-matrix-elements)
- G-021: [Coherent state expansion, eigenvalue equation, overlap and resolution](#coherent-state-expansion-eigenvalue-equation-overlap-and-resolution); [Coherent-state photon statistics](#coherent-state-photon-statistics)
- G-022: [Displacement: composition law](#displacement-composition-law); [Forced oscillator: displacement plus accumulated phase](#forced-oscillator-displacement-plus-accumulated-phase)
- G-023: [Bargmann realization and displacement as weighted translation](#bargmann-realization-and-displacement-as-weighted-translation)
- G-024: [Fock matrix elements from the coherent diagonal symbol](#fock-matrix-elements-from-the-coherent-diagonal-symbol)
- G-025: [Coherent-projector differential correspondences](#coherent-projector-differential-correspondences)
- G-026: [P representation, optical equivalence, and Q function](#p-representation-optical-equivalence-and-q-function)
- G-027: [Glauber P function of a number state](#glauber-p-function-of-a-number-state)
- G-028: [Ordered characteristic functions differ by Gaussian factors](#ordered-characteristic-functions-differ-by-gaussian-factors)
- G-029: [Normally ordered moments from the normal characteristic function](#normally-ordered-moments-from-the-normal-characteristic-function)
- G-030: [Ordered characteristic functions differ by Gaussian factors](#ordered-characteristic-functions-differ-by-gaussian-factors)
- G-031: [Complex Fourier transform between P and the normal characteristic function](#complex-fourier-transform-between-p-and-the-normal-characteristic-function)
- G-032: [Thermal characteristic functions and phase-space densities](#thermal-characteristic-functions-and-phase-space-densities)
- G-033: [Q and Wigner functions of coherent and number states](#q-and-wigner-functions-of-coherent-and-number-states)
- G-034: [Weyl operator from displaced parity, Weyl symbol and reconstruction](#weyl-operator-from-displaced-parity-weyl-symbol-and-reconstruction)
- G-035: [Normally ordered correlation functions](#normally-ordered-correlation-functions)
- G-036: [Second-order coherence and Mandel Q](#second-order-coherence-and-mandel-q)
- G-037: [Atomic transition algebra and Jaynes-Cummings excitation number](#atomic-transition-algebra-and-jaynes-cummings-excitation-number)
- G-038: [Exact resonant Jaynes-Cummings propagator](#exact-resonant-jaynes-cummings-propagator)
- G-039: [Detuned Jaynes-Cummings doublet propagator](#detuned-jaynes-cummings-doublet-propagator)
- G-040: [Dispersive Hamiltonian as conditional number rotation](#dispersive-hamiltonian-as-conditional-number-rotation)
- G-041: [Raman number-controlled atomic rotation](#raman-number-controlled-atomic-rotation)
- G-042: [Multiphoton Jaynes-Cummings couplings](#multiphoton-jaynes-cummings-couplings)
- G-043: [Passive two-port transform must be unitary](#passive-two-port-transform-must-be-unitary)
- G-044: [Hadamard conjugation series](#hadamard-conjugation-series); [Beam splitter with i-reflection: exponential and adjoint action](#beam-splitter-with-i-reflection-exponential-and-adjoint-action)
- G-045: [Real beam-splitter rotation](#real-beam-splitter-rotation)
- G-046: [SU(2) disentangling (compact case)](#su-2-disentangling-compact-case); [Schwinger SU(2) algebra, Casimir and J3 phase action](#schwinger-su-2-algebra-casimir-and-j3-phase-action)
- G-047: [Beam-splitter output of a number state (binomial state)](#beam-splitter-output-of-a-number-state-binomial-state)
- G-048: [Coherent products stay products under passive mixing](#coherent-products-stay-products-under-passive-mixing)
- G-049: [Two-photon interference and twin-Fock output at 50:50](#two-photon-interference-and-twin-fock-output-at-50-50)
- G-050: [Phase shifter and parity exponential](#phase-shifter-and-parity-exponential)
- G-051: [Single-mode squeeze: normal-ordered factorization](#single-mode-squeeze-normal-ordered-factorization); [Single-mode squeeze operator](#single-mode-squeeze-operator)
- G-052: [Single-mode squeeze: adjoint action on a](#single-mode-squeeze-adjoint-action-on-a); [Hadamard conjugation series](#hadamard-conjugation-series)
- G-053: [Single-mode su(1,1) realization and Casimir](#single-mode-su-1-1-realization-and-casimir)
- G-054: [Single-mode squeeze: normal-ordered factorization](#single-mode-squeeze-normal-ordered-factorization)
- G-055: [Squeeze operator: reversed (antinormal) Gauss factorization](#squeeze-operator-reversed-antinormal-gauss-factorization)
- G-056: [Single-mode squeezed vacuum](#single-mode-squeezed-vacuum)
- G-057: [Rotated-quadrature variance of squeezed vacuum](#rotated-quadrature-variance-of-squeezed-vacuum)
- G-058: [Quadrature covariance of squeezed vacuum](#quadrature-covariance-of-squeezed-vacuum)
- G-059: [Displacement: adjoint action on a](#displacement-adjoint-action-on-a); [Displaced squeezed state: low moments](#displaced-squeezed-state-low-moments)
- G-060: [Interchanging displacement and squeezing](#interchanging-displacement-and-squeezing)
- G-061: [Displacement: adjoint action on a](#displacement-adjoint-action-on-a); [Single-mode squeeze: normal-ordered factorization](#single-mode-squeeze-normal-ordered-factorization); [Displaced squeezed state as a Bogoliubov eigenstate](#displaced-squeezed-state-as-a-bogoliubov-eigenstate)
- G-062: [Single-mode squeezed vacuum](#single-mode-squeezed-vacuum); [Fock coefficients of a displaced squeezed vacuum (Hermite form)](#fock-coefficients-of-a-displaced-squeezed-vacuum-hermite-form)
- G-063: [Single-mode squeezed vacuum](#single-mode-squeezed-vacuum); [Coherent overlap of squeezed vacuum and Q function of a displaced squeezed state](#coherent-overlap-of-squeezed-vacuum-and-q-function-of-a-displaced-squeezed-state)
- G-064: [Wigner function of a displaced squeezed state](#wigner-function-of-a-displaced-squeezed-state)
- G-065: [Number variance of a displaced squeezed vacuum](#number-variance-of-a-displaced-squeezed-vacuum)
- G-066: [Degenerate parametric amplifier gives a squeeze of parameter 2 eta t](#degenerate-parametric-amplifier-gives-a-squeeze-of-parameter-2-eta-t)
- G-067: [Degenerate parametric amplifier gives a squeeze of parameter 2 eta t](#degenerate-parametric-amplifier-gives-a-squeeze-of-parameter-2-eta-t)
- G-068: [Two-mode squeeze: normal-ordered factorization](#two-mode-squeeze-normal-ordered-factorization); [Hadamard conjugation series](#hadamard-conjugation-series); [Two-mode squeeze operator and Bogoliubov action](#two-mode-squeeze-operator-and-bogoliubov-action)
- G-069: [Two-mode su(1,1) realization and Casimir](#two-mode-su-1-1-realization-and-casimir)
- G-070: [Two-mode squeeze: normal-ordered factorization](#two-mode-squeeze-normal-ordered-factorization)
- G-071: [Two-mode squeeze: normal-ordered factorization](#two-mode-squeeze-normal-ordered-factorization); [Two-mode squeeze as two opposite single-mode squeezes of the sum and difference modes](#two-mode-squeeze-as-two-opposite-single-mode-squeezes-of-the-sum-and-difference-modes)
- G-072: [Two-mode squeezed vacuum](#two-mode-squeezed-vacuum)
- G-073: [Rotated-quadrature variance of squeezed vacuum](#rotated-quadrature-variance-of-squeezed-vacuum); [Twin-beam number statistics, collective squeezing, and the intermode Cauchy-Schwarz test](#twin-beam-number-statistics-collective-squeezing-and-the-intermode-cauchy-schwarz-test)
- G-074: [Reduced twin-beam state is thermal; entanglement entropy](#reduced-twin-beam-state-is-thermal-entanglement-entropy)
- G-075: [Two-mode squeeze: normal-ordered factorization](#two-mode-squeeze-normal-ordered-factorization); [Nondegenerate parametric amplifier gives a two-mode squeeze of parameter eta t](#nondegenerate-parametric-amplifier-gives-a-two-mode-squeeze-of-parameter-eta-t)
- G-076: [Phase shifter and parity exponential](#phase-shifter-and-parity-exponential); [Self-Kerr evolution of a coherent state](#self-kerr-evolution-of-a-coherent-state)
- G-077: [Half-revival cat identity](#half-revival-cat-identity)
- G-078: [Normal-ordered powers as falling factorials](#normal-ordered-powers-as-falling-factorials); [Normally ordered Kerr differs by a rotation](#normally-ordered-kerr-differs-by-a-rotation)
- G-079: [Number-shift rules for functions of the number operator](#number-shift-rules-for-functions-of-the-number-operator); [Number-function conjugation and Kerr Heisenberg solution](#number-function-conjugation-and-kerr-heisenberg-solution)
- G-080: [Number-shift rules for functions of the number operator](#number-shift-rules-for-functions-of-the-number-operator); [Kerr coherent-state moments and quadrature variance](#kerr-coherent-state-moments-and-quadrature-variance)
- G-081: [Cross-Kerr conditional rotation](#cross-kerr-conditional-rotation)
- G-082: [Coherent-state cat: normalization and a^2 eigenvalue](#coherent-state-cat-normalization-and-a-2-eigenvalue)
- G-083: [Cat photon-number distribution](#cat-photon-number-distribution)
- G-084: [Cat Wigner function](#cat-wigner-function)
- G-085: [Coherent-state cat: normalization and a^2 eigenvalue](#coherent-state-cat-normalization-and-a-2-eigenvalue); [Non-formula records](#non-formula-records)
- G-086: [Beam splitter with i-reflection: exponential and adjoint action](#beam-splitter-with-i-reflection-exponential-and-adjoint-action); [Cross-Kerr gate conjugated by a 50:50 beam splitter](#cross-kerr-gate-conjugated-by-a-50-50-beam-splitter)
- G-087: [Linear quadrature exponential: normal-order conversion factor](#linear-quadrature-exponential-normal-order-conversion-factor)
- G-088: [Centered quadrature moments from normal moments](#centered-quadrature-moments-from-normal-moments)
- G-089: [Squared-field quadratures and their commutators](#squared-field-quadratures-and-their-commutators)
- G-090: [Pair coherent states](#pair-coherent-states)
- G-091: [Twin-beam number statistics, collective squeezing, and the intermode Cauchy-Schwarz test](#twin-beam-number-statistics-collective-squeezing-and-the-intermode-cauchy-schwarz-test)
- G-092: [Frequency conversion as an SU(2) rotation](#frequency-conversion-as-an-su-2-rotation)
- G-093: [Displacement: adjoint action on a](#displacement-adjoint-action-on-a); [Photon subtraction and addition](#photon-subtraction-and-addition)
- G-094: [Collective broadband commutator](#collective-broadband-commutator)
- G-095: [Low-gain down-conversion to second order](#low-gain-down-conversion-to-second-order)
- G-096: [No-jump number exponential on a coherent state](#no-jump-number-exponential-on-a-coherent-state)
- G-097: [Number-conserving exponential in normal order](#number-conserving-exponential-in-normal-order); [Normal-ordered powers as falling factorials](#normal-ordered-powers-as-falling-factorials); [Thermal state in normally ordered exponential form](#thermal-state-in-normally-ordered-exponential-form)
- G-098: [Amplitude damping: Kraus form and matrix-element solution](#amplitude-damping-kraus-form-and-matrix-element-solution)
- G-099: [Damping of coherent dyads and of the normal characteristic function](#damping-of-coherent-dyads-and-of-the-normal-characteristic-function)
- G-100: [Driven damped oscillator: coherent amplitude and steady state](#driven-damped-oscillator-coherent-amplitude-and-steady-state)
- G-101: [Displacement: normal and antinormal factorization](#displacement-normal-and-antinormal-factorization); [Number commutators and free rotation of the field](#number-commutators-and-free-rotation-of-the-field); [Recoil exponential as a displacement and its interaction-picture form](#recoil-exponential-as-a-displacement-and-its-interaction-picture-form)
- G-102: [Second-order Dyson elimination of a fast oscillating coupling](#second-order-dyson-elimination-of-a-fast-oscillating-coupling)
- G-103: [Effective Kerr Hamiltonian from an off-resonant parametric coupling](#effective-kerr-hamiltonian-from-an-off-resonant-parametric-coupling)
- G-104: [Vacuum displacement amplitude and coherent overlap](#vacuum-displacement-amplitude-and-coherent-overlap); [Reduced twin-beam state is thermal; entanglement entropy](#reduced-twin-beam-state-is-thermal-entanglement-entropy); [Coherent trace formula and Schmidt spectrum](#coherent-trace-formula-and-schmidt-spectrum)
- G-105: [Continuum CCR and first-order pair-creation kernel](#continuum-ccr-and-first-order-pair-creation-kernel)
