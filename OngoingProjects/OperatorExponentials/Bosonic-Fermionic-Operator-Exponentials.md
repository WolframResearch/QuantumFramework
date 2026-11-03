# Bosonic and Fermionic Operator Exponentials

A formula, theorem, and derivation reference from four locally supplied books.

**Compiled:** 2 November 2026. **Original reading record:** **1,483 supplied physical PDF pages**, including chapters, exercises, appendices, references, and indexes. This process record does not certify formula-by-formula coverage or correctness. **Collection:** 419 stable source entries, plus 9 cross-book atlas entries. The reading and source-display checks are documented in the book dossiers below.

**Correctness status:** the independent [correctness and coverage audit](Correctness-and-Coverage-Audit.md) records confirmed corrections, runnable checks and unresolved verification limits. A source locator, a full reading pass, and a passing numerical example provide different kinds of evidence. This reference is not an exhaustive transcription or an independently proved theorem library.

The focus is noncommutative creation–annihilation algebra and functions of its operators: BCH/Hadamard identities, ordering, displacement, single- and multimode squeezing, bosonic and fermionic canonical transformations, quadratic exponentials, Gaussian kernels and propagators, Grassmann/Pfaffian calculus, Wick and time-ordered expansions, traces, determinant phases, and useful nonlinear examples such as Kerr evolution. Material outside this focus is recorded in the coverage ledgers rather than copied indiscriminately.

## How to find and recheck a result

- Each entry has a stable ID: **B-…**, **P-… / P2-…**, **C-…**, or **G-…**. **X-…** identifies the common-convention atlas, not a book equation number.
- A source locator gives the book section, the **printed page**, the **one-based PDF page**, and the book's equation/theorem/proposition/exercise number when one exists. Roman and unnumbered preliminaries are kept separate from Arabic body pages.
- **Source** describes a statement actually in the book. **Re-expression** translates notation. **Derived / Consequence / Correction** identifies a deduction supplied in this reference; it must not be quoted as if printed in the book. Unnumbered source theorems are described as unnumbered.
- File-page offsets are copy-specific. Combescure–Robert has **nonuniform offsets**; use the checked page pair in each C entry.
- A full reading pass is distinct from independent proof checking. Important equation pages were inspected visually, selected identities were checked algebraically or with small matrices, and visible source errors are flagged. The dossiers state limits on formula transcription, formal manipulations, and unresolved interpretations rather than presenting OCR as a proof.

## The four source books

| Key | Bibliographic reference | Physical PDF pages | Arabic page mapping |
|---|---|---:|---|
| **B** | F. A. Berezin, *The Method of Second Quantization*. Translated by Nobumichi Mugibayashi and Alan Jeffrey. Academic Press, 1966. Pure and Applied Physics, vol. 24. ISBN 9780120894505. | 236 | PDF = printed + 8 |
| **P** | A. Perelomov, *Generalized Coherent States and Their Applications*. Springer, 1986. ISBN 9783540159124. | 332 | PDF = printed + 11 |
| **C** | Didier Robert and Monique Combescure, *Coherent States and Applications in Mathematical Physics*, second edition. Springer, 2021. DOI 10.1007/978-3-030-70845-0. | 582 | Nonuniform; use entry locators |
| **G** | Christopher C. Gerry and Peter L. Knight, *Introductory Quantum Optics*. Cambridge University Press. This copy states first print publication 2004 and electronic copyright 2005. ISBN 9780521820356 (hardback). | 333 | PDF = printed + 16 |

## Topic index

| Topic | Starting entries |
|---|---|
| BCH, Hadamard, and differentiation of exponentials | [P-002](#p-002), [G-002](#g-002), [G-003](#g-003), [C-002](#c-002), [C-051](#c-051), [C-094](#c-094), [B-072](#b-072) |
| Displacement ordering, Weyl products, and differential correspondences | [X-001](#x-001), [P-003](#p-003), [P-010](#p-010), [P-011](#p-011), [G-015](#g-015), [G-017](#g-017), [B-050](#b-050), [C-003](#c-003) |
| Single-mode squeezing and both Gauss orderings | [X-002](#x-002), [P-029](#p-029), [P-031](#p-031), [G-051](#g-051), [G-054](#g-054), [G-055](#g-055), [C-025](#c-025), [C-048](#c-048) |
| Two-mode squeezing and pair statistics | [X-003](#x-003), [P-032](#p-032), [G-068](#g-068), [G-070](#g-070), [G-071](#g-071), [G-072](#g-072), [P2-019](#p2-019) |
| SU(2), passive mode mixing, and fermionic pair rotations | [X-004](#x-004), [P-025](#p-025), [G-044](#g-044), [G-046](#g-046), [C-041](#c-041), [C-087](#c-087), [P2-025](#p2-025) |
| General bosonic canonical transformations and implementability | [X-009](#x-009), [P-034](#p-034), [P-035](#p-035), [B-042](#b-042), [B-046](#b-046), [B-049](#b-049), [C-024](#c-024) |
| General fermionic canonical transformations and parity | [X-009](#x-009), [P-041](#p-041), [P-042](#p-042), [B-056](#b-056), [B-058](#b-058), [B-059](#b-059), [C-081](#c-081), [P2-026](#p2-026) |
| Multimode Gaussian vacua, overlaps, and determinant branches | [X-006](#x-006), [P-037](#p-037), [P-043](#p-043), [B-048](#b-048), [B-061](#b-061), [C-082](#c-082), [P2-021](#p2-021), [P2-027](#p2-027) |
| Quadratic propagators, Riccati equations, and metaplectic/spin phases | [X-005](#x-005), [B-069](#b-069), [B-071](#b-071), [B-074](#b-074), [C-018](#c-018), [C-021](#c-021), [C-083](#c-083), [P2-024](#p2-024) |
| Forced and parametric oscillators | [G-022](#g-022), [G-066](#g-066), [C-032](#c-032), [P2-010](#p2-010), [P2-013](#p2-013), [P2-018](#p2-018) |
| Normal symbols, kernel multiplication, and normal-order exponentials | [P-015](#p-015), [B-010](#b-010), [B-014](#b-014), [B-027](#b-027), [B-054](#b-054), [B-055](#b-055), [G-097](#g-097) |
| Weyl/antinormal symbols, Moyal products, and characteristic functions | [P-019](#p-019), [P-020](#p-020), [G-028](#g-028), [G-030](#g-030), [C-016](#c-016), [C-017](#c-017), [C-075](#c-075), [C-077](#c-077) |
| Grassmann calculus, fermionic coherent displacement, and traces | [C-098](#c-098), [C-099](#c-099), [B-029](#b-029), [B-030](#b-030), [B-033](#b-033), [B-039](#b-039), [B-040](#b-040), [C-065](#c-065), [C-071](#c-071), [C-072](#c-072), [C-073](#c-073) |
| Commuting and anticommuting Gaussian integrals, Pfaffians | [B-022](#b-022), [B-023](#b-023), [B-033](#b-033), [B-036](#b-036), [C-068](#c-068), [C-069](#c-069), [C-091](#c-091) |
| Thermal exponentials, partition determinants, and trace bounds | [X-007](#x-007), [G-009](#g-009), [G-097](#g-097), [B-106](#b-106), [B-107](#b-107), [P-016](#p-016), [C-055](#c-055), [C-079](#c-079) |
| Wick contractions, Dyson series, and auxiliary fields | [B-094](#b-094), [B-096](#b-096), [B-098](#b-098), [B-099](#b-099), [B-101](#b-101), [B-103](#b-103), [B-105](#b-105), [C-063](#c-063) |
| Nonlinear exponentials: number-shift, Kerr, cats, and limitations | [G-006](#g-006), [G-076](#g-076), [G-077](#g-077), [G-079](#g-079), [G-081](#g-081), [G-087](#g-087), [G-089](#g-089), [X-008](#x-008) |
| Fock matrix elements and generating functions | [P-018](#p-018), [G-062](#g-062), [G-090](#g-090), [P2-011](#p2-011) |
| Infinite modes, singular charts, operator domains, and stability | [B-006](#b-006), [B-011](#b-011), [B-066](#b-066), [B-079](#b-079), [B-080](#b-080), [B-083](#b-083), [C-057](#c-057), [P-004](#p-004), [P-040](#p-040) |
| Lie-group decomposition, root systems, and generalized coherent states | [P-021](#p-021), [P-023](#p-023), [P-046](#p-046), [P-047](#p-047), [P-048](#p-048), [C-094](#c-094), [C-095](#c-095), [C-096](#c-096) |

## Notation: operators, scalars, and coefficient matrices

**Hats mark operators acting on states or operators.** Scalars, scalar functions, Grassmann parameters, and finite numerical coefficient matrices have no hats. This distinction is applied to the transcribed formulas throughout the reference, including the operator labels in scalar symbols and kernels.

| Example notation | Type |
|---|---|
| $\hat a_i,\hat a_i^\dagger,\hat c_i,\hat c_i^\dagger$ | Bosonic/fermionic annihilation and creation operators. |
| $\hat H,\hat N,\hat\rho,\hat D(\alpha),\hat S(\xi),\hat U(t),\hat J_\pm,\hat K_0$ | Hamiltonian, number, density, displacement, squeeze, evolution, and represented Lie-algebra operators. |
| $\hat I$ or $\hat I_{\mathcal H}$; $I_n$ | Hilbert-space identity **operator**; numerical identity **matrix**, respectively. Scalar shifts are written explicitly, e.g. $\hat N+\tfrac12\hat I$. |
| $\alpha,\beta,\xi,r,t,\hbar$; $z,\bar z$ | Scalar parameters or commuting phase-space coordinates. A physical coordinate/momentum operator is instead $\hat q,\hat p$ or $\hat Q,\hat P$. |
| $\theta,\bar\theta$ when defined as Grassmann labels | Anticommuting **scalars**, not Fock-space operators. Their parity and left/right derivative conventions still matter. |
| $U,V,Z,F_t$ when defined as canonical/flow/Gaussian coefficient matrices | Finite matrices with scalar entries. The hatted implementer $\hat U$ is a different object from the unhatted canonical block $U$. Local definitions settle overloaded letters. |
| $A_N(\bar z,z),Q_{\hat A}(z),K_{\hat A}(\bar z,w),\rho_{mn}$ | Scalar symbol, scalar covariant symbol, scalar generating kernel, and scalar density-matrix elements; $\hat A$ and $\hat\rho$ themselves are operators. |
| $\hat\Phi,\hat\Psi,\hat C$ in one-particle Hilbert-space theorems | One-particle **operators**, including the infinite-dimensional maps in Berezin. Their finite numerical matrix specializations are identified locally. |
| $\widehat{\mathcal L},\widehat{\operatorname{ad}}_{\hat H}$ | Superoperators acting on operator space. |

State vectors and bras/kets do not acquire hats. Traces, determinants, exponentials, scalar special functions, and explicit differential operations such as $\partial_z$ retain their standard notation; an identification such as $\hat a\leftrightarrow\partial_z$ makes the differential realization explicit. The scalar Gamma function $\Gamma(x)$ and a pairing matrix $\Gamma_t$ are different from the Fock-space operator $\widehat\Gamma(M)$. A zero on the right of an operator identity denotes the zero operator.

## Conventions that must be fixed before using a formula

1. **Conjugation direction:** $\hat U\hat a\hat U^{-1}$ and $\hat U^{-1}\hat a\hat U$ give inverse canonical matrices. A squeeze or displacement parameter cannot be compared until that direction is fixed.
2. **Bosons versus fermions:** bosonic pair matrices are symmetric; fermionic pair matrices are skew-symmetric. Fermionic pair annihilation uses the reverse of the creation-pair order. Even quadratic generators use ordinary commutators; odd Grassmann-labelled objects also require parity and left/right derivative conventions.
3. **Operator versus symbol:** $e^{\hat A}$, a normally ordered $:\!e^{\hat A}\!:$, a scalar kernel exponential, a covariant symbol, and a Weyl symbol are distinct objects. In particular a Gaussian symbol does not by itself specify which operator ordering is intended.
4. **Scalar phase and square-root branch:** a canonical matrix determines adjoint action but leaves a scalar implementing phase undetermined. Determinants, Pfaffians, vacuum amplitudes, connected charts, and the metaplectic/spin cover supply information that the matrix alone loses.
5. **Analytic hypotheses:** finite fermionic Fock space is finite-dimensional. Bosonic generators are generally unbounded; infinite-mode limits additionally need Hilbert–Schmidt/trace-class and domain assumptions. Time ordering and formal Wick expansions retain their qualifications in the source records.

## Frequently needed formulas in a common convention


This short atlas connects the book-specific dossiers. Hats mark operators; finite numerical coefficient matrices, scalar functions, amplitudes, and state labels remain unhatted. Here $\hat I$ is the Hilbert-space identity and $I_N$ is a numerical matrix identity. Its equations are re-expressions or explicitly identified deductions, with links to the detailed source records. It is not a replacement for their assumptions, convergence restrictions, or source-error notes.

<a id="x-001"></a>

### X-001 — Displacement: the basic linear exponential

**Source records:** [P-002](#p-002), [P-003](#p-003), [P-010](#p-010), [G-015](#g-015), [G-017](#g-017), [B-050](#b-050), [C-006](#c-006).

Fix

$$
\hat D(\alpha)=e^{\alpha \hat a^\dagger-\bar\alpha \hat a},\qquad [\hat a,\hat a^\dagger]=\hat I.
$$

Then

$$
\hat D(\alpha)=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}e^{-\bar\alpha \hat a}
=e^{|\alpha|^2/2}e^{-\bar\alpha \hat a}e^{\alpha \hat a^\dagger},
$$

$$
\hat D(\alpha)^\dagger \hat a\hat D(\alpha)=\hat a+\alpha\hat I,\qquad
\hat D(\alpha)\hat D(\beta)=e^{(\alpha\bar\beta-\bar\alpha\beta)/2}\hat D(\alpha+\beta).
$$

The scalar is the central BCH correction. For independent modes replace scalar products by sums over modes. Berezin's implementer is often oriented as $\hat U\hat a\hat U^{-1}$, which reverses the usual displacement parameter when compared to $\hat D^\dagger \hat a\hat D$.

<a id="x-002"></a>

### X-002 — Single-mode squeezing: a quadratic exponential

**Source records:** [P-028](#p-028)–[P-031](#p-031), [G-051](#g-051)–[G-061](#g-061), [B-051](#b-051), [C-048](#c-048)–[C-049](#c-049).

Fix $\xi=re^{i\theta}$, $r\ge0$, and

$$
\hat S(\xi)=\exp\!\left[\frac12\left(\bar\xi \hat a^2-\xi(\hat a^\dagger)^2\right)\right].
$$

Then

$$
\hat S(\xi)=e^{-\frac12e^{i\theta}\tanh r\,(\hat a^\dagger)^2}
e^{-\log\cosh r\,(\hat N+\tfrac12\hat I)}
e^{\frac12e^{-i\theta}\tanh r\,\hat a^2},\qquad \hat N=\hat a^\dagger \hat a,
$$

$$
\hat S^\dagger \hat a\hat S=\hat a\cosh r-e^{i\theta}\hat a^\dagger\sinh r,
\qquad
\hat S|0\rangle=(\cosh r)^{-1/2}
e^{-\frac12e^{i\theta}\tanh r\,(\hat a^\dagger)^2}|0\rangle.
$$

The $1/2$ in the number exponent is needed for the vacuum normalization. For the reversed convention $\hat S^{-1}$ every squeeze parameter changes sign. The actual middle signs printed in Perelomov (5.2.6)/(5.2.8) are corrected in [P-029](#p-029).

<a id="x-003"></a>

### X-003 — Two-mode bosonic squeezing

**Source records:** [P-032](#p-032), [G-068](#g-068)–[G-075](#g-075); common group identity [P-029](#p-029).

For independent bosons $\hat a,\hat b$, in the convention consistent with X-002,

$$
\hat S_2(\xi)=e^{\bar\xi \hat a\hat b-\xi \hat a^\dagger \hat b^\dagger}
=e^{-e^{i\theta}\tanh r\,\hat a^\dagger \hat b^\dagger}
e^{-\log\cosh r\,(\hat N_a+\hat N_b+\hat I)}
e^{e^{-i\theta}\tanh r\,\hat a\hat b},
$$

$$
\hat S_2|0,0\rangle=\operatorname{sech}r\sum_{n\ge0}
(-e^{i\theta}\tanh r)^n|n,n\rangle.
$$

The closed algebra uses $\hat K_+=\hat a^\dagger \hat b^\dagger$, $\hat K_-=\hat a\hat b$, $\hat K_0=(\hat N_a+\hat N_b+\hat I)/2$. The number difference is conserved. The inverse-sign two-mode convention in P-032 is obtained by $\xi\mapsto-\xi$.

<a id="x-004"></a>

### X-004 — Compact versus noncompact disentangling

**Source records:** [P-025](#p-025), [P-029](#p-029), [C-041](#c-041), [C-048](#c-048), [G-046](#g-046), [P-041](#p-041)–[P-043](#p-043), [P2-025](#p2-025).

For $[\hat J_0,\hat J_\pm]=\pm \hat J_\pm$, $[\hat J_+,\hat J_-]=2\hat J_0$, let $\xi=re^{i\phi}$, $\zeta=e^{i\phi}\tan r$. On a chart with $\cos r>0$,

$$
e^{\xi \hat J_+-\bar\xi \hat J_-}
=e^{\zeta \hat J_+}e^{\log(1+|\zeta|^2)\hat J_0}e^{-\bar\zeta \hat J_-}.
$$

For $[\hat K_0,\hat K_\pm]=\pm \hat K_\pm$, $[\hat K_-,\hat K_+]=2\hat K_0$, put instead $\zeta=e^{i\phi}\tanh r$:

$$
e^{\xi \hat K_+-\bar\xi \hat K_-}
=e^{\zeta \hat K_+}e^{\log(1-|\zeta|^2)\hat K_0}e^{-\bar\zeta \hat K_-}.
$$

**Derived fermion-pair specialization:** for two CAR modes, set

$$
\hat J_+=\hat c_1^\dagger \hat c_2^\dagger,\quad \hat J_-=\hat c_2\hat c_1=\hat J_+^\dagger,
\quad \hat J_0=(\hat N_1+\hat N_2-\hat I)/2.
$$

These satisfy the compact algebra, so

$$
e^{\xi \hat c_1^\dagger \hat c_2^\dagger-\bar\xi \hat c_2\hat c_1}|00\rangle
=\cos r\,|00\rangle+e^{i\phi}\sin r\,|11\rangle.
$$

The order $\hat c_2\hat c_1$ is essential. Replacing bosonic commutators by anticommutators inside a BCH formula is not a valid derivation. Quadratic **even** fermionic generators close under ordinary commutators; Grassmann-labelled **odd** coherent-state constructions use graded rules documented in the fermionic dossiers.

<a id="x-005"></a>

### X-005 — General quadratic exponentials reduce to finite matrix evolution

**Source records:** [P-034](#p-034)–[P-035](#p-035), [P-041](#p-041)–[P-042](#p-042), [B-071](#b-071)–[B-074](#b-074), [C-018](#c-018)–[C-024](#c-024).

If an operator column $\hat R$ closes under a scalar-coefficient commutator,

$$
[\hat Q,\hat R_i]=\sum_jL_{ij}\hat R_j,
$$

then the Hadamard series gives

$$
e^{t\hat Q}\hat R e^{-t\hat Q}=e^{tL}\hat R.
$$

Use $\hat R=(\hat a,\hat a^\dagger)^T$ for quadratic bosonic generators or $\hat R=(\hat c,\hat c^\dagger)^T$ for quadratic even fermionic generators. Affine generators are treated by adjoining the identity to $\hat R$. This determines the canonical transformation, but an overall implementing scalar/phase is lost in the adjoint action; vacuum normalization, determinant cocycles and continuous branch choices supply it. Time-dependent noncommuting generators require time ordering or an equivalent matrix ODE, not a plain exponential of their time integral.

<a id="x-006"></a>

### X-006 — Finite-mode bosonic and fermionic Gaussian vacua

**Source records:** [P-036](#p-036)–[P-045](#p-045), [B-048](#b-048), [B-061](#b-061); source-specific Grassmann realizations in the C dossier.

In the source convention $(\hat a+Z\hat a^\dagger)|Z\rangle=0$ or $(\hat c+Z\hat c^\dagger)|Z\rangle=0$,

$$
|Z\rangle_B=\det(I_N-ZZ^\dagger)^{1/4}
e^{-\frac12\hat a^{\dagger T}Z\hat a^\dagger}|0\rangle,
\quad Z^T=Z,\quad I_N-ZZ^\dagger>0,
$$

$$
|Z\rangle_F=\det(I_N+ZZ^\dagger)^{-1/4}
e^{-\frac12\hat c^{\dagger T}Z\hat c^\dagger}|0\rangle,
\quad Z^T=-Z.
$$

The bosonic normalizability boundary is real; the fermionic finite-mode exponential is a polynomial, but its vacuum chart can fail when $U$ is singular. These pair-coherent fermionic states are distinct from Grassmann-labelled annihilation-eigenstate coherent states. The latter are treated separately in B-029–B-040 and the fermionic chapter of C.

<a id="x-007"></a>

### X-007 — Number-conserving exponentials, normal order, and traces

**Source records:** [B-013](#b-013), [B-028](#b-028), [B-040](#b-040), [B-107](#b-107), [G-009](#g-009), [G-097](#g-097).

Here $h,M\in\mathbb C^{N\times N}$ are finite numerical coefficient matrices, while $\widehat\Gamma(M)$ acts on Fock space. For finite modes, the second-quantization identity yields the normal-ordered re-expression

$$
e^{\hat a^\dagger h \hat a}=:\!\exp[\hat a^\dagger(e^h-I_N)\hat a]\!:
$$

and the same number-conserving identity holds for fermions with CAR normal ordering. In particular,

$$
|0\rangle\langle0|=:\!e^{-\sum_i\hat N_i}\!:
$$

for a finite bosonic mode set; for fermions it equals $\prod_i(\hat I-\hat N_i)$. The thermal traces are

$$
\operatorname{Tr}_{\mathcal F_B}\widehat\Gamma(M)=\det(I_N-M)^{-1},\qquad
\operatorname{Tr}_{\mathcal F_F}\widehat\Gamma(M)=\det(I_N+M).
$$

Finite fermion dimension needs no convergence restriction. For a positive bosonic thermal trace, $0\le M<I_N$ (equivalently $\|M\|<1$ in finite dimension) is required even for a single mode; otherwise the occupation sum diverges. For the infinite-mode extension, replace the numerical matrix $M$ by the one-particle Hilbert-space operator $\hat M$ and $I_N$ by $\hat I_{\mathcal H_1}$ in the Fredholm determinants. In the usual positive thermal setting, require trace-class $\hat M$ and additionally $\|\hat M\|<1$ for bosons. Analytic continuations or nonpositive $M$ need their own hypotheses. A fermionic trace is not obtained by blindly copying a bosonic coherent integral without its sign reversal at the boundary.

<a id="x-008"></a>

### X-008 — A practical derivation route and its limits

**Source records:** [P-002](#p-002), [P-023](#p-023), [P-035](#p-035), [P-042](#p-042), [C-021](#c-021), [C-051](#c-051), [B-094](#b-094)–[B-101](#b-101), [G-076](#g-076)–[G-089](#g-089).

1. Compute the commutator closure of the proposed generators. Linear bosonic generators close with the identity; quadratic generators close as a symplectic or orthogonal Lie algebra.
2. Choose an ordered ansatz from that algebra, such as raising exponential $\times$ Cartan exponential $\times$ lowering exponential.
3. Differentiate the ansatz in its parameter, commute factors using the adjoint action, and match the generator. For a closed low-dimensional algebra this gives scalar or matrix ODEs, often a Riccati equation.
4. Fix initial conditions at the identity. Check the adjoint action, the vacuum amplitude, and a conserved quantity. Track singular charts and branches separately from a singular physical operator.
5. For cubic or quartic generators, finite quadratic closure generally fails. The books supply special alternatives—number-shift identities for Kerr evolution, Wick/time-ordered expansions, auxiliary-field formulas, and controlled quadratic approximations. A small BCH truncation is not a general exact factorization.

The finite $2\times2$ compact/noncompact factorizations and a two-mode CAR realization were checked independently while assembling this atlas. Those checks test these identities and selected source corrections; they do not constitute an independent proof of every infinite-dimensional theorem in the dossiers.

<a id="x-009"></a>

### X-009 — Explicit multimode creator–number–annihilator factorization

**Source records:** Berezin Ch.II §4.5, printed pp.98–102 / PDF106–110, (4.26)–(4.27), [B-049](#b-049); Ch.II §5.5, Theorem5.2/(5.15), printed p.124 / PDF132, full fermionic kernel [B-060](#b-060), and pp.125–126 / PDF133–134, (5.16) and its vacuum/phase discussion, [B-061](#b-061)–[B-062](#b-062); Perelomov §8.2, printed pp.104–105 / PDF115–116, (8.2.1)–(8.2.10), [P-036](#p-036), and §9.2, printed pp.113–114 / PDF124–125, [P-042](#p-042)–[P-043](#p-043). **Derived common-convention formula:** the complete coefficients below are reconstructed from the Gaussian kernel and the canonical action, rather than presented as a new numbered equation of either book.

Let $\hat{\mathbf d}$ be a column of finitely many annihilators, with $\eta=+1$ for bosons and $\eta=-1$ for fermions. Fix the direction

$$
\hat U\hat{\mathbf d}\hat U^{-1}=u\hat{\mathbf d}+v\hat{\mathbf d}^{\dagger},
\qquad uu^\dagger-\eta vv^\dagger=I_N,\quad uv^T=\eta vu^T.
$$

The unhatted $u,v$ are numerical matrices, and $\hat{\mathbf d}^{\dagger}$ is a **column** of creators in this formula. On a chart where $u$ is invertible, put

$$
Z=u^{-1}v,\qquad W=\bar v\,u^{-1},\qquad e^h=u^{-1}.
$$

Then

$$
\boxed{\displaystyle
\hat U=c_0\,
\exp\!\left[-\tfrac12\hat{\mathbf d}^{\dagger T}Z\hat{\mathbf d}^{\dagger}\right]
\exp\!\left[\hat{\mathbf d}^{\dagger T}h\hat{\mathbf d}\right]
\exp\!\left[\tfrac\eta2\hat{\mathbf d}^{T}W\hat{\mathbf d}\right],
\qquad c_0=\langle0|\hat U|0\rangle,
\quad |c_0|=\det(uu^\dagger)^{-\eta/4}.}
$$

Here $h$ is a matrix logarithm of $u^{-1}$; the middle factor is its second quantization, so its action is fixed by $e^h$. For a unitary implementer, only the modulus of $c_0$ follows from $u,v$. Its phase must be supplied by the original evolution, a vacuum matrix element, or a continuous lift; [B-071](#b-071) gives the bosonic dynamical phase with the audit correction. In particular, choosing a positive $c_0$ independently at each time can discard a physical phase.

**Derivation and checks.** For bosons the coherent kernel of the product is

$$
K_{\hat U}(x,y)=c_0\exp\!\left[-\tfrac12x^TZx+x^Tu^{-1}y+\tfrac12y^TWy\right],
$$

which is the homogeneous kernel in B-049. For fermions the last pair factor has the opposite sign; commuting the factors with $\hat{\mathbf d},\hat{\mathbf d}^{\dagger}$ using CAR gives the specified canonical matrix, and the vacuum amplitude fixes the remaining scalar. Independent checks use a complex coupled two-mode bosonic quadratic with increasing cutoffs (GEN-22), and exact finite CAR matrix exponentials for genuinely complex, noncommuting two-, three- and four-mode examples, with residuals below $1.6\times10^{-15}$.

Bosonic factor exponentials need a suitable analytic domain or their Gaussian-kernel interpretation; the individual raising and lowering factors need not be bounded operators on the whole Fock space. Bosonic $u$ is always invertible for a unitary finite-mode canonical transformation. Fermionic $u$ can be singular, so this formula describes the nonzero-vacuum chart of an even implementer. A chart failure calls for a different reference vacuum or factorization; the original unitary can remain regular. Affine bosonic displacements can be separated using X-001 and the full affine kernel/normalization in B-048–B-049.


---

<a id="book-b"></a>

## B — Berezin — The Method of Second Quantization


F. A. Berezin, *The Method of Second Quantization*, English translation by Nobumichi Mugibayashi and Alan Jeffrey, Academic Press, copyright 1966. This dossier is based exclusively on the supplied local book, not on any software implementation or other book.

There are 236 actual PDF pages; a trailing form feed in the text extraction produces an empty 237th text segment, not another book page.

### Reading and citation audit

The complete text was read sequentially: all front matter, Introduction, Chapters I–IV, the complete Wick appendix, bibliography, and index. This was a full reading, not an inference from searches. The bundled PDF-to-Markdown converter produced a complete 236-page extraction; page-delimited `pdftotext -layout` was also read, particularly for the opening material and for comparison where Markdown damaged formula layout. Important mathematics was then checked against rendered page images.

Arabic printed page **p** corresponds to one-based PDF page **p+8**, throughout this particular source. Page headers in the sequential extraction and the visually checked pages support this mapping; no missing Arabic pages were found. The eight opening PDF pages comprise publisher lists, title/copyright, and prefaces; the Roman pagination is not consecutive in the physical source, so no fictitious missing front-matter pages are counted as read.

| Material | Printed pages | One-based PDF pages | Coverage |
|---|---:|---:|---|
| Publisher/title/copyright/prefaces | unnumbered; v, vii–ix | 1–8 | complete |
| Introduction | 1–8 | 9–16 | complete |
| Chapter I, §§1–3 | 9–86 | 17–94 | complete |
| Chapter II, §§4–5 | 87–132 | 95–140 | complete |
| Chapter III, §§6–8 | 133–180 | 141–188 | complete |
| Chapter IV, §9 | 181–202 | 189–210 | complete |
| Appendix, Wick’s theorem, §§1–7 | 203–222 | 211–230 | complete |
| References | 223–226 | 231–234 | complete |
| Index | 227–228 | 235–236 | complete |

**Visually inspected printed pages:** 11, 25, 27, 32, 42, 48, 49, 55, 56, 57, 75, 76, 85, 90, 92, 95, 96, 98, 107, 115, 124, 125, 127, 132, 134, 135, 143, 145, 147, 148, 156, 157, 158, 165, 166, 167, 168, 173, 177, 190, 209, 212, 214, 215, 216, 217, 219, 220, 222. These are 49 selected pages, not a claim that all equations in the book were visually transcribed.

**Unread pages:** none. **Unreadable pages:** none. Numerous displayed formulas lose bars, primes, arrows, or line order in machine extraction; the crucial transcribed formulas below were checked visually where listed. Formula locations whose lengthy displays are only described, rather than reproduced, are explicitly identified. No independent proof audit of every theorem is claimed. The Appendix itself declares its analysis formal and omits analytical foundations (printed p.203/PDF211). Three visible source problems are recorded under B-075, B-077, and B-081.

### Conventions and scope

The reference transcription uses normalized hats for operators; source-page images retain the book’s original notation. A hat distinguishes an operator from its scalar or Grassmann symbol, including when the operator acts on the one-particle space rather than Fock space.

| Notation | Role |
|---|---|
| $\hat a,\hat a^*,\hat b,\hat H,\hat N,\hat U,\hat A,\hat P_0,\hat I$ | Operators where these act on Fock space; $\hat I$ has the space appropriate to its formula. |
| $\hat\Phi,\hat\Psi,\hat C,\hat Z,\hat W,\widehat{\mathcal L}$ and hatted Gaussian blocks | One-particle or doubled-space operators in the infinite-dimensional theorems; finite matrices of scalar coefficients are plain when explicitly presented as finite matrix specializations. |
| $a,a^*,\alpha,\alpha^*,f,f^*,z,z^*$; $K_{\hat A},N_{\hat A}$ | Commuting bosonic or anticommuting Grassmann variables and sources; scalar/Grassmann generating kernels and normal symbols of $\hat A$. Grassmann variables are anticommuting scalars, distinct from CAR operators. |
| $A(x,y),B(x,y),C(x,y),\Delta(1\mid2),R(1\mid2)$; $\widehat{\mathcal D}$ | Plain coefficient/contraction kernels; their integral maps are hatted when used as operators, e.g. $\hat A,\hat\Delta,\hat R$. The contraction map $\widehat{\mathcal D}$ acts on the symbol algebra. Explicit derivatives, $T$-ordering, $\operatorname{Tr},\det,\exp$ retain standard notation. |

- To reduce typographical ambiguity, this dossier writes the book’s scalar matrix generating functional $\widetilde A$ as $K_{\hat A}$, and its normal generating functional as $N_{\hat A}$. Neither is an anti-normal or Weyl symbol. State vectors/functionals such as $\Phi,F,P_f,Q_f$ remain plain; the same letters with hats elsewhere denote operators.
- For an operator, $\hat A^*$ is the Hilbert adjoint; $\bar{\hat A}$ is conjugation by the fixed involution; $\hat A'=\bar{\hat A}^*$ is the transpose. Stars, bars and primes retain these distinctions, with the corresponding plain notation for coefficient kernels. The bilinear pairing $fg=(f,g^*)$ explains why the smeared CCR use $(f,g^*)$, rather than the more common sesquilinear notation.
- $a^*a=\int a^*(x)a(x)\,dx$. In finite bosonic dimension, $\Pi da^*da$ means the normalized complex measure $\pi^{-n}\prod_j d\Re a_j\,d\Im a_j$. Infinite-dimensional “continual integrals” are the particular projection limits defined by the book, not translation-invariant Lebesgue measures on a Hilbert space.
- **Source** means a source formula or theorem, with harmless typographical renaming explained above. **Re-expression** means the same statement in new notation. **Derived** means a consequence supplied here; its formula is not attributed as a displayed source equation. Original equation numbers are always chapter/section numbers as printed; Appendix numbers restart at (1).

### I. Fock spaces, symbols, and operator functions

<a id="b-001"></a>

#### B-001 — State spaces and occupation basis

**Location:** Introduction §1–2, printed pp.1–6; PDF9–14; (0.1), (0.10′). **Source/re-expression.** The bosonic/fermionic state spaces are $\mathcal F_\pm(L)=\bigoplus_{n\ge0}L^{\otimes_\pm n}$, with a one-dimensional vacuum sector and dense finite-particle vectors. For $N<\infty$, the fermionic dimension is $2^N$, while bosonic Fock space remains infinite-dimensional. Occupation monomials are $\prod_j(a_j^*)^{n_j}/\sqrt{n_j!}$; CAR force $n_j\in\{0,1\}$. This finite-versus-infinite distinction underlies every later implementability assertion.

<a id="b-002"></a>

#### B-002 — Involution, transpose, and regularized determinants

**Location:** Introduction §3, printed pp.6–8; PDF14–16; (0.12)–(0.13), Theorem 0.1. **Source.** The book defines right operation by $g\hat A=\hat A'g$, so $(g\hat A)f=g(\hat Af)$. Doubled-space determinants use compressions by $\operatorname{diag}(\hat P_n,\bar{\hat P}_n)$. The limit must exist independently of the exhausting projection sequence. It agrees with the Fredholm determinant when $\widehat{\mathcal A}-\hat I$ is trace class; it also exists for blocks with trace-class diagonal deviations and Hilbert–Schmidt off-diagonal blocks. Theorem 0.1 makes a self-adjoint one-particle operator real by unitary conjugation; it does not make all operators simultaneously real.

<a id="b-003"></a>

#### B-003 — Smeared CCR and CAR

**Location:** Ch.I §1.2, printed pp.10–12; PDF18–20; (1.2)–(1.4′). **Source; p.11 visually checked.** Bosons obey $[\hat a(f),\hat a^*(g)]=(f,g^*)\hat I$, with like-type commutators zero. Fermions obey $\{\hat a(f),\hat a^*(g)\}=(f,g^*)\hat I$, with like-type anticommutators zero. Unsmeared relations replace the pairing by a delta distribution. The source warns that $\hat a^*(x)$ at a point of zero measure is generally not an operator mapping a nonzero vector into Fock space; it is an operator-valued distribution/bilinear-form notation requiring smearing.

<a id="b-004"></a>

#### B-004 — Vacuum reconstruction of arbitrary vectors

**Location:** Ch.I §1.3, printed pp.12–13; PDF20–21; (1.5). **Source.** The common annihilation kernel is one-dimensional, spanned by $\Omega$, and
$$
\Phi=\sum_{n\ge0}\frac1{\sqrt{n!}}\int K_n(x_1,\ldots,x_n)\hat a^*(x_1)\cdots\hat a^*(x_n)\,d^nx\,\Omega.
$$
The $K_n$ are symmetric or antisymmetric. No nonzero bosonic vector is killed by every creator. A fermionic vector killed by every creator can exist only for finite-dimensional $L$, reflecting the filled state.

<a id="b-005"></a>

#### B-005 — Bounded fermionic fields

**Location:** Ch.I §1.4, printed pp.13–14; PDF21–22; Theorem 1.1, Corollary 1.1, (1.6)–(1.7). **Source/re-expression.** CAR imply $\|\hat a(f)\|,\|\hat a^*(f)\|\le\|f\|$. For normalized $f$, the corresponding occupation and vacancy operators are complementary orthogonal projections. Creation and annihilation are adjoints under the book’s involution convention; their real linear combinations are self-adjoint bounded operators. This boundedness is specific to fermions and simplifies later evolution proofs.

<a id="b-006"></a>

#### B-006 — Bosonic quadratures and their domains

**Location:** Ch.I §1.4, printed pp.14–16; PDF22–24; Theorem 1.2, (1.8)–(1.10). **Source.** $\hat Q(f)=(\hat a(f)+\hat a^*(f^*))/\sqrt2$ and $\hat P(f)=(\hat a(f)-\hat a^*(f^*))/(i\sqrt2)$ are essentially self-adjoint on finite-particle vectors. The proof reduces to a Jacobi recurrence and the divergent Carleman sum $\sum_n n^{-1/2}$. Deficiency indices are ((0,0)). Their exponentials are therefore actual unitary operators; an algebraic exponential alone would not establish that fact.

<a id="b-007"></a>

#### B-007 — Number operator and a common field domain

**Location:** Ch.I §1.4, printed p.16; PDF24; (1.11)–(1.12). **Source/re-expression.** $\hat N=\int\hat a^*(x)\hat a(x)\,dx$ acts by $n$ on the $n$-particle sector. The form domain $D(\hat N^{1/2})=\{\Phi:\sum_n n\|\Phi_n\|^2<\infty\}$ is contained in the common domain of all bosonic smeared annihilators and creators. The author explicitly leaves open here whether this is the maximal common domain; do not convert that observation into a maximality theorem.

<a id="b-008"></a>

#### B-008 — Commutant and irreducibility

**Location:** Ch.I §1.5, printed pp.16–19; PDF24–27; Theorem 1.3, Corollaries 1.2–1.4, (1.13). **Source.** A closed operator commuting with every creator and annihilator and defined on the vacuum is scalar. Every bounded commuting operator is scalar. The polynomial field algebra has all bounded operators in its weak closure. Bosonic invariant-subspace statements include an essential-self-adjointness condition on restricted quadratures; fermionic boundedness makes that extra condition automatic. These domain qualifications matter for applying Schur-type arguments to unbounded objects.

<a id="b-009"></a>

#### B-009 — Uniqueness of the vacuum representation

**Location:** Ch.I §1.6, printed pp.19–21; PDF27–29; Theorem 1.4, (1.14)–(1.15). **Source.** An irreducible adjoint CCR/CAR family with a common vacuum is unitarily equivalent to the appropriate Fock representation. The book then states finite-degree-of-freedom uniqueness, citing von Neumann for bosons and Jordan–Wigner for fermions. Infinite-dimensional representations without a vacuum need not be equivalent. The source’s finite-dimensional bosonic formulation should be read with its regularity/domain framework, not as uniqueness for arbitrary formal commutator representations.

<a id="b-010"></a>

#### B-010 — What normal form means analytically

**Location:** Ch.I §1.7, printed pp.21–24; PDF29–32; (1.16B), (1.16F), (1.16). **Source.** A normal form is a series with all creators left of all annihilators, weakly convergent on a domain whose intersection with each particle sector is dense. Kernel symmetrization/antisymmetrization makes the coefficient kernels unique. The definition uses limits of regularized kernels and includes the vacuum in the operator domain. Rearranging a formal polynomial does not by itself show the resulting series represents an operator.

<a id="b-011"></a>

#### B-011 — Bilinear forms are weaker than operators

**Location:** Ch.I §1.7, printed pp.22–24; PDF30–32; discussion following (1.16F). **Source.** A distributional kernel can define a densely defined bilinear form even if it cannot define an operator with a nonzero domain. Weak convergence of forms is defined by matrix elements on a dense finite-particle set. This distinction recurs in singular quadratic Hamiltonians and the Thirring chapter; it must be retained when applying algebraic identities to continuum kernels.

<a id="b-012"></a>

#### B-012 — Bounded operators admit normal form

**Location:** Ch.I §1.10, printed pp.28–32; PDF36–40; Theorem 1.5, (1.29)–(1.37). **Source.** Every bounded operator on either Fock space is representable in the specified normal form. The proof decomposes it into particle-sector matrix blocks and inserts the vacuum projection between strings of creators and annihilators. It establishes weak convergence on finite vectors, not a universal operator-norm-convergent Taylor series.

<a id="b-013"></a>

#### B-013 — Normal symbol of the vacuum projector

**Location:** Ch.I §1.10, printed pp.28–30; PDF36–38; (1.30), (1.33)–(1.37). **Source/re-expression.** $K_{\hat P_0}=1$ while $N_{\hat P_0}=e^{-a^*a}$. Equivalently,
$$
\hat P_0=:\!e^{-\hat a^*\hat a}\!:
$$
in both statistics, with the fermionic variables and operator strings ordered as the source prescribes. Its action on a finite vector reduces to a finite sum. This is a normal-ordered exponential, not the ordinary operator $e^{-\hat N}$.

<a id="b-014"></a>

#### B-014 — Matrix kernel versus normal symbol

**Location:** Ch.I §1.10, printed p.32; PDF40; Theorem 1.6, (1.38). **Source; visually checked.**
$$
N_{\hat A}(a^*,a)=K_{\hat A}(a^*,a)e^{-a^*a},\qquad K_{\hat A}=N_{\hat A}e^{a^*a}.
$$
This relation holds for both bosons and fermions under the preceding representability assumptions. It is not an ordinary exponential map between operators; it changes two generating functionals for the same operator.

<a id="b-015"></a>

#### B-015 — Bosonic differential realization

**Location:** Ch.I §1.8, printed pp.24–26; PDF32–34; (1.17)–(1.22). **Source; p.25 visually checked.** On state functionals, $\hat a\leftrightarrow\delta/\delta a^*$, $\hat a^*\leftrightarrow a^*$. On normal symbols, left multiplication by $\hat a$ gives $(a+\delta_{a^*})N_{\hat A}$, right multiplication by $\hat a^*$ gives $(a^*+\delta_a)N_{\hat A}$, and the other two operations are multiplication. On $K_{\hat A}$, the corresponding annihilation/creation actions use pure derivatives on the relevant side. Smearing remains implicit.

<a id="b-016"></a>

#### B-016 — Fermionic derivative sides encode signs

**Location:** Ch.I §1.9, printed pp.26–27; PDF34–35; (1.23)–(1.28). **Source; p.27 visually checked.** State functionals again realize annihilation as a left derivative and creation as left multiplication. Normal symbols use a **left** derivative with respect to $a^*$ for left annihilation, and a **right** derivative with respect to $a$ for right creation. Thus $N_{\hat A\hat a}=N_{\hat A}a$, whereas $N_{\hat a^*\hat A}=a^*N_{\hat A}$. Replacing these by commuting-variable derivatives loses the parity signs.

<a id="b-017"></a>

#### B-017 — Coherent states under the name “Poisson vectors”

**Location:** Ch.I §1.11, printed p.32; PDF40; (1.39), (1.39′), (1.40). **Source; visually checked.**
$$
P_f(a^*)=e^{a^*f-\|f\|^2/2},\qquad Q_f(a^*)=e^{a^*f},\quad \|P_f\|=1,\quad\|Q_f\|^2=e^{\|f\|^2}.
$$
The one-mode particle distribution is $e^{-|f|^2}|f|^{2n}/n!$. These are standard normalized/unnormalized bosonic coherent vectors, but “coherent state” is not the book’s terminology here. Their existence requires $f\in L$.

<a id="b-018"></a>

#### B-018 — Entire analyticity and bounded-set convergence

**Location:** Ch.I §1.11, printed pp.33–34; PDF41–42; Theorems 1.7–1.8, (1.41)–(1.45). **Source.** Vector generating series converge absolutely and uniformly on every fixed ball in $L$. For bounded operators, both normal and matrix generating series converge absolutely and uniformly on bounded sets of two independent arguments $a^*,a\in L$. For unbounded operators the author only describes analogous statements on the set of coherent vectors belonging to the domain; the whole Hilbert-space conclusion cannot be carried over automatically.

<a id="b-019"></a>

#### B-019 — Symbols as coherent matrix elements

**Location:** Ch.I §1.11, printed pp.33–35; PDF41–43; (1.42), (1.46). **Source/re-expression.** $K_{\hat A}$ is the matrix element between unnormalized coherent vectors, and $N_{\hat A}(a^*,a)=\langle P_a,\hat A P_a\rangle$ on conjugate arguments (translated to the usual bracket convention). Thus the diagonal normal symbol is the coherent expectation value. This identity gives the functional meaning of the symbols and proves their convergence for bounded $\hat A$.

<a id="b-020"></a>

#### B-020 — Rank-one projection symbols

**Location:** Ch.I §1.12, printed p.35; PDF43; (1.47)–(1.48). **Source/re-expression.** For the normalized vector functional $\Phi(a^*)$,
$$
K_{|\Phi\rangle\langle\Phi|}=\Phi(a^*)\Phi^*(a),\qquad
N_{|\Phi\rangle\langle\Phi|}=\Phi(a^*)\Phi^*(a)e^{-a^*a}.
$$
The source reverses factor order under involution in the Grassmann case. The same displayed relation works in both statistics only with that order respected.

<a id="b-021"></a>

#### B-021 — Continual integrals are projection-independent limits

**Location:** Ch.I §2.1–3, printed pp.36–41; PDF44–49; (2.1)–(2.13), Theorems 2.1–2.2. **Source.** The complex integral is defined through exhausting complex finite-dimensional subspaces; the real-Hilbert-space construction is distinct and may impose stronger convergence requirements. Integration by parts requires decay at infinity and the existence of the indicated integrals. The notation is not permission to exchange arbitrary infinite sums, traces, or integrals.

<a id="b-022"></a>

#### B-022 — Real Gaussian integral

**Location:** Ch.I §2.4, printed p.41; PDF49; (2.14). **Source.** With the normalized measure,
$$
\int e^{-(\hat Af,f)}\Pi df=(\det \hat A)^{-1/2},\qquad \hat A>0.
$$
The source states existence precisely when $\hat A-\hat I$ is trace class. Finite-dimensional Gaussian formulas therefore extend only under an appropriate determinant-class hypothesis; an arbitrary positive operator on an infinite-dimensional Hilbert space is insufficient.

<a id="b-023"></a>

#### B-023 — Doubled complex Gaussian with sources

**Location:** Ch.I §2.4, printed p.42; PDF50; (2.15)–(2.17). **Source; visually checked.** For transpose-symmetric $\widehat{\mathcal A}=(\hat A_{ij})$,
$$
\int e^{-\frac12(f\ f^*)\widehat{\mathcal A}\binom f{f^*}+(\varphi_1\ \varphi_2)\binom f{f^*}}\Pi df\,df^*
=\det\!\begin{pmatrix}\hat A_{21}&\hat A_{22}\\\hat A_{11}&\hat A_{12}\end{pmatrix}^{-1/2}
e^{\frac12(\varphi_1\ \varphi_2)\widehat{\mathcal A}^{-1}\binom{\varphi_1}{\varphi_2}}.
$$
The off-diagonal deviations $\hat A_{12}-\hat I,\hat A_{21}-\hat I$ are trace class; $\hat A_{11},\hat A_{22}$ are Hilbert–Schmidt; the real part of the quadratic form has a strictly positive lower bound. These are sufficient conditions stated by the source.

<a id="b-024"></a>

#### B-024 — Gaussian moments and bosonic contractions

**Location:** Ch.I §2.4, printed p.43; PDF51; (2.18)–(2.19). **Source/re-expression.** The normalized complex Gaussian kills unequal numbers of holomorphic and antiholomorphic variables. For equal numbers, its moments sum all pairings with positive sign; for a separately symmetric kernel the result is $n!$ times the diagonal integral. This is the elementary contraction identity supporting the later symbol product and Wick theorem.

<a id="b-025"></a>

#### B-025 — Analytic Hilbert-space realization

**Location:** Ch.I §2.5, printed pp.43–46; PDF51–54; Theorem 2.3, (2.20)–(2.27). **Source/re-expression.** The bosonic inner product is the Gaussian integral of the two analytic vector functionals. An analytic functional corresponds to a Fock vector iff its Gaussian squared norm is finite, equivalently $\sum_n\int|K_n|^2<\infty$. The proof justifies sum/integral interchange by finite-dimensional and monotone limits. The remark on p.46 interprets these integrals using a Gaussian measure on an enlarged space, not on $L$ itself.

<a id="b-026"></a>

#### B-026 — Kernel action on a vector

**Location:** Ch.I §2.6, printed pp.46–47; PDF54–55; (2.28)–(2.29). **Source/re-expression.**
$$
(\hat A\Phi)(a^*)=\int K_{\hat A}(a^*,\alpha)\Phi(\alpha^*)e^{-\alpha^*\alpha}\Pi d\alpha^*d\alpha.
$$
The source first proves this for bounded operators. Expressing $K_{\hat A}=N_{\hat A}e^{a^*\alpha}$ yields the corresponding normal-symbol formula. This is the basic functional representation for applying a non-polynomial operator to a state.

<a id="b-027"></a>

#### B-027 — Product of bosonic operators

**Location:** Ch.I §2.7, printed pp.47–48; PDF55–56; (2.30)–(2.31). **Source; p.48 visually checked.**
$$
K_{\hat A\hat B}(a^*,a)=\int K_{\hat A}(a^*,\alpha)K_{\hat B}(\alpha^*,a)e^{-\alpha^*\alpha}\Pi d\alpha^*d\alpha,
$$
$$
N_{\hat A\hat B}(a^*,a)=\int N_{\hat A}(a^*,\alpha)N_{\hat B}(\alpha^*,a)e^{-(\alpha^*-a^*)(\alpha-a)}\Pi d\alpha^*d\alpha.
$$
Initially $\hat A,\hat B$ are bounded. A differential star-product expansion is a **derived** polynomial/analytic reformulation, not an independently asserted convergence theorem here.

<a id="b-028"></a>

#### B-028 — Bosonic trace

**Location:** Ch.I §2.8, printed pp.48–49; PDF56–57; (2.32)–(2.34). **Source; visually checked.** For trace-class $\hat A$,
$$
\operatorname{Tr}\hat A=\int K_{\hat A}(a^*,a)e^{-a^*a}\Pi da^*da=\int N_{\hat A}(a^*,a)\Pi da^*da.
$$
The proof factors a trace-class operator into two Hilbert–Schmidt operators. In infinite dimension, do not use this as a trace prescription for arbitrary bounded Gaussian operators; the trace-class assumption is material.

<a id="b-029"></a>

#### B-029 — Grassmann parity and nilpotence

**Location:** Ch.I §3.1, printed pp.49–50; PDF57–58; (3.1)–(3.3). **Source.** A finite Grassmann algebra with $n$ generators has dimension $2^n$, $x_i x_j=-x_jx_i$, $x_i^2=0$. Even elements commute with all elements; odd elements anticommute with each other. Every element decomposes uniquely into even and odd parts. These are auxiliary Grassmann symbols, not the CAR operators themselves.

<a id="b-030"></a>

#### B-030 — Left/right derivatives and graded Leibniz rule

**Location:** Ch.I §3.1A, printed pp.51–52; PDF59–60; (3.4)–(3.8). **Source/re-expression.** For an ordered monomial, a left derivative removes the target generator after moving it to the left, and a right derivative after moving it to the right. The signs differ by $(-1)^{p-1}$ on degree $p$. Differentiating a product with an odd first factor introduces the graded minus sign. Second derivatives anticommute. These rules, rather than ad hoc sign choices, govern fermionic Gaussian manipulations.

<a id="b-031"></a>

#### B-031 — Grassmann integration and partial-integration warning

**Location:** Ch.I §3.1B, printed pp.52–55; PDF60–63; (3.9)–(3.14′). **Source; p.55 visually checked.** $\int dx_i=0$, $\int x_i\,dx_i=1$, with ordered iterated integration. Differentials anticommute with generators and with one another. The book warns that its right-derivative integration-by-parts identity can fail for integration over only a subset of generators; an even-number-of-unintegrated-generators condition restores it. Carrying a whole-algebra formula into a partial integral without this check can reverse the answer.

<a id="b-032"></a>

#### B-032 — Inverse determinant Jacobian

**Location:** Ch.I §3.1B, printed pp.55–56; PDF63–64; (3.15). **Source; visually checked.** For a finite scalar coefficient matrix $A$, if $x=Ay$,
$$
\int f(x)\,dx_n\cdots dx_1=(\det A)^{-1}\int f(Ay)\,dy_n\cdots dy_1.
$$
The inverse matrix transforms the differentials. This inverse determinant is the fundamental difference from a commuting-variable change of coordinates; the ordering of the differential string must also remain fixed.

<a id="b-033"></a>

#### B-033 — Fermionic Gaussian is a Pfaffian

**Location:** Ch.I §3.1B, printed pp.56–57; PDF64–65; (3.16), footnote 28. **Source/re-expression; visually checked.** With the ordered measure and finite skew coefficient matrix $M$,
$$
\int e^{\frac12x^TMx}\,dx_n\cdots dx_1=\operatorname{Pf}(M),\quad M^T=-M,
$$
and it vanishes for odd $n$. Berezin writes $\det(2a_{ij})^{1/2}$ for an exponent $\sum a_{ij}x_ix_j$. This means the polynomial Pfaffian branch fixed by orientation/continuation, not a principal numerical square root of the determinant. Determinant alone does not retain the sign.

<a id="b-034"></a>

#### B-034 — Clifford operators on a Grassmann algebra

**Location:** Ch.I §3.2, printed pp.57–59; PDF65–67; (3.17)–(3.19). **Source; p.57 visually checked.** Multiplication $\hat X_j$ and left differentiation $\hat D_j$ satisfy $\{\hat X_i,\hat D_j\}=\delta_{ij}\hat I$. Then $\hat Q_j=\hat X_j+\hat D_j$, $\hat P_j=(\hat X_j-\hat D_j)/i$ satisfy $\{\hat Q_i,\hat Q_j\}=\{\hat P_i,\hat P_j\}=2\delta_{ij}\hat I$, $\{\hat P_i,\hat Q_j\}=0$. The source supplies an inner product making multiplication and differentiation adjoint. These are the fermionic counterpart of quadratures, with orthogonal rather than symplectic transformation geometry.

<a id="b-035"></a>

#### B-035 — Infinite Grassmann functionals are not ordinary functions

**Location:** Ch.I §3.3–11, printed pp.59–83; PDF67–91; (3.20)–(3.63). **Source.** The infinite-generator construction specifies graded topological spaces, generalized generators, involution, differentiability after smearing, and a projection-independent continual integral. Generating functionals of bounded operators and Fock vectors lie in a common Grassmann algebra, and a suitable completion is normed. This justifies formal-looking series in that framework; it does not identify them with complex-valued analytic functions or normalized bosonic Poisson vectors.

<a id="b-036"></a>

#### B-036 — Doubled fermionic Gaussian

**Location:** Ch.I §3.10, printed pp.75–76; PDF83–84; (3.45)–(3.49). **Corrected transcription; original determinant visually rechecked.** For block-skew-transpose $\widehat{\mathcal A}=(\hat A_{ij})$,
$$
\int e^{\frac12(\alpha\ \alpha^*)\widehat{\mathcal A}\binom\alpha{\alpha^*}}\Pi d\alpha^*d\alpha
=\det\!\begin{pmatrix}\hat A_{11}&\hat A_{12}\\\hat A_{21}&\hat A_{22}\end{pmatrix}^{1/2}.
$$
The source uses the **original block matrix**, not the column-permuted matrix in the earlier transcription. One mode $\mathcal A=\left(\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right)$ gives normalized integral 1 and determinant 1, while the old permuted determinant was $-1$; a square-root branch cannot repair that error. The square root is the continued polynomial/Pfaffian branch described in B-033. In finite block coordinate order $(\alpha_1,\ldots,\alpha_n,\alpha_1^*,\ldots,\alpha_n^*)$, with normalized paired integration $\int e^{-\alpha^*\alpha}\Pi d\alpha^*d\alpha=1$, it is explicitly
$$
(-1)^{n(n-1)/2}\operatorname{Pf}(\mathcal A),
\qquad [\operatorname{Pf}(\mathcal A)]^2=\det\mathcal A.
$$
Here $\mathcal A$ is the finite numerical specialization of the hatted one-particle block operator. The infinite-dimensional statement uses the source's generalized determinant and hypotheses, not an arbitrary determinant of an infinite array. The Gaussian has determinant power $+1/2$, versus $-1/2$ for commuting variables.

<a id="b-037"></a>

#### B-037 — Signed fermionic Gaussian moments

**Location:** Ch.I §3.10, printed pp.76–77; PDF84–85; (3.50), (3.50′), (3.51). **Source; p.76 visually checked.** Equal creation/annihilation degrees give a signed sum over permutations of delta pairings; unequal degrees vanish. The sign is the permutation parity. For separately antisymmetric kernels the integrated formula has the same $n!$-times-diagonal shape as the bosonic identity. Reversing the creator/annihilator string changes signs unless the measure and kernel ordering are adjusted together.

<a id="b-038"></a>

#### B-038 — Fermionic Fock inner product

**Location:** Ch.I §3.12, printed pp.83–84; PDF91–92; Theorem 3.1, (3.64)–(3.66′). **Source/re-expression.** Vector inner products have a Gaussian Grassmann-integral expression, and finite Gaussian squared norm characterizes the Fock functionals. The two source forms differ by ordering $\Phi_1\Phi_2^*$ versus $\Phi_2^*\Phi_1$, together with $e^{-aa^*}\Pi da\,da^*$ versus $e^{-a^*a}\Pi da^*da$. They are equivalent as complete ordered expressions, not termwise after arbitrary commuting-variable simplifications.

<a id="b-039"></a>

#### B-039 — Fermionic kernel composition

**Location:** Ch.I §3.13–14, printed pp.84–85; PDF92–93; (3.67)–(3.70). **Source/re-expression.** The vector-action and kernel-composition formulas have the same displayed form as B-026–B-027, with Grassmann variables and the source’s fixed factor/differential order. In particular (3.70) is the normal-symbol product with Gaussian $e^{-(\alpha^*-a^*)(\alpha-a)}$. These formulas give a direct route to products of exponentials of quadratic operators through Gaussian integration.

<a id="b-040"></a>

#### B-040 — Fermionic trace differs from bosonic trace

**Location:** Ch.I §3.15, printed p.85; PDF93; (3.71)–(3.72). **Source; visually checked.**
$$
\operatorname{Tr}\hat A=\int K_{\hat A}(a^*,a)e^{-aa^*}\Pi da\,da^*
=\int N_{\hat A}(a^*,a)e^{-2aa^*}\Pi da\,da^*.
$$
This is the ordinary trace of a trace-class operator. The doubled Gaussian factor in the normal-symbol form must not be dropped by copying the bosonic formula.

<a id="b-041"></a>

#### B-041 — Normalized tracial state and infinite-mode factors

**Location:** Ch.I §3.16, printed pp.85–86; PDF93–94; (3.73)–(3.77). **Source; p.85 visually checked.** The alternative trace evaluates $N_{\hat A}(a^*/\sqrt2,a/\sqrt2)$ under $e^{-aa^*}\Pi da\,da^*$. For $n$ modes it is $2^{-n}\operatorname{Tr}\hat A$. In infinite dimension the two traces have different operator classes; this is not multiplication by an infinite normalization constant. The normalized trace extends to the specified polynomial/Majorana rings and their type- II₁ factor closures; it takes value 1 on the identity.

### II. Linear canonical transformations and Gaussian implementers

<a id="b-042"></a>

#### B-042 — Bosonic affine canonical transformations

**Location:** Ch.II §4.1, printed pp.87–89; PDF95–97; (4.1)–(4.4). **Source/re-expression.**
$$
\hat b=\hat\Phi\hat a+\hat\Psi\hat a^*+f\hat I,\qquad
\hat b^*=\bar{\hat\Psi}\hat a+\bar{\hat\Phi}\hat a^*+f^*\hat I.
$$
The book begins with operators on specified dense involution-invariant one-particle domains and permits distributional $f$. “Canonical” requires reversibility and preservation of CCR. “Proper” means there exists a unitary $\hat U$ with $\hat b=\hat U\hat a\hat U^{-1}$. In finite dimension all such transformations are proper; in infinite dimension this additional requirement is substantial.

<a id="b-043"></a>

#### B-043 — Bosonic block relations and inverse

**Location:** Ch.II §4.2, printed pp.90–91; PDF98–99; (4.5)–(4.8). **Source; p.90 visually checked.** With
$$
\hat S=\begin{pmatrix}\hat\Phi&\hat\Psi\\\bar{\hat\Psi}&\bar{\hat\Phi}\end{pmatrix},\quad
\hat S^{-1}=\begin{pmatrix}\hat\Phi^*&-\hat\Psi'\\-\hat\Psi^*&\hat\Phi'\end{pmatrix},
$$
the CCR conditions include $\hat\Phi\hat\Psi'=\hat\Psi\hat\Phi'$, $\hat\Phi\hat\Phi^*-\hat\Psi\hat\Psi^*=\hat I$, and the corresponding inverse relations $\hat\Phi^*\hat\Phi-\hat\Psi'\bar{\hat\Psi}=\hat I$, $\hat\Phi^*\hat\Psi=\hat\Psi'\bar{\hat\Phi}$. For unbounded one-particle coefficients, the source interprets these as identities of Hermitian forms on the stated common domains.

<a id="b-044"></a>

#### B-044 — Symplectic geometry and group law

**Location:** Ch.II §4.2, printed pp.90–92; PDF98–100; (4.9)–(4.13). **Source/re-expression.** Composition of homogeneous canonical transformations multiplies their doubled matrices. After changing from $(\hat a,\hat a^*)$ to real quadratures, the finite-dimensional coefficient-matrix specialization is real symplectic: $M^TJM=J$, with inverse $M^{-1}=\bigl(\begin{smallmatrix}D^T&-B^T\\-C^T&A^T\end{smallmatrix}\bigr)$. The source also factors a transformation into a passive unitary one-particle change of basis and a transformation whose $\hat\Phi$ is Hermitian positive.

<a id="b-045"></a>

#### B-045 — Bosonic invertibility and squeezing contraction

**Location:** Ch.II §4.2–3, printed pp.91, 94–95, 101; PDF99, 102–103, 109; (4.19), discussion following (4.23). **Source/re-expression; p.95 visually checked.** $\hat\Phi$ has a bounded inverse. For a proper transformation, $\hat Z=\hat\Phi^{-1}\hat\Psi$ is symmetric, Hilbert–Schmidt, and $\|\hat Z\|<1$; specifically $\hat Z\hat Z^*=\hat I-(\hat\Phi^*\hat\Phi)^{-1}$. The corresponding right-chart contraction is $\hat W=\bar{\hat\Psi}\hat\Phi^{-1}$, with $\hat W^*\hat W=\hat I-(\hat\Phi\hat\Phi^*)^{-1}$ and $\|\hat W\|<1$; compare P-034 and the visually checked Perelomov (8.1.8)–(8.1.9), printed p.102/PDF113. The conjugation bar is essential: $\|\hat\Psi\hat\Phi^{-1}\|<1$ is false for general complex multimode canonical blocks. These are the normalizability conditions for the squeezed Gaussian vacuum, rather than arbitrary restrictions added to a formal Gaussian ansatz.

<a id="b-046"></a>

#### B-046 — Bosonic implementability criterion

**Location:** Ch.II §4.3, printed pp.92–96; PDF100–104; Theorem 4.1, (4.14)–(4.25). **Source; theorem visually checked.** The affine canonical transformation is proper **iff** $\hat\Psi$ is Hilbert–Schmidt and $f\in L$. Thus pair production must have finite Hilbert–Schmidt norm and the displacement finite one-particle norm. Under the criterion $\hat\Phi$ is bounded and $\hat\Phi\hat\Phi^*-\hat I=\hat\Psi\hat\Psi^*$ is trace class, so the determinant used in normalization exists. This is not a finite-volume assertion blindly extended to a continuum.

<a id="b-047"></a>

#### B-047 — Vacuum criterion for an implementer

**Location:** Ch.II §4.3, printed pp.92–93; PDF100–101; Lemma 4.1, (4.15)–(4.16). **Source; p.92 visually checked.** A canonical transformation is proper iff there is a nonzero Fock vector annihilated by every transformed annihilator. Applying the transformed creators to this vacuum gives an equivalent irreducible Fock representation. A vacuum-defined intertwiner differs from a unitary only by a scalar; this follows from the commutant theorem and uses its domain assumptions.

<a id="b-048"></a>

#### B-048 — Bosonic transformed vacuum

**Location:** Ch.II §4.3, printed pp.95–96; PDF103–104; (4.20)–(4.25). **Source; visually checked.**
$$
F_0(a^*)=c\exp\!\left[-\tfrac12a^*\hat\Phi^{-1}\hat\Psi a^*-a^*\hat\Phi^{-1}f\right].
$$
For $f=0$, $c=\theta\det(\hat\Phi\hat\Phi^*)^{-1/4}$, $|\theta|=1$. The full affine normalization (4.25), now transcribed explicitly, is
$$
c=\frac{\theta}{\det(\hat\Phi\hat\Phi^*)^{1/4}}
\exp\!\left[\frac14(f\ f^*)
\begin{pmatrix}\hat\Phi'^{-1}\hat\Psi^*&-\hat I\\-\hat I&\hat\Phi^{*-1}\hat\Psi'\end{pmatrix}
\binom f{f^*}\right].
$$
Here the primes and stars have the source meanings (transpose and adjoint), not scalar conjugation. **Independent check:** rewriting the squared Bargmann Gaussian in real coordinates gives the same positive modulus for arbitrary finite complex canonical blocks; the remaining $\theta$ is the arbitrary unit phase. See audit check GEN-07. The vacuum satisfies the first-order functional equation (4.20), and its squared norm is evaluated by the complex Gaussian formula B-023.



Source page: Berezin, printed p.96, one-based PDF p.104.


<a id="b-049"></a>

#### B-049 — Full bosonic Gaussian kernel

**Location:** Ch.II §4.5, printed pp.98–102; PDF106–110; Theorem 4.3, (4.26), (4.26′), (4.27). **Source; p.98 visually checked.**
$$
K_{\hat U}(a^*,a)=c\exp\!\left[
\tfrac12(a\ a^*)
\begin{pmatrix}\bar{\hat\Psi}\hat\Phi^{-1}&\hat\Phi'^{-1}\\\hat\Phi^{-1}&-\hat\Phi^{-1}\hat\Psi\end{pmatrix}
\binom a{a^*}
+a(f^*-\bar{\hat\Psi}\hat\Phi^{-1}f)-a^*\hat\Phi^{-1}f\right],
\quad N_{\hat U}=K_{\hat U}e^{-a^*a}.
$$
For the homogeneous case $c=\theta\det(\hat\Phi\hat\Phi^*)^{-1/4}$. The scalar for the affine case is exactly the normalization (4.26′), identical to (4.25). The kernel is uniquely determined up to phase by $\hat b\hat U=\hat U\hat a$, $\hat b^*\hat U=\hat U\hat a^*$.



Source page: Berezin, printed p.98, one-based PDF p.106.


<a id="b-050"></a>

#### B-050 — Displacement operator in the source’s conjugation direction

**Location:** Ch.II §4.1, §4.5, printed pp.88–89, 98–99; PDF96–97, 106–107; specialization of (4.26)–(4.27). **Derived, not a separately numbered source formula.** Setting $\hat\Phi=\hat I,\hat\Psi=0$,
$$
\hat U_f=e^{\hat a(f^*)-\hat a^*(f)},\quad \hat U_f\hat a\hat U_f^{-1}=\hat a+f\hat I,
\qquad N_{\hat U_f}=e^{-\|f\|^2/2+af^*-a^*f}.
$$
This follows from B-049 or the central-commutator exponential identity. It is the inverse-sign convention to the common displacement satisfying $\hat D(f)^*\hat a\hat D(f)=\hat a+f\hat I$. Infinite-mode existence requires $f\in L$, as B-046 states.

<a id="b-051"></a>

#### B-051 — One-mode squeezing as a checkable specialization

**Location:** Ch.II §4.7, printed pp.111–112; PDF119–120; example preceding (4.48), and (4.48). **Source transformation plus derived implementer/kernel.** The source takes $\hat b=\hat a\cosh\sigma+\epsilon\hat a^*\sinh\sigma$, $|\epsilon|=1$, real $\sigma$. In the source conjugation direction,
$$
\hat U_\sigma=\exp\!\left[\tfrac\sigma2(\bar\epsilon\hat a^2-\epsilon\hat a^{*2})\right],
\quad
K_{\hat U_\sigma}=(\cosh\sigma)^{-1/2}
e^{\frac12\bar\epsilon\tanh\sigma\,a^2+\operatorname{sech}\sigma\,a^*a-\frac12\epsilon\tanh\sigma\,a^{*2}}.
$$
The phase is chosen continuously from identity. The source’s (4.48) is instead a Gaussian integral transforming an arbitrary normal symbol, with a specified convergent complex contour.

<a id="b-052"></a>

#### B-052 — Projective phases and the double cover

**Location:** Ch.II §4.4–5, printed pp.97–98, 102–103; PDF105–106, 110–111; Theorem 4.2 and following (4.33). **Source.** For proper transformations with Hermitian positive $\hat\Phi$, the regularized determinant of the doubled symplectic matrix is 1. Implementers form a projective representation. On the subgroup $\hat\Phi-\hat I$ trace class, choosing $c=(\det\hat\Phi)^{-1/2}$ leaves an unavoidable sign ambiguity. The author calls this a spinor representation of the symplectic group; modern terminology is the metaplectic double cover. A square-root branch must be tracked along paths, especially around closed loops.

<a id="b-053"></a>

#### B-053 — Passive transformation of normal symbols

**Location:** Ch.II §4.6, printed p.104; PDF112; (4.34)–(4.35). **Source.** If $\hat\Psi=f=0$, $\hat\Phi$ is unitary, and
$$
N_{\hat U\hat A\hat U^{-1}}(a^*,a)=N_{\hat A}(\bar{\hat\Phi} a^*,\hat\Phi a).
$$
The bar is the source’s one-particle conjugation, not the Hilbert adjoint. This simple substitution rule applies to passive transformations. A mixing transformation generally requires contractions and a Gaussian kernel instead of merely replacing $a,a^*$ in a normal symbol.

<a id="b-054"></a>

#### B-054 — General bosonic transformation kernel for normal symbols

**Location:** Ch.II §4.6, printed pp.104–111; PDF112–119; Lemmas 4.2–4.3, Theorem 4.4, (4.38)–(4.47). **Source; p.107 visually checked.** For bounded $\hat A$, Hermitian $\hat\Phi$, invertible $\hat\Psi$, and finite-dimensional spectral projections $\hat P_n\uparrow \hat I$,
$$
N_{\hat U\hat A\hat U^{-1}}=\lim_n\int\mathcal K_n(a^*,a\mid\alpha^*,\alpha)N_{\hat A}(\alpha^*,\alpha)\,d^nt\,d^n\tau.
$$
Here $\mathcal K_n$ is the Gaussian (4.39), with matrix $\widehat{\mathcal L}=\bigl(\begin{smallmatrix}\bar{\hat\Phi}\bar{\hat\Psi}^{-1}&-\hat I\\-\hat I&\hat\Phi\hat\Psi^{-1}\end{smallmatrix}\bigr)$ from (4.40), and contours selected by Lemma 4.3 to make its exponent’s real part negative definite. The theorem is not a real-axis integral without a contour prescription. The author extends it to polynomial operators, while a still broader extension is described only tentatively on p.111.



Source page: Berezin, printed p.107, one-based PDF p.115.


<a id="b-055"></a>

#### B-055 — Gaussian normal symbols are closed under canonical conjugation

**Location:** Ch.II §4.7, printed pp.112–115; PDF120–123; (4.49)–(4.57). **Source/re-expression; p.115 visually checked.** Conjugating a Gaussian normal symbol produces a Gaussian with a transformed quadratic form and determinant prefactor. For $N_{\hat A}=e^{-\frac12(a\ a^*)\widehat{\mathcal A}\binom a{a^*}}$, the source’s (4.55) has prefactor $\det \hat C^{-1/2}$,
$$
\hat C=\begin{pmatrix}\bar{\hat\Phi}+\hat A_{11}\hat\Psi&(\hat A_{12}-\hat I)\bar{\hat\Psi}\\(\hat A_{21}-\hat I)\hat\Psi&\hat\Phi+\hat A_{22}\bar{\hat\Psi}\end{pmatrix},
\quad b=\hat\Phi a+\hat\Psi a^*+f.
$$
Here $b$ is a commuting symbol variable; $\hat C$ and $\hat A_{ij}$ are one-particle block operators. The full exponent is the explicitly ordered block product in (4.55), not reproduced here. The derivation imposes Gaussian decay (4.49), initially $\hat\Phi=\hat\Phi^*$ and $\hat\Psi^{-1}$; (4.57) handles a passive factor. Limits removing $\hat\Psi^{-1}$ require an admissible approximating family, not arbitrary algebraic substitution.



Source page: Berezin, printed p.115, one-based PDF p.123.


<a id="b-056"></a>

#### B-056 — Fermionic canonical relations

**Location:** Ch.II §5.1–2, printed pp.116–118; PDF124–126; (5.1)–(5.6). **Source/re-expression.** Fermionic transformations are homogeneous:
$$
\hat b=\hat\Phi\hat a+\hat\Psi\hat a^*,\quad
\hat\Phi\hat\Phi^*+\hat\Psi\hat\Psi^*=\hat I,\quad\hat\Phi\hat\Psi'+\hat\Psi\hat\Phi'=0,
$$
$$
\hat S^{-1}=\begin{pmatrix}\hat\Phi^*&\hat\Psi'\\\hat\Psi^*&\hat\Phi'\end{pmatrix}.
$$
The plus signs in the CAR identities replace the bosonic minus signs. Unlike the bosonic chapter, ordinary complex additive shifts are not introduced as affine CAR transformations. Grassmann external sources used in integral calculations are not scalar shifts of Hilbert-space CAR operators.

<a id="b-057"></a>

#### B-057 — Fermionic boundedness and orthogonal geometry

**Location:** Ch.II §5.2, printed pp.117–118; PDF125–126; (5.4)–(5.10). **Source.** The identity $\hat\Phi\hat\Phi^*+\hat\Psi\hat\Psi^*=\hat I$ forces $\hat\Phi,\hat\Psi$ bounded, even for improper transformations. In real Majorana coordinates the doubled canonical matrices become real orthogonal matrices. Composition again multiplies matrices. In finite dimension this is the compact orthogonal group, in contrast to the noncompact symplectic group governing bosons.

<a id="b-058"></a>

#### B-058 — Fermionic implementability, including property A

**Location:** Ch.II §5.1, §5.3, printed pp.116, 119–123; PDF124, 127–131; Theorem 5.1, Lemma 5.1, (5.11). **Source.** The book assumes **property A**, $\hat\Phi=(\hat\Phi\hat\Phi^*)^{1/2}\hat V$ with $\hat V$ unitary. Under this hypothesis, a canonical transformation is proper iff $\hat\Psi$ is Hilbert–Schmidt. Existence of a nonzero transformed vacuum is again equivalent to properness. Do not omit property A when reporting this theorem as the exact source statement; the dossier does not silently replace it by a more general later implementability theorem.

<a id="b-059"></a>

#### B-059 — Singular fermionic charts and filled sectors

**Location:** Ch.II §5.3–4, printed pp.120–124; PDF128–132; Lemmas 5.2–5.5, (5.12)–(5.14). **Source.** For $\hat\Phi=0$, properness forces the one-particle space finite-dimensional. For a proper property-A transformation, $\ker\hat\Phi$ is finite-dimensional. The source splits off this filled sector and treats the complement where $\hat\Phi^{-1}$ exists. Its limiting argument on p.124 uses an **even-dimensional** singular component; Theorem 5.3 also explicitly requires even-dimensional zero eigenspace. Therefore a formula involving $\hat\Phi^{-1}$ is a chart, not a global parametrization of all fermionic transformations or parity components.

<a id="b-060"></a>

#### B-060 — Fermionic Gaussian implementer kernel

**Location:** Ch.II §5.5, printed p.124; PDF132; Theorem 5.2, (5.15). **Source; visually checked.** For a proper transformation with two-sided $\hat\Phi^{-1}$,
$$
K_{\hat U}=\theta\det(\hat\Phi\hat\Phi^*)^{1/4}
\exp\!\left[-\tfrac12\left(a\bar{\hat\Psi}\hat\Phi^{-1}a-2a^*\hat\Phi^{-1}a+a^*\hat\Phi^{-1}\hat\Psi a^*\right)\right],
\qquad N_{\hat U}=K_{\hat U}e^{-a^*a},\quad|\theta|=1.
$$
All variables in the exponent are Grassmann variables with the displayed order. The determinant power is the opposite of the bosonic vacuum normalization power. At singular $\hat\Phi$ the source calls for limits in an admissible sector, not evaluation of an inverse that does not exist.

<a id="b-061"></a>

#### B-061 — Fermionic paired vacuum

**Location:** Ch.II §5.5, printed p.125; PDF133; (5.16). **Source; visually checked.**
$$
F(a^*)=\theta\det(\hat\Phi\hat\Phi^*)^{1/4}
\exp[-\tfrac12a^*\hat\Phi^{-1}\hat\Psi a^*].
$$
The pairing matrix is skew-symmetric. For finite modes the exponential is a finite polynomial and has even fermion parity in this invertible-$\hat\Phi$ chart. A general transformed vacuum can require an occupied singular sector, so this even-pair exponential alone is not the universal fermionic vacuum form.

<a id="b-062"></a>

#### B-062 — Fermionic projective/spin representation

**Location:** Ch.II §5.5, printed pp.125–126; PDF133–134; discussion following (5.16). **Source.** Implementers are fixed only up to phase. On the subgroup $\hat\Phi-\hat I$ trace class, the coefficient can be chosen $\det\hat\Phi^{1/2}$, leaving a double-valued representation. The source identifies it with the spinor representation of the orthogonal group. The square-root sign is not determined by the local canonical action alone; continuity along a chosen lift supplies the missing information.

<a id="b-063"></a>

#### B-063 — Fermionic transformation of normal symbols

**Location:** Ch.II §5.6–7, printed pp.126–131; PDF134–139; Theorem 5.3, (5.17)–(5.30). **Source; p.127 visually checked.** Passive conjugation is again variable substitution. With mixing, (5.19)–(5.20) give a Grassmann Gaussian integral kernel, with $\widehat{\mathcal L}=\bigl(\begin{smallmatrix}\bar{\hat\Phi}\bar{\hat\Psi}^{-1}&-\hat I\\\hat I&-\hat\Phi\hat\Psi^{-1}\end{smallmatrix}\bigr)$ and determinant factor $\det(\hat P_n\hat\Psi\hat\Psi^*\hat P_n)^{1/2}$. The theorem assumes Hermitian $\hat\Phi$, invertible $\hat\Psi$, even-dimensional singular eigenspace, and an exhausting even-dimensional spectral sequence. Section 5.7 derives the kernel directly by intertwining left/right field multiplications in finite dimension.



Source page: Berezin, printed p.126, one-based PDF p.134.


Berezin, printed p.127, one-based PDF p.135.


<a id="b-064"></a>

#### B-064 — Fermionic Gaussian-symbol conjugation and determinant class

**Location:** Ch.II §5.8, printed p.132; PDF140; (5.31)–(5.32). **Source; visually checked.** A Grassmann Gaussian normal symbol remains Gaussian after canonical conjugation; (5.32) contains $\det \hat C^{1/2}$ times an ordered quadratic exponential. Initially $\hat\Phi$ is Hermitian nonnegative and $\hat\Psi^{-1}$ exists. For a general polar factor $\hat\Phi=\hat S\hat V$, an arbitrary unitary $\hat V$ need not have a determinant. The source explicitly replaces $\hat C$ by $\hat C_1$, built from $\hat S=(\hat\Phi\hat\Phi^*)^{1/2}$ and $\hat\Psi_1=\hat\Psi\bar{\hat V}^{-1}$. An unqualified $\det\hat\Phi$ would therefore be invalid in the infinite-mode case.

### III. Quadratic Hamiltonians, exponentials, domains, and normal forms



Source page: Berezin, printed p.132, one-based PDF p.140.


<a id="b-065"></a>

#### B-065 — General quadratic expression

**Location:** Ch.III §6.1, printed pp.133–134; PDF141–142; (6.1). **Source; p.134 visually checked.**
$$
\hat H=\tfrac12\int[A(x,y)\hat a^*(x)\hat a^*(y)+B(x,y)\hat a(x)\hat a(y)+2C(x,y)\hat a^*(x)\hat a(y)]\,dxdy
+\hat a^*(f)+\hat a(g)+k\hat I.
$$
Pair kernels satisfy $A=A',B=B'$ for bosons and $A=-A',B=-B'$ for fermions. Here $A(x,y),B(x,y),C(x,y)$ are coefficient kernels; the corresponding one-particle maps used below are $\hat A,\hat B,\hat C$. The source calls the operator homogeneous when $f=g=k=0$. It distinguishes an expression of this shape from a quadratic **normal form** representing an actual operator on the prescribed domain.

<a id="b-066"></a>

#### B-066 — Vacuum-domain necessity

**Location:** Ch.III §6.1, printed p.134; PDF142; (6.2). **Source; visually checked.** A quadratic normal form must have Hilbert–Schmidt pair-creation kernel $A$ and $f\in L$, because applying it to the vacuum produces independent two-particle and one-particle components of finite norm. This is a direct domain obstruction. Pair-creation kernels proportional to a continuum delta function generally fail it; the formal quadratic degree alone does not guarantee an operator.

<a id="b-067"></a>

#### B-067 — Self-adjointness of quadratic normal forms

**Location:** Ch.III §6.2, printed pp.135–138; PDF143–146; Theorem 6.1, Lemma 6.1, (6.3)–(6.8). **Source; p.135 visually checked.** If $\hat A$ is Hilbert–Schmidt, $\hat C=\hat C^*$, and $f\in L$, the formally self-adjoint quadratic normal form is self-adjoint on its natural domain. On finite vectors in the second-quantized $\hat C$-domain it is essentially self-adjoint. The proof uses a block Jacobi/Carleman criterion and estimates pair-sector growth by $\sqrt{n(n+1)}\|\hat A\|_2$, linear-sector growth by $\sqrt n\|f\|$. Later fermionic Gaussian-evolution formulas restrict to $f=0$; this theorem alone is not a claim that odd fermionic linear terms generate homogeneous Bogoliubov flow.

<a id="b-068"></a>

#### B-068 — Exponentials of bounded perturbations

**Location:** Ch.III §6.3, printed pp.139–141; PDF147–149; Lemma 6.2, (6.9)–(6.17). **Source/re-expression.** For $\hat A=\hat S+\hat T$, $\hat S=\hat S^*$, $\hat T$ bounded, the evolution $e^{it\hat A}$ is constructed as the unique bounded solution of the evolution equation; the interaction-picture series converges in operator norm. Its bound is $\|e^{it\hat A}\|\le e^{|t|\|\hat T\|}$. If $\int_0^t e^{i\tau \hat S}\hat Te^{-i\tau \hat S}d\tau$ is Hilbert–Schmidt with locally integrable Hilbert–Schmidt norm, then $e^{it\hat A}-e^{it\hat S}$ is Hilbert–Schmidt.

<a id="b-069"></a>

#### B-069 — Integrated pairing criterion for implementable evolution

**Location:** Ch.III §6.3, printed pp.141–142; PDF149–150; Lemma 6.3, (6.18). **Source/re-expression.** For self-adjoint $\hat C$, bounded symmetric $\hat A$, and
$$
\hat F(t)=\int_0^t e^{-i\tau \hat C}\hat A e^{-i\tau\bar{\hat C}}\,d\tau,
$$
Hilbert–Schmidt $\hat F(t)$ with locally integrable norm, the doubled matrix exponential with generator $\widehat{\mathcal L}=\bigl(\begin{smallmatrix}-\hat C&-\hat A\\\bar{\hat A}&\bar{\hat C}\end{smallmatrix}\bigr)$ exists and is the matrix of a proper bosonic canonical transformation. The criterion can hold although the instantaneous pairing $\hat A$ is not Hilbert–Schmidt; this is relevant to singular quadratic expressions.

<a id="b-070"></a>

#### B-070 — Determinant-class phase condition

**Location:** Ch.III §6.3, printed pp.142–143; PDF150–151; Lemma 6.4, (6.19). **Source.** If additionally $\hat F(t)\bar{\hat A}$ is trace class with the stated local integrability of its trace norm, then $\hat\Phi(t)e^{it\hat C}-\hat I$ is trace class. This supplies the Fredholm determinant in the Gaussian propagator. The condition is on the identity deviation, not on $\hat\Phi(t)e^{it\hat C}$ itself; at $t=0$ the latter is the identity and is not trace class in infinite dimension.

<a id="b-071"></a>

#### B-071 — Bosonic quadratic propagator, including the scalar phase

**Location:** Ch.III §6.3, printed pp.143–144; PDF151–152; Theorem 6.2, (6.20)–(6.26). **Source; p.143 visually checked.** For the self-adjoint normal form,
$$
e^{it\widehat{\mathcal L}}=\begin{pmatrix}\hat\Phi&\hat\Psi\\\bar{\hat\Psi}&\bar{\hat\Phi}\end{pmatrix},\quad
\binom g{g^*}=i\int_0^t e^{i\tau\widehat{\mathcal L}}\,d\tau\binom{-f}{f^*}.
$$
The kernel is B-049 with $f$ replaced by $g(t)$, and coefficient
$$
c(t)=\det(\hat\Phi e^{it\hat C})^{-1/2}
\exp\!\left[i\int_0^t\left(\tfrac12g\hat\Phi'^{-1}\bar{\hat A}\hat\Phi^{-1}g-f^*\hat\Phi^{-1}g\right)d\tau\right].
$$
**Source correction:** printed (6.23) omits the $\tfrac12$ multiplying the term quadratic in $g$. The vacuum component of $\partial_t\hat U=i\hat H\hat U$ gives
$$
\frac{\dot c}{c}=i\!\left[\tfrac12\operatorname{Tr}(\bar{\hat A}\hat R)+\tfrac12\ell^T\bar{\hat A}\ell+f^*\ell\right],
\qquad \hat R=-\hat\Phi^{-1}\hat\Psi,\quad \ell=-\hat\Phi^{-1}g.
$$
The trace term integrates to the displayed determinant. Thus the quadratic-in-$g$ coefficient must be $1/2$. An independent one-mode Fock calculation with complex forcing, at cutoffs 30, 50 and 80, agrees with this correction within $2.2\times10^{-15}$; the printed coefficient has error $6.17\times10^{-3}$ for the same parameters (GEN-18). This is a source error, not a transcription discrepancy. The square root is fixed continuously from $c(0)=1$. This dynamical phase is more information than the canonical action alone. A constant $k$ in $\hat H$ adds the derived scalar factor $e^{itk}$.



Source page: Berezin, printed p.143, one-based PDF p.151.


<a id="b-072"></a>

#### B-072 — Heisenberg/Hadamard action closes at first degree

**Location:** Ch.III §6.4, printed pp.144–145; PDF152–153; (6.27)–(6.30). **Source/re-expression; p.145 visually checked.** For the book’s $\hat U(t)=e^{it\hat H}$,
$$
\hat U(t)\binom{\hat a}{\hat a^*}\hat U(t)^{-1}
=e^{it\widehat{\mathcal L}}\binom{\hat a}{\hat a^*}+\binom g{g^*}\hat I.
$$
The proof differentiates conjugation, computes the commutators, and solves a linear affine system. **Derived interpretation:** this resums the Hadamard series $e^{it\widehat{\operatorname{ad}}_{\hat H}}$ on the linear field span. The book does not offer a general BCH convergence theorem for arbitrary unbounded operators; its operator hypotheses are those of the preceding quadratic results.

<a id="b-073"></a>

#### B-073 — Log-determinant differentiation for the propagator

**Location:** Ch.III §6.4, printed p.146; PDF154; (6.31)–(6.36). **Source/re-expression.** The scalar coefficient is fixed by the vacuum component of the evolution equation. The source uses
$$
\frac{d}{dt}\log\det \hat M=\operatorname{Tr}(\hat M^{-1}\dot{\hat M})
$$
for the determinant-class family $\hat M=\hat\Phi e^{it\hat C}$, and relates $\hat\Phi^{-1}\dot{\hat\Phi}$ to $\hat C$ and the pairing block. This avoids losing the scalar phase when replacing a quadratic operator exponential by its induced classical matrix exponential. Determinant-class and branch assumptions from B-070 remain in force.

<a id="b-074"></a>

#### B-074 — Fermionic quadratic propagator

**Location:** Ch.III §6.5, printed pp.147–148; PDF155–156; Lemma 6.5, Theorem 6.3, (6.37)–(6.39). **Source; visually checked.** For the homogeneous fermionic normal form, $\hat A=-\hat A'$, the same block generator $\widehat{\mathcal L}$ gives $e^{it\widehat{\mathcal L}}$. The kernel is
$$
K_{e^{it\hat H}}=\det(\hat\Phi e^{it\hat C})^{1/2}
\exp\!\left[-\tfrac12(a\ a^*)
\begin{pmatrix}\bar{\hat\Psi}\hat\Phi^{-1}&\hat\Phi'^{-1}\\-\hat\Phi^{-1}&\hat\Phi^{-1}\hat\Psi\end{pmatrix}\binom a{a^*}\right].
$$
It is initially valid where $\hat\Phi^{-1}$ exists. The source invokes limits with property A at singular points. Its square root is a continuous dynamical lift, not a fresh principal root at every time.

<a id="b-075"></a>

#### B-075 — Resolvent formula requires a convergence qualification

**Location:** Ch.III §6.6, printed p.148; PDF156; (6.40). **Source display; visually checked; limitation independently derived.** For $\hat H=\hat a^*\hat C\hat a$, the book displays
$$
K_{(\hat H-z\hat I)^{-1}}=\int_0^1\exp(a^*\xi^{\hat C}a)\,\xi^{-z-1}d\xi,\qquad\Im z\ne0.
$$
**Caution:** the stated nonreal- $z$ condition does not ensure convergence of this literal integral. At $a=a^*=0$, the integrand is $\xi^{-z-1}$, convergent only for $\Re z<0$. Generally one needs a suitable spectral/convergence half-plane, and then analytic continuation for the remaining resolvent set. The source’s equation is recorded, not promoted to an unrestricted convergent integral identity.

<a id="b-076"></a>

#### B-076 — Necessary condition for a singular quadratic Hamiltonian

**Location:** Ch.III §7.1, printed pp.148–156; PDF156–164; Theorem 7.1, (7.1)–(7.11). **Source.** Assume the corresponding classical linear system has a unique solution. If the formal quadratic expression defines a self-adjoint operator on a dense Fock domain, the induced canonical transformation must be proper. The fermionic proof uses bounded fields; the bosonic proof uses bounded exponentiated quadratures before recovering the field transformation. Thus solving the classical equations is necessary but does not prove the existence of the quantum Hamiltonian on the chosen Fock space.

<a id="b-077"></a>

#### B-077 — Sufficient condition and a visible missing identity term

**Location:** Ch.III §7.2, printed pp.156–157; PDF164–165; Theorem 7.2. **Source with contextual correction explicitly identified; p.156 visually checked.** Besides proper classical evolution, the theorem requires a determinant-class condition and meaningful locally integrable scalar terms $g\hat\Phi'^{-1}\bar{\hat A}\hat\Phi^{-1}g$, $f^*\hat\Phi^{-1}g$. The displayed condition (1) says $\hat\Phi e^{it\hat C}$ “is a trace-class operator.” Read literally this fails at $t=0$ in infinite dimension. **Contextual correction:** Lemma 6.4 and determinant (6.23) require $\hat\Phi e^{it\hat C}-\hat I$ trace class. The intended conclusion constructs a unitary one-parameter group and its self-adjoint generator. The missing $-\hat I$ is not silently erased from the audit.

<a id="b-078"></a>

#### B-078 — Sufficient reduction to quadratic normal form

**Location:** Ch.III §7.3, printed pp.157–164; PDF165–172; Theorem 7.3, Corollary 7.1, (7.12)–(7.30). **Source; pp.157–158 visually checked.** Let $\widehat{\mathfrak C} \hat X=\hat C\hat X+\hat X\hat C'$, $D_1$ be the domain of the distribution $f$, and $D_2=\{\hat X:\hat A\hat X\text{ is trace class}\}$. With bounded $\hat A$, assume $(\hat C-z\hat I)^{-1}f\in D_1$, $(\widehat{\mathfrak C}-z\hat I)^{-1}\hat A\in D_2$, $\Im z\ne0$. A proper affine/Bogoliubov transformation reduces the expression to a quadratic normal form; its natural domain is then dense and self-adjoint. When $\hat C$ multiplies by $c(x)$, the source’s sufficient weighted conditions are
$$
\int\frac{|f(x)|^2}{|c(x)-z|}dx<\infty,\quad
\iint\frac{|A(x,y)|^2}{|c(x)+c(y)-z|}dxdy<\infty.
$$

<a id="b-079"></a>

#### B-079 — Weaker conditions and infinite additive constants

**Location:** Ch.III §7.3, printed p.165; PDF173; (7.31). **Source; visually checked.** Under only $(\hat C-z\hat I)^{-1}f\in L$ and $(\widehat{\mathfrak C}-z\hat I)^{-1}\hat A$ Hilbert–Schmidt, the proper transformation can still make the pair and linear coefficients regular, but a scalar term can diverge. The author defines a regularized operator by replacing that scalar by a finite real number. This does not imply that the original formal expression was already self-adjoint, or that its absolute energy origin is uniquely determined.

<a id="b-080"></a>

#### B-080 — Continuum translation-invariant pairing can fail implementability

**Location:** Ch.III §7.4(1), printed pp.165–166; PDF173–174; (7.32)–(7.33). **Source; visually checked.** For pair kernel $A(k,k')=f(k)\delta(k+k')$, the induced anomalous block has the same delta-kernel shape. A nontrivial such block is not Hilbert–Schmidt on the continuum $L^2$ space; the transformation is improper, so the expression does not define a self-adjoint Hamiltonian in the given Fock representation. The source says the natural operator domain appears trivial, while a densely defined bilinear form can remain. A finite-mode truncation does not remove the obstruction in the continuum limit.

<a id="b-081"></a>

#### B-081 — Forced oscillator field and the sign in completing the square

**Location:** Ch.III §7.4(2), printed p.166; PDF174; (7.34)–(7.38). **Source display and independently corrected consequence; visually checked.** The book considers $\hat H=\int\omega\hat a^*\hat a+g(f\hat a^*+f^*\hat a)$, $\omega>\omega_0>0$, and sets $\hat a=\hat b-gf/\omega\,\hat I$. Its (7.38) prints a **positive** scalar $+g^2\int|f|^2/\omega$. **Direct completion gives instead**
$$
\hat H=\int\omega\hat b^*\hat b\,dk- g^2\int\frac{|f|^2}{\omega}\,dk\,\hat I.
$$
This sign correction follows already for one scalar oscillator; it is not a convention change. The finiteness/domain discussion survives, but the printed energy shift must not be copied unchanged.

<a id="b-082"></a>

#### B-082 — Weighted coupling regimes and domain dependence

**Location:** Ch.III §7.4(2a–b), printed pp.166–167; PDF174–175; (7.35)–(7.40). **Source; visually checked.** If $\int|f|^2<\infty$, Theorem 6.1 applies directly. If only $\int|f|^2/\omega<\infty$, Theorem 7.3 applies and the displaced domain is described by (7.39)–(7.40). The source points out that when $f\notin L$, domains for distinct couplings $g\ne g'$ can have trivial intersection. A shared formal expression therefore need not have a shared dense domain on which perturbation in $g$ is legitimate.

<a id="b-083"></a>

#### B-083 — Implementable displacement with divergent energy shift

**Location:** Ch.III §7.4(2c), printed pp.167–168; PDF175–176; (7.41)–(7.42). **Source; visually checked.** If $\int|f|^2/\omega^2<\infty$ but $\int|f|^2/\omega=\infty$, the displacement remains proper, while the additive energy term diverges. The regularized product
$$
\hat H_{\mathrm{reg}}=\int\omega(k)\left(\hat a^*+g\frac{f^*}{\omega}\hat I\right)
\left(\hat a+g\frac f\omega\hat I\right)dk
$$
has an operator meaning on the displaced domain. The source insists the parentheses cannot simply be removed into individually undefined terms. If the displacement norm also diverges, the transformation is improper in the original Hilbert-space representation.

<a id="b-084"></a>

#### B-084 — Bosonic canonical diagonalization theorem

**Location:** Ch.III §8.1, printed pp.168–173; PDF176–181; Theorem 8.1, (8.1)–(8.16). **Source; p.168 visually checked.** For a real, homogeneous, self-adjoint normal form
$$
\hat H=\tfrac12(\hat a^*\hat T\hat a^*+\hat a\bar{\hat T}\hat a+2\hat a^*\hat S\hat a),\quad \bar{\hat T}=\hat T,\ \bar{\hat S}=\hat S,\ \|\hat S\|<\infty,
$$
and $\hat S\pm \hat T>\mu \hat I>0$, a proper canonical transformation yields $\hat H=\hat a^*\hat\Omega\hat a+k\hat I$. Normal-form hypotheses include Hilbert–Schmidt $\hat T$. This is the book’s restricted diagonalization theorem, not a general Williamson theorem covering unstable, zero-frequency, unbounded, or indefinite quadratic Hamiltonians.

<a id="b-085"></a>

#### B-085 — Siegel-domain action and Riccati fixed point

**Location:** Ch.III §8.1, printed pp.170–173; PDF178–181; (8.8)–(8.16). **Source/re-expression.** A real symplectic block operator acts on symmetric $\hat Z=\hat X+i\hat Y$, $\hat Y>0$, by $\hat Z\mapsto(\hat A\hat Z+\hat B)(\hat C\hat Z+\hat D)^{-1}$. Properness corresponds to the restricted set $\hat Z-i\hat I$ Hilbert–Schmidt. A fixed point of the Hamiltonian flow solves a matrix/operator Riccati equation. Under the theorem’s real positive hypotheses,
$$
\hat Z=i\hat\gamma^{-1/2}(\hat\gamma^{1/2}\hat\beta\hat\gamma^{1/2})^{1/2}\hat\gamma^{-1/2},\qquad
\hat\beta=\hat S+\hat T,\quad\hat\gamma=\hat S-\hat T.
$$
Moving this fixed point to $i\hat I$ produces the diagonalizing transformation.

<a id="b-086"></a>

#### B-086 — Fermionic canonical diagonalization theorem

**Location:** Ch.III §8.2, printed pp.173–177; PDF181–185; Theorem 8.2, (8.17)–(8.33). **Source; p.173 visually checked.** For the real homogeneous self-adjoint fermionic normal form, with $[\hat S,\hat T]$ Hilbert–Schmidt,
$$
\hat S>\mu \hat I>0,\qquad \hat S^2-[\hat S,\hat T]-\hat T^2>\mu \hat I>0,
$$
a proper transformation reduces $\hat H$ to $\hat a^*\hat\Omega\hat a+k\hat I$. The proof again uses fractional-linear geometry and a fixed point; its domain is an infinite-dimensional analogue of a classical symmetric space. The conditions are explicit and stronger than arbitrary finite-dimensional CAR diagonalizability.

<a id="b-087"></a>

#### B-087 — Stability of square roots in Schatten ideals

**Location:** Ch.III §8.3, printed pp.177–180; PDF185–188; Lemma 8.1, (8.34)–(8.38). **Source; p.177 visually checked.** If $\hat C=\hat C^*$, $\hat K=\hat K^*$, $\hat C>\mu \hat I>0$, $\hat C^2+\hat K>\mu_1^2\hat I>0$, then
$$
\hat X=(\hat C^2+\hat K)^{1/2}-\hat C
$$
is trace class when $\hat K$ is trace class, and Hilbert–Schmidt when $\hat K$ is Hilbert–Schmidt. The proof uses $\hat X^2+\hat C\hat X+\hat X\hat C=\hat K$, positivity, finite-rank approximations, and a contraction argument. This is a useful operator-function result independently of the Fock-space applications; the positive spectral gaps are essential hypotheses.

### IV. What the Thirring example adds to the operator calculus

<a id="b-088"></a>

#### B-088 — Inequivalent CAR representations in the continuum

**Location:** Ch.IV §9.1, printed pp.181–184; PDF189–192; (9.1)–(9.4), footnote 4. **Source.** The field-to-particle transformation (9.3) exchanges creation and annihilation on positive/negative momentum half-lines. It preserves CAR but is improper. The auxiliary field Fock vacuum is absent from the physical-particle Fock space and conversely. Thus a self-adjoint auxiliary Hamiltonian does not automatically become a self-adjoint physical-space operator after this algebraic substitution. This is a concrete application of the infinite-mode warning in B-058.

<a id="b-089"></a>

#### B-089 — A formal CAR rearrangement can prove the false identity $k=0$

**Location:** Ch.IV §9.4, printed p.190; PDF198; footnote 13. **Source; visually checked.** Berezin evaluates a commutator of two infinite fermionic sums by two formal rearrangements. Combining the answers and CAR gives a scalar sum supposedly equal to zero; evaluating its telescoping boundary contribution gives a positive integer $k$. The lesson is an actual exhibited failure of unregulated infinite rearrangement, not merely a warning about notation. Finite algebraic identities are insufficient to justify moving summation boundaries or subtracting divergent constants.

<a id="b-090"></a>

#### B-090 — Operator approximation before taking singular limits

**Location:** Ch.IV §9.4, printed pp.191–197; PDF199–205; discussion preceding (9.18), (9.20)–(9.29). **Source.** The book replaces an improper transformation by a family of proper transformations, performs operator calculations there, and then takes limits. It explicitly notes that different approximation prescriptions can produce different final answers. The current commutators in (9.22), especially the same-current relation, depend on that prescription; the footnote on p.195 supplies an alternative limiting commutator. The dossier therefore does not present those currents as a universal regulator-independent algebra.

<a id="b-091"></a>

#### B-091 — Well-defined evolved fields without a Fock-space Hamiltonian

**Location:** Ch.IV §9.3–4, printed pp.188–198; PDF196–206; (9.14)–(9.18), (9.30)–(9.31). **Source.** A form-factor condition $\int |F(\sigma)|^2/(|\sigma|+1)\,d\sigma<\infty$ ensures that exponentiated quadratic terms in the evolved fields are actual operators. Nevertheless, the formally normal-ordered Hamiltonian $\hat H'$ is not self-adjoint in the physical Fock representation, and adding a scalar cannot repair it. The source uses a momentum/vacuum argument. Existence of the conjugated field expression therefore does not imply existence of its alleged unitary conjugating operator.

<a id="b-092"></a>

#### B-092 — Scattering and renormalized fields are separate questions

**Location:** Ch.IV §9.5–6, printed pp.199–202; PDF207–210; (9.32)–(9.39). **Source.** Asymptotic in/out fields can be unitarily related even when the finite-time fields are not. The source obtains a scattering unitary $\hat S=\theta e^{-\pi igF(0)\hat Q_1'(0)\hat Q_2'(0)}$, then a Gaussian integral representation for its normal symbol. Renormalized limiting fields may instead survive only as densely defined bilinear forms. These distinctions keep “exponential,” “field solution,” “scattering operator,” and “operator-valued distribution” from being conflated in singular models.

### V. Complete Wick appendix: contractions, perturbation series, and determinants

<a id="b-093"></a>

#### B-093 — Two-body interactions as quartic normal forms

**Location:** Appendix opening, printed pp.203–204; PDF211–212; Appendix (1)–(4). **Source/re-expression.** The normal-ordered quartic interaction in position space represents a system of particles with pair potential. Restricted to the $n$-particle sector, the source obtains $\hat H_n=-\sum_i\Delta_i+2\sum_{i<j}u(x_i-x_j)\hat I$, with its own convention for the pair-potential coefficient. Symmetric or antisymmetric wavefunctions distinguish statistics. The appendix states at its outset that the following exposition is formal and that analytical foundations are omitted.

<a id="b-094"></a>

#### B-094 — Dyson expansion and time ordering

**Location:** Appendix §1, printed pp.205–206; PDF213–214; Appendix (5)–(9). **Source/re-expression.** For $\hat H=\hat H_0+\hat V$, the source uses $\hat U(t,\tau)=e^{-it\hat H_0}e^{i(t-\tau)\hat H}e^{i\tau \hat H_0}$. Its interaction-picture integral equation yields ordered integrals or equivalently
$$
\hat U(t,\tau)=\sum_{n\ge0}\frac{i^n}{n!}\int_\tau^t\cdots\int_\tau^t
T[\hat V(t_1)\cdots \hat V(t_n)]\,dt_1\cdots dt_n.
$$
Signs reflect the book’s $e^{+it\hat H}$ convention. This is a formal perturbation expansion here; no convergence theorem for arbitrary interacting $\hat V$ is supplied.

<a id="b-095"></a>

#### B-095 — Free-field time dependence and contraction kernel

**Location:** Appendix §§1–2, printed pp.206–207; PDF214–215; Appendix (10)–(14). **Source.** For $\hat H_0=\int\omega(p)\hat a^*(p)\hat a(p)dp$,
$$
\hat a(p,t)=e^{it\omega(p)}\hat a(p),\quad
\hat a^*(p,t)=e^{-it\omega(p)}\hat a^*(p),
$$
$$
D(p_1,t_1\mid p_2,t_2)=\delta(p_1-p_2)e^{i(t_1-t_2)\omega(p_1)}.
$$
The same scalar kernel $D$, multiplied by $\hat I$, is the unequal-time CCR commutator or CAR anticommutator. The author introduces abstract symbols with these relations to avoid the nonuniqueness of coefficient expansions caused by already imposing free time dependence too early.

<a id="b-096"></a>

#### B-096 — Graded time ordering and equal-time convention

**Location:** Appendix §2, printed pp.208–209; PDF216–217; Appendix (15)–(16F). **Source; p.209 visually checked.** Bosonic $T$-ordering sums time permutations without signs; fermionic $T$-ordering includes the permutation parity. Inside $T$, fermionic fields therefore anticommute. If some times coincide, the corresponding symbols are placed in **normal order by definition**. This convention makes the later contraction vanish at equal times. It is not the alternative convention $\theta(0)=1/2$, and changing it can alter contact terms and determinants.

<a id="b-097"></a>

#### B-097 — Retarded contraction and the two-field formula

**Location:** Appendix §3A, printed pp.209–210; PDF217–218; Appendix (17), (18), (18₁). **Source; p.209 visually checked.**
$$
\Delta(p_1,t_1\mid p_2,t_2)=
\begin{cases}D(p_1,t_1\mid p_2,t_2),&t_1>t_2,\\0,&t_1\le t_2.\end{cases}
$$
Then $T[\hat a^*(1)\hat a(2)]=\hat a^*(1)\hat a(2)+\Delta(2\mid1)\hat I$ in the bosonic derivation, while like-type contractions vanish. The fermionic general formula retains the same contraction kernel but obtains signs from the graded algebra and derivative sides.

<a id="b-098"></a>

#### B-098 — Wick’s theorem as a sum over contractions

**Location:** Appendix §3A–B, printed pp.210–213; PDF218–221; Appendix (19)–(21F). **Source/re-expression.** A time-ordered product equals the normal product plus every possible set of disjoint creator–annihilator contractions. Each contraction contributes $\Delta(q,\tau\mid p,t)$, the contracted fields are removed, and the remaining fields are normal ordered. For bosons all terms have positive combinatorial weight; for fermions each pairing acquires the parity of the reordering required to perform it. The $1/k!$ in the source avoids overcounting the unordered set of $k$ contractions.

<a id="b-099"></a>

#### B-099 — Functional differential form of Wick’s theorem

**Location:** Appendix §3A–B, printed pp.212–213; PDF220–221; Appendix (21B), (21F). **Source; p.212 visually checked.** Define the contraction operator
$$
\widehat{\mathcal D}=\int\frac{\delta}{\delta a(q,\tau)}\Delta(q,\tau\mid p,t)
\frac{\delta}{\delta a^*(p,t)}\,dp\,dq\,dt\,d\tau.
$$
For a normal symbol $F$, its time-ordered symbol is $e^{\widehat{\mathcal D}}F$. Bosonic derivatives are ordinary commuting derivatives. Fermionic variables anticommute and the first derivative is the right derivative, the second the left derivative, exactly as (21F) specifies. Exponentiation truncates for a polynomial of fixed degree; extension to an infinite series is formal in this appendix.

<a id="b-100"></a>

#### B-100 — Normal symbol of interaction-picture evolution

**Location:** Appendix §4, printed pp.213–214; PDF221–222; Appendix (22)–(24F). **Source; p.214 visually checked.**
$$
N_{\hat U,1}(t,\tau\mid a^*,a)=e^{\widehat{\mathcal D}}
\exp\!\left[i\int_\tau^t V(s\mid a^*,a)ds\right].
$$
Here $V(s\mid a^*,a)$ is the scalar/Grassmann normal symbol of the interaction operator $\hat V(s)$. Only **after** the derivatives/contractions are calculated does the source substitute $a(p,t)=e^{it\omega(p)}a(p)$, $a^*(p,t)=e^{-it\omega(p)}a^*(p)$ to obtain the actual normal symbol. The bosonic and fermionic formulas share this form with the derivative conventions of B-099. Sending endpoints to infinity is a scattering construction only if the relevant limits exist.

<a id="b-101"></a>

#### B-101 — Contractions acting on an exponential source

**Location:** Appendix §5, printed p.215; PDF223; Appendix (25)–(26). **Source; visually checked.** For fermionic auxiliary sources $f,f^*$ anticommuting with all $a,a^*$,
$$
e^{\widehat{\mathcal D}}e^{f^*a+a^*f}=e^{f^*\Delta f+a^*f+f^*a}.
$$
The exponential source is an eigenfunction of the contraction operator. This identity is the elementary Gaussian generating formula used to integrate more general quadratic and quartic symbols. The fermionic source order shown here must remain fixed.

<a id="b-102"></a>

#### B-102 — Determinant identity for a quadratic symbol

**Location:** Appendix §5, printed pp.215–216; PDF223–224; Appendix (27). **Source; visually checked.** In the fermionic calculation,
$$
e^{\widehat{\mathcal D}}e^{-a^*\hat Ka}
=\det(\hat I+\hat K\hat\Delta)\exp[-a^*(\hat K^{-1}+\hat\Delta)^{-1}a].
$$
Here $\hat K$ and $\hat\Delta$ are integral maps represented by the quadratic and contraction kernels. The source derives this by introducing Gaussian Grassmann sources and completing the square. It initially assumes $\hat K^{-1}$ exists. **Derived rewriting:** $(\hat K^{-1}+\hat\Delta)^{-1}=\hat K(\hat I+\hat\Delta \hat K)^{-1}$, which makes a nonsingular- $\hat K$ continuation plausible wherever the latter inverse exists. Such continuation is distinct from the literal inverse formula, and the Appendix supplies no full infinite-dimensional determinant/convergence proof.

<a id="b-103"></a>

#### B-103 — Quartic interactions through an auxiliary commuting Gaussian

**Location:** Appendix §5, printed pp.216–217; PDF224–225; identity preceding Appendix (28), footnote 5. **Source; visually checked.** Factor a quartic kernel by $\int V_1(p_1,q_1\mid s)V_2(s\mid p_2,q_2)ds=iV(p_1,p_2\mid q_1,q_2)$. The quartic exponential becomes an integral over commuting fields $z,z^*$ with weight $e^{-\int z^*z}$, coupled quadratically to $a,a^*$ through $V_1z+V_2z^*$. This is the auxiliary-field Gaussian decoupling mechanism. The factorization is not unique; it must preserve the source’s signs, the $i$, and the order of the Grassmann quadratic factors.

<a id="b-104"></a>

#### B-104 — Retarded loops vanish, so the auxiliary determinant is one

**Location:** Appendix §5, printed p.217; PDF225; Appendix (28) and following unnumbered calculation. **Source; visually checked.** For the time-local auxiliary kernel $R(p_1,t_1\mid p_2,t_2)=R_1(p_1,p_2,t_1)\delta(t_1-t_2)$, with associated integral maps $\hat R,\hat\Delta$, each $\operatorname{Tr}(\hat R\hat\Delta)^n$ vanishes. A closed loop would require a strict cyclic time inequality; $n=1$ vanishes because $\Delta(t,t)=0$. Thus the formal trace-log expansion gives $\det(\hat I+\hat R\hat\Delta)=1$. This depends on B-096’s equal-time prescription and on treating the trace-log expansion formally; it is not a universal determinant theorem for arbitrary kernels called “retarded.”

<a id="b-105"></a>

#### B-105 — Auxiliary-field propagator representation

**Location:** Appendix §§5–6, printed pp.217–219; PDF225–227; Appendix (29), (29₁), (33)–(34). **Source; pp.217, 219 visually checked.** The interaction-picture normal symbol is an auxiliary Gaussian integral $\int e^{-a^*\hat T(z^*,z)a-z^*z}\Pi dz^*dz$, with $\hat T$ built from $\hat R(\hat I+\hat\Delta \hat R)^{-1}$ and the free time phases in (29₁). Write $\hat\omega$ for the one-particle multiplication operator with dispersion $\omega(p)$. Multiplying by free evolution gives
$$
N_{e^{it\hat H}}=\int\exp\!\left[a^*\bigl(e^{it\hat\omega}(\hat I-\hat T)-\hat I\bigr)a-zz^*\right]\Pi dz\,dz^*.
$$
The source states the final representation for both statistics, under the Appendix’s formal analytical status.



Source page: Berezin, printed p.217, one-based PDF p.225.


Berezin, printed p.219, one-based PDF p.227.


<a id="b-106"></a>

#### B-106 — Grand partition functions: determinant versus reciprocal determinant

**Location:** Appendix §7, printed pp.219–220; PDF227–228; Appendix (35B)–(37). **Source; visually checked.** The source defines $\Xi(\beta,\mu)=\operatorname{Tr}e^{-\beta(\hat H+\mu \hat N)}$; its chemical-potential sign is **plus**. Writing $\hat M=e^{-\beta(\hat\omega+\mu\hat I)}(\hat I-\hat T)$,
$$
\Xi_{B}=\int\det(\hat I-\hat M)^{-1}e^{-zz^*}\Pi dz\,dz^*,\qquad
\Xi_{F}=\int\det(\hat I+\hat M)e^{-zz^*}\Pi dz\,dz^*.
$$
These are interacting auxiliary-field formulas. They require convergence and determinant-class justification beyond the formal manipulations given here; bosonic poles cannot be ignored merely because a determinant expression exists algebraically.

<a id="b-107"></a>

#### B-107 — Cauchy Green function and free thermal traces

**Location:** Appendix §7, printed pp.220–222; PDF228–230; Appendix (38)–(43F); free normal symbol (33), printed p.219/PDF227. **Source plus explicitly derived specialization.** The operator $\hat M$ in B-106 is the Cauchy propagator $\hat G(\beta)$ for
$$
\partial_tu(p,t)+(\omega(p)+\mu)u(p,t)+\sum_qR_1(p,q,t)u(q,t)=0.
$$
The source rewrites $\Xi_{B}=\int\det(\hat I-\hat G)^{-1}e^{-zz^*}$, $\Xi_{F}=\int\det(\hat I+\hat G)e^{-zz^*}$. **Derived free specialization:** for $\hat M=e^{-\beta(\hat C+\mu\hat I)}$, $\operatorname{Tr}_{\mathcal F_{B}}\widehat\Gamma(\hat M)=\det(\hat I-\hat M)^{-1}$ and $\operatorname{Tr}_{\mathcal F_{F}}\widehat\Gamma(\hat M)=\det(\hat I+\hat M)$. In finite fermionic dimension these are finite identities. For the standard positive bosonic thermal trace, $\|\hat M\|<1$ is required even in finite dimension; $\hat M=2\hat I$ in one mode gives a divergent occupation sum, not the formal value $\det(\hat I-\hat M)^{-1}=-1$. In infinite dimension require trace-class $\hat M$ as well; trace-class $\hat M$ suffices for the fermionic trace. The related normal-order identity is $e^{\hat a^*\hat h\hat a}=:\!\exp[\hat a^*(e^{\hat h}-\hat I)\hat a]\!:$ on domains where the exponential is defined; it follows from source (33) and analytic continuation where justified.

<a id="b-108"></a>

#### B-108 — Coverage limits: BCH, Weyl/anti-normal order, and what is absent

**Location:** Entire read book; especially Ch.I §§1–3, printed pp.9–86/PDF17–94; Ch.III §6.4, printed pp.144–146/PDF152–154; Appendix §§1–7, printed pp.203–222/PDF211–230. **Coverage finding, not a numbered source theorem.** The book’s central calculus is **normal order plus matrix generating kernels**. It supplies exact quadratic exponential implementers, conjugation equations, Gaussian products, traces, and Wick contractions. It does **not** develop a systematic Weyl-symbol/Moyal calculus, anti-normal symbol calculus, arbitrary-operator BCH convergence theory, or a catalog of named optical displacement/squeezing identities. B-050–B-051 and the Hadamard interpretation B-072 are clearly marked derivations from its results. The bibliography and index were read completely; references there were not themselves fetched or treated as additional read sources.

### Audit findings and use restrictions

1. **Page coverage:** all 236 physical pages read; all 228 Arabic-numbered pages present. Page references in this dossier use the verified uniform offset $+8$, with Roman/unprinted front matter kept separate.
2. **Equation transcription:** machine extraction is not reliable for bars/primes/fermionic derivative arrows. The 49-page visual list above identifies actual image checks. Long expressions only described in B-048, B-054–B-055, B-063–B-064, and B-105 remain source lookups rather than claims of a complete independently audited transcription of every block coefficient.
3. **Visible source issues:** (6.40) lacks sufficient convergence restrictions (B-075); Theorem 7.2 condition (1) appears to omit $-\hat I$ (B-077); (7.38) has the opposite scalar sign from direct square completion (B-081). These are separated from the source’s own assertions and from the dossier’s corrections.
4. **Fermionic singular limits:** property A, two-sided $\hat\Phi^{-1}$, even-dimensional spectral approximants, parity sectors, and Pfaffian/continuous-square-root branches are not interchangeable assumptions. The finite-dimensional component argument near Lemma 5.5 should not be used to bypass the explicit even-dimensional hypotheses in Theorem 5.3.
5. **Analytical status:** Chapters I–III give stated operator/domain hypotheses, but some proofs and lengthy algebra are omitted by the author. The Wick Appendix expressly works formally. This dossier documents those results and qualifications; it does not certify every proof or every infinite-dimensional exchange of limits.
6. **No unreadable source pages or unresolved page mappings.** No claims were based on a software package. The topical omissions in B-108 describe what this fully read source does not supply.

**Entry count:** 108 stable entries, B-001–B-108.





---

<a id="book-p"></a>

## P — Perelomov — Generalized Coherent States and Their Applications


Source key **P**: A. Perelomov, *Generalized Coherent States and Their Applications*, Springer, 1986. This partition covers one-based PDF pages **1–166**, corresponding to preliminaries and printed body pp. **1–155**. The second partition continues at PDF 167. For this PDF's body, **PDF page = printed page + 11**.

The entries below use $\dagger$ for the book's operator superscript $+$, a bar for complex conjugation, $T$ for matrix transpose, and $\dagger$ for matrix Hermitian conjugation. Unless stated otherwise $[\hat a_i,\hat a_j^\dagger]=\delta_{ij}\hat I$. Prose is paraphrased. A **re-expression** keeps the source result but translates notation or normalization. A **derived correction** is not a claim that the corrected expression is printed in the book.

This transcription uses normalized hats for operators; source-page images retain the book's original notation. The same letter can have different roles in different entries.

| Notation | Role |
|---|---|
| $\hat q,\hat p,\hat a,\hat b,\hat c,\hat D,\hat S,\hat\rho,\hat A,\hat H,\hat T(g)$ | Operators on the stated Hilbert, Fock or function space; bold field vectors such as $\hat{\mathbf a}$ collect operator components. |
| $\hat J_\pm,\hat J_0,\hat K_\pm,\hat K_0,\hat C,\hat X,\hat E,\hat F,\hat G$ | Represented generators and Casimir operators. A hat is retained in a matrix realization of one of these operators. |
| $\alpha,\beta,\gamma,z,\xi,\zeta,Q,P,u,v$; $Q_{\hat A},P_{\hat A},W_{\hat A},\chi_W,\chi_N,\chi_A$ | Scalar coordinates, amplitudes and displacement parameters; scalar symbols and characteristic functions. The labels $W,N,A$ on $\chi$ indicate ordering. |
| $U,V,M,L,Z,W,K$; $A_{mn},\widetilde A_{mn}$; $A,B$ in P-039/P-045 | Finite matrices of scalar coefficients, or scalar expansion coefficients; the last pair are covariant-symbol matrices. These are distinct from hatted quantum operators. |
| $\hat I$; $I_N,I_{2N}$ | Identity operator on the space of the displayed operator relation; numerical identity matrices of the indicated dimensions. Zero in an operator relation is the zero operator. |

Group and homogeneous-space names such as $G,H,SU(2),SU(1,1)$ and the Iwasawa factors $K,A,N$ remain plain. State vectors and scalar kernels, normalizations, Gamma functions, traces, determinants, exponentials, and explicit derivatives retain standard notation.

<a id="p-001"></a>

### P-001 — CCR, canonical variables, and the Heisenberg–Weyl Lie algebra

**Locator:** P §1.1.1, printed p.8 / PDF 19, Eqs. (1.1.1)–(1.1.7); multimode version §3.1, pp.48–49 / PDF 59–60, Eqs. (3.1.1)–(3.1.10).

$$
[\hat q,\hat p]=i\hbar \hat I,\qquad \hat a=\frac{\hat q+i\hat p}{\sqrt{2\hbar}},\qquad
\hat a^\dagger=\frac{\hat q-i\hat p}{\sqrt{2\hbar}},\qquad [\hat a,\hat a^\dagger]=\hat I.
$$

In the dimensionless oscillator scaling used here, linear combinations of $\hat a,\hat a^\dagger,\hat I$ close under commutation. With $N$ modes,

$$
[\hat a_i,\hat a_j^\dagger]=\delta_{ij}\hat I,\quad [\hat a_i,\hat a_j]=[\hat a_i^\dagger,\hat a_j^\dagger]=0,
\quad \hat D(\boldsymbol\alpha)=\exp\!\sum_i(\alpha_i \hat a_i^\dagger-\bar\alpha_i \hat a_i).
$$

This closure explains why linear exponentials can be disentangled with a scalar central factor. Physical position and momentum may require oscillator-length rescaling before using the displayed definition of $\hat a$.

<a id="p-002"></a>

### P-002 — Central-commutator BCH and a differential proof

**Locator:** P §1.1.1, printed p.9 / PDF 20, Eqs. (1.1.10)–(1.1.13).

If $[\hat A,[\hat A,\hat B]]=[\hat B,[\hat A,\hat B]]=0$, then

$$
e^{\hat A}e^{\hat B}=e^{[\hat A,\hat B]/2}e^{\hat A+\hat B},\qquad
e^{\hat A+\hat B}=e^{-[\hat A,\hat B]/2}e^{\hat A}e^{\hat B}.
$$

The source proves this by defining $\hat F(t)=e^{t\hat A}e^{t\hat B}e^{-t(\hat A+\hat B)}$ and obtaining $\hat F'(t)=t[\hat A,\hat B]\hat F(t)$, $\hat F(0)=\hat I$. This is a prototype for parameter differentiation and coefficient matching in more general disentangling problems. For unbounded operators, algebraic manipulation must be interpreted on an appropriate common domain; a formal commutator identity alone is not a general analytic existence theorem.

<a id="p-003"></a>

### P-003 — Weyl composition, group commutator, and geometric phase

**Locator:** P §1.1.1, printed pp.9–10 / PDF 20–21, Eqs. (1.1.14)–(1.1.20).

$$
\hat D(\alpha)\hat D(\beta)
=e^{(\alpha\bar\beta-\bar\alpha\beta)/2}\hat D(\alpha+\beta),
\qquad \hat D(\alpha)^\dagger=\hat D(-\alpha).
$$

Consequently,

$$
\hat D(\alpha)\hat D(\beta)=e^{\alpha\bar\beta-\bar\alpha\beta}\hat D(\beta)\hat D(\alpha),
\quad \hat D(\alpha)\hat D(\beta)\hat D(-\alpha)\hat D(-\beta)
=e^{\alpha\bar\beta-\bar\alpha\beta}\hat I.
$$

For an explicitly ordered product $\hat D(\alpha_n)\cdots \hat D(\alpha_1)$ the phase is
$\operatorname{Im}\sum_{j>k}\alpha_j\bar\alpha_k$ and the resulting displacement is $\hat D(\sum_j\alpha_j)$. The phase is proportional to an oriented phase-space polygon area. Reversing the product reverses the relevant orientation; this matters when checking signs.

<a id="p-004"></a>

### P-004 — Stone–von Neumann equivalence and its domain qualification

**Locator:** P §1.1.2, printed p.11 / PDF 22, unnumbered theorem and Eqs. (1.1.23)–(1.1.24′); §8.1, p.103 / PDF 114, Eq. (8.1.12).

For a fixed nonzero central character, regular irreducible unitary representations of the finite-dimensional Weyl group are unitarily equivalent. Thus finite-mode canonical transformations can be implemented unitarily. The source explicitly warns that a corresponding assertion made only for unbounded $\hat a,\hat a^\dagger$ requires additional domain assumptions. This theorem does not remove the infinite-mode implementability conditions developed in Berezin.

<a id="p-005"></a>

### P-005 — Fock basis and exponentials applied to the vacuum

**Locator:** P §1.1.3, printed p.12 / PDF 23, Eqs. (1.1.25)–(1.1.27); §3.1, pp.48–50 / PDF 59–61, Eqs. (3.1.4)–(3.1.11).

$$
\hat a|0\rangle=0,\quad |n\rangle=\frac{(\hat a^\dagger)^n}{\sqrt{n!}}|0\rangle,
\quad \hat a|n\rangle=\sqrt n|n-1\rangle,\quad
\hat a^\dagger|n\rangle=\sqrt{n+1}|n+1\rangle.
$$

The multimode basis is $|\mathbf n\rangle=\prod_j(\hat a_j^\dagger)^{n_j}|0\rangle/\sqrt{\prod_jn_j!}$. An annihilation-only exponential acts as the identity on the vacuum. This is the final simplification after normally ordering a displacement or squeeze operator.

<a id="p-006"></a>

### P-006 — Position-space factorization of displacement

**Locator:** P §1.1.3, printed p.13 / PDF 24, Eq. (1.1.34); §1.3, p.23 / PDF 34, Eqs. (1.3.24)–(1.3.25).

With scalar displacement parameters $Q,P$ and $\alpha=(Q+iP)/\sqrt{2\hbar}$,

$$
\hat D(\alpha)=e^{-iPQ/(2\hbar)}e^{iP\hat q/\hbar}e^{-iQ\hat p/\hbar},
\quad [\hat D(\alpha)\psi](q)
=e^{iP(q-Q/2)/\hbar}\psi(q-Q).
$$

Here $q$ in the wavefunction arguments and scalar phase is the position coordinate, whereas $\hat q,\hat p$ in the exponentials are operators. These are the same Weyl operator in coordinate form. The half-phase cannot be discarded when multiplying several translations.

<a id="p-007"></a>

### P-007 — General displaced fiducials, overlaps, and resolution of identity

**Locator:** P §1.2, printed pp.14–15 / PDF 25–26, Eqs. (1.2.1)–(1.2.16).

For a normalized fiducial $|\psi_0\rangle$, define $|\alpha;\psi_0\rangle=\hat D(\alpha)|\psi_0\rangle$. Then

$$
\langle\beta;\psi_0|\alpha;\psi_0\rangle
=e^{(\alpha\bar\beta-\bar\alpha\beta)/2}
\langle\psi_0|\hat D(\alpha-\beta)|\psi_0\rangle.
$$

Displacement translates the family into itself up to the Weyl phase. For the standard normalization, $\int d^2\alpha\,|\alpha;\psi_0\rangle\langle\alpha;\psi_0|/\pi=\hat I$, in the weak sense. This gives integral kernels for operator products and for extracting matrix elements of exponential operators.

<a id="p-008"></a>

### P-008 — Displacement conjugation and invariance of centered variances

**Locator:** P §1.2, printed p.16 / PDF 27, Eqs. (1.2.19)–(1.2.22).

$$
\hat D(\alpha)^\dagger \hat a\hat D(\alpha)=\hat a+\alpha\hat I,
\quad \hat D^\dagger \hat q\hat D=\hat q+Q\hat I,\quad \hat D^\dagger \hat p\hat D=\hat p+P\hat I.
$$

Centered second moments do not change under displacement. A displaced family can therefore be centered before its squeezing or minimum-uncertainty properties are calculated.

<a id="p-009"></a>

### P-009 — Minimum uncertainty and a new annihilation operator

**Locator:** P §1.2, printed pp.16–17 / PDF 27–28, Eqs. (1.2.23)–(1.2.28).

Saturation of the appropriate uncertainty inequality gives a first-order annihilation equation. A normalized linear combination

$$
\hat b=u\hat a+v\hat a^\dagger,\qquad |u|^2-|v|^2=1
$$

again satisfies $[\hat b,\hat b^\dagger]=\hat I$. Its vacuum is a squeezed vacuum; displaced eigenstates satisfy $(u\hat a+v\hat a^\dagger)|\psi\rangle=\gamma|\psi\rangle$. This is an alternative route to deriving a squeezed state without first disentangling its unitary generator.

<a id="p-010"></a>

### P-010 — Normal and antinormal displacement factorizations

**Locator:** P §1.2, printed p.18 / PDF 29, Eqs. (1.2.31)–(1.2.32). **Visually checked.**

$$
\hat D(\alpha)=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}e^{-\bar\alpha \hat a}
=e^{+|\alpha|^2/2}e^{-\bar\alpha \hat a}e^{\alpha \hat a^\dagger}.
$$

Both scalar factors follow directly from P-002. Normal and antinormal ordering must not be interchanged while retaining the same scalar multiplier.

<a id="p-011"></a>

### P-011 — Differential correspondences for multiplying displacement

**Locator:** P §1.2, printed p.18 / PDF 29, Eqs. (1.2.33)–(1.2.34). **Visually checked.** Treat $\alpha,\bar\alpha$ as independent Wirtinger variables.

$$
\begin{aligned}
\hat a^\dagger \hat D&=(\partial_\alpha+\bar\alpha/2)\hat D,&
\hat D\hat a^\dagger&=(\partial_\alpha-\bar\alpha/2)\hat D,\\
\hat a\hat D&=(-\partial_{\bar\alpha}+\alpha/2)\hat D,&
\hat D\hat a&=(-\partial_{\bar\alpha}-\alpha/2)\hat D.
\end{aligned}
$$

These turn operator multiplication and evolution into differential equations for characteristic functions. Their left/right distinction encodes the noncommutativity.

<a id="p-012"></a>

### P-012 — Coherent-state exponentials, overlaps, and measures

**Locator:** P §1.2, printed pp.18–19 / PDF 29–30, Eqs. (1.2.35)–(1.2.41); §3.1, p.50 / PDF 61, Eqs. (3.1.11)–(3.1.12).

$$
|\alpha\rangle=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}|0\rangle
=e^{-|\alpha|^2/2}\sum_{n\ge0}\frac{\alpha^n}{\sqrt{n!}}|n\rangle,
\quad \hat a|\alpha\rangle=\alpha|\alpha\rangle,
$$

$$
\langle\alpha|\beta\rangle
=e^{-(|\alpha|^2+|\beta|^2)/2+\bar\alpha\beta},\quad
|\langle\alpha|\beta\rangle|^2=e^{-|\alpha-\beta|^2},\quad
\int\frac{d^{2N}\alpha}{\pi^N}|\boldsymbol\alpha\rangle\langle\boldsymbol\alpha|=\hat I.
$$

<a id="p-013"></a>

### P-013 — Bargmann realization of operator exponentials

**Locator:** P §1.3, printed pp.19–23 / PDF 30–34, Eqs. (1.3.3)–(1.3.11), (1.3.18)–(1.3.26).

$$
\hat a\longleftrightarrow\partial_z,\qquad \hat a^\dagger\longleftrightarrow z,
\quad f(z)=\sum_n c_n\frac{z^n}{\sqrt{n!}},\quad
\|f\|^2=\int\frac{d^2z}{\pi}e^{-|z|^2}|f(z)|^2.
$$

The reproducing kernel is $e^{z\bar w}$. Displacement becomes an ordinary weighted translation:

$$
[\hat D(\alpha)f](z)=e^{-|\alpha|^2/2+\alpha z}f(z-\bar\alpha).
$$

The coordinate–Bargmann kernel is

$$
K(z,q)=(\pi\hbar)^{-1/4}
\exp\!\left(-\frac{z^2}{2}+\sqrt{\frac2\hbar}zq-\frac{q^2}{2\hbar}\right).
$$

Useful for converting operator differential equations into Gaussian generating-function calculations. Finite pointwise growth bounds alone do not ensure Bargmann square integrability; the source gives $e^{z^2/2}$ as a boundary counterexample.

<a id="p-014"></a>

### P-014 — Weyl-lattice commutation and completeness theorems

**Locator:** P §§1.4–1.5, printed pp.23–29 / PDF 34–40, Proposition 1, unnumbered lattice theorem p.24, Eqs. (1.4.4)–(1.4.10), (1.5.1)–(1.5.19); multimode extension §3.2, pp.50–53 / PDF 61–64.

Displacements on a lattice commute when their symplectic areas obey $\operatorname{Im}(\omega_1\bar\omega_2)\in\pi\mathbb Z$. Common eigenvectors generally belong to a rigged extension rather than the ordinary Fock Hilbert space. A one-mode coherent lattice with elementary complex-plane area $S$ is overcomplete for $S<\pi$, incomplete for $S>\pi$, and complete at $S=\pi$; at equality one state may be deleted while preserving completeness, but any two deletions destroy it. The source also gives formal lattice linear dependencies and theta-function realizations. These results are useful when displacement exponentials are restricted to discrete subgroups; they are not finite disentangling formulas.

<a id="p-015"></a>

### P-015 — Covariant, contravariant, normal, and antinormal symbols

**Locator:** P §1.6, printed pp.30–34 / PDF 41–45, Eqs. (1.6.1)–(1.6.10), (1.6.19)–(1.6.29).

$$
Q_{\hat A}(\alpha)=\langle\alpha|\hat A|\alpha\rangle,\qquad
\hat A=\int\frac{d^2\alpha}{\pi}P_{\hat A}(\alpha)|\alpha\rangle\langle\alpha|.
$$

The coefficients $A_{mn},\widetilde A_{mn}$ below are scalars, while $\hat A$ denotes the operator. When the ordered expansions exist,

$$
\hat A=\sum_{mn}A_{mn}(\hat a^\dagger)^m\hat a^n
\ \Longrightarrow\ Q_{\hat A}(\alpha)=\sum_{mn}A_{mn}\bar\alpha^m\alpha^n,
$$

$$
\hat A=\sum_{mn}\widetilde A_{mn}\hat a^n(\hat a^\dagger)^m
\ \Longrightarrow\ P_{\hat A}(\alpha)=\sum_{mn}\widetilde A_{mn}\bar\alpha^m\alpha^n.
$$

$$
Q_{\hat A}(\alpha)=\int\frac{d^2\beta}{\pi}e^{-|\alpha-\beta|^2}P_{\hat A}(\beta).
$$

Normal-symbol multiplication is not ordinary multiplication: coherent kernels compose with an intermediate coherent-state integral. Every bounded operator has a covariant symbol; a regular contravariant symbol need not exist. Trace and norm conditions are discussed at pp.31–32; singular symbols require distributions and domain care.

<a id="p-016"></a>

### P-016 — Trace inequalities for exponential operators

**Locator:** P §1.6, printed p.32 / PDF 43, Eqs. (1.6.11)–(1.6.18).

For self-adjoint $\hat H$ with suitable covariant/contravariant symbols and convergent integrals,

$$
\int d\mu(\alpha)e^{-tQ_{\hat H}(\alpha)}
\le\operatorname{Tr}e^{-t\hat H}
\le\int d\mu(\alpha)e^{-tP_{\hat H}(\alpha)}.
$$

The section also recalls Golden–Thompson,
$\operatorname{Tr}e^{-t(\hat A+\hat B)}\le\operatorname{Tr}(e^{-t\hat A}e^{-t\hat B})$,
and a coherent-symbol refinement. These are useful bounds when a noncommuting exponential has no convenient exact factorization. Use the positivity, self-adjointness, and trace-class hypotheses; they are not unconditional formulas for arbitrary complex generators.

<a id="p-017"></a>

### P-017 — Distributional trace of displacement and operator reconstruction

**Locator:** P §1.6, printed p.34 / PDF 45, Eqs. (1.6.30)–(1.6.33).

$$
Q_{\hat D(\gamma)}(\alpha)
=e^{-|\gamma|^2/2+\gamma\bar\alpha-\bar\gamma\alpha},\qquad
P_{\hat D(\gamma)}(\alpha)
=e^{+|\gamma|^2/2+\gamma\bar\alpha-\bar\gamma\alpha},
$$

$$
\operatorname{Tr}\hat D(\gamma)=\pi\delta^{(2)}(\gamma),\quad
\operatorname{Tr}[\hat D(\alpha)\hat D(\beta)^\dagger]
=\pi\delta^{(2)}(\alpha-\beta),
$$

$$
\hat A=\int\frac{d^2\gamma}{\pi}\operatorname{Tr}[\hat D(\gamma)\hat A]\hat D(-\gamma).
$$

These are distributional / Hilbert–Schmidt reconstruction statements. A nonzero displacement is unitary, not trace class; the delta formula is not an ordinary convergent Fock-space trace.

<a id="p-018"></a>

### P-018 — Displacement matrix elements and Laguerre polynomials

**Locator:** P §1.6, printed p.35 / PDF 46, Eqs. (1.6.34)–(1.6.36′); further derivation in Appendix B, pp.300–304 / PDF 311–315, covered in the second partition.

$$
\langle m|\hat D(\gamma)|n\rangle=e^{-|\gamma|^2/2}
\begin{cases}
\sqrt{n!/m!}\,\gamma^{m-n}L_n^{(m-n)}(|\gamma|^2),&m\ge n,\\
\sqrt{m!/n!}\,(-\bar\gamma)^{n-m}L_m^{(n-m)}(|\gamma|^2),&n\ge m.
\end{cases}
$$

In particular $\langle n|\hat D(\gamma)|n\rangle=e^{-|\gamma|^2/2}L_n(|\gamma|^2)$. The coherent generating kernel is

$$
e^{(|\alpha|^2+|\beta|^2)/2}\langle\alpha|\hat D(\gamma)|\beta\rangle
=e^{-|\gamma|^2/2+\bar\alpha\beta+\bar\alpha\gamma-\beta\bar\gamma}.
$$

<a id="p-019"></a>

### P-019 — Weyl symbols, displaced parity, and a printed notation issue

**Locator:** P §1.6, printed pp.35–37 / PDF 46–48, Eqs. (1.6.37)–(1.6.50). Page 36 visually checked.

In the normalization $d\mu(\alpha)=d^2\alpha/\pi$,

$$
\hat T(\alpha)=2\hat D(\alpha)\hat\Pi \hat D(-\alpha)=2\hat D(2\alpha)\hat\Pi,
\quad\hat\Pi=(-1)^{\hat a^\dagger \hat a},\quad\hat\Pi \hat D(\alpha)\hat\Pi=\hat D(-\alpha),
$$

$$
W_{\hat A}(\alpha)=\operatorname{Tr}[\hat T(\alpha)\hat A],\qquad
\hat A=\int d\mu(\alpha)W_{\hat A}(\alpha)\hat T(\alpha).
$$

Fourier transforms of Weyl symbols replace the classical plane waves by $\hat D(\gamma)$; this is a direct route from a classical Gaussian exponential to its Weyl-ordered operator.

**Derived normalization correction, source page 37/PDF48 visually checked:** the Gaussian smoothing kernels with this measure are $2e^{-2|\alpha-\beta|^2}$ for $P\to W$ and $W\to Q$, and $e^{-|\alpha-\beta|^2}$ for $P\to Q$. An identity symbol must map to an identity symbol. The printed Eqs. (1.6.49)–(1.6.50) omit the factors of 2 in the two narrow Gaussian convolutions. The bracketed expression in printed Eq. (1.6.46) should not be read as an ordinary trace of a commutator: the orthogonality relation uses the product $\hat T(\alpha)\hat T(\beta)$, namely $\pi^{-1}\operatorname{Tr}[\hat T(\alpha)\hat T(\beta)]=\delta^{(2)}(\alpha-\beta)$.

<a id="p-020"></a>

### P-020 — Ordered characteristic functions and moments

**Locator:** P §1.7, printed pp.37–39 / PDF 48–50, Eqs. (1.7.8a–c)–(1.7.17). The following uses the conventional argument $\gamma$ instead of the source's $-\eta$.

$$
\chi_W(\gamma)=\operatorname{Tr}[\hat\rho \hat D(\gamma)],\quad
\chi_N(\gamma)=\operatorname{Tr}[\hat\rho e^{\gamma \hat a^\dagger}e^{-\bar\gamma \hat a}],\quad
\chi_A(\gamma)=\operatorname{Tr}[\hat\rho e^{-\bar\gamma \hat a}e^{\gamma \hat a^\dagger}],
$$

$$
\chi_N=e^{|\gamma|^2/2}\chi_W,\quad
\chi_A=e^{-|\gamma|^2/2}\chi_W,\quad
\hat\rho=\int\frac{d^2\gamma}{\pi}\chi_W(\gamma)\hat D(-\gamma).
$$

$$
\langle(\hat a^\dagger)^m\hat a^n\rangle
=(-1)^n\left.\partial_\gamma^m\partial_{\bar\gamma}^n\chi_N(\gamma)\right|_0.
$$

For a number-diagonal state, $\chi_W(\gamma)=e^{-|\gamma|^2/2}\sum_n p_nL_n(|\gamma|^2)$. The source also gives positivity criteria for quantum characteristic functions; ordinary classical Bochner positivity is replaced by a criterion including the Weyl cocycle.

<a id="p-021"></a>

### P-021 — Group-orbit coherent states and why disentangling works

**Locator:** P §§2.1–2.3, printed pp.40–44 / PDF 51–55, Eqs. (2.1.1), (2.2.1)–(2.2.6), (2.3.1)–(2.3.10).

For an irreducible unitary representation $\hat T$ and a fiducial vector, the ray stabilizer $H$ satisfies $\hat T(h)|0\rangle=e^{i\theta(h)}|0\rangle$. Coherent rays are parametrized by $G/H$. A group operator maps the family into itself up to a phase. When $\int_{G/H}|\langle0|x\rangle|^2d\mu(x)$ is finite, Schur's lemma gives a resolution of identity and a reproducing kernel. Thus decomposing $g$ into triangular and stabilizer factors immediately disentangles $\hat T(g)|0\rangle$: annihilating factors drop out and stabilizer factors become scalars.

<a id="p-022"></a>

### P-022 — Maximal isotropy and quadratic/oscillator algebras

**Locator:** P §2.4, printed pp.45–47 / PDF 56–58, Eqs. (2.4.1)–(2.4.5).

The complexified Lie algebra and the fiducial vector's isotropy subalgebra determine whether the coherent-state space admits an analytic realization. The vacuum is appropriate for the Heisenberg–Weyl and oscillator algebras; highest/lowest-weight vectors are appropriate for compact groups and positive discrete series of $SU(1,1)$. The oscillator algebra includes the quadratic number Hamiltonian as well as $\hat a,\hat a^\dagger,\hat I$. This supplies a structural test for choosing an ordered exponential ansatz rather than guessing one.

<a id="p-023"></a>

### P-023 — Gaussian decomposition is local, with a singular chart boundary

**Locator:** P §4.1, printed pp.54–57 / PDF 65–68, especially Eqs. (4.1.6)–(4.1.9); generalization §11.1, p.127 / PDF 138, Eq. (11.1.2).

On the big cell of $SL(2,\mathbb C)$, a matrix can be written uniquely as upper-triangular unipotent $\times$ diagonal $\times$ lower-triangular unipotent. The decomposition fails at a vanishing pivot, even though the group operator itself can remain perfectly regular. In a representation these three factors become raising, Cartan, and lowering exponentials. This is the algebraic origin of both $SU(2)$ and $SU(1,1)$ disentangling. It generalizes to root-space decompositions of semisimple groups, on an open dense cell, not necessarily as one global parameter chart.

<a id="p-024"></a>

### P-024 — $SU(2)$ commutators, ladder action, and differential realization

**Locator:** P §4.2, printed pp.58–59 / PDF 69–70, Eqs. (4.2.1)–(4.2.7).

$$
[\hat J_0,\hat J_\pm]=\pm \hat J_\pm,\quad[\hat J_+,\hat J_-]=2\hat J_0,
\quad \hat J_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}|j,m\pm1\rangle.
$$

On polynomials of degree at most $2j$,

$$
\hat J_- =\partial_z,\qquad \hat J_+=-z^2\partial_z+2jz\hat I,\qquad \hat J_0=z\partial_z-j\hat I.
$$

Useful realizations of this algebra include number-conserving two-boson mixing and two-fermion pair rotations. The latter needs its own CAR realization; a bosonic squeezing sign cannot simply be reused.

<a id="p-025"></a>

### P-025 — $SU(2)$ exponential disentangling

**Locator:** P §4.3, printed p.61 / PDF 72, Eqs. (4.3.11)–(4.3.14). **Visually checked.**

For $\xi=re^{i\phi}$, define $\zeta=e^{i\phi}\tan r$. On a chart with $\cos r>0$,

$$
e^{\xi \hat J_+-\bar\xi \hat J_-}
=e^{\zeta \hat J_+}e^{\log(1+|\zeta|^2)\hat J_0}e^{-\bar\zeta \hat J_-}
=e^{-\bar\zeta \hat J_-}e^{-\log(1+|\zeta|^2)\hat J_0}e^{\zeta \hat J_+}.
$$

The source parametrizes the phase through polar angles; the displayed form is its re-expression in terms of $\xi$. Beyond the chosen chart a logarithm branch and possible central phase must be tracked. Parameters are representation-independent: the book recommends checking the identity in the spin-$1/2$ matrices.

<a id="p-026"></a>

### P-026 — $SU(2)$ coherent state, overlap, annihilator, and identity resolution

**Locator:** P §4.3, printed pp.61–64 / PDF 72–75, Eqs. (4.3.10), (4.3.15), (4.3.18)–(4.3.20), (4.3.23)–(4.3.30).

$$
|\zeta;j\rangle=(1+|\zeta|^2)^{-j}e^{\zeta \hat J_+}|j,-j\rangle
=(1+|\zeta|^2)^{-j}\sum_{n=0}^{2j}\sqrt{\binom{2j}{n}}\zeta^n|j,-j+n\rangle,
$$

$$
\langle\zeta;j|\eta;j\rangle
=\frac{(1+\bar\zeta\eta)^{2j}}
{(1+|\zeta|^2)^j(1+|\eta|^2)^j},\qquad
(\hat J_-+2\zeta \hat J_0-\zeta^2\hat J_+)|\zeta;j\rangle=0,
$$

$$
\frac{2j+1}{\pi}\int\frac{d^2\zeta}{(1+|\zeta|^2)^2}
|\zeta;j\rangle\langle\zeta;j|=\hat I.
$$

The annihilator is obtained by conjugating $\hat J_-$, providing a useful check on the exponential factorization.

<a id="p-027"></a>

### P-027 — Rotation exponential kernels and oscillator contraction

**Locator:** P §4.3, printed pp.64–66 / PDF 75–77, Eqs. (4.3.36)–(4.3.45), (4.3.49)–(4.3.51).

The substitutions $\hat J_\pm\sim\sqrt{2j}(\hat a^\dagger,\hat a)$ and $\zeta=\alpha/\sqrt{2j}$ contract the spin coherent family to ordinary coherent states as $j\to\infty$. Coherent matrix elements of a spin group operator are generated by a polynomial kernel of degree $2j$ in the two coherent coordinates; traces follow by coherent-state integration. Thermal spin exponentials are obtained by analytic continuation of the rotation parameter, subject to the finite-dimensional representation. This is useful for recognizing trigonometric/hyperbolic continuations rather than treating every factorization as unrelated.

<a id="p-028"></a>

### P-028 — $SU(1,1)$ algebra and positive discrete series

**Locator:** P §§5.1.1–5.1.2, printed pp.67–71 / PDF 78–82, Eqs. (5.1.1)–(5.1.22).

$$
[\hat K_0,\hat K_\pm]=\pm \hat K_\pm,\qquad[\hat K_-,\hat K_+]=2\hat K_0,\qquad
\hat C=\hat K_0^2-\tfrac12(\hat K_+\hat K_-+\hat K_-\hat K_+)=k(k-1)\hat I.
$$

$$
\hat K_0|k,n\rangle=(k+n)|k,n\rangle,\quad
\hat K_+|k,n\rangle=\sqrt{(n+1)(2k+n)}|k,n+1\rangle,
$$

$$
\hat K_-|k,n\rangle=\sqrt{n(2k+n-1)}|k,n-1\rangle.
$$

A holomorphic realization is $\hat K_- =\partial_z$, $\hat K_+=z^2\partial_z+2kz\hat I$, $\hat K_0=z\partial_z+k\hat I$. General positive $k$ may require a covering group; do not identify every oscillator representation with a single-valued representation of $SU(1,1)$ itself. The $SU(1,1)$ matrix group is noncompact and its relevant unitary representations are infinite-dimensional.

<a id="p-029"></a>

### P-029 — $SU(1,1)$ disentangling, with a verified source sign correction

**Locator:** P §5.2.1, printed pp.73–74 / PDF 84–85, Eqs. (5.2.3), (5.2.6)–(5.2.9). **Page 74 visually checked.**

For $\xi=re^{i\phi}$, $\zeta=e^{i\phi}\tanh r$,

$$
\boxed{e^{\xi \hat K_+-\bar\xi \hat K_-}
=e^{\zeta \hat K_+}e^{\log(1-|\zeta|^2)\hat K_0}e^{-\bar\zeta \hat K_-}}
$$

$$
=e^{-\bar\zeta \hat K_-}e^{-\log(1-|\zeta|^2)\hat K_0}e^{\zeta \hat K_+}.
$$

**Derived correction:** the printed Eq. (5.2.6) has $+\eta \hat K_0$, and (5.2.8) has $-\eta \hat K_0$, while (5.2.7) defines $\eta=2\log\cosh r=-\log(1-|\zeta|^2)>0$. Those two middle signs are reversed. The corrected normal middle exponent is $-\eta \hat K_0$, and the antinormal one is $+\eta \hat K_0$. This follows both from a faithful $2\times2$ algebra calculation and from the normalized lowest-weight state in the adjacent Eq. (5.2.9). Do not copy the source signs unchanged.

<a id="p-030"></a>

### P-030 — $SU(1,1)$ coherent states, overlaps, annihilator, and measure

**Locator:** P §5.2.1, printed pp.74–75 / PDF 85–86, Eqs. (5.2.9)–(5.2.18).

$$
|\zeta;k\rangle=(1-|\zeta|^2)^ke^{\zeta \hat K_+}|k,0\rangle
=(1-|\zeta|^2)^k\sum_{n\ge0}\sqrt{\frac{\Gamma(n+2k)}{n!\Gamma(2k)}}\zeta^n|k,n\rangle,
$$

$$
\langle\zeta;k|\eta;k\rangle
=\frac{(1-|\zeta|^2)^k(1-|\eta|^2)^k}{(1-\bar\zeta\eta)^{2k}},\quad
(\hat K_--2\zeta \hat K_0+\zeta^2\hat K_+)|\zeta;k\rangle=0.
$$

For $k>1/2$,

$$
\frac{2k-1}{\pi}\int_{|\zeta|<1}\frac{d^2\zeta}{(1-|\zeta|^2)^2}
|\zeta;k\rangle\langle\zeta;k|=\hat I.
$$

The $k=1/4$ oscillator sector therefore cannot use this positive area-integral resolution naively. The state itself is well defined; it is the displayed integral realization that has a restricted range.

<a id="p-031"></a>

### P-031 — Single-mode squeezing as $SU(1,1)$

**Locator:** P §5.2.1, printed p.75 / PDF 86, Eqs. (5.2.20)–(5.2.22); factorization from corrected P-029.

$$
\hat K_+=\tfrac12(\hat a^\dagger)^2,\qquad \hat K_-=\tfrac12\hat a^2,
\quad \hat K_0=\tfrac12(\hat a^\dagger \hat a+\tfrac12\hat I),\quad \hat C=-3/16\,\hat I.
$$

The even and odd Fock sectors have $k=1/4$ and $k=3/4$. In the convention
$\hat S(\xi)=\exp[\tfrac12(\bar\xi \hat a^2-\xi(\hat a^\dagger)^2)]$, $\xi=re^{i\theta}$,

$$
\hat S(\xi)=e^{-\frac12e^{i\theta}\tanh r\,(\hat a^\dagger)^2}
e^{-\log\cosh r\,(\hat a^\dagger \hat a+1/2\hat I)}
e^{\frac12e^{-i\theta}\tanh r\,\hat a^2}.
$$

This is a **derived specialization**, with the sign convention stated explicitly. It gives the vacuum amplitude $(\cosh r)^{-1/2}$ and a pair-creation exponential.

<a id="p-032"></a>

### P-032 — Two-mode squeezing and conserved number difference

**Locator:** P §5.2.1, printed p.76 / PDF 87, Eqs. (5.2.24)–(5.2.27).

$$
\hat K_+=\hat a^\dagger \hat b^\dagger,\quad \hat K_-=\hat a\hat b,\quad
\hat K_0=\tfrac12(\hat N_a+\hat N_b+\hat I),\quad
\hat C=\tfrac14[(\hat N_a-\hat N_b)^2-\hat I].
$$

A fixed number difference $d$ gives $k=(1+|d|)/2$. The pair-vacuum state has geometric pair probabilities $(1-|\zeta|^2)|\zeta|^{2n}$. Using P-029,

$$
e^{r(e^{i\phi}\hat a^\dagger \hat b^\dagger-e^{-i\phi}\hat a\hat b)}
=e^{e^{i\phi}\tanh r\,\hat a^\dagger \hat b^\dagger}
e^{-\log\cosh r\,(\hat N_a+\hat N_b+\hat I)}
e^{-e^{-i\phi}\tanh r\,\hat a\hat b}.
$$

This is a derived specialization. Number difference is conserved, so a large Fock-space calculation separates into independent $SU(1,1)$ sectors.

<a id="p-033"></a>

### P-033 — Principal-series and Lorentz coherent systems: applicability boundary

**Locator:** P §5.2.2, printed pp.76–83 / PDF 87–94, Eqs. (5.2.28)–(5.2.81); Chs.6–7, pp.84–100 / PDF 95–111.

These pages were read in full. They construct continuous/principal-series coherent systems, overlaps as spherical functions, hyperbolic plane-wave kernels, Plancherel measures, and hypergeometric matrix elements. They are useful when a chosen Lie-algebra representation is not a lowest-weight oscillator representation. They do not supply additional finite creation/annihilation disentangling identities beyond the group decomposition already recorded. In particular, a resolution of identity for a positive discrete series must not be transferred to these non-square-integrable families without changing measure and interpretation.

<a id="p-034"></a>

### P-034 — General finite-mode bosonic Bogoliubov conditions

**Locator:** P §8.1, printed pp.101–102 / PDF 112–113, Eqs. (8.1.1)–(8.1.9). Page 102 visually checked.

For finite scalar coefficient matrices $U,V$ and the operator-vector relation $\hat{\mathbf b}=U\hat{\mathbf a}+V\hat{\mathbf a}^\dagger$,

$$
UU^\dagger-VV^\dagger=I_N,\qquad UV^T=VU^T,
\quad M=\begin{pmatrix}U&V\\\bar V&\bar U\end{pmatrix},\quad
MKM^\dagger=K,\quad K=\operatorname{diag}(I_N,-I_N),
$$

$$
M^{-1}=\begin{pmatrix}U^\dagger&-V^T\\-V^\dagger&U^T\end{pmatrix},
\quad U^\dagger U-V^T\bar V=I_N,\quad U^\dagger V=V^T\bar U.
$$

Here $M,K$ are numerical $2N\times2N$ coefficient matrices. The matrix $Z=U^{-1}V$ is symmetric, $I_N-ZZ^\dagger=(U^\dagger U)^{-1}>0$, and $U$ is invertible. These matrix constraints are a much cheaper way to verify a proposed exponential's adjoint action than testing many Fock matrix elements.

<a id="p-035"></a>

### P-035 — Quadratic bosonic closure, generators, and metaplectic phase

**Locator:** P §8.1, printed pp.102–104 / PDF 113–115, Eqs. (8.1.12)–(8.1.17). Pages 103–104 visually checked.

The bilinears

$$
\hat X_{ij}=\hat a_i\hat a_j,\quad \hat X^{ij}=\hat a_i^\dagger \hat a_j^\dagger,\quad
\hat X^i{}_j=\tfrac12(\hat a_i^\dagger \hat a_j+\hat a_j\hat a_i^\dagger)
$$

close under ordinary commutators as $\mathfrak{sp}(2N,\mathbb R)$. Their full index commutator table is printed in Eq. (8.1.14), p.103 / PDF114. A matrix generator in this realization has the form

$$
L=\begin{pmatrix}C&A\\\bar A&\bar C\end{pmatrix},\qquad
C^\dagger=-C,\quad A^T=A.
$$

Thus exponentiating a quadratic Fock operator corresponds to exponentiating a finite matrix generator. The even and odd Fock sectors are invariant. The unitary realization is naturally a representation of the metaplectic double cover, so the adjoint matrix transformation alone does not determine a global sign/phase of the implementing operator.

For an explicit common-index re-expression, put $\hat E_{ij}=\hat a_i^\dagger \hat a_j+\delta_{ij}/2\,\hat I$, $\hat F_{ij}=\hat a_i\hat a_j$, $\hat G_{ij}=\hat a_i^\dagger \hat a_j^\dagger$. Then

$$
\begin{aligned}
[\hat E_{ij},\hat E_{kl}]&=\delta_{jk}\hat E_{il}-\delta_{il}\hat E_{kj},\\
[\hat E_{ij},\hat F_{kl}]&=-\delta_{ik}\hat F_{jl}-\delta_{il}\hat F_{kj},\\
[\hat E_{ij},\hat G_{kl}]&=\delta_{jk}\hat G_{il}+\delta_{jl}\hat G_{ki},\\
[\hat F_{ij},\hat G_{kl}]&=\delta_{jk}\hat E_{li}+\delta_{jl}\hat E_{ki}+\delta_{ik}\hat E_{lj}+\delta_{il}\hat E_{kj},\\
[\hat F_{ij},\hat F_{kl}]&=[\hat G_{ij},\hat G_{kl}]=0.
\end{aligned}
$$

The source's upper/lower placement is translated into the explicitly defined $\hat E,\hat F,\hat G$ symbols here.

<a id="p-036"></a>

### P-036 — General bosonic triangular factorization and vacuum chart

**Locator:** P §8.2, printed pp.104–105 / PDF 115–116, Eqs. (8.2.1)–(8.2.10). Visually checked.

The source gives the structural factorization

$$
\hat T(g)=\mathcal N\exp[-\tfrac12\xi_{ij}\hat X^{ij}]
\exp[\alpha^k{}_l\hat X^l{}_k]\exp[-\tfrac12\eta^{ij}\hat X_{ij}].
$$

The statement supplies the ordering structure, not explicit universal coefficient functions in this display. Applied to the vacuum, the annihilation factor vanishes and the Cartan factor becomes a scalar. With $Z=U^{-1}V$,

$$
|Z\rangle=\mathcal N_B(Z)e^{-\frac12\hat{\mathbf a}^{\dagger T}Z\hat{\mathbf a}^\dagger}|0\rangle,
\quad(\hat{\mathbf a}+Z\hat{\mathbf a}^\dagger)|Z\rangle=0.
$$

<a id="p-037"></a>

### P-037 — Bosonic Gaussian normalization, overlaps, and normalizability theorem

**Locator:** P §8.2, printed pp.105–107 / PDF 116–118, Eqs. (8.2.11)–(8.2.20), unnumbered lemma p.107. Pages 105 and 107 visually checked.

$$
\mathcal N_B(Z)=\det(I_N-ZZ^\dagger)^{1/4},\quad Z^T=Z,\quad I_N-ZZ^\dagger>0,
$$

$$
\langle Z|W\rangle=
\frac{[\det(I_N-ZZ^\dagger)\det(I_N-WW^\dagger)]^{1/4}}
{\det(I_N-WZ^\dagger)^{1/2}}.
$$

The entire Bargmann Gaussian $\exp(z^TZz/2)$ is in the Bargmann Hilbert space **if and only if** $I_N-ZZ^\dagger>0$. Pair exponential series have support only on even total occupation. Equation (8.2.15) supplies a Gaussian integral for their Fock expansion coefficients, while (8.2.14) is their generating overlap identity. Determinant square-root branches must be continued from the vacuum chart.

<a id="p-038"></a>

### P-038 — Fractional-linear Gaussian parameter action

**Locator:** P §§8.1–8.2, printed pp.103,107 / PDF 114,118, Eqs. (8.1.11), (8.2.21). Visually checked.

In the source's analytic representation,

$$
Z\mapsto (U+Z\bar V)^{-1}(V+Z\bar U),
$$

$$
[\hat T^{(+)}(g)f](Z)=\det(U+Z\bar V)^{-1/2}
f\!\left((U+Z\bar V)^{-1}(V+Z\bar U)\right).
$$

Keep this function-representation convention distinct from the transformation of a covariance matrix or from an active state parameter under a differently oriented conjugation. The determinant multiplier carries information absent from the fractional-linear map alone.

<a id="p-039"></a>

### P-039 — Pure bosonic Gaussian bilinear symbols and constraints

**Locator:** P §8.3, printed pp.107–108 / PDF 118–119, Eqs. (8.3.1)–(8.3.7). Page 107 visually checked.

With $[\hat a_i,\hat a_j^\dagger]=\hbar\delta_{ij}\hat I$, the source's bilinear covariant-symbol matrices are

$$
A=-\hbar Z(I_N-Z^\dagger Z)^{-1},\qquad
B=\frac\hbar2 I_N+\hbar Z(I_N-Z^\dagger Z)^{-1}Z^\dagger,
$$

$$
B^2-A\bar A=\frac{\hbar^2}{4}I_N,\quad BA=A\bar B,
\quad B^\dagger=B,\quad A^T=A.
$$

Here $A,B$ are scalar expectation-value matrices, distinct from abstract operators. The index placement follows the source's $\hat X_{ij},\hat X^{ij},\hat X^i{}_j$ definitions. These purity constraints test proposed Gaussian expectation values and expose missing vacuum half-quanta.

<a id="p-040"></a>

### P-040 — Bosonic Gaussian geometry and representation parameter restriction

**Locator:** P §8.3, printed pp.108–110 / PDF 119–121, Eqs. (8.3.8)–(8.3.23).

The allowed $Z$ form the Siegel disk, with determinant kernels and an invariant metric. The weighted holomorphic integral in Eq. (8.3.19) converges only for the displayed parameter range $k>2N$; other admissible representation values require analytic continuation or a different realization. Its invariant-measure formula must not automatically be used as a positive integral resolution for the physical oscillator representation. This is a useful warning when turning a formal Gaussian kernel identity into an integral.

<a id="p-041"></a>

### P-041 — Fermionic CAR and finite-mode Bogoliubov constraints

**Locator:** P §9.1, printed pp.111–112 / PDF 122–123, Eqs. (9.1.1)–(9.1.8). Page 112 visually checked.

$$
\{\hat c_i,\hat c_j^\dagger\}=\delta_{ij}\hat I,\quad\{\hat c_i,\hat c_j\}=\{\hat c_i^\dagger,\hat c_j^\dagger\}=0,
$$

$$
\hat{\mathbf b}=U\hat{\mathbf c}+V\hat{\mathbf c}^\dagger,\quad
UU^\dagger+VV^\dagger=I_N,\quad UV^T+VU^T=0,
$$

$$
M=\begin{pmatrix}U&V\\\bar V&\bar U\end{pmatrix},\qquad M^\dagger M=I_{2N}.
$$

On the chart where $U$ is invertible, $Z=U^{-1}V$ is skew-symmetric. In contrast to bosons, $U$ can be singular; a fully occupied state can lie outside this vacuum chart. The connected, parity-preserving quadratic transformations give $SO(2N)$ in a Majorana basis and its spin implementation. The CAR constraints alone also allow the other orthogonal component; the source's identification with $SO(2N)$ should be read with the connected/even restriction.

<a id="p-042"></a>

### P-042 — Fermionic quadratic closure, parity, and the spin cover

**Locator:** P §§9.1–9.2, printed pp.112–113 / PDF 123–124, Eqs. (9.1.9)–(9.1.12). Visually checked.

The bilinears

$$
\hat X_{ij}=\hat c_i\hat c_j,\quad \hat X^{ij}=\hat c_i^\dagger \hat c_j^\dagger,\quad
\hat X^i{}_j=\tfrac12(\hat c_i\hat c_j^\dagger-\hat c_j^\dagger \hat c_i)
$$

close under **ordinary commutators** as $\mathfrak{so}(2N)$. Equation (9.1.11), p.113/PDF124, prints the full index table, whose relative signs differ from the bosonic table. A Hermitian infinitesimal matrix has skew pairing block $A^T=-A$ and Hermitian number block $C^\dagger=C$. Quadratic even exponentials preserve parity; Fock space splits into even and odd spinor sectors. The true unitary group action is on $\operatorname{Spin}(2N)$, with a sign ambiguity if only the orthogonal adjoint matrix is retained. The source uses $\hat T^\dagger \hat c\hat T$ at (9.1.9), whereas the bosonic (8.1.12) uses $\hat T\hat a\hat T^\dagger$; conjugation orientation must be fixed before comparing their matrices.

In the common-index convention $\hat E_{ij}=\hat c_i^\dagger \hat c_j-\delta_{ij}/2\,\hat I$, $\hat F_{ij}=\hat c_i\hat c_j$, $\hat G_{ij}=\hat c_i^\dagger \hat c_j^\dagger$, the first three commutators in the $\hat E,\hat F,\hat G$ table in P-035 retain their form, while

$$
[\hat F_{ij},\hat G_{kl}]=-\delta_{jk}\hat E_{li}+\delta_{ik}\hat E_{lj}+\delta_{jl}\hat E_{ki}-\delta_{il}\hat E_{kj},
\qquad [\hat F_{ij},\hat F_{kl}]=[\hat G_{ij},\hat G_{kl}]=0.
$$

Here $\hat F_{ij}=-\hat F_{ji}$, $\hat G_{ij}=-\hat G_{ji}$, and $\hat G_{ij}=-\hat F_{ij}^\dagger$ for distinct indices. The creation-pair order therefore matters even in defining an adjoint. This re-expression was checked in a two-mode Jordan–Wigner matrix realization for all 16 index choices.

<a id="p-043"></a>

### P-043 — Fermionic pair-exponential state, normalization, and overlap

**Locator:** P §9.2, printed pp.113–114 / PDF 124–125, Eqs. (9.2.1)–(9.2.6). Visually checked.

$$
|Z\rangle=\det(I_N+ZZ^\dagger)^{-1/4}
e^{-\frac12\hat{\mathbf c}^{\dagger T}Z\hat{\mathbf c}^\dagger}|0\rangle,
\quad Z^T=-Z,\quad(\hat{\mathbf c}+Z\hat{\mathbf c}^\dagger)|Z\rangle=0,
$$

$$
\langle Z|W\rangle=
\frac{\det(I_N+WZ^\dagger)^{1/2}}
{[\det(I_N+ZZ^\dagger)\det(I_N+WW^\dagger)]^{1/4}}.
$$

Only even occupations occur. The creation exponential is a finite polynomial in finite-mode Fock space. Its coefficient generating identity is Eq. (9.2.6). The square-root overlap branch is fixed continuously in a chosen chart; a determinant alone does not choose the global spinor sign.

<a id="p-044"></a>

### P-044 — Fermionic analytic group action

**Locator:** P §9.2, printed p.114 / PDF 125, Eq. (9.2.7). Visually checked.

$$
[\hat T^{(+)}(g)f](Z)=\det(U+Z\bar V)^{1/2}
f\!\left((U+Z\bar V)^{-1}(V+Z\bar U)\right).
$$

Compare P-038: the determinant exponent changes from $-1/2$ for bosons to $+1/2$ for fermions, with the symmetry class of $Z$ changing from symmetric to skew-symmetric. The same-looking fractional-linear formula does not make the two operator algebras interchangeable.

<a id="p-045"></a>

### P-045 — Fermionic bilinear symbols and purity constraints

**Locator:** P §9.3, printed p.115 / PDF 126, Eqs. (9.3.3)–(9.3.5). Visually checked.

With $\{\hat c_i,\hat c_j^\dagger\}=\hbar\delta_{ij}\hat I$,

$$
A=-\hbar Z(I_N+Z^\dagger Z)^{-1},\quad
B=\hbar Z(I_N+Z^\dagger Z)^{-1}Z^\dagger-\frac\hbar2I_N,
$$

$$
B^2-A\bar A=\frac{\hbar^2}{4}I_N,\quad BA=A\bar B,
\quad B^\dagger=B,\quad A^T=-A.
$$

Here $A,B$ are scalar covariant-symbol matrices. The same formal constraint has a different physical interpretation because $A$ is skew-symmetric and the number-sector shift has changed. The generalized coherent-state manifold is compact, although this $Z$ chart is not global. Geometry and invariant measures are given in Eqs. (9.3.6)–(9.3.14), pp.115–116 / PDF126–127.

<a id="p-046"></a>

### P-046 — Nilpotent groups, exponential coordinates, and square integrability

**Locator:** P Ch.10, printed pp.119–125 / PDF 130–136, definitions pp.119–120, Engel theorem p.120, §10.4 theorems pp.122–124, §10.5 theorem p.125.

The descending commutator sequence defines nilpotence. For connected simply connected nilpotent groups the exponential map gives global coordinates; sufficiently high nested Lie commutators vanish, so BCH truncates algebraically. The source develops Kirillov's orbit construction and states that square integrability is equivalent to a nondegenerate orbit form on the quotient by the center. A Pfaffian of that form tests nondegeneracy. These statements explain why the Heisenberg case is exceptionally simple and why a general nonlinear polynomial in $\hat a,\hat a^\dagger$ need not have the same finite disentangling structure. They do not establish convergence of arbitrary unbounded-operator BCH series.

<a id="p-047"></a>

### P-047 — Root-space coherent exponentials and general Gaussian decomposition

**Locator:** P Ch.11, printed pp.126–133 / PDF 137–144, Eqs. (11.1.1)–(11.1.5), (11.2.1), (11.3.6)–(11.3.12).

In a Cartan–root basis $[\hat H_i,\hat E_\alpha]=\alpha_i\hat E_\alpha$, $[\hat E_\alpha,\hat E_\beta]\propto \hat E_{\alpha+\beta}$ when that root exists. Gaussian decomposition on the big cell orders positive-root, Cartan, and negative-root factors. Applied to a lowest-weight vector,

$$
|\boldsymbol\zeta\rangle=\mathcal N(\boldsymbol\zeta)
\exp\!\left(\sum_{\alpha>0}\zeta_\alpha \hat E_\alpha\right)|-\lambda\rangle,
$$

while a unitary displacement-like generator has the form

$$
\hat D(\boldsymbol\xi)=\exp\!\left(\sum_{\alpha>0}\xi_\alpha \hat E_\alpha-\bar\xi_\alpha \hat E_{-\alpha}\right).
$$

The parameter conversion is supplied by finite-dimensional group decomposition, not by simply setting $\boldsymbol\xi=\boldsymbol\zeta$. This is the route for extending the $SU(2)$/$SU(1,1)$ examples to larger closed quadratic algebras.

<a id="p-048"></a>

### P-048 — Holomorphic discrete series, determinant kernels, and admissible measures

**Locator:** P Ch.12, printed pp.134–144 / PDF 145–155, Theorems 1–3, Propositions 1–5; Eqs. (12.2.1)–(12.2.17), (12.3.6)–(12.3.31).

The source classifies bounded symmetric matrix domains and gives determinant reproducing kernels. Of direct relevance, symmetric matrix coordinates describe bosonic pair states, while skew-symmetric coordinates occur in fermionic and related representations. A normalized coherent kernel has the form

$$
\psi_Z(W)=B(Z,Z)^{-1/2}B(Z,W),
$$

and overlap follows from the reproducing property. The group action is weighted fractional-linear substitution. Kernel powers, integration weights, and allowed representation parameters are linked; one cannot vary a determinant exponent while retaining an arbitrary positive measure. The discrete-series existence theorem requires a compact Cartan subgroup. The section also uses

$$
\det A=e^{\operatorname{Tr}\log A},\quad
-\log\det(I_N-X)=\sum_{m\ge1}\frac{\operatorname{Tr}X^m}{m},\quad \rho(X)<1,
$$

to expand determinant kernels. The logarithm branch and convergence condition in the re-expression must be fixed.

<a id="p-049"></a>

### P-049 — Principal-series kernels: locator record for the remaining first-part pages

**Locator:** P Ch.13 §§13.1–13.5, printed pp.145–155 / PDF 156–166, Eqs. (13.1.1)–(13.5.19); chapter continues in the second partition.

These pages were read for scope. They develop induced group representations, Iwasawa decomposition $G=KAN$, multiplier cocycles, horospherical coherent kernels, and rank-one symmetric-space spherical functions. Their immediate contribution to the present task is the distinction between a triangular/Gaussian factorization chart and an Iwasawa factorization, and between square-integrable lowest-weight coherent systems and principal-series systems. No additional canonical boson/fermion exponential formula is transcribed here from the unrelated hyperbolic-wave spectral tables.

### First-part source audit and coverage record

| PDF pages read | Printed pages / material | Result |
|---|---|---|
| 1–11 | Cover, title, copyright, preface, complete contents, blank leaves | Bibliographic and structural reading |
| 12–17 | Introduction pp.1–3, blank leaves, Part I divider | Read |
| 18–50 | Ch.1 pp.7–39 | Read in full; displacement, symbols, ordering and characteristic-function entries above |
| 51–58 | Ch.2 pp.40–47 | Read in full |
| 59–64 | Ch.3 pp.48–53 | Read in full |
| 65–77 | Ch.4 pp.54–66 | Read in full |
| 78–94 | Ch.5 pp.67–83 | Read in full; source middle-sign error at p.74 checked visually |
| 95–103 | Ch.6 pp.84–92 | Read in full; principal/Lorentz-series scope record P-033 |
| 104–111 | Ch.7 pp.93–100 | Read in full; principal-series scope record P-033 |
| 112–121 | Ch.8 pp.101–110 | Read in full; main Gaussian formulas checked against page images |
| 122–127 | Ch.9 pp.111–116 | Read in full; main CAR/Gaussian formulas checked against page images |
| 128–129 | Part II divider and blank leaf | Read |
| 130–136 | Ch.10 pp.119–125 | Read in full |
| 137–144 | Ch.11 pp.126–133 | Read in full |
| 145–155 | Ch.12 pp.134–144 | Read in full |
| 156–166 | Ch.13 pp.145–155 | Read in full; continuation assigned to the second partition |

OCR was performed over every supplied Perelomov PDF page because the file has no embedded text. OCR is a reading aid, not mathematical evidence: bars, primes, minus signs and subscripts are frequently damaged. Main canonical transformation, state-normalization, and disentangling displays were inspected as images. The collection does not claim that every mathematical display on every out-of-scope page was individually transcribed or independently proved.

**Correction ledger:** P-029 records the visually confirmed middle-exponent sign problem in (5.2.6)/(5.2.8). P-019 records visually confirmed missing Gaussian-normalization factors in (1.6.49)/(1.6.50), and the product/commutator notation issue in (1.6.46). P-041 states the orthogonal-component restriction implicit in the source's $SO(2N)$ identification. Exact determinant branch and conjugation conventions remain explicit wherever needed. The rotated rank-one-space table on printed p.153/PDF164 was read after correcting its image orientation.

Source: Askold Perelomov, *Generalized Coherent States and Their Applications*, Springer (1986). This dossier covers **one-based PDF pages 167–332**, covering printed body pp.156–320 together with divider/blank pages and the unprinted back cover. The first partition is documented in P-001–P-049. No software-package material was consulted.

### Complete-reading and visual audit

Every page in this partition was read sequentially in the complete OCR of this local PDF, followed by selective visual checks of important formulas and damaged OCR. This was not a keyword scan. Printed page = PDF page −11 was consistent throughout observed headers.

| Material | Printed pages | PDF pages | Coverage |
|---|---:|---:|---|
| Chapter 13 remainder, §§13.5–13.7 | 156–172 | 167–183 | complete |
| Chapters 14–15 | 173–184 | 184–195 | complete |
| Chapter 16 | 185–203 | 196–214 | complete |
| Blank, Part III title, blank | 204–206 logical positions | 215–217 | inspected; blanks/title |
| Chapter 17 | 207–210 | 218–221 | complete |
| Chapter 18 | 211–230 | 222–241 | complete |
| Chapter 19 | 231–252 | 242–263 | complete |
| Chapters 20–22 | 253–269 | 264–280 | complete |
| Chapter 23 | 270–281 | 281–292 | complete |
| Chapters 24–27 | 282–295 | 293–306 | complete |
| Appendices A, B, C | 296–307 | 307–318 | complete |
| Addendum | 308–309 | 319–320 | complete |
| Blank | 310 logical position | 321 | inspected |
| References, including Addendum references | 311–317 | 322–328 | complete |
| Blank | 318 logical position | 329 | inspected |
| Subject index | 319–320 | 330–331 | complete |
| Back cover | unnumbered, logical position 321 | 332 | complete |

**Visual checks:** PDF203, 205, 223, 225, 238, 239, 247, 249, 250, 251, 252, 260, 261, 263, 282, 283, 285, 286, 287, 294, 298, 311, 312, 313, 314, 315 (26 pages). Important formulas on these pages are distinguished from longer expressions that remain source lookups. All 166 pages were readable; no unread pages or unresolved page mappings remain in this partition. Some formula-heavy OCR is severely damaged, so untranscribed displays are not claimed as fully coefficient-verified.

**Notation.** This transcription uses normalized hats for operators, including their quantum matrix representatives. The book’s plus superscript on creation operators is transcribed as $\hat a^\dagger$; numerical coefficient-matrix adjoints $A^+$ become $A^\dagger$. A superscript $T$ denotes transpose. A **Source** formula retains its source phases, ordering, and normalization conventions; **Re-expression** is a harmless notation change; **Derived correction/consequence** is supplied here and is not passed off as an original displayed formula. Parentheses such as (B.27) refer to the book’s Appendix B numbering. Historical claims are reported as statements of this 1986 source, not as current literature assessments.


| Notation | Type |
|---|---|
| $\hat a,\hat b,\hat q,\hat p,\hat N,\hat H,\hat X_k,\hat K_j,\hat J_j,\hat S_\pm,\hat U,\hat D,\hat\rho,\hat I$ | Quantum/abstract operators; states and their labels are unhatted. $\hat T_h,\hat\Delta$ in P2-006 act on scalar functions, and $\hat T$ in P2-022 is the quantum implementer. |
| $A,B,C,D,M,Z,E;\ I_p,I_N,I_{2N},I_2$ | Numerical canonical/geometry coefficient matrices and their numerical identities in P2-002 and P2-020–P2-029. The dimension in P2-002 is $p$. Their entries and determinants are scalars. |
| $A(z,\bar w),A_1,A_2;\ H(z_1,z_2),S_{mn};\ R,\rho_{mn},D_{mn};\ \chi_k,\sigma_k,W_n$ | Scalar symbols, generating kernels, matrix elements, characteristic functions, ordering constants, and probabilities; $\hat A$ and $\hat\rho$ are operators. |
| $N$ in P2-017; $D(t)$ in P2-037; $H$ in P2-048 | Scalar normalization, numerical covariance matrix, and scalar magnetic-field strength, respectively. They differ from number operator $\hat N$, displacement $\hat D$, and Hamiltonian $\hat H$. The reflection/pair parameter $\rho$ in P2-018–P2-019, amplitudes, couplings, and phase-space coordinates remain scalar. |

**Source-display reminder.** Faithful source-page images retain the book’s original notation. The transcription beside them uses hats for operators; scalar symbols, kernels, matrix elements, and coefficient matrices remain unhatted. Differential operations, traces, determinants, and scalar function names retain standard notation.

### Geometry, quantization, and dynamics

<a id="p2-001"></a>

#### P2-001 — Principal-series coherent systems need not be square integrable

**Location:** §13.5, printed pp.156–161; PDF167–172; (13.5.20)–(13.5.62). **Source.** Rank-one spherical functions are generalized orthogonal eigenfunctions with a continuous spectral parameter, not normalizable square-integrable states over the whole noncompact space. The text explicitly gives an infinite squared norm for the CS system. Plancherel expansions use the Harish-Chandra $c$-function and spectral measure proportional to $|c(\lambda)|^{-2}d\lambda$. Their “completeness” is a generalized Fourier statement, distinct from a finite constant resolution of identity for a square-integrable discrete-series orbit.

<a id="p2-002"></a>

#### P2-002 — Symplectic tube domains

**Location:** §§13.6–13.6.2, printed pp.161–167; PDF172–178; (13.6.1)–(13.6.23). **Source/re-expression.** The Siegel upper half-space is
$$
\mathfrak H_p=\{Z=X+iY:Z^T=Z,\ Y>0\}\simeq Sp(p,\mathbb R)/U(p)
$$
in the book’s group-dimension convention. Real symplectic blocks act by $Z\mapsto(AZ+B)(CZ+D)^{-1}$, with $A^TC=C^TA$, $B^TD=D^TB$, $A^TD-C^TB=I_p$. The other classical tubes use positive complex-Hermitian, quaternion-Hermitian, or Lorentz cones. These provide the geometric setting for Gaussian squeezing matrices.

<a id="p2-003"></a>

#### P2-003 — Horospherical kernels as curved plane waves

**Location:** §§13.6.3–13.7, printed pp.167–172; PDF178–183; (13.6.24)–(13.7.12). **Source.** The horospherical kernels are Laplace–Beltrami eigenfunctions, constant on horocycles; their overlaps are spherical functions. The $c$-function has a root-product formula (13.7.9), and generalized completeness/orthogonality yield Fourier-type expansions. Quaternion/octonion details and the full root-product coefficients are read but not independently transcribed here. These generalized CS are not ordinary bosonic eigenstates of annihilation operators.

<a id="p2-004"></a>

#### P2-004 — Completeness of discrete coherent subsystems

**Location:** Chapters 14–15, printed pp.173–184; PDF184–195; Proposition 1, Theorems 1–4, (14.1.1)–(15.2.4). **Source.** An orbit subsystem in the analytic discrete-series realization is incomplete iff a nonzero function in the representation Hilbert space vanishes on all its labels. For $SU(1,1)$ the source proves sufficiency of $S_\Gamma<\pi/(2k-1)$ for completeness, using its hyperbolic-area convention. Automorphic forms give sufficient incompleteness criteria $k>m_0^++1/2$; the general bounded-domain chapter gives $k>m_0+1$. Critical equalities and area normalizations must not be silently replaced by a universal density rule.

<a id="p2-005"></a>

#### P2-005 — Covariant symbols, analytic continuation, and product kernel

**Location:** §16.3.1, printed p.192; PDF203; (16.3.7)–(16.3.12). **Source; visually checked.** For the disk’s reproducing-kernel vectors,
$$
A(z,\bar z)=\frac{\langle z|\hat A|z\rangle}{\langle z|z\rangle},\qquad
A(z,\bar w)=\frac{\langle z|\hat A|w\rangle}{\langle z|w\rangle}.
$$
The off-diagonal symbol is determined by analytic continuation of the diagonal. Formula (16.3.12) reconstructs the product symbol as an integral with the squared normalized overlap kernel. A second, contravariant symbol appears in the projector integral (16.3.10); the two symbols are connected by the Berezin transform and should not be equated.

<a id="p2-006"></a>

#### P2-006 — Curved symbol products and the correspondence principle

**Location:** §16.3.2–3, printed pp.193–194; PDF204–205; (16.3.13)–(16.3.22). **Source; p.194 visually checked.** With $h=(2k)^{-1}$,
$$
\hat T_hf=f+h\hat\Delta f+o(h),\quad
\hat\Delta=(1-|z|^2)^2\partial_z\partial_{\bar z},
$$
$$
A_1\star A_2=A_1A_2+h(1-|z|^2)^2
(\partial_{\bar z}A_1)(\partial_zA_2)+o(h).
$$
Thus the scaled star commutator tends to the Poisson bracket with the source’s $i$ convention. This is an asymptotic statement, not an exact replacement of operator multiplication by pointwise multiplication at finite $h$.

<a id="p2-007"></a>

#### P2-007 — Curved analogues of Weyl quantization and compact quantization

**Location:** §§16.3.5–16.4, printed pp.196–198; PDF207–209; (16.3.28)–(16.4.9). **Source.** The author defines a disk Weyl symbol using point-inversion operators, related to covariant/contravariant symbols by invariant integral transforms. On a sphere, $h^{-1}$ must be an integer and the holomorphic Hilbert space is finite-dimensional; the reproducing kernel is $(1+z\bar w)^{1/h}$. These constructions are curved-phase-space analogues of ordering calculi, not the flat oscillator Weyl symbol with unchanged kernels.

<a id="p2-008"></a>

#### P2-008 — General Kähler quantization

**Location:** §16.5, printed pp.199–203; PDF210–214; Theorems 1–5, (16.5.1)–(16.5.27). **Source/re-expression.** For a Kähler potential $F$, holomorphic states carry the weighted norm $c(h)\int|f|^2e^{-F/h}d\mu$. The basis sum defines a reproducing kernel, operator symbols reconstruct bounded operators, and homogeneity is used to prove the correspondence principle. Projective group action compensates the change $F(gz,\overline{gz})-F(z,\bar z)$ by holomorphic plus antiholomorphic terms. Theorems 4–5 describe Berezin-transform eigenvalues via gamma functions; these lengthy OCR-damaged products remain source lookups.



Source page: Perelomov, printed p.202, one-based PDF p.213.


Perelomov, printed p.203, one-based PDF p.214.


<a id="p2-009"></a>

#### P2-009 — Linear Lie-algebra Hamiltonians preserve coherent orbits

**Location:** Chapter 17, printed pp.207–210; PDF218–221; (17.1)–(17.13). **Source.** If $\hat H(t)=\sum_kh^k(t)\hat X_k$ lies in a represented finite Lie algebra, its evolution is a group-representation operator (up to the relevant phase/cover). An initial generalized CS stays on the CS orbit. Compact-algebra Hamiltonians can be conjugated into a Cartan subalgebra for spectra; noncompact cases need orbit classification. The text discusses scattering when asymptotic limits exist and Floquet operators. A global logarithm in the same real Lie algebra is not guaranteed merely because the Floquet matrix belongs to the group.

<a id="p2-010"></a>

#### P2-010 — Forced oscillator: displacement and accumulated phase

**Location:** §18.1, printed pp.211–213; PDF222–224; (18.1.1)–(18.1.17). **Source/re-expression.** For $\hat H=\hat H_0-f(t)\hat q$ in the interaction picture, evolution is a displacement times a scalar phase. Writing $\dot\alpha=\beta(t)$, the phase in $|\psi(t)\rangle=e^{-i\phi(t)}|\alpha(t)\rangle$ obeys $\dot\phi=\operatorname{Im}(\bar\beta\alpha)$. The source relates it to symplectic area. The net displacement is proportional to the force’s resonant Fourier component. If the initial state is vacuum, the final number distribution is Poisson, $W_n=e^{-|\gamma|^2}|\gamma|^{2n}/n!$.

<a id="p2-011"></a>

#### P2-011 — Forced-oscillator transitions are displacement matrix elements

**Location:** §18.1, printed p.213; PDF224; (18.1.14)–(18.1.16). **Source/re-expression.** $W_{mn}=|\langle m|\hat D(\gamma)|n\rangle|^2$, so associated Laguerre polynomials control transitions from arbitrary initial number states. Appendix B derives the complete matrix elements with complex phases (P2-051–P2-054). The probability alone discards the scalar evolution phase, so it is insufficient to reconstruct the propagator for interferometric calculations.

<a id="p2-012"></a>

#### P2-012 — Parametric oscillator realizes two $SU(1,1)$ sectors

**Location:** §18.2, printed pp.213–214; PDF224–225; (18.2.1)–(18.2.7). **Source; p.214 visually checked.** Quadratic oscillator operators realize $\mathfrak{su}(1,1)$, with
$$
\hat K_0=\tfrac14(\hat p^2/\omega_0+\omega_0\hat q^2),\qquad \hat K_2=\tfrac14(\hat p\hat q+\hat q\hat p).
$$
Even and odd number sectors carry discrete-series parameters $k=1/4$ and $k=3/4$. Time dependence of the frequency therefore produces squeezed generalized CS in each parity sector. The full oscillator space is the direct sum of these sectors, not one irreducible representation.

<a id="p2-013"></a>

#### P2-013 — Exact Riccati dynamics and phase for squeezing

**Location:** §18.2, printed p.214; PDF225; (18.2.8)–(18.2.9). **Source Riccati equation and independently corrected phase; original page visually rechecked.** In units $\hbar=1$, for
$\hat H=i(\beta \hat K_+-\bar\beta \hat K_--i\gamma \hat K_0)$ and
$|\psi\rangle=e^{-i\phi}|\zeta\rangle$,
$$
\dot\zeta=\beta-\bar\beta\zeta^2-i\gamma\zeta,\qquad
i\dot\phi=k(\bar\beta\zeta-\beta\bar\zeta+i\gamma),\quad|\zeta|<1.
$$
**Source correction:** printed (18.2.9) has the two $\beta$ terms with opposite signs. For normalized $|\zeta;k\rangle=(1-|\zeta|^2)^k e^{\zeta\hat K_+}|k,0\rangle$, the vacuum coefficient of the Schrödinger equation fixes the phase shown above. The printed sign leaves the exact residual $4ik\operatorname{Im}(\bar\beta\zeta)$, which is generically nonzero; the corrected phase gives zero. The Riccati equation is unchanged. The classical motion is on the hyperbolic disk. The source describes its quantum evolution as exactly captured by the corresponding classical equations, with the additional accumulated phase. These formulas do not require an adiabatic approximation.

<a id="p2-014"></a>

#### P2-014 — Floquet spectral types of a quadratic oscillator

**Location:** §18.2, printed pp.215–216; PDF226–227; (18.2.12)–(18.2.21). **Source, with scope qualification.** For a Lie-algebra logarithm of the period map, the Lorentz norm $\Omega_0^2-\Omega_1^2-\Omega_2^2$ classifies the algebra element: timelike gives discrete spectra $\pm\hbar\Omega(k+n)$, spacelike gives a two-sided continuous spectrum, a nonzero null generator gives a one-sided continuous spectrum, while the zero generator has spectrum $\{0\}$. Classical instability corresponds to the noncompact cases. Quasienergies are defined modulo $2\pi\hbar/T$; logarithm branches and the metaplectic lift must be specified. The source’s direct exponential treatment is not a proof that every real symplectic period matrix has the presumed real logarithm.

<a id="p2-015"></a>

#### P2-015 — Singular oscillator and operator domains

**Location:** §18.3.1–2, printed pp.217–221; PDF228–232; (18.3.1)–(18.3.34). **Source/re-expression.** For the half-line inverse-square oscillator, raising/lowering operators generate $\mathfrak{su}(1,1)$, giving $E_n=2n+E_0$ and Laguerre wavefunctions. The book chooses the regular solution behaving as $x^a$, with $a(a-1)=2g^2$, $E_0=a+1/2$, and $k=(a+1/2)/2$. As $g\to0$ this gives the half-line oscillator with the vanishing boundary condition, not the full-line oscillator spectrum. The selected boundary condition must accompany the algebraic formulas.

<a id="p2-016"></a>

#### P2-016 — Many-body inverse-square models inherit the same algebra

**Location:** §18.3.3, printed pp.222–226; PDF233–237; (18.3.35)–(18.3.64). **Source.** After removing center-of-mass motion, homogeneous degree $-2$ interactions admit collective raising/lowering operators obeying the same $\mathfrak{su}(1,1)$ algebra. Radial excitations form equidistant ladders; new angular/homogeneous solutions generate different ladders. For the Calogero interaction, polynomial solutions yield degeneracy generating functions $\prod_{j=2}^N(1-t^j)^{-1}$. This is a selected family plus a separation construction, not the claim that one coherent orbit spans every many-body sector without additional labels.

<a id="p2-017"></a>

#### P2-017 — Visible source error in the affine Gaussian normalization

**Location:** §18.4, printed p.227; PDF238; (18.4.3)–(18.4.4). **Source display plus derived correction; visually checked.** The source writes
$|\zeta,\alpha\rangle=N e^{\zeta \hat a^{\dagger2}/2+\alpha \hat a^\dagger}|0\rangle$
with $N=(1-|\zeta|^2)^{-1/2}e^{-|\alpha|^2/2}$. This is not normalized. For the displayed state, the positive normalization is
$$
N=(1-|\zeta|^2)^{1/4}
\exp\!\left[-\frac{|\alpha|^2+\operatorname{Re}(\bar\zeta\alpha^2)}
{2(1-|\zeta|^2)}\right].
$$
At $\alpha=0$, the squeezed-vacuum factor alone disproves the printed exponent. An independent Fock-coefficient recurrence at $\zeta=.3,\alpha=.4$ gives squared norm 1.2337242593 with the source factor and 1.0000000000 with the corrected factor.

<a id="p2-018"></a>

#### P2-018 — Gaussian generating kernel for a forced parametric oscillator

**Location:** §18.4, printed pp.227–230; PDF238–241; (18.4.6)–(18.4.20). **Source; p.228 visually checked.** The amplitude-generating function $H(z_1,z_2)=\sum_{mn}S_{mn}z_1^mz_2^n/\sqrt{m!n!}$ is a Gaussian. Formula (18.4.7) parametrizes it by the classical reflection parameter $\rho=|c_2/c_1|^2<1$, a force parameter $\nu=|d|^2$, and a relative phase. Its vacuum amplitude has modulus $(1-\rho)^{1/4}e^{-\nu(1-\sqrt\rho\cos2\varphi)/2}$. Taylor coefficients are two-variable Hermite polynomials. The full long kernel is a source lookup; the normalization error in (18.4.4) is not used to derive it.



Source page: Perelomov, printed p.228, one-based PDF p.239.


Perelomov, printed p.229, one-based PDF p.240.


Perelomov, printed p.230, one-based PDF p.241.


<a id="p2-019"></a>

#### P2-019 — Bosonic pair creation: two-mode $\mathfrak{su}(1,1)$

**Location:** §19.2.1, printed pp.233–236; PDF244–247; (19.2.1)–(19.2.14). **Source; p.236 visually checked.** $\hat K_+=\hat a^\dagger \hat b^\dagger$, $\hat K_-=\hat a\hat b$, $\hat K_0=(\hat a^\dagger \hat a+\hat b^\dagger \hat b+\hat I)/2$ preserve number difference. Starting from vacuum,
$$
|\zeta\rangle=\sqrt{1-|\zeta|^2}\,e^{\zeta \hat b^\dagger \hat a^\dagger}|0,0\rangle,\quad
W_n=(1-\rho)\rho^n,\quad\rho=|\zeta|^2=|v/u|^2.
$$
The canonical coefficients obey $|u|^2-|v|^2=1$. Each momentum block is finite-dimensional at the one-particle level; global continuum implementability requires an additional integrability analysis.

<a id="p2-020"></a>

#### P2-020 — Multimode bosonic canonical blocks

**Location:** §19.2.2, printed p.237; PDF248; (19.2.15)–(19.2.19). **Source/re-expression.** For $(\hat a,\hat b^\dagger)\mapsto M(\hat a,\hat b^\dagger)$, $M=\left(\begin{smallmatrix}A&B\\C&D\end{smallmatrix}\right)$ preserves the indefinite form $E=\operatorname{diag}(I_N,-I_N)$:
$$
MEM^\dagger=E,\quad AA^\dagger-BB^\dagger=I_N,\quad
DD^\dagger-CC^\dagger=I_N,\quad AC^\dagger=BD^\dagger.
$$
Inverse relations follow from $M^{-1}=EM^\dagger E$. The source uses $SU(N,N)$; preservation alone gives $U(N,N)$, with the determinant-one restriction removing an overall phase component. Keep the group convention separate from the displayed canonical identities.

<a id="p2-021"></a>

#### P2-021 — Multimode paired vacuum and its determinant normalization

**Location:** §19.2.2, printed p.238; PDF249; (19.2.24)–(19.2.32). **Source; visually checked.**
$$
|Z\rangle=\det(I_N-Z^\dagger Z)^{1/2}
\exp\!\left(\sum_{ij}Z_{ij}\hat a_i^\dagger \hat b_j^\dagger\right)|0,0\rangle,
\quad I_N-Z^\dagger Z>0,
$$
with $Z=(A^\dagger)^{-1}C^\dagger=BD^{-1}$ in the source’s transformation convention. Number-basis amplitudes sum products of $Z_{ij}$ over all nonnegative occupation arrays with prescribed row/column sums. This is the bipartite-pair Gaussian normalization, with determinant exponent $1/2$, distinct from the general symmetric single-family squeezed vacuum exponent $1/4$.

<a id="p2-022"></a>

#### P2-022 — Bosonic overlaps and fractional-linear action

**Location:** §19.2.2, printed p.239; PDF250; (19.2.33)–(19.2.42). **Source; visually checked.**
$$
\langle Z_1|Z_2\rangle=
\frac{\det(I_N-Z_1^\dagger Z_1)^{1/2}\det(I_N-Z_2^\dagger Z_2)^{1/2}}
{\det(I_N-Z_1^\dagger Z_2)}.
$$
The orbit is a matrix ball $I_N-Z^\dagger Z>0$. Group action is fractional linear, $Z'=(AZ-B)(-CZ+D)^{-1}$ in the source convention, with a determinant phase. Equation (19.2.40) gives the full coherent matrix element and (19.2.41) the number-basis generating function. The signs depend on whether the chosen operator transforms $\hat a$ by $\hat T^\dagger \hat a\hat T$ or $\hat T\hat a\hat T^\dagger$.

<a id="p2-023"></a>

#### P2-023 — Bosonic pair-number generating determinant

**Location:** §19.2.2, printed pp.240–241; PDF251–252; (19.2.43)–(19.2.54). **Source; visually checked.**
$$
\mathcal G(t)=\sum_{n\ge0}W_nt^n
=\frac{\det(I_N-Z^\dagger Z)}{\det(I_N-tZ^\dagger Z)}.
$$
The eigenvalues $0\le\rho_j<1$ of $Z^\dagger Z$ give independent geometric-mode factors $\prod_j(1-\rho_j)/(1-t\rho_j)$. **Derived correction:** (19.2.54) prints a minus sign before $\sum_k\operatorname{Tr}[(Z^\dagger Z)^k]t^k/k$; the determinant ratio requires a **plus** sign. Its derivative at $t=0$ must be nonnegative. The determinant ratio (19.2.53), rather than the erroneous trace-log rewriting, is the reliable source formula.

<a id="p2-024"></a>

#### P2-024 — Multimode Riccati equation and vacuum phase

**Location:** §19.2.3, printed pp.241–242; PDF252–253; (19.2.55)–(19.2.62). **Source; p.241 visually checked.** For
$\hat H=A_{ij}\hat a_i^\dagger \hat b_j^\dagger+\bar A_{ij}\hat a_i\hat b_j+C_{ij}\hat a_i^\dagger \hat a_j+D_{ij}\hat b_i^\dagger \hat b_j$,
$C=C^\dagger,D=D^\dagger$,
$$
i\dot Z=A+ZA^\dagger Z+CZ+ZD^T,\qquad
\dot\phi=\tfrac12\operatorname{Tr}(A^\dagger Z+Z^\dagger A).
$$
These equations determine the normalized CS trajectory and its scalar phase. When pairing vanishes asymptotically and number-conserving blocks tend to constants, singular values of $Z$ tend to limits determining the pair-number statistics.

<a id="p2-025"></a>

#### P2-025 — Single fermionic pair is an $SU(2)$ problem

**Location:** §§19.1,19.3.1–2, printed pp.231–233, 242–247; PDF242–244, 253–258; (19.1.8)–(19.1.14), (19.3.21)–(19.3.31). **Source.** For one particle–antiparticle pair, $\hat J_+=\hat a^\dagger \hat b^\dagger$, $\hat J_-=\hat b\hat a$, $\hat J_0=(\hat a^\dagger \hat a+\hat b^\dagger \hat b-\hat I)/2$ act as spin $1/2$ on vacuum and the occupied pair. $|A|^2+|B|^2=1$ and pair probability is $|B|^2$. The normalized state is $(1+|\zeta|^2)^{-1/2}(\hat I+\zeta \hat a^\dagger \hat b^\dagger)|0\rangle$. Fixed-direction or planar external fields reduce the Dirac dynamics to one or two commuting $SU(2)$ spin equations; arbitrary directions use $SU(4)$ for a spin-$1/2$ field.

<a id="p2-026"></a>

#### P2-026 — Multimode fermionic canonical blocks

**Location:** §19.3.3, printed p.248; PDF259; (19.3.32)–(19.3.36). **Source/re-expression.** The corresponding doubled matrix is unitary:
$$
MM^\dagger=I_{2N},\quad AA^\dagger+BB^\dagger=I_N,\quad
DD^\dagger+CC^\dagger=I_N,\quad AC^\dagger=-BD^\dagger.
$$
The source uses $SU(2N)$, again with the determinant-one restriction beyond preservation alone. CAR signs replace the bosonic indefinite metric by a positive metric. Number-difference sectors carry exterior-power representations; finite fermionic occupation numbers are zero or one.

<a id="p2-027"></a>

#### P2-027 — Fermionic paired vacuum

**Location:** §19.3.3, printed pp.249–250; PDF260–261; (19.3.41)–(19.3.49). **Source; visually checked.**
$$
|Z\rangle=\det(I_N+Z^\dagger Z)^{-1/2}
e^{\sum_{ij}Z_{ij}\hat a_i^\dagger \hat b_j^\dagger}|0,0\rangle,\qquad
Z=-(A^\dagger)^{-1}C^\dagger=BD^{-1}.
$$
The exponential is a finite polynomial, and unlike the bosonic case no contraction bound on $Z$ is needed in this finite chart. Amplitudes are signed permutation sums, equivalently minors/determinants after a basis-order convention. Singular $D$ or $A$ requires another Grassmannian chart; the inverse formula is not globally valid.

<a id="p2-028"></a>

#### P2-028 — Fermionic overlaps and Grassmannian action

**Location:** §19.3.3, printed p.250; PDF261; (19.3.49)–(19.3.56). **Source; visually checked.**
$$
\langle Z_1|Z_2\rangle=
\frac{\det(I_N+Z_1^\dagger Z_2)}
{\det(I_N+Z_1^\dagger Z_1)^{1/2}\det(I_N+Z_2^\dagger Z_2)^{1/2}}.
$$
The coherent orbit is $SU(2N)/S(U(N)\times U(N))$, represented in this chart by unrestricted complex $N\times N$ matrices. Group action is $Z'=(AZ+B)(CZ+D)^{-1}$, with phase $\arg\det(CZ+D)$. Both the determinant sign in the overlap and the denominator in the action differ from the bosonic matrix ball.

<a id="p2-029"></a>

#### P2-029 — Fermionic pair-number generating determinant

**Location:** §19.3.3, printed pp.251–252; PDF262–263; (19.3.57)–(19.3.67). **Source; p.252 visually checked.**
$$
\mathcal F(t)=\sum_{n=0}^NW_nt^n
=\frac{\det(I_N+tZ^\dagger Z)}{\det(I_N+Z^\dagger Z)}.
$$
**Derived interpretation:** with $\lambda_j\ge0$ eigenvalues of $Z^\dagger Z$, this factors into Bernoulli laws with pair probabilities $p_j=\lambda_j/(1+\lambda_j)$, whereas P2-023 factors into geometric laws. These are finite-dimensional determinants. A continuum product requires summability, Fredholm determinant, and representation assumptions not established by this blockwise chapter alone.

<a id="p2-030"></a>

#### P2-030 — Representation multiplicities and finite/infinite restrictions

**Location:** §§19.2.2–19.3.3, printed pp.237–252; PDF248–263; (19.2.20)–(19.2.23), (19.3.37)–(19.3.40). **Source plus qualification.** The source fixes particle–antiparticle number difference and evaluates Casimir invariants to identify coherent representation sectors. Its finite matrix determinant formulas describe individual spin/momentum blocks. Constancy of Casimirs alone is not a general proof of irreducibility without the accompanying representation construction, and a finite-block unitary is not automatically implementable over infinitely many continuum blocks. The dossier records these limits rather than silently importing an infinite-mode theorem.

<a id="p2-031"></a>

#### P2-031 — Clebsch–Gordan coefficients from polynomial coherent states

**Location:** Chapter 20, printed pp.253–255; PDF264–266; (20.1)–(20.12). **Source.** The tensor-product lowest-weight polynomial is proportional to $(z_2-z_1)^{j_1+j_2-j}$. Applying $e^{\zeta \hat J_+}$ produces the finite polynomial generating function (20.9); its coefficients in $z_1,z_2,\zeta$ are normalized $SU(2)$ Clebsch–Gordan coefficients. This is a useful exact coefficient-generating technique, but normalization factors in the long display remain source lookups here.



Source page: Perelomov, printed p.253, one-based PDF p.264.


Perelomov, printed p.254, one-based PDF p.265.


Perelomov, printed p.255, one-based PDF p.266.


<a id="p2-032"></a>

#### P2-032 — Quasiclassical partition-function bounds

**Location:** Chapter 21, printed pp.256–259; PDF267–270; Simon theorem, (21.1)–(21.15). **Source.** For the specified compact-group highest-weight scaling and multiaffine lattice Hamiltonian, coherent lower and contravariant upper symbols bound the normalized quantum partition function:
$Z_{\rm cl}(\beta)\le Z_Q^{(N)}(\beta)\le Z_{\rm cl}(\beta(1+a/N))$.
The source distinguishes fundamental weights for the upper bound and states the lower bound more generally. The classical phase space is the appropriate coadjoint orbit, not always a sphere or a flat phase plane.

<a id="p2-033"></a>

#### P2-033 — Bilinear fermion algebras and large-$N$ mechanics

**Location:** Chapter 22, printed pp.260–269; PDF271–280; (22.1.11)–(22.3.8). **Source.** Fermionic bilinears close an orthogonal Lie algebra; coherent paired states parametrize a constrained nonlinear phase space. In the singlet sector $1/N$ serves as an effective Planck constant and operator products differ from symbol products by $O(1/N)$. The classical matrices obey quadratic constraints (22.3.7), not arbitrary unconstrained component equations. The text connects this with Gross–Neveu dynamics and regularized gap equations; it does not supply a general solution of all resulting constrained field equations.

### Gaussian channels, phase-space distributions, and operator functions

<a id="p2-034"></a>

#### P2-034 — Thermal oscillator master equation

**Location:** §23.1.1, printed pp.270–271; PDF281–282; (23.1.1)–(23.1.4). **Source/re-expression.** In the interaction picture,
$$
\dot{\hat\rho}=\gamma(\nu+1)\left(\hat a\hat\rho \hat a^\dagger-\tfrac12\{\hat a^\dagger \hat a,\hat\rho\}\right)
+\gamma\nu\left(\hat a^\dagger\hat\rho \hat a-\tfrac12\{\hat a\hat a^\dagger,\hat\rho\}\right),
\quad \nu=(e^{\hbar\omega/kT}-1)^{-1}.
$$
For the damping thermal case $\gamma\ge0,\nu\ge0$. Amplitudes decay as $e^{-\gamma t/2}$ and occupation excess as $e^{-\gamma t}$. The source also discusses an amplifying sign regime; it must not be confused with the thermal-attenuation assumptions used below.

<a id="p2-035"></a>

#### P2-035 — Normal, Weyl, and anti-normal characteristic functions

**Location:** §23.1.2, printed p.271; PDF282; (23.1.5)–(23.1.6). **Source; visually checked.**
$$
\chi_N(\eta)=\operatorname{Tr}(\hat\rho e^{\eta \hat a^\dagger}e^{-\bar\eta \hat a}),\quad
\chi_O(\eta)=\operatorname{Tr}(\hat\rho e^{\eta \hat a^\dagger-\bar\eta \hat a}),\quad
\chi_A(\eta)=\operatorname{Tr}(\hat\rho e^{-\bar\eta \hat a}e^{\eta \hat a^\dagger}),
$$
$$
\chi_k=e^{-\sigma_k|\eta|^2}\chi_N,\qquad
\sigma_N=0,\quad\sigma_O=\tfrac12,\quad\sigma_A=1.
$$
This explicitly displays the ordering-dependent Gaussian factors. The source uses $O$ for the symmetric/Weyl case. These factors are exact central-commutator consequences, not a small-parameter approximation.

<a id="p2-036"></a>

#### P2-036 — Exact attenuation solution for characteristic functions

**Location:** §23.1.2, printed p.272; PDF283; (23.1.7)–(23.1.8). **Source; visually checked.** For $p=e^{-\gamma t}$, $q=1-p$,
$$
\chi_k(\eta,t)=e^{-q(\nu+\sigma_k)|\eta|^2}
\chi_k(\sqrt p\,\eta,0).
$$
The corresponding PDE is
$\partial_t\chi_k=-\gamma[\tfrac12\eta_j\partial_{\eta_j}+(\nu+\sigma_k)|\eta|^2]\chi_k$.
This one formula solves the channel for arbitrary initial density operators admitting the characteristic function; singular quasiprobabilities can be treated distributionally rather than forced into ordinary function spaces.

<a id="p2-037"></a>

#### P2-037 — Quasiprobability and Fokker–Planck forms

**Location:** §23.1.2, printed pp.272–273; PDF283–284; (23.1.9)–(23.1.15). **Source/re-expression; p.272 visually checked.** Fourier transformation gives $W_N=P$, $W_O=W_{\rm Wigner}$, $W_A=Q/\pi$. The real phase-space drift is $-\gamma\alpha/2$ and diffusion matrix $\gamma(\nu+\sigma_k)I_2/2$. Means contract by $\sqrt p$; covariance in each representation evolves as
$D(t)=pD(0)+q(\nu+\sigma_k)I_2/2$.
Different representations differ by the fixed vacuum-noise offset; their probability interpretations are not identical.

<a id="p2-038"></a>

#### P2-038 — Visible missing minus sign in the Gaussian convolution

**Location:** §23.1.2, printed p.272; PDF283; (23.1.12). **Source display plus derived correction; visually checked.** The source’s displayed exponential is positive in
$|\alpha-\sqrt p\,\alpha'|^2/[q(\nu+\sigma_k)]$.
It must instead be
$$
W_k(\alpha,t)=\int d^2\alpha'\,W_k(\alpha',0)
\frac{\exp[-|\alpha-\sqrt p\,\alpha'|^2/(q(\nu+\sigma_k))]}
{\pi q(\nu+\sigma_k)}.
$$
The correction follows by Fourier-transforming the decaying characteristic Gaussian (23.1.8); the printed positive exponential is nonintegrable. When $q(\nu+\sigma_k)=0$, use the delta-kernel limit rather than division by zero.

<a id="p2-039"></a>

#### P2-039 — Gaussian states stay Gaussian under relaxation

**Location:** §23.1.2, printed pp.273–274; PDF284–285; (23.1.16)–(23.1.23). **Source.** A thermal initial state with mean occupation $\nu_0$ remains thermal with $\mu(t)=p\nu_0+q\nu$. A displaced thermal state contracts its displacement by $\sqrt p$ and has the same covariance evolution. Phase averaging gives an $I_0$-Bessel radial distribution; number-state input gives Laguerre-polynomial quasiprobabilities. The Gaussian formulas have smooth distributional limits at vanishing noise, but literal denominators in the source must not be evaluated at zero width.

<a id="p2-040"></a>

#### P2-040 — Analytic density kernel generates all Fock matrix elements

**Location:** §23.1.3, printed p.274; PDF285; (23.1.24)–(23.1.26). **Source; visually checked.**
$$
R(\bar\alpha,\beta)=e^{(|\alpha|^2+|\beta|^2)/2}\langle\alpha|\hat\rho|\beta\rangle,
\qquad
R(z_1,\bar z_2)=\sum_{mn}\rho_{mn}\frac{z_1^m\bar z_2^n}{\sqrt{m!n!}}.
$$
It obeys a second-order linear PDE with mixed derivative, linear drift, and quadratic multiplication terms. This is the unnormalized coherent kernel, not the diagonal $Q$ function or the normal symbol of the operator. Extracting off-diagonal matrix elements requires the two independent analytic arguments.

<a id="p2-041"></a>

#### P2-041 — Gaussian Green kernel for driven relaxation

**Location:** §23.1.3, printed pp.274–275; PDF285–286; (23.1.27)–(23.1.34). **Source; visually checked.** The forced master equation adds a linear Hamiltonian commutator. Formula (23.1.30) gives $R(t)$ as a double Bargmann-Gaussian integral against a Gaussian Green kernel (23.1.32)–(23.1.33). Its scalar parameters are $p=e^{-\gamma t}$, $q=1-p$, thermal noise $\nu$, and the driven classical amplitude. The complete lengthy exponent was read and visually checked but remains a source lookup here. The force-phase convention must be checked against the chosen Hamiltonian before copying the printed amplitude formula (23.1.34).



Source page: Perelomov, printed p.274, one-based PDF p.285.


Perelomov, printed p.275, one-based PDF p.286.


<a id="p2-042"></a>

#### P2-042 — Attenuation, additive noise, and displacement factorization

**Location:** §23.1.3, printed pp.275–276; PDF286–287; discussion following (23.1.34), (23.1.35). **Source/re-expression; visually checked.** The source decomposes the density evolution into vacuum attenuation, added Gaussian fluctuations, and a unitary displacement
$\hat\rho\mapsto \hat D(v)\hat\rho \hat D(v)^\dagger$,
$\hat D(v)=e^{v\hat a^\dagger-\bar v \hat a}$.
The first two maps are nonunitary. **Derived convention check:** if the interaction Hamiltonian is $\hat H_d=\epsilon(t)\hat a^\dagger+\bar\epsilon(t)\hat a$, the amplitude satisfies
$\dot\beta=-(\gamma/2)\beta-i\epsilon(t)$.
This directly fixes any missing $i$ or sign when translating the source’s physical-force notation.

<a id="p2-043"></a>

#### P2-043 — Displaced thermal density kernel and pure-state limit

**Location:** §23.1.3, printed pp.276–277; PDF287–288; (23.1.36)–(23.1.41). **Source; p.276 visually checked.**
$$
R(z_1,\bar z_2;t)=\frac1{1+\mu}
\exp\!\left[\frac{\mu z_1\bar z_2+\beta z_1+\bar\beta\bar z_2-|\beta|^2}
{1+\mu}\right],
\quad \mu=p\nu_0+q\nu,\quad\beta=\sqrt p\,\alpha+v(t).
$$
At zero bath temperature with coherent input, $\mu=0$ and the state stays pure coherent under attenuation and drive. For nonzero bath occupation, thermal noise makes it mixed. The source’s use of “superposition” here means adding fields/noise, not a coherent ket sum or an arbitrary convex mixture.

<a id="p2-044"></a>

#### P2-044 — Populations of a displaced thermal state

**Location:** §23.1.3, printed p.277; PDF288; (23.1.40)–(23.1.41). **Source/re-expression.**
$$
W_n=\frac{\mu^n}{(1+\mu)^{n+1}}
e^{-|\beta|^2/(1+\mu)}
L_n\!\left(-\frac{|\beta|^2}{\mu(1+\mu)}\right).
$$
For $\mu>0$ the Laguerre argument is nonpositive and $W_n\ge0$. The $\mu\to0$ limit gives $e^{-|\beta|^2}|\beta|^{2n}/n!$, even though the literal Laguerre expression has apparent inverse powers of $\mu$. This is a useful operator-exponential-to-observable bridge with an explicitly stated limiting regime.

<a id="p2-045"></a>

#### P2-045 — Number-state input and spin relaxation

**Location:** §§23.1.3–23.2, printed pp.277–281; PDF288–292; (23.1.42)–(23.2.12). **Source.** Number-state relaxation has a Laguerre generating kernel; the source gives special population formulas and generating functions. The spin analogue replaces $\hat a,\hat a^\dagger$ by $\hat S_-,\hat S_+$, preserving trace, Hermiticity, and positivity. Spin coherent symbols convert it to a Fokker–Planck equation on the sphere; its large-spin drift becomes a classical transport equation. The long exact hypergeometric solution at p.281 is read but not independently transcribed or used as a modern general channel theorem.

<a id="p2-046"></a>

#### P2-046 — Free oscillator propagator and thermal trace

**Location:** Chapter 24, printed pp.282–283; PDF293–294; (24.1)–(24.7). **Source; thermal formulas visually checked.** For $\hat H=\omega \hat a^\dagger \hat a$ with the vacuum energy omitted,
$$
\langle\alpha|e^{-it\hat H}|\beta\rangle=
\exp\!\left[\bar\alpha\beta e^{-i\omega t}
-\tfrac12(|\alpha|^2+|\beta|^2)\right],
$$
$$
\langle\alpha|e^{-\hat H/T}|\beta\rangle=
\exp[\bar\alpha\beta e^{-\omega/T}-(|\alpha|^2+|\beta|^2)/2],
\quad Z=(1-e^{-\omega/T})^{-1}.
$$
The trace uses $d^2\alpha/\pi$. Restoring vacuum energy multiplies the trace by $e^{-\omega/(2T)}$. The resolvent formula (24.3) likewise needs its convergence/analytic-continuation prescription, not just a formal time integral.

<a id="p2-047"></a>

#### P2-047 — Landau-level and determinant applications

**Location:** Chapter 24, printed pp.283–285; PDF294–296; (24.8)–(24.21). **Source.** The oscillator trace combines with momentum integration and Landau degeneracy to produce the charged-particle partition function. The author then obtains magnetization, weak-field susceptibility, and a contour representation for the degenerate Fermi gas. These applications illustrate why scalar energy shifts, degeneracies, and trace measures must accompany an oscillator exponential. Numerical coefficients beyond the visually checked thermal trace remain source lookups here, not independently validated predictions.

<a id="p2-048"></a>

#### P2-048 — Quadratic trace in the Heisenberg–Euler application is analytically continued

**Location:** Chapter 25, printed pp.286–288; PDF297–299; (25.5)–(25.17). **Source with domain qualification; p.287 visually checked.** The stable oscillator factor is
$\operatorname{Tr}e^{-seH(2\hat N+\hat I)}=[2\sinh(seH)]^{-1}$ for positive real $seH$.
The source formally diagonalizes $\hat b^2+\hat b^{\dagger2}$ using
$\hat U=e^{i\pi(\hat b^{\dagger2}-\hat b^2)/8}$ and writes a sine-denominator trace. This $\hat U$ is not a bounded unitary squeeze: its generator is Hermitian and it is a complex canonical continuation. Trace cyclicity and oscillatory traces require regularization/contour prescriptions. The source subsequently subtracts ultraviolet terms and identifies imaginary parts from proper-time poles; the intermediate manipulation is not an unrestricted trace-class identity.

<a id="p2-049"></a>

#### P2-049 — Synchrotron and entropy applications

**Location:** Chapters 26–27, printed pp.289–295; PDF300–306; (26.1)–(27.15). **Source.** Synchrotron calculations use coherent expectations of quadratic exponentials with explicit weak-field/relativistic approximations. Chapter 27 compares von Neumann entropy with the coherent-state (Wehrl) entropy using concavity and the resolution of identity. The source reports the oscillator minimum 1 and states the spin minimum as a conjecture at its publication date; this is a historical statement, not a present-day literature claim. One-particle entropy bounds distinguish Bose enhancement and fermionic holes, and require the source’s occupancy/normalization conventions.

### Complete appendices and displacement formulas

<a id="p2-050"></a>

#### P2-050 — Entire-function completeness and the critical coherent lattice

**Location:** Appendix A, printed pp.296–299; PDF307–310; (A.1)–(A.16). **Source.** Completeness reduces to absence of an entire function vanishing on all labels with finite Bargmann Gaussian norm. A regular lattice of cell area $S$ is complete for $S\le\pi$, incomplete for $S>\pi$, in the source’s dimensionless coherent coordinate. At the critical area $S=\pi$, deleting one lattice state leaves a complete system; deleting two makes it incomplete. The proof uses an adjusted Weierstrass sigma function and Gaussian growth. “Complete” here does not assert that the critical lattice is an orthonormal basis or a bounded frame.

<a id="p2-051"></a>

#### P2-051 — Displacement kernel and finite coefficient sum

**Location:** Appendix B, printed p.300; PDF311; (B.1)–(B.3). **Source; visually checked.** For $\hat D(\gamma)=e^{\gamma \hat a^\dagger-\bar\gamma \hat a}$,
$$
\sum_{mn}D_{mn}(\gamma)\frac{\bar\alpha^m\beta^n}{\sqrt{m!n!}}
=e^{-|\gamma|^2/2+\bar\alpha\beta+\bar\alpha\gamma-\beta\bar\gamma},
$$
$$
D_{mn}=e^{-|\gamma|^2/2}\sqrt{m!n!}
\sum_{r=0}^{\min(m,n)}
\frac{\gamma^{m-r}(-\bar\gamma)^{n-r}}
{r!(m-r)!(n-r)!}.
$$
The finite sum is an exact formula with all phase information; it avoids ambiguity in associated-Laguerre conventions.

<a id="p2-052"></a>

#### P2-052 — Displacement matrix elements in associated Laguerre form

**Location:** Appendix B, printed pp.300–301; PDF311–312; (B.4)–(B.9). **Source; visually checked.**
$$
D_{mn}=
\begin{cases}
\sqrt{n!/m!}\,e^{-|\gamma|^2/2}\gamma^{m-n}L_n^{m-n}(|\gamma|^2),&m\ge n,\\
\sqrt{m!/n!}\,e^{-|\gamma|^2/2}(-\bar\gamma)^{n-m}L_m^{n-m}(|\gamma|^2),&m\le n.
\end{cases}
$$
Equation (B.8) relates negative-index Laguerre forms. The source derives a Laguerre bound from $|D_{mn}|\le1$. Cases $\gamma=0$ should be evaluated by continuity or the finite sum, giving $\delta_{mn}$ without dividing by powers of $\gamma$.

<a id="p2-053"></a>

#### P2-053 — Anti-normal displacement factorization and Gaussian integral

**Location:** Appendix B, printed pp.301–302; PDF312–313; (B.10)–(B.18). **Source; visually checked.**
$$
\hat D(\gamma)=e^{|\gamma|^2/2}e^{-\bar\gamma \hat a}e^{\gamma \hat a^\dagger}.
$$
This gives the exact Bargmann integral for the matrix elements (B.10), Laguerre integrals (B.11)–(B.12), and generating functions
$$
\sum_{n\ge0}t^nL_n^k(x)=(1-t)^{-k-1}e^{-tx/(1-t)},\quad |t|<1,
$$
$$
\sum_{n\ge0}t^nL_n^{m-n}(x)=(1+t)^m e^{-tx}.
$$
The anti-normal scalar has the opposite sign to the normal-ordered displacement scalar; exchanging these factors without changing the prefactor is an error.

<a id="p2-054"></a>

#### P2-054 — Visible error in the generating function of displacement probabilities

**Location:** Appendix B, printed p.304; PDF315; (B.26)–(B.27). **Source display plus derived correction; visually checked.** The source defines $F(u,v;\gamma)=\sum_{mn}u^mv^n|D_{mn}(\gamma)|^2$ but (B.27) includes an extra $e^{|\gamma|^2}$. The corrected identity is
$$
F(u,v;\gamma)=\frac1{1-uv}
\exp\!\left[-\frac{|\gamma|^2(1-u)(1-v)}{1-uv}\right],
$$
initially in a domain of convergence such as $|u|,|v|<1$. At $u=v=0$ it gives $|D_{00}|^2=e^{-|\gamma|^2}$, whereas the source display gives 1. At $\gamma=0$ it gives $(1-uv)^{-1}$, consistent with $\hat D=\hat I$. Thus the correction is fixed independently of any convention about Laguerre polynomials.

<a id="p2-055"></a>

#### P2-055 — Semiclassical geometry and derivative formula for displacement

**Location:** Appendix B, printed pp.302–303; PDF313–314; (B.17)–(B.25), Figure B.1. **Source; visually checked.** Large-number Bargmann functions concentrate near circles of radii $\sqrt n$ and $\sqrt m$. Displacement transitions avoid exponential suppression when these shifted circles intersect:
$|\sqrt n-|\gamma||<\sqrt m<\sqrt n+|\gamma|$.
Equation (B.23) gives the oscillation-averaged classical intersection integral. The exact derivative formula is
$$
D_{mn}(\alpha)=\frac{e^{|\alpha|^2/2}}{\sqrt{m!n!}}
(-\partial_{\bar\alpha})^m(\partial_\alpha)^n e^{-|\alpha|^2}.
$$
The geometric rule is asymptotic and does not replace the exact finite-$n$ matrix element near turning boundaries.

<a id="p2-056"></a>

#### P2-056 — Jacobians of classical-domain transformations

**Location:** Appendix C, printed pp.305–307; PDF316–318; (C.1)–(C.18). **Source.** The author differentiates the fractional-linear group action to compute the complex Jacobian for each classical domain. For a type-I $p\times q$ matrix domain the determinant factor is $\det(B^Tz+D)^{-(p+q)}$ in the source’s transposed-block convention. The type-IV vector domain is treated separately with its quadratic invariant. The appendix explains the determinant powers used in analytic Hilbert-space measures; those powers are geometry-dependent and cannot be transferred unchanged among symmetric, skew, rectangular, and vector domains.

<a id="p2-057"></a>

#### P2-057 — Addendum, references, and historical scope

**Location:** Addendum and end matter, printed pp.308–320; PDF319–331; back cover PDF332. **Coverage finding.** The Addendum surveys further applications and distinguishes the author’s group-orbit CS from more general continuous representations such as Barut–Girardello states. All bibliography entries, Addendum references, and index pages in this partition were read. They were not fetched as additional sources. Blank PDF215, 217, 321, 329 and Part III title PDF216 were inspected and counted in the complete physical-page coverage, without inventing printed page content.

<a id="p2-058"></a>

#### P2-058 — Source-audit cautions for reuse

**Location:** Cross-references within this partition: (18.4.4), (19.2.54), (23.1.12), (B.27), and (25.10)–(25.11). **Audit findings.** Four literal source problems are separated from derived corrections: affine-Gaussian normalization (P2-017), bosonic trace-log sign (P2-023), diffusion-kernel exponent sign (P2-038), and displacement-probability generating prefactor (P2-054). Complex canonical continuation in the proper-time trace is separately qualified in P2-048. Long special-function expressions left as source lookups are not claims that every OCR coefficient has been audited. The $SU(1,1)$ disentangling sign correction recorded in P-029 also applies when that identity is used in these later applications.

**Final count:** 58 stable entries, P2-001–P2-058. **Coverage:** all 166 assigned PDF pages fully read, 26 selected pages visually checked, no unread or unreadable pages; formulas explicitly marked as source lookup remain untranscribed rather than invented.



---

<a id="book-c"></a>

## C — Combescure–Robert — Coherent States and Applications in Mathematical Physics


Source: Didier Robert and Monique Combescure, *Coherent States and Applications in Mathematical Physics*, second edition, Springer (2021), DOI 10.1007/978-3-030-70845-0. The supplied PDF has 582 actual pages; the 583rd form-feed extraction segment is empty. No QuantumFramework material or Wolfram Language code is used.

**Status:** complete reading of all 582 supplied PDF pages, including front matter, all 17 chapters, Appendices A–E, references, and both indexes; 97 topic-focused entries. Full extraction was followed by selective rendered-page checks. Complete coverage, the nonuniform page map, inferred omitted blank leaves, and remaining formula ambiguities are recorded at the end.

**Conventions.** Unless an entry states otherwise, $[\hat a_j,\hat a_k^\dagger]=\delta_{jk}\hat I$, $\hat a=(\hat Q+i\hat P)/\sqrt{2\hbar}$, $z=(q,p)$, $\sigma(z,z')=q\cdot p'-p\cdot q'$, $\hat T(z)=\exp[i(p\cdot \hat Q-q\cdot \hat P)/\hbar]$, $\alpha=(q+ip)/\sqrt{2\hbar}$. The Hermitian scalar product is antilinear in its first slot. `Source` means a formula/theorem actually present; `Re-expression` translates its notation; `Consequence` is an explicitly identified deduction. Unbounded-operator identities require their stated common cores. Beware that the source contains typographical and normalization inconsistencies; the errata ledger records them instead of silently propagating them.

**Hat notation.** This dossier re-expresses the source in normalized operator notation. A hat marks an operator on a Hilbert, Fock, Bargmann, or operator space, including a matrix being used as a quantum operator. Faithful source-page images retain the book's original notation. Hats do not change product order, domains, scalar phases, or any recorded correction. Zero in an operator equation means the zero operator.

| Object type | Dossier notation |
|---|---|
| Quantum operators and identities | $\hat a,\hat c,\hat Q,\hat P,\hat N,\hat H,\hat U,\hat T,\hat S,\hat\rho,\hat I$; one-particle Hilbert operators also receive hats. |
| Scalar symbols, kernels, amplitudes, states | $H(z),A_w,A_c,K(x,y),a(t),s,\psi,\Psi$ have no hats. A symbol of $\hat A$ is scalar even when its label refers to that operator. |
| Numerical coefficient matrices and identities | $F_t,A_t,B_t,U,V,Z,\Gamma_t,I_n$ stay plain when locally defined as such; a quantum operator's matrix representative retains its hat. |
| Grassmann parameters and standard functions | $\theta,\gamma,\alpha$ are anticommuting scalars in the fermionic entries; $\Gamma(x),\exp,\operatorname{Tr},\partial_z$ retain standard notation. |

### Linear exponentials, coherent states, and symbols

<a id="c-001"></a>

#### C-001 — CCR and the representation domain
**Source:** §1.2.1, printed 5–6; PDF 22–23; (1.5)–(1.7). Position and momentum act on their natural multiplication/Sobolev domains and satisfy $[\hat P_j,\hat Q_k]=-i\hbar\delta_{jk}\hat I$ on the product domain. The real linear combination $p\cdot \hat Q-q\cdot \hat P$, initially on $D(\hat Q)\cap D(\hat P)$, has a self-adjoint realization generating the unitary $\hat T(z)$. These statements concern real phase-space labels; complex labels require the analytic core in C-005.

<a id="c-002"></a>

#### C-002 — Central-commutator BCH with operator hypotheses
**Source:** §1.2.1, Lemma 1, printed 6; PDF 23; (1.8)–(1.10).
$$
e^{\hat A+\hat B}=e^{-[\hat A,\hat B]/2}e^{\hat A}e^{\hat B}.
$$
Here $\hat A,\hat B$ are anti-self-adjoint; there is a dense common core invariant under both generators and both real one-parameter groups; $[\hat A,[\hat A,\hat B]]=[\hat B,[\hat A,\hat B]]=0$ on that core; and $i[\hat A,\hat B]$ has a self-adjoint extension. The proof differentiates $e^{-t^2[\hat A,\hat B]/2}e^{t\hat A}e^{t\hat B}$. **Consequence:** reversing the factors reverses the central correction; do not apply the truncation when nested commutators survive.

<a id="c-003"></a>

#### C-003 — Weyl multiplication and central extension
**Source:** §1.2.1, printed 6–9; PDF 23–26; (1.11), (1.12), (1.17), (1.18).
$$
\hat T(z)\hat T(z')=e^{-i\sigma(z,z')/(2\hbar)}\hat T(z+z'),\qquad
\hat T(z)\hat T(z')=e^{-i\sigma(z,z')/\hbar}\hat T(z')\hat T(z),\quad \hat T(z)^{-1}=\hat T(-z).
$$
The central-coordinate multiplication is $(t,z)(t',z')=(t+t'+\sigma(z,z'),z+z')$, represented by $e^{-it/(2\hbar)}\hat T(z)$. **Consequence:** the phase for a product of many displacements is $\exp[-i\sum_{r<s}\sigma(z_r,z_s)/(2\hbar)]$; this follows inductively and is not a separately numbered source formula.

<a id="c-004"></a>

#### C-004 — Hadamard conjugation and coordinate factorization
**Source:** §1.2.1, Lemma 2 and Corollary 1, printed 7; PDF 24; (1.13), (1.14); §1.2.2, printed 10; PDF 27; (1.19).
$$
\hat T(z)\hat Q\hat T(z)^{-1}=\hat Q-q\hat I,\quad \hat T(z)\hat P\hat T(z)^{-1}=\hat P-p\hat I,
\quad \hat T(z)=e^{-iq\cdot p/(2\hbar)}e^{ip\cdot \hat Q/\hbar}e^{-iq\cdot \hat P/\hbar},
$$
$$
(\hat T(z)u)(x)=e^{ip\cdot(x-q/2)/\hbar}u(x-q).
$$
The conjugation proof differentiates an exponential conjugate, so it supplies the terminating Hadamard-series mechanism for linear CCR generators.

<a id="c-005"></a>

#### C-005 — Complex labels and a usable analytic core
**Source:** §1.3.1, Lemma 4, printed 12–13; PDF 29–30; (1.34). Let $\mathcal H_0$ be the span of $\phi_\nu(x)e^{\eta\cdot x}$, $\nu\in\mathbb N^n,\eta\in\mathbb C^n$. The map $z\mapsto \hat T(z)u$ extends analytically to complex phase space for $u\in\mathcal H_0$, preserves this core, and obeys the Weyl product with the **bilinear** complex extension of $\sigma$. This is the book's justification for BCH involving separate nonunitary creation/annihilation exponentials; it is stronger and more precise than treating them as bounded operators.

<a id="c-006"></a>

#### C-006 — Normal-ordered displacement
**Source:** §1.3.1, printed 11–13; PDF 28–30; (1.27)–(1.36), Proposition 1.
$$
\hat D(\alpha)=\hat T(z)=e^{\alpha\cdot \hat a^\dagger-\bar\alpha\cdot \hat a}
=e^{-|\alpha|^2/2}e^{\alpha\cdot \hat a^\dagger}e^{-\bar\alpha\cdot \hat a},
\quad \hat D(\alpha)\hat a\hat D(\alpha)^{-1}=\hat a-\alpha\hat I.
$$
**Consequence:** $\hat D(\alpha)=e^{+|\alpha|^2/2}e^{-\bar\alpha\cdot \hat a}e^{\alpha\cdot \hat a^\dagger}$. More generally, on the analytic core, $e^{u\cdot \hat a^\dagger+v\cdot \hat a}=e^{u\cdot v/2}e^{u\cdot \hat a^\dagger}e^{v\cdot \hat a}$, with independent complex $u,v$; the dot is bilinear.

<a id="c-007"></a>

#### C-007 — Vacuum action, eigenvalue, and number expansion
**Source:** §1.3.1, printed 12–15; PDF 29–32; (1.33), (1.36), (1.39)–(1.41), Lemmas 7–8.
$$
|\alpha\rangle=e^{-|\alpha|^2/2}e^{\alpha\cdot \hat a^\dagger}|0\rangle,
\qquad \hat a_j|\alpha\rangle=\alpha_j|\alpha\rangle,
\qquad |\alpha\rangle=e^{-|\alpha|^2/2}\sum_{\nu\in\mathbb N^n}\frac{\alpha^\nu}{\sqrt{\nu!}}|\nu\rangle.
$$
The final expression is the consistent re-expression of the surrounding derivation; the typeset (1.41) has garbled indices/factorial formatting, so cite (1.36) and Lemma 7 when exact normalization matters.

<a id="c-008"></a>

#### C-008 — Displacement expectation and coherent overlap
**Source:** §1.3.1, Proposition 2, Lemmas 5–6, printed 13–14; PDF 30–31; (1.37), (1.38).
$$
\langle0|\hat T(z)|0\rangle=e^{-|z|^2/(4\hbar)},\qquad
\langle\phi_z,\phi_{z'}\rangle=e^{i\sigma(z,z')/(2\hbar)}e^{-|z-z'|^2/(4\hbar)},
$$
$$
\hat T(z)\phi_{z'}=e^{-i\sigma(z,z')/(2\hbar)}\phi_{z+z'}.
$$
**Re-expression:** $\langle\alpha|\beta\rangle=\exp[-(|\alpha|^2+|\beta|^2)/2+\bar\alpha\cdot\beta]$.

<a id="c-009"></a>

#### C-009 — Harmonic rotation of exponentials
**Local notation:** Here $F_t$ is the finite classical flow matrix; $\hat U_t$ is its Hilbert-space propagator.
**Source:** §1.3.2, Lemma 9, Proposition 3, printed 15–16; PDF 32–33; (1.42). For $\hat H_{\rm os}=\hbar(\hat N+\tfrac n2\hat I)$, $F_t=\left(\begin{smallmatrix}\cos t&\sin t\\-\sin t&\cos t\end{smallmatrix}\right)$,
$$
\hat U_t\hat T(z)\hat U_t^*=\hat T(F_tz),\qquad \hat U_t|\alpha\rangle=e^{-int/2}|e^{-it}\alpha\rangle.
$$
This exact covariance reduces products of harmonic exponentials and linear exponentials to a rotated displacement and a scalar phase.

<a id="c-010"></a>

#### C-010 — Resolution, trace, and Hilbert–Schmidt formulas
**Source:** §1.3.3, printed 17–20; PDF 34–37; Propositions 4–6, (1.43)–(1.59).
$$
(2\pi\hbar)^{-n}\int|\phi_z\rangle\langle\phi_z|\,dz=\hat I,
\quad\operatorname{Tr}\hat B=(2\pi\hbar)^{-n}\int\langle\phi_z,\hat B\phi_z\rangle dz,
$$
$$
\|\hat B\|_{\rm HS}^2=(2\pi\hbar)^{-n}\int\|\hat B\phi_z\|^2dz.
$$
The trace formula requires trace class, and the norm formula requires Hilbert–Schmidt class. A unitary displacement or metaplectic operator in infinite dimension is generally not trace class; its formal trace needs a distributional/regularized interpretation.

<a id="c-011"></a>

#### C-011 — Bargmann realization
**Source:** §1.4, printed 20–24; PDF 37–41; Propositions 7–10, (1.60)–(1.79). The Bargmann transform realizes creation as multiplication and annihilation as differentiation (the source's displayed $\hbar$-dependent formulas are inconsistent; use its declared $\hbar=1$ realization):
$$
\hat a^\dagger_j\mapsto\zeta_j,\quad \hat a_j\mapsto\partial_{\zeta_j},\quad
\hat N\mapsto\zeta\cdot\partial_\zeta.
$$
**Consequence at $\hbar=1$:** $e^{v\cdot \hat a}f(\zeta)=f(\zeta+v)$, $e^{u\cdot \hat a^\dagger}f(\zeta)=e^{u\cdot\zeta}f(\zeta)$, and $e^{s\hat N}f(\zeta)=f(e^s\zeta)$ on analytic vectors. The reproducing kernel is proportional to $e^{\bar\eta\cdot\zeta}$, with the book's measure normalization fixed by (1.77).

<a id="c-012"></a>

#### C-012 — Mehler kernel and analytic continuation
**Source:** §1.4, Theorem 1, Corollary 5, printed 24–25; PDF 41–42; (1.80), (1.85). At $\hbar=1$,
$$
\sum_{k\in\mathbb N^n}w^{|k|}\phi_k(x)\phi_k(y)
=\pi^{-n/2}(1-w^2)^{-n/2}
\exp\!\left[-\frac{1+w^2}{2(1-w^2)}(|x|^2+|y|^2)+\frac{2w\,x\cdot y}{1-w^2}\right],\quad |w|<1.
$$
For $n=1,\Re\tau>0$, the kernel of $e^{-\tau \hat H_{\rm os}}$ is $(2\pi\sinh\tau)^{-1/2}\exp[-(x^2+y^2)/(2\tanh\tau)+xy/\sinh\tau]$. The analytic square root is continued from real positive $\tau$; continuation to real-time caustics requires the metaplectic phase, not a pointwise principal square root.

<a id="c-013"></a>

#### C-013 — Weyl symbol and translation expansion
**Local notation:** Here $\hat A=\operatorname{Op}^w_\hbar(A)$ is the operator, while $A$ and $K$ are its scalar Weyl symbol and scalar integral kernel; the quantization map $\operatorname{Op}^w_\hbar$ keeps its standard name.
**Source:** §2.1, printed 28–30; PDF 45–47; Proposition 11–12, (2.2)–(2.11). The reliable kernel convention is
$$
[\hat A\psi](x)=(2\pi\hbar)^{-n}\int A((x+y)/2,\xi)e^{i(x-y)\cdot\xi/\hbar}\psi(y)\,dy\,d\xi,
$$
$$
A(x,\xi)=\int e^{-iu\cdot\xi/\hbar}K(x+u/2,x-u/2)\,du.
$$
The operator expands in Weyl translations and its translation coefficient is computed as $\operatorname{Tr}(\hat A\hat T(-X))$ under the smoothing hypotheses of Proposition 12. Source Fourier prefactors in (2.3)–(2.6) are not mutually consistent, so the kernel convention should determine any rescaling.

<a id="c-014"></a>

#### C-014 — Symbol traces and Gaussian overlaps
**Local notation:** Here $\hat A=\operatorname{Op}^w_\hbar(A)$ and $\hat B=\operatorname{Op}^w_\hbar(B)$; $A,B,W$ in the integrals are scalar symbols.
**Source:** §2.1–2.2, Proposition 13, 14, 17, printed 30, 33, 35–36; PDF 47, 50, 52–53; (2.12), (2.13), (2.21), (2.26).
$$
\|\hat A\|_{\rm HS}^2=(2\pi\hbar)^{-n}\int|A|^2,
\quad\operatorname{Tr}(\widehat A\widehat B)=(2\pi\hbar)^{-n}\int AB,
$$
$$
\langle\psi,\widehat A\phi\rangle=(2\pi\hbar)^{-n}\int A(X)W_{\psi,\phi}(X)dX,
\quad |\langle\psi,\phi\rangle|^2=(2\pi\hbar)^{-n}\int W_\psi W_\phi.
$$
The Wigner transform here is the Weyl symbol of the rank-one projector, not a normalized probability density.

<a id="c-015"></a>

#### C-015 — Exact symplectic covariance
**Local notation:** As in C-013, $\widehat{A\circ F}=\operatorname{Op}^w(A\circ F)$; $F$ is a finite symplectic matrix and $\hat R(F)$ its quantum implementer.
**Source:** §2.1.1, Theorem 2; §2.4.4, Theorem 13, Proposition 20; printed 31–32, 51–52; PDF 48–49, 68–69; (2.18), (2.93)–(2.95).
$$
\hat R(F)^*\hat A\hat R(F)=\widehat{A\circ F},\qquad
W_{\hat R(F)\psi,\hat R(F)\eta}(z)=W_{\psi,\eta}(F^{-1}z).
$$
For each real symplectic $F$, its unitary implementer is unique up to phase. Quadratic propagator paths reduce this freedom to the metaplectic double cover (C-023).

<a id="c-016"></a>

#### C-016 — Moyal product and exact quadratic commutator
**Local notation:** Here $A,B,H,E_t$ are scalar phase-space symbols; hats denote their Weyl quantizations.
**Source:** §2.4.1, printed 43–45; PDF 60–62; Theorem 7, (2.64)–(2.67); §2.4.3, Remark 9, printed 50; PDF 67; (2.89). The Weyl symbol of a product is the Moyal product. **Re-expression using ordinary derivatives and the kernel convention C-013:**
$$
A\star B=A\exp\!\left[\frac{i\hbar}{2}\left(\overleftarrow\partial_q\cdot\overrightarrow\partial_p-\overleftarrow\partial_p\cdot\overrightarrow\partial_q\right)\right]B.
$$
For Schwartz symbols the integral product is exact; the derivative expansion is generally an asymptotic expansion under the source's temperate-weight symbol hypotheses. If $H$ is a polynomial of degree at most two, higher odd derivatives vanish and $[\widehat H,\widehat B]=i\hbar\widehat{\{H,B\}}$ with the ordinary bracket used in this re-expression. The source's (2.89) uses a conflicting bracket sign; evaluate $H=q,B=p$ to fix conventions. **Consequence:** the symbol of $e^{t\widehat A}$ solves $\partial_t E_t=A\star E_t$, $E_0=1$, and need not equal the pointwise $e^{tA}$.

<a id="c-017"></a>

#### C-017 — Wick, anti-Wick, and smoothing
**Local notation:** The functions $A_c,A^c,A_w,G_\hbar$ are scalar symbols or scalar kernels, while $\hat A$ is the operator.
**Source:** §2.6.1, printed 56–58; PDF 73–75; Definition 11, Propositions 25–27, (2.111)–(2.123). Covariant/Wick symbol means $A_c(z)=\langle\phi_z,\widehat A\phi_z\rangle$; contravariant/anti-Wick means $\widehat A=\int A^c(z)|\phi_z\rangle\langle\phi_z|d\mu(z)$. **Consistent re-expression:** with $G_\hbar(X)=(\pi\hbar)^{-n}e^{-|X|^2/\hbar}$,
$$
A_w=G_\hbar*A^c=e^{\hbar\Delta/4}A^c,\qquad A_c=G_\hbar*A_w=e^{\hbar\Delta/4}A_w.
$$
The first convolution is explicit in (2.122); the second follows from C-014 and the coherent-state Wigner Gaussian. Inverting Gaussian smoothing is an existence/domain restriction: a bounded anti-Wick symbol forces strong analyticity of its Weyl symbol. The book does not justify replacing these maps by unrestricted convergent inverse heat series.

### Quadratic exponentials, metaplectic phases, and squeezing

<a id="c-018"></a>

#### C-018 — Quadratic evolution is a finite matrix problem
**Local notation:** Here $H_t(z)$ is the scalar classical quadratic; $S_t,J,F_t$ are finite coefficient matrices, and $I_n$ is the numerical $n$-by-$n$ identity.
**Source:** §3.1, printed 63–65; PDF 80–82; Theorem 15, (3.1)–(3.8). For $H_t(z)=\tfrac12z^TS_tz$, $S_t=S_t^T$ real and continuous,
$$
\dot F_t=JS_tF_t,\quad F_0=I_{2n},\quad J=\begin{pmatrix}0&I_n\\-I_n&0\end{pmatrix},
\qquad \hat U_t^*(\hat Q,\hat P)^T\hat U_t=F_t(\hat Q,\hat P)^T.
$$
Quantize the quadratic symmetrically with the same factor $1/2$: source (3.6) omits it. For a constant $S$, $F_t=e^{tJS}$. Time-dependent matrices require the matrix fundamental solution, not generally $e^{\int JS_tdt}$.

<a id="c-019"></a>

#### C-019 — Exact Gaussian propagation
**Local notation:** All blocks of $F_t$, $\Gamma,\Gamma_t$ are finite matrices, and $a_\Gamma$ is a scalar amplitude; $\hat U_t$ and $\hat T$ act on states.
**Source:** §3.2, Theorem 16, printed 65–71; PDF 82–88; (3.10), (3.19)–(3.21), (3.30)–(3.38). If $F_t=\left(\begin{smallmatrix}A_t&B_t\\C_t&D_t\end{smallmatrix}\right)$, $\Gamma=\Gamma^T$, $\Im\Gamma>0$, then
$$
\Gamma_t=(C_t+D_t\Gamma)(A_t+B_t\Gamma)^{-1},\qquad
\hat U_t\big[a_\Gamma e^{ix^T\Gamma x/(2\hbar)}\big]
=a_\Gamma\det(A_t+B_t\Gamma)^{-1/2}e^{ix^T\Gamma_tx/(2\hbar)}.
$$
The square root follows the path from $t=0$. A displaced Gaussian satisfies $\hat U_t\hat T(z)\phi^\Gamma=\hat T(F_tz)\hat U_t\phi^\Gamma$. For $\Gamma=iI_n$, $a_\Gamma=(\pi\hbar)^{-n/4}$. This expression is global even when the coordinate kernel has caustics.

<a id="c-020"></a>

#### C-020 — Siegel-domain protection and determinant cocycle
**Local notation:** Here $A,B,C,D,F,\Gamma$ are finite coefficient matrices; $\alpha$ is the scalar-matrix action on the Siegel domain and $\delta$ its scalar determinant cocycle.
**Source:** §3.2, Lemmas 20–21, Proposition 30, printed 69–71; PDF 86–88; (3.35)–(3.38). For $\Gamma$ in the positive Siegel domain, $A+B\Gamma$ is invertible and $\Gamma_F=(C+D\Gamma)(A+B\Gamma)^{-1}$ remains in that domain. **Re-expression:**
$$
\Im\Gamma_F=(A+B\Gamma)^{-*}(\Im\Gamma)(A+B\Gamma)^{-1}.
$$
The action is transitive and obeys $\alpha(F_1F_2)=\alpha(F_1)\alpha(F_2)$. Its scalar cocycle is
$\delta(F_1F_2,\Gamma)=\delta(F_1,\alpha(F_2)\Gamma)\delta(F_2,\Gamma)$, $\delta(F,\Gamma)=\det(A+B\Gamma)$. This explains determinant factors under Gaussian composition.

<a id="c-021"></a>

#### C-021 — Riccati and amplitude equations
**Local notation:** Here $S,K,L,G,A,B,C,D,M,N,\Gamma$ are finite matrices and $a$ is a scalar amplitude; this $N$ is unrelated to the number operator $\hat N$.
**Source:** §3.2, printed 67–69; PDF 84–86; Lemmas 18–19, (3.22)–(3.34). With a consistent block choice $S=\left(\begin{smallmatrix}K&L\\L^T&G\end{smallmatrix}\right)$,
$$
\dot\Gamma=-K-L\Gamma-\Gamma L^T-\Gamma G\Gamma,
\qquad\dot a=-\tfrac12\operatorname{Tr}(L+G\Gamma)a.
$$
Writing $M=A+iB,N=C+iD$ gives $\Gamma=NM^{-1}$; $d\log\det M/dt=\operatorname{Tr}(\dot MM^{-1})$ yields the determinant amplitude. The block names change between (3.2) and (3.28); always use the block choice displayed here.

<a id="c-022"></a>

#### C-022 — Coherent matrix elements without coordinate caustics
**Local notation:** Here $A,B,C,D,F,Y,Z,W$ are finite matrices, whereas $\hat U$ is the Hilbert-space propagator.
**Source:** §3.2, Theorem 17, printed 71–73; PDF 88–90; (3.39)–(3.41). At $\hbar=1$, put
$Y=A+D+i(B-C),\ Z=A-D+i(B+C),\ W=ZY^{-1}$, and $(x,\xi)=z-Fz'$. Then
$$
\langle\phi_z,\hat U\phi_{z'}\rangle=2^{n/2}\det(Y)^{-1/2}e^{-i\sigma(Fz',z)/2}
 e^{-(|x|^2+|\xi|^2)/4}e^{-(\xi+ix)^TW(\xi+ix)/4}.
$$
The determinant branch follows the generating symplectic path. In restoring $\hbar$, divide the phase-space quadratic exponents by $\hbar$. The exact source orientation of the two coherent labels must be kept; exchanging them conjugates phases.

<a id="c-023"></a>

#### C-023 — Metaplectic double cover and parity
**Source:** §3.3, printed 73–75; PDF 90–92; Lemma 24, Propositions 31–33, (3.42)–(3.48). Equivalent implementation conditions are covariance of all Weyl observables, covariance of Weyl translations, or covariance of the linear canonical operators. A symplectic path $F_0=I_{2n},F_1=F$ has quadratic generator $S_t=-J\dot F_tF_t^{-1}$. Two such paths yield implementers differing by a sign; products obey $\hat R(F_1)\hat R(F_2)=\pm \hat R(F_1F_2)$. The sign is tracked by the winding of $\det(A_t+iB_t)$. The even and odd subspaces give two inequivalent irreducible components. Arbitrary independent principal square roots destroy the composition law.

<a id="c-024"></a>

#### C-024 — Bogoliubov blocks
**Local notation:** The Bogoliubov blocks $u,v$ are finite matrices of scalar coefficients; $\hat U$ is the implementer and $\hat a$ the column of annihilation operators.
**Source:** §3.4, printed 76–77; PDF 93–94; (3.49), (3.50); §4.1.2, Lemma 28, printed 96; PDF 113. The finite symplectic action becomes linear mixing of annihilation and creation operators. **Re-expression in the order $(\hat a,\hat a^\dagger)$:**
$$
\hat U^*\hat a\hat U=u \hat a+v \hat a^\dagger,\qquad
u=\tfrac12[A+D+i(C-B)],\quad v=\tfrac12[A-D+i(C+B)],
$$
The canonical constraints are $uu^\dagger-vv^\dagger=I_n$ and $uv^T=vu^T$. This re-expression follows directly from C-018; it avoids transposition/conjugation ambiguities caused by the source's alternating $(\hat a^\dagger,\hat a)$ order and complex coordinate $(q-ip)/\sqrt2$.

<a id="c-025"></a>

#### C-025 — Squeeze generator and bounded squeeze coordinate
**Local notation:** Here $B,W,U$ are finite coefficient matrices, including the polar unitary matrix $U$; $\hat S(B)$ is the Fock-space squeeze operator.
**Source:** §3.4, printed 77–80; PDF 94–97; (3.58)–(3.64), Lemma 25. For a complex symmetric $B$,
$$
\hat S(B)=\exp\tfrac12(\hat a^\dagger\cdot B\hat a^\dagger-\hat a\cdot B^*\hat a),\qquad \hat S(B)^{-1}=\hat S(-B).
$$
If $W=U|W|$ is symmetric with $W^*W<I_n$, set
$$
B=U\operatorname{arctanh}|W|=W\sum_{k\ge0}\frac{(W^*W)^k}{2k+1},\qquad
W=U\tanh|B|.
$$
The strict inequality is essential: it ensures convergence and a normalizable squeezed vacuum. Equality permits infinite squeezing and falls outside this construction.

<a id="c-026"></a>

#### C-026 — Squeeze conjugation
**Local notation:** The displayed block array is a numerical matrix acting on the column of operators; every $I_n$ in that array is a coefficient identity.
**Source:** §3.4, Lemma 25, printed 79; PDF 96. With the notation in C-025,
$$
\hat S(B)\binom{\hat a}{\hat a^\dagger}\hat S(-B)=
\begin{pmatrix}
(I_n-WW^*)^{-1/2}&-W(I_n-W^*W)^{-1/2}\\
-(I_n-W^*W)^{-1/2}W^*&(I_n-W^*W)^{-1/2}
\end{pmatrix}\binom{\hat a}{\hat a^\dagger}.
$$
This is the displayed source result; its proof's intermediate derivatives contain conjugation/sign errors, so use the block result or differentiate the generator directly. For one mode $B=re^{i\theta}$, this reduces to $\hat S\hat a\hat S^{-1}=\hat a\cosh r-e^{i\theta}\hat a^\dagger\sinh r$.

<a id="c-027"></a>

#### C-027 — Squeezed vacuum from an operator Riccati equation
**Local notation:** The ansatz amplitude $a(t)$ is scalar and $M_t,B,W$ are finite matrices; the $\hat a^\dagger$ in the Fock vector are creation operators.
**Source:** §3.4, printed 80–82; PDF 97–99; (3.64)–(3.66). In Bargmann space at $\hbar=1$, the ansatz $a(t)e^{\zeta^TM_t\zeta/2}$ gives
$$
\dot M=B-MB^*M,\quad M_0=0,\qquad
\dot a=-\tfrac12a\operatorname{Tr}(B^*M),\quad M_1=W.
$$
**Re-expression as a Fock vector:**
$$
\hat S(B)|0\rangle=\det(I_n-W^*W)^{1/4}
\exp\!\left(\tfrac12\hat a^\dagger\cdot W\hat a^\dagger\right)|0\rangle,
\quad\langle0|\hat S(B)|0\rangle=\det(\cosh|B|)^{-1/2}.
$$
The source Bargmann measure contributes an additional $(2\pi)^{-n/2}$ to its function representative; this is not a physical vacuum normalization factor.

<a id="c-028"></a>

#### C-028 — Gaussian coordinate form of squeezed vacuum
**Local notation:** Here $\Gamma,W,B$ are finite Gaussian-parameter matrices and $\psi_B$ is a state wavefunction.
**Source:** §3.4, Proposition 36, printed 82–83; PDF 99–100; (3.67)–(3.70). Its reliable covariance is
$$
\Gamma=i(I_n-W)(I_n+W)^{-1}.
$$
**Corrected re-expression from the preceding Gaussian integral and C-027:**
$$
\psi_B(x)=(\pi\hbar)^{-n/4}\det(I_n-W^*W)^{1/4}\det(I_n+W)^{-1/2}
 e^{ix^T\Gamma x/(2\hbar)}.
$$
The branch is analytic from $W=0$. The printed amplitude in Proposition 36 is inconsistent with both $W=0$ and unit norm; it must not be quoted uncorrected.

<a id="c-029"></a>

#### C-029 — Squeeze–rotation factorization of a quadratic propagator
**Source:** §3.4, Propositions 34–35, printed 76–78; PDF 93–95; (3.51)–(3.62); §4.1.2, Proposition 38, printed 97; PDF 114. Polar decomposition of the Bogoliubov block separates the symplectic motion into a positive squeezing part and an orthogonal symplectic rotation. Correspondingly $\hat U_t=\hat S(B_t)\hat R_t$ with a continuous lift. The vacuum is an eigenvector of the number-preserving rotation, so $\hat U_t|0\rangle$ is a scalar phase times $\hat S(B_t)|0\rangle$, and $\hat U_t|\alpha\rangle$ is displaced to $F_tz$. Source definitions of the polar unitary vary by adjoint between chapters. **Convention-fixed phase reconstruction:** use the C-018 direction $\hat U_t^*\hat{\mathbf a}\hat U_t=u_t\hat{\mathbf a}+v_t\hat{\mathbf a}^{\dagger}$, put $Y_t=A_t+D_t+i(B_t-C_t)=2\bar u_t$, and choose the squeeze factor with positive vacuum amplitude $c_S=\det(u_tu_t^\dagger)^{-1/4}$. Then
$$
\hat R_t|0\rangle=\vartheta_t|0\rangle,\qquad
\vartheta_t=\frac{2^{n/2}\det(Y_t)^{-1/2}}{\det(u_tu_t^\dagger)^{-1/4}},\qquad |\vartheta_t|=1.
$$
The square root follows the symplectic path. This follows by taking the vacuum component of $\hat U_t=\hat S(B_t)\hat R_t$ and using C-022/C-027; it fixes the scalar in this declared convention. A metaplectic lift must still be specified; an isolated endpoint matrix does not fix its sign.

<a id="c-030"></a>

#### C-030 — Weyl symbol through a Cayley transform
**Local notation:** The Weyl symbol $R_w(F;X)$ is scalar; $F,J,S$ are finite matrices and $\hat R(F)$ is the operator represented by this symbol.
**Source:** §3.5, Theorem 18, printed 83–85; PDF 100–102; (3.72), (3.81). For $\det(I_{2n}+F)\ne0$, its Weyl symbol is Gaussian in the Cayley transform. **Corrected kernel-normalized re-expression:**
$$
R_w(F;X)=\frac{2^n e^{i\pi\nu}}{|\det(I_{2n}+F)|^{1/2}}
\exp\!\left[-\frac{i}{\hbar}X^TJ(I_{2n}-F)(I_{2n}+F)^{-1}X\right].
$$
The source omits $2^n$ in (3.72); the correction is forced by $F=I_{2n}$, whose Weyl symbol is 1. The path determines $\nu$ and square-root phase. For constant $S$ and $F=e^{tJS}$, (3.81) gives a determinant of $\cosh(tJS/2)$ and a quadratic form involving $\tanh(tJS/2)$. This formula is chartwise; $\det(I_{2n}+F)=0$ calls for another chart/distributional representation.

<a id="c-031"></a>

#### C-031 — Covariant translation symbol and phase tracking

**Local notation:** The translation coefficient is scalar; $F,J,Q_F$ are finite numerical matrices. **Source:** §3.5, Theorem 18(ii), Proposition 37, printed 83–85; PDF 100–102; (3.73)–(3.80). The following is an independently derived, convention-fixed reconstruction of the structural source formula.

Use $\sigma(z,X)=z^TJX$, the translations in C-003/C-004, and the Weyl symbol $R_w$ in C-030. Define
$$
\widetilde R(F;z)=\operatorname{Tr}[\hat T(-z)\hat R(F)]
=(2\pi\hbar)^{-n}\int R_w(F;X)e^{i\sigma(z,X)/\hbar}\,dX,
\qquad
\hat R(F)=(2\pi\hbar)^{-n}\int\widetilde R(F;z)\hat T(z)\,dz.
$$
These traces/integrals are oscillatory/distributional unless additional trace-class regularization is supplied. If both $I_{2n}+F$ and $I_{2n}-F$ are invertible, set
$$
Q_F=J(I_{2n}-F)(I_{2n}+F)^{-1}=Q_F^T,\qquad c_w=R_w(F;0).
$$
A real Fresnel integral then gives the complete chart formula
$$
\widetilde R(F;z)=\frac{c_w\,e^{-i\pi\operatorname{sgn}(Q_F)/4}}{2^n|\det Q_F|^{1/2}}
\exp\!\left[-\frac{i}{4\hbar}z^TJ(I_{2n}+F)(I_{2n}-F)^{-1}z\right].
$$
Here $\operatorname{sgn}(Q_F)$ is its positive-minus-negative eigenvalue count; the phase of $c_w$ is already fixed by the metaplectic path. The modulus of the prefactor is $|\det(I_{2n}-F)|^{-1/2}$. For one oscillator $F=e^{tJ}$, $0<t<\pi$, this becomes
$$
\widetilde R(F;z)=\frac{1}{2i\sin(t/2)}\exp\!\left[\frac{i}{4\hbar}\cot(t/2)|z|^2\right],
\qquad \hat R(F)=e^{-it(\hat N+\hat I/2)}.
$$
This resolves the local Fourier normalization and scalar phase once the stated transform and lift are fixed (GEN-13). If $I+F$ is singular while $I-F$ is invertible, continue the translation coefficient from an adjacent chart along that lift. If $I-F$ is singular, use a distributional chart or the coherent integral of Proposition 37; a divergent chart prefactor does not imply a singular unitary. The book's Fourier normalization and phase-index typography should not be transferred to a different convention without this reconstruction.

<a id="c-032"></a>

#### C-032 — Driven/affine quadratic propagators
**Local notation:** Here $H(z_s,s)$ and $\delta_t$ are scalars; the hatted factors are Hilbert-space operators, and $B_t$ is a finite squeeze-parameter matrix.
**Source:** §4.1.2, Proposition 39, printed 95–98; PDF 112–115; (4.6)–(4.12). For the second-order Taylor Hamiltonian along a classical orbit $z_t=(q_t,p_t)$, define
$$
\delta_t=\int_0^t[p_s\cdot\dot q_s-H(z_s,s)]ds-\tfrac12(q_t\cdot p_t-q\cdot p).
$$
Then its exact propagator is
$$
\hat U_2(t,s)=e^{i(\delta_t-\delta_s)/\hbar}\hat T(z_t)
\hat S(B_t)\hat R_t\hat R_s^{-1}\hat S(-B_s)\hat T(-z_s).
$$
For a genuinely affine quadratic Hamiltonian this is exact for the full evolution. For a nonlinear Hamiltonian it propagates the local quadratic approximation and has the controlled remainder described in C-033.

<a id="c-033"></a>

#### C-033 — What quadratic replacement does and does not prove
**Source:** §4.1.2–4.3.4, Theorems 21–26, printed 99–117; PDF 116–134; (4.13)–(4.74). A general smooth Hamiltonian evolves a coherent state into a displaced squeezed Gaussian times a hierarchy of polynomial corrections, not exactly a single Gaussian. The leading finite-time norm error is $O(\hbar^{1/2})$; arbitrary orders require higher Taylor terms. Errors grow with the classical stability matrix and are controlled up to suitable Ehrenfest times. Analyticity permits least-term truncation with $e^{-c/\hbar}$ errors on finite time intervals. These are asymptotic operator-exponential methods, not exact BCH identities for cubic or higher generators.

<a id="c-034"></a>

#### C-034 — Propagator reconstruction and Gaussian determinant branches
**Source:** §4.1, printed 92; PDF 109; (4.1); §5.3, printed 137–139; PDF 154–156; (5.37)–(5.43), Theorem 29. Resolution of identity reconstructs a propagator kernel from all propagated coherent states. In the ensuing Gaussian integrals, one determinant square root follows the symplectic path, while the Gaussian-integral branch is chosen by continuation from positive real part (eigenvalue square roots with positive real part). Their combination produces the Maslov index. This prescription is useful for quadratic exponential composition; it is not optional scalar decoration.

<a id="c-035"></a>

#### C-035 — Periodic quadratic evolution and quasi-energy
**Local notation:** The Floquet exponent $\rho$ here is a scalar, unrelated to a density operator $\hat\rho$. Differential operations retain their standard notation.
**Source:** §3.6, Theorems 19–20, Lemma 26, printed 86–90; PDF 103–107. In stable Mathieu/Paul-trap regions, the Floquet propagator has eigenvalues of the form $e^{-i\rho\pi(k+1/2)}$. Squeezed Hermite functions supply the eigenstates; square roots of classical complex solutions are followed continuously. The equivalence between Floquet eigenvectors and the extended operator $\hat K=-i\hbar\partial_t+\hat H_t$ is exact. This is an explicit nonautonomous quadratic-exponential application.

<a id="c-036"></a>

#### C-036 — Trace formulas require a spectral cutoff
**Local notation:** Here $P_\gamma$ is the finite reduced Poincaré matrix and $d_\gamma$ its dimension; $E$ is a scalar energy.
**Source:** §5.1–5.5, printed 127–156; PDF 144–173; Theorems 28–35, especially (5.16), (5.21)–(5.23). The trace of a propagator is used after energy localization/regularization, such as $\operatorname{Tr}[\chi(\hat H)\hat A\chi(\hat H)g((E\hat I-\hat H)/\hbar)]$. Periodic-orbit Gaussian determinants yield $|\det(I_{d_\gamma}-P_\gamma)|^{-1/2}$ and Maslov phases. The chapter does not license an ordinary trace of an arbitrary infinite-dimensional unitary $e^{-it\hat H}$. Its remaining calculations concern stationary-phase asymptotics and Bohr–Sommerfeld quantization, rather than additional exact creation/annihilation exponential identities.

### Compact and noncompact rank-one Lie algebras

<a id="c-037"></a>

#### C-037 — Torus Weyl exponentials and finite matrices
**Local notation:** The clock and shift matrices $\hat Z,\hat X$ are quantum operators on the $N$-dimensional torus Hilbert space; $N$ is its scalar dimension.
**Source:** §6.3–6.6, printed 164–183; PDF 180–199; Propositions 56–72, (6.4)–(6.53). Simultaneous quasiperiodicity under unit position/momentum translations requires $2\pi\hbar N=1$. Their finite-dimensional restrictions are clock and cyclic-shift matrices with $\hat X^N=\hat Z^N=\hat I$, and a root-of-unity commutation law. The explicit matrix definitions (6.18)–(6.19) fix the orientation: $\hat Z\hat X=\omega \hat X\hat Z$, $\omega=e^{2\pi i/N}$; the source's printed order in Proposition 59(iii) is reversed. Finite Weyl translations form an orthogonal basis of all operators in dimension $N$; symplectic automorphisms have exact Egorov covariance, subject to the boundary-phase compatibility condition (6.49). This is a compact-phase-space alternative, not a finite-dimensional representation of the unexponentiated CCR.

<a id="c-038"></a>

#### C-038 — SU(2) commutators, ladder matrices, and Casimir
**Source:** §7.3.1–7.3.2, printed 198–205; PDF 213–220; (7.11)–(7.22), Propositions 79–82.
$$
[\hat J_3,\hat J_\pm]=\pm \hat J_\pm,\quad[\hat J_+,\hat J_-]=2\hat J_3,\quad
\hat J_3|j,m\rangle=m|j,m\rangle,
$$
$$
\hat J_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}|j,m\pm1\rangle,
\quad \hat J^2=j(j+1)\hat I.
$$
Here $j\in\tfrac12\mathbb N$; the representation has dimension $2j+1$, and descends from SU(2) to SO(3) precisely for integral $j$. Thus every exponential of an element of this algebra can be reduced to a finite matrix group element and then represented at arbitrary spin.

<a id="c-039"></a>

#### C-039 — Fundamental 2-by-2 exponential
**Local notation:** Here $A$ and the Pauli basis matrices $\boldsymbol\sigma$ are numerical $2$-by-$2$ Lie-algebra matrices used to compute group parameters. The unit vector is $\mathbf n=\mathbf v/|\mathbf v|$, so its normalization is not denoted by an operator hat.
**Source:** §7.3.2, Lemma 41, Proposition 80, printed 201–202; PDF 216–217; (7.15). For $A\in\mathfrak{su}(2)$, $A^2=-(\det A)I_2$. If $\det A=1$,
$$
e^{tA}=I_2\cos t+A\sin t.
$$
**Consequence:** for $A=i\mathbf v\cdot\boldsymbol\sigma$, $e^{tA}=\cos(t|\mathbf v|)I_2+i\sin(t|\mathbf v|)\mathbf n\cdot\boldsymbol\sigma$. This is the concrete matrix method used to verify disentangling factors; $|\mathbf v|=0$ is handled by continuity.

<a id="c-040"></a>

#### C-040 — Gauss decomposition as a disentangling theorem
**Local notation:** Here $g$ and its three factors are numerical group matrices; their Hilbert-space representation has hatted generators as in C-041.
**Source:** §7.4.1, Lemma 44, printed 211–212; PDF 226–227; (7.36)–(7.38). For $g=\left(\begin{smallmatrix}\alpha&\beta\\\gamma&\delta\end{smallmatrix}\right)\in\mathrm{SL}(2,\mathbb C)$, $\delta\ne0$,
$$
g=\begin{pmatrix}1&\beta/\delta\\0&1\end{pmatrix}
\begin{pmatrix}\delta^{-1}&0\\0&\delta\end{pmatrix}
\begin{pmatrix}1&0\\\gamma/\delta&1\end{pmatrix}.
$$
The decomposition is unique. In a representation of the complexified Lie algebra it becomes a raising-exponential, diagonal-exponential, lowering-exponential factorization, with logarithm branches determined locally or by the representation. The chart restriction $\delta\ne0$ is essential; vanishing $\delta$ requires another Bruhat/pivot chart, possibly inserting a Weyl element. Simply reversing the triangular ordering is insufficient when both diagonal pivots vanish.

<a id="c-041"></a>

#### C-041 — SU(2) normal factorization, corrected explicit parameters
**Source location:** §7.4.1–7.4.2, printed 212–213; PDF 227–228; unnumbered normal form following (7.42). The source states the existence of $\hat G(\xi)=e^{\xi \hat J_+-\bar\xi \hat J_-}$ and a three-exponential normal form but prints inconsistent parameters. **Re-expression derived from C-039/C-040:** for $r=|\xi|$ in the chart $\cos r\ne0$,
$$
\hat G(\xi)=e^{\zeta \hat J_+}e^{-2\log(\cos r)\hat J_3}e^{-\bar\zeta \hat J_-},\qquad
\zeta=\frac{\xi}{r}\tan r.
$$
Initially take $|r|<\pi/2$ so the logarithm is real; other charts need a continuous branch. For the lowest-weight state $\hat J_-|j,-j\rangle=0$,
$\hat G(\xi)|j,-j\rangle=(1+|\zeta|^2)^{-j}e^{\zeta \hat J_+}|j,-j\rangle$.
The source alternates highest/lowest fiducials while retaining a raising exponential; specify the chosen fiducial before using the formula.

<a id="c-042"></a>

#### C-042 — SU(2) products, overlaps, and trace resolution
**Source:** §7.4.2, Propositions 84, 87–89, printed 213–223; PDF 228–238; (7.43)–(7.63). A product of chosen sphere displacements equals another displacement times a stabilizer rotation. In a consistent stereographic chart,
$$
\langle\eta|\zeta\rangle=
\frac{(1+\bar\eta\zeta)^{2j}}{(1+|\eta|^2)^j(1+|\zeta|^2)^j},\qquad
|\langle n'|n\rangle|^2=\left(\frac{1+n\cdot n'}2\right)^{2j},
$$
$$
\frac{2j+1}{\pi}\int_{\mathbb C}\frac{|\zeta\rangle\langle\zeta|}{(1+|\zeta|^2)^2}d^2\zeta=\hat I.
$$
Overlap and product phases are geometric area/Berry phases after a gauge/shortest-arc convention is fixed. The complex overlap is the least ambiguous computational formula; the source's intermediate spherical-area argument has factors/signs that should be checked against it.

<a id="c-043"></a>

#### C-043 — Spin symbols and exponential trace inequalities
**Local notation:** The symbols $H_c,H^c,(\hat J_\ell)_c,(\hat J_\ell)^c$ are scalar functions; hats within a symbol argument identify the underlying quantum operator.
**Source:** §7.6, Theorem 43, Propositions 91–92, printed 227–231; PDF 242–246; (7.75)–(7.99).
$$
\int\chi(H_c)d\mu\le\operatorname{Tr}\chi(\widehat H)\le\int\chi(H^c)d\mu
$$
for convex $\chi$, with covariant and contravariant symbols and a resolution measure. For a single spin, $(\hat J_\ell)_c=jn_\ell$, $(\hat J_\ell)^c=(j+1)n_\ell$. In particular $\chi(s)=e^{-\beta s}$ gives bounds on traces of exponentials. The measure resolves identity and therefore has total mass $\dim\mathcal H$, even though the discussion sometimes calls it a probability measure; normalized traces require dividing by dimension.

<a id="c-044"></a>

#### C-044 — Spin-to-CCR contraction
**Source:** §7.7, printed 232–234; PDF 247–249; (7.101)–(7.112), Propositions 93–94. At scale $\zeta/\sqrt{2j}$, the spin Bargmann coherent kernel converges to $e^{\bar\eta\zeta-|\eta|^2/2}$. Appropriately rescaled ladder operators contract to $\hat a^\dagger,\hat a$ as $j\to\infty$ with $2j\epsilon^2\to1$. This supplies an asymptotic comparison of compact-group exponentials and Weyl exponentials; equality at finite spin is not claimed.

<a id="c-045"></a>

#### C-045 — SU(1,1) algebra and analytic representation domains
**Local notation:** The hats also apply in the Bergman differential realization; $z$ is the scalar coordinate and $\partial_z$ retains its standard differential notation.
**Source:** §8.2, printed 243–249; PDF 258–264; (8.11), (8.35)–(8.43), Lemma 52, Theorem 44. Use
$$
[\hat K_-,\hat K_+]=2\hat K_0,\qquad[\hat K_0,\hat K_\pm]=\pm \hat K_\pm,\qquad
\hat C=\hat K_0^2-\tfrac12(\hat K_+\hat K_-+\hat K_-\hat K_+)=k(k-1)\hat I.
$$
The holomorphic discrete series is realized in weighted Bergman spaces on $|z|<1$; $\hat K_-=\partial_z$, $\hat K_+=2kz\hat I+z^2\partial_z$, $\hat K_0=k\hat I+z\partial_z$. Polynomial vectors form a core; the closed ladders are adjoints. The source constructs genuine SU(1,1) representations for integral $2k\ge2$; other admissible $k$ may belong to covering-group representations rather than the group itself.

<a id="c-046"></a>

#### C-046 — Bosonic realizations of SU(1,1)
**Local notation:** The number difference $\nu$ labels a sector and is scalar; $\hat N_1-\hat N_2$ is its operator before that restriction.
**Source:** §8.2.6, printed 253–255; PDF 268–270; (8.45)–(8.50), Lemma 53.
$$
\hat K_+=\tfrac12(\hat a^\dagger)^2,\quad \hat K_-=\tfrac12\hat a^2,\quad \hat K_0=\tfrac12(\hat N+\tfrac12\hat I).
$$
Even and odd oscillator sectors have Bargmann indices $k=1/4,3/4$, both with Casimir $-3/16$, and give covering/projective representations. The two-mode realization is
$$
\hat K_+=\hat a_1^\dagger \hat a_2^\dagger,\quad \hat K_-=\hat a_1\hat a_2,\quad
\hat K_0=\tfrac12(\hat N_1+\hat N_2+\hat I),\quad
\hat C=\tfrac14[(\hat N_1-\hat N_2)^2-\hat I].
$$
On a fixed number-difference sector, $\hat N_1-\hat N_2=\nu\hat I$ and $k=(1+|\nu|)/2$. These realizations turn pair creation/annihilation exponentials into rank-one group exponentials.

<a id="c-047"></a>

#### C-047 — SU(1,1) coherent series and overlap
**Local notation:** The function $\Gamma(x)$ in the series is the scalar Gamma function, distinct from the Fock second-quantization operator in C-058.
**Source:** §8.3.1–8.3.2, printed 255–260; PDF 270–275; (8.51)–(8.81), Propositions 101–105. For a lowest-weight vector $\hat K_-|k,0\rangle=0$,
$$
|\zeta;k\rangle=(1-|\zeta|^2)^k e^{\zeta \hat K_+}|k,0\rangle
=(1-|\zeta|^2)^k\sum_{m\ge0}\sqrt{\frac{\Gamma(2k+m)}{\Gamma(2k)m!}}\zeta^m|k,m\rangle,
\quad|\zeta|<1,
$$
$$
\langle\eta;k|\zeta;k\rangle=
\frac{(1-|\eta|^2)^k(1-|\zeta|^2)^k}{(1-\bar\eta\zeta)^{2k}}.
$$
These are consistent re-expressions of the source, which changes conjugate-label conventions and misprints `sinh` for `tanh` in (8.51). The identity resolution uses $(2k-1)d^2\zeta/[\pi(1-|\zeta|^2)^2]$, requiring $k>1/2$ for this measure formula. It does not apply directly to the even oscillator sector $k=1/4$; the odd sector $k=3/4$ satisfies the bound and can use this resolution.

<a id="c-048"></a>

#### C-048 — Noncompact Gauss factorization and squeeze disentangling
**Source location:** §8.3.1, printed 256; PDF 271; (8.53)–(8.57), with Gauss decomposition (8.55); §8.5.1, printed 263; PDF 278. **Re-expression derived with the algebra of C-045:**
$$
e^{\xi \hat K_+-\bar\xi \hat K_-}
=e^{\zeta \hat K_+}e^{-2\log(\cosh r)\hat K_0}e^{-\bar\zeta \hat K_-},\qquad
r=|\xi|,\quad\zeta=(\xi/r)\tanh r.
$$
For one-mode pairs (C-046),
$$
\hat S(\xi)=e^{\zeta(\hat a^\dagger)^2/2}
 e^{-\log(\cosh r)(\hat N+\tfrac12\hat I)}e^{-\bar\zeta \hat a^2/2}.
$$
For two-mode pairs replace $(\hat a^\dagger)^2/2,\hat a^2/2,\hat N+\tfrac12\hat I$ by $\hat a_1^\dagger \hat a_2^\dagger,\hat a_1\hat a_2,\hat N_1+\hat N_2+\hat I$. These are normal factorization formulas, justified first on analytic vectors; the full unitary exponential has the closed extension.

<a id="c-049"></a>

#### C-049 — One-mode squeezing and Wigner ellipse
**Local notation:** Here $\xi,\omega,\delta$ and the one-mode Gaussian coordinate $W$ are scalars; the Wigner function $W(q,p)$ is also scalar.
**Source:** §8.5.1, Lemma 57, Propositions 108–109, printed 263–265; PDF 278–280. If $\xi=re^{i\theta}$, $\omega=e^{i\theta}\tanh r$,
$$
\hat S\hat a\hat S^{-1}=\frac{\hat a-\omega \hat a^\dagger}{\sqrt{1-|\omega|^2}},\quad
\hat S\hat K_0\hat S^{-1}=\cosh(2r)\hat K_0-\tfrac12\sinh(2r)(e^{i\theta}\hat K_++e^{-i\theta}\hat K_-).
$$
With $\delta=(1-\omega)/(1+\omega)$, the normalized coordinate Gaussian has exponent $-\delta x^2/(2\hbar)$. Here $\Re\delta>0$, and the Wigner symbol is
$$
W(q,p)=2\exp[-q^2\Re\delta/\hbar-(p+q\Im\delta)^2/(\hbar\Re\delta)].
$$
Rendered inspection confirms that the source correctly prints $\Re\delta>0$ and the Wigner formula above; extraction alone had obscured the real/imaginary-part glyphs. For the coordinate-wavefunction phase use the continuous-lift convention of C-028 with scalar $W=\omega$.

<a id="c-050"></a>

#### C-050 — Time-dependent SU(1,1) Hamiltonian
**Source:** §8.5.2, Proposition 110, printed 265–266; PDF 280–281; (8.95)–(8.101). For $\hat H_t=\lambda_t\hat K_++\bar\lambda_t\hat K_-+\mu_t\hat K_0$, $\mu_t$ real,
$$
i\dot\omega=\bar\lambda\omega^2+\mu\omega+\lambda,\quad\omega_0=0,\qquad
\dot\gamma=-\lambda\bar\omega-\bar\lambda\omega-\mu,\quad\gamma_0=0.
$$
**Corrected two-time form:** $\hat U(t,s)=\hat S(\beta_t)e^{i(\gamma_t-\gamma_s)\hat K_0}\hat S(-\beta_s)$, with $\omega_t=(\beta_t/|\beta_t|)\tanh|\beta_t|$. Source (8.96) incorrectly ends in $\hat S(-\beta_t)$, failing the group composition rule; the corrected endpoint agrees with Proposition 38. The disk coordinate remains $|\omega_t|<1$ for real-time self-adjoint dynamics.

<a id="c-051"></a>

#### C-051 — Derivative of a noncommuting exponential
**Source:** §8.5.2 proof of Proposition 110, printed 266; PDF 281. Duhamel's differentiation formula is
$$
\frac d{dt}e^{\hat L(t)}=\int_0^1e^{s\hat L(t)}\dot{\hat L}(t)e^{(1-s)\hat L(t)}ds.
$$
This replaces the generally false $\dot{\hat L} e^{\hat L}$ when $[\hat L,\dot{\hat L}]\ne0$. In a finite-dimensional closed Lie algebra, the adjoint action under the integral is computed by a finite matrix exponential; the source uses it to derive the squeeze-coordinate equations. For unbounded generators, invariant domains and the existence of the differentiable propagator remain hypotheses.

<a id="c-052"></a>

#### C-052 — Singular quadratic-looking generators require a boundary condition
**Source:** §8.5.2, printed 266–269; PDF 281–284; (8.102)–(8.107), Corollary 28, Lemma 58. The oscillator plus $g^2/(2\hat Q^2)$ again closes an SU(1,1) algebra, but on $L^2(\mathbb R_+)$. The construction uses the quadratic form domain $H^1_0(\mathbb R_+)\cap\{u:xu\in L^2\}$, Hardy's inequality, and the resulting self-adjoint realization. Its lowest-weight wavefunction is proportional to $x^s e^{-x^2/2}$, $s=1/2+\sqrt{1/4+g^2}$, and $k=(2s+1)/4$. The group identities do not choose a self-adjoint extension automatically; that choice is part of the operator problem.

### Additional Gaussian and infinite-boson results

<a id="c-053"></a>

#### C-053 — A consistent oscillator Weyl exponential
**Source:** §10.1.1, printed 300; PDF 313; (10.5). At $\hbar=1$, the Weyl symbol of $e^{it\hat H_{\rm os}}$, $|t|<\pi$, is
$$
(\cos(t/2))^{-n}\exp[i\tan(t/2)|X|^2].
$$
This agrees with C-030 and resolves the incorrect power in (2.15). **Consequence:** for positive real $\tau$, $[e^{-\tau \hat H_{\rm os}}]_w=(\cosh(\tau/2))^{-n}\exp[-\tanh(\tau/2)|X|^2]$. The kernel-normalized formula has determinant power $-n$, not $-n/2$. Sections 10.1.2–10.2 add commutator/matrix criteria and functional calculus, rather than new finite Lie factorizations.

<a id="c-054"></a>

#### C-054 — Exact Herman–Kluk representation for quadratic exponentials
**Local notation:** The prefactor $C_t$ is scalar; $A_t,B_t,C_t^{\rm block},D_t,F_t$ are finite classical matrices, while $\hat U_t$ acts on states.
**Source:** §11.3, Theorems 56–58 and quadratic case (11.54), printed 324–329; PDF 337–342. For a quadratic Hamiltonian the coherent-state phase-space integral is exact with no higher-order semiclassical correction. **Normalized re-expression for a homogeneous quadratic propagator, using C-010 conventions:**
$$
\hat U_t=(2\pi\hbar)^{-n}\int C_t\,|\phi_{F_tz}\rangle\langle\phi_z|\,dz,
\quad C_t=\det\!\left[\frac{A_t+D_t+i(C_t^{\rm block}-B_t)}2\right]^{1/2}.
$$
Here $C_t^{\rm block}$ denotes the lower-left block (to distinguish it from the scalar prefactor); the root starts at 1. For affine quadratics multiply the integrand by the action phase in C-032. For nonlinear subquadratic Hamiltonians this is a leading approximation with transport corrections. The source's several printed amplitude conventions in (11.40), (11.50), and (11.92) do not agree at $t=0$; the normalized expression here is explicitly a re-expression with the identity limit fixed.

<a id="c-055"></a>

#### C-055 — Gaussian density operators and thermal factors
**Local notation:** The covariances $V,V_j$ and the symplectic matrix $J$ are finite coefficient matrices; $\hat\rho,\hat H,\hat N_j$ are quantum operators.
**Source:** §13.3.2, Proposition 135 and Remark 80, printed 385–386; PDF 396–397; (13.28)–(13.30). At $\hbar=1$, a centered Gaussian Weyl symbol with covariance $V$ represents a density operator exactly when $V+iJ/2\ge0$, equivalently all symplectic eigenvalues are at least $1/2$. It is pure exactly when $2V$ is positive symplectic, and
$$
\operatorname{Tr}\hat\rho^2=2^{-n}(\det V)^{-1/2}.
$$
Up to a metaplectic conjugation it is a product of one-mode thermal states. **Consequence from the oscillator spectrum:** for $\hat H=\sum_j\hbar\omega_j(\hat N_j+\tfrac12\hat I)$, $\omega_j>0,\beta>0$,
$$
\operatorname{Tr}e^{-\beta \hat H}=\prod_j[2\sinh(\beta\hbar\omega_j/2)]^{-1},\qquad
V_j=\frac{\hbar}{2}\coth(\beta\hbar\omega_j/2)I_2
$$
in dimensionless canonical oscillator coordinates. This explicitly states the coordinates rather than copying the inconsistent frequency scalings in (13.34)/(13.66).

<a id="c-056"></a>

#### C-056 — Partial trace of a quadratic exponential evolution
**Local notation:** Here $X,Y,V_E,V_S$ are finite matrices and $\chi$ is a scalar characteristic function; the reduced quantum map is described by these scalar data.
**Source:** §13.1.2 and §13.4.1–13.4.2, printed 377, 392–399; PDF 388, 403–410; (13.10), Theorems 73–74, (13.64), (13.65). A partial trace is phase-space integration of the Weyl symbol with factor $(2\pi\hbar)^{-N}$. A joint quadratic propagator preserves Gaussianity, and its reduced action is transport followed by Gaussian noise. **Consistent Schrödinger-picture re-expression:** if the joint canonical transformation has system output $z_{\rm out}=Xz_{\rm in}+Yu_{\rm in}$, the initial state factors, and the environment is centered Gaussian with covariance $V_E$,
$$
\chi_{S,\rm out}(k)=\chi_{S,\rm in}(X^Tk)\exp[-\tfrac12k^TYV_EY^Tk].
$$
Thus $V_{S,\rm out}=XV_SX^T+YV_EY^T$. These formulas remain meaningful when $X$ is singular; the source's time-local master equation requires invertible system block up to its critical time. The full joint unitary is regular there; only the reduced inverse/chart can fail.

<a id="c-057"></a>

#### C-057 — Segal fields and a rigorous exponential domain
**Source:** §15.2–15.3, Lemmas 92–94, printed 433–439; PDF 443–449; (15.7)–(15.21). On bosonic Fock space over a one-particle Hilbert space $\mathfrak h$, $\hat a(f)$ is antilinear in $f$, $\hat a^\dagger(f)$ linear, and
$$
[\hat a(f),\hat a^\dagger(g)]=\langle f,g\rangle \hat I,\qquad
\hat\Phi(f)=\frac{\hat a(f)+\hat a^\dagger(f)}{\sqrt2}.
$$
Finite-particle vectors are analytic for $\hat\Phi(f)$, proving its essential self-adjointness by Nelson's criterion. Creation and annihilation are controlled by $(\hat N+\hat I)^{1/2}$. This supplies a rigorous foundation for infinite-mode Weyl exponentials with $f\in\mathfrak h$; point fields $\hat a(x)$ are distributions, not individual bounded operators.

<a id="c-058"></a>

#### C-058 — Second quantization commutes with one-particle exponentiation
**Local notation:** Here $\hat h$ and $\hat U$ are operators on the one-particle Hilbert space, $h_{jk}$ are matrix coefficients of $\hat h$, and $\widehat\Gamma(\cdot),\widehat{d\Gamma}(\cdot)$ are operators on Fock space. The vector $\Psi(f)$ is not an operator.
**Source:** §15.2.1, Lemma 92, Remark 93, printed 433; PDF 443; (15.6).
$$
\widehat\Gamma(e^{-it \hat h})=e^{-it\,\widehat{d\Gamma}(\hat h)},\qquad \widehat{d\Gamma}(\hat I_{\mathfrak h})=\hat N.
$$
The left side acts by tensor powers of the one-particle unitary in each symmetric or antisymmetric sector. **Consequence:** number-preserving quadratic exponentials $e^{-it\sum h_{jk}\hat a_j^\dagger \hat a_k}$ are computed from the one-particle matrix exponential; on coherent vectors, $\widehat\Gamma(\hat U)\Psi(f)=\Psi(\hat U f)$. This follows from the exponential-vector series C-059 and applies to either finite or infinite one-particle dimension when the stated self-adjoint/unitary realizations exist.

<a id="c-059"></a>

#### C-059 — Infinite-mode displacement, exponential vector, and Poisson law
**Source:** §15.3, Definition 38, Propositions 139–140, printed 439–441; PDF 449–451; (15.22)–(15.24).
$$
\hat T(f)=e^{\hat a^\dagger(f)-\hat a(f)}=e^{-\|f\|^2/2}e^{\hat a^\dagger(f)}e^{-\hat a(f)},\qquad
\Psi(f)=e^{-\|f\|^2/2}\bigoplus_{k\ge0}\frac{f^{\otimes k}}{\sqrt{k!}}.
$$
$$
\hat T(f)\hat T(g)=e^{-i\Im\langle f,g\rangle}\hat T(f+g),\qquad
\langle\Psi(f),\Psi(g)\rangle=e^{-\|f-g\|^2/2+i\Im\langle f,g\rangle}.
$$
The occupation-number distribution is Poisson with mean and variance $\|f\|^2$. Use (15.24) for the overlap phase: the proposition's earlier inline overlap has a misparenthesized factor $1/2$.

<a id="c-060"></a>

#### C-060 — Weyl shifts preserve number-domain scales
**Source:** §15.3, Proposition 140(iii), Lemma 95, printed 440–442; PDF 450–452; (15.25).
$$
\hat T(f)^*\hat a(g)\hat T(f)=\hat a(g)+\langle g,f\rangle\hat I,
\quad \hat T(f)^*\hat N\hat T(f)=\hat N+\hat a(f)+\hat a^\dagger(f)+\|f\|^2\hat I.
$$
For every real $r$, $\langle \hat N\rangle^r\hat T(f)\langle \hat N\rangle^{-r}$ extends boundedly. This makes displacement conjugations and differentiation meaningful on weighted Fock domains, not just on formal polynomials.

<a id="c-061"></a>

#### C-061 — Derivative of a moving infinite-mode displacement
**Source:** §15.3, Lemma 96, printed 442; PDF 452; (15.26). **Consistent re-expression with the antilinearity of $\hat a(f)$ explicit:**
$$
\frac d{dt}\hat T(f_t)=\hat T(f_t)\big[\hat a^\dagger(\dot f_t)-\hat a(\dot f_t)+i\Im\langle f_t,\dot f_t\rangle \hat I\big].
$$
For $f_t\in C^1(\mathbb R,\mathfrak h)$, the derivative is strong from $D(\hat N^{1/2})$ to Fock space. The scalar term is the phase correction needed when shifting a time-dependent Hamiltonian to a moving classical field.

<a id="c-062"></a>

#### C-062 — Quadratic fluctuations about an interacting boson field
**Source:** §15.4.2, Propositions 141–144, printed 444–449; PDF 454–459; (15.28)–(15.49). A displacement by $\hbar^{-1/2}\varphi_t$ around a Hartree solution removes the linear field term. The transformed generator is a quadratic fluctuation Hamiltonian plus $\hbar^{1/2}\hat H_3+\hbar \hat H_4$; its leading evolution is a Bogoliubov propagator. The strong convergence $\hat W(t,s)\to \hat U_2(t,s)$ requires the source's Hardy-potential and regularity assumptions. It is not an exact reduction of a quartic field Hamiltonian to a Gaussian exponential.

<a id="c-063"></a>

#### C-063 — Dyson expansion versus convergent exponential factorization
**Source:** §15.4.3, Theorems 82–83, printed 450–453; PDF 460–463; (15.51)–(15.57). Iterated Duhamel gives an expansion in $\kappa=\hbar^{1/2}$ for the fluctuation evolution. For bounded stable interactions and vectors in $D(e^{\beta \hat N})$, short-time analyticity and Gevrey bounds yield the stated 2-Borel summability. With singular interactions the book only asserts an asymptotic Gevrey expansion and explicitly leaves Borel summability open. This is the relevant extension beyond closed quadratic Lie algebras; a formal Dyson series alone does not prove convergence.

### Fermionic exponentials and Grassmann calculus

<a id="c-064"></a>

#### C-064 — CAR, finite dimension, and bounded generators
**Source:** §16.2, Proposition 145–146, printed 458–461; PDF 468–471; (16.1), (16.2), (16.6)–(16.11).
$$
\{\hat c_j,\hat c_k^\dagger\}=\delta_{jk}\hat I,\quad\{\hat c_j,\hat c_k\}=\{\hat c_j^\dagger,\hat c_k^\dagger\}=0,
\qquad\|\hat c(f)\|=\|\hat c^\dagger(f)\|=\|f\|.
$$
For $n$ modes the fermionic Fock space has dimension $2^n$; Jordan–Wigner strings of Pauli matrices give an explicit matrix realization. Exponentials with ordinary complex coefficients are ordinary finite-matrix exponentials. Grassmann labels introduce an additional supermodule and are a different construction.

<a id="c-065"></a>

#### C-065 — Ordinary one-mode fermionic linear exponential
**Consequence of the explicit CAR matrix model in §16.2, printed 458–459; PDF 468–469 (C-064); no separately numbered source formula.** For ordinary complex $u,v$, $\hat L=u\hat c^\dagger+v\hat c$ satisfies $\hat L^2=uv\hat I$, so
$$
e^{\hat L}=\cosh\sqrt{uv}\,\hat I+\frac{\sinh\sqrt{uv}}{\sqrt{uv}}(u\hat c^\dagger+v\hat c).
$$
Interpret both coefficients through their entire power series, including $uv=0$, so no branch choice remains. The same consequence holds for finitely many modes with $\hat L=u\cdot \hat c^\dagger+v\cdot \hat c$ and $uv$ replaced by the bilinear $u\cdot v$: CAR cancel all cross terms and again give $\hat L^2=(u\cdot v)\hat I$. A bosonic central-commutator normal-ordering factor cannot be copied: $[\hat c,\hat c^\dagger]=\hat I-2\hat c^\dagger \hat c$ is not central. This distinction is essential when the requested exponential has complex rather than Grassmann coefficients.

<a id="c-066"></a>

#### C-066 — Nilpotent functions, parity, and derivative orientation
**Local notation:** Grassmann generators and parameters are anticommuting scalars, not Fock-space operators; their functions and derivatives therefore have no operator hats.
**Source:** §16.2–16.3.2, printed 459–464; PDF 469–474; (16.4)–(16.19), Lemmas 100–102. Grassmann generators anticommute and square to zero; ordered monomials give a $2^n$-dimensional algebra. A Grassmann number is invertible exactly when its scalar body is nonzero. Functions of a nilpotent part are finite Taylor polynomials; for a nonzero body expand about that body. Even elements commute, odd elements anticommute, and
$$
\partial_\theta(fg)=(\partial_\theta f)g+(-1)^{\pi(f)}f(\partial_\theta g),\quad
f\overleftarrow\partial_\theta=(-1)^{1+\pi(f)}\overrightarrow\partial_\theta f.
$$
Do not replace a left derivative by a right derivative without the parity sign.

<a id="c-067"></a>

#### C-067 — Berezin integration is ordered coefficient extraction

**Local notation:** $A$ is a finite even scalar matrix, $\theta$ a column of Grassmann scalars, and $\widehat I_\theta$ a linear integration map on their algebra. **Source:** §16.3.2, Definition43, Lemma102, printed464–465 / PDF474–475, (16.17)–(16.19). One-generator integration is coefficient extraction: $\int1\,d\theta=0$, $\int\theta\,d\theta=1$. Under an even invertible linear map,
$$
\int f(A\theta)\,d\theta=\det A\int f(\theta)\,d\theta.
$$
The transformed odd integration element therefore carries an inverse determinant.

**Source inconsistency and declared repair:** with explicitly left derivatives, the source order (16.17) and Gaussian normalization (16.20) have opposite one-pair signs. To agree with the Gaussian/CAR formulas below, fix
$$
\widehat I_\theta\!\left[\prod_{j=1}^n\theta_j\theta_j^*\right]=1,
\qquad\widehat I_\theta[e^{-\theta^*\theta}]=1,
\qquad
\widehat I_\theta f=\partial_{\theta_1^*}\partial_{\theta_1}\cdots\partial_{\theta_n^*}\partial_{\theta_n}f.
$$
All derivatives in this composition are left derivatives and act on the expression to their right. This convention is attached here to the written paired measure $d\theta_1^*d\theta_1\cdots d\theta_n^*d\theta_n$; it reverses the literal source (16.17) result for an odd number of pairs. The paired order, rather than an ambiguous compact measure alone, fixes every subsequent sign.

<a id="c-068"></a>

#### C-068 — Fermionic Gaussian determinant and sources

**Local notation:** $B$ is a finite coefficient matrix; $\theta,\gamma$ are Grassmann scalars. **Source:** §16.3.3, Propositions147–148, printed465–466 / PDF475–476, (16.20)–(16.21). With the explicitly repaired paired convention C-067,
$$
\widehat I_\theta e^{-\theta^*B\theta}=\det B,
\qquad
\widehat I_\theta e^{-\theta^*B\theta+\gamma^*\theta+\theta^*\gamma}
=\det B\,e^{\gamma^*B^{-1}\gamma}.
$$
Odd sources anticommute with every integrated generator. The source matrix is Hermitian; the finite polynomial calculation extends the first identity to arbitrary complex and singular $B$. The second displayed inverse requires invertibility. The preceding source's derivative order is not silently assumed compatible: C-067 records the one-mode contradiction and fixes it. Exact exterior-algebra checks verify every coefficient for complex noncommuting one-, two- and three-mode examples.

<a id="c-069"></a>

#### C-069 — Pfaffian Gaussian integral

**Local notation:** $\Phi$ is a Grassmann scalar; $K,L,M,\Lambda$ are finite coefficient matrices. **Source:** §16.3.3, Proposition149, printed466 / PDF476, (16.22). Fix block coordinate order
$$
x=(\gamma_1,\ldots,\gamma_n,\gamma_1^*,\ldots,\gamma_n^*)^T,
\quad\Phi=\tfrac12x^T\Lambda x,
\quad\Lambda=\begin{pmatrix}K&-M^T\\M&L\end{pmatrix},\quad K^T=-K,\ L^T=-L.
$$
In terms of the **standard Pfaffian in this block order**, the corrected formula with the paired integration C-067 is
$$
\widehat I_\gamma e^{\Phi+\gamma^*\xi-\xi^*\gamma}
=(-1)^{n(n-1)/2}\operatorname{Pf}(\Lambda)
\exp\!\left[\tfrac12\binom{\xi^*}{\xi}^{\!T}\Lambda^{-1}\binom{\xi^*}{\xi}\right].
$$
The inverse needs nonsingularity; the integral itself is a finite polynomial. The source's Pfaffian notation does not include this explicit block-to-paired sign: either include it or define a paired Pfaffian that absorbs it. For $n=2$, $\Phi=-\gamma^*\gamma$, the normalized integral is $+1$ and the standard block Pfaffian is $-1$. In all cases $(\operatorname{Pf}\Lambda)^2=\det\Lambda$. A determinant square never determines the ordered Pfaffian sign.

<a id="c-070"></a>

#### C-070 — Fermionic Bargmann model and normal monomials
**Local notation:** The realizations are $\hat c_j^\dagger\leftrightarrow\theta_j$ and $\hat c_j\leftrightarrow\partial_{\theta_j}$; $A_{\epsilon\epsilon'}$ are scalar coefficients, possibly Grassmann-valued.
**Source:** §16.4, printed 467–470; PDF 477–480; Propositions 150–151, (16.25)–(16.31). Fermionic creation is left multiplication by $\theta_j$, annihilation is $\partial_{\theta_j}$; the holomorphic monomials are orthonormal for the Berezin Gaussian inner product with its stated measure. Every operator has a unique finite expansion
$$
\widehat A=\sum_{\epsilon,\epsilon'\in\{0,1\}^n}A_{\epsilon\epsilon'}(\hat c^\dagger)^\epsilon \hat c^{\epsilon'}.
$$
For Grassmann coefficients their position relative to odd operators matters. The book's sentence calling (16.30) a representation with annihilators first contradicts its displayed monomial; use the actual order shown.

<a id="c-071"></a>

#### C-071 — Grassmann Weyl translation and central BCH
**Local notation:** The labels $\alpha,\gamma$ are anticommuting scalars; $\hat c,\hat c^\dagger,\hat T$ act on the Grassmann-extended Fock module, whose identity is $\hat I$.
**Source:** §16.5.1, Lemma 104, printed 470–472; PDF 480–482; (16.32)–(16.36). With independent odd labels anticommuting with odd operators,
$$
\hat T(\gamma)=e^{\hat c^\dagger\cdot\gamma-\gamma^*\cdot \hat c},\qquad
\hat T(\gamma)=e^{-\gamma^*\cdot\gamma/2}e^{\hat c^\dagger\cdot\gamma}e^{-\gamma^*\cdot \hat c}.
$$
The combined generators are **even**; their ordinary commutator is central. Thus the ordinary central BCH formula applies in the enlarged Grassmann algebra. Products obey
$$
\hat T(\alpha)\hat T(\gamma)=\hat T(\alpha+\gamma)e^{-(\alpha^*\cdot\gamma+\alpha\cdot\gamma^*)/2},\quad
\hat T(\gamma)^{-1}=\hat T(-\gamma)=\hat T(\gamma)^*,
$$
$\hat T^*\hat c\hat T=\hat c+\gamma\hat I$, $\hat T^*\hat c^\dagger \hat T=\hat c^\dagger+\gamma^*\hat I$. These are Grassmann identities, unlike C-065.

<a id="c-072"></a>

#### C-072 — Fermionic coherent state, overlap, and identity resolution
**Source:** §16.5.2, Proposition 152, printed 472–474; PDF 482–484; (16.37)–(16.47).
$$
\psi_\gamma(\theta)=e^{-\gamma^*\gamma/2}e^{\theta\cdot\gamma},\qquad
\hat c_j\psi_\gamma=\gamma_j\psi_\gamma,
$$
$$
\langle\psi_\gamma,\psi_\alpha\rangle=
\exp[\gamma^*\alpha-(\gamma^*\gamma+\alpha^*\alpha)/2],\qquad
\psi=\int\langle\psi_\gamma,\psi\rangle\psi_\gamma\,d\gamma^*d\gamma.
$$
The order in the last integral is the source's order; Grassmann coefficients cannot be commuted freely through the vectors. Starting from the filled state produces creation-eigenvector coherent states with reversed measure order, (16.45)–(16.47).

<a id="c-073"></a>

#### C-073 — Trace, supertrace, and Berezinian
**Local notation:** In the first display, $\hat A,\hat B,\hat P,\hat N$ are operators on graded Fock space. In the Berezinian display, $A,B,C,D,H$ are coefficient supermatrices; $K$ is a scalar Grassmann kernel.
**Source:** §16.6.1, printed 474–477; PDF 484–487; Propositions 153–155, (16.48)–(16.60). For ordinary finite matrices,
$$\operatorname{Str}\hat A=\operatorname{Tr}(\hat P\hat A),\quad \hat P=(-1)^{\hat N},\quad
\operatorname{Str}(\hat A\hat B)=(-1)^{\pi(\hat A)\pi(\hat B)}\operatorname{Str}(\hat B\hat A)$$
for homogeneous matrices. For even invertible supermatrices, the Berezinian obeys
$$
\operatorname{Ber}\begin{pmatrix}A&B\\C&D\end{pmatrix}
=\frac{\det(A-BD^{-1}C)}{\det D},\qquad
\operatorname{Ber}(e^H)=e^{\operatorname{Str}H}.
$$
The source kernel trace has a reflected argument $K(\theta,-\theta^*)$; the supertrace uses $K(\theta,\theta^*)$, with paired Gaussian measure (16.50)–(16.51). Its subsequent rank-one/trace conventions must not be substituted for the ordinary Hilbert-space trace without checking parity and coefficient order; in particular ordinary complex rank-one matrices obey $\operatorname{Tr}|\eta\rangle\langle\phi|=\langle\phi,\eta\rangle$.

<a id="c-074"></a>

#### C-074 — Weyl translation expansion and its supertrace inversion
**Local notation:** The coefficient $A_w$ and the contravariant Weyl symbol are scalar Grassmann functions, while $\hat A,\hat T$ are operators.
**Source:** §16.6.2, Propositions 156–158, printed 477–480; PDF 487–490; (16.61)–(16.73).
$$
\widehat A=\int A_w(\xi)\hat T(-\xi)d^2\xi,\qquad
A_w(\xi)=\operatorname{Str}(\widehat A \hat T(\xi)).
$$
The contravariant Weyl symbol is the Grassmann Fourier transform of this translation coefficient; its values for $\hat I,\hat c,\hat c^\dagger$ are $1,\alpha,\alpha^*$, respectively. Supertrace integration of a contravariant symbol gives $\operatorname{Str}\widehat A$, and integrating the product of two such symbols gives $\operatorname{Str}(\hat A\hat B)$, (16.70)–(16.71). Ordinary trace has different formulas (C-076), so this distinction changes determinants in exponential traces.

<a id="c-075"></a>

#### C-075 — Normal, antinormal, and Weyl characteristic functions
**Local notation:** The symbols $P,Q$ here are scalar distributions, unrelated to the canonical operators $\hat P,\hat Q$; $s$ is a scalar ordering parameter.
**Source:** §16.6.3, Propositions 159–160, Theorem 85, printed 481–485; PDF 491–495; (16.75)–(16.113).
$$
\hat T_N(\gamma)=e^{\hat c^\dagger\gamma}e^{-\gamma^*\hat c}=\hat T(\gamma)e^{+\gamma^*\gamma/2},
\quad \hat T_A(\gamma)=e^{-\gamma^*\hat c}e^{\hat c^\dagger\gamma}=\hat T(\gamma)e^{-\gamma^*\gamma/2},
$$
$$
\hat T_s(\gamma)=\hat T(\gamma)e^{s\gamma^*\gamma/2}.
$$
The Fourier transforms of $\operatorname{Tr}(\hat\rho \hat T_s)$ give the corresponding quasiprobabilities. The source's $P,Q$ distributions reproduce normally and antinormally ordered monomial moments exactly through finite Berezin integrals. Physical density operators in this construction are assumed even. Ordering conversions are finite nilpotent multipliers, not bosonic heat smoothing.

<a id="c-076"></a>

#### C-076 — Ordinary trace of fermionic Weyl-symbol products
**Local notation:** The functions $H^w,G_w,G^w,A^w$ are scalar Grassmann symbols or translation coefficients; $\hat H,\hat G,\hat A$ are the corresponding operators.
**Source:** §16.6.3, Proposition 160, printed 485; PDF 495; (16.114), (16.115).
$$
\operatorname{Tr}(\widehat H\widehat G)=2^{-n}\widehat I_\zeta\!\left[H^w(\zeta)G_w(-2\zeta)\right].
$$
**Parity correction:** the source prints $+2\zeta$, which is valid when the relevant symbols are even. For arbitrary ordinary complex-coefficient operators, including odd ones, the required argument is $-2\zeta$, with paired integration as in C-067. In one mode take $\hat H=\hat c$, $\hat G=\hat c^\dagger$: $H^w=\alpha$, $G_w=-\alpha^*$. The printed positive argument gives $-1$, but $\operatorname{Tr}(\hat c\hat c^\dagger)=1$; reflecting the argument gives $+1$. Independent full-matrix tests with mixed even/odd coefficients verify the corrected formula in one, two and three modes. This is a restriction/correction of the source's unqualified mixed trace formula, not the earlier subscript extraction issue.

This is a **mixed-symbol** formula: $H^w$ is the contravariant Weyl symbol, while $G_w$ is the covariant translation coefficient of C-074. The source's subscript/superscript distinction is visually confirmed; replacing $G_w$ by $G^w$ would already fail for two identities. Carry over both the scaled argument and its Berezin Jacobian together. For an ordinary complex-coefficient operator, an equivalent consequence of the finite Clifford/Weyl basis is $\operatorname{Tr}\hat A=2^n A^w(0)$, where 0 means the scalar coefficient in the physical Grassmann labels. Hence $\operatorname{Tr}(\hat H\hat G)=2^n(H^w\star G^w)(0)$, not an unweighted Berezin integral of two contravariant symbols.

<a id="c-077"></a>

#### C-077 — Fermionic star product and star exponential
**Local notation:** The quantities $G,H,U^w$ in the star equation are scalar Grassmann symbols; the quantum propagator is $e^{-it\hat H}$.
**Source:** §16.6.4, Lemmas 106–107, Proposition 161, printed 485–488; PDF 495–498; (16.117)–(16.123).
$$
G\star H=G\exp\!\left[\tfrac12\left(
\overleftarrow\partial_{\alpha^*}\cdot\overrightarrow\partial_\alpha+
\overleftarrow\partial_\alpha\cdot\overrightarrow\partial_{\alpha^*}\right)\right]H.
$$
All odd derivatives carry their arrow; the plus sign is the fermionic counterpart of the bosonic symplectic difference. The derivative exponential terminates because there are finitely many odd variables. **Consequence:** the symbol of $e^{-it\widehat H}$ is the solution of $i\partial_tU^w=H\star U^w$, $U^w_0=1$; it is generally not $e^{-itH}$ in the ordinary Grassmann product. Quadratic $H$ reduces this exact finite system to Riccati and scalar amplitude equations (C-083).

<a id="c-078"></a>

#### C-078 — Fermi oscillator exponential and removable Cayley singularity
**Source:** §16.7.1, printed 488; PDF 498; (16.125), (16.126).
$$
e^{z \hat c^\dagger \hat c}=\hat I+(e^z-1)\hat c^\dagger \hat c,\qquad
(\hat c^\dagger \hat c)^w=\alpha^*\alpha+\tfrac12,
$$
$$
[e^{z \hat c^\dagger \hat c}]^w=
\frac{e^z+1}{2}\exp[2\tanh(z/2)\alpha^*\alpha]
=\frac{e^z+1}{2}+(e^z-1)\alpha^*\alpha.
$$
The polynomial expression is entire for all complex $z$. Zeros of $e^z+1$ make the exponential chart singular but not the operator or symbol. This is a simple check on determinant branches and fermionic Gaussian prefactors.

<a id="c-079"></a>

#### C-079 — Fermionic partition determinant
**Local notation:** Here $h$ is explicitly the finite matrix of scalar one-particle coefficients, and $Z(\beta,\mu)$ a scalar partition function; $\hat c$ acts on Fock space.
**Source:** §16.7.2, Lemma 108, printed 488–490; PDF 498–500; (16.127)–(16.133).
$$
\operatorname{Tr}e^{\sum_jz_j\hat c_j^\dagger \hat c_j}=\prod_j(1+e^{z_j}),
\quad Z(\beta,\mu)=\prod_j[1+e^{-\beta(\omega_j-\mu)}].
$$
The occupation mean is $\sum_j[1+e^{\beta(\omega_j-\mu)}]^{-1}$. **Consequence for finite matrices:** by diagonalization for Hermitian one-particle $h$, and analytic continuation,
$$
\operatorname{Tr}e^{\hat c^\dagger h \hat c}=\det(I_n+e^h),\qquad
\operatorname{Str}e^{\hat c^\dagger h \hat c}=\det(I_n-e^h).
$$
The second determinant follows by weighting every occupied mode by $-1$. These are ordinary finite-mode traces; no square-root sign ambiguity occurs in the number-preserving form.

<a id="c-080"></a>

#### C-080 — General quadratic fermion Hamiltonian and BdG matrix
**Local notation:** The BdG matrices $M,L,\Sigma,F_t,A_t,B_t$ have scalar entries; $\hat H,\hat c,\hat c^\dagger$ are Fock-space operators.
**Source:** §16.7.3, printed 490–491; PDF 500–501; (16.134)–(16.140). The book uses
$$
\widehat H=\tfrac12(\hat c^\dagger M\hat c+\hat c^\dagger L\hat c^\dagger-\hat c\bar L\hat c)-\tfrac14(\operatorname{Tr}M)\hat I,
\quad M=M^\dagger,\quad L=-L^T.
$$
The associated $2n$-matrix is
$$
\Sigma=\begin{pmatrix}M/2&L\\-\bar L&-\bar M/2\end{pmatrix},\qquad
F_t=e^{it\Sigma}=\begin{pmatrix}A_t&B_t\\\bar B_t&\bar A_t\end{pmatrix}.
$$
The source's chosen conjugation is $e^{-it\hat H}(\hat c,\hat c^\dagger)^Te^{it\hat H}=F_t(\hat c,\hat c^\dagger)^T$. With the opposite Heisenberg conjugation use $F_{-t}$. For time-dependent coefficients replace the simple matrix exponential by the fundamental solution. $F_t$ is unitary and preserves the symmetric CAR form; this replaces the bosonic symplectic condition.

<a id="c-081"></a>

#### C-081 — Spin implementation of fermionic canonical transformations
**Source:** §16.7.3, printed 491–492; PDF 501–502; (16.141). The connected, parity-preserving quadratic canonical transformations considered here have real Majorana matrices in $\mathrm{SO}(2n)$. They are implemented by even unitaries, unique up to phase; a continuous quadratic lift has the Spin double-cover sign ambiguity. The two parity sectors are invariant under these even implementers. General CAR-preserving real Majorana transformations can instead lie in $\mathrm{O}(2n)$ with determinant $-1$; for example $\hat c\leftrightarrow \hat c^\dagger$ in one mode is implemented by an odd operator and is outside this connected/even statement. A Grassmann translation transforms covariantly under the implementer; phase conventions must match the chosen $\hat U \hat c \hat U^{-1}$ or $\hat U^{-1}\hat c\hat U$ orientation.

<a id="c-082"></a>

#### C-082 — Fermionic squeezed vacuum and a true chart boundary
**Local notation:** The pairing $\Gamma_t$ and the canonical blocks $A_t,B_t,F_t$ are finite matrices; $\Omega,\psi_\alpha$ are state vectors.
**Source:** §16.7.3, printed 492–494; PDF 502–504; (16.142)–(16.147), Proposition 162. While $\det A_t\ne0$,
$$
\hat U_t\Omega=\det(A_t)^{1/2}\exp\!\left[\tfrac12\hat c^\dagger\Gamma_t \hat c^\dagger\right]\Omega,
\qquad\Gamma_t=B_t^T(A_t^T)^{-1},\quad\Gamma_t^T=-\Gamma_t.
$$
The root is continued from $t=0$. The source then says $\hat U_t\psi_\alpha=\hat T(F_t\alpha)\hat U_t\Omega$, but this reverses its conjugation orientation. **Corrected consequence of C-080/C-071:** $\hat U_t\hat T(\alpha)\hat U_t^{-1}=\hat T(F_t^{-1}\alpha)$, hence $\hat U_t\psi_\alpha=\hat T(F_t^{-1}\alpha)\hat U_t\Omega$. Here the doubled label $(\alpha,\alpha^*)$ is transformed together. For one number mode, the label must evolve as $e^{-itm/2}\alpha$, while its $F_t$ has upper block $e^{+itm/2}$. Unlike bosonic positive-Siegel propagation, $A_t$ can become singular at finite real time. The Fock vector remains smooth; its vacuum-overlap chart ceases to describe it. The source's two-mode example at $t=\pi/2$ becomes a fully occupied pair rather than an invertible Grassmann Gaussian.

<a id="c-083"></a>

#### C-083 — Fermionic Riccati and determinant amplitude
**Local notation:** Here $\Gamma,M,L,A_t$ are finite coefficient matrices and $s$ is the scalar vacuum amplitude.
**Source:** §16.7.3, printed 492–493; PDF 502–503; (16.144)–(16.146).
$$
i\dot\Gamma=\tfrac12(M\Gamma+\Gamma\bar M)+L+\Gamma\bar L\Gamma,\quad\Gamma_0=0,
$$
$$
i\dot s=-\tfrac12\operatorname{Tr}(\bar L\Gamma+M/2)s,\quad s_0=1,
\quad s=\det(A_t)^{1/2}.
$$
The positive determinant exponent is the fermionic counterpart of the negative exponent in bosonic Gaussian integration. It is fixed by the number-conserving test $L=0$: $A_t=e^{itM/2}$, so $s=e^{it\operatorname{Tr}M/4}$, exactly the phase from the scalar term in C-080.

<a id="c-084"></a>

#### C-084 — Fermionic Mehler kernel and coherent matrix elements

**Local notation:** $K_{\hat U}$ is a scalar Grassmann kernel; $A,B$ are finite numerical canonical blocks and $a_0$ is a scalar amplitude. **Source:** §16.7.4, Proposition163, Corollary56, printed495–496 / PDF505–506, (16.148)–(16.150); compare Berezin (5.15)–(5.16), B-060–B-061. **Explicit convention-fixed reconstruction:** for
$$
\hat U\hat{\mathbf c}\hat U^{-1}=A\hat{\mathbf c}+B\hat{\mathbf c}^{\dagger},\qquad \det A\ne0,
$$
the holomorphic kernel is
$$
K_{\hat U}(a^*,a)=a_0\exp\!\left[-\tfrac12\left(
a^T\bar B A^{-1}a-2a^{*T}A^{-1}a+a^{*T}A^{-1}Ba^*\right)\right],
\qquad a_0=\langle0|\hat U|0\rangle.
$$
All $a,a^*$ here are Grassmann scalar labels; the written order fixes the signs. In the homogeneous quadratic dynamical lift C-080, $a_0=\det A^{1/2}$ continuously from identity. General unitary implementers have $|a_0|=\det(AA^\dagger)^{1/4}$ and a separately fixed phase. To reconstruct normalized coherent matrix elements, use the corrected inverse label $F_t^{-1}\alpha$ in C-082, not the source's forward label. Every occupation-basis coefficient has been checked against independent CAR matrix exponentials in complex one-, two- and three-mode examples. At $\det A=0$, the full finite polynomial remains regular; compute it in the occupation basis or use another vacuum chart. The explicit operator factorization of this kernel is X-009.



Source page: Combescure–Robert, printed p.495, one-based PDF p.505.


Combescure–Robert, printed p.496, one-based PDF p.506.


<a id="c-085"></a>

#### C-085 — Fermionic Cayley symbol of a quadratic exponential
**Local notation:** The symbol $U_t^w$ is scalar, while $\hat U_t=e^{-it\hat H}$ is the operator; $F_t,\Sigma,\mathcal N_t$ are finite coefficient matrices.
**Source:** §16.7.4, Theorem 86, printed 496–498; PDF 506–508; (16.151)–(16.160). **As printed, with a sign correction required below:** in the chart $\det(F_t+I_{2n})\ne0$,
$$
U_t^w(\xi)=\det\!\left(\frac{F_t+I_{2n}}{2}\right)^{1/2}
\exp\!\left[
\binom{\xi^*}{\xi}\cdot(F_t-I_{2n})(F_t+I_{2n})^{-1}\binom\xi{\xi^*}
\right].
$$
**The displayed formula is the source as printed and has an orientation/sign inconsistency.** With C-080's $\hat U_t=e^{-it\hat H}$, $F_t=e^{it\Sigma}$, its derivative at zero is $+iH^w$, whereas the star equation requires $-iH^w$. A one-mode check sets $L=0,M=m$: C-078 gives
$$
[\hat U_t]^w=\cos(tm/4)\exp[-2i\tan(tm/4)\xi^*\xi],
$$
while the printed Cayley exponent has the plus sign. **Corrected re-expression in those conventions:** keep the continued determinant root and put a minus before the displayed Cayley exponent, equivalently use $F_t^{-1}$ in its Cayley transform. The source's own line immediately before (16.158) finds $\mathcal N_t=F_t^{-1}$, but its subsequent (16.158) substitutes $F_t$ without the resulting minus. The determinant root starts at 1; $\det F_t=1$, so replacing $F_t$ by its inverse leaves that continued prefactor unchanged. The derivation uses the fermionic star equation, a matrix Riccati equation and Liouville's formula. The full finite Grassmann symbol is entire even when this Gaussian chart fails; do not evaluate singular factors separately at $\det(F+I_{2n})=0$.

<a id="c-086"></a>

#### C-086 — Covariant fermionic Cayley symbol and Pfaffian branch

**Local notation:** $a_t$ and $U_{t,w}$ are scalar Grassmann data; $J_F,F_t,C_t$ are finite numerical matrices. **Source:** §16.7.4, Corollary57, printed498 / PDF508, (16.161). **Independently corrected inverse-Cayley exponent and amplitude square:** fix C-080's direction
$$
\hat U_t=e^{-it\hat H},\qquad
\hat U_t\binom{\hat c}{\hat c^\dagger}\hat U_t^{-1}=F_t\binom{\hat c}{\hat c^\dagger},\quad F_t=e^{it\Sigma},\quad
J_F=\begin{pmatrix}0&I_n\\I_n&0\end{pmatrix}.
$$
For $\det(F_t-I_{2n})\ne0$, the translation coefficient defined in C-074 is
$$
U_{t,w}(\gamma)=a_t\exp\!\left[-\tfrac14\binom\gamma{\gamma^*}^{\!T}
J_F(F_t+I_{2n})(F_t-I_{2n})^{-1}\binom\gamma{\gamma^*}\right],
\quad a_t=\operatorname{Str}\hat U_t,\quad a_t^2=(-1)^n\det(F_t-I_{2n}).
$$
The source prints a plus in the exponent and omits $(-1)^n$ in its amplitude square. Both fail a one-mode calculation; neither discrepancy is repaired by choosing another square root. On the overlap with C-085, the sign is fixed explicitly by
$$
C_t=(F_t-I_{2n})(F_t+I_{2n})^{-1},\quad
p_t=\det((F_t+I_{2n})/2)^{1/2},\quad p_0=1,
\qquad
a_t=(-1)^{n(n-1)/2}p_t\operatorname{Pf}(-2J_FC_t).
$$
The Pfaffian uses the block coordinate order of C-069 and the root follows the same quadratic lift. Continue this scalar when changing charts. In one mode $\hat H=(m/2)(\hat N-\hat I/2)$, $\vartheta=tm/4$ gives directly
$$
U_{t,w}(\gamma)=2i\sin\vartheta-\cos\vartheta\,\gamma^*\gamma
=2i\sin\vartheta\exp[(i/2)\cot\vartheta\,\gamma^*\gamma].
$$
Hence $a_t^2=-4\sin^2\vartheta$ while $\det(F_t-I_2)=4\sin^2\vartheta$. All Grassmann coefficients also match exact complex noncommuting one-, two- and three-mode CAR calculations. This chart is singular at identity, but the full translation coefficient is a finite polynomial and remains regular through an apparent Gaussian chart pole.

<a id="c-087"></a>

#### C-087 — Fermion pair rotations are compact disentanglings
**Consequence connecting the source CAR matrix model (§16.2, printed 458–459; PDF 468–469) with its SU(2) Gauss lemma (§7.4.1, printed 211–213; PDF 226–228).** For two modes,
$$
\hat J_+=\hat c_1^\dagger \hat c_2^\dagger,\quad \hat J_-=\hat c_2\hat c_1,\quad \hat J_3=(\hat N_1+\hat N_2-\hat I)/2
$$
obey SU(2), not the bosonic SU(1,1), relations. Thus C-041 applies to $e^{\xi \hat c_1^\dagger \hat c_2^\dagger-\bar\xi \hat c_2\hat c_1}$:
$$
\hat U\Omega=\cos r\,\Omega+\frac\xi r\sin r\,\hat c_1^\dagger \hat c_2^\dagger\Omega.
$$
The Gaussian vacuum chart fails at $\cos r=0$, while this finite vector is smooth. This identifies exactly why bosonic squeeze formulas have hyperbolic functions and fermionic pair formulas have trigonometric ones.

### Superalgebra and appendix tools

<a id="c-088"></a>

#### C-088 — Mixed boson–fermion displacement
**Local notation:** Here $z$ is a commuting scalar and $\gamma$ an anticommuting scalar; the hatted bosonic and fermionic generators act on their tensor-product module.
**Source:** §17.6, printed 523–525; PDF 533–535; (17.89)–(17.104).
$$
\hat T(z,\gamma)=\exp[z\hat a_B^\dagger-\bar z \hat a_B+\hat c_F^\dagger\gamma-\gamma^*\hat c_F].
$$
Bosonic and the even fermionic displacement generators commute, so this factors as their product. In consistent ordering,
$$
\hat T(z,\gamma)=e^{-(|z|^2+\gamma^*\gamma)/2}
 e^{z\hat a_B^\dagger}e^{\hat c_F^\dagger\gamma}e^{-\bar z \hat a_B}e^{-\gamma^*\hat c_F}.
$$
Its coherent vector has overlap equal to the product of the bosonic overlap C-008 and fermionic overlap C-072, and is a simultaneous annihilator eigenvector. This factorization states all odd coefficient positions explicitly.

<a id="c-089"></a>

#### C-089 — Super star product
**Local notation:** The factors $G,H$ are scalar super-symbols; the derivative operations retain their standard notation.
**Source:** §17.7, printed 525–526; PDF 535–536; (17.105)–(17.112). The super star product adds the bosonic antisymmetric derivative contraction and fermionic symmetric contraction:
$$
G\star H=G\exp\tfrac12\left[
\overleftarrow\partial_\zeta\overrightarrow\partial_{\bar\zeta}
-\overleftarrow\partial_{\bar\zeta}\overrightarrow\partial_\zeta
+\overleftarrow\partial_{\alpha^*}\overrightarrow\partial_\alpha
+\overleftarrow\partial_\alpha\overrightarrow\partial_{\alpha^*}
\right]H.
$$
Use the source's complex bosonic coordinates and its chosen $\hbar=1$ normalization. The odd derivative part terminates; the even part has the usual analytic/asymptotic symbol restrictions. Chapter 17's Dicke-model example closes only after retaining a finite matrix of bosonic operators; it does not supply a general exact scalar Gaussian factorization for the interacting model.

<a id="c-090"></a>

#### C-090 — Supercommutator and nilpotent continuation
**Local notation:** The homogeneous $\hat A,\hat B$ are operators; their parities $\pi(\hat A),\pi(\hat B)$ and the auxiliary odd parameters are scalars.
**Source:** §17.4.1–17.4.3, printed 506–511; PDF 516–521; (17.15)–(17.34). For homogeneous operators,
$$
[\hat A,\hat B]_s=\hat A\hat B-(-1)^{\pi(\hat A)\pi(\hat B)}\hat B\hat A.
$$
The source gives graded antisymmetry/Jacobi and the superderivation rule. Odd-odd brackets are anticommutators; BCH for a group element is applied to **even total generators**, after odd parameters are included. Smooth functions of an even scalar body plus even nilpotent correction use the finite Taylor formula (17.20). This is the precise setting in which boson/fermion exponential formulas can be combined.

<a id="c-091"></a>

#### C-091 — Commuting Gaussian integral and branch
**Local notation:** Here $M,A$ are finite numerical quadratic-form matrices and $j,x$ scalar-coordinate columns.
**Source:** Appendix A.1, Theorem 88, printed 531; PDF 540; (A.1). **Re-expression in a standard unambiguous convention:** if $M=M^T$ and $\Re M>0$,
$$
\int_{\mathbb R^m}e^{-x^TMx/2+j^Tx}\,dx
=(2\pi)^{m/2}\det(M)^{-1/2}e^{j^TM^{-1}j/2}.
$$
Choose the analytic square root positive on positive real $M$. Oscillatory limits with nonnegative imaginary part of $A$ follow by $M=-iA$ and regularization. For real nonsingular symmetric $A$, the Fresnel factor is $e^{i\pi\operatorname{sgn}(A)/4}|\det A|^{-1/2}$. This fixes a missing factor $1/4$ in the phase printed earlier in Lemma 16/(2.65), and the integration-variable typo in (A.1).

<a id="c-092"></a>

#### C-092 — Stationary phase is asymptotic, with a determinant on normal directions
**Source:** Appendix A.2 and A.4, printed 531–536; PDF 540–545; Theorems 29, 89, Proposition 164. If the phase has nonnegative imaginary part and its critical manifold is clean, stationary phase integrates transverse Gaussian directions and gives the inverse square root of the transverse Hessian. Analytic saddle points extend this to complex phases. These theorems justify semiclassical propagator/trace approximations; they do not make a nonquadratic exponential equal to its leading Gaussian. Appendix A.3 checks the hydrogen coherent-state Hessian, while A.5 gives the sphere/disc Kähler metric and Poisson brackets.

<a id="c-093"></a>

#### C-093 — Controlled Moyal expansion and functional calculus
**Source:** Appendix B, Theorem 90, printed 539–541; PDF 547–549; (B.1)–(B.13); Appendix C, Proposition 165, Theorem 91, printed 543–546; PDF 550–553; (C.1)–(C.7). Symbol products have a derivative-controlled $O(\hbar^{N+1})$ remainder under the stated polynomial-growth symbol hypotheses. Functional calculus is built from resolvent parametrices and Cauchy/almost-analytic integration. Hence applying a smooth function to an operator differs from applying it pointwise to its symbol by star-product corrections. The source's growth assumptions for the functional-calculus theorem must be checked before applying it to an oscillatory exponential; quadratic formulas C-030/C-053 bypass that issue directly.

<a id="c-094"></a>

#### C-094 — Adjoint exponential and local Lie coordinates
**Local notation:** Here $X,Y$ are abstract Lie-algebra elements and $\operatorname{Ad},\operatorname{ad}$ retain their standard Lie-theoretic notation. In a Hilbert-space representation their operator images receive hats.
**Source:** Appendix D.2.2, Theorem 92, printed 548–549; PDF 555–556.
$$
\operatorname{Ad}(e^X)=e^{\operatorname{ad}X},\qquad
\operatorname{ad}X(Y)=[X,Y].
$$
The exponential map is a local diffeomorphism near the identity. Thus exponential conjugation can be computed by exponentiating the matrix of $\operatorname{ad}X$ on any finite closed commutator span. **Caution:** the following source sentence claims a global diffeomorphism for every connected simply connected group; that is false (SU(2) is a counterexample). Global factorizations require their own hypotheses/charts.

<a id="c-095"></a>

#### C-095 — Schur uniqueness, square integrability, and group coherent states
**Source:** Appendix D.3–D.4, Lemma 111, Theorems 94–95, Proposition 169, printed 551–555; PDF 558–562; (D.3)–(D.11). Irreducibility makes a bounded commuting implementer scalar, establishing uniqueness up to phase. A fiducial vector and its projective stabilizer produce $G/H$-labelled coherent states with stabilizer cocycles and overlap matrix coefficients. A quotient invariant measure and the finite square-integral condition are needed for resolution of identity. Square-integrable nonunimodular representations have an admissibility operator; in the unimodular case it is scalar. This is why every formal group orbit does not automatically yield a normalized coherent-state integral calculus.

<a id="c-096"></a>

#### C-096 — Cartan polar decomposition
**Local notation:** Here $k$ and $e^X$ are numerical group matrices, with $X$ a matrix Lie-algebra element; represented quantum rotation and squeeze operators receive hats.
**Source:** Appendix D.3.3, Proposition 168, printed 554; PDF 561. For a connected linear reductive group, $(k,X)\mapsto k e^X$ is a diffeomorphism from the maximal compact subgroup times the noncompact Lie-algebra component onto the group. This is a structural explanation for rotation–squeeze decompositions. It does not imply that one arbitrary Lie-algebra exponential chart is global.

<a id="c-097"></a>

#### C-097 — Berezin symbol product on the disk
**Local notation:** The functions $A_c,B_c$ are scalar covariant symbols; $(\hat A\hat B)_c$ means the scalar symbol of the operator product.
**Source:** Appendix E, printed 557–560; PDF 563–566; (E.1), (E.2), Proposition 171. The normalized off-diagonal covariant symbol is $A_c(z,\bar w)=\langle\psi_z,\hat A\psi_w\rangle/\langle\psi_z,\psi_w\rangle$. Its exact product is
$$
(\hat A\hat B)_c(z,\bar z)=\int_D A_c(z,\bar w)B_c(w,\bar z)|\langle\psi_z,\psi_w\rangle|^2d\nu_n(w),
$$
$$
|\langle\psi_z,\psi_w\rangle|^2=
\left[\frac{(1-|z|^2)(1-|w|^2)}{|1-\bar zw|^2}\right]^n.
$$
Group conjugation transports the covariant symbol by the inverse disk action. The large-$n$ expansion is a semiclassical calculus. Proposition 170's displayed subleading scalar coefficient fails the normalization test on the constant function; retain the exact integral rather than copying (E.3) as a reliable exponential rule.

<a id="c-098"></a>

#### C-098 — Grassmann Fourier involution with an explicit measure convention

**Source:** Combescure–Robert, §16.4.3, printed p.470 / PDF480, Definition 44 and its following unnumbered properties. **Source definition, convention-fixed re-expression and derived correction.**

Let $\xi_j,\xi_j^*,\alpha_j,\alpha_j^*$ be independent odd Grassmann generators, all mutually anticommuting. A star labels the formal conjugate partner; it does not turn these symbols into commuting complex numbers. Define the paired Berezin functional $\widehat I_\xi$ by

$$
\widehat I_\xi\!\left[\prod_{j=1}^{n}\xi_j\xi_j^*\right]=1,
$$

with all other monomials missing an integrated generator mapped to zero. External generators are kept in a fixed order when extracting this coefficient. This fixes the meaning of the source's shorthand $d\xi^*d\xi$; reversing differential strings without the corresponding sign is forbidden.

Define the linear map $\widehat{\mathcal F}$ on the finite Grassmann algebra by

$$
(\widehat{\mathcal F}f)(\alpha,\alpha^*)
=\widehat I_\xi\!\left[
 e^{\xi^*\cdot\alpha-\alpha^*\cdot\xi}f(\xi,\xi^*)
\right].
$$

Then

$$
f(\xi,\xi^*)=
\widehat I_\alpha\!\left[
 e^{\alpha^*\cdot\xi-\xi^*\cdot\alpha}
 (\widehat{\mathcal F}f)(\alpha,\alpha^*)
\right],
\qquad \widehat{\mathcal F}^{2}=\hat I_{\mathcal G}.
$$

Here $\hat I_{\mathcal G}$ is the identity map on the Grassmann algebra. This is an **involution**. The source calls the displayed square-identity property “idempotent”; the equation, rather than that word, is the intended property.

For one mode the complete basis table is

$$
\widehat{\mathcal F}1=\alpha\alpha^*=-\alpha^*\alpha,
\quad
\widehat{\mathcal F}\xi=\alpha,
\quad
\widehat{\mathcal F}\xi^*=-\alpha^*,
\quad
\widehat{\mathcal F}(\xi\xi^*)=1.
$$

Consequently, for arbitrary mode number,

$$
\widehat{\mathcal F}1=
\prod_{j=1}^{n}\alpha_j\alpha_j^*.
$$

**Source correction.** The source prints $1^{\mathcal F}=\alpha^*\cdot\alpha$. That expression has the wrong one-mode sign in this paired convention and the wrong degree when $n>1$. The top-degree product above is the convention-fixed formula. This finite transform is the missing explicit operation behind C-074's statement that the Weyl symbol is a Grassmann Fourier transform.

**Independent check:** Exact exterior-algebra arithmetic verifies the transform, its square and the displayed basis table on every monomial for $n=1,2,3$. Checks `Grassmann_Fourier_involution_all_monomials_n1/n2/n3` and `Grassmann_delta_reproduction_all_monomials_n1/n2/n3`; this addition does not claim to verify all six of the source's additional Parseval/convolution/modulation statements.


The one-pair basis table proves the involution in one mode. The exponential kernel factors into even one-pair kernels, and paired integration factors in the same order; this extends the identity to arbitrary finite mode number. The listed exact finite checks test that implementation and its sign conventions. The other Parseval/convolution/modulation assertions on the source page remain source lookups.



Source page: Combescure–Robert, printed p.470, one-based PDF p.480.


<a id="c-099"></a>

#### C-099 — Holomorphic Grassmann reproducing kernel

**Source:** Combescure–Robert, §16.4.1, printed p.468 / PDF478, (16.23)–(16.24). **Source formula in the same convention-fixed notation.**

For any holomorphic polynomial $\psi(\theta_1,\ldots,\theta_n)$ in odd generators, use a fresh independent pair $\gamma,\gamma^*$ and the paired functional $\widehat I_\gamma$ defined above. Then

$$
\psi(\theta)=
\widehat I_\gamma\!\left[
\psi(\gamma)E(\theta,\gamma,\gamma^*)
\right],
$$

$$
E(\theta,\gamma,\gamma^*)=
 e^{-(\theta-\gamma)\cdot\gamma^*}
=\prod_{j=1}^{n}
\left[1-(\theta_j-\gamma_j)\gamma_j^*\right].
$$

The bracket order is explicit. This is a scalar Grassmann kernel representing the identity map $\hat I_{\mathrm{hol}}$ on holomorphic states, not a new Fock-space operator. It explains how a finite fermionic kernel calculus can reproduce vectors without a commuting-variable Dirac distribution.

**Interpretation of the following source line.** The book also writes a polynomial delta representative involving $\prod_j(\theta_j-\gamma_j)(\theta_j^*-\gamma_j^*)$. It is a different Grassmann polynomial from $E$. Equality of their reproducing action on the intended holomorphic tests must be distinguished from equality of the kernels as polynomials. The displayed exponential kernel is independently checked on every holomorphic monomial for $n=1,2,3$ in the paired convention.


### Complete reading coverage and page map

Every supplied PDF page, 1–582, was read in order for potentially relevant results, including prose, displayed formulas, appendices, bibliography, and indexes. The dossier selects operator-exponential consequences from that full reading; omission of a topic from the entries does not indicate an unread chapter. Initial tool-output truncations were recovered by smaller reads, including the slivers on PDF 27 and 71. No supplied page remains unread. This is full **text reading with selective visual mathematics checks**, not a claim to have visually inspected all 582 raster pages.

The text has 583 form-feed segments because the final segment is empty; the PDF has exactly 582 real pages.

| Supplied PDF pages (one-based) | Printed pages | Complete reading coverage / relevance |
|---|---|---|
| 1–18 | Front matter, including v–xvii | Cover, series/title/copyright pages, both complete prefaces and contents. |
| 19 | Part I divider, inferred 1 | Part divider read; no displayed Arabic page number. |
| 20–43 | 3–26 | Chapter 1: CCR, Weyl operators, canonical coherent states, Bargmann space, Mehler formula. C-001–C-012. |
| 44–79 | 27–62 | Chapter 2: Weyl/Wick/anti-Wick symbols, Wigner calculus, Moyal product, covariance, symbol/operator estimates. C-013–C-017, C-091/C-093. |
| 80–107 | 63–90 | Chapter 3: quadratic flows, Gaussian and coherent propagation, metaplectic lift, squeezes, Cayley symbols, periodic quadratic dynamics. C-018–C-031, C-035. |
| 108–143 | 91–126 | Chapter 4: propagation of coherent states, polynomial corrections, finite/Ehrenfest-time and analytic estimates. C-029, C-032–C-034. |
| 144–173 | 127–156 | Chapter 5: trace localization, stationary phase, periodic-orbit/Bohr–Sommerfeld results. C-034, C-036. |
| 174 | Part II divider, inferred 157 | Part divider read. |
| 175–207 | 159–191 | Chapter 6: torus states, finite Weyl matrices, cat maps and exact Egorov/boundary restrictions; all ergodicity subsections read. C-037. |
| 208–249 | 193–234 | Chapter 7: SU(2), representations, spin coherent states, symbol inequalities and contraction. C-038–C-044. |
| 250–287 | 235–272 | Chapter 8: SU(1,1), disk representations, oscillator realizations, squeezes, time-dependent dynamics, singular potential and wavelets. C-045–C-052. |
| 288–308 | 273–293 | Chapter 9: hydrogen coherent states, SO(4), constrained/S3 geometry and semiclassical dynamics. Read throughout; no further general CCR/CAR exponential identity selected. Its Gaussian Hessian computation is also covered in Appendix A.3. |
| 309 | Part III divider, inferred 295 | Part divider read. |
| 310–325 | 297–312 | Chapter 10: extended Weyl calculus, oscillator symbol, commutator criteria, matrix symbols, functional calculus. C-053, C-093. |
| 326–360 | 313–347 | Chapter 11: complex-phase Fourier integral operators, Herman–Kluk representation, transport and nonlinear remainder estimates. C-054. |
| 361–383 | 349–371 | Chapter 12: semiclassical measures and quantum ergodicity, with complete proofs/applications read. No additional exact creation/annihilation disentangling selected. |
| 384–415 | 373–404 | Chapter 13: partial trace, Gaussian open systems, quadratic master equation, purity/decoherence and cat-state applications. C-055–C-056. |
| 416–439 | 405–428 | Chapter 14: matrix/multicomponent adiabatic decoupling, almost projectors, Dirac example and propagated wave packets. Read throughout; operator-valued symbol asymptotics do not add an exact scalar CCR/CAR factorization. |
| 440 | Part IV divider, inferred 429 | Part divider read. |
| 441–466 | 431–456 | Chapter 15: bosonic Fock fields, analytic vectors, Weyl states, mean-field fluctuations, Duhamel/Gevrey/Borel analysis. C-057–C-063. |
| 467–508 | 457–498 | Chapter 16: CAR, Grassmann calculus/integrals, coherent states, trace/supertrace, symbols/star product and quadratic fermionic evolution. C-064–C-087. |
| 509–539 | 499–529 | Chapter 17: supergeometry, mixed boson–fermion quantization, super-Poisson/Moyal calculus, supersymmetric oscillator and anharmonic examples. C-088–C-090. |
| 540–546 | 531–537 | Appendix A, all A.1–A.5: Gaussian integration, clean/complex stationary phase, hydrogen Hessian, Kähler geometry. C-091–C-092. |
| 547–549 | 539–541 | Appendix B: entire Moyal remainder proof. C-093. |
| 550–553 | 543–546 | Appendix C: entire functional-calculus construction. C-093. |
| 554–562 | 547–555 | Appendix D: Lie groups/algebras, exponential/adjoint, Haar measure, representations, Cartan decomposition and group coherent states. C-094–C-096. |
| 563–566 | 557–560 | Appendix E: exact and asymptotic Berezin symbol calculus on the disk. C-097. |
| 567–576 | 561–570 | All bibliography pages read; references identify background, not additional formulas claimed as read from their original works. |
| 577–580 | 571–574 | Entire concepts index read and checked for missed operator-exponential topics. |
| 581–582 | 575–576 | Entire symbols index read. |

#### Piecewise conversion rule

No single `printed + offset` rule is valid. For numbered body pages, the following intervals have been checked against actual page headers/chapter openings. In each row, `PDF = printed + offset`.

| Printed interval | PDF interval | Offset |
|---|---|---|
| 3–156 | 20–173 | +17 |
| 159–191 | 175–207 | +16 |
| 193–293 | 208–308 | +15 |
| 297–347 | 310–360 | +13 |
| 349–371 | 361–383 | +12 |
| 373–428 | 384–439 | +11 |
| 431–529 | 441–539 | +10 |
| 531–537 | 540–546 | +9 |
| 539–541 | 547–549 | +8 |
| 543–555 | 550–562 | +7 |
| 557–576 | 563–582 | +6 |

The entry-specific printed/PDF pairs take precedence over any approximate global offset. Part dividers PDF 19, 174, 309, 440 have no visible Arabic folio; their positions correspond to printed 1, 157, 295, 429 respectively.

#### Missing-page report

There are **no missing supplied PDF pages** and no detected gap in substantive chapter, appendix, bibliography, or index text. Relative to the printed Arabic sequence 1–576, the supplied electronic file omits the following leaf numbers: **2, 158, 192, 294, 296, 348, 372, 430, 530, 538, 542, 556**. These are inferred blank separator/verso leaves: every gap lies at a part, chapter, or appendix boundary, and adjoining substantive text is complete. They were not supplied and therefore were not read; the claim is complete reading of the supplied PDF, not visual inspection of those absent printed leaves. The arithmetic is consistent: 576 Arabic-number positions minus 12 omitted blanks gives 564 supplied body/divider pages; adding 18 front-matter pages gives 582.

### Verification and source-ambiguity ledger

Rendered page inspection covered PDF **23, 30, 32, 38, 39, 41, 48, 61, 81, 83, 89, 94, 96, 98, 100, 102, 185, 228, 271, 279, 280, 475, 476, 481, 482, 486, 495, 497, 498, 500, 501, 502, 503, 506, 507, 508, 540, 556, 564**. The selected pages cover BCH signs, displacement/coherent factors, symbol normalizations, quadratic/squeeze matrices, rank-one disentangling, Grassmann differential order, ordinary versus supertrace, fermionic Riccati/Cayley signs, Gaussian branches and appendix global-group claims.

Checks use the identity/zero-generator limit, vacuum normalization, first derivatives at zero, known number spectra, one-mode CAR matrices, and Gaussian integrals. These are independent algebraic checks and are labeled re-expressions/consequences where they change the printed source. They do not amount to a separate proof of every long semiclassical theorem.

| Source location | Issue and treatment |
|---|---|
| §1.3, (1.41), printed 15 / PDF 32 | Number-expansion indices/factorial formatting are inconsistent; C-007 follows (1.36) and the ladder construction. |
| §1.4, printed 20–25 / PDF 37–42 | Bargmann variable and ħ scales vary among displayed formulas. C-011 fixes ħ=1 and states its multiplication/differentiation convention explicitly. |
| (2.15), printed 31 / PDF 48 | Heat Weyl prefactor has power −n/2; C-053 and the Mehler integral give −n in these coordinates. The ħ/time scaling also requires the declared convention. |
| (2.16), printed 31 / PDF 48 | The inverse in displacement conjugation conflicts with the stated shift sign. C-004 fixes it directly on $\hat Q$ and $\hat P$. |
| Lemma 16/(2.65), printed 44 / PDF 61 | Fresnel phase prints exp(iπ signature) rather than exp(iπ signature/4). C-091 derives the consistent Gaussian branch. |
| §2 Fourier/Wigner/Wick normalizations; (2.3)–(2.6), (2.23), (2.97), (2.112)–(2.113) | Several prefactors/conventions do not match the kernel normalization. C-013/C-014/C-017 declare one consistent convention; formulas transferred between conventions require explicit rescaling. |
| (2.89), printed 50 / PDF 67 | Bracket/commutator sign conflicts with $[\hat Q,\hat P]=i\hbar\hat I$. C-016 uses an explicitly declared ordinary-derivative bracket. |
| (3.6), printed 64 / PDF 81 | Missing quadratic factor 1/2; C-018 follows (3.1)/(3.5) and the linear quantum evolution. |
| (3.17), printed 66 / PDF 83 | An extra center-dependent oscillator phase contradicts exact Weyl covariance. C-009/C-019 use the exact covariance. |
| (3.54), printed 77 / PDF 94 | Block exponential fails the zero-generator identity. C-026 uses the consistent squeeze conjugation and its generator derivative. |
| Lemma 25 proof, printed 79 / PDF 96 | Intermediate derivatives have creation/annihilation/conjugation sign mistakes. C-026 records the usable block result and one-mode check. |
| Proposition 36, printed 82–83 / PDF 99–100 | Coordinate squeeze amplitude has incorrect powers and a missing determinant. C-028 derives the normalized result from C-027 and the preceding Gaussian integral. |
| (3.72), printed 83 / PDF 100 | Missing 2^n in the Weyl prefactor, detected at $F=I_{2n}$. C-030 supplies the kernel-normalized correction. |
| §3.4 / §4.1.2 polar decomposition; §3.5 covariant Fourier/Maslov phase | Adjoints and transform normalizations vary. C-029/C-031 now reconstruct the scalar phase in explicit conventions from the vacuum kernel and a Fresnel transform; a continuous lift is still required. |
| Proposition 59(iii), printed 170 / PDF 186 | Clock/shift commutation orientation conflicts with matrices (6.18)/(6.19). C-037 uses $\hat Z\hat X=\omega\hat X\hat Z$ for those displayed operator matrices. |
| §7.4.1 normal form, printed 212–213 / PDF 227–228 | Normal-form parameter/logarithm and highest/lowest fiducial conventions conflict. C-041 derives the chartwise SU(2) factors from the 2-by-2 Gauss decomposition. |
| (8.51), (8.54), printed 256 / PDF 271 | Disk coordinate prints sinh in place of tanh; the diagonal Gauss logarithm needs its minus sign. C-047/C-048 derive the bounded-disk factors. |
| Propositions 108–109, printed 264–265 / PDF 279–280 | **Resolved extraction ambiguity:** rendered Re/Im glyphs are correct. C-049 does not treat this as a source interchange error. The coordinate phase is fixed by the continuous squeeze convention C-028. |
| (8.96), printed 265 / PDF 280 | Right endpoint incorrectly uses β_t; C-050 uses β_s, as required by $\hat U(t,s)$ and identity/composition limits. |
| §11 amplitude conventions | Several Herman–Kluk prefactors disagree in the identity limit; C-054 labels its normalized quadratic re-expression. |
| (15.24) versus earlier inline overlap, printed 440 / PDF 450 | Earlier parentheses misplace the factor 1/2. C-059 uses the product-law-compatible phase. |
| (16.50)–(16.60), printed 476–477 / PDF 486–487 | Kernel trace reflects a label, supertrace does not; rank-one/supermodule coefficient conventions require parity/order care. C-073 separates these from ordinary complex Hilbert-space traces. |
| (16.114), printed 485 / PDF 495 | G_w is a subscript, not G^w. A new parity audit additionally finds that the source positive scaled argument is valid for even symbols; C-076 gives the reflected argument needed for arbitrary odd/even operators and a one-mode counterexample. |
| Proposition 162/(16.147), printed 492–493 / PDF 502–503 | The source coherent label F_t α conflicts with its $\hat U\hat c\hat U^{-1}=F_t\hat c$ orientation. C-082 derives the inverse label F_t⁻¹ α and checks it on one number mode. C-084 flags this when reconstructing the Mehler kernel. |
| Theorem 86/(16.151), (16.158), printed 496–497 / PDF 506–507 | A real sign mismatch remains in the printed Cayley exponent with $\hat U=e^{-it\hat H}$, $F=e^{it\Sigma}$. C-085 gives both the source and the corrected sign, with a one-mode first-derivative/number-spectrum test. |
| Corollary 57/(16.161), printed 498 / PDF 508 | Independent audit resolved the inverse-Cayley minus sign and the additional factor $(-1)^n$ in the amplitude square; C-086 gives the supertrace and ordered-Pfaffian amplitude explicitly. |
| (A.1), printed 531 / PDF 540 | Printed integration variable and Fourier variable are interchanged in one place. C-091 supplies an unambiguous standard Gaussian formula and branch. |
| After Theorem 92, printed 549 / PDF 556 | Global exponential diffeomorphism claim for every connected simply connected group is false. SU(2) is a counterexample; C-094 retains only the valid local claim. |
| Proposition 170/(E.3), printed 558 / PDF 564 | The displayed subleading scalar term fails the constant-function normalization check. C-097 retains the exact Berezin product (E.2). |

The independent correctness audit supplies convention-fixed phase reconstructions for C-029/C-031 and the complete corrected ordered amplitude for C-086. It also repairs C-067–C-069 integration/Pfaffian signs. These local chart results retain their path and invertibility hypotheses; source statements explicitly rejected above remain rejected. Safe routes are the kernel convention C-013, continuous Gaussian propagation C-019, squeeze normalization C-027/C-028, exact finite CAR matrices C-064/C-065, the star equation C-077, and the occupation-number determinant C-079. No QuantumFramework source, documentation, repository code, or Wolfram Language computation was used.


Source page: Combescure–Robert, printed p.468, one-based PDF p.478.




---

<a id="book-g"></a>

## G — Gerry–Knight — Introductory Quantum Optics


Source: Christopher C. Gerry and Peter L. Knight, *Introductory Quantum Optics*, Cambridge University Press, first print edition 2004 (this electronic file carries ©2005). The supplied PDF has 333 pages. This dossier uses only this book and explicitly identified algebraic deductions; no software framework or code is used.

**Citation convention.** “p.” means the number printed in the book; “PDF” means the one-based file page. Arabic pp. 1–317 are present continuously at PDF 17–333. This mapping was checked against actual headers throughout the sequential read, not inferred just from the initial offset. Equation numbers below belong to the source; deductions do not acquire invented source equation numbers. “Source” means the formula is given in the text; “Exercise” means the book asks the reader to establish it; “Derived” means the displayed result is supplied here, rather than presented as a solved result of the book. An entry can contain separately labelled source and derived parts.

**Conventions and scope.** This transcription uses normalized hat notation for operators; faithful source-page images retain the book's original notation. Unless otherwise stated, $[\hat a,\hat a^\dagger]=\hat I$, $\hat n=\hat a^\dagger \hat a$, different modes commute, and $\hat X_1=(\hat a+\hat a^\dagger)/2$, $\hat X_2=(\hat a-\hat a^\dagger)/(2i)$, so $[\hat X_1,\hat X_2]=(i/2)\hat I$. Physical evolution is $\hat U(t)=e^{-i\hat Ht/\hbar}$. The book's squeeze convention is $\hat S(\xi)=e^{(\xi^*\hat a^2-\xi \hat a^{\dagger2})/2}$, $\xi=re^{i\theta}$, with $r\ge0$. It differs by a sign from another common convention. Polynomial commutators can be read on the finite-particle core. Power-series manipulations of unbounded operators require an appropriate common domain or analytic vectors; the book treats these formally, and supplies no general operator-domain theorem. Displacement and squeeze operators themselves have their usual unitary meaning. Trace identities require the indicated products to be trace class; phase-space singular distributions must be interpreted weakly.

| Notation | Type |
|---|---|
| $\hat a,\hat b,\hat c,\hat d;\ \hat n=\hat a^\dagger\hat a;\ \hat q,\hat p,\hat X_j,\hat Y_j$ | Field, number, canonical, and quadrature operators. Fock labels $m,n,N$, moment orders, and mode subscripts remain scalar labels. |
| $\hat H,\hat U,\hat D,\hat S,\hat B,\hat R,\hat\Pi,\hat\rho,\hat\sigma_j,\hat J_j,\hat K_j,\hat M_\ell,\hat I$ | Hamiltonian, evolution, density, atomic, Lie-algebra, Kraus, and identity operators; their quantum matrix representatives keep hats. |
| $\alpha,\beta,\xi,r,\theta,\eta,\mu,\nu;\ P,Q,W,C_N,C_A,C_W;\ T,V$ | Scalar amplitudes/parameters and phase-space functions; $T$ in G-043 is a numerical mixing matrix and $V$ is a numerical covariance matrix. Elsewhere $T$ is scalar transmission. The real mixing coefficients $c,s$, Kerr coupling $K$, entropy $S$, and $\rho_{mn}$ matrix elements are scalar. |
| $\widehat\Lambda_T,\widehat{\mathcal D},\widehat{\operatorname{ad}}_{\hat A}$ | Superoperators acting on density operators. The scalar modified Bessel function $I_q$ is unhatted. Differential operations, traces, determinants, scalar functions $f,F,G$, and time ordering $\mathcal T$ keep standard notation; standalone abstract $\hat F$ is an operator. |

### Full-reading coverage ledger

The sequential full-text pass read every supplied page, including experiments, historical discussion, every exercise, references, bibliographies, and index entries. No body page was missing. Blank PDF pages 1, 4, 8, and 16 were accounted for; PDF 2 is intentionally blank. A terminal empty form-feed split numbered 334 is not a PDF page.

| Material | Printed pages | PDF pages | Full-read result relevant to this dossier |
|---|---:|---:|---|
| Cover, title, copyright, dedication, contents, acknowledgements, blanks | unnumbered; vii–xiii | 1–16 | Complete preliminaries; edition and organization checked. |
| Chapter 1, including references and further reading | 1–9 | 17–25 | Historical and scope discussion; no additional bosonic identity. |
| Chapter 2, all §§2.1–2.7, Problems 1–13, references/bibliography | 10–42 | 26–58 | CCR, Fock calculus, conjugation lemma, BCH exercise, Gibbs exponential, phase-shift obstruction. |
| Chapter 3, all §§3.1–3.8, Problems 1–14, references/bibliography | 43–73 | 59–89 | Coherent/displacement states, Weyl law, reproducing kernel, ordering, characteristic functions, differential projector identities. |
| Chapter 4, all §§4.1–4.10, complete printed problem sequence, references/bibliography | 74–114 | 90–130 | Atom–field algebra, exact operator-valued Jaynes–Cummings propagator, dispersive exponential, multiphoton/Raman exercises. The printed problem numbering skips 10; no page is absent. |
| Chapter 5, all §§5.1–5.4, Problems 1–5, references/bibliography | 115–134 | 131–150 | Normally ordered detection/correlation functions and factorial moments; first-order coherence is weaker than coherent-state factorization. |
| Chapter 6, all §§6.1–6.5, Problems 1–11, references/bibliography | 135–149 | 151–165 | Two-port CCR preservation, beam-splitter exponential, SU(2) exercise algebra and binomial/arcsine states. |
| Chapter 7, all §§7.1–7.9, Problems 1–20, references/bibliography | 150–194 | 166–210 | Single/two-mode squeezing, parametric exponents, Kerr, cat states, higher moments, pair coherent states, frequency conversion. |
| Chapter 8, all §§8.1–8.7, Problems 1–8, references/bibliography | 195–212 | 211–228 | Nonunitary number exponential, quantum jumps, damping characteristic functions, exact loss matrix-element exercise, optical balance. |
| Chapter 9, all §§9.1–9.8, Problems 1–4, references/bibliography | 213–237 | 229–253 | Down-conversion Hamiltonians and low-gain truncation; remaining HOM, erasure, Bell, tunneling, metrology material read in full. |
| Chapter 10, all §§10.1–10.8, Problems 1–7, references/bibliography | 238–262 | 254–278 | Conditional number rotations/parity, trapped-ion displacement exponential and sideband restrictions; experimental passages read in full. |
| Chapter 11, all §§11.1–11.12, Problems 1–10, references/bibliography | 263–293 | 279–309 | Parity exponential, optical phase/controlled phase, conjugated cross-Kerr gate; cryptography, teleportation, algorithms read in full. |
| Appendix A, all §§A.1–A.6, references/bibliography | 294–303 | 310–319 | Trace, partial trace, Schmidt spectrum, entropy, unitary density evolution. |
| Appendix B, complete, bibliography | 304–307 | 320–323 | Conditional projectors needed to distinguish state preparation from unitary exponentiation. |
| Appendix C, complete, references | 308–311 | 324–327 | Dyson ordering and dispersive commutator Hamiltonian; quantized parametric Kerr reduction. |
| Appendix D, complete, references | 312–313 | 328–329 | Continuum CCR, cubic interaction, first-order pair creation and phase matching. |
| Entire index | 314–317 | 330–333 | Cross-checked algebra-related locators; no additional body text hidden behind index-only references. |

**Source-display reminder.** Source images use the book's original notation; this transcription marks operators with hats and leaves scalars and coefficient matrices plain.

**Selective visual pass.** Rendered and inspected PDF 56, 65, 66, 81, 82, 119, 122, 163, 164, 168, 169, 171, 175, 177, 178, 179, 181, 190, 191, 194, 195, 198, 199, 204, 207, 208, 226, 273, 282, 304, and 326 (printed pp. 40, 49, 50, 65, 66, 103, 106, 147, 148, 152, 153, 155, 159, 161, 162, 163, 165, 174, 175, 178, 179, 182, 183, 188, 191, 192, 210, 257, 266, 288, 310). These checks recover fractions, square roots, daggers, signs, ordering, and several actual printed errors listed at the end. Other equations were legible in the page-preserving extraction; the full book was not exhaustively proofread image by image.

### Elementary CCR and functional calculus

<a id="g-001"></a>

#### G-001 — Canonical oscillator realization
**Source:** §2.1, pp. 11–12, PDF 27–28, Eqs. (2.9), (2.12)–(2.18).
$$
[\hat q,\hat p]=i\hbar \hat I,\qquad \hat a=\frac{\omega \hat q+i\hat p}{\sqrt{2\hbar\omega}},\qquad
[\hat a,\hat a^\dagger]=\hat I,\qquad \hat H=\hbar\omega(\hat n+\tfrac12\hat I).
$$
This identifies the algebra and the zero-point scalar that must be retained when an overall phase matters.

<a id="g-002"></a>

#### G-002 — Hadamard conjugation lemma
**Source:** §2.1, p. 13, PDF 29, Eqs. (2.23)–(2.25); **Exercise:** Chapter 2 Problem 3, p. 40, PDF 56 (visually checked).
$$
e^{s\hat A}\hat Be^{-s\hat A}=\sum_{k\ge0}\frac{s^k}{k!}\widehat{\operatorname{ad}}_{\hat A}^k(\hat B),
\quad\widehat{\operatorname{ad}}_{\hat A}(\hat B)=[\hat A,\hat B].
$$
The source writes $s=i\lambda$ in the problem and $s=it/\hbar$ for dynamics. Termination requires vanishing nested commutators; closure on a finite-dimensional operator span can instead permit resummation.

<a id="g-003"></a>

#### G-003 — Central-commutator BCH/disentangling
**Exercise:** Chapter 2 Problem 4, p. 40, PDF 56; **Source:** §3.2, p. 49, PDF 65, Eqs. (3.32)–(3.33), visually checked.
If $[\hat A,[\hat A,\hat B]]=[\hat B,[\hat A,\hat B]]=0$, then
$$
e^{\hat A+\hat B}=e^{\hat A}e^{\hat B}e^{-[\hat A,\hat B]/2}=e^{\hat B}e^{\hat A}e^{[\hat A,\hat B]/2}.
$$
The sufficient condition is commutation with the commutator, not $[\hat A,\hat B]=0$. The theorem is not a general factorization for arbitrary quadratic or higher bosonic polynomials.

<a id="g-004"></a>

#### G-004 — Number commutators and free rotations
**Source:** §2.1, pp. 12–13, PDF 28–29, Eqs. (2.19)–(2.25); §3.3, p. 51, PDF 67, Eq. (3.49).
$$
[\hat n,\hat a]=-\hat a,\quad[\hat n,\hat a^\dagger]=\hat a^\dagger,\quad
e^{i\omega t \hat n}\hat ae^{-i\omega t \hat n}=e^{-i\omega t}\hat a,
\quad e^{-i\omega t(\hat n+\tfrac12\hat I)}|\alpha\rangle=e^{-i\omega t/2}|\alpha e^{-i\omega t}\rangle.
$$
The omitted $1/2$ changes only the scalar phase for a fixed oscillator.

<a id="g-005"></a>

#### G-005 — Fock basis and ladder matrix elements
**Source:** §§2.1–2.2, pp. 14–15, PDF 30–31, Eqs. (2.32), (2.35), (2.38)–(2.43).
$$
\hat a|n\rangle=\sqrt n|n-1\rangle,\quad \hat a^\dagger|n\rangle=\sqrt{n+1}|n+1\rangle,
\quad |n\rangle=\frac{(\hat a^\dagger)^n}{\sqrt{n!}}|0\rangle,
\quad\sum_{n\ge0}|n\rangle\langle n|=\hat I.
$$
These formulas provide a direct coefficientwise verification of ordered expressions and exponential actions.

<a id="g-006"></a>

#### G-006 — Number-shift rules and factorial powers
**Derived from:** G-004–G-005; source anchors §2.1, pp. 13–15, PDF 29–31, Eqs. (2.20), (2.35), (2.38)–(2.40); Chapter 2 Problem 10, p. 41, PDF 57.
$$
\hat a f(\hat n)=f(\hat n+\hat I)\hat a,\quad f(\hat n)\hat a=\hat a f(\hat n-\hat I),\quad
(\hat a^\dagger)^r \hat a^r=\hat n(\hat n-\hat I)\cdots(\hat n-r\hat I+\hat I).
$$
The first two identities hold wherever the relevant functional calculus and products are defined. They are supplied here as deductions, not as numbered theorems of the book.

<a id="g-007"></a>

#### G-007 — Quadrature scaling and vacuum covariance
**Source:** §2.3, p. 17, PDF 33, Eqs. (2.52)–(2.59).
$$
[\hat X_1,\hat X_2]=\frac i2\hat I,\quad \operatorname{Var}\hat X_1\operatorname{Var}\hat X_2\ge\frac1{16},
\quad\operatorname{Var}_{|n\rangle}\hat X_j=\frac{2n+1}{4}.
$$
In particular, the vacuum reference variance is $1/4$, not $1/2$.

<a id="g-008"></a>

#### G-008 — Independent-mode CCR and tensor-product Fock states
**Source:** §2.4, pp. 22–24, PDF 38–40, Eqs. (2.99)–(2.117), (2.122).
$$
[\hat a_j,\hat a_k^\dagger]=\delta_{jk}\hat I,\quad[\hat a_j,\hat a_k]=[\hat a_j^\dagger,\hat a_k^\dagger]=0,
\quad |\{n_j\}\rangle=\prod_j\frac{(\hat a_j^\dagger)^{n_j}}{\sqrt{n_j!}}|\{0\}\rangle.
$$
For finitely many involved modes, commuting independent-mode generators permit exponential factorization. The book's free multimode Hamiltonian is $\sum_j\hbar\omega_j(\hat n_j+\tfrac12\hat I)$.

<a id="g-009"></a>

#### G-009 — Thermal exponential, partition function, and spectrum
**Source:** §2.5, pp. 25–27, PDF 41–43, Eqs. (2.132)–(2.145).
With $x=\hbar\omega/(k_BT)>0$,
$$
Z=\frac{e^{-x/2}}{1-e^{-x}},\quad \hat\rho_{\rm th}=(1-e^{-x})e^{-x\hat n}
=\frac1{1+\bar n}\left(\frac{\bar n}{1+\bar n}\right)^{\hat n},
\quad\bar n=(e^x-1)^{-1}.
$$
The last expression is an operator function of $\hat n$, with eigenvalues $\bar n^m/(1+\bar n)^{m+1}$. Trace convergence fails for $x\le0$ in the infinite oscillator space.

<a id="g-010"></a>

#### G-010 — Thermal factorial moments
**Exercise:** Chapter 2 Problem 10, p. 41, PDF 57; **Source:** §2.5, p. 27, PDF 43, Eqs. (2.147)–(2.148).
$$
\langle(\hat a^\dagger)^r \hat a^r\rangle_{\rm th}=r!\bar n^r,
\quad\langle \hat n^2\rangle=\bar n+2\bar n^2,
\quad\operatorname{Var}\hat n=\bar n(1+\bar n).
$$
**Derived solution of the exercise:** differentiate the geometric series $\sum_{m\ge0}z^m=(1-z)^{-1}$, then set $z=\bar n/(1+\bar n)$.

<a id="g-011"></a>

#### G-011 — Obstruction to a globally canonical Hermitian phase
**Source:** §2.7, pp. 34–36, PDF 50–52, Eqs. (2.184)–(2.203).
The formal polar ansatz $\hat a=e^{i\hat\phi}\sqrt{\hat n}$ and $[\hat n,\hat\phi]=i\hat I$ would require
$$
(m-n)\langle m|\hat\phi|n\rangle=i\delta_{mn},
$$
contradicting the diagonal matrix elements. The lower-bounded number spectrum also prevents a unitary lowering shift on the physical Fock space. Thus one cannot freely use $e^{i\hat\phi}$ as a unitary canonical phase translation.

<a id="g-012"></a>

#### G-012 — Susskind–Glogower one-sided shift
**Source:** §2.7, pp. 36–38, PDF 52–54, Eqs. (2.204)–(2.223).
$$
\hat E=(\hat n+\hat I)^{-1/2}\hat a=\sum_{m\ge0}|m\rangle\langle m+1|,
\quad \hat E\hat E^\dagger=\hat I,\quad \hat E^\dagger \hat E=\hat I-|0\rangle\langle0|.
$$
The formal phase kets $|\phi\rangle=\sum_{m\ge0}e^{im\phi}|m\rangle$ satisfy $\hat E|\phi\rangle=e^{i\phi}|\phi\rangle$ and $(2\pi)^{-1}\int_0^{2\pi}|\phi\rangle\langle\phi|d\phi=\hat I$, but are not normalized Hilbert-space vectors. The word “exponential” in this construction does not imply a unitary exponential of a Hermitian phase.

### Displacements, coherent states, and ordering

<a id="g-013"></a>

#### G-013 — Coherent-state eigenvector and expansion
**Source:** §3.1, pp. 43–44, PDF 59–60, Eqs. (3.2)–(3.10).
$$
\hat a|\alpha\rangle=\alpha|\alpha\rangle,\quad
|\alpha\rangle=e^{-|\alpha|^2/2}\sum_{m\ge0}\frac{\alpha^m}{\sqrt{m!}}|m\rangle,
\quad\langle\alpha|\hat a^\dagger=\alpha^*\langle\alpha|.
$$
All complex $\alpha$ are allowed. The last equality is a left-eigenvector statement, not a right eigenvector of $\hat a^\dagger$.

<a id="g-014"></a>

#### G-014 — No normalizable right eigenvector of creation
**Exercise:** Chapter 3 Problem 1, p. 71, PDF 87. **Derived solution:** if $\hat a^\dagger|\psi\rangle=\lambda|\psi\rangle$, the vacuum coefficient gives $\lambda c_0=0$; recursive comparison gives every $c_m=0$, including $\lambda=0$. Thus no nonzero normalizable right eigenvector exists in the oscillator Fock space. This is an asymmetry that must survive any algebraic simplification.

<a id="g-015"></a>

#### G-015 — Displacement exponential and its two ordered forms
**Source:** §3.2, p. 49, PDF 65, Eqs. (3.30)–(3.39), visually checked.
$$
\hat D(\alpha)=e^{\alpha \hat a^\dagger-\alpha^*\hat a}
=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}e^{-\alpha^*\hat a}
=e^{|\alpha|^2/2}e^{-\alpha^*\hat a}e^{\alpha \hat a^\dagger},
\quad \hat D^\dagger(\alpha)=\hat D(-\alpha).
$$
The scalar signs reverse when normal order is changed to antinormal order.

<a id="g-016"></a>

#### G-016 — Vacuum action of the ordered displacement
**Source:** §3.2, p. 49, PDF 65, Eqs. (3.31), (3.35)–(3.37).
$$
e^{-\alpha^*\hat a}|0\rangle=|0\rangle,\quad
\hat D(\alpha)|0\rangle=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}|0\rangle=|\alpha\rangle.
$$
The unnormalized exponential $e^{\alpha \hat a^\dagger}|0\rangle$ has norm $e^{|\alpha|^2/2}$.

<a id="g-017"></a>

#### G-017 — Weyl composition and noncommuting phases
**Source:** §3.2, p. 50, PDF 66, Eqs. (3.41)–(3.43), visually checked.
$$
\hat D(\alpha)\hat D(\beta)=e^{(\alpha\beta^*-\alpha^*\beta)/2}\hat D(\alpha+\beta),
\quad \hat D(\alpha)|\beta\rangle=e^{i\operatorname{Im}(\alpha\beta^*)}|\alpha+\beta\rangle.
$$
**Derived:** $\hat D(\alpha)\hat D(\beta)=e^{\alpha\beta^*-\alpha^*\beta}\hat D(\beta)\hat D(\alpha)$. The book calls this a semigroup relation, but the displayed operators form the projective Weyl representation of the additive phase plane. Two obvious typographical slips in the intermediate source lines are corrected here (see error ledger).

<a id="g-018"></a>

#### G-018 — Translation of operators and arbitrary ordered functions
**Source:** §7.1, pp. 155, 158, PDF 171, 174, Eqs. (7.23), (7.41).
$$
\hat D^\dagger(\alpha)\hat a\hat D(\alpha)=\hat a+\alpha\hat I,\quad
\hat D^\dagger(\alpha)\hat a^\dagger \hat D(\alpha)=\hat a^\dagger+\alpha^*\hat I.
$$
**Derived:** $\hat D^\dagger F(\hat a,\hat a^\dagger)\hat D=F(\hat a+\alpha\hat I,\hat a^\dagger+\alpha^*\hat I)$, preserving the original order of noncommuting factors. Conversely, $\hat D\hat a\hat D^\dagger=\hat a-\alpha\hat I$.

<a id="g-019"></a>

#### G-019 — Coherent overlap and reproducing measure
**Source:** §3.5, pp. 53–54, PDF 69–70, Eqs. (3.61)–(3.71).
$$
\langle\beta|\alpha\rangle=e^{-(|\alpha|^2+|\beta|^2)/2+\beta^*\alpha},
\quad |\langle\beta|\alpha\rangle|^2=e^{-|\alpha-\beta|^2},
\quad\int_{\mathbb C}|\alpha\rangle\langle\alpha|\frac{d^2\alpha}{\pi}=\hat I.
$$
These are overcomplete states, not an orthogonal basis. The overlap becomes small at large separation, but never vanishes for finite amplitudes.

<a id="g-020"></a>

#### G-020 — Normally ordered coherent matrix elements
**Derived from source:** §3.1, pp. 43–44, PDF 59–60, Eqs. (3.2)–(3.3); §3.7, p. 62, PDF 78, Eqs. (3.105)–(3.106).
$$
\langle\beta|(\hat a^\dagger)^p \hat a^q|\alpha\rangle=(\beta^*)^p\alpha^q\langle\beta|\alpha\rangle,
\quad\langle\beta|{:F(\hat a^\dagger,\hat a):}|\alpha\rangle=F(\beta^*,\alpha)\langle\beta|\alpha\rangle.
$$
For infinite expansions, convergence is an additional requirement; mere formal normal ordering does not guarantee it.

<a id="g-021"></a>

#### G-021 — Coherent statistics and factorial moments
**Source:** §3.1, pp. 45–46, PDF 61–62, Eqs. (3.16), (3.22)–(3.25); **Exercise:** Chapter 3 Problem 6, p. 71, PDF 87.
$$
P_m=e^{-|\alpha|^2}\frac{|\alpha|^{2m}}{m!},\quad
\langle \hat n\rangle=\operatorname{Var}\hat n=|\alpha|^2,\quad
\langle(\hat a^\dagger)^r \hat a^r\rangle=|\alpha|^{2r},\quad\operatorname{Var}\hat X_j=\frac14.
$$
The factorial result is the immediate eigenvector solution of the exercise.

<a id="g-022"></a>

#### G-022 — Linearly driven evolution is a displacement up to phase
**Source:** §3.4, pp. 52–53, PDF 68–69, Eqs. (3.51)–(3.59).
Successive short-time linear drive propagators are displacements. Their time-ordered product has the form
$$
\hat U(T,0)=e^{i\Phi(T)}\hat D(\alpha(T)),\qquad\alpha(T)=\int_0^T u(t)\,dt.
$$
The accumulated scalar $\Phi$ comes from the Weyl composition phase. The interaction's precise current-dependent $u$ follows the book's conventions; using the displacement amplitude itself avoids the inconsistent complex-conjugate labels in intermediate Eq. (3.55).

<a id="g-023"></a>

#### G-023 — Segal–Bargmann entire generating function
**Source:** §3.5, pp. 54–55, PDF 70–71, Eqs. (3.72)–(3.73).
$$
\langle\alpha|\psi\rangle=e^{-|\alpha|^2/2}\psi(\alpha^*),\quad
\psi(z)=\sum_{m\ge0}\psi_m\frac{z^m}{\sqrt{m!}}.
$$
**Derived realization:** $\hat a^\dagger\mapsto z$, $\hat a\mapsto\partial_z$, hence $[\partial_z,z]=\hat I$. The Gaussian inner product follows from G-019. This is an algebraic route to exponential identities; it is not a new source equation.

<a id="g-024"></a>

#### G-024 — Reconstruction from coherent diagonal symbols
**Source:** §3.5, pp. 55–56, PDF 71–72, Eqs. (3.74)–(3.83).
$$
\langle m|\hat F|n\rangle=\frac1{\sqrt{m!n!}}
\left.\partial_{\alpha^*}^{m}\partial_\alpha^n
\left(e^{|\alpha|^2}\langle\alpha|\hat F|\alpha\rangle\right)\right|_{\alpha=\alpha^*=0}.
$$
This is the coefficient-extraction identity wherever the generating expansion converges. Do not promote the source's trace bound for an arbitrary Hermitian $\hat F$ in Eq. (3.81) to a theorem: such a bound needs positivity/trace-class or other sufficient hypotheses.

<a id="g-025"></a>

#### G-025 — Coherent-projector differential correspondences
**Exercise:** Chapter 3 Problem 4, p. 71, PDF 87.
$$
\hat a^\dagger|\alpha\rangle\langle\alpha|=(\alpha^*+\partial_\alpha)|\alpha\rangle\langle\alpha|,
\quad |\alpha\rangle\langle\alpha|\hat a=(\alpha+\partial_{\alpha^*})|\alpha\rangle\langle\alpha|.
$$
Together with $\hat a|\alpha\rangle\langle\alpha|=\alpha|\alpha\rangle\langle\alpha|$ and its right-adjoint counterpart, these turn operator products into phase-space differential operators. Treat $\alpha,\alpha^*$ as independent Wirtinger variables.

<a id="g-026"></a>

#### G-026 — Normal-order notation and optical equivalence
**Source:** §3.7, pp. 59, 62–63, PDF 75, 78–79, Eqs. (3.91), (3.105)–(3.115).
$$
\hat\rho=\int P(\alpha)|\alpha\rangle\langle\alpha|d^2\alpha,
\quad\langle{:G(\hat a^\dagger,\hat a):}\rangle=\int P(\alpha)G(\alpha^*,\alpha)d^2\alpha,
\quad Q(\alpha)=\pi^{-1}\langle\alpha|\hat\rho|\alpha\rangle.
$$
Colons reorder while disregarding commutators: ${:\hat n^2:}=\hat a^{\dagger2}\hat a^2\ne \hat n^2$. Obtaining an equal operator by commuting factors is a different operation and introduces correction terms.

<a id="g-027"></a>

#### G-027 — Singular P representations must be handled weakly
**Source:** §3.7, pp. 60–62, PDF 76–78, Eqs. (3.93)–(3.104).
For $\hat\rho=|\beta\rangle\langle\beta|$, $P(\alpha)=\delta^{(2)}(\alpha-\beta)$. For $\hat\rho=|m\rangle\langle m|$,
$$
P_m(\alpha)=\frac{e^{|\alpha|^2}}{m!}\partial_\alpha^m\partial_{\alpha^*}^m\delta^{(2)}(\alpha).
$$
The latter is a distributional formula under a test-function integral. The Mehta inversion integral (3.96) may not converge as an ordinary integral; the source itself warns of this.

<a id="g-028"></a>

#### G-028 — Three ordered characteristic functions
**Source:** §3.8, pp. 65–66, PDF 81–82, Eqs. (3.128a)–(3.129), (3.132); both pages visually checked.
$$
C_W(\lambda)=\operatorname{Tr}(\hat\rho \hat D(\lambda)),\quad
C_N(\lambda)=\operatorname{Tr}(\hat\rho e^{\lambda \hat a^\dagger}e^{-\lambda^*\hat a}),\quad
C_A(\lambda)=\operatorname{Tr}(\hat\rho e^{-\lambda^*\hat a}e^{\lambda \hat a^\dagger}),
$$
$$
C_W=e^{-|\lambda|^2/2}C_N=e^{|\lambda|^2/2}C_A.
$$
The antinormal definition displayed here follows Eq. (3.132), correcting the actual printed Eq. (3.128c), which puts creation first and changes signs. This distinction is essential.

<a id="g-029"></a>

#### G-029 — Ordered moments from derivatives
**Source:** §3.8, p. 66, PDF 82, Eqs. (3.130a)–(3.130c).
$$
\langle(\hat a^\dagger)^m \hat a^n\rangle=
\left.\partial_\lambda^m\partial_{(-\lambda^*)}^nC_N\right|_0,
\quad
\langle \hat a^m(\hat a^\dagger)^n\rangle=
\left.\partial_\lambda^n\partial_{(-\lambda^*)}^mC_A\right|_0.
$$
Use $C_W$ for the fully symmetrized monomial $\{(\hat a^\dagger)^m \hat a^n\}_W$. Existence of the required moments is needed to differentiate under the trace.

<a id="g-030"></a>

#### G-030 — Cahill–Glauber s ordering
**Source:** §3.8, p. 66, PDF 82, Eq. (3.131).
$$
C(\lambda,s)=e^{s|\lambda|^2/2}\operatorname{Tr}(\hat\rho \hat D(\lambda)),
\quad s=1,0,-1\quad\text{for normal, Weyl, antinormal order.}
$$
This Gaussian scalar is the order-conversion factor for linear bosonic exponentials; it is not a prescription for factorizing every nonlinear exponential.

<a id="g-031"></a>

#### G-031 — Complex Fourier kernels and phase-space distributions
**Source:** §3.8, pp. 66–67, PDF 82–83, Eqs. (3.132)–(3.136).
$$
C_N(\lambda)=\int P(\alpha)e^{\lambda\alpha^*-\lambda^*\alpha}d^2\alpha,
\quad P(\alpha)=\frac1{\pi^2}\int C_N(\lambda)e^{\lambda^*\alpha-\lambda\alpha^*}d^2\lambda.
$$
The same inverse transform of $C_A$ gives $Q$; that of $C_W$ gives $W$. This fixes the sign and $\pi^{-2}$ normalization used throughout the dossier.

<a id="g-032"></a>

#### G-032 — Thermal Gaussian characteristic functions
**Source:** §3.8, p. 67, PDF 83, Eqs. (3.137)–(3.142); the $C_N,C_W$ expressions are immediate source-derived order conversions.
$$
C_A=e^{-(\bar n+1)|\lambda|^2},\quad C_N=e^{-\bar n|\lambda|^2},
\quad C_W=e^{-(\bar n+1/2)|\lambda|^2},
$$
$$
P(\alpha)=\frac{e^{-|\alpha|^2/\bar n}}{\pi\bar n},\quad
Q(\alpha)=\frac{e^{-|\alpha|^2/(1+\bar n)}}{\pi(1+\bar n)}.
$$
The $\bar n\to0$ limit of $P$ is a delta distribution, not an ordinary function.

<a id="g-033"></a>

#### G-033 — Coherent and Fock Q/Wigner functions
**Source:** §3.8, p. 70, PDF 86, Eqs. (3.143)–(3.146); **Exercise:** Chapter 3 Problem 12, p. 72, PDF 88.
$$
Q_{|\beta\rangle}(\alpha)=\pi^{-1}e^{-|\alpha-\beta|^2},\quad
W_{|\beta\rangle}(\alpha)=\frac2\pi e^{-2|\alpha-\beta|^2},
$$
$$
Q_{|m\rangle}(\alpha)=\frac{e^{-|\alpha|^2}|\alpha|^{2m}}{\pi m!},\quad
W_{|m\rangle}(\alpha)=\frac2\pi(-1)^m L_m(4|\alpha|^2)e^{-2|\alpha|^2}.
$$
These test both the displacement convention and the Fourier normalization.

<a id="g-034"></a>

#### G-034 — Displaced parity representation of the Wigner function
**Exercise:** Chapter 3 Problem 14, p. 72, PDF 88; source parity definition in §11.2, p. 266, PDF 282, Eq. (11.3).
$$
W_\psi(\alpha)=\frac2\pi\sum_{m\ge0}(-1)^m
\langle\psi|\hat D(\alpha)|m\rangle\langle m|\hat D^\dagger(\alpha)|\psi\rangle.
$$
**Derived compact form:** $W_{\hat\rho}(\alpha)=\frac2\pi\operatorname{Tr}[\hat\rho \hat D(\alpha)\hat\Pi \hat D^\dagger(\alpha)]$, $\hat\Pi=e^{i\pi \hat n}$; equivalently $\hat D(\alpha)\hat\Pi \hat D^\dagger(\alpha)=\hat D(2\alpha)\hat\Pi$.

<a id="g-035"></a>

#### G-035 — Photoabsorption fixes normal ordering
**Source:** §§5.2, 5.4, pp. 121–122, 130–132, PDF 137–138, 146–148, Eqs. (5.33)–(5.41), (5.87)–(5.105).
$$
G^{(r)}=\langle \hat E^{(-)}(x_1)\cdots \hat E^{(-)}(x_r)
\hat E^{(+)}(x_r)\cdots \hat E^{(+)}(x_1)\rangle.
$$
For a single free mode, $g^{(r)}=\langle(\hat a^\dagger)^r \hat a^r\rangle/\langle \hat n\rangle^r$. It is normal order, not a symmetric electric-field power, that describes ideal absorbing counts. Denominators require nonzero intensity; the vacuum does not have a defined normalized $g^{(2)}$.

<a id="g-036"></a>

#### G-036 — Factorial correlation and number variance
**Source:** §5.4, pp. 131–132, PDF 147–148, Eqs. (5.93)–(5.100); §7.5, pp. 171–172, PDF 187–188, Eqs. (7.105)–(7.107).
$$
g^{(2)}(0)=\frac{\langle \hat n(\hat n-\hat I)\rangle}{\langle \hat n\rangle^2}
=1+\frac{\operatorname{Var}\hat n-\langle \hat n\rangle}{\langle \hat n\rangle^2},
\quad Q_M=\langle \hat n\rangle(g^{(2)}(0)-1).
$$
For coherent/thermal/nonzero Fock states the values are respectively $1,2,1-1/m$. A single free mode has time-independent $g^{(2)}(\tau)$; sub-Poissonian statistics alone are not delayed antibunching. The source's claim in Eq. (5.53) that a number state's first-order correlation factors into mean fields is incorrect: its mean field vanishes although $|g^{(1)}|=1$.

### Atom–field exponentials and passive mode mixing

<a id="g-037"></a>

#### G-037 — Atomic transition algebra coupled to a boson
**Source:** §4.5, pp. 91–93, PDF 107–109, Eqs. (4.91)–(4.104).
$$
\hat\sigma_+=|e\rangle\langle g|,\quad\hat\sigma_-=|g\rangle\langle e|,\quad
[\hat\sigma_+,\hat\sigma_-]=\hat\sigma_3,\quad[\hat\sigma_3,\hat\sigma_\pm]=\pm2\hat\sigma_\pm.
$$
$$
\hat H_{\rm JC}=\hbar\omega \hat n+\frac{\hbar\omega_0}2\hat\sigma_3+
\hbar\lambda(\hat a\hat\sigma_++\hat a^\dagger\hat\sigma_-),\qquad
[\hat H_{\rm JC},\hat n+|e\rangle\langle e|]=0.
$$
The $\pm$ in the second commutator is necessary (the extracted/printed source line suppresses it). The rotating-wave approximation removes $\hat a^\dagger\hat\sigma_+$ and $\hat a\hat\sigma_-$; the exact conservation law belongs to this model.

<a id="g-038"></a>

#### G-038 — Exact resonant operator-valued exponential
**Source:** §4.7, pp. 102–103, PDF 118–119, Eqs. (4.156)–(4.167); PDF 119 visually checked.
In atomic order $(e,g)$,
$$
e^{-i\lambda t(\hat a\hat\sigma_++\hat a^\dagger\hat\sigma_-)}=
\begin{pmatrix}
\cos(\lambda t\sqrt{\hat n+\hat I})&-i\hat a\,\dfrac{\sin(\lambda t\sqrt{\hat n})}{\sqrt{\hat n}}\\
-i\hat a^\dagger\,\dfrac{\sin(\lambda t\sqrt{\hat n+\hat I})}{\sqrt{\hat n+\hat I}}&\cos(\lambda t\sqrt{\hat n})
\end{pmatrix}.
$$
The quotient at zero is defined by continuity as $\lambda t$. Factor placement is part of the formula: $\hat a f(\hat n)\ne f(\hat n)\hat a$. It follows by even/odd powers of the block matrix $\begin{psmallmatrix}0&\hat a\\\hat a^\dagger&0\end{psmallmatrix}$.

<a id="g-039"></a>

#### G-039 — Invariant doublet and detuned exponential
**Source:** §§4.5–4.6, pp. 93–95, 99–102, PDF 109–111, 115–118, Eqs. (4.108)–(4.122), (4.139)–(4.155). **Derived compact propagator:** on $(|e,m\rangle,|g,m+1\rangle)$, after removing $\hbar\omega(m+1/2)\hat I_2$,
$$
\hat h_m/\hbar=\frac12\begin{pmatrix}\Delta&2\lambda\sqrt{m+1}\\2\lambda\sqrt{m+1}&-\Delta\end{pmatrix},
\quad\Omega_m=\sqrt{\Delta^2+4\lambda^2(m+1)},
$$
$$
e^{-i\hat h_mt/\hbar}=\hat I_2\cos(\Omega_mt/2)-i\frac{2\hat h_m}{\hbar\Omega_m}\sin(\Omega_mt/2).
$$
At resonance this reproduces Eq. (4.112). The correct dressed energies are $\hbar\omega(m+1/2)\pm\hbar\Omega_m/2$; Eq. (4.145) omits the factor $1/2$, visually confirmed on printed p.100/PDF116, so its energy formula is corrected here.

<a id="g-040"></a>

#### G-040 — Dispersive interaction as conditional number rotation
**Source:** §4.8, pp. 105–107, PDF 121–123, Eqs. (4.184)–(4.192); §10.6, p. 255, PDF 271, Eq. (10.58); Appendix C, p. 310, PDF 326, Eqs. (C17)–(C18).
$$
\hat H_{\rm eff}=\hbar\chi\bigl(|e\rangle\langle e|+\hat n\hat\sigma_3\bigr),\quad\chi=\lambda^2/\Delta,
$$
$$
\hat U_{\rm eff}=|g\rangle\langle g|\,e^{i\chi t\hat n}
+|e\rangle\langle e|\,e^{-i\chi t(\hat n+\hat I)}.
$$
Thus $|g,\alpha\rangle\mapsto|g,\alpha e^{i\chi t}\rangle$ and $|e,\alpha\rangle\mapsto e^{-i\chi t}|e,\alpha e^{-i\chi t}\rangle$. The excited-state sign is fixed by the Hamiltonian and Eq. (10.58); Eq. (4.186) has the opposite sign. Large detuning is state dependent: $|\lambda|\sqrt{m+1}/|\Delta|\ll1$ on occupied sectors.

<a id="g-041"></a>

#### G-041 — Raman number-controlled atomic rotation
**Exercise:** Chapter 4 Problem 7, p. 111, PDF 127.
The source specifies $\hat H_I=\hbar\lambda \hat n(\hat\sigma_++\hat\sigma_-)$. **Derived solution:**
$$
\hat U_I(t)=\cos(\lambda t\hat n)\hat I_A-i\sin(\lambda t\hat n)\hat\sigma_1
=|+\rangle_x\langle+|e^{-i\lambda t\hat n}+|-\rangle_x\langle-|e^{i\lambda t\hat n}.
$$
Each atomic $\hat\sigma_1$ eigenstate conditions a phase rotation of a coherent field. This is a different polynomial generator from the ordinary Jaynes–Cummings interaction.

<a id="g-042"></a>

#### G-042 — Multiphoton Jaynes–Cummings blocks
**Exercises:** Chapter 4 Problems 8–9, pp. 111–112, PDF 127–128.
The source poses $\hat H=\hbar\eta(\hat a^2\hat\sigma_++\hat a^{\dagger2}\hat\sigma_-)$ and $\hat H=\hbar\eta(\hat a\hat b\hat\sigma_++\hat a^\dagger \hat b^\dagger\hat\sigma_-)$. **Derived:** the first couples $|e,m\rangle\leftrightarrow|g,m+2\rangle$ with matrix element $\hbar\eta\sqrt{(m+1)(m+2)}$; the second couples $|e,m,k\rangle\leftrightarrow|g,m+1,k+1\rangle$ with $\hbar\eta\sqrt{(m+1)(k+1)}$. Each resonant block exponentiates as a two-dimensional sine/cosine rotation. These are solved block reductions supplied here, not source solutions of the exercises.

<a id="g-043"></a>

#### G-043 — Two-port unitarity is required to preserve CCR
**Source:** §6.2, pp. 137–139, PDF 153–155, Eqs. (6.3)–(6.11).
For $\hat{\mathbf a}_{\rm out}=T\hat{\mathbf a}_{\rm in}$,
$$
[\hat a_{{\rm out},j},\hat a_{{\rm out},k}^\dagger]=(TT^\dagger)_{jk}\hat I,
\quad\text{hence }TT^\dagger=I_2.
$$
A vacant input port still contributes a mode in the vacuum. Replacing two outputs by scalar multiples of only one input fails the CCR. This gives an immediate algebraic admissibility check for any proposed passive transformation.

<a id="g-044"></a>

#### G-044 — Beam-splitter exponential and its adjoint action
**Source:** §6.2, p. 139, PDF 155, Eqs. (6.10)–(6.12); **Exercises:** Chapter 6 Problems 1–2, p. 147, PDF 163.
$$
\hat J_1=\frac{\hat a^\dagger \hat b+\hat b^\dagger \hat a}{2},\quad \hat B_1(\vartheta)=e^{i\vartheta \hat J_1},
$$
$$
\hat B_1^\dagger\begin{pmatrix}\hat a\\\hat b\end{pmatrix}\hat B_1=
\begin{pmatrix}\cos(\vartheta/2)&i\sin(\vartheta/2)\\i\sin(\vartheta/2)&\cos(\vartheta/2)\end{pmatrix}
\begin{pmatrix}\hat a\\\hat b\end{pmatrix}.
$$
The general angle is the derived solution of Problem 2; $\vartheta=\pi/2$ is the source's 50:50 example. Hadamard resummation closes on $(\hat a,\hat b)$.

<a id="g-045"></a>

#### G-045 — Real beam-splitter rotation
**Exercise:** Chapter 6 Problem 3, p. 147, PDF 163. **Derived solution:**
$$
\hat J_2=\frac{\hat a^\dagger \hat b-\hat b^\dagger \hat a}{2i},\quad \hat B_2(\vartheta)=e^{i\vartheta \hat J_2},
\quad \hat B_2^\dagger\begin{pmatrix}\hat a\\\hat b\end{pmatrix}\hat B_2=
\begin{pmatrix}c&s\\-s&c\end{pmatrix}\begin{pmatrix}\hat a\\\hat b\end{pmatrix},
\quad c=\cos(\vartheta/2),\ s=\sin(\vartheta/2).
$$
This fixes the generator sign for a real mixing convention rather than importing the $i$-reflection convention unchanged.

<a id="g-046"></a>

#### G-046 — Schwinger SU(2) algebra and fixed-number sectors
**Exercise:** Chapter 6 Problem 4, p. 148, PDF 164.
$$
\hat J_3=(\hat n_a-\hat n_b)/2,\quad \hat J_0=(\hat n_a+\hat n_b)/2,\quad
[\hat J_i,\hat J_j]=i\epsilon_{ijk}\hat J_k,\quad[\hat J_0,\hat J_i]=0,
\quad\hat{\mathbf J}^2=\hat J_0(\hat J_0+\hat I).
$$
With total photon number $N$, the representation has $j=N/2$ and dimension $N+1$. **Derived:** $e^{i\varphi \hat J_3}$ gives opposite phase shifts $\hat a\mapsto e^{i\varphi/2}\hat a$, $\hat b\mapsto e^{-i\varphi/2}\hat b$ under $\hat U^\dagger(\cdot)\hat U$. The algebra is compact SU(2), in contrast with squeezing's noncompact quadratic algebra.

<a id="g-047"></a>

#### G-047 — Binomial state from a number-state input
**Source:** §6.2, p. 142, PDF 158, Eqs. (6.19)–(6.22); **Exercise:** Chapter 6 Problem 5, p. 148, PDF 164.
$$
\hat B_1(\vartheta)|0,N\rangle=
\sum_{k=0}^N\sqrt{\binom Nk}\,[i\sin(\vartheta/2)]^k
[\cos(\vartheta/2)]^{N-k}|k,N-k\rangle.
$$
The general sine/cosine form is the solution of the exercise, equivalent on an appropriate angle branch to the book's tangent form. At 50:50 it reduces to $2^{-N/2}\sum_k i^k\sqrt{\binom Nk}|k,N-k\rangle$.

<a id="g-048"></a>

#### G-048 — Coherent products remain products under passive mixing
**Source:** §6.2, p. 140, PDF 156, Eq. (6.16); §5.3, p. 126, PDF 142, Eqs. (5.68)–(5.69); **Exercise:** Chapter 6 Problem 6, p. 148, PDF 164.
**Derived general solution:**
$$
\hat B_1(\vartheta)|\alpha,\beta\rangle=
|c\alpha+is\beta,\;is\alpha+c\beta\rangle.
$$
It follows by transforming the displacement exponent and factoring commuting different-mode terms. The output is separable for coherent-product inputs. A weak coherent state is not a one-photon Fock state.

<a id="g-049"></a>

#### G-049 — Twin-Fock interference and even outputs
**Source:** §6.2, p. 141, PDF 157, Eqs. (6.17)–(6.18); **Exercises:** Chapter 6 Problems 7–8, p. 148, PDF 164.
At 50:50,
$$
\hat B_1|1,1\rangle=\frac i{\sqrt2}(|2,0\rangle+|0,2\rangle).
$$
**Derived exercise solution, up to the overall phase $i^N$:**
$$
\hat B_1|N,N\rangle\sim\sum_{k=0}^N
\frac{\sqrt{(2k)!(2N-2k)!}}{2^N k!(N-k)!}|2k,2N-2k\rangle.
$$
Only even occupations occur, but for $N>1$ this is not a two-component NOON state. This follows from the polynomial $(\hat a^{\dagger2}+\hat b^{\dagger2})^N$, not a classical binomial distribution of $2N$ independent photons.

<a id="g-050"></a>

#### G-050 — Phase shifter and parity exponential
**Source:** §6.3, pp. 143–144, PDF 159–160, Eqs. (6.26)–(6.31); §11.2, p. 266, PDF 282, Eqs. (11.3)–(11.4).
$$
\hat R(\phi)=e^{i\phi \hat n},\quad \hat R(\phi)|\alpha\rangle=|\alpha e^{i\phi}\rangle,
\quad\hat\Pi=e^{i\pi \hat n}=(-1)^{\hat n},\quad\hat\Pi^2=\hat I,\quad\hat\Pi|\alpha\rangle=|-\alpha\rangle.
$$
Parity's last action is also a direct consequence of the number-state expansion. In the twin-single-photon interferometer the source prints $\langle\hat\Pi\rangle=\cos2\phi$, but its displayed output (11.2) gives $-\cos2\phi$; the doubled phase dependence survives this fringe-sign error.

### Single-mode and two-mode squeezing

<a id="g-051"></a>

#### G-051 — Single-mode quadratic squeeze exponential
**Source:** §7.1, pp. 152–153, PDF 168–169, Eqs. (7.10)–(7.15), visually checked.
$$
\hat S(\xi)=\exp\!\left[\tfrac12(\xi^*\hat a^2-\xi \hat a^{\dagger2})\right],
\quad \hat S^\dagger(\xi)=\hat S(-\xi),\quad|\xi\rangle=\hat S(\xi)|0\rangle.
$$
The quadratic generator preserves parity because it changes number only by two. It is unitary for finite $r$; replacing this generator by higher powers requires a separate existence/domain analysis.

<a id="g-052"></a>

#### G-052 — Single-mode Bogoliubov transform
**Source:** §7.1, p. 153, PDF 169, Eq. (7.12); p. 158, PDF 174, Eq. (7.38).
$$
\hat S^\dagger \hat a\hat S=\mu \hat a-\nu \hat a^\dagger,\quad
\hat S^\dagger \hat a^\dagger \hat S=\mu \hat a^\dagger-\nu^*\hat a,
\quad\mu=\cosh r,\quad\nu=e^{i\theta}\sinh r.
$$
The inverse has plus signs. $\mu^2-|\nu|^2=1$ ensures the transformed commutator is still $\hat I$. Nested commutators close on $(\hat a,\hat a^\dagger)$, yielding hyperbolic functions rather than a scalar-commutator BCH factor.

<a id="g-053"></a>

#### G-053 — Quadratic SU(1,1) closure
**Derived from:** §7.1, p. 152, PDF 168, Eq. (7.10); Chapter 7 Problem 10, p. 191, PDF 207, visually checked.
With canonical normalization
$$
\hat K_+=\tfrac12\hat a^{\dagger2},\quad \hat K_-=\tfrac12\hat a^2,\quad \hat K_0=\tfrac12(\hat n+\tfrac12\hat I),
\quad[\hat K_0,\hat K_\pm]=\pm \hat K_\pm,\quad[\hat K_-,\hat K_+]=2\hat K_0.
$$
The even/odd Fock sectors are invariant. The source's Problem 10 uses real quadratures twice as large as canonical $\hat K_1,\hat K_2$; its displayed normalization should not be silently treated as $[\hat K_1,\hat K_2]=-i\hat K_0$ (see G-089).

<a id="g-054"></a>

#### G-054 — Single-mode normal Gauss factorization
**Derived, not explicitly given in the book:** from the quadratic closure in G-053 and the squeeze convention anchored at §7.1, pp. 152–153, PDF 168–169, Eqs. (7.10), (7.12).
Let $\tau=e^{i\theta}\tanh r$. Then
$$
\hat S(\xi)=e^{-\tau \hat a^{\dagger2}/2}
 e^{-\ln(\cosh r)(\hat n+\tfrac12\hat I)}
 e^{\tau^*\hat a^2/2}.
$$
The factors separately are generally unbounded/nonunitary; the product is interpreted on a suitable common domain and extends to the squeeze unitary. Differentiating the three-factor ansatz and using G-053 verifies the exponent and initial value. This is a derived useful factorization, not an application of the central-commutator theorem, because $[\hat a^2,\hat a^{\dagger2}]=4\hat n+2\hat I$ is noncentral.

<a id="g-055"></a>

#### G-055 — Single-mode reversed Gauss factorization
**Derived, not explicitly given in the book:** §7.1, pp. 152–153, PDF 168–169, Eqs. (7.10), (7.12), via G-054.
$$
\hat S(\xi)=e^{\tau^*\hat a^2/2}
 e^{+\ln(\cosh r)(\hat n+\tfrac12\hat I)}
 e^{-\tau \hat a^{\dagger2}/2}.
$$
It follows from G-054 by taking the inverse of the formula for $\hat S(-\xi)$. Reversing order requires reversing the middle number-exponential sign. The common-domain qualification in G-054 applies.

<a id="g-056"></a>

#### G-056 — Squeezed-vacuum nullifier and unique normalized solution
**Source:** §7.1, pp. 157–161, PDF 173–177, Eqs. (7.35)–(7.39), (7.59)–(7.70).
$$
(\mu \hat a+\nu \hat a^\dagger)|\xi\rangle=0,
\quad |\xi\rangle=\frac1{\sqrt{\cosh r}}
\sum_{m\ge0}\frac{\sqrt{(2m)!}}{2^m m!}(-e^{i\theta}\tanh r)^m|2m\rangle.
$$
Its probabilities are $P_{2m}=(2m)!\tanh^{2m}r/[2^{2m}(m!)^2\cosh r]$, $P_{2m+1}=0$. Although the source says the recurrence has an independent odd solution, its $n=0$ equation forces $C_1=0$, so no nonzero odd solution of this zero-eigenvalue equation exists. Eq. (7.68) has the correct normalization despite a contradictory intermediate sentence.

<a id="g-057"></a>

#### G-057 — General quadrature variance and squeezing axis
**Source:** §7.1, pp. 152–155, PDF 168–171, Eqs. (7.7)–(7.9), (7.16)–(7.21), (7.28).
$$
\hat X(\vartheta)=\frac{\hat ae^{-i\vartheta}+\hat a^\dagger e^{i\vartheta}}2,
\quad\operatorname{Var}_{\xi}\hat X(\vartheta)=\frac14[\cosh2r-\sinh2r\cos(2\vartheta-\theta)].
$$
The compact angle formula is derived from the source's Cartesian covariances. Its minimum is $e^{-2r}/4$ at $\vartheta=\theta/2$; the orthogonal maximum is $e^{2r}/4$. The rotated quadratures saturate $\operatorname{Var}\hat Y_1\operatorname{Var}\hat Y_2=1/16$.

<a id="g-058"></a>

#### G-058 — Covariance matters for rotated squeezed states
**Source:** §7.1, pp. 159–160, PDF 175–176, Eqs. (7.47), (7.53)–(7.58). **Derived explicit covariance:**
$$
\operatorname{Cov}(\hat X_1,\hat X_2)=\tfrac12\langle\{\Delta \hat X_1,\Delta \hat X_2\}\rangle
=-\frac14\sinh2r\sin\theta.
$$
Thus $\det V=1/16$, while $V_{11}V_{22}$ need not equal $1/16$. The Robertson–Schrödinger bound is $\operatorname{Var}\hat A\operatorname{Var}\hat B\ge|\langle[\hat A,\hat B]\rangle|^2/4+\operatorname{Cov}(\hat A,\hat B)^2$. The book defines its $F$ as twice this covariance.

<a id="g-059"></a>

#### G-059 — Displaced squeezed state and low moments
**Source:** §7.1, p. 155, PDF 171, Eqs. (7.22)–(7.28), visually checked.
$$
|\alpha,\xi\rangle=\hat D(\alpha)\hat S(\xi)|0\rangle,
\quad\langle \hat a\rangle=\alpha,\quad\langle \hat a^2\rangle=\alpha^2-e^{i\theta}\sinh r\cosh r,
\quad\langle \hat n\rangle=|\alpha|^2+\sinh^2r.
$$
Displacement changes first moments but preserves the centered covariance matrix. The ordering $\hat D\hat S$ is part of the definition.

<a id="g-060"></a>

#### G-060 — Interchanging displacement and squeezing
**Exercise:** Chapter 7 Problem 4, pp. 190–191, PDF 206–207. **Derived solution using Eq. (7.12):**
$$
\hat S(\xi)\hat D(\beta)=\hat D(\mu\beta-\nu\beta^*)\hat S(\xi),
\quad \hat D(\alpha)\hat S(\xi)=\hat S(\xi)\hat D(\mu\alpha+\nu\alpha^*).
$$
There is no extra scalar phase in this exact conjugation. Thus $\hat S\hat D|0\rangle$ and $\hat D\hat S|0\rangle$ describe the same family after reparametrizing the displacement, but do not have the same amplitude label.

<a id="g-061"></a>

#### G-061 — Bogoliubov eigenvalue definition of the displaced squeezed state
**Source:** §7.1, pp. 158–159, PDF 174–175, Eqs. (7.40)–(7.43).
$$
(\mu \hat a+\nu \hat a^\dagger)|\alpha,\xi\rangle=\gamma|\alpha,\xi\rangle,
\qquad\gamma=\mu\alpha+\nu\alpha^*.
$$
This is the appropriate eigenproblem for a squeezed state. It also fixes the corrected rotated-quadrature eigenvalue: $(\hat Y_1+ie^{-2r}\hat Y_2)|\alpha,\xi\rangle=e^{-r-i\theta/2}\gamma|\alpha,\xi\rangle$, not the source's first expression in Eq. (7.46), which incorrectly substitutes $\alpha$ for $\gamma$.

<a id="g-062"></a>

#### G-062 — Hermite coefficients of a displaced squeezed vacuum
**Source:** §7.1, pp. 162–163, PDF 178–179, Eqs. (7.71)–(7.81), visually checked.
$$
\langle m|\alpha,\xi\rangle=
\frac{e^{-|\alpha|^2/2-\alpha^{*2}e^{i\theta}\tanh r/2}}{\sqrt{m!\cosh r}}
\left(\frac{e^{i\theta}\tanh r}{2}\right)^{m/2}
H_m\!\left(\frac{\gamma}{\sqrt{e^{i\theta}\sinh2r}}\right),
\quad\gamma=\mu\alpha+\nu\alpha^*.
$$
Square-root branches must be chosen consistently; the complete coefficient is branch independent. At $r=0$ use the continuous limit (the displayed separate factors are singular). At $\alpha=0$ it reduces to G-056.

<a id="g-063"></a>

#### G-063 — Coherent overlap with squeezed vacuum and squeezed Q function
**Source:** §7.1, p. 163, PDF 179, Eqs. (7.78)–(7.79). **Derived general overlap/Q, correcting the source's Eq. (7.82), p. 165, PDF 181:**
$$
\langle\beta|\hat S(\xi)|0\rangle=\frac1{\sqrt{\cosh r}}
 e^{-|\beta|^2/2-e^{i\theta}\beta^{*2}\tanh r/2},
$$
$$
Q_{\hat D(\alpha)\hat S(\xi)|0\rangle}(\beta)=\frac1{\pi\cosh r}
\exp\{-|\delta|^2-\tanh r\,\operatorname{Re}(e^{-i\theta}\delta^2)\},\quad\delta=\beta-\alpha.
$$
This follows by displaced overlap and G-056. It is centered at $\alpha$; the source's printed (7.82) is not the Q function of its DS-defined state for general $\alpha$.

<a id="g-064"></a>

#### G-064 — Squeezed Gaussian Wigner function with consistent scaling
**Derived correction of:** §7.1, p. 165, PDF 181, Eq. (7.83), with quadrature convention fixed by Eqs. (7.18), (7.28) and Wigner normalization (3.145).
Let $x+iy=e^{-i\theta/2}(\beta-\alpha)$. Then
$$
W(\beta)=\frac2\pi\exp[-2e^{2r}x^2-2e^{-2r}y^2].
$$
It integrates to one and has $\operatorname{Var}x=e^{-2r}/4$. The actual printed Eq. (7.83) has incompatible factors and axis signs under this convention; it is not quoted as a valid formula. A nonnegative Wigner function does not by itself imply a classical P function.

<a id="g-065"></a>

#### G-065 — Number variance of a displaced squeezed vacuum
**Source:** §7.4, pp. 170–171, PDF 186–187, Eqs. (7.100)–(7.104).
$$
\operatorname{Var}\hat n=|\mu\alpha-\nu\alpha^*|^2+2\sinh^2r\cosh^2r.
$$
For $\arg\alpha=\theta/2$, the coherent contribution is $|\alpha|^2e^{-2r}$; for $\arg\alpha=\theta/2+\pi/2$, it is $|\alpha|^2e^{2r}$. Number squeezing requires comparison with $\langle \hat n\rangle=|\alpha|^2+\sinh^2r$, not only with $|\alpha|^2$. For squeezed vacuum, $\operatorname{Var}\hat n=2\bar n(\bar n+1)$.

<a id="g-066"></a>

#### G-066 — Degenerate parametric Hamiltonian and squeeze evolution
**Source:** §7.2, pp. 165–167, PDF 181–183, Eqs. (7.84)–(7.88).
$$
\hat H=\hbar\omega \hat n_a+\hbar\omega_P\hat n_b+i\hbar\chi^{(2)}(\hat a^2\hat b^\dagger-\hat a^{\dagger2}\hat b).
$$
Under a strong, undepleted coherent pump $\hat b\to\beta e^{-i\omega_Pt}$ and resonance $\omega_P=2\omega$,
$$
\hat H_I=i\hbar(\eta^*\hat a^2-\eta \hat a^{\dagger2}),\quad
\hat U_I(t)=e^{t(\eta^*\hat a^2-\eta \hat a^{\dagger2})}=\hat S(2\eta t),\quad\eta=\chi^{(2)}\beta.
$$
The factor 2 in the squeeze parameter is required by the single-mode definition. Detuning generally restores explicit time dependence; one cannot then exponentiate a fixed generator without accounting for time ordering or changing frame.

<a id="g-067"></a>

#### G-067 — Four-wave mixing gives the same quadratic exponent after approximation
**Source:** §7.2, p. 167, PDF 183, Eq. (7.89) and following paragraph.
$$
\hat H=\hbar\omega \hat n_a+\hbar\omega \hat n_b+i\hbar\chi^{(3)}(\hat a^2\hat b^{\dagger2}-\hat a^{\dagger2}\hat b^2).
$$
The undepleted-pump substitution gives the same resonant $\hat H_I$ as G-066, with $\eta=\chi^{(3)}\beta^2$. The fully quantized generator is quartic and should not be assigned the quadratic Bogoliubov formula before the pump approximation.

<a id="g-068"></a>

#### G-068 — Two-mode squeeze exponential and Bogoliubov action
**Source:** §7.7, pp. 182–183, PDF 198–199, Eqs. (7.156)–(7.159), visually checked.
$$
\hat S_2(\xi)=e^{\xi^*\hat a\hat b-\xi \hat a^\dagger \hat b^\dagger},\quad
\hat S_2^\dagger \hat a\hat S_2=\mu \hat a-\nu \hat b^\dagger,\quad
\hat S_2^\dagger \hat b\hat S_2=\mu \hat b-\nu \hat a^\dagger.
$$
Here there is no $1/2$ in the exponent, unlike the degenerate single-mode case. The two independent input CCR imply the transformation preserves all output CCR.

<a id="g-069"></a>

#### G-069 — Two-mode SU(1,1) and number-difference invariant
**Source:** §7.7, pp. 183–184, PDF 199–200, Eqs. (7.164)–(7.170). **Derived algebra:**
$$
\hat K_+=\hat a^\dagger \hat b^\dagger,\quad \hat K_-=\hat a\hat b,\quad \hat K_0=(\hat n_a+\hat n_b+\hat I)/2,
\quad[\hat K_0,\hat K_\pm]=\pm \hat K_\pm,\quad[\hat K_-,\hat K_+]=2\hat K_0,
$$
$$
[\hat n_a-\hat n_b,\hat K_\pm]=[\hat n_a-\hat n_b,\hat S_2]=0.
$$
Thus the squeeze exponential acts independently in each number-difference sector; the vacuum lies in difference zero.

<a id="g-070"></a>

#### G-070 — Two-mode normal Gauss factorization
**Derived, not explicitly displayed in the book:** source definition §7.7, p. 182, PDF 198, Eq. (7.156), and the algebra in G-069.
$$
\hat S_2(\xi)=e^{-\tau \hat a^\dagger \hat b^\dagger}
 e^{-\ln(\cosh r)(\hat n_a+\hat n_b+\hat I)}e^{\tau^*\hat a\hat b},
\quad\tau=e^{i\theta}\tanh r.
$$
This uses noncentral SU(1,1) disentangling. Acting on vacuum gives the source's Eq. (7.169); differentiating the factors gives the defining squeeze generator. Individual factors need common-domain care, as in G-054.

<a id="g-071"></a>

#### G-071 — Two-mode reversed factorization and supermode decomposition
**Derived, not explicitly displayed in the book:** source anchors §7.7, pp. 182–183, PDF 198–199, Eqs. (7.156), (7.158)–(7.159).
$$
\hat S_2=e^{\tau^*\hat a\hat b}e^{+\ln(\cosh r)(\hat n_a+\hat n_b+\hat I)}e^{-\tau \hat a^\dagger \hat b^\dagger}.
$$
For $\hat c=(\hat a+\hat b)/\sqrt2$, $\hat d=(\hat a-\hat b)/\sqrt2$,
$$
\hat S_{2,ab}(\xi)=\hat S_c(\xi)\hat S_d(-\xi).
$$
The factors commute because $\hat c,\hat d$ are independent modes and $\hat a\hat b=(\hat c^2-\hat d^2)/2$. This does not contradict the source's statement that $\hat S_2$ is not a product of squeezers in the original $\hat a,\hat b$ modes.

<a id="g-072"></a>

#### G-072 — Twin-beam vacuum expansion and nullifiers
**Source:** §7.7, pp. 183–184, PDF 199–200, Eqs. (7.165)–(7.170).
$$
|\xi\rangle_2=\frac1{\cosh r}\sum_{m\ge0}(-e^{i\theta}\tanh r)^m|m,m\rangle,
$$
$$
(\mu \hat a+\nu \hat b^\dagger)|\xi\rangle_2=(\mu \hat b+\nu \hat a^\dagger)|\xi\rangle_2=0,
\quad(\hat n_a-\hat n_b)|\xi\rangle_2=0.
$$
Both nullifiers plus normalization identify the vacuum-generated state. One nullifier alone is not enough to select a unique two-mode state.

<a id="g-073"></a>

#### G-073 — Collective squeezed quadrature and thermal individual modes
**Source:** §7.7, pp. 182–186, PDF 198–202, Eqs. (7.158)–(7.162), (7.171)–(7.183).
For $\hat c=(\hat a+\hat b)/\sqrt2$, the collective quadratures are the ordinary $\hat c$-quadratures. At $\theta=0$, their variances are $e^{-2r}/4$ and $e^{2r}/4$. Individual modes have
$$
\langle \hat n_a\rangle=\langle \hat n_b\rangle=N=\sinh^2r,\quad
\operatorname{Var}\hat n_a=\operatorname{Var}\hat n_b=N(N+1),
\quad\operatorname{Cov}(\hat n_a,\hat n_b)=N(N+1).
$$
The source's (7.161) prints $\sin r$ where $\sinh r$ is required; the corrected general collective variance is G-057. The correlation coefficient is one for $r>0$, undefined at vacuum where both variances vanish.

<a id="g-074"></a>

#### G-074 — Reduced Gibbs state and entropy of a twin beam
**Source:** §7.7, pp. 186–187, PDF 202–203, Eqs. (7.181)–(7.185).
$$
\hat\rho_a=\operatorname{sech}^2r\sum_{m\ge0}\tanh^{2m}r|m\rangle\langle m|,
\quad e^{-\hbar\omega_a/(k_BT_{\rm eff})}=\tanh^2r,
\quad T_{\rm eff}=\frac{\hbar\omega_a}{2k_B\ln\coth r}.
$$
**Derived evaluation of the source entropy sum:** $S(\hat\rho_a)=(N+1)\ln(N+1)-N\ln N$, $N=\sinh^2r$. The pure two-mode state has mixed thermal marginals; a thermal-looking single-mode distribution does not establish a mixed global state.

<a id="g-075"></a>

#### G-075 — Nondegenerate parametric evolution and its parameter
**Source:** §7.7, pp. 186–187, PDF 202–203, Eqs. (7.186)–(7.190).
$$
\hat H=\hbar\omega_a\hat n_a+\hbar\omega_b\hat n_b+\hbar\omega_P\hat n_c
+i\hbar\chi^{(2)}(\hat a\hat b\hat c^\dagger-\hat a^\dagger \hat b^\dagger \hat c).
$$
The strong undepleted pump and $\omega_P=\omega_a+\omega_b$ give
$$
\hat H_I=i\hbar(\eta^*\hat a\hat b-\eta \hat a^\dagger \hat b^\dagger),\quad \hat U_I=\hat S_2(\eta t),\quad\eta=\chi^{(2)}\gamma.
$$
Unlike G-066, the two-mode squeeze parameter is $\eta t$, not $2\eta t$. A quantized pump creates a trilinear interaction, so the exact three-mode dynamics is not the undepleted quadratic exponential.

### Kerr generators, cat states, and conditional exponentials

<a id="g-076"></a>

#### G-076 — Self-Kerr number-squared exponential
**Source:** §7.6, p. 179, PDF 195, Eqs. (7.140)–(7.143), visually checked.
$$
\hat H_I=\hbar K \hat n^2,\quad \hat U_K(t)=e^{-iKt \hat n^2},\quad
\hat U_K(t)|\alpha\rangle=e^{-|\alpha|^2/2}\sum_{m\ge0}\frac{\alpha^m}{\sqrt{m!}}e^{-iKtm^2}|m\rangle.
$$
The period is $2\pi/K$ for $K>0$; at $t=\pi/K$ the action on a coherent state is parity. A scalar replacement $\hat n\to|\alpha|^2$ would lose the relative phases responsible for cat formation.

<a id="g-077"></a>

#### G-077 — Exact half-revival cat identity
**Source:** §7.6, p. 179, PDF 195, Eq. (7.143); Eq. (7.142) contains a sign error, visually confirmed.
$$
e^{-i\pi \hat n^2/2}=\frac{e^{-i\pi/4}}{\sqrt2}(\hat I+i\hat\Pi),
\quad e^{-i\pi \hat n^2/2}|\alpha\rangle=
\frac{e^{-i\pi/4}}{\sqrt2}(|\alpha\rangle+i|-\alpha\rangle).
$$
The operator identity is the coefficientwise derivation supplied here: the number phase is $1$ for even $m$ and $-i$ for odd $m$, whereas the intermediate source Eq. (7.142) prints $+i$. The final source Eq. (7.143) is consistent with $-i$.

<a id="g-078"></a>

#### G-078 — Normally ordered Kerr generator differs by a rotation
**Source:** §7.6, p. 179, PDF 195, Eqs. (7.144)–(7.146).
$$
\hat H_I=\hbar K \hat a^{\dagger2}\hat a^2=\hbar K \hat n(\hat n-\hat I),
\quad e^{-iKt \hat n(\hat n-\hat I)}=e^{-iKt \hat n^2}e^{iKt\hat n}.
$$
Thus its half-revival cat has amplitudes $\pm i\alpha$. **Derived:** since $m(m-1)$ is even, this unitary is exactly periodic with period $\pi/K$; its difference from $\hat n^2$ dynamics is a rotation, not a scalar phase.

<a id="g-079"></a>

#### G-079 — General number-function conjugation and Kerr Heisenberg solution
**Exercise:** Chapter 7 Problem 8(a), p. 191, PDF 207. **Derived solution, using G-006:** for real $f$,
$$
e^{itf(\hat n)}\hat ae^{-itf(\hat n)}=\hat a\,e^{-it[f(\hat n)-f(\hat n-\hat I)]}.
$$
For $f(\hat n)=K \hat n(\hat n-\hat I)$,
$$
\hat a(t)=\hat a\,e^{-2iKt(\hat n-\hat I)}=e^{-2iKt\hat n}\hat a,
\quad\dot{\hat a}(t)=-2iK \hat n\,\hat a(t).
$$
Both ordered forms are equal, but $\hat a e^{-2iKt\hat n}$ alone is missing the required scalar shift.

<a id="g-080"></a>

#### G-080 — Kerr moments and unchanged photon statistics
**Exercise:** Chapter 7 Problem 8(b)–(c), p. 191, PDF 207. **Derived solution:** for the normally ordered Kerr generator of G-078 and $N=|\alpha|^2$,
$$
\langle \hat a^m\rangle_t=\alpha^m e^{-iKtm(m-1)}
\exp\{N(e^{-2imKt}-1)\},\quad m\ge0.
$$
All number probabilities remain Poissonian because $[\hat U_K,\hat n]=0$. Quadrature variances follow from
$$
\operatorname{Var}\hat X(\vartheta)=\frac14\{1+2[\langle \hat n\rangle-|\langle \hat a\rangle|^2]
+2\operatorname{Re}[e^{-2i\vartheta}(\langle \hat a^2\rangle-\langle \hat a\rangle^2)]\}.
$$
This supplies an exact analytic route to the exercise's short-time squeezing question without asserting that every time is squeezed.

<a id="g-081"></a>

#### G-081 — Cross-Kerr conditional rotation
**Source:** §7.6, pp. 180–181, PDF 196–197, Eqs. (7.147)–(7.152).
$$
\hat H_{\rm CK}=\hbar K \hat n_a\hat n_b,\quad \hat U_{\rm CK}=e^{-iKt\hat n_a\hat n_b},
\quad \hat U_{\rm CK}|\alpha\rangle_a|m\rangle_b
=|\alpha e^{-iKtm}\rangle_a|m\rangle_b.
$$
**Derived adjoint action:** $\hat U_{\rm CK}^\dagger \hat a\hat U_{\rm CK}=\hat a e^{-iKt\hat n_b}$, and similarly $\hat b\mapsto \hat b e^{-iKt\hat n_a}$. Both number operators are invariants. This is not a linear two-mode Bogoliubov transformation.

<a id="g-082"></a>

#### G-082 — General coherent-cat normalization and eigenoperator
**Source:** §7.6, pp. 174–175, PDF 190–191, Eqs. (7.112)–(7.119).
$$
|C_\Phi\rangle=\frac{|\alpha\rangle+e^{i\Phi}|-\alpha\rangle}
{\sqrt{2(1+e^{-2|\alpha|^2}\cos\Phi)}},\qquad
 \hat a^2|C_\Phi\rangle=\alpha^2|C_\Phi\rangle.
$$
The reciprocal square root corrects the source's Eq. (7.113). Even/odd cats use $\Phi=0,\pi$; the Yurke–Stoler state uses $\Phi=\pi/2$, with normalization exactly $1/\sqrt2$ for every $\alpha$. The odd cat at $\alpha=0$ is defined only as a limit (giving $|1\rangle$, up to phase).

<a id="g-083"></a>

#### G-083 — Cat parity projectors and photon probabilities
**Source:** §7.6, pp. 175–176, PDF 191–192, Eqs. (7.115)–(7.123); §11.2, p. 266, PDF 282, Eq. (11.3). **Derived corrected compact formulas:**
$$
|C_\pm\rangle=\frac{(\hat I\pm\hat\Pi)|\alpha\rangle}{\sqrt{2(1\pm e^{-2N})}},\quad
P_m=\frac{e^{-N}N^m}{m!}
\frac{|1+e^{i\Phi}(-1)^m|^2}{2(1+e^{-2N}\cos\Phi)},\quad N=|\alpha|^2.
$$
Even/odd cats contain only even/odd number states; the Yurke–Stoler cat has a Poisson distribution despite phase-space interference. This corrected common expression avoids missing factors and inconsistent denominators in the printed (7.122)–(7.123).

<a id="g-084"></a>

#### G-084 — Cat Wigner interference with correct odd normalization
**Source:** §7.6, pp. 177–178, PDF 193–194, Eqs. (7.136)–(7.139). **Derived unified form**, for real $\alpha$, $\beta=x+iy$:
$$
W_{C_\Phi}(x,y)=\frac{e^{-2(x-\alpha)^2-2y^2}+e^{-2(x+\alpha)^2-2y^2}
+2e^{-2x^2-2y^2}\cos(4\alpha y+\Phi)}
{\pi(1+e^{-2\alpha^2}\cos\Phi)}.
$$
This reproduces source (7.137) and (7.139), and corrects the plus denominator in the odd-state (7.138) to a minus. The incoherent equal mixture contains just the two Gaussian terms with prefactor $1/\pi$. The interference term is not suppressed by a factor $e^{-2\alpha^2}$ at the phase-space midpoint.

<a id="g-085"></a>

#### G-085 — Conditional projectors produce cats but are not unitary exponentials
**Source:** §7.6, pp. 180–182, PDF 196–198, Eqs. (7.150)–(7.155); §10.5, pp. 250–251, PDF 266–267, Eqs. (10.33)–(10.41); Appendix B, pp. 305–306, PDF 321–322, Eqs. (B5)–(B6), (B12).
The book combines a conditional rotation with an atomic or optical interference measurement. At a phase shift of $\pi$, the conditional field operator is proportional to $\hat I\pm\hat\Pi$, yielding even/odd cats after normalization. For a general rank-one outcome $|q\rangle$, the unnormalized remaining state is $\langle q|\Psi\rangle$; its norm squared is the outcome probability. It is a nonunitary conditional preparation, not a factorization of a single field-only unitary.

<a id="g-086"></a>

#### G-086 — Cross-Kerr conjugation into a conditional SU(2) rotation
**Source:** §11.11, pp. 287–289, PDF 303–305, Eqs. (11.44)–(11.52), especially p. 288, PDF 304 (visually checked).
For a consistently named control mode $\hat c$, set $\hat B=e^{i\pi(\hat a^\dagger \hat b+\hat a\hat b^\dagger)/4}$ and $\hat U_K=e^{i\eta \hat n_b \hat n_c}$. **Derived from the source construction:**
$$
\hat B^\dagger \hat U_K\hat B=
 e^{i\eta \hat n_c(\hat n_a+\hat n_b)/2}
 e^{(\eta \hat n_c/2)(\hat a^\dagger \hat b-\hat b^\dagger \hat a)}.
$$
The generators commute since total target number commutes with the SU(2) rotation and with $\hat n_c$. The source's Eq. (11.46) uses $\hat d$ while Eq. (11.49) uses $\hat c$; the identity above states the required relabelling explicitly. The source's extra $i$ factors in (11.51) do not follow from its displayed real rotation alone, so they are not copied here.

### Higher-order moments, pair states, and approximation restrictions

<a id="g-087"></a>

#### G-087 — Linear quadrature exponential generates normal-order contractions
**Source:** §7.8, p. 188, PDF 204, Eqs. (7.191)–(7.193), visually checked.
With $\Delta \hat X=\hat X-\langle \hat X\rangle\hat I$ and c-number $y$,
$$
e^{y\Delta \hat X}={:e^{y\Delta \hat X}:}\,e^{y^2/8}.
$$
The $e^{y^2/8}$ factor follows from the central commutator of the creation and annihilation parts of $y\hat X$. This is an ordinary algebraic order conversion; the exponential of a high-degree field polynomial does not inherit this factor unchanged.

<a id="g-088"></a>

#### G-088 — All centered quadrature moments in terms of normal moments
**Source:** §7.8, pp. 188–189, PDF 204–205, Eqs. (7.194)–(7.200); **Exercise:** Chapter 7 Problem 12, p. 191, PDF 207. **Derived clean sum form:**
$$
\langle(\Delta \hat X)^N\rangle=
\sum_{k=0}^{\lfloor N/2\rfloor}\frac{N!}{(N-2k)!k!8^k}
\langle{: (\Delta \hat X)^{N-2k}:}\rangle.
$$
The zeroth moment is one. For a coherent state, all positive-degree centered normal moments vanish; hence $\langle(\Delta \hat X)^{2\ell}\rangle=(2\ell-1)!!/4^\ell$. Even-order Hong–Mandel squeezing compares with this value. Fourth-order squeezing requires $\langle{: (\Delta \hat X)^4:}\rangle+\tfrac32\langle{: (\Delta \hat X)^2:}\rangle<0$, not merely a negative fourth normal moment. Odd centered moments can be negative, contrary to a sentence on p. 189; they are not a positive noise benchmark.

<a id="g-089"></a>

#### G-089 — Squared-field quadratures and their noncentral commutator
**Exercise:** Chapter 7 Problem 10, p. 191, PDF 207, visually checked. The source defines
$$
\hat Q_1=(\hat a^{\dagger2}+\hat a^2)/2,\quad \hat Q_2=(\hat a^{\dagger2}-\hat a^2)/(2i),
\quad \hat K_3=(\hat n+\tfrac12\hat I)/2.
$$
**Derived commutators in exactly that normalization:**
$$
[\hat Q_1,\hat Q_2]=-4i\hat K_3,\quad[\hat K_3,\hat Q_1]=i\hat Q_2,\quad[\hat Q_2,\hat K_3]=i\hat Q_1.
$$
Thus $\operatorname{Var}\hat Q_1\operatorname{Var}\hat Q_2\ge4\langle \hat K_3\rangle^2$. Dividing $\hat Q_1,\hat Q_2$ by two yields the canonical SU(1,1) real generators associated with G-053. Squared-field squeezing uses a state-dependent commutator; it is distinct from squeezing a higher moment of the ordinary quadrature.

<a id="g-090"></a>

#### G-090 — Pair coherent states with fixed number difference
**Exercise:** Chapter 7 Problems 15, 17, p. 192, PDF 208, visually checked. The source imposes $\hat a\hat b|\eta,q\rangle=\eta|\eta,q\rangle$ and $(\hat n_a-\hat n_b)|\eta,q\rangle=q|\eta,q\rangle$. **Derived solution for integer $q\ge0$:**
$$
|\eta,q\rangle=\left[\sum_{m\ge0}\frac{|\eta|^{2m}}{m!(m+q)!}\right]^{-1/2}
\sum_{m\ge0}\frac{\eta^m}{\sqrt{m!(m+q)!}}|m+q,m\rangle.
$$
The normalization sum is $|\eta|^{-q}I_q(2|\eta|)$ for $\eta\ne0$, continued at zero. These are eigenstates of $\hat a\hat b$, whereas the squeezed vacuum is annihilated by Bogoliubov mixtures; do not identify the two families.

<a id="g-091"></a>

#### G-091 — Normally ordered intermode Cauchy–Schwarz test
**Exercise:** Chapter 7 Problem 16, p. 192, PDF 208.
For a nonnegative P distribution,
$$
\langle \hat a^{\dagger2}\hat a^2\rangle\langle \hat b^{\dagger2}\hat b^2\rangle
\ge\langle \hat n_a\hat n_b\rangle^2.
$$
**Derived check for G-072:** with $N=\sinh^2r>0$, the left side is $4N^4$, while the right side is $(2N^2+N)^2>4N^4$. Thus the twin beam violates the classical normally ordered bound. This is not a violation of the Hilbert-space Cauchy–Schwarz theorem; it is failure of the positive classical P premise.

<a id="g-092"></a>

#### G-092 — Frequency conversion is an SU(2) exponential, with a source sign caveat
**Exercise:** Chapter 7 Problem 20, p. 192, PDF 208, visually checked.
The source specifies $\hat H_I=-i\hbar \hat G$, $\hat G=\eta \hat b^\dagger \hat c-\eta^*\hat b\hat c^\dagger$, but then prints $\hat U=e^{+t\hat G}$. **Derived correction:** the stated Hamiltonian requires $\hat U=e^{-t\hat G}$. For $\eta=ge^{i\varphi}$, $g\ge0$,
$$
\hat U^\dagger \hat b\hat U=\hat b\cos gt-e^{i\varphi}\hat c\sin gt,\quad
\hat U^\dagger \hat c\hat U=\hat c\cos gt+e^{-i\varphi}\hat b\sin gt.
$$
It preserves total number, leaves vacuum unchanged, and yields binomial number splitting and product coherent outputs. If one chooses the source's printed $e^{+t\hat G}$ instead, all sine signs reverse and the Hamiltonian is $+i\hbar \hat G$.

<a id="g-093"></a>

#### G-093 — Photon subtraction and addition are conditional operations
**Exercises:** Chapter 2 Problem 7, p. 41, PDF 57; Chapter 7 Problem 5, p. 191, PDF 207. **Derived solutions:**
$$
|\psi_-\rangle=\frac{\hat a|\psi\rangle}{\sqrt{\langle \hat n\rangle}},\quad
\langle \hat n\rangle_- =\frac{\langle \hat n(\hat n-\hat I)\rangle}{\langle \hat n\rangle}
=\langle \hat n\rangle-1+\frac{\operatorname{Var}\hat n}{\langle \hat n\rangle}.
$$
For $\langle \hat n\rangle>0$, subtraction need not reduce the conditional mean by one. For coherent input,
$$
|\alpha;+\rangle=\frac{\hat a^\dagger|\alpha\rangle}{\sqrt{1+|\alpha|^2}}
=\frac{\hat D(\alpha)(|1\rangle+\alpha^*|0\rangle)}{\sqrt{1+|\alpha|^2}}.
$$
Neither operation is a unitary displacement.

<a id="g-094"></a>

#### G-094 — Broadband collective CCR and finite bandwidth
**Source:** §7.9, pp. 189–190, PDF 205–206, Eqs. (7.201)–(7.206).
For $\hat A(t)=\sum_{k\in\mathcal B}l_k \hat a_ke^{-i\omega_kt}$,
$$
[\hat A(t),\hat A^\dagger(t)]=C\hat I,\quad C=\sum_{k\in\mathcal B}|l_k|^2,\quad
[\hat X_1^{(c)},\hat X_2^{(c)}]=(iC/2)\hat I,
\quad\operatorname{Var}\hat X_i^{(c)}<C/4\ \text{signals squeezing.}
$$
The mode set/response bandwidth must make $C$ finite; a pointlike, unfiltered continuum field needs smearing or regularization. Independent-mode CCR alone do not make a formally infinite sum a well-defined operator.

<a id="g-095"></a>

#### G-095 — Low-gain down-conversion truncation and multipair warning
**Source:** §9.1, pp. 214–217, PDF 230–233, Eqs. (9.1)–(9.11); Chapter 9 Problem 1, p. 235, PDF 251.
The source Hamiltonian has the form $\hat H_I=\hbar(\eta \hat a_s^\dagger \hat a_i^\dagger+\eta^*\hat a_s\hat a_i)$. **Derived consistent expansion** for $\mu=\eta t$:
$$
\hat U|0,0\rangle=(1-|\mu|^2/2)|0,0\rangle-i\mu|1,1\rangle-\mu^2|2,2\rangle+O(t^3).
$$
The book deliberately drops the two-pair term in its weak-source discussion. That truncation is not valid when multipair probabilities or exact normalization at higher order matter. The type-II polarization Hamiltonian in Eq. (9.6) is a sum of commuting pair-creation channels; its coefficient normalization must be tracked rather than inferred from the inconsistently normalized intermediate (9.10)–(9.11).

### Loss exponentials, trace identities, and dispersive reductions

<a id="g-096"></a>

#### G-096 — Nonunitary number exponential and coherent-state normalization
**Source:** §8.3, pp. 200–202, PDF 216–218, Eqs. (8.1)–(8.11). The no-jump effective Hamiltonian is $\hat H_{\rm eff}=\hat H-i\hbar\gamma \hat n/2$. **Derived exact normalization**, for real $s\ge0$,
$$
e^{-s\hat n}|\alpha\rangle=\exp[-\tfrac12|\alpha|^2(1-e^{-2s})]|\alpha e^{-s}\rangle,
\quad\|e^{-s\hat n}|\alpha\rangle\|^2=e^{-|\alpha|^2(1-e^{-2s})}.
$$
After normalization, no-jump evolution with $s=\gamma\delta t/2$ shrinks the amplitude. The scalar prefactor must not be deleted in an unnormalized no-jump probability calculation; the source's (8.11) omits it. The source's normalized ket (8.8) also divides by the expectation of $e^{-\gamma\delta t\hat n}$ where its square root is required. With the scalar norm above, the normalized state and the no-jump probability are unambiguous. A photon jump $\hat a$ changes cat parity as in (8.9)–(8.10).

<a id="g-097"></a>

#### G-097 — Normal ordering of a number exponential and Gibbs operator
**Derived from:** §§2.1, 2.5, pp. 14–15, 25–27, PDF 30–31, 41–43, Eqs. (2.35), (2.38)–(2.40), (2.132)–(2.145); Chapter 2 Problem 10, p. 41, PDF 57; ordering convention §3.7, p. 62, PDF 78, Eq. (3.107).
$$
e^{\lambda \hat n}={:\exp[(e^\lambda-1)\hat n]:},\quad
|0\rangle\langle0|={:e^{-\hat n}:},\quad
\hat\rho_{\rm th}=\frac1{1+\bar n}{:\exp[-\hat n/(1+\bar n)]:}.
$$
Coefficientwise proof: on $|m\rangle$, ${:\hat n^r:}$ has eigenvalue $m^{\underline r}$, so the ordered exponential is the finite binomial sum $\sum_{r=0}^m\binom mr(e^\lambda-1)^r=e^{\lambda m}$. These useful formulas are deductions, not equations explicitly printed in the book.

<a id="g-098"></a>

#### G-098 — Vacuum amplitude damping and its exact matrix-element sum
**Source:** §8.3, p. 201, PDF 217, Eq. (8.7); §8.5, p. 207, PDF 223, Eq. (8.26); **Exercise:** Chapter 8 Problem 6, p. 211, PDF 227, Eq. (8.39).
For $\hat H=0$, $\dot{\hat\rho}=\gamma(\hat a\hat\rho \hat a^\dagger-\{\hat n,\hat\rho\}/2)$. Put $T=e^{-\gamma t}$. The source exercise states
$$
\rho_{mn}(t)=T^{(m+n)/2}\sum_{\ell\ge0}\frac{(1-T)^\ell}{\ell!}
\sqrt{\frac{(m+\ell)!(n+\ell)!}{m!n!}}\,\rho_{m+\ell,n+\ell}(0).
$$
**Derived Kraus form:** $\hat\rho(t)=\sum_\ell \hat M_\ell\hat\rho(0)\hat M_\ell^\dagger$, $\hat M_\ell=\sqrt{(1-T)^\ell/\ell!}\,T^{\hat n/2}\hat a^\ell$. The order of the attenuation exponential and $\hat a^\ell$ is required. The binomial identity verifies $\sum \hat M_\ell^\dagger \hat M_\ell=\hat I$.

<a id="g-099"></a>

#### G-099 — Characteristic-function damping and coherent dyads
**Source:** §8.5, pp. 207–208, PDF 223–224, Eqs. (8.27)–(8.32); §10.5, p. 252, PDF 268, Eq. (10.43).
$$
C_N(\lambda,t)=C_N(\sqrt T\lambda,0),\quad
|\alpha\rangle\langle\alpha|\mapsto|\sqrt T\alpha\rangle\langle\sqrt T\alpha|.
$$
The source uses $\kappa$ once in (8.31) but the same damping constant is $\gamma$ in the surrounding equations. **Derived general dyad action:**
$$
\widehat\Lambda_T(|\alpha\rangle\langle\beta|)=
\exp\!\left[-\frac{1-T}{2}(|\alpha|^2+|\beta|^2-2\beta^*\alpha)\right]
|\sqrt T\alpha\rangle\langle\sqrt T\beta|.
$$
For opposite amplitudes the coherence factor is $e^{-2(1-T)|\alpha|^2}$, as in the source cat solution. Its short-time decay rate is $2\gamma|\alpha|^2$, while energy decays at $\gamma$.

<a id="g-100"></a>

#### G-100 — Driven damped coherent state and a printed factor-of-two error
**Source:** §8.6, p. 210, PDF 226, Eqs. (8.33)–(8.38), visually checked.
The stated equation is $\dot{\hat\rho}=-iG[\hat a+\hat a^\dagger,\hat\rho]+\gamma(\hat a\hat\rho \hat a^\dagger-\{\hat n,\hat\rho\}/2)$. **Derived correction:** for a coherent initial state,
$$
\dot\alpha=-iG-\gamma\alpha/2,\quad
\alpha(t)=\alpha_0e^{-\gamma t/2}-\frac{2iG}{\gamma}(1-e^{-\gamma t/2}),
\quad\alpha_{\rm ss}=-2iG/\gamma.
$$
Equivalently the master equation is $\gamma\widehat{\mathcal D}[\hat a+(2iG/\gamma)\hat I]\hat\rho$ with $\widehat{\mathcal D}[\hat b]\hat\rho=\hat b\hat\rho \hat b^\dagger-\{\hat b^\dagger \hat b,\hat\rho\}/2$. The printed (8.35) puts a minus before the dissipator and (8.36)–(8.37) use $-iG/\gamma$; those do not follow from (8.34). Short-time $\hat D(-iGt)$ in (8.38) is consistent with the corrected amplitude.

<a id="g-101"></a>

#### G-101 — Recoil exponential as a displacement; high-occupation restriction
**Source:** §10.7, pp. 256–258, PDF 272–274, Eqs. (10.61)–(10.73), especially p. 257, PDF 273 (visually checked).
The exact equal-time motional factor is
$$
e^{-i\eta(\hat a+\hat a^\dagger)}=\hat D(-i\eta)=e^{-\eta^2/2}e^{-i\eta \hat a^\dagger}e^{-i\eta \hat a}.
$$
**Derived free interaction-picture factor**, with the book's $\hat U_0=e^{-i\hat H_0t/\hbar}$:
$$
\hat U_0^\dagger e^{-i\eta(\hat a+\hat a^\dagger)}\hat U_0=
\exp[-i\eta(\hat ae^{-i\nu t}+\hat a^\dagger e^{i\nu t})].
$$
The source (10.66)–(10.68) reverses these oscillator phase signs relative to its own free Hamiltonian, thereby interchanging its red/blue sideband labels in (10.70)–(10.73). A first-order Lamb–Dicke expansion requires small matrix elements, e.g. $\eta\sqrt{2m+1}\ll1$ over significantly occupied Fock levels, not simply a small numerical $\eta$. Keeping higher Taylor powers is an approximation; the full displacement remains well defined.

<a id="g-102"></a>

#### G-102 — Dispersive Hamiltonian from a time-ordered oscillatory generator
**Source:** Appendix C, pp. 308–310, PDF 324–326, Eqs. (C1)–(C16); PDF 326 visually checked.
$$
\hat H_I(t)=\hbar g(\hat Ae^{i\Delta t}+\hat A^\dagger e^{-i\Delta t}),\quad
\hat U_I=\mathcal T\exp\left[-\frac i\hbar\int_0^t\hat H_I(t')dt'\right].
$$
The Dyson expansion has a bounded rapid first-order term of order $g/\Delta$ and a secular second-order term
$$
\hat U_I\approx \hat I-\frac{ig^2t}{\Delta}[\hat A,\hat A^\dagger],\quad
\hat H_{\rm eff}=\frac{\hbar g^2}{\Delta}[\hat A,\hat A^\dagger].
$$
This is perturbative elimination of fast dynamics, not an exact ordinary exponential identity. The smallness criterion is excitation dependent; control of both $\hat A$ and $\hat A^\dagger$ on the occupied states is needed. The source (C14) gives one corresponding expectation bound. Dropping micromotion and higher terms does not establish arbitrary-time accuracy.

<a id="g-103"></a>

#### G-103 — Off-resonant quantized parametric coupling generates Kerr terms
**Source:** Appendix C, p. 310, PDF 326, Eq. (C19), visually checked; underlying trilinear Hamiltonian §7.2, p. 165, PDF 181, Eq. (7.84).
For $\hat A=i\hat a^{\dagger2}\hat b$, $g=\chi^{(2)}$, $\Delta=2\omega-\omega_P$,
$$
[\hat A,\hat A^\dagger]=\hat n_a(\hat n_a-\hat I)-(4\hat n_a+2\hat I)\hat n_b,
$$
$$
\hat H_{\rm eff}=-\frac{\hbar[\chi^{(2)}]^2}{\Delta}
\{4\hat n_a\hat n_b+2\hat n_b-\hat n_a(\hat n_a-\hat I)\}.
$$
The explicit commutator evaluation is the deduction explaining the source result. Every term is a commuting number polynomial, so its effective propagator factors into self-Kerr, cross-Kerr, and pump rotation exponentials. The factorization is exact for the effective Hamiltonian, while the effective Hamiltonian remains approximate for the original interaction.

<a id="g-104"></a>

#### G-104 — Coherent trace resolution, reduced spectra, and density evolution
**Source:** §3.5, p. 53, PDF 69, Eq. (3.63); §3.7, p. 59, PDF 75, Eq. (3.88); Appendix A, pp. 295–303, PDF 311–319, Eqs. (A9)–(A22), (A33)–(A44), (A48)–(A52).
For trace-class $\hat F$,
$$
\operatorname{Tr}\hat F=\int\langle\alpha|\hat F|\alpha\rangle\frac{d^2\alpha}{\pi},
\quad\langle \hat O\rangle=\operatorname{Tr}(\hat\rho \hat O),\quad\hat\rho(t)=\hat U\hat\rho(0)\hat U^\dagger.
$$
For a normalized pure bipartite state in **Schmidt form** $\sum_jg_j|u_j\rangle|v_j\rangle$, with each local family orthonormal, both reduced operators have nonzero eigenvalues $|g_j|^2$, so $S=-\sum_j|g_j|^2\ln|g_j|^2$. This supplies the trace/entropy justification for G-074 and prevents confusing Schmidt amplitudes $g_j$ with probabilities $|g_j|^2$, as happens in part of the source Chapter 4 notation.

<a id="g-105"></a>

#### G-105 — Continuum CCR and pair-creation kernel
**Source:** Appendix D, pp. 312–313, PDF 328–329, Eqs. (D1)–(D7) and p. 312 footnote.
$$
[\hat a(\mathbf k),\hat a^\dagger(\mathbf k')]=\delta^{(3)}(\mathbf k-\mathbf k')\hat I.
$$
The first-order down-converted component is proportional to
$$
\int d^3k_s\,d^3k_i\,
\delta(\omega_P-\omega_s-\omega_i)\delta^{(3)}(\mathbf k_P-\mathbf k_s-\mathbf k_i)
 \hat a_s^\dagger(\mathbf k_s)\hat a_i^\dagger(\mathbf k_i)|0\rangle.
$$
The delta kernels describe the idealized phase-matching limit; finite interaction times/volumes broaden them. **Derived domain clarification:** for square-integrable mode functions, $\hat a(f)=\int f^*(k)\hat a(k)dk$ has $[\hat a(f),\hat a^\dagger(g)]=\langle f,g\rangle\hat I$. Pointwise $\hat a(k)$ is an operator-valued distribution, so arbitrary continuum exponential products require specified mode functions rather than treating delta functions as finite scalars.

### Source-error and uncertainty ledger

The following distinguishes actual visually confirmed printing from correction by algebra. Corrections in the entries are labelled as such rather than falsely attributed to the source. This is a topic-focused check, not a complete errata list for every equation of the book.

| Location | Verification status | Issue and action in this dossier |
|---|---|---|
| p. 49 / PDF 65, (3.32)–(3.33), and p. 40 / PDF 56, Problem 4 | Visual | The extracted text loses the inequality sign in $[\hat A,\hat B]\ne0$. The page itself is correct; central nested-commutator hypothesis retained. |
| p. 50 / PDF 66, intermediate line of (3.42), (3.43) | Visual | Intermediate (3.42) ends with an erroneous annihilation dagger; (3.43) omits the terminal vacuum ket in its first term. Final Weyl law and coherent displacement action used. |
| p. 65 / PDF 81, (3.128c), versus p. 66 / PDF 82, (3.132) | Visual | Printed “antinormal” characteristic function has creation before annihilation with incompatible signs; use (3.132), which also agrees with (3.129). |
| p. 56 / PDF 72, (3.81) | Full text; not separately image-checked | An unrestricted Hermitian operator need not have the asserted trace upper bound. Reconstruction used only with a convergent coherent generating expansion. |
| p. 100 / PDF 116, (4.145)–(4.146) | Original image rechecked in independent audit | Dressed-energy splitting omits $1/2$; the corrected value follows directly from the doublet Hamiltonian (G-039). |
| p. 106 / PDF 122, (4.186) | Visual | Excited number-state phase has plus sign, contradicting its Hamiltonian, (4.188), and (10.58). G-040 uses minus. |
| p. 124 / PDF 140, (5.53) | Original image rechecked in independent audit | Number states have first-order normalized coherence but zero mean fields, so the stated mean-field factorization cannot apply to them. G-036 separates the notions. |
| p. 132 / PDF 148, (5.100) | Original image rechecked in independent audit | Normalized $g^{(2)}$ of vacuum has zero denominator; it is undefined, rather than a well-defined zero. G-035–G-036 restrict denominators. |
| p. 148 / PDF 164, Chapter 6 Problem 7 | Visual | Even-only twin-Fock output is true at 50:50 mixing; it is not true at arbitrary $\hat J_1$-angle (identity mixing leaves $|N,N\rangle$). G-049 specifies 50:50. |
| p. 153 / PDF 169, prose after (7.18) | Visual | Axial uncertainty saturation is at $\theta=0$ or $\pi$, not the stated $\pi/2$. G-057–G-058 use the covariance criterion. |
| p. 159 / PDF 175, (7.46), (7.51) | Visual | Their stated eigenvalue labels conflict with (7.42)–(7.43). G-061 derives the rotated eigenvalue from $\gamma=\mu\alpha+\nu\alpha^*$. |
| p. 161 / PDF 177, prose and (7.63)–(7.64) | Visual | Prose prints $C_0=\sqrt{\cosh r}$, inconsistent with normalized (7.67)–(7.68); the summation also duplicates its zeroth contribution. G-056 uses $C_0=1/\sqrt{\cosh r}$. An independent odd zero-nullifier solution is excluded by the recurrence boundary equation. |
| p. 165 / PDF 181, (7.82)–(7.83) | Visual | General squeezed Q formula is inconsistent with the defined displaced state; Wigner Gaussian has wrong axis/scaling factors. G-063–G-064 are explicitly derived corrected formulas. |
| p. 174 / PDF 190, (7.113); p. 175 / PDF 191, (7.122); p. 176 / PDF 192, (7.123) | First two visual; (7.123) full text | Cat normalization lacks reciprocal square root; even/odd photon expressions have denominator/missing-factor errors. G-082–G-083 give normalized deductions. |
| p. 178 / PDF 194, (7.138) | Visual | Odd-cat Wigner denominator requires $1-e^{-2\alpha^2}$, not plus. Corrected unified expression in G-084. |
| p. 179 / PDF 195, (7.142) | Visual | Odd-number phase is $-i$, not $+i$. Final (7.143) is consistent with the correction (G-077). |
| p. 183 / PDF 199, (7.161) | Visual | $\sin r$ must be $\sinh r$; G-073 obtains it from the verified Bogoliubov transform. |
| p. 182 / PDF 198, prose after (7.158) | Visual | A correlated two-mode P distribution does not generally factor as $P(\alpha)P(\beta)$. No such factorization is assumed here. |
| p. 189 / PDF 205, odd-moment sentence | Original image rechecked in independent audit | Odd central moments can be negative. G-088 uses positive even noise moments only. |
| p. 192 / PDF 208, Problem 20 | Visual | $\hat H=-i\hbar \hat G$ requires $\hat U=e^{-t\hat G}$; the source prints the opposite exponent sign. G-092 provides both conventions explicitly. |
| pp. 200–202 / PDF 216–218, (8.3), (8.8), (8.11) | PDF218 original image rechecked; earlier pages by text | The normalized ket (8.8) needs the square root of its displayed expectation, and (8.11) omits the unnormalized attenuation scalar. G-096 states both correctly. |
| p. 210 / PDF 226, (8.35)–(8.37) | Visual | Dissipator sign and displacement factor disagree with (8.34). G-100 solves the stated equation and supplies $\alpha_{\rm ss}=-2iG/\gamma$. |
| p. 257 / PDF 273, (10.66)–(10.71) | Visual | Oscillator rotation phases reverse relative to the defined $\hat U_0$, interchanging sidebands. G-101 gives the consistent displacement conjugation and occupation-dependent expansion restriction. |
| p. 266 / PDF 282, (11.2)–(11.4) | Visual | At zero phase, the stated output is $i|1,1\rangle$, hence parity in either output is $-1$. Its (11.4) prints $+\cos2\theta$, while that output gives $-\cos2\theta$. Phase sensitivity is unaffected by this overall fringe sign. |
| p. 288 / PDF 304, (11.46)–(11.49) | Visual | Control-mode labels change from $d$ to $c$. G-086 relabels consistently. The real rotation in (11.51) alone also cannot supply its printed $\pm i$ factors. |

**Unverified limits and omissions.** No supplied text page remained unreadable, and there are no unread chapters or absent body pages. Full-text reading is complete, but equations outside the topic were not independently rederived or visually proofread exhaustively. The extra normal/Gauss squeeze factorizations (G-054–G-055, G-070–G-071), general number-function/Kerr moment identities (G-079–G-080), corrected Gaussian/cat expressions, normal number exponential, and Kraus form are algebraic deductions anchored to the indicated source generators, not claimed as explicit theorems of this introductory book. The source supplies no general convergence, self-adjointness, or disentangling theorem for arbitrary cubic or higher creation/annihilation polynomials, and no theorem forbidding all possible unitary preparations of an even/odd cat; statements near p. 179 about difficulty of unitary generation should not be read as such a universal impossibility claim. Higher-order *moments* in §7.8 are not higher-order unitary squeeze generators. A complete rigorous treatment of arbitrary unbounded-operator exponentials must come from a source addressing those analytic questions.

**Dossier totals:** 105 stable entries, G-001–G-105; complete full-text coverage of 333/333 PDF pages; 31 selected PDF pages visually inspected; no missing body pages, no unreadable text pages. All entries carry printed and PDF page locators and source equation/problem identifiers where present.


---

## Complete entry index


### X entries (9)

- [X-001 — Displacement: the basic linear exponential](#x-001)
- [X-002 — Single-mode squeezing: a quadratic exponential](#x-002)
- [X-003 — Two-mode bosonic squeezing](#x-003)
- [X-004 — Compact versus noncompact disentangling](#x-004)
- [X-005 — General quadratic exponentials reduce to finite matrix evolution](#x-005)
- [X-006 — Finite-mode bosonic and fermionic Gaussian vacua](#x-006)
- [X-007 — Number-conserving exponentials, normal order, and traces](#x-007)
- [X-008 — A practical derivation route and its limits](#x-008)
- [X-009 — Explicit multimode creator–number–annihilator factorization](#x-009)

### B entries (108)

- [B-001 — State spaces and occupation basis](#b-001)
- [B-002 — Involution, transpose, and regularized determinants](#b-002)
- [B-003 — Smeared CCR and CAR](#b-003)
- [B-004 — Vacuum reconstruction of arbitrary vectors](#b-004)
- [B-005 — Bounded fermionic fields](#b-005)
- [B-006 — Bosonic quadratures and their domains](#b-006)
- [B-007 — Number operator and a common field domain](#b-007)
- [B-008 — Commutant and irreducibility](#b-008)
- [B-009 — Uniqueness of the vacuum representation](#b-009)
- [B-010 — What normal form means analytically](#b-010)
- [B-011 — Bilinear forms are weaker than operators](#b-011)
- [B-012 — Bounded operators admit normal form](#b-012)
- [B-013 — Normal symbol of the vacuum projector](#b-013)
- [B-014 — Matrix kernel versus normal symbol](#b-014)
- [B-015 — Bosonic differential realization](#b-015)
- [B-016 — Fermionic derivative sides encode signs](#b-016)
- [B-017 — Coherent states under the name “Poisson vectors”](#b-017)
- [B-018 — Entire analyticity and bounded-set convergence](#b-018)
- [B-019 — Symbols as coherent matrix elements](#b-019)
- [B-020 — Rank-one projection symbols](#b-020)
- [B-021 — Continual integrals are projection-independent limits](#b-021)
- [B-022 — Real Gaussian integral](#b-022)
- [B-023 — Doubled complex Gaussian with sources](#b-023)
- [B-024 — Gaussian moments and bosonic contractions](#b-024)
- [B-025 — Analytic Hilbert-space realization](#b-025)
- [B-026 — Kernel action on a vector](#b-026)
- [B-027 — Product of bosonic operators](#b-027)
- [B-028 — Bosonic trace](#b-028)
- [B-029 — Grassmann parity and nilpotence](#b-029)
- [B-030 — Left/right derivatives and graded Leibniz rule](#b-030)
- [B-031 — Grassmann integration and partial-integration warning](#b-031)
- [B-032 — Inverse determinant Jacobian](#b-032)
- [B-033 — Fermionic Gaussian is a Pfaffian](#b-033)
- [B-034 — Clifford operators on a Grassmann algebra](#b-034)
- [B-035 — Infinite Grassmann functionals are not ordinary functions](#b-035)
- [B-036 — Doubled fermionic Gaussian](#b-036)
- [B-037 — Signed fermionic Gaussian moments](#b-037)
- [B-038 — Fermionic Fock inner product](#b-038)
- [B-039 — Fermionic kernel composition](#b-039)
- [B-040 — Fermionic trace differs from bosonic trace](#b-040)
- [B-041 — Normalized tracial state and infinite-mode factors](#b-041)
- [B-042 — Bosonic affine canonical transformations](#b-042)
- [B-043 — Bosonic block relations and inverse](#b-043)
- [B-044 — Symplectic geometry and group law](#b-044)
- [B-045 — Bosonic invertibility and squeezing contraction](#b-045)
- [B-046 — Bosonic implementability criterion](#b-046)
- [B-047 — Vacuum criterion for an implementer](#b-047)
- [B-048 — Bosonic transformed vacuum](#b-048)
- [B-049 — Full bosonic Gaussian kernel](#b-049)
- [B-050 — Displacement operator in the source’s conjugation direction](#b-050)
- [B-051 — One-mode squeezing as a checkable specialization](#b-051)
- [B-052 — Projective phases and the double cover](#b-052)
- [B-053 — Passive transformation of normal symbols](#b-053)
- [B-054 — General bosonic transformation kernel for normal symbols](#b-054)
- [B-055 — Gaussian normal symbols are closed under canonical conjugation](#b-055)
- [B-056 — Fermionic canonical relations](#b-056)
- [B-057 — Fermionic boundedness and orthogonal geometry](#b-057)
- [B-058 — Fermionic implementability, including property A](#b-058)
- [B-059 — Singular fermionic charts and filled sectors](#b-059)
- [B-060 — Fermionic Gaussian implementer kernel](#b-060)
- [B-061 — Fermionic paired vacuum](#b-061)
- [B-062 — Fermionic projective/spin representation](#b-062)
- [B-063 — Fermionic transformation of normal symbols](#b-063)
- [B-064 — Fermionic Gaussian-symbol conjugation and determinant class](#b-064)
- [B-065 — General quadratic expression](#b-065)
- [B-066 — Vacuum-domain necessity](#b-066)
- [B-067 — Self-adjointness of quadratic normal forms](#b-067)
- [B-068 — Exponentials of bounded perturbations](#b-068)
- [B-069 — Integrated pairing criterion for implementable evolution](#b-069)
- [B-070 — Determinant-class phase condition](#b-070)
- [B-071 — Bosonic quadratic propagator, including the scalar phase](#b-071)
- [B-072 — Heisenberg/Hadamard action closes at first degree](#b-072)
- [B-073 — Log-determinant differentiation for the propagator](#b-073)
- [B-074 — Fermionic quadratic propagator](#b-074)
- [B-075 — Resolvent formula requires a convergence qualification](#b-075)
- [B-076 — Necessary condition for a singular quadratic Hamiltonian](#b-076)
- [B-077 — Sufficient condition and a visible missing identity term](#b-077)
- [B-078 — Sufficient reduction to quadratic normal form](#b-078)
- [B-079 — Weaker conditions and infinite additive constants](#b-079)
- [B-080 — Continuum translation-invariant pairing can fail implementability](#b-080)
- [B-081 — Forced oscillator field and the sign in completing the square](#b-081)
- [B-082 — Weighted coupling regimes and domain dependence](#b-082)
- [B-083 — Implementable displacement with divergent energy shift](#b-083)
- [B-084 — Bosonic canonical diagonalization theorem](#b-084)
- [B-085 — Siegel-domain action and Riccati fixed point](#b-085)
- [B-086 — Fermionic canonical diagonalization theorem](#b-086)
- [B-087 — Stability of square roots in Schatten ideals](#b-087)
- [B-088 — Inequivalent CAR representations in the continuum](#b-088)
- [B-089 — A formal CAR rearrangement can prove the false identity $k=0$](#b-089)
- [B-090 — Operator approximation before taking singular limits](#b-090)
- [B-091 — Well-defined evolved fields without a Fock-space Hamiltonian](#b-091)
- [B-092 — Scattering and renormalized fields are separate questions](#b-092)
- [B-093 — Two-body interactions as quartic normal forms](#b-093)
- [B-094 — Dyson expansion and time ordering](#b-094)
- [B-095 — Free-field time dependence and contraction kernel](#b-095)
- [B-096 — Graded time ordering and equal-time convention](#b-096)
- [B-097 — Retarded contraction and the two-field formula](#b-097)
- [B-098 — Wick’s theorem as a sum over contractions](#b-098)
- [B-099 — Functional differential form of Wick’s theorem](#b-099)
- [B-100 — Normal symbol of interaction-picture evolution](#b-100)
- [B-101 — Contractions acting on an exponential source](#b-101)
- [B-102 — Determinant identity for a quadratic symbol](#b-102)
- [B-103 — Quartic interactions through an auxiliary commuting Gaussian](#b-103)
- [B-104 — Retarded loops vanish, so the auxiliary determinant is one](#b-104)
- [B-105 — Auxiliary-field propagator representation](#b-105)
- [B-106 — Grand partition functions: determinant versus reciprocal determinant](#b-106)
- [B-107 — Cauchy Green function and free thermal traces](#b-107)
- [B-108 — Coverage limits: BCH, Weyl/anti-normal order, and what is absent](#b-108)

### P entries (107)

- [P-001 — CCR, canonical variables, and the Heisenberg–Weyl Lie algebra](#p-001)
- [P-002 — Central-commutator BCH and a differential proof](#p-002)
- [P-003 — Weyl composition, group commutator, and geometric phase](#p-003)
- [P-004 — Stone–von Neumann equivalence and its domain qualification](#p-004)
- [P-005 — Fock basis and exponentials applied to the vacuum](#p-005)
- [P-006 — Position-space factorization of displacement](#p-006)
- [P-007 — General displaced fiducials, overlaps, and resolution of identity](#p-007)
- [P-008 — Displacement conjugation and invariance of centered variances](#p-008)
- [P-009 — Minimum uncertainty and a new annihilation operator](#p-009)
- [P-010 — Normal and antinormal displacement factorizations](#p-010)
- [P-011 — Differential correspondences for multiplying displacement](#p-011)
- [P-012 — Coherent-state exponentials, overlaps, and measures](#p-012)
- [P-013 — Bargmann realization of operator exponentials](#p-013)
- [P-014 — Weyl-lattice commutation and completeness theorems](#p-014)
- [P-015 — Covariant, contravariant, normal, and antinormal symbols](#p-015)
- [P-016 — Trace inequalities for exponential operators](#p-016)
- [P-017 — Distributional trace of displacement and operator reconstruction](#p-017)
- [P-018 — Displacement matrix elements and Laguerre polynomials](#p-018)
- [P-019 — Weyl symbols, displaced parity, and a printed notation issue](#p-019)
- [P-020 — Ordered characteristic functions and moments](#p-020)
- [P-021 — Group-orbit coherent states and why disentangling works](#p-021)
- [P-022 — Maximal isotropy and quadratic/oscillator algebras](#p-022)
- [P-023 — Gaussian decomposition is local, with a singular chart boundary](#p-023)
- [P-024 — $SU(2)$ commutators, ladder action, and differential realization](#p-024)
- [P-025 — $SU(2)$ exponential disentangling](#p-025)
- [P-026 — $SU(2)$ coherent state, overlap, annihilator, and identity resolution](#p-026)
- [P-027 — Rotation exponential kernels and oscillator contraction](#p-027)
- [P-028 — $SU(1,1)$ algebra and positive discrete series](#p-028)
- [P-029 — $SU(1,1)$ disentangling, with a verified source sign correction](#p-029)
- [P-030 — $SU(1,1)$ coherent states, overlaps, annihilator, and measure](#p-030)
- [P-031 — Single-mode squeezing as $SU(1,1)$](#p-031)
- [P-032 — Two-mode squeezing and conserved number difference](#p-032)
- [P-033 — Principal-series and Lorentz coherent systems: applicability boundary](#p-033)
- [P-034 — General finite-mode bosonic Bogoliubov conditions](#p-034)
- [P-035 — Quadratic bosonic closure, generators, and metaplectic phase](#p-035)
- [P-036 — General bosonic triangular factorization and vacuum chart](#p-036)
- [P-037 — Bosonic Gaussian normalization, overlaps, and normalizability theorem](#p-037)
- [P-038 — Fractional-linear Gaussian parameter action](#p-038)
- [P-039 — Pure bosonic Gaussian bilinear symbols and constraints](#p-039)
- [P-040 — Bosonic Gaussian geometry and representation parameter restriction](#p-040)
- [P-041 — Fermionic CAR and finite-mode Bogoliubov constraints](#p-041)
- [P-042 — Fermionic quadratic closure, parity, and the spin cover](#p-042)
- [P-043 — Fermionic pair-exponential state, normalization, and overlap](#p-043)
- [P-044 — Fermionic analytic group action](#p-044)
- [P-045 — Fermionic bilinear symbols and purity constraints](#p-045)
- [P-046 — Nilpotent groups, exponential coordinates, and square integrability](#p-046)
- [P-047 — Root-space coherent exponentials and general Gaussian decomposition](#p-047)
- [P-048 — Holomorphic discrete series, determinant kernels, and admissible measures](#p-048)
- [P-049 — Principal-series kernels: locator record for the remaining first-part pages](#p-049)
- [P2-001 — Principal-series coherent systems need not be square integrable](#p2-001)
- [P2-002 — Symplectic tube domains](#p2-002)
- [P2-003 — Horospherical kernels as curved plane waves](#p2-003)
- [P2-004 — Completeness of discrete coherent subsystems](#p2-004)
- [P2-005 — Covariant symbols, analytic continuation, and product kernel](#p2-005)
- [P2-006 — Curved symbol products and the correspondence principle](#p2-006)
- [P2-007 — Curved analogues of Weyl quantization and compact quantization](#p2-007)
- [P2-008 — General Kähler quantization](#p2-008)
- [P2-009 — Linear Lie-algebra Hamiltonians preserve coherent orbits](#p2-009)
- [P2-010 — Forced oscillator: displacement and accumulated phase](#p2-010)
- [P2-011 — Forced-oscillator transitions are displacement matrix elements](#p2-011)
- [P2-012 — Parametric oscillator realizes two $SU(1,1)$ sectors](#p2-012)
- [P2-013 — Exact Riccati dynamics and phase for squeezing](#p2-013)
- [P2-014 — Floquet spectral types of a quadratic oscillator](#p2-014)
- [P2-015 — Singular oscillator and operator domains](#p2-015)
- [P2-016 — Many-body inverse-square models inherit the same algebra](#p2-016)
- [P2-017 — Visible source error in the affine Gaussian normalization](#p2-017)
- [P2-018 — Gaussian generating kernel for a forced parametric oscillator](#p2-018)
- [P2-019 — Bosonic pair creation: two-mode $\mathfrak{su}(1,1)$](#p2-019)
- [P2-020 — Multimode bosonic canonical blocks](#p2-020)
- [P2-021 — Multimode paired vacuum and its determinant normalization](#p2-021)
- [P2-022 — Bosonic overlaps and fractional-linear action](#p2-022)
- [P2-023 — Bosonic pair-number generating determinant](#p2-023)
- [P2-024 — Multimode Riccati equation and vacuum phase](#p2-024)
- [P2-025 — Single fermionic pair is an $SU(2)$ problem](#p2-025)
- [P2-026 — Multimode fermionic canonical blocks](#p2-026)
- [P2-027 — Fermionic paired vacuum](#p2-027)
- [P2-028 — Fermionic overlaps and Grassmannian action](#p2-028)
- [P2-029 — Fermionic pair-number generating determinant](#p2-029)
- [P2-030 — Representation multiplicities and finite/infinite restrictions](#p2-030)
- [P2-031 — Clebsch–Gordan coefficients from polynomial coherent states](#p2-031)
- [P2-032 — Quasiclassical partition-function bounds](#p2-032)
- [P2-033 — Bilinear fermion algebras and large-$N$ mechanics](#p2-033)
- [P2-034 — Thermal oscillator master equation](#p2-034)
- [P2-035 — Normal, Weyl, and anti-normal characteristic functions](#p2-035)
- [P2-036 — Exact attenuation solution for characteristic functions](#p2-036)
- [P2-037 — Quasiprobability and Fokker–Planck forms](#p2-037)
- [P2-038 — Visible missing minus sign in the Gaussian convolution](#p2-038)
- [P2-039 — Gaussian states stay Gaussian under relaxation](#p2-039)
- [P2-040 — Analytic density kernel generates all Fock matrix elements](#p2-040)
- [P2-041 — Gaussian Green kernel for driven relaxation](#p2-041)
- [P2-042 — Attenuation, additive noise, and displacement factorization](#p2-042)
- [P2-043 — Displaced thermal density kernel and pure-state limit](#p2-043)
- [P2-044 — Populations of a displaced thermal state](#p2-044)
- [P2-045 — Number-state input and spin relaxation](#p2-045)
- [P2-046 — Free oscillator propagator and thermal trace](#p2-046)
- [P2-047 — Landau-level and determinant applications](#p2-047)
- [P2-048 — Quadratic trace in the Heisenberg–Euler application is analytically continued](#p2-048)
- [P2-049 — Synchrotron and entropy applications](#p2-049)
- [P2-050 — Entire-function completeness and the critical coherent lattice](#p2-050)
- [P2-051 — Displacement kernel and finite coefficient sum](#p2-051)
- [P2-052 — Displacement matrix elements in associated Laguerre form](#p2-052)
- [P2-053 — Anti-normal displacement factorization and Gaussian integral](#p2-053)
- [P2-054 — Visible error in the generating function of displacement probabilities](#p2-054)
- [P2-055 — Semiclassical geometry and derivative formula for displacement](#p2-055)
- [P2-056 — Jacobians of classical-domain transformations](#p2-056)
- [P2-057 — Addendum, references, and historical scope](#p2-057)
- [P2-058 — Source-audit cautions for reuse](#p2-058)

### C entries (99)

- [C-001 — CCR and the representation domain](#c-001)
- [C-002 — Central-commutator BCH with operator hypotheses](#c-002)
- [C-003 — Weyl multiplication and central extension](#c-003)
- [C-004 — Hadamard conjugation and coordinate factorization](#c-004)
- [C-005 — Complex labels and a usable analytic core](#c-005)
- [C-006 — Normal-ordered displacement](#c-006)
- [C-007 — Vacuum action, eigenvalue, and number expansion](#c-007)
- [C-008 — Displacement expectation and coherent overlap](#c-008)
- [C-009 — Harmonic rotation of exponentials](#c-009)
- [C-010 — Resolution, trace, and Hilbert–Schmidt formulas](#c-010)
- [C-011 — Bargmann realization](#c-011)
- [C-012 — Mehler kernel and analytic continuation](#c-012)
- [C-013 — Weyl symbol and translation expansion](#c-013)
- [C-014 — Symbol traces and Gaussian overlaps](#c-014)
- [C-015 — Exact symplectic covariance](#c-015)
- [C-016 — Moyal product and exact quadratic commutator](#c-016)
- [C-017 — Wick, anti-Wick, and smoothing](#c-017)
- [C-018 — Quadratic evolution is a finite matrix problem](#c-018)
- [C-019 — Exact Gaussian propagation](#c-019)
- [C-020 — Siegel-domain protection and determinant cocycle](#c-020)
- [C-021 — Riccati and amplitude equations](#c-021)
- [C-022 — Coherent matrix elements without coordinate caustics](#c-022)
- [C-023 — Metaplectic double cover and parity](#c-023)
- [C-024 — Bogoliubov blocks](#c-024)
- [C-025 — Squeeze generator and bounded squeeze coordinate](#c-025)
- [C-026 — Squeeze conjugation](#c-026)
- [C-027 — Squeezed vacuum from an operator Riccati equation](#c-027)
- [C-028 — Gaussian coordinate form of squeezed vacuum](#c-028)
- [C-029 — Squeeze–rotation factorization of a quadratic propagator](#c-029)
- [C-030 — Weyl symbol through a Cayley transform](#c-030)
- [C-031 — Covariant translation symbol and phase tracking](#c-031)
- [C-032 — Driven/affine quadratic propagators](#c-032)
- [C-033 — What quadratic replacement does and does not prove](#c-033)
- [C-034 — Propagator reconstruction and Gaussian determinant branches](#c-034)
- [C-035 — Periodic quadratic evolution and quasi-energy](#c-035)
- [C-036 — Trace formulas require a spectral cutoff](#c-036)
- [C-037 — Torus Weyl exponentials and finite matrices](#c-037)
- [C-038 — SU(2) commutators, ladder matrices, and Casimir](#c-038)
- [C-039 — Fundamental 2-by-2 exponential](#c-039)
- [C-040 — Gauss decomposition as a disentangling theorem](#c-040)
- [C-041 — SU(2) normal factorization, corrected explicit parameters](#c-041)
- [C-042 — SU(2) products, overlaps, and trace resolution](#c-042)
- [C-043 — Spin symbols and exponential trace inequalities](#c-043)
- [C-044 — Spin-to-CCR contraction](#c-044)
- [C-045 — SU(1,1) algebra and analytic representation domains](#c-045)
- [C-046 — Bosonic realizations of SU(1,1)](#c-046)
- [C-047 — SU(1,1) coherent series and overlap](#c-047)
- [C-048 — Noncompact Gauss factorization and squeeze disentangling](#c-048)
- [C-049 — One-mode squeezing and Wigner ellipse](#c-049)
- [C-050 — Time-dependent SU(1,1) Hamiltonian](#c-050)
- [C-051 — Derivative of a noncommuting exponential](#c-051)
- [C-052 — Singular quadratic-looking generators require a boundary condition](#c-052)
- [C-053 — A consistent oscillator Weyl exponential](#c-053)
- [C-054 — Exact Herman–Kluk representation for quadratic exponentials](#c-054)
- [C-055 — Gaussian density operators and thermal factors](#c-055)
- [C-056 — Partial trace of a quadratic exponential evolution](#c-056)
- [C-057 — Segal fields and a rigorous exponential domain](#c-057)
- [C-058 — Second quantization commutes with one-particle exponentiation](#c-058)
- [C-059 — Infinite-mode displacement, exponential vector, and Poisson law](#c-059)
- [C-060 — Weyl shifts preserve number-domain scales](#c-060)
- [C-061 — Derivative of a moving infinite-mode displacement](#c-061)
- [C-062 — Quadratic fluctuations about an interacting boson field](#c-062)
- [C-063 — Dyson expansion versus convergent exponential factorization](#c-063)
- [C-064 — CAR, finite dimension, and bounded generators](#c-064)
- [C-065 — Ordinary one-mode fermionic linear exponential](#c-065)
- [C-066 — Nilpotent functions, parity, and derivative orientation](#c-066)
- [C-067 — Berezin integration is ordered coefficient extraction](#c-067)
- [C-068 — Fermionic Gaussian determinant and sources](#c-068)
- [C-069 — Pfaffian Gaussian integral](#c-069)
- [C-070 — Fermionic Bargmann model and normal monomials](#c-070)
- [C-071 — Grassmann Weyl translation and central BCH](#c-071)
- [C-072 — Fermionic coherent state, overlap, and identity resolution](#c-072)
- [C-073 — Trace, supertrace, and Berezinian](#c-073)
- [C-074 — Weyl translation expansion and its supertrace inversion](#c-074)
- [C-075 — Normal, antinormal, and Weyl characteristic functions](#c-075)
- [C-076 — Ordinary trace of fermionic Weyl-symbol products](#c-076)
- [C-077 — Fermionic star product and star exponential](#c-077)
- [C-078 — Fermi oscillator exponential and removable Cayley singularity](#c-078)
- [C-079 — Fermionic partition determinant](#c-079)
- [C-080 — General quadratic fermion Hamiltonian and BdG matrix](#c-080)
- [C-081 — Spin implementation of fermionic canonical transformations](#c-081)
- [C-082 — Fermionic squeezed vacuum and a true chart boundary](#c-082)
- [C-083 — Fermionic Riccati and determinant amplitude](#c-083)
- [C-084 — Fermionic Mehler kernel and coherent matrix elements](#c-084)
- [C-085 — Fermionic Cayley symbol of a quadratic exponential](#c-085)
- [C-086 — Covariant fermionic Cayley symbol and Pfaffian branch](#c-086)
- [C-087 — Fermion pair rotations are compact disentanglings](#c-087)
- [C-088 — Mixed boson–fermion displacement](#c-088)
- [C-089 — Super star product](#c-089)
- [C-090 — Supercommutator and nilpotent continuation](#c-090)
- [C-091 — Commuting Gaussian integral and branch](#c-091)
- [C-092 — Stationary phase is asymptotic, with a determinant on normal directions](#c-092)
- [C-093 — Controlled Moyal expansion and functional calculus](#c-093)
- [C-094 — Adjoint exponential and local Lie coordinates](#c-094)
- [C-095 — Schur uniqueness, square integrability, and group coherent states](#c-095)
- [C-096 — Cartan polar decomposition](#c-096)
- [C-097 — Berezin symbol product on the disk](#c-097)
- [C-098 — Grassmann Fourier involution with an explicit measure convention](#c-098)
- [C-099 — Holomorphic Grassmann reproducing kernel](#c-099)

### G entries (105)

- [G-001 — Canonical oscillator realization](#g-001)
- [G-002 — Hadamard conjugation lemma](#g-002)
- [G-003 — Central-commutator BCH/disentangling](#g-003)
- [G-004 — Number commutators and free rotations](#g-004)
- [G-005 — Fock basis and ladder matrix elements](#g-005)
- [G-006 — Number-shift rules and factorial powers](#g-006)
- [G-007 — Quadrature scaling and vacuum covariance](#g-007)
- [G-008 — Independent-mode CCR and tensor-product Fock states](#g-008)
- [G-009 — Thermal exponential, partition function, and spectrum](#g-009)
- [G-010 — Thermal factorial moments](#g-010)
- [G-011 — Obstruction to a globally canonical Hermitian phase](#g-011)
- [G-012 — Susskind–Glogower one-sided shift](#g-012)
- [G-013 — Coherent-state eigenvector and expansion](#g-013)
- [G-014 — No normalizable right eigenvector of creation](#g-014)
- [G-015 — Displacement exponential and its two ordered forms](#g-015)
- [G-016 — Vacuum action of the ordered displacement](#g-016)
- [G-017 — Weyl composition and noncommuting phases](#g-017)
- [G-018 — Translation of operators and arbitrary ordered functions](#g-018)
- [G-019 — Coherent overlap and reproducing measure](#g-019)
- [G-020 — Normally ordered coherent matrix elements](#g-020)
- [G-021 — Coherent statistics and factorial moments](#g-021)
- [G-022 — Linearly driven evolution is a displacement up to phase](#g-022)
- [G-023 — Segal–Bargmann entire generating function](#g-023)
- [G-024 — Reconstruction from coherent diagonal symbols](#g-024)
- [G-025 — Coherent-projector differential correspondences](#g-025)
- [G-026 — Normal-order notation and optical equivalence](#g-026)
- [G-027 — Singular P representations must be handled weakly](#g-027)
- [G-028 — Three ordered characteristic functions](#g-028)
- [G-029 — Ordered moments from derivatives](#g-029)
- [G-030 — Cahill–Glauber s ordering](#g-030)
- [G-031 — Complex Fourier kernels and phase-space distributions](#g-031)
- [G-032 — Thermal Gaussian characteristic functions](#g-032)
- [G-033 — Coherent and Fock Q/Wigner functions](#g-033)
- [G-034 — Displaced parity representation of the Wigner function](#g-034)
- [G-035 — Photoabsorption fixes normal ordering](#g-035)
- [G-036 — Factorial correlation and number variance](#g-036)
- [G-037 — Atomic transition algebra coupled to a boson](#g-037)
- [G-038 — Exact resonant operator-valued exponential](#g-038)
- [G-039 — Invariant doublet and detuned exponential](#g-039)
- [G-040 — Dispersive interaction as conditional number rotation](#g-040)
- [G-041 — Raman number-controlled atomic rotation](#g-041)
- [G-042 — Multiphoton Jaynes–Cummings blocks](#g-042)
- [G-043 — Two-port unitarity is required to preserve CCR](#g-043)
- [G-044 — Beam-splitter exponential and its adjoint action](#g-044)
- [G-045 — Real beam-splitter rotation](#g-045)
- [G-046 — Schwinger SU(2) algebra and fixed-number sectors](#g-046)
- [G-047 — Binomial state from a number-state input](#g-047)
- [G-048 — Coherent products remain products under passive mixing](#g-048)
- [G-049 — Twin-Fock interference and even outputs](#g-049)
- [G-050 — Phase shifter and parity exponential](#g-050)
- [G-051 — Single-mode quadratic squeeze exponential](#g-051)
- [G-052 — Single-mode Bogoliubov transform](#g-052)
- [G-053 — Quadratic SU(1,1) closure](#g-053)
- [G-054 — Single-mode normal Gauss factorization](#g-054)
- [G-055 — Single-mode reversed Gauss factorization](#g-055)
- [G-056 — Squeezed-vacuum nullifier and unique normalized solution](#g-056)
- [G-057 — General quadrature variance and squeezing axis](#g-057)
- [G-058 — Covariance matters for rotated squeezed states](#g-058)
- [G-059 — Displaced squeezed state and low moments](#g-059)
- [G-060 — Interchanging displacement and squeezing](#g-060)
- [G-061 — Bogoliubov eigenvalue definition of the displaced squeezed state](#g-061)
- [G-062 — Hermite coefficients of a displaced squeezed vacuum](#g-062)
- [G-063 — Coherent overlap with squeezed vacuum and squeezed Q function](#g-063)
- [G-064 — Squeezed Gaussian Wigner function with consistent scaling](#g-064)
- [G-065 — Number variance of a displaced squeezed vacuum](#g-065)
- [G-066 — Degenerate parametric Hamiltonian and squeeze evolution](#g-066)
- [G-067 — Four-wave mixing gives the same quadratic exponent after approximation](#g-067)
- [G-068 — Two-mode squeeze exponential and Bogoliubov action](#g-068)
- [G-069 — Two-mode SU(1,1) and number-difference invariant](#g-069)
- [G-070 — Two-mode normal Gauss factorization](#g-070)
- [G-071 — Two-mode reversed factorization and supermode decomposition](#g-071)
- [G-072 — Twin-beam vacuum expansion and nullifiers](#g-072)
- [G-073 — Collective squeezed quadrature and thermal individual modes](#g-073)
- [G-074 — Reduced Gibbs state and entropy of a twin beam](#g-074)
- [G-075 — Nondegenerate parametric evolution and its parameter](#g-075)
- [G-076 — Self-Kerr number-squared exponential](#g-076)
- [G-077 — Exact half-revival cat identity](#g-077)
- [G-078 — Normally ordered Kerr generator differs by a rotation](#g-078)
- [G-079 — General number-function conjugation and Kerr Heisenberg solution](#g-079)
- [G-080 — Kerr moments and unchanged photon statistics](#g-080)
- [G-081 — Cross-Kerr conditional rotation](#g-081)
- [G-082 — General coherent-cat normalization and eigenoperator](#g-082)
- [G-083 — Cat parity projectors and photon probabilities](#g-083)
- [G-084 — Cat Wigner interference with correct odd normalization](#g-084)
- [G-085 — Conditional projectors produce cats but are not unitary exponentials](#g-085)
- [G-086 — Cross-Kerr conjugation into a conditional SU(2) rotation](#g-086)
- [G-087 — Linear quadrature exponential generates normal-order contractions](#g-087)
- [G-088 — All centered quadrature moments in terms of normal moments](#g-088)
- [G-089 — Squared-field quadratures and their noncentral commutator](#g-089)
- [G-090 — Pair coherent states with fixed number difference](#g-090)
- [G-091 — Normally ordered intermode Cauchy–Schwarz test](#g-091)
- [G-092 — Frequency conversion is an SU(2) exponential, with a source sign caveat](#g-092)
- [G-093 — Photon subtraction and addition are conditional operations](#g-093)
- [G-094 — Broadband collective CCR and finite bandwidth](#g-094)
- [G-095 — Low-gain down-conversion truncation and multipair warning](#g-095)
- [G-096 — Nonunitary number exponential and coherent-state normalization](#g-096)
- [G-097 — Normal ordering of a number exponential and Gibbs operator](#g-097)
- [G-098 — Vacuum amplitude damping and its exact matrix-element sum](#g-098)
- [G-099 — Characteristic-function damping and coherent dyads](#g-099)
- [G-100 — Driven damped coherent state and a printed factor-of-two error](#g-100)
- [G-101 — Recoil exponential as a displacement; high-occupation restriction](#g-101)
- [G-102 — Dispersive Hamiltonian from a time-ordered oscillatory generator](#g-102)
- [G-103 — Off-resonant quantized parametric coupling generates Kerr terms](#g-103)
- [G-104 — Coherent trace resolution, reduced spectra, and density evolution](#g-104)
- [G-105 — Continuum CCR and pair-creation kernel](#g-105)


## Rechecking record

The book-specific coverage and error ledgers above are the authoritative record of which displays were checked visually and which expressions remain formal, restricted, or uncertain.
