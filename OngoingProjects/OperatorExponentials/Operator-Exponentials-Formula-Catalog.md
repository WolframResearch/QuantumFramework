# Operator Exponentials: Formula Catalog

This catalog lists every formula in the bosonic and fermionic operator-exponential atlas [Bosonic-Fermionic-Operator-Exponentials.md](Bosonic-Fermionic-Operator-Exponentials.md), each with the atlas IDs where it appears, one plain sentence saying what it states, the formula, and its conditions. The atlas's 428 IDs map to 326 distinct formulas, grouped by topic below, plus 74 records that carry no formula; when several books state the same formula, the first version is given in full and each other book's version follows under "Also stated in".

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
- [Records with no formula](#records-with-no-formula)

## BCH and Hadamard

### Hadamard series for a closed linear commutator action
Appears in: X-005 (source records cited: P-034, P-035, P-041, P-042, B-071 to B-074, C-018 to C-024)

When commuting with $\hat Q$ maps a list of operators linearly into itself, conjugation by $e^{t\hat Q}$ is the matrix exponential of that linear map.

$$ [\hat Q,\hat R_i]=\sum_jL_{ij}\hat R_j\ \Longrightarrow\ e^{t\hat Q}\hat R e^{-t\hat Q}=e^{tL}\hat R $$

Conditions: $L$ scalar; time-dependent noncommuting generators need time ordering.

### Central-commutator BCH
Appears in: P-002 (also used by P-010, P2-035); C-002 (also used in C-006, C-059, C-071); G-003

When the commutator of two operators commutes with both, the product of their exponentials equals the exponential of the sum times the exponential of half the commutator.

$$ e^{\hat A}e^{\hat B}=e^{[\hat A,\hat B]/2}e^{\hat A+\hat B},\qquad e^{\hat A+\hat B}=e^{-[\hat A,\hat B]/2}e^{\hat A}e^{\hat B} $$

Conditions: $[\hat A,[\hat A,\hat B]]=[\hat B,[\hat A,\hat B]]=0$; common domain for unbounded generators.

**Also stated in C-002 (Central-commutator BCH):**

When the commutator of two generators commutes with both, the exponential of their sum is the product of their exponentials times the scalar exponential of minus half the commutator.

$$ e^{\hat A+\hat B}=e^{-[\hat A,\hat B]/2}e^{\hat A}e^{\hat B} $$

Conditions: [A,[A,B]] = [B,[A,B]] = 0 on an invariant common core; anti-self-adjoint generators.

**Also stated in G-003 (Central-commutator BCH):**

When the commutator of A and B commutes with both, the exponential of the sum splits into the product of exponentials times a scalar exponential of half the commutator, with the sign set by the order.

$$ e^{\hat A+\hat B}=e^{\hat A}e^{\hat B}e^{-[\hat A,\hat B]/2}=e^{\hat B}e^{\hat A}e^{[\hat A,\hat B]/2} $$

Conditions: $[\hat A,[\hat A,\hat B]]=[\hat B,[\hat A,\hat B]]=0$.

### Duhamel derivative of an exponential
Appears in: C-051

The derivative of the exponential of a moving generator is an average of the derivative conjugated along the exponential, not simply L' e^L.

$$ \frac d{dt}e^{\hat L(t)}=\int_0^1e^{s\hat L(t)}\dot{\hat L}(t)e^{(1-s)\hat L(t)}ds $$

### Adjoint exponential
Appears in: C-094

Conjugation by an exponential is the exponential of the adjoint map, which computes the conjugation whenever the commutator span closes.

$$ \operatorname{Ad}(e^X)=e^{\operatorname{ad}X},\qquad \operatorname{ad}X(Y)=[X,Y] $$

### Hadamard conjugation series
Appears in: G-002 (applied in G-018, G-044, G-052, G-068)

Conjugating B by an exponential of A sums the nested commutators of A with B; the series closes in finitely many terms when the nested commutators vanish or repeat.

$$ e^{s\hat A}\hat Be^{-s\hat A}=\sum_{k\ge0}\frac{s^k}{k!}\widehat{\operatorname{ad}}_{\hat A}^k(\hat B),\quad\widehat{\operatorname{ad}}_{\hat A}(\hat B)=[\hat A,\hat B]. $$

## Displacement

### Displacement: normal and antinormal factorization
Appears in: X-001 (source records cited: P-002, P-003, P-010, G-015, G-017, B-050, C-006); P-010, P2-053; C-006; same identity in C-059 (infinite-mode form); G-015 (also G-101 at $\alpha=-i\eta$)

The displacement operator splits into a creation exponential times an annihilation exponential (or the reverse), and the only cost of reordering is a Gaussian scalar.

$$ \hat D(\alpha)=e^{\alpha \hat a^\dagger-\bar\alpha \hat a}=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}e^{-\bar\alpha \hat a}=e^{|\alpha|^2/2}e^{-\bar\alpha \hat a}e^{\alpha \hat a^\dagger},\qquad [\hat a,\hat a^\dagger]=\hat I $$

**Also stated in P-010 (Displacement: normal and antinormal factorizations):**

The displacement in normal order carries the scalar $e^{-|\alpha|^2/2}$ and in antinormal order carries $e^{+|\alpha|^2/2}$.

$$ \hat D(\alpha)=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}e^{-\bar\alpha \hat a}=e^{+|\alpha|^2/2}e^{-\bar\alpha \hat a}e^{\alpha \hat a^\dagger} $$

**Also stated in C-006 (Normal and antinormal displacement factorization, displacement conjugation):**

The displacement operator equals its normal-ordered product times e^{-|alpha|^2/2}, or its antinormal product times e^{+|alpha|^2/2}, and it shifts the annihilator by -alpha.

$$ \hat D(\alpha)=e^{\alpha\cdot \hat a^\dagger-\bar\alpha\cdot \hat a}=e^{-|\alpha|^2/2}e^{\alpha\cdot \hat a^\dagger}e^{-\bar\alpha\cdot \hat a}=e^{|\alpha|^2/2}e^{-\bar\alpha\cdot\hat a}e^{\alpha\cdot\hat a^\dagger},\quad \hat D(\alpha)\hat a\hat D(\alpha)^{-1}=\hat a-\alpha\hat I $$

**Also stated in G-015 (Displacement: normal and antinormal factorizations, adjoint):**

The displacement exponential factors into normal or antinormal order with opposite scalar factors, and its adjoint is the displacement by minus the amplitude.

$$ \hat D(\alpha)=e^{\alpha \hat a^\dagger-\alpha^*\hat a}=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}e^{-\alpha^*\hat a}=e^{|\alpha|^2/2}e^{-\alpha^*\hat a}e^{\alpha \hat a^\dagger},\quad \hat D^\dagger(\alpha)=\hat D(-\alpha). $$

### Displacement: adjoint action on a
Appears in: X-001; P-008; G-018 (used in G-059, G-061, G-093)

Conjugating the annihilator by the displacement shifts it by the c-number $\alpha$.

$$ \hat D(\alpha)^\dagger \hat a\hat D(\alpha)=\hat a+\alpha\hat I $$

Conditions: the atlas notes that Berezin's direction $\hat U\hat a\hat U^{-1}$ reverses the sign of the parameter; with $\hat U=\hat D(\alpha)$ that direction gives $\hat a-\alpha$.

**Also stated in P-008 (Displacement conjugation of a, q, p):**

Conjugating by a displacement shifts the annihilator by $\alpha$ and the quadratures by $Q$ and $P$.

$$ \hat D(\alpha)^\dagger \hat a\hat D(\alpha)=\hat a+\alpha\hat I,\quad \hat D^\dagger \hat q\hat D=\hat q+Q\hat I,\quad \hat D^\dagger \hat p\hat D=\hat p+P\hat I $$

**Also stated in G-018 (Displacement as a translation of the field):**

Conjugation by a displacement shifts a and a^dagger by the amplitude and its conjugate; the inverse conjugation shifts the other way.

$$ \hat D^\dagger(\alpha)\hat a\hat D(\alpha)=\hat a+\alpha\hat I,\quad \hat D^\dagger(\alpha)\hat a^\dagger \hat D(\alpha)=\hat a^\dagger+\alpha^*\hat I,\quad \hat D\hat a\hat D^\dagger=\hat a-\alpha\hat I. $$

### Displacement: composition law
Appears in: X-001; P-003; C-003; same law as C-059 (one-mode re-expression in terms of Im<f,g>); G-017, G-022 (phase of composed drives)

Two displacements compose to the displacement by the sum, up to a phase set by the symplectic area of the two amplitudes.

$$ \hat D(\alpha)\hat D(\beta)=e^{(\alpha\bar\beta-\bar\alpha\beta)/2}\hat D(\alpha+\beta) $$

**Also stated in P-003 (Weyl composition law and adjoint):**

Two displacements compose to the displacement by the summed amplitude times a phase set by the symplectic area, and the adjoint of a displacement is the displacement by the opposite amplitude.

$$ \hat D(\alpha)\hat D(\beta)=e^{(\alpha\bar\beta-\bar\alpha\beta)/2}\hat D(\alpha+\beta),\qquad \hat D(\alpha)^\dagger=\hat D(-\alpha) $$

**Also stated in C-003 (Weyl translation product and phase):**

Two phase-space translations compose into the translation by the summed vector times a phase set by the symplectic area between them; reversing the order reverses that phase, and the inverse is the translation by -z.

$$ \hat T(z)\hat T(z')=e^{-i\sigma(z,z')/(2\hbar)}\hat T(z+z'),\quad \hat T(z)\hat T(z')=e^{-i\sigma(z,z')/\hbar}\hat T(z')\hat T(z),\quad \hat T(z)^{-1}=\hat T(-z) $$

**Also stated in G-017 (Weyl composition law and coherent displacement phase):**

Two displacements compose to the displacement by the summed amplitude times a phase, so displacements commute up to the phase e^{alpha beta^* - alpha^* beta}.

$$ \hat D(\alpha)\hat D(\beta)=e^{(\alpha\beta^*-\alpha^*\beta)/2}\hat D(\alpha+\beta),\quad \hat D(\alpha)|\beta\rangle=e^{i\operatorname{Im}(\alpha\beta^*)}|\alpha+\beta\rangle,\quad \hat D(\alpha)\hat D(\beta)=e^{\alpha\beta^*-\alpha^*\beta}\hat D(\beta)\hat D(\alpha). $$

### Displacement in the source direction
Appears in: B-050

Berezin's displacement U_f shifts a by +f under U a U^{-1}; it is D(-f) in the common convention, and its normal symbol is a Gaussian-normalized exponential of linear terms.

$$ \hat U_f=e^{\hat a(f^*)-\hat a^*(f)},\quad \hat U_f\hat a\hat U_f^{-1}=\hat a+f\hat I, \qquad N_{\hat U_f}=e^{-\|f\|^2/2+af^*-a^*f} $$

### Canonical commutators, ladder operators and the multimode displacement
Appears in: P-001

The canonical pair fixes the ladder operators and their unit commutator, and the multimode displacement is the exponential of the linear form in the ladder operators.

$$ [\hat q,\hat p]=i\hbar \hat I,\quad \hat a=\frac{\hat q+i\hat p}{\sqrt{2\hbar}},\quad \hat a^\dagger=\frac{\hat q-i\hat p}{\sqrt{2\hbar}},\quad [\hat a,\hat a^\dagger]=\hat I,\quad \hat D(\boldsymbol\alpha)=\exp\sum_i(\alpha_i \hat a_i^\dagger-\bar\alpha_i \hat a_i) $$

### Displacement commutation phase, group commutator and polygon phase
Appears in: P-003

Swapping two displacements costs twice the composition phase, the group commutator is that phase times the identity, and an ordered product of several displacements is the total displacement times the phase $\mathrm{Im}\sum_{j>k}\alpha_j\bar\alpha_k$.

$$ \hat D(\alpha)\hat D(\beta)=e^{\alpha\bar\beta-\bar\alpha\beta}\hat D(\beta)\hat D(\alpha),\quad \hat D(\alpha)\hat D(\beta)\hat D(-\alpha)\hat D(-\beta)=e^{\alpha\bar\beta-\bar\alpha\beta}\hat I,\quad \hat D(\alpha_n)\cdots\hat D(\alpha_1)=e^{i\,\mathrm{Im}\sum_{j>k}\alpha_j\bar\alpha_k}\hat D\big(\textstyle\sum_j\alpha_j\big) $$

### Lattice commutation condition and critical coherent lattice
Appears in: P-014, P2-050

Two displacements commute exactly when their symplectic area is a multiple of $\pi$; a coherent-state lattice is complete when its cell area is at most $\pi$ and incomplete above.

$$ [\hat D(\omega_1),\hat D(\omega_2)]=0 \iff \mathrm{Im}(\omega_1\bar\omega_2)\in\pi\mathbb Z;\qquad \text{lattice of cell area } S:\ \text{complete for } S\le\pi,\ \text{incomplete for } S>\pi $$

### Displacement in position representation
Appears in: P-006

The displacement factors into a momentum kick and a position translation with a half-area phase, so on wavefunctions it shifts the argument by $Q$ and multiplies by a plane-wave phase.

$$ \hat D(\alpha)=e^{-iPQ/(2\hbar)}e^{iP\hat q/\hbar}e^{-iQ\hat p/\hbar},\quad [\hat D(\alpha)\psi](q)=e^{iP(q-Q/2)/\hbar}\psi(q-Q),\quad \alpha=\frac{Q+iP}{\sqrt{2\hbar}} $$

### Displacement conjugation and coordinate factorization
Appears in: C-004

Conjugating by a translation shifts position and momentum by the classical displacement; the translation factors into a position-space kick and a momentum-generated shift times a scalar phase, and in position space it shifts the argument and multiplies by a plane wave.

$$ \hat T(z)\hat Q\hat T(z)^{-1}=\hat Q-q\hat I,\quad \hat T(z)=e^{-iq\cdot p/(2\hbar)}e^{ip\cdot \hat Q/\hbar}e^{-iq\cdot \hat P/\hbar},\quad (\hat T(z)u)(x)=e^{ip\cdot(x-q/2)/\hbar}u(x-q) $$

### Displacement conjugation of a field and of N
Appears in: C-060

Conjugating by a Weyl operator shifts the field a(g) by <g,f> and the number operator by the linear field terms plus ||f||^2.

$$ \hat T(f)^*\hat a(g)\hat T(f)=\hat a(g)+\langle g,f\rangle\hat I,\quad \hat T(f)^*\hat N\hat T(f)=\hat N+\hat a(f)+\hat a^\dagger(f)+\|f\|^2\hat I $$

### Derivative of a moving displacement
Appears in: C-061

The time derivative of a moving Weyl operator is the operator times the field of the velocity plus the scalar i Im<f, f'>.

$$ \frac d{dt}\hat T(f_t)=\hat T(f_t)\big[\hat a^\dagger(\dot f_t)-\hat a(\dot f_t)+i\Im\langle f_t,\dot f_t\rangle\hat I\big] $$

### Vacuum action of the ordered displacement
Appears in: G-016

In the normal-ordered displacement the annihilation factor acts trivially on the vacuum, so the coherent state is a normalized exponential of a^dagger on the vacuum.

$$ e^{-\alpha^*\hat a}|0\rangle=|0\rangle,\quad \hat D(\alpha)|0\rangle=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}|0\rangle=|\alpha\rangle,\quad \|e^{\alpha\hat a^\dagger}|0\rangle\|=e^{|\alpha|^2/2}. $$

### Recoil exponential as a displacement and its interaction-picture form
Appears in: G-101

The recoil exponential is a displacement by -i eta, and the free oscillator rotation turns its generator into the time-dependent quadrature.

$$ e^{-i\eta(\hat a+\hat a^\dagger)}=\hat D(-i\eta)=e^{-\eta^2/2}e^{-i\eta \hat a^\dagger}e^{-i\eta \hat a},\qquad\hat U_0^\dagger e^{-i\eta(\hat a+\hat a^\dagger)}\hat U_0=\exp[-i\eta(\hat ae^{-i\nu t}+\hat a^\dagger e^{i\nu t})]. $$

Conditions: real $\eta$.

## Coherent states (general)

### Poisson (coherent) vectors
Appears in: B-017

The unnormalized coherent vector e^{a*f} has squared norm e^{|f|^2}; dividing by e^{|f|^2/2} gives the normalized coherent vector, whose photon-number distribution is Poisson with mean |f|^2.

$$ P_f(a^*)=e^{a^*f-\|f\|^2/2},\qquad Q_f(a^*)=e^{a^*f},\quad \|P_f\|=1,\quad\|Q_f\|^2=e^{\|f\|^2} $$

### Overlap of displaced fiducial states
Appears in: P-007

The overlap of two displaced copies of one fiducial state is the Weyl phase times the expectation of the relative displacement in the fiducial.

$$ \langle\beta;\psi_0|\alpha;\psi_0\rangle=e^{(\alpha\bar\beta-\bar\alpha\beta)/2}\langle\psi_0|\hat D(\alpha-\beta)|\psi_0\rangle $$

Conditions: normalized fiducial; resolution of identity $\int d^2\alpha\,|\alpha;\psi_0\rangle\langle\alpha;\psi_0|/\pi=\hat I$ in the weak sense.

### Coherent state expansion, eigenvalue equation, overlap and resolution
Appears in: P-012; C-007; infinite-mode form Psi(f) in C-059; G-013, G-016 (vacuum displacement), G-021 (Poisson coefficients)

A coherent state is the normalized exponential of $\alpha\hat a^\dagger$ on the vacuum, it is an eigenvector of $\hat a$, two coherent states overlap by a Gaussian, and coherent projectors resolve the identity.

$$ |\alpha\rangle=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}|0\rangle=e^{-|\alpha|^2/2}\sum_{n}\frac{\alpha^n}{\sqrt{n!}}|n\rangle,\quad \hat a|\alpha\rangle=\alpha|\alpha\rangle,\quad \langle\alpha|\beta\rangle=e^{-(|\alpha|^2+|\beta|^2)/2+\bar\alpha\beta},\quad |\langle\alpha|\beta\rangle|^2=e^{-|\alpha-\beta|^2},\quad \int\frac{d^{2N}\alpha}{\pi^N}|\boldsymbol\alpha\rangle\langle\boldsymbol\alpha|=\hat I $$

**Also stated in C-007 (Coherent state: vacuum action, eigenvalue, number expansion):**

A coherent state is the vacuum displaced by alpha; it is an eigenvector of the annihilator with eigenvalue alpha and has amplitudes alpha^n/sqrt(n!) times e^{-|alpha|^2/2} in the number basis.

$$ |\alpha\rangle=e^{-|\alpha|^2/2}e^{\alpha\cdot \hat a^\dagger}|0\rangle,\quad \hat a_j|\alpha\rangle=\alpha_j|\alpha\rangle,\quad |\alpha\rangle=e^{-|\alpha|^2/2}\sum_{\nu}\frac{\alpha^\nu}{\sqrt{\nu!}}|\nu\rangle $$

**Also stated in G-013 (Coherent state as eigenvector and its Fock series):**

A coherent state is the right eigenvector of a, obtained by displacing the vacuum, with Poisson-weighted Fock amplitudes.

$$ \hat a|\alpha\rangle=\alpha|\alpha\rangle,\quad |\alpha\rangle=e^{-|\alpha|^2/2}\sum_{m\ge0}\frac{\alpha^m}{\sqrt{m!}}|m\rangle,\quad\langle\alpha|\hat a^\dagger=\alpha^*\langle\alpha|. $$

### Bargmann realization and displacement as weighted translation
Appears in: P-013; C-011; G-023

In the Bargmann space $\hat a^\dagger$ is multiplication by $z$ and $\hat a$ is $\partial_z$, and a displacement acts as a weighted translation of the entire function.

$$ \hat a\leftrightarrow\partial_z,\quad \hat a^\dagger\leftrightarrow z,\quad \|f\|^2=\int\frac{d^2z}{\pi}e^{-|z|^2}|f(z)|^2,\quad [\hat D(\alpha)f](z)=e^{-|\alpha|^2/2+\alpha z}f(z-\bar\alpha) $$

**Also stated in C-011 (Bargmann realization):**

In the Bargmann space of entire functions, creation acts as multiplication by zeta and annihilation as the derivative, so the linear and number exponentials act as a multiplication, a shift, and a dilation.

$$ \hat a^\dagger_j\mapsto\zeta_j,\quad \hat a_j\mapsto\partial_{\zeta_j},\quad \hat N\mapsto\zeta\cdot\partial_\zeta;\quad e^{v\cdot \hat a}f(\zeta)=f(\zeta+v),\ e^{u\cdot \hat a^\dagger}f=e^{u\cdot\zeta}f,\ e^{s\hat N}f(\zeta)=f(e^s\zeta) $$

Conditions: hbar = 1, analytic vectors.

**Also stated in G-023 (Bargmann (entire-function) representation):**

Projecting a state on a coherent bra gives a Gaussian times an entire function of alpha^*, on which a^dagger acts as multiplication by z and a as d/dz.

$$ \langle\alpha|\psi\rangle=e^{-|\alpha|^2/2}\psi(\alpha^*),\quad \psi(z)=\sum_{m\ge0}\psi_m\frac{z^m}{\sqrt{m!}},\quad \hat a^\dagger\mapsto z,\ \hat a\mapsto\partial_z. $$

### Vacuum displacement amplitude and coherent overlap
Appears in: C-008; same overlap in C-059 (as exp(-|f-g|^2/2 + i Im<f,g>)); C-010; G-019, G-104 (trace formula)

The vacuum expectation of a translation is a Gaussian in the displacement, and two coherent states overlap as a Gaussian in their separation times a symplectic phase.

$$ \langle0|\hat T(z)|0\rangle=e^{-|z|^2/(4\hbar)},\quad \langle\phi_z,\phi_{z'}\rangle=e^{i\sigma(z,z')/(2\hbar)}e^{-|z-z'|^2/(4\hbar)},\quad \hat T(z)\phi_{z'}=e^{-i\sigma(z,z')/(2\hbar)}\phi_{z+z'},\quad \langle\alpha|\beta\rangle=e^{-(|\alpha|^2+|\beta|^2)/2+\bar\alpha\cdot\beta} $$

**Also stated in C-010 (Coherent resolution of identity, trace and Hilbert-Schmidt norm):**

Coherent projectors integrated over phase space with weight (2 pi hbar)^{-n} give the identity, which turns traces and Hilbert-Schmidt norms into phase-space integrals of coherent expectation values.

$$ (2\pi\hbar)^{-n}\int|\phi_z\rangle\langle\phi_z|\,dz=\hat I,\quad \operatorname{Tr}\hat B=(2\pi\hbar)^{-n}\int\langle\phi_z,\hat B\phi_z\rangle dz,\quad \|\hat B\|_{\rm HS}^2=(2\pi\hbar)^{-n}\int\|\hat B\phi_z\|^2dz $$

Conditions: trace class, resp. Hilbert-Schmidt.

**Also stated in G-019 (Coherent overlap and resolution of identity):**

Coherent states overlap with a Gaussian in their separation and resolve the identity with measure d^2 alpha/pi.

$$ \langle\beta|\alpha\rangle=e^{-(|\alpha|^2+|\beta|^2)/2+\beta^*\alpha},\quad |\langle\beta|\alpha\rangle|^2=e^{-|\alpha-\beta|^2},\quad\int_{\mathbb C}|\alpha\rangle\langle\alpha|\frac{d^2\alpha}{\pi}=\hat I. $$

### Harmonic rotation of displacements and coherent states
Appears in: C-009

The oscillator propagator rotates phase space rigidly: it conjugates a translation into the translation by the rotated vector and sends a coherent state to the coherent state labelled e^{-it} alpha, times the zero-point phase.

$$ \hat U_t\hat T(z)\hat U_t^*=\hat T(F_tz),\quad F_t=\begin{pmatrix}\cos t&\sin t\\-\sin t&\cos t\end{pmatrix},\quad \hat U_t|\alpha\rangle=e^{-int/2}|e^{-it}\alpha\rangle $$

Conditions: U_t = exp(-i t (N + n/2)).

### Canonical oscillator realization of the CCR
Appears in: G-001

Writing the annihilation operator as a complex combination of position and momentum turns the canonical commutator into the unit bosonic commutator, and the oscillator Hamiltonian becomes the number operator plus a zero-point half.

$$ [\hat q,\hat p]=i\hbar \hat I,\qquad \hat a=\frac{\omega \hat q+i\hat p}{\sqrt{2\hbar\omega}},\qquad [\hat a,\hat a^\dagger]=\hat I,\qquad \hat H=\hbar\omega(\hat n+\tfrac12\hat I). $$

### Quadrature commutator and number-state variance
Appears in: G-007

The two quadratures have commutator i/2, so their variances obey the 1/16 product bound, and a number state has variance (2n+1)/4 in both.

$$ [\hat X_1,\hat X_2]=\frac i2\hat I,\quad \operatorname{Var}\hat X_1\operatorname{Var}\hat X_2\ge\frac1{16},\quad\operatorname{Var}_{|n\rangle}\hat X_j=\frac{2n+1}{4}. $$

### No normalizable right eigenvector of the creation operator
Appears in: G-014

The creation operator has no normalizable right eigenvector: the eigenvalue equation forces every Fock component to vanish.

$$ \hat a^\dagger|\psi\rangle=\lambda|\psi\rangle\ \Rightarrow\ |\psi\rangle=0 $$

### Coherent-state photon statistics
Appears in: G-021

A coherent state has Poisson photon statistics, factorial moments |alpha|^{2r}, and vacuum quadrature noise.

$$ P_m=e^{-|\alpha|^2}\frac{|\alpha|^{2m}}{m!},\quad \langle \hat n\rangle=\operatorname{Var}\hat n=|\alpha|^2,\quad \langle(\hat a^\dagger)^r \hat a^r\rangle=|\alpha|^{2r},\quad\operatorname{Var}\hat X_j=\frac14. $$

### Photon subtraction and addition
Appears in: G-093

Subtracting a photon conditions the mean number to <n(n-1)>/<n>; adding a photon to a coherent state gives a displaced superposition of |1> and alpha^*|0>.

$$ \langle \hat n\rangle_- =\frac{\langle \hat n(\hat n-\hat I)\rangle}{\langle \hat n\rangle}=\langle \hat n\rangle-1+\frac{\operatorname{Var}\hat n}{\langle \hat n\rangle},\qquad|\alpha;+\rangle=\frac{\hat a^\dagger|\alpha\rangle}{\sqrt{1+|\alpha|^2}}=\frac{\hat D(\alpha)(|1\rangle+\alpha^*|0\rangle)}{\sqrt{1+|\alpha|^2}}. $$

## Single-mode squeezing

### Single-mode squeeze: normal-ordered factorization
Appears in: X-002 (source records cited: P-028 to P-031, G-051 to G-061, B-051, C-048, C-049); P-031; G-054; C-048 (two-mode version: Two-mode squeezing)

The squeeze operator factors into a pair-creation exponential, a power of the number operator, and a pair-annihilation exponential, with $\tanh r$ in the outer factors.

$$ \hat S(\xi)=\exp\!\left[\tfrac12\left(\bar\xi \hat a^2-\xi(\hat a^\dagger)^2\right)\right]=e^{-\frac12e^{i\theta}\tanh r\,(\hat a^\dagger)^2}e^{-\log\cosh r\,(\hat N+\tfrac12\hat I)}e^{\frac12e^{-i\theta}\tanh r\,\hat a^2} $$

Conditions: $\xi=re^{i\theta}$, $r\ge0$.

**Also stated in P-031 (Single-mode squeezing: normal-ordered factorization):**

The one-mode squeeze operator in normal order is a pair-creation exponential, a number-diagonal factor $(\mathrm{sech}\,r)^{\hat a^\dagger\hat a+1/2}$ and a pair-annihilation exponential, so its vacuum amplitude is $(\cosh r)^{-1/2}$.

$$ \hat S(\xi)=e^{\frac12(\bar\xi\hat a^2-\xi\hat a^{\dagger2})}=e^{-\frac12e^{i\theta}\tanh r\,(\hat a^\dagger)^2}e^{-\log\cosh r\,(\hat a^\dagger \hat a+1/2)}e^{\frac12e^{-i\theta}\tanh r\,\hat a^2},\quad \xi=re^{i\theta} $$

**Also stated in G-054 (Squeeze operator: normal-ordered Gauss factorization):**

The squeeze operator factors into a normal-ordered product: an exponential of a^dagger^2, a power of sech r in the number operator, and an exponential of a^2.

$$ \hat S(\xi)=e^{-\tau \hat a^{\dagger2}/2}e^{-\ln(\cosh r)(\hat n+\tfrac12\hat I)}e^{\tau^*\hat a^2/2},\qquad\tau=e^{i\theta}\tanh r. $$

**Also stated in C-048 (Noncompact Gauss factorization of squeezers):**

A squeeze exponential factors into a pair-creation exponential, a number exponential with log cosh r, and a pair-annihilation exponential with tanh r; the same holds for two-mode pairs.

$$ e^{\xi\hat K_+-\bar\xi\hat K_-}=e^{\zeta\hat K_+}e^{-2\log(\cosh r)\hat K_0}e^{-\bar\zeta\hat K_-},\ \zeta=\frac\xi r\tanh r;\qquad \hat S(\xi)=e^{\zeta\hat a^{\dagger2}/2}e^{-\log(\cosh r)(\hat N+\frac12)}e^{-\bar\zeta\hat a^2/2} $$

and the two-mode form with a1^dag a2^dag, a1 a2, N1 + N2 + 1.

### Single-mode squeeze: adjoint action on a
Appears in: X-002; G-052

Squeezing mixes the annihilator with the creator, with $\cosh r$ and $\sinh r$ weights.

$$ \hat S^\dagger \hat a\hat S=\hat a\cosh r-e^{i\theta}\hat a^\dagger\sinh r $$

Conditions: the atlas notes that $\hat S^{-1}$ flips the sign of every squeeze parameter.

**Also stated in G-052 (Single-mode Bogoliubov transformation):**

Squeezing transforms a into the Bogoliubov combination mu a - nu a^dagger with mu = cosh r and nu = e^{i theta} sinh r.

$$ \hat S^\dagger \hat a\hat S=\mu \hat a-\nu \hat a^\dagger,\quad\hat S^\dagger \hat a^\dagger \hat S=\mu \hat a^\dagger-\nu^*\hat a,\quad\mu=\cosh r,\ \nu=e^{i\theta}\sinh r. $$

### Single-mode squeezed vacuum
Appears in: X-002; G-056 (also G-062 at $\alpha=0$, G-063 overlap)

The squeezed vacuum is a normalized Gaussian superposition of even photon numbers.

$$ \hat S|0\rangle=(\cosh r)^{-1/2}e^{-\frac12e^{i\theta}\tanh r\,(\hat a^\dagger)^2}|0\rangle $$

**Also stated in G-056 (Squeezed vacuum: nullifier and Fock expansion):**

Squeezed vacuum is annihilated by the Bogoliubov combination mu a + nu a^dagger and contains only even photon numbers with the stated amplitudes.

$$ (\mu \hat a+\nu \hat a^\dagger)|\xi\rangle=0,\quad |\xi\rangle=\frac1{\sqrt{\cosh r}}\sum_{m\ge0}\frac{\sqrt{(2m)!}}{2^m m!}(-e^{i\theta}\tanh r)^m|2m\rangle,\quad P_{2m}=\frac{(2m)!\tanh^{2m}r}{2^{2m}(m!)^2\cosh r}. $$

### One-mode squeezing implementer and kernel
Appears in: B-051

The one-mode squeeze exp[(sigma/2)(conj(eps) a^2 - eps a*^2)] maps a to a cosh sigma + eps a* sinh sigma, and its kernel is a Gaussian with coefficients tanh sigma and sech sigma.

$$ \hat U_\sigma=\exp\!\left[\tfrac\sigma2(\bar\epsilon\hat a^2-\epsilon\hat a^{*2})\right], \quad K_{\hat U_\sigma}=(\cosh\sigma)^{-1/2} e^{\frac12\bar\epsilon\tanh\sigma\,a^2+\operatorname{sech}\sigma\,a^*a-\frac12\epsilon\tanh\sigma\,a^{*2}} $$

Conditions: b = a cosh sigma + epsilon a* sinh sigma, |epsilon| = 1, sigma real.

### Single-mode su(1,1) realization and Casimir
Appears in: P-031, P2-012 (parity sectors $k=1/4,3/4$); G-053

Half the squares of the ladder operators and half the shifted number operator realize su(1,1) on one mode, with Casimir $-3/16$, so the even and odd sectors carry $k=1/4$ and $3/4$.

$$ \hat K_+=\tfrac12(\hat a^\dagger)^2,\quad \hat K_-=\tfrac12\hat a^2,\quad \hat K_0=\tfrac12(\hat a^\dagger \hat a+\tfrac12\hat I),\quad \hat C=-\tfrac{3}{16}\hat I $$

**Also stated in G-053 (Quadratic SU(1,1) algebra of one mode):**

Half the squares of the ladder operators, with half of n + 1/2, close into the su(1,1) algebra.

$$ \hat K_+=\tfrac12\hat a^{\dagger2},\ \hat K_-=\tfrac12\hat a^2,\ \hat K_0=\tfrac12(\hat n+\tfrac12\hat I),\quad[\hat K_0,\hat K_\pm]=\pm \hat K_\pm,\quad[\hat K_-,\hat K_+]=2\hat K_0. $$

### Multimode squeeze operator and bounded squeeze coordinate
Appears in: C-025

The multimode squeeze operator is the exponential of a quadratic in creators minus its adjoint; its inverse flips the sign of B, and B and the bounded coordinate W are related by tanh and the arctanh series.

$$ \hat S(B)=\exp\tfrac12(\hat a^\dagger\cdot B\hat a^\dagger-\hat a\cdot B^*\hat a),\quad \hat S(B)^{-1}=\hat S(-B),\quad B=U\operatorname{arctanh}|W|=W\sum_k\frac{(W^*W)^k}{2k+1},\quad W=U\tanh|B| $$

Conditions: B complex symmetric, W^*W < I.

### Squeeze conjugation in block form
Appears in: C-026; one-mode case also in C-049 (first display)

Squeezing acts on the column (a, a^dag) by a block matrix built from W; for one mode it gives a cosh r - e^{i theta} a^dag sinh r.

$$ \hat S(B)\binom{\hat a}{\hat a^\dagger}\hat S(-B)=\begin{pmatrix}(I-WW^*)^{-1/2}&-W(I-W^*W)^{-1/2}\\-(I-W^*W)^{-1/2}W^*&(I-W^*W)^{-1/2}\end{pmatrix}\binom{\hat a}{\hat a^\dagger},\qquad \hat S\hat a\hat S^{-1}=\hat a\cosh r-e^{i\theta}\hat a^\dagger\sinh r $$

### Squeezed number operator
Appears in: C-049

Squeezing maps the annihilator to (a - omega a^dag)/sqrt(1 - |omega|^2) and rotates K0 into a hyperbolic combination of K0, K+ and K-.

$$ \hat S\hat a\hat S^{-1}=\frac{\hat a-\omega\hat a^\dagger}{\sqrt{1-|\omega|^2}},\quad \hat S\hat K_0\hat S^{-1}=\cosh(2r)\hat K_0-\tfrac12\sinh(2r)(e^{i\theta}\hat K_++e^{-i\theta}\hat K_-) $$

### Single-mode squeeze operator
Appears in: G-051

The squeeze operator is the exponential of a quadratic anti-Hermitian generator, so its adjoint is the squeeze with opposite parameter, and it squeezes the vacuum.

$$ \hat S(\xi)=\exp\!\left[\tfrac12(\xi^*\hat a^2-\xi \hat a^{\dagger2})\right],\quad \hat S^\dagger(\xi)=\hat S(-\xi),\quad|\xi\rangle=\hat S(\xi)|0\rangle. $$

### Squeeze operator: reversed (antinormal) Gauss factorization
Appears in: G-055

The squeeze operator also factors in reversed order, with the middle number exponential changing sign.

$$ \hat S(\xi)=e^{\tau^*\hat a^2/2}e^{+\ln(\cosh r)(\hat n+\tfrac12\hat I)}e^{-\tau \hat a^{\dagger2}/2}. $$

### Rotated-quadrature variance of squeezed vacuum
Appears in: G-057, G-073 (collective mode)

The quadrature at angle vartheta of squeezed vacuum has variance (cosh 2r - sinh 2r cos(2 vartheta - theta))/4.

$$ \hat X(\vartheta)=\frac{\hat ae^{-i\vartheta}+\hat a^\dagger e^{i\vartheta}}2,\quad\operatorname{Var}_{\xi}\hat X(\vartheta)=\frac14[\cosh2r-\sinh2r\cos(2\vartheta-\theta)]. $$

### Quadrature covariance of squeezed vacuum
Appears in: G-058

The two Cartesian quadratures of squeezed vacuum are correlated by -sinh 2r sin theta/4, and the covariance determinant stays at the minimum 1/16.

$$ \operatorname{Cov}(\hat X_1,\hat X_2)=\tfrac12\langle\{\Delta \hat X_1,\Delta \hat X_2\}\rangle=-\frac14\sinh2r\sin\theta,\qquad\det V=\frac1{16}. $$

### Displaced squeezed state: low moments
Appears in: G-059

Displacing squeezed vacuum shifts the first moment by alpha, adds alpha^2 to <a^2>, and adds |alpha|^2 to the mean photon number.

$$ |\alpha,\xi\rangle=\hat D(\alpha)\hat S(\xi)|0\rangle,\quad\langle \hat a\rangle=\alpha,\quad\langle \hat a^2\rangle=\alpha^2-e^{i\theta}\sinh r\cosh r,\quad\langle \hat n\rangle=|\alpha|^2+\sinh^2r. $$

### Interchanging displacement and squeezing
Appears in: G-060

Squeezing and displacement commute up to a Bogoliubov change of the displacement amplitude, with no extra phase.

$$ \hat S(\xi)\hat D(\beta)=\hat D(\mu\beta-\nu\beta^*)\hat S(\xi),\quad \hat D(\alpha)\hat S(\xi)=\hat S(\xi)\hat D(\mu\alpha+\nu\alpha^*). $$

### Displaced squeezed state as a Bogoliubov eigenstate
Appears in: G-061

A displaced squeezed state is the eigenstate of the Bogoliubov combination mu a + nu a^dagger with eigenvalue gamma = mu alpha + nu alpha^*, which fixes the rotated-quadrature eigenvalue.

$$ (\mu \hat a+\nu \hat a^\dagger)|\alpha,\xi\rangle=\gamma|\alpha,\xi\rangle,\quad\gamma=\mu\alpha+\nu\alpha^*,\qquad(\hat Y_1+ie^{-2r}\hat Y_2)|\alpha,\xi\rangle=e^{-r-i\theta/2}\gamma|\alpha,\xi\rangle. $$

### Number variance of a displaced squeezed vacuum
Appears in: G-065

The number variance of a displaced squeezed vacuum is the coherent term |mu alpha - nu alpha^*|^2 plus the squeezed-vacuum term 2 sinh^2 r cosh^2 r.

$$ \operatorname{Var}\hat n=|\mu\alpha-\nu\alpha^*|^2+2\sinh^2r\cosh^2r. $$

### Squared-field quadratures and their commutators
Appears in: G-089

The squared-field quadratures close with K_3 into su(1,1) with a state-dependent commutator, giving the bound Var Q_1 Var Q_2 >= 4 <K_3>^2.

$$ [\hat Q_1,\hat Q_2]=-4i\hat K_3,\quad[\hat K_3,\hat Q_1]=i\hat Q_2,\quad[\hat Q_2,\hat K_3]=i\hat Q_1,\qquad\operatorname{Var}\hat Q_1\operatorname{Var}\hat Q_2\ge4\langle \hat K_3\rangle^2. $$

Conditions: $\hat Q_1=(\hat a^{\dagger2}+\hat a^2)/2$, $\hat Q_2=(\hat a^{\dagger2}-\hat a^2)/(2i)$, $\hat K_3=(\hat n+\tfrac12\hat I)/2$.

## Two-mode squeezing

### Two-mode squeeze: normal-ordered factorization
Appears in: X-003 (source records cited: P-032, G-068 to G-075, P-029); P-032, P2-019; G-070, G-071 (reversed form)

The two-mode squeezer factors the same way as the single-mode one, with the pair $\hat a^\dagger\hat b^\dagger$ in place of $(\hat a^\dagger)^2/2$.

$$ \hat S_2(\xi)=e^{\bar\xi \hat a\hat b-\xi \hat a^\dagger \hat b^\dagger}=e^{-e^{i\theta}\tanh r\,\hat a^\dagger \hat b^\dagger}e^{-\log\cosh r\,(\hat N_a+\hat N_b+\hat I)}e^{e^{-i\theta}\tanh r\,\hat a\hat b} $$

Conditions: closed algebra $\hat K_+=\hat a^\dagger\hat b^\dagger$, $\hat K_-=\hat a\hat b$, $\hat K_0=(\hat N_a+\hat N_b+\hat I)/2$; $\hat N_a-\hat N_b$ conserved; P-032's convention is $\xi\mapsto-\xi$.

**Also stated in P-032 (Two-mode squeezing: normal-ordered factorization and pair state):**

The two-mode squeeze operator in normal order is a pair-creation exponential, $(\mathrm{sech}\,r)^{\hat N_a+\hat N_b+1}$ and a pair-annihilation exponential; on the vacuum it gives a pair state with geometric pair-number distribution.

$$ e^{r(e^{i\phi}\hat a^\dagger \hat b^\dagger-e^{-i\phi}\hat a\hat b)}=e^{e^{i\phi}\tanh r\,\hat a^\dagger \hat b^\dagger}e^{-\log\cosh r\,(\hat N_a+\hat N_b+\hat I)}e^{-e^{-i\phi}\tanh r\,\hat a\hat b};\qquad |\zeta\rangle=\sqrt{1-|\zeta|^2}\,e^{\zeta\hat a^\dagger\hat b^\dagger}|0,0\rangle,\ W_n=(1-\rho)\rho^n,\ \rho=|\zeta|^2 $$

**Also stated in G-070 (Two-mode squeeze: normal and reversed Gauss factorizations):**

The two-mode squeeze factors in normal order and in reversed order, the middle factor being a power of cosh r in the total number plus one.

$$ \hat S_2(\xi)=e^{-\tau \hat a^\dagger \hat b^\dagger}e^{-\ln(\cosh r)(\hat n_a+\hat n_b+\hat I)}e^{\tau^*\hat a\hat b}=e^{\tau^*\hat a\hat b}e^{+\ln(\cosh r)(\hat n_a+\hat n_b+\hat I)}e^{-\tau \hat a^\dagger \hat b^\dagger}. $$

### Two-mode squeezed vacuum
Appears in: X-003; G-072

The two-mode squeezed vacuum is a sum of twin-photon-number states with geometric weights.

$$ \hat S_2|0,0\rangle=\operatorname{sech}r\sum_{n\ge0}(-e^{i\theta}\tanh r)^n|n,n\rangle $$

**Also stated in G-072 (Twin-beam vacuum: Fock expansion and nullifiers):**

The twin-beam vacuum is a sum over equal photon pairs with weights (-e^{i theta} tanh r)^m/cosh r, annihilated by two Bogoliubov combinations and by the number difference.

$$ |\xi\rangle_2=\frac1{\cosh r}\sum_{m\ge0}(-e^{i\theta}\tanh r)^m|m,m\rangle,\quad(\mu \hat a+\nu \hat b^\dagger)|\xi\rangle_2=(\mu \hat b+\nu \hat a^\dagger)|\xi\rangle_2=0,\quad(\hat n_a-\hat n_b)|\xi\rangle_2=0. $$

### Two-mode su(1,1) realization and Casimir
Appears in: P-032, P2-019; G-069

The pair operators of two modes realize su(1,1), and the Casimir is fixed by the number difference, so a sector with difference $d$ has $k=(1+|d|)/2$.

$$ \hat K_+=\hat a^\dagger \hat b^\dagger,\quad \hat K_-=\hat a\hat b,\quad \hat K_0=\tfrac12(\hat N_a+\hat N_b+\hat I),\quad \hat C=\tfrac14[(\hat N_a-\hat N_b)^2-\hat I] $$

**Also stated in G-069 (Two-mode SU(1,1) algebra and conserved number difference):**

The pair operators a^dagger b^dagger, a b and half the total number plus one close into su(1,1), and the number difference commutes with all of them.

$$ \hat K_+=\hat a^\dagger \hat b^\dagger,\ \hat K_-=\hat a\hat b,\ \hat K_0=(\hat n_a+\hat n_b+\hat I)/2,\quad[\hat K_0,\hat K_\pm]=\pm \hat K_\pm,\ [\hat K_-,\hat K_+]=2\hat K_0,\ [\hat n_a-\hat n_b,\hat K_\pm]=0. $$

### Multimode paired vacuum normalization
Appears in: P2-021

The multimode pair-creation exponential on the vacuum has squared norm $\det(I_N-Z^\dagger Z)^{-1}$, so the stated prefactor normalizes it.

$$ |Z\rangle=\det(I_N-Z^\dagger Z)^{1/2}\exp\left(\sum_{ij}Z_{ij}\hat a_i^\dagger \hat b_j^\dagger\right)|0,0\rangle,\quad I_N-Z^\dagger Z>0 $$

### Multimode paired vacuum overlap
Appears in: P2-022

Two multimode paired vacua overlap by the product of square-root normalizations over $\det(I_N-Z_1^\dagger Z_2)$.

$$ \langle Z_1|Z_2\rangle=\frac{\det(I_N-Z_1^\dagger Z_1)^{1/2}\det(I_N-Z_2^\dagger Z_2)^{1/2}}{\det(I_N-Z_1^\dagger Z_2)},\qquad Z'=(AZ-B)(-CZ+D)^{-1} $$

### Bosonic pair-number generating determinant (trace-log sign corrected)
Appears in: P2-023

The pair-number distribution of a multimode paired vacuum has generating function $\det(I-Z^\dagger Z)/\det(I-tZ^\dagger Z)$, and its logarithm is the trace series with a plus sign.

$$ \mathcal G(t)=\sum_nW_nt^n=\frac{\det(I_N-Z^\dagger Z)}{\det(I_N-tZ^\dagger Z)},\qquad \log\mathcal G=\log\det(I_N-Z^\dagger Z)+\sum_k\frac{\operatorname{Tr}[(Z^\dagger Z)^k]t^k}{k} $$

### Bosonic realizations of SU(1,1)
Appears in: C-046

Squares of one mode and pair operators of two modes realize su(1,1); the one-mode Casimir is -3/16 and the two-mode Casimir is ((N1 - N2)^2 - 1)/4.

$$ \hat K_+=\tfrac12\hat a^{\dagger2},\ \hat K_-=\tfrac12\hat a^2,\ \hat K_0=\tfrac12(\hat N+\tfrac12),\ \hat C=-\tfrac3{16};\qquad \hat K_+=\hat a_1^\dagger\hat a_2^\dagger,\ \hat K_-=\hat a_1\hat a_2,\ \hat K_0=\tfrac12(\hat N_1+\hat N_2+1),\ \hat C=\tfrac14[(\hat N_1-\hat N_2)^2-1] $$

### Two-mode squeeze operator and Bogoliubov action
Appears in: G-068

The two-mode squeeze mixes each mode with the conjugate of the other: a -> mu a - nu b^dagger and b -> mu b - nu a^dagger.

$$ \hat S_2(\xi)=e^{\xi^*\hat a\hat b-\xi \hat a^\dagger \hat b^\dagger},\quad\hat S_2^\dagger \hat a\hat S_2=\mu \hat a-\nu \hat b^\dagger,\quad\hat S_2^\dagger \hat b\hat S_2=\mu \hat b-\nu \hat a^\dagger. $$

### Two-mode squeeze as two opposite single-mode squeezes of the sum and difference modes
Appears in: G-071

In the sum and difference modes the two-mode squeeze is a product of commuting single-mode squeezes with opposite parameters.

$$ \hat S_{2,ab}(\xi)=\hat S_c(\xi)\hat S_d(-\xi),\qquad\hat c=(\hat a+\hat b)/\sqrt2,\ \hat d=(\hat a-\hat b)/\sqrt2. $$

### Twin-beam number statistics, collective squeezing, and the intermode Cauchy-Schwarz test
Appears in: G-073, G-091

Each twin-beam mode has mean sinh^2 r and Bose-Einstein variance, the two photon numbers are perfectly correlated, the normally ordered Cauchy-Schwarz inequality fails, and the collective quadratures are squeezed by e^{-2r}.

$$ \langle \hat n_a\rangle=N=\sinh^2r,\quad\operatorname{Var}\hat n_a=\operatorname{Cov}(\hat n_a,\hat n_b)=N(N+1),\quad\langle \hat a^{\dagger2}\hat a^2\rangle\langle \hat b^{\dagger2}\hat b^2\rangle=4N^4<(2N^2+N)^2=\langle \hat n_a\hat n_b\rangle^2. $$

### Pair coherent states
Appears in: G-090

Pair coherent states are eigenstates of a b with fixed number difference q; their amplitudes follow from a one-step recurrence and their normalization is a modified Bessel function.

$$ |\eta,q\rangle=\left[\sum_{m\ge0}\frac{|\eta|^{2m}}{m!(m+q)!}\right]^{-1/2}\sum_{m\ge0}\frac{\eta^m}{\sqrt{m!(m+q)!}}|m+q,m\rangle,\qquad\sum_{m\ge0}\frac{|\eta|^{2m}}{m!(m+q)!}=|\eta|^{-q}I_q(2|\eta|). $$

### Low-gain down-conversion to second order
Appears in: G-095

At low gain the down-conversion propagator on vacuum is vacuum, minus i mu times one pair, minus mu^2 times two pairs, with the norm correction -|mu|^2/2 on vacuum.

$$ \hat U|0,0\rangle=(1-|\mu|^2/2)|0,0\rangle-i\mu|1,1\rangle-\mu^2|2,2\rangle+O(t^3),\qquad\mu=\eta t. $$

Conditions: $\hat H_I=\hbar(\eta\hat a_s^\dagger\hat a_i^\dagger+\eta^*\hat a_s\hat a_i)$.

## SU(2) and passive mode mixing

### Passive transformation of normal symbols
Appears in: B-053

A passive (number-conserving) transformation acts on normal symbols by substituting conj(Phi) a* and Phi a, because it never mixes a with a* and so keeps normal order.

$$ N_{\hat U\hat A\hat U^{-1}}(a^*,a)=N_{\hat A}(\bar{\hat\Phi} a^*,\hat\Phi a) $$

Conditions: Psi = f = 0, Phi unitary.

### Exponential of an su(2) matrix
Appears in: C-039

A 2-by-2 su(2) matrix squares to minus its determinant, so its exponential is a cosine plus a sine term.

$$ e^{tA}=I_2\cos t+A\sin t\ (\det A=1),\qquad e^{it\mathbf v\cdot\boldsymbol\sigma}=\cos(t|\mathbf v|)I_2+i\sin(t|\mathbf v|)\mathbf n\cdot\boldsymbol\sigma $$

### Second quantization of a one-particle unitary
Appears in: C-058

Second quantization turns the one-particle unitary e^{-ith} into the exponential of dGamma(h); a number-conserving quadratic exponential acts on creators by the one-particle matrix and maps Psi(f) to Psi(Uf).

$$ \widehat\Gamma(e^{-it\hat h})=e^{-it\,\widehat{d\Gamma}(\hat h)},\quad \widehat{d\Gamma}(\hat I)=\hat N,\quad \widehat\Gamma(\hat U)\Psi(f)=\Psi(\hat Uf) $$

### Atomic transition algebra and Jaynes-Cummings excitation number
Appears in: G-037

The atomic raising and lowering operators close on sigma_3, and the Jaynes-Cummings Hamiltonian conserves the total excitation number of field plus atom.

$$ [\hat\sigma_+,\hat\sigma_-]=\hat\sigma_3,\quad[\hat\sigma_3,\hat\sigma_\pm]=\pm2\hat\sigma_\pm,\quad[\hat H_{\rm JC},\hat n+|e\rangle\langle e|]=0. $$

### Detuned Jaynes-Cummings doublet propagator
Appears in: G-039

Detuning turns each doublet into a spin-1/2 precession at the generalized Rabi frequency Omega_m, because the doubled doublet Hamiltonian squares to Omega_m^2.

$$ \hat h_m/\hbar=\frac12\begin{pmatrix}\Delta&2\lambda\sqrt{m+1}\\2\lambda\sqrt{m+1}&-\Delta\end{pmatrix},\quad\Omega_m=\sqrt{\Delta^2+4\lambda^2(m+1)},\quad e^{-i\hat h_mt/\hbar}=\hat I_2\cos(\Omega_mt/2)-i\frac{2\hat h_m}{\hbar\Omega_m}\sin(\Omega_mt/2). $$

### Passive two-port transform must be unitary
Appears in: G-043

A linear two-port keeps the output commutators bosonic only if its transfer matrix is unitary.

$$ [\hat a_{{\rm out},j},\hat a_{{\rm out},k}^\dagger]=(TT^\dagger)_{jk}\hat I,\quad\text{hence }TT^\dagger=I_2. $$

### Beam splitter with i-reflection: exponential and adjoint action
Appears in: G-044, G-086 (50:50 case $\vartheta=\pi/2$)

The symmetric beam splitter exp(i theta J_1) mixes the two modes with cos(theta/2) transmission and i sin(theta/2) reflection.

$$ \hat B_1(\vartheta)=e^{i\vartheta(\hat a^\dagger \hat b+\hat b^\dagger \hat a)/2},\quad\hat B_1^\dagger\begin{pmatrix}\hat a\\\hat b\end{pmatrix}\hat B_1=\begin{pmatrix}\cos(\vartheta/2)&i\sin(\vartheta/2)\\i\sin(\vartheta/2)&\cos(\vartheta/2)\end{pmatrix}\begin{pmatrix}\hat a\\\hat b\end{pmatrix}. $$

### Real beam-splitter rotation
Appears in: G-045

The real beam splitter exp(i theta J_2) is a plain rotation of the mode pair by theta/2.

$$ \hat B_2(\vartheta)=e^{i\vartheta(\hat a^\dagger \hat b-\hat b^\dagger \hat a)/(2i)},\quad \hat B_2^\dagger\begin{pmatrix}\hat a\\\hat b\end{pmatrix}\hat B_2=\begin{pmatrix}c&s\\-s&c\end{pmatrix}\begin{pmatrix}\hat a\\\hat b\end{pmatrix},\quad c=\cos(\vartheta/2),\ s=\sin(\vartheta/2). $$

### Schwinger SU(2) algebra, Casimir and J3 phase action
Appears in: G-046

Two modes realize SU(2) through the Schwinger generators, with Casimir J_0(J_0+1); a J_3 rotation gives the two modes opposite half-phases.

$$ \hat J_3=(\hat n_a-\hat n_b)/2,\ \hat J_0=(\hat n_a+\hat n_b)/2,\ [\hat J_i,\hat J_j]=i\epsilon_{ijk}\hat J_k,\ \hat{\mathbf J}^2=\hat J_0(\hat J_0+\hat I);\quad e^{-i\varphi\hat J_3}\hat a e^{i\varphi\hat J_3}=e^{i\varphi/2}\hat a,\ \hat b\mapsto e^{-i\varphi/2}\hat b. $$

### Beam-splitter output of a number state (binomial state)
Appears in: G-047

Sending N photons into one port of a beam splitter produces the binomial superposition over the two outputs.

$$ \hat B_1(\vartheta)|0,N\rangle=\sum_{k=0}^N\sqrt{\binom Nk}\,[i\sin(\vartheta/2)]^k[\cos(\vartheta/2)]^{N-k}|k,N-k\rangle. $$

### Coherent products stay products under passive mixing
Appears in: G-048

Passive mixing maps a product of coherent states to a product of coherent states with amplitudes transformed by the beam-splitter matrix.

$$ \hat B_1(\vartheta)|\alpha,\beta\rangle=|c\alpha+is\beta,\;is\alpha+c\beta\rangle. $$

### Two-photon interference and twin-Fock output at 50:50
Appears in: G-049

At 50:50, two single photons leave together in one port (Hong-Ou-Mandel), and N photons in each port give only even occupations with the stated weights times i^N.

$$ \hat B_1|1,1\rangle=\frac i{\sqrt2}(|2,0\rangle+|0,2\rangle),\qquad \hat B_1|N,N\rangle=i^N\sum_{k=0}^N\frac{\sqrt{(2k)!(2N-2k)!}}{2^N k!(N-k)!}|2k,2N-2k\rangle. $$

Conditions: 50:50 ($\vartheta=\pi/2$).

### Cross-Kerr gate conjugated by a 50:50 beam splitter
Appears in: G-086

A cross-Kerr gate between b and c, conjugated by a 50:50 beam splitter on a and b, becomes a total-number phase times a c-controlled real rotation of a and b.

$$ \hat B^\dagger \hat U_K\hat B=e^{i\eta \hat n_c(\hat n_a+\hat n_b)/2}e^{(\eta \hat n_c/2)(\hat a^\dagger \hat b-\hat b^\dagger \hat a)},\qquad\hat B=e^{i\pi(\hat a^\dagger \hat b+\hat a\hat b^\dagger)/4},\ \hat U_K=e^{i\eta \hat n_b \hat n_c}. $$

### Frequency conversion as an SU(2) rotation
Appears in: G-092

Frequency conversion with U = e^{-tG} is a passive SU(2) rotation of modes b and c by angle g t with the pump phase on the cross terms.

$$ \hat U=e^{-t\hat G},\ \hat G=\eta\hat b^\dagger\hat c-\eta^*\hat b\hat c^\dagger,\quad\hat U^\dagger \hat b\hat U=\hat b\cos gt-e^{i\varphi}\hat c\sin gt,\quad\hat U^\dagger \hat c\hat U=\hat c\cos gt+e^{-i\varphi}\hat b\sin gt. $$

Conditions: $\hat H_I=-i\hbar\hat G$, $\eta=ge^{i\varphi}$.

## Bosonic canonical transformations

### Multimode canonical action and its constraints
Appears in: X-009 (source records cited: B-049, B-060, B-061, B-062, P-036, P-042, P-043); P-034

A Gaussian unitary maps the annihilators to linear combinations of annihilators and creators, and the two matrix constraints say that the new operators obey the same commutation (or anticommutation) relations.

$$ \hat U\hat{\mathbf d}\hat U^{-1}=u\hat{\mathbf d}+v\hat{\mathbf d}^{\dagger},\qquad uu^\dagger-\eta vv^\dagger=I_N,\quad uv^T=\eta vu^T $$

Conditions: $\eta=+1$ bosons, $-1$ fermions; $\hat{\mathbf d}^\dagger$ a column of creators.

**Also stated in P-034 (Finite-mode bosonic Bogoliubov conditions):**

A linear canonical map of $N$ bosonic modes preserves the indefinite form $K$, which fixes the block relations among $U$ and $V$, the inverse, and makes $Z=U^{-1}V$ symmetric and a strict contraction.

$$ \hat{\mathbf b}=U\hat{\mathbf a}+V\hat{\mathbf a}^\dagger:\ UU^\dagger-VV^\dagger=I_N,\ UV^T=VU^T,\ MKM^\dagger=K,\ M^{-1}=\begin{pmatrix}U^\dagger&-V^T\\-V^\dagger&U^T\end{pmatrix},\ U^\dagger U-V^T\bar V=I_N,\ U^\dagger V=V^T\bar U,\ Z=U^{-1}V=Z^T,\ I_N-ZZ^\dagger=(U^\dagger U)^{-1} $$

### Multimode creator-number-annihilator factorization
Appears in: X-009

Any Gaussian unitary with invertible $u$ factors into a pair-creation exponential, a number-conserving exponential, and a pair-annihilation exponential, times its vacuum amplitude.

$$ \hat U=c_0\exp\!\left[-\tfrac12\hat{\mathbf d}^{\dagger T}Z\hat{\mathbf d}^{\dagger}\right]\exp\!\left[\hat{\mathbf d}^{\dagger T}h\hat{\mathbf d}\right]\exp\!\left[\tfrac\eta2\hat{\mathbf d}^{T}W\hat{\mathbf d}\right],\quad Z=u^{-1}v,\ W=\bar v u^{-1},\ e^h=u^{-1},\ c_0=\langle0|\hat U|0\rangle,\ |c_0|=\det(uu^\dagger)^{-\eta/4} $$

Conditions: $u$ invertible (always so for bosons, since $uu^\dagger=I+vv^\dagger$; fermionic $u$ can be singular, e.g. the pair rotation of X-004 at $r=\pi/2$). The phase of $c_0$ is not fixed by $u,v$.

### Bosonic coherent-state kernel of the factorized implementer
Appears in: X-009

Between unnormalized coherent states the factorized implementer has a Gaussian kernel whose quadratic form is read off from $Z$, $u^{-1}$ and $W$.

$$ K_{\hat U}(x,y)=c_0\exp\!\left[-\tfrac12x^TZx+x^Tu^{-1}y+\tfrac12y^TWy\right] $$

Conditions: kernel taken between unnormalized coherent states, $x=\bar z$, $y=w$ in $\langle z|\hat U|w\rangle$ with $|w\rangle=e^{w\cdot\hat a^\dagger}|0\rangle$. The atlas does not state this convention; it is the one that makes the formula hold.

### Bosonic affine canonical transformation
Appears in: B-042

A bosonic canonical transformation sends a to a linear combination of a and a* plus a constant shift; the transformed creator is the adjoint of that expression.

$$ \hat b=\hat\Phi\hat a+\hat\Psi\hat a^*+f\hat I,\qquad \hat b^*=\bar{\hat\Psi}\hat a+\bar{\hat\Phi}\hat a^*+f^*\hat I $$

### Bosonic block relations and inverse
Appears in: B-043

The coefficient blocks of a bosonic canonical transformation satisfy Phi Phi* - Psi Psi* = I and Phi Psi' = Psi Phi', which make the doubled matrix invertible with the displayed inverse.

$$ \hat S=\begin{pmatrix}\hat\Phi&\hat\Psi\\\bar{\hat\Psi}&\bar{\hat\Phi}\end{pmatrix},\quad \hat S^{-1}=\begin{pmatrix}\hat\Phi^*&-\hat\Psi'\\-\hat\Psi^*&\hat\Phi'\end{pmatrix},\quad \hat\Phi\hat\Phi^*-\hat\Psi\hat\Psi^*=\hat I,\ \hat\Phi\hat\Psi'=\hat\Psi\hat\Phi' $$

### Symplectic form and inverse in quadratures
Appears in: B-044

In real quadrature coordinates the transformation matrix is symplectic, and its inverse is obtained by transposing and rearranging the blocks.

$$ M^TJM=J,\qquad M^{-1}=\begin{pmatrix}D^T&-B^T\\-C^T&A^T\end{pmatrix} $$

### Squeezing contraction, both charts
Appears in: B-045

The contractions Z = Phi^{-1} Psi and W = conj(Psi) Phi^{-1} that parametrize the squeezed vacuum have norm below 1, with Z symmetric; the unbarred Psi Phi^{-1} need not.

$$ \hat Z=\hat\Phi^{-1}\hat\Psi,\ \hat Z\hat Z^*=\hat I-(\hat\Phi^*\hat\Phi)^{-1};\qquad \hat W=\bar{\hat\Psi}\hat\Phi^{-1},\ \hat W^*\hat W=\hat I-(\hat\Phi\hat\Phi^*)^{-1} $$

### Full bosonic Gaussian kernel
Appears in: B-049 (affine form also used by B-071)

The matrix kernel of the implementer is the vacuum normalization c times a Gaussian whose quadratic part is built from Phi^{-1} and Psi and whose linear part carries the shift.

$$ K_{\hat U}(a^*,a)=c\exp\!\left[ \tfrac12(a\ a^*) \begin{pmatrix}\bar{\hat\Psi}\hat\Phi^{-1}&\hat\Phi'^{-1}\\\hat\Phi^{-1}&-\hat\Phi^{-1}\hat\Psi\end{pmatrix} \binom a{a^*} +a(f^*-\bar{\hat\Psi}\hat\Phi^{-1}f)-a^*\hat\Phi^{-1}f\right], \quad N_{\hat U}=K_{\hat U}e^{-a^*a} $$

### General bosonic transformation kernel for normal symbols
Appears in: B-054

For a mixing transformation, the conjugated normal symbol is a limit of Gaussian integrals of the original symbol against a kernel built from the displayed matrix L.

$$ N_{\hat U\hat A\hat U^{-1}}=\lim_n\int\mathcal K_n(a^*,a\mid\alpha^*,\alpha)N_{\hat A}(\alpha^*,\alpha)\,d^nt\,d^n\tau,\qquad \widehat{\mathcal L}=\begin{pmatrix}\bar{\hat\Phi}\bar{\hat\Psi}^{-1}&-\hat I\\-\hat I&\hat\Phi\hat\Psi^{-1}\end{pmatrix} $$

### Gaussian normal symbols stay Gaussian under conjugation
Appears in: B-055

Conjugating an operator whose normal symbol is Gaussian gives another Gaussian normal symbol, with prefactor det C^{-1/2} for the displayed block matrix C.

$$ \hat C=\begin{pmatrix}\bar{\hat\Phi}+\hat A_{11}\hat\Psi&(\hat A_{12}-\hat I)\bar{\hat\Psi}\\(\hat A_{21}-\hat I)\hat\Psi&\hat\Phi+\hat A_{22}\bar{\hat\Psi}\end{pmatrix},\qquad \text{prefactor }\det\hat C^{-1/2} $$

### Bosonic canonical diagonalization
Appears in: B-084

A stable quadratic Hamiltonian (S plus or minus T positive) is a free Hamiltonian in squeezed modes, with frequency sqrt(S^2 - T^2) in one mode and ground-energy shift (Omega - S)/2.

$$ \hat H=\tfrac12(\hat a^*\hat T\hat a^*+\hat a\bar{\hat T}\hat a+2\hat a^*\hat S\hat a)\ \to\ \hat a^*\hat\Omega\hat a+k\hat I,\qquad \hat S\pm \hat T>\mu \hat I>0 $$

### Bogoliubov combination keeps the canonical commutator
Appears in: P-009

A linear combination of $\hat a$ and $\hat a^\dagger$ is again a canonical annihilator exactly when $|u|^2-|v|^2=1$.

$$ \hat b=u\hat a+v\hat a^\dagger,\qquad [\hat b,\hat b^\dagger]=(|u|^2-|v|^2)\hat I=\hat I \text{ when } |u|^2-|v|^2=1 $$

### Quadratic bosonic generators as sp(2N,R) matrices
Appears in: P-035

A quadratic bosonic generator corresponds to a block matrix with anti-Hermitian number block and symmetric pairing block.

$$ L=\begin{pmatrix}C&A\\\bar A&\bar C\end{pmatrix},\qquad C^\dagger=-C,\quad A^T=A $$

### Bosonic quadratic commutator table (E, F, G)
Appears in: P-035

The bilinears $\hat E,\hat F,\hat G$ close under commutators with the stated structure constants.

$$ [\hat E_{ij},\hat E_{kl}]=\delta_{jk}\hat E_{il}-\delta_{il}\hat E_{kj},\ [\hat E_{ij},\hat F_{kl}]=-\delta_{ik}\hat F_{jl}-\delta_{il}\hat F_{kj},\ [\hat E_{ij},\hat G_{kl}]=\delta_{jk}\hat G_{il}+\delta_{jl}\hat G_{ki},\ [\hat F_{ij},\hat G_{kl}]=\delta_{jk}\hat E_{li}+\delta_{jl}\hat E_{ki}+\delta_{ik}\hat E_{lj}+\delta_{il}\hat E_{kj} $$

Conditions: $\hat E_{ij}=\hat a_i^\dagger\hat a_j+\delta_{ij}/2$, $\hat F_{ij}=\hat a_i\hat a_j$, $\hat G_{ij}=\hat a_i^\dagger\hat a_j^\dagger$; $[\hat F,\hat F]=[\hat G,\hat G]=0$.

### Fractional-linear action on symmetric Gaussian parameters
Appears in: P-038

Canonical maps act on the symmetric Gaussian parameter by a matrix fractional-linear map, and the analytic representation carries a determinant multiplier.

$$ Z\mapsto (U+Z\bar V)^{-1}(V+Z\bar U),\qquad [\hat T^{(+)}(g)f](Z)=\det(U+Z\bar V)^{-1/2}f\left((U+Z\bar V)^{-1}(V+Z\bar U)\right) $$

### Siegel upper half-space and symplectic fractional-linear action
Appears in: P2-002

Real symplectic matrices act on the Siegel upper half-space by a fractional-linear map, and their blocks obey three quadratic relations.

$$ \mathfrak H_p=\{Z=X+iY:Z^T=Z,\ Y>0\},\quad Z\mapsto(AZ+B)(CZ+D)^{-1},\quad A^TC=C^TA,\ B^TD=D^TB,\ A^TD-C^TB=I_p $$

### Multimode particle-antiparticle bosonic blocks
Appears in: P2-020

A canonical map mixing particles and antiparticles preserves the indefinite form $E$, which fixes the block relations and the inverse.

$$ MEM^\dagger=E,\quad AA^\dagger-BB^\dagger=I_N,\quad DD^\dagger-CC^\dagger=I_N,\quad AC^\dagger=BD^\dagger,\quad M^{-1}=EM^\dagger E $$

### Symplectic covariance of Weyl symbols
Appears in: C-015

A metaplectic operator transforms Weyl symbols and Wigner functions by composing them with the classical linear canonical map.

$$ \hat R(F)^*\hat A\hat R(F)=\widehat{A\circ F},\qquad W_{\hat R(F)\psi,\hat R(F)\eta}(z)=W_{\psi,\eta}(F^{-1}z) $$

### Metaplectic generator from a symplectic path
Appears in: C-023

The quadratic generator of a symplectic path is read off as -J F' F^{-1}, and implementers compose only up to a sign (the metaplectic double cover).

$$ S_t=-J\dot F_tF_t^{-1},\qquad \hat R(F_1)\hat R(F_2)=\pm\hat R(F_1F_2) $$

### Bogoliubov blocks of a symplectic matrix
Appears in: C-024

The Heisenberg action of a quadratic propagator on the annihilator is a Bogoliubov mixing u a + v a^dag whose blocks are fixed combinations of the classical flow blocks.

$$ \hat U^*\hat a\hat U=u\hat a+v\hat a^\dagger,\quad u=\tfrac12[A+D+i(C-B)],\quad v=\tfrac12[A-D+i(C+B)],\quad uu^\dagger-vv^\dagger=I_n,\ uv^T=vu^T $$

## Fermionic canonical transformations

### Fermion-pair realization of su(2)
Appears in: X-004 (marked by the atlas as derived); P2-025; C-087

Creating or destroying a fermion pair in two modes, together with half the shifted number, closes on the compact su(2) algebra.

$$ \hat J_+=\hat c_1^\dagger \hat c_2^\dagger,\quad \hat J_-=\hat c_2\hat c_1=\hat J_+^\dagger,\quad \hat J_0=(\hat N_1+\hat N_2-\hat I)/2 $$

**Also stated in P2-025 (Single fermionic pair as spin 1/2):**

One fermionic particle-antiparticle pair is a spin-1/2 problem: the pair operators act as $\hat J_\pm,\hat J_0$ on the vacuum and the occupied pair.

$$ \hat J_+=\hat a^\dagger\hat b^\dagger,\ \hat J_-=\hat b\hat a,\ \hat J_0=\tfrac12(\hat a^\dagger\hat a+\hat b^\dagger\hat b-\hat I),\qquad |\zeta\rangle=(1+|\zeta|^2)^{-1/2}(\hat I+\zeta\hat a^\dagger\hat b^\dagger)|0\rangle $$

**Also stated in C-087 (Fermion pair rotation as SU(2)):**

Fermion pair creation and annihilation close an su(2), so a pair rotation sends the vacuum to cos r times the vacuum plus sin r times a filled pair.

$$ \hat J_+=\hat c_1^\dagger\hat c_2^\dagger,\ \hat J_-=\hat c_2\hat c_1,\ \hat J_3=\tfrac12(\hat N_1+\hat N_2-1),\qquad \hat U\Omega=\cos r\,\Omega+\frac\xi r\sin r\,\hat c_1^\dagger\hat c_2^\dagger\Omega $$

### Fermion-pair rotation of the vacuum
Appears in: X-004

The pair rotation turns the empty state into a superposition of empty and doubly occupied, with $\cos r$ and $\sin r$ amplitudes.

$$ e^{\xi \hat c_1^\dagger \hat c_2^\dagger-\bar\xi \hat c_2\hat c_1}|00\rangle=\cos r\,|00\rangle+e^{i\phi}\sin r\,|11\rangle $$

Conditions: $|11\rangle=\hat c_1^\dagger\hat c_2^\dagger|00\rangle$ (this ordering fixes the sign of the second term). It holds for all r; no chart condition is needed for the state.

### Norm bound for fermionic fields
Appears in: B-005

A smeared fermionic field is a bounded operator whose norm is at most the norm of its test function; in fact the bound is attained.

$$ \|\hat a(f)\|,\ \|\hat a^*(f)\|\le\|f\| $$

### Fermionic Gaussian implementer kernel
Appears in: B-060

The fermionic implementer's kernel is a Grassmann Gaussian with the same block structure as the bosonic one, but with prefactor det(Phi Phi*)^{+1/4}.

$$ K_{\hat U}=\theta\det(\hat\Phi\hat\Phi^*)^{1/4} \exp\!\left[-\tfrac12\left(a\bar{\hat\Psi}\hat\Phi^{-1}a-2a^*\hat\Phi^{-1}a+a^*\hat\Phi^{-1}\hat\Psi a^*\right)\right], \qquad N_{\hat U}=K_{\hat U}e^{-a^*a} $$

Conditions: two-sided Phi^{-1}; Grassmann order as displayed.

### Fermionic transformation of normal symbols with mixing
Appears in: B-063

For a mixing fermionic transformation the conjugated normal symbol is a Grassmann Gaussian integral of the original one, with the displayed matrix and determinant factor.

$$ \widehat{\mathcal L}=\begin{pmatrix}\bar{\hat\Phi}\bar{\hat\Psi}^{-1}&-\hat I\\\hat I&-\hat\Phi\hat\Psi^{-1}\end{pmatrix},\qquad \det(\hat P_n\hat\Psi\hat\Psi^*\hat P_n)^{1/2} $$

### Fermionic quadratic propagator
Appears in: B-074

For fermions the propagator kernel is det(Phi e^{itC})^{+1/2} times a Grassmann Gaussian built from the same classical generator.

$$ K_{e^{it\hat H}}=\det(\hat\Phi e^{it\hat C})^{1/2} \exp\!\left[-\tfrac12(a\ a^*) \begin{pmatrix}\bar{\hat\Psi}\hat\Phi^{-1}&\hat\Phi'^{-1}\\-\hat\Phi^{-1}&\hat\Phi^{-1}\hat\Psi\end{pmatrix}\binom a{a^*}\right] $$

### Fermionic diagonalization conditions
Appears in: B-086

A fermionic quadratic Hamiltonian is diagonalized by a proper transformation when S and S^2 - [S, T] - T^2 are bounded below by a positive constant.

$$ \hat S>\mu \hat I>0,\qquad \hat S^2-[\hat S,\hat T]-\hat T^2>\mu \hat I>0 $$

### Fermionic Bogoliubov constraints
Appears in: P-041; B-056 (also B-057)

A fermionic canonical map keeps the anticommutators exactly when the doubled coefficient matrix is unitary, and $Z=U^{-1}V$ is then skew-symmetric.

$$ \{\hat c_i,\hat c_j^\dagger\}=\delta_{ij},\ \hat{\mathbf b}=U\hat{\mathbf c}+V\hat{\mathbf c}^\dagger:\ UU^\dagger+VV^\dagger=I_N,\ UV^T+VU^T=0,\ M^\dagger M=I_{2N},\ Z=U^{-1}V=-Z^T $$

**Also stated in B-056 (Fermionic canonical relations):**

A fermionic canonical transformation has blocks with Phi Phi* + Psi Psi* = I and Phi Psi' + Psi Phi' = 0, the plus signs replacing the bosonic minus signs.

$$ \hat b=\hat\Phi\hat a+\hat\Psi\hat a^*,\quad \hat\Phi\hat\Phi^*+\hat\Psi\hat\Psi^*=\hat I,\quad\hat\Phi\hat\Psi'+\hat\Psi\hat\Phi'=0,\qquad \hat S^{-1}=\begin{pmatrix}\hat\Phi^*&\hat\Psi'\\\hat\Psi^*&\hat\Phi'\end{pmatrix} $$

### Fermionic quadratic commutator table
Appears in: P-042

Fermionic bilinears close under commutators with structure constants that differ from the bosonic table in the $[\hat F,\hat G]$ signs.

$$ [\hat F_{ij},\hat G_{kl}]=-\delta_{jk}\hat E_{li}+\delta_{ik}\hat E_{lj}+\delta_{jl}\hat E_{ki}-\delta_{il}\hat E_{kj},\quad [\hat F,\hat F]=[\hat G,\hat G]=0,\ \text{first three bosonic relations unchanged} $$

Conditions: $\hat E_{ij}=\hat c_i^\dagger\hat c_j-\delta_{ij}/2$, $\hat F_{ij}=\hat c_i\hat c_j$, $\hat G_{ij}=\hat c_i^\dagger\hat c_j^\dagger$.

### Fermionic analytic group action
Appears in: P-044

Fermionic canonical maps act on skew $Z$ by the same fractional-linear map as bosons, with determinant exponent $+1/2$.

$$ [\hat T^{(+)}(g)f](Z)=\det(U+Z\bar V)^{1/2}f\left((U+Z\bar V)^{-1}(V+Z\bar U)\right) $$

### Multimode fermionic blocks
Appears in: P2-026

The doubled fermionic coefficient matrix is unitary, which fixes the block relations.

$$ MM^\dagger=I_{2N},\quad AA^\dagger+BB^\dagger=I_N,\quad DD^\dagger+CC^\dagger=I_N,\quad AC^\dagger=-BD^\dagger $$

### CAR and finite fermionic Fock space
Appears in: C-064

Fermionic modes satisfy canonical anticommutation relations, so n modes span a 2^n-dimensional Fock space and a smeared annihilator has operator norm ||f||.

$$ \{\hat c_j,\hat c_k^\dagger\}=\delta_{jk}\hat I,\ \{\hat c_j,\hat c_k\}=0,\ \|\hat c(f)\|=\|f\| $$

### Linear fermionic exponential with complex coefficients
Appears in: C-065

A linear combination of fermionic creators and annihilators with complex coefficients squares to the scalar u.v, so its exponential is cosh plus a sinh-over-root multiple of the generator.

$$ e^{\hat L}=\cosh\sqrt{u\cdot v}\,\hat I+\frac{\sinh\sqrt{u\cdot v}}{\sqrt{u\cdot v}}(u\cdot\hat c^\dagger+v\cdot\hat c),\qquad \hat L^2=(u\cdot v)\hat I $$

### BdG Hamiltonian and its canonical matrix
Appears in: C-080

A general quadratic fermion Hamiltonian is encoded by a 2n-by-2n BdG matrix Sigma, whose exponential gives the Heisenberg action on (c, c^dag).

$$ \hat H=\tfrac12(\hat c^\dagger M\hat c+\hat c^\dagger L\hat c^\dagger-\hat c\bar L\hat c)-\tfrac14(\operatorname{Tr}M)\hat I,\quad \Sigma=\begin{pmatrix}M/2&L\\-\bar L&-\bar M/2\end{pmatrix},\quad e^{-it\hat H}\binom{\hat c}{\hat c^\dagger}e^{it\hat H}=e^{it\Sigma}\binom{\hat c}{\hat c^\dagger} $$

### Fermionic pair vacuum
Appears in: C-082

A quadratic fermion propagator maps the vacuum to a pair Gaussian with antisymmetric Gamma_t and amplitude det(A_t)^{1/2}, and conjugates translations by the inverse canonical matrix.

$$ \hat U_t\Omega=\det(A_t)^{1/2}\exp\!\left[\tfrac12\hat c^\dagger\Gamma_t\hat c^\dagger\right]\Omega,\quad \Gamma_t=B_t^T(A_t^T)^{-1}=-\Gamma_t^T;\qquad \hat U_t\hat T(\alpha)\hat U_t^{-1}=\hat T(F_t^{-1}\alpha) $$

### Fermionic Riccati and vacuum amplitude
Appears in: C-083

The fermionic pairing matrix obeys a Riccati equation and the vacuum amplitude a linear equation with a trace, solved by det(A_t)^{1/2}.

$$ i\dot\Gamma=\tfrac12(M\Gamma+\Gamma\bar M)+L+\Gamma\bar L\Gamma,\quad i\dot s=-\tfrac12\operatorname{Tr}(\bar L\Gamma+M/2)s,\quad s=\det(A_t)^{1/2}=\langle\Omega|\hat U_t|\Omega\rangle $$

### Fermionic Mehler kernel
Appears in: C-084

The holomorphic kernel of a fermionic Gaussian unitary is the vacuum amplitude times a Grassmann Gaussian built from A^{-1} and B.

$$ K_{\hat U}(a^*,a)=a_0\exp\!\left[-\tfrac12\left(a^T\bar BA^{-1}a-2a^{*T}A^{-1}a+a^{*T}A^{-1}Ba^*\right)\right],\quad a_0=\langle0|\hat U|0\rangle $$

## Gaussian vacua and overlaps

### Bosonic multimode Gaussian vacuum
Appears in: X-006 (source records cited: P-036 to P-045, B-048, B-061); P-036; B-048

The multimode squeezed vacuum is a normalized pair-creation exponential on the vacuum, annihilated by the combinations $\hat a+Z\hat a^\dagger$.

$$ |Z\rangle_B=\det(I_N-ZZ^\dagger)^{1/4}e^{-\frac12\hat a^{\dagger T}Z\hat a^\dagger}|0\rangle,\qquad (\hat a+Z\hat a^\dagger)|Z\rangle=0 $$

Conditions: $Z^T=Z$, $I_N-ZZ^\dagger>0$.

**Also stated in P-036 (Bosonic squeezed vacuum chart and its annihilator):**

A general bosonic group operator factors into pair-creation, number-conserving and pair-annihilation exponentials; on the vacuum this leaves a pair-creation Gaussian, which is annihilated by $\hat{\mathbf a}+Z\hat{\mathbf a}^\dagger$.

$$ \hat T(g)=\mathcal N\exp[-\tfrac12\xi_{ij}\hat X^{ij}]\exp[\alpha^k{}_l\hat X^l{}_k]\exp[-\tfrac12\eta^{ij}\hat X_{ij}];\qquad |Z\rangle=\mathcal N_B(Z)e^{-\frac12\hat{\mathbf a}^{\dagger T}Z\hat{\mathbf a}^\dagger}|0\rangle,\quad(\hat{\mathbf a}+Z\hat{\mathbf a}^\dagger)|Z\rangle=0 $$

**Also stated in B-048 (Transformed vacuum (homogeneous part)):**

The vacuum of the transformed annihilators is a Gaussian in a* with pairing matrix Phi^{-1} Psi, normalized by det(Phi Phi*)^{-1/4} up to a phase.

$$ F_0(a^*)=c\exp\!\left[-\tfrac12a^*\hat\Phi^{-1}\hat\Psi a^*-a^*\hat\Phi^{-1}f\right],\qquad c=\theta\det(\hat\Phi\hat\Phi^*)^{-1/4}\ (f=0) $$

### Fermionic multimode Gaussian vacuum
Appears in: X-006; B-061

The fermionic analogue uses an antisymmetric $Z$ and a determinant with a plus sign, raised to the power $-1/4$.

$$ |Z\rangle_F=\det(I_N+ZZ^\dagger)^{-1/4}e^{-\frac12\hat c^{\dagger T}Z\hat c^\dagger}|0\rangle,\qquad (\hat c+Z\hat c^\dagger)|Z\rangle=0 $$

Conditions: $Z^T=-Z$. Exists on the chart where the canonical block $u$ is invertible (the atlas's "$U$" here is undefined in X-006; see Errors found).

**Also stated in B-061 (Fermionic paired vacuum):**

The fermionic transformed vacuum is a pair exponential with skew pairing matrix Phi^{-1} Psi, normalized by det(Phi Phi*)^{1/4}.

$$ F(a^*)=\theta\det(\hat\Phi\hat\Phi^*)^{1/4} \exp[-\tfrac12a^*\hat\Phi^{-1}\hat\Psi a^*] $$

### Transformed vacuum: affine normalization (4.25)
Appears in: B-048, B-049 ((4.26'))

With a shift f, the vacuum amplitude acquires a real Gaussian factor in f, so its modulus is det(Phi Phi*)^{-1/4} times the displayed exponential.

$$ c=\frac{\theta}{\det(\hat\Phi\hat\Phi^*)^{1/4}} \exp\!\left[\frac14(f\ f^*) \begin{pmatrix}\hat\Phi'^{-1}\hat\Psi^*&-\hat I\\-\hat I&\hat\Phi^{*-1}\hat\Psi'\end{pmatrix} \binom f{f^*}\right] $$

### Bosonic Gaussian normalization and overlap
Appears in: P-037

A pair-creation Gaussian is normalized by $\det(I-ZZ^\dagger)^{1/4}$, and two of them overlap by the normalizations over $\det(I-WZ^\dagger)^{1/2}$.

$$ \mathcal N_B(Z)=\det(I_N-ZZ^\dagger)^{1/4},\qquad \langle Z|W\rangle=\frac{[\det(I_N-ZZ^\dagger)\det(I_N-WW^\dagger)]^{1/4}}{\det(I_N-WZ^\dagger)^{1/2}} $$

Conditions: $Z^T=Z$, $I_N-ZZ^\dagger>0$; square-root branch continued from $Z=W=0$.

### Bosonic Gaussian bilinear symbols and purity constraints
Appears in: P-039

On a pure bosonic Gaussian the pairing matrix $A$ and the number matrix $B$ are rational functions of $Z$ and obey quadratic purity constraints.

$$ A=-\hbar Z(I_N-Z^\dagger Z)^{-1},\quad B=\tfrac\hbar2 I_N+\hbar Z(I_N-Z^\dagger Z)^{-1}Z^\dagger,\quad B^2-A\bar A=\tfrac{\hbar^2}{4}I_N,\quad BA=A\bar B,\quad B^\dagger=B,\quad A^T=A $$

### Fermionic pair-exponential state, normalization and overlap
Appears in: P-043

A fermionic pair-exponential state with skew $Z$ is normalized by $\det(I+ZZ^\dagger)^{-1/4}$, is annihilated by $\hat{\mathbf c}+Z\hat{\mathbf c}^\dagger$, and two such states overlap by $\det(I+WZ^\dagger)^{1/2}$ over the normalizations.

$$ |Z\rangle=\det(I_N+ZZ^\dagger)^{-1/4}e^{-\frac12\hat{\mathbf c}^{\dagger T}Z\hat{\mathbf c}^\dagger}|0\rangle,\ Z^T=-Z,\ (\hat{\mathbf c}+Z\hat{\mathbf c}^\dagger)|Z\rangle=0,\quad \langle Z|W\rangle=\frac{\det(I_N+WZ^\dagger)^{1/2}}{[\det(I_N+ZZ^\dagger)\det(I_N+WW^\dagger)]^{1/4}} $$

### Fermionic bilinear symbols and purity constraints
Appears in: P-045

On a pure fermionic Gaussian the skew pairing matrix $A$ and the number matrix $B$ are rational in $Z$ and obey the same quadratic constraints as for bosons with the signs changed.

$$ A=-\hbar Z(I_N+Z^\dagger Z)^{-1},\quad B=\hbar Z(I_N+Z^\dagger Z)^{-1}Z^\dagger-\tfrac\hbar2I_N,\quad B^2-A\bar A=\tfrac{\hbar^2}{4}I_N,\quad BA=A\bar B,\quad B^\dagger=B,\quad A^T=-A $$

### Affine Gaussian state normalization (derived correction)
Appears in: P2-017

The state $e^{\zeta\hat a^{\dagger2}/2+\alpha\hat a^\dagger}|0\rangle$ is normalized by $(1-|\zeta|^2)^{1/4}$ times a Gaussian in $\alpha$ that includes the cross term $\mathrm{Re}(\bar\zeta\alpha^2)$.

$$ |\zeta,\alpha\rangle=Ne^{\zeta\hat a^{\dagger2}/2+\alpha\hat a^\dagger}|0\rangle,\qquad N=(1-|\zeta|^2)^{1/4}\exp\left[-\frac{|\alpha|^2+\operatorname{Re}(\bar\zeta\alpha^2)}{2(1-|\zeta|^2)}\right] $$

### Fermionic paired vacuum normalization
Appears in: P2-027

The fermionic paired vacuum $e^{\sum Z_{ij}\hat a_i^\dagger\hat b_j^\dagger}|0,0\rangle$ has squared norm $\det(I_N+Z^\dagger Z)$.

$$ |Z\rangle=\det(I_N+Z^\dagger Z)^{-1/2}e^{\sum_{ij}Z_{ij}\hat a_i^\dagger \hat b_j^\dagger}|0,0\rangle $$

### Fermionic paired vacuum overlap
Appears in: P2-028

Two fermionic paired vacua overlap by $\det(I_N+Z_1^\dagger Z_2)$ over the normalizations.

$$ \langle Z_1|Z_2\rangle=\frac{\det(I_N+Z_1^\dagger Z_2)}{\det(I_N+Z_1^\dagger Z_1)^{1/2}\det(I_N+Z_2^\dagger Z_2)^{1/2}},\qquad Z'=(AZ+B)(CZ+D)^{-1} $$

### Fermionic pair-number generating determinant
Appears in: P2-029

The fermionic pair-number generating function is $\det(I+tZ^\dagger Z)/\det(I+Z^\dagger Z)$, a product of Bernoulli laws.

$$ \mathcal F(t)=\sum_{n=0}^NW_nt^n=\frac{\det(I_N+tZ^\dagger Z)}{\det(I_N+Z^\dagger Z)},\qquad p_j=\frac{\lambda_j}{1+\lambda_j} $$

### Squeezed vacuum from the Bargmann Riccati equation
Appears in: C-027

The squeezed vacuum is the Gaussian e^{a^dag W a^dag/2} on the vacuum with normalization det(I - W^*W)^{1/4}, obtained from a matrix Riccati equation in Bargmann space; its vacuum amplitude is det(cosh|B|)^{-1/2}.

$$ \dot M=B-MB^*M,\ M_0=0,\quad \dot a=-\tfrac12a\operatorname{Tr}(B^*M),\quad \hat S(B)|0\rangle=\det(I_n-W^*W)^{1/4}\exp\!\left(\tfrac12\hat a^\dagger\cdot W\hat a^\dagger\right)|0\rangle,\quad \langle0|\hat S(B)|0\rangle=\det(\cosh|B|)^{-1/2} $$

### Coordinate form of the squeezed vacuum
Appears in: C-028; one-mode case in C-049 (delta = (1-omega)/(1+omega))

In position space the squeezed vacuum is a normalized complex Gaussian whose covariance is the Cayley transform of W.

$$ \Gamma=i(I_n-W)(I_n+W)^{-1},\qquad \psi_B(x)=(\pi\hbar)^{-n/4}\det(I_n-W^*W)^{1/4}\det(I_n+W)^{-1/2}e^{ix^T\Gamma x/(2\hbar)} $$

### Reduced Gaussian channel from a joint quadratic evolution
Appears in: C-056

A joint quadratic evolution with a Gaussian environment acts on the system by a linear map of the characteristic-function argument times a Gaussian noise factor.

$$ \chi_{S,\rm out}(k)=\chi_{S,\rm in}(X^Tk)\exp[-\tfrac12k^TYV_EY^Tk],\qquad V_{S,\rm out}=XV_SX^T+YV_EY^T $$

## Quadratic propagators and Riccati

### General quadratic expression
Appears in: B-065

The general quadratic Hamiltonian has pair-creation, pair-annihilation and number-type kernels plus linear and constant terms; only the symmetric (bosons) or antisymmetric (fermions) part of the pair kernels contributes.

$$ \hat H=\tfrac12\int[A(x,y)\hat a^*(x)\hat a^*(y)+B(x,y)\hat a(x)\hat a(y)+2C(x,y)\hat a^*(x)\hat a(y)]\,dxdy +\hat a^*(f)+\hat a(g)+k\hat I $$

### Classical generator and Heisenberg action of a quadratic Hamiltonian
Appears in: B-069, B-071, B-072

A quadratic Hamiltonian evolves the fields linearly: conjugation by e^{itH} acts on (a, a*) by the matrix exponential of the classical generator L, plus a constant shift g from the linear terms.

$$ \widehat{\mathcal L}=\begin{pmatrix}-\hat C&-\hat A\\\bar{\hat A}&\bar{\hat C}\end{pmatrix},\quad e^{it\widehat{\mathcal L}}=\begin{pmatrix}\hat\Phi&\hat\Psi\\\bar{\hat\Psi}&\bar{\hat\Phi}\end{pmatrix},\quad \hat U(t)\binom{\hat a}{\hat a^*}\hat U(t)^{-1} =e^{it\widehat{\mathcal L}}\binom{\hat a}{\hat a^*}+\binom g{g^*}\hat I,\quad \binom g{g^*}=i\int_0^t e^{i\tau\widehat{\mathcal L}}\,d\tau\binom{-f}{f^*} $$

Conditions: U(t) = e^{itH}; the Hilbert-Schmidt criterion on F(t) in B-069 is a source theorem.

### Quadratic propagator vacuum coefficient (corrected (6.23))
Appears in: B-071 (B-077 cites the same scalar terms)

The vacuum amplitude of e^{itH} is det(Phi e^{itC})^{-1/2} times the exponential of a time integral of terms quadratic and linear in the shift g, the quadratic one with coefficient 1/2.

$$ c(t)=\det(\hat\Phi e^{it\hat C})^{-1/2} \exp\!\left[i\int_0^t\left(\tfrac12g\hat\Phi'^{-1}\bar{\hat A}\hat\Phi^{-1}g-f^*\hat\Phi^{-1}g\right)d\tau\right] $$

Conditions: square root continued from c(0) = 1; constant k adds e^{itk}.

### Logarithmic derivative of a determinant
Appears in: B-073

The time derivative of log det M is the trace of M^{-1} times the derivative of M (Jacobi's formula).

$$ \frac{d}{dt}\log\det \hat M=\operatorname{Tr}(\hat M^{-1}\dot{\hat M}) $$

### Riccati fixed point for the diagonalizing transformation
Appears in: B-085

The ground state of H = (q beta q + p gamma p)/2 is the Gaussian exp(i q Z q/2), with Z = iW and W the positive solution of W gamma W = beta, given in closed form by the displayed square roots.

$$ \hat Z=i\hat\gamma^{-1/2}(\hat\gamma^{1/2}\hat\beta\hat\gamma^{1/2})^{1/2}\hat\gamma^{-1/2},\qquad \hat\beta=\hat S+\hat T,\quad\hat\gamma=\hat S-\hat T $$

### Riccati equation and phase for squeezing dynamics (phase sign corrected)
Appears in: P2-013

Under a Hamiltonian linear in the su(1,1) generators, a coherent state stays coherent; its label obeys a Riccati equation and its phase obeys a first-order equation.

$$ \hat H=i(\beta \hat K_+-\bar\beta \hat K_--i\gamma \hat K_0),\ |\psi\rangle=e^{-i\phi}|\zeta;k\rangle:\quad \dot\zeta=\beta-\bar\beta\zeta^2-i\gamma\zeta,\qquad i\dot\phi=k(\bar\beta\zeta-\beta\bar\zeta+i\gamma) $$

### Floquet classification of a quadratic oscillator
Appears in: P2-014

The sign of the Lorentz norm of the generator classifies a periodic quadratic oscillator: discrete spectrum, two-sided continuum, or one-sided continuum.

$$ \Omega_0^2-\Omega_1^2-\Omega_2^2>0:\ \text{spectrum } \pm\hbar\Omega(k+n);\quad <0:\ \text{two-sided continuous};\quad =0,\ \Omega\ne0:\ \text{one-sided continuous} $$

### Multimode pair Riccati equation and vacuum phase
Appears in: P2-024

The pair matrix of a multimode paired vacuum obeys a matrix Riccati equation, and its vacuum phase obeys a trace equation.

$$ i\dot Z=A+ZA^\dagger Z+CZ+ZD^T,\qquad \dot\phi=\tfrac12\operatorname{Tr}(A^\dagger Z+Z^\dagger A) $$

Conditions: $\hat H=A_{ij}\hat a_i^\dagger\hat b_j^\dagger+\bar A_{ij}\hat a_i\hat b_j+C_{ij}\hat a_i^\dagger\hat a_j+D_{ij}\hat b_i^\dagger\hat b_j$.

### Mehler kernel and Euclidean oscillator kernel
Appears in: C-012

The sum of products of Hermite functions weighted by w^k has a closed Gaussian form; with w = e^{-tau} and the zero-point factor it becomes the imaginary-time oscillator kernel.

$$ \sum_{k}w^{|k|}\phi_k(x)\phi_k(y)=\pi^{-n/2}(1-w^2)^{-n/2}\exp\!\left[-\frac{1+w^2}{2(1-w^2)}(|x|^2+|y|^2)+\frac{2w\,x\cdot y}{1-w^2}\right],\quad |w|<1 $$

and, n = 1, kernel of e^{-tau H_os} = (2 pi sinh tau)^{-1/2} exp[-(x^2+y^2)/(2 tanh tau) + xy/sinh tau].

### Linear quadratic flow and Heisenberg action
Appears in: C-018 (constant-S case also in C-023, C-030)

For a quadratic Hamiltonian the Heisenberg evolution of (Q, P) is linear and given by the classical flow matrix solving F' = J S F, which is e^{tJS} for constant S.

$$ \dot F_t=JS_tF_t,\ F_0=I_{2n},\quad \hat U_t^*(\hat Q,\hat P)^T\hat U_t=F_t(\hat Q,\hat P)^T,\quad F_t=e^{tJS}\ (S\ \text{constant}) $$

Conditions: H = z^T S z / 2 quantized symmetrically with the factor 1/2.

### Gaussian propagation by a linear-fractional map
Appears in: C-019

A Gaussian e^{ix^T Gamma x/(2 hbar)} stays Gaussian under a quadratic propagator: its covariance moves by the linear-fractional map of the flow blocks and its amplitude picks up det(A_t + B_t Gamma)^{-1/2}.

$$ \Gamma_t=(C_t+D_t\Gamma)(A_t+B_t\Gamma)^{-1},\qquad \hat U_t\big[a_\Gamma e^{ix^T\Gamma x/(2\hbar)}\big]=a_\Gamma\det(A_t+B_t\Gamma)^{-1/2}e^{ix^T\Gamma_tx/(2\hbar)} $$

Conditions: Gamma symmetric, Im Gamma > 0; root continued from t = 0.

### Siegel-domain invariance and determinant cocycle
Appears in: C-020

The linear-fractional action of a symplectic matrix keeps a complex symmetric matrix with positive imaginary part inside that domain, composes as a group action, and its determinant factor is a cocycle.

$$ \Im\Gamma_F=(A+B\Gamma)^{-*}(\Im\Gamma)(A+B\Gamma)^{-1},\quad \alpha(F_1F_2)=\alpha(F_1)\alpha(F_2),\quad \delta(F_1F_2,\Gamma)=\delta(F_1,\alpha(F_2)\Gamma)\delta(F_2,\Gamma) $$

### Matrix Riccati and amplitude equations
Appears in: C-021

The complex covariance of a propagated Gaussian obeys a matrix Riccati equation and its amplitude a linear equation with a trace.

$$ \dot\Gamma=-K-L\Gamma-\Gamma L^T-\Gamma G\Gamma,\qquad \dot a=-\tfrac12\operatorname{Tr}(L+G\Gamma)a,\qquad S=\begin{pmatrix}K&L\\L^T&G\end{pmatrix} $$

### Coherent matrix elements of a metaplectic operator
Appears in: C-022

The coherent-state matrix element of the propagator of a quadratic Hamiltonian is a closed Gaussian in the two labels, with det(Y)^{-1/2}, Y = A + D + i(B - C), as prefactor.

$$ \langle\phi_z,\hat U\phi_{z'}\rangle=2^{n/2}\det(Y)^{-1/2}e^{-i\sigma(Fz',z)/2}e^{-(|x|^2+|\xi|^2)/4}e^{-(\xi+ix)^TW(\xi+ix)/4},\quad Y=A+D+i(B-C),\ W=ZY^{-1},\ (x,\xi)=z-Fz' $$

Conditions: hbar = 1, Z = A - D + i(B+C), root along the path.

### Squeeze-rotation factorization and its vacuum phase
Appears in: C-029

A quadratic propagator factors into a squeeze and a number-preserving rotation; the rotation multiplies the vacuum by a unit-modulus phase fixed by det(Y_t) and det(u_t u_t^dag).

$$ \hat U_t=\hat S(B_t)\hat R_t,\quad \hat R_t|0\rangle=\vartheta_t|0\rangle,\quad \vartheta_t=\frac{2^{n/2}\det(Y_t)^{-1/2}}{\det(u_tu_t^\dagger)^{-1/4}},\quad |\vartheta_t|=1,\quad Y_t=2\bar u_t $$

### Herman-Kluk representation, exact for quadratics
Appears in: C-054

For a quadratic Hamiltonian, the phase-space integral of |phi_{F_t z}><phi_z| with the prefactor det(u_t)^{1/2} reproduces the propagator exactly.

$$ \hat U_t=(2\pi\hbar)^{-n}\int C_t\,|\phi_{F_tz}\rangle\langle\phi_z|\,dz,\qquad C_t=\det\!\left[\frac{A_t+D_t+i(C_t^{\rm block}-B_t)}2\right]^{1/2} $$

## Forced and parametric oscillators

### Forced oscillator: completing the square
Appears in: B-081

A linearly forced oscillator is a free oscillator in the displaced field b = a + gf/omega, lowered by the constant g^2 |f|^2/omega.

$$ \hat H=\int\omega\hat b^*\hat b\,dk- g^2\int\frac{|f|^2}{\omega}\,dk\,\hat I,\qquad \hat a=\hat b-gf/\omega\,\hat I $$

### Regularized forced Hamiltonian
Appears in: B-083

Writing the forced Hamiltonian as the displaced number operator adds back the constant g^2 int |f|^2/omega, the term that diverges when only int |f|^2/omega^2 is finite.

$$ \hat H_{\mathrm{reg}}=\int\omega(k)\left(\hat a^*+g\frac{f^*}{\omega}\hat I\right) \left(\hat a+g\frac f\omega\hat I\right)dk $$

### Forced oscillator: displacement plus accumulated phase
Appears in: P2-010; G-022

A linearly forced oscillator stays coherent, the amplitude integrates the force, a phase accumulates as $\dot\phi=\mathrm{Im}(\bar\beta\alpha)$, and from the vacuum the number distribution is Poisson.

$$ |\psi(t)\rangle=e^{-i\phi(t)}|\alpha(t)\rangle,\quad \dot\alpha=\beta(t),\quad \dot\phi=\operatorname{Im}(\bar\beta\alpha),\quad W_n=e^{-|\gamma|^2}|\gamma|^{2n}/n! $$

Conditions: interaction picture, generator $\beta\hat a^\dagger-\bar\beta\hat a$ with no extra c-number term.

**Also stated in G-022 (Linearly driven evolution equals a displacement up to phase):**

A linearly driven oscillator, built as a product of short displacements, is one displacement by the accumulated amplitude times an accumulated phase.

$$ \hat U(T,0)=e^{i\Phi(T)}\hat D(\alpha(T)),\qquad\alpha(T)=\int_0^T u(t)\,dt. $$

### Parametric oscillator su(1,1) generators
Appears in: P2-012

The quadratic oscillator operators $\hat K_0$ and $\hat K_2$ realize su(1,1); $\hat K_0$ is half the shifted number operator.

$$ \hat K_0=\tfrac14(\hat p^2/\omega_0+\omega_0\hat q^2),\qquad \hat K_2=\tfrac14(\hat p\hat q+\hat q\hat p) $$

### Singular oscillator ground state and su(1,1) label
Appears in: P2-015

The half-line inverse-square oscillator has ground state $x^ae^{-x^2/2}$ with $a(a-1)=2g^2$, energy $a+1/2$, and an equally spaced ladder with su(1,1) label $k=E_0/2$.

$$ \psi_0\sim x^a e^{-x^2/2},\quad a(a-1)=2g^2,\quad E_0=a+\tfrac12,\quad k=\tfrac12(a+\tfrac12),\quad E_n=2n+E_0 $$

Conditions: $\hat H=\tfrac12(\hat p^2+x^2)+g^2/x^2$ on the half line, regular boundary condition.

### Forced parametric oscillator generating kernel and vacuum amplitude
Appears in: P2-018

For a forced parametric oscillator the vacuum amplitude has modulus $(1-\rho)^{1/4}e^{-\nu(1-\sqrt\rho\cos2\varphi)/2}$.

$$ |\langle0|\psi\rangle|=(1-\rho)^{1/4}e^{-\nu(1-\sqrt\rho\cos2\varphi)/2} $$

### Drive equation and displacement channel factor
Appears in: P2-042

A linear drive adds $-i\epsilon$ to the decay of the mean field, and the drive enters the state as a displacement.

$$ \hat\rho\mapsto \hat D(v)\hat\rho \hat D(v)^\dagger,\quad \hat D(v)=e^{v\hat a^\dagger-\bar v \hat a};\qquad \hat H_d=\epsilon\hat a^\dagger+\bar\epsilon\hat a\ \Rightarrow\ \dot\beta=-(\gamma/2)\beta-i\epsilon(t) $$

### Propagator of a driven quadratic Hamiltonian
Appears in: C-032

For a Hamiltonian that is quadratic plus linear, the propagator is the classical action phase times a translation along the classical orbit, the quadratic propagator, and the inverse translation at the start time.

$$ \delta_t=\int_0^t[p_s\dot q_s-H(z_s,s)]ds-\tfrac12(q_tp_t-qp),\qquad \hat U_2(t,s)=e^{i(\delta_t-\delta_s)/\hbar}\hat T(z_t)\hat S(B_t)\hat R_t\hat R_s^{-1}\hat S(-B_s)\hat T(-z_s) $$

### Periodic quadratic evolution and quasi-energies
Appears in: C-035

For a periodic quadratic Hamiltonian in a stable region, the one-period propagator has eigenvalues of the form e^{-i rho pi (k+1/2)}.

$$ \text{Floquet eigenvalues } e^{-i\rho\pi(k+1/2)} $$

### Time-dependent SU(1,1) Hamiltonian
Appears in: C-050

For a Hamiltonian linear in the su(1,1) generators with time-dependent coefficients, the propagator is a squeeze, a K0 phase rotation and an inverse squeeze, with the squeeze coordinate solving a scalar Riccati equation.

$$ i\dot\omega=\bar\lambda\omega^2+\mu\omega+\lambda,\ \dot\gamma=-\lambda\bar\omega-\bar\lambda\omega-\mu,\qquad \hat U(t,s)=\hat S(\beta_t)e^{i(\gamma_t-\gamma_s)\hat K_0}\hat S(-\beta_s),\ \omega_t=\frac{\beta_t}{|\beta_t|}\tanh|\beta_t| $$

### Degenerate parametric amplifier gives a squeeze of parameter 2 eta t
Appears in: G-066, G-067 (four-wave mixing, same exponent with $\eta=\chi^{(3)}\beta^2$)

Under an undepleted classical pump at resonance the degenerate amplifier generates a single-mode squeeze with parameter 2 eta t.

$$ \hat H_I=i\hbar(\eta^*\hat a^2-\eta \hat a^{\dagger2}),\quad\hat U_I(t)=e^{t(\eta^*\hat a^2-\eta \hat a^{\dagger2})}=\hat S(2\eta t). $$

Conditions: undepleted classical pump, resonance.

### Nondegenerate parametric amplifier gives a two-mode squeeze of parameter eta t
Appears in: G-075

Under an undepleted pump the nondegenerate amplifier generates a two-mode squeeze with parameter eta t, with no factor 2.

$$ \hat H_I=i\hbar(\eta^*\hat a\hat b-\eta \hat a^\dagger \hat b^\dagger),\quad \hat U_I=\hat S_2(\eta t). $$

### Driven damped oscillator: coherent amplitude and steady state
Appears in: G-100

A driven damped oscillator keeps coherent states coherent; the amplitude relaxes to -2iG/gamma.

$$ \dot\alpha=-iG-\gamma\alpha/2,\quad\alpha(t)=\alpha_0e^{-\gamma t/2}-\frac{2iG}{\gamma}(1-e^{-\gamma t/2}),\quad\alpha_{\rm ss}=-2iG/\gamma. $$

## Normal symbols and normal-order exponentials

### Number-conserving exponential in normal order
Appears in: X-007 (source records cited: B-013, B-028, B-040, B-107, G-009, G-097); B-107 (special cases in B-013, B-028, B-105); G-097

The exponential of a number-conserving quadratic equals the normal-ordered exponential of the same form with $h$ replaced by $e^h-I_N$.

$$ e^{\hat a^\dagger h \hat a}=:\!\exp[\hat a^\dagger(e^h-I_N)\hat a]\!: $$

Conditions: finite modes; the same identity for fermions with CAR normal ordering.

**Also stated in B-107 (Number-type exponential in normal order):**

The exponential of a quadratic number-conserving form a* h a is the normal-ordered exponential of a*(e^h - I)a.

$$ e^{\hat a^*\hat h\hat a}=:\!\exp[\hat a^*(e^{\hat h}-\hat I)\hat a]\!: $$

**Also stated in G-097 (Normal ordering of a number exponential and the vacuum projector):**

The exponential of lambda n is the normal-ordered exponential of (e^lambda - 1) n; at lambda -> -infinity this is the vacuum projector :e^{-n}:.

$$ e^{\lambda \hat n}={:\exp[(e^\lambda-1)\hat n]:},\qquad|0\rangle\langle0|={:e^{-\hat n}:}. $$

### Vacuum projector as a normal-ordered exponential
Appears in: X-007; B-013 (related: B-107 normal-order identity)

Normal ordering $e^{-\hat N}$ gives the projector onto the vacuum; for fermions the projector is the product of the single-mode empty projectors.

$$ |0\rangle\langle0|=:\!e^{-\sum_i\hat N_i}\!:\quad\text{(bosons)},\qquad |0\rangle\langle0|=\prod_i(\hat I-\hat N_i)\quad\text{(fermions)} $$

**Also stated in B-013 (Vacuum projector as a normal-ordered exponential):**

The projector onto the vacuum is the normal-ordered exponential of minus the number operator: the normal-ordered exponential :e^{c a*a}: multiplies |n> by (1 + c)^n, and c = -1 keeps only n = 0.

$$ \hat P_0=:\!e^{-\hat a^*\hat a}\!: $$

### Matrix kernel versus normal symbol
Appears in: B-014

The matrix kernel of an operator (its matrix element between unnormalized coherent vectors) equals its normal symbol times the coherent overlap e^{a*a}.

$$ N_{\hat A}(a^*,a)=K_{\hat A}(a^*,a)e^{-a^*a},\qquad K_{\hat A}=N_{\hat A}e^{a^*a} $$

### Symbols as coherent matrix elements
Appears in: B-019

The normal symbol on the diagonal a* = conj(a) is the expectation value of the operator in the normalized coherent state of amplitude a.

$$ N_{\hat A}(a^*,a)=\langle P_a,\hat A P_a\rangle $$

### Differential realization on normal symbols
Appears in: B-015

Multiplying an operator by a on the left, or by a* on the right, acts on its normal symbol as multiplication plus a derivative in the other variable, because only those two products need reordering.

$$ N_{\hat a\hat A}=(a+\partial_{a^*})N_{\hat A},\qquad N_{\hat A\hat a^*}=(a^*+\partial_{a})N_{\hat A} $$

Conditions: polynomial or analytic symbols, smearing implicit.

### Rank-one projector symbols
Appears in: B-020

The projector onto a normalized vector Phi has matrix kernel Phi(a*) Phi*(a) and normal symbol that times e^{-a*a}.

$$ K_{|\Phi\rangle\langle\Phi|}=\Phi(a^*)\Phi^*(a),\qquad N_{|\Phi\rangle\langle\Phi|}=\Phi(a^*)\Phi^*(a)e^{-a^*a} $$

### Kernel acting on a vector
Appears in: B-026

An operator acts on a vector functional by integrating its matrix kernel against the vector with the Bargmann weight e^{-alpha* alpha}.

$$ (\hat A\Phi)(a^*)=\int K_{\hat A}(a^*,\alpha)\Phi(\alpha^*)e^{-\alpha^*\alpha}\Pi d\alpha^*d\alpha $$

### Product of operators: matrix kernels
Appears in: B-027

The matrix kernel of a product is the Bargmann-weighted integral of the two kernels over the intermediate coherent label.

$$ K_{\hat A\hat B}(a^*,a)=\int K_{\hat A}(a^*,\alpha)K_{\hat B}(\alpha^*,a)e^{-\alpha^*\alpha}\Pi d\alpha^*d\alpha $$

### Product of operators: normal symbols
Appears in: B-027

The normal symbol of a product is a Gaussian convolution of the two normal symbols; for polynomials it equals the series N_A exp(left d/da, right d/da*) N_B.

$$ N_{\hat A\hat B}(a^*,a)=\int N_{\hat A}(a^*,\alpha)N_{\hat B}(\alpha^*,a)e^{-(\alpha^*-a^*)(\alpha-a)}\Pi d\alpha^*d\alpha $$

### Resolvent kernel of a number-type Hamiltonian
Appears in: B-075

The resolvent of H = a* C a is an integral over xi in (0, 1) of the kernel of xi^H weighted by xi^{-z-1}, valid only where that integral converges.

$$ K_{(\hat H-z\hat I)^{-1}}=\int_0^1\exp(a^*\xi^{\hat C}a)\,\xi^{-z-1}d\xi $$

Conditions: the atlas's added restriction: needs Re z below the spectrum of H (Re z < 0 for C >= 0), then analytic continuation.

### Covariant (Q) and contravariant (P) symbols from ordered expansions
Appears in: P-015

The covariant symbol is the coherent expectation value and is read from the normal-ordered expansion; the contravariant symbol is the weight in the coherent-projector expansion and is read from the antinormal-ordered expansion.

$$ Q_{\hat A}(\alpha)=\langle\alpha|\hat A|\alpha\rangle,\ \ \hat A=\sum A_{mn}(\hat a^\dagger)^m\hat a^n\Rightarrow Q_{\hat A}=\sum A_{mn}\bar\alpha^m\alpha^n;\quad \hat A=\int\frac{d^2\alpha}{\pi}P_{\hat A}|\alpha\rangle\langle\alpha|,\ \ \hat A=\sum \widetilde A_{mn}\hat a^n(\hat a^\dagger)^m\Rightarrow P_{\hat A}=\sum \widetilde A_{mn}\bar\alpha^m\alpha^n $$

### Off-diagonal covariant symbol by analytic continuation
Appears in: P2-005

The covariant symbol is the coherent expectation value divided by the overlap, and its off-diagonal version uses two different coherent labels.

$$ A(z,\bar z)=\frac{\langle z|\hat A|z\rangle}{\langle z|z\rangle},\qquad A(z,\bar w)=\frac{\langle z|\hat A|w\rangle}{\langle z|w\rangle} $$

### Disk star product to first order in h
Appears in: P2-006

On the disk, the Berezin transform is the identity plus $h$ times the invariant Laplacian, and the symbol product differs from the pointwise product at order $h$ by a derivative pairing.

$$ \hat T_hf=f+h\hat\Delta f+o(h),\ \hat\Delta=(1-|z|^2)^2\partial_z\partial_{\bar z},\qquad A_1\star A_2=A_1A_2+h(1-|z|^2)^2(\partial_{\bar z}A_1)(\partial_zA_2)+o(h) $$

### Normal-ordered powers as falling factorials
Appears in: G-006, G-078 ($\hat a^{\dagger2}\hat a^2=\hat n(\hat n-\hat I)$), G-097 (eigenvalue of ${:\hat n^r:}$)

Normal-ordered powers a^dagger^r a^r equal the falling factorial of the number operator.

$$ (\hat a^\dagger)^r \hat a^r=\hat n(\hat n-\hat I)\cdots(\hat n-r\hat I+\hat I) $$

### Normally ordered coherent matrix elements
Appears in: G-020

Between coherent states a normally ordered operator is replaced by its symbol with a^dagger -> beta^* and a -> alpha, times the overlap.

$$ \langle\beta|(\hat a^\dagger)^p \hat a^q|\alpha\rangle=(\beta^*)^p\alpha^q\langle\beta|\alpha\rangle,\quad\langle\beta|{:F(\hat a^\dagger,\hat a):}|\alpha\rangle=F(\beta^*,\alpha)\langle\beta|\alpha\rangle. $$

### Normally ordered correlation functions
Appears in: G-035

Ideal photodetection measures normally ordered field correlations; for one mode the normalized r-th order coherence is the r-th factorial moment over the mean to the r-th power.

$$ G^{(r)}=\langle \hat E^{(-)}(x_1)\cdots \hat E^{(-)}(x_r)\hat E^{(+)}(x_r)\cdots \hat E^{(+)}(x_1)\rangle,\qquad g^{(r)}=\langle(\hat a^\dagger)^r \hat a^r\rangle/\langle \hat n\rangle^r. $$

Conditions: $\langle\hat n\rangle\ne0$.

### Second-order coherence and Mandel Q
Appears in: G-036

g^{(2)} is the second factorial moment over the squared mean, equal to 1 plus the excess number variance; it is 1 for coherent light, 2 for thermal light and 1 - 1/m for a number state, and undefined for the vacuum.

$$ g^{(2)}(0)=\frac{\langle \hat n(\hat n-\hat I)\rangle}{\langle \hat n\rangle^2}=1+\frac{\operatorname{Var}\hat n-\langle \hat n\rangle}{\langle \hat n\rangle^2},\quad Q_M=\langle \hat n\rangle(g^{(2)}(0)-1);\quad g^{(2)}=1,\ 2,\ 1-1/m. $$

### Linear quadrature exponential: normal-order conversion factor
Appears in: G-087

The exponential of a quadrature equals its normal-ordered form times e^{y^2/8}, the Gaussian factor of the central commutator.

$$ e^{y\Delta \hat X}={:e^{y\Delta \hat X}:}\,e^{y^2/8}. $$

### Centered quadrature moments from normal moments
Appears in: G-088

Centered quadrature moments expand into normally ordered moments with Gaussian contraction weights N!/((N-2k)! k! 8^k); for a coherent state only the full contraction survives.

$$ \langle(\Delta \hat X)^N\rangle=\sum_{k=0}^{\lfloor N/2\rfloor}\frac{N!}{(N-2k)!k!8^k}\langle{: (\Delta \hat X)^{N-2k}:}\rangle,\qquad\langle(\Delta \hat X)^{2\ell}\rangle_{\rm coh}=(2\ell-1)!!/4^\ell. $$

## Weyl and antinormal symbols, characteristic functions

### Differential correspondences for multiplying a displacement
Appears in: P-011

Multiplying a displacement by a ladder operator on either side equals a first-order Wirtinger derivative of the displacement plus a linear term.

$$ \hat a^\dagger \hat D=(\partial_\alpha+\bar\alpha/2)\hat D,\quad \hat D\hat a^\dagger=(\partial_\alpha-\bar\alpha/2)\hat D,\quad \hat a\hat D=(-\partial_{\bar\alpha}+\alpha/2)\hat D,\quad \hat D\hat a=(-\partial_{\bar\alpha}-\alpha/2)\hat D $$

Conditions: $\alpha,\bar\alpha$ independent Wirtinger variables.

### Q and P symbols of a displacement
Appears in: P-017

The covariant symbol of a displacement is a Gaussian-damped plane wave in phase space, and its contravariant symbol has the opposite sign of the $|\gamma|^2/2$ term.

$$ Q_{\hat D(\gamma)}(\alpha)=e^{-|\gamma|^2/2+\gamma\bar\alpha-\bar\gamma\alpha},\qquad P_{\hat D(\gamma)}(\alpha)=e^{+|\gamma|^2/2+\gamma\bar\alpha-\bar\gamma\alpha} $$

### Distributional trace of displacement and operator reconstruction
Appears in: P-017, P-020 (the $\hat\rho=\int\chi_W\hat D(-\gamma)$ line)

A displacement has a delta-function trace, displacements are orthogonal in the Hilbert-Schmidt sense, and any operator is rebuilt from its displacement transform.

$$ \operatorname{Tr}\hat D(\gamma)=\pi\delta^{(2)}(\gamma),\quad \operatorname{Tr}[\hat D(\alpha)\hat D(\beta)^\dagger]=\pi\delta^{(2)}(\alpha-\beta),\quad \hat A=\int\frac{d^2\gamma}{\pi}\operatorname{Tr}[\hat D(\gamma)\hat A]\hat D(-\gamma) $$

Conditions: Hilbert-Schmidt / distributional sense; $d^2\gamma=d\,\mathrm{Re}\gamma\,d\,\mathrm{Im}\gamma$.

### Weyl operator from displaced parity, Weyl symbol and reconstruction
Appears in: P-019; G-034

The Weyl operator is twice a displaced parity, parity flips the sign of a displacement, the Weyl symbol is the trace against it, and these operators are delta-orthogonal.

$$ \hat T(\alpha)=2\hat D(\alpha)\hat\Pi \hat D(-\alpha)=2\hat D(2\alpha)\hat\Pi,\quad \hat\Pi \hat D(\alpha)\hat\Pi=\hat D(-\alpha),\quad W_{\hat A}(\alpha)=\operatorname{Tr}[\hat T(\alpha)\hat A],\quad \hat A=\int d\mu(\alpha)W_{\hat A}(\alpha)\hat T(\alpha),\quad \pi^{-1}\operatorname{Tr}[\hat T(\alpha)\hat T(\beta)]=\delta^{(2)}(\alpha-\beta) $$

Conditions: $d\mu=d^2\alpha/\pi$, $\hat\Pi=(-1)^{\hat a^\dagger\hat a}$.

**Also stated in G-034 (Wigner function as displaced parity):**

The Wigner function is 2/pi times the expectation of the parity operator displaced to alpha, and the displaced parity equals D(2 alpha) times the parity.

$$ W_{\hat\rho}(\alpha)=\frac2\pi\operatorname{Tr}[\hat\rho \hat D(\alpha)\hat\Pi \hat D^\dagger(\alpha)],\qquad\hat D(\alpha)\hat\Pi \hat D^\dagger(\alpha)=\hat D(2\alpha)\hat\Pi. $$

### Gaussian smoothing kernels between P, W and Q (derived normalization correction)
Appears in: P-019; P-015; C-017

With the measure $d^2\alpha/\pi$, the smoothing kernels from $P$ to $W$ and from $W$ to $Q$ are $2e^{-2|\alpha-\beta|^2}$, and from $P$ to $Q$ it is $e^{-|\alpha-\beta|^2}$.

$$ W=\int\frac{d^2\beta}{\pi}\,2e^{-2|\alpha-\beta|^2}P(\beta),\qquad Q=\int\frac{d^2\beta}{\pi}\,2e^{-2|\alpha-\beta|^2}W(\beta),\qquad Q=\int\frac{d^2\beta}{\pi}e^{-|\alpha-\beta|^2}P(\beta) $$

**Also stated in P-015 (Q symbol is the Gaussian smoothing of the P symbol):**

The covariant symbol is the contravariant symbol smoothed with the unit Gaussian.

$$ Q_{\hat A}(\alpha)=\int\frac{d^2\beta}{\pi}e^{-|\alpha-\beta|^2}P_{\hat A}(\beta) $$

**Also stated in C-017 (Wick, anti-Wick and Weyl symbols related by Gaussian smoothing):**

The Weyl symbol is the anti-Wick symbol smoothed by a Gaussian of width set by hbar, and the Wick symbol is the Weyl symbol smoothed once more by the same Gaussian.

$$ A_w=G_\hbar*A^c=e^{\hbar\Delta/4}A^c,\qquad A_c=G_\hbar*A_w=e^{\hbar\Delta/4}A_w,\qquad G_\hbar(X)=(\pi\hbar)^{-n}e^{-|X|^2/\hbar} $$

### Ordered characteristic functions differ by Gaussian factors
Appears in: P-020, P2-035; G-028, G-030 (s-ordered form $C(\lambda,s)=e^{s|\lambda|^2/2}C_W$)

Normal, Weyl and antinormal characteristic functions differ only by Gaussian factors $e^{\pm|\gamma|^2/2}$, because the ordered exponentials differ from $\hat D(\gamma)$ by those scalars.

$$ \chi_W=\operatorname{Tr}[\hat\rho\hat D(\gamma)],\ \chi_N=\operatorname{Tr}[\hat\rho e^{\gamma\hat a^\dagger}e^{-\bar\gamma\hat a}],\ \chi_A=\operatorname{Tr}[\hat\rho e^{-\bar\gamma\hat a}e^{\gamma\hat a^\dagger}];\quad \chi_N=e^{|\gamma|^2/2}\chi_W,\ \chi_A=e^{-|\gamma|^2/2}\chi_W;\quad \chi_k=e^{-\sigma_k|\gamma|^2}\chi_N,\ \sigma_N=0,\ \sigma_O=\tfrac12,\ \sigma_A=1 $$

**Also stated in G-028 (Three ordered characteristic functions and their Gaussian conversion):**

The Weyl, normal and antinormal characteristic functions are traces of the displacement in three orders; they differ only by the Gaussian factors e^{-+|lambda|^2/2}, which is the s-ordered family.

$$ C_W(\lambda)=\operatorname{Tr}(\hat\rho \hat D(\lambda)),\quad C_N(\lambda)=\operatorname{Tr}(\hat\rho e^{\lambda \hat a^\dagger}e^{-\lambda^*\hat a}),\quad C_A(\lambda)=\operatorname{Tr}(\hat\rho e^{-\lambda^*\hat a}e^{\lambda \hat a^\dagger}),\quad C_W=e^{-|\lambda|^2/2}C_N=e^{|\lambda|^2/2}C_A. $$

### Normally ordered moments from the normal characteristic function
Appears in: P-020; G-029

Normally ordered moments are mixed Wirtinger derivatives of the normal characteristic function at the origin.

$$ \langle(\hat a^\dagger)^m\hat a^n\rangle=(-1)^n\left.\partial_\gamma^m\partial_{\bar\gamma}^n\chi_N(\gamma)\right|_0 $$

**Also stated in G-029 (Ordered moments as derivatives of characteristic functions):**

Ordered moments are read off as mixed derivatives of the matching ordered characteristic function at the origin.

$$ \langle(\hat a^\dagger)^m \hat a^n\rangle=\left.\partial_\lambda^m\partial_{(-\lambda^*)}^nC_N\right|_0,\quad\langle \hat a^m(\hat a^\dagger)^n\rangle=\left.\partial_\lambda^n\partial_{(-\lambda^*)}^mC_A\right|_0. $$

### Exact attenuation solution for ordered characteristic functions
Appears in: P2-036, P2-039 (thermal $\mu(t)=p\nu_0+q\nu$)

Under thermal attenuation each ordered characteristic function is the initial one at contracted argument times a Gaussian that grows to the noise level.

$$ \chi_k(\eta,t)=e^{-q(\nu+\sigma_k)|\eta|^2}\chi_k(\sqrt p\,\eta,0),\qquad \partial_t\chi_k=-\gamma[\tfrac12\eta_j\partial_{\eta_j}+(\nu+\sigma_k)|\eta|^2]\chi_k,\quad p=e^{-\gamma t},\ q=1-p $$

### Fokker-Planck drift, diffusion and covariance evolution
Appears in: P2-037

The ordered quasiprobabilities obey Fokker-Planck equations with drift $-\gamma\alpha/2$, and their covariance relaxes linearly to the noise level plus the ordering offset.

$$ W_N=P,\ W_O=W_{\rm Wigner},\ W_A=Q/\pi;\qquad D(t)=pD(0)+q(\nu+\sigma_k)I_2/2 $$

### Gaussian convolution kernel for attenuation (sign corrected)
Appears in: P2-038

Attenuation acts on each ordered quasiprobability as convolution with a normalized Gaussian centered at $\sqrt p\,\alpha'$.

$$ W_k(\alpha,t)=\int d^2\alpha'\,W_k(\alpha',0)\frac{\exp[-|\alpha-\sqrt p\,\alpha'|^2/(q(\nu+\sigma_k))]}{\pi q(\nu+\sigma_k)} $$

### Weyl quantization kernel convention
Appears in: C-013

An operator acts through its Weyl symbol evaluated at the midpoint, and the symbol is recovered from the integral kernel by a Fourier transform in the difference variable.

$$ [\hat A\psi](x)=(2\pi\hbar)^{-n}\int A((x+y)/2,\xi)e^{i(x-y)\cdot\xi/\hbar}\psi(y)\,dy\,d\xi,\qquad A(x,\xi)=\int e^{-iu\cdot\xi/\hbar}K(x+u/2,x-u/2)\,du $$

### Symbol traces and Wigner overlaps
Appears in: C-014

With the kernel convention of C-013, Hilbert-Schmidt norms, traces of products and expectation values become phase-space integrals of symbols with the factor (2 pi hbar)^{-n}.

$$ \|\hat A\|_{\rm HS}^2=(2\pi\hbar)^{-n}\int|A|^2,\quad\operatorname{Tr}(\hat A\hat B)=(2\pi\hbar)^{-n}\int AB,\quad \langle\psi,\hat A\phi\rangle=(2\pi\hbar)^{-n}\int A\,W_{\psi,\phi},\quad |\langle\psi,\phi\rangle|^2=(2\pi\hbar)^{-n}\int W_\psi W_\phi $$

Conditions: W is the Weyl symbol of the rank-one operator (unnormalized, integral (2 pi hbar)^n).

### Moyal product and exact quadratic commutator
Appears in: C-016

The Weyl symbol of a product is the Moyal star product, and when one factor is at most quadratic the commutator is exactly i hbar times the quantized Poisson bracket.

$$ A\star B=A\exp\!\left[\frac{i\hbar}{2}\left(\overleftarrow\partial_q\cdot\overrightarrow\partial_p-\overleftarrow\partial_p\cdot\overrightarrow\partial_q\right)\right]B,\qquad [\hat H,\hat B]=i\hbar\,\widehat{\{H,B\}}\ (\deg H\le2) $$

### Weyl symbol of a metaplectic operator (Cayley form)
Appears in: C-030; oscillator case C-053

Away from det(I + F) = 0, the Weyl symbol of a metaplectic operator is a Gaussian in the Cayley transform of F, with prefactor 2^n |det(I + F)|^{-1/2} and a path-dependent phase.

$$ R_w(F;X)=\frac{2^n e^{i\pi\nu}}{|\det(I_{2n}+F)|^{1/2}}\exp\!\left[-\frac{i}{\hbar}X^TJ(I_{2n}-F)(I_{2n}+F)^{-1}X\right] $$

Conditions: det(I+F) != 0; nu and phase from the path.

### Translation coefficient of a metaplectic operator
Appears in: C-031

The translation coefficient Tr[T(-z) R(F)] of a metaplectic operator is the Fresnel transform of its Weyl symbol: a Gaussian in the inverse Cayley transform with a signature phase.

$$ \widetilde R(F;z)=\frac{c_w\,e^{-i\pi\operatorname{sgn}(Q_F)/4}}{2^n|\det Q_F|^{1/2}}\exp\!\left[-\frac{i}{4\hbar}z^TJ(I+F)(I-F)^{-1}z\right],\quad Q_F=J(I-F)(I+F)^{-1};\qquad \widetilde R=\frac{1}{2i\sin(t/2)}e^{\frac{i}{4\hbar}\cot(t/2)|z|^2}\ \text{for}\ \hat R=e^{-it(\hat N+1/2)} $$

Conditions: I + F and I - F invertible; 0 < t < pi for the example.

### Wigner function of a one-mode squeezed Gaussian
Appears in: C-049

The Wigner function of the one-mode squeezed Gaussian e^{-delta x^2/(2 hbar)} is a tilted Gaussian ellipse fixed by the real and imaginary parts of delta.

$$ W(q,p)=2\exp[-q^2\Re\delta/\hbar-(p+q\Im\delta)^2/(\hbar\Re\delta)],\qquad \psi\propto e^{-\delta x^2/(2\hbar)},\ \delta=\frac{1-\omega}{1+\omega} $$

### Weyl symbol of the oscillator exponential
Appears in: C-053 (special case of C-030)

The Weyl symbol of the oscillator propagator is (cos(t/2))^{-n} e^{i tan(t/2)|X|^2}, and in imaginary time (cosh(tau/2))^{-n} e^{-tanh(tau/2)|X|^2}.

$$ [e^{it\hat H_{\rm os}}]_w=(\cos(t/2))^{-n}\exp[i\tan(t/2)|X|^2],\qquad [e^{-\tau\hat H_{\rm os}}]_w=(\cosh(\tau/2))^{-n}\exp[-\tanh(\tau/2)|X|^2] $$

Conditions: hbar = 1, |t| < pi.

### Fermionic Weyl translation expansion
Appears in: C-074

An operator expands in fermionic translations with coefficient Str(A T(xi)).

$$ \hat A=\int A_w(\xi)\hat T(-\xi)d^2\xi,\qquad A_w(\xi)=\operatorname{Str}(\hat A\hat T(\xi)) $$

### Normal, antinormal and s-ordered fermionic translations
Appears in: C-075

Normal and antinormal fermionic translations differ from the symmetric one by the scalar Grassmann factor e^{+gamma^* gamma/2} or e^{-gamma^* gamma/2}.

$$ \hat T_N(\gamma)=e^{\hat c^\dagger\gamma}e^{-\gamma^*\hat c}=\hat T(\gamma)e^{\gamma^*\gamma/2},\quad \hat T_A(\gamma)=e^{-\gamma^*\hat c}e^{\hat c^\dagger\gamma}=\hat T(\gamma)e^{-\gamma^*\gamma/2},\quad \hat T_s=\hat T e^{s\gamma^*\gamma/2} $$

### Ordinary trace from mixed fermionic symbols
Appears in: C-076

The ordinary trace of a product of fermionic operators is a paired Berezin integral of one contravariant and one covariant symbol at a reflected, doubled argument, and the trace of one operator is 2^n times its symbol at zero.

$$ \operatorname{Tr}(\hat H\hat G)=2^{-n}\widehat I_\zeta\!\left[H^w(\zeta)G_w(-2\zeta)\right],\qquad \operatorname{Tr}\hat A=2^nA^w(0) $$

### Fermionic star product
Appears in: C-077

The fermionic star product contracts odd derivatives symmetrically, and the propagator's symbol solves the star Schrodinger equation.

$$ G\star H=G\exp\!\left[\tfrac12\left(\overleftarrow\partial_{\alpha^*}\cdot\overrightarrow\partial_\alpha+\overleftarrow\partial_\alpha\cdot\overrightarrow\partial_{\alpha^*}\right)\right]H,\qquad i\partial_tU^w=H\star U^w $$

### Fermi oscillator exponential and its Weyl symbol
Appears in: C-078

The exponential of one fermion number operator is linear in it, and its Weyl symbol is (e^z+1)/2 + (e^z-1) alpha^* alpha.

$$ e^{z\hat c^\dagger\hat c}=\hat I+(e^z-1)\hat c^\dagger\hat c,\quad (\hat c^\dagger\hat c)^w=\alpha^*\alpha+\tfrac12,\quad [e^{z\hat c^\dagger\hat c}]^w=\frac{e^z+1}2\exp[2\tanh(z/2)\alpha^*\alpha]=\frac{e^z+1}2+(e^z-1)\alpha^*\alpha $$

### Fermionic Cayley Weyl symbol (corrected sign)
Appears in: C-085

The Weyl symbol of a quadratic fermion propagator is a Grassmann Gaussian in the Cayley transform of F_t; the atlas corrects the sign of the exponent.

$$ U_t^w(\xi)=\det\!\left(\frac{F_t+I}{2}\right)^{1/2}\exp\!\left[-\binom{\xi^*}{\xi}\cdot(F_t-I)(F_t+I)^{-1}\binom\xi{\xi^*}\right]\ \text{(atlas correction: minus sign)} $$

### Super star product
Appears in: C-089

The super star product adds the antisymmetric bosonic contraction to the symmetric fermionic one.

$$ G\star H=G\exp\tfrac12\left[\overleftarrow\partial_\zeta\overrightarrow\partial_{\bar\zeta}-\overleftarrow\partial_{\bar\zeta}\overrightarrow\partial_\zeta+\overleftarrow\partial_{\alpha^*}\overrightarrow\partial_\alpha+\overleftarrow\partial_\alpha\overrightarrow\partial_{\alpha^*}\right]H $$

### Coherent-projector differential correspondences
Appears in: G-025

Multiplying a coherent projector by a^dagger on the left or a on the right acts on its matrix elements as a first-order differential operator in alpha and alpha^*.

$$ \hat a^\dagger|\alpha\rangle\langle\alpha|=(\alpha^*+\partial_\alpha)|\alpha\rangle\langle\alpha|,\quad |\alpha\rangle\langle\alpha|\hat a=(\alpha+\partial_{\alpha^*})|\alpha\rangle\langle\alpha|. $$

Conditions: $\alpha,\alpha^*$ independent Wirtinger variables.

### P representation, optical equivalence, and Q function
Appears in: G-026

A density operator is a P-weighted mixture of coherent projectors, normal-ordered expectations are P-averages of the symbol, Q is the coherent diagonal over pi, and normal ordering of n^2 drops a term n.

$$ \hat\rho=\int P(\alpha)|\alpha\rangle\langle\alpha|d^2\alpha,\quad\langle{:G(\hat a^\dagger,\hat a):}\rangle=\int P(\alpha)G(\alpha^*,\alpha)d^2\alpha,\quad Q(\alpha)=\pi^{-1}\langle\alpha|\hat\rho|\alpha\rangle,\quad {:\hat n^2:}=\hat a^{\dagger2}\hat a^2\ne\hat n^2. $$

### Glauber P function of a number state
Appears in: G-027

The P function of a coherent state is a delta at its amplitude; that of a number state is a derivative of the delta, defined only under integration against test functions.

$$ P_{|\beta\rangle}(\alpha)=\delta^{(2)}(\alpha-\beta),\qquad P_m(\alpha)=\frac{e^{|\alpha|^2}}{m!}\partial_\alpha^m\partial_{\alpha^*}^m\delta^{(2)}(\alpha). $$

### Complex Fourier transform between P and the normal characteristic function
Appears in: G-031

P and the normal characteristic function form a complex Fourier pair with kernel e^{lambda alpha^* - lambda^* alpha} and normalization pi^{-2}.

$$ C_N(\lambda)=\int P(\alpha)e^{\lambda\alpha^*-\lambda^*\alpha}d^2\alpha,\quad P(\alpha)=\frac1{\pi^2}\int C_N(\lambda)e^{\lambda^*\alpha-\lambda\alpha^*}d^2\lambda. $$

### Q and Wigner functions of coherent and number states
Appears in: G-033

Coherent states have Gaussian Q and W centered at their amplitude; number states have a Poisson-shaped Q and a Laguerre-ringed W of sign (-1)^m at the origin.

$$ Q_{|\beta\rangle}(\alpha)=\pi^{-1}e^{-|\alpha-\beta|^2},\quad W_{|\beta\rangle}(\alpha)=\frac2\pi e^{-2|\alpha-\beta|^2},\quad Q_{|m\rangle}(\alpha)=\frac{e^{-|\alpha|^2}|\alpha|^{2m}}{\pi m!},\quad W_{|m\rangle}(\alpha)=\frac2\pi(-1)^m L_m(4|\alpha|^2)e^{-2|\alpha|^2}. $$

### Coherent overlap of squeezed vacuum and Q function of a displaced squeezed state
Appears in: G-063

The coherent overlap of squeezed vacuum is a Gaussian in beta with the e^{i theta} tanh r quadratic term, and the Q function of the displaced state is that overlap squared, centered at alpha.

$$ \langle\beta|\hat S(\xi)|0\rangle=\frac1{\sqrt{\cosh r}}e^{-|\beta|^2/2-e^{i\theta}\beta^{*2}\tanh r/2},\quad Q(\beta)=\frac1{\pi\cosh r}\exp\{-|\delta|^2-\tanh r\,\operatorname{Re}(e^{-i\theta}\delta^2)\},\ \delta=\beta-\alpha. $$

### Wigner function of a displaced squeezed state
Appears in: G-064

A displaced squeezed state has a Gaussian Wigner function centered at alpha, narrowed by e^{-r} along the squeezing axis at angle theta/2 and widened by e^{r} across it.

$$ W(\beta)=\frac2\pi\exp[-2e^{2r}x^2-2e^{-2r}y^2],\qquad x+iy=e^{-i\theta/2}(\beta-\alpha). $$

### Cat Wigner function
Appears in: G-084

The cat Wigner function is two Gaussians at +-alpha plus an interference fringe cos(4 alpha y + Phi) at the origin, over the cat normalization.

$$ W_{C_\Phi}(x,y)=\frac{e^{-2(x-\alpha)^2-2y^2}+e^{-2(x+\alpha)^2-2y^2}+2e^{-2x^2-2y^2}\cos(4\alpha y+\Phi)}{\pi(1+e^{-2\alpha^2}\cos\Phi)}. $$

Conditions: real $\alpha$, $\beta=x+iy$.

## Grassmann calculus and fermionic coherent states

### Fermionic symbols: which side the derivative acts on
Appears in: B-016

For fermions, multiplying by a on the right or by a* on the left keeps normal order, so the normal symbol is just multiplied on the same side by the Grassmann variable.

$$ N_{\hat A\hat a}=N_{\hat A}a,\qquad N_{\hat a^*\hat A}=a^*N_{\hat A} $$

### Grassmann algebra relations
Appears in: B-029

Grassmann generators anticommute and square to zero, so n of them span an algebra of dimension 2^n.

$$ x_i x_j=-x_jx_i,\qquad x_i^2=0,\qquad \dim=2^n $$

### Left versus right derivatives
Appears in: B-030

On a monomial of degree p, the right and left Grassmann derivatives differ by the sign (-1)^{p-1}.

$$ \partial^{R}_{x}m=(-1)^{p-1}\,\partial^{L}_{x}m\quad\text{on a monomial }m\text{ of degree }p $$

### Berezin integration rules
Appears in: B-031

Berezin integration kills constants and gives 1 on the integration variable itself.

$$ \int dx_i=0,\qquad \int x_i\,dx_i=1 $$

### Grassmann change of variables
Appears in: B-032

Under a linear change of Grassmann variables the integral picks up the inverse determinant, the opposite of commuting variables.

$$ \int f(x)\,dx_n\cdots dx_1=(\det A)^{-1}\int f(Ay)\,dy_n\cdots dy_1 $$

### Clifford operators from multiplication and derivative
Appears in: B-034

Grassmann multiplication and left differentiation satisfy the CAR, so their sum and difference generate a Clifford algebra.

$$ \{\hat Q_i,\hat Q_j\}=\{\hat P_i,\hat P_j\}=2\delta_{ij}\hat I,\quad \{\hat P_i,\hat Q_j\}=0,\qquad \hat Q_j=\hat X_j+\hat D_j,\ \hat P_j=(\hat X_j-\hat D_j)/i $$

### Grassmann derivatives: Leibniz rule and left/right relation
Appears in: C-066

The Grassmann derivative obeys a graded Leibniz rule, and the right and left derivatives differ by a parity sign.

$$ \partial_\theta(fg)=(\partial_\theta f)g+(-1)^{\pi(f)}f(\partial_\theta g),\qquad f\overleftarrow\partial_\theta=(-1)^{1+\pi(f)}\overrightarrow\partial_\theta f $$

### Berezin integration: determinant rule and paired normalization
Appears in: C-067; the same paired functional is used in C-068, C-069, C-098, C-099

Berezin integration rescales by det A under a linear change of generators, and the paired convention normalizes e^{-theta^* theta} to 1.

$$ \int f(A\theta)\,d\theta=\det A\int f(\theta)\,d\theta,\qquad \widehat I_\theta\Big[\prod_j\theta_j\theta_j^*\Big]=1,\ \widehat I_\theta[e^{-\theta^*\theta}]=1,\ \widehat I_\theta=\partial_{\theta_1^*}\partial_{\theta_1}\cdots\partial_{\theta_n^*}\partial_{\theta_n} $$

### Normal-ordered monomial expansion of fermionic operators
Appears in: C-070

Every operator on n fermionic modes is a unique combination of normal-ordered monomials.

$$ \hat A=\sum_{\epsilon,\epsilon'\in\{0,1\}^n}A_{\epsilon\epsilon'}(\hat c^\dagger)^\epsilon\hat c^{\epsilon'} $$

### Grassmann Weyl translation
Appears in: C-071

The fermionic translation with Grassmann labels has an even generator, so the central BCH formula factorizes it and composes it up to a scalar Grassmann phase, and it shifts c by gamma.

$$ \hat T(\gamma)=e^{\hat c^\dagger\cdot\gamma-\gamma^*\cdot\hat c}=e^{-\gamma^*\cdot\gamma/2}e^{\hat c^\dagger\cdot\gamma}e^{-\gamma^*\cdot\hat c},\quad \hat T(\alpha)\hat T(\gamma)=\hat T(\alpha+\gamma)e^{-(\alpha^*\cdot\gamma+\alpha\cdot\gamma^*)/2},\quad \hat T^{-1}\hat c\hat T=\hat c+\gamma,\ \hat T^{-1}\hat c^\dagger\hat T=\hat c^\dagger+\gamma^* $$

### Fermionic coherent state, overlap, resolution
Appears in: C-072

The fermionic coherent state is an exponential of theta.gamma, an eigenvector of each c_j with Grassmann eigenvalue gamma_j; overlaps are Gaussian in the labels and the states resolve the identity.

$$ \psi_\gamma=e^{-\gamma^*\gamma/2}e^{\theta\cdot\gamma},\ \hat c_j\psi_\gamma=\gamma_j\psi_\gamma,\quad \langle\psi_\gamma,\psi_\alpha\rangle=\exp[\gamma^*\alpha-(\gamma^*\gamma+\alpha^*\alpha)/2],\quad \psi=\int\langle\psi_\gamma,\psi\rangle\psi_\gamma\,d\gamma^*d\gamma $$

### Supertrace and Berezinian
Appears in: C-073

The supertrace weights the trace by the parity operator and is graded-cyclic; the Berezinian is the super analogue of the determinant and satisfies Ber e^H = e^{Str H}.

$$ \operatorname{Str}\hat A=\operatorname{Tr}(\hat P\hat A),\ \hat P=(-1)^{\hat N},\ \operatorname{Str}(\hat A\hat B)=(-1)^{\pi(A)\pi(B)}\operatorname{Str}(\hat B\hat A),\quad \operatorname{Ber}\begin{pmatrix}A&B\\C&D\end{pmatrix}=\frac{\det(A-BD^{-1}C)}{\det D},\ \operatorname{Ber}(e^H)=e^{\operatorname{Str}H} $$

### Mixed boson-fermion displacement
Appears in: C-088

A mixed boson-fermion displacement factors into its bosonic and fermionic normal-ordered parts, since the even generators of the two parts commute.

$$ \hat T(z,\gamma)=e^{-(|z|^2+\gamma^*\gamma)/2}e^{z\hat a_B^\dagger}e^{\hat c_F^\dagger\gamma}e^{-\bar z\hat a_B}e^{-\gamma^*\hat c_F} $$

### Supercommutator
Appears in: C-090

The supercommutator is the commutator with a sign that turns it into an anticommutator when both operators are odd.

$$ [\hat A,\hat B]_s=\hat A\hat B-(-1)^{\pi(\hat A)\pi(\hat B)}\hat B\hat A $$

### Grassmann Fourier involution
Appears in: C-098

The Grassmann Fourier transform is an involution in the paired convention; it maps 1 to the ordered top product of the alpha_j alpha_j^*.

$$ (\widehat{\mathcal F}f)(\alpha,\alpha^*)=\widehat I_\xi[e^{\xi^*\cdot\alpha-\alpha^*\cdot\xi}f],\quad \widehat{\mathcal F}^2=\hat I,\quad \widehat{\mathcal F}1=\alpha\alpha^*,\ \widehat{\mathcal F}\xi=\alpha,\ \widehat{\mathcal F}\xi^*=-\alpha^*,\ \widehat{\mathcal F}(\xi\xi^*)=1,\quad \widehat{\mathcal F}1=\prod_j\alpha_j\alpha_j^* $$

### Holomorphic Grassmann reproducing kernel
Appears in: C-099

The exponential Grassmann kernel e^{-(theta - gamma).gamma^*} reproduces every holomorphic Grassmann function under paired integration.

$$ \psi(\theta)=\widehat I_\gamma[\psi(\gamma)E(\theta,\gamma,\gamma^*)],\qquad E=e^{-(\theta-\gamma)\cdot\gamma^*}=\prod_j[1-(\theta_j-\gamma_j)\gamma_j^*] $$

## Gaussian integrals and Pfaffians

### Real Gaussian integral
Appears in: B-022

The normalized real Gaussian integral of e^{-(Af, f)} is det A^{-1/2}.

$$ \int e^{-(\hat Af,f)}\Pi df=(\det \hat A)^{-1/2},\qquad \hat A>0 $$

Conditions: normalized measure pi^{-n/2} d^n f; infinite dimension needs A - I trace class.

### Doubled complex Gaussian with sources
Appears in: B-023

A complex Gaussian with both holomorphic and antiholomorphic quadratic terms and linear sources integrates to the inverse square root of a row-swapped block determinant times the exponential of half the inverse form on the sources.

$$ \int e^{-\frac12(f\ f^*)\widehat{\mathcal A}\binom f{f^*}+(\varphi_1\ \varphi_2)\binom f{f^*}}\Pi df\,df^* =\det\!\begin{pmatrix}\hat A_{21}&\hat A_{22}\\\hat A_{11}&\hat A_{12}\end{pmatrix}^{-1/2} e^{\frac12(\varphi_1\ \varphi_2)\widehat{\mathcal A}^{-1}\binom{\varphi_1}{\varphi_2}} $$

Conditions: A transpose-symmetric, A12 - I and A21 - I trace class, A11, A22 Hilbert-Schmidt, real part of the form bounded below.

### Fermionic Gaussian is a Pfaffian
Appears in: B-033

The Berezin integral of a Grassmann Gaussian with skew matrix M is the Pfaffian of M.

$$ \int e^{\frac12x^TMx}\,dx_n\cdots dx_1=\operatorname{Pf}(M),\quad M^T=-M $$

### Doubled fermionic Gaussian (corrected determinant)
Appears in: B-036

The paired Grassmann Gaussian integral equals the square root of the unpermuted block determinant, fixed in sign as (-1)^{n(n-1)/2} times the Pfaffian.

$$ \int e^{\frac12(\alpha\ \alpha^*)\widehat{\mathcal A}\binom\alpha{\alpha^*}}\Pi d\alpha^*d\alpha =\det\!\begin{pmatrix}\hat A_{11}&\hat A_{12}\\\hat A_{21}&\hat A_{22}\end{pmatrix}^{1/2} = (-1)^{n(n-1)/2}\operatorname{Pf}(\mathcal A) $$

Conditions: normalized paired integration, int e^{-alpha* alpha} = 1; block order (alpha_1..alpha_n, alpha*_1..alpha*_n).

### Determinant and log-determinant series
Appears in: P-048

The determinant is the exponential of the trace of the logarithm, and $-\log\det(I-X)$ is the trace power series, convergent for spectral radius below 1.

$$ \det A=e^{\operatorname{Tr}\log A},\qquad -\log\det(I_N-X)=\sum_{m\ge1}\frac{\operatorname{Tr}X^m}{m},\ \rho(X)<1;\qquad \psi_Z(W)=B(Z,Z)^{-1/2}B(Z,W) $$

### Fermionic Gaussian integral with sources
Appears in: C-068

The paired Berezin integral of e^{-theta^* B theta} is det B, and linear sources complete the square to det B e^{gamma^* B^{-1} gamma}.

$$ \widehat I_\theta e^{-\theta^*B\theta}=\det B,\qquad \widehat I_\theta e^{-\theta^*B\theta+\gamma^*\theta+\theta^*\gamma}=\det B\,e^{\gamma^*B^{-1}\gamma} $$

### Pfaffian Gaussian integral
Appears in: C-069

A general Grassmann Gaussian integrates to a Pfaffian with the block-to-paired sign (-1)^{n(n-1)/2}, and sources put the inverse matrix in the exponent.

$$ \widehat I_\gamma e^{\Phi+\gamma^*\xi-\xi^*\gamma}=(-1)^{n(n-1)/2}\operatorname{Pf}(\Lambda)\exp\!\left[\tfrac12\binom{\xi^*}{\xi}^{T}\Lambda^{-1}\binom{\xi^*}{\xi}\right],\quad \Phi=\tfrac12x^T\Lambda x $$

### Fermionic translation coefficient and Pfaffian amplitude
Appears in: C-086

The translation coefficient of a quadratic fermion propagator is a Grassmann Gaussian in the inverse Cayley transform, with amplitude the supertrace, whose square is (-1)^n det(F_t - I) and which equals an ordered Pfaffian.

$$ U_{t,w}(\gamma)=a_t\exp\!\left[-\tfrac14\binom\gamma{\gamma^*}^{T}J_F(F_t+I)(F_t-I)^{-1}\binom\gamma{\gamma^*}\right],\ a_t=\operatorname{Str}\hat U_t,\ a_t^2=(-1)^n\det(F_t-I),\ a_t=(-1)^{n(n-1)/2}p_t\operatorname{Pf}(-2J_FC_t) $$

### Commuting Gaussian integral and Fresnel factor
Appears in: C-091

The real Gaussian integral of e^{-x^T M x/2 + j^T x} is (2 pi)^{m/2} det(M)^{-1/2} e^{j^T M^{-1} j/2}, and its oscillatory limit carries the Fresnel phase e^{i pi sgn(A)/4}.

$$ \int_{\mathbb R^m}e^{-x^TMx/2+j^Tx}dx=(2\pi)^{m/2}\det(M)^{-1/2}e^{j^TM^{-1}j/2},\qquad \text{Fresnel factor } e^{i\pi\operatorname{sgn}(A)/4}|\det A|^{-1/2} $$

Conditions: M symmetric, Re M > 0; root positive on positive M.

## Thermal states and traces

### Traces of second-quantized one-particle operators
Appears in: X-007; B-107

The Fock-space trace of the second quantization of a one-particle matrix $M$ is the inverse determinant of $I-M$ for bosons and the determinant of $I+M$ for fermions.

$$ \operatorname{Tr}_{\mathcal F_B}\widehat\Gamma(M)=\det(I_N-M)^{-1},\qquad \operatorname{Tr}_{\mathcal F_F}\widehat\Gamma(M)=\det(I_N+M) $$

Conditions: $\widehat\Gamma(e^h)=e^{\hat a^\dagger h\hat a}$; bosons need convergence: spectral radius of $M$ below 1 (the atlas gives $0\le M<I_N$ for the positive case); infinite modes need trace class and $\|\hat M\|<1$ for bosons.

**Also stated in B-107 (Free Fock-space traces of a second-quantized one-particle map):**

The trace over bosonic Fock space of the second quantization of M is 1/det(I - M), and over fermionic Fock space det(I + M).

$$ \operatorname{Tr}_{\mathcal F_{B}}\widehat\Gamma(\hat M)=\det(\hat I-\hat M)^{-1},\qquad \operatorname{Tr}_{\mathcal F_{F}}\widehat\Gamma(\hat M)=\det(\hat I+\hat M) $$

Conditions: bosons need ||M|| < 1 (and trace class in infinite dimension); fermions need trace class only.

### Bosonic trace
Appears in: B-028

The trace of a trace-class operator is the integral of its normal symbol over the normalized complex plane (equivalently, of its kernel times e^{-a*a}).

$$ \operatorname{Tr}\hat A=\int K_{\hat A}(a^*,a)e^{-a^*a}\Pi da^*da=\int N_{\hat A}(a^*,a)\Pi da^*da $$

Conditions: A trace class.

### Fermionic trace
Appears in: B-040

The fermionic trace is a Berezin integral of the kernel with weight e^{-aa*}, or of the normal symbol with the doubled weight e^{-2aa*}.

$$ \operatorname{Tr}\hat A=\int K_{\hat A}(a^*,a)e^{-aa^*}\Pi da\,da^* =\int N_{\hat A}(a^*,a)e^{-2aa^*}\Pi da\,da^* $$

### Normalized fermionic trace
Appears in: B-041

Rescaling the symbol's arguments by 1/sqrt 2 and using the single weight gives the trace divided by 2^n, the normalized trace that is 1 on the identity.

$$ \tau(\hat A)=\int N_{\hat A}(a^*/\sqrt2,a/\sqrt2)e^{-aa^*}\Pi da\,da^*=2^{-n}\operatorname{Tr}\hat A $$

### Grand partition functions
Appears in: B-106

The grand partition function is an auxiliary-field average of 1/det(I - M) for bosons and det(I + M) for fermions, with M the one-particle Boltzmann factor dressed by the interaction.

$$ \Xi_{B}=\int\det(\hat I-\hat M)^{-1}e^{-zz^*}\Pi dz\,dz^*,\qquad \Xi_{F}=\int\det(\hat I+\hat M)e^{-zz^*}\Pi dz\,dz^*,\quad \hat M=e^{-\beta(\hat\omega+\mu\hat I)}(\hat I-\hat T) $$

### Coherent-symbol bounds on exponential traces (Berezin-Lieb, Golden-Thompson)
Appears in: P-016

The trace of a Gibbs exponential is bounded below by the integral of $e^{-tQ}$ and above by the integral of $e^{-tP}$; Golden-Thompson bounds the trace of an exponential of a sum by the trace of the product.

$$ \int d\mu(\alpha)e^{-tQ_{\hat H}(\alpha)}\le\operatorname{Tr}e^{-t\hat H}\le\int d\mu(\alpha)e^{-tP_{\hat H}(\alpha)},\qquad \operatorname{Tr}e^{-t(\hat A+\hat B)}\le\operatorname{Tr}(e^{-t\hat A}e^{-t\hat B}) $$

Conditions: self-adjoint, bounded below, trace class, convergent symbol integrals; $d\mu=d^2\alpha/\pi$.

### Weyl characteristic function of a number-diagonal state
Appears in: P-020

For a state diagonal in the number basis, the Weyl characteristic function is a Gaussian times the population-weighted sum of Laguerre polynomials.

$$ \chi_W(\gamma)=e^{-|\gamma|^2/2}\sum_np_nL_n(|\gamma|^2) $$

### Quasiclassical partition-function bounds
Appears in: P2-032

For lattice spin systems the quantum partition function is squeezed between classical ones at $\beta$ and at a slightly larger $\beta$.

$$ Z_{\rm cl}(\beta)\le Z_Q^{(N)}(\beta)\le Z_{\rm cl}(\beta(1+a/N)) $$

### Thermal oscillator master equation
Appears in: P2-034

The thermal damped-oscillator master equation makes $\langle\hat a\rangle$ decay at rate $\gamma/2$ and $\langle\hat n\rangle$ relax to $\nu$ at rate $\gamma$.

$$ \dot{\hat\rho}=\gamma(\nu+1)\left(\hat a\hat\rho \hat a^\dagger-\tfrac12\{\hat a^\dagger \hat a,\hat\rho\}\right)+\gamma\nu\left(\hat a^\dagger\hat\rho \hat a-\tfrac12\{\hat a\hat a^\dagger,\hat\rho\}\right),\quad \nu=(e^{\hbar\omega/kT}-1)^{-1} $$

### Displaced thermal density kernel
Appears in: P2-043

The coherent kernel of a displaced thermal state is a single Gaussian in the two labels, with mean occupation $\mu$ and displacement $\beta$.

$$ R(z_1,\bar z_2;t)=\frac1{1+\mu}\exp\left[\frac{\mu z_1\bar z_2+\beta z_1+\bar\beta\bar z_2-|\beta|^2}{1+\mu}\right] $$

### Populations of a displaced thermal state
Appears in: P2-044

The populations of a displaced thermal state are a geometric factor times a Gaussian times a Laguerre polynomial of negative argument.

$$ W_n=\frac{\mu^n}{(1+\mu)^{n+1}}e^{-|\beta|^2/(1+\mu)}L_n\left(-\frac{|\beta|^2}{\mu(1+\mu)}\right) $$

### Free oscillator propagator and thermal trace in coherent states
Appears in: P2-046

The coherent matrix element of $e^{-it\omega\hat a^\dagger\hat a}$ is the overlap with a rotated label, and the thermal trace is the geometric series.

$$ \langle\alpha|e^{-it\omega\hat a^\dagger\hat a}|\beta\rangle=\exp\left[\bar\alpha\beta e^{-i\omega t}-\tfrac12(|\alpha|^2+|\beta|^2)\right],\qquad Z=(1-e^{-\omega/T})^{-1} $$

### Stable oscillator proper-time trace
Appears in: P2-048

The trace of $e^{-s(2\hat N+1)}$ is $1/(2\sinh s)$ for $s>0$; the source's diagonalizing $\hat U$ is a complex continuation, not a unitary.

$$ \operatorname{Tr}e^{-s(2\hat N+\hat I)}=[2\sinh s]^{-1},\ s>0;\qquad \hat U=e^{i\pi(\hat b^{\dagger2}-\hat b^2)/8}\ \text{(not unitary)} $$

### Berezin-Lieb bounds for spin
Appears in: C-043

The trace of a convex function of a spin Hamiltonian is bounded below by the integral over its covariant symbol and above by the integral over its contravariant symbol; for J_l these symbols are j n_l and (j+1) n_l.

$$ \int\chi(H_c)d\mu\le\operatorname{Tr}\chi(\hat H)\le\int\chi(H^c)d\mu,\qquad (\hat J_\ell)_c=jn_\ell,\ (\hat J_\ell)^c=(j+1)n_\ell $$

### Gaussian density operators, purity, thermal factors
Appears in: C-055

A centered Gaussian Weyl symbol is a density operator when V + iJ/2 >= 0; its purity is 2^{-n} det(V)^{-1/2}, the oscillator partition function is a product of 1/(2 sinh) factors, and each thermal mode has covariance (hbar/2) coth(beta hbar omega/2).

$$ V+iJ/2\ge0,\quad \operatorname{Tr}\hat\rho^2=2^{-n}(\det V)^{-1/2},\quad \operatorname{Tr}e^{-\beta\hat H}=\prod_j[2\sinh(\beta\hbar\omega_j/2)]^{-1},\quad V_j=\frac\hbar2\coth(\beta\hbar\omega_j/2)I_2 $$

### Fermionic partition determinant
Appears in: C-079

Traces of fermionic quadratic exponentials are determinants: det(I + e^h) for the trace and det(I - e^h) for the supertrace.

$$ \operatorname{Tr}e^{\sum_jz_j\hat c_j^\dagger\hat c_j}=\prod_j(1+e^{z_j}),\quad \operatorname{Tr}e^{\hat c^\dagger h\hat c}=\det(I_n+e^h),\quad \operatorname{Str}e^{\hat c^\dagger h\hat c}=\det(I_n-e^h) $$

### Thermal state as a number exponential
Appears in: G-009

A thermal state is a normalized exponential of the number operator; its eigenvalues form a geometric distribution in the mean photon number, and the partition function carries the zero-point factor.

$$ Z=\frac{e^{-x/2}}{1-e^{-x}},\quad \hat\rho_{\rm th}=(1-e^{-x})e^{-x\hat n}=\frac1{1+\bar n}\left(\frac{\bar n}{1+\bar n}\right)^{\hat n},\quad\bar n=(e^x-1)^{-1}. $$

Conditions: $x=\hbar\omega/(k_BT)>0$.

### Thermal factorial moments and number variance
Appears in: G-010

In a thermal state the r-th factorial moment is r! nbar^r, which gives the Bose-Einstein number variance nbar(1 + nbar).

$$ \langle(\hat a^\dagger)^r \hat a^r\rangle_{\rm th}=r!\bar n^r,\quad\langle \hat n^2\rangle=\bar n+2\bar n^2,\quad\operatorname{Var}\hat n=\bar n(1+\bar n). $$

### Thermal characteristic functions and phase-space densities
Appears in: G-032

A thermal state has Gaussian characteristic functions with widths nbar + 1, nbar, nbar + 1/2, a Gaussian P of width nbar, and a Gaussian Q of width nbar + 1.

$$ C_A=e^{-(\bar n+1)|\lambda|^2},\quad C_N=e^{-\bar n|\lambda|^2},\quad C_W=e^{-(\bar n+1/2)|\lambda|^2},\quad P(\alpha)=\frac{e^{-|\alpha|^2/\bar n}}{\pi\bar n},\quad Q(\alpha)=\frac{e^{-|\alpha|^2/(1+\bar n)}}{\pi(1+\bar n)}. $$

### Reduced twin-beam state is thermal; entanglement entropy
Appears in: G-074, G-104 (Schmidt spectrum justification)

Tracing one twin-beam mode leaves a thermal state with tanh^2 r as the Boltzmann factor, whose entropy is the Bose-Einstein entropy of mean sinh^2 r.

$$ \hat\rho_a=\operatorname{sech}^2r\sum_{m\ge0}\tanh^{2m}r|m\rangle\langle m|,\quad e^{-\hbar\omega_a/(k_BT_{\rm eff})}=\tanh^2r,\quad S(\hat\rho_a)=(N+1)\ln(N+1)-N\ln N. $$

### Thermal state in normally ordered exponential form
Appears in: G-097

Combining the previous entry with the thermal exponential writes the thermal state as a normal-ordered Gaussian with weight 1/(1 + nbar).

$$ \hat\rho_{\rm th}=\frac1{1+\bar n}{:\exp[-\hat n/(1+\bar n)]:}. $$

### Amplitude damping: Kraus form and matrix-element solution
Appears in: G-098

Vacuum amplitude damping is a Kraus map with operators sqrt((1-T)^l/l!) T^{n/2} a^l, whose Fock components give the stated matrix-element sum and which solves the damping master equation.

$$ \rho_{mn}(t)=T^{(m+n)/2}\sum_{\ell\ge0}\frac{(1-T)^\ell}{\ell!}\sqrt{\frac{(m+\ell)!(n+\ell)!}{m!n!}}\,\rho_{m+\ell,n+\ell}(0),\qquad\hat M_\ell=\sqrt{(1-T)^\ell/\ell!}\,T^{\hat n/2}\hat a^\ell,\ \sum_\ell\hat M_\ell^\dagger\hat M_\ell=\hat I. $$

Conditions: $T=e^{-\gamma t}$, $\dot{\hat\rho}=\gamma(\hat a\hat\rho\hat a^\dagger-\{\hat n,\hat\rho\}/2)$.

### Damping of coherent dyads and of the normal characteristic function
Appears in: G-099

Damping maps a coherent dyad to the dyad of shrunken amplitudes times a coherence factor that decays with the squared distance between alpha and beta.

$$ C_N(\lambda,t)=C_N(\sqrt T\lambda,0),\qquad\widehat\Lambda_T(|\alpha\rangle\langle\beta|)=\exp\!\left[-\frac{1-T}{2}(|\alpha|^2+|\beta|^2-2\beta^*\alpha)\right]|\sqrt T\alpha\rangle\langle\sqrt T\beta|. $$

### Coherent trace formula and Schmidt spectrum
Appears in: G-104

The trace of an operator is the integral of its coherent diagonal over d^2 alpha/pi, and the entropy of a pure bipartite state is the Shannon entropy of its squared Schmidt coefficients.

$$ \operatorname{Tr}\hat F=\int\langle\alpha|\hat F|\alpha\rangle\frac{d^2\alpha}{\pi},\quad\hat\rho(t)=\hat U\hat\rho(0)\hat U^\dagger,\quad S=-\sum_j|g_j|^2\ln|g_j|^2. $$

Conditions: trace-class $\hat F$; orthonormal Schmidt families.

## Wick, Dyson and auxiliary fields

### Bound for exponentials of bounded perturbations
Appears in: B-068

If A is a self-adjoint operator plus a bounded perturbation T, then e^{itA} grows at most like e^{|t| ||T||}.

$$ \|e^{it\hat A}\|\le e^{|t|\|\hat T\|},\qquad \hat A=\hat S+\hat T,\ \hat S=\hat S^* $$

### Pair potential on the n-particle sector
Appears in: B-093

A quartic normal-ordered interaction acts on the n-particle sector as a sum over particle pairs of the pair potential; the factor 2 corresponds to the quartic term written without the usual 1/2.

$$ \hat H_n=-\sum_i\Delta_i+2\sum_{i<j}u(x_i-x_j)\hat I $$

### Dyson series with the e^{+itH} convention
Appears in: B-094

With Berezin's sign convention e^{+itH}, the interaction-picture propagator is the time-ordered exponential with powers of +i.

$$ \hat U(t,\tau)=\sum_{n\ge0}\frac{i^n}{n!}\int_\tau^t\cdots\int_\tau^t T[\hat V(t_1)\cdots \hat V(t_n)]\,dt_1\cdots dt_n $$

Conditions: U(t, tau) = e^{-itH0} e^{i(t-tau)H} e^{i tau H0}, V(t) = e^{-itH0} V e^{itH0}; formal.

### Free time dependence and contraction kernel
Appears in: B-095

In Berezin's convention the free annihilator acquires the phase e^{+it omega}, the creator e^{-it omega}, and their unequal-time commutator is the scalar D.

$$ \hat a(p,t)=e^{it\omega(p)}\hat a(p),\quad \hat a^*(p,t)=e^{-it\omega(p)}\hat a^*(p),\qquad D(p_1,t_1\mid p_2,t_2)=\delta(p_1-p_2)e^{i(t_1-t_2)\omega(p_1)} $$

### Retarded contraction and the two-field time-ordered product
Appears in: B-097

The retarded contraction is D when the annihilator's time is later and 0 otherwise, so a two-field time-ordered product is its normal product plus that contraction.

$$ \Delta(p_1,t_1\mid p_2,t_2)= \begin{cases}D(p_1,t_1\mid p_2,t_2),&t_1>t_2,\\0,&t_1\le t_2,\end{cases}\qquad T[\hat a^*(1)\hat a(2)]=\hat a^*(1)\hat a(2)+\Delta(2\mid1)\hat I $$

### Wick's theorem as an exponential of the contraction operator
Appears in: B-099

The normal symbol of a time-ordered product is obtained by applying the exponential of the second-order contraction operator to the plain product of symbols.

$$ \widehat{\mathcal D}=\int\frac{\delta}{\delta a(q,\tau)}\Delta(q,\tau\mid p,t) \frac{\delta}{\delta a^*(p,t)}\,dp\,dq\,dt\,d\tau,\qquad T\text{-symbol}=e^{\widehat{\mathcal D}}F $$

### Normal symbol of interaction-picture evolution
Appears in: B-100

Applying Wick's theorem term by term to the Dyson series gives the normal symbol of the interaction-picture propagator as e^{D} acting on the exponential of the integrated interaction symbol.

$$ N_{\hat U,1}(t,\tau\mid a^*,a)=e^{\widehat{\mathcal D}} \exp\!\left[i\int_\tau^t V(s\mid a^*,a)ds\right] $$

### Contractions acting on an exponential source
Appears in: B-101

An exponential of linear Grassmann sources is an eigenfunction of the contraction operator, picking up the factor e^{f* Delta f}.

$$ e^{\widehat{\mathcal D}}e^{f^*a+a^*f}=e^{f^*\Delta f+a^*f+f^*a} $$

### Determinant identity for a quadratic symbol
Appears in: B-102

Contracting a fermionic quadratic exponential gives det(I + K Delta) times a quadratic exponential with the dressed kernel (K^{-1} + Delta)^{-1}.

$$ e^{\widehat{\mathcal D}}e^{-a^*\hat Ka} =\det(\hat I+\hat K\hat\Delta)\exp[-a^*(\hat K^{-1}+\hat\Delta)^{-1}a] $$

### Auxiliary-field decoupling of a quartic term
Appears in: B-103

Factorizing the quartic kernel through an auxiliary label lets a commuting Gaussian integral over z, z* turn the quartic exponential into a quadratic one.

$$ \int V_1(p_1,q_1\mid s)V_2(s\mid p_2,q_2)ds=iV(p_1,p_2\mid q_1,q_2) $$

### Retarded loops vanish
Appears in: B-104

Closed loops of retarded contractions with a time-local kernel vanish, so the auxiliary determinant equals 1.

$$ \operatorname{Tr}(\hat R\hat\Delta)^n=0,\qquad \det(\hat I+\hat R\hat\Delta)=1 $$

### Auxiliary-field representation of the propagator
Appears in: B-105

The normal symbol of the full propagator is a Gaussian average over the auxiliary field of a quadratic normal-ordered exponential; with no interaction it reduces to the normal symbol of e^{it omega N}.

$$ N_{e^{it\hat H}}=\int\exp\!\left[a^*\bigl(e^{it\hat\omega}(\hat I-\hat T)-\hat I\bigr)a-zz^*\right]\Pi dz\,dz^* $$

### Second-order Dyson elimination of a fast oscillating coupling
Appears in: G-102

To second order a fast oscillating coupling leaves a secular term linear in t proportional to [A, A^dagger], which defines the effective Hamiltonian.

$$ \hat U_I\approx \hat I-\frac{ig^2t}{\Delta}[\hat A,\hat A^\dagger],\qquad\hat H_{\rm eff}=\frac{\hbar g^2}{\Delta}[\hat A,\hat A^\dagger]. $$

## Nonlinear exponentials (number shift, Kerr, cats)

### Number commutators and free rotation of the field
Appears in: G-004 (sign convention reused in G-101)

The number operator lowers and raises with a and a^dagger, so a free rotation multiplies a by a phase and carries a coherent state around a circle, the zero-point term adding only a global phase.

$$ [\hat n,\hat a]=-\hat a,\quad[\hat n,\hat a^\dagger]=\hat a^\dagger,\quad e^{i\omega t \hat n}\hat ae^{-i\omega t \hat n}=e^{-i\omega t}\hat a,\quad e^{-i\omega t(\hat n+\tfrac12\hat I)}|\alpha\rangle=e^{-i\omega t/2}|\alpha e^{-i\omega t}\rangle. $$

### Number-shift rules for functions of the number operator
Appears in: G-006 (used in G-079, G-080)

Moving a through a function of the number operator shifts the argument of the function by one.

$$ \hat a f(\hat n)=f(\hat n+\hat I)\hat a,\quad f(\hat n)\hat a=\hat a f(\hat n-\hat I) $$

### No Hermitian phase operator canonically conjugate to number
Appears in: G-011

A Hermitian phase operator canonically conjugate to the number operator cannot exist, because its diagonal number-basis elements would have to satisfy 0 = i.

$$ (m-n)\langle m|\hat\phi|n\rangle=i\delta_{mn} $$

### Susskind-Glogower one-sided shift
Appears in: G-012

The Susskind-Glogower operator shifts the number down by one with unit weight; it is an isometry from one side only, missing exactly the vacuum projector on the other.

$$ \hat E=(\hat n+\hat I)^{-1/2}\hat a=\sum_{m\ge0}|m\rangle\langle m+1|,\quad \hat E\hat E^\dagger=\hat I,\quad \hat E^\dagger \hat E=\hat I-|0\rangle\langle0|. $$

### Exact resonant Jaynes-Cummings propagator
Appears in: G-038

On resonance the Jaynes-Cummings propagator rotates each doublet {|e,m>, |g,m+1>} by the angle lambda t sqrt(m+1), which the operator-valued matrix with its square-root number functions encodes.

$$ e^{-i\lambda t(\hat a\hat\sigma_++\hat a^\dagger\hat\sigma_-)}=\begin{pmatrix}\cos(\lambda t\sqrt{\hat n+\hat I})&-i\hat a\,\dfrac{\sin(\lambda t\sqrt{\hat n})}{\sqrt{\hat n}}\\-i\hat a^\dagger\,\dfrac{\sin(\lambda t\sqrt{\hat n+\hat I})}{\sqrt{\hat n+\hat I}}&\cos(\lambda t\sqrt{\hat n})\end{pmatrix}. $$

### Dispersive Hamiltonian as conditional number rotation
Appears in: G-040

In the dispersive limit the atom picks a number-dependent rotation of the field: e^{i chi t n} in the ground state and e^{-i chi t(n+1)} in the excited state.

$$ \hat H_{\rm eff}=\hbar\chi\bigl(|e\rangle\langle e|+\hat n\hat\sigma_3\bigr),\quad\hat U_{\rm eff}=|g\rangle\langle g|\,e^{i\chi t\hat n}+|e\rangle\langle e|\,e^{-i\chi t(\hat n+\hat I)}. $$

### Raman number-controlled atomic rotation
Appears in: G-041

The Raman coupling lambda n sigma_1 rotates the atom by an angle proportional to the photon number, equivalently it rotates the field in opposite senses on the two sigma_1 eigenstates.

$$ \hat U_I(t)=\cos(\lambda t\hat n)\hat I_A-i\sin(\lambda t\hat n)\hat\sigma_1=|+\rangle_x\langle+|e^{-i\lambda t\hat n}+|-\rangle_x\langle-|e^{i\lambda t\hat n}. $$

### Phase shifter and parity exponential
Appears in: G-050, G-076 (Kerr at $t=\pi/K$ equals parity)

A phase shifter rotates a coherent amplitude; at angle pi it is the parity operator, which squares to the identity and maps |alpha> to |-alpha>.

$$ \hat R(\phi)=e^{i\phi \hat n},\quad \hat R(\phi)|\alpha\rangle=|\alpha e^{i\phi}\rangle,\quad\hat\Pi=e^{i\pi \hat n}=(-1)^{\hat n},\quad\hat\Pi^2=\hat I,\quad\hat\Pi|\alpha\rangle=|-\alpha\rangle. $$

### Self-Kerr evolution of a coherent state
Appears in: G-076

Self-Kerr evolution multiplies each Fock amplitude by e^{-i K t m^2}; at t = pi/K this is the parity operator.

$$ \hat U_K(t)=e^{-iKt \hat n^2},\quad \hat U_K(t)|\alpha\rangle=e^{-|\alpha|^2/2}\sum_{m\ge0}\frac{\alpha^m}{\sqrt{m!}}e^{-iKtm^2}|m\rangle,\quad \hat U_K(\pi/K)=\hat\Pi. $$

### Half-revival cat identity
Appears in: G-077

At the half revival the Kerr phase is 1 on even and -i on odd photon numbers, so the operator is a combination of identity and parity and turns a coherent state into a cat.

$$ e^{-i\pi \hat n^2/2}=\frac{e^{-i\pi/4}}{\sqrt2}(\hat I+i\hat\Pi),\quad e^{-i\pi \hat n^2/2}|\alpha\rangle=\frac{e^{-i\pi/4}}{\sqrt2}(|\alpha\rangle+i|-\alpha\rangle). $$

### Normally ordered Kerr differs by a rotation
Appears in: G-078

The normally ordered Kerr generator n(n-1) differs from n^2 by a linear term, so its evolution differs by a free rotation and is periodic with period pi/K.

$$ \hat a^{\dagger2}\hat a^2=\hat n(\hat n-\hat I),\quad e^{-iKt \hat n(\hat n-\hat I)}=e^{-iKt \hat n^2}e^{iKt\hat n},\quad e^{-i\pi\hat n(\hat n-\hat I)}=\hat I. $$

### Number-function conjugation and Kerr Heisenberg solution
Appears in: G-079

Conjugating a by any function of n multiplies it by the phase of the difference f(n) - f(n-1); for the normally ordered Kerr this gives a(t) = a e^{-2iKt(n-1)} = e^{-2iKtn} a.

$$ e^{itf(\hat n)}\hat ae^{-itf(\hat n)}=\hat a\,e^{-it[f(\hat n)-f(\hat n-\hat I)]},\qquad\hat a(t)=\hat a\,e^{-2iKt(\hat n-\hat I)}=e^{-2iKt\hat n}\hat a\quad(f=K\hat n(\hat n-\hat I)). $$

### Kerr coherent-state moments and quadrature variance
Appears in: G-080

Under the normally ordered Kerr Hamiltonian a coherent state's moments <a^m> acquire the phase e^{-iKt m(m-1)} and the factor exp(N(e^{-2imKt} - 1)), and quadrature variances follow from the low moments.

$$ \langle \hat a^m\rangle_t=\alpha^m e^{-iKtm(m-1)}\exp\{N(e^{-2imKt}-1)\},\quad\operatorname{Var}\hat X(\vartheta)=\frac14\{1+2[\langle \hat n\rangle-|\langle \hat a\rangle|^2]+2\operatorname{Re}[e^{-2i\vartheta}(\langle \hat a^2\rangle-\langle \hat a\rangle^2)]\}. $$

Conditions: normally ordered Kerr generator $K\hat n(\hat n-\hat I)$, $N=|\alpha|^2$.

### Cross-Kerr conditional rotation
Appears in: G-081

Cross-Kerr coupling rotates mode a by an angle proportional to the photon number in mode b.

$$ \hat U_{\rm CK}=e^{-iKt\hat n_a\hat n_b},\quad \hat U_{\rm CK}|\alpha\rangle_a|m\rangle_b=|\alpha e^{-iKtm}\rangle_a|m\rangle_b,\quad\hat U_{\rm CK}^\dagger \hat a\hat U_{\rm CK}=\hat a e^{-iKt\hat n_b}. $$

### Coherent-state cat: normalization and a^2 eigenvalue
Appears in: G-082, G-085 (parity projection prepares the even/odd cat)

A superposition of |alpha> and |-alpha> with relative phase Phi has squared norm 2(1 + e^{-2|alpha|^2} cos Phi) and is an eigenstate of a^2 with eigenvalue alpha^2.

$$ |C_\Phi\rangle=\frac{|\alpha\rangle+e^{i\Phi}|-\alpha\rangle}{\sqrt{2(1+e^{-2|\alpha|^2}\cos\Phi)}},\qquad \hat a^2|C_\Phi\rangle=\alpha^2|C_\Phi\rangle. $$

### Cat photon-number distribution
Appears in: G-083

The photon statistics of a cat are the Poisson weights times the interference factor |1 + e^{i Phi}(-1)^m|^2 over the normalization.

$$ |C_\pm\rangle=\frac{(\hat I\pm\hat\Pi)|\alpha\rangle}{\sqrt{2(1\pm e^{-2N})}},\quad P_m=\frac{e^{-N}N^m}{m!}\frac{|1+e^{i\Phi}(-1)^m|^2}{2(1+e^{-2N}\cos\Phi)}. $$

### No-jump number exponential on a coherent state
Appears in: G-096

A no-jump number exponential shrinks a coherent amplitude by e^{-s} and multiplies the state by e^{-|alpha|^2(1 - e^{-2s})/2}.

$$ e^{-s\hat n}|\alpha\rangle=\exp[-\tfrac12|\alpha|^2(1-e^{-2s})]|\alpha e^{-s}\rangle,\qquad\|e^{-s\hat n}|\alpha\rangle\|^2=e^{-|\alpha|^2(1-e^{-2s})}. $$

### Effective Kerr Hamiltonian from an off-resonant parametric coupling
Appears in: G-103

For A = i a^dagger^2 b the commutator [A, A^dagger] is a polynomial in the two number operators, so the effective Hamiltonian is a sum of self-Kerr, cross-Kerr and pump-rotation terms.

$$ [\hat A,\hat A^\dagger]=\hat n_a(\hat n_a-\hat I)-(4\hat n_a+2\hat I)\hat n_b,\qquad\hat A=i\hat a^{\dagger2}\hat b. $$

## Fock matrix elements and generating functions

### Fock basis and ladder action
Appears in: P-005; G-005; B-004

The vacuum is killed by the annihilator, number states are normalized powers of the creator on the vacuum, and the ladder operators move one step with the square-root factors.

$$ \hat a|0\rangle=0,\quad |n\rangle=\frac{(\hat a^\dagger)^n}{\sqrt{n!}}|0\rangle,\quad \hat a|n\rangle=\sqrt n|n-1\rangle,\quad \hat a^\dagger|n\rangle=\sqrt{n+1}|n+1\rangle $$

**Also stated in G-005 (Fock ladder action and resolution of identity):**

The ladder operators move one step on the number basis with square-root weights, the n-photon state is the normalized n-th power of a^dagger on the vacuum, and the number states resolve the identity.

$$ \hat a|n\rangle=\sqrt n|n-1\rangle,\quad \hat a^\dagger|n\rangle=\sqrt{n+1}|n+1\rangle,\quad |n\rangle=\frac{(\hat a^\dagger)^n}{\sqrt{n!}}|0\rangle,\quad\sum_{n\ge0}|n\rangle\langle n|=\hat I. $$

**Also stated in B-004 (Fock vector from the vacuum):**

Any Fock vector is a sum over particle number of creators applied to the vacuum, weighted by a symmetric (bosons) or antisymmetric (fermions) kernel; the 1/sqrt(n!) makes the squared norm the plain sum of integrals of |K_n|^2.

$$ \Phi=\sum_{n\ge0}\frac1{\sqrt{n!}}\int K_n(x_1,\ldots,x_n)\hat a^*(x_1)\cdots\hat a^*(x_n)\,d^nx\,\Omega $$

### Coordinate to Bargmann kernel
Appears in: P-013

The kernel that takes a position wavefunction to its Bargmann function is a Gaussian in $z$ and $q$.

$$ K(z,q)=(\pi\hbar)^{-1/4}\exp\left(-\frac{z^2}{2}+\sqrt{\frac2\hbar}zq-\frac{q^2}{2\hbar}\right) $$

### Displacement matrix elements in associated Laguerre form
Appears in: P-018, P2-052, P2-011 (transition probability $W_{mn}=|\langle m|\hat D(\gamma)|n\rangle|^2$)

Fock matrix elements of a displacement are associated Laguerre polynomials times a Gaussian and a power of the amplitude.

$$ \langle m|\hat D(\gamma)|n\rangle=e^{-|\gamma|^2/2}\begin{cases}\sqrt{n!/m!}\,\gamma^{m-n}L_n^{(m-n)}(|\gamma|^2),&m\ge n,\\ \sqrt{m!/n!}\,(-\bar\gamma)^{n-m}L_m^{(n-m)}(|\gamma|^2),&n\ge m\end{cases} $$

### Coherent generating kernel of the displacement
Appears in: P-018, P2-051 (same identity: $\sum_{mn}D_{mn}\bar\alpha^m\beta^n/\sqrt{m!n!}$ is $e^{(|\alpha|^2+|\beta|^2)/2}\langle\alpha|\hat D|\beta\rangle$)

The unnormalized coherent matrix element of a displacement is a single Gaussian exponential in the two coherent labels and $\gamma$.

$$ e^{(|\alpha|^2+|\beta|^2)/2}\langle\alpha|\hat D(\gamma)|\beta\rangle=e^{-|\gamma|^2/2+\bar\alpha\beta+\bar\alpha\gamma-\beta\bar\gamma} $$

### Displacement matrix elements as a finite sum
Appears in: P2-051

Each displacement matrix element is a finite sum over the number of contractions, with all phases explicit.

$$ D_{mn}=e^{-|\gamma|^2/2}\sqrt{m!n!}\sum_{r=0}^{\min(m,n)}\frac{\gamma^{m-r}(-\bar\gamma)^{n-r}}{r!(m-r)!(n-r)!} $$

### Analytic density kernel
Appears in: P2-040

The unnormalized coherent matrix element of a density operator is a generating function whose Taylor coefficients are the Fock matrix elements.

$$ R(\bar\alpha,\beta)=e^{(|\alpha|^2+|\beta|^2)/2}\langle\alpha|\hat\rho|\beta\rangle=\sum_{mn}\rho_{mn}\frac{\bar\alpha^m\beta^n}{\sqrt{m!n!}} $$

### Laguerre generating function in the lower index
Appears in: P2-053

The generating function of $L_n^k$ in the lower index is $(1-t)^{-k-1}e^{-tx/(1-t)}$.

$$ \sum_{n\ge0}t^nL_n^k(x)=(1-t)^{-k-1}e^{-tx/(1-t)},\quad |t|<1 $$

### Laguerre generating function along the diagonal m-n
Appears in: P2-053

Along the diagonal $m-n$ the Laguerre generating function is $(1+t)^me^{-tx}$.

$$ \sum_{n\ge0}t^nL_n^{m-n}(x)=(1+t)^m e^{-tx} $$

### Generating function of displacement transition probabilities (prefactor corrected)
Appears in: P2-054

The double generating function of the displacement transition probabilities is $(1-uv)^{-1}\exp[-|\gamma|^2(1-u)(1-v)/(1-uv)]$.

$$ F(u,v;\gamma)=\sum_{mn}u^mv^n|D_{mn}(\gamma)|^2=\frac1{1-uv}\exp\left[-\frac{|\gamma|^2(1-u)(1-v)}{1-uv}\right] $$

### Displacement matrix elements as Wirtinger derivatives
Appears in: P2-055

Each displacement matrix element is a mixed Wirtinger derivative of $e^{-|\alpha|^2}$, rescaled; transitions are classically allowed where the two number circles intersect.

$$ D_{mn}(\alpha)=\frac{e^{|\alpha|^2/2}}{\sqrt{m!n!}}(-\partial_{\bar\alpha})^m(\partial_\alpha)^n e^{-|\alpha|^2},\qquad |\sqrt n-|\gamma||<\sqrt m<\sqrt n+|\gamma|\ \text{(semiclassical allowed region)} $$

### Independent-mode CCR and multimode Fock states
Appears in: G-008

Different modes commute, each mode has its own unit commutator, and multimode number states are products of normalized creation powers on the joint vacuum.

$$ [\hat a_j,\hat a_k^\dagger]=\delta_{jk}\hat I,\quad[\hat a_j,\hat a_k]=[\hat a_j^\dagger,\hat a_k^\dagger]=0,\quad |\{n_j\}\rangle=\prod_j\frac{(\hat a_j^\dagger)^{n_j}}{\sqrt{n_j!}}|\{0\}\rangle. $$

### Fock matrix elements from the coherent diagonal symbol
Appears in: G-024

The Fock matrix elements of an operator are Taylor coefficients of its coherent-state diagonal times e^{|alpha|^2}, with alpha and alpha^* as independent variables.

$$ \langle m|\hat F|n\rangle=\frac1{\sqrt{m!n!}}\left.\partial_{\alpha^*}^{m}\partial_\alpha^n\left(e^{|\alpha|^2}\langle\alpha|\hat F|\alpha\rangle\right)\right|_{\alpha=\alpha^*=0}. $$

### Multiphoton Jaynes-Cummings couplings
Appears in: G-042

The two-photon and two-mode Jaynes-Cummings couplings connect |e,m> to |g,m+2> and |e,m,k> to |g,m+1,k+1> with the products of ladder square roots.

$$ \langle g,m+2|\hat a^{\dagger2}\hat\sigma_-|e,m\rangle=\sqrt{(m+1)(m+2)},\qquad\langle g,m+1,k+1|\hat a^\dagger\hat b^\dagger\hat\sigma_-|e,m,k\rangle=\sqrt{(m+1)(k+1)}. $$

### Fock coefficients of a displaced squeezed vacuum (Hermite form)
Appears in: G-062

The Fock amplitudes of a displaced squeezed vacuum are Hermite polynomials in gamma over sqrt(e^{i theta} sinh 2r), times a Gaussian prefactor.

$$ \langle m|\alpha,\xi\rangle=\frac{e^{-|\alpha|^2/2-\alpha^{*2}e^{i\theta}\tanh r/2}}{\sqrt{m!\cosh r}}\left(\frac{e^{i\theta}\tanh r}{2}\right)^{m/2}H_m\!\left(\frac{\gamma}{\sqrt{e^{i\theta}\sinh2r}}\right). $$

## Infinite modes and domains

### Smeared CCR and CAR
Appears in: B-003

Smeared bosonic annihilators and creators commute to the bilinear pairing of their test functions, and smeared fermionic ones anticommute to the same pairing; like-type pairs commute (bosons) or anticommute (fermions) to zero.

$$ [\hat a(f),\hat a^*(g)]=(f,g^*)\hat I,\qquad \{\hat a(f),\hat a^*(g)\}=(f,g^*)\hat I $$

Conditions: bilinear pairing fg = (f, g*); like-type (anti)commutators vanish.

### Bosonic quadratures
Appears in: B-006

The smeared position and momentum quadratures are the Hermitian combinations of a smeared field and its adjoint, normalized so that [Q, P] = i for a normalized real test function.

$$ \hat Q(f)=(\hat a(f)+\hat a^*(f^*))/\sqrt2,\qquad \hat P(f)=(\hat a(f)-\hat a^*(f^*))/(i\sqrt2) $$

### Weighted sufficient conditions for a quadratic normal form
Appears in: B-078

When C is multiplication by c(x), finiteness of these two weighted integrals suffices for the formal quadratic expression to be reduced to a quadratic normal form.

$$ \int\frac{|f(x)|^2}{|c(x)-z|}dx<\infty,\quad \iint\frac{|A(x,y)|^2}{|c(x)+c(y)-z|}dxdy<\infty $$

### Square-root perturbation stays in the ideal
Appears in: B-087

The correction X = (C^2 + K)^{1/2} - C solves X^2 + CX + XC = K, which the source uses to show X is as small (trace class, Hilbert-Schmidt) as K.

$$ \hat X=(\hat C^2+\hat K)^{1/2}-\hat C,\qquad \hat X^2+\hat C\hat X+\hat X\hat C=\hat K $$

### Scattering unitary in the Thirring model
Appears in: B-092

The Thirring model's scattering operator is a phase times the exponential of the product of the two zero-momentum charges.

$$ \hat S=\theta e^{-\pi igF(0)\hat Q_1'(0)\hat Q_2'(0)} $$

### Calogero degeneracy generating function
Appears in: P2-016

The degeneracies of the Calogero model after removing the center of mass are generated by $\prod_{j=2}^N(1-t^j)^{-1}$.

$$ \prod_{j=2}^N(1-t^j)^{-1} $$

### Canonical commutator of position and momentum
Appears in: C-001

Momentum and position along the same axis fail to commute by a fixed imaginary constant; along different axes they commute.

$$ [\hat P_j,\hat Q_k]=-i\hbar\delta_{jk}\hat I $$

Conditions: on the product domain; real phase-space labels.

### Periodic-orbit determinant in regularized traces
Appears in: C-036

Energy-localized traces of the propagator are asymptotically sums over periodic orbits, each weighted by |det(I - P_gamma)|^{-1/2} and a Maslov phase.

$$ \operatorname{Tr}[\chi(\hat H)\hat A\chi(\hat H)g((E\hat I-\hat H)/\hbar)]\ \sim\ \ldots\,|\det(I_{d_\gamma}-P_\gamma)|^{-1/2}e^{i\pi(\cdot)} $$

### Segal fields
Appears in: C-057

Smeared annihilation and creation operators satisfy the canonical commutator with the one-particle inner product, and the field operator is their real combination.

$$ [\hat a(f),\hat a^\dagger(g)]=\langle f,g\rangle\hat I,\qquad \hat\Phi(f)=\frac{\hat a(f)+\hat a^\dagger(f)}{\sqrt2} $$

### Infinite-mode displacement, exponential vector, product, overlap
Appears in: C-059 (one-mode content identical to C-003, C-006, C-007, C-008)

On Fock space over any one-particle space, the Weyl operator factors in normal order, composes with the phase e^{-i Im<f,g>}, and coherent vectors overlap as e^{-|f-g|^2/2 + i Im<f,g>}; occupation numbers are Poisson with mean ||f||^2.

$$ \hat T(f)=e^{-\|f\|^2/2}e^{\hat a^\dagger(f)}e^{-\hat a(f)},\quad \hat T(f)\hat T(g)=e^{-i\Im\langle f,g\rangle}\hat T(f+g),\quad \langle\Psi(f),\Psi(g)\rangle=e^{-\|f-g\|^2/2+i\Im\langle f,g\rangle} $$

### Collective broadband commutator
Appears in: G-094

A collective mode summed over a band keeps a bosonic commutator scaled by C = sum |l_k|^2, and its quadratures have commutator iC/2.

$$ [\hat A(t),\hat A^\dagger(t)]=C\hat I,\quad C=\sum_{k\in\mathcal B}|l_k|^2,\quad[\hat X_1^{(c)},\hat X_2^{(c)}]=(iC/2)\hat I. $$

Conditions: $C$ finite.

### Continuum CCR and first-order pair-creation kernel
Appears in: G-105

Continuum modes have a delta-function commutator, and to first order a pump photon splits into a pair constrained by energy and momentum conservation.

$$ [\hat a(\mathbf k),\hat a^\dagger(\mathbf k')]=\delta^{(3)}(\mathbf k-\mathbf k')\hat I,\qquad\int d^3k_s\,d^3k_i\,\delta(\omega_P-\omega_s-\omega_i)\delta^{(3)}(\mathbf k_P-\mathbf k_s-\mathbf k_i)\hat a_s^\dagger(\mathbf k_s)\hat a_i^\dagger(\mathbf k_i)|0\rangle. $$

## Lie-group decompositions and generalized coherent states

### SU(2) disentangling (compact case)
Appears in: X-004 (source records cited: P-025, P-029, C-041, C-048, G-046, P-041 to P-043, P2-025); P-025; C-041; fermion-pair application in C-087

A compact group element generated by $\hat J_\pm$ factors into raising, Cartan, and lowering exponentials with $\zeta=e^{i\phi}\tan r$, on the chart where $\cos r>0$.

$$ e^{\xi \hat J_+-\bar\xi \hat J_-}=e^{\zeta \hat J_+}e^{\log(1+|\zeta|^2)\hat J_0}e^{-\bar\zeta \hat J_-},\qquad \zeta=e^{i\phi}\tan r $$

Conditions: $[\hat J_0,\hat J_\pm]=\pm\hat J_\pm$, $[\hat J_+,\hat J_-]=2\hat J_0$, $\xi=re^{i\phi}$, $\cos r>0$. The chart condition is needed: for $\cos r<0$ the middle factor loses the sign of $\cos r$, which matters in half-integer representations.

**Also stated in P-025 (SU(2) disentangling):**

A compact rotation exponential factors into raising, diagonal and lowering exponentials with $\zeta=e^{i\phi}\tan r$, in either order, on the chart $\cos r>0$.

$$ e^{\xi \hat J_+-\bar\xi \hat J_-}=e^{\zeta \hat J_+}e^{\log(1+|\zeta|^2)\hat J_0}e^{-\bar\zeta \hat J_-}=e^{-\bar\zeta \hat J_-}e^{-\log(1+|\zeta|^2)\hat J_0}e^{\zeta \hat J_+},\quad \xi=re^{i\phi},\ \zeta=e^{i\phi}\tan r $$

Conditions: $\cos r>0$.

**Also stated in C-041 (SU(2) normal factorization and spin coherent state):**

The SU(2) displacement factors into a raising exponential, a J_3 exponential and a lowering exponential with tan r, and acting on the lowest weight it gives the normalized spin coherent state.

$$ e^{\xi\hat J_+-\bar\xi\hat J_-}=e^{\zeta\hat J_+}e^{-2\log(\cos r)\hat J_3}e^{-\bar\zeta\hat J_-},\ \zeta=\frac{\xi}{r}\tan r;\qquad \hat G(\xi)|j,-j\rangle=(1+|\zeta|^2)^{-j}e^{\zeta\hat J_+}|j,-j\rangle $$

Conditions: cos r != 0.

### SU(1,1) disentangling (noncompact case)
Appears in: X-004; P-029 (also referenced in P2-058)

The noncompact analogue uses $\zeta=e^{i\phi}\tanh r$ and $\log(1-|\zeta|^2)$, and needs no chart restriction.

$$ e^{\xi \hat K_+-\bar\xi \hat K_-}=e^{\zeta \hat K_+}e^{\log(1-|\zeta|^2)\hat K_0}e^{-\bar\zeta \hat K_-},\qquad \zeta=e^{i\phi}\tanh r $$

Conditions: $[\hat K_0,\hat K_\pm]=\pm\hat K_\pm$, $[\hat K_-,\hat K_+]=2\hat K_0$; no chart restriction.

**Also stated in P-029 (SU(1,1) disentangling (middle-exponent sign corrected)):**

A noncompact su(1,1) exponential factors into raising, diagonal and lowering exponentials with $\zeta=e^{i\phi}\tanh r$; the diagonal exponent is $+\log(1-|\zeta|^2)\hat K_0$ in normal order and the opposite in antinormal order.

$$ e^{\xi \hat K_+-\bar\xi \hat K_-}=e^{\zeta \hat K_+}e^{\log(1-|\zeta|^2)\hat K_0}e^{-\bar\zeta \hat K_-}=e^{-\bar\zeta \hat K_-}e^{-\log(1-|\zeta|^2)\hat K_0}e^{\zeta \hat K_+},\quad \zeta=e^{i\phi}\tanh r $$

### su(2) commutators, ladder action and holomorphic realization
Appears in: P-024; C-038

The spin generators obey the su(2) commutators, act on $|j,m\rangle$ with the standard square-root ladder factors, and are realized on polynomials of degree at most $2j$ by first-order differential operators.

$$ [\hat J_0,\hat J_\pm]=\pm \hat J_\pm,\ [\hat J_+,\hat J_-]=2\hat J_0,\ \hat J_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}|j,m\pm1\rangle;\quad \hat J_-=\partial_z,\ \hat J_+=-z^2\partial_z+2jz,\ \hat J_0=z\partial_z-j $$

**Also stated in C-038 (SU(2) commutators, ladder matrices, Casimir):**

The spin-j ladder operators satisfy the su(2) relations, raise and lower m with the standard square-root factors, and give the Casimir j(j+1).

$$ [\hat J_3,\hat J_\pm]=\pm\hat J_\pm,\ [\hat J_+,\hat J_-]=2\hat J_3,\ \hat J_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}|j,m\pm1\rangle,\ \hat J^2=j(j+1)\hat I $$

### SU(2) coherent state expansion
Appears in: P-026

The spin coherent state is the normalized exponential of $\zeta\hat J_+$ on the lowest weight, with binomial square-root amplitudes.

$$ |\zeta;j\rangle=(1+|\zeta|^2)^{-j}e^{\zeta \hat J_+}|j,-j\rangle=(1+|\zeta|^2)^{-j}\sum_{n=0}^{2j}\sqrt{\binom{2j}{n}}\zeta^n|j,-j+n\rangle $$

### SU(2) coherent overlap and annihilator
Appears in: P-026; P-026; C-042

Spin coherent states overlap by $(1+\bar\zeta\eta)^{2j}$ over the normalizations, and each is annihilated by $\hat J_-+2\zeta\hat J_0-\zeta^2\hat J_+$.

$$ \langle\zeta;j|\eta;j\rangle=\frac{(1+\bar\zeta\eta)^{2j}}{(1+|\zeta|^2)^j(1+|\eta|^2)^j},\qquad (\hat J_-+2\zeta \hat J_0-\zeta^2\hat J_+)|\zeta;j\rangle=0 $$

**Also stated in P-026 (SU(2) resolution of identity):**

Spin coherent projectors integrated with the invariant sphere measure resolve the identity.

$$ \frac{2j+1}{\pi}\int\frac{d^2\zeta}{(1+|\zeta|^2)^2}|\zeta;j\rangle\langle\zeta;j|=\hat I $$

**Also stated in C-042 (Spin coherent overlap and resolution of identity):**

Spin coherent states overlap as (1 + conj(eta) zeta)^{2j} over the norms, with modulus set by the angle between the two directions, and resolve the identity with the measure (2j+1)/pi (1+|zeta|^2)^{-2}.

$$ \langle\eta|\zeta\rangle=\frac{(1+\bar\eta\zeta)^{2j}}{(1+|\eta|^2)^j(1+|\zeta|^2)^j},\quad |\langle n'|n\rangle|^2=\left(\frac{1+n\cdot n'}2\right)^{2j},\quad \frac{2j+1}{\pi}\int\frac{|\zeta\rangle\langle\zeta|}{(1+|\zeta|^2)^2}d^2\zeta=\hat I $$

### su(1,1) algebra, Casimir, discrete-series ladder and holomorphic realization
Appears in: P-028; C-045

The su(1,1) generators close with $[\hat K_-,\hat K_+]=2\hat K_0$, the Casimir is $k(k-1)$ in the positive discrete series, and the series has explicit ladder factors and a holomorphic realization.

$$ [\hat K_0,\hat K_\pm]=\pm \hat K_\pm,\ [\hat K_-,\hat K_+]=2\hat K_0,\ \hat C=\hat K_0^2-\tfrac12(\hat K_+\hat K_-+\hat K_-\hat K_+)=k(k-1)\hat I,\ \hat K_+|k,n\rangle=\sqrt{(n+1)(2k+n)}|k,n+1\rangle,\ \hat K_-=\partial_z,\ \hat K_+=z^2\partial_z+2kz,\ \hat K_0=z\partial_z+k $$

**Also stated in C-045 (SU(1,1) algebra and Bergman realization):**

The su(1,1) generators close with [K-, K+] = 2 K0 and Casimir k(k-1), and act on holomorphic functions on the disk as first-order differential operators.

$$ [\hat K_-,\hat K_+]=2\hat K_0,\ [\hat K_0,\hat K_\pm]=\pm\hat K_\pm,\ \hat C=\hat K_0^2-\tfrac12(\hat K_+\hat K_-+\hat K_-\hat K_+)=k(k-1)\hat I;\quad \hat K_-=\partial_z,\ \hat K_+=2kz+z^2\partial_z,\ \hat K_0=k+z\partial_z $$

### SU(1,1) coherent state expansion and normalization
Appears in: P-030; P-030; C-047; overlap modulus reused in C-097

The su(1,1) coherent state is the normalized exponential of $\zeta\hat K_+$ on the lowest weight, with Pochhammer square-root amplitudes, and it is normalized for $|\zeta|<1$.

$$ |\zeta;k\rangle=(1-|\zeta|^2)^ke^{\zeta \hat K_+}|k,0\rangle=(1-|\zeta|^2)^k\sum_{n\ge0}\sqrt{\frac{\Gamma(n+2k)}{n!\Gamma(2k)}}\zeta^n|k,n\rangle $$

**Also stated in P-030 (SU(1,1) coherent overlap and annihilator):**

Two su(1,1) coherent states overlap by $(1-\bar\zeta\eta)^{-2k}$ times the normalizations, and each is annihilated by $\hat K_--2\zeta\hat K_0+\zeta^2\hat K_+$.

$$ \langle\zeta;k|\eta;k\rangle=\frac{(1-|\zeta|^2)^k(1-|\eta|^2)^k}{(1-\bar\zeta\eta)^{2k}},\qquad (\hat K_--2\zeta \hat K_0+\zeta^2\hat K_+)|\zeta;k\rangle=0 $$

**Also stated in C-047 (SU(1,1) coherent series and overlap):**

The SU(1,1) coherent state is the normalized exponential of the raising operator on the lowest-weight vector; its coefficients are square roots of Gamma-function ratios and its overlap is a power of (1 - conj(eta) zeta).

$$ |\zeta;k\rangle=(1-|\zeta|^2)^k\sum_m\sqrt{\frac{\Gamma(2k+m)}{\Gamma(2k)m!}}\zeta^m|k,m\rangle,\qquad \langle\eta;k|\zeta;k\rangle=\frac{(1-|\eta|^2)^k(1-|\zeta|^2)^k}{(1-\bar\eta\zeta)^{2k}} $$

Conditions: |zeta| < 1; resolution measure (2k-1)/pi needs k > 1/2.

### SU(1,1) resolution of identity
Appears in: P-030

su(1,1) coherent projectors integrated with the invariant disk measure resolve the identity only for $k>1/2$.

$$ \frac{2k-1}{\pi}\int_{|\zeta|<1}\frac{d^2\zeta}{(1-|\zeta|^2)^2}|\zeta;k\rangle\langle\zeta;k|=\hat I,\quad k>1/2 $$

### Root-space coherent states and displacement-like generator
Appears in: P-047

For a general semisimple group, coherent states are exponentials of positive-root generators on a lowest weight, and the unitary displacement uses positive and negative roots.

$$ |\boldsymbol\zeta\rangle=\mathcal N(\boldsymbol\zeta)\exp\left(\sum_{\alpha>0}\zeta_\alpha \hat E_\alpha\right)|-\lambda\rangle,\qquad \hat D(\boldsymbol\xi)=\exp\left(\sum_{\alpha>0}\xi_\alpha \hat E_\alpha-\bar\xi_\alpha \hat E_{-\alpha}\right) $$

### Compact Kahler kernel
Appears in: P2-007

On the sphere the reproducing kernel is $(1+z\bar w)^{1/h}$ with $1/h$ an integer.

$$ K(z,\bar w)=(1+z\bar w)^{1/h},\quad h^{-1}\in\mathbb Z $$

### Jacobians of classical-domain transformations
Appears in: P2-056

The complex Jacobian of the fractional-linear action on a type-I domain is $\det(B^Tz+D)^{-(p+q)}$.

$$ \det(B^Tz+D)^{-(p+q)} $$

### Completeness criterion for discrete SU(1,1) subsystems
Appears in: P2-004

A discrete SU(1,1) coherent subsystem is complete when the fundamental area is below $\pi/(2k-1)$ and incomplete for large $k$.

$$ S_\Gamma<\frac{\pi}{2k-1}\ \Rightarrow\ \text{complete};\qquad k>m_0^++\tfrac12\ \Rightarrow\ \text{incomplete} $$

### Clock and shift matrices on the torus
Appears in: C-037

On the quantum torus with 2 pi hbar N = 1, translations become N-by-N clock and shift matrices of order N that commute up to a root of unity.

$$ 2\pi\hbar N=1,\quad \hat X^N=\hat Z^N=\hat I,\quad \hat Z\hat X=\omega\hat X\hat Z,\ \omega=e^{2\pi i/N} $$

### Gauss decomposition of SL(2,C)
Appears in: C-040

A unimodular 2-by-2 matrix with nonzero lower-right entry factors uniquely into upper-triangular, diagonal and lower-triangular pieces.

$$ g=\begin{pmatrix}1&\beta/\delta\\0&1\end{pmatrix}\begin{pmatrix}\delta^{-1}&0\\0&\delta\end{pmatrix}\begin{pmatrix}1&0\\\gamma/\delta&1\end{pmatrix} $$

Conditions: alpha delta - beta gamma = 1, delta != 0.

### Spin-to-oscillator contraction
Appears in: C-044

As j grows, the spin coherent kernel at scale zeta/sqrt(2j) tends to the bosonic coherent kernel.

$$ \text{spin kernel at scale }\zeta/\sqrt{2j}\ \to\ e^{\bar\eta\zeta-|\eta|^2/2}\quad (j\to\infty) $$

### Singular oscillator lowest weight
Appears in: C-052

For the oscillator with an inverse-square potential on the half line, the lowest-weight state is x^s e^{-x^2/2} with s fixed by g, and its Bargmann index is k = (2s+1)/4.

$$ \psi\propto x^se^{-x^2/2},\quad s=\tfrac12+\sqrt{\tfrac14+g^2},\quad k=\frac{2s+1}{4} $$

Conditions: H = (-d^2 + x^2 + g^2/x^2)/2 on L^2(R+), chosen self-adjoint realization.

### Berezin symbol product on the disk
Appears in: C-097

On the disk, the covariant symbol of a product is an integral of the two off-diagonal symbols against the squared overlap of SU(1,1) coherent states.

$$ (\hat A\hat B)_c(z,\bar z)=\int_DA_c(z,\bar w)B_c(w,\bar z)|\langle\psi_z,\psi_w\rangle|^2d\nu_n(w),\qquad |\langle\psi_z,\psi_w\rangle|^2=\left[\frac{(1-|z|^2)(1-|w|^2)}{|1-\bar zw|^2}\right]^n $$

## Records with no formula

These atlas records state domains, theorems in prose, or coverage notes.

- X-008: a five-step derivation route (commutator closure, ordered ansatz, Riccati ODEs, initial conditions and charts, limits for cubic and quartic generators) plus a statement of what was checked when the atlas was assembled. No formula.
- B-001: Fock spaces and dimensions (2^N fermionic, infinite bosonic).
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
- P-023: Gauss decomposition on the big cell of SL(2,C), singular at a vanishing pivot (prose; the 2x2 instance is P-029).
- P-027: spin-to-oscillator contraction $\hat J_\pm\sim\sqrt{2j}(\hat a^\dagger,\hat a)$, $\zeta=\alpha/\sqrt{2j}$ (asymptotic prose; consistent with P-026, since $(1+\bar\zeta\eta)^{2j}\to e^{\bar\alpha\beta}$).
- P-033: principal and Lorentz series scope boundary (prose).
- P-040: Siegel disk geometry and the $k>2N$ integral restriction (prose).
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
- P2-050: coherent-lattice completeness theorem (stated with the lattice entry).
- P2-057: addendum and references coverage record.
- P2-058: source-audit cautions (cross-references to P2-017, P2-023, P2-038, P2-054, P2-048, P-029).
- C-005: analytic core for complex labels (domain statement).
- C-033: quadratic replacement gives O(hbar^{1/2}) and Ehrenfest-time errors; asymptotic.
- C-034: branch prescription for Gaussian determinants (Maslov index), prose.
- C-062: Bogoliubov fluctuation dynamics around a Hartree solution; convergence statement.
- C-063: Dyson/Gevrey/Borel summability statement.
- C-081: Spin double cover and O(2n) versus SO(2n) implementers, prose (consistent with C-087 and with c <-> c^dag being odd).
- C-092: stationary phase is asymptotic, prose.
- C-093: Moyal remainder and functional calculus, prose.
- C-095: Schur uniqueness and square integrability for group coherent states, prose.
- C-096: Cartan polar decomposition, prose.
- G-085: conditional measurement prepares cats; the only formula ($\propto\hat I\pm\hat\Pi$) is stated with the cat parity entry.
