# Atlas section X (common-convention formulas, X-001 to X-009)

Every check below is exact: QF's symbolic second-quantization algebra, exact finite group matrices, Jordan-Wigner fermion matrices (the complete operators), Bargmann-function realizations, exact integrals, or hand derivation. No Fock cutoff is used anywhere. Each `wl` cell assumes the loading lines `Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"]`. The cells are in `drafts/exact/X-work/cells/`, their outputs in `cells-out.txt`; the non-QF exact checks are `exactwl.wl` (output `exactwl-out.txt`) and `x009boson.wl` (output `x009boson-out.txt`).

## Errors found

None of the identities in X-001 to X-009 is wrong in sign, factor, or condition. Two wording points, not errors:

- X-006 says the fermionic vacuum chart "can fail when $U$ is singular", but no $U$ is defined in X-006. It means the canonical block $u$ of X-009 (Perelomov's $U$); the state $|Z\rangle_F$ exists only on the chart where $u$ is invertible.
- X-007 states the bosonic trace condition as $0\le M<I_N$ "equivalently $\|M\|<1$". The equivalence holds only for positive $M$ (which is the setting stated). For a general complex $M$ the trace is $\prod_i\sum_n\lambda_i^n$ over the eigenvalues of $M$, so it converges exactly when the spectral radius is below 1; $\|M\|<1$ is sufficient, not necessary (a nilpotent $M$ of any norm gives trace 1).

---

### Displacement: normal and antinormal factorization
Topic: Displacement
Appears in: X-001 (source records cited: P-002, P-003, P-010, G-015, G-017, B-050, C-006)

The displacement operator splits into a creation exponential times an annihilation exponential (or the reverse), and the only cost of reordering is a Gaussian scalar.

$$ \hat D(\alpha)=e^{\alpha \hat a^\dagger-\bar\alpha \hat a}=e^{-|\alpha|^2/2}e^{\alpha \hat a^\dagger}e^{-\bar\alpha \hat a}=e^{|\alpha|^2/2}e^{-\bar\alpha \hat a}e^{\alpha \hat a^\dagger},\qquad [\hat a,\hat a^\dagger]=\hat I $$

Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{{c01}}
```
Output: `{E^(al*SuperDagger[\[FormalA]1])**E^(-(\[FormalA]1*Conjugate[al]))/E^((al*Conjugate[al])/2), E^((al*Conjugate[al])/2)*E^(-(\[FormalA]1*Conjugate[al]))**E^(al*SuperDagger[\[FormalA]1]), <same two forms again from BosonicExpOrder>}`. Matches: both orderings are the atlas forms verbatim, from the named operator and from ordering the bare exponential.

### Displacement: adjoint action on a
Topic: Displacement
Appears in: X-001

Conjugating the annihilator by the displacement shifts it by the c-number $\alpha$.

$$ \hat D(\alpha)^\dagger \hat a\hat D(\alpha)=\hat a+\alpha\hat I $$

Conditions: the atlas notes that Berezin's direction $\hat U\hat a\hat U^{-1}$ reverses the sign of the parameter; with $\hat U=\hat D(\alpha)$ that direction gives $\hat a-\alpha$.
Audit: Correct, including the Berezin-direction remark. Checked by: QF exact.
QF code (exact):
```wl
{{c02}}
```
Output: `{\[FormalA]1 + al, \[FormalA]1 - al}`. Matches: the first is $\hat D^\dagger\hat a\hat D$, the second the reversed direction $\hat D\hat a\hat D^\dagger$.

### Displacement: composition law
Topic: Displacement
Appears in: X-001

Two displacements compose to the displacement by the sum, up to a phase set by the symplectic area of the two amplitudes.

$$ \hat D(\alpha)\hat D(\beta)=e^{(\alpha\bar\beta-\bar\alpha\beta)/2}\hat D(\alpha+\beta) $$

Audit: Correct (with the opposite phase sign the scalar would not match the QF result). Checked by: QF exact.
QF code (exact):
```wl
{{c03}}
```
Output: `{E^(-(\[FormalA]1*Conjugate[al]) - \[FormalA]1*Conjugate[be] + (-(be*Conjugate[al]) + al*Conjugate[be])/2 + al*SuperDagger[\[FormalA]1] + be*SuperDagger[\[FormalA]1]), 0}`. Matches: the merged exponent carries the scalar $(\alpha\bar\beta-\bar\alpha\beta)/2$ and the difference from the atlas right side is zero.

### Single-mode squeeze: normal-ordered factorization
Topic: Single-mode squeezing
Appears in: X-002 (source records cited: P-028 to P-031, G-051 to G-061, B-051, C-048, C-049)

The squeeze operator factors into a pair-creation exponential, a power of the number operator, and a pair-annihilation exponential, with $\tanh r$ in the outer factors.

$$ \hat S(\xi)=\exp\!\left[\tfrac12\left(\bar\xi \hat a^2-\xi(\hat a^\dagger)^2\right)\right]=e^{-\frac12e^{i\theta}\tanh r\,(\hat a^\dagger)^2}e^{-\log\cosh r\,(\hat N+\tfrac12\hat I)}e^{\frac12e^{-i\theta}\tanh r\,\hat a^2} $$

Conditions: $\xi=re^{i\theta}$, $r\ge0$.
Audit: Correct, including the $\tfrac12$ in the number exponent. Checked by: QF exact.
QF code (exact):
```wl
{{c04}}
```
Output: `{E^(-1/2*(E^(I*Arg[xi])*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]1]*Tanh[Abs[xi]]))**Sech[Abs[xi]]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^((\[FormalA]1**\[FormalA]1*Tanh[Abs[xi]])/(2*E^(I*Arg[xi]))), E^(-1/2*(E^(I*th)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]*Tanh[r]))**Sech[r]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^((\[FormalA]1**\[FormalA]1*Tanh[r])/(2*E^(I*th)))}`. Matches: $\operatorname{sech}r^{\,\hat N+1/2}=e^{-\log\cosh r(\hat N+\frac12)}$, so both the named operator and the ordered bare exponential are the atlas form.

### Single-mode squeeze: adjoint action on a
Topic: Single-mode squeezing
Appears in: X-002

Squeezing mixes the annihilator with the creator, with $\cosh r$ and $\sinh r$ weights.

$$ \hat S^\dagger \hat a\hat S=\hat a\cosh r-e^{i\theta}\hat a^\dagger\sinh r $$

Conditions: the atlas notes that $\hat S^{-1}$ flips the sign of every squeeze parameter.
Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{{c05}}
```
Output: `\[FormalA]1*Cosh[r] - E^(I*th)*Sinh[r]*SuperDagger[\[FormalA]1]`. Matches.

### Single-mode squeezed vacuum
Topic: Single-mode squeezing
Appears in: X-002

The squeezed vacuum is a normalized Gaussian superposition of even photon numbers.

$$ \hat S|0\rangle=(\cosh r)^{-1/2}e^{-\frac12e^{i\theta}\tanh r\,(\hat a^\dagger)^2}|0\rangle $$

Audit: Correct; it is X-006 with $Z=e^{i\theta}\tanh r$, since $(1-\tanh^2r)^{1/4}=(\cosh r)^{-1/2}$. Checked by: QF exact, and hand derivation: in the QF normal form above the annihilation factor acts as $\hat I$ on $|0\rangle$ and $\operatorname{sech}r^{\,\hat N+1/2}|0\rangle=(\cosh r)^{-1/2}|0\rangle$.
QF code (exact):
```wl
{{c06}}
```
Output: `{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 0}`. Matches: QF's closed-form amplitudes $\langle m|\hat S|0\rangle$ equal the formula's for $m=0,\dots,10$ with $r,\theta$ symbolic, and the formula's state has norm exactly 1 (the second entry is norm minus 1, summed in closed form over all $m$).

### Two-mode squeeze: normal-ordered factorization
Topic: Two-mode squeezing
Appears in: X-003 (source records cited: P-032, G-068 to G-075, P-029)

The two-mode squeezer factors the same way as the single-mode one, with the pair $\hat a^\dagger\hat b^\dagger$ in place of $(\hat a^\dagger)^2/2$.

$$ \hat S_2(\xi)=e^{\bar\xi \hat a\hat b-\xi \hat a^\dagger \hat b^\dagger}=e^{-e^{i\theta}\tanh r\,\hat a^\dagger \hat b^\dagger}e^{-\log\cosh r\,(\hat N_a+\hat N_b+\hat I)}e^{e^{-i\theta}\tanh r\,\hat a\hat b} $$

Conditions: closed algebra $\hat K_+=\hat a^\dagger\hat b^\dagger$, $\hat K_-=\hat a\hat b$, $\hat K_0=(\hat N_a+\hat N_b+\hat I)/2$; $\hat N_a-\hat N_b$ conserved; P-032's convention is $\xi\mapsto-\xi$.
Audit: Correct. Checked by: QF exact.
QF code (exact):
```wl
{{c07}}
```
Output: `E^(-(E^(I*th)*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]2]*Tanh[r]))**Sech[r]^(SuperDagger[\[FormalA]1]**\[FormalA]1 + SuperDagger[\[FormalA]2]**\[FormalA]2)**E^((\[FormalA]1**\[FormalA]2*Tanh[r])/E^(I*th))*Sech[r]`. Matches: the trailing `Sech[r]` is the $\hat I$ in the middle exponent.

### Two-mode squeezed vacuum
Topic: Two-mode squeezing
Appears in: X-003

The two-mode squeezed vacuum is a sum of twin-photon-number states with geometric weights.

$$ \hat S_2|0,0\rangle=\operatorname{sech}r\sum_{n\ge0}(-e^{i\theta}\tanh r)^n|n,n\rangle $$

Audit: Correct. Checked by: hand derivation: in the QF normal form above the annihilation factor and $\operatorname{sech}r^{\,\hat N_a+\hat N_b}$ act as $\hat I$ on $|0,0\rangle$, and $(\hat a^\dagger\hat b^\dagger)^n|0,0\rangle=n!\,|n,n\rangle$; and exact WL computation: in the Bargmann realization ($\hat a_i\to\partial_{z_i}$, $\hat a_i^\dagger\to z_i$) the function $\psi(s)=\operatorname{sech}(sr)\exp(-e^{i\theta}\tanh(sr)\,z_1z_2)$ satisfies $\partial_s\psi=(\bar\xi\partial_1\partial_2-\xi z_1z_2)\psi$ with $\psi(0)=1$ exactly, so it is $\hat S_2(s\xi)|0,0\rangle$.
QF code (exact): none. `BosonicMatrixElement` takes multimode Fock indices only for polynomials; `BosonicMatrixElement[{{2, 2}, {0, 0}}, Exp[x a1 ** a2 + y ad1 ** ad2]]` returns unevaluated.

### SU(2) disentangling (compact case)
Topic: Lie-group decompositions and generalized coherent states
Appears in: X-004 (source records cited: P-025, P-029, C-041, C-048, G-046, P-041 to P-043, P2-025)

A compact group element generated by $\hat J_\pm$ factors into raising, Cartan, and lowering exponentials with $\zeta=e^{i\phi}\tan r$, on the chart where $\cos r>0$.

$$ e^{\xi \hat J_+-\bar\xi \hat J_-}=e^{\zeta \hat J_+}e^{\log(1+|\zeta|^2)\hat J_0}e^{-\bar\zeta \hat J_-},\qquad \zeta=e^{i\phi}\tan r $$

Conditions: $[\hat J_0,\hat J_\pm]=\pm\hat J_\pm$, $[\hat J_+,\hat J_-]=2\hat J_0$, $\xi=re^{i\phi}$, $\cos r>0$. The chart condition is needed: for $\cos r<0$ the middle factor loses the sign of $\cos r$, which matters in half-integer representations.
Audit: Correct, including the chart remark. Checked by: exact WL computation: symbolic $r,\phi$ in the spin-1/2 and spin-1 matrices; the difference is zero for $0<r<\pi/2$; for $\pi/2<r<\pi$ the spin-1/2 right side is exactly minus the left side, while spin 1 still matches.
QF code (exact): none. QF has no abstract SU(2) disentangling.

### SU(1,1) disentangling (noncompact case)
Topic: Lie-group decompositions and generalized coherent states
Appears in: X-004

The noncompact analogue uses $\zeta=e^{i\phi}\tanh r$ and $\log(1-|\zeta|^2)$, and needs no chart restriction.

$$ e^{\xi \hat K_+-\bar\xi \hat K_-}=e^{\zeta \hat K_+}e^{\log(1-|\zeta|^2)\hat K_0}e^{-\bar\zeta \hat K_-},\qquad \zeta=e^{i\phi}\tanh r $$

Conditions: $[\hat K_0,\hat K_\pm]=\pm\hat K_\pm$, $[\hat K_-,\hat K_+]=2\hat K_0$; no chart restriction.
Audit: Correct. With $\hat K_+=\hat a^{\dagger2}/2$, $\hat K_-=\hat a^2/2$, $\hat K_0=(\hat N+\tfrac12)/2$ and $\xi\to-\xi$ it is the X-002 factorization; with the two-mode $\hat K$'s it is X-003. Checked by: exact WL computation (the 2x2 non-unitary representation $\hat K_0=\sigma_z/2$, $\hat K_+=\sigma_+$, $\hat K_-=-\sigma_-$ with symbolic $r,\phi$: algebra and identity both exact zero) and QF exact for the single-mode boson realization.
QF code (exact):
```wl
{{c10}}
```
Output: `E^((E^(I*ph)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]1], 2]*Tanh[r])/2)**Sech[r]^(1/2 + SuperDagger[\[FormalA]1]**\[FormalA]1)**E^(-1/2*(\[FormalA]1**\[FormalA]1*Tanh[r])/E^(I*ph))`. Matches: this is $e^{\zeta\hat K_+}(1-|\zeta|^2)^{\hat K_0}e^{-\bar\zeta\hat K_-}$, since $(1-\tanh^2r)^{(\hat N+1/2)/2}=\operatorname{sech}r^{\,\hat N+1/2}$.

### Fermion-pair realization of su(2)
Topic: Fermionic canonical transformations
Appears in: X-004 (marked by the atlas as derived)

Creating or destroying a fermion pair in two modes, together with half the shifted number, closes on the compact su(2) algebra.

$$ \hat J_+=\hat c_1^\dagger \hat c_2^\dagger,\quad \hat J_-=\hat c_2\hat c_1=\hat J_+^\dagger,\quad \hat J_0=(\hat N_1+\hat N_2-\hat I)/2 $$

Audit: Correct; these satisfy the compact algebra. The order $\hat c_2\hat c_1$ is essential, as the atlas says: with $\hat c_1\hat c_2$ (which equals $-\hat J_+^\dagger$) the relation $[\hat J_+,\hat J_-]=2\hat J_0$ fails. Checked by: exact fermion matrices (2 modes).
QF code (exact): none. no fermionic layer.

### Fermion-pair rotation of the vacuum
Topic: Fermionic canonical transformations
Appears in: X-004

The pair rotation turns the empty state into a superposition of empty and doubly occupied, with $\cos r$ and $\sin r$ amplitudes.

$$ e^{\xi \hat c_1^\dagger \hat c_2^\dagger-\bar\xi \hat c_2\hat c_1}|00\rangle=\cos r\,|00\rangle+e^{i\phi}\sin r\,|11\rangle $$

Conditions: $|11\rangle=\hat c_1^\dagger\hat c_2^\dagger|00\rangle$ (this ordering fixes the sign of the second term). It holds for all r; no chart condition is needed for the state.
Audit: Correct. Checked by: hand derivation: on span{|00>, |11>} the generator is the matrix ((0, -ξ̄), (ξ, 0)), whose exponential sends |00> to cos r |00> + (ξ/r) sin r |11>; and exact fermion matrices (2 modes, symbolic $r,\phi$, difference exactly zero).
QF code (exact): none. no fermionic layer.

### Hadamard series for a closed linear commutator action
Topic: BCH and Hadamard
Appears in: X-005 (source records cited: P-034, P-035, P-041, P-042, B-071 to B-074, C-018 to C-024)

When commuting with $\hat Q$ maps a list of operators linearly into itself, conjugation by $e^{t\hat Q}$ is the matrix exponential of that linear map.

$$ [\hat Q,\hat R_i]=\sum_jL_{ij}\hat R_j\ \Longrightarrow\ e^{t\hat Q}\hat R e^{-t\hat Q}=e^{tL}\hat R $$

Conditions: $L$ scalar; time-dependent noncommuting generators need time ordering.
Audit: Correct. Checked by: hand derivation: $\operatorname{ad}_{\hat Q}^k\hat R=L^k\hat R$ because $L$ is scalar, so the Hadamard series sums to $e^{tL}\hat R$; and QF exact for the general single-mode quadratic $\hat Q=p\hat a^2+q\hat a^{\dagger2}+s\hat a^\dagger\hat a$, where $[\hat Q,(\hat a,\hat a^\dagger)]$ gives $L=\begin{pmatrix}-s&-2q\\2p&s\end{pmatrix}$. The atlas remark that the adjoint action loses the overall phase is right.
QF code (exact):
```wl
{{c13}}
```
Output: `{\[FormalA]1*Cosh[Sqrt[(-4*p*qq + s^2)*t^2]] - (t*Sinh[Sqrt[(-4*p*qq + s^2)*t^2]]*(\[FormalA]1*s + 2*qq*SuperDagger[\[FormalA]1]))/Sqrt[(-4*p*qq + s^2)*t^2], 0, 0}`. Matches: the difference from the first row of $e^{tL}(\hat a,\hat a^\dagger)^T$ is zero symbolically and at a random complex point.

### Bosonic multimode Gaussian vacuum
Topic: Gaussian vacua and overlaps
Appears in: X-006 (source records cited: P-036 to P-045, B-048, B-061)

The multimode squeezed vacuum is a normalized pair-creation exponential on the vacuum, annihilated by the combinations $\hat a+Z\hat a^\dagger$.

$$ |Z\rangle_B=\det(I_N-ZZ^\dagger)^{1/4}e^{-\frac12\hat a^{\dagger T}Z\hat a^\dagger}|0\rangle,\qquad (\hat a+Z\hat a^\dagger)|Z\rangle=0 $$

Conditions: $Z^T=Z$, $I_N-ZZ^\dagger>0$.
Audit: Correct (normalization and annihilation condition). Consistent with X-009 ($Z=u^{-1}v$) and X-002. Checked by: QF exact for the annihilation condition: with $\hat X=-\frac12\hat a^{\dagger T}Z\hat a^\dagger$, QF gives $e^{\hat X}\hat a e^{-\hat X}=\hat a+Z\hat a^\dagger$, so $(\hat a+Z\hat a^\dagger)e^{\hat X}|0\rangle=e^{\hat X}\hat a|0\rangle=0$. Exact WL computation for the norm: in the Bargmann realization the squared norm is the real Gaussian integral of $\exp(-|z|^2-\operatorname{Re}z^TZz)$, whose real $2N\times2N$ form $A=\begin{pmatrix}I+X&-Y\\-Y&I-X\end{pmatrix}$ ($Z=X+iY$) has $\det A=\det(I_N-Z\bar Z)$ identically for $N=1,2,3$, so the squared norm before the prefactor is $\det(I_N-ZZ^\dagger)^{-1/2}$.
QF code (exact):
```wl
{{c14}}
```
Output: `{\[FormalA]1 + z11*SuperDagger[\[FormalA]1] + z12*SuperDagger[\[FormalA]2], \[FormalA]2 + z12*SuperDagger[\[FormalA]1] + z22*SuperDagger[\[FormalA]2]}`. Matches: these are the components of $\hat a+Z\hat a^\dagger$.

### Fermionic multimode Gaussian vacuum
Topic: Gaussian vacua and overlaps
Appears in: X-006

The fermionic analogue uses an antisymmetric $Z$ and a determinant with a plus sign, raised to the power $-1/4$.

$$ |Z\rangle_F=\det(I_N+ZZ^\dagger)^{-1/4}e^{-\frac12\hat c^{\dagger T}Z\hat c^\dagger}|0\rangle,\qquad (\hat c+Z\hat c^\dagger)|Z\rangle=0 $$

Conditions: $Z^T=-Z$. Exists on the chart where the canonical block $u$ is invertible (the atlas's "$U$" here is undefined in X-006; see Errors found).
Audit: Correct. Checked by: exact fermion matrices (4 modes, exact rational complex antisymmetric $Z$): the annihilation condition is exactly zero and the squared norm of the unnormalized state, squared, equals $\det(I_N+ZZ^\dagger)$ exactly.
QF code (exact): none. no fermionic layer.

### Number-conserving exponential in normal order
Topic: Normal symbols and normal-order exponentials
Appears in: X-007 (source records cited: B-013, B-028, B-040, B-107, G-009, G-097)

The exponential of a number-conserving quadratic equals the normal-ordered exponential of the same form with $h$ replaced by $e^h-I_N$.

$$ e^{\hat a^\dagger h \hat a}=:\!\exp[\hat a^\dagger(e^h-I_N)\hat a]\!: $$

Conditions: finite modes; the same identity for fermions with CAR normal ordering.
Audit: Correct for both statistics. Checked by: QF exact (single mode); exact WL computation for two bosonic modes in the Bargmann realization, where $\hat a^\dagger h\hat a\to\sum_{ij}z_ih_{ij}\partial_j$: the normal-ordered exponential acts on a generic cubic as $f(z)\mapsto f(z+B^Tz)$ exactly ($B=e^h-I$), and $f(e^{th^T}z)$ solves the flow equation of $e^{t\hat a^\dagger h\hat a}$ (residual zero at a random numeric point, symbolic $2\times2$ matrix exponential), so both sides send $f(z)$ to $f(e^{h^T}z)$; exact fermion matrices (3 modes, exact $h=P\,\mathrm{diag}(\log2,\log3,\log\frac15)P^{-1}$ with complex $P$, normal order $c^\dagger_{i_1}\cdots c^\dagger_{i_k}c_{j_k}\cdots c_{j_1}$): difference exactly zero.
QF code (exact):
```wl
{{c16}}
```
Output: `NormalOrdered[E^((-1 + E^lam)*SuperDagger[\[FormalA]1]**\[FormalA]1)]`. Matches (single mode).

### Vacuum projector as a normal-ordered exponential
Topic: Normal symbols and normal-order exponentials
Appears in: X-007

Normal ordering $e^{-\hat N}$ gives the projector onto the vacuum; for fermions the projector is the product of the single-mode empty projectors.

$$ |0\rangle\langle0|=:\!e^{-\sum_i\hat N_i}\!:\quad\text{(bosons)},\qquad |0\rangle\langle0|=\prod_i(\hat I-\hat N_i)\quad\text{(fermions)} $$

Audit: Correct; for fermions $:\!e^{-\sum_i\hat N_i}\!:$ also equals $\prod_i(\hat I-\hat N_i)$. Checked by: QF exact (single boson mode, all Fock indices symbolic, and coherent-state matrix elements; several modes follow because the normal-ordered exponential factorizes over modes) and exact fermion matrices (3 modes, both forms exactly $|0\rangle\langle0|$).
QF code (exact):
```wl
{{c17}}
```
Output: `{KroneckerDelta[0, n]*KroneckerDelta[m, n], E^((-Abs[w]^2 - Abs[z]^2)/2), Inactive[Sum][(-1)^k a†^k ** a^k / k!, {k, 0, Infinity}]}` (series index shortened). Matches: $\langle m|\cdot|n\rangle=\delta_{m0}\delta_{n0}$, and between normalized coherent states it is $\langle z|0\rangle\langle0|w\rangle$.

### Traces of second-quantized one-particle operators
Topic: Thermal states and traces
Appears in: X-007

The Fock-space trace of the second quantization of a one-particle matrix $M$ is the inverse determinant of $I-M$ for bosons and the determinant of $I+M$ for fermions.

$$ \operatorname{Tr}_{\mathcal F_B}\widehat\Gamma(M)=\det(I_N-M)^{-1},\qquad \operatorname{Tr}_{\mathcal F_F}\widehat\Gamma(M)=\det(I_N+M) $$

Conditions: $\widehat\Gamma(e^h)=e^{\hat a^\dagger h\hat a}$; bosons need convergence: spectral radius of $M$ below 1 (the atlas gives $0\le M<I_N$ for the positive case); infinite modes need trace class and $\|\hat M\|<1$ for bosons.
Audit: Correct. The atlas's "equivalently $\|M\|<1$" is only an equivalence for positive $M$ (see Errors found). Checked by: QF exact (single mode, with the convergence condition returned); exact WL computation for two bosonic modes: with $\langle z|\widehat\Gamma(M)|w\rangle=e^{\bar z^TMw}$ between unnormalized coherent states (the previous entry), the trace is the Gaussian integral of $e^{-\bar z^T(I-M)z}$ over $\mathbb C^2$ divided by $\pi^2$, evaluated exactly for a rational non-Hermitian $M$ with spectral radius below 1: it equals $1/\det(I_N-M)$; exact fermion matrices (3 modes, same exact $h$ as above): $\operatorname{Tr}e^{\hat c^\dagger h\hat c}-\det(I_N+e^h)$ is exactly zero.
QF code (exact):
```wl
{{c18}}
```
Output: `ConditionalExpression[-(-1 + E^lam)^(-1), E^Re[lam] < 1]`. Matches: $1/(1-e^\lambda)$ with the convergence condition $|e^\lambda|<1$.

### Multimode canonical action and its constraints
Topic: Bosonic canonical transformations (fermionic case: Fermionic canonical transformations)
Appears in: X-009 (source records cited: B-049, B-060, B-061, B-062, P-036, P-042, P-043)

A Gaussian unitary maps the annihilators to linear combinations of annihilators and creators, and the two matrix constraints say that the new operators obey the same commutation (or anticommutation) relations.

$$ \hat U\hat{\mathbf d}\hat U^{-1}=u\hat{\mathbf d}+v\hat{\mathbf d}^{\dagger},\qquad uu^\dagger-\eta vv^\dagger=I_N,\quad uv^T=\eta vu^T $$

Conditions: $\eta=+1$ bosons, $-1$ fermions; $\hat{\mathbf d}^\dagger$ a column of creators.
Audit: Correct (the two constraints are the CCR/CAR preserved by the map). Checked by: QF exact (bosons: with symbolic $u,v$, $[\hat b_i,\hat b_j^\dagger]=(uu^\dagger-vv^\dagger)_{ij}$ and $[\hat b_i,\hat b_j]=(uv^T-vu^T)_{ij}$ for $\hat b=u\hat a+v\hat a^\dagger$; and QF's exact Bogoliubov matrices of a network of two single-mode squeezers and a beam splitter satisfy both constraints exactly); exact fermion matrices (2 modes, symbolic $u,v$: $\{\hat b_i,\hat b_j^\dagger\}=(uu^\dagger+vv^\dagger)_{ij}$, $\{\hat b_i,\hat b_j\}=(uv^T+vu^T)_{ij}$; 3 and 4 modes, $u,v$ read off from $e^{-i\hat H}$ of a random complex even quadratic $\hat H$: both constraints hold to machine precision).
QF code (exact):
```wl
{{c19}}
```
Output: `{{{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}}`. Matches: the first two are the commutators minus the constraint matrices (symbolic), the last two are the constraints on QF's exact $u=\alpha^\dagger$, $v=\beta^\dagger$ (QF stores the creator action $\hat U\hat a_k^\dagger\hat U^\dagger=\sum_j\alpha_{jk}\hat a_j^\dagger+\beta_{jk}\hat a_j$).

### Multimode creator-number-annihilator factorization
Topic: Bosonic canonical transformations (fermionic case: Fermionic canonical transformations)
Appears in: X-009

Any Gaussian unitary with invertible $u$ factors into a pair-creation exponential, a number-conserving exponential, and a pair-annihilation exponential, times its vacuum amplitude.

$$ \hat U=c_0\exp\!\left[-\tfrac12\hat{\mathbf d}^{\dagger T}Z\hat{\mathbf d}^{\dagger}\right]\exp\!\left[\hat{\mathbf d}^{\dagger T}h\hat{\mathbf d}\right]\exp\!\left[\tfrac\eta2\hat{\mathbf d}^{T}W\hat{\mathbf d}\right],\quad Z=u^{-1}v,\ W=\bar v u^{-1},\ e^h=u^{-1},\ c_0=\langle0|\hat U|0\rangle,\ |c_0|=\det(uu^\dagger)^{-\eta/4} $$

Conditions: $u$ invertible (always so for bosons, since $uu^\dagger=I+vv^\dagger$; fermionic $u$ can be singular, e.g. the pair rotation of X-004 at $r=\pi/2$). The phase of $c_0$ is not fixed by $u,v$.
Audit: Correct, including the $\eta/2$ sign of the pair-annihilation factor ($Z$ and $W$ come out symmetric for bosons, antisymmetric for fermions). Checked by: hand derivation for bosons: QF gives $e^{\hat X}\hat a e^{-\hat X}=\hat a+Z\hat a^\dagger$ and $e^{\hat Y}\hat a^\dagger e^{-\hat Y}=\hat a^\dagger+W\hat a$ for the outer exponents $\hat X,\hat Y$, and the middle factor sends $\hat a\to u\hat a$, $\hat a^\dagger\to(u^{-1})^T\hat a^\dagger$; so the product sends $\hat a\to u\hat a+v\hat a^\dagger$ and $\hat a^\dagger\to W u\,\hat a+((u^{-1})^T+WuZ)\hat a^\dagger$, which is $\bar v\hat a+\bar u\hat a^\dagger$ because $u^T\bar u-v^\dagger v=I$ and $u^T\bar v=v^\dagger u$; canonical action plus the vacuum amplitude fix $\hat U$, and $|c_0|$ follows from the X-006 normalization since $\det(I-ZZ^\dagger)=\det(uu^\dagger)^{-1}$. The needed matrix identity, symmetry of $Z,W$, and $|c_0|=\det(uu^\dagger)^{-1/4}$ with QF's $c$ were confirmed on QF's exact Bogoliubov matrices of a four-element two-mode network, evaluated at 40 digits. Exact fermion matrices (3 and 4 modes): the factorization equals $e^{-i\hat H}$ to machine precision and $|c_0|=\det(uu^\dagger)^{1/4}$; with $+\tfrac12$ in place of $\eta/2$ the mismatch is O(1), so the sign is load-bearing.
QF code (exact):
```wl
{{c20}}
```
Output: `{{\[FormalA]1*Cosh[r] + E^(I*th)*Sinh[r]*SuperDagger[\[FormalA]2], \[FormalA]2*Cosh[r] + E^(I*th)*Sinh[r]*SuperDagger[\[FormalA]1]}, {{0, E^(I*th)*Tanh[r]}, {E^(I*th)*Tanh[r], 0}}, {{0, Tanh[r]/E^(I*th)}, {Tanh[r]/E^(I*th), 0}}, {{Sech[r], 0}, {0, Sech[r]}}, Sech[r], E^(-(E^(I*th)*SuperDagger[\[FormalA]1]**SuperDagger[\[FormalA]2]*Tanh[r]))**Sech[r]^(SuperDagger[\[FormalA]1]**\[FormalA]1 + SuperDagger[\[FormalA]2]**\[FormalA]2)**E^((\[FormalA]1**\[FormalA]2*Tanh[r])/E^(I*th))*Sech[r]}`. Matches for the two-mode squeezer $\hat U=\hat S_2(re^{i\theta})$: QF's canonical action gives $Z$, $W$, $u^{-1}=\operatorname{sech}r\,I$ and $\det(uu^\dagger)^{-1/4}=\operatorname{sech}r$, and the atlas product built from them, $c_0e^{-e^{i\theta}\tanh r\,\hat a_1^\dagger\hat a_2^\dagger}\operatorname{sech}r^{\,\hat N_1+\hat N_2}e^{e^{-i\theta}\tanh r\,\hat a_1\hat a_2}$, is QF's own ordering of $\hat U$.

### Bosonic coherent-state kernel of the factorized implementer
Topic: Bosonic canonical transformations
Appears in: X-009

Between unnormalized coherent states the factorized implementer has a Gaussian kernel whose quadratic form is read off from $Z$, $u^{-1}$ and $W$.

$$ K_{\hat U}(x,y)=c_0\exp\!\left[-\tfrac12x^TZx+x^Tu^{-1}y+\tfrac12y^TWy\right] $$

Conditions: kernel taken between unnormalized coherent states, $x=\bar z$, $y=w$ in $\langle z|\hat U|w\rangle$ with $|w\rangle=e^{w\cdot\hat a^\dagger}|0\rangle$. The atlas does not state this convention; it is the one that makes the formula hold.
Audit: Correct. Checked by: hand derivation: $\langle z|e^{-\frac12\hat a^{\dagger T}Z\hat a^\dagger}=e^{-\frac12\bar z^TZ\bar z}\langle z|$, $e^{\hat a^{\dagger T}h\hat a}|w\rangle=|e^hw\rangle$, $e^{\frac12\hat a^TW\hat a}|w\rangle=e^{\frac12w^TWw}|w\rangle$, and $\langle z|w'\rangle=e^{\bar z^Tw'}$; and QF exact: QF's package-private `gaussianKernel` on exact Bogoliubov matrices returns exactly the atlas blocks, and the single-mode middle line is QF's coherent matrix element of $e^{h\hat a^\dagger\hat a}$.
QF code (exact):
```wl
{{c21}}
```
Output: `{{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}, True, E^((-Abs[w]^2 - Abs[z]^2)/2 + E^h*w*Conjugate[z])}`. Matches: QF's kernel matrix minus $\begin{pmatrix}-Z&u^{-1}\\(u^{-1})^T&W\end{pmatrix}$ is exactly zero, its scalar is QF's $c$, and the normalized coherent element of $e^{h\hat a^\dagger\hat a}$ is $e^{e^h\bar zw}$ times the normalization factors.

---

## Non-formula records

- X-008: a five-step derivation route (commutator closure, ordered ansatz, Riccati ODEs, initial conditions and charts, limits for cubic and quartic generators) plus a statement of what was checked when the atlas was assembled. No formula; the claims are consistent with the checks above (the 2x2 compact and noncompact factorizations and the two-mode CAR realization all match).
