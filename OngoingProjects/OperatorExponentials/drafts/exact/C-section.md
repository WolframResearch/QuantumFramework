# C section: Combescure and Robert, Coherent States and Applications in Mathematical Physics (atlas lines 2981 to 4360)

## Errors found

No new error in the atlas text. The source misprints the atlas flags and corrects (C-018 missing 1/2, C-026 proof signs, C-028 amplitude, C-030 missing 2^n, C-047 sinh for tanh, C-050 endpoint, C-053 power, C-059 overlap parenthesis, C-067 sign, C-069 block-to-paired sign, C-076 argument sign, C-082 inverse label, C-085 Cayley sign, C-086 exponent sign and (-1)^n, C-091 Fresnel phase, C-094 global claim, C-098 Fourier of 1) were re-derived, and in each case the corrected form is the one that checks. Three corrections were tested against the uncorrected print: C-050 (the printed endpoint fails by an order-one residual in the 2-by-2 representation), C-085 (the one-mode hand check gives the minus sign) and C-098 (the exterior-algebra transform of 1 is the ordered top product, not alpha^* . alpha).

Convention that affects QF users of this section: QF's `SqueezeOperator[xi]` is exp((conj(xi) a^2 - xi a^dag^2)/2), the atlas operator S(-xi) of C-025, C-027, C-048.

Every QF cell below runs after these loading lines, each cell in a fresh session (scripts: `drafts/exact/C-work/cells/*.wl`, runner `runcells.wls`; non-QF exact checks: `exactchecks.wl`; outputs `cells_out.txt`, `exact_out.txt`):

```wl
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

### Canonical commutator of position and momentum
Topic: Infinite modes and domains
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

### Central-commutator BCH
Topic: BCH and Hadamard
Appears in: C-002 (also used in C-006, C-059, C-071)

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

### Weyl translation product and phase
Topic: Displacement
Appears in: C-003; same law as C-059 (one-mode re-expression in terms of Im<f,g>)

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

### Displacement conjugation and coordinate factorization
Topic: Displacement
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

### Normal and antinormal displacement factorization, displacement conjugation
Topic: Displacement
Appears in: C-006; same identity in C-059 (infinite-mode form)

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

### Coherent state: vacuum action, eigenvalue, number expansion
Topic: Coherent states (general)
Appears in: C-007; infinite-mode form Psi(f) in C-059

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

### Vacuum displacement amplitude and coherent overlap
Topic: Coherent states (general)
Appears in: C-008; same overlap in C-059 (as exp(-|f-g|^2/2 + i Im<f,g>))

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

### Harmonic rotation of displacements and coherent states
Topic: Coherent states (general)
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

### Coherent resolution of identity, trace and Hilbert-Schmidt norm
Topic: Coherent states (general)
Appears in: C-010

Coherent projectors integrated over phase space with weight (2 pi hbar)^{-n} give the identity, which turns traces and Hilbert-Schmidt norms into phase-space integrals of coherent expectation values.

$$ (2\pi\hbar)^{-n}\int|\phi_z\rangle\langle\phi_z|\,dz=\hat I,\quad \operatorname{Tr}\hat B=(2\pi\hbar)^{-n}\int\langle\phi_z,\hat B\phi_z\rangle dz,\quad \|\hat B\|_{\rm HS}^2=(2\pi\hbar)^{-n}\int\|\hat B\phi_z\|^2dz $$

Conditions: trace class, resp. Hilbert-Schmidt.
Audit: Correct. dq dp = 2 hbar d^2 alpha turns the prefactor into 1/pi per mode; the number-basis element <m|...|n> has an angular integral that vanishes for m != n and a radial part (2/n!) times the integral of r^{2n+1} e^{-r^2}, which exact integration gives as 1. Checked by: exact WL computation plus hand derivation.
QF code (exact): none. QF has no phase-space integration; `BosonicMatrixElement[{al, be}, 1, "Basis" -> "Coherent"]` gives only the overlaps in the integrand.

### Bargmann realization
Topic: Coherent states (general)
Appears in: C-011

In the Bargmann space of entire functions, creation acts as multiplication by zeta and annihilation as the derivative, so the linear and number exponentials act as a multiplication, a shift, and a dilation.

$$ \hat a^\dagger_j\mapsto\zeta_j,\quad \hat a_j\mapsto\partial_{\zeta_j},\quad \hat N\mapsto\zeta\cdot\partial_\zeta;\quad e^{v\cdot \hat a}f(\zeta)=f(\zeta+v),\ e^{u\cdot \hat a^\dagger}f=e^{u\cdot\zeta}f,\ e^{s\hat N}f(\zeta)=f(e^s\zeta) $$

Conditions: hbar = 1, analytic vectors.
Audit: Correct (Taylor shift, multiplication, Euler scaling). Checked by: hand derivation: e^{v d/dzeta} f(zeta) = sum_k v^k f^{(k)}(zeta)/k! = f(zeta + v), and zeta d/dzeta zeta^k = k zeta^k gives e^{sN} zeta^k = (e^s zeta)^k.
QF code (exact): none. QF has no function-space (Bargmann) realization; the operator statement behind the dilation, e^{sN} a e^{-sN} = e^{-s} a, is the C-094 call.

### Mehler kernel and Euclidean oscillator kernel
Topic: Quadratic propagators and Riccati
Appears in: C-012

The sum of products of Hermite functions weighted by w^k has a closed Gaussian form; with w = e^{-tau} and the zero-point factor it becomes the imaginary-time oscillator kernel.

$$ \sum_{k}w^{|k|}\phi_k(x)\phi_k(y)=\pi^{-n/2}(1-w^2)^{-n/2}\exp\!\left[-\frac{1+w^2}{2(1-w^2)}(|x|^2+|y|^2)+\frac{2w\,x\cdot y}{1-w^2}\right],\quad |w|<1 $$

and, n = 1, kernel of e^{-tau H_os} = (2 pi sinh tau)^{-1/2} exp[-(x^2+y^2)/(2 tanh tau) + xy/sinh tau].

Audit: Correct. Checked by: exact WL computation (n = 1, symbolic x, y: the Taylor coefficients of the closed form in w agree with the Hermite-function sum through order w^10; with w = e^{-tau} and the factor e^{-tau/2} the closed form divided by the tau kernel simplifies to 1); several dimensions factor into one-dimensional kernels.
QF code (exact): none. A position-space kernel; QF has no exact position-representation call for bosonic operators.

### Weyl quantization kernel convention
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-013

An operator acts through its Weyl symbol evaluated at the midpoint, and the symbol is recovered from the integral kernel by a Fourier transform in the difference variable.

$$ [\hat A\psi](x)=(2\pi\hbar)^{-n}\int A((x+y)/2,\xi)e^{i(x-y)\cdot\xi/\hbar}\psi(y)\,dy\,d\xi,\qquad A(x,\xi)=\int e^{-iu\cdot\xi/\hbar}K(x+u/2,x-u/2)\,du $$

Audit: Correct (the second is the Fourier inverse of the first in u = x - y). Checked by: hand derivation.
QF code (exact): none. QF has no Weyl-symbol map for operators.

### Symbol traces and Wigner overlaps
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-014

With the kernel convention of C-013, Hilbert-Schmidt norms, traces of products and expectation values become phase-space integrals of symbols with the factor (2 pi hbar)^{-n}.

$$ \|\hat A\|_{\rm HS}^2=(2\pi\hbar)^{-n}\int|A|^2,\quad\operatorname{Tr}(\hat A\hat B)=(2\pi\hbar)^{-n}\int AB,\quad \langle\psi,\hat A\phi\rangle=(2\pi\hbar)^{-n}\int A\,W_{\psi,\phi},\quad |\langle\psi,\phi\rangle|^2=(2\pi\hbar)^{-n}\int W_\psi W_\phi $$

Conditions: W is the Weyl symbol of the rank-one operator (unnormalized, integral (2 pi hbar)^n).
Audit: Correct for the C-013 kernel convention. Checked by: hand derivation.
QF code (exact): none. QF has no Weyl-symbol calculus; its `WignerFunction` takes a finite-size state and returns a normalized density.

### Symplectic covariance of Weyl symbols
Topic: Bosonic canonical transformations
Appears in: C-015

A metaplectic operator transforms Weyl symbols and Wigner functions by composing them with the classical linear canonical map.

$$ \hat R(F)^*\hat A\hat R(F)=\widehat{A\circ F},\qquad W_{\hat R(F)\psi,\hat R(F)\eta}(z)=W_{\psi,\eta}(F^{-1}z) $$

Audit: Correct (the second follows from the first through C-014). Checked by: hand derivation.
QF code (exact): none. QF has no Weyl-symbol calculus.

### Moyal product and exact quadratic commutator
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-016

The Weyl symbol of a product is the Moyal star product, and when one factor is at most quadratic the commutator is exactly i hbar times the quantized Poisson bracket.

$$ A\star B=A\exp\!\left[\frac{i\hbar}{2}\left(\overleftarrow\partial_q\cdot\overrightarrow\partial_p-\overleftarrow\partial_p\cdot\overrightarrow\partial_q\right)\right]B,\qquad [\hat H,\hat B]=i\hbar\,\widehat{\{H,B\}}\ (\deg H\le2) $$

Audit: Correct. q*p - p*q = i hbar reproduces [Q,P] = i hbar, and H = q, B = p gives i hbar {q,p} = i hbar. Checked by: hand derivation.
QF code (exact): none. QF has no Moyal product.

### Wick, anti-Wick and Weyl symbols related by Gaussian smoothing
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-017

The Weyl symbol is the anti-Wick symbol smoothed by a Gaussian of width set by hbar, and the Wick symbol is the Weyl symbol smoothed once more by the same Gaussian.

$$ A_w=G_\hbar*A^c=e^{\hbar\Delta/4}A^c,\qquad A_c=G_\hbar*A_w=e^{\hbar\Delta/4}A_w,\qquad G_\hbar(X)=(\pi\hbar)^{-n}e^{-|X|^2/\hbar} $$

Audit: Correct. The Weyl symbol of |phi_z><phi_z| is 2^n e^{-|X-z|^2/hbar}; with the (2 pi hbar)^{-n} measure this gives the kernel G; per coordinate variance hbar/2 gives e^{hbar Delta/4}. Checked by: hand derivation.
QF code (exact): none. QF's `WignerFunction`, `HusimiQFunction` and `SOrderedFunction` act on finite-size states, not on operator symbols.

### Linear quadratic flow and Heisenberg action
Topic: Quadratic propagators and Riccati
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

### Gaussian propagation by a linear-fractional map
Topic: Quadratic propagators and Riccati
Appears in: C-019

A Gaussian e^{ix^T Gamma x/(2 hbar)} stays Gaussian under a quadratic propagator: its covariance moves by the linear-fractional map of the flow blocks and its amplitude picks up det(A_t + B_t Gamma)^{-1/2}.

$$ \Gamma_t=(C_t+D_t\Gamma)(A_t+B_t\Gamma)^{-1},\qquad \hat U_t\big[a_\Gamma e^{ix^T\Gamma x/(2\hbar)}\big]=a_\Gamma\det(A_t+B_t\Gamma)^{-1/2}e^{ix^T\Gamma_tx/(2\hbar)} $$

Conditions: Gamma symmetric, Im Gamma > 0; root continued from t = 0.
Audit: Correct. Oscillator: Gamma = i stays i and det^{-1/2} = e^{-it/2} (zero-point phase); free particle: Gamma_t = i/(1+it), the known spreading packet. Checked by: hand derivation.
QF code (exact): none. A position-space Gaussian statement; QF has no exact position-representation call for bosonic states.

### Siegel-domain invariance and determinant cocycle
Topic: Quadratic propagators and Riccati
Appears in: C-020

The linear-fractional action of a symplectic matrix keeps a complex symmetric matrix with positive imaginary part inside that domain, composes as a group action, and its determinant factor is a cocycle.

$$ \Im\Gamma_F=(A+B\Gamma)^{-*}(\Im\Gamma)(A+B\Gamma)^{-1},\quad \alpha(F_1F_2)=\alpha(F_1)\alpha(F_2),\quad \delta(F_1F_2,\Gamma)=\delta(F_1,\alpha(F_2)\Gamma)\delta(F_2,\Gamma) $$

Audit: Correct. Checked by: WL computation with finite 4-by-4 symplectic matrices (two degrees of freedom, F1 and F2 as matrix exponentials, complex Gamma: all three identities hold to machine precision; Gamma_F is symmetric with positive imaginary part).
QF code (exact): none. Classical matrix algebra; QF's second-quantization layer does not apply.

### Matrix Riccati and amplitude equations
Topic: Quadratic propagators and Riccati
Appears in: C-021

The complex covariance of a propagated Gaussian obeys a matrix Riccati equation and its amplitude a linear equation with a trace.

$$ \dot\Gamma=-K-L\Gamma-\Gamma L^T-\Gamma G\Gamma,\qquad \dot a=-\tfrac12\operatorname{Tr}(L+G\Gamma)a,\qquad S=\begin{pmatrix}K&L\\L^T&G\end{pmatrix} $$

Audit: Correct. Forcing line: with M = A + B Gamma0, N = C + D Gamma0, d/dt(M;N) = JS(M;N) = (L^T M + G N; -K M - L N), and Gamma = N M^{-1} gives the Riccati equation; a ~ det(M)^{-1/2} gives the trace. Checked by: hand derivation, plus a WL computation with finite matrices (central-difference derivative of the linear-fractional map agrees to finite-difference accuracy).
QF code (exact): none. Classical matrix algebra; QF's second-quantization layer does not apply.

### Coherent matrix elements of a metaplectic operator
Topic: Quadratic propagators and Riccati
Appears in: C-022

The coherent-state matrix element of the propagator of a quadratic Hamiltonian is a closed Gaussian in the two labels, with det(Y)^{-1/2}, Y = A + D + i(B - C), as prefactor.

$$ \langle\phi_z,\hat U\phi_{z'}\rangle=2^{n/2}\det(Y)^{-1/2}e^{-i\sigma(Fz',z)/2}e^{-(|x|^2+|\xi|^2)/4}e^{-(\xi+ix)^TW(\xi+ix)/4},\quad Y=A+D+i(B-C),\ W=ZY^{-1},\ (x,\xi)=z-Fz' $$

Conditions: hbar = 1, Z = A - D + i(B+C), root along the path.
Audit: Correct. Checked by: exact WL computation (one mode, S = {{9/10, 3/10}, {3/10, 3/5}}, t = 4/5, hbar = 1: H written in su(1,1) generators, checked by QF normal ordering; U = e^{-itH} disentangled as e^{c K+} e^{d K0} e^{f K-} from the exact 2-by-2 representation; its coherent elements written in closed form and compared with the atlas formula at four label pairs at 40 digits: equal, while dropping the sigma phase leaves an order-one residual; the 2-by-2 entry whose square root is taken stays away from the negative axis along the path, so the principal root is the continuous one).
QF code (exact): none. `BosonicExpOrder[Exp[x1 a**a + x2 ad**ad + x3 ad**a]]` returns the exponential unchanged, and `BosonicMatrixElement[..., "Basis" -> "Coherent"]` leaves exponentials unevaluated.

### Metaplectic generator from a symplectic path
Topic: Bosonic canonical transformations
Appears in: C-023

The quadratic generator of a symplectic path is read off as -J F' F^{-1}, and implementers compose only up to a sign (the metaplectic double cover).

$$ S_t=-J\dot F_tF_t^{-1},\qquad \hat R(F_1)\hat R(F_2)=\pm\hat R(F_1F_2) $$

Audit: Correct (from F' = J S F and J^{-1} = -J; the sign is the metaplectic double cover). Checked by: hand derivation.
QF code (exact): none. Classical matrix statement; QF has no metaplectic-lift routine.

### Bogoliubov blocks of a symplectic matrix
Topic: Bosonic canonical transformations
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

### Multimode squeeze operator and bounded squeeze coordinate
Topic: Single-mode squeezing
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

### Squeeze conjugation in block form
Topic: Single-mode squeezing
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

### Squeezed vacuum from the Bargmann Riccati equation
Topic: Gaussian vacua and overlaps
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

### Coordinate form of the squeezed vacuum
Topic: Gaussian vacua and overlaps
Appears in: C-028; one-mode case in C-049 (delta = (1-omega)/(1+omega))

In position space the squeezed vacuum is a normalized complex Gaussian whose covariance is the Cayley transform of W.

$$ \Gamma=i(I_n-W)(I_n+W)^{-1},\qquad \psi_B(x)=(\pi\hbar)^{-n/4}\det(I_n-W^*W)^{1/4}\det(I_n+W)^{-1/2}e^{ix^T\Gamma x/(2\hbar)} $$

Audit: Correct (the atlas's corrected amplitude). One mode, W = tanh r: the generator is -r(x d/dx + 1/2), so S|0> = pi^{-1/4} e^{-r/2} e^{-e^{-2r}x^2/2}, and (1 - tanh^2)^{1/4}(1 + tanh)^{-1/2} = e^{-r/2}; complex W: the norm integral gives 1. Checked by: hand derivation.
QF code (exact): none. A position-space wavefunction; QF has no exact position-representation call for bosonic states.

### Squeeze-rotation factorization and its vacuum phase
Topic: Quadratic propagators and Riccati
Appears in: C-029

A quadratic propagator factors into a squeeze and a number-preserving rotation; the rotation multiplies the vacuum by a unit-modulus phase fixed by det(Y_t) and det(u_t u_t^dag).

$$ \hat U_t=\hat S(B_t)\hat R_t,\quad \hat R_t|0\rangle=\vartheta_t|0\rangle,\quad \vartheta_t=\frac{2^{n/2}\det(Y_t)^{-1/2}}{\det(u_tu_t^\dagger)^{-1/4}},\quad |\vartheta_t|=1,\quad Y_t=2\bar u_t $$

Audit: Correct. |det Y|^{-1/2} 2^{n/2} = |det u|^{-1/2} = det(u u^dag)^{-1/4}, so |vartheta| = 1; Y = 2 conj(u) follows from C-024. Checked by: hand derivation (consistency of C-022, C-024, C-027).
QF code (exact): none. QF's exact layer has no metaplectic factorization; `QuantumInterferometer` works on cut-off Fock spaces and is not used.

### Weyl symbol of a metaplectic operator (Cayley form)
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-030; oscillator case C-053

Away from det(I + F) = 0, the Weyl symbol of a metaplectic operator is a Gaussian in the Cayley transform of F, with prefactor 2^n |det(I + F)|^{-1/2} and a path-dependent phase.

$$ R_w(F;X)=\frac{2^n e^{i\pi\nu}}{|\det(I_{2n}+F)|^{1/2}}\exp\!\left[-\frac{i}{\hbar}X^TJ(I_{2n}-F)(I_{2n}+F)^{-1}X\right] $$

Conditions: det(I+F) != 0; nu and phase from the path.
Audit: Correct (with the atlas's added 2^n). F = I gives 1. Oscillator F = e^{tJ}: (I-F)(I+F)^{-1} = -tan(t/2) J, so the exponent is -i tan(t/2)|X|^2/hbar and the prefactor 1/cos(t/2), the known symbol of e^{-itH_os}. Checked by: hand derivation.
QF code (exact): none. QF has no Weyl-symbol calculus.

### Translation coefficient of a metaplectic operator
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-031

The translation coefficient Tr[T(-z) R(F)] of a metaplectic operator is the Fresnel transform of its Weyl symbol: a Gaussian in the inverse Cayley transform with a signature phase.

$$ \widetilde R(F;z)=\frac{c_w\,e^{-i\pi\operatorname{sgn}(Q_F)/4}}{2^n|\det Q_F|^{1/2}}\exp\!\left[-\frac{i}{4\hbar}z^TJ(I+F)(I-F)^{-1}z\right],\quad Q_F=J(I-F)(I+F)^{-1};\qquad \widetilde R=\frac{1}{2i\sin(t/2)}e^{\frac{i}{4\hbar}\cot(t/2)|z|^2}\ \text{for}\ \hat R=e^{-it(\hat N+1/2)} $$

Conditions: I + F and I - F invertible; 0 < t < pi for the example.
Audit: Correct. Direct trace: Tr[D(-alpha) e^{-it(N+1/2)}] = e^{-it/2} e^{-x/2} sum L_n(x) w^n = 1/(2i sin(t/2)) exp(i x cot(t/2)/2), w = e^{-it}, x = |z|^2/(2 hbar), via the Laguerre generating function; the general chart formula with c_w = 1/cos(t/2), Q_F = tan(t/2) I, sgn = 2 gives the same. Checked by: hand derivation.
QF code (exact): none. QF has no trace of an operator on the full Fock space; its exact Laguerre matrix elements of `DisplacementOperator[al]` are one input of the hand derivation, not a trace.

### Propagator of a driven quadratic Hamiltonian
Topic: Forced and parametric oscillators
Appears in: C-032

For a Hamiltonian that is quadratic plus linear, the propagator is the classical action phase times a translation along the classical orbit, the quadratic propagator, and the inverse translation at the start time.

$$ \delta_t=\int_0^t[p_s\dot q_s-H(z_s,s)]ds-\tfrac12(q_tp_t-qp),\qquad \hat U_2(t,s)=e^{i(\delta_t-\delta_s)/\hbar}\hat T(z_t)\hat S(B_t)\hat R_t\hat R_s^{-1}\hat S(-B_s)\hat T(-z_s) $$

Audit: Correct for an affine quadratic Hamiltonian. Checked by: exact WL computation (H = N + 1/2 + f Q with symbolic f, q, p, t, hbar = 1, B = 0, R_t = e^{-it(N+1/2)}, s = 0: in the Bargmann realization a -> d/dz, a^dag -> z, the right-hand side applied to an arbitrary analytic function solves i d_t psi = H psi identically, FullSimplify giving 0, with the classical orbit and delta_t computed symbolically and the product equal to the identity at t = 0); general s follows from U(t,s) = U(t,0) U(s,0)^{-1}.
QF code (exact): none. `BosonicExpOrder[Exp[-I t (ad**a + 1/2 + f (a + ad)/Sqrt[2])]]` returns the exponential unfactored.

### Periodic quadratic evolution and quasi-energies
Topic: Forced and parametric oscillators
Appears in: C-035

For a periodic quadratic Hamiltonian in a stable region, the one-period propagator has eigenvalues of the form e^{-i rho pi (k+1/2)}.

$$ \text{Floquet eigenvalues } e^{-i\rho\pi(k+1/2)} $$

Audit: Unchecked: needs Mathieu-equation stability data; only the constant-frequency limit (rho = omega, period pi) was checked by hand.
QF code (exact): none. QF has no Floquet or classical-stability routine.

### Periodic-orbit determinant in regularized traces
Topic: Infinite modes and domains
Appears in: C-036

Energy-localized traces of the propagator are asymptotically sums over periodic orbits, each weighted by |det(I - P_gamma)|^{-1/2} and a Maslov phase.

$$ \operatorname{Tr}[\chi(\hat H)\hat A\chi(\hat H)g((E\hat I-\hat H)/\hbar)]\ \sim\ \ldots\,|\det(I_{d_\gamma}-P_\gamma)|^{-1/2}e^{i\pi(\cdot)} $$

Audit: Unchecked: semiclassical stationary-phase asymptotics; no finite identity to test.
QF code (exact): none. QF has no semiclassical trace routine.

### Clock and shift matrices on the torus
Topic: Lie-group decompositions and generalized coherent states
Appears in: C-037

On the quantum torus with 2 pi hbar N = 1, translations become N-by-N clock and shift matrices of order N that commute up to a root of unity.

$$ 2\pi\hbar N=1,\quad \hat X^N=\hat Z^N=\hat I,\quad \hat Z\hat X=\omega\hat X\hat Z,\ \omega=e^{2\pi i/N} $$

Audit: Correct for Z|j> = omega^j |j>, X|j> = |j+1>. Checked by: exact WL computation (N = 5; the torus space is finite, so these are the complete operators). The book's (6.18)-(6.19) matrices were not seen; prior audit only for the claim that they fix this orientation.
QF code (exact): none. Finite clock and shift matrices; QF's second-quantization layer does not apply.

### SU(2) commutators, ladder matrices, Casimir
Topic: SU(2) and passive mode mixing
Appears in: C-038

The spin-j ladder operators satisfy the su(2) relations, raise and lower m with the standard square-root factors, and give the Casimir j(j+1).

$$ [\hat J_3,\hat J_\pm]=\pm\hat J_\pm,\ [\hat J_+,\hat J_-]=2\hat J_3,\ \hat J_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}|j,m\pm1\rangle,\ \hat J^2=j(j+1)\hat I $$

Audit: Correct. Checked by: exact WL computation (j = 3/2, complete 4-by-4 spin matrices).
QF code (exact): none. Finite spin matrices; QF's second-quantization layer does not apply.

### Exponential of an su(2) matrix
Topic: SU(2) and passive mode mixing
Appears in: C-039

A 2-by-2 su(2) matrix squares to minus its determinant, so its exponential is a cosine plus a sine term.

$$ e^{tA}=I_2\cos t+A\sin t\ (\det A=1),\qquad e^{it\mathbf v\cdot\boldsymbol\sigma}=\cos(t|\mathbf v|)I_2+i\sin(t|\mathbf v|)\mathbf n\cdot\boldsymbol\sigma $$

Audit: Correct (A^2 = -I). Checked by: exact WL computation (2-by-2 matrices).
QF code (exact): none. A 2-by-2 matrix identity; QF's second-quantization layer does not apply.

### Gauss decomposition of SL(2,C)
Topic: Lie-group decompositions and generalized coherent states
Appears in: C-040

A unimodular 2-by-2 matrix with nonzero lower-right entry factors uniquely into upper-triangular, diagonal and lower-triangular pieces.

$$ g=\begin{pmatrix}1&\beta/\delta\\0&1\end{pmatrix}\begin{pmatrix}\delta^{-1}&0\\0&\delta\end{pmatrix}\begin{pmatrix}1&0\\\gamma/\delta&1\end{pmatrix} $$

Conditions: alpha delta - beta gamma = 1, delta != 0.
Audit: Correct. Checked by: exact WL computation (symbolic product returns g).
QF code (exact): none. A 2-by-2 matrix identity; QF's second-quantization layer does not apply.

### SU(2) normal factorization and spin coherent state
Topic: SU(2) and passive mode mixing
Appears in: C-041; fermion-pair application in C-087

The SU(2) displacement factors into a raising exponential, a J_3 exponential and a lowering exponential with tan r, and acting on the lowest weight it gives the normalized spin coherent state.

$$ e^{\xi\hat J_+-\bar\xi\hat J_-}=e^{\zeta\hat J_+}e^{-2\log(\cos r)\hat J_3}e^{-\bar\zeta\hat J_-},\ \zeta=\frac{\xi}{r}\tan r;\qquad \hat G(\xi)|j,-j\rangle=(1+|\zeta|^2)^{-j}e^{\zeta\hat J_+}|j,-j\rangle $$

Conditions: cos r != 0.
Audit: Correct. Checked by: exact WL computation (j = 3/2 complete spin matrices, xi = 0.7 e^{0.4 i}: both identities hold to machine precision).
QF code (exact): none. Finite spin matrices; QF's second-quantization layer does not apply.

### Spin coherent overlap and resolution of identity
Topic: Lie-group decompositions and generalized coherent states
Appears in: C-042

Spin coherent states overlap as (1 + conj(eta) zeta)^{2j} over the norms, with modulus set by the angle between the two directions, and resolve the identity with the measure (2j+1)/pi (1+|zeta|^2)^{-2}.

$$ \langle\eta|\zeta\rangle=\frac{(1+\bar\eta\zeta)^{2j}}{(1+|\eta|^2)^j(1+|\zeta|^2)^j},\quad |\langle n'|n\rangle|^2=\left(\frac{1+n\cdot n'}2\right)^{2j},\quad \frac{2j+1}{\pi}\int\frac{|\zeta\rangle\langle\zeta|}{(1+|\zeta|^2)^2}d^2\zeta=\hat I $$

Audit: Correct. Checked by: exact WL computation for both overlaps (j = 3/2 complete spin matrices); resolution by hand (each diagonal element is (2j+1) C(2j,k) B(k+1, 2j+1-k) = 1).
QF code (exact): none. Finite spin matrices; QF's second-quantization layer does not apply.

### Berezin-Lieb bounds for spin
Topic: Thermal states and traces
Appears in: C-043

The trace of a convex function of a spin Hamiltonian is bounded below by the integral over its covariant symbol and above by the integral over its contravariant symbol; for J_l these symbols are j n_l and (j+1) n_l.

$$ \int\chi(H_c)d\mu\le\operatorname{Tr}\chi(\hat H)\le\int\chi(H^c)d\mu,\qquad (\hat J_\ell)_c=jn_\ell,\ (\hat J_\ell)^c=(j+1)n_\ell $$

Audit: Correct (standard Berezin-Lieb inequality and spin symbols). Checked by: hand derivation for the covariant symbol (<j,-j| rotated J_l |j,-j> = j n_l); contravariant value and the inequality from the standard result, not recomputed.
QF code (exact): none. Spin symbols; QF's second-quantization layer does not apply.

### Spin-to-oscillator contraction
Topic: Lie-group decompositions and generalized coherent states
Appears in: C-044

As j grows, the spin coherent kernel at scale zeta/sqrt(2j) tends to the bosonic coherent kernel.

$$ \text{spin kernel at scale }\zeta/\sqrt{2j}\ \to\ e^{\bar\eta\zeta-|\eta|^2/2}\quad (j\to\infty) $$

Audit: Correct as a limit ((1 + conj(eta) zeta/(2j))^{2j} -> e^{conj(eta) zeta}). Checked by: hand derivation.
QF code (exact): none. A large-spin limit; QF's second-quantization layer does not apply.

### SU(1,1) algebra and Bergman realization
Topic: Lie-group decompositions and generalized coherent states
Appears in: C-045

The su(1,1) generators close with [K-, K+] = 2 K0 and Casimir k(k-1), and act on holomorphic functions on the disk as first-order differential operators.

$$ [\hat K_-,\hat K_+]=2\hat K_0,\ [\hat K_0,\hat K_\pm]=\pm\hat K_\pm,\ \hat C=\hat K_0^2-\tfrac12(\hat K_+\hat K_-+\hat K_-\hat K_+)=k(k-1)\hat I;\quad \hat K_-=\partial_z,\ \hat K_+=2kz+z^2\partial_z,\ \hat K_0=k+z\partial_z $$

Audit: Correct. Checked by: exact WL computation (symbolic action on an arbitrary f(z): all four relations vanish).
QF code (exact): none. An abstract differential realization; the bosonic realizations are the QF call of C-046.

### Bosonic realizations of SU(1,1)
Topic: Two-mode squeezing
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

### SU(1,1) coherent series and overlap
Topic: Lie-group decompositions and generalized coherent states
Appears in: C-047; overlap modulus reused in C-097

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

### Noncompact Gauss factorization of squeezers
Topic: Single-mode squeezing
Appears in: C-048 (two-mode version: Two-mode squeezing)

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

### Squeezed number operator
Topic: Single-mode squeezing
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

### Wigner function of a one-mode squeezed Gaussian
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-049

The Wigner function of the one-mode squeezed Gaussian e^{-delta x^2/(2 hbar)} is a tilted Gaussian ellipse fixed by the real and imaginary parts of delta.

$$ W(q,p)=2\exp[-q^2\Re\delta/\hbar-(p+q\Im\delta)^2/(\hbar\Re\delta)],\qquad \psi\propto e^{-\delta x^2/(2\hbar)},\ \delta=\frac{1-\omega}{1+\omega} $$

Audit: Correct (unnormalized Weyl symbol of the projector; integral 2 pi hbar). Forcing line: psi(q+u/2) conj(psi(q-u/2)) = exp(-Re(delta) q^2 - Re(delta) u^2/4 - i Im(delta) q u), whose u-Fourier transform gives the (p + q Im delta)^2 term. Checked by: hand derivation.
QF code (exact): none. QF's `WignerFunction` takes a finite-size state; QF has no exact Weyl symbol of an operator.

### Time-dependent SU(1,1) Hamiltonian
Topic: Forced and parametric oscillators
Appears in: C-050

For a Hamiltonian linear in the su(1,1) generators with time-dependent coefficients, the propagator is a squeeze, a K0 phase rotation and an inverse squeeze, with the squeeze coordinate solving a scalar Riccati equation.

$$ i\dot\omega=\bar\lambda\omega^2+\mu\omega+\lambda,\ \dot\gamma=-\lambda\bar\omega-\bar\lambda\omega-\mu,\qquad \hat U(t,s)=\hat S(\beta_t)e^{i(\gamma_t-\gamma_s)\hat K_0}\hat S(-\beta_s),\ \omega_t=\frac{\beta_t}{|\beta_t|}\tanh|\beta_t| $$

Audit: Correct (the atlas's corrected endpoint). Checked by: WL computation in the exact 2-by-2 representation of su(1,1) (constant lambda, mu; Riccati integrated numerically: U(t,0) = S(beta_t) e^{i gamma_t K0} agrees at three times to the integration accuracy; the corrected two-time form agrees; the printed form with S(-beta_t) fails by an order-one residual).
QF code (exact): none. QF has no time-dependent disentangling routine.

### Duhamel derivative of an exponential
Topic: BCH and Hadamard
Appears in: C-051

The derivative of the exponential of a moving generator is an average of the derivative conjugated along the exponential, not simply L' e^L.

$$ \frac d{dt}e^{\hat L(t)}=\int_0^1e^{s\hat L(t)}\dot{\hat L}(t)e^{(1-s)\hat L(t)}ds $$

Audit: Correct. Checked by: WL computation with finite 3-by-3 matrices (central difference against quadrature of the Duhamel integral).
QF code (exact): none. A general matrix identity; QF's second-quantization layer does not apply.

### Singular oscillator lowest weight
Topic: Lie-group decompositions and generalized coherent states
Appears in: C-052

For the oscillator with an inverse-square potential on the half line, the lowest-weight state is x^s e^{-x^2/2} with s fixed by g, and its Bargmann index is k = (2s+1)/4.

$$ \psi\propto x^se^{-x^2/2},\quad s=\tfrac12+\sqrt{\tfrac14+g^2},\quad k=\frac{2s+1}{4} $$

Conditions: H = (-d^2 + x^2 + g^2/x^2)/2 on L^2(R+), chosen self-adjoint realization.
Audit: Correct (energy (2s+1)/2 = 2k with K0 = H/2). Checked by: exact WL computation (symbolic eigenvalue check).
QF code (exact): none. A differential operator on L^2(R+), outside QF's Fock-space layer.

### Weyl symbol of the oscillator exponential
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-053 (special case of C-030)

The Weyl symbol of the oscillator propagator is (cos(t/2))^{-n} e^{i tan(t/2)|X|^2}, and in imaginary time (cosh(tau/2))^{-n} e^{-tanh(tau/2)|X|^2}.

$$ [e^{it\hat H_{\rm os}}]_w=(\cos(t/2))^{-n}\exp[i\tan(t/2)|X|^2],\qquad [e^{-\tau\hat H_{\rm os}}]_w=(\cosh(\tau/2))^{-n}\exp[-\tanh(\tau/2)|X|^2] $$

Conditions: hbar = 1, |t| < pi.
Audit: Correct (power -n). Check: (2 pi)^{-1} times the integral of the tau symbol is 1/(2 sinh(tau/2)) = Tr e^{-tau H_os}. Checked by: hand derivation.
QF code (exact): none. QF has no Weyl-symbol calculus.

### Herman-Kluk representation, exact for quadratics
Topic: Quadratic propagators and Riccati
Appears in: C-054

For a quadratic Hamiltonian, the phase-space integral of |phi_{F_t z}><phi_z| with the prefactor det(u_t)^{1/2} reproduces the propagator exactly.

$$ \hat U_t=(2\pi\hbar)^{-n}\int C_t\,|\phi_{F_tz}\rangle\langle\phi_z|\,dz,\qquad C_t=\det\!\left[\frac{A_t+D_t+i(C_t^{\rm block}-B_t)}2\right]^{1/2} $$

Audit: Correct (the atlas's normalized re-expression; C_t = det u^{1/2}). Checked by: exact WL computation (the C-022 data: the coherent element of the Herman-Kluk integral is a two-dimensional Gaussian integral, done in closed form with positive-definite real part, and it equals the exact disentangled matrix element at four label pairs at 40 digits; coherent elements determine the operator).
QF code (exact): none. `BosonicMatrixElement[..., "Basis" -> "Coherent"]` leaves exponentials unevaluated and QF has no phase-space integration.

### Gaussian density operators, purity, thermal factors
Topic: Thermal states and traces
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

### Reduced Gaussian channel from a joint quadratic evolution
Topic: Gaussian vacua and overlaps
Appears in: C-056

A joint quadratic evolution with a Gaussian environment acts on the system by a linear map of the characteristic-function argument times a Gaussian noise factor.

$$ \chi_{S,\rm out}(k)=\chi_{S,\rm in}(X^Tk)\exp[-\tfrac12k^TYV_EY^Tk],\qquad V_{S,\rm out}=XV_SX^T+YV_EY^T $$

Audit: Correct. Forcing line: e^{ik.(Xz+Yu)} factors into system and environment parts, and the centered Gaussian environment has chi_E(m) = exp(-m^T V_E m/2). Checked by: hand derivation.
QF code (exact): none. QF has no characteristic-function calculus for bosonic channels.

### Segal fields
Topic: Infinite modes and domains
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

### Second quantization of a one-particle unitary
Topic: SU(2) and passive mode mixing
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

### Infinite-mode displacement, exponential vector, product, overlap
Topic: Infinite modes and domains
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

### Displacement conjugation of a field and of N
Topic: Displacement
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

### Derivative of a moving displacement
Topic: Displacement
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

### CAR and finite fermionic Fock space
Topic: Fermionic canonical transformations
Appears in: C-064

Fermionic modes satisfy canonical anticommutation relations, so n modes span a 2^n-dimensional Fock space and a smeared annihilator has operator norm ||f||.

$$ \{\hat c_j,\hat c_k^\dagger\}=\delta_{jk}\hat I,\ \{\hat c_j,\hat c_k\}=0,\ \|\hat c(f)\|=\|f\| $$

Audit: Correct. Checked by: exact fermion matrices (3 modes) for the CAR; hand derivation for the norm (c(f)^dag c(f) has eigenvalues 0 and ||f||^2 by the CAR).
QF code (exact): none. no fermionic layer.

### Linear fermionic exponential with complex coefficients
Topic: Fermionic canonical transformations
Appears in: C-065

A linear combination of fermionic creators and annihilators with complex coefficients squares to the scalar u.v, so its exponential is cosh plus a sinh-over-root multiple of the generator.

$$ e^{\hat L}=\cosh\sqrt{u\cdot v}\,\hat I+\frac{\sinh\sqrt{u\cdot v}}{\sqrt{u\cdot v}}(u\cdot\hat c^\dagger+v\cdot\hat c),\qquad \hat L^2=(u\cdot v)\hat I $$

Audit: Correct. Checked by: exact fermion matrices (2 modes, complex u, v).
QF code (exact): none. no fermionic layer.

### Grassmann derivatives: Leibniz rule and left/right relation
Topic: Grassmann calculus and fermionic coherent states
Appears in: C-066

The Grassmann derivative obeys a graded Leibniz rule, and the right and left derivatives differ by a parity sign.

$$ \partial_\theta(fg)=(\partial_\theta f)g+(-1)^{\pi(f)}f(\partial_\theta g),\qquad f\overleftarrow\partial_\theta=(-1)^{1+\pi(f)}\overrightarrow\partial_\theta f $$

Audit: Correct. Checked by: exact WL computation (exterior algebra, odd and even test elements).
QF code (exact): none. no fermionic layer.

### Berezin integration: determinant rule and paired normalization
Topic: Grassmann calculus and fermionic coherent states
Appears in: C-067; the same paired functional is used in C-068, C-069, C-098, C-099

Berezin integration rescales by det A under a linear change of generators, and the paired convention normalizes e^{-theta^* theta} to 1.

$$ \int f(A\theta)\,d\theta=\det A\int f(\theta)\,d\theta,\qquad \widehat I_\theta\Big[\prod_j\theta_j\theta_j^*\Big]=1,\ \widehat I_\theta[e^{-\theta^*\theta}]=1,\ \widehat I_\theta=\partial_{\theta_1^*}\partial_{\theta_1}\cdots\partial_{\theta_n^*}\partial_{\theta_n} $$

Audit: Correct as a declared convention (the atlas's repair). Checked by: exact WL computation (exterior algebra: 3-by-3 determinant rule; paired Gaussian gives 1 for n = 2).
QF code (exact): none. no fermionic layer.

### Fermionic Gaussian integral with sources
Topic: Gaussian integrals and Pfaffians
Appears in: C-068

The paired Berezin integral of e^{-theta^* B theta} is det B, and linear sources complete the square to det B e^{gamma^* B^{-1} gamma}.

$$ \widehat I_\theta e^{-\theta^*B\theta}=\det B,\qquad \widehat I_\theta e^{-\theta^*B\theta+\gamma^*\theta+\theta^*\gamma}=\det B\,e^{\gamma^*B^{-1}\gamma} $$

Audit: Correct in the C-067 convention. Checked by: exact WL computation (exterior algebra, n = 2, complex non-Hermitian B, every Grassmann coefficient matches).
QF code (exact): none. no fermionic layer.

### Pfaffian Gaussian integral
Topic: Gaussian integrals and Pfaffians
Appears in: C-069

A general Grassmann Gaussian integrates to a Pfaffian with the block-to-paired sign (-1)^{n(n-1)/2}, and sources put the inverse matrix in the exponent.

$$ \widehat I_\gamma e^{\Phi+\gamma^*\xi-\xi^*\gamma}=(-1)^{n(n-1)/2}\operatorname{Pf}(\Lambda)\exp\!\left[\tfrac12\binom{\xi^*}{\xi}^{T}\Lambda^{-1}\binom{\xi^*}{\xi}\right],\quad \Phi=\tfrac12x^T\Lambda x $$

Audit: Correct (with the atlas's block-to-paired sign). Checked by: exact WL computation (exterior algebra, n = 2, complex K, L, M: all coefficients match; the example Phi = -gamma^* gamma gives integral 1 and block Pfaffian -1; Pf^2 = det).
QF code (exact): none. no fermionic layer.

### Normal-ordered monomial expansion of fermionic operators
Topic: Grassmann calculus and fermionic coherent states
Appears in: C-070

Every operator on n fermionic modes is a unique combination of normal-ordered monomials.

$$ \hat A=\sum_{\epsilon,\epsilon'\in\{0,1\}^n}A_{\epsilon\epsilon'}(\hat c^\dagger)^\epsilon\hat c^{\epsilon'} $$

Audit: Correct (4^n monomials span the 2^n-by-2^n matrices). Checked by: hand derivation (dimension count).
QF code (exact): none. no fermionic layer.

### Grassmann Weyl translation
Topic: Grassmann calculus and fermionic coherent states
Appears in: C-071

The fermionic translation with Grassmann labels has an even generator, so the central BCH formula factorizes it and composes it up to a scalar Grassmann phase, and it shifts c by gamma.

$$ \hat T(\gamma)=e^{\hat c^\dagger\cdot\gamma-\gamma^*\cdot\hat c}=e^{-\gamma^*\cdot\gamma/2}e^{\hat c^\dagger\cdot\gamma}e^{-\gamma^*\cdot\hat c},\quad \hat T(\alpha)\hat T(\gamma)=\hat T(\alpha+\gamma)e^{-(\alpha^*\cdot\gamma+\alpha\cdot\gamma^*)/2},\quad \hat T^{-1}\hat c\hat T=\hat c+\gamma,\ \hat T^{-1}\hat c^\dagger\hat T=\hat c^\dagger+\gamma^* $$

Audit: Correct. Checked by: exact fermion matrices (one mode plus four extra Jordan-Wigner modes standing for the odd labels; all identities exact; T(gamma)^{-1} = T(-gamma)). The adjoint statement T^* = T(-gamma) needs the Grassmann conjugation and was not tested.
QF code (exact): none. no fermionic layer.

### Fermionic coherent state, overlap, resolution
Topic: Grassmann calculus and fermionic coherent states
Appears in: C-072

The fermionic coherent state is an exponential of theta.gamma, an eigenvector of each c_j with Grassmann eigenvalue gamma_j; overlaps are Gaussian in the labels and the states resolve the identity.

$$ \psi_\gamma=e^{-\gamma^*\gamma/2}e^{\theta\cdot\gamma},\ \hat c_j\psi_\gamma=\gamma_j\psi_\gamma,\quad \langle\psi_\gamma,\psi_\alpha\rangle=\exp[\gamma^*\alpha-(\gamma^*\gamma+\alpha^*\alpha)/2],\quad \psi=\int\langle\psi_\gamma,\psi\rangle\psi_\gamma\,d\gamma^*d\gamma $$

Audit: Correct for the eigenvalue. Checked by: exact fermion matrices (extended Jordan-Wigner realization). Overlap and resolution: standard results, prior audit only.
QF code (exact): none. no fermionic layer.

### Supertrace and Berezinian
Topic: Grassmann calculus and fermionic coherent states
Appears in: C-073

The supertrace weights the trace by the parity operator and is graded-cyclic; the Berezinian is the super analogue of the determinant and satisfies Ber e^H = e^{Str H}.

$$ \operatorname{Str}\hat A=\operatorname{Tr}(\hat P\hat A),\ \hat P=(-1)^{\hat N},\ \operatorname{Str}(\hat A\hat B)=(-1)^{\pi(A)\pi(B)}\operatorname{Str}(\hat B\hat A),\quad \operatorname{Ber}\begin{pmatrix}A&B\\C&D\end{pmatrix}=\frac{\det(A-BD^{-1}C)}{\det D},\ \operatorname{Ber}(e^H)=e^{\operatorname{Str}H} $$

Audit: Correct (standard). Checked by: hand derivation (odd A, B: Tr(P A B) = -Tr(P B A) because P anticommutes with odd matrices).
QF code (exact): none. no fermionic layer.

### Fermionic Weyl translation expansion
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-074

An operator expands in fermionic translations with coefficient Str(A T(xi)).

$$ \hat A=\int A_w(\xi)\hat T(-\xi)d^2\xi,\qquad A_w(\xi)=\operatorname{Str}(\hat A\hat T(\xi)) $$

Audit: Unchecked here: depends on the source's d^2 xi normalization; prior audit only (verified finite).
QF code (exact): none. no fermionic layer.

### Normal, antinormal and s-ordered fermionic translations
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-075

Normal and antinormal fermionic translations differ from the symmetric one by the scalar Grassmann factor e^{+gamma^* gamma/2} or e^{-gamma^* gamma/2}.

$$ \hat T_N(\gamma)=e^{\hat c^\dagger\gamma}e^{-\gamma^*\hat c}=\hat T(\gamma)e^{\gamma^*\gamma/2},\quad \hat T_A(\gamma)=e^{-\gamma^*\hat c}e^{\hat c^\dagger\gamma}=\hat T(\gamma)e^{-\gamma^*\gamma/2},\quad \hat T_s=\hat T e^{s\gamma^*\gamma/2} $$

Audit: Correct. Checked by: exact fermion matrices (extended Jordan-Wigner realization).
QF code (exact): none. no fermionic layer.

### Ordinary trace from mixed fermionic symbols
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-076

The ordinary trace of a product of fermionic operators is a paired Berezin integral of one contravariant and one covariant symbol at a reflected, doubled argument, and the trace of one operator is 2^n times its symbol at zero.

$$ \operatorname{Tr}(\hat H\hat G)=2^{-n}\widehat I_\zeta\!\left[H^w(\zeta)G_w(-2\zeta)\right],\qquad \operatorname{Tr}\hat A=2^nA^w(0) $$

Audit: Correct as far as checked: Tr A = 2^n A^w(0) holds by hand for I and c^dag c (symbol alpha^* alpha + 1/2 gives 1). The reflected argument -2 zeta: prior audit only.
QF code (exact): none. no fermionic layer.

### Fermionic star product
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-077

The fermionic star product contracts odd derivatives symmetrically, and the propagator's symbol solves the star Schrodinger equation.

$$ G\star H=G\exp\!\left[\tfrac12\left(\overleftarrow\partial_{\alpha^*}\cdot\overrightarrow\partial_\alpha+\overleftarrow\partial_\alpha\cdot\overrightarrow\partial_{\alpha^*}\right)\right]H,\qquad i\partial_tU^w=H\star U^w $$

Audit: Correct (one-mode consistency: alpha * alpha^* = alpha alpha^* + 1/2 and alpha^* * alpha = alpha^* alpha + 1/2 reproduce the symbols of c c^dag and c^dag c, which sum to 1). Checked by: hand derivation.
QF code (exact): none. no fermionic layer.

### Fermi oscillator exponential and its Weyl symbol
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-078

The exponential of one fermion number operator is linear in it, and its Weyl symbol is (e^z+1)/2 + (e^z-1) alpha^* alpha.

$$ e^{z\hat c^\dagger\hat c}=\hat I+(e^z-1)\hat c^\dagger\hat c,\quad (\hat c^\dagger\hat c)^w=\alpha^*\alpha+\tfrac12,\quad [e^{z\hat c^\dagger\hat c}]^w=\frac{e^z+1}2\exp[2\tanh(z/2)\alpha^*\alpha]=\frac{e^z+1}2+(e^z-1)\alpha^*\alpha $$

Audit: Correct ((e^z + 1) tanh(z/2) = e^z - 1). Checked by: hand derivation.
QF code (exact): none. no fermionic layer.

### Fermionic partition determinant
Topic: Thermal states and traces
Appears in: C-079

Traces of fermionic quadratic exponentials are determinants: det(I + e^h) for the trace and det(I - e^h) for the supertrace.

$$ \operatorname{Tr}e^{\sum_jz_j\hat c_j^\dagger\hat c_j}=\prod_j(1+e^{z_j}),\quad \operatorname{Tr}e^{\hat c^\dagger h\hat c}=\det(I_n+e^h),\quad \operatorname{Str}e^{\hat c^\dagger h\hat c}=\det(I_n-e^h) $$

Audit: Correct. Checked by: exact fermion matrices (3 modes, complex non-Hermitian h).
QF code (exact): none. no fermionic layer.

### BdG Hamiltonian and its canonical matrix
Topic: Fermionic canonical transformations
Appears in: C-080

A general quadratic fermion Hamiltonian is encoded by a 2n-by-2n BdG matrix Sigma, whose exponential gives the Heisenberg action on (c, c^dag).

$$ \hat H=\tfrac12(\hat c^\dagger M\hat c+\hat c^\dagger L\hat c^\dagger-\hat c\bar L\hat c)-\tfrac14(\operatorname{Tr}M)\hat I,\quad \Sigma=\begin{pmatrix}M/2&L\\-\bar L&-\bar M/2\end{pmatrix},\quad e^{-it\hat H}\binom{\hat c}{\hat c^\dagger}e^{it\hat H}=e^{it\Sigma}\binom{\hat c}{\hat c^\dagger} $$

Audit: Correct. Checked by: exact fermion matrices (2 modes, random Hermitian M, antisymmetric complex L; H Hermitian; Heisenberg action matches).
QF code (exact): none. no fermionic layer.

### Fermionic pair vacuum
Topic: Fermionic canonical transformations
Appears in: C-082

A quadratic fermion propagator maps the vacuum to a pair Gaussian with antisymmetric Gamma_t and amplitude det(A_t)^{1/2}, and conjugates translations by the inverse canonical matrix.

$$ \hat U_t\Omega=\det(A_t)^{1/2}\exp\!\left[\tfrac12\hat c^\dagger\Gamma_t\hat c^\dagger\right]\Omega,\quad \Gamma_t=B_t^T(A_t^T)^{-1}=-\Gamma_t^T;\qquad \hat U_t\hat T(\alpha)\hat U_t^{-1}=\hat T(F_t^{-1}\alpha) $$

Audit: Correct (the atlas's inverse label). Checked by: exact fermion matrices for the pair vacuum (2 modes, t = 0.6); hand derivation for the label in one number mode (U c U^{-1} = e^{itm/2} c, so T(alpha) maps to T(e^{-itm/2} alpha) = T(F^{-1} alpha)).
QF code (exact): none. no fermionic layer.

### Fermionic Riccati and vacuum amplitude
Topic: Fermionic canonical transformations
Appears in: C-083

The fermionic pairing matrix obeys a Riccati equation and the vacuum amplitude a linear equation with a trace, solved by det(A_t)^{1/2}.

$$ i\dot\Gamma=\tfrac12(M\Gamma+\Gamma\bar M)+L+\Gamma\bar L\Gamma,\quad i\dot s=-\tfrac12\operatorname{Tr}(\bar L\Gamma+M/2)s,\quad s=\det(A_t)^{1/2}=\langle\Omega|\hat U_t|\Omega\rangle $$

Audit: Correct. Checked by: exact fermion matrices (2 modes, the BdG data of C-080; derivatives by central differences).
QF code (exact): none. no fermionic layer.

### Fermionic Mehler kernel
Topic: Fermionic canonical transformations
Appears in: C-084

The holomorphic kernel of a fermionic Gaussian unitary is the vacuum amplitude times a Grassmann Gaussian built from A^{-1} and B.

$$ K_{\hat U}(a^*,a)=a_0\exp\!\left[-\tfrac12\left(a^T\bar BA^{-1}a-2a^{*T}A^{-1}a+a^{*T}A^{-1}Ba^*\right)\right],\quad a_0=\langle0|\hat U|0\rangle $$

Audit: Unchecked here: needs the full Grassmann kernel convention; prior audit only (finite CAR checks).
QF code (exact): none. no fermionic layer.

### Fermionic Cayley Weyl symbol (corrected sign)
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-085

The Weyl symbol of a quadratic fermion propagator is a Grassmann Gaussian in the Cayley transform of F_t; the atlas corrects the sign of the exponent.

$$ U_t^w(\xi)=\det\!\left(\frac{F_t+I}{2}\right)^{1/2}\exp\!\left[-\binom{\xi^*}{\xi}\cdot(F_t-I)(F_t+I)^{-1}\binom\xi{\xi^*}\right]\ \text{(atlas correction: minus sign)} $$

Audit: Correct as corrected by the atlas. One mode, L = 0: F = diag(e^{2i theta}, e^{-2i theta}), theta = tm/4; the printed plus sign gives +2i tan(theta) xi^* xi, while C-078 gives cos(theta) exp(-2i tan(theta) xi^* xi). Checked by: hand derivation.
QF code (exact): none. no fermionic layer.

### Fermionic translation coefficient and Pfaffian amplitude
Topic: Gaussian integrals and Pfaffians
Appears in: C-086

The translation coefficient of a quadratic fermion propagator is a Grassmann Gaussian in the inverse Cayley transform, with amplitude the supertrace, whose square is (-1)^n det(F_t - I) and which equals an ordered Pfaffian.

$$ U_{t,w}(\gamma)=a_t\exp\!\left[-\tfrac14\binom\gamma{\gamma^*}^{T}J_F(F_t+I)(F_t-I)^{-1}\binom\gamma{\gamma^*}\right],\ a_t=\operatorname{Str}\hat U_t,\ a_t^2=(-1)^n\det(F_t-I),\ a_t=(-1)^{n(n-1)/2}p_t\operatorname{Pf}(-2J_FC_t) $$

Audit: Correct. One mode: Str U = 2i sin(theta), the Grassmann coefficient of Str(U T(gamma)) is -cos(theta), the exponent gives (i/2) cot(theta) gamma^* gamma, det(F - I) = 4 sin^2(theta), and the Pfaffian form gives cos(theta) 2i tan(theta). Checked by: hand derivation.
QF code (exact): none. no fermionic layer.

### Fermion pair rotation as SU(2)
Topic: SU(2) and passive mode mixing
Appears in: C-087

Fermion pair creation and annihilation close an su(2), so a pair rotation sends the vacuum to cos r times the vacuum plus sin r times a filled pair.

$$ \hat J_+=\hat c_1^\dagger\hat c_2^\dagger,\ \hat J_-=\hat c_2\hat c_1,\ \hat J_3=\tfrac12(\hat N_1+\hat N_2-1),\qquad \hat U\Omega=\cos r\,\Omega+\frac\xi r\sin r\,\hat c_1^\dagger\hat c_2^\dagger\Omega $$

Audit: Correct. Checked by: exact fermion matrices (2 modes: su(2) relations and the vacuum image).
QF code (exact): none. no fermionic layer.

### Mixed boson-fermion displacement
Topic: Grassmann calculus and fermionic coherent states
Appears in: C-088

A mixed boson-fermion displacement factors into its bosonic and fermionic normal-ordered parts, since the even generators of the two parts commute.

$$ \hat T(z,\gamma)=e^{-(|z|^2+\gamma^*\gamma)/2}e^{z\hat a_B^\dagger}e^{\hat c_F^\dagger\gamma}e^{-\bar z\hat a_B}e^{-\gamma^*\hat c_F} $$

Audit: Correct (each part factorizes by C-006 and C-071). Checked by: hand derivation.
QF code (exact): none. the bosonic factor is QF's `DisplacementOperator[al, Infinity, "Ordering" -> "Normal"]` (C-006), but there is no fermionic layer for the other factor.

### Super star product
Topic: Weyl and antinormal symbols, characteristic functions
Appears in: C-089

The super star product adds the antisymmetric bosonic contraction to the symmetric fermionic one.

$$ G\star H=G\exp\tfrac12\left[\overleftarrow\partial_\zeta\overrightarrow\partial_{\bar\zeta}-\overleftarrow\partial_{\bar\zeta}\overrightarrow\partial_\zeta+\overleftarrow\partial_{\alpha^*}\overrightarrow\partial_\alpha+\overleftarrow\partial_\alpha\overrightarrow\partial_{\alpha^*}\right]H $$

Audit: Convention: correct when zeta is the symbol of a (then zeta * conj(zeta) - conj(zeta) * zeta = 1 = [a, a^dag]); with zeta as the symbol of a^dag the bosonic sign flips. Checked by: hand derivation.
QF code (exact): none. QF has no symbol calculus and no fermionic layer.

### Supercommutator
Topic: Grassmann calculus and fermionic coherent states
Appears in: C-090

The supercommutator is the commutator with a sign that turns it into an anticommutator when both operators are odd.

$$ [\hat A,\hat B]_s=\hat A\hat B-(-1)^{\pi(\hat A)\pi(\hat B)}\hat B\hat A $$

Audit: Correct (definition). Checked by: hand derivation.
QF code (exact): none. QF has no graded (super) algebra; no fermionic layer.

### Commuting Gaussian integral and Fresnel factor
Topic: Gaussian integrals and Pfaffians
Appears in: C-091

The real Gaussian integral of e^{-x^T M x/2 + j^T x} is (2 pi)^{m/2} det(M)^{-1/2} e^{j^T M^{-1} j/2}, and its oscillatory limit carries the Fresnel phase e^{i pi sgn(A)/4}.

$$ \int_{\mathbb R^m}e^{-x^TMx/2+j^Tx}dx=(2\pi)^{m/2}\det(M)^{-1/2}e^{j^TM^{-1}j/2},\qquad \text{Fresnel factor } e^{i\pi\operatorname{sgn}(A)/4}|\det A|^{-1/2} $$

Conditions: M symmetric, Re M > 0; root positive on positive M.
Audit: Correct (standard). Checked by: hand derivation.
QF code (exact): none. A classical integral; QF's second-quantization layer does not apply.

### Adjoint exponential
Topic: BCH and Hadamard
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

### Berezin symbol product on the disk
Topic: Lie-group decompositions and generalized coherent states
Appears in: C-097

On the disk, the covariant symbol of a product is an integral of the two off-diagonal symbols against the squared overlap of SU(1,1) coherent states.

$$ (\hat A\hat B)_c(z,\bar z)=\int_DA_c(z,\bar w)B_c(w,\bar z)|\langle\psi_z,\psi_w\rangle|^2d\nu_n(w),\qquad |\langle\psi_z,\psi_w\rangle|^2=\left[\frac{(1-|z|^2)(1-|w|^2)}{|1-\bar zw|^2}\right]^n $$

Audit: Correct (insert the resolution of identity; overlap is C-047 with n = 2k). Checked by: hand derivation.
QF code (exact): none. QF has no SU(1,1) coherent-state symbol calculus.

### Grassmann Fourier involution
Topic: Grassmann calculus and fermionic coherent states
Appears in: C-098

The Grassmann Fourier transform is an involution in the paired convention; it maps 1 to the ordered top product of the alpha_j alpha_j^*.

$$ (\widehat{\mathcal F}f)(\alpha,\alpha^*)=\widehat I_\xi[e^{\xi^*\cdot\alpha-\alpha^*\cdot\xi}f],\quad \widehat{\mathcal F}^2=\hat I,\quad \widehat{\mathcal F}1=\alpha\alpha^*,\ \widehat{\mathcal F}\xi=\alpha,\ \widehat{\mathcal F}\xi^*=-\alpha^*,\ \widehat{\mathcal F}(\xi\xi^*)=1,\quad \widehat{\mathcal F}1=\prod_j\alpha_j\alpha_j^* $$

Audit: Correct (the atlas's corrected F1). Checked by: exact WL computation (exterior algebra: n = 1 table exact; inversion formula and F^2 = identity on all monomials for n = 1, 2; F1 for n = 2 is the ordered top product).
QF code (exact): none. no fermionic layer.

### Holomorphic Grassmann reproducing kernel
Topic: Grassmann calculus and fermionic coherent states
Appears in: C-099

The exponential Grassmann kernel e^{-(theta - gamma).gamma^*} reproduces every holomorphic Grassmann function under paired integration.

$$ \psi(\theta)=\widehat I_\gamma[\psi(\gamma)E(\theta,\gamma,\gamma^*)],\qquad E=e^{-(\theta-\gamma)\cdot\gamma^*}=\prod_j[1-(\theta_j-\gamma_j)\gamma_j^*] $$

Audit: Correct. Checked by: exact WL computation (exterior algebra, all four holomorphic monomials, n = 2).
QF code (exact): none. no fermionic layer.

## Non-formula records

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
