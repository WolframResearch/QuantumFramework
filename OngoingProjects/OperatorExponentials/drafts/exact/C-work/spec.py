# Per-entry text for drafts/exact/C-section.md, in input order (90 entries).
# Each item: (sentence, audit, cell_id or None, none_reason or match_note)
NOFERM = "no fermionic layer."
S = []
def e(sent, audit, cell=None, note=""):
    S.append((sent, audit, cell, note))

e("Momentum and position along the same axis fail to commute by a fixed imaginary constant; along different axes they commute.",
  "Correct. Checked by: QF exact (one mode, hbar = 1).", "C001",
  "Matches: [P, Q] = -i at hbar = 1; different modes commute because their field variables commute.")
e("When the commutator of two generators commutes with both, the exponential of their sum is the product of their exponentials times the scalar exponential of minus half the commutator.",
  "Correct. Checked by: QF exact.", "C002",
  "Matches: with A = u a^dag, B = v a one has [A, B] = -u v, so e^{-[A,B]/2} = e^{u v/2}.")
e("Two phase-space translations compose into the translation by the summed vector times a phase set by the symplectic area between them; reversing the order reverses that phase, and the inverse is the translation by -z.",
  "Correct. With alpha = (q+ip)/sqrt(2 hbar), alpha conj(beta) - conj(alpha) beta = -i sigma/hbar, which forces the phase. Checked by: QF exact (hbar = 1).", "C003",
  "Matches: the scalar in the first product is (al conj(be) - be conj(al))/2, which the last entry shows is -i sigma(z,z')/2; the reversed product carries the opposite scalar, so their ratio is e^{-i sigma}; D(al) D(-al) = 1.")
e("Conjugating by a translation shifts position and momentum by the classical displacement; the translation factors into a position-space kick and a momentum-generated shift times a scalar phase, and in position space it shifts the argument and multiplies by a plane wave.",
  "Correct. The position action follows by hand: e^{-iqP/hbar} u(x) = u(x-q) and e^{ipQ/hbar} multiplies by e^{ipx/hbar}, so with the phase e^{-iqp/(2 hbar)} one gets e^{ip(x-q/2)/hbar} u(x-q). Checked by: QF exact (hbar = 1) for the conjugation and the factorization; hand derivation: e^{ipx} e^{-iqp/2} u(x-q) for the position action.", "C004",
  "Matches: T Q T^{-1} = Q - q, T P T^{-1} = P - p, and e^{-iqp/2} e^{ipQ} e^{-iqP} returns exp(i(p Q - q P)) = T(z).")
e("The displacement operator equals its normal-ordered product times e^{-|alpha|^2/2}, or its antinormal product times e^{+|alpha|^2/2}, and it shifts the annihilator by -alpha.",
  "Correct. Checked by: QF exact.", "C006",
  "Matches all three statements.")
e("A coherent state is the vacuum displaced by alpha; it is an eigenvector of the annihilator with eigenvalue alpha and has amplitudes alpha^n/sqrt(n!) times e^{-|alpha|^2/2} in the number basis.",
  "Correct. Checked by: QF exact for n = 0 to 5 and for the eigenvalue; hand derivation for general n: e^{alpha a^dag}|0> = sum_n alpha^n (a^dag)^n/n! |0> = sum_n alpha^n/sqrt(n!) |n>.", "C007",
  "Matches: the amplitude differences vanish, and a e^{al a^dag} = e^{al a^dag}(a + al), so a acting on e^{al a^dag}|0> returns al times the same vector.")
e("The vacuum expectation of a translation is a Gaussian in the displacement, and two coherent states overlap as a Gaussian in their separation times a symplectic phase.",
  "Correct. |alpha|^2 = |z|^2/(2 hbar) and i Im(conj(alpha) beta) = i sigma(z,z')/(2 hbar) convert one form into the other; T(z) phi_{z'} follows from C-003 applied to T(z') acting on the vacuum. Checked by: QF exact (hbar = 1) plus hand derivation.", "C008",
  "Matches: both exponents equal the phase-space forms (last entry zero).")
e("The oscillator propagator rotates phase space rigidly: it conjugates a translation into the translation by the rotated vector and sends a coherent state to the coherent state labelled e^{-it} alpha, times the zero-point phase.",
  "Correct (e^{-it} alpha corresponds to F_t z). QF gives e^{-itN} a e^{itN} = e^{it} a, so U_t (alpha a^dag - conj(alpha) a) U_t^* is the generator of D(e^{-it} alpha); with U_t|0> = e^{-int/2}|0> this gives the state formula. Checked by: QF exact (one mode) plus hand derivation.", "C009",
  "Matches: Cosh[Sqrt[-t^2]] + I t Sinh[Sqrt[-t^2]]/Sqrt[-t^2] is cos t + i sin t = e^{it} (second entry zero).")
e("Coherent projectors integrated over phase space with weight (2 pi hbar)^{-n} give the identity, which turns traces and Hilbert-Schmidt norms into phase-space integrals of coherent expectation values.",
  "Correct. dq dp = 2 hbar d^2 alpha turns the prefactor into 1/pi per mode; the number-basis element <m|...|n> has an angular integral that vanishes for m != n and a radial part (2/n!) times the integral of r^{2n+1} e^{-r^2}, which exact integration gives as 1. Checked by: exact WL computation plus hand derivation.", None,
  "QF has no phase-space integration; `BosonicMatrixElement[{al, be}, 1, \"Basis\" -> \"Coherent\"]` gives only the overlaps in the integrand.")
e("In the Bargmann space of entire functions, creation acts as multiplication by zeta and annihilation as the derivative, so the linear and number exponentials act as a multiplication, a shift, and a dilation.",
  "Correct (Taylor shift, multiplication, Euler scaling). Checked by: hand derivation: e^{v d/dzeta} f(zeta) = sum_k v^k f^{(k)}(zeta)/k! = f(zeta + v), and zeta d/dzeta zeta^k = k zeta^k gives e^{sN} zeta^k = (e^s zeta)^k.", None,
  "QF has no function-space (Bargmann) realization; the operator statement behind the dilation, e^{sN} a e^{-sN} = e^{-s} a, is the C-094 call.")
e("The sum of products of Hermite functions weighted by w^k has a closed Gaussian form; with w = e^{-tau} and the zero-point factor it becomes the imaginary-time oscillator kernel.",
  "Correct. Checked by: exact WL computation (n = 1, symbolic x, y: the Taylor coefficients of the closed form in w agree with the Hermite-function sum through order w^10; with w = e^{-tau} and the factor e^{-tau/2} the closed form divided by the tau kernel simplifies to 1); several dimensions factor into one-dimensional kernels.", None,
  "A position-space kernel; QF has no exact position-representation call for bosonic operators.")
e("An operator acts through its Weyl symbol evaluated at the midpoint, and the symbol is recovered from the integral kernel by a Fourier transform in the difference variable.",
  "Correct (the second is the Fourier inverse of the first in u = x - y). Checked by: hand derivation.", None,
  "QF has no Weyl-symbol map for operators.")
e("With the kernel convention of C-013, Hilbert-Schmidt norms, traces of products and expectation values become phase-space integrals of symbols with the factor (2 pi hbar)^{-n}.",
  "Correct for the C-013 kernel convention. Checked by: hand derivation.", None,
  "QF has no Weyl-symbol calculus; its `WignerFunction` takes a finite-size state and returns a normalized density.")
e("A metaplectic operator transforms Weyl symbols and Wigner functions by composing them with the classical linear canonical map.",
  "Correct (the second follows from the first through C-014). Checked by: hand derivation.", None,
  "QF has no Weyl-symbol calculus.")
e("The Weyl symbol of a product is the Moyal star product, and when one factor is at most quadratic the commutator is exactly i hbar times the quantized Poisson bracket.",
  "Correct. q*p - p*q = i hbar reproduces [Q,P] = i hbar, and H = q, B = p gives i hbar {q,p} = i hbar. Checked by: hand derivation.", None,
  "QF has no Moyal product.")
e("The Weyl symbol is the anti-Wick symbol smoothed by a Gaussian of width set by hbar, and the Wick symbol is the Weyl symbol smoothed once more by the same Gaussian.",
  "Correct. The Weyl symbol of |phi_z><phi_z| is 2^n e^{-|X-z|^2/hbar}; with the (2 pi hbar)^{-n} measure this gives the kernel G; per coordinate variance hbar/2 gives e^{hbar Delta/4}. Checked by: hand derivation.", None,
  "QF's `WignerFunction`, `HusimiQFunction` and `SOrderedFunction` act on finite-size states, not on operator symbols.")
e("For a quadratic Hamiltonian the Heisenberg evolution of (Q, P) is linear and given by the classical flow matrix solving F' = J S F, which is e^{tJS} for constant S.",
  "Correct. Checked by: QF exact (one mode, symbolic S, hbar = 1): i[H, (Q, P)] = J S (Q, P), and the full conjugation U^* (Q, P) U equals e^{tJS} (Q, P) in closed form.", "C018",
  "Matches: the generator is J S and the flow difference vanishes (Simplify run in the case s12^2 > s11 s22; both sides are entire in the entries of S).")
e("A Gaussian e^{ix^T Gamma x/(2 hbar)} stays Gaussian under a quadratic propagator: its covariance moves by the linear-fractional map of the flow blocks and its amplitude picks up det(A_t + B_t Gamma)^{-1/2}.",
  "Correct. Oscillator: Gamma = i stays i and det^{-1/2} = e^{-it/2} (zero-point phase); free particle: Gamma_t = i/(1+it), the known spreading packet. Checked by: hand derivation.", None,
  "A position-space Gaussian statement; QF has no exact position-representation call for bosonic states.")
e("The linear-fractional action of a symplectic matrix keeps a complex symmetric matrix with positive imaginary part inside that domain, composes as a group action, and its determinant factor is a cocycle.",
  "Correct. Checked by: WL computation with finite 4-by-4 symplectic matrices (two degrees of freedom, F1 and F2 as matrix exponentials, complex Gamma: all three identities hold to machine precision; Gamma_F is symmetric with positive imaginary part).", None,
  "Classical matrix algebra; QF's second-quantization layer does not apply.")
e("The complex covariance of a propagated Gaussian obeys a matrix Riccati equation and its amplitude a linear equation with a trace.",
  "Correct. Forcing line: with M = A + B Gamma0, N = C + D Gamma0, d/dt(M;N) = JS(M;N) = (L^T M + G N; -K M - L N), and Gamma = N M^{-1} gives the Riccati equation; a ~ det(M)^{-1/2} gives the trace. Checked by: hand derivation, plus a WL computation with finite matrices (central-difference derivative of the linear-fractional map agrees to finite-difference accuracy).", None,
  "Classical matrix algebra; QF's second-quantization layer does not apply.")
e("The coherent-state matrix element of the propagator of a quadratic Hamiltonian is a closed Gaussian in the two labels, with det(Y)^{-1/2}, Y = A + D + i(B - C), as prefactor.",
  "Correct. Checked by: exact WL computation (one mode, S = {{9/10, 3/10}, {3/10, 3/5}}, t = 4/5, hbar = 1: H written in su(1,1) generators, checked by QF normal ordering; U = e^{-itH} disentangled as e^{c K+} e^{d K0} e^{f K-} from the exact 2-by-2 representation; its coherent elements written in closed form and compared with the atlas formula at four label pairs at 40 digits: equal, while dropping the sigma phase leaves an order-one residual; the 2-by-2 entry whose square root is taken stays away from the negative axis along the path, so the principal root is the continuous one).", None,
  "`BosonicExpOrder[Exp[x1 a**a + x2 ad**ad + x3 ad**a]]` returns the exponential unchanged, and `BosonicMatrixElement[..., \"Basis\" -> \"Coherent\"]` leaves exponentials unevaluated.")
e("The quadratic generator of a symplectic path is read off as -J F' F^{-1}, and implementers compose only up to a sign (the metaplectic double cover).",
  "Correct (from F' = J S F and J^{-1} = -J; the sign is the metaplectic double cover). Checked by: hand derivation.", None,
  "Classical matrix statement; QF has no metaplectic-lift routine.")
e("The Heisenberg action of a quadratic propagator on the annihilator is a Bogoliubov mixing u a + v a^dag whose blocks are fixed combinations of the classical flow blocks.",
  "Correct. Checked by: QF exact (one mode, symbolic S, hbar = 1): u and v read off U^* a U equal (A+D+i(C-B))/2 and (A-D+i(C+B))/2 for the blocks of e^{tJS}, and |u|^2 - |v|^2 = 1; for several modes u v^T = v u^T is the symplectic condition written in these blocks (hand).", "C024",
  "Matches: both differences vanish and |u|^2 - |v|^2 = 1.")
e("The multimode squeeze operator is the exponential of a quadratic in creators minus its adjoint; its inverse flips the sign of B, and B and the bounded coordinate W are related by tanh and the arctanh series.",
  "Correct (the scalar series sums to ArcTanh; for matrices apply it on the eigenbasis of W^*W, where W sum_k (W^*W)^k/(2k+1) = U arctanh|W|). Convention: QF's `SqueezeOperator[xi]` is this S(-xi) for one mode. Checked by: QF exact.", "C025",
  "Matches: QF's squeezer is exp((conj(xi) a^2 - xi a^dag^2)/2), the atlas S(-xi); S S^{-1} = 1; the series is ArcTanh.")
e("Squeezing acts on the column (a, a^dag) by a block matrix built from W; for one mode it gives a cosh r - e^{i theta} a^dag sinh r.",
  "Correct. Checked by: QF exact (one mode: symbolic r, theta, with the atlas S(B) entered directly; two modes: QF gives the commutators of the generator with (a1, a2, a1^dag, a2^dag) exactly, the conjugation is the exponential of that 4-by-4 matrix by C-094, and it equals the block formula for a complex symmetric B at 40 digits).", "C026",
  "Matches: one mode gives a cosh r - e^{i th} sinh r a^dag; the commutator matrix is {{0, -B}, {-conj(B), 0}}; the block residual is 0.")
e("The squeezed vacuum is the Gaussian e^{a^dag W a^dag/2} on the vacuum with normalization det(I - W^*W)^{1/4}, obtained from a matrix Riccati equation in Bargmann space; its vacuum amplitude is det(cosh|B|)^{-1/2}.",
  "Correct. Forcing line for the Riccati: the generator acting on a e^{zeta M zeta/2} gives (1/2)[zeta B zeta - Tr(B^* M) - zeta M B^* M zeta] times the Gaussian. Several modes: with the Takagi form B = V D V^T (V unitary, D >= 0), M_t = V tanh(tD) V^T solves the Riccati equation because V^T conj(V) = I reduces it to D sech^2(tD) = D - tanh(tD) D tanh(tD), and the amplitude integrates to det(cosh(tD))^{-1/2}. Checked by: QF exact (one mode) and hand derivation (n modes).", "C027",
  "Matches after xi -> -xi (QF convention): the normal-ordered squeezer is the atlas Gaussian, <0|S|0> = Sech[|xi|]^(1/2) = det(cosh|B|)^{-1/2}, and <2|S|0> equals the atlas value (difference 0).")
e("In position space the squeezed vacuum is a normalized complex Gaussian whose covariance is the Cayley transform of W.",
  "Correct (the atlas's corrected amplitude). One mode, W = tanh r: the generator is -r(x d/dx + 1/2), so S|0> = pi^{-1/4} e^{-r/2} e^{-e^{-2r}x^2/2}, and (1 - tanh^2)^{1/4}(1 + tanh)^{-1/2} = e^{-r/2}; complex W: the norm integral gives 1. Checked by: hand derivation.", None,
  "A position-space wavefunction; QF has no exact position-representation call for bosonic states.")
e("A quadratic propagator factors into a squeeze and a number-preserving rotation; the rotation multiplies the vacuum by a unit-modulus phase fixed by det(Y_t) and det(u_t u_t^dag).",
  "Correct. |det Y|^{-1/2} 2^{n/2} = |det u|^{-1/2} = det(u u^dag)^{-1/4}, so |vartheta| = 1; Y = 2 conj(u) follows from C-024. Checked by: hand derivation (consistency of C-022, C-024, C-027).", None,
  "QF's exact layer has no metaplectic factorization; `QuantumInterferometer` works on cut-off Fock spaces and is not used.")
e("Away from det(I + F) = 0, the Weyl symbol of a metaplectic operator is a Gaussian in the Cayley transform of F, with prefactor 2^n |det(I + F)|^{-1/2} and a path-dependent phase.",
  "Correct (with the atlas's added 2^n). F = I gives 1. Oscillator F = e^{tJ}: (I-F)(I+F)^{-1} = -tan(t/2) J, so the exponent is -i tan(t/2)|X|^2/hbar and the prefactor 1/cos(t/2), the known symbol of e^{-itH_os}. Checked by: hand derivation.", None,
  "QF has no Weyl-symbol calculus.")
e("The translation coefficient Tr[T(-z) R(F)] of a metaplectic operator is the Fresnel transform of its Weyl symbol: a Gaussian in the inverse Cayley transform with a signature phase.",
  "Correct. Direct trace: Tr[D(-alpha) e^{-it(N+1/2)}] = e^{-it/2} e^{-x/2} sum L_n(x) w^n = 1/(2i sin(t/2)) exp(i x cot(t/2)/2), w = e^{-it}, x = |z|^2/(2 hbar), via the Laguerre generating function; the general chart formula with c_w = 1/cos(t/2), Q_F = tan(t/2) I, sgn = 2 gives the same. Checked by: hand derivation.", None,
  "QF has no trace of an operator on the full Fock space; its exact Laguerre matrix elements of `DisplacementOperator[al]` are one input of the hand derivation, not a trace.")
e("For a Hamiltonian that is quadratic plus linear, the propagator is the classical action phase times a translation along the classical orbit, the quadratic propagator, and the inverse translation at the start time.",
  "Correct for an affine quadratic Hamiltonian. Checked by: exact WL computation (H = N + 1/2 + f Q with symbolic f, q, p, t, hbar = 1, B = 0, R_t = e^{-it(N+1/2)}, s = 0: in the Bargmann realization a -> d/dz, a^dag -> z, the right-hand side applied to an arbitrary analytic function solves i d_t psi = H psi identically, FullSimplify giving 0, with the classical orbit and delta_t computed symbolically and the product equal to the identity at t = 0); general s follows from U(t,s) = U(t,0) U(s,0)^{-1}.", None,
  "`BosonicExpOrder[Exp[-I t (ad**a + 1/2 + f (a + ad)/Sqrt[2])]]` returns the exponential unfactored.")
e("For a periodic quadratic Hamiltonian in a stable region, the one-period propagator has eigenvalues of the form e^{-i rho pi (k+1/2)}.",
  "Unchecked: needs Mathieu-equation stability data; only the constant-frequency limit (rho = omega, period pi) was checked by hand.", None,
  "QF has no Floquet or classical-stability routine.")
e("Energy-localized traces of the propagator are asymptotically sums over periodic orbits, each weighted by |det(I - P_gamma)|^{-1/2} and a Maslov phase.",
  "Unchecked: semiclassical stationary-phase asymptotics; no finite identity to test.", None,
  "QF has no semiclassical trace routine.")
e("On the quantum torus with 2 pi hbar N = 1, translations become N-by-N clock and shift matrices of order N that commute up to a root of unity.",
  "Correct for Z|j> = omega^j |j>, X|j> = |j+1>. Checked by: exact WL computation (N = 5; the torus space is finite, so these are the complete operators). The book's (6.18)-(6.19) matrices were not seen; prior audit only for the claim that they fix this orientation.", None,
  "Finite clock and shift matrices; QF's second-quantization layer does not apply.")
e("The spin-j ladder operators satisfy the su(2) relations, raise and lower m with the standard square-root factors, and give the Casimir j(j+1).",
  "Correct. Checked by: exact WL computation (j = 3/2, complete 4-by-4 spin matrices).", None,
  "Finite spin matrices; QF's second-quantization layer does not apply.")
e("A 2-by-2 su(2) matrix squares to minus its determinant, so its exponential is a cosine plus a sine term.",
  "Correct (A^2 = -I). Checked by: exact WL computation (2-by-2 matrices).", None,
  "A 2-by-2 matrix identity; QF's second-quantization layer does not apply.")
e("A unimodular 2-by-2 matrix with nonzero lower-right entry factors uniquely into upper-triangular, diagonal and lower-triangular pieces.",
  "Correct. Checked by: exact WL computation (symbolic product returns g).", None,
  "A 2-by-2 matrix identity; QF's second-quantization layer does not apply.")
e("The SU(2) displacement factors into a raising exponential, a J_3 exponential and a lowering exponential with tan r, and acting on the lowest weight it gives the normalized spin coherent state.",
  "Correct. Checked by: exact WL computation (j = 3/2 complete spin matrices, xi = 0.7 e^{0.4 i}: both identities hold to machine precision).", None,
  "Finite spin matrices; QF's second-quantization layer does not apply.")
e("Spin coherent states overlap as (1 + conj(eta) zeta)^{2j} over the norms, with modulus set by the angle between the two directions, and resolve the identity with the measure (2j+1)/pi (1+|zeta|^2)^{-2}.",
  "Correct. Checked by: exact WL computation for both overlaps (j = 3/2 complete spin matrices); resolution by hand (each diagonal element is (2j+1) C(2j,k) B(k+1, 2j+1-k) = 1).", None,
  "Finite spin matrices; QF's second-quantization layer does not apply.")
e("The trace of a convex function of a spin Hamiltonian is bounded below by the integral over its covariant symbol and above by the integral over its contravariant symbol; for J_l these symbols are j n_l and (j+1) n_l.",
  "Correct (standard Berezin-Lieb inequality and spin symbols). Checked by: hand derivation for the covariant symbol (<j,-j| rotated J_l |j,-j> = j n_l); contravariant value and the inequality from the standard result, not recomputed.", None,
  "Spin symbols; QF's second-quantization layer does not apply.")
e("As j grows, the spin coherent kernel at scale zeta/sqrt(2j) tends to the bosonic coherent kernel.",
  "Correct as a limit ((1 + conj(eta) zeta/(2j))^{2j} -> e^{conj(eta) zeta}). Checked by: hand derivation.", None,
  "A large-spin limit; QF's second-quantization layer does not apply.")
e("The su(1,1) generators close with [K-, K+] = 2 K0 and Casimir k(k-1), and act on holomorphic functions on the disk as first-order differential operators.",
  "Correct. Checked by: exact WL computation (symbolic action on an arbitrary f(z): all four relations vanish).", None,
  "An abstract differential realization; the bosonic realizations are the QF call of C-046.")
e("Squares of one mode and pair operators of two modes realize su(1,1); the one-mode Casimir is -3/16 and the two-mode Casimir is ((N1 - N2)^2 - 1)/4.",
  "Correct (k = 1/4, 3/4 and k = (1+|nu|)/2 reproduce k(k-1)). Checked by: QF exact.", "C046",
  "Matches: the one-mode Casimir is -3/16 and both one-mode commutation relations give 0; the two-mode Casimir minus ((N1-N2)^2 - 1)/4 normal-orders to 0.")
e("The SU(1,1) coherent state is the normalized exponential of the raising operator on the lowest-weight vector; its coefficients are square roots of Gamma-function ratios and its overlap is a power of (1 - conj(eta) zeta).",
  "Correct. For general k the coefficients follow by hand from K+|k,m> = sqrt((m+1)(m+2k)) |k,m+1>. In the odd oscillator sector (k = 3/4, |k,m> = |2m+1>, K+ = a^dag^2/2) QF gives the coefficients for m = 0 to 4, an exact FullSimplify shows the Fock form (1/2)^m sqrt((2m+1)!)/m! equals the Gamma form for every m, and the exact sum of Gamma(2k+m)/(Gamma(2k) m!) x^m is (1-x)^{-2k}, which gives the overlap and the normalization. Checked by: QF exact plus exact WL computation.", "C047",
  "Matches: all differences are 0 and the series sums to (1 - x)^(-2 k).")
e("A squeeze exponential factors into a pair-creation exponential, a number exponential with log cosh r, and a pair-annihilation exponential with tanh r; the same holds for two-mode pairs.",
  "Correct. Checked by: exact WL computation (2-by-2 representation of su(1,1), symbolic r and theta) and QF exact (one mode with real r and complex xi, two modes; here the atlas's S(xi) is entered as an explicit exponential, so no sign flip).", "C048",
  "Matches: with Tan[Sqrt[-|xi|^2]]/Sqrt[-|xi|^2] = Tanh[r]/r and Sec[Sqrt[-|xi|^2]] = Sech[r] (last entry), the one-mode result is e^{zeta a^dag^2/2} Sech[r]^(N+1/2) e^{-conj(zeta) a^2/2} and the two-mode result e^{zeta a1^dag a2^dag} Sech[r]^(N1+N2+1) e^{-conj(zeta) a1 a2}, zeta = xi tanh(r)/r.")
e("Squeezing maps the annihilator to (a - omega a^dag)/sqrt(1 - |omega|^2) and rotates K0 into a hyperbolic combination of K0, K+ and K-.",
  "Correct. Checked by: QF exact (symbolic r > 0 and theta, atlas S(xi) entered directly).", "C049",
  "Matches: both differences are 0.")
e("The Wigner function of the one-mode squeezed Gaussian e^{-delta x^2/(2 hbar)} is a tilted Gaussian ellipse fixed by the real and imaginary parts of delta.",
  "Correct (unnormalized Weyl symbol of the projector; integral 2 pi hbar). Forcing line: psi(q+u/2) conj(psi(q-u/2)) = exp(-Re(delta) q^2 - Re(delta) u^2/4 - i Im(delta) q u), whose u-Fourier transform gives the (p + q Im delta)^2 term. Checked by: hand derivation.", None,
  "QF's `WignerFunction` takes a finite-size state; QF has no exact Weyl symbol of an operator.")
e("For a Hamiltonian linear in the su(1,1) generators with time-dependent coefficients, the propagator is a squeeze, a K0 phase rotation and an inverse squeeze, with the squeeze coordinate solving a scalar Riccati equation.",
  "Correct (the atlas's corrected endpoint). Checked by: WL computation in the exact 2-by-2 representation of su(1,1) (constant lambda, mu; Riccati integrated numerically: U(t,0) = S(beta_t) e^{i gamma_t K0} agrees at three times to the integration accuracy; the corrected two-time form agrees; the printed form with S(-beta_t) fails by an order-one residual).", None,
  "QF has no time-dependent disentangling routine.")
e("The derivative of the exponential of a moving generator is an average of the derivative conjugated along the exponential, not simply L' e^L.",
  "Correct. Checked by: WL computation with finite 3-by-3 matrices (central difference against quadrature of the Duhamel integral).", None,
  "A general matrix identity; QF's second-quantization layer does not apply.")
e("For the oscillator with an inverse-square potential on the half line, the lowest-weight state is x^s e^{-x^2/2} with s fixed by g, and its Bargmann index is k = (2s+1)/4.",
  "Correct (energy (2s+1)/2 = 2k with K0 = H/2). Checked by: exact WL computation (symbolic eigenvalue check).", None,
  "A differential operator on L^2(R+), outside QF's Fock-space layer.")
e("The Weyl symbol of the oscillator propagator is (cos(t/2))^{-n} e^{i tan(t/2)|X|^2}, and in imaginary time (cosh(tau/2))^{-n} e^{-tanh(tau/2)|X|^2}.",
  "Correct (power -n). Check: (2 pi)^{-1} times the integral of the tau symbol is 1/(2 sinh(tau/2)) = Tr e^{-tau H_os}. Checked by: hand derivation.", None,
  "QF has no Weyl-symbol calculus.")
e("For a quadratic Hamiltonian, the phase-space integral of |phi_{F_t z}><phi_z| with the prefactor det(u_t)^{1/2} reproduces the propagator exactly.",
  "Correct (the atlas's normalized re-expression; C_t = det u^{1/2}). Checked by: exact WL computation (the C-022 data: the coherent element of the Herman-Kluk integral is a two-dimensional Gaussian integral, done in closed form with positive-definite real part, and it equals the exact disentangled matrix element at four label pairs at 40 digits; coherent elements determine the operator).", None,
  "`BosonicMatrixElement[..., \"Basis\" -> \"Coherent\"]` leaves exponentials unevaluated and QF has no phase-space integration.")
e("A centered Gaussian Weyl symbol is a density operator when V + iJ/2 >= 0; its purity is 2^{-n} det(V)^{-1/2}, the oscillator partition function is a product of 1/(2 sinh) factors, and each thermal mode has covariance (hbar/2) coth(beta hbar omega/2).",
  "Correct (vacuum V = I/2 gives purity 1; thermal purity tanh(beta omega/2) = 1/(2 nbar + 1)). Checked by: QF exact (number-basis diagonals of Q^2, P^2 and (QP+PQ)/2 for symbolic n) plus exact WL sums over the thermal weights (hbar = 1): partition function, purity, and covariance; the positivity criterion by hand (standard).", "C055",
  "Matches: the diagonals are n + 1/2, n + 1/2, 0, and all four differences are 0.")
e("A joint quadratic evolution with a Gaussian environment acts on the system by a linear map of the characteristic-function argument times a Gaussian noise factor.",
  "Correct. Forcing line: e^{ik.(Xz+Yu)} factors into system and environment parts, and the centered Gaussian environment has chi_E(m) = exp(-m^T V_E m/2). Checked by: hand derivation.", None,
  "QF has no characteristic-function calculus for bosonic channels.")
e("Smeared annihilation and creation operators satisfy the canonical commutator with the one-particle inner product, and the field operator is their real combination.",
  "Correct (standard; a(f) antilinear). Checked by: QF exact (two modes, symbolic f, g).", "C057",
  "Matches: [a(f), a^dag(g)] = conj(f).g = <f, g>.")
e("Second quantization turns the one-particle unitary e^{-ith} into the exponential of dGamma(h); a number-conserving quadratic exponential acts on creators by the one-particle matrix and maps Psi(f) to Psi(Uf).",
  "Correct. Checked by: QF exact (two modes, beam-splitter generator, symbolic theta and phi): e^G a_j^dag e^{-G} = sum_k (e^g)_{kj} a_k^dag with g the one-particle matrix; since G annihilates the vacuum, e^G Psi(f) = Psi(e^g f) follows by hand.", "C058",
  "Matches: all four coefficient differences are 0.")
e("On Fock space over any one-particle space, the Weyl operator factors in normal order, composes with the phase e^{-i Im<f,g>}, and coherent vectors overlap as e^{-|f-g|^2/2 + i Im<f,g>}; occupation numbers are Poisson with mean ||f||^2.",
  "Correct (alpha conj(beta) - conj(alpha) beta = -2i Im(conj(alpha) beta)). Checked by: QF exact (one mode, symbolic f, g; several modes give the product over an orthonormal basis since different modes commute) and hand derivation for the Poisson law (|<n|Psi(f)>|^2 = e^{-|f|^2}|f|^{2n}/n!, from C-007).", "C059",
  "Matches: the product scalar equals -i Im<f,g> and the overlap exponent equals -|f-g|^2/2 + i Im<f,g> (both differences 0).")
e("Conjugating by a Weyl operator shifts the field a(g) by <g,f> and the number operator by the linear field terms plus ||f||^2.",
  "Correct (one mode: a(g) = conj(g) a, so a(g) + <g,f> is conj(g)(a + f)). Checked by: QF exact (one mode; several modes add mode by mode).", "C060",
  "Matches: D(-al) a D(al) = a + al and D(-al) N D(al) = (a^dag + conj(al))(a + al) = N + al a^dag + conj(al) a + |al|^2.")
e("The time derivative of a moving Weyl operator is the operator times the field of the velocity plus the scalar i Im<f, f'>.",
  "Correct. Checked by: QF exact (one mode): the exponent of D(-alpha) D(alpha + e beta) is linear in e with coefficient beta a^dag - conj(beta) a + (conj(alpha) beta - alpha conj(beta))/2, and (conj(alpha) beta - alpha conj(beta))/2 = i Im(conj(alpha) beta).", "C061",
  "Matches: the difference is 0 (conjugates entered as independent symbols alc, bec).")
# fermions
e("Fermionic modes satisfy canonical anticommutation relations, so n modes span a 2^n-dimensional Fock space and a smeared annihilator has operator norm ||f||.",
  "Correct. Checked by: exact fermion matrices (3 modes) for the CAR; hand derivation for the norm (c(f)^dag c(f) has eigenvalues 0 and ||f||^2 by the CAR).", None, NOFERM)
e("A linear combination of fermionic creators and annihilators with complex coefficients squares to the scalar u.v, so its exponential is cosh plus a sinh-over-root multiple of the generator.",
  "Correct. Checked by: exact fermion matrices (2 modes, complex u, v).", None, NOFERM)
e("The Grassmann derivative obeys a graded Leibniz rule, and the right and left derivatives differ by a parity sign.",
  "Correct. Checked by: exact WL computation (exterior algebra, odd and even test elements).", None, NOFERM)
e("Berezin integration rescales by det A under a linear change of generators, and the paired convention normalizes e^{-theta^* theta} to 1.",
  "Correct as a declared convention (the atlas's repair). Checked by: exact WL computation (exterior algebra: 3-by-3 determinant rule; paired Gaussian gives 1 for n = 2).", None, NOFERM)
e("The paired Berezin integral of e^{-theta^* B theta} is det B, and linear sources complete the square to det B e^{gamma^* B^{-1} gamma}.",
  "Correct in the C-067 convention. Checked by: exact WL computation (exterior algebra, n = 2, complex non-Hermitian B, every Grassmann coefficient matches).", None, NOFERM)
e("A general Grassmann Gaussian integrates to a Pfaffian with the block-to-paired sign (-1)^{n(n-1)/2}, and sources put the inverse matrix in the exponent.",
  "Correct (with the atlas's block-to-paired sign). Checked by: exact WL computation (exterior algebra, n = 2, complex K, L, M: all coefficients match; the example Phi = -gamma^* gamma gives integral 1 and block Pfaffian -1; Pf^2 = det).", None, NOFERM)
e("Every operator on n fermionic modes is a unique combination of normal-ordered monomials.",
  "Correct (4^n monomials span the 2^n-by-2^n matrices). Checked by: hand derivation (dimension count).", None, NOFERM)
e("The fermionic translation with Grassmann labels has an even generator, so the central BCH formula factorizes it and composes it up to a scalar Grassmann phase, and it shifts c by gamma.",
  "Correct. Checked by: exact fermion matrices (one mode plus four extra Jordan-Wigner modes standing for the odd labels; all identities exact; T(gamma)^{-1} = T(-gamma)). The adjoint statement T^* = T(-gamma) needs the Grassmann conjugation and was not tested.", None, NOFERM)
e("The fermionic coherent state is an exponential of theta.gamma, an eigenvector of each c_j with Grassmann eigenvalue gamma_j; overlaps are Gaussian in the labels and the states resolve the identity.",
  "Correct for the eigenvalue. Checked by: exact fermion matrices (extended Jordan-Wigner realization). Overlap and resolution: standard results, prior audit only.", None, NOFERM)
e("The supertrace weights the trace by the parity operator and is graded-cyclic; the Berezinian is the super analogue of the determinant and satisfies Ber e^H = e^{Str H}.",
  "Correct (standard). Checked by: hand derivation (odd A, B: Tr(P A B) = -Tr(P B A) because P anticommutes with odd matrices).", None, NOFERM)
e("An operator expands in fermionic translations with coefficient Str(A T(xi)).",
  "Unchecked here: depends on the source's d^2 xi normalization; prior audit only (verified finite).", None, NOFERM)
e("Normal and antinormal fermionic translations differ from the symmetric one by the scalar Grassmann factor e^{+gamma^* gamma/2} or e^{-gamma^* gamma/2}.",
  "Correct. Checked by: exact fermion matrices (extended Jordan-Wigner realization).", None, NOFERM)
e("The ordinary trace of a product of fermionic operators is a paired Berezin integral of one contravariant and one covariant symbol at a reflected, doubled argument, and the trace of one operator is 2^n times its symbol at zero.",
  "Correct as far as checked: Tr A = 2^n A^w(0) holds by hand for I and c^dag c (symbol alpha^* alpha + 1/2 gives 1). The reflected argument -2 zeta: prior audit only.", None, NOFERM)
e("The fermionic star product contracts odd derivatives symmetrically, and the propagator's symbol solves the star Schrodinger equation.",
  "Correct (one-mode consistency: alpha * alpha^* = alpha alpha^* + 1/2 and alpha^* * alpha = alpha^* alpha + 1/2 reproduce the symbols of c c^dag and c^dag c, which sum to 1). Checked by: hand derivation.", None, NOFERM)
e("The exponential of one fermion number operator is linear in it, and its Weyl symbol is (e^z+1)/2 + (e^z-1) alpha^* alpha.",
  "Correct ((e^z + 1) tanh(z/2) = e^z - 1). Checked by: hand derivation.", None, NOFERM)
e("Traces of fermionic quadratic exponentials are determinants: det(I + e^h) for the trace and det(I - e^h) for the supertrace.",
  "Correct. Checked by: exact fermion matrices (3 modes, complex non-Hermitian h).", None, NOFERM)
e("A general quadratic fermion Hamiltonian is encoded by a 2n-by-2n BdG matrix Sigma, whose exponential gives the Heisenberg action on (c, c^dag).",
  "Correct. Checked by: exact fermion matrices (2 modes, random Hermitian M, antisymmetric complex L; H Hermitian; Heisenberg action matches).", None, NOFERM)
e("A quadratic fermion propagator maps the vacuum to a pair Gaussian with antisymmetric Gamma_t and amplitude det(A_t)^{1/2}, and conjugates translations by the inverse canonical matrix.",
  "Correct (the atlas's inverse label). Checked by: exact fermion matrices for the pair vacuum (2 modes, t = 0.6); hand derivation for the label in one number mode (U c U^{-1} = e^{itm/2} c, so T(alpha) maps to T(e^{-itm/2} alpha) = T(F^{-1} alpha)).", None, NOFERM)
e("The fermionic pairing matrix obeys a Riccati equation and the vacuum amplitude a linear equation with a trace, solved by det(A_t)^{1/2}.",
  "Correct. Checked by: exact fermion matrices (2 modes, the BdG data of C-080; derivatives by central differences).", None, NOFERM)
e("The holomorphic kernel of a fermionic Gaussian unitary is the vacuum amplitude times a Grassmann Gaussian built from A^{-1} and B.",
  "Unchecked here: needs the full Grassmann kernel convention; prior audit only (finite CAR checks).", None, NOFERM)
e("The Weyl symbol of a quadratic fermion propagator is a Grassmann Gaussian in the Cayley transform of F_t; the atlas corrects the sign of the exponent.",
  "Correct as corrected by the atlas. One mode, L = 0: F = diag(e^{2i theta}, e^{-2i theta}), theta = tm/4; the printed plus sign gives +2i tan(theta) xi^* xi, while C-078 gives cos(theta) exp(-2i tan(theta) xi^* xi). Checked by: hand derivation.", None, NOFERM)
e("The translation coefficient of a quadratic fermion propagator is a Grassmann Gaussian in the inverse Cayley transform, with amplitude the supertrace, whose square is (-1)^n det(F_t - I) and which equals an ordered Pfaffian.",
  "Correct. One mode: Str U = 2i sin(theta), the Grassmann coefficient of Str(U T(gamma)) is -cos(theta), the exponent gives (i/2) cot(theta) gamma^* gamma, det(F - I) = 4 sin^2(theta), and the Pfaffian form gives cos(theta) 2i tan(theta). Checked by: hand derivation.", None, NOFERM)
e("Fermion pair creation and annihilation close an su(2), so a pair rotation sends the vacuum to cos r times the vacuum plus sin r times a filled pair.",
  "Correct. Checked by: exact fermion matrices (2 modes: su(2) relations and the vacuum image).", None, NOFERM)
e("A mixed boson-fermion displacement factors into its bosonic and fermionic normal-ordered parts, since the even generators of the two parts commute.",
  "Correct (each part factorizes by C-006 and C-071). Checked by: hand derivation.", None,
  "the bosonic factor is QF's `DisplacementOperator[al, Infinity, \"Ordering\" -> \"Normal\"]` (C-006), but there is no fermionic layer for the other factor.")
e("The super star product adds the antisymmetric bosonic contraction to the symmetric fermionic one.",
  "Convention: correct when zeta is the symbol of a (then zeta * conj(zeta) - conj(zeta) * zeta = 1 = [a, a^dag]); with zeta as the symbol of a^dag the bosonic sign flips. Checked by: hand derivation.", None,
  "QF has no symbol calculus and no fermionic layer.")
e("The supercommutator is the commutator with a sign that turns it into an anticommutator when both operators are odd.",
  "Correct (definition). Checked by: hand derivation.", None,
  "QF has no graded (super) algebra; no fermionic layer.")
e("The real Gaussian integral of e^{-x^T M x/2 + j^T x} is (2 pi)^{m/2} det(M)^{-1/2} e^{j^T M^{-1} j/2}, and its oscillatory limit carries the Fresnel phase e^{i pi sgn(A)/4}.",
  "Correct (standard). Checked by: hand derivation.", None,
  "A classical integral; QF's second-quantization layer does not apply.")
e("Conjugation by an exponential is the exponential of the adjoint map, which computes the conjugation whenever the commutator span closes.",
  "Correct; the atlas's rejection of a global diffeomorphism claim is right (on SU(2), exp(2 pi i sigma_3) = exp(0)). Checked by: QF exact for a closing case plus hand derivation.", "C094",
  "Matches: e^{sN} a e^{-sN} = e^{-s} a = e^{ad(sN)} a, since [sN, a] = -s a.")
e("On the disk, the covariant symbol of a product is an integral of the two off-diagonal symbols against the squared overlap of SU(1,1) coherent states.",
  "Correct (insert the resolution of identity; overlap is C-047 with n = 2k). Checked by: hand derivation.", None,
  "QF has no SU(1,1) coherent-state symbol calculus.")
e("The Grassmann Fourier transform is an involution in the paired convention; it maps 1 to the ordered top product of the alpha_j alpha_j^*.",
  "Correct (the atlas's corrected F1). Checked by: exact WL computation (exterior algebra: n = 1 table exact; inversion formula and F^2 = identity on all monomials for n = 1, 2; F1 for n = 2 is the ordered top product).", None, NOFERM)
e("The exponential Grassmann kernel e^{-(theta - gamma).gamma^*} reproduces every holomorphic Grassmann function under paired integration.",
  "Correct. Checked by: exact WL computation (exterior algebra, all four holomorphic monomials, n = 2).", None, NOFERM)
assert len(S) == 90, len(S)
