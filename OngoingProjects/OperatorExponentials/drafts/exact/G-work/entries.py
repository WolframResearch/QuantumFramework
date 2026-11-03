# Per-entry prose for drafts/exact/G-section.md (index = position of the ### entry in G-section.md).
QFC = "QF exact"
WLC = "exact WL computation"
E = {}
def e(i, desc, audit, checked, match=None, none=None):
    E[i] = dict(desc=desc, audit=audit, checked=checked, match=match, none=none)

e(1, "Writing the annihilation operator as a complex combination of position and momentum turns the canonical commutator into the unit bosonic commutator, and the oscillator Hamiltonian becomes the number operator plus a zero-point half.",
  "Correct", "QF exact (with q and p rebuilt from a and a^dagger, hbar and omega symbolic)",
  "Matches: [a, a^dagger] = 1, [q, p] = i hbar, and H = hbar omega (n + 1/2).")
e(2, "Conjugating B by an exponential of A sums the nested commutators of A with B; the series closes in finitely many terms when the nested commutators vanish or repeat.",
  "Correct (standard)", "hand derivation (both sides solve dX/ds = [A, X] with X(0) = B); QF exact on two instances",
  "Matches the series: e^{s n} a^dagger e^{-s n} = e^s a^dagger (rotation closure) and e^{s a^2} a^dagger e^{-s a^2} = a^dagger + 2 s a (series stops after one commutator).")
e(3, "When the commutator of A and B commutes with both, the exponential of the sum splits into the product of exponentials times a scalar exponential of half the commutator, with the sign set by the order.",
  "Correct", "QF exact on A = x a^dagger, B = y a (so [A, B] = -x y, central)",
  "Matches: normal order gives e^{xy/2} e^{x a^dagger} e^{y a} = e^A e^B e^{-[A,B]/2}, antinormal order gives e^{-xy/2} e^{y a} e^{x a^dagger} = e^B e^A e^{[A,B]/2}.")
e(4, "The number operator lowers and raises with a and a^dagger, so a free rotation multiplies a by a phase and carries a coherent state around a circle, the zero-point term adding only a global phase.",
  "Correct", "QF exact (commutators, conjugation, and coherent-basis matrix elements; the coherent-state identity holds for every bra beta, so it holds as a vector identity)",
  "Matches: [n, a] = -a, [n, a^dagger] = a^dagger, e^{i wt n} a e^{-i wt n} = e^{-i wt} a, and <beta|e^{-i wt (n + 1/2)}|alpha> = e^{-i wt/2} <beta|alpha e^{-i wt}> for all beta.")
e(5, "The ladder operators move one step on the number basis with square-root weights, the n-photon state is the normalized n-th power of a^dagger on the vacuum, and the number states resolve the identity.",
  "Correct", "QF exact (symbolic Fock indices; norm of (a^dagger)^3|0>/sqrt(3!)); hand derivation for the completeness sum",
  "Matches.")
e(6, "Moving a through a function of the number operator shifts the argument of the function by one.",
  "Correct", "QF exact for f(n) = n^3, n^2 and e^{s n}; hand derivation for general f (both sides have the Fock element f(n) sqrt(n) between <n-1| and |n>)",
  "Matches: both polynomial differences normal-order to 0, and e^{-s n} a e^{s n} = e^s a is a f(n) = f(n+1) a with f = e^{s n}.")
e(7, "Normal-ordered powers a^dagger^r a^r equal the falling factorial of the number operator.",
  "Correct", "QF exact (operator identity for r = 3, symbolic Fock element for r = 4); hand derivation for general r (the Fock element is n!/(n-r)!)",
  "Matches: the r = 3 difference is 0 and <n|a^dagger^4 a^4|n> = n(n-1)(n-2)(n-3).")
e(8, "The two quadratures have commutator i/2, so their variances obey the 1/16 product bound, and a number state has variance (2n+1)/4 in both.",
  "Correct", "QF exact (commutator and symbolic-n variances); the bound is Robertson's inequality with this commutator",
  "Matches.")
e(9, "Different modes commute, each mode has its own unit commutator, and multimode number states are products of normalized creation powers on the joint vacuum.",
  "Correct", "QF exact (two modes, and the norm of a1^dagger^2 a2^dagger |0,0>/sqrt(2! 1!))",
  "Matches.")
e(10, "A thermal state is a normalized exponential of the number operator; its eigenvalues form a geometric distribution in the mean photon number, and the partition function carries the zero-point factor.",
  "Correct", "QF exact (Fock diagonal of e^{-x n}) and exact WL sums",
  "Matches: the eigenvalues equal nbar^m/(1+nbar)^{m+1} with nbar = 1/(e^x - 1), they sum to 1, and the partition sum equals e^{-x/2}/(1 - e^{-x}).")
e(11, "In a thermal state the r-th factorial moment is r! nbar^r, which gives the Bose-Einstein number variance nbar(1 + nbar).",
  "Correct", "QF exact (symbolic Fock elements of a^dagger^r a^r) with exact WL sums over the thermal distribution",
  "Matches for r = 1..4, <n^2> = nbar(1 + 2 nbar), and Var n = nbar(1 + nbar).")
e(12, "A Hermitian phase operator canonically conjugate to the number operator cannot exist, because its diagonal number-basis elements would have to satisfy 0 = i.",
  "Correct", "hand derivation (the m = n element gives 0 = i)",
  none="an impossibility statement; there is nothing to compute.")
e(13, "The Susskind-Glogower operator shifts the number down by one with unit weight; it is an isometry from one side only, missing exactly the vacuum projector on the other.",
  "Correct", "hand derivation (E^dagger E = sum_m |m+1><m+1|) with the QF Fock element of a",
  "Matches: the only nonzero element of (n+1)^{-1/2} a is <m|.|m+1> = 1, which is the one-sided shift.")
e(14, "A coherent state is the right eigenvector of a, obtained by displacing the vacuum, with Poisson-weighted Fock amplitudes.",
  "Correct", "QF exact (D^dagger a D = a + alpha gives a D|0> = alpha D|0>; Fock amplitudes for m = 0..8); hand derivation for general m from the eigenvalue relation",
  "Matches.")
e(15, "The creation operator has no normalizable right eigenvector: the eigenvalue equation forces every Fock component to vanish.",
  "Correct", "hand derivation: the components read 0 = lambda c_0 and sqrt(m) c_{m-1} = lambda c_m; for lambda != 0 this gives c_0 = 0 and then every c_m = 0, for lambda = 0 it gives every c_{m-1} = 0",
  none="a nonexistence statement; there is nothing to compute.")
e(16, "The displacement exponential factors into normal or antinormal order with opposite scalar factors, and its adjoint is the displacement by minus the amplitude.",
  "Correct", "QF exact",
  "Matches both ordered forms; D(alpha) D(-alpha) = 1, and since the generator is anti-Hermitian D^dagger(alpha) = D(-alpha).")
e(17, "In the normal-ordered displacement the annihilation factor acts trivially on the vacuum, so the coherent state is a normalized exponential of a^dagger on the vacuum.",
  "Correct", "QF exact (Fock amplitudes m = 0..8; antinormal-to-normal reordering and vacuum element)",
  "Matches: <m|D|0> = e^{-|alpha|^2/2} <m|e^{alpha a^dagger}|0>, and <0|e^{alpha^* a} e^{alpha a^dagger}|0> = e^{|alpha|^2}, so the norm is e^{|alpha|^2/2}.")
e(18, "Two displacements compose to the displacement by the summed amplitude times a phase, so displacements commute up to the phase e^{alpha beta^* - alpha^* beta}.",
  "Correct", "QF exact",
  "Matches all three: the merged exponent, the commutation phase, and the phase (alpha beta^* - alpha^* beta)/2 = i Im(alpha beta^*).")
e(19, "Conjugation by a displacement shifts a and a^dagger by the amplitude and its conjugate; the inverse conjugation shifts the other way.",
  "Correct", "QF exact",
  "Matches.")
e(20, "Coherent states overlap with a Gaussian in their separation and resolve the identity with measure d^2 alpha/pi.",
  "Correct", "QF exact for the overlap; exact WL integrals for the resolution (the angular integral vanishes for m != n, the radial one equals 1)",
  "Matches.")
e(21, "Between coherent states a normally ordered operator is replaced by its symbol with a^dagger -> beta^* and a -> alpha, times the overlap.",
  "Correct", "QF exact for p = 2, q = 1; hand derivation for general normal-ordered F (a acts on the ket, a^dagger on the bra)",
  "Matches.")
e(22, "A coherent state has Poisson photon statistics, factorial moments |alpha|^{2r}, and vacuum quadrature noise.",
  "Correct", "QF exact (coherent-basis expectation values and Fock amplitudes m = 0..5)",
  "Matches.")
e(23, "A linearly driven oscillator, built as a product of short displacements, is one displacement by the accumulated amplitude times an accumulated phase.",
  "Correct (structure)", "QF exact on a product of three displacements, and hand derivation from the Weyl law; the time-dependent drive u is book-specific and not rederived",
  "Matches: the product is a single displacement by x1 + x2 + x3 times a scalar phase built from the pairwise Weyl phases.")
e(24, "Projecting a state on a coherent bra gives a Gaussian times an entire function of alpha^*, on which a^dagger acts as multiplication by z and a as d/dz.",
  "Correct", "QF exact (<alpha|m> for m = 0..6) and exact WL derivatives of z^m/sqrt(m!) matching the QF ladder elements for symbolic m",
  "Matches.")
e(25, "The Fock matrix elements of an operator are Taylor coefficients of its coherent-state diagonal times e^{|alpha|^2}, with alpha and alpha^* as independent variables.",
  "Correct where the series converges", "QF exact (coherent diagonal and Fock elements for a polynomial and for e^{c n}, all m, n <= 4)",
  "Matches: every coefficient difference is 0.")
e(26, "Multiplying a coherent projector by a^dagger on the left or a on the right acts on its matrix elements as a first-order differential operator in alpha and alpha^*.",
  "Correct", "QF exact coherent overlaps with alpha and alpha^* made independent, then exact WL differentiation; the identity holds between all coherent bras and kets",
  "Matches.")
e(27, "A density operator is a P-weighted mixture of coherent projectors, normal-ordered expectations are P-averages of the symbol, Q is the coherent diagonal over pi, and normal ordering of n^2 drops a term n.",
  "Correct", "QF exact for the ordering remark; hand derivation for optical equivalence (the G-020 rule inside the P integral) and the Q definition",
  "Matches: n^2 = a^dagger^2 a^2 + n, so :n^2: differs from n^2 by n.")
e(28, "The P function of a coherent state is a delta at its amplitude; that of a number state is a derivative of the delta, defined only under integration against test functions.",
  "Correct as distributions", "exact WL computation with QF coherent overlaps: integrating P_m against <beta|alpha><alpha|gamma> by parts gives (1/m!) d^m d^m [e^{alpha alpha^*} f] at 0, compared with <beta|m><m|gamma> for m = 0..6; hand derivation for the coherent case",
  "Matches: the weak action of P_m reproduces |m><m| between all coherent states.")
e(29, "The Weyl, normal and antinormal characteristic functions are traces of the displacement in three orders; they differ only by the Gaussian factors e^{-+|lambda|^2/2}, which is the s-ordered family.",
  "Correct, including the atlas's correction of the printed Eq. (3.128c)", "QF exact (the two ordered displacements hold for any rho under the trace); hand derivation that a creation-first product is the normal function, which differs from the antinormal one by e^{|lambda|^2}",
  "Matches: D = e^{-|lambda|^2/2} e^{lambda a^dagger} e^{-lambda^* a} = e^{|lambda|^2/2} e^{-lambda^* a} e^{lambda a^dagger}, giving C_W = e^{-|lambda|^2/2} C_N = e^{|lambda|^2/2} C_A.")
e(30, "Ordered moments are read off as mixed derivatives of the matching ordered characteristic function at the origin.",
  "Correct", "hand derivation (Taylor expansion of the ordered exponentials under the trace)",
  none="a Taylor-coefficient identity under a trace; no QF call produces it, and the coherent-basis BosonicMatrixElement leaves exponentials of a or a^dagger unevaluated.")
e(31, "P and the normal characteristic function form a complex Fourier pair with kernel e^{lambda alpha^* - lambda^* alpha} and normalization pi^{-2}.",
  "Correct (sign and pi^{-2})", "exact WL computation (script wl/031.wl): the inverse transform of e^{-nbar|lambda|^2} gives e^{-|alpha|^2/nbar}/(pi nbar) and the forward transform returns e^{-nbar|lambda|^2}",
  none="QF has no P function or characteristic-function transform.")
e(32, "A thermal state has Gaussian characteristic functions with widths nbar + 1, nbar, nbar + 1/2, a Gaussian P of width nbar, and a Gaussian Q of width nbar + 1.",
  "Correct", "QF exact (coherent diagonal of e^{-x n} for Q; Fock elements <m|e^{lambda a^dagger} e^{-lambda^* a}|m> = L_m(|lambda|^2) for m = 0..5) with the exact Laguerre generating sum; C_A, C_W by the conversion factors of the G-028 entry, P by the G-031 integral",
  "Matches: Q = e^{-|alpha|^2/(1+nbar)}/(pi(1+nbar)) and C_N = e^{-nbar|lambda|^2}.")
e(33, "Coherent states have Gaussian Q and W centered at their amplitude; number states have a Poisson-shaped Q and a Laguerre-ringed W of sign (-1)^m at the origin.",
  "Correct", "QF exact (coherent overlaps, parity in the coherent basis, Fock elements of D(alpha) and D(2 alpha) for m = 0..5); W through the displaced-parity identity of the next entry",
  "Matches all four.")
e(34, "The Wigner function is 2/pi times the expectation of the parity operator displaced to alpha, and the displaced parity equals D(2 alpha) times the parity.",
  "Correct (the Fock-sum form in the atlas is the same trace)", "QF exact: Pi D(-alpha) Pi^{-1} = D(alpha) and D(alpha)^2 = D(2 alpha), so D(alpha) Pi D(-alpha) = D(alpha) D(alpha) Pi",
  "Matches.")
e(35, "Ideal photodetection measures normally ordered field correlations; for one mode the normalized r-th order coherence is the r-th factorial moment over the mean to the r-th power.",
  "Correct (definition plus the ideal-absorber argument)", "prior audit only for the detection model; the single-mode g^{(2)} case is checked exactly in the next entry",
  none="a definition of a measured quantity; QF's G1Correlation and G2Coherence act on finite Fock states only.")
e(36, "g^{(2)} is the second factorial moment over the squared mean, equal to 1 plus the excess number variance; it is 1 for coherent light, 2 for thermal light and 1 - 1/m for a number state, and undefined for the vacuum.",
  "Correct", "QF exact (coherent-basis and symbolic-m Fock moments, exact thermal sums); hand derivation for the variance form and the vacuum 0/0",
  "Matches: 1, 2, (m-1)/m, and n(n-1) = a^dagger^2 a^2.")
e(37, "The atomic raising and lowering operators close on sigma_3, and the Jaynes-Cummings Hamiltonian conserves the total excitation number of field plus atom.",
  "Correct", "exact 2x2 atomic matrices with QF normal-ordered field entries",
  "Matches: the three Pauli relations hold and [H_JC, n + |e><e|] is the zero operator matrix.")
e(38, "On resonance the Jaynes-Cummings propagator rotates each doublet {|e,m>, |g,m+1>} by the angle lambda t sqrt(m+1), which the operator-valued matrix with its square-root number functions encodes.",
  "Correct, including factor placement", "exact WL MatrixExp of the 2x2 generator on each invariant doublet (symbolic m), coupling from the QF Fock element; |g,0> is invariant (hand)",
  "Matches: every entry of the claimed operator matrix equals the doublet rotation.")
e(39, "Detuning turns each doublet into a spin-1/2 precession at the generalized Rabi frequency Omega_m, because the doubled doublet Hamiltonian squares to Omega_m^2.",
  "Correct; dressed energies hbar omega(m+1/2) +- hbar Omega_m/2 (atlas's correction of printed Eq. 4.145 holds)", "exact WL computation (script wl/039.wl): (2h)^2 = Omega_m^2 I, the closed-form propagator equals MatrixExp, and the eigenvalues are +-Omega_m/2",
  none="a 2x2 atomic block with no field operator left; nothing for QF's bosonic algebra to do.")
e(40, "In the dispersive limit the atom picks a number-dependent rotation of the field: e^{i chi t n} in the ground state and e^{-i chi t(n+1)} in the excited state.",
  "Correct (excited-state sign as in the atlas)", "hand derivation (sigma_3 = -1 on g gives -chi n; +1 on e gives chi(n+1))",
  none="the propagator is diagonal in the atom-number basis; the field rotation it implies is the G-004 check.")
e(41, "The Raman coupling lambda n sigma_1 rotates the atom by an angle proportional to the photon number, equivalently it rotates the field in opposite senses on the two sigma_1 eigenstates.",
  "Correct", "exact WL computation (script wl/041.wl): on each number sector the exact 2x2 exponential equals both stated forms, symbolic m",
  none="no spin-boson layer; the operator acts as an exact 2x2 matrix in each number sector.")
e(42, "The two-photon and two-mode Jaynes-Cummings couplings connect |e,m> to |g,m+2> and |e,m,k> to |g,m+1,k+1> with the products of ladder square roots.",
  "Correct", "QF exact (symbolic Fock indices)",
  "Matches: sqrt((m+1)(m+2)) and sqrt((m+1)(k+1)).")
e(43, "A linear two-port keeps the output commutators bosonic only if its transfer matrix is unitary.",
  "Correct", "QF exact (general symbolic 2x2 T)",
  "Matches: the output commutator matrix equals T T^dagger entry by entry.")
e(44, "The symmetric beam splitter exp(i theta J_1) mixes the two modes with cos(theta/2) transmission and i sin(theta/2) reflection.",
  "Correct", "QF exact (symbolic theta)",
  "Matches both rows.")
e(45, "The real beam splitter exp(i theta J_2) is a plain rotation of the mode pair by theta/2.",
  "Correct", "QF exact (symbolic theta)",
  "Matches both rows (c, s) and (-s, c).")
e(46, "Two modes realize SU(2) through the Schwinger generators, with Casimir J_0(J_0+1); a J_3 rotation gives the two modes opposite half-phases.",
  "Correct", "QF exact",
  "Matches: all three cyclic commutators, J_0 commutes with every J_i, J^2 = J_0(J_0+1), and e^{-i phi J_3} a e^{i phi J_3} = e^{i phi/2} a, b -> e^{-i phi/2} b.")
e(47, "Sending N photons into one port of a beam splitter produces the binomial superposition over the two outputs.",
  "Correct", "QF exact: B b^dagger B^dagger = i sin(theta/2) a^dagger + cos(theta/2) b^dagger and B|0,0> = |0,0> give B|0,N> = (i s a^dagger + c b^dagger)^N/sqrt(N!)|0,0>; Fock amplitudes for N = 3, symbolic theta; hand binomial expansion for general N",
  "Matches.")
e(48, "Passive mixing maps a product of coherent states to a product of coherent states with amplitudes transformed by the beam-splitter matrix.",
  "Correct", "QF exact (conjugation of the two-mode displacement exponent, symbolic theta), with B|0,0> = |0,0> (hand)",
  "Matches: the transformed displacement has amplitudes c alpha + i s beta on mode 1 and i s alpha + c beta on mode 2.")
e(49, "At 50:50, two single photons leave together in one port (Hong-Ou-Mandel), and N photons in each port give only even occupations with the stated weights times i^N.",
  "Correct, with the overall phase i^N (the atlas writes \"up to the overall phase i^N\")", "QF exact for N = 1, 2, 3 (B a^dagger B^dagger and B b^dagger B^dagger at theta = pi/2, then two-mode Fock amplitudes); odd components vanish",
  "Matches for N = 1, 2, 3, and the |1,1> and |3,1> amplitudes are 0.")
e(50, "A phase shifter rotates a coherent amplitude; at angle pi it is the parity operator, which squares to the identity and maps |alpha> to |-alpha>.",
  "Correct. The interferometer fringe sign remark (-cos 2 phi) is Unchecked: it depends on the book's Eq. (11.2) setup, not reproduced in the atlas; prior audit only", "QF exact",
  "Matches: <beta|e^{i phi n}|alpha> = <beta|alpha e^{i phi}>, e^{i pi n} has Fock elements (-1)^n, e^{2 pi i n} is the identity on integer n, and <beta|Pi|alpha> = <beta|-alpha>.")
e(51, "The squeeze operator is the exponential of a quadratic anti-Hermitian generator, so its adjoint is the squeeze with opposite parameter, and it squeezes the vacuum.",
  "Correct; same sign convention as QF", "QF exact",
  "Matches: QF's symbolic S(xi) has the book's exponent and S(xi) S(-xi) = 1.")
e(52, "Squeezing transforms a into the Bogoliubov combination mu a - nu a^dagger with mu = cosh r and nu = e^{i theta} sinh r.",
  "Correct", "QF exact (symbolic r, theta)",
  "Matches both lines.")
e(53, "Half the squares of the ladder operators, with half of n + 1/2, close into the su(1,1) algebra.",
  "Correct", "QF exact",
  "Matches all three commutators.")
e(54, "The squeeze operator factors into a normal-ordered product: an exponential of a^dagger^2, a power of sech r in the number operator, and an exponential of a^2.",
  "Correct", "QF exact",
  "Matches: Sech[r]^(1/2 + n) = e^{-ln(cosh r)(n + 1/2)} and tau = e^{i Arg xi} Tanh[Abs xi]; the general BosonicExpOrder form reduces to it with c1 = xi^*/2, c2 = -xi/2, sqrt(c1 c2) = i r/2, Tan[i r]/(i r) = Tanh[r]/r.")
e(55, "The squeeze operator also factors in reversed order, with the middle number exponential changing sign.",
  "Correct", "QF exact",
  "Matches: Cosh[r]^(1/2 + n) = e^{+ln(cosh r)(n + 1/2)}; the general antinormal form has Sec^(1/2 - a a^dagger), and a a^dagger = n + 1 turns it into the same Cosh^(1/2 + n) after c1 = xi^*/2, c2 = -xi/2.")
e(56, "Squeezed vacuum is annihilated by the Bogoliubov combination mu a + nu a^dagger and contains only even photon numbers with the stated amplitudes.",
  "Correct (including the exclusion of an odd solution)", "QF exact (nullifier by conjugation; Fock amplitudes for 2m = 0..10 and odd ones zero) and exact WL normalization sum",
  "Matches: S^dagger (mu a + nu a^dagger) S = a, so the combination annihilates S|0>; the even amplitudes match, odd ones are 0, and the probabilities sum to 1.")
e(57, "The quadrature at angle vartheta of squeezed vacuum has variance (cosh 2r - sinh 2r cos(2 vartheta - theta))/4.",
  "Correct", "QF exact (Bogoliubov-transformed quadrature, vacuum expectation, symbolic r, theta, vartheta)",
  "Matches.")
e(58, "The two Cartesian quadratures of squeezed vacuum are correlated by -sinh 2r sin theta/4, and the covariance determinant stays at the minimum 1/16.",
  "Correct", "QF exact",
  "Matches both.")
e(59, "Displacing squeezed vacuum shifts the first moment by alpha, adds alpha^2 to <a^2>, and adds |alpha|^2 to the mean photon number.",
  "Correct", "QF exact (S^dagger D^dagger a D S = mu a - nu a^dagger + alpha, then vacuum expectations)",
  "Matches all three.")
e(60, "Squeezing and displacement commute up to a Bogoliubov change of the displacement amplitude, with no extra phase.",
  "Correct, no extra phase", "QF exact (conjugation of a displacement by a squeeze, symbolic r, theta)",
  "Matches both exponents exactly.")
e(61, "A displaced squeezed state is the eigenstate of the Bogoliubov combination mu a + nu a^dagger with eigenvalue gamma = mu alpha + nu alpha^*, which fixes the rotated-quadrature eigenvalue.",
  "Correct (the corrected rotated eigenvalue holds)", "QF exact: (D S)^dagger (mu a + nu a^dagger) D S = a + gamma, and Y_1 + i e^{-2r} Y_2 = e^{-r - i theta/2}(mu a + nu a^dagger) as an operator identity",
  "Matches.")
e(62, "The Fock amplitudes of a displaced squeezed vacuum are Hermite polynomials in gamma over sqrt(e^{i theta} sinh 2r), times a Gaussian prefactor.",
  "Correct (principal square roots used consistently)", "QF exact normal forms of D and S and the braid e^{-alpha^* a} e^{-tau a^dagger^2/2} = e^{-tau (a^dagger - alpha^*)^2/2} e^{-alpha^* a}, so D S|0> = F(a^dagger)|0> with F(z) = e^{-|alpha|^2/2} cosh(r)^{-1/2} e^{alpha z - tau (z - alpha^*)^2/2}; the Taylor coefficients of F (m = 0..6) were compared with the Hermite form at three random exact points to 30 digits, because FullSimplify leaves a root-branch residual at m = 1 that is zero for r > 0",
  "Matches: the largest difference over m = 0..6 and the three points is zero to 30 digits.")
e(63, "The coherent overlap of squeezed vacuum is a Gaussian in beta with the e^{i theta} tanh r quadratic term, and the Q function of the displaced state is that overlap squared, centered at alpha.",
  "Correct (the atlas's correction of printed Eq. 7.82)", "exact WL computation (script wl/063.wl): the Fock sum of <beta|2m><2m|S|0> with the G-056 amplitudes gives the overlap, and |<beta - alpha|S|0>|^2/pi equals the stated Q; hand derivation of D^dagger(alpha)|beta> = phase |beta - alpha>",
  none="BosonicMatrixElement[{beta, 0}, Exp[c ad ** ad], \"Basis\" -> \"Coherent\"] returns unevaluated (exponentials in the coherent basis are a known missing feature).")
e(64, "A displaced squeezed state has a Gaussian Wigner function centered at alpha, narrowed by e^{-r} along the squeezing axis at angle theta/2 and widened by e^{r} across it.",
  "Correct (the atlas's correction of printed Eq. 7.83)", "QF exact for C_W: S^dagger D(lambda) S = D(mu lambda + nu lambda^*) and <0|D(z)|0> = e^{-|z|^2/2}; exact WL Gaussian integral of the inverse transform (script wl/064.wl, standard formula for a positive quadratic form, determinant 1)",
  "Matches: the inverse transform of C_W(lambda) = e^{lambda alpha^* - lambda^* alpha} e^{-|mu lambda + nu lambda^*|^2/2} equals the stated W exactly.")
e(65, "The number variance of a displaced squeezed vacuum is the coherent term |mu alpha - nu alpha^*|^2 plus the squeezed-vacuum term 2 sinh^2 r cosh^2 r.",
  "Correct", "QF exact",
  "Matches.")
e(66, "Under an undepleted classical pump at resonance the degenerate amplifier generates a single-mode squeeze with parameter 2 eta t.",
  "Correct", "hand derivation, (xi^* a^2 - xi a^dagger^2)/2 at xi = 2 eta t is t(eta^* a^2 - eta a^dagger^2); QF exact on the exponent",
  "Matches: QF's S(2 eta t) exponent is t eta^* a^2 - t eta a^dagger^2 for real t. The quantized-pump Hamiltonians of G-066 and G-067 have no QF closed form.")
e(67, "The two-mode squeeze mixes each mode with the conjugate of the other: a -> mu a - nu b^dagger and b -> mu b - nu a^dagger.",
  "Correct", "QF exact (symbolic r, theta)",
  "Matches both.")
e(68, "The pair operators a^dagger b^dagger, a b and half the total number plus one close into su(1,1), and the number difference commutes with all of them.",
  "Correct", "QF exact",
  "Matches all five commutators.")
e(69, "The two-mode squeeze factors in normal order and in reversed order, the middle factor being a power of cosh r in the total number plus one.",
  "Correct", "QF exact (normal and antinormal BosonicExpOrder)",
  "Matches after c1 = xi^*, c2 = -xi (sqrt(c1 c2) = i r, Tan[i r]/(i r) = Tanh[r]/r, Sec[i r] = Sech[r]); in the antinormal form a a^dagger + b b^dagger = n_a + n_b + 2 turns Sec^(-(...)) Sec into Cosh^(n_a + n_b + 1).")
e(70, "In the sum and difference modes the two-mode squeeze is a product of commuting single-mode squeezes with opposite parameters.",
  "Correct", "QF exact (generator identity a b = (c^2 - d^2)/2 and commutation of the two single-mode generators)",
  "Matches: the generators add up to the two-mode generator, they commute, and c, d^dagger commute while [c, c^dagger] = 1.")
e(71, "The twin-beam vacuum is a sum over equal photon pairs with weights (-e^{i theta} tanh r)^m/cosh r, annihilated by two Bogoliubov combinations and by the number difference.",
  "Correct", "QF exact (nullifiers by conjugation, number difference commuting with the generator, and (a^dagger b^dagger)^m/m! |0,0> = |m,m>); the expansion follows from the normal form of the G-070 entry acting on the vacuum",
  "Matches.")
e(72, "Each twin-beam mode has mean sinh^2 r and Bose-Einstein variance, the two photon numbers are perfectly correlated, the normally ordered Cauchy-Schwarz inequality fails, and the collective quadratures are squeezed by e^{-2r}.",
  "Correct; the collective quadrature variances at theta = 0 are e^{-2r}/4 and e^{2r}/4", "QF exact (all moments and both collective variances through the Bogoliubov-transformed operators)",
  "Matches all eight.")
e(73, "Tracing one twin-beam mode leaves a thermal state with tanh^2 r as the Boltzmann factor, whose entropy is the Bose-Einstein entropy of mean sinh^2 r.",
  "Correct", "QF exact amplitudes |c_m|^2 = sech^2 r tanh^{2m} r (Schmidt form, so the reduced state is diagonal with these weights) and exact WL entropy sum",
  "Matches: the weights, S = (N+1) ln(N+1) - N ln N, and tanh^2 r = N/(N+1).")
e(74, "Under an undepleted pump the nondegenerate amplifier generates a two-mode squeeze with parameter eta t, with no factor 2.",
  "Correct", "hand derivation (-i H_I t/hbar = t(eta^* a b - eta a^dagger b^dagger), no factor 1/2 in S_2)",
  none="the statement is the one-line exponent identity; the exponential it names is the two-mode squeeze checked in the G-068 and G-070 entries.")
e(75, "Self-Kerr evolution multiplies each Fock amplitude by e^{-i K t m^2}; at t = pi/K this is the parity operator.",
  "Correct", "exact WL computation (script wl/075.wl): e^{-i pi m^2} = (-1)^m for integer m",
  none="BosonicMatrixElement[{m, n}, Exp[-I K t nn ** nn]] returns unevaluated (no rule for the square of a^dagger a in an exponent).")
e(76, "At the half revival the Kerr phase is 1 on even and -i on odd photon numbers, so the operator is a combination of identity and parity and turns a coherent state into a cat.",
  "Correct; confirms the atlas's correction of printed Eq. (7.142)", "exact WL computation (script wl/076.wl): e^{-i pi m^2/2} is 1 for m = 2k and -i for m = 2k+1, the same values as e^{-i pi/4}(1 + i(-1)^m)/sqrt(2); the printed +i on odd m is wrong",
  none="same missing n^2 rule as the previous entry.")
e(77, "The normally ordered Kerr generator n(n-1) differs from n^2 by a linear term, so its evolution differs by a free rotation and is periodic with period pi/K.",
  "Correct (period pi/K)", "QF exact for the ordering; exact WL for the phases on integer m",
  "Matches: n(n-1) = a^dagger^2 a^2, e^{-i pi m(m-1)} = 1, and e^{-iKt m(m-1)} = e^{-iKt m^2} e^{iKt m}.")
e(78, "Conjugating a by any function of n multiplies it by the phase of the difference f(n) - f(n-1); for the normally ordered Kerr this gives a(t) = a e^{-2iKt(n-1)} = e^{-2iKtn} a.",
  "Correct", "QF exact Fock elements with exact WL algebra for general f (both sides are compared on their only nonzero element <n-1|.|n>)",
  "Matches both ordered forms; the form a e^{-2iKt n} without the shift is off by the factor e^{2iKt}, as the atlas warns.")
e(79, "Under the normally ordered Kerr Hamiltonian a coherent state's moments <a^m> acquire the phase e^{-iKt m(m-1)} and the factor exp(N(e^{-2imKt} - 1)), and quadrature variances follow from the low moments.",
  "Correct", "QF exact Fock elements of a^m (m = 1, 2, 3) with exact WL sums over the coherent distribution; QF exact normal-ordered X(vartheta)^2 for the variance identity (hand: Var X = <X^2> - <X>^2 with that normal form)",
  "Matches for m = 1, 2, 3; X(vartheta)^2 normal-orders to (e^{-2i vartheta} a^2 + e^{2i vartheta} a^dagger^2 + 2 a^dagger a + 1)/4, which gives the stated variance.")
e(80, "Cross-Kerr coupling rotates mode a by an angle proportional to the photon number in mode b.",
  "Correct", "QF exact on each number sector of b (n_b = m), hand derivation that the cross-Kerr operator acts there as e^{-iKtm n_a}",
  "Matches: e^{iKtm n} a e^{-iKtm n} = e^{-iKtm} a and <beta|e^{-iKtm n}|alpha> = <beta|alpha e^{-iKtm}>.")
e(81, "A superposition of |alpha> and |-alpha> with relative phase Phi has squared norm 2(1 + e^{-2|alpha|^2} cos Phi) and is an eigenstate of a^2 with eigenvalue alpha^2.",
  "Correct", "QF exact (coherent overlaps and coherent-basis a^2)",
  "Matches.")
e(82, "The photon statistics of a cat are the Poisson weights times the interference factor |1 + e^{i Phi}(-1)^m|^2 over the normalization.",
  "Correct", "QF exact Fock amplitudes of |alpha> and |-alpha> for m = 0..6, symbolic Phi; hand derivation for general m; (I +- Pi)|alpha> = |alpha> +- |-alpha> by the parity entry",
  "Matches.")
e(83, "The cat Wigner function is two Gaussians at +-alpha plus an interference fringe cos(4 alpha y + Phi) at the origin, over the cat normalization.",
  "Correct (odd-cat denominator with a minus)", "QF exact coherent overlaps, Weyl phase and parity action assembled into W = (2/pi)<C|D(2 beta) Pi|C>, then exact WL simplification for real alpha, symbolic Phi, x, y",
  "Matches: the difference simplifies to 0.")
e(84, "A cross-Kerr gate between b and c, conjugated by a 50:50 beam splitter on a and b, becomes a total-number phase times a c-controlled real rotation of a and b.",
  "Correct", "QF exact (three modes: B^dagger n_b B, the exponent identity, and the commutation of the two factors)",
  "Matches: B^dagger n_b B = (n_a + n_b)/2 - (i/2)(a^dagger b - b^dagger a), the exponent equals the stated sum, and the two generators commute.")
e(85, "The exponential of a quadrature equals its normal-ordered form times e^{y^2/8}, the Gaussian factor of the central commutator.",
  "Correct", "QF exact",
  "Matches (the shift by <X> is a scalar factor on both sides).")
e(86, "Centered quadrature moments expand into normally ordered moments with Gaussian contraction weights N!/((N-2k)! k! 8^k); for a coherent state only the full contraction survives.",
  "Correct (fourth-order coefficients 3/2 and 3/16)", "QF exact operator identity for N = 1..6 (the scalar shift by <X> leaves the commutator and hence the identity unchanged, hand); exact WL for the coherent values",
  "Matches: all six differences normal-order to 0 and (2l)!/(l! 8^l) = (2l-1)!!/4^l.")
e(87, "The squared-field quadratures close with K_3 into su(1,1) with a state-dependent commutator, giving the bound Var Q_1 Var Q_2 >= 4 <K_3>^2.",
  "Correct", "QF exact for the three commutators; the bound is Robertson's inequality with the first commutator",
  "Matches all three.")
e(88, "Pair coherent states are eigenstates of a b with fixed number difference q; their amplitudes follow from a one-step recurrence and their normalization is a modified Bessel function.",
  "Correct", "QF exact (two-mode Fock elements with symbolic m, q) and exact WL Bessel sum",
  "Matches: the recurrence c_{m+1} sqrt((m+1)(m+q+1)) = eta c_m holds, n_a - n_b = q, and the sum equals |eta|^{-q} I_q(2|eta|).")
e(89, "Frequency conversion with U = e^{-tG} is a passive SU(2) rotation of modes b and c by angle g t with the pump phase on the cross terms.",
  "Correct (the atlas's sign correction of the printed e^{+tG} holds)", "QF exact (symbolic g, phi, t)",
  "Matches both lines.")
e(90, "Subtracting a photon conditions the mean number to <n(n-1)>/<n>; adding a photon to a coherent state gives a displaced superposition of |1> and alpha^*|0>.",
  "Correct", "QF exact (a^dagger n a = n(n-1), D^dagger a^dagger D = a^dagger + alpha^*, norm <alpha|a a^dagger|alpha> = 1 + |alpha|^2); hand for the variance form",
  "Matches.")
e(91, "A collective mode summed over a band keeps a bosonic commutator scaled by C = sum |l_k|^2, and its quadratures have commutator iC/2.",
  "Correct", "QF exact for three modes with symbolic weights and frequencies; hand derivation for a general finite band",
  "Matches.")
e(92, "At low gain the down-conversion propagator on vacuum is vacuum, minus i mu times one pair, minus mu^2 times two pairs, with the norm correction -|mu|^2/2 on vacuum.",
  "Correct", "QF exact (two-mode Fock amplitudes of the first- and second-order Taylor terms of e^{-i(mu a^dagger b^dagger + mu^* a b)}; the third-order term is cubic in mu by construction)",
  "Matches: first order gives -i mu on |1,1>, second order gives -|mu|^2/2 on |0,0> and -mu^2 on |2,2>.")
e(93, "A no-jump number exponential shrinks a coherent amplitude by e^{-s} and multiplies the state by e^{-|alpha|^2(1 - e^{-2s})/2}.",
  "Correct", "QF exact (coherent-basis matrix elements)",
  "Matches both.")
e(94, "The exponential of lambda n is the normal-ordered exponential of (e^lambda - 1) n; at lambda -> -infinity this is the vacuum projector :e^{-n}:.",
  "Correct", "QF exact; exact WL binomial sum for the Fock eigenvalues, and Fock elements of :e^{-n}: for m = 0..5",
  "Matches: NormalOrdered[E^((E^lam - 1) n)], sum_k binom(m,k)(e^lam - 1)^k = e^{lam m}, and :e^{-n}: has eigenvalues 1, 0, 0, 0, 0, 0.")
e(95, "Combining the previous entry with the thermal exponential writes the thermal state as a normal-ordered Gaussian with weight 1/(1 + nbar).",
  "Correct", "QF exact (BosonicExpOrder and the substitution x = ln(1 + 1/nbar))",
  "Matches: 1 - e^{-x} = 1/(1+nbar) and e^{-x} - 1 = -1/(1+nbar).")
e(96, "Vacuum amplitude damping is a Kraus map with operators sqrt((1-T)^l/l!) T^{n/2} a^l, whose Fock components give the stated matrix-element sum and which solves the damping master equation.",
  "Correct", "QF exact Fock elements (completeness for m = 0..5, <m|a^l|m+l>^2 = (m+l)!/m! for l = 1..3) and exact WL differentiation of the component coefficients against the Lindblad component equation (symbolic m, n, l)",
  "Matches: completeness, the component amplitudes, and d/dt of each coefficient equals the master-equation right side.")
e(97, "Damping maps a coherent dyad to the dyad of shrunken amplitudes times a coherence factor that decays with the squared distance between alpha and beta.",
  "Correct", "QF exact (T^{n/2}|alpha> = e^{-(1-T)|alpha|^2/2} |sqrt(T) alpha> in the coherent basis) with the exact WL sum over Kraus indices; the C_N statement follows by optical equivalence (hand)",
  "Matches.")
e(98, "A driven damped oscillator keeps coherent states coherent; the amplitude relaxes to -2iG/gamma.",
  "Correct (factor 2 versus the printed -iG/gamma)", "QF exact coherent-basis matrix elements of the Lindblad generator applied to |beta><beta| at beta = -2iG/gamma (zero between all coherent bras and kets), and exact WL DSolve for the amplitude equation",
  "Matches: the generator vanishes on the steady state, and the ODE solution and its limit are as stated.")
e(99, "The recoil exponential is a displacement by -i eta, and the free oscillator rotation turns its generator into the time-dependent quadrature.",
  "Correct", "QF exact",
  "Matches for real eta (and BosonicExpOrder gives the same normal form for any eta); the interaction-picture exponent difference is 0.")
e(100, "To second order a fast oscillating coupling leaves a secular term linear in t proportional to [A, A^dagger], which defines the effective Hamiltonian.",
  "Correct (sign of the secular term)", "exact WL computation (script wl/100.wl): the time-ordered double integrals of e^{+-i Delta(t1 - t2)} have secular parts +-i t/Delta and those of e^{+-i Delta(t1 + t2)} are bounded, giving -i g^2 t [A, A^dagger]/Delta",
  none="a time-ordered perturbative statement; QF has no Dyson-series tool.")
e(101, "For A = i a^dagger^2 b the commutator [A, A^dagger] is a polynomial in the two number operators, so the effective Hamiltonian is a sum of self-Kerr, cross-Kerr and pump-rotation terms.",
  "Correct; the resulting H_eff follows by the G-102 formula", "QF exact",
  "Matches: the difference normal-orders to 0.")
e(102, "The trace of an operator is the integral of its coherent diagonal over d^2 alpha/pi, and the entropy of a pure bipartite state is the Shannon entropy of its squared Schmidt coefficients.",
  "Correct with the stated Schmidt premise", "QF exact coherent diagonal of e^{-x n} with an exact WL integral against the Fock trace; hand derivation (G-019 resolution) for general trace-class F and for the Schmidt statement",
  "Matches: both sides equal e^x/(e^x - 1).")
e(103, "Continuum modes have a delta-function commutator, and to first order a pump photon splits into a pair constrained by energy and momentum conservation.",
  "Correct as a distributional statement; the smeared form [a(f), a^dagger(g)] = <f, g> I is right", "hand derivation",
  none="continuum modes; QF works with finitely many modes.")
