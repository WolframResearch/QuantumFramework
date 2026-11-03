# Per-entry text for drafts/exact/P-section.md. Keys are 1-based entry positions.
# desc: what the formula states. audit: the Audit line. cell: True (use cells/NN.wl) or a
# "none" sentence. out: optional shortened Output text; match: the sentence after Output.

E = {}

E[1] = dict(
    desc="The canonical pair fixes the ladder operators and their unit commutator, and the multimode displacement is the exponential of the linear form in the ladder operators.",
    audit="Correct. Checked by: QF exact (unit commutator, and $[\\hat q,\\hat p]=i\\hbar$ with $\\hat q,\\hat p$ written through $\\hat a,\\hat a^\\dagger$); hand derivation for the multimode form: the modes commute, so the exponential splits into one-mode displacements.",
    cell=True, match="Matches: the first two entries are $[\\hat a,\\hat a^\\dagger]$ and $[\\hat q,\\hat p]$, the third is $\\hat D(\\alpha)$ for one mode.")

E[2] = dict(
    desc="When the commutator of two operators commutes with both, the product of their exponentials equals the exponential of the sum times the exponential of half the commutator.",
    audit="Correct. Checked by: QF exact.",
    cell=True, match="Matches: $[x\\hat a^\\dagger,y\\hat a]=-xy$, so the result is $e^{[A,B]/2}e^{A+B}$; the second form is the same identity solved for $e^{A+B}$.")

E[3] = dict(
    desc="Two displacements compose to the displacement by the summed amplitude times a phase set by the symplectic area, and the adjoint of a displacement is the displacement by the opposite amplitude.",
    audit="Correct. Checked by: QF exact; the adjoint follows because the generator $\\alpha\\hat a^\\dagger-\\bar\\alpha\\hat a$ is anti-Hermitian, so $\\hat D(\\alpha)^\\dagger=\\hat D(-\\alpha)$, which the second output confirms is the inverse.",
    cell=True, match="Matches: the ratio of $\\hat D(\\alpha)\\hat D(\\beta)$ to $e^{(\\alpha\\bar\\beta-\\bar\\alpha\\beta)/2}\\hat D(\\alpha+\\beta)$ is 1, and $\\hat D(\\alpha)\\hat D(-\\alpha)=\\hat I$.")

E[4] = dict(
    desc="Swapping two displacements costs twice the composition phase, the group commutator is that phase times the identity, and an ordered product of several displacements is the total displacement times the phase $\\mathrm{Im}\\sum_{j>k}\\alpha_j\\bar\\alpha_k$.",
    audit="Correct (the polygon phase means the factor $e^{i\\,\\mathrm{Im}\\sum_{j>k}\\alpha_j\\bar\\alpha_k}$). Checked by: QF exact.",
    cell=True, match="Matches: swap ratio and group commutator are both $e^{\\alpha\\bar\\beta-\\bar\\alpha\\beta}$; the three-factor product times $\\hat D(-(\\alpha_1+\\alpha_2+\\alpha_3))$ divided by $e^{i\\,\\mathrm{Im}\\,s}$ is 1, with $s=\\alpha_2\\bar\\alpha_1+\\alpha_3\\bar\\alpha_1+\\alpha_3\\bar\\alpha_2$.")

E[5] = dict(
    desc="Two displacements commute exactly when their symplectic area is a multiple of $\\pi$; a coherent-state lattice is complete when its cell area is at most $\\pi$ and incomplete above.",
    audit="Correct for the commutation condition; the completeness statements are Unchecked: cited theorem, source proof via the Weierstrass sigma function. Checked by: QF exact plus hand derivation: the group commutator is $e^{\\omega_1\\bar\\omega_2-\\bar\\omega_1\\omega_2}=e^{2i\\,\\mathrm{Im}(\\omega_1\\bar\\omega_2)}$, which is 1 exactly when $\\mathrm{Im}(\\omega_1\\bar\\omega_2)\\in\\pi\\mathbb Z$. P-014 and P2-050 agree (overcomplete below $\\pi$, complete at $\\pi$ with one deletable state).",
    cell=True, match="Matches.")

E[6] = dict(
    desc="The vacuum is killed by the annihilator, number states are normalized powers of the creator on the vacuum, and the ladder operators move one step with the square-root factors.",
    audit="Correct. Checked by: QF exact, with symbolic Fock indices.",
    cell=True, match="Matches: $\\langle m|\\hat a|n\\rangle=\\sqrt n\\,\\delta_{m,n-1}$, $\\langle m|\\hat a^\\dagger|n\\rangle=\\sqrt m\\,\\delta_{m-1,n}$, and $(\\hat a^\\dagger)^3|0\\rangle/\\sqrt{3!}=|3\\rangle$.")

E[7] = dict(
    desc="The displacement factors into a momentum kick and a position translation with a half-area phase, so on wavefunctions it shifts the argument by $Q$ and multiplies by a plane-wave phase.",
    audit="Correct. Checked by: QF exact for the operator factorization ($\\hbar=1$, $\\hat q=(\\hat a+\\hat a^\\dagger)/\\sqrt2$, $\\hat p=-i(\\hat a-\\hat a^\\dagger)/\\sqrt2$); hand derivation for the action: $e^{-iQ\\hat p}$ translates by $Q$ and $e^{iP\\hat q}$ multiplies by $e^{iPq}$, which with the phase $e^{-iPQ/2}$ gives $e^{iP(q-Q/2)}\\psi(q-Q)$.",
    cell=True, match="Matches: $e^{-iPQ/2}e^{iP\\hat q}e^{-iQ\\hat p}\\,\\hat D(\\alpha)^{-1}=1$ for real $P,Q$.")

E[8] = dict(
    desc="The overlap of two displaced copies of one fiducial state is the Weyl phase times the expectation of the relative displacement in the fiducial.",
    audit="Correct for the overlap; the resolution of identity is Unchecked here (standard Schur-orthogonality result). Checked by: QF exact for the operator identity $\\hat D(\\beta)^\\dagger\\hat D(\\alpha)=\\hat D(-\\beta)\\hat D(\\alpha)=e^{(\\alpha\\bar\\beta-\\bar\\alpha\\beta)/2}\\hat D(\\alpha-\\beta)$, which gives the overlap for every fiducial.",
    cell=True, match="Matches.")

E[9] = dict(
    desc="Conjugating by a displacement shifts the annihilator by $\\alpha$ and the quadratures by $Q$ and $P$.",
    audit="Correct. Checked by: QF exact, with $\\hat q,\\hat p$ written through $\\hat a,\\hat a^\\dagger$ and symbolic $\\hbar$.",
    cell=True, match="Matches: the shifts are $\\alpha=(Q+iP)/\\sqrt{2\\hbar}$, $Q$ and $P$.")

E[10] = dict(
    desc="A linear combination of $\\hat a$ and $\\hat a^\\dagger$ is again a canonical annihilator exactly when $|u|^2-|v|^2=1$.",
    audit="Correct. Checked by: QF exact.",
    cell=True, match="Matches.")

E[11] = dict(
    desc="The displacement in normal order carries the scalar $e^{-|\\alpha|^2/2}$ and in antinormal order carries $e^{+|\\alpha|^2/2}$.",
    audit="Correct. Checked by: QF exact (the built-in ordered forms and `BosonicExpOrder` of the generator agree).",
    cell=True, match="Matches.")

E[12] = dict(
    desc="Multiplying a displacement by a ladder operator on either side equals a first-order Wirtinger derivative of the displacement plus a linear term.",
    audit="Correct. Checked by: QF exact for the two commutation steps, and hand derivation from the normal form $\\hat D=e^{-\\alpha\\bar\\alpha/2}e^{\\alpha\\hat a^\\dagger}e^{-\\bar\\alpha\\hat a}$ with $\\alpha,\\bar\\alpha$ independent: $\\partial_\\alpha\\hat D=(\\hat a^\\dagger-\\bar\\alpha/2)\\hat D$ and $\\partial_{\\bar\\alpha}\\hat D=-\\hat D(\\hat a+\\alpha/2)$ give the first and fourth lines; $e^{-\\bar\\alpha\\hat a}\\hat a^\\dagger=(\\hat a^\\dagger-\\bar\\alpha)e^{-\\bar\\alpha\\hat a}$ and $\\hat a e^{\\alpha\\hat a^\\dagger}=e^{\\alpha\\hat a^\\dagger}(\\hat a+\\alpha)$ turn them into the second and third.",
    cell=True, match="Matches: $e^{-\\bar\\alpha\\hat a}\\hat a^\\dagger e^{\\bar\\alpha\\hat a}=\\hat a^\\dagger-\\bar\\alpha$ and $e^{-\\alpha\\hat a^\\dagger}\\hat a e^{\\alpha\\hat a^\\dagger}=\\hat a+\\alpha$.")

E[13] = dict(
    desc="A coherent state is the normalized exponential of $\\alpha\\hat a^\\dagger$ on the vacuum, it is an eigenvector of $\\hat a$, two coherent states overlap by a Gaussian, and coherent projectors resolve the identity.",
    audit="Correct. Checked by: QF exact for the expansion (symbolic $n$), eigenvalue and overlap; exact WL computation for the resolution: the angular integral kills off-diagonal terms and the radial integral $\\int_0^\\infty 2r\\,e^{-r^2}r^{2n}dr/n!$ returns 1 for every $n$.",
    cell=True, match="Matches: $\\langle n|e^{\\alpha\\hat a^\\dagger}|0\\rangle=\\alpha^n/\\sqrt{n!}$; $\\hat a e^{\\alpha\\hat a^\\dagger}|0\\rangle=\\alpha e^{\\alpha\\hat a^\\dagger}|0\\rangle$; $\\langle0|e^{\\bar\\alpha\\hat a}e^{\\beta\\hat a^\\dagger}|0\\rangle=e^{\\bar\\alpha\\beta}$, which with the two normalizations is the stated overlap.")

E[14] = dict(
    desc="In the Bargmann space $\\hat a^\\dagger$ is multiplication by $z$ and $\\hat a$ is $\\partial_z$, and a displacement acts as a weighted translation of the entire function.",
    audit="Correct, including the boundary example $e^{z^2/2}$. Checked by: QF exact plus hand derivation: $\\hat D f(\\hat a^\\dagger)|0\\rangle=e^{-|\\alpha|^2/2}e^{\\alpha\\hat a^\\dagger}e^{-\\bar\\alpha\\hat a}f(\\hat a^\\dagger)|0\\rangle=e^{-|\\alpha|^2/2}e^{\\alpha\\hat a^\\dagger}f(\\hat a^\\dagger-\\bar\\alpha)|0\\rangle$, which reads $e^{-|\\alpha|^2/2+\\alpha z}f(z-\\bar\\alpha)$. Exact WL computation: the norm series of $e^{z^2/2}$, $\\sum_k(2k)!/(4^k(k!)^2)$, diverges (`SumConvergence` returns False).",
    cell=True, match="Matches: the normal form of $\\hat D$, and $e^{-\\bar\\alpha\\hat a}\\hat a^\\dagger e^{\\bar\\alpha\\hat a}=\\hat a^\\dagger-\\bar\\alpha$.")

E[15] = dict(
    desc="The kernel that takes a position wavefunction to its Bargmann function is a Gaussian in $z$ and $q$.",
    audit="Correct. Checked by: exact WL computation: $\\sum_n z^n\\psi_n(q)/\\sqrt{n!}$ with Hermite functions, summed in closed form through the Hermite generating function $\\sum_n t^nH_n(x)/n!=e^{2tx-t^2}$ at $t=z/\\sqrt2$, $x=q/\\sqrt\\hbar$, equals the stated kernel for symbolic $\\hbar$.",
    cell="none. QF's second-quantization layer has no position representation of states.")

E[16] = dict(
    desc="The covariant symbol is the coherent expectation value and is read from the normal-ordered expansion; the contravariant symbol is the weight in the coherent-projector expansion and is read from the antinormal-ordered expansion.",
    audit="Correct. Checked by: QF exact for the covariant rule; hand derivation for the contravariant rule ($\\int\\alpha^n\\bar\\alpha^m|\\alpha\\rangle\\langle\\alpha|d^2\\alpha/\\pi=\\hat a^n(\\int|\\alpha\\rangle\\langle\\alpha|d^2\\alpha/\\pi)\\hat a^{\\dagger m}$).",
    cell=True, match="Matches: $Q_{\\hat a^{\\dagger2}\\hat a}=\\bar\\alpha^2\\alpha$.")

E[17] = dict(
    desc="The covariant symbol is the contravariant symbol smoothed with the unit Gaussian.",
    audit="Correct. Checked by: QF exact plus exact WL integral, on $\\hat a\\hat a^\\dagger$, whose contravariant symbol is $|\\beta|^2$.",
    cell=True, match="Matches: the covariant symbol from QF minus the Gaussian smoothing of $|\\beta|^2$ is 0.")

E[18] = dict(
    desc="The trace of a Gibbs exponential is bounded below by the integral of $e^{-tQ}$ and above by the integral of $e^{-tP}$; Golden-Thompson bounds the trace of an exponential of a sum by the trace of the product.",
    audit="Correct. Checked by: hand derivation for the lower bound (Jensen: $e^{-t\\langle\\alpha|\\hat H|\\alpha\\rangle}\\le\\langle\\alpha|e^{-t\\hat H}|\\alpha\\rangle$, then the resolution of identity); upper bound and Golden-Thompson cited as standard theorems.",
    cell="none. QF has no operator-inequality or trace-bound function.")

E[19] = dict(
    desc="The covariant symbol of a displacement is a Gaussian-damped plane wave in phase space, and its contravariant symbol has the opposite sign of the $|\\gamma|^2/2$ term.",
    audit="Correct. Checked by: QF exact for $Q$ (vacuum element of $e^{\\bar\\alpha\\hat a}\\hat D(\\gamma)e^{\\alpha\\hat a^\\dagger}$, times $e^{-|\\alpha|^2}$); exact WL Gaussian integral showing that smoothing the stated $P$ gives the stated $Q$.",
    cell=True, match="Matches: the ratio to $e^{-|\\gamma|^2/2+\\gamma\\bar\\alpha-\\bar\\gamma\\alpha}$ is 1.")

E[20] = dict(
    desc="A displacement has a delta-function trace, displacements are orthogonal in the Hilbert-Schmidt sense, and any operator is rebuilt from its displacement transform.",
    audit="Correct. Checked by: hand derivation: $\\operatorname{Tr}\\hat D(\\gamma)=\\sum_ne^{-|\\gamma|^2/2}L_n(|\\gamma|^2)$, and $\\sum_nL_n(x)$ is $\\pi\\delta^{(2)}$ in this measure (Glauber's completeness of displacements).",
    cell="none. A distributional trace has no exact QF counterpart.")

E[21] = dict(
    desc="Fock matrix elements of a displacement are associated Laguerre polynomials times a Gaussian and a power of the amplitude.",
    audit="Correct. Checked by: QF exact: the matrix elements of the normal-ordered product $e^{-|\\gamma|^2/2}e^{\\gamma\\hat a^\\dagger}e^{-\\bar\\gamma\\hat a}$ (a route independent of QF's built-in Laguerre formula) equal the stated Laguerre form for all $m,n\\le6$, with $\\gamma,\\bar\\gamma$ symbolic.",
    cell=True, match="Matches: all 49 differences are 0.")

E[22] = dict(
    desc="The unnormalized coherent matrix element of a displacement is a single Gaussian exponential in the two coherent labels and $\\gamma$.",
    audit="Correct. Checked by: QF exact, $e^{(|\\alpha|^2+|\\beta|^2)/2}\\langle\\alpha|\\hat D(\\gamma)|\\beta\\rangle=\\langle0|e^{\\bar\\alpha\\hat a}\\hat D(\\gamma)e^{\\beta\\hat a^\\dagger}|0\\rangle$.",
    cell=True, match="Matches: the ratio to the stated kernel is 1 (`cal` stands for $\\bar\\alpha$).")

E[23] = dict(
    desc="Each displacement matrix element is a finite sum over the number of contractions, with all phases explicit.",
    audit="Correct. Checked by: QF exact, against the matrix elements of the normal-ordered product for all $m,n\\le6$.",
    cell=True, match="Matches.")

E[24] = dict(
    desc="The Weyl operator is twice a displaced parity, parity flips the sign of a displacement, the Weyl symbol is the trace against it, and these operators are delta-orthogonal.",
    audit="Correct. Checked by: QF exact for $\\hat\\Pi\\hat D(\\alpha)\\hat\\Pi=\\hat D(-\\alpha)$ and $\\hat D(\\alpha)^2=\\hat D(2\\alpha)$, which force $2\\hat D(\\alpha)\\hat\\Pi\\hat D(-\\alpha)=2\\hat D(\\alpha)\\hat D(\\alpha)\\hat\\Pi=2\\hat D(2\\alpha)\\hat\\Pi$; hand derivation for the orthogonality ($4\\operatorname{Tr}[\\hat D(2\\alpha-2\\beta)]=4\\pi\\delta^{(2)}(2(\\alpha-\\beta))=\\pi\\delta^{(2)}(\\alpha-\\beta)$).",
    cell=True, match="Matches; the third output is the Weyl symbol of $|1\\rangle\\langle1|$, $-2\\langle1|\\hat D(2\\alpha)|1\\rangle=-2e^{-2|\\alpha|^2}(1-4|\\alpha|^2)$ at $\\alpha=x+iy$.")

E[25] = dict(
    desc="With the measure $d^2\\alpha/\\pi$, the smoothing kernels from $P$ to $W$ and from $W$ to $Q$ are $2e^{-2|\\alpha-\\beta|^2}$, and from $P$ to $Q$ it is $e^{-|\\alpha-\\beta|^2}$.",
    audit="Correct (the atlas's corrected form; the printed kernels without the factor 2 give half the result). Checked by: QF exact plus exact WL integral for $W\\to Q$ on $|1\\rangle$: the QF Weyl symbol smoothed with $2e^{-2|\\alpha-\\beta|^2}$ is exactly the QF covariant symbol $|\\alpha|^2e^{-|\\alpha|^2}$. Exact WL computation for $P\\to W$ on a thermal state: smoothing $P=e^{-|\\beta|^2/\\bar n}/\\bar n$ gives $e^{-|\\alpha|^2/(\\bar n+1/2)}/(\\bar n+1/2)$, equal to $\\operatorname{Tr}[\\hat T(\\alpha)\\hat\\rho]=2e^{-2|\\alpha|^2}\\sum_np_n(-1)^nL_n(4|\\alpha|^2)$ summed exactly.",
    cell=True, match="Matches: the first output is $Q_{|1\\rangle\\langle1|}$, the second is the smoothed Weyl symbol minus it.")

E[26] = dict(
    desc="Normal, Weyl and antinormal characteristic functions differ only by Gaussian factors $e^{\\pm|\\gamma|^2/2}$, because the ordered exponentials differ from $\\hat D(\\gamma)$ by those scalars.",
    audit="Correct. P2-035 is the same identity written relative to $\\chi_N$ (source letter $O$ for Weyl, argument $\\eta$). Checked by: QF exact, at the operator level, so it holds for every density operator.",
    cell=True, match="Matches: $e^{\\gamma\\hat a^\\dagger}e^{-\\bar\\gamma\\hat a}=e^{|\\gamma|^2/2}\\hat D(\\gamma)$ and $e^{-\\bar\\gamma\\hat a}e^{\\gamma\\hat a^\\dagger}=e^{-|\\gamma|^2/2}\\hat D(\\gamma)$.")

E[27] = dict(
    desc="Normally ordered moments are mixed Wirtinger derivatives of the normal characteristic function at the origin.",
    audit="Correct. Checked by: hand derivation; each $\\partial_{\\bar\\gamma}$ brings down $-\\hat a$ on the right of $e^{\\gamma\\hat a^\\dagger}e^{-\\bar\\gamma\\hat a}$ and each $\\partial_\\gamma$ brings down $\\hat a^\\dagger$ on the left. QF exact on a coherent state: $\\chi_N=e^{\\gamma\\bar\\beta-\\bar\\gamma\\beta}$ there, and the derivative formula reproduces QF's $\\langle\\beta|\\hat a^{\\dagger2}\\hat a|\\beta\\rangle$.",
    cell=True, match="Matches: both entries are $\\bar\\beta^2\\beta$.")

E[28] = dict(
    desc="For a state diagonal in the number basis, the Weyl characteristic function is a Gaussian times the population-weighted sum of Laguerre polynomials.",
    audit="Correct. Checked by: QF exact for $\\langle n|\\hat D(\\gamma)|n\\rangle=e^{-|\\gamma|^2/2}L_n(|\\gamma|^2)$ ($n\\le8$, from the normal-ordered product; for general $n$ see the Laguerre entry), and exact WL sum for thermal $p_n=(1-x)x^n$, which gives $\\chi_W=e^{-(\\bar n+1/2)|\\gamma|^2}$ with $\\bar n=x/(1-x)$, $y=|\\gamma|^2$.",
    cell=True, match="Matches.")

E[29] = dict(
    desc="The spin generators obey the su(2) commutators, act on $|j,m\\rangle$ with the standard square-root ladder factors, and are realized on polynomials of degree at most $2j$ by first-order differential operators.",
    audit="Correct. Checked by: exact WL computation: commutators of the differential operators on a symbolic $f(z)$ vanish identically, and the exact spin matrices for $j=1/2,1,3/2,2$ satisfy the algebra with zero residual.",
    cell="none. QF's exact algebra layer is bosonic; it has no abstract su(2) generators.")

E[30] = dict(
    desc="A compact rotation exponential factors into raising, diagonal and lowering exponentials with $\\zeta=e^{i\\phi}\\tan r$, in either order, on the chart $\\cos r>0$.",
    audit="Correct on the stated chart. Checked by: exact WL computation in the faithful 2x2 matrices with symbolic $r,\\phi$: both orders equal the left side for $0<r<\\pi/2$, and for $\\pi/2<r<\\pi$ both equal minus the left side. That sign is the central element $-I$, which acts as $(-1)^{2j}$ in spin $j$, so the identity holds in every representation on the chart and fails for half-integer $j$ beyond it.",
    cell="none. QF's exact algebra layer has no su(2) generators.")

E[31] = dict(
    desc="The spin coherent state is the normalized exponential of $\\zeta\\hat J_+$ on the lowest weight, with binomial square-root amplitudes.",
    audit="Correct. Checked by: exact WL computation with exact spin matrices ($j=3/2$ and 2, symbolic $\\zeta$), and hand derivation for general $j$: $\\hat J_+^n|j,-j\\rangle=\\sqrt{n!\\,(2j)!/(2j-n)!}\\,|j,-j+n\\rangle$.",
    cell="none. No su(2) layer in QF.")

E[32] = dict(
    desc="Spin coherent states overlap by $(1+\\bar\\zeta\\eta)^{2j}$ over the normalizations, and each is annihilated by $\\hat J_-+2\\zeta\\hat J_0-\\zeta^2\\hat J_+$.",
    audit="Correct. Checked by: exact WL computation ($j=3/2$ and 2, symbolic $\\zeta,\\eta$, zero residuals); hand derivation of the annihilator: $e^{\\zeta\\hat J_+}\\hat J_-e^{-\\zeta\\hat J_+}=\\hat J_-+2\\zeta\\hat J_0-\\zeta^2\\hat J_+$.",
    cell="none. No su(2) layer in QF.")

E[33] = dict(
    desc="Spin coherent projectors integrated with the invariant sphere measure resolve the identity.",
    audit="Correct. Checked by: exact WL computation; the diagonal element $(2j+1)\\binom{2j}{n}\\int_0^\\infty t^n(1+t)^{-2j-2}dt$ simplifies to 1 for symbolic $j$ and $n$.",
    cell="none. No su(2) layer in QF.")

E[34] = dict(
    desc="The su(1,1) generators close with $[\\hat K_-,\\hat K_+]=2\\hat K_0$, the Casimir is $k(k-1)$ in the positive discrete series, and the series has explicit ladder factors and a holomorphic realization.",
    audit="Correct. Checked by: exact WL computation: the differential realization satisfies the commutators and the Casimir on a symbolic $f(z)$, and the ladder factors give $[\\hat K_-,\\hat K_+]|k,n\\rangle=2(k+n)|k,n\\rangle$ and Casimir $k(k-1)$ for symbolic $k,n$.",
    cell="none. QF has no abstract su(1,1) generators; the oscillator realizations are checked in the single-mode and two-mode realization entries.")

E[35] = dict(
    desc="A noncompact su(1,1) exponential factors into raising, diagonal and lowering exponentials with $\\zeta=e^{i\\phi}\\tanh r$; the diagonal exponent is $+\\log(1-|\\zeta|^2)\\hat K_0$ in normal order and the opposite in antinormal order.",
    audit="Correct as corrected in the atlas; the printed source signs are wrong. Checked by: QF exact in the single-mode realization ($\\hat K_+=\\hat a^{\\dagger2}/2$, $\\hat K_-=\\hat a^2/2$, $\\hat K_0=(\\hat a^\\dagger\\hat a+1/2)/2$), where the corrected middle factor is $(1-\\tanh^2r)^{\\hat K_0}=(\\mathrm{sech}\\,r)^{\\hat a^\\dagger\\hat a+1/2}$ and the printed one would be $(\\cosh r)^{\\hat a^\\dagger\\hat a+1/2}$; and exact WL computation in the faithful 2x2 matrices ($\\hat K_+=E_{12}$, $\\hat K_-=-E_{21}$, $\\hat K_0=\\mathrm{diag}(1/2,-1/2)$, commutators verified), symbolic $r>0$, $\\phi$: corrected forms give zero residual, printed forms do not. Both sides agree at $r=0$ and vary continuously, so the 2x2 identity lifts to every positive discrete series.",
    cell=True, match="Matches the corrected signs: normal order has middle factor $\\mathrm{sech}(r)^{1/2+\\hat a^\\dagger\\hat a}$; antinormal order has $\\mathrm{sech}(r)^{1/2-\\hat a\\hat a^\\dagger}=\\cosh(r)^{\\hat a^\\dagger\\hat a+1/2}$.")

E[36] = dict(
    desc="The su(1,1) coherent state is the normalized exponential of $\\zeta\\hat K_+$ on the lowest weight, with Pochhammer square-root amplitudes, and it is normalized for $|\\zeta|<1$.",
    audit="Correct. Checked by: exact WL computation: the product of ladder factors $\\prod_{j<n}\\sqrt{(j+1)(2k+j)}/n!$ equals $\\sqrt{\\Gamma(n+2k)/(n!\\Gamma(2k))}$ for symbolic $k,n$, and $\\sum_n\\Gamma(n+2k)x^n/(n!\\Gamma(2k))=(1-x)^{-2k}$ exactly.",
    cell="none. No su(1,1) layer in QF.")

E[37] = dict(
    desc="Two su(1,1) coherent states overlap by $(1-\\bar\\zeta\\eta)^{-2k}$ times the normalizations, and each is annihilated by $\\hat K_--2\\zeta\\hat K_0+\\zeta^2\\hat K_+$.",
    audit="Correct. Checked by: exact WL computation: the overlap series sums to $(1-\\bar\\zeta\\eta)^{-2k}$, and the coefficient of $|k,n\\rangle$ in the annihilator equation vanishes for symbolic $k,n$ (neighbouring amplitudes differ by $\\zeta\\sqrt{(n+2k)/(n+1)}$, so the three terms give $\\zeta(n+2k)-2\\zeta(k+n)+\\zeta n=0$).",
    cell="none. No su(1,1) layer in QF.")

E[38] = dict(
    desc="su(1,1) coherent projectors integrated with the invariant disk measure resolve the identity only for $k>1/2$.",
    audit="Correct, including the restriction $k>1/2$. Checked by: exact WL computation; the diagonal element reduces to $(2k-1)\\Gamma(n+2k)/(n!\\Gamma(2k))\\int_0^1t^n(1-t)^{2k-2}dt$, which the kernel returns as 1 under the condition $k>1/2$ (the integral diverges otherwise).",
    cell="none. No su(1,1) layer in QF.")

E[39] = dict(
    desc="Half the squares of the ladder operators and half the shifted number operator realize su(1,1) on one mode, with Casimir $-3/16$, so the even and odd sectors carry $k=1/4$ and $3/4$.",
    audit="Correct; $k(k-1)=-3/16$ for $k=1/4$ and $3/4$. Checked by: QF exact.",
    cell=True, match="Matches: $[\\hat K_0,\\hat K_+]=\\hat K_+$, $[\\hat K_-,\\hat K_+]=2\\hat K_0$, Casimir $-3/16$.")

E[40] = dict(
    desc="The one-mode squeeze operator in normal order is a pair-creation exponential, a number-diagonal factor $(\\mathrm{sech}\\,r)^{\\hat a^\\dagger\\hat a+1/2}$ and a pair-annihilation exponential, so its vacuum amplitude is $(\\cosh r)^{-1/2}$.",
    audit="Correct. Checked by: QF exact, both the built-in ordered form and `BosonicExpOrder` of the generator.",
    cell=True, match="Matches: with $\\xi=re^{i\\theta}$ the two ordered forms coincide and the vacuum element is $\\mathrm{sech}(r)^{1/2}$.")

E[41] = dict(
    desc="The pair operators of two modes realize su(1,1), and the Casimir is fixed by the number difference, so a sector with difference $d$ has $k=(1+|d|)/2$.",
    audit="Correct; $k(k-1)=(d^2-1)/4$ gives $k=(1+|d|)/2$. Checked by: QF exact.",
    cell=True, match="Matches: both commutators and the Casimir identity reduce to 0.")

E[42] = dict(
    desc="The two-mode squeeze operator in normal order is a pair-creation exponential, $(\\mathrm{sech}\\,r)^{\\hat N_a+\\hat N_b+1}$ and a pair-annihilation exponential; on the vacuum it gives a pair state with geometric pair-number distribution.",
    audit="Correct. P-032 writes the pair state with $\\hat a^\\dagger\\hat b^\\dagger$, P2-019 with $\\hat b^\\dagger\\hat a^\\dagger$ (same operator). Checked by: QF exact for the factorization; hand derivation for the pair state: $(\\hat a^\\dagger\\hat b^\\dagger)^n|0,0\\rangle=n!\\,|n,n\\rangle$ (QF confirms $n=3$), so $e^{\\zeta\\hat a^\\dagger\\hat b^\\dagger}|0,0\\rangle=\\sum_n\\zeta^n|n,n\\rangle$ with squared norm $(1-|\\zeta|^2)^{-1}$ and $W_n=(1-\\rho)\\rho^n$.",
    cell=True, match="Matches for real $r>0$.")

E[43] = dict(
    desc="A linear canonical map of $N$ bosonic modes preserves the indefinite form $K$, which fixes the block relations among $U$ and $V$, the inverse, and makes $Z=U^{-1}V$ symmetric and a strict contraction.",
    audit="Correct. Checked by: exact WL computation, $N=2$, with an exact group element obtained as the Cayley transform $(I+L/2)(I-L/2)^{-1}$ of a rational generator $L$ of the P-035 form: all ten relations hold with zero residual. Forcing line for the last one: $ZZ^\\dagger=U^{-1}(UU^\\dagger-I)U^{-\\dagger}=I-(U^\\dagger U)^{-1}$.",
    cell="none. QF has no function that tests Bogoliubov coefficient matrices.")

E[44] = dict(
    desc="A quadratic bosonic generator corresponds to a block matrix with anti-Hermitian number block and symmetric pairing block.",
    audit="Correct. Checked by: hand derivation, $LK+KL^\\dagger$ has blocks $C+C^\\dagger$ and $A^T-A$, so it vanishes exactly when $C^\\dagger=-C$ and $A^T=A$; exact WL computation confirms that block form for a general symbolic $L$.",
    cell="none. QF has no sp(2N) coefficient-matrix layer.")

E[45] = dict(
    desc="The bilinears $\\hat E,\\hat F,\\hat G$ close under commutators with the stated structure constants.",
    audit="Correct. Checked by: QF exact, all 16 index choices over two modes, six relations each.",
    cell=True, match="Matches: 96 results, all 0.")

E[46] = dict(
    desc="A general bosonic group operator factors into pair-creation, number-conserving and pair-annihilation exponentials; on the vacuum this leaves a pair-creation Gaussian, which is annihilated by $\\hat{\\mathbf a}+Z\\hat{\\mathbf a}^\\dagger$.",
    audit="Correct for the state and annihilator; the triangular factorization is structural (no coefficient functions given) and Unchecked beyond its single-mode instance. Checked by: QF exact, one mode and two modes with symmetric $Z$.",
    cell=True, match="Matches: $\\hat a_i e^{-\\frac12\\hat{\\mathbf a}^{\\dagger T}Z\\hat{\\mathbf a}^\\dagger}=e^{-\\frac12\\hat{\\mathbf a}^{\\dagger T}Z\\hat{\\mathbf a}^\\dagger}(\\hat a_i-\\sum_jZ_{ij}\\hat a_j^\\dagger)$, so on the vacuum $(\\hat a_i+\\sum_jZ_{ij}\\hat a_j^\\dagger)|Z\\rangle=0$.")

E[47] = dict(
    desc="A pair-creation Gaussian is normalized by $\\det(I-ZZ^\\dagger)^{1/4}$, and two of them overlap by the normalizations over $\\det(I-WZ^\\dagger)^{1/2}$.",
    audit="Correct. Checked by: QF exact: the one-mode vacuum element in closed form, and for two modes the vacuum element $\\langle0|e^{-\\frac12\\hat{\\mathbf a}^T\\bar Z\\hat{\\mathbf a}}e^{-\\frac12\\hat{\\mathbf a}^{\\dagger T}W\\hat{\\mathbf a}^\\dagger}|0\\rangle$ expanded in powers of $\\bar Z,W$ matches the series of $\\det(I-W\\bar Z)^{-1/2}$ through third order with symbolic symmetric $Z,W$. All orders follow from the Bargmann Gaussian integral $\\int\\prod d^2w/\\pi\\,e^{-|w|^2-\\frac12\\bar w^T\\bar Z\\bar w-\\frac12w^TWw}=\\det(I-W\\bar Z)^{-1/2}$ (standard).",
    cell=True, match="Matches: one mode gives $(1-\\bar zw)^{-1/2}$; the two-mode series differences are 0 at orders 1, 2, 3.")

E[48] = dict(
    desc="Canonical maps act on the symmetric Gaussian parameter by a matrix fractional-linear map, and the analytic representation carries a determinant multiplier.",
    audit="Correct for the map; the determinant multiplier is Unchecked (its sign and branch depend on the source's function-space convention). Checked by: exact WL computation ($N=2$, two exact Cayley group elements, rational $Z$): the image is symmetric, $I-Z'Z'^\\dagger$ has positive eigenvalues, and $g_1\\cdot(g_2\\cdot Z)=(g_2g_1)\\cdot Z$ with zero residual, which is the composition order that makes $f\\mapsto f(g\\cdot Z)$ a representation.",
    cell="none. QF has no Gaussian-parameter action.")

E[49] = dict(
    desc="On a pure bosonic Gaussian the pairing matrix $A$ and the number matrix $B$ are rational functions of $Z$ and obey quadratic purity constraints.",
    audit="Convention: on the two-mode state $|Z\\rangle$ ($\\hbar=1$), $A_{ij}=\\langle\\hat a_i\\hat a_j\\rangle$ and $B_{ij}=\\langle\\hat a_j^\\dagger\\hat a_i\\rangle+\\delta_{ij}/2$; reading $B_{ij}$ as $\\langle\\hat a_i^\\dagger\\hat a_j\\rangle+\\delta_{ij}/2$ (that is, $\\langle\\hat X^i{}_j\\rangle$) gives the transpose, which differs. The four matrix constraints are Correct. Checked by: exact WL computation with an exact rational complex symmetric $Z$: the second moments are solved exactly from the annihilator relation $\\hat a_i|Z\\rangle=-\\sum_kZ_{ik}\\hat a_k^\\dagger|Z\\rangle$ (QF exact in the previous entry) and the CCR, then compared with the formulas; the constraints hold with zero residual.",
    cell="none. `BosonicMatrixElement` takes polynomials between Fock states; it does not give moments of a multimode Gaussian state, and `CovarianceMatrix` works only on cut-off states.")

E[50] = dict(
    desc="A fermionic canonical map keeps the anticommutators exactly when the doubled coefficient matrix is unitary, and $Z=U^{-1}V$ is then skew-symmetric.",
    audit="Correct. Checked by: exact fermion matrices (3 modes), with an exact unitary group element from the Cayley transform of a rational anti-Hermitian generator of block form: the new operators satisfy the CAR exactly and all block relations hold with zero residual.",
    cell="none. no fermionic layer.")

E[51] = dict(
    desc="Fermionic bilinears close under commutators with structure constants that differ from the bosonic table in the $[\\hat F,\\hat G]$ signs.",
    audit="Correct. Checked by: exact fermion matrices (2 modes), all 16 index choices, six relations: every residual exactly 0.",
    cell="none. no fermionic layer.")

E[52] = dict(
    desc="A fermionic pair-exponential state with skew $Z$ is normalized by $\\det(I+ZZ^\\dagger)^{-1/4}$, is annihilated by $\\hat{\\mathbf c}+Z\\hat{\\mathbf c}^\\dagger$, and two such states overlap by $\\det(I+WZ^\\dagger)^{1/2}$ over the normalizations.",
    audit="Correct. Checked by: exact fermion matrices (4 modes), exact rational complex skew $Z,W$: the squared norm squared equals $\\det(I+ZZ^\\dagger)$, the annihilator residual is exactly 0, and the unnormalized overlap equals the principal square root of $\\det(I+WZ^\\dagger)$ exactly (the atlas's branch caveat stands for other parameters).",
    cell="none. no fermionic layer.")

E[53] = dict(
    desc="Fermionic canonical maps act on skew $Z$ by the same fractional-linear map as bosons, with determinant exponent $+1/2$.",
    audit="Unchecked: the determinant exponent depends on the source's function-space realization, which is not reconstructed here. The skew-symmetry of the image follows by the same algebra as the bosonic map with signs from $UV^T=-VU^T$.",
    cell="none. no fermionic layer.")

E[54] = dict(
    desc="On a pure fermionic Gaussian the skew pairing matrix $A$ and the number matrix $B$ are rational in $Z$ and obey the same quadratic constraints as for bosons with the signs changed.",
    audit="Convention: on the four-mode state of the pair-exponential entry ($\\hbar=1$), $A_{ij}=\\langle\\hat c_j\\hat c_i\\rangle=-\\langle\\hat c_i\\hat c_j\\rangle$ and $B_{ij}=\\langle\\hat c_j^\\dagger\\hat c_i\\rangle-\\delta_{ij}/2$; read as $\\langle\\hat X_{ij}\\rangle$ with $\\hat X_{ij}=\\hat c_i\\hat c_j$, the sign of $A$ is wrong. The constraints are Correct. Checked by: exact fermion matrices (4 modes), exact rational $Z$, all residuals exactly 0.",
    cell="none. no fermionic layer.")

E[55] = dict(
    desc="For a general semisimple group, coherent states are exponentials of positive-root generators on a lowest weight, and the unitary displacement uses positive and negative roots.",
    audit="Unchecked: definitions on a general semisimple group; the conversion from $\\xi$ to $\\zeta$ is not given. Rank-one instances are the SU(2) and SU(1,1) disentangling entries.",
    cell="none. QF has no root-space layer.")

E[56] = dict(
    desc="The determinant is the exponential of the trace of the logarithm, and $-\\log\\det(I-X)$ is the trace power series, convergent for spectral radius below 1.",
    audit="Correct (the kernel $\\psi_Z$ is a definition). Checked by: exact WL computation: for a fully symbolic $3\\times3$ matrix the series of $-\\log\\det(I-tX)$ agrees with $\\sum_mt^m\\operatorname{Tr}X^m/m$ through $t^5$.",
    cell="none. Not an operator identity.")

E[57] = dict(
    desc="Real symplectic matrices act on the Siegel upper half-space by a fractional-linear map, and their blocks obey three quadratic relations.",
    audit="Correct. Checked by: exact WL computation: an exact rational real symplectic $4\\times4$ matrix (Cayley transform of $JS$, $S$ symmetric) satisfies the block relations with zero residual, and the image of a point of the half-space is symmetric with positive-definite imaginary part.",
    cell="none. Not an operator identity.")

E[58] = dict(
    desc="The covariant symbol is the coherent expectation value divided by the overlap, and its off-diagonal version uses two different coherent labels.",
    audit="Correct (definitions; the continuation statement holds because the off-diagonal symbol is holomorphic in one label and antiholomorphic in the other). Checked by: hand derivation; QF exact for a flat-space instance.",
    cell=True, match="Matches: the flat off-diagonal symbol of $\\hat a^\\dagger\\hat a$ is $\\bar\\alpha\\beta$, the analytic continuation of $|\\alpha|^2$; QF has no disk coherent states.")

E[59] = dict(
    desc="On the disk, the Berezin transform is the identity plus $h$ times the invariant Laplacian, and the symbol product differs from the pointwise product at order $h$ by a derivative pairing.",
    audit="Unchecked: asymptotic statement not re-derived. The pairing $\\partial_{\\bar z}A_1\\,\\partial_zA_2$ is the one expected when $z$ labels the creation side.",
    cell="none. QF has no disk coherent states.")

E[60] = dict(
    desc="On the sphere the reproducing kernel is $(1+z\\bar w)^{1/h}$ with $1/h$ an integer.",
    audit="Correct: with $1/h=2j$ this is the unnormalized SU(2) coherent overlap. Checked by: comparison with the SU(2) overlap check (exact).",
    cell="none. No su(2) layer in QF.")

E[61] = dict(
    desc="A linearly forced oscillator stays coherent, the amplitude integrates the force, a phase accumulates as $\\dot\\phi=\\mathrm{Im}(\\bar\\beta\\alpha)$, and from the vacuum the number distribution is Poisson.",
    audit="Correct. Checked by: QF exact for the step $\\hat D(\\beta\\,dt)\\hat D(\\alpha)=e^{dt(\\beta\\bar\\alpha-\\bar\\beta\\alpha)/2}\\hat D(\\alpha+\\beta\\,dt)$, so $\\dot\\phi=-\\mathrm{Im}(\\beta\\bar\\alpha)=\\mathrm{Im}(\\bar\\beta\\alpha)$, and for the amplitude $\\langle n|\\hat D(\\gamma)|0\\rangle=e^{-|\\gamma|^2/2}\\gamma^n/\\sqrt{n!}$ (symbolic $n$), whose square is the Poisson law.",
    cell=True, match="Matches.")

E[62] = dict(
    desc="The quadratic oscillator operators $\\hat K_0$ and $\\hat K_2$ realize su(1,1); $\\hat K_0$ is half the shifted number operator.",
    audit="Convention: $\\hat K_0=\\tfrac12(\\hat a^\\dagger\\hat a+\\tfrac12)$ is Correct; $\\hat K_2=\\tfrac{i}{4}(\\hat a^{\\dagger2}-\\hat a^2)=-(\\hat K_+-\\hat K_-)/(2i)$ for the single-mode $\\hat K_\\pm$, and the source's $\\hat K_1=\\tfrac14(\\hat p^2/\\omega_0-\\omega_0\\hat q^2)=-(\\hat a^{\\dagger2}+\\hat a^2)/4$ flips sign too, so this basis matches the single-mode realization after $\\hat a\\to i\\hat a$; the algebra closes either way. Checked by: QF exact, with $\\hat q=(\\hat a+\\hat a^\\dagger)/\\sqrt{2\\omega_0}$, $\\hat p=-i\\sqrt{\\omega_0/2}(\\hat a-\\hat a^\\dagger)$, $\\hbar=1$.",
    cell=True, match="Matches: the outputs are $\\hat K_0$, $\\hat K_2$ and $\\hat K_1$ in ladder form.")

E[63] = dict(
    desc="Under a Hamiltonian linear in the su(1,1) generators, a coherent state stays coherent; its label obeys a Riccati equation and its phase obeys a first-order equation.",
    audit="Correct as corrected in the atlas; the printed phase sign is wrong. Checked by: exact WL computation in the holomorphic realization, $\\psi=e^{-i\\phi}(1-\\zeta\\bar\\zeta)^k(1-\\zeta z)^{-2k}$ with $\\zeta,\\bar\\zeta,\\phi$ functions of $t$: the Schrodinger residual divided by $\\psi$ is exactly 0 with the corrected phase, and $2ik(\\bar\\beta\\zeta-\\beta\\bar\\zeta)=-4k\\,\\mathrm{Im}(\\bar\\beta\\zeta)$ with the printed one.",
    cell="none. No su(1,1) layer in QF.")

E[64] = dict(
    desc="The sign of the Lorentz norm of the generator classifies a periodic quadratic oscillator: discrete spectrum, two-sided continuum, or one-sided continuum.",
    audit="Correct as a statement about Lie-algebra generators (timelike elements are conjugate to $\\Omega\\hat K_0$, whose spectrum is $\\Omega(k+n)$); the existence of a real logarithm of the period map is a separate assumption, as the atlas says. Checked by: hand derivation for the timelike case; other cases Unchecked.",
    cell="none. QF has no Floquet classification.")

E[65] = dict(
    desc="The half-line inverse-square oscillator has ground state $x^ae^{-x^2/2}$ with $a(a-1)=2g^2$, energy $a+1/2$, and an equally spaced ladder with su(1,1) label $k=E_0/2$.",
    audit="Correct. Checked by: exact WL computation, $\\hat H(x^ae^{-x^2/2})/(x^ae^{-x^2/2})$ simplifies to $a+1/2$ when the inverse-square coefficient is $a(a-1)/2$; $k=E_0/2$ since $\\hat K_0=\\hat H/2$.",
    cell="none. Not a Fock-space operator identity.")

E[66] = dict(
    desc="The degeneracies of the Calogero model after removing the center of mass are generated by $\\prod_{j=2}^N(1-t^j)^{-1}$.",
    audit="Correct by counting: symmetric polynomials in $N$ variables are generated by power sums of degree $1..N$; removing center-of-mass motion removes degree 1. Checked by: hand derivation.",
    cell="none. Not an operator identity.")

E[67] = dict(
    desc="The state $e^{\\zeta\\hat a^{\\dagger2}/2+\\alpha\\hat a^\\dagger}|0\\rangle$ is normalized by $(1-|\\zeta|^2)^{1/4}$ times a Gaussian in $\\alpha$ that includes the cross term $\\mathrm{Re}(\\bar\\zeta\\alpha^2)$.",
    audit="Correct as corrected; the printed $N=(1-|\\zeta|^2)^{-1/2}e^{-|\\alpha|^2/2}$ is wrong. Checked by: exact WL computation: the squared norm, the Bargmann Gaussian integral $\\int d^2w/\\pi\\,e^{-|w|^2}|e^{\\zeta w^2/2+\\alpha w}|^2$ with symbolic complex $\\zeta,\\alpha$, times the corrected $N^2$ is exactly 1; with the printed $N$ it is $(1-\\zeta^2)^{-3/2}e^{\\alpha^2\\zeta/(1-\\zeta)}\\ne1$ already for real $\\zeta,\\alpha$.",
    cell="none. `BosonicMatrixElement[{0, 0}, Exp[ca a] ** Exp[cz/2 a ** a] ** Exp[z/2 ad ** ad] ** Exp[al ad]]` returns unevaluated, and `BosonicExpOrder` of that product orders the quadratic factors but leaves the linear wings unordered.")

E[68] = dict(
    desc="For a forced parametric oscillator the vacuum amplitude has modulus $(1-\\rho)^{1/4}e^{-\\nu(1-\\sqrt\\rho\\cos2\\varphi)/2}$.",
    audit="Unchecked: the map from the force and frequency history to $(\\rho,\\nu,\\varphi)$ is a source lookup. At $\\nu=0$ it reduces to the squeezed-vacuum amplitude $(1-|\\zeta|^2)^{1/4}$ of the affine Gaussian entry, consistent with $\\rho=|\\zeta|^2$.",
    cell="none. The source parameters are not reconstructed.")

E[69] = dict(
    desc="A canonical map mixing particles and antiparticles preserves the indefinite form $E$, which fixes the block relations and the inverse.",
    audit="Correct. Checked by: hand derivation, blockwise expansion of $MEM^\\dagger=E$; same structure as the finite-mode bosonic Bogoliubov entry with $(\\hat a,\\hat b^\\dagger)$ in place of $(\\hat a,\\hat a^\\dagger)$.",
    cell="none. QF has no coefficient-matrix layer.")

E[70] = dict(
    desc="The multimode pair-creation exponential on the vacuum has squared norm $\\det(I_N-Z^\\dagger Z)^{-1}$, so the stated prefactor normalizes it.",
    audit="Correct (the relation $Z=(A^\\dagger)^{-1}C^\\dagger=BD^{-1}$ is convention-dependent and Unchecked). Checked by: QF exact, $N=2$ (four modes, symbolic complex $Z$): the vacuum element of $e^{\\sum\\bar Z_{ij}\\hat a_i\\hat b_j}e^{\\sum Z_{ij}\\hat a_i^\\dagger\\hat b_j^\\dagger}$ expanded in powers of $Z$ matches $\\det(I-Z^\\dagger Z)^{-1}$ through third order; hand derivation for all orders: in Bargmann variables the $v$ integral of $e^{-|v|^2+u^TZv+\\overline{u^TZv}}$ gives $e^{u^TZZ^\\dagger\\bar u}$, and the $u$ integral then gives $\\det(I-Z^\\dagger Z)^{-1}$ for $I-Z^\\dagger Z>0$.",
    cell=True, match="Matches: the differences at orders 1, 2, 3 are 0.")

E[71] = dict(
    desc="Two multimode paired vacua overlap by the product of square-root normalizations over $\\det(I_N-Z_1^\\dagger Z_2)$.",
    audit="Correct for the overlap; the group action is Unchecked (convention-dependent orientation, as the atlas notes). Checked by: QF exact, $N=2$ with independent symbolic $Z_1,Z_2$, through third order; all orders by the same Bargmann integral with $u^TZ_2v$ and $\\overline{u^TZ_1v}$, which gives $\\det(I-Z_1^\\dagger Z_2)^{-1}$.",
    cell=True, match="Matches: the differences at orders 1, 2, 3 are 0.")

E[72] = dict(
    desc="The pair-number distribution of a multimode paired vacuum has generating function $\\det(I-Z^\\dagger Z)/\\det(I-tZ^\\dagger Z)$, and its logarithm is the trace series with a plus sign.",
    audit="Correct as corrected; the printed minus sign before the trace sum is wrong. Checked by: hand derivation from the overlap entry: $t^{\\hat N_a}\\hat a_i^\\dagger t^{-\\hat N_a}=t\\hat a_i^\\dagger$ (QF exact) gives $t^{\\hat N_a}|Z\\rangle\\propto|tZ\\rangle$, so $\\mathcal G(t)=\\det(I-Z^\\dagger Z)\\langle0|\\ldots\\rangle_{Z,tZ}=\\det(I-Z^\\dagger Z)/\\det(I-tZ^\\dagger Z)$; exact WL computation for a symbolic $2\\times2$ matrix $M=Z^\\dagger Z$: the plus-sign trace series matches $\\log\\mathcal G$ through $t^4$, and the minus sign already fails at order $t$.",
    cell=True, match="Matches: $e^{\\lambda\\hat N}\\hat a^\\dagger e^{-\\lambda\\hat N}=e^\\lambda\\hat a^\\dagger$ (simplified for $\\lambda>0$; the raw output is even in $\\sqrt{\\lambda^2}$, so it holds for all $\\lambda$); plus sign residual 0; minus sign residual $2\\operatorname{Tr}M\\,t$.")

E[73] = dict(
    desc="The pair matrix of a multimode paired vacuum obeys a matrix Riccati equation, and its vacuum phase obeys a trace equation.",
    audit="Correct. Checked by: hand derivation. The pair-creation coefficient gives $[C_{ij}\\hat a_i^\\dagger\\hat a_j,\\,Z_{kl}\\hat a_k^\\dagger\\hat b_l^\\dagger]\\to(CZ)$, the $D$ term gives $ZD^T$, and the double commutator of $\\bar A_{ij}\\hat a_i\\hat b_j$ gives $ZA^\\dagger Z$. The vacuum coefficient gives $\\dot\\phi+i\\dot N/N=\\operatorname{Tr}(A^\\dagger Z)$; its real part is the stated phase, and its imaginary part reproduces $\\frac{d}{dt}\\tfrac12\\log\\det(I-Z^\\dagger Z)$.",
    cell="none. QF has no time-dependent Gaussian-state equations.")

E[74] = dict(
    desc="One fermionic particle-antiparticle pair is a spin-1/2 problem: the pair operators act as $\\hat J_\\pm,\\hat J_0$ on the vacuum and the occupied pair.",
    audit="Correct. Checked by: hand derivation; on $\\{|0\\rangle,\\hat a^\\dagger\\hat b^\\dagger|0\\rangle\\}$, $[\\hat J_+,\\hat J_-]=\\mathrm{diag}(-1,1)=2\\hat J_0$; exact fermion matrices (4 modes) confirm $[\\hat J_+,\\hat J_-]-2\\hat J_0$ annihilates the vacuum. The Dirac-field reductions are Unchecked.",
    cell="none. no fermionic layer.")

E[75] = dict(
    desc="The doubled fermionic coefficient matrix is unitary, which fixes the block relations.",
    audit="Correct. Checked by: hand derivation, blockwise expansion of $MM^\\dagger=I$ (same algebra as the exact fermionic Bogoliubov check).",
    cell="none. no fermionic layer.")

E[76] = dict(
    desc="The fermionic paired vacuum $e^{\\sum Z_{ij}\\hat a_i^\\dagger\\hat b_j^\\dagger}|0,0\\rangle$ has squared norm $\\det(I_N+Z^\\dagger Z)$.",
    audit="Correct (the chart relation to $A,B,C,D$ is Unchecked). Checked by: exact fermion matrices (4 modes), exact rational complex $2\\times2$ $Z$: squared norm minus $\\det(I+Z^\\dagger Z)$ is exactly 0.",
    cell="none. no fermionic layer.")

E[77] = dict(
    desc="Two fermionic paired vacua overlap by $\\det(I_N+Z_1^\\dagger Z_2)$ over the normalizations.",
    audit="Correct for the overlap; the group action is Unchecked (convention). Checked by: exact fermion matrices (4 modes), exact rational $Z_1,Z_2$: residual exactly 0.",
    cell="none. no fermionic layer.")

E[78] = dict(
    desc="The fermionic pair-number generating function is $\\det(I+tZ^\\dagger Z)/\\det(I+Z^\\dagger Z)$, a product of Bernoulli laws.",
    audit="Correct. Checked by: exact fermion matrices (4 modes), $t^{\\hat N_a}=\\prod_i(1+(t-1)\\hat n_{a_i})$ with symbolic $t$: residual exactly 0 as a polynomial in $t$. Bernoulli factorization follows from the eigenvalues of $Z^\\dagger Z$.",
    cell="none. no fermionic layer.")

E[79] = dict(
    desc="For lattice spin systems the quantum partition function is squeezed between classical ones at $\\beta$ and at a slightly larger $\\beta$.",
    audit="Unchecked: Simon's theorem cited from the source; the lower bound is the Berezin-Lieb lower bound applied to a product of spin coherent states.",
    cell="none. QF has no partition-function bounds.")

E[80] = dict(
    desc="The thermal damped-oscillator master equation makes $\\langle\\hat a\\rangle$ decay at rate $\\gamma/2$ and $\\langle\\hat n\\rangle$ relax to $\\nu$ at rate $\\gamma$.",
    audit="Correct (standard). Checked by: QF exact, the adjoint generator $\\mathcal L^\\dagger(X)=\\gamma(\\nu+1)(\\hat a^\\dagger X\\hat a-\\tfrac12\\{\\hat a^\\dagger\\hat a,X\\})+\\gamma\\nu(\\hat aX\\hat a^\\dagger-\\tfrac12\\{\\hat a\\hat a^\\dagger,X\\})$ applied to $\\hat a$ and $\\hat a^\\dagger\\hat a$.",
    cell=True, match="Matches: $\\mathcal L^\\dagger(\\hat a)=-\\tfrac\\gamma2\\hat a$ and $\\mathcal L^\\dagger(\\hat a^\\dagger\\hat a)=-\\gamma(\\hat a^\\dagger\\hat a-\\nu)$.")

E[81] = dict(
    desc="Under thermal attenuation each ordered characteristic function is the initial one at contracted argument times a Gaussian that grows to the noise level.",
    audit="Correct. Checked by: exact WL computation: with $\\eta,\\bar\\eta$ independent and an arbitrary initial function, the stated solution satisfies the PDE identically. For thermal input $\\chi_N=e^{-\\nu_0|\\eta|^2}$ the solution is thermal with $\\mu=p\\nu_0+q\\nu$ (hand).",
    cell="none. QF has no characteristic-function evolution.")

E[82] = dict(
    desc="The ordered quasiprobabilities obey Fokker-Planck equations with drift $-\\gamma\\alpha/2$, and their covariance relaxes linearly to the noise level plus the ordering offset.",
    audit="Correct. Checked by: hand derivation; a Gaussian factor $e^{-s|\\eta|^2}$ in $\\chi_k$ corresponds to variance $s/2$ per real quadrature, giving the additive $q(\\nu+\\sigma_k)/2$.",
    cell="none. QF has no quasiprobability evolution.")

E[83] = dict(
    desc="Attenuation acts on each ordered quasiprobability as convolution with a normalized Gaussian centered at $\\sqrt p\\,\\alpha'$.",
    audit="Correct as corrected; the printed positive exponent is not integrable. Checked by: exact WL computation: the corrected kernel integrates to 1, and its Fourier transform is $e^{-s|k|^2/4}$ with $s=q(\\nu+\\sigma_k)$, the Gaussian factor of the characteristic-function solution.",
    cell="none. QF has no quasiprobability evolution.")

E[84] = dict(
    desc="The unnormalized coherent matrix element of a density operator is a generating function whose Taylor coefficients are the Fock matrix elements.",
    audit="Correct (with $z_1\\leftrightarrow\\bar\\alpha$, $\\bar z_2\\leftrightarrow\\beta$). Checked by: hand derivation from the coherent expansion $e^{|\\alpha|^2/2}|\\alpha\\rangle=\\sum_n\\alpha^n|n\\rangle/\\sqrt{n!}$.",
    cell="none. The displaced thermal entry below exercises this kernel with QF.")

E[85] = dict(
    desc="The coherent kernel of a displaced thermal state is a single Gaussian in the two labels, with mean occupation $\\mu$ and displacement $\\beta$.",
    audit="Correct. Checked by: QF exact: $\\hat\\rho=\\hat D(\\beta)\\,e^{\\lambda\\hat N}\\hat D(\\beta)^\\dagger/(1+\\mu)$ with $e^\\lambda=\\mu/(1+\\mu)$; the displacements move the coherent labels with the Weyl phases (`BosonicExpSimplify`), and the coherent element of $e^{\\lambda\\hat N}$ is QF's closed form. The kernel with $z_1=\\bar\\alpha$, $\\bar z_2=\\alpha'$ divided by the stated $R$ is 1.",
    cell=True, match="Matches.")

E[86] = dict(
    desc="The populations of a displaced thermal state are a geometric factor times a Gaussian times a Laguerre polynomial of negative argument.",
    audit="Correct, including the Poisson limit $\\mu\\to0$. Checked by: exact WL computation: $n!$ times the $z_1^n\\bar z_2^n$ Taylor coefficient of the kernel $R$ (QF exact in the previous entry) equals the stated $W_n$ for $n\\le6$ with symbolic $\\mu,\\beta$.",
    cell="none. `BosonicMatrixElement` handles $\\hat D(\\beta)$ or a function of $\\hat N$ alone, not their product; the populations are read off the exact kernel of the previous entry.")

E[87] = dict(
    desc="A linear drive adds $-i\\epsilon$ to the decay of the mean field, and the drive enters the state as a displacement.",
    audit="Correct. Checked by: QF exact for the drive term, $-i[\\hat a,\\hat H_d]=-i\\epsilon$; the damping term $-\\tfrac\\gamma2\\langle\\hat a\\rangle$ is the master-equation entry.",
    cell=True, match="Matches.")

E[88] = dict(
    desc="The coherent matrix element of $e^{-it\\omega\\hat a^\\dagger\\hat a}$ is the overlap with a rotated label, and the thermal trace is the geometric series.",
    audit="Correct (and restoring vacuum energy multiplies by $e^{-\\omega/(2T)}$). Checked by: QF exact.",
    cell=True, match="Matches: with $\\lambda=-i\\omega t$ or $-\\omega/T$ the element is $e^{\\bar\\alpha\\beta e^\\lambda-(|\\alpha|^2+|\\beta|^2)/2}$, and the trace is $(1-e^{-\\omega/T})^{-1}$.")

E[89] = dict(
    desc="The trace of $e^{-s(2\\hat N+1)}$ is $1/(2\\sinh s)$ for $s>0$; the source's diagonalizing $\\hat U$ is a complex continuation, not a unitary.",
    audit="Correct, with the atlas's qualification that $\\hat U$ has a Hermitian generator and is a complex canonical continuation. Checked by: exact WL computation (the sum minus $1/(2\\sinh s)$ simplifies to 0); hand check that $i(\\hat b^{\\dagger2}-\\hat b^2)$ is Hermitian.",
    cell="none. A trace of a function of $\\hat N$ is an ordinary series; QF adds nothing.")

E[90] = dict(
    desc="The generating function of $L_n^k$ in the lower index is $(1-t)^{-k-1}e^{-tx/(1-t)}$.",
    audit="Correct. Checked by: exact WL computation, series coefficients through $t^6$ with symbolic $k$ and $x$ all equal.",
    cell="none. Not an operator identity.")

E[91] = dict(
    desc="Along the diagonal $m-n$ the Laguerre generating function is $(1+t)^me^{-tx}$.",
    audit="Correct. Checked by: exact WL computation, $m=0..6$ through order $t^{10}$, all differences 0.",
    cell="none. Not an operator identity.")

E[92] = dict(
    desc="The double generating function of the displacement transition probabilities is $(1-uv)^{-1}\\exp[-|\\gamma|^2(1-u)(1-v)/(1-uv)]$.",
    audit="Correct as corrected; the printed extra $e^{|\\gamma|^2}$ is wrong. Checked by: QF exact: the exact matrix elements $D_{mn}(\\gamma)$ with $\\gamma=g_x+ig_y$ symbolic give $\\sum_{m,n\\le5}u^mv^n|D_{mn}|^2$, which equals the double Taylor polynomial of the corrected form through $u^5v^5$; at $u=v=0$ the sum is $e^{-|\\gamma|^2}$, while the printed form gives 1.",
    cell=True, match="Matches the corrected form: difference 0, and $F(0,0)=e^{-|\\gamma|^2}$.")

E[93] = dict(
    desc="Each displacement matrix element is a mixed Wirtinger derivative of $e^{-|\\alpha|^2}$, rescaled; transitions are classically allowed where the two number circles intersect.",
    audit="Correct for the exact formula; the semiclassical region is Unchecked (asymptotic). Checked by: QF exact, symbolic with $\\alpha,\\bar\\alpha$ independent, against the normal-ordered matrix elements for $m,n\\le5$: all differences 0.",
    cell=True, match="Matches.")

E[94] = dict(
    desc="The complex Jacobian of the fractional-linear action on a type-I domain is $\\det(B^Tz+D)^{-(p+q)}$.",
    audit="Unchecked: depends on the source's transposed-block convention for type-I domains; not re-derived.",
    cell="none. Not an operator identity.")

E[95] = dict(
    desc="A discrete SU(1,1) coherent subsystem is complete when the fundamental area is below $\\pi/(2k-1)$ and incomplete for large $k$.",
    audit="Unchecked: theorem cited from the source (automorphic-form argument), area convention source-specific.",
    cell="none. QF has no su(1,1) layer.")
