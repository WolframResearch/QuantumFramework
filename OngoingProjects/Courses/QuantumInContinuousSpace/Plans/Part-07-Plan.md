# Part 7 Plan: Approximation methods

Ten questions. Class census per the class census table in `Route-Table.md`: C0 (no differential
equation) 7.1, 7.2, 7.3, 7.4, 7.10; C9 (WKB / semiclassical asymptotics) 7.5, 7.6, 7.9; C5
(truncated-basis ODE-IVP systems) 7.7, 7.8.

## Common ground

Every question in this part replaces an intractable $\hat H = \hat H_0 + V$ by a controlled
surrogate, and each method carries its own error handle. Rayleigh-Schrodinger perturbation theory:
$E_n \approx E_n^{(0)} + \langle n|V|n\rangle + \sum_{m \neq n} \frac{|\langle m|V|n\rangle|^2}{E_n^{(0)} - E_m^{(0)}}$,
with degeneracy forced through the secular problem $\det(V_{ab} - E^{(1)}\delta_{ab}) = 0$. The
variational principle bounds from above, $E_0 \le \langle\psi_T|\hat H|\psi_T\rangle / \langle\psi_T|\psi_T\rangle$;
its finite-basis form is the generalized eigenproblem $Hc = ESc$ with overlap matrix
$S_{ij} = \langle\phi_i|\phi_j\rangle$; Temple's inequality bounds from below,
$E_0 \ge \langle H\rangle - \frac{\langle H^2\rangle - \langle H\rangle^2}{E_1 - \langle H\rangle}$
for $\langle H\rangle < E_1$. Semiclassics quantizes by the action,
$\oint p\,dx = 2\pi(n + \tfrac12)$ with $p = \sqrt{2(E - V)}$, tunnels with the Gamow factor
$T \sim e^{-2\int |p|\,dx}$ across the forbidden region, and connects across a turning point
through the Airy asymptotics with the $\pi/4$ phase. Time-dependent perturbation theory yields
$\Gamma = 2\pi\,|\langle f|V|i\rangle|^2 \rho(E_f)$ (Fermi's golden rule), and a parameter ramp of
duration $T$ is bracketed by the sudden rule $P(n) = |\langle n^{\mathrm{new}}|\psi\rangle|^2$ as
$T \to 0$ and by instantaneous-eigenstate tracking as $T \to \infty$.

### 7.1 [BSc] How do I compute first- and second-order energy shifts by nondegenerate time-independent perturbation theory (the anharmonic oscillator)?

Perturb the oscillator with $\lambda x^4$ and compute $E_n^{(1)}$ and $E_n^{(2)}$ for general $n$:
get the $x^4$ matrix elements by Integrate over HermiteH eigenfunctions and again by ladder
algebra (the two must agree, and only $m \in \{n, n \pm 2, n \pm 4\}$ survive), so first order is
the known $\tfrac{3\lambda}{4}(2n^2 + 2n + 1)$ and second order a four-term Sum in closed form.
Benchmark against the truncated Fock-basis Hamiltonian (SparseArray diagonal plus banded $x^4$
block, Eigenvalues, $N$-swept): the residual against the numeric level must fall as $\lambda^3$
under $\lambda$-halving, a scaling a wrong second-order sum cannot fake. Close by watching the
orders grow with $n$ and $\lambda$: the series is asymptotic, and the coupling where it stops
helping is located, not asserted.

### 7.2 [MSc] How do I resolve a degeneracy by degenerate perturbation theory?

Split the ring's degenerate $m = \pm 1$ doublet with $\lambda\cos 2\phi$: the naive formula fails
twice (the diagonal element vanishes, the degenerate denominator blows up), so build the
$2 \times 2$ secular matrix from Integrate over $\phi \in (0, 2\pi)$ and Eigensystem it, giving
$E^{(1)} = \pm\lambda/2$ with correct zeroth-order states $\cos\phi$ and $\sin\phi$. The exact
anchor is that this ring is the Mathieu problem under $a = 2E$, $q = \lambda$:
MathieuCharacteristicA and MathieuCharacteristicB (probed live in the Route-Table C1 verdict) give
the exact split levels, and their Series in $q$ must recover $2E^{(1)} = \pm q$, a check that can
refute the secular matrix rather than re-derive it; the finite-$q$ deviation of order $\lambda^2$
shows where first order ends. Close physically: the split pair are standing waves pinned to the
minima and maxima of the potential, the onset of hindered rotation.

### 7.3 [BSc] How do I get a Rayleigh-Ritz variational upper bound on the ground-state energy from a Gaussian trial function?

Keep the Gaussian trial the question names, $e^{-bx^2/2}$ on $V = x^4/4$: Integrate for
$\langle H\rangle(b)$ under $b > 0$ assumptions, then take stationarity by Solve on
$\partial_b \langle H\rangle = 0$ built with D, confirming the second derivative positive
(Minimize returns a Piecewise on parameter signs, the known PIPELINE gotcha). Contrast with the
two-parameter trial $(1 + cx^2)\,e^{-bx^2/2}$, which strictly lowers the bound and quantifies the
Gaussian's deficit. Anchor the machinery on the harmonic oscillator, where the same trial recovers
$E_0 = \tfrac12$ exactly, and check the scaling law $E_0(\alpha x^4) = \alpha^{1/3} E_0(x^4)$ on
the functional. The refuting check is the ordering $E_{\mathrm{2par}} \le E_{\mathrm{Gauss}}$ with
both above the FD numeric $E_0$ from the C2-certified recipe (SparseArray Band Hamiltonian,
Eigenvalues with Arnoldi Shift, Richardson); an inversion refutes the integrals. Close on the
tail: the true state decays as $e^{-x^3/(3\sqrt2)}$, faster than any Gaussian, and the
two-parameter trial buys back part of exactly that mismatch.

### 7.4 [MSc] How do I run the linear variational (Ritz) method in a finite basis, reducing the problem to a generalized matrix eigenproblem?

Take the shallow-barrier quartic double well $V = (x^2 - a^2)^2/2$ at $a = 3/2$ (shallow by
design: 19.6 owns the deep-barrier instanton regime of this family) in the nonorthogonal basis of
two displaced Gaussians $e^{-b(x \mp a)^2/2}$ plus polynomial corrections
$x^k e^{-b(x \mp a)^2/2}$: every $H_{ij}$ and $S_{ij}$ closes as a Gaussian-polynomial Integrate,
and the generalized problem $Hc = ESc$ goes to Eigensystem[{hmat, smat}] (two-matrix calling
form, verify at authoring). Two bare Gaussians alone give the textbook $2 \times 2$; the
enlargement makes the overlap matrix do real work and shows Ritz convergence from above: every
level must descend as the basis grows (min-max), and one rising level refutes. Benchmark the
tunneling splitting against the C2-certified FD Richardson recipe (that verdict's two routes agree
to $9 \times 10^{-15}$ on double-well splittings at $a = 2$; the recipe transfers, the numbers
are re-measured at $a = 3/2$). Close at the edge where the enlarged basis approaches linear
dependence and $S$ goes near-singular: the conditioning limit of nonorthogonal Ritz.

### 7.5 [BSc] How do I quantize a smooth well by the WKB (Bohr-Sommerfeld) condition?

Quantize the quartic well $V = x^4/4$ by $\oint p\,dx = 2\pi(n + \tfrac12)$: raw Integrate with
turning-point limits stalls (C9 trap), so scale the turning point out with $x = (4E)^{1/4} t$; the
action closes to $S(E) = \tfrac{16}{3} K(-1) E^{3/4}$ with $K$ from EllipticK, and Solve on
$S(\mathrm{en}) = 2\pi(n + \tfrac12)$ (en, never E, which is Euler's number) gives the closed
$E_n^{\mathrm{WKB}} \propto (n + \tfrac12)^{4/3}$ (C9 verdict, probe p1). Benchmark against the
FD SparseArray plus Eigenvalues plus Richardson spectrum per the C2 recipe and measure the
verdict's error law: relative error 18.2 percent at $n = 0$ down to 0.083 percent at $n = 6$ with
$\mathrm{relerr} \cdot (n + \tfrac12)^2 \to -0.0352$; a wrong action or Maslov term destroys that
$1/n^2$ collapse, a sharper refuter than level-by-level closeness. Never trust Newton-Leibniz
antiderivative limits across a turning point (silently 0, no message), and guard any FindRoot
objective containing NIntegrate with ?NumericQ. Close on the Morse potential as the anchor where
WKB with the half-integer term is exact (its action needs its own substitution before Integrate
closes it, verify at authoring): exactness at one end of the part, a measured $1/n^2$ law at the
other.

### 7.6 [MSc] How do I compute a WKB tunneling rate through a barrier (the Gamow factor)?

Compute the Gamow factor on the inverted parabolic barrier $V = V_0 - \tfrac12 x^2$, the barrier
with an exact transmission benchmark: the tunneling exponent closes under Integrate with the
turning points scaled out, exactly $2\pi(V_0 - E)$ in the scaled units $\varepsilon = V_0 - E$
(C9 verdict, probe p3), while DSolve returns the exact ParabolicCylinderD solutions behind the
Kemble formula $T = 1/(1 + e^{2\pi\varepsilon})$. Verify the deep-tunneling limit,
$\log T + 2\pi\varepsilon \to 0$ under Limit with next order $-1$, and cross-check with the
independent NDSolveValue scattering integration, which matches Kemble to $5 \times 10^{-8}$ at
$\varepsilon \in \{2, 1, 0.5, 0, -1\}$. ParabolicCylinderD along rotated rays is an open risk
flagged in the verdict, so the connection-formula step is probed at authoring; the Kemble closed
form itself is the certified anchor. Close over the top: at $\varepsilon \le 0$ the bare Gamow
factor is meaningless while Kemble stays valid, with $T(0) = \tfrac12$ exactly at the crest.

### 7.7 [MSc] How do I compute transition rates by time-dependent perturbation theory and recover Fermi's golden rule with a continuum density of states?

Photodetach the delta-well bound state $\psi_b = \sqrt{g}\, e^{-g|x|}$ with $E_b = -g^2/2$ using
an oscillating force $F x \cos\omega t$: a genuine discrete-to-continuum transition whose pieces
all close. The golden-rule side runs entirely through Integrate (delta matrix elements never
through NIntegrate, which silently returns 0 on a point measure): the bound-free dipole element on
the exact delta-well continuum states and the free 1D density of states
$\rho(E) \propto 1/\sqrt{E}$ (normalization convention fixed at authoring) give the closed rate
$\Gamma(\omega, g, F)$ above threshold $\omega > g^2/2$. Cross-check with first-order amplitude
ODEs on a box-discretized continuum via NDSolveValue at PrecisionGoal and AccuracyGoal 10 or
higher (default-tolerance norm drift grows with integration time, C5 trap), and identify the
linear-growth window honestly: fit the total transition probability's slope only past the
bandwidth transient and before depletion, showing both failure regimes outside it. Route
agreement is blind to truncation (C5) and this generator is unprobed at scale, so box-doubling
and level-density invariance of the fitted slope carry the refutation. Close at threshold, where
$\Gamma$ switches on with the $1/\sqrt{E}$ continuum edge.

### 7.8 [MSc] How do I apply the sudden and the adiabatic approximations to a continuous system?

Double the infinite well, $L \to 2L$, suddenly and slowly on the same system so the crossover is
visible. The sudden side is closed-form Integrate: overlaps
$P(n) = |\langle \chi_n^{(2L)} | \psi_1^{(L)} \rangle|^2$, checked by the exact sum rules
$\sum_n P(n) = 1$ and $\sum_n P(n)\, E_n^{(2L)} = E_1^{(L)}$ (energy is conserved through the
quench), both refutable. The ramped side expands on the instantaneous basis, couplings
$\langle \chi_n | \partial_t \chi_m \rangle$ closed by Integrate, and propagates the truncated
amplitude system with NDSolveValue (the C5 primary for a nonautonomous generator) at
PrecisionGoal and AccuracyGoal 10 or higher, norm monitored (default-tolerance drift grows with
integration time) and truncation $N$-swept (route agreement alone is blind to it). Sweep the ramp
time $T$ from sudden to adiabatic and read the $1/T^2$ approach; when the residual transition
probability is small, benchmark by projecting on the instantaneous adiabatic basis or against an
exact finite-$T$ solution, never the bare adiabatic (Landau-Zener-type) asymptote: at finite $T$
the boundary correction of order $1/T^2$ dominates and a converged solver "fails" the asymptote
by orders of magnitude (C5 trap, probe p7). Close where the sweep lands on the sudden overlaps at
one end and on adiabatic tracking at the other.

### 7.9 [MSc] How do I match WKB solutions across a turning point with the Airy connection formulas?

Derive the Airy connection at the linear turning point $V = Fx$: the $\pi/4$ phase is read
directly off Series[AiryAi[x], {x, -Infinity, 1}], which returns
$\cos(\tfrac{\pi}{4} - \tfrac{2}{3}(-x)^{3/2}) / (\sqrt{\pi}\,(-x)^{1/4})$ (C9 verdict, probe
p4), so the connection formula is earned from asymptotics the kernel exhibits, not quoted. Charge
it against the quantum bouncer (hard wall at $x = 0$ plus linear potential): the eigenstates
$\mathrm{Ai}(2^{1/3}(x - E_n))$ have residual exactly 0 and the exact spectrum sits at the
AiryAiZero points, while WKB with the wall-plus-turning-point Maslov count $(n + \tfrac34)$ lands
within $7.6 \times 10^{-3}$ to $1.0 \times 10^{-4}$ over $n = 0..6$. The discriminating contrasts
do the refuting: the wrong Maslov $\tfrac12$ is off 24 percent at $n = 0$, and a misread $\pi/2$
phase leaves an order-one mismatch in the asymptotic cosine (error 0.89 against
$4.2 \times 10^{-4}$ at $z = 15$). Close with the phase as physics: a quarter of a wave lost at
the soft turning point and half at the hard wall, summing to the $(n + \tfrac34)$ count.

### 7.10 [MSc] How do I get a variational *lower* bound on the ground-state energy (Temple's inequality)?

Cage the quartic ground energy from both sides: the upper bound from the 7.3-style trial
redefined locally (no cross-answer symbol leak), the lower from Temple's inequality
$E_0 \ge \langle H\rangle - \frac{\langle H^2\rangle - \langle H\rangle^2}{E_1 - \langle H\rangle}$,
valid only while $\langle H\rangle < E_1$. Both moments close under Integrate ($H^2$ acts as a
symbolic fourth derivative via D), and the gap input $E_1$ plus the reference $E_0$ come from the
C2-certified FD recipe (SparseArray Band Hamiltonian, Eigenvalues with Arnoldi Shift, Richardson
$(4v_{h/2} - v_h)/3$); only eigenvalues are used, so the C2 verdict's unprobed eigenfunction
normalization is not load-bearing here. The refuting check is the sandwich
$E_{\mathrm{Temple}} \le E_0^{\mathrm{FD}} \le E_{\mathrm{var}}$, all three built from the same
$V$; anchor on the harmonic oscillator, where the exact-width Gaussian has zero energy variance
and both bounds collapse onto $\tfrac12$ exactly. Close at the edge: de-tune the trial until
$\langle H\rangle \to E_1$ and watch the lower bound degrade and then lose validity, the price of
rigor being a certified gap.
