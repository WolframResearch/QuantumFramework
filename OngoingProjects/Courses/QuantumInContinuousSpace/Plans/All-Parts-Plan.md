# All Parts Plan: Quantum in Continuous Space, Parts 4 through 23

All 130 questions of Parts 4 through 23, each with the approach planned for its worked solution.
This file is a straight concatenation of the per-part plans (`Part-NN-Plan.md`), which remain the
editable sources; regenerate it rather than editing it here. The governing standard is
`../PIPELINE.md`, the plan contract is `Shared-Brief.md`, the pinned examples are in
`Example-Ledger.md`, and every solver route cites a kernel-probed verdict in `Route-Table.md`.

## Contents

- Part 4. The harmonic oscillator
- Part 5. The time-dependent Schrodinger equation as a PDE
- Part 6. Scattering in one dimension
- Part 7. Approximation methods
- Part 8. General theorems and structural methods
- Part 9. Periodic potentials and band structure
- Part 10. Two and three dimensions: separation of variables
- Part 11. Orbital angular momentum in continuous space
- Part 12. Central potentials and the hydrogen atom
- Part 13. Electromagnetic coupling
- Part 14. Spin coupled to spatial motion
- Part 15. Identical particles in continuous space
- Part 16. Density operators, mixed states, and the Wigner function
- Part 17. Continuous-variable quantum optics and information
- Part 18. Open quantum systems in continuous space
- Part 19. Path integrals
- Part 20. Three-dimensional scattering theory
- Part 21. Nonlinear and mean-field wave mechanics
- Part 22. Relativistic wave equations
- Part 23. From one particle to fields: the second-quantization bridge

## Part 4 Plan: The harmonic oscillator

8 questions. Class census, per the `Route-Table.md` class census: C1 = {4.1, 4.5}, both entries
inheriting the C1 verdict's binding traps in their prose; C0 = {4.2, 4.3, 4.4, 4.6, 4.7, 4.8}, WL
machinery named directly per the C0 row.

### Common ground

The part rests on: the oscillator Hamiltonian $\hat H=\tfrac12(\hat p^2+\hat x^2)$ with
$\hat p=-i\,d/dx$ (natural units $\hbar=m=\omega=1$ throughout; only 4.5 restores units); the
ladder pair $a=(\hat x+i\hat p)/\sqrt2$, $a^\dagger=(\hat x-i\hat p)/\sqrt2$ with
$[a,a^\dagger]=1$, the number operator $\hat n=a^\dagger a$, and $\hat H=\hat n+\tfrac12$; the
ladder action $a\vert n\rangle=\sqrt n\,\vert n-1\rangle$,
$a^\dagger\vert n\rangle=\sqrt{n+1}\,\vert n+1\rangle$; the normalized eigenfunctions
$\psi_n(x)=H_n(x)\,e^{-x^2/2}/\sqrt{2^n\,n!\,\sqrt\pi}$ with $E_n=n+\tfrac12$; and the vacuum
quadrature variances $\Delta x^2=\Delta p^2=\tfrac12$, the reference disk against which the
coherent and squeezed states of 4.6 through 4.8 are measured.

#### 4.1 [BSc] How do I solve the oscillator's stationary equation and obtain the Hermite-function eigenstates with energies $E_n=n+\tfrac12$?

Solve the stationary ODE $-\tfrac12\psi''+\tfrac12x^2\psi=E\psi$ with DSolve at symbolic $E$ (C1
per `Route-Table.md`), getting the general solution in ParabolicCylinderD/HermiteH form; never
attach decay conditions at $\pm\infty$ (DSolve hangs), instead select the decaying branch by hand
and read the quantization off as the manual termination condition, Solve giving $E_n=n+\tfrac12$.
Take the $n=3$ state $\psi_3\propto(2x^3-3x)\,e^{-x^2/2}$ (three nodes, at $0$ and
$\pm\sqrt{3/2}$; vibrational spectra) off the general solution, close the residual
$\hat H\psi_3-\tfrac72\psi_3\to0$ with FullSimplify under integer and positivity assumptions, and
show the nonterminating branch at generic $E$ grows like $e^{+x^2/2}$ (Series at large $x$): the
check that refutes "any $E$ is allowed" and is itself the physics of quantization. Cross-check
with independent machinery twice, per the C1 verdict: the exact DEigensystem domain form
$\{x,-\infty,\infty\}$ returning $\{\tfrac12,\tfrac32,\tfrac52,\tfrac72\}$ with Hermite-Gaussians,
and shifted NDEigenvalues with a mesh sweep (shift below $\tfrac12$ per class discipline, benign
here since the spectrum is positive). Close by reading the node count off the carrier: level $n$
carries exactly $n$ nodes.

#### 4.2 [BSc] How do I build the ladder operators $a=(\hat x+i\hat p)/\sqrt2$ and $a^\dagger$ as differential operators and verify $[a,a^\dagger]=1$?

Define the ladder pair as differential operators, $af=(xf+f')/\sqrt2$ and
$a^\dagger f=(xf-f')/\sqrt2$, as pure functions built on D, and close $[a,a^\dagger]f=f$ with
Simplify on a generic unassigned $f[x]$ (PIPELINE section 7): the generic-function proof is the
check that refutes a sign error in $\hat p=-i\,d/dx$. Then act on the nodal $n=2$ carrier
$\psi_2=(2x^2-1)\,e^{-x^2/2}/\sqrt{2\sqrt\pi}$ (ladder algebra) and earn the factors
$a\psi_2=\sqrt2\,\psi_1$ and $a^\dagger\psi_2=\sqrt3\,\psi_3$ against the normalized HermiteH
states with Integrate and FullSimplify, spot-checking $\langle\psi_1,a\psi_2\rangle=\sqrt2$ with
NIntegrate reusing the same definitions. Keep this a generic-function proof rather than a matrix
identity, since 4.4 owns the truncated matrices. Close with the general law the concrete action
instantiates: $\sqrt n$ down, $\sqrt{n+1}$ up, an identity no finite matrix representation can
hold exactly, the tension 4.4 exhibits.

#### 4.3 [BSc] How do I generate the whole spectrum algebraically from $a|0\rangle=0$ and the number operator $\hat n=a^\dagger a$?

Solve $a\psi_0=0$ as a first-order ODE with DSolve, normalize, and build the whole ladder up to
$n=4$ by $\psi_{n+1}=a^\dagger\psi_n/\sqrt{n+1}$ with FoldList (algebraic spectrum): the whole
family, not one level. Verify with checks that share no machinery with the first-order
construction: the second-order residual $\hat H\psi_4-\tfrac92\psi_4=0$ under Simplify,
$\hat n\,\psi_n=n\,\psi_n$ across the family with Table, the laddered $\psi_4$ equal to the closed
HermiteH form under FullSimplify, and the orthonormality Gram matrix by Integrate with a numeric
spot check. Descend as well: four applications of $a$ return $\sqrt{4!}\,\psi_0$ and a fifth
annihilates, the ladder floor. Close with the algebraic reading: $E_n=n+\tfrac12$ follows from
$[a,a^\dagger]=1$ and the floor alone, the Hermite functions falling out of the algebra rather
than being put in.

#### 4.4 [BSc] How do I compute the matrix elements of $\hat x$ and $\hat p$ in the number basis, giving the truncated banded-matrix representation of the oscillator (with the truncation made explicit)?

Build $a$, $a^\dagger$ at truncation $N=8$ (levels $0$ through $8$, $9\times9$) with SparseArray
and Band, form $\hat x,\hat p$ by matrix algebra with entries kept as exact surds,
$\langle m|\hat x|n\rangle=(\sqrt n\,\delta_{m,n-1}+\sqrt{n+1}\,\delta_{m,n+1})/\sqrt2$, and
exhibit the banded structure (truncation diagnostics); cross-check the entries against the
independent quadrature route, Integrate of $\psi_m\,x\,\psi_n$ over HermiteH states. The honest
truncation exhibit is the commutator corner defect
$[\hat x,\hat p]=i\,(\mathbb 1-(N{+}1)\,\vert N\rangle\langle N\vert)$ with corner entry $-iN$:
the trace of any finite commutator vanishes identically, so the corner must carry exactly $-N$, a
refuting check on the whole construction; sweep $N$ with Table to show the pattern. Then
Eigenvalues of the truncated $\tfrac12(\hat x^2+\hat p^2)$: every level below the cutoff is
exactly $n+\tfrac12$ while the top one is $N/2$ instead of $N+\tfrac12$. Close on the theorem the
defect encodes: the canonical commutation relation has no finite-dimensional representation, so a
corner of this exact size is forced, not a bug.

#### 4.5 [BSc] How do I show the ground state is a Gaussian, and restore $\omega$ and $m$ to read off the oscillator length $\sqrt{\hbar/m\omega}$?

Keep the Gaussian ground state as the carrier (the question names it, so the named-trivial policy
holds): DSolve the dimensionful ODE with symbolic $E$ and symbolic $\hbar,m,\omega$ (C1: no decay
BCs, hang trap; quantization is the manual termination read-off at $n=0$), giving
$\psi_0(x)=(m\omega/\pi\hbar)^{1/4}e^{-m\omega x^2/2\hbar}$ with $E_0=\hbar\omega/2$ and
oscillator length $\ell=\sqrt{\hbar/m\omega}$; if the three-parameter symbolic solve stalls,
nondimensionalize by $x=\ell\xi$ first and restore units by substitution, since the rescaling is
the physics being asked. Close the dimensionful residual and unit norm with FullSimplify under
$\hbar,m,\omega>0$, anchor the rescaled equation on the exact DEigensystem domain form, and refute
by scaling: recompute $\ell$ from $\Delta x^2=\ell^2/2$ by Integrate and demand
$\sqrt{\hbar/m\omega}$; Limit shows $\ell\to\infty$ as $\omega\to0$, the bound state flattening
away. The discriminating contrast is two real laboratory scales via Quantity and UnitConvert: the
CO stretch ($\omega/2\pi c\approx2143\ \mathrm{cm}^{-1}$, reduced mass $\approx6.86\,\mathrm u$)
gives $\ell\approx5\ \mathrm{pm}$, while a Rb atom in a $2\pi\times100\ \mathrm{kHz}$ optical
tweezer gives $\ell\approx34\ \mathrm{nm}$, nearly four orders of magnitude apart. Close there:
one formula spans molecular vibration to cold atoms.

#### 4.6 [MSc] How do I build a coherent state $|\alpha\rangle$ as the eigenstate of $a$, and equivalently by applying $e^{\alpha a^\dagger-\alpha^* a}$ to the ground state, and show it is a minimum-uncertainty state?

Build the coherent state at complex $\alpha=2e^{i\pi/3}$, so
$x_0=\sqrt2\,\operatorname{Re}\alpha=\sqrt2$ and $p_0=\sqrt2\,\operatorname{Im}\alpha=\sqrt6$
(quantum optics), twice. Eigenstate route: DSolve the first-order ODE
$\psi'=(\sqrt2\,\alpha-x)\psi$ and normalize with Integrate plus ComplexExpand. Displaced route:
split $D(\alpha)=e^{\alpha a^\dagger-\alpha^*a}$ by BCH as
$e^{-ix_0p_0/2}\,e^{ip_0\hat x}\,e^{-ix_0\hat p}$, legitimate because $[a,a^\dagger]$ is central
(the MSc ordering point), with $e^{-ix_0\hat p}$ acting as translation on $\psi_0$. The refuting
check is exact equality first: the ratio of the two constructions must FullSimplify to a constant
of modulus 1, since a modulus-only overlap would silently pass a missing $e^{-ix_0p_0/2}$ phase,
and both must satisfy $a\psi=\alpha\psi$; NIntegrate spot-checks the overlap. Then earn minimum
uncertainty with $\alpha$ kept symbolic complex: $\Delta x^2=\Delta p^2=\tfrac12$, so
$\Delta x\,\Delta p=\tfrac12$ exactly on the Heisenberg bound, with $\alpha\to0$ recovering the
vacuum. Close with the classical reading: $\langle\hat n\rangle=\vert\alpha\vert^2$ with relative
spread $\vert\alpha\vert^{-1}$, the packet turning classical at large $\vert\alpha\vert$.

#### 4.7 [MSc] How do I build a squeezed state by applying $e^{\frac12(\xi^* a^2-\xi a^{\dagger2})}$ to the ground state and compute its unequal quadrature variances?

Apply $S(\xi)=e^{\frac12(\xi^*a^2-\xi a^{\dagger2})}$ to the vacuum with $\xi=re^{i\pi/3}$ and $r$
symbolic (quantum optics): the phase of $\xi$ kept, not $\theta=0$, so the noise ellipse is
genuinely rotated off the $x,p$ axes. Symbolic side: apply the Bogoliubov transform
$S^\dagger aS=a\cosh r-a^\dagger e^{i\theta}\sinh r$ as replacement rules, reduce the vacuum
moments, and close
$\operatorname{Var}(X_\phi)=\tfrac12\left(\cosh2r-\sinh2r\cos(2\phi-\theta)\right)$ with
FullSimplify: principal axes at $\phi=\theta/2$ carrying $e^{-2r}/2$ and $e^{+2r}/2$, product
exactly $\tfrac14$ there and larger on any other quadrature pair. Numeric side, on independent
machinery: truncated ladder matrices via SparseArray, MatrixExp of the truncated generator applied
to the vacuum vector; truncation and exponentiation do not commute, so sweep the truncation $N$
and compare against the exact even-$n$ amplitudes
$c_{2n}=(-e^{i\theta}\tanh r)^n\sqrt{(2n)!}/(2^n\,n!\,\sqrt{\cosh r})$, a check that refutes
truncation error and a sign-convention error at once; the $(\tanh r)^n$ tail slows convergence at
large $r$. Close with the limits: $r\to0$ recovers the vacuum disk, and large $r$ pushes one
rotated quadrature below shot noise while the conjugate diverges at fixed product.

#### 4.8 [MSc] How do I compute the photon-number statistics of a number state and a coherent state (the Poisson distribution)?

Set the number state $n=4$ beside the coherent state $\alpha=2e^{i\pi/3}$, so
$\vert\alpha\vert^2=4$: identical mean energy $\langle\hat H\rangle=\tfrac92$, so only
the statistics discriminate (photon counting). Compute the amplitudes
$c_n=e^{-\vert\alpha\vert^2/2}\alpha^n/\sqrt{n!}$, probabilities with Abs under ComplexExpand,
and close the symbolic Sums $\sum_n P(n)=1$, $\bar n=\vert\alpha\vert^2$, and
$\Delta n^2=\vert\alpha\vert^2$: the check that refutes normalization or weight errors. Compare
against PDF of PoissonDistribution with Mean and Variance as the independent route, and watch the
partial-sum tail numerically. Mandel $Q=(\Delta n^2-\bar n)/\bar n$ lands at $0$ for the coherent
state and at $-1$ for the Fock state, the sub-Poissonian floor, and the large-$\bar n$ Gaussian
envelope of the Poisson is the classical limit. Close on the phase: the phase of $\alpha$ cancels
in every $P(n)$, so photon counting is blind to exactly what the quadratures of 4.6 and 4.7
resolve.

## Part 5 Plan: The time-dependent Schrodinger equation as a PDE

10 questions. Class census per the class-census table in `Route-Table.md`: C3 (time-dependent
linear Schrodinger PDE): 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.9, 5.10; C0 (no differential equation):
5.1, 5.8.

### Common ground

Everything in this part rests on the time-dependent Schrodinger equation
$i\,\partial_t\psi = -\tfrac12\,\partial_{xx}\psi + V(x)\,\psi$ and four of its consequences: a
stationary state evolves only by the phase $e^{-iEt}$; a general state propagates in the energy
eigenbasis as $\psi(x,t)=\sum_n c_n e^{-iE_n t}\psi_n(x)$ with $c_n=\int\bar\psi_n\psi_0\,dx$;
the same evolution is carried by the propagator, $\psi(x,t)=\int K(x,t;x')\,\psi(x',0)\,dx'$;
the norm $\int|\psi|^2\,dx$ is conserved through the continuity equation
$\partial_t\rho+\partial_x j=0$ with $j=\operatorname{Im}(\bar\psi\,\partial_x\psi)$; and
expectation values obey Ehrenfest's relations $\frac{d\langle x\rangle}{dt}=\langle p\rangle$ and
$\frac{d\langle p\rangle}{dt}=-\langle V'(x)\rangle$, which reduce to classical mechanics only
when $\langle V'(x)\rangle=V'(\langle x\rangle)$.

The solver contract below is authoring guidance for this plan, not material for the Part opening.
Seven of the eight C3 entries, all but 5.5, share the certified contract of the C3 verdict in
`Route-Table.md`: primary `NDSolveValue` with `Method -> {"MethodOfLines",
"SpatialDiscretization" -> {"TensorProductGrid", "MinPoints" -> n, "MaxPoints" -> n,
"DifferenceOrder" -> "Pseudospectral"}}`, periodic identification $\psi(t,L_1)=\psi(t,L_2)$, and
`AccuracyGoal -> 10, PrecisionGoal -> 10`. Pinning `MinPoints` and `MaxPoints` to the same $n$ is
what makes $n$ the swept knob rather than a solver choice, and every resolution sweep below moves
that one number over the verdict's measured ladder $n=65,129,257,513$; the tight goals are
load-bearing, since defaults are silently wrong at the $0.25\%$ level with no message and
refinement at default goals is non-monotone. Cross-check independently with a hand-built Strang
split-step Fourier propagator on `Fourier`/`InverseFourier` (`FourierParameters -> {1,-1}`,
packed arrays). Size the box from the state's own energy content and hold it fixed across every
sweep, because norm conservation does not detect wall reflection. Take moments and probabilities
as $dx$-weighted grid sums, never `NIntegrate` on an interpolant. Entry 5.5 leaves this contract
and states its own discretization.

#### 5.1 [BSc] How do I show that a stationary state evolves only by the phase $e^{-iEt}$, while a superposition acquires genuine time dependence?

Work in the Poschl-Teller family $V=-\tfrac{\lambda(\lambda+1)}{2}\operatorname{sech}^2x$, whose
bound levels are $E_n=-\tfrac{(\lambda-n)^2}{2}$, and specialize to $\lambda=2$, that is
$V=-3\operatorname{sech}^2x$ with exactly two states, $\psi_0\propto\operatorname{sech}^2x$ at
$E_0=-2$ and $\psi_1\propto\operatorname{sech}x\tanh x$ at $E_1=-\tfrac12$: a finite spectrum
makes the beat a single closed-form frequency with no truncation question. This is C0 machinery
per the `Route-Table.md` C0 row: earn both levels with `D` + `FullSimplify` residuals
$\hat H\psi_n-E_n\psi_n\to0$, then reduce the densities of $e^{-iE_0t}\psi_0$ and of the
equal-weight superposition with `ComplexExpand` under `Element[{x, t}, Reals]`. The stationary
density is $t$-free identically; the superposition carries a cross term $\cos(\Delta E\,t)$ with
$\Delta E=\lambda-\tfrac12$, so the whole family beats, and $\lambda=2$ gives $\Delta E=\tfrac32$
with period $4\pi/3$, refutable by a numeric spot value of the density at one $(x,t)$. Close by
taking `Limit` of that same cross term as $\Delta E\to0$: the density goes static, so genuine
time dependence lives in energy differences, not in superposition itself.

#### 5.2 [BSc] How do I integrate the time-dependent Schrodinger equation directly as a PDE with `NDSolve` for a given initial $\psi(x,0)$?

Release the cusp packet $\psi_0\propto e^{-|x-x_0|}$, $x_0=3$, in the harmonic trap $V=x^2/2$:
its momentum amplitude falls only as $(1+k^2)^{-1}$, a genuine grid stress a smooth default
packet would never apply. Integrate under the shared C3 contract, with a default-goals twin run
alongside to exhibit the silent-error trap live. Size the box from the state's own energy, edge
$\gtrsim\sqrt{2(\langle H\rangle+4\sigma_H)}$ plus several cusp decay lengths (the tail falls as
$e^{-|x-x_0|}$, so a handful of units suffices), and hold that box fixed while $n$ doubles, so
the convergence law is measured on one geometry rather than on a box that moves with the grid.
The trap supplies exact anchors that reuse $\psi_0$: the density at $t=\pi$ must equal
$\rho_0(-x)$, the recurrence overlap $|\langle\psi_0|\psi(2\pi)\rangle|$ must return to 1, and
the grid-sum drift of $\langle H\rangle$ refutes independently, with the split-step propagator as
the cross-check. Close with the measured convergence law in $n$, tail-limited by the cusp, set
against the spectral convergence a smooth packet would show.

#### 5.3 [BSc] How do I follow a free Gaussian wave packet's spreading and group velocity, both analytically and by `NDSolve`?

The question names the free Gaussian, so it is the carrier, with width $\sigma_0$ and drift $k_0$
kept symbolic. Get the closed form analytically first: the C3 verdict certifies `DSolve` on the
free IVP only for the fixed packet whose $t=0$ form is $e^{-x^2/4}e^{2ix}$ ($\sigma_0=1$,
$k_0=2$), so try symbolic-parameter `DSolve` and tag it (verify at authoring), with the transform
route as the named fallback: `FourierTransform` of $\psi_0$, multiply by $e^{-ik^2t/2}$,
transform back under `Assumptions -> {sigma0 > 0, Element[{x, t}, Reals]}`, then residual-verify
the resulting closed form against the TDSE. Run the shared C3 contract with box edge
$\gtrsim k_0T+8\sigma_T$, compare pointwise, and check the moment laws $\langle x\rangle(t)=k_0t$
and $\sigma^2(t)=\sigma_0^2+t^2/(4\sigma_0^2)$ from grid sums; then exhibit the norm-blind
boundary trap deliberately, an undersized box reflecting at norm 1. The discriminating contrast
is a same-width sech packet whose rescaled density fails to collapse onto its initial shape:
self-similar spreading is a Gaussian privilege, not a free-particle law. Close with the ballistic
limit $\sigma(t)\to t/(2\sigma_0)$.

#### 5.4 [MSc] How do I implement the split-step Fourier propagator from scratch and benchmark it against `NDSolve`?

Build the Strang map $e^{-iV\,dt/2}e^{-iK\,dt}e^{-iV\,dt/2}$ from scratch on
`Fourier`/`InverseFourier` (`FourierParameters -> {1,-1}`, packed arrays, the $k$-grid convention
checked before use) and benchmark it against the shared `NDSolveValue` contract on a packet in
the symmetric quartic double well $V=(x^2-a^2)^2/2$, started as a local-harmonic-width Gaussian
displaced onto the left minimum. Pin $a=1.5$, where the barrier $a^4/2\approx2.5$ is comparable
to the local zero-point energy $\omega/2=a=1.5$: the doublet splitting is then of order $10^{-1}$
and the tunneling period $\pi/\Delta E$ of order $10^2$, inside reach of the integrator, whereas
the C2-measured deep cases $a=2$ (splitting $7.58\times10^{-4}$, period about $8\times10^3$) and
$a=2.2$ ($2.9\times10^{-5}$, about $2\times10^5$) sit three to five orders beyond the probed
windows and would be dominated by drift. Confirm the splitting with a quick eigenvalue check (the
C2 recipe) before committing run time. Verify that the pointwise disagreement between the two
integrators shrinks under $dt$ halving and goal tightening at fixed $n$, and measure the Strang
order in a coarse-$dt$ window, since the verdict flags the order ratio as contaminated at the
$10^{-13}$ floor. The discriminator is per-step norm, measured for this run rather than
transplanted: the verdict's $10^{-13}$ split-step and $10^{-5}$ `NDSolveValue` drifts are
free-packet references at $T\lesssim6$, not predictions for a longer double-well run. Close with
the $V\to0$ limit, which must reproduce the free-Gaussian closed form.

#### 5.5 [MSc] How do I evolve a state by expanding it in energy eigenstates and propagating term by term, and observe wave-packet revivals?

Expand the off-center triangular packet peaked at $x_0=L/3$ in the infinite well of width $L$:
the kink populates many modes with $c_n$ decaying only as $1/n^2$. Certify the spectrum before
leaning on it, since the whole revival argument needs $E_n$ exactly quadratic in $n$: substitute
$\psi_n=\sqrt{2/L}\sin(n\pi x/L)$ with symbolic $n$ into the stationary equation, `FullSimplify`
the residual to $0$ under `Element[n, Integers]`, and check $\psi_n(0)=\psi_n(L)=0$; then
$E_n=n^2\pi^2/(2L^2)$ gives $E_nT_{\mathrm{rev}}=2\pi n^2$ at $T_{\mathrm{rev}}=4L^2/\pi$. The
primary is the exact eigenbasis per the `Route-Table.md` cross-class hand-off: $c_n$ by
`Integrate` in closed form with `Assumptions -> Element[n, Integers]`, a truncated `Sum` bounded
by the tail mass $\sum_{n>N}|c_n|^2$. The revival is then exact by construction, so the refutable
content is the mirror identity $\psi(x,T_{\mathrm{rev}}/2)=-\psi_0(L-x)$, compared pointwise
between the truncated series and the reflected initial function, plus an independent numeric
evolution whose overlap $|\langle\psi_0|\psi(T_{\mathrm{rev}})\rangle|$ must return to 1. That
run leaves the periodic-pseudospectral primary because these walls are physical, and the
verdict's only Dirichlet evidence is a free packet striking an artificial wall, which does not
transfer; pseudospectral on a non-periodic grid with Dirichlet conditions is gated in no verdict,
so take `"DifferenceOrder" -> 4` with `MinPoints` and `MaxPoints` pinned and swept, and tag a
pseudospectral variant (verify at authoring) if the kink's high modes demand it. As a long-time
member it inherits the flagged risks, default-goal drift growing with integration length and
pseudospectral cost above $n=513$ unprobed, so the numeric leg carries a goal-and-resolution
sweep. Close with the quarter-revival two-image cat at $T_{\mathrm{rev}}/4$.

#### 5.6 [MSc] How does an oscillator coherent state evolve, so that $\langle x\rangle(t)$ traces the classical oscillation while the packet stays minimal and does not spread?

Evolve the coherent state $\alpha=2$ ($x_0=2\sqrt2$, $p_0=0$, width
$\sigma_{\mathrm{gs}}^2=\tfrac12$) beside a squeezed-displaced Gaussian with the same $(x_0,p_0)$
but $\sigma_0^2=\tfrac18$ ($r=\ln 2$) in the trap $V=x^2/2$: the same-centroid partner isolates
what "coherent" adds to "displaced", turning "does not spread" into a statement about the width
alone. Get the exact side first: for quadratic $V$ the moment system
$\langle x\rangle,\langle p\rangle,\langle x^2\rangle,\langle xp{+}px\rangle,\langle p^2\rangle$
closes and `DSolve` solves it, giving $\langle x\rangle(t)=2\sqrt2\cos t$ for both states and
$\sigma^2(t)=\sigma_0^2\cos^2t+\tfrac{1}{4\sigma_0^2}\sin^2t$. Run both under the shared C3
contract with moments as $dx$-weighted grid sums, and refute or confirm on: the coherent width
holding $\tfrac12$ for all $t$ while the squeezed width breathes between $\tfrac18$ and $2$ at
frequency 2, the two centroids agreeing pointwise, and the exact recurrence
$|\langle\psi_0|\psi(2\pi)\rangle|=1$. Close with the $r\to0$ limit collapsing the pair onto one
state, noting the edge that a strong squeeze under-resolves the narrow phase of the breath.

#### 5.7 [MSc] How do I watch a wave packet scatter off a barrier in real time and see tunneling and reflection split the packet?

Boost the compact-support raised-cosine packet
$\psi_0\propto\left(1+\cos\frac{\pi(x-x_c)}{w}\right)e^{ik_0x}$ on $|x-x_c|<w$, with $x_c=-20$,
$w=5$, $k_0=2$, at the Eckart barrier $V_0\operatorname{sech}^2(x/a)$ with $V_0=2$ and $a=1$.
That places the mean energy $E_0=k_0^2/2=2$ mid-ramp on the exact anchor
$T(E)=\sinh^2(\pi ka)\,/\,[\sinh^2(\pi ka)+\cosh^2(\tfrac\pi2\sqrt{8V_0a^2-1})]$ with
$k=\sqrt{2E}$, valid for $8V_0a^2>1$: here $\sinh^2(2\pi)=7.17\times10^4$ against
$\cosh^2(\tfrac\pi2\sqrt{15})=4.83\times10^4$, so $T=0.60$ and both fragments are large enough to
see and to weigh. The barrier vanishes nowhere, so compact support does not give zero initial
overlap; it gives a quantitatively negligible one, $V/V_0=\operatorname{sech}^2(d/a)\approx
4e^{-2d/a}=3.7\times10^{-13}$ at the nearest support edge $d=15$, an energy shift far below
tolerance, while the real payoff is an unambiguous late-time mass split and a clean initial
momentum amplitude. Run the shared C3 contract with box edge $\gtrsim k_0T+8\sigma_T$ on both
sides, since both fragments travel, and the split-step cross-check. The refuting check is the
late-time transmitted mass (grid sum beyond the barrier) against
$\int T(k)\,|\tilde\psi_0(k)|^2\,dk$ built from the same $\psi_0$ via `FourierTransform`, with
$R+T=1$ from the same sums. The C3 verdict flags 5.7's ringing as its unmeasured open risk: the
smooth Eckart barrier retires the potential-side Gibbs ringing while keeping the exact anchor,
but the packet's $C^1$ support seams still seed high-$k$ content, so a resolution sweep in $n$ of
the seam ringing is mandatory at authoring. Close with $T\to1$ at $E\gg V_0$ set against the
deep-tunneling suppression.

#### 5.8 [MSc] How do I construct the propagator (kernel) $K(x,t;x',0)$ for the free particle and for the oscillator (the Mehler kernel)?

Construct the free kernel $K_0=(2\pi it)^{-1/2}e^{i(x-x')^2/(2t)}$ from the momentum integral
$\frac{1}{2\pi}\int e^{ik(x-x')-ik^2t/2}\,dk$ by `Integrate` with a $t\to t-i\epsilon$ regulator
and `Limit`, and the Mehler kernel
$K_{\mathrm{osc}}=(2\pi i\sin t)^{-1/2}\exp\!\left\{\tfrac{i}{2\sin t}\left[(x^2+x'^2)\cos t-2xx'\right]\right\}$
by closing the Wick-rotated eigen-sum $\sum_n e^{-(n+1/2)\tau}\psi_n(x)\psi_n(x')$ with `Sum`
(`HermiteH`) and continuing $\tau\to it$; if `Sum` refuses the bilinear Hermite series, posit the
closed form and certify it instead by the imaginary-time residual
$\partial_\tau K=\tfrac12\partial_{xx}K-\tfrac12x^2K$ and the $\tau\to0^+$ delta limit. At real
$t$ the eigen-sum oscillates without converging pointwise, which is why the comparison runs
through imaginary time (C0 per the `Route-Table.md` C0 row). Each kernel earns its name three
ways, any failure refuting it: the TDSE residual in $(x,t)$ is identically $0$ under `D` +
`FullSimplify`; the composition law $\int K(x,t_1;x'')\,K(x'',t_2;x')\,dx''=K(x,t_1{+}t_2;x')$
closes; and the $t\to0^+$ action on a test packet returns the packet, the $\delta(x-x')$ limit
taken through a test function, never as a bare limit. Add a numeric leg: propagate one packet by
explicit kernel quadrature at fixed $t$ and compare against a locally rebuilt split-step Fourier
map. Then let the free kernel act on the Berry-Balazs Airy beam (`AiryAi`), which propagates as
$\operatorname{Ai}(x-t^2/4)\,e^{it(x-t^2/6)/2}$, an accelerating profile out of a free kernel.
Close at the caustics $t=n\pi$, where the Mehler prefactor diverges and the kernel degenerates to
mirrored delta transport, with the $\sin t\to t$ free-kernel limit as the bridge between the two.

#### 5.9 [BSc] How do I verify Ehrenfest's theorem numerically, that $\langle x\rangle$ and $\langle p\rangle$ obey the classical equations of motion?

Evolve the asymmetric two-lobe packet, the normalized unequal-weight sum
$e^{-(x-1.2)^2}+\tfrac12 e^{-(x+0.6)^2}$, in the quartic well $V=x^4/4$, where Ehrenfest's actual
content is the difference between $\frac{d\langle p\rangle}{dt}=-\langle V'(x)\rangle=
-\langle x^3\rangle$ and the classical $-V'(\langle x\rangle)=-\langle x\rangle^3$. That gap is
$\langle x^3\rangle-\langle x\rangle^3=3\langle x\rangle\operatorname{Var}(x)+\mu_3$, so what the
example must supply is a large variance at nonzero $\langle x\rangle$ together with a nonzero
third central moment $\mu_3$; skewness alone is not the mechanism, and even a symmetric packet
displaced off the origin already opens the first term. Run the shared C3 contract with the box
beyond the packet's classical turning points, moments as $dx$-weighted grid sums and
$\langle p\rangle$ from the spectral derivative on the $k$-grid. The machinery's refuter is the
harmonic control, the same packet in $V=x^2/2$, where quadratic $V$ closes Ehrenfest on
$\langle x\rangle$ exactly, so any gap there is an error. On the quartic run both residuals
$\frac{d\langle x\rangle}{dt}-\langle p\rangle$ and
$\frac{d\langle p\rangle}{dt}+\langle x^3\rangle$ must vanish to tolerance, allowing for the
finite-difference noise floor of $d/dt$ on a sampled moment series, while the classical curve
departs. Close with the narrow-packet limit shrinking both terms of the gap: classicality is a
property of the state's spread, not of the theorem.

#### 5.10 [BSc] How do I confirm that norm and probability current are conserved along an `NDSolve` evolution (a numerical-fidelity check)?

Launch the boosted supergaussian $\psi_0\propto e^{-(x-x_0)^4}e^{ik_0x}$ from $x_0=-8$ across the
reflectionless well $V=-\operatorname{sech}^2x$, the $\lambda=1$ Poschl-Teller well with the
single bound state $\psi_b=\operatorname{sech}(x)/\sqrt2$ at $E_b=-\tfrac12$: $|t(k)|=1$ for
every momentum component, so any reflected lobe above tolerance is guaranteed spurious. Derive
the current symbolically before the run, from the TDSE for generic $\psi$ and real $V$, writing
$\rho=\psi\bar\psi$ (never differentiate $|\psi|^2$ literally, since `D[Abs[..]]` yields no
delta, a C3-cited trap) to get $\partial_t\rho+\partial_x j=0$ with
$j=\operatorname{Im}(\bar\psi\,\partial_x\psi)$. Run the shared C3 contract with box edge
$\gtrsim k_0T+8\sigma_T$ and grade three checks of increasing strength on the same run: the norm
grid-sum drift; the continuity residual on the grid, shrinking under refinement in $n$; and
station-flux bookkeeping, $\int_0^T j(x_s,t)\,dt$ against the probability mass transferred past
$x_s$. Reflectionless does not mean everything transmits: the well captures its bound component
permanently, so predict $T=1-|\langle\psi_b|\psi_0\rangle|^2$ by `Integrate` from the same
definitions, a truthful and strictly stronger anchor than $T\to1$. Exhibit the norm-blindness
lesson live, an undersized box reflecting at norm 1, so norm alone certifies nothing about walls.
Close on that ranking: a fidelity check that cannot fail teaches nothing.

## Part 6 Plan: Scattering in one dimension

Questions: 6 (6.1 through 6.6). Class census: C4 for all six questions, per the `Route-Table.md`
class census row C4 and the full C4 verdict block; no other class appears in this part.

### Common ground

Every question here is a fixed-energy scattering boundary-value problem: the stationary Schrodinger
equation at a continuum energy $E$, with connection amplitudes, not eigenvalues, as the unknowns.
The probability current $j=\operatorname{Im}(\bar\psi\,\partial_x\psi)$ defines the physical
coefficients as current ratios, $R=|r|^2$ and $T=(k'/k)\,|t|^2$, with $k=\sqrt{2E}$ the incident and
$k'$ the far-side asymptotic wavevector; $|t|^2$ alone is not $T$ whenever the two media differ. All
dynamics enters through matching: continuity of $\psi$ and $\psi'$ at every interface of a
piecewise-constant potential, replaced for a point interaction $\lambda\,\delta(x)$ by the jump
condition $\psi'(0^{+})-\psi'(0^{-})=2\lambda\,\psi(0)$ with $\psi$ continuous. The scattering
matrix collects the amplitudes, outgoing $=S\,\cdot$ incoming with
$S=\begin{pmatrix} r & t' \\ t & r' \end{pmatrix}$; current conservation is unitarity
$S^{\dagger}S=\mathbb{1}$, and the continuation $k\to i\kappa$ turns poles of $S$ into bound states.
Standing C4 inheritances (Route-Table.md, C4 verdict): answers are phrased region by region because
DSolve on a Piecewise coefficient silently echoes its input (probe 3); the primary route is symbolic
matching, Solve on the continuity system plus ComplexExpand plus FullSimplify under positivity
assumptions with wavevectors substituted last; the transfer matrix is the independent cross-check;
and the standing verification is current-based unitarity, $R+T=1$ derived from currents, never from
amplitude moduli alone.

### Per-question entries

#### 6.1 [BSc] How do I compute the reflection and transmission coefficients at a potential step?

Set up the step $V=V_0\,\theta(x)$ (a heterojunction) region by region, because DSolve refuses the
piecewise potential silently: for $E>V_0$ match $\psi$ and $\psi'$ at $x=0$ with transmitted
wavevector $k'=\sqrt{2(E-V_0)}$, Solve for $r,t$, and reduce with ComplexExpand and FullSimplify
under $k,k'>0$. The point of the pinned two-regime example is that only currents survive as
probabilities: near threshold $|t|^2=4k^2/(k+k')^2\to4$ while the current-based
$T=(k'/k)|t|^2\to0$, so build $R$ and $T$ from $j=\operatorname{Im}(\bar\psi\,\partial_x\psi)$ on
the solved amplitudes and verify $R+T=1$ there, a check any amplitude-only definition fails for
every $E>V_0$. Take the limits $E\to V_0^{+}$ ($T\to0$) and $E\gg V_0$ ($T\to1$), spot-check
numerically with the single-interface transfer matrix, then rerun the matching for $E<V_0$: the
evanescent branch $e^{-\kappa x}$ with $\kappa=\sqrt{2(V_0-E)}$ carries zero current, so $R=1$
exactly with penetration depth $1/\kappa$. Close on the contrast the ledger pins: at a step,
amplitude moduli mislead and currents do not.

#### 6.2 [BSc] How do I compute the tunneling transmission $T(E)$ through a rectangular barrier?

Set up plane waves in the three regions of a barrier of height $V_0$ and width $L$, match $\psi$
and $\psi'$ at both walls, and Solve the linear system for the transmission amplitude; $T=|t|^2$
comes out in closed form via ComplexExpand and FullSimplify under $k,\kappa,L>0$ (region by
region, because DSolve refuses the piecewise potential silently), with $k=\sqrt{2E}$ and
$\kappa=\sqrt{2(V_0-E)}$ substituted last, and the result provably equal to the textbook
$[1+\tfrac{V_0^2\sinh^2\kappa L}{4E(V_0-E)}]^{-1}$ (probe 1: difference 0). Verify $R+T=1$ from
the probability currents, recover the over-barrier resonances by the continuation $\kappa\to iq$
(with $T=1$ exactly at $\sin qL=0$), and cross-check the whole thing against a transfer-matrix
product (probe 2: difference 0, $\det M=1$). Close with the opaque limit
$T\approx16\,\tfrac{E(V_0-E)}{V_0^2}\,e^{-2\kappa L}$ by Series, and the threshold edge
$E\to V_0$, where Limit gives $T\to(1+V_0L^2/2)^{-1}$.

#### 6.3 [MSc] How do I exhibit transmission resonances in a well (the Ramsauer-Townsend effect)?

Take the rectangular well of depth $V_0$ and width $a$ and continue the barrier matching by
$\kappa\to ik'$ with inside wavevector $k'=\sqrt{2(E+V_0)}$ (probe 1 of the C4 verdict: $T=1$
exactly at $\sin k'a=0$). Predict the resonance ladder first, Solve/Reduce under positivity giving
$E_n=\tfrac{n^2\pi^2}{2a^2}-V_0$, with $V_0$ and $a$ placed so the $n=1$ resonance sits at low
energy (Ramsauer-Townsend); then substitute $E_n$ back into the closed-form $T(E)$, reused and
never re-typed, and FullSimplify $T-1$ to 0, with strict $T<1$ off resonance. Verify $R+T=1$ from
the currents throughout and bracket the predicted $E_1$ with transfer-matrix numeric values of $T$.
Close at the two ends, the strong-scattering dip as $E\to0$ between resonances and the deep-well
limit where the resonances sharpen, against the barrier of 6.2, whose closed form has no $T=1$
point below $V_0$.

#### 6.4 [MSc] How do I assemble the transfer-matrix method for a piecewise-constant potential?

Build the interface and free-propagation $2\times2$ matrices for a symmetric double barrier
(heights $V_0$, widths $a$, gap $b$: a resonant-tunneling diode) and assemble the total $M$ with
Inverse and Dot; the transfer matrix is primary here per the C4 verdict, because the cell structure
is the physics, and the symbolic product stays tractable (probe 2: LeafCount 9775 simplifies to 186
in under a minute). Read $T=1/|M_{11}|^2$ off the product and locate the inter-barrier resonance,
invisible in either factor: on resonance $T=1$ for the symmetric structure, off resonance
$T\sim T_1^2$, far sharper than either barrier alone. Verify $\det M=1$ per factor and for the
product (flux conservation, refutable factor by factor), confirm the single-barrier factor
reproduces 6.2's closed form and that $R+T=1$ holds from the currents, and pin the numerics against
the direct 8-unknown Solve matching on the same geometry (probe 2: both give
$T=0.2965025061851854$). Close with the opaque limit, where the resonance width collapses toward a
bound-state-like line.

#### 6.5 [MSc] How do I build the one-dimensional scattering matrix, verify its unitarity, and identify bound states as poles?

For the finite square well of depth $V_0$ and width $a$, run the symbolic matching for both
incidence directions (Solve, ComplexExpand, FullSimplify under $k,k'>0$, region by region as
everywhere in this part) and assemble $S=\begin{pmatrix} r & t \\ t & r \end{pmatrix}$; then
FullSimplify $S^{\dagger}S-\mathbb{1}$ to the zero matrix, which is the current-conservation
statement itself, its diagonal the standing $R+T=1$. Continue $k\to i\kappa$ and show the poles of
$S$ land on the even and odd conditions $\kappa=k'\tan(k'a/2)$ and $\kappa=-k'\cot(k'a/2)$ with
$k'=\sqrt{2(E+V_0)}$; this pole step is algebra on the probed closed form, and the C4 verdict's
open risks flag it as not itself probed, so the continuation must be shown explicitly in a cell,
never cited. The refuting check is independence: match the bound sector directly
(evanescent-outside ansatz, Solve) and locate both conditions with FindRoot at pinned parameters;
the two routes must give the same energies, and bound states enter as poles of the continued
denominator, never through an eigensolver. Close with the limits: $V_0\to0$ gives
$S\to\sigma_x=\begin{pmatrix} 0 & 1 \\ 1 & 0\end{pmatrix}$, pure transmission with $r=0$ and
$t=1$, since on this convention the identity would mean $r=1$, a hard wall rather than free
space; and a pole arriving at threshold $\kappa\to0$ announces a new bound state.

#### 6.6 [MSc] How do I check Levinson's theorem relating the phase shift to the number of bound states?

Give the attractive delta well $V=\lambda\,\delta(x)$, $\lambda=-g<0$ (a contact interaction),
entirely to the symbolic jump condition $\psi'(0^{+})-\psi'(0^{-})=2\lambda\,\psi(0)$ with $\psi$
continuous: nothing numeric may touch the potential, since NIntegrate returns 0. on a point measure
and `D[Abs[x],{x,2}]` yields no delta (C4 traps (c) and (d)). Solve the jump system for $r,t$,
extract the even-channel phase $\delta_{+}(k)=\arctan(g/k)$ via Arg/ArcTan, and take Limit at both
ends, $\delta\to0$ as $k\to\infty$ and the Levinson value as $k\to0$; the odd channel is exactly
unscattered ($\delta_{-}\equiv0$), the discriminating contrast the zero-range interaction gives for
free. The bound-state count comes from exact delta-well algebra, not a C4 route (the verdict's
members-sanity note): the continued $t=ik/(ik-\lambda)$ has its single pole at $k=ig$, hence
exactly one bound state at $E=-g^2/2$. Cross-check the phase as the eigenphase of $S$ (the
symmetric eigenvector $(1,1)$ gives $r+t$, the antisymmetric $(1,-1)$ gives $r-t$, which is the
free value $-1$ here, matching $\delta_{-}\equiv0$) against the jump-condition $\delta_{+}$ at a
pinned $k$, verify $R+T=1$ from the
currents, and close by cashing the theorem: the phase winding at threshold counts the spectrum.
Concern: the pinned $\delta(0)-\delta(\infty)=\pi$ is the winding of the even eigenphase
$2\delta_{+}$ of $S$; the phase shift itself obeys the half-corrected 1D Levinson
$\delta_{+}(0)-\delta_{+}(\infty)=\pi\,(n_{+}-\tfrac12)=\pi/2$ for one bound state and no
zero-energy resonance, so authoring must pin the convention before quoting $\pi$.

## Part 7 Plan: Approximation methods

Ten questions. Class census per the class census table in `Route-Table.md`: C0 (no differential
equation) 7.1, 7.2, 7.3, 7.4, 7.10; C9 (WKB / semiclassical asymptotics) 7.5, 7.6, 7.9; C5
(truncated-basis ODE-IVP systems) 7.7, 7.8.

### Common ground

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

#### 7.1 [BSc] How do I compute first- and second-order energy shifts by nondegenerate time-independent perturbation theory (the anharmonic oscillator)?

Perturb the oscillator with $\lambda x^4$ and compute $E_n^{(1)}$ and $E_n^{(2)}$ for general $n$:
get the $x^4$ matrix elements by Integrate over HermiteH eigenfunctions and again by ladder
algebra (the two must agree, and only $m \in \{n, n \pm 2, n \pm 4\}$ survive), so first order is
the known $\tfrac{3\lambda}{4}(2n^2 + 2n + 1)$ and second order a four-term Sum in closed form.
Benchmark against the truncated Fock-basis Hamiltonian (SparseArray diagonal plus banded $x^4$
block, Eigenvalues, $N$-swept): the residual against the numeric level must fall as $\lambda^3$
under $\lambda$-halving, a scaling a wrong second-order sum cannot fake. Close by watching the
orders grow with $n$ and $\lambda$: the series is asymptotic, and the coupling where it stops
helping is located, not asserted.

#### 7.2 [MSc] How do I resolve a degeneracy by degenerate perturbation theory?

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

#### 7.3 [BSc] How do I get a Rayleigh-Ritz variational upper bound on the ground-state energy from a Gaussian trial function?

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

#### 7.4 [MSc] How do I run the linear variational (Ritz) method in a finite basis, reducing the problem to a generalized matrix eigenproblem?

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

#### 7.5 [BSc] How do I quantize a smooth well by the WKB (Bohr-Sommerfeld) condition?

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

#### 7.6 [MSc] How do I compute a WKB tunneling rate through a barrier (the Gamow factor)?

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

#### 7.7 [MSc] How do I compute transition rates by time-dependent perturbation theory and recover Fermi's golden rule with a continuum density of states?

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

#### 7.8 [MSc] How do I apply the sudden and the adiabatic approximations to a continuous system?

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

#### 7.9 [MSc] How do I match WKB solutions across a turning point with the Airy connection formulas?

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

#### 7.10 [MSc] How do I get a variational *lower* bound on the ground-state energy (Temple's inequality)?

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

## Part 8 Plan: General theorems and structural methods

Six questions. Class census (per the class census table in `Route-Table.md`): C0, no differential
equation, for 8.1, 8.2, 8.3, 8.5, 8.6; C9, WKB / semiclassical asymptotics, for 8.4.

### Common ground

The part rests on identities every stationary state must obey, and on operator constructions that
generate spectra without solving a differential equation. On a normalized eigenstate $H\psi=E\psi$
with $\langle O\rangle=\int\bar\psi\,O\psi$, the virial theorem fixes
$2\langle T\rangle=\langle\vec r\cdot\nabla V\rangle$; a parameter in $H_\lambda$ moves the energy
by the Hellmann-Feynman theorem $\partial E/\partial\lambda=\langle\partial H_\lambda/\partial\lambda\rangle$;
the double commutator $[\hat x,[\hat H,\hat x]]=1$ (natural units) has the Thomas-Reiche-Kuhn sum
rule $\sum_n f_{kn}=1$, $f_{kn}=2(E_n-E_k)\vert\langle n\vert\hat x\vert k\rangle\vert^2$, as its
diagonal matrix element; at large quantum number the locally averaged stationary density approaches
the classical dwell-time density $\rho_{\mathrm{cl}}(x)=2/(T_{\mathrm{cl}}\,v(x))$, for the
oscillator the arcsine law $1/(\pi\sqrt{A^2-x^2})$; the factorization $H=A^\dagger A+E_0$ with
$A=\tfrac{1}{\sqrt2}(\tfrac{d}{dx}+W)$ and $W=-\psi_0'/\psi_0$ builds partner potentials
$V_\pm=\tfrac12(W^2\pm W')+E_0$ sharing every level except the annihilated ground state; and when
the partner reproduces its own shape, $V_+(x;a_0)=V_-(x;a_1)+R(a_1)$, the spectrum telescopes to
$E_n=\sum_{k=1}^{n}R(a_k)$ with no equation solved.

#### 8.1 [BSc] How do I verify the virial theorem $2\langle T\rangle=\langle \vec r\cdot\nabla V\rangle$ on a stationary state?

Compute $\langle T\rangle$, $\langle V\rangle$, and $\langle\vec r\cdot\nabla V\rangle$ as three
independent computations on two power-law systems, so the theorem's exponent dependence is what
gets tested (C0 per `Route-Table.md`): the hydrogen $2p$ reduced radial
$u_{21}(r)=\tfrac{1}{2\sqrt6}\,r^2e^{-r/2}$ via exact Integrate,
$\langle T\rangle=\int_0^\infty(\tfrac12u'^2+\tfrac{1}{r^2}u^2)\,dr=\tfrac18$ and
$\langle V\rangle=-\tfrac14$ with $\langle T\rangle+\langle V\rangle=E_2=-\tfrac18$ as a second
identity, and the ground state of the quartic well $V=x^4/4$, which has no closed form but stays
inside C0 as a truncated harmonic number basis: ladder matrices $\hat x=(a+a^\dagger)/\sqrt2$,
$\hat p=i(a^\dagger-a)/\sqrt2$ through SparseArray and Eigensystem, the three expectations as
separate quadratic forms swept in the truncation $N$. The ratios $2\langle T\rangle=-\langle V\rangle$
and $2\langle T\rangle=4\langle V\rangle$ must then emerge, never be imposed, with an NIntegrate
spot-check on the Coulomb integrals. Close by evaluating the same combination on a 50/50
superposition of two hydrogen levels, where it fails: stationarity is load-bearing.

#### 8.2 [MSc] How do I apply the Hellmann-Feynman theorem to get $\partial E/\partial\lambda$ from an expectation value?

Treat $l$ as a continuous Hellmann-Feynman parameter in the radial family
$H_l=-\tfrac12\partial_r^2+\tfrac{l(l+1)}{2r^2}-\tfrac1r$ (C0 per `Route-Table.md`), stating the
continuation honestly: at fixed node count $n_r$ the normalizable eigenfunctions exist for real
$l>-\tfrac12$ with $E=-\tfrac{1}{2(n_r+l+1)^2}$, the physical integer-$l$ spectrum is the
restriction of this analytic family, and the derivative is taken at fixed $n_r$, never fixed $n$.
D in $l$ then gives $\langle\tfrac{2l+1}{2r^2}\rangle=\tfrac{1}{(n_r+l+1)^3}$, hence
$\langle 1/r^2\rangle=\tfrac{1}{n^3(l+1/2)}$ with no integral computed; close the loop by an
independent Integrate of $\int_0^\infty u_{n_r,l}^2\,r^{-2}\,dr$ on the C1-verdict radial forms
$u_{n,l}=r^{l+1}e^{-r/n}L_{n-l-1}^{2l+1}(2r/n)$ (LaguerreL, FullSimplify under $r>0$,
$l>-\tfrac12$): fully symbolic in $l$ on the nodeless family $u_{0,l}\propto r^{l+1}e^{-r/(l+1)}$,
where both sides come out $(l+1)^{-3}$, and exact at the nodal instance $(n_r,l)=(1,1)$, a fixed
instance because the symbolic-$(n,l)$ norm is recorded unprobed in the C1 verdict. Spot-check one
case with NIntegrate and close with the large-$l$ circular-orbit limit,
$\langle 1/r^2\rangle\to n^{-4}$ (verify at authoring).

#### 8.3 [MSc] How do I verify the Thomas-Reiche-Kuhn oscillator-strength sum rule?

Verify the rule where convergence is genuine content: the infinite well on $(0,L)$ from the ground
state (C0 per `Route-Table.md`), with dipole elements $x_{1n}$ from Integrate, symbolic in $n$,
nonzero only for even $n$ by parity and falling as $n^{-3}$. Build
$f_{1n}=2(E_n-E_1)\vert x_{1n}\vert^2$ reusing those elements and watch the partial sums
$S_N=\sum_{n\le N}f_{1n}$ from Sum converge to 1, led by $f_{12}=\tfrac{256}{27\pi^2}\approx0.961$,
with the tail $f_{1n}=O(n^{-4})$ and the exact closure $\lim_{N\to\infty}S_N=1$ attempted
symbolically (verify at authoring). Contrast with the oscillator, where saturation must fall out
rather than be asserted: Hermite-Gaussian Integrate gives $f_{0n}=0$ for every $n\neq1$ and
$f_{01}=1$ exactly, and from an excited state $k$ the negative emission strength $f_{k,k-1}=-k$
joins $f_{k,k+1}=k+1$ to keep the sum at 1. Close on that pair: one system spreads the strength
over an infinite even-parity ladder, the other concentrates it in a single line.

#### 8.4 [MSc] How do I exhibit the correspondence principle, the classical limit of stationary states at large quantum number?

Make the limit quantitative on the oscillator $n=30$ (C9 per `Route-Table.md`): the exact
eigenfunction via HermiteH with $E_{30}=\tfrac{61}{2}$, against the classical arcsine density
$\rho_{\mathrm{cl}}(x)=\tfrac{1}{\pi\sqrt{A^2-x^2}}$ with $A=\sqrt{2E_{30}}=\sqrt{61}$, normalized
exactly on $(-A,A)$; both legs are certified machinery in the C9 verdict, whose cross-check, the
C2-certified FD recipe (SparseArray Band tridiagonal, Eigenvalues with Arnoldi shift, Richardson),
independently reproduces the $n=30$ eigenpair. The C9 verdict flags the local-averaging window as
an open authoring decision with no kernel evidence: candidates are the integral mean
$\bar\rho(x)=\tfrac{1}{2w}\int_{x-w}^{x+w}\vert\psi_{30}\vert^2$ or a grid MovingAverage, with $w$
tied to the local de Broglie wavelength $2\pi/p_{\mathrm{cl}}(x)$; record the choice and a
$w$-sensitivity check, never silently. The principle itself is the refutable object: the windowed
error $\max\vert\bar\rho-\rho_{\mathrm{cl}}\vert$ must shrink from $n=10$ to $n=30$ with the same
definitions reused at both $n$, while the raw density does not converge, its oscillation staying
comparable to the density itself. Close at the turning points, where $\rho_{\mathrm{cl}}$ diverges
while the quantum density stays finite with tunneling tails beyond: the edge where the
correspondence stops.

#### 8.5 [MSc] How do I factorize a Hamiltonian as $H=A^\dagger A+E_0$ and build its supersymmetric partner potential?

Factorize the infinite well on $(0,L)$ by pure algebra (C0 per `Route-Table.md`):
$W=-\psi_1'/\psi_1=-\tfrac{\pi}{L}\cot\tfrac{\pi x}{L}$, $A=\tfrac{1}{\sqrt2}(\partial_x+W)$,
$A^\dagger=\tfrac{1}{\sqrt2}(-\partial_x+W)$. Expand $A^\dagger A$ on a generic $f(x)$ and
FullSimplify: $-\tfrac12f''-E_1f$ with $E_1=\tfrac{\pi^2}{2L^2}$ falls out rather than being
imposed, so $H_1=A^\dagger A+E_1$, and the reversed product $AA^\dagger$ yields the partner
$V_2=\tfrac{\pi^2}{L^2}\csc^2\tfrac{\pi x}{L}$, which looks nothing like a box yet must carry the
box spectrum minus its ground level. Check isospectrality against exact solvability, not against
$H_1$ itself: $V_2$ is the trigonometric Poschl-Teller well
$\tfrac{\lambda(\lambda-1)}{2}\tfrac{\pi^2}{L^2}\csc^2\tfrac{\pi x}{L}$ at $\lambda=2$ with exact
spectrum $\tfrac{\pi^2(n+2)^2}{2L^2}$ (verify at authoring), equal to the well levels for $n\ge2$;
certify the mapped eigenfunctions by the residual $H_2(A\psi_n)-E_n(A\psi_n)=0$ with
$\psi_n=\sqrt{2/L}\sin(n\pi x/L)$ at symbolic integer $n$ (fixed-$n$ instances if the kernel
stalls), with the C2-certified FD recipe on $V_2$ under Dirichlet walls as the numeric fallback.
Close on the deleted level: $A\psi_1=0$, the annihilated ground state, is exactly the state the
partner cannot have.

#### 8.6 [MSc] How do I use shape invariance to obtain the oscillator and Poschl-Teller spectra purely algebraically?

Obtain both spectra from the recursion alone, no differential equation anywhere (C0 per
`Route-Table.md`): oscillator $W=x$ with constant shift $R=1$, and Poschl-Teller
$W=\lambda\tanh x$ with $V_-=\tfrac12\bigl(\lambda^2-\lambda(\lambda+1)\operatorname{sech}^2x\bigr)$
and partner parameter $\lambda\to\lambda-1$. Verify the shape-invariance identity
$V_+(x;a_0)=V_-(x;a_1)+R(a_1)$ as a zero FullSimplify residual at symbolic $x$ and $\lambda$, then
telescope $E_n=\sum_{k=1}^{n}R(a_k)$ with Sum: $E_n=n$ above the oscillator ground energy, and
$E_n=\tfrac12\bigl(\lambda^2-(\lambda-n)^2\bigr)$ on the $\lambda$ ladder. The recursion output
must match independently known exact spectra: the oscillator anchor
$\{\tfrac12,\tfrac32,\tfrac52,\tfrac72\}$ is kernel-evidenced in the C1 verdict (DEigensystem
domain form), and the $\lambda=2$ ladder must reproduce the exactly solvable pair $E_0=-2$,
$E_1=-\tfrac12$ of $V=-3\operatorname{sech}^2x$ after the asymptote shift, with the C2-certified
FD recipe on $V_-$ as the numeric fallback. Close on termination: the ladder stops when
$\lambda-n\le0$, so at $\lambda=2$ it predicts exactly two bound states, a count the algebra must
get right, not just spacings; pairing an infinite ladder with a terminating one shows the
recursion carries the bound-state count too.

## Part 9 Plan: Periodic potentials and band structure

Part 9, "Periodic potentials and band structure", 6 questions. Class census (`Route-Table.md`):
C1: 9.1; C4: 9.2, 9.3; C0: 9.4, 9.5; C3: 9.6.

A potential with period $a$, $V(x+a)=V(x)$, commutes with the lattice translation, so energy
eigenstates take the Bloch form $\psi_k(x)=e^{ikx}u_k(x)$ with periodic $u_k(x+a)=u_k(x)$ and
quasimomentum $k$ in the Brillouin zone $(-\pi/a,\pi/a]$; the spectrum organizes into bands
$E_n(k)$ separated by gaps. For a piecewise cell the one-cell transfer matrix $M(E)$ with
$\det M=1$ carries the whole spectrum through $\cos ka=\tfrac12\operatorname{Tr}M(E)$, allowed
bands being $\lvert\operatorname{Tr}M\rvert\le2$. The density of states per site is
$g(E)=\frac{a}{\pi}\lvert dE/dk\rvert^{-1}$ with $\int_{\mathrm{band}}g\,dE=1$ and van Hove
$1/\sqrt{\lvert E-E_{\mathrm{edge}}\rvert}$ edges; Wannier functions
$w_m(x)=\frac{a}{2\pi}\int_{\mathrm{BZ}}e^{-ikma}\psi_k(x)\,dk$ localize when the Bloch gauge is
smooth; a deep lattice gives the tight-binding band $E(k)=\varepsilon-2J\cos ka$, and a static
tilt $Fx$ replaces transport by the Wannier-Stark ladder with spacing and Bloch frequency
$\omega_B=Fa$.

#### 9.1 [BSc] How do I state and use Bloch's theorem, writing eigenstates as $e^{ikx}u_k(x)$ with periodic $u_k$?

Pinned example: the cosine lattice $V(x)=V_0\cos2x$ (period $\pi$, the optical-lattice standard)
at $V_0=1$; the substitution $\psi''+(2E-2V_0\cos2x)\psi=0$ maps it onto the Mathieu equation
with $a=2E$, $q=V_0$, so the entry rides the probed C1 Mathieu gate: band energies
$E_n(k)=\tfrac12\,$`MathieuCharacteristicA`/`MathieuCharacteristicB` (certified against a
periodic-cell `NDEigenvalues` solve with `PeriodicBoundaryCondition` at the $3\times10^{-6}$
level) and Bloch states from the `MathieuC`/`MathieuS` Floquet pair (certified against direct
`NDSolve` integration near $10^{-11}$). State the theorem, then verify the form instead of
asserting it: assemble $\psi_k$ at an interior $k$ (say $k=0.3$ in zone-edge units), extract
$u_k(x)=e^{-ikx}\psi_k(x)$, and check $u_k(x+\pi)-u_k(x)$ at sample points, a check that fails
outright if the Floquet assembly is wrong. Exact anchor and refuting limit reusing the same band
function: at $q=0$ the characteristic reduces to $\nu^2$, so `FullSimplify`/`Limit` must return
the free dispersion $E(k)=k^2/2$ exactly. Close at the zone edge: a small-$q$ `Series` of the
first characteristic gap gives $\Delta E\approx V_0$, the nearly-free-electron gap opening
linear in lattice depth.

#### 9.2 [MSc] How do I solve the Kronig-Penney model and plot the allowed and forbidden energy bands?

Pinned example: the rectangular lattice, wells of width $a=1$ between barriers of height $V_0=8$
and width $b=1/4$, period $d=a+b$. Work region by region with $e^{\pm iqx}$ in wells and
$e^{\pm\kappa x}$ under barriers ($q=\sqrt{2E}$, $\kappa=\sqrt{2(V_0-E)}$ substituted last),
because `DSolve` on the `Piecewise` potential silently echoes unevaluated (C4 trap): build the
one-cell matrix $M(E)$ by `Dot` from the probed interface matrices, confirm `Det`$[M]=1$, and
read off $\cos kd=\tfrac12\operatorname{Tr}M(E)$, simplified with `ComplexExpand` and
`FullSimplify` under $q,\kappa,a,b>0$. The Bloch trace step is flagged unprobed algebra in the
C4 verdict, so certify it independently: the identical condition must emerge from the direct
continuity-plus-Bloch-phase `Solve` on one cell, and the `FullSimplify` of the difference to $0$
can refute either derivation. `Plot` the trace function with the $[-1,1]$ strip marking allowed
bands, invert to the band diagram $E_n(k)$ with `FindRoot` on a $k$-grid, and pin the first two
bands and their gap from the edge roots $\lvert\tfrac12\operatorname{Tr}M\rvert=1$ (the
$E>V_0$ branch via the continuation $\kappa\to iq'$). Refuting limits reusing the trace:
$V_0\to0$ collapses onto the folded free parabola $E=k^2/2$, and $\kappa b\gg1$ flattens each
band onto an isolated-well level with width $\sim e^{-2\kappa b}$. Close with that opaque-barrier
scaling: bandwidth is tunneling.

#### 9.3 [MSc] How do I treat a Dirac comb (a lattice of delta potentials) and find its band edges?

Pinned example: the attractive Dirac comb $V(x)=-g\sum_{n}\delta(x-nd)$ with $g=1$, $d=2$. The
delta never touches numerics (`NIntegrate` silently returns $0.$ on a point measure, and `D` of
`Abs` yields no `DiracDelta`, both binding here per C4): it enters only through the symbolic
jump condition $\psi'(nd^{+})-\psi'(nd^{-})=-2g\,\psi(nd)$, derived once by `Integrate` across
the point. Plane waves between sites plus continuity, jump, and Bloch phase go through `Solve`
to the band condition $\cos kd=\cos qd-\frac{g}{q}\sin qd$ for $E=q^2/2>0$, continued by
$q\to i\kappa$ to $\cos kd=\cosh\kappa d-\frac{g}{\kappa}\sinh\kappa d$ for the negative band;
band edges are the `FindRoot`/`Reduce` roots at $\cos kd=\pm1$. Cross-checks that could refute:
the delta-interface transfer matrix (`Dot` per cell, `Det` $=1$) must reproduce the same trace,
and the ledger's connecting limit, a rectangular-cell trace re-derived inside this entry for
attractive cells and collapsed by `Limit` $b\to0$ at fixed area $V_0b=-g$, must land exactly on
the comb condition. Anchor: as $d\to\infty$ the negative band pinches onto the isolated
delta-well level $E=-g^2/2$ with exponentially small width. Close there: a band is bound-state
hybridization, its width the tunneling overlap of neighboring wells.

#### 9.4 [MSc] How do I compute the density of states within a band?

Pinned example: the lowest band of the cosine lattice $V_0\cos2x$ at $V_0=1$, dispersion
$E(k)=\tfrac12\,$`MathieuCharacteristicA`$[k,V_0]$ for $k\in(0,1)$ in zone-edge units, a C0
entry riding already-probed Mathieu evaluations with no new equation solved. The density of
states per site is $g(E)=\frac{a}{\pi}\lvert dE/dk\rvert^{-1}$, normalized so
$\int_{\mathrm{band}}g\,dE=1$; the symbolic $\nu$-derivative of `MathieuCharacteristicA` is
uncertain (verify at authoring), with the fallback a fine `Table` $k$-grid and centered
`Differences`, keeping energies named `en` since `E` is Euler's number (C0 row). Van Hove
content: `Series` of $E(k)$ at $k=0$ and $k=1$ is quadratic in $k-k_{\mathrm{edge}}$, so
$g\sim1/\sqrt{\lvert E-E_{\mathrm{edge}}\rvert}$; exhibit the $-\tfrac12$ slope of $\log g$
against $\log\lvert E-E_{\mathrm{edge}}\rvert$ near both edges. Verification by state counting
reusing the dispersion, able to refute a wrong $g$: sample $E(k)$ on a uniform $N$-point
$k$-grid, bin with `BinCounts`, and compare bin fractions to $\int g\,dE$ over the same bins
(`NIntegrate`, safe here, no point measure), plus the sum-rule total $\int g\,dE=1$. Close on
the contrast between the divergent edges and the flat band interior: the $1/\sqrt{}$ spikes are
what a tunneling or absorption measurement actually resolves.

#### 9.5 [MSc] How do I take the tight-binding limit and build Wannier functions from Bloch states?

Pinned example: the deep cosine lattice $V_0\cos2x$ at $V_0=10$ (cold-atom regime), lowest band.
Assemble Bloch states $\psi_k$ on a uniform $k$-grid from the `MathieuC`/`MathieuS` Floquet pair
at `MathieuCharacteristicA` (C1-probed evaluations), and fix the gauge per $k$, say $u_k(0)$
real and positive, before summing; without a smooth gauge the Wannier sum delocalizes, which is
the one nontrivial choice in the construction. Build
$w_0(x)=\frac{a}{2\pi}\int_{\mathrm{BZ}}\psi_k(x)\,dk$ as a uniform $k$-sum (`Table`, `Total`;
the integrand is periodic in $k$ so the uniform sum converges fast) and $w_1(x)=w_0(x-a)$ by
translation. Extract $J$ two independent ways: from the band,
$J=\big(E(\pi/a)-E(0)\big)/4$ off `MathieuCharacteristicA`/`MathieuCharacteristicB`, with the
residual of $E(k)$ against $\varepsilon-2J\cos ka$ quantifying the neglected next-nearest
hopping; and from the states, $J=-\int w_0\,\hat H\,w_1\,dx$ with
$\hat Hw=-\tfrac12w''+V_0\cos(2x)\,w$ via `NIntegrate` (delta-free, safe). Agreement of the two
$J$ values is the refuting check; disagreement convicts the gauge or the $k$-sum. Show
exponential localization as linear $\log\lvert w_0\rvert$ tails steepening across
$V_0\in\{5,10,15\}$. Close with the experimental knob: $J$ falls roughly exponentially in
$\sqrt{V_0}$, which is how a lattice-depth ramp freezes tunneling in cold-atom experiments.

#### 9.6 [MSc] How do I exhibit Bloch oscillations (a Wannier-Stark ladder) under a static field?

Pinned example: the tilted deep lattice in its tight-binding frame,
$H=-J\sum_{n}(\vert n{+}1\rangle\langle n\vert+\mathrm{h.c.})+Fa\sum_{n}n\,\vert n\rangle\langle n\vert$
on $N=101$ sites via `SparseArray` with `Band`, pinned at $J=1$, $Fa=0.4$ so $2J/(Fa)=5$ spans
several sites; the ladder is primary because the C3 verdict flags the continuum run as an open
conflict (the unbounded tilt $Fx$ is incompatible with the periodic pseudospectral grid as
written, and the long-time cost is unprobed). Sorted `Eigenvalues` show constant bulk spacing
$Fa$ to machine precision, the Wannier-Stark ladder itself, with chain-edge deviations as the
finite-size failure edge; the ladder amplitudes are `BesselJ`$[n-m,2J/(Fa)]$ (verify at
authoring). Evolve a broad $k_0=0$ packet with the `MatrixExp` action form
`MatrixExp[-I H t, v]`: $\langle x\rangle(t)=a\sum_{n}n\lvert\psi_n\rvert^2$ oscillates at
$\omega_B=Fa$ with peak-to-peak excursion $4J/F$, the semiclassical
$x(t)=x_0+\frac{2J}{F}\big(\cos(Fat)-1\big)$, and no DC transport, the counterintuitive point.
Refuting checks reusing the evolution: the equally spaced spectrum forces exact periodicity,
$1-\lvert\langle\psi(0)\vert\psi(T_B)\rangle\rvert$ at machine zero for $T_B=2\pi/(Fa)$, and the
measured excursion must equal $4J/F$. Flagged cross-check: one Bloch period of the continuum C3
recipe (`NDSolveValue`, tight goals 10) on the $V_0=10$ lattice tilted by $Fx$ in a large box
with explicit Dirichlet walls sized several times $4J/F$ (a time-dependent gauge sweeping
$k\to k-Ft$ would restore periodicity but is unprobed); norm conservation cannot detect wall
reflection (walls reflect at norm 1), so certify by doubling the box and demanding
$\langle x\rangle(t)$ unchanged, with $dx$-weighted grid sums, never `NIntegrate` on the
oscillatory interpolant. Keep $Fa$ well below the band gap or Landau-Zener leakage breaks the
single-band picture. Close with the experiment: cold atoms in a vertical optical lattice read
the local gravitational force off $\omega_B=Fa$.

## Part 10 Plan: Two and three dimensions: separation of variables

Questions: 7 (10.1 through 10.7). Class census: C7 for 10.1, 10.3, 10.5, 10.7 and C0 for 10.2,
10.4, 10.6, per the `Route-Table.md` class census (C7 and C0 rows) and the full C7 verdict block.

### Common ground

A product ansatz turns one eigenvalue PDE into ordinary ones: substituting $\psi=X(x)Y(y)$ into
$-\tfrac12\nabla^2\psi+V\psi=E\psi$ and dividing by $XY$ leaves additive pieces each depending on a
single coordinate, so each is constant; eigenvalues add, $E_{n_x n_y}=E_{n_x}+E_{n_y}$,
eigenfunctions multiply, and degeneracy is the number of index tuples sharing one sum. The two-body
kinetic energy separates in center-of-mass and relative coordinates with total mass $M=m_1+m_2$ and
reduced mass $\mu=m_1 m_2/(m_1+m_2)$; a central potential separates as
$\psi=R(r)\,Y_{lm}(\theta,\phi)$, and $u=rR$ turns the radial problem into the one-dimensional
$-\tfrac12 u''+\bigl(V(r)+\tfrac{l(l+1)}{2r^2}\bigr)u=Eu$. Natural units $\hbar=m=1$ throughout.
Standing C7 inheritances (`Route-Table.md`, C7 verdict): the derivation is the substitute-and-divide
computation, Expand on the divided equation plus a per-term FreeQ certificate of the split;
DEigensystem confirms exact spectra on exact regions but refuses symbolic side lengths, so symbolic
spectra come from the derivation alone; the independent cross-check is full-dimensional FEM
NDEigensystem with explicit DirichletCondition (omitting BCs is silent Neumann, probed in 2D) and a
mesh sweep; degenerate shells are compared by projector or multiplicity, never by individual
eigenfunctions (the FEM basis inside a shell is a genuinely mixed unitary rotation, probed); and
NDEigensystem eigenfunctions are InterpolatingFunction expressions whose pointwise evaluation nests
silently unevaluated and looks like a hang, so evaluation goes through the Head or substitution
idiom.

### Per-question entries

#### 10.1 [BSc] How do I separate variables in a two- or three-dimensional box and form the product eigenfunctions?

Work the pinned rectangle with $L_x/L_y=\sqrt2$ (a quantum dot; take $L_x=\sqrt2$, $L_y=1$), where
the incommensurate ratio isolates separability from degeneracy accidents. Substitute $X(x)Y(y)$,
Expand the divided equation, certify the split with the per-term FreeQ check, and assemble the 1D
box spectra into $E_{n_x n_y}=\tfrac{\pi^2}{2}\bigl(n_x^2/L_x^2+n_y^2/L_y^2\bigr)$ with $L_x,L_y$
symbolic, from the derivation alone, since DEigensystem refuses symbolic side lengths; run the same
mechanics once in 3D, where the certificate is probed too. Then confirm exactly: DEigensystem with
DirichletCondition on the exact region returns the first six levels
$\{\tfrac{3\pi^2}{4},\tfrac{3\pi^2}{2},\tfrac{9\pi^2}{4},\tfrac{11\pi^2}{4},3\pi^2,\tfrac{17\pi^2}{4}\}$
with product sines (probe 2), identical to the enumerated $E_{n_x n_y}$. Cross-check with FEM
NDEigenvalues on ToElementMesh over the Rectangle, explicit DirichletCondition (omit it and the
$\tfrac{3\pi^2}{4}$ ground state silently vanishes into the Neumann spectrum, probe 4), and a
MaxCellMeasure sweep (monotone, $2.4\times10^{-3}\to1.5\times10^{-5}$, probe 3); the refuting
comparison is the interleaved level ordering $(1,1),(2,1),(1,2),(3,1),(2,2)$, which any wrong
quantum-number assignment scrambles. Close on the contrast with the commensurate square, where
swapped index pairs collide and the one-to-one level labeling is lost.

#### 10.2 [BSc] How do I find and explain the degeneracies of a square or cubic box?

The pinned exhibit is the square-box shell $s=n_x^2+n_y^2=50=1^2+7^2=5^2+5^2$, three states
$(1,7),(7,1),(5,5)$ at $E=\tfrac{\pi^2}{2L^2}\,50$. Enumerate levels with Tuples and Select up to a
cap $s\le s_{\max}$ and GatherBy the exact integer $s$ (integer arithmetic, no floating
comparison), reading multiplicities off with Tally; count the same shells number-theoretically with
SquaresR[2,s] and PowersRepresentations[s,2,2], reducing the signed lattice count to ordered
positive pairs (drop zero-component representations, divide by 4 for signs: $r_2(50)=12$ collapses
to 3). The refuting check runs both routes over every shell up to the cap and demands equality
level by level, which any mishandled zero or equal-component case breaks; the same Tuples
enumeration and SquaresR[3,s] extend the count to the cubic box. Explain the split: $(1,7)$ against
$(7,1)$ is the square's exchange symmetry, while their coincidence with $(5,5)$ is arithmetic,
forced by the number 50 and not by geometry. Close on the contrast with an incommensurate
rectangle, where every such collision disappears: degeneracy lives in the spectrum's arithmetic,
not in dimensionality.

#### 10.3 [BSc] How do I solve the two-dimensional isotropic oscillator and count its degeneracies?

Work the pinned $N=2$ shell at $E=3$: the Cartesian triplet $\psi_{2,0},\psi_{1,1},\psi_{0,2}$
built from HermiteH 1D factors, against the polar states $(n_r,m)=(0,\pm2),(1,0)$ built as
$e^{im\phi}\,r^{|m|}L_{n_r}^{|m|}(r^2)\,e^{-r^2/2}$ with LaguerreL. Separate in Cartesian
coordinates by substitute-and-divide with the FreeQ certificate, read off
$E_{n_x n_y}=n_x+n_y+1$, then compute the explicit $3\times3$ unitary between the two bases by
Integrate in polar coordinates ($x=r\cos\phi$, $y=r\sin\phi$), FullSimplify certifying
$UU^\dagger=\mathbb 1$ exactly: a wrong normalization or wrong polar radial factor refutes it, and
so does a nonzero residual of any polar state in the reused eigenvalue equation. Cross-check with
FEM NDEigensystem on a box (explicit DirichletCondition, mesh sweep, plus a box-size sweep, which
the probes held fixed): compare the shell multiplicities $\{1,2,3\}$ and the $N=2$ shell
projector, never individual eigenfunctions, since the probed FEM shell basis is a genuinely mixed
rotation ($|M_{ij}|_{\max}=0.958$ while $MM^{T}=\mathbb 1$ to $10^{-6}$), and evaluate
eigenfunctions only through the Head or substitution idiom. Close by widening the count: solutions
of $2n_r+|m|=N$ number exactly $N+1$, the same ladder the Cartesian pairs give.

#### 10.4 [BSc] How do I separate the two-body problem into center-of-mass and relative motion, reducing it to a one-body problem with the reduced mass $\mu$?

Apply the two-body kinetic operator $-\tfrac{1}{2m_1}\nabla_1^2-\tfrac{1}{2m_2}\nabla_2^2$ to the
composed $\Psi\bigl(\tfrac{m_1\vec r_1+m_2\vec r_2}{m_1+m_2},\,\vec r_1-\vec r_2\bigr)$ with D (the
chain rule is automatic) and Simplify with $m_1,m_2$ symbolic: the identity against
$-\tfrac{1}{2M}\nabla_R^2-\tfrac{1}{2\mu}\nabla_r^2$, $M=m_1+m_2$, $\mu=m_1m_2/M$, must simplify
to zero, the refuting check, since any wrong mass combination leaves a surviving cross term.
Confirm $\lim_{m_2\to\infty}\mu=m_1$ with Limit, and split off the free center of mass with the
FreeQ certificate on the product $\Phi(R)\phi(r)$. Then cash into the pinned spectroscopy: the
Rydberg scales as $R_M=R_\infty/(1+m_e/M)$, so the hydrogen and deuterium Balmer $n=3\to2$ lines
near $656\,\mathrm{nm}$ differ by
$\Delta\lambda\approx\lambda_\infty(m_e/m_p-m_e/m_d)\approx0.18\,\mathrm{nm}$; compute the mass
ratios with Quantity and UnitConvert on the "ElectronMass", "ProtonMass", "DeuteronMass" units
(verify at authoring; typed CODATA ratios otherwise) and compare with the measured
$0.179\,\mathrm{nm}$ isotope shift. Close there: the reduced mass is a measured line splitting, the
one that revealed deuterium.

#### 10.5 [BSc] How do I separate the Schrodinger equation in spherical coordinates into a radial and an angular equation?

Keep the potential a generic symbolic $V(r)$: the point is that the $l(l+1)$ term emerges for every
central potential, not for one example. Write the operator with Laplacian in the "Spherical" chart,
substitute $R(r)\,Y(\theta,\phi)$, multiply by $2r^2/(RY)$, Expand, and certify the split per term
with FreeQ (the 3D mechanics probed in the C7 verdict); name the separation constant $\lambda$ and
reduce the radial side with $u=rR$ to $-\tfrac12u''+\bigl(V(r)+\tfrac{\lambda}{2r^2}\bigr)u=Eu$.
Pin $\lambda=l(l+1)$ by substituting SphericalHarmonicY over the whole family $l=0,\dots,3$, all
$m$, with Table and FullSimplify: every ratio must return exactly $l(l+1)$, a family check that
refutes the sign or $1/\sin^2\theta$ slip a single $(l,m)$ could mask. Then reassemble: impose the
radial equation on symbolic $R$ and confirm the full 3D residual simplifies to zero with $V$ still
generic, reusing the definitions. The full-3D FEM cross-check on a Ball is unprobed for this class;
if attempted at authoring it is gated (explicit DirichletCondition, mesh sweep), never assumed.
Close on the universality: the angular problem is solved once for all central potentials, and
everything $V$-specific lives in one radial dimension.

#### 10.6 [BSc] How do I build the effective radial potential and read off the centrifugal barrier?

Build $V_{\mathrm{eff}}(r)=-\tfrac1r+\tfrac{l(l+1)}{2r^2}$ for the pinned hydrogen family
$l=0,1,2$ with Table, and extract its geometry symbolically with D and Solve: the zero crossing at
$r=l(l+1)/2$, the well minimum $-\tfrac{1}{2l(l+1)}$ at $r=l(l+1)$, and the classical turning
points $r_{\pm}=n^2\bigl(1\pm\sqrt{1-l(l+1)/n^2}\bigr)$ at $E_n=-\tfrac{1}{2n^2}$; Plot the three
curves with the $E_n$ lines. The refuting check reuses the turning points: Reduce must show they
are real exactly when $l\le n-1$, so the classical geometry of $V_{\mathrm{eff}}$ reproduces the
hydrogen selection rule $l_{\max}=n-1$, and a wrong centrifugal coefficient (an $l^2$, a missing
factor 2) breaks the count. Contrast $l=0$, which has no wall at all, with the inner turning point
receding from the origin as $l$ grows. Close in the lab: only $s$ states keep finite density at the
nucleus, which is why contact physics (hyperfine structure, electron capture) selects $l=0$.
Concern: the attractive-Coulomb $V_{\mathrm{eff}}$ has no local maximum, so the ledger's "barrier
heights" name no finite number; the entry reads the barrier quantitatively as the zero crossing,
well depth, and inner turning point instead.

#### 10.7 [MSc] How do I separate variables in parabolic coordinates (the natural frame for the Stark and Coulomb problems)?

Work the pinned hydrogen problem in $\xi=r+z$, $\eta=r-z$, $\phi$, where
$\nabla^2=\tfrac{4}{\xi+\eta}\bigl[\partial_\xi(\xi\,\partial_\xi)+\partial_\eta(\eta\,\partial_\eta)\bigr]+\tfrac{1}{\xi\eta}\,\partial_\phi^2$
and $-1/r=-2/(\xi+\eta)$; derive that form by explicit change of variables with D rather than a
chart name (CoordinateChartData's "Parabolic" convention differs by squares, verify at authoring if
used). Substitute $f_1(\xi)f_2(\eta)e^{im\phi}$, divide by the product, multiply by $(\xi+\eta)$
(the $m^2$ term splits via $\tfrac{1}{\xi\eta}=\tfrac{1}{\xi+\eta}(\tfrac1\xi+\tfrac1\eta)$),
Expand, and certify with per-term FreeQ; this is the C7 verdict's flagged weakest fit, because the
separation constants couple through the Coulomb term as $\beta_1+\beta_2=1$ instead of splitting
freely. Read the termination condition of each confluent-hypergeometric 1D equation,
$\beta_i=\sqrt{-2E}\bigl(n_i+\tfrac{|m|+1}{2}\bigr)$, and Solve $\beta_1+\beta_2=1$ for
$E=-\tfrac{1}{2n^2}$ with $n=n_1+n_2+|m|+1$. Verification stays symbolic, per the verdict: the
reassembly residual (impose both separated ODEs and simplify the full PDE to zero), and the
recount, enumerating $(n_1,n_2,m)$ at fixed $n$ with Tuples and Select to exactly $n^2$ states,
matching the known spherical degeneracy; the recount refutes any misplaced $\tfrac{|m|+1}{2}$. A
numeric cross-check is authoring-gated: the spectrum is negative, so the magnitude-ordering trap is
load-bearing (shift mandatory) and FEM near $r=0$ is unprobed. Close on why the frame exists: a
uniform field adds $Fz=F(\xi-\eta)/2$ and the equation stays separable in exactly these
coordinates, the groundwork for the Stark effect.

## Part 11 Plan: Orbital angular momentum in continuous space

Six questions. Class census (per the class census table in `Route-Table.md`): C0 for all of 11.1
through 11.6, differential-operator algebra and quadrature on the sphere, no solver anywhere.

### Common ground

The part rests on one operator family and its representation theory. Orbital angular momentum is
$\vec L=\vec r\times\vec p$ with $\vec p=-i\nabla$, which in spherical coordinates
($x=r\sin\theta\cos\phi$, $y=r\sin\theta\sin\phi$, $z=r\cos\theta$) closes into first-order
differential operators in the angles alone; the algebra $[L_i,L_j]=i\epsilon_{ijk}L_k$ with
Casimir $L^2=L_x^2+L_y^2+L_z^2$ and ladders $L_\pm=L_x\pm i\,L_y$; the sphere inner product
$\langle f,g\rangle=\int_0^{2\pi}\!\!\int_0^\pi\bar f\,g\,\sin\theta\,d\theta\,d\phi$ under which
the joint eigenfunctions obey $L^2Y_{lm}=l(l+1)\,Y_{lm}$ and $L_zY_{lm}=m\,Y_{lm}$; rotations act
on a multiplet by $R\,Y_{lm}=\sum_{m'}D^l_{m'm}\,Y_{lm'}$ without leaving it; and products couple
by Clebsch-Gordan coefficients $C^{LM}_{l_1m_1,l_2m_2}$. Every entry runs on D, Integrate,
Simplify and FullSimplify under the assumptions $r>0$, $0<\theta<\pi$, $\phi$ real (the C0
assumption discipline), with the built-ins SphericalHarmonicY, WignerD, and ClebschGordan serving
as anchors that constructions are checked against, never as sources the derivations quote.

#### 11.1 [BSc] How do I write the orbital angular-momentum operators $L_x,L_y,L_z$ as differential operators on the angles?

Earn all three forms by explicit change of variables on the pinned generic $f(\theta,\phi)$, never
by quotation (C0 per `Route-Table.md`): form the composite $f(\theta(x,y,z),\phi(x,y,z))$ with
$\theta=\arccos\bigl(z/\sqrt{x^2+y^2+z^2}\bigr)$ and the two-argument ArcTan for $\phi$, let D
chain-rule it through $-i(x\,\partial_y-y\,\partial_x)$ and the two cyclic partners of
$\vec r\times\vec p$, then substitute $x=r\sin\theta\cos\phi$, $y=r\sin\theta\sin\phi$,
$z=r\cos\theta$ and Simplify under $r>0$, $0<\theta<\pi$, $\phi$ real, reading off
$L_z=-i\,\partial_\phi$, $L_x=i(\sin\phi\,\partial_\theta+\cot\theta\cos\phi\,\partial_\phi)$,
$L_y=i(-\cos\phi\,\partial_\theta+\cot\theta\sin\phi\,\partial_\phi)$ with every power of $r$
cancelling in the output. The read-off from the simplified generic result is where a transcription
slip hides, so verify refutably there: apply the Cartesian and the angular representation of each
component to the concrete carrier $(x+iy)^2z/r^3$ (a $Y_{32}$ shape, alive under all three
operators) and Simplify the three differences to zero, reusing both definitions. Note the
$\cot\theta$ poles at $\theta=0,\pi$ as chart artifacts of operators that are regular in Cartesian
form: the edge regime of the coordinate change. Close on what cancelled: no $L_i$ contains $r$ or
$\partial_r$, so orbital angular momentum generates pure rotations and the rest of the part lives
on the unit sphere.

#### 11.2 [BSc] How do I verify $[L_i,L_j]=i\epsilon_{ijk}L_k$ and that $L^2$ is the Casimir with eigenvalue $l(l+1)$?

Redefine the three angular operators inside the answer as pure functions acting on expressions in
$(\theta,\phi)$ and run the algebra on the pinned generic $f(\theta,\phi)$ (C0 per
`Route-Table.md`): Simplify $[L_x,L_y]f-i\,L_zf$ and its two cyclic partners to zero in a single
Table over the index triples, so the whole $\epsilon_{ijk}$ structure is exercised rather than one
instance; assemble $L^2f=L_x(L_xf)+L_y(L_yf)+L_z(L_zf)$ and Simplify $[L^2,L_i]f$ to zero for all
three $i$, the Casimir property proved on a generic function rather than a lucky state. Then the
eigenvalue on the pinned concrete family: over a Table of SphericalHarmonicY spanning the full
$l=2$ multiplet $m=-2,\dots,2$, FullSimplify $L^2Y_{2m}-6\,Y_{2m}$ to zero for all five members,
so $l(l+1)=6$ is earned on the whole block. Discriminate rather than illustrate: on the
cross-multiplet mixture $Y_{11}+Y_{22}$ the same $L^2$ fails to act as a scalar (the ratio
$L^2\psi/\psi$ does not Simplify to a constant), which refutes the reading of $l(l+1)$ as a
property of the operator instead of the multiplet. Close: commutators plus Casimir spectrum are
the complete rotation fingerprint that every later construction in the part must carry.

#### 11.3 [BSc] How do I build the spherical harmonics $Y_{lm}$ as the simultaneous eigenfunctions of $L^2$ and $L_z$ and check their orthonormality?

Construct the pinned $l=2$ family from its highest weight with $L_z$ and $L_\pm=L_x\pm i\,L_y$
redefined in their angular forms (C0 per `Route-Table.md`): $L_zY_{22}=2\,Y_{22}$ is a first-order
equation in $\phi$ whose DSolve solution forces $Y_{22}=g(\theta)\,e^{2i\phi}$, and $L_+Y_{22}=0$
then collapses to the one-variable ODE $g'(\theta)=2\cot\theta\,g(\theta)$, DSolve giving
$g\propto\sin^2\theta$; normalize by Integrate of ComplexExpand of $Y\bar Y$ against the sphere
measure $\sin\theta\,d\theta\,d\phi$. Ladder down with NestList of $L_-$, normalizing each rung
the same way, and confirm a fifth application annihilates the bottom state, the multiplet closing
from within. The refuting check reuses the constructed family: the $5\times5$ Gram matrix by
Outer and Integrate must equal IdentityMatrix, and any wrong rung breaks its row. Only then
compare with the built-in: each ratio of constructed $Y_{2m}$ to SphericalHarmonicY must Simplify
to a $(\theta,\phi)$-independent unimodular constant, the comparison holding up to a fixed phase
per $m$ because the ladder fixes relative phases while Condon-Shortley fixes the built-in's
absolute ones (verify at authoring which $m$ land on phase exactly $1$). Close: one first-order
PDE plus four ladder steps produce the entire multiplet; the built-in is a lookup table for what
the algebra just constructed.

#### 11.4 [MSc] How do the raising and lowering operators $L_\pm$ act on the $Y_{lm}$?

Take the built-in $l=3$ family as the pinned carrier, ladder operators redefined in their angular
forms (C0 per `Route-Table.md`): over a Table in $m=-3,\dots,3$, FullSimplify
$L_+Y_{3m}-\sqrt{12-m(m+1)}\,Y_{3,m+1}$ and $L_-Y_{3m}-\sqrt{12-m(m-1)}\,Y_{3,m-1}$ to zero, the
whole multiplet at once with SphericalHarmonicY supplying every carrier, and the edge
annihilations $L_+Y_{33}=0$, $L_-Y_{3,-3}=0$ computed as the exact zeros of the factor where
$m(m\pm1)=l(l+1)$. A single-rung identity can hide a convention error, so triangulate magnitude
and phase separately, reusing the definitions: the two-step diagonal
$L_-L_+Y_{3m}=(12-m(m+1))\,Y_{3m}$ is phase-insensitive, and the norm
$\int\vert L_+Y_{3m}\vert^2\,d\Omega=12-m(m+1)$ by Integrate and ComplexExpand checks magnitude
alone; the positive real $\sqrt{l(l+1)-m(m\pm1)}$ one-step factors then hold exactly in the
Condon-Shortley phase the built-in carries, and an $m$-dependent sign in the one-step residuals is
precisely what a convention mismatch would look like. Close on termination: the factor vanishes
only at the edges, so the ladder stops itself after $2l+1=7$ states; the multiplet dimension is a
computed consequence, not an input.

#### 11.5 [MSc] How do I rotate a wavefunction and represent the rotation on the $Y_{lm}$ by the Wigner $D$-matrix?

Rotate the pinned $l=2$ multiplet through its argument (C0 per `Route-Table.md`): each $Y_{2m}$ is
a homogeneous quadratic in the unit vector, so substitute the inverse rotation of $(x,y,z)$ built
from RotationMatrix in $zyz$ Euler order with $\beta$ symbolic about the $y$ axis; the flanking
$z$ rotations act as computed $\phi$ shifts, pure $e^{-im\alpha}$ and $e^{-im'\gamma}$ phases, so
the middle $\beta$ block is the honest content. Rewrite the rotated functions in angles and
project back on the basis: coefficients $c_{m'm}(\beta)$ by Integrate of $\bar Y_{2m'}$ times the
rotated $Y_{2m}$ over the sphere measure, symbolic $\beta$ surviving the quadrature. The
$5\times5$ coefficient matrix must match the WignerD entries $D^2_{m'm}$: FullSimplify each
entrywise difference to zero at symbolic $\beta$; WignerD's Euler-angle ordering, its
active-versus-passive sense, and its conjugation convention vary across references, so fix the
dictionary on one nonvanishing entry and hold it for the whole matrix (verify at authoring).
Refuting checks that reuse the machinery: the matrix is unitary at symbolic $\beta$ (ComplexExpand
under $\beta$ real), collapses to IdentityMatrix at $\beta=0$, and leaves the addition-theorem
invariant $\sum_{m}\vert Y_{2m}\vert^2=5/(4\pi)$ untouched. Close: no $l=1$ or $l=3$ component
ever appears in the expansion; the five-dimensional block is rotationally closed, which is what
irreducible means in practice.

#### 11.6 [MSc] How do I couple two angular momenta and compute the Clebsch-Gordan coefficients (with the three-$Y_{lm}$ Gaunt integral as the orbital instance), the change of basis reused later for adding spin?

Compute the pinned $1\otimes1\to0,1,2$ decomposition twice on independent machinery (C0 per
`Route-Table.md`). First, ClebschGordan tabulated with Table over the physical grid only,
$M=m_1+m_2$: off-shell argument triples emit a message while returning zero, so never generate
them (the no-Quiet rule). Second, the Gaunt integrals
$\int Y_{1m_1}Y_{1m_2}\overline{Y_{LM}}\,d\Omega$ by Integrate over the sphere measure with
SphericalHarmonicY and Conjugate, compared entry by entry against the closed Gaunt form
$\sqrt{9/(4\pi(2L+1))}\,C^{L0}_{10,10}\,C^{LM}_{1m_1,1m_2}$ assembled from the same ClebschGordan
values: FullSimplify each difference to zero over all $(L,M,m_1,m_2)$. The $L=0,2$ channels pin
the coefficients up to the computed $m$-independent reduced factor; the $L=1$ channel is the edge,
both sides vanishing identically since $C^{10}_{10,10}=0$ by parity ($1+1+1$ odd), equivalently
the $L=1$ combination is antisymmetric under $m_1\leftrightarrow m_2$ while a one-sphere product
$Y_{1m_1}Y_{1m_2}$ is symmetric. Certify the odd channel, and every other, with the part's own
operators on two spheres: assemble
$\Psi_{LM}=\sum C^{LM}_{1m_1,1m_2}\,Y_{1m_1}(\theta_1,\phi_1)\,Y_{1m_2}(\theta_2,\phi_2)$ and
verify $L_z^{\mathrm{tot}}\Psi_{LM}=M\,\Psi_{LM}$ and
$(L^{\mathrm{tot}})^2\Psi_{LM}=L(L+1)\,\Psi_{LM}$ for all nine states, the total operators built
from the per-sphere angular forms; one wrong coefficient in any channel breaks it. Close in the
lab: the same Gaunt zero is the parity selection rule forbidding electric-dipole transitions
between two $p$ states, and the $9=1+3+5$ change of basis is exactly the recoupling matrix reused
when spin is added.

## Part 12 Plan: Central potentials and the hydrogen atom

Questions: 9 (12.1 through 12.9). Class census, per the class census table in `Route-Table.md`:
C1: 12.1, 12.2; C0: 12.3, 12.6; C2: 12.4, 12.7, 12.8; C7: 12.5; C9: 12.9. Both probed class
representatives live in this part (12.1 for C1, 12.7 for C2), so their entries cite gates verbatim;
the Route-Table C0 kernel facts bind throughout (`E` is Euler's number, use `en` for energies;
`Simplify` closes residuals and norms only under explicit integer, positivity, and reality
assumptions).

### Common ground

Rotational invariance reduces every question here to one radial dimension:
$\psi_{nlm}=\tfrac{u_{nl}(r)}{r}\,Y_{lm}(\theta,\phi)$, with the reduced radial equation
$-\tfrac12 u_{nl}''+\left(V(r)+\tfrac{l(l+1)}{2r^2}\right)u_{nl}=E\,u_{nl}$ on $(0,\infty)$ and
$u_{nl}(0)=0$; the centrifugal term defines the effective potential
$V_{\mathrm{eff}}(r)=V(r)+l(l+1)/(2r^2)$. Normalization splits, $\int_0^\infty u_{nl}^2\,dr=1$ with
the angular unit carried by $\int|Y_{lm}|^2\,d\Omega=1$, and the radial distribution is
$P(r)=u_{nl}^2=r^2R_{nl}^2$. For hydrogen, $V=-1/r$ in natural units $\hbar=m=e=1$ (Bohr radius
and Hartree both 1), the spectrum is $E_n=-1/(2n^2)$ with degeneracy
$\sum_{l=0}^{n-1}(2l+1)=n^2$: larger than rotational symmetry alone explains, which is the thread
running from the exact solution through the numerics to the dynamical symmetry and the
semiclassical limit.

### Per-question entries

#### 12.1 [BSc] How do I solve the Coulomb radial equation and obtain the bound-state energies $E_n=-1/(2n^2)$?

C1 representative, gates cited verbatim from the Route-Table C1 verdict. DSolve the reduced radial
equation $-\tfrac12 u''+(\tfrac{l(l+1)}{2r^2}-\tfrac1r)u=E\,u$ with symbolic $E$ and $l$: the
general solution arrives as WhittakerM and WhittakerW in under a second; never hand DSolve both
decay boundary conditions, it hangs with no message. Quantization is the manual termination
read-off: the normalizable branch terminates exactly when $1/\sqrt{-2E}=n$ with $n$ a positive
integer, and Solve turns that into $E_n=-1/(2n^2)$, tabulated for $n=1..4$. Convert the complex
Whittaker form to the real Laguerre form through the explicit branch substitution
$en\to-\kappa^2/2$ with $\kappa>0$, made per level rather than left to automatic simplification
(the C1 verdict flags this bookkeeping as an open risk), assemble
$u_{nl}=r^{l+1}e^{-r/n}L_{n-l-1}^{2l+1}(2r/n)$ with LaguerreL, close the residual to 0 with
FullSimplify under the integer and positivity assumptions ($n,l$ integers, $n\ge l+1$, $l\ge0$,
$r>0$), and take exact Integrate norms. Cross-check on independent machinery: NDEigenvalues with
explicit DirichletCondition at both ends and Arnoldi Shift below the ground state, since the
unshifted call silently returns only near-zero values and misses every bound level; take the
domain literally at $r=0$, which C2 probe p2 measured as best, because an origin cutoff installs
an error floor $\Delta E\approx2\varepsilon$ that mesh refinement cannot pass (C1 probe 7:
$1.8\times10^{-2}$, $2.0\times10^{-4}$, $2.0\times10^{-6}$ at $\varepsilon=10^{-2},10^{-4},10^{-6}$)
and is actively harmful once $l\ge1$. Close with the Rydberg accumulation $E_n\to0^-$: infinitely
many levels crowd under the continuum, the structural fact any truncated numeric box must answer
for.

#### 12.2 [BSc] How do I build the radial wavefunctions from the associated Laguerre polynomials and normalize them?

Carrier $R_{31}$: $l=1$ with exactly one radial node ($n-l-1=1$). State the general form
$R_{nl}(r)=\sqrt{(2/n)^3\,\tfrac{(n-l-1)!}{2n\,(n+l)!}}\,(2r/n)^{l}e^{-r/n}L_{n-l-1}^{2l+1}(2r/n)$
in prose first, then build it with LaguerreL. The symbolic-$(n,l)$ closed-form norm is an open
risk in the C1 verdict while fixed-$(n,l)$ exact norms are certified: attempt
$\int_0^\infty R_{nl}^2\,r^2\,dr=1$ with Integrate under the integer and positivity assumptions,
and whatever the general integral returns, certify the normalization at fixed pairs, at least
$(3,1)$, $(2,0)$, $(4,2)$, with an NIntegrate spot-check of the $(3,1)$ norm reusing the same
defined $R_{31}$ as the independent numeric reference. Verification that can refute: the
radial-equation residual of the defined $R_{31}$ FullSimplifies to 0 under the same assumptions,
and the orthogonality $\int_0^\infty R_{21}R_{31}\,r^2\,dr=0$ reuses both definitions against the
same measure. The limiting regime is the origin: Series or Limit on $u_{nl}=rR_{nl}$ as $r\to0$
must give the centrifugal suppression $u_{nl}\sim r^{l+1}$, so the $l=1$ carrier vanishes
quadratically where an $s$ state vanishes only linearly. Close by reading the node count $n-l-1$
off the Laguerre degree and locating the $R_{31}$ node exactly at $r=6$ from $L_1^{3}(2r/3)=0$.

#### 12.3 [BSc] How do I assemble the full $\psi_{nlm}$, plot probability densities and the radial distribution, and compute $\langle r\rangle$?

Assemble $\psi_{321}=R_{32}(r)\,Y_2^1(\theta,\phi)$, a $3d$, $m=1$ state, from LaguerreL and
SphericalHarmonicY, defining $R_{32}$ inside the entry. Normalization splits: the spherical
harmonic is orthonormal on the sphere, so $\int|\psi_{321}|^2\,d^3r$ reduces to an exact radial
Integrate; the density $|\psi_{321}|^2$ is $\phi$-independent because the $e^{i\phi}$ phase drops
under ComplexExpand, so DensityPlot on the $(x,z)$ half-plane plus Plot of the radial distribution
$P(r)=r^2R_{32}^2$ carry the visualization (a 3D cut via SliceDensityPlot3D is optional and tagged
(verify at authoring)). Compute $\langle r\rangle$ exactly with Integrate and set it against the
closed formula $\langle r\rangle_{nl}=\tfrac12\left(3n^2-l(l+1)\right)=\tfrac{21}{2}$; an
independent NIntegrate spot-check reusing the defined $P(r)$ can refute both. Close with the
circular-orbit reading: for $l=n-1$ the most probable radius, from Solve on $P'(r)=0$, sits
exactly at the Bohr value $r=n^2=9$, while $\langle r\rangle=\tfrac{21}{2}$ sits above it because
$P(r)$ is skewed outward.

#### 12.4 [BSc] How do I find the bound states of a spherical square well numerically?

Spherical finite well $V=-V_0$ for $r<a$, zero beyond, at $V_0=10$, $a=1$: deep enough to bind one
level each for $l=0$ and $l=1$, since the first $s$ level appears at $\sqrt{2V_0}\,a=\pi/2$ and
the first $p$ level at $\sqrt{2V_0}\,a=\pi$, while $\sqrt{20}\approx4.47$ clears both. Numeric
side per the C2 recipe:
NDEigenvalues on the reduced $u$ with explicit DirichletCondition at both ends, Arnoldi Shift
below $-V_0$ then Sort, MaxCellMeasure 0.01, the domain taken literally at $r=0$ (an $\epsilon$
cutoff is actively harmful for $l\ge1$), and the truncation radius fixed by an $R$-doubling sweep.
Benchmark side: exact matching, region by region because DSolve silently echoes a Piecewise
potential: with $k=\sqrt{2(E+V_0)}$ and $\kappa=\sqrt{-2E}$, the $l=0$ condition is
$k\cot ka=-\kappa$, and $l=1$ matches the log-derivative of $r\,j_1(kr)$ (SphericalBesselJ) to the
elementary decaying exterior $e^{-\kappa r}(1+1/(\kappa r))$, FindRoot solving each condition.
Certify the closed forms before comparing them: substitute each region's expression back into the
reduced radial equation and Simplify the residual to 0 under $k,\kappa,r>0$, so a mismatch cannot
be blamed on a mistyped Bessel form. Verification: the FEM levels against the FindRoot roots,
independent machinery on both sides, so any disagreement beyond the mesh floor refutes one route.
Close by sweeping $V_0$ downward through the two thresholds: at $\sqrt{2V_0}\,a=\pi$ the $l=1$
level is expelled, and only below $\sqrt{2V_0}\,a=\pi/2$ does the well hold no bound state at all,
the sharpest contrast with one dimension, where any attractive well binds.

#### 12.5 [MSc] How do I solve the three-dimensional isotropic oscillator and reconcile its Cartesian and spherical degeneracy counts?

The carrier is the $N=2$ shell of the 3D isotropic oscillator at $E=\tfrac72$, but the count is
earned before the shell is inspected. The separation is the route, per C7's primary: substitute
$\psi=X(x)Y(y)Z(z)$ into $\hat H\psi=E\psi$, divide by $\psi$, Expand, and certify the additive
split with a per-term FreeQ test (the 3D mechanics are probed in C7 probe1), leaving three Hermite
problems and $E=n_x+n_y+n_z+\tfrac32$; the spherical substitution $\psi=\tfrac{u(r)}{r}Y_{lm}$
reduces the same operator to the radial equation whose solutions
$R_{n_r l}\propto r^{l}e^{-r^2/2}L_{n_r}^{l+1/2}(r^2)$ give $E=2n_r+l+\tfrac32$. Equating the two
labelings turns the degeneracy into an identity rather than an observation: the Cartesian shell
$n_x+n_y+n_z=N$ holds $(N+1)(N+2)/2$ states, the spherical shell runs $l=N,N-2,\dots$ down to 0 or
1 with $n_r=(N-l)/2$, and $\sum_{l=N,N-2,\dots}(2l+1)=(N+1)(N+2)/2$ must hold for every $N$, closed
symbolically with Sum and exhibited at $N=2$ ($6=5+1$) and $N=3$ ($10=7+3$), the odd shell being
the one where a single $l$ cannot carry the count. The explicit map at $N=2$: rewrite each
spherical state as polynomial times Gaussian in Cartesian coordinates and compute the $6\times6$
overlap matrix by exact Integrate; verification that can refute: the matrix must be exactly
unitary and both bases must give residual $\hat H\psi=\tfrac72\psi$ under Simplify, either failure
exposing a wrong radial or angular factor. Degeneracy is reconciled by shell projector or
multiplicities, never by individual FEM eigenfunctions: inside a degenerate shell an eigensolver
returns an arbitrarily mixed basis
(largest single overlap 0.958 in the C7 verdict's 2D shell probe) while the shell projector is
reproduced. Any 3D NDEigensystem confirmation is authoring-gated, since 3D FEM cost is unprobed
in the C7 verdict: gate it with its own mesh sweep, or let the multiplicity count $\{5,1\}$ at the
clustered eigenvalue carry the numeric side. Close with the algebraic reading: the reconciliation
exhibits the $SU(3)$ degeneracy of the oscillator, the analogue of the Coulomb degeneracy that a
conserved Runge-Lenz vector explains.

#### 12.6 [MSc] How do I expose the extra Coulomb degeneracy through the conserved Runge-Lenz vector (the dynamical symmetry)?

Runge-Lenz $\hat{\vec A}=\tfrac12(\hat{\vec p}\times\hat{\vec L}-\hat{\vec L}\times\hat{\vec p})-\hat{\vec r}/r$,
with the symmetrization kept because $\hat{\vec p}\times\hat{\vec L}\ne-\hat{\vec L}\times\hat{\vec p}$
as operators. Realize everything as explicit differential operators on a generic test function:
$\hat p_j f=-i\,D[f,x_j]$, $\hat L$ composed from position and D, $r=\sqrt{x^2+y^2+z^2}$ written
out, then compute $[\hat H,\hat A_z]f=\hat H(\hat A_z f)-\hat A_z(\hat H f)$ and close it to
exactly 0 with Together and Simplify under $r>0$; this plain-D route needs nothing beyond core,
and any noncommutative-algebra package used instead is tagged (verify at authoring). The generic
$f(x,y,z)$ is what makes the check refuting: dropping the symmetrization or misordering $1/r$
against the derivatives leaves a visible nonzero remainder, so run the deliberately unsymmetrized
commutator beside the honest one as the contrast. Then the degeneracy: apply $\hat A_z$ to the
closed-form $\psi_{200}$ and Simplify the output to an explicit multiple of $\psi_{210}$: the
conserved vector rotates $2s$ into $2p$ inside the $n=2$ shell, which is exactly why states of
different $l$ share one energy. Close with the broken case: any screening deformation of $1/r$
destroys the commutator and with it the $l$ degeneracy, the mechanism behind alkali quantum
defects.

#### 12.7 [MSc] How do I compute central-potential eigenstates by applying `NDEigensystem` to the radial equation?

C2 representative, gates cited verbatim from the Route-Table C2 verdict. Hydrogen reduced radial
equation on $(0,R)$, benchmark $E_n=-1/(2n^2)$, with all four load-bearing elements of the recipe:
explicit `DirichletCondition[u[r]==0, True]` at both ends (omitting BCs silently yields a ground
level 24 percent wrong, no message); `Method -> {"Eigensystem" -> {"Arnoldi", "Shift" -> -0.6}}`
then Sort (unshifted, NDEigenvalues on this exact problem returns
$\{-0.0111, 0.0157, -0.0306, 0.0501, -0.0555, 0.0916\}$, a near-zero Rydberg tail plus box
continuum spanning $(-0.0555, 0.0916)$, the ground state absent with no message, C2 probe
p1-hydrogen-l0.wls); MaxCellMeasure 0.01; the domain taken literally at $r=0$, since an $\epsilon$
cutoff is actively harmful for $l\ge1$ ($r_{\min}=0$ reaches $5.5\times10^{-13}$ on $n=2$ where
$\epsilon=0.1$ costs $5.4\times10^{-5}$). Truncation is the teaching exhibit: failure invades from
the top as the Rydberg accumulation outgrows the box ($R=20$ puts $0.133$ on $n=5$; $R=160$ holds
all five levels at or below $9.4\times10^{-11}$), so $R$ is fixed by the doubling sweep; the
certified cross-check, a hand-built FD tridiagonal (SparseArray Band, Eigenvalues with Arnoldi
Shift, Richardson $(4v_{h/2}-v_h)/3$), reproduces the same truncation floor to five digits, so a
cross-check sharing $R$ is blind to truncation and only the sweep detects it. Evaluate
eigenfunctions by argument substitution or the function head, never by applying the $u[r]$-form
output to a number (an enormous unevaluated echo). Close with the paired-detector moral: each
silent failure mode (ordering, boundary conditions, truncation) is caught by its own loud
detector (Shift plus Sort, explicit Dirichlet, $R$-doubling).

#### 12.8 [MSc] How do I compute the quantum defect of an alkali-like screened-Coulomb potential?

Sodium-like screened Coulomb $V(r)=-\tfrac1r\left(1+(Z-1)e^{-r/a}\right)$ with $Z=11$ and $a$
tuned at authoring so the valence $s$ defect lands near sodium's $\delta_0\approx1.35$. Levels per
$l$ channel by the C2 recipe (explicit DirichletCondition at both ends, Arnoldi Shift below the
deepest level then Sort, MaxCellMeasure 0.01, domain literally at $r=0$); high-$n$ levels sit near
the Rydberg accumulation point, so the $R$-doubling sweep is mandatory per level, with $R$ well
beyond $2n^2$ of the highest wanted state ($R=160$ for $n=5$ at floor accuracy). The sorted output
is the channel spectrum counted from its bottom, so the principal quantum number is assigned, not
assumed: the $k$-th state of the sorted $l$ channel ($k$ from 0) carries $k$ radial nodes, hence
$n=k+l+1$. At $Z=11$ the bottom of each channel is a deeply bound core level near $-Z^2/2$, which
is not a perturbed Coulomb state at all, so the defect table covers only the valence series (in
this model $n\ge3$ in all three channels, the boundary located by where the levels rejoin the
Rydberg pattern) and the core levels stand beside it as the contrast. Extract
$\delta_l=n-1/\sqrt{-2E_{nl}}$ from eigenvalues alone: eigenfunction normalization is unprobed in
the C2 verdict, and the eigenvalue-only extraction avoids it entirely (any wavefunction-based
matrix element would have to check its norm explicitly first). Verification that can refute: turn
the screening off, $(Z-1)\to0$ on the same grid, and every extracted $\delta_l$ must collapse to
zero with the levels back on $-1/(2n^2)$; the Rydberg-formula content is the near-constancy of
$\delta_l$ across the valence $n$, shown both as the per-level table and as a one-parameter
FindFit of those same levels to $-1/(2(n-\delta_l)^2)$, whose residual measures the leftover
$n$-dependence the single-parameter formula cannot absorb. Close with the hierarchy
$\delta_0>\delta_1>\delta_2\approx0$ and with why the core levels sit outside the table: only
penetrating orbits see the unscreened core, and a state bound inside it is no longer Rydberg-like.

#### 12.9 [MSc] How do I apply WKB to the radial equation with the Langer correction $l(l+1)\to(l+\tfrac12)^2$ and recover the Coulomb and oscillator spectra?

Radial WKB action $S(E)=\int_{r_-}^{r_+}\sqrt{2E+\tfrac2r-\tfrac{\lambda^2}{r^2}}\,dr$ run twice,
$\lambda=l+\tfrac12$ (Langer) against $\lambda=\sqrt{l(l+1)}$ (uncorrected), quantized by
$S=\pi(n_r+\tfrac12)$. The Kepler substitution $r=a(1-e\cos u)$ is mandatory: raw Integrate over
the exact turning points stalls (70 s timeout in the C9 verdict), and the Newton-Leibniz
antiderivative with Limit across a turning point silently returns 0 with no message; under the
substitution the action closes exactly to $S=\pi/\sqrt{-2E}-\pi\lambda$, and Solve gives
$E=-1/(2(n_r+l+1)^2)$ exactly with the Langer $\lambda$, while the uncorrected form is off by 6.0
percent at $(n_r,l)=(1,1)$: the on/off contrast is the load-bearing exhibit. The oscillator side
follows the same discipline: substitute the turning points out of
$S(E)=\int\sqrt{2E-r^2-\lambda^2/r^2}\,dr$ via $s=r^2$, closing to $\tfrac{\pi}{2}(E-\lambda)$ and
hence $E=2n_r+l+\tfrac32$ again only with the Langer $\lambda$ (this quadrature is outside the
probed set, so its closing is tagged (verify at authoring) under the same trap discipline).
Verification reusing the definitions: guarded 30-digit NIntegrate of the same integrand between
numeric turning points at several $(n_r,l)$ pairs must return $\pi(n_r+\tfrac12)$ (the probed
Coulomb family agrees to about $10^{-28}$), and the uncorrected $\lambda$ visibly fails the same
test. Close by naming the accident: with the Langer correction radial WKB is exact for both
Coulomb and oscillator, the two central potentials whose classical orbits close.

## Part 13 Plan: Electromagnetic coupling

Six questions, all MSc. Class census per `Route-Table.md`: C1: 13.2 (exactly solvable stationary
ODE, reached only after a hand reduction); C0: 13.1, 13.3, 13.4, 13.5, 13.6 (no differential
equation; WL machinery named per entry, C0 kernel facts binding).

Everything in this part flows from minimal coupling for a unit charge in natural units: the gauged
Schrodinger equation $i\partial_t\psi=[\tfrac12(\hat p-A)^2+\phi]\psi$, invariant under
$A\to A+\nabla\chi$, $\phi\to\phi-\partial_t\chi$, $\psi\to e^{i\chi}\psi$; the flux quantum
$\Phi_0=2\pi$, which turns enclosed flux into the phase $e^{2\pi i\Phi/\Phi_0}$; degenerate
perturbation theory, which replaces the naive first-order shift by the secular problem of
diagonalizing the perturbation inside the degenerate manifold; and the Berry phase
$\gamma=\oint i\langle\psi\vert\nabla_R\psi\rangle\cdot dR$, the loop holonomy left after the
dynamical phase is removed. The C0 row's kernel facts bind throughout (in code `E` is Euler's
number, so field strengths and energies need other names; `Simplify` and `Integrate` need explicit
reality and positivity assumptions to close residuals); 13.2 alone touches a differential equation
and cites the C1 route.

#### 13.1 [MSc] How do I implement minimal coupling $\hat p\to\hat p-A$ and show how $\psi$ transforms under a gauge change?

Work the uniform electric field $F$ in two gauges: scalar gauge $\phi=-Fx$, $A=0$, and velocity
gauge $\phi=0$, $A(t)=-Ft$, connected by the gauge function $\chi(x,t)=-Ftx$ with
$\psi'=e^{i\chi}\psi$. Prove the transformation generically first: substitute $e^{i\chi}\psi$ into
the velocity-gauge Schrodinger equation with `D`, eliminate $\partial_t\psi$ through the
scalar-gauge equation as a replacement rule, and `Simplify` the residual to zero for arbitrary
$\psi(x,t)$; since $\chi$ is real, `ComplexExpand` shows $\vert e^{i\chi}\psi\vert^2=\vert\psi\vert^2$
at the same generic level. Then carry a concrete state from strong-field physics, the accelerating
Volkov plane wave $\psi_k=\exp[i(k+Ft)x-\tfrac{i}{2}\int_0^t(k+Fs)^2\,ds]$ in the scalar gauge and
its partner $e^{i\chi}\psi_k$ in the velocity gauge, and let `D` plus `Simplify` under
`Element[{x, t, k}, Reals]` drive both Schrodinger residuals to zero, the check that could refute
either gauge's Hamiltonian. Close on what is not invariant: $\langle\hat p\rangle$ differs between
the gauges while the mechanical momentum $\langle\hat p-A\rangle$ agrees, so only the latter is the
measurable velocity.

#### 13.2 [MSc] How do I derive the Landau levels of a charged particle in a uniform magnetic field?

In the Landau gauge $A=(0,Bx)$ the ansatz $e^{iky}f(x)$ reduces
$H=\tfrac12[p_x^2+(p_y-Bx)^2]$ to a shifted oscillator of frequency $B$ centered at $x_0=k/B$ by
completing the square; per the C1 verdict this reduction is hand algebra, so present it as such in
prose and certify the operator identity with a `Simplify` residual before handing the shifted
oscillator to the C1 route (`DSolve` plus the manual termination read-off, `HermiteH`
eigenfunctions, `FullSimplify` residual zero), giving $E_n=B(n+\tfrac12)$ with no $k$ dependence.
In the symmetric gauge $A=\tfrac{B}{2}(-y,x)$ define the lowest-level family
$\psi_m\propto(x+iy)^m e^{-B(x^2+y^2)/4}$ (the chirality of $m$ against the sign of $B$: verify at
authoring) and let a `Table` of `Simplify` residuals show $(H-\tfrac{B}{2})\psi_m=0$ for
$m=0,\dots,6$: one spectrum from two gauges, plus an explicit $m$ degeneracy, is the refuting
check. State the degeneracy per area $B/2\pi$ with its flux-counting argument (Landau-gauge strip:
$k=2\pi j/L_y$ with guiding centers $k/B\in[0,L_x]$ gives $N=BL_xL_y/2\pi=\Phi/\Phi_0$), noting
the finite $m$ family exhibits but cannot exhaust the infinite degeneracy. Close with the quantum
Hall reading: the filling factor counts electrons per flux quantum.

#### 13.3 [MSc] How do I compute the Aharonov-Bohm phase acquired around a flux line?

Put the unit charge on a ring of radius $R$ threaded by flux $\Phi$:
$H=\tfrac{1}{2R^2}(-i\partial_\phi-\Phi/\Phi_0)^2$ with $\Phi_0=2\pi$. Single-valued eigenfunctions
$e^{im\phi}$, $m\in\mathbb Z$, give $E_m=\tfrac{1}{2R^2}(m-\Phi/\Phi_0)^2$ by a `D` plus `Simplify`
residual, so the Aharonov-Bohm phase is not bookkeeping but a spectrum. The refuting check is
twofold: exact periodicity, with `Simplify` showing $E_m(\Phi+\Phi_0)-E_{m+1}(\Phi)=0$ so the
spectrum as a set is $\Phi_0$-periodic, and gauge equivalence, removing $A$ by
$\chi=(\Phi/\Phi_0)\phi$ at the price of the twisted boundary condition
$\psi(2\pi)=e^{2\pi i\Phi/\Phi_0}\psi(0)$, whose free-particle spectrum must reproduce the same
$E_m$. Then the persistent current $I_m=-\partial E_m/\partial\Phi$ via `D`: for the ground state a
sawtooth in $\Phi$, flipping sign at half flux where two levels cross. Close on the laboratory:
equilibrium currents in mesoscopic normal-metal rings, periodic in exactly one flux quantum.

#### 13.4 [MSc] How do I compute the normal (orbital) Zeeman splitting of hydrogen perturbatively from the $-\tfrac12 B L_z$ coupling?

On the fourfold degenerate $n=2$ manifold $\{2s,2p_{-1},2p_0,2p_1\}$, built from explicit
`LaguerreL` radial parts and `SphericalHarmonicY` angular parts, degenerate perturbation theory
demands the secular matrix of the perturbation inside the manifold. Assemble the $4\times4$ matrix
of $-\tfrac12BL_z$ with $L_z=-i\partial_\phi$ acting on the wavefunctions and raw `Integrate` over
$(r,\theta,\phi)$ with the volume element $r^2\sin\theta$ and assumptions $B>0$; the refuting check
is that this raw-`Integrate` matrix comes out exactly diagonal (compare against `DiagonalMatrix`)
with entries $-\tfrac12Bm$, so `Eigenvalues` returns the shifts $\{0,\tfrac12B,0,-\tfrac12B\}$ and
the triplet pattern of three line components at $0,\pm\tfrac12B$. Say why the secular problem is
trivial here, the perturbation commutes with $L_z$ and the basis is already adapted, so the
degenerate manifold splits without mixing. Close by restoring units to the Bohr magneton scale
$\mu_B B$ for a laboratory field, with the spin-free caveat stated: real hydrogen shows the
anomalous pattern once spin enters.

#### 13.5 [MSc] How do I compute the linear Stark effect on the degenerate $n=2$ hydrogen manifold?

With a field $E$ along $z$ the perturbation is $Ez=Er\cos\theta$ (in code name the field `F0`:
`E` is Euler's number in WL, per the C0 row). Build the $4\times4$ secular matrix on
$\{2s,2p_0,2p_{-1},2p_1\}$ from raw `Integrate` matrix elements of the explicit hydrogen
wavefunctions under $E>0$-style positivity assumptions on the field symbol: parity and the
$\Delta m=0$ selection rule leave the single element $\langle2s\vert z\vert2p_0\rangle=-3$ in
atomic units. The refuting check: `Eigenvalues` of that raw-`Integrate` matrix is exactly
$\{3E,0,0,-3E\}$, a splitting linear in the field. Then define the parabolic pair
$(\vert2s\rangle\mp\vert2p_0\rangle)/\sqrt2$ inside the entry (the separation frame 10.7
establishes) and confirm with a matrix `Dot` that they are the $\pm3E$ eigenvectors while
$2p_{\pm1}$ stay unshifted. Close on the surprise: linearity means the degenerate atom responds as
if it carried a permanent dipole of magnitude $3$ atomic units, where any nondegenerate level
starts at $O(E^2)$.

#### 13.6 [MSc] How do I compute Berry's geometric phase for a state carried adiabatically around a closed loop in parameter space?

Drag the displaced trap $H(x_0,p_0)=\tfrac12(p-p_0)^2+\tfrac12(x-x_0)^2$ around the ellipse
$x_0=a\cos s$, $p_0=b\sin s$, $s\in[0,2\pi]$, whose instantaneous ground state is the coherent
state $\psi=\pi^{-1/4}\exp[-\tfrac12(x-x_0)^2+ip_0x]$ (the phase convention, for example an extra
$e^{-ip_0x_0/2}$, moves the connection but not the holonomy: verify at authoring). Compute the
Berry connection $A(s)=i\langle\psi\vert\partial_s\psi\rangle$ by symbolic `Integrate` over $x$
under reality assumptions, then $\gamma=\int_0^{2\pi}A(s)\,ds$ by a second `Integrate` over the
loop parameter, landing on the enclosed phase-space area $\pi ab$ (overall sign is
convention-dependent: verify at authoring). Cross-check by adiabatic numeric evolution in the
truncated Fock basis with C5 machinery: `NDSolveValue` on the coefficient ODE system at
`PrecisionGoal` and `AccuracyGoal` 10 with an $N$-sweep for truncation, traverse the loop in time
$T$, subtract the dynamical phase $\int_0^T\langle H\rangle\,dt$ from
$\arg\langle\psi(0)\vert\psi(T)\rangle$, and watch the remainder approach $\gamma$ as $T$ grows;
per the C5 verdict norm drift grows with $T$ at default goals, and the finite-$T$ deviation is a
physical finite-rate correction, not solver error, so the honest exhibit is the deviation-versus-$T$
sweep rather than one number. Close on the reading: the holonomy is the symplectic area, the
continuous-variable analogue of the spin solid angle, and the working principle of trapped-ion
geometric phase gates.

## Part 14 Plan: Spin coupled to spatial motion

Six questions, all MSc. Class census per `Route-Table.md`: C0: 14.1, 14.2, 14.3, 14.6 (no
differential equation; WL machinery named per entry, C0 kernel facts binding); C3: 14.4
(time-dependent Schrodinger PDE, the probed spinor gate); C1: 14.5 (exactly solvable stationary
ODE, reached only after a hand reduction).

Everything in this part lives in the tensor product $L^2(\mathbb R)\otimes\mathbb C^2$: a state is
a two-component spinor $\Psi=(\psi_\uparrow,\psi_\downarrow)$ with norm
$\sum_\sigma\int\vert\psi_\sigma\vert^2\,dx=1$, and tracing out position leaves the spin density
matrix $(\rho_s)_{\sigma\sigma'}=\int\psi_\sigma\overline{\psi_{\sigma'}}\,dx$, whose purity
$\operatorname{Tr}\rho_s^2$ quantifies position-spin entanglement. Spin couples to motion through
three recurring structures: the identity $\vec L\cdot\vec S=\tfrac12(J^2-L^2-S^2)$, which converts
spin-orbit coupling into $j$-dependent numbers; the Clebsch-Gordan change of basis
$\vert j,m_j\rangle=\sum_{m_l,m_s}C^{j m_j}_{l m_l; s m_s}\vert l,m_l\rangle\vert s,m_s\rangle$;
and the Pauli Hamiltonian $H=\tfrac12(\hat p-A)^2-\tfrac{g}{2}\vec B\cdot\vec S$ with $g=2$, whose
gradient-field limit is Stern-Gerlach and whose uniform-field spectrum is Landau plus Zeeman. The
Zeeman sign is fixed by the family convention of a unit positive charge (`Shared-Brief.md`,
"Charge and sign conventions"): writing the kinetic term as $(\hat p-A)$ forces
$-\tfrac{g}{2}\vec B\cdot\vec S$. The C0 row's kernel facts bind throughout (`E` is Euler's number,
so energies and field strengths need other names; `Integrate` and `Simplify` need explicit reality
and positivity assumptions to close norms and residuals); 14.4 evolves a PDE and cites the C3
route, 14.5 cites the C1 route.

#### 14.1 [MSc] How do I build a two-component spinor wavefunction $\psi_\sigma(x)$ by tensoring a spin-1/2 with the spatial state?

Build the pinned entangled spinor as a one-parameter family: normalize each branch on its own,
$\hat\psi_\uparrow=\operatorname{sech}(x-a)/\sqrt2$ (an even nodeless profile, the single bound
state of $V=-\operatorname{sech}^2x$) and
$\hat\psi_\downarrow=\sqrt{3/2}\,\operatorname{sech}(x+a)\tanh(x+a)$ (an odd one-node profile, the
excited state of the deeper well $V=-3\operatorname{sech}^2x$), then mix them with a symbolic angle,
$\Psi=\cos\theta\,\hat\psi_\uparrow\vert\uparrow\rangle+\sin\theta\,\hat\psi_\downarrow\vert\downarrow\rangle$,
whose $\theta=\pi/6$ member is the ledger's pinned state. Every spin observable follows from
`Integrate` under `Assumptions -> {a > 0, Element[x, Reals]}`: the populations are
$p_\uparrow=\cos^2\theta$, $p_\downarrow=\sin^2\theta$ with no $a$ dependence at all (both
$\int\operatorname{sech}^2u\,du=2$ and $\int\operatorname{sech}^2u\tanh^2u\,du=2/3$ are
translation invariant), so the whole geometry sits in the single real overlap
$s(a)=\int\hat\psi_\uparrow\hat\psi_\downarrow\,dx$ (closed form in $a$: verify that `Integrate`
closes it) and the purity is
$\operatorname{Tr}\rho_s^2=1-\tfrac12\sin^2(2\theta)\,[1-s(a)^2]$. The refuting contrast is the
product state $\phi(x)(\alpha\vert\uparrow\rangle+\beta\vert\downarrow\rangle)$, where the identical
pipeline must return purity exactly 1, with $\operatorname{Tr}\rho_s=1$ standing on both and an
`NIntegrate` spot check of $s$ at $a=1$ reusing the definitions. Two limits carry the physics:
$s$ vanishes both at $a=0$ (odd integrand, exact parity zero) and as $a\to\infty$ (disjoint
supports), so the entanglement is maximal at zero displacement and at infinite displacement alike
and is weakest in between, and the floor is $\tfrac12$ only at $\theta=\pi/4$, while the pinned
$\theta=\pi/6$ member floors at $\tfrac58$. Close on the measurement reading: an apparatus probing
spin alone sees the mixed $\rho_s$, and the missing purity is exactly the which-branch information
stored in position, which parity can supply as completely as distance does.

#### 14.2 [MSc] How do I add the spin-orbit term $\propto \vec L\cdot\vec S$ and compute the hydrogen fine structure?

Assemble the pinned hydrogen $2p$ fine structure from its radial and angular factors. Radially,
define $R_{nl}(r)$ from `LaguerreL` (the associated-Laguerre argument convention and the resulting
normalization constant: verify at authoring) and earn the norms by `Integrate` with the $r^2$
measure; then compute $\langle r^{-3}\rangle_{nl}$ from those same definitions over the family
$(n,l)\in\{(2,1),(3,1),(3,2),(4,1)\}$ and check that all four reproduce
$1/[n^3\,l(l+\tfrac12)(l+1)]$, a definition-reusing test that exposes a wrong radial function or a
wrong normalization at once and hands over the $n^{-3}$ scaling of the whole fine structure;
$(2,1)$ gives $\langle r^{-3}\rangle_{2p}=\tfrac1{24}$. Angularly, build the $6\times6$
$\vec L\cdot\vec S=\tfrac12(J^2-L^2-S^2)$ from `KroneckerProduct` of the $l=1$ matrices and
`PauliMatrix[i]/2`, and let `Eigenvalues` return $\tfrac12$ four times and $-1$ twice, matching
$\tfrac12[j(j+1)-l(l+1)-s(s+1)]$ at $j=\tfrac32,\tfrac12$ with multiplicities $2j+1$; the two
trace identities $\operatorname{Tr}(\vec L\cdot\vec S)=0$ and
$\operatorname{Tr}[(\vec L\cdot\vec S)^2]=3$ pin the assembly independently of the eigensolver.
With $\xi(r)=\alpha^2/(2r^3)$ the splitting is
$\tfrac32\cdot\tfrac{\alpha^2}{2}\langle r^{-3}\rangle_{2p}=\alpha^2/32$ hartree. Note the edge the
family exposes: at $l=0$ the radial factor diverges while $\langle\vec L\cdot\vec S\rangle$
vanishes, so the product is indeterminate and the $2s$ shift is carried by the Darwin term, not by
spin-orbit. Close in the lab: $\alpha^2/32$ hartree is about $10.9\,\mathrm{GHz}$, the measured
fine-structure doubling of the Lyman-$\alpha$ line.

#### 14.3 [MSc] How do I form the total angular momentum $\vec J=\vec L+\vec S$ and change to the coupled basis?

Construct $J_i=L_i\otimes\mathbb 1_2+\mathbb 1_3\otimes S_i$ with `KroneckerProduct`, the $l=1$
matrices built from the ladder elements $\sqrt{l(l+1)-m(m\pm1)}$ via `SparseArray` and
$S_i=$ `PauliMatrix[i]/2`, and certify the algebra first: `Simplify` of $[J_x,J_y]-iJ_z$ to the
zero matrix. Then form the $6\times6$ change-of-basis matrix $U$ whose rows are
`ClebschGordan[{1, ml}, {1/2, ms}, {j, mj}]` coefficients (argument order and Condon-Shortley
phase: verify at authoring) and let exact linear algebra deliver the payoff: $U.U^\dagger=\mathbb 1$,
and $U.J^2.U^\dagger$, $U.J_z.U^\dagger$ exactly diagonal, symbolically, with entries $\tfrac{15}4$
four times and $\tfrac34$ twice beside the matching $m_j$, every off-diagonal an exact 0. Two
independent refutations: build the same basis by ladder algebra, starting from the stretched state
$\vert\tfrac32,\tfrac32\rangle=\vert1,1\rangle\vert\uparrow\rangle$, stepping down with $J_-$ and
renormalizing, then taking the orthogonal complement inside each $m_j$ block for the $j=\tfrac12$
pair, and compare row by row (a phase-convention slip shows up as a sign mismatch); and check the
general $s=\tfrac12$ closed form
$\vert l\pm\tfrac12,m_j\rangle=\pm\sqrt{\tfrac{l\pm m_j+1/2}{2l+1}}\vert m_j-\tfrac12\rangle\vert\uparrow\rangle+\sqrt{\tfrac{l\mp m_j+1/2}{2l+1}}\vert m_j+\tfrac12\rangle\vert\downarrow\rangle$
(overall sign convention: verify at authoring) against `ClebschGordan` over a `Table` of
$l\in\{1,2,3\}$ and every $m_j$, since `ClebschGordan` at symbolic $l,m_j$ may simply not evaluate.
Close with the forward pointer: these coupled columns are the weak-field eigenvectors 14.6
recovers, and their multiplet structure is what 14.2's eigenvalue count already displayed.

#### 14.4 [MSc] How do I model Stern-Gerlach as a spin-dependent spatial deflection of a spinor wave packet in a field gradient?

Evolve the pinned two-component sech packet $\Psi(x,0)=\tfrac{1}{\sqrt2}\psi_0(x)(1,1)$, with
$\psi_0=\operatorname{sech}(x)/\sqrt2$, under
$i\partial_t\Psi=[-\tfrac12\partial_x^2+\lambda x\,\sigma_z]\Psi$ at the probed gradient
$\lambda=0.2$ to $T=8$, using the C3 spec in full: `NDSolveValue` on the coupled pair with
`Method -> {"MethodOfLines", "SpatialDiscretization" -> {"TensorProductGrid", "MinPoints" -> n, "MaxPoints" -> n, "DifferenceOrder" -> "Pseudospectral"}}`
and `AccuracyGoal -> 10, PrecisionGoal -> 10`, the point pins giving a fixed $dx$ for the
$dx$-weighted grid sums that must replace `NIntegrate` on oscillatory interpolants, and tight goals
being load-bearing (defaults are silently wrong at the $10^{-3}$ level and refinement at default
goals is non-monotone). Size the box in both position and momentum: on $[-40,40]$ with $n=513$ the
walls sit far beyond the deflected centers $\mp\lambda t^2/2=\mp6.4$ plus the dispersive spread,
while $k_{\max}=\pi/dx\approx20$ clears the momentum the tilt pumps in, $\vert\langle p\rangle_\sigma\vert=\lambda T=1.6$,
plus several $k$-widths of $\tilde\psi_0\propto\operatorname{sech}(\pi k/2)$; the linear potential
is not periodic, so the periodic identification is honest only while $\psi$ is negligible at both
walls, and the C3 wrap-around contamination is structural and invisible to norm conservation, so
carry an $n$-doubling and $T$-lengthening sweep (the probe evidence reaches only $t\le2$ at this
$\lambda$). The deliverable is the measurement record, the coherence
$c(t)=\int\psi_\uparrow\overline{\psi_\downarrow}\,dx$ and the spin purity falling from 1 while
both populations hold at $\tfrac12$, so the check must bite on shape, not on centroids: the tilt is
a momentum translation, giving the exact momentum-space form
$c(t)=\int\tilde\psi_0(k-\lambda t)\overline{\tilde\psi_0(k+\lambda t)}\,e^{ik\lambda t^2}\,dk$
(closed form for the sech pair: verify that `Integrate` closes it), and the C3 cross-check route,
a hand-built split-step Fourier evolution of the same coupled pair, must reproduce $\vert c(t)\vert$
itself. The probed $\Omega=0$ Ehrenfest anchor
($\langle x\rangle_\uparrow(2)=-0.39999$ against the exact $-0.4$, total norm drift
$5\times10^{-8}$) stays as a cheap standing sanity, but it is exact for any initial packet in a
linear potential and therefore cannot see a wrong width, wrong dispersion, or an aliased tail.
Close on measurement as dynamics: the coherence dies mainly through momentum separation, so once
$\vert c\vert$ is negligible the beams are which-path records and no downstream spin rotation
revives the interference, the projection postulate emerging from unitary evolution plus a
traced-out position.

#### 14.5 [MSc] How do I solve the Pauli equation for a nonrelativistic spin in an electromagnetic field?

Take uniform $\vec B=B\hat z$ in the Landau gauge $A=(0,Bx,0)$, where the Pauli Hamiltonian
$H=\tfrac12(\hat p-A)^2-\tfrac{g}{2}BS_z$ (the family sign convention for a unit positive charge)
is diagonal in $S_z$: split into the two spin blocks and translate along $y$ with the ansatz
$e^{iky}f(x)$. Per the C1 verdict this member is covered only after the hand reduction: complete the
square to a shifted oscillator of frequency $B$ centered at $x_0=k/B$ as hand algebra, certify the
reduction with a `Simplify` residual on the operator identity, and only then hand the oscillator to
the C1 route, `DSolve` with symbolic energy followed by the manual series-termination read-off and
`HermiteH` eigenfunctions with a `FullSimplify` residual of zero; the read-off is not a stylistic
choice, `DSolve` handed the full decay conditions hangs silently rather than failing (the C1 gate
timed out at 60 s with no message), and the whole-line `DEigensystem` domain form is a free exact
confirmation only if it accepts a symbolic $B$ (the gate probed it at numeric coefficients: verify
at authoring). The spectrum $E_{n,m_s}=B(n+\tfrac12)-\tfrac{g}{2}Bm_s$ becomes at $g=2$ the two
interleaved ladders $E_{n,\uparrow}=Bn$ and $E_{n,\downarrow}=B(n+1)$, so `Simplify` of
$E_{n,\downarrow}-E_{n+1,\uparrow}$ is 0 symbolically in $B$ and $n$, and the discriminating
refutation is that `Solve` on $E_{n,\downarrow}=E_{n+1,\uparrow}$ for $g$ returns exactly $g=2$:
the coincidence is a property of the gyromagnetic ratio, not of magnetism. Add one numeric
confirmation at $B=1$, an `NDEigenvalues` run per spin block on a truncated box with explicit
`DirichletCondition` (omitting boundary conditions is a silent Neumann-zero answer, wrong at the
percent level on a ground state, per the probed gates), whose sorted merge is
$\{0,1,1,2,2,\dots\}$: the coincidence appears as numerical degeneracy. The infinite per-level
Landau degeneracy was not probed and its flux counting belongs to 13.2; here only the $k$
independence of $E_{n,m_s}$, read directly off the closed form, records it. Close at the bottom of
the ladder: at $g=2$ the state $(n=0,m_s=+\tfrac12)$ sits at exactly zero energy, the zero-point
cost cancelled by the Zeeman gain, the nonrelativistic shadow of the Dirac value of $g$ that 22.3
derives.

#### 14.6 [MSc] How do I compute the anomalous Zeeman effect and the Paschen-Back crossover, where spin-orbit and the magnetic field compete?

On the hydrogen $2p$ manifold build the exact $6\times6$
$H=\xi\,\vec L\cdot\vec S+\tfrac{B}{2}(L_z+2S_z)$ from `KroneckerProduct` with $\xi$ and $B$
symbolic (rebuild the $l=1$ and Pauli matrices in-entry; the answer stays self-contained). Certify
the block structure first, `Simplify` of $[H,J_z]$ to the zero matrix, then make the blocks explicit
rather than hoping the eigensolver finds them: order the uncoupled basis by $m_j=m_l+m_s$ with
`Ordering` (or conjugate by the coupled-basis matrix of 14.3) so $H$ becomes block diagonal with
sizes 1, 2, 2, 1, and take `Eigenvalues` on the $1\times1$ and $2\times2$ blocks, where the
Breit-Rabi square roots are guaranteed closed forms in $\xi$ and $B$; a direct `Eigenvalues` on the
full two-parameter $6\times6$ may return `Root` objects instead (verify at authoring). The two
limits are the refuting checks because they rest on independent physics: `Series` in $B$ at fixed
$\xi$ must reproduce the anomalous pattern $E_{j,m_j}\approx E_j+g_j m_j B/2$ with Lande factors
$g_{3/2}=\tfrac43$ and $g_{1/2}=\tfrac23$, matching the coupled-basis diagonal elements of
$(L_z+2S_z)/2$, while the opposite `Series` in $\xi$ at fixed $B$ must land on the Paschen-Back
pattern $B(m_l+2m_s)/2$ with its characteristic degeneracies; and the coupling-resolved identity
$\operatorname{Tr}(H^2)=3\xi^2+\tfrac52B^2$, exact and symbolic in both, catches an assembly error
that a traceless sum rule would let through. `Plot` the six branches against $B$ at fixed $\xi$ to
exhibit the crossover near $B\sim\xi$: branches of different $m_j$ cross freely while the paired
same-$m_j$ branches avoid crossing, two-level repulsion inside each block. Close in the
spectroscope: weak-field lines splitting by $g_j$ rather than by 1 is the historical anomalous
Zeeman effect, and the field where the pattern reorganizes, Zeeman energy comparable to the
$\alpha^2/32$ splitting of 14.2, is about $0.8\,\mathrm{T}$ for hydrogen.

## Part 15 Plan: Identical particles in continuous space

Five questions, all MSc. Class census (`Route-Table.md`): C0 for 15.1, 15.2, 15.3, 15.4 (no
differential equation, symbolic `Integrate` and matrix algebra); C6 for 15.5 (mean-field SCF).

### Common ground

Two identical particles admit only states of definite exchange symmetry,
$\psi(x_2,x_1)=\pm\psi(x_1,x_2)$: from distinct orbitals $\phi_a,\phi_b$ the normalized pair states
are $\psi_{\pm}=(\phi_a(x_1)\phi_b(x_2)\pm\phi_b(x_1)\phi_a(x_2))/\sqrt{2}$, and for $N$ fermions the
antisymmetrizer is the Slater determinant $\Psi=\det[\phi_i(x_j)]/\sqrt{N!}$. Statistics becomes
observable through the two-body density $\rho_2(x_1,x_2)=|\psi(x_1,x_2)|^2$ at coincidence, and
through the energy via the direct and exchange integrals
$J=\iint|\phi_a(\mathbf r_1)|^2\,v_{12}\,|\phi_b(\mathbf r_2)|^2$ and
$K=\iint\phi_a^{*}(\mathbf r_1)\phi_b(\mathbf r_1)\,v_{12}\,\phi_b^{*}(\mathbf r_2)\phi_a(\mathbf r_2)$,
which split spatially symmetric from antisymmetric pairs as $J\pm K$; when the pair interaction is too
hard to diagonalize exactly, a mean field closes it self-consistently, $F[\phi]\,\phi=\varepsilon\phi$.

#### 15.1 [MSc] How do I build a two-particle wavefunction $\psi(x_1,x_2)$ and symmetrize or antisymmetrize it?

Work in the infinite well of width $L$ with orbitals $\chi_n(x)=\sqrt{2/L}\,\sin(n\pi x/L)$, $n=1,2$,
distinct so exchange has something to act on. Define the product $\chi_1(x_1)\chi_2(x_2)$ and the pair
$\psi_{\pm}=(\chi_1(x_1)\chi_2(x_2)\pm\chi_2(x_1)\chi_1(x_2))/\sqrt{2}$; earn the normalization with
`Integrate` of $|\psi_{\pm}|^2$ over $[0,L]^2$ under $L>0$ (exactly 1, riding on orbital
orthonormality), and earn the symmetry by applying the swap rule `{x1 -> x2, x2 -> x1}` and
confirming with `Simplify` that $\psi_{+}$ is invariant and $\psi_{-}$ flips sign. The fermion node
line is `Simplify` of $\psi_{-}$ at $x_2\to x_1$ returning 0 identically; the boson coincidence
density $|\psi_{+}(x,x)|^2=2\,\chi_1(x)^2\chi_2(x)^2$, exactly twice the distinguishable product's
coincidence value, quantifies the enhancement. All three checks reuse the defined $\psi_{\pm}$ and can
refute: a wrong prefactor fails the norm, a dropped sign fails the swap or the node. Close on the
equal-orbital edge: with both orbitals $\chi_1$ the antisymmetric combination vanishes identically
(exclusion falls out of the algebra), while the naive symmetric $1/\sqrt{2}$ state has norm
$\sqrt{2}$, so the general prefactor is $1/\sqrt{2(1+\delta_{ab})}$.

#### 15.2 [MSc] How do I exhibit the exchange hole for fermions and the bunching for bosons (Pauli exclusion versus enhancement)?

Use the oscillator pair $\phi_0(x)=\pi^{-1/4}e^{-x^2/2}$ and
$\phi_1(x)=\sqrt{2}\,\pi^{-1/4}\,x\,e^{-x^2/2}$ (built from `HermiteH`), the same orbitals for both
statistics so only the sign differs. Form $\rho_2^{\pm}=|\psi_{\pm}(x_1,x_2)|^2$ and reduce with
`Simplify` under `Element[{x1, x2}, Reals]`: the fermion density collapses to the closed form
$\rho_2^{-}=(x_1-x_2)^2\,e^{-x_1^2-x_2^2}/\pi$, whose quadratic zero along $x_1=x_2$ is the exchange
hole (exhibit the $s^2$ opening with `Series` in $s=x_1-x_2$), while the boson coincidence value
$\rho_2^{+}(x,x)=2\,\phi_0(x)^2\phi_1(x)^2$ doubles the distinguishable reference
$\rho_2^{\mathrm{d}}=\tfrac12\left(|\phi_0(x_1)\phi_1(x_2)|^2+|\phi_1(x_1)\phi_0(x_2)|^2\right)$:
factor 0 against factor 2, both closed forms. For the pair correlation compute the marginal
$\rho(x)=\int\rho_2^{\pm}\,dx_2$ by `Integrate` and form
$g(x_1,x_2)=\rho_2^{\pm}/(\rho(x_1)\rho(x_2))$; the refuting check is that the marginal must equal
$\tfrac12(\phi_0(x)^2+\phi_1(x)^2)$ for both signs (orthogonality kills the cross term), so any
normalization slip in $\psi_{\pm}$ surfaces immediately. Close in the lab: coincidence counters on
cold-atom pairs read exactly these two numbers, a doubled accidental rate for bosons and an
antibunching null for fermions, the atomic Hanbury Brown-Twiss experiment.

#### 15.3 [MSc] How do I build a Slater determinant for $N$ fermions in a well or oscillator?

Fill the well orbitals $\chi_n(x)=\sqrt{2/L}\,\sin(n\pi x/L)$, $n=1,2,3$, and build
$\Psi(x_1,x_2,x_3)=\det[\chi_n(x_j)]/\sqrt{3!}$ with `Det` on the `Table`-built $3\times3$ orbital
matrix, a genuine $N=3$ object rather than the two-body shortcut. Verify antisymmetry under one
explicit transposition, `Simplify` of $(\Psi$ at `{x1 -> x2, x2 -> x1}`$)+\Psi$ to 0, and the
coincidence node, `Simplify` of $\Psi$ at $x_2\to x_1$ to 0. The load-bearing check is the
determinant identity for orthonormal orbitals: the one-body density
$\rho(x)=3\iint|\Psi(x,x_2,x_3)|^2\,dx_2\,dx_3$ must equal $\sum_{n=1}^{3}\chi_n(x)^2$; compute the
left side by nested `Integrate` under $0<x<L$ and `Simplify` the difference against
`Sum[chi[n, x]^2, {n, 3}]` to 0, which refutes a wrong $1/\sqrt{3!}$ or a non-orthonormal orbital
choice while stating the physics: filled levels add in density with no interference. Close by
occupying $\chi_1$ twice in the matrix: `Det` with two equal rows returns 0 identically, Pauli
exclusion as a determinant fact rather than a postulate.

#### 15.4 [MSc] How do I compute the exchange energy and the para/ortho splitting of a two-electron (helium-like) model?

Use hydrogenic orbitals with symbolic nuclear charge $Z$: $\phi_{1s}=Z^{3/2}e^{-Zr}/\sqrt{\pi}$ and
$\phi_{2s}=Z^{3/2}(2-Zr)\,e^{-Zr/2}/(4\sqrt{2\pi})$, orthonormality earned first by `Integrate` with
$Z>0$. Route both Coulomb integrals through the multipole expansion
$1/r_{12}=\sum_{l}r_{<}^{l}\,r_{>}^{-l-1}\,P_l(\cos\theta_{12})$: the angular integrals of s orbitals
kill every $l>0$ term, so $J$ and $K$ collapse to nested radial quadratures with kernel $1/r_{>}$,
done as two explicit `Integrate` pieces (inner $\int_0^{r_1}dr_2$ weighted by $1/r_1$ plus
$\int_{r_1}^{\infty}dr_2$ weighted by $1/r_2$, then the outer $r_1$ integral, assumptions $Z>0$,
$r_1>0$); expected closed forms $J=\tfrac{17}{81}Z$ and $K=\tfrac{16}{729}Z$ (verify at authoring),
splitting $2K=\tfrac{32}{729}Z$. Cross-check with `NIntegrate` on the same nested radial integrand at
$Z=2$, reusing the orbital definitions; the physics check is the sign chain $J>K>0$ ($K$ is the
Coulomb self-energy of the overlap density $\phi_{1s}\phi_{2s}$, hence positive), which places the
ortho triplet at $J-K$ below the para singlet at $J+K$ from a spin-independent Hamiltonian. Close
with the number: at $Z=2$, $2K\approx0.088$ Hartree $\approx2.4$ eV, overshooting the measured helium
$1s2s$ splitting near $0.8$ eV because the unscreened $2s$ orbital sits too deep; screening is
exactly what a self-consistent field supplies.

#### 15.5 [MSc] How do I solve the Hartree-Fock self-consistent-field equations for a model atom numerically?

This is the C6 SCF member (`Route-Table.md`, C6 verdict, probe p6): softened-Coulomb one-dimensional
helium on a grid, one-body Hamiltonian by finite differences via `SparseArray` with `Band`, the
electron-electron kernel a dense matrix on the $dx$-weighted grid; the softening parameters and the
closed-shell reduction are re-derived in the part, showing that for the doubly occupied orbital the
exchange term equals the direct self-term, so RHF collapses to the Hartree operator
$F[\phi]=h+J[\phi]$. Iterate with `NestList`, dense `Eigensystem` per step, linear mixing: the probed
pattern converges geometrically with contraction $\approx0.135$ at plain mixing, but convergence is
declared only on the orbital change norm $\|\Delta\phi\|$ (`Norm`), never from an eigenvalue plateau
(the C6 refutation: a plateau can mask a limit cycle); the converged orbital must pass the stationary
residual $\|(h+J[\phi])\phi-\varepsilon\phi\|$, two mixings $\alpha=1$ and $\alpha=0.5$ must land on
one fixed point, and the probe anchors $\varepsilon\approx-0.75032$,
$E_{\mathrm{HF}}=2\varepsilon-U\approx-2.22450$ with a parity-symmetric orbital (energies named `en`,
never `E`). Grade against the exact benchmark: the two-electron grid Hamiltonian
$h\otimes\mathbb{1}+\mathbb{1}\otimes h+V_{\mathrm{ee}}$ assembled with `KroneckerProduct` on
`SparseArray`, ground state by the C2-certified recipe (`Eigenvalues` with Arnoldi and `"Shift"`
below the ground energy, since unshifted magnitude ordering silently misses the negative spectrum)
plus a grid-doubling sweep. The honest close is the correlation energy
$E_{\mathrm{corr}}=E_{\mathrm{exact}}-E_{\mathrm{HF}}<0$, the part of the physics no mean field can
reach, pinned to a number by the exact grid.

## Part 16 Plan: Density operators, mixed states, and the Wigner function

Six questions, all MSc. Class census (`Route-Table.md`): C0 for all of 16.1 through 16.6 (kernels,
quadrature, and integral transforms; no differential-equation solver).

### Common ground

A mixed state enters as a kernel $\rho(x,x')=\sum_{k}p_k\,\psi_k(x)\,\overline{\psi_k(x')}$ with
unit trace $\int\rho(x,x)\,dx=1$ and purity
$\operatorname{Tr}\rho^2=\iint|\rho(x,x')|^2\,dx\,dx'$, equal to 1 exactly when the state is pure;
a subsystem of a pure pair state is the partial trace
$\rho_A(x,x')=\int\Psi(x,y)\,\overline{\Psi(x',y)}\,dy$. Phase space sees $\rho$ through the
Wigner-Weyl transform $W(x,p)=\frac{1}{2\pi}\int e^{-ipy}\,\rho(x+\tfrac{y}{2},x-\tfrac{y}{2})\,dy$,
whose marginals are the position and momentum densities and which obeys $|W|\le 1/\pi$; smoothing
orders the representations, with $Q(\alpha)=\frac{1}{\pi}\langle\alpha|\rho|\alpha\rangle\ge0$ at
the regular end and $P$, defined by $\rho=\int P(\alpha)\,|\alpha\rangle\langle\alpha|\,d^{2}\alpha$,
at the singular end. Evolution is the Moyal equation
$\partial_{t}W=-p\,\partial_{x}W+V'(x)\,\partial_{p}W-\tfrac{1}{24}V'''(x)\,\partial_{p}^{3}W+\ldots$,
which truncates to the classical Liouville flow whenever $V$ is quadratic. Oscillator entries set
$\omega=1$, so the coherent state $|\alpha\rangle$ is a vacuum-width Gaussian centered at
$(x_0,p_0)=(\sqrt{2}\,\operatorname{Re}\alpha,\sqrt{2}\,\operatorname{Im}\alpha)$.

#### 16.1 [MSc] How do I represent a density operator as a kernel $\rho(x,x')$ and compute its purity $\iint|\rho(x,x')|^2\,dx\,dx'$?

Build the two levels of the Poschl-Teller well $V=-3\operatorname{sech}^2x$ in closed form,
$\psi_0(x)=\tfrac{\sqrt{3}}{2}\operatorname{sech}^2x$ at $E_0=-2$ and
$\psi_1(x)=\sqrt{3/2}\,\operatorname{sech}x\tanh x$ at $E_1=-\tfrac{1}{2}$, earning both by a
`FullSimplify` residual of the stationary equation and unit norms by `Integrate` (energies named
`en`, never `E`). Form the discriminating pair: the equal-population mixture
$\rho_{\mathrm{mix}}(x,x')=\tfrac{1}{2}\psi_0(x)\psi_0(x')+\tfrac{1}{2}\psi_1(x)\psi_1(x')$ against
the coherent superposition $\psi=(\psi_0+i\psi_1)/\sqrt{2}$ with kernel
$\rho_{\mathrm{sup}}(x,x')=\psi(x)\overline{\psi(x')}$; the relative phase $i$ makes the two
diagonals coincide pointwise (`FullSimplify` of the difference to 0), so no position histogram
separates the pair. Purity does: `Integrate` the double integral $\iint|\rho|^2\,dx\,dx'$ under
reality assumptions to exactly $\tfrac{1}{2}$ for the mixture and $1$ for the superposition, and
cross-check both numbers as $\sum_{k}p_k^2$ from the `Integrate`-verified orthonormality of the
levels, reusing the same definitions, so a wrong kernel fails one route or the other. Close on
where the difference lives: `FullSimplify` $\rho_{\mathrm{sup}}-\rho_{\mathrm{mix}}$ down to the
pure off-diagonal cross term $\propto i\,(\psi_1(x)\psi_0(x')-\psi_0(x)\psi_1(x'))$, the coherence
a position measurement never sees and decoherence destroys.

#### 16.2 [MSc] How do I obtain a reduced density matrix by tracing out one particle (integrating over its coordinate)?

Take the coupled pair $H=\tfrac{1}{2}(p_1^2+p_2^2)+\tfrac{1}{2}(x_1^2+x_2^2)+\tfrac{1}{2}g(x_1-x_2)^2$
with normal modes $x_{\pm}=(x_1\pm x_2)/\sqrt{2}$ at $\omega_{+}=1$, $\omega_{-}=\sqrt{1+2g}$, and
write the exactly Gaussian ground state
$\Psi(x_1,x_2)=(\omega_{+}\omega_{-}/\pi^2)^{1/4}\,e^{-\omega_{+}x_{+}^2/2-\omega_{-}x_{-}^2/2}$,
earned by the 2D stationary residual (`D` plus `FullSimplify`, energy
$\tfrac{1}{2}(\omega_{+}+\omega_{-})$). Trace out the partner by `Integrate` over $y$ with
`Assumptions` $\omega_{\pm}>0$: $\rho_A(x,x')=\int\Psi(x,y)\,\overline{\Psi(x',y)}\,dy$ is again
Gaussian, $\propto e^{-\gamma(x^2+x'^2)/2+\beta xx'}$ with $\gamma,\beta$ read off by `Coefficient`.
Purity comes from the Gaussian double `Integrate` of $|\rho_A|^2$, a
$2\sqrt{\omega_{+}\omega_{-}}/(\omega_{+}+\omega_{-})$-type closed form (verify at authoring); the
independent route that could refute it computes $\langle x^2\rangle$ and $\langle p^2\rangle$ from
the same kernel (`D` across the diagonal before setting $x'=x$, then `Integrate`), forms the
symplectic eigenvalue $\nu=\sqrt{\langle x^2\rangle\langle p^2\rangle}$, and demands purity
$=1/(2\nu)$, after which the entanglement entropy follows as
$S=(\nu+\tfrac{1}{2})\ln(\nu+\tfrac{1}{2})-(\nu-\tfrac{1}{2})\ln(\nu-\tfrac{1}{2})$ (verify at
authoring). Anchor with `Limit` at $g\to0$: purity 1 and $S=0$ exactly; a large-$g$ `Series` shows
the purity falling like $(\omega_{+}/\omega_{-})^{1/2}$ toward zero. Close: the global state stays
pure at every $g$, yet each oscillator alone is strictly mixed once $g\neq0$; mixedness of a
subsystem is the operational face of entanglement.

#### 16.3 [MSc] How do I compute the Wigner quasi-probability function by the Wigner-Weyl transform?

Define the $n=1$ Fock state $\psi_1(x)=\sqrt{2}\,\pi^{-1/4}\,x\,e^{-x^2/2}$ and the transform
$W(x,p)=\frac{1}{2\pi}\int e^{-ipy}\,\psi_1(x+\tfrac{y}{2})\,\overline{\psi_1(x-\tfrac{y}{2})}\,dy$,
then evaluate by one direct `Integrate` under `Element[{x, p}, Reals]`: the closed form
$W_1(x,p)=\frac{1}{\pi}\,(2x^2+2p^2-1)\,e^{-x^2-p^2}$ lands with $W_1(0,0)=-1/\pi$ exactly, the
extremal value permitted by the bound $|W|\le1/\pi$. Three checks reuse the definitions and could
each refute the result: `Integrate` over $p$ must return $|\psi_1(x)|^2$, `Integrate` over $x$ must
return $|\tilde{\psi}_1(p)|^2$ with $\tilde{\psi}_1$ obtained independently by `FourierTransform`,
and the full phase-plane `Integrate` must give 1. `Solve` locates the negative region as the disk
$x^2+p^2<\tfrac{1}{2}$ around the origin. Close: a genuine probability density is nowhere negative,
so the clean $-1/\pi$ certifies the Fock state as non-classical in the strongest phase-space sense,
and no state can dip deeper.

#### 16.4 [MSc] How do I compute the Wigner function of coherent, number, and cat states and exhibit its negativity?

Normalize the even cat at $\alpha=2$: with $x_0=\sqrt{2}\,\alpha$,
$\psi_{\mathrm{cat}}(x)=N\,(e^{-(x-x_0)^2/2}+e^{-(x+x_0)^2/2})$ with $N$ fixed by `Integrate`
($N^{-2}=2\sqrt{\pi}\,(1+e^{-x_0^2})$), and push it through the Wigner kernel of the common ground
beside the coherent state $|\alpha\rangle$ (a displaced vacuum Gaussian, positive everywhere) and
the $n=1$ number state (the negative disk): one `Integrate` under reality assumptions produces the
two blobs at $(\pm x_0,0)$ plus the interference term $\propto e^{-x^2-p^2}\cos(2x_0p)$ at the
origin, fringes of period $\pi/x_0$ in $p$ set by the blob separation $2x_0$. The discriminating
contrast is the incoherent mixture
$\tfrac{1}{2}(|\alpha\rangle\langle\alpha|+|{-\alpha}\rangle\langle{-\alpha}|)$: its kernel through
the same transform must give exactly the two blobs with no fringe term, so the fringe visibility
(interference amplitude over blob amplitude at the origin, unity for the pure cat up to
$e^{-2\alpha^2}$ normalization corrections) is the check that separates coherence from classical
ignorance, with the marginal and normalization checks of the transform repeated on the cat. Close
in the lab: the momentum marginal $\propto e^{-p^2}\cos^{2}(x_0p)$ shows the fringes to a momentum
measurement while the position histogram sees two dead bumps, and the fringe period shrinking as
$1/x_0$ is the seed of why large cats are the first to decohere.

#### 16.5 [MSc] How do I compute the Husimi-$Q$ and Glauber-Sudarshan-$P$ representations?

Run the same trio, the coherent state $\alpha_0=2$, the $n=1$ number state, and the even cat
$\alpha=2$ (each redefined here), through both representations. For Husimi, take
$Q(\alpha)=\tfrac{1}{\pi}|\langle\alpha|\psi\rangle|^2$ with the coherent-state wavefunction
$\phi_{\alpha}(x)=\pi^{-1/4}\,e^{-(x-\sqrt{2}\operatorname{Re}\alpha)^2/2+i\sqrt{2}\operatorname{Im}\alpha\,x}$
(phase convention fixed in-entry) and overlaps by `Integrate` under real
$\operatorname{Re}\alpha,\operatorname{Im}\alpha$: the coherent state gives
$\tfrac{1}{\pi}e^{-|\alpha-\alpha_0|^2}$, the number state $\tfrac{|\alpha|^2}{\pi}e^{-|\alpha|^2}$,
the cat two bumps with exponentially suppressed interference, all manifestly nonnegative. The check
that could refute: $Q$ is the Wigner function smoothed by the vacuum Gaussian, so convolving each
state's independently computed $W$ with $\tfrac{1}{\pi}e^{-x^2-p^2}$ by `Integrate` must reproduce
each $Q$ exactly (convolution normalization to be pinned, verify at authoring). For $P$: the
coherent state is $P=\delta^{(2)}(\alpha-\alpha_0)$, exercised only through `Integrate` against
`DiracDelta` (never `NIntegrate`, which silently returns `0.` on a point measure); the thermal
state at mean occupation $\bar{n}$ joins as the regular anchor,
$P=\tfrac{1}{\pi\bar{n}}\,e^{-|\alpha|^2/\bar{n}}$, a true density; the number and cat states admit
no density, so exhibit them honestly through a regularized representation, the $n=1$ case as the
delta-derivative form $e^{|\alpha|^2}\partial_{\alpha}\partial_{\bar{\alpha}}\delta^{(2)}(\alpha)$
(verify at authoring), stating that the limit exists only as a distribution. Verify all of them by
optical equivalence: `Integrate` each $P$ against $|\alpha|^2$ must return the independently
computed $\langle a^{\dagger}a\rangle$ ($|\alpha_0|^2$, $\bar{n}$, $1$, and the cat's closed form).
Close on the ladder: one Gaussian smoothing separates $P$ from $W$ and another separates $W$ from
$Q$, so singular $P$, negative $W$, positive $Q$ are three verdicts on the same state, and
"classical" means the state whose $P$ already sits at the top rung as a density. Concern: the
ledger's $P$ triple (delta, Gaussian, singular) does not map onto the pinned trio, since the number
state's $P$ is a delta derivative rather than a Gaussian; the Gaussian $P$ belongs to the thermal
state, kept here as the regular anchor beside the trio.

#### 16.6 [MSc] How do I build the thermal (Gibbs) oscillator state, compute its Wigner function, and evolve a Wigner function by the Moyal equation?

State the Mehler-form Gibbs kernel
$\rho_{\beta}(x,x')=\sqrt{\tfrac{\tanh(\beta/2)}{\pi}}\,\exp\left(-\tfrac{\tanh(\beta/2)}{4}(x+x')^2-\tfrac{\coth(\beta/2)}{4}(x-x')^2\right)$
with $Z=1/(2\sinh(\beta/2))$ and earn it rather than assert it: the Bloch-equation residual
$-\partial_{\beta}\rho_{\beta}=(-\tfrac{1}{2}\partial_x^2+\tfrac{1}{2}x^2)\rho_{\beta}$ must
`FullSimplify` to 0 under $\beta>0$, with unit trace by `Integrate` (whether `Sum` closes the
Fock-ladder Mehler sum symbolically: verify at authoring). Displace by $\alpha$: the kernel picks
up the shift $x\to x-x_0$, $x'\to x'-x_0$ and the phase $e^{ip_0(x-x')}$ with
$(x_0,p_0)=(\sqrt{2}\operatorname{Re}\alpha,\sqrt{2}\operatorname{Im}\alpha)$. One Wigner
`Integrate` under $\beta>0$ then gives the displaced Gaussian
$W(x,p)=\tfrac{1}{\pi\coth(\beta/2)}\,\exp\left(-\tfrac{(x-x_0)^2+(p-p_0)^2}{\coth(\beta/2)}\right)$
of width $\tfrac{1}{2}\coth(\beta/2)$: `Limit` at $\beta\to\infty$ recovers the coherent state's
vacuum width $\tfrac{1}{2}$, a small-$\beta$ `Series` the classical equipartition width $1/\beta$,
and the refuting moment check
$\iint\tfrac{1}{2}(x^2+p^2)\,W\,dx\,dp=\tfrac{1}{2}\coth(\beta/2)+\tfrac{1}{2}(x_0^2+p_0^2)$ must
match the independently known $\langle H\rangle$, reusing the kernel. For the evolution write the
Moyal equation through its first correction,
$\partial_{t}W=-p\,\partial_{x}W+V'(x)\,\partial_{p}W-\tfrac{1}{24}V'''(x)\,\partial_{p}^{3}W$: for
$V=\tfrac{1}{2}x^2$ the correction vanishes identically, so the flow is exactly classical
Liouville, and substituting the rigidly rotating Gaussian
$W_0(x\cos t-p\sin t,\;p\cos t+x\sin t)$ must `FullSimplify` the residual to 0; give the check
teeth by repeating the substitution with $V=x^4/4$, where the retained $V'''$ term leaves a
visibly nonzero residual. Close: a displaced thermal blob orbits the trap rigidly at the
oscillator frequency with no spreading at any temperature; what shears a phase-space distribution
is anharmonicity, never thermal width, and that is the exact content of Moyal reducing to
Liouville for quadratic $H$.

## Part 17 Plan: Continuous-variable quantum optics and information

Seven questions, all MSc. Class census per `Route-Table.md`: C0: 17.1, 17.2, 17.3, 17.4, 17.5,
17.7 (no differential equation; WL machinery named per entry, C0 kernel facts binding); C5: 17.6
(truncated-basis ODE-IVP: `MatrixExp` action form primary, route agreement provably blind to
truncation, $N$-sweep mandatory).

This part lives on one or two bosonic modes, and in 17.6 on one mode times a two-level atom, with
$[a,a^\dagger]=1$ and the quadratures $X=(a+a^\dagger)/\sqrt2$, $P=-i(a-a^\dagger)/\sqrt2$ obeying
$[X,P]=i$ with vacuum variances $\tfrac12$ each, the rotated one being
$X_\theta=(ae^{-i\theta}+a^\dagger e^{i\theta})/\sqrt2$. The Gaussian unitaries are the
displacement $D(\alpha)=e^{\alpha a^\dagger-\alpha^*a}$; the squeeze
$S(\xi)=e^{(\xi^*a^2-\xi a^{\dagger2})/2}$, the convention pinned by 4.7, under which $S(r)$ with
$r>0$ squeezes $X$ to $\langle X^2\rangle=e^{-2r}/2$ and stretches $P$; the beam splitter
$B(\theta)=e^{\theta(a^\dagger b-ab^\dagger)}$; and the two-mode squeezer
$S_2(r)=e^{r(a^\dagger b^\dagger-ab)}$, written creation-minus-annihilation, the opposite overall
sign from the single-mode pattern, which is exactly what makes $X_1-X_2$ and $P_1+P_2$ the
squeezed pair (each scaling by $e^{-r}$) rather than the anti-squeezed one. Where the
beam-splitter phase sits is a genuine convention (verify at authoring): the coincidence null does
not depend on it, the relative sign between $\vert2,0\rangle$ and $\vert0,2\rangle$ does. The atom
enters through the Jaynes-Cummings Hamiltonian
$H=\omega_c a^\dagger a+\tfrac{\omega_a}{2}\sigma_z+g(a\sigma_++a^\dagger\sigma_-)$ with detuning
$\Delta=\omega_a-\omega_c$; in the frame rotating at $\omega_c$ for both mode and atom this becomes
$\tfrac{\Delta}{2}\sigma_z+g(a\sigma_++a^\dagger\sigma_-)$, where the exact two-amplitude solutions
below are read. Phase space carries the Wigner function $W(x,p)$ and the Husimi
$Q(\alpha)=\vert\langle\alpha\vert\psi\rangle\vert^2/\pi$. The working representation throughout is
truncated Fock: `SparseArray` ladder matrices, `KroneckerProduct` for two modes or mode times atom,
`MatrixExp` for the unitaries, with the truncation edge exhibited and swept rather than hidden. The
C0 kernel facts bind (in code `E` is Euler's number, so energies and fields need other names;
`Simplify` and `Integrate` need explicit reality and positivity assumptions such as $r>0$ to close
Gaussian moments).

#### 17.1 [MSc] How do I define the quadrature operators and the optical phase space of a single mode?

Build the ladder matrix $a$ as a `SparseArray` with superdiagonal $\sqrt1,\dots,\sqrt N$ at $N=12$,
form $X=(a+a^\dagger)/\sqrt2$ and $P=-i(a-a^\dagger)/\sqrt2$ with `ConjugateTranspose`, and read the
optical phase space off them: $[X,P]=i$ follows from $[a,a^\dagger]=1$ alone, the vacuum unit vector
gives $\langle X\rangle=\langle P\rangle=0$ and $\langle X^2\rangle=\langle P^2\rangle=\tfrac12$ by
`Dot`, so the vacuum saturates $\Delta X\,\Delta P=\tfrac12$ and occupies a disk of radius
$\sqrt{1/2}$ in the $(x,p)$ plane, and a `Table` of $\langle X_\theta^2\rangle$ over $\theta$ comes
out flat at $\tfrac12$: the vacuum is isotropic, which is what makes $\tfrac12$ one shot-noise
floor rather than a phase-dependent one. On truncated matrices the commutator carries the corner
defect $[X,P]=i(\mathbb 1-(N{+}1)\vert N\rangle\langle N\vert)$, exact in the interior, the same
edge 4.4 exhibits for $\hat x,\hat p$. The refuting cross-check runs on machinery the matrices never
touch: `Integrate` of $x^2\vert\psi_0\vert^2$ for $\psi_0(x)=\pi^{-1/4}e^{-x^2/2}$, and of
$p^2\vert\tilde\psi_0\vert^2$ after `FourierTransform`, must both return $\tfrac12$, where a
misplaced $\sqrt2$ in either quadrature definition would land on $1$ or $\tfrac14$ instead. Close
operationally: $\tfrac12$ is the noise power a balanced detector reads with nothing in the signal
port, and every Gaussian state in this part is this one disk translated, rotated, or deformed at
fixed area.

#### 17.2 [MSc] How do I revisit the displacement and squeeze operators of Part 4 as phase-space transformations and read off their action on the quadratures?

Carry Part 4's parameters $\alpha=2e^{i\pi/3}$ and $\xi=re^{i\pi/3}$ so both phases stay live, and
conjugate by BCH: for $D(\alpha)$ the commutator series terminates after one step, since the
commutator of $X$ with $\alpha a^\dagger-\alpha^*a$ is a c-number, giving
$D^\dagger(\alpha)XD(\alpha)=X+\sqrt2\operatorname{Re}\alpha$ and $P$ shifted by
$\sqrt2\operatorname{Im}\alpha$; for $S(\xi)$ the nested commutators resum to hyperbolic mixing,
the quadratures along the squeeze axes scaling as $e^{\mp r}$ and mixing $X$ with $P$ whenever
$\phi\neq0$, with the overall sign fixed by 4.7's convention rather than left open ($S(r)$ shrinks
$X$). Verify twice, on pictures that must agree: conjugate the truncated matrices of the 17.1
construction by `MatrixExp` of the truncated generators at $N=40$ and compare the leading
$20\times20$ block against the closed forms, demanding `Max@Abs` of the difference below $10^{-8}$
and raising $N$ until it is, with the deviation climbing toward the corner as the diagnostic rather
than a failure; and transform the vacuum Wigner Gaussian, letting `Integrate` under reality
assumptions put the ellipse center at
$(\sqrt2\operatorname{Re}\alpha,\sqrt2\operatorname{Im}\alpha)$ and scale the principal variances by
$e^{\mp2r}$ at fixed area. Close with the shared picture: Gaussian unitaries act affinely on phase
space, displacement translating the vacuum disk, squeezing deforming it at constant area, which is
why neither can produce Wigner negativity.

#### 17.3 [MSc] How do I model a beam splitter on two modes and exhibit Hong-Ou-Mandel interference?

The beam splitter is $B(\theta)=e^{\theta(a^\dagger b-ab^\dagger)}$, realized by `MatrixExp` on the
`KroneckerProduct` two-mode Fock space truncated at two photons per mode; the generator conserves
total photon number, so the two-photon sector
$\{\vert2,0\rangle,\vert1,1\rangle,\vert0,2\rangle\}$ closes exactly and truncation is not an
approximation here. At the balanced point $\theta=\pi/4$ act on $\vert1,1\rangle$: the substitution
this convention produces, $a^\dagger\to(a^\dagger+b^\dagger)/\sqrt2$ and
$b^\dagger\to(b^\dagger-a^\dagger)/\sqrt2$ in $a^\dagger b^\dagger\vert0,0\rangle$, gives
$(\vert0,2\rangle-\vert2,0\rangle)/\sqrt2$ with the $\vert1,1\rangle$ terms cancelling bosonically,
so the coincidence amplitude is exactly zero; the matrix route returns the same state up to an
overall sign, a global phase carrying no physics, so the comparison is on the coincidence weight and
the relative sign, never the global one. Unitarity of $B$ and output probabilities summing to 1 are
the structural gates. Then the pinned distinguishability sweep, the point of the example: give the
second photon the internal mode $s\,e_1+\sqrt{1-s^2}\,e_2$ on a four-mode (two spatial times two
internal) space, derive the closed form $P_{\mathrm{coinc}}(s)=\tfrac12(1-s^2)$, and evaluate the
four-mode matrix computation at intermediate $s$, notably $s=1/\sqrt2$, where the rival forms
$(1-s)/2$, $(1-s^3)/2$ and $(1-s^2)^2/2$ that agree with it at both endpoints separate from it by
tens of percent: the endpoints alone certify nothing. Close in the laboratory: the dip visibility
$s^2$ is the standard meter of single-photon indistinguishability, and a classical field caps the
dip at half the distinguishable plateau, so a deeper dip certifies quantum interference.

#### 17.4 [MSc] How do I build the two-mode squeezed vacuum and exhibit its EPR (position-momentum) correlations?

Apply `MatrixExp` in its action form to the two-mode vacuum with the generator
$r(a^\dagger b^\dagger-ab)$ on a `KroneckerProduct` space at $N=24$ per mode, keeping $r$ symbolic
in the closed forms and $r=1$ for numbers. The state must be the Schmidt-diagonal
$\sqrt{1-\lambda^2}\sum_n\lambda^n\vert n,n\rangle$ with $\lambda=\tanh r$: `ArrayReshape` the
output into its coefficient matrix $c$ and let `SingularValueList` return a geometric sequence of
ratio $\tanh r$, while the reduced state $cc^\dagger$ is thermal with $\bar n=\sinh^2r$, the first
quantification of the entanglement; the second is the EPR pair
$\Delta(X_1-X_2)^2=\Delta(P_1+P_2)^2=e^{-2r}$ against the vacuum value 1, computed with per-mode
quadrature matrices on the same state. That variance identity is also the sign discriminator, and
the only check here that is one: `SingularValueList` cannot see a sign (the flipped generator gives
$(-\lambda)^n$, the same singular values) and neither can the truncation tail, whereas the flipped
generator sends $\Delta(X_1-X_2)^2$ to $e^{+2r}$ and squeezes $X_1+X_2$ instead. Truncation is
checked by an $N$-sweep with its tail predicted in advance, the missing weight
$\lambda^{2(N+1)}\approx10^{-6}$ at $r=1$, $N=24$, and pushing $r$ up at fixed $N$ must visibly
break the variance identity, the honest failure edge. Close on the limit: as $r\to\infty$ the pair
approaches the ideal EPR state, the correlations sharpening without bound while each mode alone
heats toward a featureless thermal state.

#### 17.5 [MSc] How do I detect continuous-variable entanglement with the Duan separability criterion?

Evaluate the Duan combination $\Delta(X_1-X_2)^2+\Delta(P_1+P_2)^2$, bounded below by 2 for every
separable state in the $[X,P]=i$ normalization (bound normalization: verify at authoring). On the
two-mode squeezed vacuum, rebuilt inside the entry from the $r(a^\dagger b^\dagger-ab)$ generator
at $N=24$ so the answer stays self-contained, the Bogoliubov action sends both combinations to
$e^{-r}$ times themselves, so the sum is $2e^{-2r}$, below the bound for every $r>0$; on the
coherent product $D(\alpha)\otimes D(\beta)\vert0,0\rangle$ with $\alpha,\beta$ kept symbolic,
displacement moves means and not variances, so the sum is exactly 2 for all $\alpha,\beta$: one
violating state and one boundary state, both in closed form. The truncated-matrix evaluation must
reproduce both numbers, the $r\to0$ limit of the squeezed sum must land exactly on the
coherent-product value 2, gluing the two computations together, and a sign slip in the two-mode
generator would squeeze the complementary combinations and push the sum to $2e^{+2r}>2$, reporting
no violation at any $r$, so the check discriminates conventions as sharply as it discriminates
states. Close on tightness: a product state sitting exactly on the boundary shows the constant 2
cannot be improved, and this same combination is the witness measured on twin beams in the
laboratory.

#### 17.6 [MSc] How do I solve the Jaynes-Cummings model (one mode plus a two-level atom, truncated Fock) and exhibit vacuum Rabi oscillations?

C5 per `Route-Table.md`. On resonance the initial state $\vert e,0\rangle$ closes into the
two-dimensional block $\{\vert e,0\rangle,\vert g,1\rangle\}$, and `DSolveValue` on the
two-amplitude system returns the exact pair $\{\cos gt,\,-i\sin gt\}$ (the probed C5 gate), so
$P_e(t)=\cos^2gt$: vacuum Rabi oscillation at frequency $2g$ out of a field containing zero
photons. Generalize with detuning through the same `DSolveValue`,
$P_e(t)=1-\frac{4g^2}{4g^2+\Delta^2}\sin^2(\tfrac12\sqrt{4g^2+\Delta^2}\,t)$, oscillating at
$\sqrt{4g^2+\Delta^2}$ with peak transfer probability $4g^2/(4g^2+\Delta^2)$, swept over $\Delta$.
The certified numeric route is the C5 primary: assemble the full truncated Jaynes-Cummings
Hamiltonian with `SparseArray` and `KroneckerProduct` (mode times atom) and evolve with the
`MatrixExp` action form. From $\vert e,0\rangle$ the rotating-wave coupling conserves excitation
number, so the state never leaves the one-excitation manifold, the truncated ladder is never asked
for $a^\dagger\vert N\rangle$, and the $N$-sweep is exactly flat: keep it, but label it the
assembly check it actually is, since all it can detect is an excitation-non-conserving error. The
truncation trap is exercised by a second run that does leave that manifold, an initial coherent
field $\vert g,\alpha\rangle$ with $\vert\alpha\vert^2\approx9$ showing collapse and revival near
$t\approx2\pi\sqrt{\bar n}/g$, where $N$ must exceed $\vert\alpha\vert^2$ plus several
$\vert\alpha\vert$ and the sweep genuinely moves; per the C5 verdict route agreement never
certifies truncation (two routes agreed to $1.1\times10^{-7}$ while both were wrong by
$5.6\times10^{-3}$, a blindness factor $5.1\times10^4$), so only the sweep or an exact benchmark
can. Close in the cavity: the $2g$ oscillation is the single-atom vacuum Rabi frequency of cavity
QED, the revivals are the discreteness of the photon number made visible in an atomic signal, and
the dispersive limit $\Delta\gg g$ shrinks the transfer to $4g^2/\Delta^2$, the regime qubit
readout lives in.

#### 17.7 [MSc] How do I model balanced homodyne and heterodyne detection, measuring a quadrature and reconstructing the Husimi-$Q$ distribution operationally?

Take the squeezed vacuum $S(r)\vert0\rangle$ at $r=1$ on a signal mode truncated at $N_a=40$, where
the neglected weight sits near $10^{-4}$; at $N=12$ it is concentrated at high $n$ and the
anti-squeezed variance falls about 4 percent short of $e^{2r}/2$, ten times the statistical
tolerance quoted below. Derive by `Integrate` (assumptions $r>0$, $\theta$ real) the closed-form
marginal of the state's Wigner Gaussian along $X_\theta$, variance
$\sigma_\theta^2=\tfrac12(e^{-2r}\cos^2\theta+e^{2r}\sin^2\theta)$, sub-vacuum at $\theta=0$ under
4.7's convention, where a flipped squeeze sign would show up loudly as $\langle X^2\rangle=e^{2r}/2$
instead; that marginal is the benchmark, never the sampler. Homodyne is built as apparatus: put a
local oscillator in the coherent state $\beta=\vert\beta\vert e^{i\theta}$ on mode $b$, combine it
with the signal on 17.3's $B(\pi/4)$ by `MatrixExp` in action form, and detect photon number in
both output ports; the difference $n_c-n_d$ pulled back through the splitter is
$n_-=a^\dagger b+ab^\dagger$, whose scaled version $n_-/(\sqrt2\vert\beta\vert)$ tends to $X_\theta$
as the local oscillator grows classical. Sampling is then genuinely operational, `RandomChoice` on
the output Fock weights $\vert\Psi_{\mathrm{out}}\vert^2$ followed by binning
$(n_c-n_d)/(\sqrt2\vert\beta\vert)$, with no appeal to the answer. The refuting check is
quantitative and exact: the scaled variance is $\sigma_\theta^2+\sinh^2r/(2\vert\beta\vert^2)$, so
a sweep over $\vert\beta\vert\in\{2,4,8\}$ (sizing $N_b$ beyond $\vert\beta\vert^2$ plus several
$\vert\beta\vert$, and $N_a$ beyond the $\vert\beta\vert^2/2$ photons each output port then carries)
must show the excess falling as $1/\vert\beta\vert^2$ with that coefficient while the histogram
converges onto the closed marginal, and an independent $N$-sweep at $\theta=\pi/2$, where
truncation bites hardest, must leave the anti-squeezed variance stationary. Heterodyne adds the
vacuum port: splitting the signal against vacuum and reading both quadratures samples
$Q=W\ast W_{\mathrm{vac}}$ with covariance $\sigma_\theta^2+\tfrac12$, cross-checked against
$Q(\alpha)=\vert\langle\alpha\vert\psi\rangle\vert^2/\pi$ evaluated directly from the truncated Fock
amplitudes. Close on the contrast the question asks for: histograms of detector clicks converge to
the formal objects (the Wigner marginal, the $Q$ function), heterodyne adds one full vacuum unit of
variance per quadrature ($\tfrac12$ in these units, half a photon per mode) and so doubles the
vacuum noise, the 3 dB penalty for reading both knobs at once, and the $\theta$-resolved homodyne
variances are the raw data of quantum state tomography.

## Part 18 Plan: Open quantum systems in continuous space

5 questions. Class census per the class-census table in `Route-Table.md`: C5 (truncated-basis
ODE-IVP systems): all of 18.1 through 18.5; 18.1 is the class's probed representative, so entries
cite the C5 gates directly.

### Common ground

Everything in this part rests on the Lindblad master equation
$\dot\rho=-i[H,\rho]+\sum_k\left(L_k\rho L_k^\dagger-\tfrac12\{L_k^\dagger L_k,\rho\}\right)$ for
the reduced density operator of a mode in contact with an environment, and four of its
consequences: trace preservation and complete positivity, which hold for every Lindblad generator
and fail for non-Lindblad approximations such as Caldeira-Leggett; the exact damped-coherent
benchmark, that under $H=\hat n$, $L=\sqrt\gamma\,a$ a coherent state stays coherent with
$\alpha(t)=\alpha e^{-(\gamma/2+i)t}$ and Poisson populations of mean $|\alpha|^2e^{-\gamma t}$;
the equivalence, for quadratic systems, of the operator master equation to a phase-space
drift-diffusion (Fokker-Planck) flow with exactly closing Gaussian moments; and the unravelling
identity, that the deterministic $\rho(t)$ is the ensemble average of stochastic pure-state
trajectories punctuated by quantum jumps.

All five entries share the certified contract of the C5 verdict in `Route-Table.md`: the Fock
space is truncated at $N$, the ladder operator is a `SparseArray` with $a_{n,n+1}=\sqrt{n+1}$, and
the vectorized Liouvillian is assembled as a `SparseArray` of `KroneckerProduct` terms under the
row-major convention $\mathrm{vec}(A\rho B)=(A\otimes B^{T})\,\mathrm{vec}\,\rho$; evolution goes
through the `MatrixExp[L t, v]` action form, exact in time, so truncation is the only error and it
decays with the Poisson tail mass; the independent cross-check is `NDSolveValue` on the
matrix-valued ODE $\rho'(t)=\mathcal L[\rho(t)]$, with the probed caveat that route agreement is
blind to truncation (agreement $1.1\times10^{-7}$ while both routes were wrong by
$5.6\times10^{-3}$ at $N=8$), so the exact benchmark and the $N$-sweep are mandatory wherever they
exist; default-tolerance norm and trace drift grows with integration length, and the positivity
floor at default tolerances is $\sim-10^{-6}$, so every positivity assertion uses a
tolerance-aware threshold.

#### 18.1 [MSc] How do I solve the Lindblad master equation for a damped oscillator?

Damp the coherent state $\alpha=2$ under $H=\hat n$, $L=\sqrt\gamma\,a$, with $\gamma$ symbolic in
the assembly and pinned to $\gamma=\tfrac25$ for the runs: this is the C5 verdict's own probed
representative, and the coherent state is the pointer-like carrier whose stability under damping
is the physics. Assemble the vectorized Liouvillian per the shared contract, evolve by
`MatrixExp[L t, v]`, and cross-check with `NDSolveValue` on the matrix ODE. The refuting check is
the exact benchmark, reusing the stored $a$ and $\rho(t)$: $\mathrm{Tr}[a\,\rho(t)]$ must equal
$\alpha e^{-(\gamma/2+i)t}$ and the diagonal populations must stay Poisson with mean
$|\alpha|^2e^{-\gamma t}$, so the energy $\mathrm{Tr}[\hat n\,\rho(t)]$ decays at exactly
$\gamma$. Because route agreement is blind to truncation, run the $N$-sweep from $N=8$ to $40$
against the benchmark and confirm the super-exponential tail-mass decay the verdict records
(population error $5.4\times10^{-3}$ at $N=8$ falling to $10^{-13}$ by $N=40$), while trace error
and the minimum eigenvalue of $\rho$ (`Eigenvalues`) are monitored at tolerance-aware thresholds.
Close on the pointer reading: damping sends $|\alpha\rangle$ to another coherent state with purity
intact, the stability that 18.3 shows failing spectacularly for superpositions of two of them.

#### 18.2 [MSc] How do I model quantum Brownian motion (Caldeira-Leggett) on a truncated basis?

Evolve a displaced position-squeezed packet ($x_0=2$, squeezing $r=1$, so $\Delta x$ sits well
below the thermal length; the ledger pins "Caldeira-Leggett packet", and this specification sits
inside that class) in the trap $H=\tfrac12(\hat p^2+\hat x^2)$ under the high-temperature
Caldeira-Leggett generator
$\dot\rho=-i[H,\rho]-i\gamma[\hat x,\{\hat p,\rho\}]-2\gamma T[\hat x,[\hat x,\rho]]$ with
$\gamma=\tfrac1{10}$, $T=10$, $N=30$ (prefactor conventions verify at authoring), building
$\hat x=(a+a^\dagger)/\sqrt2$ and $\hat p=i(a^\dagger-a)/\sqrt2$ from the shared `SparseArray`
ladder, vectorizing each commutator term by the row-major convention, evolving by `MatrixExp`
action with the `NDSolveValue` cross-check. The class trap inverts here, per the C5 members-sanity
note: this generator is NOT completely positive, so the transient negativity of $\rho$ is physics,
not solver error, and the positivity monitor must not be used as a correctness gate; instead
exhibit the early-time minimum eigenvalue of $\rho(t)$ dipping far below the $-10^{-6}$ tolerance
floor and let the honest failure regime carry the lesson. The refuting check moves to the moments:
for quadratic $H$ the first and second moments close as a small linear ODE system that
`DSolveValue` solves exactly (damped $\langle p\rangle$, diffusion feeding $\Delta p^2$ toward its
thermal value), and the truncated-basis moments $\mathrm{Tr}[\hat x\,\rho(t)]$,
$\mathrm{Tr}[\hat p\,\rho(t)]$, reusing the stored operators, must land on those closed forms with
the trace pinned at 1 and the $N$-sweep guarding truncation. Close with the repair: adding the
minimal $[\hat p,[\hat p,\rho]]$ diffusion term (coefficient verify at authoring) restores
Lindblad form, and the same monitor that indicted Caldeira-Leggett now certifies positivity,
separating approximation from artifact.

#### 18.3 [MSc] How does dephasing select a pointer basis and decohere a spatial superposition (einselection)?

Damp the superposition $\mathcal N(|\alpha_1\rangle+|\alpha_2\rangle)$ of two coherent states with
$\alpha_1=2$, $\alpha_2=-1$ (an asymmetric pair, distinct from 16.4's even cat) under $H=\hat n$,
$L=\sqrt\gamma\,a$, $\gamma=\tfrac15$, $N\geq30$, built and evolved per the shared contract with
the `NDSolveValue` cross-check and $N$-sweep: the channel that drains energy through $a$ acts on
this state as dephasing in the coherent-state basis it selects. Track two rates on the same run,
reusing the stored $\rho(t)$: the energy $\mathrm{Tr}[\hat n\,\rho(t)]$ relaxing at $\gamma$, and
the off-diagonal coherence weight $|\langle\beta_2(t)|\rho(t)|\beta_1(t)\rangle|$ between the
drifting pointer components $\beta_i(t)=\alpha_i e^{-(\gamma/2+i)t}$, whose suppression factor is
$\exp[-\tfrac12|\alpha_1-\alpha_2|^2(1-e^{-\gamma t})]$ with short-time rate
$\Gamma_{\mathrm{dec}}=\tfrac\gamma2|\alpha_1-\alpha_2|^2$ (closed forms verify at authoring). The
refuting check is the scaling law itself: sweep the separation $d=|\alpha_1-\alpha_2|$ over
$\{\tfrac32,2,3,4\}$ at fixed midpoint with `Table`, extract each rate from the initial slope of
the log-coherence, and the rates must fall on $\tfrac\gamma2 d^2$ while the energy rate stays
pinned at $\gamma$; that quadratic-over-constant contrast is einselection, and a wrong jump of
machinery anywhere would bend the fitted exponent off 2. Close by cashing to the lab: at
macroscopic separation the same law kills a cat's coherence orders of magnitude before its energy
moves, which is why the world looks classical in the basis the damping selects.

#### 18.4 [MSc] How do I cast the optical master equation as a Fokker-Planck equation in the Wigner representation?

Take the optical master equation for a mode in thermal contact, $H=\hat n$ with
$L_1=\sqrt{\gamma(\bar n+1)}\,a$ (emission) and $L_2=\sqrt{\gamma\bar n}\,a^\dagger$ (absorption),
at $\gamma=\tfrac15$, $\bar n=2$. The derivation is symbolic algebra, since per the C5
members-sanity note this member is probed only through the small-system `DSolve` gate and the
Fokker-Planck framing itself is unprobed: apply the Wigner correspondence rules
($a\rho\leftrightarrow(\alpha+\tfrac12\partial_{\bar\alpha})W$ and companions, verify at
authoring) with `D` on a generic $W(x,p,t)$ and `Simplify`, and read off the drift and diffusion
coefficients of the resulting form: drift is the classical rotation from $H$ plus an isotropic
contraction at $\gamma/2$, diffusion is isotropic with $D=\tfrac\gamma2(\bar n+\tfrac12)$ (exact
coefficients verify at authoring). Then close the Gaussian moment ODEs exactly, the probed gate:
$\langle x\rangle$, $\langle p\rangle$ and the three second moments form a small linear system
`DSolveValue` solves in closed form, relaxing to the thermal Gaussian of variance
$\bar n+\tfrac12$, and residual-check that this Gaussian annihilates the derived right-hand side
(`D` plus `Simplify` to 0, reusing the derived drift and diffusion, so a wrong coefficient cannot
survive). The refuting cross-check is the truncated-basis route: the same moments as
$\mathrm{Tr}[\hat x\,\rho(t)]$ and companions from the shared `SparseArray` plus `MatrixExp`
machinery must land on the `DSolveValue` closed forms, with the $N$-sweep guarding truncation.
Close on the reading: the mode performs a classical Ornstein-Uhlenbeck relaxation whose noise
floor $\bar n+\tfrac12$ never reaches zero, the $\tfrac12$ being vacuum fluctuation that survives
at $\bar n=0$.

#### 18.5 [MSc] How do I unravel a damped mode into quantum trajectories (the quantum-jump picture)?

Unravel the damped Fock state $|2\rangle$ under $H=\hat n$, $L=\sqrt\gamma\,a$, $\gamma=\tfrac12$:
the dynamics is confined to $\mathrm{span}\{|0\rangle,|1\rangle,|2\rangle\}$, so truncation is
exact and any ensemble discrepancy indicts the jump logic itself. Run the probed functional loop:
precompute the non-Hermitian step `MatrixExp[-I Heff dt]` with
$H_{\mathrm{eff}}=\hat n-\tfrac{i\gamma}2\,a^\dagger a$, then `FoldList` the unnormalized state
over the time grid, firing a jump (apply $a$, renormalize, redraw the `RandomReal` threshold)
whenever the squared norm crosses the threshold. The probed warning binds: a coherent-state
ensemble test is VACUOUS, $|\alpha\rangle$ being an eigenstate of the jump operator (probed
trajectory spread $3\times10^{-6}$), which is exactly why the ledger pins the Fock start. The
ensemble-vs-master-equation check has teeth here: the mean of $\langle n\rangle(t)$ over roughly
300 trajectories must match $2e^{-\gamma t}$ from the shared vectorized-Liouvillian `MatrixExp`
route within the Monte Carlo error bar; per trajectory, the integer staircase $2\to1\to0$ (a Fock
state stays Fock between jumps; probed at $4.4\times10^{-16}$) and the waiting-time statistics,
the first jump exponential at rate $2\gamma$ because the no-jump norm of $|n\rangle$ decays as
$e^{-n\gamma t}$, with mean jump count $2(1-e^{-\gamma t})$ by photon bookkeeping. Cross-check a
single trajectory on the verified `WhenEvent` plus `DiscreteVariables` machinery inside `NDSolve`
(probed clean), and sweep $dt$ at fixed ensemble since jump-time resolution enters at $O(dt)$, a
flagged open risk. Close on the punchline: each trajectory is what a photodetector records, clicks
and staircases, and the smooth exponential of 18.1 is nothing but their average.

## Part 19 Plan: Path integrals

Part 19, "Path integrals": 6 questions, all MSc.
Class census (`Route-Table.md`): C0 for 19.1, 19.2, 19.3, 19.5 (no differential equation, C0 row); C9 for 19.4, 19.6 (WKB / semiclassical asymptotics, C9 verdict).

Everything in this part rests on the time-sliced propagator
$K(x,x';t)=\langle x\vert e^{-i\hat{H}t}\vert x'\rangle=\lim_{N\to\infty}\int\prod_{j=1}^{N-1}dx_{j}\,(2\pi i\epsilon)^{-N/2}e^{iS_{N}}$
with step $\epsilon=t/N$ and discretized action $S_{N}=\sum_{j}\bigl[(x_{j+1}-x_{j})^{2}/(2\epsilon)-\epsilon V(x_{j})\bigr]$;
the Wick rotation $t\to-i\beta$ turns the weight into $e^{-S_{E}}$ and the trace over periodic paths into the
partition function $Z(\beta)=\operatorname{Tr}e^{-\beta\hat{H}}$; stationary phase $\delta S=0$ selects the classical
path with action $S_{\mathrm{cl}}$, and the Gaussian fluctuations around it give the Van Vleck form
$K_{\mathrm{sc}}=\sqrt{-\partial^{2}S_{\mathrm{cl}}/(\partial x\,\partial x')/(2\pi i)}\;e^{iS_{\mathrm{cl}}}$, exact
whenever the action is quadratic. Every kernel-checkable claim in the part reduces to three facts: a propagator
solves the time-dependent Schrodinger equation (TDSE), collapses to $\delta(x-x')$ as $t\to0^{+}$, and composes,
$\int K(x,y;t_{1})K(y,x';t_{2})\,dy=K(x,x';t_{1}+t_{2})$.

#### 19.1 [MSc] How do I derive the free-particle propagator from the Gaussian path integral?

Slice the free evolution into $N$ steps of width $\epsilon=t/N$, each carrying the short-time kernel
$K_{\epsilon}(x,y)=e^{i(x-y)^{2}/(2\epsilon)}/\sqrt{2\pi i\epsilon}$, and let a single Integrate do the two-slice
composition $\int K_{\epsilon_{1}}(x,y)K_{\epsilon_{2}}(y,x')\,dy$ with symbolic widths under assumptions
$\epsilon_{1}>0$, $\epsilon_{2}>0$, $x,x'$ real (the Fresnel Gaussian stalls or picks a wrong branch without them;
if the oscillatory form still refuses, evaluate the Wick-rotated Gaussian and continue back, recording the swap):
the output must be $K_{\epsilon_{1}+\epsilon_{2}}(x,x')$, which is the induction step. A symbolic Nest of the same
composition at fixed $N=4$ collapses the sliced integral to $K_{4\epsilon}$, and with the induction established the
general $N$-slice form is $K_{N\epsilon}$, so Limit under $\epsilon=t/N$, $N\to\infty$ lands on
$K=e^{i(x-x')^{2}/(2t)}/\sqrt{2\pi i t}$, the limit taken, not asserted. Verify with checks that can refute and
reuse the derived $K$: FullSimplify of the TDSE residual $i\partial_{t}K+\tfrac{1}{2}\partial_{x}^{2}K$ to zero
under $t>0$, and the $t\to0^{+}$ delta limit taken as $\lim_{t\to0^{+}}\int K(x,x';t)f(x')\,dx'=f(x)$ on an explicit
test function through Integrate then Limit, never NIntegrate, which silently returns $0.$ on an emerging point
measure (C0 row). Close by reading the exponent as the classical action of the straight-line path,
$S_{\mathrm{cl}}=(x-x')^{2}/(2t)$: the free kernel is its own stationary-phase form.

#### 19.2 [MSc] How do I get the oscillator propagator (Mehler) from the path integral?

Run the same $N$-slice discretization on the oscillator ($\omega=1$): integrating the $N-1$ interior points of
$S_{N}=\sum_{j}\bigl[(x_{j+1}-x_{j})^{2}/(2\epsilon)-\epsilon x_{j}^{2}/2\bigr]$ leaves a Gaussian whose prefactor
carries the determinant $d_{N-1}$ of the tridiagonal fluctuation matrix, obeying the Gelfand-Yaglom recursion
$d_{j+1}=(2-\epsilon^{2})d_{j}-d_{j-1}$, $d_{0}=1$, $d_{1}=2-\epsilon^{2}$; solve it with RSolve (verify at
authoring; the fallback is the explicit closed form $d_{j}=\sin((j+1)\theta)/\sin\theta$ with
$\cos\theta=1-\epsilon^{2}/2$, certified by FullSimplify of the recursion residual), cross-check the recursion at
fixed $N=3$ against the brute-force sliced Integrate, and take Limit of $\epsilon\,d_{N-1}\to\sin t$ under
$\epsilon=t/N$, $N\to\infty$; the boundary quadratic form closes through cofactors of the same tridiagonal matrix,
landing on the Mehler kernel $K=e^{i[(x^{2}+x'^{2})\cos t-2xx']/(2\sin t)}/\sqrt{2\pi i\sin t}$. Verify by the
refutable pair reusing the derived $K$: the TDSE residual with $V=x^{2}/2$ FullSimplified to zero under $0<t<\pi$,
and the $t\to0^{+}$ delta limit against a test function via Integrate and Limit. The honest edge is the caustic at
$t=\pi$: $\sin t\to0$, the prefactor diverges, and the kernel collapses onto $\delta(x+x')$; state the derived
branch as $0<t<\pi$ and name the Maslov phase beyond it as out of scope. Close on the caustic as the time where
every classical path refocuses regardless of $x'$.

#### 19.3 [MSc] How do I use the imaginary-time path integral to compute the partition function of a finite system?

Compute the oscillator $Z(\beta)$ three ways and make them collide: the spectrum sum
$\sum_{n\geq0}e^{-\beta(n+1/2)}$ closed by Sum; the exact $1/(2\sinh(\beta/2))$, their identity certified by
FullSimplify under $\beta>0$; and the discretized imaginary-time trace as a finite matrix product: grid $[-L,L]$ by
Subdivide, symmetrized Euclidean short-time matrix
$T_{jk}=\sqrt{1/(2\pi\epsilon)}\,e^{-(x_{j}-x_{k})^{2}/(2\epsilon)-\epsilon(V(x_{j})+V(x_{k}))/2}\,\delta x$ with
$V=x^{2}/2$ built by Outer on packed machine arrays (apply N to the grid first so MatrixPower stays packed), then
$Z_{N}=\operatorname{Tr}T^{N}$ via MatrixPower and Tr at $\epsilon=\beta/N$. The refuting check: the deviation of
$Z_{N}$ from the closed form at fixed $\beta$ must fall fourfold per doubling of $N$ (the symmetrized Trotter error
is $O(\epsilon^{2})$); a plateau above that law exposes grid truncation, so sweep $L$ as well. Anchor the classical
limit: Series of $1/(2\sinh(\beta/2))$ at $\beta\to0$ gives $1/\beta$, and the phase-space integral
$\int dx\,dp\,e^{-\beta(p^{2}+x^{2})/2}/(2\pi)$ closed by Integrate gives the same $1/\beta$, equipartition
regained. Close on the two ends of $\beta$: quantum ground-state dominance $e^{-\beta/2}$ against the classical
$1/\beta$.

#### 19.4 [MSc] How do I apply the stationary-phase (semiclassical) approximation and identify the classical action?

State the Van Vleck stationary-phase propagator
$K_{\mathrm{sc}}=\sqrt{-\partial^{2}S_{\mathrm{cl}}/(\partial x\,\partial x')/(2\pi i)}\,e^{iS_{\mathrm{cl}}}$ and
identify the classical action by computing it: DSolve the classical boundary problem $\ddot{x}=-x$, $x(0)=x'$,
$x(t)=x$, then $S_{\mathrm{cl}}=\int_{0}^{t}(\dot{x}^{2}-x^{2})/2\,ds$ by Integrate, closing to
$[(x^{2}+x'^{2})\cos t-2xx']/(2\sin t)$; D gives the mixed derivative $-1/\sin t$, and the assembled kernel is
exact for the oscillator, confirmed by the refutable pair reusing the assembled $K_{\mathrm{sc}}$: FullSimplify of
the TDSE residual to zero under $0<t<\pi$ and the $t\to0^{+}$ delta limit (the C9 verdict rates this member
low-risk pure symbolic Gaussian quadrature). Then locate the boundary of stationary phase: with
$V=x^{2}/2+\lambda x^{4}$, expand the path-integral weight about the classical path; the quartic fluctuation term
contributes at $O(\lambda)$ through the Gaussian moments $\int_{0}^{t}x_{\mathrm{cl}}(s)^{2}G(s,s)\,ds$ and
$\int_{0}^{t}G(s,s)^{2}\,ds$ with $G(s,s)=\sin s\,\sin(t-s)/\sin t$ the Dirichlet fluctuation propagator, each
closed by Integrate and neither vanishing, so the correction enters at exactly $O(\lambda)$ and at no lower order.
Close: for a quadratic action stationary phase is not an approximation at all, and $\lambda x^{4}$ is precisely
where it becomes one.

#### 19.5 [MSc] How do I evaluate a discretized path integral numerically for a simple potential?

Get the quartic-well ground energy from pure matrix multiplication: on a Subdivide grid over $[-L,L]$ build the
symmetrized Euclidean transfer matrix
$T_{jk}=\sqrt{1/(2\pi\epsilon)}\,e^{-(x_{j}-x_{k})^{2}/(2\epsilon)-\epsilon(V(x_{j})+V(x_{k}))/2}\,\delta x$ with
$V=x^{4}/4$ by Outer on packed machine arrays, and read the ground energy from the large-$\beta$ decay of the
trace: $-\epsilon^{-1}\ln\bigl(\operatorname{Tr}T^{N+1}/\operatorname{Tr}T^{N}\bigr)$, computed with MatrixPower
and Tr (repeated Dot of $T$ onto a vector doubles as an internal power-iteration consistency check), flattens to a
plateau once $e^{-\beta(E_{1}-E_{0})}$ is negligible; store the estimate as en0, never E, which is Euler's number
(C0 row). The Trotter step $\epsilon$ is the swept knob: halve it at fixed $\beta$ and the plateau must converge as
$O(\epsilon^{2})$; sweep $L$ too, since a short box lifts $E_{0}$ silently. The refuting benchmark is independent
machinery: the C2-certified FD recipe (SparseArray Band tridiagonal Hamiltonian, Eigenvalues with the Arnoldi
shift, Richardson in the grid step) pins $E_{0}$ far below the transfer-matrix error, and the extrapolated plateau
must land on it within its own $O(\epsilon^{2})$ bar or the discretization is wrong. Close spectrally:
$Z(\beta)=e^{-\beta E_{0}}(1+e^{-\beta(E_{1}-E_{0})}+\dots)$, so the pre-plateau bend of the same data measures the
gap $E_{1}-E_{0}$ as the leading correction.

#### 19.6 [MSc] How do I compute the double-well tunneling splitting from an instanton (imaginary-time bounce)?

Deep-barrier quartic double well $V=(x^{2}-a^{2})^{2}/2$: compute the instanton action
$S_{0}=\int_{-a}^{a}\sqrt{2V(x)}\,dx$ by the C9 quadrature discipline, scaling the turning points out with $x=au$
before Integrate (raw Integrate with symbolic turning-point limits stalls, and the antiderivative-plus-Limit
Newton-Leibniz route silently returns $0$ when the limit direction lands in the forbidden region: the C9 traps
that bite exactly this integral), which closes to $S_{0}=a^{3}\int_{-1}^{1}(1-u^{2})\,du=4a^{3}/3$. The prediction
$\Delta E\propto e^{-S_{0}}$ is asymptotic-only per the C9 verdict (the one-loop prefactor has no kernel
evidence), so never demand equality: measure the splitting $\Delta E=E_{1}-E_{0}$ with the C2-certified FD recipe
(SparseArray Band + Eigenvalues with the Arnoldi shift + Richardson; certified at the $3\times10^{-5}$ splitting
for $a=2.2$; below roughly $10^{-7}$ switch to the parity-resolved half-domain solves the C2 verdict prescribes)
over a sweep of $a$ from about $1.6$ to $2.2$, and compare $\ln\Delta E$ trends: the local slope
$d\ln\Delta E/da$ extracted with Differences must track $-dS_{0}/da=-4a^{2}$ increasingly well as $a$ grows, the
trend that refutes or confirms the exponent while leaving the prefactor alone. Close: $S_{0}$ is the Euclidean
action of the path that crosses the forbidden barrier in imaginary time, a number invisible to every order of
perturbation theory around either minimum.

## Part 20 Plan: Three-dimensional scattering theory

Questions: 6 (20.1 through 20.6). Class census: C4 for 20.1, 20.4, 20.5 (per the `Route-Table.md`
C4 verdict, radial members); C0 for 20.2, 20.3 (the C0 row); C1 for 20.6 (the C1 verdict,
continuum member).

### Common ground

Every question lives in the stationary continuum at energy $E=k^2/2$, where the scattering state
obeys the asymptotics $\psi(\vec r)\to e^{ikz}+f(\theta)\,e^{ikr}/r$ and the amplitude carries all
observables through $d\sigma/d\Omega=|f(\theta)|^2$. A central potential decouples into partial
waves, each radial channel ending in one number, the phase shift defined by
$u_l(r)\to\sin(kr-l\pi/2+\delta_l)$, so that
$f(\theta)=\tfrac1k\sum_l(2l+1)e^{i\delta_l}\sin\delta_l\,P_l(\cos\theta)$ and
$\sigma_{tot}=\tfrac{4\pi}{k^2}\sum_l(2l+1)\sin^2\delta_l$; unitarity compresses into the optical
theorem $\sigma_{tot}=\tfrac{4\pi}{k}\operatorname{Im}f(0)$. Weak scatterers admit the Born
amplitude $f_B(q)=-\tfrac2q\int_0^\infty r\,V(r)\sin(qr)\,dr$ at momentum transfer
$q=2k\sin(\theta/2)$, and the Coulomb $1/r$ tail is the exception that never reaches free
asymptotics: its radial waves are the Coulomb wave functions with logarithmically distorted phases.
Standing inheritances: the C4 radial members use the certified fixed-$E$ recipe (NDSolveValue
outward integration, log-derivative matching to SphericalBesselJ/SphericalBesselY, hard-sphere
anchor $\delta_0=-ka$ extracted to $10^{-9}$ in a clean sweep), matching away from nodes of $u$
where $\beta=u'/u$ diverges (or matching the pair $(u,u')$ linearly) and evaluating interpolants
only at explicit numbers, since trap (e) silently echoes symbolic arguments and extrapolates
out-of-domain; the C0 members name their machinery directly and keep the symbol `E` out of code
(`en` for energies); the C1 member leans only on the probed fact that CoulombF/CoulombG evaluate
in WL 15, with the Rutherford phase and normalization conventions unprobed.

### Per-question entries

#### 20.1 [MSc] How do I expand a scattering state in partial waves and extract the phase shifts $\delta_l$?

Expand the plane wave as $e^{ikz}=\sum_l(2l+1)\,i^l j_l(kr)\,P_l(\cos\theta)$ and let the potential
act channel by channel; the pinned pair is the hard sphere of radius $a$ as the exact anchor with a
finite well of depth $V_0$ and the same radius beside it as the numeric generalization. For the
hard sphere, $u_l(a)=0$ forces the exterior combination $\cos\delta_l\,j_l-\sin\delta_l\,y_l$ to
vanish at the surface, so Solve and FullSimplify give $\tan\delta_l=j_l(ka)/y_l(ka)$ exactly
(SphericalBesselJ, SphericalBesselY), the probed C4 anchor. For the well, run the certified
fixed-$E$ recipe: NDSolveValue outward from the regular start $u_l\sim r^{l+1}$, form the
log-derivative $\beta$ at a matching radius $r_m$ outside the well, and extract
$\tan\delta_l=\frac{k\,j_l'(kr_m)-\beta\,j_l(kr_m)}{k\,y_l'(kr_m)-\beta\,y_l(kr_m)}$; place $r_m$
away from zeros of the interpolant (the log-derivative diverges at nodes of $u$, and matching the
pair $(u,u')$ linearly is the fallback), and evaluate the interpolant only at explicit numeric
points, since trap (e) bites exactly at computed matching points. Verify per the probe: the
PrecisionGoal sweep $\{4,6,8,10\}$ must drop the hard-sphere error monotonically toward $10^{-9}$
against $\delta_0=-ka$, the extraction must not move when $r_m$ does, and the well's numeric
$\delta_0$ must land on its exact region-matching closed form from Solve on the continuity system
(a disagreement refutes the recipe). Close on the low-energy reading: $\delta_0=-ka$ says the
hard sphere's scattering length is its geometric radius, and $\delta_l$ collapses for $l\gg ka$,
the centrifugal barrier screening the distant waves.

#### 20.2 [MSc] How do I assemble the scattering amplitude and the differential cross section?

Recompute the hard-sphere shifts in-entry from $\tan\delta_l=j_l(ka)/y_l(ka)$ (SphericalBesselJ,
SphericalBesselY; each answer is self-contained, nothing imported from 20.1), then assemble
$f(\theta)=\tfrac1k\sum_{l=0}^{l_{\max}}(2l+1)e^{i\delta_l}\sin\delta_l\,P_l(\cos\theta)$ with
LegendreP over a Table of shifts and read off $d\sigma/d\Omega=|f(\theta)|^2$ via ComplexExpand.
Sweep the truncation $l_{\max}$ at $ka\sim1$ until the amplitude stops moving (the shifts die once
$l\gtrsim ka$, and the sweep is the honest convergence exhibit). The refuting check computes
$\sigma_{tot}$ two ways from the same shift table: Integrate of $|f(\theta)|^2$ over the sphere,
which Legendre orthogonality collapses, against the partial-wave sum
$\tfrac{4\pi}{k^2}\sum(2l+1)\sin^2\delta_l$, with the $ka\to0$ end anchored at
$\sigma_{tot}\to4\pi a^2$ (pure $l=0$, four times the geometric cross section). Then push $ka$
high with $l_{\max}$ scaled up and watch the forward peak sharpen while $\sigma_{tot}\to2\pi a^2$;
close on that shadow doubling: the classical $\pi a^2$ is doubled by forward diffraction that a
detector at any finite angle misses.

#### 20.3 [MSc] How do I compute a cross section in the Born approximation?

Take the Yukawa potential $V(r)=-g\,e^{-\mu r}/r$ (screened Coulomb, the Born staple) and evaluate
the central-potential Born amplitude in one Integrate:
$f_B(q)=-\tfrac2q\int_0^\infty r\,V(r)\sin(qr)\,dr=\tfrac{2g}{q^2+\mu^2}$ under Assumptions
$\mu>0$, $q>0$, $g>0$ (the assumptions are load-bearing for the clean closed form, per the C0
row). Substitute $q=2k\sin(\theta/2)$ for $d\sigma/d\Omega=|f_B|^2$, and let a second Integrate
close the total Born cross section $\sigma_B=\tfrac{16\pi g^2}{\mu^2(\mu^2+4k^2)}$ over the solid
angle. The refuting check is the screening-off limit: Limit $\mu\to0$ must turn $d\sigma/d\Omega$
into exactly the Rutherford form $(\eta/2k)^2\csc^4(\theta/2)$ with $\eta=g/k$ (any stray factor
refutes the amplitude convention), while $\sigma_B$ itself diverges in the same limit, the forward
signature of an unscreened $1/r$ tail. Close on the bridge to 20.6: for the Coulomb potential the
Born modulus is not merely a first-order estimate, the exact treatment with Coulomb wave functions
returns the same Rutherford cross section.

#### 20.4 [MSc] How do I extract the low-energy scattering length and verify the optical theorem?

Define $a_s=-\lim_{k\to0}\delta_0/k$ and work the finite well of radius $a$ as a function of depth
$V_0$, class C4. Get $\delta_0(k)$ at a ladder of small $k$ by the certified fixed-$E$ recipe
(NDSolveValue outward, log-derivative matched to SphericalBesselJ/SphericalBesselY away from
interpolant zeros, or the pair $(u,u')$ matched linearly, arguments numeric per trap (e)), and
beside it the exact route: Solve the region-matching continuity system for $\tan\delta_0$ in
closed form and Limit $k\to0$ to $a_s(V_0)=a\left(1-\tan(\kappa_0 a)/(\kappa_0 a)\right)$ with
$\kappa_0=\sqrt{2V_0}$. The refuting structure is the depth sweep: $a_s$ must diverge exactly at
$\kappa_0 a=(2n+1)\pi/2$, the well's known thresholds for a new $l=0$ bound state, and the numeric
extraction must reproduce each pole location (a shifted pole refutes the small-$k$ extraction).
Verify the optical theorem at a moderate $k$ on the same shift set: assemble
$f(0)=\tfrac1k\sum(2l+1)e^{i\delta_l}\sin\delta_l$ and check
$\sigma_{tot}=\tfrac{4\pi}{k}\operatorname{Im}f(0)$ against the partial-wave sum
$\tfrac{4\pi}{k^2}\sum(2l+1)\sin^2\delta_l$; the identity holds term by term because
$\operatorname{Im}(e^{i\delta_l}\sin\delta_l)=\sin^2\delta_l$, so any residual exposes a wrongly
assembled amplitude rather than roundoff. Close on the ultracold reading: near each threshold
$a_s$ swings through $\pm\infty$ and the low-energy cross section $4\pi a_s^2$ dwarfs the
geometric size, the zero-range universality behind Feshbach-tuned atomic gases.

#### 20.5 [MSc] How do I fit a resonance to the Breit-Wigner form and check Levinson's theorem in three dimensions?

Build the well-plus-barrier $V(r)=-V_0$ for $r<a$, $V_b>0$ for $a<r<b$, zero beyond, and choose
$(V_0,V_b,a,b)$ so a single $l=0$ shape resonance sits inside the barrier window $0<E_R<V_b$. The
potential is piecewise constant, so the C4 primary applies region by region (DSolve silently
echoes a Piecewise coefficient): Solve the continuity system for $\tan\delta_0(E)$ in closed form,
tabulate $\delta_0$ on an energy grid through the window with the arctan branch unwrapped so the
phase is continuous, and cross-check a few grid points with the fixed-$E$ radial recipe, matching
away from nodes of $u$ and keeping interpolant arguments numeric (trap (e) bites at computed
matching points). Fit $\delta_0(E)=\delta_{bg}+\arctan\tfrac{\Gamma/2}{E_R-E}$ with
NonlinearModelFit (verify at authoring; flagged unprobed in the C4 verdict), seeded at the
steepest phase rise; the refuting checks are that the fitted $E_R$ must sit on the near-$\pi$ jump
of the tabulated curve, and that a $b$-sweep must shrink $\Gamma$ as the barrier thickens, or the
shape-resonance reading is wrong. Re-run Levinson in three dimensions on the same closed form:
$\delta_0(0)-\delta_0(\infty)=n_b\pi$ with $n_b$ the well's actual bound-state count from its
exact quantization conditions, the resonance contributing its rapid $\pi$ rise at $E_R$ without
changing the count. Close on the lifetime $\tau=1/\Gamma$ growing exponentially with barrier
thickness: the fit parameters carry the decay physics.

#### 20.6 [MSc] How do I treat Coulomb scattering (Rutherford) with the Coulomb wave functions?

Take the repulsive Coulomb potential $V=\alpha/r$ with Sommerfeld parameter $\eta=\alpha/k$, class
C1: continuum physics with no quantization, and the verdict certifies exactly that CoulombF and
CoulombG evaluate in WL 15, while every phase-shift and normalization convention beyond that is
unprobed, so each convention-dependent equality below carries (verify at authoring). Exhibit the
$1/r$ obstruction first: the radial waves never reach free asymptotics, CoulombF with arguments
$l$, $\eta$, $\rho=kr$ oscillating as $\sin(\rho-\eta\ln2\rho-l\pi/2+\sigma_l)$ with the Coulomb
phase $\sigma_l=\arg\Gamma(l+1+i\eta)$ computed from Gamma and Arg. The refuting check reuses
those definitions: evaluate CoulombF at large $\rho$ numerically and compare against the predicted
distorted sinusoid built from the computed $\sigma_l$ (a constant phase offset refutes the
convention, which is precisely the unprobed item). Then assemble the closed Coulomb amplitude
$f_C(\theta)=-\tfrac{\eta}{2k}\csc^2(\theta/2)\,e^{-i\eta\ln\sin^2(\theta/2)+2i\sigma_0}$ (verify
at authoring) and close on the exact Rutherford
$d\sigma/d\Omega=|f_C(\theta)|^2=(\eta/2k)^2\csc^4(\theta/2)$: the logarithmic phases drop out of
the modulus, the quantum result coincides with the classical formula at every angle, and the
forward $\csc^4$ divergence is the unscreened long-range tail that 20.3's Yukawa regulator cuts
off.

## Part 21 Plan: Nonlinear and mean-field wave mechanics

Questions: 4 (21.1 through 21.4). Class census: C6 for all four questions, per the `Route-Table.md`
class census row C6 and the full C6 verdict block; 21.1 is the probed class representative, so
entries cite its gates verbatim.

### Common ground

Every question here lives on the nonlinear Schrodinger equation
$i\,\partial_t\psi = -\tfrac12\,\partial_{xx}\psi + V\psi + g|\psi|^2\psi$, where superposition
fails, so every linear-spectral route is blocked as primary (C6 classification). The stationary
problem is a nonlinear eigenproblem, $\mu\,\phi = -\tfrac12\phi'' + V\phi + g|\phi|^2\phi$ at
fixed norm $\|\phi\|^2=1$, with the chemical potential $\mu$ the Lagrange multiplier of the energy
functional $E[\phi]=\int(\tfrac12|\phi'|^2+V|\phi|^2+\tfrac{g}{2}|\phi|^4)\,dx$ at fixed norm, so
$E = \mu - \tfrac{g}{2}\int|\phi|^4\,dx \ne \mu$ whenever $g>0$: conflating the two is a silent
physics error every entry guards against. Imaginary time turns the equation into a norm-projected
gradient flow whose fixed point is the ground state; the healing length $\xi = 1/\sqrt{2gn}$ at
local density $n$ sets the scale over which the condensate can bend. Standing C6 inheritances
(Route-Table.md, C6 verdict): the primary route is imaginary-time split-step Fourier with per-step
renormalization on the C3-certified Fourier/InverseFourier stack (FourierParameters -> {1,-1},
Nest, packed arrays), with $\mu$ read off the spectral quadratic form; the scheme is first order
in $dt$ in imaginary time because the frozen-nonlinearity substep breaks Strang (probe p3), so
$\mu$ is always $dt$-swept and extrapolated; the independent cross-check is finite differences
plus Newton on the bordered stationary system; SCF linear mixing on the GPE is refuted outright
(probe p4: limit cycles whose eigenvalue plateau fakes convergence); and DSolve does not hand over
soliton closed forms (probe p1: implicit JacobiSN family with `DSolve::ifun`), so closed forms are
certified by symbolic residual.

### Per-question entries

#### 21.1 [MSc] How do I find the ground state of the Gross-Pitaevskii equation by imaginary-time evolution?

Take the pinned trapped repulsive condensate, $V=\tfrac12 x^2$ at $g=100$ with $g=20$ alongside
for the sweep, the probed class representative. Build the imaginary-time split-step evolution on a
Fourier grid of $n=256$ to $512$ points (Fourier/InverseFourier with FourierParameters -> {1,-1},
packed arrays), renormalize after every step, iterate with Nest, and read
$\mu=\langle\phi|H_{GP}[\phi]|\phi\rangle$ off the spectral quadratic form. The scheme is first
order in $dt$ (probe p3: the frozen-nonlinearity substep breaks Strang in imaginary time,
residual-floor ratios 2.000), so a $dt$-halving sweep with extrapolation of $\mu$ is mandatory,
anchored twice: at $g=0$ the solver must return the exact Gaussian ground state with
$\mu-\tfrac12$ at $10^{-15}$, and at large $g$ it must approach the symbolic
$\mu_{TF}=(3g/(4\sqrt2))^{2/3}$ from above. The refuting cross-check is independent machinery end
to end: a SparseArray Band finite-difference Hamiltonian plus Newton on the bordered system
$\{H_0\phi + g\phi^3 - \mu\phi = 0,\ \|\phi\|^2=1\}$ with an ArrayFlatten Jacobian and
LinearSolve, seeded by the Thomas-Fermi profile (probe p4b: agreement $1.2\times10^{-8}$ on $\mu$
after the respective extrapolations). Never SCF linear mixing on the GPE: refuted at mixing
strengths $\alpha\in\{1,0.5,0.08\}$, all limit cycles, and the $\alpha=0.08$ eigenvalue plateau
fakes convergence while $\|\Delta\phi\|$ stalls (probe p4). Close with both numbers the class
demands, $\mu = 14.1343$ and $E = 8.5085$ at $g=100$ since $E=\mu-\tfrac{g}{2}\int|\phi|^4$, and
with the healing length at the trap center against the cloud radius: the condensate is large and
flat except where it must bend.

#### 21.2 [MSc] How do I build the bright and dark solitons of the nonlinear Schrodinger equation?

Give both pinned species their closed forms and certify each by symbolic residual, because DSolve
will not produce them (probe p1: implicit JacobiSN two-constant family with `DSolve::ifun`): the
bright soliton $\psi = A\,\mathrm{sech}(Ax)\,e^{iA^2t/2}$ of the attractive equation and the dark
tanh notch on the repulsive background, each with a FullSimplify residual of exactly 0 (probed).
Then evolve the Galilean-boosted bright soliton, $A=1$, $v=0.5$, out to $T=10$, with the
real-time split-step Fourier stack: in real time $|\psi|$ is invariant in the nonlinear substep,
so Strang stays second order (probe p5: $dt$-halving ratio 4.013, max density error
$6.1\times10^{-8}$ at $dt=10^{-3}$, norm drift $10^{-12}$ over $10^4$ steps), and the refuting
check compares the evolved density against the exact translated profile built from the same
closed-form definition, never a re-typed literal. Trap: the soliton norm is $2A$, not 1, so only
densities and normalized overlaps are honest metrics; an unnormalized overlap sold as a fidelity
is meaningless. The dark soliton's evolution is the verdict's largest unprobed item: the tanh
background carries a $\pi$ phase jump incompatible with a plain periodic grid, so plan the
phase-ramped box or a Dirichlet finite-difference evolution and flag the dark-dynamics cell
authoring-gated. Close on the contrast the pinned pair exists for: attractive self-focusing
balancing dispersion in a localized lump, versus a repulsive phase kink that only exists riding a
finite background.

#### 21.3 [MSc] How do I apply the Thomas-Fermi approximation to a trapped condensate?

Derive the profile symbolically: dropping the kinetic term from the stationary GPE leaves
$|\phi_{TF}(x)|^2 = (\mu - V(x))/g$ where positive, the inverted parabola for $V=\tfrac12 x^2$;
Integrate the norm over the support (the probed $\tfrac{4\sqrt2\,\mu^{3/2}}{3g}$) and Solve the
normalization for $\mu_{TF}=(3g/(4\sqrt2))^{2/3}$, an identity FullSimplify closes to zero (probe
p1). Compare against the full GPE numeric at large $g$: rerun the imaginary-time split-step solver
(redefined inside the answer, per self-containment) at $g=20$ and $g=100$, overlay $|\phi|^2$ on
the inverted parabola, and let the edge healing layer be the visible failure: the TF cloud ends
sharply at $R_{TF}=\sqrt{2\mu_{TF}}$ while the true density bends smoothly through a layer of
width set by the local healing length. The refuting check is quantitative and directional: the
$dt$-extrapolated $\mu$ must approach $\mu_{TF}$ from above, $\mu/\mu_{TF}=1.0094$ at $g=20$ and
$1.0013$ at $g=100$ (probed), and an approach from below or a non-monotone trend refutes either
the numeric or the derivation. Close on that measured approach: the kinetic energy TF discards
lives only in the edge layer, and the ratio's march toward 1 is the approximation earning its
regime.

#### 21.4 [MSc] How do I solve the Hartree mean-field equation as a self-consistent nonlinear Schrodinger problem?

Take the pinned pair of bosons in a harmonic trap with a softened contact interaction, the
normalized Gaussian $V_{int}(u)=\tfrac{g}{\sqrt{2\pi}\,s}\,e^{-u^2/(2s^2)}$ with $g$ and $s$ of
order one (exact values an authoring choice), and solve the Hartree equation
$\varepsilon\,\phi = -\tfrac12\phi'' + \tfrac12 x^2\phi + (V_{int}*|\phi|^2)\,\phi$ by plain
fixed-point iteration: finite-difference Hamiltonian, dense Eigensystem for the lowest orbital,
the mean-field potential by discrete convolution (ListConvolve, alignment conventions checked at
authoring, or the Fourier stack), iterated with NestList. Plain iteration is a contraction for
Hartree (probe p6: contraction 0.135, with $\alpha=0.5$ reaching the same fixed point to
$10^{-12}$), but convergence is monitored on $\|\Delta\phi\|$, never the eigenvalue alone, because
the class refutation showed an eigenvalue plateau masking a limit cycle. Grade the mean field
against the exact two-body solution: KroneckerProduct assembles the 1D pieces into a SparseArray
Hamiltonian on the $(x_1,x_2)$ grid, lowest eigenvalue by the C2-certified sparse eigensolver
machinery (Arnoldi with an explicit shift below the spectrum; the Route-Table hands this reference
off to the linear classes). The refuting pair reuses the converged orbital: the Hartree energy
$E_{H}=2\varepsilon-U$, with the mean-field interaction energy $U$ removed once for double
counting, against the exact $E_0$, and the product state $\phi(x_1)\phi(x_2)$'s normalized overlap
against the exact symmetric ground state. Close with the mean-field error as the honest close:
both gaps grow with $g$, and their growth measures exactly the interparticle correlation a product
ansatz cannot carry.

## Part 22 Plan: Relativistic wave equations

Questions: 5 (22.1 through 22.5). Class census per `Route-Table.md`: C0 (no differential
equation) for 22.1, 22.2, 22.3; C8 (coupled singular radial system, sole member) for 22.4, per
the full C8 verdict as amended by revision R1 in its log; C4 (fixed-energy scattering BVP) for
22.5, per the C4 verdict and its members-sanity note.

### Common ground

The part works in $\hbar=c=1$ with the mass $m$ kept explicit, because the rest energy $m$ is the
protagonist: every pathology and every limit in these five questions happens on the scale
$E\sim m$. The Klein-Gordon equation $(\partial_t^2-\partial_x^2+m^2)\phi=0$ quantizes the
relativistic dispersion $E^2=k^2+m^2$, so both roots $E=\pm\sqrt{k^2+m^2}$ solve it, and its
conserved density $\rho=\tfrac{i}{2m}(\bar\phi\,\partial_t\phi-\phi\,\partial_t\bar\phi)$ is not
sign-definite. The Dirac equation $(i\gamma^\mu\partial_\mu-m)\psi=0$ rests on the Clifford
algebra $\{\gamma^\mu,\gamma^\nu\}=2\eta^{\mu\nu}$ with
$\eta=\operatorname{diag}(1,-1,-1,-1)$; it buys the positive density $\rho=\psi^\dagger\psi$ and
current $j^\mu=\bar\psi\gamma^\mu\psi$ ($\bar\psi=\psi^\dagger\gamma^0$) at the price of keeping
the negative-energy branch, which returns as the Klein paradox and zitterbewegung. Its exact
hydrogen spectrum $E_{nj}$ is the fine structure that 14.2 obtains perturbatively. The notation
is crowded, so code names carry the distinctions the TeX cannot: the gamma matrices are
`gmat[mu]` with `alphaD` and `betaD` for the $1+1$ pair, while the fine-structure combination is
`zalpha` ($Z\alpha$) and the Dirac-Coulomb exponent is `gam` ($\gamma$); energies are `en`
(`E` is Euler's number); the step potential is `HeavisideTheta` while 22.1's mixing angle is
`th`. Standing kernel facts binding all five entries (Route-Table.md, C0 row and C8 verdict):
densities and currents are written as explicit $\phi\bar\phi$ products before differentiation,
never through Abs; FullSimplify closes residuals only under explicit assumptions ($m>0$, reality
of $x,t,k$, $0<Z\alpha<1$); and every DSolve on a coupled or piecewise system is time-boxed,
because the failure modes in this part are silent hangs and silent echoes, not messages.
Concern: PIPELINE fixes $\hbar=m=1$, and this part departs from it, keeping $\hbar=c=1$ with $m$
symbolic, because $m=1$ would erase the rest energy that every question here measures against
($V_0>2m$, $g=2$, $E_{nj}\to m$, the $1/(2m)$ tremor); the departure is deliberate and part-wide,
and the coordinator should confirm it before authoring.

### Per-question entries

#### 22.1 [MSc] How do I solve the Klein-Gordon equation, write its plane-wave modes, and confront the negative-energy and indefinite-density problem?

Write both branches $\phi_{\pm}=e^{i(kx\mp\omega t)}$ with $\omega=\sqrt{k^2+m^2}$ and earn them
twice: the residual of $(\partial_t^2-\partial_x^2+m^2)\phi_{\pm}$ FullSimplifies to 0, and Solve
on the dispersion polynomial returns both roots $\pm\omega$, the negative branch arriving
uninvited. The pinned example is a two-mode superposition with negative-energy admixture, a
positive-energy rest mode plus a negative-energy moving mode,
$\phi=\cos\theta\,e^{-imt}+\sin\theta\,e^{i(\sqrt3\,mx+2mt)}$ at $\theta=\pi/6$. Build
$\rho=\tfrac{i}{2m}(\bar\phi\,\partial_t\phi-\phi\,\partial_t\bar\phi)$ and
$j=\tfrac{1}{2im}(\bar\phi\,\partial_x\phi-\phi\,\partial_x\bar\phi)$ with the conjugate written
out and reduced by ComplexExpand under reality of $x,t,\theta$ and $m>0$ (never differentiate an
Abs square), then FullSimplify the continuity $\partial_t\rho+\partial_x j$ to 0, a check any
sign slip refutes. Compute the indefiniteness instead of narrating it: $\rho$ reduces to
$\cos^2\theta-2\sin^2\theta-\sin\theta\cos\theta\,\cos\Theta$ with travelling phase
$\Theta=3mt+\sqrt3\,mx$, so Minimize over one period gives the closed-form minimum
$(1-\sqrt3)/4<0$ at the pinned angle, and one Solve of
$\cos^2\theta-2\sin^2\theta-\sin\theta\cos\theta=0$ widens that single case into the general
threshold $\tan\theta_c=\tfrac12$, above which $\rho$ goes negative somewhere, while each branch
alone gives the constant $\pm\omega/m$. Two plane waves are not square integrable, so this is a
pointwise statement; the normalizable version, worth one Integrate on a mixed-branch $L^2$
profile, is that the total $\int\rho\,dx$ itself can be negative. Anchor with the nonrelativistic
limit: substitute $\phi=e^{-imt}\psi(x,t)$ with generic $\psi$ and Series in $1/m$ to recover
$\rho=|\psi|^2+O(1/m)$, the Schrodinger density. Close on the reinterpretation: $\rho$ survives
as a charge density, and the sign it refuses to fix belongs to the antiparticle.

#### 22.2 [MSc] How do I build the Dirac equation from gamma matrices and write the free spinor solutions?

Construct the Dirac-basis gammas as explicit $4\times4$ matrices, $\gamma^0=\sigma_3\otimes\mathbb 1$
and $\gamma^i=i\sigma_2\otimes\sigma_i$ via KroneckerProduct and PauliMatrix with
$\eta=$ DiagonalMatrix of $\{1,-1,-1,-1\}$, and verify the Clifford algebra exhaustively: Table
over all sixteen $(\mu,\nu)$ pairs of $\gamma^\mu\gamma^\nu+\gamma^\nu\gamma^\mu-2\eta^{\mu\nu}\mathbb 1$,
Flatten and Union collapsing to $\{0\}$, an exact matrix computation with nothing sampled. The
pinned example keeps the momentum nontrivial: $\vec p=p\,(\sin\vartheta,0,\cos\vartheta)$
symbolic in the $x$-$z$ plane, helicity two-spinors $\chi_{\pm}$ in half-angle form checked as
eigenvectors of $\vec\sigma\cdot\hat p$, then
$u_{\pm}=\sqrt{E+m}\,\bigl(\chi_{\pm},\pm\tfrac{p}{E+m}\chi_{\pm}\bigr)^T$ and $v_{\pm}$ with the
blocks swapped, all simplification under $p>0$, $m>0$, $0<\vartheta<\pi$ with
$E=\sqrt{p^2+m^2}$ (`en` in code). Verify per spinor, reusing the defined gammas: the Dirac
residuals $(\gamma^\mu p_\mu-m)u_{\pm}=0$ and $(\gamma^\mu p_\mu+m)v_{\pm}=0$ by FullSimplify,
the normalizations $\bar uu=2m$ and $\bar vv=-2m$ with $\bar u=u^\dagger\gamma^0$ via
ConjugateTranspose, and helicity $\vec\Sigma\cdot\hat p\,u_{\pm}=\pm u_{\pm}$ with
$\Sigma_i=\mathbb 1\otimes\sigma_i$; each equation fails under any wrong sign or normalization.
Take both limits with Series on those same expressions: at $p\ll m$, $u_{\pm}\to\sqrt{2m}
(\chi_{\pm},0)^T$ with the lower block $O(p/m)$, which is precisely the small component 22.3
eliminates, and at $m\to0$ the normalization $\bar uu=2m$ degenerates to 0 while
$\vec\Sigma\cdot\hat p$ and $\gamma^5$ agree on $u_{\pm}$, helicity collapsing onto chirality.
The genuine edge is $p\to0$ exactly, where $\hat p$ is undefined and with it the helicity
spinors: at rest only the spin projection survives, and the $\pm$ labels stop meaning anything.
Close with completeness: $\sum_s u_s\bar u_s-(\gamma^\mu p_\mu+m)$ FullSimplifies to the zero
matrix, tying the spinors back to the algebra that built them.

#### 22.3 [MSc] How do I take the nonrelativistic limit to the Pauli equation and recover the $g=2$ prediction?

Couple the Dirac equation minimally, $\vec p\to\vec\pi=\vec p-q\vec A$ with a fully symbolic
$\vec A(x,y,z)$ and scalar $\Phi$, charge-sign convention tagged (verify at authoring), and let
every operator act on a generic two-spinor test function $f$ through D, so the noncommutativity
of $\vec p$ and $\vec A$ is handled by differentiation, never by symbol shuffling. First earn the
key identity on generic $f$: $(\vec\sigma\cdot\vec\pi)^2f-(\vec\pi^2-q\,\vec\sigma\cdot\vec B)f$
FullSimplifies to $\{0,0\}$ with $\vec B=$ Curl of $\vec A$ and PauliMatrix carrying the algebra.
The identity fixes the coefficient of $\vec\sigma\cdot\vec B$, nothing more; $g$ enters one step
later, when $-\tfrac{q}{2m}\vec\sigma\cdot\vec B$ is read against $-\vec\mu\cdot\vec B$ with
$\vec\mu=g\,\tfrac{q}{2m}\vec S$ and $\vec S=\tfrac12\vec\sigma$, giving $g=2$. Then reduce:
split $\psi=e^{-imt}(\varphi,\chi)^T$, stripping the rest phase before any expansion (Series in
$1/m$ chokes on the oscillatory factor otherwise), solve the lower component as
$\chi=\tfrac{\vec\sigma\cdot\vec\pi}{2m}\varphi+O(1/m^2)$, substitute into the upper equation,
and a Series at $m\to\infty$ lands on
$i\partial_t\varphi=[\tfrac{\vec\pi^2}{2m}-\tfrac{q}{2m}\vec\sigma\cdot\vec B+q\Phi]\varphi$.
Coefficient on that expansion confirms the $\vec\sigma\cdot\vec B$ term, but it and the identity
share the elimination $\chi=\tfrac{\vec\sigma\cdot\vec\pi}{2m}\varphi$, so a factor-of-2 slip
there moves $g$ and both readings still pass; the independent refuter is the exact Dirac Landau
spectrum in uniform $\vec B=B\hat z$, obtained by squaring the Dirac operator with the same
$\vec A$ and Solve, $E^2=m^2+p_z^2+qB(2n+1-\sigma)$ with $\sigma=\pm1$, whose Series at
$m\to\infty$ returns $E\approx m+\tfrac{p_z^2}{2m}+\tfrac{qB}{2m}(2n+1-\sigma)$ with no
elimination assumed anywhere, and whose $(n,\sigma=+1)$ and $(n+1,\sigma=-1)$ coincidence (the
degeneracy 14.5 exhibits) happens only at $g=2$. Close in the lab: the measured electron value
2.00232 sits a tenth of a percent above the Dirac prediction, and that excess belongs to
radiative corrections outside this equation.

#### 22.4 [MSc] How do I obtain the Dirac fine-structure spectrum of hydrogen?

State the coupled radial system $F'=-\tfrac{\kappa}{r}F+(E+m+\tfrac{Z\alpha}{r})G$,
$G'=\tfrac{\kappa}{r}G-(E-m+\tfrac{Z\alpha}{r})F$ (the convention the C8 probes
residual-verified) and never hand it to DSolve: the coupled system hangs even without boundary
conditions, and naive elimination returns only DifferentialRoot (C8 gates; time-box any DSolve at
authoring). Fix the quantum-number dictionary first, because it carries the structure of the
answer: $\kappa=-(l+1)$ for $j=l+\tfrac12$ and $\kappa=+l$ for $j=l-\tfrac12$, hence
$|\kappa|=j+\tfrac12$ and $\gamma=\sqrt{(j+\tfrac12)^2-(Z\alpha)^2}$ (code `gam`, `zalpha`), the
identity from which the spectrum's independence of $l$ follows rather than being asserted. The
probed backbone is the squared Biedenharn route: DSolve the effective Coulomb ODE
$-u''+(\tfrac{s(s+1)}{r^2}-\tfrac{2EZ\alpha}{r})u=(E^2-m^2)u$ with non-integer
$s\in\{\gamma-1,\gamma\}$, which returns WhittakerM with symbolic energy (`en`, and
$0<Z\alpha<1$, $m>0$ explicit in every Simplify); quantization is the manual termination
read-off $EZ\alpha/\sqrt{m^2-E^2}=n_r+s+1$, and Solve gives the exact
$E_{nj}=m\,[1+(\tfrac{Z\alpha}{n_r+\gamma})^2]^{-1/2}$ with $n=n_r+|\kappa|$, the $s=\gamma-1$
and $s=\gamma$ branches paired by hand and $n_r=0$ existing only for $\kappa<0$. Certify every
exhibited state by substituting the assembled $(F,G)$ into the first-order system and
FullSimplifying the residual to $\{0,0\}$ (the probes certify general-$n$ closed forms only for
the ground state, so the residual runs per state), and assert the normalization
$\int_0^\infty(F^2+G^2)\,dr=1$ by Gamma integrals. Cross-check against 14.2's perturbative fine
structure: Series of $E_{nj}$ in $Z\alpha$ reproduces
$m\,[1-\tfrac{(Z\alpha)^2}{2n^2}-\tfrac{(Z\alpha)^4}{2n^4}(\tfrac{n}{j+1/2}-\tfrac34)]$ exactly,
remainder scaling as $(Z\alpha)^6$ (probed). The numeric refuter is 32-digit NDSolve shooting on
the coupled system from $r_0=10^{-8}$ with the $E$-independent indicial ratio
$G/F=(\gamma+\kappa)/(Z\alpha)$, the balance of the $r^{\gamma-1}$ terms, which reduces to the
probed $(\gamma-1)/(Z\alpha)$ only at $\kappa=-1$; seeding a $\kappa=+1$ state with that special
value picks the irregular branch and silently disarms the discriminator, whereas the general
seed blows the far tail up by a factor near $10^7$ for a $10^{-2}$ energy error (probed), with
node counting and $r_{\max}\sim1/\sqrt{m^2-E^2}$ for excited states. Exhibit the degeneracy at
fixed $n$ as $\kappa\to-\kappa$: $2s_{1/2}$ ($\kappa=-1$) and $2p_{1/2}$ ($\kappa=+1$) coincide
exactly, since $E_{nj}$ sees $\kappa$ only through $|\kappa|=j+\tfrac12$, while $2p_{3/2}$
($\kappa=-2$) stays split off by $\alpha^4m/32$ in the same Series. The edges are where the
formula ends: as $Z\alpha\to|\kappa|^{-}$ the exponent $\gamma\to0$ and Limit on
$E_{1,1/2}=m\gamma$ sends the ground level to 0 with divergent slope near $Z\approx137$, the fall
to the center; and at the origin the components go as $r^\gamma$ with $\gamma<1$, finite in
amplitude but with divergent derivative, still normalizable since $2\gamma>-1$, against the
Schrodinger $1s$ reduced radial that leaves the origin linearly as $r$. Close on that cusp and
its missing partner: the relativistic origin behavior is qualitative, and the splitting the
Dirac equation refuses to produce, $2s_{1/2}$ against $2p_{1/2}$, is the Lamb shift.

#### 22.5 [MSc] How do I exhibit the Klein paradox and zitterbewegung numerically?

Work in the $1+1$ Dirac representation $\alpha=\sigma_1$, $\beta=\sigma_3$ (code `alphaD`,
`betaD`; representation choice verify at authoring), so spinors are two-component. Klein step:
the pinned supercritical step $V=V_0\,$HeavisideTheta$(x)$ with $V_0>E+m$ (hence $V_0>2m$),
matched region by region because DSolve on a Piecewise coefficient silently echoes (C4 verdict;
the probes covered the Schrodinger case only, so time-box the Dirac matching at authoring): the
equation is first order, so only $\psi$ itself is continuous at $x=0$, a two-component Solve, and
the transmitted momentum $q=\sqrt{(E-V_0)^2-m^2}$, real again in the Klein zone, has its sign
fixed by the group velocity $dE/dq$ (verify at authoring). Transmission comes from the current
ratio $j=\psi^\dagger\alpha\psi$, never from $|t|^2$ (the C4 members-sanity note); FullSimplify
$R+T-1$ to 0 from the currents, then sweep $T(V_0)$ at fixed $E$ across all three windows:
ordinary transmission for $V_0<E-m$, a closed window of width exactly $2m$ for
$E-m<V_0<E+m$ where $q$ is imaginary and $T=0$, and the Klein zone $V_0>E+m$ where $T$ reopens
and stays finite however high the step is raised. The same sweep on the Schrodinger step decays
monotonically and never reopens, so the discriminator is the reopening itself, and the closed
window being exactly $2m$ wide is the pair-creation threshold in disguise. Zitterbewegung is a
C0 mode sum (Route-Table C0 row), no PDE solver: pin
$\psi_0(x)=\operatorname{sech}^2(x/a)\,e^{ik_0x}\,(1,0)^T$ with $k_0=m/2$ and $a=4/m$, the flat
spinor chosen because it is not a $\pm$ eigenvector, and report $\|\Lambda_{-}(k_0)\psi_0\|^2$
as the number that sets the tremor amplitude (of order $1/(2m)$, half the Compton wavelength);
starting from $u_{+}(k_0)$ instead would send it to zero and leave the checks with nothing to
measure. The momentum amplitudes stay symbolic through
$\int\operatorname{sech}^2(x)e^{-ikx}dx=\pi k/\sinh(\pi k/2)$ by Integrate, rescaled for width
$a$, so the only discretization is the $k$-quadrature (Subdivide); split each mode with the
projectors $\Lambda_{\pm}(k)=\tfrac12(\mathbb 1\pm H(k)/E_k)$, $H(k)=\alpha k+\beta m$, rebuild
$\psi(x,t)$ by Total over the grid, and take $\langle x\rangle(t)$ as a $dx$-weighted grid sum:
a drift with a superposed tremor. Three refuters reuse those same definitions: reconstructing
$\psi(x,0)$ from the discrete mode sum must reproduce the pinned profile and keep improving
under a $k$-grid doubling, or the quadrature is inventing the signal; the tremor period must
match $\pi/E_{k_0}=2.81/m$; and rerunning the identical sum on $\Lambda_{+}$-projected data must
flatten the tremor while leaving the drift. State the validity window with the period check,
since it is the reason $a$ is pinned wide: $a=4/m$ gives $\Delta k\approx m/4$, hence
$\Delta(2E_k)\approx0.22m$ and dephasing near $4.5/m$, about one and a half clean oscillations
(and $a=8/m$ doubles that), while at the literal Compton width $a=1/m$ the tremor dephases in
$1.1/m$ and the period is unmeasurable; the drift over one period, $\approx1.26/m$ at
$k_0=m/2$, stays comparable to the tremor amplitude, so both live on one plot. Close where the
two computations meet: both effects are $\pm E$ interference, and both dissolve only in the
many-particle reinterpretation the supercritical step forces.

## Part 23 Plan: From one particle to fields: the second-quantization bridge

Questions: 6 (23.1 through 23.6), all MSc. Class census: C0 for 23.1 through 23.5 (no differential
equation, per the `Route-Table.md` C0 row); C5 for 23.6 (time-dependent generator, per the C5
verdict block).

### Common ground

Everything in this part is one move repeated: promote the mode amplitudes of a one-particle problem
to ladder operators. The single ladder obeys $[a,a^{\dagger}]=1$ with number operator
$\hat n=a^{\dagger}a$; the multimode Fock space is the tensor product of per-mode ladders, with
occupation basis $\vert n_1,n_2,\ldots\rangle$ and operators of different modes commuting; the
field operator is the mode sum $\hat\psi(x)=\sum_n a_n\varphi_n(x)$ over an orthonormal set, with
$[\hat\psi(x),\hat\psi^{\dagger}(x')]=\delta(x-x')$ holding only in the untruncated limit; the
second-quantized Hamiltonian is
$\hat H=\sum_{ij}h_{ij}\,a_i^{\dagger}a_j+\tfrac12\sum_{ijkl}V_{ijkl}\,a_i^{\dagger}a_j^{\dagger}a_l a_k$
with $h_{ij}$ and $V_{ijkl}$ first-quantized matrix elements; the vacuum energy is the mode sum
$E_0=\tfrac12\sum_n\omega_n$; and a time-dependent quadratic Hamiltonian mixes $a$ with
$a^{\dagger}$ through a Bogoliubov transformation $b=\alpha a+\beta^{*}a^{\dagger}$ constrained by
$\vert\alpha\vert^{2}-\vert\beta\vert^{2}=1$, whose $\vert\beta\vert^{2}$ counts particles created
from vacuum. Every construction lives on an explicitly truncated space (a per-mode quantum cap, a
mode-number cap $N$, a site count, a frequency cutoff), and each answer states what its truncation
costs and where the cap is exact. Standing C0 inheritances (`Route-Table.md`, C0 row): truncated
matrices via SparseArray and KroneckerProduct with Eigenvalues/Eigensystem; delta-function matrix
elements go through Integrate only, never NIntegrate; `E` is Euler's number, energies are `en`;
Simplify needs explicit positivity, reality, and integer assumptions to close sums and residuals.

### Per-question entries

#### 23.1 [MSc] How do I build multimode occupation-number (Fock) space from copies of the oscillator ladder?

Build one truncated ladder $a$ as a SparseArray with $\sqrt{1},\ldots,\sqrt{n_{\max}}$ on the
superdiagonal at $n_{\max}=4$ (five levels per mode), then embed three copies as
$a_1,a_2,a_3$ by KroneckerProduct with `IdentityMatrix[5, SparseArray]` into the $125$-dimensional
three-mode space, whose product basis is the occupation basis $\vert n_1,n_2,n_3\rangle$. Verify
the algebra where it can refute the construction: each per-mode commutator equals the truncated
identity with its corner defect, $[a_j,a_j^{\dagger}]=\mathbb 1-(n_{\max}{+}1)\,P_j$ with $P_j$
projecting mode $j$ onto its top level, exactly (the defect 4.4 exhibits for $[\hat x,\hat p]$;
demanding the ideal $\mathbb 1$ would be a false test on any finite matrix), while every cross-mode
commutator $[a_i,a_j^{\dagger}]$ with $i\ne j$ is exactly the zero matrix, defect-free. Then do
occupation arithmetic: the total number operator $\hat N=\sum_j a_j^{\dagger}a_j$ must come out
diagonal with entries $n_1+n_2+n_3$ read against the product-basis index arithmetic, and
$a_2^{\dagger}a_1$ applied to the explicit unit vector for $\vert 1,0,2\rangle$ must return exactly
the unit vector for $\vert 0,1,2\rangle$, a check a wrong tensor-slot order fails loudly. Close on
what the checks just earned: cross-mode commutativity is exact even when truncated because it is
tensor-product structure, while each mode's own ladder carries the finite-dimensional corner defect.

#### 23.2 [MSc] How do I define the field operators $\hat\psi(x)$ as mode sums and impose the (anti)commutation relations?

Assemble $\hat\psi(x)=\sum_{n=1}^{N}a_n\varphi_n(x)$ on the first $N$ infinite-well modes
$\varphi_n(x)=\sqrt{2/L}\,\sin(n\pi x/L)$, and impose $[a_n,a_m^{\dagger}]=\delta_{nm}$
symbolically with KroneckerDelta so that Sum reduces the operator commutator to the c-number kernel
$K_N(x,x')=\sum_{n=1}^{N}\varphi_n(x)\varphi_n(x')$: the content of the question is what this
kernel does as $N$ grows. Get the exact partial sum in closed form, a Dirichlet-type kernel: the
product-to-sum split writes $K_N$ as a difference of Dirichlet kernels in $x-x'$ and $x+x'$, and
Sum with Simplify under $0<x<L$, $0<x'<L$ should close it (verify at authoring). The refuting
check is weak convergence, never a pointwise delta claim: smear against the smooth test function
$f(x')=x'(L-x')$, whose exact Integrate against $K_N$ at the pinned interior point $x=3L/10$ must
converge to $f(3L/10)$ as $N$ doubles through $4,8,16,32$ (Fourier coefficients of $f$ decay as
$n^{-3}$, so the error should fall roughly as $N^{-2}$), while the diagonal $K_N(x,x)$ grows
without bound like $N/L$, exhibited beside it. For fermions the identical c-number kernel follows
from $\{c_n,c_m^{\dagger}\}=\delta_{nm}$, one KroneckerDelta line, so the statistics are invisible
at this level. Close on the statement the pair of computations earns: $\hat\psi(x)$ is an
operator-valued distribution, $K_N\to\delta(x-x')$ only against test functions, and the divergent
diagonal is that same fact seen pointwise.

#### 23.3 [MSc] How do I represent a free scalar field as infinitely many oscillators and read off the vacuum energy?

Massless scalar field in a box of length $L$ with Dirichlet walls: mode frequencies
$\omega_n=n\pi/L$, vacuum energy $E_0=\tfrac12\sum_n\omega_n$. Exhibit the divergence with an
explicit sharp cutoff first: Sum gives $E_0(n_{\max})=\tfrac{\pi}{4L}n_{\max}(n_{\max}{+}1)$,
growing as the square of the cutoff. Then regulate smoothly with $e^{-\omega_n/\Lambda}$: Sum
closes geometrically to $\tfrac{\pi}{2L}\,e^{\pi/(L\Lambda)}/(e^{\pi/(L\Lambda)}-1)^{2}$ (verify at
authoring), and Series in $1/\Lambda$ splits it as $\tfrac{L\Lambda^{2}}{2\pi}-\tfrac{\pi}{24L}
+O(\Lambda^{-2})$: a bulk divergence proportional to $L$ plus a finite geometry-dependent term. The
Casimir-flavored difference is two geometries at one cutoff: partition a box of fixed total length
$D$ at position $L$ and form $\Delta E(L)=E_0(L)+E_0(D{-}L)-2E_0(D/2)$, reusing the same regulated
summand; the bulk terms cancel identically and Limit as $\Lambda\to\infty$ leaves
$-\tfrac{\pi}{24}\left(\tfrac1L+\tfrac1{D-L}-\tfrac4D\right)$. State plainly what is and is not
claimed: this is mode-sum regularization of a one-dimensional scalar toy, not the electromagnetic
Casimir experiment; the claim is only that the difference of two divergent sums has a finite,
cutoff-independent limit. The refuting check is exactly that independence: sweep $\Lambda$ over a
decade numerically with the same summand and watch $\Delta E$ settle on the Series value; residual
$\Lambda$ dependence refutes the regularization. Close with the force: $-\partial_L\Delta E$ pulls
the partition toward the nearer wall, approaching the attraction $-\pi/(24L^{2})$ as $D\to\infty$.

#### 23.4 [MSc] How do I show the second-quantized many-body Hamiltonian equals the symmetrized first-quantized one?

Take the first $M=3$ well modes on $(0,1)$, one-body part diagonal with $E_i=i^{2}\pi^{2}/2$, and
the contact pair interaction $g\,\delta(x_1-x_2)$, whose tensor
$V_{ijkl}=g\int_0^1\varphi_i\varphi_j\varphi_k\varphi_l\,dx$ is exact rationals times $g$ by
Integrate only (NIntegrate returns 0. on a point measure, the C0 binding fact). First-quantized
side: the six symmetrized pair states
$\vert ij\rangle_S=(\vert ij\rangle+\vert ji\rangle)/\sqrt{2(1+\delta_{ij})}$ and the $6\times6$
matrix of $h\otimes\mathbb 1+\mathbb 1\otimes h+g\,\delta(x_1-x_2)$ by exact Integrate, the delta
collapsing the double integral to one dimension. Second-quantized side: rebuild per-mode ladders
capped at two quanta (a $27$-dimensional space) via SparseArray and KroneckerProduct, assemble
$\hat H=\sum_i E_i\,a_i^{\dagger}a_i+\tfrac12\sum_{ijkl}V_{ijkl}\,a_i^{\dagger}a_j^{\dagger}a_l a_k$,
and restrict to the total-number-two eigenspace of $\hat N$ in the normalized pair basis
$a_i^{\dagger}a_j^{\dagger}\vert 0\rangle/\sqrt{1+\delta_{ij}}$; the two-quanta cap represents the
two-boson sector exactly, since normal-ordered $\hat H$ never crosses the truncation edge inside
the sector, so equality can be demanded exactly rather than approximately. The interacting term is
the content: the free parts agree trivially, so the refuting check is Simplify of the entrywise
difference of the two $6\times6$ matrices to the zero matrix with $g$ symbolic, where a wrong
$\tfrac12$, a dropped exchange term, or an index-order slip in $V_{ijkl}$ shifts the pair-changing
off-diagonal entries. Close by reading one physical number off the shared matrix: the
$\vert 11\rangle$ diagonal entry $\pi^{2}+\tfrac{3g}{2}$ (verify at authoring), the first-order
shift of two bosons in one orbital, contrasted with the pair-changing entries a mean-field
treatment discards.

#### 23.5 [MSc] How do I put a bosonic field on a finite lattice (truncated) and compute its dispersion relation?

Harmonic chain of $N=24$ sites, unit masses and springs, periodic boundary conditions: the
dynamical matrix is the circulant SparseArray with $2$ on the diagonal, $-1$ on the neighbors, and
$-1$ in the corners. Take the exact route first: substitute the Fourier mode $u_j=e^{ikj}$ with
$k=2\pi m/N$ into the eigenvalue equation and Simplify to $2-2\cos k=4\sin^{2}(k/2)$, so
$\omega(k)=2\vert\sin(k/2)\vert$ exactly, doubly degenerate in $\pm m$. The refuting check reuses
the matrix, never a re-typed formula: Sqrt of the numeric Eigenvalues of the SparseArray, sorted,
against the sorted `Table` of $2\vert\sin(\pi m/N)\vert$ over $m=0,\ldots,N{-}1$, agreeing to
machine zero; an open-chain slip or sign error yields the shifted sine family and the wrong
degeneracy pattern, failing loudly. Then quantize: in normal coordinates
$\hat H=\sum_k\omega(k)\,(a_k^{\dagger}a_k+\tfrac12)$, the lattice field is $N$ independent
ladders, the multimode structure of 23.1 realized by a mechanical system. Quantify the continuum
limit with Series: $2\sin(k/2)=k-k^{3}/24+O(k^{5})$, the linear phonon branch $\omega\to\vert
k\vert$ with unit sound speed, and at the zone edge $\omega(\pi)=2$ against the linear
extrapolation $\pi$, a deficit of $1-2/\pi\approx36\%$, with the group velocity $\cos(k/2)$ (by D)
vanishing there. Close on the physics: the chain emulates a massless continuum field only at
wavelengths long against the spacing, and the zone-edge mode is a standing wave that transports
nothing.

#### 23.6 [MSc] How does a time-dependent boundary or parameter create particles from the vacuum (a dynamical-Casimir-type capstone, with the truncation made explicit)?

C5 member (`Route-Table.md`, C5 verdict): the time-dependent generator goes through NDSolveValue at
tight goals, norm drift grows with integration length, and parametric-resonance stiffness was
probed only by proxy (Landau-Zener), flagged. Single mode with modulated frequency
$\omega(t)=\omega_0(1+\epsilon\sin 2\omega_0 t)$, $\omega_0=1$, $\epsilon=1/20$: the mode function
obeys $\ddot f+\omega(t)^{2}f=0$ with vacuum data $f(0)=1/\sqrt{2\omega_0}$,
$\dot f(0)=-i\sqrt{\omega_0/2}$; integrate the equivalent first-order complex system with
NDSolveValue at PrecisionGoal and AccuracyGoal 10 or higher (WorkingPrecision above machine only
with exact rationalized input, C5 trap). The truncation is explicit and exact here: the field is
reduced to one mode, legitimate because a parametric frequency drive does not couple modes, whereas
a genuinely moving boundary would mix them and demand the multimode machinery of 23.1. The
adiabatic particle number is phase-free,
$n(t)=\bigl(\vert\dot f\vert^{2}+\omega(t)^{2}\vert f\vert^{2}\bigr)/\bigl(2\omega(t)\bigr)
-\tfrac12=\vert\beta(t)\vert^{2}$, and the refuting invariant is the Wronskian
$f\dot f^{*}-f^{*}\dot f=i$, equivalently $\vert\alpha\vert^{2}-\vert\beta\vert^{2}=1$, monitored
along the whole trajectory: drift beyond the goal tolerance refutes the run before any physics is
read. On the resonant drive the number grows exponentially, $\vert\beta\vert\sim\sinh(\lambda t)$
with parametric rate $\lambda=\epsilon\omega_0/2$ for this modulation depth (verify at authoring):
fit the late-time slope of $\log n(t)$ out to $t=200$ with LinearModelFit and compare rates in the
asymptotic window only, never pointwise at early times (the C5 wrong-benchmark trap at finite
time); because the stiffness was probed only by proxy, a goals sweep at authoring must show the
fitted rate solver-independent. The contrast is the detuned drive
$\omega(t)=\omega_0(1+\epsilon\sin(2.3\,\omega_0 t))$, outside the resonance tongue of width of
order $\epsilon\omega_0$: $n(t)$ stays bounded and oscillatory. Close on the capstone reading:
modulating a cavity parameter at $2\omega_0$ converts vacuum into pairs, the dynamical Casimir
mechanism, and it is resonance, not modulation as such, that creates particles secularly.

