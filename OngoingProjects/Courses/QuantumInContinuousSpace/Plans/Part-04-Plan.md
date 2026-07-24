# Part 4 Plan: The harmonic oscillator

8 questions. Class census, per the `Route-Table.md` class census: C1 = {4.1, 4.5}, both entries
inheriting the C1 verdict's binding traps in their prose; C0 = {4.2, 4.3, 4.4, 4.6, 4.7, 4.8}, WL
machinery named directly per the C0 row.

## Common ground

The part rests on: the oscillator Hamiltonian $\hat H=\tfrac12(\hat p^2+\hat x^2)$ with
$\hat p=-i\,d/dx$ (natural units $\hbar=m=\omega=1$ throughout; only 4.5 restores units); the
ladder pair $a=(\hat x+i\hat p)/\sqrt2$, $a^\dagger=(\hat x-i\hat p)/\sqrt2$ with
$[a,a^\dagger]=1$, the number operator $\hat n=a^\dagger a$, and $\hat H=\hat n+\tfrac12$; the
ladder action $a\vert n\rangle=\sqrt n\,\vert n-1\rangle$,
$a^\dagger\vert n\rangle=\sqrt{n+1}\,\vert n+1\rangle$; the normalized eigenfunctions
$\psi_n(x)=H_n(x)\,e^{-x^2/2}/\sqrt{2^n\,n!\,\sqrt\pi}$ with $E_n=n+\tfrac12$; and the vacuum
quadrature variances $\Delta x^2=\Delta p^2=\tfrac12$, the reference disk against which the
coherent and squeezed states of 4.6 through 4.8 are measured.

### 4.1 [BSc] How do I solve the oscillator's stationary equation and obtain the Hermite-function eigenstates with energies $E_n=n+\tfrac12$?

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

### 4.2 [BSc] How do I build the ladder operators $a=(\hat x+i\hat p)/\sqrt2$ and $a^\dagger$ as differential operators and verify $[a,a^\dagger]=1$?

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

### 4.3 [BSc] How do I generate the whole spectrum algebraically from $a|0\rangle=0$ and the number operator $\hat n=a^\dagger a$?

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

### 4.4 [BSc] How do I compute the matrix elements of $\hat x$ and $\hat p$ in the number basis, giving the truncated banded-matrix representation of the oscillator (with the truncation made explicit)?

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

### 4.5 [BSc] How do I show the ground state is a Gaussian, and restore $\omega$ and $m$ to read off the oscillator length $\sqrt{\hbar/m\omega}$?

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

### 4.6 [MSc] How do I build a coherent state $|\alpha\rangle$ as the eigenstate of $a$, and equivalently by applying $e^{\alpha a^\dagger-\alpha^* a}$ to the ground state, and show it is a minimum-uncertainty state?

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

### 4.7 [MSc] How do I build a squeezed state by applying $e^{\frac12(\xi^* a^2-\xi a^{\dagger2})}$ to the ground state and compute its unequal quadrature variances?

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

### 4.8 [MSc] How do I compute the photon-number statistics of a number state and a coherent state (the Poisson distribution)?

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
