# Part 12 Plan: Central potentials and the hydrogen atom

Questions: 9 (12.1 through 12.9). Class census, per the class census table in `Route-Table.md`:
C1: 12.1, 12.2; C0: 12.3, 12.6; C2: 12.4, 12.7, 12.8; C7: 12.5; C9: 12.9. Both probed class
representatives live in this part (12.1 for C1, 12.7 for C2), so their entries cite gates verbatim;
the Route-Table C0 kernel facts bind throughout (`E` is Euler's number, use `en` for energies;
`Simplify` closes residuals and norms only under explicit integer, positivity, and reality
assumptions).

## Common ground

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

## Per-question entries

### 12.1 [BSc] How do I solve the Coulomb radial equation and obtain the bound-state energies $E_n=-1/(2n^2)$?

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

### 12.2 [BSc] How do I build the radial wavefunctions from the associated Laguerre polynomials and normalize them?

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

### 12.3 [BSc] How do I assemble the full $\psi_{nlm}$, plot probability densities and the radial distribution, and compute $\langle r\rangle$?

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

### 12.4 [BSc] How do I find the bound states of a spherical square well numerically?

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

### 12.5 [MSc] How do I solve the three-dimensional isotropic oscillator and reconcile its Cartesian and spherical degeneracy counts?

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

### 12.6 [MSc] How do I expose the extra Coulomb degeneracy through the conserved Runge-Lenz vector (the dynamical symmetry)?

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

### 12.7 [MSc] How do I compute central-potential eigenstates by applying `NDEigensystem` to the radial equation?

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

### 12.8 [MSc] How do I compute the quantum defect of an alkali-like screened-Coulomb potential?

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

### 12.9 [MSc] How do I apply WKB to the radial equation with the Langer correction $l(l+1)\to(l+\tfrac12)^2$ and recover the Coulomb and oscillator spectra?

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
