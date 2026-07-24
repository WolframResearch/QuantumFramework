# Part 10 Plan: Two and three dimensions: separation of variables

Questions: 7 (10.1 through 10.7). Class census: C7 for 10.1, 10.3, 10.5, 10.7 and C0 for 10.2,
10.4, 10.6, per the `Route-Table.md` class census (C7 and C0 rows) and the full C7 verdict block.

## Common ground

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

## Per-question entries

### 10.1 [BSc] How do I separate variables in a two- or three-dimensional box and form the product eigenfunctions?

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

### 10.2 [BSc] How do I find and explain the degeneracies of a square or cubic box?

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

### 10.3 [BSc] How do I solve the two-dimensional isotropic oscillator and count its degeneracies?

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

### 10.4 [BSc] How do I separate the two-body problem into center-of-mass and relative motion, reducing it to a one-body problem with the reduced mass $\mu$?

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

### 10.5 [BSc] How do I separate the Schrodinger equation in spherical coordinates into a radial and an angular equation?

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

### 10.6 [BSc] How do I build the effective radial potential and read off the centrifugal barrier?

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

### 10.7 [MSc] How do I separate variables in parabolic coordinates (the natural frame for the Stark and Coulomb problems)?

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
