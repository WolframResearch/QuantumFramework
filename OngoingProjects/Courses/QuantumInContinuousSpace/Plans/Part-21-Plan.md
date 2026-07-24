# Part 21 Plan: Nonlinear and mean-field wave mechanics

Questions: 4 (21.1 through 21.4). Class census: C6 for all four questions, per the `Route-Table.md`
class census row C6 and the full C6 verdict block; 21.1 is the probed class representative, so
entries cite its gates verbatim.

## Common ground

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

## Per-question entries

### 21.1 [MSc] How do I find the ground state of the Gross-Pitaevskii equation by imaginary-time evolution?

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

### 21.2 [MSc] How do I build the bright and dark solitons of the nonlinear Schrodinger equation?

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

### 21.3 [MSc] How do I apply the Thomas-Fermi approximation to a trapped condensate?

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

### 21.4 [MSc] How do I solve the Hartree mean-field equation as a self-consistent nonlinear Schrodinger problem?

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
