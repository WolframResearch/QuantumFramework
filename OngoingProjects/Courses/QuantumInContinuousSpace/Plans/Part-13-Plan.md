# Part 13 Plan: Electromagnetic coupling

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

### 13.1 [MSc] How do I implement minimal coupling $\hat p\to\hat p-A$ and show how $\psi$ transforms under a gauge change?

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

### 13.2 [MSc] How do I derive the Landau levels of a charged particle in a uniform magnetic field?

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

### 13.3 [MSc] How do I compute the Aharonov-Bohm phase acquired around a flux line?

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

### 13.4 [MSc] How do I compute the normal (orbital) Zeeman splitting of hydrogen perturbatively from the $-\tfrac12 B L_z$ coupling?

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

### 13.5 [MSc] How do I compute the linear Stark effect on the degenerate $n=2$ hydrogen manifold?

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

### 13.6 [MSc] How do I compute Berry's geometric phase for a state carried adiabatically around a closed loop in parameter space?

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
