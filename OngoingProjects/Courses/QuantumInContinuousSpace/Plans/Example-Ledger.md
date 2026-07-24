# Example Ledger: pinned primary examples, Parts 4 through 23

One primary example per question, pinned at plan time so authoring inherits a conflict-free ledger.
Rules: no state reused as a crutch across answers (PIPELINE Appendix A); a transformation question
may take a known packet as its object without that counting as reuse; an example the kernel refuses
at authoring time may be replaced within its stated class, with the swap recorded here and in the
plan. Hydrogen recurs by scope in Parts 8, 10, 12, 13, 14, 22: it is the subject of those questions,
not a crutch. The quartic double well appears in 5.4, 7.4, 19.6 deliberately, in three different
regimes (benchmark dynamics, shallow-barrier Ritz, deep-barrier instanton); the ledger treats this
as one intentional family, not a repeat.

Reserved by Parts 1 through 3 (do not reuse as primary carriers): Morse ground state; Lorentzian
$1/(x^2+a^2)$; plane wave; sinc; hydrogen $1s$ reduced radial; displaced Hermite $n=1$;
infinite-well two-level superposition; Berry-Balazs Airy beam; tent trial state; chirped sech;
boosted Hermite $n=1$; boosted Poschl-Teller excited state $\operatorname{sech}(x)\tanh(x)e^{ikx}$;
infinite-well ground state; sech-vs-Gaussian contrast pair.

Class codes cite `Route-Table.md`. C0 = no differential equation (operator algebra, quadrature,
transforms, truncated linear algebra).

## Part 4. The harmonic oscillator

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 4.1 | BSc | stationary ODE solved; the $n=3$ state (three nodes) read off the general solution, quantization from decay at $\pm\infty$ | vibrational spectra | quantization earned from normalizability, not asserted | C1 |
| 4.2 | BSc | $a,a^\dagger$ act on generic $f(x)$ for $[a,a^\dagger]=1$; concrete action on the $n=2$ state showing $\sqrt n$, $\sqrt{n+1}$ | ladder algebra | generic-function proof plus a nodal concrete carrier | C0 |
| 4.3 | BSc | $a\psi_0=0$ solved as a first-order ODE, ladder up to $n=4$; residual of $\hat H\psi_4=\tfrac92\psi_4$ | algebraic spectrum | the whole ladder built, not one level | C0 |
| 4.4 | BSc | truncated $N=8$ matrices of $\hat x,\hat p$; banded structure; commutator corner defect $[\hat x,\hat p]=i\,(\mathbb 1-(N{+}1)\,\vert N\rangle\langle N\vert)$ | truncation diagnostics | the corner defect is the honest truncation exhibit | C0 |
| 4.5 | BSc | named-trivial: the Gaussian ground state with $\hbar,m,\omega$ restored; contrast: CO stretch vs a Rb atom in an optical tweezer, oscillator lengths orders of magnitude apart | molecular and cold-atom scales | the contrast is two real laboratory scales | C1 |
| 4.6 | MSc | coherent state with complex $\alpha=2e^{i\pi/3}$; the two constructions (eigenstate of $a$; displaced vacuum) shown equal | quantum optics | complex $\alpha$, not real; two constructions must agree | C0 |
| 4.7 | MSc | squeezed vacuum with $\xi=r e^{i\pi/3}$: rotated variance ellipse, unequal quadratures | quantum optics | the phase of $\xi$ kept, not $\theta=0$ | C0 |
| 4.8 | MSc | number state $n=4$ vs coherent $\vert\alpha\vert^2=4$: same mean energy, deterministic vs Poisson statistics, Mandel $Q$ | photon counting | a same-mean discriminating pair | C0 |

## Part 5. The time-dependent Schrodinger equation as a PDE

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 5.1 | BSc | Poschl-Teller well $V=-3\operatorname{sech}^2x$ (exactly two levels $E_0=-2$, $E_1=-\tfrac12$): phase-only stationary evolution vs the two-level beat | exactly solvable well | non-box system; beat frequency $\Delta E$ is closed form | C0 |
| 5.2 | BSc | cusp packet $e^{-\vert x-x_0\vert}$ released off-center in a harmonic trap | trapped-atom quench | the cusp radiates high momenta immediately; a genuine grid stress | C3 |
| 5.3 | BSc | named-trivial: free Gaussian with drift $k_0$; contrast: a same-width sech packet does not spread self-similarly | matter-wave optics | the contrast exposes what Gaussian self-similarity hides | C3 |
| 5.4 | MSc | split-step Fourier vs NDSolve on a packet oscillating in a symmetric quartic double well | tunneling dynamics | two independent integrators, measurable disagreement | C3 |
| 5.5 | MSc | off-center triangular packet in the infinite well; exact revival at $T_{\mathrm{rev}}=4L^2/\pi$ | box revivals | non-smooth packet, many modes; revival exact for the quadratic spectrum | C3 |
| 5.6 | MSc | coherent state $\alpha=2$ vs a squeezed-displaced state: identical $\langle x\rangle(t)$, constant vs breathing width | trapped-ion optics | isolates what "coherent" adds to "displaced" | C3 |
| 5.7 | MSc | compact-support raised-cosine packet, boost $k_0$, on an Eckart barrier $V_0\operatorname{sech}^2(x/a)$ with exact $T(E)$ | reaction barriers | zero initial barrier overlap; exact transmission anchor | C3 |
| 5.8 | MSc | free kernel and Mehler kernel; free kernel applied to the Airy beam reproduces its accelerating propagation (transformation-object reuse, allowed) | propagator theory | kernel earns the $\delta(x-x')$ limit, TDSE residual, composition law | C0 |
| 5.9 | BSc | packet in the quartic well $V=x^4/4$: $\tfrac{d\langle p\rangle}{dt}=-\langle V'(x)\rangle\neq-V'(\langle x\rangle)$ exhibited | Ehrenfest fine print | anharmonicity makes the theorem's actual content visible | C3 |
| 5.10 | BSc | boosted supergaussian $e^{-x^4}$ packet crossing the reflectionless well $V=-\operatorname{sech}^2x$: norm and current conserved while transmission is total | integrable scattering | the fidelity check rides a striking exact property | C3 |

## Part 6. Scattering in one dimension

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 6.1 | BSc | potential step: both $E>V_0$ (partial reflection) and $E<V_0$ (evanescent, $R=1$ with penetration); current-verified $R+T=1$ | heterojunction | both regimes; current-based, not amplitude-based | C4 |
| 6.2 | BSc | rectangular barrier: closed-form $T(E)$; opaque limit $\sim e^{-2\kappa a}$ | tunneling junctions | class representative; symbolic form plus limit | C4 |
| 6.3 | MSc | rectangular well transmission: $T=1$ at $k'a=n\pi$, parameters placed so a resonance sits at low energy | Ramsauer-Townsend | resonance positions predicted, then found | C4 |
| 6.4 | MSc | double barrier (resonant-tunneling diode): transfer-matrix product; inter-barrier resonance sharper than either barrier alone | mesoscopic transport | assembly from interface matrices does real work | C4 |
| 6.5 | MSc | finite-well S-matrix: unitarity symbolically; poles at imaginary $k$ reproduce the even/odd bound-state conditions | scattering theory | the pole-bound-state identity computed, not asserted | C4 |
| 6.6 | MSc | attractive delta well: exact phase shift, Levinson $\delta(0)-\delta(\infty)=\pi$ with exactly one bound state | contact interactions | exact Levinson count; delta handled by Integrate, never NIntegrate | C4 |

## Part 7. Approximation methods

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 7.1 | BSc | quartic perturbation $\lambda x^4$ on the oscillator: $E_n^{(1)}$, $E_n^{(2)}$ for general $n$, compared to numerics | anharmonic vibration | the whole family, not one level; second order hints at divergence | C0 |
| 7.2 | MSc | particle on a ring: degenerate $m=\pm1$ doublet split by $\lambda\cos2\phi$ | quantum rotor | the naive nondegenerate formula fails; the secular matrix is forced | C0 |
| 7.3 | BSc | named-trivial: Gaussian trial on $V=x^4/4$; contrast: the two-parameter trial $(1+cx^2)e^{-bx^2/2}$ strictly lowers the bound | variational method | the contrast quantifies the Gaussian's deficit | C0 |
| 7.4 | MSc | shallow-barrier quartic double well in a nonorthogonal basis of two displaced Gaussians plus polynomial corrections: genuinely generalized eigenproblem, splitting estimate | molecular inversion | the overlap matrix is the point of the question | C0 |
| 7.5 | BSc | WKB quantization of the quartic well vs its numeric spectrum; Morse anchor where WKB is exact | semiclassics | exactness anchor plus an honest error-vs-$n$ curve | C9 |
| 7.6 | MSc | inverted parabolic barrier: WKB Gamow exponent vs the exact parabolic-cylinder transmission | barrier penetration | an exact benchmark exists for this barrier | C9 |
| 7.7 | MSc | photodetachment of the delta-well bound state by an oscillating force: golden-rule rate with continuum DOS; numeric amplitude-ODE cross-check | model photoionization | discrete-to-continuum with closed-form matrix elements | C5 |
| 7.8 | MSc | infinite well suddenly vs slowly doubled, $L\to2L$: sudden overlaps $P(n)$ vs adiabatic tracking; a rate sweep spans the crossover | quench dynamics | both approximations on one system, crossover visible | C5 |
| 7.9 | MSc | Airy connection derived at a linear turning point $V=Fx$; verified against the quantum-bouncer spectrum (Airy zeros) | semiclassical matching | the $\pi/4$ phase earned, then checked against exact zeros | C9 |
| 7.10 | MSc | Temple lower bound plus the 7.3 upper bound sandwiching the quartic ground energy around the numeric value | rigorous bounds | a two-sided cage, not a single estimate | C0 |

## Part 8. General theorems and structural methods

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 8.1 | BSc | virial on hydrogen $2p$ ($2\langle T\rangle=-\langle V\rangle$) and on the quartic ground state ($2\langle T\rangle=4\langle V\rangle$) | power-law potentials | two exponents discriminate the theorem's content | C0 |
| 8.2 | MSc | Hellmann-Feynman on hydrogen with $l$ as the parameter: $\langle1/r^2\rangle$ from $\partial E/\partial l$; Integrate cross-check | atomic expectation values | a nonobvious expectation value without an integral | C0 |
| 8.3 | MSc | TRK sum rule in the infinite well: partial sums of $f_{1n}$ converge to 1; the oscillator saturates in one term | sum rules | convergence watched, not asserted | C0 |
| 8.4 | MSc | oscillator $n=30$: locally averaged $\vert\psi\vert^2$ vs the classical arcsine density $1/(\pi\sqrt{A^2-x^2})$ | correspondence principle | quantitative local-average comparison, not a picture | C9 |
| 8.5 | MSc | infinite well factorized: superpartner $V_2\propto\csc^2(\pi x/L)$; isospectrality minus the ground state verified on the exactly solvable partner | SUSY QM | the partner is a different-looking exactly solvable well | C0 |
| 8.6 | MSc | shape invariance: Poschl-Teller $\lambda$ ladder and the oscillator spectrum from the recursion alone | algebraic spectra | spectra with no ODE solving at all | C0 |

## Part 9. Periodic potentials and band structure

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 9.1 | BSc | cosine lattice: Bloch form $e^{ikx}u_k(x)$ verified on Mathieu band states, band energy from the Mathieu characteristic | optical lattices | exactly solvable lattice; Bloch form verified, not asserted | C1 |
| 9.2 | MSc | Kronig-Penney rectangular lattice: band condition, first two bands and the gap | solid-state textbook | the transcendental condition solved honestly | C4 |
| 9.3 | MSc | attractive Dirac comb: band edges; collapse limit $b\to0$, $V_0b$ fixed, recovers 9.2 | point-scatterer lattice | the limit connects the two lattices; delta via Integrate only | C4 |
| 9.4 | MSc | density of states of the lowest cosine-lattice band: $g(E)=1/\vert dE/dk\vert$, van Hove edges | band thermodynamics | the $1/\sqrt{}$ singularity is the content; checked against state counting | C0 |
| 9.5 | MSc | deep cosine lattice: $E(k)\approx\varepsilon-2J\cos ka$; Wannier function exponential localization vs depth; $J$ extracted two independent ways | cold atoms in lattices | band fit and hopping integral must agree | C0 |
| 9.6 | MSc | tilted deep lattice: $\langle x\rangle(t)$ oscillates at $\omega_B=Fa$; Wannier-Stark ladder | Bloch oscillations | oscillation instead of transport is the counterintuitive point | C3 |

## Part 10. Two and three dimensions: separation of variables

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 10.1 | BSc | rectangle with $L_x/L_y=\sqrt2$: product eigenstates; a 2D numeric eigensolve reproduces $E_{n_xn_y}$ | quantum dots | incommensurate ratio isolates separability from degeneracy accidents | C7 |
| 10.2 | BSc | square box: degeneracy of $E\propto50=1^2+7^2=5^2+5^2$ (three states), counted number-theoretically | spectral counting | degeneracy from representations as sums of two squares | C0 |
| 10.3 | BSc | 2D isotropic oscillator $N=2$ shell: Cartesian triplet vs polar $(n_r,m)$ states with the explicit unitary between them | 2D traps | same shell, two bases, explicit map | C7 |
| 10.4 | BSc | hydrogen vs deuterium: reduced mass in the Rydberg, the measured Balmer isotope shift | spectroscopy | the two-body reduction cashes into a measured number | C0 |
| 10.5 | BSc | generic central $V(r)$: full separation; $Y_{lm}$ substitution leaves the radial equation with $l(l+1)$ | formalism | derived for generic $V$, not one example | C7 |
| 10.6 | BSc | hydrogen effective radial potentials for $l=0,1,2$: barrier heights and turning points | atomic structure | the family, not a single $l$ | C0 |
| 10.7 | MSc | hydrogen in parabolic coordinates: separation constants, $n=n_1+n_2+\vert m\vert+1$ | Stark groundwork | re-counts the degeneracy in the frame 13.5 will use | C7 |

## Part 11. Orbital angular momentum in continuous space

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 11.1 | BSc | $L_x,L_y,L_z$ derived from $\vec r\times\vec p$ by change of variables, acting on generic $f(\theta,\phi)$ | formalism | full derivation, not quoted forms | C0 |
| 11.2 | BSc | commutators and Casimir on generic $f$; eigenvalue check on the full $l=2$ multiplet | formalism | generic proof plus a concrete family | C0 |
| 11.3 | BSc | $Y_{2m}$ built from $L_+Y_{22}=0$ (a first-order PDE) and laddered down; orthonormality integrals; compared to the built-in | representation construction | constructed, then compared, not quoted | C0 |
| 11.4 | MSc | $L_\pm$ on the full $l=3$ multiplet: $\sqrt{l(l+1)-m(m\pm1)}$ for all $m$, edge annihilation at $m=\pm3$ | ladder algebra | the whole multiplet including the edges | C0 |
| 11.5 | MSc | rotate $Y_{2m}$ by Euler angles with $\beta$ symbolic: expansion coefficients equal the Wigner $D^2$ entries | rotation theory | symbolic $\beta$ where the kernel allows | C0 |
| 11.6 | MSc | $1\otimes1\to0,1,2$: ClebschGordan vs the Gaunt triple-$Y$ integrals, two independent computations | angular-momentum coupling | the same coefficients from two machineries | C0 |

## Part 12. Central potentials and the hydrogen atom

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 12.1 | BSc | Coulomb radial ODE solved; quantization from normalizability; $E_n=-1/(2n^2)$ for $n=1..4$ | hydrogen | class representative; spectrum earned from the ODE | C1 |
| 12.2 | BSc | $R_{31}$ (one node, $l=1$) from the associated Laguerre form; general $(n,l)$ normalization | hydrogen radial functions | nodal case plus the general normalization | C1 |
| 12.3 | BSc | full $\psi_{321}$ (a $3d$, $m=1$ state): density, radial distribution, $\langle r\rangle$ vs the closed formula | orbital structure | a $d$ orbital, not $1s$ | C0 |
| 12.4 | BSc | spherical finite well: $l=0$ and $l=1$ bound states numerically, checked against Bessel matching conditions | nuclear-well toy | the $l>0$ case exceeds the textbook $l=0$ | C2 |
| 12.5 | MSc | 3D isotropic oscillator $N=2$ shell: Cartesian count 6 equals spherical $5_{(l=2)}+1_{(l=0)}$ with the explicit map | shell structure | degeneracy reconciled, not just counted | C7 |
| 12.6 | MSc | Runge-Lenz: $[\hat H,\hat{\vec A}]=0$ with operator ordering done honestly; $\hat{\vec A}$ rotates $2s$ into $2p$ | dynamical symmetry | ordering is the hard part; degeneracy explained | C0 |
| 12.7 | MSc | hydrogen levels from a numeric radial eigensolve on $(0,R)$: negative-spectrum shift trap live, truncation artifact vs $R$ | numerical atomic structure | class representative; the trap is teaching content | C2 |
| 12.8 | MSc | sodium-like screened Coulomb: quantum defects $\delta_l$ from numeric levels, Rydberg-formula fit | alkali spectroscopy | numerics land on the empirical defect | C2 |
| 12.9 | MSc | radial WKB with and without the Langer correction: only the corrected form reproduces $-1/(2n^2)$ | semiclassical atoms | the on/off contrast makes the correction load-bearing | C9 |

## Part 13. Electromagnetic coupling

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 13.1 | MSc | uniform electric field in two gauges (scalar $-Ex$ vs time-dependent $A(t)$): gauge function, phase transformation, identical $\vert\psi\vert^2$ | gauge invariance | the same physics computed in both gauges | C0 |
| 13.2 | MSc | Landau gauge reduction to a shifted oscillator; symmetric-gauge $m$ degeneracy; degeneracy per area $B/2\pi$ | quantum Hall groundwork | two gauges, one spectrum | C1 |
| 13.3 | MSc | charged particle on a flux-threaded ring: $E_m\propto(m-\Phi/\Phi_0)^2$, periodicity in $\Phi_0$, persistent current | Aharonov-Bohm rings | the phase becomes a measurable spectral shift | C0 |
| 13.4 | MSc | hydrogen $n=2$: orbital Zeeman triplet from $-\tfrac12BL_z$ | Zeeman spectroscopy | the degenerate manifold handled properly | C0 |
| 13.5 | MSc | linear Stark on the $n=2$ manifold: the $4\times4$ secular problem; parabolic states diagonalize; splitting $\pm3E$ | Stark spectroscopy | linearity in the field is the surprise; uses 10.7's frame | C0 |
| 13.6 | MSc | harmonic trap center dragged around a closed loop in $(x_0,p_0)$: geometric phase equals the enclosed phase-space area; adiabatic evolution cross-check | geometric phase | a continuous-variable Berry phase, not the spin cliche | C0 |

## Part 14. Spin coupled to spatial motion

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 14.1 | MSc | position-spin entangled spinor $\operatorname{sech}(x-a)\vert\uparrow\rangle+\operatorname{sech}(x+a)\tanh(x+a)\vert\downarrow\rangle$: reduced spin purity below 1 | spinor wave mechanics | different spatial structure per branch; entanglement quantified | C0 |
| 14.2 | MSc | hydrogen $2p$ fine structure: $\xi(r)\vec L\cdot\vec S$ with $\langle1/r^3\rangle$ closed form; $j=3/2,1/2$ splitting | fine structure | full radial-angular assembly, a real number out | C0 |
| 14.3 | MSc | $l=1\otimes s=\tfrac12\to j=\tfrac32,\tfrac12$ via CG; $J^2,J_z$ verified diagonal | recoupling | the basis change later parts reuse | C0 |
| 14.4 | MSc | Stern-Gerlach: two-component sech packet in a field gradient; spatial split entangles with spin; overlap decay as measurement | Stern-Gerlach | measurement as dynamics, not postulate | C3 |
| 14.5 | MSc | Pauli equation in uniform $B$: Landau levels plus Zeeman with $g=2$; the $(n,\uparrow)$ and $(n{+}1,\downarrow)$ degeneracy coincidence | Pauli theory | the accidental degeneracy is the exhibit | C1 |
| 14.6 | MSc | hydrogen $2p$ with spin-orbit and $B$: $6\times6$ diagonalization; weak and strong limits; Paschen-Back crossover | anomalous Zeeman | the full crossover with both limits recovered | C0 |

## Part 15. Identical particles in continuous space

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 15.1 | MSc | two particles in the infinite well, orbitals $n=1,2$: symmetrized vs antisymmetrized; fermion node line $x_1=x_2$ | two-body statistics | different orbitals so exchange bites | C0 |
| 15.2 | MSc | oscillator orbitals $n=0,1$: pair correlation $g(x_1,x_2)$; fermion hole vs boson bunching at coincidence | quantum statistics | same orbitals, opposite statistics, quantitative | C0 |
| 15.3 | MSc | $N=3$ Slater determinant in the well ($n=1,2,3$): total density and nodal structure | few-fermion systems | genuine $N=3$, not the $N=2$ shortcut | C0 |
| 15.4 | MSc | helium-like $1s2s$: direct $J$ and exchange $K$ Coulomb integrals in closed form; para/ortho splitting $2K$ | helium spectroscopy | both integrals computed, the splitting earned | C0 |
| 15.5 | MSc | softened-Coulomb one-dimensional helium: Hartree-Fock SCF on a grid vs exact two-electron grid diagonalization | computational atomic physics | self-consistency graded against an exact benchmark | C6 |

## Part 16. Density operators, mixed states, and the Wigner function

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 16.1 | MSc | equal-population mixture vs coherent superposition of the two Poschl-Teller levels: same populations, purity and off-diagonal kernel discriminate | decoherence groundwork | the discriminating pair, not one state | C0 |
| 16.2 | MSc | two coupled oscillators' ground state: reduced single-mode state, purity vs coupling, entanglement entropy | bipartite Gaussian states | partial trace with a tunable knob, closed forms | C0 |
| 16.3 | MSc | Wigner function of the $n=1$ Fock state: negativity $W(0,0)=-1/\pi$ | phase-space QM | negativity is the non-classical smoking gun | C0 |
| 16.4 | MSc | even cat state $\alpha=2$ beside number and coherent states: blobs plus interference fringes | cat states | the fringes are the physics | C0 |
| 16.5 | MSc | the same trio under Husimi $Q$ (all positive) and Glauber $P$ (delta, Gaussian, singular) | phase-space hierarchy | the singularity ladder across representations | C0 |
| 16.6 | MSc | displaced thermal state: Wigner Gaussian with $\coth(\beta/2)$ width; Moyal evolution is exact rigid rotation for quadratic $H$ | thermal quantum optics | Moyal reduces to Liouville exactly here, and is verified | C0 |

## Part 17. Continuous-variable quantum optics and information

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 17.1 | MSc | single-mode quadratures: commutator, vacuum disk, shot-noise variance | quantum optics | operational definitions with variances | C0 |
| 17.2 | MSc | $D^\dagger XD$, $S^\dagger XS$ by BCH; the Wigner ellipse transformation alongside | Gaussian optics | operator and phase-space pictures must agree | C0 |
| 17.3 | MSc | Hong-Ou-Mandel on $\vert1,1\rangle$: coincidence null at 50:50; dip visibility vs mode-overlap parameter | photonic interference | the distinguishability sweep, not just the null | C0 |
| 17.4 | MSc | two-mode squeezed vacuum: Schmidt form, thermal reduced state, $\Delta(X_1-X_2)^2=e^{-2r}$ | EPR correlations | entanglement quantified two ways | C0 |
| 17.5 | MSc | Duan criterion: TMSV violates for all $r>0$; a two-mode coherent product sits exactly on the boundary | CV entanglement detection | a passing and a boundary state, both exact | C0 |
| 17.6 | MSc | Jaynes-Cummings from $\vert e,0\rangle$: vacuum Rabi at $2g$, detuning sweep, truncation checked | cavity QED | truncation explicit; detuning generalizes the resonance | C5 |
| 17.7 | MSc | homodyne on the squeezed state: quadrature histograms vs LO phase; heterodyne as $Q$ sampling | quantum measurement | operational reconstruction, not formal definition | C0 |

## Part 18. Open quantum systems in continuous space

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 18.1 | MSc | damped coherent state: stays coherent with $\alpha(t)=\alpha e^{-\gamma t/2}$; energy decay | cavity damping | class representative; pointer-like stability is the physics | C5 |
| 18.2 | MSc | Caldeira-Leggett packet on a truncated basis: damping plus diffusion; the short-time positivity violation exhibited honestly | quantum Brownian motion | the known artifact is part of the lesson | C5 |
| 18.3 | MSc | superposition of two coherent states under damping: decoherence rate scaling as $\vert\alpha_1-\alpha_2\vert^2$ | einselection | the scaling law is pointer-basis selection | C5 |
| 18.4 | MSc | damped-oscillator Wigner function: drift-diffusion form; Gaussian moment ODEs close exactly; relaxation to thermal | phase-space open systems | the PDE cross-checked by exact moments | C5 |
| 18.5 | MSc | quantum-jump unravelling of a damped $n=2$ Fock state: trajectories vs ensemble average; waiting-time statistics | quantum trajectories | the stochastic average earns the deterministic equation | C5 |

## Part 19. Path integrals

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 19.1 | MSc | free propagator from $N$ time slices: Gaussian recursion in closed form, then $N\to\infty$ | path integrals | the limit taken, not asserted | C0 |
| 19.2 | MSc | Mehler kernel from the same slicing on the oscillator (determinant recursion) | path integrals | the determinant recursion in closed form | C0 |
| 19.3 | MSc | oscillator partition function: discretized imaginary-time trace vs $1/(2\sinh(\beta/2))$ vs the spectrum sum | thermal QM | three routes agree | C0 |
| 19.4 | MSc | Van Vleck stationary-phase propagator: exact for the oscillator; the first anharmonic correction located | semiclassics | exactness for quadratic actions placed precisely | C9 |
| 19.5 | MSc | quartic-well $E_0$ from numeric imaginary-time transfer matrices: large-$\beta$ decay rate | lattice path integrals | ground energy from pure matrix multiplication | C0 |
| 19.6 | MSc | deep-barrier quartic double well: instanton action in $e^{-S_0}$ vs the numeric splitting | instantons | the exponent predicted, then measured | C9 |

## Part 20. Three-dimensional scattering theory

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 20.1 | MSc | hard sphere: $\delta_l$ exact from Bessel matching; a finite well numerically beside it | partial waves | exact anchor plus numeric generalization | C4 |
| 20.2 | MSc | hard-sphere cross section from the $\delta_l$ sum: forward shadow, total $\sigma\to2\pi a^2$ | cross sections | the shadow-doubling surprise | C0 |
| 20.3 | MSc | Yukawa Born amplitude in closed form; the screening-off limit recovers Rutherford | Born approximation | the limit connects to 20.6 | C0 |
| 20.4 | MSc | finite-well scattering length vs depth: divergence at each new bound state; optical-theorem check | ultracold scattering | the resonance structure of $a$, not one number | C4 |
| 20.5 | MSc | well-plus-barrier shape resonance: Breit-Wigner fit of $\delta_0(E)$; Levinson recount with the resonance | resonance physics | fit parameters carry meaning | C4 |
| 20.6 | MSc | Coulomb scattering with the full Coulomb wave functions: Rutherford $d\sigma/d\Omega$ exactly | Rutherford scattering | the exceptional long-range case done exactly | C1 |

## Part 21. Nonlinear and mean-field wave mechanics

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 21.1 | MSc | trapped repulsive condensate: GPE ground state by imaginary-time evolution; healing length, chemical potential | BEC | class representative | C6 |
| 21.2 | MSc | bright soliton (attractive) and dark tanh notch (repulsive): closed forms residual-verified; moving-soliton stability under evolution | nonlinear matter waves | both species, residual plus dynamics | C6 |
| 21.3 | MSc | Thomas-Fermi profile vs the full GPE numeric at large $g$: inverted parabola with edge healing | strongly interacting BEC | the approximation's edge failure shown | C6 |
| 21.4 | MSc | two softened-contact bosons: Hartree self-consistent orbital vs the exact two-body grid solution | mean-field validity | mean field graded against exact | C6 |

## Part 22. Relativistic wave equations

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 22.1 | MSc | Klein-Gordon: $\pm E$ branches; a packet with negative-energy admixture makes the density indefinite | relativistic waves | the pathology computed, not narrated | C0 |
| 22.2 | MSc | gamma algebra verified; free spinors $u,v$ with $\bar uu=2m$; helicity | Dirac formalism | full spinor construction with checks | C0 |
| 22.3 | MSc | nonrelativistic reduction of Dirac in a field: Pauli equation with $g=2$ emerging | Dirac to Pauli | $g=2$ derived, not quoted | C0 |
| 22.4 | MSc | Dirac hydrogen: the coupled radial system, exact $E_{nj}$; compared to 14.2's perturbative fine structure | relativistic hydrogen | class representative; exact vs perturbative | C8 |
| 22.5 | MSc | supercritical step $V_0>2m$: Klein transmission; zitterbewegung of a packet with $\pm E$ interference | Klein paradox | both classic pathologies in one framework | C4 |

## Part 23. From one particle to fields: the second-quantization bridge

| q | tag | pinned example | field | nontrivial because / contrast | class |
|---|---|---|---|---|---|
| 23.1 | MSc | three-mode Fock space: occupation arithmetic, per-mode ladders commuting across modes | second quantization | multimode from the start | C0 |
| 23.2 | MSc | field operator on the first $N$ well modes: the truncated $[\hat\psi(x),\hat\psi^\dagger(x')]$ kernel converging to $\delta(x-x')$ | field construction | the delta emerges under truncation refinement | C0 |
| 23.3 | MSc | scalar field in a box: mode-sum vacuum energy, cutoff scaling, a Casimir-flavored regularized difference | vacuum energy | the divergence structure made explicit | C0 |
| 23.4 | MSc | two-boson sector: second-quantized $\hat H$ matrix elements equal the symmetrized first-quantized ones, on an interacting term | equivalence theorem | checked on interaction, not the free part | C0 |
| 23.5 | MSc | harmonic chain of $N$ sites: dispersion $\omega(k)=2\vert\sin(k/2)\vert$, continuum limit, zone-edge deviation | lattice field theory | the continuum limit quantified | C0 |
| 23.6 | MSc | parametrically driven mode frequency $\omega(t)$: Bogoliubov ODEs, particles from vacuum, resonance at $2\omega$ | dynamical Casimir | pair creation computed from mode equations | C5 |
