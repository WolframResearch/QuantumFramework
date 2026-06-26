# Quantum in Continuous Space: A Comprehensive Question List (BSc through MSc, done in the Wolfram Language)

A curriculum-as-questions for learning to *do* wave mechanics by computing it: the quantum
theory of a particle (and a few particles, and a field as its limit) living in a continuous,
infinite-dimensional Hilbert space $L^2$. The spine is the Schrodinger equation read two ways:
as an ordinary differential equation whose eigenvalues are the energy spectrum (the
time-independent problem), and as a partial differential equation that propagates a state in
time (the time-dependent problem). Each question is a thing you build, solve, or verify in the
Wolfram Language. Questions are tagged **[BSc]** (foundational, expected in a strong
undergraduate course) or **[MSc]** (advanced, graduate).

This is the continuous-space companion to the finite-dimensional list in
`../QuantumInDiscreteSpace/Question-List.md`. The two are complementary by construction: that
list deliberately excluded everything infinite-dimensional (wavefunctions, the harmonic
oscillator, coherent states, the hydrogen atom, scattering, fields); this list is exactly that
excluded material.

**Conventions.** Natural units $\hbar=m=1$ throughout (and the oscillator frequency $\omega=1$
where natural), so the formulas stay clean and the symbolic and numerical solvers are
unobstructed. A handful of items deliberately restore $\hbar$, $m$, or $\omega$ to exhibit a
physical scale (the oscillator length, the Bohr radius); these are noted where they occur.
Every driven or coupled system uses classical c-number electromagnetic fields until the final
part, where the field itself is quantized.

**Tooling.** Strictly pure Wolfram Language: `DSolve` and `DSolveValue` for closed-form ODE
and PDE solutions, `NDEigensystem` and `NDSolve` for numerical spectra and propagation,
`FourierTransform` and `Fourier` for the position-momentum bridge, `Integrate`/`NIntegrate`
for normalization and matrix elements, `FindRoot` for transcendental quantization conditions,
and hand-built schemes (Numerov, split-step Fourier, transfer matrix, imaginary-time descent)
where the physics is the algorithm. No external quantum package is used; where a continuous
object must become a finite array to be computed (a truncated number basis, a spatial grid),
the truncation is made explicit and its convergence is checked, not hidden.

**How verification works here, and why it differs from the finite-dimensional companion.** In
finite dimensions a claim can be settled by *exact* equality of two finite arrays. In
continuous space there is in general no exact finite matrix, the spectrum can be continuous,
and an eigenstate can be non-normalizable (delta-normalized). So a result is confirmed in one
of three ways, in order of preference: (1) a closed form when the special function exists
(Hermite, Laguerre, Airy, spherical harmonics, parabolic-cylinder, Coulomb wave functions),
checked by substitution back into the differential equation; (2) agreement between two
independent methods on the same quantity (an analytic spectrum against `NDEigensystem`, a
`DSolve` propagation against split-step Fourier); (3) numerical convergence as a grid is
refined or a basis is enlarged. The representation and discretization machinery is therefore
itself early content, not a footnote.

The order is a valid prerequisite chain: a question never needs a concept or skill that a
later question introduces. You can work straight down. This is why the wavefunction, the
position and momentum operators, and the Fourier bridge come first; the harmonic oscillator
(which underlies coherent states, the Wigner function, quantum optics, and the field bridge)
comes early; angular momentum and central potentials precede three-dimensional scattering;
and spin coupled to spatial motion and electromagnetic coupling precede the relativistic wave
equations.

Each question is one concept. Where a single computation has several inseparable facets (an
inner product and its modulus and transition probability) those stay together; genuinely
distinct concepts are separate questions.

---

## Part 1. The wavefunction and the state space $L^2$
- 1.1 [BSc] How do I represent a one-dimensional wavefunction $\psi(x)$ as a Wolfram Language function and plot its probability density $|\psi(x)|^2$?
- 1.2 [BSc] How do I normalize a wavefunction with $\int_{-\infty}^{\infty}|\psi|^2\,dx=1$, and recognize a state that cannot be normalized?
- 1.3 [BSc] How do I compute the Born-rule probability $\int_a^b|\psi|^2\,dx$ of finding the particle in a region?
- 1.4 [BSc] How do I compute the inner product $\langle\phi|\psi\rangle=\int\phi^*\psi\,dx$, its modulus, and the transition probability $|\langle\phi|\psi\rangle|^2$?
- 1.5 [BSc] How do I compute the position expectation $\langle x\rangle$ and the spread $\Delta x=\sqrt{\langle x^2\rangle-\langle x\rangle^2}$?
- 1.6 [BSc] How do I compute the probability current density $j(x,t)$ and verify the continuity equation $\partial_t|\psi|^2+\partial_x j=0$?
- 1.7 [BSc] How do I represent a state on a discrete spatial grid so that integrals become sums, building the bridge from $L^2$ to the finite vector a computer can hold?
- 1.8 [BSc] How do I show that an overall constant phase is unobservable, while the local phase gradient $\partial_x\arg\psi$ carries the current (the velocity field)?
- 1.9 [BSc] How does a Galilean boost transform a wave packet (the boost phase factor), and why is $|\psi|^2$ unchanged in shape?

## Part 2. Operators, position and momentum, and the Fourier bridge
- 2.1 [BSc] How do I represent the position operator $\hat x$ as multiplication and the momentum operator $\hat p=-i\,d/dx$ as differentiation, acting on a test function?
- 2.2 [BSc] How do I verify the canonical commutator $[\hat x,\hat p]=i$ by letting it act on an arbitrary function?
- 2.3 [BSc] How do I compute the momentum expectation $\langle p\rangle$ and spread $\Delta p$ in the position representation?
- 2.4 [BSc] How do I show that $\hat p$ is Hermitian by integration by parts, and see that the boundary terms vanish for an $L^2$ state?
- 2.5 [BSc] How do I write the momentum eigenfunctions $e^{ipx}$, see that they are not square-integrable, and normalize them to a Dirac delta?
- 2.6 [BSc] How do I obtain the momentum-space wavefunction $\tilde\psi(p)$ by Fourier transform and compute $\langle p^n\rangle$ as a moment in $p$?
- 2.7 [BSc] How do I verify the position-momentum uncertainty relation $\Delta x\,\Delta p\ge 1/2$ and show that a Gaussian saturates it?
- 2.8 [MSc] How do I build a function $f(\hat x,\hat p)$ of noncommuting operators, and resolve the ordering ambiguity by Weyl symmetric ordering?
- 2.9 [MSc] How do I distinguish self-adjointness from formal Hermiticity through the operator domain and boundary conditions (the box versus the half-line)?

## Part 3. The time-independent Schrodinger equation as an ODE eigenvalue problem
- 3.1 [BSc] How do I write the stationary Schrodinger equation $-\tfrac12\psi''+V\psi=E\psi$, and what distinguishes a normalizable bound state from a scattering state?
- 3.2 [BSc] How do I solve the infinite square well by imposing boundary conditions on the ODE, and read off the spectrum and eigenfunctions?
- 3.3 [BSc] How do I verify the orthonormality and completeness of the well eigenfunctions?
- 3.4 [BSc] How do I find the finite-square-well bound states from the even/odd transcendental quantization conditions with `FindRoot`?
- 3.5 [BSc] How do I compute numerical bound states with `NDEigensystem` on a finite interval and compare them to the analytic answer?
- 3.6 [BSc] How do I implement the shooting and Numerov methods from scratch to find a bound state?
- 3.7 [BSc] How do I show that the eigenstates of a symmetric potential have definite parity?
- 3.8 [MSc] How do I verify the node theorem, that the $n$-th bound state has exactly $n$ nodes?
- 3.9 [BSc] How do I find the single bound state of an attractive delta-function potential?
- 3.10 [MSc] How do I solve the linear (triangular) potential in terms of Airy functions, the quantum bouncer?

## Part 4. The harmonic oscillator
- 4.1 [BSc] How do I solve the oscillator's stationary equation and obtain the Hermite-function eigenstates with energies $E_n=n+\tfrac12$?
- 4.2 [BSc] How do I build the ladder operators $a=(\hat x+i\hat p)/\sqrt2$ and $a^\dagger$ as differential operators and verify $[a,a^\dagger]=1$?
- 4.3 [BSc] How do I generate the whole spectrum algebraically from $a|0\rangle=0$ and the number operator $\hat n=a^\dagger a$?
- 4.4 [BSc] How do I compute the matrix elements of $\hat x$ and $\hat p$ in the number basis, giving the truncated banded-matrix representation of the oscillator (with the truncation made explicit)?
- 4.5 [BSc] How do I show the ground state is a Gaussian, and restore $\omega$ and $m$ to read off the oscillator length $\sqrt{\hbar/m\omega}$?
- 4.6 [MSc] How do I build a coherent state $|\alpha\rangle$ as the eigenstate of $a$, and equivalently by applying $e^{\alpha a^\dagger-\alpha^* a}$ to the ground state, and show it is a minimum-uncertainty state?
- 4.7 [MSc] How do I build a squeezed state by applying $e^{\frac12(\xi^* a^2-\xi a^{\dagger2})}$ to the ground state and compute its unequal quadrature variances?
- 4.8 [MSc] How do I compute the photon-number statistics of a number state and a coherent state (the Poisson distribution)?

## Part 5. The time-dependent Schrodinger equation as a PDE
- 5.1 [BSc] How do I show that a stationary state evolves only by the phase $e^{-iEt}$, while a superposition acquires genuine time dependence?
- 5.2 [BSc] How do I integrate the time-dependent Schrodinger equation directly as a PDE with `NDSolve` for a given initial $\psi(x,0)$?
- 5.3 [BSc] How do I follow a free Gaussian wave packet's spreading and group velocity, both analytically and by `NDSolve`?
- 5.4 [MSc] How do I implement the split-step Fourier propagator from scratch and benchmark it against `NDSolve`?
- 5.5 [MSc] How do I evolve a state by expanding it in energy eigenstates and propagating term by term, and observe wave-packet revivals?
- 5.6 [MSc] How does an oscillator coherent state evolve, so that $\langle x\rangle(t)$ traces the classical oscillation while the packet stays minimal and does not spread?
- 5.7 [MSc] How do I watch a wave packet scatter off a barrier in real time and see tunneling and reflection split the packet?
- 5.8 [MSc] How do I construct the propagator (kernel) $K(x,t;x',0)$ for the free particle and for the oscillator (the Mehler kernel)?
- 5.9 [BSc] How do I verify Ehrenfest's theorem numerically, that $\langle x\rangle$ and $\langle p\rangle$ obey the classical equations of motion?
- 5.10 [BSc] How do I confirm that norm and probability current are conserved along an `NDSolve` evolution (a numerical-fidelity check)?

## Part 6. Scattering in one dimension
- 6.1 [BSc] How do I compute the reflection and transmission coefficients at a potential step?
- 6.2 [BSc] How do I compute the tunneling transmission $T(E)$ through a rectangular barrier?
- 6.3 [MSc] How do I exhibit transmission resonances in a well (the Ramsauer-Townsend effect)?
- 6.4 [MSc] How do I assemble the transfer-matrix method for a piecewise-constant potential?
- 6.5 [MSc] How do I build the one-dimensional scattering matrix, verify its unitarity, and identify bound states as poles?
- 6.6 [MSc] How do I check Levinson's theorem relating the phase shift to the number of bound states?

## Part 7. Approximation methods
- 7.1 [BSc] How do I compute first- and second-order energy shifts by nondegenerate time-independent perturbation theory (the anharmonic oscillator)?
- 7.2 [MSc] How do I resolve a degeneracy by degenerate perturbation theory?
- 7.3 [BSc] How do I get a Rayleigh-Ritz variational upper bound on the ground-state energy from a Gaussian trial function?
- 7.4 [MSc] How do I run the linear variational (Ritz) method in a finite basis, reducing the problem to a generalized matrix eigenproblem?
- 7.5 [BSc] How do I quantize a smooth well by the WKB (Bohr-Sommerfeld) condition?
- 7.6 [MSc] How do I compute a WKB tunneling rate through a barrier (the Gamow factor)?
- 7.7 [MSc] How do I compute transition rates by time-dependent perturbation theory and recover Fermi's golden rule with a continuum density of states?
- 7.8 [MSc] How do I apply the sudden and the adiabatic approximations to a continuous system?
- 7.9 [MSc] How do I match WKB solutions across a turning point with the Airy connection formulas?
- 7.10 [MSc] How do I get a variational *lower* bound on the ground-state energy (Temple's inequality)?

## Part 8. General theorems and structural methods
- 8.1 [BSc] How do I verify the virial theorem $2\langle T\rangle=\langle \vec r\cdot\nabla V\rangle$ on a stationary state?
- 8.2 [MSc] How do I apply the Hellmann-Feynman theorem to get $\partial E/\partial\lambda$ from an expectation value?
- 8.3 [MSc] How do I verify the Thomas-Reiche-Kuhn oscillator-strength sum rule?
- 8.4 [MSc] How do I exhibit the correspondence principle, the classical limit of stationary states at large quantum number?
- 8.5 [MSc] How do I factorize a Hamiltonian as $H=A^\dagger A+E_0$ and build its supersymmetric partner potential?
- 8.6 [MSc] How do I use shape invariance to obtain the oscillator and Poschl-Teller spectra purely algebraically?

## Part 9. Periodic potentials and band structure
- 9.1 [BSc] How do I state and use Bloch's theorem, writing eigenstates as $e^{ikx}u_k(x)$ with periodic $u_k$?
- 9.2 [MSc] How do I solve the Kronig-Penney model and plot the allowed and forbidden energy bands?
- 9.3 [MSc] How do I treat a Dirac comb (a lattice of delta potentials) and find its band edges?
- 9.4 [MSc] How do I compute the density of states within a band?
- 9.5 [MSc] How do I take the tight-binding limit and build Wannier functions from Bloch states?
- 9.6 [MSc] How do I exhibit Bloch oscillations (a Wannier-Stark ladder) under a static field?

## Part 10. Two and three dimensions: separation of variables
- 10.1 [BSc] How do I separate variables in a two- or three-dimensional box and form the product eigenfunctions?
- 10.2 [BSc] How do I find and explain the degeneracies of a square or cubic box?
- 10.3 [BSc] How do I solve the two-dimensional isotropic oscillator and count its degeneracies?
- 10.4 [BSc] How do I separate the two-body problem into center-of-mass and relative motion, reducing it to a one-body problem with the reduced mass $\mu$?
- 10.5 [BSc] How do I separate the Schrodinger equation in spherical coordinates into a radial and an angular equation?
- 10.6 [BSc] How do I build the effective radial potential and read off the centrifugal barrier?
- 10.7 [MSc] How do I separate variables in parabolic coordinates (the natural frame for the Stark and Coulomb problems)?

## Part 11. Orbital angular momentum in continuous space
- 11.1 [BSc] How do I write the orbital angular-momentum operators $L_x,L_y,L_z$ as differential operators on the angles?
- 11.2 [BSc] How do I verify $[L_i,L_j]=i\epsilon_{ijk}L_k$ and that $L^2$ is the Casimir with eigenvalue $l(l+1)$?
- 11.3 [BSc] How do I build the spherical harmonics $Y_{lm}$ as the simultaneous eigenfunctions of $L^2$ and $L_z$ and check their orthonormality?
- 11.4 [MSc] How do the raising and lowering operators $L_\pm$ act on the $Y_{lm}$?
- 11.5 [MSc] How do I rotate a wavefunction and represent the rotation on the $Y_{lm}$ by the Wigner $D$-matrix?
- 11.6 [MSc] How do I couple two angular momenta and compute the Clebsch-Gordan coefficients (with the three-$Y_{lm}$ Gaunt integral as the orbital instance), the change of basis reused later for adding spin?

## Part 12. Central potentials and the hydrogen atom
- 12.1 [BSc] How do I solve the Coulomb radial equation and obtain the bound-state energies $E_n=-1/(2n^2)$?
- 12.2 [BSc] How do I build the radial wavefunctions from the associated Laguerre polynomials and normalize them?
- 12.3 [BSc] How do I assemble the full $\psi_{nlm}$, plot probability densities and the radial distribution, and compute $\langle r\rangle$?
- 12.4 [BSc] How do I find the bound states of a spherical square well numerically?
- 12.5 [MSc] How do I solve the three-dimensional isotropic oscillator and reconcile its Cartesian and spherical degeneracy counts?
- 12.6 [MSc] How do I expose the extra Coulomb degeneracy through the conserved Runge-Lenz vector (the dynamical symmetry)?
- 12.7 [MSc] How do I compute central-potential eigenstates by applying `NDEigensystem` to the radial equation?
- 12.8 [MSc] How do I compute the quantum defect of an alkali-like screened-Coulomb potential?
- 12.9 [MSc] How do I apply WKB to the radial equation with the Langer correction $l(l+1)\to(l+\tfrac12)^2$ and recover the Coulomb and oscillator spectra?

## Part 13. Electromagnetic coupling
- 13.1 [MSc] How do I implement minimal coupling $\hat p\to\hat p-A$ and show how $\psi$ transforms under a gauge change?
- 13.2 [MSc] How do I derive the Landau levels of a charged particle in a uniform magnetic field?
- 13.3 [MSc] How do I compute the Aharonov-Bohm phase acquired around a flux line?
- 13.4 [MSc] How do I compute the normal (orbital) Zeeman splitting of hydrogen perturbatively from the $-\tfrac12 B L_z$ coupling?
- 13.5 [MSc] How do I compute the linear Stark effect on the degenerate $n=2$ hydrogen manifold?
- 13.6 [MSc] How do I compute Berry's geometric phase for a state carried adiabatically around a closed loop in parameter space?

## Part 14. Spin coupled to spatial motion
- 14.1 [MSc] How do I build a two-component spinor wavefunction $\psi_\sigma(x)$ by tensoring a spin-1/2 with the spatial state?
- 14.2 [MSc] How do I add the spin-orbit term $\propto \vec L\cdot\vec S$ and compute the hydrogen fine structure?
- 14.3 [MSc] How do I form the total angular momentum $\vec J=\vec L+\vec S$ and change to the coupled basis?
- 14.4 [MSc] How do I model Stern-Gerlach as a spin-dependent spatial deflection of a spinor wave packet in a field gradient?
- 14.5 [MSc] How do I solve the Pauli equation for a nonrelativistic spin in an electromagnetic field?
- 14.6 [MSc] How do I compute the anomalous Zeeman effect and the Paschen-Back crossover, where spin-orbit and the magnetic field compete?

## Part 15. Identical particles in continuous space
- 15.1 [MSc] How do I build a two-particle wavefunction $\psi(x_1,x_2)$ and symmetrize or antisymmetrize it?
- 15.2 [MSc] How do I exhibit the exchange hole for fermions and the bunching for bosons (Pauli exclusion versus enhancement)?
- 15.3 [MSc] How do I build a Slater determinant for $N$ fermions in a well or oscillator?
- 15.4 [MSc] How do I compute the exchange energy and the para/ortho splitting of a two-electron (helium-like) model?
- 15.5 [MSc] How do I solve the Hartree-Fock self-consistent-field equations for a model atom numerically?

## Part 16. Density operators, mixed states, and the Wigner function
- 16.1 [MSc] How do I represent a density operator as a kernel $\rho(x,x')$ and compute its purity $\iint|\rho(x,x')|^2\,dx\,dx'$?
- 16.2 [MSc] How do I obtain a reduced density matrix by tracing out one particle (integrating over its coordinate)?
- 16.3 [MSc] How do I compute the Wigner quasi-probability function by the Wigner-Weyl transform?
- 16.4 [MSc] How do I compute the Wigner function of coherent, number, and cat states and exhibit its negativity?
- 16.5 [MSc] How do I compute the Husimi-$Q$ and Glauber-Sudarshan-$P$ representations?
- 16.6 [MSc] How do I build the thermal (Gibbs) oscillator state, compute its Wigner function, and evolve a Wigner function by the Moyal equation?

## Part 17. Continuous-variable quantum optics and information
- 17.1 [MSc] How do I define the quadrature operators and the optical phase space of a single mode?
- 17.2 [MSc] How do I revisit the displacement and squeeze operators of Part 4 as phase-space transformations and read off their action on the quadratures?
- 17.3 [MSc] How do I model a beam splitter on two modes and exhibit Hong-Ou-Mandel interference?
- 17.4 [MSc] How do I build the two-mode squeezed vacuum and exhibit its EPR (position-momentum) correlations?
- 17.5 [MSc] How do I detect continuous-variable entanglement with the Duan separability criterion?
- 17.6 [MSc] How do I solve the Jaynes-Cummings model (one mode plus a two-level atom, truncated Fock) and exhibit vacuum Rabi oscillations?
- 17.7 [MSc] How do I model balanced homodyne and heterodyne detection, measuring a quadrature and reconstructing the Husimi-$Q$ distribution operationally?

## Part 18. Open quantum systems in continuous space
- 18.1 [MSc] How do I solve the Lindblad master equation for a damped oscillator?
- 18.2 [MSc] How do I model quantum Brownian motion (Caldeira-Leggett) on a truncated basis?
- 18.3 [MSc] How does dephasing select a pointer basis and decohere a spatial superposition (einselection)?
- 18.4 [MSc] How do I cast the optical master equation as a Fokker-Planck equation in the Wigner representation?
- 18.5 [MSc] How do I unravel a damped mode into quantum trajectories (the quantum-jump picture)?

## Part 19. Path integrals
- 19.1 [MSc] How do I derive the free-particle propagator from the Gaussian path integral?
- 19.2 [MSc] How do I get the oscillator propagator (Mehler) from the path integral?
- 19.3 [MSc] How do I use the imaginary-time path integral to compute the partition function of a finite system?
- 19.4 [MSc] How do I apply the stationary-phase (semiclassical) approximation and identify the classical action?
- 19.5 [MSc] How do I evaluate a discretized path integral numerically for a simple potential?
- 19.6 [MSc] How do I compute the double-well tunneling splitting from an instanton (imaginary-time bounce)?

## Part 20. Three-dimensional scattering theory
- 20.1 [MSc] How do I expand a scattering state in partial waves and extract the phase shifts $\delta_l$?
- 20.2 [MSc] How do I assemble the scattering amplitude and the differential cross section?
- 20.3 [MSc] How do I compute a cross section in the Born approximation?
- 20.4 [MSc] How do I extract the low-energy scattering length and verify the optical theorem?
- 20.5 [MSc] How do I fit a resonance to the Breit-Wigner form and check Levinson's theorem in three dimensions?
- 20.6 [MSc] How do I treat Coulomb scattering (Rutherford) with the Coulomb wave functions?

## Part 21. Nonlinear and mean-field wave mechanics
- 21.1 [MSc] How do I find the ground state of the Gross-Pitaevskii equation by imaginary-time evolution?
- 21.2 [MSc] How do I build the bright and dark solitons of the nonlinear Schrodinger equation?
- 21.3 [MSc] How do I apply the Thomas-Fermi approximation to a trapped condensate?
- 21.4 [MSc] How do I solve the Hartree mean-field equation as a self-consistent nonlinear Schrodinger problem?

## Part 22. Relativistic wave equations
- 22.1 [MSc] How do I solve the Klein-Gordon equation, write its plane-wave modes, and confront the negative-energy and indefinite-density problem?
- 22.2 [MSc] How do I build the Dirac equation from gamma matrices and write the free spinor solutions?
- 22.3 [MSc] How do I take the nonrelativistic limit to the Pauli equation and recover the $g=2$ prediction?
- 22.4 [MSc] How do I obtain the Dirac fine-structure spectrum of hydrogen?
- 22.5 [MSc] How do I exhibit the Klein paradox and zitterbewegung numerically?

## Part 23. From one particle to fields: the second-quantization bridge
- 23.1 [MSc] How do I build multimode occupation-number (Fock) space from copies of the oscillator ladder?
- 23.2 [MSc] How do I define the field operators $\hat\psi(x)$ as mode sums and impose the (anti)commutation relations?
- 23.3 [MSc] How do I represent a free scalar field as infinitely many oscillators and read off the vacuum energy?
- 23.4 [MSc] How do I show the second-quantized many-body Hamiltonian equals the symmetrized first-quantized one?
- 23.5 [MSc] How do I put a bosonic field on a finite lattice (truncated) and compute its dispersion relation?
- 23.6 [MSc] How does a time-dependent boundary or parameter create particles from the vacuum (a dynamical-Casimir-type capstone, with the truncation made explicit)?

---

## Coverage map (dependency-correct order, with the BSc/MSc split)

| Part | Theme | Items | BSc | MSc |
|------|-------|-------|-----|-----|
| 1 | The wavefunction and $L^2$ | 9 | 9 | 0 |
| 2 | Operators, position/momentum, Fourier bridge | 9 | 7 | 2 |
| 3 | TISE as an ODE eigenproblem | 10 | 8 | 2 |
| 4 | The harmonic oscillator | 8 | 5 | 3 |
| 5 | TDSE as a PDE | 10 | 5 | 5 |
| 6 | Scattering in one dimension | 6 | 2 | 4 |
| 7 | Approximation methods | 10 | 3 | 7 |
| 8 | General theorems and structural methods | 6 | 1 | 5 |
| 9 | Periodic potentials and band structure | 6 | 1 | 5 |
| 10 | Two and three dimensions: separation | 7 | 6 | 1 |
| 11 | Orbital angular momentum | 6 | 3 | 3 |
| 12 | Central potentials and hydrogen | 9 | 4 | 5 |
| 13 | Electromagnetic coupling | 6 | 0 | 6 |
| 14 | Spin coupled to spatial motion | 6 | 0 | 6 |
| 15 | Identical particles in continuous space | 5 | 0 | 5 |
| 16 | Density operators, mixed states, Wigner | 6 | 0 | 6 |
| 17 | Continuous-variable optics and information | 7 | 0 | 7 |
| 18 | Open systems in continuous space | 5 | 0 | 5 |
| 19 | Path integrals | 6 | 0 | 6 |
| 20 | Three-dimensional scattering theory | 6 | 0 | 6 |
| 21 | Nonlinear and mean-field wave mechanics | 4 | 0 | 4 |
| 22 | Relativistic wave equations | 5 | 0 | 5 |
| 23 | Second-quantization bridge | 6 | 0 | 6 |
| **Total** | | **158** | **54** | **104** |

A two-semester reading worked strictly in order. The 54 [BSc] items (Parts 1 through 12, the
single-particle wave mechanics core: the wavefunction, the position and momentum operators,
the time-independent equation and its bound states, the harmonic oscillator, the
time-dependent equation and wave packets, one-dimensional scattering, the elementary
approximation methods, separation of variables, orbital angular momentum, and the hydrogen
atom) form a self-contained first course, each answerable from earlier items. The 104 [MSc]
items extend it to a graduate course: band structure, electromagnetic coupling, spin and
identical particles, mixed states and the Wigner function, continuous-variable optics, open
systems, path integrals, three-dimensional scattering, nonlinear mean-field theory, the
relativistic wave equations, and the bridge to quantum field theory.
