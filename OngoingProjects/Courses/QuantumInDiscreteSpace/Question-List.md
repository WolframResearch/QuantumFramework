# Quantum in Finite Dimensions: A Comprehensive Question List (BSc through MSc, done in the Wolfram Language)

A curriculum-as-questions for learning to *do* quantum mechanics and quantum information by
computing it, with one hard constraint: every object lives in a finite-dimensional Hilbert space.
Each question is a thing you build, compute, or verify in the Wolfram Language. Questions are
tagged **[BSc]** (foundational, expected in a strong undergraduate course) or **[MSc]** (advanced,
graduate).

The order is a valid prerequisite chain: a question never needs a concept or skill that a later
question introduces. You can work straight down. This is why the density operator, the tensor
product, and elementary circuits appear early (they are used everywhere), the commutator sits with
the observables, and quantum channels precede the entropic information theory that needs them.

Each question is one concept. Where a single computation has several facets (an inner product and
its modulus and square) or one idea is shown across systems (degrees of freedom of qubit, qudit,
and mixed state), those stay together; genuinely distinct concepts are separate questions.

Genuinely infinite-dimensional material (bosonic Fock space, coherent and squeezed states of light,
continuous-variable systems, the harmonic oscillator) is deliberately excluded. Its
finite-dimensional counterparts, spin coherent states and the discrete Wigner function, appear in
its place. Every driven system uses a classical c-number field, never a quantized mode.

---

## Part 1. Pure states and the Born rule
- 1.1 [BSc] How do I represent a qubit as a state vector $\{\alpha,\beta\}$ in $\mathbb{C}^2$?
- 1.2 [BSc] How do I read the Born-rule probabilities of a state?
- 1.3 [BSc] How do I check that a state is normalized, and normalize one that is not?
- 1.4 [BSc] How do I show that a global phase is physically unobservable (two vectors differing by a phase give identical predictions)?
- 1.5 [BSc] How do I count the real parameters (degrees of freedom) of a pure qubit and a pure qudit in $\mathbb{C}^d$?
- 1.6 [BSc] How do I compute an inner product $\langle\phi|\psi\rangle$, its modulus, and the transition probability $|\langle\phi|\psi\rangle|^2$?
- 1.7 [BSc] How do I write a state both as a column vector of amplitudes in the computational basis and as a Dirac-notation superposition $\sum_i c_i\,|i\rangle$?
- 1.8 [BSc] How do I get the Bloch vector of a qubit, and the Bloch-sphere angles?
- 1.9a [BSc] How do I prepare a qudit in the computational basis state $|0\rangle$?
- 1.9b [BSc] How do I prepare a uniform superposition $\tfrac1{\sqrt d}\sum_k|k\rangle$ in $\mathbb{C}^d$?

## Part 2. Observables, spectra, and approximation
- 2.1 [BSc] How do I compute a commutator and an anticommutator of two operators?
- 2.2 [BSc] How do I build the Pauli operators and verify their algebra (Levi-Civita product rule)?
- 2.3a [BSc] How do I test that an operator is Hermitian?
- 2.3b [BSc] Why must an observable be Hermitian (its eigenvalues, the measurement outcomes, must be real)?
- 2.4 [BSc] How do I compute eigenvalues and eigenvectors and form the spectral decomposition?
- 2.5 [BSc] How do I evaluate a function of an operator $f(\hat A)$ by functional calculus?
- 2.6 [BSc] How do I compute an expectation value and a variance on a pure state as $\langle\psi|\hat A|\psi\rangle$?
- 2.7 [BSc] How do I write the resolution of the identity in an orthonormal basis?
- 2.8 [BSc] How do I decide whether two observables are compatible (commute) and find a simultaneous eigenbasis?
- 2.9a [MSc] How do I construct a complete set of commuting observables (CSCO)?
- 2.9b [MSc] How do I label basis states by the joint spectrum of a CSCO?
- 2.10a [MSc] How do I define the Hilbert-Schmidt operator inner product $\langle A,B\rangle=\mathrm{Tr}[A^\dagger B]$ and show the Pauli operators are orthogonal?
- 2.10b [MSc] How do I expand an operator in the Pauli basis?
- 2.11 [BSc] How do I bound the ground-state energy from above using a trial state (the Rayleigh-Ritz variational method)?
- 2.12 [MSc] How do I compute first- and second-order energy shifts by non-degenerate time-independent perturbation theory?
- 2.13 [MSc] How do I resolve an avoided crossing by degenerate perturbation theory, obtaining the split levels $E_\pm=\bar E\pm\sqrt{\delta^2+|V|^2}$ (level repulsion)?

## Part 3. States as density operators
- 3.1a [BSc] How do I build a density matrix from a classical ensemble?
- 3.1b [BSc] How do I compute the purity $\mathrm{Tr}[\rho^2]$ of a density matrix?
- 3.2 [BSc] How do I test whether a state is pure or mixed ($\mathrm{Tr}[\rho^2]=1$, rank one, $\rho^2=\rho$)?
- 3.3 [BSc] How do I build the maximally mixed state $I/d$ and read its purity?
- 3.4 [BSc] How do I compute an expectation value and a variance of an observable on a general (mixed) state as $\mathrm{Tr}[\hat A\rho]$?
- 3.5 [BSc] How do I count the real parameters (degrees of freedom) of a $d$-dimensional mixed state (density matrix)?
- 3.6 [BSc] How do I map a qubit density matrix to its Bloch vector and read off positivity (radius $\le 1$) and the pure-state boundary (radius $=1$)?
- 3.7 [MSc] How do I show that the state space is convex, with pure states as its extreme points (the Bloch ball as the $d=2$ case)?
- 3.8 [BSc] How do I compute the von Neumann entropy of a density matrix?
- 3.9 [BSc] How do I test whether a matrix is a valid density operator ($\rho=\rho^\dagger$, $\mathrm{Tr}\,\rho=1$, $\rho\succeq0$)?

## Part 4. Composite systems: tensor product and partial trace
- 4.1 [BSc] How do I form the tensor product of states and of operators?
- 4.2 [BSc] How do I take a partial trace of a two-party state and obtain a reduced density matrix (mixed when the input is entangled)?

## Part 5. Projective measurement
- 5.1 [BSc] How do I perform a projective measurement and read outcomes, probabilities, and the mean?
- 5.2 [BSc] How do I measure in a non-computational basis and get the post-measurement (collapsed) states?
- 5.3 [BSc] How do I simulate finite-shot statistics and watch the empirical frequencies approach the Born rule?
- 5.4 [BSc] How do I apply a non-selective projective measurement (the Lüders channel $\rho\mapsto\sum_k P_k\rho P_k$), which decoheres a state without selecting an outcome?

## Part 6. Uncertainty and incompatibility
- 6.1 [MSc] How do I verify the Robertson-Schrodinger uncertainty relation for a given state?
- 6.2 [MSc] How do I find a minimum-uncertainty state that saturates the Robertson bound for a chosen pair of observables?
- 6.3 [BSc] How do I compute the classical Shannon entropy $H(p)=-\sum_i p_i\log_2 p_i$ of a measurement distribution?
- 6.4 [MSc] How do I construct a maximal set of mutually unbiased bases for a qubit (bases with every cross-overlap $|\langle e_i|f_j\rangle|^2=1/d$; maximal = no further basis can be added that stays unbiased to every member)?
- 6.4b [BSc] How do I verify that given bases are pairwise mutually unbiased, i.e. every cross-overlap satisfies $|\langle e_i|f_j\rangle|^2 = 1/d$?
- 6.5 [MSc] How do I verify the Maassen-Uffink entropic uncertainty relation $H(A)+H(B)\ge-\log_2 c^2$ (with $H(A)$ the Shannon entropy of the outcome distribution of measuring $A$, and $c=\max_{j,k}|\langle a_j|b_k\rangle|$ the largest overlap between an eigenstate $|a_j\rangle$ of $A$ and an eigenstate $|b_k\rangle$ of $B$) for two observables?
- 6.6 [MSc] How do I find the states that saturate the Maassen-Uffink entropic bound $H(X)+H(Z)\ge 1$ bit for the complementary $X,Z$ (maximal complementarity)?

## Part 7. Unitary dynamics and pictures
- 7.1a [BSc] How do I build the time-evolution operator $U(t)=e^{-i\hat H t}$ (working in units where $\hbar=1$)?
- 7.1b [BSc] How do I solve the Schrodinger equation for a time-independent Hamiltonian, $|\psi(t)\rangle=e^{-i\hat H t}|\psi_0\rangle$?
- 7.2 [BSc] How do I find stationary states and show they evolve only by a phase?
- 7.3 [BSc] How do I compute Rabi oscillations of a two-level system under a resonant classical drive?
- 7.4 [BSc] How do I compute the Larmor precession of a spin in a static magnetic field?
- 7.5a [BSc] How do I evolve an observable in the Heisenberg picture, $A_H(t)=U^\dagger(t)\,A\,U(t)$?
- 7.5b [BSc] How do I check Ehrenfest's theorem, $\tfrac{d}{dt}\langle A\rangle=i\langle[\hat H,A]\rangle$, for a time-independent observable?
- 7.6 [MSc] How do I transform between the Schrodinger, Heisenberg, and interaction pictures and show they agree on expectation values?
- 7.7 [MSc] How do I derive the Mandelstam-Tamm quantum speed limit $t\geq\tfrac{\hbar\,\arccos|\langle\psi_0|\psi_t\rangle|}{\Delta E}$ from the time-energy uncertainty relation?
- 7.8 [MSc] How do I solve a classically driven two-level system in the rotating frame (the Rabi problem)?
- 7.9a [MSc] How do I compute the first-order transition probability under a periodic perturbation $H'(t)=V\cos(\omega t)$, $P_{i\to f}(t)\approx|V_{fi}|^2\,t^2\,\mathrm{sinc}^2\!\big(\tfrac{(\omega_{fi}-\omega)t}{2}\big)$?
- 7.9b [MSc] How do I recover Fermi's golden rule rate $\Gamma=2\pi|V_{fi}|^2\rho(E_f)$ from a finite dense band of final states in the dense-band limit?
- 7.10a [MSc] How do I find a symmetry operator $S$ that commutes with $\hat H$ and use its eigenvalue as a conserved quantum number?
- 7.10b [MSc] How do I block-diagonalize $\hat H$ into symmetry sectors using the eigenspaces of a commuting symmetry $S$?
- 7.11 [MSc] How do I implement a Trotter-Suzuki decomposition to simulate $e^{-i\hat H t}$ for a sum of noncommuting terms?
- 7.12 [MSc] How do I verify the adiabatic theorem's gap and rate condition (the transition probability vanishing when the sweep is slow compared to the inverse gap squared, $T \gg 1/\Delta^2$)?
- 7.13 [MSc] How do I compute Berry's geometric phase around a closed loop as the discrete Bargmann invariant $\gamma=-\arg\prod_{k}\langle\psi_k|\psi_{k+1}\rangle$ with $|\psi_{N}\rangle\equiv|\psi_0\rangle$?
- 7.14 [MSc] How do I simulate a quantum quench and compute the Loschmidt echo?
- 7.15 [MSc] How do I exhibit the quantum Zeno effect, freezing evolution under a Hamiltonian by frequent projective measurement?
- 7.16 [BSc] How do I encode each quantum postulate as a one-line WL operation (state $=$ unit vector or density operator, observable $=$ Hermitian matrix, evolution $=$ unitary, measurement $=$ projector and Born rule, composite $=$ tensor product)?

## Part 8. Spin and angular momentum
- 8.1 [BSc] How do I model a spin-1/2 in a magnetic field and reproduce the Stern-Gerlach spin-projection outcomes?
- 8.2 [BSc] How do I build the spin-$j$ angular-momentum operators and verify $[J_x,J_y]=iJ_z$ and the Casimir?
- 8.3 [BSc] How do I rotate a spin state and read the rotation off the expectation values?
- 8.4 [MSc] How do I construct the Wigner $D$-matrix element $D^j_{m'm}(\alpha,\beta,\gamma)=\langle j\,m'|e^{-i\alpha J_z}e^{-i\beta J_y}e^{-i\gamma J_z}|j\,m\rangle$ for a spin-$j$ rotation?
- 8.5 [MSc] How do I add two angular momenta and compute the Clebsch-Gordan coefficients?
- 8.6 [MSc] How do I build the singlet and triplet states and identify them by total spin?
- 8.7 [MSc] How do I construct a spin coherent state and place it on the generalized Bloch sphere?

## Part 9. Single-qubit operations and SU(2)
- 9.1a [BSc] How do I build the standard named single-qubit gates $X, Y, Z, H, S, T$?
- 9.1b [BSc] How do I build the axis-angle rotation gates $R_{\hat n}(\theta)=e^{-i\theta\,\hat n\cdot\vec\sigma/2}$?
- 9.2 [BSc] How do I derive the closed-form axis-angle exponential $\cos(\alpha/2)I - i\sin(\alpha/2)\,\hat n\cdot\vec\sigma$?
- 9.3 [BSc] How do I verify the unitarity and unit-determinant conditions that define SU(2)?
- 9.4 [BSc] How do I decompose an arbitrary single-qubit unitary into Euler ($Z$-$Y$-$Z$) angles?
- 9.5 [MSc] How do I show that SU(2) double-covers SO(3) and that a $2\pi$ rotation is not the identity?
- 9.6 [MSc] How do I realize a single-qubit unitary to within a target error $\epsilon$ as a sequence of $\{H, T\}$ gates (the Solovay-Kitaev idea)?

## Part 10. Elementary multi-qubit gates and circuits
- 10.1 [BSc] How do I build the standard multi-qubit gates (CNOT, CZ, SWAP, Toffoli, controlled-$U$)?
- 10.2a [BSc] How do I compose a circuit from gates and apply it to an input state?
- 10.2b [BSc] How do I read off the overall unitary of a composed circuit?
- 10.3 [BSc] How do I compute the depth and the parallel layers of a circuit?

## Part 11. Composite systems and entanglement
- 11.1 [BSc] How do I build the Bell, GHZ, and W states?
- 11.2 [BSc] How do I measure a Bell pair jointly and exhibit the perfect correlations?
- 11.3 [BSc] How do I compute the Schmidt decomposition and Schmidt rank of a bipartite pure state?
- 11.4 [BSc] How do I quantify pure-state entanglement by the entropy of entanglement (the von Neumann entropy of a reduced state)?
- 11.5 [BSc] How do I verify no-signalling through the invariance of reduced states?
- 11.6 [MSc] How do I test a mixed state for entanglement by the positive-partial-transpose (Peres-Horodecki) criterion?
- 11.7 [MSc] How do I compute the concurrence and the negativity of a mixed two-qubit state?
- 11.8 [MSc] How do I construct an entanglement witness $W$ for a given entangled state $\rho$, with $\mathrm{Tr}[W\sigma]\geq0$ for all separable $\sigma$ but $\mathrm{Tr}[W\rho]<0$?
- 11.9 [MSc] How do I decide whether pure state $|\psi\rangle$ can be converted to $|\phi\rangle$ by LOCC, using Nielsen's theorem that it is possible iff the Schmidt-coefficient vector of $|\psi\rangle$ is majorized by that of $|\phi\rangle$?
- 11.10 [MSc] How do I verify monogamy of entanglement (the CKW inequality) using the three-qubit tangle?
- 11.11 [MSc] How do I distinguish the GHZ and W classes of tripartite entanglement under SLOCC by computing the three-tangle $\tau_3$ (nonzero for GHZ, zero for W)?
- 11.12a [MSc] How do I construct a $3\times3$ (or $2\times4$) PPT-entangled (bound) state?
- 11.12b [MSc] Why can no bound-entangled state exist for two qubits (the Horodecki $2\times2$/$2\times3$ theorem)?
- 11.13 [MSc] How do I build the symmetric and antisymmetric subspaces of identical particles and project onto them (exchange symmetry, bosons and fermions in a finite mode set)?

## Part 12. Mixed states: distinguishability and thermal states
- 12.1 [BSc] How do I build a Werner state and locate its entanglement threshold?
- 12.2 [BSc] How do I purify a mixed state and recover it by partial trace?
- 12.3 [MSc] How do I show that different ensembles can give the same density matrix (the ensemble ambiguity / unitary freedom)?
- 12.4 [BSc] How do I build the Gibbs (thermal) state $e^{-\beta\hat H}/Z$ of a finite system and compute thermal averages?
- 12.5 [BSc] How do I compare two states by trace distance and (pure-state) fidelity?
- 12.6 [MSc] How do I compute the Uhlmann fidelity $F(\rho,\sigma)=(\mathrm{Tr}\sqrt{\sqrt{\rho}\,\sigma\sqrt{\rho}})^2$ between two mixed states?
- 12.7 [MSc] How do I find the optimal state-discrimination measurement and the Helstrom error probability $P_{\mathrm{err}}=\tfrac12\big(1-\|p_1\rho_1-p_2\rho_2\|_1\big)$?
- 12.8a [MSc] How do I verify Gleason's frame-function constraint for a qutrit ($d\geq3$), that every frame function is $\langle\psi|\rho|\psi\rangle$ for some density operator $\rho$?
- 12.8b [MSc] Why does the qubit ($d=2$) fail Gleason and instead need the POVM (Busch) version?

## Part 13. Quantum channels and open-system dynamics
- 13.1 [BSc] How do I build a channel from Kraus operators and apply it to a state?
- 13.2 [BSc] How do I apply the standard noise channels (bit flip, phase flip, depolarizing, amplitude damping)?
- 13.3 [BSc] How do I get the Choi matrix and certify complete positivity?
- 13.4 [BSc] How do I verify trace preservation ($\sum_k K_k^\dagger K_k = I$)?
- 13.5 [MSc] How do I build the Stinespring dilation (an isometry into system plus environment)?
- 13.6 [MSc] How do I convert between the Kraus, Choi, Stinespring, and superoperator representations?
- 13.7 [MSc] How do I construct the complementary channel $\mathcal{N}^c(\rho)=\mathrm{Tr}_{\mathrm{sys}}[V\rho V^\dagger]$ by tracing the system out of the Stinespring isometry, leaving the environment output?
- 13.8 [MSc] How do I test whether a channel is unital, and what does that mean for the Bloch ball?
- 13.9 [MSc] How do I compose channels and test CP-divisibility ($\mathcal{N}_{t_2,t_0}=\mathcal{N}_{t_2,t_1}\circ\mathcal{N}_{t_1,t_0}$ with each intermediate map completely positive) as the criterion for Markovianity?
- 13.10 [BSc] How do I solve the Lindblad master equation for $T_1$ relaxation and $T_2$ dephasing?
- 13.11 [MSc] How do I solve the optical Bloch equations for a driven, damped two-level atom (interaction picture plus Lindblad damping)?
- 13.12 [BSc] How do I find the steady state of a dissipative open system?
- 13.13 [MSc] How do I unravel the Lindblad equation into quantum trajectories (the quantum-jump picture)?
- 13.14 [MSc] How does local noise destroy entanglement, and can it die at finite time (sudden death)?
- 13.15 [MSc] How do I bound the distinguishability of two channels by the diamond norm $\|\mathcal{N}_1-\mathcal{N}_2\|_\diamond=\sup_\rho\|(\mathcal{N}_1\otimes\mathrm{id})(\rho)-(\mathcal{N}_2\otimes\mathrm{id})(\rho)\|_1$?
- 13.16 [MSc] How does a dephasing channel select a pointer basis and suppress off-diagonal coherences (decoherence and einselection)?

## Part 14. Quantum information and entropy
- 14.1 [MSc] How do I compute and compare the Shannon and von Neumann entropies?
- 14.2 [MSc] How do I compute the quantum relative entropy and check its nonnegativity?
- 14.3 [MSc] How do I compute joint, conditional, and mutual quantum entropies (and find a negative conditional entropy)?
- 14.4 [MSc] How do I compute the coherent information $I_c(\rho,\mathcal{N})=S(\mathcal{N}(\rho))-S((\mathcal{N}\otimes\mathrm{id})(|\psi_\rho\rangle\langle\psi_\rho|))$ of a state through a channel?
- 14.5 [MSc] How do I verify the data-processing inequality under a channel?
- 14.6 [MSc] How do I compute the Holevo bound $\chi=S\!\big(\sum_i p_i\rho_i\big)-\sum_i p_i\,S(\rho_i)$ on the accessible information of an ensemble?
- 14.7 [MSc] How do I assemble a channel-capacity estimate for a simple channel, the classical (Holevo) capacity or the quantum (coherent-information) capacity, from the entropic quantities above?
- 14.8 [MSc] How do I quantify the coherence of a state in a fixed basis by the relative-entropy $C_{\mathrm{rel}}(\rho)=S(\rho_{\mathrm{diag}})-S(\rho)$ or $\ell_1$ measure $C_{\ell_1}(\rho)=\sum_{i\neq j}|\rho_{ij}|$?

## Part 15. Two-qubit decomposition and universality
- 15.1 [MSc] How do I decompose a two-qubit gate (KAK / Cartan) and count the CNOTs it needs?
- 15.2 [MSc] How do I confirm a gate set is universal by exactly recompiling a random unitary?
- 15.3 [MSc] How do I express a Clifford$+T$ circuit and count its non-Clifford ($T$) cost (see also magic in 19.6)?

## Part 16. Generalized measurement and dilations
- 16.1 [MSc] How do I build the qubit (tetrahedral) SIC-POVM, check its effects are positive and complete, and apply it?
- 16.2 [MSc] How do I realize a POVM as a projective measurement on a larger space (Naimark dilation)?
- 16.3 [MSc] How do I describe a measurement by its Kraus operators (a quantum instrument) and get outcome-conditioned states?
- 16.4 [MSc] How do I compute a weak value as a pre- and post-selected matrix element (with a finite two-level pointer) and see it exceed the eigenvalue range?

## Part 17. Tomography, estimation, and metrology
- 17.1 [MSc] How do I reconstruct an unknown density matrix by state tomography?
- 17.2 [MSc] How do I reconstruct an unknown channel by process tomography?
- 17.3 [MSc] How do I compute the quantum Fisher information and the Cramer-Rao bound for a phase-estimation problem?
- 17.4 [MSc] How do I estimate average gate fidelity by randomized benchmarking, fitting the survival probability to $A p^m + B$ and reading the depolarizing rate from the decay $p$?

## Part 18. Qudits and higher-dimensional systems
- 18.1 [BSc] How do I build the clock and shift (Weyl-Heisenberg) operators and verify the Weyl relation?
- 18.2 [BSc] How do I build the qutrit entangling gate (the SUM gate, a qudit CNOT)?
- 18.3 [BSc] How do I build the qudit Fourier transform?
- 18.4 [MSc] How do I construct the generalized Pauli (Heisenberg-Weyl) group for a qudit?
- 18.5 [BSc] How do I build mutually unbiased bases for a general qudit (via the Fourier basis)?
- 18.6 [MSc] How do I construct a SIC-POVM in general dimension $d$ and verify its equiangular symmetry?

## Part 19. The stabilizer formalism
- 19.1 [MSc] How do I represent a stabilizer state by its generators (tableau)?
- 19.2 [MSc] How do I describe the Pauli and Clifford groups as groups (order and generators), with the Clifford group as the normalizer of the Pauli group, $C P C^\dagger=P$?
- 19.3 [MSc] How do I apply Clifford gates and watch the generators transform?
- 19.4 [MSc] How do I simulate a large Clifford circuit efficiently (Gottesman-Knill)?
- 19.5 [MSc] How do I build a graph (cluster) state and read its stabilizers?
- 19.6 [MSc] How do I detect the non-stabilizerness (magic) of the state produced when a non-Clifford $T$ gate acts on a stabilizer state?

## Part 20. Quantum error correction
- 20.1 [BSc] How do I encode the three-qubit bit-flip code and correct an error by syndrome measurement?
- 20.2 [MSc] How do I state and check the Knill-Laflamme error-correction conditions for a code?
- 20.3 [MSc] How do I build the nine-qubit Shor code and correct an arbitrary single-qubit error?
- 20.4 [MSc] How do I build a classical linear code from its generator and parity-check matrices ($G$, $H$ with $H G^\top = 0$ over $\mathrm{GF}(2)$) and compute its minimum distance, as the foundation for CSS codes?
- 20.5 [MSc] How do I build a CSS code from two classical codes, with the Steane $[[7,1,3]]$ as the worked instance?
- 20.6 [MSc] How do I build the five-qubit perfect code's stabilizers and correct an arbitrary single-qubit error with it (the smallest, non-CSS, code, saturating the quantum Hamming bound)?
- 20.7 [MSc] How do I compute a code's distance?
- 20.8 [MSc] How do I identify a code's logical operators?
- 20.9 [MSc] How do I lay out a surface-code patch and read its stabilizers and logical qubit?
- 20.10a [MSc] How do I demonstrate a transversal (fault-tolerant) logical gate by applying $U^{\otimes n}$ and verifying it acts as the logical gate on codewords?
- 20.10b [MSc] What does the threshold theorem guarantee (arbitrarily reliable computation once the physical error rate is below a threshold, $p<p_{\mathrm{th}}$)?

## Part 21. Discrete phase space and magic as a resource
- 21.1 [MSc] How do I compute the discrete Wigner function of a qudit state for odd $d$?
- 21.2 [MSc] Why does even $d$ (in particular a qubit) need a different phase-space construction (Gibbons-Hoffman-Wootters), and how do I build it?
- 21.3 [MSc] How do I detect Wigner negativity as a marker of non-classicality?
- 21.4 [MSc] How do I compute the mana, a magic monotone, from the Wigner negativity?
- 21.5 [MSc] How do I estimate the stabilizer rank $\chi$ of a magic state, the least number of stabilizer states in a superposition equal to it?
- 21.6 [MSc] How do I run one round of the 5-to-1 (Bravyi-Kitaev) magic-state distillation protocol and watch the output fidelity improve?

## Part 22. Finite-dimensional many-body systems
- 22.1 [MSc] How do I build the transverse-field Ising Hamiltonian on a finite chain and diagonalize it?
- 22.2 [MSc] How do I map a spin-1/2 chain to free fermions by the Jordan-Wigner transformation $\sigma_j^+ = c_j^\dagger\prod_{k<j}(1-2c_k^\dagger c_k)$, making the transverse-field Ising and XY chains exactly solvable?
- 22.3 [MSc] How do I find the ground state and the energy gap of a Heisenberg or XY chain?
- 22.4 [MSc] How do I locate a finite-size quantum phase transition by the gap and an order parameter?
- 22.5 [MSc] How do I compute the entanglement entropy of a block and observe area-law scaling?
- 22.6 [MSc] How do I write a small state as a matrix product state and read its bond dimension?
- 22.7 [MSc] How do I evolve a spin chain and observe a Lieb-Robinson light cone?

## Part 23. Quantum algorithms
- 23.1 [BSc] Deutsch: how do I decide a one-bit function's type in a single query?
- 23.2 [BSc] Deutsch-Jozsa: how do I tell constant from balanced in one query?
- 23.3 [BSc] Bernstein-Vazirani: how do I recover a hidden string in one query?
- 23.4 [MSc] Simon: how do I find a hidden period and finish the classical post-processing?
- 23.5 [BSc] Grover: how do I amplitude-amplify a marked item, and how many iterations are optimal?
- 23.6 [MSc] How do I generalize Grover to amplitude amplification on an arbitrary initial state $A|0\rangle$ using $Q=-A S_0 A^\dagger S_\chi$?
- 23.7 [MSc] How do I estimate a marked fraction by amplitude estimation (phase estimation on the Grover operator)?
- 23.8 [BSc] How do I build and verify the quantum Fourier transform?
- 23.9 [MSc] How do I run quantum phase estimation and read out an eigenphase?
- 23.10 [MSc] Shor: how do I find the order of an element modulo $N$ by phase estimation?
- 23.11 [MSc] HHL: how do I solve $A\vec x=\vec b$ by eigenvalue inversion (phase estimation plus a controlled rotation)?
- 23.12 [MSc] How do I implement a discrete-time quantum walk on a finite graph by alternating a coin operator on the coin register with a shift operator $S=\sum_{c,v}|c\rangle\langle c|\otimes|v{+}c\rangle\langle v|$?
- 23.13 [MSc] VQE: how do I find a ground-state energy variationally?
- 23.14 [MSc] QAOA: how do I solve Max-Cut, and how does adding layers improve it?
- 23.15 [MSc] How do I compute a variational gradient exactly by the parameter-shift rule?
- 23.16 [MSc] How do I exhibit a barren plateau, the gradient variance $\mathrm{Var}[\partial_\theta\langle H\rangle]$ decaying exponentially, $\sim 2^{-n}$, with qubit number $n$?

## Part 24. Quantum communication protocols
- 24.1 [BSc] How do I teleport an unknown qubit?
- 24.2 [BSc] How do I send two classical bits with one qubit (superdense coding)?
- 24.3 [BSc] How do I entangle two qubits that never interacted (entanglement swapping)?
- 24.4 [MSc] How do I distill a better Bell pair from several noisy ones?
- 24.5 [MSc] How do I run the BB84 key-distribution protocol and detect an eavesdropper?

## Part 25. Foundations: nonlocality and contextuality
- 25.1 [BSc] No-cloning: how do I show that no unitary can copy an unknown state?
- 25.2 [BSc] No-deleting: how do I show that no unitary can erase one of two copies of an unknown state?
- 25.3 [MSc] No-broadcasting: how do I show that no channel can broadcast an unknown mixed state (the generalization of no-cloning)?
- 25.4 [BSc] CHSH: how do I compute the quantum violation and reach the Tsirelson bound?
- 25.5 [MSc] PR boxes: how do I show they beat quantum mechanics on CHSH ($4>2\sqrt2>2$), placing quantum correlations above the classical bound but below the no-signalling maximum?
- 25.6 [BSc] GHZ: how do I exhibit an all-or-nothing, single-shot refutation of local realism?
- 25.7 [MSc] How do I exhibit state-independent (Kochen-Specker-type) contextuality with the Mermin-Peres magic square?
- 25.8 [MSc] How do I verify a Kochen-Specker uncolorable ray set (for example the Peres 24-ray configuration)?
- 25.9 [MSc] How do I violate a Leggett-Garg (temporal Bell) inequality?
- 25.10 [MSc] Hardy: how do I exhibit nonlocality without inequalities?
- 25.11 [MSc] PBR: how do I build the antidistinguishing measurement on independent copies of non-orthogonal states that (under preparation independence) rules out a $\psi$-epistemic, knowledge-only, model?

---

## Coverage map (dependency-correct order, with the BSc/MSc split)

| Part | Theme | Items | BSc | MSc |
|------|-------|-------|-----|-----|
| 1 | Pure states and the Born rule | 10 | 10 | 0 |
| 2 | Observables, spectra, approximation | 16 | 10 | 6 |
| 3 | States as density operators | 10 | 9 | 1 |
| 4 | Composite systems: tensor product, partial trace | 2 | 2 | 0 |
| 5 | Projective measurement | 4 | 4 | 0 |
| 6 | Uncertainty and incompatibility | 7 | 2 | 5 |
| 7 | Unitary dynamics and pictures | 20 | 8 | 12 |
| 8 | Spin and angular momentum | 7 | 3 | 4 |
| 9 | Single-qubit operations, SU(2) | 7 | 5 | 2 |
| 10 | Elementary multi-qubit gates and circuits | 4 | 4 | 0 |
| 11 | Composite systems and entanglement | 14 | 5 | 9 |
| 12 | Mixed states: distinguishability, thermal | 9 | 4 | 5 |
| 13 | Channels and open systems | 16 | 6 | 10 |
| 14 | Quantum information and entropy | 8 | 0 | 8 |
| 15 | Two-qubit decomposition, universality | 3 | 0 | 3 |
| 16 | Generalized measurement and dilations | 4 | 0 | 4 |
| 17 | Tomography, estimation, metrology | 4 | 0 | 4 |
| 18 | Qudits | 6 | 4 | 2 |
| 19 | Stabilizer formalism | 6 | 0 | 6 |
| 20 | Error correction | 11 | 1 | 10 |
| 21 | Discrete phase space, magic | 6 | 0 | 6 |
| 22 | Finite many-body | 7 | 0 | 7 |
| 23 | Algorithms | 16 | 5 | 11 |
| 24 | Communication | 5 | 3 | 2 |
| 25 | Foundations | 11 | 4 | 7 |
| **Total** | | **213** | **89** | **124** |

A two-semester reading worked strictly in order. The 89 [BSc] items form a self-contained first
course in finite-dimensional quantum theory and quantum information, each answerable from earlier
items; the 124 [MSc] items extend it to a graduate course covering open systems, error correction,
many-body, resource theories, and the foundations.
