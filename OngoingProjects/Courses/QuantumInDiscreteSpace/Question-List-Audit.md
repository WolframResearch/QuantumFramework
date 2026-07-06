# Audit Review of `Question-List.md`

Source file fully read: `Question-List.md`  
Audit date: 2026-07-02  
Scope requested: BSc/MSc finite-dimensional quantum mechanics and quantum information, with questions that are self-contained, unambiguous, and centered on one concept.

## Executive Verdict

The list is ambitious and broadly valuable. It covers far more than a standard BSc quantum course and reaches a serious MSc quantum information/computation curriculum: density matrices, entanglement, channels, entropy, stabilizers, error correction, phase space, many-body systems, algorithms, communication, and foundations.

It is not yet ready under its own stated standard. The main issue is not missing advanced topics; it is that many entries are too compressed, some are not self-contained, and several are technically over-strong or ambiguous. A student could use the list as a map, but not yet as a clean sequence of atomic, answerable questions.

The coverage table is internally consistent: I count 207 bullet questions, with 85 tagged BSc and 122 tagged MSc. The stronger claim that the order is a valid prerequisite chain is only partly true.

## Main Strengths

- The finite-dimensional constraint is mostly respected. The file avoids harmonic oscillators, continuous-variable Hilbert spaces, and quantized modes, while retaining finite spin, qudit, channel, and many-body examples.
- The density-operator and tensor-product machinery appears early enough to support later quantum information topics.
- The list correctly includes modern MSc-level material that many older curricula omit: Choi matrices, Stinespring dilations, complementary channels, coherent information, stabilizer formalism, surface codes, magic, tomography, randomized benchmarking, and contextuality.
- The BSc/MSc split is useful as a first pass, even though several individual tags need correction.
- The list is written in an active computational style, which fits the stated Wolfram Language goal.

## Major Issues

### 1. Missing Foundation Layer

The list begins with qubit state vectors but never explicitly establishes the finite-dimensional linear-algebra conventions that all later questions depend on. For a self-contained BSc/MSc list, add a short Part 0 before Part 1.

Recommended additions:

- 0.1 [BSc] Given a finite orthonormal basis `{|0>,...,|d-1>}`, how do I represent kets, bras, and operators as vectors and matrices?
- 0.2 [BSc] How do I compute the adjoint of a vector or matrix, and what inner-product convention am I using?
- 0.3 [BSc] How do I test whether vectors form an orthonormal basis, and how do I orthonormalize a linearly independent set?
- 0.4 [BSc] How do I construct the projector onto a subspace and verify `P^2=P=P^\dagger`?
- 0.5 [BSc] How do I test whether a matrix is unitary or an isometry, and how do I interpret it as a change of basis?
- 0.6 [BSc] How do I embed a local operator into a composite Hilbert space using a declared tensor-factor order?
- 0.7 [BSc] What global conventions are used for `hbar`, logarithm bases, tensor ordering, computational basis ordering, and numerical tolerances?

Without these, later prompts such as spectral decomposition, partial trace, channels, circuits, and tomography depend on unspoken choices.

### 2. The "One Central Concept" Rule Is Violated Often

Many questions bundle construction, proof, simulation, interpretation, and a second theorem into one prompt. This is the largest structural problem.

High-priority examples to split:

| Item | Problem | Suggested split |
|------|---------|-----------------|
| 2.2 | Builds Pauli matrices and verifies the full algebra. | One item for construction; one for commutators/anticommutators/product rule. |
| 2.8 | Compatibility test plus simultaneous eigenbasis construction. | One item for commutation/compatibility; one for constructing a joint eigenbasis. |
| 3.6 | Bloch-vector map plus positivity plus pure-state boundary. | One item for the Bloch map; one for positivity radius; one for pure boundary. |
| 4.2 | Partial trace plus entanglement-induced mixedness. | One item for partial trace; later item for reduced mixedness of an entangled pure state. |
| 7.1 | Builds `exp(-iHt)` and solves the Schrodinger equation. | One item for time-evolution operator; one for solving state evolution. |
| 7.5 | Heisenberg evolution plus Ehrenfest theorem. | One item for Heisenberg-picture observables; one for Ehrenfest. |
| 7.9 | Finite-time transition probability plus Fermi golden rule. | One item for first-order transition probability; one for dense-spectrum golden-rule limit. |
| 7.10 | Symmetry, conserved quantum number, and block diagonalization. | One item for commutation/conservation; one for symmetry-sector block diagonalization. |
| 9.1 | All standard gates plus arbitrary-axis rotations. | Separate standard gates, phase gates, and rotations. |
| 10.2 | Circuit composition, application, and unitary extraction. | Separate state application from overall unitary construction. |
| 11.1 | Bell, GHZ, and W states are three different ideas. | Separate bipartite Bell states, GHZ states, and W states. |
| 11.7 | Concurrence and negativity are different entanglement monotones. | Separate them. |
| 12.5 | Trace distance and fidelity are distinct state distinguishability notions. | Separate them. |
| 13.6 | Converts among four channel representations. | Use one item per conversion or one carefully bounded pairwise conversion. |
| 13.9 | Channel composition, CP-divisibility, and Markovianity. | Separate composition from CP-divisibility test. |
| 14.3 | Joint, conditional, mutual entropy, and negative conditional entropy. | Separate each quantity; then add an example with negative conditional entropy. |
| 16.1 | SIC construction, positivity/completeness, and application. | Separate construction/verification from measurement application. |
| 20.4 | Classical linear codes plus CSS foundation. | Separate classical code arithmetic from CSS construction. |
| 22.2 | Jordan-Wigner map plus exact solvability of two chains. | Separate the map from diagonalizing a concrete finite chain. |
| 23.5 | Grover amplification plus optimal iteration count. | Separate amplitude amplification dynamics from optimal iteration formula. |
| 23.14 | QAOA construction plus optimization-depth behavior. | Separate p=1 construction from numerical depth comparison. |

### 3. Several Questions Are Not Self-Contained

A self-contained question should specify the input object, the convention, the output, and the success criterion. Several entries use phrases like "a given state", "a simple channel", "arbitrary", "general", "watch", "exhibit", "locate", or "estimate" without enough data.

Examples:

| Item | Ambiguity | Revision direction |
|------|-----------|--------------------|
| 1.2 | Born-rule probabilities for which measurement? | Say "in the computational basis" or supply a PVM `{P_k}`. |
| 1.8 | Bloch angles depend on phase and coordinate conventions. | Specify `rho=(I+r.sigma)/2` and angle convention. |
| 2.8 | Degeneracy complicates "find a simultaneous eigenbasis". | Specify finite Hermitian matrices and handle degenerate eigenspaces. |
| 4.2 | Partial trace depends on tensor ordering and subsystem dimensions. | Include declared dimensions `{d_A,d_B}` and ordering. |
| 5.1 | "Projective measurement" needs projectors and outcome labels. | Given `{P_k,a_k}`, compute probabilities, selective states, and mean. |
| 6.2 | A saturating state may not exist for an arbitrary finite pair in a useful form. | Use a chosen pair, e.g. Pauli `X,Z`, or ask for numerical minimization. |
| 7.8 | Rotating-frame Rabi problem depends on Hamiltonian and rotating-wave approximation. | State the driven two-level Hamiltonian and whether RWA is used. |
| 12.7 | Helstrom measurement needs priors and two specified density matrices. | Give `p1,p2,rho1,rho2`; output projectors onto positive/negative parts. |
| 14.7 | "Channel-capacity estimate" is not a single well-defined computation. | Choose one finite channel and one lower/upper bound to compute. |
| 17.1 | Tomography needs an informationally complete measurement and estimator. | Specify Pauli tomography or SIC tomography and a reconstruction rule. |
| 21.5 | Exact stabilizer rank is hard and dimension-dependent. | Restrict to a named one- or two-qubit magic state, or ask for an upper-bound decomposition. |

### 4. Some Technical Statements Need Correction

These are the most important correctness issues.

| Item | Issue | Suggested correction |
|------|-------|----------------------|
| 7.1 | Uses `e^{-iHt}` without declaring `hbar=1`; later 7.7 includes `hbar`. | Add a global `hbar=1` convention or write `e^{-iHt/hbar}` consistently. |
| 7.9 | Fermi's golden rule is a continuum/dense-spectrum limit, not literally finite-dimensional dynamics. | Phrase it as a finite dense-band approximation and compare to the continuum rate in the dense-band limit. |
| 7.12 | `T >> 1/Delta^2` is incomplete and dimensionally schematic. | Use the standard finite-dimensional adiabatic condition involving `max_s ||dH/ds|| / (T Delta_min^2)` or a Landau-Zener example. |
| 12.8a | "Verify Gleason's theorem" is not a finite computational exercise in full generality. | Ask students to reconstruct a density matrix from a finite informationally complete frame function, then state the theorem separately. |
| 12.8b | Purely conceptual "Why does the qubit fail Gleason" violates the compute/verify style. | Reframe as constructing a dispersion-free frame function on qubit projectors and showing why POVMs remove the exception. |
| 13.10 | Lindblad `T1/T2` dynamics is usually not BSc; it also needs generator conventions. | Retag as MSc or split into a BSc matrix-exponential semigroup example and an MSc Lindblad item. |
| 15.2 | "Exactly recompiling a random unitary" with a universal gate set is generally false for finite gate sets. | Use "approximate to error epsilon" for finite universal sets, or restrict exact synthesis to a continuous gate set. |
| 18.5 | "MUBs for a general qudit via the Fourier basis" is false if it means a maximal set in all dimensions; Fourier gives one mutually unbiased partner to the computational basis. | Say "build the computational and Fourier MUB pair in dimension d" or "build a complete set in prime-power dimension". |
| 18.6 | "Construct a SIC-POVM in general dimension d" is not known constructively for all d. | Restrict to `d=2` or `d=3`, or say "given candidate vectors in dimension d, verify the SIC conditions". |
| 20.10b | "What does the threshold theorem guarantee" is conceptual, not a WL operation. | Replace with a finite concatenated-code threshold toy model and put the theorem statement in exposition. |
| 21.2 | The even-dimension phase-space issue is subtle; Gibbons-Hoffman-Wootters is a finite-field construction for prime powers, not just "the even-d workaround". | Specify qubit `d=2` and the exact phase-point operators being used. |
| 23.11 | HHL does not output the classical vector `x`; it prepares a quantum state proportional to `A^{-1}|b>`. | Say "prepare `|x>` proportional to `A^{-1}|b>` and compare amplitudes to the normalized classical solution". |
| 23.14 | "Adding layers improve it" is not guaranteed for a practical local optimizer, although the variational optimum is nonincreasing for energy minimization with nested ansatz families. | Ask to compare the optimized objective for `p=1,2,3` on a fixed graph and optimizer protocol. |

### 5. BSc/MSc Tags Need Rebalancing

The BSc track is strong but uneven. A BSc student can reasonably do pure states, observables, projective measurement, simple unitary dynamics, spin-1/2, basic gates, Bell states, teleportation, Deutsch/Deutsch-Jozsa/Grover basics, and no-cloning.

The following should probably be MSc unless the course is explicitly quantum-information-heavy and mathematically advanced:

- 12.1 Werner-state entanglement threshold.
- 13.3 Choi matrix and complete positivity.
- 13.10 Lindblad `T1/T2` relaxation and dephasing.
- 18.5 MUBs for general qudits.
- 9.4 arbitrary single-qubit Euler decomposition, if the BSc course is physics rather than quantum computing.

Some MSc items could remain MSc but need stronger prerequisites inserted:

- 11.9 Nielsen theorem requires majorization to be introduced first.
- 13.15 diamond norm requires trace norm, ancilla optimization, and preferably an SDP formulation.
- 17.3 quantum Fisher information requires classical Fisher information and the symmetric logarithmic derivative.
- 20.5 CSS codes require a stabilizer-code construction item before Steane.

## Coverage Gaps

The list is broad, but the following topics are important for a finite-dimensional BSc/MSc curriculum and are either missing or only implicit.

### Core Finite-Dimensional Quantum Mechanics

- Change of basis and active/passive transformations.
- Projectors, subspaces, and degeneracy as first-class objects.
- Normal operators and why Hermitian/unitary matrices have orthonormal eigenbases.
- Degenerate projective measurements and incomplete measurements.
- Sequential measurements and measurement disturbance.
- Wigner theorem in finite dimension, at least as a statement about unitary/antiunitary symmetry transformations.
- Time-reversal as an antiunitary finite-dimensional symmetry, if foundations/symmetry is a goal.

### Composite Systems

- Product states versus separable mixed states versus entangled states before advanced entanglement measures.
- Local operators and identity padding.
- Swap/permutation operators for tensor factors.
- Conditional states after local measurement.

### Quantum Information Theory

- Entropy inequalities: subadditivity, Araki-Lieb, and strong subadditivity.
- Fannes/Audenaert continuity bound or at least continuity of entropy in trace distance.
- Purified distance or Bures distance if fidelity is used heavily.
- Majorization and Schur concavity before Nielsen's theorem.
- Resource-theory definitions: free states, free operations, monotones, at least once before coherence/magic/entanglement resources.

### Quantum Computation

- Reversible classical computation and oracle construction.
- Phase kickback as its own BSc/MSc bridge item.
- Resource scaling: gate count, query count, qubit count, and precision.
- Controlled arithmetic for phase estimation/Shor/HHL.
- Solovay-Kitaev should be framed as an approximation theorem, not as a practical exact compiler.

### Error Correction and Fault Tolerance

- Stabilizer code from an abelian Pauli subgroup, code space, syndrome, normalizer, logical Pauli group.
- Syndrome table construction before Shor/Steane/five-qubit examples.
- Fault-tolerant syndrome extraction as a separate item from code construction.
- Eastin-Knill/no-universal-transversal-gates as an MSc concept if transversal gates are included.

### Finite Many-Body Physics

- Correlation functions and connected correlations.
- Translation symmetry and momentum sectors.
- Exact diagonalization using symmetry sectors.
- Finite-size scaling and limits of diagnosing "phase transitions" in finite systems.
- Thermalization/ETH in finite spectra, if quench and many-body dynamics are included.

## Suggested Global Style Rule

Use this template for every question:

> Given `[explicit finite-dimensional input]` under `[declared convention]`, how do I compute or verify `[one output or property]`?

Examples:

- Weak: "How do I read the Born-rule probabilities of a state?"
- Strong: "Given a normalized state `psi=sum_k c_k |k>` in the computational basis, how do I compute the probability vector `p_k=|c_k|^2` and verify `sum_k p_k=1`?"

- Weak: "How do I take a partial trace?"
- Strong: "Given `rho_AB` on `C^dA tensor C^dB` with basis order `|a>|b>`, how do I compute `rho_A=Tr_B rho_AB` and verify `Tr rho_A=1`?"

- Weak: "How do I confirm a gate set is universal?"
- Strong: "Given the finite gate set `{H,T,CNOT}` and a target one-qubit unitary `U`, how do I synthesize a circuit whose operator-norm error is at most `epsilon`?"

## Concrete Revision Examples

These examples show how to rewrite existing items without changing the overall curriculum.

### Part 1

- Replace 1.2 with: "Given `|psi>=sum_k c_k |k>` in the computational basis, how do I compute `p_k=|c_k|^2` and verify normalization of the probability distribution?"
- Split 1.8 into:
  - "Given a pure qubit `|psi>`, how do I compute its Bloch vector `r=<sigma>`?"
  - "Given a nonzero Bloch vector `r`, how do I compute polar and azimuthal angles under a declared convention?"
- Clarify 1.9a/1.9b by saying whether "prepare" means "construct the state vector" or "give a circuit that prepares it."

### Part 2

- Add an item before 2.4: "Given a Hermitian operator with degenerate eigenvalues, how do I construct its spectral projectors?"
- Replace 2.8 with two items:
  - "Given two Hermitian matrices `A,B`, how do I test whether `[A,B]=0`?"
  - "Given commuting Hermitian matrices `A,B`, how do I construct a simultaneous eigenbasis, including degenerate eigenspaces?"
- Move the Hilbert-Schmidt inner product earlier if it is used later for tomography, Pauli expansions, Choi matrices, or operator bases.

### Parts 4 and 5

- Add: "Given subsystem dimensions and tensor ordering, how do I reshape a vector or matrix to expose subsystem indices?"
- Add: "Given a local operator `A` on subsystem `A`, how do I form `A tensor I_B` and compute its expectation value on `rho_AB`?"
- Replace 5.1 with: "Given a PVM `{P_k}` and outcome labels `{a_k}`, how do I compute `p_k=Tr(P_k rho)`, the conditional post-measurement states, and the mean `sum_k a_k p_k`?"
- Add a separate item for degenerate measurements and Luders updates.

### Part 7

- State globally whether `hbar=1`. If not, revise all dynamics formulas.
- Replace 7.9 with:
  - "Given `H(t)=H0+lambda V cos(omega t)`, how do I compute the first-order transition probability between two nondegenerate eigenstates of `H0`?"
  - "Given a finite dense band of final states, how do I numerically observe the approach to Fermi's golden-rule scaling?"
- Replace 7.12 with a concrete two-level avoided crossing or Landau-Zener model.

### Parts 13 and 14

- Retag Choi/complete positivity and Lindblad equations as MSc unless the course is explicitly advanced undergraduate quantum information.
- Split 13.6 into pairwise conversions:
  - Kraus to Choi.
  - Choi to Kraus by diagonalization.
  - Kraus to superoperator.
  - Stinespring isometry to Kraus by environment projection.
- Add strong subadditivity before data processing, or explicitly state data processing as a theorem being numerically checked in examples.

### Parts 18 and 21

- Replace 18.5 with: "How do I build the computational basis and Fourier basis in dimension `d` and verify they are mutually unbiased?"
- Add a separate MSc item: "For prime-power dimension `d`, how do I construct a complete set of `d+1` mutually unbiased bases using finite-field structure?"
- Replace 18.6 with: "In `d=2` or `d=3`, how do I construct a known SIC-POVM and verify positivity, completeness, and equiangularity?"
- For general `d`, use: "Given a candidate set of `d^2` normalized vectors, how do I verify the SIC-POVM equations?"

### Parts 19 and 20

- Add before 20.1: "Given commuting independent Pauli generators, how do I compute the stabilizer code space dimension?"
- Add before 20.2: "Given a stabilizer code, how do I compute its syndrome table for a specified finite error set?"
- Replace 20.10b with a computational toy model: "Given a recursive logical error model `p_{l+1}=A p_l^2`, how do I compute the threshold and logical error after `L` concatenation levels?"

### Part 23

- Add: "How do I implement phase kickback for a Boolean oracle?"
- Add: "How do I convert a classical reversible function into a unitary permutation matrix?"
- Revise HHL to: "Given a small Hermitian matrix `A` with known condition number and a normalized `|b>`, how do I prepare a state proportional to `A^{-1}|b>` and compare it with the normalized classical solution?"
- Revise QAOA to use a fixed graph, fixed depth range, and fixed optimizer protocol.

## Recommended Revision Plan

1. Add Part 0 for conventions and finite-dimensional linear algebra.
2. Enforce the question template: explicit input, convention, one computation, one success criterion.
3. Split every item that contains "and" joining two independent tasks, especially in Parts 7, 13, 14, 20, and 23.
4. Correct the technically over-strong items: 15.2, 18.5, 18.6, 20.10b, 23.11.
5. Rebalance the BSc/MSc tags after adding prerequisites.
6. Add the missing bridge questions for oracle construction, stabilizer codes, entropy inequalities, and local operators.
7. Regenerate the coverage table after splitting, because the item count will likely increase from 207 to roughly 240-270.

## Bottom Line

As a topic map, the file is strong. As a self-contained question curriculum, it needs another pass. The biggest improvements are mechanical and achievable: add explicit conventions, split overloaded prompts, fix a few technical overclaims, and turn vague "show/why/exhibit/watch" entries into finite, checkable computations.
