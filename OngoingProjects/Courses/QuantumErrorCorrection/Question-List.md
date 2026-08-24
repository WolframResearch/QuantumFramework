# Quantum Error Correction: A Pedagogical Question Collection

From classical redundancy and the first quantum codes to fault-tolerant architectures and the
research frontier, pinned to 2026-08-23.

This collection contains 188 atomic questions in dependency order. It is designed to be answered
with a mixture of exact derivation, finite computation, circuit analysis, decoder implementation,
stochastic benchmarking, and primary-paper reconstruction. The binding authoring and verification
rules are in [`PIPELINE.md`](PIPELINE.md).

## Tags and evidence modes

- `[Bridge]`: prerequisite or transition material.
- `[Core]`: durable QEC knowledge.
- `[Advanced]`: graduate or research-facing structure.
- `[Frontier-2026]`: a dated result that must be re-audited when the frontier changes.
- `[Q0]`: exact finite-field or Pauli algebra.
- `[Q1]`: exact finite-dimensional channel or recovery calculation.
- `[Q2]`: circuit fault propagation and bounded fault enumeration.
- `[Q3]`: decoder construction on a fixed syndrome model.
- `[Q4]`: stochastic simulation, threshold estimation, or resource analysis.
- `[Q5]`: proof or conceptual argument.
- `[Q6]`: reconstruction of a dated primary result.

## Part 0. Algebraic bridge and conventions (7 questions)

- 0.1 `[Bridge] [Q0]` How do I row-reduce a binary matrix over $\mathrm{GF}(2)$, compute its rank and null space, and verify the result without accidentally using real arithmetic?
- 0.2 `[Bridge] [Q0]` How do I represent an $n$-qubit Pauli operator modulo global phase by a binary vector $(\mathbf{x}\mid\mathbf{z})$?
- 0.3 `[Bridge] [Q0]` How does the binary symplectic product decide whether two Pauli strings commute or anticommute?
- 0.4 `[Bridge] [Q0]` How do Clifford gates transform binary Pauli vectors, and how do I verify that the transformation preserves the symplectic product?
- 0.5 `[Bridge] [Q1]` How do I convert a small quantum channel among Kraus, superoperator, and Choi representations and verify complete positivity and trace preservation?
- 0.6 `[Bridge] [Q1]` How do I represent a quantum measurement together with its classical outcome register, rather than treating the outcome as an unexplained collapse label?
- 0.7 `[Bridge] [Q5]` What is the operational difference among a physical qubit, data qubit, ancilla qubit, gauge qubit, logical qubit, code block, stabilizer, syndrome, recovery, decoder, and Pauli frame?

## Part 1. Classical codes as the inference model (8 questions)

- 1.1 `[Bridge] [Q1]` How do I encode one bit in the three-bit repetition code, send it through a binary symmetric channel, decode by majority vote, and compute the exact failure probability?
- 1.2 `[Bridge] [Q0]` Given a binary generator matrix $G$, how do I construct a compatible parity-check matrix $H$ and verify $H G^{\mathsf{T}} = 0$ over $\mathrm{GF}(2)$?
- 1.3 `[Bridge] [Q0]` How does a syndrome identify an error coset rather than a unique physical error?
- 1.4 `[Bridge] [Q0]` How do I compute the minimum distance of a small linear code by exhaustive enumeration and by dependence among columns of $H$?
- 1.5 `[Bridge] [Q3]` How does maximum-likelihood decoding change when different bit locations have different error probabilities?
- 1.6 `[Bridge] [Q3]` Why are known erasure locations easier to decode than unknown flips, and how do I solve a fixed erasure pattern by linear algebra?
- 1.7 `[Advanced] [Q3]` How does belief propagation decode a classical code on a cycle-free factor graph, and which exact marginal provides the verification anchor?
- 1.8 `[Bridge] [Q4]` For the concatenated repetition-code recursion $p_{\ell+1} = 3p_{\ell}^{2} - 2p_{\ell}^{3}$, how do I find the unstable fixed point and distinguish a threshold from finite-level improvement?

## Part 2. Why quantum error correction is possible (9 questions)

- 2.1 `[Core] [Q5]` Why does the no-cloning theorem forbid protecting an unknown state by making independent copies but not by encoding it into an entangled subspace?
- 2.2 `[Core] [Q1]` How do I encode $\alpha\lvert 0\rangle + \beta\lvert 1\rangle$ in the three-qubit bit-flip code, extract a syndrome without learning $\alpha$ or $\beta$, and recover from one $X$ error?
- 2.3 `[Core] [Q1]` How does conjugation by Hadamard turn the bit-flip code into a phase-flip code that corrects one $Z$ error?
- 2.4 `[Core] [Q5]` Why is the syndrome distribution independent of the logical amplitudes for a correctable error set?
- 2.5 `[Core] [Q1]` How does expanding a one-qubit operator in the Pauli basis reduce correction of a continuum of errors to correction of a finite operator basis?
- 2.6 `[Core] [Q1]` For a coherent $X$-rotation on one physical qubit of the three-qubit bit-flip code, what logical channel remains after ideal syndrome measurement and recovery, both conditionally and after the syndrome outcome is forgotten?
- 2.7 `[Core] [Q1]` How does the nine-qubit Shor code combine phase-flip and bit-flip protection to correct an arbitrary single-qubit error?
- 2.8 `[Core] [Q5]` What operationally distinguishes error detection with postselection, active recovery, and a Pauli-frame update?
- 2.9 `[Core] [Q1]` How do I compose encoding, noise, syndrome extraction, recovery, and decoding into one logical channel and exhibit its first uncorrectable error pattern?

## Part 3. Noise models and performance observables (11 questions)

- 3.1 `[Core] [Q1]` How do I write a one-qubit Pauli channel, verify that its probabilities define a completely positive trace-preserving map, and identify its action on the Bloch vector?
- 3.2 `[Core] [Q1]` How does pure dephasing suppress off-diagonal coherence while preserving populations, and which repetition basis is adapted to it?
- 3.3 `[Core] [Q1]` Why is amplitude damping not a Pauli channel, and how does that difference appear in its Kraus operators and fixed state?
- 3.4 `[Core] [Q1]` How can a coherent over-rotation and a stochastic Pauli channel have the same average infidelity but produce different logical scaling under repeated use?
- 3.5 `[Advanced] [Q1]` What does Pauli twirling preserve, which phase-sensitive or off-diagonal process information does it discard, and when can the discarded structure change a QEC prediction?
- 3.6 `[Core] [Q5]` What is the operational difference among leakage, loss, and erasure, and exactly which location information is available to the decoder in each case?
- 3.7 `[Core] [Q1]` How do I parameterize biased Pauli noise and quantify the bias without hiding the total physical error probability?
- 3.8 `[Advanced] [Q4]` How can a rare correlated burst produce a logical-error floor that an independent-error model cannot predict?
- 3.9 `[Core] [Q2]` What locations and conditional dependencies must be specified before a noise model deserves to be called circuit-level?
- 3.10 `[Core] [Q4]` Given repeated logical-memory trials, how do I estimate logical error per cycle with an uncertainty interval and distinguish it from physical gate infidelity?
- 3.11 `[Advanced] [Q1]` How do average gate infidelity, entanglement fidelity, unitarity, diamond distance, and logical failure probability differ, and which of them can support a fault-tolerance claim?

## Part 4. Error-correction conditions and recovery maps (11 questions)

- 4.1 `[Core] [Q5]` How do I derive the Knill-Laflamme condition $P E_{a}^{\dagger} E_{b} P = c_{ab}P$ from the requirement that the environment learn nothing about the encoded state?
- 4.2 `[Core] [Q1]` How do I check the Knill-Laflamme condition explicitly for the three-qubit bit-flip code and the error set $\{I,X_{1},X_{2},X_{3}\}$?
- 4.3 `[Core] [Q5]` Why does correcting an operator basis imply correction of every channel whose Kraus operators lie in its linear span?
- 4.4 `[Core] [Q1]` How does degeneracy appear when two distinct physical errors act identically on the code space?
- 4.5 `[Advanced] [Q5]` How is correctability of a known erasure region related to decoupling of that region from a reference system entangled with the logical state?
- 4.6 `[Core] [Q1]` Given orthogonal error subspaces, how do I construct an explicit syndrome-and-recovery channel and verify it on a logical state entangled with a reference?
- 4.7 `[Core] [Q1]` How do I compute entanglement fidelity of an encoded channel without testing only the logical computational-basis states?
- 4.8 `[Advanced] [Q1]` How do approximate Knill-Laflamme conditions quantify residual logical information in the environment when exact correction is impossible?
- 4.9 `[Advanced] [Q5]` How does operator-algebra or subsystem QEC generalize protection of a subspace to protection of a logical algebra while allowing gauge information to change?
- 4.10 `[Advanced] [Q1]` For a fixed encoding and noise channel, how do the transpose or Petz-type recovery and an SDP-optimized recovery compare in entanglement fidelity, and what exact small instance anchors the optimization?
- 4.11 `[Advanced] [Q5]` How does regularized coherent information characterize quantum channel capacity, how can a one-shot coherent-information calculation supply an achievable rate, and why are capacity, a memory threshold, and experimental break-even different claims?

## Part 5. Stabilizer and binary symplectic methods (12 questions)

- 5.1 `[Core] [Q0]` Given $r$ independent commuting Pauli generators on $n$ qubits, how do I derive the stabilized subspace dimension $2^{n-r}$?
- 5.2 `[Core] [Q0]` How do I test whether a proposed generator list is independent, mutually commuting, and free of the forbidden element $-I$?
- 5.3 `[Core] [Q0]` How do I row-reduce a stabilizer tableau without changing the stabilized code space?
- 5.4 `[Core] [Q0]` Given a stabilizer and a finite Pauli error set, how do I construct the complete syndrome table and identify syndrome degeneracies?
- 5.5 `[Core] [Q0]` How do I compute the normalizer and choose logical $\overline{X}_{i},\overline{Z}_{i}$ representatives with the required commutation relations?
- 5.6 `[Core] [Q0]` How does the binary condition $H_{X}H_{Z}^{\mathsf{T}} = 0$ encode commutation for a CSS stabilizer?
- 5.7 `[Core] [Q0]` How do I update a stabilizer tableau under $H$, $S$, and CNOT gates and verify the result by direct conjugation?
- 5.8 `[Core] [Q2]` How does an ancilla circuit measure a Pauli stabilizer without measuring any logical observable?
- 5.9 `[Core] [Q1]` How do I construct the projector onto a stabilizer code as an average over the stabilizer group and verify its rank?
- 5.10 `[Advanced] [Q0]` How do I synthesize an encoding Clifford circuit from a stabilizer tableau and verify that it maps standard generators to the target stabilizer?
- 5.11 `[Advanced] [Q1]` How does stabilizer syndrome measurement process a non-Pauli error without implying that the original physical noise was stochastic Pauli noise?
- 5.12 `[Advanced] [Q0]` How do generalized Pauli operators and the trace-symplectic form extend stabilizer codes from qubits to prime-power qudits?

## Part 6. Code parameters, distance, and bounds (9 questions)

- 6.1 `[Core] [Q0]` How do I determine $n$, $k$, and $d$ for a stabilizer code from its generators and normalizer?
- 6.2 `[Core] [Q0]` How do I compute distance by searching for the minimum-weight Pauli in the normalizer but outside the stabilizer?
- 6.3 `[Advanced] [Q0]` When a CSS code has unequal $X$- and $Z$-distances, how do I compute $d_{X}$ and $d_{Z}$ separately and match them to biased noise?
- 6.4 `[Core] [Q5]` Why can degeneracy allow different physical errors to share a syndrome without reducing the code's logical distance?
- 6.5 `[Core] [Q5]` How is the quantum Hamming bound derived for a nondegenerate $[[n,k,2t+1]]$ code, and why is degeneracy a necessary caveat?
- 6.6 `[Advanced] [Q5]` How does the quantum Singleton bound $n-k \geq 2(d-1)$ follow from erasure correctability and no-cloning?
- 6.7 `[Core] [Q5]` Why can a code of distance $d$ correct $d-1$ erasures but only $\lfloor(d-1)/2\rfloor$ errors at unknown locations?
- 6.8 `[Advanced] [Q2]` How can a syndrome-extraction circuit have circuit-level distance smaller than the distance of its abstract code?
- 6.9 `[Advanced] [Q5]` How do the cleaning lemma and locality-rate-distance trade-offs constrain two-dimensional geometrically local codes?

## Part 7. Canonical block, CSS, and subsystem codes (11 questions)

- 7.1 `[Core] [Q0]` Given classical codes $C_2\subseteq C_1$, how do I construct the associated CSS stabilizer and derive its encoded dimension?
- 7.2 `[Core] [Q0]` How do I construct the Steane $[[7,1,3]]$ code from the classical Hamming code and identify its logical Pauli operators?
- 7.3 `[Core] [Q0]` How do I build the full single-qubit Pauli syndrome table for the Steane code and verify correction modulo its stabilizer?
- 7.4 `[Core] [Q0]` How do I express the Shor code as a stabilizer code and recover its bit- and phase-syndrome structure from the generators?
- 7.5 `[Core] [Q0]` How do I verify that the five-qubit $[[5,1,3]]$ code is non-CSS, corrects every single-qubit Pauli error, and saturates the quantum Hamming bound?
- 7.6 `[Core] [Q4]` Under one fixed circuit-level noise model, how do the Shor, Steane, and five-qubit memories differ in logical error and syndrome-extraction overhead?
- 7.7 `[Advanced] [Q0]` How do the gauge generators of the Bacon-Shor code produce a lower-weight measurement scheme while preserving the protected logical subsystem?
- 7.8 `[Advanced] [Q2]` How does gauge fixing change the stabilizer and available logical operations without measuring the logical state?
- 7.9 `[Advanced] [Q1]` How does a four-qubit channel-adapted code satisfy the correction conditions for leading-order amplitude-damping errors?
- 7.10 `[Advanced] [Q0]` How do I search exhaustively for a small stabilizer code with requested $[[n,k,d]]$ parameters while quotienting out redundant generator descriptions?
- 7.11 `[Advanced] [Q5]` Given a stated noise bias, connectivity graph, measurement budget, and target logical operation, what evidence is sufficient to justify choosing one block code over another?

## Part 8. Fault-tolerant syndrome extraction, preparation, and measurement (14 questions)

- 8.1 `[Core] [Q2]` How can one fault on a reused syndrome ancilla propagate into an uncorrectable multi-qubit data error?
- 8.2 `[Core] [Q2]` How does the ordering of controlled gates create or suppress hook errors in a weight-four stabilizer measurement?
- 8.3 `[Core] [Q3]` Why do noisy stabilizer measurements in standard surface-code syndrome extraction require repeated rounds and decoding in space-time rather than majority voting each stabilizer independently?
- 8.4 `[Advanced] [Q2]` How does verification of a Shor cat ancilla prevent one preparation fault from becoming a high-weight data error?
- 8.5 `[Advanced] [Q2]` How does Steane error correction extract an encoded syndrome transversally, and which ancilla faults must be rejected?
- 8.6 `[Advanced] [Q2]` How does teleportation-based Knill error correction combine syndrome extraction, recovery, and movement of the logical state?
- 8.7 `[Advanced] [Q2]` How can a flag qubit make a dangerous correlated fault distinguishable without using a full verified cat state?
- 8.8 `[Core] [Q2]` How do I exhaustively verify the circuit-level distance of a small extraction circuit by enumerating all faults up to a specified order?
- 8.9 `[Advanced] [Q2]` How can leakage persist across QEC cycles, and how does a leakage-removal unit alter both the fault model and the schedule?
- 8.10 `[Advanced] [Q4]` Under identical hardware and noise assumptions, how do I compare two syndrome-extraction schedules using logical failure, depth, ancilla count, and decoder information?
- 8.11 `[Core] [Q2]` How do I prepare a logical $\lvert 0_{L}\rangle$ fault-tolerantly and verify that one preparation fault cannot become an undetected logical error?
- 8.12 `[Core] [Q2]` How do I measure a logical Pauli fault-tolerantly when physical readout errors require redundancy or repeated syndrome information?
- 8.13 `[Advanced] [Q5]` What redundant syndrome structure and confinement property allow single-shot QEC to tolerate measurement errors without repeating a number of extraction rounds that grows with code distance?
- 8.14 `[Advanced] [Q2]` How do rectangles, extended rectangles, malignant fault sets, and level reduction turn gadget-level fault enumeration into a recursive threshold proof?

## Part 9. Fault-tolerant logical operations (12 questions)

- 9.1 `[Core] [Q2]` Why does a transversal gate limit the spread of one physical fault within each code block without automatically providing a universal gate set?
- 9.2 `[Core] [Q0]` Which Clifford gates act transversally on the Steane code, and how do I verify their logical action on stabilizers and logical Paulis?
- 9.3 `[Advanced] [Q5]` What does the Eastin-Knill theorem rule out, which assumptions carry the theorem, and which fault-tolerant constructions evade those assumptions without contradicting it?
- 9.4 `[Core] [Q2]` How does gate teleportation implement a logical operation using entanglement, measurement, and a classically controlled frame update?
- 9.5 `[Core] [Q2]` How does magic-state injection turn a verified non-stabilizer resource state into a logical $T$ gate?
- 9.6 `[Advanced] [Q1]` In one fixed magic-state distillation protocol, how do input error, acceptance probability, and output error transform in one round?
- 9.7 `[Advanced] [Q5]` How does the Clifford hierarchy organize gates that propagate Pauli corrections into progressively more complex corrections?
- 9.8 `[Advanced] [Q2]` How can code switching or gauge fixing expose a logical gate unavailable transversally in the original code?
- 9.9 `[Core] [Q2]` How does lattice surgery implement a logical CNOT through joint logical parity measurements and frame updates?
- 9.10 `[Advanced] [Q4]` How do I compare transversal, teleportation, code-switching, and lattice-surgery implementations of one logical gate under a common error and resource model?
- 9.11 `[Advanced] [Q5]` How do locality-preserving gate classifications strengthen Eastin-Knill for topological stabilizer codes, and why are geometrically local constant-depth logical gates in two dimensions restricted to Clifford operations?
- 9.12 `[Advanced] [Q4]` Under one circuit-level model and target output error, how do conventional magic-state distillation and magic-state cultivation differ in acceptance, space-time volume, and failure mechanisms?

## Part 10. Thresholds, break-even, and resource overhead (9 questions)

- 10.1 `[Core] [Q4]` Given a recursive logical-error model, how do I identify its threshold and compute the suppression after a fixed number of concatenation levels?
- 10.2 `[Core] [Q5]` What is the difference between an asymptotic threshold and a finite-size pseudo-threshold?
- 10.3 `[Advanced] [Q5]` Which locality, independence, fresh-ancilla, and decoder assumptions enter a fault-tolerance threshold theorem?
- 10.4 `[Core] [Q4]` Why do code-capacity, phenomenological, and circuit-level noise models produce different threshold values for the same code family?
- 10.5 `[Advanced] [Q4]` How do I estimate a threshold and critical scaling exponent from finite-distance crossing data without mistaking a finite-size drift for convergence?
- 10.6 `[Core] [Q4]` How do I compute the error-suppression factor $\Lambda$ between code distances and state precisely what experiment duration and decoder it describes?
- 10.7 `[Core] [Q4]` How do I test whether an encoded memory is beyond break-even relative to the best relevant unencoded memory rather than an arbitrarily chosen physical average?
- 10.8 `[Advanced] [Q4]` How do splitting, importance-sampling, or related rare-event methods estimate logical failure probabilities below the reach of direct Monte Carlo without biasing the result?
- 10.9 `[Advanced] [Q4]` Given an algorithmic logical-error budget, how do I estimate data-block, syndrome, routing, and magic-state-factory resources without hiding time overhead?

## Part 11. Topological, surface, and colour codes (14 questions)

- 11.1 `[Core] [Q0]` How do I construct toric-code star and plaquette stabilizers on a periodic square lattice and verify all commutators?
- 11.2 `[Core] [Q0]` How do Pauli error strings create syndrome defects only at their endpoints, and how does this support an anyon interpretation?
- 11.3 `[Core] [Q0]` How do noncontractible loop operators become logical Paulis, and how does their shortest weight determine toric-code distance?
- 11.4 `[Advanced] [Q0]` How does the chain-complex relation $\partial_{1}\partial_{2} = 0$ encode stabilizer commutation and identify logical operators with homology classes?
- 11.5 `[Core] [Q0]` How do rough and smooth boundaries turn the toric code into a planar surface-code patch with one logical qubit?
- 11.6 `[Core] [Q0]` How do I lay out a rotated surface-code patch of distance $d$, count its data and check qubits, and verify its logical distance?
- 11.7 `[Core] [Q2]` How do I schedule one round of surface-code syndrome extraction and identify the hook-error orientation required to preserve distance?
- 11.8 `[Core] [Q3]` How do repeated surface-code measurements turn data and measurement faults into detection events on a three-dimensional space-time graph?
- 11.9 `[Advanced] [Q2]` How does code deformation move a boundary or defect while preserving the encoded logical information?
- 11.10 `[Advanced] [Q2]` How do merge and split operations realize logical parity measurements in lattice surgery?
- 11.11 `[Advanced] [Q0]` How do I construct a small triangular colour code, identify its logical operators, and verify its transversal Clifford action?
- 11.12 `[Advanced] [Q4]` Under biased Pauli noise, how does an XZZX surface-code layout change the decoding geometry and logical suppression relative to the standard layout?
- 11.13 `[Advanced] [Q2]` How can a periodic schedule of low-weight Pauli measurements dynamically create, move, and protect logical qubits in a Floquet code even when no fixed stabilizer code describes the whole cycle?
- 11.14 `[Advanced] [Q2]` How does a three-dimensional topological cluster state convert one spatial direction into syndrome time, and how are logical operations represented by measurements and boundary conditions?

## Part 12. Decoding and logical-error estimation (14 questions)

- 12.1 `[Core] [Q5]` Why is quantum decoding inference over logical equivalence classes rather than identification of the exact physical error?
- 12.2 `[Core] [Q3]` How do I build and verify a complete lookup decoder for a small stabilizer code under a specified Pauli channel?
- 12.3 `[Core] [Q3]` How do I map surface-code detection events to a minimum-weight perfect-matching problem and recover the logical class?
- 12.4 `[Core] [Q3]` How are matching-edge weights derived from an explicit noise model rather than assigned by geometric distance alone?
- 12.5 `[Advanced] [Q3]` How can degeneracy-aware decoding prefer a more probable error class even when its most likely representative is not the minimum-weight error?
- 12.6 `[Advanced] [Q3]` How does union-find decoding grow and merge clusters, and on which erasure or Pauli-noise benchmark can I verify it?
- 12.7 `[Advanced] [Q3]` How do I formulate belief propagation for a CSS quantum LDPC code using its factor graph?
- 12.8 `[Advanced] [Q3]` Why can short cycles and quantum degeneracy prevent ordinary belief propagation from converging to a useful correction?
- 12.9 `[Advanced] [Q3]` How does ordered-statistics post-processing repair a failed belief-propagation result, and what runtime cost does it add?
- 12.10 `[Core] [Q3]` How does a decoder use known erasure locations, and what improvement is lost if those locations are discarded?
- 12.11 `[Advanced] [Q3]` How can analogue measurement confidence and correlated fault probabilities be incorporated without double-counting information?
- 12.12 `[Advanced] [Q4]` How do I train and test a learned decoder without leaking test labels, device identity, or privileged noise-model parameters into the evaluation?
- 12.13 `[Core] [Q4]` How do I compare decoder accuracy, uncertainty, latency, memory use, and scaling on exactly the same event stream?
- 12.14 `[Advanced] [Q3]` How do logical-coset probabilities become partition functions or tensor-network contractions, and how can an exact or controlled approximate maximum-likelihood decoder be verified on a small surface code?

## Part 13. Quantum low-density parity-check codes (13 questions)

- 13.1 `[Advanced] [Q5]` What bounded check weight and bounded qubit degree make a stabilizer family quantum LDPC, and which desirable property is not guaranteed by sparsity alone?
- 13.2 `[Advanced] [Q0]` How do I construct a hypergraph-product CSS code from two classical parity-check matrices and verify stabilizer commutation?
- 13.3 `[Advanced] [Q0]` How do kernels and cokernels of the classical matrices determine the encoded dimension of a hypergraph-product code?
- 13.4 `[Advanced] [Q5]` How does the product of chain complexes explain the checks, logical operators, and distance bounds of a product code?
- 13.5 `[Advanced] [Q5]` What does an asymptotically good quantum LDPC family mean in terms of rate and relative distance, and why was establishing one a qualitative breakthrough?
- 13.6 `[Advanced] [Q5]` How do lifted and balanced product constructions evade the distance limitations of earlier product codes at the level of their underlying complexes?
- 13.7 `[Advanced] [Q5]` Which expansion and local-testability ideas enter quantum Tanner codes, and which of their asymptotic guarantees survive at practical finite size?
- 13.8 `[Advanced] [Q0]` How do I construct a small bivariate bicycle code from two binary circulant-polynomial operators and verify $H_{X}H_{Z}^{\mathsf{T}} = 0$?
- 13.9 `[Advanced] [Q0]` How do I compute the parameters and logical operators of a fixed bivariate bicycle instance rather than inferring them from its nominal family?
- 13.10 `[Advanced] [Q2]` How can nonlocal weight-six checks of a bivariate bicycle code be scheduled on a stated long-range connectivity graph, and which faults limit circuit-level distance?
- 13.11 `[Advanced] [Q4]` Under one circuit-level noise model, how do I compare a finite qLDPC memory with a surface-code memory in logical error, physical qubits, connectivity, extraction depth, and decoder latency?
- 13.12 `[Advanced] [Q3]` How does the small-set-flip decoder use expansion to correct a quantum expander code, and how does its guarantee differ from the empirical performance of BP-OSD?
- 13.13 `[Advanced] [Q2]` How can fault-tolerant logical Pauli measurements and addressable Clifford operations be implemented on a finite-rate qLDPC block without destroying its encoding advantage?

## Part 14. Bosonic quantum error correction (11 questions)

- 14.1 `[Advanced] [Q1]` How do photon loss, gain, dephasing, and small phase-space displacements act differently on an oscillator encoding?
- 14.2 `[Advanced] [Q1]` How do even and odd cat-code states encode a qubit, and which physical error changes photon-number parity?
- 14.3 `[Advanced] [Q4]` How does increasing cat size suppress one logical error while enhancing another, producing a tunable noise bias rather than free protection?
- 14.4 `[Advanced] [Q1]` How does a binomial code satisfy the Knill-Laflamme conditions for a specified finite set of photon-loss and dephasing errors?
- 14.5 `[Advanced] [Q0]` How do the displacement stabilizers and logical translations of an ideal GKP code define a qubit lattice in phase space?
- 14.6 `[Advanced] [Q1]` How does finite squeezing convert ideal GKP correction into a nonzero logical displacement-error probability?
- 14.7 `[Advanced] [Q1]` How does modular-quadrature syndrome extraction infer and correct a small displacement without revealing the encoded logical value?
- 14.8 `[Advanced] [Q1]` How can engineered dissipation or repeated autonomous feedback stabilize a bosonic code space, and what logical noise remains after stabilization?
- 14.9 `[Advanced] [Q4]` How does concatenating biased cat qubits or GKP modes with an outer qubit code change the total logical-error and hardware-overhead budget?
- 14.10 `[Advanced] [Q5]` Under a fixed physical noise model, what operational evidence distinguishes the useful regimes of cat, binomial, and GKP encodings?
- 14.11 `[Advanced] [Q2]` How do I design and fault-enumerate a bosonic logical gate or recovery gadget so that ancilla faults, photon loss, and finite-energy shifts do not spread into an uncorrectable logical error?

## Part 15. Hardware-adapted and autonomous QEC (9 questions)

- 15.1 `[Advanced] [Q2]` How do probabilistic photonic fusion failures and photon loss become located faults in a fusion-based architecture, and what residual Pauli errors must the decoder still model?
- 15.2 `[Advanced] [Q2]` How can an experiment convert loss or leakage into a flagged erasure without introducing a larger unflagged error channel?
- 15.3 `[Advanced] [Q3]` How do I modify a stabilizer decoder to combine erasure flags with ordinary syndrome information?
- 15.4 `[Advanced] [Q2]` What does it mean for a physical entangling gate to preserve noise bias, and how do I test the claim at the logical level?
- 15.5 `[Advanced] [Q4]` Under the same biased circuit-level noise, when does a repetition-based or XZZX encoding outperform a symmetric code?
- 15.6 `[Advanced] [Q2]` How does a leakage-reduction or swap-reset circuit alter the effective syndrome graph and logical failure mechanisms?
- 15.7 `[Advanced] [Q4]` How can syndrome histories calibrate a data-driven noise model, and which correlations remain unidentifiable from syndrome data alone?
- 15.8 `[Advanced] [Q2]` How can mid-circuit erasure detection, adaptive ancilla selection, reset, and qubit reuse be incorporated without breaking fault-tolerance assumptions?
- 15.9 `[Advanced] [Q4]` Given hardware connectivity, native gates, noise bias, measurement speed, and decoder latency, how do I perform a reproducible code-noise-decoder co-design comparison?

## Part 16. Frontier reconstructions, snapshot 2026-08-23 (14 questions)

The primary sources, claim statuses, bounded deliverables, and limitations for these questions are
pinned in [`Frontier-Snapshot-2026.md`](Frontier-Snapshot-2026.md).

- 16.1 `[Frontier-2026] [Q6]` How do I reconstruct the reported distance scaling of a below-threshold superconducting surface-code memory and distinguish the measured suppression factor from an asymptotic threshold claim?
- 16.2 `[Frontier-2026] [Q6]` How do I audit whether a reported real-time surface-code decoder meets the latency requirement imposed by the physical syndrome-cycle time?
- 16.3 `[Frontier-2026] [Q6]` How do I compare the demonstrated memory scaling and logical-operation advantages of colour and surface codes on the same superconducting platform?
- 16.4 `[Frontier-2026] [Q6]` Which experimentally demonstrated ingredients of a neutral-atom logical processor support scalable universal fault-tolerant computation, and which architectural conclusions remain extrapolations?
- 16.5 `[Frontier-2026] [Q6]` What did low-overhead bivariate-bicycle-code experiments directly demonstrate about stabilizer extraction and logical memory, and what did they not yet demonstrate?
- 16.6 `[Frontier-2026] [Q6]` How do theoretical low-overhead qLDPC projections depend on connectivity, circuit-level noise, syndrome extraction, and decoder assumptions?
- 16.7 `[Frontier-2026] [Q6]` How do I reconstruct the logical-error budget of a distance-five repetition code built from biased bosonic cat qubits?
- 16.8 `[Frontier-2026] [Q6]` How was beyond-break-even correction established for GKP-encoded logical qutrits and ququarts, and what was the physical comparison baseline?
- 16.9 `[Frontier-2026] [Q6]` How does a metastable-atom experiment use erasure conversion and adaptive logical circuits, and which Pauli errors remain outside the demonstrated correction capability?
- 16.10 `[Frontier-2026] [Q6]` How do I reproduce a learned decoder's comparison with matching and tensor-network baselines while separating accuracy from training cost and real-time latency?
- 16.11 `[Frontier-2026] [Q6]` How do I reconstruct the reported break-even performance of trapped-ion qLDPC memories and compare their code, connectivity, lifetime baseline, and per-logical error with the earlier superconducting bivariate-bicycle experiment?
- 16.12 `[Frontier-2026] [Q6]` Which ingredients of a measurement-free universal fault-tolerant demonstration were executed coherently on encoded trapped-ion qubits, and which scaling assumptions remain untested?
- 16.13 `[Frontier-2026] [Q6]` How do I reconstruct an experimental logical magic-state-distillation protocol and separate improved output fidelity from postselection, finite-distance, and state-preparation effects?
- 16.14 `[Frontier-2026] [Q6]` How do I build an evidence table that compares current QEC platforms without ranking incomparable claims such as memory lifetime, below-threshold scaling, gate fidelity, logical-qubit count, and projected overhead?

## Coverage and count

| phase | Parts | questions |
|---|---:|---:|
| Algebraic and classical bridge | 0-1 | 15 |
| First quantum codes, noise, and general recovery | 2-4 | 31 |
| Stabilizers, parameters, and canonical codes | 5-7 | 32 |
| Fault tolerance, logical gates, and thresholds | 8-10 | 35 |
| Topological codes and decoding | 11-12 | 28 |
| qLDPC, bosonic, and hardware-adapted QEC | 13-15 | 33 |
| Dated frontier reconstructions | 16 | 14 |
| **Total** | **0-16** | **188** |
