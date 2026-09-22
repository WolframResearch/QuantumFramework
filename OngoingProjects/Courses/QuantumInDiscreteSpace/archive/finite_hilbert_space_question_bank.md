# Quantum Questions in Finite-Dimensional Hilbert Space

This file consolidates questions extracted and transformed from the six local PDFs in this folder. The wording is intentionally not a verbatim reproduction of textbook exercises. Each question is rewritten to target one central mathematical or physical idea in finite-dimensional quantum theory.

PDF page numbers refer to the page number in the PDF file, not necessarily the printed page number in the book.

## Source Abbreviations

- CT1: Cohen-Tannoudji, Diu, Laloe, *Quantum Mechanics, Volume 1*. Main finite-dimensional regions: PDF pages 96-211, 215-383, 388-480, 643-771.
- CT2: Cohen-Tannoudji, Diu, Laloe, *Quantum Mechanics 2*. Main finite-dimensional regions: PDF pages 68-97, 99-193, 196-309, 387-469, 473-556.
- G: Griffiths, *Introduction to Quantum Mechanics*. OCR-derived text. Main finite-dimensional regions: PDF pages 105-142, 172-210, 213-260, 269-278, 352-378, 432-444.
- M: Mermin, *Quantum Computer Science*. Main finite-dimensional regions: most of PDF pages 18-175 and appendices 177-235.
- S: Sakurai and Napolitano, *Modern Quantum Mechanics*. Main finite-dimensional regions: PDF pages 6-70, 72-154, 164-270, 312-390, 459-501.
- Sh: Shankar, *Principles of Quantum Mechanics*. OCR-derived text. Main finite-dimensional regions: PDF pages 15-86, 127-162, 270-284, 290-313, 316-360, 380-432.

## Review Rules Applied

- Keep only questions that can be formulated inside a finite-dimensional complex Hilbert space, or inside a finite-dimensional invariant subspace of a larger problem.
- Rewrite every question so it has one central idea.
- Remove duplicate questions across books unless the variant changes the mathematical content.
- Avoid hidden infinite-dimensional assumptions such as unbounded operators, continuous spectra, or position eigenstates unless the question explicitly asks for a finite-dimensional analogue.
- Prefer precise terms: ray, density operator, Hermitian observable, unitary, projector, POVM element, tensor product, partial trace, degeneracy, irreducible representation.

## A. Finite Hilbert Spaces, Bases, and Dimension

1. In an n-dimensional complex Hilbert space, what algebraic properties distinguish a vector space from a Hilbert space? [CT1, S, Sh, G]
2. What conditions must a set of n vectors satisfy to be an orthonormal basis of an n-dimensional Hilbert space? [CT1, Sh, G]
3. How does the expansion of a vector in one orthonormal basis determine its column-vector representation? [CT1, S, Sh, G]
4. How do the coordinates of a state vector change under a change of orthonormal basis? [CT1, S, Sh, G]
5. What is the distinction between an active unitary transformation of a state and a passive change of basis? [CT1, S, Sh]
6. Why is the dimension of a spin-s Hilbert space equal to 2s+1? [CT1, CT2, G, S, Sh]
7. In a finite-dimensional Hilbert space, how can completeness of an orthonormal basis be expressed as a resolution of the identity? [CT1, S, Sh, G]
8. How does one construct the dual vector associated with a ket, and how does the inner product appear in bra-ket notation? [CT1, S, Sh, M]
9. What is the difference between a subspace, an invariant subspace, and an eigenspace of an operator? [CT1, S, Sh]
10. How does one project an arbitrary vector onto a given subspace using an orthogonal projector? [CT1, S, Sh, G]
11. What information is basis-independent in a finite-dimensional quantum state, and what information is representation-dependent? [CT1, S, Sh, M]
12. How can a finite-dimensional Hilbert space model a two-level atom, a spin-1/2 particle, or a qubit without reference to position wavefunctions? [CT1, CT2, G, M, S, Sh]

## B. States, Rays, Superposition, and Phases

13. Why do normalized vectors that differ only by a nonzero global phase represent the same pure quantum state? [CT1, G, S, M]
14. What is a ray in Hilbert space, and why is it the correct mathematical object for a pure state? [CT1, S, Sh, G]
15. How does the superposition principle differ from classical probabilistic mixing? [CT1, CT2, G, M, S]
16. How can relative phase between two amplitudes affect measurement statistics even when the absolute global phase cannot? [CT1, CT2, G, M, S]
17. In a two-dimensional Hilbert space, how can every normalized pure state be parameterized up to global phase? [CT1, CT2, G, M, S]
18. What physical information is encoded in the Bloch-vector representation of a qubit or spin-1/2 state? [CT1, CT2, G, M, S]
19. How does a basis-state preparation differ from a measurement that merely finds the system in that basis state? [CT1, G, S]
20. What does it mean for two pure states to be orthogonal, and why does orthogonality matter operationally? [CT1, G, M, S, Sh]
21. How does one compute the transition probability between two normalized pure states? [CT1, G, M, S]
22. When are two apparently different superpositions physically indistinguishable? [CT1, G, M, S]
23. What is the relation between a pure state vector and the rank-one projector that represents the same pure state? [CT1, G, M, S]
24. How can the state of a finite system be specified without choosing a privileged measurement basis? [CT1, S, Sh, M]

## C. Operators, Observables, Spectra, and Commutators

25. What conditions make a linear operator a valid observable in a finite-dimensional Hilbert space? [CT1, G, S, Sh]
26. Why must the eigenvalues of a Hermitian operator be real? [CT1, G, S, Sh]
27. Why are eigenvectors of a Hermitian operator with distinct eigenvalues orthogonal? [CT1, G, S, Sh]
28. How does the spectral decomposition of a Hermitian operator encode its possible measurement outcomes? [CT1, G, S, Sh]
29. What changes in the spectral decomposition when an observable has degenerate eigenvalues? [CT1, G, S, Sh]
30. How is a projector onto a degenerate eigenspace different from a projector onto a single eigenvector? [CT1, G, S, Sh]
31. What is the physical meaning of an operator's matrix elements in a chosen basis? [CT1, S, Sh]
32. How can an operator be reconstructed from its action on a complete basis? [CT1, S, Sh]
33. What does it mean for two observables to commute, and why does this imply simultaneous diagonalizability in finite dimensions? [CT1, G, S, Sh]
34. What is a complete set of commuting observables in a finite-dimensional Hilbert space? [CT1, G, S]
35. How does degeneracy prevent one observable from uniquely labeling basis states? [CT1, G, S, Sh]
36. How does the commutator [A,B] quantify incompatibility of two observables? [CT1, G, S, Sh]
37. In what sense is the uncertainty relation in finite dimensions a statement about variances and commutators, not about measurement error? [CT1, G, S, Sh]
38. How does one define a function of a diagonalizable operator in finite dimensions? [CT1, Sh, S]
39. What is the difference between Hermitian, unitary, normal, and projection operators? [CT1, S, Sh, M]
40. How can a unitary operator be expressed through the exponential of a Hermitian generator? [CT1, CT2, G, S, Sh]

## D. Measurement, Probability, and State Update

41. Given a pure state and a Hermitian observable with nondegenerate spectrum, how are the probabilities of the possible outcomes computed? [CT1, G, S, Sh]
42. How is the probability rule modified when the measured eigenspace is degenerate? [CT1, G, S, Sh]
43. What is the post-measurement state after an ideal projective measurement with a nondegenerate result? [CT1, G, S]
44. What is the post-measurement state after an ideal projective measurement with a degenerate result? [CT1, G, S]
45. How can repeated measurement of the same observable yield certainty after the first ideal measurement? [CT1, G, S]
46. Why can a measurement of Sx disturb a previously prepared Sz eigenstate? [CT1, CT2, G, S]
47. How do sequential Stern-Gerlach experiments reveal noncommutativity of spin components? [CT1, CT2, G, S, Sh]
48. How is the expectation value of an observable computed from a pure state? [CT1, G, S, Sh]
49. How is the variance of an observable computed, and what does zero variance imply about the state? [CT1, G, S, Sh]
50. How can incompatible observables both have well-defined expectation values but not simultaneous sharp values? [CT1, G, S, Sh]
51. What is the distinction between preparing an ensemble and measuring a single system? [CT1, G, S]
52. How does a projective measurement in one basis transform coherence relative to another basis? [CT1, CT2, G, M, S]
53. How can a finite-dimensional measurement be represented by a family of mutually orthogonal projectors summing to the identity? [CT1, G, S, M]
54. What additional freedom is introduced by a generalized measurement described by POVM elements? [M, S]
55. How can two different measurement procedures have the same outcome probabilities but different post-measurement states? [M, S]
56. What is the operational meaning of distinguishability for two nonorthogonal pure states? [M, G, S]

## E. Unitary Dynamics and Two-Level Time Evolution

57. What condition must a time-evolution operator satisfy to preserve normalization in a finite-dimensional Hilbert space? [CT1, G, S, Sh]
58. How does a time-independent finite-dimensional Hamiltonian generate unitary evolution? [CT1, G, S, Sh]
59. How does one solve the Schrodinger equation by diagonalizing a time-independent Hamiltonian matrix? [CT1, G, S, Sh]
60. What physical information is contained in the relative phases accumulated by energy eigenstates? [CT1, G, S, Sh]
61. How does a two-level system exhibit oscillations when the Hamiltonian is not diagonal in the preparation basis? [CT1, CT2, G, M, S]
62. How does the energy splitting of a two-level Hamiltonian determine the oscillation frequency of transition probabilities? [CT1, CT2, G, S]
63. What finite-dimensional Hamiltonian describes a spin-1/2 magnetic moment in a uniform magnetic field? [CT1, CT2, G, S, Sh]
64. How does the Bloch vector rotate under that Hamiltonian, and which expectation values are conserved? [CT1, CT2, G, S]
65. Under what condition does a sinusoidal perturbation drive resonant transitions in a two-level system? [CT1, CT2, G]
66. How does the rotating-wave or near-resonance approximation simplify driven two-level dynamics? [CT2, G]
67. What distinguishes exact Rabi oscillations from first-order time-dependent perturbation theory? [CT2, G]
68. How does adiabatic evolution of a nondegenerate eigenstate differ from sudden evolution? [CT2, G, Sh]
69. What is Berry phase for a cyclic adiabatic evolution in a finite-dimensional eigenspace? [G, Sh, S]
70. How can a unitary gate be interpreted as a controlled finite-time Hamiltonian evolution? [M, S]

## F. Spin-1/2, Pauli Matrices, and Qubits

71. How do the Pauli matrices represent the Cartesian components of spin-1/2? [CT1, CT2, G, M, S, Sh]
72. How can one derive the eigenvectors of Sx and Sy from the Sz eigenbasis? [CT1, CT2, G, S]
73. What measurement probabilities result when an Sz-up spin is measured along an arbitrary direction n? [CT1, CT2, G, S]
74. How does the Bloch sphere encode the expectation values of the three Pauli observables? [CT1, CT2, G, M, S]
75. What is the relation between SU(2) transformations of spinors and SO(3) rotations of physical axes? [CT1, CT2, G, S, Sh]
76. Why can a 2pi rotation change the sign of a spin-1/2 state vector without changing the physical state? [CT1, CT2, G, S, Sh]
77. What experimental consequences distinguish spin-1/2 from a classical vector with continuously variable projection? [CT1, CT2, G, S]
78. How does a Stern-Gerlach apparatus function as both a state preparer and a measurement device? [CT1, CT2, G, S, Sh]
79. How does one construct the spin operator S.n for an arbitrary unit vector n? [CT1, CT2, G, S]
80. What are the eigenvalues and eigenvectors of S.n for spin-1/2? [CT1, CT2, G, S]
81. How does a magnetic field gradient couple spin to spatial path in a Stern-Gerlach measurement? [CT2, G, S, Sh]
82. How can a spin measurement along one axis be represented as a projective measurement in the eigenbasis of another spin operator? [CT1, CT2, G, S, Sh]
83. What is the difference between a spinor and an ordinary three-dimensional vector under rotations? [CT1, CT2, G, S, Sh]
84. How do Pauli matrix commutation and anticommutation relations constrain spin-1/2 observables? [CT1, CT2, G, M, S]
85. How does one express an arbitrary 2x2 Hermitian Hamiltonian as a scalar term plus an effective magnetic field term? [CT1, CT2, M, S]
86. How does the qubit formalism generalize the spin-1/2 formalism beyond physical spin? [M, CT1, G]

## G. Angular Momentum and Rotations

87. What algebraic commutation relations define angular momentum operators in quantum mechanics? [CT1, CT2, G, S, Sh]
88. How do ladder operators determine the allowed m values for a fixed angular momentum j? [CT1, CT2, G, S, Sh]
89. Why must an irreducible angular-momentum representation have dimension 2j+1? [CT1, CT2, G, S, Sh]
90. How does one construct the standard basis |j,m> for a finite angular-momentum irrep? [CT1, CT2, G, S, Sh]
91. How are the matrix elements of Jx and Jy obtained from J+ and J-? [CT1, CT2, G, S, Sh]
92. What is the difference between orbital angular momentum representations and abstract spin representations? [CT1, CT2, G, S, Sh]
93. Why do orbital angular momenta have integer l while spin may have half-integer j? [CT1, CT2, G, S, Sh]
94. How does a rotation operator act on a finite-dimensional spin state? [CT1, CT2, G, S, Sh]
95. How can rotational invariance of a Hamiltonian imply angular momentum conservation? [CT1, CT2, G, S, Sh]
96. What does it mean for a Hamiltonian to commute with J^2 and Jz but not separately with Lz and Sz? [CT2, G, S, Sh]
97. How does spin-orbit coupling make total angular momentum a better quantum number than spin or orbital angular momentum separately? [CT2, G, S, Sh]
98. How are Clebsch-Gordan coefficients used to change from an uncoupled product basis to a coupled angular-momentum basis? [CT2, G, S, Sh]
99. What total angular momenta can arise from adding j1 and j2? [CT2, G, S, Sh]
100. How does the singlet-triplet decomposition arise from adding two spin-1/2 systems? [CT2, G, S, Sh]
101. How does rotational symmetry restrict the matrix elements of vector operators? [G, S, Sh]
102. What is the finite-dimensional content of the Wigner-Eckart theorem? [S, Sh]
103. How do selection rules follow from angular-momentum conservation and tensor-operator structure? [G, S, Sh]
104. How can accidental degeneracy be distinguished from degeneracy enforced by angular momentum symmetry? [G, S, Sh]

## H. Composite Systems, Tensor Products, and Entanglement

105. Why is the state space of two finite quantum systems the tensor product, not the Cartesian product, of their individual Hilbert spaces? [CT1, CT2, G, M, S, Sh]
106. How does the dimension of a composite Hilbert space depend on the subsystem dimensions? [CT1, M, S, Sh]
107. What distinguishes a product state from an entangled state? [CT1, CT2, G, M, S]
108. How can one decide whether a pure bipartite state is entangled by examining its coefficient matrix? [M, S, Sh]
109. What is the Schmidt decomposition of a finite-dimensional bipartite pure state? [M, S]
110. How do Schmidt coefficients quantify entanglement for a pure bipartite state? [M, S]
111. How does measurement on one subsystem update the conditional state of another subsystem in an entangled state? [G, M, S]
112. What is the difference between classical correlation and quantum entanglement? [G, M, S]
113. Why can entanglement generate correlations that cannot be explained by predetermined local spin values? [G, M, S]
114. How does the Bell singlet state encode perfect anticorrelation along any common spin-measurement axis? [G, M, S]
115. How does the tensor-product structure make it possible for local observables on different subsystems to commute? [CT1, CT2, M, S]
116. What is the partial trace, and how does it produce a subsystem density operator? [M, S]
117. How can a subsystem of a pure entangled state be in a mixed state? [M, S]
118. What does it mean for a unitary operation to be local with respect to a tensor-product decomposition? [M, S]
119. How can an entangling unitary transform product states into entangled states? [M, S]
120. Why is CNOT an entangling gate even though it maps computational basis states to computational basis states? [M]
121. How does exchange symmetry constrain the state space of identical bosons or fermions? [CT2, G, S, Sh]
122. Why must the total state of identical fermions be antisymmetric under particle exchange? [CT2, G, S, Sh]
123. How does spin symmetry affect the allowed spatial states of two identical spin-1/2 fermions? [CT2, G, S, Sh]
124. How do singlet and triplet spin states determine symmetric or antisymmetric spatial wavefunctions for two electrons? [CT2, G, S, Sh]

## I. Density Operators, Mixed States, and General Measurements

125. What properties characterize a valid density operator on a finite-dimensional Hilbert space? [M, S]
126. How does a rank-one density operator represent a pure state? [M, S, CT1]
127. How can the purity Tr(rho^2) distinguish pure states from mixed states? [M, S]
128. How does a classical ensemble of pure states define a density operator? [M, S]
129. Why can different ensembles represent the same density operator? [M, S]
130. How are expectation values computed using the trace rule Tr(rho A)? [M, S]
131. How does unitary time evolution act on a density operator? [M, S]
132. How does a projective measurement update a density operator when the outcome is known? [M, S, CT1]
133. How does a projective measurement update a density operator when the outcome is ignored? [M, S]
134. What is the difference between decoherence in a basis and ignorance about a pure state? [M, S]
135. How can a POVM describe measurement probabilities without specifying a unique post-measurement state? [M, S]
136. How can a finite-dimensional quantum channel be represented by Kraus operators? [M, S]
137. What constraints must Kraus operators satisfy for a channel to preserve trace? [M, S]
138. How do depolarizing, dephasing, and bit-flip channels act on a qubit density matrix? [M]

## J. Degeneracy, Perturbations, and Effective Finite Subspaces

139. Why does ordinary nondegenerate perturbation theory fail inside a degenerate eigenspace? [CT2, G, S, Sh]
140. How does degenerate perturbation theory reduce first-order corrections to diagonalizing a finite matrix in the degenerate subspace? [CT2, G, S, Sh]
141. How should one choose the zeroth-order basis when a perturbation breaks a degeneracy? [CT2, G, S, Sh]
142. What does it mean for a perturbation to lift degeneracy partially rather than completely? [CT2, G, S, Sh]
143. How can symmetry predict which matrix elements of a perturbation vanish? [CT2, G, S, Sh]
144. How does a two-level avoided crossing arise from a 2x2 effective Hamiltonian? [CT2, G, S, Sh]
145. How does coupling between two nearly degenerate states modify the energy eigenstates? [CT2, G, S, Sh]
146. How can the Zeeman effect be formulated as diagonalization in a finite angular-momentum multiplet? [CT2, G, S, Sh]
147. How does fine or hyperfine structure select a coupled angular-momentum basis? [CT2, G, S, Sh]
148. When is a finite-dimensional truncation of a larger Hilbert space physically justified? [CT2, G, S, Sh]
149. What errors can arise when a truncated finite subspace is treated as exactly closed under the Hamiltonian? [CT2, G, S]
150. How can time-dependent perturbation theory be interpreted as amplitude flow among a finite set of basis states? [CT2, G]

## K. Quantum Information, Gates, Protocols, and Algorithms

151. What makes a qubit a two-dimensional quantum system rather than a classical bit with probabilistic values? [M, CT1, G]
152. How is a one-qubit gate represented by a 2x2 unitary matrix? [M]
153. How is an n-qubit pure state represented in a 2^n-dimensional Hilbert space? [M]
154. Why must reversible closed-system quantum computation be unitary? [M]
155. How can the Hadamard gate create coherent superposition from a computational basis state? [M]
156. How do Pauli X, Y, and Z gates act on computational basis states and phases? [M]
157. How does a controlled gate act differently on basis states and superpositions? [M]
158. Why is a two-qubit entangling gate necessary for universal quantum computation? [M]
159. How can arbitrary single-qubit rotations and one entangling two-qubit gate generate universal quantum circuits? [M]
160. How does phase kickback encode information in a relative phase rather than a classical output bit? [M]
161. How does the quantum Fourier transform act on computational basis states? [M]
162. Why is periodicity extraction central to Shor's factoring algorithm? [M]
163. How does measurement of an output register affect the state used for period finding in Shor's algorithm? [M]
164. What is amplitude amplification in Grover search, viewed as rotations in a two-dimensional subspace? [M]
165. How does Grover's algorithm use two reflections to rotate amplitude toward the marked subspace? [M]
166. Why is the no-cloning theorem a consequence of linearity and inner-product preservation? [G, M]
167. How does quantum teleportation transfer an unknown qubit state using entanglement and two classical bits? [M]
168. How does superdense coding use one shared Bell pair to transmit two classical bits by sending one qubit? [M]
169. What distinguishes entanglement-assisted protocols from protocols using only classical shared randomness? [M]
170. How does a stabilizer measurement reveal error information without measuring the encoded quantum state directly? [M]
171. What is the difference between bit-flip, phase-flip, and combined Pauli errors on a qubit? [M]
172. How can a repetition code protect against one type of qubit error but not arbitrary errors? [M]
173. Why must quantum error correction preserve coherence among logical basis states? [M]
174. How do syndrome subspaces partition the Hilbert space in a quantum error-correcting code? [M]
175. How can fault-tolerant reasoning be framed in terms of finite-dimensional subspaces and projectors? [M]

## L. Conceptual and Foundational Questions in Finite Dimension

176. What does a finite-dimensional version of the measurement problem ask about state update and unitary evolution? [G, S, M]
177. How do Bell-type spin experiments rule out assigning predetermined values to all spin components under locality assumptions? [G, M, S]
178. What is the finite-dimensional content of contextuality: why can outcome assignments depend on the measurement context? [S, M]
179. How does the EPR argument appear when formulated with two spin-1/2 systems? [G, M, S]
180. What does the no-signaling principle say about local measurements on one part of an entangled finite-dimensional system? [G, M, S]
181. How can collapse be understood operationally as state-update of information without enabling faster-than-light signaling? [G, M, S]
182. What distinguishes a hidden-variable explanation of spin correlations from the quantum density-operator account? [G, M, S]
183. How does the quantum Zeno effect appear in repeated projective measurements of a finite-dimensional system? [G, S]
184. How can the classical limit of a large spin differ from the behavior of a single spin-1/2 system? [CT2, G, S, Sh]
185. Which features of finite-dimensional quantum theory are purely kinematic, and which require a Hamiltonian or dynamics? [CT1, CT2, G, M, S, Sh]

## Consolidation Notes

- Repeated textbook problems about computing a specific spin probability were consolidated into the general questions on arbitrary-axis spin measurement, sequential Stern-Gerlach experiments, and Pauli-matrix spectra.
- Repeated perturbation exercises were consolidated into questions about diagonalizing perturbations in degenerate finite subspaces.
- Infinite-dimensional wavefunction problems were excluded unless they reduced to a finite invariant subspace, a finite angular-momentum multiplet, a two-level approximation, or a finite spin/composite subsystem.
- Quantum-computation material from Mermin was retained broadly because it is natively finite-dimensional.
- OCR-derived books were included after full-page OCR: Griffiths and Shankar.
