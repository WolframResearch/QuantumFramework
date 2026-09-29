---
Template: Default
---

# A formal semantics for QuantumFramework

This document describes QuantumFramework (QF) 2.1.1 (kernel `b13a9077`, 2026-09-27). Sections 1 to 11 cover the finite-dimensional core: for each object, what QF stores, and what QF does with the stored information for each operation. The appendices hold the full nesting and return-head tables (A), the constructor edge cases (B), the phase-space formulas (C), and the placement against the published literature (D, the only part that cites).

One idea runs through the whole account: **QF represents states and transformations together with their bases and subsystem structure, and then uses that information to compose them, discard subsystems, and read out measurements.** A single running example carries it: one object whose meaning comes from its declared basis, not from its stored coefficient array alone.

## 1. States and operators

A QuantumFramework state, operator, channel, or measurement is a coefficient array together with the information needed to read it. That information is what makes the array mean something, so it comes first: the array alone does not fix what it represents.

Take six complex numbers. The same six can be

- one six-level system, `QuantumState[{c1, ..., c6}, 6]`, a qudit of dimension six;
- a qubit together with a qutrit, `QuantumState[{c1, ..., c6}, {2, 3}]`;
- a linear map from a qutrit to a qubit, the $2 \times 3$ matrix with rows $(c_1, c_2, c_3)$ and $(c_4, c_5, c_6)$, `QuantumOperator[{{c1, c2, c3}, {c4, c5, c6}}, {{1}, {2}}]`.

Nothing in the numbers chooses among these. The subsystem dimensions and the split into what the object takes in and what it puts out decide it, and QF stores exactly that: for each object, the dimensions of its subsystems, an output side, an input side, and the basis on each side.

The third reading, the $2 \times 3$ map, lowers dimension from three levels to two. One physical operation has exactly this shape: removing a photon applies the annihilation operator $a$, with $a|n\rangle = \sqrt{n}\,|n-1\rangle$, which on inputs of at most two photons is the $2 \times 3$ matrix sending $\{|0\rangle, |1\rangle, |2\rangle\}$ to $\{|0\rangle, |1\rangle\}$. This [photon subtraction](https://www.science.org/doi/abs/10.1126/science.1146204) is heralded, conditional on detecting the removed photon, so the map does not preserve the norm: it is the unnormalized operator applied to the state and renormalized after the herald.

For each kind of object, three things: its ordinary quantum meaning, what QF stores, and how QF uses what it stores.

A **state** is a ket $|\psi\rangle$ or a density operator $\rho$, and an **operator** $A$ is a linear map. An operator acts in three ways:

$$ |\psi\rangle \mapsto A|\psi\rangle, \qquad \rho \mapsto A\rho A^\dagger, \qquad \rho \mapsto \mathcal E(\rho). $$

$A|\psi\rangle$ is an operator on kets; $A\rho A^\dagger$ is the same operator acting on a density matrix by conjugation; $\mathcal E(\rho)$ is a general linear map on density matrices, which need not be a single conjugation, such as a channel or the generator of an open-system evolution. QF represents all three with `QuantumOperator`: $A|\psi\rangle$ by its matrix on kets, $A\rho A^\dagger$ by that matrix acting through conjugation, and $\mathcal E(\rho)$ by a matrix on the operator space.

Mathematically a density operator is an operator, $\rho = \rho^\dagger \ge 0$ with $\operatorname{Tr}\rho = 1$; the two differ in how QF stores them. QF stores a ket as its amplitudes with no input subsystem. It stores an ordinary operator by flattening its matrix into a single coefficient vector over the output and input subsystems together. It stores a density operator as a `QuantumState` carrying matrix coefficients, a ket index and its conjugate, with no input subsystem.

So the storage records two independent things: whether the coefficients are a flattened vector or a matrix, and how the subsystems split into output and input. A density operator and an operator are the same kind of mathematical object stored two different ways; conflating the storage with the object, or calling both "a vector," hides the split QF actually tracks.

A **channel** and a **measurement** are operators with one extra output subsystem set aside, and a **circuit** is a list of these placed on shared wires. The public objects and their meanings:

| Object | Quantum meaning |
|---|---|
| `QuantumState` | a ket, or a density operator |
| `QuantumOperator` | a linear map, with an output side and an input side |
| `QuantumChannel` | a completely positive map, stored as a dilation with an environment set aside |
| `QuantumMeasurementOperator` | a general measurement; a projective observable is stored as the Hermitian operator itself, while a basis, a named measurement, or a list of POVM effects is stored as a dilation that keeps the outcome index rather than tracing it, and `"SuperOperator"` is that dilation either way |
| `QuantumMeasurement` | the result of a measurement: outcome probabilities and post-measurement states |
| `QuantumCircuitOperator` | a list of the above on numbered wires, a diagram |

These heads are layers, each adding one thing to the one it contains:

| Layer | What it contributes |
|---|---|
| `QuditBasis`, `QuantumBasis` | the basis elements and subsystem dimensions on an output side and an input side, together with the picture and the parameters: everything needed to read an array as a state or an operator |
| `QuantumState` | a coefficient array together with that basis |
| `QuantumOperator` | a `QuantumState` placed on specified output and input wires |
| `QuantumChannel` | an operator with one extra output wire, an environment to be traced |
| `QuantumMeasurementOperator` | an operator whose extra output wire carries the outcome when it is present; a projective observable adds none until `"SuperOperator"` supplies it |
| `QuantumCircuitOperator` | operations and their wiring, from which a method computes a result |

The exhaustive nesting is in Appendix A. A head names a representation, not a physical object. A `QuantumState` need only have the right shape; it need not be normalized or positive, and the same head holds the trace-zero output of a Liouvillian as readily as a physical state. Whether an object is a valid state, a unitary, or a trace-preserving channel is a separate question the constructor does not ask.

Throughout what follows the running example is one operator on a qubit on wire 1, introduced next.

## 2. Coefficients and basis: stored numbers versus the represented map

The coefficients QF stores represent an operator in whatever basis is declared; the stored matrix can differ from the operator's matrix in the computational basis. Write $C$ for the stored coefficient matrix and $A$ for the same operator's computational-basis matrix. If the columns of $F_{\mathrm{out}}$ and $F_{\mathrm{in}}$ are the declared output and input basis vectors, then

$$A = F_{\mathrm{out}}\, C\, F_{\mathrm{in}}^{-1}.$$

This needs only that the basis matrices $F_{\mathrm{out}}$ and $F_{\mathrm{in}}$ are invertible; it does not need them orthonormal. A concrete case: store $C = \operatorname{diag}(1, -1)$ and declare the basis to be the eigenbasis of $X$, with columns $|{+}\rangle, |{-}\rangle$. Then $A = |{+}\rangle\langle{+}| - |{-}\rangle\langle{-}| = X$: the stored coefficients $\operatorname{diag}(1, -1)$ represent the operator $X$ in the $X$ basis, though the same $C$ in the computational basis would be $Z$. In QF this is $R = $ `QuantumOperator[{{1, 0}, {0, -1}}, "PauliX"]`; its stored coefficient matrix $C$ is $\operatorname{diag}(1, -1)$ and its computational matrix $A$ is $X$. $R$ is the running operator.

Two constructions use the basis for opposite purposes, and QF keeps them apart. Declaring a basis fixes what supplied numbers mean: `QuantumOperator[{{1, 0}, {0, -1}}, "PauliX"]` reads the numbers $\operatorname{diag}(1, -1)$ as coefficients in the $X$ basis, and so builds $X$. Converting an existing object to another basis preserves the object and recomputes its numbers: `QuantumOperator[QuantumOperator["X"], "PauliX"]` takes the operator $X$, already built, and re-expresses it in the $X$ basis, where its stored matrix becomes $\operatorname{diag}(1, -1)$. The first turns numbers into an operator; the second turns an operator into different numbers for the same operator. Both reach $R$.

A change of coordinates preserves the operator $A$; applying a gate changes the state. This is why QF reconciles bases before it composes: coefficients on a shared subsystem cannot be contracted until both sides are written in one basis. The `QuantumBasis` also records a picture, the convention that fixes how the stored coefficients are read. This account treats the Schrödinger, Heisenberg, and phase-space pictures; the kernel also stores an interaction picture, which is outside this account.

## 3. Maps on density matrices

The map $\rho \mapsto \mathcal E(\rho)$ on density matrices is the foundation for channels, evolution, and the phase-space picture. It cannot always be written as a single conjugation $\rho \mapsto A\rho A^\dagger$: a channel is a sum $\sum_m M_m \rho M_m^\dagger$, and the generator of a continuous evolution is a commutator plus dissipators. QF represents such a map the same way it represents an operator, as a matrix of coefficients. The step that makes this possible is to turn $\rho$ into a vector.

For a system of dimension $d$, vectorization arranges the $d^2$ elements of its $d \times d$ density matrix into a vector of length $d^2$, and a map on density matrices then becomes a $d^2 \times d^2$ matrix acting on that vector. That matrix is again a `QuantumOperator`, on the operator space of $d \times d$ matrices. It is how the map is represented; applying the map need not build it in full, since a tensor contraction can carry out the same action. Left and right multiplication are `QuantumOperator["Left"[A]]` $= A \otimes I$ and `QuantumOperator["Right"[A]]` $= I \otimes A^{\mathsf T}$; a conjugation $\rho \mapsto A\rho A^\dagger$ combines left multiplication by $A$ with right multiplication by $A^\dagger$; the generator of the von Neumann equation is `QuantumOperator["Liouvillian"[H]]`, the matrix of $-i[H, \cdot]$. The Liouvillian of $X$ applied to the density matrix $|0\rangle\langle 0|$ returns the commutator, a `QuantumState` of trace zero:

$$ \texttt{QuantumOperator["Liouvillian"["X"]]}\,\big[\,|0\rangle\langle 0|\,\big] = -i\,[X, |0\rangle\langle 0|] = \begin{pmatrix} 0 & i \\ -i & 0 \end{pmatrix}. $$

The stored coefficients and this matrix are the same numbers read with different indices. To build the matrix that acts on the vectorized density operator, QF reshapes the stored coefficient array and swaps its two index groups (the property `"StateMatrix"`); the numbers do not change, only their arrangement. Two steps are worth separating: `"Bend"` forms the density matrix, from a ket as $|\psi\rangle\langle\psi|$ and from a stored matrix as itself, then vectorizes it into the `"DensityVector"` of length $d^2$; `"Double"` additionally interleaves the axes so each subsystem sits beside its conjugate, the arrangement the maps above act on.

Three things now share one model. A channel is the map $\rho \mapsto \sum_m M_m \rho M_m^\dagger$, an operator-space matrix that QF stores as a dilation with an environment rather than as the matrix itself. A continuous evolution integrates the map generated by the Liouvillian, or by the Lindbladian when there are dissipators. The phase-space picture is a change of coordinates of this same matrix. Each is the same operator-space representation read or generated differently, not a separate mechanism.

## 4. Combining objects

Apply, compose, tensor, adjoint, and discard, each with what QF does to perform it.

#| colwidths: 0.3 0.1 0.5
| Operation | Meaning | What QF adds |
|---|---|---|
| Apply to a ket | $A\lvert\psi\rangle$ | reads the state on the input wires, carries untouched wires unchanged, reconciles the basis |
| Apply to a density matrix | $A\rho A^\dagger$ | the same, on the stored density coefficients |
| Compose ($B$ then $A$) | $AB$ | contracts the shared wire, aligning bases first |
| Tensor | $A \otimes B$ | places on disjoint wires, renumbering to avoid clashes |
| Adjoint | $A^\dagger$ | swaps input and output and conjugates; the physical adjoint only on orthonormal bases |
| Discard a subsystem | $\operatorname{Tr}_B \rho$ | contracts the discarded subsystem's row and column indices |

**Applying.** An operator acts on the wires it names and leaves the rest alone. Pauli $X$ on qubit 1 of a two-qubit register sends $|00\rangle$ to $|10\rangle$: QF puts an identity on qubit 2 and matches the dimension on qubit 1. If the state has fewer qubits than the operator reaches, QF pads it with $|0\rangle$ on the missing wires. On a density matrix the same operator gives $X\rho X^\dagger$. The result carries a basis of its own, so its stored amplitudes need not be its computational amplitudes; applying $R$ to $|0\rangle$ gives $|1\rangle$, stored as $(1, -1)/\sqrt2$ in the $X$ basis.

Applying runs as a short sequence, and the other operations share its shape. QF reads which wires the operator and the state have in common and their dimensions; writes both on a common basis, changing coordinates when the shared wires carry different ones; takes the amplitude form for a ket and the vectorized form for a density matrix; contracts the shared wires; and rebuilds the result as a new object carrying its own basis and wire orders. Compose, tensor, and trace are the same steps with a different contraction, and a method runs this for a whole circuit at once.

**Composing.** Applying $B$ then $A$ is the matrix product $AB$,

$$(AB)_{ik} = \sum_j A_{ij}\, B_{jk}.$$

The summed index $j$ runs over the basis of the shared subsystem, the one $A$ takes as input and $B$ produces as output. QF contracts that subsystem, keeps every subsystem the contraction does not touch (the remaining outputs and inputs of both operands), and on the shared subsystem writes both operators in one basis before summing. When the two operators overlap only partially, this is the general rule: composing a two-qubit gate on wires 1 and 2 with a one-qubit gate on wire 2 contracts wire 2 and leaves wire 1's output and the two-qubit gate's other input in place. Maps on density matrices compose by the same contraction on the operator space, so a channel after a channel, or a gate after a channel, is one more contraction.

That basis step is where the stored coefficients matter. Compose $Z$ after $R$, whose stored coefficient matrix is $\operatorname{diag}(1, -1)$ although the operator is $X$. Contracting the stored numbers as they sit would give $Z\operatorname{diag}(1,-1) = I$, the wrong operator. QF reconciles the bases and returns $ZX$,

$$ZX = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}.$$

**Adjoint.** $A^\dagger$ swaps the input and output sides and conjugate-transposes the coefficients. This is the physical adjoint only when the bases are orthonormal. In general, with Gram matrices $G = F^\dagger F$ on each side, the adjoint of the represented map has coefficients $G_{\mathrm{in}}^{-1}\, C^\dagger\, G_{\mathrm{out}}$, and QF omits the Gram factors, so on a non-orthonormal basis one first changes to the computational basis.

**Tensor product.** $A \otimes B$ places the two operators on disjoint wires, renumbering $B$'s wires so none clashes; the pictures are not checked for compatibility, and the first operand's picture is adopted.

**Discarding.** The partial trace has one rule per object. On a density matrix, discarding a subsystem sums over the matching row and column indices of that subsystem, $\operatorname{Tr}_B\rho$. On an operator, a trace pair contracts one output wire against one input wire; the two wires must carry the same basis, or the result is the trace of the stored coefficients rather than the trace of the map. On a circuit, the trace is drawn as a cap: a wire that closes the two ends into a loop.

**Circuits.** A circuit is the composition of its listed operations on numbered wires, applied in order. QF composes them by the same contraction, so the circuit is one object; a method turns the diagram into a value.

**The correctness requirement for composition.** Convert the composed object to its computational matrix and you get the product of the two operators' computational matrices, once both are placed on a common register with identities on the untouched subsystems and their bases are aligned. This holds for any invertible basis matrices on the shared subsystems, and it rests entirely on the change-of-basis step.

Composition needs the shared subsystems to carry the same dimensions. A dimension mismatch is refused. A mismatch between the Schrödinger and Heisenberg pictures is not refused; it is silently resolved to the Schrödinger picture. A phase-space operand instead fails on a dimension mismatch.

## 5. Channels and measurements

A channel and a measurement both come from a family of operators, but the family means different things and the list constructors read it differently. That difference fixes the physical object the user creates.

A **channel** is the map $\rho \mapsto \sum_m M_m \rho M_m^\dagger$, the $M_m$ its Kraus operators. `QuantumChannel[{M1, M2, ...}]` reads the list as Kraus operators directly. QF stores the family as one operator with an extra output subsystem,

$$V = \sum_m |m\rangle \otimes M_m,$$

the Kraus operators stacked along that index; discarding the index, the environment, and tracing it out returns the map. This is the Stinespring form: a dilation followed by a trace. $V$ is an isometry, and the map preserves the trace, exactly when $\sum_m M_m^\dagger M_m = I$, which QF does not check at construction and exposes as `"TracePreservingQ"`; a channel built from an arbitrary list is completely positive but need not preserve the trace. On the running qubit, the bit-flip channel `QuantumChannel[{Sqrt[3/4] I, Sqrt[1/4] X}]` sends $|0\rangle\langle 0|$ to $\operatorname{diag}(3/4, 1/4)$.

A **measurement** is read branch by branch. With measurement operators $M_m$ acting on a state $\rho$ by

$$\mathcal I_m(\rho) = M_m \rho M_m^\dagger, \qquad p_m = \operatorname{Tr}\mathcal I_m(\rho), \qquad \rho_m = \frac{\mathcal I_m(\rho)}{p_m}\ (p_m > 0),$$

four things are distinct and QF keeps them distinct: the unnormalized branch $\mathcal I_m(\rho)$, its probability $p_m$, the conditional state $\rho_m$, and the state $\sum_m \mathcal I_m(\rho)$ obtained by ignoring the outcome. The probabilities alone come from the effects $E_m = M_m^\dagger M_m$, a positive operator-valued measure (POVM), through the Born rule $p_m = \operatorname{Tr}(E_m \rho)$; the POVM fixes the probabilities but not the disturbance, since many measurement operators share one effect.

The list constructors read a measurement's family as effects, not as Kraus operators. `QuantumMeasurementOperator[{E1, E2, ...}]` takes the $E_m$ to be POVM effects and forms each measurement operator as the matrix square root $M_m = \sqrt{E_m}$, the canonical choice. The stored measurement operators are these $\sqrt{E_m}$. `"POVMElements"` returns the effects, rescaled so their sum has unit mean diagonal; it returns them verbatim only when the supplied effects already sum that way. An observable is its own construction: `QuantumMeasurementOperator[obs]` stores the Hermitian operator itself, with no outcome wire, and measures it in its eigenbasis; its outcomes are the operator's eigenvalues, where a basis measurement instead labels outcomes by basis index. The dilation, and the extra outcome wire, appear through `"SuperOperator"`, and the two encodings give the same outcome probabilities. So the same list of matrices builds a different object under `QuantumChannel`, where it is Kraus operators, than under `QuantumMeasurementOperator`, where it is effects.

The measurement stores its family as the same $V = \sum_m |m\rangle \otimes M_m$, now with the extra index kept to carry the outcome. This is where the QF-specific semantics lives: the extra index can be kept as a live wire, a bare isometry with nothing measured yet; discarded, the channel above; or read, the measurement. Kept, the branches stay coherent in $V$ and record no classical outcome; the outcome becomes classical only once the index is read.

Read, applying the measurement returns a `QuantumMeasurement` exposing `"Probabilities"`, `"States"`, and `"PostMeasurementState"`. `"Probabilities"` is the Born distribution, an association keyed by outcome; the bare vector is `"ProbabilitiesList"`, and returning either is not sampling a single outcome. `"States"` returns the branch states on the full register, not the normalized conditional states $\rho_m$. For a pure input each branch comes back as the unnormalized ket $M_m|\psi\rangle$, of norm $\sqrt{p_m}$. Its density operator is the branch $\mathcal I_m(\rho) = M_m|\psi\rangle\langle\psi|M_m^\dagger$, of trace $p_m$. Measuring one qubit of a Bell state returns two branch kets of norm $1/\sqrt2$, each the unnormalized half of the split with probability $1/2$. `"PostMeasurementState"` reduces to the target subsystems, tracing out the outcome index and every non-target subsystem, so it is not the full-register object `"States"` returns.

`R["Matrix"]` is the stored $\operatorname{diag}(1, -1)$ and `R["MatrixRepresentation"]` is its computational matrix $X$; applying it, `R[QuantumState["0"]]` is $|1\rangle$, stored $(1, -1)/\sqrt2$ in the $X$ basis. Measuring the $X$ observable on that state, `QuantumMeasurementOperator["X"][R[QuantumState["0"]]]`, returns a `QuantumMeasurement` with probabilities $\{1/2, 1/2\}$ and two branch kets of norm $1/\sqrt2$, the split of $|1\rangle = (|{+}\rangle - |{-}\rangle)/\sqrt2$.

## 6. Parameters and evolution

A parameterized operator is a family indexed by the values of its parameters, $\theta \mapsto A(\theta)$; QF records the parameter names and their ranges alongside the coefficients, and composing two families composes them value by value on the union of their names. Substituting a parameter value is evaluation at a point, $A(\theta_0)$. Because the host is a computer-algebra system, an amplitude may itself stay symbolic, so a value can enter either as a declared parameter or as a symbolic amplitude; which one to use is a modeling choice.

Continuous-time evolution is `QuantumEvolve`, and it acts on the representations already built. From a Hamiltonian $H$ it integrates the Schrödinger equation on the state; from a Hamiltonian, a list of Lindblad operators, and rates it integrates the map generated by the Lindbladian, so the evolved object is a density-matrix `QuantumState`. With no initial state it returns the propagator as an operator, and from an initial observable it returns the Heisenberg-picture evolution as an operator. The equation is solved when `QuantumEvolve` is called, since the solver runs inside the constructor; only the evaluation of the stored solution at a chosen time is deferred. On the running qubit, `QuantumEvolve[QuantumOperator["X"], {QuantumOperator["Z"]}, ρ, {t, 0, 1}]` returns a state whose density matrix carries the time $t$.

## 7. Coordinate changes: phase space

A phase-space picture is a change of coordinates of the operator-space representation. Until now a basis has been the state vectors that expand a ket; here it is an operator basis, a set of matrices $B_a$ that expand a density matrix as $\rho = \sum_a c_a B_a$. For one qubit the four matrices $I, X, Y, Z$ are such a basis. Which picture the `QuantumBasis` records fixes what the coefficients are: a ket's amplitudes, a density matrix's matrix elements in the Schrödinger picture, and in a phase-space picture the operator-basis coefficients $c_a$, which are quasi-probabilities; the Heisenberg picture instead carries the evolution on operators rather than states. QF reaches a phase-space picture with `QuantumPhaseSpaceTransform`: read the operator-space coefficients in an operator basis with $d^2$ elements per subsystem. Because it is a change of coordinates of the operator-space matrices, it preserves composition and application: transforming a product gives the product of the transforms. When both the density matrix and the operator basis are Hermitian, as the Wigner, Pauli, and Gell-Mann bases are, the coefficients $c_a$ are real, though the vectorized matrix is complex.

A complete operator basis always gives coordinates, its expansion coefficients fixing the density matrix; a measurement's outcome probabilities fix it only when the measurement is informationally complete, which an arbitrary POVM need not be. The Wigner basis (displaced parities in odd dimension, and a half-phase variant in even dimension) gives the discrete Wigner function; a Hermitian operator basis gives real Bloch-like components; an informationally complete POVM, such as the symmetric informationally complete (SIC) POVM of the named `"QBismSIC"` basis, gives its outcome probabilities as coordinates, built from the POVM's dual operators. The inverse, `QuantumWeylTransform`, takes square roots of the dimensions and returns to the Schrödinger picture; it recovers a density matrix exactly but forgets the global phase of a ket or a unitary, because the density matrix $\rho = \psi\psi^\dagger$ does not carry it. The rate matrix `HamiltonianTransitionRate[H]` is the Liouvillian $-i[H,\cdot]$ expressed directly in these coordinates. The even-dimension convention, the weighted trace-preservation condition, and the sign-sector transform are in Appendix C.

## 8. The stabilizer representation

A qubit stabilizer state is a different data structure, a tableau: the group of Pauli operators that fixes the state, stored as bits, with one recorded global phase. It is not a tensor of amplitudes, and QF joins it to the tensor representation by two conversions. `QuantumState[ps]` builds the state vector from the tableau and multiplies in the recorded global phase; its cost is set by the construction it uses, a $2^n \times 2^n$ operator and its null space, not by the $2^n$ amplitudes it returns. `PauliStabilizer[qs]` converts a state to a tableau, succeeds only on stabilizer states, and records the global phase it finds. Because the phase is recorded, converting a state to a tableau and back returns the same physical state including its global phase; the reverse round trip returns the same stabilizer group, though not necessarily the same generator list. The conversions keep the phase; only the Clifford gate updates inside the tableau act up to global phase.

The boundary is sharp. Clifford gates, including the $\pi/2$ phase gate $S$, keep the tableau. A non-Clifford phase gate such as $T$ leaves it for a sum of stabilizer states. A non-Clifford non-phase unitary is refused.

A separate result connects the two representations in odd dimension; it concerns qudit states, not the qubit tableau above. A pure qudit state has a nonnegative discrete Wigner function exactly when it is a stabilizer state. So Wigner negativity appears only once a state leaves the stabilizer set. It appears only on states that a non-Clifford gate actually moves: a gate diagonal in the computational basis fixes its own stabilizer eigenstates, and those stay nonnegative. Qubits are even-dimensional, and there Wigner nonnegativity is not this criterion.

## 9. Execution methods

A circuit is a diagram until a method turns it into a value, and the methods do not all return the same kind of value. The default method contracts the circuit as one network of tensors joined at its wires, in a chosen order. The Schrödinger method instead multiplies the gates one at a time. The stabilizer method runs in the tableau. These are exact simulators: each must return, after decoding, the state the circuit's operator produces, so on any one circuit they must agree. An external method hands the circuit to a simulator or to hardware, and there the return is not a state but what the backend gives back: a table of shot counts, or a job to be queried. The exact methods agree on the output state; a shot-based or hardware run returns samples of it, not the state.

## 10. Distances and entanglement

Two families read a number off the density-matrix representation rather than build a new object.

A distance compares two states. `QuantumDistance[a, b]` defaults to the `"Fidelity"` measure. Despite the name this measure is a distance, $1 - \operatorname{Re}\operatorname{Tr}\sqrt{\rho_a \rho_b}$, zero on identical states and not the Uhlmann fidelity itself. A third argument selects another measure instead: trace, Bures, Bures-angle, Hilbert-Schmidt, relative-entropy, relative-purity, or Bloch. `QuantumSimilarity` is the companion quantity $1 - d$.

An entanglement monotone reads one state across a bipartition. `QuantumEntanglementMonotone[a, part]` defaults to the concurrence, taken from the reduced state for a pure input and from the concurrence vector for a mixed one, and exact only for two qubits or a pure state. Its alternatives are negativity, log-negativity, the entanglement and Renyi entropies, realignment, mutual information, and discord. `QuantumEntangledQ` tests entanglement and defaults to the realignment criterion.

All of these are functionals on the state already represented, so they inherit its shape assumptions. Nothing here checks that the input is normalized or positive.

## 11. Assumptions and limitations

QF checks shapes, not physics. A constructor accepts a vector or a square matrix of the right dimension; it does not check normalization, positivity, unitarity, or trace preservation, which are assumptions of the model, readable afterward but not enforced. So a `QuantumState` need not be normalized, and a `QuantumChannel` need not preserve the trace.

Two operations are computed directly on the coefficient matrix in the declared basis, and each matches the physical operation only under a condition. The adjoint conjugate-transposes the coefficient matrix, which is the Hermitian adjoint $A^\dagger$, fixed by $\langle\phi|A\psi\rangle = \langle A^\dagger\phi|\psi\rangle$, only in an orthonormal basis; in a non-orthonormal basis the overlaps of the basis vectors, collected in the Gram matrix, enter, and the bare conjugate-transpose omits them. The partial trace contracts a subsystem's output index with its input index, and because the trace of an operator does not depend on the basis, this reproduces the physical partial trace $\operatorname{Tr}_S(\cdot) = \sum_i \langle i|\,\cdot\,|i\rangle$ whenever the subsystem's two sides carry the same basis; when they carry different bases the contraction pairs vectors that are not dual to each other, and the result is not the trace. In both cases QF still returns a well-defined object, but off these conditions it is not the physical adjoint or partial trace.

Two identifications must not be over-read. Two kets that differ by a phase are the same physical state, and two normalized density matrices equal as operators are the same state; the freedom to rescale an unnormalized positive operator to unit trace is a chosen normalization, not a general equality of QF objects. In particular a subnormalized measurement branch is not identified up to an arbitrary positive factor, because its trace is its probability. And the finite-dimensional compositional core described here does not exhaust the kernel: the bosonic and second-quantization layers carry their own objects and are outside this account.

## Appendix A. Objects and their nesting

The heads wrap one another; there is no second data structure in the chain. A `QuantumOperator` contains a `QuantumState`; a `QuantumChannel` and a `QuantumMeasurementOperator` contain a `QuantumOperator`; a `QuantumMeasurement` contains a `QuantumMeasurementOperator` after application. The stabilizer heads (`PauliStabilizer`, `StabilizerFrame`, `CliffordChannel`, `GraphState`) are the one other representation, joined to the tensor representation by two conversions.

The head that comes back from each call, run against the working tree. A lowercase name stands for any object of the matching head: `qs` a `QuantumState`, `qo` a `QuantumOperator`, `qmo` a `QuantumMeasurementOperator`, `qm` a `QuantumMeasurement`, `chan` a `QuantumChannel`, `qco` and `qc` a `QuantumCircuitOperator`, and `ps` a `PauliStabilizer`; `op` is any operation a circuit can hold.

| Call | Head of result |
|---|---|
| `qo[qs]` | `QuantumState` |
| `qo1[qo2]` | `QuantumOperator` |
| `qo[qmo]`, `qmo[qo]` | `QuantumMeasurementOperator` |
| `qo[chan]` | `QuantumChannel` |
| `chan[qs]` | `QuantumState` |
| `qmo[qs]` | `QuantumMeasurement` |
| `qco[qs]`, `qco[]` without a measurement | `QuantumState` |
| `qco[qs]`, `qco[]` with a measurement | `QuantumMeasurement` |
| `qco[op]`, `qc1 /* qc2` | `QuantumCircuitOperator` |
| `QuantumPartialTrace` of a state, operator, or circuit | same head as the argument |
| `SuperDagger[qo]`, `qo1 + qo2`, `QuantumTensorProduct[qo1, qo2]` | `QuantumOperator` |
| `QuantumEvolve[...]` | `QuantumState` or `QuantumOperator` |

## Appendix B. Constructor edge cases and wire numbers

The integer on each wire records its role by its sign: a positive number is a live quantum wire; a non-positive number is a trace subsystem inside a `QuantumChannel` and an outcome subsystem inside a `QuantumMeasurementOperator`, and carries no role on a bare operator. The role is assigned by the wrapper, so the same coefficients serve as a unitary, a channel, or a measurement, and the role is invisible in a bare wire list.

When an operator reaches wires the state does not fill, the state is padded with $|0\rangle$ on the missing wires, up to the full input order of the one-element circuit, gaps included; the padding is by qubits, so a missing qutrit fails with a dimension message rather than being prepared. A coefficient array whose dimension factors, such as six, is split into subsystems by default (`QuantumState[Range[6]]` is a qubit and a qutrit, not one six-level qudit); pass the dimension explicitly to override. A tableau is stored bit-packed, sixty-two qubits to a machine word, with one recorded global phase (the sixty-two-rows-per-word layout is a transient packing used only by the compiled gate fold).

## Appendix C. Phase-space formulas

In odd dimension the Wigner basis element at $(q,p)$ is a displaced parity, built from three operators on a $d$-level system with $\omega = e^{2\pi i/d}$: the Fourier operator $F\lvert k\rangle = \tfrac{1}{\sqrt d}\sum_j \omega^{jk}\lvert j\rangle$, the clock $Z\lvert k\rangle = \omega^k\lvert k\rangle$, and the shift $X\lvert k\rangle = \lvert k-1 \bmod d\rangle$. The parity is $F^2$, with $F^2\lvert k\rangle = \lvert{-k}\bmod d\rangle$, and the displaced parity is $A(q,p) = e^{i\pi q p/d}\, F^2\, Z^q X^p$, which the kernel builds directly in this form rather than through a separate displacement operator. The Wigner basis element at $(q,p)$ is $A(2q, 2p)$, reaching every grid point because $2$ is invertible mod $d$; the elements are Hermitian, unit trace, mutually orthogonal with $\operatorname{Tr}[A_a A_b] = d\,\delta_{ab}$, summing to $d\cdot I$, and the coordinates are $W(q,p) = \operatorname{Tr}[A(2q,2p)\,\rho]/d$. In even dimension $2$ is not invertible, the parity orbit is too small, and the element is $2\,A(q,p)$ instead, with Gram matrix $4d\cdot I$; for a qubit these are $2\{I, X, Z, -Y\}$. The elements no longer share a trace, so $\operatorname{Tr}\rho = \sum_a \tau_a c_a$ with $\tau_a = \operatorname{Tr} A_a$, which for a qubit is $4c_0$. The inverse `QuantumWeylTransform` takes the square root of each dimension ($9 \to 3$ for a qutrit).

The rate matrix `HamiltonianTransitionRate[H]` is the matrix of $-i[H,\cdot]$ in these coordinates with no further factor. Trace preservation is the weighted condition $\tau^T R = 0$ with $\tau_a = \operatorname{Tr} A_a$, since $\operatorname{Tr}\rho = \tau^T c$; it holds for a valid generator. The ordinary column sums $\mathbf 1^T R$, with $\mathbf 1$ the all-ones vector, are a different quantity. When all basis elements share one nonzero trace, $\tau \propto \mathbf 1$, so $\mathbf 1^T R$ is proportional to $\tau^T R$ and vanishes, as for the odd-dimension displaced parities (zero for the qutrit); with unequal traces it need not, and the even-dimension qubit basis has $\mathbf 1^T R = (0, 0, 2h, -2h)$ for $H = hX$. The sign-sector transform `QuantumPositiveTransform` splits the coordinate axis into a nonnegative pair $(w^+, w^-)$ with $w = w^+ - w^-$; on a matrix it yields the block form $[[T^+, T^-],[T^-, T^+]]$, and the decoding $\Delta = (I, -I)$ satisfies $\Delta\,P(T) = T\,\Delta$.

## Appendix D. Related work

The nearest published model for QF is a composite: a categorical process model with an environment and classical structure as the meaning, and a Proto-Quipper-style typed core as the language. The categorical line [Abr04, Sel07, Coe06, Coe10, Coe14] is the closest reading of the channels and measurements: the doubling `"Bend"` corresponds to the completely positive maps (CPM) construction, the cup and cap to a dagger-compact structure, the spiders to basis structures, and a measurement to a morphism to a classical object. The programming-language line [Rio17, Fu22, Pay17, Sel04, Sel04b, Ros26] is the closest reading of the circuit shape: a circuit as a morphism in a free category with a functor for its meaning. The operator-algebraic line [Cho14, Jia21, Lin25c, Ran21] is where the picture structure lives, the Schrödinger description evolving states and the Heisenberg description evolving observables. Phase space is Wootters and Gross [Woo87, Gro06, Gib04], the rate matrix is Braasch and Wootters [Bra20], and the tableau is Aaronson and Gottesman [Aar04, Got98, Nes09].

Against that composite, QF combines four things no surveyed system combines: a basis stored with every state and operator, with an implicit change of basis rather than a rejection (Qwerty [Ada24] and Díaz-Caro [Dia25b, Dia25] put the basis in the type instead; Díaz-Caro's core carries an implicit change of basis of its own, but in the type rather than stored with the object); three pictures in one object model, connected by transforms, with the phase-space picture unrepresented in any surveyed system; channels stored as a Stinespring operator with the environment as trace subsystems; and a second representation, the stabilizer tableau, joined by conversions that carry the global phase. Yao.jl [Luo19] shares the unified block type but has none of the four; the idealized cores lambda-Q# [Sin22], QIR [Luo23], and Qunity [Voi22] formalize a designed subset rather than a shipping framework.

This comparison records what each cited work builds into its published model, not a limit of its approach: a feature a work does not present may be a matter of scope rather than reach, since for a designed core or a categorical axiomatization the construction could often be added within the same framework. It is a comparison against published formal models; coexisting statevector and stabilizer representations with conversions between them are also routine in circuit simulators such as Stim, Qiskit Aer, and Cirq, so what is distinctive in QF is not that the two coexist but that they form one object model, with the phase-carrying round trip stated as a property.

A vocabulary note, since this appendix speaks the language of the works it cites. A *functor* is a representation-to-representation map that preserves composition, of the kind this account calls a transform; a *morphism* is a linear map, channel, or state; a *dagger-compact category* is finite-dimensional states and operators with tensor product and adjoint; the *CPM construction* is the doubling; a *Frobenius algebra*, or *spider*, is the copy-and-compare structure of a basis.

The deep search behind this survey is at https://app.undermind.ai/projects/bcedc042-7a0d-4f3d-b4ca-4e1c138571a3?path=/Formal%20semantics%20and%20object%20models%20of%20quantum%20software%20frameworks .

- [Aar04] S. Aaronson, D. Gottesman, "Improved simulation of stabilizer circuits," Phys. Rev. A 70 (2004) 052328. doi:10.1103/PhysRevA.70.052328
- [Abr04] S. Abramsky, B. Coecke, "A categorical semantics of quantum protocols," LICS 2004. doi:10.1109/LICS.2004.1319636
- [Ada24] A. J. Adams, S. Khan, J. S. Young, T. Conte, "Qwerty: A Basis-Oriented Quantum Programming Language," IEEE QCE 2025, pp. 804-815. doi:10.1109/QCE65121.2025.00093
- [Ber18] V. Bergholm et al., "PennyLane: Automatic differentiation of hybrid quantum-classical computations," arXiv:1811.04968.
- [Bra20] W. F. Braasch Jr., W. K. Wootters, "Transition probabilities and transition rates in discrete phase space," Phys. Rev. A 102 (2020) 052204. doi:10.1103/PhysRevA.102.052204
- [Cho14] K. Cho, "Semantics for a Quantum Programming Language by Operator Algebras," New Generation Computing 34 (2016) 25-68. doi:10.1007/s00354-016-0204-3
- [Cle22] A. Clément, N. Heurtel, S. Mansfield, S. Perdrix, B. Valiron, "A Complete Equational Theory for Quantum Circuits," LICS 2023. doi:10.1109/LICS56636.2023.10175801
- [Cle23b] A. Clément, N. Delorme, S. Perdrix, R. Vilmart, "Quantum Circuit Completeness: Extensions and Simplifications," CSL 2024. doi:10.4230/LIPIcs.CSL.2024.20
- [Coe06] B. Coecke, D. Pavlovic, "Quantum measurements without sums," in Mathematics of Quantum Computation and Quantum Technology (2007). doi:10.1201/9781584889007.ch16
- [Coe09] B. Coecke, R. Duncan, "Interacting quantum observables: categorical algebra and diagrammatics," New J. Phys. 13 (2011) 043016. doi:10.1088/1367-2630/13/4/043016
- [Coe10] B. Coecke, S. Perdrix, "Environment and Classical Channels in Categorical Quantum Mechanics," Log. Methods Comput. Sci. 8(4:14) (2012). doi:10.2168/LMCS-8(4:14)2012
- [Coe13] B. Coecke, C. Heunen, A. Kissinger, "Categories of quantum and classical channels," Quantum Inf. Process. 15 (2016) 5179-5209. doi:10.1007/s11128-014-0837-4
- [Coe14] B. Coecke, C. Heunen, A. Kissinger, "Categories of Quantum and Classical Channels (extended abstract)," QPL 2012, EPTCS 158, pp. 1-14. doi:10.4204/EPTCS.158.1
- [Cro21] A. W. Cross et al., "OpenQASM 3: A Broader and Deeper Quantum Assembly Language," ACM TQC 3(3) (2022). doi:10.1145/3505636
- [Dia25] A. Díaz-Caro, N. A. Monzon, "A Quantum-Control Lambda-Calculus with Multiple Measurement Bases," APLAS 2025, pp. 151-170. doi:10.1007/978-981-95-3585-9_8
- [Dia25b] A. Díaz-Caro, O. Malherbe, R. Romero, "Basis-Sensitive Quantum Typing via Realisability," arXiv:2510.18542. doi:10.48550/arXiv.2510.18542
- [Fu20b] P. Fu, K. Kishida, P. Selinger, "Linear Dependent Type Theory for Quantum Programming Languages," LICS 2020. doi:10.1145/3373718.3394765
- [Fu22] P. Fu, K. Kishida, N. J. Ross, P. Selinger, "Proto-Quipper with Dynamic Lifting," Proc. ACM Program. Lang. (POPL 2023) 309-334. doi:10.1145/3571204
- [Fu24] P. Fu, K. Kishida, N. J. Ross, P. Selinger, "Proto-Quipper with Reversing and Control," QPL 2024, EPTCS 426. doi:10.4204/EPTCS.426.1
- [Gib04] K. S. Gibbons, M. J. Hoffman, W. K. Wootters, "Discrete phase space based on finite fields," Phys. Rev. A 70 (2004) 062101. doi:10.1103/PhysRevA.70.062101
- [Got98] D. Gottesman, "The Heisenberg Representation of Quantum Computers," arXiv:quant-ph/9807006 (1998).
- [Gre13] A. S. Green, P. L. Lumsdaine, N. J. Ross, P. Selinger, B. Valiron, "Quipper: a scalable quantum programming language," PLDI 2013. doi:10.1145/2499370.2462177
- [Gro06] D. Gross, "Hudson's theorem for finite-dimensional quantum systems," J. Math. Phys. 47 (2006) 122107. doi:10.1063/1.2393152
- [Jia21] X. Jia, A. Kornell, B. Lindenhovius, M. Mislove, V. Zamdzhiev, "Semantics for variational quantum programming," Proc. ACM Program. Lang. (POPL 2022) 1-31. doi:10.1145/3498687
- [Lee21] D. Lee, V. Perrelle, B. Valiron, Z. Xu, "Concrete Categorical Model of a Quantum Circuit Description Language with Measurement," FSTTCS 2021. doi:10.4230/LIPIcs.FSTTCS.2021.51
- [Lin25c] B. Lindenhovius, V. Zamdzhiev, "Operator Spaces, Linear Logic and the Heisenberg-Schrödinger Duality of Quantum Theory," LICS 2025, pp. 870-883. doi:10.1109/LICS65433.2025.00071
- [Luo19] X.-Z. Luo, J.-G. Liu, P. Zhang, L. Wang, "Yao.jl: Extensible, Efficient Framework for Quantum Algorithm Design," Quantum 4 (2020) 341. doi:10.22331/q-2020-10-11-341
- [Luo23] J. Luo, J. Zhao, "Formalization of Quantum Intermediate Representations for Code Safety," J. Syst. Softw. (2023) 112236. doi:10.48550/arXiv.2303.14500
- [Nes09] M. Van den Nest, "Simulating quantum computers with probabilistic methods," Quantum Inf. Comput. 11 (2011) 784-812. doi:10.26421/QIC11.9-10-5
- [Pay17] J. Paykin, R. Rand, S. Zdancewic, "QWIRE: a core language for quantum circuits," POPL 2017. doi:10.1145/3009837.3009894
- [Ran21] R. Rand, A. Sundaram, K. Singhal, B. Lackey, "Gottesman Types for Quantum Programs," QPL 2020, EPTCS 340, pp. 279-290. doi:10.4204/EPTCS.340.14
- [Rio17] F. Rios, P. Selinger, "A categorical model for a quantum circuit description language," QPL 2017, EPTCS 266, pp. 164-178. doi:10.4204/EPTCS.266.11
- [Ros26] N. J. Ross, S. Wesley, "Parameterized Quantum Circuit Semantics Through Enriched Categories," arXiv:2607.16114. doi:10.48550/arXiv.2607.16114
- [Sel04] P. Selinger, "Towards a quantum programming language," Math. Struct. Comput. Sci. 14 (2004) 527-586. doi:10.1017/S0960129504004256
- [Sel04b] P. Selinger, B. Valiron, "A lambda calculus for quantum computation with classical control," Math. Struct. Comput. Sci. 16 (2006) 527-552. doi:10.1017/S0960129506005238
- [Sel07] P. Selinger, "Dagger Compact Closed Categories and Completely Positive Maps," ENTCS 170 (2007) 139-163. doi:10.1016/j.entcs.2006.12.018
- [Sin22] K. Singhal, K. Hietala, S. Marshall, R. Rand, "Q# as a Quantum Algorithmic Language," QPL 2022, EPTCS 394, pp. 170-191. doi:10.4204/EPTCS.394.10
- [Sta15] S. Staton, "Algebraic Effects, Linearity, and Quantum Programming Languages," POPL 2015. doi:10.1145/2775051.2676999
- [Voi22] F. Voichick, L. Li, R. Rand, M. Hicks, "Qunity: A Unified Language for Quantum and Classical Computing," Proc. ACM Program. Lang. (POPL 2023) 921-951. doi:10.1145/3571225
- [Woo87] W. K. Wootters, "A Wigner-function formulation of finite-state quantum mechanics," Ann. Phys. 176 (1987) 1-21. doi:10.1016/0003-4916(87)90176-X
