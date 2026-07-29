## Part 10. Elementary multi-qubit gates and circuits

### 10.1 [BSc] How do I build the standard multi-qubit gates (CNOT, CZ, SWAP, Toffoli, controlled-$U$)?

The standard entangling gates are controlled operations: CNOT and CZ apply $X$ or $Z$ to a target when the
control is $|1\rangle$, SWAP exchanges two qubits, Toffoli is a doubly-controlled $X$, and the general
controlled-$U$ is $|0\rangle\langle0|\otimes I + |1\rangle\langle1|\otimes U$. Every controlled gate has
this projector form, which is its definition.

**WL** : each controlled gate is a sum of projector-tensor-target terms. Here is the general
controlled-$U$ with an arbitrary single-qubit $U$, the pattern the controlled gates follow.

```wl
MatrixForm@With[{p0 = {{1, 0}, {0, 0}}, p1 = {{0, 0}, {0, 1}}, u = Array[Subscript[\[FormalU], ##] &, {2, 2}]}, KroneckerProduct[p0, IdentityMatrix[2]] + KroneckerProduct[p1, u]]
```

The control leaves the $|0\rangle$ subspace alone and applies $U$ on the $|1\rangle$ subspace, so the
matrix is block diagonal, $\mathrm{diag}(I, U)$. One can therefore also assemble it straight from its two
diagonal blocks.

```wl
BlockDiagonalMatrix[{IdentityMatrix[2], Array[Subscript[\[FormalU], ##] &, {2, 2}]}]
```

This returns a structured array rather than an explicit list of rows, so display it to see the entries.

```wl
% // MatrixForm
```

The named gates are this pattern with $U = X$ (CNOT) and $U = Z$ (CZ). Toffoli is the doubly-controlled
$X$: with $P_{11} = |11\rangle\langle11|$ the projector onto both controls being set, it is
$\mathrm{TOF} = (I_4 - P_{11})\otimes I + P_{11}\otimes X$, the controlled-$U$ pattern carried over to two
controls, equivalently $\mathrm{TOF} = I_8 + P_{11}\otimes(X - I)$, the identity plus a correction acting
only on the $|11\rangle$ control subspace. SWAP is not a controlled gate, so it falls outside the projector
pattern. It is the exchange operator
$\tfrac12\sum_{k=0}^{3}\sigma_k\otimes\sigma_k = \tfrac12\big(I\otimes I + X\otimes X + Y\otimes Y + Z\otimes Z\big)$
(with $\sigma_0 = I$), equivalently the permutation that swaps the two qubits.

All four are nevertheless block diagonal, because each fixes some basis states outright and acts
non-trivially only inside the complementary subspace, so each can be assembled from its invariant blocks
alone: $\mathrm{CNOT} = \mathrm{diag}(I_2, X)$, $\mathrm{CZ} = \mathrm{diag}(I_2, Z)$,
$\mathrm{TOF} = \mathrm{diag}(I_6, X)$, and $\mathrm{SWAP} = \mathrm{diag}(1, X, 1)$, the last because SWAP
fixes $|00\rangle$ and $|11\rangle$ and exchanges $|01\rangle \leftrightarrow |10\rangle$.

```wl
MatrixForm /@ <|"CNOT" -> BlockDiagonalMatrix[{IdentityMatrix[2], PauliMatrix[1]}], "CZ" -> BlockDiagonalMatrix[{IdentityMatrix[2], PauliMatrix[3]}], "SWAP" -> BlockDiagonalMatrix[{IdentityMatrix[1], PauliMatrix[1], IdentityMatrix[1]}], "Toffoli" -> BlockDiagonalMatrix[{IdentityMatrix[6], PauliMatrix[1]}]|>
```

These block forms assume the standard wire ordering, controls before target; permuting the wires conjugates
the matrix by a permutation and destroys the pattern.

**QF** : these are named `QuantumOperator`s; read their matrices off the objects.

```wl
AssociationMap[MatrixForm[QuantumOperator[#]["Matrix"]] &, {"CNOT", "CZ", "SWAP", "Toffoli"}]
```

Both give the same gates: CNOT and CZ flip $X$ or phase $Z$ on the target when the control is set, SWAP
permutes the two qubits, and Toffoli flips its target only when both controls are $|1\rangle$. The QF names
carry the object's structure (arity, control and target wires, unitarity), so they drop straight into a
circuit; the projector form on the WL side shows *why* each is the gate it is.

The main advantage of the QF form is that you decide *where* the operator acts by setting its target wires
as a second argument, something a bare matrix cannot carry. CNOT with qubit 2 as the control and qubit 4 as
the target:

```wl
MatrixForm[QuantumOperator["CNOT", {2, 4}]["Matrix"]]
```

The matrix that comes back is still $4\times4$: the wires $\{2,4\}$ sit in the operator's order, not in its
entries. A bare matrix has nowhere to keep that, so in plain WL the placement has to be built by hand,
padding every spectator wire with an `IdentityMatrix` factor. On a four-qubit register that means identity
on wires 1 and 3, the control projectors on wire 2, and the target on wire 4:

```wl
MatrixForm@With[{p0 = {{1, 0}, {0, 0}}, p1 = {{0, 0}, {0, 1}}, id = IdentityMatrix[2]}, KroneckerProduct[id, p0, id, id] + KroneckerProduct[id, p1, id, PauliMatrix[1]]]
```

That one is $16\times16$: writing the placement down forces a commitment to the width of the register, and
adding a fifth qubit means another `IdentityMatrix` factor and a $32\times32$ matrix. The QF operator stays
two-qubit and carries $\{2,4\}$ as structure, so the same object drops into a circuit of any width.

Or SWAP between qubit 1 and qubit 3:

```wl
MatrixForm[QuantumOperator["SWAP", {1, 3}]["Matrix"]]
```

The matrix printed is the bare two-wire gate, not padded with identities on the wires in between; the wire
labels live in the object's `"Order"`. So when this CNOT meets a register it flips qubit 4 exactly when
qubit 2 is set (and the order within the pair fixes which wire is the control), padding the intervening
wires with the identity automatically. Placing the same gate in native WL means tensoring those identities
in by hand and reordering for the non-adjacent wires.

### 10.2a [BSc] How do I compose a circuit from gates and apply it to an input state?

A circuit is an ordered list of gates on named wires; applying it evolves an input state through them in
order. Compose the two-gate Bell circuit, $H$ on qubit 1 then CNOT from 1 to 2, and apply it to
$|00\rangle$.

**WL** : the circuit's action is the ordered matrix product acting on the input vector, CNOT after
$H\otimes I$.

```wl
With[{cnot = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}}}, cnot . KroneckerProduct[HadamardMatrix[2], IdentityMatrix[2]] . {1, 0, 0, 0}]
```

**QF** : a `QuantumCircuitOperator` built from the gate list applies to the `QuantumState` directly.

```wl
QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}][QuantumState["00"]]
```

Both produce the Bell state $|\Phi^+\rangle = \tfrac1{\sqrt2}(|00\rangle + |11\rangle)$: the Hadamard makes
a superposition on the control, and the CNOT correlates the target with it, so a product input comes out
entangled. The QF circuit is a reusable object that also knows its wires and its overall unitary (10.2b).

### 10.2b [BSc] How do I read off the overall unitary of a composed circuit?

The whole circuit is itself one unitary, the ordered product of its gate matrices; reading it off lets you
treat a subcircuit as a single black-box gate. Take the same Bell circuit.

**WL** : multiply the gate matrices in circuit order (with each embedded on its wires), CNOT after
$H\otimes I$.

```wl
With[{cnot = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}}}, cnot . KroneckerProduct[HadamardMatrix[2], IdentityMatrix[2]]]
```

**QF** : the whole circuit collapses to one operator whose overall unitary you read off directly with
`"Matrix"`.

```wl
Normal[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]["Matrix"]]
```

Both give the same $4\times4$ unitary: multiplying the embedded gate matrices in circuit order and reading
`"Matrix"` off the composed operator are two routes to the one matrix that represents the whole circuit.

### 10.3 [BSc] How do I compute the depth and the parallel layers of a circuit?

The depth of a circuit is the number of layers when gates are packed as early as possible: two gates share
a layer only if they act on disjoint qubits, so depth measures the longest chain of gates that must run in
sequence, the circuit's runtime on parallel hardware. Take $H$ on qubits 1 and 2, then CNOT from 1 to 2,
CNOT from 2 to 3, and finally $S$ on qubit 1 (a bare gate name defaults to the first wire): five gates on
three qubits. Here is that circuit, drawn one column per layer, so the packing the rest of the answer
computes is visible from the start.

```wl
QuantumCircuitOperator[{"H" -> 1, "H" -> 2, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}, "S"}]["Diagram"]
```

**WL** : schedule each gate into the earliest layer its qubits allow (one past the latest layer any of its
qubits is currently busy), tracking a per-qubit clock with one entry per wire; the depth is the maximum
over qubits.

```wl
With[{gates = {{1}, {2}, {1, 2}, {2, 3}, {1}}}, Max[Fold[Function[{clock, g}, ReplacePart[clock, Thread[g -> Max[clock[[g]]] + 1]]], ConstantArray[0, Max[gates]], gates]]]
```

**QF** : a `QuantumCircuitOperator` reports its `"Depth"` and its `"Layers"` (the list of parallel slices)
directly.

```wl
With[{qc = QuantumCircuitOperator[{"H" -> 1, "H" -> 2, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}, "S"}]}, {qc["Depth"], Length[qc["Layers"]]}]
```

Both give depth $3$. The two Hadamards share layer 1 (disjoint qubits); the first CNOT waits for both of
them in layer 2; then in layer 3 the second CNOT (on qubits 2 and 3) runs alongside the $S$ on qubit 1,
since those touch disjoint wires. So five gates pack into three layers: it is the depth, not the gate count
of $5$, that sets how long the circuit takes on hardware that fires disjoint gates simultaneously, so
scheduling into layers is the first step in estimating a circuit's real cost.
