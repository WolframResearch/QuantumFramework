---
Template: Default
---

# Quantum in Finite Dimensions: Answers in WL and in QuantumFramework (Part 10)

A companion answer key to `Question-List.md`. For each question it gives **two** worked answers:

1. **WL** : native Wolfram Language only (plain vectors and matrices, `PauliMatrix`, `KroneckerProduct`,
   `HadamardMatrix`, and friends). Nothing beyond the built-in language.
2. **QF** : the same task done through the **QuantumFramework** paclet, leaning on its objects and
   property downvalues (`QuantumCircuitOperator[...]["..."]`) rather than rebuilding matrices by hand. A
   QF result is a rich object, a gate or a circuit, not a bare array: it can act on a target, compose with
   others, report its overall unitary, its depth and its layer structure, and more. The answers keep
   results as objects and treat a matrix as just one of the object's many representations, not the thing
   itself.

Two habits run through the answers. They write the computation **explicitly** rather than hiding it in a
helper function, because the code is itself part of the explanation: a controlled gate is built from its
definition $|0\rangle\langle0|\otimes I + |1\rangle\langle1|\otimes U$, not pasted in as a fixed matrix.
And they prefer **built-in features** (a named gate, a circuit's own `"Depth"`) over hand-rolled
constructions. Each cell does **one** thing and is preceded by a sentence saying what. Run the WL answers
in a bare kernel; the QF answers need the paclet loaded once.

## Setup

The QF answers load the paclet once:

```wl
Needs["Wolfram`QuantumFramework`"]
```

---

## Part 10. Elementary multi-qubit gates and circuits

### 10.1 [BSc] How do I build the standard multi-qubit gates (CNOT, CZ, SWAP, Toffoli, controlled-$U$)?

The standard entangling gates are controlled operations: CNOT and CZ apply $X$ or $Z$ to a target when the
control is $|1\rangle$, SWAP exchanges two qubits, Toffoli is a doubly-controlled $X$, and the general
controlled-$U$ is $|0\rangle\langle0|\otimes I + |1\rangle\langle1|\otimes U$. Every controlled gate has
this projector form, which is its definition.

**WL** : each controlled gate is a sum of projector-tensor-target terms; SWAP is the permutation matrix.
Here is the general controlled-$U$ with an arbitrary single-qubit $U$, the pattern all the others follow.

```wl
With[{p0 = {{1, 0}, {0, 0}}, p1 = {{0, 0}, {0, 1}}}, KroneckerProduct[p0, IdentityMatrix[2]] + KroneckerProduct[p1, u]]
```

The named gates are this pattern with $U = X, Z$, plus the two permutation gates SWAP and Toffoli.

```wl
With[{p0 = {{1, 0}, {0, 0}}, p1 = {{0, 0}, {0, 1}}, p11 = DiagonalMatrix[{0, 0, 0, 1}]}, <|"CNOT" -> KroneckerProduct[p0, IdentityMatrix[2]] + KroneckerProduct[p1, PauliMatrix[1]], "CZ" -> KroneckerProduct[p0, IdentityMatrix[2]] + KroneckerProduct[p1, PauliMatrix[3]], "SWAP" -> {{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}}, "Toffoli" -> IdentityMatrix[8] + KroneckerProduct[p11, PauliMatrix[1] - IdentityMatrix[2]]|>]
```

**QF** : these are named `QuantumOperator`s; read their matrices off the objects.

```wl
AssociationMap[Normal[QuantumOperator[#]["Matrix"]] &, {"CNOT", "CZ", "SWAP", "Toffoli"}]
```

Both give the same gates: CNOT and CZ flip $X$ or phase $Z$ on the target when the control is set, SWAP
permutes the two qubits, and Toffoli flips its target only when both controls are $|1\rangle$. The QF names
carry the object's structure (arity, control and target wires, unitarity), so they drop straight into a
circuit; the projector form on the WL side shows *why* each is the gate it is.

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

**QF** : the circuit collapses to a single `QuantumOperator` via `"QuantumOperator"`; here is its matrix.

```wl
Normal[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]["QuantumOperator"]["Matrix"]]
```

Both give the same $4\times4$ Bell-preparation unitary, whose first column is $|\Phi^+\rangle$ (its action
on $|00\rangle$, as in 10.2a). Compressing a circuit to its overall unitary is how a multi-gate block is
reused as a named subroutine, at the cost of hiding the internal parallelism that 10.3 exposes.

### 10.3 [BSc] How do I compute the depth and the parallel layers of a circuit?

The depth of a circuit is the number of layers when gates are packed as early as possible: two gates share
a layer only if they act on disjoint qubits, so depth measures the longest chain of gates that must run in
sequence, the circuit's runtime on parallel hardware. Take $H$ on qubit 1, $H$ on qubit 2, CNOT from 1 to
2, then $X$ on qubit 1.

**WL** : schedule each gate into the earliest layer its qubits allow (one past the latest layer any of its
qubits is currently used), tracking a per-qubit clock; the depth is the maximum over qubits.

```wl
With[{gates = {{1}, {2}, {1, 2}, {1}}}, Max[Fold[Function[{clock, g}, ReplacePart[clock, Thread[g -> Max[clock[[g]]] + 1]]], {0, 0}, gates]]]
```

**QF** : a `QuantumCircuitOperator` reports its `"Depth"` and its `"Layers"` (the list of parallel slices)
directly.

```wl
With[{qc = QuantumCircuitOperator[{"H" -> 1, "H" -> 2, "CNOT" -> {1, 2}, "X" -> 1}]}, {qc["Depth"], Length[qc["Layers"]]}]
```

Both give depth $3$: the two Hadamards run together in layer 1 (disjoint qubits), the CNOT waits for both
in layer 2, and the final $X$ waits for qubit 1 in layer 3. Depth, not gate count, sets how long a circuit
takes on hardware that can fire disjoint gates simultaneously, so scheduling into layers is the first step
in estimating a circuit's real cost.
