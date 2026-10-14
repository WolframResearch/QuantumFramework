---
Template: Default
Title: Primitive Gates
Context: Wolfram`QuantumFramework`
---

# Primitive Gates

The atomic gates every circuit is built from. These are commodity, so this page
is a compact reference: one cell per distinct gate, with trivial name aliases
noted rather than repeated.

## Pauli gates: `X`, `Y`, `Z`

The bit flip $X$, swapping $|0\rangle \leftrightarrow |1\rangle$:

```wl
QuantumCircuitOperator[{"X"}]["Diagram"]
```

---

The $Y$ gate, a combined bit-and-phase flip sending $|0\rangle \to i|1\rangle$:

```wl
QuantumCircuitOperator[{"Y"}]["Diagram"]
```

---

The sign flip $Z$, negating $|1\rangle$:

```wl
QuantumCircuitOperator[{"Z"}]["Diagram"]
```

The aliases `PauliX`, `PauliY`, `PauliZ` build the same gates, as does `NOT` for `X`.

## Hadamard: `H`

The Hadamard, mapping $|0\rangle \to |+\rangle$ and $|1\rangle \to |-\rangle$ (alias `Hadamard`):

```wl
QuantumCircuitOperator[{"H"}]["Diagram"]
```

---

Acting on $|0\rangle$ it gives the equal superposition $(|0\rangle + |1\rangle)/\sqrt2$:

```wl
Normal[QuantumCircuitOperator[{"H"}][QuantumState["0"]]["StateVector"]]
```

## Phase gates: `S`, `T`

The `S` gate, a $\pi/2$ phase on $|1\rangle$, $\mathrm{diag}(1, i)$:

```wl
QuantumCircuitOperator[{"S"}]["Diagram"]
```

---

The `T` gate, a $\pi/4$ phase on $|1\rangle$, $\mathrm{diag}(1, e^{i\pi/4})$:

```wl
QuantumCircuitOperator[{"T"}]["Diagram"]
```

## Square root of NOT: `V`

The square root of `NOT`, with $V^2 = X$ (note: bare `SX` is the two-qubit $S \otimes X$, while `SX[]` resolves to the single-qubit root):

```wl
QuantumCircuitOperator[{"V"}]["Diagram"]
```

## Rotations: `RX`, `RY`, `RZ`

Rotation about the X axis by an angle:

```wl
QuantumCircuitOperator[{"RX"[Pi/2]}]["Diagram"]
```

---

Rotation about the Y axis:

```wl
QuantumCircuitOperator[{"RY"[Pi/3]}]["Diagram"]
```

---

Rotation about the Z axis:

```wl
QuantumCircuitOperator[{"RZ"[Pi/4]}]["Diagram"]
```

## General single-qubit gate: `U`

The general single-qubit gate from three Euler angles; here $U[\pi/2, 0, \pi]$ is the Hadamard:

```wl
QuantumCircuitOperator[{"U"[Pi/2, 0, Pi]}]["Diagram"]
```

## Phase gate: `P`

A phase gate $\mathrm{diag}(1, e^{i\theta})$ by an explicit angle (the integer form `PhaseShift[k]` gives a $2\pi/2^k$ phase for positive $k$, and $2\pi\,\mathrm{Sign}(k)/2^{|k|}$ in general):

```wl
QuantumCircuitOperator[{"P"[Pi/4]}]["Diagram"]
```

## Two-qubit gates: `CNOT`, `CZ`, `SWAP`

Controlled NOT (alias `CX`):

```wl
QuantumCircuitOperator[{"CNOT"}]["Diagram"]
```

---

Controlled Z:

```wl
QuantumCircuitOperator[{"CZ"}]["Diagram"]
```

---

The swap of two wires:

```wl
QuantumCircuitOperator[{"SWAP"}]["Diagram"]
```
