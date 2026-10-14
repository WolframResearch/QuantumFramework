---
Template: Default
Title: Fourier Transform and Phase Estimation
Context: Wolfram`QuantumFramework`
---

# Fourier Transform and Phase Estimation

The quantum Fourier transform moves a register between the computational and the
Fourier basis; phase estimation uses its inverse to read the eigenphase of a
unitary into a counting register. The Fourier transform and its inverse also
exist as single `QuantumOperator` matrices, so under the circuit form here the
diagram is the gate *decomposition*, Hadamards and controlled phase rotations,
not an atomic box.

## Quantum Fourier transform

The transform on three qubits: a Hadamard and a ladder of controlled phase shifts on each wire, then the wire-reversing swap:

```wl
QuantumCircuitOperator[{"Fourier"[3]}]["Diagram"]
```

---

A second argument offsets the transform up by that many leading wires:

```wl
QuantumCircuitOperator[{"Fourier"[3, 1]}]["Diagram"]
```

---

The inverse transform, the dagger of the forward one:

```wl
QuantumCircuitOperator[{"InverseFourier"[3]}]["Diagram"]
```

## Phase estimation

The eigenphase of the `T` gate read into a three-qubit counting register; the constructor prepares the target in $|1\rangle$, selecting its eigenphase 1/8 (both `T` and `S` have eigenphase 0 on $|0\rangle$), and ends in an inverse Fourier transform and measurement:

```wl
QuantumCircuitOperator[{"PhaseEstimation"[QuantumOperator["T"], 3]}]["Diagram"]
```

---

Running it, the counting register reads `001` = 1/8, the eigenphase of `T`, with certainty:

```wl
QuantumCircuitOperator[{"PhaseEstimation"[QuantumOperator["T"], 3]}][]["Probabilities"]
```

---

The same for the `S` gate (phase 1/4), a different target unitary:

```wl
QuantumCircuitOperator[{"PhaseEstimation"[QuantumOperator["S"], 3]}]["Diagram"]
```

---

A larger counting register, five qubits: more resolution for a general phase, though the dyadic 1/8 here is already exact with three:

```wl
QuantumCircuitOperator[{"PhaseEstimation"[QuantumOperator["T"], 5]}]["Diagram"]
```

---

The `PowerExpand` option replaces each controlled `U^(2^k)` by `2^k` repeated controlled-`U` gates:

```wl
QuantumCircuitOperator[{"PhaseEstimation"[QuantumOperator["T"], 3, "PowerExpand" -> True]}]["Diagram"]
```
