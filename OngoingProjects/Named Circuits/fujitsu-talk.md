---
Template: Default
Title: QuantumFramework Selected Circuits
Context: Wolfram`QuantumFramework`
---

# QuantumFramework Selected Circuits

Selected examples assembled from the named-circuit pages. Every cell is
`QuantumCircuitOperator[{"Name"[...]}]`: a bare string builds the default
instance, `"Name"[...]` passes arguments, and the same object draws a diagram,
returns a matrix, and acts on a state.

## Grover

The search step marking the states that satisfy a boolean formula:

```wl
QuantumCircuitOperator["Grover"[\[FormalX][1] && \[FormalX][2] && ! \[FormalX][3]]]["Diagram"]
```

---

Two phase-Grover iterations on the uniform superposition concentrate the probability on the marked $|111\rangle$:

```wl
Nest[QuantumCircuitOperator["GroverPhase"[\[FormalX][1] && \[FormalX][2] && \[FormalX][3]]], QuantumState["+++"], 2]["Probabilities"]
```

## Fourier transform and phase estimation

The quantum Fourier transform moves a register between the computational and the
Fourier basis; phase estimation uses its inverse to read the eigenphase of a
unitary into a counting register. The Fourier transform and its inverse also
exist as single `QuantumOperator` matrices, so under the circuit form here the
diagram is the gate *decomposition*, Hadamards and controlled phase rotations,
not an atomic box.

The transform on three qubits: a Hadamard and a ladder of controlled phase shifts on each wire, then the wire-reversing swap:

```wl
QuantumCircuitOperator[{"Fourier"[3]}]["Diagram"]
```

---

The eigenphase of the `T` gate read into a three-qubit counting register; the constructor prepares the target in $|1\rangle$, selecting its eigenphase 1/8 (both `T` and `S` have eigenphase 0 on $|0\rangle$), and ends in an inverse Fourier transform and measurement:

```wl
QuantumCircuitOperator[{"PhaseEstimation"[QuantumOperator["T"], 3]}]["Diagram"]
```

---

Running it, the counting register reads `001` = 1/8, the eigenphase of `T`, with certainty:

```wl
QuantumCircuitOperator[{"PhaseEstimation"[QuantumOperator["T"], 3]}][]["Probabilities"]
```

## Hamiltonian simulation: `Trotterization`

A second-order Trotter step for the sum `X + Z`, symmetrically ordered:

```wl
QuantumCircuitOperator[{"Trotterization"[{"X", "Z"}, 2]}]["Diagram"]
```

## Non-unitary elements

Not every circuit element is a unitary gate. A circuit can also measure a wire,
reset it, trace it out, or discard it. Each of these is a named element whose
action on a state is described below, not just drawn as a diagram.

An unrecorded projective measurement on wire 1; applied to $|+\rangle$ it returns the maximally mixed state $I/2$:

```wl
QuantumCircuitOperator[{"Measurement"[1]}]["Diagram"]
```

---

Acting on $|+\rangle$, the density matrix is $I/2$: the coherence is gone:

```wl
Normal[QuantumCircuitOperator[{"Measurement"[1]}][QuantumState["+"]]["DensityMatrix"]]
```

---

A completely positive map applied directly in the circuit (the strength nests inside the name, `"Depolarizing"[p]`); the depolarizing channel at 3/10 partially mixes a pure input:

```wl
QuantumCircuitOperator[{"Channel"["Depolarizing"[3/10]]}]["Diagram"]
```

---

Acting on $|+\rangle$ it damps the coherence rather than erasing it, giving $\begin{pmatrix}1/2 & 7/20\\ 7/20 & 1/2\end{pmatrix}$:

```wl
Normal[QuantumCircuitOperator[{"Channel"["Depolarizing"[3/10]]}][QuantumState["+"]]["DensityMatrix"]]
```

## Exact linear algebra: `Multiplexer`

Its matrix is block-diagonal: $X$ on the control-0 block, $H$ on the control-1 block:

```wl
Normal[QuantumCircuitOperator[{"Multiplexer"["X", "H"]}]["Matrix"]]
```

## Tensor-network atoms

These are the wire-level building blocks of the ZX and tensor-network calculus:
bending a wire into a state or an effect, crossing two wires, copying a basis
state, and the phase spider. Each caption states what the atom's matrix does,
verified in the kernel.

A one-in one-out Z-spider carries a phase on the $|1\rangle$ leg; at phase $\pi$ it is $\mathrm{diag}(1, -1) = Z$:

```wl
QuantumCircuitOperator[{"ZSpider"[Pi]}]["Diagram"]
```

---

The braid crosses two wires: a swap carrying a $-1$ phase on $|11\rangle$:

```wl
QuantumCircuitOperator[{"Braid"[2]}]["Diagram"]
```

---

The cup bends a wire from two outputs, giving the unnormalized Bell state $|00\rangle + |11\rangle$:

```wl
QuantumCircuitOperator[{"Cup"[2]}]["Diagram"]
```

## Qudits and angular momentum

Beyond qubits, the library carries the clock-and-shift operators of a
$d$-level qudit, the qudit entangling gate, and angular-momentum operators whose
dimension is set by the spin. Each caption states the matrix or the action.

The clock shift on a three-level qudit, cycling $|j\rangle \to |j+1 \bmod 3\rangle$:

```wl
QuantumCircuitOperator[{"Shift"[3]}]["Diagram"]
```

---

Acting on $|0\rangle$ it advances the qutrit to $|1\rangle$, amplitudes $(0, 1, 0)$:

```wl
Normal[QuantumCircuitOperator[{"Shift"[3]}][QuantumState[{1, 0, 0}, 3]]["StateVector"]]
```

---

The qudit analogue of CNOT, adding the control to the target modulo 3:

```wl
QuantumCircuitOperator[{"SUM"[3]}]["Diagram"]
```
