---
Template: Default
Title: Non-Unitary Elements
Context: Wolfram`QuantumFramework`
---

# Non-Unitary Elements

Not every circuit element is a unitary gate. A circuit can also measure a wire,
reset it, trace it out, or discard it. Each of these is a named element whose
action on a state is described below, not just drawn as a diagram.

## Measurement: `Measurement`

An unrecorded projective measurement on wire 1; applied to $|+\rangle$ it returns the maximally mixed state $I/2$:

```wl
QuantumCircuitOperator[{"Measurement"[1]}]["Diagram"]
```

---

Acting on $|+\rangle$, the density matrix is $I/2$: the coherence is gone:

```wl
Normal[QuantumCircuitOperator[{"Measurement"[1]}][QuantumState["+"]]["DensityMatrix"]]
```

## Reset: `Reset`

Resets a wire to $|0\rangle$ regardless of its input, so $|1\rangle$ becomes $|0\rangle\langle0|$:

```wl
QuantumCircuitOperator[{"Reset"[]}]["Diagram"]
```

## Partial trace: `Trace`

The partial trace over a wire, keeping the reduced state of the rest; tracing one half of a Bell pair leaves the other maximally mixed:

```wl
QuantumCircuitOperator[{"Trace"[2]}]["Diagram"]
```

## Discard: `Discard`

The effect $\langle0| + \langle1|$ that contracts a wire's amplitudes, not a trace-preserving removal: it sends a Bell pair's kept wire to $|+\rangle$ but maps $|-\rangle$ to zero, so it carries a normalization or postselection caveat:

```wl
QuantumCircuitOperator[{"Discard"[2]}]["Diagram"]
```

## Channel: `Channel`

A completely positive map applied directly in the circuit (the strength nests inside the name, `"Depolarizing"[p]`); the depolarizing channel at 3/10 partially mixes a pure input:

```wl
QuantumCircuitOperator[{"Channel"["Depolarizing"[3/10]]}]["Diagram"]
```

---

Acting on $|+\rangle$ it damps the coherence rather than erasing it, giving $\begin{pmatrix}1/2 & 7/20\\ 7/20 & 1/2\end{pmatrix}$:

```wl
Normal[QuantumCircuitOperator[{"Channel"["Depolarizing"[3/10]]}][QuantumState["+"]]["DensityMatrix"]]
```

## Generators and the diagonal map

Two of these names are the one open-system generator in its two conventions.
`Hamiltonian[H, {jumps}]` is the form `QuantumEvolve` consumes: the commutator
$[H, \cdot]$ for a closed system, and $i\mathcal L$, with $\mathcal L$ the
Lindblad generator, once jump operators are given. `Liouvillian[H, {jumps}]` is
$\mathcal L$ itself, the bare generator of $\dot\rho = \mathcal L[\rho]$, equal
to $-i\,$`Hamiltonian[H, {jumps}]`; reach for it when $\mathcal L$ is wanted (its
spectrum, say) rather than for `QuantumEvolve`, which exponentiates $-i$ times
its argument. Evolving with `Hamiltonian[H, {jumps}]` and the separate
Hamiltonian-and-jump-operator form `QuantumEvolve[H, {jumps}, ...]` give the same
result. A third, `Decohere`, is not a generator at all but a
dimension-reducing *representation* map: the $2\times4$ operator
$\{\{1,0,0,0\},\{0,0,0,1\}\}$ that pulls the diagonal from a doubled
representation (4 amplitudes in, 2 out). Because its input and output dimensions
differ it cannot serve as an evolution generator, and it is not a
trace-preserving dephasing channel.
