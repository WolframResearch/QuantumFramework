---
Template: Default
Title: Gate Gadgets and Encoders
Context: Wolfram`QuantumFramework`
---

# Gate Gadgets and Encoders

These are the small fixed circuits that recur as building blocks: entangled-state
preparations, the three-qubit permutation gates, integer encoders, graph states,
and the general state-preparation compiler.

## Bell pair: `Bell`

A Hadamard then a CNOT, preparing the Bell state $(|00\rangle + |11\rangle)/\sqrt2$:

```wl
QuantumCircuitOperator[{"Bell"}]["Diagram"]
```

---

Acting on $|00\rangle$ it produces the Bell state, amplitudes $(1, 0, 0, 1)/\sqrt2$:

```wl
Normal[QuantumCircuitOperator[{"Bell"}][QuantumState["00"]]["StateVector"]]
```

## GHZ state: `GHZ`

The n-qubit GHZ state $(|0\cdots0\rangle + |1\cdots1\rangle)/\sqrt2$, here on four wires:

```wl
QuantumCircuitOperator[{"GHZ"[4]}]["Diagram"]
```

## Doubly-controlled NOT: `Toffoli`

The doubly-controlled NOT, shown as its 15-gate Clifford+T decomposition:

```wl
QuantumCircuitOperator[{"Toffoli"}]["Diagram"]
```

## Controlled-SWAP: `Fredkin`

The controlled-SWAP, shown as its elementary-gate decomposition:

```wl
QuantumCircuitOperator[{"Fredkin"}]["Diagram"]
```

## Magic basis: `Magic`

The two-qubit magic-basis circuit, `S` on each wire then `H` and `CNOT`:

```wl
QuantumCircuitOperator[{"Magic"}]["Diagram"]
```

## Integer encoding: `Number`

Encodes an integer in the computational basis, here 5 as the bit string `101` with an `X` on each 1:

```wl
QuantumCircuitOperator[{"Number"[5]}]["Diagram"]
```

---

A second argument pads the encoding to a fixed number of qubits:

```wl
QuantumCircuitOperator[{"Number"[5, 4]}]["Diagram"]
```

## Fourier-basis integer: `PhaseNumber`

Encodes an integer in the Fourier basis, Hadamards followed by phase shifts rather than bit flips:

```wl
QuantumCircuitOperator[{"PhaseNumber"[5]}]["Diagram"]
```

## Graph state: `Graph`

The graph state of a given graph: a Hadamard on every vertex and a controlled-Z on every edge, here a four-vertex path:

```wl
QuantumCircuitOperator[{"Graph"[PathGraph[Range[4]]]}]["Diagram"]
```

## State preparation: `QuantumState`

Compiles a circuit that prepares a given pure state (mixed inputs are rejected), here the Bell state:

```wl
QuantumCircuitOperator[{"QuantumState"[QuantumState["Bell"]]}]["Diagram"]
```

---

The same constructor under its other name, `StatePreparation`, here preparing a GHZ state:

```wl
QuantumCircuitOperator[{"StatePreparation"[QuantumState["GHZ"]]}]["Diagram"]
```
