---
Template: Default
Title: Composition and Simulation
Context: Wolfram`QuantumFramework`
---

# Composition and Simulation

These constructions take operators or whole circuits as input and assemble a
larger circuit: controlling a circuit on extra wires, multiplexing a choice of
gate on a control register, routing two unitaries in indefinite order, solving a
linear system, and compiling a Hamiltonian into a product-formula circuit.

## Controlling a circuit: `Controlled`

A whole circuit controlled on an extra wire, here a two-qubit Fourier transform on wires 2,3 fired by wire 1 (the transform is offset to 2,3 so it does not collide with the control):

```wl
QuantumCircuitOperator[{"Controlled"[QuantumCircuitOperator["Fourier"[2, 1]], {1}]}]["Diagram"]
```

---

Two control wires (3 and 4) gating a GHZ preparation on wires 1,2, so the pair is prepared only when both controls fire:

```wl
QuantumCircuitOperator[{"Controlled"[QuantumCircuitOperator["GHZ"[2]], {3, 4}]}]["Diagram"]
```

## Multiplexing a gate: `Multiplexer`

A gate chosen by the control register, one target gate per control value, here three gates on a two-wire control (shown as its gate decomposition):

```wl
QuantumCircuitOperator[{"Multiplexer"["X", "Y", "Z"]}]["Diagram"]
```

---

Two gates on a single control wire:

```wl
QuantumCircuitOperator[{"Multiplexer"["X", "H"]}]["Diagram"]
```

---

Its matrix is block-diagonal: $X$ on the control-0 block, $H$ on the control-1 block:

```wl
Normal[QuantumCircuitOperator[{"Multiplexer"["X", "H"]}]["Matrix"]]
```

## Indefinite causal order: `Switch`

The quantum switch of two unitaries, applying them in an order set by an ancilla; preparing that ancilla in superposition puts the order in superposition (shown as its gate decomposition):

```wl
QuantumCircuitOperator[{"Switch"[QuantumOperator["X"], QuantumOperator["Z"]]}]["Diagram"]
```

## Linear systems: `ControlledMultiplexer`

The variational linear-solver construction for `A x = b` (coefficient-state preparation, an ansatz, multiplexing, and an ancilla contraction, returned before the fidelity cost and optimization), here a diagonal `A` and a two-component `b`:

```wl
QuantumCircuitOperator[{"ControlledMultiplexer"[{{1, 0}, {0, 2}}, {1, 1}]}]["Diagram"]
```

## Hamiltonian simulation: `Trotterization`

A second-order Trotter step for the sum `X + Z`, symmetrically ordered:

```wl
QuantumCircuitOperator[{"Trotterization"[{"X", "Z"}, 2]}]["Diagram"]
```

---

A first-order step for `X + Y + Z` repeated three times, which shrinks the step error:

```wl
QuantumCircuitOperator[{"Trotterization"[{"X", "Y", "Z"}, 1, 3]}]["Diagram"]
```
