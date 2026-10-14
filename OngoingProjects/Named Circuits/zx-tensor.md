---
Template: Default
Title: Tensor-Network Atoms
Context: Wolfram`QuantumFramework`
---

# Tensor-Network Atoms

These are the wire-level building blocks of the ZX and tensor-network calculus:
bending a wire into a state or an effect, crossing two wires, copying a basis
state, and the phase spider. Each caption states what the atom's matrix does,
verified in the kernel.

## Cup: `Cup`

The cup bends a wire from two outputs, giving the unnormalized Bell state $|00\rangle + |11\rangle$:

```wl
QuantumCircuitOperator[{"Cup"[2]}]["Diagram"]
```

---

It has no inputs; its output state has amplitudes $(1, 0, 0, 1)$:

```wl
Normal[QuantumCircuitOperator[{"Cup"[2]}][]["StateVector"]]
```

## Cap: `Cap`

The cap is the matching effect $\langle00| + \langle11|$, a wire bent into two inputs that contracts a pair:

```wl
QuantumCircuitOperator[{"Cap"[2]}]["Diagram"]
```

## Braid: `Braid`

The braid crosses two wires: a swap carrying a $-1$ phase on $|11\rangle$:

```wl
QuantumCircuitOperator[{"Braid"[2]}]["Diagram"]
```

## Copy: `Copy`

The copy spider clones a basis state onto two wires, $|0\rangle \mapsto |00\rangle$ and $|1\rangle \mapsto |11\rangle$:

```wl
QuantumCircuitOperator[{"Copy"[2]}]["Diagram"]
```

## Phase spider: `ZSpider`

A one-in one-out Z-spider carries a phase on the $|1\rangle$ leg; at phase $\pi$ it is $\mathrm{diag}(1, -1) = Z$:

```wl
QuantumCircuitOperator[{"ZSpider"[Pi]}]["Diagram"]
```

---

At phase $\pi/2$ the same spider is $\mathrm{diag}(1, i) = S$:

```wl
QuantumCircuitOperator[{"ZSpider"[Pi/2]}]["Diagram"]
```
