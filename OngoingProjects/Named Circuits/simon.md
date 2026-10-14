---
Template: Default
Title: Simon's Algorithm
Context: Wolfram`QuantumFramework`
---

# Simon's Algorithm

Simon's algorithm finds a hidden nonzero mask *s* of a function promised to be
two-to-one with that period: `f(x) = f(x')` exactly when `x' = x` or `x' = x ⊕ s`.
(The bare period `f(x) = f(x ⊕ s)` alone is weaker: a constant function satisfies
it for every *s*.) It uses exponentially fewer queries than any classical method.
The named circuit is the *sampling* step: oracle, Hadamard layers on the input
register, and a measurement that returns one string `y` with `y · s = 0`;
recovering *s* needs repeated samples and classical linear algebra on top. The
mask is given as a bit string, a bit list, or left to a random draw. One caveat
verified in the kernel: the generated function respects the period `f(x) = f(x ⊕ s)`
but is not guaranteed strictly two-to-one (the outputs are drawn at random), so
it realizes the period rather than the full Simon promise.

## Sampling circuit: `Simon`

The mask given as a bit string, on a 3-wire input register feeding a 3-wire output register:

```wl
QuantumCircuitOperator[{"Simon"["101"]}]["Diagram"]
```

---

The mask given as a bit list:

```wl
QuantumCircuitOperator[{"Simon"[{1, 0}]}]["Diagram"]
```

---

A random mask, defaulting to three bits (so the register is fixed at three input and three output wires):

```wl
QuantumCircuitOperator[{"Simon"[]}]["Diagram"]
```

## Oracle alone: `SimonOracle`

The two-register oracle for a chosen mask, mapping a 3-wire input register into a 3-wire output register:

```wl
QuantumCircuitOperator[{"SimonOracle"["110"]}]["Diagram"]
```

---

The same mask given as a bit list:

```wl
QuantumCircuitOperator[{"SimonOracle"[{1, 1, 0}]}]["Diagram"]
```
