---
Template: Default
Title: Deutsch and Deutsch-Jozsa
Context: Wolfram`QuantumFramework`
---

# Deutsch and Deutsch-Jozsa

Deutsch-Jozsa decides, for a function *promised* to be either constant or
balanced, which of the two it is in a single query. The named circuit assembles
the oracle, the surrounding Hadamard layers, and the final measurement. The
function argument does not enforce the promise: a formula or an integer indexes
a general truth table, and `Automatic` draws from all Boolean functions (for
three inputs only 72 of 256 are constant or balanced), so a random draw can be a
valid circuit on an invalid problem instance. The explicit examples below pick
balanced or constant functions; the `Automatic` example does not enforce the
promise. Each algorithm comes in a *boolean* form (bit-flip oracle plus a
$|-\rangle$ ancilla) and a *phase* form (phase oracle, no ancilla); both read the
input register out in the X basis.

## Deutsch: one qubit

The balanced case, with the oracle `f(x) = x`:

```wl
QuantumCircuitOperator[{"Deutsch"[\[FormalX][1]]}]["Diagram"]
```

---

The constant case, selected by the integer index (here `f ≡ 0`):

```wl
QuantumCircuitOperator[{"Deutsch"[1]}]["Diagram"]
```

---

The phase form, which drops the ancilla and measures in the X basis (width 2):

```wl
QuantumCircuitOperator[{"DeutschPhase"[\[FormalX][1]]}]["Diagram"]
```

## Deutsch-Jozsa: n qubits

A balanced function on three inputs, here the parity `x1 ⊻ x2 ⊻ x3`:

```wl
QuantumCircuitOperator[{"DeutschJozsa"[Xor[\[FormalX][1], \[FormalX][2], \[FormalX][3]]]}]["Diagram"]
```

---

A constant function on three inputs, selected by the integer index:

```wl
QuantumCircuitOperator[{"DeutschJozsa"[1, 3]}]["Diagram"]
```

---

The decision, read out: the parity function sends all the weight onto the all-minus outcome, so the all-plus string that would certify a constant function has probability zero (the constant case puts probability 1 there instead):

```wl
QuantumCircuitOperator[{"DeutschJozsa"[Xor[\[FormalX][1], \[FormalX][2], \[FormalX][3]]]}][]["Probabilities"]
```

---

A second argument fixes the number of input qubits when the function is left to a random draw:

```wl
QuantumCircuitOperator[{"DeutschJozsa"[Automatic, 3]}]["Diagram"]
```

---

The phase form on two inputs, here the balanced equivalence `¬(x1 ⊻ x2)`: phase oracle, no ancilla:

```wl
QuantumCircuitOperator[{"DeutschJozsaPhase"[! Xor[\[FormalX][1], \[FormalX][2]]]}]["Diagram"]
```

## The oracles alone

The bit-flip oracle for a balanced function, flipping the ancilla where the function is true:

```wl
QuantumCircuitOperator[{"DeutschJozsaBooleanOracle"[Xor[\[FormalX][1], \[FormalX][2]]]}]["Diagram"]
```

---

The phase oracle for the same function, marking by sign with no ancilla:

```wl
QuantumCircuitOperator[{"DeutschJozsaPhaseOracle"[Xor[\[FormalX][1], \[FormalX][2]]]}]["Diagram"]
```

---

The constant oracle on three inputs, selected by the integer index:

```wl
QuantumCircuitOperator[{"DeutschJozsaBooleanOracle"[1, 3]}]["Diagram"]
```
