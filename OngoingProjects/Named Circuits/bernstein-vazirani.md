---
Template: Default
Title: Bernstein-Vazirani
Context: Wolfram`QuantumFramework`
---

# Bernstein-Vazirani

Bernstein-Vazirani recovers a hidden bit string *a* in a single query, where the
oracle computes the dot product `f(x) = a · x (mod 2)`. The named circuit builds
the oracle (a `CNOT` from each wire whose secret bit is 1 onto the ancilla), the
surrounding Hadamard layers, and the measurement, whose outcome reads the secret
back directly.

## The algorithm

The secret given as a bit string:

```wl
QuantumCircuitOperator[{"BernsteinVazirani"["101"]}]["Diagram"]
```

---

Running the full circuit, a single measurement returns the secret `101` with certainty:

```wl
QuantumCircuitOperator[{"BernsteinVazirani"["101"]}][]["Probabilities"]
```

---

The secret given as a bit list:

```wl
QuantumCircuitOperator[{"BernsteinVazirani"[{1, 1, 0, 1}]}]["Diagram"]
```

---

A second argument shifts every wire up by one, leaving wire 1 idle with the register on wires 2-5 and the ancilla on 6:

```wl
QuantumCircuitOperator[{"BernsteinVazirani"["1101", 1]}]["Diagram"]
```

## The oracle alone

The oracle for a chosen secret: one `CNOT` onto the ancilla from each wire whose bit is 1:

```wl
QuantumCircuitOperator[{"BernsteinVaziraniOracle"["1011"]}]["Diagram"]
```

---

The same oracle from a bit list:

```wl
QuantumCircuitOperator[{"BernsteinVaziraniOracle"[{1, 0, 1, 1}]}]["Diagram"]
```
