---
Template: Default
Title: Named Quantum Circuits
Context: Wolfram`QuantumFramework`
---

# Named Quantum Circuits

`QuantumCircuitOperator` accepts a library of named constructions. A bare string
builds the default instance; `"Name"[...]` passes arguments. This page walks the
**Grover** family through its main signatures, one line each, showing the
circuit diagram.

The names split into two layers: the diffusion block alone, and the full search
step $G = D\,O$ that wraps a diffusion block around an oracle.

## The diffusion block

The diffusion reflection on *n* wires, $I - 2|{+}\cdots{+}{-}\rangle\langle{+}\cdots{+}{-}|$: every input wire in $|+\rangle$ and the flag wire in $|-\rangle$:

```wl
QuantumCircuitOperator["GroverDiffusion"[4]]["Diagram"]
```

---

The block placed on a chosen (possibly non-contiguous) set of wires, marking onto the highest:

```wl
QuantumCircuitOperator["GroverDiffusion"[{1, 3}]]["Diagram"]
```

---

A supplied operator sets the marking gate: its wire (2) becomes the target, the other wires its controls:

```wl
QuantumCircuitOperator["GroverDiffusion"[3, QuantumOperator["Z", {2}]]]["Diagram"]
```

---

The same reflection built from zero-controls instead of an X-sandwich: an identical operator in a shallower circuit (depth 3 rather than 5):

```wl
QuantumCircuitOperator["GroverDiffusion0"[3]]["Diagram"]
```

---

The phase-oracle diffusion, $I - 2|{+}\cdots{+}\rangle\langle{+}\cdots{+}|$: a reflection about the all-plus state with no flag wire, distinct from `GroverDiffusion`:

```wl
QuantumCircuitOperator["GroverPhaseDiffusion"[3]]["Diagram"]
```

---

The same phase reflection from zero-controls: an identical operator to `GroverPhaseDiffusion` in a shallower circuit:

```wl
QuantumCircuitOperator["GroverPhaseDiffusion0"[3]]["Diagram"]
```

---

The same block under its amplitude-named synonym, which builds only in call form (even `[]`), not as the bare string:

```wl
QuantumCircuitOperator["GroverAmplification"[3]]["Diagram"]
```

## The full search operator

The default search step over the built-in three-variable boolean function, with one oracle ancilla:

```wl
QuantumCircuitOperator["Grover"]["Diagram"]
```

---

The search step marking the states that satisfy a boolean formula:

```wl
QuantumCircuitOperator["Grover"[\[FormalX][1] && \[FormalX][2] && ! \[FormalX][3]]]["Diagram"]
```

---

A variable-to-wire list fixing which wire each variable occupies (both this and the bare two-variable form use three wires, with wire 3 the target):

```wl
QuantumCircuitOperator["Grover"[\[FormalX][1] && \[FormalX][2], {1, 2}]]["Diagram"]
```

---

A second argument that is an integer sets the oracle's target wire; choosing 1 collides with a variable, so the variables shift to 2,3 and the target lands on 1:

```wl
QuantumCircuitOperator["Grover"[\[FormalX][1] && \[FormalX][2], 1]]["Diagram"]
```

---

A `BooleanFunction` object accepted directly in place of a symbolic formula:

```wl
QuantumCircuitOperator["Grover"[BooleanFunction[7, 2]]]["Diagram"]
```

---

A prebuilt phase oracle wrapped with the matching phase diffusion, no formula compilation (here `CZ`, which sign-flips the marked $|11\rangle$; pairing a phase oracle with the plain `Grover` instead would mismatch the diffusion):

```wl
QuantumCircuitOperator["GroverPhase"[QuantumOperator["CZ"]]]["Diagram"]
```

---

The phase-oracle search step, which marks by sign flip on the data wires and so needs no ancilla:

```wl
QuantumCircuitOperator["GroverPhase"[\[FormalX][1] && \[FormalX][2]]]["Diagram"]
```

---

The search step built with the shallower zero-control diffusion in place of the X-sandwich form (the same reflection):

```wl
QuantumCircuitOperator["Grover0"[\[FormalX][1] && \[FormalX][2]]]["Diagram"]
```

## Does it search?

Two phase-Grover iterations on the uniform superposition concentrate the probability on the marked $|111\rangle$:

```wl
Nest[QuantumCircuitOperator["GroverPhase"[\[FormalX][1] && \[FormalX][2] && \[FormalX][3]]], QuantumState["+++"], 2]["Probabilities"]
```
