---
Template: Default
Title: Oracles
Context: Wolfram`QuantumFramework`
---

# Oracles

An oracle is the black box the search and query algorithms compile to: it marks
the basis states that satisfy a Boolean condition. The named oracles differ in
*how* they mark: flip a target qubit, flip a sign, rotate a target, or load a
whole function of the input. The first three take a Boolean formula and
variable-placement arguments; `GrayOracle` instead takes a function or angle
list, a precision, a target, and a rotation axis.

## Flag oracle: `BooleanOracle`

Flipping the target on the states that satisfy a formula, here a single multi-controlled gate for the one product term:

```wl
QuantumCircuitOperator[{"BooleanOracle"[\[FormalX][1] && \[FormalX][2] && \[FormalX][3]]}]["Diagram"]
```

---

A sum of exclusive-or terms compiles to one controlled gate per term (the ESOP form), so the gate count tracks the formula:

```wl
QuantumCircuitOperator[{"BooleanOracle"[Xor[\[FormalX][1], \[FormalX][2]]]}]["Diagram"]
```

---

A variable-to-wire list placing the variables on wires 1 and 3, leaving wire 2 idle and the target on wire 4:

```wl
QuantumCircuitOperator[{"BooleanOracle"[\[FormalX][1] || \[FormalX][2], {1, 3}]}]["Diagram"]
```

---

An integer in the second slot sets the target wire, here wire 5, leaving the variables on 1,2 and wires 3,4 idle:

```wl
QuantumCircuitOperator[{"BooleanOracle"[\[FormalX][1] && ! \[FormalX][2], 5]}]["Diagram"]
```

---

An operator in the second slot becomes the marking gate, turning the oracle into a controlled-*U* that applies an arbitrary unitary to the target on the satisfying states:

```wl
QuantumCircuitOperator[{"BooleanOracle"[! \[FormalX][1] && \[FormalX][2], QuantumOperator["RandomUnitary", "Label" -> "U"]]}]["Diagram"]
```

## Phase oracle: `PhaseOracle`

Flipping the sign of the satisfying states with no ancilla (width 2, not 3):

```wl
QuantumCircuitOperator[{"PhaseOracle"[Nand[\[FormalX][1], \[FormalX][2]]]}]["Diagram"]
```

---

On the uniform state $|{+}{+}\rangle$, the oracle for $x_1 \wedge x_2$ flips the sign of the marked $|11\rangle$ amplitude alone:

```wl
Normal[QuantumCircuitOperator[{"PhaseOracle"[\[FormalX][1] && \[FormalX][2]]}][QuantumState["++"]]["StateVector"]]
```

---

A variable-to-wire list fixing the wire of each variable:

```wl
QuantumCircuitOperator[{"PhaseOracle"[Implies[\[FormalX][1], \[FormalX][2]], {1, 2}]}]["Diagram"]
```

---

An integer in the second slot offsets every wire up by *m*, here placing the variables on wires 3,4 with 1,2 idle:

```wl
QuantumCircuitOperator[{"PhaseOracle"[\[FormalX][1] && \[FormalX][2], 2]}]["Diagram"]
```

## Rotation oracle: `BooleanOracleR`

Marking the states that satisfy a formula by a target rotation (`RZ` by π) rather than a flip:

```wl
QuantumCircuitOperator[{"BooleanOracleR"[\[FormalX][1] && \[FormalX][2]]}]["Diagram"]
```

---

A variable-to-wire list fixing each variable's wire:

```wl
QuantumCircuitOperator[{"BooleanOracleR"[Xor[\[FormalX][1], \[FormalX][2]], {1, 2}]}]["Diagram"]
```

---

A fourth argument sets the axis and angle of the marking rotation:

```wl
QuantumCircuitOperator[{"BooleanOracleR"[! \[FormalX][1] && ! \[FormalX][2], Automatic, Automatic, "RY"[Pi/2]]}]["Diagram"]
```

## Gray-code angle oracle: `GrayOracle`

A pure function loaded as target rotations $\theta_x = f(x)/2^{\text{prec}}$, sampled over the input values (default precision 4, so $R_Y(x^2/16)$ here):

```wl
QuantumCircuitOperator[{"GrayOracle"[(#^2 &)]}]["Diagram"]
```

---

An explicit list of angles in place of a function, padded with zeros to the $2^{\text{prec}}$ the default four input wires expect (16 angles):

```wl
QuantumCircuitOperator[{"GrayOracle"[{0, 1, 2, 3}]}]["Diagram"]
```

---

A precision argument sets the number of input wires, here 3 instead of 4, so the angles become $f(x)/8$ and the ladder halves:

```wl
QuantumCircuitOperator[{"GrayOracle"[(#^2 &), 3]}]["Diagram"]
```

---

A rotation-axis argument sets which rotation the angles drive, here `RZ`:

```wl
QuantumCircuitOperator[{"GrayOracle"[(#^2 &), 3, Automatic, "RZ"]}]["Diagram"]
```
