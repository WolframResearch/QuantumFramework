---
Template: Default
Title: Qudits and Angular Momentum
Context: Wolfram`QuantumFramework`
---

# Qudits and Angular Momentum

Beyond qubits, the library carries the clock-and-shift operators of a
$d$-level qudit, the qudit entangling gate, and angular-momentum operators whose
dimension is set by the spin. Each caption states the matrix or the action.

## Clock shift: `Shift`

The clock shift on a three-level qudit, cycling $|j\rangle \to |j+1 \bmod 3\rangle$:

```wl
QuantumCircuitOperator[{"Shift"[3]}]["Diagram"]
```

---

Acting on $|0\rangle$ it advances the qutrit to $|1\rangle$, amplitudes $(0, 1, 0)$:

```wl
Normal[QuantumCircuitOperator[{"Shift"[3]}][QuantumState[{1, 0, 0}, 3]]["StateVector"]]
```

## Clock phase: `ShiftPhase`

The companion clock phase, a diagonal of the cube roots of unity on a qutrit:

```wl
QuantumCircuitOperator[{"ShiftPhase"[3]}]["Diagram"]
```

## Qudit sum gate: `SUM`

The qudit analogue of CNOT, adding the control to the target modulo 3:

```wl
QuantumCircuitOperator[{"SUM"[3]}]["Diagram"]
```

## Sign flip: `FlipSign`

Flips the sign of the all-ones basis state on three qutrits, leaving every other state fixed:

```wl
QuantumCircuitOperator[{"FlipSign"[3]}]["Diagram"]
```

## Angular momentum: `JZ`

The spin-1 angular-momentum z-component; in its stored eigenbasis $\mathrm{diag}(-1, 0, 1)$ (the computational representation is $\mathrm{diag}(1, 0, -1)$):

```wl
QuantumCircuitOperator[{"JZ"[1]}]["Diagram"]
```

## Raising operator: `J+`

The spin-1 raising operator, $J_+|m\rangle = \sqrt{(j-m)(j+m+1)}\,|m+1\rangle$ with $J_+|j,j\rangle = 0$; for spin 1 the nonzero coefficients are $\sqrt2$:

```wl
QuantumCircuitOperator[{"J+"[1]}]["Diagram"]
```

## Spin rotation: `WignerD`

A spin rotation by Euler angles, here a $\pi/2$ z-rotation on a spin-1/2:

```wl
QuantumCircuitOperator[{"WignerD"[{Pi/2, 0, 0}]}]["Diagram"]
```
