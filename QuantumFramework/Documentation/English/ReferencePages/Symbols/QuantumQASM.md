---
Template: Symbol
Name: QuantumQASM
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/QuantumQASM
Keywords: [OpenQASM, QASM, qiskit, transpile, transpiler, native gates, coupling map, hardware, quantum circuit]
SeeAlso: [QiskitTarget, QuantumCircuitOperator, QiskitCircuit, ImportQASMCircuit]
RelatedGuides: [WolframQuantumComputationFramework]
---

## Usage

<code>[QuantumQASM]()[*circ*]</code> gives the OpenQASM 3 string for the quantum circuit *circ*, serialized exactly as built.

<code>[QuantumQASM]()[*circ*, {$g_1$, $g_2$, …}]</code> transpiles *circ* to the native gate set $\{g_1, g_2, \dots\}$ and gives the resulting native OpenQASM.

<code>[QuantumQASM]()[*circ*, *basis*, *opts*]</code> applies qiskit transpiler options *opts* such as `"CouplingMap"`, `"OptimizationLevel"` and `"Target"`.

<code>[QuantumQASM]()[*circ*, "Version" -> 2]</code> gives OpenQASM 2 instead of the default OpenQASM 3.

## Details & Options

- *circ* can be a [QuantumCircuitOperator]() or a [QiskitCircuit]() (the value of `qco["Qiskit"]`).
- `QuantumQASM` is **faithful by default**: with no native gate set and no `"Target"`, `"CouplingMap"` or other transpiler option, *circ* is serialized exactly as built, with no qubit routing or optimization. Supplying any native gate set, target, coupling map or transpiler option makes `QuantumQASM` transpile *circ* first, then serialize the result.
- `QuantumQASM` returns an OpenQASM `String`. On any error it returns a [Failure]() instead of a partial program (see Possible Issues). It has no further property or accessor calls; the device model it accepts through `"Target"`, a [QiskitTarget]() object, is what carries accessors such as `["NumQubits"]` and `["OperationNames"]` (shown under Options).
- The native gate set is a list of qiskit *standard gate names* such as `{"x", "sx", "rz", "cz"}`; it can equivalently be given through the `"BasisGates"` option. The names are defined by the installed qiskit (54 in qiskit 2.4.1, e.g. `"x"`, `"sx"`, `"rz"`, `"cz"`, `"h"`, `"cx"`, `"ecr"`, `"swap"`, `"ccx"`); this set can differ between qiskit versions, so the complete, version-correct list is qiskit's standard gate-name mapping, obtainable with `qiskit.circuit.library.get_standard_gate_name_mapping` or from qiskit's [circuit-library documentation](https://docs.quantum.ibm.com/api/qiskit/circuit_library).
- Every option other than `"Version"` is forwarded to qiskit's transpiler. Option names are CamelCase and map to qiskit's snake_case parameter names (`"CouplingMap"` to `coupling_map`); the raw snake_case name is accepted too. Any current or future qiskit transpile parameter therefore works without a change to `QuantumQASM`.
- A `"CouplingMap"` is a list of directed pairs `{{c, t}, …}` of physical qubit indices that may interact; two-qubit gates are emitted in that calibrated direction.
- A qiskit-free alternative is the pure-Wolfram emitter [QuantumCircuitOperator]() property `"QASM"`, which emits a generic gate set; `QuantumQASM` instead routes through qiskit and can decompose to a hardware-native basis.
- `QuantumQASM` builds the qiskit circuit through `qco["Qiskit"]` and runs the qiskit transpiler in the framework's Python session; the qiskit package is installed automatically on first use.

| Option | Default value | Description |
|---|---|---|
| `"Version"` | 3 | the OpenQASM version to emit, 2 or 3 |
| `"BasisGates"` | Automatic | the native gate set (the same as the second argument) |
| `"CouplingMap"` | Automatic | directed pairs of physical qubits that may interact |
| `"OptimizationLevel"` | Automatic | the qiskit optimization level, 0 through 3 |
| `"Target"` | None | a [QiskitTarget]() object or an Association device spec |

Any other qiskit `transpile` parameter (for example `"InitialLayout"`, `"RoutingMethod"`, `"SeedTranspiler"`) is also accepted.

The string value of a *method* option — `"LayoutMethod"`, `"RoutingMethod"`, `"TranslationMethod"`, `"SchedulingMethod"`, `"InitMethod"`, `"OptimizationMethod"`, `"UnitarySynthesisMethod"` — is the name of a qiskit transpiler *stage plugin*. These names are defined by the installed qiskit, not by `QuantumQASM`, and vary by version; the authoritative list is qiskit's [transpilation defaults and configuration options](https://quantum.cloud.ibm.com/docs/en/guides/defaults-and-configuration-options) guide, and they can be enumerated from the installed qiskit through `qiskit.transpiler.preset_passmanagers.plugin.list_stage_plugins`. For the qiskit 2.4.1 bundled with the framework the values are:

| Method option | qiskit 2.4.1 plugin names |
|---|---|
| `"LayoutMethod"` | `"trivial"`, `"dense"`, `"sabre"`, `"default"` |
| `"RoutingMethod"` | `"basic"`, `"lookahead"`, `"sabre"`, `"none"`, `"default"` |
| `"TranslationMethod"` | `"translator"`, `"synthesis"`, `"default"` (plus backend-specific `"ibm_…"` variants) |
| `"SchedulingMethod"` | `"asap"`, `"alap"`, `"default"` |
| `"UnitarySynthesisMethod"` | `"default"`, `"aqc"`, `"clifford"`, `"gridsynth"`, `"sk"` |
| `"InitMethod"`, `"OptimizationMethod"` | `"default"` |

The table above is a snapshot of qiskit 2.4.1. A different qiskit version may add, rename, or remove these plugin names, so the list could differ from what is shown here — the authoritative set is always whatever the installed qiskit reports (queryable through `list_stage_plugins`), not this page. `QuantumQASM` itself does not depend on the table: it forwards whatever value you give straight to qiskit, so a value the installed qiskit does not recognize comes back as a `Failure` carrying qiskit's own message.

`"OptimizationLevel"` takes an integer 0 through 3; `"SeedTranspiler"` takes any integer.

## Basic Examples

The OpenQASM 3 of a circuit, serialized exactly as built:

```wl
QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}]]
```

<!-- => OPENQASM 3.0; include "stdgates.inc"; bit[2] c; qubit[2] q; h q[0]; cx q[0], q[1]; c[0] = measure q[0]; c[1] = measure q[1]; -->

---

Transpile the same circuit to a native gate set:

```wl
QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], {"x", "sx", "rz", "cz"}]
```

<!-- => OPENQASM 3.0; ... rz(pi/2) q[0]; sx q[0]; rz(pi/2) q[0]; ... (native sx/rz/cz, no h/cx) -->

---

Route onto a device coupling map, emitting two-qubit gates in the calibrated direction:

```wl
QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], {"x", "sx", "rz", "cz"}, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}, "OptimizationLevel" -> 3]
```

<!-- => OPENQASM 3.0; ... physical qubits $0/$1; cz $1, $0; (directed) -->

## Scope

A [QiskitCircuit]() (the value of `qco["Qiskit"]`) is accepted directly:

```wl
QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}]["Qiskit"], {"x", "sx", "rz", "cz"}]
```

<!-- => OPENQASM 3.0; ... native sx/rz/cz ... -->

---

Option names are case-insensitive: CamelCase and qiskit's snake_case spellings are interchangeable:

```wl
With[{qc = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], cm = {{1, 0}, {1, 3}, {2, 0}, {2, 3}}, b = {"x", "sx", "rz", "cz"}},
  QuantumQASM[qc, b, "CouplingMap" -> cm, "OptimizationLevel" -> 1] === QuantumQASM[qc, b, "couplingMap" -> cm, "optimization_level" -> 1]
]
```

<!-- => True -->

---

Any qiskit transpiler parameter is forwarded, here the routing method and seed:

```wl
QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], {"x", "sx", "rz", "cz"}, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}, "RoutingMethod" -> "sabre", "SeedTranspiler" -> 42]
```

<!-- => OPENQASM 3.0; ... native routed program ... -->

## Options

### "Version"

Emit OpenQASM 2 instead of the default OpenQASM 3:

```wl
QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], "Version" -> 2]
```

<!-- => OPENQASM 2.0; include "qelib1.inc"; qreg q[2]; creg c[2]; h q[0]; cx q[0],q[1]; measure q[0] -> c[0]; measure q[1] -> c[1]; -->

### "CouplingMap"

Restrict two-qubit gates to a device's connected qubit pairs, given as directed pairs:

```wl
QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], {"x", "sx", "rz", "cz"}, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}]
```

<!-- => OPENQASM 3.0; ... cz $1, $0; (directed onto the coupling map) -->

### "OptimizationLevel"

Higher optimization levels can shrink the program; compare the OpenQASM length at levels 0 and 3:

```wl
StringLength @ QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], {"x", "sx", "rz", "cz"}, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}, "OptimizationLevel" -> #] & /@ {0, 3}
```

<!-- => {196, 181} -->

### "Target"

A [QiskitTarget]() bundles a native basis with directed connectivity; transpile against it:

```wl
QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], "Target" -> QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}|>]]
```

<!-- => OPENQASM 3.0; ... native program on the target ... -->

---

A [QiskitTarget]() carries device accessors, including the auto-added `measure` and `reset`:

```wl
With[{tgt = QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}|>]},
  {tgt["NumQubits"], tgt["OperationNames"]}
]
```

<!-- => {4, {cz, measure, reset, rz, sx, x}} -->

---

The target can also be given inline as an Association:

```wl
QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], "Target" -> <|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}|>]
```

<!-- => OPENQASM 3.0; ... native program ... -->

## Applications

Take an algorithm built in the framework and emit native, routed OpenQASM 3 ready for a hardware compiler — decomposed to the device gate set and mapped onto its connectivity:

```wl
QuantumQASM[QuantumCircuitOperator["GHZ"[3]], {"x", "sx", "rz", "cz"}, "CouplingMap" -> {{0, 1}, {1, 2}}, "OptimizationLevel" -> 3]
```

<!-- => OPENQASM 3.0; ... three-qubit native program routed on the line {0,1,2} ... -->

## Possible Issues

An option qiskit does not recognize gives a [Failure]() that names the unknown key and the accepted set:

```wl
QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], {"x", "sx", "rz", "cz"}, "NoSuchParam" -> 5]
```

<!-- => Failure["QuantumQASM", <|... "Unknown" -> {"no_such_param"}, "Accepted" -> {...}, "QiskitVersion" -> "2.4.1"|>] : Unknown transpile option(s): {no_such_param}. -->

---

A non-universal native gate set that qiskit cannot compile to gives a transpile Failure rather than a partial program:

```wl
QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], {"x", "cz"}]
```

<!-- => Failure["QuantumQASM", <|... "PythonError" -> "TranspilerError: Unable to translate ..."|>] : qiskit transpile failed. -->

---

Passing the same option under two spellings is rejected as ambiguous:

```wl
QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], {"x", "sx", "rz", "cz"}, "CouplingMap" -> {{1, 0}}, "coupling_map" -> {{1, 0}}]
```

<!-- => Failure["QuantumQASM", ...] : Ambiguous duplicate option after key normalization (for example CouplingMap and coupling_map together). -->

---

Giving the native gates both positionally and through the `"BasisGates"` option is an error:

```wl
QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], {"x", "sx", "rz", "cz"}, "BasisGates" -> {"x"}]
```

<!-- => Failure["QuantumQASM", ...] : Native gates given both positionally and as the BasisGates option. -->

---

A gate the chosen OpenQASM version cannot express makes the export fail; here a global phase in OpenQASM 2:

```wl
QuantumQASM[QuantumCircuitOperator[{"GlobalPhase"[0.3] -> 1}], "Version" -> 2]
```

<!-- => Failure["QuantumQASM", <|... "PythonError" -> ...|>] : OpenQASM2 export failed. -->
