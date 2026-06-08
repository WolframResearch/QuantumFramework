---
Template: Symbol
Name: QiskitTarget
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/QiskitTarget
Keywords: [qiskit, Target, transpiler, device, coupling map, basis gates, backend, connectivity]
SeeAlso: [QuantumQASM, QuantumCircuitOperator, QiskitCircuit]
RelatedGuides: [WolframQuantumComputationFramework]
---

## Usage

<code>[QiskitTarget]()[*spec*]</code> builds a qiskit transpiler target from the device specification *spec*, an Association of qiskit `Target.from_configuration` parameters.

<code>*tgt*["*prop*"]</code> gives the property *prop* of a `QiskitTarget` object *tgt*, one of `"NumQubits"`, `"OperationNames"`, `"CouplingMap"` or `"Spec"`.

## Details & Options

- A `QiskitTarget` describes a device for transpilation: its native gate set, qubit connectivity, and (optionally) per-instruction errors and durations. It is the device argument of [QuantumQASM]()'s `"Target"` option.
- *spec* is an Association whose keys are qiskit `Target.from_configuration` parameters: `"BasisGates"`, `"NumQubits"`, `"CouplingMap"`, `"InstructionDurations"`, `"ConcurrentMeasurements"`, `"Dt"`, `"TimingConstraints"`, `"CustomNameMapping"`. Keys are CamelCase and map to qiskit's snake_case names; the raw snake_case name is accepted too. The accepted key set is read from the installed qiskit, so it tracks the [qiskit Target.from_configuration documentation](https://docs.quantum.ibm.com/api/qiskit/qiskit.transpiler.Target#from_configuration) of whichever version is installed.
- Each value takes qiskit's own form: `"NumQubits"` an integer; `"CouplingMap"` a list of directed integer pairs `{{c, t}, …}`; `"Dt"` a number of seconds; `"InstructionDurations"`, `"ConcurrentMeasurements"`, `"TimingConstraints"` and `"CustomNameMapping"` the structures qiskit's `Target.from_configuration` expects.
- `"BasisGates"` is a list of qiskit *standard gate names* — the names are defined by the installed qiskit, not by `QiskitTarget`, and vary by version. In the bundled qiskit 2.4.1 there are 54, including the single-qubit `"x"`, `"y"`, `"z"`, `"h"`, `"s"`, `"sx"`, `"t"`, `"rx"`, `"ry"`, `"rz"`, `"p"`, `"u"`, `"u3"`, the two-qubit `"cx"`, `"cz"`, `"cy"`, `"ch"`, `"swap"`, `"iswap"`, `"ecr"`, `"rzz"`, `"rxx"`, and the three-qubit `"ccx"`, `"ccz"`, `"cswap"`. The names listed here are those of qiskit 2.4.1; a different qiskit version may add, rename, or remove gate names, so the list could differ from what is shown. The complete, version-correct list is qiskit's standard gate-name mapping, obtainable with `qiskit.circuit.library.get_standard_gate_name_mapping` or from qiskit's [circuit-library documentation](https://docs.quantum.ibm.com/api/qiskit/circuit_library). A name outside that set must be supplied through `"CustomNameMapping"`.
- The keys are validated against the installed qiskit at call time, so an unknown key gives a [Failure]() (see Possible Issues) and any current or future qiskit `Target` parameter is accepted without a change to `QiskitTarget`.
- `QiskitTarget` stores the validated specification rather than a serialized qiskit `Target`; the live target is rebuilt from the spec when it is used (inside [QuantumQASM]()), so the same object works across qiskit sessions and versions.
- `measure` and `reset` are added on every qubit automatically — they are universal but `Target.from_configuration` omits them for a bare unitary basis — so a target built from `{"x", "sx", "rz", "cz"}` still transpiles circuits that contain measurements.
- A `QiskitTarget` displays as a summary panel showing its qubit count, native gate names, and number of coupling edges.

| Property | Result |
|---|---|
| `"NumQubits"` | the number of qubits |
| `"OperationNames"` | the sorted native operation names, including `measure` and `reset` |
| `"CouplingMap"` | the directed coupling pairs |
| `"Spec"` | the normalized specification Association |

## Basic Examples

A device target from a native gate set, a qubit count, and a coupling map:

```wl
QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}|>]
```

<!-- => QiskitTarget summary panel — Qubits: 4, Gates: {cz, measure, reset, rz, sx, x}, Couplings: 4 -->

---

A larger, linearly connected device:

```wl
QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 5, "CouplingMap" -> {{0, 1}, {1, 2}, {2, 3}, {3, 4}}|>]
```

<!-- => QiskitTarget summary panel — Qubits: 5, Gates: {cz, measure, reset, rz, sx, x}, Couplings: 4 -->

---

Once built, query a property — here the qubit count:

```wl
QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}|>]["NumQubits"]
```

<!-- => 4 -->

## Scope

The native operation names, including the automatically added `measure` and `reset`:

```wl
QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}|>]["OperationNames"]
```

<!-- => {cz, measure, reset, rz, sx, x} -->

---

The directed coupling pairs:

```wl
QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}|>]["CouplingMap"]
```

<!-- => {{1, 0}, {1, 3}, {2, 0}, {2, 3}} -->

---

The normalized device specification, with keys in qiskit's snake_case form:

```wl
QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}|>]["Spec"]
```

<!-- => <|basis_gates -> {x, sx, rz, cz}, num_qubits -> 4, coupling_map -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}|> -->

### Case-insensitive specification keys

CamelCase and qiskit's snake_case spellings of the spec keys are interchangeable:

```wl
With[{cm = {{1, 0}, {1, 3}, {2, 0}, {2, 3}}, b = {"x", "sx", "rz", "cz"}},
  QiskitTarget[<|"BasisGates" -> b, "NumQubits" -> 4, "CouplingMap" -> cm|>]["Spec"] ===
    QiskitTarget[<|"basis_gates" -> b, "num_qubits" -> 4, "coupling_map" -> cm|>]["Spec"]
]
```

<!-- => True -->

### Use with QuantumQASM

A `QiskitTarget` is the device argument of [QuantumQASM](); the circuit is transpiled to the target's basis and connectivity:

```wl
QuantumQASM[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], "Target" -> QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}|>]]
```

<!-- => OPENQASM 3.0; ... native program on the target ... -->

## Applications

Build a device target once and reuse it to transpile several circuits to native OpenQASM:

```wl
With[{tgt = QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}|>]},
  StringStartsQ[QuantumQASM[#, "Target" -> tgt], "OPENQASM 3.0"] & /@ {QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], QuantumCircuitOperator["GHZ"[2]]}
]
```

<!-- => {True, True} -->

## Properties and Relations

Passing a `QiskitTarget` is equivalent to passing its specification Association inline to [QuantumQASM]()'s `"Target"` option — the object simply validates and carries the spec:

```wl
With[{circ = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], spec = <|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}|>},
  QuantumQASM[circ, "Target" -> QiskitTarget[spec]] === QuantumQASM[circ, "Target" -> spec]
]
```

<!-- => True -->

## Possible Issues

A specification key that qiskit's `Target` does not recognize gives a [Failure]() naming the unknown key and the accepted set:

```wl
QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "Bogus" -> 1|>]
```

<!-- => Failure["QiskitTarget", <|... "Unknown" -> {"bogus"}, ...|>] : Unknown Target option(s): {bogus}. -->

---

`QiskitTarget` also returns a `Failure` if the qiskit `Target` API is unavailable in the framework's Python session (for example a qiskit too old to provide `Target.from_configuration`); in a working session, where qiskit is installed automatically, this does not arise.
