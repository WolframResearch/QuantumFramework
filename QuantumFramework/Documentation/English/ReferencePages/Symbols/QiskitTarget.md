---
Template: Symbol
Name: QiskitTarget
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/QiskitTarget
Keywords: [qiskit, Target, transpiler, device, coupling map, basis gates, backend, connectivity]
SeeAlso: [QuantumQASM, QuantumCircuitOperator, QiskitCircuit]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: [IBMQuantumErrorMap]
---

## Usage

<code>[QiskitTarget]()[*spec*]</code> builds a qiskit transpiler target from the device specification *spec*, an Association of qiskit `Target.from_configuration` parameters.

<code>[QiskitTarget]()["FromBackend" -> *name*, "Provider" -> *p*]</code> builds an error-aware target by reading the live backend *name*'s own qiskit `Target`, including its per-instruction error and duration; *p* selects the provider (`"IBMProvider"`, `"AWSProvider"`, or `"GenericV2"` for a local fake device).

## Details & Options

- A `QiskitTarget` describes a device for transpilation: its native gate set, qubit connectivity, and (optionally) per-instruction errors and durations. It is the device argument of [QuantumQASM]()'s `"Target"` option.
- *spec* is an Association whose keys are qiskit `Target.from_configuration` parameters: `"BasisGates"`, `"NumQubits"`, `"CouplingMap"`, `"InstructionDurations"`, `"ConcurrentMeasurements"`, `"Dt"`, `"TimingConstraints"`, `"CustomNameMapping"`. Keys are CamelCase and map to qiskit's snake_case names; the raw snake_case name is accepted too. The accepted key set is read from the installed qiskit, so it tracks the [qiskit Target.from_configuration documentation](https://docs.quantum.ibm.com/api/qiskit/qiskit.transpiler.Target#from_configuration) of whichever version is installed.
- Each value takes qiskit's own form: `"NumQubits"` an integer; `"CouplingMap"` a list of directed integer pairs `{{c, t}, …}`; `"Dt"` a number of seconds; `"InstructionDurations"`, `"ConcurrentMeasurements"`, `"TimingConstraints"` and `"CustomNameMapping"` the structures qiskit's `Target.from_configuration` expects.
- `"BasisGates"` is a list of qiskit *standard gate names*. The names are defined by the installed qiskit, not by `QiskitTarget`, and vary by version. In the bundled qiskit 2.4.1 there are 54, including the single-qubit `"x"`, `"y"`, `"z"`, `"h"`, `"s"`, `"sx"`, `"t"`, `"rx"`, `"ry"`, `"rz"`, `"p"`, `"u"`, `"u3"`, the two-qubit `"cx"`, `"cz"`, `"cy"`, `"ch"`, `"swap"`, `"iswap"`, `"ecr"`, `"rzz"`, `"rxx"`, and the three-qubit `"ccx"`, `"ccz"`, `"cswap"`. The names listed here are those of qiskit 2.4.1; a different qiskit version may add, rename, or remove gate names, so the list could differ from what is shown. The complete, version-correct list is qiskit's standard gate-name mapping, obtainable with `qiskit.circuit.library.get_standard_gate_name_mapping` or from qiskit's [circuit-library documentation](https://docs.quantum.ibm.com/api/qiskit/circuit_library). A name outside that set must be supplied through `"CustomNameMapping"`.
- The keys are validated against the installed qiskit at call time, so an unknown key gives a [Failure]() (see Possible Issues) and any current or future qiskit `Target` parameter is accepted without a change to `QiskitTarget`.
- `QiskitTarget` stores the validated specification rather than a serialized qiskit `Target`; the live target is rebuilt from the spec when it is used (inside [QuantumQASM]()), so the same object works across qiskit sessions and versions.
- `measure` and `reset` are added on every qubit automatically (they are universal but `Target.from_configuration` omits them for a bare unitary basis), so a target built from `{"x", "sx", "rz", "cz"}` still transpiles circuits that contain measurements.
- `"InstructionProperties"` is a reserved spec key that makes transpilation **error-aware**: a list of rows `{*gate*, {*qubits*…}, <|"Error" -> *e*, "Duration" -> *d*|>}` (either field may be omitted). After the target is built, each row is applied with qiskit's `Target.update_instruction_properties`, so `optimization_level` 1 and above bias layout and routing toward low-error qubits and pairs. A row whose gate or qubits are not in the target is skipped (recorded on the transpiled circuit's metadata, surfaced by `QuantumQASM::skippedProps`), never applied blindly; an `error` of `measure` rows also steers the layout away from high-readout-error qubits.
- `[QiskitTarget]()["FromBackend" -> *name*, …]` reads `*name*`'s populated `Target` and emits exactly this `"InstructionProperties"` spec, so a circuit transpiled against the result lands on the device's lowest-error qubits. The same spec can be produced for non-qiskit devices (for example by `QUALink`'s `QUAConnect[…]["TargetSpec"]` for the IQCC `arbel` QPU), so one error-aware-target mechanism serves every backend.
- A `QiskitTarget` displays as a summary panel showing its qubit count, native gate names, number of coupling edges, and whether it carries per-instruction error data.
- A property of a `QiskitTarget` object *tgt* is read with <code>*tgt*["*prop*"]</code>, where *prop* is one of:

| Property | Result |
|---|---|
| `"NumQubits"` | the number of qubits |
| `"OperationNames"` | the sorted native operation names, including `measure` and `reset` |
| `"CouplingMap"` | the directed coupling pairs |
| `"InstructionProperties"` | the per-instruction error/duration rows, or `Missing` if none |
| `"Spec"` | the normalized specification Association |

## Basic Examples

A device target from a native gate set, a qubit count, and a coupling map:

```wl
QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}|>]
```

<!-- => QiskitTarget summary panel, Qubits: 4, Gates: {cz, measure, reset, rz, sx, x}, Couplings: 4 -->

---

A larger, linearly connected device:

```wl
QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 5, "CouplingMap" -> {{0, 1}, {1, 2}, {2, 3}, {3, 4}}|>]
```

<!-- => QiskitTarget summary panel, Qubits: 5, Gates: {cz, measure, reset, rz, sx, x}, Couplings: 4 -->

---

Once built, query a property, here the qubit count:

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

### Error-aware targets

The reserved `"InstructionProperties"` key attaches a per-instruction error to the target, which makes the layout and routing passes prefer low-error qubits and pairs. Build a 4-qubit ring in which one `cz` edge, `{1, 2}`, is far cleaner than the rest:

```wl
ring  = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3}, {3, 2}, {3, 0}, {0, 3}};
props = Join[
   ({"cz", #, <|"Error" -> 0.2|>} &) /@ {{0, 1}, {1, 0}, {2, 3}, {3, 2}, {3, 0}, {0, 3}},
   ({"cz", #, <|"Error" -> 0.001|>} &) /@ {{1, 2}, {2, 1}}];
tgt = QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> ring, "InstructionProperties" -> props|>];
Length @ tgt["InstructionProperties"]
```

<!-- => 8 -->

At `optimization_level` 3 the Bell's `cz` is routed onto that low-error edge:

```wl
qc = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}];
First @ StringCases[QuantumQASM[qc, "Target" -> tgt, "OptimizationLevel" -> 3, "SeedTranspiler" -> 0],
   "cz $" ~~ a : DigitCharacter .. ~~ ", $" ~~ b : DigitCharacter .. :> {FromDigits @ a, FromDigits @ b}]
```

<!-- => {1, 2} -->

Drop the error rows and the same circuit no longer prefers that edge, it places by connectivity alone:

```wl
First @ StringCases[QuantumQASM[qc, "Target" -> QiskitTarget[KeyDrop[tgt["Spec"], "instruction_properties"]], "OptimizationLevel" -> 3, "SeedTranspiler" -> 0],
   "cz $" ~~ a : DigitCharacter .. ~~ ", $" ~~ b : DigitCharacter .. :> {FromDigits @ a, FromDigits @ b}]
```

<!-- => {0, 1} -->

### Reading a backend's error model

`["FromBackend" -> *name*, "Provider" -> *p*]` reads a live backend's own populated `Target` instead of building one by hand. `"GenericV2"` is a local fake device (no credentials, no network) useful for testing the path; `"IBMProvider"` and `"AWSProvider"` reach real hardware (see Applications). Here the fake device's target arrives already error-populated:

```wl
gen = QiskitTarget["FromBackend" -> 5, "Provider" -> "GenericV2"];
{gen["NumQubits"], Length @ gen["InstructionProperties"]}
```

<!-- => {5, 45} -->

Each row carries the backend's measured error and duration:

```wl
SelectFirst[gen["InstructionProperties"], MatchQ[#, {"cx", _, _}] &]
```

<!-- => {"cx", {0, 1}, <|"Error" -> 0.00278, "Duration" -> 2.662*^-7|>} -->

## Applications

Build a device target once and reuse it to transpile several circuits to native OpenQASM:

```wl
With[{tgt = QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> {{1, 0}, {1, 3}, {2, 0}, {2, 3}}|>]},
  StringStartsQ[QuantumQASM[#, "Target" -> tgt], "OPENQASM 3.0"] & /@ {QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], QuantumCircuitOperator["GHZ"[2]]}
]
```

<!-- => {True, True} -->

### Error-aware routing on a real IBM device

On real hardware every qubit and every two-qubit pair has a different error rate, so *where* a circuit is placed sets its fidelity. `["FromBackend"]` reads the device's own calibrated error model and lets [QuantumQASM]() choose that placement automatically. (Outputs below are recorded from a live IBM backend; `"FromBackend"` authenticates through the same stored IBM credentials as [IBMJobSubmit]().)

Read the live backend's target, already populated, no assembly needed:

```wl
#| eval: false
tgt = QiskitTarget["FromBackend" -> Automatic, "Provider" -> "IBMProvider"];
{tgt["NumQubits"], Length @ tgt["InstructionProperties"]}
```

<!-- => {156, 1288}  (a 156-qubit device, with 1288 per-instruction error rows) -->

The errors span a wide range across the chip, so qubit choice matters: `cz` error runs from 0.12% to a stuck 100%, and readout from 0.2% to 32%:

```wl
#| eval: false
{MinMax @ Cases[tgt["InstructionProperties"], {"cz", {_, _}, p_} :> p["Error"]],
 MinMax @ Cases[tgt["InstructionProperties"], {"measure", {_}, p_} :> p["Error"]]}
```

<!-- => {{0.0012, 1.}, {0.0023, 0.3221}} -->

Transpile a Bell against the error-aware target and against a connectivity-only one (the same coupling map with the error rows dropped), and read off the physical pair each chose:

```wl
#| eval: false
qc      = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}];
naive   = QiskitTarget[KeyDrop[tgt["Spec"], "instruction_properties"]];
pairOf[t_] := First @ StringCases[QuantumQASM[qc, "Target" -> t, "OptimizationLevel" -> 3],
   "cz $" ~~ a : DigitCharacter .. ~~ ", $" ~~ b : DigitCharacter .. :> {FromDigits @ a, FromDigits @ b}];
{pairOf[tgt], pairOf[naive]}
```

<!-- => {{25, 37}, {0, 1}}  (error-aware picks {25, 37}; connectivity-only takes the {0, 1} corner) -->

Scoring each chosen pair by its total expected error, that is its `cz` error plus the two readout errors, the error-aware placement is about four times cleaner, **1.11%** versus **4.88%**, decided automatically from the device's own published numbers with no hand-screening of qubits:

```wl
#| eval: false
budget[t_, {a_, b_}] := With[{ip = t["InstructionProperties"]},
   SelectFirst[ip, MatchQ[#, {"cz", {a, b} | {b, a}, _}] &][[3]]["Error"] +
    SelectFirst[ip, MatchQ[#, {"measure", {a}, _}] &][[3]]["Error"] +
    SelectFirst[ip, MatchQ[#, {"measure", {b}, _}] &][[3]]["Error"]];
{budget[tgt, {25, 37}], budget[tgt, {0, 1}]}
```

<!-- => {0.0111, 0.0488} -->

The emitted QASM is ordinary native, physical-qubit OpenQASM 3, so the rest of the workflow is unchanged: the same error-aware circuit can be submitted through [IBMJobSubmit]() (or, for the IQCC `arbel` QPU, through QUALink's gate-level pipeline). The one-line change, `"FromBackend"` instead of a hand-built target, places the circuit on the device's best qubits for free.

## Properties and Relations

Passing a `QiskitTarget` is equivalent to passing its specification Association inline to [QuantumQASM]()'s `"Target"` option. The object simply validates and carries the spec:

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
