---
Template: Symbol
Name: QiskitCircuit
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/QiskitCircuit
Keywords: [qiskit, circuit, OpenQASM, QASM, transpile, QPY, IBM, interoperability, pickled circuit]
SeeAlso: [QuantumCircuitOperator, QuantumQASM, QiskitTarget, IBMJobSubmit]
RelatedGuides: [WolframQuantumComputationFramework]
---

## Usage

<code>[QiskitCircuit]()[*qco*]</code> builds a qiskit circuit from the [QuantumCircuitOperator]() *qco*.

<code>[QiskitCircuit]()["*qasm*"]</code> builds a qiskit circuit from OpenQASM source *qasm*.

<code>*qkc*["*prop*"]</code> gives the property *prop* of a `QiskitCircuit` object *qkc*.

<code>*qkc*[]</code> runs *qkc* on a simulator or backend, and <code>*qkc*[*state*]</code> runs it on the initial [QuantumState]() *state*.

## Details & Options

- A `QiskitCircuit` is the bridge between [QuantumFramework]() and qiskit: it wraps a pickled qiskit `QuantumCircuit` as a [ByteArray]() and exposes it through a property interface.
- A `QiskitCircuit` is produced by [QuantumCircuitOperator]()'s `"Qiskit"` property, by [QuantumQASM]() import, and by the transform properties below; <code>[QuantumCircuitOperator]()[*qkc*]</code> converts one back to a [QuantumCircuitOperator]().
- <code>[QiskitCircuit]()["*qasm*"]</code> accepts OpenQASM 2.0 or 3.0 source: the source is read with [QuantumCircuitOperator]() and re-emitted to qiskit, so it is equivalent to <code>[QuantumCircuitOperator]()["*qasm*"]["Qiskit"]</code>. A string that is not OpenQASM gives a [Failure]() (see Possible Issues).
- Property evaluation runs in the framework's qiskit Python session, which is set up automatically. The hardware-facing properties `"Transpile"`, `"Layout"`, `"Validate"`, `"QPY"`, and running *qkc* on a provider additionally need that provider configured (an IBM or AWS account); the default for `qkc[]` is a local Aer simulator.
- A `QiskitCircuit` displays as a summary panel showing its qubit count, depth, and gate counts.

| Property | Result |
|---|---|
| `"Qubits"` | the number of qubits |
| `"Depth"` | the circuit depth |
| `"Ops"` | the gate counts, an Association of gate name to count |
| `"Clbits"` | the classical-bit indices |
| `"QuantumCircuitOperator"` | convert to a [QuantumCircuitOperator]() (also `"QuantumCircuit"`) |
| `"QuantumOperator"` | the circuit unitary as a [QuantumOperator]() |
| `"Matrix"` | the unitary matrix |
| `"Diagram"` | qiskit's Matplotlib drawing of the circuit |
| `"Graphics"` | the circuit drawn through LaTeX |
| `"QASM"`, `"QASM2"` | OpenQASM 2.0 source |
| `"QASM3"` | OpenQASM 3.0 source |
| `"QPY"` | compressed qiskit QPY bytes |
| `"Decompose"`, *n* | decompose each gate *n* times, giving a `QiskitCircuit` |
| `"Transpile"`, *basis* | transpile to a native *basis* or backend, giving a `QiskitCircuit` |
| `"Layout"` | a plot of the transpiled qubit layout on a backend |
| `"Validate"` | validate the circuit on an IBM backend |
| `"Bytes"` | the underlying pickled qiskit circuit, a [ByteArray]() |

## Basic Examples

A circuit in [QuantumFramework](), drawn with its [QuantumCircuitOperator]() `"Diagram"`:

```wl
QuantumCircuitOperator["Toffoli"]["Diagram"]
```

<!-- => QuantumFramework circuit diagram of the Toffoli gate -->

---

[QiskitCircuit]() builds the same circuit as a qiskit object, shown as a summary panel of its qubit count, depth, and gate counts:

```wl
QiskitCircuit[QuantumCircuitOperator["Toffoli"]]
```

<!-- => QiskitCircuit summary panel — Qubits: 3, Depth: 11, Ops: <|cx -> 6, t -> 4, tdg -> 3, h -> 2|> -->

---

The `"Diagram"` property draws it the way qiskit does, so the two diagrams above and below are the same circuit in each framework's notation:

```wl
QiskitCircuit[QuantumCircuitOperator["Toffoli"]]["Diagram"]
```

<!-- => qiskit's drawing of the same circuit -->

## Scope

### Convert between QuantumCircuitOperator and QiskitCircuit

A `QiskitCircuit` converts back to a [QuantumCircuitOperator]() with [QuantumCircuitOperator]()[*qkc*]:

```wl
QuantumCircuitOperator[QiskitCircuit[QuantumCircuitOperator["Toffoli"]]]
```

<!-- => QuantumCircuitOperator summary box -->

---

Drawn again as a [QuantumFramework]() circuit, it matches the original, completing the round trip:

```wl
With[{qkc = QiskitCircuit[QuantumCircuitOperator["Toffoli"]]}, QuantumCircuitOperator[qkc]["Diagram"]]
```

<!-- => QuantumFramework diagram of the round-tripped circuit -->

### Import OpenQASM source

A string of OpenQASM source builds a qiskit circuit directly:

```wl
QiskitCircuit["OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\nh q[0];\ncx q[0], q[1];\n"]
```

<!-- => QiskitCircuit summary panel — Qubits: 2 -->

### Circuit structure

The circuit depth:

```wl
QiskitCircuit[QuantumCircuitOperator["Toffoli"]]["Depth"]
```

<!-- => 11 -->

---

The gate counts, an Association of qiskit gate name to count:

```wl
QiskitCircuit[QuantumCircuitOperator["Toffoli"]]["Ops"]
```

<!-- => <|cx -> 6, t -> 4, tdg -> 3, h -> 2|> -->

### Decompose

`"Decompose"` expands each gate one level, giving another `QiskitCircuit`:

```wl
QiskitCircuit[QuantumCircuitOperator["Toffoli"]]["Decompose"]
```

<!-- => QiskitCircuit summary panel (decomposed) -->

---

Its diagram shows the lower-level gates the decomposition introduces:

```wl
QiskitCircuit[QuantumCircuitOperator["Toffoli"]]["Decompose"]["Diagram"]
```

<!-- => qiskit diagram of the decomposed circuit -->

### Unitary

The circuit unitary as a [QuantumOperator]():

```wl
QiskitCircuit[QuantumCircuitOperator["Toffoli"]]["QuantumOperator"]
```

<!-- => QuantumOperator summary box (3 qubits) -->

### Export OpenQASM

Export the circuit as OpenQASM 3.0 source:

```wl
QiskitCircuit[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["QASM3"]
```

<!-- => "OPENQASM 3.0; ..." -->

### Run on a simulator

Running a `QiskitCircuit` on the default local simulator gives a [QuantumState]():

```wl
QiskitCircuit[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]][]
```

<!-- => QuantumState summary box (Bell state) -->

## Possible Issues

A string that is not OpenQASM source cannot be read as a circuit and gives [$Failed]() with a message:

```wl
QiskitCircuit["not openqasm source"]
```

<!-- => QiskitCircuit::badarg message; $Failed -->
