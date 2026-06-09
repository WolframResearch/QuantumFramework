---
Template: Symbol
Name: IBMJob
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/IBMJob
Keywords: [IBM Quantum, QPU, job, result, counts, measurement, SamplerV2, quantum seconds, ServiceConnect, hardware]
SeeAlso: [IBMJobSubmit, QuantumMeasurement, QuantumCircuitOperator, QuantumQASM]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: [SendingQueriesToIBMQPUs]
---

## Usage

<code>[IBMJob]()[*assoc*]</code> represents a job submitted to an IBM Quantum processing unit, carrying its identifier, backend and the raw service responses, with curated accessors layered on top.

<code>[IBMJob]()[*assoc*][*"prop"*]</code> gives the value of the property *"prop"*, such as `"Status"`, `"Counts"` or `"QuantumSeconds"`.

<code>[IBMJob]()[*assoc*][]</code> gives the primary result, the [QuantumMeasurement]() reconstructed from the hardware counts.

## Details & Options

- An `IBMJob` is normally produced by [IBMJobSubmit]() or by the <code>[QuantumCircuitOperator]()[*circ*][Method -> {"Qiskit", "Provider" -> "IBMProvider", …}]</code> hardware path, not constructed by hand.
- The object is a **lossless** carrier: the raw service responses are kept verbatim under `"Raw"`, and every accessor reads from there, so no field the service returns is discarded.
- Accessors **fetch lazily**: a property reads from the stored snapshot when present, otherwise queries the service object. With no reachable connection an accessor gives `Missing["NoConnection"]`; before the job is `"Completed"` the result accessors give `Missing["JobNotComplete", `*status*`]`.
- `"Status"` is one of `"Queued"`, `"Validating"`, `"Initializing"`, `"Running"`, `"Completed"`, `"Cancelled"` or `"Failed"`; the terminal statuses are `"Completed"`, `"Cancelled"` and `"Failed"`.
- The hardware counts decode into the Wolfram Language qubit order (ascending qubit) using the `"MeasuredQubits"` map (classical bit to original qubit) captured at submission time, so a permuted or partial measurement and any backend layout decode correctly.
- `IBMJob[`*assoc*`]["Refresh"]` re-queries the service and gives a **new** `IBMJob` with an updated `"Raw"` snapshot; once the job is `"Completed"` it also pulls the results and metrics. `IBMJob[`*assoc*`]["Cancel"]` requests cancellation.
- An unknown property name gives a [Failure]() listing the valid properties (see Possible Issues).
- The summary box shows the status (with a colored indicator), the backend and the job id; once `"Completed"` it also shows the qubit count, shot count and quantum-processing time.

The full property list (`IBMJob[`*assoc*`]["Properties"]`):

| Property | Result |
|---|---|
| `"ID"` | the job identifier string |
| `"Backend"` | the backend name |
| `"ServiceObject"` | the `ServiceObject` the job is bound to |
| `"Status"` | the current job status |
| `"UserID"` | the submitting user id |
| `"ProgramID"` | the primitive program, `"sampler"` or `"estimator"` |
| `"Cost"` | the reserved cost, in seconds |
| `"EstimatedRunTime"` | the estimated running time, in seconds |
| `"QuantumSeconds"` | the quantum-processing time used, in seconds |
| `"SubmittedCircuit"` | the submitted ISA OpenQASM program |
| `"Options"` | the `params.options` block sent with the job |
| `"Timestamps"` | an Association of lifecycle `DateObject`s |
| `"Duration"` | the running-to-finished wall-clock duration |
| `"Samples"` | the per-shot classical-register readouts |
| `"MeasuredQubits"` | the classical-bit to original-qubit decode map |
| `"NumBits"`, `"Qubits"` | the number of measured qubits |
| `"Shots"` | the total number of shots |
| `"Measurement"` | the [QuantumMeasurement]() of the outcome |
| `"Counts"` | exact integer counts keyed by the bit-list outcome |
| `"Probabilities"` | the outcome probabilities |
| `"ExecutedCircuit"` | the server-side transpiled circuit (see Possible Issues) |
| `"ExecutionSpans"` | the execution time spans reported by the service |
| `"Refresh"` | a new `IBMJob` with a re-queried snapshot |
| `"Cancel"` | requests cancellation of the job |
| `"Raw"` | the raw service responses; `["Raw", `*"sect"*`]` gives one section |
| `"Properties"` | the list of property names |

## Basic Examples

A handle returned by [IBMJobSubmit]() displays as a status summary box:

```wl
#| eval: false
job = IBMJobSubmit[QuantumCircuitOperator[{"GHZ", "Fourier"[3], {1, 2, 3}}], "ibm_fez"]
```

<!-- => IBMJob[<| Status: Queued, Backend: ibm_fez, Job ID: d8jhf29e8nrc73bjk610 |>] -->

---

Check the status, refreshing the snapshot from the service:

```wl
#| eval: false
job = job["Refresh"]; job["Status"]
```

<!-- => Completed -->

---

The default output is the measurement reconstructed from the hardware counts:

```wl
#| eval: false
job[]
```

<!-- => QuantumMeasurement[3 qubits, 4096 shots] -->

## Scope

### Results

The outcome probabilities, computed from the measured frequencies and aligned to the Wolfram Language qubit order:

```wl
#| eval: false
job["Probabilities"]
```

<!-- => <|{0,0,0}->0.255, {0,0,1}->0.197, {0,1,0}->0.149, {0,1,1}->0.054, {1,0,0}->0.019, {1,0,1}->0.035, {1,1,0}->0.120, {1,1,1}->0.171|> -->

---

Exact integer counts keyed by the bit-list outcome:

```wl
#| eval: false
job["Counts"]
```

<!-- => <|{0,0,0}->1043, {0,0,1}->807, ... |> -->

---

The number of measured qubits and the total number of shots:

```wl
#| eval: false
{job["NumBits"], job["Shots"]}
```

<!-- => {3, 4096} -->

---

The classical-bit to original-qubit map used to decode the samples into qubit order:

```wl
#| eval: false
job["MeasuredQubits"]
```

<!-- => {0, 1, 2} -->

---

The raw per-shot classical-register readouts:

```wl
#| eval: false
Short[job["Samples"], 2]
```

<!-- => {"0x5", "0x0", "0x7", ...} -->

### Job metadata

The job identifier and backend:

```wl
#| eval: false
{job["ID"], job["Backend"]}
```

<!-- => {d8jhf29e8nrc73bjk610, ibm_fez} -->

---

The quantum-processing time consumed, as a `Quantity` in seconds:

```wl
#| eval: false
job["QuantumSeconds"]
```

<!-- => Quantity[3.2, "Seconds"] -->

---

The lifecycle timestamps and the running-to-finished duration:

```wl
#| eval: false
{job["Timestamps"], job["Duration"]}
```

<!-- => {<|created->DateObject[...], running->DateObject[...], finished->DateObject[...]|>, Quantity[..., "Seconds"]} -->

---

The submitted ISA OpenQASM program, as actually sent to the QPU:

```wl
#| eval: false
StringTake[job["SubmittedCircuit"], 40]
```

<!-- => OPENQASM 2.0; include "qelib1.inc"; qreg q -->

### Lifecycle

`"Refresh"` re-queries the service and gives a new handle with an updated snapshot:

```wl
#| eval: false
job["Refresh"]["Status"]
```

<!-- => Completed -->

---

`"Cancel"` requests cancellation of a queued or running job:

```wl
#| eval: false
job["Cancel"]
```

<!-- => <|...|> -->

### Raw responses

The complete raw service responses are kept verbatim; ask for a single section by name:

```wl
#| eval: false
Keys @ job["Raw"]
```

<!-- => {Details, Results, Metrics} -->

```wl
#| eval: false
job["Raw", "Metrics"]
```

<!-- => <|"usage" -> <|"quantum_seconds" -> 3.2, ...|>, "timestamps" -> <|...|>|> -->

## Possible Issues

An unknown property gives a [Failure]() listing the valid property names:

```wl
#| eval: false
job["Frobnicate"]
```

<!-- => Failure["IBMJob", <|"MessageTemplate" -> "Unknown property `1`. Use one of: `2`.", "MessageParameters" -> {"Frobnicate", {"ID", "Backend", ...}}|>] -->

---

`"ExecutedCircuit"` downloads the server-side transpiled circuit, which IBM stores only when the Runtime transpiles server-side. A circuit submitted already in ISA form to the V2 primitives stores none, so the request gives a [Failure]() rather than a circuit:

```wl
#| eval: false
job["ExecutedCircuit"]
```

<!-- => Failure["IBMExecutedCircuit", <|"MessageTemplate" -> "IBM stored no transpiled circuit for job `1`. ..."|>] -->

---

When the transpiled-circuit object does exist, IBM serves it from a Cloud Object Storage direct endpoint (host starting `s3.direct`) reachable only from inside IBM Cloud, so the download from elsewhere gives a [Failure]() instead of timing out silently:

```wl
#| eval: false
job["ExecutedCircuit"]
```

<!-- => Failure["IBMExecutedCircuit", <|"MessageTemplate" -> "Could not reach `1` to download the transpiled circuit. ..."|>] -->
