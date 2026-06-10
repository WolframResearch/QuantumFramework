---
Template: Symbol
Name: IBMJob
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/IBMJob
Keywords: [IBM Quantum, QPU, job, result, counts, measurement, expectation value, SamplerV2, EstimatorV2, quantum seconds, ServiceConnect, hardware]
SeeAlso: [IBMJobSubmit, QuantumMeasurement, QuantumCircuitOperator, QuantumQASM]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: [SendingQueriesToIBMQPUs]
---

## Usage

<code>[IBMJob]()[*assoc*]</code> represents a job submitted to an IBM Quantum processing unit, carrying its identifier, backend and the raw service responses, with curated accessors layered on top.

<code>[IBMJob]()[*assoc*][*"prop"*]</code> gives the value of the property *"prop"*, such as `"Status"`, `"Counts"` or `"QuantumSeconds"`.

<code>[IBMJob]()[*assoc*][]</code> gives the primary result: the [QuantumMeasurement]() reconstructed from the hardware counts for a sampler job, or the observable's expectation value for an estimator job.

## Details & Options

- An `IBMJob` is normally produced by [IBMJobSubmit]() or by the <code>[QuantumCircuitOperator]()[*circ*][Method -> {"Qiskit", "Provider" -> "IBMProvider", …}]</code> hardware path, not constructed by hand.
- The object is a **lossless** carrier: the raw service responses are kept verbatim under `"Raw"`, and every accessor reads from there, so no field the service returns is discarded.
- Accessors **fetch lazily**: a property reads from the stored snapshot when present, otherwise queries the service object. With no reachable connection an accessor gives `Missing["NoConnection"]`; before the job is `"Completed"` the result accessors give `Missing["JobNotComplete", `*status*`]`.
- `"Status"` is one of `"Queued"`, `"Validating"`, `"Initializing"`, `"Running"`, `"Completed"`, `"Cancelled"` or `"Failed"`; the terminal statuses are `"Completed"`, `"Cancelled"` and `"Failed"`.
- The hardware counts decode into the Wolfram Language qubit order (ascending qubit) using the `"MeasuredQubits"` map (classical bit to original qubit) captured at submission time, so a permuted or partial measurement and any backend layout decode correctly.
- The handle exposes two kinds of names. **Properties** (`IBMJob[`*assoc*`]["Properties"]`) are local accessors that read the cached snapshot and never touch the network, so a bulk map over them is fast and has no side effects. **Actions** (`IBMJob[`*assoc*`]["Actions"]`) each require a live connection and are kept out of the property list, so that map never re-queries, cancels or blocks on a remote fetch; each is still callable by name. `IBMJob[`*assoc*`]["Refresh"]` re-queries the service and gives a **new** `IBMJob` with an updated `"Raw"` snapshot (once `"Completed"` it also pulls the results and metrics, and for an estimator job caches the qiskit-decoded expectation values); `IBMJob[`*assoc*`]["Cancel"]` requests cancellation; `IBMJob[`*assoc*`]["ExecutedCircuit"]` downloads the server-side transpiled circuit (see Possible Issues).
- The result accessors are primitive-specific. For a **sampler** job the live accessors are `"Counts"`, `"Measurement"`, `"Probabilities"`, `"Samples"`, `"NumBits"` and `"Shots"`; for an **estimator** job they are `"ExpectationValue"`, `"ExpectationValues"` and `"StandardErrors"`. An accessor that does not apply to the job's primitive gives `Missing["NotApplicable", `*primitive*`]`.
- An unknown property name gives a [Failure]() listing the valid properties (see Possible Issues).
- The summary box shows the status (with a colored indicator), the backend and the job id; once `"Completed"` a sampler job also shows the qubit count, shot count and quantum-processing time, and an estimator job shows the expectation value and quantum-processing time.

The data properties (`IBMJob[`*assoc*`]["Properties"]`), each read from the cached snapshot with no network access:

| Property | Result |
|---|---|
| `"ID"` | the job identifier string |
| `"Backend"` | the backend name |
| `"ServiceObject"` | the `ServiceObject` the job is bound to |
| `"Status"` | the current job status |
| `"Primitive"` | the IBM Runtime primitive the job was submitted with, `"sampler"` or `"estimator"` |
| `"UserID"` | the submitting user id |
| `"ProgramID"` | the primitive program, `"sampler"` or `"estimator"` |
| `"Cost"` | the reserved cost, in seconds |
| `"EstimatedRunTime"` | the estimated running time, in seconds |
| `"QuantumSeconds"` | the quantum-processing time used, in seconds |
| `"SubmittedCircuit"` | the submitted circuit, as stored in the job's `params` |
| `"Options"` | the `params.options` block sent with the job |
| `"Timestamps"` | an Association of lifecycle `DateObject`s |
| `"Duration"` | the running-to-finished wall-clock duration |
| `"Samples"` | the per-shot outcomes, each the integer value of the classical register for that shot |
| `"MeasuredQubits"` | the classical-bit to original-qubit decode map |
| `"NumBits"`, `"Qubits"` | the number of measured qubits |
| `"Shots"` | the total number of shots |
| `"Measurement"` | the [QuantumMeasurement]() of the outcome (sampler) |
| `"Counts"` | exact integer counts keyed by the bit-list outcome (sampler) |
| `"Probabilities"` | the outcome probabilities (sampler) |
| `"ExpectationValue"` | the observable's expectation value (estimator); a single value for one observable, else a list |
| `"ExpectationValues"` | the list of expectation values, one per observable (estimator) |
| `"StandardErrors"` | the per-observable standard error of the expectation value (estimator) |
| `"ExecutionSpans"` | the execution time spans reported by the service |
| `"Raw"` | the raw service responses; `["Raw", `*"sect"*`]` gives one section |
| `"Properties"` | the list of data property names |
| `"Actions"` | the list of action names |

The actions (`IBMJob[`*assoc*`]["Actions"]`), each requiring a live connection and so kept out of `"Properties"`:

| Action | Result |
|---|---|
| `"Refresh"` | re-queries the service and gives a new `IBMJob` with an updated snapshot |
| `"Cancel"` | requests cancellation of the job |
| `"ExecutedCircuit"` | downloads the server-side transpiled circuit (see Possible Issues) |

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

The raw per-shot outcomes, each the integer value of the classical register for that shot:

```wl
#| eval: false
Short[job["Samples"], 2]
```

<!-- => {5, 0, 7, 0, 4, 7, 2, ...} -->

### Estimator results

For a job submitted with `"Primitive" -> "estimator"`, the result is the observable's expectation value and its standard error rather than a measurement; the sampler accessors report `Missing["NotApplicable", "estimator"]`:

```wl
#| eval: false
estJob = IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez",
  "Primitive" -> "estimator", "Observable" -> "ZZZ", "Wait" -> True];
{estJob["ExpectationValue"], estJob["StandardErrors"]}
```

<!-- => {0.03..., {0.01...}} -->

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

The submitted circuit, as stored in the job's `params` (qiskit serializes it before submission, so it is not OpenQASM text):

```wl
#| eval: false
job["SubmittedCircuit"]
```

<!-- => (the submitted circuit, as serialized in the job's params) -->

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

### Properties and actions

`"Properties"` lists the local accessors; mapping over them reads only the cached snapshot, with no network access:

```wl
#| eval: false
AssociationMap[job, job["Properties"]]
```

<!-- => <|ID -> d8jhf29e8nrc73bjk610, Backend -> ibm_fez, Status -> Completed, ...|> -->

---

`"Actions"` lists the names that need a live connection, kept out of `"Properties"` so a bulk map never re-queries, cancels or blocks on a remote fetch:

```wl
#| eval: false
job["Actions"]
```

<!-- => {Refresh, Cancel, ExecutedCircuit} -->

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

An unknown property gives a [Failure]() listing the valid data properties and actions:

```wl
#| eval: false
job["Frobnicate"]
```

<!-- => Failure["IBMJob", <|"MessageTemplate" -> "Unknown property `1`. Data properties: `2`. Actions (need a live connection): `3`.", "MessageParameters" -> {"Frobnicate", {"ID", "Backend", ...}, {"Refresh", "Cancel", "ExecutedCircuit"}}|>] -->

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
