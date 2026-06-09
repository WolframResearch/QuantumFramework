---
Template: Symbol
Name: IBMJobSubmit
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/IBMJobSubmit
Keywords: [IBM Quantum, QPU, job, sampler, estimator, SamplerV2, OpenQASM, ServiceConnect, hardware, asynchronous]
SeeAlso: [IBMJob, QuantumQASM, QuantumCircuitOperator, QuantumMeasurement]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: [SendingQueriesToIBMQPUs]
---

## Usage

<code>[IBMJobSubmit]()[*circ*]</code> submits the quantum circuit *circ* to the default backend of the active IBM Quantum Platform connection and gives an [IBMJob]() handle.

<code>[IBMJobSubmit]()[*circ*, *"backend"*]</code> submits *circ* to the named backend.

<code>[IBMJobSubmit]()[*circ*, *backend*, *opts*]</code> submits with the options *opts*, such as the primitive, the number of shots and the pass-through primitive option tree.

## Details & Options

- *circ* is a [QuantumCircuitOperator](). It is transpiled to an ISA OpenQASM program against the chosen backend through [QuantumQASM]() (the `"Provider" -> "IBMProvider"` path) before submission.
- `IBMJobSubmit` requires an **active** service connection: it consults <code>ServiceFramework`GetDefaultServiceObject["IBMQuantumPlatform"]</code> and never creates a connection or prompts for credentials. Run <code>[ServiceConnect]()["IBMQuantumPlatform"]</code> first (see the [Sending Queries to IBM QPUs]() tech note). With no active connection it returns a [Failure]() and never submits.
- `IBMJobSubmit` returns immediately with an asynchronous [IBMJob]() handle whose `"Status"` is whatever the service reports (typically `"Queued"`). Query the handle later, or pass `"Wait" -> True` to block until the job reaches a terminal status (`"Completed"`, `"Cancelled"` or `"Failed"`).
- The submission carries the per-classical-bit to original-qubit map captured at transpile time into the [IBMJob]() as `"MeasuredQubits"`, so the returned samples decode into ascending-qubit order regardless of how the backend lays out and routes the circuit.
- *backend* defaults to `Automatic`, the first backend in `so["Backends"]` for the active connection *so*; give a string such as `"ibm_fez"` to target a specific QPU.
- The job is submitted through the connection's low-level `"RawJobRun"` request, which gives full control over the `params.options` block; the user's `"PrimitiveOptions"` tree is deep-merged over the defaults `<|"default_shots" -> shots, "dynamical_decoupling" -> <|"enable" -> True|>|>`.
- Option keys in `"PrimitiveOptions"` are written in CamelCase and converted to the service's snake_case keys recursively (`"DynamicalDecoupling"` to `dynamical_decoupling`), so the whole IBM Runtime options schema is reachable without a curated key list.

| Option | Default value | Description |
|---|---|---|
| `"Primitive"` | `"sampler"` | the IBM Runtime primitive, `"sampler"` or `"estimator"` |
| `"Shots"` | 4096 | the number of shots (`default_shots`) |
| `"Observable"` | `"Z"` | the observable, used only when `"Primitive"` is `"estimator"` |
| `"Version"` | 2 | the OpenQASM version of the submitted ISA circuit, 2 or 3 |
| `"Wait"` | False | whether to block until the job reaches a terminal status |
| `"PrimitiveOptions"` | `<\|\|>` | an Association of IBM Runtime primitive options, deep-merged over the defaults |

## Basic Examples

With an active IBM Quantum Platform connection, submit a circuit and get back a job handle (the call returns immediately, while the job sits in the queue):

```wl
#| eval: false
job = IBMJobSubmit[QuantumCircuitOperator[{"GHZ", "Fourier"[3], {1, 2, 3}}], "ibm_fez"]
```

<!-- => IBMJob[<| Status: Queued, Backend: ibm_fez, Job ID: d8jhf29e8nrc73bjk610 |>] -->

---

Once the job completes, the handle's default output is the measurement reconstructed from the hardware counts, in the Wolfram Language qubit order:

```wl
#| eval: false
job["Refresh"][]
```

<!-- => QuantumMeasurement[3 qubits, 4096 shots] -->

## Scope

Submit to the active connection's default backend by omitting the backend argument:

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}]]
```

<!-- => IBMJob[<| Status: Queued, Backend: ibm_fez, Job ID: ... |>] -->

---

Block until the job finishes with `"Wait" -> True`; the returned handle is already `"Completed"`:

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez", "Wait" -> True]["Status"]
```

<!-- => Completed -->

---

Set the number of shots:

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez", "Shots" -> 1024]
```

<!-- => IBMJob[<| Status: Queued, Backend: ibm_fez, Job ID: ... |>] -->

---

Reach into the IBM Runtime options schema through `"PrimitiveOptions"`, given in CamelCase and converted to the service's snake_case keys (here turning off dynamical decoupling and selecting a resilience level):

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez",
  "PrimitiveOptions" -> <|"DynamicalDecoupling" -> <|"Enable" -> False|>, "ResilienceLevel" -> 1|>]
```

<!-- => IBMJob[<| Status: Queued, Backend: ibm_fez, Job ID: ... |>] -->

---

Submit an `"estimator"` job with an observable instead of a sampler job:

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez", "Primitive" -> "estimator", "Observable" -> "ZZZ"]
```

<!-- => IBMJob[<| Status: Queued, ProgramID: estimator, Backend: ibm_fez, Job ID: ... |>] -->

## Options

### "Shots"

The number of shots is forwarded as the primitive's `default_shots`:

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez", "Shots" -> 8192]["Shots"]
```

<!-- => 8192 -->

### "Version"

Emit the submitted ISA circuit as OpenQASM 3 instead of the default OpenQASM 2:

```wl
#| eval: false
StringStartsQ[IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez", "Version" -> 3]["SubmittedCircuit"], "OPENQASM 3.0"]
```

<!-- => True -->

### "Wait"

`"Wait" -> True` polls the service until the job is terminal before returning the handle:

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez", "Wait" -> True]["QuantumSeconds"]
```

<!-- => Quantity[..., "Seconds"] -->

### "Primitive"

Choose between the `"sampler"` and `"estimator"` primitives; `"Observable"` applies to the estimator:

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez", "Primitive" -> "estimator", "Observable" -> "ZZZ"]["ProgramID"]
```

<!-- => estimator -->

## Possible Issues

With no active IBM Quantum Platform connection, `IBMJobSubmit` returns a [Failure]() and submits nothing:

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez"]
```

<!-- => Failure["IBMJobSubmit", <|"MessageTemplate" -> "No active IBM Quantum connection. Run ServiceConnect[\"`1`\"] first.", "MessageParameters" -> {"IBMQuantumPlatform"}|>] -->

---

When the backend is `Automatic` and the active connection exposes no backends, the default backend cannot be determined:

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}]]
```

<!-- => Failure["IBMJobSubmit", <|"MessageTemplate" -> "Could not determine a default backend."|>] -->

---

A circuit carrying any non-qubit (higher-dimensional) system has no OpenQASM representation, so it is rejected before submission with the same [Failure]() [QuantumQASM]() gives:

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"H"[3] -> 1, {1}}], "ibm_fez"]
```

<!-- => Failure["QuantumQASM", <|... "NonQubitDimensions" -> {3}|>] : QuantumQASM supports qubit (2-dimensional) systems only. -->

---

If the service accepts the request but returns no job id, the response is surfaced in the [Failure]():

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez"]
```

<!-- => Failure["IBMJobSubmit", <|"MessageTemplate" -> "Job submission did not return a job id.", "Response" -> <|...|>|>] -->
