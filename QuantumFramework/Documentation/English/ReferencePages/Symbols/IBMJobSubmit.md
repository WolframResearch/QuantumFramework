---
Template: Symbol
Name: IBMJobSubmit
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/IBMJobSubmit
Keywords: [IBM Quantum, QPU, job, sampler, estimator, SamplerV2, EstimatorV2, observable, ServiceConnect, hardware, asynchronous]
SeeAlso: [IBMJob, QuantumQASM, QuantumCircuitOperator, QuantumMeasurement]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: [QPUServiceConnect]
---

## Usage

<code>[IBMJobSubmit]()[*circ*]</code> submits the quantum circuit *circ* to the default backend of the active IBM Quantum Platform connection and gives an [IBMJob]() handle.

<code>[IBMJobSubmit]()[*circ*, *"backend"*]</code> submits *circ* to the named backend.

<code>[IBMJobSubmit]()[*circ*, *backend*, *opts*]</code> submits with the options *opts*, such as the primitive, the number of shots and the pass-through primitive option tree.

## Details & Options

- *circ* is a [QuantumCircuitOperator](). It is transpiled against the chosen backend and submitted through qiskit's own `SamplerV2` / `EstimatorV2` primitive, so qiskit owns every wire format: the circuit serialization, the estimator observable's layout, and the primitive options schema. The job is launched without blocking on its result, so the handle returns while the job is still queued.
- `IBMJobSubmit` requires an **active** service connection: it consults <code>ServiceFramework`GetDefaultServiceObject["IBMQuantumPlatform"]</code> and never creates a connection or prompts for credentials. Run <code>[ServiceConnect]()["IBMQuantumPlatform"]</code> first (see the [Sending Queries to IBM QPUs]() tech note). With no active connection it returns a [Failure]() and never submits.
- `IBMJobSubmit` returns immediately with an asynchronous [IBMJob]() handle whose `"Status"` is whatever the service reports (typically `"Queued"`). Query the handle later, or pass `"Wait" -> True` to block until the job reaches a terminal status (`"Completed"`, `"Cancelled"` or `"Failed"`).
- For a sampler job the submission carries the per-classical-bit to original-qubit map captured at transpile time into the [IBMJob]() as `"MeasuredQubits"`, so the returned samples decode into ascending-qubit order regardless of how the backend lays out and routes the circuit.
- *backend* defaults to `Automatic`, the first backend in `so["Backends"]` for the active connection *so*; give a string such as `"ibm_fez"` to target a specific QPU.
- For the `"estimator"` primitive the `"Observable"` is built into a qiskit `SparsePauliOp` and mapped onto the transpiled circuit's layout, which both reorders and **widens** it to the backend's physical register. Submitting the bare logical-width observable is what an IBM Runtime estimator rejects (error 1501, "the number of qubits of the circuit does not match the number of qubits of the observable").
- The `"PrimitiveOptions"` tree is applied onto the qiskit primitive's own options object over the defaults `default_shots -> shots` and `dynamical_decoupling.enable -> True`, so qiskit validates the option schema. Keys are written in CamelCase and converted to qiskit's snake_case attributes recursively (`"DynamicalDecoupling"` to `dynamical_decoupling`), so the whole options schema is reachable without a curated key list.

| Option | Default value | Description |
|---|---|---|
| `"Primitive"` | `"sampler"` | the IBM Runtime primitive, `"sampler"` or `"estimator"` |
| `"Shots"` | 4096 | the number of shots (`default_shots`) |
| `"Observable"` | `"Z"` | the estimator observable: a Pauli string (`"ZZZ"`), a list of Pauli strings, or a list of `{pauliString, coefficient}` pairs; used only when `"Primitive"` is `"estimator"` |
| `"Wait"` | False | whether to block until the job reaches a terminal status |
| `"PrimitiveOptions"` | `<\|\|>` | an Association of IBM Runtime primitive options, applied onto the qiskit primitive's options object over the defaults |
| `"OptimizationLevel"` | Automatic | the transpiler optimization level (0–3) for the submitted circuit. Submission always transpiles against the backend's own `Target` (per-instruction error and duration), so layout and routing are error-aware; this controls how aggressively. Automatic uses qiskit's preset default. |

## Basic Examples

With an active IBM Quantum Platform connection, submit a circuit and get back a job handle (the call returns immediately, while the job sits in the queue):

```wl
#| eval: false
job = IBMJobSubmit[QuantumCircuitOperator[{"GHZ", "Fourier"[3], {1, 2, 3}}], "ibm_fez"]
```

<!-- => IBMJob[<| Status: Queued, Backend: ibm_fez, Job ID: d8kcdur2d42s73ca1r40 |>] -->

```wl
#| eval: false
job["Status"]
```

<!-- => "Queued" -->

---

Once the job completes, the handle's default output is the measurement reconstructed from the hardware counts, in the Wolfram Language qubit order:

```wl
#| eval: false
job["Refresh"]
```

<!-- => IBMJob[<| Status: Completed, Backend: ibm_fez, Job ID: d8kcdur2d42s73ca1r40 |>] -->

Get the measurement results from the QPU:

```wl
#| eval: false
qpu = job["Refresh"][]
```

<!-- => QuantumMeasurement[3 qubits, 4096 shots] -->

Compute the corresponding measurement results in the Wolfram Language:

```wl
#| eval: false
wl = QuantumCircuitOperator[{"GHZ", "Fourier"[3], {1, 2, 3}}][];
```

Compare the exact Wolfram Language probabilities with the QPU estimates:

```wl
#| eval: false
BarChart[Transpose[Values @* KeySort /@ {wl["Probabilities"], qpu["Probabilities"]}],
  AspectRatio -> 1/2, Frame -> True, ChartLegends -> {"Exact WL", "QPU"},
  PlotLabels -> Keys[wl["Probabilities"]]]
```

<!-- => (bar chart: exact-WL vs QPU probability per basis state) -->

## Scope

Submit to the active connection's default backend by omitting the backend argument:

```wl
#| eval: false
job = IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}]]
```

<!-- => IBMJob[<| Status: Queued, Backend: ibm_fez, Job ID: d8kcgfg32u0s73f8jbu0 |>] -->

Once the job completes, plot the probabilities the QPU measured:

```wl
#| eval: false
job["Refresh"][]["ProbabilityPlot", AspectRatio -> 1/3]
```

<!-- => (probability bar chart over the measured basis states) -->

---

Block until the job finishes with `"Wait" -> True`; the returned handle is already `"Completed"`:

```wl
#| eval: false
job = IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez", "Wait" -> True]
```

<!-- => IBMJob[<| Status: Completed, Backend: ibm_fez, Job ID: d8kcguj2d42s73ca1u3g |>] -->

```wl
#| eval: false
job["Status"]
```

<!-- => "Completed" -->

---

Set the number of shots, then read the hardware counts back:

```wl
#| eval: false
job = IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez", "Shots" -> 1024]
```

<!-- => IBMJob[<| Status: Queued, Backend: ibm_fez, Job ID: d8kch2032u0s73f8jch0 |>] -->

```wl
#| eval: false
job["Refresh"]["Counts"]
```

<!-- => <|{1, 1, 1} -> 454, {1, 0, 0} -> 54, {0, 1, 1} -> 68, {0, 0, 0} -> 420, {1, 0, 1} -> 9, {0, 1, 0} -> 3, {0, 0, 1} -> 6, {1, 1, 0} -> 10|> -->

---

Reach into the IBM Runtime options schema through `"PrimitiveOptions"`, given in CamelCase and converted to the service's snake_case keys (here turning off dynamical decoupling):

```wl
#| eval: false
job = IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez",
  "PrimitiveOptions" -> <|"DynamicalDecoupling" -> <|"Enable" -> False|>|>]
```

<!-- => IBMJob[<| Status: Queued, Backend: ibm_fez, Job ID: d8kci4rnn5bs738qmb50 |>] -->

```wl
#| eval: false
N[job["Refresh"][]["Probabilities"]]
```

<!-- => <||000⟩ -> 0.451172, |001⟩ -> 0.003174, |010⟩ -> 0.004150, |011⟩ -> 0.070313, |100⟩ -> 0.057373, |101⟩ -> 0.009277, |110⟩ -> 0.003174, |111⟩ -> 0.401367|> -->

---

Submit an `"estimator"` job with an observable instead of a sampler job. The observable is laid out onto the transpiled circuit, so a 3-qubit `"ZZZ"` is accepted against a many-qubit backend:

```wl
#| eval: false
job = IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez", "Primitive" -> "estimator", "Observable" -> "ZZZ"]
```

<!-- => IBMJob[<| Status: Queued, Backend: ibm_fez, Job ID: d8kcilr2d42s73ca1vo0 |>] -->

```wl
#| eval: false
job["Refresh"]
```

<!-- => IBMJob[<| Status: Completed, Backend: ibm_fez, Job ID: d8kcilr2d42s73ca1vo0 |>] -->

Once an estimator job completes, its default output and `"ExpectationValue"` are the qiskit-decoded expectation value of the observable, with `"StandardErrors"` giving the per-observable standard error:

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez", "Primitive" -> "estimator", "Observable" -> "ZZZ", "Wait" -> True]["ExpectationValue"]
```

<!-- => 0.055201698513800426 -->

The exact Wolfram Language value of the same expectation is zero, so the hardware estimate sits near zero up to shot noise:

```wl
#| eval: false
With[{qs = QuantumCircuitOperator[{"GHZ"}][]}, qs["Dagger"][QuantumOperator["ZZZ"][qs]]["Scalar"]]
```

<!-- => 0 -->

## Options

### "Shots"

The number of shots is forwarded as the primitive's `default_shots`:

```wl
#| eval: false
job = IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez", "Shots" -> 8192]
```

<!-- => IBMJob[<| Status: Queued, Backend: ibm_fez, Job ID: d8kei2o32u0s73f8lml0 |>] -->

```wl
#| eval: false
job["Refresh"]["Shots"]
```

<!-- => 8192 -->

### "Wait"

`"Wait" -> True` polls the service until the job is terminal before returning the handle:

```wl
#| eval: false
job = IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez", "Wait" -> True]
```

<!-- => IBMJob[<| Status: Completed, Backend: ibm_fez, Job ID: d8keibjqv2lc73856j40 |>] -->

```wl
#| eval: false
job["QuantumSeconds"]
```

<!-- => Quantity[3, "Seconds"] -->

### "Primitive"

Choose between the `"sampler"` and `"estimator"` primitives; `"Observable"` applies to the estimator:

```wl
#| eval: false
job = IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez", "Primitive" -> "estimator", "Observable" -> "ZZZ"]
```

<!-- => IBMJob[<| Status: Queued, Backend: ibm_fez, Job ID: d8keikjnn5bs738qol20 |>] -->

```wl
#| eval: false
job["ProgramID"]
```

<!-- => "estimator" -->

### "OptimizationLevel"

Submission always transpiles the circuit against the backend's own `Target`, which carries per-instruction error and duration, so layout and routing are **error-aware**: the circuit is placed on the backend's lowest-error qubits and pairs automatically. `"OptimizationLevel"` (0–3) controls how hard the transpiler works at that; `Automatic` uses qiskit's preset default. Raise it for a deeper search on a noisy device:

```wl
#| eval: false
job = IBMJobSubmit[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}], "ibm_fez", "OptimizationLevel" -> 3]
```

<!-- => IBMJob[<| Status: Queued, Backend: ibm_fez, Job ID: ... |>] -->

The fidelity payoff of this error-aware placement is quantified on the [QiskitTarget]() page (reading the same backend's error model, the chosen embedding's total expected error was about 4x lower than a connectivity-only placement). To inspect or pre-select that placement before submitting, build the target explicitly with `QiskitTarget["FromBackend" -> "ibm_fez", "Provider" -> "IBMProvider"]` and transpile with [QuantumQASM]().

## Possible Issues

With no active IBM Quantum Platform connection, `IBMJobSubmit` returns a [Failure]() and submits nothing:

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez"]
```

<!-- => Failure["IBMJobSubmit", <|"MessageTemplate" -> "No active IBM Quantum connection. Run ServiceConnect[\"`1`\"] first.", "MessageParameters" -> {"IBMQuantumPlatform"}|>] -->

---

When the backend is `Automatic`, the first backend the active connection exposes is used. If the connection exposes none, the default backend cannot be determined and the submission fails before any request is sent:

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}]]
```

<!-- => Failure["IBMJobSubmit", <|"MessageTemplate" -> "Could not determine a default backend."|>] -->

---

A circuit carrying any non-qubit (higher-dimensional) system is not a qubit circuit a QPU can run, so it is rejected before submission with the same [Failure]() [QuantumQASM]() gives:

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"H"[3], "M"[QuantumBasis[3]]}], "ibm_fez"]
```

<!-- => Failure["QuantumQASM", <|"MessageTemplate" -> "QuantumQASM supports qubit (2-dimensional) systems only; the given `1` has non-qubit qudit dimension(s) `2`. OpenQASM has no representation for higher-dimensional qudits.", "MessageParameters" -> {"circuit", {3}}, "NonQubitDimensions" -> {3}|>] -->

---

If the service accepts the request but the response carries no job id, the raw response is surfaced under the `"Response"` key of the [Failure]():

```wl
#| eval: false
IBMJobSubmit[QuantumCircuitOperator[{"GHZ", {1, 2, 3}}], "ibm_fez"]
```

<!-- => Failure["IBMJobSubmit", <|"MessageTemplate" -> "Job submission did not return a job id.", "Response" -> <|"status" -> "accepted"|>|>] -->
