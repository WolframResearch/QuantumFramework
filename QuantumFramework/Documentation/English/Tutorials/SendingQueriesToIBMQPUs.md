---
Template: TechNote
Name: SendingQueriesToIBMQPUs
Title: Sending Queries to IBM QPUs
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/SendingQueriesToIBMQPUs
Keywords: [IBM Quantum, QPU, OpenQASM, Qiskit, ServiceConnect, IBMJobSubmit, IBMJob, sampler, SamplerV2, transpile, quantum hardware]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: [QPUServiceConnect, GettingStarted]
---

In the [Wolfram Quantum Framework](https://www.wolfram.com/quantum-computation-framework/) a quantum circuit is a single symbolic object: you can simulate it exactly and send that same object to a real [IBM Quantum](https://www.ibm.com/quantum) processing unit (QPU). This tutorial walks the full path: connect to the IBM Quantum Platform, build and check a circuit, inspect it as [OpenQASM](https://openqasm.com/) (the portable text format QPUs read), submit it to hardware with <code>[IBMJobSubmit]()</code>, and hold the hardware counts up against the exact result.

Building and checking run locally with no account. The hardware cells need an IBM Quantum account and are marked so you can read the tutorial without running them. For help, write to quantum@wolfram.com.

```wl
<< Wolfram`QuantumFramework`
```

## Connecting to the IBM Quantum Platform

Create an account at [quantum.cloud.ibm.com](https://quantum.cloud.ibm.com/). From the dashboard you need two things: an **IAM API key** (the long secret string shown once when you create the key, not the key's display name) and the **instance CRN** (a `crn:v1:…` string identifying your instance).

The API key is a secret and is stored in your encrypted system credential store. The CRN is read from a local symbol that the connection consults when it builds the request headers, so set both before connecting:

```wl
#| eval: false
api = "XXX";
LocalSymbol["ibm_crn"] = "XXX";
```

Loading the Wolfram Quantum Framework registers the IBM Quantum Platform service automatically, so [ServiceConnect]() recognizes the `"IBMQuantumPlatform"` name directly:

```wl
#| eval: false
ibm = ServiceConnect["IBMQuantumPlatform", "New", Authentication -> {"apikey" -> api}]
```

<!-- => ServiceObject["IBMQuantumPlatform", …] -->

A connected service object lists the QPUs your account can reach through `"Backends"`. Submitting and retrieving jobs is handled for you by [IBMJobSubmit]() and the [IBMJob]() handle, shown below. List the available backends:

```wl
#| eval: false
ibm["Backends"]
```

<!-- => {ibm_fez, ibm_marrakesh, ibm_kingston} -->

## Building and checking a circuit

Build the circuit in the Wolfram Language first, so you have an exact reference to compare the hardware against. This example prepares a three-qubit GHZ state $|\mathrm{GHZ}\rangle = \tfrac{1}{\sqrt{2}}(|000\rangle + |111\rangle)$, applies a quantum Fourier transform and measures all three qubits. The transform turns the GHZ entanglement into a definite interference pattern across the eight outcomes, and that structure is exactly what device noise erodes, which makes the circuit a sensitive benchmark:

```wl
qc = QuantumCircuitOperator[{"GHZ", "Fourier"[3], {1, 2, 3}}];
qc["Diagram"]
```

Because the circuit is small, plot its exact output distribution: the noiseless ground truth a QPU run approximates, and the reference we overlay with the hardware counts:

```wl
qc[]["ProbabilityPlot"]
```

Before spending QPU time, you can sample the circuit on a local noiseless simulator. [QuantumCircuitOperator]()'s `"Qiskit"` interface runs through a local Qiskit Aer simulator when no provider is given, returning a [QuantumMeasurement]() whose statistics converge to the exact distribution as the shot count grows (this needs the Python Qiskit stack):

```wl
#| eval: false
qc["Qiskit"]["Shots" -> 4096]["ProbabilityPlot"]
```

## Inspecting the circuit as OpenQASM

A QPU does not read a Wolfram Language expression: the circuit is sent as OpenQASM. <code>[QuantumQASM]()</code> is the single hub for it. Called on a circuit with no extra arguments it emits OpenQASM **3.0** natively, in pure Wolfram Language with no external dependency, serializing the circuit exactly as built. Measurements become classical-bit assignments:

```wl
QuantumQASM[qc]
```

<!-- =>
OPENQASM 3.0;
qubit[3] q;
bit[3] c;
U(1.5707963267948966, 0., 3.141592653589793) q[0];
ctrl(1) @ negctrl(0) @ U(3.141592653589793, 0., 3.141592653589793) q[0] q[1];
ctrl(1) @ negctrl(0) @ U(3.141592653589793, 0., 3.141592653589793) q[1] q[2];
U(1.5707963267948966, 0., 3.141592653589793) q[0];
ctrl(1) @ negctrl(0) @ U(0., 0., 1.5707963267948966) q[1] q[0];
ctrl(1) @ negctrl(0) @ U(0., 0., 0.7853981633974483) q[2] q[0];
U(1.5707963267948966, 0., 3.141592653589793) q[1];
ctrl(1) @ negctrl(0) @ U(0., 0., 1.5707963267948966) q[2] q[1];
U(1.5707963267948966, 0., 3.141592653589793) q[2];
swap q[0] q[2];
c[0] = measure q[0];
c[1] = measure q[1];
c[2] = measure q[2];
-->

The hub is also an importer: handing OpenQASM source (or a `.qasm` file) back to [QuantumCircuitOperator]() reconstructs the circuit, again with no external dependency. The round trip recovers a circuit equal to the original:

```wl
QuantumCircuitOperator[QuantumQASM[qc]] == qc
```

<!-- => True -->

OpenQASM 3.0 is the native and recommended target. For older toolchains that still require the legacy `qelib1.inc` gate set, pass `"Version" -> 2` to emit OpenQASM 2.0. This path renders through Qiskit, so it needs the Python Qiskit stack configured:

```wl
#| eval: false
QuantumQASM[qc, "Version" -> 2]
```

<!-- =>
OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg c[3];
h q[0];
cx q[0],q[1];
…
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
-->

## Running on hardware

With an active connection, [IBMJobSubmit]() is the one call that runs a circuit on a QPU. It transpiles the circuit against the backend's own error-aware `Target` (per-instruction error and duration), submits it through the `SamplerV2` primitive, and returns immediately with an asynchronous [IBMJob]() handle while the job sits in the queue. There is no manual OpenQASM, no raw request and no hex decoding of the classical register: the handle carries the per-bit to qubit map captured at transpile time, so results come back already in the Wolfram Language qubit order.

```wl
#| eval: false
job = IBMJobSubmit[qc, "ibm_fez"]
```

<!-- => IBMJob[<| Status: Queued, Backend: ibm_fez, Job ID: d8kcdur2d42s73ca1r40 |>] -->

The handle reports the job's status. Poll it while the job is queued, or pass `"Wait" -> True` to [IBMJobSubmit]() to block until the job reaches a terminal status:

```wl
#| eval: false
job["Status"]
```

<!-- => "Queued" -->

Refresh the handle to pull a fresh snapshot from the service. Once the status is `"Completed"`, applying the handle returns the measurement reconstructed from the hardware counts:

```wl
#| eval: false
qpu = job["Refresh"][]
```

<!-- => QuantumMeasurement[3 qubits, 4096 shots] -->

The same handle keeps every response IBM returns under `"Raw"` and exposes curated accessors on top, so the run's metadata is one query away. For example, the billed quantum runtime:

```wl
#| eval: false
job["QuantumSeconds"]
```

<!-- => Quantity[…, "Seconds"] -->

## Comparing exact and hardware results

`qpu` is the hardware measurement returned by the handle; the circuit `qc` gives the exact, noiseless reference. Take the corresponding Wolfram Language measurement:

```wl
#| eval: false
wl = qc[];
```

Putting the two distributions side by side shows how much of the interference pattern survives the hardware. Each [QuantumMeasurement]() exposes its outcome distribution through `"Probabilities"`; sorting both by outcome and overlaying them, the dominant `000`, `001` and `111` peaks and the all-but-forbidden `100` still come through, while decoherence erodes the contrast and pulls the rest toward a uniform distribution:

```wl
#| eval: false
BarChart[Transpose[Values @* KeySort /@ {wl["Probabilities"], qpu["Probabilities"]}],
    AspectRatio -> 1/2, Frame -> True, ChartLegends -> {"Exact WL", "QPU ibm_fez"},
    PlotLabels -> Keys[wl["Probabilities"]]]
```

<!-- => (bar chart: exact-WL vs QPU probability per basis state) -->

The differences are not in the queries but in device noise: each shot passes through a deep transpiled circuit on superconducting hardware, and the noise varies from job to job.
