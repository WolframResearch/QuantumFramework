---
Template: TechNote
Name: SendingQueriesToIBMQPUs
Title: Sending Queries to IBM QPUs
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/SendingQueriesToIBMQPUs
Keywords: [IBM Quantum, QPU, OpenQASM, Qiskit, ServiceConnect, sampler, transpile]
RelatedGuides: [QuantumComputation]
RelatedTutorials: [QPUServiceConnect]
---

This tutorial shows the recommended pipeline for running a Wolfram Language quantum circuit on a real IBM Quantum processing unit (QPU). You build and check a circuit symbolically, inspect it as OpenQASM with <code>[QuantumQASM]()</code>, run it on hardware, and compare the hardware counts against the exact result. IBM's cloud API changed substantially in 2025 (a new IBM Quantum Platform on IBM Cloud, with `qiskit-ibm-runtime` and the `SamplerV2` primitive); the steps below target that platform. For help, write to quantum@wolfram.com.

## Connecting to IBM Quantum

Create an account at [quantum.cloud.ibm.com](https://quantum.cloud.ibm.com/). From the dashboard you need two things: an **IAM API key** (the long secret string shown once when you create the key, not the key's display name) and the **instance CRN** (a `crn:v1:…` string identifying your instance).

The API key is a secret and is stored in your encrypted system credential store. The CRN is read from a local symbol that the connection consults when it builds the request headers, so set both before connecting:

```wl
#| eval: false
api = "XXX";                          (* your IAM API key secret *)
LocalSymbol["ibm_crn"] = "XXX";       (* your instance CRN, as a string *)
```

Loading the Wolfram Quantum Framework registers the IBM Quantum Platform service connection automatically, so <code>[ServiceConnect]()</code> recognizes the `"IBMQuantumPlatform"` name directly: no separate package load and no <code>[Quiet]()</code> wrapper are needed.

```wl
#| eval: false
ibm = ServiceConnect["IBMQuantumPlatform", "New", Authentication -> {"apikey" -> api}]
<!-- => ServiceObject["IBMQuantumPlatform", …] -->
```

A connected service object exposes three requests: list the available QPUs with `"Backends"`, submit a job with `"JobRun"`, and fetch results with `"JobResults"`:

```wl
#| eval: false
ibm["Backends"]
<!-- => {ibm_fez, ibm_marrakesh, ibm_kingston} -->
```

The exact set of backends depends on which instances your account can reach.

## Building and checking a circuit

Build the circuit in the Wolfram Language first, so you have an exact reference to compare the hardware against. This example prepares a three-qubit GHZ state, applies a quantum Fourier transform, and measures all three qubits:

```wl
qc = QuantumCircuitOperator[{"GHZ", "Fourier"[3], {1, 2, 3}}];
qc["Diagram"]
```

Because the circuit is small, compute its exact output distribution symbolically. This is the noiseless ground truth that a QPU run approximates:

```wl
FullSimplify /@ qc[]["Probabilities"] // Dataset
<!-- => <|{0,0,0}->1/4, {0,0,1}->(2+Sqrt[2])/16, {0,1,0}->1/8, {0,1,1}->(2-Sqrt[2])/16,
          {1,0,0}->0, {1,0,1}->(2-Sqrt[2])/16, {1,1,0}->1/8, {1,1,1}->(2+Sqrt[2])/16|> -->
```

Before spending QPU time, sample the circuit on a local noiseless simulator. <code>[QuantumCircuitOperator]()</code>'s `"Qiskit"` interface runs through a local Qiskit Aer simulator when no provider is given, returning a <code>[QuantumMeasurement]()</code> whose statistics converge to the exact distribution as the shot count grows:

```wl
#| eval: false
qc["Qiskit"][Shots -> 4096]["ProbabilityPlot"]
```

## Inspecting the circuit as OpenQASM

<code>[QuantumQASM]()</code> is the single hub for OpenQASM. Called on a circuit with no extra arguments it emits OpenQASM **3.0** natively, in pure Wolfram Language with no Qiskit dependency, serializing the circuit exactly as built. Measurements become classical-bit assignments:

```wl
QuantumQASM[qc]
<!-- =>
OPENQASM 3.0;
qubit[3] q;
bit[3] c;
U(1.5708, 0., 3.14159) q[0];
ctrl(1) @ negctrl(0) @ U(3.14159, 0., 3.14159) q[0] q[1];
…
c[0] = measure q[0];
c[1] = measure q[1];
c[2] = measure q[2];
-->
```

Pass `"Version" -> 2` to get OpenQASM 2.0 instead (this path uses Qiskit to render the legacy `qelib1.inc` gate set):

```wl
#| eval: false
QuantumQASM[qc, "Version" -> 2]
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
```

The hub is also an importer: handing OpenQASM source (or a `.qasm` file) back to <code>[QuantumCircuitOperator]()</code> reconstructs the circuit, with no Qiskit needed for the import:

```wl
QuantumCircuitOperator[QuantumQASM[qc]] // Head
<!-- => QuantumCircuitOperator -->
```

## Running on hardware

The simplest and most robust way to run on a QPU is to set the circuit's evaluation `Method` to Qiskit with an IBM provider and a backend. This transpiles the circuit to the backend's native gate set and connectivity, runs it through the `SamplerV2` primitive, and returns a <code>[QuantumMeasurement]()</code> with the qubit ordering already aligned to the Wolfram Language convention, so no manual bit-string bookkeeping is needed:

```wl
#| eval: false
result = qc[Method -> {"Qiskit", "Provider" -> "IBMProvider", "Backend" -> "ibm_fez"}]
<!-- => QuantumMeasurement[…] -->
```

This call blocks until the job finishes. To avoid tying up the kernel while the job sits in the queue, submit it asynchronously in a separate kernel with <code>[LocalSubmit]()</code> and store the result through a handler:

```wl
#| eval: false
LocalSubmit[
    Needs["Wolfram`QuantumFramework`"];
    qc = QuantumCircuitOperator[{"GHZ", "Fourier"[3], {1, 2, 3}}];
    qc[Method -> {"Qiskit", "Provider" -> "IBMProvider", "Backend" -> "ibm_fez"}],
    HandlerFunctions -> <|"TaskFinished" :> Function[result = #EvaluationResult]|>
]
<!-- => TaskObject[…] -->
```

The probabilities are computed from the measured frequencies (here over 4096 shots):

```wl
#| eval: false
result["Probabilities"]
<!-- => <|{0,0,0}->0.255, {0,0,1}->0.197, {0,1,0}->0.149, {0,1,1}->0.054,
          {1,0,0}->0.019, {1,0,1}->0.035, {1,1,0}->0.120, {1,1,1}->0.171|> -->
```

## Submitting OpenQASM directly

For full control you can transpile and submit the OpenQASM yourself. IBM's sampler requires an ISA circuit (mapped to the backend's native gate set and qubit connectivity). Giving <code>[QuantumQASM]()</code> a `"Provider"` and a `"Backend"` transpiles the circuit against that backend and emits the ISA OpenQASM in a single call, entirely on the IBM provider:

```wl
#| eval: false
qasm = QuantumQASM[qc, "Provider" -> "IBMProvider", "Backend" -> "ibm_fez", "Version" -> 2];
<!-- =>
OPENQASM 2.0;
include "qelib1.inc";
qreg q[156];
creg c[3];
rz(pi/2) q[142];
… ISA gates rz / sx / cz, laid out on the ibm_fez 156-qubit device …
measure q[…] -> c[…];
-->
```

Submit the OpenQASM string as a `sampler` job and capture the returned job id:

```wl
#| eval: false
job = ibm["JobRun", {"QASM" -> qasm, "ProgramID" -> "sampler", "Backend" -> "ibm_fez"}]
<!-- => <|id -> d8jhf29e8nrc73bjk610, backend -> ibm_fez|> -->
```

Poll the job until it reports `Completed`, then retrieve the raw results:

```wl
#| eval: false
ibm["RawJobDetails", "JobID" -> job["id"]]["state"]
<!-- => <|status -> Completed|> -->
```

```wl
#| eval: false
raw = ibm["RawJobResults", "JobID" -> job["id"]];
samples = raw[["results", 1, "data", "c", "samples"]];
```

The samples are hexadecimal readouts of the classical register. Going through the high-level `Method` path above avoids this manual decoding; if you do decode the raw register yourself, mind that IBM orders classical bits least-significant first, the reverse of the Wolfram Language convention.

## Comparing exact and hardware results

Putting the exact distribution next to the hardware run shows the structure surviving the device noise: the dominant `000`, `001`, and `111` outcomes and the near-zero `100` all come through, with the noise flattening the distribution. The hardware values below are the verified `ibm_fez` run from the previous section:

```wl
exact = N[Values[KeySort[FullSimplify /@ qc[]["Probabilities"]]]];
hardware = Values @ KeySort @ <|
    {0, 0, 0} -> 0.255, {0, 0, 1} -> 0.197, {0, 1, 0} -> 0.149, {0, 1, 1} -> 0.054,
    {1, 0, 0} -> 0.019, {1, 0, 1} -> 0.035, {1, 1, 0} -> 0.120, {1, 1, 1} -> 0.171|>;
BarChart[Transpose[{exact, hardware}],
    ChartLegends -> {"Exact WL", "QPU ibm_fez"},
    ChartLabels -> {Tuples[{0, 1}, 3], None},
    Frame -> True, AspectRatio -> 1/2]
```

The differences are not in the queries but in device noise: each shot passes through a deep transpiled circuit on superconducting hardware, and the noise varies from job to job.
