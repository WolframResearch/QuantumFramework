# Amazon Braket Support in QuantumFramework: Design Plan

**Status:** Design only. No kernel edits, no implementation. This document is the plan an implementer follows.

| Field | Value |
|---|---|
| QuantumFramework | `Wolfram/QuantumFramework 2.0.0`, repo HEAD `c3f21a43` (Qiskit/QASM/IBM kernel files unchanged since audit anchor `a601a187`) |
| AWSClient paclet | `2.1.1` at `/Applications/Wolfram 15.0.app/.../Links/AWSClient`; Braket API def `Resources/APIDefinitions/braket/2019-09-01/service.json` (`apiVersion 2019-09-01`, `rest-json`, SigV4) |
| Template to mirror | `IBMJobSubmit` / `IBMJob` in `QuantumFramework/Kernel/QuantumCircuitOperator/IBMQuantum.m` |
| Sibling investigation | "QF native QPU QASM" report (forthcoming) at `OngoingProjects/AWS Braket/QF-Native-QPU-QASM-Investigation.md`; owns the compilation-engine decision |
| Date | 2026-06-18 |

---

## 1. Big picture

A quantum circuit built in the Wolfram Language is a symbolic object. To run it on Amazon Braket you have to do three physical things: translate the circuit into a text the device's compiler accepts (OpenQASM), hand that text to a named device through a cloud API, and read back the measurement statistics from where the device deposited them. Braket is a marketplace, not a single machine: one API fronts several vendors (IonQ, Rigetti, IQM, QuEra) plus on-demand simulators, and each vendor speaks a slightly different gate dialect and returns results in a slightly different shape.

The goal of this work is to give a physicist the same clean experience for Braket that already exists for IBM. On the IBM side you write `job = IBMJobSubmit[circuit, backend]`, get back a handle immediately, and ask it physical questions later: `job["Status"]`, `job["Measurement"]`, `job["Probabilities"]`, `job["Counts"]`, `job["Cancel"]`. The handle is *lossless*: it keeps every raw service response and exposes curated accessors on top, so nothing the cloud returned is thrown away. We want `BraketJobSubmit` / `BraketJob` to be that, for Braket.

The constraint that shapes the whole design: this must **not be a thin wrapper** around today's blog snippet. The current documented Braket pipeline is dead (every line of it fails against AWSClient 2.1.1, see Section 3), and it is dead precisely because it hardcoded the gate set, the calling convention, and the result shape of one moment in time. A wrapper that re-hardcodes those will rot the same way the moment a new device appears, a vendor changes its dialect, AWSClient changes how it takes parameters, or Braket revises its result schema. So the design is organized around **discovering** what a device can do (rather than assuming it) and around **single insulation points** for each thing that drifts.

The physics payoff is concrete and identical to IBM's: the same `QuantumCircuitOperator` you simulate exactly in WL can be sent to real hardware, and the device's finite-shot, noisy outcome distribution can be overlaid on the exact distribution QF computes, so the gap *is* the device noise. Braket makes this a cross-technology comparison: trapped-ion (IonQ) versus superconducting (Rigetti, IQM) versus an exact simulator (SV1), all from one circuit object.

---

## 2. Design principles (the "not a thin wrapper" contract)

1. **Discover, never hardcode.** Nothing about a device's gate set, modifiers, qubit connectivity, shot limits, native pragmas, cost, or availability is written into the kernel. It all comes from Braket `GetDevice` at runtime (Section 7B). A device that ships next year works without a code change.
2. **One seam per drift axis.** Each thing that can change independently gets exactly one place in the code that knows about it:
   - AWSClient calling convention -> the transport shim (7A).
   - Braket device-capability schema -> the capability parser (7B).
   - Vendor gate dialects / verbatim rules -> the compilation layer (7C).
   - Braket result schema (SV1 vs QPU) and bit order -> the result decoder (7E).
   Everything else calls *through* these seams and stays ignorant of the quirk.
3. **Mirror IBM's object model.** Same accessor vocabulary (`Status`/`Refresh`/`Measurement`/`Probabilities`/`Counts`/`Cancel`/`Raw`/...), same Properties-vs-Actions split (local reads are bulk-mappable and side-effect-free; anything needing a live call is an Action). A physicist who learned `IBMJob` already knows `BraketJob`.
4. **Async-first.** Braket tasks are genuinely asynchronous (a QPU task can sit in a queue for minutes to hours, gated by the vendor's execution window). The primary surface returns a handle immediately, exactly as `IBMJobSubmit` does; blocking is opt-in.
5. **Prefer WL-native transport and decode.** Unlike the IBM path, which leans on qiskit for submission and result deserialization, Braket submission and S3 retrieval go through AWSClient (pure WL), and result decoding is pure WL. Python (qiskit + `qiskit_braket_provider`) appears, if at all, only as a transitional *compilation* fallback behind the compilation seam, never on the transport or decode path. This unifies credentials under one mechanism (`ServiceConnect["AWS"]`) instead of also requiring AWS keys in a boto3 environment.
6. **Lossless.** `BraketJob` keeps the verbatim `CreateQuantumTask` / `GetQuantumTask` responses and the raw `results.json` under `"Raw"`, with curated accessors on top, so no field Braket returns is lost.

**Non-goals.** No decoder / error-mitigation layer (the device returns raw counts). No Braket Hybrid Jobs (`CreateJob`) or analog/AHS programs in this phase: scope is gate-model `CreateQuantumTask`. No attempt to bill or reconcile against the AWS invoice; `"Cost"` is a pre-submission estimate (7D).

---

## 3. The problem: eight concrete breakages this design must erase

All eight were verified live on 2026-06-18 against AWSClient 2.1.1 and the device APIs. The current tutorial (`Documentation/English/Tutorials/QPUServiceConnect.md`, "Amazon Braket" section, lines 187-336) demonstrates the dead pattern: `braket["CreateQuantumTask", "DeviceArn" -> ..., "Action" -> ..., ...]` with bare trailing rules and bare-symbol action keys, `qc["Qiskit"]["QASM", "Provider" -> "AWSBraket"]`, and a hand-rolled `GetQuantumTask` -> S3 `GetObject` -> parse.

| # | Breakage (verified) | Theme | Seam that erases it |
|---|---|---|---|
| 1 | AWSClient 2.1.1 takes request params as **one list**: `braket["CreateQuantumTask", {"DeviceArn"->..., "Shots"->..., "RequestRegion"->...}]`. The blog's bare trailing rules error `Unknown option Shots`. | Transport | 7A transport shim |
| 2 | `ClientToken` is **required**; the idempotency auto-fill does not fire. Must pass `"ClientToken" -> CreateUUID[]`. | Transport | 7A |
| 3 | The `action` member carries the AWS `jsonvalue` trait: it must be a **string-keyed Association** (AWSClient JSON-encodes it). A pre-stringified JSON double-encodes and the server rejects it: `Provided action is not valid`. | Transport | 7A |
| 4 | Native `QuantumQASM[qc]` (which emits `U` and `ctrl @ negctrl @ U` modifier syntax) is **rejected by Braket's parser**. | Compilation | 7C compiler |
| 5 | `qc["Qiskit"]["QASM", "Provider"->"AWSBraket"]` needs region + creds in the **Python/boto3 env** and a separate `braket` Python environment (`qiskit_braket_provider`). | Compilation / creds | 7C + 8 |
| 6 | The AWSBraket backend handling in `qiskitInitBackend` (`Qiskit.m:436-441`) is **broken**: it hardcodes `get_backend('basic_simulator')` for *any* named backend, so you cannot transpile to a real QPU's gate set. | Compilation (kernel bug) | Stage 0 fix |
| 7 | Consequently QF cannot emit QPU-runnable QASM: IonQ rejects `cphaseshift`; native `gpi`/`gpi2` need a **verbatim box**; `gphase` is rejected and must be **stripped**. | Compilation | 7C |
| 8 | Result retrieval is **100% manual**: `GetQuantumTask` -> `outputS3Directory` -> S3 `GetObject` `results.json` -> parse, and shapes **differ by device** (SV1 returns per-shot `measurements`; IonQ returns aggregated `measurementProbabilities`), with bit-order bookkeeping. | Retrieval | 7E decoder |

Read top to bottom, the breakages cluster into the three layers the architecture is built around: a transport layer that must own the AWSClient quirks (1-3), a compilation layer that must produce device-runnable QASM from discovered capabilities (4-7), and a retrieval layer that must normalize heterogeneous result shapes (8).

---

## 4. Architecture: layered seams

```
  user  -->  BraketJobSubmit[qco, device, opts]            (D: public surface)
                 |
                 |  picks bucket/region, assembles action
                 v
  C  BraketCompile[qco, capabilityRecord, mode] -> qasm     (compilation seam)
        engine = native-WL  | qiskit fallback (behind seam)
                 |
                 v
  A  braketRequest["CreateQuantumTask", <|...|>]            (transport seam)
        owns: one-list form, ClientToken, jsonvalue action
                 |                              ^
                 v                              |
            AWS Braket  ----(S3 results.json)---+
                 |
                 v
  D  BraketJob[<|"TaskARN"->..., "Device"->rec, "Raw"->...|>]
        "Refresh" (Action) --> A GetQuantumTask + S3 GetObject
                            --> E braketDecodeResult (decode seam)
                                  owns: SV1 vs QPU shape, bit order
        accessors: Status/Measurement/Counts/Probabilities/Cost/QueuePosition/...

  B  BraketDevice / capability cache  <-- GetDevice (jsonvalue deviceCapabilities)
        owns: supportedOperations, supportedModifiers, supportedPragmas,
              connectivity, shotsRange, deviceCost, executionWindows, status
        consumed by C (gate set) and D (region, shots bounds, cost)
```

Five internal layers (A transport, B capability, C compilation, D submission/handle, E decode) plus a provider-abstraction seam above D. The dependency arrows only point one way: the handle and the decoder never reach into AWSClient directly, the compiler never reaches into the service object, and capability discovery is the single source of truth that the compiler and the handle both read. Each lettered subsection of Section 7 specifies one layer.

A useful contrast with IBM: there, qiskit owns the wire format end to end (circuit serialization, submission via `SamplerV2`, result deserialization), and QF mostly forwards. Here, **QF owns the wire format**: AWSClient is a transport, OpenQASM is the interchange, and the capability record is QF's own normalized device model. That is more code, but it is the code that makes Braket support robust instead of a wrapper, and it removes the second credential path.

---

## 5. Public surface: recommendation

**Primary (recommended): `BraketJobSubmit` / `BraketJob`.** Direct mirror of the IBM pair.

```
job = BraketJobSubmit[qco, device, opts]   (* device = ARN string | BraketDevice | Automatic *)
job["Status"]            (* "Queued" | "Running" | "Completed" | ... *)
job["Measurement"]       (* QuantumMeasurement once Completed *)
job["Probabilities"]     (* normalized over outcomes *)
job["Counts"]            (* exact integer counts, ascending-qubit order *)
job["QueuePosition"]
job["Cost"]              (* estimate, Quantity[..., "USDollars"] *)
job["Cancel"]            (* Action *)
job["Refresh"]           (* Action: re-query + fetch+decode S3 when done *)
job["Raw"]               (* lossless service + results.json payloads *)
```

**Secondary (optional, for symmetry): `qco[Method -> {"Braket", "Device" -> arn, "Shots" -> n, ...}]`.** A thin alias that calls `BraketJobSubmit` and returns the handle (or, with `"Wait" -> True`, blocks and returns the decoded `QuantumMeasurement`, the way the qiskit Method path returns one through `iMethodIBMJob`).

**Recommendation: ship `BraketJobSubmit`/`BraketJob` as the primary, documented surface, and add the `Method -> {"Braket", ...}` alias as a convenience.** Rationale and trade-offs:

- **Async is the honest model.** A QPU `CreateQuantumTask` returns an ARN in milliseconds; the result lands in S3 minutes to hours later, gated by the vendor's `executionWindows`. A `Method` path that *looks* like a local apply (`qco[Method -> ...]` returning a `QuantumMeasurement`) invites the user to block a kernel for an hour. `BraketJobSubmit` returning a handle makes the latency visible and lets the user poll, walk away, or reconnect to the ARN later. This is exactly why the IBM design made `IBMJobSubmit` the async front door and relegated blocking to the Method path.
- **The handle is the natural home** for everything Braket-specific that has no place in a bare `QuantumMeasurement`: queue position, S3 location, cost estimate, device record, failure reason, the raw payloads. Lossless retention requires an object; a `Method` that returns only a measurement discards all of it.
- **Symmetry is a real usability asset.** `BraketJob` and `IBMJob` sharing an accessor vocabulary means tutorials, downstream analysis code, and the user's muscle memory transfer directly. This also sets up the provider-abstraction seam (Section 7F).
- **Keep the `Method` alias** because (a) it matches how every other QF execution backend is reached (`Method -> "Schrodinger"`, `"TensorNetwork"`, `"Qiskit"`), so omitting it is a surprising gap, and (b) for SV1 a short blocking call (`Method -> {"Braket", "Device" -> "SV1", "Wait" -> True}`) is genuinely convenient and cheap. The alias must route to the same `BraketJobSubmit` engine, not duplicate it.
- **Do not** extend the existing python `qiskitApply` Method path (`Method -> "Qiskit"`, `Provider -> "AWSBraket"`) to carry Braket execution. That path is broken (breakage #6), python-bound, and credential-duplicating. The new `"Braket"` Method is a parallel, WL-native route; the old `Provider -> "AWSBraket"` stays only as a compilation fallback behind the 7C seam.

New public heads to register in `PacletInfo.wl` `Symbols` (with `::usage` written in the doc page Usage section, per the auto-generated-`Usage.m` convention): `BraketJobSubmit`, `BraketJob`, and `BraketDevice` (the capability handle, 7B). Whether `BraketDevice` is public or PackageScope is an open question (Section 11); it is useful enough as a discovery/inspection object to justify exporting.

---

## 6. The active connection

A single `braketActiveConnection[]`, mirroring `ibmActiveConnection[]` (`IBMQuantum.m:29`): it returns the active AWS `ServiceObject` or `$Failed`, and **never creates or prompts**. From it, two sub-services are obtained the way the AWS connector exposes them:

```
aws    = ServiceFramework`GetDefaultServiceObject["AWS"]   (* or the user's stored AWS connection *)
braket = ServiceExecute[aws, "GetService", {"Name" -> "Braket"}]
s3     = ServiceExecute[aws, "GetService", {"Name" -> "S3"}]
```

`braketActiveConnection[]` caches `{aws, braket, s3}` together. If there is no active AWS connection, every public entry point returns the same shaped failure `IBMJobSubmit` gives ("No active AWS connection. Run `ServiceConnect[\"AWS\"]` first."). The handle stores the AWS `ServiceObject` (under `"ServiceObject"`, as `IBMJob` does) so Actions can re-query.

---

## 7. Layer-by-layer design

### 7A. Transport shim: `braketRequest` (erases breakages 1, 2, 3)

**One private function is the only code in the kernel that calls a raw Braket service operation.** It centralizes every AWSClient 2.1.1 calling-convention quirk, so a future AWSClient that changes the convention again touches exactly this function.

```
braketRequest[op_String, params_Association, opts___] := Module[{svc, payload},
   svc = (* braket sub-service from braketActiveConnection[] *);
   payload = prepareBraketParams[op, params];     (* the quirk-handling *)
   svc[op, payload]                                (* AWSClient 2.1.1: ONE list *)
]
```

`prepareBraketParams` owns:

1. **One-list form (breakage 1).** Parameters are assembled into a single `{ "Key" -> val, ... }` list and passed as the *one* argument after the operation name. Never bare trailing rules. (Verified: bare rules error `Unknown option Shots`.)
2. **ClientToken auto-fill (breakage 2).** For the idempotent mutating operations (`CreateQuantumTask`, `CancelQuantumTask`), if the caller did not supply `"ClientToken"`, inject `CreateUUID[]`. A fresh UUID per logical submission; on a `Refresh`-driven retry the *same* token must be reused so AWS treats it idempotently rather than creating a duplicate task. (So the token is generated at `BraketJobSubmit` time and stored on the handle, not regenerated inside `braketRequest` on every call. `prepareBraketParams` only fills it when absent.)
3. **`jsonvalue` members as string-keyed Associations (breakage 3).** `action` and `deviceParameters` carry the AWS `jsonvalue` trait (confirmed in `service.json`: `CreateQuantumTaskRequest.action -> JsonValue [jsonvalue]`, `deviceParameters -> ... [jsonvalue]`). AWSClient JSON-encodes them itself. So they must be handed over as **Associations with string keys**, e.g.

   ```
   "action" -> <|
       "braketSchemaHeader" -> <|"name" -> "braket.ir.openqasm.program", "version" -> "1"|>,
       "source" -> qasmString
   |>
   ```

   never as a pre-built JSON string (that double-encodes -> `Provided action is not valid`). The dead blog used bare *symbols* (`braketSchemaHeader`, `name`, `source`) as keys, which is wrong twice over (undefined symbols, and not strings). `prepareBraketParams` is the one place that builds these Associations; callers pass the QASM string and the shim wraps it in the schema-headed envelope.
4. **`RequestRegion`.** Braket devices are region-pinned (IonQ in `us-east-1`, Rigetti in `us-west-1`, IQM in `eu-north-1`, ...). Every call carries `"RequestRegion" -> region`, where `region` is derived from the device (7B). This is centralized here so no caller forgets it.

**Why this is the drift seam:** AWSClient already changed its convention once (the blog was written for the pre-2.1.1 bare-rule form). When it changes again, the symptom will be a transport failure from exactly one function. Optionally, `prepareBraketParams` can branch on `PacletObject["AWSClient"]["Version"]` if two conventions must coexist, but the simpler contract is: code to the installed version, isolate it here.

### 7B. Device capability model: `BraketDevice` + capability cache (the robustness engine)

This layer is what makes "robust to new devices and vendors" true rather than aspirational. **Everything the rest of the system needs to know about a device is read from `GetDevice` at runtime; nothing is hardcoded.**

`GetDevice` (`GET /device/{deviceArn}`) returns (from `service.json`): `deviceArn`, `deviceName`, `providerName`, `deviceType` (`QPU` | `SIMULATOR`), `deviceStatus` (`ONLINE` | `OFFLINE` | `RETIRED`), `deviceQueueInfo`, and the prize: **`deviceCapabilities`, a `jsonvalue` field** that AWSClient decodes into a (deeply nested) string-keyed Association. The Braket device-capability schema inside it carries, per device:

- **`action` / supported program types** (e.g. `braket.ir.openqasm.program`) and, under the OpenQASM action, **`supportedOperations`** (the gate set the device accepts), **`supportedModifiers`** (ctrl/negctrl/pow/inv and their limits), **`supportedPragmas`** / `forbiddenPragmas` (including whether `verbatim` is available), and `supportedResultTypes`.
- **`paradigm`**: `qubitCount`, and **`connectivity`** (`fullyConnected` or an explicit coupling graph).
- **`deviceParameters`** schema (what the optional `deviceParameters` action member may carry).
- **`service`**: **`shotsRange`** (`{minShots, maxShots}`), **`deviceCost`** (`{price, unit}`, e.g. price per task and per shot), **`executionWindows`** (when the QPU runs), `deviceLocation` (region), and update timestamps.
- **`provider`**: vendor-specific native-gate detail (e.g. IonQ `gpi`, `gpi2`, `ms`; the native basis used in verbatim mode).

The design defines a **normalized capability record** `braketCapabilityRecord[deviceArn]` (an Association) that flattens the parts the kernel uses into stable QF keys, computed by a **schema-tolerant parser** (`Lookup` with `Missing` fallbacks down each path, never positional indexing into the vendor JSON), validated against the embedded `braketSchemaHeader`:

```
<|
  "DeviceArn" -> ..., "DeviceName" -> ..., "Provider" -> ..., "DeviceType" -> "QPU"|"SIMULATOR",
  "Status" -> "ONLINE"|"OFFLINE"|"RETIRED",
  "Region" -> "us-east-1",                          (* for RequestRegion + bucket default *)
  "SupportedGates" -> {"h","rx","ry","rz","cnot",...},   (* logical gate set *)
  "NativeGates"    -> {"gpi","gpi2","ms",...},            (* verbatim-mode basis *)
  "SupportedModifiers" -> {...}, "VerbatimAvailable" -> True|False,
  "Connectivity" -> "FullyConnected" | <coupling graph>,
  "QubitCount" -> 32,
  "ShotsRange" -> {1, 100000},
  "DeviceCost" -> <|"PerTask" -> Quantity[0.3,"USDollars"], "PerShot" -> Quantity[0.01,"USDollars"]|>,
  "ExecutionWindows" -> {...},
  "SchemaVersion" -> "...",                          (* the header that was parsed *)
  "Raw" -> <the full deviceCapabilities Association>  (* lossless *)
|>
```

`BraketDevice[arn]` is the public face: it runs `GetDevice`, builds the record, and is the object the user can inspect (`dev["SupportedGates"]`, `dev["ShotsRange"]`, `dev["Status"]`) and pass to `BraketJobSubmit`. A `SearchDevices` (`POST /devices`) helper lists what is available, mirroring the blog's discovery step but returning `BraketDevice` handles.

**Robustness properties this buys:**

- A **new device or vendor** is usable the day it appears on Braket: its gate set, connectivity, shot bounds, and cost are read, not assumed.
- The **compiler (7C)** asks the record for the gate set instead of hardcoding `cnot`/`cphaseshift`/`gpi`. When IonQ's accepted gates change, the record changes, the compiler follows.
- **Pre-flight validation** is free and friendly: before submitting, `BraketJobSubmit` can check `Status === "ONLINE"`, `shots` within `ShotsRange`, and the circuit width against `QubitCount`, failing in milliseconds with a clear message instead of after a slow rejected submission.
- The parser tolerates **schema drift**: an added field is ignored, a missing field degrades to `Missing[...]` and a documented fallback, and the parsed `braketSchemaHeader` is recorded so a future schema bump is visible rather than silently mis-read.

A short-lived in-session cache (keyed by ARN) avoids re-fetching `GetDevice` on every submission; `deviceStatus` and `executionWindows` are the only volatile parts and are re-read on submission.

### 7C. Compilation strategy (erases breakages 4, 5, 7; depends on the Stage 0 fix for 6)

This is the layer where the "logical-vs-native" and "qiskit-vs-WL" decisions live, behind one interface:

```
BraketCompile[qco_, capRec_, mode_, opts___] -> qasmString | Failure
   mode = "Logical" (default) | "Native"
```

Everything above this interface (submission, the handle) is engine-agnostic; everything below is swappable. That swappability is the point: it lets us ship a working path now and replace the engine later without touching the rest.

**Two compilation modes, chosen from the capability record, not hardcoded:**

- **Logical / recompiled mode (default).** Emit a portable OpenQASM 3 program in the device's *accepted logical* gate set (`capRec["SupportedGates"]`), and let Braket's own compiler lower it to native gates and route it onto the connectivity. This is the right default because it is robust (the device owns the hard part of native compilation and calibration-aware routing) and avoids the verbatim box. The circuit is decomposed to the discovered logical basis: e.g. `cphaseshift` (which IonQ rejects, breakage 7) is rewritten in terms of gates the record lists. `U` and `ctrl @ negctrl @ U` modifier syntax (which Braket's parser rejects, breakage 4) is never emitted: the lowering targets named gates only.
- **Native / verbatim mode (opt-in, or forced when the device does not recompile).** Emit the device's *native* gates (`capRec["NativeGates"]`, e.g. IonQ `gpi`/`gpi2`/`ms`) wrapped in a verbatim box,

  ```
  #pragma braket verbatim
  box { gpi(...) $0; gpi2(...) $1; ms(...) $0, $1; ... }
  ```

  which Braket requires for native-gate circuits (breakage 7), and only when `capRec["VerbatimAvailable"]` is true. This mode gives the user exact control of the pulses-as-gates the device runs, at the cost of doing the native decomposition in QF.

**In both modes, `gphase` / global phase is stripped** (Braket rejects `gphase`, breakage 7). A QF circuit can carry a `GlobalPhase` / `"GlobalPhase"[...]` element; the emitter drops it. This is physically safe for sampling: a global phase is unobservable in measurement statistics, so removing it cannot change any count or probability. (This must be a deliberate, documented step in the emitter, not an accident: note it in the code as "Braket rejects gphase; global phase is unobservable in sampling.")

**Engine choice (the consequential decision), with trade-offs:**

- **`qiskit_braket_provider` (the existing path).** *Pros:* already wired (`qiskitQASM`, `Qiskit.m:648-652`, calls `convert_qiskit_to_braket_circuit(...).to_ir(OPENQASM)`); it already knows the verbatim box, the native-gate conversion, and gphase handling. *Cons:* drags in a heavy Python dependency and a *separate* `braket` Python environment (`PythonTools.m:16`); needs AWS region + credentials in the **boto3** environment, duplicating the `ServiceConnect["AWS"]` credentials (breakage 5); and its backend selection is broken in QF (breakage 6) so it cannot currently target a real QPU's gate set at all. It also couples QF to two external release cadences (qiskit and the provider).
- **Native WL transpiler.** *Pros:* QF already has a native, dependency-free OpenQASM emitter (`qasmEmitCircuit`, `QuantumQASM.m:336`) and decomposition machinery (KAK / `Decompose`); a gate-set-directed WL lowering removes Python from the path entirely, unifies credentials under `ServiceConnect["AWS"]`, and is immune to qiskit/provider drift. This is the only engine that fully honors design principle 5. *Cons:* it must be built: gate-set-directed decomposition to a discovered basis, the verbatim-box emitter, and gphase stripping are real work. That work is exactly what the sibling "QF native QPU QASM" investigation is scoping.

**Recommendation.** Make the **native WL transpiler the strategic target** and the **default once it lands**, with the qiskit path as a clearly-isolated transitional fallback reachable only behind the `BraketCompile` seam. Concretely:

- Stage 0 fixes breakage 6 so the qiskit fallback at least *works* for QPUs while the native path is built.
- `BraketCompile` is written interface-first, with a pluggable `"Engine"` option (`Automatic` | `"WL"` | `"Qiskit"`). `Automatic` prefers the native WL engine when it can express the target gate set, and falls back to qiskit otherwise. Submission and the handle never know which engine ran.
- The native engine's design is owned by the sibling investigation; this plan commits only to the *seam* and the *mode semantics* (logical vs native, gphase stripping, capability-driven gate set), so the two efforts compose. When the sibling report lands, its emitter becomes the `"WL"` engine behind this interface.

This keeps the consequential engine decision *reversible*: shipping on qiskit first does not lock QF in, because nothing above `BraketCompile` depends on it.

### 7D. Submission and the `BraketJob` handle (mirrors `IBMJob`)

**`BraketJobSubmit[qco, device : _String | _BraketDevice | Automatic, opts]`** (`Options`: `"Shots" -> 1000`, `"Mode" -> "Logical"`, `"Engine" -> Automatic`, `"OutputS3Bucket" -> Automatic`, `"OutputS3KeyPrefix" -> "tasks"`, `"DeviceParameters" -> <||>`, `"Wait" -> False`, `"Tags" -> <||>`) does, in order:

1. `braketActiveConnection[]`; fail clearly if none.
2. Reject non-qubit circuits up front with the same message `QuantumQASM` gives (`qasmQubitsQ` / `qasmNonQubitFailure`, `QuantumQASM.m:321`), exactly as `IBMJobSubmit` does (`IBMQuantum.m:155`). OpenQASM models qubits only.
3. Resolve the device to a capability record (7B). `Automatic` picks a sensible default (SV1 simulator). Pre-flight checks: `Status === "ONLINE"`, `MinShots <= shots <= MaxShots`, circuit width `<= QubitCount`.
4. `qasm = BraketCompile[qco, capRec, mode, "Engine" -> ...]` (7C). Surface its `Failure` verbatim.
5. Resolve bucket + region (Section 8). Generate the `ClientToken` here (one per submission) and keep it for idempotent retries/cancel.
6. `braketRequest["CreateQuantumTask", <| "DeviceArn" -> arn, "Shots" -> shots, "Action" -> <schema-headed action>, "OutputS3Bucket" -> bucket, "OutputS3KeyPrefix" -> prefix, "ClientToken" -> token, "DeviceParameters" -> devParams (* if any *), "RequestRegion" -> region |>]` (7A wraps `action`/`deviceParameters` and the one-list form).
7. Read `quantumTaskArn` from the response; build the handle:

   ```
   BraketJob[<|
      "TaskARN" -> arn, "Device" -> capRec, "ServiceObject" -> aws,
      "Shots" -> shots, "OutputS3Bucket" -> bucket, "OutputS3Directory" -> Missing[],
      "ClientToken" -> token, "Region" -> region, "Action" -> action
   |>]
   ```
8. `If[waitQ, braketWaitFor[job], job["Refresh"]]` (one initial status fetch so the summary box has a status), exactly as `IBMJobSubmit` ends (`IBMQuantum.m:183`).

**`BraketJob` accessors**, mapped to the IBM vocabulary so they read identically:

| `BraketJob` | Meaning | Source |
|---|---|---|
| `"TaskARN"` | the Braket task id (Braket's analog of `IBMJob["ID"]`) | local |
| `"Device"` / `"DeviceName"` / `"Provider"` | capability record fields | local |
| `"Status"` | normalized status (see below) | local if cached, else `GetQuantumTask` |
| `"Shots"` | requested shots | local |
| `"NumSuccessfulShots"` | shots actually returned | `GetQuantumTask.numSuccessfulShots` |
| `"QueuePosition"` | position in the device queue | `GetQuantumTask.queueInfo.position` |
| `"Cost"` | **estimate** `PerTask + PerShot * shots` from `deviceCost` | capability record |
| `"CreatedAt"` / `"EndedAt"` / `"Duration"` | timestamps | `GetQuantumTask.createdAt/endedAt` |
| `"FailureReason"` | why a `FAILED` task failed | `GetQuantumTask.failureReason` |
| `"Counts"` | exact integer counts, ascending-qubit order | `results.json` via 7E |
| `"Measurement"` | `QuantumMeasurement` | 7E |
| `"Probabilities"` | normalized over outcomes | 7E |
| `"OutputS3Bucket"` / `"OutputS3Directory"` | where results live | local / `GetQuantumTask` |
| `"Raw"` | lossless: task response + `results.json` | stored |
| `"Refresh"` (**Action**) | re-query `GetQuantumTask`; on `COMPLETED`, fetch + decode S3 | live |
| `"Cancel"` (**Action**) | `CancelQuantumTask` (needs the stored `ClientToken`) | live |
| `"Properties"` / `"Actions"` | the two accessor sets | local |

**Status normalization.** Braket's enum is `{CREATED, QUEUED, RUNNING, COMPLETED, FAILED, CANCELLING, CANCELLED}` (`service.json`). Map to the IBM-aligned display vocabulary so the summary box and any cross-provider code are uniform: `CREATED`/`QUEUED` -> `"Queued"`, `RUNNING` -> `"Running"`, `COMPLETED` -> `"Completed"`, `FAILED` -> `"Failed"`, `CANCELLING`/`CANCELLED` -> `"Cancelled"`. Terminal set (stop polling) = `{COMPLETED, FAILED, CANCELLED}`, the analog of `$ibmTerminal` (`IBMQuantum.m:263`). Reuse the `$ibmStatusColor` scheme for the summary box.

**Properties-vs-Actions split** (the IBM contract, `IBMQuantum.m:518/532`): `"Properties"` are local, network-free, side-effect-free, and bulk-mappable (`AssociationMap[job, job["Properties"]]` never touches the network or duplicates a task). `"Actions"` (`"Refresh"`, `"Cancel"`) each need a live call and are kept out of `"Properties"` so a "dump everything" never re-queries, cancels, or blocks. An unknown-property catch-all returns a `Failure` listing both sets, as `IBMJob` does (`IBMQuantum.m:541`).

**Cost is an estimate, and the accessor must say so.** Unlike `IBMJob["QuantumSeconds"]`, which IBM reports back, Braket's `GetQuantumTask` returns no dollar figure. `"Cost"` is computed from the device's published `deviceCost` times shots and returned as `Quantity[..., "USDollars"]`, with the understanding (documented, and ideally flagged in the value's provenance) that the authoritative number is on the AWS bill. This is honest and still useful for "should I run this on a QPU?" budgeting.

**Lossless `"Raw"`** holds the verbatim `CreateQuantumTask` response, the latest `GetQuantumTask` response, and the parsed `results.json`, so any field Braket returns that QF does not yet surface is still reachable, the way `IBMJob["Raw"]` works (`IBMQuantum.m:512`).

### 7E. Result decoder: `braketDecodeResult` (erases breakage 8)

**One function normalizes every device's result shape and bit order.** `"Refresh"` on a `COMPLETED` task does: `outputS3Directory = GetQuantumTask[...].outputS3Directory`, then S3 `GetObject` `results.json` from `OutputS3Bucket`/that directory, then `braketDecodeResult[json, capRec]` -> `Counts`. All of that is encapsulated in the handle's Action; the user never sees S3.

`braketDecodeResult` handles the two known shapes (and is written so a third degrades gracefully, not crashes):

- **Simulators (SV1, DM1, TN1): per-shot `measurements`.** A list of per-shot bit lists (plus `measuredQubits`). Decode directly into `Counts` over the bit-list outcomes, one entry per shot.
- **QPUs (IonQ and others): aggregated `measurementProbabilities`.** A string-keyed map `bitstring -> probability`. Reconstruct `Counts` by `Round[probability * shots]` per outcome (when `measurementCounts` is present, use it directly; some payloads carry both). Keep the probabilities too.

**Bit order is the subtle part and must be pinned by a live known-state test, not assumed.** Braket bitstrings and IBM/qiskit bitstrings use opposite endianness conventions, and the decoder must reorder Braket outcomes into QF's ascending-qubit convention (qubit 0 first), the same convention `IBMJob` decodes to via `qiskitReorderByQubit` (`Qiskit.m:484`) and the `measuredQubits` map. The decoder reads Braket's `measuredQubits` to map result columns to circuit qubits. The single verification that locks the convention down: submit a circuit preparing a known asymmetric basis state (e.g. an `X` on qubit 0 only, expecting outcome `|100...>` in QF's ordering) to SV1 and confirm the decoded `Counts` key matches the QF-simulated `QuantumMeasurement` of the same circuit. Until that test passes, treat the bit order as unverified.

The decoder also stashes Braket's `taskMetadata` / `additionalMetadata` (which carries the compiled program and device-side info) into `"Raw"` for losslessness.

### 7F. Provider-abstraction seam (future vendors)

`IBMJob` and `BraketJob` already converge on one shape: submit -> handle; handle answers `Status`/`Refresh`/`Measurement`/`Counts`/`Probabilities`/`Cancel`/`Raw` with a Properties/Actions split. The seam is to **name that shape** as an internal protocol so a third backend (a future direct-vendor connector, an on-prem device, another marketplace) plugs in by implementing it:

- a `submit[circuit, deviceSpec, shots, opts] -> handle` for the provider, and
- a handle answering the shared accessor set.

The `qco[Method -> {provider, ...}]` dispatch then routes on the provider name to the right `submit`. This plan does *not* build a heavy abstraction layer now (premature for two providers); it (a) keeps `BraketJob`'s accessor names and Properties/Actions split identical to `IBMJob`'s so the protocol is implicit and real, and (b) routes the new `"Braket"` Method through the same dispatch point the `"Qiskit"`/`"IBM"` methods use, so formalizing the protocol later is a refactor, not a redesign. That is the seam the explicit "non-Braket backends could plug in" requirement asks for, sized to two providers rather than over-engineered for ten.

---

## 8. Credentials, region, and the S3 bucket

- **Credentials: one mechanism.** The active `ServiceConnect["AWS"]` connection (SigV4 under the hood) is the single credential path for *both* submission and S3 retrieval. The Braket and S3 sub-services come from the same connection (Section 6). This is the design's credential win over the qiskit path, which additionally needs AWS keys in a boto3 env (breakage 5); the WL-native transport + decode path removes that second path entirely. (The qiskit compilation fallback, if used, still needs boto3 creds *only because the provider insists on constructing an `AWSBraketProvider`*; the native WL engine removes even that.)
- **Region: derived from the device, set per call.** A device is pinned to a region (`capRec["Region"]`, read from `deviceCapabilities.service.deviceLocation` or the ARN). `braketRequest` sets `"RequestRegion"` to it on every call (7A). The output bucket must be reachable from that region. Centralizing region in the transport shim and the capability record means no caller hardcodes `us-east-1` (the blog's bug-in-waiting).
- **S3 bucket: explicit, then canonical default, then a clear failure.** Resolution order in `BraketJobSubmit`: (1) the user's `"OutputS3Bucket"` option, if given; (2) the canonical Braket default bucket `amazon-braket-{accountId}-{region}` (account id from `STS GetCallerIdentity` via the AWS connection, or from the device ARN's account; region from the capability record); (3) if neither resolves to an existing bucket, return a clear `Failure` telling the user to create one (link to the S3 console), *not* a raw AWS error. **Do not auto-create the bucket silently:** creating an S3 bucket is an outward-facing, billable side effect, so it is offered (an explicit `"CreateBucket" -> True` opt-in) but never default. The default-bucket *name* is constructed and probed (a `HeadBucket`/list check), which is read-only.

---

## 9. Staged implementation plan

Ordered so each stage is independently testable and the riskiest external dependency (the broken qiskit path) is decoupled early.

**Stage 0: Kernel fix, backend selection (erases breakage 6).** In `qiskitInitBackend` (`Qiskit.m:436-441`), the `$AWSProvider` branch hardcodes `provider.get_backend('basic_simulator')` for any named backend. Change it to `provider.get_backend(backend_name)` when a name is given (mirroring the IBM branch at `Qiskit.m:423-426`), so the qiskit fallback can target a real QPU's gate set. Smallest, highest-leverage change; unblocks the transitional compilation engine. Guard: this stage is the *only* kernel edit before the new files; it is a one-line-shaped fix with a regression test (transpile a circuit against a named Braket device and assert the basis matches the device, run offline against `GenericV2` where possible). Hand this stage to the adversarial WL verification loop (`/wl-verify`) since it touches a `.m` file.

**Stage 1: Transport shim + capability model (Layers A, B).** New file `Kernel/QuantumCircuitOperator/Braket.m` (or a `Kernel/Braket/` subdir if it grows): `braketActiveConnection`, `braketRequest` + `prepareBraketParams` (one-list, ClientToken, jsonvalue), `BraketDevice`, `braketCapabilityRecord`, `SearchDevices` helper. Fully testable offline against recorded `GetDevice` JSON fixtures (SV1, IonQ, Rigetti) and a mocked service object that records the payload it was handed (so the one-list form, the presence of `ClientToken`, and the string-keyed `action` Association are asserted without a network call).

**Stage 2: Compilation layer (Layer C).** `BraketCompile[qco, capRec, mode, "Engine" -> ...]` with the logical/native modes, gphase stripping, and the engine plug. Initial `"Engine" -> "Qiskit"` wrapping the (now-fixed) `qiskitQASM` path; design the interface for the `"WL"` engine the sibling investigation will deliver. Coordinate the gate-set-directed lowering and verbatim-box emission with `QF-Native-QPU-QASM-Investigation.md` so the `"WL"` engine drops in behind the seam.

**Stage 3: Submission + handle + decoder (Layers D, E).** `BraketJobSubmit`, `BraketJob` with the full accessor set and Properties/Actions split, `braketWaitFor`, `braketDecodeResult`. Bucket/region resolution (Section 8). Decoder tested offline against recorded `results.json` fixtures for SV1 (per-shot `measurements`) and IonQ (`measurementProbabilities`), including the bit-order known-state fixture.

**Stage 4: `Method` alias + provider seam.** Wire `qco[Method -> {"Braket", ...}]` to `BraketJobSubmit` at the same dispatch point the `"Qiskit"` method uses; keep `BraketJob`'s accessor names identical to `IBMJob`'s so the implicit provider protocol holds.

**Stage 5: Documentation.** Rewrite the AWS section of the tutorial (Section 10).

**Public surface bookkeeping** (Stages 1, 3): register `BraketDevice`, `BraketJobSubmit`, `BraketJob` in `PacletInfo.wl` `Symbols`; write each `::usage` in the doc-page Usage section (never hand-edit the auto-generated `Usage.m`).

---

## 10. Test strategy

- **Offline (the bulk; no creds, runs in CI):** a `Tests/Braket.wlt` `VerificationTest` suite over recorded fixtures.
  - *Transport:* assert `prepareBraketParams` produces the one-list form, injects `ClientToken` when absent and preserves it when given, and renders `action`/`deviceParameters` as string-keyed Associations (never JSON strings). A mock service object captures and returns the payload.
  - *Capability parser:* against SV1 / IonQ / Rigetti `GetDevice` fixtures, assert the normalized record fields (gate set, shots range, region, cost, status) and graceful `Missing[...]` on a fixture with a field removed (schema-drift simulation).
  - *Compiler:* assert `gphase` is stripped, logical mode emits only named gates from the device's set (no `U`/modifier syntax, no `cphaseshift` for IonQ), and native mode wraps in `#pragma braket verbatim box { ... }` only when `VerbatimAvailable`.
  - *Decoder:* against SV1 (`measurements`) and IonQ (`measurementProbabilities`) `results.json` fixtures, assert identical normalized `Counts`; the known-asymmetric-state fixture locks bit order.
- **Live (gated behind credentials + an explicit opt-in; not in default CI):**
  - *SV1 round-trip* (cheap, noiseless): the integration smoke test. Submit, poll to `COMPLETED`, decode, and confirm the finite-shot distribution matches the QF-simulated exact distribution within sampling error. Also the bit-order ground truth (Section 7E).
  - *IonQ task* (costs money, behind explicit opt-in): the QPU validation. Confirm the full async lifecycle (queue position visible, `COMPLETED`, decoded counts) and that the overlay against the exact WL distribution shows device noise. Tag the task so it is identifiable on the AWS bill.
  - Credentials are requested from the user at run time; any key shared earlier in this engagement is considered compromised and must be rotated. (None are needed for this design document or for the offline suite.)

---

## 11. Documentation plan

The tutorial `Documentation/English/Tutorials/QPUServiceConnect.md` (a TechNote, authored as literate Markdown and built with md2nb) currently demonstrates the **dead** blog pipeline in its "Amazon Braket" section (lines 187-336). Rewrite that section around the new API, mirroring the IBM section's structure:

1. `ServiceConnect["AWS"]` -> active connection (keep, it is correct).
2. **Device discovery** via `SearchDevices` / `BraketDevice`, showing `dev["SupportedGates"]`, `dev["ShotsRange"]`, `dev["Status"]` instead of the manual JSON decode.
3. **Submission** via `BraketJobSubmit[qc, dev, "Shots" -> n]` returning a `BraketJob`, replacing the hand-built `CreateQuantumTask` with bare-symbol action keys.
4. **Result retrieval** via `job["Refresh"]` / `job["Measurement"]` / `job["Counts"]`, replacing the manual `GetQuantumTask` -> S3 `GetObject` -> parse.
5. Keep the honest physics payoff: overlay SV1's finite-shot sample on the exact WL distribution (the sampling spread), and note that a QPU run adds device noise on top, making it a cross-technology comparison.

Live cloud cells are marked non-evaluating (`eval:false` in the md twin, as the IBM live cells are); offline cells must be verifier-clean. The deleted blog code's gotchas (bare trailing rules, missing `ClientToken`, bare-symbol `action`, `"Provider" -> "AWSBraket"`, manual S3) do not reappear anywhere. Per the doc-build conventions already recorded for this tutorial: `<!-- => -->` output hints go *after* the fence; `Off[General::shdw]` for the QF doc build; verify the built `.nb` with `Get`, and round-trip through `NotebookToMarkdown`.

A separate short Symbol page each for `BraketJobSubmit`, `BraketJob`, `BraketDevice` (md2nb `Template: Symbol`), following the `IBMJobSubmit` / `IBMJob` pages as the model.

---

## 12. Open questions and coordination

- **Compilation engine timing.** The native WL engine (the strategic default) depends on the sibling `QF-Native-QPU-QASM-Investigation.md`. Until it lands, ship on the qiskit fallback (post Stage 0). The `BraketCompile` seam makes the swap a contained change. Confirm with that investigation that the mode semantics here (logical vs native, gphase stripping, capability-driven gate set, verbatim box) match what its emitter will produce.
- **`BraketDevice` public or scoped.** Recommended public (it is a useful inspection/discovery object), but it could start PackageScope and be promoted. Decide before registering symbols in `PacletInfo`.
- **Cost fidelity.** `"Cost"` as a pre-submission estimate from `deviceCost` is the honest best available (Braket returns no dollar figure). If a more precise post-hoc number is ever wanted, it would come from AWS Cost Explorer / Budgets, out of scope here.
- **Default device.** `Automatic` -> SV1 is the safe, cheap default; confirm this is the desired behavior versus requiring an explicit device for any submission.
- **DeviceParameters.** Optional `deviceParameters` (also `jsonvalue`) carries vendor-specific run options (e.g. error mitigation, gate-level control). The transport shim already handles its encoding; surfacing specific options is a follow-on once a concrete need appears.

---

## Appendix A: Braket API reference (from `service.json`, `2019-09-01`)

Operations: `CreateQuantumTask`, `GetQuantumTask`, `CancelQuantumTask`, `GetDevice`, `SearchDevices`, `SearchQuantumTasks` (plus Hybrid-Jobs and spending-limit ops, out of scope).

**`CreateQuantumTask`** (`POST /quantum-task`, 201). Required: `clientToken` (`String64`), `deviceArn`, `shots` (`Long`), `outputS3Bucket`, `outputS3KeyPrefix`, **`action`** (`JsonValue`, **`jsonvalue`**). Optional: **`deviceParameters`** (`jsonvalue`), `tags`, `jobToken`, `associations`. Returns `quantumTaskArn`.

**Corrected `action` payload** (string keys, passed as an Association, never pre-stringified):
```
"action" -> <|
   "braketSchemaHeader" -> <|"name" -> "braket.ir.openqasm.program", "version" -> "1"|>,
   "source" -> qasmString
|>
```

**`GetQuantumTask`** (`GET /quantum-task/{quantumTaskArn}`, 200). Returns `status`, `failureReason`, `deviceArn`, `deviceParameters`, `shots`, `outputS3Bucket`, **`outputS3Directory`**, `createdAt`, `endedAt`, `queueInfo` (`queue`, **`position`**, `queuePriority`, `message`), `numSuccessfulShots`, `actionMetadata`, `tags`, `jobArn`.

**`CancelQuantumTask`** (`PUT /quantum-task/{quantumTaskArn}/cancel`, 200). Required: `quantumTaskArn`, **`clientToken`** (required again). Returns `cancellationStatus` (`CANCELLING` | `CANCELLED`).

**`GetDevice`** (`GET /device/{deviceArn}`, 200). Returns `deviceName`, `providerName`, `deviceType` (`QPU` | `SIMULATOR`), `deviceStatus` (`ONLINE` | `OFFLINE` | `RETIRED`), **`deviceCapabilities`** (`jsonvalue`: the nested capability schema with `supportedOperations`, `supportedModifiers`, `supportedPragmas`, `connectivity`, `shotsRange`, `deviceCost`, `executionWindows`, native-gate detail), `deviceQueueInfo`.

`QuantumTaskStatus` enum: `CREATED`, `QUEUED`, `RUNNING`, `COMPLETED`, `FAILED`, `CANCELLING`, `CANCELLED`.

## Appendix B: IBM template anchors (the mirror)

`IBMQuantum.m`: `IBMJobSubmit` (`:126`), `ibmActiveConnection` (`:29`), `iMethodIBMJob` blocking constructor (`:231`), `Refresh` (`:268`), `ibmRaw`/`ibmFetch` lossless raw fetch (`:247-256`), `$ibmTerminal` (`:263`), `$ibmProperties` (`:518`), `$ibmActions` (`:532`), unknown-property catch-all (`:541`), summary box (`:569`). `Qiskit.m`: `qiskitInitBackend` AWS branch with the backend-selection bug (`:436-441`), AWSBraket OpenQASM conversion (`:648-652`), `qiskitReorderByQubit` bit-order map (`:484`). `QuantumQASM.m`: native emitter `qasmEmitCircuit` (`:336`), qubit-only guard (`:321`), provider-export routing (`:367`).
