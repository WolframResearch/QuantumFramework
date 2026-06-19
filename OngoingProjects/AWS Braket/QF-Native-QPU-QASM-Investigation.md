# Why QuantumFramework Cannot Yet Emit QPU-Runnable OpenQASM for AWS Braket, and What It Would Take

**Investigation only. No kernel changes were made; this is a written report.**

| Field | Value |
|---|---|
| Date | 2026-06-18 |
| Repo | `/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework` |
| HEAD | `c3f21a43` (audit anchor `a601a187`; the QASM/Qiskit/IBM kernel files are **untouched** between the two, so all audit citations below hold at HEAD: `git diff a601a187..HEAD -- QuantumFramework/Kernel/QuantumCircuitOperator/` is empty) |
| Paclet version | `Wolfram/QuantumFramework 2.0.0` |
| Verified offline | native QASM emission (live kernel), the 1- and 2-qubit decomposition machinery, the qiskit-bridge logical-gate conversion via QF's own `"braket"` Python env |
| Needs live AWS (NOT done here) | end-to-end IonQ/SV1 acceptance, `GetDevice` capability JSON. The Braket parse error, the `cphaseshift` rejection, and the IonQ `supportedOperations` list quoted below are from the live AWS audit of 2026-06-18, not reproduced here. |
| Security | Any AWS access key shared in an earlier session should be treated as compromised and rotated. |

---

## Big picture

A quantum processor does not need much: it accepts a circuit written in OpenQASM, in (or reducible to) the gate alphabet its hardware actually implements. The claim under investigation, "QuantumFramework cannot currently produce QPU-runnable code for AWS Braket," is **true today**, but the obstruction is shallow plumbing, not a missing capability.

QF has two ways to turn a circuit into OpenQASM, and both miss Braket for distinct, fixable reasons:

1. **The native Wolfram Language emitter** writes a *dialect* of OpenQASM 3 built around the primitive gate $U(\theta,\phi,\lambda)$, the gate **modifiers** `ctrl @` / `negctrl @`, and a bare `gphase`. Braket's OpenQASM parser is a restricted subset of the language and rejects this syntax outright, before any device or gate set is even considered. This is a *grammar* mismatch.

2. **The qiskit bridge** can produce exactly the named-gate OpenQASM Braket accepts (verified offline below), but a one-line backend-selection defect routes every named device to a nonexistent simulator, and the only path that does run targets the SV1 statevector simulator's basis (which contains `cphaseshift`), a gate real IonQ hardware rejects.

Neither failure is fundamental. The physics content of "compile to a device's native set" is unitary synthesis, and QF already owns the relevant machinery: it decomposes any 1-qubit gate to an Euler $ZYZ$ form and any 2-qubit gate to a Cartan/KAK circuit over $\{U, \mathrm{CX}\}$ (both verified live). What is missing is (a) an emitter that speaks Braket's accepted named-gate dialect, (b) a step that rewrites the synthesized $\{U,\mathrm{CX}\}$ into the *device's* advertised gate set, and, only for hardware with restricted qubit connectivity, (c) a router. For an all-to-all device such as IonQ, items (a) and a trivial $U \to \{r_z, r_y\}$ rewrite are essentially the whole job.

The remainder of this report traces the exact code paths, separates "the minimal patch" from "the real fix," and assesses the feasibility of a native WL transpiler that bypasses the qiskit/boto3 round-trip entirely.

---

## 1. The two QASM paths, and exactly where each one dies

### 1a. Native WL path: `QuantumQASM[qco]` produces a dialect Braket cannot parse

The default export is pure Wolfram Language (no Python, no qiskit):

- `QuantumQASM[qco_QuantumCircuitOperator]` $\to$ `qasmEmitCircuit` (`QuantumQASM.m:336`, worker at `QuantumCircuitOperator/Properties.m:476`).
- The header is fixed at `QuantumCircuitOperator/Properties.m:518`: `OPENQASM 3.0;` then `qubit[N] q;` and `bit[M] c;`.
- Each gate is serialized by `qasmEmitOperator` / `qasmEmitSimple` (`QuantumOperator/Properties.m:944-1068`).

The emitter's gate vocabulary is the OpenQASM-3 *primitive* $U$ gate plus modifiers:

- **Single-qubit gate** $\to$ `U(a, b, c) q[i];` (`QuantumOperator/Properties.m:966-970`, applied at `:1067-1068`). The angles are the Euler $ZYZ$ angles of the gate; a leftover global phase is dropped.
- **Controlled gate** $\to$ `ctrl(n) @ negctrl(m) @ U(...) q[...] q[...];` followed, when the relative phase is nonzero, by `ctrl(n) @ negctrl(m) @ gphase(...) ...;` (`QuantumOperator/Properties.m:1000-1047`, specifically the prefix at `:1034` and the gate/phase lines at `:1036-1040`).
- **Global phase** $\to$ `gphase(phi);` (`QuantumOperator/Properties.m:1050-1051`).
- **SWAP / permutation** $\to$ `swap q[i] q[j];` (`:972-977`, `:1055-1065`).

Verified live against the working tree (Bell circuit `{H, CNOT}`):

```
OPENQASM 3.0;
qubit[2] q;
bit[0] c;
U(1.5707963267948966, 0., 3.141592653589793) q[0];
ctrl(1) @ negctrl(0) @ U(3.141592653589793, 0., 3.141592653589793) q[0] q[1];
```

and for a circuit carrying a `T` and a `CZ`:

```
OPENQASM 3.0;
qubit[3] q;
bit[0] c;
U(1.5707963267948966, 0., 3.141592653589793) q[0];
U(0., 0., 0.7853981633974483) q[0];
ctrl(1) @ negctrl(0) @ U(0., 0., 3.141592653589793) q[0] q[1];
ctrl(1) @ negctrl(0) @ U(3.141592653589793, 0., 3.141592653589793) q[1] q[2];
```

**Why Braket rejects this.** From the live AWS audit, feeding this string to Braket's parser yields:

```
Parsing error: ... missing ';' at 'q'
```

Braket implements a restricted subset of OpenQASM 3. It does not accept the bare primitive-$U$ gate call, the `ctrl @` / `negctrl @` gate modifiers, or `gphase` in a plain (non-`verbatim`) program. The IonQ device capabilities confirm the modifier half directly: `supportedModifiers: {}` (none). So the parse failure is a **dialect** failure: it happens regardless of which device you target and regardless of the gate alphabet. The native path never even reaches the gate-set question.

Two further facts about this path, both correct by design and not the problem here:

- It is **qubit-only**. A non-2-dimensional wire returns a clean `Failure["QuantumQASM", ...]` (`QuantumQASM.m:321-331`); verified live on a qutrit. OpenQASM has no qudit representation, so this is the right behavior, not a bug.
- When it cannot serialize a gate it emits a `// Unimplemented` marker that `qasmEmitCircuit` turns into a `Failure` naming the gate and suggesting a native gate set (`QuantumCircuitOperator/Properties.m:502-515`).

### 1b. Qiskit-bridge path: `qc["Qiskit"]["QASM", "Provider" -> "AWSBraket"]` is broken by a backend-selection defect and a wrong basis

Any *extra argument* to `QuantumQASM` routes export through qiskit (`QuantumQASM.m:347-349`). A `"Provider"` / `"Backend"` request specifically engages the hardware path (`QuantumQASM.m:367-372`), which calls `qiskitQASM` (`Qiskit.m:639`). The chain is:

`qiskitQASM` (`Qiskit.m:639`) $\to$ `qiskitInitBackend` (`Qiskit.m:355`, boots the provider/backend) $\to$ `transpile(circuit, backend, optimization_level=...)` then, for an AWS provider, `convert_qiskit_to_braket_circuit(circuit).to_ir(IRType.OPENQASM).source` (`Qiskit.m:644-661`, conversion at `:648-652`).

The provider name routing is at `Qiskit.m:350-351`:

```wolfram
$IBMProvider = "IBM" | "IBMQ" | "IBMProvider" | "IBMRuntime"
$AWSProvider = "AWS" | "AWSBraket" | "Braket"
```

**The defect.** Inside `qiskitInitBackend`, the per-provider backend selection (`Qiskit.m:420-466`) handles the AWS branch at `Qiskit.m:435-440`:

```python
    $AWSProvider,
"
if backend_name is None:
    backend = provider.get_backend('SV1')
else:
    backend = provider.get_backend('basic_simulator')
",
```

The `else` branch hardcodes the string literal `'basic_simulator'` and **discards `backend_name` entirely**. So `"Backend" -> "Forte Enterprise 1"` (or any real device) becomes a request for a backend literally named `basic_simulator`, which does not exist in the AWS Braket provider, producing "No backend matches the criteria." Contrast the sibling branches, all of which honor `backend_name`:

- IBM (`Qiskit.m:421-426`): `backend = provider.backend(backend_name)`.
- GenericV2 fake device (`Qiskit.m:446-449`): `GenericBackendV2(num_qubits=int(backend_name) ...)`.
- Default/Aer (`Qiskit.m:451-465`): uses the circuit, no named backend.

This is a copy-paste defect: the AWS `else` line should read `backend = provider.get_backend(backend_name)`. (The qiskit_braket_provider uses space-separated device names, e.g. `provider.get_backend('Forte Enterprise 1')`, so the surface name must be passed through verbatim, which it is, once the line stops ignoring it.)

**The consequence chain.** Because every named backend fails, the *only* AWS path that runs is the no-backend path, `get_backend('SV1')`. `transpile` then compiles to SV1's basis, which includes `cphaseshift`. Submitting that to IonQ hardware gives (live AWS audit):

```
uses a gate: cphaseshift which is ... not supported by the device
```

So even the one path that produces a string produces the wrong gate set for hardware. The IonQ Forte-Enterprise-1 `supportedOperations` (from `GetDevice` $\to$ `DeviceCapabilities.action."braket.ir.openqasm.program".supportedOperations`) are:

```
{x, y, z, h, s, si, t, ti, v, vi, rx, ry, rz, cnot, swap, xx, yy, zz}
```

with `supportedModifiers: {}` and `supportedPragmas` including `verbatim`. `cphaseshift`, `gphase`, gate modifiers, and a bare `U` outside a `verbatim` box are all absent.

**A third obstruction:** instantiating `AWSBraketProvider()` (`Qiskit.m:402-406`) boots boto3 and requires a region plus credentials (`AWS_DEFAULT_REGION`, `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`) *merely to emit text*. There is no offline AWS QASM path through this bridge.

---

## 2. Is the minimal fix just correcting the backend name?

Correcting `Qiskit.m:440` to `provider.get_backend(backend_name)` is **necessary but not sufficient.**

What it buys: `"Backend" -> "<device>"` would reach the real device, and `transpile(circuit, backend, optimization_level=...)` (`Qiskit.m:646`) would then be genuinely device-aware, compiling against that device's `Target` (its basis gates, coupling map, and, since the error-aware-transpile work shipped at `f370f568`, its per-instruction error/duration). For an all-to-all device with a logical basis, this is most of the way to runnable QASM.

What it still does not handle:

- **Residual global phase.** The `convert_qiskit_to_braket_circuit` output still carries a trailing `gphase(...)` line (verified offline; see §4). IonQ rejects it ("GPhase is not supported"). The phase is physically irrelevant to measurement, so the fix is to strip it, but the backend-name patch alone does not.
- **Native vs logical gates.** Transpiling against a real IonQ `Target` can yield the vendor's *native* gates (`gpi`, `gpi2`, `ms`/`zz`), which are only legal inside a `#pragma braket verbatim` box; or it can yield *logical* gates Braket recompiles for you. Which one you get depends on how qiskit_braket_provider configures the target. The robust, verbatim-free route is to transpile to an explicit logical basis (`basis_gates = {h, cx, rz, ry, rx, swap}`), all of which are in IonQ's `supportedOperations`.
- **Credentials.** The provider still boots boto3; you still need live AWS to emit anything.

So the one-line patch unblocks the device-aware qiskit route (and, as a bonus, `QiskitTarget["FromBackend"]` for AWS, see §5), but a runnable IonQ program also needs gphase stripping, a logical-basis choice, and credentials. And it keeps QF permanently tethered to the qiskit + boto3 stack for hardware export.

---

## 3. The deeper question: a native WL transpiler with no qiskit/boto3 round-trip

**Verdict: feasible, and for all-to-all devices it is a small project.** The hard part of "compile to a device's gate set" is unitary synthesis, and QF already has the pieces.

### 3.1 What QF already owns (verified live)

- **1-qubit synthesis.** `qo["ZYZ"]` (`QuantumOperator/Properties.m:691`) returns a $U(\theta,\phi,\lambda)$ gate (Euler $ZYZ$) plus a recovered `GlobalPhase`. The Euler-angle engine is `UnitaryEulerAngles` / `UnitaryAnglesWithPhase` (`:665-684`).
- **2-qubit synthesis.** `qo["KAK"]` (`QuantumOperator/Properties.m:869`, numeric core `TwoQubitKAK` at `:752`) factors a 4x4 unitary in the magic (Bell) basis and returns a `QuantumCircuitOperator` of $U$ gates and CX gates plus a `GlobalPhase`. Verified live on `CNOT` and on a random $SU(4)$: both produce a circuit whose gate labels are `U[...]` and $\mathrm{C}_X$. The gate-emitting helper is `kakUGate` (`:837`), which itself calls `UnitaryAnglesWithPhase`.
- **The umbrella.** `qo["Decompose", gateset]` (`QuantumOperator/Properties.m:938-942`) routes 1-qubit $\to$ `ZYZ`, 2-qubit $\to$ `KAK`. **Caveat:** the `gateset` argument (default `{"U","CX"}`) is *accepted but ignored*; the property always emits $\{U, \mathrm{CX}\}$. So "decompose to an arbitrary basis" is not yet real; a native transpiler would take the $\{U, \mathrm{CX}\}$ output and rewrite to the device set itself. ($n \ge 3$-qubit gates, the Quantum Shannon Decomposition, are unimplemented and fall through; in practice circuits are built from 1- and 2-local gates, and QF already ships multi-qubit named gates like Toffoli as pre-decomposed 1-/2-qubit circuits.)
- **A named-gate vocabulary.** The native importer's `$qasmGateTable` (`OpenQASMImport.m:160-200`) already maps `h, x, y, z, s, t, sx, rx, ry, rz, p, u, swap, cx/cnot, cz, ...` to QF operators. The *inverse* of this table (QF label $\to$ named QASM gate) is most of a named-gate emitter; the table proves QF already knows the names.

### 3.2 What a native Braket emitter needs, and the effort

1. **A Braket-dialect string emitter.** Emit *named* gates with comma-separated operands (`h q[0];`, `cnot q[0], q[1];`, `rz(theta) q[0];`), with **no** bare $U$, **no** `ctrl/negctrl` modifiers, **no** `gphase`. This is a formatting layer parallel to the existing `qasmEmit*` functions, reusing the inverted `$qasmGateTable`. The one real rewrite is $U(\theta,\phi,\lambda) \to r_z(\phi)\, r_y(\theta)\, r_z(\lambda)$ (up to global phase), since `ZYZ`/`KAK` produce $U$ gates and $\{r_y, r_z\}$ are in IonQ's set. **Effort: low.**

2. **Device-gate-set targeting.** Take the device's `supportedOperations` (from `GetDevice`, §5) and decompose to it. For an IonQ *logical* basis $\{r_x, r_y, r_z, \mathrm{cnot}, \ldots\}$: `KAK` + `ZYZ` already land in $\{U, \mathrm{CX}\}$; rewrite $U \to r_z r_y r_z$ and $\mathrm{CX} \to \mathrm{cnot}$. **Effort: medium.** For the *native* basis $\{\mathrm{gpi}, \mathrm{gpi2}, \mathrm{ms}/\mathrm{zz}\}$: needs a gpi/gpi2 Euler synthesis (a $\mathrm{gpi2}$-$r_z$-$\mathrm{gpi2}$ sandwich realizes any 1-qubit rotation) and MS/ZZ as the 2-qubit primitive. Textbook, but new code. **Effort: medium-high.**

3. **Global-phase stripping.** The 0-qudit `GlobalPhase` operators (`qasmEmitOperator` at `QuantumOperator/Properties.m:1050-1051`) are simply omitted, or emitted only when a controlled context makes the relative phase physical. **Effort: trivial.**

4. **Verbatim box.** Native vendor gates are legal only inside a `#pragma braket verbatim` / `box { ... }`. The emitter would wrap the native-gate sequence and emit the pragma, gated on the device advertising `verbatim` in `supportedPragmas` (IonQ does). Needed *only* for the device-native route; the logical route needs no verbatim. **Effort: low, once the native gates exist.**

5. **Connectivity / coupling.** This is the one genuinely missing capability: QF has **no router** (no SABRE-style layout/SWAP-insertion). For all-to-all hardware (IonQ) it is a non-issue; any logical qubit talks to any other. For limited-connectivity superconducting devices (Rigetti, IQM) it requires SWAP insertion against the coupling map. This is exactly the piece the qiskit bridge gives for free. **Effort: medium-high; can be deferred** by only offering the native route for all-to-all devices (or for `verbatim` where the user has pre-mapped qubits), and falling back to the qiskit bridge otherwise.

**Net.** For IonQ-class all-to-all devices via the logical route, a native WL Braket emitter is items 1 + 3 plus the already-existing $\{U,\mathrm{CX}\}$ decompose. That is a small, self-contained addition with no Python dependency. The native-in-verbatim route and the connectivity router are larger and device-class-specific, and can come later.

---

## 4. Offline confirmation of the logical-gate bypass (no AWS)

To anchor the recommendation in evidence rather than in qiskit's documentation, I exercised the conversion offline through QF's own `"braket"` Python environment (`PythonTools.m:16`, which provisions `qiskit-braket-provider` on a Python 3.11 session). No provider is booted and no AWS call is made: this is a pure local `transpile` to a logical basis followed by the circuit-format conversion.

For `QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "T" -> 1}]`, transpiled to `basis_gates = {h, cx, rz, ry, rx, swap}` and converted with `convert_qiskit_to_braket_circuit(tc).to_ir(IRType.OPENQASM).source`:

```
OPENQASM 3.0;
bit[2] b;
qubit[2] q;
h q[0];
cnot q[0], q[1];
rz(0.7853981633974483) q[0];
gphase(0.39269908169872414);
b[0] = measure q[0];
b[1] = measure q[1];
```

Three things this confirms:

- **The accepted shape uses named gates with comma-separated operands** (`h q[0];`, `cnot q[0], q[1];`, `rz(...) q[0];`), in direct contrast to the native WL emitter's `U(...)` / `ctrl(1) @ negctrl(0) @ U(...)` dialect from §1a. This is the target a native WL Braket emitter must reproduce.
- **A `gphase(...)` line survives the conversion** and must be stripped for IonQ. This independently reproduces the user's finding 5 (strip the trailing `gphase`).
- **Every emitted gate (`h`, `cnot`, `rz`) is in IonQ's `supportedOperations`.** With the `gphase` line removed, this circuit is, by inspection, in IonQ's logical alphabet. (End-to-end acceptance by the live device still needs credentials; not done here.)

The user's verified end-to-end recipe (finding 5) is the same idea: in the `"braket"` env, `transpile(circuit, basis_gates=['h','cx','rz','ry','rx','swap'], optimization_level=1)`, then `convert_qiskit_to_braket_circuit(...).to_ir(IRType.OPENQASM).source`, then strip the trailing `gphase(...)` line. That bypasses QF's emitter and its provider plumbing entirely, which is why it works.

---

## 5. Discovering the device gate set dynamically (so it survives new devices/vendors)

Hardcoding `{x, y, z, h, ...}` would rot the moment a vendor adds a gate or AWS onboards a new device. The gate set must be read from the device.

**Source of truth.** AWS Braket `GetDevice` returns `DeviceCapabilities` as JSON; the OpenQASM action carries everything needed:

- `action."braket.ir.openqasm.program".supportedOperations` $\to$ the basis gates to decompose to.
- `... .supportedModifiers` $\to$ whether `ctrl`/`negctrl` may be emitted (empty for IonQ, so do not).
- `... .supportedPragmas` $\to$ whether `verbatim` is available (it is, for IonQ).
- `paradigm.connectivity` (`fullyConnected` flag or an explicit adjacency) $\to$ feeds the router / decides whether routing is even needed.

**Getting it without the qiskit round-trip.** Two options:

1. **A direct WL `GetDevice` client.** The capability blob is plain JSON behind an authenticated Braket control-plane API call. A small WL client (boto3 via `ExternalEvaluate`, or `URLExecute` with SigV4 signing) returns it, parsed natively into an Association like `<|"BasisGates" -> {...}, "Modifiers" -> {...}, "Pragmas" -> {...}, "Connectivity" -> ...|>`. This still needs AWS credentials (it is authenticated), but it needs **no qiskit and no transpile**, and it is the canonical source.

2. **Via the qiskit `Target` (already in QF, but currently broken for AWS).** `QiskitTarget["FromBackend" -> name, "Provider" -> "AWSBraket"]` (`Qiskit.m:671`) reads any backend's populated `Target` into a QF spec carrying `BasisGates`, `CouplingMap`, and `InstructionProperties` (`QuantumQASM.m:237-243`). It is the natural discovery hook, **but it calls `qiskitInitBackend` with `"Backend" -> name` (`Qiskit.m:673`), so it hits the identical AWS backend-selection defect from §1b**: for a named AWS device the bug substitutes the wrong backend. The offline `"GenericV2"` provider works (it honors `backend_name`); a named AWS device does not, until `Qiskit.m:440` is fixed. Note also (from the traps investigation, carried in the audit) that `QiskitTarget["FromBackend" -> "GenericV2"]` *without* a `"Provider"` silently yields an idealized Aer-style target rather than erroring.

**Recommendation.** Build option 1 (a thin WL `GetDevice` capability client), cache the result per device-ARN, feed `BasisGates` into the native decompose target and `Connectivity` into the (future) router. Hardcode nothing; keep a static per-vendor table only as an offline fallback. Once `Qiskit.m:440` is fixed, option 2 becomes a second, qiskit-backed discovery path that reuses the shipped `QiskitTarget` machinery.

---

## 6. Compact call-graph (for reference)

Native path (the dialect Braket rejects):

```
QuantumQASM[qco]                                  QuantumQASM.m:336
  -> qasmEmitCircuit[qco]                          QuantumCircuitOperator/Properties.m:476
       header "OPENQASM 3.0; qubit[N] q; ..."      :518
       per operator -> op["QASM"] -> qasmEmitOperator / qasmEmitSimple
                                                    QuantumOperator/Properties.m:944-1068
         1q gate   -> "U(a,b,c) q[i];"              :966-970, :1067-1068
         ctrl gate -> "ctrl(n) @ negctrl(m) @ U(...) ...;"  :1000-1047 (prefix :1034)
         phase     -> "gphase(phi);"                :1050-1051
  => Braket parser: "Parsing error: ... missing ';' at 'q'"  (restricted grammar; no U / modifiers / gphase)
```

Qiskit-bridge path (backend defect + wrong basis):

```
QuantumQASM[qc_QiskitCircuit, ... "Provider"->"AWSBraket" ...]   QuantumQASM.m:347-372
  -> qiskitQASM[qc, "Provider"->..., "Backend"->...]              Qiskit.m:639
       -> qiskitInitBackend                                       Qiskit.m:355
            $AWSProvider -> AWSBraketProvider()  (needs region+creds)  :400-406
            backend select (AWS branch)                           :435-440
               backend_name is None -> get_backend('SV1')         :438   (only path that runs)
               else                 -> get_backend('basic_simulator')  :440  *** BUG: ignores backend_name ***
       -> transpile(circuit, backend, optimization_level=1)       :646
       -> convert_qiskit_to_braket_circuit(circuit).to_ir(OPENQASM).source   :648-652
  => named device: "No backend matches the criteria"
  => SV1 path: basis includes cphaseshift -> IonQ "uses a gate: cphaseshift which is ... not supported"
  => output still carries a gphase(...) line  (IonQ "GPhase is not supported")
```

---

## 7. Verdict and recommended sequence

The claim is true today, for two independent and shallow reasons (a wrong-dialect native emitter; a one-line backend defect plus a wrong-basis fallback in the qiskit bridge). A QPU does indeed only need OpenQASM in its native set, and QF has the synthesis machinery to produce it. Recommended order, smallest lever first:

1. **One-line fix** `Qiskit.m:440` $\to$ `provider.get_backend(backend_name)`. Unblocks device-aware qiskit transpile *and* `QiskitTarget["FromBackend"]` for AWS. (Low effort, high leverage.)
2. **Strip `gphase`** in the AWS export (in `qiskitQASM`, or as a post-filter on the converted source). Makes IonQ accept the qiskit-route output.
3. **Native Braket-dialect emitter**: named gates, comma operands, no modifiers/gphase, with a $U \to r_z r_y r_z$ rewrite, targeting a supplied `basis_gates`. Removes the boto3/qiskit dependency for the logical route on all-to-all devices. (Builds on the existing `ZYZ`/`KAK` decompose and the `$qasmGateTable` names.)
4. **WL `GetDevice` capability client** for dynamic gate-set + connectivity discovery, cached per device-ARN.
5. **Optional, larger**: native-gate (`gpi`/`gpi2`/`ms`) synthesis + `verbatim` wrapping for the device-native route, and a connectivity router (SWAP insertion) for limited-connectivity hardware.

## What still needs AWS (not done in this investigation)

- End-to-end acceptance of a QF-produced string by a live IonQ (and SV1) device.
- The exact `verbatim` box syntax and the per-device native gate list, which should be confirmed against a live `GetDevice` rather than assumed.
- The Braket parse error, the `cphaseshift` rejection, and the IonQ `supportedOperations` list are quoted from the live AWS audit of 2026-06-18; everything in sections 1a, 3, and 4 was reproduced offline against the working-tree paclet.
