# Design proposal: a version-resilient Qiskit transpile + faithful QASM export layer for QuantumFramework

> ## IMPLEMENTED (2026-06-05) — final form supersedes the QiskitCircuit-method design below
>
> The shipped feature is a single new **top-level function `QuantumQASM`** (plus a
> `QiskitTarget` device-spec object), authored in its own kernel file, NOT methods on
> `QiskitCircuit`. The existing `QiskitCircuit["QASM3"|"QASM2"|"Transpile"]` methods and
> their bugs were left untouched (the change is purely additive; `QuantumQASM` bypasses
> them with its own clean transpile + faithful dump).
>
> - **Surface:** `QuantumQASM[qco, basisGates : Automatic, opts___Rule] -> OpenQASM String | Failure`.
>   Faithful dump by default; transpiles to a native basis only when a basis / `Target` /
>   `CouplingMap` / any transpile option is given. `"Version" -> 3|2` selects the OpenQASM
>   version. Also accepts a `QiskitCircuit` as input.
> - **`QiskitTarget[spec_Association]`:** validated, normalized device-spec carrier (NOT a
>   pickled Target); the live qiskit `Target` is rebuilt from the spec at point-of-use
>   (no cross-version pickle), entry points existence-guarded, `measure`/`reset`
>   auto-added per qubit. Accessors `["Spec"|"NumQubits"|"OperationNames"|"CouplingMap"]`
>   + a summary box.
> - **Version resilience:** options are forwarded as `**kwargs` to `qiskit.transpile`,
>   validated at call time against the live `inspect.signature`; CamelCase keys are mapped
>   to qiskit snake_case by an *algorithmic* transform (case-tolerant; raw snake_case also
>   passes through); only `coupling_map`/`target` are coerced. Not a thin wrapper.
> - **Files:** `QuantumFramework/Kernel/QuantumCircuitOperator/QuantumQASM.m` (auto-loaded),
>   tests in `Tests/Qiskit.wlt`.
> - **Verified:** `Tests/Qiskit.wlt` 24/24 (incl. string-for-string cross-check vs raw
>   qiskit 2.4.1 and a unitary-equivalence correctness test); existing
>   `Tests/QuantumCircuitOperator.wlt` 37/37 (no regression).
>
> The sections below are the original investigation/design that led here; the
> mechanism (pass-through kwargs + live-signature validation + algorithmic CamelCase +
> coercion) is what was implemented, on a top-level function rather than a method.

Status (original): investigation + design. Prototype in `prototype.wl`, validated 9/9
by `validate.wl` against the live paclet and qiskit 2.4.1.

**Decisions made (2026-06-04), see §7:** (1) `["QASM3"]`/`["QASM2"]` flip to
faithful-dump default. (2) Option surface is **full CamelCase WL-style**, mapped to
qiskit snake_case by an *algorithmic* transform (so resilience is preserved, see
§3.1/§3.2). (3) This session stays design-only; no kernel implementation yet.

Author context: prompted by the verified end-to-end "QF circuit ->
OpenQASM3 -> hardware" path (IQCC `arbel`, `openqasm2qua` compiler). The
transpile step (decompose to a native gate set, route onto a coupling map) is
the piece that belongs in QF. The reference investigation lives at
`QuantumMachines/QUALink/Claude/PLAN-QF-to-Arbel-Gate-Pipeline.md`.

Repo state read for this design: `main` @ `1a810bfc`. The Qiskit bridge has not
changed since the audit anchor (`bfddf5ed`); every line citation below was
re-verified against live HEAD.

---

## 0. TL;DR (the recommendation in five lines)

1. **Mechanism for version resilience:** forward a WL `Association` of qiskit
   `transpile` parameters straight through as Python `**kwargs`, and validate the
   keys at call time against `inspect.signature(transpile)` read live from the
   installed qiskit. Nothing about the option list is hardcoded in WL. A small
   coercion layer wraps the handful of rich-typed keys (`coupling_map`, `target`).
2. **Shape:** compositional, not monolithic. `qc["Qiskit"]["Transpile", optsAssoc]`
   returns a `QiskitCircuit`, so `["Ops"]`, `["Depth"]`, `["Diagram"]`,
   `["QuantumCircuitOperator"]`, `["QASM3"]` all compose on the result. Option keys
   are **CamelCase WL-style** (`"BasisGates"`, `"CouplingMap"`, `"OptimizationLevel"`),
   converted to qiskit snake_case by an algorithmic CamelCase->snake_case transform so
   no per-param dictionary is maintained and future qiskit params auto-map.
3. **Faithful export:** make `QiskitCircuit["QASM3"]` / `["QASM2"]` dump the
   circuit **as-is** (no transpile-against-AerSimulator), so an already-native
   circuit stays native. Re-engage the backend path only when a `Backend`/`Provider`
   is explicitly given (preserves the AWS-Braket and IBM workflows).
4. **Target model:** add a `QiskitTarget[spec]` compound object that forwards to
   `Target.from_configuration` with the same kwargs-forwarding philosophy; it is the
   version-stable, highest-value device interface and is reusable across circuits.
5. **Migration:** fix `["Transpile"]` and `["QASM3"]`/`["QASM2"]` in place (both
   are bugfixes); keep the old behavior reachable behind explicit `Backend`/`Provider`.
   Add `VerificationTest`s + a CHANGELOG note flagging the `["QASM3"]` default change.

---

## 1. Current-state audit (what exists, and the four gaps)

File: `QuantumFramework/Kernel/QuantumCircuitOperator/Qiskit.m` (680 LOC). Python
session plumbing in `Kernel/PythonTools.m`. Pure-WL QASM emitter in
`Kernel/QuantumCircuitOperator/Properties.m:420-422`.

### 1.1 The object and its generic escape hatch

- `QuantumCircuitOperatorToQiskit[qco]` / `qco["Qiskit"]` (`Qiskit.m:63-143`)
  builds a Python `QuantumCircuit` and returns `QiskitCircuit[bytes]` (a pickled
  blob). `shortcutToGate` (`:11-61`) is the gate-translation map; unknown gates
  fall back to a `UnitaryGate` at **single precision** (`ComplexReal32`, `:57`,`:59`).
- `QiskitCircuit` is a compound object with a **generic method dispatcher**:
  `qc["Eval", attr, args, kwargs]` and `qc["EvalBytes", ...]` (`:149-153`) call
  *any* method on the pickled circuit by name, via `PythonEvalAttr`
  (`PythonTools.m:60-83`). Concrete accessors built on it: `["Qubits"]`,
  `["Depth"]`, `["Ops"]` (= `count_ops`), `["Clbits"]`, `["Diagram"]`,
  `["QuantumCircuitOperator"]` (roundtrip back to QF), `["Matrix"]`, `["Decompose", n]`.

`["Eval"]` is already a version-resilient primitive **for instance methods**.
But `transpile` is a *module-level* function `qiskit.transpile(circuit, **kwargs)`,
not a method on the circuit, so `["Eval"]` cannot reach it. That gap is the whole
reason a dedicated transpile method is needed.

### 1.2 Gap A: `qiskitInitBackend` injects `save_statevector` (the clbit-less trap)

`qiskitInitBackend` (`Qiskit.m:314`, options `:312`) builds an IBM / AWS-Braket /
Aer backend. For the default Aer path, when the circuit has no clbits it injects
`qc.save_statevector()` (or `save_density_matrix()`) at `:404-410`. That Aer
instruction cannot be translated to a discrete native basis, so any downstream
`transpile(qc, backend)` against that circuit fails. It also has **no way to accept
a coupling map or a qiskit `Target`.** Every method that routes through
`qiskitInitBackend` inherits this trap.

### 1.3 Gap B: `["QASM3"]` re-transpiles and de-natures

`qiskitQASM` (`Qiskit.m:530-550`), exposed as `["QASM"|"QASM2"]` (`:552`) and
`["QASM3"]` (`:554`), **always** runs `transpile(qc, backend)` first (`:533-534`)
against the backend from `qiskitInitBackend` (the default AerSimulator). On an
already-native circuit this **re-fuses** the native `sx`/`rz` back into `u2`/`u3`.

This is not hypothetical. Live, on a circuit already transpiled to `{x,sx,rz,cz}`
(real `count_ops = {rz, sx, cz, measure}`):

```
transpiled["Ops"]    ->  <|rz -> 5, sx -> 3, measure -> 2, cz -> 1|>   (correct)
transpiled["QASM3"]  ->  ... u2(0, -pi) q[0]; u2(pi/2, -pi) q[1]; ...   (de-natured!)
```

A direct `qasm3.dumps(circuit)` on the same bytes gives the correct native gates.
So `["QASM3"]` actively destroys the native-basis property the user transpiled for.

### 1.4 Gap C: `["Transpile"]` is option-starved and backend-coupled

`qiskitTranspile` (`Qiskit.m:586-608`), exposed as `["Transpile", basisGates, opts]`
(`:610`), accepts only `basis_gates` (positional), `"InitialLayout"`, and
`"OptimizationLevel"`. It exposes **none** of: `coupling_map`, `target`,
`layout_method`, `routing_method`, `translation_method`, `scheduling_method`,
`approximation_degree`, `seed_transpiler`, `unitary_synthesis_method`, ... (qiskit
2.4.1 `transpile` has **24** parameters; this exposes 3). It also passes
`backend=backend` from `qiskitInitBackend` (`:600`), so it inherits the
save_statevector trap **and** the backend's target can override/conflict with the
requested `basis_gates`. Net effect for the native-hardware path: it cannot route
(no coupling map) and often fails outright on a clbit-less circuit.

### 1.5 Gap D: the pure-WL emitter is not native

`qco["QASM"]` (`Properties.m:420-422`) is a pure-WL OpenQASM-3.0 emitter that
riffles each operator's `["QASM"]`. It emits generic `U(theta,phi,lambda)` and
`ctrl @` modifiers, not a hardware native basis, and is independent of qiskit. It
is fine as a generic serialization but is not a transpile path.

### 1.6 The known API trap (do not try to "fix" the parse)

`qco["Qiskit", "Transpile", ...]` (comma form) is parsed as one property `"Qiskit"`
with extra args; it matches no rule, fires `QuantumCircuitOperator::undefprop`, and
returns the circuit unchanged. The chained form `qco["Qiskit"]["Transpile", ...]` is
the correct surface and is what this design extends. (Documenting the trap is worth
a one-line note; changing property-arg parsing is out of scope and risky.)

---

## 2. Qiskit API audit (live, qiskit 2.4.1)

Read live via the QF Python bridge (`inspect.signature`), not hardcoded:

```
qiskit.__version__ = 2.4.1

transpile parameters (24):
  circuits, backend, basis_gates, coupling_map, initial_layout, layout_method,
  routing_method, translation_method, scheduling_method, dt, approximation_degree,
  seed_transpiler, optimization_level, callback, output_name,
  unitary_synthesis_method, unitary_synthesis_plugin_config, target, hls_config,
  init_method, optimization_method, ignore_backend_supplied_default_methods,
  num_processes, qubits_initially_zero

Target.from_configuration parameters (8):
  basis_gates, num_qubits, coupling_map, instruction_durations,
  concurrent_measurements, dt, timing_constraints, custom_name_mapping

generate_preset_pass_manager parameters (19):
  optimization_level, backend, target, basis_gates, coupling_map, initial_layout,
  layout_method, routing_method, translation_method, scheduling_method,
  approximation_degree, seed_transpiler, unitary_synthesis_method,
  unitary_synthesis_plugin_config, hls_config, init_method, optimization_method,
  dt, qubits_initially_zero, _skip_target
```

Key observations that shape the design:

- **`transpile` is a flat kwargs surface.** Almost every parameter is a scalar /
  string / list that marshals directly through `ExternalEvaluate` (WL Association ->
  Python dict, confirmed). Only a few carry rich Python objects: `coupling_map`
  (-> `CouplingMap`), `target` (-> `Target`), `backend`, `initial_layout`,
  `hls_config`, `unitary_synthesis_plugin_config`. This is exactly the situation a
  "pass-through kwargs + coerce a known short list" design is built for.
- **`transpile`, `generate_preset_pass_manager`, and `Target.from_configuration`
  share almost the same vocabulary.** A single option model covers all three.
- **A `Target` bundles basis + directed connectivity + per-instruction
  errors/durations.** Passing a directed coupling map (or a Target built from one)
  makes qiskit emit two-qubit gates in the calibrated direction automatically. This
  is verified below and matches how a real backend (`arbel`'s `transpiler_target`)
  describes itself.

### 2.1 Robustness reality (qiskit decides, not WL) — verified live

| Input | Result | Surfaced as |
|---|---|---|
| `basis_gates = {x, sx, rz, cz}` + directed `coupling_map`, `opt=3` | native, `$0`/`$1`, `cz $1, $0` directed | `QiskitCircuit` |
| non-universal basis `{x, cz}` | `TranspilerError: Unable to translate ...` | `Failure["QiskitTranspile", ...]` |
| `Target.from_configuration` with only `{x,sx,rz,cz}` (no `measure`) | `HighLevelSynthesis is unable to synthesize "measure"` | `Failure["QiskitTranspile", ...]` |
| unknown option key `no_such_param` | caught **before** the call | `Failure["QiskitTranspileOption", ...]` |

The third row is a genuine design lesson: a `Target` built from
`from_configuration` must declare `measure`/`reset` (via `custom_name_mapping` or by
fetching a real device target, which already carries them) or measurements fail at
HighLevelSynthesis. The function must surface all of these as a `Failure`, never a
fake QASM string.

The canonical verified recipe (the thing the kernel should do):

```
take qc["Qiskit"]["Bytes"]                            # clean circuit, no save_statevector
transpile(circuit, basis_gates=..., coupling_map=CouplingMap(directed pairs),
          optimization_level=...)                     # DIRECT call, no backend
qasm3.dumps(out)                                       # faithful, native
```

Live output of that recipe on a QF Bell with measurements:

```
OPENQASM 3.0;
include "stdgates.inc";
bit[2] c;
rz(pi/2) $0; sx $0; rz(pi/2) $0;
rz(pi/2) $1; sx $1; rz(pi) $1;
cz $1, $0;
sx $1; rz(pi/2) $1;
c[0] = measure $0;
c[1] = measure $1;
```

Native, physical-qubit, directed. This compiled on `arbel` in the PLAN
investigation.

---

## 3. The recommended design

### 3.1 Version-resilience mechanism (the centerpiece)

**Decision: pass-through `**kwargs` + live signature validation + a small coercion
layer.** Justification, weighing the two candidates the task names:

- *(A) Forward an `Association` straight as `**kwargs`.* This is the load-bearing
  choice. The WL bridge already renders an `Association` as a Python dict usable for
  `**` (`PythonEvalAttr` does exactly this at `PythonTools.m:71-72`). A new qiskit
  release that adds a parameter (e.g. a future `placement_method`) works with **zero
  QF change**: the user just puts the new key in the Association. A removed/renamed
  parameter needs no QF change either; it simply becomes an "unknown option" caught
  by (B). This is strictly more future-proof than any WL-side `Options[...]` list,
  which would have to be hand-edited every qiskit release (the same manual-lockstep
  problem the TN version pin already suffers, per `qf-tn-integration.md`).

- *(B) Runtime-introspect `inspect.signature(transpile)`.* Use it as a **guard**, not
  as the transport. Before the call, compute `set(opts) - set(signature.parameters)`.
  If non-empty, return a clean `Failure` that names the bad keys **and** lists the
  accepted set together with `qiskit.__version__`. This converts an opaque Python
  `TypeError: transpile() got an unexpected keyword argument` into an actionable QF
  message, and it is the mechanism that keeps the surface honest across versions. It
  can also drive FE autocompletion (populate special-arg completion from the live
  signature instead of a static list).

- *(C) Algorithmic CamelCase mapping (the chosen option-key style).* The WL surface
  uses CamelCase option names; an algorithmic `CamelCase -> snake_case` transform
  (`"BasisGates" -> "basis_gates"`, `"SeedTranspiler" -> "seed_transpiler"`)
  produces the qiskit key, and validation (B) runs on the translated keys. This is
  the crucial detail that keeps the CamelCase choice from defeating resilience: it is
  **not** a hand-maintained per-param dictionary, so a future qiskit parameter
  (`PlacementMethod` -> `placement_method`) auto-works with no QF change. All 24
  current `transpile` params round-trip cleanly. Raw snake_case keys also pass
  through untouched, so qiskit's own names remain usable as an escape hatch (verified
  by the `snake-case-passthrough` test).

- *Coercion layer.* The only WL-side knowledge that is qiskit-shaped is the short
  list of rich-typed keys. The prototype coerces:
  - `coupling_map`: a list of integer pairs -> `CouplingMap([tuple(e) for e in ...])`
    (pass a `CouplingMap` through untouched).
  - `target`: a `QiskitTarget` byte-blob -> `pickle.loads`; a pre-built `Target`
    passes through.
  Everything else passes verbatim. This list is short, stable, and additive: if a
  future key needs coercion, add one `if`. It does **not** grow with the qiskit option
  count.

Why not "introspect and rebuild a typed WL options interface"? Because that
re-introduces the hardcoded list under a different name and re-creates the
lockstep-maintenance burden. The pass-through is the point.

### 3.2 Surface and home

All methods live on `QiskitCircuit` (in `Qiskit.m`), the bytes-level object where
composition already happens. `qco["Qiskit"]` remains the gateway from a
`QuantumCircuitOperator`.

```wolfram
qc["Qiskit"]["Transpile", optsAssoc]        (* -> QiskitCircuit | Failure *)
qc["Qiskit"]["Transpile", k1 -> v1, ...]    (* flat-rule sugar, same thing *)

tcirc = qc["Qiskit"]["Transpile", <|
    "BasisGates" -> {"x", "sx", "rz", "cz"},
    "CouplingMap" -> {{1,0},{1,3},{2,0},{2,3}},
    "OptimizationLevel" -> 3 |>];

tcirc["Ops"]                                  (* count_ops, reject deep results *)
tcirc["Depth"]
tcirc["QASM3"]                                (* faithful native dump *)
tcirc["QuantumCircuitOperator"]               (* roundtrip back to QF *)

target = QiskitTarget[<| "BasisGates" -> {...}, "NumQubits" -> 21,
                         "CouplingMap" -> directedPairs |>];
qc["Qiskit"]["Transpile", <| "Target" -> target, "OptimizationLevel" -> 3 |>];
```

- **Option keys are CamelCase WL-style** (`"BasisGates"`, `"CouplingMap"`,
  `"OptimizationLevel"`, `"SeedTranspiler"`, ...), converted to qiskit snake_case by
  the algorithmic `CamelCase -> snake_case` transform of §3.1(C). Because the mapping
  is algorithmic (not a curated per-param dictionary), the surface needs no
  per-release maintenance and any future qiskit parameter auto-maps. Raw snake_case
  keys also pass through untouched, so qiskit's own names stay usable as an escape
  hatch (verified by the `snake-case-passthrough` test). All 24 current `transpile`
  params round-trip cleanly through the transform.

### 3.3 Compositional vs monolithic — decision

**Compositional.** `["Transpile"]` returns a `QiskitCircuit`. Rationale:

- It composes with the entire existing `QiskitCircuit` surface for free, including
  the exact thing the task asks for: `["Ops"]` (count_ops) and `["Depth"]` so callers
  can reject results that exploded the two-qubit-gate count.
- It mirrors the existing `["Decompose", n] -> QiskitCircuit` and the (current)
  `["Transpile"] -> QiskitCircuit` shape; minimal surprise, no new mental model.
- The end-to-end "QF -> native QASM3" is then a one-liner by composition:
  `qc["Qiskit"]["Transpile", opts]["QASM3"]`. A monolithic helper is unnecessary; if
  desired it is trivial sugar on top (`qco["NativeQASM", opts] :=
  qco["Qiskit"]["Transpile", opts]["QASM3"]`), but it should not be the primitive.

### 3.4 Faithful QASM export — decision

**Make `QiskitCircuit["QASM3"]` / `["QASM2"]` dump the circuit as-is** (call
`qasm3.dumps` / `qasm2.dumps` directly on the pickled bytes), with an output guard:
the result must be a string starting with `OPENQASM`, else a
`Failure["QiskitQASMExport", ...]`. Re-engage the current transpile-against-backend
path **only when a `"Backend"` or `"Provider"` option is explicitly supplied** (this
is also what the AWS-Braket conversion branch in `qiskitQASM:536-540` needs). Rule:

- no `Backend`/`Provider` -> faithful dump of exactly the circuit you hold.
- `Backend`/`Provider` given -> transpile against it, then dump (today's behavior,
  preserved for the IBM/Braket flows).

This makes the common case correct (an already-native circuit stays native) while
keeping the hardware-provider workflows intact. It is a behavior change to the
default and therefore a documented bugfix (see §4).

### 3.5 `QiskitTarget` — the device model

Add a `QiskitTarget[spec_Association]` constructor returning a `QiskitTargetObject`
(pickled `Target`), parallel to `QiskitCircuit`. It forwards `spec` to
`Target.from_configuration` with the same validate-against-live-signature +
coerce-`coupling_map` mechanism. Why a first-class object:

- A `Target` is **reusable across many circuits** (fetch once per device per session;
  `arbel`'s coupling map is a device property, not per-circuit).
- It is the **version-stable, highest-value** interface: basis + directed
  connectivity + per-instruction errors/durations in one object, enabling
  error-aware layout/routing and correct directed two-qubit gates.
- It maps onto how real backends describe themselves. A device's
  `transpiler_target` (e.g. fetched read-only from IQCC) is exactly a
  `from_configuration`-shaped spec.

`["Transpile"]` accepts the target three ways, all coerced uniformly: `"target" ->
QiskitTarget[...]` (a `QiskitTargetObject`), `"target" -> <|spec|>` (built inline),
or the decomposed `"basis_gates"` + `"coupling_map"` (no Target object). Device
constants stay out of QF entirely (goal 4): the user supplies the spec.

### 3.6 Error handling

Following the `Qiskit.m` idiom (`Enclose @ Block[...]`, `Confirm`):

- **`Enclose`/`Confirm` -> `Failure`** at the WL boundary; `ConfirmBy` asserts the
  Python return is `_QiskitCircuit | _Failure` (or `_String | _Failure` for dumps).
- **Precondition message** (unknown option): `Failure["QiskitTranspileOption",
  <|"MessageTemplate" -> "Unknown transpile option(s): \`1\`. Accepted (qiskit \`2\`):
  \`3\`.", "MessageParameters" -> {unknown, version, accepted}, "Unknown" -> ...|>]`.
  Accepted set + version read live.
- **Transpile failure** (non-universal basis, coupling too small, target missing
  `measure`, exotic-gate blowup that raises): catch the Python exception, return
  `Failure["QiskitTranspile", <|..., "PythonError" -> str(e)|>]`. **Never a fake
  circuit.**
- **QASM output guard**: non-`OPENQASM` result or a dump exception ->
  `Failure["QiskitQASMExport", ...]`.
- **Depth/count surfacing**: the returned `QiskitCircuit` already answers `["Ops"]`
  and `["Depth"]`; callers reject deep results themselves. (Optional sugar:
  `["TranspileInfo"]` returning `<|"Ops" -> ..., "Depth" -> ..., "Qubits" -> ...|>`.)

### 3.7 Reuse of existing infrastructure

- `QuantumCircuitOperatorToQiskit` / `qco["Qiskit"]["Bytes"]` (`:63-143`, `:147`) is
  the clean-circuit source — unchanged. (Caveat from `mistakes.md`: unknown gates
  fall back to single-precision `UnitaryGate`; acceptable for hardware-bound work but
  documented.)
- `PythonEvaluate` (`PythonTools.m:43-57`) is the transport. **Idiom to respect**
  (learned the hard way during prototyping): the `<* $x *>` template resolves the
  symbol name *literally* in the passed context, so these helpers must scope their
  interpolated globals with **`Block`, not `Module`** (`Module` lexically renames
  `$x` -> `$x$nnn` and the template then finds an unassigned global). Every existing
  python-interpolation helper in `Qiskit.m` already uses `Block` for this reason.
- `["Eval"]` / `PythonEvalAttr` (`:149-153`, `PythonTools.m:60-83`) stays the generic
  escape hatch for any *instance method* not given a named wrapper.
- Existing `QiskitCircuit["Ops"|"Depth"|"QuantumCircuitOperator"|"Diagram"]`
  (`:199-292`) compose on the transpiled result with no change.

---

## 4. Migration / compatibility

### `["Transpile"]` — fix in place (capability superset)

The new option-Association form is a strict superset. Keep the positional
`basisGates` first argument and the `"InitialLayout"` / `"OptimizationLevel"` /
`"CouplingMap"` CamelCase options working via the alias map, so existing call sites
are unaffected. The one behavioral change is that it **no longer routes through
`qiskitInitBackend`** (no backend built, no `save_statevector` injected) unless a
`"Backend"`/`"Provider"` is explicitly passed. Since the old backend-coupled path
*failed* on clbit-less circuits and silently constrained the basis to the
AerSimulator target, this is a bugfix, not a regression. Risk: a caller who relied on
`["Transpile", basis, "Backend" -> ...]` building an Aer target — preserved, because
passing `"Backend"`/`"Provider"` still engages that path.

### `["QASM3"]` / `["QASM2"]` — fix in place (documented behavior change)

Default becomes faithful (dump-as-is). This is **technically breaking but a bugfix**:
the old default de-natured already-native circuits (§1.3, verified). The old
transpile-first behavior is reachable two ways: explicitly pass
`"Backend"`/`"Provider"` (re-engages transpile-against-backend, and is required for
the AWS-Braket conversion branch anyway), or compose `["Transpile", ...]["QASM3"]`.
Ship with:

- `VerificationTest`s asserting (a) faithful dump of an un-transpiled circuit is
  generic (`h`/`cx`), (b) faithful dump of a `{x,sx,rz,cz}`-transpiled circuit is
  native and directed, (c) `["QASM3", "Backend" -> Automatic]` still transpiles.
- A CHANGELOG note: "`QiskitCircuit["QASM3"]`/`["QASM2"]` now export the circuit
  faithfully instead of transpiling against the default AerSimulator first; pass a
  `Backend`/`Provider` option for the previous transpile-then-dump behavior."

### Leave alone

- The comma-form parse trap `qc["Qiskit", "Transpile", ...]` — document it, do not
  change property-arg parsing.
- The pure-WL `qco["QASM"]` emitter (`Properties.m:420`) — orthogonal; keep as the
  generic, qiskit-free serialization.

---

## 5. Alternatives considered (and why rejected)

1. **A static WL `Options[QiskitTranspile]` mirroring the 24 params.** Rejected:
   re-creates the manual per-release lockstep (the exact pain the TN version pin
   shows) and silently breaks when qiskit evolves — the explicit anti-goal.
2. **Reuse `["Eval"]` for transpile.** Rejected: `["Eval"]` does `getattr` on the
   *circuit instance*; `transpile` is a module-level function. `["Eval"]` stays the
   tool for instance methods.
3. **Monolithic `qco -> QASM3 string` function.** Rejected as the primitive (loses
   composition, hides `count_ops`/`depth`); fine as one-line sugar on top.
4. **Build a `PassManager` from `generate_preset_pass_manager` instead of calling
   `transpile`.** Deferred. `transpile` is the stable, documented, highest-level
   entry; its kwargs are a near-superset of `generate_preset_pass_manager`'s. A
   `QiskitPassManager` object (same kwargs-forwarding philosophy) is a clean *future*
   extension for callers who want to run a custom pass pipeline or reuse one pass
   manager across circuits, but it is not needed for the native-hardware path.
5. **Vendor/clone qiskit into the repo to read pass internals.** Not needed: live
   `inspect.signature` is the authoritative, version-correct source for the installed
   qiskit, and it is exactly what the resilience mechanism consumes at runtime.

---

## 6. Prototype (validated)

`prototype.wl` implements the mechanism as standalone Global-context functions
(`qTranspile`, `qDumpQASM`, `QiskitTarget`) **without editing the kernel**; they map
1:1 onto the proposed `QiskitCircuit` methods. `validate.wl` exercises them:

```
PASS transpile-returns-QiskitCircuit      (CamelCase opts Association -> composable QiskitCircuit)
PASS native-count-ops                     (count_ops = {cz, measure, rz, sx})
PASS faithful-native-qasm3                (directed "cz $1, $0", no u2/u3)
PASS faithful-generic-qasm3               (un-transpiled dump stays h/cx)
PASS snake-case-passthrough               (raw qiskit names also forwarded; resilience escape hatch)
PASS unknown-option-failure               (bad key caught pre-call, named in Failure)
PASS non-universal-basis-failure          ({x,cz} -> Failure, not fake QASM)
PASS target-builds                        (QiskitTarget via from_configuration)
PASS end-to-end-oneliner                  (qc -> Transpile -> QASM3 native string)
9 / 9 passed
```

The prototype uses the chosen CamelCase surface (`"BasisGates"`, `"CouplingMap"`,
`"OptimizationLevel"`) via the algorithmic transform, and the `snake-case-passthrough`
test confirms raw qiskit names still work.

Run: `wolframscript -file OngoingProjects/QiskitTranspile/validate.wl` (loads the
paclet from the repo, boots qiskit 2.4.1 via the QF Python session).

### Proposed kernel placement when greenlit

- `Qiskit.m`: add `PackageExport["QiskitTarget"]`; `qiskitTranspileKwargs[qc, opts]`
  (the `Block`-scoped pass-through+coerce+validate helper), wired as
  `qc_QiskitCircuit["Transpile", opts:OptionsPattern[]|_Association]`; rewrite
  `qiskitQASM` to dump faithfully unless `Backend`/`Provider` is set;
  `QiskitTarget[spec_Association]` + `QiskitTargetObject` formatting.
- `Tests/QuantumCircuitOperator.wlt` (or a new `Tests/Qiskit.wlt`): the nine
  assertions above, plus the `["QASM3"]` behavior-change tests from §4.
- `CHANGELOG`: the `["QASM3"]`/`["QASM2"]` faithful-default note.

---

## 7. Decisions

Resolved 2026-06-04:

1. **Option-key style — RESOLVED: full CamelCase WL-style.** Implemented via the
   algorithmic `CamelCase -> snake_case` transform (§3.1(C), §3.2) so resilience is
   preserved (no per-param dictionary; future qiskit params auto-map; raw snake_case
   still passes through). Prototype + `snake-case-passthrough` test updated.
2. **`["QASM3"]` default flip — RESOLVED: flip to faithful-dump default.** Old
   transpile-first reachable via an explicit `Backend`/`Provider` option. Ship with
   the §4 `VerificationTest`s + CHANGELOG note.
3. **Proceed — RESOLVED: design-only for now.** No kernel implementation this
   session.

Still open (not blocking, raise when implementation is greenlit):

4. **Scope of `QiskitTarget` at implementation time.** First-class `QiskitTarget`
   object in v1 (recommended; prototyped), or inline `"Target" -> <|spec|>` coercion
   only to start?
5. **`QiskitPassManager` (alt. 4).** In scope later, or explicitly dropped?
6. **Placement of this deliverable.** Currently `OngoingProjects/QiskitTranspile/`
   (tracked). Move to a `Documentation/` ROADMAP-style folder (like the Stabilizer
   subproject) if you prefer.
7. **Greenlight to implement.** When you want it, I will implement §6's kernel
   changes + tests + CHANGELOG on a branch off `main` (CamelCase surface + faithful
   `["QASM3"]` default, per the resolved decisions above).
