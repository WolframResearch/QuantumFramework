# Idioms & Performance


When QF offers multiple paths to the same answer, this records which is correct, which is fast, and which is wrong but tempting.

## Stabilizer simulation: complexity map (`Kernel/Stabilizer/`, updated through `02cc92d9`)

For Clifford / stabilizer circuits, the `PauliStabilizer` tableau path is exponentially cheaper than the dense state vector, **but only as long as you stay in the tableau**. The moment you materialize, you pay O(2ⁿ) or O(4ⁿ). **Perf note (June 2026 rework):** the hot paths now run on a packed machine-word representation (`Packed.m`, 62 tableau bits/word) and, for whole circuits, a `Compile`-to-C bulk fold (`Compiled.m`). Asymptotics below are unchanged; the constants dropped by 1-3 orders of magnitude, putting QF at 2.5-5.3× of Stim and ahead of QuantumClifford.jl at $n=1000$ (full benchmarks: `OngoingProjects/Platform Comparison/QF-Stabilizer-Optimization-Report.md`).

| Operation | Cost | Where | Note |
|---|---|---|---|
| **Whole Clifford circuit** `ps["ApplyCircuit", specs]` | **~0.2-0.3 µs/gate**, flat in n | `Compiled.m` | the route to quote: one compiled-C fold over the encoded circuit; reached automatically by `qc[Method -> "Stabilizer"]` / `PauliStabilizerApply` when every gate is an encodable Clifford and signs are concrete; returns a canonical (unpacked) object |
| Clifford gate update (H/S/CNOT/…), single | O(n) per gate, packed words | `GateUpdates.m` → `Packed.m` | direct tableau mutation; drops global phase by design; returns a **packed** object (compare via `["Tableau"]`, not `SameQ`) |
| `PauliStabilizer[n]` register ctor | closed-form (was O(n²·2n) pattern-match) | `Constructors.m` | direct `Signs`/`Tableau` assembly; n=1000 355 → 4.9 ms |
| Composition `left[right]` | O(n³) | `Compose.m` | symplectic mod-2 matmul + `phaseLookup`; still O(n³) but packed reconstruction made it ~55× faster (n=100 2577 → 46 ms); a crash on all-zero/≤1-generator rows was fixed |
| Single-qubit / Pauli-string measurement | O(n²) (packed; was rank-3 rebuild per call) | `Measurement.m`, `PauliMeasure.m` | AG destabilizer trick on packed words; single-qubit n=500 16.5 → 0.19 ms, non-det Pauli-string n=500 42.3 s → 59 ms |
| Bipartite entropy `["Entropy", A]` | **O(n³), no OOM** | `Entropy.m` | Fattal closed form `rank_{F₂} - |A|`; independent of `|A|` |
| Inner product, `Method -> "ClosedForm"` | O(n³) | `InnerProduct.m` | **magnitude + ±1 sign only**, concrete signs required |
| Inner product, `Method -> "Direct"` (default) | **O(2ⁿ)** | `InnerProduct.m` | materializes vectors; full complex phase; n ≤ ~8 |
| `ps["State"]` / `["QuantumState"]` | **O(4ⁿ)** | `Conversions.m` | materializes 2ⁿ vector; n ≤ ~10 |
| `PauliStabilizer[qs]` (tomography) | **O(4ⁿ)** | `Constructors.m` | every Pauli expectation; n ≤ ~8 |
| `ps["Circuit"]` (AG synthesis) | O(n³) gates | `Conversions.m` | NOT length-minimized (Reid24/Winderl23 are TODO) |
| `qmo[ps]` / `qc[ps]` hybrid interop | O(n²) fast path, else O(2ⁿ) | `HybridInterop.m` | stays in tableau for Pauli-basis QMO / named Pauli channel; else `::nonpaulibasis` materialize |

**Idiom**: keep stabilizer work tableau-resident. For long gate streams use **`ps["ApplyCircuit", specs]`** (or `PauliStabilizerApply`), not a one-gate-at-a-time fold of objects (the per-gate route still pays object wrap/unwrap; ~66× over the original but ~25× the compiled route) and not a `QuantumCircuitOperator` built from 10⁴ gate objects (host-framework construction is the bottleneck then). Use `["Entropy", A]` (not the Schmidt path) and `Method -> "ClosedForm"` inner products for large n; only call `["State"]`/`["QuantumState"]`/`["Operator"]` at the very end, for n small enough to materialize. The non-Clifford P/T boundary returns a `StabilizerFrame` whose term count **doubles per non-Clifford gate**, fine for a few T gates, exponential for T-rich circuits. Its dense readout (`["StateVector"]`/`["State"]`) is **phase-coherent** as of `61bc0e39`: a **gate-built** frame carries a relating Pauli per component, so materialization builds the reference component once and relates the rest by those Paulis, so a single-state `["StateVector"]` or `⟨O⟩` is amplitude-exact up to one global phase (same contract as a bare `PauliStabilizer`), at a cost of ≈ +15% `ByteCount` / +8-12% build time per frame. Only the inner-product **phase** between two distinct, separately-gauged frames stays gauge-limited (magnitude correct); hand-built / `Plus`-combined frames keep the prior independent per-component readout. **Storage-form trap:** single gates return packed objects, `ApplyCircuit` returns canonical ones; both are valid everywhere, but `SameQ` across the two forms fails, so compare `["Tableau"]` / `["Stabilizers"]`. **CompiledFunction trap:** the compiled fold silently falls back to interpreted (~100× slower, no warning) if any argument arrives unpacked, so the kernel forces `Developer`ToPackedArray` at every boundary.

## KAK decomposition: numeric vs symbolic (`QuantumOperator/Properties.m`)

`qo["KAK"]` (and `["Decompose"]`) has two internal paths, auto-selected by whether the matrix is numeric:
- **Numeric** (`TwoQubitKAK`): magic-basis SVD/eigensolve, fast, emits `U`-gate + CNOT gadgets, recovers the global phase by trace-projection. **Output equals the operator up to a global phase only** (see `mistakes.md`) — compare up-to-phase.
- **Structured-symbolic** (`TwoQubitKAKSymbolic`): exact arithmetic under a real-parameter assumption, **time-bounded at `$kakSymbolicTimeBudget = 60` s**, with a probe-substitution self-check that returns the operator *undecomposed* (`::kaksymbolic`) rather than a silently-wrong circuit. A generic dense symbolic 2-qubit unitary is a Root blowup and will hit the budget — pass a *structured/parametric* gate (e.g. RXX(θ), a controlled-phase), not a fully generic symbolic matrix.

## QASM: native WL vs qiskit path (`QuantumQASM.m`, `OpenQASMImport.m`)

`QuantumQASM` import and the **default** export are **pure Wolfram Language — no Python session, no qiskit dependency**, immune to qiskit version churn (depend only on the OpenQASM spec). Reach for the qiskit path (extra argument: a native gate set, `QiskitTarget`/`CouplingMap`, a transpiler option, or a `"Provider"`/`"Backend"`) only when you need transpilation to a device basis. The qiskit transpile core validates options against the **live** `inspect.signature(transpile)` (version-resilient) rather than a hard-coded dictionary. Qubit-only (non-2-dim wire → `Failure`).

The native emitter's coverage and fidelity were extended after `ac32bda8` (`541a104b`, `d7b8b3fd`, `a44c1eea`):
- **Full-precision reals**: `qasmReal` (`QuantumOperator/Properties.m:952`) emits angles via `ToString[N[x], InputForm]` with `*^`→`e`, replacing `NumberForm` (which truncated to 6 digits and rendered tiny/huge magnitudes as a line-wrapping `×10` superscript). Round-trip fidelity is no longer capped at 6 digits.
- **Single-qubit + controlled gates carry the right phase**: `qasmEmitSimple` (`:960`) is guarded by `qo["SquareQ"]` (so a non-square state-injection tensor that also reports `Dimensions === {2,2}`, like a `Cup`, no longer reaches `UnitaryAngles`) and uses `UnitaryAnglesWithPhase`. `UnitaryAngles` alone folds `Arg[mat[[1,1]]]` into the φ/λ angles (a Z-conjugation, not a global phase), giving the wrong operator when `mat[[1,1]]` is complex. The controlled path emits the leftover phase as a controlled `gphase`.
- **State injection** (`QuantumCircuitOperator/Properties.m:427-461`): a no-input/≥1-output operator (a `|+⟩` ket or `Cup`) is lowered to `reset` + a preparation gate sequence rather than failing. The circuit's qubit count is now `qco["Max"]` (highest wire), not `qco["Arity"]`, so a circuit that *starts* with state prep (input arity 0) still declares all its wires.
- **Non-computational measurement rebasing** (`:464-475`): an X/Y (or any non-Z) measurement basis `B` is lowered to `Inverse[B]` (decomposed via `QuantumShortcut`) followed by a computational measurement, the same trick the qiskit path uses. A computational measurement is left byte-for-byte unchanged; a genuine POVM (non-square basis) is left as-is. This was a silent-drop bug before.

### Error-aware transpilation (working tree, UNCOMMITTED at `f9dc1cdc`)

In progress in `Qiskit.m` / `IBMQuantum.m` / `QuantumQASM.m`; verified by source read only, not committed, not benchmarked. The qiskit-bridge transpile is becoming error-weighted:
- An `"OptimizationLevel"` option (Automatic → 1) on `IBMJobSubmit`, `qiskitQASM`, and `qiskitPrimitiveSubmit`; transpilation runs against the backend's own qiskit `Target` (per-instruction error + duration, coupling map), so layout/routing are error-aware. Users may raise to 2/3.
- `QiskitTarget["FromBackend" -> name, "Provider" -> p]` reads any qiskit-ecosystem backend's populated `Target` into a QF `QiskitTarget` (one implementation through the shared `qiskitInitBackend`: IBM, AWS Braket, and a local `GenericBackendV2` fake device that needs no credentials/network), filtering `operation_names` to real basis gates and carrying `InstructionProperties` (error/duration).
- `QuantumQASM`'s `build_target` applies the per-instruction properties via `update_instruction_properties` (stable across qiskit 1.x/2.x), recording rather than silently dropping rows for unknown gates/loci.
This stays bridge-based (no native SABRE router). Re-confirm against the working tree before relying on it; it may land or change.

## Quiet-tightening (no-Quiet hygiene, `ac32bda8`; cache `Quiet` removed `f9dc1cdc`)

Several broad `Quiet[...]` / `ConfirmQuiet` calls were narrowed to tag-specific suppression: `Utilities.m` (`Eigensystem::eivec0`, `Inverse::sing`, `FindIntegerNullVector`), `QuantumOptimization.m` (`ExternalEvaluate::*`), `QuEST.m` (`Needs::*`), `QuantumState`/`QuantumMeasurement` Formatting (`ConfirmQuiet`→`Confirm`). **The last broad cache `Quiet` is now gone:** the property-cache `Set`-on-non-symbol (`Quiet[HeadProp[...] = result, Rule::rhs]`, with its `(* TODO *)`) was replaced at `f9dc1cdc` by the `cacheProperty` store, which only commits a pattern-free key and so never needs to suppress `Rule::rhs` (see the Caching section). Net: fewer silently-swallowed messages; if a previously-`Quiet`ed warning now surfaces, it is the kernel telling you something real.

## Caching (two-tier)

QF has **TWO independent global caches**:

| Cache | Where | What it caches | Bypass |
|---|---|---|---|
| `$QuantumFrameworkPropCache` | `QuantumFramework.m:15` | Property lookups on `QuantumState`/`QuantumOperator`/`QuantumBasis`/`QuditBasis`. | Set to `False`. State cache also bypassed when `ParameterArity > 0`. |
| `$QuditBasisCache` | `QuditBasis/NamedBases.m:27` | Construction of named QuditBases. So `QuditBasis["GellMann"[d]]` is nearly free after first call. | Set `$QuditBasisCaching = False`. Built-in bypass for `RandomMIC \| RandomHaarMIC \| RandomBlochMIC`. |

- `$QuantumFrameworkPropCache = False` does NOT clear `$QuditBasisCache`. Long sessions should periodically clear: ``Wolfram`QuantumFramework`PackageScope`$QuditBasisCache = <||>``.
- **`Memoize[f]`** (`Utilities.m`) is the QF-internal memoization for arbitrary functions. **Rewritten `5b696969`:** results now live in a private `Association` keyed on `Hash[Hold[args]]` (held key kept for collision safety), not appended `f[args] -> res` DownValues, so a patterned argument can no longer trip `Rule::rhs`. New `Memoize[f, "CacheQ" -> pred]` option restricts which calls cache (pred runs once on a miss, never on a hit). Clear by `Clear[f]`. Used by `Init.m` on `FromOperatorShorthand`.
- **No automatic invalidation on object mutation.** Building `qs1`, querying `qs1["VonNeumannEntropy"]`, then re-using the same expression after substitution may return a stale cached value. Pin intermediates as named symbols to force re-eval.

### Property-cache machinery repair (`dfc741dc`/`f9dc1cdc`, 2026-06-14): the per-hit tax is gone

Before this, every property access (even a cache **hit**) paid a fixed tax: the cache-eligibility guard evaluated `obj["Basis"]["ParameterArity"]`, and `obj["Basis"]` had no direct rule, so it fell through ~124 `QuantumOperatorProp` clauses into the generic `AllProperties`-intersection delegation, ~0.4 ms, never cached (`"Basis"` is excluded by design). A warm hit on the trivial getter `op["Order"]` cost ~400 µs; the apply path made ~66k such accesses per circuit application. Three coordinated fixes removed it:

1. **Direct `"Basis"` getter** on every wrapped head (the basis lives in the wrapped state): `QuantumOperator`, `QuantumChannel`, `QuantumHamiltonianOperator`, `QuantumMeasurementOperator`, `QuantumMeasurement` (`QuantumState` already had one). So the guard no longer triggers the fall-through.
2. **Widened cache-exclusion lists** to the $O(1)$ structural getters: `QuantumOperator` now also excludes `Order`/`InputOrder`/`OutputOrder`; `QuantumState` excludes `State`; each wrapper head adds `Basis`; `QuditBasis` stops caching `Representations`/`Properties`. Excluded props short-circuit the eligibility guard, so those reads get cheaper too (caching a structural extract was wasted memory anyway).
3. **`cacheProperty` store** (`Utilities.m`, `HoldFirst`) replaces `Quiet[HeadProp[obj,prop,args] = result, Rule::rhs]`: it commits the cache DownValue only when the held key is `FreeQ` of pattern constructs, so the latent over-matching DownValue (which the `Quiet` was masking) is impossible. Patterned/meta keys are left uncached.

**Effect (M2 Pro): warm 105-gate 12-qubit apply 836 → 273 ms, cold 2956 → 783 ms, `op["Basis"]`×105 32 → 2 ms, output bit-identical.** Practical consequence: the old advice "warm up, then set `$QuantumFrameworkPropCache = False` for repeated numeric apply" is now nearly moot (default-on 273 ms is within ~1.2× of the cache-off floor 222 ms). The two-tier cache table above is unchanged in structure; only the per-hit machinery got cheap. Full decomposition: `OngoingProjects/Platform Comparison/QF-Apply-Path-Deep-Audit.md` (status banner) + master comparison §12.2.

## Profiling

- `$QuantumFrameworkProfile = True` enables `EchoTiming` profile output via the internal `profile[label]` wrapper (`Utilities.m:362`). Use to find slow paths in your own QF computations.

## Sparse vs dense

- **`pauliMatrix`, `spinMatrix`, `identityMatrix`, `fanoMatrix` return `SparseArray`** (`Utilities.m:120, 201–239`). Algebraic chains of Paulis preserve sparsity.
- **`MatrixPartialTrace` is dense-only** (`Utilities.m:132`). Uses `ArrayReshape` + `TensorContract`. For large mixed states, this is a memory bottleneck.
- **`SparseArrayFlatten` only kicks in when `Length[dims] > 11`** (`Utilities.m:275`). For typical few-qubit tensors, falls back to `Flatten` (which densifies).

## Symbolic vs numeric

- **`SymbolicQ[expr]`** (`Utilities.m:65`) flags any free symbol that's not `Rational | Complex | NumericFunction | Constant`. Used internally to gate symbolic vs numeric paths.
- **`eigensystem` calls `Simplify`** on output (`Utilities.m:151`). For symbolic matrices with parameters, `Simplify` is the dominant cost.
- For pure-numeric eigendecomposition where speed matters, prefer the underlying `Eigensystem[matrix, Method -> "Hermitian"]` directly — QF's wrapper is generic.

## Eigendecomposition

- QF wraps `Eigensystem` generically — **no Hermitian-specialized path** even for density matrices and Hamiltonians.
- Options on `eigensystem`: `"Sort" -> True | False | <comparator>`, `"Normalize" -> True | False`, `Chop -> True | False`.
- `eigenvalues` and `eigenvectors` (`Utilities.m:185, 189`) are thin wrappers — same generic path.

## Hard dimension limits (entropy / purity)

`VonNeumannEntropy` and `Purity` for `QuantumState` have a **hardcoded upper bound** of `Dimension < 2^10 = 1024` (`QuantumState/Properties.m:455, 475`):

```
ConfirmAssert[qs["Dimension"] < 2^10]
```

Above 10 qubits, the assertion fails inside an `Enclose` and the property returns `Indeterminate` **silently** (no warning). For ≥10-qubit work:
- Reduce to a tractable bipartition first (`QuantumPartialTrace[qs, complement]`), then compute entropy on the reduction.
- Use `QuantumMPS` for tensor-network routing.
- Use `qs["TraceNorm"]` (`Properties.m:375`) — singular-value sum, no dimension limit.

Display limits (cosmetic, not computational):
- `Dimension < 2^9 = 512`: full `ComplexArrayPlot` icon in summary.
- `Dimension >= 2^9`: sparse-array placeholder icon.
- Purity/Entropy in summary: `TimeConstrained[..., 1]` (1 sec); `$Aborted` if it doesn't finish.

Bloch-related properties only work for `Dimension == 2` (`Properties.m:929, 962, 987`). For general dimensions, use `qs["BlochVector"]` which computes `Tr[ρ G_i]` against `GellMannMatrices[d]`.

## Property dispatch & lazy materialization

`QuantumState` property dispatch (`Properties.m:51-62`) caches **conditionally**:
- `$QuantumFrameworkPropCache == True` (default).
- AND the prop is not `"Properties"`, `"AllProperties"`, or `"Basis"`.
- AND the basis has no free parameters (`ParameterArity == 0`).

### `"Properties"` is now O(1); the union moved to `"AllProperties"` (post-`ac32bda8`)

On every wrapped head (`QuantumState`, `QuantumOperator`, `QuantumChannel`,
`QuantumHamiltonianOperator`, `QuantumMeasurementOperator`, `QuantumMeasurement`,
`QuantumCircuitOperator`), `obj["Properties"]` now returns just the head's stored static list
`$Quantum…Properties`, a constant-time lookup, cheap to render in summaries. The old behavior
(union with the wrapped object's properties, which forces materializing the wrapped object's full
catalog) is now `obj["AllProperties"]`, and it is on each head's cache-prevention list (recomputed
every call). For "does this respond to X?" introspection use `"AllProperties"`; for a quick cheap
header list use `"Properties"`. See `architecture.md` §2 and the `mistakes.md` HEAD-resync entry.

Parametric states bypass cache. Substitute parameters before heavy property extraction:

```mathematica
qs  = QuantumState[{a, b}, basis];   (* parametric, never caches *)
qs2 = qs[{a -> 1, b -> 0}];          (* concrete, caches *)
qs2["VonNeumannEntropy"]             (* cached on second call *)
```

Lazy materialization costs to know:
| Property | Cost | Why |
|---|---|---|
| `qs["DensityMatrix"]` on Vector-type state | N → N² memory | `KroneckerProduct[v, Conjugate[v]]` (`Properties.m:415`). |
| `qs["StateVector"]` on Matrix-type pure state | O(d³) time | Eigendecomposes density matrix (`Properties.m:130-145`). |
| `qs["Type"]` | O(d³) | `PositiveSemidefiniteMatrixQ[N @ DensityMatrix]` (`:481-497`). |
| `qs["SchmidtBasis"]` | O(d³) | `SingularValueDecomposition` (`:524-548`). |
| `qs["VonNeumannEntropy"]` | O(d³) + Simplify | Eigenvalues + Total[-p Log[p]] + Simplify (`:449-460`). |

Pin the cheap representation early (`vec = qs["StateVector"]` once) to avoid repeated lazy work.

## Validity caching (symbol-level)

`QuantumState[state, basis]` validates once via `quantumStateQ` and then marks itself with `System\`Private\`HoldSetValid` (`QuantumState/QuantumState.m:21`). Subsequent calls to `QuantumStateQ` short-circuit to `HoldValidQ`. **You don't pay the validation cost twice for the same expression.**

Implication: Returning the same `QuantumState[...]` expression structure is faster than rebuilding from raw matrices. Hold and reuse.

## Eigenvalues vs Eigensystem inconsistency

`qs["Eigenvalues"]` (`Properties.m:438`) calls bare `Eigenvalues[qs["DensityMatrix"]]` — **NOT** the QF-wrapped `eigenvalues` helper. So no `Chop`, `Sort`, or `Normalize` options.

`qs["Eigenvectors"]` (`:440`) and `qs["Eigensystem"]` (`:442`) DO use the wrapped versions and accept options.

For consistent behavior across all spectral queries, prefer `qs["Eigensystem", opts]` and take `First` for eigenvalues.

## Inversion

- **`MatrixInverse` silently falls back to `PseudoInverse`** if `Inverse` fails (`Utilities.m:303`). No warning emitted.
- For non-invertible inputs, you get a pseudo-inverse without knowing. Caller responsibility to check `MatrixRank` first if it matters.
- Change-of-basis (`QuantumState/QuantumState.m:145-177`) uses `MatrixInverse` on `Output["ReducedMatrix"]` and `Input["ReducedMatrix"]` — non-invertible bases produce silent pseudo-inverse results.

## Operator ordering (hidden subtlety)

- `autoOrderQ[order]` (`Utilities.m:84`) accepts a single order, `Automatic`, OR a `{out, in}` 2-tuple.
- **Operators can have asymmetric input/output orders.** `op[qs]` may rearrange qubit ordering. Verify with `op["Order"]` if surprised.

## Property dispatch — two forms

```mathematica
qs["Purity"]               (* string form *)
qs["VonNeumannEntropy", 2]  (* {name, args...} form *)
```
Both valid (`Utilities.m:76`).

## State representation duality

- `stateQ` (`Utilities.m:80`) accepts BOTH a vector (pure state) AND a square matrix (density matrix).
- A `QuantumState[vector]` and `QuantumState[densityMatrix]` are both valid; they're different *purities* of the same conceptual object.

## Tensor flatten threshold

- Below 11 dims: `SparseArrayFlatten` falls through to plain `Flatten`. **Densifies sparse tensors** for typical few-qubit work.
- Workaround for many-mode bosonic/fermionic problems where dim > 11: stay in `SparseArray` form explicitly; avoid intermediate `Flatten`.

## Auto-installed dependencies

- First load of `Wolfram/QuantumFramework` triggers PacletInstall of `IBMQuantumPlatform 0.0.2` (from local asset) and `Wolfram/TensorNetworks 1.0.5` (from cloud; pin bumped from 1.0.4 in `468f2dac`) if either is missing (`QuantumFrameworkMain.m:6-11`).
- This is silent and one-shot; subsequent loads skip it. Be aware of network dependency on first run for `TensorNetworks`.

## Anti-patterns *(seeded; expand from object Properties.m)*

- ❌ Repeatedly building `QuantumState` from the same matrix in a loop without caching. Use `Module[{qs = QuantumState[...]}, ...]`.
- ❌ Calling `QuantumPartialTrace` on a sparse N-qubit density matrix expecting it to stay sparse — it densifies.
- ❌ Using `Inverse[op["Matrix"]]` and assuming it gives a true inverse — `MatrixInverse` falls back to pseudo-inverse silently.
- ❌ Constructing a `QuantumCircuitOperator` and immediately calling `["Matrix"]` for an N-qubit circuit with N large — that materializes the 2^N × 2^N gate matrix. Compose, then evaluate at the end (or use `QuantumMPS` routing).
- ❌ Calling `qc["Matrix"]` and then doing `qs["DensityMatrix"]` operations downstream when an idiomatic property like `qc[qs]["Purity"]` would route through cached state form.
- ❌ Using `Do`/`For`/`AppendTo` to build `QuantumCircuitOperator` lists. Use `Table` / `Map` and pass the list once. (User preference per `wl-native-style` skill.)

## QuantumOperator-specific notes

### Composition has a direct fast path (since `a321a6f2`, 2026-05-08)

`qo1[qo2]` now dispatches to a **direct `TensorContract` fast path** that bypasses `QuantumCircuitOperator` and `TensorNetwork` whenever the operators are index-compatible (`QuantumOperator/QuantumOperator.m:344-410`). Two specialized rules:

1. **Matrix-form / channel path** (`MatrixQ` on either side, same Picture, compatible orders): `qo1["Bend", shift][qo2["Bend", shift]]["Unbend"]` — doubles the qudits, applies the vector rule, then unbends back. Recovers ~30× over the legacy route for cases like `CNOTmix @ CNOTmix`.
2. **Vector-form path** (both `VectorQ`, same Picture, index compatibility on `OutputOrder` / `InputOrder`): joins the two `StateTensor`s with `TensorContract` on the shared-qudit axis pairs, then permutes to canonical `[outputs … inputs]` order. ~2× for small aligned vector ops.

Both paths bail out to the **legacy fallback** `QuantumOperator @ QuantumCircuitOperator[{qo2, qo1}]` (`:412`) when the index pattern doesn't match (different Picture, disjoint orders that overlap weirdly, etc.). Same routing remains for `qo[qs]`, `qo[qmo]`, `qo[qc]`, `qo[qm]` (no fast path there).

**Practical implication:** the old advice "always batch into `QuantumCircuitOperator[{...}]` because `op1[op2[op3[...]]]` traverses N circuits" is **partly obsolete**. For nested `qo[qo]` composition with aligned indices, the nested form is now fine (and sometimes faster than building a 1-step circuit). The circuit form is still preferred when:
- You're mixing `QuantumOperator`s with `QuantumState`s / `QuantumMeasurement`s — those still route through `QuantumCircuitOperator`.
- You want lazy evaluation / circuit-level optimization (TN routing, Schrödinger FullSimplify, Qiskit export).
- Your operators have non-aligned orders that would trigger the fallback anyway.

`QuantumState[qs1][qs2]` got the **same treatment** (`QuantumState/QuantumState.m:302-389`): a Kraus-style fast path for `Schrodinger`-picture vector states with prefix-matching dimensions, plus a `Double`-rep fallback for matrix-form `qs1`. `phiBra[rho2]`: 150 ms → 1.7 ms.

### Symbolic vs numeric: `HermitianQ` ≠ `UnitaryQ` symmetry

| Property | Implementation | Symbolic-safe? |
|---|---|---|
| `qo["HermitianQ"]` | `HermitianMatrixQ[qo["Matrix"]]` (`Properties.m:443`) | ✅ yes — no `N` |
| `qo["UnitaryQ"]` | `UnitaryMatrixQ[N[qo["Matrix"]]]` (`:445`) | ❌ no — applies `N` first |

For symbolic operators (e.g. `R_y(θ)`), `UnitaryQ` can return `False` even when the operator is genuinely unitary, because `N` doesn't simplify enough. Verify symbolically: `Simplify[qo["Matrix"].qo["Dagger"]["Matrix"] == IdentityMatrix[d]]`.

### `qo` spectral properties are consistent (unlike `qs`)

For QuantumOperator, all three spectral keys use the wrapped helpers from `Utilities.m` with `"Sort" -> False, "Normalize" -> True` defaults:

```
qo["Eigenvalues"]  = eigenvalues[qo["MatrixRepresentation"], opts, "Sort" -> False, "Normalize" -> True]
qo["Eigenvectors"] = eigenvectors[..., same defaults]
qo["Eigensystem"]  = eigensystem[..., same defaults]
```

(`Properties.m:447-451`). All accept option overrides. **Different from QuantumState**, where `qs["Eigenvalues"]` calls bare `Eigenvalues`.

### `PauliDecompose` two-path performance

- Qubit operators (`Dimensions == {2..}`): tensor-reshape recursion, `O(4^n)` (`Properties.m:626-637`). Fast for n ≤ ~10 qubits.
- Non-qubit qudits: `SolveValues` over `Tuples[{"I","X","Y","Z"}, n]` Pauli operator basis (`:639-649`). Slow even for moderate dim.

For qudits, decompose against Gell-Mann directly: `Tr[op["Matrix"] . #] & /@ GellMannMatrices[d]`.

### Order vs Label semantics

Crucial: control structure (CNOT, Toffoli, etc.) is encoded in the Label of the operator, NOT in the Order. The Order specifies *which qudits the operator touches*; the Label specifies *which of those are control vs target*.

```
qo["Order"]        (* {{outputQudits}, {inputQudits}} *)
qo["ControlOrder1"]  (* extracted by pattern-match on Label: Subscript["C",_][c1,_] *)
qo["TargetOrder"]   (* InputOrder minus ControlOrder *)
```

Reordering a CNOT does NOT change which qudit is control. To change control, build a new operator with a different label or use the named-operator constructor.

### `Diagonalize` collapses the Order

`qo["Diagonalize"]` (`Properties.m:475-482`) replaces the operator with its eigendecomposition AND reduces the Order to a single qudit (`Take[#, UpTo[1]] & /@ qo["Order"]`). The eigenstates become the new basis. Useful for measurement-operator construction; **destructive of multi-qudit structure**.

### `MatrixRepresentation` vs `Matrix`

- `qo["Matrix"]` = `qo["Sort"]["StateMatrix"]` (`Properties.m:180`) — sorted by Order, in the operator's stored basis.
- `qo["MatrixRepresentation"]` = `qo["Computational"]["Matrix"]` (`:182`) — in the computational basis.

For algebraic comparisons across operators in different bases, use `MatrixRepresentation`. For "what does this matrix look like in MY basis", use `Matrix`.

## Circuit application & compilation

### Method dispatch (default = TensorNetwork)

`qco[qs]` calls `quantumCircuitApply` (`QuantumCircuitOperator/QuantumCircuitOperator.m:87-102`). The `Method` option selects the engine:

| Method | Path | Notes |
|---|---|---|
| `Automatic` | `TensorNetwork` | **Default.** |
| `"TensorNetwork"` or `{"TensorNetwork", subOpts}` | `TensorNetworkApply` | Uses the Wolfram/TensorNetworks paclet. `subOpts` are `FilterRules`-forwarded to `TensorNetworkApply`/`TensorNetworkCompile`; this is the **only** channel for TN options now (e.g. `{"TensorNetwork", "Path" -> path}` picks the contraction order). The old top-level option leak was removed at `77b9ccba`. |
| `"Schrodinger"` | `Fold` of gates with `FullSimplify` per gate | Slow for deep symbolic circuits. Has `EvaluateWithProgress` live progress bar. Since `f9dc1cdc` each gate is applied with `Method -> {"TensorNetwork", "Computational" -> False}` (`QuantumCircuitOperator.m:99`), the post-`77b9ccba` sub-option contract; the prior bare `"Computational" -> False` was dropped + emitted `OptionValue::nodef` per gate. Output bit-identical. Benefited from the cache refactor: n=12 105-gate fold 14.3 → 6.6 s. |
| `"QuEST"` | `QuESTApply` via QuEST kernel | **Qubits only** (Dimensions == {2..}, `QuEST.m:58`). |
| `"Qiskit"` or `{"Qiskit", subOpts}` | Qiskit roundtrip | External bridge. |
| `"Stabilizer"` | `PauliStabilizerApply` | For Clifford/stabilizer circuits. |

The same `Method` sub-option threading applies to `quantumCircuitCompile`
(`QuantumCircuitOperator/Properties.m:148-162`) for `qco["QuantumOperator", Method -> …]` /
`qco["Compile", …]`. `Options[quantumCircuitApply]`/`Options[quantumCircuitCompile]` are now just
`{Method -> Automatic}`; they no longer Join the TN option set, so a bare TN option at top level is
dropped. New contraction-path option: `"Path"` on `TensorNetworkCompile` (`TensorNetwork.m:182`),
forwarded as `TensorNetworkContract[net, OptionValue["Path"], …]` (`:215`).

### `TensorNetworkCompile` auto-conversions

The compile path makes two automatic conversions you should know about:

1. **Auto-bend** (`TensorNetwork.m:203-211`): if any operator is in matrix form OR the circuit has trace qudits AND we're not in PhaseSpace picture, the circuit is replaced with `circuit["Bend"]`. This doubles the operator count via Stinespring dilation. Hidden cost.

2. **Computational-basis conversion** (`TensorNetwork.m:30-35, :100-105`): with `"Computational" -> True` (default), non-computational input/output bases get extra `MatrixInverse[ReducedMatrix]` and `ReducedMatrix` tensors prepended/appended. Extra nodes in the network.

### Circuit cache exclusions

Circuit properties that always bypass the cache (`$QuantumCircuitPreventCache`, `Properties.m:30-33`):

`Association`, `Elements`, `Operators`, `FullOperators`, `FullElements`, `Options`, `Diagram`, `Icon`, `Qiskit`, `QiskitCircuit`, `QuantumOperator`, `Flatten`, `Double`, `Bend`, `DiscardExtraQudits`, `ToggleExpand`, `ExpandElements`, `Picture`, `Parameters`, `ParameterArity`, `AllProperties`.

For tight loops over these, pin manually with `Module[{ops = qco["Operators"]}, ...]`.

### Composition shortcuts

`Composition[qco1, qco2, ...]` and `RightComposition[qco1, qco2, ...]` are upvalues (`QuantumCircuitOperator.m:151-159`) that **flatten** the elements. So `qco1 @* qco2` is structurally equivalent to a single circuit with all gates concatenated. Cleaner than nested circuit-of-circuits.

### `Layers` for parallelism

`qco["Layers"]` (`Properties.m:188-203`) partitions operators into layers where each layer has no qudit overlap. Useful as a model of circuit depth and for parallel hardware mapping. `qco["Depth"]` is `Length @ qco["Layers"]`.

### Trotterization gate-count growth

`QuantumCircuitOperator[{"Trotterization", ops, order, reps}]` (`NamedCircuits.m:631-644`) uses Suzuki recursion for even orders ≥ 4 (`:606-614`):
- order 1: linear in `ops × reps` (Lie-Trotter).
- order 2: 2× linear (Strang).
- order 4: ~5× order-2 ≈ 10× linear.
- order 6: ~25× linear.
- order 8: ~125× linear.

For accuracy, prefer increasing `reps` over increasing `order` past 2 — same convergence rate, no exponential gate blowup.

### Qiskit roundtrip precision

`QuantumCircuitOperatorToQiskit` falls back to `Unitary` gate with `NumericArray[Normal[N[matrix]], "ComplexReal64"]` for unknown shortcuts (`Qiskit.m:57, :59, :105`). **Double precision (64-bit complex)** — `N[matrix]` is `MachinePrecision`, so arbitrary-precision input is collapsed to ~15-16 decimal digits (not the ~7 of single precision).

Stay in the WL kernel (`Method -> "Schrodinger"` or `"TensorNetwork"`) only when you need precision *above* machine double. For ordinary numeric work the roundtrip preserves full double precision. *(Corrected 2026-06-05: was `"ComplexReal32"`/single precision; live kernel at HEAD `fb16ffc2` uses `"ComplexReal64"`.)*

### `QiskitCircuit["QuantumOperator"]` materializes 2^N × 2^N

`Qiskit.m:518-521` calls Python `Operator(qc).data` — full unitary. Same memory limits as Qiskit's own statevector simulator.

For large Qiskit circuits, prefer `qiskitCircuit["QuantumCircuitOperator"]` (gate-list form) and apply via the QF tensor-network method.

### Circuit drawing: dynamic mode is interactive

`CircuitDraw[qco]` (`Draw.m:962-987`) returns a static `Graphics`. `CircuitDraw[qco, "Dynamic" -> True]` (`:962-976`) returns a `DynamicModule[Dynamic[Graphics[...], TrackedSymbols :> {qc}]]`. **Sub-circuits become click-to-expand**: each is wrapped in `EventHandler[..., {"MouseClicked" :> (qc = qc["ToggleExpand", path])}]` (`:824-843`).

For static publication-quality output, default `CircuitDraw[qco]` is correct. For exploration of nested circuits, dynamic mode is far more useful.

### Wire thickness encodes dimension

`Draw.m:53` — `defaultWireThickness[dim] := Directive[AbsoluteThickness[Log2[dim/4 + 3/2]], Opacity[Min[0.3 Log2[dim], 1]]]`. So qubit (d=2) wires are thin, qutrit thicker, etc. To override, pass `"DimensionWires" -> False` to `CircuitDraw` — flat thickness 1.

### `Expand` recursion levels

`"Expand"` option for `CircuitDraw` controls nested sub-circuit unfolding (`Draw.m:743-755`):
- `True` → `Infinity` (fully expand all nested)
- `False` → `-Infinity` (collapse all nested to single gates)
- `None` → 0 (default-level)
- Integer → exact recursion depth

## Hamiltonian evolution paths (numeric vs symbolic)

`QuantumHamiltonianOperator` does not load on `main` (see `function-index.md` / `activity-status.md`), so do Hamiltonian evolution on `QuantumOperator`; both paths live there:

| Property | Solver | When to use |
|---|---|---|
| `qo["NEvolutionOperator"]` (`QuantumOperator/Properties.m:586-618`) | `NDSolveValue` | Numeric H, single parameter. **Default for performance.** |
| `qo["EvolutionOperator"]` (`QuantumOperator/Properties.m:584` → `QuantumEvolve[qo, None]`) | symbolic `DSolveValue` | Closed-form symbolic evolution. **Fragile for complex H** — same fragility the (non-loading) `qho["EvolutionOperator"]` had. |

For a parametric `H(θ)` where you want `U(t)` numerically: build a `QuantumOperator[H, "ParameterSpec" -> {t, t0, t1}]` and call `["NEvolutionOperator"]`. For symbolic insight on a small example: `qo["EvolutionOperator"]` and inspect the `DSolveValue` output. Don't use the symbolic path for production numeric work — `DSolveValue` will hang.

For `MatrixExp[-I H t]` directly: this is the cheapest path when H has a known closed-form matrix exponential.

## Partial trace path selection

`QuantumPartialTrace` has **two paths with different memory profiles** (`QuantumPartialTrace.m`):

| Form | Path | Cost |
|---|---|---|
| `QuantumPartialTrace[qs, {q, q', ...}]` (flat list) | calls `MatrixPartialTrace[qs["DensityMatrix"], qudits, qs["Dimensions"]]` (`:30`) | **Dense** — materializes density matrix even for pure states. |
| `QuantumPartialTrace[qs, {{out, in}, ...}]` (pair list) | calls `TensorContract[qs["StateTensor"], ...]` directly (`:36-46`) | Works on the state tensor; cheaper for vector-form pure states. |
| `QuantumPartialTrace[qco, ...]` | builds a Cup-Cap circuit (`:72-87`) | Categorical implementation — evaluate via tensor-network. |

For pure-state partial trace where in/out indices are symmetric, prefer `QuantumPartialTrace[qs, {{q, q}, ...}]` to avoid materializing the density matrix.

## Concurrence: vector vs matrix asymmetry

`qs["Concurrence"]` (or `QuantumEntanglementMonotone[qs, partition, "Concurrence"]`) has two paths in `QuantumEntanglement.m:45-51`:
- **Vector (pure state)**: closed form `Sqrt[Max[0, 2(1 - PartialTrace^2.Norm)]]` — fast, O(d²).
- **Matrix (mixed state)**: full `ConcurrenceVector` (`:27-39`) — iterates over `O(d₁²·d₂²)` Y-operator pairs, each requiring `SingularValueList` of a square root. **O(d⁴)** for d-dim subsystem.

For mixed states with d > 4, this is the dominant cost in any entanglement-quantification workflow. Prefer `Negativity` (2-qubit) or `Realignment` (general) when the state is mixed and Concurrence isn't strictly needed.

## QuantumStateEstimate cost profile

`QuantumStateEstimate` (`QuantumStateEstimation.m`) is **not a single computation** — it's a 4-stage pipeline:
1. Linear inversion via `PseudoInverse` (cheap).
2. Physical projection (eigendecomposition).
3. Maximum likelihood: **300 iterations** by default (`:236`), each involving a likelihood evaluation over all POVM elements.
4. Bayesian sampling: 150-step burn-in + 150-step adaptive (`:294-296`).

For quick state inversion only, call `inversion[result]` (`:67-82`) directly. For full Bayesian estimate, expect multi-second to multi-minute runtime depending on dim and POVM count.

## `MergeInterpolatingFunctions` for chained evolution

When you `QuantumEvolve[H, ...]` and the resulting state has `InterpolatingFunction` entries, then evolve again, QF's `MergeInterpolatingFunctions` (`QuantumEvolution.m:213-249`) composes the InterpolatingFunctions into a single linearly-interpolated `SparseInterpolationFunction`. This lets you chain numeric evolutions without re-discretizing.

Activate via the `"MergeInterpolatingFunctions" -> True` option (default).

## Adiabatic gap analysis (QHO non-loading)

The dedicated `qho["Gap"]` / `qho["MinimumGap"]` properties exist only in the non-loading `QuantumHamiltonianOperator` source on `main`, so they are not usable. Do it user-side on a parametric `QuantumOperator`:

```mathematica
gap[v_] := Differences[Sort[qo["Eigenvalues"][[1 ;; 2]]] /. param -> v][[1]]

(* single point *)              gap[0.3]
(* explicit n-point sweep *)    gap /@ Range[t0, tf, (tf - t0)/n]
(* minimum *)                   FindMinimum[gap[param], qo["ParameterSpec"][[1]]]
```

Each evaluation of `gap[v]` computes the two smallest eigenvalues of `H(v)` — `O((2^N)^3)` for an N-qubit Hamiltonian. The old `qho["Gap"]` defaulted to 100 points, which was a known perf trap; the user-side form forces explicit choice of resolution.

## Tensor network thresholds *(further details in TensorNetwork.m read)*

- Above how many qubits should `QuantumMPS` routing replace full-matrix simulation?
- For `QuantumCircuitOperator` with N gates and Q qubits, what's the bond-dimension cost vs. dense?

## Channel form trade-offs *(pending — read `Kernel/QuantumChannel/QuantumChannel.m` and `Properties.m`)*

- Which form does QF use internally for `channel[state]`? Kraus, Choi, Stinespring, natural?
- Cost asymmetry of converting between forms.

## VQE / QAOA cost-function patterns (from `QuantumOptimization.nb`)

The canonical symbolic cost-function pattern in QF:

```
ansatz = QuantumCircuitOperator[{...gates...}, "Parameters" -> {θ1, θ2, ...}];
|φ⟩ = ansatz[];
⟨φ| = |φ⟩^†;
cost[θ1_, θ2_, ...] = ("Scalar" //. ⟨φ|@H@|φ⟩) // ComplexExpand // FullSimplify;
```

This is faster than building the cost in a `Module`: QF's parameter machinery reuses the symbolic `QuantumState` and only substitutes parameter values at call time. The `["Scalar"]` extraction is constant-time after symbolic construction; `FullSimplify` once is much cheaper than `Simplify` per call.

For optimization: `MinValue[cost[θ1,θ2,...], {θ1,θ2,...}]` gives exact ground-state energy for analytically-tractable Hamiltonians (single-qubit Pauli sums). For numeric, `NMinimize` typically beats `FindMinimum` on these surfaces; `FindMinimum` often diverges (e.g., `QuantumOptimization.nb` VQLS section returns `θ ~ 10^10` from `FindMinimum` while `NMinimize` converges).

For *parameterized* states, prefer `vqs[<|θ1->v1, θ2->v2|>]["Formula"]` over closing over Module-bound symbols — the framework propagates parameters through the entire pipeline.

## QAOA layer scaling

For Max-Cut and similar combinatorial problems, single-layer QAOA gives **approximate** solutions only. `QuantumOptimization.nb` 4-edge-graph example shows `θ_opt = ±π/4` gives `C = 3` but exact answer is `4`. Two layers (`(qcost@mixer) /* (qcost@mixer)`) is enough to reach the optimum on small graphs.

Pattern: `QuantumCircuitOperator[{stateprep, qcost[γ1], mixer[θ1] -> qubits, qcost[γ2], mixer[θ2] -> qubits}]` — note that `qcost[γ1]` is shorthand for setting the `γ` parameter via the Parameters substitution.

## QuantumLinearSolve recovery: NMinimize vs gradient descent

`QuantumLinearSolve[m, b]` does internal hybrid optimization. If you build VQLS manually, two pitfalls:
1. `FindMinimum` typically fails (cost-function landscape has wide barren plateaus). Use `NMinimize`.
2. Even after optimization, `|ψ(ω_opt)⟩ = e^{iφ}|b⟩` only up to a global phase. To recover amplitudes, compute `globalphase = optψ / (b/Norm[b])` (componentwise; ratios should agree), then `x_sol = (1/Mean[globalphase]) * (Norm[b]/Norm[m]) * x(ω_opt)`.

## Quantum Natural Gradient Descent vs vanilla GD

Use QNGD over GD when the cost-function landscape has near-singular regions (`fubini["MatrixForm"]` near-degenerate). QNGD computes `θ → θ − η · F^+(θ) · ∇cost(θ)` where `F` is the Fubini-Study metric tensor (use the **variational state** `|ψ⟩`, not the ansatz, as input — the variational state captures the full geometry).

QNGD on the same problem can converge in ~2 steps where vanilla GD needs 50+ (`QuantumOptimization.nb` final HHL/VQLS example).

## Phase-space dimension parity rule

`QuantumWignerTransform` and the related state properties have a **parity-dependent dimension**:
- Odd `d`: `qs["PhaseSpace"]` is `d × d`. `qs["QuasiProbability"] == qs["PhaseSpace"]` for single qudit.
- Even `d`: `qs["PhaseSpace"]` is `2d × 2d`. The first quadrant is the actual Wigner state vector; the other 3 are sign-flipped copies. Take every 2nd element of column/row sums to recover position/momentum probabilities.
- Multi-qudit: `qs["PhaseSpace"]` is a tensor of dim list, with each pair `{2 d_i, 2 d_i}` (even) or `{d_i, d_i}` (odd). `qs["QuasiProbability"]` is the *flattened total-dim* matrix.

Marginalizing the joint quasi-probability matrix recovers `Probabilities` in computational basis (column sum) and momentum basis (row sum).

## Beam-splitter ordering and method

For `BeamSplitterOperator[{θ, φ}, ...]`:
- Default `Method -> "MatrixExp"` is faster for **numeric** angles (compiled matrix exponential).
- `Method -> "Recurrence"` uses an analytic recurrence; faster for **exact** angles (e.g., `π/4`, `π/2`) where it avoids floating-point.
- Pure symbolic with `Method -> "Recurrence"` returns clean trig expressions.

Same pattern for `DisplacementOperator` and `SqueezeOperator`'s `"Ordering"` option: `"Normal"` for large `|α|, |ξ|` and small truncation; `"Weak"` for small `|α|` and high accuracy.

## ServiceConnect for quantum hardware

`ServiceConnect["IBMQ"]` (or `"AWS"`) returns a `ServiceObject` with `["Requests"]` listing methods (`"Backend"`, `"BackendQueue"`, `"RunCircuit"`, `"JobResults"`, `"JobStatus"`). Two patterns to run a circuit:

1. **Synchronous (kernel blocked):** `qc[Method -> {"Qiskit", "Provider" -> "IBMProvider", "Backend" -> ibm_brisbane}]`. Returns `QuantumMeasurement` after the IBM job completes.
2. **Asynchronous via LocalSubmit:** wrap in `LocalSubmit[Needs[...]; qc[...], HandlerFunctions -> <|"TaskFinished" -> ...|>]`. Returns `TaskObject[]` immediately; the kernel is free.

For >8 qubits using Qiskit method, the result is an `Association` of bitstring counts, not a `QuantumMeasurement` (`QuantumCircuitOperator.nb:1620`).

## EntanglementLayer strategies

| Strategy | Connections | Use case |
|---|---|---|
| `Automatic` | All-pairs | Maximally expressive, deepest |
| `"Linear"` | `(1,2), (2,3), ...` | NISQ-friendly, lowest depth |
| `"ReverseLinear"` | `(n-1, n), (n-2, n-1), ...` | Backward chain |
| `"Pairwise"` | even-odd pair layers | Two-layer alternating |
| `"Circular"` | All-pairs + wrap | Periodic boundary, ring topology |

`"Linear"` is the canonical hardware-efficient choice; `"Pairwise"` is preferred for short-depth Trotterization-like circuits.

## Doc reference pages are stubs — use Usage.m + tutorials

Most of `Documentation/.../Symbols/*.nb` reference pages have `XXXX` placeholders for "Tech Notes", "Related Demonstrations", "Related Links", etc. Sections that actually contain content: signature templates, named-symbol lists, and "Examples" / "More Examples". The "Properties & Relations", "Possible Issues", "Neat Examples" sections are usually empty.

Don't trust a missing example to mean a feature is unimplemented. Cross-reference:
1. `Kernel/Usage.m` — signatures (the contract).
2. `Documentation/.../Tutorials/*.nb` — usage (contains workflows, not just examples).
3. `Kernel/<file>.m` — implementation.

## Notebook-confirmed idioms (Phase 4: 90+ user notebooks read)

### Composition operators

The canonical "compose two circuits" pattern is the right-composition operator `/*`:

```
qc = qc1 /* qc2 /* qc3            (* read left-to-right *)
qc = QuantumMeasurementOperator[{1,2}] @* qc1   (* same as qc1 then measure *)
```

`/*` flattens into a single `QuantumCircuitOperator`. `@*` (built-in `Composition`) is right-to-left like usual `f @* g`. Pre-composing a measurement with `@*` is equivalent to appending it to the circuit.

### Expectation value: 4 idiomatic forms

From `Notebooks/Examples Landing Page/Expectation-value.nb`:

```
(* 1. Density-matrix trace — works for any state *)
Tr[op["Matrix"] . qs["DensityMatrix"]]

(* 2. Bend + Cap circuit — categorical form *)
QuantumCircuitOperator[{qs["Bend"], op,
  Sequence @@ Table["Cap" -> {i, i + n}, {i, n}]}][]["Scalar"]

(* 3. Operator-on-state via partial trace *)
QuantumPartialTrace[op[qs["Operator"]]]["Scalar"]

(* 4. Measurement Mean *)
QuantumMeasurementOperator[op][qs]["Mean"]
```

All four return the same scalar. Path 4 is the most idiomatic for QF; path 1 is the textbook formula; path 2 is useful when you want to keep things as a circuit (e.g., for visualization).

### Inner product / fidelity via Dagger composition

```
(qc /* {ψ^†})[]["Norm"]^2          (* fidelity |⟨ψ|ψf⟩|² *)
QuantumPartialTrace[(qc /* {ψ^†})[], {1,2,3}]["Scalar"]   (* analytical *)
```

`QuantumState["..."]["Dagger"]` (or `^†`) at the end of a circuit collapses all qudits to a scalar. `(qc /* {ψ^†})[]` returns a `QuantumState` with no qudits whose `["Scalar"]` is the inner product, and `["Number"]` returns the underlying complex number.

### Coherent matrix-on-vector via Multiplexer

Pattern from `Examples Landing Page/Coherent_applying_of_a_matrix_on_a_vector.nb`:

```
σ = QuantumOperator[A]["PauliDecompose"];
mp = QuantumCircuitOperator["Multiplexer" @@ Keys[σ]];
ancillary = QuantumState[Sqrt @ Values[σ], "Label" -> "Sqrt(c)"];
qc = QuantumCircuitOperator[{
  ancillary, QuantumState[x] -> Range[L+1, L+1+n], mp,
  (ancillary["Conjugate"])^†
}];
qc[]["AmplitudesList"]   (* === A.x *)
```

This is the building block for `QuantumLinearSolve` — implements `A.x` coherently using Pauli decomposition + multiplexer + ancilla amplitudes.

### Lindblad evolution: 4 equivalent forms

From `Examples Landing Page/LindbladEq.nb`:

```
(* Form 1: jump operators + rates separate *)
QuantumEvolve[H, {L1, L2}, {γ1, γ2}, ρ0, {t, 0, T}]

(* Form 2: rates folded into operators *)
QuantumEvolve[H, {Sqrt[γ1] L1, Sqrt[γ2] L2}, ρ0, {t, 0, T}]

(* Form 3: Liouvillian explicit *)
ℒ = QuantumOperator["Liouvillian"[H, {Ls}, {γs}]];
expℒ = QuantumEvolve[I ℒ, None, {t, 0, T}];
ρt = expℒ[t][ρ0]

(* Form 4: time-independent shortcut (when ℒ is t-independent) *)
ρt = Exp[ℒ t][ρ0]
```

All four forms rebase the Hamiltonian and every jump operator into the Hamiltonian's
output basis (computational when `H === None`) since 2026-06-11, so named jump operators
(`"J-"`, `"J+"`, ...) keep their operator meaning and every route returns its state in the
same basis. See `mistakes.md` "Named-basis jump operators".

### Symbolic Bell-correlation pattern

From `Bellstheorem.nb` and `Examples Landing Page/CHSH_game.nb`:

```
(* Compute Bell expectation as a Wolfram scalar *)
Φ12 = QuantumState["PsiMinus"]["Dagger"][
  QuantumTensorProduct[σ[θ1, φ1], σ[θ2, φ2]][QuantumState["PsiMinus"]]
]["Number"]

(* CHSH winning probability via MultivariateDistribution *)
NProbability[BitAnd[x, y] == BitXor[a, b],
  {a, x, y, b} ∼ QuantumCircuitOperator["CHSH"][]["MultivariateDistribution"]]
```

`["MultivariateDistribution"]` on a `QuantumMeasurement` returns a `MultivariateDistribution` you can use in `Probability[]`/`NProbability[]` directly — far cleaner than manually summing over outcomes.

### Counts simulation via SimulatedMeasurement

```
counts = Counts[exp["SimulatedMeasurement", 1024]]    (* <|"01"->498, ...|> *)
```

`["SimulatedMeasurement", n]` simulates `n` shots from the *exact* probability distribution. Useful for matching expected QPU statistics. Note: this is sampled from the exact distribution, not a noisy simulation.

### Parametric Hamiltonians via "ParameterSpec"

For complex parametric Hamiltonians (transmon, Rabi, QuEra Aquila):

```
H[t_] := QuantumOperator[Cos[ω t] QuantumOperator["X"] +
  Sin[ω t] QuantumOperator["Y"], "ParameterSpec" -> t];
ψf = QuantumEvolve[H, ψ0, {t, 0, T}];
state = ψf[1.5]   (* substitute time *)
```

The `"ParameterSpec" -> t` makes time a hidden parameter. `QuantumEvolve` returns an interpolating function — call `ψf[time]` for the state at that time.

### Hardware backend chain

```
qpy = BaseEncode @ qc["Qiskit"]["QPY",
   "Provider" -> "IBMProvider", "Backend" -> ibm_brisbane];
ibmq["RunCircuit", {"QPY" -> qpy, "Backend" -> ibm_brisbane}]
```

The QPY route bypasses Wolfram's default Aer simulator and submits directly to IBMQuantumPlatform. For AWS Braket: use `qc["Qiskit"]["QASM", "Provider" -> "AWSBraket"]` then `braket["CreateQuantumTask", "DeviceArn" -> ..., "Action" -> qasm_program]`. For QuEra Aquila (analog/Rydberg): construct a custom `BraketAHS` association — `AquilaHamiltonian[sites, {Ω, Δ, φ}]` builds the analog Rydberg Hamiltonian (see `Notebooks/QuEra_WhitePaper.nb`).

### Stabilizer formalism (Quantum_Error_Correction.nb)

```
(* 5-qubit code *)
s1 = QuantumOperator["XZZXI"];  s2 = QuantumOperator["IXZZX"];
s3 = QuantumOperator["XIXZZ"];  s4 = QuantumOperator["ZXIXZ"];

(* Code-subspace projector *)
projector = QuantumCircuitOperator[{
  (1/2)(("I" -> Range[5]) + "XZZXI"),
  (1/2)(("I" -> Range[5]) + "IXZZX"),
  ...
}];
```

Pauli-string operators with the `"XZZXI"` shorthand build N-qubit Pauli operators from one string. Combined with `+`/`-` and rational coefficients, you can build code projectors algebraically.

### Tensor network: ContractTensorNetwork beats default for large circuits

```
net = circuit["TensorNetwork"];
ContractTensorNetwork[net]                       (* optimized via EinsteinSummation *)
ContractTensorNetwork[net, Method -> "Naive"]    (* edge-order, slower *)
```

The default uses `EinsteinSummation` to find an optimal contraction order. `Method -> "Naive"` is `~5x` slower on a 10-gate circuit. To start from a non-register state: `circuit[ψ0]["Tensor"]` or `InitializeTensorNetwork[net, ψ0["Tensor"]]`.

### QAOA-in-QAOA divide-and-conquer

For Max-Cut on graphs > 6 nodes, hierarchical QAOA:

```
(* Phase 1 — Divide *)
subgraphs = partition[graph];
subQAOA = (build QAOA per subgraph);
solution = NMaximize per subgraph;

(* Phase 2 — Conquer *)
bigGraph = (subgraphs as nodes, weighted by inter-subgraph edges);
bigQAOA = (QAOA on the coarse-grained graph);

(* Phase 3 — Flip & correct *)
flip = Keys[Select[bigSolution, MatchQ[#, Rule[_, 1]] &]];
finalSolution = correctedSubgraphSolution
```

Documented in `OngoingProjects/Courses/QuantumOptimization/Lesson3_QAOA.nb` and `Quantum_Hackaton_LATAM_OpenQuantumInstitute.nb`.

### Phase oracle vs Boolean oracle

```
QuantumCircuitOperator["BooleanOracle"[bf]]   (* needs |q⟩ ancillary, |x⟩|q⟩ → |x⟩|q⊕f(x)⟩ *)
QuantumCircuitOperator["PhaseOracle"[bf]]     (* no ancilla, |x⟩ → (-1)^{f(x)} |x⟩ *)
```

For Grover, the phase oracle variant is shorter and consumes one fewer qubit. Use `"GroverPhase"[bf]` over `"Grover"[bf]` when ancilla budget matters.

### Programmatic Hamiltonian construction for Max-Cut

```
index = MapApply[List, EdgeRules[g]] /. x_Integer :> {x, x};
op[I] = StringRepeat["I", VertexCount[g]];
sum = Total[(1/2)(op[I] - StringReplacePart[op[I], "Z", #]) & /@ index];
Hc = QuantumOperator[sum]
```

Builds the Max-Cut cost Hamiltonian from a `Graph` programmatically. Generalizes to weighted graphs (multiply by edge weights), QUBO problems, and any sum-of-Pauli-strings Hamiltonian.

### Notebook-only QF idioms (not in docs)

- `qs["Number"]` — when state is a 0-qudit "scalar" object, returns the underlying complex number (vs `["Scalar"]` which returns a Wolfram-Lang scalar expression).
- `qs["BlochVector"]` returns 3-vector for qubit; `qs["BlochCartesianCoordinates"]` is a synonym.
- `Composition @@ Reap[stateEvolution[qs]][[-1, 1]]` — using `Sow`/`Reap` to capture an emitted circuit during a `FoldList` traversal.
- `qc[<|θ -> π/3|>]` — substitute parameters via `Association`. Equivalent to `qc[π/3]` if there's only one parameter; required when there are multiple.
- `qc[Method -> "Schrodinger"]` — explicitly pick the Schrodinger application path (vs default tensor-network); useful when the TN path is slow for very small circuits or for fidelity comparison.
- `qc[]["Reverse"]` — reverse the qudit order of an output measurement (matches IBM's little-endian convention).
- `QuantumLabelName[qc]` — returns the textual label list of a circuit's gates (used to extract a `QuantumShortcut`).
- Sub-circuit slicing: `qc[[2;;-2]]` works on `QuantumCircuitOperator` and returns a sliced circuit (drops the wrapping measurements/state-prep).
- `QuantumState["..."]["Operator"]` returns the conventional QuantumOperator (vs `["Projector"]` which returns the higher-dim "Pure map" form).
- For state preparation: `QuantumCircuitOperator["StatePreparation"[state]]` or `QuantumCircuitOperator["QuantumState"[state]]` builds a circuit whose `[]` output equals `state`.

### `QuantumChannel` symbolic algebra

```
QuantumChannel["BitFlip"[p], {2, 3}]                                 (* on qubits 2, 3 *)
QuantumChannel["BitFlip"[p], {2}] @ QuantumChannel["BitFlip"[p], {3}]   (* parallel *)
```

Channels accept symbolic parameters and compose via `@`. The result is another `QuantumChannel` (uses the Stinespring dilation internally). Symbolic parameters in the rate are supported and propagate through `QuantumEvolve`.

### `qc[<|param -> v|>]["AmplitudesList"]` for fast QPU comparison

```
quantum_pred = qc[]["ProbabilitiesList"];
qpu_results = ibmq["JobResults", {"JobID" -> jobid}];
counts = qpu_results["results", 1, "data", "counts"];
qpu_pred = N @ Normalize[Values @ Normal @ counts, Total];
BarChart[Transpose[{quantum_pred, qpu_pred}]]
```

Standard pattern for comparing Wolfram's exact prediction against IBMQ/AWS hardware results. The `Normalize[..., Total]` makes counts into a probability distribution.

### State preparation via recursive multiplexer (state-prepV3)

```
RZY[vec_] := Module[{a, b, φ, ψ}, {{a, φ}, {b, ψ}} = AbsArg[vec];
  {{"RZ", φ - ψ}, {"RY", 2 ArcSin[(a - b)/(Sqrt[2 (a^2 + b^2)])]}}];

multiplexer[qs_, permute_] := ...;
qc = Composition @@ Reap[stateEvolution[qs]][[-1, 1]];
qc[] == qs    (* recovers the original state *)
```

For a state vector of length 2^n, builds an n-qubit circuit whose action on `|0...0⟩` produces that state. Pattern from `Notebooks/state-prepV3.nb`. More efficient than `QuantumCircuitOperator["StatePreparation"[state]]` for non-symbolic states.

### Lindblad shortcut via `Ls -> γs` Rule

From `OngoingProjects/Courses/QuantumOptimization/Tutorials/2025-02_Training_QF3.nb` (and 2025-03 mirror):

```
(* Compact: jump operators + rates as a Rule, no explicit "Liouvillian" wrap *)
ρt = QuantumEvolve[H, Ls -> γs, ρ0, {t, 0, 1}]

(* Or with implicit unit rates *)
ρt = QuantumEvolve[H, Ls, ρi, t]   (* Ls is a flat list, no rates given *)

(* T1/T2 relaxation pattern *)
Ls = {Sqrt[1/T1] QuantumOperator["J+"], Sqrt[1/T1] QuantumOperator["J-"],
      Sqrt[1/T2] QuantumOperator["JZ"]};
ρt = QuantumEvolve[Ω0 QuantumOperator["JZ"], Ls, ρi, t]
```

This is a 5th, even more compact, equivalent of the 4 Lindblad forms above.

### `"Hamiltonian"[H, Ls, γs]` operator (master-equation-ready)

From QuantumOptics dissipative-Jaynes-Cummings:

```
ℋ = QuantumOperator["Hamiltonian"[H, {decayAtom, lossField}, {γ1, γ2}]];
evolved = QuantumEvolve[ℋ, ψ0, {t, 0, 2}]
```

`QuantumOperator["Hamiltonian"[...]]` packs a coherent `H` plus jump operators + rates into a single object that `QuantumEvolve` understands without needing the `ⅈ` factor (in contrast to `Liouvillian`-form, which does).

### Two-mode `["Bend"]` shortcut

From `OngoingProjects/Courses/QuantumOptics/Mads_Edits/BookChapter6.nb`:

```
{a1, b1} = Array[AnnihilationOperator[{#}] &, 2];
a1["Bend"] == a1 @ b1     (* same operator on both modes simultaneously *)
```

`["Bend"]` on a single-mode operator constructs the symmetric two-mode action — useful for two-mode squeezing / beam-splitter algebra.

### Stabilizer subsystem (post-2026-05-06 audit)

The PauliStabilizer / StabilizerFrame / GraphState / CliffordChannel / hybrid-interop kernel files in `QuantumFramework/Kernel/Stabilizer/` follow these idioms after the 2026-05-06 WL-idiom audit:

- **`Internal``Bag`` instead of `AppendTo`** for sequential greedy algorithms (e.g., `agPickIndependent` in `Constructors.m`, the AG-canonical circuit synthesis `gateBag` in `Conversions.m`'s `ps["Circuit"]`). Amortized O(1) per stuff vs O(n) per `AppendTo`.
- **`Fold` instead of `Do` + accumulator mutation** when the algorithm is genuinely sequential but each step is a function of prior state (e.g., `agExtendToSymplecticBasis` in `Constructors.m`, `stabilizerRowSumAGPhase` in `CliffordChannel.m`). Reads as a state-monad transformation; lets the kernel decide when to specialize.
- **`MapThread` for vectorized per-qubit AG g-function** in `rowsum` (`Measurement.m`) and the row-sum phase tracking (`CliffordChannel.m`). Replaces `Total @ Table[..., {j, n}]` with a single MapThread over the X / Z bit vectors. Faster on packed arrays.
- **Per-row XOR via `BitXor` over `SparseArray`** in `Compose.m` and `GateUpdates.m`. The Phase 1 design choice (`synthesis §9`).
- **No `Quiet`** anywhere in the subsystem. Per the user's `feedback_no_quiet` memory rule; if a message is expected, declare it in a `VerificationTest`'s third arg, never suppress.
- **`Catch` + `Throw` instead of `Break`** for short-circuit search (e.g., `agPickIndependent` early termination, `stabilizerPauliFromMatrix` first-match in `HybridInterop.m`). `Break` doesn't compose with nested constructs; `Throw` does.
- **Pattern-matched dispatch** (e.g., the named-code constructors in `NamedCodes.m`, `["Random", n]` named pattern, `op_ -> order_` rewrite in `PauliStabilizer.m:80`) instead of `If`/`Switch` cascades.
- **Sequential `Do` is acceptable** for state-machine algorithms: AG canonical circuit synthesis (`Conversions.m:ps["Circuit"]`) and the AG measurement primitive (`Measurement.m:rowsum`-based `["M"]`). Documented in the file headers.

When extending the Stabilizer kernel, run `wolframscript -file Tests/RunTests.wls Stabilizer` after every change. The test suite recurses into subdirectories (`Tests/Stabilizer/AuditMatrix.wlt` is the canonical coverage matrix).


### Polarization × path joint basis (PBS)

```
polBasis  = QuantumBasis[{"H", "V"}];
pathBasis = QuantumBasis[{e_0, e_1}];
PBS = QuantumOperator[<<4×4 matrix>>, {polBasis, pathBasis}];
```

`QuantumOperator[matrix, {basis1, basis2}]` constructs an operator on a *list* of bases that aren't computational. Lets you build a polarizing beam splitter where horizontal-vs-vertical polarization couples to different path modes — the resulting `QuantumOperator` then composes naturally with `QuantumState["+", polBasis] // QuantumState["0", pathBasis]`.
