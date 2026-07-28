# Generic array module + QuantumState container generalization - design and plan

**Goal** (from PROMPT.md): replace `QuantumState`'s hard-coded `SparseArray` amplitude
container with a generic array-like contract, backed by **one domain-neutral array
module** usable anywhere in the framework - not a state-specific layer.

Supported container classes:

1. **Explicit** - `SparseArray`, packed arrays and plain `List`s, `NumericArray`,
   `StructuredArray`s (`SymmetrizedArray` etc.): element-addressable, materializable.
2. **Lazy parametric** - an array-valued expression evaluated on demand, first target
   `InterpolatingFunction[...][t]` as returned whole by NDSolve, so `QuantumEvolve` no
   longer splits one multivalued InterpolatingFunction into 2^n (4^n for Lindblad)
   scalar ones.
3. **Symbolic** (to a workable degree) - `VectorSymbol` / `MatrixSymbol` / `ArraySymbol`,
   atomic symbols with `$Assumptions`-registered dimensions, and inactive/structural
   expression trees (`Inactive[TensorProduct]`, `Verbatim[Transpose][...]`,
   `TensorContract[...]`) whose shape is computable without materialization.

Plus: extensive tests for the new design; profile the new `QuantumEvolve`.

**Precedent to follow:**
`/Users/swish/src/wolfram/TensorNetworks/TensorNetworks/Kernel/IndexArray/ArrayUtilities.wl`
(the `` Wolfram`TensorNetworks`IndexArray` `` context). It already implements the symbolic
tier: `ArrayDimensions` recurses structurally through `Inactive[D]`, `Transpose`, `Plus`,
`Inactive[TensorProduct]`, `TensorContract`; `ArrayName`/`ArrayPart`/`setDimensions`
handle `VectorSymbol | MatrixSymbol | ArraySymbol` and assumption-registered atomic
symbols (rewriting `$Assumptions` in place); `ArrayTranspose`/`ArrayContract` compose
permutations and keep results in simplified inactive form (`SimplifyArray`);
`ZeroArrayQ`/`squareMatrixQ` are shape-level predicates. The QF module mirrors this API
shape and naming, and adds what the precedent lacks: container classification, the
SparseArray explicit-values accessor tier, the lazy-parametric tier, and
pack/materialize conversions.

**Decisions taken (user-confirmed):**
- **Transparent container + generic array module.** No new wrapper head:
  `QuantumState[state, basis]` keeps holding the raw array; the module supplies generic
  dispatch for every array operation. `qs["State"]` returns whatever container the state
  holds. The module itself knows nothing about quantum states.
- **Preserve the given container.** The constructor stops coercing to `SparseArray`;
  whatever supported container arrives is stored as-is (unpacked numeric lists may
  normalize to packed where free). Internal producers keep emitting what they emit.
- **Generic by design, precedent-aligned.** Named and shaped after the TensorNetworks
  `ArrayUtilities` precedent so the two can converge later (see open question).
- **Separate paclet, installed from the cloud.** The module is NOT a QF Kernel file:
  it is its own paclet `Wolfram/Arrays` in its own repo at
  `~/src/wolfram/Arrays`, laid out TensorNetworks-style with the paclet nested
  one level deep (`Arrays/PacletInfo.wl`, `Kernel/`, `Documentation/`; repo root holds
  `Tests/`, README, build scripts). QF installs it from a public cloud object exactly as
  it installs `Wolfram/TensorNetworks`, and imports `` Wolfram`Arrays` `` per file.
- **Lazy tier is a head registry, not an InterpolatingFunction special case.** Survey of
  built-in candidates (live-kernel verified on WL 15.0, see
  `Container-Candidates-Survey.md` in this folder): v1 admits InterpolatingFunction
  applied forms (+ their Derivative forms), ParametricFunction (all three arities - the
  two-stage flagship: parameters then time, kernel-cached solves), unapplied `Function`
  closures (analytic gates, probe-based shape discovery), and array-valued `Piecewise`
  (the one form whose unevaluated `Dimensions` just works). v2 registry candidates:
  BSplineFunction/BezierFunction (real-only, arity+domain guards), TimeSeries as
  accept-and-convert interop, LinearSolveFunction if a lazy-OPERATOR tier is added.
  Rejected with ingest adapters: Tabular/Dataset (Normal+pack), ByteArray
  (NumericArray lift), QuantityArray (QuantityMagnitude strip), FunctionInterpolation
  (corrupts array input in 15.0). Cross-cutting v1 requirements from the survey:
  intercept Part/Map/Dimensions on every inert applied form (they operate on the
  expression tree, not the array); document the NumericQ evaluation trigger (exact
  values like Pi do NOT stay lazy, only non-numeric symbols); one-substitution
  semantics (a single `/.` = one whole-array eval); repack unpacked eval outputs.

**QML arc (net-backed lazy arrays, design in progress):** a further lazy class where an
array is an unevaluated NetGraph with NetArrayLayer leaves, and TensorNetworkContract
gains a backend emitting an Inactive NetGraph (Transpose/Reshape/Dot layer motifs per
pairwise contraction step) - making circuit tensor networks differentiable/trainable.
The net framework is real-only, and that is NOT a blocker: **QML runs in the phase-space
picture**, where the whole network is real - kernel-verified in QF: the Bell state's
Wigner vector and the CNOT and RZ(0.7) phase-space kernels all have max|Im| = 0 exactly
(QuantumWignerTransform), and TensorNetworkCompile already has the phaseSpaceQ branch,
so a NetGraph backend inherits phase-space routing for free. Complex-number support in
nets (or the real/imag doubling lift) becomes a later nicety, not a prerequisite.
Feasibility probe results and the full design section land separately.

**Parallel arc, explicitly out of scope now:** the packed-array TN fast paths
("PackTensors" / "PackedFold" per the ethereal-cat notes). This module makes those
easier later (its pack/materialize primitives are exactly what the TN boundary needs,
and a packed state stays packed instead of paying the sparse tax), but no TensorNetwork
bypass shortcuts are implemented here.

---

## 1. Where the SparseArray contract lives today (verified inventory)

Two chokepoints define the contract:

1. **Validity predicate** - `QuantumState/QuantumState.m:15`
   `quantumStateQ[QuantumState[state_SparseArray ? stateQ, qb_]] := Length[state] === qb["Dimension"]`
   The head is baked into validity. Anything else fails `QuantumStateQ`, and since every
   property dispatch is guarded on it, a non-SparseArray container makes the whole
   property API silently inert (returns unevaluated, no message).
2. **Coercion rule** - `QuantumState/QuantumState.m:138-139`
   `QuantumState[state_ ? stateQ, basis_] /; ! SparseArrayQ[state] && state =!= {} :=
    ... SparseArray[state, Length[state]] ...`
   The sole intake funnel; ~80-100 `QuantumState[<computed array>, basis]` reconstruction
   sites across the kernel all pass through it.

Supporting contract surface:

- `stateQ` (`Utilities.m:84`): `VectorQ[state] && Length[state] >= 0 || SquareMatrixQ[state]` -
  the rank test; neither `InterpolatingFunction[...][t]` nor a `MatrixSymbol` satisfies it.
- `"StateType"` (`QuantumState/Properties.m:82-88`): `VectorQ`/`SquareMatrixQ` on the
  container; a lazy or symbolic container falls to `"UnknownType"` and poisons the
  `"Type"` cascade.
- `"NumericQ"`/`"NumberQ"` (`Properties.m:98-100`): `Replace[..., sa_ ? SparseArrayQ :> sa["ExplicitValues"]]`.
- **RAW SparseArray subfield reads on state data** (break silently under any other
  container): `Properties.m:125, 129` ("Amplitude"), `:357, 361` ("GraphRule"),
  `:609` ("Disentangle"), `:621, 639` ("Decompose"), `:695` ("Compress");
  `QuantumOperator/Properties.m:393-409` ("OrderedFormula"). Plus WRAPPED sites
  (re-wrap in `SparseArray` first - already container-tolerant): `Properties.m:173,
  347-349, 398, 675-677`, `QuantumMeasurement.m:211-214`, QSVT/Experimental/Optimization.
- `SparseArrayFlatten` (`Utilities.m:288-313`): 25+ call sites in 10 files; already a
  3-clause head dispatch (CSR fast path for rank > 11, scalar, generic `Flatten`).
  **This is the existing seam for structural ops.**
- `"State"` getter (`Properties.m:74`): the single structural read of slot 1 - a one-line
  interception point; ~90 `["State"]` call sites all route through it.
- Element-wise maps that densify: parameter substitution
  (`QuantumState.m:427`: `Map[ReplaceAll[rules], state, {1|2}]`), `"Simplify"` family
  (`Properties.m:368`), degenerate scan in `"Type"` (`Properties.m:489`).
- History: the container was `Developer`ToPackedArray` originally, switched wholesale to
  `SparseArray` in f07e5825 ("trying to optimize stuff", Oct 2021). There is exactly one
  `NumericArray` entry point (`QuantumOperator.m:227`) and it immediately `Normal`s it.

QuantumEvolve specifics (`QuantumEvolution.m`):

- NDSolve returns ONE array-valued `InterpolatingFunction`; `ExpandInterpolatingFunction`
  (`:254-262`, invoked at `:180-182`) explodes it into one scalar `Interpolation` per
  amplitude, each duplicating the full time grid: O(d^2 * nGrid) memory for Lindblad.
  It exists only because the container contract cannot hold `if[t]`.
- The inverse already exists: `MergeInterpolatingFunctions` (`:216-251`) collapses
  per-entry IFs into a single array-valued closure to feed NDSolve. The new container
  is essentially "store what MergeInterpolatingFunctions produces".
- Parameters live in the basis (`"ParameterSpec"`); `qs[t0]` does element-wise
  `ReplaceAll` + full reconstruction; the property cache is disabled entirely for
  parametric objects (`ParameterArity == 0` guards).
- Ranked profiling targets: the IF explosion; `ProbabilitiesList` embedding the full
  `Total` in each entry (O(d^2) scalar IF calls per plot sample); the disabled cache;
  per-timestep basis rebuild in `qs[t0]`; `Piecewise` overhead in the merge path.
- No `QuantumEvolve.wlt` exists; only `Tests/Liouvillian.wlt` covers evolution routes.

---

## 2. Design

### 2.1 The module: paclet `Wolfram/Arrays` (domain-neutral)

The module is a standalone paclet (repo `~/src/wolfram/Arrays`, paclet nested at
`Arrays/Arrays/`, installed into QF from the cloud). **It imports nothing
from the quantum layer and mentions no quantum head** - a self-contained abstraction
over WL array containers. Naming follows the TensorNetworks precedent (`Array*`,
CamelCase, PackageExport with ::usage strings, documented via
create-wolfram-documentation). Every function head-dispatches with a generic fallback,
so a new container type is one clause per function, not a kernel-wide hunt; the lazy
tier specifically dispatches on head-of-head through a registry table
({InterpolatingFunction, ParametricFunction, Function, Piecewise} in v1 - see the
survey doc for the verified per-head protocols).

Intended consumers beyond `QuantumState` (why this must stay generic):

- `QuditBasis` representations - the same SparseArray coercion gate exists at
  `QuditBasis/QuditBasis.m:58-59` and the same `ExplicitValues` idioms at
  `QuditBasis/Properties.m:355, 357`; migrating it later is the same one-line swap.
- `QuantumOperator`'s `NumericArray` intake (`QuantumOperator.m:227`), and the
  Qiskit/QuEST/Classiq FFI wire paths that produce/consume `NumericArray`.
- The TensorNetwork boundary (the ethereal-cat arc): its pack/densify decisions are
  `ArrayMaterialize`/`ArrayPack` calls on tensors.
- Symbolic derivations: a `QuantumState` over `MatrixSymbol`/`VectorSymbol` amplitudes
  with known shape, usable with `$Assumptions`-driven simplification.
- Any future container (a generating function, an out-of-core store) becomes one
  dispatch clause per operation. GPUArray is already surveyed (survey section F.3):
  explicit-tier compute-native, ~10x on complex Dot at 512^2, but Metal compute is
  fp32-only (real storage Real32, complex Dot error ~3e-5), so it enters via the v2
  registry behind an availability guard and an explicit precision-consent guard;
  recognition must use GPUArrayQ, never the head (malformed calls stay inert with
  head GPUArray).

**Admission criterion (revised, user-confirmed):** a type is admissible to the explicit
tier iff it has a Dimensions equivalent (shape introspection without materializing) AND
a materialization path. Compute-nativeness (Dot/arithmetic working in-container) is NOT
an admission gate - it is a per-container capability flag; storage-only containers
(NumericArray, QuantityArray, TabularColumn, Tabular, ByteArray, ...) get
materialize-before-compute semantics, and where a reconstruction path exists the module
re-wraps results to preserve the container round-trip. Numericity is judged at module
level (numeric-with-units is numeric; unit-stripping is consumer policy). Genuine
rejection is reserved for types with NO shape or NO materialization (e.g. Compress
strings). Per-type verdicts under this criterion: see Container-Candidates-Survey.md
(revision in progress).

**Container classification**

```
ArrayContainerQ[a]    supported container? (any tier)
ArrayExplicitQ[a]     tier 1: shape-introspectable + materializable; spans
                      compute-native (SparseArray | packed | List | QuantityArray |
                      TabularColumn - both re-probe-verified with full native
                      arithmetic) and storage-only (NumericArray | Tabular | Dataset |
                      ByteArray | EventSeries | DataStructure "DynamicArray" |
                      StructuredArray atoms) sub-classes
ArrayComputeNativeQ[a] capability flag: in-container arithmetic/Dot verified for this
                      head; storage-only containers materialize before compute
                      (QuantityMagnitude for QuantityArray - Normal is the trap;
                      Normal[ev["Values"]] for EventSeries; snapshot-on-ingest for
                      DataStructure reference semantics)
ArrayLazyQ[a]         tier 2: array-valued inert applied form with >= 1 non-numeric
                      arg; head registry, v1 = InterpolatingFunction[__][args__]
                      (+ Derivative forms), ParametricFunction (3 arities),
                      unapplied Function, array-valued Piecewise (per survey)
ArraySymbolicQ[a]     tier 3: VectorSymbol | MatrixSymbol | ArraySymbol, atomic symbols
                      with Vectors|Matrices|Arrays assumptions, and structural
                      inactive trees (per precedent's ArrayDimensions clauses)
```

**Shape (never materializes)**

```
ArrayDimensions[a]    explicit -> Dimensions/tensorDimensions;
                      lazy IF  -> Head[a]["OutputDimensions"];
                      symbolic -> the precedent's structural recursion
                      (TensorDimensions, Transpose/Plus/TensorProduct/TensorContract
                      clauses, assumption lookup)
ArrayRank[a]          Length @ ArrayDimensions
ZeroArrayQ[a]         MemberQ[ArrayDimensions[a], 0]
ArrayNumericQ[a] / ArrayNumberQ[a]
                      SparseArray -> ExplicitValues test (as today); packed/NumericArray
                      -> True fast path; lazy/symbolic -> False
```

The state-domain rank contract ("Vector" | "Matrix" | "UnknownType", the
`Length[state] === qb["Dimension"]` check) stays in the QuantumState layer, built on
`ArrayDimensions` - the module classifies containers and shapes, never quantum
semantics.

**Accessors (fill the ExplicitValues/Positions API gap; explicit tier only,
Missing/fallback elsewhere)**

```
ArrayExplicitValues[a]     SparseArray native; dense -> nonzero values (wrap on demand)
ArrayExplicitPositions[a]  likewise
ArrayExplicitLength[a]     likewise
ArrayMaterialize[a]        Normal for explicit; lazy -> ExpandInterpolatingFunction
                           (the ONLY place it remains); symbolic -> unevaluated
ArrayPack[a]               best-effort packed conversion (2-arg ToPackedArray, Complex),
                           for numeric hot paths and the future TN arc
```

**Structural ops** (subsume the `SparseArrayFlatten` seam - it becomes a thin alias so
its 25+ call sites keep working)

```
ArrayFlatten*, ArrayReshape*, ArrayTranspose, ArrayContract, ArrayPart,
ArrayPad, ArrayConjugate, ArrayChop, ArrayMap[f, a, level], ArrayAllZeroQ[a],
SimplifyArray[a]
```

(* names colliding with System` or TensorNetworks exports get resolved at
implementation time - see open question 1; the precedent already owns
ArrayTranspose/ArrayContract/ArrayPart/SimplifyArray semantics and those are adopted
as-is, including inactive-form results for symbolic inputs. *)

Per-tier behavior:
- **Explicit**: pass through to the native op, preserving container where the op does
  (SparseArray stays sparse, packed stays packed; NumericArray/StructuredArray convert
  on ops WL does not support natively, result stored as produced).
- **Lazy**: materialize through the ONE choke (`ArrayMaterialize`), documented per op;
  structural ops on lazy containers are rare (the common path substitutes the
  parameter first).
- **Symbolic**: stay symbolic - shape-tracked inactive forms per the precedent
  (`ArrayTranspose` composes permutations, `ArrayContract` wraps in inactive
  `TensorContract`, `SimplifyArray` cleans up); ops without a sound symbolic form
  return unevaluated rather than guessing.

**The lazy fast path (the QuantumEvolve win)**

```
ArrayReplaceAll[a, rules]   substitution:
    lazy      -> a /. rules   (ONE array-valued IF evaluation for the whole array)
    SparseArray -> MapSparseArray-style over explicit values
    symbolic  -> ReplaceAll on the expression tree
    other     -> ReplaceAll / Map at level
```

`qs[t0]` on an evolved state becomes one IF call returning a packed array, instead of
2^n scalar IF calls - and the result is then an ordinary explicit (packed) state.

### 2.2 Chokepoint edits (the actual kernel diff, small by design)

1. `quantumStateQ` (`QuantumState.m:15`):
   `state_SparseArray ? stateQ` -> `state_ ? ArrayContainerQ` plus the state-side shape
   check, with the dimension read via `ArrayDimensions[state]`.
2. Coercion (`QuantumState.m:138-139`) -> normalization:
   supported containers pass through unchanged; plain lists normalize (packed where
   numeric, and the existing 0/1-dimension padding cases keep their current SparseArray
   behavior so edge-case semantics do not move).
3. `stateQ` (`Utilities.m:84`): shape test routed through `ArrayDimensions` so lazy and
   symbolic containers pass (rank contract itself unchanged, still state-side).
4. `"StateType"`, `"NumericQ"`, `"NumberQ"` (`Properties.m:82-100`) -> module calls.
5. The ~10 RAW `ExplicitValues`/`ExplicitPositions`/`ExplicitLength` sites -> accessors.
   (WRAPPED sites keep working as-is but migrate opportunistically.)
6. Parameter substitution (`QuantumState.m:427`) -> `ArrayReplaceAll`.
7. `"Simplify"` maps (`Properties.m:368-371`) and the `"Type"` degenerate scan
   (`Properties.m:489`) -> `ArrayMap` / `ArrayAllZeroQ` (this also resolves the noted
   rule-ordering ambiguity between the generic map and the "works on SparseArrays
   directly" fast clause).
8. `QuantumEvolve` (`QuantumEvolution.m:180-182`): stop calling
   `ExpandInterpolatingFunction`; store the solver's array-valued
   `InterpolatingFunction[...][parameter]` directly as the (lazy) container.
   `ExpandInterpolatingFunction` stays available for explicit materialization
   (used by `ArrayMaterialize`), and an option (`"Expand" -> True`) restores old
   behavior.

Not migrated in this arc (but now one-line swaps later): the `QuditBasis`
representation coercion, TensorNetwork tensor packing, Stabilizer/* (independent packed
world, stays untouched).

### 2.3 Behavior contract for parametric states (compat)

- `psi[t0]` -> one IF evaluation -> ordinary explicit state; everything downstream
  unchanged. This is the dominant usage (docs, README, Liouvillian tests).
- Symbolic-t reads (`Plot[Evaluate @ Normal @ psi["ProbabilitiesList"], ...]`) route
  through `ArrayMaterialize` - same results as today, no API change. A faster
  whole-array symbolic form (`Abs[if[t]]^2` un-expanded) is a documented follow-up, not
  in v1, to keep return shapes stable.
- Formatting on lazy states short-circuits via `"Type" === "Parametric"` before any
  materializing property (Purity/Entropy under TimeConstrained today).
- The property cache stays disabled for parametric states (existing guard, unchanged).

Symbolic containers get a minimal, honest contract in v1: construction, shape/type
properties, conjugate/transpose/tensor-product in inactive form, `Simplify`, and
substitution of concrete arrays via `ArrayReplaceAll`; materializing properties
(probabilities, plots, FFI) return unevaluated or fail loudly, never silently.

---

## 3. Phasing

**Phase 1 - paclet + chokepoints, zero behavior change.**
The `Wolfram/Arrays` paclet (classification, shape, accessors,
structural ops, its own standalone test suite and documentation) is built in its own
repo, installed into QF from the cloud, and QF edits 1-4 land importing it. SparseArray
path bit-identical (it dispatches to exactly today's code). Full QF suite must stay
green with internal producers still emitting SparseArray everywhere. The paclet's
suite pins the survey's unlisted-API dependencies (IF "OutputDimensions"/"Domain",
TimeSeries "ValueDimensions", Piecewise Dimensions-threading, structured-array
atomicity) so version breaks surface as test failures, not silent misclassification.

**Phase 2 - container preservation + RAW-site migration.**
Edits 5-7; constructor preserves packed/NumericArray/StructuredArray inputs; tests that
each property class works on each explicit container (construct, change basis, tensor
product, partial trace, measure, compare materialized results across containers), plus
symbolic-container construction/shape/inactive-op tests.

**Phase 3 - lazy container + QuantumEvolve.**
Recognition + shape-from-IF + `ArrayReplaceAll` fast path; edit 8. New
`Tests/QuantumEvolve.wlt`: all `Liouvillian.wlt` equivalence routes repeated against the
lazy container; `psi[t0]` equality lazy-vs-expanded; `ProbabilitiesList` parity;
`ParameterSpec` propagation; `"Expand" -> True` compatibility; phase-space evolve.

**Phase 4 - profiling.**
Instrument with the existing `profile[...]` hook (currently only 3 sites): change of
basis, tensor product, evolve solve/wrap, `qs[t0]`. New benchmark script
`OngoingProjects/Platform Comparison/scripts/evolve01_lazy_if.wls` on the
`bench_qf_common.wl` harness (`tm` = min-over-7, `RESULT|label|value` lines):
old-vs-new memory of an evolved state, `qs[t0]` latency, `Plot` sampling latency,
end-to-end `QuantumEvolve` wall time at n = 2..8 qubits (Schrodinger + Lindblad).

**Test conventions** (per existing suite): `VerificationTest` + kebab TestIDs,
`BeginTestSection` per topic, materialize before comparison, compare-to-zero with
`Chop`/`Simplify`, fully-qualified PackageScope names where needed. The runner
auto-discovers new `.wlt` files recursively; no registration needed.

---

## 4. Risks

1. **Validity relaxation has kernel-wide blast radius** (operators/channels/measurements
   all wrap states). Mitigation: Phase 1 changes no producer, so every stored container
   is still a SparseArray until Phase 2; the suite gates each phase.
2. **Silently-inert property API** when a predicate misses (existing failure mode, seen
   with the `{}` escape hatch today). The module predicates get their own unit tests
   first, and `ArrayContainerQ` failures in the constructor raise a message instead of
   passing through.
3. **Lazy-container ops that materialize by accident** defeat the point. All
   materialization funnels through `ArrayMaterialize`; a test asserts `psi[t0]` on an
   evolved state performs no expansion (e.g. by counting `Interpolation` calls via a
   hook).
4. **NumericArray/StructuredArray op gaps** in WL (e.g. SymmetrizedArray de-optimizing
   `Dot`, already noted in `Stabilizer/RandomClifford.m:36`). Contract: storage is
   preserved; compute may convert; results are stored as produced. Documented per op.
5. **`SparseArray[state, Length[state]]` edge-case semantics** (0/1-dimension padding
   normalization) must be preserved for list inputs - pinned by tests before the
   constructor changes.
6. **Symbolic tier scope creep.** "To some degree" means shape + structural ops +
   substitution; anything requiring element enumeration is out. The v1 contract in 2.3
   draws that line explicitly so reviewers can push back on it rather than discover it.

## 5. Open questions

1. **Relationship to TensorNetworks' `ArrayUtilities`.** QF already depends on
   `` Wolfram`TensorNetworks` `` for contraction. Options: (a) QF keeps its own module
   mirroring the API (no new coupling for core state storage; some duplication),
   (b) upstream the container/lazy/accessor tiers INTO
   `TensorNetworks`IndexArray`ArrayUtilities` and import it (single source of truth;
   makes QF's core depend on TensorNetworks). Recommendation: (a) now, converge later
   once the API stabilizes; resolve name collisions (`ArrayTranspose`, `ArrayContract`,
   `ArrayPart`, `SimplifyArray`, `ZeroArrayQ`) at implementation time by adopting the
   precedent's semantics verbatim so a future merge is a deletion, not a migration.
2. Should `QuantumOperator`'s `NumericArray` intake (`QuantumOperator.m:227`) stop
   `Normal`izing in Phase 2 (preserve NumericArray through operator storage too)?
   Recommendation: yes, it is the same one-line swap.

## 6. Deliverables checklist

- [ ] `Kernel/ArrayUtilities.m` with full standalone unit coverage
      (`Tests/ArrayUtilities.wlt`) across explicit / lazy / symbolic tiers
- [ ] Chokepoints 1-8 migrated; full suite green at every phase
- [ ] Container-preservation semantics for NumericArray / packed / List /
      StructuredArray / SparseArray
- [ ] Lazy `InterpolatingFunction[...][t]` container end-to-end through QuantumEvolve
      (no per-amplitude splitting), with `"Expand"` compatibility option
- [ ] Symbolic-container v1 contract (shape, inactive structural ops, substitution)
- [ ] `Tests/QuantumEvolve.wlt` (new; today zero direct coverage)
- [ ] `evolve01_lazy_if.wls` benchmark + recorded before/after numbers
- [ ] `profile[...]` instrumentation at the four hot sites
