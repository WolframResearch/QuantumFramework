# Architecture


The mental model of QuantumFramework — how objects relate, how their properties are computed, where the cost lives. **Load this first** before writing any QF code.

## 0. Package boot

- Entry: `Kernel/QuantumFrameworkMain.m` — loader; auto-installs `IBMQuantumPlatform` (**≥0.0.4** as of `ac32bda8`, was 0.0.2) + `Wolfram/TensorNetworks` (≥1.0.5) paclets if missing; **now also `Needs["IBMQuantumPlatform`"]` at load** so `ServiceConnect["IBMQuantumPlatform"]` works with no extra Needs from the user; sets `$ContextAliases["H`"] = "WolframInstitute`Hypergraph`"`.
- Main package: `Kernel/QuantumFramework.m:1-25` — defines `QuantumFrameworkOperatorQ`, `$QuantumFrameworkPropCache`, FE special-arg autocompletion.
- Init: `Kernel/Init.m:7` — memoizes `FromOperatorShorthand`.
- Internal helpers: `Kernel/Utilities.m` — ~400 lines of math primitives.

External dependencies:
- `WolframInstitute`Hypergraph` (via `H``  context alias)
- `Wolfram/TensorNetworks` paclet (auto-installed; provides MPS support)
- `IBMQuantumPlatform` paclet (auto-installed from local asset)
- `DocumentationSearch` (built-in)

## 1. Object map

*(Per-object detail to be filled as object subdir kernel files are read.)*

### Top-level constructors (from `PacletInfo.wl` `Symbols` list)
- `QuditName` — `QuditName[name, "Dual" -> False]` wrapper. Tracks dual flag (Bra vs Ket). `splitQuditName` flattens nested QuditNames; `groupQuditName` re-groups consecutive same-dual entries. Tensor product via `QuantumTensorProduct` upvalue. `$QuditZero = CircleTimes[]`, `$QuditIdentity = \[FormalCapitalI]`.
- `QuditBasis` — `QuditBasis[Association[{QuditName, qudit-index} -> Tensor]]`. Each entry maps a (name, index) pair to its tensor representation. `Dimension = Times @@ Dimensions`. ~36 named bases (`$QuditBasisNames`) covering Computational, Pauli (eigenbases), Spin/J, Bell, Fourier, Schwinger, Dirac, GellMann, Wigner/MIC, Bloch family, Wootters, SICs (Tetrahedron, QBismSIC, HesseSIC, HoggarSIC), Ivanovic MUB. **Has its own construction cache** `$QuditBasisCache` separate from property cache.
- `QuantumBasis` — `QuantumBasis[Association[Input, Output, Picture, Label, ParameterSpec]]`. Wraps two QuditBases (Output + Input). `Picture ∈ {Schrodinger, Heisenberg, Interaction, PhaseSpace}` is part of basis identity (PhaseSpace gets special handling in `QuantumState` constructor for Wigner-style transforms). `MatrixRepresentation = KroneckerProduct[OutputMatrix, InputMatrix]`. `Dimension = OutputDimension × InputDimension`.
- `QuantumState` — `QuantumState[SparseArray, QuantumBasis]`. State is always a SparseArray (vector for pure, matrix for mixed; both validated by `stateQ`). Property dispatch caches per-state with parameter-arity guard (parametric states bypass cache).
- `QuantumOperator` — `QuantumOperator[QuantumState, {outputOrder, inputOrder}]`. Thin wrapper: an operator is a state plus a 2-tuple of qudit orders. `NHoldRest`. Composition routes through `QuantumCircuitOperator` (`qo1[qo2] = QuantumOperator @ QuantumCircuitOperator[{qo2, qo1}]`).
- `QuantumMeasurementOperator` — `QuantumMeasurementOperator[QuantumOperator, targets]` (`:9`). The wrapped operator's NonPositive output orders represent the eigenvalue register. Type auto-detected: `"Projection"` (no NonPositive output, equal in/out dims), `"POVM"` (NonPositive output present), `"Unknown"`. 9 named POVMs (RandomHermitian, WignerMICPOVM, GellMannMICPOVM, Tetrahedron/QBism/Hesse/Hoggar SIC POVMs, RandomPOVM Haar/Bloch). 5 cache-bypass properties.
- `QuantumChannel` — `QuantumChannel[QuantumOperator]` in **Stinespring dilation form**. The wrapped operator has `OutputQudits >= InputQudits`; extra output qudits, encoded with **non-positive** order indices (`0` and below), are the environment to be traced. So a single channel internally lifts to an isometric operator + a deferred partial trace. Constructor accepts a list of Kraus operators (stacked into one isometry) or 8 named forms (BitFlip/PhaseFlip/Depolarizing/AmplitudeDamping/etc.). Channel-on-anything composition currently routes through `QuantumCircuitOperator` — direct channel composition code is commented out.
- `QuantumCircuitOperator` — `QuantumCircuitOperator[Association]` with required `"Elements"` key (`:13-14`). Elements are `QuantumFrameworkOperatorQ` instances OR Barriers (`BarrierQ` matches `"Barrier" | "Barrier"[order, ...] | "Barrier"[Span, ...]`). Default application `qco[qs]` routes through `TensorNetworkApply` (`:87-102`). Method options: `Automatic | "TensorNetwork" | "Schrodinger" | "QuEST" | "Qiskit" | "Stabilizer"`. Composition via `Composition`/`RightComposition` upvalues flattens nested circuits (`:151-159`). 17 properties bypass cache (`Properties.m:28-31`). Recent commit `cbbe9368` added `"Parameters"` option simplification.
- `QuantumHamiltonianOperator` — `QuantumHamiltonianOperator[QuantumOperator]` (`:9`). A QuantumOperator with explicit `ParameterSpec` for time evolution. Default parameter `\[FormalT]`. Specialized properties for adiabatic computation: `Gap[n=100]` sweeps 100 parameter values, `MinimumGap` via `FindMinimum`, `EvolutionOperator` solves Schrödinger equation symbolically via `DSolve`. Constructors: `H1 -> H2` linear interpolation, `EvolutionOperator[U(t), spec]` inverse via `Reduce`, `Ising[J]`, `TransverseIsing[J]`. Filename typo (`QuantumHamiltonialOperator.m`) doesn't affect runtime. **Still on disk at HEAD `bfddf5ed` but does NOT load**: the constructor file's `Package[...]` is commented out (`:1`), so `QuantumHamiltonianOperator` resolves to `` Global` `` with 0 DownValues and stays unevaluated. Not in PacletInfo Symbols list; not reachable via context either. See `activity-status.md` and `function-index.md` for the corrected status.

### Stabilizer subsystem (`Kernel/Stabilizer/`, 6 public heads at `ac32bda8`; 20-file subdir since the June 2026 perf rounds)
A self-contained Aaronson-Gottesman tableau simulator, rebuilt from the old monolithic `Kernel/PauliStabilizer.m` into a subdir (17 files at `ac32bda8`; **20 now**, after `Packed.m`, `Compiled.m` and the perf-round splits). Design choice that shapes everything: integration with the rest of QF is via **UpValues on the Stabilizer heads** (`QuantumState[ps]`, `qo[ps]`, `qmo[ps]`, `qc[ps]`), NOT a Picture flag or QuantumBasis wrapper — the latter would route through `KroneckerProduct[Output,Input]` + `MatrixInverse` and pay O(2ⁿ)…O(8ⁿ), defeating the formalism's O(n²)–O(n³) advantage. Heads:
- `PauliStabilizer[<|"Signs", "Tableau"|>]` — canonical `{2,n,2n}` tableau (X/Z blocks; n destabilizer rows then n stabilizer rows). Signs may be symbolic F₂ polynomials (SymPhase). The Clifford generators are direct tableau mutations (`GateUpdates.m`); composition is symplectic mod-2 matrix multiply with a precomputed rank-4 `phaseLookup` SparseArray (`Compose.m`). Materialization to a `QuantumState`/`QuantumOperator` is O(4ⁿ) and carries a recovered `"GlobalPhase"`. **Second representation (June 2026, `Packed.m`):** a *packed* form `<|"PackedX","PackedZ","Signs","Qubits"|>` storing 62 tableau bits per machine word (chunk-major), with `psConcreteFastQ` gating the fast path; the hot gate kernels (`packedGateH/S/...`), single-qubit + Pauli-string measurement (`Measurement.m`/`PauliMeasure.m`), and formatting all operate directly on packed words, and `ps["Tableau"]` lazily reconstructs the canonical array for any other property. Single gates return packed objects; **`SameQ` across packed vs canonical fails; compare `["Tableau"]`/`["Stabilizers"]`.** Symbolic-sign / `StabilizerFrame` states fail `psConcreteFastQ` and stay on the canonical path.
- `StabilizerFrame[<|"Components" -> {{c, ps}, ...}|>]` — `Σ cᵢ|sᵢ⟩`; the closure object for non-Clifford P/T gates (which double the term count). **Gate-built frames also carry a `"Paulis"` key** (`61bc0e39`): one relating Pauli per component `{xbits, zbits, coeff∈{1,-1,I,-I}}` under the contract `|componentᵢ⟩ = matrix[Pᵢ]·|reference⟩`, maintained O(n)/gate by Clifford conjugation (`P→G·P·G†`) and `Z`-multiplication. The dense `["StateVector"]` branches on this key: with it, materialize the reference component once and relate the rest (phase-coherent up to one global phase, same contract as a bare `PauliStabilizer`); without it (hand-built or `Plus`-combined, where components need not share a tableau) sum each component independently. Only the inner-product phase between two distinct gauged frames stays gauge-limited.
- `CliffordChannel[<|"UA","UB","c",...|>]` — Yashin25 Boolean-tableau CPTP map; composition via F₂ null-space + AG/|Φ⁺⟩ phase tracking.
- `GraphState[<|"Graph","VOPs"|>]` + `LocalComplement` — AndBri05 graph-state representation (VOP tracking stubbed at identity).
- `StabilizerStateQ` — predicate.

`PauliStabilizerApply[qco, ps]` (in `Stabilizer/PauliStabilizer.m`) is the `Method -> "Stabilizer"` circuit dispatcher: returns `$Failed` (one `::nonclifford`) on a non-Clifford gate, or a `StabilizerFrame`/`Plus` for the P/T boundary. **Since the June 2026 round 2 (`Compiled.m`)** it routes an all-Clifford concrete circuit through a single `Compile`-to-C **bulk gate fold** (`ps["ApplyCircuit", specs]`) over a generator-major packed tableau (~0.2-0.3 µs/gate, flat in n; auto-selects `"C"`/`"WVM"` by compiler availability), falling back to the per-gate fold for anything unencodable. `ApplyCircuit` returns a **canonical** (rank-3) object, unlike single gates which return packed ones. The closed-form `PauliStabilizer[n]` register constructor and the packed measurement kernels are the other round-2 pieces; `RandomClifford.m` was hardened (`\mathbb{F}_2` inverse + packed matmul; exact-rational Mallows ratio) so `PauliStabilizer["Random"[n]]` no longer blows up near n≈500 or underflows past n≈538. Two pre-existing AG measurement bugs (destabilizer rowsum omitted; non-abelian Pauli-string group) were fixed in passing (`1c524432`/`52d3ea89`/`730abcec`). Full story: `OngoingProjects/Platform Comparison/QF-Stabilizer-Optimization-Report.md`; benchmarks 2.5-5.3× of Stim, ahead of QC.jl at n=1000.

### QASM hub + IBM jobs (new at `ac32bda8`)
- **QASM hub** (`QuantumCircuitOperator/QuantumQASM.m`): `QuantumQASM` is the single import/export entry point. **Import is native WL** (`OpenQASMImport.m`, depends only on the OpenQASM spec, not on any qiskit release); **export default is native WL emission** (`qasmEmit*` workers in the two `Properties.m` files). Any extra arg routes export through qiskit (faithful dump or transpile). The public `["QASM"]` properties on operator/circuit/QiskitCircuit all funnel here. `QiskitTarget` carries a validated device-spec, rebuilt at point-of-use (no cross-version pickle).
- **IBM jobs** (`QuantumCircuitOperator/IBMQuantum.m`): `IBMJobSubmit` submits async through qiskit's SamplerV2/EstimatorV2 (`qiskitPrimitiveSubmit` in `Qiskit.m`); `IBMJob` is a lossless handle storing raw service responses under `"Raw"` with WL-native accessors layered on top. Sampler counts are reordered into ascending-qubit order via a captured classical-bit→original-qubit map (`qiskitReorderByQubit`), so a permuted/transpiled measurement decodes correctly. Since `ac32bda8` three changes (`a44c1eea`):
  - **Pre-submit option validation** (`ibmValidatePrimitiveOptions`, `IBMQuantum.m:52`, called at `158`): the user's `PrimitiveOptions` are checked against qiskit's own `SamplerOptions`/`EstimatorOptions` dataclass schema **before** the slow transpile + submit, so a bad/mis-nested key fails in milliseconds with a message naming the offending key, the primitive, the valid keys at that level, and (for the common sampler↔estimator mix-up) a cross-primitive hint.
  - **Native result decoding** (`ibmSamples`, `IBMQuantum.m:357`): the completed-job sampler path now decodes the RuntimeEncoder-serialized `PrimitiveResult` envelope directly in WL. `ibmUnwrap` strips `<|"__type__","__value__"|>` layers, `ibmNpyDataBytes` inflates a zlib'd `.npy` (`Developer`RawUncompress` + header skip), `ibmBitArraySamples` unpacks a `BitArray`, and per-register outcomes are bit-packed into one integer per shot. Replaces the old `Interpreter["HexInteger"]` hex-string path.
  - **Properties vs Actions split** (`IBMQuantum.m:514,528`): `IBMJob["Properties"]` (`$ibmProperties`) lists only **local, network-free** accessors, so `AssociationMap[job, …]` never re-queries; the three accessors that need a live connection (`"Refresh"`, `"Cancel"`, `"ExecutedCircuit"`) moved to the new `IBMJob["Actions"]` set (`$ibmActions`). `"ExecutedCircuit"` is thus no longer in `"Properties"` (it downloads a server-transpiled circuit from a presigned VPC-private COS URL that is slow or 404s). Each Action is still callable by name; the unknown-property message now distinguishes the two sets.

### Two-qubit KAK (Cartan) decomposition (`QuantumOperator/Properties.m`)
`qo["KAK"]` / `qo["Decompose"]` factor a 4×4 unitary in the magic (Bell) basis where local SU(2)⊗SU(2) becomes real SO(4). Numeric (`TwoQubitKAK`) and structured-symbolic (`TwoQubitKAKSymbolic`, time-bounded + probe-substitution self-check) paths. **Output is correct up to a global phase** (the residual-phase projection does not always null it — see `mistakes.md`).

### Dimension semantics (subtle)

For both `QuantumBasis` and `QuditBasis`:
- `qb["Dimension"]` = total *element-space* dimension. For QuantumBasis: `OutputDimension × InputDimension`. For state-space (input dim 1) it equals OutputDimension; for operator-space (in==out) it's the square.
- `qb["OutputDimension"]` / `qb["InputDimension"]` — per-side.
- `qb["Dimensions"]` — list of per-qudit dimensions (e.g. `{2, 2, 3}` for two qubits + a qutrit).
- `qb["MatrixDimensions"]` = `{ElementDimension, Dimension}` — shape of the matrix when flattened.
- `qb["Qudits"]` = `Length[qb["Dimensions"]]`. Excludes 1-dim qudits.

### Operations (top-level, in `PacletInfo.wl`)
- `QuantumTensorProduct` — combine objects
- `QuantumPartialTrace` — reduce by tracing out indices
- `QuantumMeasurement` — non-demolishing measurement result wrapper
- `QuantumMeasurementSimulation` — repeated-shot simulation
- `QuantumStateEstimate` — tomography from `<|qmo → result, ...|>`
- `QuantumEvolve` — Hamiltonian/Liouvillian time-evolution (8 forms; Lindblad jumps with damping rates; passes `NDSolve`/`DSolve` opts)
- `QuantumEntangledQ`, `QuantumEntanglementMonotone` — entanglement detection/quantification
- `QuantumDistance`, `QuantumSimilarity` — state metrics *(QuantumSimilarity is undocumented — see `mistakes.md`)*
- `QuantumWignerTransform`, `QuantumPhaseSpaceTransform` — phase-space representations

### Secondary symbols (in `Usage.m` but NOT in `PacletInfo` Symbols list)
- `QuantumCircuitMultiwayGraph` — multiway graph of a circuit
- `QuantumMPS` — Matrix-Product-State / Operator → `QuantumCircuitOperator`
- `QuantumShortcut` — shorthand
- `QuantumWignerMICTransform` — minimal informationally complete basis

### "Operator-like" predicate
`QuantumFrameworkOperatorQ` (PackageScope, `QuantumFramework.m:10-13`) returns `True` for any of: `QuantumOperator`, `QuantumMeasurementOperator`, `QuantumChannel`, `QuantumCircuitOperator`. **`QuantumState` is NOT considered an "operator"** despite being callable.

## 2. Property dispatch

*(Detail per object pending — read from `*/Properties.m` files.)*

QF property keys take two forms:
- A bare string: `qs["Purity"]`
- A list with arguments: `qs["VonNeumannEntropy", base]`

### Structural getters short-circuit dispatch (`dfc741dc`/`f9dc1cdc`)

The $O(1)$ "structural" accessors that read a field straight off the wrapped expression now have **direct rules and bypass the cache machinery**. Each wrapped head has a one-line `"Basis"` getter (`QuantumOperatorProp[QuantumOperator[state_, _], "Basis"] := state["Basis"]`, and analogous getters on `QuantumChannel`/`QuantumHamiltonianOperator`/`QuantumMeasurementOperator`/`QuantumMeasurement` delegating to the wrapped operator; `QuantumState` already had one), and the cache-eligibility guard's exclusion list was widened to cover them (`QuantumOperator`: + `Order`/`InputOrder`/`OutputOrder`; `QuantumState`: + `State`; wrapper heads: + `Basis`; `QuditBasis`: stops caching `Representations`/`Properties`). This matters because the eligibility guard itself reads `obj["Basis"]["ParameterArity"]` on every access; before the direct getter, that resolved through ~124 clauses into the `AllProperties`-intersection delegation (§2 below), taxing every property read. The cache **store** is now `cacheProperty` (§5), not `Quiet[... = result, Rule::rhs]`.

### `"Properties"` vs `"AllProperties"` (changed across the wrapped heads, post-`ac32bda8`)

Every head that wraps another object (`QuantumOperator` over a `QuantumState`;
`QuantumChannel`/`QuantumHamiltonianOperator`/`QuantumMeasurementOperator` over a
`QuantumOperator`; `QuantumMeasurement` over a `QuantumMeasurementOperator`;
`QuantumCircuitOperator` over its last operator) now splits its property catalog in two:

- **`obj["Properties"]`** returns the head's **own** static list `$Quantum…Properties` only. It no
  longer recurses into the wrapped object, so it is O(1) (a stored list) and side-effect-free.
- **`obj["AllProperties"]`** is the old `"Properties"` behavior: `Union @ Join[ownList,
  wrapped["AllProperties"]]`, the full dynamically-computed set of everything the object answers
  to (this is what the delegated accessors actually dispatch against).

`obj["Properties"]` is a **strict subset** of `obj["AllProperties"]` (live: a `QuantumOperator`
gives 67 vs 209). The accessor-delegation conditions throughout (`MemberQ[Intersection[wrapped[…],
self[…]], prop]`) were all repointed from `"Properties"` to `"AllProperties"`, so delegation to the
wrapped object is unchanged. Per-head definitions: `QuantumState/Properties.m:64,69`,
`QuantumOperator/Properties.m:44,46`, `QuantumChannel/Properties.m:37,39`,
`QuantumHamiltonianOperator/Properties.m:24,26`, `QuantumMeasurementOperator/Properties.m:52,54`,
`QuantumMeasurement.m:94,96`, `QuantumCircuitOperator/Properties.m:17,19`. `"AllProperties"` is
added to every head's cache-prevention list (it must recompute each call). The behavioral catch
(callers that relied on `"Properties"` returning the full union) is logged in `mistakes.md`.

Internal predicates (`Utilities.m:76-78`):
```
propQ[prop_]    := MatchQ[prop, _String | {_String, ___}]
propName[prop_] := Replace[prop, name_String | {name_String, ___} :> name]
```

### State / order validation predicates (`Utilities.m:80-90`)

| Predicate | Accepts | Implication |
|---|---|---|
| `stateQ` | `VectorQ` (pure) **or** `SquareMatrixQ` (density) | States are dual-form. |
| `orderQ` | duplicate-free integer vector | |
| `autoOrderQ` | `orderQ`, `Automatic`, **or** `{out, in}` 2-tuple | Operators can have asymmetric input/output orders. |
| `targetQ` / `targetsQ` | positive-integer vector / list of those | |
| `measurementReprQ` | tensor of rank 2 or 3 | |
| `SymbolicQ` | any free symbol other than `Rational`/`Complex`/`NumericFunction`/`Constant` | Used internally to gate symbolic vs numeric paths. |

## 3. Computation paths

*(Hot-path traces pending — extract from `*/Properties.m`.)*

### Internal math primitives (from `Utilities.m`)

| Helper | Definition | Notes |
|---|---|---|
| `eigensystem` | `Utilities.m:151` — `Eigensystem` + optional `Chop`/`Sort`/`Normalize` + `Simplify` | **Generic** (no Hermitian specialization). Falls back to plain `Eigensystem` on `Eigensystem::eivec0`. |
| `kroneckerProduct` | `Utilities.m:126` — `KroneckerProduct` with scalar fallback | Smart: collapses to `Times` for scalars. |
| `MatrixPartialTrace` | `Utilities.m:132-137` — `ArrayReshape` + `TensorContract` | **Dense path only.** No sparse short-circuit. |
| `pauliMatrix`, `spinMatrix`, `fanoMatrix` | `Utilities.m:193-239` | Return `SparseArray`. Preserve sparsity through Pauli chains. |
| `identityMatrix` | `Utilities.m:120` — `IdentityMatrix[n, SparseArray]` | Sparse default. Edge-case: `0` or `{_,0}` → `{{}}`. |
| `MatrixInverse` | `Utilities.m:303` — `Inverse` then `PseudoInverse` fallback | **Quietly** falls back if `Inverse` fails — caller doesn't know precision was lost. |
| `SparseArrayFlatten` | `Utilities.m:275` — custom flatten preserving sparsity | Only triggers when `Length[dims] > 11`; smaller tensors use plain `Flatten`. |
| `matrixFunction` | `Utilities.m:306` — dispatches by `f` | `Plus/Minus/Times/Conjugate` direct; `Power` → `MatrixPower`; else `ResourceFunction["ComputeMatrixFunction"]`. |
| `GellMannMatrices`, `RegularSimplex` | `Utilities.m:254, 265` | SU(N) Gell-Mann basis; vertex-of-simplex utility. |
| `GramMatrix`, `GramDual` | `Utilities.m:268-270` | Dual-frame basis tools. |

## 4. Bottleneck catalog *(seeded; expand from `*/Properties.m`)*

- **`MatrixPartialTrace` is always dense.** For an N-qubit density matrix, materializes 2^N × 2^N. No sparse path even if input is sparse.
- **`eigensystem` is generic.** No Hermitian-specialized path; symbolic mode applies `Simplify` (can be slow). For large Hermitian matrices, `Eigensystem[matrix, Method -> "Hermitian"]` would be faster.
- **`MatrixInverse` silent fallback.** Failed `Inverse` → `PseudoInverse` quietly. Non-invertible operators give pseudo-inverse without warning.
- **`SparseArrayFlatten` threshold = 11 dims.** Below that, sparsity is lost on flatten.

## 5. Caching & profiling

- **`$QuantumFrameworkPropCache = True`** (`QuantumFramework.m:15`) — global flag. When `True`, property results are cached. **Process-wide; no automatic invalidation.**
- **Cache store = `cacheProperty`** (`Utilities.m`, PackageScope, `HoldFirst`; **new at `f9dc1cdc`**). Each head's property dispatcher (`QuantumOperatorProp`, `QuantumStateProp`, `QuantumBasisProp`, `QuditBasisProp`, and the wrapper-head `*Prop`s) commits a hit through `cacheProperty[HeadProp[obj, prop, args], result]` instead of the former `Quiet[HeadProp[...] = result, Rule::rhs]`. `cacheProperty` only executes the `Set` when the held key is `FreeQ[_, _Pattern | _Blank | _BlankSequence | _BlankNullSequence | _Optional]`, so a patterned key (symbolic/meta call) cannot create an over-matching DownValue and no `Rule::rhs` ever fires; such keys are simply left uncached. **The per-hit tax is gone:** direct `"Basis"` getters on every head + adding the $O(1)$ structural getters to the cache-exclusion lists mean an eligibility check no longer triggers the `AllProperties`-intersection fall-through. See §2 and `idioms-and-performance.md` (Caching).
- **`Memoize[f]`** (`Utilities.m`) — function-level memoization. **Rewritten `5b696969`:** stores results in a private `Association` keyed on `Hash[Hold[args]]` (not appended literal DownValues), with an optional `"CacheQ" -> pred` filter. Used by `Init.m` on `FromOperatorShorthand`; now also fit for the `*Prop` accessors. Clear by `Clear[f]`.
- **`$QuantumFrameworkProfile = False`** (`Utilities.m:362`) — set `True` for `EchoTiming` profile output via the `profile[label]` wrapper.

## 6. Performance map

*(Empirical table pending — to be calibrated with WL kernel benchmarks during synthesis.)*

Planned benchmarks:
- Partial trace: QF dense vs MPS-routed for 6/8/10/12-qubit density matrices.
- Hermitian eigensolve: QF generic vs `Method -> "Hermitian"` for d ∈ {16, 64, 256, 1024}.
- Property cache hit rate: warm vs cold access of `["VonNeumannEntropy"]`, `["Purity"]`, `["Eigenvalues"]`.
- Circuit composition: symbolic-symbolic vs numeric-numeric for N ∈ {5, 10, 15, 20} gates.
