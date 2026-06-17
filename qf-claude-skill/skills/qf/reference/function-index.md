# Function Index


> **`"Properties"` vs `"AllProperties"` (all wrapped heads, post-`ac32bda8`).** On `QuantumState`, `QuantumOperator`, `QuantumChannel`, `QuantumHamiltonianOperator`, `QuantumMeasurementOperator`, `QuantumMeasurement`, `QuantumCircuitOperator`, the bare `obj["Properties"]` now returns only the head's **own** static list (`$Quantum…Properties`); the **new** `obj["AllProperties"]` returns the full union including the wrapped object's catalog (the old `"Properties"` value). Use `"AllProperties"` to test whether an object answers to an inherited property. Per-head defs: `QuantumState/Properties.m:64,69` · `QuantumOperator/Properties.m:44,46` · `QuantumChannel/Properties.m:37,39` · `QuantumHamiltonianOperator/Properties.m:24,26` · `QuantumMeasurementOperator/Properties.m:52,54` · `QuantumMeasurement.m:94,96` · `QuantumCircuitOperator/Properties.m:17,19`. (Live: `QuantumOperator` = 67 `"Properties"` vs 209 `"AllProperties"`.)

Every public function: signature(s), kernel location, canonical example, pitfalls.

## Top-level objects (in `PacletInfo` `Symbols`)

### `QuditName[names]`
Convenient wrapper around qudit names with special formatting.
- Defined in `Kernel/QuditBasis/QuditName.m`.
- Signatures: `Usage.m:611-617`.
- Example: TODO — read `Documentation/.../Symbols/QuditName.nb`.

### `QuditBasis[names | dim]`
Single-qudit (or composite) basis.
- Forms (`Usage.m:591-607`):
  - `QuditBasis[<|name -> v, ...|>]` — element-by-element.
  - `QuditBasis[dim]` — N-dimensional canonical basis.
- Defined in `Kernel/QuditBasis/QuditBasis.m`. Named bases in `Kernel/QuditBasis/NamedBases.m`.

### `QuantumBasis[...]` — 5 forms
Full quantum basis with input/output dimensions and named bases (`Usage.m:7-82`):
- `QuantumBasis["name"]` — named.
- `QuantumBasis[<|name1 -> b1, ...|>]` — element-by-element.
- `QuantumBasis[{n1, n2, ...}]` — composite (list of dims).
- `QuantumBasis[n, m]` — n^m dimensional (m qudits, each n-dim).
- `QuantumBasis[{{n1, n2}, {m1, m2}}]` — explicit output × input dims.
- Defined in `Kernel/QuantumBasis/`.

### `QuantumState[qs[, qb]]` — 4 forms
Pure vector or density matrix in a basis (`Usage.m:508-555`):
- `QuantumState[qs, qb]` — vector/density matrix in a quantum basis.
- `QuantumState[qs]` — in computational basis (default).
- `QuantumState[<|name -> amp, ...|>, qb]` — by named amplitudes.
- `QuantumState["name"]` — named state.
- `QuantumState[QuantumState[..., qb1], qb2]` — basis change.
- Defined in `Kernel/QuantumState/`. Properties dispatcher in `Kernel/QuantumState/Properties.m` (1008 LOC — biggest property file).
- **`QuantumState::padded` (guard, `7e12e8f0`)**: `QuantumState[vector, dim]` warns when zero-padding rounds the qudit count up (`multiplicity > 1`, i.e. the length is not a power of the qudit dimension) — the `QuantumState[Table[0,{12}],2]`-is-a-4-qubit-norm-0-state trap. The single-qudit amplitude-prefix pad (`QuantumState[{1}, d]`) stays silent. The object is still produced; use `QuantumState["Register"[n]]` for a register.

### `QuantumOperator[...]` — 3 forms
Matrix/tensor + order + basis (`Usage.m:402-442`):
- `QuantumOperator[rep, order, qb]`
- `QuantumOperator["name", order, qb]`
- `QuantumOperator[qo, basis]` — basis change.
- Defined in `Kernel/QuantumOperator/QuantumOperator.m` (639 LOC).
- **`QuantumOperator::broadcast` (guard, `7e12e8f0`)**: when an operator is placed on an order longer than its qudit count, the multiplicity-broadcast rule (`QuantumOperator.m:~254`) tensors `LCM[...]` copies. It now refuses if `qo["Dimension"]^multiplicity > $QuantumOperatorBroadcastLimit` (PackageScope, `2^24`), emitting the message + `$Failed` instead of OOM-killing the kernel (the apply-path `{psi -> Range[12]}`-with-4-qubit-`psi` → `16^12`-amplitude trap). Legitimate small broadcasts (`"H" -> {1,2,3}`) are unaffected; a huge **named**-gate order (`QuantumOperator["H", hugeOrder]`) takes the Fourier-kronecker path, not this rule, so it still builds eagerly. Direct structural getter `QuantumOperatorProp[QuantumOperator[state_, _], "Basis"] := state["Basis"]` added (`dfc741dc`); see `architecture.md` §2.

**Decomposition properties** (`QuantumOperator/Properties.m`, added at `ac32bda8`):
- `qo["ZYZ"]` — 1-qubit Euler `U`-gate decomposition (with `GlobalPhase` when phase ≠ 0).
- `qo["KAK"]` — **two-qubit Cartan/KAK decomposition** (4×4 / `{2,2,2,2}` only). Numeric path `TwoQubitKAK` (magic-basis factorization → `e^{iφ}(a₁⊗a₂)·exp(i(c₁XX+c₂YY+c₃ZZ))·(b₁⊗b₂)`, emitted as `U`-gate + CNOT gadgets); structured-symbolic path `TwoQubitKAKSymbolic` (exact arithmetic, time-bounded `$kakSymbolicTimeBudget = 60` s, with a probe-substitution self-check that returns the operator undecomposed rather than a silently-wrong circuit, emitting `::kaksymbolic`). PackageScope: `TwoQubitKAK`, `TwoQubitKAKSymbolic`, `$MagicBasis`.
  - ⚠️ **Reconstructs up to a global phase only** (verified: `QuantumOperator["CNOT"]["KAK"]` matches CNOT up-to-phase but not exactly). See `mistakes.md`.
- `qo["Decompose", gateset:{"U","CX"}]` — umbrella: 1-qubit → `ZYZ`, 2-qubit → `KAK`. n ≥ 3 (QSD) is future work and falls through to the undefined-property path.

### `QuantumMeasurementOperator[...]`
Internal: `QuantumMeasurementOperator[QuantumOperator, targets]` (`QuantumMeasurementOperator/QuantumMeasurementOperator.m:9`) where `targets` matches `targetsQ` (list of `targetQ` integer-vector specs). `NHoldRest` attribute (`:24`).

Constructor forms (`Usage.m:313-365` + `QuantumMeasurementOperator.m:26-161`):
- `QuantumMeasurementOperator[matrix, target, basis]` — from a matrix + target + basis.
- `QuantumMeasurementOperator[basis, target]` — from a `QuantumBasis` (computational measurement).
- `QuantumMeasurementOperator[basis -> eigenvalues, target]` — basis with custom eigenvalues. Adds an `0` extra output order qudit for the eigenvalue register.
- `QuantumMeasurementOperator[QuantumOperator[...]["Diagonalize"], ...]` — from a diagonalized operator.
- `QuantumMeasurementOperator[QuantumOperator]` — defaults target to `InputOrder`.
- `QuantumMeasurementOperator[ops_List, target]` — from list of Kraus-like operators. Stacked via `StackQuantumOperators[Sqrt /@ ops]`.
- `QuantumMeasurementOperator[tensor /; TensorRank == 3, target]` — from POVM tensor. Splits by first index.
- `QuantumMeasurementOperator[Integer]` — singleton-list shorthand.

Composition (`:205-276`):
- `qmo[qs]` / `qmo[qm]` / `qmo[op]` — routes through `QuantumCircuitOperator`. Returns `QuantumMeasurement` if outcome shape matches.
- Direct `qmo[qmo']` is **commented out** (`:217-273`).

Properties (`Properties.m:5-19`, ~40):
- Targets: `Targets`, `Target`, `TargetCount`, `TargetIndex`, `TargetDimensions`, `TargetDimension`, `TargetBasis`.
- Eigen-side: `Eigenqudits`, `Eigendimensions`, `Eigenbasis`, `Eigenvalues`, `EigenvalueVectors`, `Eigenvectors`, `Eigenorder`, `Eigenindex`, `ExtraQudits`.
- State-side: `StateDimensions`, `StateQudits`, `StateBasis`.
- Type: `Type` (`Projection | POVM | Unknown`), `ProjectionQ`, `POVMQ`, `POVMElements`, `Operators`.
- Forms: `SuperOperator` (Stinespring dilation), `POVM`, `Canonical`, `CanonicalBasis`, `Sort`, `SortTarget`, `ReverseEigenQudits`, `DiscardExtraQudits`.
- `"Computational"` is now a **no-op for a non-computational eigenbasis** (`Properties.m:328`, `d7b8b3fd`): `qmo["Computational"] /; ! TrueQ[qmo["Eigenbasis"]["ComputationalQ"]] := qmo`. Rebasing an X/Y measurement into the computational frame stripped the eigenvalue-tagged pointer labels and silently corrupted the outcome statistics; for an already-computational basis it routes through the operator's `"Computational"` unchanged.

Type detection (`Properties.m:181-188`):
- `"Projection"`: no NonPositive output orders AND `OutputDimensions == InputDimensions`.
- `"POVM"`: ≥ 1 NonPositive output order (eigenvalue register).
- `"Unknown"`: neither.

Cache bypass: `{"Properties", "AllProperties", "QuantumOperator", "Operator", "SuperOperator", "DiscardExtraQudits"}` (`Properties.m:33`).

POVM normalization (`Properties.m:194`): `# / Mean[Diagonal[Total[#]]] & @ (# . # & /@ Through[Operators["MatrixRepresentation"]])` — non-standard normalization. See `mistakes.md`.

Named POVMs (`$QuantumMeasurementOperatorNames`, `NamedMeasurementOperators.m:7`):
- `M[basisArgs]` — generic measurement on a basis (default 2-dim).
- `RandomHermitian[basisArgs]` — random Hermitian via `(m + m†)/2` on `RandomComplex[1+I, ...]`.
- `WignerMICPOVM[basisArgs]` — Wigner-MIC POVM.
- `GellMannMICPOVM[d, s]` — d²-element Gell-Mann.
- `RandomPOVM[d, Method -> "Haar" | "Bloch"]` (`:44-56`) — random; default Haar.
- `TetrahedronSICPOVM[angles]` — 4-element rotated tetrahedron.
- `QBismSICPOVM[d]` — d²-element QBism.
- `HesseSICPOVM[]` — 9-element 3-dim.
- `HoggarSICPOVM[]` — 64-element 8-dim.

### ⚠️ NOT removed on `main` — `QuantumHamiltonianOperator` head and subdir
The whole `Kernel/QuantumHamiltonianOperator/` subdir (4 files, ~281 LOC) is **still present** on `main` (verified at HEAD `bfddf5ed`, 2026-05-27), including the typo filename `QuantumHamiltonialOperator.m`. The deprecation described in earlier passes of this audit was carried out only on the **unmerged** `qho-decapitation-shim` branch and never landed; there is no `deprecate/quantum-hamiltonian-operator` branch on this repo.

**The subsystem does not load.** Verified in a fresh kernel: `QuantumHamiltonianOperator` resolves to context `` Global` `` with **0 DownValues / 0 SubValues**. The constructor file `QuantumHamiltonialOperator.m:1` has a commented-out `Package[...]`, so its `PackageExport["QuantumHamiltonianOperator"]` never registers in the framework context. The head is therefore a non-functional orphan: `QuantumHamiltonianOperator[...]` stays unevaluated.

**`"Ising"` / `"TransverseIsing"` were NOT migrated to `QuantumOperator`.** They remain QHO-only (`QuantumHamiltonianOperator/NamedOperators.m:50, :86`, in the v2.0 `"Name"[args]` form). Verified: `QuantumOperator["Ising"[J], {1, 2}]` returns `Failure["InvalidName", ...]` (neither name is in `$QuantumOperatorNames`).

To get this functionality on `main` today, build it from `QuantumOperator` directly:
- Parameter wrapping → `QuantumOperator[H, "ParameterSpec" -> {t, t0, tf}]`.
- Lindblad packing → `QuantumOperator["Hamiltonian"[H, Ls, γs]]` (`NamedOperators.m:917`).
- Symbolic evolution operator → `qo["EvolutionOperator"]` (`QuantumOperator/Properties.m:584`, routes through `QuantumEvolve[qo, None]` / `DSolveValue`).
- Numeric evolution operator → `qo["NEvolutionOperator"]` (`:586-618`, uses `NDSolveValue`).
- Spectral gap → user-side `Differences[Sort[qo["Eigenvalues"][[1 ;; 2]]] /. param -> v][[1]]` and `FindMinimum[..., qo["ParameterSpec"][[1]]]`.
- Eigenvalue plot → `qo["EigenvaluePlot"]` (`QuantumOperator/Properties.m:620-623`).

Live QHO source (for reference; not loaded): `"Ising"`/`"TransverseIsing"` at `QuantumHamiltonianOperator/NamedOperators.m:50, :86`; `["Gap"]` defaults to a 100-point sweep at `Properties.m:56`; `["EvolutionOperator"]` uses symbolic `DSolve` at `Properties.m:86, :105`. `"LinearInterpolation"` and the `"EvolutionOperator"[U]` inverse also still live there.

### `QuantumChannel[...]`
Internal: `QuantumChannel[QuantumOperator]` in **Stinespring form** — the operator has more output qudits than input; the extra output qudits (encoded as non-positive orders: `0` and below) represent the environment to be traced out (`QuantumChannel/QuantumChannel.m:9-10`).

Constructor forms (`Usage.m:86-112` + `QuantumChannel.m:21-43`):
- `QuantumChannel[{op1, op2, ...}, args]` — Kraus form: list of `QuantumOperator`s. Stacked via `StackQuantumOperators` (`Utilities`).
- `QuantumChannel["name"[params], ...]` — named (see catalog below).
- `QuantumChannel[QuantumMeasurementOperator | QuantumMeasurement]` — converts via `["POVM"]`.
- `QuantumChannel[QuantumOperator, args]` — wraps an operator after sorting and reordering.

Named channels (`$QuantumChannelNames`, `NamedChannels.m:7-11`) — all 2-dim by default:
- `BitFlip[p]`, `PhaseFlip[p]`, `BitPhaseFlip[p]` — 2-Kraus probabilistic Pauli.
- `Depolarizing[p]` — 4-Kraus uniform.
- `AmplitudeDamping[γ]`, `PhaseDamping[λ]` — 2-Kraus.
- `GeneralizedAmplitudeDamping[γ, p]` — 4-Kraus thermal.
- `ResetError[p0, p1]` — 5-Kraus reset model.

Properties (`QuantumChannel/Properties.m:5-11`): channel-specific are `QuantumOperator`, `Operator`, `TraceOrder`, `TraceQudits`, `TraceDimensions`, `Trace`, `TracePreservingQ`. **Inherits all `QuantumOperator` properties**.

Notable properties:
- `qc["TracePreservingQ"]` — checks `Trace[E^†·E] == I` numerically; returns `Undefined` for symbolic.
- `qc["Adjoint"]` — Stinespring dual via index-permutation.
- `qc["EvolutionChannel"]` / `["NEvolutionChannel"]` — exponentiate each Kraus as if it were a Hamiltonian (`Properties.m:59-65`). See `mistakes.md`.
- `qc["DiscardExtraQudits"]` — collapses environment qudits explicitly.

### `QuantumCircuitOperator[{obj1, obj2, ...}]`
Internal: `QuantumCircuitOperator[Association[Elements, Label, ...CircuitDraw opts]]`. Elements are `QuantumFrameworkOperatorQ` objects OR Barriers (`QuantumCircuitOperator.m:11, :13-14`).

Constructor forms (`Usage.m:126-137` + `QuantumCircuitOperator.m:44-77`):
- `QuantumCircuitOperator[{op1, op2, ...}, opts]` — list of operators (or shorthands; parsed via `FromCircuitOperatorShorthand`).
- `QuantumCircuitOperator[op]` — single op → wrapped.
- `QuantumCircuitOperator[qco, opts]` — re-wrap with new opts.
- `QuantumCircuitOperator[{params...}]` — list of parameter sets → operators.
- `QuantumCircuitOperator[Graph]`: converted via `TensorNetworkQuantumCircuit[ToTensorNetworkGraph[g, Method -> "Random"]]` (`QuantumCircuitOperator.m:81`; repointed from `GraphTensorNetwork` at `4f2e2034`).

Application (`QuantumCircuitOperator.m:82-148`):
- `qco[qs, Method -> ...]` — apply to a state via `quantumCircuitApply`.
- `qco[opts]` — apply to default `|0...0⟩` register (`:129-142`).
- `qco[qm]` — apply to a measurement.
- `qco[op]` — prepend operator to elements.

Method dispatch table:
| Method | Engine | Constraint |
|---|---|---|
| `Automatic`, `"TensorNetwork"` | `TensorNetworkApply` (Wolfram/TensorNetworks paclet) | default |
| `"Schrodinger"` | `Fold` + `FullSimplify` per step + progress bar | symbolic-friendly but slow |
| `"QuEST"` | external QuESTLink | qubits only (`{2..}` dimensions) |
| `"Qiskit"` | Qiskit roundtrip | external bridge |
| `"Stabilizer"` | `PauliStabilizerApply` (`Kernel/Stabilizer/PauliStabilizer.m`) | Clifford/stabilizer circuits (see Stabilizer subsystem section) |

Sub-files (`Kernel/QuantumCircuitOperator/`):
- `QuantumCircuitOperator.m` (229 LOC) — constructor, application dispatch
- `Properties.m` (461 LOC) — ~50 properties; 17 bypass cache
- `TensorNetwork.m` (280 LOC) — `TensorNetworkApply`, `TensorNetworkCompile`, `QuantumTensorNetwork`, `QuantumTensorNetworkGraph`, `TensorNetworkQuantumCircuit` (inverse)
- `Draw.m` (988 LOC) — `CircuitDraw`, `drawGate`, `drawMeasurement`, `drawWires`, `drawBarrier`, `circuitDraw`, `expandLevel`. Public API: `CircuitDraw[qco, "Dynamic" -> True]` for click-to-expand interactive circuits. Sub-circuits are `EventHandler`-wrapped (`:824-843`) for `qc["ToggleExpand", path]` mutations. Wire thickness is dimension-aware: `Log2[dim/4 + 3/2]` (`:53`). Style maps `$GateDefaultBoundaryStyle` and `$GateDefaultBackgroundStyle` (`:13-35`) — pattern-based color assignment by label.
- `NamedCircuits.m` (795 LOC) — see catalog below.
- `Qiskit.m` (869 LOC) — `QiskitCircuit`, `QuantumCircuitOperatorToQiskit`, `ImportQPY`, `qiskitQASM`, `qiskitPrimitiveSubmit` (async SamplerV2/EstimatorV2), `qiskitReorderByQubit` (layout-aware bit decode), `qiskitNamedOperator`/`qiskitControlledOperator` (emit v2.0 `"Name"[args]`). `ImportQASMCircuit` **deprecated** → forwards to `QuantumQASM`. See QASM-hub section below.
- `QuantumQASM.m` (337 LOC) — **the single QASM hub**. `QuantumQASM` (import + export), `QiskitTarget` (validated device-spec carrier). `qasmEmit*` PackageScope workers. See QASM-hub section.
- `OpenQASMImport.m` (505 LOC) — native WL OpenQASM 2/3 importer (`ImportOpenQASM`, `qasmStringQ`, `qasmFileQ`). **No Python/qiskit.** `$qasmGateTable` (~40 gate specs). See QASM-hub section.
- `IBMQuantum.m` (474 LOC) — `IBMJobSubmit` + `IBMJob` (lossless result handle). `iMethodIBMJob` (blocking-Method constructor). See IBM-jobs section.
- `QuEST.m` (107 LOC) — `FromQuESTLink`, `ToQuESTLink`, `QuESTCompile`, `QuESTApply`. Qubits only. (`FromQuESTLink[name[args]]` now emits the v2.0 call form.)
- `Classiq.m` (42 LOC) — `ClassiqQuantumState[qs]`. Routes its returned QASM through `QuantumQASM` (was `ImportQASMCircuit`).
- `Formatting.m` (29 LOC) — summary box; >32 gates → fake placeholder icon.

Notable circuit properties:
- `qco["Gates"]` / `["GateCount"]` / `["OperatorCount"]` — flattened operator count.
- `qco["Depth"]` — circuit depth (max layer count per qudit).
- `qco["Width"]` — qudit-range span.
- `qco["Layers"]` — partition into commuting blocks; useful for parallelism.
- `qco["Topology"]` — graph of operator-label connectivity.
- `qco["TensorNetwork"]` / `["TensorNetworkGraph"]` — tensor network representation.
- `qco["Hypergraph"]` — via WolframInstitute/Hypergraph.
- `qco["ZXTensorNetwork"]` — ZX-calculus tensor network.
- `qco["QASM"]` — OpenQASM 3.0 export.
- `qco["Bend"]`, `qco["Double"]`, `qco["Dagger"]`, `qco["Conjugate"]`, `qco["Transpose"]` — superoperator/dual transforms.
- `qco["Flatten"[lvl]]` — recursively unwrap nested circuits.

Recent commit (`cbbe9368`): `"Parameters"` option simplification (each operator inherits `Parameters` from the circuit options).

### ⚠️ `QuantumHamiltonianOperator` — present but non-loading
- Subdir still on disk; the head does not load (`` Global` ``, 0 DownValues). See the QHO entry above for status and build-it-yourself targets.

## Stabilizer subsystem (`Kernel/Stabilizer/`, 6 new public symbols)

Added to `PacletInfo` `Symbols` at HEAD `ac32bda8`: `PauliStabilizer`, `StabilizerStateQ`, `StabilizerFrame`, `CliffordChannel`, `GraphState`, `LocalComplement`. The old 493-line `Kernel/PauliStabilizer.m` monolith was **deleted** and rebuilt as a 17-file `Kernel/Stabilizer/` subdir. All verified loading in a fresh 2.0.0 kernel.

### `PauliStabilizer[...]` — Aaronson-Gottesman tableau (`Stabilizer/PauliStabilizer.m`, `Constructors.m`)
Internal: `PauliStabilizer[<|"Signs" -> {...}, "Tableau" -> t|>]` where `t` is a `{2, n, 2n}` array (X-block, Z-block; rows = `n` destabilizers then `n` stabilizers). `Signs` may be symbolic (F₂ polynomials in `\[FormalS][k]` outcome symbols — FangYing23 SymPhase); `ConcretePauliStabilizerQ` is the strict `{-1,1}` predicate used where symbolic signs aren't supported. `StabilizerStateQ[ps]` ⟺ `PauliStabilizerQ[ps]`; also true for a rank-1 `StabilizerFrame`. Rejects bare `QuantumState` (detection needs 4ⁿ tomography).

Constructor forms (`Constructors.m`, `PauliStabilizer.m`):
- `PauliStabilizer[n_Integer]` — `|0…0⟩` register. **Closed-form since the June 2026 rework** (`Constructors.m`): assembles `Signs`/`Tableau` directly (destabilizers Xᵢ, stabilizers Zᵢ, all signs +1) instead of routing a known tableau through the generic validating constructor; n=1000 355 → 4.9 ms.
- `PauliStabilizer[{"XZZXI", ...}]` / `[{stabStrings}, {destabStrings}]` — signed Pauli strings (`:228, :238`).
- `PauliStabilizer[qs_QuantumState]` — **4ⁿ tomography** (`:89`): `<ψ|P|ψ>` over all Paulis, greedy F₂-rank generator pick, symplectic Gram-Schmidt destabilizers (AarGot04). Stores a `"GlobalPhase"` key so `ps["State"] === qs` exactly. Cost O(4ⁿ), practical n ≤ ~8.
- `PauliStabilizer[qo_QuantumOperator]` (`:168`) — AG decomposition via `QuantumOperatorTableau` (`Conversions.m:51`), with recovered `"GlobalPhase"`.
- `PauliStabilizer[qco_QuantumCircuitOperator]` (`:189`) — fold over normal operators.
- `PauliStabilizer["name"]` — named codes (`NamedCodes.m`, `$PauliStabilizerNames`): `5QubitCode(1)`, `SteaneCode`/`7QubitCode(1)`/`SteaneCode1`, `9QubitCode(1)`, `"Random"[n:5]` (Bravyi-Maslov/Koenig-Smolin Mallows sampler, `RandomClifford.m`). The `n`-th generator is logical Z-bar (states are `|0_L⟩`/`|1_L⟩`). **`"Random"[n]` hardened (June 2026):** the symplectic `Inverse` is taken over `\mathbb{F}_2` (`Modulus -> 2`, was over `\mathbb{Z}` where unit-triangular inverses grow exponentially) and intermediate matmuls stay on packed integers (was a structured-array route with tens-of-GB intermediates), fixing the n≈500 blow-up; the Mallows acceptance ratio uses exact rationals, fixing `General::munfl` at n≳538 (seeded output bit-identical). Viable to n≈538; past that still emits `General::munfl`.
- `PauliStabilizer["shortcut"]` — any `QuantumCircuitOperator` shortcut (string / rule / list).

Properties (`Properties.m`, ~30): `"Qubits"`/`"Qudits"`, `"GeneratorCount"`, `"Signs"`, `"Phase"`, `"X"`/`"Z"`/`"Tableau"`, `"Matrix"` (flat 2n×2n), `"Stabilizer(X/Z/Signs/Tableau)"`, `"Destabilizer(...)"`, `"PauliForm"`/`"Generators"`/`"Stabilizers"`/`"Destabilizers"`/`"PauliStrings"`, `"TableauForm"`, `"State"`/`"QuantumState"` (O(4ⁿ) materialize, `Conversions.m:69`), `"Circuit"`/`"QuantumCircuit"` (AG greedy synthesis, O(n³) gates, NOT length-minimized, `Conversions.m:91`), `"Operator"`/`"QuantumOperator"`.

**Performance representation (June 2026 rework, `Packed.m` + `Compiled.m`):** concrete stabilizers carry a *packed* form (62 tableau bits per machine word, chunk-major) gated by `psConcreteFastQ`; single gates, single-qubit + Pauli-string measurement, and formatting run directly on packed words, and `ps["Tableau"]` lazily reconstructs the canonical array. **`ps["ApplyCircuit", specs]`** (`Compiled.m`, also reached by `PauliStabilizerApply` / `qc[Method -> "Stabilizer"]`) runs a whole Clifford circuit through one `Compile`-to-C bulk fold over a generator-major packed tableau (~0.2-0.3 µs/gate, flat in n; auto-selects `"C"`/`"WVM"`); `specs` accepts `op -> q`, `op -> {q...}`, the `"C"[g -> t, c, _]` controlled form, and legacy `{op, q...}`; any unencodable gate routes to the per-gate fold (so non-Clifford P/T circuits still return a `StabilizerFrame`). **Storage-form trap:** single gates return packed objects, `ApplyCircuit` and constructors return canonical ones; `SameQ` across the two fails, so compare `["Tableau"]`/`["Stabilizers"]`. **CompiledFunction trap:** the bulk fold silently interprets (~100× slower, no warning) on unpacked args, so the kernel forces `Developer\`ToPackedArray` at every boundary.

Methods (gate updates, all return a new `PauliStabilizer`): `ps["H"|"S"|SuperDagger["S"]|"X"|"Y"|"Z"|"CNOT"|"CX"|"CZ"|"SWAP"|"V"|SuperDagger["V"], j (,k)]` (`GateUpdates.m`, packed kernels in `Packed.m` when `psConcreteFastQ`); `ps["ApplyCircuit", specs]` (compiled bulk fold, above); `ps["Permute"|"PermuteQudits", perm]`; `ps["Dagger"|"Inverse"]` (`GateUpdates.m:114`, emits `::singular` on a singular symplectic matrix); `ps["PadRight"|"PadLeft", n]`; composition `left[right]` via symplectic mod-2 product with `phaseLookup` corrections (`Compose.m`, O(n³); ~55× faster via packed reconstruction, and a crash on all-zero / ≤1-generator rows was fixed); `QuantumTensorProduct[a, b]` block-diagonal merge (`Compose.m:56`).
- **Non-Clifford boundary**: `ps["P"[θ], j]`, `ps["T", j]`, `ps[SuperDagger["T"], j]` return a **`StabilizerFrame`** (doubles the term count; GarMar15).
- **Measurement** (`PauliMeasure.m`, `Measurement.m`): `ps["M"|"Measure", q]` (single Z, AarGot04 destabilizer trick) → `<|outcome -> ps'|>` (1 key deterministic, `{0->,1->}` random); `ps["M", {q1,...}]` recursive; `ps["M", "XZZXI"]` joint Pauli-string measurement. `ps[q__Integer]` / `ps[{q...}]` / `ps[]` are `"M"` shorthands. **June 2026 rework:** both branches now operate on the packed representation (single-qubit n=500 16.5 → 0.19 ms; non-deterministic Pauli-string n=500 42.3 s → 59 ms). Two pre-existing AG bugs were fixed in passing: the row-sum update omitted the destabilizer rows (`1c524432`; corrupted *later* measurement correlations on some branches, invisible to single-measurement tests) and the Pauli-string branch left a non-abelian group (`52d3ea89`); plus a `packedMeasureZ` destabilizer-phase `Boole` mismatch (`730abcec`). Cross-validated against dense simulation + Stim; branch outputs unchanged for users.
- **Symbolic measurement** (`SymbolicMeasure.m`, FangYing23): `ps["SymbolicMeasure", q]` (fresh `\[FormalS][k]` symbol), `ps["SubstituteOutcomes", rules]`, `ps["SampleOutcomes", n]`. `$StabilizerSymbolCounter` is the session counter. **Limitation**: a deterministic measurement after a symbolic one is physically correct but its outcome polynomial is not stamped into the signs (see `mistakes.md`).
- **Inner product / expectation** (`InnerProduct.m`): `ps["InnerProduct", other, Method -> "Direct"|"ClosedForm"]` — Direct materializes vectors (O(2ⁿ), full complex phase); ClosedForm GarMarCro12 §3 is O(n³) but **magnitude + ±1 sign only** (complex phase TODO), requires concrete signs. `ps["Expectation", "XZZXI"]` closed-form AG i-power tracking → `±1` or `0`.
- **Entropy** (`Entropy.m`): `ps["Entropy", partition]` / `[q]` — Fattal-et-al closed form `rank_{F₂}(restricted gens) - |A|`, O(n³), no OOM at large n. Out-of-range → `::partition`.
- **UpValues** (`Conversions.m:142`, `HybridInterop.m`): `QuantumState[ps]`, `QuantumOperator[ps]`, `QuantumCircuitOperator[ps]`, `qo[ps]`, `qmo[ps]` (Pauli-string fast path else `::nonpaulibasis` fallback), `qc[ps]` (named Pauli-channel → probabilistic-mixture list).

### `StabilizerFrame[...]` (`Stabilizer/StabilizerFrame.m`)
Superpositions `Σ cᵢ|sᵢ⟩` of stabilizer states (GarMar15 Quipu). `StabilizerFrame[{{c, ps}, ...}]` or `[ps]`, internally `StabilizerFrame[<|"Components" -> {{c, ps}, ...}, "Paulis" -> ...|>]`. Properties: `"Components"`, `"Coefficients"`, `"Stabilizers"`, `"Length"`, `"Qubits"`, `"StateVector"`/`"State"` (materialize), `"InnerProduct"`. Clifford gates distribute over components; `"P"[θ]`/`"T"`/`SuperDagger["T"]` double the frame. `Plus`/`Times`/`Equal` upvalues defined.

**Relating Paulis + phase-coherent materialization (`61bc0e39`, line numbers in this file shifted with the fix, so cite handlers by rule name).** A **gate-built** frame additionally carries a `"Paulis"` key: one relating Pauli per component, encoded `{xbits, zbits, coeff}` with `coeff ∈ {1, -1, I, -I}` (so `{1,1}` on a qubit is Y), component 1's relating Pauli the identity, under the contract `|componentᵢ⟩ = matrix[Pᵢ] . |reference⟩`. Justification: every component of a gate-built frame shares one stabilizer tableau (a Z update leaves it unchanged, a Clifford update transforms it identically), so each component is a definite Pauli on a single reference component; the relative phases are physical but **not** recoverable from the final frame data alone, hence the bookkeeping. Maintained O(n)/gate by the `sfGate1`/`sfGate2`/`sfTimesZ` lookup tables: the Clifford-distribute rule (`f_StabilizerFrame[gate_String, args___]`) conjugates every relating Pauli `P -> G.P.G†` via `sfConjugatePauli`; the `"P"[θ]`/`"T"` doubling rule left-multiplies the new branch by `Z_q` via `sfTimesZPauli`; the `Times` upvalue carries `"Paulis"` through. `f["StateVector"]` (the `"StateVector"` handler) now **branches**: with a `"Paulis"` key it materializes the reference component **once** and relates the rest by `sfApplyPauli` (recovering all relative phases, phase-coherent up to a single global phase, same contract as a bare `PauliStabilizer`); without the key (hand-built `StabilizerFrame[{{c, ps}, ...}]`, or any `Plus`-combined frame, where components need not share a tableau) it keeps the prior independent per-component sum, so nothing regresses. **Verified live (working-tree paclet):** `T|+⟩ = {1, e^{iπ/4}}/√2` exactly (the case the old per-component readout swapped to fidelity 0.5); amplitude-exact vs dense `Method -> "Schrodinger"` across 216 random Clifford+T circuits (worst phase-aligned amplitude error 0); `⟨Y⟩` on `T|+⟩` is now `+0.7071` (was sign-flipped to `-0.7071`); manual frames unchanged. Cost ≈ +15% `ByteCount`, +8-12% build time per frame; the `2^k` component-doubling-per-non-Clifford-gate law and the lack of rank compression are unchanged.

### `CliffordChannel[...]` (`Stabilizer/CliffordChannel.m`, Yashin25)
CPTP map stabilizer→stabilizer as a Boolean tableau `[U_A | U_B | c]`. Constructors: `CliffordChannel[<|"UA","UB","c","InputQubits","OutputQubits"|>]`, `["Identity", n]`, `[ps]` (state-prep Choi), `[qc_QuantumChannel]` (deterministic single-Pauli only, else `::stochastic`). Accessors: `"Rank"`, `"Tableau"`, `"UA"`/`"UB"`/`"c"`, `"InputQubits"`/`"OutputQubits"`, `"Source"`. Composition `cc1[cc2]` via F₂ null-space with AG row-sum phase + |Φ⁺⟩ contraction-sign correction (the X·Z→Y phase that simple c-bit XOR misses). `cc[ps]` evolves a stabilizer state (identity / state-prep / matching-input cases; else `::stateevol`).

### `GraphState[...]`, `LocalComplement[...]` (`Stabilizer/GraphState.m`, AndBri05)
`GraphState[<|"Graph"->g, "VOPs"->{0..}|>]`; `GraphState[g_Graph]` (VOPs all-identity). `GraphState[ps_PauliStabilizer]` only for **graph-form** stabilizers (Xᵢ on diagonal, Z/I elsewhere) else `::nongraph`. Properties: `"Graph"`, `"VOPs"`, `"Vertices"`/`"Edges"`/`"VertexCount"`/`"EdgeCount"`, `"AdjacencyMatrix"`, `"Stabilizers"` (Kᵢ = Xᵢ ∏_{j∈N(i)} Zⱼ), `"Qubits"`, `"PauliStabilizer"`. `LocalComplement[g|gs, v]` toggles edges among v's neighbors (AndBri05 Def 1). **VOP tracking is stubbed** (only index 0 = identity; richer VOPs and VOP-tracked LC are on the roadmap).

## QASM hub: `QuantumQASM` + native importer (`QuantumQASM.m`, `OpenQASMImport.m`)

Secondary symbols (exported from the framework context, NOT in `PacletInfo` `Symbols`).

- **`QuantumQASM[...]`** — the single QASM entry point.
  - Import: `QuantumQASM[src_String]` (OpenQASM source, detected by `OPENQASM` header), `QuantumQASM[File[f]]`, `QuantumQASM["*.qasm"]` → `QuantumCircuitOperator`, **native WL, no Python** (`OpenQASMImport.m`). `QuantumCircuitOperator[qasmString | File[...] | "*.qasm"]` now routes here.
  - Export (default, native WL): `QuantumQASM[qco]` / `[qco, "WL"]`, `QuantumQASM[qo]` / `[qo, "Simple"]`. Qubit-only: a non-2-dim wire → `Failure["QuantumQASM", ...]` with `"NonQubitDimensions"`.
  - Export (qiskit path): any extra arg (a native gate set, `QiskitTarget`/`CouplingMap`, `"Version"`, transpiler option) routes through qiskit — faithful dump or transpile-to-basis. `"Provider"`/`"Backend"` engages the hardware-provider path (`qiskitQASM`). Options are CamelCase→snake_case normalized and validated against the live `transpile` signature (version-resilient). **(Shipped `f370f568`):** `qiskitQASM` (`Qiskit.m:639`) gains `"OptimizationLevel"` (Automatic → 1) and now transpiles against the backend's `Target` (per-instruction error + duration), so layout/routing are error-aware; `build_target` (`QuantumQASM.m:58`) applies `instruction_properties` via `update_instruction_properties` (recording skipped rows under a `qualink_skipped_props` metadata key, not silently dropping, for unknown gates/loci). The `"Target"` option is now `Switch`-validated (`QuantumQASM.m:381`): a `QiskitTarget` → its `"Spec"`, an `Association` → normalized keys, anything else (a failed `QiskitTarget[spec]`, `$Failed`, a stray value) → a `Failure` naming the bad head instead of an opaque Python `build_coupling_map` error. Still bridge-based, no native router.
- **`ImportOpenQASM`** (PackageScope) — the native parser. Supports: version header, include, qreg/creg + qubit/bit, standard + parametrized gates, custom `gate` defs (→ named subcircuit), inv/pow/ctrl/negctrl modifiers, measure (v2 + v3), reset, barrier, gphase. Anything else → categorized `Failure` (never a half-built circuit). `$qasmGateTable` (~40 gate name→spec entries). Import refinements since `ac32bda8` (`OpenQASMImport.m`): **physical qubits `$N`** resolve to absolute wire N+1 (`qasmResolveQ`, the form qiskit emits for transpiled circuits); **controlled gphase** (`inv/pow/ctrl/negctrl @ gphase`) imports as a multi-controlled phase (`qasmGPhase` hosts a 1-qubit phase gate on one control and controls it with the rest) instead of a bare `GlobalPhase`; **scientific-notation angles** (`1.5e-3`) are normalized to WL `*^` before parsing (`qasmAngle`); **`reset` at circuit start is dropped** to identity (a QF wire with no explicit initial state is already `|0⟩`; operands are still resolved so an out-of-range qubit is reported).
- **`QiskitTarget[spec_Association]`** — validated, normalized device-spec carrier; rebuilds the live qiskit `Target` from the spec at point-of-use (no cross-version pickle). Properties: `"Spec"`, `"NumQubits"`, `"CouplingMap"`, `"OperationNames"`, and (`f370f568`) `"InstructionProperties"` (`QuantumQASM.m:242`): the per-instruction error/duration rows `{{gate, {qubits...}, <|"Error"->_, "Duration"->_|>}, ...}` or `Missing["NotAvailable"]` when none were supplied; the summary box adds an "Error-aware:" item. The spec's `instruction_properties` key is **reserved** (`QuantumQASM.m:208`, `reserved = {'instruction_properties'}`): it is a `QiskitTarget`-only carrier, excluded from the `from_configuration` unknown-key check and applied afterward by `build_target` via `update_instruction_properties`. **(Shipped `f370f568`):** a new extractor `QiskitTarget["FromBackend" -> name, "Provider" -> p, opts]` (`Qiskit.m:671`) reads any qiskit-ecosystem backend's already-populated `Target` (per-instruction error + duration, coupling map) into a QF `QiskitTarget` through the shared `qiskitInitBackend` (so IBM / AWS Braket / the local credential-free `GenericBackendV2` fake device (provider `"GenericV2"`, `Qiskit.m:446`, `GenericBackendV2(num_qubits, seed=42)`, no credentials/network) all work via one code path); `operation_names` is filtered to real basis gates (via `get_standard_gate_name_mapping`), `measure` kept for the instruction-properties walk; the carried `InstructionProperties` make downstream transpiles error-aware. **Verified live (working-tree paclet):** `QiskitTarget[<|...,"instruction_properties"->{{"cx",{0,1},<|"Error"->0.01,"Duration"->3.5*^-7|>}}|>]` constructs and `["InstructionProperties"]` returns the rows; a spec with none returns `Missing["NotAvailable"]`.
- **`ImportQASMCircuit`** — **deprecated** (`Qiskit.m`): emits `::deprecated`, forwards to `QuantumQASM`.

The public `qco["QASM"]` / `qo["QASM"]` / `qo["SimpleQASM"]` / `qc["QASM"|"QASM2"|"QASM3"]` properties all route through this hub (`qasmEmitCircuit` in `QuantumCircuitOperator/Properties.m:477`, `qasmEmitOperator`/`qasmEmitSimple` in `QuantumOperator/Properties.m:960`). The circuit emitter threads a global classical-bit counter so successive measurements write `c[0], c[1], …`; an unserializable gate → `Failure` naming the gate and suggesting a native gate set. New PackageScope helpers `qasmQubitsQ`, `qasmNonQubitFailure` (`QuantumQASM.m`). Emit refinements since `ac32bda8`:
  - **Full-precision real literals** via `qasmReal` (`QuantumOperator/Properties.m:952`, `ToString[N[x], InputForm]` with `*^`→`e`), replacing 6-digit `NumberForm`.
  - **Single-qubit emit** (`qasmEmitSimple`, `:960`) guarded by `qo["SquareQ"]` and using `UnitaryAnglesWithPhase` (not `UnitaryAngles`, which mis-folds a complex `mat[[1,1]]` into the angles); the controlled-gate path emits the leftover phase as a controlled `gphase`.
  - **State injection** (`qasmStateInjectionQ`/`qasmStatePrepLines`, `QuantumCircuitOperator/Properties.m:427-461`): a no-input/≥1-output operator (ket, `Cup`) lowers to `reset` + prep gates; qubit count uses `qco["Max"]`, not `qco["Arity"]`.
  - **Non-computational measurement rebasing** (`qasmRebaseMeasurement`, `:464-475`): an X/Y measurement basis `B` lowers to `Inverse[B]` (via `QuantumShortcut`) + computational measurement; computational measurements and genuine POVMs are left untouched.

## IBM Quantum jobs (`IBMQuantum.m`)

Secondary symbols (framework context). Require an active `ServiceConnect["IBMQuantumPlatform"]` (auto-`Needs` at load now).

- **`IBMJobSubmit[qco, backend:Automatic, opts]`**: async submit through qiskit's own SamplerV2/EstimatorV2 (`qiskitPrimitiveSubmit`); `run()` is **not** awaited, so a handle returns immediately. Options: `"Primitive" -> "sampler"|"estimator"`, `"Shots" -> 4096`, `"Observable" -> "Z"` (Pauli string / list / `{pauli, coeff}` pairs), `"Wait" -> False`, `"PrimitiveOptions" -> <||>`. **(Shipped `f370f568`):** an `"OptimizationLevel" -> Automatic` option (Automatic → 1, `IBMQuantum.m:123`) forwarded to `qiskitPrimitiveSubmit`/`generate_preset_pass_manager`, transpiling against the backend's error-aware `Target`; users may raise to 2/3. Qubit-only (same guard as `QuantumQASM`). The `"PrimitiveOptions"` are **validated before transpile/submit** by `ibmValidatePrimitiveOptions` (PackageScope, `IBMQuantum.m:52`) against qiskit's own `SamplerOptions`/`EstimatorOptions` schema: a bad/mis-nested key fails in ms with a message naming the key, the primitive, the valid keys there, and a sampler/estimator cross-hint. Returns an `IBMJob`.
- **`IBMJob[<|...|>]`**: lossless result handle; carries raw service responses under `"Raw"`. Its accessors are split into two sets:
  - **`IBMJob["Properties"]`** (`$ibmProperties`, `IBMQuantum.m:514`): **local, network-free** accessors only, safe to bulk-map (`AssociationMap[job, job["Properties"]]`): `"ID"`, `"Backend"`, `"ServiceObject"`, `"Status"`, `"Primitive"`, `"UserID"`, `"ProgramID"`, `"Cost"`, `"EstimatedRunTime"`, `"SubmittedCircuit"`, `"Options"`, `"QuantumSeconds"`, `"Timestamps"`, `"Duration"`, `"Samples"`, `"MeasuredQubits"`, `"NumBits"`, `"Qubits"`, `"Shots"`, `"Measurement"`, `"Counts"`, `"Probabilities"`, `"ExpectationValue(s)"`, `"StandardErrors"`, `"ExecutionSpans"`, `"Raw"`, `"Properties"`, `"Actions"`.
  - **`IBMJob["Actions"]`** (`$ibmActions`, `IBMQuantum.m:528`): **need a live connection**, kept out of `"Properties"`: `"Refresh"` (re-query + cache decoded results), `"Cancel"`, `"ExecutedCircuit"`. Each still callable by name; the unknown-property message lists both sets. **`"ExecutedCircuit"` moved here from `"Properties"`** (`a44c1eea`).
  - Sampler counts are decoded little-endian and reordered into ascending-qubit order via `qiskitReorderByQubit` + the captured `"MeasuredQubits"` (classical-bit to original-qubit) map.
- **`iMethodIBMJob`** (PackageScope): constructor used by the blocking `Method -> "Qiskit"` path; stores a precomputed `"Measurement"` + `"Counts"`.
- **Native sampler-result decoding** (`ibmSamples`, `IBMQuantum.m:357`, with helpers `ibmUnwrap`/`ibmNpyDataBytes`/`ibmBitArraySamples`): the completed-job path now decodes the RuntimeEncoder-serialized `PrimitiveResult` envelope directly in WL. Unwrap `<|"__type__","__value__"|>` layers, inflate the zlib'd `.npy` (`Developer`RawUncompress` + header skip), unpack each register's `BitArray`, and bit-pack per-register outcomes into one integer per shot. Replaces the old `Interpreter["HexInteger"]` hex path.
- **`IBMJob[...]["ExecutedCircuit"]`**: server-transpiled circuit via presigned QPY URL (now an **Action**, see above). Distinguishes HTTP 404 (no server-side transpile happened: already-ISA circuits to V2 primitives store none) from an unreachable `s3.direct` VPC-private endpoint (see `mistakes.md`).

## Operations

### `QuantumTensorProduct[objects]`
Tensor product of QF objects (`Usage.m:559-569`).
- Defined in `Kernel/QuantumTensorProduct.m`.

### `QuantumPartialTrace[qs | qb, s]` — 2 forms
(`Usage.m:446-467`):
- `QuantumPartialTrace[qs, s]` — trace out qubits indexed by `s`.
- `QuantumPartialTrace[qb, s]` — same on a basis.
- Defined in `Kernel/QuantumPartialTrace.m`.
- Uses `MatrixPartialTrace` (`Utilities.m:132`). **Dense-only path** — see `idioms-and-performance.md`.

### `QuantumMeasurement[...]`
Non-demolishing measurement result wrapper (`Usage.m:302-309`).
- Defined in `Kernel/QuantumMeasurement.m`.

### `QuantumMeasurementSimulation[state, qmo, counts]`
Repeated-shot measurement simulation (`Usage.m:369-383`).

### `QuantumStateEstimate[<|qmo1 -> result1, ...|>]`
Tomography from measurement outcomes (`Usage.m:488-504`).
- Defined in `Kernel/QuantumStateEstimation.m`.

### `QuantumEvolve[op[, ...]]` — 8 forms (heaviest dispatch)
Symbolic and numeric Hamiltonian/Liouvillian time-evolution (`Usage.m:196-298`):
- `QuantumEvolve[op]` — symbolic from initial register state, with `t` as time parameter.
- `QuantumEvolve[op, qs, t]` — symbolic from `qs`.
- `QuantumEvolve[op, qs, {t, ti, tf}]` — numeric over interval.
- `QuantumEvolve[op, None, ...]` — evolution operator only.
- `QuantumEvolve[op, {L1, L2, ...}, ...]` — Lindblad jump operators.
- `QuantumEvolve[op, {L1, ...} -> {γ1, γ2, ...}, ...]` — with damping rates.
- `QuantumEvolve[..., "AdditionalEquations" -> spec]` — piecewise/hybrid systems.
- `QuantumEvolve[..., "ReturnEquations" -> True]` — returns DEs without solving.
- Plus `NDSolve` / `DSolve` options pass-through.
- Defined in `Kernel/QuantumEvolution.m`.

### `QuantumEntangledQ[qs, biPartition, method:"Realignment"]`
True/False for entanglement on bipartition (`Usage.m:156-173`, `QuantumEntanglement.m:14-15`). **Default method is `"Realignment"`**, NOT Concurrence.

### `QuantumEntanglementMonotone[qs, biPartition, method:"Concurrence"]`
9 measures (`$QuantumEntanglementMonotones`, `:9-12`):
| Measure | Path | File:line | Constraints |
|---|---|---|---|
| `"Concurrence"` | vector: closed form; matrix: O(d₁²·d₂²) Y-operator iteration | `:45-51` | works generally, expensive for mixed |
| `"ConcurrenceVector"` | full vector | `:43` | |
| `"Negativity"` | `(\|ρ^T_B\|_1 - 1)/2` | `:54-55` | **2-qubit only** |
| `"LogNegativity"` | `Log_2 \|ρ^T_B\|_1` | `:58-59` | **2-qubit only** |
| `"EntanglementEntropy"` | Schmidt for pure; partial-trace VonNeumann for mixed | `:62-69` | |
| `"RenyiEntropy"` ≡ `"RenyiEntanglementEntropy"` | `(1/(1-α)) Log_2 Tr[ρ_A^α]` | `:74-83` | default α = 1/2 (`:71-72`) |
| `"Realignment"` | `Σ σᵢ - 1` (CCNR criterion) | `:85-88` | general, default for entangled-Q |
| `"MutualInformationI"` | `S(ρ_A) + S(ρ_B) - S(ρ)` | `:93-97, :109-110` | |
| `"MutualInformationJ"` | `S(ρ_B) - Σ pᵢ S(ρ_B\|i)` | `:99-104, :112-113` | requires explicit `QuantumMeasurementOperator` |
| `"Discord"` | `MutualInformationI - MutualInformationJ` | `:106-107, :115-116` | requires QMO |

### `QuantumPartialTrace[obj, qudits | {{out,in}, ...}]` — defined in `Kernel/QuantumPartialTrace.m` (98 LOC)
- `QuantumPartialTrace[qb_QuditBasis, qudits]` (`:9`): `qb["Delete", qudits]`.
- `QuantumPartialTrace[qb_QuantumBasis, {{out, in}, ...}]` (`:12-19`): traces matched output/input pairs.
- `QuantumPartialTrace[qs_QuantumState, qudits]` (`:28-34`): **uses `MatrixPartialTrace` (dense path)**.
- `QuantumPartialTrace[qs_QuantumState, {{out, in}, ...}]` (`:36-46`): **uses `TensorContract` directly (faster)**.
- `QuantumPartialTrace[qo_QuantumOperator, qudits]` (`:51-62`): order-mapped trace.
- **`QuantumPartialTrace[qco_QuantumCircuitOperator, qudits]`** (`:72-87`): builds a Cup-Cap categorical circuit.
- `QuantumPartialTrace[qm_QuantumMeasurement, qudits]` (`:67-69`): wraps state-level trace.
- Default qudits: full trace (`:48, :64, :89`).

### `QuantumWignerTransform[obj, opts]`, `QuantumWeylTransform[obj]` — defined in `Kernel/QuantumWignerTransform.m` (98 LOC)
- `WignerBasis[qb, opts]` (`:12-35`): builds the Wigner basis. Even-dim has 2× prefactor (`:32`); odd vs even-dim use different label patterns.
- Options: `"Exact" -> True`, `"Decompose" -> True`, `"Basis" -> "Pauli"`.
- `QuantumWeylTransform[obj]` (`:73-97`): **inverse** transform (PhaseSpace → Schrödinger). Requires `Sqrt[Dimensions]` integer.
- Both dispatch over basis/state/operator/channel/measurement/circuit.

### `QuantumPhaseSpaceTransform[obj, basis]` — defined in `Kernel/QuantumPhaseSpaceTransform.m` (94 LOC)
- Default basis is `Wigner` (`:12, :15`).
- Companion: `QuantumPositiveTransform[obj]` (`:74-93`) — splits quasi-probabilities into positive/negative parts via `Ramp`. Useful for displaying Wigner functions as actual probabilities.
- Dispatches over basis/state/operator/channel/measurement/circuit.

### `QuantumTensorProduct[objects...]` — defined in `Kernel/QuantumTensorProduct.m` (137 LOC)
Type-aware dispatch:
- `QuditBasis × QuditBasis` (`:16-20`): increments qudit index of second by `qb1["FullQudits"]`.
- `QuantumBasis × QuantumBasis` (`:25-43`): merges Output/Input/Picture (PhaseSpace wins)/Label.
- `QuantumState × QuantumState` (`:48-70`): mixed-state path uses `Double`/`Undouble`; pure path uses `KroneckerProduct`.
- `QuantumOperator × QuantumOperator` (`:75-88`): order-shifts on overlap.
- `QuantumChannel × QuantumChannel` (`:90-99`): preserves trace-qudit (NonPositive) order structure.
- `QuantumMeasurementOperator × QuantumMeasurementOperator` (`:102-120`): asymmetric POVM/projection auto-promotes to all-POVM (`:102-103`).
- Cross-type fallbacks (`:126-128`).
- `QuantumCircuit × QuantumCircuit` (`:133-134`): `Join[ops1, shift(ops2)]`.
- Final fallback returns `Failure[QuantumTensorProduct, "Undefined"]` (`:137`).

### `QuantumPartialTranspose[obj, args]` — defined in `Kernel/QuatumPartialTranspose.m` (filename typo, 7 LOC)
**One-line wrapper**: `QuantumPartialTranspose[qobj_, args___] /; QuantumStateQ[qobj] || QuantumOperatorQ[qobj] := qobj["Transpose", args]` (`:6`).

### `QuantumEvolve[op, lindblad, observable, ...]` — defined in `Kernel/QuantumEvolution.m` (267 LOC)
Top-level evolution dispatcher.
- 8 forms (`Usage.m:196-298`).
- Method branches on `defaultParameter` shape (`:50-58`):
  - `{symbol, t0, t1}` (numeric range) → `NDSolveValue` path with `Progress\`EvaluateWithProgress`.
  - bare symbol or `Automatic` (symbolic) → `DSolveValue` path.
- Lindblad form at `:91`: `dρ/dt = matrix.s - s.matrix + Σ(L.s.L† - ½(L†.L.s + s.L†.L))`.
- **PhaseSpace special path** (`:46, :63-69`) — Hamiltonian becomes a `HamiltonianTransitionRate`.
- `TrigToExp` applied to Hamiltonian and Lindblad operators (`:62, :74`) to assist symbolic simplification.
- **`MergeInterpolatingFunctions`** (`:213-249`) — composes prior numeric solutions into linearly-interpolated `SparseInterpolationFunction` for further evolution.
- Companion exports: `HamiltonianTransitionRate`, `LindbladTransitionRates` (`:264-266`).

### `QuantumMeasurement[qmo]` — defined in `Kernel/QuantumMeasurement.m` (367 LOC)
Wrapper `QuantumMeasurement[QuantumMeasurementOperator]` (`:12`). Distinct from `QuantumMeasurementOperator` — represents the outcome of a measurement.

Constructor forms (`:22-51`):
- `QuantumMeasurement[qs, target]` — wraps via `QuantumMeasurementOperator[qs["Operator"], target]`.
- `QuantumMeasurement[<|outcome -> probability, ...|>]` — multi-outcome measurement, default basis is computational with `Log2[Length]` qubits.
- `QuantumMeasurement[<|...|>, {qs1, qs2, ...}]` — outcomes paired with explicit states.

14 properties (`$QuantumMeasurementProperties`, `:56-66`):
- Outcomes: `Outcomes`, `Probability`, `Probabilities`, `ProbabilityList`, `ProbabilitiesList`.
- Distributions: `Distribution`, `NDistribution`, `MultivariateDistribution`, `NMultivariateDistribution`.
- States: `States`, `MixedStates`, `StateAssociation`, `StateAmplitudes`, `StateProbabilities`.
- Stats: `Mean`, `MeanState`, `Entropy`. **`"Entropy"` is now a closed-form Shannon entropy** (`:235-245`, `da372b3f`), mirroring `QuantumState["VonNeumannEntropy"]`: `-Σ p Log[b, p]` over the nonzero `Normalize`d outcome probabilities, `Chop @ Simplify`, `Enclose` → `Indeterminate` fallback, capped at `< 2^10` outcomes. Bare `qm["Entropy"]` returns `Quantity[_, "Bits"]` (base 2); the `logBase` argument `qm["Entropy", b]` and list-form `qm[{"Entropy", b}]` (which `QuantumState` already had) return a bare number in that base. Previously routed through `Information[CategoricalDistribution[outcomes, probs]] / Log[2]`, which fell back to Monte-Carlo sampling and spammed `RandomVariate::unsdst` + `Extract::psl1` on any symbolic/parametric outcome probabilities (and on plain display, via the summary-box `N @ qm["Entropy"]`).
- Eigen: `Eigenvalues`, `EigenvalueVectors`, `Eigenvectors`, `Projectors`, `Eigenstate`.
- Visualization: `ProbabilityChart`, `ProbabilitiesChart`.
- Simulation: `SimulatedMeasurement`, `SimulatedCounts`, `SimulatedStateMeasurement`.
- Cache excludes `Simulated*` properties (`:80`).
- `MixedStates` branches on `LabelHead` (`:129-141`): `Computational | Automatic` → computational form; `Eigen` → eigenbasis; else canonical.
- `SimulatedCounts[n_:100]` (`:253`) uses `MultinomialDistribution` — much faster than n iterations.
- **Symbolic-outcome guard (`fc430744`)**: `numericProbabilitiesQ` + the `HoldRest` `whenNumericProbabilities` wrapper (`:182-185`, `Enclose` + `ConfirmAssert`, no `Quiet`) gate every property whose ordering/Monte-Carlo path needs numeric probabilities (`"DistributionInformation"`, `"TopProbabilities"`, `"SimulatedMeasurement"`, `"SimulatedCounts"`, `"SimulatedStateMeasurement"`, `"TopStateProbabilities"`), returning `Indeterminate` (mirroring `QuantumState["VonNeumannEntropy"]`) instead of leaking `RandomVariate::unsdst`/`MultinomialDistribution::vprobprm`/`CategoricalDistribution::elmntavsl`/`Extract::psl1`/`KeyMap::invak`. `"Categories"`/`"Probabilities"`/`"ProbabilityTable"`/`"ProbabilityArray"` route straight to the distribution and **keep** their symbolic values; `"Distribution"`/`"NDistribution"` still return a symbolic `CategoricalDistribution`.

### `QuantumStateEstimate[measurementResult, opts]` — defined in `Kernel/QuantumStateEstimation.m` (1072 LOC, by Mihai Vidrighin)
Sophisticated quantum state tomography. Pipeline:
1. **Linear inversion** (`:67-82`): `PseudoInverse` of stacked POVM elements.
2. **Physical projection** (`:285-286`): pure-state via `EigenPrune` or general via `Physical` property.
3. **Maximum likelihood** (`:215-253`): adaptive step size, default 300 iterations (`:236`).
4. **Bayesian sampling** (`:294-296`): 150-step burn-in + 150-step adaptive on Bures manifold.

Companion: `QuantumMeasurementSimulation[state, qmos, n]` (`:59-60`) — sample from each QMO via `MultinomialDistribution`.

Returns `QuantumStateEstimation[<|14 keys|>]` including:
- `InversionState`, `MaximumLikelihoodState` (both `QuantumState`).
- `Invertible`, `PhysicalInversion` (booleans).
- `MaximumLikelihoodAcceptanceRatio`, `BayesianAcceptanceRatio`.
- `MaximumLikelihoodConvergence`, `MaximumLikelihoodConvergenceStepSize`.
- `Dimension`, `MeasurementOutcomes`, `TotalCounts`, `CountsPerMeasurement`, `MetropolisStepSize`.
- `BayesianSampler` — callable: `bayesianSampler[n, "StepSize" -> ..., "Burn" -> 0]` returns `n` `QuantumState`s sampled from the posterior.

### `QuantumDistance[qs1, qs2, t:"Fidelity"]`
State distance with measure `t` (`Usage.m:141-152`). Defined in `Kernel/QuantumDistance.m` (80 LOC).

8 measures (`$QuantumDistances`, `:10`):
| Measure | Formula | File:line | Notes |
|---|---|---|---|
| `"Fidelity"` | `1 - Re[Tr[MatrixPower[ρ₁ρ₂, 1/2]]]` | `:15-16` | **Returns `1 - F`, NOT F.** Use `QuantumSimilarity` for F. |
| `"RelativeEntropy"` | `Tr[ρ₁(log ρ₁ - log ρ₂)] / Log[2]` (Quantity Bits) | `:18-34` | Uses `MatrixLog[..., Method -> "Jordan"]`. |
| `"RelativePurity"` | `1 - Tr[ρ₁ ρ₂]` | `:36-41` | |
| `"Trace"` | `Re[Tr[MatrixPower[(ρ₁-ρ₂)†(ρ₁-ρ₂), 1/2]]]/2` | `:44-48` | |
| `"Bures"` | `Sqrt[2 × Fidelity]` | `:50-51` | |
| `"BuresAngle"` | `Re[ArcCos[1 - Fidelity]]` | `:53-54` | |
| `"HilbertSchmidt"` | `Norm[ρ₁-ρ₂, "Frobenius"]` | `:56-57` | |
| `"Bloch"` | `Re[EuclideanDistance[bloch₁, bloch₂]]/2` | `:59-60` | **`Dimension == 2` only.** |

Recent change (`23d8c92b`): updated.

### `QuantumSimilarity[qs1, qs2, distance:"Fidelity"]`
Implementation: `Kernel/QuantumDistance.m:63-79`. Maps each distance measure to a similarity in `[0, 1]`:
- `"Fidelity" | "Trace" | "Bloch" | "RelativePurity"`: `1 - d`
- `"Bures" | "HilbertSchmidt"`: `1 - d/Sqrt[2]`
- `"BuresAngle"`: `1 - d/(Pi/2)`
- `"RelativeEntropy"`: `2^(-d)` (exponential decay for unbounded metric)
- Default measure: `"Fidelity"`.

**Note**: No `Usage` message in `Usage.m` and no doc page, but symbol IS exported in `PacletInfo.wl:42` and fully implemented. See `mistakes.md`.

### `QuantumWignerTransform[obj]`
Phase-space (Wigner) representation (`Usage.m:582-587`).
- Defined in `Kernel/QuantumWignerTransform.m`.

### `QuantumPhaseSpaceTransform[obj, basis]`
Generic phase-space transform (`Usage.m:471-476`).
- Defined in `Kernel/QuantumPhaseSpaceTransform.m`.

## Secondary symbols (in `Usage.m`, NOT in PacletInfo `Symbols`)

These need explicit context loading or are accessed via primary symbols' methods:

- `QuantumCircuitMultiwayGraph[qc]` (`Usage.m:116-122`) — `Kernel/Multiway.m`.
- `QuantumMPS[QuantumState[...]]` / `QuantumMPS[QuantumOperator[...]]` (`Usage.m:387-398`) — returns `QuantumCircuitOperator`. In `Kernel/QuantumCircuitOperator/TensorNetwork.m`.
- `QuantumShortcut[]` (`Usage.m:480-484`) — shorthand factory.
- `QuantumWignerMICTransform[obj]` (`Usage.m:573-578`) — minimal informationally complete (MIC) basis. In `Kernel/MIC.m`.
- `QuantumAdiabaticEvolve[Hb, Hp, opts]` / `QuantumAdiabaticEvolution[<|...|>]` (`Kernel/QuantumOptimization.m`, added `ac32bda8`) — adiabatic-evolution analysis. `QuantumAdiabaticEvolve` builds `H(s) = (1-s)Hb + s Hp` (or a custom `"Schedule"`), eigensolves under `0 ≤ s ≤ 1`, and returns a `QuantumAdiabaticEvolution` carrier with `"MinimalEnergyGap"`, `"MaxAdiabaticCoupling"`, `"AdiabaticTimeEstimate"` (= `"TimeScaling" · ξ/g²_min`, default scaling 10), `"Energies"`, `"Eigensystem"`, plus `"EnergySpectrumPlot"`/`"SpectralGapPlot"`/`"AdiabaticCouplingPlot"`/`"AdiabaticPathPlot"`. Emits `::gap` on zero spectral gap, `::num` on Root-object eigenvalues. Options: `"Parameter" -> \[FormalS]`, `"TimeScaling" -> 10`, `"Schedule" -> None`. **Like the QHO head, `PackageExport` from the `QuantumOptimization`` context but resolves to `` Global` `` on a bare `Needs["Wolfram`QuantumFramework`"]`** — not in `PacletInfo` `Symbols`.

## Named operator catalog (`$QuantumOperatorNames`)

`Kernel/QuantumOperator/NamedOperators.m:10-35` lists ~110 named operators dispatchable via `QuantumOperator["Name"]` (defaults) or `QuantumOperator["Name"[args, ...], order, opts]` (the v2.0 form since `fa0d0113`). The legacy `QuantumOperator[{"Name", args, ...}]` form is **removed** — passing it now returns `Failure["InvalidName", ...]`. List form is kept only for genuine sequences (tensor products, gate sequences inside a circuit). `"SimpleSpider"` and `"HeisenbergWeyl"` were added to this list in `bfddf5ed` (2026-05-27); they had definitions before but were not in `$QuantumOperatorNames`. Grouped (all line refs re-verified at HEAD `bfddf5ed`):

### Identity / structural
- `Identity`, `I` — identity (`:113-149`); `QuantumOperator["I"[dims]]` for arbitrary dimensions.
- `Permutation` — qudit permutation gate (`:635-653`); supports cycles or list form.
- `Curry`, `Uncurry` — basis-bending operators (`:657-659`).
- `Zero`, `O` — zero operator (`:153`).

### Fourier / Hadamard
- `Fourier`, `InverseFourier` (`:486-497`).
- `Hadamard`, `H` — generalized DFT (`:584`); equivalent to Fourier.

### Rotations
- `XRotation`, `YRotation`, `ZRotation`, `RX`, `RY`, `RZ` (`:156-172`); built from `Exp[-I angle/2 Pauli]` via `MatrixExp`.
- `R[angle, ops...]` — generic rotation about an arbitrary product of operators (`:180`). Powerful.
- `U`, `U3`, `U2`, `U1` — IBM/qiskit-style 3-, 2-, 1-parameter unitaries (`:190-198`; `U1` shares the `Phase` rule at `:206`).
- `Phase`, `P` — diagonal phase gate (`:206`).

### Phase / shift
- `PhaseShift` — `Phase` with quantum-circuit-style step (`:216`).
- `GlobalPhase` (`:222`).
- `Diagonal` — diagonal operator from a list (`:231-236`).
- `FlipSign` — phase-flips a specific basis state (`:250-270`).

### Pauli / NOT
- `X`, `Y`, `Z`, `PauliX`, `PauliY`, `PauliZ` — 2-dim Paulis; `{name, dim}` for higher-dim generalizations using shift/clock matrices (`:549-559`).
- `Shift` ≡ `PauliX`, `ShiftPhase` ≡ `PauliZ` (`:549, :559`).
- `NOT` — `X` with `"Label" -> "NOT"` (`:568`).
- `RootNOT`, `V`, `SX` — `Sqrt[X]` (`V`/`SX` at `:277`, `RootNOT` at `:571`).
- `0`, `1` — projector-like Paulis (`-Z`, `Z`) (`:564-566`).
- `S`, `T` — phase gates `Phase[Pi/2]`, `Phase[Pi/4]` (`:273-275`).

### Swap family
- `SWAP` — `Permutation` with `Cycles[{{1,2}}]` (`:516`).
- `RootSWAP` — `Sqrt[SWAP]` (`:519`).
- `Braid` — Bell-basis or sign-permutation braid (`:523`).
- `CSWAP`, `Fredkin`, `C0SWAP`, `Fredkin0` (`:595-597`).
- `SUM` — qudit SUM gate (`:530`).

### Controlled
- `C`, `Controlled` — generic controlled-`op`, with `control1` and optionally `control0` lists (`:318-451`). Many overloads.
- `C0`, `Controlled0` — controlled-on-0 (`:339-451`).
- `CX`, `CY`, `CZ`, `CH`, `CT`, `CS`, `CPHASE`, `CP` — qubit-specific controlled gates (`:285-315`).
- `CNOT`, `C0NOT` (`:280-282`).
- `Toffoli` — `Controlled` with two controls (`:591`).
- `Deutsch` — controlled `i * Rx(θ)` (`:786`).
- **Note**: control structure is encoded in the Label, not Order — see `mistakes.md`.

### Random
- `RandomUnitary` — `CircularUnitaryMatrixDistribution` for square dim, otherwise `RandomPure` state (`:601-621`).
- `RandomHermitian` — random unitary + its dagger (`:625-631`).

### ZX-calculus / categorical
- `Spider`, `SimpleSpider` — generic phase spider (`Spider` at `:691`, `SimpleSpider` at `:710`; shared no-basis/arg forms `:687-689`).
- `XSpider`, `YSpider`, `ZSpider` — Pauli-flavored spiders (`:674`).
- `WSpider` — W-state spider (`:728`).
- `Measure`, `Encode`, `Copy`, `Decohere` — categorical primitives via spiders (`:753-761`).
- `Marginal`, `Discard`, `Trace`, `Reset` — discarding/reset (`:763-769`).
- `Cup`, `Cap` — bend/unbend categorical wires (`:777-781`).

### Quantum control / multiplexing
- `Switch[a, b]` — quantum-controlled switch via CSWAP (`:796-810`).
- `Multiplexer`, `BlockDiagonal` — block-diagonal of multiple operators (`:465-483`).

### Angular momentum
- `WignerD` — Wigner D-functions for spin-`j` (`:843-845`).
- `JX`, `JY`, `JZ`, `AngularMomentumX`/`Y`/`Z` (`:849-853`).
- `J+`, `J-`, `JX+`, `JX-`, `JY+`, `JY-`, `JZ+`, `JZ-` — raising/lowering (`:855-857`).
- `+`, `-`, `Up`, `Down`, `I+`, `I-` — bosonic raising/lowering (`:859-861`).

### Heisenberg-Weyl
- `HeisenbergWeyl` — generalized displacement-like operator (`:814`). The dispatch typo (`QuantumOperartor`) was fixed in `c74f3baf`; the name was added to `$QuantumOperatorNames` in `bfddf5ed`, so `QuantumOperator["HeisenbergWeyl"[d]]` now resolves (verified).

### Operator transforms
- `Double` — calls `["Double"]` on the wrapped operator (`:864`).
- `Left`, `Right` — left/right multiplication by an operator (vectorized form) (`:866-873`).

### Open systems / time evolution
- `Hamiltonian[args]` — `I * Liouvillian` (`:927`).
- `Liouvillian[H, Ls, Gammas]` — full Lindblad superoperator (`:900-923`). Supports vector or matrix form for damping rates. Internal helpers: `HamiltonianMixedOperator` (`:881`), `LindbladMixedOperator` (`:886`). **Common-basis rebase (2026-06-11):** before building, `H` and every jump operator are rebased into the Hamiltonian's output basis (computational when `H === None`), the same convention as `QuantumEvolve`'s direct path (`QuantumEvolution.m:80`). Named operators carrying their own basis (`"J-"`, `"J+"`, `"JX"`, ...) enter with their operator meaning, and the superoperator (and states evolved through it) come back in the Hamiltonian's basis rather than in whichever jump basis `Plus` happened to keep. See `mistakes.md` "Named-basis jump operators".

### String shortcuts
- `QuantumOperator["IXYZ"]` — string composed of `{I, X, Y, Z, H, S, T, V, P}` decodes as tensor product (`:920`).
- `QuantumOperator["cnot"]` — auto-uppercases to `"CNOT"` (`:927-930`).
- `QuantumOperator[{op1 -> {1}, op2 -> {1, 2}}, opts]` — list of rules → `QuantumCircuitOperator` → single operator (`:940`).

### Failure forms
- `QuantumOperator["InvalidName"[...]]` (or any bare `"InvalidName"`) returns `Failure["InvalidName", ...]` (`:948-953`). Use `Enclose`/`Confirm` to handle. (Same dispatch added across all object heads in `fa0d0113`.)

### Operator-shorthand parser
- `FromOperatorShorthand` (`:56-109`) — central parser; recognized by `$ShorthandOperatorPattern` (`:48-54`). Memoized in `Init.m:7`. Supports rules, strings, symbols, dagger, named-op patterns, numeric-function applications.

## Named QuditBasis catalog (`$QuditBasisNames`)

`Kernel/QuditBasis/NamedBases.m:8-23` — ~36 named bases. `QuantumBasis["name", ...]` and `QuditBasis["name", ...]` both dispatch here.

### Computational / identity
- `Computational`, `Identity`, `I` — n-dim computational basis (`:40-46`).

### Pauli (single-qubit eigenbases)
- `PauliX` ≡ `X` — eigenbasis `{|+⟩, |−⟩}` of σ_x (`:67-83`).
- `PauliY` ≡ `Y` — eigenbasis of σ_y.
- `PauliZ` ≡ `Z` — eigenbasis of σ_z (= computational for d=2).
- **All three return EIGENBASES, not the matrix-form Pauli operators** — see `mistakes.md`.

### Pauli (4-element σ-basis)
- `Pauli` — `{σ_0, σ_1, σ_2, σ_3}` = `{I, X, Y, Z}` matrices as basis elements (`:140-147`). 4 elements in 2×2 matrix space. Used for operator decomposition / process tomography.

### Spin / angular momentum
- `JX`, `JY`, `JZ`, `J`, `JI` — spin-`j` eigenbases of `spinMatrix[i, 2j+1]` (`:85-104`). Default `j = 1/2`.

### Bell
- `Bell` — 4 Bell states (`:48-58`).

### Fourier / phase
- `Fourier` — DFT of an underlying basis (`:107-121`). Default 2-dim.
- `Schwinger` — clock+shift Schwinger basis for `d × d` (`:128-138`).

### Dirac
- `Dirac` — 16-element gamma-matrix basis (`:149-175`).

### Phase-space
- `Wigner` — Wigner basis for d-dim (`:180-185`).
- `WignerMIC` — Wigner with MIC normalization (`:188`).

### MIC (minimal informationally complete)
- `Bloch`, `BlochSphere` — Bloch-vector bases (`:203-211`).
- `GellMann` — d²-element Gell-Mann SU(d) basis (`:191-194`).
- `GellMannMIC`, `GellMannBloch`, `GellMannBlochMIC` (`:196-221`).
- `RandomMIC`, `RandomHaarMIC`, `RandomBlochMIC` — random POVMs (`:232-246`). **Skipped by `$QuditBasisCache`** (each call gives a fresh basis).

### Wootters / Feynman
- `Wootters`, `Feynman` — finite-phase-space bases (`:224-229`).

### SICs (Symmetric Informationally Complete POVMs)
- `Tetrahedron` (`:273-278`).
- `QBismSIC` — d²-element QBism SIC (`:280-283`).
- `HesseSIC` — Hesse SIC (`:285-290`).
- `HoggarSIC` — Hoggar lines SIC (`:292-296`).

### MUB (Mutually Unbiased Bases)
- `Ivanovic` — d-prime MUB construction (`:251-264`); generalizes to composite via prime factorization (`:266-268`).

### String shortcuts
- `QuditBasis["XYZI"]` — chain of single-char Pauli names (`:302-308`). **Only accepts `{I, X, Y, Z}`** (narrower than QuantumOperator's `{I, X, Y, Z, H, S, T, V, P}`).
- `QuditBasis["XYZ"[d]]` — same chain at dim `d` (`:306`).
- `QuditBasis["XBasis"]` — `"Basis"` suffix is auto-stripped (`:310-311`).

### Phase-space subset
`$QuditPhaseSpaceBasisNames` (`:25`) = `{Wigner, WignerMIC, Pauli, GellMann, GellMannMIC, Bloch, BlochSphere, GellMannBloch, GellMannBlochMIC, Wootters, Feynman, Tetrahedron, RandomMIC, RandomHaarMIC, RandomBlochMIC, QBismSIC, HesseSIC, HoggarSIC}`. When constructing `QuantumBasis["name"]` with one of these, `Picture -> "PhaseSpace"` is set automatically (`QuantumBasis.m:173`).

## Named circuit catalog (`$QuantumCircuitOperatorNames`)

`Kernel/QuantumCircuitOperator/NamedCircuits.m` — **44** named circuit constructions (verified `Length[$QuantumCircuitOperatorNames] == 44` at `ac32bda8`). Since the prior audit, new entries: `Trotterization[ops, order, steps, time]` (Suzuki-Trotter), `LeggettGarg`, `WignerCHSH`, `ControlledMultiplexer`, `StatePreparation`, `QuantumState`.

**State preparation** (`NamedCircuits.m`): `QuantumCircuitOperator[qs_QuantumState]` / `["QuantumState"[qs]]` / `["StatePreparation"[qs]]` now routes through `QuantumStatePreparation[qs, Method -> Automatic|"Multiplexer"|"BlockDiagonal"|"Classiq"]`. Automatic sends all-qubit `{2..}` states to the qubit `QuantumStateMultiplexer` (RZY amplitude-pair gates) and **any qudit state to the new `QuantumStateBlockDiagonal`** (block-diagonal uniformly-controlled gates per qudit, `QuditPrepGates`). `"Multiplexer"` is qubit-only (`::qubitonly`); non-pure input → `::nonpure`; unknown method → `::badmethod`.

**v2.0 list-form regressions fixed** (`797cf7b4` + others): `GrayOracle`, `Fourier`, `Multiplexer` dagger, `Grover`/`DeutschJozsa` dispatchers all migrated from `{Name, args}` to `Name[args]` call form.

### Algorithmic circuits
- `Bell` (`:459-460`) — `{H, CNOT}`. **Prep circuit, NOT the Bell state.**
- `Toffoli`, `Fredkin` (`:462-477`) — explicit 15+/18+ gate decompositions.
- `Magic` (`:646`) — `{S→1, S→2, H→2, CNOT→{2,1}}`.
- `GHZ[n]` (`:44-45`) — `{H, CNOT chain}`.
- `Graph[g, ...]` (`:34-41`) — graph-state preparation (`H` on each + entangling along edges).

### Grover family
- `GroverDiffusion[xs, gate]` / `GroverDiffusion0` / `GroverPhaseDiffusion` / `GroverPhaseDiffusion0` (`:48-118`) — diffusion operators (4 variants for ±0/+phase combinations).
- `Grover` / `Grover0` / `GroverPhase` / `GroverPhase0` (`:127-163`) — full Grover operator (oracle + diffusion). Accepts a Boolean formula or operator.

### Boolean / phase oracles
- `BooleanOracle[formula, vars, n, gate]` (`:203-241`) — ESOP-decomposed oracle. Picks smaller of positive vs negative form.
- `BooleanOracleR` (`:243-288`) — rotation-based oracle (RX/RY/RZ).
- `PhaseOracle` (`:334-363`) — phase-kickback variant.
- `GrayOracle` (`:290-315`) — Gray-code optimized angle decomposition. Helpers: `GrayMatrix`, `GrayOrders`, `BooleanGrayAngles` (`:317-332`).

### Quantum algorithms
- `Deutsch`, `DeutschPhase`, `DeutschJozsa`, `DeutschJozsaPhase`, plus oracle variants (`:366-417`).
- `SimonOracle`, `Simon` (`:422-454`) — random-secret default; accepts string `"01101"` or list `{0,1,1,0,1}`.
- `BernsteinVaziraniOracle`, `BernsteinVazirani` (`:484-509`) — secret string default `"101"`.
- `Fourier[n, m]` / `InverseFourier[n, m]` (`:533-547`) — QFT and inverse, `n` qubits, optional shift `m`.
- `PhaseEstimation[op, n, m]` (`:550-590`) — standard QPE template. **Qubit-only** (`MatchQ[op["OutputDimensions"], {2 ..}]`). `"PowerExpand" -> True` option to flatten the controlled-U^k.
- `Number[n, qubits]` / `PhaseNumber[n, qubits, h]` (`:512-531`) — encode integer `n` as a basis state.

### Composite / structural
- `Controlled[qc, control1, control0]` (`:592-596`) — wraps a circuit in controlled form by mapping each operator/sub-circuit to its controlled version.
- `Multiplexer[ops...]` / `Multiplexer[ops...] -> defaultK` (`:649-668`) — mux of `ops` controlled by binary index.
- `ControlledMultiplexer` (`:670-671`) — uses `QuantumLinearSolve`.
- `Switch[a, b]` (`:479-482`) — quantum-controlled switch via Cup/CSWAP/Cap.
- `Trotterization[ops, order, reps, const]` (`:631-644`) — order-1/2/even Suzuki recursion. **Order > 2 explodes** as `~5^(order/2 - 1)`.

### State preparation
- `QuantumState[qs]` / `StatePreparation[qs]` (`:739-744, :746`) — circuit that prepares `qs` (qubits only). Methods:
  - `Method -> Automatic` → `QuantumStateMultiplexer` (recursive RZY decomposition; `:719-737`).
  - `Method -> "Classiq"` → calls Classiq's `prepare_state` via Python (`Classiq.m`).

### Foundations / non-locality tests
- `CHSH[θ]` (`:749-761`) — Bell-test circuit with measurement angle.
- `LeggettGarg[θ]` (`:764-771`) — Leggett-Garg inequality circuit.
- `WignerCHSH[θ]` (`:773-791`) — CHSH in WignerMIC basis.

### Pauli-string shortcut
- `QuantumCircuitOperator["IXYHS"]` (`:599-601`) — chains of single-char ops in `{I, X, Y, Z, H, S, T, V, P}` decode as a tensor-product gate. Same alphabet as `QuantumOperator`.

## Qiskit integration (`QuantumCircuitOperator/Qiskit.m`)

### Forward direction (QF → Qiskit)
- `QuantumCircuitOperatorToQiskit[qco]` / `qco["Qiskit"]` / `qco["QiskitCircuit"]` (`:63-145`) — builds a Python Qiskit `QuantumCircuit`, returns `QiskitCircuit[bytes]` (pickled).
- Gate translation map `shortcutToGate` (`:11-61`): X|NOT → XGate, Y → YGate, Z|"1" → ZGate, H → HGate, S → SGate, T → TGate, V → SXGate, SWAP → SwapGate, RX/RY/RZ → respective Gates, U2 → U2Gate, U → U3Gate, Permutation → PermutationGate, GlobalPhase → GlobalPhaseGate, P → PhaseGate, PhaseShift[k] → PhaseGate[Sign[k] 2π/2^|k|], C/Controlled → Control, Dagger → Dagger, Barrier passthrough.
- **Unknown gates fallback**: `Unitary` with `NumericArray[Normal[N[matrix]], "ComplexReal32"]` — **single precision**.

### `QiskitCircuit[bytes]` properties
- `["Bytes"]` — raw pickled Python.
- `["Qubits"]`, `["Depth"]`, `["Ops"]`, `["Clbits"]` — basic info via `Eval`.
- `["Diagram", "Scale" -> n]` — matplotlib via `qiskitDiagram` (`:175-192`).
- `["Graphics"]` — LaTeX via MaTeX paclet (auto-installs).
- `["QuantumCircuit" | "QuantumCircuitOperator"]` (`:215-292`) — convert back to QF.
- `["Matrix"]` / `["QuantumOperator"]` (`:518-523`) — full unitary via `Operator(qc).data`. **Materializes 2^N × 2^N**.
- `["Decompose", n]` (`:516`) — `qc.decompose(reps=n)`.
- `["QASM"]` / `["QASM2"]` / `["QASM3"]` (`:552-554`) — QASM 2/3 via Qiskit's `qasm2.dumps` / `qasm3.dumps`. Uses `transpile` first.
- `["QPY"]` (`:557-570`) — QPY-format export, zlib-compressed bytes.
- `["Transpile", basisGates, opts]` (`:586-610`) — Qiskit transpile with optional basis gates and optimization level.
- `["Layout", opts]` (`:612-625`) — `plot_circuit_layout` via matplotlib.
- `["Validate"]` (`:628-639`) — FireOpal validation (qctrl).
- `qc[opts]` / `qc[qs, opts]` (`:525-527`) — apply circuit. Provider options: `"Provider" -> "IBM" | "AWS" | None`, `"Backend" -> name | Automatic`, `"Shots" -> 1024`, `"DensityMatrix" -> False`, `"FireOpal" -> False`. Returns: `QuantumState` (numeric statevector), `QuantumMeasurement` (counts dict ≤ 8 qubits), or raw result.

### QASM I/O
- `ImportQASMCircuit[file_or_string, backend, basisGates]` (`:642-660`) — parse QASM → `QiskitCircuit`. Uses Qiskit's `transpile` for backend conversion.
- `ImportQPY[bytes]` (`:572-582`) — load QPY (zlib-decompressed).

### Provider patterns (`:309-310`)
- `$IBMProvider = "IBM" | "IBMQ" | "IBMProvider" | "IBMRuntime"` — uses `qiskit_ibm_runtime.QiskitRuntimeService`. Auth via `SystemCredential["IBMQuantumPlatform_APIKEY"]` or explicit `"Token"`.
- `$AWSProvider = "AWS" | "AWSBraket" | "Braket"` — uses `qiskit_braket_provider.AWSBraketProvider`.
- Default — local `BasicProvider` or `AerSimulator` (with `method='statevector'` or `'density_matrix'`).

## Sub-package contexts (from `PacletInfo.wl`)

- `Wolfram`QuantumFramework`Init` ` — boot
- `Wolfram`QuantumFramework`Gates` ` — gate library (`Kernel/Gates.m`)
- `Wolfram`QuantumFramework`Experimental` ` — research-track (`Kernel/Experimental.m`)
- `Wolfram`QuantumFramework`ExampleRepository` ` — example library (`Kernel/ExampleRepository.m`)
- `Wolfram`QuantumFramework`DiagramPlot` ` — diagram plotting (`Kernel/DiagramPlot.m`)
- `Wolfram`QuantumFramework`SecondQuantization` ` — boson/fermion (`Kernel/SecondQuantization.m`)
- `Wolfram`QuantumFramework`QuantumOptimization` ` — VQE/QAOA/QNGD/etc. (`Kernel/QuantumOptimization.m`)

*(Public functions per sub-context to be enumerated in Phase 2 reads.)*

## PackageScope helpers (from `Kernel/Utilities.m`)

Not user-facing but valuable when building inside QF:
- Predicates: `SymbolicQ`, `nameQ`, `propQ`, `stateQ`, `orderQ`, `autoOrderQ`, `targetQ`, `measurementReprQ`, `emptyTensorQ`
- Tensor helpers: `tensorToVector`, `toTensor`, `tensorDimensions`, `tensorRank`, `identityMatrix`, `kroneckerProduct`, `projector`, `MatrixPartialTrace`, `blockDiagonalMatrix`
- Spectral: `eigensystem`, `eigenvalues`, `eigenvectors`
- Operator bases: `pauliMatrix`, `spinMatrix`, `fanoMatrix`, `GellMannMatrices`
- Geometry/frame: `RegularSimplex`, `GramMatrix`, `GramDual`
- Performance: `SparseArrayFlatten`, `MatrixInverse`, `matrixFunction`, `SetPrecisionNumeric`, `TranscendentalRecognize`
- Profiling / caching: `$QuantumFrameworkProfile`, `profile`, `Memoize` (rewritten `5b696969`: `Hash[Hold[args]]`-keyed `Association` + `"CacheQ" -> pred` option), **`cacheProperty`** (new `f9dc1cdc`, `HoldFirst`: the property-cache store. `cacheProperty[HeadProp[obj,prop,args], result]` commits the cache DownValue only when the held key is `FreeQ` of pattern constructs, replacing `Quiet[... = result, Rule::rhs]` in every head's dispatcher)

## Sub-package public API (from doc reads, 2026-04-29)

### `Wolfram`QuantumFramework`SecondQuantization`` (`Kernel/SecondQuantization.m`)

Doc: `Documentation/.../Tutorials/SecondQuantization.nb` (1245 LOC). Marked "New in: 1.4.4".

**State constructors:**
- `FockState[n[, size]]` — single-mode Fock `|n⟩`. Size defaults to `$FockSize=16`. Raises `FockState::clip` if `n ∉ [0, size-1]`; clips silently to nearest.
- `FockState[{n1, n2, ...}[, size]]` — multimode. Same size for every mode (truncation does NOT change per mode).
- `CoherentState[size][α]` — parametric coherent state. Use `CoherentState[][α]` with default size.
- `ThermalState[nbar[, size]]` — mixed thermal state with mean photon number.
- `CatState[size][α, φ]` — `𝒩(|α⟩ + e^{iφ}|-α⟩)`.

**Operators (all return `QuantumOperator`, with optional `order` arg = list of mode positions):**
- `AnnihilationOperator[[size,] [order]]` — bosonic `a`. Memoized (`SecondQuantization.m:122`). `a^†` via `Dagger` or `^†`.
- `DisplacementOperator[α[, size,] [order], "Ordering" -> "Normal"|"Weak"|"Antinormal"]` — `Exp[α a^† − α^* a]`. Default ordering is `"Normal"`. `"Weak"` gives the most accurate `UnitaryQ` for moderate truncation.
- `SqueezeOperator[ξ[, size,] [order], "Ordering" -> "Normal"|"Weak"|"Antinormal"]` — `Exp[(1/2)(ξ^* a^2 − ξ a^{†2})]`. Default `"Normal"` is most stable for large `|ξ|` at moderate truncation.
- `PhaseShiftOperator[θ[, size,] [order]]` — diagonal `e^{i θ a^† a}`. At `θ = 2π/d`, equals `QuantumOperator["Z"[d]]["Conjugate"]`.
- `BeamSplitterOperator[{θ, φ}[, size,] [order], Method -> "MatrixExp"|"Recurrence"]` — two-mode `Exp[θ(e^{iφ} a_1 a_2^† − e^{-iφ} a_1^† a_2)]`. `Method->"Recurrence"` is faster for **exact** angles (default `"MatrixExp"` is faster for numeric).
- `QuadratureOperators[[size,] [order]]` — returns `{X1, X2}` where `X1 = (a + a^†)/2`, `X2 = (a − a^†)/(2i)`. The commutator at index `size-1` deviates due to truncation.

**Diagnostics:**
- `OperatorVariance[ψ, op]` — `⟨op^2⟩ − ⟨op⟩^2`.
- `G2Coherence[ψ]` — `⟨a^† a^† a a⟩ / ⟨a^† a⟩^2`. Coherent → 1, Fock |1⟩ → 0, thermal → 2 (truncation bias toward < 2).
- `G1Correlation[ψ, {{r1,t1}, {r2,t2}}]` — first-order field correlation.

**Phase-space representations (numeric grid, return `InterpolatingFunction`):**
- `WignerRepresentation[ψ, {xmin,xmax}, {pmin,pmax}, opts]`
- `HusimiQRepresentation[ψ, {xmin,xmax}, {pmin,pmax}, opts]`
- Options: `"GaussianScaling" -> Sqrt[2]`, `"GridSize" -> 100`. Increase `"GridSize"` for higher Fock-number states; default 100 is too coarse for `|n⟩` with `n ≥ 2` (off by ~0.07 at origin for Fock |2|).

**Globals:**
- `$FockSize` — current default truncation. Initialized to 16. Memoizes implicitly (changing it after building operators leaves stale `AnnihilationOperator` definitions). Reset with `SetFockSpaceSize[]`.
- `SetFockSpaceSize[n]` / `SetFockSpaceSize[]` — mutate / reset.

### `Wolfram`QuantumFramework`QuantumOptimization`` (`Kernel/QuantumOptimization.m`)

Doc: `Documentation/.../Tutorials/QuantumOptimization.nb` (4810 LOC).

**Variational helpers:**
- `GenerateParameters[n, m, "Symbol" -> "θ"]` — list of `n*m` symbols `θ1..θ(n*m)`. Use for an `n`-qubit, `m`-layer ansatz.
- `ParametrizedLayer[gate, [qubits,] index, "Symbol" -> "θ"]` — `Sequence[gate[θi]->{q}, ...]`. With only `index`, qubits default to `Range[Length[index]]`.
- `EntanglementLayer[cgate, qubits, "Entanglement" -> Automatic|"Linear"|"Circular"|"ReverseLinear"|"Pairwise"]` — emits `Sequence[cgate->{i,j}, ...]`. Linear chain, all-pairs (Automatic), ring (Circular), pairwise even/odd layers.

**Optimizers:**
- `GradientDescent[f, init, "LearningRate" -> 0.8, "MaxIterations" -> 50, "Gradient" -> None]` — vanilla GD. Returns full trajectory (list of param vectors, one per step).
- `QuantumNaturalGradientDescent[f, metric, "InitialPoint" -> Automatic, "LearningRate" -> 0.8, "MaxIterations" -> 50, "Gradient" -> None]` — same trajectory shape, but step `θ_{t+1} = θ_t − η · g^+(θ) · ∇f(θ)` using pseudo-inverse of `metric`.
- `FubiniStudyMetricTensor[QuantumState[..., "Parameters" -> {...}], property]` — Fubini-Study metric on parameterized state. Properties: `"Matrix"`, `"MatrixForm"`, `"Parameters"`, `"SparseArray"`, `All`. Default returns a `QuantumOperator`.
- `SPSRGradientValues[generator[θ], "PauliX"|"PauliY"|"PauliZ"]` — exact stochastic-parameter-shift gradient values.
- `ASPSRGradientValues[generator[θ], pauli, H]` — approximate SPSR (faster, with H = the non-θ part).

**Linear solver (top-level — exposed in `PacletInfo` `Symbols`):**
- `QuantumLinearSolve[m, b[, prop|All, opts]]` — multiplexer-based VQLS. Works on numeric matrices, requires `Length[b]` a power of 2.
  - **Properties (`prop`):** `"Result"` (default), `"Ansatz"` (return parametric `QuantumState`), `"CircuitOperator"` (return parametric `QuantumCircuitOperator`), `"GlobalPhase"` (the recovered `e^{iφ}`), `"OptimizedParameters"`, `"Parameters"`. `All` returns an Association of all.
  - **Options:** `"Ansatz" -> Automatic | QuantumState | QuantumCircuitOperator`, `"GlobalPhaseAccuracy" -> 10^-5`, `AccuracyGoal -> Automatic`, `MaxIterations -> Automatic`, `Method -> Automatic | "NelderMead" | "DifferentialEvolution" | "SimulatedAnnealing" | "RandomSearch" | "Couenne"`, `PrecisionGoal -> Automatic`, `WorkingPrecision -> MachinePrecision`.
  - Internally uses `QuantumOperator[m]["PauliDecompose"]` then `QuantumCircuitOperator["Multiplexer" @@ Keys[σ]]`. Cost is `1 - QuantumDistance[..., "Fidelity"]`.

**Python integration:**
- `ClassiqSetup[prop, "CheckDependencies" -> {...}, "InstallPackages" -> {...}]` — boots a Python session and installs Classiq. Properties: `"Evaluators"`, `"Session"` (returns `ExternalSessionObject`), `"ClassiqVersion"`, `All`.

### `Wolfram`QuantumFramework`Experimental`` — `QuantumWignerTransform`

Doc explicitly tags `QuantumWignerTransform` as **[EXPERIMENTAL]** (`Documentation/.../Symbols/QuantumWignerTransform.nb:19`). API may change. Equivalent to basis-changing `QuantumState[obj["Double"], QuantumBasis["Wigner"[d]]]`.

Related state properties (also experimental):
- `qs["PhaseSpace"]` — tensor `𝒲_{p1,q1,p2,q2,...}`. Dimensions: `d_i × d_i` per qudit if `d_i` odd, `2 d_i × 2 d_i` if even.
- `qs["QuasiProbability"]` — flattened matrix `𝒲_{p_t, q_t}` over total dimension. **Same as PhaseSpace only for single qudit.**
- `qs["TransitionPhaseSpace"]` / `qs["TransitionQuasiProbability"]` — 4-quadrant `2 d × 2 d` tensor/matrix. Quadrants encode position↔position, position↔momentum, etc.
- `qs["TransitionGraph"]` — `Graph` of phase-space transitions weighted by transition quasi-probability.
- Norm: pure states on an n-sphere of radius `1/Sqrt[d]` (odd `d`), `1/(2 Sqrt[d])` (even `d`). Mixed states inside.
- Overlap: `Tr[ρ1 . ρ2] = d * Total[W1 * W2, 2]` (odd `d`). Useful for fast inner products in phase space.

### Top-level confirmed API additions (from docs, beyond `Usage.m`)

- `QuantumOperator["PauliDecompose"]` (tutorial-only) — Pauli decomposition Association, used by `QuantumLinearSolve` for the multiplexer encoding.
- `QuantumCircuitOperator["QASM"]` / `["TensorNetwork"]` / `["QiskitCircuit"]` / `["Topology"]` / `["CircuitOperator"]` / `["Width"]` / `["Shift", n]` — properties documented in QuantumCircuitOperator ref page that aren't called out in `Usage.m`.
- `QuantumState["Operator"]` returns `QuantumOperator` (projector as operator). `QuantumState["Projector"]` is the higher-dim "Pure map" form.
- `QuantumState["SchmidtBasis"]["Disentangle"]` returns the list of `c_i |u_i v_i⟩` terms; `["Decompose"]` recurses Schmidt until indecomposable.
- `QuantumState["Bipartition", spec]` accepts: a pair `{{i1,...}, {j1,...}}`, a single side `{i, j, ...}` (complement is the other side), or a target dimension.
- `QuantumState["Bend"]` / `["Unbend"]` — vectorize density matrix into a higher-dim pure state and back. `Length[ψ_bent] = d^2`.
- `QuantumState["Permute", Cycles[...]]` — qudit permutation.

### Doc-confirmed named states (extension of `QuantumState["name"]`)

`Documentation/.../Symbols/QuantumState.nb:75-103`:
`"0"`, `"1"`, `"00"`...`"111"` (any binary string), `"+"`, `"-"`, `"L"`, `"R"`, `"Plus"`, `"Minus"`, `"Left"`, `"Right"`, `"Bell"`, `"PhiPlus"`, `"PhiMinus"`, `"PsiPlus"`, `"PsiMinus"`, `"BlochVector"[{rx,ry,rz}]`, `"Register"[s]` / `"Register"[s,i]`, `"UniformSuperposition"[s]`, `"UniformMixture"[s]`, `"RandomPure"` / `"RandomPure"[s]` / `"RandomPure",{d1,d2,...}`, `"RandomMixed"`, `"GHZ"` / `"GHZ"[s]`, `"W"` / `"W"[s]`, `"Werner"[p, d]`, `"Dicke"[n,k]` / `"Dicke"[{k1,k2,...}]`, `"Graph"[g]`.

Mixed-format strings combine: `QuantumState["+0L"]` = `+ ⊗ 0 ⊗ L`.

### Doc-confirmed named bases

`Documentation/.../Symbols/QuantumBasis.nb:31-72`:
- Hilbert: `"I"`/`"Computational"[d]`, `"Bell"`, `"Pauli"`, `"X"[d]` / `"Y"[d]` / `"Z"[d]` (generalized Pauli), `"PauliX"`/`"PauliY"`/`"PauliZ"`, `"JX"[j]` / `"JY"[j]` / `"JZ"[j]` (spin), `"Fourier"[d]`, `"Schwinger"[d]`, `"Dirac"`.
- IC phase space: `"WignerMIC"[d]`, `"GellMannMIC"[d]`, `"Tetrahedron"` / `"Tetrahedron"[a,b,c]`, `"HasseSIC"`, `"HoggarSIC"` (8-dim), `"QBismSIC"[d]` (numeric SIC up to `d=151`, **downloads from GitHub**), `"RandomHaarMIC"`, `"RandomBlochMIC"`.
- Non-IC phase space: `"GellMann"[d]`, `"Wigner"[d]`, `"Wootters"[p]` (`p` prime), `"Bloch"`.

### Doc-confirmed named circuits (`QuantumCircuitOperator["name"]`)

`Documentation/.../Symbols/QuantumCircuitOperator.nb:73-100`:
`"Bell"`, `"Cup"` / `"Cap"`, `"Toffoli"`, `"CHSH"`, `"Magic"` (computational → Bell basis transformer), `"Graph"[g]`, `"Grover"[f]` / `"Grover0"[f]` / `"GroverPhase"[f]` / `"GroverPhase0"[f]` / `"GroverDiffusion"[n]` / `"GroverDiffusion0"[n]` / `"GroverPhaseDiffusion"[n]` / `"GroverPhaseDiffusion0"[n]`, `"BooleanOracle"[f]` / `"BooleanOracleR"[f]` / `"PhaseOracle"[f]`, `"Fourier"[n]` / `"InverseFourier"[n]`, `"BernsteinVazirani"[s]` / `"BernsteinVaziraniOracle"[s]`, `"Simon"[s]` / `"SimonOracle"[s]`, `"PhaseEstimation"[u, n]`, `"DeutschJozsa"[k, n]` / `"DeutschJozsaPhase"[k, n]` / `"DeutschJozsaBooleanOracle"[k, n]` / `"DeutschJozsaPhaseOracle"[k, n]`, `"Trotterization"[ops, order, steps, time]`, `"PhaseNumber"[n, m, boo]`, `"LeggettGarg"`, `"GrayOracle"[angles, prec, target]`, `"Multiplexer"[op1, op2, ...]`, `"QuantumState"[state]` / `"StatePreparation"[state]`.

**Option notes** (`8da732a2`, post-v2.0):
- `"BooleanOracleR"[f]` accepts a `rotationGate` option of form `("RX"|"RY"|"RZ")[_ ? NumericQ]`; default is `"RZ"[Pi]`. The v2.0 named-call form replaces the old `{"RZ", Pi}` default.
- `"Graph"[g]` accepts an optional third arg `gate` defaulting to `"C"["1"]` (was `{"C", "1"}` pre-`a5958975`) — a controlled-by-`"1"` gate spec.

### Doc-confirmed named operators (extension of `QuantumOperator["name"]`)

`Documentation/.../Symbols/QuantumOperator.nb:47-122` — most relevant additions to the catalog already in the function-index:
- `"Switch"[op1, op2]` — quantum switch (commute/anti-commute probe).
- `"Curry"[basis|bases]` / `"Uncurry"[basis|bases]` — pack/unpack qudits without changing total dimension.
- `"Cup"[d]` / `"Cap"[d]` — identity-with-no-input / identity-with-no-output (categorical caps).
- `"Spider"[basis, angles] -> inputs -> outputs` — generic spider; `"ZSpider"`, `"XSpider"`, `"Copy"`, `"Decohere"`, `"Encode"`, `"Measure"`, `"Discard"`, `"Marginal"` — specializations.
- `"Trace"[d]` — partial-trace operator. `"Liouvillian"[H, {L1,...}, {γ1,...}]` and `"Hamiltonian"[H, {Ls}, {γs}]` (already in index).
- `"Right"[op]` / `"Left"[op]` — superoperators acting on the right/left of a doubled state.
