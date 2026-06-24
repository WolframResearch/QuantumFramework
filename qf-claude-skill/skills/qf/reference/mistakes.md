# Recurring Mistakes


## HEAD resync gotchas (2026-06-19, `b4fcdb31`)

### `PauliStabilizer` now takes whole circuits as `ps[{spec, ...}]` / `ps[spec, ...]` (multi-gate bracket forms)
- **What changed** (`c3f21a43`; `PauliStabilizer.m:141-168`, `StabilizerFrame.m:263-266`, catalogs in `Compiled.m:57-58`): `ps[{"H", "CNOT"}]` (list) and `ps["H", "CNOT"]` (braceless variadic) apply a whole circuit, routing to the compiled `ps["ApplyCircuit", ...]` engine. Specs are arrows (`"H" -> 1`), bare names with a default target (one-qubit → qubit 1, two-qubit → `{1, 2}`), or call-forms (`"P"[θ]`, `SuperDagger["S"|"V"|"T"]`); they can be mixed in one list.
- **The `ps[{"X"}]` vs `ps["X"]` distinction (the trap that motivated the design):** a **single-gate list** `ps[{"X"}]` *applies* the X gate; the **single string** `ps["X"]` is the property accessor and returns the X-block of the tableau (a `{_List, _List}`, **not** a new stabilizer). They never collide because one argument is a `_List` and the other a `_String`. Likewise `ps[{1, 2}]` stays a **measurement** (the more specific `{___Integer}` rule), not a two-gate circuit.
- **A single bare gate stays inert:** the variadic form requires `Length ≥ 2` (`stabilizerGateSpecSeqQ`), so `ps["H"]` does **not** apply H. it stays unevaluated as `PauliStabilizer[...]["H"]`. To apply one gate use the explicit-target form `ps["H", 1]` or the single-gate list `ps[{"H"}]`. An unrecognized token (`ps["Wobble", "Splork"]`) is rejected by the `stabilizerGateSpecQ` membership guard (the dispatch rule does not fire; the call stays unevaluated), **not** via a message. `PauliStabilizer::badgate` exists but is unreachable from the public bracket path (it fires only if `stabilizerDefaultTarget` is called directly on a non-string head).
- **Catalog hygiene:** the dispatch recognizes `$stabilizerCliffordNames` = `{H,S,X,Y,Z,V,T,CNOT,CX,CZ,SWAP}`; this is a **third** gate-name list, deliberately distinct from `$stabilizerGateCodes` (the C-kernel encoding, which omits the non-Clifford `V`/`T` that decompose / route to a frame) and from `sfGate1`/`sfGate2` (the frame Pauli-conjugation tables keyed by normalized names like `"Sdg"`). A new gate must be added to whichever of the three its dispatch / engine / frame path actually supports; they are not interchangeable.

### `StabilizerFrame` single-gate rule was silently corrupting `f["H", "S"]` (FIXED)
- **What changed** (`c3f21a43`; `StabilizerFrame.m:201, 212`): the Clifford-distribute and `SuperDagger` rules were `f_StabilizerFrame[gate_String, args___]` / `[SuperDagger[gate_String], args___]`. With the unrestricted `args___`, a call like `f["H", "S"]` captured the **second gate name** `"S"` as `args` and handed it to `sfConjugatePauli` as a **qubit index**, silently corrupting the relating Paulis (no message, wrong amplitudes). The rules are now `args : ___Integer`, so a trailing gate name no longer matches and instead falls through to the new multi-gate variadic rule. Recall that a frame is what you get from a non-Clifford gate (`PauliStabilizer[2]["T", 1]` is a `StabilizerFrame`), so this bug was reachable from ordinary `T`-circuit code, not just hand-built frames.
- **The mistake** (now fixed): assuming `sf["H", "S"]` chained two gates pre-`c3f21a43`. it did not. On the working tree it now equals `sf["H", 1]["S", 1]` via the variadic→list→fold path. Frame multi-gate forms `sf[{spec, ...}]` / `sf[spec, ...]` fold per gate (no compiled engine for a frame) and track the relating Paulis coherently; verified to amplitude 0 against the explicit-target chain (`AuditMatrix.wlt` `Audit-Multi-Frame-{List,Variadic}-Coherent`).

## HEAD resync gotchas (2026-06-17, `a601a187`)

### A gate-built `StabilizerFrame`'s dense `["StateVector"]` is now amplitude-exact (narrowed caveat)
- **What changed** (`61bc0e39`): the old dense readout summed each component's `+1` eigenvector materialized **independently**, each with an arbitrary global phase from its `NullSpace` representative, so the coherent sum mixed incoherent phases: `T|+⟩` came out with its two amplitudes **swapped** (fidelity 0.5) and observables read from the frame were sign-flipped. A **gate-built** frame now carries a relating Pauli per component (`"Paulis"` key), so `f["StateVector"]` materializes one reference component and relates the rest by the actual Pauli operators, recovering all relative phases.
- **The (narrowed) caveat**: as of HEAD, a single gate-built frame's `["StateVector"]` and any single-state `⟨ψ|O|ψ⟩` are **amplitude-exact up to one overall global phase**, the same contract as a bare `PauliStabilizer`, not a weaker one. Do **not** still warn that the frame is "only a combinatorial tool, correct up to a per-component phase, not for amplitude readout"; that was the pre-`61bc0e39` reality. Verified: `T|+⟩ = {1, e^{iπ/4}}/√2` exactly, `⟨Y⟩` on `T|+⟩` `= +0.7071` (was `-0.7071`), amplitude-exact vs dense across 216 random Clifford+T circuits.
- **The one residual numeric caveat**: the **phase of an inner product between two distinct, separately-gauged stabilizer states**. `Stabilizer/InnerProduct.m`'s frame inner product still routes through `["StateVector"]`, so its **magnitude** is correct but the **relative phase** is gauge-limited (each state fixes its own global phase independently). Magnitudes, single-state expectations, and `["StateVector"]` are fine.
- **Frames without relating Paulis are unchanged**: a hand-built `StabilizerFrame[{{c, ps}, ...}]` or any `Plus`-combined frame (components need not share a tableau) keeps the old independent per-component readout. The `2^k` component-doubling-per-non-Clifford-gate law and the lack of rank compression are unchanged; cost of the fix is ≈ +15% `ByteCount`, +8-12% build time per frame. Full account: `OngoingProjects/Stabilizer/magic-and-pathsums/StabilizerFrame-NonClifford-Blowup.md` §4.

### Symbolic / parametric `QuantumMeasurement` no longer spams messages on access or display (FIXED)
- **What changed** (`da372b3f` + `fc430744`): `QuantumMeasurement["Entropy"]` used to route through `Information[CategoricalDistribution[outcomes, probs]]`, which falls back to Monte-Carlo sampling and emitted `RandomVariate::unsdst` + `Extract::psl1` whenever the outcome probabilities were symbolic (any parametric circuit), including on **plain display**, because the summary box called `N @ qm["Entropy"]`. `"Entropy"` is now a closed-form Shannon entropy (`-Σ p Log[b, p]`, `logBase` arg + list-form), and a `HoldRest` `whenNumericProbabilities` guard (no `Quiet`) returns `Indeterminate` for the sibling Monte-Carlo properties (`"DistributionInformation"`, `"TopProbabilities"`, `"SimulatedMeasurement"`, `"SimulatedCounts"`, `"SimulatedStateMeasurement"`, `"TopStateProbabilities"`) instead of leaking `MultinomialDistribution::vprobprm`/`CategoricalDistribution::elmntavsl`/`KeyMap::invak`.
- **The mistake** (now obsolete): expecting message spam from a parametric measurement and wrapping it in `Quiet`. Unnecessary now; verified clean (no messages on `qm["Entropy"]`, `qm["Probabilities"]`, `qm["SimulatedCounts", n]`, or display for a `QuantumCircuitOperator[{"RY"[θ]->1, {1}}, "Parameters"->{θ}][]` measurement). `"Categories"`/`"Probabilities"`/`"ProbabilityTable"`/`"ProbabilityArray"` keep their symbolic values; `"Distribution"`/`"NDistribution"` still return a symbolic `CategoricalDistribution`. Shipped on `main`.

### Error-aware Qiskit transpilation is now **shipped**, not working-tree
- **What changed** (`f370f568`): the `"OptimizationLevel"` option (`IBMJobSubmit`/`qiskitQASM`/`qiskitPrimitiveSubmit`, Automatic → 1), `QiskitTarget["FromBackend" -> name, "Provider" -> p]`, the offline `"GenericV2"` fake-device provider, `QiskitTarget["InstructionProperties"]`, and the `instruction_properties` carrier (applied via `update_instruction_properties`) are all committed. Any prior note saying these are "uncommitted in the working tree at `f9dc1cdc`" is stale.
- **Subtleties to keep in mind** (carried over from the traps investigation): `QiskitTarget["FromBackend" -> "GenericV2"]` with **no** `"Provider"` silently yields an idealized Aer-style target rather than erroring, and `build_target`'s `qualink_skipped_props` is a write-only metadata key (recorded, never read back by QF), so a Target row dropped for an unknown gate/locus is logged into circuit metadata but not surfaced to the user. The `"Target"` option is now `Switch`-validated: a failed `QiskitTarget[spec]` / `$Failed` / stray value is rejected with a `Failure` naming the head instead of an opaque Python error.

## HEAD resync gotchas (2026-06-14, `f9dc1cdc`)

### `$QuantumFrameworkPropCache = False` is no longer a meaningful speed lever
- **What changed** (`dfc741dc`/`f9dc1cdc`): the property-cache machinery repair removed the per-access "bookkeeping tax" (direct `"Basis"` getters + widened cache exclusions so $O(1)$ structural getters skip the eligibility guard; `cacheProperty` store). The default-on warm apply now costs about what the cache-off path costs.
- **The mistake**: carrying forward the old advice "warm up, then set `$QuantumFrameworkPropCache = False` for repeated numeric circuit evaluation" expecting a ~3.7× win. Measured now: warm 105-gate 12-qubit apply 273 ms default-on vs 222 ms cache-off (~1.2×). The flag still exists as an escape hatch, but it is no longer the lever it was. Output was, and remains, bit-identical with the flag on or off.
- **Right way**: just apply; the default is fast. Reach for `False` only if you suspect a stale cached value across a substitution (the no-auto-invalidation entry still stands).

### `Rule::rhs` no longer fires from heavy property access (don't `Quiet` it yourself)
- **What changed** (`f9dc1cdc`): the cache store is now `cacheProperty`, which commits a cache DownValue only when the held key is free of pattern constructs, so the latent over-matching DownValue (which the old `Quiet[... = result, Rule::rhs]` masked) cannot form. `Memoize` was likewise rewritten to a `Hash`-keyed `Association` (`5b696969`).
- **The mistake**: wrapping QF property-heavy code in `Quiet[..., Rule::rhs]` "to be safe". Unnecessary now, and it would mask a *real* `Rule::rhs` from your own code. If you see `Rule::rhs` from current QF, it is a genuine bug worth reporting.

### The misshapen-state OOM and the `Table[0,{n}]` register misparse are now guarded (but still wrong inputs)
- **What changed** (`7e12e8f0`): `QuantumOperator::broadcast` + `$Failed` when the order-driven multiplicity broadcast would exceed `$QuantumOperatorBroadcastLimit` ($2^{24}$): the apply-path kernel-killer (`{psi -> Range[12]}` with a 4-qubit `psi` building a `16^12`-amplitude object) is now a fast `Failure`. And `QuantumState[vector, dim]` emits `QuantumState::padded` when zero-padding rounds the qudit count up.
- **The mistake**: still passing `QuantumState[Table[0, {n}], 2]` as an n-qubit register. It is a 4-qubit **norm-0** state; you now get a `QuantumState::padded` warning, but the object is still wrong. Use `QuantumState["Register"[n]]`. Likewise a state→order size mismatch now errors instead of OOMing, but it is still a mismatch to fix. (Note: named `QuantumOperator["H", hugeOrder]` takes the Fourier-kronecker path, not the broadcast rule, so a huge **named**-gate order still builds eagerly.)

### Stabilizer: `SameQ` across packed and canonical storage fails; use `ApplyCircuit` for big circuits
- **What changed** (June 2026 perf rounds): single Clifford gates return a **packed** `PauliStabilizer` (`Packed.m`), while `ps["ApplyCircuit", specs]` and constructors return the **canonical** rank-3 object. Both are valid everywhere (`PauliStabilizerQ` handles both).
- **The mistake**: comparing two stabilizer objects with `===`/`SameQ` across storage forms (returns `False` even for the same state), or folding a long gate stream one object-gate at a time (pays object wrap/unwrap; ~25× the compiled route) or building a `QuantumCircuitOperator` from 10⁴ gate objects (host-framework construction dominates).
- **Right way**: compare `ps1["Tableau"] === ps2["Tableau"]` (or `["Stabilizers"]`). For long Clifford streams call `ps["ApplyCircuit", specs]` (or `qc[Method -> "Stabilizer"]`, which routes there). `PauliStabilizer["Random"[n]]` is now viable to n≈538 (the `\mathbb{F}_2`-inverse and Mallows-underflow fixes); past that it still emits `General::munfl`.

## HEAD resync gotchas (2026-06-11, `77b9ccba`)

### `obj["Properties"]` no longer returns the full property set: use `"AllProperties"`
- **What changed** (post-`ac32bda8`): on every wrapped head (`QuantumState`, `QuantumOperator`, `QuantumChannel`, `QuantumHamiltonianOperator`, `QuantumMeasurementOperator`, `QuantumMeasurement`, `QuantumCircuitOperator`), `obj["Properties"]` now returns only the head's **own** static list (`$Quantum…Properties`); it no longer unions in the wrapped object's properties. The full set moved to the **new** `obj["AllProperties"]`.
- **The mistake**: `MemberQ[qo["Properties"], "Eigenvalues"]` (or any property inherited from the wrapped object) now returns `False`, even though `qo["Eigenvalues"]` still works. Live: a `QuantumOperator` has 67 `"Properties"` but 209 `"AllProperties"`; the operator answers to all 209.
- **The right way**: to test "does this object respond to property X?", query `obj["AllProperties"]`, not `obj["Properties"]`. The split is deliberate: `"Properties"` is now an O(1) stored list (cheap to display), `"AllProperties"` recomputes the union (and is on every head's cache-prevention list). Accessor delegation is unaffected; only code that *introspected* the old union breaks.

### `IBMJob["ExecutedCircuit"]` is now an **Action**, not a **Property**
- **What changed** (`a44c1eea`): `IBMJob["Properties"]` lists only local, network-free accessors. `"Refresh"`, `"Cancel"`, and `"ExecutedCircuit"` moved to the new `IBMJob["Actions"]` set because each needs a live connection.
- **The mistake**: `AssociationMap[job, job["Properties"]]` no longer triggers `"ExecutedCircuit"` (so it won't 404 on the presigned VPC-private COS URL, won't re-query, won't cancel). That's the intended fix, but code that walked `"Properties"` expecting `"ExecutedCircuit"` in the dump won't find it. Call the three Actions explicitly by name; they're still valid keys.

### TN contraction tuning goes through `Method -> {"TensorNetwork", subOpts}`, not top-level options
- **What changed** (`77b9ccba`): `Options[quantumCircuitApply]`/`Options[quantumCircuitCompile]` no longer Join in `Options[TensorNetworkApply]`/`Options[TensorNetworkCompile]`. A TN option (e.g. the new `"Path"`) passed as a bare top-level option to `qco[state, …]` / `qco["QuantumOperator", …]` is now silently dropped (or errors as an unknown option).
- **The right way**: nest it under the `Method` sub-option list: `qco[state, Method -> {"TensorNetwork", "Path" -> path}]`. See `qf-tn-integration.md` §4.

### Named-basis jump operators ("J-"/"J+") and the Liouvillian basis convention
- **The trap**: `QuantumOperator["J-"]` lives in the spin-J basis, whose element ordering is *reversed* relative to the computational basis: `["Matrix"]` (raw) is `{{0,1},{0,0}}` but the faithful computational representation is the transpose `{{0,0},{1,0}}`. So `"J-"` lowers the spin projection $m$, mapping $|0\rangle = m{=}{+}1/2$ to $|1\rangle = m{=}{-}1/2$, the opposite orientation of the qubit-community $\sigma_- = |0\rangle\langle 1|$ read off the raw matrix. With unequal rates on a `"J-"`/`"J+"` pair, an orientation flip swaps the steady-state populations.
- **What was wrong (fixed 2026-06-11)**: `QuantumEvolve`'s direct path rebases jumps into the Hamiltonian's basis (`QuantumEvolution.m:80`), but `QuantumOperator["Liouvillian"[...]]` built the superoperator from raw matrices and attached the jump's own basis. Both were *faithful* (same physical dynamics, verified to 1e-9 in trace distance), but the superoperator route returned states in the hidden J basis, so raw `DensityMatrix` / `ProbabilitiesList` / `BlochVector` reads came out swapped vs the direct route; the QuantumEvolve doc page's four "equivalent" Lindblad routes visibly disagreed.
- **The fix**: the `"Liouvillian"`/`"Hamiltonian"` constructors now rebase `H` and every jump into the Hamiltonian's output basis (computational when `H === None`) before building (`NamedOperators.m:900`). All routes return states in the same basis; raw reads agree. Regression: `Tests/Liouvillian.wlt` "Named-J-*" tests.
- **The convention to remember**: `"J-"` at rate $\gamma$ relaxes toward $|1\rangle$ (the ground state of $H \propto +Z$, since `"JZ"` is $+\sigma_z/2$ computationally). For the textbook $\sigma_-$ decaying toward $|0\rangle$, use the explicit matrix `QuantumOperator[{{0,1},{0,0}}]`.
- **Date**: 2026-06-11.

## New-surface gotchas (2026-06-09 re-audit, `ac32bda8`)

### KAK / `["Decompose"]` reconstructs only up to a global phase
- **The mistake**: assuming `QuantumOperator[op["KAK"]]` equals `op` exactly. **Verified**: `QuantumOperator["CNOT"]["KAK"]` matches CNOT up-to-phase (`Abs[Tr[U^\[Dagger] V]] == 4`) but `Norm[U - V] =!= 0`. The numeric path computes `residualPhase = Arg[Tr[ConjugateTranspose[circMat] . mat]/4]` and prepends a `GlobalPhase` to fix this, but in the CNOT case the residual is not nulled.
- **The right way**: compare KAK / `["Decompose"]` output up to a global phase (e.g. `Abs[Tr[U^\[Dagger] V]] == Length[U]`, or `QuantumDistance`/fidelity), not by exact matrix equality. The decomposition's *canonical coordinates* and *local factors* are correct; only the overall phase may be off. (`QuantumOperator/Properties.m`, numeric `TwoQubitKAK` path.)

### A registered name with a bad call shape now returns `Failure["InvalidArguments"]`, not an unevaluated expression
- **What changed** (`ac32bda8`): every Named* head (`QuantumState`, `QuantumOperator`, `QuantumMeasurementOperator`, `QuantumChannel`, `QuditBasis`, `QuantumCircuitOperator`) gained a last-resort `name_String[args___] /; MemberQ[$...Names, name]` rule emitting `::invalidArgs`. Previously such a call could fall through unevaluated, recurse, or hit `$RecursionLimit`.
- **The right way**: wrap constructor calls whose arguments come from data in `Enclose`/`Confirm` and handle `Failure["InvalidArguments", ...]` (distinct from `Failure["InvalidName", ...]`, which fires for an *unregistered* name).

### `QuantumQASM` / `["QASM"]` is qubit-only
- **The mistake**: exporting a circuit/operator with any non-2-dim wire to OpenQASM. OpenQASM models qubit registers only, so QF returns `Failure["QuantumQASM", <|..., "NonQubitDimensions" -> {...}|>]` rather than emitting malformed QASM. The dimension guard runs per operator over `qco["Flatten"]["Operators"]` (measurement eigen-registers / channel trace-registers at non-positive orders are excluded).
- **Note**: import (`QuantumQASM[src]`) and the default export are **native WL (no Python)**; only an extra-argument export (gate set / Target / option / provider) needs the Python/qiskit session. `ImportQASMCircuit` is deprecated — use `QuantumQASM`.

### Stabilizer gate updates drop the global phase by design
- **The mistake**: expecting `PauliStabilizer[qs][gate]["State"] === qs[gate]` exactly. The Clifford gate rules in `GateUpdates.m` deliberately do **not** forward `"GlobalPhase"` (recovering it would need O(2ⁿ) state materialization per gate, defeating the AG advantage). The contract is **correct up to a global phase**.
- **The right way**: for exact equality, re-run the constructor on the post-gate *state*: `PauliStabilizer[qs[gate]]["State"] === qs[gate]` (the constructor path recovers and stores the phase). The `PauliStabilizer[qs]`/`[qo]` constructors *do* store a `"GlobalPhase"` so a freshly-built object round-trips exactly.

### `ps["InnerProduct", other, Method -> "ClosedForm"]` is magnitude + sign only, and needs concrete signs
- The O(n³) closed form (GarMarCro12) returns the **magnitude** `2^(-s/2)` and a ±1 sign for perfectly-aligned states; the full **complex phase is TODO**. It also requires `ConcretePauliStabilizerQ` inputs (numeric signs) — symbolic signs emit `::expectationdim` and `$Failed`. The default `Method -> "Direct"` materializes vectors (O(2ⁿ), full phase) — use it for n ≤ ~8, ClosedForm for larger n when only the magnitude is needed.

### `SymbolicMeasure`: a deterministic measurement after a symbolic one is not stamped into the signs
- A deterministic outcome that follows a prior symbolic measurement (e.g. Bell `ZZ` correlation `m₂ = m₁`) is computed correctly by AG (the post-state is right) but the outcome polynomial is **not** written into the stabilizer signs, so `SampleOutcomes` cannot recover the correlation from signs alone. The kernel tracks these via an `"Outcomes"` map in the association, applied iteratively by `SubstituteOutcomes`/`SampleOutcomes`. For workflows needing explicit outcome correlation, prefer `ps["M", q]` (Association branching) or a `StabilizerFrame`. (`Stabilizer/SymbolicMeasure.m`.)

### `CliffordChannel[qc_QuantumChannel]` handles only deterministic single-Pauli channels
- A stochastic named channel (BitFlip/PhaseFlip/Depolarizing…) → `::stochastic` + `$Failed`. For tableau-level stochastic application use `qc[ps]` (returns a `{probability, ps_after_pauli}` mixture list, `HybridInterop.m`). `cc[ps]` itself only implements identity / state-prep / matching-input cases (else `::stateevol`).

### `GraphState[ps]` only accepts graph-form stabilizers; `LocalComplement` VOPs are stubbed
- `GraphState[ps_PauliStabilizer]` requires each stabilizer row to be Xᵢ on the diagonal and Z/I elsewhere (no Y, no off-diagonal X) → else `::nongraph` + `$Failed`. A general stabilizer needs a local-Clifford conversion (AndBri05 Lemma 1) first, which is **not yet implemented**. `LocalComplement` toggles edges but does **not** update VOPs (only identity VOP supported) — the returned graph state is correct as a graph but the local-Clifford correction is dropped.

### `PauliStabilizer[qs]` asserts `qs` is an actual stabilizer state
- The 4ⁿ-tomography constructor `ConfirmAssert`s that exactly 2ⁿ Paulis have `|<P>| ≈ 1`; a non-stabilizer state fails the Enclose with a clear message. It is O(4ⁿ) — practical only to n ≈ 8.

### `QuantumAdiabaticEvolve` resolves to `` Global` `` on a bare load
- Like the QHO head, `QuantumAdiabaticEvolve`/`QuantumAdiabaticEvolution` are `PackageExport` from the `QuantumOptimization`` context but are **not** in `PacletInfo` `Symbols`, so a bare `Needs["Wolfram`QuantumFramework`"]` leaves them in `` Global` `` (unevaluated). `IBMJob`/`IBMJobSubmit`/`QuantumQASM`/`QiskitTarget`, by contrast, *are* in the framework context and resolve normally.

Indexed by function/topic.

## Format
- **The mistake** — concrete code or pattern.
- **Why it fails / why it's bad** — root cause with file:line.
- **The right way** — what to do instead.
- **Date observed**.

## Seeded entries (from synthesis observations)

### ⚰️ RESOLVED 2026-05-04 — ~~`QuantumSimilarity` — undocumented public symbol~~
- **The mistake**: Calling `QuantumSimilarity[qs1, qs2]` based on a guess of the signature.
- **Why**: The symbol is exported in `PacletInfo.wl:42` but has **no `Usage` message** in `Kernel/Usage.m` and **no reference page** in `Documentation/.../Symbols/`. Behavior is not contractually documented.
- **The right way**: Read `Kernel/QuantumDistance.m` (likely co-located) to determine the actual signature before calling. Until then, prefer `QuantumDistance` with an explicit measure.
- **Date**: 2026-04-29.
- **Resolution**: fixed in commit `edb226fb` (2026-05-04). `Usage.m` now has `QuantumSimilarity::usage` (alphabetically between `QuantumShortcut` and `QuantumStateEstimate`); `Documentation/.../Symbols/QuantumSimilarity.nb` created from the `QuantumDistance.nb` template, with the Notes section enumerating the 4 distance→similarity transformation rules from `QuantumDistance.m:63-79`.

### ⚰️ RESOLVED 2026-05-04 — ~~`QuatumPartialTranspose.m` filename typo~~
- **The mistake**: Spelling the symbol `QuantumPartialTranspose` (correct English) when calling.
- **Why**: The kernel file is `Kernel/QuatumPartialTranspose.m` (missing 'n'). The exported symbol may be misspelled to match the file.
- **The right way**: Verify with ``?Wolfram`QuantumFramework`*PartialTranspose*`` before using.
- **Date**: 2026-04-29.
- **Resolution**: file renamed to `Kernel/QuantumPartialTranspose.m` on 2026-05-04 (commit `f6124ea8`). The exported symbol was always correctly spelled `QuantumPartialTranspose`; only the file was wrong. Two stale notebook file-manager references (`QuantumFramework/ResourceDefinition.nb`, `Notebooks/ResourceDefinition_Jul2.nb`) updated to the new path.

### `MatrixInverse` silent pseudo-inverse fallback
- **The mistake**: Calling `Inverse[op["Matrix"]]` or trusting `MatrixInverse[m]` and assuming the result is a true inverse.
- **Why**: `Utilities.m:303` — `Inverse` is `Quiet`-tried via `Check`, then falls back to `PseudoInverse` on failure. **No warning emitted.**
- **The right way**: For potentially singular matrices, call `MatrixRank` first or use `Inverse[m, ZeroTest -> ...]` explicitly to detect singularity.
- **Date**: 2026-04-29.

### `QuantumPartialTrace` densifies sparse density matrices
- **The mistake**: Building a sparse N-qubit density matrix, then calling `QuantumPartialTrace` and expecting the output to stay sparse / fit in memory.
- **Why**: `MatrixPartialTrace` (`Utilities.m:132-137`) uses `ArrayReshape` + `TensorContract` — always dense. The 2^N × 2^N intermediate is materialized.
- **The right way**: For large N, partial-trace via tensor-network paths (`QuantumMPS`) or compute the reduced state from the bipartition directly. *(Concrete recipe to be filled from `QuantumCircuitOperator/TensorNetwork.m` read.)*
- **Date**: 2026-04-29.

### `eigensystem` is generic (no Hermitian path)
- **The mistake**: Expecting `qs["Eigenvalues"]` on a Hermitian density matrix to use the fast Hermitian eigensolver.
- **Why**: `Utilities.m:151-181` — calls plain `Eigensystem` with `Simplify` and `Chop` wrappers. **No `Method -> "Hermitian"` is ever set.**
- **The right way**: For large pure-numeric Hermitian inputs where this matters, drop to `Eigensystem[qs["DensityMatrix"], Method -> "Hermitian"]`.
- **Date**: 2026-04-29.

### Secondary symbols need full context
- **The mistake**: Calling `QuantumMPS[...]` directly after ``Needs["Wolfram`QuantumFramework`"]``.
- **Why**: `QuantumMPS`, `QuantumCircuitMultiwayGraph`, `QuantumShortcut`, `QuantumWignerMICTransform` have `Usage` messages but are NOT in `PacletInfo.wl`'s `Symbols` list. They may not be auto-exported into the user namespace.
- **The right way**: Verify symbol resolution after loading; use ``Wolfram`QuantumFramework`QuantumMPS`` if unqualified resolution fails.
- **Date**: 2026-04-29.

### `SparseArrayFlatten` densifies for ≤11 dims
- **The mistake**: Assuming sparse representation survives a flatten in any QF tensor pipeline.
- **Why**: `Utilities.m:275-301` — the custom sparse-preserving flatten only triggers when `Length[dims] > 11`. Below that, falls through to plain `Flatten` which produces a dense list.
- **The right way**: For workflows where sparsity is critical, explicitly stay in `SparseArray` form and avoid `Flatten` until the very end.
- **Date**: 2026-04-29.

### Property cache has no auto-invalidation
- **The mistake**: Building `qs = QuantumState[v, qb]`, querying `qs["VonNeumannEntropy"]`, substituting into `v`, and expecting a fresh result.
- **Why**: `$QuantumFrameworkPropCache = True` (`QuantumFramework.m:15`) caches results process-wide. There's no auto-invalidation when the underlying expression changes via substitution; only structural identity matters.
- **The right way**: Pin intermediate states as named symbols (`qs2 = qs /. v -> v2`) so the cache key changes; or set ``Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = False`` for the duration of the substitution-heavy code.
- **Date**: 2026-04-29.

### `qs["VonNeumannEntropy"]` and `qs["Purity"]` silently fail above 10 qubits, but only for MIXED states
- **The mistake**: Computing `qs["VonNeumannEntropy"]` or `qs["Purity"]` on a **mixed** state with `Dimension >= 1024` (10+ qubits) and getting `Indeterminate` with no warning. A **pure** state of any size is fine (see below).
- **Why**: `QuantumState/Properties.m:460, :480` — `ConfirmAssert[qs["Dimension"] < 2^10]` inside an `Enclose`, but on the **non-pure branch only**. Both handlers fast-path pure states first: entropy is `If[qs["PureStateQ"], 0, ConfirmAssert[…]; …]` and purity is `If[qs["StateType"] === "Vector", 1, ConfirmAssert[…]; …]`. So a pure state returns `0` / `1` regardless of size; only a mixed state of `Dimension >= 1024` trips the assert, and `Enclose` then returns `Indeterminate`. **No warning emitted.** Verified live: `QuantumState["GHZ"[10]]` gives entropy `0` / purity `1`; a mixed dimension-1024 state gives `Indeterminate`; a mixed 9-qubit (dimension-512) state computes normally.
- **The right way**: This only bites mixed states `>= 10` qubits. Reduce via `QuantumPartialTrace[qs, complement]` first (entropy of a small subsystem is fine), or route through `QuantumMPS`. As a purity proxy without the limit, `qs["TraceNorm"]` (`Properties.m:375`) computes `Total[SingularValueList[DensityMatrix]]`.
- **Date**: 2026-04-29.

### `qs["Eigenvalues"]` is NOT the same as the `eigenvalues` helper
- **The mistake**: Expecting `qs["Eigenvalues"]` to use the wrapped `eigenvalues` helper from `Utilities.m` (with `Chop`, `Sort`, `Normalize` options).
- **Why**: `QuantumState/Properties.m:438` calls bare `Eigenvalues[qs["DensityMatrix"]]`. But `Properties.m:440-442` for `"Eigenvectors"` and `"Eigensystem"` DO use the wrapped helpers (`eigenvectors`, `eigensystem`). **Inconsistent.**
- **The right way**: For wrapped behavior with options, call `qs["Eigensystem", "Sort" -> True, "Normalize" -> True]` and take `First`. Or use the helper directly: ``eigenvalues[qs["DensityMatrix"], opts]`` after `Needs["Wolfram`QuantumFramework`PackageScope`"]`.
- **Date**: 2026-04-29.

### `qs["DensityMatrix"]` on pure states allocates N² memory
- **The mistake**: Calling `qs["DensityMatrix"]` on a pure-state vector and not realizing the cost.
- **Why**: `QuantumState/Properties.m:415` builds it as `KroneckerProduct[stateVector, Conjugate[stateVector]]`. A 2^N-entry vector becomes a 2^(2N) matrix.
- **The right way**: For pure-state computations, use `qs["StateVector"]` directly. Materialize the density matrix only when a property genuinely needs it.
- **Date**: 2026-04-29.

### `qs["StateVector"]` on density-matrix-given pure states triggers eigendecomp
- **The mistake**: Calling `qs["StateVector"]` and assuming it's a free read of the underlying state.
- **Why**: `QuantumState/Properties.m:130-145` — if the state is stored as a density matrix but happens to be pure, QF eigendecomposes the density matrix to recover the vector. O(d³) work.
- **The right way**: Build pure states from vectors (`QuantumState[vector, basis]`), not from density matrices. If you must work with a density-matrix-given pure state, pin the result: `vec = qs["StateVector"]` once.
- **Date**: 2026-04-29.

### `qs["Type"]` is O(d³) due to `PositiveSemidefiniteMatrixQ`
- **The mistake**: Calling `qs["Type"]` repeatedly (e.g. inside `Cases` or a loop predicate).
- **Why**: `QuantumState/Properties.m:481-497` — branches through `PositiveSemidefiniteMatrixQ[N @ qs["DensityMatrix"]]`, which runs an eigendecomposition. **The cache helps for non-parametric states**, but parametric states bypass cache (see next entry).
- **The right way**: Pin: `Module[{type = qs["Type"]}, ...]` if branching on Type multiple times in a non-cached path.
- **Date**: 2026-04-29.

### Property cache skips parametric states
- **The mistake**: Expecting `qs["VonNeumannEntropy"]` to be cached and instant on repeated calls when `qs` has free parameters.
- **Why**: `QuantumState/Properties.m:51-62` — caching is gated on `qs["Basis"]["ParameterArity"] == 0`. **Parametric states never cache, regardless of `$QuantumFrameworkPropCache`.**
- **The right way**: Substitute parameters once before extraction: `qs2 = qs[{p -> 0.5}]`. The substituted state is non-parametric and will cache normally.
- **Date**: 2026-04-29.

### `qs["BlochPlot"]`, `BlochSphericalCoordinates`, `BlochCartesianCoordinates` only work for `Dimension == 2`
- **The mistake**: Calling Bloch-related properties on a multi-qubit state.
- **Why**: `QuantumState/Properties.m:929, 962, 987` — guarded by `qs["Dimension"] == 2`. The property catalog (`:65-66`) also strips `Bloch*` from properties for `Dimension != 2`.
- **The right way**: Use `qs["BlochVector"]` (`Properties.m:953`) for general dimensions — it computes `Tr[ρ G_i]` against `GellMannMatrices[d]` for any d.
- **Date**: 2026-04-29.

### Display icon falls back to placeholder above `Dimension == 512`
- **The mistake**: Surprised that a `QuantumState[...]` summary box doesn't show the matrix plot for a moderately-sized state.
- **Why**: `QuantumState/Formatting.m:23` — `ComplexArrayPlot` icon only renders if `qs["Dimension"] < 2^9 = 512`. Above that, falls back to `$SparseArrayBox` placeholder. Purity/Entropy in summary use `TimeConstrained[..., 1]` and become `$Aborted` if 1 sec isn't enough.
- **The right way**: Just be aware — this is a display-only fallback, not a computational issue. The state itself is fine.
- **Date**: 2026-04-29.

### `qo["UnitaryQ"]` is numeric-only and lies for symbolic unitaries
- **The mistake**: Calling `qo["UnitaryQ"]` on a symbolic operator (e.g., a parametric rotation) and trusting the False result.
- **Why**: `QuantumOperator/Properties.m:445` — `UnitaryMatrixQ[N[qo["Matrix"]]]`. `N` on a symbolic matrix often doesn't simplify enough to satisfy `UnitaryMatrixQ`'s tolerance, so a genuinely unitary `R_y(θ)` can return `False`.
- **The right way**: For symbolic operators, verify by hand: `Simplify[qo["Matrix"] . qo["Dagger"]["Matrix"] - IdentityMatrix[qo["OutputDimension"]]] === 0`. Or check `qo["HermitianQ"]` for Hamiltonians (it IS symbolic-friendly — `HermitianMatrixQ[qo["Matrix"]]` without `N`).
- **Date**: 2026-04-29.

### Control structure of an operator lives in the Label, not the Order
- **The mistake**: Trying to add controls by manipulating `qo["Order"]` or `qo["InputOrder"]`.
- **Why**: `QuantumOperator/Properties.m:127-136` — `ControlOrder1`, `ControlOrder0`, `ControlOrder`, `TargetOrder` are extracted by pattern-matching the operator's `Label`: `Subscript["C", _][control1_, control0_]`. Reordering qudits does NOT add or remove controls.
- **The right way**: Build controlled operators with the named-operator system (e.g., `QuantumOperator["CNOT"]`) or by setting the `"Label"` option to `Subscript["C", innerLabel][{control1Qubits}, {control0Qubits}]`. The Order specifies WHICH qudits — Label specifies WHICH ARE CONTROL.
- **Date**: 2026-04-29.

### Operator composition routes through `QuantumCircuitOperator` (overhead)
- **The mistake**: Assuming `qo1[qo2]` is a direct tensor contraction.
- **Why**: `QuantumOperator/QuantumOperator.m:346` — `qo1[qo2] := QuantumOperator @ QuantumCircuitOperator[{qo2, qo1}]`. The chain wraps in a circuit then unwraps. There's a disabled direct EinsteinSummation path (`:348-435`) commented out. Same overhead applies to `qo[qs]`, `qo[qmo]`, `qo[qc]`, `qo[qm]`.
- **The right way**: For chains of >5 operators, build a `QuantumCircuitOperator[{op1, op2, ...}]` once and apply at the end — that's the single circuit traversal, not N nested ones. For a single multiply, the overhead is small but real.
- **Date**: 2026-04-29.

### `qo["PauliDecompose"]` has fast/slow paths
- **The mistake**: Calling `qo["PauliDecompose"]` on a high-dimension qudit operator (say a 6-level qudit operator) and waiting forever.
- **Why**: `QuantumOperator/Properties.m:626-649` — for `Dimensions == {2..}` (qubits), uses a fast tensor-reshape recursion. For non-qubit dimensions, builds a system of equations and calls `SolveValues` — exponential in size for symbolic input.
- **The right way**: `PauliDecompose` is intended for qubit operators. For qudits, decompose into `GellMannMatrices` directly: `Tr[op . #] & /@ GellMannMatrices[d]`.
- **Date**: 2026-04-29.

### Operator equality silently coerces to matrix form
- **The mistake**: Comparing two `QuantumOperator`s and assuming the picture/basis matters when the matrix entries match.
- **Why**: `QuantumOperator/QuantumOperator.m:597` — `Equal` first checks `Picture` equality, then `Chop @ SetPrecisionNumeric @ SparseArrayFlatten @ Sort["MatrixRepresentation"]`. **Promotes vector form to matrix via `ToMatrix` if any operand is in matrix form.**
- **The right way**: For exact picture/basis equality (not matrix equivalence), compare `qo["Picture"]`, `qo["Basis"]` and `qo["Order"]` separately. `Equal` is matrix-equivalence, not structural equivalence.
- **Date**: 2026-04-29.

### ⚰️ RESOLVED 2026-04-29 — ~~KERNEL BUG: `QuantumOperator["HeisenbergWeyl"]` returns an undefined symbol~~
- **The mistake**: Calling `QuantumOperator["HeisenbergWeyl"]` and getting `QuantumOperartor[{"HeisenbergWeyl", 2}, ...]` back, unevaluated.
- **Why**: `QuantumOperator/NamedOperators.m:855` has a typo — `QuantumOperartor` (missing letter) instead of `QuantumOperator`. Real bug at `cbbe9368`.
- **The right way**: Until the typo is fixed, use the explicit form: `QuantumOperator[{"HeisenbergWeyl", p, i, a}, order]` (`:857`). Or build the matrix directly: `Total @ Table[Exp[I 2 Pi a l / d] KroneckerProduct[UnitVector[d, Mod[i + l, d] + 1], UnitVector[d, l + 1]], {l, 0, d - 1}]`.
- **Date**: 2026-04-29.
- **Resolution**: fixed in commit `c74f3baf` (2026-04-29) — `QuantumOperartor` → `QuantumOperator` at NamedOperators.m:855.

### ⚰️ RESOLVED 2026-04-29 — ~~KERNEL BUG: One `{"Controlled0", non-operator-args, control0}, target` dispatch path is broken~~
- **The mistake**: Building a controlled-0 operator from a non-operator name + target: `QuantumOperator[{"Controlled0", "X", {1}}, {2}]`.
- **Why**: `QuantumOperator/NamedOperators.m:379` has a typo — `QuantumOperartorQ` instead of `QuantumOperatorQ`. The `ConfirmBy` fails because the predicate symbol is undefined.
- **The right way**: Pass an explicit `QuantumOperator` instead of a name + args: `QuantumOperator[{"Controlled0", QuantumOperator["X", {2}], {1}}]`. This skips the typoed dispatch.
- **Date**: 2026-04-29.
- **Resolution**: fixed in commit `c74f3baf` (2026-04-29) — `QuantumOperartorQ` → `QuantumOperatorQ` at NamedOperators.m:379.

### `QuantumOperator["IXYZ"]` is a Pauli-string shortcut — but ONLY for `{I, X, Y, Z, H, S, T, V, P}`
- **The mistake**: Expecting `QuantumOperator["XYZ"]` to mean a single-qubit named gate, or `QuantumOperator["CNOT"]` to be decoded as a string chain.
- **Why**: `NamedOperators.m:963-966` decodes strings composed *entirely* of `{"I", "X", "Y", "Z", "H", "S", "T", "V", "P"}` as a tensor product. `"CNOT"` is NOT in that set, so it falls through to the named-operator dispatch and yields a CNOT gate. But `"IXX"` becomes `I⊗X⊗X` and `"HHH"` becomes `H⊗H⊗H`.
- **The right way**: For tensor products of Paulis/Hadamards, the string shortcut works. For named gates (CNOT, Toffoli, Fredkin, etc.), use the full name (the dispatch falls through correctly). Don't mix concepts.
- **Date**: 2026-04-29.

### `QuantumOperator["H"[d]]` is the d-dim **Fourier**, not a textbook block-Hadamard
- **The mistake**: Expecting `QuantumOperator["H"[d]]` to be a Hadamard-style operator with structure for `d > 2`.
- **Why**: `NamedOperators.m:613-620` — `QuantumOperator[("Hadamard" | "H")[dim], order]` rewrites to `QuantumOperator["Fourier"[dim] -> order]` with `"Label" -> "H"`. Hadamard IS the discrete Fourier transform. For d = 2 they're equivalent (`F_2 = H`); for d > 2 it's the generalized DFT.
- **The right way**: For a qudit DFT, `"H"` is correct. For a literal `H ⊗ ... ⊗ I`, build the matrix explicitly or use the string-chain shortcut `QuantumOperator["HII..."]`.
- **Date**: 2026-04-29 (syntax refreshed for v2.0 named constructors, 2026-05-22).

### Operator names are auto-uppercased — `QuantumOperator["cnot"]` works
- **The mistake**: Doubting whether lowercase names resolve.
- **Why**: `NamedOperators.m:970-976` — auto-uppercase mapping `$upperCasesOperatorNames`. Lowercase names map to canonical `$QuantumOperatorNames` form. So `"cnot" → "CNOT"`, `"hadamard" → "Hadamard"`, etc.
- **The right way**: Both forms work; for clarity, prefer the canonical case from `$QuantumOperatorNames`.
- **Date**: 2026-04-29.

### Invalid named operators return `Failure`, not error
- **The mistake**: Wrapping a possibly-invalid named-operator call in `Quiet`/`Check` expecting a thrown error.
- **Why**: `NamedOperators.m:993-995` — invalid `QuantumOperator["name"[args]]` (and the bare `QuantumOperator["name"]`) returns `Failure["InvalidName", ...]` rather than throwing. Same dispatch holds for `QuantumState`, `QuantumChannel`, `QuantumCircuitOperator`, `QuantumBasis`, `QuditBasis`, `QuantumMeasurementOperator` after `fa0d0113`. Caller must check with `FailureQ` or pattern-match.
- **The right way**: `Enclose @ With[{op = ConfirmBy[QuantumOperator["name"[args]], QuantumOperatorQ]}, ...]`. Use `Confirm` to convert `Failure` into a caught exception. The legacy `QuantumOperator[{"name", args}]` form returns the same `Failure["InvalidName", ...]` now (it's no longer a valid surface).
- **Date**: 2026-04-29 (refreshed 2026-05-22 for v2.0 named-constructor switch).

### `QuditBasis["PauliX"]` is the X-eigenbasis, NOT the σ_x matrix
- **The mistake**: Calling `QuditBasis["PauliX"]` and expecting `{{0, 1}, {1, 0}}` as the basis matrix.
- **Why**: `QuditBasis/NamedBases.m:67-83` — `PauliX/Y/Z` builds the eigensystem of `pauliMatrix[i, dim]` with `Sort -> True`. So `QuditBasis["PauliX"]` returns the basis `{|+⟩, |−⟩}` with eigenvalue-labeled names. This is what `Properties.m:208` then KroneckerProducts as the basis matrix.
- **The right way**: For the σ_x matrix as a *named* basis element, use `QuditBasis["Pauli"]` (gives 4-element {I, X, Y, Z} basis), or pull `pauliMatrix[1, 2]` from `Utilities.m`. For the X-eigenbasis, `"PauliX"` is correct.
- **Date**: 2026-04-29.

### `QuditBasis["Pauli"]` is a 4-element basis (σ_0..σ_3), NOT a single qubit
- **The mistake**: Treating `QuditBasis["Pauli"]` as if it were a 2-element single-qubit basis.
- **Why**: `QuditBasis/NamedBases.m:140-147` — `"Pauli"` returns `AssociationThread[{σ_0, σ_1, σ_2, σ_3}, {I, X, Y, Z matrices}]`. **4 elements**, in 2×2 matrix space. Used for operator decomposition / process tomography, not state representation.
- **The right way**: For a single-qubit basis, use `"PauliX" | "PauliY" | "PauliZ"` (each 2-element eigenbasis) or `"Computational"` (default). Reserve `"Pauli"` for tomography/superoperator work.
- **Date**: 2026-04-29.

### `QuditBasis["XYZ"]` string-chain accepts only `{I, X, Y, Z}` — narrower than QuantumOperator's
- **The mistake**: Trying `QuditBasis["XH"]` or `QuditBasis["IS"]` and getting unexpected results.
- **Why**: `QuditBasis/NamedBases.m:302-308` — string decoding gates on `ContainsOnly[chars, {"I", "X", "Y", "Z"}]`. **Strictly stricter** than `QuantumOperator`'s 9-char alphabet (`I, X, Y, Z, H, S, T, V, P`). `QuditBasis["XH"]` falls through to single-name dispatch, doesn't decode as a chain.
- **The right way**: For Pauli-only chains as a basis, the shortcut works. For mixed-gate chains, build with `QuantumTensorProduct[QuditBasis["X"], QuditBasis[QuantumOperator["H"]["Output"]], ...]` or build operators directly.
- **Date**: 2026-04-29.

### Two independent QF caches; both global, both can grow unboundedly
- **The mistake**: Assuming `$QuantumFrameworkPropCache = False` clears all QF caching.
- **Why**: There are TWO global caches:
  1. **`$QuantumFrameworkPropCache`** (`QuantumFramework.m:15`) — caches *property* lookups on `QuantumState`/`QuantumOperator`/`QuantumBasis`/`QuditBasis`. Conditional on parameter-arity for state.
  2. **`$QuditBasisCache`** (`QuditBasis/NamedBases.m:27-33`) — caches *named-basis construction* via `Lookup[$QuditBasisCache, Key[{args}], ...]`. Bypassed only for `Random*MIC`. Independent of `$QuantumFrameworkPropCache`.
- **The right way**: Long sessions should periodically clear `Wolfram\`QuantumFramework\`PackageScope\`$QuditBasisCache = <||>` to release memory. Setting `$QuantumFrameworkPropCache = False` doesn't touch the construction cache.
- **Date**: 2026-04-29.

### `qb["Dimension"]` is `OutputDimension × InputDimension`, not just one side
- **The mistake**: Building a 4-qudit basis and being surprised that `qb["Dimension"]` is 16 (4 output × 4 input), not 4.
- **Why**: `QuantumBasis/Properties.m:167` — `Dimension = OutputDimension × InputDimension`. For a state-space (input dim 1) basis it equals OutputDimension. For an operator-space basis (equal in/out), it's the square.
- **The right way**: For "the dimension I think of as Hilbert space," use `qb["OutputDimension"]`. For matrix-element count, `qb["Dimension"]`. See `qb["MatrixDimensions"]` for `{ElementDimension, Dimension}` shape.
- **Date**: 2026-04-29.

### `QuantumBasis` and `QuditBasis` are `NHoldAll` — `N[basis, prec]` only works through limited paths
- **The mistake**: Calling `N[qb, 30]` and seeing nothing change.
- **Why**: `QuantumBasis/QuantumBasis.m:200` and `QuditBasis/QuditBasis.m:157` set `NHoldAll`. The custom `N` upvalue (`QuantumBasis/QuantumBasis.m:196`) only fires if not all representations are already inexact-numeric.
- **The right way**: For numeric basis evaluation, build with explicit numeric matrices, or use `qb["Output"]` and `N` the underlying tensors via `MapAt`. Don't expect bare `N[qb, 30]` to coerce.
- **Date**: 2026-04-29.

### `qc["TracePreservingQ"]` returns `Undefined` for symbolic channels
- **The mistake**: Building a parametric `QuantumChannel` and asking `qc["TracePreservingQ"]` to verify CPTP.
- **Why**: `QuantumChannel/Properties.m:55-57` — `If[MatrixQ[m, NumericQ], m == IdentityMatrix[...], Undefined]`. **Symbolic matrices return `Undefined`, not `True`/`False`.**
- **The right way**: For symbolic verification, expand and `Simplify[qc["Trace"]["MatrixRepresentation"] - IdentityMatrix[d]] === ConstantArray[0, ...]`. Or substitute parameters before checking.
- **Date**: 2026-04-29.

### Channel auxiliary qudits are encoded as non-positive integers in `Order`
- **The mistake**: Trying to access channel environment qudits via positive indices, or assuming `qc["FullOutputOrder"]` returns only positive integers.
- **Why**: `QuantumChannel/QuantumChannel.m:9-10` — channel validity requires `OutputQudits >= InputQudits`, with the extra outputs given non-positive orders. `Properties.m:45` — `qc["TraceOrder"] := Select[qc["FullOutputOrder"], NonPositive]`. The 0 / negative orders ARE the environment.
- **The right way**: Use `qc["TraceOrder"]`, `qc["TraceQudits"]`, `qc["TraceDimensions"]` to access environment qudits explicitly. Don't filter on positive only.
- **Date**: 2026-04-29.

### `qc["EvolutionChannel"]` treats Kraus operators as Hamiltonian generators (surprising)
- **The mistake**: Calling `qc["EvolutionChannel"]` on a Kraus-form channel and getting back `Exp[I K_1] ⊗ Exp[I K_2] ⊗ ...` instead of an evolution-equivalent representation.
- **Why**: `QuantumChannel/Properties.m:59-61` — `(I #)["EvolutionOperator"] & /@ qc["UnstackOutput", 1]`. Each Kraus operator is treated as a Hamiltonian: `K → exp(I·K·t)`. This is a specific physical model, not a general Kraus-to-Lindblad conversion.
- **The right way**: For genuine Kraus-to-Lindblad evolution, use `QuantumEvolve[H, {L1, L2}, ...]` directly with the Lindblad form. `EvolutionChannel` is a different operation suited for specific use cases (gate-as-channel exponentiation).
- **Date**: 2026-04-29.

### Channel composition routes through `QuantumCircuitOperator`
- **The mistake**: Assuming `qc1[qc2]` is a direct channel composition / Stinespring concatenation.
- **Why**: `QuantumChannel/QuantumChannel.m:50` — `(qc_QuantumChannel)[args___] := QuantumCircuitOperator[qc][args]`. The direct channel-on-channel implementation (lines 52-67) is **commented out**. Like operator composition, channels build a circuit and unwrap.
- **The right way**: For correctness this works. For performance with many channels, build `QuantumCircuitOperator[{c1, c2, c3, ...}]` once and apply at the end.
- **Date**: 2026-04-29.

### Default `qco[qs]` is TensorNetwork, not Schrödinger
- **The mistake**: Reasoning about `qco[qs]` as a sequence of `qs |-> gate1[qs] |-> gate2[qs] |-> ...` Schrödinger applications.
- **Why**: `QuantumCircuitOperator/QuantumCircuitOperator.m:87-102` — `quantumCircuitApply` defaults to `Automatic` which routes to `TensorNetworkApply[qco["Flatten"], qs, ...]`. The Schrödinger fold path is opt-in via `Method -> "Schrodinger"`.
- **The right way**: For symbolic, hand-traceable evolution, pass `Method -> "Schrodinger"`. For numeric large-circuit performance, the default tensor-network path is correct. Be explicit about which you want — cost characteristics differ.
- **Date**: 2026-04-29.

### `qco[qs, Method -> "Schrodinger"]` calls `FullSimplify` after every gate
- **The mistake**: Selecting `Method -> "Schrodinger"` for a symbolic 30-gate circuit and watching it grind.
- **Why**: `QuantumCircuitOperator/QuantumCircuitOperator.m:92` — `Fold[(n++; FullSimplify @ #2[#1, "Computational" -> False]) &, qs, qco["NormalOperators"]]`. **Every intermediate state is `FullSimplify`-d.** For symbolic states with deep parameter dependencies, this can be the dominant cost.
- **The right way**: For symbolic-but-shallow circuits, this is fine and gives clean output. For deep symbolic circuits, use `Method -> "TensorNetwork"` (default) and apply `FullSimplify` once at the end. Or wrap parameters before applying.
- **Date**: 2026-04-29.

### TensorNetworkCompile silently auto-bends for matrix operators / trace orders
- **The mistake**: Building a circuit with one mixed-state QuantumChannel among many pure operators and being surprised by 2× memory.
- **Why**: `QuantumCircuitOperator/TensorNetwork.m:203` — `bendQ = (AnyTrue[circuit["Operators"], #["MatrixQ"] &] || circuit["TraceQudits"] > 0) && ! phaseSpaceQ`. When `bendQ` is True, the circuit is replaced with `circuit["Bend"]` (`:207`), which doubles the operator count via the Stinespring dilation. **Hidden cost.**
- **The right way**: For mostly-pure circuits, isolate the mixed/channel operations and apply at the end. For genuinely mixed-state evolution, this is the correct cost; just be aware it's automatic.
- **Date**: 2026-04-29.

### TensorNetworkCompile inserts basis-conversion tensors for non-computational bases
- **The mistake**: Building operators in a custom basis (e.g., Pauli eigenbasis) and being surprised by extra tensor nodes in `qco["TensorNetworkGraph"]`.
- **Why**: `QuantumCircuitOperator/TensorNetwork.m:30-35, :100-105` — when `"Computational" -> True` (default) and an operator's `Input` or `Output` basis isn't computational, QF prepends `MatrixInverse[Input["ReducedMatrix"]]` and appends `Output["ReducedMatrix"]` as extra `QuantumOperator` nodes. Node count grows.
- **The right way**: For symbolic-basis circuits where you want to inspect the network, pass `"Computational" -> False`. For numeric performance, the conversion is usually beneficial.
- **Date**: 2026-04-29.

### Circuit cache **always bypasses** 17 properties
- **The mistake**: Setting `$QuantumFrameworkPropCache = True` and expecting `qco["Operators"]` or `qco["Flatten"]` to be cached on second access.
- **Why**: `QuantumCircuitOperator/Properties.m:28-31` — `$QuantumCircuitPreventCache = {"Association", "Elements", "Operators", "FullOperators", "FullElements", "Options", "Diagram", "Icon", "Qiskit", "QiskitCircuit", "QuantumOperator", "Flatten", "Double", "Bend", "DiscardExtraQudits", "Picture", "Parameters", "ParameterArity"}`. The cache check explicitly excludes these 17 (`:36`).
- **The right way**: For repeated access to these in tight loops, pin manually: `Module[{ops = qco["Operators"]}, ...]`.
- **Date**: 2026-04-29.

### `qco["Bend"]` doubles the operator count (not just superoperator wrap)
- **The mistake**: Calling `qco["Bend"]` and assuming it just changes the basis to superoperator form.
- **Why**: `QuantumCircuitOperator/Properties.m:389-399` — for each non-matrix operator, splices in `{op, op["Conjugate"]["Shift", qco["Width"]]}`. So a 10-gate pure circuit becomes a 20-gate bent circuit with a shifted-width copy. Channels and measurement operators get wrapped via Bend instead. **Operator count doubles.**
- **The right way**: Use `Bend` deliberately when you need the duplicated form (for instance, to compute `E[ρ]` via state-on-bent-operator pattern). For ordinary application of pure circuits, don't bend.
- **Date**: 2026-04-29.

### `ToQuESTLink` only handles 2-dim qudits
- **The mistake**: Building a qutrit circuit and asking `qco[qs, Method -> "QuEST"]` to evaluate it.
- **Why**: `QuantumCircuitOperator/QuEST.m:58` — `ToQuESTLink[qo_QuantumOperator] /; MatchQ[qo["Dimensions"], {2 ..}]`. Non-qubit operators don't match the pattern; the catch-all `:86` returns `Failure["Unknown", ...]`. Unrecognized labels in `:64-83` return `Missing["NotImplemented", label]`.
- **The right way**: Use QuEST only for qubit circuits. For higher-dim qudits use the default tensor-network path or `"Schrodinger"`.
- **Date**: 2026-04-29.

### Display icon is fake for circuits with > 32 gates
- **The mistake**: Looking at the summary box for a 100-gate circuit and seeing what looks like a small Fourier circuit, then assuming the circuit IS a Fourier transform.
- **Why**: `QuantumCircuitOperator/Formatting.m:14` — when `qco["GateCount"] > 32`, the icon is replaced with `QuantumCircuitOperator[{"Fourier"[3]}]["Icon", "GateBackgroundStyle" -> _ -> LightGray, "GateBoundaryStyle" -> _ -> Gray]`. **It's a generic placeholder, not the real circuit.**
- **The right way**: Use `qco["Diagram"]` or `qco["Diagram", ImageSize -> Large]` to see the actual structure for any circuit. The summary icon is purely a visual cue.
- **Date**: 2026-04-29.

### `QuantumCircuitOperator["Bell"]` is `{H, CNOT}`, not the Bell state
- **The mistake**: Calling `QuantumCircuitOperator["Bell"]` and treating it as if it were `QuantumState["Bell"]`.
- **Why**: `QuantumCircuitOperator/NamedCircuits.m:459-460` — `QuantumCircuitOperator["Bell"[n]]` returns `QuantumCircuitOperator[{"H", "CNOT"}]["Shift", n]`. It's the **preparation circuit**, not the state. (List form `{"H", "CNOT"}` is the gate-sequence form, preserved by the v2.0 refactor.)
- **The right way**: For the Bell state itself, `QuantumState["Bell"]` (`QuantumState/NamedStates.m:130`). For the prep circuit applied to `|00⟩`, `QuantumCircuitOperator["Bell"][QuantumState["00"]]` produces the same state as `QuantumState["Bell"]` (in computational basis).
- **Date**: 2026-04-29 (syntax refreshed for v2.0, 2026-05-22).

### Toffoli/Fredkin are 14+ gates, not multi-controlled
- **The mistake**: Building `QuantumCircuitOperator["Toffoli"]` and assuming `qco["GateCount"] == 1`.
- **Why**: `QuantumCircuitOperator/NamedCircuits.m:462-477` — Toffoli is a hardcoded sequence: `{H -> 3, CNOT -> {2,3}, T† -> 3, CNOT -> {1,3}, T -> 3, CNOT -> {2,3}, T† -> 3, CNOT -> {1,3}, T -> 2, T -> 3, H -> 3, CNOT -> {1,2}, T -> 1, T† -> 2, CNOT -> {1,2}}`. 15 gates. Fredkin is similar (18 gates).
- **The right way**: For the multi-controlled abstract operator, use `QuantumOperator["Toffoli"]` (single 3-qubit operator). For an executable decomposition, the circuit form. Choose deliberately.
- **Date**: 2026-04-29.

### Trotterization order > 2 explodes via Suzuki recursion
- **The mistake**: Calling `QuantumCircuitOperator[{"Trotterization", {ops...}, order, reps}]` with order=4, 6, 8 expecting linear growth.
- **Why**: `QuantumCircuitOperator/NamedCircuits.m:606-614` — `trotterCoeffs[l, n_? EvenQ, c]` uses Suzuki recursion: 5 sub-blocks of order n-2 per recursion level. Gate count is `~5^(n/2 - 1) × Length[ops] × reps`. Order 4 ≈ 5×, order 6 ≈ 25×, order 8 ≈ 125×.
- **The right way**: Order 1 (Lie-Trotter) and order 2 (Strang) are linear in `Length[ops]`. For higher accuracy without explosion, increase `reps` instead of `order`.
- **Date**: 2026-04-29.

### `PhaseEstimation` is qubit-only
- **The mistake**: Building `QuantumCircuitOperator[{"PhaseEstimation", QuantumOperator[..., qutrit_basis]}]` and waiting for it to work.
- **Why**: `QuantumCircuitOperator/NamedCircuits.m:556, :575` — guarded by `MatchQ[op["OutputDimensions"], {2 ..}]`. Non-qubit operators don't match the dispatch rule and fall through.
- **The right way**: For qudit phase estimation, decompose to qubits manually or build a custom circuit with controlled qudit operators + qudit Fourier.
- **Date**: 2026-04-29.

### `QiskitCircuit` unknown-gate injection collapses to machine (double) precision
- **The mistake**: Sending a high-precision (e.g., 30-digit) symbolic-numeric state through `qco["Qiskit"]` and being surprised by precision loss.
- **Why**: `QuantumCircuitOperator/Qiskit.m:57, :59` — fallback `Unitary` gate is built with `NumericArray[Normal @ N @ arr, "ComplexReal64"]` — **double precision (64-bit complex)**. `N[arr]` is already `MachinePrecision`, so an arbitrary-precision (e.g. 30-digit) input is collapsed to ~15-16 decimal digits.
- **The right way**: For arbitrary-precision needs (`WorkingPrecision` above machine), stay in the WL kernel (Schrödinger or TensorNetwork methods). For ordinary numeric work the Qiskit roundtrip is now full double precision, so there is no longer a single-precision floor to worry about.
- **Date**: 2026-04-29; corrected 2026-06-05 (was `ComplexReal32`/32-bit; live kernel at HEAD `fb16ffc2` uses `"ComplexReal64"`).

### `QiskitCircuit["QuantumOperator"]` materializes the full matrix
- **The mistake**: Calling `qiskitCircuit["QuantumOperator"]` on a 20-qubit Qiskit circuit and waiting for hours.
- **Why**: `QuantumCircuitOperator/Qiskit.m:518-521` — calls Python `Operator(qc).data` which builds the full 2^N × 2^N unitary. Same complexity as Qiskit's own statevector simulator memory limits.
- **The right way**: For large circuits, convert via `qiskitCircuit["QuantumCircuitOperator"]` (returns the gate-list form, no matrix), then apply with the tensor-network method. Or work with the Qiskit circuit directly.
- **Date**: 2026-04-29.

### IBM auth uses `SystemCredential` silently
- **The mistake**: Running `qco[Method -> "Qiskit", "Provider" -> "IBM"]` and getting an opaque Python error about credentials.
- **Why**: `QuantumCircuitOperator/Qiskit.m:332-336` — auth token is looked up via `SystemCredential["IBMQuantumPlatform_APIKEY"]`. If neither token is passed nor stored in the system credential vault, the Python side fails (the error is unhelpful).
- **The right way**: Set up once via ``SystemCredential["IBMQuantumPlatform_APIKEY"] = "..."``, or pass `"Provider" -> {"IBMProvider", "Token" -> "..."}` explicitly.
- **Date**: 2026-04-29.

### Qiskit roundtrip relabels: `X` ↔ `NOT`
- **The mistake**: Building `QuantumOperator["X"]`, exporting via Qiskit, importing back, and being surprised the label is `"NOT"` instead of `"X"`.
- **Why**: `QuantumCircuitOperator/Qiskit.m:241` — `gate_to_QuantumOperator` maps Qiskit's `x` gate to `QuantumOperator['NOT', order]` rather than `'X'`. Asymmetric with the `shortcutToGate` mapping (`:18`) where both `X` and `NOT` map to `XGate`.
- **The right way**: Don't rely on label preservation through Qiskit roundtrip. The matrix and behavior are correct; only the display label may flip. Use `qo["MatrixRepresentation"]` for equivalence checking, not labels.
- **Date**: 2026-04-29.

### `save_statevector` from a Qiskit import is silently dropped
- **The mistake**: Importing a Qiskit circuit that ends with `qc.save_statevector()` and being surprised the QF circuit doesn't have a "save state" gate.
- **Why**: `QuantumCircuitOperator/Qiskit.m:259` — `if gate.name == 'save_statevector': return wl.Nothing`. Drops the instruction silently.
- **The right way**: To capture an intermediate state in QF, use the explicit Schrödinger evaluation or build a `QuantumMeasurementOperator` at the desired position. Qiskit's `save_*` instructions don't roundtrip.
- **Date**: 2026-04-29.

### ⚠️ NOT resolved — `QuantumHamiltonianOperator` subsystem still on `main`
The removal described in earlier passes was done only on the unmerged `qho-decapitation-shim` branch; on `main` (HEAD `bfddf5ed`, verified 2026-05-27) the `Kernel/QuantumHamiltonianOperator/` subdir is still present, so the four pitfalls below remain **live in source**. They are moot in practice only because the subsystem fails to load: `QuantumHamiltonianOperator` resolves to `` Global` `` with 0 DownValues / 0 SubValues (verified in a fresh kernel), so the head stays unevaluated rather than executing buggy code. See `audit/qho-deprecation-audit.md` for the (still-accurate) feature inventory.

#### KERNEL FILENAME: `QuantumHamiltonialOperator.m` (typo, missing 'n')
- **The mistake**: Looking for `QuantumHamiltonianOperator.m` (correct English) and not finding it.
- **Why**: `Kernel/QuantumHamiltonianOperator/QuantumHamiltonialOperator.m` — a filename typo (the `QuatumPartialTranspose.m` typo was since renamed in `f6124ea8`, and the two `QuantumOperartor` typos in NamedOperators.m fixed in `c74f3baf`). The exported symbol `QuantumHamiltonianOperator` is spelled correctly; only the filename is wrong.
- **Status**: still present on `main` at HEAD `bfddf5ed`.

#### `QuantumHamiltonialOperator.m` has a commented-out `Package` declaration
- **The mistake**: Editing this file and being confused why standard `PackageScope`/`PackageExport` declarations seem to behave inconsistently.
- **Why**: Line 1 is `(* Package["Wolfram\`QuantumFramework\`"] *)` — commented out, while other files in the subdir start with `Package[...]` directly. As a result `PackageExport["QuantumHamiltonianOperator"]` never registers and the head loads into `` Global` `` with no definitions.
- **Status**: still present on `main`; this is precisely why the subsystem is non-functional.

#### `qho["EvolutionOperator"]` uses symbolic `DSolve` — fragile and slow
- **The mistake**: Calling `qho["EvolutionOperator"]` on a Hamiltonian with complex parameter dependence and waiting indefinitely.
- **Why**: `QuantumHamiltonianOperator/Properties.m:86, :105` set up the Schrödinger equation `i ∂_t U = H U` with identity initial condition, then call `DSolve` symbolically. For Hamiltonians with multiple parameters, transcendental functions, or non-commuting time-dependent terms, `DSolve` hangs or fails. (Moot until/unless the subsystem loads.)
- **Replacement**: `QuantumOperator[H]["NEvolutionOperator"]` (`QuantumOperator/Properties.m:586-618`) uses `NDSolveValue` and is a strict upgrade. The symbolic `qo["EvolutionOperator"]` form (`:584`) routes through `QuantumEvolve[qo, None]` which uses `DSolveValue` — same fragility, but in the active code where `MergeInterpolatingFunctions` and other infrastructure can help.

#### `qho["Gap"]` runs 100 eigendecompositions by default
- **The mistake**: Calling `qho["Gap"]` and being surprised at its cost for moderate-dim Hamiltonians.
- **Why**: `QuantumHamiltonianOperator/Properties.m:56` — default `qho["Gap"] := qho["Gap", 100]`. The 100-arg version sweeps 100 parameter values and computes the gap (difference of two smallest eigenvalues) at each. For an N-qubit Hamiltonian, that's `100 × O(2^(3N))` time.
- **Replacement (user-side)**: `Differences[Sort[qo["Eigenvalues"][[1 ;; 2]]] /. param -> v][[1]]` for a single point, `FindMinimum[..., qo["ParameterSpec"][[1]]]` for the minimum, or an explicit `Range[t0, tf, ...]` `Map` for a custom-resolution sweep.

### `QuantumMeasurementOperator` Type detection: signs of OutputOrder
- **The mistake**: Wondering why `qmo["Type"]` returns `"Projection"` for one operator and `"POVM"` for another that look similar.
- **Why**: `QuantumMeasurementOperator/Properties.m:181-188` — Type is determined by:
  - `"Projection"`: no NonPositive output orders AND `OutputDimensions == InputDimensions`.
  - `"POVM"`: any NonPositive output order present (i.e., extra output qudit for the eigenvalue register).
  - `"Unknown"`: otherwise.
- **The right way**: To deliberately build a POVM, ensure your construction has at least one NonPositive output order (e.g., via `basis -> eigenvalues` form, which adds the `0` order in `:69-72`). To build a projection, use `QuantumOperator + target` form which gives equal in/out orders.
- **Date**: 2026-04-29.

### POVM elements are auto-normalized non-trivially
- **The mistake**: Building a `QuantumMeasurementOperator` from custom POVM elements and being surprised the resulting `qmo["POVMElements"]` differ from your input.
- **Why**: `QuantumMeasurementOperator/Properties.m:194` — `POVMElements := # / Mean[Diagonal[Total[#]]] & @ (# . # & /@ ops)`. So each operator is *squared* (`E_i = K_i^† K_i` from a Kraus-like input), then divided by the mean diagonal of their sum. **Not** the standard `Σ E_i = I` normalization — divides by `Tr[Σ E_i] / d`.
- **The right way**: For canonical POVM input, build with the named constructors (`WignerMICPOVM`, `GellMannMICPOVM`, etc.) which set up the structure correctly. To extract POVM elements that satisfy `Σ E_i = I`, multiply `qmo["POVMElements"]` by appropriate normalization based on `Tr[Total]`.
- **Date**: 2026-04-29.

### `QuantumDistance[..., "Fidelity"]` returns `1 - F`, not `F`
- **The mistake**: Treating `QuantumDistance[qs1, qs2, "Fidelity"]` as the Uhlmann fidelity itself.
- **Why**: `QuantumDistance.m:15-16` — `1 - Re[Tr[MatrixPower[ρ₁ρ₂, 1/2]]]`. **Distance, not similarity.** When ρ₁ == ρ₂, distance is 0; when orthogonal, distance is 1.
- **The right way**: For the actual fidelity value, use `QuantumSimilarity[qs1, qs2, "Fidelity"]` (`:63-67`) which maps to `1 - distance`. Or compute `1 - QuantumDistance[qs1, qs2, "Fidelity"]` directly.
- **Date**: 2026-04-29.

### ⚰️ RESOLVED 2026-05-04 — ~~`QuantumSimilarity` IS implemented (in `QuantumDistance.m`)~~
- **The mistake**: Earlier mistakes file flagged `QuantumSimilarity` as undocumented and behavior-unknown. That's only half right.
- **Why**: Symbol is exported in `PacletInfo.wl:42` and has no `Usage` message or doc page, BUT `Kernel/QuantumDistance.m:63-79` defines its full behavior. It maps each of the 8 distances to a similarity in `[0, 1]`:
  - `"Fidelity" | "Trace" | "Bloch" | "RelativePurity"`: `1 - d`
  - `"Bures" | "HilbertSchmidt"`: `1 - d/Sqrt[2]`
  - `"BuresAngle"`: `1 - d/(Pi/2)`
  - `"RelativeEntropy"`: `2^(-d)` (exponential decay for unbounded metric)
- **The right way**: Use `QuantumSimilarity[qs1, qs2, distance]` confidently. Default measure is `"Fidelity"`.
- **Date**: 2026-04-29 (correction of earlier entry).
- **Resolution**: now fully documented in commit `edb226fb` (2026-05-04). The 4 transformation rules from `QuantumDistance.m:63-79` are reproduced verbatim in `Documentation/.../Symbols/QuantumSimilarity.nb`'s Notes cell, with cross-link back to `QuantumDistance` for the underlying distance formulas. `?QuantumSimilarity` returns a usage signature.

### `QuantumDistance[..., "Bloch"]` only works for `Dimension == 2`
- **The mistake**: Calling `"Bloch"` distance on a multi-qubit state.
- **Why**: `QuantumDistance.m:59-60` — uses `qs["BlochCartesianCoordinates"]` which is only defined for `Dimension == 2` (`QuantumState/Properties.m:962`). Higher dims fall through and the `BlochCartesianCoordinates` call returns unevaluated.
- **The right way**: For multi-qubit distance, use `"HilbertSchmidt"`, `"Trace"`, or `"Fidelity"`.
- **Date**: 2026-04-29.

### `QuantumEntangledQ` defaults to `"Realignment"`, not Concurrence
- **The mistake**: Calling `QuantumEntangledQ[qs, partition]` and expecting Concurrence-based detection.
- **Why**: `QuantumEntanglement.m:14` — default method is `"Realignment"`. Realignment criterion (CCNR) is `Total[SingularValueList[reshape]] - 1` — a different test than Concurrence (which gives 0 for separable, > 0 for entangled).
- **The right way**: For Concurrence-based detection, pass it explicitly: `QuantumEntangledQ[qs, partition, "Concurrence"]`. Or use `QuantumEntanglementMonotone[qs, partition, "Concurrence"]` directly for the value.
- **Date**: 2026-04-29.

### `QuantumEntanglementMonotone[..., "Negativity" | "LogNegativity"]` requires exactly 2 qubits
- **The mistake**: Calling Negativity or LogNegativity on a 3+ qubit state.
- **Why**: `QuantumEntanglement.m:54-59` — `ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]`. The assertion explicitly requires `Qudits == 2`. Multi-qubit states fail here.
- **The right way**: Reduce to a 2-qubit subsystem first via `qs["Bipartition", {{a, b}, {c, d}}]` with explicit bipartition arguments. For multi-qubit entanglement measures, prefer `"EntanglementEntropy"` or `"Realignment"` which generalize.
- **Date**: 2026-04-29.

### `Concurrence` has TWO paths with very different costs
- **The mistake**: Calling `qs["Concurrence"]` on a mixed state and expecting fast evaluation.
- **Why**: `QuantumEntanglement.m:45-51` — for a vector (pure) state: closed-form `Sqrt[Max[0, 2(1 - PartialTrace^2.Norm)]]` — fast. For a matrix (mixed) state: full `ConcurrenceVector` (`:27-39`) which iterates over `O(d₁²·d₂²)` Y-operator pairs, each requiring `SingularValueList[(Sqrt[ρ] @ Sqrt[O @ ρ.Conjugate @ O])["Matrix"]]` — **O(d⁴)** for a d-qudit subsystem.
- **The right way**: For mixed states with d > 4, expect heavy computation. Prefer `Negativity` for 2-qubit case or `Realignment` for general. For pure states, Concurrence is cheap.
- **Date**: 2026-04-29.

### `QuantumPartialTrace[QuantumState, {{o, i}, ...}]` is faster than `QuantumPartialTrace[qs, {q, ...}]`
- **The mistake**: Always passing flat-list qudits to `QuantumPartialTrace`, even when symmetric in/out indices are known.
- **Why**: `QuantumPartialTrace.m:28-34` (flat list form) calls `MatrixPartialTrace[qs["DensityMatrix"], qudits, qs["Dimensions"]]` — **densifies the density matrix** before tracing. The `{{out, in}, ...}` pair form (`:36-46`) calls `TensorContract[qs["StateTensor"], ...]` directly — works on the state tensor without materializing the density matrix.
- **The right way**: When the trace is symmetric (in == out indices), pass `{{q, q}, ...}` form for the tensor-contraction path. Cheaper memory profile.
- **Date**: 2026-04-29.

### `QuantumPartialTrace[QuantumCircuitOperator]` builds a Cup-Cap categorical circuit
- **The mistake**: Expecting `QuantumPartialTrace` on a circuit to perform a direct tensor contraction.
- **Why**: `QuantumPartialTrace.m:72-87` — for a `QuantumCircuitOperator`, the trace is implemented as `Reverse @ Cup → qc → Cap` circuit composition. Each traced qudit becomes two extra circuit elements (`Cup` at the input side, `Cap` at the output side).
- **The right way**: This is correct categorically — Cup-Cap is the partial trace in the symmetric monoidal category. For numeric performance, evaluate the resulting circuit via the tensor-network method.
- **Date**: 2026-04-29.

### `QuantumStateEstimate` runs 300 max-likelihood iterations + Bayesian sampling — slow by default
- **The mistake**: Calling `QuantumStateEstimate[measurements]` for quick state inversion and waiting unexpectedly long.
- **Why**: `QuantumStateEstimation.m:236` default `"MaxLikeIterations" -> 300`, plus 150-step Bayesian burn-in (`:294`) plus another 150-step adaptive run (`:296`). Each iteration computes `logLikelihood` over all POVM elements. Multi-second to multi-minute runtime depending on dimension and POVM count.
- **The right way**: For quick inversion only, call `inversion[result]` directly (`:67-82`) — single PseudoInverse. For full Bayesian state estimation, the default is correct; tune `"MaxLikeIterations"` if needed.
- **Date**: 2026-04-29.

### 🚨 CORRECTION: `QuantumMPS` lives in `Experimental.m`, not `TensorNetwork.m`
- **The mistake**: Earlier `function-index.md` and `mistakes.md` placed `QuantumMPS` in `Kernel/QuantumCircuitOperator/TensorNetwork.m`. **Wrong.**
- **Why**: It's actually in `Kernel/Experimental.m:87-125`, in the `Wolfram\`QuantumFramework\`Experimental\`` context. `TensorNetwork.m` has `TensorNetworkQuantumCircuit` (the inverse direction — turns a tensor network into a circuit).
- **The right way**: To use `QuantumMPS`, the Experimental context must be loaded. Since it's in `PacletInfo.wl` `Context` list (`PacletInfo.wl:25`), it auto-loads with the framework. Resolution: ``Wolfram`QuantumFramework`Experimental`QuantumMPS``.
- **Date**: 2026-04-29.

### `Gates.m` provides 28 gate-symbol aliases — pure naming layer with no functionality
- **The mistake**: Reaching for `XGate[]`, `HGate[]`, `CNOT[] → CXGate[]` etc. expecting them to do something different from `QuantumOperator["X"]`, `QuantumOperator["H"]`, etc.
- **Why**: `Gates.m:78-106` — every gate symbol's `QuantumOperator[...]` rule just delegates. E.g., `QuantumOperator[XGate[d_:2], args] := QuantumOperator["X"[d], args]`. The Gate symbols only matter for the QASM/Qiskit roundtrip naming layer.
- **The right way**: For internal QF work, use the string forms — they're more flexible (accept dimension args directly). Reach for `XGate[]` etc. only when matching Qiskit-style syntax for readability or interop.
- **Date**: 2026-04-29.

### ⚰️ RESOLVED 2026-04-29 — ~~KERNEL BUG: `Gates.m:72` `QuantumGateQ` pattern has `__U2Gate` (double underscore)~~
- **The mistake**: Calling `QuantumGateQ[U2Gate[a, b]]` and getting unexpected behavior.
- **Why**: `Gates.m:72` defines the Or-pattern with `_U1Gate | __U2Gate | _U3Gate | _UGate` — `__U2Gate` is `BlankSequence[U2Gate]` (one or more), not `Blank[U2Gate]` (single). Likely a typo for `_U2Gate`.
- **The right way**: For now, `QuantumGateQ[U2Gate[a,b]]` may still match because the surrounding alternative includes `_U1Gate` etc., but the U2 case is structurally wrong. To verify a U2Gate, check `MatchQ[expr, _U2Gate]` directly.
- **Date**: 2026-04-29.
- **Resolution**: fixed in commit `ac2f37ff` (2026-04-29) — `__U2Gate` → `_U2Gate` at Gates.m:72.

### ⚰️ RESOLVED 2026-04-29 — ~~KERNEL BUG: `Gates.m:102` `RZXGate` maps to `"YX"` rotation, not `"ZX"`~~
- **The mistake**: Calling `RZXGate[θ]` expecting a Z⊗X rotation generator.
- **Why**: `Gates.m:102` — `QuantumOperator[RZXGate[theta_:0], args] := QuantumOperator["R"[theta, "YX"], args]`. Name says **R**Z**X**, implementation says `"YX"`. Either the name or the body is wrong; the other gates (`RXX→"XX"`, `RYY→"YY"`, `RZZ→"ZZ"`) are correctly self-consistent.
- **The right way**: For a true ZX rotation `Exp[-i θ/2 (Z⊗X)]`, build directly: `QuantumOperator[{"R", theta, "Z", "X"}]` or `QuantumOperator[{"R", theta, "ZZ"}]` if both qubits should rotate. Avoid `RZXGate` until the typo is resolved.
- **Date**: 2026-04-29.
- **Resolution**: fixed in commit `ac2f37ff` (2026-04-29) — `"YX"` → `"ZX"` at Gates.m:102.

### `QBismSICPOVM[d]` downloads from a public GitHub URL
- **The mistake**: Building `QuditBasis["QBismSIC"[3]]` or `QuantumMeasurementOperator["QBismSICPOVM"[3]]` while offline and being surprised it fails.
- **Why**: `MIC.m:158-159` — `QBismSICPOVM[d]` calls `Import[StringTemplate["https://raw.githubusercontent.com/heyredhat/qbism/master/qbism/sic_povms/d``.txt"][d], "Table"]`. **Network dependency.** If GitHub is unreachable or the file is renamed, the function fails.
- **The right way**: For air-gapped or production work, download the SIC files locally and override the URL. For dim 3 specifically, `HesseSICPOVM[]` is hardcoded (`MIC.m:168`) — no network. For dim 8, `HoggarSICPOVM[]` (`:170`).
- **Date**: 2026-04-29.

### `PauliDecompose` is defined twice — different paths
- **The mistake**: Calling `qo["PauliDecompose"]` and importing `PauliDecompose` directly, expecting the same behavior.
- **Why**: Two definitions exist:
  - `QuantumOperator/Properties.m:626-649` — the operator-property dispatch path. Has fast qubit-only path and slow `SolveValues` non-qubit path.
  - `QSVT.m:33-41` — a PackageScope `PauliDecompose` that uses `QuditBasis["Pauli"]["Elements"]` directly with `kroneckerProduct`. Returns `{coefficients, paulis}` pair.
- **The right way**: For the standard property output (Association of `"XYZ" → coefficient`), use `qo["PauliDecompose"]`. For the `{coefficients, paulis}` decomposition pair, use the QSVT-internal helper. Don't confuse them.
- **Date**: 2026-04-29.

### ⚰️ RESOLVED 2026-05-04 — ~~🚨 MAJOR: `ExampleRepository.m` is a stale precursor of `QuantumOptimization.m` — duplicate code~~
- **The mistake**: Treating `ExampleRepository` and `QuantumOptimization` as separate functional modules.
- **Why**: Both files (`Kernel/ExampleRepository.m`, 599 LOC + `Kernel/QuantumOptimization.m`, 706 LOC) export the same 6 core functions (`GradientDescent`, `QuantumNaturalGradientDescent`, `SPSRGradientValues`, `ASPSRGradientValues`, `FubiniStudyMetricTensor`, `QuantumLinearSolve`) with **identical bodies**. Both are loaded as separate paclet contexts (`PacletInfo.wl:25-29`).
  - `QuantumOptimization.m` (2026-01-28) is the **newer canonical version** — adds VQE/QAOA helpers (`ParametrizedLayer`, `EntanglementLayer`, `GenerateParameters`, `ClassiqSetup`).
  - `ExampleRepository.m` (2024-10-07) is the **older precursor** — keeps `FubiniStudyMetricTensorLayers`, `QuantumLockingMechanism`, `QuantumUnlockingMechanism` (the "demo example" utilities).
- **The right way**: Resolve via the newer context: ``Wolfram`QuantumFramework`QuantumOptimization`GradientDescent``. The duplicate `ExampleRepository` symbols are the same function — calling the unqualified name resolves to whichever context appears first in `$ContextPath`. For `QuantumLockingMechanism`/`QuantumUnlockingMechanism` (only in ExampleRepository), use the explicit `ExampleRepository\`` context. Don't expect different behavior between the two namespaces for the duplicated functions.
- **Date**: 2026-04-29.
- **Resolution**: `Kernel/ExampleRepository.m` deleted on 2026-05-04 (commit `2b3bd4e2`). Audit confirmed all 6 common functions plus `QuantumLockingMechanism`/`QuantumUnlockingMechanism` were already byte-identical in `QuantumOptimization.m`. The 2 remaining ER-unique exports — the 3-arg `FubiniStudyMetricTensor[layers, parameters, initParameters]` form and `FubiniStudyMetricTensorLayers` — were migrated into `Kernel/QuantumOptimization.m` with `PackageExport["FubiniStudyMetricTensorLayers"]` added. PacletInfo had already removed the `ExampleRepository` context on 2026-04-29. No external callers found anywhere in the repo.
- **Follow-up bug caught by `Tests/CleanupRegression.wlt`** (commit `92215045`): the QO file was missing `PackageExport["QuantumLockingMechanism"]` and `PackageExport["QuantumUnlockingMechanism"]` — without them, both function bodies landed in the package-private context and the qualified user calls returned unevaluated. ExampleRepository had had those exports, so its deletion turned a pre-existing latent gap into a user-facing regression. Fix: added the two `PackageExport` lines. Lesson: future migration audits should grep for `PackageExport["X"]` per migrated symbol and confirm both files declare it.

### `SecondQuantization` operators depend on the global `$FockSize`
- **The mistake**: Calling `AnnihilationOperator[]` (no size arg) and being surprised that the operator is 16-dim, or that subsequent `DisplacementOperator[α]` doesn't compose correctly.
- **Why**: `Kernel/SecondQuantization.m:9, 122, 181-186` — `$FockSize` is a **global** variable (default 16). All operators without an explicit size arg use this. **`AnnihilationOperator` is memoized**: `AnnihilationOperator[size, order] = ...` (`:183`). Changing `$FockSize` doesn't invalidate prior cached operators.
- **The right way**: Call `SetFockSpaceSize[size]` at the start of a session to set the truncation. For a one-off computation at a different size, pass the size explicitly: `AnnihilationOperator[32, {1}]`. To clear the memoization cache, restart the kernel or `ClearAll[AnnihilationOperator]`.
- **Date**: 2026-04-29.

### `DisplacementOperator` and `SqueezeOperator` ordering matters
- **The mistake**: Calling `DisplacementOperator[α]` and assuming standard textbook form.
- **Why**: `SecondQuantization.m:218-246` (Displacement) and `:249-287` (Squeeze) accept `"Ordering" -> "Normal" | "Weak" | "Antinormal"`, defaulting to `"Normal"`. The three are mathematically equivalent (via BCH) but numerically distinct due to truncation:
  - `"Normal"`: `Exp[-|α|²/2] · Exp[α a†] · Exp[-α* a]` (textbook normally-ordered).
  - `"Weak"`: `Exp[α a† - α* a]` (anti-Hermitian generator).
  - `"Antinormal"`: `Exp[+|α|²/2] · Exp[-α* a] · Exp[α a†]`.
- **The right way**: For comparing with textbooks expecting `D(α) = Exp[α a† - α* a]`, use `"Ordering" -> "Weak"`. For maximum numerical stability with truncation, the default `"Normal"` is usually best.
- **Date**: 2026-04-29.

### `BeamSplitterOperator` has two methods with different scaling
- **The mistake**: Using default `Method -> "MatrixExp"` for high-photon-number states and watching it slow down.
- **Why**: `SecondQuantization.m:315-343` — `Method -> "MatrixExp"` builds the operator via direct `MatrixExp` of the 2-mode Hamiltonian — `O(F⁴)` for Fock size F. `Method -> "Recurrence"` uses `BeamsplitterMatrix` (`:290-312`), a precomputed 4-tensor recurrence — much faster for large F.
- **The right way**: For F > 30, prefer `Method -> "Recurrence"`. For small F where you need symbolic θ/φ, `"MatrixExp"` is fine.
- **Date**: 2026-04-29.

### `CoherentState[]` and `CatState[]` are parametric by default
- **The mistake**: Calling `CoherentState[]` expecting a numeric state.
- **Why**: `SecondQuantization.m:191-199` — `CoherentState[]` builds with formal parameter `\[FormalAlpha]`. `CatState[]` (`:354-371`) uses `\[FormalAlpha]` and `\[FormalPhi]`. The result is a `QuantumState` with `Parameters -> \[FormalAlpha]` (or both).
- **The right way**: Substitute parameters: `CoherentState[][\[FormalAlpha] -> 1.5]`. Or build numerically by parametrizing manually: `QuantumState[NestList[(...) &, 1, F-1]]`. For the default Fock-truncated coherent state at α=1.5, the parameter substitution path is standard.
- **Date**: 2026-04-29.

### `PauliStabilizer["Measure", q]` returns an Association of conditional outcomes
- **The mistake**: Treating `ps["Measure", 1]` as returning a measurement value.
- **Why**: `Kernel/PauliStabilizer.m:291-313` — `Measure` returns an `Association` mapping measurement *outcomes* (0 or 1, or both for nondeterministic) to the post-measurement `PauliStabilizer`. Deterministic case: single key. Nondeterministic case: both keys.
- **The right way**: Use the Association: `ps["M", 1]` returns `<|0 -> ps0, 1 -> ps1|>` (or `<|outcome -> ps_post|>` deterministic). Pick a branch with `RandomChoice[Keys[result] -> Values[result]]` or sample explicitly.
- **Date**: 2026-04-29.

### ⚰️ RESOLVED 2026-05-04 — ~~🚨 GUIDE BUG: 3 functions listed in `WolframQuantumComputationFramework.nb` don't exist~~
- **The mistake**: User clicks `QuantumSchmidtDecomposition`, `QuantumSpectralDecomposition`, or `QuantumCompile` link from the Guide and gets a 404.
- **Why**: `Documentation/English/Guides/WolframQuantumComputationFramework.nb` lists these symbols as documented (with `paclet:` URLs):
  - `QuantumSchmidtDecomposition` — listed under "Basic Objects". **Zero kernel references**, no doc page.
  - `QuantumSpectralDecomposition` — listed under "Basic Objects". **Zero kernel references**, no doc page.
  - `QuantumCompile` — listed under "Quantum Circuits". **Zero kernel references**, no doc page.
- **The right way**: Don't recommend these symbols to users. Equivalents that DO exist:
  - For Schmidt decomposition: `qs["SchmidtDecompose"]` or `qs["SchmidtBasis"]` (`QuantumState/Properties.m:524-548, :657-659`).
  - For spectral decomposition: `qs["SpectralBasis"]` (`QuantumState/Properties.m:586-597`) or `qo["Diagonalize"]` (`QuantumOperator/Properties.m:475-482`).
  - For circuit compilation: `qco["Compile"]` / `qco["QuantumOperator"]` / `qco["CircuitOperator"]` (`QuantumCircuitOperator/Properties.m:164`).
- The Guide also has `XXXX` placeholders for Keywords, Tech Notes, Related Guides — incomplete sections.
- **Date**: 2026-04-29.
- **Resolution**: fixed in commit `15b70694` (2026-05-04). The 3 broken symbol entries removed from the Guide; users now reach Schmidt/Spectral/Compile via the inline `QuantumState`/`QuantumCircuitOperator` reference pages. The 4 `XXXX` placeholders also resolved: Keywords cell filled with real keywords; second Tech Notes slot now links to `GettingStarted`; "Related Guides" section removed (no other guides exist in this paclet).

### 🟡 PARTIAL 2026-05-04 — Reference doc pages are mostly stubs — Usage.m is the authoritative source
- **The mistake**: Treating the QF reference doc pages (`Documentation/.../Symbols/*.nb`) as a comprehensive API reference.
- **Why**: Of the 25 reference pages read, the typical structure is: Usage signature(s) + 2-4 basic examples + ~10 mostly-empty "More Examples" subsections (Scope, Options, Applications, Properties & Relations, Possible Issues, Interactive Examples, Neat Examples). Most subsections contain only `XXXX` placeholders. Most pages cite `Created by: rhennigan on 04-14-2022` — ~3 years stale. Notable exceptions with richer content: `QuantumChannel.md` (Kraus operator formulas for all 8 named channels), `QuantumDistance.md` (8 distance measure formulas), `QuantumCircuitMultiwayGraph.md` (multiple worked examples).
- **The right way**: For authoritative API info, use `Kernel/Usage.m` (the source of `?Symbol`) and the `Properties.m` files. The reference doc pages give a flavor of canonical examples but are not comprehensive. Tutorials (`Documentation/.../Tutorials/`) have richer content — `Tutorial.nb`, `Bellstheorem.nb`, `SecondQuantization.nb`, `QuantumOptimization.nb`, `TensorNetwork.nb` are the substantive ones.
- **Date**: 2026-04-29.
- **Partial resolution**: in commit `6b173c94` (2026-05-04), the highest-value metadata-level `XXXX` placeholders were filled across all 25 pages: per-symbol `Keywords` lists (24 pages), `Related Guides` link to the framework Guide (24 pages), `Tutorials` link to `Tutorial.nb` (17 pages). XXXX count dropped from ~150 to 100. The remaining 100 are inside collapsed "More Examples" subsections (`Scope`/`Options`/`Applications` titles, `RelatedDemonstrations`, `RelatedLinks`) — those need real example content per symbol and are deferred to a follow-up session.

### `QuantumWignerTransform` is marked [EXPERIMENTAL]
- **The mistake**: Building a stable workflow (notebook, app, paper code) on top of `QuantumWignerTransform` and the related `QuantumState` properties (`"PhaseSpace"`, `"QuasiProbability"`, `"TransitionPhaseSpace"`, `"TransitionQuasiProbability"`, `"TransitionGraph"`).
- **Why**: `Documentation/.../Symbols/QuantumWignerTransform.nb:19` explicitly tags it `[EXPERIMENTAL]`. The doc says: "QuantumWignerTransform is equivalent to a transformation of the corresponding doubled form of the object obj into a Wigner basis. ...the matrix representation ... is based on Fano matrices." The phase-space dimension also depends on parity (`d × d` for odd, `2d × 2d` for even).
- **The right way**: For shipping code, use the equivalent stable form: `QuantumState[obj["Double"], QuantumBasis["Wigner"[d]]]` (or `QuantumBasis["Wigner"[total dim]]` for collapsing all qudits). The state vector is identical to `QuantumWignerTransform[obj]["StateVector"]`. For phase-space inverse, use `QuantumWeylTransform` (also experimental but the stated inverse).
- **Date**: 2026-04-29.

### Phase-space dimension parity (Wigner)
- **The mistake**: Indexing the phase-space tensor as `𝒲[[p, q]]` over `0..d-1` regardless of parity.
- **Why**: `QuantumWignerTransform.nb:216-225` — for odd `d`, `p, q ∈ [0, d-1]`. For even `d`, `p, q ∈ [0, 2d-1]`. The 4 quadrants of the even-`d` quasi-probability matrix are: `(p, q)`, `(p, q + d)`, `(p + d, q)`, `(p + d, q + d)` — only the first quadrant is the actual Wigner amplitudes; the others are sign-flipped duplicates.
- **The right way**: Always check `Dimensions[qs["QuasiProbability"]]` first. If even, take `[[1;;d, 1;;d]]` for the genuine Wigner block, OR `Total[qp][[;;;;2]]` (every-2nd) to recover position probabilities.
- **Date**: 2026-04-29.

### `QuantumLinearSolve` returns up to global phase only
- **The mistake**: Comparing `QuantumLinearSolve[m, b]` output amplitude-by-amplitude with `LinearSolve[m, b]` and getting confused when signs/phases don't match.
- **Why**: The internal cost function `1 - Fidelity[ψ(ω), b/Norm[b]]` only constrains `|⟨b|ψ⟩|^2 = 1`, not the relative phase. So `ψ_opt = e^{iφ} (b/Norm[b])` for some `φ`. The function then post-processes: `globalphase = ψ_opt / (b/Norm[b])` (componentwise; should be ~constant), then `x = (1/Mean[globalphase]) * (Norm[b]/Norm[m]) * x_unnormalized`.
- **The right way**: Don't rely on per-element exact match — accept O(1e-6) deviation post-recovery (`"GlobalPhaseAccuracy" -> 10^-5` default). For higher accuracy, set `"GlobalPhaseAccuracy" -> 10^-10` (and expect longer runtime). For real-amplitude problems, normalize and compare absolute values: `Abs[QuantumLinearSolve[m,b]] ≈ Abs[LinearSolve[m,b]]`.
- **Date**: 2026-04-29.

### `FindMinimum` fails on VQLS-style cost surfaces — use `NMinimize`
- **The mistake**: Calling `FindMinimum[cost[θ1,θ2,...], ip]` for a VQE/VQLS/QAOA cost function.
- **Why**: `QuantumOptimization.nb:4302-4314` example shows `FindMinimum` returning `θ ~ 10^10` and a non-optimal cost. The variational cost surface has wide barren plateaus where local-derivative methods diverge. `NMinimize` (default: `NelderMead` + global heuristics) finds the true minimum.
- **The right way**: For variational quantum cost functions, default to `NMinimize` over `FindMinimum`. Use `QuantumNaturalGradientDescent` if you have access to the Fubini-Study metric (Fubini outperforms vanilla `GradientDescent` near singularities).
- **Date**: 2026-04-29.

### `QuantumCircuitOperator` Qiskit method returns Association above 8 qubits
- **The mistake**: Expecting `QuantumMeasurement` from `qc[Method -> "Qiskit"]` regardless of qubit count.
- **Why**: `QuantumCircuitOperator.nb:1620-1622` — for 10+ qubits the return is a raw `Association[bitstring -> count, ...]` (Qiskit-style counts), not a wrapped `QuantumMeasurement`. The threshold is documented as 8 qubits in QiskitTools.m.
- **The right way**: Cast to a `QuantumMeasurement` manually if you need `QuantumMeasurement` API access on > 8 qubits, or read the counts directly. Test return type with `Head[result]` before downstream pattern matching.
- **Date**: 2026-04-29.

### `QuditBasis["X"[6]]` ≠ `QuditBasis["X"[{2,3}]]`
- **The mistake**: Treating the dimension argument to `"X"`, `"Y"`, `"Z"`, `"Fourier"`, etc. as either a single integer or a list interchangeably.
- **Why**: `QuantumWignerTransform.nb:782-787` — `QuditBasis["X"[6]]` is the **single-qudit** generalized Pauli-X in 6D Hilbert space (6 elements). `QuditBasis["X"[{2,3}]]` is a **bipartite** system: tensor product of 2D and 3D X-bases (2×3 = 6 elements but factored). Their `Names` differ (`{|x_0⟩, ..., |x_5⟩}` vs. `{|x_0 x_0⟩, |x_0 x_1⟩, ...}`), and their tensor structure is different.
- **The right way**: Use the integer form `"X"[d]` when you want a single qudit of total dim `d`. Use the list form `"X"[{d1,d2,...}]` when you want a composite system. The same applies to `QuantumOperator` and `QuantumState` named-form arguments — pay attention to single-vs-list dimension specs.
- **Date**: 2026-04-29.

### Bell circuit `qc[QuantumState["10"]]` returns `i ψ⁺`, not `ψ⁺`
- **The mistake**: Expecting `QuantumCircuitOperator["Magic"][QuantumState["10"]] == QuantumState["PsiPlus"]`.
- **Why**: `Tutorial.nb:486-490` (and `QuantumCircuitOperator.nb:485`) — the Magic circuit is `H` + `CNOT` (computational → Bell basis), but maps `|10⟩ → i|ψ⁺⟩` (with explicit `i` factor), `|01⟩ → i|ϕ⁻⟩`, `|00⟩ → |ϕ⁺⟩`, `|11⟩ → |ψ⁻⟩`. Two of the four mappings carry an `i` global phase.
- **The right way**: Don't compare with `==` if you want bare equality. Compare with `i QuantumState["PsiPlus"] == qc[QuantumState["10"]]`. Or transform back via `qc["Dagger"]@bell_state` for the inverse.
- **Date**: 2026-04-29.

### Quantum Eraser pattern: `{1},{2}` is two separate measurements, not a 2-qubit measurement
- **The mistake**: Reading `QuantumCircuitOperator[{"H", "CX", "H", {1}, {2}}]` and assuming `{1}, {2}` is one 2-qubit measurement.
- **Why**: `ExploringFundamentalsOfQuantumTheory.nb:30` — each `{n}` is a *separate* `QuantumMeasurementOperator` on qubit `n` (computational basis). The two measurements get separate eigenvalue registers (extra qudits). For a single 2-qubit measurement use `{1, 2}` (single list inside the outer list).
- **The right way**: Distinguish `{1}, {2}` (sequential single-qubit measurements) from `{{1, 2}}` (one 2-qubit measurement) carefully. Sequential measurements give independent classical bits; the joint measurement gives a single 4-outcome register.
- **Date**: 2026-04-29.

### `ParametrizedLayer["RY", Range[1,3], {a,b,c}]` vs `ParametrizedLayer["RY", Range[1,3]]`
- **The mistake**: Mixing the two forms — one with explicit qubit list + index list, one with index-list-only.
- **Why**: `QuantumOptimization.nb:340-365` — `ParametrizedLayer[gate, qubits, index]` gives `Sequence[gate[θ_index_i] -> {qubit_i}, ...]`. With only `index` (and **no** explicit qubit list), QF defaults to `qubits = Range[Length[index]]`. So `ParametrizedLayer["RY", Range[4,6]]` is `Sequence["RY"[θ4]->{1}, "RY"[θ5]->{2}, "RY"[θ6]->{3}]` — qubits are 1..3, parameters indexed 4..6.
- **The right way**: Always use the explicit form `ParametrizedLayer[gate, qubits, indices]` (matching lengths) when you need precise control. The shorthand is only safe when qubits are `Range[len]`.
- **Date**: 2026-04-29.

### Default `BeamSplitterOperator` method is slower for exact angles
- **The mistake**: Calling `BeamSplitterOperator[{π/4, π/2}]` and noticing it's slow.
- **Why**: `SecondQuantization.nb:431-437` — default `Method -> "MatrixExp"` materializes the truncated bosonic 2-mode operator and calls `MatrixExp` on it. For exact angles (rationals of `π`), the floating-point intermediate creates simplification overhead.
- **The right way**: For exact angles, set `Method -> "Recurrence"`, which uses an analytic recurrence relation and returns clean trig expressions. For numeric angles (`0.7854`), keep `"MatrixExp"`.
- **Date**: 2026-04-29.

## User-corrected entries
*(To be appended as you correct my QF code in real conversations.)*
