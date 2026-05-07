# Stabilizer subsystem roadmap

> Status tracker for items that are **partial**, **deferred**, or have **latent bugs** in the Stabilizer subsystem (Phases 1–6). Each entry has a concrete next step (file path, signature, algorithm sketch, test to add). Updated: 2026-05-06 (Phase 6 API consolidation — see "A.0 — API consolidation" below). Branch: `stabilizer-phases-1-4`.

> **Audit-doc context.** This document complements [`synthesis-implementation.md`](synthesis-implementation.md) (the *what works today* tour) and [`API.md`](API.md) (per-function reference) by recording *what doesn't yet work, and exactly how to finish it*. The source synthesis is at [`OngoingProjects/Stabilizer/package-design-synthesis.md`](../package-design-synthesis.md).

## Overall status

- **Tests:** 249 PauliStabilizer + 32 QuantumDistance + 20 Roundtrips + 34 HybridInterop + 41 CliffordChannel + 41 Correctness + 44 CrossPackage_Stim + 15 CrossPackage_QuantumClifford = **476 / 476 passing**.
- **Phases done:** 1 (refactor + tests), 2 (hygiene), 3 (symbolic phases), 4 (frame + inner products + Pauli measurement), 5a (graph state + LC), 5b (companion MD), 5c (QO/QS round-trip + post-mortem + cross-module audit), 6 + 6.5 (API consolidation 10 → 4 public symbols), 7.1 (hybrid interop UpValues, Pauli fast path), 7.2 (named Pauli channels routed via tableau), **7.3 (extended Pauli-label detection: Times[-1, ...]/Superscript forms)**, 8.1 (CliffordChannel scaffolding per Yashin25 §2.3), **8.2 (CliffordChannel composition via Boolean null space + cc[ps] state evolution + QC->CC round-trip)**.
- **Public surface:** 5 top-level symbols (`PauliStabilizer`, `StabilizerFrame`, `CliffordChannel`, `GraphState`, `LocalComplement`) + 6 method-grade operations on `PauliStabilizer` (`SymbolicMeasure`, `SubstituteOutcomes`, `SampleOutcomes`, `InnerProduct`, `Expectation`, plus the `["Random", n]` named pattern). `StabilizerFrame` carries `["InnerProduct", other]` as well. `CliffordChannel` is the new Phase 8.1 head; once Phase 8.2 (composition) lands, it becomes the underlying tableau type for `PauliStabilizer` / `StabilizerFrame` / measurement operations (Yashin25 unification).
- **Branches:** `stabilizer-phases-1-4` (Phases 1–6.5), `claude/phase7-hybrid-interop` (7.1+7.2 stacked off 6.5), `claude/phase8-clifford-channel` (8.1 stacked off 7.2).
- **Open items:** 19 (12 partial + 7 deferred). A.9 reclassified as **inherent trade-off**. A.11/A.12/A.13 added 2026-05-04 from refpage-validation findings. B.4 (`CliffordTableau` distinct head) **dropped 2026-05-06**: superseded by the proposed Phase 8 (Yashin25 Choi-tableau unifier — see B.2).
- **Synthesis priority hits**: 5 / 10 ✅; 2 / 10 ⚠️ partial; 3 / 10 ⏸ deferred (per `package-design-synthesis.md` §11).

> **Lesson (2026-05-04).** Phase 5c surfaced a structural test-writing miss: TIER 1.4 was labeled "Round-trips" but contained no `A[C[x]] === x` exact-equality test for any constructor-accessor pair. The Y round-trip bug (`PauliStabilizer[QuantumOperator["Y"]]["QuantumOperator"] = i·Y`) was therefore invisible to a 185-test suite. Full root-cause + design rationale: [`post-mortem-phase-5c.md`](post-mortem-phase-5c.md). Process rule for future test-writing: [`feedback_user_facing_roundtrip_first.md`](~/.claude/projects/-Users-mohammadb-Documents-GitHub-QuantumFramework/memory/feedback_user_facing_roundtrip_first.md).

> **Lesson (2026-05-06, Phase 6).** Phases 3–5 added 10 public symbols organically, one per literature primitive (FangYing23, GarMar15, AndBri05, KoeSmo14). A design review with N. Murzin flagged the surface as excessive. Phase 6 demoted four (`StabilizerMeasure`, `SubstituteOutcomes`, `SampleOutcomes`, `RandomClifford`) to method-grade operations, no semantic change. Process rule for future phases: end every feature-growth phase with a "consolidate the API surface" pass before merging.

---

## Phase 5c — DONE (2026-05-04)

`QuantumOperator` ↔ `PauliStabilizer` and `QuantumState` ↔ `PauliStabilizer` now round-trip exactly on the *first hop* (constructor → accessor). Six commits on `stabilizer-phases-1-4`:

| Commit | Step | What |
|---|---|---|
| `93edc400` | (1) Lock the spec | New TIER 1.4a/1.4b/1.4c blocks in [`Tests/PauliStabilizer.wlt`](../../../Tests/PauliStabilizer.wlt) (42 new tests). 21 red on commit. |
| `6f8637ee` | (2) State tomography | Rewrite `PauliStabilizer[qs_QuantumState]` ([`Stabilizer/Constructors.m`](../../../QuantumFramework/Kernel/Stabilizer/Constructors.m)). The old `ResourceFunction["RowSpace"]` path canonicalized away generator signs — `PauliStabilizer[QuantumState[{0,1}]]` returned `\|0⟩`'s tableau, Bell phi+ returned phi-'s tableau. New path: 4ⁿ Pauli expectations → keep `\|⟨P⟩\| ≈ 1` as stabilizer group → greedy F₂ rank growth picks `n` independent generators with their signs → symplectic Gram-Schmidt extends to destabilizers. Adds `"GlobalPhase"` association key (default 1) so `["State"]` replays the overall phase. 21 → 6 red. |
| `f7ebfd56` | (3) Operator phase capture | `PauliStabilizer[qo_QuantumOperator]` ([`Stabilizer/Constructors.m`](../../../QuantumFramework/Kernel/Stabilizer/Constructors.m)) computes `phase = qo["Matrix"] / recovered["Matrix"]` at the first nonzero entry and stores under `"GlobalPhase"`. `["QuantumOperator"]` ([`Stabilizer/Conversions.m`](../../../QuantumFramework/Kernel/Stabilizer/Conversions.m)) multiplies the canonical AG-decomposed operator by that phase. Fixes Y, YY, XY, YX, YZ, ZY round-trip. 6 → 0 red. |
| `57d8b53f` | (4) Doc updates | ROADMAP / API.md / synthesis-implementation.md reflect 5c. |
| `f9e2d711` | (5) Post-mortem + A.9 contract + TIER 1.4d | New [`post-mortem-phase-5c.md`](post-mortem-phase-5c.md) writes down the structural reasons the Y bug was missable, including the explicit phase-aware vs phase-oblivious design rationale. Comment block added to [`Stabilizer/GateUpdates.m`](../../../QuantumFramework/Kernel/Stabilizer/GateUpdates.m) documenting that `"GlobalPhase"` is intentionally dropped on gate updates (A.9 is an inherent trade-off, not a closeable bug — see §A.9 below). New TIER 1.4d block with up-to-phase tests + 2 "escape hatch" exact-equality tests showing the user-facing workaround. |
| `cbf51576` | (6) Cross-module audit | New [`Tests/Roundtrips.wlt`](../../../Tests/Roundtrips.wlt) probes the two QF pairs the audit flagged as untested for exact-equality round-trip: `QuantumOperator ↔ QuantumChannel` and `QuantumOperator ↔ QuantumCircuitOperator`. 20/20 PASS — no new bugs surfaced. The original Y bug was unique to the AG-tableau projection. |

Removed: `PauliStabilizerTableau` (the old broken tomography helper, 14 LOC).

The fix surfaces one new partial item — **A.10** (cosmetic, generator-order canonicalization) — and reclassifies **A.9** as an inherent trade-off with a documented contract (see below).

### Phase 5c — follow-up commits (2026-05-04)

After `cb043441` (which produced the 5c content above), four follow-up commits polished the docs and tests without changing the underlying design:

| Commit | What |
|---|---|
| `070eb336` | `verify-API.wls` block reorder to match `API.md` doc order. No content change. |
| `3e0539f7` | `API.md` documents the `"GlobalPhase"` association key + the round-trip contract (`A[C[x]] === x` on first hop; up-to-phase under gate updates per §A.9). |
| `90639b6c` | New TIER 1.4e block in [`Tests/PauliStabilizer.wlt`](../../../Tests/PauliStabilizer.wlt): coverage matrix asserting every accessor in `_PauliStabilizer["Properties"]` is exercised by at least one test (+15 tests, 235 → 250). Closes the structural gap that motivated the post-mortem. |
| `abc51f3a` | TIER 1.1 (`multiplication-via-symplectic`) rewritten to be derivation-driven, then embedded verbatim in `API.md` so the doc and the test cannot drift. No test-count change. |
| `48072598` + (next) | **Reference page drafts** for all 10 Stabilizer subsystem public exports (`PauliStabilizer`, `RandomClifford`, `StabilizerMeasure`, `SubstituteOutcomes`, `SampleOutcomes`, `StabilizerInnerProduct`, `StabilizerExpectation`, `StabilizerFrame`, `GraphState`, `LocalComplement`). Pipeline: 4 parallel `general-purpose` subagents wrote markdown to `/tmp/refpage-*.md` against [`refpage-template.md`](~/.claude/skills/documentation-writing/references/refpage-template.md); `wolframscript`-validated 110 code blocks (107 pass, 3 fail); `nb-writer-v2` (kernel-driven JSON→NB serializer) produced structurally valid `.nb` files. The 3 failing blocks surfaced kernel bugs (now ROADMAP §A.11/§A.12/§A.13), not doc bugs. The 10 `.nb` files now live in [`OngoingProjects/Stabilizer/Documentation/Examples for doc pages/`](Examples%20for%20doc%20pages/) as draft material — ready for promotion to `QuantumFramework/Documentation/English/ReferencePages/Symbols/` once the §A.11–§A.13 bugs are fixed and the pages are polished. |

Two process artifacts also produced:

- `~/.claude/projects/-Users-mohammadb-Documents-GitHub-QuantumFramework/memory/feedback_user_facing_roundtrip_first.md` — global memory rule: "for any QF constructor-accessor pair `C`/`A`, the FIRST test in the tier must be `A[C[x]] === x`". Will load in future Claude sessions on this repo.
- This document's **Lesson** callout above and the post-mortem are cross-linked from API.md and synthesis-implementation.md.

---

## Phase 6 — DONE (2026-05-06): API consolidation

Design review with N. Murzin flagged the post-Phase-5c public surface (10 symbols) as excessive. Phase 6 hard-removes 4 of them and replaces with method-grade operations on `PauliStabilizer`. **No semantic change**; all functionality is preserved as `ps[...]` methods. Net test count: 250 → 249 (one redundant `RandomClifford::usage` test dropped).

| Old top-level | New form | Rationale |
|---|---|---|
| `RandomClifford[n]` | `PauliStabilizer["Random", n]` | Already existed as a named-pattern dispatch ([API.md:118-128](API.md)); the standalone alias was redundant. |
| `StabilizerMeasure[ps, q]` | `ps["SymbolicMeasure", q]` | Is a method on a `PauliStabilizer`; symmetric with `ps["M", q]`. |
| `SubstituteOutcomes[ps, rules]` | `ps["SubstituteOutcomes", rules]` | Same — single-receiver method. |
| `SampleOutcomes[ps, n]` | `ps["SampleOutcomes", n]` | Same. |

**Files touched** (1 commit set):
- [PacletInfo.wl](../../../QuantumFramework/PacletInfo.wl) — Symbols list 10 → 6 (lines 51–54 dropped).
- [Kernel/Usage.m](../../../QuantumFramework/Kernel/Usage.m) — drop 4 standalone usage messages; extend `PauliStabilizer::usage` to mention the new method forms.
- [Kernel/Stabilizer/RandomClifford.m](../../../QuantumFramework/Kernel/Stabilizer/RandomClifford.m) — `PackageExport[RandomClifford]` → `PackageScope[RandomClifford]`. The bare `RandomClifford[n]` definition is unchanged; only the export status changes. The named-pattern dispatch `PauliStabilizer["Random", n] := RandomClifford[n]` continues to resolve within the package.
- [Kernel/Stabilizer/SymbolicMeasure.m](../../../QuantumFramework/Kernel/Stabilizer/SymbolicMeasure.m) — three exports → `PackageScope[symbolicMeasure]`, `PackageScope[substituteOutcomes]`, `PackageScope[sampleOutcomes]`. Internal callsites renamed accordingly.
- [Kernel/Stabilizer/Properties.m](../../../QuantumFramework/Kernel/Stabilizer/Properties.m) — added 5 method dispatches: `ps["SymbolicMeasure", q_Integer]`, `ps["SymbolicMeasure", qudits_List]`, `ps["SubstituteOutcomes", rules_]`, `ps["SampleOutcomes"]`, `ps["SampleOutcomes", n_Integer ? Positive]`.
- [Kernel/Stabilizer/PauliStabilizer.m](../../../QuantumFramework/Kernel/Stabilizer/PauliStabilizer.m) — `_PauliStabilizer["Properties"]` registry extended with `"SymbolicMeasure"`, `"SubstituteOutcomes"`, `"SampleOutcomes"`.
- [Tests/PauliStabilizer.wlt](../../../Tests/PauliStabilizer.wlt) — TIER 5 + TIER 6 callsites rewritten; `Phase2-RandomClifford-IsPublic` renamed to `Phase2-PauliStabilizerRandom-Valid`; `Phase2-RandomClifford-Callable` renamed to `Phase2-PauliStabilizerRandom-Callable`; `Phase2-RandomClifford-UsageMessage` removed.
- [Documentation/API.md](API.md) — public-API table 10 → 6 rows + new "Method-grade operations" subsection. Standalone sections for the four demoted symbols folded into one "Methods (Symbolic measurement)" + "Methods (Random Clifford)" pair under `# PauliStabilizer`.

### Phase 6.5 — DONE (2026-05-06): borderline cases demoted

After Phase 6 landed, two more public symbols flagged in the design review (the borderline cases) were demoted to methods. Public surface 6 → 4. No semantic change.

| Old top-level | New form | Files |
|---|---|---|
| `StabilizerInnerProduct[a, b]` | `a["InnerProduct", b]` (and `frame["InnerProduct", b]`) | [Stabilizer/InnerProduct.m](../../../QuantumFramework/Kernel/Stabilizer/InnerProduct.m) → `PackageScope[stabilizerInnerProduct]`; method dispatches in [Stabilizer/Properties.m](../../../QuantumFramework/Kernel/Stabilizer/Properties.m) and [Stabilizer/StabilizerFrame.m](../../../QuantumFramework/Kernel/Stabilizer/StabilizerFrame.m). |
| `StabilizerExpectation[ps, pauli]` | `ps["Expectation", pauli]` | [Stabilizer/InnerProduct.m](../../../QuantumFramework/Kernel/Stabilizer/InnerProduct.m) → `PackageScope[stabilizerExpectation]`. Message tag renamed `StabilizerExpectation::dim` → `PauliStabilizer::expectationdim`. |

[`_PauliStabilizer["Properties"]`](../../../QuantumFramework/Kernel/Stabilizer/PauliStabilizer.m) extended with `"InnerProduct"` and `"Expectation"`; [`_StabilizerFrame["Properties"]`](../../../QuantumFramework/Kernel/Stabilizer/StabilizerFrame.m) extended with `"InnerProduct"`. [`PauliStabilizer::usage`](../../../QuantumFramework/Kernel/Usage.m) and [`StabilizerFrame::usage`](../../../QuantumFramework/Kernel/Usage.m) updated.

**Phase 7 — Hybrid interop (in progress; separate branch / PR):**

Phase 7 is split into two sub-phases. Phase 7.1 (the scaffolding + Pauli-basis fast path) lives on branch `claude/phase7-hybrid-interop` (off `stabilizer-phases-1-4`).

- **Phase 7.1 — DONE (2026-05-06).** UpValues `qmo[ps]`, `qmo[sf]`, `qc[ps]`, `qc[sf]` attached to `PauliStabilizer` / `StabilizerFrame` ([Stabilizer/HybridInterop.m](../../../QuantumFramework/Kernel/Stabilizer/HybridInterop.m)). When the QMO's operator label is a Pauli string, the dispatcher routes to the existing `ps["M", pauli]` AG fast path (`O(n²)`, stays in tableau). When the basis is non-Pauli, the dispatcher emits `PauliStabilizer::nonpaulibasis` and falls back to the legacy `PauliStabilizerApply[QuantumCircuitOperator[qmo], ps]` generic path (which itself was previously the unconditional behavior at `Conversions.m:140`, removed in this commit). Tests in [Tests/HybridInterop.wlt](../../../Tests/HybridInterop.wlt) (8 tests, all passing). 309 / 309 tests across the whole suite. Design rationale (UpValues over `Picture` flag over `QuantumBasis`-wrap) documented in [post-mortem-phase-5c.md](post-mortem-phase-5c.md) Phase 6 footnote.
- **Phase 7.2 — DONE (2026-05-06)** *(channel half).* `qc[ps_PauliStabilizer]` now detects the four named Pauli channels (`BitFlip`, `PhaseFlip`, `BitPhaseFlip`, `Depolarizing`) by Label and returns a probabilistic-mixture list `{{prob, ps_after_pauli}, ...}` where each branch is the original `ps` with the Kraus-corresponding Pauli gate applied via tableau update (`O(n)` per branch). Non-Pauli channels (`AmplitudeDamping`, `PhaseDamping`, `GeneralizedAmplitudeDamping`, `ResetError`) still fall back to dense materialization — their Kraus operators are not all Paulis. New helper [`stabilizerCliffordChannelMixture`](../../../QuantumFramework/Kernel/Stabilizer/HybridInterop.m) does the detection. Tests: 7 new in [Tests/HybridInterop.wlt](../../../Tests/HybridInterop.wlt) (BitFlip num-branches + identity-branch + probability-sum, PhaseFlip / BitPhaseFlip branches, Depolarizing 4-branch + probability-sum, AmplitudeDamping fallback). 15 / 15 in HybridInterop, 316 / 316 across the full suite.
- **Phase 7.3 — DONE (2026-05-06)** *(label-detection extension).* Extended `stabilizerPauliLabelFromQMO` to recognize three more label forms: `Times[-1, Superscript[X|Y|Z|I, CircleTimes[m]]]` → `"-XXXX..."`; `Superscript[X|Y|Z|I, CircleTimes[m]]` → `"XXXX..."`; `Times[-1, str]` for Pauli string `str` → `"-" <> str`. This routes `qmo[ps]` through the AG tableau fast path for QMOs built as `QuantumOperator[-"XX"]` and similar tensor-power Pauli expressions. Tests added in [Tests/HybridInterop.wlt](../../../Tests/HybridInterop.wlt) TIER J/K (5 new).
- **Phase 7.4 — PARTIAL DONE (2026-05-06)** *(matrix-iteration detector half).* Extended `stabilizerPauliLabelFromQMO` with a `stabilizerPauliFromMatrix` fallback that, for `n ≤ $stabilizerPauliMatrixSearchMaxQubits` (default 4) qubits, iterates `4ⁿ × {+1, -1}` candidate Pauli strings against `qmo["Operator"]["MatrixRepresentation"]` and returns the matching string. This catches QMOs built as `QuantumOperator[matrix]` (where the symbolic Label is `None`) and routes them through the AG tableau fast path. Cap is PackageScope-tunable via `Block`. Tests added in [Tests/HybridInterop.wlt](../../../Tests/HybridInterop.wlt) TIER L (15 new). 52 / 52 HybridInterop tests passing. Remaining sub-task: real `StabilizerFrame` decomposition for arbitrary low-rank non-Pauli QMO basis vectors (deferred to Phase 7.5; depends on richer mixed-state stabilizer support).

**Phase 8 — Yashin25 Choi-tableau unifier (architectural, separate PR):** promote [B.2](#b2--clifford-channel-via-choi-tableau-yashin25) from "deferred" to "next milestone". Implements `CliffordChannel[<|"UA", "UB", "c"|>]` per Yashin25 §2.3. After Phase 8, `PauliStabilizer` / `StabilizerFrame` / the `SymbolicMeasure` family become facade methods over `CliffordChannel`. Public API unchanged. ~400 LOC + 15 tests.

---

## A. Partial implementations (work-but-incomplete)

### A.1 — `ps["InnerProduct", other]` closed-form (`O(n³)`)
| | |
|---|---|
| **Current** | Direct vector materialization, `O(2ⁿ)` memory + time. Works for `n ≤ ~8`. |
| **Source** | [`QuantumFramework/Kernel/Stabilizer/InnerProduct.m`](../../../QuantumFramework/Kernel/Stabilizer/InnerProduct.m) lines 26–35 |
| **Why deferred** | Phase 4 v1 prioritized correctness over performance. Closed-form requires Gaussian elimination over `𝔽₂` plus phase-factor accumulation. |
| **Reference** | GarMarCro12 §3 (arxiv:1210.6646); deSilSalYin23 §3 (faster amplitude path) |
| **Next step** | Implement `Method -> "ClosedForm"` option. Algorithm: (a) compute generator-set intersection `G_ψ ∩ G_φ` via Gaussian elimination over `𝔽₂`; (b) check sign-disagreement (return 0 if found); (c) return `2^(-s/2)` where `s = n − dim(G_ψ ∩ G_φ)`. |
| **File** | extend [`Stabilizer/InnerProduct.m`](../../../QuantumFramework/Kernel/Stabilizer/InnerProduct.m) |
| **Tests to add** | n = 10, 12 cases (where direct-vector OOMs); cross-check vs. direct-vector for n ≤ 8 with seeded `RandomClifford` |
| **Effort** | ~150 LOC kernel + ~10 tests |

### A.2 — `ps["Expectation", pauli]` AG-phase i-factor tracking
| | |
|---|---|
| **Current** | When `P ∈ ⟨g_i⟩` (commutes with all stabilizers but is in their span), falls back to direct vector for `<P>` because the `𝔽₂`-decomposition misses the `i`-factor from `Y = iXZ`. |
| **Source** | [`Stabilizer/InnerProduct.m`](../../../QuantumFramework/Kernel/Stabilizer/InnerProduct.m) lines 53–64; symptom: `<Bell\|YY\|Bell> = -1` requires fallback (synthesis-implementation.md §2.5) |
| **Why deferred** | The `agPhase[x1, z1, x2, z2]` function in `Measurement.m` already handles the per-pair phase. Need to thread it through a sequence of binary multiplications recovered by `LinearSolve`. |
| **Reference** | AarGot04 §3 `g`-function; Yashin25 Eq 4 (cocycle) |
| **Next step** | Replace the direct-vector fallback at `InnerProduct.m:62-64` with: iterate the `recoverable` coefficient list `{c_1, …, c_n}`, multiply generators in order using `Mod[BitXor[currentVec, generator_i], 2]` for the `𝔽₂` part and `Mod[totalPhase + agPhase[…], 4]` for the `i`-factor accumulator. Final `<P>` sign = `(-1)^(totalPhase / 2)` × `targetSign` × `Π signs[i]^c_i`. |
| **File** | [`Stabilizer/InnerProduct.m`](../../../QuantumFramework/Kernel/Stabilizer/InnerProduct.m); reuse `agPhase` from [`Stabilizer/Measurement.m`](../../../QuantumFramework/Kernel/Stabilizer/Measurement.m) |
| **Tests to add** | `Phase6-Expectation-YY-Closed-Form-NoFallback` — no direct-vector fallback issued. Run on n = 8, 10 cases. |
| **Effort** | ~50 LOC + ~5 tests |

### A.3 — `["Dagger"]["Dagger"]` infinite recursion (latent bug from Phase 1)
| | |
|---|---|
| **Current** | `ps["Dagger"]["Dagger"]` does NOT terminate; hits `TerminatedEvaluation["IterationLimit"]`. Discovered during Phase 1 integration check (Tests/PauliStabilizer.wlt does not currently exercise `Dagger ∘ Dagger` — only single `["Dagger"]`). |
| **Source** | [`Stabilizer/GateUpdates.m`](../../../QuantumFramework/Kernel/Stabilizer/GateUpdates.m) lines 95–105 (definition of `["Dagger" \| "Inverse"]`) |
| **Why** | The recursive `ps[PauliStabilizer[<|"Matrix" -> mat, "Signs" -> ps["Signs"]\|>]]["Phase"]` call constructs an inner PauliStabilizer without a paired tableau-padding, then composes with itself, triggering infinite recursion. |
| **Reference** | None — kernel-level bug |
| **Next step** | Refactor to compute the inverse phase directly without composing the inverse PauliStabilizer with the original. The right algorithm: invert the matrix mod 2 to get `M⁻¹`, then compute the new phase as `Phase(M⁻¹)` using the symplectic phase formula directly (no recursive composition). |
| **File** | [`Stabilizer/GateUpdates.m`](../../../QuantumFramework/Kernel/Stabilizer/GateUpdates.m) |
| **Tests to add** | `Roundtrip-DaggerInvolution` — `ps["Dagger"]["Dagger"] === ps` for Bell, GHZ-3, 5Q (use the canonical-form equality of the post-Dagger-Dagger tableau against the original). |
| **Effort** | ~30 LOC + 3 tests |

### A.4 — `StabilizerMeasure` deterministic-outcome correlation
| | |
|---|---|
| **Current** | `StabilizerMeasure` returns the post-state but DROPS the deterministic outcome polynomial. `SampleOutcomes` cannot recover Bell ZZ correlation from stabilizer signs alone. Tracked in `Tier 6 KNOWN LIMITATIONS` of `Tests/PauliStabilizer.wlt` (`Phase3-LIMITATION-DeterministicOutcomeNotStamped`). |
| **Source** | [`Stabilizer/SymbolicMeasure.m`](../../../QuantumFramework/Kernel/Stabilizer/SymbolicMeasure.m) lines 50–82 |
| **Why deferred** | Phase 3 implementation predates `StabilizerFrame`. Fix requires explicit per-measurement outcome record, which is exactly what frames provide. |
| **Reference** | FangYing23 §3 (the right way to track outcome polynomials via fresh symbols + outcome-record matrix) |
| **Next step** | Extend `StabilizerMeasure` to return either `<\|"State" -> ps, "Outcome" -> polynomial\|>` or wrap in `StabilizerFrame` with per-measurement components. Specifically: track an "Outcomes" key on a new wrapper head. |
| **File** | extend [`Stabilizer/SymbolicMeasure.m`](../../../QuantumFramework/Kernel/Stabilizer/SymbolicMeasure.m) and possibly extend `StabilizerFrame` |
| **Tests to add** | Flip `Phase3-LIMITATION-DeterministicOutcomeNotStamped` from `{{0, 1}, {0}}` to `{{0, 1}, {0, 1}}` (m₂ mirrors m₁ across all 20 random Bell samples). |
| **Effort** | ~80 LOC + 4 tests |

### A.5 — `["Circuit"]` synthesis is greedy AG only (not optimal)
| | |
|---|---|
| **Current** | `ps["Circuit"]` produces a Clifford circuit whose dagger equals `ps`, via the AG greedy algorithm. Circuit length is **not minimized** — see audit `paulistabilizer-source-audit.md §8.2`. |
| **Source** | [`Stabilizer/Conversions.m`](../../../QuantumFramework/Kernel/Stabilizer/Conversions.m) lines 87–115 |
| **Why deferred** | Reid24 / Winderl23 algorithms are substantial (≥ 200 LOC each). Phase 1 prioritized correctness. |
| **Reference** | Reid24 (arxiv:2404.19408 — pivot-based 2Q-gate minimization); Winderl23 (arxiv:2309.08972 — Steiner-tree-based connectivity-aware) |
| **Next step** | Add `Method -> "Reid"` and `Method -> "Winderl"` options. Both algorithms iteratively pivot a qubit, sanitize destabilizer + stabilizer rows, then remove inter-qubit interactions. Winderl additionally restricts CNOT to a Steiner tree over a connectivity graph. |
| **File** | new [`Stabilizer/Synthesis.m`](../../../QuantumFramework/Kernel/Stabilizer/Synthesis.m) |
| **Tests to add** | Gate-count comparison vs. AG (Reid should be ≤); connectivity-respecting verification for Winderl (no CNOT outside the connectivity graph). Cross-check that resulting circuits, when applied to `\|0⟩⊗ⁿ`, still produce the right stabilizer state. |
| **Effort** | ~250 LOC + 8 tests |

### A.6 — `LocalComplement` does not update VOPs
| | |
|---|---|
| **Current** | `LocalComplement[gs, v]` toggles edges among `v`'s neighbors but copies VOPs unchanged. AndBri05 Theorem 1 specifies the VOP update rule (Eq 8 of the paper), which we don't apply. |
| **Source** | [`Stabilizer/GraphState.m`](../../../QuantumFramework/Kernel/Stabilizer/GraphState.m) lines 130–139 |
| **Why deferred** | Requires the 24-element LocalClifford table (item B.3) to interpret VOP indices and compose them. |
| **Reference** | AndBri05 §2 footnote, Eq 8 |
| **Next step** | (1) Build the 24×24 LocalClifford composition table per item B.3. (2) For each neighbor `n` of `v` after LC, update `VOPs[n] ← compose(VOPs[n], √(-iX))`. For `v` itself, `VOPs[v] ← compose(VOPs[v], √(iZ))`. |
| **File** | extend [`Stabilizer/GraphState.m`](../../../QuantumFramework/Kernel/Stabilizer/GraphState.m) (after item B.3 lands) |
| **Tests to add** | LC followed by `["PauliStabilizer"]` reproduces the same stabilizer state up to a known local Clifford (cross-check against direct construction). |
| **Effort** | ~30 LOC after item B.3 + 4 tests |

### A.7 — `StabilizerEntanglement[ψ, partition]` via Schmidt rank only
| | |
|---|---|
| **Current** | Bipartite entanglement entropy computed via `MatrixRank @ ArrayReshape[ps["State"]["StateVector"], {2^|A|, 2^|B|}]`. Cost `O(2ⁿ)` memory. Works for `n ≤ ~10`. |
| **Source** | Used in `Tier 3.D` of `Tests/PauliStabilizer.wlt` |
| **Why deferred** | Closed-form formula for stabilizer entropy via the rank of `[X_A \| Z_A]` mod 2 is more nuanced than I expected (synthesis §3.5 references it; Fattal et al.). |
| **Reference** | Fattal, Cubitt, Yamamoto, Bravyi, Chuang 2004 (entanglement entropy of stabilizer states); deSilSalYin23 |
| **Next step** | Implement `StabilizerEntropy[ps, partition]` using the formula `S(A) = |A| − rank_F2([gens restricted to A])` where the restriction is the X+Z bits projected onto A's qubit indices. Closed-form, polynomial-time. |
| **File** | new `Stabilizer/Entropy.m` (or extend `Properties.m`) |
| **Tests to add** | Cross-check vs. Schmidt rank for n ≤ 6 (5 cases including Bell, GHZ, cluster); verify n = 12 case where Schmidt rank OOMs but closed-form succeeds. |
| **Effort** | ~80 LOC + 6 tests |

### A.8 — Coverage table partial markers in `synthesis-implementation.md`
Three items in the §4–§6 menus are explicitly listed as `partial` in [`synthesis-implementation.md`](synthesis-implementation.md):
- §4 `StabilizerStateQ[ψ]` — needs a top-level public symbol; currently use `MatchQ[..., _PauliStabilizer]`. **Effort:** ~5 LOC + 2 tests.
- §5 `CliffordTableau[U]` — needs a distinct head from `PauliStabilizer` (which is a state's tableau, not a gate's). **Effort:** ~100 LOC + 5 tests. Tied to item B.4.
- §6 `StabilizerRankDecomposition[ρ]` — currently relies on user manually constructing a `StabilizerFrame`. Need an automatic decomposition (Bravyi 2016). **Effort:** ~200 LOC + 5 tests.

### A.9 — `"GlobalPhase"` does not propagate through gate updates **(inherent trade-off; documented as contract)**
| | |
|---|---|
| **Current** | Phase 5c added a `"GlobalPhase"` association key set by `PauliStabilizer[qs_QuantumState]` and `PauliStabilizer[qo_QuantumOperator]`, so the *first* round-trip is exact. Gate updates (`ps["H", 1]`, `ps["CNOT" -> {1, 2}]`, etc.) intentionally do not copy `"GlobalPhase"` — see the comment block at the top of [`Stabilizer/GateUpdates.m`](../../../QuantumFramework/Kernel/Stabilizer/GateUpdates.m). As a result `PauliStabilizer[qs]["H", 1]["State"]` is correct *up to* a global phase. |
| **Source** | [`Stabilizer/GateUpdates.m`](../../../QuantumFramework/Kernel/Stabilizer/GateUpdates.m) — every gate-update rule constructs a fresh association without `"GlobalPhase"`. |
| **Why this is inherent (not "deferred")** | Investigated 2026-05-04 (commit `f9e2d711`). The phase factor a Clifford gate `U` picks up when applied to a stabilizer state depends not just on `U` but on the input state's tableau. Concrete counterexamples: `Z\|0⟩ = \|0⟩` (no phase) vs `Z\|1⟩ = −\|1⟩` (sign flip) — same `Z` gate, different state, different phase. Same for `S`, `Y`, etc. So no constant per-gate "phase factor" rule exists. Recovering the new `"GlobalPhase"` exactly requires materializing the state vector at every gate update — `O(2ⁿ)` per gate, which defeats the AG `O(n²)` complexity advantage that motivates the tableau representation in the first place. |
| **Reference** | Post-mortem [`OngoingProjects/Stabilizer/Documentation/post-mortem-phase-5c.md`](post-mortem-phase-5c.md) §2 (phase-aware vs phase-oblivious design). |
| **Documented contract** | `PauliStabilizer[qs]["gate", q]["State"]` equals `gate @ qs` *up to global phase*. For exact equality, re-run the constructor on the post-gate state: `PauliStabilizer[gate @ qs]["State"] === gate @ qs`. Tests at [`Tests/PauliStabilizer.wlt`](../../../Tests/PauliStabilizer.wlt) TIER 1.4d (8 tests including 2 "escape hatch" exact-equality tests). |
| **Future option** | If exact phase tracking through gate updates is wanted, a separate kernel option (`Method -> "PhaseAware"` or similar) could opt into the `O(2ⁿ)` materialization on every update. Out of scope for the current AG-tableau design. |
| **Effort to upgrade to opt-in** | ~80 LOC + 15 tests, *if* the design choice is made. |

### A.10 — Generator order differs between default-register and QS-derived paths
| | |
|---|---|
| **Current** | `QuantumCircuitOperator[circ][Method -> "Stabilizer"]["Stabilizers"]` returns generators in one order; `QuantumCircuitOperator[circ][qs_QuantumState, Method -> "Stabilizer"]["Stabilizers"]` returns the same group in a different order when `qs` has fewer qubits than the circuit and gets auto-padded. The mathematical contract (which group) is honored but the list-equality of the result depends on insertion order during pad. The `Integration-MethodStabilizer-ExplicitState` test (Tier 4) had to be made set-wise with `Sort` to accommodate this. |
| **Source** | the auto-pad path inside [`Stabilizer/GateUpdates.m`](../../../QuantumFramework/Kernel/Stabilizer/GateUpdates.m) (gates that operate on a qubit beyond the current register size) |
| **Why deferred** | Cosmetic / API consistency. The group is correct — only the human-readable generator order differs from the Phase 1 convention. |
| **Reference** | None — internal convention. |
| **Next step** | Either (a) make the auto-pad insert new-qubit rows in canonical AG order matching the integer-constructor path, or (b) document a single canonical `Sort`-able form for `["Stabilizers"]` and apply it on output. (b) is simpler; (a) is more aesthetic but touches the gate-update code. |
| **File** | [`Stabilizer/GateUpdates.m`](../../../QuantumFramework/Kernel/Stabilizer/GateUpdates.m) (option a) or [`Stabilizer/Properties.m`](../../../QuantumFramework/Kernel/Stabilizer/Properties.m) (option b) |
| **Tests to add** | revert `Integration-MethodStabilizer-ExplicitState` to list-equality once order is canonical |
| **Effort** | ~50 LOC + 0 new tests |

### A.11 — `pauliStringMatrix` fails on single-qubit input
| | |
|---|---|
| **Current** | `StabilizerExpectation[PauliStabilizer[1], "Z"]` (and `"X"`, `"Y"`) on a 1-qubit input falls into the direct-vector fallback path, where `pauliStringMatrix` calls `KroneckerProduct @@ {single}`. With a one-element list, `KroneckerProduct` throws `KroneckerProduct::argmu` and the result becomes an unevaluated `Re[…KroneckerProduct[…]…]` instead of `1`. |
| **Source** | [`Stabilizer/InnerProduct.m:104-113`](../../../QuantumFramework/Kernel/Stabilizer/InnerProduct.m) |
| **Why deferred** | Surfaced 2026-05-04 during refpage validation; not blocking multi-qubit usage but breaks the simplest 1-qubit example. |
| **Reference** | None — kernel implementation detail. |
| **Next step** | Wrap the single-element list case: `If[Length[mats] == 1, First[mats], KroneckerProduct @@ mats]`. Or use `Apply[Dot, mats]` style fallback. |
| **File** | [`Stabilizer/InnerProduct.m:104-113`](../../../QuantumFramework/Kernel/Stabilizer/InnerProduct.m) |
| **Tests to add** | `Phase4-Expectation-OneQubit-Z` and similar 1-qubit cases (need to be added; previously only multi-qubit tested). |
| **Effort** | ~5 LOC + 6 tests |

### A.12 — `ps["Stabilizers"]` formatter throws `StringJoin::string` on symbolic phases
| | |
|---|---|
| **Current** | After a `StabilizerMeasure`, the resulting `PauliStabilizer` has phase entries like `\[FormalS][k]` or `1 - 2 \[FormalS][k]`. Reading `ps["Stabilizers"]` then calls `PauliForm`, which tries `StringJoin[1 - 2 \[FormalS][k], "Z"]` and fails with `StringJoin::string`. The user can still read `ps["Phase"]` but the human-readable string list is broken until `SubstituteOutcomes` is applied. |
| **Source** | [`Stabilizer/Formatting.m`](../../../QuantumFramework/Kernel/Stabilizer/Formatting.m) `PauliForm` — sign-prefix replacement `Replace[signs, {1 -> "", -1 -> "-"}, {1}]` doesn't handle symbolic signs. |
| **Why deferred** | Surfaced 2026-05-04 during refpage validation. The `["Phase"]` accessor is the documented escape hatch for symbolic states; doc pages now show that path. |
| **Reference** | None — formatter detail. |
| **Next step** | Extend the sign-prefix replacement to handle non-numeric signs: produce `ToString[s]` or wrap symbolic signs in a `Row[{s, "·", paulistring}]` form. |
| **File** | [`Stabilizer/Formatting.m`](../../../QuantumFramework/Kernel/Stabilizer/Formatting.m) |
| **Tests to add** | `Phase3-Stabilizers-FormatsSymbolic` — assert `ps["Stabilizers"]` does not throw `StringJoin::string` after a `StabilizerMeasure`. |
| **Effort** | ~10 LOC + 2 tests |

### A.13 — `GraphState[ps_PauliStabilizer]` silently returns edgeless graph for non-graph-form input
| | |
|---|---|
| **Current** | `GraphState[PauliStabilizer[{"XXX", "ZZI", "IZZ"}]]` (a GHZ stabilizer set) silently returns an edgeless 3-vertex graph. The constructor only handles graph-form stabilizers (`X_i ⊗ Π_{j ∈ N(i)} Z_j`); other inputs return an empty edge set without warning. Users discovering this think the conversion succeeded. |
| **Source** | [`Stabilizer/GraphState.m`](../../../QuantumFramework/Kernel/Stabilizer/GraphState.m) — `GraphState[ps_PauliStabilizer]` constructor. |
| **Why deferred** | Surfaced 2026-05-04 during refpage validation. The user-facing escape is "use `LocalComplement` to convert non-graph stabilizer sets to graph form first" but that's not documented and not always possible. |
| **Reference** | AndBri05 §2 (graph-state stabilizer form). |
| **Next step** | Either (a) emit a `GraphState::nongraph` message when the input isn't graph-form, returning `$Failed`; or (b) attempt graph-form conversion via local Cliffords (AndBri05 Lemma 1). (a) is the immediate-fix path. |
| **File** | [`Stabilizer/GraphState.m`](../../../QuantumFramework/Kernel/Stabilizer/GraphState.m) |
| **Tests to add** | `Phase5-GraphState-NonGraphForm-Fails` — assert `GraphState[PauliStabilizer["GHZ-stabilizers"]]` issues `::nongraph` and returns `$Failed`. |
| **Effort** | ~15 LOC + 2 tests |

---

## B. Deferred features (not started)

### B.1 — Quadratic-form triple `(V, Q, ℓ)` (DehMoo03)
| | |
|---|---|
| **Reference** | DehMoo03 §4 Theorem 3 (arxiv:quant-ph/0304125); HosDehMoo04 §V Theorem 1 (qudit generalization); deSilSalYin23 §3 (fast algorithms) |
| **Why deferred** | Niche representation; useful for symbolic computation but not on the v1 critical path. |
| **What to build** | `QuadraticForm[<\|"Subspace" -> V, "Quadratic" -> Q, "Linear" -> ℓ\|>]` head + bidirectional conversion `StabilizerToQuadraticForm[ps]` ↔ `QuadraticFormToStabilizer[…]`. |
| **File** | new `Stabilizer/QuadraticForm.m` |
| **Tests** | Round-trip closure for n = 3, 5; symbolic Q after Clifford action; closed-form amplitudes match direct computation for n ≤ 8. |
| **Effort** | ~300 LOC + 12 tests |

### B.2 — Clifford channel via Choi tableau (Yashin25) — Phase 8
| | |
|---|---|
| **Reference** | Yashin25 §2.3 (arxiv:2504.14101) |
| **Status** | **Phase 8.1 + 8.2 + 8.3 DONE (2026-05-06)**: head, predicate, basic constructors, accessors, **composition via Boolean null space + AG row-sum phase tracking + Φ⁺ contraction sign** (Yashin25 §3.2/§3.3), `cc[ps_PauliStabilizer]` state evolution for identity / state-prep / general dim-matched channels, `CliffordChannel[qc_QuantumChannel]` for deterministic-Pauli channels (stochastic emits notice). |
| **Files** | [`Stabilizer/CliffordChannel.m`](../../../QuantumFramework/Kernel/Stabilizer/CliffordChannel.m) (head, predicate, ctors, composition with full phase tracking, state evolution, QC interop). |
| **Tests** | [`Tests/CliffordChannel.wlt`](../../../Tests/CliffordChannel.wlt) — 76 tests: TIER A-D (predicate + basic ctors), E-I (extended ctors / properties / invariants), J-O (Phase 8.2 composition: identity-on-state, associativity, state-evolution physics, QC roundtrip, random-Clifford-rank), TIER P-T (Phase 8.3 row-sum AG phase + contraction sign + correctness on S/H/X/Y compositions including S² = Z, S³ = S⁻¹ X→-Y, H S H Y → ±Y). |
| **Algorithm (Phase 8.2/8.3 composition)** | For `cc1[cc2]` (compose `cc1 ∘ cc2`), stack `cc2.UB` (k1 rows) above `cc1.UA` (k2 rows) into a `(k1+k2) × 2nB` matrix. Compute the *left* null space `λ s.t. λ · stack = 0` via `NullSpace[Transpose[stack], Modulus -> 2]`. Each kernel vector `λ = (λ_A, λ_C)` produces a composition row `[λ_A · cc2.UA, λ_C · cc1.UB, c_new]` mod 2 where `c_new = (λ · cConcat) ⊕ ((rowSumPhase + contractionPhase) / 2)`. The row-sum phase tracks the AG g-function across F₂-summing `λ_A · UA`, `λ_C · UC`, and both `λ · UB` sides; the contraction phase adds `(-1)^{∑_q x_q z_q}` for the combined u_B Pauli (Y_B contributes -1 since Y^T = -Y). Tested on `cc_S² = cc_Z` (X-row picks up sign), `cc_S³ = cc_S⁻¹` (X→-Y), Heisenberg picture H S H Y → -Y. |
| **Helpers** | `stabilizerRowSumAGPhase[U, lambda, n]` (PackageScope) — tracks cumulative AG i-power mod 4 of the F₂-summed Pauli rows. `stabilizerContractionPhase[uB, n]` (PackageScope) — returns `(-1)^{∑_q x_q z_q} → 0 or 2 mod 4` for the Φ⁺ contraction sign. |

### B.3 — 24-element LocalClifford group (AndBri05)
| | |
|---|---|
| **Reference** | AndBri05 §2 footnote; And05 (the supplementary table — see external-packages-audit.md) |
| **Why deferred** | Required for richer `GraphState` algorithms (item A.6) but not blocking v1 cluster-state demos. |
| **What to build** | `LocalCliffordGroup[]` (list of 24 2×2 unitaries over `ℤ[i]`); `LocalCliffordIndex[m]` (matrix → index 0..23); `LocalCliffordCompose[a, b]` (24×24 composition table). |
| **File** | new `Stabilizer/LocalClifford.m` |
| **Tests** | Group closure (`g_i · g_j ∈ {g_k}`); identity = index 0; each element has an inverse; `LocalCliffordCompose[a, LocalCliffordIndex[Inverse[…]]]` = 0. |
| **Effort** | ~250 LOC + 8 tests |

### B.4 — ~~`CliffordTableau` head distinct from `PauliStabilizer`~~ (DROPPED 2026-05-06)

**Superseded by Phase 8 (Yashin25 Choi tableau, [B.2](#b2--clifford-channel-via-choi-tableau-yashin25)).** The original B.4 motivation was to disambiguate "state's tableau" (`PauliStabilizer`) from "gate's tableau" (proposed `CliffordTableau`). Phase 6 design review surfaced that adding *another* head goes the wrong direction — the principled answer is the Yashin25 unifier, which subsumes states, gates, channels, and measurements into a single Choi-tableau type. Keep the historical entry for trace; do not implement.

### B.5 — `IndexClifford` (KoeSmo14 §3.3 inverse map)
| | |
|---|---|
| **Reference** | KoeSmo14 §3.3 (arxiv:1406.2170) |
| **Why deferred** | Useful for hashing and dedup; not on critical path. |
| **What to build** | `IndexClifford[ct]` = integer index `0 ≤ i < |C_n|` such that `RandomClifford[n][i]` reproduces `ct`. Inverse of the Mallows construction. |
| **File** | extend [`Stabilizer/RandomClifford.m`](../../../QuantumFramework/Kernel/Stabilizer/RandomClifford.m) |
| **Tests** | `IndexClifford[RandomClifford[3]]` is in `[0, 92897280)`; `RandomClifford[3, IndexClifford[ct]] === ct`. |
| **Effort** | ~100 LOC + 4 tests |

### B.6 — Pauli tracking (Paler14, RuhDev25)
| | |
|---|---|
| **Reference** | Paler14 Tables 1–3 (arxiv:1401.5872); RuhDev25 §3 (arxiv:2405.03970) |
| **Why deferred** | Required for QF MBQC pipeline. Complex: per-qubit `{I, X, Z, XZ}` frame + measurement-induced byproduct propagation + scheduling constraints. |
| **What to build** | `PauliFrame[{...}]` head + `PauliTrack[circuit, frame]` propagation per Paler14 update tables. |
| **File** | new `Stabilizer/PauliTracking.m` |
| **Tests** | Surface-code patch trace; MBQC measurement-order partial-order verification. |
| **Effort** | ~300 LOC + 10 tests |

### B.7 — Symbolic Clifford parameters + DLA (Mueller26)
| | |
|---|---|
| **Reference** | Mueller26 §3 (arxiv:2601.02233) |
| **Why deferred** | Variational ansatz support; useful for VQE/QAOA differentiability but separate concern from stabilizer-state simulation. |
| **What to build** | `ParametricPauliRotation[P, θ]` symbolic operator + parameter-shift gradient + nested-commutator DLA Lie closure with hash-set-based duplicate detection. |
| **File** | new `Stabilizer/Parametric.m` |
| **Tests** | DLA closure for `{ZZ, XX}` is `{ZZ, XX, YY, II}` (4 elements); parameter-shift gradient matches finite-difference for a small ansatz. |
| **Effort** | ~400 LOC + 12 tests |

### B.8 — Distance / nearest neighbors (GarMarCro12 §5)
| | |
|---|---|
| **Reference** | GarMarCro12 §5 |
| **Why deferred** | Sanity-check utility; not v1 critical. |
| **What to build** | `StabilizerNearestNeighbors[ps]` returning the `4(2ⁿ − 1)` nearest-neighbor stabilizer states with `\|⟨ψ\|φ⟩\| = 2^(-1/2)`. |
| **File** | new `Stabilizer/NearestNeighbors.m` |
| **Tests** | Count = `4(2ⁿ − 1)` for n ∈ {1, 2, 3}; |⟨ψ\|φ⟩| = 2^(-1/2) for each returned state. |
| **Effort** | ~120 LOC + 6 tests |

### B.9 — `EncoderCircuit[code]` from `H_q` (MonPar23)
| | |
|---|---|
| **Reference** | MonPar23 §IV Algorithm 1 (arxiv:2309.11793) |
| **Why deferred** | The current named codes ([5,1,3], Steane, Shor) hard-code stabilizers; an `EncoderCircuit` would synthesize the encoder for *any* user-supplied stabilizer code. |
| **What to build** | `EncoderCircuit[code]` extracting the parity-check matrix `H_q`, putting it in standard form `[I_1, A_1, A_2 \| B, C_1, C_2; 0, 0, 0 \| D, I_2, E]`, then producing the encoding circuit. |
| **File** | new `Stabilizer/EncoderSynthesis.m` |
| **Tests** | `EncoderCircuit[PauliStabilizer["5QubitCode"]]` reproduces a known [5,1,3] encoding circuit. Cross-check against Got97 §6.4 for Steane, Shor, and 5-qubit codes. |
| **Effort** | ~250 LOC + 6 tests |

### B.10 — Fast classical interconversion (deSilSalYin23)
| | |
|---|---|
| **Reference** | deSilSalYin23 §3, §4.3 (arxiv:2311.10357) |
| **Why deferred** | The 10 fast algorithms span amplitudes ↔ quadratic-form ↔ check-matrix interconversion in `O(N · n)` instead of `O(N⁴)`. Implementing them well requires the quadratic-form representation (item B.1). |
| **What to build** | `CliffordTableauFromMatrix[U]` in `O(n · 2ⁿ)` (deSilSalYin23 §4.3) — exponentially faster than the current `4ⁿ` tomography path. |
| **File** | new `Stabilizer/FastInterconversion.m` |
| **Tests** | Cross-check vs. existing 4ⁿ paths for n ≤ 8; verify n = 10, 12 work where 4ⁿ OOMs. |
| **Effort** | ~300 LOC + 10 tests |

### B.11 — Cluster-state graph rule book (PatGuh26)
| | |
|---|---|
| **Reference** | PatGuh26 §3, §4 (arxiv:2312.02377) |
| **Why deferred** | Patil & Guha derive a *graphical* (not tableau) rule book for X/Y/Z measurements on cluster states, plus fusion operations. Phase 4's `ps["M", "XZZXI"]` works at the tableau level, but the graph-rewrite forms (PatGuh26 §3, §4) are not yet a separate `ClusterStateRules.wl` module. |
| **What to build** | Graph-rewrite rules for X, Y, Z measurements on cluster-state vertices, with explicit graph transformations + sign updates. |
| **File** | new `Stabilizer/ClusterStateRules.m` |
| **Tests** | Each measurement type's rule reproduces the underlying tableau measurement on a small cluster state. |
| **Effort** | ~400 LOC + 12 tests |

### B.12 — Wigner ↔ tableau equivalence (KocHuaLov17, odd `d`)
| | |
|---|---|
| **Reference** | KocHuaLov17 (arxiv:1703.04630) |
| **Why deferred** | Ties to qudit unification (item B.13). For odd qudit dimensions, AarGot04 simulation = discrete Wigner phase-space simulation. |
| **What to build** | `QuantumWignerTransform` overload for `PauliStabilizer` that returns the tableau's discrete Wigner representation. |
| **File** | extend `QuantumWignerTransform.m` (in the broader QF kernel) |
| **Tests** | Deferred until qudit unification lands. |
| **Effort** | ~150 LOC + 5 tests, requires B.13 first |

### B.13 — Symbolic qudits (HosDehMoo04, Beaudrap11, WinPay24b) — v2 unification
| | |
|---|---|
| **Reference** | HosDehMoo04, Beaudrap11, WinPay24b |
| **Why deferred** | Biggest "comparative advantage" in synthesis §3.6 but requires a substantial qudit-aware refactor of every Pauli operation. |
| **What to build** | Refactor `Stabilizer/` to operate on `Z_d` (or `Z_{2d}` for even `d`) instead of hard-coding `Z_2`. Use WinPay24b's "fake ½" `Z_d`-module trick to unify even/odd parity. The internal `PauliRow` already takes a `d` parameter — most other code paths need generalization. |
| **File** | refactor across `Stabilizer/` |
| **Tests** | Qutrit Bell-state stabilizers; qudit Clifford composition; `<Bell\|YY\|Bell>` analog for `d = 3`. |
| **Effort** | ~1000 LOC refactor + 30 tests |

---

## C. Tracking convention

Each item above has:
- **ID** (e.g., `A.1`, `B.7`).
- **Current state** — what works today.
- **Source pointer** — file:line if applicable.
- **Why** — reason for the gap.
- **Reference** — paper anchor.
- **Next step** — concrete algorithm sketch.
- **File** — where the new code goes.
- **Tests to add** — what tier and what to assert.
- **Effort** — rough LOC + test count.

When an item is implemented:
1. Update the relevant section of `synthesis-implementation.md` (move from `partial`/`⏸` to `✅`).
2. Re-run `verify-synthesis-implementation.wls` and update embedded outputs.
3. Add the new tests to `Tests/PauliStabilizer.wlt`.
4. Mark the corresponding item in this document as **DONE** with a commit hash.
5. Promote any new public symbols in `PacletInfo.wl` and add `Usage.m` entries.

## D. Cross-references

- Companion: [`synthesis-implementation.md`](synthesis-implementation.md) — the *what works* document (capability tour).
- Companion: [`API.md`](API.md) — the per-function reference (49 verified examples).
- Test suite: [`Tests/PauliStabilizer.wlt`](../../../Tests/PauliStabilizer.wlt) — 4 top-level tiers, 32+ sub-tiers, 250 tests.
- Verifier (synthesis): [`verify-synthesis-implementation.wls`](verify-synthesis-implementation.wls) — re-runnable with `wolframscript`.
- Verifier (API): [`verify-API.wls`](verify-API.wls) — re-runnable with `wolframscript`.
- Original synthesis: [`OngoingProjects/Stabilizer/package-design-synthesis.md`](../package-design-synthesis.md) — distilled from 28 papers.
- Plan: `/Users/mohammadb/.claude/plans/audit-this-users-mohammadb-documents-git-robust-russell.md` (local-only).
