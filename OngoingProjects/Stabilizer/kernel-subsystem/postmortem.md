# Phase 5c post-mortem: the Y round-trip miss

> **TL;DR.** TIER 1.4 in `Tests/Stabilizer/PauliStabilizer.wlt` was labeled "Round‑trips" but did not contain a single QO→PS→QO or QS→PS→QS exact‑equality test. `PauliStabilizer[QuantumOperator["Y"]]["QuantumOperator"]` was returning `i·Y`, not `Y`, and 185/185 tests went green. The miss was structural, not careless: an unstated design choice was inherited from the stabilizer-simulator tradition.
>
> **Resolution.** Phase 5c made the design choice explicit (track global phase) and finished implementing it. This document records the *why* so the same class of mistake doesn't recur.

Anchor: branch `stabilizer-phases-1-4`, commits `93edc400` (red tests), `6f8637ee` (state tomography), `f7ebfd56` (operator phase capture), `57d8b53f` (5c docs), plus the gate-update propagation commit closing §A.9.

> **Phase 6 footnote (2026-05-06).** After Phase 5c shipped, a design review with N. Murzin flagged the public API surface (10 symbols across Phases 3–5) as excessive. Phase 6 demoted four of them (`RandomClifford`, `StabilizerMeasure`, `SubstituteOutcomes`, `SampleOutcomes`) to method-grade operations on `PauliStabilizer`. **Phase 6.5** then demoted the two borderline cases (`StabilizerInnerProduct`, `StabilizerExpectation`) on the same grounds — both are single-receiver operations that read more naturally as `ps["InnerProduct", other]` and `ps["Expectation", pauli]`. Combined surface: 10 → 4 (`PauliStabilizer`, `StabilizerFrame`, `GraphState`, `LocalComplement`). No semantic change. The same review argued for routing hybrid stabilizer ↔ Schrödinger interop through cross-head **UpValues** rather than `Picture -> "Stabilizer"` or `QuantumBasis`-wrapping (both would force `O(2ⁿ)` materialization through the basis pipeline and defeat the formalism's `O(n²)` advantage). See `~/.claude/plans/i-am-conviced-of-vast-kettle.md` for the design plan and the audit findings (Picture is a guard, not metadata; `QuantumBasis` and `QuantumState` data shapes don't admit side-car keys without breaking validation; UpValues are already the proven cross-type pattern at [Stabilizer/Conversions.m:139-142](../../../QuantumFramework/Kernel/Stabilizer/Conversions.m)). Phase 7 (UpValues for `qmo[ps]`, `qc[ps]`, etc.) and Phase 8 (Yashin25 Choi-tableau unifier per [ROADMAP §B.2](roadmap.md#b2--clifford-channel-via-choi-tableau-yashin25)) are the next milestones, scoped to separate PRs. Process rule going forward: every feature-growth phase ends with an "API surface consolidation" pass.

---

## 1. What was missed

Before Phase 5c, the only tests labeled "Round‑trips" were five `Roundtrip‑*` entries at [`Tests/Stabilizer/PauliStabilizer.wlt:224‑263`](../../../Tests/Stabilizer/PauliStabilizer.wlt):

| Test | What it actually tested |
|---|---|
| `Roundtrip-CircuitToBellStabilizers` | one-way: QCO → stabilizer-string list |
| `Roundtrip-PhaseSigns` | internal formula `Signs = 1 − 2·Phase` |
| `Roundtrip-MatrixIdempotent` | `<\|"Matrix" -> m\|>` constructor idempotency on `m` |
| `Roundtrip-HSquaredIsIdentity` | `H · H = I` gate-closure |
| `Roundtrip-SFourthIsIdentity` | `S⁴ = I` gate-closure |

None of these is `A[C[x]] === x` for any constructor `C` paired with any accessor `A`. The label looked thorough; the contents didn't exercise the round‑trip the user would actually try.

The four most obvious user examples that should have been there but weren't:

1. `PauliStabilizer[QuantumOperator["Y"]]["QuantumOperator"] === QuantumOperator["Y"]` — was returning `i·Y` (off by a factor of `i` because the AG decomposition recovers `Z·X = i·Y`).
2. `PauliStabilizer[QuantumOperator["YY"]]["QuantumOperator"] === QuantumOperator["YY"]` — was returning `−Y⊗Y` (`i² = −1`).
3. `PauliStabilizer[QuantumState[{0, 1}]]["State"] === QuantumState[{0, 1}]` — was returning `\|0⟩` (the state-tomography path lost the stabilizer sign and reconstructed `+Z`'s eigenstate instead of `−Z`'s).
4. `PauliStabilizer[QuantumState["Bell"]]["Stabilizers"] === {"XX", "ZZ"}` — was returning `{"−XX", "YY"}` (the same `RowSpace`-canonicalization sign loss returned the stabilizer set of `β₁₀`, not `β₀₀`).

---

## 2. Design rationale: phase-aware vs phase-oblivious

### Pure stabilizer simulators

Phase-oblivious is the right call. The whole point of [Stim](https://github.com/quantumlib/Stim), Aaronson–Gottesman, [QuantumClifford.jl](https://github.com/QuantumSavory/QuantumClifford.jl) is to simulate measurement statistics on millions of qubits in `O(n²)` memory. Measurement statistics are phase-invariant: `|ψ⟩` and `e^{iα}|ψ⟩` give identical outcome probabilities for any observable. Tracking phase is wasted bytes against an unobservable quantity. The canonical view in this tradition is "the Clifford group acts on the Pauli group by conjugation" — under conjugation, `Y` and `i·Y` are the same Clifford gate.

### QuantumFramework

Phase-aware is the right call. Three concrete reasons:

1. **`===` in WL means exact identity.** Every other QF module — `QuantumOperator`, `QuantumState`, `QuantumChannel` — compares as literal matrices/vectors. A user writes `PauliStabilizer[qo]["QuantumOperator"] === qo` and expects a boolean answer that respects the actual matrix. Treating `i·Y` and `Y` as the same operator violates the contract every other module honors.
2. **`PauliStabilizer` is one method among many.** Users compose `someNonCliffordOp @ pauliStabilizerOp`. The relative phase of the stabilizer factor is observable in that combined object, even when it would be unobservable for the stabilizer factor alone.
3. **Didactic / book context.** *Quantum Computing by Model and Simulation* makes notational distinctions between `Y` and `i·X·Z`; users following along expect the framework to preserve them. "Up to global phase" is fine in a paper; it isn't fine in a pedagogy where every formula round-trips through the kernel.

The cost: one complex number per `PauliStabilizer` object (stored under the `"GlobalPhase"` association key). Negligible against the `O(n²)` tableau advantage.

**Phase 5c made this choice explicit.** Past sessions had inherited the simulator convention without flagging that QF wants the stricter contract; tests were therefore written under the simulator convention.

---

## 3. Why the miss was structurally easy

Three reasons, each with a code anchor:

1. **Mislabeled tier.** [`Tests/Stabilizer/PauliStabilizer.wlt:224‑263`](../../../Tests/Stabilizer/PauliStabilizer.wlt) had a tier called "Round‑trips" that didn't contain round-trips. Reviewers seeing "TIER 1.4 — Round-trips ✅ all green" had no signal that the round-trip a user would actually try wasn't covered.
2. **Inherited convention.** Past test-writing reproduced the stabilizer-simulator convention of phase-oblivious comparison (see §2 above) without flagging it. There was no statement anywhere that QF's contract was stricter than the convention.
3. **The verifier scripts (now retired) reproduced the same blind spot.** The legacy `verify-API.wls` (deleted 2026-05-07) explicitly avoided exact equality, using `Chop @ Abs[Conjugate[vec1] . vec2] − 1` with the comment *"String equality of stabilizer generators does NOT hold"*. The legacy `verify-synthesis-implementation.wls` (also deleted) tested exact equality but only for stabilizer-set roundtrips that stay inside `PauliStabilizer`; it never crossed the constructor↔accessor type boundary. Their replacement is `Tests/Stabilizer/AuditMatrix.wlt` plus the existing `VerificationTest`-based `.wlt` files, which DO cross type boundaries (TIER 15 in `AuditMatrix.wlt`).

The miss was structural. A single more-careful pass wouldn't have caught it; a different test-writing rule would have.

---

## 4. What changed in Phase 5c

| Commit | Step | Outcome |
|---|---|---|
| `93edc400` | TIER 1.4a/b/c added (42 tests) | 21 red — locked the spec |
| `6f8637ee` | State tomography rewrite (`Stabilizer/Constructors.m`); 4ⁿ Pauli expectations + greedy F₂ rank growth + symplectic Gram-Schmidt; introduced `"GlobalPhase"` association key | 21 → 6 red |
| `f7ebfd56` | Operator global-phase capture (`Stabilizer/Constructors.m`, `Stabilizer/Conversions.m`) | 6 → 0 red |
| `57d8b53f` | Doc updates: roadmap / api.md reflect 5c | docs aligned |
| `f9e2d711` | This post-mortem; A.9 contract documented (`Stabilizer/GateUpdates.m` comment block); TIER 1.4d added (up-to-phase + escape-hatch tests) | 227 → 235 |
| `cbf51576` | Phase C cross-module audit — `Tests/Stabilizer/Roundtrips.wlt` probes `QuantumChannel ↔ QuantumOperator` and `QuantumCircuitOperator ↔ QuantumOperator` | 20/20 PASS, no new bugs |

Test count: 185 → 235 on `PauliStabilizer.wlt`; 217 → 287 across the suite (235 + 32 + 20).

---

## 5. What's still partial after Phase 5c

- **§A.9** — `"GlobalPhase"` does not propagate through *gate updates*. Phase 5c implemented exact phase tracking on the *first hop* (constructor → accessor); gate-update propagation was investigated 2026-05-04 and found to require `O(2ⁿ)` materialization at every gate update because the phase factor a Clifford gate picks up depends on the input state's tableau, not just on the gate. (Concrete counterexample: `Z\|0⟩ = \|0⟩` vs `Z\|1⟩ = −\|1⟩`.) That cost defeats the AG `O(n²)` complexity advantage. **Reclassified as an inherent trade-off with a documented contract**, not a closeable bug. The contract: `PauliStabilizer[qs]["gate", q]["State"]` equals `gate @ qs` *up to global phase*; for exact equality, re-run the constructor on the post-gate state. TIER 1.4d in [`Tests/Stabilizer/PauliStabilizer.wlt`](../../../Tests/Stabilizer/PauliStabilizer.wlt) documents this with up-to-phase + escape-hatch tests.
- **§A.10** — generator order differs between default-register and QS-derived paths after auto-pad. The mathematical group is correct; only the list order differs (which forced `Integration-MethodStabilizer-ExplicitState` to use `Sort` on both sides). Cosmetic; not on the v1 critical path.
- **Cross-module audit (Phase C)** — [`Tests/Stabilizer/Roundtrips.wlt`](../../../Tests/Stabilizer/Roundtrips.wlt) probed `QuantumChannel ↔ QuantumOperator` and `QuantumCircuitOperator ↔ QuantumOperator`. **20/20 PASS — no new bugs found.** The Phase 5c Y bug was unique to the AG-tableau projection.

---

## 6. Lessons (mirrors the global memory rule)

1. **User-facing exact-equality round-trips first.** When designing tests for any constructor-accessor pair `C[x]` / `A[ps]`, the *first* test must be `A[C[x]] === x` for representative simple `x`. Phase-oblivious / set-equivalence comparisons come *after*, opt-in via clearly-named helpers.
2. **State the design choice explicitly.** When inheriting a convention from a tradition (stabilizer simulators, MPS canonical forms, etc.), check whether the convention matches QF's contract. If it doesn't, write down the divergence in the relevant module doc — not just in a code comment.
3. **Verifiers must cross type boundaries.** A verifier that only checks invariants *inside* a type ("the stabilizer-set roundtrips through `Circuit`") will not catch contract violations *across* types ("the `QuantumOperator` recovered from a `PauliStabilizer` is not the original `QuantumOperator`"). At least one test per accessor-pair must compare against the input as the user would.

---

## 7. References

- ROADMAP entry: [Phase 5c — DONE](roadmap.md#phase-5c--done-2026-05-04) (commit hashes).
- Global memory rule: `~/.claude/projects/-Users-mohammadb-Documents-GitHub-QuantumFramework/memory/feedback_user_facing_roundtrip_first.md`.
- Cross-module probe: [`Tests/Stabilizer/Roundtrips.wlt`](../../../Tests/Stabilizer/Roundtrips.wlt) (Phase C of the plan).
- Original audit: [`paulistabilizer-source-audit.md`](paulistabilizer-source-audit.md).
- Source synthesis: [`package-design-synthesis.md`](package-design-synthesis.md).
