# Phase 5c post-mortem: the Y round-trip miss

> **TL;DR.** TIER 1.4 in `Tests/PauliStabilizer.wlt` was labeled "Round‑trips" but did not contain a single QO→PS→QO or QS→PS→QS exact‑equality test. `PauliStabilizer[QuantumOperator["Y"]]["QuantumOperator"]` was returning `i·Y`, not `Y`, and 185/185 tests went green. The miss was structural, not careless: an unstated design choice was inherited from the stabilizer-simulator tradition.
>
> **Resolution.** Phase 5c made the design choice explicit (track global phase) and finished implementing it. This document records the *why* so the same class of mistake doesn't recur.

Anchor: branch `stabilizer-phases-1-4`, commits `93edc400` (red tests), `6f8637ee` (state tomography), `f7ebfd56` (operator phase capture), `57d8b53f` (5c docs), plus the gate-update propagation commit closing §A.9.

---

## 1. What was missed

Before Phase 5c, the only tests labeled "Round‑trips" were five `Roundtrip‑*` entries at [`Tests/PauliStabilizer.wlt:224‑263`](../../../Tests/PauliStabilizer.wlt):

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

1. **Mislabeled tier.** [`Tests/PauliStabilizer.wlt:224‑263`](../../../Tests/PauliStabilizer.wlt) had a tier called "Round‑trips" that didn't contain round-trips. Reviewers seeing "TIER 1.4 — Round-trips ✅ all green" had no signal that the round-trip a user would actually try wasn't covered.
2. **Inherited convention.** Past test-writing reproduced the stabilizer-simulator convention of phase-oblivious comparison (see §2 above) without flagging it. There was no statement anywhere that QF's contract was stricter than the convention.
3. **The verifier scripts reproduced the same blind spot.** [`OngoingProjects/Stabilizer/Documentation/verify-API.wls:242‑256`](verify-API.wls) explicitly avoids exact equality (uses `Chop @ Abs[Conjugate[vec1] . vec2] − 1`, comment on line 245: *"String equality of stabilizer generators does NOT hold"*). [`verify-synthesis-implementation.wls:98‑108, 324‑338`](verify-synthesis-implementation.wls) does test exact equality but only for stabilizer-set roundtrips that stay inside `PauliStabilizer`; it never crosses the constructor↔accessor type boundary.

The miss was structural. A single more-careful pass wouldn't have caught it; a different test-writing rule would have.

---

## 4. What changed in Phase 5c

| Commit | Step | Outcome |
|---|---|---|
| `93edc400` | TIER 1.4a/b/c added (42 tests) | 21 red — locked the spec |
| `6f8637ee` | State tomography rewrite (`Stabilizer/Constructors.m`); 4ⁿ Pauli expectations + greedy F₂ rank growth + symplectic Gram-Schmidt; introduced `"GlobalPhase"` association key | 21 → 6 red |
| `f7ebfd56` | Operator global-phase capture (`Stabilizer/Constructors.m`, `Stabilizer/Conversions.m`) | 6 → 0 red |
| `57d8b53f` | Doc updates: ROADMAP/API/synthesis-implementation reflect 5c | docs aligned |
| (next commit) | Phase E — gate-update propagation closes ROADMAP §A.9 | gate-update round-trips also exact |

Test count: 185 → 227+ on `PauliStabilizer.wlt`; 217 → 274+ across the suite.

---

## 5. What's still partial

After Phase 5c + this plan's gate-update commit:

- **§A.10** — generator order differs between default-register and QS-derived paths after auto-pad. The mathematical group is correct; only the list order differs (which forced `Integration-MethodStabilizer-ExplicitState` to use `Sort` on both sides). Cosmetic; not on the v1 critical path.
- **Cross-module audit findings** (Phase C of the plan) — any new partial items surfaced by [`Tests/Roundtrips.wlt`](../../../Tests/Roundtrips.wlt) on the QS↔QO and QChannel↔QO pairs are filed in their respective module's status docs.

---

## 6. Lessons (mirrors the global memory rule)

1. **User-facing exact-equality round-trips first.** When designing tests for any constructor-accessor pair `C[x]` / `A[ps]`, the *first* test must be `A[C[x]] === x` for representative simple `x`. Phase-oblivious / set-equivalence comparisons come *after*, opt-in via clearly-named helpers.
2. **State the design choice explicitly.** When inheriting a convention from a tradition (stabilizer simulators, MPS canonical forms, etc.), check whether the convention matches QF's contract. If it doesn't, write down the divergence in the relevant module doc — not just in a code comment.
3. **Verifiers must cross type boundaries.** A verifier that only checks invariants *inside* a type ("the stabilizer-set roundtrips through `Circuit`") will not catch contract violations *across* types ("the `QuantumOperator` recovered from a `PauliStabilizer` is not the original `QuantumOperator`"). At least one test per accessor-pair must compare against the input as the user would.

---

## 7. References

- ROADMAP entry: [Phase 5c — DONE](ROADMAP.md#phase-5c--done-2026-05-04) (commit hashes).
- Global memory rule: `~/.claude/projects/-Users-mohammadb-Documents-GitHub-QuantumFramework/memory/feedback_user_facing_roundtrip_first.md`.
- Cross-module probe: [`Tests/Roundtrips.wlt`](../../../Tests/Roundtrips.wlt) (Phase C of the plan).
- Original audit: [`paulistabilizer-source-audit.md`](../paulistabilizer-source-audit.md).
- Source synthesis: [`package-design-synthesis.md`](../package-design-synthesis.md).
