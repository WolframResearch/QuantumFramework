# Stabilizer subsystem roadmap

> Status tracker for items that are **partial**, **deferred**, or have **latent bugs** in the Stabilizer subsystem (Phases 1–5). Each entry has a concrete next step (file path, signature, algorithm sketch, test to add). Updated: 2026-04-30 (last content update); moved 2026-05-02 to `OngoingProjects/Stabilizer/Documentation/`. Branch: `stabilizer-phases-1-4`.

> **Audit-doc context.** This document complements [`synthesis-implementation.md`](synthesis-implementation.md) (the *what works today* tour) and [`API.md`](API.md) (per-function reference) by recording *what doesn't yet work, and exactly how to finish it*. The source synthesis is at [`OngoingProjects/Stabilizer/package-design-synthesis.md`](../package-design-synthesis.md).

## Overall status

- **Tests:** 185 / 185 PauliStabilizer + 32 / 32 QuantumDistance = 217 / 217 passing.
- **Phases done:** 1 (refactor + tests), 2 (hygiene), 3 (symbolic phases), 4 (frame + inner products + Pauli measurement), 5a (graph state + LC), 5b (companion MD).
- **Open items:** 16 (8 partial + 7 deferred + 1 latent bug).
- **Synthesis priority hits**: 5 / 10 ✅; 2 / 10 ⚠️ partial; 3 / 10 ⏸ deferred (per `package-design-synthesis.md` §11).

---

## A. Partial implementations (work-but-incomplete)

### A.1 — `StabilizerInnerProduct` closed-form (`O(n³)`)
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

### A.2 — `StabilizerExpectation` AG-phase i-factor tracking
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

### B.2 — Clifford channel via Choi tableau (Yashin25)
| | |
|---|---|
| **Reference** | Yashin25 §2.3 (arxiv:2504.14101) |
| **Why deferred** | Most ambitious unification: subsumes pure states, mixed states, measurements, post-selection into a single `[U_A \| U_B \| c]` Boolean matrix. Composition via vector-space intersection. |
| **What to build** | `CliffordChannel[<\|"UA" -> ..., "UB" -> ..., "c" -> ...\|>]` head + `ChannelCompose[Φ, Ψ]` via Gaussian elimination + interop with existing `QuantumChannel`. |
| **File** | new `Stabilizer/CliffordChannel.m` |
| **Tests** | (a) Round-trip `CliffordChannel[QuantumChannel["BitFlip", p]]` recovers the Choi matrix. (b) `ChannelCompose[Φ_BitFlip[p], Φ_BitFlip[q]] == Φ_BitFlip[p + q − 2 p q]`. (c) Pure stabilizer `CliffordChannel[ps]` produces the right `U_A`. |
| **Effort** | ~400 LOC + 15 tests |

### B.3 — 24-element LocalClifford group (AndBri05)
| | |
|---|---|
| **Reference** | AndBri05 §2 footnote; And05 (the supplementary table — see external-packages-audit.md) |
| **Why deferred** | Required for richer `GraphState` algorithms (item A.6) but not blocking v1 cluster-state demos. |
| **What to build** | `LocalCliffordGroup[]` (list of 24 2×2 unitaries over `ℤ[i]`); `LocalCliffordIndex[m]` (matrix → index 0..23); `LocalCliffordCompose[a, b]` (24×24 composition table). |
| **File** | new `Stabilizer/LocalClifford.m` |
| **Tests** | Group closure (`g_i · g_j ∈ {g_k}`); identity = index 0; each element has an inverse; `LocalCliffordCompose[a, LocalCliffordIndex[Inverse[…]]]` = 0. |
| **Effort** | ~250 LOC + 8 tests |

### B.4 — `CliffordTableau` head distinct from `PauliStabilizer`
| | |
|---|---|
| **Reference** | Synthesis §5 (the "Clifford operations" menu) |
| **Why deferred** | Conceptually clear (`PauliStabilizer` = state's tableau; `CliffordTableau` = gate's tableau) but the current architecture conflates them since `PauliStabilizer[qo_QuantumOperator]` happens to be how Clifford gates are tableau-encoded. |
| **What to build** | New head `CliffordTableau` with an `Apply[ct, ps]` method that lifts a gate's tableau onto a state's tableau. |
| **File** | new `Stabilizer/CliffordTableau.m` |
| **Tests** | `CliffordTableau[QuantumCircuitOperator[{"H" -> 1}]][PauliStabilizer[1]]` produces the same result as `PauliStabilizer[1]["H", 1]`. |
| **Effort** | ~150 LOC + 6 tests |

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
- Test suite: [`Tests/PauliStabilizer.wlt`](../../../Tests/PauliStabilizer.wlt) — 8 tiers, 185 tests.
- Verifier (synthesis): [`verify-synthesis-implementation.wls`](verify-synthesis-implementation.wls) — re-runnable with `wolframscript`.
- Verifier (API): [`verify-API.wls`](verify-API.wls) — re-runnable with `wolframscript`.
- Original synthesis: [`OngoingProjects/Stabilizer/package-design-synthesis.md`](../package-design-synthesis.md) — distilled from 28 papers.
- Plan: `/Users/mohammadb/.claude/plans/audit-this-users-mohammadb-documents-git-robust-russell.md` (local-only).
