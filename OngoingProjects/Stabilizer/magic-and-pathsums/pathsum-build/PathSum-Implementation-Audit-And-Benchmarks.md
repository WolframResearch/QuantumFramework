---
Title: 'Path-sum implementation: audit and fresh benchmarks'
Author: investigation in the Wolfram QuantumFramework
Description: 'A concise, precise audit of the symbolic path-sum build in this folder (pathsum.wl, elimination.wl, diagonal.wl, rankcap.wl, harness.wl, benchmarks/, tests/): what it implements, the correctness guarantee, fresh re-run of all four stage tests and the benchmark suite, the limitations and benchmark caveats found, and what is deliberately not built. Every number below was produced this session against the working-tree paclet.'
---

# Path-sum implementation: audit and fresh benchmarks

Audit of the path-sum build in this folder, re-run from scratch on 2026-06-24 against the working-tree paclet `Wolfram/QuantumFramework 2.0.0` (loaded with `PacletDirectoryLoad`), machine `mohammadb3maclap`, 12 cores, Wolfram Language 15.0. All four test suites pass; the benchmark numbers below are this session's, not the ones quoted in [`PathSum-WL-Note.md`](PathSum-WL-Note.md). Conventions: amplitudes are exact in $\mathbb{Z}[\omega, 1/\sqrt2]$, $\omega = e^{i\pi/4}$; equality means `RootReduce[a-b] === 0`; the engine uses the full-angle convention $T = \mathrm{diag}(1, \omega)$.

## 1. What is implemented

| File | Provides | Status |
|---|---|---|
| `pathsum.wl` | The carrier: each qubit value an $\mathbb{F}_2$ multilinear poly in path vars `x[i]`; phase a degree-$\le 3$ $\mathbb{Z}_8$ pseudo-Boolean poly; gates $H, T, T^\dagger, S, Z, CZ, CCZ, CNOT$; the brute-force Gauss-sum amplitude `amp` (the exact oracle and fallback); the single-variable $[\mathrm{HH}]$ sign-cofactor elimination `reduce`. | complete, exact |
| `elimination.wl` | The rewrite engine: classify each internal var as `const` / `hh` ($C=4Q$, $Q$ linear) / `s` ($C\equiv 2 \bmod 4$) / `skip` (magic), fire the proven identity, complete the residual by brute force; degree-$\le 3$ guard; `engineReduce`; confluence-by-replay `confluentQ`; the `Inactive[Sum]` formal carrier on a `HoldAll`, `NonThreadable` head. | complete, exact |
| `diagonal.wl` | CNOT-dihedral (no-$H$) compressor: linear $\mathbb{F}_2$ output map, one `LinearSolve[M, y, Modulus->2]`, one closed-form term $2^{-n/2}\omega^{\psi(y)}$, rank $1$; isolated single-variable Gauss-sum lemmas ($\mathbb{Z}_2$ / $\mathbb{Z}_4$ / odd / quadratic completion). | complete, exact |
| `rankcap.wl` | Exact span-dimension cap on a gate-built `StabilizerFrame`: materialize components coherently via relating Paulis, keep a maximal independent subset, re-express by exact `LinearSolve`. A span-dimension witness ($\le 2^n$), not the optimal stabilizer rank. | complete, exact |
| `harness.wl` | QF dense oracle (`Method -> "Schrodinger"`), exact differential matchers (`RootReduce`), randomized circuit generator. | complete |
| `benchmarks/` , `tests/` | `bench.wl` + `run-benchmarks.wls`; `stage0..3-verify.wls`. | see Sections 2-3 |

The engine handles the gate set $\{H, T, T^\dagger, S, Z, CZ, CCZ, CNOT\}$ only. The phase is carried as an *integer* $\mathbb{Z}_8$ polynomial, so the core engine does **not** carry free symbolic gate angles; the symbolic-angle capability lives in a separate builder, [`../pathsum-blowup-alternatives-demos.wls`](../pathsum-blowup-alternatives-demos.wls), not here.

## 2. Correctness audit (all four suites green, this session)

The build's correctness does **not** depend on the rewriter's completeness: every reduction is checked against both the un-reduced brute-force amplitude and QF dense, exactly. Re-run results:

- **Stage 0** (carrier + reduction): $T$-convention exact; the four regression demos reproduce `{vars before, after, reduced===QF, scale=2^k}` exactly (`H(ZH)^8`: $9\to1$; 2q deep: $10\to2$; magic `H(TH)^8`: $9\to9$; mixed: $5\to3$); builder $\equiv$ QF on 3 circuits; randomized differential suite $9/9$.
- **Stage 1** (rank cap): $T\lvert0\rangle$ $2\to1$, interleaved $8\to2$, diagonal $T_1T_2CZ$ $4\to4$ (no over-collapse), 4-$T$ two-qubit $16\to4$, each reproducing the state exactly ($\lvert\Delta\rvert = 0$).
- **Stage 2** (diagonal compressor + lemmas): `diagAmp === QF` and `psi === diagAmp` on $d_1, d_2, d_3$; rank-1 vs eager frame $2^3$; randomized diagonal $n\le 6$ suite $10/10$; single-variable lemmas $12/12$, with the $\mathbb{Z}_2$ / $\mathbb{Z}_4$ / odd-$T$ classes all correct.
- **Stage 3** (full engine): `engineMatchesQ` on the 5-circuit battery and $10/10$ random; degree-$\le 3$ invariant $15/15$; confluence-by-replay $8/8$; Clifford ladder telescopes (internal residual $\to 0$); magic residue grows exactly with $\#T$; `Inactive[Sum]` carrier with `Activate` $\equiv$ engine $\equiv$ QF.

**Verdict: the implemented engine is exact on everything it is tested against, and the test battery is the right one** (differential against an independent dense oracle, plus invariants and confluence replay).

## 3. Fresh benchmarks

**M1, diagonal single-amplitude** (`diagAmp` vs QF dense one entry, seconds; the genuine speed win, exact, flat in $n$):

| $n$ | gates | `diagAmp` (s) | QF dense (s) |
|---|---|---|---|
| 4 | 12 | 0.0023 | 0.566 |
| 6 | 18 | 0.0033 | 0.915 |
| 8 | 24 | 0.0046 | 1.309 |
| 10 | 30 | 0.0059 | 1.795 |
| 12 | 36 | 0.0072 | 2.902 |
| 14 | 42 | 0.0086 | 6.847 |

**M2, Clifford telescoping** (H/S/CZ ladder, $n=6$; internal residual collapses to $0$, steps linear in depth):

| depth | vars before | vars after | internal residual | steps | reduce (s) |
|---|---|---|---|---|---|
| 1 | 12 | 6 | 0 | 4 | 0.003 |
| 2 | 18 | 6 | 0 | 8 | 0.009 |
| 3 | 24 | 6 | 0 | 11 | 0.014 |
| 4 | 30 | 6 | 0 | 15 | 0.023 |
| 6 | 42 | 6 | 0 | 22 | 0.049 |
| 8 | 54 | 0 | 0 | 32 | 0.091 |

**M3, parametric capability**: one symbolic traversal $0.009$ s vs dense over 17 angles $1.84$ s. Caveat in Section 4.

**Narrow-deep magic, the honest loss** (`H(TH)^k` on one qubit; residual $= k$ magic vars, engine brute-forces $2^k$; from `benchNarrow`, re-run this session):

| $k = \#T$ | internal residual | `engineAmp` (s) | QF dense (s) |
|---|---|---|---|
| 4 | 4 | 0.0016 | 0.345 |
| 8 | 8 | 0.013 | 0.530 |
| 12 | 12 | 0.248 | 0.766 |
| 14 | 14 | 1.175 | 0.888 |
| 16 | 16 | 5.465 | 1.033 |
| 18 | 18 | 25.28 | 1.140 |

**The crossover is at $k \approx 14$** on this machine (engine $1.17$ s vs dense $0.89$ s), not $\approx 16$ as the note states; the engine already loses at $k=14$ and loses by $5\times$ at $k=16$. This is the precise honest statement: when the magic paths $2^k$ exceed the state-space dimension, dense wins.

**M5, representational over-count**: the eager `StabilizerFrame` keeps $2^{\#T}$ components ($2, 4, 8, 16, 32, 64$ for $\#T = 1..6$) where the diagonal compressor is rank $1$.

## 4. Findings and caveats (precise)

**Correctness: no defects found.** The engine is exact wherever tested, and the design makes correctness independent of the rewriter (brute-force fallback). The two items below are about *the rewriter's reach* and *the benchmark harness*, not about wrong answers.

1. **`engineReduce` can stop early on a degree-violating $S$ step.** `step` selects the *first* eligible variable in `ps["vars"]`; if that variable classifies as `s` but firing it would push the phase past degree $3$, `fire` refuses and returns the state unchanged, which `engineReduce` reads as a fixed point and halts, leaving later still-reducible variables un-eliminated. Correctness is preserved (brute force completes the residual), but the reduction is not maximal and is variable-order dependent. A retry-next-variable loop would fix the reach.

2. **`rankcap` picks its basis by numeric rank.** `capFrameRank` chooses the maximal independent subset with `MatrixRank[N[...]]` (numeric), then re-expresses exactly with `LinearSolve`. Robust in practice for $\mathbb{Z}[\omega]$ component vectors and exact on every tested case ($\lvert\Delta\rvert = 0$), but the independence test carries a floating-point tolerance risk; an exact `MatrixRank` over `RootReduce`d entries would remove it.

3. **M3 mixes two angle conventions.** The "engine" column times `PauliStabilizer[1][{"H","P"[\[Theta]],"H"}]["StateVector"]`, which is the `StabilizerFrame` *half-angle* `"P"` ($\mathrm{diag}(1, e^{i\theta/2})$, closed form $\cos(\theta/2)$), while the "dense" column rebuilds the *full-angle* `QuantumCircuitOperator["P"[\[Theta]]]` ($\cos\theta$). The benchmark only times, never compares values, so it does not error, but the two columns compute different functions and the reported closed form is the frame's, not the engine's. The timing illustration (symbolic once vs grid) stands; the clean, convention-consistent capability demonstration is [`../pathsum-blowup-alternatives-demos.wls`](../pathsum-blowup-alternatives-demos.wls) ($\langle Z\rangle(\theta)=\cos\theta$, exact vs full-angle dense).

4. **`run-benchmarks.wls` does not exercise the honest loss.** Its M4 calls `benchMagic[3, d, ...]`, and for $n=3$ those random circuits fully reduced (internal residual $0$), so the engine "wins" there; that row is not representative. The genuine loss regime is `benchNarrow` (Section 3 table), which the runner omits. Recommend adding `benchNarrow` to the runner.

## 5. What is deliberately not built (scope boundary)

- **No confluent $\mathbb{Z}_8$ system.** Only the `const` / $[\mathrm{HH}]$ / $S$ rules; no quadratic completion or multi-variable constraints. Confluence for the $\mathbb{Z}_8$ phase is observed empirically (shuffled-order replay) but not proven; Amy-Stinchcombe prove it only for the $\mathbb{F}_2$ sign fragment.
- **No closed-form Clifford residual.** The quadratic (Clifford) part of the residual is brute-forced, not evaluated by the Guan-Regan $\mathbb{F}_2$-rank Gauss sum.
- **No competitive scaling.** The residual is completed at $2^{\lvert V\rvert}$; the de Colnet rank-width fixed-parameter contraction is not implemented.
- **No GroebnerBasis / variety point-count path.** Named as one call away, not wired into the engine.
- **No frame merging.** `rankcap` caps at the span dimension ($\le 2^n$), not the NP-hard optimal stabilizer rank; no Garcia-Markov term merging.
- **No free symbolic angles in the core engine** (integer $\mathbb{Z}_8$ phases only), and **no general rotations** ($RZ$, $RY$); the symbolic-angle capability is a separate builder.

## 6. Verdict

The implemented path-sum build is a small, exact, well-tested engine: a correct carrier, a sound degree-bounded rewriter with a brute-force safety net, a genuine win on the diagonal fragment (flat in $n$, milliseconds where dense grows) and on Clifford telescoping, a verified rank cap, and an honest, reproduced loss on narrow-deep magic with the crossover at $k\approx 14$. Correctness is solid; the open work is reach and scaling (maximal $S$-reduction, closed-form Clifford residual, rank-width contraction, the variety path), and two harness fixes (the M3 convention mix and adding `benchNarrow` to the runner).

*Re-run: `tests/stage0..3-verify.wls` and `benchmarks/run-benchmarks.wls` (plus a one-off `benchNarrow` for Section 3), this session, against `Wolfram/QuantumFramework 2.0.0`.*
