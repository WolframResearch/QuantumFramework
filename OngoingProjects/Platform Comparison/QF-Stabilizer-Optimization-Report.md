# QF stabilizer optimization, round 2: compiled fold + full hot-path rework

**Date:** 2026-06-12 (revision; original 2026-06-11. Round 1 = `QF-Stabilizer-Packed-FastPath.md`, plan = `QF-Stabilizer-Bottleneck-Audit.md`). **Still current at HEAD `f9dc1cdc` (2026-06-14):** the mid-June property-cache refactor (`dfc741dc`/`f9dc1cdc`) touches the `QuantumOperator`/`QuantumState`/`QuantumBasis` property dispatch and the state-vector apply path, not `Stabilizer/`; the numbers below reproduce within run-to-run noise (e.g. `ApplyCircuit` $n=1000$ 62.7 ms, ctor 5.08 ms).
**Repo:** WolframResearch/QuantumFramework (paclet 2.0.0). Round 1 shipped as `95566e93`, round 2 as `e37b7fad`; this revision re-measures everything at HEAD `6f953ad5`, which includes three post-ship stabilizer fixes (`1c524432`, `52d3ea89`, `6f953ad5`; see "Kernel updates since the round-2 ship").
**Machine:** Apple M2 Pro, Mathematica 15.0.0. Externals re-measured fresh at this revision: Stim 1.16.0, QuantumClifford.jl 0.11 (Julia 1.12.6).
**Scope:** `QuantumFramework/Kernel/Stabilizer/` only; nothing outside the stabilizer subsystem was touched.

## The physics, in one paragraph

A stabilizer state is specified by the Pauli group that fixes it, stored as a binary
symplectic tableau ($2n$ generator rows over $n$ qubits, plus a $\pm 1$ sign per row), so a
Clifford gate is a constant-size rewrite of one or two qubit columns. Round 1 packed the
tableau bits into machine words, fixing the $O(n^2)$-per-gate scaling. This round removes
the remaining cost: the Wolfram interpreter's fixed $\sim 1\,\mu\text{s}$ tax per primitive
operation, paid a dozen times per gate plus object wrap/unwrap. The cure is a single
compiled (C) routine that takes the whole encoded circuit and the packed tableau and runs
the entire gate loop natively, plus surgical fixes to every other hot path that was
re-traversing the tableau element-by-element (measurement, formatting, construction).

## Headline result: full Clifford-circuit simulation, identical op stream

Deterministic H/S/CNOT stream (the cross-platform benchmark from the master comparison);
total wall-clock ms, min-of-3, `ClearSystemCache[]` each rep, constructor included. All
columns re-measured 2026-06-12 at HEAD `6f953ad5` except "QF original" (the pre-rework
implementation, from the master comparison):

| $n$ (gates) | QF original | QF round 1 (per-gate packed) | **QF round 2 (`ApplyCircuit`)** | Stim | QuantumClifford.jl |
|---:|---:|---:|---:|---:|---:|
| 100 (2 000) | 681 | 69.2 | **2.56** | 0.047 | 1.08 |
| 500 (10 000) | 15 836 | 412 | **26.2** | 0.330 | 23.1 |
| 1000 (20 000) | 63 469 | 951 | **58.0** | 0.953 | 108.4 |

> **Correction (2026-07-20).** Stim column restated: the harness had been issuing one Python call per gate, which times the pybind11 boundary rather than the tableau engine (per-gate cost flat at ~0.5 us across all $n$, where an $O(n)$ update cannot be flat). Batched through one `sim.do(circuit)`: 23.6/33.0/47.6 ns per gate. QF and QC.jl columns are unaffected and reproduce.

- **Cumulative speedup vs the original implementation: $266\times$ / $606\times$ / $1095\times$.**
- **Gap to Stim: quote phases, not one number** (2026-07-20 re-measurement, `scripts/bench_phases_*`, after two opposite-direction single-ratio errors). Compiled kernel vs Stim engine call: $7.3$/$5.4$/$4.4\times$ at $n=100$/$500$/$1000$ (narrowing with $n$; both builds scalar). End-to-end vs Stim fast-I/O pipeline: $11$-$16\times$, dominated by QF's encode/pack/materialize boundaries. The pre-2026-07-20 "$2.5$/$5.1$/$5.3\times$" mixed boundaries over a Stim harness timing Python call overhead.
- **QF beats QuantumClifford.jl at $n = 1000$** (58.0 vs 108.4 ms, $1.9\times$ faster) and roughly matches it at $n = 500$.
- Per-gate cost: 1.3-2.9 $\mu\text{s}$/gate (was $\sim 340$-$3170\,\mu\text{s}$ originally).
- The round-2 numbers are stable across the post-ship fixes: the 2026-06-11 measurements (2.63 / 26.3 / 58.6 ms) reproduce within noise at HEAD.

The compiled kernel alone (state already packed, circuit already encoded) runs at
0.22-0.27 $\mu\text{s}$/gate. (Earlier editions called this "Stim's own regime"; it is not, Stim runs the same stream end to end at 0.024-0.048 $\mu\text{s}$/gate.) The remaining wall-clock is
circuit encoding ($\sim 0.7\,\mu\text{s}$/gate) and the pack/unpack boundary (one-time per
circuit).

## What was implemented (file by file)

### 1. Compiled bulk gate fold: `Stabilizer/Compiled.m` (new file)

The core of round 2. One `Compile[..., CompilationTarget -> "C"]` kernel
(`compiledCliffordFold`) executes an entire encoded Clifford circuit on a
**generator-major** packed tableau: for each qubit, its X (resp. Z) bits across all $2n$
generators are packed into $\lceil 2n/62 \rceil$ machine words, and the $2n$ sign bits are
one word vector. A single-qubit gate then reads/writes the few words of its qubit's column,
and the sign update is one vectorized XOR, the same data layout Stim and
QuantumClifford.jl use; that layout is what makes Clifford updates word-parallel.

- Gate set: H, S, S†, CNOT/CX, SWAP, X, Y, Z natively; CZ expands to H·CNOT·H. The
  integer encoding lives in **one** association (`$stabilizerGateCodes`) injected into both
  the compiled body and the encoder via `With`; no duplicated literals.
- Circuit encoder (`encodeStabilizerGates`): a single `Dispatch`-table `Replace` pass over
  the spec list ($\sim 0.7\,\mu\text{s}$/spec; the first, layered version cost
  $4.5\,\mu\text{s}$/spec). Accepts `op -> q`, `op -> {q...}`, the `QuantumShortcut`
  controlled form `"C"[g -> t, c, _]`, and the legacy `{op, q...}` form. Any unrecognized
  gate makes the whole encode return `$Failed`, which routes the caller to the per-gate
  fold, so non-Clifford circuits (`T` and `P[θ]`, which return a `StabilizerFrame`)
  behave exactly as before.
- Compilation target auto-selects: `"C"` when a C compiler is present
  (`CCompilerDriver`CCompilers[]`), `"WVM"` bytecode otherwise (still removes the
  per-gate interpreter dispatch). The kernel is built lazily once and memoized.
- Entry points: `ps["ApplyCircuit", specs]` (public method) and `PauliStabilizerApply`
  ([PauliStabilizer.m:123](../../QuantumFramework/Kernel/Stabilizer/PauliStabilizer.m)), so
  `qc[Method -> "Stabilizer"]` circuits take the compiled path automatically when every
  gate is encodable and the state is concrete.

**Two packed-array traps cost a day and are worth recording.** (a) A compiled function
silently falls back to *interpreted* evaluation when any argument is an unpacked array,
$\sim 100\times$ slower with zero warning; every boundary now forces
`Developer`ToPackedArray`. (b) The word values reach $2^{61}$, and several construction
routes (notably `Dot` results fed through other ops) come back unpacked; this is invisible
unless you check `Developer`PackedArrayQ` at each step.

### 2. Closed-form integer constructor: `Constructors.m:259`

`PauliStabilizer[n]` (the $|0\rangle^{\otimes n}$ register) was routing a known-in-closed-form
tableau through the generic validating constructor: `ArrayQ[..., 3, MatchQ[0|1|-1]]` over
a $\{2, n, 2n\}$ array, i.e. millions of per-element pattern matches. It now assembles
`Signs`/`Tableau` directly (destabilizers $X_i$, stabilizers $Z_i$, all signs $+1$).

| | $n=100$ | $n=500$ | $n=1000$ |
|---|---:|---:|---:|
| before | 4.0 ms | 89.5 ms | 355 ms |
| after (at HEAD) | **0.078 ms** | **1.17 ms** | **5.1 ms** |

($\sim 50$-$77\times$; verified `===`-identical to the old route for $q = 1..12$, test
`Tier10-IntegerCtor-ClosedForm`.)

### 3. Packed AG measurement: `Measurement.m`

`ps["M", a]` was reconstructing the full rank-3 array from the packed state on every call
and folding interpreted `rowsum`s. It now measures **directly on the packed
representation**: the anticommuting-row scan is one `BitAnd`+`Sign` over the qubit's chunk
column, and the AG phase of row products is the carry-save two-parity-bit accumulation
(`packedRowSumPhase`), the audit-§14.1 trap, validated against the canonical `agPhase`
sum on 3000 random row pairs and against the canonical measurement Association on 500
random states before integration.

| single Z-measurement | $n=100$ | $n=500$ |
|---|---:|---:|
| before | 6.85 ms | 16.5 ms |
| after (at HEAD) | **0.101 ms** | **0.228 ms** |

($\sim 68\times$.) Deterministic/random branching, post-states, and correlations (Bell
`{00,11}`, GHZ `{000,111}`) are exactly preserved; an $n$-qubit measurement sweep dropped
from $O(n^3)$ to $O(n^2/62)$.

### 4. Pauli-string measurement deterministic branch: `PauliMeasure.m:55`

The deterministic outcome was computed by materializing the $2^n$ state vector and
sandwiching the $2^n \times 2^n$ Pauli matrix. It now calls the existing AG closed-form
`stabilizerExpectation` ($O(n^2)$, i-factor tracking), so the exponential wall is gone
entirely on the deterministic branch. (The *non-deterministic* string branch was
subsequently rewritten for correctness in `52d3ea89`, which made it the one slow path left
in the subsystem, and then moved onto the packed representation in `730abcec`; see "Round
3" below.)

### 5. Vectorized Pauli-string formatting: `Formatting.m`

`ps["Stabilizers"]` (and every summary box) was running an element-by-element pattern
`Replace` over the transposed tableau. The Pauli letter is now read by arithmetic indexing
$x + 2z$ into a letter table, one vectorized `Part` per row: at $n=100$, from 128.6 ms to
**2.03 ms** ($63\times$).

### 6. Packed-aware dimension handlers: `Properties.m:9`

`Qubits`/`Qudits`/`GeneratorCount` on a packed object read the stored count instead of
reconstructing the tableau to take its `Dimensions`; this was a hidden $O(n^2)$ rebuild
inside every measurement and formatting call.

### 7. Composition crash fix: `Compose.m:35`

The phase loop called `Rest`/`Most` on the picked-generator list, which threw on any row
selecting 0 or 1 generators (all-zero Pauli rows; mixed-size compositions hit it
reliably). Rows selecting $\leq 1$ generators contribute no consecutive-product phase, so
they now short-circuit to 0. Verified: 200 random mixed-size compositions, zero messages.
Incidentally, composition also got $\sim 55\times$ faster (at $n=100$, from 2577 ms to **46.1 ms**)
because the packed states now reconstruct to *packed* arrays, which the phase loop's
`Pick`/`FoldList` consume far more efficiently; the algorithm is still the $O(n^3)$
symplectic multiply, unchanged and correct.

## Kernel updates since the round-2 ship (2026-06-11 evening through 2026-06-12)

Four commits touched the stabilizer subsystem after `e37b7fad`; all are in the numbers
above.

### `1c524432`: AG measurement destabilizer rowsum + register sizing (correctness)

A final end-to-end check, sequential full-register measurement of random Clifford circuits
compared against dense simulation, failed on 3 of 15 circuits. Arbitration with Stim showed
the **dense result was right and the stabilizer measurement wrong**: a *pre-existing* bug
(present before any of this work; verified on the pre-rework tree) in the non-deterministic
Z-measurement branch, which rowsummed only the anticommuting rows *above* the pivot,
omitting the destabilizers. AarGot04 §3 requires $\mathrm{rowsum}(i, p)$ for *every*
$i \neq p$ with $x_{ia} = 1$; the stale destabilizers break the symplectic pairing and
corrupt **later deterministic outcomes** on a subset of branches, invisible to
single-measurement tests, to the gate-level cross-package suites, and to the state-vector
materialization (which reads only the stabilizer half). The packed measurement had
faithfully reproduced the canonical algorithm, bug included. Fixed in both branches. A
second pre-existing hole (`PauliStabilizerApply` sizing the register by wire *count*
instead of highest wire *index*, garbling circuits like `{"H" -> 2}`) was fixed in the same
commit. 15 permanent TIER 11 regression tests pin sequential measurement to the
materialized state-vector support, including the original failing seed.

The lesson for the test suite: the original tests checked single measurements, gate
conjugations, and textbook states, but never compared *multi-qubit measurement enumeration*
against an independent oracle. That class of check now exists.

### `52d3ea89`: Pauli-string measurement correctness + `RandomClifford` viability

Two independent fixes:

- **`PauliMeasure.m`**: the non-deterministic Pauli-*string* branch had the same disease
  `1c524432` cured for single-qubit Z (it overwrote the pivot row without rowsumming the
  other anticommuting rows, destabilizers included, leaving a non-abelian "stabilizer"
  group and no destabilizer update at all). Rewritten to mirror the corrected algorithm:
  rowsum every other anticommuting row into the pivot, promote the pivot's old row to its
  destabilizer slot, install $(-1)^b P$ as the new stabilizer. Cross-validated against
  dense simulation and Stim; regression tests `Tier11-PauliStringMeasure-*`.
- **`RandomClifford.m`**: two latent blow-ups that made `PauliStabilizer["Random"[n]]`
  die around $n \approx 500$. (a) Structured-array wrappers
  (`SymmetrizedArray`/`LowerTriangularMatrix`) pushed the downstream `Dot` off the packed
  matmul path onto a generic tensor route whose intermediates reach tens of GB; plain
  packed integer matrices are now used throughout. (b) `Inverse` of a unit lower-triangular
  0/1 matrix was taken over the integers, where entries grow exponentially with $n$, turning
  the final $2n \times 2n$ `Dot` into a big-integer matmul; the tableau lives over
  $\mathbb{F}_2$, so it is now `Inverse[..., Modulus -> 2]`.

### `aa45cd0e`: test-suite alignment (no kernel behavior change in the stabilizer subsystem)

Stabilizer tests updated to the current fallback semantics (non-Pauli-basis measurement
falls back to the dense `QuantumMeasurement` path with `PauliStabilizer::nonpaulibasis`;
stochastic Pauli channels return `$Failed` with `CliffordChannel::stochastic`;
non-frame-decomposable gates like `RX[0.3]` emit one `PauliStabilizer::nonclifford` and
return `$Failed`, while `T` returns a `StabilizerFrame` silently). Plus two cosmetic
pattern-variable renames outside the subsystem (`Labels.m`, `QuantumState/Properties.m`).

### `6f953ad5`: named-circuit call-form shortcut dispatch (constructor surface)

`PauliStabilizer["GHZ"[5]]` fell through unevaluated: the circuit-shortcut constructor rule
matched bare strings, rules, and lists, but never the v2.0 `"Name"[args]` call form. The
pattern is widened, bare shorthand rules are list-wrapped (the circuit constructor rejects
them bare), and an unrecognized name propagates a `Failure` via `Enclose`/`ConfirmBy`
instead of returning an unevaluated expression. Separately, bare
`PauliStabilizer["Random"]` was being captured by the shortcut (returning a corrupt
object) because WL does not reorder the composite `Optional`-argument rule from a
later-loading file; literal rules for `"Random"` and `"Random"[n]` were added and
registered `$PauliStabilizerNames` are excluded from the shortcut, so dispatch is
load-order independent. Nine TIER-1 regression tests in `AuditMatrix.wlt`.

## Round 3 (`730abcec`, 2026-06-12): packed non-deterministic Pauli-string measurement

The `52d3ea89` correctness rewrite of `ps["M", "ZZ...Z"]` (non-deterministic branch) traded
speed for exactness: essentially its entire cost was the `Fold` of the interpreted `rowsum`
over the anticommuting rows (97 of them at $n=100$), each step copying the full rank-3
$\{2, n, 2n\}$ tableau, $O(n^3)$ data motion with interpreter constants. (The symplectic
product `Mod[ps["Matrix"] . omega . pVec, 2]` and the `Signs`/`Tableau` reads were
immaterial, $\sim 2$ ms combined at $n=100$.) `730abcec` gives the branch the same packed
treatment the single-qubit measurement got in round 2 (`packedMeasurePauliString`,
`PauliMeasure.m`):

- The anticommuting-row scan encodes the measured Pauli's symplectic vector into packed
  X/Z word-vectors (one machine word per 62-qubit chunk) and takes, per chunk, the
  popcount parity of `BitAnd[x_row, z_P]` and `BitAnd[z_row, x_P]` with a vectorized
  `DigitCount` over all $2n$ rows, replacing the dense matrix product.
- The AG row-clearing loop reuses `packedRowSumPhase` (the carry-save two-parity-bit
  i-power accumulation), so each cleared row costs $O(n/62)$ machine-word XORs instead of
  a rank-3 tableau copy.
- Pivot promotion and the $(-1)^b P$ install are per-chunk column writes. Branch states
  come back *packed*, so chained measurements stay on the fast path (a second string
  measurement on a branch state at $n=100$: 9.1 ms end to end).
- Symbolic-sign states (non-Clifford P/T residues, SymPhase) fail the `psConcreteFastQ`
  gate and keep the canonical rank-3 path, unchanged.

Measured on the benchmark random-Clifford state (`SeedRandom[3]`, all-Z string, min-of-3,
`ClearSystemCache` per rep, M2 Pro):

| `ps["M", "Z..Z"]`, non-deterministic | $n=100$ | $n=500$ | $n=1000$ |
|---|---:|---:|---:|
| round 2 (`e37b7fad`, **incorrect** branch) | 3.4 ms | (not measured) | |
| `6f953ad5` (correct, interpreted rank-3 `Fold`) | 393 ms | 42.3 s | (not measured) |
| `730abcec` (correct, packed) | **4.9 ms** | **58.7 ms** | **211 ms** |

That is $80\times$ at $n=100$ and $720\times$ at $n=500$, with the deterministic branch
(`stabilizerExpectation`, $O(n^2)$) untouched.

### Destabilizer-phase fix in the round-2 `packedMeasureZ` (same commit)

The exact packed-vs-canonical equivalence gate caught a latent bug in the *shipped*
round-2 single-qubit path as well: when a cleared row is a **destabilizer**, the AG
i-power of the row product can be odd ($\pm i$; destabilizer rows carry no Hermitian sign
constraint), and the packed phase update `Mod[exponent, 4]/2` then produced a fractional
phase bit, i.e. a sign of $0$, where the canonical `rowsum` collapses an odd exponent to
sign $+1$ (`Boole[exponent == 2]`). Destabilizer signs are not physical, so no outcome
statistics or post-states were ever wrong, but the $0$ broke exact equality with the
canonical path and, worse, de-concretized the branch's `"Signs"`, silently knocking every
subsequent operation on that branch off the packed fast path (22 of 200 random single-qubit
measurements at $n=5$ were affected). Both the new string path and `packedMeasureZ` now use
the canonical `Boole` convention; the deterministic scratch-row accumulation (stabilizer
products only, provably even exponent) was aligned for uniformity.

## Verification

Everything is gated on exact equivalence, not statistics.

- **Stabilizer suite at `730abcec`: 965 / 965, zero failures** (re-run 2026-06-12,
  fresh kernel; was 945/945 at `6f953ad5`, plus the 20 new round-3 pins below). The 6
  message-expectation failures that trailed the suite through round 2 were resolved by the
  `aa45cd0e` test alignment. The suite has also grown and been reorganized since the
  round-2 report (was 921 tests): `AuditMatrix.wlt` 203, `CliffordChannel.wlt` 76,
  `Connections.wlt` 52 (new), `Correctness_TextbookResults.wlt` 42,
  `CrossPackage_QuantumClifford.wlt` 15, `CrossPackage_Stim.wlt` 44, `Formulas.wlt` 132
  (new), `HybridInterop.wlt` 52, `PauliStabilizer.wlt` 329, `Roundtrips.wlt` 20.
- **Round-3 pins** (`PauliStabilizer.wlt`, TIER 11): 9
  `Tier11-PauliStringMeasure-PackedCanonical-*` and 9 `Tier11-Measure-PackedCanonical-*`
  tests compare the packed paths against the canonical rank-3 path (forced by Blocking the
  `psConcreteFastQ` gate) for **exact** outcome-key / `"Signs"` / `"Tableau"` equality plus
  sign concreteness, at $n \in \{5, 20, 40\}$ with signed random strings; 2
  `Tier11-PauliStringMeasure-Packed-Perf-*` tests pin the budget ($n=100$ under 50 ms,
  $n=500$ under 2 s, $\sim 10\times$ headroom over the measured 4.9 / 58.7 ms).
- **Round-3 cross-validation** (throwaway harnesses, 2026-06-12): packed-vs-canonical
  exact equality on 435 random cases up to $n=70$; branch post-states vs the dense
  projection $(1 + (1-2b)P)/2\,|\psi\rangle$ up to global phase and deterministic outcomes
  vs the dense $\langle P \rangle$ on 300 random cases at $n \le 5$; deterministic/random
  classification and deterministic outcome sign vs Stim 1.16
  (`peek_observable_expectation`) on 120 random cases up to $n=20$, all clean.
- **Full framework suite at HEAD: 1493 / 1493** (re-run 2026-06-12, `Tests/RunTests.wls`,
  fresh kernel, 100%).
- **Cross-package suites: `CrossPackage_Stim.wlt` 44/44 and
  `CrossPackage_QuantumClifford.wlt` 15/15.** The packed/compiled paths reproduce both
  external simulators' documented gate-conjugation rules exactly.
- **New test tiers** in `Tests/Stabilizer/PauliStabilizer.wlt` (all passing):
  - TIER 9 (8 tests): interpreted packed gates vs a self-contained canonical
    transcription, 400 random Cliffords each, $n \in \{3,...,130\}$ spanning the 62/63
    chunk boundary.
  - TIER 10 (10 tests): compiled `ApplyCircuit` vs the per-gate fold (6 sizes, all 8
    gates + CZ), empty-circuit identity, CZ symmetry, single-gate equivalence per gate
    type, non-Clifford fallback to `StabilizerFrame`, and the constructor closed-form
    equivalence.
  - TIER 11 (12 + 15 tests): packed measurement against an *independent* projective
    measurement of the materialized state vector ($P_b = (1 \pm Z_a)/2$ on $|\psi\rangle$,
    matched up to global phase) at $n = 2..5$, determinism/correlation structure (Bell,
    GHZ, stabilizer vs anticommuting Pauli strings), the out-of-range message, plus the
    `1c524432` sequential-measurement regression pins and the `Tier11-PauliStringMeasure-*`
    pins from `52d3ea89`.
- **Pre-integration validation** (throwaway harnesses, then promoted into the tiers):
  carry-save phase vs canonical on 3000 random row pairs; packed measurement vs canonical
  Association on 500 random states; compiled fold vs interpreted fold bit-for-bit at 7
  sizes; encoder output for every accepted spec shape.

## Caveats

- `PauliStabilizer["Random"[1000]]` emits `General::munfl` (machine underflow) during
  construction, a pre-existing Mallows-sampler issue independent of the measurement paths
  ($n \le 500$ is clean; the constructed state measures correctly).
- The `"WVM"` compilation-target branch (no C compiler installed) is selected by a guard
  that is straightforward, but this machine has clang, so the WVM path itself has not been
  exercised end-to-end here.
- `ApplyCircuit` returns a **canonical** (rank-3) object, unlike single gates which return
  packed ones; both forms are accepted everywhere (`PauliStabilizerQ` handles both), but
  `SameQ` across storage forms remains the round-1 caveat: compare `["Tableau"]` /
  `["Stabilizers"]`, not raw objects.
- Composition and `["Circuit"]` synthesis remain algorithmically $O(n^3)$ (correct,
  crash-fixed, and now fast enough at the benchmark sizes); a packed rewrite of the
  composition *phase convention* was prototyped but could not be verified against the
  existing `phaseLookup` semantics, so it was **not** shipped; the verified crash fix was.
- `QuantumCircuitOperator` *construction* for $10^4$-gate circuits is slow in the host
  framework (minutes); that is outside the stabilizer subsystem and out of scope. Use
  `ps["ApplyCircuit", specs]` or `PauliStabilizerApply` for big gate streams.

## Reproduce

```
# Op streams to /tmp first:  cp scripts/stab_ops_*.json /tmp/
# QF (repo working tree):    wolframscript -file scripts/qf_final_all.wl
# Stim (venv):               /tmp/qcbench/bin/python scripts/bench_stim_full.py
# QuantumClifford.jl:        JULIA_DEPOT_PATH=/tmp/qc_depot julia scripts/bench_qc_full.jl
# Stabilizer test suite:     wolframscript -file scripts/run_stab_tests.wl
# Full framework suite:      wolframscript -file Tests/RunTests.wls   (repo root)
```

Op streams and all scripts are persisted under `OngoingProjects/Platform Comparison/scripts/`.
