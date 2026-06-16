# Bit-packed fast path for `PauliStabilizer` Clifford gate updates

**Date:** 2026-06-11
**Repo:** WolframResearch/QuantumFramework (paclet 2.0.0), branch worktree off `main` @ `77b9ccba`
**Machine:** Apple M2 Pro, Mathematica 15.0.0
**Scope:** `QuantumFramework/Kernel/Stabilizer/` Clifford gate updates (the hot path)

> **This is the round-1 snapshot (historical).** Rounds 2 and 3 (the compiled bulk fold and packed measurement) superseded these per-gate numbers; see `QF-Stabilizer-Optimization-Report.md` for the shipped end state ($266$ to $1095\times$ over the original) and `QF-Stabilizer-Bottleneck-Audit.md` for the plan-vs-outcome board. The 2026-06-13/14 property-cache refactor (`dfc741dc`/`f9dc1cdc`) does not touch the stabilizer subsystem, so nothing here changed at the current HEAD `f9dc1cdc`.

## The physics, in one paragraph

A stabilizer state on $n$ qubits is specified not by its $2^n$ amplitudes but by the
group of Pauli operators that fix it. The Aaronson-Gottesman (AG) tableau stores a
generating set as a binary symplectic array: for each of $2n$ rows (the $n$ stabilizer
generators plus $n$ destabilizers) it records, per qubit, whether an $X$ and/or $Z$ appears,
together with a $\pm 1$ sign. A Clifford gate conjugates the Paulis, so simulating it is
just a rewrite of this array, with no exponential state vector. The cost of a single gate
should therefore be $O(n)$: a Clifford generator touches only its one or two qubit
columns. The previous implementation instead paid $O(n^2)$ per gate because it rebuilt the
entire rank-3 array (`ReplacePart` over all $2\cdot n \cdot 2n$ entries) and re-scanned the
length-$2n$ sign list element by element through pattern-matched property dispatch. Packing
each row's $X$ and $Z$ bits into 62-bit machine words and operating on whole columns at once
restores the $O(n)$ scaling and removes the per-element interpreter overhead.

## What changed and where

The canonical rank-3 integer array (`ps["Tableau"]`, dimensions $\{2, n, 2n\}$) and the
length-$2n$ `ps["Signs"]` list remain the display/serialization representation and the
fallback for symbolic-phase work. A second, *packed* representation is introduced and the
hot gate kernels route through it.

### New file: `QuantumFramework/Kernel/Stabilizer/Packed.m`

Packed storage is **chunk-major**: a packed `PauliStabilizer` carries the keys

```
"PackedX" -> { chunk_1, chunk_2, ... },   (* each chunk is a length-2n packed Int64 vector *)
"PackedZ" -> { ... },
"Signs"   -> { (-1|1) ... },              (* length 2n *)
"Qubits"  -> n
```

and **no** `"Tableau"` key. Bit $b$ of chunk $c$ of row $r$ holds qubit $(c\cdot 62 + b)$ of
that tableau row, so each chunk packs up to 62 qubits and stays inside a machine integer
($2^{62} < 2^{63}-1$). Because qubits within a chunk share one Int64 column vector, a gate
on qubit $j$ touches only the single chunk column $\lceil j/62 \rceil$ (CNOT/SWAP touch the
two columns of $j$ and $k$). Each gate is then a handful of `BitAnd` / `BitXor` / `Sign`
operations on length-$2n$ **packed arrays of machine integers**, which the kernel executes
as vectorized Listable ops with no arbitrary-precision overhead.

Key symbols (`Packed.m`):
- `packChunks` / `unpackChunks` / `tableauFromPacked` ([Packed.m:71](../../QuantumFramework/Kernel/Stabilizer/Packed.m), :76, :81): convert between the rank-3 array and chunk-major form.
- `psPackedQ` ([Packed.m:79](../../QuantumFramework/Kernel/Stabilizer/Packed.m)), `psConcreteFastQ` ([Packed.m:86](../../QuantumFramework/Kernel/Stabilizer/Packed.m)): the strict numeric predicate that gates the fast path (packed already, or concrete binary tableau with concrete $\pm 1$ signs). Symbolic-phase / `StabilizerFrame` states fail it and take the canonical path.
- `psGetPacked` ([Packed.m:101](../../QuantumFramework/Kernel/Stabilizer/Packed.m)) / `psFromPacked` ([Packed.m:110](../../QuantumFramework/Kernel/Stabilizer/Packed.m)): tuple `{n, packedX, packedZ, signs}` ⇄ object.
- `ps["Tableau"]` lazy reconstruction handler ([Packed.m:116](../../QuantumFramework/Kernel/Stabilizer/Packed.m)): so every downstream property (`State`, `Measure`, `Stabilizers`, formatting, ...) works unchanged on a packed object.
- `packedGateH/S/Sdg/CNOT/SWAP/X/Y/Z` ([Packed.m:134](../../QuantumFramework/Kernel/Stabilizer/Packed.m)-:220): the bitwise kernels, each mirroring its canonical AG rule bit-for-bit.

### `QuantumFramework/Kernel/Stabilizer/GateUpdates.m`

Each Clifford generator rule now dispatches through the packed kernel when
`psConcreteFastQ[ps]` holds, with the original rank-3-array body kept verbatim as the
`else` branch:
H ([GateUpdates.m:30](../../QuantumFramework/Kernel/Stabilizer/GateUpdates.m)),
S (:41), S† (:52), CNOT/CX (:63), X (:77), Y (:83), Z (:89), SWAP (:98).
`["Dagger"]/["Inverse"]` ([GateUpdates.m:164](../../QuantumFramework/Kernel/Stabilizer/GateUpdates.m))
now **preserves the input's storage form** (packed in → packed out) so structural
(`===`) comparisons stay consistent. `CZ` composes from packed `H`+`CNOT` automatically;
`V`/`V†` go through composition (see "Not packed" below).

### `QuantumFramework/Kernel/Stabilizer/PauliStabilizer.m`

`PauliStabilizerQ` and `ConcretePauliStabilizerQ` were rewritten as a single `Which`-body
that accepts both the canonical and packed key sets
([PauliStabilizer.m:32](../../QuantumFramework/Kernel/Stabilizer/PauliStabilizer.m), :49). This
also removed a latent fragility where the old `Tableau`-form pattern condition
(`Dimensions[tableau][[3]]`) misfired when probed against a keyless association.

### Tests

- `Tests/Stabilizer/PauliStabilizer.wlt`: the test-internal `psValidQ` helper now validates
  the *reconstructed* binary tableau (accepting both forms); a new **TIER 9** section adds
  8 packed-vs-canonical equivalence tests that fold 400 random Cliffords (H, S, S†, X, Y, Z,
  SWAP, CNOT) and assert `kernel["Tableau"] === canonical-reference` and same `["Signs"]`,
  at $n \in \{3,8,40,62,63,70,100,130\}$ (single-chunk, chunk-boundary, multi-chunk). The
  reference functions are a self-contained transcription of the pre-packing rank-3 rules.

## The packed gate algebra (for the record)

Writing $m = 2^{(j-1) \bmod 62}$ for the qubit-$j$ mask inside its chunk, and `Sign` of a
masked column as the per-row 0/1 indicator:

| Gate | $X$ update | $Z$ update | sign flip where |
|---|---|---|---|
| $H_j$ | swap bit $j$ with $Z$ | swap bit $j$ with $X$ | $x_j = z_j = 1$ |
| $S_j$ | - | $z_j \mathrel{\oplus}= x_j$ | $x_j = z_j = 1$ |
| $S^\dagger_j$ | - | $z_j \mathrel{\oplus}= x_j$ | $x_j = 1,\ z_j = 0$ |
| $\mathrm{CNOT}_{j\to k}$ | $x_k \mathrel{\oplus}= x_j$ | $z_j \mathrel{\oplus}= z_k$ | $x_j = z_k = 1 \ \wedge\ x_k = z_j$ |
| $\mathrm{SWAP}_{jk}$ | exchange bits $j,k$ | exchange bits $j,k$ | (none) |
| $X_j,Z_j,Y_j$ | - | - | $z_j{=}1$ / $x_j{=}1$ / $x_j \neq z_j$ |

The sign update is a single Listable multiply `s (1 - 2*flip)` over the length-$2n$ sign
vector. This is the AG rule set; the two-parity-bit carry-save phase of the general Pauli
multiply (audit §4) is not needed here because single Clifford generators only ever flip a
sign, never accrue an $i$.

## Benchmark

Identical deterministic op stream across all three simulators (0-indexed; `+1` for WL):
`i mod 3 == 0 → H(i mod n)`; `== 1 → S((7i+1) mod n)`; `== 2 → CNOT(a,b)` with
`a=i mod n, b=(i+(i mod 5)+1) mod n` (bumped if `b==a`). Same JSON op files fed to QF, Stim,
and QuantumClifford.jl. QF timed by `AbsoluteTiming`, min-of-3, `ClearSystemCache[]` each run
(`SetAttributes[tm, HoldFirst]`); Stim `perf_counter` min-of-5; QC.jl `@belapsed`.

### Total wall-clock (ms)

| $n$ | gates | **QF before** | **QF after** | **speedup** | Stim 1.16.0 | QuantumClifford.jl 0.11 |
|---:|---:|---:|---:|---:|---:|---:|
| 100  | 2 000  | 681.2   | **74.6**  | **9.1×**  | 1.04  | 1.18   |
| 500  | 10 000 | 15 836  | **489.1** | **32.4×** | 5.51  | 23.76  |
| 1000 | 20 000 | 63 469  | **1 254.8** | **50.6×** | 11.33 | 110.78 |

### Per-gate (ms/gate)

| $n$ | QF before | QF after | Stim | QC.jl |
|---:|---:|---:|---:|---:|
| 100  | 0.341 | 0.037 | 0.00052 | 0.00059 |
| 500  | 1.584 | 0.049 | 0.00055 | 0.00238 |
| 1000 | 3.173 | 0.063 | 0.00057 | 0.00554 |

### Reading the numbers

- **Scaling fixed.** Before, per-gate cost grew linearly in $n$ (0.34 → 1.58 → 3.17 ms/gate),
  the signature of $O(n^2)$ work per gate. After, it grows only mildly (0.037 → 0.049 → 0.063
  ms/gate): the residual growth is the length-$2n$ sign-vector multiply and the
  $\lceil n/62 \rceil$-element chunk-list `ReplacePart`, both small. Hence the speedup *grows*
  with system size: 9× at $n=100$, 32× at $n=500$, 51× at $n=1000$.
- **Gap to Stim.** The factor over Stim collapsed from roughly $650$-$5600\times$ (before) to
  $72$-$111\times$ (after). Because both now scale near-linearly per gate, this remaining gap is
  an almost-constant factor, not a widening asymptotic one.
- **Versus QuantumClifford.jl.** QC.jl uses the same chunked-UInt64 idea but in compiled Julia
  with SIMD. After packing, QF closes from $\sim\!570\times$ (before, $n{=}1000$) to $\sim\!11\times$
  QC.jl at $n=1000$.

## Correctness

- All `Tests/Stabilizer/*.wlt` run with the working tree loaded: **878 / 884 pass**. The 6
  failures are **pre-existing on the unmodified baseline** (verified by stashing the change):
  they are message-expectation tests unrelated to stabilizer arithmetic
  (`Phase2-NonClifford-Message`, `Phase2-PauliStabilizer-UsageMessage`,
  `Audit-Hybrid-QMO-NonPauli-Fallback`, `Phase8.2-CCfromQC-BitFlip-StochasticNotice`,
  `Phase7-QMO-NonPauliBasis-Fallback-EmitsMessage`, `Conn5-qmoComputational-FallsBack`).
  **Zero regressions.**
- The two external cross-validation suites are the anchor: `CrossPackage_Stim.wlt` (44/44)
  and `CrossPackage_QuantumClifford.wlt` (15/15) both pass: the packed path reproduces the
  documented Stim and QuantumClifford.jl gate-conjugation rules exactly.
- New TIER 9 equivalence tests (8/8) assert packed `===` canonical across single-chunk,
  chunk-boundary ($n=62,63$), and multi-chunk sizes under `SeedRandom`.

### Caveats

- **Storage form is observable via `===`.** A concrete `PauliStabilizer` that has passed
  through a gate is stored packed (different keys from a freshly constructed canonical one),
  so `gateChain === bareConstructor` can be `False` even when the states are equal. All
  in-repo equality tests compare expressions that share a prefix (hence the same form), and
  `Dagger` was made form-preserving to keep the involution tests green, but external callers
  who compare a gate result to a hand-built canonical association with `SameQ` should compare
  `ps["Tableau"]`/`ps["Stabilizers"]` instead.
- **62-bit chunks**, not 64: an Int64 sign bit is reserved, so the mask $2^{61}$ and sums like
  $2^{61}+2^{61}=2^{62}$ stay representable. No effect on results; a negligible ~3% more chunks
  than a full-64 packing.

## Not packed (canonical fallback, correct but not yet accelerated)

These route through the packed object via lazy `ps["Tableau"]` reconstruction and then run
the existing canonical code. They are correct and fully test-green (the Stim/QC suites
exercise measurement), but are not yet bit-packed:

- **Composition** `left[right]` (`Compose.m`): still the dense $O(n^3)$ symplectic
  matrix-mod-2 multiply (measured 179 / 2578 / 10085 ms at $n=50/100/200$). Profiling shows
  the dominant cost is the per-row phase-correction loop (`phaseLookup` + `Pick`/`FoldList`),
  not the matmul. The packed plan (audit §4.2): pack the $2n\times 2n$ symplectic matrix
  row-major and accumulate the product by XOR of selected rows, tracking phase with the
  two-parity-bit carry-save scheme. `V`/`V†` ride on this. *(A latent edge bug was noted in
  `Compose.m`'s phase code: an all-zero Pauli row makes `Rest[Pick[...]]` take `Rest[{}]`;
  worth a separate fix.)*
- **Measurement** (`Measurement.m` single-qubit Z, `PauliMeasure.m` Pauli-string): run AG
  `rowsum` on the reconstructed rank-3 array. The packed kernel is the audit §4 / §7 Pauli
  multiply with destabilizer-trick projection; deferred because `rowsum`'s two-parity-bit
  phase is the single most correctness-sensitive primitive in the subsystem and is not on
  the gate-simulation hot path the benchmark targets.

## Reproduce

```
# QF (before = stash the kernel change; after = working tree):
wolframscript -file /tmp/qf_bench_full.wl          # n = 100, 500, 1000

# Stim (venv with `pip install stim`):
/tmp/qcbench/bin/python /tmp/bench_stim_full.py

# QuantumClifford.jl:
JULIA_DEPOT_PATH=/tmp/qc_depot julia /tmp/bench_qc_full.jl

# Full stabilizer test suite:
wolframscript -file /tmp/run_stab_tests.wl
```

Op streams: `OngoingProjects/Platform Comparison/scripts/stab_ops_{100,500,1000}.json`.
