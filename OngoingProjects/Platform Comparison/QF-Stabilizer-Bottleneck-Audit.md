# QF stabilizer subsystem: full performance bottleneck audit

**Date:** 2026-06-11 · **Machine:** Apple M2 Pro, Mathematica 15.0.0 · **Scope:** all 18 files in `QuantumFramework/Kernel/Stabilizer/`
**State:** packed gate fast path already landed (see `QF-Stabilizer-Packed-FastPath.md`). QF is still ~70-110× Stim on gate-heavy Clifford simulation. This audit finds *why* and ranks the fixes.

## TL;DR

The packing removed the **algorithmic** $O(n^2)$-per-gate cost. What remains is a **constant-factor interpreter tax**: every gate is ~12 separate Listable calls plus object wrap/unwrap, each costing ~1 µs of Wolfram-kernel dispatch *regardless of array length*. At 31-56 µs/gate that tax, not arithmetic, is now the wall.

**The one fix that matters: compile the gate fold.** A `Compile`-to-C bulk kernel, validated below, runs the identical bit algebra at **1.2 µs/gate, a 40× speedup over the current packed path**, landing within ~2× of Stim. Everything else in this report is secondary.

## Method

`AbsoluteTiming`, `ClearSystemCache[]` before each run, min-of-3, `SetAttributes[…, HoldFirst]` on the timer (the one trap that silently reports 0 ms). Per-gate decomposition isolates each helper by calling it 1000× on a warmed packed object.

## 1. Where a single gate's time goes (the headline)

Per-gate cost of the *current* packed path, decomposed (µs):

| component | n=100 | n=500 | n=1000 | what it is |
|---|---:|---:|---:|---|
| property dispatch | 4.4 | 4.4 | 5.1 | matching `ps["H",j]` through the downvalue table |
| `psConcreteFastQ` | 2.1 | 2.1 | 2.1 | the fast-path predicate |
| `psGetPacked` | 3.0 | 2.9 | 5.0 | object → `{n,px,pz,s}` tuple |
| **`packedGateH` kernel** | **15.3** | **18.8** | **30.6** | the actual bit algebra |
| `psFromPacked` | 6.4 | 6.1 | 7.5 | tuple → object |
| **total H gate** | **31.0** | **34.7** | **52.1** | |
| total CNOT gate | 46.4 | 47.5 | 55.8 | |

Two equally-large culprits:

- **Object overhead ≈ 16 µs/gate** (dispatch + predicate + getPacked + fromPacked). This is paid *per gate* only because the public API folds `PauliStabilizer` objects one gate at a time. A bulk path that converts once amortizes it to ~0.
- **Kernel ≈ 15-30 µs/gate.** `packedGateH` is ~12 Listable ops (`BitAnd`, `Sign`, `BitXor`, subtract/add, `ReplacePart`, the `s(1-2…)` sign multiply) over the length-$2n$ vectors. Each Listable call carries ~0.5-1 µs of interpreter dispatch *independent of length*; 12 of them ≈ the observed cost. It grows with $n$ (200→2000 elements) but the floor is the call count, not the element count.

This is the structural ceiling of an interpreted packed path. You cannot Listable your way under ~5-10 µs/gate.

## 2. The fix, validated: compile the fold

A single `Compile[…, CompilationTarget -> "C"]` function that takes the packed `Int64` arrays + an integer-encoded gate list and loops over all gates **in compiled code** (no per-gate interpreter dispatch, no object wrapping):

```
n = 50, 2000-gate stream (H/S/CNOT):
  current interpreted packed fold : 95.0 ms   (47 µs/gate)
  compiled-to-C bulk kernel       :  2.3 ms   ( 1.2 µs/gate)   -> 41× faster
  correctness: compiled X, Z, signs  ===  interpreted packed path  (exact)
```

The compiled kernel reproduces the canonical tableau and signs **bit-for-bit** (verified against the interpreted packed path on a 2000-gate random Clifford stream). 1.2 µs/gate is ~2× Stim's 0.5 µs/gate — the remaining gap is SIMD and Stim's hand-tuned inner loop, not algorithm.

**Why this works:** the bit algebra is already pure integer/bitwise ops (`BitAnd`/`BitXor`/`BitShiftLeft`/`Sign`), all of which `Compile` lowers to native C. The interpreter tax (both the ~16 µs object overhead and the ~1 µs-per-Listable-call kernel tax) vanishes because the entire gate loop becomes one C function call.

### Concrete design

1. **Encode** the circuit as an `Integer` matrix `gates` with rows `{code, j, k}` (`code`: 0=H, 1=S, 2=Sdg, 3=CNOT, 4=SWAP, 5/6/7=X/Y/Z; `j,k` 0-based bit positions). Single pass over `QuantumShortcut[qco]`.
2. **Convert** the input state to packed once.
3. **Call** the compiled kernel once. For $n \le 62$ it is the rank-1 `Int64` version above. For $n > 62$ use the rank-2 `{nchunks, 2n}` version: same body, the chunk index $\lceil j/62 \rceil$ selects the column; per-gate work is $O(\text{nchunks})$ in compiled C (still ~1-3 µs/gate).
4. **Rebuild** the object once.
5. **Wire** it into `PauliStabilizerApply` (`PauliStabilizer.m:112`) so *every* `qc[Method -> "Stabilizer"]` and circuit-fold benefits, and expose a direct `ps["ApplyCircuit", gates]` for the object-fold microbenchmark path.

Fallback: if no C compiler is present, `CompilationTarget -> "WVM"` (bytecode) still removes the per-call dispatch and should give ~5-10×. Guard with `compiledFunctionLoadableQ`; fall back to the current interpreted packed path otherwise. LibraryLink (the audit's §11 option) would close the last 2× to Stim but is not needed to satisfy "not happy yet."

Projected end-to-end (scaling 1.2 µs/gate, + multi-chunk factor):

| n (gates) | QF now | QF compiled (proj.) | Stim |
|---:|---:|---:|---:|
| 100 (2k) | 75 ms | ~3-5 ms | 1.0 ms |
| 500 (10k) | 489 ms | ~20-40 ms | 5.5 ms |
| 1000 (20k) | 1255 ms | ~50-90 ms | 11.3 ms |

i.e. from ~70-110× Stim down to ~3-8× Stim.

## 3. Secondary bottlenecks (ranked)

Measured at n=100 (random Clifford state), with the dominant cause:

| # | operation | n=100 cost | complexity | root cause | fix |
|---|---|---:|---|---|---|
| 3a | **`ps["Stabilizers"]` / `PauliForm`** | **128 ms** | $O(n^2)$ | `Formatting.m:24` builds $n$ strings of length $n$ via per-element `Replace[Transpose[...], {...}, {2}]` then `StringJoin@@@`; the rank-3 `Transpose` and elementwise `Replace` are the cost | map the 4 symbol codes with one `Part`-index into a char array, `StringJoin` per row once; or decode straight from packed bits. ~10-50× |
| 3b | **composition `left[right]`** | **2577 ms** | $O(n^3)$ | `Compose.m:35` — the per-row phase loop (`phaseLookup` + `Pick`/`FoldList`/`Extract`), *not* the matmul, dominates; also a latent `Rest[{}]` crash on an all-zero Pauli row | packed row-XOR product + carry-save phase (audit §4.2); at minimum vectorize/compile the phase loop. Fixes `V`/`V†`, `Dagger` inner, explicit compose |
| 3c | **single-qubit measurement** | **6.85 ms** | $O(n^2)$/call | `Measurement.m` reconstructs the full rank-3 array from packed *every call* (`ps["Tableau"]`, `ps["Signs"]`, `ps["GeneratorCount"]`) then runs interpreted `rowsum` folds | packed/compiled measurement kernel (AG destabilizer trick on chunks, no reconstruction). A measurement-per-qubit sweep is currently $O(n^3)$ |
| 3d | Pauli-string measurement | 35 ms | $O(n^2)$ | `PauliMeasure.m`: builds `ArrayFlatten[ω]` ($2n\times2n$) and `gens.ω.pVec` per call; deterministic branch even materializes the $2^n$ state vector (`PauliMeasure.m:59`) | reuse the symplectic product on packed rows; avoid the $2^n$ fallback (use AG i-factor tracking, already done in `InnerProduct.m`) |
| 3e | `ps["Circuit"]` synthesis | $O(n^3)$ | $O(n^3)$ | `Conversions.m:91` greedy AG: ~$n^2$ gate applications, each a full object gate-update | runs on top of the compiled gate path automatically once 1.2 is wired; otherwise unchanged-correct |
| 3f | `ps["Matrix"]`, `ps["GeneratorCount"]`, `ps["Qudits"]` on packed | 3.6 ms (Matrix) | $O(n^2)$/call | each reconstructs the tableau from packed every call; no memoization | cache the reconstructed `"Tableau"` into the packed association on first access (the head is an immutable carrier; store it once). Removes repeated $O(n^2)$ rebuilds inside measurement/composition/formatting |

Not bottlenecks (inherent or already polynomial-good): `["State"]` (`Conversions.m:69`, $O(2^n)$ — fundamental, n≤10), `InnerProduct` closed-form (`InnerProduct.m:50`, $O(n^3)$, fine), `Entropy` (`Entropy.m`, $O(n^3)$ rank, 5 ms), `RandomClifford` ($O(n^3)$, one-shot), `CliffordChannel` composition (F₂ nullspace, channel formalism — separate concern), `GraphState`/`LocalComplement` (graph ops). `["Tableau"]` reconstruction itself is cheap (0.98 ms at n=100); the cost is calling it *repeatedly*.

## 4. The through-line

Every secondary item shares the root cause of #1: **interpreted per-element / per-call work that should be one compiled or vectorized pass.** Formatting re-scans element-by-element; measurement reconstructs and folds in the interpreter; composition loops the phase row-by-row; properties rebuild instead of caching. The packing fixed the *data layout*; the remaining cost is the *interpreter touching that data too many times*.

So the priority order is unambiguous:

1. **Compile the gate fold** (§2). Validated 40×. Fixes the headline and, transitively, Circuit synthesis (3e).
2. **Memoize the reconstructed tableau** in the packed association (3f). One-line-ish; removes repeated $O(n^2)$ rebuilds feeding measurement/formatting/composition.
3. **Vectorize `PauliForm`** (3a). 128 ms → few ms; pure display path, low risk.
4. **Compiled/packed measurement kernel** (3c/3d). Needed if the workload measures a lot.
5. **Packed composition + the `Rest[{}]` fix** (3b). Needed for `V`-heavy or composition-heavy use.

Items 1-3 are high-value/low-risk and would take QF from "~100× Stim" to "single-digit× Stim" on the workloads the benchmark targets. I can implement #1 (the compiled kernel, with the WVM fallback and `PauliStabilizerApply` wiring) next if you want — the prototype is already correctness-validated.

## Appendix: raw per-gate decomposition data

```
n=100  H 31.0us  CNOT 46.4us  | dispatch 4.4  fastQ 2.1  getPacked 3.0  kernelH 15.3  fromPacked 6.4
n=500  H 34.7us  CNOT 47.5us  | dispatch 4.4  fastQ 2.1  getPacked 2.9  kernelH 18.8  fromPacked 6.1
n=1000 H 52.1us  CNOT 55.8us  | dispatch 5.1  fastQ 2.1  getPacked 5.0  kernelH 30.6  fromPacked 7.5
compiled-to-C bulk fold, n=50, 2000 gates: 2.3 ms (1.2 us/gate) vs 95 ms interpreted; X/Z/signs exact match
other ops @ n=100: Stabilizers 128ms  single-measure 6.85ms  Pauli-string-measure 35ms  Matrix 3.6ms  Entropy 5.2ms  Tableau-recon 0.98ms  compose 2577ms
compose scaling: n=50 179ms, n=100 2577ms, n=200 10085ms  (O(n^3))
```
