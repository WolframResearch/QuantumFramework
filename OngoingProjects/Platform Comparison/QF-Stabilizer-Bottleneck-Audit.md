# QF stabilizer subsystem: performance bottleneck audit, with resolution status

**Original audit:** 2026-06-11, written after round 1 (the packed gate fast path, `QF-Stabilizer-Packed-FastPath.md`), when QF was still ~70-110x Stim on gate-heavy Clifford simulation. The audit found *why* and ranked the fixes.
**This revision:** 2026-06-12, at repo HEAD `02cc92d9`. The audit's recommendations have since shipped as rounds 2 and 3 (`e37b7fad`, `730abcec`; see `QF-Stabilizer-Optimization-Report.md`), so this document now carries a **status marker on every item** and compares each projection against the measured outcome. The analysis sections are kept intact: they explain *why* the kernel is built the way it now is. Numbers re-verified at HEAD on the revision date are marked as such. **Machine:** Apple M2 Pro, Mathematica 15.0.0. **Scope:** `QuantumFramework/Kernel/Stabilizer/` (now 20 files, including the round-2 additions `Packed.m` and `Compiled.m`).

## TL;DR (as resolved)

The audit's diagnosis: round-1 packing removed the **algorithmic** $O(n^2)$-per-gate cost, and what remained was a **constant-factor interpreter tax**, every gate being ~12 separate vectorized calls plus object wrap/unwrap, each call costing ~1 µs of Wolfram-kernel dispatch *regardless of array length*.

The audit's one recommendation that mattered: **compile the gate fold** into a single native (C) routine. That shipped as `Compiled.m` / `ps["ApplyCircuit"]` and performed *better* than projected: the compiled kernel runs at $0.22$ to $0.27\,\mu\mathrm{s}$/gate (projection said $1.2$), and end-to-end QF landed at $2.5$ to $5.3\times$ Stim (projection said $3$ to $8\times$). The secondary items (formatting, measurement, composition crash, property reconstruction) also shipped, except the full packed-composition rewrite, which was deliberately dropped (status board below).

## Status board (2026-06-12, HEAD `02cc92d9`)

| Audit item | Recommendation | Status | Measured outcome |
|---|---|---|---|
| §2 headline | compile the gate fold to C, wire into `PauliStabilizerApply`, expose `ps["ApplyCircuit"]` | **SHIPPED** (round 2, `Compiled.m`) | $0.22$-$0.27\,\mu\mathrm{s}$/gate compiled kernel; end-to-end $2.6/26.2/58.0\,\mathrm{ms}$ at $n=100/500/1000$ (re-verified at HEAD: $3.6/25.5/57.0$) |
| §2 fallback | `"WVM"` bytecode target when no C compiler | **SHIPPED** (auto-select via `CCompilerDriver`) | guard in place; not exercised end-to-end on this machine (clang present) |
| 3a `PauliForm` | vectorize the Pauli-letter decode | **SHIPPED** (round 2, `Formatting.m`) | $128.6 \to 2.0\,\mathrm{ms}$ at $n=100$ (re-verified: $2.2$) |
| 3b composition | packed symplectic product + carry-save phase; fix the `Rest[{}]` crash | **CRASH FIX SHIPPED; rewrite dropped** | crash fixed (`Compose.m:35`); $2577 \to 46\,\mathrm{ms}$ at $n=100$ anyway (re-verified: $45.4$), see note below |
| 3c single-qubit measurement | packed/compiled measurement kernel, no tableau reconstruction | **SHIPPED** (round 2, `Measurement.m`) | $6.85 \to 0.10\,\mathrm{ms}$ at $n=100$, $16.5 \to 0.23$ at $n=500$ (re-verified: $0.10/0.19$); plus two correctness fixes the work surfaced (below) |
| 3d Pauli-string measurement | symplectic product on packed rows; kill the $2^n$ deterministic fallback | **SHIPPED** (round 2 deterministic, round 3 non-deterministic, `PauliMeasure.m`) | deterministic branch now closed-form $O(n^2)$; non-deterministic branch $393\,\mathrm{ms} \to 4.9$ at $n=100$, $42.3\,\mathrm{s} \to 58.7\,\mathrm{ms}$ at $n=500$ (re-verified: $4.1/50.3$) |
| 3e `ps["Circuit"]` synthesis | rides the compiled gate path automatically | **TRANSITIVE** | unchanged algorithm ($O(n^3)$ gates), now fed by fast gate updates |
| 3f property reconstruction | memoize the reconstructed tableau on the packed object | **SUPERSEDED** | not memoized; instead `Qubits`/`Qudits`/`GeneratorCount` read the stored count (`Properties.m:9`) and measurement/formatting operate **directly on packed words**, so the repeated $O(n^2)$ rebuilds this item targeted no longer happen on the hot paths |

**Note on 3b.** The packed phase-convention rewrite was prototyped but could not be verified against the existing `phaseLookup` semantics, so per the exact-equality discipline it was **not shipped**; only the verified crash fix was. The $55\times$ improvement happened anyway, because packed states now reconstruct to *packed* arrays, which the existing phase loop's `Pick`/`FoldList` consume far more efficiently. The algorithm is still the $O(n^3)$ symplectic multiply, unchanged and correct.

**Correctness dividends.** The measurement work (3c/3d) was gated on exact equivalence with the canonical path and on dense-simulation/Stim cross-checks, which caught two **pre-existing** kernel bugs (Aaronson-Gottesman row-sum omitting the destabilizers, `1c524432`; the same omission in the Pauli-string branch leaving a non-abelian "stabilizer" group, `52d3ea89`) and one latent round-2 bug (a destabilizer-phase convention mismatch that silently knocked branches off the fast path, fixed in `730abcec`). Plain-language account: master comparison, Section 12.1. Stabilizer suite at HEAD: 965/965; full framework suite 1493/1493.

---

The sections below are the original audit (2026-06-11), kept as the rationale for the design that shipped, with status notes inline.

## Method

`AbsoluteTiming`, `ClearSystemCache[]` before each run, min-of-3, `SetAttributes[..., HoldFirst]` on the timer (the one trap that silently reports 0 ms). Per-gate decomposition isolates each helper by calling it 1000x on a warmed packed object.

## 1. Where a single gate's time goes (the headline)

Per-gate cost of the round-1 packed path, decomposed (µs). *Status note: this path still exists as the single-gate route (`ps["H", j]` etc.) and these numbers still hold at HEAD (re-verified 2026-06-12: $37/40/48\,\mu\mathrm{s}$/gate on the benchmark stream at $n=100/500/1000$); the bulk route of §2 is what removed the tax for circuits.*

| component | n=100 | n=500 | n=1000 | what it is |
|---|---:|---:|---:|---|
| property dispatch | 4.4 | 4.4 | 5.1 | matching `ps["H",j]` through the downvalue table |
| `psConcreteFastQ` | 2.1 | 2.1 | 2.1 | the fast-path predicate |
| `psGetPacked` | 3.0 | 2.9 | 5.0 | object to `{n,px,pz,s}` tuple |
| **`packedGateH` kernel** | **15.3** | **18.8** | **30.6** | the actual bit algebra |
| `psFromPacked` | 6.4 | 6.1 | 7.5 | tuple to object |
| **total H gate** | **31.0** | **34.7** | **52.1** | |
| total CNOT gate | 46.4 | 47.5 | 55.8 | |

Two equally-large culprits:

- **Object overhead ~16 µs/gate** (dispatch + predicate + getPacked + fromPacked). This is paid *per gate* only because the public API folds `PauliStabilizer` objects one gate at a time. A bulk path that converts once amortizes it to ~0.
- **Kernel ~15-30 µs/gate.** `packedGateH` is ~12 Listable ops (`BitAnd`, `Sign`, `BitXor`, subtract/add, `ReplacePart`, the `s(1-2...)` sign multiply) over the length-$2n$ vectors. Each Listable call carries ~0.5-1 µs of interpreter dispatch *independent of length*; 12 of them is the observed cost. It grows with $n$ (200 to 2000 elements) but the floor is the call count, not the element count.

This is the structural ceiling of an interpreted packed path. You cannot Listable your way under ~5-10 µs/gate.

## 2. The fix, validated: compile the fold. **SHIPPED as `Compiled.m` (round 2)**

A single `Compile[..., CompilationTarget -> "C"]` function that takes the packed `Int64` arrays + an integer-encoded gate list and loops over all gates **in compiled code** (no per-gate interpreter dispatch, no object wrapping). Prototype validation at audit time:

```
n = 50, 2000-gate stream (H/S/CNOT):
  current interpreted packed fold : 95.0 ms   (47 µs/gate)
  compiled-to-C bulk kernel       :  2.3 ms   ( 1.2 µs/gate)   -> 41x faster
  correctness: compiled X, Z, signs  ===  interpreted packed path  (exact)
```

The compiled kernel reproduces the canonical tableau and signs **bit-for-bit** (verified against the interpreted packed path on a 2000-gate random Clifford stream).

**Why this works:** the bit algebra is already pure integer/bitwise ops (`BitAnd`/`BitXor`/`BitShiftLeft`/`Sign`), all of which `Compile` lowers to native C. The interpreter tax (both the ~16 µs object overhead and the ~1 µs-per-Listable-call kernel tax) vanishes because the entire gate loop becomes one C function call.

### Concrete design (all five points implemented; deltas noted)

1. **Encode** the circuit as an integer gate table. *Shipped as `encodeStabilizerGates`, a single `Dispatch`-table pass at ~0.7 µs/spec; the gate codes live in one association (`$stabilizerGateCodes`) injected into both the encoder and the compiled body. CZ is additionally accepted (expanded to H·CNOT·H); any unrecognized gate fails the whole encode, routing the caller to the per-gate fold so non-Clifford circuits behave exactly as before.*
2. **Convert** the input state to packed once. *Shipped; the shipped layout is generator-major (per-qubit columns packed across all $2n$ generators), the same data layout Stim and QC.jl use.*
3. **Call** the compiled kernel once, multi-chunk for $n > 62$. *Shipped; measured $0.22$-$0.27\,\mu\mathrm{s}$/gate, flat in $n$, beating the $1.2\,\mu\mathrm{s}$ prototype because the shipped kernel avoids the prototype's per-gate array copies.*
4. **Rebuild** the object once. *Shipped; note `ApplyCircuit` returns a canonical (rank-3) object while single gates return packed ones; both are accepted everywhere, but compare via `["Tableau"]`/`["Stabilizers"]`, not `SameQ` on raw objects.*
5. **Wire** it into `PauliStabilizerApply` so every `qc[Method -> "Stabilizer"]` benefits, and expose `ps["ApplyCircuit", gates]`. *Shipped (`PauliStabilizer.m:123`). The same commit train also fixed a pre-existing register-sizing bug in `PauliStabilizerApply` (wire count vs highest wire index, `1c524432`).*

Fallback: if no C compiler is present, `CompilationTarget -> "WVM"` (bytecode) still removes the per-call dispatch. *Shipped as an auto-select guard; the WVM branch has not been exercised end-to-end on this machine (clang present). A packed-array trap recorded during integration: a compiled function silently falls back to interpreted evaluation (~100x slower, zero warning) if any argument arrives unpacked, so every boundary forces `Developer`ToPackedArray`.*

Projected vs actual end-to-end (total wall-clock ms, identical op stream; actual measured at `6f953ad5`, re-verified at HEAD `02cc92d9`):

| n (gates) | QF at audit time | QF projected | **QF actual** | Stim |
|---:|---:|---:|---:|---:|
| 100 (2k) | 75 ms | ~3-5 ms | **2.6 ms** | 1.03 ms |
| 500 (10k) | 489 ms | ~20-40 ms | **26.2 ms** | 5.14 ms |
| 1000 (20k) | 1255 ms | ~50-90 ms | **58.0 ms** | 10.9 ms |

i.e. the projection "from ~70-110x Stim down to ~3-8x Stim" resolved to **2.5-5.3x Stim**, and ahead of QuantumClifford.jl at $n=1000$ ($58.0$ vs $108.4\,\mathrm{ms}$). The remaining gap to Stim is SIMD lane width (256-bit lanes vs 62-bit words) and circuit-encoding overhead, not algorithm. LibraryLink (the §11 option of the original kernel audit) would close part of the last factor and remains unimplemented and not currently needed.

## 3. Secondary bottlenecks (ranked); statuses in the board above

Measured at n=100 (random Clifford state), with the dominant cause:

| # | operation | n=100 cost at audit time | complexity | root cause | fix (as proposed; see status board) |
|---|---|---:|---|---|---|
| 3a | **`ps["Stabilizers"]` / `PauliForm`** | **128 ms** | $O(n^2)$ | `Formatting.m:24` built $n$ strings of length $n$ via per-element `Replace[Transpose[...], {...}, {2}]` then `StringJoin@@@`; the rank-3 `Transpose` and elementwise `Replace` were the cost | map the 4 symbol codes with one `Part`-index into a char array (a Pauli letter is just $x + 2z$). **Shipped: 2.0 ms** |
| 3b | **composition `left[right]`** | **2577 ms** | $O(n^3)$ | `Compose.m:35`: the per-row phase loop (`phaseLookup` + `Pick`/`FoldList`/`Extract`), *not* the matmul, dominated; also a latent `Rest[{}]` crash on any row selecting 0 or 1 generators | packed row-XOR product + carry-save phase; at minimum fix the crash. **Crash fix shipped; rewrite dropped (unverifiable vs `phaseLookup`); 46 ms anyway via packed reconstruction** |
| 3c | **single-qubit measurement** | **6.85 ms** | $O(n^2)$/call | `Measurement.m` reconstructed the full rank-3 array from packed *every call* then ran interpreted `rowsum` folds | packed measurement kernel (AG destabilizer trick on chunks, no reconstruction). **Shipped: 0.10 ms; an $n$-qubit sweep went from $O(n^3)$ to $O(n^2/62)$** |
| 3d | Pauli-string measurement | 35 ms | $O(n^2)$ | `PauliMeasure.m` built `ArrayFlatten[omega]` ($2n \times 2n$) and `gens.omega.pVec` per call; the deterministic branch even materialized the $2^n$ state vector | reuse the symplectic product on packed rows; use AG i-factor tracking instead of the $2^n$ fallback. **Shipped (deterministic round 2, non-deterministic round 3): 4.9 ms.** *Caveat: the 35 ms figure was measured on the pre-`52d3ea89` branch, which was producing wrong post-states; the corrected interpreted branch cost 393 ms before round 3 packed it.* |
| 3e | `ps["Circuit"]` synthesis | $O(n^3)$ | $O(n^3)$ | `Conversions.m:91` greedy AG: ~$n^2$ gate applications, each a full object gate-update | rides the compiled gate path automatically. **Transitive** |
| 3f | `ps["Matrix"]`, `ps["GeneratorCount"]`, `ps["Qudits"]` on packed | 3.6 ms (Matrix) | $O(n^2)$/call | each reconstructed the tableau from packed on every call; no memoization | cache the reconstructed `"Tableau"` on first access. **Superseded: dimension handlers read the stored count (`Properties.m:9`) and the hot paths (measurement, formatting) now work directly on packed words, so the repeated rebuilds are gone where they mattered** |

Not bottlenecks (inherent or already polynomial-good): `["State"]` (`Conversions.m`, $O(2^n)$: fundamental, $n \lesssim 10$), `InnerProduct` closed-form ($O(n^3)$, fine), `Entropy` ($O(n^3)$ rank, 5 ms), `RandomClifford` ($O(n^3)$, one-shot; *its $n \approx 500$ blow-ups and the $n \gtrsim 538$ underflow were separate correctness items, fixed in `52d3ea89` and `02cc92d9`*), `CliffordChannel` composition ($\mathbb{F}_2$ nullspace, channel formalism: separate concern), `GraphState`/`LocalComplement` (graph ops). `["Tableau"]` reconstruction itself is cheap (0.98 ms at n=100); the cost was calling it *repeatedly*.

## 4. The through-line

Every secondary item shares the root cause of #1: **interpreted per-element / per-call work that should be one compiled or vectorized pass.** Formatting re-scanned element-by-element; measurement reconstructed and folded in the interpreter; composition loops the phase row-by-row; properties rebuilt instead of reading stored data. The packing fixed the *data layout*; the remaining cost was the *interpreter touching that data too many times*.

The priority order the audit set, with resolutions:

1. **Compile the gate fold** (§2). *Shipped round 2; validated 41x prototype became ~$10^3\times$ end-to-end against the original implementation.*
2. **Stop the repeated tableau rebuilds** (3f). *Superseded by packed-native hot paths and stored-count handlers.*
3. **Vectorize `PauliForm`** (3a). *Shipped, 63x.*
4. **Compiled/packed measurement kernel** (3c/3d). *Shipped rounds 2-3, 68x and 80-720x, plus three correctness fixes.*
5. **Packed composition + the `Rest[{}]` fix** (3b). *Crash fix shipped; phase-convention rewrite consciously dropped; 55x landed anyway.*

### What remains open (the honest residue, 2026-06-12)

- **Gap to Stim ($2.5$ to $5.3\times$):** SIMD lane width and the ~0.7 µs/gate circuit-encoding pass. A LibraryLink kernel could close part of it; not currently justified.
- **Constant-cost Pauli-frame bulk sampling** (the master comparison's Section 3.2 gap): out of this audit's scope, still absent.
- **Packed composition phase rewrite** (3b): would need a verified mapping to `phaseLookup` semantics first.
- **WVM fallback path:** implemented but untested end-to-end (no compiler-free machine at hand).
- **`QuantumCircuitOperator` construction overhead** for $10^4$-gate circuits (minutes): a host-framework issue, outside the stabilizer subsystem; use `ps["ApplyCircuit", specs]` / `PauliStabilizerApply` for long streams.

## Appendix: raw data

Per-gate decomposition (2026-06-11, round-1 path):

```
n=100  H 31.0us  CNOT 46.4us  | dispatch 4.4  fastQ 2.1  getPacked 3.0  kernelH 15.3  fromPacked 6.4
n=500  H 34.7us  CNOT 47.5us  | dispatch 4.4  fastQ 2.1  getPacked 2.9  kernelH 18.8  fromPacked 6.1
n=1000 H 52.1us  CNOT 55.8us  | dispatch 5.1  fastQ 2.1  getPacked 5.0  kernelH 30.6  fromPacked 7.5
compiled-to-C bulk fold prototype, n=50, 2000 gates: 2.3 ms (1.2 us/gate) vs 95 ms interpreted; X/Z/signs exact match
other ops @ n=100: Stabilizers 128ms  single-measure 6.85ms  Pauli-string-measure 35ms  Matrix 3.6ms  Entropy 5.2ms  Tableau-recon 0.98ms  compose 2577ms
compose scaling: n=50 179ms, n=100 2577ms, n=200 10085ms  (O(n^3))
```

Re-verification at HEAD `02cc92d9` (2026-06-12, `scripts/qf_final_all.wl`, min-of-3, `ClearSystemCache[]` per rep):

```
ApplyCircuit (bulk, shipped §2)   n=100: 3.58 ms   n=500: 25.5 ms   n=1000: 57.0 ms
per-gate fold (round-1 path, §1)  n=100: 74.2 ms   n=500: 402 ms    n=1000: 955 ms
PauliStabilizer[n] constructor    n=100: 0.071 ms  n=500: 1.24 ms   n=1000: 4.9 ms
ps["M", 1]                        n=100: 0.103 ms  n=500: 0.187 ms
ps["M", "Z..Z"]                   n=100: 4.07 ms   n=500: 50.3 ms
ps["Stabilizers"]                 n=100: 2.2 ms
composition a[b]                  n=100: 45.4 ms
```

Re-confirmed at HEAD `f9dc1cdc` (2026-06-14): the property-cache refactor of 2026-06-13/14 (`dfc741dc`/`f9dc1cdc`, the state-vector apply-path fix) touches `QuantumOperator`/`QuantumState`/`QuantumBasis` property dispatch, **not** the `Stabilizer/` subsystem, and the stabilizer numbers reproduce within run-to-run noise:

```
ApplyCircuit (bulk)               n=100: 2.81 ms   n=500: 25.9 ms   n=1000: 62.7 ms
per-gate fold (round-1 path)      n=100: 71.1 ms   n=500: 406 ms    n=1000: 948 ms
PauliStabilizer[n] constructor    n=100: 0.062 ms  n=500: 1.51 ms   n=1000: 5.08 ms
ps["M", 1]                        n=100: 0.099 ms  n=500: 0.177 ms
ps["M", "Z..Z"]                   n=100: 3.98 ms   n=500: 49.9 ms
ps["Stabilizers"]                 n=100: 1.92 ms
composition a[b]                  n=100: 48.7 ms
```
