# Audit: QuantumClifford.jl and Stim against QuantumFramework's `PauliStabilizer`

**Date:** 2026-04-30
**Subject:** Two reference implementations of Clifford-circuit / stabilizer-state simulation, audited as **upgrade material** for the existing `QuantumFramework/Kernel/PauliStabilizer.m` (494 LOC, paclet 1.6.5, SHA cbbe9368, last modified 2024-05-29).
**Repos audited (external):**
- `External Packages/QuantumClifford.jl` — pure-Julia, Krastanov / QuantumSavory.
- `External Packages/Stim` — Google Quantum AI / Gidney, C++ with pybind11.

**Anchor for the in-tree comparison:** `audit/Stabilizer/paulistabilizer-source-audit.md` — deep audit of `PauliStabilizer.m` with file:line citations for every claim made about current behavior. Read alongside this document.

## 0. Executive summary

QuantumFramework already ships a working stabilizer simulator. The kernel implements the Aaronson–Gottesman tableau formalism, all standard single- and two-qubit Cliffords, symplectic-mod-2 composition, tensor product, AG canonical-circuit synthesis, Z-basis measurement (deterministic and random branches), Mallows-distribution random Clifford sampling, and bidirectional conversion to/from `QuantumState`, `QuantumOperator`, and `QuantumCircuitOperator`. By the standards of either reference implementation it is roughly a feature-complete v0.1.

The gaps are **specific and well-localized**, and they cluster into four bands:

1. **Bulk-sampling throughput** — Pauli-frame / reference-frame sampling is the single feature both QC.jl and Stim treat as the reason their packages exist; `PauliStabilizer.m` has no analogue. Per-shot tableau measurement is the only sampling path. This is the largest functional gap.
2. **Measurement surface** — only Z-basis on a single qubit is supported; arbitrary Pauli-string measurement (the syndrome-extraction primitive) and expectation values `⟨P⟩` for arbitrary Pauli `P` are absent. Both are short additions on top of the existing AG `rowsum` machinery.
3. **Symbolic extension** — the entire data structure commits to integer entries (`PauliTableauQ`, line 12) and integer signs `∈ {-1, +1}` (line 13). Non-Clifford rotations (`"P"[θ]`, `"T"`) escape the formalism by returning a sum of two `PauliStabilizer` objects with complex coefficients (line 321) — but the sum is not closed under further gate application, so it is a one-shot expansion rather than a real symbolic mode. This is the area the project owner has explicitly flagged as priority.
4. **Performance path** — composition (line 389) is dense matrix-mod-2 multiplication; the rank-3 tableau stores `{0, 1, -1}` integers, not packed UInt64s. Both QC.jl and Stim spend essentially all of their compute on bit-packed XOR + popcount over packed Pauli rows; `PauliStabilizer.m` does not have this fast path.

The QC.jl/Stim audit material below is reorganized to map directly onto these four bands. For each algorithmic technique I identify whether `PauliStabilizer.m` already has it (in which case the reference implementations are useful for sanity-checking and corner cases), partially has it (where the reference shows the next refactor), or lacks it (where the reference is the implementation blueprint).

The most important strategic fact: the existing kernel's tableau layout is `{2, n, 2n}` with destabilizers in rows `1..n` and stabilizers in rows `n+1..2n` — i.e., already a destabilizer-augmented form, semantically equivalent to QC.jl's `MixedDestabilizer` with `rank == n` always. So the deterministic-measurement path is *already* O(n), and the inverted-tableau trick that Stim relies on does not need to be ported. This matters because it concentrates the upgrade work on items 1–4 above and rules out a class of structural rewrites that would otherwise look attractive.

## 1. Coverage and methodology

The two external repos contain ~870 source files combined. I read directly:
- `QuantumClifford.jl/src/`: the top module file, `pauli_operator.jl`, `mul_leftright.jl`, `canonicalization.jl`, `project_trace_reset.jl`, `pauli_frames.jl`, `symbolic_cliffords.jl`, plus the `docs/src/index.md` and `README.md`.
- `Stim/src/stim/`: `stabilizers/pauli_string.h`, `stabilizers/tableau.h`, `simulators/tableau_simulator.h`, plus the `README.md`.

Two parallel Explore agents mapped the rest of each repo (data layouts, kernels, gate dispatch, DEM, search, public API). All algorithmic claims below are cross-checked against either a direct read or an Explore citation; each citation is `path:line`.

For the in-tree comparison I rely on `audit/Stabilizer/paulistabilizer-source-audit.md`, which is a line-by-line audit of `QuantumFramework/Kernel/PauliStabilizer.m`. Where this document cites a `line N`, the reference is to that file unless a different path is given. The deep audit is the source of truth for what is currently implemented; this document is the analysis of how the QC.jl and Stim techniques map onto the existing code.

**Not covered** (deliberately, since they don't bear on the upgrade path):
- `lib/QECCore/`, `src/ecc/` (~100 files of code-family generators — Reed-Muller, BCH, Tanner, lifted-product, Gallager, color codes, surface, toric, Steane, Shor, etc.).
- `Stim/glue/` (Cirq, ZX, Crumble JS editor, Sinter/sample subprojects).
- `Stim/src/stim/diagram/`, `gen/`, `cmd/` (rendering, generators, command-line driver).
- The various Julia ext packages (`ext/QuantumCliffordHeckeExt`, `…OscarExt`, `…PyQDecodersExt`, `…PyTesseractDecoderExt`, `…JuMPExt`, `…MakieExt`, `…QOpticsExt`, `…QuantikzExt`).
- GPU paths (`ext/QuantumCliffordGPUExt`, `ext/QuantumCliffordKAExt`).
- Test directories on both sides.

If any of these turn out to matter for QuantumFramework (e.g. you want code-family generators, a Quantikz exporter, or LDPC decoder bindings), they're a separate audit. The 28-paper SDK literature in `Planning for research/Stabilizer/` covers most of these capabilities at the algorithmic level.

## 2. Three design philosophies, side by side

| Axis | Stim | QuantumClifford.jl | **`PauliStabilizer.m` (current)** |
|---|---|---|---|
| Primary state object | `TableauSimulator` holds one `Tableau<W> inv_state` (the **inverse** stabilizer tableau) | Four parallel types: `Stabilizer`, `Destabilizer`, `MixedStabilizer`, `MixedDestabilizer` | **One** type: `PauliStabilizer[<\|"Signs"->…, "Tableau"->…\|>]` with an always-present destabilizer half — semantically a `MixedDestabilizer` with `rank == n` always |
| Tableau layout | column-major, two `simd_bits<W>` arrays `xs` and `zs`, padded to W=256 | rank-2 chunked array `xzs` with X-chunks first, then Z-chunks; phase as a 0-d `UInt8` array | **rank-3 array of dimensions `{2, n, 2n}`**: first axis splits X/Z, second is qubit, third is row index (destabilizers `1..n`, stabilizers `n+1..2n`) |
| Phase storage | `bool sign` per row — only ±1; ±i is a transient `bool *imag` during multiplication | `UInt8` per row, two bits encoding {+1, +i, −1, −i} | **integer `Signs ∈ {-1, +1}` list of length 2n** (line 13); equivalent ±1 only, no ±i |
| Phase-bit access | one bit per row | two bits per row | one full integer per row (overspecified — three states `{−1, 0, +1}` is wider than needed but harmless) |
| Gate set | ~80 hand-written `prepend_*` and `do_*` kernels | ~30 symbolic gate types via `@qubitop1` / `@qubitop2` macros | **~12 named gate updates** (H, S, S†, X, Y, Z, CNOT, CZ, SWAP, V, V†, Permute, PermuteQudits) plus dagger and pad — all direct tableau mutations on lines 216–319 |
| Bulk sampling | `FrameSimulator` (reference-frame trick) | `PauliFrame` + `pftrajectories` (multithreaded) | **none** — measurement returns a conditional `Association`; user must sample by hand |
| Circuit IR | First-class: typed `Circuit` with monotonic buffers, REPEAT blocks, ASCII text round-trip | Generic vector of operations | Defers to `QuantumCircuitOperator`; integration via `PauliStabilizerApply` (line 426) |
| Detector / error model | `DetectorErrorModel` — separate compiled error-graph format | Implicit; PauliFrame trajectories are the substrate | not present |
| Random Clifford | uniform symplectic, not Bravyi–Maslov | uniform invertible symplectic, not Bravyi–Maslov | **Bravyi–Maslov / Koenig–Smolin Mallows-distribution sampler** (line 165) — `PauliStabilizer.m` is **ahead of both** here |
| API style | Imperative methods on simulator | Functional, `@P_str`, `@S_str`, `@C_str` literal macros | Functional, property-based: `ps[prop_String]` dispatches; `ps["H", j]` updates; `ps["M", q]` measures |
| Lines that matter | ~20K C++ in `src/stim/` core | ~8.4K Julia in `src/` | **494 lines of WL** in `Kernel/PauliStabilizer.m` |

The structural conclusion: `PauliStabilizer.m` is closest to QC.jl's `MixedDestabilizer` in semantics (always-present destabilizer half, rank tracking implicit), to Stim in API discipline (one type, properties as the user surface), and **ahead of both** on random-Clifford sampling (already uses Mallows). It lags both on bulk sampling, measurement breadth, and symbolic extensibility.

## 3. Pauli operator: the core data type

### 3.1 How `PauliStabilizer.m` stores Paulis today

The deep audit confirms (lines 11–13):

> "Tableau" in this file means the **rank-3 array** `Dimensions == {2, n, 2n}` — first axis splits X/Z bits, second axis is qubit index (n qubits), third axis is row index (2n rows: n destabilizers + n stabilizers, in that order).
> "Signs" is a length-`2n` integer list `∈ {-1, +1}` that tracks the ±1 phase per row.

Validation (line 12) requires entries to be `0 | 1 | -1` integers. The tableau is **not bit-packed** and entries are **not symbolic**; both are real constraints on what kinds of refactors are possible.

### 3.2 How QuantumClifford.jl stores Paulis

`External Packages/QuantumClifford.jl/src/pauli_operator.jl:42–48`:

```julia
struct PauliOperator{P <: AbstractArray{<: Unsigned, 0}, XZ <: AbstractVector{<: Unsigned}}
    phase::P            # 0-d array of UInt8 holding ∈ {0x0, 0x1, 0x2, 0x3}
    nqubits::Int
    xz::XZ              # X bits in xz[1:end÷2], Z bits in xz[end÷2+1:end]
end
```

The `xz` field is a **single concatenated array** with X chunks in the first half and Z chunks in the second. Cache-locality argument: `mul_left!`'s inner loop reads `l[i]`, `l[i+len]`, `r[i]`, `r[i+len]` — adjacent indices into the same array. The `phase` is a *zero-dimensional array* (not a scalar) so a GPU kernel can broadcast over it.

### 3.3 How Stim stores Paulis

`External Packages/Stim/src/stim/stabilizers/pauli_string.h:71–79`:

```cpp
template <size_t W>
struct PauliString {
    size_t num_qubits;
    bool sign;          // true = -1, false = +1; imaginary phase NOT permitted
    simd_bits<W> xs, zs;  // two SEPARATE simd_bits, each padded to W=256
};
```

Two key differences from QC.jl:
1. **Separate `xs` and `zs` arrays**, not concatenated. Each is independently 256-bit aligned for AVX.
2. **`sign` is a single bool** — ±1 only. Imaginary phase ±i is forbidden in the persistent type.

**Why Stim can get away with one bit:** A stabilizer generator is Hermitian, so its eigenvalues are ±1, never ±i. Tableaux that arise from conjugating Hermitian generators by Cliffords stay Hermitian — therefore the phase column is always ±1. Imaginary phase only enters as a side computation (`X·Z = −iY`), and gets cancelled before storage. QC.jl, in contrast, lets you build *arbitrary* Paulis with `P"-iXYZ"` — the full Pauli group — at the cost of one extra bit per row.

### 3.4 X/Y/Z encoding

All three implementations use the standard symplectic encoding:

| Pauli | (x, z) |
|---|---|
| I | (0, 0) |
| X | (1, 0) |
| Z | (0, 1) |
| Y | (1, 1) |

Stim provides explicit conversion helpers between this `xz` encoding and the alternative XYZ-as-{0,1,2,3} encoding for parsing/printing (`pauli_string.h:41–61`). `PauliStabilizer.m` uses the `xz` encoding throughout; the XYZ encoding only appears at the string-parsing boundary in the `PauliStabilizer[{"XZZXI", ...}]` constructor (line 109).

### 3.5 What this means for upgrades

`PauliStabilizer.m`'s rank-3 `{2, n, 2n}` integer tableau is **expressive** — every entry is a top-level value the user can inspect via `ps["Tableau"]` — but **not performant**: each entry is a full WL integer, and the per-bit arithmetic in gate updates is `Mod[..., 2]` rather than `BitXor`. Two distinct upgrade paths exist:

1. **Add a packed view alongside the current tableau.** Keep the rank-3 array as the canonical representation; provide a `ps["PackedXZ"]` accessor that returns a pair of `IntegerString`-backed packed UInt64 lists, computed lazily and cached. Hot kernels (composition, measurement) take a fast path that operates on the packed view; correctness kernels and display fall back to the rank-3 array. This preserves the user surface and is the lowest-risk path.

2. **Loosen the type contract to admit symbolic entries.** Replace `PauliTableauQ`'s `MatchQ[0 | 1 | -1]` constraint (line 12) with something that accepts symbolic expressions, gate by gate audit each update rule (lines 216–319) for whether `Mod[..., 2]` and the sign-flip predicates work symbolically. This is the path the project owner flagged as priority. It is **incompatible with bit-packing** in the obvious way, so the two paths bifurcate the implementation.

A reasonable architecture: keep the symbolic tableau as the default (extending the current code), and treat the packed-UInt64 path as a separate optimized representation triggered when entries are confirmed numeric (`Developer\`PackedArrayQ` or `MatchQ[..., 0 | 1]`). The dispatch is done at the boundary of expensive operations (composition, measurement, AG canonicalization) and is invisible to user code.

## 4. The multiplication kernel — the one that everything else rides on

This is the inner loop that the QC.jl README's claim "15μs to multiply two 1M-qubit Paulis" describes. It's also the only place where SIMD really matters.

### 4.1 The algebra

Two Paulis multiply via:
- X-bits: `x_out = x_l ⊕ x_r`
- Z-bits: `z_out = z_l ⊕ z_r`
- Phase: ±i for each qubit where the two inputs anticommute, summed mod 4.

The phase contribution per qubit is determined by `anti_comm = (x2 & z1) ⊕ (x1 & z2)` and a refinement bit. Both reference packages implement the same trick: **carry-save addition over two parity bits**, then `popcount` at the end.

QC.jl `External Packages/QuantumClifford.jl/src/mul_leftright.jl:11–31`, the non-vectorised reference:

```julia
function _mul_ordered_nonvec!(r, l; phases=true)
    if !phases
        r .⊻= l                    # XOR-only path when phase doesn't matter (Pauli frames)
        return (0, 0)
    end
    cnt1 = zero(T); cnt2 = zero(T)
    len = length(l) ÷ 2
    @inbounds @simd for i in 1:len
        x1, x2, z1, z2 = l[i], r[i], l[i+len], r[i+len]
        r[i]      = newx1 = x1 ⊻ x2
        r[i+len]  = newz1 = z1 ⊻ z2
        x1z2      = x1 & z2
        anti_comm = (x2 & z1) ⊻ x1z2
        cnt2     ⊻= (cnt1 ⊻ newx1 ⊻ newz1 ⊻ x1z2) & anti_comm   # carry
        cnt1     ⊻= anti_comm                                    # parity
    end
    rcnt1 = count_ones(cnt1)
    rcnt2 = count_ones(cnt2)
    return rcnt1, rcnt2
end
```

Phase combination at the call site:
```julia
phase_out = (rcnt1 ⊻ (rcnt2 << 1)) & 0x3   # mod 4
```

Stim does the exact same algebra in `src/stim/stabilizers/pauli_string_ref.inl`.

### 4.2 What `PauliStabilizer.m` does today

Two distinct code paths:

**(a) Pauli multiplication via the AG `g` and `rowsum` primitives** (lines 280, 282):
- `g[x1, z1, x2, z2]` is the standard AG symplectic phase function.
- `rowsum[{r, t}, h, i]` replaces row h with row i ⊕ row h, updating the sign according to `g`.

This is the per-row primitive used inside measurement (`ps["Measure", a]`, line 291) and exposed as `ps["RowSum", h, i]` (line 288). It works element-by-element on the rank-3 integer tableau — no bit-packing.

**(b) Composition `left[right]`** (line 389): full symplectic-mod-2 matrix multiplication. The deep audit notes:

> Performance: O(n³) for n qubits (matrix multiply dominated). No specialization for sparse tableaux.

This path uses `Inverse[Matrix, Modulus -> 2]` (line 263 for dagger) and dense `Dot` over the {0,1}-block matrix. Phase tracking is via `phaseLookup`, a precomputed sparse rank-4 array (line 380).

### 4.3 The `phases=false` fast path

Both QC.jl and Stim have a **conditional** that skips phase tracking entirely. In QC.jl, the `phases=false` branch above collapses the entire body to one XOR over the array. In Stim, the equivalent lives in `frame_simulator.cc` (the reference-frame sampling path, §9).

This is the single most important micro-optimisation in the Pauli-frame path: phase tracking is the expensive part, but Pauli frames don't need phase (they only track *which* Pauli error happened, not its global phase). Skipping the carry-save accounting roughly **halves** sampling cost per gate.

`PauliStabilizer.m` has no analogue, because it has no Pauli-frame path (§9). When that gets added, the `phases=false` branch is part of the spec.

### 4.4 What translates to WL

For a packed-UInt64 fast path on top of the existing rank-3 tableau:
- The XOR-only path (`phases=false`): one `BitXor` over packed lists, trivially fast in WL.
- The phased path: `BitAnd`, `BitXor`, and `DigitCount[..., 2, 1]` (popcount) — all native, called once per chunk. For `n ≤ 64` this is one chunk and runs in microseconds; for `n = 10⁶` it's ~16K chunks — the loop overhead in WL becomes the bottleneck unless compiled.
- LibraryLink option: ~30 lines of C exposing `mul_pauli_packed(uint64_t* l, uint64_t* r, size_t nchunks, uint8_t* phase_out)` would close the gap for large n.

**Recommendation:** don't write the LibraryLink shim until the WL-native packed path is shipped, profiled, and a real workload demands it. The current dense matrix-mod-2 multiply at line 389 is fine for the n ≲ 50 regime that fits a Wolfram notebook example; the packed path matters only when users start running surface-code-distance-7+ circuits.

## 5. State types: one is enough, but it should expose more

QC.jl ships **four** state types where most simulators ship one:

| Type | Storage | Adds | Use when… |
|---|---|---|---|
| `Stabilizer` | n rows × n qubits, single sign per row | nothing — bare list of generators | You only need to *track* a state, not measure it efficiently |
| `Destabilizer` | 2n rows | destabilizer half — measurements that anticommute become O(n) instead of O(n²) | You measure a lot; state is pure (rank n) |
| `MixedStabilizer` | up to n rows + explicit `rank::Int` | rank tracking — supports states of rank < n | You're modelling a code's syndrome readout where the rank changes |
| `MixedDestabilizer` | 2n rows + `rank` + logical X / logical Z slots | logicals — direct access to logical operators of a code | You're doing QEC and need logical Pauli observables |

The cost ramp is real: `Stabilizer` projection is **O(n³)** in the worst case (commuting measurement requires a full canonicalisation), `MixedDestabilizer` projection is **O(n²)**.

### 5.1 What `PauliStabilizer.m` is already

The current type is **structurally** a `MixedDestabilizer` with `rank == n` always:
- 2n rows always present (deep audit §1: "A valid object always has **2n rows**…destabilizers first, stabilizers second", lines 102-107).
- Single-row constructors auto-pad with the destabilizer half.
- Deterministic measurement is O(n) because the destabilizer half is there to be XOR-walked.

What it does **not** do:
- **Explicit rank tracking.** `rank < n` (mixed states, codes mid-readout) is not representable.
- **Logical operator slots.** Logicals are not exposed as a property; the user must canonicalise and read them off.
- **Type lattice.** A pure-state user pays for the destabilizer half they don't need, but the cost is one constant factor and probably not worth introducing a separate `Stabilizer` type.

### 5.2 Recommendation

Keep the single `PauliStabilizer` head; it is the right level of abstraction for a user-facing WL package. Add two things on top:

1. **A `"Rank"` property.** For now, return `n`. Refactor later to track rank explicitly when measurements rank-decrease the state (relevant for syndrome readout where some stabilizers get destroyed).
2. **`"LogicalX"` and `"LogicalZ"` properties.** The information is already in the destabilizer half once the state is in Gottesman canonical form. Expose it. This is the single most-asked-for QEC primitive that the current API doesn't have.

Both are short additions on top of the existing structure — no breaking changes to the data layout.

## 6. Canonicalisation — three flavors in QC.jl, one in current

QC.jl ships three canonicalisation routines (`canonicalization.jl`):

1. **`canonicalize!`** (line 92) — standard reduced row-echelon over GF(2), two-pass: first sweep X-pivots column by column, then Z-pivots. O(r²·n) where r = rank, n = qubits. Returns `(rx, rz)`, the X-pivot and Z-pivot ranks.
2. **`canonicalize_rref!`** (line 150) — RREF along a *subset* of columns, indexing in *reverse*. Specialised for `traceout!` (partial trace).
3. **`canonicalize_gott!`** (line 221) — Gottesman canonical form with column permutations. Used internally to compute logicals when constructing `MixedDestabilizer`. Also supports a `backtrack=true` mode that records the sequence of `rowswap`/`mul_left` ops so they can be undone.

Stim folds canonicalisation into `Tableau::stabilizers(canonical=true)` (`tableau.h:240`); the internal logic in `tableau.inl` is the same Gaussian-elimination algorithm.

### 6.1 What `PauliStabilizer.m` has

The `["Circuit"]` property (line 344) implements the **greedy AG canonicalization**: for each qubit, set destabilizer-X column to a single 1, then zero out the destabilizer-Z column, then zero out the stabilizer-X column. The deep audit notes this is "Standard, well-tested algorithm. Performance: `O(n³)` gates in the worst case; circuit length is **not minimized**" — better synthesis exists per [Reid24](https://arxiv.org/abs/2404.19408) and [Winderl23](https://arxiv.org/abs/2309.08972).

That's the only canonicalization-shaped routine. There is no `canonicalize_rref!` analogue (so partial trace, if added, will need it) and no separate Gottesman form with logical extraction (so the `LogicalX`/`LogicalZ` recommendation above will need a routine of its own).

**WL port note for what's missing:** `RowReduce[m, Modulus -> 2]` is a one-liner, but it doesn't track phases — you need to either fold a "phase column" into the matrix or implement the Gottesman row-op loop by hand. The latter is ~50 lines and is what the current `["Circuit"]` already does internally.

## 7. Projection / measurement — Z-basis only, single-qubit primary

This is the single most-used non-trivial routine. Both reference packages implement the standard Aaronson–Gottesman projection algorithm:

1. Test whether the measurement Pauli `M` commutes with all stabilizer generators.
2. If **yes** (deterministic outcome): the result is determined by the current sign data. Requires either a generation step (`generate!`) — O(n²) — or, with destabilizer info, a single XOR walk — O(n).
3. If **no** (random outcome): pick the first anticommuting generator, multiply every *other* anticommuting generator by it (so they now commute), replace the picked generator with `M`, and randomise its sign.

QC.jl `project_trace_reset.jl:276–308`, the `Stabilizer` version:

```julia
function _project!(stabilizer, pauli; keep_result, phases)
    anticommutes = 0
    for i in 1:n
        if comm(pauli, stabilizer, i) != 0x0
            anticommutes = i; break
        end
    end
    if anticommutes == 0           # deterministic branch
        if keep_result
            _, _, r = _canonicalize!(stabilizer; phases, ranks=Val(true))   # O(n³)
            gen = _generate!(copy(pauli), stabilizer; phases)               # O(n²)
            ...
        end
    else                            # random branch
        for i in anticommutes+1:n
            if comm(pauli, stabilizer, i) != 0
                mul_left!(stabilizer, i, anticommutes; phases)
            end
        end
        stabilizer[anticommutes] = pauli
        result = nothing
    end
end
```

The `MixedDestabilizer` version (lines 375–439) is structurally similar but cheaper: with destabilizers around, the deterministic-outcome value is computed via XORing in the destabilizers that anticommute with `pauli` — no canonicalisation needed, O(n) instead of O(n³).

Stim's equivalent is `TableauSimulator::measure_kickback_z` at `tableau_simulator.h:206`. Returns `(measurement_bit, kickback_pauli)`. The kickback is empty if deterministic, non-empty if random.

### 7.1 What `PauliStabilizer.m` has

`ps["Measure", a]` / `ps["M", a]` (line 291) implements exactly the AG projection algorithm for **single-qubit Pauli-Z measurement** on qubit `a`, using the destabilizer half for O(n) deterministic measurement (already efficient). Returns an `Association`:
- Deterministic case: `<|0 -> ps_post|>` or `<|1 -> ps_post|>`.
- Non-deterministic case: `<|0 -> ps0, 1 -> ps1|>` — the conditional branches.

`ps["M", {q1, q2, …}]` (line 314) recurses, returning a flat `Association` keyed by tuples of outcomes.

Convenience overloads `ps["M", q1, q2, …]`, `ps[q1, q2, …]`, and `ps[]` (lines 329–332) cover the common variadic and "measure all" cases.

### 7.2 What's missing — and how to add it

The deep audit calls out three gaps (§6.5):

1. **Pauli-X / Pauli-Y / arbitrary Pauli-string measurement.** To measure `XZZXI` (e.g., a 5-qubit-code stabilizer), the user must rotate first. This is a **major limitation for QEC stabilizer-syndrome workflows** — every stabilizer code defines its parity checks as Pauli strings, not single-qubit Z's. Adding it is a thin wrapper:

   ```mathematica
   ps["M", pauliString_String] := Module[{rotated, m},
       (* basis-rotation Cliffords: H for X-positions, H S H for Y-positions *)
       rotated = applyBasisRotation[ps, pauliString];
       m = rotated["M", firstNonIdentityQubit[pauliString]];
       Map[#["Dagger" /@ basisRotationCircuit] &, m]
   ]
   ```
   ~30 lines including the basis-rotation helper. The single biggest functional addition the kernel needs.

2. **Expectation value `⟨P⟩` for arbitrary Pauli `P`.** The AG framework supports this in O(n²) via the same commute-with-stabilizers test used in measurement, but without the projection step:
   - If `P` commutes with all stabilizers: `⟨P⟩ ∈ {±1}` based on the sign data.
   - If `P` anticommutes with any stabilizer: `⟨P⟩ = 0`.
   ~20 lines, reuses the existing `g`/`comm` primitives.

3. **Sampling.** `ps["M", q]` returns conditional outcomes, not a sampled outcome. The user has to write `RandomChoice[Keys[result] -> Values[result]]`. There is no built-in batched sampler — this is what §9 (Pauli frames) addresses.

### 7.3 Avoid the O(n³) trap

The deep audit confirms `["Circuit"]` is O(n³) and that the dense composition path at line 389 is also O(n³). For measurement, the existing single-qubit Z-measure correctly uses the destabilizer half and is O(n). When adding Pauli-string measurement, **do not** route through `Circuit . measureZ . Circuit^†` — that materializes the rotation as Cliffords and triggers the O(n³) composition. Instead, apply the basis rotations directly as tableau mutations (using the existing `["H", j]`, `["S", j]` updates, which are O(n) each).

## 8. The inverted-tableau trick — not needed

The Stim README articulates this clearly:

> When doing general stabilizer simulation, Stim tracks the inverse of the stabilizer tableau that was historically used. This has the unexpected benefit of making measurements that commute with the current stabilizers take linear time instead of quadratic time. This is beneficial in error correcting codes, because the measurements they perform are usually redundant and so commute with the current stabilizers.

Confirmed at `tableau_simulator.h:41`: the simulator's only state field is `Tableau<W> inv_state`. Gates conjugate the inverse (`prepend_*` on the inverse equals an `append` on the forward). Measurements use `inverse_x_output(qubit)` / `inverse_z_output(qubit)` (`tableau.h:228–238`) — these read out a single column of the inverse tableau in O(n) time.

QC.jl does **not** use the inverted tableau; it uses the destabilizer trick instead, which buys the same O(n) deterministic measurement at the cost of doubling the row count.

### 8.1 PauliStabilizer.m is already on the right side of this

The kernel's rank-3 tableau stores 2n rows always — destabilizers `1..n`, stabilizers `n+1..2n`. This is the destabilizer-augmented form, semantically equivalent to QC.jl's `MixedDestabilizer`. Deterministic measurement is therefore already O(n) (line 291 walks the destabilizer half on the deterministic branch).

**Recommendation: do not port the inverted-tableau trick.** Stim's inversion saves a constant factor *on top of* the destabilizer trick; for WL where the inner loop is already 10–100× slower than C, that constant factor isn't where the win is. Porting it would also require plumbing the `prepend` vs `append` distinction through the gate layer (lines 216–319), which is a real cost and zero benefit for this codebase's regime.

## 9. Pauli-frame / reference-frame sampling — the largest functional gap

Both reference packages implement the same idea:

> **Reference Frame Sampling.** When bulk sampling, Stim only uses a general stabilizer simulator for an initial reference sample. After that, it cheaply derives as many samples as needed by propagating simulated errors diffed against the reference. This simple trick is *ridiculously* cheaper than the alternative: constant cost per gate, instead of linear cost or even quadratic cost. (Stim README)

Mechanically, a "Pauli frame" is just a Pauli operator that you propagate alongside the reference state. When a noisy gate fires, you XOR a random error into the frame (instead of running the noise through the full stabilizer machinery). When a measurement fires, you XOR the frame's bit on that qubit into the reference outcome to get the actual outcome.

QC.jl `pauli_frames.jl:9–12`:

```julia
struct PauliFrame{T,S} <: AbstractQCState
    frame::T            # Tableau, one row per trajectory
    measurements::S     # Bool matrix [trajectory × measurement_index]
end
```

Each row of `frame` is one shot's worth of accumulated Pauli error. Gates are applied with `_apply!(frame, op; phases=Val(false))` (line 59) — the `phases=false` is the throughput knob from §4.3. Measurements (lines 85–97) are a single bit-test against the frame:

```julia
should_flip = !iszero(xzs[ibig, f] & ismallm)
frame.measurements[f, op.bit] = should_flip
```

Noise (lines 136–155) is one `rand()` per qubit per shot:

```julia
if r < p
    xzs[ibig, f] ⊻= ismallm                       # X error
elseif r < 2p
    xzs[end÷2 + ibig, f] ⊻= ismallm               # Z error
elseif r < 3p
    xzs[ibig, f] ⊻= ismallm                       # Y = X·Z
    xzs[end÷2 + ibig, f] ⊻= ismallm
end
```

**This is the bulk of the runtime in any QEC simulation, and it's all XOR.** Stim's `FrameSimulator` (in `src/stim/simulators/frame_simulator.{h,cc}`) is the same algorithm with SIMD lanes wider than QC.jl's.

### 9.1 What `PauliStabilizer.m` has

Nothing. The deep audit is explicit (§6.5, §13.3):

> Sampling — `ps["M", q]` returns conditional outcomes, not a sampled outcome. The user has to do `RandomChoice[Keys[result] -> Values[result]]` themselves. There is no built-in batched sampler — Stim's headline feature ([arxiv:2103.02202](https://arxiv.org/abs/2103.02202)) is not implemented.

This is the single biggest functional gap in the kernel relative to either reference implementation.

### 9.2 What to add

Implement the QC.jl `PauliFrame` interface directly. A WL sketch:

```mathematica
PauliFrame[ps_PauliStabilizer, nShots_Integer] := <|
  "ReferenceState" -> ps,
  "Frame"  -> ConstantArray[0, {nShots, 2 ps["Qubits"]}],  (* X-bits then Z-bits, all zero initially *)
  "Outcomes" -> {}
|>;

ApplyToFrame[frame_, "H", q_] := MapAt[
  Function[bits, swapAtPosition[bits, q, q + nQubits]],   (* H swaps X and Z *)
  frame,
  "Frame"
];

(* Noise: each shot, with probability p, XOR an X error into the frame *)
NoisyApply[frame_, gate_, p_] := AssociateTo[frame, "Frame" ->
  MapIndexed[
    With[{r = RandomReal[]},
      Which[
        r < p, BitXor[#1, xPositionMask[gate]],
        r < 2 p, BitXor[#1, zPositionMask[gate]],
        r < 3 p, BitXor[#1, BitOr[xPositionMask[gate], zPositionMask[gate]]],
        True, #1
      ]
    ] &,
    frame["Frame"]
  ]
];

(* Measurement: one XOR + bit-test per shot *)
MeasureFrame[frame_, q_] := With[{
    refOutcome = First[Keys[frame["ReferenceState"]["M", q]]],
    flips = Map[BitGet[#[[q]], 0] &, frame["Frame"]]
  },
  refOutcome ⊕ flips
];
```

The implementation is ~150 lines of WL. The throughput target is "constant per-gate cost across all shots" — i.e., 1000 shots run in roughly the same time as 10 shots, modulo `Length`-of-list overhead.

**Recommendation:** Implement Pauli frames as the *primary* sampling path. Even at WL speeds, the constant cost per gate (vs. linear in n) keeps it usable for distance ~7–13 surface codes with thousands of trajectories — i.e., the regime that fits a textbook chapter or a research note. The kernel's missing this is what blocks the QEC chapter of the book from running on real-scale examples.

## 10. Symbolic gates: already done, but only numeric

This is the section where Wolfram is naturally *better* than C++.

QC.jl `symbolic_cliffords.jl:91–118`:

```julia
macro qubitop1(name, kernel, inv_kernel)
    prefixname = Symbol(:s, name)
    quote
        struct $(esc(prefixname)) <: AbstractSingleQubitOperator
            q::Int
        end
        @inline $(esc(:qubit_kernel))(::$prefixname, x, z) = $kernel
        @inline $(esc(:inv_qubit_kernel))(::$prefixname, x, z) = $inv_kernel
    end
end

@qubitop1 Hadamard     (z   , x   , x!=0 && z!=0) (z   , x   , x!=0 && z!=0)
@qubitop1 Phase        (x   , x⊻z , x!=0 && z!=0) (x   , x⊻z , x!=0 && z==0)
@qubitop1 X            (x   , z   , z!=0)         (x   , z   , z!=0)
...
```

Stim's equivalent — `prepend_H_XZ`, `prepend_H_YZ`, `prepend_SQRT_XX`, ~80 of them — are hand-coded one-by-one in `tableau_specialized_prepend.cc`. QC.jl's macro generates the same set in ~50 lines; Stim's hand-coded set is several thousand lines of nearly-identical code.

### 10.1 What `PauliStabilizer.m` does already

Lines 216–319 implement gate updates as **direct tableau mutations** — exactly the kernel-association pattern, written by hand rather than via a macro. From the deep audit (§4):

| Gate | Property | Line | Update rule |
|---|---|---|---|
| H | `ps["H", j]` | 216 | Swap `Tableau[[1,j]]` ↔ `Tableau[[2,j]]`; flip sign where `x[j] == z[j] == 1` |
| S | `ps["S", j]` | 222 | `Z[j] ← Z[j] ⊕ X[j]`; flip sign where `x[j] == z[j] == 1` |
| S† | `ps[SuperDagger["S"], j]` | 228 | Same Z-update; sign flip on `x[j] == 1, z[j] == 0` |
| CNOT/CX | `ps["CNOT", j, k]` | 234 | `X[k] ← X[k] ⊕ X[j]`; `Z[j] ← Z[j] ⊕ Z[k]`; AG sign rule |
| X / Y / Z | lines 243–245 | direct phase update via `Phase ← Phase ⊕ Z[j]` etc. |
| CZ / SWAP | 246–247 | decomposed (`H_k . CNOT_jk . H_k`; `PermuteQudits`) |

This is **architecturally identical** to QC.jl's symbolic kernels — the WL natural form is direct pattern-dispatched function definitions, which is even cleaner than Julia's macro-generated structs. The current code is a direct translation of the same algebraic content.

### 10.2 What's not there

1. **Symbolic entries.** Every kernel assumes integer `{0, 1}` bits and integer `{-1, +1}` signs. Symbolic `θ`, `α`, etc. as tableau entries would require auditing every gate's sign-flip predicate — most are written as `Pick[..., x_row * z_row, 1]` or equivalent, which depends on multiplication being numeric.

2. **CY.** The deep audit (§4) flags: "`"CY"` is **NOT** a Clifford generator here. There is no built-in CY rule." Trivial to add (CY = (I⊗S†)·CNOT·(I⊗S)), but reflects that the gate set was assembled by hand rather than by closure under conjugation.

3. **Parametric Cliffords.** No `Rx[θ]`, `Ry[θ]`, `Rz[θ]` — even at θ ∈ {π/2, π, 3π/2} where they reduce to combinations of H and S. The user composes them manually.

4. **Closure of the non-Clifford output.** The `ps["P"[θ], j]` (line 321) and `ps["T", j]` (line 324) handlers return a *symbolic linear combination* of two `PauliStabilizer` objects:

   ```mathematica
   ps_PauliStabilizer["P"[phase_], j_Integer] := With[{c = Exp[I phase / 2]},
       (1 + c) / 2 ps + (1 - c) / 2 ps["Z", j]
   ]
   ```

   This is **not closed** under further gate application — `(...)["H", j]` does not distribute over the sum. The deep audit calls this "a token gesture, not real symbolic-phase support." Users will hit this edge the first time they try to compose a T with anything else.

### 10.3 The path forward for symbolic capability

This is the area the project owner has flagged as priority. Three options, in increasing order of work:

**Option A: Distribute non-Clifford outputs over further gate applications.** Add a `pSum /: pSum["H", j]` upvalue (or equivalent) so that the sum form survives composition. This makes T-gate compounds *expressible* but the representation grows exponentially with the number of T's — `k` T-gates give `2^k` terms. Useful for hand-checked small cases (k ≤ 4 or so).

**Option B: Dyadic-rational sign tracking.** Loosen `Signs ∈ {-1, +1}` to `Signs ∈ {Exp[I θ_k]}` where `θ_k` is a dyadic rational of π. This admits T-gate phases natively at the cost of replacing the integer sign list with a list of phase exponents. Most existing gate updates port directly — the predicates that flip signs become predicates that add π to the phase exponent.

**Option C: Full symbolic phase per [SymPhase / FangYing23](https://arxiv.org/abs/2311.03906).** The paper introduces a "stabilizer formalism with symbolic phase" — every sign is a symbolic expression over the group of phases generated by the Clifford+T gate set. Composition does symbolic phase arithmetic; measurement reduces the phase via a normal form. This is the most expressive and the largest project (~500 lines including the phase-simplification primitives).

**Recommendation:** Skip Option A (the exponential blow-up makes it a footgun). Implement Option B as a 0.2-style refactor — it's the cheapest path that gets symbolic T-gates into the formalism and stays closed under further composition. Layer Option C on top for an 0.3 release if the symbolic-phase regime turns out to be the dominant use case.

The implication for `PauliTableauQ` (line 12): the validation should be split — keep the strict `{0, 1, -1}` predicate as the fast-path check, and add a permissive `Symbolic` predicate for the symbolic-phase mode. Operations dispatch on whichever matches.

## 11. Random Clifford sampling — `PauliStabilizer.m` is ahead

Both reference packages implement uniform sampling on the symplectic group, **not** Bravyi–Maslov / Mallows-distribution sampling:
- QC.jl has `random_invertible_gf2`, `random_pauli`, and `random_clifford` — the last samples a uniformly random invertible symplectic matrix and converts (`src/randoms.jl`).
- Stim's `Tableau<W>::random` (declared at `tableau.h:86`) is the same approach.

`PauliStabilizer.m` (line 165) implements the **Bravyi–Maslov / Koenig–Smolin Mallows-distribution sampler** ([KoeSmo14, arxiv:1406.2170](https://arxiv.org/abs/1406.2170)). `SampleMallows[n]` (line 147) is the underlying permutation-with-Mallows-bias sampler. Exposed as `PauliStabilizer["Random", n]` (line 184); needs O(n²) random bits and O(n³) matrix work for the binary inverse.

This is the one place where the WL kernel is ahead of both reference implementations on algorithmic substance, not just expressiveness. The Mallows sampler is mathematically the right uniform-on-the-Clifford-group sampler; the symplectic-uniform approaches in QC.jl and Stim are uniform modulo Pauli right-multiplication, which is *almost* uniform on Cliffords but not quite. For pedagogical correctness this matters.

**Recommendation:** Keep the existing implementation. No upgrade needed.

## 12. API surface — what's there, what should be there

QC.jl exports ~80 functions; the user-facing core is much smaller:

```
@P_str, @S_str, @T_str, @C_str               # literal macros
PauliOperator, Stabilizer, MixedDestabilizer # types
apply!, project!, projectrand!               # state evolution + measurement
canonicalize!                                # normal form
mul_left!, mul_right!                        # Pauli algebra
sH, sS, sCNOT, sCPHASE, sSWAP                # symbolic gates (~30)
PauliFrame, pftrajectories                   # bulk sampling
nqubits, stabilizerview, destabilizerview    # accessors
random_pauli, random_clifford, random_stabilizer
graphstate                                    # graph-state form
expect, traceout!, ptrace                    # observables, partial trace
```

Stim's user surface (Python via pybind11) is even tighter:
```
stim.PauliString, stim.Tableau               # algebra
stim.Circuit                                  # IR
stim.TableauSimulator                         # interactive simulation
stim.FrameSimulator                           # bulk sampling
stim.DetectorErrorModel                       # downstream decoder input
```

### 12.1 What `PauliStabilizer.m` exposes

The deep audit notes (§1):

> **Properties contract** (line 17): `{"Qubits", "Generators", "Matrix", "Phase", "TableauForm", "State", "Operator", "Circuit"}`. This list is **incomplete vs. what's actually defined**: at least 25 more property names dispatch through the property handlers.

The actual property surface is much wider:
- Storage / structure: `"Qubits"`, `"Qudits"`, `"GeneratorCount"`, `"Signs"`, `"Tableau"`, `"X"`, `"Z"`, `"Phase"`, `"Matrix"`, `"p"`, `"TableauPhase"`
- Splits: `"Stabilizer"`, `"StabilizerTableau"`, `"StabilizerX"`, `"StabilizerZ"`, `"StabilizerSigns"`, `"Destabilizer"`, `"DestabilizerTableau"`, `"DestabilizerX"`, `"DestabilizerZ"`, `"DestabilizerSigns"`
- Display: `"PauliForm"`, `"Generators"`, `"Stabilizers"`, `"Destabilizers"`, `"PauliStrings"`, `"PauliSymbols"`, `"StabilizerTableauForm"`, `"TableauForm"`
- Conversion: `"State"`, `"Operator"`, `"Circuit"`
- Operations as properties: `"H"`, `"S"`, `"CNOT"`, `"CZ"`, `"SWAP"`, `"V"`, `"Permute"`, `"PermuteQudits"`, `"PadLeft"`, `"PadRight"`, `"Dagger"`, `"Inverse"`, `"M"`, `"Measure"`, `"RowSum"`, `"P"[θ]`, `"T"`

### 12.2 What's missing from the API

| Capability | Why it matters | Cost to add |
|---|---|---|
| `ps["M", pauliString]` | Pauli-string measurement; the syndrome-extraction primitive | ~30 lines on top of existing primitives |
| `ps["Expect", pauliString]` | Expectation values; the most-asked observable | ~20 lines, reuses `g`/`comm` |
| `ps["Sample", n]` or `Sample[circuit, shots]` | Pauli-frame backed bulk sampling | ~150 lines (§9) |
| `ps["LogicalX", k]` / `ps["LogicalZ", k]` | Logical operator access for QEC | ~50 lines, reuses `["Circuit"]` canonicalization machinery |
| `ps["GraphState"]` | Canonical form per [AndBri05](https://arxiv.org/abs/quant-ph/0504117) | ~80 lines |
| `ps["InnerProduct", other]` | `⟨ps_a\|ps_b⟩` via [GarMarCro12](https://arxiv.org/abs/1210.6646) | ~100 lines |
| `ps["TraceOut", qs]` / `ps["PartialTrace", qs]` | Standard QM operation; needs `canonicalize_rref!`-style routine | ~80 lines |
| `Rx[θ]`, `Ry[θ]`, `Rz[θ]` as parametric Cliffords | Reduce to H/S compositions at θ ∈ {π/2, π, 3π/2} | ~30 lines |

### 12.3 Properties contract maintenance

Three small cleanups from the deep audit (§13.7):
- The `"Properties"` contract at line 17 is out of date — it lists 8 of the 25+ actual properties. Keep it in sync, or compute it programmatically from the dispatch rules.
- `"QuantumSttate"` (line 335) — typo, likely intended `"QuantumState"`. Either fix or remove.
- `$PauliStabilizerNames` (line 9) lists `"SteaneCode"` twice. Cosmetic but visible.

## 13. Concrete recommendations for upgrading `PauliStabilizer.m`

Ordered by priority. This replaces the original from-scratch v0.1/v0.2/v0.3 plan; the kernel already has the v0.1 capability set.

### Tier 1 — Functional gaps blocking real use

1. **Pauli-string measurement.** `ps["M", "XZZXI"]`, with sign-aware basis rotations. The single most important addition for QEC workflows. ~30 lines.

2. **Expectation value.** `ps["Expect", "XZZXI"]`, returns ±1 (commuting case) or 0 (anticommuting case). ~20 lines.

3. **Pauli-frame bulk sampling.** Implement `PauliFrame[ps, nShots]` per §9. Required for any non-trivial QEC simulation; without it, the kernel cannot scale past ~50-shot pedagogical demos. ~150 lines.

4. **Tests file.** `Tests/PauliStabilizer.wlt` does not currently exist (deep audit §13.5). The minimum suite should lock down: every gate's update on a known initial state, the measurement outcomes for the 5/7/9-qubit codes, the AG canonicalization round-trip `ps == PauliStabilizer[ps["Circuit"]]`, the Mallows random-Clifford reproducibility under `SeedRandom`, and the integration UpValues with `QuantumState`/`QuantumOperator`. ~300 lines of test code.

5. **Documentation.** No `Usage.m` entry, no documentation page (deep audit §13.6). The 3 demo notebooks (`Stabilizers.nb`, `StabilizerError.nb`, `PauliGroupTheory.nb`) reference `PauliStabilizer` zero times — they use raw `QuantumOperator["XZZXI"]` instead. Either rewrite the notebooks or add a focused tutorial that exercises the main API.

### Tier 2 — Symbolic extensions (the project owner's flagged priority)

Per §10.3, the proposed path is:

6. **Loosen sign storage to dyadic-rational phases.** Replace `Signs ∈ {-1, +1}` with `Signs` as a list of phase exponents over `Z[1/2^k]·π`. T-gate composition stays in the formalism, exponentially in `k` rather than as a sum-over-Cliffords. Concretely, `PauliStabilizerQ` accepts the new form; gate updates that previously did `signs[[i]] *= -1` now do `signs[[i]] += π`; phase normalization happens at export.

7. **Distribute or close `"P"[θ]`/`"T"`.** With Tier 2.6 in place, `ps["P"[θ], j]` no longer needs to escape to a sum. The phase update is `Signs[[stabilizers carrying Z on j]] += θ`. Direct, closed under further gates.

8. **Parametric Clifford gates.** `Rx[θ]`, `Ry[θ]`, `Rz[θ]` for θ ∈ {π/2, π, 3π/2, 2π}; auto-reduce to H/S compositions. `RX[π/2]` = `H S H`. ~30 lines.

9. **Symbolic measurement outcomes.** The `Bernoulli[p]` form for the random branch — output `<|0 -> ps0, 1 -> ps1, "Probability" -> 1/2|>` so users can compute expectation values of measurement-conditioned observables symbolically.

### Tier 3 — Performance

10. **Bit-packed UInt64 storage as a parallel representation.** Add a `ps["PackedXZ"]` accessor that returns packed UInt64 lists, computed lazily and cached. Hot kernels (composition, measurement) take a fast path on the packed view; the rank-3 array remains canonical. ~200 lines of WL.

11. **LibraryLink shim for `mul_pauli_packed`.** Only after Tier 3.10 ships and is profiled, and only if a real workload demands it. ~150 lines of C, gets within ~3× of Stim.

12. **Caching of computed properties.** `["Matrix"]`, `["Phase"]`, etc. are recomputed every call (deep audit §13.4). Cache via the `data` Association inside the head.

### Tier 4 — Scope expansion

13. **Logical operator access.** `ps["LogicalX", k]`, `ps["LogicalZ", k]`. The information is in the destabilizer half once the state is in Gottesman canonical form. ~50 lines.

14. **Graph-state extraction.** `ps["GraphState"]` per [AndBri05](https://arxiv.org/abs/quant-ph/0504117). Run `canonicalize_gott!` and read off the off-diagonal pattern. ~80 lines.

15. **Inner product.** `InnerProduct[ps_a, ps_b]` per [GarMarCro12](https://arxiv.org/abs/1210.6646). O(n³). ~100 lines.

16. **Partial trace.** `ps["PartialTrace", qubits]` — needs the `canonicalize_rref!`-style subset-canonicalization routine that QC.jl ships separately. ~80 lines.

17. **Qudit support.** `PauliRow` (line 58) already takes a `d` parameter; promote it through the API. Per [HosDehMoo04](https://arxiv.org/abs/quant-ph/0408190), [Beaudrap11](https://arxiv.org/abs/1102.3354). The string parser, the named-codes dispatch, and the `["State"]` materialization (`PauliMatrix[0..3]` hard-coded at line 335) all need generalization.

18. **Better synthesis.** Replace the greedy AG canonicalization in `["Circuit"]` (line 344) with [Reid24](https://arxiv.org/abs/2404.19408) or [Winderl23](https://arxiv.org/abs/2309.08972) for shorter circuits. Optional; the current output is correct, just not minimal.

### What to skip

- **Inverted tableau (Stim).** The destabilizer form (already present) gives the same O(n) deterministic measurement asymptotic. Stim's inversion saves a constant factor that doesn't translate to WL.
- **Bravyi–Maslov sampling** *as a separate addition*. Already implemented at line 165.
- **Full DEM (detector error model, Stim).** It's a downstream-decoder format; if QuantumFramework doesn't ship a decoder, it doesn't need DEM.
- **ECC code-family generators.** `lib/QECCore/` in QC.jl is ~100 files (Reed-Muller, BCH, lifted-product, etc.) — a separate package.
- **GPU paths.** Both reference packages have them; not appropriate scope here.

## 14. Correctness traps observed during the audit

These are subtle points where both reference packages had to get something exactly right that a fresh implementation will likely get wrong. For each, I note the corresponding location in `PauliStabilizer.m`.

1. **Phase tracking under multiplication uses TWO parity bits, not one.** A naïve "every anticommute contributes ±i" loop is off by half-integer phases. Both QC.jl and Stim implement carry-save addition over two parity bits. **In `PauliStabilizer.m`:** the AG `g` function (line 280) implements this directly, term-by-term over the integer tableau. The phase carry logic is **algebraically equivalent** to QC.jl's `cnt2 ⊻= (cnt1 ⊻ ...) & anti_comm; cnt1 ⊻= anti_comm`, just expressed as an integer formula. Already correct; verify before changing.

2. **`mul_left!` on a `Destabilizer` flips index order.** In QC.jl `mul_leftright.jl:200–206`, the Destabilizer overload calls `mul_left!(t, j, i; phases=Val(false))` instead of `(t, i, j)` to preserve commutation between stabilizers and destabilizers. **In `PauliStabilizer.m`:** the corresponding invariant lives in the `["Circuit"]` canonicalization (line 344) — the order of operations within the greedy loop matters. Worth a focused test.

3. **`Stabilizer` projection requires a full canonicalisation in the deterministic branch — but `MixedDestabilizer` doesn't.** **In `PauliStabilizer.m`:** the kernel is structurally a MixedDestabilizer, and the Z-measure path (line 291) correctly uses the destabilizer half. Don't accidentally route a future upgrade through a forward-tableau-only path.

4. **Pauli-frame initialisation needs a uniform-random Z component.** In QC.jl `pauli_frames.jl:47–55`, `initZ!` randomises the Z bits of every frame at start. Without it, all trajectories start identical and you get systematically biased samples. **In `PauliStabilizer.m`:** when implementing Tier 1.3 (Pauli frames), the initial frame must have `RandomInteger[1, {nShots, n}]` filled into the Z half, not zeros.

5. **Stim's "no imaginary phase in storage" is an algebraic claim, not a convenience.** **In `PauliStabilizer.m`:** the integer signs `{-1, +1}` are inherently ±1-only, so this trap doesn't apply. But if Tier 2.6 lands and signs become symbolic phases, the symbolic-phase representation must enforce or verify the same invariant — half-integer angles within a single stabilizer row will silently corrupt the tableau.

6. **Tableau orientation matters for prepend vs append.** Stim's `Tableau` is column-major (`tableau.h:47–49`); QC.jl's `Tableau.xzs` is chunk-major-then-row. **In `PauliStabilizer.m`:** the rank-3 `{2, n, 2n}` array has the "row-of-2n" axis last, which is fine for WL's natural `[[All, All, k]]` slice extraction. Don't refactor the orientation; the gate updates assume it.

7. **Signs validation rejects symbolic outputs from `"P"[θ]`.** **In `PauliStabilizer.m`:** the current behavior is that `ps["P"[θ], j]` returns `<expr> + <expr>` outside the `PauliStabilizer` head, so `PauliStabilizerQ` is never asked. If Tier 2.7 lands and the result stays a `PauliStabilizer`, `PauliStabilizerQ` (line 13) will reject the symbolic phases — so it must be loosened in lockstep.

8. **`["State"]` and `["Operator"]` materialize 2ⁿ×2ⁿ matrices.** Practical limit is ~n=10 (deep audit §13.4). The TraditionalForm summary (line 480) defers to `ToBoxes[ps["State"], TraditionalForm]` — so calling `TraditionalForm[ps]` with n > 10 will hang or OOM the kernel. Either guard the formatter at line 480 with the `GeneratorCount < 32` check used elsewhere, or never auto-materialize.

9. **`PauliStabilizerApply` silently passes through unknown gates** (deep audit §10, line 426). A gate that isn't in the dispatch table simply returns unevaluated, with no error. This is brittle: a typo in a circuit description corrupts the simulation invisibly. Add an explicit guard.

## 15. File map (citation index)

For when something in the recommendations needs to be cross-checked back to source.

### `PauliStabilizer.m` (in-tree, the upgrade target)

| Topic | Lines |
|---|---|
| `PauliTableauQ` validation (rank-3, integer entries) | 12 |
| `PauliStabilizerQ` validation (signs ∈ {-1,+1}) | 13 |
| Properties contract list (incomplete) | 17 |
| `PauliStabilizerTableau` (4ⁿ Pauli-expectation tomography) | 22 |
| `PauliRow` (the qudit-capable inner kernel) | 58 |
| Constructor: `QuantumOperator → PauliStabilizer` | 88–92 |
| Constructor: `QuantumCircuitOperator → PauliStabilizer` | 94–97 |
| Constructor: rank-3 array, single-row pad | 99–107 |
| Constructor: Pauli-string list | 109–119 |
| Constructor: explicit destab+stab string halves | 121–129 |
| Constructor: `q_Integer` ⇒ `\|0…0⟩` | 131 |
| Named codes dispatch | 133–140 |
| `SampleMallows` (Mallows permutation sampler) | 147 |
| `RandomClifford[n]` (Bravyi–Maslov / Koenig–Smolin) | 165 |
| `PauliStabilizer["Random", n]` | 184 |
| Property dispatch (direct) | 192–214 |
| Gate update: H | 216 |
| Gate update: S | 222 |
| Gate update: S† | 228 |
| Gate update: CNOT/CX | 234 |
| Gate updates: X / Y / Z | 243–245 |
| Gate updates: CZ / SWAP / Permute / PermuteQudits | 246–256 |
| Dagger / Inverse | 263 |
| Pad operations | 276–277 |
| AG `g` function | 280 |
| `rowsum` primitive | 282 |
| `ps["M", a]` / `ps["Measure", a]` | 291 |
| Multi-qubit measurement | 314 |
| Gate updates: V / V† | 318–319 |
| `ps["P"[θ], j]` (non-Clifford, returns sum-of-PS) | 321 |
| `ps["T", j]` and `ps[SuperDagger["T"], j]` | 324–325 |
| Order-argument rewriting | 327 |
| Convenience measurement overloads | 329–332 |
| `ps["State"]` (materializes 2ⁿ×2ⁿ) | 335 |
| `ps["Circuit"]` (greedy AG canonicalization) | 344 |
| `ps["Operator"]` | 375 |
| `phaseLookup` (precomputed sparse rank-4 array) | 380 |
| Composition `left[right]` (dense matrix mod 2) | 389 |
| `QuantumTensorProduct` | 408 |
| `PauliStabilizerApply` (the dispatcher) | 426 |
| UpValues for QF integration | 432–436 |
| `PauliForm` display | 441 |
| `TableauForm` display | 453 |
| Property names list | 469–475 |
| TraditionalForm summary (auto-materializes!) | 480 |
| SummaryBox with truncation | 484 |

### QuantumClifford.jl

| Topic | File | Lines |
|---|---|---|
| Pauli operator type, F(2,2) layout | `src/pauli_operator.jl` | 42–48 |
| `@P_str` macro, phase dictionary | `src/pauli_operator.jl` | 97–105 |
| In-place multiplication, non-vectorised | `src/mul_leftright.jl` | 6–31 |
| Vectorised multiplication via SIMD.jl | `src/mul_leftright.jl` | 75–121 |
| `mul_left!` for Destabilizer (index flip) | `src/mul_leftright.jl` | 200–206 |
| Top-level types, exports | `src/QuantumClifford.jl` | 1–110 |
| `Tableau` struct | `src/QuantumClifford.jl` | 145–151 |
| `canonicalize!` (Echelon) | `src/canonicalization.jl` | 92–133 |
| `canonicalize_rref!` (subset, reverse) | `src/canonicalization.jl` | 150–180 |
| `canonicalize_gott!` (with column perms) | `src/canonicalization.jl` | 221–272 |
| `_project!(Stabilizer)` | `src/project_trace_reset.jl` | 276–308 |
| `_project!(MixedDestabilizer)` | `src/project_trace_reset.jl` | 375–439 |
| `projectX/Y/Z!` specialisations | `src/project_trace_reset.jl` | 447–504 |
| `traceout!`, `ptrace` | `src/project_trace_reset.jl` | 642–924 |
| `expect` | `src/project_trace_reset.jl` | 744–752 |
| `PauliFrame` type | `src/pauli_frames.jl` | 9–37 |
| Pauli-frame gate application | `src/pauli_frames.jl` | 58–117 |
| Pauli-frame noise injection | `src/pauli_frames.jl` | 136–174 |
| Symbolic gate macro `@qubitop1` | `src/symbolic_cliffords.jl` | 91–103 |
| Single-qubit gate kernels | `src/symbolic_cliffords.jl` | 105–118 |
| `SingleQubitOperator` general form | `src/symbolic_cliffords.jl` | 168–249 |

### Stim

| Topic | File | Lines |
|---|---|---|
| `PauliString` type (xs, zs, sign) | `src/stim/stabilizers/pauli_string.h` | 71–142 |
| XZ ↔ XYZ encoding helpers | `src/stim/stabilizers/pauli_string.h` | 41–61 |
| `TableauHalf`, `Tableau` types | `src/stim/stabilizers/tableau.h` | 31–56 |
| `prepend_*` specialised gates (~40) | `src/stim/stabilizers/tableau.h` | 175–216 |
| Constant-time output queries | `src/stim/stabilizers/tableau.h` | 222–238 |
| `TableauSimulator`, `inv_state` field | `src/stim/simulators/tableau_simulator.h` | 39–45 |
| `do_*` gate dispatchers (~70) | `src/stim/simulators/tableau_simulator.h` | 109–180 |
| `measure_kickback_z/y/x` | `src/stim/simulators/tableau_simulator.h` | 206–208 |
| `peek_x/y/z`, postselect | `src/stim/simulators/tableau_simulator.h` | 186–197 |
| `peek_observable_expectation` | `src/stim/simulators/tableau_simulator.h` | 271 |
| `collapse_qubit_z` (random branch) | `src/stim/simulators/tableau_simulator.h` | 229 |

## 16. Coverage caveats

The Explore agents reported the following; I did not directly read the source for these and so they are second-hand. Nothing in the recommendations depends on these being exactly right, but flag if you decide to lean on them:

- Stim's `gen/`, `cmd/`, `diagram/` directories — only know structure from the `find` listing.
- DEM construction details (`src/stim/dem/`) — algorithm summary is from the agent.
- Stim's `search/graphlike/algo.h` shortest-error finder — the Explore call returned a short summary.
- QC.jl `ext/QuantumCliffordGPUExt/` GPU adapters — read the file list, not contents.
- The non-vectorised assertions about Stim's SIMD behaviour (256-bit alignment, "100 billion Pauli multiplications/second") come from Stim's own README, which I trust but did not benchmark.
- The Bravyi–Maslov non-implementation claim for QC.jl and Stim — the Explore agents say neither uses it, but the random-sampling code in `src/randoms.jl` and `tableau.inl` was not directly read. (`PauliStabilizer.m` definitely does use Mallows/Bravyi–Maslov; that's read directly at line 165.)
- The deep audit's claim that `Tests/PauliStabilizer.wlt` does not exist (§13.5) — confirmed by the deep audit but not re-verified for this report.

## 17. Bottom line

`PauliStabilizer.m` is a respectable v0.1 of a stabilizer simulator: Cliffords, deterministic and random Z-measurement, symplectic-mod-2 composition, AG canonical synthesis, Mallows random sampling, tensor product, and full integration with the rest of QuantumFramework via `QuantumState`/`QuantumOperator`/`QuantumCircuitOperator` UpValues. Architecturally it is a `MixedDestabilizer` with `rank == n`, written as a 494-line direct WL translation of the AG algorithms.

The four upgrade bands, in priority order:

1. **Bulk-sampling throughput** (Pauli-frame engine) is the single biggest functional gap; without it the kernel cannot run any QEC simulation past pedagogical scale. Tier 1.3 in §13.
2. **Measurement breadth** (Pauli-string measurement, expectation values) is a small, mechanical addition that unlocks the QEC-syndrome workflow that the kernel was designed for. Tier 1.1–1.2.
3. **Symbolic capability** (the project owner's flagged priority) requires loosening the integer-only contracts on `PauliTableauQ` and `Signs`, and either distributing or closing the non-Clifford `"P"[θ]`/`"T"` outputs. Tier 2.6–2.9. The dyadic-rational-phase representation (Option B in §10.3) is the cheapest path that gets symbolic T-gates into the formalism.
4. **Performance** (bit-packed UInt64, optional LibraryLink) is a follow-up that only matters once Tiers 1–2 expose a real workload — `n ≲ 50` is fine on the existing rank-3 integer tableau. Tier 3.

The Stim "inverted tableau" trick should not be ported — the destabilizer-augmented form already gives O(n) deterministic measurement. The QC.jl four-type lattice (`Stabilizer` / `Destabilizer` / `MixedStabilizer` / `MixedDestabilizer`) should not be ported as four separate heads — one `PauliStabilizer` head with rank tracking and logical-operator accessors layered on is the cleaner WL fit.

The expected total work, end-to-end through Tier 3, is roughly +1500 lines of WL on top of the existing 494, plus ~300 lines of tests and a documentation pass. Performance after Tier 3.10 (packed-UInt64 fast path) should land within ~10× of QC.jl single-threaded for `n ≤ 64`; closing the remaining gap requires LibraryLink (Tier 3.11) and is only justified if a real workload demands it.
